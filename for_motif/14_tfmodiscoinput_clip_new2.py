import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader
from torch import nn
import os
import datetime
import random
import math
from path_explain import PathExplainerTorch
import sys

args = sys.argv
looptime=0

samplename=str(args[1])
chunknum=1
enhaclusfile=str(args[2])
cellclusfile=str(args[3])
celltype=int(str(args[4]))
prefix=str(args[5])

ATACmatrix_norm=np.load(samplename+"/ATAC_pred_fold0.npy")
ATACmatrix_norm=ATACmatrix_norm[:,:,0]
ATACmatrix_norm=torch.from_numpy(ATACmatrix_norm)
ATACmatrix_norm=ATACmatrix_norm.to(torch.float32)
ATACmatrix_norm=(ATACmatrix_norm/ATACmatrix_norm.sum(axis=0))*(ATACmatrix_norm.sum(axis=0).mean())
ATACmatrix_norm=nn.functional.normalize(ATACmatrix_norm,p=2.0, dim=1, eps=1e-12, out=None)

RNAmatrix=pd.read_csv(samplename+"/log1praw.csv",sep=",")
RNAmatrix=RNAmatrix.to_numpy()
RNAmatrix=torch.from_numpy(RNAmatrix)
RNAmatrix=RNAmatrix.to(torch.float32)
#RNAmatrix=torch.log(RNAmatrix)
RNAmatrix_b=((RNAmatrix.transpose(0,1)-RNAmatrix.mean(dim=1))/RNAmatrix.std(dim=1)).transpose(0,1)

predRNA=np.load(samplename+"/predrna_t2_raw_2_fold0.npy")
predRNA=torch.from_numpy(predRNA)
predRNA=predRNA.to(torch.float32)

peaklist=pd.read_csv(samplename+"/peaks.bed",sep="\t",header=None)
peaklist=peaklist[[1,2]].to_numpy()
peaklist=torch.from_numpy(peaklist)
peaklist=peaklist.to(torch.float32)

RNAembed=pd.read_csv(samplename+"/pca_50dim_rnaembed_raw.txt",sep=" ",header=None)
RNAembed=RNAembed.to_numpy()
RNAembed=torch.from_numpy(RNAembed)
RNAembed=RNAembed.to(torch.float32)

pairlist_prom=pd.read_csv(samplename+"/pair_promoter.csv",header=None)
pairlist_prom=pairlist_prom[[1,2,3]].to_numpy()
pairlist_prom=torch.from_numpy(pairlist_prom)
pairlist_prom=pairlist_prom.to(torch.float32)

pairlist=pd.read_csv(samplename+"/pair_300000.csv",header=None,names=("gene","start","end","pos"))
pairlist_train=pd.read_csv(samplename+"/pair_300000_train_fold0.csv",header=None,names=("gene","start","end","pos"))
pairlist_test=pd.read_csv(samplename+"/pair_300000_test_fold0.csv",header=None,names=("gene","start","end","pos"))

traingenelist=pairlist["gene"].isin(pairlist_train["gene"]).values.tolist()
testgenelist=pairlist["gene"].isin(pairlist["gene"]).values.tolist()

pairlist=pairlist[["start","end","pos"]].to_numpy()
pairlist=torch.from_numpy(pairlist)
pairlist=pairlist.to(torch.float32)

peakreadlist=np.load(samplename+"/peak_extend_read.npy",allow_pickle=True)
peakreadlist=torch.from_numpy(peakreadlist)
peakreadlist=peakreadlist.to(torch.float32) #peaks,length,4
peakreadlist=peakreadlist.transpose(1,2) #peaks,4,length

embedlist=pd.read_csv(samplename+"/embed_fold0.txt",sep=" ",header=None)
embedlist=embedlist.to_numpy()
embedlist=torch.from_numpy(embedlist)
embedlist=embedlist.to(torch.float32)

cellclus=pd.read_csv(cellclusfile,header=None)
cellclus=cellclus.to_numpy()
cellclus=cellclus[:,0]
cellclus=torch.from_numpy(cellclus)
cellclus=cellclus.long()

print(cellclus.shape)

traincellnum=np.load(samplename+"/traincellnum.npy")

attr_list_full = []
input_list_full = []

enhaclus=np.load(enhaclusfile)
usecellmtx=np.load(samplename+"/enhausetagmtx.npy")


fullcell=RNAembed.shape[0]
traincell=traincellnum #10000 11830 11740
testcell=fullcell #fullcell-traincell #1893
#fullgenenum=traingenenum+testgenenum
fullgenenum=pairlist.shape[0]
traintag=0
testtag=0 #traincell

peaknum=ATACmatrix_norm.shape[0]
peaknum_train=round(peaknum*0.9)
trainpeaktag=0
testpeaktag=peaknum_train
trainpeaknum=peaknum_train
testpeaknum=peaknum-trainpeaknum

batchsize=1

trainpeaklist=list(range(trainpeaktag,trainpeaktag+trainpeaknum))
testpeaklist=list(range(testpeaktag,testpeaktag+testpeaknum))
fullpeaklist=list(range(peaknum))

max_len=int((pairlist[:,1]-pairlist[:,0]).max()+1)
PAD_IDX=0

b_epoch=0
epoch=1
epoch_final=0

fdim=512
edim=32
lamd=0
llam=0.001
alam=0.1

basenum=1500
basenuminput=1344
gradual=3
taun=0.3
mtk1=0.1
mtk2=0.1
alpha = 0.95

opt_p_b=1e-2
opt_p_f=1e-3

opt_e=1e-5
opt_t=1e-5
opt_r=1e-3

opt_e_f=1e-6
opt_t_f=1e-6
opt_r_f=1e-4

flipr=1
revr=0.5

traingenenum=1 #3599 #11881
testgenenum=math.ceil(fullgenenum/chunknum)

traingenelist=torch.from_numpy(np.array(range(traingenenum)))
if looptime!=(chunknum-1):
  testgenelist=torch.from_numpy(np.array(range(testgenenum*looptime,testgenenum*(looptime+1))))
else:
  testgenelist=torch.from_numpy(np.array(range(testgenenum*looptime,fullgenenum)))
  testgenenum=len(testgenelist)

cellclus_test=cellclus[testtag:testcell+testtag]



os.system("free -h")
torch.backends.cudnn.benchmark = True

def dataprep(train=True):
  inputmtx=[]
  prommtx=[]
  atacmtx=[]
  rnamtx=[]
  peakttmtx=[]
  posmtx=[]
  if train==True:
    inputmtx = [inputcal(traingenelist[j]) for j in range(traingenenum) ]
    prommtx = [promcal(traingenelist[j]) for j in range(traingenenum) ]
    atacmtx = [ataccal(traingenelist[j],traintag,traincell+traintag) for j in range(traingenenum) ]
    rnamtx= [rnacal(traingenelist[j],traintag,traincell+traintag) for j in range(traingenenum) ]
    posmtx = [poscal(traingenelist[j]) for j in range(traingenenum) ]
  else:
    inputmtx = [inputcal(testgenelist[j]) for j in range(testgenenum) ]
    prommtx = [promcal(testgenelist[j]) for j in range(testgenenum) ]
    atacmtx = [ataccal(testgenelist[j],testtag,testcell+testtag) for j in range(testgenenum) ]
    rnamtx= [rnacal(testgenelist[j],testtag,testcell+testtag) for j in range(testgenenum) ]
    posmtx = [poscal(testgenelist[j]) for j in range(testgenenum) ]
  return inputmtx,prommtx,atacmtx,rnamtx,posmtx


def inputcal(j):
  peakstart=int(pairlist[j,0].item())
  peakend=int(pairlist[j,1].item())
  prompos=int(pairlist_prom[j,0].item())
  enha=list(range(peakstart,(peakend+1)))
  enha.remove(prompos)
  genepos=int(pairlist[j,2].item())
  peaknum_gene=peakend-peakstart+1
  input = torch.zeros((max_len,peakreadlist.shape[1],peakreadlist.shape[2]))
  embedvec=peakreadlist[enha,:,:]
  input[0,:,:]=peakreadlist[prompos,:,:]
  input[1:peaknum_gene,:,:]=embedvec
  return input


def promcal(j):
  peakstart=int(pairlist_prom[j,0].item())
  embedvec=peakreadlist[peakstart,:,:]
  return embedvec


def ataccal(j,cellstart,cellend):
  cellnum=cellend-cellstart
  atacvec_return = torch.zeros((cellnum,max_len))
  peakstart=int(pairlist[j,0].item())
  peakend=int(pairlist[j,1].item())
  peaknum_gene=peakend-peakstart+1
  prompos=int(pairlist_prom[j,0].item())
  enha=list(range(peakstart,(peakend+1)))
  enha.remove(prompos)
  atacvec=ATACmatrix_norm[enha,cellstart:cellend]
  atacvec = atacvec.transpose(1,0)
  atacvec_return[:,0] = ATACmatrix_norm[prompos,cellstart:cellend]
  atacvec_return[:,1:peaknum_gene]=atacvec
  return atacvec_return



def rnacal(j,cellstart,cellend):
  rnavec=RNAembed[cellstart:cellend,:]
  #rnavec = rnavec.unsqueeze(-1)
  return rnavec


def poscal(j):
  peakstart=int(pairlist[j,0].item())
  peakend=int(pairlist[j,1].item())
  prompos=int(pairlist_prom[j,0].item())
  enha=list(range(peakstart,(peakend+1)))
  enha.remove(prompos)
  genepos=int(pairlist[j,2].item())
  peaknum_gene=peakend-peakstart+1
  input = torch.zeros((max_len))
  posvec=((((peaklist[enha,0]+peaklist[enha,1])/2)-genepos)/300000)
  input[0]=((((peaklist[prompos,0]+peaklist[prompos,1])/2)-genepos)/300000)
  input[1:peaknum_gene]=posvec
  return input



traininput,trainprom,trainatac,trainrna,trainpos=dataprep(True)
testinput,testprom,testatac,testrna,testpos=dataprep(False)

def peak_b_cal(j):
  input = peakreadlist[j,:,:]
  return input

def atac_b_cal(j):
  atacvec = ATACmatrix_norm[j,:]
  return atacvec


class CustomDataset(torch.utils.data.Dataset):
  def __init__(self,train=True):
    self.genes=[]
    self.trmd=[]
    self.atac=[]
    self.seqe=[]
    if train==True:
      self.genes=[traingenelist[j] for j in range(traingenenum)]
      self.trmd=[0 for j in range(traingenenum)]
      self.mode=[0 for j in range(traingenenum) ]
    else:
      self.genes=[testgenelist[j] for j in range(testgenenum)]
      self.trmd=[1 for j in range(testgenenum)]
      self.mode=[1 for j in range(testgenenum) ]
  def __getitem__(self,index):
    if self.trmd[index]==0:
      input_b=traininput[index]
      promvec=trainprom[index]
      atacvec=trainatac[index]
      rnavec=trainrna[index]
      posvec=trainpos[index]
    else:
      input_b=testinput[index]
      promvec=testprom[index]
      atacvec=testatac[index]
      rnavec=testrna[index]
      posvec=testpos[index]
    input=input_b.clone().detach()
    input[:,:,0:78]=0.25
    input[:,:,1422:1500]=0.25
    input=input[:,:,78:1422]
    input_o=torch.zeros((10,input.shape[0],input.shape[1],input.shape[2]))
    shlist=list(range(1422-78))
    for rs in range(10):
      random.shuffle(shlist)
      input_o[rs,:,:,:]=input[:,:,shlist].clone().detach()
    return input,input_o,promvec,atacvec,rnavec,posvec,self.genes[index]
  def __len__(self):
    return len(self.genes)



class CustomDataset_b(torch.utils.data.Dataset):
  def __init__(self,train=True):
    self.atac=[]
    self.seqe=[]
    if train==True:
      self.atac=[atac_b_cal(j) for j in trainpeaklist ]
      self.seqe=[peak_b_cal(j) for j in trainpeaklist ]
      self.mode=[0 for j in trainpeaklist ]
    else:
      self.atac=[atac_b_cal(j) for j in testpeaklist ]
      self.seqe=[peak_b_cal(j) for j in testpeaklist ]
      self.mode=[1 for j in testpeaklist ]
  def __getitem__(self,index):
    input_b=self.seqe[index]
    input=input_b
    input[:,0:78]=0.25
    input[:,1422:1500]=0.25
    if self.mode[index]==0:
      if torch.rand(1)<flipr:
        idx=torch.randint(low=75, high=82, size=(1,))
        input=input[:,idx:(idx+1344)]
      else:
        input=input[:,78:1422]
    else:
        input=input[:,78:1422]
    if self.mode[index]==0:
      if torch.rand(1)<revr:
        input=torch.flip(input,dims=[1])
        input_o=input.clone().detach()
        input_o[0,:]=input[3,:]
        input_o[1,:]=input[0,:]
        input_o[2,:]=input[1,:]
        input_o[3,:]=input[2,:]
      else:
        input_o=input.clone().detach()
    else:
        input_o=input.clone().detach()
    atacvec=self.atac[index]
    return input_o,input_o,atacvec
  def __len__(self):
    return len(self.atac)


b_train_dataset=CustomDataset_b(train=True)
b_test_dataset=CustomDataset_b(train=False)

train_dataset=CustomDataset(train=True)
test_dataset=CustomDataset(train=False)

dt_now = datetime.datetime.now()
print(dt_now)
os.system("free -h")

train_batch = DataLoader(
    dataset=train_dataset,  # データセットの指定
    batch_size=batchsize,  # バッチサイズの指定
    shuffle=True,  # シャッフルするかどうかの指定
    num_workers=2,
    pin_memory=True)  # コアの数

test_batch = DataLoader(
    dataset=test_dataset,
    batch_size=batchsize,
    shuffle=False,
    num_workers=2,
    pin_memory=True)

b_train_batch = DataLoader(
    dataset=b_train_dataset,
    batch_size=128,
    shuffle=True,
    num_workers=2,
    pin_memory=True)

b_test_batch = DataLoader(
    dataset=b_test_dataset,
    batch_size=128,
    shuffle=False,
    num_workers=2,
    pin_memory=True)

# ニューラルネットワークの定義

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
#device = torch.device("cpu")

class PeakembedModel(nn.Module):
    def __init__(self):
        super().__init__()
        self.seqembed= nn.Sequential(
            nn.GELU(),
            nn.Conv1d(4, 288, kernel_size=17, stride=1, bias=False, padding=8),
            nn.BatchNorm1d(288,momentum=0.01,eps=1e-03),
            nn.MaxPool1d(3, stride=3),
            nn.GELU(),
            nn.Conv1d(288, 288, kernel_size=5, stride=1, bias=False, padding=2),
            nn.BatchNorm1d(288,momentum=0.01,eps=1e-03),
            nn.MaxPool1d(2, stride=2),
            nn.GELU(),
            nn.Conv1d(288, 323, kernel_size=5, stride=1, bias=False, padding=2),
            nn.BatchNorm1d(323,momentum=0.01,eps=1e-03),
            nn.MaxPool1d(2, stride=2),
            nn.GELU(),
            nn.Conv1d(323, 363, kernel_size=5, stride=1, bias=False, padding=2),
            nn.BatchNorm1d(363,momentum=0.01,eps=1e-03),
            nn.MaxPool1d(2, stride=2),
            nn.GELU(),
            nn.Conv1d(363, 407, kernel_size=5, stride=1, bias=False, padding=2),
            nn.BatchNorm1d(407,momentum=0.01,eps=1e-03),
            nn.MaxPool1d(2, stride=2),
            nn.GELU(),
            nn.Conv1d(407, 456, kernel_size=5, stride=1, bias=False, padding=2),
            nn.BatchNorm1d(456,momentum=0.01,eps=1e-03),
            nn.MaxPool1d(2, stride=2),
            nn.GELU(),
            nn.Conv1d(456, 512, kernel_size=5, stride=1, bias=False, padding=2),
            nn.BatchNorm1d(512,momentum=0.01,eps=1e-03),
            nn.MaxPool1d(2, stride=2),
            nn.GELU(),
            nn.Conv1d(512, 256, kernel_size=1, stride=1, bias=False),
            nn.BatchNorm1d(256,momentum=0.01,eps=1e-03),
            nn.GELU(),
          )
        torch.nn.init.kaiming_normal_(self.seqembed[1].weight, mode="fan_in", nonlinearity="relu")
        torch.nn.init.kaiming_normal_(self.seqembed[5].weight, mode="fan_in", nonlinearity="relu")
        torch.nn.init.kaiming_normal_(self.seqembed[9].weight, mode="fan_in", nonlinearity="relu")
        torch.nn.init.kaiming_normal_(self.seqembed[13].weight, mode="fan_in", nonlinearity="relu")
        torch.nn.init.kaiming_normal_(self.seqembed[17].weight, mode="fan_in", nonlinearity="relu")
        torch.nn.init.kaiming_normal_(self.seqembed[21].weight, mode="fan_in", nonlinearity="relu")
        torch.nn.init.kaiming_normal_(self.seqembed[25].weight, mode="fan_in", nonlinearity="relu")
        torch.nn.init.kaiming_normal_(self.seqembed[29].weight, mode="fan_in", nonlinearity="relu")
        self.seqfc = nn.Sequential(
            nn.Linear(256*7, 32,bias=False),
            nn.BatchNorm1d(32,momentum=0.01,eps=1e-03),
            nn.Dropout(0.2),
            nn.GELU(),
          )
        torch.nn.init.kaiming_normal_(self.seqfc[0].weight, mode="fan_in", nonlinearity="relu")
        self.atacpred = nn.Sequential(
          nn.Linear(32, fullcell),
          nn.Sigmoid(),
          )
        torch.nn.init.kaiming_normal_(self.atacpred[0].weight, mode="fan_in", nonlinearity="relu")
    def forward(self, src, prom, mode):
        #src:B,L,4,bp prom:B,4,bp atac:B,C,L
        if mode==1:
          src=src.view(-1,mode,4,basenuminput) #pre:B,1,4,bp tr:B,L,4,bp
          src_n=src.view(-1,4,basenuminput)
          embseq=self.seqembed(src_n)
          embseq=embseq.view(-1,256*7)
          embseq=self.seqfc(embseq)
          embseq=embseq.view(-1,mode,32) #B,L,dim
          embprom=embseq[:,0,:]
          mask=src.sum(dim=(2,3)) #B,L
          mask = (mask != PAD_IDX) #B,L
          atacpredmtx=self.atacpred(embseq) #B,L,dim to B,L,C
          atacpredmtx=atacpredmtx.transpose(1,2) #B,C,L
          atacpredmtx=atacpredmtx.transpose(0,1) #C,B,L
          atacpredmtx_o=torch.mul(atacpredmtx,mask).transpose(0,1) #B,C,L
        else:
          batchsize=src.shape[0]
          src=src.view(-1,mode,4,basenuminput) #pre:B,1,4,bp tr:B,L,4,bp
          src_n=src.view(-1,4,basenuminput)
          src_ns=src_n[src_n.sum(dim=(1,2))!=0,:,:]
          embseq_in=self.seqembed(src_ns)
          embseq_in=embseq_in.view(-1,256*7)
          embseq_out=self.seqfc(embseq_in)
          embseq=torch.zeros((batchsize*mode,32)).to(device)
          embseq[src_n.sum(dim=(1,2))!=0,:]=embseq_out
          embseq=embseq.view(-1,mode,32) #B,L,dim
          embprom=embseq[:,0,:]
          mask=src.sum(dim=(2,3)) #B,L
          mask = (mask != PAD_IDX) #B,L
          atacpredmtx=self.atacpred(embseq) #B,L,dim to B,L,C
          atacpredmtx=atacpredmtx.transpose(1,2) #B,C,L
          atacpredmtx=atacpredmtx.transpose(0,1) #C,B,L
          atacpredmtx_o=torch.mul(atacpredmtx,mask).transpose(0,1) #B,C,L
        return embprom, embseq, mask, atacpredmtx_o



class TransformerModel(nn.Module):
    def __init__(self):
        super().__init__()
        self.prom_q = nn.Linear(32, 32*fdim)
        self.prom_sig = nn.Linear(32, 1)
        self.enha_q = nn.Linear(33, 32*fdim)
        self.layer_norm = nn.LayerNorm([max_len-1],elementwise_affine=False)
        self.encoder = nn.Sequential(
          nn.Linear(fdim, 128),
          nn.ReLU(),
          nn.Dropout(0.25),
          nn.Linear(128, 64),
          )
    def forward(self, src, prom, atac, mask):
        #src:B,L,D prom:B,D atac:B,C,L
        prom_qm=self.prom_q(prom).view(-1,fdim,32)
        prom_qm=prom_qm.unsqueeze(2) #B,MH,L(1),D
        enha_km=self.enha_q(src).view(-1,max_len,32,fdim)
        enha_km=enha_km.transpose(1,3) #B,MH,D,L
        QK=torch.matmul(prom_qm, enha_km).squeeze(dim=2) #B,MH,L
        QK=QK.transpose(0,1) #MH,B,L
        QK=torch.mul(QK,mask).transpose(0,1) #B,MH,L
        QKs=QK[:,:,1:max_len]
        QKs=self.layer_norm(QKs)
        prom_qms=self.prom_sig(prom_qm) #B,MH,L(1),1
        prom_qms=prom_qms.squeeze(dim=3) #B,MH,1
        atac=atac.transpose(1,2) #B,L,C
        atac_s=atac[:,1:max_len,:]
        atac_p=atac[:,0,:] #B,C
        atac_p=atac_p.unsqueeze(1) #B,1,C
        QKs=torch.matmul(QKs, atac_s) #B,MH,C
        QKp=torch.matmul(prom_qms, atac_p) #B,MH,C
        QKall=torch.add(QKs,QKp) #B,MH,C
        QKall=QKall.transpose(1,2) #B,C,MH
        QKall=self.encoder(QKall)
        return QKall



class RNAModel(nn.Module):
    def __init__(self):
        super().__init__()
        self.decoder = nn.Sequential(
          nn.Linear(64, 5120),
          nn.ReLU(),
          nn.Dropout(0.25),
          nn.Linear(5120, 50),
          )
    def forward(self, src):
        output2_e = self.decoder(src)
        return output2_e




class FullModel(nn.Module):
    def __init__(self,pos, genet):
        super().__init__()
        self.seqembed= nn.Sequential(
            nn.GELU(),
            nn.Conv1d(4, 288, kernel_size=17, stride=1, bias=False, padding=8),
            nn.BatchNorm1d(288,momentum=0.01,eps=1e-03),
            nn.MaxPool1d(3, stride=3),
            nn.GELU(),
            nn.Conv1d(288, 288, kernel_size=5, stride=1, bias=False, padding=2),
            nn.BatchNorm1d(288,momentum=0.01,eps=1e-03),
            nn.MaxPool1d(2, stride=2),
            nn.GELU(),
            nn.Conv1d(288, 323, kernel_size=5, stride=1, bias=False, padding=2),
            nn.BatchNorm1d(323,momentum=0.01,eps=1e-03),
            nn.MaxPool1d(2, stride=2),
            nn.GELU(),
            nn.Conv1d(323, 363, kernel_size=5, stride=1, bias=False, padding=2),
            nn.BatchNorm1d(363,momentum=0.01,eps=1e-03),
            nn.MaxPool1d(2, stride=2),
            nn.GELU(),
            nn.Conv1d(363, 407, kernel_size=5, stride=1, bias=False, padding=2),
            nn.BatchNorm1d(407,momentum=0.01,eps=1e-03),
            nn.MaxPool1d(2, stride=2),
            nn.GELU(),
            nn.Conv1d(407, 456, kernel_size=5, stride=1, bias=False, padding=2),
            nn.BatchNorm1d(456,momentum=0.01,eps=1e-03),
            nn.MaxPool1d(2, stride=2),
            nn.GELU(),
            nn.Conv1d(456, 512, kernel_size=5, stride=1, bias=False, padding=2),
            nn.BatchNorm1d(512,momentum=0.01,eps=1e-03),
            nn.MaxPool1d(2, stride=2),
            nn.GELU(),
            nn.Conv1d(512, 256, kernel_size=1, stride=1, bias=False),
            nn.BatchNorm1d(256,momentum=0.01,eps=1e-03),
            nn.GELU(),
          )
        self.seqfc = nn.Sequential(
            nn.Linear(256*7, 32,bias=False),
            nn.BatchNorm1d(32,momentum=0.01,eps=1e-03),
            nn.Dropout(0.2),
            nn.GELU(),
          )
        self.atacpred = nn.Sequential(
          nn.Linear(32, fullcell),
          nn.Sigmoid(),
          )
        self.makeQ = nn.Linear(33, 33*edim)
        self.makeK = nn.Linear(33, 33*edim)
        self.layer_norm = nn.LayerNorm([max_len-1,edim])
        self.softm = nn.Softmax(dim=3)
        self.decoder = nn.Sequential(
          nn.Linear(edim, 512),
          nn.BatchNorm1d(512),
          nn.ReLU(),
          nn.Dropout(0.25),
          nn.Linear(512, 1),
          )
        self.prom_q = nn.Linear(32, 32*fdim)
        self.prom_sig = nn.Linear(32, 1)
        self.enha_q = nn.Linear(33, 32*fdim)
        self.layer_norm2 = nn.LayerNorm([max_len-1],elementwise_affine=False)
        self.encoder = nn.Sequential(
          nn.Linear(fdim, 128),
          nn.ReLU(),
          nn.Dropout(0.25),
          nn.Linear(128, 64),
          )
        self.decoder2 = nn.Sequential(
          nn.Linear(64, 5120),
          nn.ReLU(),
          nn.Dropout(0.25),
          nn.Linear(5120, 50),
          )
        self.pos=pos
        self.genet=genet
        #self.atac=atac
    def forward(self, src):
        #src:B,L,4,bp prom:B,4,bp atac:B,C,L
        #prom, pos, genet
        mode=max_len
        pos=self.pos
        genet=self.genet
        fbatchsize=src.shape[0]
        src=src.view(-1,mode,4,basenuminput) #pre:B,1,4,bp tr:B,L,4,bp
        src_n=src.reshape(src.shape[0]*src.shape[1],4,basenuminput)
        src_ns=src_n[src_n.sum(dim=(1,2))!=0,:,:]
        embseq_in=self.seqembed(src_ns)
        embseq_in=embseq_in.view(-1,256*7)
        embseq_out=self.seqfc(embseq_in)
        embseq=torch.zeros((fbatchsize*mode,32)).to(device)
        embseq[src_n.sum(dim=(1,2))!=0,:]=embseq_out
        embseq=embseq.view(-1,mode,32) #B,L,dim
        embprom=embseq[:,0,:]
        mask=src.sum(dim=(2,3)) #B,L
        mask = (mask != PAD_IDX) #B,L
        atacpredmtx=self.atacpred(embseq) #B,L,dim to B,L,C
        atacpredmtx=atacpredmtx.transpose(1,2) #B,C,L
        atacpredmtx=atacpredmtx.transpose(0,1) #C,B,L
        atacpredmtx_o=torch.mul(atacpredmtx,mask).transpose(0,1) #B,C,L
        inputvec=torch.zeros((embseq.shape[0],embseq.shape[1],embseq.shape[2]+1))
        inputvec[:,:,0:32]=embseq
        inputvec[:,:,32]=pos
        inputvec=inputvec.to(device)
        #def forward(self, src, prom, atac, mask):
        #src:B,L,D prom:B,D atac:B,C,L
        prom_qm=self.prom_q(embprom).view(-1,fdim,32)
        prom_qm=prom_qm.unsqueeze(2) #B,MH,L(1),D
        enha_km=self.enha_q(inputvec).view(-1,max_len,32,fdim)
        enha_km=enha_km.transpose(1,3) #B,MH,D,L
        QK2=torch.matmul(prom_qm, enha_km).squeeze(dim=2) #B,MH,L
        QK2=QK2.transpose(0,1) #MH,B,L
        QK2=torch.mul(QK2,mask).transpose(0,1) #B,MH,L
        QKs2=QK2[:,:,1:max_len]
        QKs2=self.layer_norm2(QKs2)
        prom_qms=self.prom_sig(prom_qm) #B,MH,L(1),1
        prom_qms=prom_qms.squeeze(dim=3) #B,MH,1
        atac_o=atacpredmtx_o.transpose(1,2) #B,L,C
        atac_s2=atac_o[:,1:max_len,:]
        atac_p=atac_o[:,0,:] #B,C
        atac_p=atac_p.unsqueeze(1) #B,1,C
        QKs2=torch.matmul(QKs2, atac_s2) #B,MH,C
        QKp2=torch.matmul(prom_qms, atac_p) #B,MH,C
        QKall2=torch.add(QKs2,QKp2) #B,MH,C
        QKall2=QKall2.transpose(1,2) #B,C,MH
        QKall2=self.encoder(QKall2)
        #return QKall2
        #def forward(self, src):
        output2_e = self.decoder2(QKall2)
        #return output2_e
        output2_e=output2_e.view(-1,50)
        pred_o = torch.matmul(output2_e, Vs)
        pred_o=pred_o.view(-1,testcell,fullgenenum)
        pred_o=pred_o[range(pred_o.shape[0]),:,genet] #B,C
        pred_os=pred_o[:,cellclus_test==celltype]
        usecellmtx_s=usecellmtx[:,cellclus_test==celltype]
        pred_oss=pred_os*usecellmtx_s[genet,:]
        pred_om=pred_oss.sum(axis=1)
        pred_om=pred_om.unsqueeze(1)
        return pred_om

print(torch.cuda.memory_allocated())

transformernet = TransformerModel()
rnanet = RNAModel()
peakembed = PeakembedModel()

peakembed.load_state_dict(torch.load(samplename+"/DNA_embed_fold0.pth"))
transformernet.load_state_dict(torch.load(samplename+"/transformer_rnapred_transformer_t2_raw_2_fold0.pth"))
rnanet.load_state_dict(torch.load(samplename+"/transformer_rnapred_rna_t2_raw_2_fold0.pth"))
Vs=torch.load(samplename+"/pca_50dim_matrix_raw.pt").detach().clone().to(device)

usecellmtx=torch.from_numpy(usecellmtx).detach().clone().to(device)

# デバイスの確認
print("Device: {}".format(device))

#criterion = losscal()  # 損失関数（平均二乗誤差: MSE）
criterion = nn.MSELoss()
kl_loss = nn.KLDivLoss(reduction="batchmean", log_target=True)

BCEL = torch.nn.BCELoss()
BCEL_each = torch.nn.BCELoss(reduction="none")

smx = nn.Softmax(dim=2)

def losscal(pred_e,pred_u,rna):
  pred_e_f=pred_e.view(-1,2)
  rna=rna.flatten()
  pred_e=smx(pred_e)
  pred_u=smx(pred_u)
  pred_u=pred_u.unsqueeze(dim=3)
  pred_e=pred_e.unsqueeze(dim=2)
  c=torch.matmul(pred_u, pred_e).sum(dim=1)
  d=c[:,:,0].log_softmax(dim=1)
  e=c[:,:,1].log_softmax(dim=1)
  loss=criterion(pred_e_f,rna)+lamd*kl_loss(d,e)
  return loss


def totallosscal(pred_e,rna,qkm):
  loss=criterion(pred_e,rna)+llam*torch.abs(qkm).sum()
  return loss


def rnalosscal(pred_e,rna):
  loss=criterion(pred_e,rna)
  return loss


def ataclosscal(atacpred,atac):
  loss=BCEL(atacpred,atac)
  #loss=msel2(atacpred,atac)
  return loss


def atac_totallosscal(atacpred,atac,pred_e,rna,qkm,peaktt):
  aloss=BCEL_each(atacpred,atac) #B,C,L
  aloss_s=(torch.mul(aloss.transpose(0,1),peaktt)).transpose(0,1)
  rloss=criterion(pred_e,rna)+llam*torch.abs(qkm).sum()
  loss=rloss+alam*(aloss_s.mean())
  return loss

# 最適化関数の定義


# 損失を保存するリストを作成
train_loss_list = []  # 学習損失
test_loss_list = []  # 評価損失

train_loss_list_rna = []  # 学習損失
test_loss_list_rna = []  # 評価損失

dt_now = datetime.datetime.now()
print(dt_now)

attr_list = []
input_list = []


niter=20
tmpi=0

for data, rdata, prom, atac, rna, pos, genet in test_batch:
    print(tmpi)
    tmpi=tmpi+1
    tmpenha=(enhaclus[genet,:]).astype(int)
    if (tmpenha==1).sum()>0:
      tmp=atac.sum(dim=1)
      niter=20
      if ((tmp!=0).sum(dim=1)>50).sum()>0:
        niter=10
      # GPUにTensorを転送
      data = data.to(device)
      pos = pos.to(device)
      genet = genet.to(device)
      data=data.detach().requires_grad_(True)
      transformernet = TransformerModel()
      rnanet = RNAModel()
      peakembed = PeakembedModel()
      peakembed.load_state_dict(torch.load(samplename+"/DNA_embed_fold0.pth"))
      transformernet.load_state_dict(torch.load(samplename+"/transformer_rnapred_transformer_t2_raw_2_fold0.pth"))
      rnanet.load_state_dict(torch.load(samplename+"/transformer_rnapred_rna_t2_raw_2_fold0.pth"))
      Vs=torch.load(samplename+"/pca_50dim_matrix_raw.pt").detach().clone().to(device)
      transformernet = transformernet.to(device)
      rnanet = rnanet.to(device)
      peakembed = peakembed.to(device)
      fullnet=FullModel(pos,genet)
      fullnet.seqembed = peakembed.seqembed
      fullnet.seqfc = peakembed.seqfc
      fullnet.atacpred = peakembed.atacpred
      fullnet.prom_q = transformernet.prom_q
      fullnet.prom_sig = transformernet.prom_sig
      fullnet.enha_q = transformernet.enha_q
      fullnet.layer_norm2 = transformernet.layer_norm
      fullnet.encoder = transformernet.encoder
      fullnet.decoder2 = rnanet.decoder
      del peakembed
      del transformernet
      del rnanet
      fullnet.to(device)
      fullnet.eval()
      # データを入力して予測値を計算（順伝播）
      #for rs in range(10):
        #print(data.shape)
        #print(rdata.shape)
      rdata=torch.zeros((1,data.shape[1],data.shape[2],data.shape[3]))
      explainer = PathExplainerTorch(fullnet)
      interactions = explainer.attributions(input_tensor = data,
                                      baseline=rdata,
                                      num_samples=niter,
                                      use_expectation=False)
      if interactions.shape[0]==1:
        tmpenha=tmpenha[np.newaxis,:]
      interactions_out=interactions.to('cpu').detach().numpy().copy()
      data_out=data.to('cpu').detach().numpy().copy()
      interactions_out=interactions_out[tmpenha==1,:,:]
      data_out=data_out[tmpenha==1,:,:]
      attr_list.append(interactions_out)
      input_list.append(data_out)
      del data_out
      del interactions_out
      del tmpenha
      del fullnet
      del data
      del atac
      del interactions
      del explainer
      del rdata
      del Vs
      torch.cuda.empty_cache()
      print(torch.cuda.memory_allocated())



inputm_full_s=np.concatenate(input_list, 0) #B,RS,L,4,1344
fname = samplename+"/"+prefix+"inputmtx.npy"
np.save(fname, inputm_full_s)

featurem_full_s=np.concatenate(attr_list,0) 
fname = samplename+"/"+prefix+"contribution.npy"
np.save(fname, featurem_full_s)