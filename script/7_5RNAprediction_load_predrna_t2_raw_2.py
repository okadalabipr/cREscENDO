import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch
from torch.utils.data import TensorDataset, DataLoader
from torch import nn
import torch.nn.functional as F
from torch import optim
import os
import cv2
import datetime
import pickle
from tqdm import tqdm
from scipy import stats
from sync_batchnorm import convert_model, DataParallelWithCallback
import sys

args = sys.argv
samplename=str(args[1])
chunkid=int(float(args[2]))

#celllist=pd.read_csv(samplename+"/celllist_scVI_afterjoin.csv")
#celllist=celllist.to_numpy()
#celllist_lis=celllist[:,3]-1
#usecell=celllist_lis.tolist()

peaklist=pd.read_csv(samplename+"/peaks.bed",sep="\t",header=None)
peaklist=peaklist[[1,2]].to_numpy()
peaklist=torch.from_numpy(peaklist)
peaklist=peaklist.to(torch.float32)

#sig = nn.Sigmoid()
ATACmatrix=np.load(samplename+"/ATAC_pred_fold"+str(chunkid)+".npy")
ATACmatrix=ATACmatrix[:,:,0]
ATACmatrix=torch.from_numpy(ATACmatrix)
ATACmatrix=ATACmatrix.to(torch.float32)
#ATACmatrix=ATACmatrix[:,usecell]
#ATACmatrix=sig(ATACmatrix)
ATACmatrix=nn.functional.normalize(ATACmatrix,p=2.0, dim=1, eps=1e-12, out=None)

print(ATACmatrix.shape)

RNAmatrix=pd.read_csv(samplename+"/log1praw.csv",sep=",",header=None)
RNAmatrix=RNAmatrix.to_numpy()
RNAmatrix=torch.from_numpy(RNAmatrix)
RNAmatrix=RNAmatrix.to(torch.float32)
#RNAmatrix=RNAmatrix+1
#RNAmatrix=torch.log(RNAmatrix)
RNAmatrix_b=((RNAmatrix.transpose(0,1)-RNAmatrix.mean(dim=1))/RNAmatrix.std(dim=1)).transpose(0,1)

RNAembed=pd.read_csv(samplename+"/pca_50dim_rnaembed_raw.txt",sep=" ",header=None)
RNAembed=RNAembed.to_numpy()
RNAembed=torch.from_numpy(RNAembed)
RNAembed=RNAembed.to(torch.float32)

print(RNAembed.shape)

pairlist_prom=pd.read_csv(samplename+"/pair_promoter.csv",header=None)
pairlist_prom=pairlist_prom[[1,2,3]].to_numpy()
pairlist_prom=torch.from_numpy(pairlist_prom)
pairlist_prom=pairlist_prom.to(torch.float32)

pairlist=pd.read_csv(samplename+"/pair_300000.csv",header=None,names=("gene","start","end","pos"))
pairlist_train=pd.read_csv(samplename+"/pair_300000_train_fold"+str(chunkid)+".csv",header=None,names=("gene","start","end","pos"))
pairlist_test=pd.read_csv(samplename+"/pair_300000_test_fold"+str(chunkid)+".csv",header=None,names=("gene","start","end","pos"))

traingenelist=pairlist["gene"].isin(pairlist_train["gene"]).values.tolist()
testgenelist=pairlist["gene"].isin(pairlist["gene"]).values.tolist()

pairlist=pairlist[["start","end","pos"]].to_numpy()
pairlist=torch.from_numpy(pairlist)
pairlist=pairlist.to(torch.float32)

embedlist=pd.read_csv(samplename+"/embed_fold"+str(chunkid)+".txt",sep=" ",header=None)
embedlist=embedlist.to_numpy()
embedlist=torch.from_numpy(embedlist)
embedlist=embedlist.to(torch.float32)

traincellnum=np.load(samplename+"/traincellnum.npy")

fullcell=RNAembed.shape[0]
traincell=traincellnum #10000 11830 11740
#testcell=fullcell-traincell #1893
testcell=RNAembed.shape[0] #1893
traingenenum=sum(traingenelist) #3599 #11881
testgenenum=sum(testgenelist)
#fullgenenum=traingenenum+testgenenum
fullgenenum=testgenenum
traintag=0
testtag=0

max_len=int((pairlist[:,1]-pairlist[:,0]).max()+1)
print(max_len)

PAD_IDX=0
epoch=100
fdim=512
edim=32
lamd=0
llam=0.001

traingenelist=torch.from_numpy(np.array(range(fullgenenum))[traingenelist])
testgenelist=torch.from_numpy(np.array(range(fullgenenum))[testgenelist])

os.system("free -h")
torch.backends.cudnn.benchmark = True

def dataprep(train=True):
  inputmtx=[]
  prommtx=[]
  atacmtx=[]
  rnamtx=[]

  if train==True:
    inputmtx = [inputcal(traingenelist[j]) for j in range(traingenenum) ]
    prommtx = [promcal(traingenelist[j]) for j in range(traingenenum) ]
    atacmtx = [ataccal(traingenelist[j],traintag,traincell+traintag) for j in range(traingenenum) ]
    rnamtx= [rnacal(traingenelist[j],traintag,traincell+traintag) for j in range(traingenenum) ]
  else:
    inputmtx = [inputcal(testgenelist[j]) for j in range(testgenenum) ]
    prommtx = [promcal(testgenelist[j]) for j in range(testgenenum) ]
    atacmtx = [ataccal(testgenelist[j],testtag,testcell+testtag) for j in range(testgenenum) ]
    rnamtx= [rnacal(testgenelist[j],testtag,testcell+testtag) for j in range(testgenenum) ]
  return inputmtx,prommtx,atacmtx,rnamtx

def inputcal(j):
  peakstart=int(pairlist[j,0].item())
  peakend=int(pairlist[j,1].item())
  prompos=int(pairlist_prom[j,0].item())
  enha=list(range(peakstart,(peakend+1)))
  enha.remove(prompos)

  genepos=int(pairlist[j,2].item())
  peaknum_gene=peakend-peakstart+1
  input = torch.zeros((max_len,33))

  embedvec=embedlist[enha]
  posvec=((((peaklist[enha,0]+peaklist[enha,1])/2)-genepos)/300000)

  input[0,0:32]=embedlist[prompos]
  input[0,32]=((((peaklist[prompos,0]+peaklist[prompos,1])/2)-genepos)/300000)

  input[1:peaknum_gene,0:32]=embedvec
  input[1:peaknum_gene,32]=posvec
  return input

def promcal(j):
  peakstart=int(pairlist_prom[j,0].item())
  embedvec=embedlist[peakstart,:]
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

  atacvec=ATACmatrix[enha,cellstart:cellend]
  atacvec = atacvec.transpose(1,0)
  atacvec_return[:,0] = ATACmatrix[prompos,cellstart:cellend]
  atacvec_return[:,1:peaknum_gene]=atacvec
  return atacvec_return

def rnacal(j,cellstart,cellend):
  #rnavec=RNAembed[cellstart:cellend,:]
  rnavec=RNAmatrix_b[j,cellstart:cellend]
  #rnavec = rnavec.unsqueeze(-1)
  return rnavec

traininput,trainprom,trainatac,trainrna=dataprep(True)
testinput,testprom,testatac,testrna=dataprep(False)

class CustomDataset(torch.utils.data.Dataset):
  def __init__(self,train=True):
    self.genes=[]
    self.trmd=[]

    if train==True:
      self.genes=[traingenelist[j] for j in range(traingenenum)]
      self.trmd=[0 for j in range(traingenenum)]
    else:
      self.genes=[testgenelist[j] for j in range(testgenenum)]
      self.trmd=[1 for j in range(testgenenum)]
  def __getitem__(self,index):
    if self.trmd[index]==0:
      input=traininput[index]
      promvec=trainprom[index]
      atacvec=trainatac[index]
      rnavec=trainrna[index]
    else:
      input=testinput[index]
      promvec=testprom[index]
      atacvec=testatac[index]
      rnavec=testrna[index]

    return input,promvec,atacvec,rnavec,self.genes[index]
  def __len__(self):
    return len(self.genes)



full_dataset=CustomDataset(train=False)

dt_now = datetime.datetime.now()
print(dt_now)
os.system("free -h")

full_batch = DataLoader(
    dataset=full_dataset,
    batch_size=8,
    shuffle=False,
    num_workers=2,
    pin_memory=True)

# ニューラルネットワークの定義



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
    def forward(self, src, prom, atac):
        #src:B,L,D prom:B,D atac:B,C,L
        mask=src.transpose(1,2) #B,D,L
        mask=mask.transpose(0,1) # D,B,L
        #src_padding_mask_m = (mask[0,:,:] != PAD_IDX) #B,L
        mask = (mask[0,:,:] != PAD_IDX) #B,L
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


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

transformernet = TransformerModel()
rnanet = RNAModel()

transformernet.load_state_dict(torch.load(samplename+"/transformer_rnapred_transformer_t2_raw_2_fold"+str(chunkid)+".pth"))
rnanet.load_state_dict(torch.load(samplename+"/transformer_rnapred_rna_t2_raw_2_fold"+str(chunkid)+".pth"))

transformernet = convert_model(transformernet).to(device)
transformernet = DataParallelWithCallback(transformernet)

rnanet = convert_model(rnanet).to(device)
rnanet = DataParallelWithCallback(rnanet)

Vs=torch.load(samplename+"/pca_50dim_matrix_raw.pt").detach().clone().to(device)

# デバイスの確認
print("Device: {}".format(device))

#criterion = losscal()  # 損失関数（平均二乗誤差: MSE）
criterion = nn.MSELoss()
kl_loss = nn.KLDivLoss(reduction="batchmean", log_target=True)

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

def rnalosscal_grad(pred_e,rna):
  pred_e=pred_e.flatten()
  rna=rna.flatten()
  loss=(pred_e-rna).sum()
  return loss

# 最適化関数の定義
optimizer = optim.Adam(transformernet.parameters(),lr=1e-5,betas=(0.95, 0.9995),eps=1e-07)
rnaoptimizer = optim.Adam(rnanet.parameters(),lr=1e-3,betas=(0.95, 0.9995),eps=1e-07)


# 損失を保存するリストを作成
train_loss_list = []  # 学習損失
test_loss_list = []  # 評価損失


train_loss_list_rna = []  # 学習損失
test_loss_list_rna = []  # 評価損失

dt_now = datetime.datetime.now()
print(dt_now)



m = nn.Softmax(dim=1)
relu = nn.ReLU()

transformernet.eval()
rnanet.eval()

pred_i_lis = []
true_x_lis = []
pred_u_lis = []
true_u_lis = []

qk_lis = []
grad_lis = []
prom_list = []
atac_list = []
rna_list = []

with torch.no_grad():
  for data, prom, atac, rna, genet in full_batch:
      # GPUにTensorを転送
      data = data.to(device)
      prom = prom.to(device)
      atac = atac.to(device)
      rna = rna.to(device)
      genet = genet.to(device)

      # データを入力して予測値を計算（順伝播）
      output=transformernet(data,prom,atac)
      pred_e=rnanet(output)

      pred_e=pred_e.view(-1,50)
      pred_o = torch.matmul(pred_e, Vs)
      pred_o=pred_o.view(-1,testcell,fullgenenum)
      pred_o=pred_o[range(pred_o.shape[0]),:,genet]

      rna_list.append(pred_o)

rna_m=torch.cat(rna_list, dim=0) #gene,cell,len
outn = rna_m.to("cpu").detach().numpy().copy()
np.save(samplename+"/predrna_t2_raw_2_fold"+str(chunkid)+".npy",outn)