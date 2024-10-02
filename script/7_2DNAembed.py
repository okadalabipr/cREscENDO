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
import copy
import random
import sys

args = sys.argv
samplename=str(args[1])
chunkid=int(float(args[2]))

ATACmatrix=np.load(samplename+"/atac_count.npy")
ATACmatrix[ATACmatrix>0]=1
#ATACmatrix=ATACmatrix>0.2
ATACmatrix=torch.from_numpy(ATACmatrix)
ATACmatrix=ATACmatrix.to(torch.float32)
print(ATACmatrix.shape)

peakreadlist=np.load(samplename+"/peak_extend_read.npy",allow_pickle=True)
peakreadlist=torch.from_numpy(peakreadlist)
peakreadlist=peakreadlist.to(torch.float32) #peaks,length,4
peakreadlist=peakreadlist.transpose(1,2) #peaks,4,length

peaknum=ATACmatrix.shape[0]

c_list=np.array(range(peaknum))
splitlist=np.percentile(c_list, q=[0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]).astype(int)
peaknum_train_start=splitlist[chunkid]
peaknum_train_end=splitlist[chunkid+1]

trainpeaklist=list(range(0,peaknum_train_start))+list(range(peaknum_train_end+1,peaknum))
testpeaklist=list(range(peaknum_train_start,peaknum_train_end+1))
fullpeaklist=list(range(peaknum))

PAD_IDX=0
b_epoch=1000
fullcell=ATACmatrix.shape[1]
basenum=1500
basenuminput=1344

opt_p_b=1e-2
flipr=1
revr=0.5

os.system("free -h")
torch.backends.cudnn.benchmark = True


def peak_b_cal(j):
  input = peakreadlist[j,:,:]
  return input

def atac_b_cal(j):
  atacvec = ATACmatrix[j,:]
  return atacvec


  
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
        input_o[3,:]=input[0,:]
        input_o[2,:]=input[1,:]
        input_o[1,:]=input[2,:]
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

dt_now = datetime.datetime.now()
print(dt_now)
os.system("free -h")

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

class CustomDataset_b_all(torch.utils.data.Dataset):
  def __init__(self,train=True):
    self.atac=[]
    self.seqe=[]

    self.atac=[atac_b_cal(j) for j in fullpeaklist ]
    self.seqe=[peak_b_cal(j) for j in fullpeaklist ]
    self.mode=[1 for j in fullpeaklist ]
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



b_all_dataset=CustomDataset_b_all(train=False)

b_all_batch = DataLoader(
    dataset=b_all_dataset,
    batch_size=128,
    shuffle=False,
    num_workers=2,
    pin_memory=True)

# ニューラルネットワークの定義

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')


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

        return embprom, embseq, mask, atacpredmtx_o

peakembed = PeakembedModel()
peakembed = peakembed.to(device)

# デバイスの確認
print("Device: {}".format(device))
#criterion = losscal()  # 損失関数（平均二乗誤差: MSE）

BCEL = torch.nn.BCELoss()

def ataclosscal(atacpred,atac):
  loss=BCEL(atacpred,atac)
  return loss

# 最適化関数の定義
optimizer_b = optim.Adam(peakembed.parameters(),lr=opt_p_b,betas=(0.95, 0.9995),eps=1e-07)

# 損失を保存するリストを作成

dt_now = datetime.datetime.now()
print(dt_now)


##################################

for i in range(b_epoch):

    # エポックの進行状況を表示
    print('---------------------------------------------')
    print("Epoch: {}/{}".format(i+1, b_epoch))

    # 損失の初期化
    train_loss_atac = 0  # 学習損失
    test_loss_atac = 0  # 学習損失


    # ---------学習パート--------- #
    # ニューラルネットワークを学習モードに設定
    peakembed.train()

    # ミニバッチごとにデータをロードし学習

    for data, prom, atac in b_train_batch:

        # GPUにTensorを転送
        data = data.to(device)
        prom = prom.to(device)
        atac = atac.to(device)


        # データを入力して予測値を計算（順伝播）
        torch.cuda.empty_cache()
        _, _, _, atacpred1=peakembed(data,prom,1)
        atacpred1=atacpred1.squeeze(dim=2)

        optimizer_b.zero_grad()
        atacloss = ataclosscal(atacpred1, atac)
        atacloss.backward()
        optimizer_b.step()


        # ミニバッチごとの損失を蓄積
        train_loss_atac += atacloss.item()

    # ミニバッチの平均の損失を計算
    batch_train_loss_atac = train_loss_atac / len(b_train_batch)

    # ---------学習パートはここまで--------- #

    # ---------評価パート--------- #
    # ニューラルネットワークを評価モードに設定

    peakembed.eval()

    # 評価時の計算で自動微分機能をオフにする
    with torch.no_grad():
        for data, prom, atac in b_test_batch:
            # GPUにTensorを転送
            data = data.to(device)
            prom = prom.to(device)
            atac = atac.to(device)

            # データを入力して予測値を計算（順伝播）
            torch.cuda.empty_cache()
            _, _, _, atacpred1=peakembed(data,prom,1)
            atacpred1=atacpred1.squeeze(dim=2)

            # 損失（誤差）を計算
            atacloss = ataclosscal(atacpred1, atac)
            # ミニバッチごとの損失を蓄積
            test_loss_atac += atacloss.item()
    
    batch_test_loss_atac = test_loss_atac / len(b_test_batch)

    # ---------評価パートはここまで--------- #

    # エポックごとに損失を表示
    print("Train_Loss: {:.2E} Test_Loss: {:.2E}".format(
        batch_train_loss_atac, batch_test_loss_atac))
    
torch.save(peakembed.state_dict(), samplename+"/DNA_embed_fold"+str(chunkid)+".pth")

dt_now = datetime.datetime.now()
print(dt_now)


peakembed.eval()
##################################
with torch.no_grad():
    pred_i_lis = []
    atac_i_lis = []
    for data, prom, atac in b_all_batch:
        # GPUにTensorを転送
        data = data.to(device)
        prom = prom.to(device)
        atac = atac.to(device)

        # データを入力して予測値を計算（順伝播）
        torch.cuda.empty_cache()
        _,embseq,_,atacpred=peakembed(data,prom,1)
        pred_i_lis.append(embseq)
        atac_i_lis.append(atacpred)

featurem=torch.cat(pred_i_lis, dim=0)
featurem_o=featurem.squeeze(dim=1).to('cpu').detach().numpy().copy()
np.savetxt(samplename+"/embed_fold"+str(chunkid)+".txt",featurem_o)

featurem=torch.cat(atac_i_lis, dim=0)
featurem_o=featurem.to('cpu').detach().numpy().copy()
np.save(samplename+"/ATAC_pred_fold"+str(chunkid)+".npy",featurem_o)