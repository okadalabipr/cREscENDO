import numpy as np
import pandas as pd
import torch
from torch import nn
import sys

args = sys.argv
samplename=str(args[1])

loss = nn.MSELoss()

RNAmatrix=pd.read_csv(samplename+"/log1praw.csv",sep=",",header=None)
RNAmatrix=RNAmatrix.to_numpy()
RNAmatrix=torch.from_numpy(RNAmatrix)
RNAmatrix=RNAmatrix.to(torch.float32)
#RNAmatrix=RNAmatrix+1
#RNAmatrix=torch.log(RNAmatrix)
RNAmatrix_b=((RNAmatrix.transpose(0,1)-RNAmatrix.mean(dim=1))/RNAmatrix.std(dim=1)).transpose(0,1)
RNAmatrix_np=RNAmatrix_b.transpose(0,1) ## Cell Gene

traincellnum=np.load(samplename+"/traincellnum.npy")

print(RNAmatrix_np.shape)
U, S, V = torch.pca_lowrank(RNAmatrix_np[0:traincellnum,], q=50)
Vs = V.detach()
projection = torch.matmul(RNAmatrix_np, V)
approx = torch.matmul(projection, Vs.T)

mse = loss(approx,RNAmatrix_np)

print("Loss")
print(mse)

np.savetxt(samplename+"/pca_50dim_rnaembed_raw.txt",projection.to("cpu").detach().numpy().copy())
np.savetxt(samplename+"/RNAmatrix_raw.txt",RNAmatrix_np.to("cpu").detach().numpy().copy())
np.savetxt(samplename+"/RNAmatrix_reconst_embed_raw.txt",approx.to("cpu").detach().numpy().copy())
torch.save(Vs.T,samplename+"/pca_50dim_matrix_raw.pt") 