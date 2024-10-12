import scanpy as sc
import pandas as pd
import numpy as np
import sys

args = sys.argv
samplename=str(args[1])


fname=samplename+"/Deeplift_full_ver2_all.npy"
enhamtx=np.load(fname)
grad=np.load(samplename+"/allgrad_ssep_max.npy")

fullgenename=pd.read_csv(samplename+"/pair_300000.csv",header=None,names=("gene","start","end","pos"))
fullgenename=fullgenename["gene"]
fullgenename=fullgenename.to_numpy()

pairlist=pd.read_csv(samplename+"/pair_300000.csv",header=None,names=("gene","start","end","pos"))
pairlist=pairlist[["start","end","pos"]].to_numpy()

pairlist_prom=pd.read_csv(samplename+"/pair_promoter.csv",header=None)
pairlist_prom=pairlist_prom[[1,2,3]].to_numpy()

ATAC_use=np.load(samplename+"/ATAC_pred_fold0.npy")
ATAC_use=ATAC_use[:,:,0]

ATAC_corr=np.zeros(grad.shape)
for j in range(ATAC_corr.shape[0]):
  print(j)
  peakstart=int(pairlist[j,0])
  peakend=int(pairlist[j,1])
  peaknum_gene=peakend-peakstart+1
  prompos=int(pairlist_prom[j,0])
  enha=list(range(peakstart,(peakend+1)))
  enha.remove(prompos)
  ATAC_tmp=np.zeros((enhamtx.shape[1],enhamtx.shape[2]))
  ATAC_tmp[:,0]=ATAC_use[prompos,:]
  ATAC_tmp[:,1:peaknum_gene]=ATAC_use[enha,:].transpose(1,0)
  for i in range(peaknum_gene):
    ATAC_corr[j,i]=np.corrcoef([ATAC_tmp[:,i],enhamtx[j,:,i]])[0,1]


thv=np.percentile(grad[grad!=0], [40])
sthv=np.percentile(grad[grad!=0], [85])

promenhatag=np.zeros(grad.shape)
promenhatag[(grad>thv)&(ATAC_corr>0)]=2
promenhatag[(grad>sthv)&(ATAC_corr>0)]=4
promenhatag[:,0]=3



enhathr=np.percentile(enhamtx.transpose(0,2,1),75,axis=2) #B,L
enhathrmtx=(enhamtx.transpose(1,0,2)>enhathr).transpose(1,0,2) #B,C,L

enhausemtx=((enhathrmtx.transpose(1,0,2)*((promenhatag==4)|(promenhatag==2))).sum(axis=2))/((promenhatag==4)|(promenhatag==2)).sum(axis=1)

enhausetagmtx=enhausemtx.transpose(1,0)>0.7

np.save(samplename+'/enhausetagmtx.npy',enhausetagmtx)
np.save(samplename+'/promenhatag.npy',promenhatag)


################



RNAembed=pd.read_csv(samplename+"/pca_50dim_rnaembed_raw.txt",sep=" ",header=None)
RNAembed=RNAembed.to_numpy()
fullcell=RNAembed.shape[0]

pd.DataFrame(np.ones((fullcell))).to_csv(samplename+"/celltag1.txt",sep="\t",header=False, index=False)

###################

