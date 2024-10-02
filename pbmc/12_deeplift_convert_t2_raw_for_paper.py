import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import math
import subprocess



args = sys.argv
samplename=str(args[1])


fname=samplename+"/Deeplift_full_ver2_all.npy"
enhamtx=np.load(fname)

fullgenename=pd.read_csv(samplename+"/pair_300000.csv",header=None,names=("gene","start","end","pos"))
fullgenename=fullgenename["gene"]

pairlist=pd.read_csv(samplename+"/pair_300000.csv",header=None,names=("gene","start","end","pos"))
pairlist=pairlist[["start","end","pos"]].to_numpy()

pairlist_prom=pd.read_csv(samplename+"/pair_promoter.csv",header=None)
pairlist_prom=pairlist_prom[[1,2,3]].to_numpy()

max_len=int((pairlist[:,1]-pairlist[:,0]).max()+1)

RNAembed=pd.read_csv(samplename+"/pca_50dim_rnaembed_raw.txt",sep=" ",header=None)
RNAembed=RNAembed.to_numpy()
traincellnum=np.load(samplename+"/traincellnum.npy")

fullcell=RNAembed.shape[0]
traincell=traincellnum #10000 11830 11740
testcell=fullcell-traincell #1893

cellnamelist=pd.read_csv(samplename+"/celllist_scVI.csv",header=None)
cellnamelist=cellnamelist.to_numpy()
#cellnamelist=cellnamelist[:,1]
print(cellnamelist[0])

ATACmatrix=np.load(samplename+"/ATAC_pred.npy")
#ATACmatrix=ATACmatrix[:,:,0]

ATAC_corr=np.zeros((enhamtx.shape[0],enhamtx.shape[2]))
for j in range(ATAC_corr.shape[0]):
  print(j)
  peakstart=int(pairlist[j,0])
  peakend=int(pairlist[j,1])
  peaknum_gene=peakend-peakstart+1
  prompos=int(pairlist_prom[j,0])
  enha=list(range(peakstart,(peakend+1)))
  enha.remove(prompos)
  ATAC_tmp=np.zeros((enhamtx.shape[1],enhamtx.shape[2]))
  ATAC_tmp[:,0]=ATACmatrix[prompos,:]
  ATAC_tmp[:,1:peaknum_gene]=ATACmatrix[enha,:].transpose(1,0)
  for i in range(peaknum_gene):
    ATAC_corr[j,i]=np.corrcoef([ATAC_tmp[:,i],enhamtx[j,:,i]])[0,1]


np.save(samplename+"/ATAC_corr.npy",ATAC_corr)
ATAC_corr=np.load(samplename+"/ATAC_corr.npy")


enhalist=enhamtx.max(axis=1)
enhalist[ATAC_corr<0]=0
enhalist=enhalist.flatten()

enhamtx=enhamtx.transpose(0,2,1) #B,L,C
enhamtx=enhamtx.reshape(enhamtx.shape[0]*enhamtx.shape[1],enhamtx.shape[2]) #B*L,C

enhause=np.zeros((pairlist.shape[0],max_len))
enhaname=np.zeros((pairlist.shape[0],max_len),dtype=object)
for j in range(pairlist.shape[0]):
  peakstart=int(pairlist[j,0])
  peakend=int(pairlist[j,1])
  peaknum_gene=peakend-peakstart+1
  enhause[j,0:peaknum_gene]=1
  enhaname[j,0]=fullgenename[j]+"_promoter"
  for k in range(peaknum_gene-1):
    enhaname[j,k+1]=fullgenename[j]+"_enhancer_"+str(k+1)

enhamtx_s=enhamtx[enhause.flatten()==1,:]
enhaname_s=enhaname.flatten()[enhause.flatten()==1]
enhalist_s=enhalist[enhause.flatten()==1]

plt.figure()
plt.hist(enhalist_s, bins=100, range=(0.00001, 0.3))
plt.savefig(samplename+"/hist_t2_raw.png")
plt.show()

enhamtx_se=enhamtx_s[enhalist_s>0.05,:]
enhaname_se=enhaname_s.flatten()[enhalist_s>0.05]


print(enhaname)
print(cellnamelist)
print(enhaname_se)
print(enhamtx_se.shape)

df = pd.DataFrame(enhamtx_se,columns=cellnamelist,index=enhaname_se)
dft = df.T
dft.to_csv(samplename+"/Deeplift_enhamtx_t2_raw.csv",header=True, index=True)