import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt
import scipy.stats as stats
import sys

args = sys.argv
samplename=str(args[1])


peaklist=pd.read_csv(samplename+"/peaks.bed",sep="\t",header=None)
peaklist=peaklist[[1,2]].to_numpy()
pairlist=pd.read_csv(samplename+"/pair_300000.csv",header=None)
pairlist=pairlist[[1,2,3]].to_numpy()

pairlist_prom=pd.read_csv(samplename+"/pair_promoter.csv",header=None)
pairlist_prom=pairlist_prom[[1,2,3]].to_numpy()

max_len=int((pairlist[:,1]-pairlist[:,0]).max()+1)


peakidmtx=np.zeros((pairlist.shape[0],max_len))
for j in range(pairlist.shape[0]):
    peakstart=int(pairlist[j,0].item())
    peakend=int(pairlist[j,1].item())
    prompos=int(pairlist_prom[j,0].item())
    enha=list(range(peakstart,(peakend+1)))
    enha.remove(prompos)
    peaknum_gene=peakend-peakstart+1
    peakidmtx[j,0]=prompos
    peakidmtx[j,1:peaknum_gene]=enha

posmtx=np.zeros((pairlist.shape[0],max_len))
posmtx=posmtx-1
for j in range(pairlist.shape[0]):
    peakstart=int(pairlist[j,0].item())
    peakend=int(pairlist[j,1].item())
    prompos=int(pairlist_prom[j,0].item())
    enha=list(range(peakstart,(peakend+1)))
    enha.remove(prompos)
    genepos=int(pairlist[j,2].item())
    peaknum_gene=peakend-peakstart+1
    posvec=((((peaklist[enha,0]+peaklist[enha,1])/2)-genepos))
    posmtx[j,1:peaknum_gene]=posvec
    posmtx[j,0]=0

peaklist_full=pd.read_csv(samplename+"/peaks_extend.bed",sep="\t",header=None)
peaklist_full=peaklist_full.to_numpy()


enha_bed=pd.read_csv(samplename+"/F5.hg38.enhancers.bed",sep="\t",header=None)
enha_bed=enha_bed.to_numpy()

enha_matrix=np.zeros((pairlist.shape[0],max_len))
enha_matrix=enha_matrix-1
for gn in range(pairlist.shape[0]):
    print(gn)
    peaknum=(posmtx[gn,:]!=(-1)).sum()
    prom_id=peakidmtx[gn,0].astype(int)
    tmphic=enha_bed[enha_bed[:,0]==peaklist_full[prom_id,0],:].copy()
    for p1 in range(peaknum):
        p1_id=peakidmtx[gn,p1].astype(int)
        p1_in=(tmphic[:,2]>=peaklist_full[p1_id,1])&(tmphic[:,1]<=peaklist_full[p1_id,2])
        maemuki=tmphic[(p1_in),:]
        enha_matrix[gn,p1]=int(maemuki.shape[0]>0)

np.save(samplename+"/enha.npy",enha_matrix)
