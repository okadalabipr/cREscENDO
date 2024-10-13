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

hic=pd.read_csv(samplename+"/PCHiC_peak_matrix_cutoff5.tsv",sep="\t")
hic_s=hic.iloc[:,[0,1,2,5,6,7]]
hic_s["signal"]=hic.iloc[:,11:28].to_numpy().max(axis=1)
hic_s["baitChr"]="chr"+hic_s["baitChr"].astype(str)
hic_s["oeChr"]="chr"+hic_s["oeChr"].astype(str)

bait=hic_s.iloc[:,0:3]
oe=hic_s.iloc[:,3:6]

bait["tag"]=bait["baitChr"]+str("-")+bait["baitStart"].astype("str")+str("-")+bait["baitEnd"].astype("str")
bait.to_csv(samplename+"/bait_hg19_tag.bed", sep="\t",header=False, index=False)
oe["tag"]=oe["oeChr"]+str("-")+oe["oeStart"].astype("str")+str("-")+oe["oeEnd"].astype("str")
oe.to_csv(samplename+"/oe_hg19_tag.bed", sep="\t",header=False, index=False)

