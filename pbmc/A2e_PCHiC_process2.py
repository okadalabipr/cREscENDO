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

bait_hg19=pd.read_csv(samplename+"/bait_hg19_tag.bed",sep="\t",header=None,names=("chr","start","end","tag"))
bait_hg38=pd.read_csv(samplename+"/bait_hg38_tag.bed",sep="\t",header=None,names=("chr","start","end","tag"))
chrs = ["chr"+str(i) for i in range(1,23)] + ["chrX", "chrY"]
bait_hg38_s=bait_hg38[bait_hg38["chr"].isin(chrs)]

oe_hg19=pd.read_csv(samplename+"/oe_hg19_tag.bed",sep="\t",header=None,names=("chr","start","end","tag"))
oe_hg38=pd.read_csv(samplename+"/oe_hg38_tag.bed",sep="\t",header=None,names=("chr","start","end","tag"))
chrs = ["chr"+str(i) for i in range(1,23)] + ["chrX", "chrY"]
oe_hg38_s=oe_hg38[oe_hg38["chr"].isin(chrs)]

bait_hg38_sd=bait_hg38_s.drop_duplicates(subset="tag")
bait_hg19_m=pd.merge(bait_hg19, bait_hg38_sd, on="tag", how="left")
oe_hg38_sd=oe_hg38_s.drop_duplicates(subset="tag")
oe_hg19_m=pd.merge(oe_hg19, oe_hg38_sd, on="tag", how="left")

bait_hg19_to_38=bait_hg19_m.iloc[:,4:7]
bait_hg19_to_38.columns=["baitChr","baitStart","baitEnd"]
oe_hg19_to_38=oe_hg19_m.iloc[:,4:7]
oe_hg19_to_38.columns=["oeChr","oeStart","oeEnd"]

hg19_m=pd.concat([bait_hg19_to_38, oe_hg19_to_38],axis=1)
hg19_m["signal"]=hic.iloc[:,11:28].to_numpy().max(axis=1)
hg19_mdn=hg19_m.dropna(how="any")

hic_sn=hg19_mdn.to_numpy()

pchic_matrix=np.zeros((pairlist.shape[0],max_len))
for gn in range(pairlist.shape[0]):
    print(gn)
    peaknum=(posmtx[gn,:]!=(-1)).sum()
    prom_id=peakidmtx[gn,0].astype(int)
    tmphic=hic_sn[hic_sn[:,0]==peaklist_full[prom_id,0],:].copy()
    prom_in=(tmphic[:,2]>=peaklist_full[prom_id,1])&(tmphic[:,1]<=peaklist_full[prom_id,2])
    tmphic_in=tmphic[(prom_in),:]
    for p1 in range(peaknum):
        p1_id=peakidmtx[gn,p1].astype(int)
        p1_in=(tmphic_in[:,5]>=peaklist_full[p1_id,1])&(tmphic_in[:,4]<=peaklist_full[p1_id,2])
        maemuki=tmphic_in[(p1_in),:]
        pchic_matrix[gn,p1]=maemuki.shape[0]

np.save(samplename+"/pchic_matrix.npy",pchic_matrix)



