import numpy as np
import pandas as pd
import os

import sys
args = sys.argv
work_dir=str(args[1])

work_dir="../scenicplus8"


peaklist=pd.read_csv('peaks.bed',sep="\t",header=None)
peaklist=peaklist[[1,2]].to_numpy()
pairlist=pd.read_csv('pair_300000.csv',header=None)
pairlist=pairlist[[1,2,3]].to_numpy()

pairlist_prom=pd.read_csv('pair_promoter.csv',header=None)
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


genelist=pd.read_csv('pair_300000.csv',header=None)
genelist=genelist[[0]].to_numpy()

R2G=pd.read_csv(os.path.join(work_dir, 'scenicplus_R2G.txt'),sep="\t")
R2G=R2G.iloc[:,1:4]
R2G["Chr"]=R2G["Region"].apply(lambda s: s.split(":")[0])
R2G["tag"]=R2G["Region"].apply(lambda s: s.split(":")[1])
R2G["Start"]=R2G["tag"].apply(lambda s: s.split("-")[0])
R2G["End"]=R2G["tag"].apply(lambda s: s.split("-")[1])
R2G=R2G[["Chr","Start","End","Gene","R2G_importance"]]
R2G["Start"]=R2G["Start"].astype(float)
R2G["End"]=R2G["End"].astype(float)
R2G=R2G.to_numpy()

peaklist_full=pd.read_csv('peaks_extend.bed',sep="\t",header=None)
peaklist_full=peaklist_full.to_numpy()

R2G_matrix=np.zeros((pairlist.shape[0],max_len))
for gn in range(pairlist.shape[0]):
    print(gn)
    peaknum=(posmtx[gn,:]!=(-1)).sum()
    for p1 in range(peaknum):
        p1_id=peakidmtx[gn,p1].astype(int)
        tmphic=R2G[(R2G[:,0]==peaklist_full[p1_id,0])&(R2G[:,3]==genelist[gn]),:].copy()
        p1_in=(tmphic[:,2]>=peaklist_full[p1_id,1])&(tmphic[:,1]<=peaklist_full[p1_id,2])
        maemuki=tmphic[(p1_in),:]
        if len(maemuki)>0:
            R2G_matrix[gn,p1]=maemuki[:,4].max()

np.save(os.path.join(work_dir, "R2G_matrix.npy"),R2G_matrix)


#############


R2G=pd.read_csv(os.path.join(work_dir, 'scenicplus_R2G_rho.txt'),sep="\t")
R2G=R2G.iloc[:,1:4]
R2G["Chr"]=R2G["Region"].apply(lambda s: s.split(":")[0])
R2G["tag"]=R2G["Region"].apply(lambda s: s.split(":")[1])
R2G["Start"]=R2G["tag"].apply(lambda s: s.split("-")[0])
R2G["End"]=R2G["tag"].apply(lambda s: s.split("-")[1])
R2G=R2G[["Chr","Start","End","Gene","R2G_importance_x_rho"]]
R2G["Start"]=R2G["Start"].astype(float)
R2G["End"]=R2G["End"].astype(float)
R2G=R2G.to_numpy()

peaklist_full=pd.read_csv('peaks_extend.bed',sep="\t",header=None)
peaklist_full=peaklist_full.to_numpy()

R2G_rho_matrix=np.zeros((pairlist.shape[0],max_len))
for gn in range(pairlist.shape[0]):
    print(gn)
    peaknum=(posmtx[gn,:]!=(-1)).sum()
    for p1 in range(peaknum):
        p1_id=peakidmtx[gn,p1].astype(int)
        tmphic=R2G[(R2G[:,0]==peaklist_full[p1_id,0])&(R2G[:,3]==genelist[gn]),:].copy()
        p1_in=(tmphic[:,2]>=peaklist_full[p1_id,1])&(tmphic[:,1]<=peaklist_full[p1_id,2])
        maemuki=tmphic[(p1_in),:]
        if len(maemuki)>0:
            R2G_rho_matrix[gn,p1]=maemuki[:,4].max()

np.save(os.path.join(work_dir, "R2G_rho_matrix.npy"),R2G_rho_matrix)


#################


R2G=pd.read_csv(os.path.join(work_dir, 'scenicplus_R2G_absrho.txt'),sep="\t")
R2G=R2G.iloc[:,1:4]
R2G["Chr"]=R2G["Region"].apply(lambda s: s.split(":")[0])
R2G["tag"]=R2G["Region"].apply(lambda s: s.split(":")[1])
R2G["Start"]=R2G["tag"].apply(lambda s: s.split("-")[0])
R2G["End"]=R2G["tag"].apply(lambda s: s.split("-")[1])
R2G=R2G[["Chr","Start","End","Gene","R2G_importance_x_abs_rho"]]
R2G["Start"]=R2G["Start"].astype(float)
R2G["End"]=R2G["End"].astype(float)
R2G=R2G.to_numpy()

peaklist_full=pd.read_csv('peaks_extend.bed',sep="\t",header=None)
peaklist_full=peaklist_full.to_numpy()

R2G_absrho_matrix=np.zeros((pairlist.shape[0],max_len))
for gn in range(pairlist.shape[0]):
    print(gn)
    peaknum=(posmtx[gn,:]!=(-1)).sum()
    for p1 in range(peaknum):
        p1_id=peakidmtx[gn,p1].astype(int)
        tmphic=R2G[(R2G[:,0]==peaklist_full[p1_id,0])&(R2G[:,3]==genelist[gn]),:].copy()
        p1_in=(tmphic[:,2]>=peaklist_full[p1_id,1])&(tmphic[:,1]<=peaklist_full[p1_id,2])
        maemuki=tmphic[(p1_in),:]
        if len(maemuki)>0:
            R2G_absrho_matrix[gn,p1]=maemuki[:,4].max()

np.save(os.path.join(work_dir, "R2G_absrho_matrix.npy"),R2G_absrho_matrix)

