import numpy as np
import pandas as pd

peaklist=pd.read_csv('pbmc_fold1/peaks.bed',sep="\t",header=None)
peaklist=peaklist[[1,2]].to_numpy()
pairlist=pd.read_csv('pbmc_fold1/pair_300000.csv',header=None)
pairlist=pairlist[[1,2,3]].to_numpy()

pairlist_prom=pd.read_csv('pbmc_fold1/pair_promoter.csv',header=None)
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


peaklist_full=pd.read_csv('pbmc_fold1/peaks_extend.bed',sep="\t",header=None)
peaklist_full=peaklist_full.to_numpy()

cicero=pd.read_csv("cicero_for_paper/try1/PBMC_cicero_connections.txt",sep="\t")

cicero_p1=cicero["Peak1"].str.split(':',expand=True)
cicero_p2=cicero["Peak2"].str.split(':',expand=True)

cicero_p1=cicero_p1.to_numpy()
cicero_p2=cicero_p2.to_numpy()
cicero_v=cicero["coaccess"].to_numpy()
cicero_all=np.concatenate([cicero_p1,cicero_p2,cicero_v[:,np.newaxis]],1)
cicero_all[:,1:3]=cicero_all[:,1:3].astype(int)
cicero_all[:,4:6]=cicero_all[:,4:6].astype(int)

cicero_matrix=np.zeros((pairlist.shape[0],max_len))
for gn in range(pairlist.shape[0]):
    print(gn)
    peaknum=(posmtx[gn,:]!=(-1)).sum()
    prom_id=peakidmtx[gn,0].astype(int)
    tmphic=cicero_all[cicero_all[:,0]==peaklist_full[prom_id,0],:].copy()
    prom_in=(tmphic[:,2]>=peaklist_full[prom_id,1])&(tmphic[:,1]<=peaklist_full[prom_id,2])
    tmphic_in_m=tmphic[(prom_in),:]
    prom_in=(tmphic[:,5]>=peaklist_full[prom_id,1])&(tmphic[:,4]<=peaklist_full[prom_id,2])
    tmphic_in_u=tmphic[(prom_in),:]
    for p1 in range(peaknum):
        p1_id=peakidmtx[gn,p1].astype(int)
        p1_in=(tmphic_in_m[:,5]>=peaklist_full[p1_id,1])&(tmphic_in_m[:,4]<=peaklist_full[p1_id,2])
        maemuki=tmphic_in_m[(p1_in),:]
        p1_in=(tmphic_in_u[:,2]>=peaklist_full[p1_id,1])&(tmphic_in_u[:,1]<=peaklist_full[p1_id,2])
        ushiromuki=tmphic_in_u[(p1_in),:]
        alltmp=np.concatenate([maemuki,ushiromuki],0)
        if alltmp.shape[0]>0:
            cicero_matrix[gn,p1]=alltmp[:,6].max()

np.save("cicero_for_paper/try1/cicero_matrix.npy",cicero_matrix)