import numpy as np
import pandas as pd

peaklist=pd.read_csv('pbmc_fold1/peaks.bed',sep="\t",header=None)
peaklist=peaklist[[1,2]].to_numpy()
pairlist=pd.read_csv('pbmc_fold1/pair_300000.csv',header=None)
pairlist=pairlist[[1,2,3]].to_numpy()

pairlist_prom=pd.read_csv('pbmc_fold1/pair_promoter.csv',header=None)
pairlist_prom=pairlist_prom[[1,2,3]].to_numpy()

genelist=pd.read_csv('pbmc_fold1/pair_promoter.csv',header=None)
genelist=genelist[0].to_numpy()

max_len=int((pairlist[:,1]-pairlist[:,0]).max()+1)

peaklist_full=pd.read_csv('pbmc_fold1/peaks_extend.bed',sep="\t",header=None)
peaklist_full=peaklist_full.to_numpy()


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

###########
    
for i in range(1,11):
    fname="archr_for_paper/try"+str(i)+"/PBMC_ArchR_links_genes.txt"
    archr=pd.read_csv(fname,sep="\t",header=None)
    archr_tmp=archr[1].str.split('_',expand=True)
    archr_all = pd.concat([archr[0], archr_tmp,archr[2]], axis=1)
    archr_all=archr_all.to_numpy()
    archr_all[:,2:4]=archr_all[:,2:4].astype(int)
    archr_matrix=np.zeros((pairlist.shape[0],max_len))
    for gn in range(pairlist.shape[0]):
        print(gn)
        peaknum=(posmtx[gn,:]!=(-1)).sum()
        genetmp=genelist[gn]
        tmphic=archr_all[archr_all[:,0]==genetmp,:].copy()
        for p1 in range(peaknum):
            p1_id=peakidmtx[gn,p1].astype(int)
            p1_in=(tmphic[:,3]>=peaklist_full[p1_id,1])&(tmphic[:,2]<=peaklist_full[p1_id,2])
            maemuki=tmphic[(p1_in),:]
            if maemuki.shape[0]>0:
                archr_matrix[gn,p1]=maemuki[:,4].max()
    fname="archr_for_paper/try"+str(i)+"/archr_matrix.npy"
    np.save(fname,archr_matrix)