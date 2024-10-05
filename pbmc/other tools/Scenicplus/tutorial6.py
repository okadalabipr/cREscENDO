import dill
import scanpy as sc
import os
import warnings
warnings.filterwarnings("ignore")
import pandas
import pyranges
# Set stderr to null to avoid strange messages from ray
import sys
_stderr = sys.stderr
null = open(os.devnull,'wb')

import sys
args = sys.argv
work_dir=str(args[1])

work_dir="scenicplus/scenicplus8"

scplus_obj = dill.load(open(os.path.join(work_dir, 'scplus_obj.pkl'),'rb'))


scplus_obj["direct_gene_based_AUC"]

scout=scplus_obj.uns['eRegulon_metadata'][["Region","Gene","R2G_importance"]]
scout.to_csv(os.path.join(work_dir, "scenicplus_R2G.txt"),sep="\t")

scout=scplus_obj.uns['eRegulon_metadata'][["Region","Gene","R2G_importance_x_rho"]]
scout.to_csv(os.path.join(work_dir, "scenicplus_R2G_rho.txt"),sep="\t")

scout=scplus_obj.uns['eRegulon_metadata'][["Region","Gene","R2G_importance_x_abs_rho"]]
scout.to_csv(os.path.join(work_dir, "scenicplus_R2G_absrho.txt"),sep="\t")

scout=scplus_obj.uns['eRegulon_metadata'][["Region_signature_name","Region"]]
scout2=scplus_obj.uns['eRegulon_AUC']
scout3=scout2["Region_based"]
scout3.index=scout3.index.str.replace('-10x_pbmc', '')
scout3["Cellid"]=scout3.index


np.unique(scout["Region_signature_name"])

import numpy as np
import pandas as pd

cellclus=pd.read_csv("pbmc_after/cellcluster_all.txt",sep=",",header=None)
#cellclus=cellclus.to_numpy()
#cellclus=cellclus[:,0]

cellid=pd.read_csv("pbmc_after/celllist_scVI.csv",sep="\t",header=None)


cellmtx=pd.concat([cellid,cellclus],axis=1)
cellmtx.columns=["Cellid","Cluster"]

scout4=pd.merge(scout3, cellmtx, on='Cellid', how='left')

scout4.to_csv(os.path.join(work_dir, "scout4.csv"))


scout41=scout4[scout4["Cluster"]==1]
scout42=scout4[(scout4["Cluster"]==1)|(scout4["Cluster"]==2)]
scout43=scout4[scout4["Cluster"]==3]
scout44=scout4[scout4["Cluster"]==4]

eregn=scout4.shape[1]-2


eregmtx=pd.concat([scout41.iloc[:,:eregn].mean(axis=0),scout42.iloc[:,:eregn].mean(axis=0),scout43.iloc[:,:eregn].mean(axis=0),scout44.iloc[:,:eregn].mean(axis=0)],axis=1)

eregmtx.to_csv(os.path.join(work_dir, "eregmtx.csv"))

scout5=scplus_obj.uns["RSS"]
scout6=scout5["GEX_celltype_region_based"]
scout7=scout6.T
scout8=scout7[["CD14+_Monocytes","CD4_T_cells","CD8_T_cells","B_cells"]]

scout8.to_csv(os.path.join(work_dir, "scout8.csv"))

##########################

import dill
import scanpy as sc
import os
import warnings
warnings.filterwarnings("ignore")
import pandas
import pyranges
# Set stderr to null to avoid strange messages from ray
import sys
import numpy as np
import pandas as pd
_stderr = sys.stderr
null = open(os.devnull,'wb')



work_dir="scenicplus/scenicplus0"
scplus_obj = dill.load(open(os.path.join(work_dir, 'scplus_obj.pkl'),'rb'))
cellclus=pd.read_csv("pbmc_after/cellcluster_all.txt",sep=",",header=None)
cellid=pd.read_csv("pbmc_after/celllist_scVI.csv",sep="\t",header=None)
cellmtx=pd.concat([cellid,cellclus],axis=1)
cellmtx.columns=["Cellid","Cluster"]
scout2=scplus_obj.uns['eRegulon_AUC']
scout3=scout2["Region_based"]
scout3.index=scout3.index.str.replace('-10x_pbmc', '')
scout3["Cellid"]=scout3.index
scout4=pd.merge(scout3, cellmtx, on='Cellid', how='left')
scout4.to_csv(os.path.join(work_dir, "scout4.csv"))

work_dir="scenicplus/scenicplus4"
scplus_obj = dill.load(open(os.path.join(work_dir, 'scplus_obj.pkl'),'rb'))
cellclus=pd.read_csv("pbmc_after/cellcluster_all.txt",sep=",",header=None)
cellid=pd.read_csv("pbmc_after/celllist_scVI.csv",sep="\t",header=None)
cellmtx=pd.concat([cellid,cellclus],axis=1)
cellmtx.columns=["Cellid","Cluster"]
scout2=scplus_obj.uns['eRegulon_AUC']
scout3=scout2["Region_based"]
scout3.index=scout3.index.str.replace('-10x_pbmc', '')
scout3["Cellid"]=scout3.index
scout4=pd.merge(scout3, cellmtx, on='Cellid', how='left')
scout4.to_csv(os.path.join(work_dir, "scout4.csv"))

work_dir="scenicplus/scenicplus7"
scplus_obj = dill.load(open(os.path.join(work_dir, 'scplus_obj.pkl'),'rb'))
cellclus=pd.read_csv("pbmc_after/cellcluster_all.txt",sep=",",header=None)
cellid=pd.read_csv("pbmc_after/celllist_scVI.csv",sep="\t",header=None)
cellmtx=pd.concat([cellid,cellclus],axis=1)
cellmtx.columns=["Cellid","Cluster"]
scout2=scplus_obj.uns['eRegulon_AUC']
scout3=scout2["Region_based"]
scout3.index=scout3.index.str.replace('-10x_pbmc', '')
scout3["Cellid"]=scout3.index
scout4=pd.merge(scout3, cellmtx, on='Cellid', how='left')
scout4.to_csv(os.path.join(work_dir, "scout4.csv"))

work_dir="scenicplus/scenicplus2"
scplus_obj = dill.load(open(os.path.join(work_dir, 'scplus_obj.pkl'),'rb'))
cellclus=pd.read_csv("pbmc_after/cellcluster_all.txt",sep=",",header=None)
cellid=pd.read_csv("pbmc_after/celllist_scVI.csv",sep="\t",header=None)
cellmtx=pd.concat([cellid,cellclus],axis=1)
cellmtx.columns=["Cellid","Cluster"]
scout2=scplus_obj.uns['eRegulon_AUC']
scout3=scout2["Region_based"]
scout3.index=scout3.index.str.replace('-10x_pbmc', '')
scout3["Cellid"]=scout3.index
scout4=pd.merge(scout3, cellmtx, on='Cellid', how='left')
scout4.to_csv(os.path.join(work_dir, "scout4.csv"))

work_dir="scenicplus/scenicplus8"
scplus_obj = dill.load(open(os.path.join(work_dir, 'scplus_obj.pkl'),'rb'))
cellclus=pd.read_csv("pbmc_after/cellcluster_all.txt",sep=",",header=None)
cellid=pd.read_csv("pbmc_after/celllist_scVI.csv",sep="\t",header=None)
cellmtx=pd.concat([cellid,cellclus],axis=1)
cellmtx.columns=["Cellid","Cluster"]
scout2=scplus_obj.uns['eRegulon_AUC']
scout3=scout2["Region_based"]
scout3.index=scout3.index.str.replace('-10x_pbmc', '')
scout3["Cellid"]=scout3.index
scout4=pd.merge(scout3, cellmtx, on='Cellid', how='left')
scout4.to_csv(os.path.join(work_dir, "scout4.csv"))

#########################



peaklist=pd.read_csv('pbmc_after/peaks.bed',sep="\t",header=None)
peaklist=peaklist[[1,2]].to_numpy()
pairlist=pd.read_csv('pbmc_after/pair_300000.csv',header=None)
pairlist=pairlist[[1,2,3]].to_numpy()

pairlist_prom=pd.read_csv('pbmc_after/pair_promoter.csv',header=None)
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


genelist=pd.read_csv('pbmc_after/pair_300000.csv',header=None)
genelist=genelist[[0]].to_numpy()


R2G=pd.read_csv(os.path.join(work_dir, 'scenicplus_R2G_absrho.txt'),sep="\t")
R2G=R2G.iloc[:,1:4]

scout.index=R2G.index
((R2G["Region"]==scout["Region"])==0).sum()
R2G=pd.concat([R2G, scout["Region_signature_name"]], axis=1)

R2G["Chr"]=R2G["Region"].apply(lambda s: s.split(":")[0])
R2G["tag"]=R2G["Region"].apply(lambda s: s.split(":")[1])
R2G["Start"]=R2G["tag"].apply(lambda s: s.split("-")[0])
R2G["End"]=R2G["tag"].apply(lambda s: s.split("-")[1])
R2G=R2G[["Chr","Start","End","Gene","R2G_importance_x_abs_rho","Region_signature_name"]]
R2G["Start"]=R2G["Start"].astype(float)
R2G["End"]=R2G["End"].astype(float)
R2G=R2G.to_numpy()

peaklist_full=pd.read_csv('pbmc_after/peaks_extend.bed',sep="\t",header=None)
peaklist_full=peaklist_full.to_numpy()

R2G_absrho_matrix=np.zeros((pairlist.shape[0],max_len),dtype="object")
for gn in range(pairlist.shape[0]):
    print(gn)
    peaknum=(posmtx[gn,:]!=(-1)).sum()
    for p1 in range(peaknum):
        p1_id=peakidmtx[gn,p1].astype(int)
        tmphic=R2G[(R2G[:,0]==peaklist_full[p1_id,0])&(R2G[:,3]==genelist[gn]),:].copy()
        p1_in=(tmphic[:,2]>=peaklist_full[p1_id,1])&(tmphic[:,1]<=peaklist_full[p1_id,2])
        maemuki=tmphic[(p1_in),:]
        if len(maemuki)>0:
            R2G_absrho_matrix[gn,p1]=maemuki[np.argmax(maemuki[:,4]),5]

np.save(os.path.join(work_dir, "eregulonmtx.npy"),R2G_absrho_matrix)

