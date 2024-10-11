import scanpy as sc
import pandas as pd
import numpy as np
import sys

args = sys.argv
samplename=str(args[1])
cellcluspath=str(args[2])


cellclus=pd.read_csv(cellcluspath,sep=",",header=None)
cellclus=cellclus.to_numpy()
cellclus=cellclus[:,0]


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



enhamtx_n=enhamtx.copy()
enhatmp=enhamtx.transpose(0,2,1)[(promenhatag==4)]
enhatmp=(enhatmp-enhatmp.mean(axis=0))/enhatmp.std(axis=0)
enhamtx_n.transpose(0,2,1)[(promenhatag==4)]=enhatmp

enhamtxuse=enhamtx_n.transpose(0,2,1)/10
enhamtxuse=enhamtxuse[(promenhatag==4)]
enhamtxusepow=pow(2, enhamtxuse.transpose(1,0))
enhaall = sc.AnnData(enhamtxusepow)
sc.tl.pca(enhaall, svd_solver="arpack")
sc.pp.neighbors(enhaall, n_neighbors=10, n_pcs=40)
sc.tl.umap(enhaall)
enhaall.obs["leiden"] = pd.Categorical(list(cellclus))
sc.tl.rank_genes_groups(enhaall, "leiden", method="wilcoxon")

DAEmtxall=pd.DataFrame(enhaall.uns["rank_genes_groups"]["names"])
pvalall=pd.DataFrame(enhaall.uns["rank_genes_groups"]["pvals_adj"])
fcall=pd.DataFrame(enhaall.uns["rank_genes_groups"]["logfoldchanges"])


enhaclusgenemtxall_allenha=np.zeros((DAEmtxall.shape[1],promenhatag.shape[0],promenhatag.shape[1]))

for i in range(DAEmtxall.shape[1]):
	enhaname=DAEmtxall.iloc[:,i]
	usetag=enhaname[((pvalall.iloc[:,i]<0.05)&(fcall.iloc[:,i]>0))].to_numpy()
	for j in range(usetag.shape[0]):
		tmp2=usetag[j].split("_")
		xtag=np.where(fullgenename == tmp2[0])[0][0]
		ytag=int(tmp2[2])
		enhaclusgenemtxall_allenha[i,xtag,ytag]=1


np.save(samplename+"/enha_clus_mtx_deeplift.npy",enhaclusgenemtxall_allenha)

################



