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

fullgenename=pd.read_csv(samplename+"/pair_300000.csv",header=None,names=("gene","start","end","pos"))
fullgenename=fullgenename["gene"]
fullgenename=fullgenename.to_numpy()

promenhatag=np.load(samplename+'/promenhatag.npy')


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

