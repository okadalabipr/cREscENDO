import scanpy as sc
import pandas as pd
import numpy as np
import sys

args = sys.argv
samplename=str(args[1])
cellcluspath=str(args[2])

peaknum=1000

pairlist=pd.read_csv(samplename+"/pair_300000.csv",header=None,names=("gene","start","end","pos"))
pairlist=pairlist[["start","end","pos"]].to_numpy()
max_len=int((pairlist[:,1]-pairlist[:,0]).max()+1)

RNAembed=pd.read_csv(samplename+"/pca_50dim_rnaembed_raw.txt",sep=" ",header=None)
RNAembed=RNAembed.to_numpy()
traincellnum=np.load(samplename+"/traincellnum.npy")

fullcell=RNAembed.shape[0]
traincell=traincellnum #10000 11830 11740
testcell=fullcell #1893
testtag=0

df=pd.read_csv(samplename+"/Deeplift_enhamtx_t2_raw.csv", header=0 , index_col=0)
df.columns.str.contains("promoter")
df2=df.iloc[:,~df.columns.str.contains("promoter")]

enha = sc.AnnData(df2, df2.index.to_frame(), df2.columns.to_frame())
sc.tl.pca(enha, svd_solver="arpack")
sc.pp.neighbors(enha, n_neighbors=10, n_pcs=40)
sc.tl.umap(enha)

sc.pl.umap(enha, frameon=False, save="pbmc_for_paper_enhancer.png")
cellclus=pd.read_csv(cellcluspath,header=None)
cellclus=cellclus.to_numpy()
cellclus=cellclus[:,0]-1
cellclus=cellclus[testtag:testtag+testcell]
enha.obs["rnacluster"] = pd.Categorical(list(cellclus))

sc.pl.umap(enha, frameon=False, save=samplename+"/myplot_deeplift_out.png",color="rnacluster")
sc.pl.umap(enha, frameon=False, save=samplename+"/myplot_deeplift.png",color="rnacluster",add_outline=True, legend_loc="on data")
sc.tl.rank_genes_groups(enha, "rnacluster", method="wilcoxon")

sc.pl.rank_genes_groups(enha, n_genes=25, sharey=False, save="DAE_deeplift_pbmc_for_paper.png")

enha.uns["rank_genes_groups"]["names"]

DAEmtx=pd.DataFrame(enha.uns["rank_genes_groups"]["names"])
pvalall=pd.DataFrame(enha.uns["rank_genes_groups"]["pvals_adj"])
fcall=pd.DataFrame(enha.uns["rank_genes_groups"]["logfoldchanges"])

DAEmtx.to_csv(samplename+"/DAE_mtx_deeplift.csv",header=True, index=True)

dae=pd.read_csv(samplename+"/DAE_mtx_deeplift.csv", header=0 , index_col=0)
dae=dae[0:peaknum]
dae_n=dae.to_numpy()
fullgenename=pd.read_csv(samplename+"/pair_300000.csv",header=None,names=("gene","start","end","pos"))
fullgenename=fullgenename["gene"]
fullgenename=fullgenename.to_numpy()


grad=np.load(samplename+"/allgrad_ssep_max.npy")

thv=np.percentile(grad[grad!=0], [40])
sthv=np.percentile(grad[grad!=0], [85])

promenhatag=np.zeros(grad.shape)
promenhatag[(grad>thv)&(ATAC_corr>0)]=2
promenhatag[(grad>sthv)&(ATAC_corr>0)]=4
promenhatag[:,0]=3


fname=samplename+"/Deeplift_full_ver2_all.npy"
enhamtx=np.load(fname)

enhathr=np.percentile(enhamtx.transpose(0,2,1),75,axis=2) #B,L
enhathrmtx=(enhamtx.transpose(1,0,2)>enhathr).transpose(1,0,2) #B,C,L

enhausemtx=((enhathrmtx.transpose(1,0,2)*((promenhatag==4)|(promenhatag==2))).sum(axis=2))/((promenhatag==4)|(promenhatag==2)).sum(axis=1)

enhausetagmtx=enhausemtx.transpose(1,0)>0.7

np.save(samplename+'/enhausetagmtx.npy',enhausetagmtx)
np.save(samplename+'/promenhatag.npy',promenhatag)


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
