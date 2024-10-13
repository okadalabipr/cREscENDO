import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams["pdf.fonttype"] = 42
import sys


args = sys.argv

samplename=str(args[1])


fname=samplename+"/Deeplift_full_ver2.npy"
enhamtx=np.load(fname)

genelist=pd.read_csv(samplename+"/pair_300000.csv",header=None)
genelist=genelist[[0]].to_numpy()
genelist=genelist[:,0]


grad=enhamtx.max(axis=1)
useidx=grad!=0

pairlist_prom=pd.read_csv(samplename+"/pair_promoter.csv",header=None)
pairlist_prom=pairlist_prom[[1,2,3]].to_numpy()

pairlist=pd.read_csv(samplename+"/pair_300000.csv",header=None)
pairlist=pairlist[[1,2,3]].to_numpy()

promenhatag=np.load(samplename+"/promenhatag.npy")

cellclus=pd.read_csv(samplename+"/cellcluster.txt",sep=",",header=None)
cellclus=cellclus.to_numpy()
cellclus=cellclus[:,0]


#################################



RNAmatrix=pd.read_csv(samplename+"/log1praw.csv",sep=",")
RNAmatrix=RNAmatrix.to_numpy()
rnanormfac=RNAmatrix.mean(axis=0)
RNAmatrix_n=RNAmatrix/rnanormfac
RNAmatrix_b=((RNAmatrix_n.transpose(1,0)-RNAmatrix_n.mean(axis=1))/RNAmatrix_n.std(axis=1)).transpose(1,0)



enhamtx_n=enhamtx.copy()
enhatmp=enhamtx.transpose(0,2,1)[(promenhatag==4)]
enhatmp=(enhatmp-enhatmp.mean(axis=0))/enhatmp.std(axis=0)
enhamtx_n.transpose(0,2,1)[(promenhatag==4)]=enhatmp

enhamtxuse=enhamtx.transpose(0,2,1)
enhamtxuse=enhamtxuse[(promenhatag==4)]
enhavar=enhamtxuse[:,cellclus==0].var(axis=1)

enhatumornorm=((enhamtxuse-enhamtxuse.mean(axis=0))/enhamtxuse.std(axis=0))
enhamtxnorm=enhamtx.transpose(0,2,1)
normfac=np.abs(enhamtxnorm).sum(axis=(0,1))
enhamtxnorm=enhamtxnorm/normfac
enhamtxuse=enhamtxnorm[(promenhatag==2)|(promenhatag==4)]


peakidmtx=np.zeros((enhamtx.shape[0],enhamtx.shape[2]))
for j in range(enhamtx.shape[0]):
    peakstart=int(pairlist[j,0].item())
    peakend=int(pairlist[j,1].item())
    prompos=int(pairlist_prom[j,0].item())
    enha=list(range(peakstart,(peakend+1)))
    enha.remove(prompos)
    peaknum_gene=peakend-peakstart+1
    peakidmtx[j,0]=prompos
    peakidmtx[j,1:peaknum_gene]=enha

ATACr=np.load(samplename+"/atac_count.npy")
cno=ATACr.sum(axis=0)
cnomin,cnomax=np.percentile(cno, [10,90])

normfac=np.abs(enhamtx).sum(axis=(0,2))
cnotag=(cno>cnomin)&(cno<cnomax)

tfidxr=np.load(samplename+"/enhamotif.npy")
tfidx=np.load(samplename+"/enhamotif_mix.npy")

tfact=pd.read_csv(samplename+"/motif_cell_activity.csv",index_col=0)
tfact=tfact.to_numpy()
tfactnorm=tfact
zeb1min,zeb1max=np.percentile(tfactnorm[(cellclus==0)&(cnotag==1),1], [10,90])
cremmin,cremmax=np.percentile(tfactnorm[(cellclus==0)&(cnotag==1),2], [10,90])

###############
enhamtxuse=enhamtx_n.transpose(0,2,1)/10
enhamtxuse=enhamtxuse[(promenhatag==4)]
enhamtxusepow=pow(2, enhamtxuse.transpose(1,0))

enha = sc.AnnData(enhamtxusepow[(cellclus==0)&(cnotag==1)])
sc.tl.pca(enha, svd_solver="arpack")
sc.pp.neighbors(enha, n_neighbors=10, n_pcs=40)
sc.tl.umap(enha)
sc.tl.tsne(enha)

sc.tl.leiden(
    enha,
    resolution=0.2
)

cellleiden=np.zeros((cellclus.shape[0]))
cellleiden=cellleiden-1
cellleiden[(cellclus==0)&(cnotag==1)]=enha.obs["leiden"].to_numpy()
np.save(samplename+"/cellleiden.npy",cellleiden)
pd.DataFrame(cellleiden).to_csv(samplename+"/cellleiden.txt",sep="\t",header=False, index=False)

enha.obs["leiden"] = pd.Categorical(list(cellleiden[cellclus==0]))
enha.obs["ZEB1"] = tfactnorm[(cellclus==0)&(cnotag==1),1]
enha.obs["CREM"] = tfactnorm[(cellclus==0)&(cnotag==1),2]

sc.pl.umap(enha, frameon=False, color="ZEB1",save="glioma_enhancer3zeb1.pdf",vmin=zeb1min, vmax=zeb1max)
sc.pl.umap(enha, frameon=False, color="CREM",save="glioma_enhancer3crem.pdf",vmin=cremmin, vmax=cremmax)
sc.pl.umap(enha, frameon=False,legend_loc="on data", color="leiden",save="glioma_enhancer2.pdf")



DAEmtx=pd.DataFrame(enha.uns["rank_genes_groups"]["names"])
pval=pd.DataFrame(enha.uns["rank_genes_groups"]["pvals_adj"])
fc=pd.DataFrame(enha.uns["rank_genes_groups"]["logfoldchanges"])


genelist=pd.read_csv(samplename+"/pair_300000.csv",header=None)
genelist=genelist[[0]].to_numpy()
genelist=genelist[:,0]

enhaclusgenemtxall=np.zeros((DAEmtx.shape[1],promenhatag.shape[0],promenhatag.shape[1]))

for i in range(DAEmtx.shape[1]):
   tmp=np.zeros(enhamtxuse.shape[0])
   tmp=tmp
   enhaname=DAEmtx.iloc[:,i]
   usetag=enhaname[((pval.iloc[:,i]<0.0001)&(fc.iloc[:,i]>0.01))].to_numpy()
   usetag=usetag.astype(int)
   tmp[usetag]=1
   enhaclusgenemtx=np.zeros(promenhatag.shape)
   enhaclusgenemtx=enhaclusgenemtx
   enhaclusgenemtx[(promenhatag==4)]=tmp
   enhaclusgenemtxall[i]=enhaclusgenemtx

np.save(samplename+"/enhaclusgenemtxall.txt",enhaclusgenemtxall)

peaklistout=pd.read_csv(samplename+"/peaks_extend.bed",sep="\t",header=None)
sox2cre=peaklistout.iloc[np.unique(peakidmtx[(tfidxr[:,:,0]==1)&(enhaclusgenemtxall[0]==1)]).astype(int),:]
sox2noncre=peaklistout.iloc[np.unique(peakidmtx[(tfidxr[:,:,0]==1)&((promenhatag!=2)&(promenhatag!=4))]).astype(int),:]

sox2cre.to_csv(samplename+"/sox2cre.bed",sep="\t",header=False, index=False)
sox2noncre.to_csv(samplename+"/sox2noncre.bed",sep="\t",header=False, index=False)
