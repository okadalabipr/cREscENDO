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

pd.DataFrame(genelist).to_csv(samplename+"/genelistfull.txt",sep="\t",header=False, index=False)
strogene=genelist[(promenhatag==4).sum(axis=1)>0]
pd.DataFrame(strogene).to_csv(samplename+"/strogenelist.txt",sep="\t",header=False, index=False)

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


###############################################

enhamtxuse=enhamtx_n.transpose(0,2,1)/10
enhamtxuse=enhamtxuse[(promenhatag==4)]
enhamtxusepow=pow(2, enhamtxuse.transpose(1,0))
enhaall = sc.AnnData(enhamtxusepow)
sc.tl.pca(enhaall, svd_solver="arpack")
sc.pp.neighbors(enhaall, n_neighbors=10, n_pcs=40)
sc.tl.umap(enhaall)

enhaall.obs["leiden"] = pd.Categorical(list(cellclus))

enhaall.obs["celltype"] = ""
enhaall.obs["celltype"][cellclus==0]="Tumor"
enhaall.obs["celltype"][cellclus==1]="OPC"
enhaall.obs["celltype"][cellclus==2]="Astrocyte"
enhaall.obs["celltype"][cellclus==3]="Myeloid"

sc.pl.umap(enhaall, frameon=False,legend_loc="on data", color="celltype",save="glioma_enhancer_all3.pdf")


sox2minall,sox2maxall=np.percentile(tfactnorm[:,0], [10,90])
enhaall.obs["SOX2"] = tfactnorm[:,0]
sc.pl.umap(enhaall, frameon=False, color="SOX2",save="glioma_enhancer3sox2_all.pdf",vmin=sox2minall, vmax=sox2maxall)

sox2rna=pd.read_csv(samplename+"/sox2exp.txt",sep="\t",header=None)
sox2rna=sox2rna.to_numpy()
sox2rna=sox2rna[:,0]
sox2rnaminall,sox2rnamaxall=np.percentile(sox2rna, [10,90])
enhaall.obs["SOX2_RNA"] = sox2rna
sc.pl.umap(enhaall, frameon=False, color="SOX2_RNA",save="glioma_enhancer3sox2rna_all.pdf",vmin=sox2rnaminall, vmax=sox2rnamaxall)



#################

gexall=sc.read_h5ad(samplename+"/gexm_allcell.h5ad")
sc.pp.normalize_per_cell(gexall, counts_per_cell_after=1e4)
sc.pp.log1p(gexall)
sc.tl.pca(gexall, svd_solver="arpack")
sc.pp.neighbors(gexall, n_neighbors=10, n_pcs=40)
sc.tl.umap(gexall)
sc.tl.tsne(gexall)

gexall.obs["celltype"] = ""
gexall.obs["celltype"][cellclus==0]="Tumor"
gexall.obs["celltype"][cellclus==1]="OPC"
gexall.obs["celltype"][cellclus==2]="Astrocyte"
gexall.obs["celltype"][cellclus==3]="Myeloid"

sc.pl.umap(gexall, frameon=False,legend_loc="on data", color="celltype",save="glioma_gex_all3.pdf")
