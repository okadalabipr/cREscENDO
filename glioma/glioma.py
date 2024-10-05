import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
import sys
import math
import subprocess

args = sys.argv

samplename=str(args[1])


fname=samplename+"/Deeplift_full_ver2_mergenorm.npy"
enhamtx=np.load(fname)

genelist=pd.read_csv(samplename+"/pair_promoter.csv",sep=",",header=None)


grad=enhamtx.max(axis=1)
useidx=grad!=0

pairlist_prom=pd.read_csv(samplename+'/pair_promoter.csv',header=None)
pairlist_prom=pairlist_prom[[1,2,3]].to_numpy()

pairlist=pd.read_csv(samplename+'/pair_300000.csv',header=None)
pairlist=pairlist[[1,2,3]].to_numpy()

ATAC_use=np.load(samplename+"/ATAC_pred.npy")
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


np.save(samplename+"/ATAC_corr.npy",ATAC_corr)
ATAC_corr=np.load(samplename+"/ATAC_corrnorm.npy")

thv=np.percentile(grad[grad!=0], [40])
sthv=np.percentile(grad[grad!=0], [85])

promenhatag=np.zeros(grad.shape)
promenhatag[(grad>thv)&(ATAC_corr>0)]=2
promenhatag[(grad>sthv)&(ATAC_corr>0)]=4
promenhatag[(grad>thv)&(ATAC_corr<0)]=5
promenhatag[(grad>sthv)&(ATAC_corr<0)]=6
promenhatag[:,0]=3

cellclus=pd.read_csv(samplename+"/cellcluster.txt",sep=",",header=None)
cellclus=cellclus.to_numpy()
cellclus=cellclus[:,0]


#################################

strogene=genelist[(promenhatag==4).sum(axis=1)>0]
pd.DataFrame(strogene).to_csv('glioma/strogenelist.txt',sep="\t",header=False, index=False)

##############################
import scanpy as sc
import numpy as np
import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns


RNAmatrix=pd.read_csv("glioma/log1praw.csv",sep=",")
RNAmatrix=RNAmatrix.to_numpy()
rnanormfac=RNAmatrix.mean(axis=0)
RNAmatrix_n=RNAmatrix/rnanormfac
RNAmatrix_b=((RNAmatrix_n.transpose(1,0)-RNAmatrix_n.mean(axis=1))/RNAmatrix_n.std(axis=1)).transpose(1,0)

genelist=pd.read_csv('neuro/pair_300000.csv',header=None)
genelist=genelist[[0]].to_numpy()
genelist=genelist[:,0]
pd.DataFrame(genelist).to_csv('glioma/genelistfull.txt',sep="\t",header=False, index=False)

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



gsea=pd.read_csv('glioma/ssGSEA.csv',index_col=0)
tfidxr=np.load("glioma/enhamotif.npy")
tfidx=np.load("glioma/enhamotif_mix.npy")
tfactprom=pd.read_csv('glioma/motif_cell_activity_prom.csv',index_col=0)
tfactprom=tfactprom.to_numpy()

tfact=pd.read_csv('glioma/motif_cell_activity.csv',index_col=0)
tfact=tfact.to_numpy()
#tfactnorm=(tfact.transpose(1,0)/normfac).transpose(1,0)
tfactnorm=tfact
zeb1min,zeb1max=np.percentile(tfactnorm[(cellclus==0)&(cnotag==1),3], [10,90])
cremmin,cremmax=np.percentile(tfactnorm[(cellclus==0)&(cnotag==1),15], [10,90])

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
np.save('glioma/cellleiden.npy',cellleiden)
pd.DataFrame(cellleiden).to_csv("glioma/cellleiden.txt",sep="\t",header=False, index=False)
enha.obs["leiden"] = pd.Categorical(list(cellleiden[cellclus==0]))

enha.obs["ZEB1"] = tfactnorm[(cellclus==0)&(cnotag==1),3]
enha.obs["CREM"] = tfactnorm[(cellclus==0)&(cnotag==1),15]

sc.pl.umap(enha, frameon=False, color="ZEB1",save="glioma_enhancer3zeb1.pdf",vmin=zeb1min, vmax=zeb1max)
sc.pl.umap(enha, frameon=False, color="CREM",save="glioma_enhancer3crem.pdf",vmin=cremmin, vmax=cremmax)
sc.pl.umap(enha, frameon=False,legend_loc="on data", color="leiden",save="glioma_enhancer2.pdf")



from statannotations.Annotator import Annotator

enhatmpuse=enhamtx_n.transpose(0,2,1)[useidx==1]
enhatmpuse1=enhatmpuse[:,cellleiden==1]
enhatmpuse2=enhatmpuse[:,cellleiden==2]
enhatmpuse3=enhatmpuse[:,cellleiden==3]

enhamean1=enhatmpuse1.mean(axis=1)
enhamean2=enhatmpuse2.mean(axis=1)
enhamean3=enhatmpuse3.mean(axis=1)



for i in range(15):
   print(stats.ttest_ind(tfactnorm[cellleiden==1,i], tfactnorm[cellleiden==0,i]))

cnomin2,cnomax2=np.percentile(cno, [25,75])

cnotag2=(cno>cnomin2)&(cno<cnomax2)

from scipy.stats import spearmanr



enhamtxusetumor=enhamtxuse[:,cellclus==0].mean(axis=1)
enhamtxuseother=enhamtxuse[:,cellclus!=0].mean(axis=1)

tumorspe=enhamtxusetumor-enhamtxuseother

sc.tl.rank_genes_groups(enha, "leiden", method="wilcoxon")
sc.pl.rank_genes_groups(enha, n_genes=25, sharey=False, save="DAE_deeplift_pbmc_for_paper.png")


DAEmtx=pd.DataFrame(enha.uns["rank_genes_groups"]["names"])
pval=pd.DataFrame(enha.uns["rank_genes_groups"]["pvals_adj"])
fc=pd.DataFrame(enha.uns["rank_genes_groups"]["logfoldchanges"])

((pval.iloc[:,0]<0.0001)&(fc.iloc[:,0]>0)).sum()

((pval.iloc[:,0]>0.0001)&(fc.iloc[:,0]>0)).sum()

genelist=pd.read_csv('neuro/pair_300000.csv',header=None)
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
   tgtgene=genelist[(enhaclusgenemtx==1).sum(axis=1)>0]
   pd.DataFrame(tgtgene).to_csv("glioma/clus"+str(i)+"genenew.txt",sep="\t",header=False, index=False)

tgtgene=genelist[(enhaclusgenemtxall==1).sum(axis=(0,2))>0]
pd.DataFrame(tgtgene).to_csv("glioma/clusallgenenew.txt",sep="\t",header=False, index=False)

np.save('glioma/enhaclusgenemtxall.txt',enhaclusgenemtxall)

peaklistout=pd.read_csv('neuro/peaks_extend.bed',sep="\t",header=None)
clus0peak=peaklistout.iloc[np.unique(peakidmtx[enhaclusgenemtxall[0]==1]).astype(int),:]
clus1peak=peaklistout.iloc[np.unique(peakidmtx[enhaclusgenemtxall[1]==1]).astype(int),:]
clus2peak=peaklistout.iloc[np.unique(peakidmtx[enhaclusgenemtxall[2]==1]).astype(int),:]
clus3peak=peaklistout.iloc[np.unique(peakidmtx[enhaclusgenemtxall[3]==1]).astype(int),:]


clus0peak.to_csv("glioma/clus0peak.bed",sep="\t",header=False, index=False)
clus1peak.to_csv("glioma/clus1peak.bed",sep="\t",header=False, index=False)
clus2peak.to_csv("glioma/clus2peak.bed",sep="\t",header=False, index=False)
clus3peak.to_csv("glioma/clus3peak.bed",sep="\t",header=False, index=False)


sox2cre=peaklistout.iloc[np.unique(peakidmtx[(tfidxr[:,:,1]==1)&(enhaclusgenemtxall[0]==1)]).astype(int),:]
sox2noncre=peaklistout.iloc[np.unique(peakidmtx[(tfidxr[:,:,1]==1)&((promenhatag!=2)&(promenhatag!=4))]).astype(int),:]

sox2cre.to_csv("glioma/sox2cre.bed",sep="\t",header=False, index=False)
sox2noncre.to_csv("glioma/sox2noncre.bed",sep="\t",header=False, index=False)

##############

##code for bedtools
#bedtools intersect -wa -u -f 1 -b glioma/sox2cre.bed -a glioma/GSM6008250_SMNB19_SOX2_3_broadPeak.bed > glioma/tmp30b.bed
#bedtools intersect -wa -u -f 1  -b glioma/sox2noncre.bed -a glioma/GSM6008250_SMNB19_SOX2_3_broadPeak.bed > glioma/tmp31b.bed

############

tmp30b=pd.read_csv("glioma/tmp30b.bed",sep="\t",header=None)
tmp31b=pd.read_csv("glioma/tmp31b.bed",sep="\t",header=None)

tmp30b["Type"]="Tumor"
tmp31b["Type"]="Non-enhancer"

tmp30bout=tmp30b.iloc[:,[6,9]]
tmp31bout=tmp31b.iloc[:,[6,9]]

tmp30bout.columns=["fc","Type"]
tmp31bout.columns=["fc","Type"]

df_concat = pd.concat([tmp30bout,tmp31bout], axis=0)

from statannotations.Annotator import Annotator

x = "Type"
y = "fc"
order = ["Tumor","Non-enhancer"]
my_pal = {"Tumor": "blue", "Non-enhancer": "green"}
plt.clf()
plt.figure(figsize=(6.5,8))
plt.rcParams.update({'font.size': 18})
ax = sns.boxplot(data=df_concat, x=x, y=y, order=order,palette=my_pal,showfliers= False)
pairs = [("Tumor","Non-enhancer")]
annotator = Annotator(ax, pairs, data=df_concat, x=x, y=y, order=order)
annotator.configure(test='t-test_ind', text_format='star', loc='inside')
annotator.apply_and_annotate()
ax.set_xlabel("Enhancer", fontsize=18)
ax.set_ylabel("log2FC of SOX2 signal from background (Cut-n-Tag)", fontsize=18)
plt.title("SOX2 binding in glioblastoma",fontsize=24)
plt.xticks(fontsize=18,rotation=90)
plt.yticks(fontsize=18)
plt.subplots_adjust(left=0.20)
plt.subplots_adjust(bottom=0.40)
plt.savefig("glioma/sox2bind.pdf",format="pdf")



###############################################

enhamtxuse=enhamtx_n.transpose(0,2,1)/10
#enhamtxuse=enhamtxuse/normfac
enhamtxuse=enhamtxuse[(promenhatag==4)]
enhamtxusepow=pow(2, enhamtxuse.transpose(1,0))

enhaall = sc.AnnData(enhamtxusepow)
#enha = sc.AnnData(df2, df2.index.to_frame(), df2.columns.to_frame())
sc.tl.pca(enhaall, svd_solver="arpack")
sc.pp.neighbors(enhaall, n_neighbors=10, n_pcs=40)
sc.tl.umap(enhaall)
sc.tl.tsne(enhaall)

sc.pl.umap(enhaall, frameon=False, save="glioma_enhancer_all.png")
enhaall.obs["leiden"] = pd.Categorical(list(cellclus))
sc.pl.umap(enhaall, frameon=False,legend_loc="on data", color="leiden",save="glioma_enhancer_all2.pdf")
sc.pl.tsne(enhaall, frameon=False,legend_loc="on data", color="leiden",save="glioma_enhancer_all2t.pdf")

enhaall.obs["celltype"] = ""
enhaall.obs["celltype"][cellclus==0]="Tumor"
enhaall.obs["celltype"][cellclus==1]="OPC"
enhaall.obs["celltype"][cellclus==2]="Astrocyte"
enhaall.obs["celltype"][cellclus==3]="Myeloid"

sc.pl.umap(enhaall, frameon=False,legend_loc="on data", color="celltype",save="glioma_enhancer_all3.pdf")


sox2minall,sox2maxall=np.percentile(tfactnorm[:,1], [10,90])
enhaall.obs["SOX2"] = tfactnorm[:,1]
sc.pl.umap(enhaall, frameon=False, color="SOX2",save="glioma_enhancer3sox2_all.pdf",vmin=sox2minall, vmax=sox2maxall)

sox2rna=pd.read_csv(samplename+"/sox2exp.txt",sep="\t",header=None)
sox2rna=sox2rna.to_numpy()
sox2rna=sox2rna[:,0]
#emtrna=np.log2(emtrna)
sox2rnaminall,sox2rnamaxall=np.percentile(sox2rna, [10,90])
enhaall.obs["SOX2_RNA"] = sox2rna
sc.pl.umap(enhaall, frameon=False, color="SOX2_RNA",save="glioma_enhancer3sox2rna_all.pdf",vmin=sox2rnaminall, vmax=sox2rnamaxall)

sc.tl.rank_genes_groups(enhaall, "leiden", method="wilcoxon")
DAEmtxall=pd.DataFrame(enhaall.uns["rank_genes_groups"]["names"])
pvalall=pd.DataFrame(enhaall.uns["rank_genes_groups"]["pvals_adj"])
fcall=pd.DataFrame(enhaall.uns["rank_genes_groups"]["logfoldchanges"])


genelist=pd.read_csv('neuro/pair_300000.csv',header=None)
genelist=genelist[[0]].to_numpy()
genelist=genelist[:,0]
enhaclusgenemtxall_allenha=np.zeros((DAEmtxall.shape[1],promenhatag.shape[0],promenhatag.shape[1]))

for i in range(DAEmtxall.shape[1]):
   tmp=np.zeros(enhamtxuse.shape[0])
   tmp=tmp
   enhaname=DAEmtxall.iloc[:,i]
   usetag=enhaname[((pvalall.iloc[:,i]<0.05)&(fcall.iloc[:,i]>0))].to_numpy()
   usetag=usetag.astype(int)
   tmp[usetag]=1
   enhaclusgenemtx=np.zeros(promenhatag.shape)
   enhaclusgenemtx=enhaclusgenemtx
   enhaclusgenemtx[(promenhatag==4)]=tmp
   enhaclusgenemtxall_allenha[i]=enhaclusgenemtx
   tgtgene=genelist[(enhaclusgenemtx==1).sum(axis=1)>0]
   pd.DataFrame(tgtgene).to_csv("glioma/clus"+str(i)+"genenew_allenha.txt",sep="\t",header=False, index=False)

tgtgene=genelist[(enhaclusgenemtxall_allenha==1).sum(axis=(0,2))>0]
pd.DataFrame(tgtgene).to_csv("glioma/clusallgenenewall.txt",sep="\t",header=False, index=False)
np.save('glioma/enhaclusgenemtxall_allenha.npy',enhaclusgenemtxall_allenha)

peaklistout=pd.read_csv('neuro/peaks_extend.bed',sep="\t",header=None)
clus0peak=peaklistout.iloc[np.unique(peakidmtx[enhaclusgenemtxall[0]==1]).astype(int),:]
clus1peak=peaklistout.iloc[np.unique(peakidmtx[enhaclusgenemtxall[1]==1]).astype(int),:]
clus2peak=peaklistout.iloc[np.unique(peakidmtx[enhaclusgenemtxall[2]==1]).astype(int),:]
clus3peak=peaklistout.iloc[np.unique(peakidmtx[enhaclusgenemtxall[3]==1]).astype(int),:]

clus0peak.to_csv("glioma/clus0peak.bed",sep="\t",header=False, index=False)
clus1peak.to_csv("glioma/clus1peak.bed",sep="\t",header=False, index=False)
clus2peak.to_csv("glioma/clus2peak.bed",sep="\t",header=False, index=False)
clus3peak.to_csv("glioma/clus3peak.bed",sep="\t",header=False, index=False)



#################


gexall=sc.read_h5ad("glioma/gexm_allcell.h5ad")
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
