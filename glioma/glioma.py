import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
import sys
import math
import subprocess



samplename="glioma"

fname=samplename+"/Deeplift_full_ver2_mergenorm.npy"
#fname="pbmc_fold1/allgrad_ssep_max.npy"
enhamtx=np.load(fname)


genelist=pd.read_csv(samplename+"/pair_promoter.csv",sep=",",header=None)


clus1gene=pd.read_csv(samplename+"/cluster0_deg_glioma.txt",sep="\t",header=None)
clus2gene=pd.read_csv(samplename+"/cluster1_deg_glioma.txt",sep="\t",header=None)
clus3gene=pd.read_csv(samplename+"/cluster2_deg_glioma.txt",sep="\t",header=None)
clus4gene=pd.read_csv(samplename+"/cluster3_deg_glioma.txt",sep="\t",header=None)



geneclus=np.zeros((genelist.shape[0],4))

geneclus[:,0]=genelist[0].isin(clus1gene[0])
geneclus[:,1]=genelist[0].isin(clus2gene[0])
geneclus[:,2]=genelist[0].isin(clus3gene[0])
geneclus[:,3]=genelist[0].isin(clus4gene[0])



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


ATAC_corrt=np.zeros(grad.shape)
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
    ATAC_corrt[j,i]=np.corrcoef([ATAC_tmp[cellclus==0,i],enhamtx[j,cellclus==0,i]])[0,1]

np.save(samplename+"/ATAC_corrnorm.npy",ATAC_corr)
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

np.save("glioma/promenhatag.npy",promenhatag)

enhamtx1=enhamtx[geneclus[:,0]==1]
enhamtx2=enhamtx[geneclus[:,1]==1]
enhamtx3=enhamtx[geneclus[:,2]==1]
enhamtx4=enhamtx[geneclus[:,3]==1]


enhamtx1n=((enhamtx1.transpose(1,0,2)-enhamtx1.mean(axis=1))/enhamtx1.std(axis=1)).transpose(1,2,0)
enhamtx2n=((enhamtx2.transpose(1,0,2)-enhamtx2.mean(axis=1))/enhamtx2.std(axis=1)).transpose(1,2,0)
enhamtx3n=((enhamtx3.transpose(1,0,2)-enhamtx3.mean(axis=1))/enhamtx3.std(axis=1)).transpose(1,2,0)
enhamtx4n=((enhamtx4.transpose(1,0,2)-enhamtx4.mean(axis=1))/enhamtx4.std(axis=1)).transpose(1,2,0)


promenhatag1=promenhatag[geneclus[:,0]==1]
promenhatag2=promenhatag[geneclus[:,1]==1]
promenhatag3=promenhatag[geneclus[:,2]==1]
promenhatag4=promenhatag[geneclus[:,3]==1]



cellclus=pd.read_csv(samplename+"/cellcluster.txt",sep=",",header=None)
cellclus=cellclus.to_numpy()
cellclus=cellclus[:,0]

enhanorm=np.zeros((3,4,4))



tmp=enhamtx1n[(promenhatag1==3)]
for i in range(4):
   enhanorm[0,0,i]=tmp[:,cellclus==(i)].mean()


tmp=enhamtx2n[(promenhatag2==3)]
for i in range(4):
   enhanorm[0,1,i]=tmp[:,cellclus==(i)].mean()


tmp=enhamtx3n[(promenhatag3==3)]
for i in range(4):
   enhanorm[0,2,i]=tmp[:,cellclus==(i)].mean()


tmp=enhamtx4n[(promenhatag4==3)]
for i in range(4):
   enhanorm[0,3,i]=tmp[:,cellclus==(i)].mean()





tmp=enhamtx1n[(promenhatag1==2)]
for i in range(4):
   enhanorm[1,0,i]=tmp[:,cellclus==(i)].mean()


tmp=enhamtx2n[(promenhatag2==2)]
for i in range(4):
   enhanorm[1,1,i]=tmp[:,cellclus==(i)].mean()


tmp=enhamtx3n[(promenhatag3==2)]
for i in range(4):
   enhanorm[1,2,i]=tmp[:,cellclus==(i)].mean()


tmp=enhamtx4n[(promenhatag4==2)]
for i in range(4):
   enhanorm[1,3,i]=tmp[:,cellclus==(i)].mean()





tmp=enhamtx1n[(promenhatag1==4)]
for i in range(4):
   enhanorm[2,0,i]=tmp[:,cellclus==(i)].mean()


tmp=enhamtx2n[(promenhatag2==4)]
for i in range(4):
   enhanorm[2,1,i]=tmp[:,cellclus==(i)].mean()


tmp=enhamtx3n[(promenhatag3==4)]
for i in range(4):
   enhanorm[2,2,i]=tmp[:,cellclus==(i)].mean()


tmp=enhamtx4n[(promenhatag4==4)]
for i in range(4):
   enhanorm[2,3,i]=tmp[:,cellclus==(i)].mean()



import seaborn as sns

x_axis_labels = ["Glioma","OPC","Astrocyte","Myeloid"] # labels for x-axis
y_axis_labels = ["Promoter Glioma","Promoter OPC","Promoter Astrocyte","Promoter Myeloid","Weak enhancer Glioma","Weak enhancer OPC","Weak enhancer Astrocyte","Weak enhancer Myeloid","Strong enhancer Glioma","Strong enhancer OPC","Strong enhancer Astrocyte","Strong enhancer Myeloid"] # labels for y-axis

plt.figure()
sns.heatmap(enhanorm.reshape(12,4),xticklabels=x_axis_labels, yticklabels=y_axis_labels)
plt.subplots_adjust(left=0.35)
plt.title('Model prediction')
plt.savefig(samplename+"/model_glioma.png")




fantom1=np.load(samplename+"/alldonor1.npy")
fantom2=np.load(samplename+"/alldonor2.npy")
fantom3=np.load(samplename+"/alldonor3.npy")

for i in range(4):
    fantom1[:,:,i]=fantom1[:,:,i]/fantom1[:,:,i].sum()
    fantom2[:,:,i]=fantom2[:,:,i]/fantom2[:,:,i].sum()
    fantom3[:,:,i]=fantom3[:,:,i]/fantom3[:,:,i].sum()

fantomall=np.zeros((fantom1.shape[0],fantom1.shape[1],fantom1.shape[2]*3))
fantomall[:,:,[0,3,6,9]]=fantom1
fantomall[:,:,[1,4,7,10]]=fantom2
fantomall[:,:,[2,5,8,11]]=fantom3

geta=fantomall[fantomall>0].min()
fantomalllog=np.log(fantomall+geta)


fantom11=fantomall[geneclus[:,0]==1]
fantom12=fantomall[(geneclus[:,1]==1)|(geneclus[:,2]==1)]
fantom13=fantomall[geneclus[:,3]==1]
fantom14=fantomall[geneclus[:,4]==1]

fantom11n=((fantom11.transpose(2,0,1)-fantom11.mean(axis=2))/fantom11.std(axis=2)).transpose(1,2,0)
fantom12n=((fantom12.transpose(2,0,1)-fantom12.mean(axis=2))/fantom12.std(axis=2)).transpose(1,2,0)
fantom13n=((fantom13.transpose(2,0,1)-fantom13.mean(axis=2))/fantom13.std(axis=2)).transpose(1,2,0)
fantom14n=((fantom14.transpose(2,0,1)-fantom14.mean(axis=2))/fantom14.std(axis=2)).transpose(1,2,0)


fantomalln=((fantomalllog.transpose(2,0,1)-fantomalllog.mean(axis=2))/fantomalllog.std(axis=2)).transpose(1,2,0)
fantomcell=np.zeros(fantom1.shape)

fantomcell[:,:,0]=fantomalln[:,:,0:3].mean(axis=2)
fantomcell[:,:,1]=fantomalln[:,:,3:6].mean(axis=2)
fantomcell[:,:,2]=fantomalln[:,:,6:9].mean(axis=2)
fantomcell[:,:,3]=fantomalln[:,:,9:12].mean(axis=2)

fantomcellstd=fantomcell.std(axis=2)


fantomnorm=np.zeros((3,4,12))


fantomnorm[0,0,:]=np.nanmean(fantom11n[(promenhatag1==3)],axis=0)
fantomnorm[0,1,:]=np.nanmean(fantom12n[(promenhatag2==3)],axis=0)
fantomnorm[0,2,:]=np.nanmean(fantom13n[(promenhatag3==3)],axis=0)
fantomnorm[0,3,:]=np.nanmean(fantom14n[(promenhatag4==3)],axis=0)

fantomnorm[1,0,:]=np.nanmean(fantom11n[(promenhatag1==2)],axis=0)
fantomnorm[1,1,:]=np.nanmean(fantom12n[(promenhatag2==2)],axis=0)
fantomnorm[1,2,:]=np.nanmean(fantom13n[(promenhatag3==2)],axis=0)
fantomnorm[1,3,:]=np.nanmean(fantom14n[(promenhatag4==2)],axis=0)

fantomnorm[2,0,:]=np.nanmean(fantom11n[(promenhatag1==4)],axis=0)
fantomnorm[2,1,:]=np.nanmean(fantom12n[(promenhatag2==4)],axis=0)
fantomnorm[2,2,:]=np.nanmean(fantom13n[(promenhatag3==4)],axis=0)
fantomnorm[2,3,:]=np.nanmean(fantom14n[(promenhatag4==4)],axis=0)



#################################

##############################
import scanpy as sc
import numpy as np
import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import yeojohnson
from pygam import s, LinearGAM
import anndata
from scipy.spatial import distance
import seaborn as sns
from scipy.stats import spearmanr
import random
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import umap
from scipy import stats
import os

RNAmatrix=pd.read_csv("glioma/log1praw.csv",sep=",")
RNAmatrix=RNAmatrix.to_numpy()
rnanormfac=RNAmatrix.mean(axis=0)
RNAmatrix_n=RNAmatrix/rnanormfac
RNAmatrix_b=((RNAmatrix_n.transpose(1,0)-RNAmatrix_n.mean(axis=1))/RNAmatrix_n.std(axis=1)).transpose(1,0)

genelist=pd.read_csv('neuro/pair_300000.csv',header=None)
genelist=genelist[[0]].to_numpy()
genelist=genelist[:,0]

np.abs(enhamtx[:,cellleiden==0,:]).sum(axis=(0,2)).mean()
np.abs(enhamtx[:,cellleiden==1,:]).sum(axis=(0,2)).mean()
np.abs(enhamtx[:,cellleiden==2,:]).sum(axis=(0,2)).mean()
np.abs(enhamtx[:,cellleiden==3,:]).sum(axis=(0,2)).mean()

(enhamtx[:,cellleiden==0,:]>0).sum(axis=(0,2)).mean()
(enhamtx[:,cellleiden==1,:]>0).sum(axis=(0,2)).mean()
(enhamtx[:,cellleiden==2,:]>0).sum(axis=(0,2)).mean()
(enhamtx[:,cellleiden==3,:]>0).sum(axis=(0,2)).mean()

(enhamtx[:,cellleiden==0,:]<0).sum(axis=(0,2)).mean()
(enhamtx[:,cellleiden==1,:]<0).sum(axis=(0,2)).mean()
(enhamtx[:,cellleiden==2,:]<0).sum(axis=(0,2)).mean()
(enhamtx[:,cellleiden==3,:]<0).sum(axis=(0,2)).mean()


cno[cellleiden==0].mean()
cno[cellleiden==1].mean()
cno[cellleiden==2].mean()
cno[cellleiden==3].mean()


np.percentile(enhavar, [90])

enhamtx_n=enhamtx.copy()
#normfac=np.abs(enhamtx).sum(axis=(0,2))
#normfac1=(np.abs(enhamtx)*(enhamtx>0)).sum(axis=(0,2))/(enhamtx>0).sum(axis=(0,2))
#normfac2=(np.abs(enhamtx)*(enhamtx<0)).sum(axis=(0,2))/(enhamtx<0).sum(axis=(0,2))
enhatmp=enhamtx.transpose(0,2,1)[(promenhatag==4)]
enhatmp=(enhatmp-enhatmp.mean(axis=0))/enhatmp.std(axis=0)
enhamtx_n.transpose(0,2,1)[(promenhatag==4)]=enhatmp
#enhamtx_n=(enhamtx_n.transpose(0,2,1)/normfac).transpose(0,2,1)
#enhamtx_n=((enhamtx*(enhamtx>0)).transpose(0,2,1)/normfac1).transpose(0,2,1)+((enhamtx*(enhamtx<0)).transpose(0,2,1)/normfac2).transpose(0,2,1)

enhamtxuse=enhamtx.transpose(0,2,1)
enhamtxuse=enhamtxuse[(promenhatag==4)]


enhavar=enhamtxuse[:,cellclus==0].var(axis=1)

varmtx=np.zeros((enhamtx.shape[0],enhamtx.shape[2]))
varmtx[promenhatag==4]=enhavar
varmtxidx=(varmtx>=np.percentile(enhavar, [90]))

np.save("glioma/varmtxidx.npy",varmtxidx)

vargene=genelist[(varmtxidx.sum(axis=1)>0)]

pd.DataFrame(vargene).to_csv('glioma/vargene.txt',sep="\t",header=False, index=False)

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

((cno<=cnomin)&(cellleiden==3)).sum()
((cno>cnomin)&(cellleiden==3)).sum()

enhamtxall=enhamtx_n.transpose(0,2,1)
enhamtxall=enhamtxall[(promenhatag==4)]
stronormfac=np.abs(enhamtxall).sum(axis=0)#/normfac
#stronormfac=cno
stronormfacmin,stronormfacmax=np.percentile(stronormfac[cellclus==0], [10,90])
normfacst=np.abs(enhamtxall).sum(axis=0)

normfac=np.abs(enhamtx).sum(axis=(0,2))

cnotag=(cno>cnomin)&(cno<cnomax)

enhamtxall=enhamtx_n.transpose(0,2,1)[(promenhatag==4)]
stronormfac=(enhamtxall).sum(axis=0)
stronormfacmin,stronormfacmax=np.percentile(stronormfac[cellclus==0], [10,90])



gsea=pd.read_csv('glioma/ssGSEA.csv',index_col=0)
emtgseamin,emtgseamax=np.percentile(gsea["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"], [10,90])
nfkbgseamin,nfkbgseamax=np.percentile(gsea["HALLMARK_TNFA_SIGNALING_VIA_NFKB"], [10,90])

krasdngseamin,krasdngseamax=np.percentile(gsea["HALLMARK_KRAS_SIGNALING_DN"], [10,90])
infgseamin,infgseamax=np.percentile(gsea["HALLMARK_INFLAMMATORY_RESPONSE"], [10,90])

uvdngseamin,uvdngseamax=np.percentile(gsea["HALLMARK_UV_RESPONSE_DN"], [10,90])
hhgseamin,hhgseamax=np.percentile(gsea["HALLMARK_HEDGEHOG_SIGNALING"], [10,90])

tfidxr=np.load("glioma/enhamotif.npy")
tfidx=np.load("glioma/enhamotif_mix.npy")


tfactprom=pd.read_csv('glioma/motif_cell_activity_prom.csv',index_col=0)
tfactprom=tfactprom.to_numpy()

tfact=pd.read_csv('glioma/motif_cell_activity.csv',index_col=0)
tfact=tfact.to_numpy()
#tfactnorm=(tfact.transpose(1,0)/normfac).transpose(1,0)
tfactnorm=tfact
nficmin,nficmax=np.percentile(tfactnorm[(cellclus==0)&(cnotag==1),0], [10,90])
sox2min,sox2max=np.percentile(tfactnorm[(cellclus==0)&(cnotag==1),1], [10,90])
smad3min,smad3max=np.percentile(tfactnorm[(cellclus==0)&(cnotag==1),2], [10,90])
zeb1min,zeb1max=np.percentile(tfactnorm[(cellclus==0)&(cnotag==1),3], [10,90])
ap1min,ap1max=np.percentile(tfactnorm[(cellclus==0)&(cnotag==1),4], [10,90])
smadap1min,smadap1max=np.percentile(tfactnorm[(cellclus==0)&(cnotag==1),5], [10,90])
nfkbmin,nfkbmax=np.percentile(tfactnorm[(cellclus==0)&(cnotag==1),6], [10,90])
egr1min,egr1max=np.percentile(tfactnorm[(cellclus==0)&(cnotag==1),7], [10,90])
nfixmin,nfixmax=np.percentile(tfactnorm[(cellclus==0)&(cnotag==1),8], [10,90])
nfictlx1min,nfictlx1max=np.percentile(tfactnorm[(cellclus==0)&(cnotag==1),9], [10,90])
smad3zeb1min,smad3zeb1max=np.percentile(tfactnorm[(cellclus==0)&(cnotag==1),10], [10,90])
sox2zeb1min,sox2zeb1max=np.percentile(tfactnorm[(cellclus==0)&(cnotag==1),11], [10,90])
sox2ap1min,sox2ap1max=np.percentile(tfactnorm[(cellclus==0)&(cnotag==1),12], [10,90])
ascl1min,ascl1max=np.percentile(tfactnorm[(cellclus==0)&(cnotag==1),13], [10,90])
nfibmin,nfibmax=np.percentile(tfactnorm[(cellclus==0)&(cnotag==1),14], [10,90])
cremmin,cremmax=np.percentile(tfactnorm[(cellclus==0)&(cnotag==1),15], [10,90])

for i in range(15):
   tmp1=tfidx[(promenhatag==4),i]
   tmp2=tfidx[(promenhatag==4),4]
   tmp3=(((tmp1==1)&(tmp2==1)).sum()*((tmp1==0)&(tmp2==0)).sum())/(((tmp1==0)&(tmp2==1)).sum()*((tmp1==1)&(tmp2==0)).sum())
   print(tmp3)


for i in range(15):
   tmp1=tfidx[(promenhatag==4)&(emtmtx==1),i]
   tmp2=tfidx[(promenhatag==4)&(emtmtx==1),1]
   tmp3=tfidx[(promenhatag==4),i]
   tmp4=tfidx[(promenhatag==4),1]
   tmp5=(((tmp1==1)&(tmp2==1)).sum()/((tmp2==1)).sum())/(((tmp3==1)&(tmp4==1)).sum()/((tmp4==1)).sum())
   print(tmp5)

for i in range(15):
   tmp1=tfidx[(promenhatag==4)&(nfkbmtx==1),i]
   tmp3=tfidx[(promenhatag==4),i]
   tmp5=tmp1.mean()/tmp3.mean()
   print(tmp5)

tmp1=tfidx[(promenhatag==4),3]
tmp2=tfidx[(promenhatag==4),6]
tmp3=(((tmp1==1)&(tmp2==1)).sum()*((tmp1==0)&(tmp2==0)).sum())/(((tmp1==0)&(tmp2==1)).sum()*((tmp1==1)&(tmp2==0)).sum())
print(tmp3)

emtgene=pd.read_csv(samplename+"/emtgene.txt",sep="\t",header=None)
emttag=np.isin(genelist,emtgene)
emtmtx=np.ones(promenhatag.shape)
emtmtx=(emtmtx.transpose(1,0)*(emttag.astype(int))).transpose(1,0)

krasgene=pd.read_csv(samplename+"/krasgene.txt",sep="\t",header=None)
krastag=np.isin(genelist,krasgene)
krasmtx=np.ones(promenhatag.shape)
krasmtx=(krasmtx.transpose(1,0)*(krastag.astype(int))).transpose(1,0)

tgfbgene=pd.read_csv(samplename+"/tgfbgene.txt",sep="\t",header=None)
tgfbtag=np.isin(genelist,tgfbgene)
tgfbmtx=np.ones(promenhatag.shape)
tgfbmtx=(tgfbmtx.transpose(1,0)*(tgfbtag.astype(int))).transpose(1,0)

tmptag=((tfidx[:,:,2]==1).sum(axis=1)>0)
tmptagmtx=np.ones(promenhatag.shape)
tmptagmtx=(tmptagmtx.transpose(1,0)*(tmptag.astype(int))).transpose(1,0)


#enhamtxuseemt=enhamtx.transpose(0,2,1)
#enhamtxuseemt=enhamtxuseemt/normfacst
enhamtxuseemt=enhamtx_n.transpose(0,2,1)[((promenhatag==4))&(emtmtx==1)]
emtsignal=enhamtxuseemt.sum(axis=0)
emtmin,emtmax=np.percentile(emtsignal[(cellclus==0)&(cnotag==1)0], [10,90])

#enhamtxusekras=enhamtx.transpose(0,2,1)
#enhamtxusekras=enhamtxusekras/normfacst
enhamtxusekras=enhamtx_n.transpose(0,2,1)[((promenhatag==4))&(krasmtx==1)]
krassignal=enhamtxusekras.sum(axis=0)
krasmin,krasmax=np.percentile(krassignal[(cellclus==0)&(cnotag==1)], [10,90])

#enhamtxusetgfb=enhamtx.transpose(0,2,1)
#enhamtxusetgfb=enhamtxusetgfb/normfacst
enhamtxusetgfb=enhamtx_n.transpose(0,2,1)[((promenhatag==4))&(tgfbmtx==1)]
#enhamtxusetgfb=enhamtx_n.transpose(0,2,1)[(promenhatag==4)&(tgfbmtx==1)&(tfidx[:,:,2]==1)&(tmptagmtx==1)]
tgfbsignal=enhamtxusetgfb.sum(axis=0)
tgfbmin,tgfbmax=np.percentile(tgfbsignal[(cellclus==0)&(cnotag==1)], [10,90])

#tfidx[(promenhatag==4)&(tgfbmtx==1)]
#tfidx[(promenhatag==4)&(tgfbmtx==1)].mean(axis=0)/tfidx[(promenhatag==4)&(tgfbmtx==0)].mean(axis=0)
#nfmin,nfmax=np.percentile(normfac[cellclus==0], [10,90])


hypogene=pd.read_csv(samplename+"/hypogene.txt",sep="\t",header=None)
hypotag=np.isin(genelist,hypogene)
hypomtx=np.ones(promenhatag.shape)
hypomtx=(hypomtx.transpose(1,0)*(hypotag.astype(int))).transpose(1,0)
#enhamtxusehypo=enhamtx.transpose(0,2,1)
#enhamtxusehypo=enhamtxusehypo/normfacst
enhamtxusehypo=enhamtx_n.transpose(0,2,1)[((promenhatag==4))&(hypomtx==1)]
hyposignal=enhamtxusehypo.sum(axis=0)
hypomin,hypomax=np.percentile(hyposignal[(cellclus==0)&(cnotag==1)], [10,90])


hghggene=pd.read_csv(samplename+"/hghggene.txt",sep="\t",header=None)
hghgtag=np.isin(genelist,hghggene)
hghgmtx=np.ones(promenhatag.shape)
hghgmtx=(hghgmtx.transpose(1,0)*(hghgtag.astype(int))).transpose(1,0)
#enhamtxusehghg=enhamtx.transpose(0,2,1)
#enhamtxusehghg=enhamtxusehghg/normfacst
enhamtxusehghg=enhamtx_n.transpose(0,2,1)[((promenhatag==4))&(hghgmtx==1)]
hghgsignal=enhamtxusehghg.sum(axis=0)
hghgmin,hghgmax=np.percentile(hghgsignal[(cellclus==0)&(cnotag==1)0], [10,90])


nfkbgene=pd.read_csv(samplename+"/nfkbgene.txt",sep="\t",header=None)
nfkbtag=np.isin(genelist,nfkbgene)
nfkbmtx=np.ones(promenhatag.shape)
nfkbmtx=(nfkbmtx.transpose(1,0)*(nfkbtag.astype(int))).transpose(1,0)
#enhamtxusenfkb=enhamtx.transpose(0,2,1)
#nhamtxusenfkb=enhamtxusenfkb/normfacst
enhamtxusenfkb=enhamtx_n.transpose(0,2,1)[((promenhatag==4))&(nfkbmtx==1)]
nfkbsignal=enhamtxusenfkb.sum(axis=0)
nfkbactmin,nfkbactmax=np.percentile(nfkbsignal[(cellclus==0)&(cnotag==1)], [10,90])


emtsignalnorm=(emtsignal-emtsignal.mean())/emtsignal.std()


emtrna=RNAmatrix_b[emttag].mean(axis=0)
emtrnamin,emtrnamax=np.percentile(emtrna[(cellclus==0)&(cnotag==1)], [10,90])

emtrna=pd.read_csv(samplename+"/emtexp.txt",sep="\t",header=None)
emtrna=emtrna.to_numpy()
emtrna=emtrna[:,0]
#emtrna=np.log2(emtrna)
emtrnamin,emtrnamax=np.percentile(emtrna[(cellclus==0)&(cnotag==1)], [10,90])

nfkbrna=pd.read_csv(samplename+"/nfkbexp.txt",sep="\t",header=None)
nfkbrna=nfkbrna.to_numpy()
nfkbrna=nfkbrna[:,0]
#emtrna=np.log2(emtrna)
nfkbrnamin,nfkbrnamax=np.percentile(nfkbrna[(cellclus==0)&(cnotag==1)], [10,90])

[:,cellclus==0]



pca = PCA()
pca.fit(enhamtxuse.transpose(1,0))
embedding_comb_pca = pca.transform(enhamtxuse[:,cno>1000].transpose(1,0))
embedding_comb = umap.UMAP().fit_transform(embedding_comb_pca[:,0:30])

embedding_comb = umap.UMAP().fit_transform(enhamtxuse.transpose(1,0))
kmeans_embed = KMeans(n_clusters=3,max_iter=30, init="random")
cluster_kmeans = kmeans_embed.fit_predict(enhatumornorm.transpose(1,0))


fig, ax = plt.subplots()
ax.scatter(embedding_comb[:,0], embedding_comb[:,1],s=1, c=(normfac).astype(float))
plt.xlabel('UMAP1')
plt.ylabel("UMAP2")
plt.legend(fontsize=20)
fname=samplename+"/umap_kmeans_promenha_new.png"
plt.savefig(fname)
plt.show()


fig, ax = plt.subplots()
ax.scatter(embedding_comb_pca[:,0], embedding_comb_pca[:,1],s=1, c=(normfac[cellclus==0]).astype(float))
plt.xlabel('UMAP1')
plt.ylabel("UMAP2")
plt.legend(fontsize=20)
fname=samplename+"/umap_kmeans_promenha_new.png"
plt.savefig(fname)
plt.show()



fig, ax = plt.subplots()
ax.scatter(embedding_comb_pca[:,1], normfac,s=1, c=(np.log(cno)).astype(float))
fname=samplename+"/umap_kmeans_promenha_normfac.png"
plt.savefig(fname)
plt.show()



enhamtxuse=enhamtx.transpose(0,2,1)
enhamtxuse=enhamtxuse[(promenhatag==4)]
embedding_comb = umap.UMAP().fit_transform(enhamtxuse[:,cellclus==0].transpose(1,0))

fig, ax = plt.subplots()
ax.scatter(embedding_comb[:,0], embedding_comb[:,1],s=1, c=(emtsignal[cellclus==0]).astype(float),vmin=emtmin, vmax=emtmax)
plt.xlabel('UMAP1')
plt.ylabel("UMAP2")
plt.legend(fontsize=20)
fname=samplename+"/umap_kmeans_promenha_newemt.png"
plt.savefig(fname)
plt.show()

fig, ax = plt.subplots()
ax.scatter(embedding_comb[:,0], embedding_comb[:,1],s=1, c=(krassignal[cellclus==0]).astype(float),vmin=krasmin, vmax=krasmax)
plt.xlabel('UMAP1')
plt.ylabel("UMAP2")
plt.legend(fontsize=20)
fname=samplename+"/umap_kmeans_promenha_newkras.png"
plt.savefig(fname)
plt.show()

fig, ax = plt.subplots()
ax.scatter(embedding_comb[:,0], embedding_comb[:,1],s=1, c=(normfac[cellclus==0]).astype(float),vmin=nfmin, vmax=nfmax)
plt.xlabel('UMAP1')
plt.ylabel("UMAP2")
plt.legend(fontsize=20)
fname=samplename+"/umap_kmeans_promenha_new.png"
plt.savefig(fname)
plt.show()

fig, ax = plt.subplots()
ax.scatter(embedding_comb[:,0], embedding_comb[:,1],s=1, c=(cellleiden[cellclus==0]).astype(float))
plt.xlabel('UMAP1')
plt.ylabel("UMAP2")
plt.legend(fontsize=20)
fname=samplename+"/umap_kmeans_promenha_newclus.png"
plt.savefig(fname)
plt.show()


############################
enhagsea = sc.AnnData(gsea)
sc.tl.pca(enhagsea, svd_solver="arpack")
sc.pp.neighbors(enhagsea, n_neighbors=10, n_pcs=40)
sc.tl.umap(enhagsea)

sc.pl.umap(enhagsea, frameon=False, save="glioma_enhagsea.png")



###############
enhamtxuse=enhamtx_n.transpose(0,2,1)/10
#enhamtxuse=enhamtxuse/normfac
enhamtxuse=enhamtxuse[(promenhatag==4)]
enhamtxusepow=pow(2, enhamtxuse.transpose(1,0))

enha = sc.AnnData(enhamtxusepow[(cellclus==0)&(cnotag==1)])
#enha = sc.AnnData(df2, df2.index.to_frame(), df2.columns.to_frame())
sc.tl.pca(enha, svd_solver="arpack")
sc.pp.neighbors(enha, n_neighbors=10, n_pcs=40)
sc.tl.umap(enha)
sc.tl.tsne(enha)

sc.pl.umap(enha, frameon=False, save="glioma_enhancer.png")


sc.tl.leiden(
    enha,
    resolution=0.2
)

sc.pl.umap(enha, frameon=False,legend_loc="on data", color="leiden",save="glioma_enhancer2.png")
sc.pl.tsne(enha, frameon=False,legend_loc="on data", color="leiden",save="glioma_enhancer2t.png")


cellleiden=pd.read_csv("glioma/cellleiden.txt",header=None)
cellleiden=cellleiden.to_numpy()
cellleiden=cellleiden[:,0]
enha.obs["leiden"] = pd.Categorical(list(cellleiden[cellclus==0]))

enha.obs["emt"] = emtsignal[(cellclus==0)&(cnotag==1)]
enha.obs["kras"] = krassignal[(cellclus==0)&(cnotag==1)]
enha.obs["tgfb"] = tgfbsignal[(cellclus==0)&(cnotag==1)]
enha.obs["hypo"] = hyposignal[(cellclus==0)&(cnotag==1)]
enha.obs["hghg"] = hghgsignal[(cellclus==0)&(cnotag==1)]
enha.obs["NFKBact"] = nfkbsignal[(cellclus==0)&(cnotag==1)]
enha.obs["NFIC"] = tfactnorm[(cellclus==0)&(cnotag==1),0]
enha.obs["SOX2"] = tfactnorm[(cellclus==0)&(cnotag==1),1]
enha.obs["SMAD3"] = tfactnorm[(cellclus==0)&(cnotag==1),2]
enha.obs["ZEB1"] = tfactnorm[(cellclus==0)&(cnotag==1),3]
enha.obs["AP1"] = tfactnorm[(cellclus==0)&(cnotag==1),4]
enha.obs["SMAD3_AP1"] = tfactnorm[(cellclus==0)&(cnotag==1),5]
enha.obs["NFKB"] = tfactnorm[(cellclus==0)&(cnotag==1),6]
enha.obs["EGR1"] = tfactnorm[(cellclus==0)&(cnotag==1),7]
enha.obs["NFIX"] = tfactnorm[(cellclus==0)&(cnotag==1),8]
enha.obs["NFIC_TLX1"] = tfactnorm[(cellclus==0)&(cnotag==1),9]
enha.obs["SMAD3_ZEB1"] = tfactnorm[(cellclus==0)&(cnotag==1),10]
enha.obs["SOX2_ZEB1"] = tfactnorm[(cellclus==0)&(cnotag==1),11]
enha.obs["SOX2_AP1"] = tfactnorm[(cellclus==0)&(cnotag==1),12]
enha.obs["ASCL1"] = tfactnorm[(cellclus==0)&(cnotag==1),13]
enha.obs["NFIB"] = tfactnorm[(cellclus==0)&(cnotag==1),14]
enha.obs["CREM"] = tfactnorm[(cellclus==0)&(cnotag==1),15]

enha.obs["EMT_RNA"] = emtrna[(cellclus==0)&(cnotag==1)]
enha.obs["NFKB_RNA"] = nfkbrna[(cellclus==0)&(cnotag==1)]
enha.obs["Strongnormfac"] = stronormfac[(cellclus==0)&(cnotag==1)]

enha.obs["EMT_GSEA"] = gsea["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"].to_numpy()
enha.obs["NFKB_GSEA"] = gsea["HALLMARK_TNFA_SIGNALING_VIA_NFKB"].to_numpy()
enha.obs["KRAS_DOWN_GSEA"] = gsea["HALLMARK_KRAS_SIGNALING_DN"].to_numpy()
enha.obs["INFLAMMATORY_GSEA"] = gsea["HALLMARK_INFLAMMATORY_RESPONSE"].to_numpy()
enha.obs["UVresponse_down_GSEA"] = gsea["HALLMARK_UV_RESPONSE_DN"].to_numpy()
enha.obs["Hedgehog_GSEA"] = gsea["HALLMARK_HEDGEHOG_SIGNALING"].to_numpy()



cellleiden=np.zeros((cellclus.shape[0]))
cellleiden=cellleiden-1
cellleiden[(cellclus==0)&(cnotag==1)]=enha.obs["leiden"].to_numpy()

np.save('glioma/cellleiden.npy',cellleiden)
pd.DataFrame(cellleiden).to_csv("glioma/cellleiden.txt",sep="\t",header=False, index=False)

sc.pl.umap(enha, frameon=False, color="emt",save="glioma_enhancer3emt.png",vmin=emtmin, vmax=emtmax)
sc.pl.umap(enha, frameon=False, color="kras",save="glioma_enhancer3kras.png",vmin=krasmin, vmax=krasmax)
sc.pl.umap(enha, frameon=False, color="tgfb",save="glioma_enhancer3tgfb.png",vmin=tgfbmin, vmax=tgfbmax)
sc.pl.umap(enha, frameon=False, color="hypo",save="glioma_enhancer3hypo.png",vmin=hypomin, vmax=hypomax)
sc.pl.umap(enha, frameon=False, color="hghg",save="glioma_enhancer3hghg.png",vmin=hghgmin, vmax=hghgmax)
sc.pl.umap(enha, frameon=False, color="NFKBact",save="glioma_enhancer3nfkbact.png",vmin=nfkbactmin, vmax=nfkbactmax)
sc.pl.umap(enha, frameon=False, color="NFIC",save="glioma_enhancer3nfic.png",vmin=nficmin, vmax=nficmax)
sc.pl.umap(enha, frameon=False, color="SOX2",save="glioma_enhancer3sox2.png",vmin=sox2min, vmax=sox2max)
sc.pl.umap(enha, frameon=False, color="SMAD3",save="glioma_enhancer3smad3.png",vmin=smad3min, vmax=smad3max)
sc.pl.umap(enha, frameon=False, color="ZEB1",save="glioma_enhancer3zeb1.png",vmin=zeb1min, vmax=zeb1max)
sc.pl.umap(enha, frameon=False, color="AP1",save="glioma_enhancer3ap1.png",vmin=ap1min, vmax=ap1max)
sc.pl.umap(enha, frameon=False, color="SMAD3_AP1",save="glioma_enhancer3smad3ap1.png",vmin=smadap1min, vmax=smadap1max)
sc.pl.umap(enha, frameon=False, color="NFKB",save="glioma_enhancer3nfkb.png",vmin=nfkbmin, vmax=nfkbmax)
sc.pl.umap(enha, frameon=False, color="EGR1",save="glioma_enhancer3egr1.png",vmin=egr1min, vmax=egr1max)
sc.pl.umap(enha, frameon=False, color="NFIX",save="glioma_enhancer3nfix.png",vmin=nfixmin, vmax=nfixmax)
sc.pl.umap(enha, frameon=False, color="NFIC_TLX1",save="glioma_enhancer3nfictlx1.png",vmin=nfictlx1min, vmax=nfictlx1max)
sc.pl.umap(enha, frameon=False, color="SMAD3_ZEB1",save="glioma_enhancer3smad3zeb1.png",vmin=smad3zeb1min, vmax=smad3zeb1max)
sc.pl.umap(enha, frameon=False, color="SOX2_ZEB1",save="glioma_enhancer3sox2zeb1.png",vmin=sox2zeb1min, vmax=sox2zeb1max)
sc.pl.umap(enha, frameon=False, color="SOX2_AP1",save="glioma_enhancer3sox2ap1.png",vmin=sox2ap1min, vmax=sox2ap1max)
sc.pl.umap(enha, frameon=False, color="ASCL1",save="glioma_enhancer3ascl1.png",vmin=ascl1min, vmax=ascl1max)
sc.pl.umap(enha, frameon=False, color="NFIB",save="glioma_enhancer3nfib.png",vmin=nfibmin, vmax=nfibmax)
sc.pl.umap(enha, frameon=False, color="CREM",save="glioma_enhancer3crem.png",vmin=cremmin, vmax=cremmax)


sc.pl.umap(enha, frameon=False, color="Strongnormfac",save="glioma_enhancer3Strongnormfac.png",vmin=stronormfacmin, vmax=stronormfacmax)

sc.pl.umap(enha, frameon=False, color="EMT_RNA",save="glioma_enhancer3emtrna.png",vmin=emtrnamin, vmax=emtrnamax)
sc.pl.umap(enha, frameon=False, color="NFKB_RNA",save="glioma_enhancer3nfkbrna.png",vmin=emtrnamin, vmax=emtrnamax)

sc.pl.umap(enha, frameon=False, color="EMT_GSEA",save="glioma_enhancer3emtgsea.png",vmin=emtgseamin, vmax=emtgseamax)
sc.pl.umap(enha, frameon=False, color="NFKB_GSEA",save="glioma_enhancer3nfkbgsea.png",vmin=nfkbgseamin, vmax=nfkbgseamax)
sc.pl.umap(enha, frameon=False, color="KRAS_DOWN_GSEA",save="glioma_enhancer3krasdngsea.png",vmin=krasdngseamin, vmax=krasdngseamax)
sc.pl.umap(enha, frameon=False, color="INFLAMMATORY_GSEA",save="glioma_enhancer3infgsea.png",vmin=infgseamin, vmax=infgseamax)
sc.pl.umap(enha, frameon=False, color="UVresponse_down_GSEA",save="glioma_enhancer3uvdngsea.png",vmin=uvdngseamin, vmax=uvdngseamax)
sc.pl.umap(enha, frameon=False, color="Hedgehog_GSEA",save="glioma_enhancer3hhgsea.png",vmin=hhgseamin, vmax=hhgseamax)


sc.pl.umap(enha, frameon=False, color="ZEB1",save="glioma_enhancer3zeb1.pdf",vmin=zeb1min, vmax=zeb1max)
sc.pl.umap(enha, frameon=False, color="CREM",save="glioma_enhancer3crem.pdf",vmin=cremmin, vmax=cremmax)



emtgseamin,emtgseamax=np.percentile(gsea["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"], [5,95])
nfkbgseamin,nfkbgseamax=np.percentile(gsea["HALLMARK_TNFA_SIGNALING_VIA_NFKB"], [5,95])


sc.pl.umap(enha, frameon=False, color="EMT_GSEA",save="glioma_enhancer3emtgsea.pdf",vmin=emtgseamin, vmax=emtgseamax)
sc.pl.umap(enha, frameon=False, color="NFKB_GSEA",save="glioma_enhancer3nfkbgsea.pdf",vmin=nfkbgseamin, vmax=nfkbgseamax)

sc.pl.umap(enha, frameon=False,legend_loc="on data", color="leiden",save="glioma_enhancer2.pdf")

gsea["cellclus"]="Cluster 0"
gsea["cellclus"].iloc[cellleiden[cellleiden>=0]==1]="Cluster 1"



from statannotations.Annotator import Annotator

x = "cellclus"
y = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
order = ["Cluster 0","Cluster 1"]
my_pal = {"Cluster 0": "blue", "Cluster 1": "green"}
plt.clf()
plt.figure(figsize=(6.5,8))
plt.rcParams.update({'font.size': 18})
ax = sns.violinplot(data=gsea, x=x, y=y, order=order,palette=my_pal,showfliers= False)
pairs = [('Cluster 0', 'Cluster 1')]
annotator = Annotator(ax, pairs, data=gsea, x=x, y=y, order=order)
annotator.configure(test='t-test_ind', text_format='star', loc='inside')
annotator.apply_and_annotate()
ax.set_xlabel("Cluster", fontsize=18)
ax.set_ylabel("ssGSEA score", fontsize=18)
plt.ylim(-0.05,0.48)
plt.title("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",fontsize=24)
plt.xticks(fontsize=18,rotation=90)
plt.yticks(fontsize=18)
plt.subplots_adjust(left=0.20)
plt.subplots_adjust(bottom=0.40)
plt.savefig("glioma/HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_varplot.pdf",format="pdf")



x = "cellclus"
y = "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
order = ["Cluster 0","Cluster 1"]
my_pal = {"Cluster 0": "blue", "Cluster 1": "green"}
plt.clf()
plt.figure(figsize=(6.5,7))
plt.rcParams.update({'font.size': 18})
ax = sns.violinplot(data=gsea, x=x, y=y, order=order,palette=my_pal,showfliers= False)
pairs = [('Cluster 0', 'Cluster 1')]
annotator = Annotator(ax, pairs, data=gsea, x=x, y=y, order=order)
annotator.configure(test='t-test_ind', text_format='star', loc='inside')
annotator.apply_and_annotate()
plt.ylim(-0.28,0.34)
ax.set_xlabel("Cluster", fontsize=18)
ax.set_ylabel("ssGSEA score", fontsize=18)
plt.title("HALLMARK_TNFA_SIGNALING_VIA_NFKB",fontsize=24)
plt.xticks(fontsize=18,rotation=90)
plt.yticks(fontsize=18)
plt.subplots_adjust(left=0.20)
plt.subplots_adjust(bottom=0.40)
plt.savefig("glioma/HALLMARK_TNFA_SIGNALING_VIA_NFKB_varplot.pdf",format="pdf")









tfactnorm[cellleiden==0,0].mean()
tfactnorm[cellleiden==1,0].mean()
tfactnorm[cellleiden==2,0].mean()
tfactnorm[cellleiden==3,0].mean()

tfactnorm[cellleiden==0,1].mean()
tfactnorm[cellleiden==1,1].mean()
tfactnorm[cellleiden==2,1].mean()
tfactnorm[cellleiden==3,1].mean()

tfactnorm[cellleiden==0,2].mean()
tfactnorm[cellleiden==1,2].mean()
tfactnorm[cellleiden==2,2].mean()
tfactnorm[cellleiden==3,2].mean()

tfactnorm[cellleiden==0,3].mean()
tfactnorm[cellleiden==1,3].mean()
tfactnorm[cellleiden==2,3].mean()
tfactnorm[cellleiden==3,3].mean()

emtsignal[cellleiden==0].mean()
emtsignal[cellleiden==1].mean()
emtsignal[cellleiden==2].mean()
emtsignal[cellleiden==3].mean()

enhatmpuse=enhamtx_n.transpose(0,2,1)[useidx==1]
enhatmpuse1=enhatmpuse[:,cellleiden==1]
enhatmpuse2=enhatmpuse[:,cellleiden==2]
enhatmpuse3=enhatmpuse[:,cellleiden==3]

enhamean1=enhatmpuse1.mean(axis=1)
enhamean2=enhatmpuse2.mean(axis=1)
enhamean3=enhatmpuse3.mean(axis=1)

np.abs(enhamean1).mean()
np.abs(enhamean2).mean()
np.abs(enhamean3).mean()

np.abs(enhamean1).mean()/np.abs(enhamean2).mean()
np.abs(enhamean3).mean()/np.abs(enhamean2).mean()


np.abs(enhamean1[enhamean1>0]).mean()/np.abs(enhamean2[enhamean2>0]).mean()
np.abs(enhamean3[enhamean3>0]).mean()/np.abs(enhamean2[enhamean2>0]).mean()

np.abs(enhamean1[enhamean1<0]).mean()/np.abs(enhamean2[enhamean2<0]).mean()
np.abs(enhamean3[enhamean3<0]).mean()/np.abs(enhamean2[enhamean2<0]).mean()

(enhamean1).mean()
(enhamean2).mean()
(enhamean3).mean()

(enhamean1).mean()/(enhamean2).mean()
(enhamean3).mean()/(enhamean2).mean()

fig, ax = plt.subplots()
ax.hist(enhatmpuse1.mean(axis=1),bins=500,range=(-0.00002,0.00002))
plt.savefig(samplename+"/histclus1.png")
plt.show()

fig, ax = plt.subplots()
ax.hist(enhatmpuse2.mean(axis=1),bins=500,range=(-0.00002,0.00002))
plt.savefig(samplename+"/histclus2.png")
plt.show()

fig, ax = plt.subplots()
ax.hist(enhatmpuse3.mean(axis=1),bins=500,range=(-0.00002,0.00002))
plt.savefig(samplename+"/histclus3.png")
plt.show()

stats.ttest_ind(emtrna[cellleiden==1],emtrna[cellleiden==0])
stats.ttest_ind(emtrna[cellleiden==2],emtrna[cellleiden==1])
stats.ttest_ind(emtrna[cellleiden==2],emtrna[cellleiden==3])

stats.ttest_ind(tfactnorm[cellleiden==1,3], tfactnorm[cellleiden==0,3])

for i in range(15):
   print(stats.ttest_ind(tfactnorm[cellleiden==1,i], tfactnorm[cellleiden==0,i]))

cnomin2,cnomax2=np.percentile(cno, [25,75])

cnotag2=(cno>cnomin2)&(cno<cnomax2)

from scipy.stats import spearmanr
np.corrcoef([tfactnorm[(cellleiden>=0)&(cellleiden<=2),3],emtrna[(cellleiden>=0)&(cellleiden<=2)]])

np.corrcoef([tfactnorm[cellclus==0,3],emtrna[cellclus==0]])

np.corrcoef([tfactnorm[(cellclus==0)&(cnotag==1),3],gsea["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"].to_numpy()])
np.corrcoef([tfactnorm[(cellclus==0)&(cnotag==1),15],gsea["HALLMARK_TNFA_SIGNALING_VIA_NFKB"].to_numpy()])


spearmanr(tfactnorm[(cellclus==0)&(cnotag==1),3],gsea["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"].to_numpy())
spearmanr(tfactnorm[(cellclus==0)&(cnotag==1),4],gsea["HALLMARK_TNFA_SIGNALING_VIA_NFKB"].to_numpy())

for i in range(tfactnorm.shape[1]):
   print(np.corrcoef([tfactnorm[(cellclus==0)&(cnotag==1),i],gsea["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"].to_numpy()])[0,1])


for i in range(tfactnorm.shape[1]):
   print(np.corrcoef([tfactnorm[(cellclus==0)&(cnotag==1),i],gsea["HALLMARK_TNFA_SIGNALING_VIA_NFKB"].to_numpy()])[0,1])

np.corrcoef([tfactnorm[cellclus==0,10],emtrna[cellclus==0]])

np.corrcoef([tfactnorm[(cellclus==0)&(cnotag==1),3],emtrna[(cellclus==0)&(cnotag==1)]])
np.corrcoef([tfactnorm[(cellclus==0)&(cnotag==1),15],nfkbrna[(cellclus==0)&(cnotag==1)]])

np.corrcoef([tfactnorm[(cellclus==0)&(cnotag2==1),3],emtrna[(cellclus==0)&(cnotag2==1)]])

correlation, pvalue = spearmanr(tfactnorm[(cellclus==0)&(cnotag==1),3],emtrna[(cellclus==0)&(cnotag==1)])

spearmanr(tfactnorm[(cellclus==0)&(cnotag==1),10],emtrna[(cellclus==0)&(cnotag==1)])

spearmanr(tfactnorm[(cellclus==0)&(cnotag==1),0],emtrna[(cellclus==0)&(cnotag==1)])
spearmanr(tfactnorm[(cellclus==0)&(cnotag==1),1],emtrna[(cellclus==0)&(cnotag==1)])
spearmanr(tfactnorm[(cellclus==0)&(cnotag==1),2],emtrna[(cellclus==0)&(cnotag==1)])

spearmanr(tfactnorm[(cellclus==0)&(cnotag==1),11],emtrna[(cellclus==0)&(cnotag==1)])

spearmanr(tfactnorm[(cellclus==0)&(cnotag==1),2]+tfactnorm[(cellclus==0)&(cnotag==1),3],emtrna[(cellclus==0)&(cnotag==1)])

spearmanr(tfactnorm[(cellclus==0)&(cnotag==1),3],emtrna[(cellclus==0)&(cnotag==1)])

fig, ax = plt.subplots()
ax.scatter(tfactnorm[(cellclus==0)&(cnotag==1),3],emtrna[(cellclus==0)&(cnotag==1)],s=1,c=(cellleiden[(cellclus==0)&(cnotag==1)]).astype(float))
#plt.xlabel('UMAP1')
#plt.ylabel("UMAP2")
plt.legend(fontsize=20)
fname=samplename+"/scatter_zeb1_emtrna.png"
plt.savefig(fname)
plt.show()


fig, ax = plt.subplots()
ax.scatter(tfactnorm[(cellclus==0)&(cnotag==1),3],gsea["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"].to_numpy(),s=1,c=(cellleiden[(cellclus==0)&(cnotag==1)]).astype(float))
#plt.xlabel('UMAP1')
#plt.ylabel("UMAP2")
plt.xlim(-0.05,-0.03)
plt.legend(fontsize=20)
fname=samplename+"/scatter_zeb1_emtgsea.png"
plt.savefig(fname)
plt.show()


np.corrcoef([tfactnorm[(cellclus==0)&(cnotag==1),3],gsea["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"].to_numpy()])



for i in range(15):
   tmp1=((tfidx[:,:,i]==1)&(nfkbmtx==1)&(promenhatag==4)).sum()/((nfkbmtx==1)&(promenhatag==4)).sum()
   tmp2=((tfidx[:,:,i]==1)&(promenhatag==4)).sum()/((promenhatag==4)).sum()
   print(tmp1/tmp2)


for i in range(15):
   tmp1=((tfidx[:,:,i]==1)&(emtmtx==1)&(promenhatag==4)).sum()/((emtmtx==1)&(promenhatag==4)).sum()
   tmp2=((tfidx[:,:,i]==1)&(promenhatag==4)).sum()/((promenhatag==4)).sum()
   print(tmp1)

for i in range(15):
   tmp1=((tfidxr[:,:,i]==1)&(emtmtx==1)&(promenhatag==4)).sum()/((tfidxr[:,:,i]==1)&(emtmtx==1)&(useidx==1)).sum()
   print(tmp1)

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



emtpeak=peaklistout.iloc[np.unique(peakidmtx[(emtmtx==1)&((promenhatag==2)|(promenhatag==4))]).astype(int),:]
nfkbpeak=peaklistout.iloc[np.unique(peakidmtx[(nfkbmtx==1)&((promenhatag==2)|(promenhatag==4))]).astype(int),:]

emtpeak.to_csv("glioma/emtpeak.bed",sep="\t",header=False, index=False)
nfkbpeak.to_csv("glioma/nfkbpeak.bed",sep="\t",header=False, index=False)

sox2cre=peaklistout.iloc[np.unique(peakidmtx[(tfidxr[:,:,1]==1)&(enhaclusgenemtxall[0]==1)]).astype(int),:]

sox2noncre=peaklistout.iloc[np.unique(peakidmtx[(tfidxr[:,:,1]==1)&((promenhatag!=2)&(promenhatag!=4))]).astype(int),:]

sox2cre.to_csv("glioma/sox2cre.bed",sep="\t",header=False, index=False)
sox2noncre.to_csv("glioma/sox2noncre.bed",sep="\t",header=False, index=False)

#####
fastaFromBed -fi Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa -bed glioma/clus0peak.bed -fo glioma/clus0peak.fasta
fastaFromBed -fi Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa -bed glioma/clus1peak.bed -fo glioma/clus1peak.fasta
fastaFromBed -fi Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa -bed glioma/clus2peak.bed -fo glioma/clus2peak.fasta
fastaFromBed -fi Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa -bed glioma/clus3peak.bed -fo glioma/clus3peak.fasta
####
fastaFromBed -fi Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa -bed glioma/upout.bed -fo glioma/upout.fasta
fastaFromBed -fi Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa -bed glioma/downout.bed -fo glioma/downout.fasta
####

fastaFromBed -fi Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa -bed glioma/emtpeak.bed -fo glioma/emtpeak.fasta
fastaFromBed -fi Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa -bed glioma/nfkbpeak.bed -fo glioma/nfkbpeak.fasta
####

bedtools intersect -wa -u -a glioma/clus0peak.bed -b glioma/GSE58345.bed > glioma/tmp0.bed
bedtools intersect -wa -u -a glioma/clus1peak.bed -b glioma/GSE58345.bed > glioma/tmp1.bed
bedtools intersect -wa -u -a glioma/clus2peak.bed -b glioma/GSE58345.bed > glioma/tmp2.bed
bedtools intersect -wa -u -a glioma/clus3peak.bed -b glioma/GSE58345.bed > glioma/tmp3.bed

wc -l glioma/tmp0.bed
wc -l glioma/tmp1.bed
wc -l glioma/tmp2.bed
wc -l glioma/tmp3.bed

wc -l glioma/clus0peak.bed
wc -l glioma/clus1peak.bed
wc -l glioma/clus2peak.bed
wc -l glioma/clus3peak.bed


bedtools intersect -wa -u -F 1 -a glioma/sox2cre.bed -b glioma/GSM6008250_SMNB19_SOX2_3_broadPeak.bed > glioma/tmp0.bed
bedtools intersect -wa -u -F 1  -a glioma/sox2noncre.bed -b glioma/GSM6008250_SMNB19_SOX2_3_broadPeak.bed > glioma/tmp1.bed

bedtools intersect -wa -u -a glioma/sox2cre.bed -b glioma/GSM6008250_SMNB19_SOX2_3_broadPeak_e3.bed > glioma/tmp0.bed
bedtools intersect -wa -u  -a glioma/sox2noncre.bed -b glioma/GSM6008250_SMNB19_SOX2_3_broadPeak_e3.bed > glioma/tmp1.bed


bedtools intersect -wa -u -f 1 -b glioma/sox2cre.bed -a glioma/GSM6008250_SMNB19_SOX2_3_broadPeak.bed > glioma/tmp0.bed
bedtools intersect -wa -u -f 1  -b glioma/sox2noncre.bed -a glioma/GSM6008250_SMNB19_SOX2_3_broadPeak.bed > glioma/tmp1.bed

wc -l glioma/tmp0.bed
wc -l glioma/tmp1.bed

wc -l glioma/sox2cre.bed
wc -l glioma/sox2noncre.bed

bedtools intersect -a glioma/clus0peak.bed -b glioma/clus1peak.bed > glioma/tmp.bed


######
peakssox=pd.read_csv("glioma/GSM6008248_SMNB19_SOX2_1_broadPeak.bed",sep="\t",header=None)
peakssox[peakssox[8]>3].to_csv("glioma/GSM6008248_SMNB19_SOX2_1_broadPeak_e3.bed",sep="\t",header=False, index=False)

peakssox=pd.read_csv("glioma/GSM6008249_SMNB19_SOX2_2_broadPeak.bed",sep="\t",header=None)
peakssox[peakssox[8]>3].to_csv("glioma/GSM6008249_SMNB19_SOX2_2_broadPeak_e3.bed",sep="\t",header=False, index=False)

peakssox=pd.read_csv("glioma/GSM6008250_SMNB19_SOX2_3_broadPeak.bed",sep="\t",header=None)
peakssox[peakssox[8]>3].to_csv("glioma/GSM6008250_SMNB19_SOX2_3_broadPeak_e3.bed",sep="\t",header=False, index=False)

bedtools intersect -wa -u -a glioma/sox2cre.bed -b glioma/GSM6008248_SMNB19_SOX2_1_broadPeak_e3.bed > glioma/tmp10.bed
bedtools intersect -wa -u  -a glioma/sox2noncre.bed -b glioma/GSM6008248_SMNB19_SOX2_1_broadPeak_e3.bed > glioma/tmp11.bed

bedtools intersect -wa -u -a glioma/sox2cre.bed -b glioma/GSM6008249_SMNB19_SOX2_2_broadPeak_e3.bed > glioma/tmp20.bed
bedtools intersect -wa -u  -a glioma/sox2noncre.bed -b glioma/GSM6008249_SMNB19_SOX2_2_broadPeak_e3.bed > glioma/tmp21.bed

bedtools intersect -wa -u -a glioma/sox2cre.bed -b glioma/GSM6008250_SMNB19_SOX2_3_broadPeak_e3.bed > glioma/tmp30.bed
bedtools intersect -wa -u  -a glioma/sox2noncre.bed -b glioma/GSM6008250_SMNB19_SOX2_3_broadPeak_e3.bed > glioma/tmp31.bed

wc -l glioma/tmp10.bed

wc -l glioma/tmp11.bed

wc -l glioma/tmp20.bed

wc -l glioma/tmp21.bed

wc -l glioma/tmp30.bed

wc -l glioma/tmp31.bed

wc -l glioma/sox2cre.bed

wc -l glioma/sox2noncre.bed

#####


bedtools intersect -wa -u -f 1 -b glioma/sox2cre.bed -a glioma/GSM6008248_SMNB19_SOX2_1_broadPeak.bed > glioma/tmp10b.bed
bedtools intersect -wa -u -f 1  -b glioma/sox2noncre.bed -a glioma/GSM6008248_SMNB19_SOX2_1_broadPeak.bed > glioma/tmp11b.bed


bedtools intersect -wa -u -f 1 -b glioma/sox2cre.bed -a glioma/GSM6008249_SMNB19_SOX2_2_broadPeak.bed > glioma/tmp20b.bed
bedtools intersect -wa -u -f 1  -b glioma/sox2noncre.bed -a glioma/GSM6008249_SMNB19_SOX2_2_broadPeak.bed > glioma/tmp21b.bed


bedtools intersect -wa -u -f 1 -b glioma/sox2cre.bed -a glioma/GSM6008250_SMNB19_SOX2_3_broadPeak.bed > glioma/tmp30b.bed
bedtools intersect -wa -u -f 1  -b glioma/sox2noncre.bed -a glioma/GSM6008250_SMNB19_SOX2_3_broadPeak.bed > glioma/tmp31b.bed

tmp10b=pd.read_csv("glioma/tmp10b.bed",sep="\t",header=None)
tmp11b=pd.read_csv("glioma/tmp11b.bed",sep="\t",header=None)
tmp20b=pd.read_csv("glioma/tmp20b.bed",sep="\t",header=None)
tmp21b=pd.read_csv("glioma/tmp21b.bed",sep="\t",header=None)
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



####
ATAC_tmp=ATAC_use[np.unique(peakidmtx[((nfkbmtx==1)&(promenhatag==4))]).astype(int),:]

tmp0=ATAC_tmp[:,cellleiden==0].mean(axis=1)
tmp1=ATAC_tmp[:,cellleiden==1].mean(axis=1)

enhatmp=enhamtx.transpose(0,2,1)[((emtmtx==1)&(promenhatag==4))]

tmp2=enhatmp[:,cellleiden==0].mean(axis=1)
tmp3=enhatmp[:,cellleiden==1].mean(axis=1)


enhathr=np.percentile(enhatmp,75,axis=1)
enhatmp2=(enhatmp.transpose(1,0)>enhathr).transpose(1,0)
enhatmp2[:,cellleiden==0].mean()
enhatmp2[:,cellleiden==1].mean()


enhathr=np.percentile(enhamtx.transpose(0,2,1),75,axis=2) #B,L
enhathrmtx=(enhamtx.transpose(1,0,2)>enhathr).transpose(1,0,2) #B,C,L

enhausemtx=((enhathrmtx.transpose(1,0,2)*((promenhatag==4)|(promenhatag==2))).sum(axis=2))/((promenhatag==4)|(promenhatag==2)).sum(axis=1)

enhausetagmtx=enhausemtx.transpose(1,0)>0.7

np.save('glioma/enhathrmtx.npy',enhathrmtx)
np.save('glioma/enhausetagmtx.npy',enhausetagmtx)

ATAC_corrt[enhaclusgenemtxall.sum(axis=0)>0]


useidx=np.unique(peakidmtx[enhaclusgenemtxall.sum(axis=0)>0]).astype(int)
pd.DataFrame(useidx).to_csv("glioma/tumorpeakidx.txt",sep="\t",header=False, index=False)

##################################


enhaname=DAEmtx.iloc[:,4]
usetag=enhaname[((pval.iloc[:,4]<0.001)&(fc.iloc[:,4]>0))].to_numpy()
usetag=usetag.astype(int)

strogene=genelist[(promenhatag==4).sum(axis=1)>0]

clus0gene=genelist[(enhaclusgenemtxall[0]==1).sum(axis=1)>0]
clus1gene=genelist[(enhaclusgenemtxall[1]==1).sum(axis=1)>0]


pd.DataFrame(genelist).to_csv('glioma/genelistfull.txt',sep="\t",header=False, index=False)
pd.DataFrame(clus0gene).to_csv('glioma/clus0genenew2.txt',sep="\t",header=False, index=False)
pd.DataFrame(clus1gene).to_csv('glioma/clus1genenew2.txt',sep="\t",header=False, index=False)

pd.DataFrame(strogene).to_csv('glioma/strogenelist.txt',sep="\t",header=False, index=False)

ascl1gene=genelist[(tfidx[:,:,13]==1).sum(axis=1)>0]
pd.DataFrame(ascl1gene).to_csv('glioma/ascl1gene.txt',sep="\t",header=False, index=False)


nfibgene=genelist[(tfidx[:,:,14]==1).sum(axis=1)>0]
pd.DataFrame(nfibgene).to_csv('glioma/nfibgene.txt',sep="\t",header=False, index=False)


ap1gene=genelist[(tfidx[:,:,4]==1).sum(axis=1)>0]
pd.DataFrame(ap1gene).to_csv('glioma/ap1gene.txt',sep="\t",header=False, index=False)

smad3ap1gene=genelist[(tfidx[:,:,5]==1).sum(axis=1)>0]
pd.DataFrame(smad3ap1gene).to_csv('glioma/smad3ap1gene.txt',sep="\t",header=False, index=False)

tfidx
tfidx

enhaclusgenemtxall.sum(axis=0)>0

enhatumor=enhamtx.transpose(0,2,1)
enhatumor=enhatumor[(promenhatag==4)&(enhaclusgenemtxall.sum(axis=0)>0)]
enhatumorpow=pow(2, enhatumor.transpose(1,0))

pca = PCA()
embedding_comb = umap.UMAP().fit_transform(enhatumor.transpose(1,0))
embedding_comb = pca.fit_transform(enhatumor[:,cellclus==0].transpose(1,0))
kmeans_embed = KMeans(n_clusters=4,max_iter=30, init="random")
cluster_kmeans = kmeans_embed.fit_predict(enhatumor.transpose(1,0))

fig, ax = plt.subplots()
ax.scatter(embedding_comb[:,0],embedding_comb[:,1],s=1,c=(cellclus).astype(float))
plt.xlabel('UMAP1')
plt.ylabel("UMAP2")
plt.legend(fontsize=20)
fname=samplename+"/umap_kmeans_promenha_new_tumor.png"
plt.savefig(fname)
plt.show()



embedding_comb = umap.UMAP().fit_transform(enhatumor[:,cellclus==0].transpose(1,0))
kmeans_embed = KMeans(n_clusters=4,max_iter=30, init="random")
cluster_kmeans = kmeans_embed.fit_predict(enhatumor[:,cellclus==0].transpose(1,0))
fig, ax = plt.subplots()
ax.scatter(embedding_comb[:,0],embedding_comb[:,1],s=1,c=(cluster_kmeans).astype(float))
plt.xlabel('UMAP1')
plt.ylabel("UMAP2")
plt.legend(fontsize=20)
fname=samplename+"/umap_kmeans_promenha_new_tumor2.png"
plt.savefig(fname)
plt.show()
####################

ATACuse_select=ATACr[useidx]
ATACuse_select=ATACr

normatacr=ATACr.mean(axis=0)

ATACuse_select=ATACr/normatacr

emtgene=pd.read_csv(samplename+"/emtgene.txt",sep="\t",header=None)
emttag=np.isin(genelist,emtgene)
emtmtx=np.ones(promenhatag.shape)
emtmtx=(emtmtx.transpose(1,0)*(emttag.astype(int))).transpose(1,0)
ATACremt=ATACr/normatacr
ATACremt=ATACremt[np.unique(peakidmtx[(promenhatag==4)&(emtmtx==1)]).astype(int),:]
emtsignalatacr=ATACremt.sum(axis=0)
emtminatacr,emtmaxatacr=np.percentile(emtsignalatacr[cellclus==0], [10,90])

krasgene=pd.read_csv(samplename+"/krasgene.txt",sep="\t",header=None)
krastag=np.isin(genelist,krasgene)
krasmtx=np.ones(promenhatag.shape)
krasmtx=(krasmtx.transpose(1,0)*(krastag.astype(int))).transpose(1,0)
ATACrkras=ATACr/normatacr
ATACrkras=ATACrkras[np.unique(peakidmtx[(promenhatag==4)&(krasmtx==1)]).astype(int),:]
krassignalatacr=ATACrkras.sum(axis=0)
krasminatacr,krasmaxatacr=np.percentile(krassignalatacr[cellclus==0], [10,90])

tgfbgene=pd.read_csv(samplename+"/tgfbgene.txt",sep="\t",header=None)
tgfbtag=np.isin(genelist,tgfbgene)
tgfbmtx=np.ones(promenhatag.shape)
tgfbmtx=(tgfbmtx.transpose(1,0)*(tgfbtag.astype(int))).transpose(1,0)
ATACrtgfb=ATACr/normatacr
ATACrtgfb=ATACrtgfb[np.unique(peakidmtx[(promenhatag==4)&(tgfbmtx==1)]).astype(int),:]
tgfbsignalatacr=ATACrtgfb.sum(axis=0)
tgfbminatacr,tgfbmaxatacr=np.percentile(tgfbsignalatacr[cellclus==0], [10,90])

hypogene=pd.read_csv(samplename+"/hypogene.txt",sep="\t",header=None)
hypotag=np.isin(genelist,hypogene)
hypomtx=np.ones(promenhatag.shape)
hypomtx=(hypomtx.transpose(1,0)*(hypotag.astype(int))).transpose(1,0)
ATACrhypo=ATACr/normatacr
ATACrhypo=ATACrhypo[np.unique(peakidmtx[(promenhatag==4)&(hypomtx==1)]).astype(int),:]
hyposignalatacr=ATACrhypo.sum(axis=0)
hypominatacr,hypomaxatacr=np.percentile(hyposignalatacr[cellclus==0], [10,90])

hghggene=pd.read_csv(samplename+"/hghggene.txt",sep="\t",header=None)
hghgtag=np.isin(genelist,hghggene)
hghgmtx=np.ones(promenhatag.shape)
hghgmtx=(hghgmtx.transpose(1,0)*(hghgtag.astype(int))).transpose(1,0)
ATACrhghg=ATACr/normatacr
ATACrhghg=ATACrhghg[np.unique(peakidmtx[(promenhatag==4)&(hghgmtx==1)]).astype(int),:]
hghgsignalatacr=ATACrhghg.sum(axis=0)
hghgminatacr,hghgmaxatacr=np.percentile(hghgsignalatacr[cellclus==0], [10,90])


pd.DataFrame(np.unique(peakidmtx[(promenhatag==4)&(emtmtx==1)]).astype(int)).to_csv('glioma/emtpeaks.txt',sep="\t",header=False, index=False)
pd.DataFrame(np.unique(peakidmtx[(promenhatag==4)&(krasmtx==1)]).astype(int)).to_csv('glioma/kraspeaks.txt',sep="\t",header=False, index=False)
pd.DataFrame(np.unique(peakidmtx[(promenhatag==4)&(tgfbmtx==1)]).astype(int)).to_csv('glioma/tgfbpeaks.txt',sep="\t",header=False, index=False)
pd.DataFrame(np.unique(peakidmtx[(promenhatag==4)&(hypomtx==1)]).astype(int)).to_csv('glioma/hypopeaks.txt',sep="\t",header=False, index=False)
pd.DataFrame(np.unique(peakidmtx[(promenhatag==4)&(hghgmtx==1)]).astype(int)).to_csv('glioma/hghgpeaks.txt',sep="\t",header=False, index=False)


enhar = sc.AnnData(ATACuse_select[:,cellclus==0].transpose(1,0))
#enha = sc.AnnData(df2, df2.index.to_frame(), df2.columns.to_frame())
sc.tl.pca(enhar, svd_solver="arpack")
sc.pp.neighbors(enhar, n_neighbors=10, n_pcs=40)
sc.tl.umap(enhar)
sc.tl.tsne(enhar)
sc.pl.umap(enhar, frameon=False, save="glioma_enhancer_rawatac.png")


sc.tl.leiden(
    enhar,
    resolution=0.1
)

enhar.obs["emt"] = emtsignalatacr[cellclus==0]
enhar.obs["kras"] = krassignalatacr[cellclus==0]
enhar.obs["tgfb"] = tgfbsignalatacr[cellclus==0]
enhar.obs["hypo"] = hyposignalatacr[cellclus==0]
enhar.obs["hghg"] = hghgsignalatacr[cellclus==0]

sc.pl.umap(enhar, frameon=False,legend_loc="on data", color="leiden",save="glioma_enhancer2_rawatac.png")
sc.pl.tsne(enhar, frameon=False,legend_loc="on data", color="leiden",save="glioma_enhancer2t_rawatac.png")
enhar.obs["rnacluster"] = pd.Categorical(list(cellleiden[cellclus==0]))

sc.pl.umap(enhar, frameon=False, color="rnacluster",save="glioma_enhancer3_rawatac.png")

sc.pl.umap(enhar, frameon=False, color="emt",save="glioma_enhancer3emt_rawatac.png",vmin=emtminatacr, vmax=emtmaxatacr)
sc.pl.umap(enhar, frameon=False, color="kras",save="glioma_enhancer3kras_rawatac.png",vmin=krasminatacr, vmax=krasmaxatacr)
sc.pl.umap(enhar, frameon=False, color="tgfb",save="glioma_enhancer3tgfb_rawatac.png",vmin=tgfbminatacr, vmax=tgfbmaxatacr)
sc.pl.umap(enhar, frameon=False, color="hypo",save="glioma_enhancer3hypo_rawatac.png",vmin=hypominatacr, vmax=hypomaxatacr)
sc.pl.umap(enhar, frameon=False, color="hghg",save="glioma_enhancer3hghg_rawatac.png",vmin=hghgminatacr, vmax=hghgmaxatacr)





ATACuse_select=ATAC_use[useidx]

enhap = sc.AnnData(ATACuse_select[:,cellclus==0].transpose(1,0))
#enha = sc.AnnData(df2, df2.index.to_frame(), df2.columns.to_frame())
sc.tl.pca(enhap, svd_solver="arpack")
sc.pp.neighbors(enhap, n_neighbors=10, n_pcs=40)
sc.tl.umap(enhap)
sc.tl.tsne(enhap)
sc.pl.umap(enhap, frameon=False, save="glioma_enhancer_predatac.png")


sc.tl.leiden(
    enhap,
    resolution=0.1
)

sc.pl.umap(enhap, frameon=False,legend_loc="on data", color="leiden",save="glioma_enhancer2_predatac.png")
sc.pl.tsne(enhap, frameon=False,legend_loc="on data", color="leiden",save="glioma_enhancer2t_predatac.png")
enhap.obs["rnacluster"] = pd.Categorical(list(cellleiden[cellclus==0]))

sc.pl.umap(enhap, frameon=False, color="rnacluster",save="glioma_enhancer3_predatac.png")

ATAC_use

######################


pairlist_prom=pd.read_csv('neuro/pair_promoter.csv',header=None)
pairlist_prom=pairlist_prom[[1,2,3]].to_numpy()

pairlist=pd.read_csv('neuro/pair_300000.csv',header=None)
pairlist=pairlist[[1,2,3]].to_numpy()


peaklist=pd.read_csv('neuro/peaks_extend.bed',sep="\t",header=None)
peaklist=peaklist[[1,2]].to_numpy()

embedlist=pd.read_csv("neuro/embed.txt",sep=" ",header=None)
embedlist=embedlist.to_numpy()



posmtx=np.zeros((grad.shape[0],grad.shape[2]))
for j in range(grad.shape[0]):
    peakstart=int(pairlist[j,0].item())
    peakend=int(pairlist[j,1].item())
    prompos=int(pairlist_prom[j,0].item())
    enha=list(range(peakstart,(peakend+1)))
    enha.remove(prompos)
    genepos=int(pairlist[j,2].item())
    peaknum_gene=peakend-peakstart+1
    posvec=((((peaklist[enha,0]+peaklist[enha,1])/2)-genepos))
    posmtx[j,1:peaknum_gene]=posvec

peakidmtx=np.zeros((grad.shape[0],grad.shape[2]))
for j in range(grad.shape[0]):
    peakstart=int(pairlist[j,0].item())
    peakend=int(pairlist[j,1].item())
    prompos=int(pairlist_prom[j,0].item())
    enha=list(range(peakstart,(peakend+1)))
    enha.remove(prompos)
    peaknum_gene=peakend-peakstart+1
    peakidmtx[j,0]=prompos
    peakidmtx[j,1:peaknum_gene]=enha



comb_enha_embed=np.zeros((promenhatag.shape[0],promenhatag.shape[1],32,2))
for i in range(comb_enha_embed.shape[0]):
    for j in range(comb_enha_embed.shape[1]):
      idx=peakidmtx[i,j].astype(int)
      pidx=peakidmtx[i,0].astype(int)
      embed_tmp=embedlist[idx,:]
      pembe_tmp=embedlist[pidx,:]
      comb_enha_embed[i,j,:,0]=embed_tmp
      comb_enha_embed[i,j,:,1]=pembe_tmp


comb_enha_embed_s=comb_enha_embed[(promenhatag==2)|(promenhatag==4)]




comb_enha_embed_re=comb_enha_embed_s.reshape(comb_enha_embed_s.shape[0],comb_enha_embed_s.shape[1]*comb_enha_embed_s.shape[2])
promemblist=pd.DataFrame(comb_enha_embed[:,:,1]).drop_duplicates().to_numpy()
comb_enha_embed_comb=np.concatenate([comb_enha_embed[:,:,0],promemblist])

from sklearn.decomposition import PCA 

pca = PCA()
pca.fit(comb_enha_embed_re)
embedding_comb_pca = pca.transform(comb_enha_embed_re)


fig, ax = plt.subplots()
ax.scatter(embedding_comb_pca[:,0],embedding_comb_pca[:,1],s=3,c=(tumorspe>0).astype(float),alpha=0.3)
plt.xlabel('UMAP1')
plt.ylabel("UMAP2")
plt.legend(fontsize=20)
fname=samplename+"/tumor_embed_pca.png"
plt.savefig(fname)
plt.show()




pca = PCA()
pca.fit(comb_enha_embed_s[:,:,0])
embedding_comb_pca = pca.transform(comb_enha_embed_s[:,:,0])


fig, ax = plt.subplots()
ax.scatter(embedding_comb_pca[:,0],embedding_comb_pca[:,1],s=3,c=(tumorspe>0).astype(float),alpha=0.3)
plt.xlabel('UMAP1')
plt.ylabel("UMAP2")
plt.legend(fontsize=20)
fname=samplename+"/tumor_embed_pca2.png"
plt.savefig(fname)
plt.show()


################

import seaborn as sns

enhatumor=enhamtxuse
normfac=np.abs(enhatumor).sum(axis=(0))
enhatumornorm=((enhamtxuse.transpose(1,0)-enhamtxuse.mean(axis=1))/enhamtxuse.std(axis=1)).transpose(1,0)


#sns.set(font_scale=1.2)
fig, ax= plt.subplots()
g2=sns.clustermap(pd.DataFrame(enhatumornorm))
#g2.ax_heatmap.set_xticklabels(g2.ax_heatmap.get_xmajorticklabels(), fontsize = 12,rotation=0)
#g2.ax_cbar.set_position([0.05,0.8,0.05,0.15])
#g2.subplots_adjust(right=0.8)
#g2.ax_cbar.set_position((0.85, .2, .03, .4))
fname=samplename+"/tumor_clustermap.png"
plt.savefig(fname)


kmeans_embed = KMeans(n_clusters=12,max_iter=30, init="random")
cluster_kmeans = kmeans_embed.fit_predict(enhatumornorm)


clusmtx=np.zeros((12,7))

cluscellmtx=np.zeros((12,enhatumornorm.shape[1]))

for i in range(12):
   tmp=enhatumornorm[cluster_kmeans==i].mean(axis=0)
   for j in range(7):
      clusmtx[i,j]=tmp[cellclus==j].mean()

for i in range(12):
   cluscellmtx[i,:]=enhatumornorm[cluster_kmeans==i].mean(axis=0)


embedding_comb = umap.UMAP().fit_transform(cluscellmtx.transpose(1,0))
fig, ax = plt.subplots()
ax.scatter(embedding_comb[:,0],embedding_comb[:,1],s=3,c=(normfac).astype(float),alpha=0.3)
plt.xlabel('UMAP1')
plt.ylabel("UMAP2")
plt.legend(fontsize=20)
fname=samplename+"/tumor_embed_pca5.png"
plt.savefig(fname)
plt.show()




kmeans_embed2 = KMeans(n_clusters=9,max_iter=30, init="random")
cluster_cell = kmeans_embed2.fit_predict(cluscellmtx.transpose(1,0))


normclusmtx=np.zeros((9,7))

cluscell=np.zeros((9,7))
for i in range(9):
   for j in range(7):
      cluscell[i,j]=((cluster_cell==i)&(cellclus==j)).sum()
      normclusmtx[i,j]=normfac[((cluster_cell==i)&(cellclus==j))].mean()

pd.DataFrame(cluscell)
pd.DataFrame(normclusmtx)

clusclusmtx=np.zeros((12,9))

for i in range(12):
   tmp=enhatumornorm[cluster_kmeans==i].mean(axis=0)
   for j in range(9):
      clusclusmtx[i,j]=tmp[cluster_cell==j].mean()

pd.DataFrame(clusclusmtx)


fig, ax= plt.subplots()
g2=sns.heatmap(pd.DataFrame(clusclusmtx))
#g2.ax_heatmap.set_xticklabels(g2.ax_heatmap.get_xmajorticklabels(), fontsize = 12,rotation=0)
#g2.ax_cbar.set_position([0.05,0.8,0.05,0.15])
#g2.subplots_adjust(right=0.8)
#g2.ax_cbar.set_position((0.85, .2, .03, .4))
fname=samplename+"/tumor_clusmtx.png"
plt.savefig(fname)


enhaclusgenemtx=np.zeros(promenhatag.shape)
enhaclusgenemtx=enhaclusgenemtx-1
enhaclusgenemtx[(promenhatag==2)|(promenhatag==4)]=cluster_kmeans



genelist=pd.read_csv('neuro/pair_300000.csv',header=None)
genelist=genelist[[0]].to_numpy()
genelist=genelist[:,0]


clus0gene=genelist[(enhaclusgenemtx==0).sum(axis=1)>0]
clus1gene=genelist[(enhaclusgenemtx==1).sum(axis=1)>0]
clus2gene=genelist[(enhaclusgenemtx==2).sum(axis=1)>0]
clus3gene=genelist[(enhaclusgenemtx==3).sum(axis=1)>0]
clus4gene=genelist[(enhaclusgenemtx==4).sum(axis=1)>0]
clus5gene=genelist[(enhaclusgenemtx==5).sum(axis=1)>0]
clus6gene=genelist[(enhaclusgenemtx==6).sum(axis=1)>0]
clus7gene=genelist[(enhaclusgenemtx==7).sum(axis=1)>0]

pd.DataFrame(genelist).to_csv('glioma/genelistfull.txt',sep="\t",header=False, index=False)
pd.DataFrame(clus0gene).to_csv('glioma/clus0gene.txt',sep="\t",header=False, index=False)
pd.DataFrame(clus1gene).to_csv('glioma/clus1gene.txt',sep="\t",header=False, index=False)
pd.DataFrame(clus2gene).to_csv('glioma/clus2gene.txt',sep="\t",header=False, index=False)
pd.DataFrame(clus3gene).to_csv('glioma/clus3gene.txt',sep="\t",header=False, index=False)
pd.DataFrame(clus4gene).to_csv('glioma/clus4gene.txt',sep="\t",header=False, index=False)
pd.DataFrame(clus5gene).to_csv('glioma/clus5gene.txt',sep="\t",header=False, index=False)
pd.DataFrame(clus6gene).to_csv('glioma/clus6gene.txt',sep="\t",header=False, index=False)
pd.DataFrame(clus7gene).to_csv('glioma/clus7gene.txt',sep="\t",header=False, index=False)


clus0geneonly=genelist[((enhaclusgenemtx==0).sum(axis=1)>0)&((enhaclusgenemtx>0).sum(axis=1)==0)]

pd.DataFrame(clus0geneonly).to_csv('glioma/clus0geneonly.txt',sep="\t",header=False, index=False)



###########




embedding_comb = umap.UMAP().fit_transform(enhatumornorm[cluster_kmeans==5].transpose(1,0))
#kmeans_embed = KMeans(n_clusters=6,max_iter=30, init="random")
#cluster_kmeans = kmeans_embed.fit_predict(cluscellmtx.transpose(1,0))

fig, ax = plt.subplots()
ax.scatter(embedding_comb[:,0], embedding_comb[:,1],s=1, c=(cluster_cell==5).astype(float))
plt.xlabel('UMAP1')
plt.ylabel("UMAP2")
plt.legend(fontsize=20)
fname=samplename+"/umap_kmeans_promenha_new4.png"
plt.savefig(fname)
plt.show()

######################
from scipy.stats import yeojohnson

grad_yj=np.zeros(enhamtx.shape)
for i in range(grad.shape[0]):
    print(i)
    tmpdf=enhamtx[i,:,:].copy()
    if (len(tmpdf[tmpdf!=0])>0):
        xt, _ = yeojohnson(tmpdf[tmpdf!=0])
        xt_norm=stats.zscore(xt)
        grad_yj[i,tmpdf!=0]=xt_norm

np.save(samplename+"/grad_yj.npy",grad_yj)






enhamtxuse=grad_yj.transpose(0,2,1)
normfac=grad_yj.sum(axis=(0,2))
enhamtxuse=enhamtxuse[(promenhatag==4)]


enhamtxnorm=grad_yj.transpose(0,2,1)
normfac=np.abs(enhamtxnorm).sum(axis=(0,1))
enhamtxnorm=enhamtxnorm/normfac
enhamtxuse=enhamtxnorm[(promenhatag==2)|(promenhatag==4)]



[:,cellclus==0]



pca = PCA()
pca.fit(enhamtxuse.transpose(1,0))
embedding_comb_pca = pca.transform(enhamtxuse.transpose(1,0))
embedding_comb = umap.UMAP().fit_transform(embedding_comb_pca[:,1:11])

embedding_comb = umap.UMAP().fit_transform(enhamtxuse.transpose(1,0))
kmeans_embed = KMeans(n_clusters=3,max_iter=30, init="random")
cluster_kmeans = kmeans_embed.fit_predict(enhamtxuse.transpose(1,0))


fig, ax = plt.subplots()
ax.scatter(embedding_comb[:,0], embedding_comb[:,1],s=1, c=(np.log(normfac)).astype(float))
plt.xlabel('UMAP1')
plt.ylabel("UMAP2")
plt.legend(fontsize=20)
fname=samplename+"/umap_kmeans_promenha_new.png"
plt.savefig(fname)
plt.show()

##################################




enhamtxuse=enhamtx.transpose(0,2,1)
normfac=enhamtx.sum(axis=(0,2))
enhamtxuse=enhamtxuse[(promenhatag==4)]


pca = PCA()
pca.fit(enhamtxuse.transpose(1,0))
embedding_comb_pca = pca.transform(enhamtxuse.transpose(1,0))


pcacorr=np.zeros((embedding_comb_pca.shape[1]))

for i in range(pcacorr.shape[0]):
   pcacorr[i]=np.corrcoef(normfac,embedding_comb_pca[:,i])[0,1]



fig, ax = plt.subplots()
ax.hist(pcacorr,bins=50)
plt.savefig(samplename+"/pcacorr.png")
plt.show()

embedding_comb_pcas=embedding_comb_pca[:,np.abs(pcacorr)<0.15]

embedding_comb = umap.UMAP().fit_transform(embedding_comb_pcas[:,:30])



fig, ax = plt.subplots()
ax.scatter(embedding_comb[:,0], embedding_comb[:,1],s=1, c=(np.log(normfac)).astype(float))
plt.xlabel('UMAP1')
plt.ylabel("UMAP2")
plt.legend(fontsize=20)
fname=samplename+"/umap_kmeans_promenha_new.png"
plt.savefig(fname)
plt.show()


###################




pairlist=pd.read_csv('glioma/pair_300000.csv',header=None,names=("gene","start","end","pos"))
pairlist_train=pd.read_csv('glioma/pair_300000_train.csv',header=None,names=("gene","start","end","pos"))
pairlist_test=pd.read_csv('glioma/pair_300000_test.csv',header=None,names=("gene","start","end","pos"))

traingenelist=pairlist["gene"].isin(pairlist_train["gene"]).values.tolist()
testgenelist=pairlist["gene"].isin(pairlist_test["gene"]).values.tolist()

traincellnum=np.load("glioma/traincellnum.npy")

tmp=enhamtx[testgenelist]
tmp=enhamtx
enhamtx_test=tmp[:,traincellnum:,:]
promenhatag_s=promenhatag[testgenelist]


enhamtxuse=enhamtx_test.transpose(0,2,1)
enhamtxuse=enhamtxuse[(promenhatag==4)]
enhamtxusepow=pow(2, enhamtxuse.transpose(1,0))

#[cellclus[traincellnum:]==0]

enha = sc.AnnData(enhamtxusepow)
#enha = sc.AnnData(df2, df2.index.to_frame(), df2.columns.to_frame())
sc.tl.pca(enha, svd_solver="arpack")
sc.pp.neighbors(enha, n_neighbors=10, n_pcs=40)
sc.tl.umap(enha)
sc.tl.tsne(enha)

sc.pl.umap(enha, frameon=False, save="glioma_enhancer.png")


sc.tl.leiden(
    enha,
    resolution=0.1
)

sc.pl.umap(enha, frameon=False,legend_loc="on data", color="leiden",save="glioma_enhancer2.png")
sc.pl.tsne(enha, frameon=False,legend_loc="on data", color="leiden",save="glioma_enhancer2t.png")


cellclus=pd.read_csv("glioma/cellcluster.txt",header=None)
cellclus=cellclus.to_numpy()
cellclus=cellclus[:,0]
enha.obs["rnacluster"] = pd.Categorical(list(cellclus[traincellnum:]))

cellleiden=np.zeros((cellclus.shape[0]))
cellleiden=cellleiden-1
cellleiden[cellclus==0]=enha.obs["leiden"].to_numpy()

np.save('glioma/cellleiden.npy',cellleiden)
pd.DataFrame(cellleiden).to_csv("glioma/cellleiden.txt",sep="\t",header=False, index=False)

sc.pl.umap(enha, frameon=False, color="rnacluster",save="glioma_enhancer3.png")

###################


cellclus_s=cellclus[traincellnum:]

enha = sc.AnnData(enhamtxusepow[cellclus_s==0])
#enha = sc.AnnData(df2, df2.index.to_frame(), df2.columns.to_frame())
sc.tl.pca(enha, svd_solver="arpack")
sc.pp.neighbors(enha, n_neighbors=10, n_pcs=40)
sc.tl.umap(enha)
sc.tl.tsne(enha)

sc.pl.umap(enha, frameon=False, save="glioma_enhancer.png")


sc.tl.leiden(
    enha,
    resolution=0.1
)

sc.pl.umap(enha, frameon=False,legend_loc="on data", color="leiden",save="glioma_enhancer2.png")
sc.pl.tsne(enha, frameon=False,legend_loc="on data", color="leiden",save="glioma_enhancer2t.png")


enha.obs["rnacluster"] = pd.Categorical(list(cellclus_s]))

cellleiden=np.zeros((cellclus_s.shape[0]))
cellleiden=cellleiden-1
cellleiden[cellclus_s==0]=enha.obs["leiden"].to_numpy()
np.save('glioma/cellleiden.txt',cellleiden)
pd.DataFrame(cellleiden).to_csv("glioma/cellleiden.txt",sep="\t",header=False, index=False)

sc.pl.umap(enha, frameon=False, color="rnacluster",save="glioma_enhancer3.png")


enhamtxusetumor=enhamtxuse[:,cellclus_s==0].mean(axis=1)
enhamtxuseother=enhamtxuse[:,cellclus_s!=0].mean(axis=1)

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
   pthr=np.percentile(pval.iloc[:,i], [30])
   usetag=enhaname[((pval.iloc[:,i]<pthr.item())&(fc.iloc[:,i]>0))].to_numpy()
   usetag=usetag.astype(int)
   tmp[usetag]=1
   enhaclusgenemtx=np.zeros(promenhatag.shape)
   enhaclusgenemtx=enhaclusgenemtx
   enhaclusgenemtx[(promenhatag==4)]=tmp
   enhaclusgenemtxall[i]=enhaclusgenemtx
   tgtgene=genelist[(enhaclusgenemtx==1).sum(axis=1)>0]
   pd.DataFrame(tgtgene).to_csv("glioma/clus"+str(i)+"genenewtest.txt",sep="\t",header=False, index=False)


np.save('glioma/enhaclusgenemtxall.npy',enhaclusgenemtxall)

tgtgene0=genelist[(enhaclusgenemtxall[0]==1).sum(axis=1)>0]
tgtgene1=genelist[(enhaclusgenemtxall[1]==1).sum(axis=1)>0]
tgtcommon=np.intersect1d(tgtgene0, tgtgene1)



ATAC_corrt[enhaclusgenemtxall.sum(axis=0)>0]


useidx=np.unique(peakidmtx[enhaclusgenemtxall.sum(axis=0)>0]).astype(int)
pd.DataFrame(useidx).to_csv("glioma/tumorpeakidx.txt",sep="\t",header=False, index=False)

enhaname=DAEmtx.iloc[:,4]
usetag=enhaname[((pval.iloc[:,4]<0.001)&(fc.iloc[:,4]>0))].to_numpy()
usetag=usetag.astype(int)

clus0gene=genelist[(enhaclusgenemtx==0).sum(axis=1)>0]
clus1gene=genelist[(enhaclusgenemtx==1).sum(axis=1)>0]
clus2gene=genelist[(enhaclusgenemtx==2).sum(axis=1)>0]
clus3gene=genelist[(enhaclusgenemtx==3).sum(axis=1)>0]
clus4gene=genelist[(enhaclusgenemtx==4).sum(axis=1)>0]
clus5gene=genelist[(enhaclusgenemtx==5).sum(axis=1)>0]
clus6gene=genelist[(enhaclusgenemtx==6).sum(axis=1)>0]
clus7gene=genelist[(enhaclusgenemtx==7).sum(axis=1)>0]

pd.DataFrame(genelist).to_csv('glioma/genelistfull.txt',sep="\t",header=False, index=False)
pd.DataFrame(clus2gene).to_csv('glioma/clus2genenew.txt',sep="\t",header=False, index=False)


###########################




ATACuse_select=ATACr[useidx]
ATACuse_select=ATACuse_select[:,traincellnum:]

enhar = sc.AnnData(ATACuse_select[:,cellclus_s==0].transpose(1,0))
#enha = sc.AnnData(df2, df2.index.to_frame(), df2.columns.to_frame())
sc.tl.pca(enhar, svd_solver="arpack")
sc.pp.neighbors(enhar, n_neighbors=10, n_pcs=40)
sc.tl.umap(enhar)
sc.tl.tsne(enhar)
sc.pl.umap(enhar, frameon=False, save="glioma_enhancer_rawatac.png")


sc.tl.leiden(
    enhar,
    resolution=0.1
)

sc.pl.umap(enhar, frameon=False,legend_loc="on data", color="leiden",save="glioma_enhancer2_rawatac.png")
sc.pl.tsne(enhar, frameon=False,legend_loc="on data", color="leiden",save="glioma_enhancer2t_rawatac.png")
enhar.obs["rnacluster"] = pd.Categorical(list(cellleiden[cellclus_s==0]))

sc.pl.umap(enhar, frameon=False, color="rnacluster",save="glioma_enhancer3_rawatac.png")






ATACuse_select=ATAC_use[useidx]
ATACuse_select=ATACuse_select[:,traincellnum:]

enhap = sc.AnnData(ATACuse_select[:,cellclus_s==0].transpose(1,0))
#enha = sc.AnnData(df2, df2.index.to_frame(), df2.columns.to_frame())
sc.tl.pca(enhap, svd_solver="arpack")
sc.pp.neighbors(enhap, n_neighbors=10, n_pcs=40)
sc.tl.umap(enhap)
sc.tl.tsne(enhap)
sc.pl.umap(enhap, frameon=False, save="glioma_enhancer_predatac.png")


sc.tl.leiden(
    enhap,
    resolution=0.1
)

sc.pl.umap(enhap, frameon=False,legend_loc="on data", color="leiden",save="glioma_enhancer2_predatac.png")
sc.pl.tsne(enhap, frameon=False,legend_loc="on data", color="leiden",save="glioma_enhancer2t_predatac.png")
enhap.obs["rnacluster"] = pd.Categorical(list(cellleiden[cellclus_s==0]))

sc.pl.umap(enhap, frameon=False, color="rnacluster",save="glioma_enhancer3_predatac.png")



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

sox2minallprom,sox2maxallprom=np.percentile(tfactprom[:,1], [10,90])
enhaall.obs["SOX2prom"] = tfactprom[:,1]
sc.pl.umap(enhaall, frameon=False, color="SOX2prom",save="glioma_enhancer3sox2_prom_all.pdf",vmin=sox2minallprom, vmax=sox2maxallprom)



sox2rna=pd.read_csv(samplename+"/sox2exp.txt",sep="\t",header=None)
sox2rna=sox2rna.to_numpy()
sox2rna=sox2rna[:,0]
#emtrna=np.log2(emtrna)
sox2rnaminall,sox2rnamaxall=np.percentile(sox2rna, [10,90])

enhaall.obs["SOX2_RNA"] = sox2rna
sc.pl.umap(enhaall, frameon=False, color="SOX2_RNA",save="glioma_enhancer3sox2rna_all.pdf",vmin=sox2rnaminall, vmax=sox2rnamaxall)



stats.ttest_ind(tfactnorm[cellclus==0,1], tfactnorm[cellclus!=0,1])

sc.tl.rank_genes_groups(enhaall, "leiden", method="wilcoxon")
sc.pl.rank_genes_groups(enhaall, n_genes=25, sharey=False, save="DAE_deeplift_pbmc_for_paper.png")


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