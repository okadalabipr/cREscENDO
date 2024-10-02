import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
import sys
import math
import subprocess
from scipy.stats import spearmanr
import scipy.stats as stats
import seaborn as sns
from statannotations.Annotator import Annotator

samplename="pbmc_after"

fname="pbmc_fold1/Deeplift_full_ver2_all.npy"
#fname="pbmc_fold1/allgrad_ssep_max.npy"
enhamtx=np.load(fname)


genelist=pd.read_csv(samplename+"/pair_promoter.csv",sep=",",header=None)


clus1gene=pd.read_csv(samplename+"/clus1gene.txt",sep="\t",header=None)
clus2gene=pd.read_csv(samplename+"/clus2gene.txt",sep="\t",header=None)
clus3gene=pd.read_csv(samplename+"/clus3gene.txt",sep="\t",header=None)
clus4gene=pd.read_csv(samplename+"/clus4gene.txt",sep="\t",header=None)
clus5gene=pd.read_csv(samplename+"/clus5gene.txt",sep="\t",header=None)
clus6gene=pd.read_csv(samplename+"/clus6gene.txt",sep="\t",header=None)


geneclus=np.zeros((genelist.shape[0],6))

geneclus[:,0]=genelist[0].isin(clus1gene[0])
geneclus[:,1]=genelist[0].isin(clus2gene[0])
geneclus[:,2]=genelist[0].isin(clus3gene[0])
geneclus[:,3]=genelist[0].isin(clus4gene[0])
geneclus[:,4]=genelist[0].isin(clus5gene[0])
geneclus[:,5]=genelist[0].isin(clus6gene[0])


grad=np.load("pbmc_fold1/allgrad_ssep_max.npy")


pairlist_prom=pd.read_csv('pbmc_fold1/pair_promoter.csv',header=None)
pairlist_prom=pairlist_prom[[1,2,3]].to_numpy()

pairlist=pd.read_csv('pbmc_fold1/pair_300000.csv',header=None)
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
      #ATAC_corr[j,i]=np.corrcoef([ATAC_tmp[:,i],enhamtx[j,:,i]])[0,1]
      ATAC_corr[j,i], _ = spearmanr(ATAC_tmp[:,i],enhamtx[j,:,i])

np.save("pbmc_fold1/ATAC_corr.npy",ATAC_corr)
np.save("pbmc_fold1/ATAC_corr_sp.npy",ATAC_corr)

ATAC_corr=np.load("pbmc_fold1/ATAC_corr.npy")

thv=np.percentile(grad[grad!=0], [40])
sthv=np.percentile(grad[grad!=0], [85])

promenhatag=np.zeros(grad.shape)
promenhatag[(grad>thv)&(ATAC_corr>0)]=2
promenhatag[(grad>sthv)&(ATAC_corr>0)]=4
promenhatag[:,0]=3

#################


promenhatag1=promenhatag[geneclus[:,0]==1]
promenhatag2=promenhatag[(geneclus[:,1]==1)|(geneclus[:,2]==1)]
promenhatag3=promenhatag[geneclus[:,3]==1]
promenhatag4=promenhatag[geneclus[:,4]==1]


cellclus=pd.read_csv(samplename+"/cellcluster_all.txt",sep=",",header=None)
cellclus=cellclus.to_numpy()
cellclus=cellclus[:,0]


##############


tmp=np.load(samplename+"/encodemtx.npy")

encode=tmp.copy()
encode[:,:,0]=tmp[:,:,3]
encode[:,:,1]=tmp[:,:,0]
encode[:,:,2]=tmp[:,:,1]
encode[:,:,3]=tmp[:,:,2]

encodelog=np.log(encode+encode[encode!=0].min())
encoden=encodelog

tmp1=encoden[promenhatag==3]
encoden[promenhatag==3]=(tmp1-tmp1.mean())
tmp1=encoden[promenhatag==4]
encoden[promenhatag==4]=(tmp1-tmp1.mean())
tmp1=encoden[promenhatag==2]
encoden[promenhatag==2]=(tmp1-tmp1.mean())

enhamtxtmp=enhamtx.transpose(0,2,1)
enhamtxm=np.zeros((enhamtxtmp.shape[0],enhamtxtmp.shape[1],4))
enhamtxm[:,:,0]=enhamtxtmp[:,:,cellclus==(1)].mean(axis=2)
enhamtxm[:,:,1]=enhamtxtmp[:,:,((cellclus==(2))|(cellclus==(3)))].mean(axis=2)
enhamtxm[:,:,2]=enhamtxtmp[:,:,cellclus==(4)].mean(axis=2)
enhamtxm[:,:,3]=enhamtxtmp[:,:,cellclus==(5)].mean(axis=2)


#####33


encoden1=encoden[geneclus[:,0]==1]
encoden2=encoden[(geneclus[:,1]==1)|(geneclus[:,2]==1)]
encoden3=encoden[geneclus[:,3]==1]
encoden4=encoden[geneclus[:,4]==1]

enhamtxm1=enhamtxm[geneclus[:,0]==1]
enhamtxm2=enhamtxm[(geneclus[:,1]==1)|(geneclus[:,2]==1)]
enhamtxm3=enhamtxm[geneclus[:,3]==1]
enhamtxm4=enhamtxm[geneclus[:,4]==1]




###########

encoden1m=encoden1[(promenhatag1==3)]
encoden2m=encoden2[(promenhatag2==3)]
encoden3m=encoden3[(promenhatag3==3)]
encoden4m=encoden4[(promenhatag4==3)]

enhamtxm1m=enhamtxm1[(promenhatag1==3)]
enhamtxm2m=enhamtxm2[(promenhatag2==3)]
enhamtxm3m=enhamtxm3[(promenhatag3==3)]
enhamtxm4m=enhamtxm4[(promenhatag4==3)]

tmp13e=encoden1m[:,0]-encoden1m[:,[1,2,3]].mean(axis=1)
tmp23e=encoden2m[:,1]-encoden2m[:,[0,2,3]].mean(axis=1)
tmp33e=encoden3m[:,2]-encoden3m[:,[0,1,3]].mean(axis=1)
tmp43e=encoden4m[:,3]-encoden4m[:,[0,1,2]].mean(axis=1)

tmp13m=enhamtxm1m[:,0]-enhamtxm1m[:,[1,2,3]].mean(axis=1)
tmp23m=enhamtxm2m[:,1]-enhamtxm2m[:,[0,2,3]].mean(axis=1)
tmp33m=enhamtxm3m[:,2]-enhamtxm3m[:,[0,1,3]].mean(axis=1)
tmp43m=enhamtxm4m[:,3]-enhamtxm4m[:,[0,1,2]].mean(axis=1)

##


encoden1m=encoden1[(promenhatag1==2)]
encoden2m=encoden2[(promenhatag2==2)]
encoden3m=encoden3[(promenhatag3==2)]
encoden4m=encoden4[(promenhatag4==2)]

enhamtxm1m=enhamtxm1[(promenhatag1==2)]
enhamtxm2m=enhamtxm2[(promenhatag2==2)]
enhamtxm3m=enhamtxm3[(promenhatag3==2)]
enhamtxm4m=enhamtxm4[(promenhatag4==2)]

tmp12e=encoden1m[:,0]-encoden1m[:,[1,2,3]].mean(axis=1)
tmp22e=encoden2m[:,1]-encoden2m[:,[0,2,3]].mean(axis=1)
tmp32e=encoden3m[:,2]-encoden3m[:,[0,1,3]].mean(axis=1)
tmp42e=encoden4m[:,3]-encoden4m[:,[0,1,2]].mean(axis=1)

tmp12m=enhamtxm1m[:,0]-enhamtxm1m[:,[1,2,3]].mean(axis=1)
tmp22m=enhamtxm2m[:,1]-enhamtxm2m[:,[0,2,3]].mean(axis=1)
tmp32m=enhamtxm3m[:,2]-enhamtxm3m[:,[0,1,3]].mean(axis=1)
tmp42m=enhamtxm4m[:,3]-enhamtxm4m[:,[0,1,2]].mean(axis=1)

##


encoden1m=encoden1[(promenhatag1==4)]
encoden2m=encoden2[(promenhatag2==4)]
encoden3m=encoden3[(promenhatag3==4)]
encoden4m=encoden4[(promenhatag4==4)]

enhamtxm1m=enhamtxm1[(promenhatag1==4)]
enhamtxm2m=enhamtxm2[(promenhatag2==4)]
enhamtxm3m=enhamtxm3[(promenhatag3==4)]
enhamtxm4m=enhamtxm4[(promenhatag4==4)]

tmp14e=encoden1m[:,0]-encoden1m[:,[1,2,3]].mean(axis=1)
tmp24e=encoden2m[:,1]-encoden2m[:,[0,2,3]].mean(axis=1)
tmp34e=encoden3m[:,2]-encoden3m[:,[0,1,3]].mean(axis=1)
tmp44e=encoden4m[:,3]-encoden4m[:,[0,1,2]].mean(axis=1)

tmp14m=enhamtxm1m[:,0]-enhamtxm1m[:,[1,2,3]].mean(axis=1)
tmp24m=enhamtxm2m[:,1]-enhamtxm2m[:,[0,2,3]].mean(axis=1)
tmp34m=enhamtxm3m[:,2]-enhamtxm3m[:,[0,1,3]].mean(axis=1)
tmp44m=enhamtxm4m[:,3]-enhamtxm4m[:,[0,1,2]].mean(axis=1)

##

import scipy.stats as stats

_, encode_p_stroweak = stats.ttest_ind(tmp14e, tmp13e)
_, encode_p_stroprom = stats.ttest_ind(tmp14e, tmp12e)
_, encode_p_weakprom = stats.ttest_ind(tmp13e, tmp12e)
_, enha_p_stroweak = stats.ttest_ind(tmp14m, tmp13m)
_, enha_p_stroprom = stats.ttest_ind(tmp14m, tmp12m)
_, enha_p_weakprom = stats.ttest_ind(tmp13m, tmp12m)

encode_p_stroweak='{:.2e}'.format(encode_p_stroweak)
encode_p_stroprom='{:.2e}'.format(encode_p_stroprom)
encode_p_weakprom='{:.2e}'.format(encode_p_weakprom)
enha_p_stroweak='{:.2e}'.format(enha_p_stroweak)
enha_p_stroprom='{:.2e}'.format(enha_p_stroprom)
enha_p_weakprom='{:.2e}'.format(enha_p_weakprom)

pvalue1=[encode_p_stroweak,encode_p_weakprom,encode_p_stroprom]
pvalue2=[enha_p_stroweak,enha_p_weakprom,enha_p_stroprom]


tmp14epd=pd.DataFrame(tmp14e)
tmp14epd.columns=["H3K27Ac"]
tmp14epd["Class"]="Strong"
tmp13epd=pd.DataFrame(tmp13e)
tmp13epd.columns=["H3K27Ac"]
tmp13epd["Class"]="Weak"
tmp12epd=pd.DataFrame(tmp12e)
tmp12epd.columns=["H3K27Ac"]
tmp12epd["Class"]="Promoter"

Q1 = tmp14epd['H3K27Ac'].quantile(0.25)
Q3 = tmp14epd['H3K27Ac'].quantile(0.75)
IQR = Q3 - Q1
pmin1=Q1 - 1.5 * IQR
pmax1=Q3 + 1.5 * IQR
Q1 = tmp13epd['H3K27Ac'].quantile(0.25)
Q3 = tmp13epd['H3K27Ac'].quantile(0.75)
IQR = Q3 - Q1
pmin2=Q1 - 1.5 * IQR
pmax2=Q3 + 1.5 * IQR
Q1 = tmp12epd['H3K27Ac'].quantile(0.25)
Q3 = tmp12epd['H3K27Ac'].quantile(0.75)
IQR = Q3 - Q1
pmin3=Q1 - 1.5 * IQR
pmax3=Q3 + 1.5 * IQR

pmin=np.array([pmin1,pmin2,pmin3]).min()
pmax=np.array([pmax1,pmax2,pmax3]).max()

median=['{:.3f}'.format(np.median(tmp14e)),'{:.3f}'.format(np.median(tmp13e)),'{:.3f}'.format(np.median(tmp12e))]

df_concat = pd.concat([tmp14epd, tmp13epd, tmp12epd], axis=0)
x = "Class"
y = "H3K27Ac"
order = ['Strong', 'Weak', 'Promoter']
my_pal = {"Strong": "orange", "Weak": "magenta", "Promoter":"gray"}

plt.clf()
plt.figure(figsize=(6.5,6))
plt.rcParams.update({'font.size': 18})
ax = sns.boxplot(data=df_concat, x=x, y=y, order=order,palette=my_pal,showfliers= False)
pairs = [(0, 1),
         (1, 2),
         (0, 2)]

annotatevalue=pvalue1

h=(pmax-pmin)*0.05

for i in range(len(pairs)):
   x1, x2 = pairs[i]
   yh, col = pmax+(h*2)+(i*h*2.5), 'k'
   plt.plot([x1, x1, x2, x2], [yh, yh+h, yh+h, yh], lw=1.5, c=col)
   ax.text((x1+x2)*.5, yh+h, annotatevalue[i], ha='center', va='bottom', color=col)

for xtick in ax.get_xticks():
    ax.text(xtick,pmax+(h*0),median[xtick], 
            horizontalalignment='center',size='x-small',color='k',weight='semibold')

ax.set_xlabel("Class", fontsize=18)
ax.set_ylabel("LogFC of H3K27Ac\n(Monocyte / others)", fontsize=18)
plt.title("Monocyte",fontsize=24)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.ylim(pmin-(h*0.3),pmax+10*h)
plt.subplots_adjust(left=0.20)
plt.subplots_adjust(bottom=0.15)
plt.savefig(samplename+"/barplot_encode_Monocyte_paper.png")


tmp14mpd=pd.DataFrame(tmp14m)
tmp14mpd.columns=["Predicted Activitiy"]
tmp14mpd["Class"]="Strong"
tmp13mpd=pd.DataFrame(tmp13m)
tmp13mpd.columns=["Predicted Activitiy"]
tmp13mpd["Class"]="Weak"
tmp12mpd=pd.DataFrame(tmp12m)
tmp12mpd.columns=["Predicted Activitiy"]
tmp12mpd["Class"]="Promoter"

Q1 = tmp14mpd['Predicted Activitiy'].quantile(0.25)
Q3 = tmp14mpd['Predicted Activitiy'].quantile(0.75)
IQR = Q3 - Q1
pmin1=Q1 - 1.5 * IQR
pmax1=Q3 + 1.5 * IQR
Q1 = tmp13mpd['Predicted Activitiy'].quantile(0.25)
Q3 = tmp13mpd['Predicted Activitiy'].quantile(0.75)
IQR = Q3 - Q1
pmin2=Q1 - 1.5 * IQR
pmax2=Q3 + 1.5 * IQR
Q1 = tmp12mpd['Predicted Activitiy'].quantile(0.25)
Q3 = tmp12mpd['Predicted Activitiy'].quantile(0.75)
IQR = Q3 - Q1
pmin3=Q1 - 1.5 * IQR
pmax3=Q3 + 1.5 * IQR

pmin=np.array([pmin1,pmin2,pmin3]).min()
pmax=np.array([pmax1,pmax2,pmax3]).max()


median=['{:.3f}'.format(np.median(tmp14m)),'{:.3f}'.format(np.median(tmp13m)),'{:.3f}'.format(np.median(tmp12m))]

df_concat = pd.concat([tmp14mpd, tmp13mpd, tmp12mpd], axis=0)
x = "Class"
y = "Predicted Activitiy"
order = ['Strong', 'Weak', 'Promoter']
my_pal = {"Strong": "orange", "Weak": "magenta", "Promoter":"gray"}

plt.clf()
plt.figure(figsize=(6.5,6))
plt.rcParams.update({'font.size': 18})
ax = sns.boxplot(data=df_concat, x=x, y=y, order=order,palette=my_pal,showfliers= False)

pairs = [(0, 1),
         (1, 2),
         (0, 2)]

annotatevalue=pvalue2

h=(pmax-pmin)*0.05

for i in range(len(pairs)):
   x1, x2 = pairs[i]
   yh, col = pmax+(h*3)+(i*h*2.5), 'k'
   plt.plot([x1, x1, x2, x2], [yh, yh+h, yh+h, yh], lw=1.5, c=col)
   ax.text((x1+x2)*.5, yh+h, annotatevalue[i], ha='center', va='bottom', color=col)

for xtick in ax.get_xticks():
    ax.text(xtick,pmax+(h*1),median[xtick], 
            horizontalalignment='center',size='x-small',color='k',weight='semibold')

ax.set_xlabel("Class", fontsize=18)
ax.set_ylabel("Difference of predicted activity\n(Monocyte - others))", fontsize=18)
plt.title("Monocyte",fontsize=24)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.ylim(pmin-(h*0.3),pmax+11*h)
plt.subplots_adjust(left=0.20)
plt.subplots_adjust(bottom=0.15)
plt.savefig(samplename+"/barplot_enha_Monocyte_paper.png")

##


_, encode_p_stroweak = stats.ttest_ind(tmp24e, tmp23e)
_, encode_p_stroprom = stats.ttest_ind(tmp24e, tmp22e)
_, encode_p_weakprom = stats.ttest_ind(tmp23e, tmp22e)
_, enha_p_stroweak = stats.ttest_ind(tmp24m, tmp23m)
_, enha_p_stroprom = stats.ttest_ind(tmp24m, tmp22m)
_, enha_p_weakprom = stats.ttest_ind(tmp23m, tmp22m)

encode_p_stroweak='{:.2e}'.format(encode_p_stroweak)
encode_p_stroprom='{:.2e}'.format(encode_p_stroprom)
encode_p_weakprom='{:.2e}'.format(encode_p_weakprom)
enha_p_stroweak='{:.2e}'.format(enha_p_stroweak)
enha_p_stroprom='{:.2e}'.format(enha_p_stroprom)
enha_p_weakprom='{:.2e}'.format(enha_p_weakprom)

pvalue1=[encode_p_stroweak,encode_p_weakprom,encode_p_stroprom]
pvalue2=[enha_p_stroweak,enha_p_weakprom,enha_p_stroprom]


tmp24epd=pd.DataFrame(tmp24e)
tmp24epd.columns=["H3K27Ac"]
tmp24epd["Class"]="Strong"
tmp23epd=pd.DataFrame(tmp23e)
tmp23epd.columns=["H3K27Ac"]
tmp23epd["Class"]="Weak"
tmp22epd=pd.DataFrame(tmp22e)
tmp22epd.columns=["H3K27Ac"]
tmp22epd["Class"]="Promoter"

Q1 = tmp24epd['H3K27Ac'].quantile(0.25)
Q3 = tmp24epd['H3K27Ac'].quantile(0.75)
IQR = Q3 - Q1
pmin1=Q1 - 1.5 * IQR
pmax1=Q3 + 1.5 * IQR
Q1 = tmp23epd['H3K27Ac'].quantile(0.25)
Q3 = tmp23epd['H3K27Ac'].quantile(0.75)
IQR = Q3 - Q1
pmin2=Q1 - 1.5 * IQR
pmax2=Q3 + 1.5 * IQR
Q1 = tmp22epd['H3K27Ac'].quantile(0.25)
Q3 = tmp22epd['H3K27Ac'].quantile(0.75)
IQR = Q3 - Q1
pmin3=Q1 - 1.5 * IQR
pmax3=Q3 + 1.5 * IQR

pmin=np.array([pmin1,pmin2,pmin3]).min()
pmax=np.array([pmax1,pmax2,pmax3]).max()

median=['{:.3f}'.format(np.median(tmp24e)),'{:.3f}'.format(np.median(tmp23e)),'{:.3f}'.format(np.median(tmp22e))]

df_concat = pd.concat([tmp24epd, tmp23epd, tmp22epd], axis=0)
x = "Class"
y = "H3K27Ac"
order = ['Strong', 'Weak', 'Promoter']
my_pal = {"Strong": "orange", "Weak": "magenta", "Promoter":"gray"}

plt.clf()
plt.figure(figsize=(6.5,6))
plt.rcParams.update({'font.size': 18})
ax = sns.boxplot(data=df_concat, x=x, y=y, order=order,palette=my_pal,showfliers= False)
pairs = [(0, 1),
         (1, 2),
         (0, 2)]

annotatevalue=pvalue1

h=(pmax-pmin)*0.05

for i in range(len(pairs)):
   x1, x2 = pairs[i]
   yh, col = pmax+(h*3)+(i*h*2.5), 'k'
   plt.plot([x1, x1, x2, x2], [yh, yh+h, yh+h, yh], lw=1.5, c=col)
   ax.text((x1+x2)*.5, yh+h, annotatevalue[i], ha='center', va='bottom', color=col)

for xtick in ax.get_xticks():
    ax.text(xtick,pmax+(h*1),median[xtick], 
            horizontalalignment='center',size='x-small',color='k',weight='semibold')

ax.set_xlabel("Class", fontsize=18)
ax.set_ylabel("LogFC of H3K27Ac\n(CD4 T cell / others)", fontsize=18)
plt.title("CD4 T cell",fontsize=24)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.ylim(pmin-(h*0.3),pmax+11*h)
plt.subplots_adjust(left=0.20)
plt.subplots_adjust(bottom=0.15)
plt.savefig(samplename+"/barplot_encode_CD4Tcell_paper.png")


tmp24mpd=pd.DataFrame(tmp24m)
tmp24mpd.columns=["Predicted Activitiy"]
tmp24mpd["Class"]="Strong"
tmp23mpd=pd.DataFrame(tmp23m)
tmp23mpd.columns=["Predicted Activitiy"]
tmp23mpd["Class"]="Weak"
tmp22mpd=pd.DataFrame(tmp22m)
tmp22mpd.columns=["Predicted Activitiy"]
tmp22mpd["Class"]="Promoter"

Q1 = tmp24mpd['Predicted Activitiy'].quantile(0.25)
Q3 = tmp24mpd['Predicted Activitiy'].quantile(0.75)
IQR = Q3 - Q1
pmin1=Q1 - 1.5 * IQR
pmax1=Q3 + 1.5 * IQR
Q1 = tmp23mpd['Predicted Activitiy'].quantile(0.25)
Q3 = tmp23mpd['Predicted Activitiy'].quantile(0.75)
IQR = Q3 - Q1
pmin2=Q1 - 1.5 * IQR
pmax2=Q3 + 1.5 * IQR
Q1 = tmp22mpd['Predicted Activitiy'].quantile(0.25)
Q3 = tmp22mpd['Predicted Activitiy'].quantile(0.75)
IQR = Q3 - Q1
pmin3=Q1 - 1.5 * IQR
pmax3=Q3 + 1.5 * IQR

pmin=np.array([pmin1,pmin2,pmin3]).min()
pmax=np.array([pmax1,pmax2,pmax3]).max()


median=['{:.3f}'.format(np.median(tmp24m)),'{:.3f}'.format(np.median(tmp23m)),'{:.3f}'.format(np.median(tmp22m))]

df_concat = pd.concat([tmp24mpd, tmp23mpd, tmp22mpd], axis=0)
x = "Class"
y = "Predicted Activitiy"
order = ['Strong', 'Weak', 'Promoter']
my_pal = {"Strong": "orange", "Weak": "magenta", "Promoter":"gray"}

plt.clf()
plt.figure(figsize=(6.5,6))
plt.rcParams.update({'font.size': 18})
ax = sns.boxplot(data=df_concat, x=x, y=y, order=order,palette=my_pal,showfliers= False)

pairs = [(0, 1),
         (1, 2),
         (0, 2)]

annotatevalue=pvalue2

h=(pmax-pmin)*0.05

for i in range(len(pairs)):
   x1, x2 = pairs[i]
   yh, col = pmax+(h*3)+(i*h*2.5), 'k'
   plt.plot([x1, x1, x2, x2], [yh, yh+h, yh+h, yh], lw=1.5, c=col)
   ax.text((x1+x2)*.5, yh+h, annotatevalue[i], ha='center', va='bottom', color=col)

for xtick in ax.get_xticks():
    ax.text(xtick,pmax+(h*1),median[xtick], 
            horizontalalignment='center',size='x-small',color='k',weight='semibold')

ax.set_xlabel("Class", fontsize=18)
ax.set_ylabel("Difference of predicted activity\n(CD4 Tcell - others))", fontsize=18)
plt.title("CD4 T cell",fontsize=24)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.ylim(pmin-(h*0.3),pmax+11*h)
plt.subplots_adjust(left=0.20)
plt.subplots_adjust(bottom=0.15)
plt.savefig(samplename+"/barplot_enha_CD4Tcell_paper.png")



##



_, encode_p_stroweak = stats.ttest_ind(tmp34e, tmp33e)
_, encode_p_stroprom = stats.ttest_ind(tmp34e, tmp32e)
_, encode_p_weakprom = stats.ttest_ind(tmp33e, tmp32e)
_, enha_p_stroweak = stats.ttest_ind(tmp34m, tmp33m)
_, enha_p_stroprom = stats.ttest_ind(tmp34m, tmp32m)
_, enha_p_weakprom = stats.ttest_ind(tmp33m, tmp32m)

encode_p_stroweak='{:.2e}'.format(encode_p_stroweak)
encode_p_stroprom='{:.2e}'.format(encode_p_stroprom)
encode_p_weakprom='{:.2e}'.format(encode_p_weakprom)
enha_p_stroweak='{:.2e}'.format(enha_p_stroweak)
enha_p_stroprom='{:.2e}'.format(enha_p_stroprom)
enha_p_weakprom='{:.2e}'.format(enha_p_weakprom)

pvalue1=[encode_p_stroweak,encode_p_weakprom,encode_p_stroprom]
pvalue2=[enha_p_stroweak,enha_p_weakprom,enha_p_stroprom]


tmp34epd=pd.DataFrame(tmp34e)
tmp34epd.columns=["H3K27Ac"]
tmp34epd["Class"]="Strong"
tmp33epd=pd.DataFrame(tmp33e)
tmp33epd.columns=["H3K27Ac"]
tmp33epd["Class"]="Weak"
tmp32epd=pd.DataFrame(tmp32e)
tmp32epd.columns=["H3K27Ac"]
tmp32epd["Class"]="Promoter"

Q1 = tmp34epd['H3K27Ac'].quantile(0.25)
Q3 = tmp34epd['H3K27Ac'].quantile(0.75)
IQR = Q3 - Q1
pmin1=Q1 - 1.5 * IQR
pmax1=Q3 + 1.5 * IQR
Q1 = tmp33epd['H3K27Ac'].quantile(0.25)
Q3 = tmp33epd['H3K27Ac'].quantile(0.75)
IQR = Q3 - Q1
pmin2=Q1 - 1.5 * IQR
pmax2=Q3 + 1.5 * IQR
Q1 = tmp32epd['H3K27Ac'].quantile(0.25)
Q3 = tmp32epd['H3K27Ac'].quantile(0.75)
IQR = Q3 - Q1
pmin3=Q1 - 1.5 * IQR
pmax3=Q3 + 1.5 * IQR

pmin=np.array([pmin1,pmin2,pmin3]).min()
pmax=np.array([pmax1,pmax2,pmax3]).max()

median=['{:.3f}'.format(np.median(tmp34e)),'{:.3f}'.format(np.median(tmp33e)),'{:.3f}'.format(np.median(tmp32e))]

df_concat = pd.concat([tmp34epd, tmp33epd, tmp32epd], axis=0)
x = "Class"
y = "H3K27Ac"
order = ['Strong', 'Weak', 'Promoter']
my_pal = {"Strong": "orange", "Weak": "magenta", "Promoter":"gray"}

plt.clf()
plt.figure(figsize=(6.5,6))
plt.rcParams.update({'font.size': 18})
ax = sns.boxplot(data=df_concat, x=x, y=y, order=order,palette=my_pal,showfliers= False)
pairs = [(0, 1),
         (1, 2),
         (0, 2)]

annotatevalue=pvalue1

h=(pmax-pmin)*0.05

for i in range(len(pairs)):
   x1, x2 = pairs[i]
   yh, col = pmax+(h*3)+(i*h*2.5), 'k'
   plt.plot([x1, x1, x2, x2], [yh, yh+h, yh+h, yh], lw=1.5, c=col)
   ax.text((x1+x2)*.5, yh+h, annotatevalue[i], ha='center', va='bottom', color=col)

for xtick in ax.get_xticks():
    ax.text(xtick,pmax+(h*1),median[xtick], 
            horizontalalignment='center',size='x-small',color='k',weight='semibold')

ax.set_xlabel("Class", fontsize=18)
ax.set_ylabel("LogFC of H3K27Ac\n(CD8 T cell / others)", fontsize=18)
plt.title("CD8 T cell",fontsize=24)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.ylim(pmin-(h*0.3),pmax+11*h)
plt.subplots_adjust(left=0.20)
plt.subplots_adjust(bottom=0.15)
plt.savefig(samplename+"/barplot_encode_CD8Tcell_paper.png")


tmp34mpd=pd.DataFrame(tmp34m)
tmp34mpd.columns=["Predicted Activitiy"]
tmp34mpd["Class"]="Strong"
tmp33mpd=pd.DataFrame(tmp33m)
tmp33mpd.columns=["Predicted Activitiy"]
tmp33mpd["Class"]="Weak"
tmp32mpd=pd.DataFrame(tmp32m)
tmp32mpd.columns=["Predicted Activitiy"]
tmp32mpd["Class"]="Promoter"

Q1 = tmp34mpd['Predicted Activitiy'].quantile(0.25)
Q3 = tmp34mpd['Predicted Activitiy'].quantile(0.75)
IQR = Q3 - Q1
pmin1=Q1 - 1.5 * IQR
pmax1=Q3 + 1.5 * IQR
Q1 = tmp33mpd['Predicted Activitiy'].quantile(0.25)
Q3 = tmp33mpd['Predicted Activitiy'].quantile(0.75)
IQR = Q3 - Q1
pmin2=Q1 - 1.5 * IQR
pmax2=Q3 + 1.5 * IQR
Q1 = tmp32mpd['Predicted Activitiy'].quantile(0.25)
Q3 = tmp32mpd['Predicted Activitiy'].quantile(0.75)
IQR = Q3 - Q1
pmin3=Q1 - 1.5 * IQR
pmax3=Q3 + 1.5 * IQR

pmin=np.array([pmin1,pmin2,pmin3]).min()
pmax=np.array([pmax1,pmax2,pmax3]).max()


median=['{:.3f}'.format(np.median(tmp34m)),'{:.3f}'.format(np.median(tmp33m)),'{:.3f}'.format(np.median(tmp32m))]

df_concat = pd.concat([tmp34mpd, tmp33mpd, tmp32mpd], axis=0)
x = "Class"
y = "Predicted Activitiy"
order = ['Strong', 'Weak', 'Promoter']
my_pal = {"Strong": "orange", "Weak": "magenta", "Promoter":"gray"}

plt.clf()
plt.figure(figsize=(6.5,6))
plt.rcParams.update({'font.size': 18})
ax = sns.boxplot(data=df_concat, x=x, y=y, order=order,palette=my_pal,showfliers= False)

pairs = [(0, 1),
         (1, 2),
         (0, 2)]

annotatevalue=pvalue2

h=(pmax-pmin)*0.05

for i in range(len(pairs)):
   x1, x2 = pairs[i]
   yh, col = pmax+(h*3)+(i*h*2.5), 'k'
   plt.plot([x1, x1, x2, x2], [yh, yh+h, yh+h, yh], lw=1.5, c=col)
   ax.text((x1+x2)*.5, yh+h, annotatevalue[i], ha='center', va='bottom', color=col)

for xtick in ax.get_xticks():
    ax.text(xtick,pmax+(h*1),median[xtick], 
            horizontalalignment='center',size='x-small',color='k',weight='semibold')

ax.set_xlabel("Class", fontsize=18)
ax.set_ylabel("Difference of predicted activity\n(CD8 Tcell - others))", fontsize=18)
plt.title("CD8 T cell",fontsize=24)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.ylim(pmin-(h*0.3),pmax+11*h)
plt.subplots_adjust(left=0.20)
plt.subplots_adjust(bottom=0.15)
plt.savefig(samplename+"/barplot_enha_CD8Tcell_paper.png")





##




_, encode_p_stroweak = stats.ttest_ind(tmp44e, tmp43e)
_, encode_p_stroprom = stats.ttest_ind(tmp44e, tmp42e)
_, encode_p_weakprom = stats.ttest_ind(tmp43e, tmp42e)
_, enha_p_stroweak = stats.ttest_ind(tmp44m, tmp43m)
_, enha_p_stroprom = stats.ttest_ind(tmp44m, tmp42m)
_, enha_p_weakprom = stats.ttest_ind(tmp43m, tmp42m)

encode_p_stroweak='{:.2e}'.format(encode_p_stroweak)
encode_p_stroprom='{:.2e}'.format(encode_p_stroprom)
encode_p_weakprom='{:.2e}'.format(encode_p_weakprom)
enha_p_stroweak='{:.2e}'.format(enha_p_stroweak)
enha_p_stroprom='{:.2e}'.format(enha_p_stroprom)
enha_p_weakprom='{:.2e}'.format(enha_p_weakprom)

pvalue1=[encode_p_stroweak,encode_p_weakprom,encode_p_stroprom]
pvalue2=[enha_p_stroweak,enha_p_weakprom,enha_p_stroprom]


tmp44epd=pd.DataFrame(tmp44e)
tmp44epd.columns=["H3K27Ac"]
tmp44epd["Class"]="Strong"
tmp43epd=pd.DataFrame(tmp43e)
tmp43epd.columns=["H3K27Ac"]
tmp43epd["Class"]="Weak"
tmp42epd=pd.DataFrame(tmp42e)
tmp42epd.columns=["H3K27Ac"]
tmp42epd["Class"]="Promoter"

Q1 = tmp44epd['H3K27Ac'].quantile(0.25)
Q3 = tmp44epd['H3K27Ac'].quantile(0.75)
IQR = Q3 - Q1
pmin1=Q1 - 1.5 * IQR
pmax1=Q3 + 1.5 * IQR
Q1 = tmp43epd['H3K27Ac'].quantile(0.25)
Q3 = tmp43epd['H3K27Ac'].quantile(0.75)
IQR = Q3 - Q1
pmin2=Q1 - 1.5 * IQR
pmax2=Q3 + 1.5 * IQR
Q1 = tmp42epd['H3K27Ac'].quantile(0.25)
Q3 = tmp42epd['H3K27Ac'].quantile(0.75)
IQR = Q3 - Q1
pmin3=Q1 - 1.5 * IQR
pmax3=Q3 + 1.5 * IQR

pmin=np.array([pmin1,pmin2,pmin3]).min()
pmax=np.array([pmax1,pmax2,pmax3]).max()

median=['{:.3f}'.format(np.median(tmp44e)),'{:.3f}'.format(np.median(tmp43e)),'{:.3f}'.format(np.median(tmp42e))]

df_concat = pd.concat([tmp44epd, tmp43epd, tmp42epd], axis=0)
x = "Class"
y = "H3K27Ac"
order = ['Strong', 'Weak', 'Promoter']
my_pal = {"Strong": "orange", "Weak": "magenta", "Promoter":"gray"}

plt.clf()
plt.figure(figsize=(6.5,6))
plt.rcParams.update({'font.size': 18})
ax = sns.boxplot(data=df_concat, x=x, y=y, order=order,palette=my_pal,showfliers= False)
pairs = [(0, 1),
         (1, 2),
         (0, 2)]

annotatevalue=pvalue1

h=(pmax-pmin)*0.05

for i in range(len(pairs)):
   x1, x2 = pairs[i]
   yh, col = pmax+(h*3)+(i*h*2.5), 'k'
   plt.plot([x1, x1, x2, x2], [yh, yh+h, yh+h, yh], lw=1.5, c=col)
   ax.text((x1+x2)*.5, yh+h, annotatevalue[i], ha='center', va='bottom', color=col)

for xtick in ax.get_xticks():
    ax.text(xtick,pmax+(h*1),median[xtick], 
            horizontalalignment='center',size='x-small',color='k',weight='semibold')

ax.set_xlabel("Class", fontsize=18)
ax.set_ylabel("LogFC of H3K27Ac\n(B cell / others)", fontsize=18)
plt.title("B cell",fontsize=24)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.ylim(pmin-(h*0.3),pmax+11*h)
plt.subplots_adjust(left=0.20)
plt.subplots_adjust(bottom=0.15)
plt.savefig(samplename+"/barplot_encode_Bcell_paper.png")


tmp44mpd=pd.DataFrame(tmp44m)
tmp44mpd.columns=["Predicted Activitiy"]
tmp44mpd["Class"]="Strong"
tmp43mpd=pd.DataFrame(tmp43m)
tmp43mpd.columns=["Predicted Activitiy"]
tmp43mpd["Class"]="Weak"
tmp42mpd=pd.DataFrame(tmp42m)
tmp42mpd.columns=["Predicted Activitiy"]
tmp42mpd["Class"]="Promoter"

Q1 = tmp44mpd['Predicted Activitiy'].quantile(0.25)
Q3 = tmp44mpd['Predicted Activitiy'].quantile(0.75)
IQR = Q3 - Q1
pmin1=Q1 - 1.5 * IQR
pmax1=Q3 + 1.5 * IQR
Q1 = tmp43mpd['Predicted Activitiy'].quantile(0.25)
Q3 = tmp43mpd['Predicted Activitiy'].quantile(0.75)
IQR = Q3 - Q1
pmin2=Q1 - 1.5 * IQR
pmax2=Q3 + 1.5 * IQR
Q1 = tmp42mpd['Predicted Activitiy'].quantile(0.25)
Q3 = tmp42mpd['Predicted Activitiy'].quantile(0.75)
IQR = Q3 - Q1
pmin3=Q1 - 1.5 * IQR
pmax3=Q3 + 1.5 * IQR

pmin=np.array([pmin1,pmin2,pmin3]).min()
pmax=np.array([pmax1,pmax2,pmax3]).max()


median=['{:.3f}'.format(np.median(tmp44m)),'{:.3f}'.format(np.median(tmp43m)),'{:.3f}'.format(np.median(tmp42m))]

df_concat = pd.concat([tmp44mpd, tmp43mpd, tmp42mpd], axis=0)
x = "Class"
y = "Predicted Activitiy"
order = ['Strong', 'Weak', 'Promoter']
my_pal = {"Strong": "orange", "Weak": "magenta", "Promoter":"gray"}

plt.clf()
plt.figure(figsize=(6.5,6))
plt.rcParams.update({'font.size': 18})
ax = sns.boxplot(data=df_concat, x=x, y=y, order=order,palette=my_pal,showfliers= False)

pairs = [(0, 1),
         (1, 2),
         (0, 2)]

annotatevalue=pvalue2

h=(pmax-pmin)*0.05

for i in range(len(pairs)):
   x1, x2 = pairs[i]
   yh, col = pmax+(h*3)+(i*h*2.5), 'k'
   plt.plot([x1, x1, x2, x2], [yh, yh+h, yh+h, yh], lw=1.5, c=col)
   ax.text((x1+x2)*.5, yh+h, annotatevalue[i], ha='center', va='bottom', color=col)

for xtick in ax.get_xticks():
    ax.text(xtick,pmax+(h*1),median[xtick], 
            horizontalalignment='center',size='x-small',color='k',weight='semibold')

ax.set_xlabel("Class", fontsize=18)
ax.set_ylabel("Difference of predicted activity\n(Bcell - others))", fontsize=18)
plt.title("B cell",fontsize=24)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.ylim(pmin-(h*0.3),pmax+11*h)
plt.subplots_adjust(left=0.20)
plt.subplots_adjust(bottom=0.15)
plt.savefig(samplename+"/barplot_enha_Bcell_paper.png")




########################################
############################################
###############################################






##########################
from sklearn.cluster import KMeans
from sklearn.metrics.cluster import adjusted_rand_score
ATAC_corr=np.load("pbmc_fold1/ATAC_corr.npy")

ari_origact=np.zeros((5))
ari_scenicact=np.zeros((5))

fnamelist1=["pbmc_fold1","pbmc_fold2","pbmc_fold3","pbmc_fold4","pbmc_fold5"]
fnamelist2=["scenicplus0","scenicplus4","scenicplus7","scenicplus2","scenicplus8"]

for i in range(5):
   print(i)
   fnametmp1=fnamelist1[i]
   fnametmp2=fnamelist2[i]
   enhamtx=np.load(fnametmp1+"/Deeplift_full_ver2_all.npy")
   grad=np.load("pbmc_fold1/allgrad_ssep_max.npy")
   scout4=pd.read_csv("scenicplus/"+fnametmp2+"/scout4.csv")
   thv=np.percentile(grad[grad!=0], [40])
   sthv=np.percentile(grad[grad!=0], [85])
   promenhatag=np.zeros(grad.shape)
   promenhatag[(grad>thv)&(ATAC_corr>0)]=2
   promenhatag[(grad>sthv)&(ATAC_corr>0)]=4
   promenhatag[:,0]=3
   cellclus=pd.read_csv("pbmc_after/cellcluster_all.txt",sep=",",header=None)
   cellid=pd.read_csv("pbmc_after/celllist_scVI.csv",sep="\t",header=None)
   cellmtx=pd.concat([cellid,cellclus],axis=1)
   cellmtx.columns=["Cellid","Cluster"]
   cellmtx["tag"]=np.array(range(cellmtx.shape[0]))
   scout4s=scout4[scout4["Cellid"].isin(cellmtx["Cellid"])]
   scouttag=pd.merge(scout4s, cellmtx, on='Cellid', how='left')
   scouttag=scouttag["tag"].to_numpy().astype(int)
   cellmtxs=cellmtx.iloc[scouttag,:]
   cellcluss=cellmtxs["Cluster"].to_numpy()
   enhamtxuse=enhamtx[:,scouttag,:].transpose(0,2,1)
   enhamtxusestro=enhamtxuse[promenhatag==4].transpose(1,0)
   enhadf=pd.DataFrame(enhamtxusestro)
   kmeans_model = KMeans(n_clusters = 9).fit(enhadf) 
   kmeans_modelsc = KMeans(n_clusters = 9).fit(scout4s.iloc[:,1:(scout4.shape[1]-2)]) 
   ari_origact[i]=adjusted_rand_score(cellcluss, kmeans_model.labels_)
   ari_scenicact[i]=adjusted_rand_score(cellcluss, kmeans_modelsc.labels_)

np.save("pbmc_after/ari_origact.npy",ari_origact)
np.save("pbmc_after/ari_scenicact.npy",ari_scenicact)



##########################
from sklearn.cluster import KMeans
from sklearn.metrics.cluster import adjusted_rand_score


ari_origact=np.zeros((5))
ari_scenicact=np.zeros((5))

fnamelist1=["pbmcnew/try1","pbmcnew/try2","pbmcnew/try3","pbmcnew/try4","pbmcnew/try5"]
fnamelist2=["scenicplus0","scenicplus4","scenicplus7","scenicplus2","scenicplus8"]

for i in range(5):
   print(i)
   fnametmp1=fnamelist1[i]
   fnametmp2=fnamelist2[i]
   enhamtxusestro=np.load(fnametmp1+"/enhamtxusestro.npy")
   scout4=pd.read_csv("scenicplus/"+fnametmp2+"/scout4.csv")
   cellclus=pd.read_csv("pbmc_after/cellcluster_all.txt",sep=",",header=None)
   cellid=pd.read_csv("pbmc_after/celllist_scVI.csv",sep="\t",header=None)
   cellmtx=pd.concat([cellid,cellclus],axis=1)
   cellmtx.columns=["Cellid","Cluster"]
   cellmtx["tag"]=np.array(range(cellmtx.shape[0]))
   scout4s=scout4[scout4["Cellid"].isin(cellmtx["Cellid"])]
   scouttag=pd.merge(scout4s, cellmtx, on='Cellid', how='left')
   scouttag=scouttag["tag"].to_numpy().astype(int)
   cellmtxs=cellmtx.iloc[scouttag,:]
   cellcluss=cellmtxs["Cluster"].to_numpy()
   enhamtxusestro=enhamtxusestro[scouttag,:]
   enhadf=pd.DataFrame(enhamtxusestro)
   kmeans_model = KMeans(n_clusters = 9).fit(enhadf) 
   kmeans_modelsc = KMeans(n_clusters = 9).fit(scout4s.iloc[:,1:(scout4.shape[1]-2)]) 
   ari_origact[i]=adjusted_rand_score(cellcluss, kmeans_model.labels_)
   ari_scenicact[i]=adjusted_rand_score(cellcluss, kmeans_modelsc.labels_)

np.save("pbmc_after/ari_origact.npy",ari_origact)
np.save("pbmc_after/ari_scenicact.npy",ari_scenicact)

##########################
from sklearn.cluster import KMeans
from sklearn.metrics.cluster import adjusted_rand_score


ari_origact_weak=np.zeros((5))

fnamelist1=["pbmcnew/try1","pbmcnew/try2","pbmcnew/try3","pbmcnew/try4","pbmcnew/try5"]
fnamelist2=["scenicplus0","scenicplus4","scenicplus7","scenicplus2","scenicplus8"]

for i in range(5):
   print(i)
   fnametmp1=fnamelist1[i]
   fnametmp2=fnamelist2[i]
   enhamtxusestro=np.load(fnametmp1+"/enhamtxuseweak.npy")
   scout4=pd.read_csv("scenicplus/"+fnametmp2+"/scout4.csv")
   cellclus=pd.read_csv("pbmc_after/cellcluster_all.txt",sep=",",header=None)
   cellid=pd.read_csv("pbmc_after/celllist_scVI.csv",sep="\t",header=None)
   cellmtx=pd.concat([cellid,cellclus],axis=1)
   cellmtx.columns=["Cellid","Cluster"]
   cellmtx["tag"]=np.array(range(cellmtx.shape[0]))
   scout4s=scout4[scout4["Cellid"].isin(cellmtx["Cellid"])]
   scouttag=pd.merge(scout4s, cellmtx, on='Cellid', how='left')
   scouttag=scouttag["tag"].to_numpy().astype(int)
   cellmtxs=cellmtx.iloc[scouttag,:]
   cellcluss=cellmtxs["Cluster"].to_numpy()
   enhamtxusestro=enhamtxusestro[scouttag,:]
   enhadf=pd.DataFrame(enhamtxusestro)
   kmeans_model = KMeans(n_clusters = 9).fit(enhadf) 
   ari_origact_weak[i]=adjusted_rand_score(cellcluss, kmeans_model.labels_)

np.save("pbmc_after/ari_origact_weak.npy",ari_origact_weak)


#####

from sklearn.cluster import KMeans
from sklearn.metrics.cluster import adjusted_rand_score


ari_origact_prom=np.zeros((5))

fnamelist1=["pbmcnew/try1","pbmcnew/try2","pbmcnew/try3","pbmcnew/try4","pbmcnew/try5"]
fnamelist2=["scenicplus0","scenicplus4","scenicplus7","scenicplus2","scenicplus8"]

for i in range(5):
   print(i)
   fnametmp1=fnamelist1[i]
   fnametmp2=fnamelist2[i]
   enhamtxusestro=np.load(fnametmp1+"/enhamtxuseprom.npy")
   scout4=pd.read_csv("scenicplus/"+fnametmp2+"/scout4.csv")
   cellclus=pd.read_csv("pbmc_after/cellcluster_all.txt",sep=",",header=None)
   cellid=pd.read_csv("pbmc_after/celllist_scVI.csv",sep="\t",header=None)
   cellmtx=pd.concat([cellid,cellclus],axis=1)
   cellmtx.columns=["Cellid","Cluster"]
   cellmtx["tag"]=np.array(range(cellmtx.shape[0]))
   scout4s=scout4[scout4["Cellid"].isin(cellmtx["Cellid"])]
   scouttag=pd.merge(scout4s, cellmtx, on='Cellid', how='left')
   scouttag=scouttag["tag"].to_numpy().astype(int)
   cellmtxs=cellmtx.iloc[scouttag,:]
   cellcluss=cellmtxs["Cluster"].to_numpy()
   enhamtxusestro=enhamtxusestro[scouttag,:]
   enhadf=pd.DataFrame(enhamtxusestro)
   kmeans_model = KMeans(n_clusters = 9).fit(enhadf) 
   ari_origact_prom[i]=adjusted_rand_score(cellcluss, kmeans_model.labels_)

np.save("pbmc_after/ari_origact_prom.npy",ari_origact_prom)

#####


stats.ttest_ind(ari_origact, ari_origact_weak)

ari_origact=np.load("pbmc_after/ari_origact.npy")


######

cellclus=pd.read_csv(samplename+"/cellcluster_all.txt",sep=",",header=None)
cellclus=cellclus.to_numpy()
cellclus=cellclus[:,0]



grad_archr=np.load("grad_archr.npy")
grad_cicero=np.load("grad_cicero.npy")
grad_direct=np.load("grad_direct.npy")
grad_scenic=np.load("grad_scenic.npy")
grad_orig=np.load("grad_orig.npy")


ATAC_use=np.load("pbmc_after/atac_count.npy")

ATAC_mtxuse=np.zeros((grad_orig.shape[0],grad_orig.shape[1],ATAC_use.shape[1]))
for j in range(ATAC_mtxuse.shape[0]):
  print(j)
  peakstart=int(pairlist[j,0])
  peakend=int(pairlist[j,1])
  peaknum_gene=peakend-peakstart+1
  prompos=int(pairlist_prom[j,0])
  enha=list(range(peakstart,(peakend+1)))
  enha.remove(prompos)
  tmp1=ATAC_use[prompos,:]
  tmp2=ATAC_use[enha,:]
  ATAC_mtxuse[j,0,:]=tmp1
  ATAC_mtxuse[j,1:peaknum_gene,:]=tmp2


enhatag=np.load("pbmc_fold1/enha.npy")
enhatag[:,0]=(-1)


ATACcluster_orig=np.zeros((cellclus.shape[0],5))
ATACcluster_archr=np.zeros((cellclus.shape[0],5))
ATACcluster_cicero=np.zeros((cellclus.shape[0],5))
ATACcluster_scenic=np.zeros((cellclus.shape[0],5))
ATACcluster_direct=np.zeros((cellclus.shape[0],5))

for i in range(5):
   enhath_orig=np.percentile(grad_orig[enhatag!=(-1),i], [85])
   enhath_archr=np.percentile(grad_archr[enhatag!=(-1),i], [85])
   enhath_cicero=np.percentile(grad_cicero[enhatag!=(-1),i], [85])
   enhath_scenic=np.percentile(grad_scenic[enhatag!=(-1),i], [85])
   enhath_direct=np.percentile(grad_direct[enhatag!=(-1),i], [85])
   enhatag_orig=grad_orig[:,:,i]>enhath_orig
   enhatag_archr=grad_archr[:,:,i]>enhath_archr
   enhatag_cicero=grad_cicero[:,:,i]>enhath_cicero
   enhatag_scenic=grad_scenic[:,:,i]>enhath_scenic
   enhatag_direct=grad_direct[:,:,i]>enhath_direct
   enhatag_orig[:,0]=0
   enhatag_archr[:,0]=0
   enhatag_cicero[:,0]=0
   enhatag_scenic[:,0]=0
   enhatag_direct[:,0]=0
   ATACtmp_orig=ATAC_mtxuse[enhatag_orig==1].transpose(1,0)
   ATACtmp_archr=ATAC_mtxuse[enhatag_archr==1].transpose(1,0)
   ATACtmp_cicero=ATAC_mtxuse[enhatag_cicero==1].transpose(1,0)
   ATACtmp_scenic=ATAC_mtxuse[enhatag_scenic==1].transpose(1,0)
   ATACtmp_direct=ATAC_mtxuse[enhatag_direct==1].transpose(1,0)
   ATACtmp_orig_df=pd.DataFrame(ATACtmp_orig)
   ATACtmp_archr_df=pd.DataFrame(ATACtmp_archr)
   ATACtmp_cicero_df=pd.DataFrame(ATACtmp_cicero)
   ATACtmp_scenic_df=pd.DataFrame(ATACtmp_scenic)
   ATACtmp_direct_df=pd.DataFrame(ATACtmp_direct)
   ATACtmp_orig_df_model = KMeans(n_clusters = 9).fit(ATACtmp_orig_df) 
   ATACtmp_archr_df_model = KMeans(n_clusters = 9).fit(ATACtmp_archr_df) 
   ATACtmp_cicero_df_model = KMeans(n_clusters = 9).fit(ATACtmp_cicero_df) 
   ATACtmp_scenic_df_model = KMeans(n_clusters = 9).fit(ATACtmp_scenic_df) 
   ATACtmp_direct_df_model = KMeans(n_clusters = 9).fit(ATACtmp_direct_df)
   ATACcluster_orig[:,i]=ATACtmp_orig_df_model.labels_
   ATACcluster_archr[:,i]=ATACtmp_archr_df_model.labels_
   ATACcluster_cicero[:,i]=ATACtmp_cicero_df_model.labels_
   ATACcluster_scenic[:,i]=ATACtmp_scenic_df_model.labels_
   ATACcluster_direct[:,i]=ATACtmp_direct_df_model.labels_


np.save("pbmc_after/ATACcluster_orig.npy",ATACcluster_orig)
np.save("pbmc_after/ATACcluster_archr.npy",ATACcluster_archr)
np.save("pbmc_after/ATACcluster_cicero.npy",ATACcluster_cicero)
np.save("pbmc_after/ATACcluster_scenic.npy",ATACcluster_scenic)
np.save("pbmc_after/ATACcluster_direct.npy",ATACcluster_direct)


ATACtmp_all=ATAC_mtxuse[enhatag!=(-1)].transpose(1,0)
ATACtmp_all_df=pd.DataFrame(ATACtmp_all)
ATACtmp_all_df_model = KMeans(n_clusters = 9).fit(ATACtmp_all_df) 
ari_all=adjusted_rand_score(cellclus, ATACtmp_all_df_model.labels_)


ari_orig=np.zeros((5))
ari_archr=np.zeros((5))
ari_cicero=np.zeros((5))
ari_scenic=np.zeros((5))
ari_direct=np.zeros((5))

for i in range(5):
   ari_orig[i]=adjusted_rand_score(cellclus, ATACcluster_orig[:,i])
   ari_archr[i]=adjusted_rand_score(cellclus, ATACcluster_archr[:,i])
   ari_cicero[i]=adjusted_rand_score(cellclus, ATACcluster_cicero[:,i])
   ari_scenic[i]=adjusted_rand_score(cellclus, ATACcluster_scenic[:,i])
   ari_direct[i]=adjusted_rand_score(cellclus, ATACcluster_direct[:,i])

ari_orig=adjusted_rand_score(cellclus, ATACtmp_orig_df_model.labels_)
ari_archr=adjusted_rand_score(cellclus, ATACtmp_archr_df_model.labels_)
ari_cicero=adjusted_rand_score(cellclus, ATACtmp_cicero_df_model.labels_)
ari_scenic=adjusted_rand_score(cellclus, ATACtmp_scenic_df_model.labels_)
ari_direct=adjusted_rand_score(cellclus, ATACtmp_direct_df_model.labels_)

left = np.array(["Our model", "Scenic+", "DIRECT-NET", "ArchR", "Cicero"])
height = np.array([ari_orig, ari_scenic,ari_direct,ari_archr,ari_cicero])

points = (ari_orig, ari_scenic,ari_direct,ari_archr,ari_cicero)
fig, ax = plt.subplots(figsize=(6, 6.5))
ax.bar(left,height)
ax.set_xticklabels(["Our model", "Scenic+", "DIRECT-NET", "ArchR", "Cicero"], fontsize=16,rotation=45,fontweight='bold')
plt.tick_params(labelsize=16,width=2)
ax.spines["top"].set_linewidth(2)
ax.spines["left"].set_linewidth(2)
ax.spines["bottom"].set_linewidth(2)
ax.spines["right"].set_linewidth(2)
plt.title("Clutering", fontsize=28)
#plt.ylim(0.43,0.68)
ax.set_xlabel("Method", fontsize=20)
ax.set_ylabel("Adjusted Rand Index", fontsize=20)
plt.subplots_adjust(left=0.18)
plt.subplots_adjust(bottom=0.25)
plt.savefig("pbmc_after/ARI.png")
plt.show()





points = (ari_orig, ari_scenic,ari_direct,ari_archr,ari_cicero)
fig, ax = plt.subplots(figsize=(6, 6.5))
bp = ax.boxplot(points)
ax.set_xticklabels(["Our model", "Scenic+", "DIRECT-NET", "ArchR", "Cicero"], fontsize=16,rotation=45,fontweight='bold')
plt.tick_params(labelsize=16,width=2)
ax.spines["top"].set_linewidth(2)
ax.spines["left"].set_linewidth(2)
ax.spines["bottom"].set_linewidth(2)
ax.spines["right"].set_linewidth(2)
plt.title("Clustering", fontsize=28)
plt.ylim(0.43,0.68)
ax.set_xlabel("Method", fontsize=20)
ax.set_ylabel("Adjusted Rand Index", fontsize=20)
plt.subplots_adjust(left=0.18)
plt.subplots_adjust(bottom=0.25)
plt.savefig("pbmc_after/ARI.png")
plt.show()


ari_origactpd=pd.DataFrame(ari_origact)
ari_origactpd.columns=["ARI"]
ari_origactpd["Method"]="Our model(Activity)"
ari_scenicactpd=pd.DataFrame(ari_scenicact)
ari_scenicactpd.columns=["ARI"]
ari_scenicactpd["Method"]="Scenic+(eRegulon)"
ari_origpd=pd.DataFrame(ari_orig)
ari_origpd.columns=["ARI"]
ari_origpd["Method"]="Our model(ATAC-seq)"
ari_scenicpd=pd.DataFrame(ari_scenic)
ari_scenicpd.columns=["ARI"]
ari_scenicpd["Method"]="Scenic+(ATAC-seq)"
ari_directpd=pd.DataFrame(ari_direct)
ari_directpd.columns=["ARI"]
ari_directpd["Method"]="DIRECT-NET"
ari_archrpd=pd.DataFrame(ari_archr)
ari_archrpd.columns=["ARI"]
ari_archrpd["Method"]="ArchR"
ari_ciceropd=pd.DataFrame(ari_cicero)
ari_ciceropd.columns=["ARI"]
ari_ciceropd["Method"]="Cicero"
ari_allpd = pd.DataFrame([[ari_all,"All peaks"]], columns=["ARI","Method"])


df_concat = pd.concat([ari_origactpd,ari_scenicactpd,ari_origpd, ari_scenicpd, ari_directpd, ari_archrpd, ari_ciceropd,ari_allpd,ari_allpd,ari_allpd,ari_allpd,ari_allpd], axis=0)

medians = df_concat.groupby(['Method'])['ARI'].median()
vertical_offset = df_concat['ARI'].median() * 0.05 # offset from median for display

x = "Method"
y = "ARI"
order = ["Our model(Activity)","Scenic+(eRegulon)","Our model(ATAC-seq)", "Scenic+(ATAC-seq)", "DIRECT-NET", "ArchR", "Cicero","All peaks"]
my_pal = {"Our model(Activity)": "blue", "Scenic+(eRegulon)": "green","Our model(ATAC-seq)": "blue", "Scenic+(ATAC-seq)": "green", "DIRECT-NET":"yellow", "ArchR": "cyan", "Cicero":"brown", "All peaks":"gray"}


plt.clf()
plt.figure(figsize=(6.5,7))
plt.rcParams.update({'font.size': 18})
ax = sns.boxplot(data=df_concat, x=x, y=y, order=order,palette=my_pal,showfliers= False)
#for xtick in ax.get_xticks():
#    ax.text(xtick,medians[xtick] + vertical_offset,medians[xtick], 
#            horizontalalignment='center',size='x-small',color='w',weight='semibold')

pairs = [('Our model(Activity)', 'Scenic+(eRegulon)'),
         ('Our model(ATAC-seq)', 'Scenic+(ATAC-seq)'),
         ('Our model(ATAC-seq)', 'DIRECT-NET'),
         ('Our model(ATAC-seq)', 'ArchR'),
         ('Our model(ATAC-seq)', 'Cicero'),
         ('Our model(ATAC-seq)', 'All peaks')]

annotator = Annotator(ax, pairs, data=df_concat, x=x, y=y, order=order)
annotator.configure(test='t-test_ind', text_format='star', loc='inside')
annotator.apply_and_annotate()
ax.set_xlabel("Method", fontsize=18)
ax.set_ylabel("Adjusted Rand Index", fontsize=18)
plt.title("Clustering",fontsize=24)
plt.xticks(fontsize=18,rotation=90)
plt.yticks(fontsize=18)
plt.subplots_adjust(left=0.20)
plt.subplots_adjust(bottom=0.40)
plt.savefig("barplot_ARI.pdf",format="pdf")

