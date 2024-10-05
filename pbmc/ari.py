import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
import scipy.stats as stats
import seaborn as sns
from statannotations.Annotator import Annotator
import sys

args = sys.argv
samplename=str(args[1])

fname=samplename+"/Deeplift_full_ver2_all.npy"
enhamtx=np.load(fname)

genelist=pd.read_csv(samplename+"/pair_promoter.csv",sep=",",header=None)
grad=np.load(samplename+"/allgrad_ssep_max.npy")


pairlist_prom=pd.read_csv(samplename+"/pair_promoter.csv",header=None)
pairlist_prom=pairlist_prom[[1,2,3]].to_numpy()

pairlist=pd.read_csv(samplename+"/pair_300000.csv",header=None)
pairlist=pairlist[[1,2,3]].to_numpy()


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
   cellclus=pd.read_csv(samplename+"/cellcluster_all.txt",sep=",",header=None)
   cellid=pd.read_csv(samplename+"/celllist_scVI.csv",sep="\t",header=None)
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

np.save(samplename+"/ari_origact.npy",ari_origact)
np.save(samplename+"/ari_scenicact.npy",ari_scenicact)



######

cellclus=pd.read_csv(samplename+"/cellcluster_all.txt",sep=",",header=None)
cellclus=cellclus.to_numpy()
cellclus=cellclus[:,0]

grad_archr=np.load(samplename+"/grad_archr.npy")
grad_cicero=np.load(samplename+"/grad_cicero.npy")
grad_direct=np.load(samplename+"/grad_direct.npy")
grad_scenic=np.load(samplename+"/grad_scenic.npy")
grad_orig=np.load(samplename+"/grad_orig.npy")


ATAC_use=np.load(samplename+"/atac_count.npy")

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


enhatag=np.load(samplename+"/enha.npy")
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


np.save(samplename+"/ATACcluster_orig.npy",ATACcluster_orig)
np.save(samplename+"/ATACcluster_archr.npy",ATACcluster_archr)
np.save(samplename+"/ATACcluster_cicero.npy",ATACcluster_cicero)
np.save(samplename+"/ATACcluster_scenic.npy",ATACcluster_scenic)
np.save(samplename+"/ATACcluster_direct.npy",ATACcluster_direct)


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

