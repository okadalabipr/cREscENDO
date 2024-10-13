import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt
import scipy.stats as stats
import sys

args = sys.argv
samplename=str(args[1])


peaklist=pd.read_csv(samplename+"/peaks.bed",sep="\t",header=None)
peaklist=peaklist[[1,2]].to_numpy()
pairlist=pd.read_csv(samplename+"/pair_300000.csv",header=None)
pairlist=pairlist[[1,2,3]].to_numpy()

pairlist_prom=pd.read_csv(samplename+"/pair_promoter.csv",header=None)
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

peaklist_full=pd.read_csv(samplename+"/peaks_extend.bed",sep="\t",header=None)
peaklist_full=peaklist_full.to_numpy()


enha=np.load(samplename+"/enha.npy")
enha_ans=enha.copy()
enha_ans[:,0]=-1
pchic_matrix=np.load(samplename+"/pchic_matrix.npy")
c_array = np.percentile(pchic_matrix[enha_ans!=(-1)], q=[90])
pchic_matrix_th=(pchic_matrix>c_array).astype(float)
pchic_matrix_th[enha_ans==(-1)]=-1

#####################

grad_orig=np.zeros((pairlist.shape[0],max_len,5))
for i in range(1,6):
    fname=samplename+"/try"+str(i)+"/allgrad_ssep_max.npy"
    grad_tmp=np.load(fname)
    grad_orig[:,:,i-1]=grad_tmp

np.save(samplename+"/grad_orig.npy",grad_orig)



auroc_orig=np.zeros((5))
for i in range(5):
    enha_ans=enha.copy()
    enha_ans[:,0]=-1
    grad_tmp=grad_orig[:,:,i]
    gradmm_s=grad_tmp[enha_ans!=-1]
    enha_s=enha_ans[enha_ans!=-1]
    auroc_orig[i]=roc_auc_score(enha_s, gradmm_s)

auroc_orig_hic=np.zeros((5))
for i in range(5):
    enha_ans=enha.copy()
    enha_ans[:,0]=-1
    grad_tmp=grad_orig[:,:,i]
    gradmm_s=grad_tmp[enha_ans!=-1]
    enha_s=pchic_matrix_th[enha_ans!=-1]
    auroc_orig_hic[i]=roc_auc_score(enha_s, gradmm_s)


##########

grad_archr=np.zeros((pairlist.shape[0],max_len,5))
for i in range(2,7):
    fname=samplename+"/archr_for_paper/try"+str(i)+"/archr_matrix.npy"
    grad_tmp=np.load(fname)
    grad_archr[:,:,i-2]=grad_tmp

np.save(samplename+"/grad_archr.npy",grad_archr)

auroc_archr=np.zeros((5))
for i in range(5):
    enha_ans=enha.copy()
    enha_ans[:,0]=-1
    grad_tmp=grad_archr[:,:,i]
    gradmm_s=grad_tmp[(enha_ans!=-1)&(np.abs(posmtx)<=250000)]
    enha_s=enha_ans[(enha_ans!=-1)&(np.abs(posmtx)<=250000)]
    auroc_archr[i]=roc_auc_score(enha_s, gradmm_s)

auroc_archr_hic=np.zeros((5))
for i in range(5):
    enha_ans=enha.copy()
    enha_ans[:,0]=-1
    grad_tmp=grad_archr[:,:,i]
    gradmm_s=grad_tmp[(enha_ans!=-1)&(np.abs(posmtx)<=250000)]
    enha_s=pchic_matrix_th[(enha_ans!=-1)&(np.abs(posmtx)<=250000)]
    auroc_archr_hic[i]=roc_auc_score(enha_s, gradmm_s)

##########
grad_cicero=np.zeros((pairlist.shape[0],max_len,5))
for i in range(1,6):
    fname=samplename+"/cicero_for_paper/try"+str(i)+"/cicero_matrix.npy"
    grad_tmp=np.load(fname)
    grad_cicero[:,:,i-1]=grad_tmp

np.save(samplename+"/grad_cicero.npy",grad_cicero)

auroc_cicero=np.zeros((5))
for i in range(5):
    enha_ans=enha.copy()
    enha_ans[:,0]=-1
    grad_tmp=grad_cicero[:,:,i]
    gradmm_s=grad_tmp[enha_ans!=-1]
    enha_s=enha_ans[enha_ans!=-1]
    auroc_cicero[i]=roc_auc_score(enha_s, gradmm_s)

auroc_cicero_hic=np.zeros((5))
for i in range(5):
    enha_ans=enha.copy()
    enha_ans[:,0]=-1
    grad_tmp=grad_cicero[:,:,i]
    gradmm_s=grad_tmp[enha_ans!=-1]
    enha_s=pchic_matrix_th[enha_ans!=-1]
    auroc_cicero_hic[i]=roc_auc_score(enha_s, gradmm_s)

##########
grad_direct=np.zeros((pairlist.shape[0],max_len,5))
for i in range(1,6):
    fname=samplename+"/direct_for_paper/try"+str(i)+"/"+"peakenha_direct.txt"
    grad_tmp=np.loadtxt(fname)
    grad_direct[:,:,i-1]=grad_tmp

np.save(samplename+"/grad_direct.npy",grad_direct)

auroc_direct=np.zeros((5))
for i in range(5):
    enha_ans=enha.copy()
    enha_ans[:,0]=-1
    grad_tmp=grad_direct[:,:,i]
    gradmm_s=grad_tmp[(enha_ans!=-1)&(np.abs(posmtx)<=250000)]
    enha_s=enha_ans[(enha_ans!=-1)&(np.abs(posmtx)<=250000)]
    auroc_direct[i]=roc_auc_score(enha_s, gradmm_s)

auroc_direct_hic=np.zeros((5))
for i in range(5):
    enha_ans=enha.copy()
    enha_ans[:,0]=-1
    grad_tmp=grad_direct[:,:,i]
    gradmm_s=grad_tmp[(enha_ans!=-1)&(np.abs(posmtx)<=250000)]
    enha_s=pchic_matrix_th[(enha_ans!=-1)&(np.abs(posmtx)<=250000)]
    auroc_direct_hic[i]=roc_auc_score(enha_s, gradmm_s)
##########

grad_scenic=np.zeros((pairlist.shape[0],max_len,5))

for i in range(1,6):
    fname=samplename+"/scenicplus_for_paper/try"+str(i)+"/R2G_absrho_matrix.npy"
    grad_tmp=np.load(fname)
    grad_scenic[:,:,i-1]=grad_tmp

np.save(samplename+"/grad_scenic.npy",grad_scenic)

auroc_scenic=np.zeros((5))
for i in range(5):
    enha_ans=enha.copy()
    enha_ans[:,0]=-1
    grad_tmp=grad_scenic[:,:,i]
    gradmm_s=grad_tmp[(enha_ans!=-1)&(np.abs(posmtx)<=150000)]
    enha_s=enha_ans[(enha_ans!=-1)&(np.abs(posmtx)<=150000)]
    auroc_scenic[i]=roc_auc_score(enha_s, gradmm_s)

auroc_scenic_hic=np.zeros((5))
for i in range(5):
    enha_ans=enha.copy()
    enha_ans[:,0]=-1
    grad_tmp=grad_scenic[:,:,i]
    gradmm_s=grad_tmp[(enha_ans!=-1)&(np.abs(posmtx)<=150000)]
    enha_s=pchic_matrix_th[(enha_ans!=-1)&(np.abs(posmtx)<=150000)]
    auroc_scenic_hic[i]=roc_auc_score(enha_s, gradmm_s)

###########

import seaborn as sns
from statannotations.Annotator import Annotator

auroc_origpd=pd.DataFrame(auroc_orig)
auroc_origpd.columns=["AUROC"]
auroc_origpd["Method"]="Our model"
auroc_scenicpd=pd.DataFrame(auroc_scenic)
auroc_scenicpd.columns=["AUROC"]
auroc_scenicpd["Method"]="Scenic+"
auroc_directpd=pd.DataFrame(auroc_direct)
auroc_directpd.columns=["AUROC"]
auroc_directpd["Method"]="DIRECT-NET"
auroc_archrpd=pd.DataFrame(auroc_archr)
auroc_archrpd.columns=["AUROC"]
auroc_archrpd["Method"]="ArchR"
auroc_ciceropd=pd.DataFrame(auroc_cicero)
auroc_ciceropd.columns=["AUROC"]
auroc_ciceropd["Method"]="Cicero"

df_concat = pd.concat([auroc_origpd,auroc_scenicpd,auroc_directpd, auroc_archrpd, auroc_ciceropd], axis=0)

x = "Method"
y = "AUROC"
order = ["Our model","Scenic+","DIRECT-NET", "ArchR", "Cicero"]
my_pal = {"Our model": "blue", "Scenic+": "green","DIRECT-NET":"yellow", "ArchR": "cyan", "Cicero":"brown", "All peaks":"gray"}

plt.clf()
plt.figure(figsize=(6,8.5))
plt.rcParams.update({"font.size": 18})
ax = sns.boxplot(data=df_concat, x=x, y=y, order=order,palette=my_pal,showfliers= False)
#for xtick in ax.get_xticks():
#    ax.text(xtick,medians[xtick] + vertical_offset,medians[xtick], 
#            horizontalalignment="center",size="x-small",color="w",weight="semibold")

pairs = [("Our model", "Scenic+"),
         ("Our model", "DIRECT-NET"),
         ("Our model", "ArchR"),
         ("Our model", "Cicero")]

annotator = Annotator(ax, pairs, data=df_concat, x=x, y=y, order=order)
annotator.configure(test="t-test_ind", text_format="star", loc="inside")
annotator.apply_and_annotate()
ax.set_xlabel("Method", fontsize=18)
ax.set_ylabel("AUROC", fontsize=18)
plt.title("FANTOM5 enhancers",fontsize=24)
plt.xticks(fontsize=18,rotation=90)
plt.yticks(fontsize=18)
plt.subplots_adjust(left=0.20)
plt.subplots_adjust(bottom=0.40)
plt.savefig(samplename+"/FANTOM5_for_paper.pdf",format="pdf")



#####


auroc_origpd=pd.DataFrame(auroc_orig_hic)
auroc_origpd.columns=["AUROC"]
auroc_origpd["Method"]="Our model"
auroc_scenicpd=pd.DataFrame(auroc_scenic_hic)
auroc_scenicpd.columns=["AUROC"]
auroc_scenicpd["Method"]="Scenic+"
auroc_directpd=pd.DataFrame(auroc_direct_hic)
auroc_directpd.columns=["AUROC"]
auroc_directpd["Method"]="DIRECT-NET"
auroc_archrpd=pd.DataFrame(auroc_archr_hic)
auroc_archrpd.columns=["AUROC"]
auroc_archrpd["Method"]="ArchR"
auroc_ciceropd=pd.DataFrame(auroc_cicero_hic)
auroc_ciceropd.columns=["AUROC"]
auroc_ciceropd["Method"]="Cicero"


df_concat = pd.concat([auroc_origpd,auroc_scenicpd,auroc_directpd, auroc_archrpd, auroc_ciceropd], axis=0)


x = "Method"
y = "AUROC"
order = ["Our model","Scenic+","DIRECT-NET", "ArchR", "Cicero"]
my_pal = {"Our model": "blue", "Scenic+": "green","DIRECT-NET":"yellow", "ArchR": "cyan", "Cicero":"brown", "All peaks":"gray"}

plt.clf()
plt.figure(figsize=(6,8.5))
plt.rcParams.update({"font.size": 18})
ax = sns.boxplot(data=df_concat, x=x, y=y, order=order,palette=my_pal,showfliers= False)


pairs = [("Our model", "Scenic+"),
         ("Our model", "DIRECT-NET"),
         ("Our model", "ArchR"),
         ("Our model", "Cicero")]

annotator = Annotator(ax, pairs, data=df_concat, x=x, y=y, order=order)
annotator.configure(test="t-test_ind", text_format="star", loc="inside")
annotator.apply_and_annotate()
ax.set_xlabel("Method", fontsize=18)
ax.set_ylabel("AUROC", fontsize=18)
plt.title("PCHiC",fontsize=24)
plt.xticks(fontsize=18,rotation=90)
plt.yticks(fontsize=18)
plt.subplots_adjust(left=0.20)
plt.subplots_adjust(bottom=0.40)
plt.savefig(samplename+"/PCHiC_for_paper.pdf",format="pdf")



###############

startpos=0
endpos=250000
binnum=20
dist=(endpos-startpos)/(binnum-1)

auroc_dist_orig=np.zeros((binnum,2,5))
for k in range(5):
    for i in range(binnum):
        tag=dist*i+startpos
        tag2=dist*(i+1)+startpos
        auroc_dist_orig[i,0,k]=(tag+tag2)/2
        tmp=grad_orig[((np.abs(posmtx)>tag)&(np.abs(posmtx)<=tag2)),k]
        tmpe=enha[(np.abs(posmtx)>tag)&(np.abs(posmtx)<=tag2)]
        tmp=tmp[tmpe!=-1]
        tmpe=tmpe[tmpe!=-1]
        if len(tmpe)>0:
            auroc_dist_orig[i,1,k]=roc_auc_score(tmpe,tmp)



auroc_dist_orig_hic=np.zeros((binnum,2,5))
for k in range(5):
    for i in range(binnum):
        tag=dist*i+startpos
        tag2=dist*(i+1)+startpos
        auroc_dist_orig_hic[i,0,k]=(tag+tag2)/2
        tmp=grad_orig[((np.abs(posmtx)>tag)&(np.abs(posmtx)<=tag2)),k]
        tmpe=enha[(np.abs(posmtx)>tag)&(np.abs(posmtx)<=tag2)]
        tmphic=pchic_matrix_th[(np.abs(posmtx)>tag)&(np.abs(posmtx)<=tag2)]
        tmp=tmp[tmpe!=-1]
        tmphic=tmphic[tmpe!=-1]
        if len(tmpe)>0:
            auroc_dist_orig_hic[i,1,k]=roc_auc_score(tmphic,tmp)

auroc_dist_orig_mtx=np.zeros((binnum,3))
auroc_dist_orig_mtx[:,0]=auroc_dist_orig[:,0,0]
auroc_dist_orig_mtx[:,1]=np.mean(auroc_dist_orig[:,1,:],axis=1)
auroc_dist_orig_mtx[:,2]=np.std(auroc_dist_orig[:,1,:],axis=1)

auroc_dist_orig_hic_mtx=np.zeros((binnum,3))
auroc_dist_orig_hic_mtx[:,0]=auroc_dist_orig_hic[:,0,0]
auroc_dist_orig_hic_mtx[:,1]=np.mean(auroc_dist_orig_hic[:,1,:],axis=1)
auroc_dist_orig_hic_mtx[:,2]=np.std(auroc_dist_orig_hic[:,1,:],axis=1)


#########
            
auroc_dist_archr=np.zeros((binnum,2,5))
for k in range(5):
    for i in range(binnum):
        tag=dist*i+startpos
        tag2=dist*(i+1)+startpos
        auroc_dist_archr[i,0,k]=(tag+tag2)/2
        tmp=grad_archr[((np.abs(posmtx)>tag)&(np.abs(posmtx)<=tag2)),k]
        tmpe=enha[(np.abs(posmtx)>tag)&(np.abs(posmtx)<=tag2)]
        tmp=tmp[tmpe!=-1]
        tmpe=tmpe[tmpe!=-1]
        if len(tmpe)>0:
            auroc_dist_archr[i,1,k]=roc_auc_score(tmpe,tmp)

auroc_dist_archr_hic=np.zeros((binnum,2,5))
for k in range(5):
    for i in range(binnum):
        tag=dist*i+startpos
        tag2=dist*(i+1)+startpos
        auroc_dist_archr_hic[i,0,k]=(tag+tag2)/2
        tmp=grad_archr[((np.abs(posmtx)>tag)&(np.abs(posmtx)<=tag2)),k]
        tmpe=enha[(np.abs(posmtx)>tag)&(np.abs(posmtx)<=tag2)]
        tmphic=pchic_matrix_th[(np.abs(posmtx)>tag)&(np.abs(posmtx)<=tag2)]
        tmp=tmp[tmpe!=-1]
        tmphic=tmphic[tmpe!=-1]
        if len(tmpe)>0:
            auroc_dist_archr_hic[i,1,k]=roc_auc_score(tmphic,tmp)

auroc_dist_archr_mtx=np.zeros((binnum,3))
auroc_dist_archr_mtx[:,0]=auroc_dist_archr[:,0,0]
auroc_dist_archr_mtx[:,1]=np.mean(auroc_dist_archr[:,1,:],axis=1)
auroc_dist_archr_mtx[:,2]=np.std(auroc_dist_archr[:,1,:],axis=1)

auroc_dist_archr_hic_mtx=np.zeros((binnum,3))
auroc_dist_archr_hic_mtx[:,0]=auroc_dist_archr_hic[:,0,0]
auroc_dist_archr_hic_mtx[:,1]=np.mean(auroc_dist_archr_hic[:,1,:],axis=1)
auroc_dist_archr_hic_mtx[:,2]=np.std(auroc_dist_archr_hic[:,1,:],axis=1)



#########
            
auroc_dist_cicero=np.zeros((binnum,2,5))
for k in range(5):
    for i in range(binnum):
        tag=dist*i+startpos
        tag2=dist*(i+1)+startpos
        auroc_dist_cicero[i,0,k]=(tag+tag2)/2
        tmp=grad_cicero[((np.abs(posmtx)>tag)&(np.abs(posmtx)<=tag2)),k]
        tmpe=enha[(np.abs(posmtx)>tag)&(np.abs(posmtx)<=tag2)]
        tmp=tmp[tmpe!=-1]
        tmpe=tmpe[tmpe!=-1]
        if len(tmpe)>0:
            auroc_dist_cicero[i,1,k]=roc_auc_score(tmpe,tmp)

auroc_dist_cicero_hic=np.zeros((binnum,2,5))
for k in range(5):
    for i in range(binnum):
        tag=dist*i+startpos
        tag2=dist*(i+1)+startpos
        auroc_dist_cicero_hic[i,0,k]=(tag+tag2)/2
        tmp=grad_cicero[((np.abs(posmtx)>tag)&(np.abs(posmtx)<=tag2)),k]
        tmpe=enha[(np.abs(posmtx)>tag)&(np.abs(posmtx)<=tag2)]
        tmphic=pchic_matrix_th[(np.abs(posmtx)>tag)&(np.abs(posmtx)<=tag2)]
        tmp=tmp[tmpe!=-1]
        tmphic=tmphic[tmpe!=-1]
        if len(tmpe)>0:
            auroc_dist_cicero_hic[i,1,k]=roc_auc_score(tmphic,tmp)

auroc_dist_cicero_mtx=np.zeros((binnum,3))
auroc_dist_cicero_mtx[:,0]=auroc_dist_cicero[:,0,0]
auroc_dist_cicero_mtx[:,1]=np.mean(auroc_dist_cicero[:,1,:],axis=1)
auroc_dist_cicero_mtx[:,2]=np.std(auroc_dist_cicero[:,1,:],axis=1)

auroc_dist_cicero_hic_mtx=np.zeros((binnum,3))
auroc_dist_cicero_hic_mtx[:,0]=auroc_dist_cicero_hic[:,0,0]
auroc_dist_cicero_hic_mtx[:,1]=np.mean(auroc_dist_cicero_hic[:,1,:],axis=1)
auroc_dist_cicero_hic_mtx[:,2]=np.std(auroc_dist_cicero_hic[:,1,:],axis=1)



#########
            
auroc_dist_direct=np.zeros((binnum,2,5))
for k in range(5):
    for i in range(binnum):
        tag=dist*i+startpos
        tag2=dist*(i+1)+startpos
        auroc_dist_direct[i,0,k]=(tag+tag2)/2
        tmp=grad_direct[((np.abs(posmtx)>tag)&(np.abs(posmtx)<=tag2)),k]
        tmpe=enha[(np.abs(posmtx)>tag)&(np.abs(posmtx)<=tag2)]
        tmp=tmp[tmpe!=-1]
        tmpe=tmpe[tmpe!=-1]
        if len(tmpe)>0:
            auroc_dist_direct[i,1,k]=roc_auc_score(tmpe,tmp)

auroc_dist_direct_hic=np.zeros((binnum,2,5))
for k in range(5):
    for i in range(binnum):
        tag=dist*i+startpos
        tag2=dist*(i+1)+startpos
        auroc_dist_direct_hic[i,0,k]=(tag+tag2)/2
        tmp=grad_direct[((np.abs(posmtx)>tag)&(np.abs(posmtx)<=tag2)),k]
        tmpe=enha[(np.abs(posmtx)>tag)&(np.abs(posmtx)<=tag2)]
        tmphic=pchic_matrix_th[(np.abs(posmtx)>tag)&(np.abs(posmtx)<=tag2)]
        tmp=tmp[tmpe!=-1]
        tmphic=tmphic[tmpe!=-1]
        if len(tmpe)>0:
            auroc_dist_direct_hic[i,1,k]=roc_auc_score(tmphic,tmp)


auroc_dist_direct_mtx=np.zeros((binnum,3))
auroc_dist_direct_mtx[:,0]=auroc_dist_direct[:,0,0]
auroc_dist_direct_mtx[:,1]=np.mean(auroc_dist_direct[:,1,:],axis=1)
auroc_dist_direct_mtx[:,2]=np.std(auroc_dist_direct[:,1,:],axis=1)

auroc_dist_direct_hic_mtx=np.zeros((binnum,3))
auroc_dist_direct_hic_mtx[:,0]=auroc_dist_direct_hic[:,0,0]
auroc_dist_direct_hic_mtx[:,1]=np.mean(auroc_dist_direct_hic[:,1,:],axis=1)
auroc_dist_direct_hic_mtx[:,2]=np.std(auroc_dist_direct_hic[:,1,:],axis=1)


#########
            
auroc_dist_scenic=np.zeros((binnum,2,5))
for k in range(5):
    for i in range(binnum):
        tag=dist*i+startpos
        tag2=dist*(i+1)+startpos
        auroc_dist_scenic[i,0,k]=(tag+tag2)/2
        tmp=grad_scenic[((np.abs(posmtx)>tag)&(np.abs(posmtx)<=tag2)),k]
        tmpe=enha[(np.abs(posmtx)>tag)&(np.abs(posmtx)<=tag2)]
        tmp=tmp[tmpe!=-1]
        tmpe=tmpe[tmpe!=-1]
        if len(tmpe)>0:
            auroc_dist_scenic[i,1,k]=roc_auc_score(tmpe,tmp)

auroc_dist_scenic_hic=np.zeros((binnum,2,5))
for k in range(5):
    for i in range(binnum):
        tag=dist*i+startpos
        tag2=dist*(i+1)+startpos
        auroc_dist_scenic_hic[i,0,k]=(tag+tag2)/2
        tmp=grad_scenic[((np.abs(posmtx)>tag)&(np.abs(posmtx)<=tag2)),k]
        tmpe=enha[(np.abs(posmtx)>tag)&(np.abs(posmtx)<=tag2)]
        tmphic=pchic_matrix_th[(np.abs(posmtx)>tag)&(np.abs(posmtx)<=tag2)]
        tmp=tmp[tmpe!=-1]
        tmphic=tmphic[tmpe!=-1]
        if len(tmpe)>0:
            auroc_dist_scenic_hic[i,1,k]=roc_auc_score(tmphic,tmp)

auroc_dist_scenic_mtx=np.zeros((binnum,3))
auroc_dist_scenic_mtx[:,0]=auroc_dist_scenic[:,0,0]
auroc_dist_scenic_mtx[:,1]=np.mean(auroc_dist_scenic[:,1,:],axis=1)
auroc_dist_scenic_mtx[:,2]=np.std(auroc_dist_scenic[:,1,:],axis=1)

auroc_dist_scenic_hic_mtx=np.zeros((binnum,3))
auroc_dist_scenic_hic_mtx[:,0]=auroc_dist_scenic_hic[:,0,0]
auroc_dist_scenic_hic_mtx[:,1]=np.mean(auroc_dist_scenic_hic[:,1,:],axis=1)
auroc_dist_scenic_hic_mtx[:,2]=np.std(auroc_dist_scenic_hic[:,1,:],axis=1)

tmp=(auroc_dist_scenic_mtx[:,0]-150000)
tmp[tmp>0]=-150000
sceidx=np.argmax(tmp)

fig, ax = plt.subplots(figsize=(7,7))
ax.errorbar(auroc_dist_orig_mtx[:,0], auroc_dist_orig_mtx[:,1], auroc_dist_orig_mtx[:,2],color="b",label="Our model")
ax.errorbar(auroc_dist_scenic_mtx[:sceidx,0], auroc_dist_scenic_mtx[:sceidx,1], auroc_dist_scenic_mtx[:sceidx,2],color="g",label="Scenic+")
ax.errorbar(auroc_dist_direct_mtx[:,0], auroc_dist_direct_mtx[:,1], auroc_dist_direct_mtx[:,2],color="y",label="DIRECT-NET")
ax.errorbar(auroc_dist_archr_mtx[:,0], auroc_dist_archr_mtx[:,1], auroc_dist_archr_mtx[:,2],color="c",label="Archr")
ax.errorbar(auroc_dist_cicero_mtx[:,0], auroc_dist_cicero_mtx[:,1], auroc_dist_cicero_mtx[:,2],color="k",label="Cicero")
plt.tick_params(labelsize=16,width=2)
ax.spines["top"].set_linewidth(2)
ax.spines["left"].set_linewidth(2)
ax.spines["bottom"].set_linewidth(2)
ax.spines["right"].set_linewidth(2)
plt.title("FANTOM5 enhancer", fontsize=28)
ax.set_xlim(startpos,endpos)
ax.set_xlabel("Absolute Distance from TSS(bp)", fontsize=20)
ax.set_ylabel("AUROC", fontsize=20)
plt.subplots_adjust(left=0.18)
plt.subplots_adjust(bottom=0.25)
ax.legend(fontsize=18,loc="center left", bbox_to_anchor=(1., .5))
plt.savefig(samplename+"/AUROC_FANTOM5_dist.pdf", bbox_inches="tight")
plt.show()

fig, ax = plt.subplots(figsize=(7,7))
ax.errorbar(auroc_dist_orig_hic_mtx[:,0], auroc_dist_orig_hic_mtx[:,1], auroc_dist_orig_hic_mtx[:,2],color="b",label="Our model")
ax.errorbar(auroc_dist_scenic_hic_mtx[:sceidx,0], auroc_dist_scenic_hic_mtx[:sceidx,1], auroc_dist_scenic_hic_mtx[:sceidx,2],color="g",label="Scenic+")
ax.errorbar(auroc_dist_direct_hic_mtx[:,0], auroc_dist_direct_hic_mtx[:,1], auroc_dist_direct_hic_mtx[:,2],color="y",label="DIRECT-NET")
ax.errorbar(auroc_dist_archr_hic_mtx[:,0], auroc_dist_archr_hic_mtx[:,1], auroc_dist_archr_hic_mtx[:,2],color="c",label="Archr")
ax.errorbar(auroc_dist_cicero_hic_mtx[:,0], auroc_dist_cicero_hic_mtx[:,1], auroc_dist_cicero_hic_mtx[:,2],color="k",label="Cicero")
plt.tick_params(labelsize=16,width=2)
ax.spines["top"].set_linewidth(2)
ax.spines["left"].set_linewidth(2)
ax.spines["bottom"].set_linewidth(2)
ax.spines["right"].set_linewidth(2)
plt.title("PCHiC", fontsize=28)
ax.set_xlim(startpos,endpos)
ax.set_xlabel("Absolute Distance from TSS(bp)", fontsize=20)
ax.set_ylabel("AUROC", fontsize=20)
plt.subplots_adjust(left=0.18)
plt.subplots_adjust(bottom=0.25)
ax.legend(fontsize=18,loc="center left", bbox_to_anchor=(1., .5))
plt.savefig(samplename+"/AUROC_PCHiC_dist.pdf", bbox_inches="tight")
plt.show()




#####################
##########

eqtl=np.load(samplename+"/eqtl_mtx_all_cav.npy")
eqtl_sig=np.load(samplename+"/eqtl_mtx_sig05.npy")

enha_ans=enha.copy()
enha_ans[:,0]=-1

grad_orig_th=np.zeros(grad_orig.shape)
for i in range(5):
    c_array = np.percentile(grad_orig[enha!=(-1),i], q=[90])
    grad_orig_th[:,:,i]=(grad_orig[:,:,i]>c_array).astype(float)

grad_direct_th=np.zeros(grad_direct.shape)
for i in range(5):
    c_array = np.percentile(grad_direct[enha!=(-1),i], q=[90])
    grad_direct_th[:,:,i]=(grad_direct[:,:,i]>c_array).astype(float)

grad_archr_th=np.zeros(grad_archr.shape)
for i in range(5):
    c_array = np.percentile(grad_archr[enha!=(-1),i], q=[90])
    grad_archr_th[:,:,i]=(grad_archr[:,:,i]>c_array).astype(float)

grad_cicero_th=np.zeros(grad_cicero.shape)
for i in range(5):
    c_array = np.percentile(grad_cicero[enha!=(-1),i], q=[90])
    grad_cicero_th[:,:,i]=(grad_cicero[:,:,i]>c_array).astype(float)

grad_scenic_th=np.zeros(grad_scenic.shape)
for i in range(5):
    c_array = np.percentile(grad_scenic[enha!=(-1),i], q=[90])
    grad_scenic_th[:,:,i]=(grad_scenic[:,:,i]>c_array).astype(float)




eqtl_s=eqtl[enha_ans!=(-1)]
eqtl_sig_s=eqtl_sig[enha_ans!=(-1)]
enha_s=enha[enha_ans!=(-1)]
eqtl_ratio_fantom=(eqtl_sig_s[enha_s==1].mean()/eqtl_s[enha_s==1].mean())/(eqtl_sig_s.mean()/eqtl_s.mean())

enha_s=pchic_matrix_th[enha_ans!=(-1)]
eqtl_ratio_pchic=(eqtl_sig_s[enha_s==1].mean()/eqtl_s[enha_s==1].mean())/(eqtl_sig_s.mean()/eqtl_s.mean())

eqtl_ratio_orig=np.zeros(5)
for i in range(5):
    enha_s=grad_orig_th[enha_ans!=(-1),i]
    eqtl_ratio_orig[i]=(eqtl_sig_s[enha_s==1].mean()/eqtl_s[enha_s==1].mean())/(eqtl_sig_s.mean()/eqtl_s.mean())

eqtl_ratio_archr=np.zeros(5)
for i in range(5):
    enha_s=grad_archr_th[enha_ans!=(-1),i]
    eqtl_ratio_archr[i]=(eqtl_sig[(grad_archr_th[:,:,i]==1)&(np.abs(posmtx)<=250000)&(enha_ans!=(-1))].mean()/eqtl[(grad_archr_th[:,:,i]==1)&(np.abs(posmtx)<=250000)&(enha_ans!=(-1))].mean())/(eqtl_sig[(np.abs(posmtx)<=250000)&(enha_ans!=(-1))].mean()/eqtl[(np.abs(posmtx)<=250000)&(enha_ans!=(-1))].mean())

eqtl_ratio_direct=np.zeros(5)
for i in range(5):
    enha_s=grad_direct_th[enha_ans!=(-1),i]
    eqtl_ratio_direct[i]=(eqtl_sig[(grad_direct_th[:,:,i]==1)&(np.abs(posmtx)<=250000)&(enha_ans!=(-1))].mean()/eqtl[(grad_direct_th[:,:,i]==1)&(np.abs(posmtx)<=250000)&(enha_ans!=(-1))].mean())/(eqtl_sig[(np.abs(posmtx)<=250000)&(enha_ans!=(-1))].mean()/eqtl[(np.abs(posmtx)<=250000)&(enha_ans!=(-1))].mean())

eqtl_ratio_cicero=np.zeros(5)
for i in range(5):
    enha_s=grad_cicero_th[enha_ans!=(-1),i]
    eqtl_ratio_cicero[i]=(eqtl_sig_s[enha_s==1].mean()/eqtl_s[enha_s==1].mean())/(eqtl_sig_s.mean()/eqtl_s.mean())

eqtl_ratio_scenic=np.zeros(5)
for i in range(5):
    enha_s=grad_scenic_th[enha_ans!=(-1),i]
    eqtl_ratio_scenic[i]=(eqtl_sig[(grad_scenic_th[:,:,i]==1)&(np.abs(posmtx)<=150000)&(enha_ans!=(-1))].mean()/eqtl[(grad_scenic_th[:,:,i]==1)&(np.abs(posmtx)<=150000)&(enha_ans!=(-1))].mean())/(eqtl_sig[(np.abs(posmtx)<=150000)&(enha_ans!=(-1))].mean()/eqtl[(np.abs(posmtx)<=150000)&(enha_ans!=(-1))].mean())


#####


points = (eqtl_ratio_fantom,eqtl_ratio_pchic,eqtl_ratio_orig, eqtl_ratio_scenic,eqtl_ratio_direct,eqtl_ratio_archr,eqtl_ratio_cicero)




auroc_fantompd = pd.DataFrame([[eqtl_ratio_fantom,"FANTOM5"]], columns=["eqtl","Method"])
auroc_pchicpd = pd.DataFrame([[eqtl_ratio_pchic,"PCHiC"]], columns=["eqtl","Method"])

auroc_origpd=pd.DataFrame(eqtl_ratio_orig)
auroc_origpd.columns=["eqtl"]
auroc_origpd["Method"]="Our model"
auroc_scenicpd=pd.DataFrame(eqtl_ratio_scenic)
auroc_scenicpd.columns=["eqtl"]
auroc_scenicpd["Method"]="Scenic+"
auroc_directpd=pd.DataFrame(eqtl_ratio_direct)
auroc_directpd.columns=["eqtl"]
auroc_directpd["Method"]="DIRECT-NET"
auroc_archrpd=pd.DataFrame(eqtl_ratio_archr)
auroc_archrpd.columns=["eqtl"]
auroc_archrpd["Method"]="ArchR"
auroc_ciceropd=pd.DataFrame(eqtl_ratio_cicero)
auroc_ciceropd.columns=["eqtl"]
auroc_ciceropd["Method"]="Cicero"


df_concat = pd.concat([auroc_fantompd,auroc_pchicpd,auroc_origpd,auroc_scenicpd,auroc_directpd, auroc_archrpd, auroc_ciceropd], axis=0)


x = "Method"
y = "eqtl"
order = ["FANTOM5","PCHiC","Our model","Scenic+","DIRECT-NET", "ArchR", "Cicero"]
my_pal = {"FANTOM5": "black","PCHiC": "black","Our model": "blue", "Scenic+": "green","DIRECT-NET":"yellow", "ArchR": "cyan", "Cicero":"brown", "All peaks":"gray"}

plt.clf()
plt.figure(figsize=(6,8.5))
plt.rcParams.update({"font.size": 18})
ax = sns.boxplot(data=df_concat, x=x, y=y, order=order,palette=my_pal,showfliers= False)
#for xtick in ax.get_xticks():
#    ax.text(xtick,medians[xtick] + vertical_offset,medians[xtick], 
#            horizontalalignment="center",size="x-small",color="w",weight="semibold")

pairs = [("Our model", "Scenic+"),
         ("Our model", "DIRECT-NET"),
         ("Our model", "ArchR"),
         ("Our model", "Cicero")]

annotator = Annotator(ax, pairs, data=df_concat, x=x, y=y, order=order)
annotator.configure(test="t-test_ind", text_format="star", loc="inside")
annotator.apply_and_annotate()
ax.set_xlabel("Method", fontsize=18)
ax.set_ylabel("Enrichment ratio of eQTL", fontsize=18)
plt.title("eQTL enrichment",fontsize=24)
plt.xticks(fontsize=18,rotation=90)
plt.yticks(fontsize=18)
plt.subplots_adjust(left=0.20)
plt.subplots_adjust(bottom=0.40)
plt.savefig(samplename+"/eqtl_for_paper.pdf",format="pdf")

