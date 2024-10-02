import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import f1_score
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import scipy.stats as stats

peaklist=pd.read_csv('pbmc_fold1/peaks.bed',sep="\t",header=None)
peaklist=peaklist[[1,2]].to_numpy()
pairlist=pd.read_csv('pbmc_fold1/pair_300000.csv',header=None)
pairlist=pairlist[[1,2,3]].to_numpy()

pairlist_prom=pd.read_csv('pbmc_fold1/pair_promoter.csv',header=None)
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

peaklist_full=pd.read_csv('pbmc_fold1/peaks_extend.bed',sep="\t",header=None)
peaklist_full=peaklist_full.to_numpy()

hic=pd.read_csv('PCHiC_peak_matrix_cutoff5.tsv',sep="\t")
hic_s=hic.iloc[:,[0,1,2,5,6,7]]
hic_s["signal"]=hic.iloc[:,11:28].to_numpy().max(axis=1)
hic_s["baitChr"]="chr"+hic_s["baitChr"].astype(str)
hic_s["oeChr"]="chr"+hic_s["oeChr"].astype(str)

hic_s.iloc[:,0:3].to_csv("bait.bed",header=False, index=False)
hic_s.iloc[:,3:6].to_csv("oe.bed",header=False, index=False)

bait=hic_s.iloc[:,0:3]
oe=hic_s.iloc[:,3:6]

bait["tag"]=bait["baitChr"]+str("-")+bait["baitStart"].astype('str')+str("-")+bait["baitEnd"].astype('str')
bait.to_csv("bait_hg19_tag.bed", sep="\t",header=False, index=False)
oe["tag"]=oe["oeChr"]+str("-")+oe["oeStart"].astype('str')+str("-")+oe["oeEnd"].astype('str')
oe.to_csv("oe_hg19_tag.bed", sep="\t",header=False, index=False)

####
liftOver bait_hg19_tag.bed hg19ToHg38.over.chain bait_hg38_tag.bed bait_unlifted.bed
liftOver oe_hg19_tag.bed hg19ToHg38.over.chain oe_hg38_tag.bed oe_unlifted.bed
####


bait_hg19=pd.read_csv("bait_hg19_tag.bed",sep="\t",header=None,names=("chr","start","end","tag"))
bait_hg38=pd.read_csv("bait_hg38_tag.bed",sep="\t",header=None,names=("chr","start","end","tag"))
chrs = ["chr"+str(i) for i in range(1,23)] + ["chrX", "chrY"]
bait_hg38_s=bait_hg38[bait_hg38["chr"].isin(chrs)]

oe_hg19=pd.read_csv("oe_hg19_tag.bed",sep="\t",header=None,names=("chr","start","end","tag"))
oe_hg38=pd.read_csv("oe_hg38_tag.bed",sep="\t",header=None,names=("chr","start","end","tag"))
chrs = ["chr"+str(i) for i in range(1,23)] + ["chrX", "chrY"]
oe_hg38_s=oe_hg38[oe_hg38["chr"].isin(chrs)]

bait_hg38_sd=bait_hg38_s.drop_duplicates(subset='tag')
bait_hg19_m=pd.merge(bait_hg19, bait_hg38_sd, on='tag', how='left')
oe_hg38_sd=oe_hg38_s.drop_duplicates(subset='tag')
oe_hg19_m=pd.merge(oe_hg19, oe_hg38_sd, on='tag', how='left')

bait_hg19_to_38=bait_hg19_m.iloc[:,4:7]
bait_hg19_to_38.columns=["baitChr","baitStart","baitEnd"]
oe_hg19_to_38=oe_hg19_m.iloc[:,4:7]
oe_hg19_to_38.columns=["oeChr","oeStart","oeEnd"]

hg19_m=pd.concat([bait_hg19_to_38, oe_hg19_to_38],axis=1)
hg19_m["signal"]=hic.iloc[:,11:28].to_numpy().max(axis=1)
hg19_mdn=hg19_m.dropna(how='any')

hic_sn=hg19_mdn.to_numpy()

pchic_matrix=np.zeros((pairlist.shape[0],max_len))
for gn in range(pairlist.shape[0]):
    print(gn)
    peaknum=(posmtx[gn,:]!=(-1)).sum()
    prom_id=peakidmtx[gn,0].astype(int)
    tmphic=hic_sn[hic_sn[:,0]==peaklist_full[prom_id,0],:].copy()
    prom_in=(tmphic[:,2]>=peaklist_full[prom_id,1])&(tmphic[:,1]<=peaklist_full[prom_id,2])
    tmphic_in=tmphic[(prom_in),:]
    for p1 in range(peaknum):
        p1_id=peakidmtx[gn,p1].astype(int)
        p1_in=(tmphic_in[:,5]>=peaklist_full[p1_id,1])&(tmphic_in[:,4]<=peaklist_full[p1_id,2])
        maemuki=tmphic_in[(p1_in),:]
        pchic_matrix[gn,p1]=maemuki.shape[0]

np.save("pbmc_fold1/pchic_matrix.npy",pchic_matrix)




##################

enha_bed=pd.read_csv('F5.hg38.enhancers.bed',sep="\t",header=None)
enha_bed=enha_bed.to_numpy()

enha_matrix=np.zeros((pairlist.shape[0],max_len))
enha_matrix=enha_matrix-1
for gn in range(pairlist.shape[0]):
    print(gn)
    peaknum=(posmtx[gn,:]!=(-1)).sum()
    prom_id=peakidmtx[gn,0].astype(int)
    tmphic=enha_bed[enha_bed[:,0]==peaklist_full[prom_id,0],:].copy()
    for p1 in range(peaknum):
        p1_id=peakidmtx[gn,p1].astype(int)
        p1_in=(tmphic[:,2]>=peaklist_full[p1_id,1])&(tmphic[:,1]<=peaklist_full[p1_id,2])
        maemuki=tmphic[(p1_in),:]
        enha_matrix[gn,p1]=int(maemuki.shape[0]>0)

np.save("pbmc_fold1/enha.npy",enha_matrix)

enha=np.load("pbmc_fold1/enha.npy")
enha_ans=enha.copy()
enha_ans[:,0]=-1
pchic_matrix=np.load("pbmc_fold1/pchic_matrix.npy")
c_array = np.percentile(pchic_matrix[enha_ans!=(-1)], q=[90])
pchic_matrix_th=(pchic_matrix>c_array).astype(float)
pchic_matrix_th[enha_ans==(-1)]=-1


#####################

grad_orig=np.zeros((pairlist.shape[0],max_len,5))
for i in range(1,6):
    fname="pbmcnew/try"+str(i)+"/allgrad_ssep_max.npy"
    grad_tmp=np.load(fname)
    grad_orig[:,:,i-1]=grad_tmp

enhath=np.percentile(grad_orig[enha!=(-1),0], [40])
#enhath=np.percentile(grad_orig[enha!=(-1),0], [10,20,30,35,40,45,50,60,70,85])
#enhath=np.percentile(grad_orig[enha!=(-1),0], [85])

np.save("grad_orig.npy",grad_orig)

(grad_orig[enha!=(-1),0]>enhath).sum()

grad_orig2=grad_orig.copy()
grad_orig2[grad_orig2<enhath]=0

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


auroc_orig2=np.zeros((5))
for i in range(5):
    enha_ans=enha.copy()
    enha_ans[:,0]=-1
    grad_tmp=grad_orig2[:,:,i]
    gradmm_s=grad_tmp[(enha_ans!=-1)&(np.abs(posmtx)<=150000)]
    enha_s=enha_ans[(enha_ans!=-1)&(np.abs(posmtx)<=150000)]
    auroc_orig2[i]=roc_auc_score(enha_s, gradmm_s)

auroc_orig_hic2=np.zeros((5))
for i in range(5):
    enha_ans=enha.copy()
    enha_ans[:,0]=-1
    grad_tmp=grad_orig2[:,:,i]
    gradmm_s=grad_tmp[(enha_ans!=-1)&(np.abs(posmtx)<=150000)]
    enha_s=pchic_matrix_th[(enha_ans!=-1)&(np.abs(posmtx)<=150000)]
    auroc_orig_hic2[i]=roc_auc_score(enha_s, gradmm_s)
##########

grad_archr=np.zeros((pairlist.shape[0],max_len,5))
for i in range(2,7):
    fname="archr_for_paper/try"+str(i)+"/archr_matrix.npy"
    grad_tmp=np.load(fname)
    grad_archr[:,:,i-2]=grad_tmp

enhath=np.percentile(grad_archr[enha!=(-1),0], [40,85])
#enhath=np.percentile(grad_archr[enha!=(-1),0], [10,20,30,35,40,45,50,60,70,85])

np.save("grad_archr.npy",grad_archr)

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
    fname="cicero_for_paper/try"+str(i)+"/cicero_matrix.npy"
    grad_tmp=np.load(fname)
    grad_cicero[:,:,i-1]=grad_tmp

enhath=np.percentile(grad_cicero[enha!=(-1),0], [10,20,30,40,50,60,70,85])

np.save("grad_cicero.npy",grad_cicero)

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
    fname="peakenha_direct_"+str(1)+".txt"
    grad_tmp=np.loadtxt(fname)
    grad_direct[:,:,i-1]=grad_tmp

enhath=np.percentile(grad_direct[enha!=(-1),0], [10,20,30,35,40,45,50,60,70,85])

np.save("grad_direct.npy",grad_direct)

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

fname="scenicplus/scenicplus0/R2G_absrho_matrix.npy"
grad_tmp=np.load(fname)
grad_scenic[:,:,0]=grad_tmp

fname="scenicplus/scenicplus4/R2G_absrho_matrix.npy"
grad_tmp=np.load(fname)
grad_scenic[:,:,1]=grad_tmp

fname="scenicplus/scenicplus7/R2G_absrho_matrix.npy"
grad_tmp=np.load(fname)
grad_scenic[:,:,2]=grad_tmp

fname="scenicplus/scenicplus2/R2G_absrho_matrix.npy"
grad_tmp=np.load(fname)
grad_scenic[:,:,3]=grad_tmp

fname="scenicplus/scenicplus8/R2G_absrho_matrix.npy"
grad_tmp=np.load(fname)
grad_scenic[:,:,4]=grad_tmp


np.save("grad_scenic.npy",grad_scenic)
enhath=np.percentile(grad_scenic[enha!=(-1),1], [10,20,30,35,40,45,50,60,70,85])

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
fantom_t_scenic, fantom_p_scenic = stats.ttest_ind(auroc_orig, auroc_scenic)
fantom_t_direct, fantom_p_direct = stats.ttest_ind(auroc_orig, auroc_direct)
fantom_t_archr, fantom_p_archr = stats.ttest_ind(auroc_orig, auroc_archr)
fantom_t_cicero, fantom_p_cicero = stats.ttest_ind(auroc_orig, auroc_cicero)

points = (auroc_orig, auroc_scenic,auroc_direct,auroc_archr,auroc_cicero)
fig, ax = plt.subplots(figsize=(6, 6.5))
bp = ax.boxplot(points)
ax.set_xticklabels(["Our model", "Scenic+", "DIRECT-NET", "ArchR", "Cicero"], fontsize=16,rotation=45,fontweight='bold')
plt.tick_params(labelsize=16,width=2)
ax.spines["top"].set_linewidth(2)
ax.spines["left"].set_linewidth(2)
ax.spines["bottom"].set_linewidth(2)
ax.spines["right"].set_linewidth(2)
plt.title("FANTOM5 enhancer", fontsize=28)
plt.ylim(0.43,0.68)
ax.set_xlabel("Method", fontsize=20)
ax.set_ylabel("AUROC", fontsize=20)
plt.subplots_adjust(left=0.18)
plt.subplots_adjust(bottom=0.25)
plt.savefig("FANTOM5_for_paper.pdf")
plt.show()


hic_t_scenic, hic_p_scenic = stats.ttest_ind(auroc_orig_hic, auroc_scenic_hic)
hic_t_direct, hic_p_direct = stats.ttest_ind(auroc_orig_hic, auroc_direct_hic)
hic_t_archr, hic_p_archr = stats.ttest_ind(auroc_orig_hic, auroc_archr_hic)
hic_t_cicero, hic_p_cicero = stats.ttest_ind(auroc_orig_hic, auroc_cicero_hic)

points = (auroc_orig_hic, auroc_scenic_hic,auroc_direct_hic,auroc_archr_hic,auroc_cicero_hic)
fig, ax = plt.subplots(figsize=(6, 6.5))
bp = ax.boxplot(points)
ax.set_xticklabels(["Our model", "Scenic+", "DIRECT-NET", "ArchR", "Cicero"], fontsize=16,rotation=45,fontweight='bold')
plt.tick_params(labelsize=16,width=2)
ax.spines["top"].set_linewidth(2)
ax.spines["left"].set_linewidth(2)
ax.spines["bottom"].set_linewidth(2)
ax.spines["right"].set_linewidth(2)
plt.title("PCHiC", fontsize=28)
plt.ylim(0.51,0.595)
ax.set_xlabel("Method", fontsize=20)
ax.set_ylabel("AUROC", fontsize=20)
plt.subplots_adjust(left=0.18)
plt.subplots_adjust(bottom=0.25)
plt.savefig("PCHiC_for_paper.png")
plt.show()

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
plt.rcParams.update({'font.size': 18})
ax = sns.boxplot(data=df_concat, x=x, y=y, order=order,palette=my_pal,showfliers= False)
#for xtick in ax.get_xticks():
#    ax.text(xtick,medians[xtick] + vertical_offset,medians[xtick], 
#            horizontalalignment='center',size='x-small',color='w',weight='semibold')

pairs = [('Our model', 'Scenic+'),
         ('Our model', 'DIRECT-NET'),
         ('Our model', 'ArchR'),
         ('Our model', 'Cicero')]

annotator = Annotator(ax, pairs, data=df_concat, x=x, y=y, order=order)
annotator.configure(test='t-test_ind', text_format='star', loc='inside')
annotator.apply_and_annotate()
ax.set_xlabel("Method", fontsize=18)
ax.set_ylabel("AUROC", fontsize=18)
plt.title("FANTOM5 enhancers",fontsize=24)
plt.xticks(fontsize=18,rotation=90)
plt.yticks(fontsize=18)
plt.subplots_adjust(left=0.20)
plt.subplots_adjust(bottom=0.40)
plt.savefig("FANTOM5_for_paper.pdf",format="pdf")




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
plt.rcParams.update({'font.size': 18})
ax = sns.boxplot(data=df_concat, x=x, y=y, order=order,palette=my_pal,showfliers= False)
#for xtick in ax.get_xticks():
#    ax.text(xtick,medians[xtick] + vertical_offset,medians[xtick], 
#            horizontalalignment='center',size='x-small',color='w',weight='semibold')

pairs = [('Our model', 'Scenic+'),
         ('Our model', 'DIRECT-NET'),
         ('Our model', 'ArchR'),
         ('Our model', 'Cicero')]

annotator = Annotator(ax, pairs, data=df_concat, x=x, y=y, order=order)
annotator.configure(test='t-test_ind', text_format='star', loc='inside')
annotator.apply_and_annotate()
ax.set_xlabel("Method", fontsize=18)
ax.set_ylabel("AUROC", fontsize=18)
plt.title("PCHiC",fontsize=24)
plt.xticks(fontsize=18,rotation=90)
plt.yticks(fontsize=18)
plt.subplots_adjust(left=0.20)
plt.subplots_adjust(bottom=0.40)
plt.savefig("PCHiC_for_paper.pdf",format="pdf")



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
ax.legend(fontsize=18,loc='center left', bbox_to_anchor=(1., .5))
plt.savefig("AUROC_FANTOM5_dist.pdf", bbox_inches='tight')
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
ax.legend(fontsize=18,loc='center left', bbox_to_anchor=(1., .5))
plt.savefig("AUROC_PCHiC_dist.pdf", bbox_inches='tight')
plt.show()




#####################
##########

eqtl=np.load("eQTL_for_paper/eqtl_mtx_all_cav.npy")
eqtl_sig=np.load("eQTL_for_paper/eqtl_mtx_sig05.npy")

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

#eqgene=eqtl_sig.sum(axis=1)
#enha_ans[eqgene==0,:]=-1



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



eqtl_t_scenic, eqtl_p_scenic = stats.ttest_ind(eqtl_ratio_orig, eqtl_ratio_scenic)
eqtl_t_direct, eqtl_p_direct = stats.ttest_ind(eqtl_ratio_orig, eqtl_ratio_direct)
eqtl_t_archr, eqtl_p_archr = stats.ttest_ind(eqtl_ratio_orig, eqtl_ratio_archr)
eqtl_t_cicero, eqtl_p_cicero = stats.ttest_ind(eqtl_ratio_orig, eqtl_ratio_cicero)

from scipy.stats import ttest_1samp
ttest_1samp(eqtl_ratio_orig, popmean=eqtl_ratio_fantom)
ttest_1samp(eqtl_ratio_orig, popmean=eqtl_ratio_pchic)

colors = ["k","k","b","g","y","c","k"]
points = (eqtl_ratio_fantom,eqtl_ratio_pchic,eqtl_ratio_orig, eqtl_ratio_scenic,eqtl_ratio_direct,eqtl_ratio_archr,eqtl_ratio_cicero)
fig, ax = plt.subplots(figsize=(6, 6.5))
bp = ax.boxplot(points)
plt.tick_params(labelsize=16,width=2)
ax.spines["top"].set_linewidth(2)
ax.spines["left"].set_linewidth(2)
ax.spines["bottom"].set_linewidth(2)
ax.spines["right"].set_linewidth(2)
ax.set_xticklabels(["FANTOM5","PCHiC","Our model", "Scenic+", "DIRECT-NET", "ArchR", "Cicero"], fontsize=16,rotation=45,fontweight='bold')
plt.title("eQTL enrichment", fontsize=28)
plt.ylim(0.6,2.0)
ax.set_xlabel("Method", fontsize=20)
ax.set_ylabel("Enrichment ratio of eQTL", fontsize=20)
plt.subplots_adjust(left=0.18)
plt.subplots_adjust(bottom=0.25)
plt.savefig("eQTL_for_paper.png")
plt.show()


######################


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
plt.rcParams.update({'font.size': 18})
ax = sns.boxplot(data=df_concat, x=x, y=y, order=order,palette=my_pal,showfliers= False)
#for xtick in ax.get_xticks():
#    ax.text(xtick,medians[xtick] + vertical_offset,medians[xtick], 
#            horizontalalignment='center',size='x-small',color='w',weight='semibold')

pairs = [('Our model', 'Scenic+'),
         ('Our model', 'DIRECT-NET'),
         ('Our model', 'ArchR'),
         ('Our model', 'Cicero')]

annotator = Annotator(ax, pairs, data=df_concat, x=x, y=y, order=order)
annotator.configure(test='t-test_ind', text_format='star', loc='inside')
annotator.apply_and_annotate()
ax.set_xlabel("Method", fontsize=18)
ax.set_ylabel("Enrichment ratio of eQTL", fontsize=18)
plt.title("eQTL enrichment",fontsize=24)
plt.xticks(fontsize=18,rotation=90)
plt.yticks(fontsize=18)
plt.subplots_adjust(left=0.20)
plt.subplots_adjust(bottom=0.40)
plt.savefig("eqtl_for_paper.pdf",format="pdf")




################################

points = (eqtl_ratio_fantom,eqtl_ratio_pchic,eqtl_ratio_orig, eqtl_ratio_scenic,eqtl_ratio_direct,eqtl_ratio_archr,eqtl_ratio_cicero)
fig, ax = plt.subplots()
bp = ax.boxplot(points)
ax.set_xticklabels(["FANTOM5","PCHiC","Our model", "Scenic+", "DIRECT-NET", "ArchR", "Cicero"], fontsize=10)
#plt.title("The prediction of FANTOM5 enhancer", fontsize=18)
ax.set_xlabel("Method", fontsize=14)
ax.set_ylabel("Enrichment ratio of eQTL", fontsize=14)
plt.savefig("eQTL_for_paper_test.png")
plt.show()

#########################


eqtl=np.load("eQTL_for_paper/eqtl_mtx_all_cav.npy")
eqtl_sig=np.load("eQTL_for_paper/eqtl_mtx_sig05.npy")

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

#eqgene=eqtl_sig.sum(axis=1)
#enha_ans[eqgene==0,:]=-1


eqtl_ratio_fantom=((eqtl_sig[:,1:]*(enha[:,1:]>0)).sum(axis=1)/(eqtl[:,1:]*(enha[:,1:]>0)).sum(axis=1))/(eqtl_sig[:,1:].sum(axis=1)/eqtl[:,1:].sum(axis=1))
eqtl_ratio_fantom=eqtl_ratio_fantom[~(np.isnan(eqtl_ratio_fantom)|np.isinf(eqtl_ratio_fantom))].mean()

eqtl_ratio_pchic=((eqtl_sig[:,1:]*(pchic_matrix_th[:,1:]>0)).sum(axis=1)/(eqtl[:,1:]*(pchic_matrix_th[:,1:]>0)).sum(axis=1))/(eqtl_sig[:,1:].sum(axis=1)/eqtl[:,1:].sum(axis=1))
eqtl_ratio_pchic=eqtl_ratio_pchic[~(np.isnan(eqtl_ratio_pchic)|np.isinf(eqtl_ratio_pchic))].mean()


eqtl_ratio_orig=np.zeros(5)
for i in range(5):
    enha_tmp=grad_orig_th[:,:,i]
    tmp=((eqtl_sig[:,1:]*enha_tmp[:,1:]).sum(axis=1)/(eqtl[:,1:]*enha_tmp[:,1:]).sum(axis=1))/(eqtl_sig[:,1:].sum(axis=1)/eqtl[:,1:].sum(axis=1))
    eqtl_ratio_orig[i]=tmp[~(np.isnan(tmp)|np.isinf(tmp))].mean()

eqtl_ratio_archr=np.zeros(5)
for i in range(5):
    enha_tmp=grad_archr_th[:,:,i]
    tmp=((eqtl_sig[:,1:]*enha_tmp[:,1:]).sum(axis=1)/(eqtl[:,1:]*enha_tmp[:,1:]).sum(axis=1))/(eqtl_sig[:,1:].sum(axis=1)/eqtl[:,1:].sum(axis=1))
    eqtl_ratio_archr[i]=tmp[~(np.isnan(tmp)|np.isinf(tmp))].mean()


eqtl_ratio_direct=np.zeros(5)
for i in range(5):
    enha_tmp=grad_direct_th[:,:,i]
    tmp=((eqtl_sig[:,1:]*enha_tmp[:,1:]).sum(axis=1)/(eqtl[:,1:]*enha_tmp[:,1:]).sum(axis=1))/(eqtl_sig[:,1:].sum(axis=1)/eqtl[:,1:].sum(axis=1))
    eqtl_ratio_direct[i]=tmp[~(np.isnan(tmp)|np.isinf(tmp))].mean()


eqtl_ratio_cicero=np.zeros(5)
for i in range(5):
    enha_tmp=grad_cicero_th[:,:,i]
    tmp=((eqtl_sig[:,1:]*enha_tmp[:,1:]).sum(axis=1)/(eqtl[:,1:]*enha_tmp[:,1:]).sum(axis=1))/(eqtl_sig[:,1:].sum(axis=1)/eqtl[:,1:].sum(axis=1))
    eqtl_ratio_cicero[i]=tmp[~(np.isnan(tmp)|np.isinf(tmp))].mean()


eqtl_ratio_scenic=np.zeros(5)
for i in range(5):
    enha_tmp=grad_scenic_th[:,:,i]
    tmp=((eqtl_sig[:,1:]*enha_tmp[:,1:]).sum(axis=1)/(eqtl[:,1:]*enha_tmp[:,1:]).sum(axis=1))/(eqtl_sig[:,1:].sum(axis=1)/eqtl[:,1:].sum(axis=1))
    eqtl_ratio_scenic[i]=tmp[~(np.isnan(tmp)|np.isinf(tmp))].mean()



((eqtl_sig*enha_tmp).sum(axis=1)/

eqtl_ratio_orig=np.zeros(5)
for i in range(5):
    enha_tmp=grad_orig_th[:,1:,i]
    tmp=(eqtl_sig[:,1:]*enha_tmp).sum(axis=1)
    eqtl_ratio_orig[i]=tmp[~(np.isnan(tmp)|np.isinf(tmp))].mean()

eqtl_ratio_scenic=np.zeros(3)
for i in range(3):
    enha_tmp=grad_scenic_th[:,1:,i]
    tmp=(eqtl_sig*enha_tmp).sum(axis=1)
    eqtl_ratio_scenic[i]=tmp[~(np.isnan(tmp)|np.isinf(tmp))].mean()

eqtl_ratio_orig2=np.zeros(5)
for i in range(5):
    enha_tmp=grad_orig_th[:,1:,i]
    tmp=(eqtl*enha_tmp).sum(axis=1)
    eqtl_ratio_orig2[i]=tmp[~(np.isnan(tmp)|np.isinf(tmp))].mean()

eqtl_ratio_scenic2=np.zeros(5)
for i in range(5):
    enha_tmp=grad_scenic_th[:,1:,i]
    tmp=(eqtl*enha_tmp).sum(axis=1)
    eqtl_ratio_scenic2[i]=tmp[~(np.isnan(tmp)|np.isinf(tmp))].mean()

enha_ans


points = (eqtl_ratio_fantom,eqtl_ratio_pchic,eqtl_ratio_orig, eqtl_ratio_scenic,eqtl_ratio_direct,eqtl_ratio_archr,eqtl_ratio_cicero)
fig, ax = plt.subplots()
bp = ax.boxplot(points)
ax.set_xticklabels(["FANTOM5","PCHiC","Our model", "Scenic+", "DIRECT-NET", "ArchR", "Cicero"], fontsize=10)
#plt.title("The prediction of FANTOM5 enhancer", fontsize=18)
ax.set_xlabel("Method", fontsize=14)
ax.set_ylabel("Enrichment ratio of eQTL", fontsize=14)
plt.savefig("eQTL_for_paper.png")
plt.show()


points = (eqtl_ratio_fantom,eqtl_ratio_pchic,eqtl_ratio_orig, eqtl_ratio_scenic,eqtl_ratio_direct,eqtl_ratio_archr,eqtl_ratio_cicero)
fig, ax = plt.subplots()
bp = ax.boxplot(points)
ax.set_xticklabels(["FANTOM5","PCHiC","Our model", "Scenic+", "DIRECT-NET", "ArchR", "Cicero"], fontsize=10)
#plt.title("The prediction of FANTOM5 enhancer", fontsize=18)
ax.set_xlabel("Method", fontsize=14)
ax.set_ylabel("Enrichment ratio of eQTL", fontsize=14)
plt.savefig("eQTL_for_paper_test2.png")
plt.show()



###################

eqtl_ratio_orig_tmp=np.zeros(eqtl_sig.shape[0])
for i in range(eqtl_sig.shape[0]):
    eqtl_tmp=eqtl_sig[i,:]
    grad_tmp=grad_orig[i,:,0]
    tagtmp=enha_ans[i,:]
    tmp1=eqtl_tmp[tagtmp!=(-1)]
    tmp2=grad_tmp[tagtmp!=(-1)]
    eqtl_ratio_orig_tmp[i]=tmp2[tmp1==1].sum()/tmp2[tmp1==0].sum()

eqtl_ratio_orig_tmp[np.isnan(eqtl_ratio_orig_tmp)]=0
eqtl_ratio_orig_tmp[np.isinf(eqtl_ratio_orig_tmp)]=0


eqtl_sig_t1=np.load("eQTL_for_paper/eqtl_mtx_sig_t1.npy")
eqtl_sig_t2=np.load("eQTL_for_paper/eqtl_mtx_sig_t2.npy")


eqtl_sig_s=eqtl_sig_t2[enha_ans!=(-1)]
enha_s=enha[enha_ans!=(-1)]
(eqtl_sig_s[enha_s==1].mean())/(eqtl_sig_s.mean())

pchic_matrix_th
hic_ans=pchic_matrix_th.copy()
hic_ans[:,0]=-1


eqtl_s=eqtl[hic_ans!=(-1)]
eqtl_sig_s=eqtl_sig[hic_ans!=(-1)]
enha_s=hic_ans[hic_ans!=(-1)]
(eqtl_sig_s[enha_s==1].mean()/eqtl_s[enha_s==1].mean())/(eqtl_sig_s.mean()/eqtl_s.mean())

#####
###################

peaklist_full_signal=pd.read_csv("encode_for_paper/peaks_extend_signal.txt",sep="\t")
peaklist_full_signal=peaklist_full_signal.to_numpy()

chip_mtx=np.zeros((pairlist.shape[0],max_len,4))
chip_mtx=chip_mtx-1
for j in range(pairlist.shape[0]):
    peakstart=int(pairlist[j,0].item())
    peakend=int(pairlist[j,1].item())
    prompos=int(pairlist_prom[j,0].item())
    enha=list(range(peakstart,(peakend+1)))
    enha.remove(prompos)
    genepos=int(pairlist[j,2].item())
    peaknum_gene=peakend-peakstart+1
    chip_mtx[j,1:peaknum_gene,:]=peaklist_full_signal[enha,4:8]
    chip_mtx[j,0,:]=peaklist_full_signal[prompos,4:8]


grad_all=np.load("pbmc_fold1/Deeplift_full_ver2_all.npy")


cellclus=pd.read_csv("pbmc_fold1/cellcluster_all.txt",header=None)
cellclus=cellclus.to_numpy()
cellclus=cellclus[:,0]

grad_cellclus=np.zeros((grad_all.shape[0],grad_all.shape[2],5))
for i in range(5):
    grad_cellclus[:,:,i]=grad_all[:,cellclus==(i+1),:].mean(axis=1)

np.save("grad_cellclus.npy",grad_cellclus)

grad_orig[:,:,0]

tmpchip=chip_mtx[enha_ans!=(-1),0]
tmpgrad=grad_cellclus[enha_ans!=(-1),0]

normfac=chip_mtx[enha_ans!=(-1),:].mean(axis=0)
chip_mtx_norm=np.zeros(chip_mtx.shape)
for i in range(normfac.shape[0]):
    chip_mtx_norm[:,:,i]=chip_mtx[:,:,i]/normfac[i]


grad_mtx_norm=np.zeros(chip_mtx.shape)
grad_mtx_norm[:,:,0]=grad_cellclus[:,:,0]
grad_mtx_norm[:,:,1]=grad_cellclus[:,:,3]
grad_mtx_norm[:,:,2]=grad_cellclus[:,:,4]
grad_mtx_norm[:,:,3]=grad_cellclus[:,:,1]

tmpchip=chip_mtx_norm[enha_ans!=(-1),:]
tmpgrad=grad_mtx_norm[enha_ans!=(-1),:]

corr_chip_grad=np.zeros((tmpchip.shape[0]))
for i in range(tmpchip.shape[0]):
    corr_chip_grad[i]=np.corrcoef([tmpchip[i,:],tmpgrad[i,:]])[1,0]

np.nanmean(corr_chip_grad[enha_ans[enha_ans!=(-1)]==1])
np.nanmean(corr_chip_grad[enha_ans[enha_ans!=(-1)]==0])

np.nanmean(np.abs(corr_chip_grad[enha_ans[enha_ans!=(-1)]==1]))
np.nanmean(np.abs(corr_chip_grad[enha_ans[enha_ans!=(-1)]==0]))

tmpchip[enha_ans[enha_ans!=(-1)]==1,:]
tmpgrad[enha_ans[enha_ans!=(-1)]==1,:]

grad_orig[enha_ans!=(-1),0]


roc_auc_score(tmpchip.flatten()>2,tmpgrad.flatten())

roc_auc_score(tmpgrad.flatten()>0.8,tmpchip.flatten())

fig, ax = plt.subplots()
ax.scatter(tmpchip.flatten(),tmpgrad.flatten(),s=1)
plt.savefig("scatter_chip_vs_grad_alll.png")
plt.show()

from sklearn.metrics import confusion_matrix

tmpchip.argmax(axis=1)
tmpgrad.argmax(axis=1)

cm = confusion_matrix(tmpchip.argmax(axis=1), tmpgrad.argmax(axis=1))

tmp1=cm.sum(axis=1)/cm.sum()
tmp2=cm.sum(axis=0)
predmtx=np.matmul(tmp1[:,np.newaxis],tmp2[np.newaxis,:])
cm/predmtx


np.nanmean(corr_chip_grad[grad_orig[enha_ans!=(-1),0]>1])

np.corrcoef([tmpchip,tmpgrad])

fig, ax = plt.subplots()
ax.scatter(tmpchip,tmpgrad,s=1)
plt.savefig("scatter_chip_vs_grad.png")
plt.show()

###########

grad_r2g=np.load('R2G_matrix.npy')
#grad_r2g=grad_r2g[1090:2090,:]
enha=np.loadtxt('peak_enha_pbmc.txt')
enha[:,0]=-1
grad_r2g_s=grad_r2g[enha!=-1]
enha_s=enha[enha!=-1]
roc_auc_score(enha_s, grad_r2g_s)

precision, recall, thresholds = precision_recall_curve(enha_s, grad_r2g_s)
numerator = 2 * recall * precision
denom = recall + precision
f1_scores = np.divide(numerator, denom, out=np.zeros_like(denom), where=(denom!=0))
max_f1_r2g = np.max(f1_scores)
max_f1_thresh_r2g = thresholds[np.argmax(f1_scores)]
print(max_f1_r2g)

################

grad_r2gr=np.load('R2G_rho_matrix.npy')
grad_r2gr=grad_r2gr[1090:2090,:]
enha=np.loadtxt('peak_enha_pbmc.txt')
enha[:,0]=-1
grad_r2gr_s=grad_r2gr[enha!=-1]
enha_s=enha[enha!=-1]
roc_auc_score(enha_s, grad_r2gr_s)

precision, recall, thresholds = precision_recall_curve(enha_s, grad_r2gr_s)
numerator = 2 * recall * precision
denom = recall + precision
f1_scores = np.divide(numerator, denom, out=np.zeros_like(denom), where=(denom!=0))
max_f1_r2gr = np.max(f1_scores)
max_f1_thresh_r2gr = thresholds[np.argmax(f1_scores)]
print(max_f1_r2gr)

##########

grad_r2g_ar=np.load('R2G_absrho_matrix.npy')
#grad_r2g_ar=grad_r2g_ar[1090:2090,:]
enha=np.loadtxt('peak_enha_pbmc.txt')
enha[:,0]=-1
grad_r2g_ar_s=grad_r2g_ar[enha!=(-1)]
enha_s=enha[enha!=-1]
roc_auc_score(enha_s, grad_r2g_ar_s)
precision_recall_curve(enha_s, grad_r2g_ar_s)

precision, recall, thresholds = precision_recall_curve(enha_s, grad_r2g_ar_s)
numerator = 2 * recall * precision
denom = recall + precision
f1_scores = np.divide(numerator, denom, out=np.zeros_like(denom), where=(denom!=0))
max_f1_r2g_ar_s = np.max(f1_scores)
max_f1_thresh_r2g_ar_s= thresholds[np.argmax(f1_scores)]
print(max_f1_r2g_ar_s)


###############


grad_r2g_ar=np.load('R2G_absrho_matrix_0.npy')
#grad_r2g_ar=grad_r2g_ar[1090:2090,:]
enha=np.loadtxt('peak_enha_pbmc.txt')
enha[:,0]=-1
grad_r2g_ar_s=grad_r2g_ar[enha!=(-1)]
enha_s=enha[enha!=-1]
roc_auc_score(enha_s, grad_r2g_ar_s)
precision_recall_curve(enha_s, grad_r2g_ar_s)

precision, recall, thresholds = precision_recall_curve(enha_s, grad_r2g_ar_s)
numerator = 2 * recall * precision
denom = recall + precision
f1_scores = np.divide(numerator, denom, out=np.zeros_like(denom), where=(denom!=0))
max_f1_r2g_ar_s = np.max(f1_scores)
max_f1_thresh_r2g_ar_s= thresholds[np.argmax(f1_scores)]
print(max_f1_r2g_ar_s)

###############
grad_dire=pd.read_csv("peakenha_direct.txt",sep="\t",header=None).to_numpy()
grad_dire=pd.read_csv("peakenha_direct_prto.txt",sep="\t",header=None).to_numpy()
grad_dire=grad_dire[1090:2090,:]
enha_st=np.loadtxt('peak_enha_pbmc.txt')
enha_st=enha_st[1090:2090,:]
enha_st[:,0]=-1
grad_dire_s=grad_dire[enha_st!=-1]
enha_sts=enha_st[enha_st!=-1]
roc_auc_score(enha_sts, grad_dire_s)

precision, recall, thresholds = precision_recall_curve(enha_s, grad_dire_s)
numerator = 2 * recall * precision
denom = recall + precision
f1_scores = np.divide(numerator, denom, out=np.zeros_like(denom), where=(denom!=0))
max_f1_dire_s = np.max(f1_scores)
max_f1_thresh_dire_s= thresholds[np.argmax(f1_scores)]
print(max_f1_dire_s)

###############
grad_dire=pd.read_csv("peakenha_direct_prto.txt",sep="\t",header=None).to_numpy()
enha=np.loadtxt('peak_enha_pbmc.txt')
enha[:,0]=-1
grad_dire_s=grad_dire[enha!=-1]
enha_s=enha[enha!=-1]
roc_auc_score(enha_s, grad_dire_s)

precision, recall, thresholds = precision_recall_curve(enha_s, grad_dire_s)
numerator = 2 * recall * precision
denom = recall + precision
f1_scores = np.divide(numerator, denom, out=np.zeros_like(denom), where=(denom!=0))
max_f1_dire_s = np.max(f1_scores)
max_f1_thresh_dire_s= thresholds[np.argmax(f1_scores)]
print(max_f1_dire_s)

###############
grad_dire=pd.read_csv("peakenha_direct_0.txt",sep="\t",header=None).to_numpy()
enha=np.loadtxt('peak_enha_pbmc.txt')
enha[:,0]=-1
grad_dire_s=grad_dire[enha!=-1]
enha_s=enha[enha!=-1]
roc_auc_score(enha_s, grad_dire_s)

precision, recall, thresholds = precision_recall_curve(enha_s, grad_dire_s)
numerator = 2 * recall * precision
denom = recall + precision
f1_scores = np.divide(numerator, denom, out=np.zeros_like(denom), where=(denom!=0))
max_f1_dire_s = np.max(f1_scores)
max_f1_thresh_dire_s= thresholds[np.argmax(f1_scores)]
print(max_f1_dire_s)

#########


startpos=0
endpos=280000
binnum=20
dist=(endpos-startpos)/(binnum-1)
pairbin_abspos_r_all_se=np.zeros((binnum,2))

for i in range(binnum):
    tag=dist*i+startpos
    tag2=dist*(i+1)+startpos
    pairbin_abspos_r_all_se[i,0]=(tag+tag2)/2
    tmp=gradmm_all[(np.abs(posmtx)>=tag)&(np.abs(posmtx)<tag2)]
    tmpe=enha[(np.abs(posmtx)>=tag)&(np.abs(posmtx)<tag2)]
    tmp=tmp[tmpe!=-1]
    tmpe=tmpe[tmpe!=-1]
    if len(tmpe)>0:
        pairbin_abspos_r_all_se[i,1]=roc_auc_score(tmpe,tmp)


startpos=0
endpos=280000
binnum=20
dist=(endpos-startpos)/(binnum-1)
pairbin_abspos_r_all_se_r2g=np.zeros((binnum,2))

for i in range(binnum):
    tag=dist*i+startpos
    tag2=dist*(i+1)+startpos
    pairbin_abspos_r_all_se_r2g[i,0]=(tag+tag2)/2
    tmp=grad_r2g_ar[(np.abs(posmtx)>=tag)&(np.abs(posmtx)<tag2)]
    tmpe=enha[(np.abs(posmtx)>=tag)&(np.abs(posmtx)<tag2)]
    tmp=tmp[tmpe!=-1]
    tmpe=tmpe[tmpe!=-1]
    if len(tmpe)>0:
        pairbin_abspos_r_all_se_r2g[i,1]=roc_auc_score(tmpe,tmp)


posmtx_s=posmtx[1090:2090,:]

startpos=0
endpos=280000
binnum=20
dist=(endpos-startpos)/(binnum-1)
pairbin_abspos_r_all_se_dire=np.zeros((binnum,2))

for i in range(binnum):
    tag=dist*i+startpos
    tag2=dist*(i+1)+startpos
    pairbin_abspos_r_all_se_dire[i,0]=(tag+tag2)/2
    tmp=grad_dire[(np.abs(posmtx_s)>=tag)&(np.abs(posmtx_s)<tag2)]
    tmpe=enha_st[(np.abs(posmtx_s)>=tag)&(np.abs(posmtx_s)<tag2)]
    tmp=tmp[tmpe!=-1]
    tmpe=tmpe[tmpe!=-1]
    if len(tmpe)>0:
        pairbin_abspos_r_all_se_dire[i,1]=roc_auc_score(tmpe,tmp)

fig, ax = plt.subplots()
ax.plot(pairbin_abspos_r_all_se[:,0],pairbin_abspos_r_all_se[:,1],color="b")
ax.plot(pairbin_abspos_r_all_se_r2g[:,0],pairbin_abspos_r_all_se_r2g[:,1],color="g")
ax.plot(pairbin_abspos_r_all_se_dire[:,0],pairbin_abspos_r_all_se_dire[:,1],color="c")
ax.set_xlim(startpos,endpos)
#ax.set_ylim(-0.08,-0.04)
plt.savefig("AUROC_model_scenic_directnet_dist.png")
plt.show()


#############
from sklearn.metrics import f1_score
f1_score(pchic_matrix, gradmm)
pchic_matrix

pchic_matrix[gradmm>np.mean(gradmm)].mean()
pchic_matrix[gradmm<np.mean(gradmm)].mean()

pchic_matrix[grad_r2g_ar>np.mean(grad_r2g_ar)].mean()
pchic_matrix[grad_r2g_ar<np.mean(grad_r2g_ar)].mean()

c_array = np.percentile(pchic_matrix[enha!=(-1)], q=[90])
pchic_matrix_th=(pchic_matrix>c_array).astype(float)

gradmm_all_s=gradmm_all[enha!=-1]
gradmm_s=gradmm[enha!=-1]
pchic_matrix_th_s=pchic_matrix_th[enha!=-1]
roc_auc_score(pchic_matrix_th_s, gradmm_s)
roc_auc_score(pchic_matrix_th_s, gradmm_all_s)

grad_r2g_ar_s=grad_r2g_ar[enha!=-1]
roc_auc_score(pchic_matrix_th_s, grad_r2g_ar_s)

grad_r2g_s=grad_r2g[enha!=-1]
roc_auc_score(pchic_matrix_th_s, grad_r2g_s)

grad_dire_s=grad_dire[enha!=-1]
roc_auc_score(pchic_matrix_th_s, grad_dire_s)
###########

enha=np.loadtxt('peak_enha_pbmc.txt')
enha[:,0]=-1

startpos=0
endpos=280000
binnum=20
dist=(endpos-startpos)/(binnum-1)
pairbin_abspos_r_all_se_hic=np.zeros((binnum,2))

for i in range(binnum):
    tag=dist*i+startpos
    tag2=dist*(i+1)+startpos
    pairbin_abspos_r_all_se_hic[i,0]=(tag+tag2)/2
    tmp=gradmm_all[(np.abs(posmtx)>=tag)&(np.abs(posmtx)<tag2)]
    tmpe=pchic_matrix_th[(np.abs(posmtx)>=tag)&(np.abs(posmtx)<tag2)]
    enhatmp=enha[(np.abs(posmtx)>=tag)&(np.abs(posmtx)<tag2)]
    tmp=tmp[enhatmp!=-1]
    tmpe=tmpe[enhatmp!=-1]
    if len(tmpe)>0:
        pairbin_abspos_r_all_se_hic[i,1]=roc_auc_score(tmpe,tmp)


startpos=0
endpos=280000
binnum=20
dist=(endpos-startpos)/(binnum-1)
pairbin_abspos_r_all_se_r2g_hic=np.zeros((binnum,2))

for i in range(binnum):
    tag=dist*i+startpos
    tag2=dist*(i+1)+startpos
    pairbin_abspos_r_all_se_r2g_hic[i,0]=(tag+tag2)/2
    tmp=grad_r2g_ar[(np.abs(posmtx)>=tag)&(np.abs(posmtx)<tag2)]
    tmpe=pchic_matrix_th[(np.abs(posmtx)>=tag)&(np.abs(posmtx)<tag2)]
    enhatmp=enha[(np.abs(posmtx)>=tag)&(np.abs(posmtx)<tag2)]
    tmp=tmp[enhatmp!=-1]
    tmpe=tmpe[enhatmp!=-1]
    if len(tmpe)>0:
        pairbin_abspos_r_all_se_r2g_hic[i,1]=roc_auc_score(tmpe,tmp)


posmtx_s=posmtx[1090:2090,:]
pchic_matrix_th_st=pchic_matrix_th[1090:2090,:]

startpos=0
endpos=280000
binnum=20
dist=(endpos-startpos)/(binnum-1)
pairbin_abspos_r_all_se_dire_hic=np.zeros((binnum,2))

for i in range(binnum):
    tag=dist*i+startpos
    tag2=dist*(i+1)+startpos
    pairbin_abspos_r_all_se_dire_hic[i,0]=(tag+tag2)/2
    tmp=grad_dire[(np.abs(posmtx_s)>=tag)&(np.abs(posmtx_s)<tag2)]
    tmpe=pchic_matrix_th_st[(np.abs(posmtx_s)>=tag)&(np.abs(posmtx_s)<tag2)]
    enhatmp_s=enha_st[(np.abs(posmtx_s)>=tag)&(np.abs(posmtx_s)<tag2)]
    tmp=tmp[enhatmp_s!=-1]
    tmpe=tmpe[enhatmp_s!=-1]
    if len(tmpe)>0:
        pairbin_abspos_r_all_se_dire_hic[i,1]=roc_auc_score(tmpe,tmp)

fig, ax = plt.subplots()
ax.plot(pairbin_abspos_r_all_se_hic[:,0],pairbin_abspos_r_all_se_hic[:,1],color="b")
ax.plot(pairbin_abspos_r_all_se_r2g_hic[:,0],pairbin_abspos_r_all_se_r2g_hic[:,1],color="g")
ax.plot(pairbin_abspos_r_all_se_dire_hic[:,0],pairbin_abspos_r_all_se_dire_hic[:,1],color="c")
ax.set_xlim(startpos,endpos)
#ax.set_ylim(-0.08,-0.04)
plt.savefig("AUROC_model_scenic_directnet_dist_hic.png")
plt.show()

###############
enha=np.loadtxt('peak_enha_pbmc.txt')
enha[:,0]=-1

startpos=0
endpos=280000
binnum=20
dist=(endpos-startpos)/(binnum-1)
pchic_matrix_normth=np.zeros(pchic_matrix.shape)

for i in range(binnum):
    tag=dist*i+startpos
    tag2=dist*(i+1)+startpos
    pairbin_abspos_r_all_se_hic[i,0]=(tag+tag2)/2
    tmp=gradmm[(np.abs(posmtx)>=tag)&(np.abs(posmtx)<tag2)]
    tmpe=pchic_matrix[(np.abs(posmtx)>=tag)&(np.abs(posmtx)<tag2)]
    enhatmp=enha[(np.abs(posmtx)>=tag)&(np.abs(posmtx)<tag2)]
    c_array = np.percentile(tmpe[enhatmp!=(-1)], q=[90])
    pchic_matrix_normth[(np.abs(posmtx)>=tag)&(np.abs(posmtx)<tag2)]=(tmpe>c_array).astype(float)


gradmm_s=gradmm[enha!=-1]
pchic_matrix_normth_s=pchic_matrix_normth[enha!=-1]
roc_auc_score(pchic_matrix_normth_s, gradmm_s)

grad_r2g_ar_s=grad_r2g_ar[enha!=-1]
roc_auc_score(pchic_matrix_normth_s, grad_r2g_ar_s)


#############
enha=np.loadtxt('peak_enha_pbmc.txt')
enha[:,0]=-1


startpos=0
endpos=280000
binnum=20
dist=(endpos-startpos)/(binnum-1)
pairbin_abspos_r_all_se_hic=np.zeros((binnum,2))

for i in range(binnum):
    tag=dist*i+startpos
    tag2=dist*(i+1)+startpos
    pairbin_abspos_r_all_se_hic[i,0]=(tag+tag2)/2
    tmp=gradmm[(np.abs(posmtx)>=tag)&(np.abs(posmtx)<tag2)]
    tmpe=pchic_matrix[(np.abs(posmtx)>=tag)&(np.abs(posmtx)<tag2)]
    enhatmp=enha[(np.abs(posmtx)>=tag)&(np.abs(posmtx)<tag2)]
    c_array = np.percentile(tmpe[enhatmp!=(-1)], q=[90])
    tmpe_th=(tmpe>c_array).astype(float)
    tmp=tmp[enhatmp!=-1]
    tmpe_th=tmpe_th[enhatmp!=-1]
    if len(tmpe_th)>0:
        pairbin_abspos_r_all_se_hic[i,1]=roc_auc_score(tmpe_th,tmp)


startpos=0
endpos=280000
binnum=20
dist=(endpos-startpos)/(binnum-1)
pairbin_abspos_r_all_se_r2g_hic=np.zeros((binnum,2))

for i in range(binnum):
    tag=dist*i+startpos
    tag2=dist*(i+1)+startpos
    pairbin_abspos_r_all_se_r2g_hic[i,0]=(tag+tag2)/2
    tmp=grad_r2g_ar[(np.abs(posmtx)>=tag)&(np.abs(posmtx)<tag2)]
    tmpe=pchic_matrix[(np.abs(posmtx)>=tag)&(np.abs(posmtx)<tag2)]
    enhatmp=enha[(np.abs(posmtx)>=tag)&(np.abs(posmtx)<tag2)]
    c_array = np.percentile(tmpe[enhatmp!=(-1)], q=[90])
    tmpe_th=(tmpe>c_array).astype(float)
    tmp=tmp[enhatmp!=-1]
    tmpe_th=tmpe_th[enhatmp!=-1]
    if len(tmpe)>0:
        pairbin_abspos_r_all_se_r2g_hic[i,1]=roc_auc_score(tmpe_th,tmp)


posmtx_s=posmtx[1090:2090,:]
pchic_matrix_th_st=pchic_matrix_th[1090:2090,:]

startpos=0
endpos=280000
binnum=20
dist=(endpos-startpos)/(binnum-1)
pairbin_abspos_r_all_se_dire_hic=np.zeros((binnum,2))

for i in range(binnum):
    tag=dist*i+startpos
    tag2=dist*(i+1)+startpos
    pairbin_abspos_r_all_se_dire_hic[i,0]=(tag+tag2)/2
    tmp=grad_dire[(np.abs(posmtx_s)>=tag)&(np.abs(posmtx_s)<tag2)]
    tmpe=pchic_matrix_th_st[(np.abs(posmtx_s)>=tag)&(np.abs(posmtx_s)<tag2)]
    tmp=tmp[tmpe!=-1]
    tmpe=tmpe[tmpe!=-1]
    if len(tmpe)>0:
        pairbin_abspos_r_all_se_dire_hic[i,1]=roc_auc_score(tmpe,tmp)

fig, ax = plt.subplots()
ax.plot(pairbin_abspos_r_all_se_hic[:,0],pairbin_abspos_r_all_se_hic[:,1],color="b")
ax.plot(pairbin_abspos_r_all_se_r2g_hic[:,0],pairbin_abspos_r_all_se_r2g_hic[:,1],color="g")
ax.plot(pairbin_abspos_r_all_se_dire_hic[:,0],pairbin_abspos_r_all_se_dire_hic[:,1],color="c")
ax.set_xlim(startpos,endpos)
#ax.set_ylim(-0.08,-0.04)
plt.savefig("AUROC_model_scenic_directnet_dist_hic.png")
plt.show()



#############


peaklist=pd.read_csv('peaks.bed',sep="\t",header=None)
peaklist=peaklist[[1,2]].to_numpy()
pairlist=pd.read_csv('peak_gene_pair_promoter300000_magic.csv',header=None)
pairlist=pairlist[[1,2,3]].to_numpy()
#pairlist=pairlist[1090:,:]
pairlist_prom=pd.read_csv('peak_gene_pair_promoter_magic.csv',header=None)
pairlist_prom=pairlist_prom[[1,2,3]].to_numpy()
#pairlist_prom=pairlist_prom[1090:,:]

peakidmtx=np.zeros((2090,80))
for j in range(2090):
    peakstart=int(pairlist[j,0].item())
    peakend=int(pairlist[j,1].item())
    prompos=int(pairlist_prom[j,0].item())
    enha=list(range(peakstart,(peakend+1)))
    enha.remove(prompos)
    peaknum_gene=peakend-peakstart+1
    peakidmtx[j,0]=prompos
    peakidmtx[j,1:peaknum_gene]=enha

posmtx=np.zeros((2090,80))
for j in range(2090):
    peakstart=int(pairlist[j,0].item())
    peakend=int(pairlist[j,1].item())
    prompos=int(pairlist_prom[j,0].item())
    enha=list(range(peakstart,(peakend+1)))
    enha.remove(prompos)
    genepos=int(pairlist[j,2].item())
    peaknum_gene=peakend-peakstart+1
    posvec=((((peaklist[enha,0]+peaklist[enha,1])/2)-genepos))
    posmtx[j,1:peaknum_gene]=posvec



############################
promenhatag=np.load("promenhatag_from_synanalysis2.npy")
grad=np.load('DeepLift_ver16_3_18_correct3_50dim_pca_raw_log1p_fullgene.npy')
cellclus=np.loadtxt("cellcluster.txt")


grad_clus=np.zeros((grad.shape[0],5,grad.shape[2]))
for i in range(5):
    grad_clus[:,i,:]=grad[:,cellclus==i,:].mean(axis=1)

fantom_cd4_matrix=np.load("fantom_cd4_matrix.npy")
fantom_cd8_matrix=np.load("fantom_cd8_matrix.npy")
fantom_b_matrix=np.load("fantom_b_matrix.npy")
fantom_mono_matrix=np.load("fantom_mono_matrix.npy")

fantom_clus=np.zeros((grad.shape[0],5,grad.shape[2]))
fantom_clus[:,0,:]=fantom_mono_matrix
fantom_clus[:,1,:]=fantom_cd4_matrix
fantom_clus[:,2,:]=fantom_cd4_matrix
fantom_clus[:,3,:]=fantom_cd8_matrix
fantom_clus[:,4,:]=fantom_b_matrix

fantom_grad_corr_matrix=np.zeros((2090,80))
for gn in range(2090):
    print(gn)
    peaknum=(promenhatag[gn,:]!=0).sum()
    for p1 in range(peaknum):
        fantom_grad_corr_matrix[gn,p1]=np.corrcoef([fantom_clus[gn,:,p1],grad_clus[gn,:,p1]])[0,1]


fantom_grad_corr_matrix_sp=np.zeros((2090,80))
for gn in range(2090):
    print(gn)
    peaknum=(promenhatag[gn,:]!=0).sum()
    for p1 in range(peaknum):
        correlation, pvalue =spearmanr(fantom_clus[gn,:,p1],grad_clus[gn,:,p1])
        fantom_grad_corr_matrix_sp[gn,p1]=correlation


fantom_tmp=fantom_clus.transpose(0,2,1)[promenhatag!=0].flatten()
grad_tmp=grad_clus.transpose(0,2,1)[promenhatag!=0].flatten()

np.corrcoef([fantom_tmp,grad_tmp])

fantom_grad_corr_matrix[np.isnan(fantom_grad_corr_matrix)]=0


fantom_grad_corr_matrix[promenhatag==4].mean()

np.nanmean(fantom_grad_corr_matrix_sp[promenhatag==4])

np.nanmean(fantom_grad_corr_matrix_sp[(promenhatag==4)|(promenhatag==2)])


##########

fantom_grad_corr_matrix=np.zeros((2090,80))
for gn in range(2090):
    print(gn)
    peaknum=(promenhatag[gn,:]!=0).sum()
    for p1 in range(peaknum):
        fantom_grad_corr_matrix[gn,p1]=np.corrcoef([fantom_clus[gn,:,p1],grad_clus[gn,:,p1]])[0,1]


###########

ATAC_use=np.load("ATAC_use_full.npy")

grad_attr_corr_matrix=np.zeros((2090,80))
for gn in range(2090):
    print(gn)
    peaknum=(promenhatag[gn,:]!=0).sum()
    for p1 in range(peaknum):
        grad_attr_corr_matrix[gn,p1]=np.corrcoef([ATAC_use[gn,:,p1],grad[gn,:,p1]])[0,1]


for gn in range(2090):
    print(gn)
    peaknum=(promenhatag[gn,:]!=0).sum()
    for p1 in range(peaknum):
        if promenhatag[gn,p1]>4:
            fig, ax = plt.subplots()
            ax.scatter(ATAC_use[gn,:,p1],grad[gn,:,p1],s=1)
            fname="enha_test/"+str(gn)+"_"+str(p1)+".png"
            plt.savefig(fname)
            plt.show()

        grad_attr_corr_matrix[gn,p1]=np.corrcoef([ATAC_use[gn,:,p1],grad[gn,:,p1]])[0,1]

grad_attr_corr_matrix=np.zeros((2090,80))
for gn in range(2090):
    print(gn)
    peaknum=(promenhatag[gn,:]!=0).sum()
    for p1 in range(peaknum):
        correlation, pvalue =spearmanr(ATAC_use[gn,:,p1],grad[gn,:,p1])
        grad_attr_corr_matrix[gn,p1]=correlation

correlation, pvalue =spearmanr(fantom_clus[gn,:,p1],grad_clus[gn,:,p1])

grad_attr_corr_matrix=np.zeros((1000,80))
for gn in range(1000):
    print(gn)
    peaknum=(promenhatag[gn,:]!=0).sum()
    for p1 in range(peaknum):
        grad_attr_corr_matrix[gn,p1]=np.corrcoef([ATAC_use[(1090+gn),:,p1],grad[gn,:,p1]])[0,1]


(fantom_grad_corr_matrix[(promenhatag==4)&(grad_attr_corr_matrix<0)]<0).mean()
(fantom_grad_corr_matrix[(promenhatag==4)&(grad_attr_corr_matrix>0)]<0).mean()

grad_attr_corr_matrix[]

###########################

promenhatag=np.load("promenhatag_from_synanalysis2.npy")

peaklist=pd.read_csv('peaks.bed',sep="\t",header=None)
peaklist=peaklist[[1,2]].to_numpy()
pairlist=pd.read_csv('peak_gene_pair_promoter300000_magic.csv',header=None)
pairlist=pairlist[[1,2,3]].to_numpy()
#pairlist=pairlist[1090:,:]
pairlist_prom=pd.read_csv('peak_gene_pair_promoter_magic.csv',header=None)
pairlist_prom=pairlist_prom[[1,2,3]].to_numpy()
#pairlist_prom=pairlist_prom[1090:,:]

peakidmtx=np.zeros((2090,80))
for j in range(2090):
    peakstart=int(pairlist[j,0].item())
    peakend=int(pairlist[j,1].item())
    prompos=int(pairlist_prom[j,0].item())
    enha=list(range(peakstart,(peakend+1)))
    enha.remove(prompos)
    peaknum_gene=peakend-peakstart+1
    peakidmtx[j,0]=prompos
    peakidmtx[j,1:peaknum_gene]=enha


fantom_data=pd.read_csv('CD4_fantom_ctss.bed',sep="\t",header=None)
fantom_data=fantom_data.iloc[:,[0,1,2,4]]
fantom_data=fantom_data.drop_duplicates()
fantom_data=fantom_data.to_numpy()

peaklist_full=pd.read_csv('peaks_extend.bed',sep="\t",header=None)
peaklist_full=peaklist_full.to_numpy()

fantom_cd4_matrix=np.zeros((2090,80))
for gn in range(2090):
    print(gn)
    peaknum=(promenhatag[gn,:]!=0).sum()
    for p1 in range(peaknum):
        p1_id=peakidmtx[gn,p1].astype(int)
        tmphic=fantom_data[fantom_data[:,0]==peaklist_full[p1_id,0],:].copy()
        p1_in=(tmphic[:,1]>=peaklist_full[p1_id,1])&(tmphic[:,2]<=peaklist_full[p1_id,2])
        maemuki=tmphic[(p1_in),:]
        fantom_cd4_matrix[gn,p1]=maemuki[:,3].sum()

np.save("fantom_cd4_matrix.npy",fantom_cd4_matrix)



fantom_data=pd.read_csv('CD8_fantom_ctss.bed',sep="\t",header=None)
fantom_data=fantom_data.iloc[:,[0,1,2,4]]
fantom_data=fantom_data.drop_duplicates()
fantom_data=fantom_data.to_numpy()

peaklist_full=pd.read_csv('peaks_extend.bed',sep="\t",header=None)
peaklist_full=peaklist_full.to_numpy()

fantom_cd8_matrix=np.zeros((2090,80))
for gn in range(2090):
    print(gn)
    peaknum=(promenhatag[gn,:]!=0).sum()
    for p1 in range(peaknum):
        p1_id=peakidmtx[gn,p1].astype(int)
        tmphic=fantom_data[fantom_data[:,0]==peaklist_full[p1_id,0],:].copy()
        p1_in=(tmphic[:,1]>=peaklist_full[p1_id,1])&(tmphic[:,2]<=peaklist_full[p1_id,2])
        maemuki=tmphic[(p1_in),:]
        fantom_cd8_matrix[gn,p1]=maemuki[:,3].sum()

np.save("fantom_cd8_matrix.npy",fantom_cd8_matrix)


fantom_data=pd.read_csv('CD19_fantom_ctss.bed',sep="\t",header=None)
fantom_data=fantom_data.iloc[:,[0,1,2,4]]
fantom_data=fantom_data.drop_duplicates()
fantom_data=fantom_data.to_numpy()

peaklist_full=pd.read_csv('peaks_extend.bed',sep="\t",header=None)
peaklist_full=peaklist_full.to_numpy()

fantom_b_matrix=np.zeros((2090,80))
for gn in range(2090):
    print(gn)
    peaknum=(promenhatag[gn,:]!=0).sum()
    for p1 in range(peaknum):
        p1_id=peakidmtx[gn,p1].astype(int)
        tmphic=fantom_data[fantom_data[:,0]==peaklist_full[p1_id,0],:].copy()
        p1_in=(tmphic[:,1]>=peaklist_full[p1_id,1])&(tmphic[:,2]<=peaklist_full[p1_id,2])
        maemuki=tmphic[(p1_in),:]
        fantom_b_matrix[gn,p1]=maemuki[:,3].sum()

np.save("fantom_b_matrix.npy",fantom_b_matrix)




fantom_data=pd.read_csv('CD14_fantom_ctss.bed',sep="\t",header=None)
fantom_data=fantom_data.iloc[:,[0,1,2,4]]
fantom_data=fantom_data.drop_duplicates()
fantom_data=fantom_data.to_numpy()

peaklist_full=pd.read_csv('peaks_extend.bed',sep="\t",header=None)
peaklist_full=peaklist_full.to_numpy()

fantom_mono_matrix=np.zeros((2090,80))
for gn in range(2090):
    print(gn)
    peaknum=(promenhatag[gn,:]!=0).sum()
    for p1 in range(peaknum):
        p1_id=peakidmtx[gn,p1].astype(int)
        tmphic=fantom_data[fantom_data[:,0]==peaklist_full[p1_id,0],:].copy()
        p1_in=(tmphic[:,1]>=peaklist_full[p1_id,1])&(tmphic[:,2]<=peaklist_full[p1_id,2])
        maemuki=tmphic[(p1_in),:]
        fantom_mono_matrix[gn,p1]=maemuki[:,3].sum()

np.save("fantom_mono_matrix.npy",fantom_mono_matrix)