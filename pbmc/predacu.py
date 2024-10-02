import numpy as np
import pandas as pd
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import f1_score
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import scipy.stats as stats
import seaborn as sns


samplename="pbmcnew"

genelist=pd.read_csv(samplename+"/pair_300000.csv",sep=",",header=None)


clus1gene=pd.read_csv(samplename+"/clus1gene.txt",sep="\t",header=None)
clus2gene=pd.read_csv(samplename+"/clus2gene.txt",sep="\t",header=None)
clus3gene=pd.read_csv(samplename+"/clus3gene.txt",sep="\t",header=None)
clus4gene=pd.read_csv(samplename+"/clus4gene.txt",sep="\t",header=None)
clus5gene=pd.read_csv(samplename+"/clus5gene.txt",sep="\t",header=None)
clus6gene=pd.read_csv(samplename+"/clus6gene.txt",sep="\t",header=None)


genein=genelist[0].isin(clus1gene[0])+genelist[0].isin(clus2gene[0])+genelist[0].isin(clus3gene[0])+genelist[0].isin(clus4gene[0])+genelist[0].isin(clus5gene[0])+genelist[0].isin(clus6gene[0])

pred=np.load("pbmcnew/predrna_t2_raw_2_fold0.npy")


RNAmatrix=pd.read_csv(samplename+"/log1praw.csv",sep=",",header=None)
RNAmatrix=RNAmatrix.to_numpy()
RNAmatrix_b=((RNAmatrix.transpose(1,0)-RNAmatrix.mean(axis=1))/RNAmatrix.std(axis=1)).transpose(1,0)




obser=pd.read_csv("pbmcnew/RNAmatrix_raw.txt",header=None,sep=" ")
obser=obser.to_numpy().transpose(1,0)

from sklearn.metrics import r2_score
from scipy.stats import spearmanr

corrmtx=np.zeros((pred.shape[0]))
for i in range(corrmtx.shape[0]):
    #corrmtx[i]=np.corrcoef([pred[i,:],RNAmatrix_b[i,:]])[1,0]
    corrmtx[i],_=spearmanr(pred[i,:],RNAmatrix_b[i,:])
    #corrmtx[i]=r2_score(pred[i,:],RNAmatrix_b[i,:])

corrmtx[genein==0].mean()
corrmtx[genein==1].mean()


from scipy import stats
def fisher_z_transform(r):
    return 0.5 * np.log((1 + r) / (1 - r))


tmpa=fisher_z_transform(corrmtx[genein==1])
tmpb=fisher_z_transform(corrmtx[genein==0])

stats.ttest_ind(tmpa, tmpb)

corrmtxpd1=pd.DataFrame(corrmtx[genein==1])
corrmtxpd1.columns=["Correlation"]
corrmtxpd1["Type"]="Celltype"
corrmtxpd2=pd.DataFrame(corrmtx[genein==0])
corrmtxpd2.columns=["Correlation"]
corrmtxpd2["Type"]="noncelltype"

df_concat = pd.concat([corrmtxpd1,corrmtxpd2], axis=0)

from statannotations.Annotator import Annotator

x = "Type"
y = "Correlation"
order = ["Celltype","noncelltype"]
my_pal = {"Celltype": "blue", "noncelltype": "green"}
plt.clf()
plt.figure(figsize=(6.5,8))
plt.rcParams.update({'font.size': 18})
ax = sns.boxplot(data=df_concat, x=x, y=y, order=order,palette=my_pal,showfliers= False)
pairs = [("Celltype","noncelltype")]
annotator = Annotator(ax, pairs, data=df_concat, x=x, y=y, order=order)
annotator.configure(test='t-test_ind', text_format='star', loc='inside')
annotator.apply_and_annotate()
ax.set_xlabel("Type", fontsize=18)
ax.set_ylabel("Correlation coefficient", fontsize=18)
plt.title("Prediction accuracy",fontsize=24)
plt.xticks(fontsize=18,rotation=90)
plt.yticks(fontsize=18)
plt.subplots_adjust(left=0.20)
plt.subplots_adjust(bottom=0.40)
plt.savefig("pbmcnew/predacu.pdf",format="pdf")



enha=np.load("pbmc_fold1/enha.npy")
enha_ans=enha.copy()
enha_ans[:,0]=-1
pchic_matrix=np.load("pbmc_fold1/pchic_matrix.npy")
c_array = np.percentile(pchic_matrix[enha_ans!=(-1)], q=[90])
pchic_matrix_th=(pchic_matrix>c_array).astype(float)
pchic_matrix_th[enha_ans==(-1)]=-1


np.nanmean((enha_ans[genein==1]==1).sum(axis=1)/(enha_ans[genein==1]>=0).sum(axis=1))
np.nanmean((enha_ans[genein==0]==1).sum(axis=1)/(enha_ans[genein==0]>=0).sum(axis=1))


np.nanmean((enha_ans[genein==1]==1).sum(axis=1))
np.nanmean((enha_ans[genein==0]==1).sum(axis=1))


np.nanmean((pchic_matrix_th[genein==1]==1).sum(axis=1))
np.nanmean((pchic_matrix_th[genein==0]==1).sum(axis=1))


np.nanmean((pchic_matrix_th[genein==1]==1).sum(axis=1)/(pchic_matrix_th[genein==1]>=0).sum(axis=1))
np.nanmean((pchic_matrix_th[genein==0]==1).sum(axis=1)/(pchic_matrix_th[genein==0]>=0).sum(axis=1))





corrmtxpd1=pd.DataFrame((enha_ans[genein==1]==1).sum(axis=1))
corrmtxpd1.columns=["Correlation"]
corrmtxpd1["Type"]="Celltype"
corrmtxpd2=pd.DataFrame((enha_ans[genein==0]==1).sum(axis=1))
corrmtxpd2.columns=["Correlation"]
corrmtxpd2["Type"]="noncelltype"

df_concat = pd.concat([corrmtxpd1,corrmtxpd2], axis=0)

from statannotations.Annotator import Annotator

x = "Type"
y = "Correlation"
order = ["Celltype","noncelltype"]
my_pal = {"Celltype": "blue", "noncelltype": "green"}
plt.clf()
plt.figure(figsize=(6.5,8))
plt.rcParams.update({'font.size': 18})
ax = sns.boxplot(data=df_concat, x=x, y=y, order=order,palette=my_pal,showfliers= False)
pairs = [("Celltype","noncelltype")]
annotator = Annotator(ax, pairs, data=df_concat, x=x, y=y, order=order)
annotator.configure(test='t-test_ind', text_format='star', loc='inside')
annotator.apply_and_annotate()
ax.set_xlabel("Type", fontsize=18)
ax.set_ylabel("The number of FANTOM5 enhancers", fontsize=18)
plt.title("FANTOM5",fontsize=24)
plt.xticks(fontsize=18,rotation=90)
plt.yticks(fontsize=18)
plt.subplots_adjust(left=0.20)
plt.subplots_adjust(bottom=0.40)
plt.savefig("pbmcnew/predacu_fantom5.pdf",format="pdf")







corrmtxpd1=pd.DataFrame((pchic_matrix_th[genein==1]==1).sum(axis=1))
corrmtxpd1.columns=["Correlation"]
corrmtxpd1["Type"]="Celltype"
corrmtxpd2=pd.DataFrame((pchic_matrix_th[genein==0]==1).sum(axis=1))
corrmtxpd2.columns=["Correlation"]
corrmtxpd2["Type"]="noncelltype"

df_concat = pd.concat([corrmtxpd1,corrmtxpd2], axis=0)

from statannotations.Annotator import Annotator

x = "Type"
y = "Correlation"
order = ["Celltype","noncelltype"]
my_pal = {"Celltype": "blue", "noncelltype": "green"}
plt.clf()
plt.figure(figsize=(6.5,8))
plt.rcParams.update({'font.size': 18})
ax = sns.boxplot(data=df_concat, x=x, y=y, order=order,palette=my_pal,showfliers= False)
pairs = [("Celltype","noncelltype")]
annotator = Annotator(ax, pairs, data=df_concat, x=x, y=y, order=order)
annotator.configure(test='t-test_ind', text_format='star', loc='inside')
annotator.apply_and_annotate()
ax.set_xlabel("Type", fontsize=18)
ax.set_ylabel("The number of PCHiC interacting regions", fontsize=18)
plt.title("PCHiC",fontsize=24)
plt.xticks(fontsize=18,rotation=90)
plt.yticks(fontsize=18)
plt.subplots_adjust(left=0.20)
plt.subplots_adjust(bottom=0.40)
plt.savefig("pbmcnew/predacu_pchic.pdf",format="pdf")

