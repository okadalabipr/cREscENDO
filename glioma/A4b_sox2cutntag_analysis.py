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


tmp30b=pd.read_csv(samplename+"/tmp30b.bed",sep="\t",header=None)
tmp31b=pd.read_csv(samplename+"/tmp31b.bed",sep="\t",header=None)

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
plt.rcParams.update({"font.size": 18})
ax = sns.boxplot(data=df_concat, x=x, y=y, order=order,palette=my_pal,showfliers= False)
pairs = [("Tumor","Non-enhancer")]
annotator = Annotator(ax, pairs, data=df_concat, x=x, y=y, order=order)
annotator.configure(test="t-test_ind", text_format="star", loc="inside")
annotator.apply_and_annotate()
ax.set_xlabel("Enhancer", fontsize=18)
ax.set_ylabel("log2FC of SOX2 signal from background (Cut-n-Tag)", fontsize=18)
plt.title("SOX2 binding in glioblastoma",fontsize=24)
plt.xticks(fontsize=18,rotation=90)
plt.yticks(fontsize=18)
plt.subplots_adjust(left=0.20)
plt.subplots_adjust(bottom=0.40)
plt.savefig(samplename+"/sox2bind.pdf",format="pdf")

