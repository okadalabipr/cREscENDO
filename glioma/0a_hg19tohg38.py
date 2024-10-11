import numpy as np
import pandas as pd
import sys

args = sys.argv
samplename=str(args[1])

peaklist_hg19=pd.read_csv(samplename+"/peaks_hg19_tag.bed",sep="\t",header=None,names=("chr","start","end","tag"))
peaklist_hg38=pd.read_csv(samplename+"/peaks_hg38_tag.bed",sep="\t",header=None,names=("chr","start","end","tag"))
chrs = ['chr'+str(i) for i in range(1,23)] + ['chrX', 'chrY']
peaklist_hg38_s=peaklist_hg38[peaklist_hg38["chr"].isin(chrs)]
ataccount=np.load(samplename+"/atac_count_hg19.npy")
ataccount_s=ataccount[peaklist_hg19["tag"].isin(peaklist_hg38_s["tag"]),:]
np.save(samplename+"/atac_count.npy",ataccount_s)
peaklist_hg38_s[["chr","start","end"]].to_csv(samplename+"/peaks.bed", sep="\t",header=False, index=False)