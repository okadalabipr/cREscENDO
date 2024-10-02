import pandas as pd
import numpy as np
import sys

args = sys.argv
samplename=str(args[1])


peaks=pd.read_table(samplename+"/peaks.bed",delimiter="\t",header=None,names=["chr","start","end"])
peaks["center"]=(peaks["start"]+peaks["end"])/2
peaks["center"]=peaks["center"].astype('int')
peaks["start"]=peaks["center"]-750
peaks["end"]=peaks["center"]+750
peaks.to_csv(samplename+"/peaks_extend.bed",sep="\t",header=False, index=False)
