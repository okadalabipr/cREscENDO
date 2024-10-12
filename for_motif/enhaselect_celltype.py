import numpy as np
import pandas as pd
import sys

args = sys.argv

samplename=str(args[1])
celltype=int(str(args[2]))


promenhatag=np.load(samplename+"/promenhatag.npy")

dar=np.load(samplename+"/enha_clus_mtx_deeplift.npy") #Cls,B,L
dar=dar[celltype]
genelist=pd.read_csv(samplename+"/pair_300000.csv",header=None)
genelist=genelist[[0]].to_numpy()
genelist=genelist[:,0]

fname=samplename+"/clus"+str(celltype)+"gene.txt"
tgtgene=pd.read_csv(fname,sep="\t",header=None)
tgttag=np.isin(genelist,tgtgene)
tgtmtx=np.ones(promenhatag.shape)
tgtmtx=(tgtmtx.transpose(1,0)*(tgttag.astype(int))).transpose(1,0)

enhaclus=(tgtmtx==1)&(dar==1)

outpath=samplename+"/enhaclus"+str(celltype)+".npy"

np.save(outpath,enhaclus)

