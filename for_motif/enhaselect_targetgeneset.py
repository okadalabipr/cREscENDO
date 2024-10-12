import numpy as np
import pandas as pd
import sys

args = sys.argv

samplename=str(args[1])
genesetpath=str(args[2])
outprefix=str(args[3])


promenhatag=np.load(samplename+"/promenhatag.npy")

genelist=pd.read_csv(samplename+"/pair_300000.csv",header=None)
genelist=genelist[[0]].to_numpy()
genelist=genelist[:,0]

tgtgene=pd.read_csv(genesetpath,sep="\t",header=None)
tgttag=np.isin(genelist,tgtgene)
tgtmtx=np.ones(promenhatag.shape)
tgtmtx=(tgtmtx.transpose(1,0)*(tgttag.astype(int))).transpose(1,0)

enhaclus=(tgtmtx==1)&((promenhatag==4)|(promenhatag==2))

outpath=samplename+"/enhaclus_"+outprefix+".npy"

np.save(outpath,enhaclus)

