import numpy as np
import pandas as pd
import os
import sys

args = sys.argv
samplename=str(args[1])

filelis=os.listdir(samplename)

alllist=[]
i=0
print(i)
tmplis=[]
tagname="Deeplift_full_ver2_fold"+str(i)
l_in = [s for s in filelis if tagname in s]
for j in range(len(l_in)):
  fname=samplename+"/"+"Deeplift_full_ver2_fold"+str(i)+"_"+str(j)+".npy"
  tmpdf=np.load(fname)
  tmplis.append(tmpdf)
concdf=np.concatenate(tmplis,1)
alllist.append(concdf)

allgrad=np.concatenate(alllist,0)
fname=samplename+"/"+"Deeplift_full_ver2_all.npy"
np.save(fname,allgrad)


