import numpy as np
import pandas as pd
import torch
import sys
import math
import subprocess

args = sys.argv
samplename=str(args[1])
eachcell=int(float(args[2])) #1000
foldid=int(float(args[3]))


RNAembed=pd.read_csv(samplename+"/pca_50dim_rnaembed_raw.txt",sep=" ",header=None)
RNAembed=RNAembed.to_numpy()
RNAembed=torch.from_numpy(RNAembed)
RNAembed=RNAembed.to(torch.float32)

fullcell=RNAembed.shape[0] #11754
chanknum=math.ceil(fullcell/eachcell) #12

for i in range(chanknum):
  celltag=i*eachcell
  cellnum=eachcell
  if i==chanknum-1:
    cellnum=fullcell-celltag
  command = ["python","script/7_7DeepLift_full_ver2_10foldsplit.py", samplename, str(cellnum), str(celltag), str(i), str(foldid)]
  subprocess.call(command)

