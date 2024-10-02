import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch
from torch.utils.data import TensorDataset, DataLoader
from torch import nn
import torch.nn.functional as F
from torch import optim
import os
#import cv2
import datetime
import pickle
from tqdm import tqdm
from scipy import stats
import copy
import random
from sync_batchnorm import convert_model, DataParallelWithCallback
import sys
import math
import subprocess

args = sys.argv
samplename=str(args[1])
eachcell=int(float(args[2])) #1000
foldid=int(float(args[3]))

#samplename="/users/ken.murakami/workspace/pbmcnew/try1"
#eachcell=1000

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
  command = ["python","7_7DeepLift_full_ver2.py", samplename, str(cellnum), str(celltag), str(i), str(foldid)]
  subprocess.call(command)


# python 7_7DeepLift_full_ver2.py samplename 1000 celltag i foldid(kotei)