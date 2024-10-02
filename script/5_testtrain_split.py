import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch
from torch.utils.data import TensorDataset, DataLoader
from torch import nn
import torch.nn.functional as F
from torch import optim
import os
import cv2
import datetime
import pickle
from tqdm import tqdm
from scipy import stats
import copy
import random
import sys

args = sys.argv
samplename=str(args[1])

ATACmatrix=np.load(samplename+"/atac_count.npy")
pairlist=pd.read_csv(samplename+"/pair_300000.csv",header=None,names=("gene","start","end","pos"))
peaknum=ATACmatrix.shape[0]

c_list=np.array(range(peaknum))
splitlist=np.percentile(c_list, q=[0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]).astype(int)

for i in range(10):
    peaknum_train_start=splitlist[i]
    peaknum_train_end=splitlist[i+1]
    gene_train=pairlist[(pairlist.end<peaknum_train_start)|(pairlist.start>peaknum_train_end)]
    gene_test=pairlist[(pairlist.end>=peaknum_train_start)&(pairlist.start<=peaknum_train_end)]

    gene_train.to_csv(samplename+"/pair_300000_train_fold"+str(i)+".csv", header=False, index=False)
    gene_test.to_csv(samplename+"/pair_300000_test_fold"+str(i)+".csv", header=False, index=False)
