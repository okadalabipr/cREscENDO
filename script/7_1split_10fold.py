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
from sync_batchnorm import convert_model, DataParallelWithCallback
import sys
import math
import subprocess

args = sys.argv
samplename=str(args[1])
chanknum=10

for i in range(chanknum):
  command = ["python","7_2DNAembed.py",samplename,str(i)]
  subprocess.call(command)
  command = ["python","7_3RNAprediction_t2_raw_2.py",samplename,str(i)]
  subprocess.call(command)
  command = ["python","7_4RNAprediction_load_t2_raw_2.py",samplename,str(i)]
  subprocess.call(command)
  command = ["python","7_5RNAprediction_load_predrna_t2_raw_2.py",samplename,str(i)]
  subprocess.call(command)
  command = ["python","7_6DeepLift_full_split.py",samplename,str(1000),str(i)]
  subprocess.call(command)