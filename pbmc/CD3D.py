import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
args = sys.argv
samplename=str(args[1])

fname=samplename+"/Deeplift_full_ver2_all.npy"
enhamtx=np.load(fname)
CD3D=enhamtx[4854,:,:]
pd.DataFrame(CD3D).to_csv(samplename+"/CD3Dactivity.csv",sep="\t",header=False, index=False)