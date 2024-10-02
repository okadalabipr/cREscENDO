import numpy as np
import pandas as pd
import os
import sys

args = sys.argv
samplename=str(args[1])
#samplename="/users/ken.murakami/workspace/pbmcnew/try5"

filelis=os.listdir(samplename)

alllist=[]
for i in range(10):
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


##########
import numpy as np
import pandas as pd
import os
import sys

#samplename="/users/ken.murakami/workspace/pbmcnew/try5"
fname=samplename+"/"+"Deeplift_full_ver2_all.npy"
allgrad=np.load(fname)

pairlistall=pd.read_csv(samplename+"/pair_300000.csv",header=None)
pairlistall=pairlistall[[0]].to_numpy()

pairlist=pd.read_csv(samplename+"/pair_300000.csv",header=None)
pairlist=pairlist[[0]].to_numpy()

pairlist_list=[]
for i in range(10):
  fname=samplename+"/pair_300000_test_fold"+str(i)+".csv"
  pairlist_tmp=pd.read_csv(fname,header=None)
  pairlist_tmp=pairlist_tmp[[0]].to_numpy()
  pairlist_list.append(pairlist_tmp)

pairlist_tgt=np.concatenate(pairlist_list,0)
allgrad_s=allgrad[pd.DataFrame(pairlist_tgt).drop_duplicates().index,:,:]


allgrad_s_max=allgrad_s.max(axis=1)
allgrad_s_max_pd=pd.DataFrame(allgrad_s_max.transpose(1,0),columns= pd.DataFrame(pairlist_tgt).drop_duplicates().to_numpy()[:,0])
allgrad_ssep_max=allgrad_s_max_pd[pairlist[:,0]].T.to_numpy()

fname=samplename+"/"+"allgrad_ssep_max.npy"
np.save(fname,allgrad_ssep_max)

tag=np.array(range(allgrad_s_max_pd.shape[1]))
tag=tag[np.newaxis,:]
tag=pd.DataFrame(tag,columns= pd.DataFrame(pairlist_tgt).drop_duplicates().to_numpy()[:,0])
tag=tag[pairlist[:,0]].T.to_numpy()
tag=tag[:,0]
allgrad_so=allgrad_s[tag]

fname=samplename+"/"+"Deeplift_full_ver2_all.npy"
np.save(fname,allgrad_so)