import numpy as np
import pandas as pd
import sys

args = sys.argv
samplename=str(args[1])

peaklist=pd.read_csv(samplename+"/peaks_extend.bed",sep="\t",header=None, names=("chr","start","end","center"))
peaklist["sequence_name"]=peaklist["chr"]+":"+peaklist["start"].astype(str)+"-"+peaklist["end"].astype(str)

fname="glioma/SOX2/fimo.tsv"
motiflist=pd.read_csv(fname,sep="\t")
#motiflist=motiflist.dropna(how='any')
motiflist=motiflist.sort_values(by=['p-value'], ascending=True)
motiflist=motiflist.drop_duplicates(subset='sequence_name')
pthv=motiflist["p-value"].iloc[5000]
motiflist=motiflist[motiflist["p-value"]<=pthv]
genename=motiflist["motif_alt_id"].iloc[0]
motiflist[genename]=True
tmpdf=motiflist[["sequence_name",genename]]
peaklist=pd.merge(peaklist, tmpdf, on='sequence_name', how='left',copy=False)
#peaklist=peaklist.join(tmpdf, on='sequence_name')
peaklist=peaklist.fillna(False)

fname=samplename+"/ZEB1/fimo.tsv"
motiflist=pd.read_csv(fname,sep="\t")
#motiflist=motiflist.dropna(how="any")
motiflist=motiflist.sort_values(by=["p-value"], ascending=True)
motiflist=motiflist.drop_duplicates(subset="sequence_name")
pthv=motiflist["p-value"].iloc[5000]
motiflist=motiflist[motiflist["p-value"]<=pthv]
genename=motiflist["motif_alt_id"].iloc[0]
motiflist[genename]=True
tmpdf=motiflist[["sequence_name",genename]]
peaklist=pd.merge(peaklist, tmpdf, on="sequence_name", how="left",copy=False)
#peaklist=peaklist.join(tmpdf, on="sequence_name")
peaklist=peaklist.fillna(False)

fname=samplename+"/CREM/fimo.tsv"
motiflist=pd.read_csv(fname,sep="\t")
#motiflist=motiflist.dropna(how="any")
motiflist=motiflist.sort_values(by=["p-value"], ascending=True)
motiflist=motiflist.drop_duplicates(subset="sequence_name")
pthv=motiflist["p-value"].iloc[5000]
motiflist=motiflist[motiflist["p-value"]<=pthv]
genename=motiflist["motif_alt_id"].iloc[0]
motiflist[genename]=True
tmpdf=motiflist[["sequence_name",genename]]
peaklist=pd.merge(peaklist, tmpdf, on="sequence_name", how="left",copy=False)
#peaklist=peaklist.join(tmpdf, on="sequence_name")
peaklist=peaklist.fillna(False)

peaklist.to_csv(samplename+"/motif_matrix.csv")


#####################


motifmtx=pd.read_csv(samplename+"/motif_matrix.csv")
enhamtxr=np.load(samplename+"/Deeplift_full_ver2.npy") #B,C,L

motifna=motifmtx.iloc[:,6:].to_numpy()

promenhatag=np.load(samplename+"/promenhatag.npy")

enhamtx=enhamtxr.copy()
enhatmp=enhamtxr.transpose(0,2,1)[(promenhatag==4)]
enhatmp=(enhatmp-enhatmp.mean(axis=0))/enhatmp.std(axis=0)
enhamtx.transpose(0,2,1)[(promenhatag==4)]=enhatmp


pairlist=pd.read_csv(samplename+"/pair_300000.csv",header=None,names=("gene","start","end","pos"))
pairlist=pairlist[["start","end","pos"]].to_numpy()

max_len=int((pairlist[:,1]-pairlist[:,0]).max()+1)


enhause=np.zeros((pairlist.shape[0],max_len))
enhamotif=np.zeros((pairlist.shape[0],max_len,motifna.shape[1]))

for j in range(pairlist.shape[0]):
  peakstart=int(pairlist[j,0])
  peakend=int(pairlist[j,1])
  peaknum_gene=peakend-peakstart+1
  enhause[j,0:peaknum_gene]=1
  enhamotif[j,0:peaknum_gene,:]=motifna[peakstart:(peakend+1),:]



enhamotif_mix=(enhamotif.transpose(2,0,1)*((promenhatag==4))).transpose(1,2,0)

np.save(samplename+"/enhamotif.npy",enhamotif)
np.save(samplename+"/enhamotif_mix.npy",enhamotif_mix)

enhamotif_re=enhamotif_mix.reshape(enhamotif_mix.shape[0]*enhamotif_mix.shape[1],enhamotif_mix.shape[2]) #B,L,T to B*L,T
enhamotif_re_norm=enhamotif_re/enhamotif_re.sum(axis=0)

enhamtx_re=enhamtx.transpose(0,2,1)
enhamtx_re=enhamtx_re.reshape(enhamtx_re.shape[0]*enhamtx_re.shape[1],enhamtx_re.shape[2])
enhamtx_re=enhamtx_re.transpose(1,0) #C,B*L

motifact=np.matmul(enhamtx_re, enhamotif_re_norm) #C,T
cellnamelist=pd.read_csv(samplename+"/celllist_scVI.csv",header=None)
motifact_df=pd.DataFrame(motifact,columns=motifmtx.columns[6:],index=np.array(cellnamelist[0]))
motifact_df.to_csv(samplename+"/motif_cell_activity.csv")


