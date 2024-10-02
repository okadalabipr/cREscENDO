mkdir glioma/AP1
fimo -oc glioma/AP1 --text -motif MA0490.3 --thresh 0.001 JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt glioma/peaks_extend.fasta > glioma/AP1/fimo.tsv

fimo -oc glioma/AP1 --text -motif MA0490.3 --thresh 0.001 JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt glioma/peaks_extend.fasta > glioma/AP1/fimo.tsv

fimo -oc glioma/NFKB --text -motif MA0105.4 --thresh 0.001 JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt glioma/peaks_extend.fasta > glioma/NFKB/fimo.tsv

fimo -oc glioma/EGR1 --text -motif MA0162.5 --thresh 0.001 JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt glioma/peaks_extend.fasta > glioma/EGR1/fimo.tsv

mkdir glioma/NFIX
fimo -oc glioma/NFIX --text -motif MA0671.2 --thresh 0.001 JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt glioma/peaks_extend.fasta > glioma/NFIX/fimo.tsv

mkdir glioma/NFICTLX1
fimo -oc glioma/NFICTLX1 --text -motif MA0119.1 --thresh 0.001 JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt glioma/peaks_extend.fasta > glioma/NFICTLX1/fimo.tsv

mkdir glioma/ASCL1
fimo -oc glioma/ASCL1 --text -motif MA1100.3 --thresh 0.001 JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt glioma/peaks_extend.fasta > glioma/ASCL1/fimo.tsv

mkdir glioma/NFIB
fimo -oc glioma/NFIB --text -motif MA1643.2 --thresh 0.001 JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt glioma/peaks_extend.fasta > glioma/NFIB/fimo.tsv

mkdir glioma/CREM
fimo -oc glioma/CREM --text -motif MA0609.3 --thresh 0.001 JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt glioma/peaks_extend.fasta > glioma/CREM/fimo.tsv




#############3




import numpy as np
import pandas as pd

peaklist=pd.read_csv('glioma/peaks_extend.bed',sep="\t",header=None, names=("chr","start","end","center"))
peaklist["sequence_name"]=peaklist["chr"]+":"+peaklist["start"].astype(str)+"-"+peaklist["end"].astype(str)



fname="glioma/NFIC/fimo.tsv"
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

fname="glioma/SMAD3/fimo.tsv"
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

fname="glioma/ZEB1/fimo.tsv"
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


fname="glioma/AP1/fimo.tsv"
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


peaklist["SMAD3_AP1"]=(peaklist["SMAD3"]*peaklist["JUNB"])



fname="glioma/NFKB/fimo.tsv"
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

fname="glioma/EGR1/fimo.tsv"
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


fname="glioma/NFIX/fimo.tsv"
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

fname="glioma/NFICTLX1/fimo.tsv"
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

peaklist["SMAD3_ZEB1"]=(peaklist["SMAD3"]*peaklist["ZEB1"])
peaklist["SOX2_ZEB1"]=(peaklist["SOX2"]*peaklist["ZEB1"])
peaklist["SOX2_AP1"]=(peaklist["SOX2"]*peaklist["JUNB"])

fname="glioma/ASCL1/fimo.tsv"
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

fname="glioma/NFIB/fimo.tsv"
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


fname="glioma/CREM/fimo.tsv"
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

peaklist.to_csv('glioma/motif_matrix.csv')


#####################


import numpy as np
import pandas as pd

motifmtx=pd.read_csv('glioma/motif_matrix.csv')
enhamtxr=np.load("glioma/Deeplift_full_ver2_mergenorm.npy") #B,C,L

motifna=motifmtx.iloc[:,6:].to_numpy()

promenhatag=np.load("glioma/promenhatag.npy")

enhamtx=enhamtxr.copy()
enhatmp=enhamtxr.transpose(0,2,1)[(promenhatag==4)]
enhatmp=(enhatmp-enhatmp.mean(axis=0))/enhatmp.std(axis=0)
enhamtx.transpose(0,2,1)[(promenhatag==4)]=enhatmp
#normfac1=(np.abs(enhamtxr)*(enhamtxr>0)).sum(axis=(0,2))/(enhamtxr>0).sum(axis=(0,2))
#normfac2=(np.abs(enhamtxr)*(enhamtxr<0)).sum(axis=(0,2))/(enhamtxr<0).sum(axis=(0,2))
#enhamtx=((enhamtxr*(enhamtxr>0)).transpose(0,2,1)/normfac1).transpose(0,2,1)+((enhamtxr*(enhamtxr<0)).transpose(0,2,1)/normfac2).transpose(0,2,1)

pairlist=pd.read_csv('glioma/pair_300000.csv',header=None,names=("gene","start","end","pos"))
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
#enhamotif_mix=(enhamotif.transpose(2,0,1)).transpose(1,2,0)
enhamotif_mixprom=(enhamotif.transpose(2,0,1)*((promenhatag==3))).transpose(1,2,0)

np.save("glioma/enhamotif.npy",enhamotif)
np.save("glioma/enhamotif_mix.npy",enhamotif_mix)

enhamotif_re=enhamotif_mix.reshape(enhamotif_mix.shape[0]*enhamotif_mix.shape[1],enhamotif_mix.shape[2]) #B,L,T to B*L,T
enhamotif_re_norm=enhamotif_re/enhamotif_re.sum(axis=0)

enhamotif_reprom=enhamotif_mixprom.reshape(enhamotif_mixprom.shape[0]*enhamotif_mixprom.shape[1],enhamotif_mixprom.shape[2]) #B,L,T to B*L,T
enhamotif_reprom_norm=enhamotif_reprom/enhamotif_reprom.sum(axis=0)

enhamtx_re=enhamtx.transpose(0,2,1)
enhamtx_re=enhamtx_re.reshape(enhamtx_re.shape[0]*enhamtx_re.shape[1],enhamtx_re.shape[2])
enhamtx_re=enhamtx_re.transpose(1,0) #C,B*L

motifact=np.matmul(enhamtx_re, enhamotif_re_norm) #C,T
motifactprom=np.matmul(enhamtx_re, enhamotif_reprom_norm) #C,T

#motifact_norm=motifact/enhause.sum()

cellnamelist=pd.read_csv('glioma/celllist_scVI.csv',header=None)

motifact_df=pd.DataFrame(motifact,columns=motifmtx.columns[6:],index=np.array(cellnamelist[0]))
motifact_dfprom=pd.DataFrame(motifactprom,columns=motifmtx.columns[6:],index=np.array(cellnamelist[0]))

motifact_df.to_csv('glioma/motif_cell_activity.csv')
motifact_dfprom.to_csv('glioma/motif_cell_activity_prom.csv')


