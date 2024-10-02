from Bio import SeqIO
import numpy as np
import pickle
import sys

args = sys.argv
samplename=str(args[1])


maxlen=0
lennum=0

for record in SeqIO.parse(samplename+"/peaks_extend.fasta", 'fasta'):
    seq = record.seq
    if maxlen<len(seq):
        maxlen=len(seq)
    lennum=lennum+1

print(lennum)
print(maxlen)

df=np.zeros((lennum,maxlen,4))

cp=0
for record in SeqIO.parse(samplename+"/peaks_extend.fasta", 'fasta'):
    seq = record.seq

    for j in range(len(seq)):
        if(seq[j]=="A"):
            df[cp,j,0]=1
        if(seq[j]=="C"):
            df[cp,j,1]=1
        if(seq[j]=="G"):
            df[cp,j,2]=1
        if(seq[j]=="T"):
            df[cp,j,3]=1
    cp=cp+1

np.save(samplename+"/peak_extend_read.npy",df)

#f = open('peak_read.obj', 'wb')
#pickle.dump(read_lis, f)

#f = open("peak_read.obj","rb")
#read_lis = pickle.load(f)