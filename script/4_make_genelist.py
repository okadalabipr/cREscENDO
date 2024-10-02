import scanpy as sc
import numpy as np
import pandas as pd
import sys

args = sys.argv
samplename=str(args[1])

adata = sc.read_h5ad(samplename+"/gexm_allcell.h5ad")
genelist=adata.var.index

peaks=pd.read_table(samplename+"/peaks.bed",delimiter="\t",header=None,names=["chr","start","end"])
peaks["center"]=(peaks["start"]+peaks["end"])/2
peaks["center"]=peaks["center"].astype('int')

TSSlist=pd.read_table(samplename+"/TSS_list.txt",delimiter="\t",header=None,names=["ensembl_gene_id", "chromosome_name","start_position", "end_position","external_gene_name","strand"])
TSSlist=TSSlist[TSSlist["external_gene_name"].isin(genelist)]
TSSlist_f=TSSlist[TSSlist["strand"]==1]
TSSlist_r=TSSlist[TSSlist["strand"]==-1]

distance=500

Promoter_f=TSSlist_f[["chromosome_name","start_position", "end_position","external_gene_name"]].copy()
Promoter_r=TSSlist_r[["chromosome_name","start_position", "end_position","external_gene_name"]].copy()
Promoter_f["start_position"]=TSSlist_f["start_position"]-distance
Promoter_f["end_position"]=TSSlist_f["start_position"]+distance
Promoter_f["TSS_position"]=TSSlist_f["start_position"]
Promoter_r["start_position"]=TSSlist_r["end_position"]-distance
Promoter_r["end_position"]=TSSlist_r["end_position"]+distance
Promoter_r["TSS_position"]=TSSlist_r["end_position"]

Promoter_range= pd.concat([Promoter_f, Promoter_r])
Promoter_range["chromosome_name"]=str("chr")+Promoter_range["chromosome_name"]
Promoter_range=Promoter_range.drop_duplicates(subset="external_gene_name")

peak_pos = pd.DataFrame(index=range(len(Promoter_range)), columns=["gene_name", "peakstart","peakend", "TSS"])
for i in range(len(Promoter_range)):
    tmppeaks=peaks[(peaks["chr"]==Promoter_range["chromosome_name"].iloc[i]) & (peaks["center"]>Promoter_range["start_position"].iloc[i]) & (peaks["center"]<Promoter_range["end_position"].iloc[i])]                         
    if len(tmppeaks)!=0:
        minnum=abs(tmppeaks["center"]-Promoter_range["TSS_position"].iloc[i]).idxmin()
        maxnum=abs(tmppeaks["center"]-Promoter_range["TSS_position"].iloc[i]).idxmax()
        peak_pos["gene_name"].iloc[i]=Promoter_range["external_gene_name"].iloc[i]
        peak_pos["peakstart"].iloc[i]=minnum
        peak_pos["peakend"].iloc[i]=minnum
        peak_pos["TSS"].iloc[i]=Promoter_range["TSS_position"].iloc[i]

peak_pos = peak_pos.dropna()
peak_pos_select=peak_pos[peak_pos["gene_name"].isin(genelist)]
peak_pos_select.to_csv(samplename+"/pair_promoter.csv",header=False, index=False)

distance=300000

Promoter_f=TSSlist_f[["chromosome_name","start_position", "end_position","external_gene_name"]].copy()
Promoter_r=TSSlist_r[["chromosome_name","start_position", "end_position","external_gene_name"]].copy()
Promoter_f["start_position"]=TSSlist_f["start_position"]-distance
Promoter_f["end_position"]=TSSlist_f["start_position"]+distance
Promoter_f["TSS_position"]=TSSlist_f["start_position"]
Promoter_r["start_position"]=TSSlist_r["end_position"]-distance
Promoter_r["end_position"]=TSSlist_r["end_position"]+distance
Promoter_r["TSS_position"]=TSSlist_r["end_position"]

Promoter_range= pd.concat([Promoter_f, Promoter_r])
Promoter_range["chromosome_name"]=str("chr")+Promoter_range["chromosome_name"]
Promoter_range=Promoter_range.drop_duplicates(subset="external_gene_name")

peak_pos_enhancer = pd.DataFrame(index=range(len(Promoter_range)), columns=["gene_name", "peakstart","peakend", "TSS"])
for i in range(len(Promoter_range)):
    tmppeaks=peaks[(peaks["chr"]==Promoter_range["chromosome_name"].iloc[i]) & (peaks["center"]>Promoter_range["start_position"].iloc[i]) & (peaks["center"]<Promoter_range["end_position"].iloc[i])]                         
    if len(tmppeaks)!=0:
        minnum=min(tmppeaks.index)
        maxnum=max(tmppeaks.index)
        peak_pos_enhancer["gene_name"].iloc[i]=Promoter_range["external_gene_name"].iloc[i]
        peak_pos_enhancer["peakstart"].iloc[i]=minnum
        peak_pos_enhancer["peakend"].iloc[i]=maxnum
        peak_pos_enhancer["TSS"].iloc[i]=Promoter_range["TSS_position"].iloc[i]

peak_pos_enhancer_select=peak_pos_enhancer[peak_pos_enhancer["gene_name"].isin(peak_pos_select["gene_name"])]
peak_pos_enhancer_select.to_csv(samplename+"/pair_300000.csv",header=False, index=False)