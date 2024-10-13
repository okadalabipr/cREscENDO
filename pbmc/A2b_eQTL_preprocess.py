import numpy as np
import pandas as pd
import sys

args = sys.argv
samplename=str(args[1])

eqtl_tmp=pd.read_csv(samplename+"/CAVIAR_Results_v8_GTEx_LD_ALL_NOCUTOFF.txt",sep="\t")
eqtl_tmp_sep=eqtl_tmp[eqtl_tmp["TISSUE"]=="Whole_Blood"]
eqtl_tmp_sep.to_csv(samplename+"/CAVIAR_Results_v8_GTEx_LD_ALL_NOCUTOFF_Whole_Blood.txt",sep="\t")

peaklist=pd.read_csv(samplename+"/peaks.bed",sep="\t",header=None)
peaklist=peaklist[[1,2]].to_numpy()
pairlist=pd.read_csv(samplename+"/pair_300000.csv",header=None)
pairlist=pairlist[[1,2,3]].to_numpy()

pairlist_prom=pd.read_csv(samplename+"/pair_promoter.csv",header=None)
pairlist_prom=pairlist_prom[[1,2,3]].to_numpy()

peaklist_full=pd.read_csv(samplename+"/peaks_extend.bed",sep="\t",header=None)
peaklist_full=peaklist_full.to_numpy()

genelist=pd.read_csv(samplename+"/pair_300000.csv",header=None)
genelist=genelist[[0]].to_numpy()
genelist=genelist[:,0]

ensemblelist=pd.read_csv(samplename+"/ensemblelist.csv")
ensemblelist.columns = ["tag", "gene_id", "name"]

max_len=int((pairlist[:,1]-pairlist[:,0]).max()+1)

peakidmtx=np.zeros((genelist.shape[0],max_len))
for j in range(genelist.shape[0]):
    peakstart=int(pairlist[j,0].item())
    peakend=int(pairlist[j,1].item())
    prompos=int(pairlist_prom[j,0].item())
    enha=list(range(peakstart,(peakend+1)))
    enha.remove(prompos)
    peaknum_gene=peakend-peakstart+1
    peakidmtx[j,0]=prompos
    peakidmtx[j,1:peaknum_gene]=enha

eqtl_tmp_sep=pd.read_csv(samplename+"/CAVIAR_Results_v8_GTEx_LD_ALL_NOCUTOFF_Whole_Blood.txt",sep="\t",index_col=0)

tgtgene=eqtl_tmp_sep["GENE"].str.split(".",expand=True)
tgtgene.columns=["gene_id","tag"]
tgtgene_mg=pd.merge(tgtgene, ensemblelist, how="left", on = "gene_id")

eqtl_mdn=eqtl_tmp_sep
snps=eqtl_mdn[["CHROM","POS","Probability"]]

snps.index=tgtgene_mg.index
snps_mg = pd.concat([snps, tgtgene_mg], axis=1)
snps_mg=snps_mg.dropna(how="any")



snpdfa=np.array(snps_mg[["CHROM", "POS", "name", "Probability"]])
snpdfa[:,1]=snpdfa[:,1].astype(int)

eqtl_mtx_all_cav=np.zeros((genelist.shape[0],max_len))
for gn in range(genelist.shape[0]):
    print(gn)
    peakstart=int(pairlist[gn,0].item())
    peakend=int(pairlist[gn,1].item())
    peaknum=peakend-peakstart+1
    tmphic_full=snpdfa[snpdfa[:,2]==genelist[gn],:].copy()
    full_in=(tmphic_full[:,1].astype(int)>=peaklist_full[peakstart,1])&(tmphic_full[:,1].astype(int)<=peaklist_full[peakend,2])
    tmphic_full_s=tmphic_full[(full_in),:].copy()
    for p1 in range(peaknum):
        p1_id=peakidmtx[gn,p1].astype(int)
        p1_in=(tmphic_full_s[:,1].astype(int)>=peaklist_full[p1_id,1])&(tmphic_full_s[:,1].astype(int)<=peaklist_full[p1_id,2])
        maemuki=tmphic_full_s[(p1_in),:]
        if maemuki.shape[0]>0:
            eqtl_mtx_all_cav[gn,p1]=maemuki.shape[0]

np.save(samplename+"/eqtl_mtx_all_cav.npy",eqtl_mtx_all_cav)



snps_mg_sep=snps_mg[snps_mg["Probability"]>0.5]


snpdfa=np.array(snps_mg_sep[["CHROM", "POS", "name", "Probability"]])
snpdfa[:,1]=snpdfa[:,1].astype(int)

eqtl_mtx_sig=np.zeros((genelist.shape[0],max_len))
for gn in range(genelist.shape[0]):
    print(gn)
    peakstart=int(pairlist[gn,0].item())
    peakend=int(pairlist[gn,1].item())
    peaknum=peakend-peakstart+1
    tmphic_full=snpdfa[snpdfa[:,2]==genelist[gn],:].copy()
    full_in=(tmphic_full[:,1].astype(int)>=peaklist_full[peakstart,1])&(tmphic_full[:,1].astype(int)<=peaklist_full[peakend,2])
    tmphic_full_s=tmphic_full[(full_in),:].copy()
    for p1 in range(peaknum):
        p1_id=peakidmtx[gn,p1].astype(int)
        p1_in=(tmphic_full_s[:,1].astype(int)>=peaklist_full[p1_id,1])&(tmphic_full_s[:,1].astype(int)<=peaklist_full[p1_id,2])
        maemuki=tmphic_full_s[(p1_in),:]
        if maemuki.shape[0]>0:
            eqtl_mtx_sig[gn,p1]=maemuki.shape[0]

np.save(samplename+"/eqtl_mtx_sig05.npy",eqtl_mtx_sig)