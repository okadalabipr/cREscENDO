import scanpy as sc
import numpy as np
import pandas as pd
import sys

args = sys.argv
samplename=str(args[1])
inputname=str(args[2])

adata = sc.read_10x_h5(inputname,gex_only=False)
gex=adata[:,adata.var.feature_types=="Gene Expression"]
atac=adata[:,adata.var.feature_types=="Peaks"]
thres = int(adata.shape[0]*0.05)
del adata

sc.pp.filter_genes(gex, min_cells=thres)
sc.pp.filter_cells(gex, min_genes=500)
sc.pp.normalize_total(gex, target_sum=1e4)
sc.pp.log1p(gex)
gex.raw = gex
gex.var_names_make_unique()

genelist=pd.read_csv(samplename+"/pair_promoter.csv",sep=",",header=None)
RNAmatrix_sep=gex.T[genelist[0]].X.toarray()
RNAmatrix_sep=pd.DataFrame(RNAmatrix_sep)
RNAmatrix_sep.to_csv(samplename+"/log1praw.csv", header=False, index=False)

with open(samplename+"/celllist_scVI.csv", "w") as f:
    f.writelines([d+"\n" for d in gex.obs.index.tolist()])

traincellnum=int(RNAmatrix_sep.shape[1]*0.5)
np.save(samplename+"/traincellnum.npy",traincellnum)