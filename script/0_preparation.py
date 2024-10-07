import scanpy as sc
import numpy as np
import pandas as pd
import numpy as np
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

cellnum=gex.shape[0]
train_cellnum=int(cellnum*0.9)
gex_train=gex[0:train_cellnum,:]

gex_train.write_h5ad(samplename+"/gexm_traincell.h5ad")
del gex_train

gexcelllis=gex.obs.index
gex.write_h5ad(samplename+"/gexm_allcell.h5ad")
del gex


sc.pp.filter_genes(atac, min_cells=thres)
atac=atac[gexcelllis,:]
atac.write_h5ad(samplename+"/gexm_allcell_atac.h5ad")

peakbed=np.empty((atac.var.index.shape[0],3), dtype=object)
for i in range(atac.var.index.shape[0]):
	tmplist=str(atac.var.index[i]).split(":")
	peakbed[i,0]=tmplist[0]
	peakbed[i,1:3]=tmplist[1].split("-")
	
np.savetxt(samplename+"/peaks_before.bed", peakbed, delimiter='\t', fmt='%s')

######################################
# scBasset preparation

bed_file=samplename+"/peaks_before.bed"
peak = pd.read_csv(bed_file, sep='\t', names=['chr','start','end'])

atac.var['chr'] = peak['chr'].values
atac.var['start'] = peak['start'].values
atac.var['end'] = peak['end'].values

atac.raw=atac
atac.__dict__['_raw'].__dict__['_var'] = atac.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
del atac.raw

chrs = ['chr'+str(i) for i in range(1,23)] + ['chrX', 'chrY']
atac_s = atac[:, atac.var['chr'].isin(chrs)]
atac_s.write(samplename+"/atac_ad.h5ad")

m = atac_s.X.tocoo().transpose().toarray()
np.save(samplename+"/atac_count.npy",m)

peakbed_s=peakbed[(atac.var['chr'].isin(chrs)),:]
np.savetxt(samplename+"/peaks.bed", peakbed_s, delimiter='\t', fmt='%s')