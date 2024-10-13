setwd("D:/scbasset_test/multiome_human_neuro/P-1694_S-1694/for_paper")
#rm(list =ls())

library(dplyr)
library(Seurat)


h5file_dir_path=""
working_dir_path=""


count<-Read10X(paste0(h5file_dir_path,"/datamtx"))

celllist<-read.csv(paste0(working_dir_path,"/celllist_scVI.csv"),header=F)
pairlist<-read.csv(paste0(working_dir_path,"/pair_300000.csv"),header=F)


count_s<-count$"Gene Expression"[,celllist$V1]
count_s<-count_s[pairlist$V1,]
count_g<-count$"Gene Expression"[,celllist$V1]

gexm <- CreateSeuratObject(counts = count_s, project = "pbmc10k")
gexm <- NormalizeData(gexm)
gexm <- FindVariableFeatures(gexm, selection.method = "vst", nfeatures = 2000)
gexm <- ScaleData(gexm, verbose = FALSE)
gexm <- RunPCA(gexm, npcs = 30, verbose = FALSE)
gexm <- RunUMAP(gexm, reduction = "pca", dims = 1:30)
gexm <- FindClusters(gexm, resolution = 0.1)

write.table(gexm@meta.data$seurat_clusters,paste0(working_dir_path,"/cellcluster.txt"),quote=F,sep="\t",col.names=F,row.names=F)


all.markersrna <- FindAllMarkers(object = gexm)

gene0=all.markersrna[(all.markersrna$cluster==0) & (all.markersrna$p_val_adj<0.05) & (all.markersrna$avg_log2FC>0.3),]
gene1=all.markersrna[(all.markersrna$cluster==1) & (all.markersrna$p_val_adj<0.05) & (all.markersrna$avg_log2FC>0.3),]
gene2=all.markersrna[(all.markersrna$cluster==2) & (all.markersrna$p_val_adj<0.05) & (all.markersrna$avg_log2FC>0.3),]
gene3=all.markersrna[(all.markersrna$cluster==3) & (all.markersrna$p_val_adj<0.05) & (all.markersrna$avg_log2FC>0.3),]

write.table(gene0$gene,paste0(working_dir_path,"/cluster0_deg_glioma.txt"),sep="\t",quote=F,row.names=F,col.names=F)
write.table(gene1$gene,paste0(working_dir_path,"/cluster1_deg_glioma.txt"),sep="\t",quote=F,row.names=F,col.names=F)
write.table(gene2$gene,paste0(working_dir_path,"/cluster2_deg_glioma.txt"),sep="\t",quote=F,row.names=F,col.names=F)
write.table(gene3$gene,paste0(working_dir_path,"/cluster3_deg_glioma.txt"),sep="\t",quote=F,row.names=F,col.names=F)


gexmg <- CreateSeuratObject(counts = count_g, project = "pbmc10k")
gexmg <- NormalizeData(gexmg,normalization.method="LogNormalize")
gexmg <- ScaleData(gexmg, verbose = FALSE)
log1pgs=(gexmg@assays$RNA@scale.data)
tmp=log1pgs["SOX2",]
write.table(tmp,paste0(working_dir_path,"/sox2exp.txt"),sep="\t",quote=F,row.names=F,col.names=F)
