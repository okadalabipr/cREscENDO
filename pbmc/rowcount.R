rm(list =ls())

library(dplyr)
library(Seurat)
library(Signac)
library(Matrix)
library(SeuratData)
library(SeuratDisk)
library(ggplot2)

count<-Read10X_h5("pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
dim(count$`Gene Expression`)
dim(count$`Peaks`)

celllist_raw<-read.csv("celllist_scVI.csv",header=F)
genelist<-read.csv("pair_promoter.csv",header=F)

gexm <- CreateSeuratObject(counts = count$`Gene Expression`, project = "pbmc10k",min.cells = 3, min.features = 200)
gexm <- NormalizeData(gexm)
gexm <- FindVariableFeatures(gexm, selection.method = "vst", nfeatures = 2000)
gexm <- ScaleData(gexm, verbose = FALSE)
gexm <- RunPCA(gexm, npcs = 30, verbose = FALSE)
gexm <- RunUMAP(gexm, reduction = "pca", dims = 1:30)
gexm <- FindNeighbors(gexm, reduction = "pca", dims = 1:30)
gexm <- FindClusters(gexm, resolution = 0.1)
DimPlot(gexm, reduction = "umap",group.by="seurat_clusters",label = T)

gexm@meta.data$seurat_clusters

colnames(celllist_raw)<-"cellname"
gexmmeta=as.data.frame(cbind(rownames(gexm@meta.data),gexm@meta.data$seurat_clusters))
colnames(gexmmeta)<-c("cellname","cluster")

gexmmeta_merge=dplyr::left_join(celllist_raw , gexmmeta ,by="cellname")
write.table(gexmmeta_merge$cluster,"cellcluster_all.txt",sep="\t",quote=F,col.names=F,row.names = F)


############################

clusterall<-read.csv("cellcluster_all.txt",header=F)

rna=count$`Gene Expression`
tmp=rna[rownames(rna) %in% genelist$V1,]
rnas=tmp[,colnames(rna) %in% celllist_raw$V1]

gexm <- CreateSeuratObject(counts = rnas, project = "pbmc10k",min.cells = 3, min.features = 200)
gexm <- NormalizeData(gexm)
gexm <- FindVariableFeatures(gexm, selection.method = "vst", nfeatures = 2000)
gexm <- ScaleData(gexm, verbose = FALSE)
gexm <- RunPCA(gexm, npcs = 30, verbose = FALSE)

gexm@meta.data$clusterall<-clusterall

Idents(gexm)<-gexm@meta.data$clusterall
all.markers <- FindAllMarkers(object = gexm)

clus1m=all.markers[all.markers$cluster=="1",]
clus2m=all.markers[all.markers$cluster=="2",]
clus3m=all.markers[all.markers$cluster=="3",]
clus4m=all.markers[all.markers$cluster=="4",]
clus5m=all.markers[all.markers$cluster=="5",]
clus6m=all.markers[all.markers$cluster=="6",]

clus1gene=clus1m[(clus1m$p_val<0.05)&(clus1m$avg_log2FC>0),"gene"]
clus2gene=clus2m[(clus2m$p_val<0.05)&(clus2m$avg_log2FC>0),"gene"]
clus3gene=clus3m[(clus3m$p_val<0.05)&(clus3m$avg_log2FC>0),"gene"]
clus4gene=clus4m[(clus4m$p_val<0.05)&(clus4m$avg_log2FC>0),"gene"]
clus5gene=clus5m[(clus5m$p_val<0.05)&(clus5m$avg_log2FC>0),"gene"]
clus6gene=clus6m[(clus6m$p_val<0.05)&(clus6m$avg_log2FC>0),"gene"]


write.table(clus1gene,"clus1gene.txt",sep="\t",quote=F,row.names=F,col.names=F)
write.table(clus2gene,"clus2gene.txt",sep="\t",quote=F,row.names=F,col.names=F)

write.table(clus3gene,"clus3gene.txt",sep="\t",quote=F,row.names=F,col.names=F)
write.table(clus4gene,"clus4gene.txt",sep="\t",quote=F,row.names=F,col.names=F)

write.table(clus5gene,"clus5gene.txt",sep="\t",quote=F,row.names=F,col.names=F)
write.table(clus6gene,"clus6gene.txt",sep="\t",quote=F,row.names=F,col.names=F)