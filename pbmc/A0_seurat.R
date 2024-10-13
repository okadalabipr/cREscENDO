rm(list =ls())

library(dplyr)
library(Seurat)
#library(Matrix)
#library(SeuratData)
#library(SeuratDisk)

h5file_dir_path=""
working_dir_path=""

count<-Read10X_h5(paste0(h5file_dir_path,"pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)

celllist_raw<-read.csv(paste0(working_dir_path,"/celllist_scVI.csv"),header=F)
genelist<-read.csv(paste0(working_dir_path,"../../../papertest/pair_promoter.csv"),header=F)

gexm <- CreateSeuratObject(counts = count$`Gene Expression`,min.cells = 3, min.features = 200)
gexm <- NormalizeData(gexm)
gexm <- FindVariableFeatures(gexm, selection.method = "vst", nfeatures = 2000)
gexm <- ScaleData(gexm, verbose = FALSE)
gexm <- RunPCA(gexm, npcs = 30, verbose = FALSE)
gexm <- RunUMAP(gexm, reduction = "pca", dims = 1:30)
gexm <- FindNeighbors(gexm, reduction = "pca", dims = 1:30)
gexm <- FindClusters(gexm, resolution = 0.1)

DimPlot(gexm, reduction = "umap",group.by="seurat_clusters",label = T)


colnames(celllist_raw)<-"cellname"
gexmmeta=as.data.frame(cbind(rownames(gexm@meta.data),gexm@meta.data$seurat_clusters))
colnames(gexmmeta)<-c("cellname","cluster")
gexmmeta$cluster=as.integer(gexmmeta$cluster)-1

gexmmeta_merge=dplyr::left_join(celllist_raw , gexmmeta ,by="cellname")
write.table(gexmmeta_merge$cluster,paste0(working_dir_path,"cellcluster_all.txt"),sep="\t",quote=F,col.names=F,row.names = F)


############################

clusterall<-read.csv(paste0(working_dir_path,"cellcluster_all.txt"),header=F)

rna=count$`Gene Expression`
tmp=rna[rownames(rna) %in% genelist$V1,]
rnas=tmp[,colnames(rna) %in% celllist_raw$cellname]

gexm <- CreateSeuratObject(counts = rnas,min.cells = 3, min.features = 200)
gexm <- NormalizeData(gexm)
gexm <- FindVariableFeatures(gexm, selection.method = "vst", nfeatures = 2000)
gexm <- ScaleData(gexm, verbose = FALSE)
gexm <- RunPCA(gexm, npcs = 30, verbose = FALSE)

gexm@meta.data$clusterall<-clusterall

Idents(gexm)<-gexm@meta.data$clusterall
all.markers <- FindAllMarkers(object = gexm)

clus0m=all.markers[all.markers$cluster=="0",]
clus1m=all.markers[all.markers$cluster=="1",]
clus2m=all.markers[all.markers$cluster=="2",]
clus3m=all.markers[all.markers$cluster=="3",]
clus4m=all.markers[all.markers$cluster=="4",]
clus5m=all.markers[all.markers$cluster=="5",]

clus0gene=clus0m[(clus0m$p_val<0.05)&(clus0m$avg_log2FC>0),"gene"]
clus1gene=clus1m[(clus1m$p_val<0.05)&(clus1m$avg_log2FC>0),"gene"]
clus2gene=clus2m[(clus2m$p_val<0.05)&(clus2m$avg_log2FC>0),"gene"]
clus3gene=clus3m[(clus3m$p_val<0.05)&(clus3m$avg_log2FC>0),"gene"]
clus4gene=clus4m[(clus4m$p_val<0.05)&(clus4m$avg_log2FC>0),"gene"]
clus5gene=clus5m[(clus5m$p_val<0.05)&(clus5m$avg_log2FC>0),"gene"]


write.table(clus0gene,paste0(working_dir_path,"clus0gene.txt"),sep="\t",quote=F,row.names=F,col.names=F)
write.table(clus1gene,paste0(working_dir_path,"clus1gene.txt"),sep="\t",quote=F,row.names=F,col.names=F)
write.table(clus2gene,paste0(working_dir_path,"clus2gene.txt"),sep="\t",quote=F,row.names=F,col.names=F)
write.table(clus3gene,paste0(working_dir_path,"clus3gene.txt"),sep="\t",quote=F,row.names=F,col.names=F)
write.table(clus4gene,paste0(working_dir_path,"clus4gene.txt"),sep="\t",quote=F,row.names=F,col.names=F)
write.table(clus5gene,paste0(working_dir_path,"clus5gene.txt"),sep="\t",quote=F,row.names=F,col.names=F)