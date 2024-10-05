setwd("D:/scbasset_test/pbmc_for_paper")
rm(list =ls())

library(dplyr)
library(Seurat)
library(Signac)
library(Matrix)
library(SeuratData)
library(SeuratDisk)
library(ggplot2)

count<-Read10X_h5("pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
celllist<-read.csv("celllist_scVI.csv",header=F)
count_s<-count$"Gene Expression"[,celllist$V1]


gexm <- CreateSeuratObject(counts = count_s, project = "pbmc10k",min.cells = 3)
gexm <- NormalizeData(gexm)
gexm <- FindVariableFeatures(gexm, selection.method = "vst", nfeatures = 2000)
gexm <- ScaleData(gexm, verbose = FALSE)
gexm <- RunPCA(gexm, npcs = 30, verbose = FALSE)
gexm <- RunUMAP(gexm, reduction = "pca", dims = 1:30)
gexm <- FindNeighbors(gexm, reduction = "pca", dims = 1:30)
gexm <- FindClusters(gexm, resolution = 0.1)
DimPlot(gexm, reduction = "umap",group.by="seurat_clusters",label = T)


g<-DimPlot(gexm, reduction = "umap",group.by="seurat_clusters")
print(g)
ggsave("pbmc_cluster.pdf",width = 8, height = 6)


CD3D<-read.csv("CD3Dactivity.csv",sep="\t",header=F)
gexm@meta.data$CD3D_promoter<-CD3D[,1]
gexm@meta.data$CD3D_chr11_118334494_118335388<-CD3D[,27]
gexm@meta.data$CD3D_chr11_118241249_118242001<-CD3D[,19]



degma <- FindMarkers(gexm, ident.1 = 5, ident.2 = 0)

png("CD3D_rna.png", width = 700, height = 500)
g<-FeaturePlot(gexm, features="CD3D")&labs(color = "Observed\nExpression")&theme(axis.text=element_text(size=20),
       axis.title = element_text(size=30),
       legend.text=element_text(size=15),
       legend.title=element_text(size=20),title=element_text(size=30),panel.background = element_rect(fill = "white", colour = "black", size = 1.2))
print(g)
dev.off()

ggsave("CD3D_rna.pdf",width = 8, height = 6)

png("CD3D_promoter.png", width = 700, height = 500)
FeaturePlot(gexm, features="CD3D_promoter",min.cutoff = 0,max.cutoff = 0.003)&labs(color = "Predicted\nActivity")&theme(axis.text=element_text(size=20),
                                                                                    axis.title = element_text(size=30),
                                                                                    legend.text=element_text(size=15),
                                                                                    legend.title=element_text(size=20),title=element_text(size=30),panel.background = element_rect(fill = "white", colour = "black", size = 1.2))
dev.off()

g<-FeaturePlot(gexm, features="CD3D_promoter",min.cutoff = 0,max.cutoff = 0.002)&labs(color = "Predicted\nActivity")&theme(axis.text=element_text(size=20),
                                                                                                                        axis.title = element_text(size=30),
                                                                                                                        legend.text=element_text(size=15),
                                                                                                                        legend.title=element_text(size=20),title=element_text(size=30),panel.background = element_rect(fill = "white", colour = "black", size = 1.2))



ggsave("CD3D_promoter.pdf",width = 8, height = 6)


png("CD3D_chr11_118334494_118335388.png", width = 700, height = 500)
FeaturePlot(gexm, features="CD3D_chr11_118334494_118335388",min.cutoff = 0.05,max.cutoff = 0.125)&labs(color = "Predicted\nActivity")&theme(axis.text=element_text(size=20),
                                                                                                                                        axis.title = element_text(size=30),
                                                                                                                                        legend.text=element_text(size=15),
                                                                                                                                        legend.title=element_text(size=20),title=element_text(size=30),panel.background = element_rect(fill = "white", colour = "black", size = 1.2))
dev.off()


g<-FeaturePlot(gexm, features="CD3D_chr11_118334494_118335388",min.cutoff = 0.05,max.cutoff = 0.125)&labs(color = "Predicted\nActivity")&theme(axis.text=element_text(size=20),
                                                                                                                                            axis.title = element_text(size=30),
                                                                                                                                            legend.text=element_text(size=15),
                                                                                                                                            legend.title=element_text(size=20),title=element_text(size=30),panel.background = element_rect(fill = "white", colour = "black", size = 1.2))

ggsave("CD3D_chr11_118334494_118335388.pdf",width = 8, height = 6)



png("CD3D_chr11_118241249_118242001.png", width = 700, height = 500)
FeaturePlot(gexm, features="CD3D_chr11_118241249_118242001",min.cutoff = 0.00,max.cutoff = 0.005)&labs(color = "Predicted\nActivity")&theme(axis.text=element_text(size=20),
                                                                                                                                       axis.title = element_text(size=30),
                                                                                                                                       legend.text=element_text(size=15),
                                                                                                                                       legend.title=element_text(size=20),title=element_text(size=30),panel.background = element_rect(fill = "white", colour = "black", size = 1.2))
dev.off()

g<-FeaturePlot(gexm, features="CD3D_chr11_118241249_118242001",min.cutoff = 0.0,max.cutoff = 0.005)&labs(color = "Predicted\nActivity")&theme(axis.text=element_text(size=20),
                                                                                                                                        axis.title = element_text(size=30),
                                                                                                                                        legend.text=element_text(size=15),
                                                                                                                                        legend.title=element_text(size=20),title=element_text(size=30),panel.background = element_rect(fill = "white", colour = "black", size = 1.2))


ggsave("CD3D_chr11_118241249_118242001.pdf",width = 8, height = 6)
