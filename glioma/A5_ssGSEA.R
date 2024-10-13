setwd("D:/scbasset_test/multiome_human_neuro/P-1694_S-1694/for_paper")
#rm(list =ls())

library(dplyr)
library(Seurat)
library(Matrix)


h5file_dir_path=""
working_dir_path=""


count<-Read10X(paste0(h5file_dir_path,"/datamtx"))

celllist<-read.csv(paste0(working_dir_path,"/celllist_scVI.csv"),header=F)
pairlist<-read.csv(paste0(working_dir_path,"/pair_300000.csv"),header=F)

count_s<-count$"Gene Expression"[,celllist$V1]
count_s<-count_s[pairlist$V1,]

gexm <- CreateSeuratObject(counts = count_s, project = "pbmc10k")
gexm <- NormalizeData(gexm)
gexm <- FindVariableFeatures(gexm, selection.method = "vst", nfeatures = 2000)
gexm <- ScaleData(gexm, verbose = FALSE)
gexm <- RunPCA(gexm, npcs = 30, verbose = FALSE)
gexm <- RunUMAP(gexm, reduction = "pca", dims = 1:30)

######

cellclus<-read.csv(paste0(working_dir_path,"/cellcluster.txt"),header=F)
cellleiden<-read.csv(paste0(working_dir_path,"/cellleiden.txt"),header=F)
gexm$subclus<-cellleiden$V1

gexm_sub=subset(gexm, subset = subclus>=0)
gexm_sub <- FindVariableFeatures(gexm_sub, selection.method = "vst", nfeatures = 7805)
gexm_sub <- ScaleData(gexm_sub, verbose = FALSE)
gexm_sub <- RunPCA(gexm_sub, npcs = 30, verbose = FALSE)
gexm_sub <- RunUMAP(gexm_sub, reduction = "pca", dims = 1:30)


gexm_sub@meta.data["subclus"]<-cellleiden[cellclus==0]


library(escape)
gene.sets3 <- getGeneSets(library = "H")

ES.seurat1 <- enrichIt(obj = gexm_sub, 
                       gene.sets = gene.sets3, 
                       groups = 1000, cores = 2)

write.table(ES.seurat1,paste0(working_dir_path,"/ssGSEA.csv"),sep=",",quote=F,row.names=T,col.names=T)


ES.seurat1<-read.csv(paste0(working_dir_path,"/ssGSEA.csv"))

wdiff=matrix(0, nrow = 50, ncol = 2)

ESo=apply(ES.seurat1[leidensub==0,],2,mean)-apply(ES.seurat1[leidensub==1,],2,mean)

leidensub=gexm_sub@meta.data$subclus

for (i in 1:50){
  tmptest=wilcox.test(ES.seurat1[leidensub==0,i], ES.seurat1[leidensub==1,i], paired = FALSE)
  wdiff[i,1]=tmptest$p.value
  wdiff[i,2]=ESo[i]
}

rownames(wdiff)=colnames(ES.seurat1)
write.table(wdiff,paste0(working_dir_path,"/wdiff.csv"),sep=",",quote=F,row.names=T,col.names=T)

