setwd("D:/scbasset_test/multiome_human_neuro/P-1694_S-1694/for_paper")
#rm(list =ls())

library(dplyr)
library(Seurat)
library(Signac)
library(Matrix)
library(SeuratData)
library(SeuratDisk)
library(ggplot2)
library(stringr)
library(ggrepel)


library(extrafont)
font_import(pattern="arial")
windowsFonts(sans="arial")
loadfonts(device="win")
names(windowsFonts())


count<-Read10X("../datamtx")

celllist<-read.csv("../celllist_scVI.csv",header=F)
pairlist<-read.csv("../pair_300000.csv",header=F)

cellclus<-read.csv("../cellcluster.txt",header=F)
tumorclus<-read.csv("cellleiden.txt",header=F)

tumorpeaklist<-read.csv("tumorpeakidx.txt",header=F)
tumorpeaklist$V1=tumorpeaklist$V1+1

emtgene<-read.csv("emtgene.txt",header=F)
nfkbgene<-read.csv("nfkbgene.txt",header=F)

tfact<-read.csv("motif_cell_activity.csv")

thr=as.integer(dim(X)[2]*0.05)

count_s<-count$"Gene Expression"[,celllist$V1]
count_s<-count_s[pairlist$V1,]
count_g<-count$"Gene Expression"[,celllist$V1]


gexm <- CreateSeuratObject(counts = count_s, project = "pbmc10k")
gexm <- NormalizeData(gexm)
log1praw=as.matrix(gexm@assays$RNA@data)
write.table(log1praw,"log1praw.csv",sep=",",quote=F,col.names=T,row.names=T)

gexm <- FindVariableFeatures(gexm, selection.method = "vst", nfeatures = 2000)
gexm <- ScaleData(gexm, verbose = FALSE)
gexm <- RunPCA(gexm, npcs = 30, verbose = FALSE)
gexm <- RunUMAP(gexm, reduction = "pca", dims = 1:30)


cellcluster<-read.csv("../cellcluster.txt",header=F)
gexm$cluster<-cellcluster$V1
Idents(object = gexm) <- "cluster"
all.markersrna <- FindAllMarkers(object = gexm)

gene0=all.markersrna[(all.markersrna$cluster==0) & (all.markersrna$p_val_adj<0.05) & (all.markersrna$avg_log2FC>0.3),]
gene1=all.markersrna[(all.markersrna$cluster==1) & (all.markersrna$p_val_adj<0.05) & (all.markersrna$avg_log2FC>0.3),]
gene2=all.markersrna[(all.markersrna$cluster==2) & (all.markersrna$p_val_adj<0.05) & (all.markersrna$avg_log2FC>0.3),]
gene3=all.markersrna[(all.markersrna$cluster==3) & (all.markersrna$p_val_adj<0.05) & (all.markersrna$avg_log2FC>0.3),]

write.table(gene0$gene,"cluster0_deg_glioma.txt",sep="\t",quote=F,row.names=F,col.names=F)
write.table(gene1$gene,"cluster1_deg_glioma.txt",sep="\t",quote=F,row.names=F,col.names=F)
write.table(gene2$gene,"cluster2_deg_glioma.txt",sep="\t",quote=F,row.names=F,col.names=F)
write.table(gene3$gene,"cluster3_deg_glioma.txt",sep="\t",quote=F,row.names=F,col.names=F)


cellleiden<-read.csv("cellleiden.txt",header=F)
gexm$subclus<-cellleiden$V1

gexm_sub=subset(gexm, subset = subclus>=0)
gexm_sub <- FindVariableFeatures(gexm_sub, selection.method = "vst", nfeatures = 7805)
gexm_sub <- ScaleData(gexm_sub, verbose = FALSE)
gexm_sub <- RunPCA(gexm_sub, npcs = 30, verbose = FALSE)
gexm_sub <- RunUMAP(gexm_sub, reduction = "pca", dims = 1:30)

Idents(object = gexm_sub) <- "subclus"
all.markers <- FindAllMarkers(object = gexm_sub,logfc.threshold = 0.05)
deg=all.markers[(all.markers$cluster==1)&(all.markers$p_val_adj<0.05),]

posdeg=deg[deg$avg_log2FC>0,]
negdeg=deg[deg$avg_log2FC<0,]
write.table(posdeg$gene,"posdeggeneclus1.txt",sep="\t",quote=F,col.names=F,row.names=F)
write.table(negdeg$gene,"negdeggeneclus1.txt",sep="\t",quote=F,col.names=F,row.names=F)

cellleiden<-read.csv("cellleiden.txt",header=F)
gexm_sub@meta.data["subclus"]<-cellleiden[cellcluster==0]

gexmg <- CreateSeuratObject(counts = count_g, project = "pbmc10k")
gexmg <- NormalizeData(gexmg,normalization.method="LogNormalize")
gexmg <- ScaleData(gexmg, verbose = FALSE)

gexmg@meta.data["subclus"]<-cellleiden
gexmg@meta.data["cellclus"]<-cellclus

log1pg=(gexmg@assays$RNA@data)
log1pgs=(gexmg@assays$RNA@scale.data)



library(escape)
gene.sets3 <- getGeneSets(library = "H")

ES.seurat1 <- enrichIt(obj = gexm_sub, 
                       gene.sets = gene.sets3, 
                       groups = 1000, cores = 2)

write.table(ES.seurat1,"ssGSEA.csv",sep=",",quote=F,row.names=T,col.names=T)

ES.seurat1<-read.csv("ssGSEA.csv")

wdiff=matrix(0, nrow = 50, ncol = 2)

ESo=apply(ES.seurat1[leidensub==0,],2,mean)-apply(ES.seurat1[leidensub==1,],2,mean)

leidensub=gexm_sub@meta.data$subclus

for (i in 1:50){
  tmptest=wilcox.test(ES.seurat1[leidensub==0,i], ES.seurat1[leidensub==1,i], paired = FALSE)
  wdiff[i,1]=tmptest$p.value
  wdiff[i,2]=ESo[i]
}

rownames(wdiff)=colnames(ES.seurat1)
write.table(wdiff,"wdiff.csv",sep=",",quote=F,row.names=T,col.names=T)


hdiff=matrix(0, nrow = 50, ncol = 1)

for (i in 1:50){
  tmpgene=as.data.frame(msigdbr_list_h[i])
  tmpgene=tmpgene[,1]
  tmp=log1pgs[rownames(log1pgs) %in% tmpgene,]
  tmp20=apply(tmp[,cellleiden$V1==0],1,mean)
  tmp21=apply(tmp[,cellleiden$V1==1],1,mean)
  tmp4=tmp20-tmp21
  hdiff[i,]=mean(tmp4[is.finite(tmp4)])
}


msigdbr_list_h0<-msigdbr_list_h[hdiff>0]
msigdbr_list_h1<-msigdbr_list_h[hdiff<(-0.07)]

tmp=log1pgs["SOX2",]
write.table(tmp,"sox2exp.txt",sep="\t",quote=F,row.names=F,col.names=F)


Idents(object = gexm_sub) <- "cluster"
all.markers <- FindAllMarkers(object = gexm_sub,assay="RNA",only.pos=TRUE)

gene0=all.markers[(all.markers$cluster==0) & (all.markers$p_val_adj<0.05) & (all.markers$avg_log2FC>0.3),]
gene1=all.markers[(all.markers$cluster==1) & (all.markers$p_val_adj<0.05) & (all.markers$avg_log2FC>0.3),]
gene2=all.markers[(all.markers$cluster==2) & (all.markers$p_val_adj<0.05) & (all.markers$avg_log2FC>0.3),]
gene3=all.markers[(all.markers$cluster==3) & (all.markers$p_val_adj<0.05) & (all.markers$avg_log2FC>0.3),]

genelist<-rownames(gexm_sub)

genelist_sub<-genelist[!genelist %in% gene0$gene]
genelist_sub<-genelist_sub[!genelist_sub %in% gene1$gene]
genelist_sub<-genelist_sub[!genelist_sub %in% gene2$gene]
genelist_sub<-genelist_sub[!genelist_sub %in% gene3$gene]


write.table(genelist_sub,"non_deg_glioma.txt",sep="\t",quote=F,row.names=F,col.names=F)
write.table(gene0$gene,"cluster0_deg_glioma.txt",sep="\t",quote=F,row.names=F,col.names=F)
write.table(gene1$gene,"cluster1_deg_glioma.txt",sep="\t",quote=F,row.names=F,col.names=F)
write.table(gene2$gene,"cluster2_deg_glioma.txt",sep="\t",quote=F,row.names=F,col.names=F)
write.table(gene3$gene,"cluster3_deg_glioma.txt",sep="\t",quote=F,row.names=F,col.names=F)



