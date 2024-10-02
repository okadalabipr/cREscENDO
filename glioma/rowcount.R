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
hypogene<-read.csv("hypogene.txt",header=F)
tgfbgene<-read.csv("tgfbgene.txt",header=F)
krasgene<-read.csv("krasgene.txt",header=F)
nfkbgene<-read.csv("nfkbgene.txt",header=F)
andrgene<-read.csv("andrgene.txt",header=F)
krasdngene<-read.csv("krasdngene.txt",header=F)

tfact<-read.csv("motif_cell_activity.csv")


count$"Gene Expression"["NFIC",]

dim(count$Peaks)
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

FeaturePlot(gexm, features="PCNA")

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
DimPlot(object = gexm_sub,label = T)

Idents(object = gexm_sub) <- "subclus"
all.markers <- FindAllMarkers(object = gexm_sub,logfc.threshold = 0.05)

deg=all.markers[(all.markers$cluster==1)&(all.markers$p_val_adj<0.05),]

#deg2 <- FindMarkers(gexm_sub, ident.1 = 1, ident.2 = 0)
posdeg=deg[deg$avg_log2FC>0,]
negdeg=deg[deg$avg_log2FC<0,]

write.table(posdeg$gene,"posdeggeneclus1.txt",sep="\t",quote=F,col.names=F,row.names=F)
write.table(negdeg$gene,"negdeggeneclus1.txt",sep="\t",quote=F,col.names=F,row.names=F)

selectgeneact<-gexm_sub@assays$RNA@data

gexm_sub@meta.data["emt_acitivty"]<-apply(selectgeneact[rownames(selectgeneact) %in% emtgene$V1,],2,mean)
gexm_sub@meta.data["hypo_acitivty"]<-apply(selectgeneact[rownames(selectgeneact) %in% hypogene$V1,],2,mean)

cellleiden<-read.csv("cellleiden.txt",header=F)
gexm_sub@meta.data["subclus"]<-cellleiden[cellcluster==0]

DimPlot(object = gexm_sub, group.by="subclus",label = T)
FeaturePlot(gexm_sub, features="emt_acitivty")
FeaturePlot(gexm_sub, features="hypo_acitivty")

FeaturePlot(gexm_sub, features="NFIX")
FeaturePlot(gexm_sub, features="CDH1")



count_g

tmp=count_g["ZEB2",]
mean(tmp[cellleiden$V1==0])/mean(count_g[,cellleiden$V1==0])
mean(tmp[cellleiden$V1==1])/mean(count_g[,cellleiden$V1==1])
mean(tmp[cellleiden$V1==2])/mean(count_g[,cellleiden$V1==2])
mean(tmp[cellleiden$V1==3])/mean(count_g[,cellleiden$V1==3])

tmp=count_g[rownames(count_g) %in% emtgene$V1,]
mean(tmp[,cellleiden$V1==0])/mean(count_g[,cellleiden$V1==0])
mean(tmp[,cellleiden$V1==1])/mean(count_g[,cellleiden$V1==1])
mean(tmp[,cellleiden$V1==2])/mean(count_g[,cellleiden$V1==2])
mean(tmp[,cellleiden$V1==3])/mean(count_g[,cellleiden$V1==3])

mean(count_g[,cellleiden$V1==0])
mean(count_g[,cellleiden$V1==1])
mean(count_g[,cellleiden$V1==2])
mean(count_g[,cellleiden$V1==3])



gexmg <- CreateSeuratObject(counts = count_g, project = "pbmc10k")
gexmg <- NormalizeData(gexmg,normalization.method="RC")
gexmg <- NormalizeData(gexmg,normalization.method="LogNormalize")
gexmg <- ScaleData(gexmg, verbose = FALSE)

gexmg@meta.data["subclus"]<-cellleiden
gexmg@meta.data["cellclus"]<-cellclus

Idents(object = gexmg) <- "cellclus"
deg2 <- FindMarkers(gexmg, ident.1 = 0, ident.2 = 1)
deg22 <- FindMarkers(gexmg, ident.1 = 0, ident.2 = 2)
deg23 <- FindMarkers(gexmg, ident.1 = 0, ident.2 = 3)

deg3<-deg2[(deg2$p_val_adj<0.00001)&(deg2$avg_log2FC>0),]
deg32<-deg22[(deg22$p_val_adj<0.00001)&(deg22$avg_log2FC>0),]
deg33<-deg23[(deg23$p_val_adj<0.00001)&(deg23$avg_log2FC>0),]

tmpdef=intersect(rownames(deg3),rownames(deg32))
deg4=intersect(tmpdef,rownames(deg33))

write.table(rownames(deg3),"deg3gene.txt",sep="\t",quote=F,row.names=F,col.names=F)
write.table(deg4,"deg4gene.txt",sep="\t",quote=F,row.names=F,col.names=F)

log1pg=(gexmg@assays$RNA@data)
log1pgs=(gexmg@assays$RNA@scale.data)

gexmg_sub=subset(gexmg, subset = subclus>=0)
Idents(object = gexmg_sub) <- "subclus"
deg2 <- FindMarkers(gexm_sub, ident.1 = 1, ident.2 = 0)




mean(log1pg[,cellleiden$V1==0])
mean(log1pg[,cellleiden$V1==1])
mean(log1pg[,cellleiden$V1==2])
mean(log1pg[,cellleiden$V1==3])

tmp=log1pg[rownames(count_g) %in% emtgene$V1,]
mean(tmp[,cellleiden$V1==0])
mean(tmp[,cellleiden$V1==1])
mean(tmp[,cellleiden$V1==2])
mean(tmp[,cellleiden$V1==3])


tmp=log1pg[rownames(count_g) %in% tgfbgene$V1,]
mean(tmp[,cellleiden$V1==0])
mean(tmp[,cellleiden$V1==1])
mean(tmp[,cellleiden$V1==2])
mean(tmp[,cellleiden$V1==3])

tmp=log1pg[rownames(count_g) %in% krasgene$V1,]
mean(tmp[,cellleiden$V1==0])
mean(tmp[,cellleiden$V1==1])
mean(tmp[,cellleiden$V1==2])
mean(tmp[,cellleiden$V1==3])

tmp=log1pg[rownames(count_g) %in% krasdngene$V1,]
mean(tmp[,cellleiden$V1==0])
mean(tmp[,cellleiden$V1==1])


tmp=log1pgs[rownames(count_g) %in% nfkbgene$V1,]
mean(tmp[,cellleiden$V1==0])
mean(tmp[,cellleiden$V1==1])


tmp=log1pgs[rownames(count_g) %in% krasdngene$V1,]
mean(tmp[,cellleiden$V1==0])
mean(tmp[,cellleiden$V1==1])


tmp=log1pg[rownames(count_g) %in% nfkbgene$V1,]
mean(tmp[,cellleiden$V1==0])
mean(tmp[,cellleiden$V1==1])

mean(tmp[,cellleiden$V1==0])/mean(log1pg[,cellleiden$V1==0])
mean(tmp[,cellleiden$V1==1])/mean(log1pg[,cellleiden$V1==1])
mean(tmp[,cellleiden$V1==2])/mean(log1pg[,cellleiden$V1==2])
mean(tmp[,cellleiden$V1==3])/mean(log1pg[,cellleiden$V1==3])


tmp=log1pg["CDH1",]
mean(tmp[cellleiden$V1==0])/mean(log1pg[,cellleiden$V1==0])
mean(tmp[cellleiden$V1==1])/mean(log1pg[,cellleiden$V1==1])
mean(tmp[cellleiden$V1==2])/mean(log1pg[,cellleiden$V1==2])
mean(tmp[cellleiden$V1==3])/mean(log1pg[,cellleiden$V1==3])


tmp=log1pg["PDGFRA",]
mean(tmp[cellclus$V1==0])
mean(tmp[cellclus$V1==1])
mean(tmp[cellclus$V1==2])
mean(tmp[cellclus$V1==3])
mean(tmp[cellclus$V1==4])
mean(tmp[cellclus$V1==5])
mean(tmp[cellclus$V1==6])


tmp=log1pg["GFAP",]
mean(tmp[cellclus$V1==0])
mean(tmp[cellclus$V1==1])
mean(tmp[cellclus$V1==2])
mean(tmp[cellclus$V1==3])
mean(tmp[cellclus$V1==4])
mean(tmp[cellclus$V1==5])
mean(tmp[cellclus$V1==6])

tmp=log1pg["NRF1",]
mean(tmp[cellclus$V1==0])
mean(tmp[cellclus$V1==1])
mean(tmp[cellclus$V1==2])
mean(tmp[cellclus$V1==3])
mean(tmp[cellclus$V1==4])
mean(tmp[cellclus$V1==5])
mean(tmp[cellclus$V1==6])

tmp=log1pg["S100A8",]
mean(tmp[cellclus$V1==0])
mean(tmp[cellclus$V1==1])
mean(tmp[cellclus$V1==2])
mean(tmp[cellclus$V1==3])
mean(tmp[cellclus$V1==4])
mean(tmp[cellclus$V1==5])
mean(tmp[cellclus$V1==6])

tmp=log1pg["FGF1",]
mean(tmp[cellclus$V1==0])
mean(tmp[cellclus$V1==1])
mean(tmp[cellclus$V1==2])
mean(tmp[cellclus$V1==3])
mean(tmp[cellclus$V1==4])
mean(tmp[cellclus$V1==5])
mean(tmp[cellclus$V1==6])


tmp=log1pg["ATF3",]
mean(tmp[cellclus$V1==0])
mean(tmp[cellclus$V1==1])
mean(tmp[cellclus$V1==2])
mean(tmp[cellclus$V1==3])
mean(tmp[cellclus$V1==4])
mean(tmp[cellclus$V1==5])
mean(tmp[cellclus$V1==6])

tmp=log1pgs["SOX2",]
mean(tmp[cellclus$V1==0])
mean(tmp[cellclus$V1==1])
mean(tmp[cellclus$V1==2])
mean(tmp[cellclus$V1==3])
mean(tmp[cellclus$V1==4])
mean(tmp[cellclus$V1==5])
mean(tmp[cellclus$V1==6])

tmp=log1pgs["CDK1",]
mean(tmp[cellclus$V1==0])
mean(tmp[cellclus$V1==1])
mean(tmp[cellclus$V1==2])
mean(tmp[cellclus$V1==3])
mean(tmp[cellclus$V1==4])
mean(tmp[cellclus$V1==5])
mean(tmp[cellclus$V1==6])

tmp=log1pg["CREM",]
mean(tmp[cellleiden$V1==0])
mean(tmp[cellleiden$V1==1])

library(fgsea)

clusFC<-log2(apply(log1pg[,cellleiden$V1==0],1,mean)/apply(log1pg[,cellleiden$V1==1],1,mean))
original_gene_list <- clusFC
names(original_gene_list) <- rownames(log1pg)
gene_list<-original_gene_list[is.finite(original_gene_list)]
gene_list = sort(gene_list, decreasing = TRUE)

msigdbr_df = msigdbr(species = "human", category = "H")
msigdbr_list = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)
fgseaRes <- fgsea(msigdbr_list, gene_list, minSize=15, maxSize = 500, nperm=1000)
head(fgseaRes[order(padj, -abs(NES)), ], n=20)
fgseaRes[ES > 0][head(order(padj, -abs(NES)), n=30), pathway]
fgseaRes[ES < 0][head(order(padj, -abs(NES)), n=30), pathway]


library(escape)
gene.sets3 <- getGeneSets(library = "H")

ES.seurat1 <- enrichIt(obj = gexm_sub, 
                       gene.sets = gene.sets3, 
                       groups = 1000, cores = 2)

write.table(ES.seurat1,"ssGSEA.csv",sep=",",quote=F,row.names=T,col.names=T)

ES.seurat1<-read.csv("ssGSEA.csv")

wdiff=matrix(0, nrow = 50, ncol = 2)

ESo=apply(ES.seurat1[leidensub==0,],2,mean)-apply(ES.seurat1[leidensub==1,],2,mean)
ESo[ESo<(-0.02)]

leidensub=gexm_sub@meta.data$subclus

for (i in 1:50){
  tmptest=wilcox.test(ES.seurat1[leidensub==0,i], ES.seurat1[leidensub==1,i], paired = FALSE)
  wdiff[i,1]=tmptest$p.value
  wdiff[i,2]=ESo[i]
}



rownames(wdiff)=colnames(ES.seurat1)


write.table(wdiff,"wdiff.csv",sep=",",quote=F,row.names=T,col.names=T)



tmp2=tmp[cellclus==0]
tmp3=tfact["ZEB1"][cellclus==0]

cor(tmp2[tmp2!=0],tmp3[tmp2!=0])
plot(tmp2[tmp2!=0],tmp3[tmp2!=0])

tmpemt=apply(log1pgs[rownames(count_g) %in% emtgene$V1,],2,mean)
cor(tmpemt[cellclus==0],tmp3)

plot(tmpemt[cellclus==0],tmp3)


mean(log1pgs[,cellleiden$V1==0])
mean(log1pgs[,cellleiden$V1==1])
mean(log1pgs[,cellleiden$V1==2])
mean(log1pgs[,cellleiden$V1==3])

tmp=log1pgs[rownames(log1pgs) %in% emtgene$V1,]
mean(tmp[,cellleiden$V1==0])
mean(tmp[,cellleiden$V1==1])
mean(tmp[,cellleiden$V1==2])
mean(tmp[,cellleiden$V1==3])


tmp=log1pgs[rownames(log1pgs) %in% krasgene$V1,]
mean(tmp[,cellleiden$V1==0])
mean(tmp[,cellleiden$V1==1])
mean(tmp[,cellleiden$V1==2])
mean(tmp[,cellleiden$V1==3])

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

for (i in 1:50){
  tmpgene=as.data.frame(msigdbr_list_h[i])
  tmpgene=tmpgene[,1]
  tmp=log1pg[rownames(log1pgs) %in% tmpgene,]
  tmp20=apply(tmp[,cellleiden$V1==0],1,mean)
  tmp21=apply(tmp[,cellleiden$V1==1],1,mean)
  tmp30=log2(tmp20)
  tmp31=log2(tmp21)
  tmp4=tmp30-tmp31
  hdiff[i,]=mean(tmp4[is.finite(tmp4)])
}

dim(tmp)

msigdbr_list_h0<-msigdbr_list_h[hdiff>0]
msigdbr_list_h1<-msigdbr_list_h[hdiff<(-0.07)]

tmp=log1pgs[rownames(count_g) %in% emtgene$V1,]
emtexp=apply(tmp,2,mean)
write.table(emtexp,"emtexp.txt",sep="\t",quote=F,row.names=F,col.names=F)

tmp=log1pgs[rownames(count_g) %in% nfkbgene$V1,]
nfkbexp=apply(tmp,2,mean)
write.table(nfkbexp,"nfkbexp.txt",sep="\t",quote=F,row.names=F,col.names=F)

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

genelist

FeaturePlot(gexm, features="SOX4")
FeaturePlot(gexm, features="KCNIP4")
FeaturePlot(gexm, features="IGF2BP2")


####################


peaks_s=count$Peaks[,celllist$V1]

peakbed=read.table("../peaks_hg38_tag.bed",sep="\t",header=F)
peaktag=peakbed$V4
peaktag=sub("-", ":", peaktag)
selectpeaks=peaks_s[peaktag,]

selectpeaks2=selectpeaks[,cellclus$V1==0]

dim(selectpeaks)

dim(selectpeaks[tumorpeaklist$V1,])


tumorpeaks=selectpeaks[tumorpeaklist$V1,]
tumorpeaks2=tumorpeaks[,cellclus$V1==0]

peaks_t=peaks_s[,cellleiden$V1>=0]

peaks_assay <- CreateAssayObject(counts = peaks_t)

pbmc <- CreateSeuratObject(
  counts = peaks_assay,
  assay = "peaks",
)

pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
DimPlot(object = pbmc, label = TRUE) + NoLegend()
pbmc@meta.data["originalclus"]=cellleiden[cellleiden>=0]
DimPlot(object = pbmc, group.by="originalclus",label = T)

Idents(object = pbmc) <- "originalclus"

da_peaks <- FindMarkers(
  object = pbmc,
  ident.1 = "1",
  ident.2 = "0",
  test.use = 'LR',
  logfc.threshold = 0.05,
  latent.vars = 'nCount_peaks'
)

uppeaks=da_peaks[(da_peaks$avg_log2FC>0)&(da_peaks$p_val_adj<0.05),]
downpeaks=da_peaks[(da_peaks$avg_log2FC<0)&(da_peaks$p_val_adj<0.05),]

rownames(uppeaks)


uptmp=str_split(rownames(uppeaks), pattern = ":", simplify = TRUE)
uptmp2=str_split(uptmp[,2], pattern = "-", simplify = TRUE)
upout=cbind(uptmp[,1],uptmp2)

downtmp=str_split(rownames(downpeaks), pattern = ":", simplify = TRUE)
downtmp2=str_split(downtmp[,2], pattern = "-", simplify = TRUE)
downout=cbind(downtmp[,1],downtmp2)

uptag=paste(uptmp[,1],uptmp[,2],sep="-")
upout=peakbed[peakbed$V4 %in% uptag,]

downtag=paste(downtmp[,1],downtmp[,2],sep="-")
downout=peakbed[peakbed$V4 %in% downtag,]

write.table(upout,"upout.bed",sep="\t",quote=F,row.names=F,col.names=F)
write.table(downout,"downout.bed",sep="\t",quote=F,row.names=F,col.names=F)

mean(peaks_t["chr7:106808702-106809572",])


pbmc



peaks_assay <- CreateAssayObject(counts = selectpeaks)

pbmc <- CreateSeuratObject(
  counts = peaks_assay,
  assay = "peaks",
)

pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
DimPlot(object = pbmc, label = TRUE) + NoLegend()
pbmc@meta.data["originalclus"]=cellclus
DimPlot(object = pbmc, group.by="originalclus",label = T)


rownames(pbmc@meta.data)[1]
celllist$V1[1]

peaks_assay2 <- CreateAssayObject(counts = tumorpeaks)

tumor <- CreateSeuratObject(
  counts = peaks_assay2,
  assay = "peaks",
)

tumor <- RunTFIDF(tumor)
tumor <- FindTopFeatures(tumor, min.cutoff = 'q0')
tumor <- RunSVD(tumor)
tumor <- RunUMAP(object = tumor, reduction = 'lsi', dims = 2:30)
tumor <- FindNeighbors(object = tumor, reduction = 'lsi', dims = 2:30)
tumor <- FindClusters(object = tumor, verbose = FALSE, algorithm = 3)
DimPlot(object = tumor, label = TRUE) + NoLegend()
tumor@meta.data["tumorclus"]=cellclus
DimPlot(object = tumor, group.by="tumorclus",label = T)







peaks_assay <- CreateAssayObject(counts = selectpeaks2)

pbmc <- CreateSeuratObject(
  counts = peaks_assay,
  assay = "peaks",
)

pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
DimPlot(object = pbmc, label = TRUE) + NoLegend()
pbmc@meta.data["originalclus"]=tumorclus[cellclus==0]
DimPlot(object = pbmc, group.by="originalclus",label = T)


rownames(tumorpeaks2)[8000]



peaks_assay2 <- CreateAssayObject(counts = tumorpeaks2)

tumor <- CreateSeuratObject(
  counts = peaks_assay2,
  assay = "peaks",
)

tumor <- RunTFIDF(tumor)
tumor <- FindTopFeatures(tumor, min.cutoff = 'q0')
tumor <- RunSVD(tumor)
tumor <- RunUMAP(object = tumor, reduction = 'lsi', dims = 2:30)
tumor <- FindNeighbors(object = tumor, reduction = 'lsi', dims = 2:30)
tumor <- FindClusters(object = tumor, verbose = FALSE, algorithm = 3)
DimPlot(object = tumor, label = TRUE) + NoLegend()
tumor@meta.data["tumorclus"]=tumorclus[cellclus==0]
DimPlot(object = tumor, group.by="tumorclus",label = T)




#############

peaks_assay <- CreateAssayObject(counts = selectpeaks)

#gexm_cut <- CreateSeuratObject(counts = count_s, project = "neuro",min.cells = 500)
peaks_assay <- CreateAssayObject(counts = peaks_s)
gexm_cut[["ATAC"]] <- peaks_assay
#gexm_cut <- subset(gexm_cut, subset = nFeature_RNA > 500)

gexm_cut[["percent.mt"]] <- PercentageFeatureSet(gexm_cut, pattern = "^MT-")
VlnPlot(gexm_cut, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(gexm_cut, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(gexm_cut, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

DefaultAssay(gexm_cut) <- "RNA"
gexm_cut <- NormalizeData(gexm_cut)
gexm_cut <- FindVariableFeatures(gexm_cut, selection.method = "vst", nfeatures = 2000)
gexm_cut <- ScaleData(gexm_cut, verbose = FALSE)
gexm_cut <- RunPCA(gexm_cut, npcs = 30, verbose = FALSE)
gexm_cut <- RunUMAP(gexm_cut, reduction = "pca", dims = 1:30)
gexm_cut <- FindNeighbors(gexm_cut, reduction = "pca", dims = 1:30)
gexm_cut <- FindClusters(gexm_cut, resolution = 0.1)
DimPlot(gexm_cut, reduction = "umap",group.by="rna_cluster",label = T)

DefaultAssay(gexm_cut) <- "ATAC"
gexm_cut <- RunTFIDF(gexm_cut)
gexm_cut <- FindTopFeatures(gexm_cut, min.cutoff = 'q50')
gexm_cut <- RunSVD(gexm_cut)
gexm_cut <- RunUMAP(object = gexm_cut, reduction = 'lsi', dims = 2:30)
gexm_cut <- FindNeighbors(object = gexm_cut, reduction = 'lsi', dims = 2:30)
gexm_cut <- FindClusters(object = gexm_cut, verbose = FALSE, algorithm = 3, resolution = 0.1)
DimPlot(object = gexm_cut, label = TRUE,group.by="RNA_snn_res.0.1")

motif_cell<-read.csv("motif_cell_activity.csv",header=T,row.names=1)
motif_assay <- CreateAssayObject(counts = t(motif_cell))
gexm_cut[["MOTIF"]] <- motif_assay

celllist_100<-read.csv("celllist_pca100.txt",header=F)
gexm_cut$pc100<-(rownames(gexm_cut@meta.data) %in% celllist_100$V1)
DimPlot(gexm_cut, reduction = "umap",group.by="pc100",label = T)

celllist_tumor<-read.csv("cellcluster_tumor_celltag.txt",header=F)
celllist_tumor_clus<-read.csv("cellcluster_tumor.txt",header=F)

tumorclus=as.data.frame(cbind(celllist_tumor,celllist_tumor_clus))
gexm_cut$tumorclus<-6
gexm_cut$tumorclus[tumorclus[,1]]<-tumorclus[,2]
DimPlot(gexm_cut, reduction = "umap",group.by="tumorclus",label = T)

gexm_cut$tumorclus_0<-(gexm_cut$tumorclus==0)
DimPlot(gexm_cut, reduction = "umap",group.by="tumorclus_0",label = T)

gexm_cut$tumorclus_1<-(gexm_cut$tumorclus==1)
DimPlot(gexm_cut, reduction = "umap",group.by="tumorclus_1",label = T)

gexm_cut$tumorclus_2<-(gexm_cut$tumorclus==2)
DimPlot(gexm_cut, reduction = "umap",group.by="tumorclus_2",label = T)

gexm_cut$tumorclus_3<-(gexm_cut$tumorclus==3)
DimPlot(gexm_cut, reduction = "umap",group.by="tumorclus_3",label = T)

gexm_cut$tumorclus_4<-(gexm_cut$tumorclus==4)
DimPlot(gexm_cut, reduction = "umap",group.by="tumorclus_4",label = T)

DefaultAssay(gexm_cut) <- "RNA"
DimPlot(gexm_cut, reduction = "umap",group.by="RNA_snn_res.0.1",label = T)

gexm_cut_t <- subset(gexm_cut, subset = tumorclus!=6)

DefaultAssay(gexm_cut_t) <- "RNA"
gexm_cut_t <- NormalizeData(gexm_cut_t)
gexm_cut_t <- FindVariableFeatures(gexm_cut_t, selection.method = "vst", nfeatures = 2000)
gexm_cut_t <- ScaleData(gexm_cut_t, verbose = FALSE)
gexm_cut_t <- RunPCA(gexm_cut_t, npcs = 30, verbose = FALSE)
gexm_cut_t <- RunUMAP(gexm_cut_t, reduction = "pca", dims = 1:30)
DimPlot(object = gexm_cut_t, label = TRUE,group.by="tumorclus")

gexm_cut@meta.data
gexm_cut_t@meta.data

gexm_cut_t$seurat_clusters<-gexm_cut_t$tumorclus

Idents(object = gexm_cut_t) <- "tumorclus"

fgene_t <- FindAllMarkers(gexm_cut_t, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.05)
fgene_st<-fgene_t %>%
  group_by(cluster) %>%
  slice_max(n = 500, order_by = avg_log2FC)


clus0t=fgene_t[fgene_t$cluster==0,]
clus0st=clus0t[(clus0t$p_val_adj<0.1)&(clus0t$avg_log2FC>0.05),]
enhalist<-read.csv("enhancer_clus0_tumor.csv",header=F)
targetgene=str_split(enhalist$V1, pattern = "_", simplify = TRUE)[,1]
gene_list<-unique(targetgene)
tumorgene_0st=gene_list[(gene_list %in% clus0st$gene)]


clus1t=fgene_t[fgene_t$cluster==1,]
clus1st=clus1t[(clus1t$p_val_adj<0.1)&(clus1t$avg_log2FC>0.05),]
enhalist<-read.csv("enhancer_clus1_tumor.csv",header=F)
targetgene=str_split(enhalist$V1, pattern = "_", simplify = TRUE)[,1]
gene_list<-unique(targetgene)
tumorgene_1st=gene_list[(gene_list %in% clus1st$gene)]


clus3t=fgene_t[fgene_t$cluster==3,]
clus3st=clus3t[(clus3t$p_val_adj<0.1)&(clus3t$avg_log2FC>0.05),]
enhalist<-read.csv("enhancer_clus3_tumor.csv",header=F)
targetgene=str_split(enhalist$V1, pattern = "_", simplify = TRUE)[,1]
gene_list<-unique(targetgene)
tumorgene_3st=gene_list[(gene_list %in% clus3st$gene)]




DefaultAssay(gexm_cut_t) <- "ATAC"
gexm_cut_t <- RunTFIDF(gexm_cut_t)
gexm_cut_t <- FindTopFeatures(gexm_cut_t, min.cutoff = 'q50')
gexm_cut_t <- RunSVD(gexm_cut_t)
gexm_cut_t <- RunUMAP(object = gexm_cut_t, reduction = 'lsi', dims = 2:30)
gexm_cut_t <- FindNeighbors(object = gexm_cut_t, reduction = 'lsi', dims = 2:30)
gexm_cut_t <- FindClusters(object = gexm_cut_t, verbose = FALSE, algorithm = 3, resolution = 0.1)

DimPlot(object = gexm_cut_t, label = TRUE,group.by="tumorclus")

DimPlot(gexm_cut_t, reduction = "umap",group.by="tumorclus_0",label = T)
DimPlot(gexm_cut_t, reduction = "umap",group.by="tumorclus_1",label = T)
DimPlot(gexm_cut_t, reduction = "umap",group.by="tumorclus_2",label = T)
DimPlot(gexm_cut_t, reduction = "umap",group.by="tumorclus_3",label = T)
DimPlot(gexm_cut_t, reduction = "umap",group.by="tumorclus_4",label = T)




gexm_cut@meta.data

mean(gexm_cut$nCount_ATAC[gexm_cut$pc100==TRUE])
mean(gexm_cut$nCount_ATAC[gexm_cut$pc100!=TRUE])

fgene <- FindAllMarkers(gexm_cut, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
fgene_s<-fgene %>%
  group_by(cluster) %>%
  slice_max(n = 500, order_by = avg_log2FC)

write.table(unique(fgene_s$gene),"rna_vargene.txt",quote=F,sep="\t",col.names=F,row.names=F,eol = "\n")


gene01<-FindMarkers(gexm_cut, ident.1 = "0", ident.2 = "1", group.by = "rna_cluster", min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE)
gene02<-FindMarkers(gexm_cut, ident.1 = "0", ident.2 = "2", group.by = "rna_cluster", min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE)


DefaultAssay(gexm_cut) <- "RNA"
gene012<-FindMarkers(gexm_cut, ident.1 = "0", ident.2 = c("1","2"), group.by = "rna_cluster", min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE)
#enhancer_act<-read.csv("enhancer_nonorm_meanact_clus0_to12.csv",header=T)
enhancer_act<-read.csv("enhancer_nonorm_meanact_clus0_to12_adj.csv",header=T)
colnames(enhancer_act)=c("enhancer","activity_change")
targetgene=str_split(enhancer_act$enhancer, pattern = "_", simplify = TRUE)[,1]
enhancer_act$target<-targetgene

enhancer_act
gene012df<-as.data.frame(cbind(rownames(gene012),gene012))
colnames(gene012df)[1]<-"target"

enhancer_df<-inner_join(enhancer_act,gene012df,by="target")

enhancer_calc<-as.data.frame(table(enhancer_df$target))
colnames(enhancer_calc)<-c("target","Freq")
enhancer_calc$meanact<-0
enhancer_calc$meanact_nonabs<-0

"SOX10" %in% rownames(gene012)
"SOX10" %in% enhancer_act$target

for(i in 1:length(enhancer_calc$target)){
  enhancer_calc$meanact[i]=mean(abs(enhancer_df$activity_change[enhancer_df$target==enhancer_calc$target[i]]))
  enhancer_calc$meanact_nonabs[i]=mean((enhancer_df$activity_change[enhancer_df$target==enhancer_calc$target[i]]))
}


enhancer_df[enhancer_df$activity_change>0.5,]

hist(enhancer_df$activity_change,breaks=100,range=c(-1,1))

hist(enhancer_calc$meanact_nonabs)

plot(enhancer_df$avg_log2FC,enhancer_df$activity_change)
plot(enhancer_calc$Freq,enhancer_calc$meanact)
plot(enhancer_calc$Freq,enhancer_calc$meanact_nonabs)

enhancer_df_tmp<-enhancer_df[enhancer_df$activity_change>0.3,]

library(plotly)
fig <- plot_ly(x = -1*log(enhancer_df$p_val_adj), y = enhancer_df$activity_change)
fig2 <- subplot(
  fig %>% add_markers(alpha = 0.2),
  fig %>% add_histogram2d(nbinsx=100,nbinsy=100)
)
fig2

fig <- plot_ly(x = -1*log(enhancer_df$p_val_adj), y = enhancer_df$avg_log2FC)
fig2 <- subplot(
  fig %>% add_markers(alpha = 0.2),
  fig %>% add_histogram2d(nbinsx=100,nbinsy=100)
)
fig2


fig <- plot_ly(x = enhancer_calc$Freq, y = enhancer_calc$meanact_nonabs)
fig2 <- subplot(
  fig %>% add_markers(alpha = 0.2),
  fig %>% add_histogram2d(nbinsx=100,nbinsy=100)
)
fig2

fig <- plot_ly(x = enhancer_calc$Freq, y = log(enhancer_calc$meanact))
fig2 <- subplot(
  fig %>% add_markers(alpha = 0.2),
  fig %>% add_histogram2d(nbinsx=100,nbinsy=100)
)
fig2

enhancer_calc[log(enhancer_calc$meanact)>0,]

gene0102<-rownames(gene01)[rownames(gene01) %in% rownames(gene02)]



"NFKB1" %in% rownames(gene01)

fgene$cluster

clus0=fgene[fgene$cluster==0,]

clus0s=clus0[(clus0$p_val_adj<0.001)&(clus0$avg_log2FC>0.5),]
clus0sc=clus0[(clus0$p_val_adj<0.001)&(clus0$avg_log2FC>1),]

enhalist<-read.csv("enhancer_clus0_all.csv",header=F)
targetgene=str_split(enhalist$V1, pattern = "_", simplify = TRUE)[,1]
gene_list<-unique(targetgene)
tumorgene=gene_list[(gene_list %in% clus0s$gene)]
tumorgene_n=gene_list[(gene_list %in% clus0sc$gene)]

hvar<-read.csv("high_variance_gene.txt",header=F)
tumorgene_h=gene_list[(gene_list %in% hvar$V1)]

enhalist_hv<-read.csv("DAE_mtx_deeplift_high_variance.csv",header=T)
enhalist_hv=enhalist_hv[,-1]
targetgene=str_split(enhalist_hv[1:1000,1], pattern = "_", simplify = TRUE)[,1]
gene_list<-unique(targetgene)
tumorgene=gene_list[(gene_list %in% clus0s$gene)]

enhalist<-read.csv("enhancer_hv_clus0_to12.csv",header=F)
targetgene=str_split(enhalist$V1, pattern = "_", simplify = TRUE)[,1]
gene_list<-unique(targetgene)
tumorgene=gene_list[(gene_list %in% gene0102)]
write.table(tumorgene,"cluster0_12_enha_genelist_withRNA.txt",sep="\t",quote=F,col.names=F,row.names=F)

"NFKB1" %in% gene0102
"NFKB1" %in% gene_list
length(gene_list)
length(gene0102)

length(gene_list)
length(hvar$V1)
length(tumorgene_h)

write.table(tumorgene,"cluster0_all_enha_genelist_withRNA.txt",sep="\t",quote=F,col.names=F,row.names=F)

cluster0.markers <- FindMarkers(gexm_cut, ident.1 = 0, ident.2 = 1, min.pct = 0.25)

cluster0.markersp<-cluster0.markers[cluster0.markers$avg_log2FC>0,]

cluster0.markerspm<-cbind(rownames(cluster0.markersp),cluster0.markersp)


DefaultAssay(gexm_cut) <- "RNA"
FeaturePlot(object = gexm_cut, features = 'FABP7')
FeaturePlot(object = gexm_cut, features = 'AQP4')
FeaturePlot(object = gexm_cut, features = 'PDGFRA')
FeaturePlot(object = gexm_cut, features = 'S100B')
FeaturePlot(object = gexm_cut, features = 'TUBB')
FeaturePlot(object = gexm_cut, features = 'CCR5')
FeaturePlot(object = gexm_cut, features = 'FOXJ1')
FeaturePlot(object = gexm_cut, features = 'OLIG1')
FeaturePlot(object = gexm_cut, features = 'SOX15')
FeaturePlot(object = gexm_cut, features = 'NFKB1')
FeaturePlot(object = gexm_cut, features = 'SOX8')

FeaturePlot(object = gexm_cut, features = 'CD55')
FeaturePlot(object = gexm_cut, features = 'SIRT2')
FeaturePlot(object = gexm_cut, features = 'MRRF')

DefaultAssay(gexm_cut) <- "MOTIF"
FeaturePlot(object = gexm_cut, features = 'SOX18')

count_cell=t(count_s)
count_cell=t(adjmtx)
count_cell=adjmtx

sum(colnames(count_cell)=="SOX10")
sum(colnames(motif_cell)=="SOX10")

colnames(motif_cell)=toupper(colnames(motif_cell))

unique(colnames(motif_cell))

count_cell_s=count_cell[,colnames(count_cell) %in% colnames(motif_cell)]
count_cell_smtx=as.data.frame(as.matrix(count_cell_s))
motif_cell_s=motif_cell[,colnames(motif_cell) %in% colnames(count_cell_smtx)]

motif_cell_su=motif_cell_s[,unique(colnames(motif_cell_s))]
count_cell_smtx_u=count_cell_smtx[,unique(colnames(count_cell_smtx))]
count_cell_smtx_u=count_cell_smtx_u[,colnames(motif_cell_su)]

motifnum=dim(motif_cell_su)[2]
cormotif=matrix(0, motifnum)
rownames(cormotif)=colnames(motif_cell_su)

for(i in 1:motifnum){
  cormotif[i]=cor(count_cell_smtx_u[,i],motif_cell_su[,i])
}

sub_motif_act=apply(motif_cell_su[cellcluster==0,],2,mean)-apply(motif_cell_su[cellcluster!=0,],2,mean)



plot(cormotif,sub_motif_act)

sub_motif_act["SOX10"]
cormotif["SOX10"]

length(cormotif)

plotdf=as.data.frame(cbind(cormotif,sub_motif_act))
colnames(plotdf)=c("Correlation","Motifact")

tmpn=rownames(cormotif)

tmpn[cormotif<0.1]<-""
tmpn[sub_motif_act<0.006]<-""

g <- ggplot(plotdf,aes(x=Correlation, y=Motifact))
g <- g + geom_point(size = 8, colour = "black")
g <- g +theme_bw()
g <- g + ggtitle("Enriched TF motifs in tumor-specific enhancers")
g <- g + xlab("Correlation between TF activity vs expression")
g <- g + ylab("Tumor Specificity of TF acitivity")
g <- g +theme(text = element_text("Arial"),
              legend.position = "none",
              legend.title=element_text(size = 24),
              legend.text=element_text(size = 24),
              axis.title.x=element_text(size = 24),
              axis.title.y=element_text(size = 24),
              axis.text.x=element_text(size = 24),
              axis.text.y=element_text(size = 24),
              plot.title = element_text(size = 40,hjust = 0.5),
              panel.grid = element_blank()
)
g <- g+ ggrepel::geom_text_repel(aes(label=tmpn),size=8,point.size = 15,fontface = "bold",segment.color = "black",segment.size=1,min.segment.length = 0.1)
g

png(paste0("TF_tumor.png"), width = 1000, height = 700)  # 描画デバイスを開く
plot(g)
dev.off()

sub_motif_act[(sub_motif_act>0.0025) & (cormotif>0.1)]
sub_motif_act[(sub_motif_act>0.0005) & (cormotif>0.1)]
sub_motif_act["SOX8"]
cormotif[rownames(cormotif)=="SOX10"]

colnames(count_cell)

count_cell_smtx=t(count_cell_smtx)
motif_cell=t(motif_cell)


for(i in tumorgene){
  png(paste0("tumorgene/",i,"_vlnplot.png"), width = 700, height = 500)  
  print(VlnPlot(object = gexm_cut, features = i))
  dev.off()
  
  png(paste0("tumorgene/",i,"_featureplot.png"), width = 700, height = 500)  
  print(FeaturePlot(object = gexm_cut, features = i))
  dev.off()
}
tumorgene

VlnPlot(object = gexm_cut, features = 'CCDC85B')

dim(gexm_cut@meta.data)
write.table(gexm_cut@meta.data$seurat_clusters,"cellcluster.txt",quote=F,sep="\t",col.names=F,row.names=F)
