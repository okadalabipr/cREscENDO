library(DIRECTNET)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(patchwork)
library(dplyr)
library(ggplot2)
options(stringsAsFactors = FALSE)

genelist<-read.csv("pair_promoter_for_directnet.csv",header=F)


library(hdf5r)
inputdata.10x <- Read10X_h5("pbmc_10k/filtered_feature_bc_matrix.h5")
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks
genome.info <- read.table(file = "pbmc_10k/hg38.promoter.regions.txt")
names(genome.info) <- c("Chrom","Starts","Ends","genes")
genes <- lapply(genome.info$genes, function(x) strsplit(x,"[|]")[[1]][1])
genes <- lapply(genes, function(x) strsplit(x,"[.]")[[1]][1])
genes <- unlist(genes)
genome.info$genes <- genes
unik <- !duplicated(genes)# filter out different transcript
genome.info <- genome.info[unik,]

pbmc <- CreateSeuratObject(counts = rna_counts)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]

frag.file <- "pbmc_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(counts = atac_counts,sep = c(":", "-"),fragments = frag.file,min.cells = 10)
pbmc[["ATAC"]] <- chrom_assay
pbmc <- subset(x = pbmc,subset = nCount_ATAC < 7e4 &nCount_ATAC > 5e3 &nCount_RNA < 25000 &nCount_RNA > 1000 &percent.mt < 20)

DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

DefaultAssay(pbmc) <- "ATAC"
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(pbmc, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
pbmc <- FindMultiModalNeighbors(pbmc, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pbmc <- FindClusters(pbmc, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
pbmc <- FindSubCluster(pbmc, cluster = 6, graph.name = "wsnn", algorithm = 3)

pbmc <- Run_DIRECT_NET(pbmc, peakcalling = FALSE, k_neigh = 50, atacbinary = TRUE, max_overlap=0.5, size_factor_normalize = FALSE, genome.info = genome.info, focus_markers = genelist$V1)
direct.net_result <- Misc(pbmc, slot = 'direct.net')
saveRDS(direct.net_result, file = "direct_net_enhancerlist_0.obj")