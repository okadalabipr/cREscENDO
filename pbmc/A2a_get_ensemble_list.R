library(dplyr)
library(biomaRt)

working_dir_path=""

db <- useMart("ensembl",host="https://jan2020.archive.ensembl.org")
hg <- useDataset("hsapiens_gene_ensembl", mart = db)
res <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),mart = hg)
chr<-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")

genelist<-read.csv(paste0(working_dir_path,"/pair_300000.csv"),header=F)
genelist<-genelist[,1]

res_s_s=res[res$external_gene_name %in% genelist,]
res_s_s=distinct(res_s_s,external_gene_name,.keep_all = T)

write.csv(res_s_s[,c("ensembl_gene_id","external_gene_name")],paste0(working_dir_path,"/ensemblelist.csv"))
