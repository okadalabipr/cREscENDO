setwd("D:/scbasset_test/multiome_human_neuro/P-1694_S-1694/for_paper")

library(msigdbr)
library(fgsea)




genelist<-read.csv("genelistfull.txt",header=F)
genelist<-genelist$V1



####################


genetag<-read.table("strogenelist.txt",header=F)
gene_list<-genetag$V1

msigdbr_df_c2 = msigdbr(species = "human", category = "C2")
msigdbr_list_c2 = split(x = msigdbr_df_c2$gene_symbol, f = msigdbr_df_c2$gs_name)

msigdbr_df_c3 = msigdbr(species = "human", category = "C3")
msigdbr_list_c3 = split(x = msigdbr_df_c3$gene_symbol, f = msigdbr_df_c3$gs_name)

msigdbr_df_c5 = msigdbr(species = "human", category = "C5")
msigdbr_list_c5 = split(x = msigdbr_df_c5$gene_symbol, f = msigdbr_df_c5$gs_name)

msigdbr_df_h = msigdbr(species = "human", category = "H")
msigdbr_list_h = split(x = msigdbr_df_h$gene_symbol, f = msigdbr_df_h$gs_name)

fgseaRes <- fora(msigdbr_list_c2, gene_list, genelist)
out2<-head(fgseaRes[order(padj), ], n=30)

fgseaRes <- fora(msigdbr_list_c3, gene_list, genelist)
out3_01<-head(fgseaRes[order(padj), ], n=30)

fgseaRes <- fora(msigdbr_list_c5, gene_list, genelist)
out5<-head(fgseaRes[order(padj), ], n=30)

fgseaRes <- fora(msigdbr_list_h, gene_list, genelist)
outh<-head(fgseaRes[order(padj), ], n=30)

write.csv(out2[,1:5],"strogenelist_c2.csv")
write.csv(out5[,1:5],"strogenelist_c5.csv")
write.csv(outh[,1:5],"strogenelist_h.csv")

