setwd("D:/scbasset_test/multiome_human_neuro/P-1694_S-1694/for_paper")

library(gt)
library(ggplot2)
library(png)



clus0<-read.csv("clus0genenew_allenha_h.csv")
clus1<-read.csv("clus1genenew_allenha_h.csv")
clus2<-read.csv("clus2genenew_allenha_h.csv")
clus3<-read.csv("clus3genenew_allenha_h.csv")

clus0$overlap<-clus0$overlap/clus0$size
clus1$overlap<-clus1$overlap/clus1$size
clus2$overlap<-clus2$overlap/clus2$size
clus3$overlap<-clus3$overlap/clus3$size

clus0<-clus0[c("pathway","pval","padj","overlap")]
clus1<-clus1[c("pathway","pval","padj","overlap")]
clus2<-clus2[c("pathway","pval","padj","overlap")]
clus3<-clus3[c("pathway","pval","padj","overlap")]

clus0<-clus0[clus0$padj<0.001,]
clus1<-clus1[clus1$padj<0.001,]
clus2<-clus2[clus2$padj<0.001,]
clus3<-clus3[clus3$padj<0.001,]

colnames(clus0)<-c("pathway","cl0_pval","cl0_padj","cl0_overlap")
colnames(clus1)<-c("pathway","cl1_pval","cl1_padj","cl1_overlap")
colnames(clus2)<-c("pathway","cl2_pval","cl2_padj","cl2_overlap")
colnames(clus3)<-c("pathway","cl3_pval","cl3_padj","cl3_overlap")


clusall<-dplyr::full_join(clus0,clus1)
clusall<-dplyr::full_join(clusall,clus2)
clusall<-dplyr::full_join(clusall,clus3)



g<-ggplot(df, aes(x = x, y = y, size = size)) +
  geom_point() +
  scale_size(name = "Size", range = c(3, 15))

######################



library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(stringr)

data0<-read.csv("strogenelist_h.csv",row.names=1)


data0["logPadj"]=-1*log(data0["padj"])
data0["GeneRatio"]=data0["overlap"]/data0["size"]


data0$pathway<-str_replace(data0$pathway, "HALLMARK_", "")
data00 <- transform(data0[data0$padj<0.01,], pathway= factor(pathway, levels = rev(data0$pathway[data0$padj<0.05])))


g <- ggplot(data00, aes(x = pathway, y = logPadj, fill = GeneRatio))+ theme(text = element_text(size = 18),axis.title = element_text(size = 20),axis.text.x  = element_text(size = 18),plot.title = element_text(size = 24))
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()+ ggtitle("Strong Enhancer")+
  theme(plot.title = element_text(hjust = 0.5,size = 24),text = element_text(face="bold",size = 20))
g <- g + theme(panel.border = element_rect(fill = NA, size=1))
g <- g + scale_y_continuous(labels = scales::number_format(accuracy = 1))
#g <- g + scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6))
ggsave("stroenha.pdf",width = 9, height = 6)


data00

#######################

stro<-read.csv("strogenelist_h.csv",row.names=1)
stro<-stro[stro$padj<0.05,]

gseadiff<-read.csv("wdiff.csv")
colnames(gseadiff)<-c("Pval","Score_diff")
gseadiff<-gseadiff[gseadiff$`Pval`<0.0001,]
gseadiff<-gseadiff[rownames(gseadiff) %in% stro$pathway,]
gseadiff$pathway<-rownames(gseadiff)
gseadiff$logp<-(-1)*log10(gseadiff$Pval)

gseadiff<-gseadiff[order(gseadiff$Score_diff,decreasing = TRUE),]


gseadiff <- transform(gseadiff, pathway= factor(pathway, levels = gseadiff$pathway))

g <- ggplot(gseadiff, aes(x = pathway, y = Score_diff, fill = logp))+ theme(text = element_text(size = 18),axis.title = element_text(size = 20),axis.text.x  = element_text(size = 18),plot.title = element_text(size = 24))
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()+ ggtitle("ssGSEA score")+
  theme(plot.title = element_text(hjust = 0.5,size = 24),text = element_text(face="bold",size = 20))
g <- g + theme(panel.border = element_rect(fill = NA, size=1))
#g <- g + scale_y_continuous(labels = scales::number_format(accuracy = 1))
#g <- g + scale_color_gradient(low="blue", high="red")
g

ggsave("ssGSEA.pdf",width = 12, height = 6)
#################

library(gt)
library(ggplot2)
library(png)

df=data.frame(matrix(ncol = 3,nrow = 2))
colnames(df)=c("Matched TF","Motif","Q value")
df[1,3]=0.000289
df[2,3]=0.000504
df[1,1]="SOX2"
df[2,1]="CTCF"
df[1,2]=0
df[2,2]=1


tableA<-text_transform(gt(df),locations = cells_body(columns = Motif),
                   fn = function(x) {
                     c(local_image(filename = "glioma_modisco_clus0_new_new2_modiscoreport/MA0143.5.png"),
                       local_image(filename = "glioma_modisco_clus0_new_new2_modiscoreport/MA0139.2.png"))
                   }
)%>%cols_align(align = "center", columns = c("Matched TF","Motif","Q value"))%>%tab_header(
  title = md("Glioma cells")
)%>%tab_options(table.font.size = 14,table.font.weight="bold",heading.title.font.weight="bold",column_labels.font.weight = 'bold')


gtsave(tableA, "./glioma_clus0.pdf")


####

library(gt)
library(ggplot2)
library(png)

df=data.frame(matrix(ncol = 3,nrow = 2))
colnames(df)=c("Matched TF","Motif","Q value")
df[1,3]=0.000357
df[2,3]=0.0265
df[1,1]="NFIC"
df[2,1]="SOX6"
df[1,2]=0
df[2,2]=1


tableB<-text_transform(gt(df),locations = cells_body(columns = Motif),
                       fn = function(x) {
                         c(local_image(filename = "glioma_modisco_clus1_new_new2_modiscoreport/MA01527.2.png"),
                           local_image(filename = "glioma_modisco_clus1_new_new2_modiscoreport/MA0515.1.png"))
                       }
)%>%cols_align(align = "center", columns = c("Matched TF","Motif","Q value"))%>%tab_header(
  title = md("OPC")
)%>%tab_options(table.font.size = 14,table.font.weight="bold",heading.title.font.weight="bold",column_labels.font.weight = 'bold')


gtsave(tableB, "glioma_clus1.pdf")

#####################



df=data.frame(matrix(ncol = 5,nrow = 2))
colnames(df)=c("Forward","Reverse","Matched TF","Motif","Q value")
df[1,5]=0.000338
df[2,5]=0.004793
df[1,3]="FOXJ2::ELF1"
df[2,3]="CTCFL"
df[1,1]=0
df[2,1]=1
df[1,2]=0
df[2,2]=1
df[1,4]=0
df[2,4]=1


tableA<-text_transform(gt(df),
                       locations = cells_body(columns = Forward),
                       fn = function(x) {
                         c(local_image(filename = "pbmc_modisco_type2ig_modiscoreport/trimmed_logos/pos_patterns.pattern_1.cwm.fwd.png"),
                           local_image(filename = "pbmc_modisco_type2ig_modiscoreport/trimmed_logos/pos_patterns.pattern_2.cwm.fwd.png"))
                       }
)%>%text_transform(locations = cells_body(columns = Reverse),
                   fn = function(x) {
                     c(local_image(filename = "pbmc_modisco_type2ig_modiscoreport/trimmed_logos/pos_patterns.pattern_1.cwm.rev.png"),
                       local_image(filename = "pbmc_modisco_type2ig_modiscoreport/trimmed_logos/pos_patterns.pattern_2.cwm.rev.png"))
                   }
)%>%text_transform(locations = cells_body(columns = Motif),
                   fn = function(x) {
                     c(local_image(filename = "pbmc_modisco_type2ig_modiscoreport/MA1952.2.png"),
                       local_image(filename = "pbmc_modisco_type2ig_modiscoreport/MA1102.3.png"))
                   }
)%>%cols_align(align = "center", columns = c("Forward","Reverse","Matched TF","Motif","Q value"))%>%tab_header(
  title = md("CD4 T cell (Cluster 2)")
)%>%tab_options(table.font.size = 14,table.font.weight="bold",heading.title.font.weight="bold",column_labels.font.weight = 'bold')



gtsave(tableA, "Fig2_clus2.png")



####################




df=data.frame(matrix(ncol = 5,nrow = 2))
colnames(df)=c("Forward","Reverse","Matched TF","Motif","Q value")
df[1,5]=0.006185
df[2,5]=0.000849
df[1,3]="FOXO1::ELK3"
df[2,3]="GABPA"
df[1,1]=0
df[2,1]=1
df[1,2]=0
df[2,2]=1
df[1,4]=0
df[2,4]=1


tableA<-text_transform(gt(df),
                       locations = cells_body(columns = Forward),
                       fn = function(x) {
                         c(local_image(filename = "pbmc_modisco_type3ig_modiscoreport/trimmed_logos/pos_patterns.pattern_0.cwm.fwd.png"),
                           local_image(filename = "pbmc_modisco_type3ig_modiscoreport/trimmed_logos/pos_patterns.pattern_5.cwm.fwd.png"))
                       }
)%>%text_transform(locations = cells_body(columns = Reverse),
                   fn = function(x) {
                     c(local_image(filename = "pbmc_modisco_type3ig_modiscoreport/trimmed_logos/pos_patterns.pattern_0.cwm.rev.png"),
                       local_image(filename = "pbmc_modisco_type3ig_modiscoreport/trimmed_logos/pos_patterns.pattern_5.cwm.rev.png"))
                   }
)%>%text_transform(locations = cells_body(columns = Motif),
                   fn = function(x) {
                     c(local_image(filename = "pbmc_modisco_type3ig_modiscoreport/MA1955.2.png"),
                       local_image(filename = "pbmc_modisco_type3ig_modiscoreport/MA0062.4.png"))
                   }
)%>%cols_align(align = "center", columns = c("Forward","Reverse","Matched TF","Motif","Q value"))%>%tab_header(
  title = md("CD4 T cell (Cluster 3)")
)%>%tab_options(table.font.size = 14,table.font.weight="bold",heading.title.font.weight="bold",column_labels.font.weight = 'bold')



gtsave(tableA, "Fig2_clus3.png")





####################




df=data.frame(matrix(ncol = 5,nrow = 2))
colnames(df)=c("Forward","Reverse","Matched TF","Motif","Q value")
df[1,5]=0.020701
df[2,5]=0.000152
df[1,3]="ZEB1"
df[2,3]="GABPA"
df[1,1]=0
df[2,1]=1
df[1,2]=0
df[2,2]=1
df[1,4]=0
df[2,4]=1


tableA<-text_transform(gt(df),
                       locations = cells_body(columns = Forward),
                       fn = function(x) {
                         c(local_image(filename = "pbmc_modisco_type4ig_modiscoreport/trimmed_logos/pos_patterns.pattern_0.cwm.fwd.png"),
                           local_image(filename = "pbmc_modisco_type4ig_modiscoreport/trimmed_logos/pos_patterns.pattern_6.cwm.fwd.png"))
                       }
)%>%text_transform(locations = cells_body(columns = Reverse),
                   fn = function(x) {
                     c(local_image(filename = "pbmc_modisco_type4ig_modiscoreport/trimmed_logos/pos_patterns.pattern_0.cwm.rev.png"),
                       local_image(filename = "pbmc_modisco_type4ig_modiscoreport/trimmed_logos/pos_patterns.pattern_6.cwm.rev.png"))
                   }
)%>%text_transform(locations = cells_body(columns = Motif),
                   fn = function(x) {
                     c(local_image(filename = "pbmc_modisco_type4ig_modiscoreport/MA0103.4.png"),
                       local_image(filename = "pbmc_modisco_type4ig_modiscoreport/MA0062.4.png"))
                   }
)%>%cols_align(align = "center", columns = c("Forward","Reverse","Matched TF","Motif","Q value"))%>%tab_header(
  title = md("CD8 T cell (Cluster 4)")
)%>%tab_options(table.font.size = 14,table.font.weight="bold",heading.title.font.weight="bold",column_labels.font.weight = 'bold')



gtsave(tableA, "Fig2_clus4.png")





####################




df=data.frame(matrix(ncol = 5,nrow = 2))
colnames(df)=c("Forward","Reverse","Matched TF","Motif","Q value")
df[1,5]=0.000137
df[2,5]=0.032144
df[1,3]="SPIB"
df[2,3]="BCL11A"
df[1,1]=0
df[2,1]=1
df[1,2]=0
df[2,2]=1
df[1,4]=0
df[2,4]=1


tableA<-text_transform(gt(df),
                       locations = cells_body(columns = Forward),
                       fn = function(x) {
                         c(local_image(filename = "pbmc_modisco_type5ig_modiscoreport/trimmed_logos/pos_patterns.pattern_0.cwm.fwd.png"),
                           local_image(filename = "pbmc_modisco_type5ig_modiscoreport/trimmed_logos/pos_patterns.pattern_1.cwm.fwd.png"))
                       }
)%>%text_transform(locations = cells_body(columns = Reverse),
                   fn = function(x) {
                     c(local_image(filename = "pbmc_modisco_type5ig_modiscoreport/trimmed_logos/pos_patterns.pattern_0.cwm.rev.png"),
                       local_image(filename = "pbmc_modisco_type5ig_modiscoreport/trimmed_logos/pos_patterns.pattern_1.cwm.rev.png"))
                   }
)%>%text_transform(locations = cells_body(columns = Motif),
                   fn = function(x) {
                     c(local_image(filename = "pbmc_modisco_type5ig_modiscoreport/MA0081.3.png"),
                       local_image(filename = "pbmc_modisco_type5ig_modiscoreport/MA2324.1.png"))
                   }
)%>%cols_align(align = "center", columns = c("Forward","Reverse","Matched TF","Motif","Q value"))%>%tab_header(
  title = md("B cell (Cluster 5)")
)%>%tab_options(table.font.size = 14,table.font.weight="bold",heading.title.font.weight="bold",column_labels.font.weight = 'bold')



gtsave(tableA, "Fig2_clus5.png")



########################

df=data.frame(matrix(ncol = 5,nrow = 2))
colnames(df)=c("Forward","Reverse","Matched TF","Motif","Q value")
df[1,5]=0.000363
df[2,5]=0.003558
df[1,3]="Spi1"
df[2,3]="CEBPA"
df[1,1]=0
df[2,1]=1
df[1,2]=0
df[2,2]=1
df[1,4]=0
df[2,4]=1


tableA<-text_transform(gt(df),
                       locations = cells_body(columns = Forward),
                       fn = function(x) {
                         c(local_image(filename = "pbmc_modisco_type6ig_modiscoreport/trimmed_logos/pos_patterns.pattern_0.cwm.fwd.png"),
                           local_image(filename = "pbmc_modisco_type6ig_modiscoreport/trimmed_logos/pos_patterns.pattern_3.cwm.fwd.png"))
                       }
)%>%text_transform(locations = cells_body(columns = Reverse),
                   fn = function(x) {
                     c(local_image(filename = "pbmc_modisco_type6ig_modiscoreport/trimmed_logos/pos_patterns.pattern_0.cwm.rev.png"),
                       local_image(filename = "pbmc_modisco_type6ig_modiscoreport/trimmed_logos/pos_patterns.pattern_3.cwm.rev.png"))
                   }
)%>%text_transform(locations = cells_body(columns = Motif),
                   fn = function(x) {
                     c(local_image(filename = "pbmc_modisco_type6ig_modiscoreport/MA0080.7.png"),
                       local_image(filename = "pbmc_modisco_type6ig_modiscoreport/MA0102.5.png"))
                   }
)%>%cols_align(align = "center", columns = c("Forward","Reverse","Matched TF","Motif","Q value"))%>%tab_header(
  title = md("Macrophage (Cluster 6)")
)%>%tab_options(table.font.size = 14,table.font.weight="bold",heading.title.font.weight="bold",column_labels.font.weight = 'bold')



gtsave(tableA, "Fig2_clus6.png")



