library(ggplot2)


stro<-read.csv(paste0(working_dir_path,"/strogenelist_h.csv"),row.names=1)
stro<-stro[stro$padj<0.05,]

gseadiff<-read.csv(paste0(working_dir_path,"/wdiff.csv"))
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
g

ggsave(paste0(working_dir_path,"/ssGSEA.pdf"),width = 12, height = 6)
