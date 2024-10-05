library(stringr)
library(dplyr)

enhancer_list<-readRDS("direct_net_enhancerlist_0.obj")
enha_mtx<-as.data.frame(cbind(enhancer_list$gene,enhancer_list$Peak2))

peaks<-read.table("pbmc_fold1/peaks.bed",sep="\t",header=F)
peaks$tag=paste(peaks$V1,peaks$V2,peaks$V3,sep="_")

genelist<-read.csv("pbmc_fold1/pair_300000.csv",header=F)
genelist_prom<-read.csv("pbmc_fold1/pair_promoter.csv",header=F)
genelist_prom_s=genelist_prom

maxlen=max(genelist$V3-genelist$V2)+1


peakenha_direct=as.data.frame(matrix(-1, ncol = maxlen, nrow = dim(genelist)[1]))
for(i in 1:dim(genelist)[1]){
  genename<-genelist$V1[i]
  tmplist=enhancer_list$Peak2[enhancer_list$gene==genename]
  tmplist_Importance=enhancer_list$Importance[enhancer_list$gene==genename]
  
  pstart=genelist$V2[i]+1
  pend=genelist$V3[i]+1
  pnum=pend-pstart+1
  prompos=genelist_prom_s$V2[i]+1
  ppos=prompos-pstart+1
  tmp=peaks$tag[pstart:pend]
  tmp=tmp[-c(ppos)]
  
  tmpdf=as.data.frame(cbind(tmplist,tmplist_Importance))
  colnames(tmpdf)=c("tag","score")
  tmpdf_d=as.data.frame(tmp)
  colnames(tmpdf_d)=c("tag")
  tmpdf_j<-left_join(tmpdf_d,tmpdf,by="tag")
  tmpdf_j$score[is.na(tmpdf_j$score)]=0
  
  peakenha_direct[i,1]=0
  if(pnum>=2){peakenha_direct[i,2:pnum]=tmpdf_j$score}
}

write.table(peakenha_direct,"peakenha_direct_0.txt",sep="\t",quote=F,col.names=F,row.names=F)
