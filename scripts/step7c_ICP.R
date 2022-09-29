################################################################
### step7c: compare immune checkpoints genes between UHR and CHR group (Fig6C)
################################################################
rm(list=ls())
library(RColorBrewer);
library(ggplot2)
library(reshape2)
mycol = colorRampPalette(brewer.pal(9, "Set1"))(10)

setwd('/local/yanzijun/CRU/NB_FN/')
Receptor <- c('ADORA2A','BTLA','CD27','CD40','CTLA4','EDNRB','HAVCR2','ICOS','IL2RA','LAG3','PDCD1',
              'TIGIT','TLR4','TNFRSF14','TNFRSF18','TNFRSF4','TNFRSF9')
Ligand <- c('CCL5','CD40LG','CD70','CX3CL1','CX3CL10','CX3CL9','IFNG','IL10','IL12A','IL1B','TGFB1',
            'TNF','TNFSF4','TNFSF9','VEGFA','VEGFB')
Costim <- c('CD28','CD80','ICOSLG')
Adhension <- c('ICAM1','ITGB2','SELP')
Antigen <- paste('HLA-',c('A','B','C','DMA','DMB','DOA','DOB','DPA1','DPB1',
                          'DQA1','DQA2','DQB1','DQB2','DRA','DRB1','DRB5','E','F','G'),sep='')
Coinhi <- c('BTN3A1','BTN3A2','CD247','CD246','PDCD1LG2','SLAMF7','VTCN1')
Other <- c('ENTPD1','GZMA','HMGB1','IDO1','PRF1')


GSE_ID='GSE62564'
if(GSE_ID=='GSE62564'){
  EXP <- readRDS('/local/yanzijun/CRU/NB/data/GSE62564/EXP_HR.RDS')
  model.res <- readRDS('res_scMK/Lasso/model.res.RDS')
}else if(GSE_ID=='TARGET'){
  EXP <- readRDS('/local/yanzijun/CRU/NB/data/TARGET/RNAseq/fromGDC/EXP_HR.RDS')
  model.res <- readRDS('res_scMK/Lasso/TARGET.res.RDS')
}

## plot 
Lst <- list(Receptor=Receptor,Ligand=Ligand,Costimulator=Costim,Adhension,Antigen=Antigen,
                Coinhibitor=Coinhi,Other=Other)

for(i in 1:length(Lst)){
  type=names(Lst)[i]
  geneLst=Lst[[i]]
  print(length(geneLst))
  
  sub.EXP <- as.data.frame(t(EXP[rownames(EXP) %in% geneLst,]))
  dim(sub.EXP)
  sub.EXP  <- tibble::rownames_to_column(sub.EXP,'sampleID')
  
  RS <- model.res$clin_data
  RS$sampleID=rownames(RS)
  RS$type='high risk'
  RS$type[RS$rs < quantile(RS$rs,0.5)]='low risk'
  RS[1:2,]
  
  
  df <- merge(sub.EXP,RS[,c('sampleID','type')])
  df$type=factor(df$type,levels = c('low risk','high risk'))
  df[1:3,]
  melt.df <- reshape2::melt(df)
  melt.df[1:4,]
  colnames(melt.df)[3:4] <- c('gene','exp')
  
  p <- ggboxplot(melt.df, x="gene", y="exp", color = "black",
                 palette = mycol[2:1], #add = "jitter"
                 add = FALSE,fill = "type")+
    stat_compare_means(aes(group=type),label = "p.signif",paired = F)+
    theme(axis.text = element_text(size=8,colour = 'black'),
          axis.text.x = element_text(angle = 30,hjust = 1,vjust = 1,size = 8),
          legend.position = 'none')+labs(y=type)
  print(p)
  
  if(length(geneLst)<10){
    width=5
  }else{
    width=10
  }
  ggsave(paste('FIG_scMK/ICP_',type,'_',GSE_ID,'.pdf',sep=''),plot = last_plot(),
         width = width,height = 5,units = 'cm')
}
