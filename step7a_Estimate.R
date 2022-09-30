################################################################
### step7a: Estimate analysis between UHR and CHR group (Fig6A)
################################################################
rm(list=ls())
library(ggplot2)
library(estimate)
library(RColorBrewer);
library(scater)
mycol = colorRampPalette(brewer.pal(9, "Set1"))(9)

setwd('/local/yanzijun/CRU/NB_FN/')
SEQC.exp <- readRDS('../NB/data/GSE62564/EXP_HR.RDS')
TARGET.exp <- readRDS('../NB/data/TARGET/RNAseq/fromGDC/EXP_HR.RDS')

GSE_ID='GSE62564'
if(GSE_ID=='GSE62564'){
  EXP=SEQC.exp
  model.res <- readRDS('res_scMK/Lasso/model.res.RDS')
  RS <- model.res$clin_data
  write.table(EXP,'../NB/data/GSE62564/EXP_HR.txt',row.names = T,col.names = T,sep='\t',quote = F)
}else if(GSE_ID=='TARGET'){
  EXP=TARGET.exp
  model.res <- readRDS('res_scMK/Lasso/TARGET.res.RDS')
  RS <- model.res$clin_data
  write.table(EXP,'../NB/data/TARGET/RNAseq/fromGDC/EXP_HR.txt',row.names = T,col.names = T,sep='\t',quote = F)
}

RS$sampleID=rownames(RS)
RS$type='high risk'
RS$type[RS$rs < quantile(RS$rs,0.5)]='low risk'
RS$rank=NULL
head(RS)

### Estimate
if(GSE_ID=='GSE62564'){
  exp.txt <- '/local/yanzijun/CRU/NB/data/GSE62564/EXP_HR.txt'
  exp.gct <- '/local/yanzijun/CRU/NB/data/GSE62564/EXP_HR.gct'
  PL="illumina"
}else{
  exp.txt <- '/local/yanzijun/CRU/NB/data/TARGET/RNAseq/fromGDC/EXP_HR.txt'
  exp.gct <-'/local/yanzijun/CRU/NB/data/TARGET/RNAseq/fromGDC/EXP_HR.gct'
  PL="illumina"
}

score.gct <- paste('/local/yanzijun/CRU/NB_FN/res_scMK/Estimate/',GSE_ID,'/Estimate.gct',sep='')
score.txt <- paste('/local/yanzijun/CRU/NB_FN/res_scMK/Estimate/',GSE_ID,'/Estimate.txt',sep='')

filterCommonGenes(input.f=exp.txt, output.f=exp.gct, id="GeneSymbol")
estimateScore(input.ds=exp.gct, output.ds=score.gct, platform=PL)

scores=read.table(score.gct,skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
scores[1:2,]
write.table(scores,score.txt,sep='\t',col.names = T,row.names = T,quote = F)

### plot
scores.df <- as.data.frame(scores)
scores.df <- tibble::rownames_to_column(scores.df,'sampleID')

if(GSE_ID=='TARGET'){
  scores.df$sampleID <- apply(as.matrix(scores.df$sampleID),1,function(x) gsub('\\.','-',x))
}

## plot boxplot
boxplot.lst <- list() 
variables <- c('StromalScore','ImmuneScore','ESTIMATEScore')
for(i in 1:length(variables)){
  var=variables[i]
  print(var)
  df <- as.data.frame(merge(RS,scores.df[,c('sampleID',var)],by='sampleID'))
  head(df)
  colnames(df)[ncol(df)] <- 'score' 
  df <- na.omit(df)
  df$var=var
  df$type=factor(df$type,levels = c('low risk','high risk'))
  head(df)
  
  plots <- ggboxplot(df, x = 'type', y = "score",
                     color = 'var', fill = 'var', 
                     add=NULL)+
    stat_compare_means(label = "p.format")+
    theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1 ,colour = 'black'),
          legend.position = 'none')+
    scale_color_manual(values = 'black')+
    scale_fill_manual(values = mycol[i])+xlab('')+ylab(var) 
  boxplot.lst[[i]] <- plots 
}
boxplots.p <- ggarrange(plotlist = boxplot.lst,ncol = length(variables))
boxplots.p
ggsave(paste('FIG_scMK/Estimate_boxplot_',GSE_ID,'.pdf',sep=''),boxplots.p,width = 7,height = 3)