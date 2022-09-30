################################################################
### step7b: CIBERSORT analysis between UHR and CHR group (Fig6B)
################################################################
rm(list=ls())
setwd('/local/yanzijun/CRU/NB_FN/')

source('/local/yanzijun/CRU/NB_FN/res_scMK/CIBERSORT/CIBERSORT.R')
get_CIBERSORT <- function(expFile,resFile,sigresFile){
  res_cibersort <- CIBERSORT('/local/yanzijun/CRU/NB_FN/res_scMK/CIBERSORT/LM22.txt',
                             expFile,
                             perm = 1000)
  res_cibersort <- as.data.frame(res_cibersort,stringsAsFactors = F)
  print(dim(res_cibersort))
  saveRDS(res_cibersort,resFile)
  
  sigres_cibersort <- res_cibersort[which(res_cibersort$`P-value` < 0.05),]
  print(dim(sigres_cibersort))
  saveRDS(sigres_cibersort,sigresFile)
  
  res <- list(res_cibersort=res_cibersort,sigres_cibersort=sigres_cibersort)
  return(res)
}


## make input: EXP
SEQC.exp <- readRDS('../NB/data/GSE62564/EXP_HR.RDS')
TARGET.exp <- readRDS('../NB/data/TARGET/RNAseq/fromGDC/EXP_HR.RDS')
LM22 <- read.table('/local/yanzijun/CRU/NB_FN/res_scMK/CIBERSORT/LM22.txt',header = T,sep='\t',stringsAsFactors = F)[,1]

GSE_ID='GSE62564'

if(GSE_ID=='GSE62564'){
  EXP=SEQC.exp
}else if(GSE_ID=='TARGET'){
  EXP=TARGET.exp
}

EXP=EXP[rownames(EXP) %in% LM22,]
EXP <- tibble::rownames_to_column(EXP,'symbol')
write.table(EXP,paste('res_scMK/CIBERSORT/',GSE_ID,'/EXP.txt',sep=''),row.names = F,col.names = T,sep='\t',quote = F)


## CIBERSORT
expFile=paste('/local/yanzijun/CRU/NB_FN/res_scMK/CIBERSORT/',GSE_ID,'/EXP.txt',sep='')
resFile=paste('/local/yanzijun/CRU/NB_FN/res_scMK/CIBERSORT/',GSE_ID,'/res_cibersort.RDS',sep='')
sigresFile=paste('/local/yanzijun/CRU/NB_FN/res_scMK/CIBERSORT/',GSE_ID,'/sigres_cibersort.RDS',sep='')
CIBER <- get_CIBERSORT(expFile,resFile,sigresFile)


##########
## 1. Fig6B: Boxplot of TICC of each celltype between HR and LR
##########
rm(list=ls())
library(ggplot2)
library(RColorBrewer);
library(reshape2)
mycol = colorRampPalette(brewer.pal(9, "Set1"))(9)

GSE_ID='GSE62564' #TARGET
CIBER <- readRDS(paste('res_scMK/CIBERSORT/',GSE_ID,'/res_cibersort.RDS',sep=''))
CIBER <- tibble::rownames_to_column(CIBER,'sampleID')

if(GSE_ID=='GSE62564'){
  model.res <- readRDS('res_scMK/Lasso/model.res.RDS')
}else if(GSE_ID=='TARGET'){
  model.res <- readRDS('res_scMK/Lasso/TARGET.res.RDS')
}

RS <- model.res$clin_data
RS$sampleID=rownames(RS)
RS$type='UHR'
RS$type[RS$rs <= quantile(RS$rs,0.5)]='CHR'
RS[1:2,]

df <- merge(CIBER[,1:23],RS[,c('sampleID','type')])
df[1:3,]
melt.df <- melt(df)
melt.df[1:4,]
colnames(melt.df)[3:4] <- c('CellType','Prop')
melt.df$type <- factor(melt.df$type,levels = c('CHR','UHR'))
p <- ggboxplot(melt.df, x="CellType", y="Prop", color = "type",
               palette = mycol[2:1], add = NULL)+
  stat_compare_means(aes(group=type),label = "p.signif",paired = F)+
  theme(axis.text.x = element_text(angle = 30,hjust = 1,vjust = 1))
p

ggsave(paste('FIG_scMK/CIBER_boxplot_',GSE_ID,'.pdf',sep=''),plot = last_plot(),
       width = 20,height = 10,units = 'cm')
