################################################################
### step5: Assocation between risk score and different clinical features (Fig4/sFig4)
################################################################
rm(list=ls())
library(ggplot2)
library(dplyr)
library("survival")
.libPaths(c("/local/yzj/R/x86_64-pc-linux-gnu-library/4.0" ,
            "/local/yanzijun/R/x86_64-pc-linux-gnu-library/4.0" ))
library("survminer")
library(survivalROC)
library(RColorBrewer);
library(scater)
library(Rmisc)
mycol = colorRampPalette(brewer.pal(9, "Set1"))(9)

setwd('/local/yanzijun/CRU/NB_FN/')
SEQC.clin <-  read.csv('../NB/data/GSE62564/clin.info_new.csv')
TARGET.clin <- read.csv('../NB/data/TARGET/RNAseq/fromGDC/clin/clin.info_new.csv')

### preprocess clinical information
SEQC.clin$age_b <- '>5'
SEQC.clin$age_b[SEQC.clin$age<= 5] <- '<=5'
SEQC.clin$age_b<- factor(SEQC.clin$age_b,levels = c('<=5','>5'))

SEQC.clin$sex_raw[SEQC.clin$sex_raw=='F']='Female'
SEQC.clin$sex_raw[SEQC.clin$sex_raw=='M']='Male'
SEQC.clin$sex_raw <- factor(SEQC.clin$sex_raw,levels = c('Female','Male'))
table(SEQC.clin$sex_raw)

SEQC.clin$mycn[SEQC.clin$mycn==0]='nonAmp'
SEQC.clin$mycn[SEQC.clin$mycn==1]='Amp'
SEQC.clin$mycn  <- factor(SEQC.clin$mycn ,levels = c('nonAmp','Amp'))

SEQC.clin$stage_b <- 'Stage 4'
SEQC.clin$stage_b[SEQC.clin$stage_raw %in% c('1','2','3')]='Stage 1_3'
SEQC.clin$stage_b<-factor(SEQC.clin$stage_b,levels = c('Stage 1_3','Stage 4'))
table(SEQC.clin$stage_b)

SEQC.clin$OS.status  <- factor(SEQC.clin$OS.status ,levels = c(0,1))
SEQC.clin$EFS.status  <- factor(SEQC.clin$EFS.status ,levels = c(0,1))


TARGET.clin$age_b <- '>5'
TARGET.clin$age_b[TARGET.clin$age<= 5] <- '<=5'
TARGET.clin$age_b<- factor(TARGET.clin$age_b,levels = c('<=5','>5'))
table(TARGET.clin$age_b)

TARGET.clin$sex_raw <- factor(TARGET.clin$sex_raw,levels = c('Female','Male'))

TARGET.clin$mycn[TARGET.clin$mycn==0]='nonAmp'
TARGET.clin$mycn[TARGET.clin$mycn==1]='Amp'
TARGET.clin$mycn  <- factor(TARGET.clin$mycn ,levels = c('nonAmp','Amp'))

TARGET.clin$stage_b <- 'Stage 4'
TARGET.clin$stage_b[TARGET.clin$stage_raw %in% c('Stage 3','Stage 2b')] <- 'Stage 1_3'
TARGET.clin$stage_b <- factor(TARGET.clin$stage_b,levels = c('Stage 1_3','Stage 4'))
table(TARGET.clin$stage_b)

TARGET.clin$OS.status  <- factor(TARGET.clin$OS.status ,levels = c(0,1))
TARGET.clin$EFS.status  <- '0'
TARGET.clin$EFS.status[TARGET.clin$EFS_raw.status=='Death'] <- 1
TARGET.clin$EFS.status  <- factor(TARGET.clin$OS.status ,levels = c(0,1))
table(TARGET.clin$EFS.status)

saveRDS(SEQC.clin,'../NB/data/GSE62564/clin.info_new_binary.RDS')
saveRDS(TARGET.clin,'../NB/data/TARGET/RNAseq/fromGDC/clin/clin.info_new_binary.RDS')


datSet='training'
if(datSet=='training'){
  clin=SEQC.clin
  model.res <- readRDS('res_scMK/Lasso/model.res.RDS')
  RS <- model.res$clin_data
}else if(datSet=='TARGET'){
  clin=TARGET.clin
  model.res <- readRDS('res_scMK/Lasso/TARGET.res.RDS')
  RS <- model.res$clin_data
}else if(datSet=='All'){
  clin=rbind(SEQC.clin[,c('sampleID','age_b','sex_raw','mycn','OS.status','stage_b','EFS.status')],
             TARGET.clin[,c('sampleID','age_b','sex_raw','mycn','OS.status','stage_b','EFS.status')])
  model2 <- readRDS('res_scMK/Lasso/TARGET.res.RDS')
  model1 <- readRDS('res_scMK/Lasso/model.res.RDS')
  RS <- rbind(model1$clin_data,model2$clin_data)
}


RS$sampleID=rownames(RS)
RS$type='UHR'
RS$type[RS$rs < quantile(RS$rs,0.5)]='CHR'
RS$type <- factor(RS$type,levels = c('UHR','CHR'))
RS$rank=NULL
head(RS)

## plot boxplot
boxplot.lst <- list()
pdf(paste('FIG_scMK/clinRS_KMplot_',datSet,'.pdf',sep=''),width = 10,height = 5) 
variables <- c('age_b','sex_raw','mycn','OS.status','stage_b','EFS.status')

for(i in 1:length(variables)){
  var=variables[i]
  print(var)
  df <- as.data.frame(merge(RS,clin[,c('sampleID',var)],by='sampleID'))
  dim(df)
  colnames(df)[ncol(df)] <- 'group' 
  df <- na.omit(df)
  df$var=var
  head(df)
  
  plots <- ggboxplot(df, x = 'group', y = "rs",
                     color = 'var', fill = 'var', 
                     add=NULL)+
    stat_compare_means(label = "p.format")+
    theme(#axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1 ,colour = 'black'),
      legend.position = 'none')+
    scale_color_manual(values = 'black')+
    scale_fill_manual(values = mycol[i])+xlab(var)+ylab('Risk score') 
  boxplot.lst[[i]] <- plots 
  
  ## KM curve
  kmplot.lst <- c()
  groups <- levels(df$group)
  for(j in 1:length(groups)){
    g=groups[j]
    sub.df <- df[df$group==g,]
    if(nrow(sub.df)>1){
      fit=survfit(Surv(time,status)~type, data = sub.df)
      p <- ggsurvplot(fit, pval = TRUE,palette ='Set1',
                      ggtheme = theme_survminer(base_size = 10),
                      title=paste(g,': KM of ',datSet,sep='')) ## days
      #surv_pvalue(fit)
      kmplot.lst[j] <- p
    }
  }
  if(length(groups)==2){
    multiplot(kmplot.lst[1],kmplot.lst[2], cols=2)
  }else if(length(groups)==3){
    multiplot(kmplot.lst[1],kmplot.lst[2],kmplot.lst[3], cols=3)
  }else if(length(groups)==4){
    multiplot(kmplot.lst[1],kmplot.lst[2],kmplot.lst[3],kmplot.lst[4],cols=4)
  }
}

boxplots.p <- ggarrange(plotlist = boxplot.lst,ncol = length(variables))
ggsave(paste('FIG_scMK/clinRS_boxplot_',datSet,'.pdf',sep=''),boxplots.p,width = 12,height = 3)
dev.off()


## varibale with p<0.05: plot KM curve with risk table
variable <- c('age_b','sex_raw','stage_b')
groups=c('<=5','Male','Stage 4')

for(i in 1:length(variable)){
  var=variable[i]
  group=groups[i]
  df <- as.data.frame(merge(RS,clin[,c('sampleID',var)],by='sampleID'))
  colnames(df)[ncol(df)] <- 'group' 
  df <- na.omit(df)
  head(df)
  sub.df <- df[df$group==group,]
  fit=survfit(Surv(time,status)~type, data = sub.df)
  p <- ggsurvplot(fit, pval = TRUE,palette ='Set1',
                  ggtheme = theme_survminer(base_size = 8),
                  risk.table = TRUE) ## days
  p
  ggsave(paste('FIG_scMK/clinRS_KMplot_',var,'_',datSet,'_table.pdf',sep=''),plot = last_plot(),
         width = 4.5,height = 1.5)
}

## sankey plot
library(ggalluvial)
library(tidyverse)
library(hrbrthemes)
library(wesanderson)

#input data
df <- as.data.frame(merge(RS,clin[clin$sampleID %in% RS$sampleID,c('sampleID','mycn','OS.status')],by='sampleID'))
dim(df)
df$mycn <- factor(df$mycn,levels = c('Amp','nonAmp'))
df[1:2,]
tmp <- df[,c('type','mycn','OS.status')]
tmp1 <- table(tmp$type,tmp$mycn,tmp$OS.status)

## chisq.test
tableR <- as.matrix(table(df$type,df$mycn,df$OS.status))
chisq.test(tableR)

## plot
plot.df <- as.data.frame(tmp1)
head(plot.df)
plot.df$Var3 <- as.character(plot.df$Var3)
plot.df$Var3 <- factor(plot.df$Var3,levels = c('1','0'))

is_alluvia_form(plot.df, axes = 1:3, silent = TRUE)
alluvia.p <- ggplot(plot.df,
                    aes(y = Freq, axis1 = Var1, axis2 = Var2)) +
  geom_alluvium(aes(fill = Var3), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Risk Group", "MYCN"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = 'Set1') +
  ggtitle("distribution of patients, by risk, status and mycn status")+
  theme_classic2()
alluvia.p
pdf(paste('FIG_scMK/clinRS_alluvia_',datSet,'_status.pdf',sep=''),width = 5.1,height = 4)
alluvia.p
dev.off()