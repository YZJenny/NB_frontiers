################################################################
### step10: univariable and multivariable analysis in clinical features and MMGS (Table1-2)
################################################################
rm(list=ls())
library(reshape2)
library(ggplot2)
library(forestplot)
library(plotrix)
library(RColorBrewer)
library("survival")
getPalette = colorRampPalette(brewer.pal(9, "Set1"))(9)

get_Unires <- function(Cox_df){
  covariates <- colnames(Cox_df)[4:ncol(Cox_df)]
  univ_formulas <- sapply(covariates,
                          function(x) as.formula(paste('Surv(time, status)~', x)))
  univ_models <- lapply(univ_formulas, function(x){coxph(x, data = Cox_df)})
  univ_results <- lapply(univ_models,
                         function(x){ 
                           x <- summary(x)
                           p.value<-signif(x$sctest["pvalue"], digits=4)
                           sc.test<-signif(x$sctest["test"], digits=4)
                           beta<-signif(x$coef[1], digits=4);#coeficient beta
                           HR <-signif(x$coef[2], digits=4);#exp(beta)
                           HR.confint.lower <- signif(x$conf.int[,"lower .95"], 4)
                           HR.confint.upper <- signif(x$conf.int[,"upper .95"],4)
                           HR.confint <- paste0(HR, " (", 
                                                HR.confint.lower, "-", HR.confint.upper, ")")
                           res<-c(HR, HR.confint.lower, HR.confint.upper, HR.confint, p.value)
                           names(res)<-c( 'HR','lower','upper',"HR (95% CI for HR)", "p.value")
                           return(res)
                         })
  res <- t(as.data.frame(univ_results, check.names = FALSE))
  res <- as.data.frame(res,stringsAsFactors = F)
  res$Gene <- rownames(res)
  res <- as.matrix(res)
  res <- rbind(colnames(res),res)
  return(res)
  
}

#### survival file
GSE_ID='GSE62564' #'GSE62564',TARGET
setwd(paste('/local/yanzijun/CRU/NB/data/',GSE_ID,sep=''))

#### candiGene file
if(GSE_ID=='GSE62564'){
  EXP <- readRDS('/local/yanzijun/CRU/NB/data/GSE62564/EXP_HR.RDS')
  Surv.lst <- readRDS('/local/yanzijun/CRU/NB/data/GSE62564/surv.info.RDS')
  clin.info <-  read.csv('/local/yanzijun/CRU/NB/data/GSE62564/clin.info_new.csv')
  model.res <- readRDS('/local/yanzijun/CRU/NB_FN/res_scMK/Lasso/model.res.RDS')
}else{
  EXP <- readRDS('/local/yanzijun/CRU/NB/data/TARGET/RNAseq/fromGDC/EXP_HR.RDS')
  Surv.lst <- readRDS('/local/yanzijun/CRU/NB/data/TARGET/RNAseq/fromTARGET/clin/surv.info.RDS')
  clin.info <- read.csv('/local/yanzijun/CRU/NB/data/TARGET/RNAseq/fromGDC/clin/clin.info_new.csv')
  model.res <- readRDS('/local/yanzijun/CRU/NB_FN/res_scMK/Lasso/TARGET.res.RDS')
}

## 年龄/stage变成二分类
clin.info[1:2,]
clin.info$Age <- 1
clin.info$Age[clin.info$age < 5] <- 0
table(clin.info$Age)

clin.info$Stage <- 1
if(GSE_ID=='GSE62564'){
  clin.info$Stage[clin.info$stage_raw <= 2] <- 0
}else{
  clin.info$Stage[clin.info$stage_raw %in% c('Stage 2b')] <- 0
}

table(clin.info$Stage)

type='OS'
if(type=='OS'){
  Surv <- Surv.lst$OS
}else if(type=='EFS'){
  Surv <- Surv.lst$EFS
}else if(type=='PFS'){
  Surv <- Surv.lst$DFS
}


########
## 1. univariable analysis
########
#### candiGene file
RS <- model.res$clin_data
RS[1:2,]
RS <- tibble::rownames_to_column(RS,'sampleID')
RS$MMGS <- 1
RS$MMGS[RS$rs <= quantile(RS$rs,0.5)] <- 0
table(RS$MMGS)

clin.info <- clin.info[,c('sampleID','Age','sex','mycn')]
#clin.df <- merge(RS[,c('sampleID','status','time','rs')],clin.info)
clin.df <- merge(RS[,c('sampleID','status','time','MMGS')],clin.info)
clin.df[1:2,]
dim(clin.df)

Unires <- get_Unires(clin.df)
print(Unires)

########
## 2. multivaribale analysis
########
.libPaths(c("/local/yzj/R/x86_64-pc-linux-gnu-library/4.0",
            "/local/yanzijun/R/x86_64-pc-linux-gnu-library/4.0"))
library(survival)
library(survminer)
y='Surv(time, status)'
x <- paste(colnames(clin.df)[4:ncol(clin.df)],collapse = '+')
print(x)
surv <<- as.formula(paste(y, "~", x))
Model <- coxph(surv , data =  clin.df)
print(Model)
print(summary(Model))
