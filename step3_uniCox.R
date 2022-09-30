################################################################
### step3: get prognostic genes from SEQC HR patients 
################################################################
############
## 1. processing SEQC/TARGET RNA-seq data and clinical information
############
rm(list=ls())
library(GEOquery)

### 1.1 SEQC cohort (GSE62564 re-anlaysis for GSE49711)
GSE_ID='GSE62564'
setwd(paste('/local/yanzijun/CRU/NB/data/',GSE_ID,sep=''))
EXP <- read.table('GSE62564_SEQC_NB_RNA-Seq_log2RPM.txt',sep='\t',header = T)
gene2NM <- readRDS('/local/yanzijun/public/annotation/hg19/gene2NM.hg19.refGene.RDS')
dim(gene2NM)
ids=gene2NM
colnames(ids) <- c('probe_id','symbol')
head(ids)
ids=na.omit(ids)
dim(ids)

median=apply(EXP[,2:ncol(EXP)],1,median)
df_median <- data.frame(probe_id=EXP$RefSeqID,median=median)

ids <- merge(ids,df_median,by='probe_id')
ids=ids[order(ids$symbol,ids$median,decreasing = T),]
ids=ids[!duplicated(ids$symbol),1:2]
colnames(ids)[1] <- 'RefSeqID'

EXP=merge(ids,EXP,by='RefSeqID')
EXP$RefSeqID=EXP$gene <- NULL

dim(EXP)
dim(na.omit(EXP))
saveRDS(EXP,'EXP.RDS')

## clinical info
library(dplyr)
clin.info1 <- read.csv('/local/yanzijun/CRU/NB/data/GSE62564/clin.info.csv')
clin.df1 <- clin.info1[,c('title','age.ch1','Sex.ch1','high.risk.ch1',
                          'os.bin.ch1','os.day.ch1','efs.bin.ch1','efs.day.ch1')]
colnames(clin.df1) <- c('sampleID','age','sex_raw','risk_raw','OS.status','OS.time','EFS.status','EFS.time')
clin.df1$sampleID <- apply(as.matrix(clin.df1$sampleID),1,function(x) unlist(strsplit(x,' '))[1])
clin.df1$age <- round((clin.df1$age)/365,2)

clin.info2 <- read.csv('/local/yanzijun/CRU/NB/data/SEQC/GSE49711/clin.info.csv')
clin.df2 <- clin.info2[,c('title','inss.stage.ch1','mycn.status.ch1')]
colnames(clin.df2) <- c('sampleID','stage_raw','mycn')

clin.df <- merge(clin.df1,clin.df2)
col.order <- c('sampleID','age','sex_raw','mycn','risk_raw','stage_raw','OS.status','OS.time','EFS.status','EFS.time')
clin.df <- select(clin.df,col.order)

## to binary
table(clin.df$sex_raw)
clin.df$sex <- NA
clin.df$sex[clin.df$sex_raw=='F']<- 0
clin.df$sex[clin.df$sex_raw=='M']<- 1
table(clin.df$sex)

table(clin.df$risk_raw)
clin.df$risk <- 0
clin.df$risk[clin.df$risk_raw=='HR']<- 1
table(clin.df$risk)
clin.df[1:2,]
write.csv(clin.df,'/local/yanzijun/CRU/NB/data/SEQC/GSE49711/clin.info_new.csv')
write.csv(clin.df,'/local/yanzijun/CRU/NB/data/SEQC/GSE62564/clin.info_new.csv')

## to surv info
OS.info <- clin.df[,c('sampleID','OS.status','OS.time')]
EFS.info <- clin.df[,c('sampleID','EFS.status','EFS.time')]
colnames(OS.info)[2:3]=colnames(EFS.info)[2:3]=c('status','time')
surv.lst=list(OS=OS.info,EFS=EFS.info)
saveRDS(surv.lst,'/local/yanzijun/CRU/NB/data/SEQC/GSE49711/surv.info.RDS')
saveRDS(surv.lst,'/local/yanzijun/CRU/NB/data/SEQC/GSE62564/surv.info.RDS')

## HR expression
EXP <- readRDS('EXP.RDS')
GSM.HR <- clin.df$sampleID[clin.df$risk=='1']
GSM.LR <- clin.df$sampleID[clin.df$risk=='0']

rownames(EXP) <- EXP$symbol;EXP$symbol=NULL
EXP_HR <- select(EXP,GSM.HR)
saveRDS(EXP_HR,'EXP_HR.RDS')
write.csv(EXP_HR,'SEQC_HRNB_expression.csv')


## 1.2 TARGET cohort
rm(list=ls())
library(xlsx)
setwd('/local/yanzijun/CRU/NB/data/TARGET/RNAseq')

clinTARGET_dis <- read.xlsx('fromTARGET/clin/TARGET_NBL_ClinicalData_Discovery.xlsx',sheetIndex = 1)
clinTARGET_val <- read.xlsx('fromTARGET/clin/TARGET_NBL_ClinicalData_Validation.xlsx',sheetIndex = 1)

clinGDC <- read.table('fromGDC/clin/gdc_sample_sheet.2022-06-02.tsv',header = T,sep='\t')
length(intersect(clinGDC$Case.ID,clinTARGET_dis$TARGET.USI))
length(intersect(clinGDC$Case.ID,clinTARGET_val$TARGET.USI))

ol <- intersect(clinGDC$Case.ID,clinTARGET_dis$TARGET.USI)
clin.df <- clinTARGET_dis[clinTARGET_dis$TARGET.USI %in% ol,
                          c('TARGET.USI','Age.at.Diagnosis.in.Days','Gender','MYCN.status',
                            'INSS.Stage','COG.Risk.Group',
                            'Vital.Status','Overall.Survival.Time.in.Days',
                            'First.Event','Event.Free.Survival.Time.in.Days'
                          )]
dim(clin.df)
clin.df[1:2,]
colnames(clin.df) <- c('sampleID','age','sex_raw','mycn_raw','stage_raw','risk_raw',
                       'OS_raw.status','OS.time','EFS_raw.status','EFS.time')
clin.df[1:2,]
clin.df$age <- round(clin.df$age/365,2)

table(clin.df$sex_raw)
clin.df$sex <- NA
clin.df$sex[clin.df$sex_raw=='Male']<- 1
clin.df$sex[clin.df$sex_raw=='Female']<- 0
table(clin.df$sex)

table(clin.df$mycn_raw)
clin.df$mycn <- NA
clin.df$mycn[clin.df$mycn_raw=='Amplified']<- 1
clin.df$mycn[clin.df$mycn_raw=='Not Amplified']<- 0
table(clin.df$mycn)

table(clin.df$risk_raw)
clin.df$risk <- NA
clin.df$risk[clin.df$risk_raw=='High Risk']<- 1
clin.df$risk[clin.df$risk_raw=='Low Risk']<- 0
table(clin.df$risk)
clin.df[1:2,]

table(clin.df$OS_raw.status)
clin.df$OS.status <- NA
clin.df$OS.status[clin.df$OS_raw.status=='Dead']<- 1
clin.df$OS.status[clin.df$OS_raw.status=='Alive']<- 0
table(clin.df$OS.status)
write.csv(clin.df,'fromGDC/clin/clin.info_new.csv',row.names = F)
write.csv(clin.df,'fromTARGET/clin/clin.info_new.csv',row.names = F)

## to surv info
OS.info <- clin.df[,c('sampleID','OS.status','OS.time')]
colnames(OS.info)[2:3]=c('status','time')
surv.lst=list(OS=OS.info)
saveRDS(surv.lst,'fromGDC/clin/surv.info.RDS')
saveRDS(surv.lst,'fromTARGET/clin/surv.info.RDS')

## exp
tmp <- clinGDC[clinGDC$Case.ID %in%ol,c('File.ID','Case.ID')]
name2ID=tmp$Case.ID
names(name2ID) <- tmp$File.ID

files <- names(name2ID)
head(files)

i=1
file=files[i]
exp_file <- list.files(paste('fromGDC/exp/',file,sep=''))
EXP <- read.table(paste('fromGDC/exp/',file,'/',exp_file,sep=''),header = T,skip = 1,sep='\t')
EXP <- EXP[-c(1:4),c('gene_name','fpkm_uq_unstranded')]
colnames(EXP)[2] <- name2ID[file]
EXP <- EXP[!duplicated(EXP$gene_name),]
dim(EXP)
EXP[1:2,]

for(i in 2:length(files)){
  print(i)
  file=files[i]
  exp_file <- list.files(paste('fromGDC/exp/',file,sep=''))
  df <- read.table(paste('fromGDC/exp/',file,'/',exp_file,sep=''),header = T,skip = 1,sep='\t')
  df <- df[-c(1:4),c('gene_name','fpkm_uq_unstranded')]
  colnames(df)[2] <- name2ID[file]
  df <- df[!duplicated(df$gene_name),]
  EXP <- merge(EXP,df,by='gene_name',all=FALSE)
}
colnames(EXP)[1] <- 'symbol'
tmp <- EXP
rownames(tmp) <- tmp$symbol;tmp$symbol=NULL
max(tmp)
if(max(tmp)>50){
  tmp <- log(tmp+1,2)
}
tmp <- tibble::rownames_to_column(tmp,'symbol')

saveRDS(tmp,'fromGDC/EXP.RDS')
saveRDS(tmp,'fromTARGET/EXP.RDS')

write.csv(tmp,'fromGDC/TARGET_NB_expression.csv',col.names = TRUE,row.names = F)
write.csv(tmp,'fromTARGET/TARGET_NB_expression.csv',col.names = TRUE,row.names = F)

## HR expression
tmp <- readRDS('fromTARGET/EXP.RDS')
clin.df <- read.csv('fromTARGET/clin/clin.info_new.csv')
HR.sample <- clin.df$sampleID[clin.df$risk=='1']

rownames(tmp)=tmp$symbol;tmp$symbol=NULL
EXP_HR <- tmp[,colnames(tmp)%in%HR.sample]
saveRDS(EXP_HR,'fromTARGET/EXP_HR.RDS')
saveRDS(EXP_HR,'fromGDC/EXP_HR.RDS')

write.csv(EXP_HR,'fromTARGET/TARGET_HRNB_expression.csv')
write.csv(EXP_HR,'fromGDC/TARGET_HRNB_expression.csv')

############
## 2. uniCox
############
rm(list=ls())
library("survival")

## uni Cox（OS) function
get_Unires <- function(Cox_df){
  covariates <- colnames(Cox_df)[4:ncol(Cox_df)]
  univ_formulas <- sapply(covariates,
                          function(x) as.formula(paste('Surv(time, status)~', x)))
  univ_models <- lapply(univ_formulas, function(x){coxph(x, data = Cox_df)})
  univ_results <- lapply(univ_models,
                         function(x){ 
                           x <- summary(x)
                           p.value<-signif(x$sctest["pvalue"], digits=2)
                           sc.test<-signif(x$sctest["test"], digits=2)
                           beta<-signif(x$coef[1], digits=2);#coeficient beta
                           HR <-signif(x$coef[2], digits=2);#exp(beta)
                           HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                           HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                           HR.confint <- paste0(HR, " (", 
                                                HR.confint.lower, "-", HR.confint.upper, ")")
                           # res<-c(beta, HR.confint, sc.test, p.value)
                           # names(res)<-c("beta", "HR (95% CI for HR)", "sc.test", "p.value")
                           res<-c(HR, HR.confint.lower, HR.confint.upper, HR.confint, p.value)
                           names(res)<-c( 'HR','lower','upper',"HR (95% CI for HR)", "p.value")
                           return(res)
                         })
  
  res <- t(as.data.frame(univ_results, check.names = FALSE))
  res <- as.data.frame(res,stringsAsFactors = F)
  return(res)
}


## survival file
GSE_ID='GSE62564' 
setwd(paste('/local/yanzijun/CRU/NB/data/',GSE_ID,sep=''))

if(GSE_ID=='GSE62564'){
  Surv.lst <- readRDS('surv.info.RDS')
}else{
  Surv.lst <- readRDS('RNAseq/fromTARGET/clin/surv.info.RDS')
}

type='OS'
if(type=='OS'){
  Surv <- Surv.lst$OS
}else if(type=='EFS'){
  Surv <- Surv.lst$EFS
}else if(type=='PFS'){
  Surv <- Surv.lst$DFS
}


#### gene exp file
## only for HighRisk patients
if(GSE_ID=='GSE62564'){
  EXP <-  readRDS('EXP_HR.RDS')
}else{
  EXP <- readRDS('RNAseq/fromGDC/EXP_HR.RDS')
}


ol.samples <- intersect(colnames(EXP),Surv$sampleID)
Surv <- Surv[Surv$sampleID %in% ol.samples,]
dim(Surv)
EXP <- dplyr::select(EXP,Surv$sampleID)
dim(EXP)
EXP[1:4,1:3]

geneLst <- rownames(EXP)
sub.EXP <- EXP[rownames(EXP) %in% geneLst,]
sub.EXP <- sub.EXP[-c(grep('-',rownames(sub.EXP)),
                      grep('\\/',rownames(sub.EXP)),
                      grep('@',rownames(sub.EXP)),
                      grep('\\.',rownames(sub.EXP)),
                      grep('^\\d',rownames(sub.EXP)),
                      grep('NOP5',rownames(sub.EXP))),]

dim(sub.EXP)

sub.EXP <- as.data.frame(t(sub.EXP))
print(all(rownames(sub.EXP)==Surv$ID))

Cox_df <- cbind(Surv,sub.EXP)
dim(Cox_df)

Unires <- get_Unires(Cox_df)
print(Unires[1:2,])

sig.Unires <- Unires[Unires$p.value<0.05,]
dim(sig.Unires)

system(paste('mkdir -p /local/yanzijun/CRU/NB_FN/res/GeneSet/UniCox/',GSE_ID,sep=''))
saveRDS(Unires,paste('/local/yanzijun/CRU/NB_FN/res/GeneSet/UniCox/',GSE_ID,'/Unires_',type,'_HR.RDS',sep=""))
saveRDS(sig.Unires,paste('/local/yanzijun/CRU/NB_FN/res/GeneSet/UniCox/',GSE_ID,'/sigUnires_',type,'_HR.RDS',sep=''))

favorG=rownames(sig.Unires)[sig.Unires$HR<1]
unfavorG=rownames(sig.Unires)[sig.Unires$HR>1]
progG.lst=list(favorG=favorG,unfavorG=unfavorG)
saveRDS(progG.lst,paste('/local/yanzijun/CRU/NB_FN/res/GeneSet/UniCox/',GSE_ID,'/progGeneLst_',type,'_HR.RDS',sep=''))