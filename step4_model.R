################################################################
### step4: construction model (Fig2/Fig3/sFig3)
################################################################
########################
### 1. get candidate gene set
########################
rm(list=ls())
## 1.1 marker gene set of malignant cells
MK <- read.csv("/local/yanzijun/CRU/NB_FN/res/GeneSet/scMK/markers_H.csv")
table(MK$cluster)
MK.T <- MK$gene[MK$cluster=='tumor']
print(length(MK.T)) 

## 1.2 METcor gene set
sigcor.df <- readRDS('/local/yanzijun/CRU/NB_FN/res/GeneSet/Cor/METcor.RDS')
print(nrow(sigcor.df)) 
## MK-METcor gene set
corMK.T <- intersect(MK.T,sigcor.df$gene)
print(length(corMK.T)) 

## 1.3 prognostic-related gene set
progG.lst <- readRDS('/local/yanzijun/CRU/NB_FN/res/GeneSet/UniCox/GSE62564/progGeneLst_OS_HR.RDS')
riskG <- c(progG.lst$unfavorG,progG.lst$favorG)


########################
### 2. training model
########################
library(glmnet)
library(dplyr)
library(RColorBrewer)
mycol = colorRampPalette(brewer.pal(9, "Set1"))(9)
setwd('/local/yanzijun/CRU/NB_FN/')

### 2.0 input genelist
geneLst <- intersect(corMK.T,riskG)
print(length(geneLst)) 

EXP <- readRDS('/local/yanzijun/CRU/NB/data/GSE62564/EXP_HR.RDS')
Surv.lst <- readRDS('/local/yanzijun/CRU/NB/data/GSE62564/surv.info.RDS')
datTraits <- Surv.lst$OS

## clinical time can't be 0
clin_info=datTraits
rownames(clin_info) <- clin_info$sampleID;clin_info$sampleID=NULL
clin_info  <- clin_info[clin_info$time>0,]
clin_info=as.matrix(clin_info)

ol.samples=intersect(colnames(EXP),rownames(clin_info))
clin_info=clin_info[rownames(clin_info) %in% ol.samples,]
datExpr <- as.data.frame(t(select(EXP,rownames(clin_info))))
X=datExpr[,which(colnames(datExpr)%in%geneLst)]
X=as.matrix(X)
dim(X)
print(all(rownames(clin_info)==rownames(X)))


### 2.1 feature selection: lasso_cox
set.seed(123)
fit = glmnet(X, clin_info, family = "cox",alpha=1)
cvfit = cv.glmnet(X, clin_info, family = "cox",nfolds=5,alpha=1)

print(fit$dev.ratio[which(fit$lambda==cvfit$lambda.min)])
coef.min = as.matrix(coef(cvfit, s = "lambda.min"))

lasso.risk_factor <- coef.min[which(coef.min[,1]!=0),1]
print(length(lasso.risk_factor)) 
print(lasso.risk_factor)

lasso.df <- as.data.frame(cbind(clin_info,X[,colnames(X)%in%names(lasso.risk_factor)]))
dim(lasso.df)

## sFig3A and B
pdf('FIG_scMK/lasso_cox.pdf',width = 6,height = 6)
plot(fit)
plot(cvfit)
dev.off()

## multi-stepwise
library(My.stepwise)
if(length(lasso.risk_factor)>10){
My.stepwise.coxph(Time = "time", Status = "status", variable.list = names(lasso.risk_factor),
                  data = lasso.df,sle = 0.25,sls = 0.25)
}


if(length(lasso.risk_factor)>10){
  stepwise.risk_factor <- c(0.02401,-0.07335,0.27887,0.27786,-0.34653 ,-0.10144)
  names(stepwise.risk_factor) <- c('NPW','MAPT','RAMP1','CDT1','C1QTNF4','MEG3')
  lasso.obj <- list(fit=fit,cvfit=cvfit,coef.min=coef.min,lasso.df=lasso.df,
                    lasso.risk_factor=lasso.risk_factor,stepwise.risk_factor=stepwise.risk_factor)
}else{
  lasso.obj <- list(fit=fit,cvfit=cvfit,coef.min=coef.min,lasso.df=lasso.df,
                    lasso.risk_factor=lasso.risk_factor)
}
saveRDS(lasso.obj,'res_scMK/Lasso/lassoCox_stepwise.RDS')


## Fig2A: barplot of coefficient
library(ggplot2)
lasso.obj <- readRDS('res_scMK/Lasso/lassoCox_stepwise.RDS')

risk_factor=lasso.obj$lasso.risk_factor
if(length(risk_factor)>10){
  risk_factor=lasso.obj$stepwise.risk_factor
}
gene <- names(risk_factor)
print(length(gene))

data <- data.frame(factor=factor(gene,levels=gene[order(risk_factor)]),coef=risk_factor)
data$group='Protective'
data$group[data$coef>0]='Unprotective'
p <- ggplot(data,aes(factor,coef,fill=group))+
  geom_bar(stat='identity',width = 0.9,colour='black')+
  coord_flip() +
  labs(x="",y="individual coefficients")+
  theme(panel.grid=element_blank(),
        panel.background = element_blank(),
        axis.ticks.y= element_blank())+ #delete backgroud
  guides(fill=FALSE)+ #delete legend
  theme(axis.text.x = element_text(size = 10,colour = 'black'),
        axis.text.y = element_text(size = 10,colour = 'black'), 
        axis.title.x = element_text(size = 10,colour = 'black'), 
        axis.title.y = element_text(size = 10,colour = 'black'))+
  scale_fill_manual(values = mycol[2:1])
p
ggsave(filename='FIG_scMK/Coef_barplot.pdf',plot=p,width=6,height=3)


### 2.2 construct model
X.model <- dplyr::select(as.data.frame(X),names(risk_factor))
dim(X.model)
print(all(colnames(X.model)==names(risk_factor)))

RS=c()
for(i in 1:nrow(X.model)){
  RS[i]=(risk_factor[1]*X.model[i,1])+(risk_factor[2]*X.model[i,2])+(risk_factor[3]*X.model[i,3])
  +(risk_factor[4]*X.model[i,4])+(risk_factor[5]*X.model[i,5])+(risk_factor[6]*X.model[i,6])
}

names(RS)=rownames(X.model)
saveRDS(RS,'res_scMK/Lasso/RS_training.RDS')

## Fig2B: plot RS
rs_value <- sort(RS)
left <- rs_value[1:floor(length(rs_value)/2)]
right <- rs_value[floor(length(rs_value)/2+1):length(rs_value)]
pdf('FIG_scMK/RS_training_plot.pdf',width=8,height=6)
plot(x=1:length(left),y=left,col=mycol[2],pch=19,xlim=c(0,length(RS)),ylim=c(min(rs_value),max(rs_value)),
     xlab="",ylab="",cex=1.2,cex.lab=1.3,cex.axis=1.2,xaxt='n')
par(new=TRUE)
plot(x=(length(left)+1):length(RS),y=right,col=mycol[1],pch=19,xlim=c(0,length(RS)),ylim=c(min(rs_value),max(rs_value)),
     xaxt='n',yaxt='n',xlab="",ylab="Risk Score",cex=1.2,cex.lab=1.3,cex.axis=1.2)
abline(v=round(length(RS)/2),col = "gray", lwd = 3.5, lty = 2,xaxt='n')
dev.off()


## Fig2B: plot follow time
print(all(names(RS)==rownames(clin_info)))
data=cbind(rs=RS,status=clin_info)
data=data[order(data[,1]),]
data=as.data.frame(cbind(data,rank=1:nrow(data)))

alive=data[which(data$status==0),]
dead=data[which(data$status==1),]
pdf('FIG_scMK/FollowTime_training_plot.pdf',width=8,height=6)
plot(x=dead$rank,y=dead$time,col=mycol[1],pch=17,xlim=c(0,nrow(data)),ylim=c(min(data$time),max(data$time)),
     xlab="",ylab="",cex=1.2,cex.lab=1.3,cex.axis=1.2)
par(new=TRUE)
plot(x=alive$rank,y=alive$time,col=mycol[2],pch=17,xlim=c(0,nrow(data)),ylim=c(min(data$time),max(data$time)),
     xaxt='n',yaxt='n',xlab="Patients cohort",ylab="Follow-up time(days)",cex=1.2,cex.lab=1.3,cex.axis=1.2)
abline(v=round(length(RS)/2)+0.5,col = "gray", lwd = 3.5, lty = 2)
dev.off()

## Fig2C: plot heatmap of 6 risk gene
library(pheatmap)
library(dplyr)

heatmap.df <- as.data.frame(t(X.model))
heatmap.df <- dplyr::select(heatmap.df,rownames(data))
print(all(colnames(heatmap.df)==rownames(data)))

if(ncol(heatmap.df)%%2==1){
  annotation_col = data.frame(Type = factor(c(rep("lowrisk", floor(length(RS)/2)),
                                              rep('highrisk',floor(length(RS)/2)+1))))
}else{
  annotation_col = data.frame(Type = factor(c(rep("lowrisk", floor(length(RS)/2)),
                                              rep('highrisk',floor(length(RS)/2)))))
}
rownames(annotation_col) = colnames(heatmap.df)
p.heatmap <- pheatmap(heatmap.df,cluster_rows =T,cluster_cols = F,annotation_col = annotation_col,
                      legend = TRUE, annotation_legend = TRUE, show_colnames =  F,
                      color= colorRampPalette(c("navy", "white", "firebrick3"))(50),border_color=NA,
                      annotation_colors =list(Type=c(lowrisk=mycol[2],highrisk=mycol[1])),
                      fontsize=10,fontsize_row=16)
ggsave(filename = 'FIG_scMK/Coef_heatmap.pdf',p.heatmap,width = 8,height = 3)

model.res <- list(RS.train=RS,X.model=X.model,clin_data=data,heatmap.df=heatmap.df)
saveRDS(model.res,'res_scMK/Lasso/model.res.RDS')


### Fig2D: survival analysis
setwd('/local/yanzijun/CRU/NB_FN/')
.libPaths(c("/local/yzj/R/x86_64-pc-linux-gnu-library/4.0",
            "/local/yanzijun/R/x86_64-pc-linux-gnu-library/4.0"))
library(ggplot2)
library(pheatmap)
library(dplyr)
library("survival")
library("survminer")
library(survivalROC)
model.res <- readRDS("res_scMK/Lasso/model.res.RDS")
RS=model.res$RS.train
clin_info=model.res$clin_data
clin_info$rank=NULL

cut_off=quantile(RS,0.5)
sur_df <- clin_info
sur_df$time=round(sur_df$time/365,3)
sur_df$surtype <- 'UHR'
sur_df$surtype[sur_df$rs < cut_off]='CHR'
sur_df$surtype <- factor(sur_df$surtype,levels = c('UHR','CHR'))

fit <- survfit(Surv(time =time, event =status) ~ surtype,data=sur_df)
p <- ggsurvplot(fit, pval = TRUE,palette ='Set1',risk.table = TRUE,
                ggtheme = theme_survminer(base_size = 10),
                title='Survival curve of training set') 
p
print(surv_pvalue(fit)$pval)

pdf('FIG_scMK/survCurve_training.pdf',width = 5,height = 5.5)
print(p)
dev.off()




####################################
#### 3.Validation model
####################################
rm(list=ls())
.libPaths(c("/local/yzj/R/x86_64-pc-linux-gnu-library/4.0",
            "/local/yanzijun/R/x86_64-pc-linux-gnu-library/4.0"))
library(ggplot2)
library(pheatmap)
library(dplyr)
library("survival")
library("survminer")
library(survivalROC)
library(RColorBrewer)
mycol = colorRampPalette(brewer.pal(9, "Set1"))(9)


lasso.res <- readRDS("res_scMK/Lasso/lassoCox_stepwise.RDS")
risk_factor=lasso.res$lasso.risk_factor
if(length(risk_factor)>10){
  risk_factor=lasso.res$stepwise.risk_factor
}
print(length(risk_factor))


validSet='TARGET'

datExpr=readRDS('../NB/data/TARGET/RNAseq/fromGDC/EXP_HR.RDS')
print(max(datExpr))

X=dplyr::select(as.data.frame(t(datExpr)),names(risk_factor))
dim(X)

if(max(datExpr)>50){
  X=log(X+1,2)
}
print(max(X))

Surv.lst=readRDS('/local/yanzijun/CRU/NB/data/TARGET/RNAseq/fromGDC/clin/surv.info.RDS')
clin_info=Surv.lst$OS
clin_info <- clin_info[!duplicated(clin_info$sampleID),]
rownames(clin_info)=clin_info$sampleID;clin_info$sampleID=NULL

ol.samples <- intersect(rownames(X),rownames(clin_info))
print(length(ol.samples))

clin_info <- clin_info[rownames(clin_info)%in% ol.samples,]
X= as.data.frame(t(dplyr::select(as.data.frame(t(X)),rownames(clin_info))))

print(all(rownames(X)==rownames(clin_info)))
print(all(colnames(X)==names(risk_factor)))

RS_Valid=c()
for(i in 1:nrow(X)){
  RS_Valid[i]=(risk_factor[1]*X[i,1])+(risk_factor[2]*X[i,2])+(risk_factor[3]*X[i,3])
  +(risk_factor[4]*X[i,4])+(risk_factor[5]*X[i,5])+(risk_factor[6]*X[i,6])
}
names(RS_Valid)=rownames(X)


## Fig3A: plot RS
rs_value <- sort(RS_Valid)
left <- rs_value[1:floor(length(rs_value)/2)]
right <- rs_value[floor(length(rs_value)/2+1):length(rs_value)]

pdf(paste('FIG_scMK/RS_',validSet,'_plot.pdf',sep=''),width=8,height=6)
plot(x=1:length(left),y=left,col=mycol[2],pch=19,xlim=c(0,length(RS_Valid)),ylim=c(min(rs_value),max(rs_value)),
     xlab="",ylab="",cex=1.2,cex.lab=1.3,cex.axis=1.2,xaxt='n')
par(new=TRUE)
plot(x=(length(left)+1):length(RS_Valid),y=right,col=mycol[1],pch=19,xlim=c(0,length(RS_Valid)),ylim=c(min(rs_value),max(rs_value)),
     xaxt='n',yaxt='n',xlab="",ylab="Risk Score",cex=1.2,cex.lab=1.3,cex.axis=1.2)
abline(v=length(rs_value)/2+0.5,col = "gray", lwd = 3.5, lty = 2,xaxt='n')
title(main=paste(validSet,' (validation set)',sep=''))
dev.off()


##Fig3A: plot follow time
print(all(names(RS_Valid)==rownames(clin_info)))
data=cbind(rs=RS_Valid,status=clin_info)
data=data[order(data[,1]),]
data=as.data.frame(cbind(data,rank=1:nrow(data)))
colnames(data) <- c('rs','status','time','rank')

alive=data[which(data$status==0),]
dead=data[which(data$status==1),]

pdf(paste('FIG_scMK/FollowTime_',validSet,'_plot.pdf',sep=''),width=8,height=6)
plot(x=dead$rank,y=dead$time,col=mycol[1],pch=17,xlim=c(0,nrow(data)),ylim=c(min(data$time),max(data$time)),
     xlab="",ylab="",cex=1.2,cex.lab=1.3,cex.axis=1.2)
par(new=TRUE)
plot(x=alive$rank,y=alive$time,col=mycol[2],pch=17,xlim=c(0,nrow(data)),ylim=c(min(data$time),max(data$time)),
     xaxt='n',yaxt='n',xlab="Patients cohort",ylab="Follow-up time(days)",cex=1.2,cex.lab=1.3,cex.axis=1.2)
abline(v=length(rs_value)/2+0.5,col = "gray", lwd = 3.5, lty = 2)
dev.off()


##Fig3B: survival analysis
cut_off=quantile(RS_Valid,0.5)
sur_df <- data
sur_df$time=round(sur_df$time/365,3)
sur_df$surtype <- 'UHR'
sur_df$surtype[sur_df$rs < cut_off]='CHR'
sur_df$surtype <- factor(sur_df$surtype,levels = c('UHR','CHR'))
table(sur_df$surtype)


fit <- survfit(Surv(time =time, event =status) ~ surtype,data=sur_df)
p <- ggsurvplot(fit, pval = TRUE,palette ='Set1',risk.table = TRUE,
                ggtheme = theme_survminer(base_size = 10),
                title=paste('Survival curve of ',validSet,' (validation set)',sep='') )## month
p

pdf(paste('FIG_scMK/survCurve_',validSet,'.pdf',sep=''),width = 5,height = 5.5)
print(p)
dev.off()

print(surv_pvalue(fit)$pval)

valid.res <- list(RS=RS_Valid,X.model=X,clin_data=data)
saveRDS(valid.res,paste('res_scMK/Lasso/',validSet,'.res.RDS',sep = ''))


######
### 4. Fig2E/3C: tROC analysis
######
rm(list=ls())
.libPaths(c("/local/yzj/R/x86_64-pc-linux-gnu-library/4.0" ,
            "/local/yanzijun/R/x86_64-pc-linux-gnu-library/4.0" ))
library(ggplot2)
library(dplyr)
library("survminer")
library(survivalROC)
library(rmda)
library(RColorBrewer);
mycol = colorRampPalette(brewer.pal(9, "Set1"))(9)

setwd('/local/yanzijun/CRU/NB_FN/')

datSet='GSE85047' #TARGET
if(datSet=='GSE85047'){
  model.res <- readRDS('res_scMK/Lasso/model.res.RDS')
}else if(datSet=='TARGET'){
  model.res <- readRDS('res_scMK/Lasso/TARGET.res.RDS')
}

clin_info <- model.res$clin_data
clin_info$sampleID=rownames(clin_info);

df <- clin_info
head(df)


ROC_1 <- survivalROC(Stime = df$time,status = df$status,marker = df$rs,predict.time = 365*1,method = 'KM')
ROC_3 <- survivalROC(Stime = df$time,status = df$status,marker = df$rs,predict.time = 365*3,method = 'KM')
ROC_5 <- survivalROC(Stime = df$time,status = df$status,marker = df$rs,predict.time = 365*5,method = 'KM')

pdf(paste('FIG_scMK/ROC_',datSet,'.pdf',sep=''),width = 5,height = 5)
plot(ROC_1$FP,ROC_1$TP,type='l',col='#BC3C29FF',xlim=c(0,1),ylim=c(0,1),
     xlab='False positive rate',
     ylab='True positive rate', main=paste('Time-dependent ROC curve\n',datSet,sep=''),lwd=2.5)
abline(0,1,col='gray',lwd=2,lty=2)
lines(ROC_3$FP,ROC_3$TP,type='l',col='#0072B5FF',xlim=c(0,1),ylim=c(0,1),lwd=2.5,lty=1)
lines(ROC_5$FP,ROC_5$TP,type='l',col='#20854EFF',xlim=c(0,1),ylim=c(0,1),lwd=2.5,lty=1)
legend(0.5,0.2,c(paste('AUC at 1 year:',round(ROC_1$AUC,3)),
                 paste('AUC at 3 year :',round(ROC_3$AUC,3)),
                 paste('AUC at 5 year :',round(ROC_5$AUC,3))),
       x.intersp=1, y.intersp=0.8,lty=1, lwd=2, col=mycol[1:3],bty='n',seg.len = 1,cex=0.8)
dev.off()
