################################################################
### step2: get METcor genes based on bulk expression and methylation data from 56 HR patients (sFig2)
################################################################
############
## 1. processing bulk RNA-seq data
############
rm(list=ls())
library(GEOquery)
GSE_ID='GSE73517'
setwd(paste('/local/yanzijun/CRU/NB/data/',GSE_ID,sep=''))

# load series and platform data from GEO
gset <- getGEO(GSE_ID, GSEMatrix =TRUE, getGPL=FALSE)
print(names(gset))

if(length(gset)==1){
  GPL_ID=gset[[paste(GSE_ID,'_series_matrix.txt.gz',sep='')]]@annotation
}

if (length(gset) > 1) idx <- grep(GPL_ID, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
saveRDS(gset,'gset.RDS')

#### EXP mtx
EXP <- exprs(gset)

# log2 transform
qx <- as.numeric(quantile(EXP, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { EXP[which(EXP <= 0)] <- NaN
EXP <- log2(EXP) }
dim(EXP)

GPL <- getGEO(GPL_ID, destdir='/local/yzj/publicData/GPL')
print(names(Table(GPL)))
tmp=Table(GPL)
ids=tmp[,c('ID','GeneSymbol')]

colnames(ids) <- c('probe_id','symbol')
head(ids)
ids=na.omit(ids)
dim(ids)

EXP <- tibble::rownames_to_column(as.data.frame(EXP),'ID_REF')
median=apply(EXP[,2:ncol(EXP)],1,median)
df_median <- data.frame(probe_id=EXP$ID_REF,median=median)

ids <- merge(ids,df_median,by='probe_id')
ids=ids[order(ids$symbol,ids$median,decreasing = T),]
ids=ids[!duplicated(ids$symbol),1:2]
colnames(ids)[1] <- 'ID_REF'

EXP=merge(ids,EXP,by='ID_REF')
EXP$ID_REF <- NULL

dim(EXP)
dim(na.omit(EXP))
EXP[is.na(EXP)] <- 0
saveRDS(EXP,'EXP.RDS')

#### clin info
meta <- gset@phenoData@data
colnames(meta)
table(apply(as.matrix(tmp),1,function(x) unlist(strsplit(x,','))[3]))
write.csv(meta,'clin.info.csv')


### get HR/LR Expression
GSM.HR <- na.omit(meta$geo_accession[meta$current.risk.category.ch1=='high-risk'])
GSM.LR <- na.omit(meta$geo_accession[meta$current.risk.category.ch1=='low-risk'])
print(length(GSM.HR))
print(length(GSM.LR))

rownames(EXP) <- EXP$symbol;EXP$symbol=NULL
EXP <- dplyr::select(EXP,c(GSM.LR,GSM.HR))
EXP_HR <- dplyr::select(EXP,GSM.HR)
EXP_LR <- dplyr::select(EXP,GSM.LR)
saveRDS(EXP,'EXP_new_risk.RDS')
saveRDS(EXP_HR,'EXP_HR.RDS')
saveRDS(EXP_LR,'EXP_LR.RDS')

EXP_HR <- readRDS('EXP_HR.RDS')
write.csv(EXP_HR,'GSE73517_HRNB_expression.csv')

############
## 2. processing bulk methylation data
############
rm(list=ls())
library(GEOquery)
GSE_ID="GSE73515"
setwd(paste('/local/yanzijun/CRU/NB/data/',GSE_ID,sep=''))

# load series and platform data from GEO
gset <- getGEO(GSE_ID, GSEMatrix =TRUE, getGPL=FALSE)
print(names(gset))

if(length(gset)==1){
  GPL_ID=gset[[paste(GSE_ID,'_series_matrix.txt.gz',sep='')]]@annotation
}
print(GPL_ID)

if (length(gset) > 1) idx <- grep(GPL_ID, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
saveRDS(gset,'gset.RDS')

####  mtx
gset <- readRDS('gset.RDS')
EXP <- exprs(gset)

GPL <- getGEO('GPL13534', destdir='/local/yzj/publicData/GPL')
print(names(Table(GPL)))
ids=Table(GPL)
saveRDS(ids,'map.info.RDS')
saveRDS(EXP,'EXP.RDS')

#### clin info
meta <- gset@phenoData@data
colnames(meta)
write.csv(meta,'clin.info.csv')

GSM.HR <-meta$geo_accession[meta$current.risk.category.ch1=='high-risk']
GSM.LR <- meta$geo_accession[meta$current.risk.category.ch1=='low-risk']
length(GSM.HR)

### filter probes
## delete probed with NA values 
dim(na.omit(EXP))
## delete probes in chrX/Y
probes <- ids$ID[!ids$CHR %in% c('X','Y')]
sub.EXP <- EXP[rownames(EXP) %in% probes,]
HR.EXP <- dplyr::select(as.data.frame(sub.EXP),GSM.HR)
dim(HR.EXP)
saveRDS(HR.EXP,'EXP_HR.RDS')
write.csv(MET,'GSE73515_HRNB_methylation.csv')

### change probeID to gene symbol
rm(list=ls())
library(dplyr)
library(tidyr)
RNA <- readRDS('/local/yanzijun/CRU/NB/data/GSE73517/EXP_HR.RDS')
dim(RNA)
MET <- readRDS('/local/yanzijun/CRU/NB/data/GSE73515/EXP_HR.RDS')
dim(MET)

RNA.clin <- read.csv('/local/yanzijun/CRU/NB/data/GSE73517/clin.info.csv')
MET.clin <- read.csv('/local/yanzijun/CRU/NB/data/GSE73515/clin.info.csv')
print(all(RNA.clin$current.risk.category.ch1==MET.clin$current.risk.category.ch1))
RNA.clin$sampleID <- paste('SampleID',1:nrow(RNA.clin),sep='')
MET.clin$sampleID <- paste('SampleID',1:nrow(MET.clin),sep='')

RNA.map <- RNA.clin$sampleID
names(RNA.map) <- RNA.clin$geo_accession
MET.map <- MET.clin$sampleID
names(MET.map) <- MET.clin$geo_accession

colnames(RNA) <- RNA.map[colnames(RNA)]
colnames(MET) <- MET.map[colnames(MET)]

RNA <- dplyr::select(RNA,colnames(MET))
print(all(colnames(RNA)==colnames(MET)))

ids.map <- readRDS('/local/yanzijun/CRU/NB/data/GSE73515/map.info.RDS')
ids <- ids.map[,c('ID','UCSC_RefGene_Name')]
MET_new <- MET
MET_new <- tibble::rownames_to_column(MET_new,'ID')
MET_new <- merge(ids,MET_new)
MET_new <- MET_new[which(MET_new$UCSC_RefGene_Name != ""),]
dim(MET_new)
MET_new <- MET_new %>%
  dplyr::arrange(UCSC_RefGene_Name)


MET_new <- aggregate(x = MET_new[,3:ncol(MET_new)],
                             by = list(MET_new$UCSC_RefGene_Name),
                             FUN = mean)

MET_new2 <- c()
for(i in 1:nrow(MET_new)){
  genes <- unique(unlist(strsplit(MET_new$Group.1[i],';')))
  if(length(genes)==1){
    MET_new2 <- rbind(MET_new2,cbind(data.frame(gene=genes),MET_new[i,2:ncol(MET_new)]))
  }else if(length(genes)>1){
    for(j in 1:length(genes)){
      MET_new2 <- rbind(MET_new2,c(data.frame(gene=genes[j]),MET_new[i,2:ncol(MET_new)]))
    }
  }
}
saveRDS(MET_new2,'/local/yanzijun/CRU/NB/data/GSE73515/EXP_tmp.RDS')

MET_new3 <- aggregate(x = MET_new2[,2:ncol(MET_new2)],
                     by = list(MET_new2$gene),
                     FUN = mean)
dim(MET_new3)
colnames(MET_new3)[1] <- 'gene'
saveRDS(MET_new3,'/local/yanzijun/CRU/NB/data/GSE73515/EXP_HR_gene.RDS')


############ 
## 3. calculate Pearson correlation
############
rm(list=ls())
library(dplyr)
library(tidyr)

## input data
RNA <- readRDS('/local/yanzijun/CRU/NB/data/GSE73517/EXP_HR.RDS')
dim(RNA)
MET <- readRDS('/local/yanzijun/CRU/NB/data/GSE73515/EXP_HR.RDS')
dim(MET)

RNA.clin <- read.csv('/local/yanzijun/CRU/NB/data/GSE73517/clin.info.csv')
MET.clin <- read.csv('/local/yanzijun/CRU/NB/data/GSE73515/clin.info.csv')
print(all(RNA.clin$current.risk.category.ch1==MET.clin$current.risk.category.ch1))
RNA.clin$sampleID <- paste('SampleID',1:nrow(RNA.clin),sep='')
MET.clin$sampleID <- paste('SampleID',1:nrow(MET.clin),sep='')

RNA.map <- RNA.clin$sampleID
names(RNA.map) <- RNA.clin$geo_accession
MET.map <- MET.clin$sampleID
names(MET.map) <- MET.clin$geo_accession
map.lst <- list(RNA.map=RNA.map,MET.map=MET.map)

colnames(RNA) <- RNA.map[colnames(RNA)]
colnames(MET) <- MET.map[colnames(MET)]

RNA <- dplyr::select(RNA,colnames(MET))
print(all(colnames(RNA)==colnames(MET)))

MET_new3 <- readRDS('/local/yanzijun/CRU/NB/data/GSE73515/EXP_HR_gene.RDS')

MET_MK <- MET_new3
rownames(MET_MK) <- MET_MK$gene
MET_MK$gene <- NULL

RNA_MK <- RNA
dim(RNA_MK)
dim(MET_MK)

commonG <- intersect(rownames(MET_MK),rownames(RNA_MK))
length(commonG)

MET_MK <- MET_MK[match(commonG,rownames(MET_MK)),]
RNA_MK <- RNA_MK[match(commonG,rownames(RNA_MK)),]

print(all(rownames(MET_MK)==rownames(RNA_MK)))
print(all(colnames(MET_MK)==colnames(RNA_MK)))

print('Run correlation!')

pvalues <- c()
cors <- c()
for(i in 1:length(commonG)){
  print(i)
  res <- cor.test(as.numeric(MET_MK[i,]),as.numeric(RNA_MK[i,]))
  pvalues[i] <- as.numeric(res$p.value)
  cors[i] <- as.numeric(res$estimate)
}
res <- data.frame(gene=commonG,cor=cors,pvalue=pvalues)
saveRDS(res,'/local/yanzijun/CRU/NB_FN/res/GeneSet/Cor/Pearson_all.RDS')

## get METcor genes
cor.df <- readRDS('/local/yanzijun/CRU/NB_FN/res/GeneSet/Cor/Pearson_all.RDS')
pcutoff=0.001
sigcor.df <- cor.df[cor.df$p < pcutoff & abs(cor.df$cor) >=0.4,]
print(nrow(sigcor.df))
saveRDS(sigcor.df,'/local/yanzijun/CRU/NB_FN/res/GeneSet/Cor/METcor.RDS')


############ 
## 4. plot
############
rm(list=ls())
library(ggplot2)
library(RColorBrewer)
mycol = colorRampPalette(brewer.pal(9, "Set1"))(9)

## plot METcor genes distribution
cor.df <- readRDS('/local/yanzijun/CRU/NB_FN/res/GeneSet/Cor/Pearson_all.RDS')
cor.df$log10pvalue <- -log(cor.df$pvalue,10)
cor.df$Type=rep('nonSig',nrow(cor.df))

pcutoff=0.001
cor.df$Type[cor.df$pvalue < pcutoff & cor.df$cor >=0.4] <- 'pos.METcor'
cor.df$Type[cor.df$pvalue < pcutoff & cor.df$cor <=-0.4] <- 'neg.METcor'

## sFig2A: correlation
p.cor <- ggplot(data=cor.df, aes(x=cor,y=log10pvalue,color=Type)) + 
  geom_point(alpha=0.7) +
  theme_classic(base_size = 9)+
  scale_color_manual(values = c(mycol[1],"grey",mycol[2]))+
  geom_hline(yintercept = -log(pcutoff,10),lty=3,col="black",lwd=0.5)+
  geom_vline(xintercept = 0.4,lty=3,col="black",lwd=0.5)+
  geom_vline(xintercept = -0.4,lty=3,col="black",lwd=0.5)+
  labs(x='Pearson correlation',y = expression(paste("-log"[10], "(", italic("P"), "-value)")))
p.cor
ggsave('/local/yanzijun/CRU/NB_FN/FIG_scMK/METcor_CorPlot.pdf',p.cor,width = 4.5,height = 3.5)

## sFig2B: gene distribution
sigcor.df <- cor.df[cor.df$Type !='nonSig',]
corMK.T <- data.frame(gene=sigcor.df$gene)

gtf <- read.table('/mdshare/node9/yzj/publicData/annotation/hg19/gencode_v19_gene_pos.txt',
                  header = F,sep='\t',stringsAsFactors = F)
colnames(gtf) <- c("gene","chr","start","end")
dim(gtf)

chr.num <- as.data.frame(table(gtf$chr))
colnames(chr.num) <- c('chr','num')

corMK.T.gtf <- merge(corMK.T,gtf,by='gene',all.x=TRUE)
dim(corMK.T.gtf)

corMK.T.num <- as.data.frame(table(corMK.T.gtf$chr))
colnames(corMK.T.num) <- c('chr','corMK.T.num')

num.df <- merge(corMK.T.num,chr.num)
dim(num.df)
num.df$freq <- round(num.df$corMK.T.num/num.df$num,3)
num.df$chr <- factor(num.df$chr,levels = rev(paste('chr',c(1:22,'X','Y'),sep='')))
num.p <- ggplot() +  geom_bar(data = num.df, 
                              aes(x = chr, y = freq),
                              stat = "identity",
                              width = 0.8, fill = mycol[2], colour="white",
                              position = position_dodge(width = 1))+
  theme_classic()+
  theme(axis.text=element_text(size=9,colour = 'black'),
        # axis.text.x = element_text(angle = 45, hjust = 0.5,vjust = 0.5)
  )+
  labs(x='',y = "Frequency")+coord_flip()
num.p
ggsave('/local/yanzijun/CRU/NB_FN/FIG_scMK/METcor_chrnum.pdf',num.p,height = 6.5, width = 3.5)


## sFig2C: type distribution
ids.map <- readRDS('/local/yanzijun/CRU/NB/data/GSE73515/map.info.RDS')
sub.ids.map <- ids.map[,c('UCSC_RefGene_Name','Relation_to_UCSC_CpG_Island')]

corMK.T.map <- sub.ids.map[sub.ids.map$UCSC_RefGene_Name %in% corMK.T$gene,]
corMK.T.map$Relation_to_UCSC_CpG_Island[corMK.T.map$Relation_to_UCSC_CpG_Island==""] <- 'Other'
type.df <- as.data.frame(table(corMK.T.map$Relation_to_UCSC_CpG_Island))
colnames(type.df) <- c('type','num')
type.df$freq <- type.df$num/sum(type.df$num)

freq.p <- ggplot() +  geom_bar(data = type.df, 
                               aes(x = type, y = freq),
                               stat = "identity",
                               width = 0.8, fill = mycol[2], colour="white",
                               position = position_dodge(width = 1))+
  theme_classic()+
  theme(axis.text=element_text(size=9,colour = 'black'),
        #axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)
  )+
  labs(x='',y = "Frequency")+coord_flip()
freq.p
ggsave('/local/yanzijun/CRU/NB_FN/FIG_scMK/METcor_Islandfreq.pdf',freq.p,height = 2.5, width = 4)


## SFig 3D: six-signature genes: correlation between expression and methylation
rm(list=ls())
geneLst <- c('RAMP1','CDT1','MEG3','NPW','MAPT','C1QTNF4')
res <- readRDS('/local/yanzijun/CRU/NB_FN/res/GeneSet/Cor/Pearson.RDS')

map.lst <- readRDS('/local/yanzijun/CRU/NB/data/GSE73515/sample2ID.RDS')
RNA.map <- map.lst$RNA.map
RNA <- readRDS('/local/yanzijun/CRU/NB/data/GSE73517/EXP_HR.RDS')
colnames(RNA) <- RNA.map[colnames(RNA)]
MET <- readRDS('/local/yanzijun/CRU/NB/data/GSE73515/EXP_HR_gene.RDS')
rownames(MET) <- MET$gene;MET$gene=NULL


subRNA <- RNA[match(geneLst,rownames(RNA)),]
subMET <- MET[match(geneLst,rownames(MET)),]
all(rownames(subRNA)==rownames(subMET))
all(colnames(subRNA)==colnames(subMET))

plots <- list()
for(i in 1:length(geneLst)){
  g=geneLst[i]
  df <- data.frame(Expression=as.numeric(subRNA[rownames(subRNA)==g,]),
                   Methylation=as.numeric(subMET[rownames(subMET)==g,]))
  p <- ggplot(data=df, aes(x=Expression, y=Methylation))+geom_point(size=1)+stat_smooth(method="lm")+
    theme_bw()+
    theme(axis.text =element_text(size=9,colour = 'black'),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(x=paste(g,'expression',sep=' '),y=paste(g,'methylation',sep=' '))
  plots[[i]] <- p
}
library(ggpubr)
cor.plot <- ggarrange(plotlist = plots,nrow = 2,ncol = 3)
ggsave('/local/yanzijun/CRU/NB_FN/FIG_scMK/corPlot.pdf',width = 20,height = 12,units = 'cm')