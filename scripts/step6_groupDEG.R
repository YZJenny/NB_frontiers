################################################################
### step6: DEG/GSEA analysis between UHR and CHR group (Fig5/sFig5)
################################################################
rm(list=ls())
library(readxl)
library(limma)
library(dplyr)
setwd('/local/yanzijun/CRU/NB_FN/')
SEQC.exp <- readRDS('../NB/data/GSE62564/EXP_HR.RDS')
TARGET.exp <- readRDS('../NB/data/TARGET/RNAseq/fromGDC/EXP_HR.RDS')

GSE_ID='GSE62564'
if(GSE_ID=='GSE62564'){
  EXP=SEQC.exp
  model.res <- readRDS('res_scMK/Lasso/model.res.RDS')
  RS <- model.res$clin_data
}else if(GSE_ID=='TARGET'){
  EXP=TARGET.exp
  model.res <- readRDS('res_scMK/Lasso/TARGET.res.RDS')
  RS <- model.res$clin_data
}


GSM.HR=rownames(RS)[RS$rs >= quantile(RS$rs,0.5)]
GSM.LR=rownames(RS)[RS$rs < quantile(RS$rs,0.5)]


####################    
##### 1. DEG analysis: Limma
#################### 
print(length(GSM.HR))
print(length(GSM.LR))

## 执行limma的DEG
label <- c(rep('LR',length(GSM.LR)),rep('HR',length(GSM.HR)))

EXP <- dplyr::select(EXP,c(GSM.LR,GSM.HR))
EXP_HR <- dplyr::select(EXP,GSM.HR)
EXP_LR <- dplyr::select(EXP,GSM.LR)

Group=data.frame(sampleID=colnames(EXP),Group=label)

group_list <- factor(label,levels=c('LR','HR'))
design <- model.matrix(~0+group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(EXP)

contrast.matrix <- makeContrasts(HR-LR, levels=design)
fit <- lmFit(EXP, design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit,trend=T,robust=T) #FPKM

output <- topTable(fit, coef=1,n=Inf,adjust.method='fdr')
output <- tibble::rownames_to_column(output,'gene') 

system(paste('mkdir -p /local/yanzijun/CRU/NB_FN/res_scMK/Enrich/groupDEG/',GSE_ID,sep=''))
write.table(output,paste('/local/yanzijun/CRU/NB_FN/res_scMK/Enrich/groupDEG/',GSE_ID,'/Limma_output.txt',sep=''),sep='\t',row.names = F,quote = F)

## DEGlst
for(FC in c(1,1.5,2)){
  print(FC)
  up.Genes <- output$gene[output$logFC > log(FC,2) & output$adj.P.Val < 0.05]
  dn.Genes <- output$gene[output$logFC < -log(FC,2) & output$adj.P.Val < 0.05]
  print(length(up.Genes))
  print(length(dn.Genes))
  GeneLst <- list(up.Genes=up.Genes,dn.Genes=dn.Genes)
  saveRDS(GeneLst,paste('/local/yanzijun/CRU/NB_FN/res_scMK/Enrich/groupDEG/',GSE_ID,'/GeneLst_FC',FC,'.RDS',sep=''))
}

################
#### 2. sFig5: plot heatmap
################
rm(list=ls())
library(readxl)
library(limma)
library(dplyr)
library(RColorBrewer)
mycol = colorRampPalette(brewer.pal(9, "Set1"))(9)

setwd('/local/yanzijun/CRU/NB_FN/')
EXP <- readRDS('../NB/data/GSE62564/EXP_HR.RDS')
model.res <- readRDS('res_scMK/Lasso/model.res.RDS')
RS <- model.res$clin_data
GSM.HR=rownames(RS)[RS$rs >= quantile(RS$rs,0.5)]
GSM.LR=rownames(RS)[RS$rs < quantile(RS$rs,0.5)]

output <- read.table('/local/yanzijun/CRU/NB_FN/res_scMK/Enrich/groupDEG/GSE62564/Limma_output.txt',
                     sep='\t',header = TRUE)
output <- output[output$logFC !=0 & output$adj.P.Val < 0.05,]
write.csv(output, '/local/yanzijun/CRU/NB_FN/res_scMK/Enrich/groupDEG/GSE62564/sigLimma_output_FC1.csv',
          row.names = F)
dim(output)
output <- output[order(output$logFC,decreasing = T),]
topgene <- c(output$gene[1:20],tail(output$gene,20))

heatmap.df <- EXP[match(topgene,rownames(EXP)),]
heatmap.df <- dplyr::select(heatmap.df,c(GSM.HR,GSM.LR))

annotation_col = data.frame(Type = factor(c(rep("CHR", length(GSM.LR)),
                                            rep('UHR',length(GSM.HR)))))
rownames(annotation_col) = colnames(heatmap.df)

heatmap.df[heatmap.df>8]=8
heatmap.df[heatmap.df<-8]=-8
p.heatmap <- pheatmap(heatmap.df,cluster_rows =F,cluster_cols = F,annotation_col = annotation_col,
                      legend = TRUE, annotation_legend = TRUE, show_colnames =  F,
                      color= colorRampPalette(c("navy", "white", "firebrick3"))(50),border_color=NA,
                      annotation_colors =list(Type=c(CHR=mycol[1],UHR=mycol[2])),
                      fontsize=10,fontsize_row=10)
ggsave(filename = 'FIG_scMK/DEG_heatmap.pdf',p.heatmap,width = 14,height = 14,units = 'cm')


####################
##### 3. GSEA
####################
rm(list=ls())
#.libPaths(c("/local/yzj/R/x86_64-pc-linux-gnu-library/4.0",'/local/yanzijun/tools/R-4.0.0/library'))
library(Biobase)
library(genefilter)
library(RColorBrewer)
library(AnnotationHub)
library(org.Hs.eg.db) 
library(clusterProfiler)
library(DOSE)
library(dplyr)
library(tidyverse)
library(res_scMKhape2)
library(enrichplot)
library(RColorBrewer);
library(scater)
mycol = colorRampPalette(brewer.pal(9, "Set1"))(10)

setwd('/local/yanzijun/CRU/NB_FN/res_scMK/Enrich/groupDEG/')

GSE_ID='GSE62564'
output <- read.table(paste('/local/yanzijun/CRU/NB_FN/res_scMK/Enrich/groupDEG/',GSE_ID,'/Limma_output.txt',sep=''),sep='\t',
                     header = T)
output[1:2,]
df=output[,c('gene','logFC')]

## get ENTREZID
symbol2id=mapIds(org.Hs.eg.db,df$gene,"ENTREZID",'SYMBOL')
symbol2id=symbol2id[which(symbol2id!='')] 
symbol2id <- data.frame(gene=names(symbol2id),id=symbol2id)

df <- merge(df,symbol2id)

## sort FC
sortdf<-df[order(df$logFC, decreasing = T),]
head(sortdf)
gene.expr = sortdf$logFC
names(gene.expr) <- sortdf$id
head(gene.expr)


### GO
go <- gseGO(gene.expr, OrgDb = org.Hs.eg.db,ont = 'BP',keyType = "ENTREZID")
dim(go)

sortgo <- go[order(go$enrichmentScore, decreasing = T),]
head(sortgo$Description)
dim(sortgo)
write.table(sortgo,paste(GSE_ID,'/gseaGO.txt',sep=''),sep = "\t",quote = F,col.names = T,row.names = F)

GO.p <- gseaplot2(go,
                  row.names(sortgo)[1:10],
                  title = "GO BP (Ultra group)",
                  base_size = 7,
                  color = mycol[1:length(terms)],
                  pvalue_table = FALSE,
                  ES_geom="line")
ggsave(paste(GSE_ID,'/gseaGO.pdf',sep=''),GO.p,width = 6,height = 5)


#### Fig5A-B: KEGG
kk <- gseKEGG(gene.expr, organism = "hsa")
dim(kk)

sortkk <- kk[order(kk$enrichmentScore, decreasing = T),]
head(sortkk$Description,10)
dim(sortkk)
write.table(sortkk,paste(GSE_ID,'/gseaKEGG.txt',sep=''),sep = "\t",quote = F,col.names = T,row.names = F)

library(enrichplot)
N=5
type='DN'
if(type=='UP'){
  terms <- sortkk$ID[1:N]
  tt="Enriched in ultra group"
}else{
  terms <- tail(sortkk$ID,N)
  tt="Enriched in common group"
}

KEGG.p <- gseaplot2(kk,
                    terms,
                    title = tt,
                    base_size = 8,
                    color = mycol[1:length(terms)],
                    pvalue_table = FALSE,
                    ES_geom="line")
KEGG.p
ggsave(paste(GSE_ID,'/gseaKEGG_',type,'.pdf',sep='l'),KEGG.p,
       width = 10,height = 9,units = 'cm')



#### Fig5C-D: Hallmarks
geneset <- read.gmt('/mdshare/node9/yzj/publicData/GMT/h.all.v6.2.symbols.gmt') 
gene.expr2 = sortdf$logFC
names(gene.expr2) <- sortdf$gene
head(gene.expr2)
hm <- GSEA(gene.expr2, TERM2GENE=geneset, verbose=FALSE)
HM.df <- hm@result
HM.df <- HM.df[order(HM.df$NES,decreasing = T),]
dim(HM.df)

N=5
type='DN'
if(type=='UP'){
  terms <- HM.df$Description[1:N]
  tt="Enriched in Ultra group"
}else{
  terms <- tail(HM.df$Description,N)
  tt="Enriched in Common group"
}

HM.p <- gseaplot2(hm,
                    terms,
                    title = tt,
                    base_size = 10,
                    color = mycol[1:length(terms)],
                    pvalue_table = FALSE,
                    ES_geom="line")
HM.p
ggsave(paste(GSE_ID,'/gseaHM_',type,'.pdf',sep=''),HM.p,width = 6,height = 5)