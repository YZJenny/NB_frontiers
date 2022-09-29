################################################################
### step8: inflammatory activity analysis between UHR and CHR group (Fig7)
################################################################
rm(list=ls())
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
mycol = colorRampPalette(brewer.pal(9, "Set1"))(9)

IgG <- unique(c('POU2AF1','IGLC2','IGHM','IGH@','IGHG1','IGHM','IGHG3','LOC388078','IGHM','IGHG1','IGHG1',
                'IGKC','IGLJ3','LOC91316','LOC440871','IGKC','IGHD','IGLC2','IGLL1','IGKV1D‐13',
                'IGLJ3','IGHG1','IGLC2','IGLJ3','IGLC2','IGKC','IGL','IGHG1','LOC391427','LOC339562','IGKC'))
HCK <- unique(c('IFI30','LAPTM5','ITGB2','C1QB','SLCO2B1','CD163','TYROBP','FCER1G','SLC7A7','CCR1',
                'TFEC','HCK','NCF2','LAIR1','CD86','CD163','C1QA','MS4A4A','MNDA','AIF1',
                'LST1','AIF1','DOCK2','RNASE6','AIF1','MS4A6A'))
MHCII <- unique(c('HLA‐DPB1','PRG1','CTSS','HLA‐DMB','HLA‐DRB1','LCP2','PTPRC','HLA‐DRA','CD74',
                  'HLA‐DRA','HLA‐DPA1','HLA‐DQA1','HLA‐DRB1','HLA‐DMA'))
MHCI <- unique(c('HLA-A','HLA-B','HLA-C','HLA-F','HLA-G','HLA-J'))
LCK <- unique(c('CCL5','IL2RG','CD48','SELL','SCYA5','LCK','GZMA','IL7R','KLRK1','CD2','STAT4','TNFRSF7',
                'SLAMF1','CCR7','GZMK','CCR2','LTB','PRKCB1','CD3Z','SH2D1A','TRBC1','ITK','TRBC1',
                'CD3D','GIMAP5','PLAC8','GIMAP4','GIMAP5','PRG1','HCLS1','INPP5D','CD53','SLA',
                'PIK3CD','IRF8','GMFG','FGL2','IL10RA','CXorf9','CSF2RB','LCP2','CORO1A','HEM1',
                'SELPLG','EVI2B','PTPRC','RAC2','LPXN','ARHGAP15','SAMSN1','KIAA0053'))
STAT1 <- unique(c('STAT1','GBP1','TAP1','IRF1','CXCL9','PSMB9','CXCL10','STAT1','INDO',
                  'CXCL11','IFIH1'))
Interferon <- unique(c('MX1','IFI27','OAS1','IFIT1','G1P3','IFI44L','IFIT3','OAS2','G1P2','OAS1','RSAD2',
                       'IFI44','OAS3','FLJ20035'))

geneLst <- unique(c(IgG,HCK,MHCI,MHCII,LCK,STAT1,Interferon))

gene2Type <- rep(c('IgG','HCK','MHCI','MHCII',
                   'LCK','STAT1','Interferon'),
                 c(length(IgG),length(HCK),length(MHCI),length(MHCII),
                   length(LCK),length(STAT1),length(Interferon)))
names(gene2Type) <- geneLst


GSE_ID='GSE62564'
if(GSE_ID=='GSE62564'){
  EXP <- readRDS('/local/yanzijun/CRU/NB/data/GSE62564/EXP_HR.RDS')
  model.res <- readRDS('res_scMK/Lasso/model.res.RDS')
}else if(GSE_ID=='TARGET'){
  EXP <- readRDS('/local/yanzijun/CRU/NB/data/TARGET/RNAseq/fromGDC/EXP_HR.RDS')
  model.res <- readRDS('res_scMK/Lasso/TARGET.res.RDS')
}

sub.EXP <- EXP[match(geneLst,rownames(EXP)),] 
dim(sub.EXP)

RS <- model.res$clin_data

RS$type='UHR'
RS$type[RS$rs < quantile(RS$rs,0.5)]='CHR'
RS[1:2,]
CHRID <- rownames(RS)[RS$type=='CHR']
UHRID <- rownames(RS)[RS$type=='UHR']
ID2Type <- RS$type
names(ID2Type) <- rownames(RS)

sub.EXP <- dplyr::select(sub.EXP,c(CHRID,UHRID))
dim(sub.EXP)
sub.EXP <- na.omit(sub.EXP)
dim(sub.EXP)

### Fig7A: pheatmap
annotation_col = data.frame(Group = factor(ID2Type[colnames(sub.EXP)]))
rownames(annotation_col) = colnames(sub.EXP)

annotation_row = data.frame(IC = factor(gene2Type[rownames(sub.EXP)]))
rownames(annotation_row) = rownames(sub.EXP)


annotation_colors =list(Group=c(CHR=mycol[2],UHR=mycol[1]),
                        IC=c(IgG=mycol[3],HCK=mycol[4],MHCI=mycol[5],MHCII=mycol[6],
                             LCK=mycol[7],STAT1=mycol[8],Interferon=mycol[9]))


table(gene2Type[rownames(sub.EXP)])
p.heatmap <- pheatmap(sub.EXP,cluster_rows =F,cluster_cols = F,
                      annotation_col = annotation_col,annotation_row = annotation_row,
                      scale = 'row',
                      legend = TRUE, annotation_legend = TRUE, show_colnames =  F,
                      color= colorRampPalette(c("navy", "white", "firebrick3"))(50),border_color=NA,
                      annotation_colors = annotation_colors,
                      gaps_row = c(2,21,26,29,61,67), 
                      gaps_col = c(length(CHRID)),
                      fontsize=8,fontsize_row=8)

ggsave(paste('FIG_scMK/Infla_heatmap_',GSE_ID,'.pdf',sep=''),plot = p.heatmap,
       width = 13,height = 20,units = 'cm')


### Fig7B: ssGSEA
library(pheatmap)
library(GSVA)
library(GSVAdata)
library(clusterProfiler)
library(RColorBrewer);
library(scater)
mycol = colorRampPalette(brewer.pal(9, "Set1"))(9)

gene.set <- data.frame(term=rep('Inflammatory',length(geneLst)),gene=geneLst)

SEQC.exp <- readRDS('../NB/data/GSE62564/EXP_HR.RDS')
TARGET.exp <- readRDS('../NB/data/TARGET/RNAseq/fromGDC/EXP_HR.RDS')

datSet='training'
if(datSet=='training'){
  EXP=SEQC.exp
  model.res <- readRDS('res_scMK/Lasso/model.res.RDS')
  RS <- model.res$clin_data
}else if(datSet=='TARGET'){
  EXP=TARGET.exp
  model.res <- readRDS('res_scMK/Lasso/TARGET.res.RDS')
  RS <- model.res$clin_data
}

RS <- tibble::rownames_to_column(RS,'sampleID')
RS$type='UHR'
RS$type[RS$rs < quantile(RS$rs,0.5)]='CHR'
table(RS$type)


res.ssgsea <- gsva(as.matrix(EXP), gene.set, method = "ssgsea", kcdf = "Gaussian", min.sz = 10)

score.df <- data.frame(t(res.ssgsea))
score.df <- tibble::rownames_to_column(score.df,'sampleID')
head(score.df)
plot.df <- merge(RS[,c('sampleID','type')],score.df)
colnames(plot.df)[3] <- 'score'
head(plot.df)
plot.df$type=factor(plot.df$type,levels = c('CHR','UHR'))

p <- ggboxplot(plot.df, x = 'type', y = "score",
               color = 'type', add=NULL)+
  stat_compare_means(label = "p.format")+
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1 ,colour = 'black'),
        legend.position = 'none', )+
  scale_color_manual(values = mycol[2:1])+
  labs(title="Inflammtory activity", x="", y="ssGSEA score")+
  theme(plot.title = element_text(hjust = 0.5))
p

ggsave(paste('FIG_scMK/GSVA_Inflamm_',datSet,'.pdf',sep=''),p,width = 3,height = 4)


### Fig7C: correlation 
library(ggplot2)
library(ggsci)
library(pheatmap)
library(reshape2)
library(ggpubr)
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))(9)

df <- sub.EXP
df$group <- gene2Type[rownames(df)]

i=1
sub.df <- df[,c(i,ncol(df))]
colnames(sub.df)[1] <- 'SEQC'
plot.df <- aggregate(sub.df$SEQC,list(sub.df$group),mean)
colnames(plot.df) <- c('group',colnames(sub.EXP)[i])

for(i in 2:ncol(sub.EXP)){
  sub.df <- df[,c(i,ncol(df))]
  colnames(sub.df)[1] <- 'SEQC'
  mean.df <- aggregate(sub.df$SEQC,list(sub.df$group),mean)
  colnames(mean.df) <- c('group',colnames(sub.EXP)[i])
  plot.df <- merge(plot.df,mean.df,by='group')
}
rownames(plot.df) <- plot.df$group;plot.df$group=NULL

cor.df <- as.data.frame(t(plot.df))
print(all(rownames(cor.df)==rownames(RS)))
cor.df$riskscore=RS$rs

cormat <- round(cor(cor.df), 1) 
res1 <- cor.mtest(cor.df, conf.level = .95) 

col2 <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                               "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                               "#4393C3", "#2166AC", "#053061")))

pdf('/local/yanzijun/CRU/NB_FN/FIG_scMK/Inflam_corplot.pdf',width = 6,height = 8)
corrplot(cormat, p.mat = res1$p,insig = "label_sig",col = col2(200),
         sig.level = c(.001, .01, .05), pch.cex = .9, pch.col = "white") 
dev.off()