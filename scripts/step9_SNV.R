################################################################
### step9: mutaion landscape between UHR and CHR group (sFig6)
################################################################
rm(list=ls())
.libPaths(c("/local/yzj/R/x86_64-pc-linux-gnu-library/4.0",
            "/local/yanzijun/R/x86_64-pc-linux-gnu-library/4.0"))
library(RColorBrewer);
mycol = colorRampPalette(brewer.pal(9, "Set1"))(9)
library(maftools)
library(ggpubr)
library(ggplot2)

setwd('/local/yanzijun/CRU/NB/data/TARGET/SNV/')
clinSNV <- read.table('/local/yanzijun/CRU/NB/data/TARGET/SNV/clin/gdc_sample_sheet.2022-06-09.tsv',header = T,sep='\t')
clinSNV$Tumor_Sample_Barcode <- apply(as.matrix(clinSNV$Case.ID),1,function(x) unique(unlist(strsplit(x,', '))))

model.res <- readRDS('/local/yanzijun/CRU/NB_FN/res_scMK/Lasso/TARGET.res.RDS')
RS <- model.res$clin_data
RS<- tibble::rownames_to_column(RS,'sampleID')
RS$type='high risk'
RS$type[RS$rs < quantile(RS$rs,0.5)]='low risk'
head(RS)

ol <- intersect(RS$sampleID,clinSNV$Tumor_Sample_Barcode)
length(ol)

clin=RS[RS$sampleID %in% ol,]
sub.clinSNV <- clinSNV[clinSNV$Tumor_Sample_Barcode %in% ol,]

########
## 1. sFig 6A: plot top mutaion
########
maf_path <- paste('maf/',paste(sub.clinSNV$File.ID,sub.clinSNV$File.Name,sep='/'),sep='')
maf=maftools:::merge_mafs(maf=maf_path, verbose = TRUE)

## Tumor_Sample_Barcode
sample.sum <- getSampleSummary(maf)
tmp <- apply(as.matrix(sample.sum$Tumor_Sample_Barcode),1,
             function(x) paste0(unlist(strsplit(x,'-'))[1:3],collapse ='-'))
ID <- sample.sum$Tumor_Sample_Barcode
names(ID)=tmp
clin$Tumor_Sample_Barcode <- ID[clin$sampleID]

maf=maftools:::merge_mafs(maf=maf_path, verbose = TRUE,clinicalData=clin)
getClinicalData(maf)
getGeneSummary(maf)
getFields(maf)

vc_cols = mycol[1:8]
#vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)

Typecolors = mycol[1:2]
names(Typecolors) = c("high risk","low risk")
Typecolors = list(type = Typecolors)

pdf('/local/yanzijun/CRU/NB_FN/FIG_scMK/oncoplot_all.pdf',width = 10,height = 6)
oncoplot(
  maf = maf,
  #colors = vc_cols,
  clinicalFeatures = 'type',
  sortByAnnotation = TRUE,
  annotationColor = Typecolors
)
dev.off()

########
## 2. sFig6B-C: association between ALK mutaion/WT and survival
########
clin.HR <- subset(clin, type=="high risk")$Tumor_Sample_Barcode
clin.LR <- subset(clin, type=="low risk")$Tumor_Sample_Barcode

maf.HR <- subsetMaf(maf=maf, tsb=clin.HR, isTCGA=FALSE)
maf.LR <- subsetMaf(maf=maf, tsb=clin.LR, isTCGA=FALSE)

gene='ALK'
pdf(paste('/local/yanzijun/CRU/NB_FN/FIG_scMK/oncoplot_',gene,'.pdf',sep=''),width = 4,height = 4)
mafSurvival(maf=maf, genes=gene, time="time", Status="status", isTCGA=FALSE)
mafSurvival(maf=maf.HR, genes=gene, time="time", Status="status", isTCGA=FALSE)
mafSurvival(maf=maf.LR, genes=gene, time="time", Status="status", isTCGA=FALSE)
dev.off()
