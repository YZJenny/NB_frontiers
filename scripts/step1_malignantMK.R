################################################################
### step1: get marker genes of malignant cells based on scRNA-seq data from 11 HR patients (Fig1)
################################################################
rm(list=ls())
library(gdata)
library(GEOquery)
library(data.table)
library(Seurat)
library(dplyr)
library(future)
library(future.apply)

##########
## 1. preprocessing NB scRNA-seq
##########
work_dir='/local/yanzijun/CRU/NB/data/CancerCell_2020/GSE137804_RAW/'
files <- list.files(work_dir)

meta.info <- read.csv('/local/yanzijun/CRU/NB/data/CancerCell_2020/GSE137804_tumor_dataset_annotation.csv')
rownames(meta.info)=meta.info$cellname

for(i in 1:length(files)){
  f=files[i]
  ID <- unlist(strsplit(f,'_'))[2]
  if(grepl('^T',ID,ignore.case = FALSE)){
    if(nchar(ID)==4){
      dat <- fread(paste(work_dir,f,sep=''))
      colnames(dat)[1] <- 'Symbol'
    }else if(nchar(ID)==3){
      dat <- fread(paste(work_dir,f,sep=''),header = T)
      dat$Gene_ID <- NULL
    }
    colnames(dat) <- c('Symbol',paste(ID,colnames(dat)[-1],sep='_'))
    dat <- as.data.frame(dat[!duplicated(dat$Symbol),])
    rownames(dat) <- dat$Symbol
    dat$Symbol <- NULL
    
    comcells <- intersect(meta.info$cellname,colnames(dat)[-1])
    print(length(comcells))
    meta = meta.info[meta.info$cellname %in% comcells,]
    print(dim(meta))
    dat <- select(dat,comcells)
    print(all(colnames(dat)==meta$cellname))
    
    pbmc=CreateSeuratObject(counts = dat, min.cells = 0, min.features = 0,
                            project = ID,meta.data = meta)
    saveRDS(pbmc,paste('/local/yanzijun/CRU/NB/res/RDS/pbmc_',ID,'.RDS',sep=''))
    print(f)
  }
}


pbmc_T10 <- readRDS('/local/yanzijun/CRU/NB/res/RDS/pbmc_T10.RDS')
pbmc_T19 <- readRDS('/local/yanzijun/CRU/NB/res/RDS/pbmc_T19.RDS')
pbmc_T27 <- readRDS('/local/yanzijun/CRU/NB/res/RDS/pbmc_T27.RDS')
pbmc_T34 <- readRDS('/local/yanzijun/CRU/NB/res/RDS/pbmc_T34.RDS')
pbmc_T40 <- readRDS('/local/yanzijun/CRU/NB/res/RDS/pbmc_T40.RDS')
pbmc_T44 <- readRDS('/local/yanzijun/CRU/NB/res/RDS/pbmc_T44.RDS')
pbmc_T69 <- readRDS('/local/yanzijun/CRU/NB/res/RDS/pbmc_T69.RDS')
pbmc_T71 <- readRDS('/local/yanzijun/CRU/NB/res/RDS/pbmc_T71.RDS')
pbmc_T75 <- readRDS('/local/yanzijun/CRU/NB/res/RDS/pbmc_T75.RDS')
pbmc_T92 <- readRDS('/local/yanzijun/CRU/NB/res/RDS/pbmc_T92.RDS')
pbmc_T162 <- readRDS('/local/yanzijun/CRU/NB/res/RDS/pbmc_T162.RDS')
pbmc_T174 <- readRDS('/local/yanzijun/CRU/NB/res/RDS/pbmc_T174.RDS')
pbmc_T188 <- readRDS('/local/yanzijun/CRU/NB/res/RDS/pbmc_T188.RDS')
pbmc_T200 <- readRDS('/local/yanzijun/CRU/NB/res/RDS/pbmc_T200.RDS')
pbmc_T214 <- readRDS('/local/yanzijun/CRU/NB/res/RDS/pbmc_T214.RDS')
pbmc_T230 <- readRDS('/local/yanzijun/CRU/NB/res/RDS/pbmc_T230.RDS')

clin.info <- data.frame(sample=c('T10','T19','T27','T34','T40','T44','T69','T71','T75','T92','T162',
                                 'T174','T188','T200','T214','T230'),
                        gender=c('M','F','M','F','M','F','M','M','F','F','M','M','M','F','F','M'),
                        hist=c(rep('GNB',2),rep('NB',14)),
                        risk=c('H','I','H','H','I','L','H','H','H','H','H','L','L','H','H','H'),
                        MYCN=c(rep('No',10),'Yes','No','No','Yes','No','Yes'))
saveRDS(clin.info,'/local/yanzijun/CRU/NB/data/CancerCell_2020/RDS/clin.info.RDS')

pbmc <- merge(pbmc_T10,
              y = list(pbmc_T19,pbmc_T27,pbmc_T34,pbmc_T40,pbmc_T44,
                       pbmc_T69,pbmc_T71,pbmc_T75,pbmc_T92,pbmc_T162,
                       pbmc_T174,pbmc_T188,pbmc_T200,pbmc_T214,pbmc_T230),
              add.cell.ids = clin.info$sample, project = "NB")


####
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 5000)
pbmc <- ScaleData(object = pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunPCA(object = pbmc, seed.use=123, npcs=50,
               features = VariableFeatures(object = pbmc), ndims.print=1,nfeatures.print=1)
pbmc <- RunTSNE(pbmc, dims = 1:50, seed.use = 123,n.components=2)
pbmc <- RunUMAP(pbmc, dims = 1:50, seed.use = 123,n.components=2)
saveRDS(pbmc,'/local/yanzijun/CRU/NB/res/RDS/pbmc.RDS')

###########
## 2. select High Risk samples
###########
select.ID <- clin.info$sample[clin.info$risk=='H']
Idents(pbmc) <- pbmc$sample
sub.pbmc <- subset(pbmc,idents=select.ID)
saveRDS(sub.pbmc,'/local/yanzijun/CRU/NB/data/CancerCell_2020/RDS/pbmc_H.RDS')

# sub.pbmc <- readRDS('/local/yanzijun/CRU/NB/data/CancerCell_2020/RDS/pbmc_H.RDS')
# EXP <- sub.pbmc@assays$RNA@counts
# write.csv(EXP,'/local/yanzijun/CRU/NB/data/CancerCell_2020/RDS/GSE137804_HRNB_expression.csv')


##########
## 3. Find marker genes
##########
plan("multisession", workers = 16) ###set the compute core
options(future.globals.maxSize = 60000 * 1024^2)
Idents(sub.pbmc) <- sub.pbmc$celltype
table(sub.pbmc$celltype)

MK <- FindAllMarkers(sub.pbmc,only.pos = T)
write.csv(MK,file="/local/yanzijun/CRU/NB_FN/res/GeneSet/scMK/markers_H.csv")


##########
##4. plot
##########
rm(list=ls())
library(Seurat)
library(ggplot2)
color.cell <- c("#1F78B4","#E31A1C","#A6CEE3" ,"#CAB2D6","#33A02C",
                "#FF7F00" ,"#6A3D9A","#B96B63")

## Fig1A: plot MK from HR scRNA-seq
pbmc <- readRDS('/local/yanzijun/CRU/NB/data/CancerCell_2020/RDS/pbmc_H.RDS')
pbmc$celltype <- factor(pbmc$celltype,levels = unique(pbmc$celltype))
Idents(pbmc) <- pbmc$celltype

umap.CT <- DimPlot(pbmc,group.by = 'celltype',cols = color.cell,raster=FALSE)
pdf('/local/yanzijun/CRU/NB_FN/FIG_scMK/umap_CT.pdf')
print(umap.CT)
dev.off()


## Fig1B: heatmap of MK 
MK <- read.csv("/local/yanzijun/CRU/NB_FN/res/GeneSet/scMK/markers_H.csv")
MK[1:2,]
table(MK$cluster)
MK.T <- MK$gene[MK$cluster=='tumor' & MK$p_val_adj<0.05 ]

library(dplyr)
sel.genes <- MK %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

tmp <- AverageExpression(pbmc, return.seurat = TRUE)
tmp@active.ident <- factor(tmp@active.ident, levels = levels(pbmc$celltype))

plot.heatmap <- 
  DoHeatmap(tmp, features = sel.genes$gene,draw.lines = FALSE, 
            label = T, group.colors = color.cell, size = 3.5,
            slot = 'scale.data', disp.min = -2, disp.max = 2) +
  guides(color = F) + 
  labs(fill = 'Row Z-score') + 
  scale_fill_gradientn(colors = c("navy", "white", "firebrick3"), 
                       breaks = c(-1, 0, 1)) +
  theme(#axis.text = element_blank(),
    legend.text = element_text(size = 8,colour = 'black'),
    legend.title = element_text(size = 8,colour = 'black'),
    legend.position = 'bottom', legend.direction = 'horizontal',
    legend.key.size = unit(0.2, "cm"),
    legend.key.width = unit(0.2,"cm") ,
    plot.margin = unit(c(0.5,0,0,0), "cm"))

ggsave('/local/yanzijun/CRU/NB_FN/FIG_scMK/Heatmap_celltype.pdf',plot.heatmap,
       height = 12, width = 10,units = 'cm')



## Fig 1C/D: Enrichment
data(c2BroadSets)
library(Biobase)
library(genefilter)
library(limma)
library(AnnotationHub)
library(org.Hs.eg.db)  
library(clusterProfiler)
library(DOSE)
library(dplyr)
library(tidyverse)
library(reshape2)
library(xlsx)

symbol2id=mapIds(org.Hs.eg.db,MK.T,"ENTREZID",'SYMBOL')
id=symbol2id[which(symbol2id!='')]

#GO BP#
ego <- enrichGO(OrgDb="org.Hs.eg.db", gene = id, ont = "BP", pvalueCutoff = 0.05, readable= TRUE) #GO富集分析
ego_sim <- clusterProfiler::simplify(ego,cutoff=0.6,by="p.adjust",select_fun=min,measure="Wang")
ego_res <- as.data.frame(ego_sim)
dim(ego_res)
if(nrow(ego_res)>0){
  write.xlsx(ego_res,"/local/yanzijun/CRU/NB_FN/res_scMK/Enrich/scMK/GO_BP_all.xlsx",row.names = FALSE,sheetName = 'GOBP',append = TRUE)
}    

#KEGG#
ekk <- enrichKEGG(gene= id,organism  = 'hsa')	 
ekk_res <- as.data.frame(ekk)
if(nrow(ekk_res)>0){
  write.xlsx(ekk_res,"/local/yanzijun/CRU/NB_FN/res_scMK/Enrich/scMK/KEGG.xlsx",row.names = FALSE,sheetName = 'KEGG',append = TRUE)
}

### plot
plot.df <- ego_res
print(plot.df$Description)
plot.df <- plot.df[-9,]
plot.df$log10Qval <- -log(plot.df$qvalue,10)
plot.df <- plot.df[order(plot.df$log10Qval,decreasing = F),]
dim(plot.df)
if(nrow(plot.df)>20){
  plot.df <- tail(plot.df,20)
}

plot.df$Description <- factor(plot.df$Description,levels = plot.df$Description)
GO.p <- ggplot(plot.df, aes(x = Description, y = log10Qval)) + 
  geom_bar(stat = 'identity', color = "white", fill = "#D6221E",width = 0.9,alpha = 0.3) + 
  theme_classic() + coord_flip() +
  labs(x='',y = expression(paste("-log"[10], "(", italic("Q"), "-value)"))) +
  theme(axis.title = element_text(size = 8, colour = 'black'), 
        axis.text.y = element_text(size = 8, colour = 'black'), 
        axis.text.x = element_text(size = 8, colour = 'black'))
GO.p
ggsave('/local/yanzijun/CRU/NB_FN/FIG_scMK/GO_scMK.pdf',GO.p,height = 9, width = 10,units = 'cm')


plot.df <- ekk_res[c(4,8,20,21,24,25,28,33,34,37,38,43,48,49,55:56,60,64),]
print(plot.df$Description)
plot.df$log10Qval <- -log(plot.df$qvalue,10)
plot.df <- plot.df[order(plot.df$log10Qval,decreasing = F),]
dim(plot.df)
if(nrow(plot.df)>20){
  plot.df <- tail(plot.df,20)
}

plot.df$Description <- factor(plot.df$Description,levels = plot.df$Description)
KEGG.p <- ggplot(plot.df, aes(x = Description, y = log10Qval)) + 
  geom_bar(stat = 'identity', color = "white", fill = "#148CEF",width = 0.9,alpha = 0.4) + 
  theme_classic() + coord_flip() +
  labs(x='',y = expression(paste("-log"[10], "(", italic("Q"), "-value)"))) +
  theme(axis.title = element_text(size = 8, colour = 'black'), 
        axis.text.y = element_text(size = 8, colour = 'black'), 
        axis.text.x = element_text(size = 8, colour = 'black'))
KEGG.p
ggsave('/local/yanzijun/CRU/NB_FN/FIG_scMK/KEGG_scMK.pdf',KEGG.p,height = 9, width = 9.5,units = 'cm')


## sFig 3C: expression of six-signature genes in each celltype
sub.pbmc <- readRDS('/local/yanzijun/CRU/NB/data/CancerCell_2020/RDS/pbmc_H.RDS')
dim(sub.pbmc)
Idents(sub.pbmc) <- sub.pbmc$celltype
table(sub.pbmc$celltype)
geneLst <- c('RAMP1','CDT1','MEG3','NPW','MAPT','C1QTNF4')
p <- DotPlot(sub.pbmc,features = geneLst,scale = TRUE)
ggsave('/local/yanzijun/CRU/NB_FN/FIG_scMK/Dotplot.pdf',p,height = 6,width = 18,units = 'cm')