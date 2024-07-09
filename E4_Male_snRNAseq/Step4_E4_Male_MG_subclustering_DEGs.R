#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)

setwd("/Users/lifan/Desktop/data_analysis/LG31E3E4/male/integration_E4")
male_MG <- readRDS("male_MG_subset.rds")

male_MG <- male_MG[!grepl("^mt-", rownames(male_MG)), ]
DefaultAssay(male_MG) <- 'integrated'
male_MG <- ScaleData(male_MG, verbose = FALSE)
male_MG <- RunPCA(male_MG, features = VariableFeatures(object = male_MG), verbose = FALSE)
ElbowPlot(male_MG)
male_MG <- FindNeighbors(male_MG, dims = 1:15)
male_MG <- FindClusters(male_MG, resolution = 0.2)
male_MG <- RunUMAP(male_MG, dims = 1: 15)
# rename cluster
n <- dim(table(male_MG@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
male_MG@active.ident <- plyr::mapvalues(x = male_MG@active.ident, from = current.cluster.ids, to = new.cluster.ids)
male_MG@active.ident <- factor(male_MG@active.ident, levels=1:n)

# 
write.csv(table(male_MG$seurat_clusters, male_MG$Sample_Name), "male_MG_subcluster_cell_counts.csv")

df <- FindAllMarkers(male_MG, logfc.threshold = 0.25, min.pct = 0.25, only.pos = T,test.use = "MAST")
df <- df[df$p_val_adj < 0.05,]
write.csv(df, "male_MG_markers.csv")

# remove cluster 6 (154 cells, 2% of total MG, mainly from 2 samples of M_E4_P301S_WT), cluster 7 (73 cells, 1% of total MG, expressing neuronal genes)

male_MG <- subset(male_MG, idents=c("6","7"), invert=T)

Idents(male_MG) <- "seurat_clusters"
# rename cluster
n <- dim(table(male_MG@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
male_MG@active.ident <- plyr::mapvalues(x = male_MG@active.ident, from = current.cluster.ids, to = new.cluster.ids)
male_MG@active.ident <- factor(male_MG@active.ident, levels=1:n)

# Figure S4E
pdf("male_MG_umap_annotation_Condition_new_color_new.pdf", width=7, height=2.66)
DimPlot(male_MG, reduction = "umap", split.by = "Condition", label = T, cols = c('1' = '#00BFC4', '3' = '#B79F00','2'='#00BA38',
                                                                               '4'='#619CFF','5'='#F564E3'))
dev.off()

# Figure S4G
DefaultAssay(male_MG) <- "RNA"
pdf("male_MG_DotPlot_Condition_nomt_with_Siglech_Hexb.pdf.pdf", width=10, height=2.5)
VlnPlot(male_MG, features = c("P2ry12","Myo1e","Cst3","Mrc1","Skap1"), pt.size = 0, ncol = 5, cols = c('1' = '#00BFC4', '3' = '#B79F00','2'='#00BA38',
                                                                                                     '4'='#619CFF','5'='#F564E3'))
dev.off()

DefaultAssay(male_MG) <- "RNA"
Idents(male_MG) <- "Condition"

# Figure S4H
pdf("male_MG_IFN_DotPlot_Condition_new.pdf", width=7, height=3)
DotPlot(male_MG, features = c( "Pik3ap1", "B2m", "Trim14", "Rnf213", "Nlrc5", "Sp100", "Sp110", "Usp25",
                          "Samd9l", "Ifi204", "Gm4951", "Adar"
)) + RotatedAxis() + scale_colour_gradient2(low = "darkblue", mid = "white", high = "darkred")
dev.off()

#.	DEG list for E4-R47H-P301S vs. E4-WT-P301S in male mice

df <- FindMarkers(male_MG, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = F, ident.1 = "M_E4_P301S_R47H",
                  ident.2 = "M_E4_P301S_WT")
df <- df[df$p_val_adj < 0.05,]
write.csv(df, "male_E4_MG_M_E4_P301S_R47H_vs_M_E4_P301S_WT_DEGs.csv")





