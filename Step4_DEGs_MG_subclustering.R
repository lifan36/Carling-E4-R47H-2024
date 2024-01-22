#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)
setwd("/athena/ganlab/scratch/lif4001/LG31E3E4/integration")
Cluster_MG <- readRDS("LG31E3E4_MG_subset.rds")
Cluster_OL <- readRDS("LG31E3E4_OL_subset.rds")
Cluster_AST <- readRDS("LG31E3E4_AST_subset.rds")

Idents(Cluster_MG) <- "Condition"
Idents(Cluster_AST) <- "Condition"
Idents(Cluster_OL) <- "Condition"

setwd("/athena/ganlab/scratch/lif4001/LG31E3E4/integration/DEGs")
df <- FindMarkers(Cluster_MG, ident.1 = "E3_P301S_R47H", ident.2 = "E3_P301S_WT", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
df <- df[df$p_val_adj < 0.05,]
write.csv(df, "MG_E3_P301S_R47H_vs_E3_P301S_WT_DEGs.csv")
df <- FindMarkers(Cluster_AST, ident.1 = "E3_P301S_R47H", ident.2 = "E3_P301S_WT", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
df <- df[df$p_val_adj < 0.05,]
write.csv(df, "AST_E3_P301S_R47H_vs_E3_P301S_WT_DEGs.csv")
df <- FindMarkers(Cluster_OL, ident.1 = "E3_P301S_R47H", ident.2 = "E3_P301S_WT", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
df <- df[df$p_val_adj < 0.05,]
write.csv(df, "OL_E3_P301S_R47H_vs_E3_P301S_WT_DEGs.csv")
df <- FindMarkers(Cluster_MG, ident.1 = "E4_P301S_R47H", ident.2 = "E4_P301S_WT", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
df <- df[df$p_val_adj < 0.05,]
write.csv(df, "MG_E4_P301S_R47H_vs_E4_P301S_WT_DEGs.csv")
df <- FindMarkers(Cluster_AST, ident.1 = "E4_P301S_R47H", ident.2 = "E4_P301S_WT", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
df <- df[df$p_val_adj < 0.05,]
write.csv(df, "AST_E4_P301S_R47H_vs_E4_P301S_WT_DEGs.csv")
df <- FindMarkers(Cluster_OL, ident.1 = "E4_P301S_R47H", ident.2 = "E4_P301S_WT", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
df <- df[df$p_val_adj < 0.05,]
write.csv(df, "OL_E4_P301S_R47H_vs_E4_P301S_WT_DEGs.csv")

#E4 subset and subclustering
E4_MG <- subset(Cluster_MG, idents=c("E4_NTG_WT","E4_P301S_WT","E4_P301S_R47H"))
E4_MG <- E4_MG[!grepl("^mt-", rownames(E4_MG)), ]
DefaultAssay(E4_MG) <- 'integrated'
E4_MG <- ScaleData(E4_MG, verbose = FALSE)
E4_MG <- RunPCA(E4_MG, features = VariableFeatures(object = E4_MG), verbose = FALSE)
ElbowPlot(E4_MG)
E4_MG <- FindNeighbors(E4_MG, dims = 1:15)
E4_MG <- FindClusters(E4_MG, resolution = 0.3)
E4_MG <- RunUMAP(E4_MG, dims = 1: 15)
# rename cluster
n <- dim(table(E4_MG@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
E4_MG@active.ident <- plyr::mapvalues(x = E4_MG@active.ident, from = current.cluster.ids, to = new.cluster.ids)
E4_MG@active.ident <- factor(E4_MG@active.ident, levels=1:n)

pdf("E4_MG_umap_Condition.pdf", width=8.2, height=2.7)
DimPlot(E4_MG, reduction = "umap", split.by = "Condition", label = T, ncol = 3)
dev.off()

DefaultAssay(E4_MG) <- "RNA"
pdf("E4_MG_VlnPlot.pdf", width=10, height=2)
VlnPlot(E4_MG, features = c("P2ry12","Stat1","mt-Co1","Mrc1","Skap1"), pt.size = 0, ncol = 5)
dev.off()

DefaultAssay(E4_MG) <- 'RNA'
LG31E3E4_MG_markers <- FindAllMarkers(E4_MG, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
test <- LG31E3E4_MG_markers
test <- test[test$p_val_adj < 0.05,]
test <- test[order(-test$avg_log2FC),]
df <- createWorkbook()
cluster_save <- function(cluster=7, data=test){
  for (i in seq(1, cluster)) {
    addWorksheet(df, sheetName = paste(i, sep=""))
    writeData(df, sheet = (paste(i, sep="")), x = (data %>% filter(cluster == i)))
  }
}
cluster_save(cluster=7, data=test)
saveWorkbook(df, "E4_MG_markers.xlsx")

E4_MG <- RenameIdents(E4_MG, `1` = "Homeostatic", `2`="DAM", `3`="DAM transition", `4`="Interferon", `5`="Macrophage", `6`="T cell")
E4_MG$state <- Idents(E4_MG)
E4_MG$state <- factor(x = E4_MG$state, levels = c("Homeostatic","DAM transition","DAM","Interferon","Macrophage","T cell"))
# Figure 2C
pdf("E4_MG_umap.pdf", width=4.5, height=2.7)
DimPlot(E4_MG, reduction = "umap", label = T)
dev.off()
# Figure 2D
pdf("E4_MG_umap_annotation_Condition.pdf", width=8, height=2.7)
DimPlot(E4_MG, reduction = "umap", split.by = "Condition", label = F)
dev.off()

Idents(E4_MG) <- "state"
DefaultAssay(E4_MG) <- "RNA"
pdf("E4_MG_HeatMap_1.pdf", width=12, height=6)
DoHeatmap(E4_MG, features = c("P2ry12","Cx3cr1","Tmem119","Selplg","Siglech",
                                 "Cd83","Myo1e","Arhgap24","Cacna1a","Fgr",
                                 "Apoe","Spp1","Csf1","Lpl","Cd74",
                                 "Oasl2","Ifi204","Rnf213","Stat1","Trim30a",
                                 "F13a1","Mrc1","Dab2","Rbpj","Cd163",
                                 "Skap1","Themis","Grap2","Cd247","Stat4"
                                 
)) + NoLegend()
dev.off()

DefaultAssay(E4_MG) <- "RNA"
pdf("E4_MG_DotPlot.pdf", width=10, height=4)
DotPlot(E4_MG, features = c("P2ry12","Cx3cr1","Tmem119","Selplg","Siglech",
                               "Cd83","Myo1e","Arhgap24","Cacna1a","Fgr",
                               "Apoe","Spp1","Csf1","Lpl","Cd74",
                               "Oasl2","Ifi204","Rnf213","Stat1","Trim30a",
                               "F13a1","Mrc1","Dab2","Rbpj","Cd163",
                               "Skap1","Themis","Grap2","Cd247","Stat4"
                               
)) + RotatedAxis()
dev.off()

# Figure 2E
pdf("E4_MG_VlnPlot.pdf", width=12, height=3.2)
VlnPlot(E4_MG, features = c("P2ry12","Myo1e","Stat1","Mrc1","Skap1"), pt.size = 0, ncol = 5, cols = c('Homeostatic' = '#00BFC4', 'DAM transition' = '#B79F00','DAM'='#00BA38',
                                                                                                         'Interferon' = '#F8766D','Macrophage'='#619CFF','T cell'='#F564E3'))
dev.off()

# Source data for Figure 2F
write.csv(table(E4_MG$seurat_clusters, E4_MG$Sample_Name), "E4_MG_subcluster_cell_counts.csv")

DefaultAssay(E4_MG) <- 'RNA'
LG31E3E4_MG_markers <- FindAllMarkers(E4_MG, logfc.threshold = 0.1, test.use = "MAST", only.pos = F)
write.csv(LG31E3E4_MG_markers, "LG31E3E4_MG_DEGs.csv")
test <- LG31E3E4_MG_markers
test <- test[test$p_val_adj < 0.05,]
test <- test[order(-test$avg_log2FC),]
df <- createWorkbook()
cluster_save <- function(cluster=7, data=test){
  for (i in seq(1, cluster)) {
    addWorksheet(df, sheetName = paste(i, sep=""))
    writeData(df, sheet = (paste(i, sep="")), x = (data %>% filter(cluster == i)))
  }
}
cluster_save(cluster=7, data=test)
saveWorkbook(df, "E4_MG_DEGs.xlsx")


Idents(E4_MG) <- "orig.ident"
MG_gene_expressions <- AverageExpression(E4_MG, features = c("Trim30d", "Trim30a", "Stat1","Sp100", "Slfn8", "Rnf213", "Parp14",
                                                                "Nlrc5", "Herc6", "Eif2ak2", "Ddx60","Gm4951", "Apobec3","Cd300lf"), assays = 'RNA')
write.csv(MG_gene_expressions$RNA, file = 'average_expression_of_interferon_genes_in_MG_by_orig.ident.csv')


saveRDS(E4_MG, file = "LG31E3E4_E4_MG_03132023.rds")

#E3 subset and subclustering
E3_MG <- subset(Cluster_MG, idents=c("E3_NTG_WT","E3_P301S_WT","E3_P301S_R47H"))
E3_MG <- E3_MG[!grepl("^mt-", rownames(E3_MG)), ]
DefaultAssay(E3_MG) <- 'integrated'
E3_MG <- ScaleData(E3_MG, verbose = FALSE)
E3_MG <- RunPCA(E3_MG, features = VariableFeatures(object = E3_MG), verbose = FALSE)
ElbowPlot(E3_MG)
E3_MG <- FindNeighbors(E3_MG, dims = 1:11)
E3_MG <- FindClusters(E3_MG, resolution = 0.15)
E3_MG <- RunUMAP(E3_MG, dims = 1: 11)
# rename cluster
n <- dim(table(E3_MG@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
E3_MG@active.ident <- plyr::mapvalues(x = E3_MG@active.ident, from = current.cluster.ids, to = new.cluster.ids)
E3_MG@active.ident <- factor(E3_MG@active.ident, levels=1:n)

pdf("E3_MG_umap_Condition.pdf", width=8.2, height=2.7)
DimPlot(E3_MG, reduction = "umap", split.by = "Condition", label = T, ncol = 3)
dev.off()

DefaultAssay(E3_MG) <- "RNA"
pdf("E3_MG_VlnPlot.pdf", width=10, height=2)
VlnPlot(E3_MG, features = c("P2ry12","Stat1","mt-Co1","Mrc1","Skap1"), pt.size = 0, ncol = 5)
dev.off()

FeaturePlot(E3_MG, features = "Stat1")

DefaultAssay(E3_MG) <- "RNA"
Idents(E3_MG) <- "Condition"
pdf("E3_MG_DotPlot_Condition.pdf", width=8, height=3)
DotPlot(E3_MG, features = c("Trim30d", "Trim30a", "Stat1","Sp100", "Slfn8", "Rnf213", "Parp14",
                               "Nlrc5", "Herc6", "Eif2ak2", "Ddx60","Gm4951", "Apobec3","Cd300lf"
                               
)) + RotatedAxis() + scale_colour_gradient2(low = "darkblue", mid = "white", high = "darkred")
dev.off()

Idents(E3_MG) <- "seurat_clusters"
# rename cluster
n <- dim(table(E3_MG@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
E3_MG@active.ident <- plyr::mapvalues(x = E3_MG@active.ident, from = current.cluster.ids, to = new.cluster.ids)
E3_MG@active.ident <- factor(E3_MG@active.ident, levels=1:n)

DefaultAssay(E3_MG) <- 'RNA'
LG31E3E4_MG_markers <- FindAllMarkers(E3_MG, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
test <- LG31E3E4_MG_markers
test <- test[test$p_val_adj < 0.05,]
test <- test[order(-test$avg_log2FC),]
df <- createWorkbook()
cluster_save <- function(cluster=7, data=test){
  for (i in seq(1, cluster)) {
    addWorksheet(df, sheetName = paste(i, sep=""))
    writeData(df, sheet = (paste(i, sep="")), x = (data %>% filter(cluster == i)))
  }
}
cluster_save(cluster=7, data=test)
saveWorkbook(df, "E3_MG_markers.xlsx")

# Figure S3A
pdf("E3_MG_umap_annotation_Condition_new_color.pdf", width=8, height=2.7)
DimPlot(E3_MG, reduction = "umap", split.by = "Condition", label = F, cols = c('1' = '#00BFC4', '3' = '#B79F00','2'='#00BA38',
                                                                            '4'='#619CFF','5'='#F564E3'))
dev.off()

# Figure S3C
pdf("E3_MG_VlnPlot.pdf", width=12, height=3.2)
VlnPlot(E3_MG, features = c("P2ry12","Myo1e","Cst3","Mrc1","Skap1"), pt.size = 0, ncol = 5, cols = c('1' = '#00BFC4', '3' = '#B79F00','2'='#00BA38',
                                                                                                  '4'='#619CFF','5'='#F564E3'))
dev.off()

# Source data for Figure S3B
write.csv(table(E3_MG$seurat_clusters, E3_MG$Sample_Name), "E3_MG_subcluster_cell_counts.csv")

saveRDS(E3_MG, file = "LG31E3E4_E3_MG_03132023.rds" )






