#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)

setwd("/athena/ganlab/scratch/lif4001/LG31E3E4/DF_2ndRound")
E3_NTG_WT_1 <- readRDS(file = "E3_NTG_WT_1_singlets_PCA.rds")
E3_NTG_WT_2 <- readRDS(file = "E3_NTG_WT_2_singlets_PCA.rds")
E3_NTG_WT_3 <- readRDS(file = "E3_NTG_WT_3_singlets_PCA.rds")

E3_P301S_WT_1 <- readRDS(file = "E3_P301S_WT_1_singlets_PCA.rds")
E3_P301S_WT_2 <- readRDS(file = "E3_P301S_WT_2_singlets_PCA.rds")
E3_P301S_WT_3 <- readRDS(file = "E3_P301S_WT_3_singlets_PCA.rds")

E3_P301S_R47H_1 <- readRDS(file = "E3_P301S_R47H_1_singlets_PCA.rds")
E3_P301S_R47H_2 <- readRDS(file = "E3_P301S_R47H_2_singlets_PCA.rds")
E3_P301S_R47H_3 <- readRDS(file = "E3_P301S_R47H_3_singlets_PCA.rds")

E4_NTG_WT_1 <- readRDS(file = "E4_NTG_WT_1_singlets_PCA.rds")
E4_NTG_WT_2 <- readRDS(file = "E4_NTG_WT_2_singlets_PCA.rds")
E4_NTG_WT_3 <- readRDS(file = "E4_NTG_WT_3_singlets_PCA.rds")

E4_P301S_WT_1 <- readRDS(file = "E4_P301S_WT_1_singlets_PCA.rds")
E4_P301S_WT_2 <- readRDS(file = "E4_P301S_WT_2_singlets_PCA.rds")
E4_P301S_WT_3 <- readRDS(file = "E4_P301S_WT_3_singlets_PCA.rds")

E4_P301S_R47H_1 <- readRDS(file = "E4_P301S_R47H_1_singlets_PCA.rds")
E4_P301S_R47H_2 <- readRDS(file = "E4_P301S_R47H_2_singlets_PCA.rds")
E4_P301S_R47H_3 <- readRDS(file = "E4_P301S_R47H_3_singlets_PCA.rds")

setwd("/athena/ganlab/scratch/lif4001/LG31E3E4/integration_4more")
E3_NTG_WT <- c(E3_NTG_WT_1, E3_NTG_WT_2, E3_NTG_WT_3)
anchors_E3_NTG_WT <- FindIntegrationAnchors(object.list = E3_NTG_WT, dims = 1:30)
E3_NTG_WT_integrated <- IntegrateData(anchorset = anchors_E3_NTG_WT, dims = 1:30)
rm(E3_NTG_WT_1, E3_NTG_WT_2, E3_NTG_WT_3, E3_NTG_WT)

E3_P301S_WT <- c(E3_P301S_WT_1, E3_P301S_WT_2, E3_P301S_WT_3)
anchors_E3_P301S_WT <- FindIntegrationAnchors(object.list = E3_P301S_WT, dims = 1:30)
E3_P301S_WT_integrated <- IntegrateData(anchorset = anchors_E3_P301S_WT, dims = 1:30)
rm(E3_P301S_WT_1, E3_P301S_WT_2, E3_P301S_WT_3, E3_P301S_WT)

E3_P301S_R47H <- c(E3_P301S_R47H_1, E3_P301S_R47H_2, E3_P301S_R47H_3)
anchors_E3_P301S_R47H <- FindIntegrationAnchors(object.list = E3_P301S_R47H, dims = 1:30)
E3_P301S_R47H_integrated <- IntegrateData(anchorset = anchors_E3_P301S_R47H, dims = 1:30)
rm(E3_P301S_R47H_1, E3_P301S_R47H_2, E3_P301S_R47H_3, E3_P301S_R47H)

E4_NTG_WT <- c(E4_NTG_WT_1, E4_NTG_WT_2, E4_NTG_WT_3)
anchors_E4_NTG_WT <- FindIntegrationAnchors(object.list = E4_NTG_WT, dims = 1:30)
E4_NTG_WT_integrated <- IntegrateData(anchorset = anchors_E4_NTG_WT, dims = 1:30)
rm(E4_NTG_WT_1, E4_NTG_WT_2, E4_NTG_WT_3, E4_NTG_WT)

E4_P301S_WT <- c(E4_P301S_WT_1, E4_P301S_WT_2, E4_P301S_WT_3)
anchors_E4_P301S_WT <- FindIntegrationAnchors(object.list = E4_P301S_WT, dims = 1:30)
E4_P301S_WT_integrated <- IntegrateData(anchorset = anchors_E4_P301S_WT, dims = 1:30)
rm(E4_P301S_WT_1, E4_P301S_WT_2, E4_P301S_WT_3, E4_P301S_WT)

E4_P301S_R47H <- c(E4_P301S_R47H_1, E4_P301S_R47H_2, E4_P301S_R47H_3)
anchors_E4_P301S_R47H <- FindIntegrationAnchors(object.list = E4_P301S_R47H, dims = 1:30)
E4_P301S_R47H_integrated <- IntegrateData(anchorset = anchors_E4_P301S_R47H, dims = 1:30)
rm(E4_P301S_R47H_1, E4_P301S_R47H_2, E4_P301S_R47H_3, E4_P301S_R47H)

LG31E3E4 <- c(E3_NTG_WT_integrated, E3_P301S_WT_integrated, E3_P301S_R47H_integrated, E4_NTG_WT_integrated, E4_P301S_WT_integrated, E4_P301S_R47H_integrated)
anchors_LG31E3E4 <- FindIntegrationAnchors(object.list = LG31E3E4, dims = 1:30)
LG31E3E4_integrated <- IntegrateData(anchorset = anchors_LG31E3E4, dims = 1:30)
rm(E3_NTG_WT_integrated, E3_P301S_WT_integrated, E3_P301S_R47H_integrated, E4_NTG_WT_integrated, E4_P301S_WT_integrated, E4_P301S_R47H_integrated, LG31E3E4)


#saveRDS(LG31E3E4_integrated, file = "LG31E3E4_integrated.rds")

DefaultAssay(LG31E3E4_integrated) <- 'integrated'

# LG31E3E4_integrated <- NormalizeData(LG31E3E4_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
# LG31E3E4_integrated <- FindVariableFeatures(LG31E3E4_integrated, selection.method = "vst", nfeatures = 3000)

LG31E3E4_integrated <- ScaleData(LG31E3E4_integrated, verbose = FALSE)
LG31E3E4_integrated <- RunPCA(LG31E3E4_integrated, features = VariableFeatures(object = LG31E3E4_integrated), verbose = FALSE)

LG31E3E4_integrated <- FindNeighbors(LG31E3E4_integrated, dims = 1:15)
LG31E3E4_integrated <- FindClusters(LG31E3E4_integrated, resolution = 0.1)
LG31E3E4_integrated <- RunUMAP(LG31E3E4_integrated, dims = 1: 15)

DefaultAssay(LG31E3E4_integrated) <- 'RNA'
LG31E3E4_integrated <- NormalizeData(LG31E3E4_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
LG31E3E4_integrated <- ScaleData(LG31E3E4_integrated, features = rownames(LG31E3E4_integrated))

#saveRDS(LG31E3E4_integrated, file = 'LG31E3E4_integrated_PCA_0.1.rds')
#LG31E3E4_integrated <- readRDS(file = "LG31E3E4_integrated_PCA_0.1.rds")

LG31E3E4_integrated$Condition <- factor(x = LG31E3E4_integrated$Condition, levels = c("E3_NTG_WT","E3_P301S_WT","E3_P301S_R47H","E4_NTG_WT","E4_P301S_WT","E4_P301S_R47H"))
LG31E3E4_integrated$Sample_Name <- factor(x = LG31E3E4_integrated$Sample_Name, levels = c("E3_NTG_WT_1","E3_NTG_WT_2","E3_NTG_WT_3",
                                                                                    "E3_P301S_WT_1","E3_P301S_WT_2","E3_P301S_WT_3",
                                                                                    "E3_P301S_R47H_1","E3_P301S_R47H_2","E3_P301S_R47H_3",
                                                                                    "E4_NTG_WT_1","E4_NTG_WT_2","E4_NTG_WT_3",
                                                                                    "E4_P301S_WT_1","E4_P301S_WT_2","E4_P301S_WT_3",
                                                                                    "E4_P301S_R47H_1","E4_P301S_R47H_2","E4_P301S_R47H_3"))

pdf("LG31E3E4_QC.pdf", width=9, height=4)
Idents(LG31E3E4_integrated) <- "Condition"
VlnPlot(object = LG31E3E4_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()

# Figure S2E,F,G
Idents(LG31E3E4_integrated) <- "Sample_Name"
pdf("LG31E3E4_QC_Sample.pdf", width=12, height=4)

VlnPlot(object = LG31E3E4_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()

Idents(LG31E3E4_integrated) <- "seurat_clusters"
pdf("LG31E3E4_integrated_umap.pdf", width=5, height=4)
DimPlot(LG31E3E4_integrated, reduction = 'umap', label = T)
dev.off()
pdf("LG31E3E4_integrated_umap_split_individual.pdf", width=15, height=12)
DimPlot(LG31E3E4_integrated, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 5)
dev.off()
pdf("LG31E3E4_integrated_umap_split_Condition.pdf", width=9, height=6)
DimPlot(LG31E3E4_integrated, reduction = "umap", split.by = "Condition", label = T, ncol = 3)
dev.off()

write.csv(table(LG31E3E4_integrated$seurat_clusters, LG31E3E4_integrated$Sample_Name), "LG31E3E4_cell_counts_cluster_by_sample.csv")

DefaultAssay(LG31E3E4_integrated) <- 'RNA'

LG31E3E4_markers <- FindAllMarkers(LG31E3E4_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1, test.use = "MAST")
write.csv(LG31E3E4_markers, "LG31E3E4_markers.csv")


saveRDS(LG31E3E4_integrated, file = 'LG31E3E4_integrated_PCA_0.1.rds')

LG31E3E4_markers <- read.csv(file = "LG31E3E4_markers.csv", header=T,row.names =1)
top5 <- LG31E3E4_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5$gene <- as.character(top5$gene)
pdf("LG31E3E4_HeatMapTop5_0.1_new.pdf", width=24, height=16)
DoHeatmap(LG31E3E4_integrated, features = top5$gene) + NoLegend()
dev.off()

#Add marker genes

sig_EN<-c("Snap25","Slc17a7", "Nrgn","Gad1", "Gad2","Plp1", "Mbp", "Mobp","Sntn","Aqp4", "Clu", "Aldoc", "Pla2g7","Cx3cr1", "P2ry12", "Csf1r",
          "Pdgfra", "Vcan", "Flt1","Vtn", "Igfbp7")
markers.to.plot <- as.matrix(sig_EN)
pdf("LG31E3E4_annotation_combine.pdf", width=10, height=5)
DotPlot(object = LG31E3E4_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()
dev.off()

