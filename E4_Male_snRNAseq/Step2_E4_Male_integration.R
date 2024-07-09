
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)

#load in data from Cell Ranger or other counts data ====

#for loading Cell Ranger counts:
setwd("/athena/ganlab/scratch/lif4001/LG31E3E4_male/DF_2ndRound")

LG31E4_R47H_78 <- readRDS(file = "LG31E4_R47H_78_singlets_PCA.rds")
LG31E4_R47H_80 <- readRDS(file = "LG31E4_R47H_80_singlets_PCA.rds")
LG31E4_R47H_94 <- readRDS(file = "LG31E4_R47H_94_singlets_PCA.rds")
LG31E4_R47H_60 <- readRDS(file = "LG31E4_R47H_60_singlets_PCA.rds")
LG31E4_R47H_61 <- readRDS(file = "LG31E4_R47H_61_singlets_PCA.rds")
LG31E4_R47H_82 <- readRDS(file = "LG31E4_R47H_82_singlets_PCA.rds")
LG31E4_R47H_32 <- readRDS(file = "LG31E4_R47H_32_singlets_PCA.rds")
LG31E4_R47H_68 <- readRDS(file = "LG31E4_R47H_68_singlets_PCA.rds")
LG31E4_R47H_101 <- readRDS(file = "LG31E4_R47H_101_singlets_PCA.rds")

LG31E4_R47H_78[["Condition"]] = c('M_E4_NTG_WT')
LG31E4_R47H_78[["Sample_Name"]] = c('M_E4_NTG_WT_1')
LG31E4_R47H_80[["Condition"]] = c('M_E4_NTG_WT')
LG31E4_R47H_80[["Sample_Name"]] = c('M_E4_NTG_WT_2')
LG31E4_R47H_94[["Condition"]] = c('M_E4_NTG_WT')
LG31E4_R47H_94[["Sample_Name"]] = c('M_E4_NTG_WT_3')

LG31E4_R47H_60[["Condition"]] = c('M_E4_P301S_WT')
LG31E4_R47H_60[["Sample_Name"]] = c('M_E4_P301S_WT_1')
LG31E4_R47H_61[["Condition"]] = c('M_E4_P301S_WT')
LG31E4_R47H_61[["Sample_Name"]] = c('M_E4_P301S_WT_2')
LG31E4_R47H_82[["Condition"]] = c('M_E4_P301S_WT')
LG31E4_R47H_82[["Sample_Name"]] = c('M_E4_P301S_WT_3')

LG31E4_R47H_32[["Condition"]] = c('M_E4_P301S_R47H')
LG31E4_R47H_32[["Sample_Name"]] = c('M_E4_P301S_R47H_1')
LG31E4_R47H_68[["Condition"]] = c('M_E4_P301S_R47H')
LG31E4_R47H_68[["Sample_Name"]] = c('M_E4_P301S_R47H_2')
LG31E4_R47H_101[["Condition"]] = c('M_E4_P301S_R47H')
LG31E4_R47H_101[["Sample_Name"]] = c('M_E4_P301S_R47H_3')

setwd("/athena/ganlab/scratch/lif4001/LG31E3E4_male/integration_E4")

M_E4_NTG_WT <- c(LG31E4_R47H_78, LG31E4_R47H_80, LG31E4_R47H_94)
anchors_M_E4_NTG_WT <- FindIntegrationAnchors(object.list = M_E4_NTG_WT, dims = 1:30)
M_E4_NTG_WT_integrated <- IntegrateData(anchorset = anchors_M_E4_NTG_WT, dims = 1:30)
rm(LG31E4_R47H_78, LG31E4_R47H_80, LG31E4_R47H_94, M_E4_NTG_WT)

M_E4_P301S_WT <- c(LG31E4_R47H_60, LG31E4_R47H_61, LG31E4_R47H_82)
anchors_M_E4_P301S_WT <- FindIntegrationAnchors(object.list = M_E4_P301S_WT, dims = 1:30)
M_E4_P301S_WT_integrated <- IntegrateData(anchorset = anchors_M_E4_P301S_WT, dims = 1:30)
rm(LG31E4_R47H_60, LG31E4_R47H_61, LG31E4_R47H_82, M_E4_P301S_WT)

M_E4_P301S_R47H <- c(LG31E4_R47H_32, LG31E4_R47H_68, LG31E4_R47H_101)
anchors_M_E4_P301S_R47H <- FindIntegrationAnchors(object.list = M_E4_P301S_R47H, dims = 1:30)
M_E4_P301S_R47H_integrated <- IntegrateData(anchorset = anchors_M_E4_P301S_R47H, dims = 1:30)
rm(LG31E4_R47H_32, LG31E4_R47H_68, LG31E4_R47H_101, M_E4_P301S_R47H)


male <- c(M_E4_NTG_WT_integrated, M_E4_P301S_WT_integrated, M_E4_P301S_R47H_integrated)
anchors_male <- FindIntegrationAnchors(object.list = male, dims = 1:30)
male_integrated <- IntegrateData(anchorset = anchors_male, dims = 1:30)
rm(M_E4_NTG_WT_integrated, M_E4_P301S_WT_integrated, M_E4_P301S_R47H_integrated, male)

#saveRDS(male_integrated, file = "male_integrated.rds")

DefaultAssay(male_integrated) <- 'integrated'

# male_integrated <- NormalizeData(male_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
# male_integrated <- FindVariableFeatures(male_integrated, selection.method = "vst", nfeatures = 3000)

male_integrated <- ScaleData(male_integrated, verbose = FALSE)
male_integrated <- RunPCA(male_integrated, features = VariableFeatures(object = male_integrated), verbose = FALSE)

male_integrated <- FindNeighbors(male_integrated, dims = 1:15)
male_integrated <- FindClusters(male_integrated, resolution = 0.1)
male_integrated <- RunUMAP(male_integrated, dims = 1: 15)

DefaultAssay(male_integrated) <- 'RNA'
male_integrated <- NormalizeData(male_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
male_integrated <- ScaleData(male_integrated, features = rownames(male_integrated))

#saveRDS(male_integrated, file = 'male_integrated_PCA_0.1.rds')
#male_integrated <- readRDS(file = "male_integrated_PCA_0.1.rds")

male_integrated$Condition <- factor(x = male_integrated$Condition, levels = c("M_E4_NTG_WT","M_E4_P301S_WT","M_E4_P301S_R47H"))
male_integrated$Sample_Name <- factor(x = male_integrated$Sample_Name, levels = c("M_E4_NTG_WT_1","M_E4_NTG_WT_2","M_E4_NTG_WT_3",
                                                                                  "M_E4_P301S_WT_1","M_E4_P301S_WT_2","M_E4_P301S_WT_3",
                                                                                  "M_E4_P301S_R47H_1","M_E4_P301S_R47H_2","M_E4_P301S_R47H_3"))

pdf("male_QC.pdf", width=9, height=4)
Idents(male_integrated) <- "Condition"
VlnPlot(object = male_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()

# Figure S5E, F, G
Idents(male_integrated) <- "Sample_Name"
pdf("male_QC_Sample.pdf", width=22, height=4)

VlnPlot(object = male_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()

Idents(male_integrated) <- "seurat_clusters"
pdf("male_integrated_umap.pdf", width=5, height=4)
DimPlot(male_integrated, reduction = 'umap', label = T)
dev.off()
pdf("male_integrated_umap_split_individual.pdf", width=9, height=9)
DimPlot(male_integrated, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 3)
dev.off()
pdf("male_integrated_umap_split_Condition.pdf", width=9, height=3)
DimPlot(male_integrated, reduction = "umap", split.by = "Condition", label = T, ncol = 3)
dev.off()

write.csv(table(male_integrated$seurat_clusters, male_integrated$Sample_Name), "male_cell_counts_cluster_by_sample.csv")

DefaultAssay(male_integrated) <- 'RNA'

male_markers <- FindAllMarkers(male_integrated, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 1, test.use = "MAST")
write.csv(male_markers, "male_markers.csv")


saveRDS(male_integrated, file = 'male_integrated_PCA_0.1.rds')

#Add marker genes

sig_EN<-c("Snap25","Slc17a7", "Nrgn","Gad1", "Gad2","Plp1", "Mbp", "Mobp","Sntn","Aqp4", "Clu", "Aldoc", "Pla2g7","Cx3cr1", "P2ry12", "Csf1r",
          "Pdgfra", "Vcan", "Flt1","Vtn", "Igfbp7")
markers.to.plot <- as.matrix(sig_EN)
pdf("male_annotation_combine.pdf", width=10, height=5)
DotPlot(object = male_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()
dev.off()
