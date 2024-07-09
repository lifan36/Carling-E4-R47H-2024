#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)
setwd("/athena/ganlab/scratch/lif4001/LG31E3E4/integration_4more")
LG31E3E4_integrated <- readRDS("LG31E3E4_integrated_PCA_0.1.rds")

Idents(LG31E3E4_integrated) <- "seurat_clusters"
LG31E3E4_integrated <- RenameIdents(LG31E3E4_integrated,
                                    `0` = "oligodendrocytes", `1`="astrocytes", `2`="excitatory neurons", `3`="excitatory neurons",
                                    `4`="microglia", `5`="excitatory neurons", `6`="excitatory neurons", `7`="excitatory neurons",
                                    `8`="OPCs", `9`="inhibitory neurons", `10`="inhibitory neurons", `11`="excitatory neurons",
                                    `12`="vascular cells", `13`="excitatory neurons", `14`="CHOR", `15`="vascular cells"
)

# Figure S2A
pdf("LG31E3E4_integrated_umap_annotation.pdf", width=6, height=3.8)
DimPlot(LG31E3E4_integrated, reduction = 'umap', label = T)
dev.off()

LG31E3E4_integrated$celltype.orig.ident <- paste(Idents(LG31E3E4_integrated), LG31E3E4_integrated$orig.ident, sep = "_")
LG31E3E4_integrated$celltype <- Idents(LG31E3E4_integrated)

saveRDS(LG31E3E4_integrated, file = "LG31E3E4_integrated_Annotation.rds")


Idents(LG31E3E4_integrated) <- "celltype"
# Figure S2B
pdf("annotation_1.pdf", width=10.5, height=3.2)
DotPlot(LG31E3E4_integrated, features = c("Plp1", "Mbp", "Mobp", "Clu", "Plpp3",
                                          "Pla2g7","Slc17a7", "Nrgn","Cx3cr1", "P2ry12", "Csf1r","Vcan", "Pdgfra", "Gad1", "Gad2", "Bmp6", "Adam12", 
                                          "Cped1","Clic6","Ttr")) + RotatedAxis() + scale_colour_gradient2(low = "darkblue", mid = "white", high = "darkred")
dev.off()


# Figure S2H,I
Idents(LG31E3E4_integrated) <- "Sample_Name"
pdf("correlation_1_temp.pdf", width=12, height=3.5)
plot1 <- FeatureScatter(LG31E3E4_integrated, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(LG31E3E4_integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()

Idents(LG31E3E4_integrated) <- "Sample_Name"
pdf("correlation_2_temp.pdf", width=12, height=8)
plot1 <- FeatureScatter(LG31E3E4_integrated, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(LG31E3E4_integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()

# Figure S2C
# calculate ratio of each genotype in each cell type cluster
a<-as.data.frame(table(LG31E3E4_integrated$Condition,LG31E3E4_integrated$celltype))
colnames(a)<-c("clusters","cell.type","cell.no")
agg<-aggregate(cell.no~clusters,a,sum)
a$cluster.total <- agg$cell.no[match(a$clusters,agg$clusters)]
a$ratio<-a$cell.no/a$cluster.total

write.csv(a, "LG31E3E4_by_Condition_ratio.csv")

ggplot(a,aes(x=clusters, y=ratio, fill=cell.type))+
  geom_bar(stat="identity")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("Genotype")+
  ylab("Cell type ratio per genotype") + RotatedAxis()

ggsave("genotype_celltype_distribution_new.pdf",plot=last_plot(),
       width=4,height=3,units="in")

# Figure S2D
# calculate ratio of each sample in each cell type cluster
a<-as.data.frame(table(LG31E3E4_integrated$Sample_Name,LG31E3E4_integrated$celltype))
colnames(a)<-c("clusters","cell.type","cell.no")
agg<-aggregate(cell.no~clusters,a,sum)
a$cluster.total <- agg$cell.no[match(a$clusters,agg$clusters)]
a$ratio<-a$cell.no/a$cluster.total

write.csv(a, "LG31E3E4_by_Sample_ratio.csv")

ggplot(a,aes(x=clusters, y=ratio, fill=cell.type))+
  geom_bar(stat="identity")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("Sample")+
  ylab("Cell type ratio per sample") + RotatedAxis()

ggsave("sample_celltype_distribution_new.pdf",plot=last_plot(),
       width=6,height=3,units="in")


Idents(LG31E3E4_integrated) <- "celltype"
DefaultAssay(LG31E3E4_integrated) <- 'RNA'
pdf("LG31E3E4_integrated_umap_annotation_noLabel.pdf", width=6, height=4)
DimPlot(LG31E3E4_integrated, reduction = 'umap', label = F)
dev.off()

Cluster_EN <- subset(LG31E3E4_integrated, idents = "excitatory neurons")
Cluster_IN <- subset(LG31E3E4_integrated, idents = "inhibitory neurons")
Cluster_MG <- subset(LG31E3E4_integrated, idents = "microglia")
Cluster_AST <- subset(LG31E3E4_integrated, idents = "astrocytes")
Cluster_OL <- subset(LG31E3E4_integrated, idents = "oligodendrocytes")
Cluster_OPC <- subset(LG31E3E4_integrated, idents = "OPCs")
Cluster_VC <- subset(LG31E3E4_integrated, idents = "vascular cells")
Cluster_CHOR <- subset(LG31E3E4_integrated, idents = "CHOR")


saveRDS(Cluster_EN, file = "LG31E3E4_EN_subset.rds")
saveRDS(Cluster_IN, file = "LG31E3E4_IN_subset.rds")
saveRDS(Cluster_MG, file = "LG31E3E4_MG_subset.rds")
saveRDS(Cluster_AST, file = "LG31E3E4_AST_subset.rds")
saveRDS(Cluster_OL, file = "LG31E3E4_OL_subset.rds")
saveRDS(Cluster_OPC, file = "LG31E3E4_OPC_subset.rds")
saveRDS(Cluster_VC, file = "LG31E3E4_VC_subset.rds")
saveRDS(Cluster_CHOR, file = "LG31E3E4_CHOR_subset.rds")
