
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)
setwd("/athena/ganlab/scratch/lif4001/LG31E3E4_male/integration_E4")
male_integrated <- readRDS("male_integrated_PCA_0.1.rds")
#Remove cluster 15 - Choroid plexus epithelial cells, 13 - mixed neuron
male_integrated <- subset(male_integrated, idents=c("13","15"), invert=T)
Idents(male_integrated) <- "seurat_clusters"
male_integrated <- RenameIdents(male_integrated,
                                 `0` = "oligodendrocytes", `1`="astrocytes", `2`="excitatory neurons", `3`="microglia",
                                 `4`="excitatory neurons", `5`="excitatory neurons", `6`="excitatory neurons", `7`="excitatory neurons",
                                 `8`="inhibitory neurons",`9`="vascular cells", `10`="inhibitory neurons", `11`="excitatory neurons", `12`="OPCs",
                                `14`="excitatory neurons", `16`="vascular cells"
)

# Figure S5A
pdf("male_integrated_umap_annotation.pdf", width=6, height=3.8)
DimPlot(male_integrated, reduction = 'umap', label = T)
dev.off()

male_integrated$celltype.orig.ident <- paste(Idents(male_integrated), male_integrated$orig.ident, sep = "_")
male_integrated$celltype <- Idents(male_integrated)

saveRDS(male_integrated, file = "male_integrated_Annotation.rds")

# Figure S5B
#markers for annotation
pdf("annotation_1.pdf", width=10.5, height=2.7)
DotPlot(data, features = c("Plp1", "Mbp", "Mobp","Slc17a7", "Nrgn", "Clu", "Plpp3",
                           "Pla2g7", "Cx3cr1", "P2ry12", "Csf1r","Gad1", "Gad2","Vcan", "Pdgfra", "Bmp6", "Adam12",
                           "Cped1","Clic6")) + RotatedAxis()
dev.off()


# Figure S5C
data <- male_integrated
# calculate ratio of each genotype in each cell type cluster
a<-as.data.frame(table(data$Condition,data$celltype))
colnames(a)<-c("clusters","cell.type","cell.no")
agg<-aggregate(cell.no~clusters,a,sum)
a$cluster.total <- agg$cell.no[match(a$clusters,agg$clusters)]
a$ratio<-a$cell.no/a$cluster.total

ggplot(a,aes(x=clusters, y=ratio, fill=cell.type))+
  geom_bar(stat="identity")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("Genotype")+
  ylab("Cell type ratio per genotype") + RotatedAxis()

ggsave("genotype_celltype_distribution.pdf",plot=last_plot(),
       width=4,height=4,units="in")


# Figure S5D
data <- male_integrated
# calculate ratio of each sample in each cell type cluster
a<-as.data.frame(table(data$Sample_Name,data$celltype))
colnames(a)<-c("clusters","cell.type","cell.no")
agg<-aggregate(cell.no~clusters,a,sum)
a$cluster.total <- agg$cell.no[match(a$clusters,agg$clusters)]
a$ratio<-a$cell.no/a$cluster.total

ggplot(a,aes(x=clusters, y=ratio, fill=cell.type))+
  geom_bar(stat="identity")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("Sample")+
  ylab("Cell type ratio per sample") + RotatedAxis()

ggsave("sample_celltype_distribution.pdf",plot=last_plot(),
       width=6,height=4,units="in")


# Figure S5H, I
Idents(male_integrated) <- "Sample_Name"
pdf("correlation_1_temp.pdf", width=12, height=3.5)
plot1 <- FeatureScatter(male_integrated, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(male_integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()

Cluster_MG <- subset(male_integrated, idents = "microglia")

saveRDS(Cluster_MG, file = "male_MG_subset.rds")


