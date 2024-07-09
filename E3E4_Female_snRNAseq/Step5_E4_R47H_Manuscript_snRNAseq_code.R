# Single-Nuclei RNA-Seq analysis of mouse hippocampal tissue (all female mice, n=3 per genotype)
# Genotypes: E4/E4_mTrem2+/+_NTG, E4/E4_mTrem2+/+_P301S, E4/E4_R47H/+_P301S, E3/E3_mTrem2+/+_NTG, E3/E3_mTrem2+/+_P301S, E3/E3_R47H/+_P301S

# Packages used
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(pheatmap)
library(stringr)
library(Seurat)
library(dplyr)
library(cowplot)
library(reshape2)
library(openxlsx)
library(smplot2)
library(devtools)

# Directory
setwd("~/Documents/Weill Cornell/Gan Lab/Sequencing/2023 R47H-APOE-Tau Reintegrated")

# Figure 2G - IPA upstream regulator dot plot for MG4 cluster ####
ipa <- read.csv('USRs_MG4.csv')
ggplot(data=ipa, aes(x=reorder(name,z.score), y= (z.score), fill=-log(p.value))) +
  geom_dotplot(position ="identity",binaxis='y', stackdir='center')+ coord_flip()+ scale_fill_continuous(low = 'blue', high='red')+theme(legend.position="top")+theme(axis.line = element_line(colour = "black",size = 0.5, linetype = "solid"),
                                                                                                                                                                      panel.background = element_rect(fill = "transparent",colour = NA),
                                                                                                                                                                      plot.background = element_rect(fill = "transparent",colour = NA))+labs(x=' ', y ="Activation Z-Score")

# Figure 2I - APOE4 microglia IFN-stimulated gene dot plot ####
MG <- readRDS("LG31E3E4_E4_MG_03132023.rds")
DefaultAssay(MG) <- "RNA"
Idents(MG) <- "Condition"

DotPlot(MG, features = c( "Pik3ap1", "B2m", "Trim14", "Rnf213", "Nlrc5", "Sp100", "Sp110", "Usp25",
                          "Samd9l", "Ifi204", "Gm4951", "Adar"
)) + RotatedAxis() + scale_colour_gradient2(low = "darkblue", mid = "white", high = "darkred")

# Supplementary Figure S3D - APOE3 microglia IFN-stimulated gene dot plot ####
MG <- readRDS("LG31E3E4_E3_MG_03132023.rds")
DefaultAssay(MG) <- "RNA"
Idents(MG) <- "Condition"

DotPlot(MG, features = c( "Pik3ap1", "B2m", "Trim14", "Rnf213", "Nlrc5", "Sp100", "Sp110", "Usp25",
                          "Samd9l", "Ifi204", "Gm4951", "Adar"
)) + RotatedAxis() + scale_colour_gradient2(low = "darkblue", mid = "white", high = "darkred")

# Figure 3A - APOE4 oligodendrocyte volcano plot ####
options(ggrepel.max.overlaps = 30)

data <- read.csv("OL_E4_P301S_R47H_vs_E4_P301S_WT_DEGs.csv",header=T,row.names=1)
data$color[data$avg_log2FC > 0] <- "#D7301F"
data$color[data$avg_log2FC < 0] <- "#2B8CBE"
data$color[data$p_val_adj > 0.05] <- "grey"
labels <- subset(data, -log(p_val_adj) > 100 | avg_log2FC < -0.5 | avg_log2FC > 0.5)

ggplot(data, aes(x = avg_log2FC, y = -log(p_val_adj), color = color))+
  geom_point(shape=19, alpha=1, size=2)+
  theme_bw()+
  theme(panel.grid.major.x  = element_blank(),
        panel.grid.major.y  = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.text.x = element_text(colour = "black", size=15),
        axis.text.y = element_text(colour = "black", size=15),
        axis.title.x = element_text(colour = "black", size=15),
        axis.title.y = element_text(colour = "black", size=15),
        plot.title = element_text(size = 15, face = "bold")) +
  theme(aspect.ratio = 1) +
  scale_x_continuous(limits=c(-1, 1)) +
  xlab("log2FoldChange")+
  ylab("-log(padj)")+
  ggtitle("OL: E4-R47H-Tau vs. E4-CV-Tau")+
  scale_color_manual(values = c("#D7301F" = "#D7301F", "#2B8CBE" = "#2B8CBE","grey"="grey"),
                     name = "DEG",
                     breaks = c("#2B8CBE","#D7301F","grey"),
                     labels = c("Downregulated","Upregulated","No Change"))+
  geom_text_repel(labels, mapping = aes(label = row.names(labels)))

# Figure 3B - APOE4 disease-associated oligodendrocyte gene dot plot ####
OL <- readRDS("E4_OL_reclusted_res0.15.rds")
DefaultAssay(OL) <- "RNA"
Idents(OL) <- "Condition"

DotPlot(OL, features = c("C4b", "H2-D1", "H2-K1", "Plin4", "B2m", "Cdkn1a"
), scale.max = 30, scale.min = 0) + RotatedAxis() + scale_colour_gradient2(low = "darkblue", mid = "white", high = "darkred")

# Figure 3C - APOE4 astrocyte volcano plot ####
options(ggrepel.max.overlaps = 25)

data <- read.csv("AST_E4_P301S_R47H_vs_E4_P301S_WT_DEGs.csv",header=T,row.names=1)
data$color[data$avg_log2FC > 0] <- "#D7301F"
data$color[data$avg_log2FC < 0] <- "#2B8CBE"
data$color[data$p_val_adj > 0.05] <- "grey"
labels <- subset(data, -log(p_val_adj) > 50 | avg_log2FC < -3 | avg_log2FC > 3)

ggplot(data, aes(x = avg_log2FC, y = -log(p_val_adj), color = color))+
  geom_point(shape=19, alpha=1, size=2)+
  theme_bw()+
  theme(panel.grid.major.x  = element_blank(),
        panel.grid.major.y  = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.text.x = element_text(colour = "black", size=15),
        axis.text.y = element_text(colour = "black", size=15),
        axis.title.x = element_text(colour = "black", size=15),
        axis.title.y = element_text(colour = "black", size=15),
        plot.title = element_text(size = 15, face = "bold")) +
  theme(aspect.ratio = 1) +
  scale_x_continuous(limits=c(-1, 1)) +
  xlab("log2FoldChange")+
  ylab("-log(padj)")+
  ggtitle("AST: E4-R47H-Tau vs. E4-CV-Tau")+
  scale_color_manual(values = c("#D7301F" = "#D7301F", "#2B8CBE" = "#2B8CBE","grey"="grey"),
                     name = "DEG",
                     breaks = c("#2B8CBE","#D7301F","grey"),
                     labels = c("Downregulated","Upregulated","No Change"))+
  geom_text_repel(labels, mapping = aes(label = row.names(labels)))

# Figure 3D - APOE4 reactive astrocyte gene dot plot ####
AST <- readRDS("E4_AST_reclusted_res0.15.rds")
DefaultAssay(AST) <- 'RNA'
Idents(AST) <- "Condition"

DotPlot(AST, features = c("Gfap", "C4b", "Stat3", "H2-D1", "Ggta1", "Fkbp5"
)) + RotatedAxis() + scale_colour_gradient2(low = "darkblue", mid = "white", high = "darkred")

# Figure 3E - APOE4 oligodendrocyte gene set enrichment analysis (hallmark pathways) ####
GO <- read.csv("Hallmark_OL_Combined.csv",header=T)
colnames(GO) <- c("name", "genes_in_gsea", "description", "genes_in_data", "k/K", "p-value", "FDR", "Type")
GO$FDR <- as.numeric(GO$FDR)
GO$logP <- -log10(GO$FDR)
GO$Name <-  gsub("HALLMARK_", "", GO$name)
GO$Name <- gsub("*_", " ", GO$Name)
GO$Name <-  factor(GO$Name, levels=rev(GO$Name))

ggplot(data=GO, aes(x=Name  , y=logP, fill=Type)) +
  theme_classic() +
  ylab("-Log(FDR)") + xlab(NULL) +
  geom_bar(stat="Identity",  width=0.7, alpha=0.8) +
  coord_flip() + 
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue"))+
  theme(aspect.ratio = 1.5)

# Figure 3F - APOE4 oligodendrocyte interferon gamma hallmark pathway heatmap ####
OL <- readRDS("E4_OL_reclusted_res0.15.rds")
DefaultAssay(OL) <- "RNA"
Idents(OL) <- "Condition"

Genes <- c("Jak2", "B2m", "Arid5b", "Nfkb1", "Cdkn1a", "Nfkbia",
           "Stat3", "Herc6", "Znfx1")

cluster.averages <- AverageExpression(OL, return.seurat = TRUE)
cluster.averages <- cluster.averages$RNA
cluster.averages <- cluster.averages[rownames(cluster.averages) %in% Genes]
cluster.averages <- as.matrix(cluster.averages)

mypal = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256))

p=pheatmap::pheatmap(t(cluster.averages),
                     scale="column", clustering_method="ward.D", color = mypal,border_color = NA,
                     angle_col = 90, fontsize = 11, main='OL: Interferon Gamma Response', legend_labels = 'Average expression',cluster_rows=FALSE)

# Figure 3G - APOE4 astrocyte gene set enrichment analysis (hallmark pathways) ####
GO <- read.csv("Hallmark_AST_Combined.csv",header=T)
colnames(GO) <- c("name", "genes_in_gsea", "description", "genes_in_data", "k/K", "p-value", "FDR", "Type")
GO$FDR <- as.numeric(GO$FDR)
GO$logP <- -log10(GO$FDR)
GO$Name <-  gsub("HALLMARK_", "", GO$name)
GO$Name <- gsub("*_", " ", GO$Name)
GO$Name <-  factor(GO$Name, levels=rev(GO$Name))

ggplot(data=GO, aes(x=Name  , y=logP, fill=Type)) +
  theme_classic() +
  ylab("-Log(FDR)") + xlab(NULL) +
  geom_bar(stat="Identity",  width=0.7, alpha=0.8) +
  coord_flip() + 
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue"))+
  theme(aspect.ratio = 1.5)

# Figure 3H - APOE4 astrocyte interferon gamma and alpha hallmark pathways heatmap ####
AST <- readRDS("E4_AST_reclusted_res0.15.rds")
DefaultAssay(AST) <- 'RNA'
Idents(AST) <- "Condition"

Genes <- c('Ifih1','Rtp4','Irf2','Herc6','Ddx60',
           'Stat3','Nlrc5','Trim5')

cluster.averages <- AverageExpression(AST, return.seurat = TRUE)
cluster.averages <- cluster.averages$RNA
cluster.averages <- cluster.averages[rownames(cluster.averages) %in% Genes]
cluster.averages <- as.matrix(cluster.averages)

mypal = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256))

p=pheatmap::pheatmap(t(cluster.averages),
                     scale="column", clustering_method="ward.D", color = mypal,border_color = NA,
                     angle_col = 90, fontsize = 11, main='AST: Interferon Response', legend_labels = 'Average expression',cluster_rows=FALSE)

# Figure 4A - APOE4 microglia cGAS dot plot ####
MG <- readRDS("LG31E3E4_E4_MG_03132023.rds")
DefaultAssay(MG) <- "RNA"
Idents(MG) <- "Condition"

DotPlot(MG, features = c("Cgas"
), scale.max = 12, scale.min = 0) + RotatedAxis() + scale_colour_gradient2(low = "darkblue", mid = "white", high = "darkred")

# Figure 6A - APOE4 microglia senescence gene heatmap ####
MG <- readRDS("LG31E3E4_E4_MG_03132023.rds")
DefaultAssay(MG) <- 'RNA'
Idents(MG) <- "Condition"

Senescence_genes <- c('Akt3','Anapc1','Anapc15', 'Atm','Ep300',
                      'Hbp1','Ikbkb','Itsn2','Mapkapk3','Mtor','Pik3ca',
                      'Pik3cb','Pik3cd','Raf1','Rap1b','Rb1', 'B2m', 'Tnrc6a')

cluster.averages <- AverageExpression(MG, return.seurat = TRUE)
cluster.averages <- cluster.averages$RNA
cluster.averages <- cluster.averages[rownames(cluster.averages) %in% Senescence_genes]
cluster.averages <- as.matrix(cluster.averages)

mypal = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256))

p=pheatmap::pheatmap(t(cluster.averages),
                     scale="column", clustering_method="ward.D", color = mypal,border_color = NA,
                     angle_col = 90, fontsize = 11, main='Senescence-associated genes', legend_labels = 'Average expression',cluster_rows=FALSE)


