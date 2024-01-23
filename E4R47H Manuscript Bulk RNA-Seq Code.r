# Bulk RNA-Seq analysis of primary microglia with humanized APOE and TREM2, at baseline (unstimulated) or with 0N4R tau fibril stimulation
# Genotypes: E4/E4-CV/CV, E4/E4-R47H/CV, E3/E3-CV/CV, E3/E3-R47H/CV

# Packages used
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(VennDiagram)
library(ggrepel)
library(data.table)
library(reshape2)
library(circlize)
library(stringr)
library(dplyr)
library(patchwork)

# Directory
setwd("~/Documents/Weill Cornell/Gan Lab/Sequencing/2023 R47HCV-APOE Bulk-Seq")

# DESeq2 and DE analysis ####
# Remove duplicate genes from readcounts
data <- read.csv("readcounts.csv", header=T)
data <- subset(data, !duplicated(data$gene_name))
row.names(data) <- data$gene_name
data <- data[,-1]
write.csv(data, "readcounts_noduplicates.csv")

#DESeq2 analysis
data <- read.csv("readcounts_noduplicates.csv", header=T, row.names=1)
filtered = data[rowSums(data)>15,]
sample <- colnames(data)
genotype <- c(rep("E3_CV",6),rep("E3_R47H",6),rep("E4_CV",6),rep("E4_R47H",6))
treatment <- c(rep("Baseline",3),rep("Tau",3),rep("Baseline",3),rep("Tau",3),rep("Baseline",3),rep("Tau",3),rep("Baseline",3),rep("Tau",3))
genotype_treatment <- paste(genotype,treatment,sep="_")
meta <- data.frame(sample=sample, genotype=genotype, treatment=treatment, genotype_treatment=genotype_treatment)

dds <- DESeqDataSetFromMatrix(countData = filtered, colData = meta, design = ~genotype_treatment)
rld <- rlog(dds, blind = T)

plotPCA(rld, intgroup = "genotype_treatment", ntop = 500)+
  geom_text_repel(aes(label = genotype_treatment))+
  theme_classic()

rld_matrix <- assay(rld)
rld_cor <- cor(rld_matrix)
pheatmap(rld_cor)

dds <- DESeq(dds)

# Pairwise Comparisons

# APOE4 R47H vs. CV (Baseline)
contrast_R47HvsCV_E4 <- c("genotype_treatment","E4_R47H_Baseline","E4_CV_Baseline")
res_R47HvsCV_E4_unshrunken <- results(dds,contrast=contrast_R47HvsCV_E4,alpha=0.05)
res_R47HvsCV_E4 <- lfcShrink(dds,contrast=contrast_R47HvsCV_E4,res=res_R47HvsCV_E4_unshrunken, type="normal")
write.csv(res_R47HvsCV_E4, "DE_R47HvsCV_E4.csv")

# APOE4 R47H vs. CV (Tau)
contrast_R47HvsCV_E4_tau <- c("genotype_treatment","E4_R47H_Tau","E4_CV_Tau")
res_R47HvsCV_E4_tau_unshrunken <- results(dds,contrast=contrast_R47HvsCV_E4_tau,alpha=0.05)
res_R47HvsCV_E4_tau <- lfcShrink(dds,contrast=contrast_R47HvsCV_E4_tau,res=res_R47HvsCV_E4_tau_unshrunken, type="normal")
write.csv(res_R47HvsCV_E4_tau, "DE_R47HvsCV_E4_Tau.csv")

# Figure 5B-C - Gene set enrichment analysis (hallmark pathways) ####

# Baseline (unstimulated)
GO <- read.csv("Hallmark_E4_R47HvsCV_Baseline.csv",header=T)
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

# 0N4R tau fibril stimulation
GO <- read.csv("Hallmark_E4_R47HvsCV_Tau.csv",header=T)
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

# Figure 5D - Volcano plot ####
options(ggrepel.max.overlaps = 20)

data <- read.csv("DE_R47HvsCV_E4_Tau.csv",header=T,row.names=1)
data$color[data$log2FoldChange > 0] <- "#D7301F"
data$color[data$log2FoldChange < 0] <- "#2B8CBE"
data$color[data$padj > 0.05] <- "grey"
labels <- subset(data, -log(padj) > 25 | log2FoldChange < -3 | log2FoldChange > 3)

ggplot(data, aes(x = log2FoldChange, y = -log(padj), color = color))+
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
  scale_x_continuous(limits=c(-8, 8)) +
  xlab("log2FoldChange")+
  ylab("-log(padj)")+
  ggtitle("E4-R47H-Tau vs. E4-CV-Tau")+
  scale_color_manual(values = c("#D7301F" = "#D7301F", "#2B8CBE" = "#2B8CBE","grey"="grey"),
                     name = "DEG",
                     breaks = c("#2B8CBE","#D7301F","grey"),
                     labels = c("Downregulated","Upregulated","No Change"))+
  geom_text_repel(labels, mapping = aes(label = row.names(labels)))

# Figure 5E - IPA upstream regulator dot plot ####
ipa <- read.csv('USRs_R47HvsCV_E4_Tau.csv')
ggplot(data=ipa, aes(x=reorder(Upstream.Regulator,Activation.z.score), y= (Activation.z.score), fill=-log(p.value.of.overlap))) +
  geom_dotplot(position ="identity",binaxis='y', stackdir='center')+ coord_flip()+ scale_fill_continuous(low = 'blue', high='red')+theme(legend.position="top")+theme(axis.line = element_line(colour = "black",size = 0.5, linetype = "solid"),
                                                                                                                                                                      panel.background = element_rect(fill = "transparent",colour = NA),
                                                                                                                                                                      plot.background = element_rect(fill = "transparent",colour = NA))+labs(x=' ', y ="Activation Z-Score")

# Supplementary Figure S5 - quality control assessment of bulk RNA-Seq for APOE4 only ####

# Figure S5A - Venn diagram
Baseline <- read.csv("DE_R47HvsCV_E4.csv",header=T)
Baseline = filter(Baseline, padj < 0.05&(log2FoldChange>=0.1|log2FoldChange<=-0.1))
Tau <- read.csv("DE_R47HvsCV_E4_Tau.csv",header=T)
Tau = filter(Tau, padj < 0.05&(log2FoldChange>=0.1|log2FoldChange<=-0.1))

set1 <- Baseline$X
set2 <- Tau$X
myCol <- brewer.pal(3, "Pastel2")

grid.newpage()
venn<-venn.diagram(
  x = list(set1,set2),
  category.names = c("Baseline","Tau"),
  filename = NULL,
  output=FALSE,
  show.plot=TRUE,
  
  height = 3600, 
  width = 1800 , 
  resolution = 300,
  #compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill =c('blue','red'),
  
  # Numbers
  cex = 1.5,
  # Set names
  cat.cex = 1,
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
)
grid.newpage()
grid.draw(venn)
library(grDevices)
pdf(file="venn.pdf")
grid.draw(venn)
dev.off()

# Figure S5B-C QC
data <- read.csv("readcounts_noduplicates_APOE4.csv", header=T, row.names=1)
filtered = data[rowSums(data)>15,]
sample <- colnames(data)
genotype <- c(rep("E4_CV",6),rep("E4_R47H",6))
treatment <- c(rep("Baseline",3),rep("Tau",3),rep("Baseline",3),rep("Tau",3))
genotype_treatment <- paste(genotype,treatment,sep="_")
meta <- data.frame(sample=sample, genotype=genotype, treatment=treatment, genotype_treatment=genotype_treatment)

all(colnames(filtered) %in% meta$sample)

dds <- DESeqDataSetFromMatrix(countData = filtered, colData = meta, design = ~genotype_treatment)
rld <- rlog(dds, blind = T)

# Figure S5B - PCA plot
plotPCA(rld, intgroup = "genotype_treatment", ntop = 500)+
  geom_text_repel(aes(label = genotype_treatment))+
  theme_classic()

# Figure S5C - sample-sample correlation heatmap
rld_matrix <- assay(rld)
rld_cor <- cor(rld_matrix)
pheatmap(rld_cor)

# Figure 6F, 6H - senescence gene expression heatmaps ####
counts <- counts(dds,normalized = T)
write.csv(counts, "counts_normalized.csv")

mypal = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256))
counts <- read.csv("counts_normalized.csv", header=T, row.names=1)

# Figure 6F - SASP gene expression heatmap
senescence <- read.csv("SASP_Genes.csv", header=T)
row.names(senescence) <- senescence$SASP
overlap.senescence <- subset(counts, row.names(counts) %in% row.names(senescence))

p=pheatmap::pheatmap(overlap.senescence[,c(13:24)],
                     scale="row", clustering_method="ward.D", color = mypal,border_color = NA,
                     angle_col = 90, fontsize = 11, main='SASP Genes', legend_labels = 'Average expression',cluster_rows=FALSE)


# Figure 6H - proliferation gene expression heatmap
senescence <- read.csv("Proliferation_Genes.csv", header=T)
row.names(proliferation) <- proliferation$Proliferation
overlap.proliferation <- subset(counts, row.names(counts) %in% row.names(proliferation))

p=pheatmap::pheatmap(overlap.proliferation[,c(13:24)],
                     scale="row", clustering_method="ward.D", color = mypal,border_color = NA,
                     angle_col = 90, fontsize = 11, main='Proliferation Genes', legend_labels = 'Average expression',cluster_rows=FALSE)

