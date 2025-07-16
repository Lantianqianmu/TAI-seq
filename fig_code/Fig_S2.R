library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(ggnewscale)
library(ggforce)
source("/home/zeemeeuw/YangLab/JoINT-seq/QC/code/utils.R")
setwd("/home/zeemeeuw/YangLab/JoINT-seq/Manuscript/fig_data/fig_S2_data")

# S2A
{
  # 3T3
  data <- read.delim("source_fig_S2A_NIH3T3.txt", header = TRUE)
  data$temp <- factor(data$temp, levels = c("37","50"))
  one <- data %>% dplyr::filter(temp == 37)
  two <- data %>% dplyr::filter(temp == 50)
  ggplot(data = data, aes(x = temp, y = log10(counts))) +
    theme_classic() +
    geom_violin(aes(fill = temp), adjust = 1, trim = FALSE, width = 0.8, size = 0.5, scale = "width", show.legend = TRUE) +
    geom_boxplot(aes(group = temp), width = 0.1, alpha = 0.7, show.legend = FALSE, outlier.shape = NA, notch = FALSE, fill = "lightGrey", size = 0.5) +
    scale_fill_manual(values = c('#fd5f00', '#900c3f')) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 13, color = c("black", "black")),
          axis.text.y = element_text(size = 13, color = "black"),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15)) +
    labs(title = "3T3 ATAC 37 vs 50") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylim(c(3, 5.5)) +
    ylab(expression(log[10](fragments))) +
    annotate(geom = "text", x = 1, y = 5.4, angle = 30, label = round(median(one$counts)), size = 3.5) +
    annotate(geom = "text", x = 2, y = 5.4, angle = 30, label = round(median(two$counts)), size = 3.5) +
    coord_fixed(1.2)
  
  
  # GM12878
  data <- read.delim("source_fig_S2A_GM12878.txt", header = TRUE)
  data$temp <- factor(data$temp, levels = c("37","50"))
  one <- data %>% dplyr::filter(temp == 37)
  two <- data %>% dplyr::filter(temp == 50)
  ggplot(data = data, aes(x = temp, y = log10(counts))) +
    theme_classic() +
    geom_violin(aes(fill = temp), adjust = 1, trim = FALSE, width = 0.8, size = 0.5, scale = "width", show.legend = TRUE) +
    geom_boxplot(aes(group = temp), width = 0.1, alpha = 0.7, show.legend = FALSE, outlier.shape = NA, notch = FALSE, fill = "lightGrey", size = 0.5) +
    scale_fill_manual(values = c('#fd5f00', '#900c3f')) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 13, color = c("black", "black")),
          axis.text.y = element_text(size = 13, color = "black"),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15)) +
    labs(title = "3T3 ATAC 37 vs 50") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylim(c(3, 5.5)) +
    ylab(expression(log[10](fragments))) +
    annotate(geom = "text", x = 1, y = 5.4, angle = 30, label = round(median(one$counts)), size = 3.5) +
    annotate(geom = "text", x = 2, y = 5.4, angle = 30, label = round(median(two$counts)), size = 3.5) +
    coord_fixed(1.2)
}


# S2B
{
  # 3T3
  data <- read.delim("source_fig_S2B_NIH3T3.txt", header = TRUE)
  data$temp <- factor(data$temp, levels = c("37","50"))
  one <- data %>% dplyr::filter(temp == 37)
  two <- data %>% dplyr::filter(temp == 50)
  ggplot(data = data, aes(x = temp, y = log10(counts))) +
    theme_classic() +
    geom_violin(aes(fill = temp), adjust = 1, trim = FALSE, width = 0.8, size = 0.5, scale = "width", show.legend = TRUE) +
    geom_boxplot(aes(group = temp), width = 0.1, alpha = 0.7, show.legend = FALSE, outlier.shape = NA, notch = FALSE, fill = "lightGrey", size = 0.5) +
    scale_fill_manual(values = c('#fd5f00', '#900c3f')) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 13, color = c("black", "black")),
          axis.text.y = element_text(size = 13, color = "black"),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15)) +
    labs(title = "3T3 umis 37 vs 50") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylim(c(3, 5.5)) +
    ylab(expression(log[10](umis))) +
    annotate(geom = "text", x = 1, y = 5.4, angle = 30, label = round(median(one$counts)), size = 3.5) +
    annotate(geom = "text", x = 2, y = 5.4, angle = 30, label = round(median(two$counts)), size = 3.5) +
    coord_fixed(1.2)
  
  # GM12878
  data <- read.delim("source_fig_S2B_GM12878.txt", header = TRUE)
  data$temp <- factor(data$temp, levels = c("37","50"))
  one <- data %>% dplyr::filter(temp == 37)
  two <- data %>% dplyr::filter(temp == 50)
  
  ggplot(data = data, aes(x = temp, y = log10(counts))) +
    theme_classic() +
    geom_violin(aes(fill = temp), adjust = 1, trim = FALSE, width = 0.8, size = 0.5, scale = "width", show.legend = TRUE) +
    geom_boxplot(aes(group = temp), width = 0.1, alpha = 0.7, show.legend = FALSE, outlier.shape = NA, notch = FALSE, fill = "lightGrey", size = 0.5) +
    scale_fill_manual(values = c('#0278ae', '#8bd156')) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 13, color = c("black", "black")),
          axis.text.y = element_text(size = 13, color = "black"),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15)) +
    labs(title = "GM12878 umis 37 vs 50") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylim(c(3, 5.5)) +
    ylab(expression(log[10](umis))) +
    annotate(geom = "text", x = 1, y = 5.4, angle = 30, label = round(median(one$counts)), size = 3.5) +
    annotate(geom = "text", x = 2, y = 5.4, angle = 30, label = round(median(two$counts)), size = 3.5) +
    coord_fixed(1.2)
}

# S2C
{
  # 3T3
  data <- read.delim("source_fig_S2C_NIH3T3.txt", header = TRUE)
  data$temp <- factor(data$temp, levels = c("37","50"))
  one <- data %>% dplyr::filter(temp == 37)
  two <- data %>% dplyr::filter(temp == 50)
  
  ggplot(data = data, aes(x = temp, y = log10(counts))) +
    theme_classic() +
    geom_violin(aes(fill = temp), adjust = 1, trim = FALSE, width = 0.8, size = 0.5, scale = "width", show.legend = TRUE) +
    geom_boxplot(aes(group = temp), width = 0.1, alpha = 0.7, show.legend = FALSE, outlier.shape = NA, notch = FALSE, fill = "lightGrey", size = 0.5) +
    scale_fill_manual(values = c('#fd5f00', '#900c3f')) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 13, color = c("black", "black")),
          axis.text.y = element_text(size = 13, color = "black"),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15)) +
    labs(title = "3T3 genes 37 vs 50") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylim(c(2, 4.5)) +
    ylab(expression(log[10](genes))) +
    annotate(geom = "text", x = 1, y = 4.5, angle = 30, label = round(median(one$counts)), size = 3.5) +
    annotate(geom = "text", x = 2, y = 4.5, angle = 30, label = round(median(two$counts)), size = 3.5) +
    coord_fixed(1.2)
  
  # GM12878
  data <- read.delim("source_fig_S2C_GM12878.txt", header = TRUE)
  data$temp <- factor(data$temp, levels = c("37","50"))
  one <- data %>% dplyr::filter(temp == 37)
  two <- data %>% dplyr::filter(temp == 50)
  
  ggplot(data = data, aes(x = temp, y = log10(counts))) +
    theme_classic() +
    geom_violin(aes(fill = temp), adjust = 1, trim = FALSE, width = 0.8, size = 0.5, scale = "width", show.legend = TRUE) +
    geom_boxplot(aes(group = temp), width = 0.1, alpha = 0.7, show.legend = FALSE, outlier.shape = NA, notch = FALSE, fill = "lightGrey", size = 0.5) +
    scale_fill_manual(values = c('#0278ae', '#8bd156')) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 13, color = c("black", "black")),
          axis.text.y = element_text(size = 13, color = "black"),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15)) +
    labs(title = "GM12878 genes 37 vs 50") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylim(c(2, 4.5)) +
    ylab(expression(log[10](genes))) +
    annotate(geom = "text", x = 1, y = 4.5, angle = 30, label = round(median(one$counts)), size = 3.5) +
    annotate(geom = "text", x = 2, y = 4.5, angle = 30, label = round(median(two$counts)), size = 3.5) +
    coord_fixed(1.2)
  
}

# S2D
{
  library(VennDiagram)
  Chevalier1 <- c("#355243", "#fbca50", "#c9d5d4", "#baa28a")
  FantasticFox1 <- c("#d37a20", "#dbcb09", "#3a9cbc", "#dd7208", "#a30019")
  Moonrise3 <- c("#75cbdc", "#f0a4af", "#8a863a", "#c2b479", "#f8d068")
  Cavalcanti1 <- c("#ceab0d", "#083215", "#919562", "#6f997a", "#831e11")
  Darjeeling2 <- c("#e6c09e", "#0d5888", "#cb8b3e", "#9cd6d6", "#000000")
  Darjeeling1 <- c("#fb0007", "#139177", "#ed9e08", "#f56f08", "#4caecc")
  Royal2 <- c("#e4c9b2", "#f1c2a5", "#f49d98", "#fcd68f", "#629076")
  IsleofDogs2 <- c("#e4c9b2", "#998273", "#a6723d", "#2b2523", "#151213")
  
  # 37 temp
  atac_1 <- read.csv("G3-ata_hg38_1_bc_whitelist.txt", header = TRUE, row.names = 1)
  atac_2 <- read.csv("G3-ata_mm10_1_bc_whitelist.txt", header = TRUE, row.names = 1)
  
  rna_1 <- read.csv("G3-rna_hg38_1_bc_whitelist.txt", header = TRUE, row.names = 1)$x
  rna_2 <- read.csv("G3-rna_mm10_1_bc_whitelist.txt", header = TRUE, row.names = 1)$x
  
  rna_1 <- rna_1[grepl("\\-GCA$", rna_1)]
  rna_1 <- substr(rna_1, 1, 23)
  rna_1 <- unique(rna_1)
  
  rna_2 <- rna_2[grepl("\\-GCA$", rna_2)]
  rna_2 <- substr(rna_2, 1, 23)
  rna_2 <- unique(rna_2)
  
  
  veendata <- list("atac" = atac_1$x,
                   "rna" = rna_1)
  venn.diagram(veendata, filename = "G3-multione_GM12878_37.tiff", imagetype = "tiff",
               lwd = 0.5, cex = 0.5, cat.cex = 0, width = 2000, height = 2000,
               col="white", fill = Darjeeling1[c(5, 1)],
               disable.logging = TRUE, alpha = 0.5,
               fontfamily = "Arial", cat.fontfamily = "Arial")
  
  veendata <- list("atac" = atac_2$x,
                   "rna" = rna_2)
  venn.diagram(veendata, filename = "G3-multione_3T3_37.tiff", imagetype = "tiff",
               ext.text = FALSE, cex = 0.5, cat.cex = 0, width = 2000, height = 2000,
               col="white", fill = Darjeeling1[c(5, 1)],
               disable.logging = TRUE, alpha = 0.5,
               fontfamily = "Arial", cat.fontfamily = "Arial")
  
  # 50 temp
  atac_1 <- read.csv("G3-ata_hg38_2_bc_whitelist.txt", header = TRUE, row.names = 1)
  atac_2 <- read.csv("/home/zeemeeuw/YangLab/JoINT-seq/QC/data/G3-ata/G3-ata_mm10_2_bc_whitelist.txt", header = TRUE, row.names = 1)
  
  rna_1 <- read.csv("G3-rna_hg38_2_bc_whitelist.txt", header = TRUE, row.names = 1)$x
  rna_2 <- read.csv("G3-rna_mm10_2_bc_whitelist.txt", header = TRUE, row.names = 1)$x
  
  rna_1 <- rna_1[grepl("\\-TCC$", rna_1)]
  rna_1 <- substr(rna_1, 1, 23)
  rna_1 <- unique(rna_1)
  
  rna_2 <- rna_2[grepl("\\-TCC$", rna_2)]
  rna_2 <- substr(rna_2, 1, 23)
  rna_2 <- unique(rna_2)
  
  veendata <- list("atac" = atac_1$x,
                   "rna" = rna_1)
  venn.diagram(veendata, filename = "G3-multione_GM12878_50_55_552_24.tiff", imagetype = "tiff",
               ext.text = FALSE, cex = 0.5, cat.cex = 0, width = 2000, height = 2000,
               col="white", fill = Darjeeling1[c(5, 1)],
               disable.logging = TRUE, alpha = 0.5,
               fontfamily = "Arial", cat.fontfamily = "Arial")
  
  veendata <- list("atac" = atac_2$x,
                   "rna" = rna_2)
  venn.diagram(veendata, filename = "G3-multione_3T3_50_84_1291_31.tiff", imagetype = "tiff",
               ext.text = FALSE, cex = 0.5, cat.cex = 0, width = 2000, height = 2000,
               col="white", fill = Darjeeling1[c(5, 1)],
               disable.logging = TRUE, alpha = 0.5,
               fontfamily = "Arial", cat.fontfamily = "Arial")
}








library(ggplot2)
library(ArchR)
library(Seurat)
library(patchwork)
library(scCustomize)
library(scattermore)
library(scplotter)
library(readr)

RDSdir <- "/home/zeemeeuw/data/GJx/RNA/" # output RDS file dir
ATACdir <- "/home/zeemeeuw/data/GJx/ATAC_overlapped/" # ATAC root dir
seu <- readRDS(paste0(RDSdir, "GJx_overlapped.rds"))
loadArchRProject(ATACdir)

# S2E
{
  GJx_colors <- c("GM12878" = "#3971ff",
                  "Jurkat" = "#ff882e")
  meta <- read_tsv("source_fig_S2E.txt")

  p <- ggplot(meta, aes(x = sample)) +
    geom_violin(aes(fill = sample, y = log10(nFeature_RNA)), alpha = 0.5, trim = TRUE, scale = "width") +
    geom_boxplot(aes(fill = sample, y = log10(nFeature_RNA)), alpha = 1, show.legend = FALSE, width = 0.3, outlier.shape = NA) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.ticks.x = element_blank()) +
    ylim(1, 4) +
    scale_fill_manual(values = GJx_colors) +
    ggtitle("genes")
  print(p)
  p <- ggplot(meta, aes(x = sample)) +
    geom_violin(aes(fill = sample, y = log10(nCount_RNA)), alpha = 0.5, trim = TRUE, scale = "width") +
    geom_boxplot(aes(fill = sample, y = log10(nCount_RNA)), alpha = 1, show.legend = FALSE, width = 0.3, outlier.shape = NA) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.ticks.x = element_blank()) +
    ylim(1, 5) +
    scale_fill_manual(values = GJx_colors) +
    ggtitle("umis")
  print(p)
  p <- ggplot(meta, aes(x = sample)) +
    geom_violin(aes(fill = sample, y = log10(nFrags)), alpha = 0.5, trim = TRUE, scale = "width") +
    geom_boxplot(aes(fill = sample, y = log10(nFrags)), alpha = 1, show.legend = FALSE, width = 0.3, outlier.shape = NA) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.ticks.x = element_blank()) +
    ylim(1, 5.5) +
    scale_fill_manual(values = GJx_colors) +
    ggtitle("nfrags")
  print(p)
}


# S2G
{
  
  CellDimPlot(seu, reduction = "UMAP_LSI_Combined",  theme = "theme_blank", legend.position = "right", group_by = "cell_module", raster = FALSE, pt_size = 0.5, palcolor = GJx_colors)
  
  CellDimPlot(seu, reduction = "UMAP_LSI_RNA",  theme = "theme_blank", legend.position = "right", group_by = "cell_module", raster = FALSE, pt_size = 0.5, palcolor = GJx_colors)
  
  CellDimPlot(seu, reduction = "UMAP_LSI_ATAC",  theme = "theme_blank", legend.position = "right", group_by = "cell_module", raster = FALSE, pt_size = 0.5, palcolor = GJx_colors)
  
}

# S2F
{
  Stacked_VlnPlot(seurat_object = seu, features = c("CD3E", "PAX5"), 
                  x_lab_rotate = TRUE, colors_use = GJx_colors, plot_spacing = 0.3,
                  vln_linewidth = 0.5)
}




