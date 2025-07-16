library(ggplot2)
library(ArchR)
library(Seurat)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Mmusculus.v79) # mm10
library(EnsDb.Hsapiens.v86) # hg38
library(patchwork)
library(scplotter)
library(scCustomize)
library(scattermore)
library(future)
library(readxl)
library(openxlsx)
options(future.globals.maxSize = 1000*1024^2)
source("/home/zeemeeuw/YangLab/JoINT-seq/QC/code/utils.R")

RDSdir <- "/home/zeemeeuw/data/HLX/RNA/rds/" # output RDS file dir
ATACdir <- "/home/zeemeeuw/data/HLX/ATAC_overlapped/" # ATAC root dir
seu <- readRDS(paste0(RDSdir, "HLX_overlapped.rds"))
proj <- loadArchRProject(ATACdir)
setwd("/home/zeemeeuw/YangLab/JoINT-seq/Manuscript/fig_data/fig_S3_data")

hs_colors <- c("C1_Mono_LYZ_CEBPA" = "#bbdefb",
               "C2_B_IGHM_STAT2" = "#3361A5",
               "C3_B_TNFRSF13B_POU2F2" = "#3971ff",
               "C4_CD4T_TCF7_TCF7" = "#48b352",
               "C5_CD4T_CCL5_EOMES" = "#ffb669",
               "C6_CD4T_CCR4_TCF7" = "#1f6933",
               "C7_CD8T_TCF7_E2F6" = "#a9df91",
               "C8_CD8T_DKK3_BACH1" = "#FFE060",
               "C9_CD8T_GZMB_EOMES" = "#ff882e",
               "C10_CD8T_KLRC1_BACH1" = '#E280BA',
               "C11_CD8T_IKZF2_EOMES" = "#e53a46",
               "C12_NK_NCAM1_MGA" = "#9e79b6",
               "C13_NKT_KLRF1_MGA" = "#5C0F71",
               "C14_MAIT_SLC4A10_RORA" = "#d6d6d6")

# S3A
{
  CellDimPlot(seu, reduction = "UMAP_LSI_Combined",  theme = "theme_blank", legend.position = "right", group_by = "cell_module", raster = FALSE, pt_size = 0.5, palcolor = hs_colors)

  CellDimPlot(seu, reduction = "UMAP_LSI_RNA",  theme = "theme_blank", legend.position = "right", group_by = "cell_module", raster = FALSE, pt_size = 0.5, palcolor = hs_colors)

  CellDimPlot(seu, reduction = "UMAP_LSI_ATAC",  theme = "theme_blank", legend.position = "right", group_by = "cell_module", raster = FALSE, pt_size = 0.5, palcolor = hs_colors)

}

# S3B
{
  plotFragmentSizes(ArchRProj = proj, pal = "#ff882e")
  plotTSSEnrichment(ArchRProj = proj, pal = "#ff882e")
  
}

# S3C
{
  pdata <- data.frame(umis = seu$nCount_RNA, genes = seu$nFeature_RNA, nfrags = proj@cellColData[colnames(seu),]$nFrags, frip = proj@cellColData[colnames(seu),]$FRIP)
  
  ggplot(pdata, aes(x = 1)) +   
    geom_violin(aes(y = frip), alpha = 0.5, trim = TRUE, scale = "width", fill = "#ff882e") +
    geom_boxplot(aes(y = frip), alpha = 1, show.legend = FALSE, width = 0.3, outlier.shape = NA, fill = "#ff882e") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, size = 12, hjust = 1, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.ticks.x = element_blank()) +
    ggtitle("FRiP")

  ggplot(pdata, aes(x = 1)) +   
    geom_violin(aes(y = log10(umis)), alpha = 0.5, trim = TRUE, scale = "width", fill = "#ff882e") +
    geom_boxplot(aes(y = log10(umis)), alpha = 1, show.legend = FALSE, width = 0.3, outlier.shape = NA, fill = "#ff882e") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, size = 12, hjust = 1, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.ticks.x = element_blank()) +
    ylim(c(1, 5)) +
    ggtitle("umis")

  ggplot(pdata, aes(x = 1)) +   
    geom_violin(aes(y = log10(genes)), alpha = 0.5, trim = TRUE, scale = "width", fill = "#ff882e") +
    geom_boxplot(aes(y = log10(genes)), alpha = 1, show.legend = FALSE, width = 0.3, outlier.shape = NA, fill = "#ff882e") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, size = 12, hjust = 1, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.ticks.x = element_blank()) +
    ylim(c(1, 4)) +
    ggtitle("genes")

}

# S3D
{
  #  write_tsv(pdata_joint, "source_fig_S3D.txt")
  pdata <- read_tsv("source_fig_S3D.txt")
  
  length(intersect(pdata_10X$CDR3, pdata_joint$CDR3))
  length(pdata_10X$CDR3)
  length(pdata_joint$CDR3)
}

# S3E-F
{

  tcr_merged_circ <- read_tsv("source_fig_S3E_circ.txt")
  tcr_merged_vmix <- read_tsv("source_fig_S3E_vmix.txt")
  
  # for TRA gene
  TRA_circ <- str_split_i(tcr_merged_circ$chainA, "\\+", 1)
  TRA_circ <- TRA_circ[TRA_circ != "NA"]
  TRA_circ <- as.data.frame(table(TRA_circ))
  colnames(TRA_circ) <- c("gene", "count_circ")
  TRA_circ$method <- "circ"
  TRA_circ$prop_circ <- TRA_circ$count_circ/sum(TRA_circ$count_circ)
  
  TRA_vmix <- str_split_i(tcr_merged_vmix$chainA, "\\+", 1)
  TRA_vmix <- TRA_vmix[TRA_vmix != "NA"]
  TRA_vmix <- as.data.frame(table(TRA_vmix))
  colnames(TRA_vmix) <- c("gene", "count_vmix")
  TRA_vmix$method <- "vmix"
  TRA_vmix$prop_vmix <- TRA_vmix$count_vmix/sum(TRA_vmix$count_vmix)
  
  TRA_data <- dplyr::full_join(TRA_circ, TRA_vmix, by = "gene")
  TRA_data[is.na(TRA_data)] <- 0
  fit <- summary(lm(count_circ~count_vmix, data = TRA_data))
  sqrt(fit$r.squared)
  
  
  ratio <- (max(TRA_data$prop_vmix)-min(TRA_data$prop_vmix))/(max(TRA_data$prop_circ)-min(TRA_data$prop_circ))
  ggplot(TRA_data, aes(x = prop_vmix, y = prop_circ)) +
    geom_abline(slope = 1 , intercept = 0, size = 1, linetype ="dashed", alpha = 0.7, color = "red") +
    geom_point(aes(fill = gene), show.legend = FALSE, size = 5, alpha = 0.7, shape = 21, color = "black") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 0, size = 10, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.title.x = element_text(size = 13),
          axis.title.y = element_text(size = 13),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 13)) +
    scale_fill_manual(values = choose_colorset("bubble", n = nrow(TRA_data))) +
    xlab("Proportion in Vmix") +
    ylab("Proportion in Circ") +
    ggtitle("TRA gene usage") +
    coord_fixed(ratio) +
    ggrepel::geom_label_repel(aes(label = gene),  max.overlaps = 1, box.padding = 0.5)

  
  
  # for TRB gene
  TRB_circ <- str_split_i(tcr_merged_circ$chainB, "\\+", 1)
  TRB_circ <- TRB_circ[TRB_circ != "NA"]
  TRB_circ <- as.data.frame(table(TRB_circ))
  colnames(TRB_circ) <- c("gene", "count_circ")
  TRB_circ$method <- "circ"
  TRB_circ$prop_circ <- TRB_circ$count_circ/sum(TRB_circ$count_circ)
  
  TRB_vmix <- str_split_i(tcr_merged_vmix$chainB, "\\+", 1)
  TRB_vmix <- TRB_vmix[TRB_vmix != "NA"]
  TRB_vmix <- as.data.frame(table(TRB_vmix))
  colnames(TRB_vmix) <- c("gene", "count_vmix")
  TRB_vmix$method <- "vmix"
  TRB_vmix$prop_vmix <- TRB_vmix$count_vmix/sum(TRB_vmix$count_vmix)
  
  TRB_data <- dplyr::full_join(TRB_circ, TRB_vmix, by = "gene")
  TRB_data[is.na(TRB_data)] <- 0
  fit <- summary(lm(count_circ~count_vmix, data = TRB_data))
  sqrt(fit$r.squared)
  
  ratio <- (max(TRB_data$prop_vmix)-min(TRB_data$prop_vmix))/(max(TRB_data$prop_circ)-min(TRB_data$prop_circ))
  ggplot(TRB_data, aes(x = prop_vmix, y = prop_circ)) +
    geom_abline(slope = 1 , intercept = 0, size = 1, linetype ="dashed", alpha = 0.7, color = "red") +
    geom_point(aes(fill = gene), show.legend = FALSE, size = 5, alpha = 0.7, shape = 21, color = "black") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 0, size = 10, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.title.x = element_text(size = 13),
          axis.title.y = element_text(size = 13),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 13)) +
    scale_fill_manual(values = choose_colorset("bubble", n = nrow(TRB_data))) +
    xlab("Proportion in Vmix") +
    ylab("Proportion in Circ") +
    ggtitle("TRB gene usage") +
    coord_fixed(ratio) +
    ggrepel::geom_label_repel(aes(label = gene),  max.overlaps = 1, box.padding = 0.5)

  
  # S3F
  tcr_merged_ob <- dplyr::full_join(tcr_merged_circ[c("barcode", "origin")], tcr_merged_vmix[c("barcode", "origin")], by = "barcode")
  
  tcr_merged_ob$`TCR detected by` <- "Both"
  tcr_merged_ob$`TCR detected by`[is.na(tcr_merged_ob$origin.x)] <- "Vmix only"
  tcr_merged_ob$`TCR detected by`[is.na(tcr_merged_ob$origin.y)] <- "Circ only"
  ggplot(tcr_merged_ob, aes(x = 1, fill = `TCR detected by`)) +
    theme_bw() +
    geom_bar(stat="count",width=0.5,position='fill') +
    ylab("proportion") + xlab(NULL) +
    theme(axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line = element_blank(),
          panel.grid = element_blank()) +
    scale_fill_manual(values = c("#1B46C0FF", "#1e88e5", "#bbdefb"))
  
  
}

# S3G
# GEX
{
  markersGenes <- getMarkerFeatures(ArchRProj = proj, useMatrix = "GeneExpressionMatrix", groupBy = "cell_module", 
                                    bias = c("Gex_nUMI"), testMethod = "wilcoxon")
  markerList <- getMarkers(markersGenes, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")
  markerList2 <- lapply(markerList, function(x) {
    df <- tibble(as.data.frame(x))
    df %>% arrange(desc(Log2FC)) %>% slice_head(n = 8) %>% dplyr::select(name)
  })
  
  
  
  genes <- unique(unlist(markerList2))
  
  # genes <- c("Tcf7", "Sell", "Il7r", "Mki67", "Gzma", "Gzmb", "Tox", "Havcr2", "Klrk1", "Klrc1", "Klra5", "Klra9", "Klra3", "Ikzf2", "Hdac9", "Myb")
  
  # plotMarkerHeatmap(seMarker = markersGenes, cutOff = "FDR <= 0.05 & Log2FC >= 1", transpose = FALSE, returnMatrix = FALSE, clusterCols = FALSE)
  hm <- plotMarkerHeatmap(seMarker = markersGenes, cutOff = "FDR <= 1 & Log2FC >= 0.0", transpose = FALSE, returnMatrix = TRUE, clusterCols = FALSE, binaryClusterRows = TRUE)
  hm <- binarysort_heatmap(hm[genes,], cutOff = 3.0, clusterCols = FALSE, scale = FALSE, invert = FALSE)$mat
  pheatmap::pheatmap(hm[,],
                     show_colnames = TRUE,
                     show_rownames = TRUE,
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     border_color = NA,
                     treeheight_col = 0,
                     treeheight_row = 6,
                     clustering_method = "ward.D2",
                     scale = "none",
                     color = choose_colorset("gothic"),
                     fontsize_row = 10,
                     angle_col = 45,
                     width = 2.5, height = 4)

  
}



# GSC
{
  markersGenes <- getMarkerFeatures(ArchRProj = proj, useMatrix = "GeneScoreMatrix", groupBy = "cell_module", 
                                    bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
  markerList <- getMarkers(markersGenes, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")
  markerList2 <- lapply(markerList, function(x) {
    df <- tibble(as.data.frame(x))
    df %>% arrange(desc(Log2FC)) %>% slice_head(n = 10) %>% dplyr::select(name)
  })
  genes <- unique(unlist(markerList2))
  
  # genes <- c("Tcf7", "Sell", "Il7r", "Mki67", "Gzma", "Gzmb", "Tox", "Havcr2", "Klrk1", "Klrc1", "Klra5", "Klra9", "Klra3", "Ikzf2", "Hdac9", "Myb")
  
  # plotMarkerHeatmap(seMarker = markersGenes, cutOff = "FDR <= 0.05 & Log2FC >= 1", transpose = FALSE, returnMatrix = FALSE, clusterCols = FALSE)
  hm <- plotMarkerHeatmap(seMarker = markersGenes, cutOff = "FDR <= 1 & Log2FC >= 0.0", transpose = FALSE, returnMatrix = TRUE, clusterCols = FALSE, binaryClusterRows = TRUE)
  hm <- binarysort_heatmap(hm[genes,], cutOff = 3, clusterCols = FALSE, scale = FALSE, invert = FALSE)$mat
  pheatmap::pheatmap(hm[,],
                     show_colnames = TRUE,
                     show_rownames = TRUE,
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     border_color = NA,
                     treeheight_col = 0,
                     treeheight_row = 6,
                     clustering_method = "ward.D2",
                     scale = "none",
                     color = choose_colorset("wolfgang_extra"),
                     fontsize_row = 10,
                     width = 2.5, height = 4)

  
}



# Motif
{
  
  markersMotif <- getMarkerFeatures(ArchRProj = proj, useMatrix = "MotifMatrix", groupBy = "cell_module", 
                                    bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon", useSeqnames = "z")
  
  markerList <- getMarkers(markersMotif, cutOff = "FDR <= 0.05 & MeanDiff > 0.3")
  
  markerList2 <- lapply(markerList, function(x) {
    df <- tibble(as.data.frame(x))
    df %>% arrange(desc(MeanDiff)) %>% slice_head(n = 10) %>% dplyr::select(idx)
  })
  
  idx <- unique(unlist(markerList2))
  
  
  hm <- assay(markersMotif)[idx,]
  rownames(hm) <- rowData(markersMotif)$name[idx]
  hm <- as.matrix(hm[,])
  hm <- sweep(hm - rowMeans(hm), 1, matrixStats::rowSds(hm), `/`)
  hm <- binarysort_heatmap(hm, scale = TRUE, clusterCols = FALSE, invert = FALSE, cutOff = 1)$mat
  # hm <- plotMarkerHeatmap(seMarker = markersMotif, cutOff = "FDR <= 0.05 & MeanDiff > 0.03", transpose = FALSE, returnMatrix = TRUE, clusterCols = FALSE, binaryClusterRows = TRUE)
  pheatmap::pheatmap(hm[,],
                     show_colnames = TRUE,
                     show_rownames = TRUE,
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     border_color = NA,
                     treeheight_col = 0,
                     treeheight_row = 6,
                     clustering_method = "ward.D2",
                     scale = "none",
                     color = choose_colorset("solar_extra"),
                     fontsize_row = 8,
                     width = 2.4, height = 4)

  
}










