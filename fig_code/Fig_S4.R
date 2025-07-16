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
options(future.globals.maxSize = 1000*1024^2)
source("/home/zeemeeuw/YangLab/JoINT-seq/QC/code/utils.R")
setwd("/home/zeemeeuw/YangLab/JoINT-seq/Manuscript/fig_data/fig_S4_data")

RDSdir <- "/home/zeemeeuw/data/mmJoint_WTx/RNA/rds/" # output RDS file dir
ATACdir <- "/home/zeemeeuw/data/mmJoint_WTx/ATAC/" # ATAC root dir
seu <- readRDS(paste0(RDSdir, "WTx_overlapped.rds"))
proj <- loadArchRProject(ATACdir)


wt_colors <- c("C1_Il7r_Tcf7" = "#bbdefb", 
               "C2_Klra3_Tcf7" = "#3361A5",
               "C3_Klrk1_Eomes" = "#a9df91", 
               "C4_Mki67_Mga" = "#48b352", 
               "C5_Havcr2_Nfatc1" = "#ff882e", 
               "C6_Myb_Smarcc1" = "#e53a46")
time_colors <-  c("WT" = "#ABD8E6", "d3" = "#35a153", "d5" = "#1f6933", "d8" = "#f26a11", "d14" = "#B22222" , "d30" = "#B17BA6", "R5" = "#85B0FF", "Y5" = "#3648A2")
tissue_colors <-  c("Spleen" = "#36A23D", "Ln" = "#2884E7", "Liver" = "#817ab9", "Blood" = "#EE3B00")



# S4A
{
  CellDimPlot(seu, reduction = "UMAP_LSI_Combined",  theme = "theme_blank", legend.position = "right", group_by = "cell_module", raster = FALSE, pt_size = 0.5, palcolor = wt_colors)

  CellDimPlot(seu, reduction = "UMAP_LSI_RNA",  theme = "theme_blank", legend.position = "right", group_by = "cell_module", raster = FALSE, pt_size = 0.5, palcolor = wt_colors)

  CellDimPlot(seu, reduction = "UMAP_LSI_ATAC",  theme = "theme_blank", legend.position = "right", group_by = "cell_module", raster = FALSE, pt_size = 0.5, palcolor = wt_colors)
  
}

# S4B
{
  CellDimPlot(seu, reduction = "UMAP_LSI_Combined",  theme = "theme_blank", legend.position = "right", group_by = "tissue", raster = FALSE, pt_size = 0.5, palcolor = tissue_colors)
}

# S4C
{
  proj$tissue <- str_extract(proj$Sample, "(.*)-", group = 1)
  plotFragmentSizes(ArchRProj = proj, pal = tissue_colors, groupBy = "tissue")
  plotTSSEnrichment(ArchRProj = proj, pal = tissue_colors, groupBy = "tissue")
}

# S4D
{
  pdata <- data.frame(umis = seu$nCount_RNA, genes = seu$nFeature_RNA, tissue = seu$tissue)
  pdata$tissue <- factor(pdata$tissue, levels = c("Ln", "Spleen", "Liver", "Blood"))
  ggplot(pdata, aes(x = tissue)) +   
    geom_violin(aes(fill = tissue, y = log10(genes)), alpha = 0.5, trim = TRUE, scale = "width") +
    geom_boxplot(aes(fill = tissue, y = log10(genes)), alpha = 1, show.legend = FALSE, width = 0.3, outlier.shape = NA) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, size = 12, hjust = 1, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.ticks.x = element_blank()) +
    ylim(c(1, 4)) +
    scale_fill_manual(values = c("Spleen" = "#36A23D", "Ln" = "#2884E7", "Liver" = "#817ab9", "Blood" = "#EE3B00"))
}

# S4F

{
  markersGenes <- getMarkerFeatures(ArchRProj = proj, useMatrix = "GeneExpressionMatrix", groupBy = "cell_module", 
                                    bias = c("Gex_nUMI"), testMethod = "wilcoxon")
  markerList <- getMarkers(markersGenes, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")
  markerList2 <- lapply(markerList, function(x) {
    df <- tibble(as.data.frame(x))
    df %>% arrange(desc(Log2FC)) %>% slice_head(n = 8) %>% dplyr::select(name)
  })
  markersdf <- rbind(markerList2$`C1_Il7r_Tcf7`, markerList2$`C2_Klra3_Tcf7`, markerList2$`C3_Klrk1_Eomes`, markerList2$`C4_Mki67_Mga`, markerList2$`C5_Havcr2_Nfatc1`, markerList2$`C6_Myb_Smarcc1`)
  genes <- unique(markersdf$name)
  
  # genes <- c("Tcf7", "Sell", "Il7r", "Mki67", "Gzma", "Gzmb", "Tox", "Havcr2", "Klrk1", "Klrc1", "Klra5", "Klra9", "Klra3", "Ikzf2", "Hdac9", "Myb")
  
  # plotMarkerHeatmap(seMarker = markersGenes, cutOff = "FDR <= 0.05 & Log2FC >= 1", transpose = FALSE, returnMatrix = FALSE, clusterCols = FALSE)
  hm <- plotMarkerHeatmap(seMarker = markersGenes, cutOff = "FDR <= 1 & Log2FC >= 0.0", transpose = FALSE, returnMatrix = TRUE, clusterCols = FALSE, binaryClusterRows = TRUE)
  hm <- binarysort_heatmap(hm[genes,], cutOff = 1.0, clusterCols = FALSE, scale = FALSE, invert = FALSE)$mat
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
                     width = 2.5, height = 4)

  
}

# S4E
{
  markersGenes <- getMarkerFeatures(ArchRProj = proj, useMatrix = "GeneScoreMatrix", groupBy = "cell_module", 
                                    bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
  markerList <- getMarkers(markersGenes, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")
  markerList2 <- lapply(markerList, function(x) {
    df <- tibble(as.data.frame(x))
    df %>% arrange(desc(Log2FC)) %>% slice_head(n = 10) %>% dplyr::select(name)
  })
  markersdf <- rbind(markerList2$`C1_Il7r_Tcf7`, markerList2$`C2_Klra3_Tcf7`, markerList2$`C3_Klrk1_Eomes`, markerList2$`C4_Mki67_Mga`, markerList2$`C5_Havcr2_Nfatc1`, markerList2$`C6_Myb_Smarcc1`)
  genes <- unique(markersdf$name)
  
  # genes <- c("Tcf7", "Sell", "Il7r", "Mki67", "Gzma", "Gzmb", "Tox", "Havcr2", "Klrk1", "Klrc1", "Klra5", "Klra9", "Klra3", "Ikzf2", "Hdac9", "Myb")
  
  # plotMarkerHeatmap(seMarker = markersGenes, cutOff = "FDR <= 0.05 & Log2FC >= 1", transpose = FALSE, returnMatrix = FALSE, clusterCols = FALSE)
  hm <- plotMarkerHeatmap(seMarker = markersGenes, cutOff = "FDR <= 1 & Log2FC >= 0.0", transpose = FALSE, returnMatrix = TRUE, clusterCols = FALSE, binaryClusterRows = TRUE)
  hm <- binarysort_heatmap(hm[genes,], cutOff = 1.0, clusterCols = FALSE, scale = FALSE, invert = FALSE)$mat
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


# S4G
{
  motifs <- c("Tcf7", "Zeb1", "Foxp3")
  
  markerMotifs <- getFeatures(proj, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
  markerMotifs <- grep("z:", markerMotifs, value = TRUE)
  markerMotifs <- str_remove(markerMotifs, "z:")
  
  
  
  markersMotif <- getMarkerFeatures(ArchRProj = proj, useMatrix = "MotifMatrix", groupBy = "cell_module", 
                                    bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon", useSeqnames = "z")
  
  idx1 <- unlist(getMarkers(markersMotif, cutOff = "FDR <= 0.05 & MeanDiff > 1"))$idx
  idx2 <- rowData(markersMotif)$idx[which(rowData(markersMotif)$name %in% markerMotifs)]
  idx <- unique(c(idx1, idx2))
  
  
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





# S4H
{
  data_longer <- readxl::read_excel("Lm-TCR4.57-clonotypes-paired in rna.xlsx")
  data <- data_longer %>% 
    dplyr::filter(timepoint == "WT") %>% 
    dplyr::mutate(tissue = str_extract(Tissue_time, "(.*)-", group = 1)) %>% 
    group_by(CDRa_CDRb, timepoint, tissue, barcode) %>% 
    distinct(CDRa_CDRb, timepoint, tissue, barcode, .keep_all = TRUE) %>% 
    dplyr::filter(!str_detect(CDR3a, "_") & !str_detect(CDR3b, "_")) %>% 
    ungroup()
  data_sum <- data %>% left_join(seu@meta.data[,c("cellnames", "cell_module")], by = join_by("barcode" == "cellnames")) %>% 
    group_by(CDRa_CDRb, cell_module) %>% 
    dplyr::mutate(n = n()) %>% 
    ungroup() 
  exp_cells <- data_sum %>% dplyr::filter(CDRa_CDRb == "CAVSTNYGSSGNKLIF++CASGRTGGAYAEQFF" | CDRa_CDRb == "CALRGNNNNAPRF++CASGRTGGAYAEQFF") %>%  pull(barcode)
  seu$tcr_module <- case_when(seu$cellnames %in% exp_cells ~ "expanded", .default = "others")
  proj$tcr_module <- seu@meta.data[proj$cellNames,]$tcr_module
  
  markersGenes <- getMarkerFeatures(ArchRProj = proj, useMatrix = "GeneExpressionMatrix", groupBy = "tcr_module", 
                                    bias = c("Gex_nUMI"), testMethod = "ttest")
  markerList <- getMarkers(markersGenes, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")
  markerList2 <- lapply(markerList, function(x) {
    df <- tibble(as.data.frame(x))
    df %>% arrange(desc(Log2FC)) %>% slice_head(n = 30) %>% dplyr::select(name)
  })
  genes <- unique(unlist(markerList2))
  hm <- plotMarkerHeatmap(seMarker = markersGenes, cutOff = "FDR <= 0.05 & Log2FC >= 0.5", transpose = FALSE, returnMatrix = TRUE, clusterCols = FALSE, binaryClusterRows = TRUE, plotLog2FC = TRUE)
  hm <- binarysort_heatmap(hm[genes,], cutOff = 1.0, clusterCols = FALSE, scale = FALSE, invert = FALSE)$mat
  pheatmap::pheatmap(hm[genes,],
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
                     width = 2.5, height = 4)
  
}












