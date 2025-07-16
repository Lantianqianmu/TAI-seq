options(repos="http://mirrors.tuna.tsinghua.edu.cn/CRAN")
options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor")
# devtools::install_local("/home/zeemeeuw/Setup/ArchR-1.0.2.zip", ref = "master")
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
library(Cairo)
library(openxlsx)
library(readxl)
library(circlize)
source("/home/zeemeeuw/YangLab/JoINT-seq/QC/code/utils.R")

mmjoint_colors <- c("Quiescent" = "#bbdefb", 
                    "Priming_1" =  "#519ec6", 
                    "Priming_2" = "#3971ff",
                    "Activation_1" = "#a9df91", 
                    "Activation_2" = "#7eca72", 
                    "Activation_3" = "#48b352", 
                    "Effector_1" = "#ffb669", 
                    "Effector_2" = "#ff882e", 
                    "Effector-memory_transition_1" = "#cda7c7",
                    "Effector-memory_transition_2" = "#9e79b6",
                    "Effector-memory_transition_3" = "#7c4d9f",
                    "Memory_1" = "#ef999c", 
                    "Memory_2" = "#eb5e63",
                    "Memory_3" = "#e53a46")
time_colors <-  c("d3" = "#35a153", "d5" = "#1f6933", "d8" = "#f26a11", "d14" = "#B22222" , "d30" = "#f0f0f0", "R5" = "#85B0FF", "Y5" = "#3648A2")

# set workingdir
RDSdir <- "/home/zeemeeuw/data/mmJoint/RNA/data/lm_ym/" # output RDS file dir
ATACdir <- "/home/zeemeeuw/data/mmJoint/ATAC/overlapped_harmony/" # ATAC root dir


seu <- readRDS(paste0(RDSdir, "lm_merged_overlapped_harmony.rds"))
proj <- loadArchRProject(ATACdir)

setwd("/home/zeemeeuw/YangLab/JoINT-seq/Manuscript/fig_data/fig_S5_data")

# S5B
{
  meta <- readxl::read_excel("joint.Lm-overlapped-clean.seurat.meta.xlsx")
  tcr_paired <- readxl::read_excel("Lm-TCR4.57-clonotypes-paired in rna.xlsx") %>% 
    dplyr::mutate(tissue = str_extract(Tissue_time, "(.*)-", group = 1)) %>% 
    group_by(CDRa_CDRb, timepoint, tissue, barcode) %>% 
    distinct(CDRa_CDRb, timepoint, tissue, barcode, .keep_all = TRUE) %>% 
    ungroup()
  tcr_all <- readxl::read_excel("Lm-TCR4.57-clonotypes-all in rna.xlsx")
  
  bc_paired <- tcr_paired$barcode
  bc_a <- tcr_all$barcode[is.na(tcr_all$chainB)]
  bc_b <- tcr_all$barcode[is.na(tcr_all$chainA)]
  
  meta <- meta %>% dplyr::mutate(has_paired_tcr = if_else(cellnames %in% tcr_paired$barcode, "Y", "N"),
                                 has_all_tcr = if_else(cellnames %in% tcr_all$barcode, "Y", "N"),
                                 tcr = case_when(cellnames %in% bc_paired ~ "paired",
                                                 cellnames %in% bc_a ~ "A only",
                                                 cellnames %in% bc_b ~ "B only",
                                                 .default = "undetected"))
  
  meta_summary <- meta %>% 
    group_by(sample) %>% 
    dplyr::summarise(paired_tcr_prop = sum(has_paired_tcr == "Y")/n(), 
                     all_tcr_prop = sum(has_all_tcr == "Y")/n(), 
                     cell_number = n(),
                     sum_paired = sum(tcr == "paired"),
                     prop_paired = sum(tcr == "paired")/n(),
                     prop_a_only = sum(tcr == "A only")/n(),
                     prop_b_only = sum(tcr == "B only")/n(),
                     undetected = 1 - prop_paired - prop_a_only - prop_b_only)
  
  
  samples <- c("Spleen-WTx", "Spleen-d3SL", "Spleen-d3SL2", "Spleen-d5SL", "Spleen-d5SL2", "Spleen-d8SL", "Spleen-d8SL2", "Spleen-d14SL", "Spleen-d14SL2", "Spleen-d14SL3", "Spleen-d30SL", "Spleen-d30SL2", "Spleen-R5SL", "Spleen-R5SL2", "Spleen-R5r-SLL1", "Spleen-R5r3-SLL3", "Spleen-Y5SL", "Spleen-Y5SLr3",
               "Ln-WTx", "Ln-d3SL", "Ln-d3SL2", "Ln-d5SL", "Ln-d5SL2", "Ln-d8SL", "Ln-d8SL2", "Ln-d14SL", "Ln-d14SL2", "Ln-d14SL3", "Ln-d30SL", "Ln-d30SL2",  "Ln-R5SL", "Ln-R5SL2", "Ln-R5r-SLL1", "Ln-R5r3-SLL3", "Ln-Y5SL", "Ln-Y5SLr3", 
               "Blood-WTx", "Blood-d8x", "Blood-d8x2", "Blood-d14x", "Blood-d14x2", "Blood-d30BV", "Blood-d30BV2", "Blood-R5BV", "Blood-R5BV2", "Blood-Y5BV",  
               "Liver-WTx", "Liver-d5L", "Liver-d5L2",  "Liver-d8x", "Liver-d8x2", "Liver-d14x", "Liver-d14x2", "Liver-d14r-SLL1", "Liver-d14r3-SLL3", "Liver-d30BV", "Liver-d30BV2", "Liver-R5BV", "Liver-R5BV2", "Liver-R5r-YR", "Liver-Y5BV",  "Liver-Y5r-YR")
  
  
  meta_summary$sample <- factor(meta_summary$sample, levels = samples)
  meta_summary <- meta_summary %>% arrange(sample)
  meta_summary$id <- consecutive_id(meta_summary$sample)
  meta_summary$x <- meta_summary$id - 0.5
  
  
  
  
  # set color
  col_nfrags <- colorRampPalette(c("#B6BFEB", "#8694D8", "#5668C3", "#3648A2", "#1C2B77"))(62)
  col_TSS <- colorRampPalette(c("#B6DDEB", "#86C2D8", "#56A6C3", "#3685A2", "#1C5F77"))(62)
  col_numis <- colorRampPalette(c("#F7B4B4", "#ED8080", "#DF4A4A", "#BE2A2A", "#8C1616"))(62)
  col_ngenes <- colorRampPalette(c("#E4EBB6", "#CDD886", "#B4C356", "#93A236", "#6B771C"))(62)
  col_cellnum <- colorRampPalette(c("#B6EBBA", "#86D88B", "#56C35D", "#36A23D", "#1C7722"))(62)
  col_TCR <- colorRampPalette(c("#FFC782", "#FFA12F", "#FF851D", "#ff6200", "#FF4E0A", "#CF3D00"))(62)
  
  
  
  
  # Initialize Circos plot
  height_sample_names <- 0.05
  height_nfrags <- 0.08
  height_TSS <- 0.08
  height_numis <- 0.08
  height_ngenes <- 0.08
  height_cellnum <- 0.08
  height_TCR <- 0.15
  
  
  start_degree <- 90  
  circos.par(gap.degree = 0, 
             start.degree = start_degree, 
             track.margin = c(0.01, 0.01),
             canvas.xlim =  c(-1.5, 1.5),
             canvas.ylim =  c(-1.5, 1.5))
  circos.initialize(sectors = c("sample", "blank"), xlim = c(0, 62), sector.width = c(0.75, 0.25))
  
  # sample names
  circos.track(sectors = factor(c("sample", "blank"), levels = c("sample", "blank")),
               bg.border = NA,
               bg.col = c(sample = NA, blank = NA),
               track.height = height_sample_names,
               ylim = c(0, 1),
               panel.fun = function(x, y) {
                 if(CELL_META$sector.index == "sample"){
                   for(pos in seq(0.5, 61.5, by = 1)) {
                     value = meta_summary$sample[which(pos == meta_summary$x)]
                     circos.text(pos, 0, labels = value, 
                                 facing = "reverse.clockwise",
                                 adj = c(1, 0.5),
                                 cex = 0.8,
                                 niceFacing = TRUE,
                                 col = "black") # col = col_nfrags[pos+0.5]
                     
                   }
                 }
               })
  
  # nFrags
  circos.track(sectors = factor(c("sample", "blank"), levels = c("sample", "blank")),
               bg.border = NA,
               bg.col = c(sample = "#F2F4FE", blank = NA),
               track.height = height_nfrags,
               ylim = c(0, 4),
               panel.fun = function(x, y) {
                 if(CELL_META$sector.index == "sample"){
                   for(pos in seq(0.5, 61.5, by = 1)) {
                     value = log10(meta$nFrags[meta$sample == meta_summary$sample[which(pos == meta_summary$x)]])
                     circos.violin(value, pos, border = NA, col = col_nfrags[pos+0.5], violin_width = 1,  cex = 0.2)
                     
                   }
                   circos.yaxis(
                     side = "left",
                     at = c(0, 4),
                     labels = TRUE,
                     labels.cex = 0.8,
                     tick.length = 0.2
                   )
                 }
               })
  
  # TSS
  circos.track(sectors = factor(c("sample", "blank"), levels = c("sample", "blank")),
               bg.border = NA,
               bg.col = c(sample = "#F2FBFE", blank = NA),
               track.height = height_TSS,
               ylim = c(0, 30),
               panel.fun = function(x, y) {
                 if(CELL_META$sector.index == "sample"){
                   for(pos in seq(0.5, 61.5, by = 1)) {
                     value = meta$TSSEnrichment[meta$sample == meta_summary$sample[which(pos == meta_summary$x)]]
                     circos.violin(value, pos, border = NA, col = col_TSS[pos+0.5], violin_width = 1,  cex = 0.2)
                     
                   }
                   circos.yaxis(
                     side = "left",
                     at = c(0, 20),
                     labels = TRUE,
                     labels.cex = 0.8,
                     tick.length = 0.2
                   )
                 }
               })
  
  # nCount_RNA
  circos.track(sectors = factor(c("sample", "blank"), levels = c("sample", "blank")),
               track.height = height_numis,
               bg.border = NA,
               bg.col = c(sample = "#FFEEEE", blank = NA),
               ylim = c(0, 4), 
               panel.fun = function(x, y) {
                 if(CELL_META$sector.index == "sample"){
                   for(pos in seq(0.5, 61.5, by = 1)) {
                     value = log10(meta$nCount_RNA[meta$sample == meta_summary$sample[which(pos == meta_summary$x)]])
                     circos.violin(value, pos, border = NA, col = col_numis[pos+0.5], violin_width = 1, cex = 0.2)
                   }
                   circos.yaxis(
                     side = "left",
                     at = c(0, 4),
                     labels = TRUE,
                     labels.cex = 0.8,
                     tick.length = 0.2
                   )
                 }
                 
               })
  
  # nFeature_RNA
  circos.track(sectors = factor(c("sample", "blank"), levels = c("sample", "blank")),
               track.height = height_ngenes,
               bg.border = NA,
               bg.col = c(sample = "#FCFEF2",  blank = NA),
               ylim = c(0, 4), 
               panel.fun = function(x, y) {
                 if(CELL_META$sector.index == "sample"){
                   for(pos in seq(0.5, 61.5, by = 1)) {
                     value = log10(meta$nFeature_RNA[meta$sample == meta_summary$sample[which(pos == meta_summary$x)]])
                     circos.violin(value, pos, border = NA, col = col_ngenes[pos+0.5], violin_width = 1, cex = 0.2)
                   }
                   circos.yaxis(
                     side = "left",
                     at = c(0, 4),
                     labels = TRUE,
                     labels.cex = 0.8,
                     tick.length = 0.2
                   )
                 }
                 
               })
  
  # cell numbers
  circos.track(sectors = factor(c("sample", "blank"), levels = c("sample", "blank")),
               track.height = height_cellnum,
               bg.border = NA,
               bg.col = c(sample = "#F2FEF3", blank = NA),
               ylim = c(0, 5), panel.fun = function(x, y) {
                 if(CELL_META$sector.index == "sample"){
                   for(pos in seq(0.5, 61.5, by = 1)) {
                     value = log10(meta_summary$cell_number[meta_summary$sample == meta_summary$sample[which(pos == meta_summary$x)]])
                     circos.barplot(value, pos, border = NA, col = col_cellnum[pos+0.5])
                   }
                   circos.yaxis(
                     side = "left",
                     at = c(0, 5),
                     labels = TRUE,
                     labels.cex = 0.8,
                     tick.length = 0.2
                   )
                 }
               })
  
  # TCR numbers
  circos.track(sectors = factor(c("sample", "blank"), levels = c("sample", "blank")),
               track.height = height_TCR,
               bg.border = NA,
               bg.col = c(sample = "#FBE6BF", blank = NA),
               ylim = c(0, 4), panel.fun = function(x, y) {
                 if(CELL_META$sector.index == "sample"){
                   for(pos in seq(0.5, 61.5, by = 1)) {
                     value = log10(meta_summary$sum_paired[meta_summary$sample == meta_summary$sample[which(pos == meta_summary$x)]])
                     circos.barplot(value, pos, border = NA, col = col_TCR[pos+0.5])
                   }
                   circos.yaxis(
                     side = "left",
                     at = c(0, 4),
                     labels = TRUE,
                     labels.cex = 0.8,
                     tick.length = 0.2
                   )
                 }
               })
  
  
  
  # 结束绘图
  circos.clear()
  
  
}

# S5C
{
  CellDimPlot(seu, reduction = "UMAP_LSI_Combined_batch3.3",  theme = "theme_blank", legend.position = "right", 
              group_by = "cell_type", raster = TRUE, pt_size = 0.1)
}


# S5D
{
  CellDimPlot(seu, reduction = "UMAP_LSI_Combined_batch3.3",  theme = "theme_blank", legend.position = "right", 
              group_by = "timepoint", raster = TRUE, pt_size = 0.1)
  
  CellDimPlot(seu, reduction = "UMAP_LSI_Combined_batch3.3",  theme = "theme_blank", legend.position = "right", 
              group_by = "tissue", raster = TRUE, pt_size = 0.1, palcolor = tissue_colors)
}

# S5E Ro/Re
{
  seu$group <- paste0(seu$timepoint, "_", seu$tissue)
  data <- as.data.frame.array(table(seu$cell_type, seu$group))
  chi_sq_test <- chisq.test(data)
  exp <- chi_sq_test$expected
  Roe <- data/exp
  Roe[Roe > 5] <- 5
  
  pheatmap::pheatmap(Roe[,c("WT_Ln", "WT_Spleen", "WT_Liver", "WT_Blood",
                            "d3_Ln", "d3_Spleen",
                            "d5_Ln", "d5_Spleen", "d5_Liver", 
                            "d8_Ln", "d8_Spleen", "d8_Liver", "d8_Blood",
                            "d14_Ln", "d14_Spleen", "d14_Liver", "d14_Blood", 
                            "d30_Ln", "d30_Spleen", "d30_Liver", "d30_Blood", 
                            "R5_Ln", "R5_Spleen", "R5_Liver", "R5_Blood", 
                            "Y5_Ln", "Y5_Spleen", "Y5_Liver", "Y5_Blood")],
                     border_color = NA,
                     show_colnames = TRUE,
                     show_rownames = TRUE,
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     fontsize_row = 8,
                     fontsize_col = 8,
                     # legend_breaks = c(0, 1),
                     # breaks = seq(0, 1, length.out = 256),
                     # filename = "All_Ro_Re_celltype.pdf",
                     width = 8, height = 8,
                     scale = "none",
                     angle_col = 45,
                     display_numbers = TRUE, number_color = "black", fontsize_number = 6,
                     color = choose_colorset("oranges", 256))
}






