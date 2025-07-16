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
library(gghalves)
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
tissue_colors <-  c("Spleen" = "#36A23D", "Ln" = "#2884E7", "Liver" = "#817ab9", "Blood" = "#EE3B00")

# set workingdir
RDSdir <- "/home/zeemeeuw/data/mmJoint/RNA/data/lm_ym/" # output RDS file dir
ATACdir <- "/home/zeemeeuw/data/mmJoint/ATAC/overlapped_harmony/" # ATAC root dir
weightList <- list(paste0(ATACdir, "ImputeWeights/Impute-Weights-Rep-1"),
                   paste0(ATACdir, "ImputeWeights/Impute-Weights-Rep-2"))

seu <- readRDS(paste0(RDSdir, "lm_merged_overlapped_harmony.rds"))
proj <- loadArchRProject(ATACdir)

setwd("/home/zeemeeuw/YangLab/JoINT-seq/Manuscript/fig_data/fig_S8_data")


# S8J
{
  load(paste0(RDSdir, "lm_merged_overlapped_harmony_DORC3.RData"))
  gpcorr.use <- gpcorr.use_raw[gpcorr.use_raw$rObs > 0,]
  gpcorr.use <- gpcorr.use[gpcorr.use$pvalZ <= 0.05,]
  dorcGenes <- sort(unique(gpcorr.use$Gene))
  
  load(paste0(RDSdir, "lm_merged_overlapped_harmony_mat_smoothed.RData"))
  
  
  # Ln_re
  trajectory <- "Ln.trajectory_re"
  cells <- proj$cellNames[!is.na(getCellColData(proj, select = trajectory, drop = TRUE))]
  pseudotime <- getCellColData(proj, select = trajectory, drop = TRUE)[!is.na(getCellColData(proj, select = trajectory, drop = TRUE))]
  archdata <- data.frame(cells = proj$cellNames, pseudotime = getCellColData(proj, select = trajectory, drop = TRUE))
  rownames(archdata) <- archdata$cells
  data <- data.frame(seu@reductions$UMAP_LSI_Combined_batch3.3@cell.embeddings) 
  data$barcode <- rownames(data)
  data <- cbind(data, seu@meta.data)
  data$pseudotime <- archdata$pseudotime[rownames(data)]
  data_sub <- data[cells,]
  
  
  idx <- which(rowMeans(dorcMat_sm_norm - rnaMat_sm_norm) > 0)
  norm_list <- filter_matrix(dmat = dorcMat_sm[idx,cells], rmat = rnaMat_sm[idx,cells], scale_max_abs_val = 10, gene_corr_cutoff = 0, scaling = 'minmax', ncore = 12)
  ndmat <- norm_list[["ndmat"]]
  nrmat <- norm_list[["nrmat"]]
  
  embedding <- seu@reductions$UMAP_LSI_Combined_batch3.3@cell.embeddings[cells,]
  embedding_arr_length_list <- calc_velocity(ndmat, nrmat, embedding, max_per_cell = 10, metric = 'cosine', cutoff = 1)
  embedding_arr_length <- embedding_arr_length_list[[1]]
  arr_length <- embedding_arr_length_list[[2]]
  
  arrow_data <- smooth_arrows(embedding, embedding_arr_length, smooth_w = 30, min_count = 5, coef = 2, draw_all = FALSE)
  
  xrange <- max(data$UMAPCombinedbatch33_1) - min(data$UMAPCombinedbatch33_1)
  yrange <- max(data$UMAPCombinedbatch33_2) - min(data$UMAPCombinedbatch33_2)
  data_sub$exp <- colMeans(dorcMat_sm_norm[idx,cells])
  data_sub$obs <- colMeans(rnaMat_sm_norm[idx,cells])
  data_sub$residual <- data_sub$exp - data_sub$obs
  data_sub$residual <- ifelse(data_sub$residual > 0, data_sub$residual, 0)
  
  data_sub$norm_arrow_length <- arr_length
  ggplot(data, aes(x = UMAPCombinedbatch33_1, y = UMAPCombinedbatch33_2)) +
    geom_scattermore(aes(x = UMAPCombinedbatch33_1, y = UMAPCombinedbatch33_2), color = "grey90", show.legend = FALSE) +
    theme_classic() +
    geom_scattermore(data = data_sub, aes(x = UMAPCombinedbatch33_1, y = UMAPCombinedbatch33_2, color = cell_module), pointsize = 3, show.legend = FALSE) +
    geom_segment(data = arrow_data, aes(x = x_start, y = y_start, xend = x_end, yend = y_end), arrow = arrow(length = unit(3 , 'pt'), type = "closed"), lwd = 0.3, color = 'black') +
    coord_fixed(ratio = xrange/yrange) +
    scale_color_manual(values = mmjoint_colors)
  
  
  # Sp_re
  trajectory <- "Spleen.trajectory_re"
  cells <- proj$cellNames[!is.na(getCellColData(proj, select = trajectory, drop = TRUE))]
  pseudotime <- getCellColData(proj, select = trajectory, drop = TRUE)[!is.na(getCellColData(proj, select = trajectory, drop = TRUE))]
  archdata <- data.frame(cells = proj$cellNames, pseudotime = getCellColData(proj, select = trajectory, drop = TRUE))
  rownames(archdata) <- archdata$cells
  data <- data.frame(seu@reductions$UMAP_LSI_Combined_batch3.3@cell.embeddings) 
  data$barcode <- rownames(data)
  data <- cbind(data, seu@meta.data)
  data$pseudotime <- archdata$pseudotime[rownames(data)]
  data_sub <- data[cells,]
  
  
  idx <- which(rowMeans(dorcMat_sm_norm - rnaMat_sm_norm) > 0)
  norm_list <- filter_matrix(dmat = dorcMat_sm[idx,cells], rmat = rnaMat_sm[idx,cells], scale_max_abs_val = 10, gene_corr_cutoff = 0, scaling = 'minmax', ncore = 12)
  ndmat <- norm_list[["ndmat"]]
  nrmat <- norm_list[["nrmat"]]
  
  embedding <- seu@reductions$UMAP_LSI_Combined_batch3.3@cell.embeddings[cells,]
  embedding_arr_length_list <- calc_velocity(ndmat, nrmat, embedding, max_per_cell = 10, metric = 'cosine', cutoff = 1)
  embedding_arr_length <- embedding_arr_length_list[[1]]
  arr_length <- embedding_arr_length_list[[2]]
  
  arrow_data <- smooth_arrows(embedding, embedding_arr_length, smooth_w = 30, min_count = 5, coef = 2, draw_all = FALSE)
  
  xrange <- max(data$UMAPCombinedbatch33_1) - min(data$UMAPCombinedbatch33_1)
  yrange <- max(data$UMAPCombinedbatch33_2) - min(data$UMAPCombinedbatch33_2)
  data_sub$exp <- colMeans(dorcMat_sm_norm[idx,cells])
  data_sub$obs <- colMeans(rnaMat_sm_norm[idx,cells])
  data_sub$residual <- data_sub$exp - data_sub$obs
  data_sub$residual <- ifelse(data_sub$residual > 0, data_sub$residual, 0)
  
  data_sub$norm_arrow_length <- arr_length
  ggplot(data, aes(x = UMAPCombinedbatch33_1, y = UMAPCombinedbatch33_2)) +
    geom_scattermore(aes(x = UMAPCombinedbatch33_1, y = UMAPCombinedbatch33_2), color = "grey90", show.legend = FALSE) +
    theme_classic() +
    geom_scattermore(data = data_sub, aes(x = UMAPCombinedbatch33_1, y = UMAPCombinedbatch33_2, color = cell_module), pointsize = 3, show.legend = FALSE) +
    geom_segment(data = arrow_data, aes(x = x_start, y = y_start, xend = x_end, yend = y_end), arrow = arrow(length = unit(3 , 'pt'), type = "closed"), lwd = 0.3, color = 'black') +
    coord_fixed(ratio = xrange/yrange) +
    scale_color_manual(values = mmjoint_colors)
  
  
  # Liver_re
  trajectory <- "Liver.trajectory_re"
  cells <- proj$cellNames[!is.na(getCellColData(proj, select = trajectory, drop = TRUE))]
  pseudotime <- getCellColData(proj, select = trajectory, drop = TRUE)[!is.na(getCellColData(proj, select = trajectory, drop = TRUE))]
  archdata <- data.frame(cells = proj$cellNames, pseudotime = getCellColData(proj, select = trajectory, drop = TRUE))
  rownames(archdata) <- archdata$cells
  data <- data.frame(seu@reductions$UMAP_LSI_Combined_batch3.3@cell.embeddings) 
  data$barcode <- rownames(data)
  data <- cbind(data, seu@meta.data)
  data$pseudotime <- archdata$pseudotime[rownames(data)]
  data_sub <- data[cells,]
  
  idx <- which(rowMeans(dorcMat_sm_norm - rnaMat_sm_norm) > 0)
  norm_list <- filter_matrix(dmat = dorcMat_sm[idx,cells], rmat = rnaMat_sm[idx,cells], scale_max_abs_val = 10, gene_corr_cutoff = 0, scaling = 'minmax', ncore = 12)
  ndmat <- norm_list[["ndmat"]]
  nrmat <- norm_list[["nrmat"]]
  
  embedding <- seu@reductions$UMAP_LSI_Combined_batch3.3@cell.embeddings[cells,]
  embedding_arr_length_list <- calc_velocity(ndmat, nrmat, embedding, max_per_cell = 10, metric = 'cosine', cutoff = 1)
  embedding_arr_length <- embedding_arr_length_list[[1]]
  arr_length <- embedding_arr_length_list[[2]]
  
  arrow_data <- smooth_arrows(embedding, embedding_arr_length, smooth_w = 30, min_count = 5, coef = 2, draw_all = FALSE)
  
  xrange <- max(data$UMAPCombinedbatch33_1) - min(data$UMAPCombinedbatch33_1)
  yrange <- max(data$UMAPCombinedbatch33_2) - min(data$UMAPCombinedbatch33_2)
  data_sub$exp <- colMeans(dorcMat_sm_norm[idx,cells])
  data_sub$obs <- colMeans(rnaMat_sm_norm[idx,cells])
  data_sub$residual <- data_sub$exp - data_sub$obs
  data_sub$residual <- ifelse(data_sub$residual > 0, data_sub$residual, 0)
  
  data_sub$norm_arrow_length <- arr_length
  ggplot(data, aes(x = UMAPCombinedbatch33_1, y = UMAPCombinedbatch33_2)) +
    geom_scattermore(aes(x = UMAPCombinedbatch33_1, y = UMAPCombinedbatch33_2), color = "grey90", show.legend = FALSE) +
    theme_classic() +
    geom_scattermore(data = data_sub, aes(x = UMAPCombinedbatch33_1, y = UMAPCombinedbatch33_2, color = cell_module), pointsize = 3, show.legend = FALSE) +
    geom_segment(data = arrow_data, aes(x = x_start, y = y_start, xend = x_end, yend = y_end), arrow = arrow(length = unit(3 , 'pt'), type = "closed"), lwd = 0.3, color = 'black') +
    coord_fixed(ratio = xrange/yrange) +
    scale_color_manual(values = mmjoint_colors)
  
  
  # Blood_re
  trajectory <- "Blood.trajectory_re"
  cells <- proj$cellNames[!is.na(getCellColData(proj, select = trajectory, drop = TRUE))]
  pseudotime <- getCellColData(proj, select = trajectory, drop = TRUE)[!is.na(getCellColData(proj, select = trajectory, drop = TRUE))]
  archdata <- data.frame(cells = proj$cellNames, pseudotime = getCellColData(proj, select = trajectory, drop = TRUE))
  rownames(archdata) <- archdata$cells
  data <- data.frame(seu@reductions$UMAP_LSI_Combined_batch3.3@cell.embeddings) 
  data$barcode <- rownames(data)
  data <- cbind(data, seu@meta.data)
  data$pseudotime <- archdata$pseudotime[rownames(data)]
  data_sub <- data[cells,]
  
  
  idx <- which(rowMeans(dorcMat_sm_norm - rnaMat_sm_norm) > 0)
  norm_list <- filter_matrix(dmat = dorcMat_sm[idx,cells], rmat = rnaMat_sm[idx,cells], scale_max_abs_val = 10, gene_corr_cutoff = 0, scaling = 'minmax', ncore = 12)
  ndmat <- norm_list[["ndmat"]]
  nrmat <- norm_list[["nrmat"]]
  
  embedding <- seu@reductions$UMAP_LSI_Combined_batch3.3@cell.embeddings[cells,]
  embedding_arr_length_list <- calc_velocity(ndmat, nrmat, embedding, max_per_cell = 10, metric = 'cosine', cutoff = 1)
  embedding_arr_length <- embedding_arr_length_list[[1]]
  arr_length <- embedding_arr_length_list[[2]]
  
  arrow_data <- smooth_arrows(embedding, embedding_arr_length, smooth_w = 30, min_count = 5, coef = 2, draw_all = FALSE)
  
  xrange <- max(data$UMAPCombinedbatch33_1) - min(data$UMAPCombinedbatch33_1)
  yrange <- max(data$UMAPCombinedbatch33_2) - min(data$UMAPCombinedbatch33_2)
  data_sub$exp <- colMeans(dorcMat_sm_norm[idx,cells])
  data_sub$obs <- colMeans(rnaMat_sm_norm[idx,cells])
  data_sub$residual <- data_sub$exp - data_sub$obs
  data_sub$residual <- ifelse(data_sub$residual > 0, data_sub$residual, 0)
  
  data_sub$norm_arrow_length <- arr_length
  ggplot(data, aes(x = UMAPCombinedbatch33_1, y = UMAPCombinedbatch33_2)) +
    geom_scattermore(aes(x = UMAPCombinedbatch33_1, y = UMAPCombinedbatch33_2), color = "grey90", show.legend = FALSE) +
    theme_classic() +
    geom_scattermore(data = data_sub, aes(x = UMAPCombinedbatch33_1, y = UMAPCombinedbatch33_2, color = cell_module), pointsize = 3, show.legend = FALSE) +
    geom_segment(data = arrow_data, aes(x = x_start, y = y_start, xend = x_end, yend = y_end), arrow = arrow(length = unit(3 , 'pt'), type = "closed"), lwd = 0.3, color = 'black') +
    coord_fixed(ratio = xrange/yrange) +
    scale_color_manual(values = mmjoint_colors)
}


