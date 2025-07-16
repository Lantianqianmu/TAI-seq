# figure 3
library(ggplot2)
library(ggforce)
library(ggalluvial)
library(patchwork)
library(Seurat)
library(ArchR)
library(scplotter)
library(pheatmap)
library(scattermore)
library(BSgenome.Mmusculus.UCSC.mm10)
source("/home/zeemeeuw/YangLab/JoINT-seq/QC/code/utils.R")
setwd("/home/zeemeeuw/YangLab/JoINT-seq/Manuscript/fig_data/fig_3_data")

# set workingdir
RDSdir <- "/home/zeemeeuw/data/mmJoint/RNA/data/lm_ym/" # output RDS file dir
ATACdir <- "/home/zeemeeuw/data/mmJoint/ATAC/overlapped_harmony/" # ATAC root dir
seu <- readRDS(paste0(RDSdir, "lm_merged_overlapped_harmony.rds"))
proj <- loadArchRProject(ATACdir)

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

# Fig 3A: supervised trajectory
{
  proj$cell_state <- NA
  proj$cell_state[proj$cell_module %in% c("Quiescent")] <- "Quiescent" # nolint
  proj$cell_state[proj$cell_module %in% c("Priming_1" ,"Priming_2") & !(proj$timepoint %in% c("R5", "Y5"))] <- "pri_priming" 
  proj$cell_state[proj$cell_module %in% c("Activation_1", "Activation_2", "Activation_3") & !(proj$timepoint %in% c("R5", "Y5"))] <- "pri_Activation" 
  proj$cell_state[proj$cell_module %in% c("Effector_1", "Effector_2")  & !(proj$timepoint %in% c("R5", "Y5"))] <- "pri_Effector"
  proj$cell_state[proj$cell_module %in% c("Effector-memory_transition_1", "Effector-memory_transition_2", "Effector-memory_transition_3") & !(proj$timepoint %in% c("R5", "Y5"))] <- "pri_EMT"
  proj$cell_state[proj$cell_module %in% c("Memory_1", "Memory_2", "Memory_3") & !(proj$timepoint %in% c("R5", "Y5"))] <- "pri_Memory"
  
  proj$cell_state[proj$cell_module %in% c("Priming_1" ,"Priming_2") & proj$timepoint %in% c("R5", "Y5")] <- "re_priming"
  proj$cell_state[proj$cell_module %in% c("Activation_1", "Activation_2", "Activation_3") & proj$timepoint %in% c("R5", "Y5")] <- "re_Activation"
  proj$cell_state[proj$cell_module %in% c("Effector_1", "Effector_2")  & proj$timepoint %in% c("R5", "Y5")] <- "re_Effector"
  proj$cell_state[proj$cell_module %in% c("Effector-memory_transition_1", "Effector-memory_transition_2", "Effector-memory_transition_3") & proj$timepoint %in% c("R5", "Y5")] <- "re_EMT"
  proj$cell_state[proj$cell_module %in% c("Memory_1", "Memory_2", "Memory_3") & proj$timepoint %in% c("R5", "Y5")] <- "re_Memory"
  
  proj <- addTrajectory(proj, name = "all.pri.trajectory", groupBy = "cell_state",
                        trajectory = c("Quiescent", "pri_priming", "pri_Activation", "pri_Effector", "pri_EMT", "pri_Memory"),
                        embedding = "UMAP_Combined_batch3.3", force = TRUE)
  proj <- addTrajectory(proj, name = "all.re.trajectory", groupBy = "cell_state",
                        trajectory = c("re_priming", "re_Activation", "re_Effector", "re_EMT", "re_Memory"),
                        embedding = "UMAP_Combined_batch3.3", force = TRUE)
  
  p <- plotTrajectory(proj, embedding = "UMAP_Combined_batch3.3", colorBy = "cellColData", trajectory = "all.pri.trajectory", name = "all.pri.trajectory")
  p[[1]]
  p <- plotTrajectory(proj, embedding = "UMAP_Combined_batch3.3", colorBy = "cellColData", trajectory = "all.re.trajectory", name = "all.re.trajectory")
  p[[1]]
  
}

# Fig 3B: pseudotime heatmap
{
  # primary
  traj_pri_gsc  <- getTrajectory(ArchRProj = proj, name = "all.pri.trajectory", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
  plotTrajectoryHeatmap(traj_pri_gsc, pal = paletteContinuous(set = "whiteBlue"))
  traj_pri_gex  <- getTrajectory(ArchRProj = proj, name = "all.pri.trajectory", useMatrix = "GeneExpressionMatrix", log2Norm = TRUE)
  plotTrajectoryHeatmap(traj_pri_gex, pal = paletteContinuous(set = "blueYellow"))
  
  # recall
  traj_re_gsc  <- getTrajectory(ArchRProj = proj, name = "all.re.trajectory", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
  plotTrajectoryHeatmap(traj_pri_gsc, pal = paletteContinuous(set = "whiteBlue"))
  traj_re_gex  <- getTrajectory(ArchRProj = proj, name = "all.re.trajectory", useMatrix = "GeneExpressionMatrix", log2Norm = TRUE)
  plotTrajectoryHeatmap(traj_pri_gex, pal = paletteContinuous(set = "blueYellow"))
  
}

# Fig 3C: correlated pseudotime heatmaps of TFs
{
  # primary
  traj_pri_gex  <- getTrajectory(ArchRProj = proj, name = "all.pri.trajectory", useMatrix = "GeneExpressionMatrix", log2Norm = TRUE)
  traj_pri_tf  <- getTrajectory(ArchRProj = proj, name = "all.pri.trajectory", useMatrix = "MotifMatrix", log2Norm = TRUE)
  cor_gex_tf <- correlateTrajectories(traj_pri_gex, traj_pri_tf)
  idxToRemove <- grep(pattern = "deviations", x = cor_gex_tf[["correlatedMappings"]]$name2)
  if(length(idxToRemove > 1)){
    cor_gex_tf[["correlatedMappings"]] <- cor_gex_tf[["correlatedMappings"]][-idxToRemove,]
  }
  traj_pri_gex2 <- traj_pri_gex[cor_gex_tf[["correlatedMappings"]]$name1, ]
  traj_pri_tf2 <- traj_pri_tf[cor_gex_tf[["correlatedMappings"]]$name2, ]
  
  trajCombined <- traj_pri_gex2
  assay(trajCombined, withDimnames = FALSE) <- t(apply(assay(traj_pri_gex2), 1, scale)) + t(apply(assay(traj_pri_tf2), 1, scale))
  
  combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
  rowOrder <- match(rownames(combinedMat), rownames(traj_pri_gex2))
  
  ht1 <- plotTrajectoryHeatmap(traj_pri_gex2,  pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0, rowOrder = rowOrder)
  ht2 <- plotTrajectoryHeatmap(traj_pri_tf2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)
  ComplexHeatmap::draw(ht1 + ht2)
  
  # recall
  traj_re_gex  <- getTrajectory(ArchRProj = proj, name = "all.re.trajectory", useMatrix = "GeneExpressionMatrix", log2Norm = TRUE)
  traj_re_tf  <- getTrajectory(ArchRProj = proj, name = "all.re.trajectory", useMatrix = "MotifMatrix", log2Norm = TRUE)
  cor_gex_tf <- correlateTrajectories(traj_re_gex, traj_re_tf)
  idxToRemove <- grep(pattern = "deviations", x = cor_gex_tf[["correlatedMappings"]]$name2)
  if(length(idxToRemove > 1)){
    cor_gex_tf[["correlatedMappings"]] <- cor_gex_tf[["correlatedMappings"]][-idxToRemove,]
  }
  traj_re_gex2 <- traj_re_gex[cor_gex_tf[["correlatedMappings"]]$name1, ]
  traj_re_tf2 <- traj_re_tf[cor_gex_tf[["correlatedMappings"]]$name2, ]
  
  trajCombined <- traj_re_gex2
  assay(trajCombined, withDimnames = FALSE) <- t(apply(assay(traj_re_gex2), 1, scale)) + t(apply(assay(traj_re_tf2), 1, scale))
  
  combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
  rowOrder <- match(rownames(combinedMat), rownames(traj_re_gex2))
  
  ht1 <- plotTrajectoryHeatmap(traj_re_gex2,  pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0, rowOrder = rowOrder)
  ht2 <- plotTrajectoryHeatmap(traj_re_tf2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)
  ComplexHeatmap::draw(ht1 + ht2)
}

# Fig 3D: chromatin potential
{
  
  load(paste0(RDSdir, "lm_merged_overlapped_harmony_DORC3.RData"))
  gpcorr.use <- gpcorr.use_raw[gpcorr.use_raw$rObs > 0,]
  gpcorr.use <- gpcorr.use[gpcorr.use$pvalZ <= 0.05,]
  dorcGenes <- sort(unique(gpcorr.use$Gene))
  
  load(paste0(RDSdir, "lm_merged_overlapped_harmony_mat_smoothed.RData"))
  
  
  time_colors <-  c("WT" = "#ABD8E6", "d3" = "#35a153", "d5" = "#1f6933", "d8" = "#f26a11", "d14" = "#B22222" , "d30" = "#B17BA6", "R5" = "#85B0FF", "Y5" = "#3648A2")
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
  
  # Ln
  trajectory <- "Ln.trajectory"
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

  
  # Sp
  trajectory <- "Spleen.trajectory"
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
  data_sub$exp <- colMeans(dorcMat_sm_norm[idx,cells, drop = FALSE])
  data_sub$obs <- colMeans(rnaMat_sm_norm[idx,cells, drop = FALSE])
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
  
  
  # Liver
  trajectory <- "Liver.trajectory"
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
  ggsave("chrom_potential_Liver_2.pdf", width = 4, height = 4)
  
  # Blood
  trajectory <- "Blood.trajectory"
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
  ggsave("chrom_potential_Blood_2.pdf", width = 4, height = 4)
  

}

# Fig 3E: residual classification
{
  load(paste0(RDSdir, "Spleen_DORC_RNA_residual.RData"))
  
  cols_to_keep <- colMaxs(as.matrix(res_dorc)) > 0.5 & colMaxs(as.matrix(res_rna)) > 0.5
  res <- res_dorc[,cols_to_keep] - res_rna[,cols_to_keep]
  res_positive <- t(res[, colSums(res <= 0) != 100])
  res_dorc_used <- t(res_dorc[,rownames(res_positive)])
  res_rna_used <- t(res_rna[,rownames(res_positive)])
  line_colors <- c("DORC" = "#314A8C",
                   "Motif" = "#D03542",
                   "RNA" = "#F6CE48")
  checkpoint <- list(c(1, 10), c(26, 35), c(51, 60), c(82, 91))
  aggr_range <- list(c(1, 25), c(26, 50), c(51, 75), c(76, 100))
  color_fill <- c("#B6BFEB", "#C7E9B4", "#FBE6BF", "#F7B4B4")
  trajectory <- "Spleen.trajectory"
  genes_predi <- list("checkpoint_1" = c("Azin1", "Brip1", "Casp3", "Dock9", "Fasl", "Havcr2", "Hist1h1d", 
                                         "Hist1h1e", "Hist1h3b", "Igf2bp3", "Lamc1", "Lgals1", "Ly6e", 
                                         "Mcm6", "Myb", "Ptma", "Ptpn13", "S100a6", "Tbc1d2", "Txn1", 
                                         "Zfp281"),
                      "checkpoint_2" = c("Acss2", "Adap1", "Alcam", "Anxa2", "As3mt", "Asf1b", "Atad5", 
                                         "Atp2b4", "Atxn1", "Atxn7l1", "Birc5", "Brip1", "Car5b", "Casp3", 
                                         "Ccdc88c", "Ccl4", "Ccl5", "Ccnf", "Ccr2", "Ccr5", "Cd48", "Cdc20b", 
                                         "Cdca2", "Chsy1", "Cldnd2", "Clspn", "Cmklr1", "Cpped1", "Crip1", 
                                         "Crot", "Cx3cr1", "Cxcr3", "Cyba", "Cyth4", "Dock5", "Dtl", "Dusp5", 
                                         "E2f7", "Ehbp1l1", "Entpd1", "Epas1", "Fam129b", "Fasl", "Fbxl2", 
                                         "Fcgr2b", "Fhl2", "Fkbp5", "Gas7", "Gata3", "Gna15", "Gnai2", 
                                         "Gnptab", "Gsap", "Gzma", "Gzmk", "Havcr2", "Hist1h1b", "Hist1h1c", 
                                         "Hist1h1d", "Hist1h1e", "Hist1h2ab", "Hist1h2ae", "Hist1h3b", 
                                         "Hist1h4d", "Hmgb2", "Hopx", "Id2", "Ifitm10", "Igf2bp3", "Igf2r", 
                                         "Il18r1", "Il18rap", "Inpp4a", "Iqgap2", "Itga1", "Itga2", "Itga4", 
                                         "Itgal", "Itgax", "Itgb1", "Itgb2", "Kcnab2", "Kif13b", "Kif18b", 
                                         "Kif23", "Klrc1", "Klrc2", "Klrg1", "Klrk1", "Lamc1", "Lfng", 
                                         "Lgals1", "Lgals3", "Lmnb1", "Lrrk1", "Map7d1", "Mast2", 
                                         "Mcm6", "Mki67", "Mknk2", "N4bp1", "Ncapd2", "Ncapg2", "Ndc80", 
                                         "Neat1", "Nod1", "Pcna", "Pkm", "Plxnc1", "Pmaip1", "Prc1", "Prdm1", 
                                         "Prex1", "Ptpn13", "Ptpre", "Ptprj", "Rac2", "Rad51b", "Rbm47", 
                                         "Ric1", "Runx1", "Runx2", "S100a10", "Serpinb6b", "Serpinb9", 
                                         "Sh3bgrl3", "Slamf7", "Slc20a1", "Smpdl3b", "Sntb2", "Snx10", 
                                         "Ssh1", "Stk32c", "Stmn1", "Stx11", "Sun1", "Tbc1d10b", "Tbc1d2", 
                                         "Tbkbp1", "Tbx21", "Tcf3", "Tmem154", "Tmem163", "Tmpo", "Ttc39c", 
                                         "Ttc7b", "Ubash3b", "Uhrf1", "Xdh", "Zeb2", "Zgrf1", "Zmiz1"),
                      "checkpoint_3" = c("As3mt", "Bcl6", "Chd3", "Dennd4a", "Etf1", "Etv3", "Ifrd1", "Junb", 
                                         "Tbkbp1", "Tnfaip3" ,"Vezf1"),
                      "checkpoint_4" = c("Ccr9", "Hmgb1", "E2f7",
                                         "Hsp90aa1", "Hspa5", "Ifi27l2a", "Iigp1", "Lrig1", "Ly6a", "Nap1l1", 
                                         "Ncapd2", "Nedd4l", "Pcgf5", "Pcna", "Pdk1", "Pmepa1", "Ptma", 
                                         "Rangap1", "Rapgef4", "Rps11", "Rps26", 
                                         "Socs3", "Tnfrsf26", "Ttc28", "Tuba1b", 
                                         "Ulbp1", "Zbp1"))
  genes_maint <- list("checkpoint_1" = c("Acss2", "Actn1", "Acvr1b", "Adgre5", "Afap1", "Ahnak", "Als2cl", 
                                         "Ap1ar", "Arpc5", "Atp1b1", "B4galnt1", "Boll", "Ccm2", "Ccr7", 
                                         "Ccr9", "Cd27", "Cmah", "Cnp", "Dapl1", "Dcaf17", "Dip2b", "Dym", 
                                         "Eif4a2", "Elovl5", "Fam169b", "Fam78a", "Fbxo32", "Galnt6", 
                                         "Gimap3", "Gne", "Gtf2i", "Hdac7", "Hlcs", "Hspbap1", "Hvcn1", 
                                         "Ifitm10", "Ifnar1", "Ift80", "Igf1r", "Ikbkb", "Il4ra", "Il6ra", 
                                         "Iqgap2", "Itga6", "Itpkb", "Klf3", "Klhl6", "Lfng", "Lrig1", 
                                         "Lyst", "Maml3", "Map4k4", "Map7d1", "Mettl8", "Mgat5", "Mrps28", 
                                         "Mtss1", "N4bp2", "Nfatc3", "Nsg2", "Pacs2", "Pdk1", "Qprt", 
                                         "R3hdm1", "Rab3ip", "Ralgps2", "Ramp1", "Rapgef4", "Rasgrp2", 
                                         "Rfx1", "Rgcc", "Rgs10", "Rhobtb2", "Rnf144a", "Rpl14", "Rpl17", 
                                         "Rpl21", "Rpl34", "Rpl41", "Rplp0", "Rps26", "Rps6ka5", "Sesn1", 
                                         "Sgms1", "Siah1a", "Sipa1l2", "Ski", "Skp1a", "Sla2", "Slamf6", 
                                         "Slc12a7", "Slc26a11", "Slc44a1", "Slc6a19", "Smad7", "Smc4", 
                                         "Spag9", "Ssbp3", "St6gal1", "St8sia1", "Stk38", "Sun1", "Tanc1", 
                                         "Tbc1d10b", "Tcf3", "Tdrp", "Tet1", "Tet3", "Tgfbr2", "Tgfbr3", 
                                         "Tmem154", "Tmie", "Tmsb10", "Tnfrsf14", "Tnfrsf26", "Tnrc6c", 
                                         "Trat1", "Trib2", "Trio", "Tspan13", "Ttc28", "Ttc3", "Ubash3a", 
                                         "Ucp2", "Usp28", "Vamp1", "Vezf1", "Zdhhc14", "Zfp664", "Zmiz1"),
                      "checkpoint_2"  = c("Airn", "Arf4", "Arhgap25", "Batf", "Bcas3", "BE692007", "Boll", 
                                          "Cd8b1", "Cers6", "Chd3", "Crim1", "Cyb5a", "Dpp4", "Epsti1", 
                                          "Etf1", "Ets1", "Etv3", "Fas", "Fchsd2", "Gbp7", "Ggact", "Gramd4", 
                                          "Hbs1l", "Hif1a", "Hnrnpd", "Hsd11b1", "Hspa4", "Jmy", "Jun", 
                                          "Junb", "Kat6b", "Klf3", "Klhdc2", "Klhl6", "Kmt2e", "Lncpint", 
                                          "Lyst", "Map3k8", "March3", "Myb", "Naa16", "Nap1l1", "Ncl", 
                                          "Nfatc3", "Nfe2l2", "Nfkbiz", "Pacs1", "Phc2", "Plac8", "Rab2a", 
                                          "Relb", "Rinl", "Rps11", "Rps26", "Rsbn1l", "Slc26a11", "Spag9", 
                                          "Spred2", "Srgn", "Ssbp3", "Stat2", "Stat4", "Tbl1xr1", "Tnfrsf26", 
                                          "Tom1l2", "Trat1", "Ubr2", "Vezf1", "Vps54", "Zc3h12d", "Zfp664", 
                                          "Zfp827"),
                      "checkpoint_3" = c("Actn4", "Adap1", "Agfg2", "Alcam", "Anxa2", "Asf1b", "Atp5b", 
                                         "Azin1", "Baz1b", "Birc5", "Brip1", "Casp3", "Cbx5", "Ccnf", 
                                         "Cd2bp2", "Cdca2", "Cit", "Clic1", "Clspn", "Cnn2", "Crip1", 
                                         "Ctcf", "Cxcr3", "Dip2c", "Dtl", "E2f3", "E2f7", "Eef1d", "Ehbp1l1", 
                                         "Epas1", "Fam129b", "Fhl2", "Fkbp5", "Gas7", "Gnptab", "Gzmb", 
                                         "Gzmk", "Havcr2", "Hif1a", "Hmgb1", "Hmgb2", "Hnrnpd", "Hnrnpl", 
                                         "Hsp90aa1", "Hspa4", "Hspa5", "Id2", "Igf2bp3", "Il12rb2", "Kif18b", 
                                         "Kif23", "Lgals1", "Lgals3", "Lockd", "Map7d1", "Mast2", "Memo1", 
                                         "Mki67", "Nap1l1", "Nasp", "Ncapg2", "Ndc80", "Nedd4l", "Nkg7", 
                                         "Nucks1", "Nusap1", "Pcgf5", "Pcna", "Pkm", "Plxnc1", "Pmaip1", 
                                         "Prc1", "Ptch1", "Ptma", "Rad51b", "Ralgps2", "Rangap1", "Rfwd3", 
                                         "Rpl35", "Rps11", "Rps6ka5", "Serpinb6b", "Serpinb9", "Sh3bgrl3", 
                                         "Slc38a2", "Spred2", "Srsf2", "Srsf3", "Srsf7", "Ssrp1", "Stil", 
                                         "Stk32c", "Stmn1", "Sun1", "Taf15", "Tbc1d10b", "Tbc1d2", "Tcf3", 
                                         "Tex2", "Tmpo", "Top2a", "Tpm4", "Traf2", "Trib2", "Ttc28", "Tuba1b", 
                                         "Txn1", "Ubb", "Uhrf1", "Ulbp1", "Vim", "Zfp664", "Zgrf1"),
                      "checkpoint_4"  = c("Acss2", "Arhgap25", "Btg2", "Cdc20b", "Cmklr1", "Cx3cr1", 
                                          "Fbxl2", "Gna15", "Gnptab", "Gzma", "Igf2bp3", "Ikzf3", "Klrc1", 
                                          "Klrc2", "Klrg1", "Lamc1", "Malt1", "N4bp1", "Osbpl3", "Smpdl3b", 
                                          "Snx10", "Srgap3", "Sun1", "Tbkbp1", "Tbx21", "Tcf3"))
  for(i in seq_along(checkpoint)){
    time_1 <- checkpoint[[i]][1]
    time_2 <- checkpoint[[i]][2]
    genes_cp <- (res_positive[,time_2] - res_positive[,time_1] > 0.015) & (res_positive[,time_2] > 0)
    
    mat_cp <- res_positive[genes_cp,]
    gene_predi <- genes_predi[[i]]
    n_predi <- length(gene_predi)
    gene_maint <- genes_maint[[i]]
    n_maint <- length(gene_maint)
    
    mat_cp <- mat_cp[order(mat_cp[,time_2], decreasing = TRUE),]
    mat_predi <- mat_cp[gene_predi,]
    mat_predi <- mat_predi[order(mat_predi[,time_2], decreasing = TRUE),]
    mat_maint <- mat_cp[gene_maint,]
    mat_maint <- mat_maint[order(mat_maint[,time_2], decreasing = TRUE),]
    mat_cp <- rbind(mat_predi, mat_maint)
    # heatmap
    pheatmap::pheatmap(mat_cp,
                       cluster_rows = FALSE, 
                       cluster_cols = FALSE,
                       show_rownames = FALSE,
                       treeheight_col = 0,
                       treeheight_row = 10,
                       border_color = NA,
                       clustering_method = "average",
                       scale = "none",
                       gaps_row = n_predi,
                       color = choose_colorset("solar_extra"),
                       filename = paste0("Spleen_residual_checkpoint_", i, "_heatmap_predi_", n_predi, "_maint_", n_maint, ".pdf"),
                       breaks = seq(from = -0.3, to = 0.3, length.out = 256),
                       # fontsize_row = 8,
                       width = 4, height = 6)

    # aggregated line plot, predisposition
    range <- aggr_range[[i]]
    aggr_predi_dorc <- colSums(res_dorc_used[gene_predi,])/length(gene_predi)
    aggr_predi_rna <- colSums(res_rna_used[gene_predi,])/length(gene_predi)
    data_dorc <- tibble(score = aggr_predi_dorc, pseudotime = 1:100, type = "DORC")
    data_rna <- tibble(score = aggr_predi_rna, pseudotime = 1:100, type = "RNA")
    data <- bind_rows(data_dorc, data_rna)
    ratio <- 45/(max(data$score)-min(data$score))
    pdf(paste0("Spleen_checkpoint_", i, "_lineplot_aggregated.pdf"), width = 4, height = 3)
    p1 <- ggplot(data) +
      theme_bw() +
      geom_rect(xmin = time_1, xmax = time_2, ymin = 0, ymax = 1, fill = color_fill[i], color = "NA", alpha = 0.8) +
      geom_line(aes(x = pseudotime, y = score, color = type), linewidth = 1) +
      scale_color_manual(values = line_colors) +
      scale_x_continuous(limits = c(0, 100), expand = c(0.025, 0.025)) +
      scale_y_continuous(limits = c(NA, NA)) +
      theme(panel.grid = element_blank(),
            axis.text.x = element_text(size = 10, color = "black"),
            axis.text.y = element_text(size = 10, color = "black"),
            axis.title.x = element_text(size = 13),
            axis.title.y = element_text(size = 13)) +
      ggtitle("Predisposition") +
      theme(plot.title = element_text(hjust = 0.5)) +
      coord_fixed(ratio) +
      ylab("Normalized counts") +
      xlab("Pseudotime")
    print(p1)
    data <- data %>% dplyr::filter(pseudotime <= range[2] & pseudotime > range[1])
    ratio <- 45/ (max(data$score)-min(data$score))
    p2 <- ggplot(data) +
      theme_bw() +
      geom_rect(xmin = time_1, xmax = time_2, ymin = 0, ymax = 1, fill = color_fill[i], color = "NA", alpha = 0.8) +
      geom_line(aes(x = pseudotime, y = score, color = type), linewidth = 1) +
      scale_color_manual(values = line_colors) +
      scale_x_continuous(limits = range, expand = c(0.025, 0.025), breaks = range, labels = range) +
      scale_y_continuous(limits = c(NA, NA)) +
      theme(panel.grid = element_blank(),
            axis.text.x = element_text(size = 10, color = "black"),
            axis.text.y = element_text(size = 10, color = "black"),
            axis.title.x = element_text(size = 13),
            axis.title.y = element_text(size = 13)) +
      ggtitle("Predisposition") +
      theme(plot.title = element_text(hjust = 0.5)) +
      coord_fixed(ratio) +
      ylab("Normalized counts") +
      xlab("Pseudotime")
    print(p2)
    # aggregated line plot, maintenance
    aggr_maint_dorc <- colSums(res_dorc_used[gene_maint,])/length(gene_maint)
    aggr_maint_rna <- colSums(res_rna_used[gene_maint,])/length(gene_maint)
    data_dorc <- tibble(score = aggr_maint_dorc, pseudotime = 1:100, type = "DORC")
    data_rna <- tibble(score = aggr_maint_rna, pseudotime = 1:100, type = "RNA")
    data <- bind_rows(data_dorc, data_rna)
    ratio <- 45/ (max(data$score)-min(data$score))
    p1 <- ggplot(data) +
      theme_bw() +
      geom_rect(xmin = time_1, xmax = time_2, ymin = 0, ymax = 1, fill = color_fill[i], color = "NA", alpha = 0.8) +
      geom_line(aes(x = pseudotime, y = score, color = type), linewidth = 1) +
      scale_color_manual(values = line_colors) +
      scale_x_continuous(limits = c(0, 100), expand = c(0.025, 0.025)) +
      scale_y_continuous(limits = c(NA, NA)) +
      theme(panel.grid = element_blank(),
            axis.text.x = element_text(size = 10, color = "black"),
            axis.text.y = element_text(size = 10, color = "black"),
            axis.title.x = element_text(size = 13),
            axis.title.y = element_text(size = 13)) +
      ggtitle("Maintenance") +
      theme(plot.title = element_text(hjust = 0.5)) +
      coord_fixed(ratio) +
      ylab("Normalized counts") +
      xlab("Pseudotime")
    print(p1)
    data <- data %>% dplyr::filter(pseudotime <= range[2] & pseudotime > range[1])
    ratio <- 45/ (max(data$score)-min(data$score))
    p2 <- ggplot(data) +
      theme_bw() +
      geom_rect(xmin = time_1, xmax = time_2, ymin = 0, ymax = 1, fill = color_fill[i], color = "NA", alpha = 0.8) +
      geom_line(aes(x = pseudotime, y = score, color = type), linewidth = 1) +
      scale_color_manual(values = line_colors) +
      scale_x_continuous(limits = range, expand = c(0.025, 0.025), breaks = range, labels = range) +
      scale_y_continuous(limits = c(NA, NA)) +
      theme(panel.grid = element_blank(),
            axis.text.x = element_text(size = 10, color = "black"),
            axis.text.y = element_text(size = 10, color = "black"),
            axis.title.x = element_text(size = 13),
            axis.title.y = element_text(size = 13)) +
      ggtitle("Maintenance") +
      theme(plot.title = element_text(hjust = 0.5)) +
      coord_fixed(ratio) +
      ylab("Normalized counts") +
      xlab("Pseudotime")
    print(p2)
    dev.off()
  }

}

# Fig 3F: line plots
{
  load(paste0(RDSdir, "Spleen_DORC_RNA_residual.RData"))
  line_colors <- c("DORC" = "#314A8C",
                   "Motif" = "#D03542",
                   "RNA" = "#F6CE48")
  cols_to_keep <- colMaxs(as.matrix(res_dorc)) > 0.5 & colMaxs(as.matrix(res_rna)) > 0.5
  res <- res_dorc[,cols_to_keep] - res_rna[,cols_to_keep]
  res_positive <- t(res[, colSums(res <= 0) != 100])
  res_dorc_used <- t(res_dorc[,rownames(res_positive)])
  res_rna_used <- t(res_rna[,rownames(res_positive)])
  trajectory <- "Spleen.trajectory"
  cells <- proj$cellNames[!is.na(getCellColData(proj, select = trajectory, drop = TRUE))]
  pseudotime <- getCellColData(proj, select = trajectory, drop = TRUE)[!is.na(getCellColData(proj, select = trajectory, drop = TRUE))]
  genes <- c("Ccr7", "Ly6e", "St6gal1", "Myb", "Klf3", "Klrg1", "Tbx21", "Gzmb", "Mki67", "Etv3", "Ccr9")
  plot_list <- list()
  for(gene in genes){
    data_dorc <- data.frame(score = res_dorc_used[gene,], pseudotime = 1:100, type = "DORC")
    data_rna <- data.frame(score = res_rna_used[gene,], pseudotime = 1:100, type = "RNA")
    data <- rbind(data_dorc, data_rna)
    data$type <- factor(data$type, levels = c("DORC", "RNA"))
    
    p <- ggplot(data) +
      theme_bw() +
      geom_line(aes(x = pseudotime, y = score, color = type), linewidth = 1) +
      scale_color_manual(values = line_colors) +
      scale_x_continuous(limits = c(0, 100), expand = c(0.025, 0.025)) +
      scale_y_continuous(limits = c(0, 1), expand = c(0.025, 0.025)) +
      theme(panel.grid = element_blank(),
            axis.text.x = element_text(size = 10, color = "black"),
            axis.text.y = element_text(size = 10, color = "black"),
            axis.title.x = element_text(size = 13),
            axis.title.y = element_text(size = 13)) +
      ggtitle(gene) +
      theme(plot.title = element_text(hjust = 0.5)) +
      coord_fixed(50) +
      ylab("Normalized counts") +
      xlab("Pseudotime")
    plot_list[[gene]] <- p
  }
  wrap_plots(plot_list, ncol = 2)
}

# Fig 3G: TF regulators of each chromatin checkpoint
{
  load(paste0(RDSdir, "Spleen_DORC_RNA_residual_checkpoint.RData"))
  load(paste0(RDSdir, "lm_merged_overlapped_harmony_GRN_list_by_tissue.RData"))
  grn <- grn_list$grn_sp
  
  # lollipop
  groups <- rep(c("predisposition", "maintenance"), 4)
  checkpoints <- c(1, 1, 2, 2, 3, 3, 4, 4)
  data <- lapply(seq_along(1:8), function(i){
    group <- groups[i]
    checkpoint <- checkpoints[i]
    genes <- rownames(cp_list[[i]]$residual)
    grn_used <- grn %>% dplyr::filter(DORC %in% genes) %>% 
      dplyr::filter(abs(Score) >= 0.8) %>% 
      group_by(Motif) %>% 
      dplyr::summarise(corr = mean(Corr), n = n()) %>% 
      dplyr::arrange(desc(n), .by_group = TRUE) %>% 
      dplyr::mutate(group = group, checkpoint = checkpoint)
  }) %>% Reduce("rbind", .) %>% 
    dplyr::mutate_at(vars(checkpoint), ~as.integer(.))
  cp1_pre <- c("Bach2", "Nfe2l2", "Atf6", "Nfkb2", "Tbx21")
  cp1_mai <- c("Zeb1", "Tcf3", "Zfp281", "Bhlhe40", "Tcf7")
  cp2_pre <-c("Bach2", "Tbx21", "Nfe2l2", "Arid5b", "Klf12")
  cp2_mai <-c("Tcf3", "Zeb2", "Zeb1", "Tcf7", "Nfkb2")
  cp3_pre <-c("E2f7", "Foxo3", "Foxp1", "Stat3", "Tcf7")
  cp3_mai <-c("Bach2", "Nfe2l2", "Arid5b", "Atf6", "Tet1")
  cp4_pre <-c("Tcf3", "Zeb1", "Nfkb2", "Zeb2", "Bach2")
  cp4_mai <-c("Tbx21", "Bach2", "Klf12", "Rara", "Trps1")
  
  pdata <- dplyr::bind_rows(dplyr::filter(data, ((group == "predisposition") & (checkpoint == 1) & (Motif %in% cp1_pre))),
                            dplyr::filter(data, ((group == "maintenance") & (checkpoint == 1) & (Motif %in% cp1_mai))),
                            dplyr::filter(data, ((group == "predisposition") & (checkpoint == 2) & (Motif %in% cp2_pre))),
                            dplyr::filter(data, ((group == "maintenance") & (checkpoint == 2) & (Motif %in% cp2_mai))),
                            dplyr::filter(data, ((group == "predisposition") & (checkpoint == 3) & (Motif %in% cp3_pre))),
                            dplyr::filter(data, ((group == "maintenance") & (checkpoint == 3) & (Motif %in% cp3_mai))),
                            dplyr::filter(data, ((group == "predisposition") & (checkpoint == 4) & (Motif %in% cp4_pre))),
                            dplyr::filter(data, ((group == "maintenance") & (checkpoint == 4) & (Motif %in% cp4_mai)))) %>% 
    dplyr::group_by(group, checkpoint) %>% 
    dplyr::arrange(dplyr::desc(n), .by_group = TRUE) %>% 
    group_by(group) %>% 
    dplyr::mutate(x = row_number(),
                  y = if_else(group == "predisposition", n, -n)) %>% 
    dplyr::mutate(ytext = if_else(group == "predisposition", n+10, -n-10)) %>% 
    dplyr::mutate(color = case_when(y > 30 ~ 30, y < -30 ~ -30, .default = y))
  motif_selected <- pdata %>% pull("Motif")
  data <- pdata
  write_tsv(data, "source_fig_3G.tsv", na = "")
  ggplot(data, aes(x = x, y = y, color = corr)) +
    theme_bw() +
    ylim(-80, 80) +
    ylab("Number of regulated genes, + = predi, - = maint") +
    geom_vline(xintercept = c(5.5, 10.5, 15.5), linetype = 2, color = 'grey90') +
    geom_hline(yintercept = 0, linetype = 1, color = 'grey90') +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(hjust = 0.5, vjust = 0.5, size = 8, color = "black"),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 10, color = "black"),
          panel.grid = element_blank()) +
    geom_segment(aes(x = x, xend = x, y = 0, yend = y)) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(aes(label = Motif), color = "black", size = 2, direction = "y") +
    scale_color_gradientn(colors = choose_colorset("imagine"), limits = c(-1, 1))
}


# Fig 3H: TF regulatory networks
{
  load(paste0(RDSdir, "Spleen_DORC_RNA_residual_checkpoint.RData"))
  load(paste0(RDSdir, "lm_merged_overlapped_harmony_GRN_list_by_tissue.RData"))
  grn <- grn_list$grn_sp
  library(ggraph)
  library(igraph)
  conn <- grn %>% 
    dplyr::filter(Motif %in% motif_selected) %>% 
    dplyr::mutate(Motif = paste0("Motif:", Motif)) %>% 
    dplyr::filter(abs(Score) >= 0.8) %>% 
    dplyr::mutate(group = if_else(Corr > 0, "Positive", "Negative")) %>% 
    dplyr::select(DORC, Motif, group)
  
  vertices_DORC <- conn %>%
    ungroup() %>% 
    dplyr::select(DORC) %>% 
    group_by(DORC) %>% 
    summarise(n = n()) %>% 
    rename_with(~c("name", "n")) %>% 
    dplyr::mutate(group = "DORC")
  
  vertices_TF <- conn %>%
    ungroup() %>% 
    dplyr::select(Motif) %>% 
    group_by(Motif) %>% 
    summarise(n = n()) %>% 
    rename_with(~c("name", "n")) %>% 
    dplyr::mutate(group = "Motif")
  vertices <- bind_rows(vertices_DORC, vertices_TF)
  
  conn <- conn %>% dplyr::select(DORC, Motif, group) %>% 
    rename_with(~c("from", "to", "group")) %>% 
    dplyr::mutate(value = 1)
  gr <- graph_from_data_frame(conn, vertices = vertices)
  
  set.seed(61)
  ggraph(gr, layout = "igraph", algorithm = 'fr') +
    geom_edge_link2(aes(edge_color = group), edge_alpha = 0.5, edge_width = 1, angle_calc = "along", linejoin = "round", check_overlap = TRUE) + 
    scale_edge_color_manual(values = c("Positive" = "#C09ADB", "Negative" = "#86d9ac")) +
    geom_node_point(aes(fill = group, size = n), stroke = NA, shape = 21, alpha = 0.5) +
    scale_size_continuous(range = c(4, 10)) +
    geom_node_text(aes(label = ifelse(group == "Motif", str_extract(as.character(name), ":(.*)", group = 1), "")), size = 5, color = "black") +
    theme_void() +
    scale_fill_manual(values = c(DORC = "#1563AA", "Motif" = "#A31D1D")) +
    guides(fill = guide_legend(override.aes = list(size = 5))) +
    expand_limits(x = c(-3, 3), y = c(-3, 3)) +
    coord_equal() +
    geom_node_text(aes(label = ifelse(group == "DORC", name, "")), size = 1, color = "black") 
}

# Fig 3J: TF RNA-DORC-MOtif lineplot
{
  line_colors <- c("DORC" = "#314A8C",
                   "Motif" = "#D03542",
                   "RNA" = "#F6CE48")
  load(paste0(RDSdir, "lm_merged_overlapped_harmony_mat_smoothed.RData"))
  plot_list <- list()
  cells <- proj$cellNames[!is.na(proj$Spleen.trajectory)]
  pseudotime <- proj$Spleen.trajectory[!is.na(proj$Spleen.trajectory)]
  for(gene in c("Bach2", "Tbx21", "Tcf7", "Bhlhe40")){
    data_motif <- data.frame(score = motifMat_sm_norm[gene, cells], pseudotime = pseudotime, type = "Motif")
    data_dorc <- data.frame(score = dorcMat_sm_norm[gene, cells], pseudotime = pseudotime, type = "DORC")
    data_rna <- data.frame(score = rnaMat_sm_norm[gene, cells], pseudotime = pseudotime, type = "RNA")
    
    Motif_fit <- loess(score ~ pseudotime, data = data_motif, span = 0.2)
    data_motif$smooth <- Motif_fit$fitted
    DORC_fit <- loess(score ~ pseudotime, data = data_dorc, span = 0.2)
    data_dorc$smooth <- DORC_fit$fitted
    RNA_fit <- loess(score ~ pseudotime, data = data_rna, span = 0.2)
    data_rna$smooth <- RNA_fit$fitted
    
    data <- rbind(data_motif, data_dorc, data_rna)
    data$type <- factor(data$type, levels = c("Motif", "DORC", "RNA"))
    
    p <- ggplot(data) +
      theme_bw() +
      # geom_point(aes(x = pseudotime, y = score, color = type), size = 0.5, alpha = 0.1, shape = 16) +
      geom_line(aes(x = pseudotime, y = smooth, color = type), size = 1) +
      # geom_smooth(aes(x = pseudotime, y = smooth, color = type), method = "loess", span = 0.2, size = 1, se = FALSE) +
      scale_color_manual(values = line_colors) +
      scale_x_continuous(limits = c(0, 100), expand = c(0.025, 0.025)) +
      scale_y_continuous(limits = c(0, 1), expand = c(0.025, 0.025)) +
      theme(panel.grid = element_blank(),
            axis.text.x = element_text(size = 10, color = "black"),
            axis.text.y = element_text(size = 10, color = "black"),
            axis.title.x = element_text(size = 13),
            axis.title.y = element_text(size = 13)) +
      ggtitle(gene) +
      theme(plot.title = element_text(hjust = 0.5)) +
      coord_fixed(50) +
      ylab("Normalized counts") +
      xlab("Pseudotime")
    plot_list[[gene]] <- p
  }
  wrap_plots(plot_list, ncol = 2)
}

# Fig 3K: Bach2 footpringing
{
  # proj <- addGroupCoverages(ArchRProj = proj, groupBy = "cell_module", force = TRUE)
  motifPositions <- getPositions(proj)
  tf_abbr <- c("Bach2")
  markerMotifs <- unlist(lapply(tf_abbr, function(x) grep(x, names(motifPositions), value = TRUE)))
  markerMotifs
  seFoot <- getFootprints(
    ArchRProj = proj, 
    useGroups = c("Quiescent", "Priming_1", "Priming_2", "Activation_1", "Activation_2", 
                  "Activation_3", "Effector_1", "Effector_2", "Effector-memory_transition_1", 
                  "Effector-memory_transition_2", "Effector-memory_transition_3", 
                  "Memory_1", "Memory_2", "Memory_3"),
    positions = motifPositions[markerMotifs], 
    groupBy = "cell_module"
  )
  
  plotFootprints(
    seFoot = seFoot,
    pal = mmjoint_colors,
    ArchRProj = proj, 
    flank = 200,
    flankNorm = 30,
    normMethod = "Subtract",
    plotName =  "Bach2",
    addDOC = FALSE,
    smoothWindow = 5,
    height = 4, width = 3,
    force = TRUE
  )
  
}


# Fig 3L: TF driver plot
{
  load(paste0(RDSdir, "lm_merged_overlapped_harmony_GRN_list_by_tissue.RData"))
  grn <- grn_list$grn_sp
  d <- grn %>% dplyr::filter(DORC == "Bach2") %>% 
    dplyr::mutate(isSig = if_else(abs(Score) >= 0.8 & abs(Corr) >= 0.7, "Yes", "No"))
  d$Label <- d$Motif
  d$Label[d$isSig %in% "No"] <- ""
  write_tsv(d, "source_fig_3L.tsv", na = "")
  
  ggplot(data = d, aes(x = Corr, y = Enrichment.log10P, color = isSig)) + 
    geom_hline(yintercept = -log10(0.01), color = "gray60", linetype = "dashed") + 
    geom_vline(xintercept = c(-0.7, 0.7), color = "gray60", linetype = "dashed") + 
    geom_point(aes(size = isSig)) + 
    theme_classic() + 
    scale_color_manual(values = c("gray66", "firebrick3")) + 
    scale_size_manual(values = c(1, 3)) +
    labs(y = "Enrichment log10 P", x = "DORC-TF Correlation", title = "Bach2") + 
    ylim(-10, 15) + 
    xlim(-1, 1) + 
    theme(legend.position = "none", 
          axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 10, color = "black"),
          axis.text.y = element_text(hjust = 0.5, vjust = 0.5, size = 10, color = "black"),
          plot.title = element_text(hjust = 0.5, face = "italic"), 
          panel.background = element_rect(fill = NA)) + 
    coord_fixed(2/25) +
    ggrepel::geom_text_repel(aes(label = Label), size = 5, color = "black", vjust = 0.5, hjust = 0.5, max.overlaps = 13)
}







