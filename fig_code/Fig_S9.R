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

setwd("/home/zeemeeuw/YangLab/JoINT-seq/Manuscript/fig_data/fig_S9_data")




# S9A
{
  library(UpSetR)
  load(paste0(RDSdir, "Spleen_DORC_RNA_residual_checkpoint.RData"))
  gname <- lapply(cp_list, function(l){rownames(l$dorc)})
  df <- tibble(gene = gname$cp1_predi, count_1 = 1)
  for (i in 2:8) {
    cname <- paste0("count_", i)
    mdf <- tibble(gene = gname[[i]], !!cname := 1)
    df <- df %>% dplyr::full_join(mdf, by = "gene")
  }
  df <- df %>% replace(is.na(.), 0)
  df <- df %>% tibble::column_to_rownames(var = "gene")
  colnames(df) <- c(
    "checkpoint_1_predisposition", "checkpoint_1_maintenance",
    "checkpoint_2_predisposition", "checkpoint_2_maintenance", 
    "checkpoint_3_predisposition", "checkpoint_3_maintenance",  
    "checkpoint_4_predisposition", "checkpoint_4_maintenance"
  )
  
  upset(df,
        nsets = 8,
        sets = rev(colnames(df) ),
        nintersects= 100,
        keep.order = TRUE,
        main.bar.color = "#8f1f1e",
        matrix.color="black", 
        sets.bar.color= "#85B0FF",
        set_size.show = TRUE, 
        mb.ratio = c(0.5, 0.5))
}


# S9B
{
  # ---------------- trajectory Ln ---------------- #
  pdf("Ln.trajectory_dorc_gene_lineplot.pdf", width = 4, height = 3)
  cells <- proj$cellNames[!is.na(proj$Ln.trajectory)]
  pseudotime <- proj$Ln.trajectory[!is.na(proj$Ln.trajectory)]
  for(gene in c("Ccr7", "St6gal1", "Klf3", "Gzmb", "MKi67", "Klrg1", "Tbx21", "Ly6e", "Myb", "Etv3", "Ccr9")){
    data_dorc <- data.frame(score = dorcMat_sm_norm[gene, cells], pseudotime = pseudotime, type = "DORC")
    data_rna <- data.frame(score = rnaMat_sm_norm[gene, cells], pseudotime = pseudotime, type = "RNA")
    
    DORC_fit <- loess(score ~ pseudotime, data = data_dorc, span = 0.2)
    data_dorc$smooth <- DORC_fit$fitted
    RNA_fit <- loess(score ~ pseudotime, data = data_rna, span = 0.2)
    data_rna$smooth <- RNA_fit$fitted
    
    data <- rbind(data_dorc, data_rna)
    data$type <- factor(data$type, levels = c("DORC", "RNA"))
    
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
    print(p)
  }
  dev.off()
  
  
  
  # ---------------- trajectory Liver ---------------- #
  pdf("Liver.trajectory_dorc_gene_lineplot.pdf", width = 4, height = 3)
  cells <- proj$cellNames[!is.na(proj$Liver.trajectory)]
  pseudotime <- proj$Liver.trajectory[!is.na(proj$Liver.trajectory)]
  for(gene in c("Gzmb", "Mki67", "Ifrd1", "Tnfaip3")){
    data_dorc <- data.frame(score = dorcMat_sm_norm[gene, cells], pseudotime = pseudotime, type = "DORC")
    data_rna <- data.frame(score = rnaMat_sm_norm[gene, cells], pseudotime = pseudotime, type = "RNA")
    
    DORC_fit <- loess(score ~ pseudotime, data = data_dorc, span = 0.2)
    data_dorc$smooth <- DORC_fit$fitted
    RNA_fit <- loess(score ~ pseudotime, data = data_rna, span = 0.2)
    data_rna$smooth <- RNA_fit$fitted
    
    data <- rbind(data_dorc, data_rna)
    data$type <- factor(data$type, levels = c("DORC", "RNA"))
    
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
    print(p)
  }
  dev.off()
  
  # ---------------- trajectory Blood ---------------- #
  pdf("Blood.trajectory_dorc_gene_lineplot.pdf", width = 4, height = 3)
  cells <- proj$cellNames[!is.na(proj$Blood.trajectory)]
  pseudotime <- proj$Blood.trajectory[!is.na(proj$Blood.trajectory)]
  for(gene in c("Gzmb", "Mki67", "Ifrd1", "Tnfaip3")){
    data_dorc <- data.frame(score = dorcMat_sm_norm[gene, cells], pseudotime = pseudotime, type = "DORC")
    data_rna <- data.frame(score = rnaMat_sm_norm[gene, cells], pseudotime = pseudotime, type = "RNA")
    
    DORC_fit <- loess(score ~ pseudotime, data = data_dorc, span = 0.2)
    data_dorc$smooth <- DORC_fit$fitted
    RNA_fit <- loess(score ~ pseudotime, data = data_rna, span = 0.2)
    data_rna$smooth <- RNA_fit$fitted
    
    data <- rbind(data_dorc, data_rna)
    data$type <- factor(data$type, levels = c("DORC", "RNA"))
    
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
    print(p)
  }
  dev.off()
  
  
  
  # ---------------- trajectory Ln_re ---------------- #
  pdf("Ln.trajectory_re_dorc_gene_lineplot.pdf", width = 4, height = 3)
  cells <- proj$cellNames[!is.na(proj$Ln.trajectory_re)]
  pseudotime <- proj$Ln.trajectory_re[!is.na(proj$Ln.trajectory_re)]
  for(gene in c("Gzmb", "Mki67", "Ifrd1", "Tnfaip3")){
    data_dorc <- data.frame(score = dorcMat_sm_norm[gene, cells], pseudotime = pseudotime, type = "DORC")
    data_rna <- data.frame(score = rnaMat_sm_norm[gene, cells], pseudotime = pseudotime, type = "RNA")
    
    DORC_fit <- loess(score ~ pseudotime, data = data_dorc, span = 0.2)
    data_dorc$smooth <- DORC_fit$fitted
    RNA_fit <- loess(score ~ pseudotime, data = data_rna, span = 0.2)
    data_rna$smooth <- RNA_fit$fitted
    
    data <- rbind(data_dorc, data_rna)
    data$type <- factor(data$type, levels = c("DORC", "RNA"))
    
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
    print(p)
  }
  dev.off()
  
  
  # ---------------- trajectory Spleen_re ---------------- #
  pdf("Spleen.trajectory_re_dorc_gene_lineplot.pdf", width = 4, height = 3)
  cells <- proj$cellNames[!is.na(proj$Spleen.trajectory_re)]
  pseudotime <- proj$Spleen.trajectory_re[!is.na(proj$Spleen.trajectory_re)]
  for(gene in c("Gzmb", "Mki67", "Ifrd1", "Tnfaip3")){
    data_dorc <- data.frame(score = dorcMat_sm_norm[gene, cells], pseudotime = pseudotime, type = "DORC")
    data_rna <- data.frame(score = rnaMat_sm_norm[gene, cells], pseudotime = pseudotime, type = "RNA")
    
    DORC_fit <- loess(score ~ pseudotime, data = data_dorc, span = 0.2)
    data_dorc$smooth <- DORC_fit$fitted
    RNA_fit <- loess(score ~ pseudotime, data = data_rna, span = 0.2)
    data_rna$smooth <- RNA_fit$fitted
    
    data <- rbind(data_dorc, data_rna)
    data$type <- factor(data$type, levels = c("DORC", "RNA"))
    
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
    print(p)
  }
  dev.off()
  
  
  # ---------------- trajectory Liver_re ---------------- #
  pdf("Liver.trajectory_re_dorc_gene_lineplot.pdf", width = 4, height = 3)
  cells <- proj$cellNames[!is.na(proj$Liver.trajectory_re)]
  pseudotime <- proj$Liver.trajectory_re[!is.na(proj$Liver.trajectory_re)]
  for(gene in c("Gzmb", "Mki67", "Ifrd1", "Tnfaip3")){
    data_dorc <- data.frame(score = dorcMat_sm_norm[gene, cells], pseudotime = pseudotime, type = "DORC")
    data_rna <- data.frame(score = rnaMat_sm_norm[gene, cells], pseudotime = pseudotime, type = "RNA")
    
    DORC_fit <- loess(score ~ pseudotime, data = data_dorc, span = 0.2)
    data_dorc$smooth <- DORC_fit$fitted
    RNA_fit <- loess(score ~ pseudotime, data = data_rna, span = 0.2)
    data_rna$smooth <- RNA_fit$fitted
    
    data <- rbind(data_dorc, data_rna)
    data$type <- factor(data$type, levels = c("DORC", "RNA"))
    
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
    print(p)
  }
  dev.off()
  
  
  # ---------------- trajectory Blood_re ---------------- #
  pdf("Blood.trajectory_re_dorc_gene_lineplot.pdf", width = 4, height = 3)
  cells <- proj$cellNames[!is.na(proj$Blood.trajectory_re)]
  pseudotime <- proj$Blood.trajectory_re[!is.na(proj$Blood.trajectory_re)]
  for(gene in c("Gzmb", "Mki67", "Ifrd1", "Tnfaip3")){
    data_dorc <- data.frame(score = dorcMat_sm_norm[gene, cells], pseudotime = pseudotime, type = "DORC")
    data_rna <- data.frame(score = rnaMat_sm_norm[gene, cells], pseudotime = pseudotime, type = "RNA")
    
    DORC_fit <- loess(score ~ pseudotime, data = data_dorc, span = 0.2)
    data_dorc$smooth <- DORC_fit$fitted
    RNA_fit <- loess(score ~ pseudotime, data = data_rna, span = 0.2)
    data_rna$smooth <- RNA_fit$fitted
    
    data <- rbind(data_dorc, data_rna)
    data$type <- factor(data$type, levels = c("DORC", "RNA"))
    
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
    print(p)
  }
  dev.off()
}
