# figure 1
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
source("/home/zeemeeuw/YangLab/JoINT-seq/QC/code/utils.R")
setwd("/home/zeemeeuw/YangLab/JoINT-seq/Manuscript/fig_data/fig_1_data")


# Fig 1c: comparison of ATAC fragments
{
  joint_3t3 <- read.delim("G3MERGE-ata_mm10_1_fragments_count.txt", header = FALSE, stringsAsFactors = FALSE, col.names = c("barcode", "counts"))
  rownames(joint_3t3) <- joint_3t3$barcode
  whitelist_1 <- read.csv("G3-ata_mm10_1_bc_whitelist.txt")$x
  joint_3t3 <- joint_3t3[joint_3t3$barcode %in% whitelist_1,]
  joint_3t3$type <- "JoINT-seq 37"
  joint_3t3$assay <- "JoINT-seq"
  # joint_3t3 <- joint_3t3[joint_3t3$counts >= 3000,]
  joint_3t3$x <- 1
  
  
  joint_3t32 <- read.delim("G3MERGE-ata_mm10_2_fragments_count.txt", header = FALSE, stringsAsFactors = FALSE, col.names = c("barcode", "counts"))
  rownames(joint_3t32) <- joint_3t32$barcode
  whitelist_2 <- read.csv("G3-ata_mm10_2_bc_whitelist.txt")$x
  joint_3t32 <- joint_3t32[joint_3t32$barcode %in% whitelist_2,]
  joint_3t32$type <- "JoINT-seq 50"
  joint_3t32$assay <- "JoINT-seq"
  # joint_3t32 <- joint_3t32[joint_3t32$counts >= 3000,]
  joint_3t32$x <- 2
  
  multiome_GM <- read.delim("multiome-GM12878-atac.csv", header = TRUE, sep = ",", col.names = c("n", "project_cellnames", "counts", "barcode"))
  multiome_GM <- multiome_GM[,c("barcode", "counts")]
  multiome_GM$type <- "10X multiome"
  multiome_GM$assay <- "10X multiome"
  multiome_GM$x <- 3
  
  
  share_GM2 <- read.delim("share-GM12878rep1.fragments.count.txt", header = FALSE, col.names = c("barcode", "counts"))
  # share_GM2$type <- "GM12878_rep1"
  share_GM2$type <- "SHARE-seq"
  share_GM2$assay <- "SHARE-seq"
  share_GM2 <- share_GM2[share_GM2$counts >= 3000,]
  share_GM2$x <- 4
  
  paired_293T <- read.csv("paire-293T-fragment.count.csv", header = TRUE, row.names = 1, col.names = c("barcode", "counts", "identity"))
  paired_293T <- paired_293T[,c(1,2)]
  paired_293T$type <- "Paired-seq"
  paired_293T$assay <- "Paired-seq"
  paired_293T <- paired_293T[paired_293T$counts >= 500,]
  paired_293T$x <- 5
  
  sciCAR_293T <- read.csv("sciCAR-293T-fragment.count.csv", header = TRUE, row.names = 1, col.names = c("barcode", "counts", "identity"))
  sciCAR_293T <- sciCAR_293T[,c(1,2)]
  sciCAR_293T$type <- "sci-CAR-seq"
  sciCAR_293T$assay <- "sci-CAR-seq"
  sciCAR_293T <- sciCAR_293T[sciCAR_293T$counts >= 100,]
  sciCAR_293T$x <- 6
  
  SNARE_GM <- read.csv("SNARE-GM12878-fragment.count.csv", header = TRUE, row.names = 1, col.names = c("barcode", "counts", "identity"))
  SNARE_GM <- SNARE_GM[,c(1,2)]
  SNARE_GM$type <- "SNARE-seq"
  SNARE_GM$assay <- "SNARE-seq"
  SNARE_GM$x <- 7
  
  data <- rbind(joint_3t3, joint_3t32, multiome_GM, share_GM2, paired_293T, sciCAR_293T, SNARE_GM)
  data$type <- factor(data$type, levels = c("JoINT-seq 37", "JoINT-seq 50", "10X multiome", "SHARE-seq", "Paired-seq", "sci-CAR-seq", "SNARE-seq"))
  data$assay <- factor(data$assay, levels = c("JoINT-seq", "10X multiome", "SHARE-seq", "Paired-seq", "sci-CAR-seq", "SNARE-seq"))
  data$y <- log10(data$counts)
  
  data$type <- factor(data$type, levels = c("JoINT-seq 37", "JoINT-seq 50", "10X multiome", "SHARE-seq", "Paired-seq", "sci-CAR-seq", "SNARE-seq"))
  data$assay <- factor(data$assay, levels = c("JoINT-seq", "10X multiome", "SHARE-seq", "Paired-seq", "sci-CAR-seq", "SNARE-seq"))
  
  write_tsv(data, "source_fig_1c.tsv")
  
  ggplot(data = data, aes(x = x, y = y)) +
    theme_classic() +
    geom_boxplot(aes(group = x), width = 0.3, alpha = 0.7, show.legend = TRUE, outlier.shape = NA, notch = FALSE, fill = "lightGrey", size = 0.3) +
    geom_violin(aes(color = type), adjust = 1, trim = TRUE, width = 0.6, fill = NA, size = 0.5, scale = "width") +
    theme(axis.text.x = element_text(size = 10, angle = 30, hjust = 1, color = c("red", "red", "black", "black", "black", "black", "black")),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.title.y = element_text(size = 10)) +
    guides(fill = "none", color = "none") +
    scale_color_manual(values = c('#364070', colorRampPalette(c('#364070', '#3A64BC', '#4e8fcf', '#2ed978'))(6))) +
    scale_x_continuous(breaks = c(1,2,3,4,5,6,7), labels = c("JoINT-seq 37", "JoINT-seq 50", "10X multiome", "SHARE-seq", "Paired-seq", "sci-CAR-seq", "SNARE-seq")) +
    scale_y_continuous(limits = c(2, 6.1)) +
    xlab("") +
    ylab(expression(log[10]~(fragments~per~cell))) +
    annotate(geom = "text", x = 1, y = 6, colour = "red", size = 3.5, angle = 30, label = round(median(joint_3t3$counts))) +
    annotate(geom = "text", x = 2, y = 6, colour = "red", size = 3.5, angle = 30, label = round(median(joint_3t32$counts))) +
    annotate(geom = "text", x = 3, y = 6, colour = "black", size = 3.5, angle = 30, label = round(median(multiome_GM$counts))) +
    annotate(geom = "text", x = 4, y = 6, colour = "black", size = 3.5, angle = 30, label = round(median(share_GM2$counts))) +
    annotate(geom = "text", x = 5, y = 6, colour = "black", size = 3.5, angle = 30, label = round(median(paired_293T$counts))) +
    annotate(geom = "text", x = 6, y = 6, colour = "black", size = 3.5, angle = 30, label = round(median(sciCAR_293T$counts))) +
    annotate(geom = "text", x = 7, y = 6, colour = "black", size = 3.5, angle = 30, label = round(median(SNARE_GM$counts)))
}

# Fig 1d: comparison of umis
{
  joint_3t3_umis <- read.delim("G3MERGE-rna_mm10_1_umis_count.txt", header = FALSE, stringsAsFactors = FALSE, col.names = c("barcode", "counts"))
  joint_3t3_genes <- read.delim("G3MERGE-rna_mm10_1_genes_count.txt", header = FALSE, stringsAsFactors = FALSE, col.names = c("barcode", "genes"))
  joint_3t3 <- merge(joint_3t3_umis, joint_3t3_genes, by = "barcode")
  rownames(joint_3t3) <- joint_3t3$barcode
  whitelist_1 <- read.csv("G3-rna_mm10_1_bc_whitelist.txt")$x
  joint_3t3 <- joint_3t3[joint_3t3$barcode %in% whitelist_1, ]
  joint_3t3$type <- "JoINT-seq 37"
  joint_3t3$assay <- "JoINT-seq"
  # joint_3t3 <- joint_3t3[joint_3t3$counts >= 3000,] # 2400 for unstranded featureCounts
  joint_3t3$x <- 1
  
  joint_3t3_umis2 <- read.delim("G3MERGE-rna_mm10_2_umis_count.txt", header = FALSE, stringsAsFactors = FALSE, col.names = c("barcode", "counts"))
  joint_3t3_genes2 <- read.delim("G3MERGE-rna_mm10_2_genes_count.txt", header = FALSE, stringsAsFactors = FALSE, col.names = c("barcode", "genes"))
  joint_3t32 <- merge(joint_3t3_umis2, joint_3t3_genes2, by = "barcode")
  rownames(joint_3t32) <- joint_3t32$barcode
  whitelist_2 <- read.csv("G3-rna_mm10_2_bc_whitelist.txt")$x
  joint_3t32 <- joint_3t32[joint_3t32$barcode %in% whitelist_2, ]
  joint_3t32$type <- "JoINT-seq 50"
  joint_3t32$assay <- "JoINT-seq"
  # joint_3t3 <- joint_3t3[joint_3t3$counts >= 3000,] # 2400 for unstranded featureCounts
  joint_3t32$x <- 2
  
  multiome_GM <- read.delim("multiome-GM12878-RNA.csv", header = TRUE, sep = ",")
  multiome_GM <- multiome_GM[c(2,4,6,8,10)]
  colnames(multiome_GM) <- c( "barcode", "project_cellnames", "counts", "genes", "percent.mt")
  multiome_GM <- multiome_GM[,c("barcode", "counts", "genes")]
  multiome_GM$type <- "10X multiome"
  multiome_GM$assay <- "10X multiome"
  multiome_GM$x <- 3

  share_GM2 <- read.csv("share-GM12878rep1-umi-gene-cell.count.csv", header = TRUE, col.names = c("barcode", "counts", "genes"))
  # share_GM2$type <- "GM12878_rep1"
  share_GM2$type <- "SHARE-seq"
  share_GM2$assay <- "SHARE-seq"
  share_GM2 <- share_GM2[share_GM2$counts >= 3000,]
  share_GM2$x <- 4
  
  paired_293T <- read.csv("paire-293T-umi-gene.count.csv", header = TRUE, row.names = 1, col.names = c("barcode", "counts", "genes", "identity"))
  paired_293T <- paired_293T[,c(1,2,3)]
  paired_293T$type <- "Paired-seq"
  paired_293T$assay <- "Paired-seq"
  paired_293T <- paired_293T[paired_293T$counts >= 300,]
  paired_293T$x <- 5
  
  sciCAR_293T <- read.csv("sciCAR-293T-umi-gene-cell.count.csv", header = TRUE, col.names = c("barcode", "counts", "genes", "identity"))
  sciCAR_293T <- sciCAR_293T[,c(1,2,3)]
  sciCAR_293T$type <- "sci-CAR-seq"
  sciCAR_293T$assay <- "sci-CAR-seq"
  sciCAR_293T <- sciCAR_293T[sciCAR_293T$counts >= 800,]
  sciCAR_293T$x <- 6
  
  SNARE_GM <- read.csv("SNARE-GM12878-umi-gene.count.csv", header = TRUE, col.names = c("barcode", "counts", "genes", "identity"))
  SNARE_GM <- SNARE_GM[,c(1,2,3)]
  SNARE_GM$type <- "SNARE-seq"
  SNARE_GM$assay <- "SNARE-seq"
  SNARE_GM$x <- 7
  
  data <- rbind(joint_3t3, joint_3t32, multiome_GM, share_GM2, paired_293T, sciCAR_293T, SNARE_GM)
  data$type <- factor(data$type, levels = c("JoINT-seq 37", "JoINT-seq 50", "10X multiome", "ISSAAC-seq", "SHARE-seq", "Paired-seq", "sci-CAR-seq", "SNARE-seq"))
  data$assay <- factor(data$assay, levels = c("JoINT-seq", "10X multiome", "SHARE-seq", "Paired-seq", "sci-CAR-seq", "SNARE-seq"))
  data$y <- log10(data$counts)
  
  data$type <- factor(data$type, levels = c("JoINT-seq 37", "JoINT-seq 50", "10X multiome", "ISSAAC-seq", "SHARE-seq", "Paired-seq", "sci-CAR-seq", "SNARE-seq"))
  data$assay <- factor(data$assay, levels = c("JoINT-seq", "10X multiome", "SHARE-seq", "Paired-seq", "sci-CAR-seq", "SNARE-seq"))
  
  write_tsv(data, "source_fig_1d.tsv")
  
  ggplot(data = data, aes(x = x, y = y)) +
    theme_classic() +
    geom_boxplot(aes(group = x), width = 0.3, alpha = 0.7, show.legend = TRUE, outlier.shape = NA, notch = FALSE, fill = "lightGrey", size = 0.3) +
    geom_violin(aes(color = type), adjust = 1, trim = TRUE, width = 0.6, fill = NA, size = 0.5, scale = "width") +
    ylim(2, 5.1) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 10, color = c("red", "red", "black", "black", "black", "black", "black")),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.title.y = element_text(size = 10)) +
    guides(fill = "none", color = "none") +
    scale_color_manual(values = c('#800505', '#800505', "#cf0e0e", '#e04b19', '#f08b3e', '#f0ac3e', '#f0d83e')) +
    scale_x_continuous(breaks = c(1,2,3,4,5,6,7), labels = c("JoINT-seq 37", "JoINT-seq 50", "10X multiome", "SHARE-seq", "Paired-seq", "sci-CAR-seq", "SNARE-seq")) +
    xlab("") +
    ylab(expression(log[10]~(umis~per~cell))) +
    ylim(2, 5.45) +
    annotate(geom = "text", x = 1, y = 5.4, colour = "red", size = 3.5, angle = 30, label = round(median(joint_3t3$counts))) +
    annotate(geom = "text", x = 2, y = 5.4, colour = "red", size = 3.5, angle = 30, label = round(median(joint_3t32$counts))) +
    annotate(geom = "text", x = 3, y = 5.4, colour = "black", size = 3.5, angle = 30, label = round(median(multiome_GM$counts))) +
    annotate(geom = "text", x = 4, y = 5.4, colour = "black", size = 3.5, angle = 30, label = round(median(share_GM2$counts))) +
    annotate(geom = "text", x = 5, y = 5.4, colour = "black", size = 3.5, angle = 30, label = round(median(paired_293T$counts))) +
    annotate(geom = "text", x = 6, y = 5.4, colour = "black", size = 3.5, angle = 30, label = round(median(sciCAR_293T$counts))) +
    annotate(geom = "text", x = 7, y = 5.4, colour = "black", size = 3.5, angle = 30, label = round(median(SNARE_GM$counts)))
}


# Fig 1e: TPM of full-length TCR/BCR
{
  imgt_colors <- c(FR1 = "#1f6933", 
                   CDR1 = "#66c992", 
                   FR2 = "#248AF3", 
                   CDR2 = "#88CEEF", 
                   FR3 = "#CF3D00", 
                   CDR3 = "#FF851D", 
                   FR4 = "#FFC782")
  
  pdata_tra <- read_delim('GJx-rC1_alignments_TRA_noumi_coverage_inner.tsv')
  pdata_trb <- read_delim('GJx-rC1_alignments_TRB_noumi_coverage_inner.tsv')
  pdata_igh <- read_delim('GJx-B1_alignments_IGH_noumi_coverage_inner.tsv')
  pdata_igl <- read_delim('GJx-B1_alignments_IGL_noumi_coverage_inner.tsv')
  
  pdata <-bind_rows(pdata_tra, pdata_trb, pdata_igh, pdata_igl) %>% 
    dplyr::mutate_at(vars(group), ~factor(.x, levels = c('FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4'))) %>% 
    dplyr::mutate_at(vars(chain), ~factor(.x, levels = c('TRA', 'TRB', 'IGH', 'IGL'))) %>% 
    arrange(percentile)
  
  ggplot() +
    geom_line(data = pdata_tra, aes(color = group, x = percentile, y = y), linewidth = 0.8, lineend = "round") +
    geom_line(data = pdata_trb, aes(color = group, x = percentile, y = y), linewidth = 0.8, lineend = "round") +
    geom_line(data = pdata_igh, aes(color = group, x = percentile, y = y), linewidth = 0.8, lineend = "round") +
    geom_line(data = pdata_igl, aes(color = group, x = percentile, y = y), linewidth = 0.8, lineend = "round") +
    ylim(c(0, 7)) +
    theme_bw() +
    scale_color_manual(values = imgt_colors) +
    ylab('log10(reads)') +
    xlab('Position percentile (%)') +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(size = 10, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.title.x = element_text(size = 13),
          axis.title.y = element_text(size = 13))
  
  
  
  
}

# Fig 1f: Detection rate of BCR and TCR
{

   
  ### TCR
  meta <- read_tsv("source_fig_1f_TCR.tsv") %>% 
    dplyr::mutate_at(vars(group), ~factor(., levels = c("undetected", "tra_only", "trb_only", "paired")))
  ggplot(meta, aes(x = sample)) +
    geom_col(aes(fill = group, y = prop), position = "fill") +
    theme_bw() +
    scale_fill_manual(values = c("undetected" = "grey90", "tra_only" = "#FCD48F", "trb_only" = "#FF9F01", "paired" = "#EE3B00")) +
    theme(axis.text.x = element_text(size = 12, color = "black"),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank()) +
    annotate("text", size = 4, color = "white", fontface = "bold", x = 1, y = meta$prop[meta$group == "paired"] - 0.03, label = paste0(round(100 * meta$prop[meta$group == "paired"], digits = 0), "%" )) +
    annotate("text", size = 4, color = "black", fontface = "bold", x = 1, y = sum(meta$prop[meta$group != "undetected"]) + 0.03, label = paste0(round(100 * sum(meta$prop[meta$group != "undetected"]) + 0.03, digits = 0), "%" )) +
    coord_fixed(2.5)
  
  ### BCR
  meta <- read_tsv("source_fig_1f_BCR.tsv") %>% 
    dplyr::mutate_at(vars(group), ~factor(., levels = c("undetected", "igl_only", "igh_only", "paired")))
  ggplot(meta, aes(x = sample)) +
    geom_col(aes(fill = group, y = prop), position = "fill") +
    theme_bw() +
    scale_fill_manual(values = c("undetected" = "grey90", "igl_only" = "#FCD48F", "igh_only" = "#FF9F01", "paired" = "#EE3B00")) +
    theme(axis.text.x = element_text(size = 12, color = "black"),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank()) +
    annotate("text", size = 4, color = "white", fontface = "bold", x = 1, y = meta$prop[meta$group == "paired"] - 0.03, label = paste0(round(100 * meta$prop[meta$group == "paired"], digits = 0), "%" )) +
    annotate("text", size = 4, color = "black", fontface = "bold", x = 1, y = sum(meta$prop[meta$group != "undetected"]) + 0.03, label = paste0(round(100 * sum(meta$prop[meta$group != "undetected"]) + 0.03, digits = 0), "%" )) +
    coord_fixed(2.5)
}


# Fig 1h: PBMC clustering
{

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
  
  data <- read_tsv("source_fig_1h.tsv") %>% 
    dplyr::mutate_at(vars(cell_module), ~factor(., levels = names(hs_colors)))
  ggplot(data, aes(x = UMAPLSICombined_1, y = UMAPLSICombined_2, color = cell_module)) +
    theme_classic() +
    geom_point(size = 0.1) +
    scale_color_manual(values = hs_colors) +
    guides(color = guide_legend(override.aes = list(size = 3)))
    
  
}

# Fig 1i: TCR projection
{

  # project TCR
  library(ggplot2)
  library(readxl)

  data2 <- read_tsv("source_fig_1i.tsv")
  data_tcr <- data2 %>% dplyr::filter(!is.na(n_clonotype)) %>% 
    dplyr::filter(!str_detect(cell_module, "_B_")) %>% 
    dplyr::filter(!str_detect(cell_module, "_Mono_"))
  xrange <- max(data2$UMAPLSICombined_1) - min(data2$UMAPLSICombined_1)
  yrange <- max(data2$UMAPLSICombined_2) - min(data2$UMAPLSICombined_2)
  ggplot(data2) +
    theme_void() +
    coord_fixed(xrange/yrange) +
    geom_point(data = data_tcr, aes(x = UMAPLSICombined_1, y = UMAPLSICombined_2, size = n_clonotype, color = cell_module), alpha = 1) +
    scale_size_continuous(range = c(0.1, 4), breaks = c(1, 5, 10)) +
    scale_color_manual(values = hs_colors)
  
  # MAIT TCR composition
  data <- read_excel("Fig1-TCRm457-clonotypes.xlsx") %>% 
    group_by(CDRa_CDRb, barcode) %>% 
    distinct(CDRa_CDRb, barcode, .keep_all = TRUE) %>% 
    dplyr::filter(!is.na(chainA) & !is.na(chainB)) %>% 
    ungroup() %>% 
    dplyr::mutate(TRAV = str_extract(chainA, "(.*)\\+(.*)\\+(.*)", group = 1),
                  TRAJ = str_extract(chainA, "(.*)\\+(.*)\\+(.*)", group = 2),
                  CDR3a = str_extract(chainA, "(.*)\\+(.*)\\+(.*)", group = 3),
                  TRBV = str_extract(chainB, "(.*)\\+(.*)\\+(.*)", group = 1),
                  TRBJ = str_extract(chainB, "(.*)\\+(.*)\\+(.*)", group = 2),
                  CDR3b = str_extract(chainB, "(.*)\\+(.*)\\+(.*)", group = 3),
                  comb_a = paste0(TRAV, "_", TRAJ),
                  comb_b = paste0(TRBV, "_", TRBJ))
  
  data_sum <- data %>% left_join(data2[,c("cellnames", "cell_module")], by = join_by("barcode" == "cellnames")) %>% 
    dplyr::filter(!is.na(cell_module)) %>% 
    dplyr::mutate(MAIT = if_else(cell_module == "C14_MAIT_SLC4A10_RORA", "MAIT", "others")) %>% 
    group_by(comb_a, MAIT) %>% 
    dplyr::summarise(count_comb_a = n()) %>% 
    arrange(desc(count_comb_a)) %>% 
    group_by(MAIT) %>% 
    dplyr::mutate(prop_comb_a = count_comb_a/sum(count_comb_a)) %>% 
    ungroup() %>% 
    arrange(desc(prop_comb_a))
  comb_a_plot <- data_sum %>% slice_head(n = 40) %>% pull("comb_a") %>% unique()
  data_sum <- data_sum %>% dplyr::filter(comb_a %in% comb_a_plot) %>% 
    dplyr::mutate_at(vars(comb_a), ~factor(., levels = unique(comb_a)))
  data_sum2 <- full_join(dplyr::filter(data_sum, MAIT == "MAIT"), dplyr::filter(data_sum, MAIT != "MAIT"), by = "comb_a")
  df1 <- data_sum2 %>% dplyr::select(comb_a, ends_with("x")) %>% 
    rename_with(~c("comb_a", "MAIT", "count_comb_a", "prop_comb_a")) %>% 
    dplyr::mutate(MAIT = "MAIT")
  df2 <- data_sum2 %>% dplyr::select(comb_a, ends_with("y")) %>% 
    rename_with(~c("comb_a", "MAIT", "count_comb_a", "prop_comb_a")) %>% 
    dplyr::mutate(MAIT = "others") %>% 
    dplyr::mutate(prop_comb_a = -prop_comb_a)
  data_sum2 <- bind_rows(df1, df2) %>% 
    replace(is.na(.), 0) 
  write_tsv(data_sum2, "source_fig_1i_MAIT.tsv", na = "")
  ggplot(data_sum2, aes(x = comb_a, y = prop_comb_a, fill = MAIT)) +
    theme_classic() +
    geom_col() +
    theme(axis.text.x = element_text(angle = 45, size = 5, hjust = 1, color = "black"),
          axis.text.y = element_text(size = 5, color = "black"),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(t = 5, l = 40)) +
    ylim(-0.1, 0.45) +
    scale_fill_manual(values = c("#DF4A4A", "#5668C3"))
}

# Fig 1J: TCR correlation with 10X
{

  
  pdata <- read_tsv("source_fig_1J.tsv")
  fit <- summary(lm(n.y~n.x, data = pdata))
  sqrt(fit$r.squared)
  ratio <- (max(pdata$n.x)-min(pdata$n.x))/(max(pdata$n.y)-min(pdata$n.y))
  ggplot(pdata, aes(x = n.x, y = n.y)) +
    geom_point(aes(fill = CDRa_CDRb), show.legend = FALSE, size = 5, alpha = 0.7, shape = 21, color = "black") +
    geom_line(stat="smooth", method = "lm", formula = y ~ x,
              size = 1, linetype ="dashed", alpha = 0.5, color = "blue") +
    theme_classic() +
    ylab("joint") +
    xlab("10X") +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 10, color = c("black", "black")),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.title.x = element_text(size = 13),
          axis.title.y = element_text(size = 13),
          legend.text = element_text(size = 10)) +
    coord_equal(ratio) +
    scale_fill_manual(values = choose_colorset("bubble", 72))
}

# Fig 1K: Master regulators of TCR clones
{

  
  pdata <- read_tsv("source_fig_1K.tsv")
  ggplot(pdata, aes(x = cluster, y = mean_diff, color = fill, group = cluster)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 8, color = "black"),
          axis.text.y = element_text(hjust = 0.5, vjust = 0.5, size = 6, color = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 8, color = "black"),
          panel.grid = element_blank()) +
    geom_jitter(aes(size = `-log10p`), alpha = 0.8, position = "jitter") +
    geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), linetype = 2, color = '#CAC3F2') +
    scale_size_continuous(range = c(1, 3)) +
    scale_color_gradientn(colors = colorRampPalette(colors = c("#14B3FF", "#88CEEF", "#C1D5DC", "#EAD397", "#FDB31A", "#E42A2A"))(256), limits = c(-0.1, 0.1)) +  
    ggrepel::geom_text_repel(aes(label = Motif), size = 2, color = "black", fontface = "italic", vjust = 0.5, hjust = 0.5, max.overlaps = 100, direction = "both",
                             force = 1) +
    ylab("Mean difference compared to background") +
    labs(color = 'Mean\ndifference')
}



# Fig 1M: FRiP and umis of mouse WT libraries
{
  
  pdata <- read_tsv("source_fig_1M_frip.tsv")
  pdata$tissue <- factor(pdata$tissue, levels = c("Ln", "Spleen", "Liver", "Blood"))
  ggplot(pdata, aes(x = tissue)) +   
    geom_violin(aes(fill = tissue, y = frip), alpha = 0.5, trim = TRUE, scale = "width") +
    geom_boxplot(aes(fill = tissue, y = frip), alpha = 1, show.legend = FALSE, width = 0.3, outlier.shape = NA) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, size = 12, hjust = 1, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.ticks.x = element_blank()) +
    scale_fill_manual(values = c("Spleen" = "#36A23D", "Ln" = "#2884E7", "Liver" = "#817ab9", "Blood" = "#EE3B00")) +
    ggtitle("FRiP")


  pdata <- read_tsv("source_fig_1M_umi.tsv")
  pdata$tissue <- factor(pdata$tissue, levels = c("Ln", "Spleen", "Liver", "Blood"))
  ggplot(pdata, aes(x = tissue)) +   
    geom_violin(aes(fill = tissue, y = log10(umis)), alpha = 0.5, trim = TRUE, scale = "width") +
    geom_boxplot(aes(fill = tissue, y = log10(umis)), alpha = 1, show.legend = FALSE, width = 0.3, outlier.shape = NA) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, size = 12, hjust = 1, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.ticks.x = element_blank()) +
    ylim(c(1, 5)) +
    scale_fill_manual(values = c("Spleen" = "#36A23D", "Ln" = "#2884E7", "Liver" = "#817ab9", "Blood" = "#EE3B00")) +
    ggtitle("UMI")
}

# Fig 1N: Detection rate of mouse TCRs
{

  
  meta <- read_tsv("source_fig_1N.tsv")
  meta$tissue <- factor(meta$tissue, levels = c("Ln", "Spleen", "Liver", "Blood"))
  meta$group <- factor(meta$group, levels = c("undetected", "alpha_only", "beta_only", "paired"))
  meta_ln <- meta %>% dplyr::filter(tissue == "Ln")
  meta_sp <- meta %>% dplyr::filter(tissue == "Spleen") 
  meta_liver <- meta %>% dplyr::filter(tissue == "Liver")
  meta_blood <- meta %>% dplyr::filter(tissue == "Blood")
  ggplot(meta, aes(x = tissue)) +
    geom_col(aes(fill = group, y = prop), position = "fill") +
    theme_bw() +
    scale_fill_manual(values = c("undetected" = "grey90", "alpha_only" = "#FCD48F", "beta_only" = "#FF9F01", "paired" = "#EE3B00")) +
    theme(axis.text.x = element_text(size = 12, color = "black"),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank()) +
    annotate("text", size = 4, color = "white", fontface = "bold", x = 1, y = meta_ln$prop[meta_ln$group == "paired"] - 0.03, label = paste0(round(100 * meta_ln$prop[meta_ln$group == "paired"], digits = 0), "%" )) +
    annotate("text", size = 4, color = "white", fontface = "bold", x = 2, y = meta_sp$prop[meta_sp$group == "paired"]- 0.03, label = paste0(round(100 * meta_sp$prop[meta_sp$group == "paired"], digits = 0), "%" )) +
    annotate("text", size = 4, color = "white", fontface = "bold", x = 3, y = meta_liver$prop[meta_liver$group == "paired"]- 0.03, label = paste0(round(100 * meta_liver$prop[meta_liver$group == "paired"], digits = 0), "%" )) +
    annotate("text", size = 4, color = "white", fontface = "bold", x = 4, y = meta_blood$prop[meta_blood$group == "paired"]- 0.03, label = paste0(round(100 * meta_blood$prop[meta_blood$group == "paired"], digits = 0), "%" )) +
    annotate("text", size = 4, color = "black", fontface = "bold", x = 1, y = sum(meta_ln$prop[meta_ln$group != "undetected"]) + 0.03, label = paste0(round(100 * sum(meta_ln$prop[meta_ln$group != "undetected"]) + 0.03, digits = 0), "%" )) +
    annotate("text", size = 4, color = "black", fontface = "bold", x = 2, y = sum(meta_sp$prop[meta_sp$group != "undetected"]) + 0.03, label = paste0(round(100 * sum(meta_sp$prop[meta_sp$group != "undetected"]) + 0.03, digits = 0), "%" )) +
    annotate("text", size = 4, color = "black", fontface = "bold", x = 3, y = sum(meta_liver$prop[meta_liver$group != "undetected"]) + 0.03, label = paste0(round(100 * sum(meta_liver$prop[meta_liver$group != "undetected"]) + 0.03, digits = 0), "%" )) +
    annotate("text", size = 4, color = "black", fontface = "bold", x = 4, y = sum(meta_blood$prop[meta_blood$group != "undetected"]) + 0.03, label = paste0(round(100 * sum(meta_blood$prop[meta_blood$group != "undetected"]) + 0.03, digits = 0), "%" ))
  
}

# Fig 1O: UMAP of mouse data
{

  wt_colors <- c("C1_Il7r_Tcf7" = "#bbdefb", 
                 "C2_Klra3_Tcf7" = "#3361A5",
                 "C3_Klrk1_Eomes" = "#a9df91", 
                 "C4_Mki67_Mga" = "#48b352", 
                 "C5_Havcr2_Nfatc1" = "#ff882e", 
                 "C6_Myb_Smarcc1" = "#e53a46")

  
  data <- read_tsv("source_fig_1O.tsv")
  data_tcr <- data %>% dplyr::filter(!is.na(n_clonotype))
  xrange <- max(data$UMAPLSICombined_1) - min(data$UMAPLSICombined_1)
  yrange <- max(data$UMAPLSICombined_2) - min(data$UMAPLSICombined_2)
  ggplot(data) +
    theme_void() +
    coord_fixed(xrange/yrange) +
    geom_point(data = data_tcr, aes(x = UMAPLSICombined_1, y = UMAPLSICombined_2, size = n_clonotype, color = cell_module), alpha = 1) +
    scale_size_continuous(range = c(0.1, 4), breaks = c(1, 5, 10)) +
    scale_color_manual(values = wt_colors)
  
}

# Fig 1P: TF enrichment of the largest TCR clone
{

  
  pdata <- read_tsv("source_fig_1P.tsv")
  xrange <- max(abs(pdata$mean_diff))*2
  yrange <- max(pdata$`-log10p`) - 0
  ggplot(pdata, aes(x = mean_diff, y = `-log10p`, color = mean_diff)) +
    geom_point(aes(color = fill), alpha = 1, size = 2, show.legend = TRUE) +
    geom_point(data = pdata %>% dplyr::filter(top10), aes(fill = fill), color = "black", alpha = 1, size = 4, show.legend = FALSE, shape = 21) +
    scale_fill_gradientn(colors = colorRampPalette(c("#66c992", "#86d9ac", "#d8dcc1", "#FFC782", "#FFA12F"))(100), limits = c(-0.3, 0.3)) +
    scale_color_gradientn(colors = colorRampPalette(c("#66c992", "#86d9ac", "#d8dcc1", "#FFC782", "#FFA12F"))(100), , limits = c(-0.3, 0.3)) +
    # scale_color_manual(values = c("Upregulated" = unname(mmjoint_colors[c1]), "Downregulated" = unname(mmjoint_colors[c2]), "Not significant" = "#f0f0f0")) +
    ggrepel::geom_text_repel(data = dplyr::filter(pdata, top10), aes(label = name), size = 3, color = "black", vjust = 0.5, hjust = 0.5, max.overlaps = 13) +
    # geom_vline(xintercept = c(-fc_line, fc_line), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    labs(x = "Mean Difference",
         y = "-log10(p)",
         color = "mean_diff") +
    theme_bw() +
    xlim(-max(abs(pdata$mean_diff)), max(abs(pdata$mean_diff))) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(color = "black")) +
    coord_fixed(ratio = xrange/yrange) +
    ggtitle("Others <-  -> C1")
}






