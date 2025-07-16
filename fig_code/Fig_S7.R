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

setwd("/home/zeemeeuw/YangLab/JoINT-seq/Manuscript/fig_data/fig_S7_data")

# S7A
{
  data_longer <- readxl::read_excel("Lm-TCR4.57-clonotypes-paired in rna.xlsx")
  
  data <- data_longer %>% 
    dplyr::mutate(tissue = str_extract(Tissue_time, "(.*)-", group = 1)) %>% 
    group_by(CDRa_CDRb, timepoint, tissue, barcode) %>% 
    distinct(CDRa_CDRb, timepoint, tissue, barcode, .keep_all = TRUE) %>% 
    dplyr::filter(!str_detect(CDR3a, "_") & !str_detect(CDR3b, "_")) %>% 
    ungroup()
  
  data_sum <- data %>% left_join(seu@meta.data[,c("cellnames", "cell_module")], by = join_by("barcode" == "cellnames")) %>% 
    group_by(CDRa_CDRb, cell_module) %>% 
    dplyr::mutate(n = n()) %>% 
    ungroup() %>% 
    arrange(desc(n)) %>% 
    dplyr::select(barcode, n) %>% 
    group_by(barcode) %>% 
    slice_head(n = 1)
  
  dft <- seu@meta.data %>% left_join(data_sum, by = join_by(cellnames == barcode))
  
  seu$n_clonotype <- dft$n
  
  data <- data.frame(seu@reductions$UMAP_LSI_Combined_batch3.3@cell.embeddings) 
  data$barcode <- rownames(data)
  data <- cbind(data, seu@meta.data)
  data_tcr <- data %>% dplyr::filter(!is.na(n_clonotype))
  
  ggplot(data_tcr) +
    geom_boxplot(aes(x = cell_module, fill = cell_module, y = log10(n_clonotype + 1)), show.legend = FALSE) +
    scale_fill_manual(values = mmjoint_colors) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
          axis.text.y = element_text(size = 10, color = "black")) 
  
}

# S7B
{
  seu$group <- paste0(seu$timepoint, "_", seu$tissue)
  data <- as.data.frame.array(table(seu$cell_module, seu$group))
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
                     # filename = "All_Ro_Re.pdf",
                     width = 8, height = 4,
                     scale = "none",
                     angle_col = 45,
                     display_numbers = TRUE, number_color = "black", fontsize_number = 6,
                     color = choose_colorset("oranges", 256))
}

# S7C
{
  files1 <- c("Activation_1-v-Priming_1-markerTF.xlsx", "Effector_1-v-Activation_1-markerTF.xlsx", 
              "Effector_1-v-Effector_2-markerTF.xlsx", "Effector-memory_transition_2-v-Effector_1-markerTF.xlsx", 
              "Effector-memory_transition_2-v-Effector_2-markerTF.xlsx", "Effector-memory_transition_2-v-Effector-memory_transition_3-markerTF.xlsx", 
              "Memory_1-v-Effector-memory_transition_2-markerTF.xlsx", "Memory_1-v-Quiescent-markerTF.xlsx", 
              "Memory_2-v-Effector-memory_transition_2-markerTF.xlsx", "Memory_2-v-Quiescent-markerTF.xlsx", 
              "Priming_1-v-Quiescent-markerTF.xlsx")
  
  tf_highlight <- c("Tcf7", "Lef1", "Fli1", "Tcf7l2", "Tcf7l2", "Fos", "Junb", "Smarcc1", "Bach1", "Bach2", "Irf1", "Nfkb1", "Nfkb2", "Mga", 
                    "Eomes", "Pou3f1", "Irf2", "Ets2", "Ctcf", "Id2", "Id3", "Id4", "Zeb1", "Zeb2", "Mnt", "Snai2", "Prdm1")
  
  for(f in files1){
    outname <- paste0(str_extract(f, "(.*)\\.xlsx", group = 1), ".pdf")
    c1 <- str_extract(f, "(.*)-v-", group = 1)
    c2 <- str_extract(f, "-v-(.*)-markerTF", group = 1)
    data <- read_excel(f)
    
    fc_line <- 0
    df <- data %>%
      dplyr::mutate(
        significance = case_when(
          FDR < 0.05 & MeanDiff > fc_line ~ "Upregulated",
          FDR < 0.05 & MeanDiff < -fc_line ~ "Downregulated",
          .default = "Not significant"), 
        `-log10FDR` = -log10(FDR)) %>% 
      dplyr::mutate(label = str_extract(name, "(.*)_", group = 1)) %>% 
      arrange(FDR) %>%
      group_by(significance) %>% 
      dplyr::mutate(id = consecutive_id(name)) %>% 
      dplyr::mutate(top10 = (id <= 20 & significance %in% c("Upregulated", "Downregulated")) | 
                      ((label %in% tf_highlight | str_detect(label, "^Klf|^Etv")) & significance %in% c("Upregulated", "Downregulated"))) %>% 
      dplyr::mutate(fill = case_when(MeanDiff >= 0.1 ~ 0.1, MeanDiff <= -0.1 ~ -0.1, .default = MeanDiff))
    
    
    
    # 绘制火山图
    xrange <- max(abs(df$MeanDiff))*2
    yrange <- max(df$`-log10FDR`) - 0
    ggplot(df, aes(x = MeanDiff, y = `-log10FDR`, color = significance)) +
      geom_point(aes(color = fill), alpha = 1, size = 2, show.legend = TRUE) +
      geom_point(data = df %>% dplyr::filter(top10), aes(fill = fill), color = "black", alpha = 1, size = 4, show.legend = FALSE, shape = 21) +
      scale_fill_gradientn(colors = colorRampPalette(c("#66c992", "#86d9ac", "#d8dcc1", "#FFC782", "#FFA12F"))(100), , limits = c(-0.1, 0.1)) +
      ggtitle(str_extract(f, "(.*)\\.xlsx", group = 1)) +
      scale_color_gradientn(colors = colorRampPalette(c("#66c992", "#86d9ac", "#d8dcc1", "#FFC782", "#FFA12F"))(100), , limits = c(-0.1, 0.1)) +
      # scale_color_manual(values = c("Upregulated" = unname(mmjoint_colors[c1]), "Downregulated" = unname(mmjoint_colors[c2]), "Not significant" = "#f0f0f0")) +
      ggrepel::geom_text_repel(data = dplyr::filter(df, top10), aes(label = label), size = 3, color = "black", vjust = 0.5, hjust = 0.5, max.overlaps = 13) +
      geom_vline(xintercept = c(-fc_line, fc_line), linetype = "dashed", color = "black") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
      labs(x = "Mean Difference",
           y = "-log10(FDR)",
           color = "MeanDiff") +
      theme_bw() +
      xlim(-max(abs(df$MeanDiff)), max(abs(df$MeanDiff))) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(color = "black"),
            axis.text.y = element_text(color = "black")) +
      coord_fixed(ratio = xrange/yrange)
    ggsave(outname, width = 4, height = 4)
    
  }
}




# S7D
{
  load(paste0(RDSdir, "calculated_seurat_gene_module.Rdata"))
  pdata <- seu@meta.data %>% 
    dplyr::select(cell_module) %>% 
    bind_cols(module_score_cd8) %>% 
    group_by(cell_module) %>% 
    dplyr::summarise_all(~mean(.x)) %>% 
    tibble::column_to_rownames("cell_module")
  pdata <- t(pdata)
  
  pheatmap::pheatmap(pdata[!(rownames(pdata) %in% c("Stress response", "Anergy", "Exhaustion", "Senescence")),], 
                     border_color = NA,
                     show_colnames = TRUE,
                     show_rownames = TRUE,
                     cluster_rows = TRUE,
                     cluster_cols = FALSE,
                     legend_breaks = c(-1, 0, 1),
                     clustering_method = "ward.D",
                     scale = "row",
                     fontsize_row = 10,
                     fontsize_col = 10,
                     angle_col = 45,
                     # color = choose_colorset("wolfgang_extra", 101),
                     color = rev(choose_colorset("RdBu", 101))[5:95],
                     # filename = "cell_module_module_score_cd8_heatmap.pdf",
                     width = 6, height = 4.5)
}


# S7E
{
  load(paste0(RDSdir, "calculated_seurat_gene_module.Rdata"))
  seu$timepoint <- factor(seu$timepoint, levels = c("WT", "d3", "d5", "d8", "d14", "d30", "R5", "Y5"))
  pdata <- seu@meta.data %>% 
    dplyr::select(tissue, timepoint) %>% 
    bind_cols(module_score_cd8) %>% 
    dplyr::filter(tissue == "Spleen") %>% 
    group_by(timepoint) %>% 
    dplyr::select(-tissue) %>% 
    dplyr::summarise_all(~mean(.x)) %>% 
    tibble::column_to_rownames("timepoint")
  pdata <- t(pdata)
  
  pheatmap::pheatmap(pdata[!(rownames(pdata) %in% c("Stress response", "Anergy", "Exhaustion", "Senescence")),c("WT", "d3", "d5", "d8", "d14", "d30", "R5", "Y5")], 
                     border_color = NA,
                     show_colnames = TRUE,
                     show_rownames = TRUE,
                     cluster_rows = TRUE,
                     cluster_cols = FALSE,
                     legend_breaks = c(-1, 0, 1),
                     clustering_method = "ward.D",
                     scale = "row",
                     fontsize_row = 10,
                     fontsize_col = 10,
                     angle_col = 45,
                     # color = choose_colorset("wolfgang_extra", 101),
                     color = rev(choose_colorset("RdBu", 101))[5:95],
                     # filename = "Spleen_timepoint_module_score_cd8_heatmap.pdf",
                     width = 5.5, height = 4.5)
}



# S7G
{
  # get matrix
  gs_mat <- getMatrixFromProject(proj, "GeneScoreMatrix")
  ge_mat <- getMatrixFromProject(proj, "GeneExpressionMatrix")
  
  
  gene <- "Havcr2"
  mat_gs <- assay(gs_mat[which(rowData(gs_mat)$name == gene),])
  mat_ge <- assay(ge_mat[which(rowData(ge_mat)$name == gene),])
  
  plot_data_Havcr2 <- data.frame(
    cell_id = proj$cellNames,
    cell_type = proj$cell_module,
    GeneScore = calculated_smoothed_ge(weightList, mat = mat_gs)[proj$cellNames],
    GeneExpression = calculated_smoothed_ge(weightList, mat = mat_ge)[proj$cellNames]
  ) 
  rt <- mean(plot_data_Havcr2$GeneScore)/mean(plot_data_Havcr2$GeneExpression)
  # rt <- 1
  pdata <- plot_data_Havcr2 %>% dplyr::mutate(GeneExpression = GeneExpression*rt)
  pdata <- pdata %>%
    pivot_longer(cols = c("GeneScore", "GeneExpression"),
                 names_to = "data_type",
                 values_to = "value")
  pdata$cell_type <- factor(pdata$cell_type, levels = c("Quiescent", "Priming_1", "Priming_2", "Activation_1", "Activation_2", 
                                                        "Activation_3", "Effector_1", "Effector_2", "Effector-memory_transition_1", 
                                                        "Effector-memory_transition_2", "Effector-memory_transition_3", 
                                                        "Memory_1", "Memory_2", "Memory_3"))
  
  ggplot()+
    geom_half_violin(data = pdata %>% filter(data_type == "GeneScore"),
                     aes(x = cell_type, y = log1p(value), fill = cell_type, color = cell_type), adjust = 1.2, 
                     side = "l", scale = "width", width = 0.9, alpha = 0.8, show.legend = FALSE) +
    scale_fill_manual(values = mmjoint_colors) +
    labs(title = gene, 
         x = "cell module", 
         y = "Log1p(smoothed genescore)") +
    scale_y_continuous(limits = c(0, 1*log1p(max(pdata$value))),
                       sec.axis = sec_axis(~./(log1p(mean(plot_data_Havcr2$GeneScore))/log1p(mean(plot_data_Havcr2$GeneExpression))), name = 'Log1p(smoothed gene expression)')) +
    geom_half_violin(data = pdata %>% filter(data_type == "GeneExpression"),
                     aes(x = cell_type, y = log1p(value), color = cell_type), fill = "white", adjust = 1.2,
                     side = "r", scale = "width", width = 0.9, alpha = 0.3, show.legend = FALSE) +
    scale_color_manual(values = mmjoint_colors) +
    geom_point(data = pdata, aes(x = cell_type, y = log1p(value), fill = data_type),
               stat = 'summary', fun = median, size = 0.5,
               position = position_dodge(width = -0.8),
               show.legend = FALSE) +
    stat_summary(data = pdata, aes(x = cell_type, y = log1p(value), fill = data_type),
                 fun.min = function(x){quantile(x)[2]},
                 fun.max = function(x){quantile(x)[4]},
                 geom = 'errorbar', color='black',
                 width = 0.01, linewidth = 0.3,
                 position = position_dodge(width = -0.8),
                 show.legend = FALSE) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
          axis.text.y = element_text(color = "black"),
          legend.box = "vertical",
          panel.grid = element_blank())
  
  
  
  
  
  gene <- "Klrg1"
  mat_gs <- assay(gs_mat[which(rowData(gs_mat)$name == gene),])
  mat_ge <- assay(ge_mat[which(rowData(ge_mat)$name == gene),])
  
  plot_data_Klrg1 <- data.frame(
    cell_id = proj$cellNames,
    cell_type = proj$cell_module,
    GeneScore = calculated_smoothed_ge(weightList, mat = mat_gs)[proj$cellNames],
    GeneExpression = calculated_smoothed_ge(weightList, mat = mat_ge)[proj$cellNames]
  ) 
  rt <- mean(plot_data_Klrg1$GeneScore)/mean(plot_data_Klrg1$GeneExpression)
  # rt <- 1
  pdata <- plot_data_Klrg1 %>% dplyr::mutate(GeneExpression = GeneExpression*rt)
  pdata <- pdata %>%
    pivot_longer(cols = c("GeneScore", "GeneExpression"),
                 names_to = "data_type",
                 values_to = "value")
  pdata$cell_type <- factor(pdata$cell_type, levels = c("Quiescent", "Priming_1", "Priming_2", "Activation_1", "Activation_2", 
                                                        "Activation_3", "Effector_1", "Effector_2", "Effector-memory_transition_1", 
                                                        "Effector-memory_transition_2", "Effector-memory_transition_3", 
                                                        "Memory_1", "Memory_2", "Memory_3"))
  
  
  ggplot()+
    geom_half_violin(data = pdata %>% filter(data_type == "GeneScore"),
                     aes(x = cell_type, y = log1p(value), fill = cell_type, color = cell_type), adjust = 1.2, 
                     side = "l", scale = "width", width = 0.9, alpha = 0.8, show.legend = FALSE) +
    scale_fill_manual(values = mmjoint_colors) +
    labs(title = gene, 
         x = "cell module", 
         y = "Log1p(smoothed genescore)") +
    scale_y_continuous(limits = c(0, 1*log1p(max(pdata$value))),
                       sec.axis = sec_axis(~./(log1p(mean(plot_data_Klrg1$GeneScore))/log1p(mean(plot_data_Klrg1$GeneExpression))), name = 'Log1p(smoothed gene expression)')) +
    geom_half_violin(data = pdata %>% filter(data_type == "GeneExpression"),
                     aes(x = cell_type, y = log1p(value), color = cell_type), fill = "white", adjust = 1.2,
                     side = "r", scale = "width", width = 0.9, alpha = 0.3, show.legend = FALSE) +
    scale_color_manual(values = mmjoint_colors) +
    geom_point(data = pdata, aes(x = cell_type, y = log1p(value), fill = data_type),
               stat = 'summary', fun = median, size = 0.5,
               position = position_dodge(width = -0.8),
               show.legend = FALSE) +
    stat_summary(data = pdata, aes(x = cell_type, y = log1p(value), fill = data_type),
                 fun.min = function(x){quantile(x)[2]},
                 fun.max = function(x){quantile(x)[4]},
                 geom = 'errorbar', color='black',
                 width = 0.01, linewidth = 0.3,
                 position = position_dodge(width = -0.8),
                 show.legend = FALSE) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
          axis.text.y = element_text(color = "black"),
          legend.box = "vertical",
          panel.grid = element_blank())
  
  
  
  
  
  gene <- "Bcl2"
  mat_gs <- assay(gs_mat[which(rowData(gs_mat)$name == gene),])
  mat_ge <- assay(ge_mat[which(rowData(ge_mat)$name == gene),])
  
  plot_data_Bcl2 <- data.frame(
    cell_id = proj$cellNames,
    cell_type = proj$cell_module,
    GeneScore = calculated_smoothed_ge(weightList, mat = mat_gs)[proj$cellNames],
    GeneExpression = calculated_smoothed_ge(weightList, mat = mat_ge)[proj$cellNames]
  ) 
  rt <- mean(plot_data_Bcl2$GeneScore)/mean(plot_data_Bcl2$GeneExpression)
  # rt <- 1
  pdata <- plot_data_Bcl2 %>% dplyr::mutate(GeneExpression = GeneExpression*rt)
  pdata <- pdata %>%
    pivot_longer(cols = c("GeneScore", "GeneExpression"),
                 names_to = "data_type",
                 values_to = "value")
  pdata$cell_type <- factor(pdata$cell_type, levels = c("Quiescent", "Priming_1", "Priming_2", "Activation_1", "Activation_2", 
                                                        "Activation_3", "Effector_1", "Effector_2", "Effector-memory_transition_1", 
                                                        "Effector-memory_transition_2", "Effector-memory_transition_3", 
                                                        "Memory_1", "Memory_2", "Memory_3"))
  
  ggplot()+
    geom_half_violin(data = pdata %>% filter(data_type == "GeneScore"),
                     aes(x = cell_type, y = log1p(value), fill = cell_type, color = cell_type), adjust = 1.2, 
                     side = "l", scale = "width", width = 0.9, alpha = 0.8, show.legend = FALSE) +
    scale_fill_manual(values = mmjoint_colors) +
    labs(title = gene, 
         x = "cell module", 
         y = "Log1p(smoothed genescore)") +
    scale_y_continuous(limits = c(0, 1*log1p(max(pdata$value))),
                       sec.axis = sec_axis(~./(log1p(mean(plot_data_Bcl2$GeneScore))/log1p(mean(plot_data_Bcl2$GeneExpression))), name = 'Log1p(smoothed gene expression)')) +
    geom_half_violin(data = pdata %>% filter(data_type == "GeneExpression"),
                     aes(x = cell_type, y = log1p(value), color = cell_type), fill = "white", adjust = 1.2,
                     side = "r", scale = "width", width = 0.9, alpha = 0.3, show.legend = FALSE) +
    scale_color_manual(values = mmjoint_colors) +
    geom_point(data = pdata, aes(x = cell_type, y = log1p(value), fill = data_type),
               stat = 'summary', fun = median, size = 0.5,
               position = position_dodge(width = -0.8),
               show.legend = FALSE) +
    stat_summary(data = pdata, aes(x = cell_type, y = log1p(value), fill = data_type),
                 fun.min = function(x){quantile(x)[2]},
                 fun.max = function(x){quantile(x)[4]},
                 geom = 'errorbar', color='black',
                 width = 0.01, linewidth = 0.3,
                 position = position_dodge(width = -0.8),
                 show.legend = FALSE) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
          axis.text.y = element_text(color = "black"),
          legend.box = "vertical",
          panel.grid = element_blank())
  
}


