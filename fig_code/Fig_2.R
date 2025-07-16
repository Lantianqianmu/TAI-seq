# figure 2
library(ggplot2)
library(ggforce)
library(ggalluvial)
library(patchwork)
library(Seurat)
library(ArchR)
library(scplotter)
library(pheatmap)
source("/home/zeemeeuw/YangLab/JoINT-seq/QC/code/utils.R")
setwd("/home/zeemeeuw/YangLab/JoINT-seq/Manuscript/fig_data/fig_2_data")

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
time_colors <-  c("WT" = "#ABD8E6", "d3" = "#35a153", "d5" = "#1f6933", "d8" = "#f26a11", "d14" = "#B22222" , "d30" = "#B17BA6", "R5" = "#85B0FF", "Y5" = "#3648A2")
tissue_colors <-  c("Spleen" = "#36A23D", "Ln" = "#2884E7", "Liver" = "#817ab9", "Blood" = "#EE3B00")

# fig 2C: cell module UMAP
{
  data <- as_tibble(seu@reductions$UMAP_LSI_Combined_batch3.3@cell.embeddings)
  data <- data %>% dplyr::bind_cols(seu@meta.data) %>% 
    dplyr::select(UMAPCombinedbatch33_1, UMAPCombinedbatch33_2, cell_module, tissue, timepoint)
  write_tsv(data, "source_fig_2C.tsv", na = "")
  CellDimPlot(seu, reduction = "UMAP_LSI_Combined_batch3.3",  theme = "theme_blank", legend.position = "right", group_by = "cell_module", raster = TRUE, pt_size = 0.1, palcolor = mmjoint_colors)
}

# fig 2D: Projecting TCR on UMAP
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
  xrange <- max(data$UMAPCombinedbatch33_1) - min(data$UMAPCombinedbatch33_1)
  yrange <- max(data$UMAPCombinedbatch33_2) - min(data$UMAPCombinedbatch33_2)
  ggplot(data) +
    theme_void() +
    coord_fixed(xrange/yrange) +
    geom_point(data = data_tcr, aes(x = UMAPCombinedbatch33_1, y = UMAPCombinedbatch33_2, size = log10(n_clonotype + 1), color = cell_module), alpha = 1) +
    scale_size_continuous(range = c(0.1, 4)) +
    scale_color_manual(values = mmjoint_colors)
  data <- data %>% 
    dplyr::select(UMAPCombinedbatch33_1, UMAPCombinedbatch33_2, cell_module, tissue, timepoint, n_clonotype) %>% 
    as_tibble()
  write_tsv(data, "source_fig_2D.tsv", na = "")
}

# fig 2E: gene module score
{
  load("calculated_seurat_gene_module.Rdata")
  seu$timepoint <- factor(seu$timepoint, levels = c("WT", "d3", "d5", "d8", "d14", "d30", "R5", "Y5"))
  pdata <- seu@meta.data %>% 
    dplyr::select(tissue, timepoint) %>% 
    bind_cols(module_score_cd8) %>% 
    dplyr::filter(tissue == "Ln") %>% 
    group_by(timepoint) %>% 
    dplyr::select(-tissue) %>% 
    dplyr::summarise_all(~mean(.x)) %>% 
    tibble::column_to_rownames("timepoint")
  pdata <- t(pdata)
  pdata <- pdata[!(rownames(pdata) %in% c("Stress response", "Anergy", "Exhaustion", "Senescence")),c("WT", "d3", "d5", "d8", "d14", "d30", "R5", "Y5")]
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
                     color = rev(choose_colorset("RdBu", 101))[5:95],
                     width = 5.5, height = 4.5)
  pdata <- data.frame(pdata) %>% 
    tibble::rownames_to_column(var = "module")
  write.table(pdata, "source_fig_2E.tsv", sep = "\t", quote = FALSE, na = "", row.names = FALSE, col.names = TRUE)
  
}



# fig 2H: sample distribution on umap
{
  data <- as_tibble(seu@reductions$UMAP_LSI_Combined_batch3.3@cell.embeddings)
  data <- data %>% dplyr::bind_cols(seu@meta.data) %>% 
    dplyr::select(UMAPCombinedbatch33_1, UMAPCombinedbatch33_2, cell_module, tissue, timepoint, Tissue_time)
  write_tsv(data, "source_fig_2H.tsv", na = "")
  CellDimPlot(seu, reduction = "UMAP_LSI_Combined_batch3.3",  theme = "theme_blank", legend.position = "none", group_by = "cell_module", facet_by = "Tissue_time",
              raster = TRUE, bg_color = "grey80", palcolor = mmjoint_colors, highlight = TRUE)
}

# fig 2I: alluvial plots of cluster proportion per tissue
{
  meta <- seu@meta.data %>% 
    group_by(timepoint, tissue, cell_module) %>% 
    dplyr::summarise(count = n()) %>% 
    group_by(timepoint, tissue) %>% 
    dplyr::mutate(freq = count/sum(count))
  meta$timepoint <- factor(meta$timepoint, levels = c("WT", "d3", "d5", "d8", "d14", "d30", "R5", "Y5"))
  meta$tissue <- factor(meta$tissue, levels = c("Ln", "Spleen", "Liver", "Blood"))
  write_tsv(meta, "source_fig_2I.tsv", na = "")

  ggplot(meta, aes(x = timepoint, y = freq, stratum = cell_module, alluvium = cell_module, fill = cell_module)) +
    geom_col(position = "stack", width = 0.5, alpha = 1) +
    geom_flow(width = 0.5, knot.pos = 0.3, alpha = 1) +
    theme_classic() +
    scale_fill_manual(values = mmjoint_colors) +
    xlab("Timepoint") +
    ylab("Proportion") +
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 10, color = "black"),
          axis.text.y = element_text(hjust = 0.5, vjust = 0.5, size = 10, color = "black"),
          strip.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 10, color = "black")) +
    facet_wrap(~tissue)
}

# fig 2J: pie chart memory
{
  meta <- seu@meta.data %>% 
    group_by(timepoint, tissue, cell_module) %>% 
    dplyr::summarise(count = n()) %>% 
    group_by(timepoint, tissue) %>% 
    dplyr::mutate(freq = count/sum(count))
  meta$timepoint <- factor(meta$timepoint, levels = c("WT", "d3", "d5", "d8", "d14", "d30", "R5", "Y5"))
  meta$tissue <- factor(meta$tissue, levels = c("Ln", "Spleen", "Liver", "Blood"))
  
  
  meta_memory <- meta %>% 
    dplyr::filter(timepoint != "WT") %>% 
    dplyr::filter(str_detect(cell_module, "^Memory")) %>% 
    group_by(timepoint, tissue) %>% 
    dplyr::summarise(prop = sum(freq)) %>% 
    dplyr::mutate(focus = if_else(timepoint %in% c("d3", "d5", "d8", "d14"), 0.1, 0)) %>% 
    group_by(tissue) %>% 
    dplyr::mutate(freq = prop/sum(prop)) %>% 
    dplyr::mutate(start = -(2 * pi - cumsum(lag(freq, default = 0)) * 2 * pi), 
                  end = -(2 * pi - cumsum(freq) * 2 * pi),
                  mid = (start + end) / 2,
                  label = sprintf("%.1f%%", freq * 100)) 
  
  
  p1 <- ggplot(meta_memory %>% dplyr::filter(tissue == "Ln")) +
    ggtitle("Ln") +
    geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 1, amount = prop, fill = timepoint, explode = focus), show.legend = FALSE, color = "white", linewidth = 0.3, stat = 'pie') +
    scale_fill_manual(values = time_colors) +
    theme_void() +
    coord_fixed() +
    geom_text(aes(x = 0.7 * sin(mid), y = 0.7 * cos(mid), label = label, angle = -mid * 180 / pi + 90), size = 4, color = "black") 
  p2 <- ggplot(meta_memory %>% dplyr::filter(tissue == "Spleen")) +
    ggtitle("Spleen") +
    geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 1, amount = prop, fill = timepoint, explode = focus), show.legend = FALSE, color = "white", linewidth = 0.3, stat = 'pie') +
    scale_fill_manual(values = time_colors) +
    theme_void() +
    coord_fixed() +
    geom_text(aes(x = 0.7 * sin(mid), y = 0.7 * cos(mid), label = label, angle = -mid * 180 / pi + 90), size = 4, color = "black") 
  p3 <- ggplot(meta_memory %>% dplyr::filter(tissue == "Liver")) +
    ggtitle("Liver") +
    geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 1, amount = prop, fill = timepoint, explode = focus), show.legend = FALSE, color = "white", linewidth = 0.3, stat = 'pie') +
    scale_fill_manual(values = time_colors) +
    theme_void() +
    coord_fixed() +
    geom_text(aes(x = 0.7 * sin(mid), y = 0.7 * cos(mid), label = label, angle = -mid * 180 / pi + 90), size = 4, color = "black") 
  p4 <- ggplot(meta_memory %>% dplyr::filter(tissue == "Blood")) +
    ggtitle("Blood") +
    geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 1, amount = prop, fill = timepoint, explode = focus), show.legend = FALSE, color = "white", linewidth = 0.3, stat = 'pie') +
    scale_fill_manual(values = time_colors) +
    theme_void() +
    coord_fixed() +
    geom_text(aes(x = 0.7 * sin(mid), y = 0.7 * cos(mid), label = label, angle = -mid * 180 / pi + 90), size = 4, color = "black") 
  wrap_plots(list(p1, p2, p3, p4), ncol = 2)
  
  write_tsv(meta_memory, "source_fig_2J.tsv", na = "")
}
