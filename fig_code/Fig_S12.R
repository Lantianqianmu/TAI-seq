library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(ggplot2)
library(Biostrings)
library(ggpubr)
library(ggforce)
library(patchwork)
library(ggnewscale)
library(readxl)
library(openxlsx)
library(scattermore)
options(repos="http://mirrors.tuna.tsinghua.edu.cn/CRAN")
options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor")
source("/home/zeemeeuw/YangLab/JoINT-seq/QC/code/utils.R")

setwd("/home/zeemeeuw/YangLab/JoINT-seq/Manuscript/fig_data/fig_S12_data")

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

load(paste0(RDSdir, "calculated_seurat_gene_module.Rdata"))

# A
{
  dft <- seu@meta.data %>% as_tibble() %>% 
    dplyr::select(timepoint, tissue, cell_module, chain1, chain2) 
  d1 <- dft %>% dplyr::select(timepoint, tissue, cell_module, chain1) %>% rename_with(~c("timepoint", "tissue", "cell_module", "chain"))
  d2 <- dft %>% dplyr::select(timepoint, tissue, cell_module, chain2) %>% rename_with(~c("timepoint", "tissue", "cell_module", "chain"))
  dft <- bind_rows(d1, d2) %>% 
    dplyr::filter(!is.na(chain)) %>% 
    group_by(timepoint, chain) %>% 
    dplyr::mutate(count_timepoint = n())
  
  # d8 
  dft_d8 <- dft %>% 
    dplyr::filter(timepoint == "d8") %>% 
    dplyr::filter(count_timepoint > 30) %>% 
    group_by(chain) %>% 
    dplyr::mutate(MEratio = sum(cell_module %in% c("Effector-memory_transition_1", "Effector-memory_transition_2", "Effector-memory_transition_3")) / sum(cell_module %in% c("Effector_1", "Effector_2"))) %>% 
    dplyr::mutate(group = case_when(sum(cell_module %in% c("Effector_1", "Effector_2"))/n() > 0.7 ~ "Eff biased",
                                    sum(cell_module %in% c("Effector-memory_transition_1", "Effector-memory_transition_2", "Effector-memory_transition_3"))/n() > 0.7 ~ "EMT biased",
                                    .default = "no bias")) %>% 
    dplyr::mutate(has_memory = case_when(sum(cell_module %in% c("Memory_1", "Memory_2", "Memory_3")) > 0 ~ "Yes",
                                         .default = "No"))
  dft_d8$group <- factor(dft_d8$group, levels = c("Eff biased", "no bias", "EMT biased")) 
  dft_d8 <- dft_d8 %>% arrange(MEratio)
  dft_d8$chain <- factor(dft_d8$chain, levels = unique(dft_d8$chain))
  
  p1 <- ggplot(dft_d8) +
    geom_bar(aes(x = chain, fill = tissue), position = "fill") +
    scale_fill_manual(values = tissue_colors) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"),
          plot.margin = margin(1, 1, 1, 30)) +
    ggtitle("d8 clonotypes (>30)") +
    coord_fixed(ratio = 0.3*n_distinct(dft_d8$chain))
  
  p2 <- ggplot(dft_d8) +
    geom_bar(aes(x = chain, fill = cell_module), position = "fill") +
    scale_fill_manual(values = mmjoint_colors) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"),
          plot.margin = margin(1, 1, 1, 30)) +
    ggtitle("d8 clonotypes (>30)") +
    coord_fixed(ratio = 0.3*n_distinct(dft_d8$chain))
  
  wrap_plots(list(p1, p2), ncol = 2)
}


# B
{
  dft <- seu@meta.data %>% as_tibble() %>% 
    dplyr::select(timepoint, tissue, cell_module, chain1, chain2) 
  d1 <- dft %>% dplyr::select(timepoint, tissue, cell_module, chain1) %>% rename_with(~c("timepoint", "tissue", "cell_module", "chain"))
  d2 <- dft %>% dplyr::select(timepoint, tissue, cell_module, chain2) %>% rename_with(~c("timepoint", "tissue", "cell_module", "chain"))
  dft <- bind_rows(d1, d2) %>% 
    dplyr::filter(!is.na(chain)) %>% 
    group_by(timepoint, chain) %>% 
    dplyr::mutate(count_timepoint = n())
  
  dft_d14 <- dft %>% 
    dplyr::filter(timepoint == "d14") %>% 
    dplyr::filter(count_timepoint > 10) %>% 
    group_by(chain) %>% 
    dplyr::mutate(MEratio = sum(cell_module %in% c("Effector-memory_transition_1", "Effector-memory_transition_2", "Effector-memory_transition_3")) / sum(cell_module %in% c("Effector_1", "Effector_2"))) %>% 
    dplyr::mutate(group = case_when(sum(cell_module %in% c("Effector_1", "Effector_2"))/n() > 0.7 ~ "Eff biased",
                                    sum(cell_module %in% c("Effector-memory_transition_1", "Effector-memory_transition_2", "Effector-memory_transition_3"))/n() > 0.7 ~ "EMT biased",
                                    .default = "no bias"))  %>% 
    dplyr::mutate(has_memory = case_when(sum(cell_module %in% c("Memory_1", "Memory_2", "Memory_3")) > 0 ~ "Yes", .default = "No"))
  dft_d14$group <- factor(dft_d14$group, levels = c("Eff biased", "no bias", "EMT biased")) 
  dft_d14 <- dft_d14 %>% arrange(MEratio)
  dft_d14$chain <- factor(dft_d14$chain, levels = unique(dft_d14$chain))
  dft_d14$chain <- factor(dft_d14$chain, levels = c("CALSDPNTGKLTF++CASRRDRNQDTQYF", "CATVASSGSWQLIF++CASSMWGYEQYF", 
                                                    "CAMREGQGGSAKLIF++CASSPDREGNLNTEVFF", "CAASRPSGSWQLIF++CAWSLTLANTGQLYF", 
                                                    "CAADSSSGSWQLIF++CASSQDPTGGKSQNTLYF", "CAASRTSSGQKLVF++CASSRTGGASF", 
                                                    "CALSDQGVTGNTGKLIF++CASSRTGEQYF", 
                                                    "CAVNGYSNNRLTL++CASSRPRTGGSYEQYF", 
                                                    "CAAPPGTGGYKVVF++CASSRSSNSDYTF", "CALWELASSGQKLVF++CASSRPRTGGSYEQYF", 
                                                    "CASRNNYAQGLTF++CASGDAGTVNF", 
                                                    "CALSAGSWQLIF++CASSGDTYAEQFF", "CAASRTSSGQKLVF++CASSRTGEQYF", 
                                                    "CAVSRTGGYKVVF++CASSPRQTDYTF", "CAANSGTYQRF++CASGDAQ_KDTQYF", 
                                                    "CALGGTGGYKVVF++CASSDSSNERLFF", "CAASESSGTYQRF++CASSPGAATGQLYF", 
                                                    "CAIGAGNTGKLIF++CASGDAGAPLF", 
                                                    "CAVSDPSSGQKLVF++CASRDTGQLYF", 
                                                    "CALWELESNNRIFF++CASSDSSNERLFF", # yes
                                                    "CAVRYPRNNNNAPRF++CASSDSSNERLFF", # yes
                                                    "CALGVWGQLIF++CASSLSGQSQNTLYF", # yes
                                                    "CAANSGTYQRF++CASSRDWGLLSQNTLYF", # yes
                                                    "CAVNTGGLSGKLTF++CASSRQGDTGQLYF", # yes
                                                    "CAVSMNSNNRIFF++CASGDAGAPLF", # yes
                                                    "CAVSMTGGYKVVF++CASSPRQSDYTF", 
                                                    "CAMRDYGSSGNKLIF++CASSLRDNNSGNTLYF", 
                                                    "CALGGTGGYKVVF++CASSPRGDEQYF", "CALRNSGTYQRF++CASGDRGNTEVFF", 
                                                    "CAMREGGGGSNYKLTF++CASSLDRGGNLNERLFF", "CVLGGYAQGLTF++CASGRAQDTQYF", 
                                                    "CALGGTGGYKVVF++CASSPRDSDYTF", 
                                                    "CAASENAYKVIF++CASSDRSSAETLYF"))
  
  
  p1 <- ggplot(dft_d14) +
    geom_bar(aes(x = chain, fill = tissue), position = "fill") +
    scale_fill_manual(values = tissue_colors) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"),
          plot.margin = margin(1, 1, 1, 30)) +
    ggtitle("d14 clonotypes (>10)") +
    coord_fixed(ratio = 0.3*n_distinct(dft_d14$chain))
  
  p2 <- ggplot(dft_d14) +
    geom_bar(aes(x = chain, fill = cell_module), position = "fill") +
    scale_fill_manual(values = mmjoint_colors) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"),
          plot.margin = margin(1, 1, 1, 30)) +
    ggtitle("d14 clonotypes (>10)") +
    coord_fixed(ratio = 0.3*n_distinct(dft_d14$chain))
  
  wrap_plots(list(p1, p2), ncol = 2)
}




# Fig C: gene module score of d8 clones
{
  tcr_paired <- readxl::read_excel("Lm-TCR4.57-clonotypes-paired in rna.xlsx") %>%
    dplyr::mutate(tissue = str_extract(Tissue_time, "(.*)-", group = 1)) %>%
    group_by(CDRa_CDRb, timepoint, tissue, barcode) %>%
    distinct(CDRa_CDRb, timepoint, tissue, barcode, .keep_all = TRUE) %>%
    group_by(CDRa_CDRb, timepoint, tissue) %>%
    dplyr::mutate(counts = n()) %>%
    arrange(dplyr::desc(counts))
  
  dft <- seu@meta.data %>% as_tibble() %>% 
    dplyr::select(timepoint, tissue, cell_module, chain1, chain2) 
  d1 <- dft %>% dplyr::select(timepoint, tissue, cell_module, chain1) %>% rename_with(~c("timepoint", "tissue", "cell_module", "chain"))
  d2 <- dft %>% dplyr::select(timepoint, tissue, cell_module, chain2) %>% rename_with(~c("timepoint", "tissue", "cell_module", "chain"))
  dft <- bind_rows(d1, d2) %>% 
    dplyr::filter(!is.na(chain)) %>% 
    group_by(timepoint, chain) %>% 
    dplyr::mutate(count_timepoint = n())
  
  dft_d8 <- dft %>% 
    dplyr::filter(timepoint == "d8") %>% 
    dplyr::filter(count_timepoint > 30) %>% 
    group_by(chain) %>% 
    dplyr::mutate(MEratio = sum(cell_module %in% c("Effector-memory_transition_1", "Effector-memory_transition_2", "Effector-memory_transition_3")) / sum(cell_module %in% c("Effector_1", "Effector_2"))) %>% 
    dplyr::mutate(group = case_when(sum(cell_module %in% c("Effector_1", "Effector_2"))/n() > 0.7 ~ "Eff biased",
                                    sum(cell_module %in% c("Effector-memory_transition_1", "Effector-memory_transition_2", "Effector-memory_transition_3"))/n() > 0.7 ~ "EMT biased",
                                    .default = "no bias"))
  dft_d8$group <- factor(dft_d8$group, levels = c("Eff biased", "no bias", "EMT biased")) 
  dft_d8 <- dft_d8 %>% arrange(MEratio)
  dft_d8$chain <- factor(dft_d8$chain, levels = unique(dft_d8$chain))
  d8_tcr <- levels(dft_d8$chain)
  
  pdata <- lapply(seq_along(d8_tcr), function(i){
    tcr_i <- d8_tcr[i]
    bc_i <- tcr_paired %>% 
      dplyr::filter(timepoint == "d8") %>% 
      dplyr::filter(CDRa_CDRb == tcr_i) %>% 
      pull(barcode)
    colMeans(module_score_cd8[bc_i,])
  }) %>% Reduce("cbind", .)
  colnames(pdata) <- d8_tcr
  pheatmap::pheatmap(pdata[c("Lipid metabolism", "OXPHOS", "Glycolysis",  "Cytotoxicity", "TCR Signaling", "Chemokine/Chemokine receptor", "Cytokine/Cytokine receptor", 
                             "IFN Response",   "Activation/Effector function", "Adhesion", "Naive", "NFKB Signaling", "MAPK Signaling", 
                             "Pro-apoptosis", "Anti-apoptosis"),], 
                     border_color = NA,
                     show_colnames = TRUE,
                     show_rownames = TRUE,
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     legend_breaks = c(-1, 0, 1),
                     clustering_method = "ward.D",
                     scale = "row",
                     fontsize_row = 10,
                     fontsize_col = 10,
                     angle_col = 45,
                     # color = choose_colorset("wolfgang_extra", 101),
                     color = rev(choose_colorset("RdBu", 101))[5:95],
                     # filename = "d8_tcr_desc_module_score_cd8_heatmap.pdf",
                     width = 13, height = 6)
  
  pdata <- data.frame(pdata) %>% 
    tibble::rownames_to_column(var = "module")
  write.table(pdata, "source_fig_S12C.tsv", sep = "\t", quote = FALSE, na = "", row.names = FALSE, col.names = TRUE)
}

# Fig D: gene module score of d14 clones
{
  tcr_paired <- readxl::read_excel("Lm-TCR4.57-clonotypes-paired in rna.xlsx") %>%
    dplyr::mutate(tissue = str_extract(Tissue_time, "(.*)-", group = 1)) %>%
    group_by(CDRa_CDRb, timepoint, tissue, barcode) %>%
    distinct(CDRa_CDRb, timepoint, tissue, barcode, .keep_all = TRUE) %>%
    group_by(CDRa_CDRb, timepoint, tissue) %>%
    dplyr::mutate(counts = n()) %>%
    arrange(dplyr::desc(counts))
  
  dft <- seu@meta.data %>% as_tibble() %>% 
    dplyr::select(timepoint, tissue, cell_module, chain1, chain2) 
  d1 <- dft %>% dplyr::select(timepoint, tissue, cell_module, chain1) %>% rename_with(~c("timepoint", "tissue", "cell_module", "chain"))
  d2 <- dft %>% dplyr::select(timepoint, tissue, cell_module, chain2) %>% rename_with(~c("timepoint", "tissue", "cell_module", "chain"))
  dft <- bind_rows(d1, d2) %>% 
    dplyr::filter(!is.na(chain)) %>% 
    group_by(timepoint, chain) %>% 
    dplyr::mutate(count_timepoint = n())
  
  dft_d14 <- dft %>% 
    dplyr::filter(timepoint == "d14") %>% 
    dplyr::filter(count_timepoint > 10) %>% 
    group_by(chain) %>% 
    dplyr::mutate(MEratio = sum(cell_module %in% c("Effector-memory_transition_1", "Effector-memory_transition_2", "Effector-memory_transition_3")) / sum(cell_module %in% c("Effector_1", "Effector_2"))) %>% 
    dplyr::mutate(group = case_when(sum(cell_module %in% c("Effector_1", "Effector_2"))/n() > 0.7 ~ "Eff biased",
                                    sum(cell_module %in% c("Effector-memory_transition_1", "Effector-memory_transition_2", "Effector-memory_transition_3"))/n() > 0.7 ~ "EMT biased",
                                    .default = "no bias"))
  dft_d14$group <- factor(dft_d14$group, levels = c("Eff biased", "no bias", "EMT biased")) 
  dft_d14 <- dft_d14 %>% arrange(MEratio)
  dft_d14$chain <- factor(dft_d14$chain, levels = unique(dft_d14$chain))
  d14_tcr <- levels(dft_d14$chain)
  d14_tcr <- c("CALSDPNTGKLTF++CASRRDRNQDTQYF", "CATVASSGSWQLIF++CASSMWGYEQYF", 
               "CAMREGQGGSAKLIF++CASSPDREGNLNTEVFF", "CAASRPSGSWQLIF++CAWSLTLANTGQLYF", 
               "CAADSSSGSWQLIF++CASSQDPTGGKSQNTLYF", "CAASRTSSGQKLVF++CASSRTGGASF", 
               "CALSDQGVTGNTGKLIF++CASSRTGEQYF", 
               "CAVNGYSNNRLTL++CASSRPRTGGSYEQYF", 
               "CAAPPGTGGYKVVF++CASSRSSNSDYTF", "CALWELASSGQKLVF++CASSRPRTGGSYEQYF", 
               "CASRNNYAQGLTF++CASGDAGTVNF", 
               "CALSAGSWQLIF++CASSGDTYAEQFF", "CAASRTSSGQKLVF++CASSRTGEQYF", 
               "CAVSRTGGYKVVF++CASSPRQTDYTF", "CAANSGTYQRF++CASGDAQ_KDTQYF", 
               "CALGGTGGYKVVF++CASSDSSNERLFF", "CAASESSGTYQRF++CASSPGAATGQLYF", 
               "CAIGAGNTGKLIF++CASGDAGAPLF", 
               "CAVSDPSSGQKLVF++CASRDTGQLYF", 
               "CALWELESNNRIFF++CASSDSSNERLFF", # yes
               "CAVRYPRNNNNAPRF++CASSDSSNERLFF", # yes
               "CALGVWGQLIF++CASSLSGQSQNTLYF", # yes
               "CAANSGTYQRF++CASSRDWGLLSQNTLYF", # yes
               "CAVNTGGLSGKLTF++CASSRQGDTGQLYF", # yes
               "CAVSMNSNNRIFF++CASGDAGAPLF", # yes
               "CAVSMTGGYKVVF++CASSPRQSDYTF", 
               "CAMRDYGSSGNKLIF++CASSLRDNNSGNTLYF", 
               "CALGGTGGYKVVF++CASSPRGDEQYF", "CALRNSGTYQRF++CASGDRGNTEVFF", 
               "CAMREGGGGSNYKLTF++CASSLDRGGNLNERLFF", "CVLGGYAQGLTF++CASGRAQDTQYF", 
               "CALGGTGGYKVVF++CASSPRDSDYTF", 
               "CAASENAYKVIF++CASSDRSSAETLYF")
  
  # gene module score
  pdata <- lapply(seq_along(d14_tcr), function(i){
    tcr_i <- d14_tcr[i]
    bc_i <- tcr_paired %>% 
      dplyr::filter(timepoint == "d14") %>% 
      dplyr::filter(CDRa_CDRb == tcr_i) %>% 
      pull(barcode)
    colMeans(module_score_cd8[bc_i,])
  }) %>% Reduce("cbind", .)
  colnames(pdata) <- d14_tcr
  pheatmap::pheatmap(pdata[c("Lipid metabolism", "OXPHOS", "Glycolysis",  "Cytotoxicity", "TCR Signaling", "Chemokine/Chemokine receptor", "Cytokine/Cytokine receptor", 
                             "IFN Response",   "Activation/Effector function", "Adhesion", "Naive", "NFKB Signaling", "MAPK Signaling", 
                             "Pro-apoptosis", "Anti-apoptosis"),], 
                     border_color = NA,
                     show_colnames = TRUE,
                     show_rownames = TRUE,
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     legend_breaks = c(-1, 0, 1),
                     clustering_method = "complete",
                     scale = "row",
                     fontsize_row = 10,
                     fontsize_col = 10,
                     angle_col = 45,
                     # color = choose_colorset("wolfgang_extra", 101),
                     color = rev(choose_colorset("RdBu", 101))[5:95],
                     # filename = "d14_tcr_desc_module_score_cd8_heatmap.pdf",
                     width = 11, height = 6)
  
  pdata <- data.frame(pdata) %>% 
    tibble::rownames_to_column(var = "module")
  write.table(pdata, "source_fig_S12D.tsv", sep = "\t", quote = FALSE, na = "", row.names = FALSE, col.names = TRUE)
}


# E
{
  data <- read_tsv("source_fig_S12.tsv")
  ggplot(data, aes(x = checkpoint_3_predisposition, y = checkpoint_3_maintenance, fill = timepoint)) +
    geom_point(size = 3, shape = 21) +
    scale_fill_manual(values = time_colors) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(size = 10, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.title.x = element_text(size = 13),
          axis.title.y = element_text(size = 13),
          strip.background = element_rect(fill = "grey90"),
          strip.text = element_text(size = 10, color = "black")) +
    coord_fixed(ratio = (max(data$checkpoint_3_predisposition) - min(data$checkpoint_3_predisposition))/
                  (max(data$checkpoint_3_maintenance) - min(data$checkpoint_3_maintenance)))
}


