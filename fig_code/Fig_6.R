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

setwd("/home/zeemeeuw/YangLab/JoINT-seq/Manuscript/fig_data/fig_6_data")

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

color_tcr <- c("jurkat"= "#F2F2F2",
               "213"= "#F70F33",
               "214"= "#F2CCDD",
               "215"= "#99D25C",
               "216"= "#8CC23E",
               "217"= "#F96566",
               "218"= "#F01485",
               "219"= "#F264B0",
               "220"= "#D0ECEE",
               "221"= "#6AC7D9",
               "222"= "#F01788",
               "223"= "#0D3333",
               "224"= "#0F6837",
               "225"= "#34A344",
               "226"= "#168C43",
               "227"= "#F5CCE4",
               "228"= "#FA0F0C",
               "229"= "#F40F67",
               "230"= "#EF1385")
highlight <- c("213" = "CAVRYPRNNNNAPRF++CASSDSSNERLFF", # d14 high aff
               "214" = "CALWELESNNRIFF++CASSDSSNERLFF",
               "215" = "CAASEMAGGSNAKLTF++CASSGASAETLYF",
               "216" = "CVVGAGSNYQLIW++CASSGASAETLYF",
               "217" = "CALGEASSGQKLVF++CAWSLGLGNYAEQFF",
               "218" = "CALGGTGGYKVVF++CASSPRDSDYTF",
               "219" = "CAASGYNQGKLIF++CASRQGARDTQYF",
               "220" = "CAASLFYQGGRALIF++CASSQDWGSSAETLYF",
               "221" = "CALRTGGYKVVF++CASSQRGSQNTLYF",
               "222" = "CAASGTGNTGKLIF++CAWIPTVANTGQLYF", # d8 catch
               "223" = "CAIPNNYAQGLTF++CASRDTGQLYF",
               "224" = "CAIPNNYAQGLTF++CAWQTFNNQAPLF",
               "225" = "CVVGAGSNYQLIW++CGAKNRDRWEVFF",
               "226" = "CAASEMAGGSNAKLTF++CGAKNRDRWEVFF",
               "227" = "CAANDGSSGNKLIF++CASSRANSDYTF",
               "228" = "CAVRDSNYQLIW++CASSFSLNYAEQFF", # d8 high aff
               "229" = "CIVTDIRTGGYKVVF++CGARDQGNNNQAPLF",
               "230" = "CAASFSGNEKITF++CTCSAGGAGNTLYF")

load(paste0(RDSdir, "calculated_seurat_gene_module.Rdata"))
load(paste0(RDSdir, "gene_module.Rdata"))
load(paste0(RDSdir, "lm_merged_overlapped_harmony_mat_smoothed.RData"))
# load(paste0(RDSdir, "chromatin_potential_gene_module.Rdata"))
load(paste0(RDSdir, "checkpoint_mean_sd_manual_2.RData"))


tcr_paired <- readxl::read_excel("Lm-TCR4.57-clonotypes-paired in rna.xlsx") %>% 
  dplyr::mutate(tissue = str_extract(Tissue_time, "(.*)-", group = 1)) %>% 
  group_by(CDRa_CDRb, timepoint, tissue, barcode) %>% 
  distinct(CDRa_CDRb, timepoint, tissue, barcode, .keep_all = TRUE) %>% 
  ungroup() 

tcr_d8 <- tcr_paired %>% dplyr::filter(timepoint == "d8" & tissue %in% c("Spleen"))
bc_217 <- tcr_d8[tcr_d8$CDRa_CDRb == highlight["217"],]$barcode # slip
bc_218 <- tcr_d8[tcr_d8$CDRa_CDRb == highlight["218"],]$barcode # slip
bc_219 <- tcr_d8[tcr_d8$CDRa_CDRb == highlight["219"],]$barcode # catch
bc_222 <- tcr_d8[tcr_d8$CDRa_CDRb == highlight["222"],]$barcode # catch
tcr_d14 <- tcr_paired %>% dplyr::filter(timepoint == "d14" & tissue %in% c("Spleen", "Ln"))
bc_214 <- tcr_d14[tcr_d14$CDRa_CDRb == highlight["214"],]$barcode # catch


# 6A
{

  module_list <- gm_cd8
  gene_list_names <- names(gm_cd8)
  
  bc_slip <- c(bc_217, bc_218)
  bc_catch <- c(bc_219, bc_222)
  df <- lapply(seq_along(module_list), function(i) {
    gene_list_name <- gene_list_names[i]
    residual_slip <- module_score_cd8[bc_slip,][[gene_list_name]]
    data_slip <- tibble(potential = residual_slip, TCR = "slip")
    residual_catch <- module_score_cd8[bc_catch,][[gene_list_name]]
    data_catch <- tibble(potential = residual_catch, TCR = "catch")
    residual_214 <- module_score_cd8[bc_214,][[gene_list_name]]
    data_214 <- tibble(potential = residual_214, TCR = "214")
    data <- bind_rows(data_slip, data_catch, data_214)
    data$TCR <- factor(data$TCR, levels = c('slip', 'catch', '214'))
    t1 <- wilcox.test(data$potential[data$TCR == 'slip'], data$potential[data$TCR == 'catch'])
    t2 <- wilcox.test(data$potential[data$TCR == 'slip'], data$potential[data$TCR == '214'])
    t3 <- wilcox.test(data$potential[data$TCR == 'catch'], data$potential[data$TCR == '214'])
    tibble(group = gene_list_name, "slip_catch" = t1$p.value, "slip_214" = t2$p.value, "catch_214" = t3$p.value)
  }) %>% Reduce("rbind", .)


  for (i in seq_along(module_list)) {
    gene_list_name <- gene_list_names[i]
    residual_217 <- module_score_cd8[bc_217,][[gene_list_name]]
    data_217 <- tibble(potential = residual_217, TCR = "217")
    residual_218 <- module_score_cd8[bc_218,][[gene_list_name]]
    data_218 <- tibble(potential = residual_218, TCR = "218")
    residual_219 <- module_score_cd8[bc_219,][[gene_list_name]]
    data_219 <- tibble(potential = residual_219, TCR = "219")
    residual_222 <- module_score_cd8[bc_222,][[gene_list_name]]
    data_222 <- tibble(potential = residual_222, TCR = "222")
    residual_214 <- module_score_cd8[bc_214,][[gene_list_name]]
    data_214 <- tibble(potential = residual_214, TCR = "214")
    data <- bind_rows(data_217, data_218, data_219, data_222, data_214)
    data$TCR <- factor(data$TCR, levels = c('217', '218', '222', '219', '214'))
    
    p <- ggplot(data, aes(x = TCR, y = potential, color = TCR)) +
      theme_classic() +
      geom_jitter(size = 2, alpha = 0.5, width = 0.3) +
      geom_boxplot(width = 0.5, fill = NA, color = "black", outlier.shape = NA) +
      geom_segment(x = 1.5, xend = 3.5, y = 0.7, yend = 0.7, color = "#000000", linewidth = 0.3) +
      geom_segment(x = 3.5, xend = 5, y = 0.8, yend = 0.8, color = "#000000", linewidth = 0.3) +
      geom_segment(x = 1.5, xend = 5, y = 0.9, yend = 0.9, color = "#000000", linewidth = 0.3) +
      annotate("text", size = 2, color = "#000000",  x = 2.5, y = 0.73, label = df$slip_catch[i]) +
      annotate("text", size = 2, color = "#000000",  x = 4.25, y = 0.83, label = df$catch_214[i]) +
      annotate("text", size = 2, color = "#000000",  x = 3.25, y = 0.93, label = df$slip_214[i]) +
      ylim(-0.5, 1) +
      ylab(gene_list_name) +
      theme(panel.grid = element_blank(),
            axis.text.x = element_text(size = 10, color = "black"),
            axis.text.y = element_text(size = 10, color = "black"),
            axis.title.x = element_text(size = 13),
            axis.title.y = element_text(size = 13)) +
      scale_color_manual(values = color_tcr)
    print(p)
    
    
  }
}

# 6B
{
  

  fname <- c("checkpoint_1_predisposition", "checkpoint_1_maintenance",
             "checkpoint_2_predisposition", "checkpoint_2_maintenance", 
             "checkpoint_3_predisposition", "checkpoint_3_maintenance", 
             "checkpoint_4_predisposition", "checkpoint_4_maintenance")
  genes_predi_score <- list("checkpoint_1" = c("Dock9", "Ly6e", "Ptpn13", "Zfp281"),
                            # "checkpoint_2" = c("Ccl5", "Ccr2", "Ccr5", "Cx3cr1", "Fasl", "Gsap", "Gzma", "Il18r1", "Il18rap", "Itga4", "Itgb1", "Klrc1", "Klrk1", "Lamc1", "Sntb2", "Tbx21", "Tcf3", "Ttc39c", "Zeb2"),
                            "checkpoint_2" = c("Cyba", "Cyth4", "Dock5", "Fkbp5", "Itga4", "Itgb1", "Lamc1"),
                            "checkpoint_3" = c("Bcl6", "Chd3", "Dennd4a", "Etf1", "Etv3", "Ifrd1", "Junb", "Vezf1"),
                            "checkpoint_4" = c("Iigp1", "Ly6a", "Nap1l1", "Ncapd2", "Nedd4l", "Pcgf5"))
  
  genes_maint_score <- list("checkpoint_1" = c("Ccm2", "Ccr7", "Ccr9", "Dcaf17", 
                                               "Hspbap1", "Hvcn1", "Il4ra", "Il6ra", "Itga6", "Mettl8", "Skp1a", "Smc4",  "St6gal1", 
                                               "Tet1", "Trib2", "Trio", "Ttc28", "Vamp1", "Zdhhc14"),
                            "checkpoint_2"  = c("Arhgap25", "Batf", "Cers6", "Fchsd2", "Ggact", "Klf3", "Naa16", "Spag9"),
                            "checkpoint_3" = c("Anxa2", "Gzmb", "Gzmk", "Mki67", "Nkg7", "Ssrp1", "Ubb"),
                            "checkpoint_4"  = c("Gnptab", "Ikzf3", "N4bp1", "Osbpl3"))
  g1 <- lapply(seq_along(genes_predi_score), function(i){
    out <- tibble(checkpoint = names(genes_predi_score)[i], type = "predisposition", genes = genes_predi_score[[i]])
    out$group <- paste0(out$checkpoint, "_", out$type)
    return(out)
  }) %>% Reduce("rbind", .)
  g2 <- lapply(seq_along(genes_maint_score), function(i){
    out <- tibble(checkpoint = names(genes_maint_score)[i], type = "maintenance", genes = genes_maint_score[[i]])
    out$group <- paste0(out$checkpoint, "_", out$type)
    return(out)
  }) %>% Reduce("rbind", .)
  genes_chrom_potential <- rbind(g1, g2)
  genes_chrom_potential_2 <- genes_chrom_potential %>% dplyr::mutate(group = case_match(checkpoint, "checkpoint_1" ~ "checkpoint_1_all",
                                                                                        "checkpoint_2" ~ "checkpoint_2_all",
                                                                                        "checkpoint_3" ~ "checkpoint_3_all",
                                                                                        "checkpoint_4" ~ "checkpoint_4_all"))
  genes_chrom_potential <- genes_chrom_potential %>% bind_rows(genes_chrom_potential_2)
  fname <- c("checkpoint_1_predisposition", "checkpoint_1_maintenance", "checkpoint_1_all",
             "checkpoint_2_predisposition", "checkpoint_2_maintenance", "checkpoint_2_all", 
             "checkpoint_3_predisposition", "checkpoint_3_maintenance", "checkpoint_3_all",
             "checkpoint_4_predisposition", "checkpoint_4_maintenance", "checkpoint_4_all")
  

  bc_slip <- c(bc_217, bc_218)
  bc_catch <- c(bc_219, bc_222)
  df <- lapply(seq_along(fname), function(i) {
    gp <- fname[i]
    mean_gp <- mean_checkpoint[[gp]]
    sd_gp <- sd_checkpoint[[gp]]
    genes <- dplyr::filter(genes_chrom_potential, group == gp) %>% pull(genes)
    
    residual_slip <- colMeans(dorcMat_sm_norm[genes, bc_slip] - rnaMat_sm_norm[genes, bc_slip])
    residual_slip <- (residual_slip - mean_gp)/(sd_gp)
    data_slip <- tibble(potential = residual_slip, TCR = "slip")
    residual_catch <- colMeans(dorcMat_sm_norm[genes, bc_catch] - rnaMat_sm_norm[genes, bc_catch])
    residual_catch <- (residual_catch - mean_gp)/(sd_gp)
    data_catch<- tibble(potential = residual_catch, TCR = "catch")
    residual_214 <- colMeans(dorcMat_sm_norm[genes, bc_214] - rnaMat_sm_norm[genes, bc_214])
    residual_214 <- (residual_214 - mean_gp)/(sd_gp)
    data_214 <- tibble(potential = residual_214, TCR = "214")
    data <- bind_rows(data_slip, data_catch, data_214)
    data$TCR <- factor(data$TCR, levels = c('slip', 'catch', '214'))
    t1 <- wilcox.test(data$potential[data$TCR == 'slip'], data$potential[data$TCR == 'catch'])
    t2 <- wilcox.test(data$potential[data$TCR == 'slip'], data$potential[data$TCR == '214'])
    t3 <- wilcox.test(data$potential[data$TCR == 'catch'], data$potential[data$TCR == '214'])
    tibble(group = gp, "slip_catch" = t1$p.value, "slip_214" = t2$p.value, "catch_214" = t3$p.value)
  }) %>% Reduce("rbind", .)
  
  
  
  for (i in seq_along(fname)) {
    gp <- fname[i]
    mean_gp <- mean_checkpoint[[gp]]
    sd_gp <- sd_checkpoint[[gp]]
    genes <- dplyr::filter(genes_chrom_potential, group == gp) %>% pull(genes)
    
    residual_217 <- colMeans(dorcMat_sm_norm[genes, bc_217] - rnaMat_sm_norm[genes, bc_217])
    residual_217 <- (residual_217 - mean_gp)/(sd_gp)
    data_217 <- tibble(potential = residual_217, TCR = "217")
    residual_218 <- colMeans(dorcMat_sm_norm[genes, bc_218] - rnaMat_sm_norm[genes, bc_218])
    residual_218 <- (residual_218 - mean_gp)/(sd_gp)
    data_218 <- tibble(potential = residual_218, TCR = "218")
    residual_219 <- colMeans(dorcMat_sm_norm[genes, bc_219] - rnaMat_sm_norm[genes, bc_219])
    residual_219 <- (residual_219 - mean_gp)/(sd_gp)
    data_219 <- tibble(potential = residual_219, TCR = "219")
    residual_222 <- colMeans(dorcMat_sm_norm[genes, bc_222] - rnaMat_sm_norm[genes, bc_222])
    residual_222 <- (residual_222 - mean_gp)/(sd_gp)
    data_222 <- tibble(potential = residual_222, TCR = "222")
    residual_214 <- colMeans(dorcMat_sm_norm[genes, bc_214] - rnaMat_sm_norm[genes, bc_214])
    residual_214 <- (residual_214 - mean_gp)/(sd_gp)
    data_214 <- tibble(potential = residual_214, TCR = "214")
    
    
    data <- bind_rows(data_217, data_218, data_219, data_222, data_214)
    data$TCR <- factor(data$TCR, levels = c('217', '218', '222', '219', '214'))
    
    p <- ggplot(data, aes(x = TCR, y = potential, color = TCR)) +
      theme_classic() +
      geom_jitter(size = 2, alpha = 0.5, width = 0.3) +
      geom_boxplot(width = 0.5, fill = NA, color = "black", outlier.shape = NA) +
      geom_segment(x = 1.5, xend = 3.5, y = 3, yend = 3, color = "#000000", linewidth = 0.3) +
      geom_segment(x = 3.5, xend = 5, y = 3.5, yend = 3.5, color = "#000000", linewidth = 0.3) +
      geom_segment(x = 1.5, xend = 5, y = 4, yend = 4, color = "#000000", linewidth = 0.3) +
      annotate("text", size = 2, color = "#000000",  x = 2.5, y = 3.2, label = df$slip_catch[i]) +
      annotate("text", size = 2, color = "#000000",  x = 4.25, y = 3.7, label = df$catch_214[i]) +
      annotate("text", size = 2, color = "#000000",  x = 3.25, y = 4.2, label = df$slip_214[i]) +
      ylim(-6, 5) +
      ylab(gp) +
      theme(panel.grid = element_blank(),
            axis.text.x = element_text(size = 10, color = "black"),
            axis.text.y = element_text(size = 10, color = "black"),
            axis.title.x = element_text(size = 13),
            axis.title.y = element_text(size = 13)) +
      scale_color_manual(values = color_tcr)
    print(p)
  }
  
}


# Fig 6C: Clonotype projection, d8
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
  
  plot_list <- list()
  for(chain in levels(dft_d8$chain)){
    data <- data.frame(seu@reductions$UMAP_LSI_Combined_batch3.3@cell.embeddings) 
    data$barcode <- rownames(data)
    data <- cbind(data, seu@meta.data)
    
    xrange <- max(data$UMAPCombinedbatch33_1) - min(data$UMAPCombinedbatch33_1)
    yrange <- max(data$UMAPCombinedbatch33_2) - min(data$UMAPCombinedbatch33_2)
    p <- ggplot(data) +
      geom_scattermore(aes(x = UMAPCombinedbatch33_1, y = UMAPCombinedbatch33_2), color = "grey90") +
      theme_void() +
      coord_fixed(xrange/yrange) +
      geom_point(data = data[which(data$timepoint == "d8" & (data$chain1 ==  chain | data$chain2 == chain)),], aes(x = UMAPCombinedbatch33_1, y = UMAPCombinedbatch33_2, fill = cell_module), 
                 shape = 21, size = 3, show.legend = FALSE) +
      scale_fill_manual(values = mmjoint_colors) +
      ggtitle(chain)
    plot_list[[chain]] <- p
  }
  p2 <- list(plot_list[["CAASGTSGSWQLIF++CAWSLGVANTGQLYF"]],
             plot_list[["CALGEEGGRALIF++CASSPQGRTGQLYF"]],
             plot_list[["CALGEASSGQKLVF++CAWSLGLGNYAEQFF"]],
             plot_list[["CALGGTGGYKVVF++CASSPRDSDYTF"]],
             plot_list[["CAPGNSNNRIFF++CASRGYTEVFF"]])
  wrap_plots(p2, ncol = 5)
}

# Fig 6D: Clonotype projection, d14
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
  
  plot_list <- list()
  for(chain in levels(dft_d14$chain)){
    data <- data.frame(seu@reductions$UMAP_LSI_Combined_batch3.3@cell.embeddings) 
    data$barcode <- rownames(data)
    data <- cbind(data, seu@meta.data)
    
    xrange <- max(data$UMAPCombinedbatch33_1) - min(data$UMAPCombinedbatch33_1)
    yrange <- max(data$UMAPCombinedbatch33_2) - min(data$UMAPCombinedbatch33_2)
    p <- ggplot(data) +
      geom_scattermore(aes(x = UMAPCombinedbatch33_1, y = UMAPCombinedbatch33_2), color = "grey90") +
      theme_void() +
      coord_fixed(xrange/yrange) +
      geom_point(data = data[which(data$timepoint == "d14" & (data$chain1 ==  chain | data$chain2 == chain)),], aes(x = UMAPCombinedbatch33_1, y = UMAPCombinedbatch33_2, fill = cell_module), 
                 shape = 21, size = 3, show.legend = FALSE) +
      scale_fill_manual(values = mmjoint_colors) +
      ggtitle(chain)
    plot_list[[chain]] <- p
  }
  p2 <- list(plot_list[["CAADSSSGSWQLIF++CASSQDPTGGKSQNTLYF"]],
             plot_list[["CALWELESNNRIFF++CASSDSSNERLFF"]],
             plot_list[["CAANSGTYQRF++CASSRDWGLLSQNTLYF"]],
             plot_list[["CALGGTGGYKVVF++CASSPRDSDYTF"]],
             plot_list[["CAVSMNSNNRIFF++CASGDAGAPLF"]])
  wrap_plots(p2, ncol = 5)
}


# Fig 6E: d8 cell module composition of each TCR clone
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
  write_tsv(dft_d8, "source_fig_6E.tsv", na = "")
  
  
}

# Fig 6F: d14 cell module composition of each TCR clone
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
  write_tsv(dft_d14, "source_fig_6F.tsv", na = "")
}



# Fig 6G: chromatin checkpoint score of d8 clones
{
  genes_predi_score <- list("checkpoint_1" = c("Dock9", "Ly6e", "Ptpn13", "Zfp281"),
                            # "checkpoint_2" = c("Ccl5", "Ccr2", "Ccr5", "Cx3cr1", "Fasl", "Gsap", "Gzma", "Il18r1", "Il18rap", "Itga4", "Itgb1", "Klrc1", "Klrk1", "Lamc1", "Sntb2", "Tbx21", "Tcf3", "Ttc39c", "Zeb2"),
                            "checkpoint_2" = c("Cyba", "Cyth4", "Dock5", "Fkbp5", "Itga4", "Itgb1", "Lamc1"),
                            "checkpoint_3" = c("Bcl6", "Chd3", "Dennd4a", "Etf1", "Etv3", "Ifrd1", "Junb", "Vezf1"),
                            "checkpoint_4" = c("Iigp1", "Ly6a", "Nap1l1", "Ncapd2", "Nedd4l", "Pcgf5"))
  
  genes_maint_score <- list("checkpoint_1" = c("Ccm2", "Ccr7", "Ccr9", "Dcaf17", 
                                               "Hspbap1", "Hvcn1", "Il4ra", "Il6ra", "Itga6", "Mettl8", "Skp1a", "Smc4",  "St6gal1", 
                                               "Tet1", "Trib2", "Trio", "Ttc28", "Vamp1", "Zdhhc14"),
                            "checkpoint_2"  = c("Arhgap25", "Batf", "Cers6", "Fchsd2", "Ggact", "Klf3", "Naa16", "Spag9"),
                            "checkpoint_3" = c("Anxa2", "Gzmb", "Gzmk", "Mki67", "Nkg7", "Ssrp1", "Ubb"),
                            "checkpoint_4"  = c("Gnptab", "Ikzf3", "N4bp1", "Osbpl3"))
  
  g1 <- lapply(seq_along(genes_predi_score), function(i){
    out <- tibble(checkpoint = names(genes_predi_score)[i], type = "predisposition", genes = genes_predi_score[[i]])
    out$group <- paste0(out$checkpoint, "_", out$type)
    return(out)
  }) %>% Reduce("rbind", .)
  g2 <- lapply(seq_along(genes_maint_score), function(i){
    out <- tibble(checkpoint = names(genes_maint_score)[i], type = "maintenance", genes = genes_maint_score[[i]])
    out$group <- paste0(out$checkpoint, "_", out$type)
    return(out)
  }) %>% Reduce("rbind", .)
  g3 <- lapply(seq_along(genes_maint_score), function(i){
    out <- tibble(checkpoint = names(genes_maint_score)[i], type = "all", genes = c(genes_predi_score[[i]], genes_maint_score[[i]]))
    out$group <- paste0(out$checkpoint, "_", out$type)
    return(out)
  }) %>% Reduce("rbind", .)
  genes_chrom_potential <- rbind(g1, g2, g3)
  
  fname <- c("checkpoint_1_predisposition", "checkpoint_1_maintenance", "checkpoint_1_all", 
             "checkpoint_2_predisposition", "checkpoint_2_maintenance", "checkpoint_2_all", 
             "checkpoint_3_predisposition", "checkpoint_3_maintenance", "checkpoint_3_all", 
             "checkpoint_4_predisposition", "checkpoint_4_maintenance", "checkpoint_4_all")
  
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
  
  pdata <- mclapply(seq_along(fname), function(n){
    gp <- fname[n]
    mean_gp <- mean_checkpoint[gp]
    sd_gp <- sd_checkpoint[gp]
    genes <- dplyr::filter(genes_chrom_potential, group == gp) %>% pull(genes)
    tcr_data <- lapply(seq_along(d8_tcr), function(i){
      tcr_i <- d8_tcr[i]
      tsub <- tcr_paired %>% dplyr::filter(timepoint == "d8") %>% dplyr::filter(CDRa_CDRb == tcr_i)
      bc_i <- tsub$barcode
      residual <- mean(colMeans(dorcMat_sm_norm[genes, bc_i, drop = FALSE] - rnaMat_sm_norm[genes, bc_i, drop = FALSE]))
      (residual-mean_gp)/sd_gp
    }) %>% unlist()
    tibble(gp = tcr_data)
  }, mc.cores = 8) %>% Reduce("cbind", .)
  colnames(pdata) <- fname
  pdata <- t(pdata)
  colnames(pdata) <- d8_tcr
  
  pheatmap::pheatmap(pdata[c("checkpoint_1_predisposition", "checkpoint_1_maintenance",
                             "checkpoint_2_predisposition", "checkpoint_2_maintenance",
                             "checkpoint_3_predisposition", "checkpoint_3_maintenance", 
                             "checkpoint_4_predisposition", "checkpoint_4_maintenance"),], 
                     border_color = NA,
                     show_colnames = TRUE,
                     show_rownames = TRUE,
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     breaks = seq(-2, 2, length.out = 91),
                     clustering_method = "ward.D",
                     scale = "none",
                     fontsize_row = 10,
                     fontsize_col = 10,
                     angle_col = 45,
                     # color = choose_colorset("wolfgang_extra", 101),
                     color = rev(choose_colorset("PuOr", 101))[5:95],
                     # filename = "d8_tcr_desc_sorting_checkpoint_all_heatmap_cut.pdf",
                     width = 13, height = 3.5)
  pdata <- data.frame(pdata) %>% 
    tibble::rownames_to_column(var = "checkpoint")
  write.table(pdata, "source_fig_6G.tsv", sep = "\t", quote = FALSE, na = "", row.names = FALSE, col.names = TRUE)
}



# Fig 6H: chromatin checkpoint score of d14 clones
{
  # load(paste0(RDSdir, "chromatin_potential_gene_module.Rdata"))
  load(paste0(RDSdir, "lm_merged_overlapped_harmony_mat_smoothed.RData"))
  load(paste0(RDSdir, "checkpoint_mean_sd_manual_2.RData"))
  
  genes_predi_score <- list("checkpoint_1" = c("Dock9", "Ly6e", "Ptpn13", "Zfp281"),
                            # "checkpoint_2" = c("Ccl5", "Ccr2", "Ccr5", "Cx3cr1", "Fasl", "Gsap", "Gzma", "Il18r1", "Il18rap", "Itga4", "Itgb1", "Klrc1", "Klrk1", "Lamc1", "Sntb2", "Tbx21", "Tcf3", "Ttc39c", "Zeb2"),
                            "checkpoint_2" = c("Cyba", "Cyth4", "Dock5", "Fkbp5", "Itga4", "Itgb1", "Lamc1"),
                            "checkpoint_3" = c("Bcl6", "Chd3", "Dennd4a", "Etf1", "Etv3", "Ifrd1", "Junb", "Vezf1"),
                            "checkpoint_4" = c("Iigp1", "Ly6a", "Nap1l1", "Ncapd2", "Nedd4l", "Pcgf5"))
  
  genes_maint_score <- list("checkpoint_1" = c("Ccm2", "Ccr7", "Ccr9", "Dcaf17", 
                                               "Hspbap1", "Hvcn1", "Il4ra", "Il6ra", "Itga6", "Mettl8", "Skp1a", "Smc4",  "St6gal1", 
                                               "Tet1", "Trib2", "Trio", "Ttc28", "Vamp1", "Zdhhc14"),
                            "checkpoint_2"  = c("Arhgap25", "Batf", "Cers6", "Fchsd2", "Ggact", "Klf3", "Naa16", "Spag9"),
                            "checkpoint_3" = c("Anxa2", "Gzmb", "Gzmk", "Mki67", "Nkg7", "Ssrp1", "Ubb"),
                            "checkpoint_4"  = c("Gnptab", "Ikzf3", "N4bp1", "Osbpl3"))
  
  g1 <- lapply(seq_along(genes_predi_score), function(i){
    out <- tibble(checkpoint = names(genes_predi_score)[i], type = "predisposition", genes = genes_predi_score[[i]])
    out$group <- paste0(out$checkpoint, "_", out$type)
    return(out)
  }) %>% Reduce("rbind", .)
  g2 <- lapply(seq_along(genes_maint_score), function(i){
    out <- tibble(checkpoint = names(genes_maint_score)[i], type = "maintenance", genes = genes_maint_score[[i]])
    out$group <- paste0(out$checkpoint, "_", out$type)
    return(out)
  }) %>% Reduce("rbind", .)
  g3 <- lapply(seq_along(genes_maint_score), function(i){
    out <- tibble(checkpoint = names(genes_maint_score)[i], type = "all", genes = c(genes_predi_score[[i]], genes_maint_score[[i]]))
    out$group <- paste0(out$checkpoint, "_", out$type)
    return(out)
  }) %>% Reduce("rbind", .)
  genes_chrom_potential <- rbind(g1, g2, g3)
  
  
  fname <- c("checkpoint_1_predisposition", "checkpoint_1_maintenance", "checkpoint_1_all", 
             "checkpoint_2_predisposition", "checkpoint_2_maintenance", "checkpoint_2_all", 
             "checkpoint_3_predisposition", "checkpoint_3_maintenance", "checkpoint_3_all", 
             "checkpoint_4_predisposition", "checkpoint_4_maintenance", "checkpoint_4_all")
  
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
    dplyr::mutate(major_module = str_extract(cell_module, "(.*)_", group = 1)) %>%
    dplyr::mutate_at(vars(major_module), ~factor(.x, levels = c("Quiescent", "Priming", "Activation", "Effector", "Effector-memory_transition", "Memory"))) %>%
    dplyr::mutate(ratio_Quiescent = sum(major_module %in% c("Quiescent")) / n()) %>%
    dplyr::mutate(ratio_Priming = sum(major_module %in% c("Priming")) / n()) %>%
    dplyr::mutate(ratio_Activation = sum(major_module %in% c("Activation")) / n()) %>%
    dplyr::mutate(ratio_Effector = sum(major_module %in% c("Effector")) / n()) %>%
    dplyr::mutate(ratio_EMT = sum(major_module %in% c("Effector-memory_transition")) / n()) %>%
    dplyr::mutate(ratio_Memory = sum(major_module %in% c("Memory")) / n()) %>%
    dplyr::mutate(poly = n_distinct(major_module)) %>%
    dplyr::mutate(MEratio = ratio_Effector/(ratio_EMT + ratio_Effector)) %>%
    arrange(desc(ratio_Quiescent), desc(ratio_Priming), desc(ratio_Activation), desc(ratio_Effector), desc(ratio_EMT), desc(ratio_Memory)) 
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
  
  fname <- c("checkpoint_1_predisposition", "checkpoint_1_maintenance",
             "checkpoint_2_predisposition", "checkpoint_2_maintenance",
             "checkpoint_3_predisposition", "checkpoint_3_maintenance",
             "checkpoint_4_predisposition", "checkpoint_4_maintenance",
             "checkpoint_1_all", "checkpoint_2_all", "checkpoint_3_all", "checkpoint_4_all")
  
  pdata <- mclapply(seq_along(fname), function(n){
    gp <- fname[n]
    mean_gp <- mean_checkpoint[gp]
    sd_gp <- sd_checkpoint[gp]
    genes <- dplyr::filter(genes_chrom_potential, group == gp) %>% pull(genes)
    tcr_data <- lapply(seq_along(d14_tcr), function(i){
      tcr_i <- d14_tcr[i]
      tsub <- tcr_paired %>% dplyr::filter(timepoint == "d14") %>% dplyr::filter(CDRa_CDRb == tcr_i)
      bc_i <- tsub$barcode
      residual <- mean(colMeans(dorcMat_sm_norm[genes, bc_i, drop = FALSE] - rnaMat_sm_norm[genes, bc_i, drop = FALSE]))
      (residual-mean_gp)/sd_gp
    }) %>% unlist()
    tibble(gp = tcr_data)
  }, mc.cores = 8) %>% Reduce("cbind", .)
  colnames(pdata) <- fname
  pdata <- t(pdata)
  colnames(pdata) <- d14_tcr
  score_d14 <- as.data.frame.array(t(pdata))
  score_d14$timepoint <- "d14"
  score_d14$MEratio <- distinct(dft_d14 %>% dplyr::select(chain, MEratio))$MEratio
  
  # d14 double gradient heatmap, cut to min-max
  pheatmap::pheatmap(pdata[c("checkpoint_1_predisposition", "checkpoint_1_maintenance",
                             "checkpoint_2_predisposition", "checkpoint_2_maintenance",
                             "checkpoint_3_predisposition", "checkpoint_3_maintenance", 
                             "checkpoint_4_predisposition", "checkpoint_4_maintenance"),], 
                     border_color = NA,
                     show_colnames = TRUE,
                     show_rownames = TRUE,
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     breaks = seq(-2, 2, length.out = 91),
                     clustering_method = "ward.D",
                     scale = "none",
                     fontsize_row = 10,
                     fontsize_col = 10,
                     angle_col = 45,
                     # color = choose_colorset("wolfgang_extra", 101),
                     color = rev(choose_colorset("PuOr", 101))[5:95],
                     # filename = "d14_tcr_desc_sorting_checkpoint_manual_2_heatmap.pdf",
                     width = 10.6, height = 4.5)
  pdata <- data.frame(pdata) %>% 
    tibble::rownames_to_column(var = "checkpoint")
  write.table(pdata, "source_fig_6H.tsv", sep = "\t", quote = FALSE, na = "", row.names = FALSE, col.names = TRUE)
}


# 6I 6J 6K
{
  data <- read_tsv("source_fig_6I.tsv")
  
  # 6I
  ggplot(data, aes(x = checkpoint_3_predisposition, y = checkpoint_3_maintenance, fill = potential)) +
    geom_point(size = 3, shape = 21) +
    scale_fill_gradientn(colors = choose_colorset("green_blue")) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(size = 10, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.title.x = element_text(size = 13),
          axis.title.y = element_text(size = 13),
          strip.background = element_rect(fill = "grey90"),
          strip.text = element_text(size = 10, color = "black")) +
    guides(fill = guide_colourbar(ticks = FALSE)) +
    coord_fixed(ratio = (max(data$checkpoint_3_predisposition) - min(data$checkpoint_3_predisposition))/
                  (max(data$checkpoint_3_maintenance) - min(data$checkpoint_3_maintenance)))
  
  # 6J
  ggplot(data, aes(x = checkpoint_3_predisposition, y = checkpoint_3_maintenance, fill = checkpoint_3_all)) +
    geom_point(size = 3, shape = 21) +
    scale_fill_gradientn(colors = choose_colorset("solar_extra")) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(size = 10, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.title.x = element_text(size = 13),
          axis.title.y = element_text(size = 13),
          strip.background = element_rect(fill = "grey90"),
          strip.text = element_text(size = 10, color = "black")) +
    guides(fill = guide_colourbar(ticks = FALSE, title = "potential")) +
    coord_fixed(ratio = (max(data$checkpoint_3_predisposition) - min(data$checkpoint_3_predisposition))/
                  (max(data$checkpoint_3_maintenance) - min(data$checkpoint_3_maintenance)))
  
  # 6K
  model <- lm(potential ~ checkpoint_3_all, data = data)
  s <- summary(model)
  sqrt(s$r.squared)
  ggplot(data, aes(x = checkpoint_3_all, y = potential)) +
    geom_smooth(se = FALSE, formula = y ~ x, method = lm, fullrange = FALSE) +
    geom_point(aes(fill = checkpoint_3_all), size = 3, shape = 21) +
    scale_fill_gradientn(colors = choose_colorset("solar_extra")) +
    ylim(0, 1) +
    # scale_fill_gradientn(colors = choose_colorset("nostelgia")) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(size = 10, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.title.x = element_text(size = 13),
          axis.title.y = element_text(size = 13),
          strip.background = element_rect(fill = "grey90"),
          strip.text = element_text(size = 10, color = "black")) +
    guides(fill = guide_colourbar(ticks = FALSE, title = "potential")) +
    coord_fixed(ratio = (max(data$checkpoint_3_all) - min(data$checkpoint_3_all))/
                  (max(data$potential) - min(data$potential))) +
    annotate(geom = "text", x = -1, y = 0.9, colour = "black", size = 5, angle = 0, label = paste0("R = 0.84"))
}




