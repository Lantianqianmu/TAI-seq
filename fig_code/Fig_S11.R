#----------------------------- sc96 -----------------------------#
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(ggplot2)
library(Biostrings)
library(ggpubr)
library(ggforce)
library(patchwork)
library(pheatmap)
options(repos="http://mirrors.tuna.tsinghua.edu.cn/CRAN")
options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor")
source("/home/zeemeeuw/YangLab/JoINT-seq/QC/code/utils.R")


setwd("/home/zeemeeuw/YangLab/JoINT-seq/Manuscript/fig_data/fig_S11_data")


# S11A
{
  tcr_paired <- readxl::read_excel("Lm-TCR4.57-clonotypes-paired in rna.xlsx") %>% 
    dplyr::mutate(tissue = str_extract(Tissue_time, "(.*)-", group = 1)) %>% 
    group_by(CDRa_CDRb, timepoint, tissue, barcode) %>% 
    distinct(CDRa_CDRb, timepoint, tissue, barcode, .keep_all = TRUE) %>% 
    ungroup()
  
  pdata <- as.data.frame.array(table(tcr_paired$TRAV, tcr_paired$timepoint))
  pdata <- pdata[,c("WT", "d3", "d5", "d8", "d14", "d30", "R5", "Y5")]
  pdata <- t(t(pdata)/colSums(pdata))
  pdata[pdata > 0.2] <- 0.2
  pheatmap::pheatmap(pdata,
                     show_colnames = TRUE,
                     show_rownames = TRUE,
                     cluster_rows = TRUE,
                     cluster_cols = FALSE,
                     treeheight_col = 0,
                     treeheight_row = 10,
                     border_color = NA,
                     clustering_method = "centroid",
                     scale = "none",
                     color = choose_colorset("wolfgang_basic"),
                     # filename = "TRAV_gene_heatmap.pdf",
                     fontsize_row = 8,
                     width = 3, height = 10)
}


# S11C
{
  tcr_paired <- readxl::read_excel("mmJoint96_filtered_tcr_longer.xlsx") %>% 
    dplyr::mutate(TRBV = str_extract(TRBV, "(.*)\\*", group = 1),
                  TRBJ = str_extract(TRBJ, "(.*)\\*", group = 1))
  
  pdata <- as.data.frame.array(table(tcr_paired$TRBV, tcr_paired$timepoint))
  pdata <- pdata[, 
                 c("D5", "D8", "D14", "D21", "R5")]
  pdata <- t(t(pdata)/colSums(pdata))
  pdata[pdata > 0.3] <- 0.3
  
  pheatmap::pheatmap(pdata,
                     show_colnames = TRUE,
                     show_rownames = TRUE,
                     cluster_rows = TRUE,
                     cluster_cols = FALSE,
                     treeheight_col = 0,
                     treeheight_row = 10,
                     border_color = NA,
                     clustering_method = "average",
                     scale = "none",
                     color = choose_colorset("wolfgang_basic"),
                     # filename = "sc96_TRBV_gene_heatmap.pdf",
                     fontsize_row = 8,
                     width = 2.08, height = 3.4)
}

# S11E-F
{
  # ----------------------------- mapping 213-230 to pie plots -----------------------------#
  
  highlight <- c("213" = "CAVRYPRNNNNAPRF++CASSDSSNERLFF",
                 "214" = "CALWELESNNRIFF++CASSDSSNERLFF",
                 "215" = "CAASEMAGGSNAKLTF++CASSGASAETLYF",
                 "216" = "CVVGAGSNYQLIW++CASSGASAETLYF",
                 "217" = "CALGEASSGQKLVF++CAWSLGLGNYAEQFF",
                 "218" = "CALGGTGGYKVVF++CASSPRDSDYTF",
                 "219" = "CAASGYNQGKLIF++CASRQGARDTQYF",
                 "220" = "CAASLFYQGGRALIF++CASSQDWGSSAETLYF",
                 "221" = "CALRTGGYKVVF++CASSQRGSQNTLYF",
                 "222" = "CAASGTGNTGKLIF++CAWIPTVANTGQLYF",
                 "223" = "CAIPNNYAQGLTF++CASRDTGQLYF",
                 "224" = "CAIPNNYAQGLTF++CAWQTFNNQAPLF",
                 "225" = "CVVGAGSNYQLIW++CGAKNRDRWEVFF",
                 "226" = "CAASEMAGGSNAKLTF++CGAKNRDRWEVFF",
                 "227" = "CAANDGSSGNKLIF++CASSRANSDYTF",
                 "228" = "CAVRDSNYQLIW++CASSFSLNYAEQFF",
                 "229" = "CIVTDIRTGGYKVVF++CGARDQGNNNQAPLF",
                 "230" = "CAASFSGNEKITF++CTCSAGGAGNTLYF")
  
  
  tcr_paired <- readxl::read_excel("Lm-TCR4.57-clonotypes-paired in rna.xlsx")
  
  tcr_paired_pie <- tcr_paired %>% 
    dplyr::mutate(tissue = str_extract(Tissue_time, "(.*)-", group = 1)) %>% 
    group_by(CDRa_CDRb, timepoint, tissue, barcode) %>% 
    distinct(CDRa_CDRb, timepoint, tissue, barcode, .keep_all = TRUE) %>% 
    group_by(CDRa_CDRb, timepoint, tissue) %>% 
    dplyr::mutate(counts = n()) %>% 
    dplyr::filter(counts > 1) %>% 
    group_by(timepoint, tissue) %>% 
    dplyr::mutate(prop = counts/n()) %>% 
    arrange(desc(prop)) %>% 
    ungroup() %>% 
    dplyr::select(CDRa_CDRb, timepoint, tissue, counts, prop) %>% 
    distinct()
  
  
  
  for(num in names(highlight)){
    plots <- list()
    n <- 0
    for(time in c("d3", "d5", "d8", "d14", "d30", "R5", "Y5")){
      for(ti in c("Ln", "Spleen", "Liver", "Blood")){
        pdata <- tcr_paired_pie %>% 
          dplyr::filter(timepoint == time & tissue == ti) %>% 
          dplyr::mutate(highlight = if_else(CDRa_CDRb == highlight[num], 1, 0)) %>% 
          arrange(desc(highlight), desc(prop)) %>% 
          dplyr::mutate(
            start = 2 * pi - cumsum(lag(prop, default = 0)) * 2 * pi, 
            end = 2 * pi - cumsum(prop) * 2 * pi                       
          ) 
        
        p <- ggplot() +
          ggforce::stat_arc_bar(data = pdata, aes(x0 = 0, y0 = 0, r0 = 0, r = 1, fill = prop, start = start, end = end), 
                                show.legend = FALSE, color = NA, linewidth = 0) +
          scale_fill_gradientn(colors = colorRampPalette(c("#F2F4FE", "#B6BFEB"))(1000), limits = c(0, 1)) +
          theme_void() +
          coord_fixed()
        if(max(pdata$highlight == 1)){
          p <- p + ggforce::stat_arc_bar(data = pdata %>% dplyr::filter(highlight == 1), aes(x0 = 0, y0 = 0, r0 = 0, r = 1, start = start, end = end),  fill = "#DC050A",
                                         show.legend = FALSE, color = NA, linewidth = 0)
        }
        n <- n + 1
        plots[[n]] <- p
      }
    }
    
    pp <- wrap_plots(plots, ncol = 4, nrow = 8)
    print(pp)
    ggsave(filename = paste0("pie_chart_", num, ".pdf"), pp, width = 4, height = 8)
  }
  
  
  
  
  
  
  
  
  
  # 227-230
  data_longer <- read_excel("mmJoint96_filtered_tcr_longer.xlsx") 
  
  highlight <- c("227" = "CAANDGSSGNKLIF__CASSRANSDYTF",
                 "228" = "CAVRDSNYQLIW__CASSFSLNYAEQFF",
                 "229" = "CIVTDIRTGGYKVVF__CGARDQGNNNQAPLF",
                 "230" = "CAASFSGNEKITF__CTCSAGGAGNTLYF")
  
  tcr_paired_pie <- data_longer %>% 
    group_by(CDR3a_CDR3b, timepoint) %>% 
    dplyr::mutate(counts = n()) %>% 
    dplyr::filter(counts > 1) %>% 
    group_by(timepoint) %>% 
    dplyr::mutate(prop = counts/n()) %>% 
    arrange(desc(prop)) %>% 
    ungroup() %>% 
    dplyr::select(CDR3a_CDR3b, timepoint, counts, prop) %>% 
    distinct()
  
  
  
  for(num in names(highlight)){
    plots <- list()
    n <- 0
    for(time in c("D5", "D8", "D14", "D21", "R5")){
      pdata <- tcr_paired_pie %>% 
        dplyr::filter(timepoint == time) %>% 
        dplyr::mutate(highlight = if_else(CDR3a_CDR3b == highlight[num], 1, 0)) %>% 
        arrange(desc(highlight), desc(prop)) %>% 
        dplyr::mutate(
          start = 2 * pi - cumsum(lag(prop, default = 0)) * 2 * pi, 
          end = 2 * pi - cumsum(prop) * 2 * pi                       
        ) 
      
      p <- ggplot() +
        ggforce::stat_arc_bar(data = pdata, aes(x0 = 0, y0 = 0, r0 = 0, r = 1, fill = prop, start = start, end = end), 
                              show.legend = FALSE, color = NA, linewidth = 0) +
        scale_fill_gradientn(colors = colorRampPalette(c("#F2F4FE", "#B6BFEB"))(1000), limits = c(0, 1)) +
        theme_void() +
        coord_fixed()
      if(max(pdata$highlight == 1)){
        p <- p + ggforce::stat_arc_bar(data = pdata %>% dplyr::filter(highlight == 1), aes(x0 = 0, y0 = 0, r0 = 0, r = 1, start = start, end = end),  fill = "#DC050A",
                                       show.legend = FALSE, color = NA, linewidth = 0)
      }
      n <- n + 1
      plots[[n]] <- p
      
    }
    
    pp <- wrap_plots(plots, ncol = 5)
    print(pp)
    ggsave(filename = paste0("pie_chart_", num, ".pdf"), pp, width = 4, height = 8)
  }
  
  
}














