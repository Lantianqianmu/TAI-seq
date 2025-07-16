# figure 6
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
library(igraph)
library(ggraph)
library(ggsignif)
options(repos="http://mirrors.tuna.tsinghua.edu.cn/CRAN")
options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor")
source("/home/zeemeeuw/YangLab/JoINT-seq/QC/code/utils.R")

setwd("/home/zeemeeuw/YangLab/JoINT-seq/Manuscript/fig_data/fig_5_data")

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

# Fig 6A: gliph motif network
{
  cdr3b_highlight <- str_extract(highlight, "\\+\\+(.*)", group = 1)
  
  # merge with scTCR
  data <- read_csv("mmJoint_GLIPH2_merge_cluster.csv") %>% 
    dplyr::select(index, pattern, Fisher_score, TcRb, V, J, TcRa, Sample, Freq) %>% 
    rename_with(~c("index", "pattern", "Fisher_score", "CDR3b", "Vb", "Jb", "CDR3a", "Sample", "Freq")) %>% 
    dplyr::mutate(tissue = str_extract(Sample, "(.*):(.*)", group = 1),
                  timepoint = str_extract(Sample, "(.*):(.*)", group = 2)) %>% 
    filter(pattern != "single") %>% 
    unite("CDR3a_CDR3b", CDR3a, CDR3b, sep = "__", remove = FALSE) %>% 
    dplyr::mutate(highlight_tcr = CDR3a_CDR3b %in% highlight) %>% 
    group_by(pattern) %>% 
    dplyr::mutate(highlight_motif = sum(highlight_tcr) > 0) %>% 
    filter(index < 500| highlight_motif)
  
  
  conn <- data %>% 
    dplyr::select(CDR3b, pattern, Freq) %>% 
    group_by(pattern, CDR3b) %>% 
    dplyr::mutate(counts = sum(Freq), value = 1) %>% 
    arrange(desc(counts)) %>% 
    slice_head(n = 1)
  
  vertices_CDR3 <- conn %>%
    ungroup() %>% 
    dplyr::select(CDR3b, counts) %>% 
    group_by(CDR3b) %>% 
    dplyr::summarise(n = sum(counts)) %>% 
    rename_with(~c("name", "n")) %>% 
    dplyr::mutate(group = "CDR3β")
  vertices_motif <- conn %>%
    ungroup() %>% 
    dplyr::select(pattern, counts) %>% 
    group_by(pattern) %>% 
    dplyr::summarise(n = sum(counts)) %>% 
    rename_with(~c("name", "n")) %>% 
    dplyr::mutate(group = "Motif")
  vertices <- bind_rows(vertices_CDR3, vertices_motif)
  vertices <- vertices %>% dplyr::mutate(group = if_else(name %in% cdr3b_highlight, "highlight", group))
  vertices <- vertices %>% dplyr::mutate(n = if_else(n > 100, 100, n))
  vertices <- dplyr::bind_rows(vertices[vertices$group != "highlight",], vertices[vertices$group == "highlight",])
  
  conn <- conn %>% dplyr::select(CDR3b, pattern, value) %>% 
    rename_with(~c("from", "to", "value"))
  gr <- graph_from_data_frame(conn, vertices = vertices)
  
  ggraph(gr, layout = "igraph", algorithm = 'mds') +
    geom_edge_link(edge_color = "#FBE6BF", edge_alpha = 1, edge_width = 0.5, angle_calc = "along", linejoin = "round", check_overlap = TRUE, show.legend = FALSE) + 
    geom_node_point(aes(fill = group, size = n), stroke = 0.1, shape = 21, alpha = 1, position = "jitter", show.legend = FALSE) +
    scale_size_continuous(range = c(1, 5)) +
    geom_node_text(aes(label = ifelse(name %in% cdr3b_highlight, as.character(name), "")), size = 4, color = "black") +
    theme_void() +
    scale_fill_manual(values = c(CDR3β = "#bcbddd", "Motif" = "#FFE060", highlight = "#DC050A")) +
    guides(fill = guide_legend(override.aes = list(size = 5))) +
    expand_limits(x = c(-3, 3), y = c(-3, 3)) +
    coord_equal()
}
