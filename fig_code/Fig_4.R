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
options(repos="http://mirrors.tuna.tsinghua.edu.cn/CRAN")
options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor")
source("/home/zeemeeuw/YangLab/JoINT-seq/QC/code/utils.R")

setwd("/home/zeemeeuw/YangLab/JoINT-seq/Manuscript/fig_data/fig_4_data")

# Fig 4A: pie plot of TCR expanstion
{
  tcr_paired <- readxl::read_excel("Lm-TCR4.57-clonotypes-paired in rna.xlsx")
  
  tcr_paired_pie <- tcr_paired %>% 
    dplyr::mutate(tissue = str_extract(Tissue_time, "(.*)-", group = 1)) %>% 
    group_by(CDRa_CDRb, timepoint, tissue, barcode) %>% 
    distinct(CDRa_CDRb, timepoint, tissue, barcode, .keep_all = TRUE) %>% 
    group_by(CDRa_CDRb, timepoint, tissue) %>% 
    dplyr::mutate(counts = n()) %>% 
    group_by(timepoint, tissue) %>% 
    dplyr::mutate(prop = counts/n()) %>% 
    arrange(desc(prop)) %>% 
    ungroup() %>% 
    dplyr::select(CDRa_CDRb, timepoint, tissue, counts, prop) %>% 
    distinct()
  
  generate_color_grad <- function(gap = c(0.1, 0.02), maximum = 0.5){
    c1 <- colorRampPalette(c("#E9BFCA", "#F7B4B4", "#ED8080", "#DF4A4A", "#BE2A2A"))(160)
    c2 <- colorRampPalette(rev(c("#DCCAE0", '#CFD5F7', "#B6BFEB", "#8694D8")))(36)      
    c3 <- colorRampPalette(c("#f0f0f0", "#DEFADC","#BEEBD9", "#ABCED8", "#98B1D8"))(4)
    return(c(c3, c2, c1))
  }
  
  plots <- list()
  n <- 0
  for(time in c("WT", "d3", "d5", "d8", "d14", "d30", "R5", "Y5")){
    for(ti in c("Ln", "Spleen", "Liver", "Blood")){
      pdata <- tcr_paired_pie %>% dplyr::filter(timepoint == time & tissue == ti) %>% 
        dplyr::mutate(fill = if_else(counts == 1, 0, prop))
      pdata <- pdata %>% dplyr::mutate(
        start = 2 * pi - cumsum(lag(prop, default = 0)) * 2 * pi, 
        end = 2 * pi - cumsum(prop) * 2 * pi                       
      )
      pdata_sum <- pdata %>% dplyr::mutate(expansion_level = case_when(prop >= 0.1 ~ "high", prop < 0.1 & prop >= 0.01 ~ "medium", prop < 0.01 & counts > 1 ~ "low", .default = "no expansion")) %>% 
        group_by(expansion_level) %>% 
        dplyr::summarise(prop_sum = sum(prop))
      pdata_sum$expansion_level <- factor(pdata_sum$expansion_level, levels = c("high", "medium", "low", "no expansion"))
      pdata_sum <- pdata_sum %>% 
        arrange(expansion_level) %>%
        dplyr::mutate(
          start = 2 * pi - cumsum(lag(prop_sum, default = 0)) * 2 * pi, 
          end = 2 * pi - cumsum(prop_sum) * 2 * pi                       
        )
      
      pdata_sum <- pdata_sum 
      
      p <- ggplot() +
        ggforce::stat_arc_bar(data = pdata, aes(x0 = 0, y0 = 0, r0 = 0, r = 0.7, fill = fill, start = start, end = end), 
                              show.legend = FALSE, color = NA, linewidth = 0) +
        scale_fill_gradientn(colors = generate_color_grad(), limits = c(0, 0.5)) +
        new_scale_fill() +
        ggforce::stat_arc_bar(data = pdata_sum, aes(x0 = 0, y0 = 0, r0 = 0.8, r = 1, fill = expansion_level, start = start, end = end), 
                              show.legend = FALSE,  color = NA, linewidth = 0) +
        # scale_fill_manual(values = c("high" = "#3361A5", "medium" = "#248AF3", "low" = "#A5C8F0", "no expansion" = "#f0f0f0")) +
        # scale_fill_manual(values = c("high" = '#C61A2C', "medium" = "#FA8E24", "low" = "#FFEC77", "no expansion" = "#f0f0f0")) +
        # scale_fill_manual(values = c("high" =  "#ABDEB6", "medium" = "#E3F4B1", "low" = "#FCFED3", "no expansion" = "#f0f0f0")) +
        scale_fill_manual(values = c("high" =  "#7AB5F6", "medium" = "#a0c9e5", low = '#D1E5E7', "no expansion" = "#f0f0f0")) +
        theme_void() +
        coord_fixed()
      n <- n + 1
      plots[[n]] <- p
    }
  }
  pp <- wrap_plots(plots, ncol = 4, nrow = 8)
  print(pp)
  
  
  # generate source data
  pdata_list <- list()
  pdata_sum_list <- list()
  for(time in c("WT", "d3", "d5", "d8", "d14", "d30", "R5", "Y5")){
    for(ti in c("Ln", "Spleen", "Liver", "Blood")){
      pdata <- tcr_paired_pie %>% dplyr::filter(timepoint == time & tissue == ti) %>% 
        dplyr::mutate(fill = if_else(counts == 1, 0, prop))
      pdata <- pdata %>% dplyr::mutate(
        start = 2 * pi - cumsum(lag(prop, default = 0)) * 2 * pi, 
        end = 2 * pi - cumsum(prop) * 2 * pi                       
      )
      pdata_sum <- pdata %>% dplyr::mutate(expansion_level = case_when(prop >= 0.1 ~ "high", prop < 0.1 & prop >= 0.01 ~ "medium", prop < 0.01 & counts > 1 ~ "low", .default = "no expansion")) %>% 
        group_by(expansion_level) %>% 
        dplyr::summarise(prop_sum = sum(prop))
      pdata_sum$expansion_level <- factor(pdata_sum$expansion_level, levels = c("high", "medium", "low", "no expansion"))
      pdata_sum <- pdata_sum %>% 
        arrange(expansion_level) %>%
        dplyr::mutate(
          start = 2 * pi - cumsum(lag(prop_sum, default = 0)) * 2 * pi, 
          end = 2 * pi - cumsum(prop_sum) * 2 * pi                       
        )
      
      pdata_sum <- pdata_sum %>% 
        dplyr::mutate(timepoint = time, tissue = ti)
      
      n <- n + 1
      pdata_list[[n]] <- pdata
      pdata_sum_list[[n]] <- pdata_sum
    }
  }
  write_tsv(pdata_list %>% Reduce("rbind", .), "source_fig_4A_inner.tsv", na = "")
  write_tsv(pdata_sum_list %>% Reduce("rbind", .), "source_fig_4A_outer.tsv", na = "")
}

# Fig 4B: TCR diversity and evenness
{
  tcr_paired <- readxl::read_excel("Lm-TCR4.57-clonotypes-paired in rna.xlsx") %>% 
    dplyr::mutate(tissue = str_extract(Tissue_time, "(.*)-", group = 1)) %>% 
    group_by(CDRa_CDRb, timepoint, tissue, barcode) %>% 
    distinct(CDRa_CDRb, timepoint, tissue, barcode, .keep_all = TRUE) %>% 
    ungroup()
  
  tcr_paired2 <- tcr_paired %>% 
    dplyr::mutate(tissue = str_extract(Tissue_time, "(.*)-", group = 1)) %>% 
    group_by(CDRa_CDRb, timepoint, tissue) %>% 
    dplyr::mutate(counts = n()) %>% 
    group_by(timepoint, tissue) %>% 
    dplyr::mutate(prop = counts/n()) %>% 
    arrange(desc(prop)) %>% 
    ungroup() %>% 
    dplyr::select(CDRa_CDRb, timepoint, tissue, counts, prop) %>% 
    distinct()
  
  tcr_paired_chao1 <- tcr_paired2 %>% 
    group_by(timepoint, tissue) %>% 
    dplyr::summarise(chao1 = n() + (sum(counts == 1)*(-1 + sum(counts == 1)))/(2*(1 + sum(counts == 2))),
                     shannon = -sum(prop * log(prop))) 
  tcr_paired_chao1$timepoint <- factor(tcr_paired_chao1$timepoint, levels = c("WT", "d3", "d5", "d8", "d14", "d30", "R5", "Y5"))
  
  lim <- 30000
  tcr_paired_chao1 <- tcr_paired_chao1 %>% dplyr::mutate(Chao1 = if_else(chao1 > lim, lim, chao1))
  
  
  
  tcr_paired_chao1 <- tcr_paired_chao1 %>% dplyr::mutate(cluster = case_when(shannon < 2.5 & Chao1 < 5000 ~ "1", 
                                                                             shannon > 2.5 & shannon < 5 & Chao1 < 5000 ~ "2",
                                                                             shannon < 7.5 & shannon > 5 & Chao1 < 15000 ~ "3",
                                                                             shannon < 5 & Chao1 > 10000 ~ "4", .default = "5"))
  tcr_paired_chao1$cluster <- factor(tcr_paired_chao1$cluster)
  
  library(ggnewscale)
  library(ggalt)
  tcr_paired_chao1 <- tcr_paired_chao1 %>% dplyr::mutate(Chao1 = if_else(chao1 > lim, lim, chao1))
  xrange <- max(tcr_paired_chao1$Chao1) - min(tcr_paired_chao1$Chao1)
  yrange <- max(tcr_paired_chao1$shannon) - min(tcr_paired_chao1$shannon)
  ggplot(tcr_paired_chao1, aes(x = Chao1, y = shannon)) +
    geom_point(aes(fill = timepoint, shape = tissue,  color = timepoint), size = 3.5, alpha = 0.8) +
    scale_fill_manual(values =   c(WT = "#A55628", 
                                   d3 = "#4CAD49",
                                   d5 = "#984e9d",
                                   d8 = "#f77f11",
                                   d14 = "#E31B1B",
                                   d30 = "#347eb4",
                                   R5 = "#f7f732",
                                   Y5 = "#F280B6")) +
    scale_color_manual(values =   c(WT = "#A55628", 
                                    d3 = "#4CAD49",
                                    d5 = "#984e9d",
                                    d8 = "#f77f11",
                                    d14 = "#E31B1B",
                                    d30 = "#347eb4",
                                    R5 = "#f7f732",
                                    Y5 = "#F280B6")) +
    scale_shape_manual(values = c(21, 22, 23, 24)) +
    theme_bw() +
    new_scale_fill() +
    new_scale_color() +
    # geom_mark_rect(aes(colour = cluster), linetype = 2, alpha = 0.1, size = 0.8, expand = 0.04) +
    scale_color_manual(values = c('#C9B2E6', '#5F12B9', "#F1BC31", "#8E2C00", "#86d9ac")) +
    # stat_ellipse(aes(color = cluster, fill = cluster), linetype = 2, alpha = 0.1, geom = "polygon", level = 0.8) +
    xlim(0, lim) +
    xlab("Chao1 (diversity)") +
    ylab("Shannon index (diversity + evenness)") +
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 10, color = "black"),
          axis.text.y = element_text(hjust = 0.5, vjust = 0.5, size = 10, color = "black"),
          axis.title = element_text(size = 13, color = "black")) +
    coord_fixed(ratio = xrange/yrange)
  # ggsave(filename = "/home/zeemeeuw/data/mmJoint/cleanTCR/chao1_shannon_new.pdf", width = 5, height = 4)
  write_tsv(tcr_paired_chao1, "source_fig_4B.tsv", na = "")
}

# Fig 4C: TRBV gene usage
{
  tcr_paired <- readxl::read_excel("Lm-TCR4.57-clonotypes-paired in rna.xlsx") %>% 
    dplyr::mutate(tissue = str_extract(Tissue_time, "(.*)-", group = 1)) %>% 
    group_by(CDRa_CDRb, timepoint, tissue, barcode) %>% 
    distinct(CDRa_CDRb, timepoint, tissue, barcode, .keep_all = TRUE) %>% 
    ungroup()
  pdata <- as.data.frame.array(table(tcr_paired$TRBV, tcr_paired$timepoint))

  pdata <- pdata[  c("TRBV12-1", "TRBV13-3", "TRBV3", "TRBV20", "TRBV31", "TRBV14", "TRBV13-2", "TRBV16", "TRBV8", "TRBV12-3", "TRBV22", "TRBV18", "TRBV25", "TRBV21", "TRBV30", "TRBV23", 
                     "TRBV24",   "TRBV17",  "TRBV4", "TRBV15", "TRBV5", "TRBV2", "TRBV26", "TRBV1",  "TRBV12-2", "TRBV19", "TRBV13-1", "TRBV29"), 
                 c("WT", "d3", "d5", "d8", "d14", "d30", "R5", "Y5")]
  pdata <- t(t(pdata)/colSums(pdata))
  pdata[pdata > 0.2] <- 0.2
  pheatmap::pheatmap(pdata,
                     show_colnames = TRUE,
                     show_rownames = TRUE,
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     treeheight_col = 0,
                     treeheight_row = 10,
                     border_color = NA,
                     clustering_method = "ward.D",
                     scale = "none",
                     color = choose_colorset("wolfgang_basic"),
                     # filename = "TRBV_gene_heatmap.pdf",
                     fontsize_row = 8,
                     width = 2.5, height = 4)
  
  pdata <- data.frame(pdata) %>% 
    tibble::rownames_to_column(var = "TRBV")
  write.table(pdata, "source_fig_4C.tsv", sep = "\t", quote = FALSE, na = "", row.names = FALSE, col.names = TRUE)
}

# Fig 4D: line plots of shared TCR
{
  data_longer <- readxl::read_excel("Lm-TCR4.57-clonotypes-paired in rna.xlsx")
  data <- data_longer %>% 
    dplyr::mutate(tissue = str_extract(Tissue_time, "(.*)-", group = 1)) %>% 
    group_by(CDRa_CDRb, timepoint, tissue, barcode) %>% 
    distinct(CDRa_CDRb, timepoint, tissue, barcode, .keep_all = TRUE) %>% 
    dplyr::filter(!str_detect(CDR3a, "_") & !str_detect(CDR3b, "_")) %>% 
    ungroup() %>% 
    dplyr::select(CDRa_CDRb, timepoint) %>% 
    group_by(CDRa_CDRb) %>% 
    dplyr::mutate(count = n()) %>% 
    group_by(CDRa_CDRb, timepoint) %>% 
    dplyr::mutate(count_timepoint = n()) %>%  
    ungroup() %>% 
    distinct() %>% 
    arrange(CDRa_CDRb)
  
  data_hmp <- data %>%
    dplyr::select(CDRa_CDRb, timepoint, count, count_timepoint) %>%
    group_by(timepoint) %>%
    dplyr::mutate(prop = count_timepoint/sum(count_timepoint)) %>%
    dplyr::select(CDRa_CDRb, timepoint, prop) %>%
    group_by(CDRa_CDRb) %>%
    dplyr::mutate(n = n()) %>%
    dplyr::filter(n > 1) %>%
    pivot_wider(names_from = timepoint, values_from = prop) %>%
    replace_na(list( d3 = 0,  d5 = 0, d8 = 0, d14 = 0, d30 = 0, R5 = 0, Y5 = 0))
  mat <- as.matrix(data_hmp[, c("d3", "d5", "d8", "d14", "d30", "R5", "Y5")])
  rownames(mat) <- data_hmp$CDRa_CDRb
  
  data <- data.frame(mat)
  data <- data.frame(data/rowSums(data)*100)
  data$CDR3a_CDR3b <- rownames(mat)
  data <- data %>% dplyr::mutate(identity = case_when(d5 >= d8 & d5 >= d14 & d5 >= d30 & d5 >= R5 & d5 >= Y5 ~ "d5",
                                                      d5 < d8 & d8 >= d14 & d8 >= d30 & d8 >= R5 & d8 >= Y5 ~ "d8",
                                                      d5 < d14 & d8 < d14 & d14 >= d30 & d14 >= R5 & d14 >= Y5 ~ "d14",
                                                      d5 < d30 & d8 < d30 & d14 < d30 & d30 >= R5 & d30 >= Y5 ~ "d30",
                                                      d5 < R5 & d8 < R5 & d14 < R5 & d30 < R5 & R5 >= Y5 ~ "R5",
                                                      .default = "Y5"))
  
  pdata <- data %>% pivot_longer(d3:Y5, names_to = "timepoint", values_to = "prop")
  pdata$timepoint <- factor(pdata$timepoint, levels = c("d3", "d5", "d8", "d14", "d30", "R5", "Y5"))

  time_color <- c(WT = "#A55628", 
    d3 = "#4CAD49",
    d5 = "#984e9d",
    d8 = "#f77f11",
    d14 = "#E31B1B",
    d30 = "#347eb4",
    R5 = "#f7f732",
    Y5 = "#F280B6")

  
  # d8
  ppdata1 <- pdata %>% dplyr::filter(identity == "d8")
  ppdata1_od <- ppdata1 %>% dplyr::filter(timepoint == "d8") %>%
    arrange(prop)
  ppdata1$CDR3a_CDR3b <- factor(ppdata1$CDR3a_CDR3b, levels = ppdata1_od$CDR3a_CDR3b)
  p1 <- ggplot(ppdata1) +
    geom_line(aes(x = timepoint, y = prop, group = CDR3a_CDR3b, color = CDR3a_CDR3b), show.legend = FALSE, linewidth = 1) +
    theme_classic() +
    ylab("Relative proportion (%)") +
    # scale_color_manual(values = colorRampPalette(c("#F7B4B4", "#ED8080", "#DF4A4A", "#BE2A2A", "#8C1616", "#590707"))(n_distinct(ppdata1$CDR3a_CDR3b))) +
    scale_color_manual(values = colorRampPalette(c("#FAC491", "#F77F11", "#633205"))(n_distinct(ppdata1$CDR3a_CDR3b))) +
    theme(axis.text.x = element_text(size = 10, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"))  +
    ylim(0, 100)  +
    facet_wrap(~identity)
  # d14
  ppdata2 <- pdata %>% dplyr::filter(identity == "d14")
  ppdata2_od <- ppdata2 %>% dplyr::filter(timepoint == "d14") %>%
    arrange(prop)
  ppdata2$CDR3a_CDR3b <- factor(ppdata2$CDR3a_CDR3b, levels = ppdata2_od$CDR3a_CDR3b)
  p2 <- ggplot(ppdata2) +
    geom_line(aes(x = timepoint, y = prop, group = CDR3a_CDR3b, color = CDR3a_CDR3b), show.legend = FALSE, linewidth = 1) +
    theme_classic() +
    ylab("Relative proportion (%)") +
    # scale_color_manual(values = colorRampPalette(c("#B6DDEB", "#86C2D8", "#56A6C3", "#3685A2", "#1C5F77", "#093B4D"))(n_distinct(ppdata2$CDR3a_CDR3b))) +
    scale_color_manual(values = colorRampPalette(c("#F29E9E", "#E31B1B", "#570707"))(n_distinct(ppdata2$CDR3a_CDR3b))) +
    theme(axis.text.x = element_text(size = 10, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"))  +
    ylim(0, 100)  +
    facet_wrap(~identity)
  # r5
  ppdata3 <- pdata %>% dplyr::filter(identity == "R5")
  ppdata3_od <- ppdata3 %>% dplyr::filter(timepoint == "R5") %>%
    arrange(prop)
  ppdata3$CDR3a_CDR3b <- factor(ppdata3$CDR3a_CDR3b, levels = ppdata3_od$CDR3a_CDR3b)
  p3 <- ggplot(ppdata3) +
    geom_line(aes(x = timepoint, y = prop, group = CDR3a_CDR3b, color = CDR3a_CDR3b), show.legend = FALSE, linewidth = 1) +
    theme_classic() +
    ylab("Relative proportion (%)") +
    scale_color_manual(values = colorRampPalette(c("#f7f732", "#A1A13D", "#45450A"))(n_distinct(ppdata3$CDR3a_CDR3b))) +
    theme(axis.text.x = element_text(size = 10, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"))  +
    ylim(0, 100)  +
    facet_wrap(~identity)
  # d5
  ppdata4 <- pdata %>% dplyr::filter(identity == "d5")
  ppdata4_od <- ppdata4 %>% dplyr::filter(timepoint == "d5") %>%
    arrange(prop)
  ppdata4$CDR3a_CDR3b <- factor(ppdata4$CDR3a_CDR3b, levels = ppdata4_od$CDR3a_CDR3b)
  p4 <- ggplot(ppdata4) +
    geom_line(aes(x = timepoint, y = prop, group = CDR3a_CDR3b, color = CDR3a_CDR3b), show.legend = FALSE, linewidth = 1) +
    theme_classic() +
    ylab("Relative proportion (%)") +
    # scale_color_manual(values = colorRampPalette(c("#B6EBD2", "#86D8B2", "#56C390", "#36A26F", "#1C774D", "#094D2D"))(n_distinct(ppdata4$CDR3a_CDR3b))) +
    scale_color_manual(values = colorRampPalette(c("#DEBCE0", "#984E9D", "#4A104D"))(n_distinct(ppdata4$CDR3a_CDR3b))) +
    theme(axis.text.x = element_text(size = 10, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"))  +
    ylim(0, 100)  +
    facet_wrap(~identity)
  # d30
  ppdata5 <- pdata %>% dplyr::filter(identity == "d30")
  ppdata5_od <- ppdata5 %>% dplyr::filter(timepoint == "d30") %>%
    arrange(prop)
  ppdata5$CDR3a_CDR3b <- factor(ppdata5$CDR3a_CDR3b, levels = ppdata5_od$CDR3a_CDR3b)
  p5 <- ggplot(ppdata5) +
    geom_line(aes(x = timepoint, y = prop, group = CDR3a_CDR3b, color = CDR3a_CDR3b), show.legend = FALSE, linewidth = 1) +
    theme_classic() +
    ylab("Relative proportion (%)") +
    # scale_color_manual(values = colorRampPalette(c("#B6BFEB", "#8694D8", "#5668C3", "#3648A2", "#1C2B77", "#09154D"))(n_distinct(ppdata5$CDR3a_CDR3b))) +
    scale_color_manual(values = colorRampPalette(c("#347eb4"))(n_distinct(ppdata5$CDR3a_CDR3b))) +
    theme(axis.text.x = element_text(size = 10, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"))  +
    ylim(0, 100)  +
    facet_wrap(~identity)
  
  ppdata6 <- pdata %>% dplyr::filter(identity == "Y5")
  ppdata6_od <- ppdata6 %>% dplyr::filter(timepoint == "Y5") %>%
    arrange(prop)
  ppdata6$CDR3a_CDR3b <- factor(ppdata6$CDR3a_CDR3b, levels = ppdata6_od$CDR3a_CDR3b)
  p6 <- ggplot(ppdata6) +
    geom_line(aes(x = timepoint, y = prop, group = CDR3a_CDR3b, color = CDR3a_CDR3b), show.legend = FALSE, linewidth = 1) +
    theme_classic() +
    ylab("Relative proportion (%)") +
    scale_color_manual(values = colorRampPalette(c("#F5C1D9", "#F280B6", "#66183B"))(n_distinct(ppdata6$CDR3a_CDR3b))) +
    theme(axis.text.x = element_text(size = 10, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"))  +
    ylim(0, 100)  +
    facet_wrap(~identity)
  
  wrap_plots(list(p4, p1, p2, p5, p3, p6), ncol = 1)
  # ggsave("/home/zeemeeuw/data/mmJoint/cleanTCR/mmjoint_TCR_module_lineplot.pdf", width = 3, height = 10)

  
  data <- bind_rows(ppdata1, ppdata2, ppdata3, ppdata4, ppdata5)
  write_tsv(data, "source_fig_4E.tsv", na = "")
}

# Fig 4F: Heatmap of clonotypes
{
  
  tcr_paired <- readxl::read_excel("mmJoint96_filtered_tcr_longer.xlsx") %>% 
    dplyr::mutate(TRBV = str_extract(TRBV, "(.*)\\*", group = 1),
                  TRBJ = str_extract(TRBJ, "(.*)\\*", group = 1))
  
  pdata <- as.data.frame.array(table(tcr_paired$CDR3a_CDR3b, tcr_paired$timepoint))
  pdata <- pdata[, 
                 c("D5", "D8", "D14", "D21", "R5")]
  mat <- t(t(pdata)/colSums(pdata))
  # pdata <- pdata[rowSums(pdata>0) > 1,]
  # pdata[pdata > 0.03] <- 0.03
  mat_scale <- sweep(mat - rowMeans(mat), 1, matrixStats::rowSds(mat), `/`)
  mat_scale <- data.frame(mat_scale)
  rownames(mat_scale) <- rownames(mat)
  mat2 <- binarysort_heatmap(mat_scale, scale = FALSE, cutOff = -0.44, clusterCols = FALSE, invert = FALSE)$mat
  mat2 <- data.frame(mat2)
  mat2[mat[rownames(mat2),] == 0] <- NA
  pheatmap::pheatmap(mat2,
                     show_colnames = TRUE,
                     show_rownames = FALSE,
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     treeheight_col = 0,
                     treeheight_row = 10,
                     border_color = NA,
                     clustering_method = "ward.D2",
                     scale = "none",
                     na_col = "white",
                     color = choose_colorset("fireworks"),
                     breaks = seq(-1.5, 1.5, length.out = 256),
                     legend_breaks = c(-1.5, 0, 1.5),
                     # filename = "sc96_shared_clone_scale_allTCR_heatmap.pdf",
                     fontsize_row = 8,
                     width = 1.8, height = 3.4)
}


# Fig 4G: lineplot, 96 well
{
  data_longer <- read_excel(path = "mmJoint96_filtered_tcr_longer.xlsx") 
  # extract shared clone
  data <- data_longer %>% 
    filter(!is.na(TRA_CDR3_aa) & !is.na(TRB_CDR3_aa)) %>% 
    filter(!str_detect(TRA_CDR3_aa, "_") & !str_detect(TRB_CDR3_aa, "_")) %>% 
    dplyr::select(CDR3a_CDR3b, timepoint) %>% 
    group_by(CDR3a_CDR3b) %>% 
    dplyr::filter(n_distinct(timepoint) > 1) %>% 
    dplyr::mutate(count = n()) %>% 
    group_by(CDR3a_CDR3b, timepoint) %>% 
    dplyr::mutate(count_timepoint = n()) %>%  
    ungroup() %>% 
    distinct() %>% 
    arrange(CDR3a_CDR3b)
  
  data_hmp <- data %>%
    group_by(timepoint) %>%
    dplyr::mutate(prop = count_timepoint/sum(count_timepoint)) %>%
    dplyr::select(CDR3a_CDR3b, timepoint, prop) %>% 
    pivot_wider(names_from = timepoint, values_from = prop) %>%
    replace_na(list(D5 = 0, D8 = 0, D14 = 0, D21 = 0, R5 = 0))
  # line plot, clone dominance at each timepoint
  data_dom <- data_hmp %>% 
    mutate(dominance = names(dplyr::select(., D5, D8, D14, D21, R5))[max.col(dplyr::select(., D5, D8, D14, D21, R5), "first")]) %>% 
    pivot_longer(cols = c(D5, D8, D14, D21, R5), names_to = "timepoint", values_to = "prop") %>% 
    dplyr::mutate(freq = 100 * prop) %>% 
    group_by(CDR3a_CDR3b) %>% 
    dplyr::mutate(norm_freq = prop/sum(prop)*100) %>% 
    dplyr::mutate_at(vars(timepoint), ~factor(.x, levels = c("D5", "D8", "D14", "D21", "R5"))) %>% 
    dplyr::mutate_at(vars(dominance), ~factor(.x, levels = c("D5", "D8", "D14", "D21", "R5")))
  
  tcr_d5 <- data_dom %>% dplyr::filter(dominance == "D5") %>% arrange(desc(prop)) %>% pull(CDR3a_CDR3b) %>% unique()
  tcr_d8 <- data_dom %>% dplyr::filter(dominance == "D8") %>% arrange(desc(prop)) %>% pull(CDR3a_CDR3b) %>% unique()
  tcr_d14 <- data_dom %>% dplyr::filter(dominance == "D14") %>% arrange(desc(prop)) %>% pull(CDR3a_CDR3b) %>% unique()
  tcr_d21 <- data_dom %>% dplyr::filter(dominance == "D21") %>% arrange(desc(prop)) %>% pull(CDR3a_CDR3b) %>% unique()
  tcr_R5 <- data_dom %>% dplyr::filter(dominance == "R5") %>% arrange(desc(prop)) %>% pull(CDR3a_CDR3b) %>% unique()
  
  color_d5 <- rev(colorRampPalette(c("#DEBCE0", "#984E9D", "#4A104D"))(n_distinct(tcr_d5)))
  names(color_d5) <- tcr_d5
  color_d8 <- rev(colorRampPalette(c("#FAC491", "#F77F11", "#703A07"))(n_distinct(tcr_d8)))
  names(color_d8) <- tcr_d8
  color_d14 <- rev(colorRampPalette(c("#F29E9E", "#E31B1B", "#570707"))(n_distinct(tcr_d14)))
  names(color_d14) <- tcr_d14
  color_d21 <- rev(colorRampPalette(c("#B1D7F5", "#347EB4", "#15354D"))(n_distinct(tcr_d21)))
  names(color_d21) <- tcr_d21
  color_R5 <- rev(colorRampPalette(c("#f7f732", "#A1A13D", "#45450A"))(n_distinct(tcr_R5)))
  names(color_R5) <- tcr_R5
  colors_used <- c(color_d5, color_d8, color_d14, color_d21, color_R5)
  
  # relative proportion
  ggplot(data_dom) +
    geom_line(aes(x = timepoint, y = norm_freq, group = CDR3a_CDR3b, color = CDR3a_CDR3b), show.legend = FALSE, linewidth = 0.5) +
    theme_classic() +
    scale_color_manual(values = colors_used) +
    theme(axis.text.x = element_text(size = 10, color = "black"),
          axis.text.y = element_text(size = 10, color = "black")) +
    facet_wrap(~dominance, ncol = 1, scales = "free") +
    ylab("Relative proportion %") +
    ylim(0, 100)
}


# Fig 4H-I: TCR classification
{
  # 4I
  data_longer <- read_excel(path = "mmJoint96_filtered_tcr_longer.xlsx") 
  # extract shared clone
  data <- data_longer %>% 
    filter(!is.na(TRA_CDR3_aa) & !is.na(TRB_CDR3_aa)) %>% 
    filter(!str_detect(TRA_CDR3_aa, "_") & !str_detect(TRB_CDR3_aa, "_")) %>% 
    dplyr::select(CDR3a_CDR3b, timepoint) %>% 
    group_by(CDR3a_CDR3b) %>% 
    dplyr::filter(n_distinct(timepoint) > 1) %>% 
    dplyr::mutate(count = n()) %>% 
    group_by(CDR3a_CDR3b, timepoint) %>% 
    dplyr::mutate(count_timepoint = n()) %>%  
    ungroup() %>% 
    distinct() %>% 
    arrange(CDR3a_CDR3b)
  
  data_hmp <- data %>%
    group_by(timepoint) %>%
    dplyr::mutate(prop = count_timepoint/sum(count_timepoint)) %>%
    dplyr::select(CDR3a_CDR3b, timepoint, prop) %>% 
    pivot_wider(names_from = timepoint, values_from = prop) %>%
    replace_na(list(D5 = 0, D8 = 0, D14 = 0, D21 = 0, R5 = 0))
  
  data_grid <- data_hmp %>% 
    rowwise() %>% 
    dplyr::mutate(d8 = D14 <= D8 & D8 >= D21, 
                  d14 = D8 <= D14 & D14 >= D21,
                  d21 = D8 <= D21 & D21 >= D14) %>% 
    dplyr::mutate(increase = R5 > max(c(D8, D14, D21)),
                  decrease = R5 < max(c(D8, D14, D21)) & R5 >= min(c(D8, D14, D21)),
                  decay = R5 < min(c(D8, D14, D21))) %>% 
    dplyr::mutate(group = case_when(d8 & increase ~ "d8_increase", d8 & decrease ~ "d8_decrease", d8 & decay ~ "d8_decay",
                                    d14 & increase ~ "d14_increase", d14 & decrease ~ "d14_decrease", d14 & decay ~ "d14_decay",
                                    d21 & increase ~ "d21_increase", d21 & decrease ~ "d21_decrease", d21 & decay ~ "d21_decay")) %>% 
    dplyr::mutate_at(vars(group), ~factor(.x, levels = c("d8_increase", "d8_decrease", "d8_decay", "d14_increase", "d14_decrease", "d14_decay", "d21_increase", "d21_decrease", "d21_decay"))) %>% 
    dplyr::select(D8, D14, D21, R5, group, CDR3a_CDR3b) %>% 
    pivot_longer(cols = c(D8, D14, D21, R5), names_to = "timepoint", values_to = "prop") %>% 
    dplyr::mutate(freq = 100 * prop) %>% 
    group_by(CDR3a_CDR3b) %>% 
    dplyr::mutate(norm_freq = prop/sum(prop)*100) %>% 
    dplyr::mutate_at(vars(timepoint), ~factor(.x, levels = c("D8", "D14", "D21", "R5")))
  
  tcr_d8 <- data_grid %>% dplyr::filter(group %in% c("d8_increase", "d8_decrease", "d8_decay")) %>% arrange(desc(prop)) %>% pull(CDR3a_CDR3b) %>% unique()
  tcr_d14 <- data_grid %>% dplyr::filter(group %in% c("d14_increase", "d14_decrease", "d14_decay")) %>% arrange(desc(prop)) %>% pull(CDR3a_CDR3b) %>% unique()
  tcr_d21 <- data_grid %>% dplyr::filter(group %in% c("d21_increase", "d21_decrease", "d21_decay")) %>% arrange(desc(prop)) %>% pull(CDR3a_CDR3b) %>% unique()
  
  color_d8 <- rev(colorRampPalette(c("#FAC491", "#F77F11", "#703A07"))(n_distinct(tcr_d8)))
  names(color_d8) <- tcr_d8
  color_d14 <- rev(colorRampPalette(c("#F29E9E", "#E31B1B", "#570707"))(n_distinct(tcr_d14)))
  names(color_d14) <- tcr_d14
  color_d21 <- rev(colorRampPalette(c("#B1D7F5", "#347EB4", "#15354D"))(n_distinct(tcr_d21)))
  names(color_d21) <- tcr_d21
  colors_used <- c(color_d8, color_d14, color_d21)
  
  # relative proportion
  ggplot(data_grid) +
    geom_line(aes(x = timepoint, y = norm_freq, group = CDR3a_CDR3b, color = CDR3a_CDR3b), show.legend = FALSE, linewidth = 0.5) +
    theme_classic() +
    scale_color_manual(values = colors_used) +
    theme(axis.text.x = element_text(size = 10, color = "black"),
          axis.text.y = element_text(size = 10, color = "black")) +
    facet_wrap(~group, ncol = 3, scales = "free") +
    ylab("Relative proportion %") +
    ylim(0, 100)
  
  # 4H
  data_summary <- data_grid %>% dplyr::mutate(pattern = str_extract(group, "_(.*)", group = 1)) %>% 
    dplyr::mutate_at(vars(pattern), ~factor(.x, levels = c("increase", "decrease", "decay"))) %>% 
    group_by(timepoint, pattern) %>% 
    dplyr::summarise(prop = sum(prop), freq = sum(freq))
  ggplot(data_summary) +
    geom_line(aes(x = timepoint, y = freq, group = pattern, color = pattern), show.legend = FALSE, linewidth = 1) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 10, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"),
          strip.text = element_text(size = 10),
          strip.background = element_rect(fill = c("#EBE5E1"))) +
    ylim(0, 70) +
    scale_color_manual(values = c("#E31B1B", "#56C35D", "#347EB4")) +
    ylab("Aggregated proportion in shared clonotypes %")
}










