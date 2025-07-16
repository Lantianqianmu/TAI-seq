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
setwd("/home/zeemeeuw/YangLab/JoINT-seq/Manuscript/fig_data/fig_S10_data")

RDSdir <- "/home/zeemeeuw/data/mmJoint_WTx/RNA/rds/" # output RDS file dir

load(paste0(RDSdir, "lm_merged_overlapped_harmony_GRN_list_by_tissue.RData"))

grn <- grn_list$grn_ln
# ln DORC driver plot: Klrg1
d <- grn %>% dplyr::filter(DORC == "Klrg1") %>% 
  dplyr::mutate(isSig = if_else(abs(Score) >= 0.8 & abs(Corr) >= 0.7, "Yes", "No"))
d$Label <- d$Motif
d$Label[d$isSig %in% "No"] <- ""

ggplot(data = d, aes(x = Corr, y = Enrichment.log10P, color = isSig)) + 
  geom_hline(yintercept = -log10(0.01), color = "gray60", linetype = "dashed") + 
  geom_vline(xintercept = c(0), color = "gray60", linetype = "dashed") + 
  geom_point(aes(size = isSig)) + 
  theme_classic() + 
  scale_color_manual(values = c("gray66", "firebrick3")) + 
  scale_size_manual(values = c(1, 3)) +
  labs(y = "Enrichment log10 P", x = "DORC-TF Correlation", title = "Klrg1") + 
  ylim(-10, 15) + 
  xlim(-1, 1) + 
  theme(legend.position = "none", 
        axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 10, color = "black"),
        axis.text.y = element_text(hjust = 0.5, vjust = 0.5, size = 10, color = "black"),
        plot.title = element_text(hjust = 0.5, face = "italic"), 
        panel.background = element_rect(fill = NA)) + 
  coord_fixed(2/25) +
  ggrepel::geom_text_repel(aes(label = Label), size = 5, color = "black", vjust = 0.5, hjust = 0.5, max.overlaps = 50)

# spleen DORC driver plot: Bcl2
grn <- grn_list$grn_sp
d <- grn %>% dplyr::filter(DORC == "Bcl2") %>% 
  dplyr::mutate(isSig = if_else(abs(Score) >= 0.8 & abs(Corr) >= 0.7, "Yes", "No"))
d$Label <- d$Motif
d$Label[d$isSig %in% "No"] <- ""

ggplot(data = d, aes(x = Corr, y = Enrichment.log10P, color = isSig)) + 
  geom_hline(yintercept = -log10(0.01), color = "gray60", linetype = "dashed") + 
  geom_vline(xintercept = c(-0.7, 0.7), color = "gray60", linetype = "dashed") + 
  geom_point(aes(size = isSig)) + 
  theme_classic() + 
  scale_color_manual(values = c("gray66", "firebrick3")) + 
  scale_size_manual(values = c(1, 3)) +
  labs(y = "Enrichment log10 P", x = "DORC-TF Correlation", title = "Bcl2") + 
  ylim(-10, 15) + 
  xlim(-1, 1) + 
  theme(legend.position = "none", 
        axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 10, color = "black"),
        axis.text.y = element_text(hjust = 0.5, vjust = 0.5, size = 10, color = "black"),
        plot.title = element_text(hjust = 0.5, face = "italic"), 
        panel.background = element_rect(fill = NA)) + 
  coord_fixed(2/25) +
  ggrepel::geom_text_repel(aes(label = Label), size = 5, color = "black", vjust = 0.5, hjust = 0.5, max.overlaps = 50)
ggsave("Spleen_Bcl2_driver.pdf", width = 4, height = 4)




# DORC driver plot: Havcr2, all tissue GRN
load(paste0(RDSdir, "lm_merged_overlapped_harmony_GRN_2.RData"))

d <- grn %>% dplyr::filter(DORC == "Havcr2") %>% 
  # dplyr::mutate(isSig = if_else(abs(Score) >= 0.8 & abs(Corr) >= 0.7, "Yes", "No"))
  dplyr::mutate(isSig = if_else(Enrichment.log10P >= 2, "Yes", "No"))
d$Label <- d$Motif
d$Label[d$isSig %in% "No"] <- ""

ggplot(data = d, aes(x = Corr, y = Enrichment.log10P, color = isSig)) + 
  geom_hline(yintercept = -log10(0.01), color = "gray60", linetype = "dashed") + 
  geom_vline(xintercept = c(-0.7, 0.7), color = "gray60", linetype = "dashed") + 
  geom_point(aes(size = isSig)) + 
  theme_classic() + 
  scale_color_manual(values = c("gray66", "firebrick3")) + 
  scale_size_manual(values = c(1, 3)) +
  labs(y = "Enrichment log10 P", x = "DORC-TF Correlation", title = "Havcr2") + 
  ylim(-10, 15) + 
  xlim(-1, 1) + 
  theme(legend.position = "none", 
        axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 10, color = "black"),
        axis.text.y = element_text(hjust = 0.5, vjust = 0.5, size = 10, color = "black"),
        plot.title = element_text(hjust = 0.5, face = "italic"), 
        panel.background = element_rect(fill = NA)) + 
  coord_fixed(2/25) +
  ggrepel::geom_text_repel(aes(label = Label), size = 5, color = "black", vjust = 0.5, hjust = 0.5, max.overlaps = 1000)
