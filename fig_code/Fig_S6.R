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
library(circlize)
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


seu <- readRDS(paste0(RDSdir, "lm_merged_overlapped_harmony.rds"))
proj <- loadArchRProject(ATACdir)

setwd("/home/zeemeeuw/YangLab/JoINT-seq/Manuscript/fig_data/fig_S6_data")


# S6B
{
  
  CellDimPlot(seu, reduction = "UMAP_LSI_ATAC",  theme = "theme_blank", legend.position = "right", group_by = "cell_module", raster = TRUE, pt_size = 0.01, palcolor = mmjoint_colors)

  CellDimPlot(seu, reduction = "UMAP_LSI_RNA",  theme = "theme_blank", legend.position = "right", group_by = "cell_module", raster = TRUE, pt_size = 0.01, palcolor = mmjoint_colors)

}


# S6E
{
  
  n_label <- 20
  plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = FALSE, n = n_label)
  plotVarDev <- data.frame(plotVarDev)
  ggplot(plotVarDev, aes(x = rank, y = combinedVars)) +
    theme_bw() +
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 10, color = "black"),
          axis.text.y = element_text(hjust = 0.5, vjust = 0.5, size = 10, color = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12, color = "black"),
          panel.grid = element_blank()) +
    geom_point(aes(color = combinedVars), size = 1) +
    ggrepel::geom_label_repel(data = plotVarDev[rev(seq_len(n_label)), 
    ], aes(x = rank, y = combinedVars, label = name), 
    size = 1.5, color = "black", nudge_x = 2, max.overlaps = 2000) +
    scale_color_gradientn(colors = choose_colorset("comet"))
}














