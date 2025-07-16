library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(ggnewscale)
library(ggforce)
assignIdents <- function(x, min_frags_human, min_frags_mouse = min_frags_human){
  human <- as.integer(x[2])
  mouse <- as.integer(x[3])
  if(human < min_frags_human && mouse < min_frags_mouse)(return("low yield"))
  if(human/mouse > 4){
    return("human")
  }else if(mouse/human > 4){
    return("mouse")
  }else{
    return("mix")
  }
}
read_peaks_merge <- function(files, width = 250, non_overlapping = TRUE){
  cn <- c("chr", "start", "end", "name", "score", "strand", "fc", "pval", "qval", "summit")
  beds <- lapply(files, function(x) read.delim(file = x, header = FALSE, sep = "\t", stringsAsFactors = FALSE, col.names = cn))
  beds <- lapply(beds, function(x) {x[, "summit"] <- x[, "start"] + x[, "summit"]; return(x)})
  bed <- do.call(rbind, beds)
  bed <- as(bed, "DataFrame")
  
  bed <- makeGRangesFromDataFrame(bed[, c("chr", "summit", "score", "qval", "name")], start.field = "summit", end.field = "summit", keep.extra.columns = TRUE)
  bed <- resize(bed, width = width, fix = "center")
  
  if (non_overlapping) {
    bed <- sortSeqlevels(bed)
    bed <- sort(bed)
    keep_peaks <- seq_along(bed)
    while (!(isDisjoint(bed[keep_peaks]))) {
      chr_names <- as.character(seqnames(bed[keep_peaks]))
      starts <- start(bed[keep_peaks])
      ends <- end(bed[keep_peaks])
      overlap_first <- intersect(which(chr_names[seq_len(length(keep_peaks) - 1)] == chr_names[seq_len(length(keep_peaks) - 1) + 1]),
                                 which(ends[seq_len(length(keep_peaks) - 1)] >= starts[seq_len(length(keep_peaks) - 1) + 1]))
      overlap_second <- overlap_first + 1
      overlap_comparison <- bed[keep_peaks[overlap_second]]$qval > bed[keep_peaks[overlap_first]]$qval
      discard <- keep_peaks[c(overlap_second[!overlap_comparison],
                              overlap_first[overlap_comparison])]
      keep_peaks <- keep_peaks[!(keep_peaks %in% discard)]
    }
    bed <- bed[keep_peaks]
  }
  return(bed)
}
get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x, y, n = 100)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
count_matrix_from_fragments <- function(file, whitelist = NULL, peaks = NULL, min_fragments = NULL){
  # read fragment file
  fragments <- readr::read_delim(file, col_names = FALSE)
  colnames(fragments) <- c("chr", "start", "end", "cell", "counts")
  
  if(!is.null(min_fragments)){
    message("Filtering cells according to min_fragments...")
    fragments <- fragments %>% group_by(cell) %>% mutate(n = n())
    fragments <- fragments[fragments$n >= min_fragments,]
    fragments <- fragments[c("chr", "start", "end", "cell", "counts")]
  }
  
  if(is.null(whitelist)){
    message("Skip whitelist filtering.")
    fragments <- fragments %>% group_by(cell)
  }else{
    message("Filtering cells according to barcode whitelist...")
    fragments <- fragments %>% mutate(valid = cell %in% whitelist)
    fragments <- fragments[fragments$valid,]
    fragments <- fragments[c("chr", "start", "end", "cell", "counts")]
    fragments <- fragments %>% group_by(cell)
  }
  
  message("Splitting cell data...")
  fragment_grouped <- fragments %>% group_split()
  
  message("Creating cell*peak matrix...")
  results <- lapply(fragment_grouped, function(x, peaks) {
    fragment_ranges <- GRanges(x$chr, ranges = IRanges(x$start, x$end))
    return(list(counts = countOverlaps(peaks, fragment_ranges, type = "any", ignore.strand = TRUE), fragment_counts = length(fragment_ranges)))
  }, peaks = peaks)
  
  mat <- Matrix::Matrix(vapply(results, function(x) x[["counts"]], rep(0, length(peaks))))
  colnames(mat) <- (fragments %>% group_keys())$cell
  fragment_counts <- vapply(results, function(x) x[["fragment_counts"]], 0)
  
  colData <- data.frame(fragment_counts = fragment_counts)
  
  out <- SummarizedExperiment(assays = list(counts = mat),
                              rowRanges = peaks, colData = colData)
  return(out)
}
source("/home/zeemeeuw/YangLab/JoINT-seq/QC/code/utils.R")
setwd("/home/zeemeeuw/YangLab/JoINT-seq/Manuscript/fig_data/fig_S1_data")

# Fig S1C: species mixing 37, ATAC
{
  hg38 <- read.delim("G3MERGE-ata_hg38_1_fragments_count.txt",
                     header = FALSE, stringsAsFactors = FALSE, col.names = c("barcode", "count.hg38"))
  mm10 <- read.delim("G3MERGE-ata_mm10_1_fragments_count.txt",
                     header = FALSE, stringsAsFactors = FALSE, col.names = c("barcode", "count.mm10"))
  
  minFrags_human <- 3000
  minFrags_mouse <- 3000
  data <- merge(hg38, mm10, by = 'barcode', all = TRUE)
  data[is.na(data)] <- 0
  data$ident <- apply(data, 1, assignIdents, min_frags_human = minFrags_human, min_frags_mouse = minFrags_mouse)
  data$ident <- factor(data$ident, levels = c("low yield", "human", "mouse", "mix"))
  median_hg38 <- median(data$count.hg38[data$ident == "human"])
  median_mm10 <- median(data$count.mm10[data$ident == "mouse"])
  drate <- round(100*sum(data$ident == "mix")/sum(data$ident != "low yield"), 3)
  data$Type <- data$ident
  data <- data %>% arrange(factor(data$ident, levels = c("mix", "low yield", "human", "mouse")))
  data$x <- data$count.hg38/1000
  data$y <- data$count.mm10/1000
  
  x_limit <- 60
  y_limit <- 60
  data_sum <- data.frame(table(data$Type))[2:4,]
  colnames(data_sum) <- c("Type", "Counts")
  data_sum$Type <- factor(data_sum$Type, levels = c("mouse", "human", "mix"))
  data_sum$x <- x_limit*0.7
  data_sum$y <- seq(from = y_limit*0.7-y_limit/10, to = y_limit*0.7+y_limit/10, by = y_limit/10)
  
  data$Type <- factor(data$Type, levels = c("low yield", "human", "mouse", "mix"))
  ggplot(data, aes(x = x, y = y)) +
    geom_point(size = 0.2, aes(colour = Type)) +
    scale_color_manual(values = c("#AAAAAA", "#4169E1", "#DC143C", "#c9c9c9")) +
    theme_classic() +
    # labs(title = paste0("ATAC | min ", minFrags_human, " ", minFrags_mouse, " | ", sum(data$ident != "low yield"), " cells | ", drate, "%")) +
    labs(title = "ATAC") +
    scale_x_continuous(limits = c(0, x_limit)) +
    scale_y_continuous(limits = c(0, y_limit)) +
    geom_text(data = data_sum, aes(label = paste0(Type, " (", Counts, ")"), x = x, y = y), size = 4) +
    xlab(bquote(paste('hg38 unique fragments (×', 10^3, ')', ' (', .(median_hg38), ')'))) +
    ylab(bquote(paste('mm10 unique fragments (×', 10^3, ')', ' (', .(median_mm10), ')'))) +
    xlab(bquote(paste('hg38 unique fragments (×', 10^3, ')'))) +
    ylab(bquote(paste('mm10 unique fragments (×', 10^3, ')'))) +
    theme(plot.title = element_text(hjust = 0.5)) +
    coord_fixed(ratio = x_limit/y_limit) +
    guides(color = guide_legend(override.aes = list(size = 2)))

}


# Fig S1D: species mixing 37, RNA
{

  hg38 <- read.delim("G3MERGE-rna_hg38_1_umis_count.txt",
                     header = FALSE, stringsAsFactors = FALSE, col.names = c("barcode", "count.hg38"))
  mm10 <- read.delim("G3MERGE-rna_mm10_1_umis_count.txt",
                     header = FALSE, stringsAsFactors = FALSE, col.names = c("barcode", "count.mm10"))
  minFrags_human <- 1000
  minFrags_mouse <- 3000
  data <- merge(hg38, mm10, by = 'barcode', all = TRUE)
  data[is.na(data)] <- 0
  data$ident <- apply(data, 1, assignIdents, min_frags_human = minFrags_human, min_frags_mouse = minFrags_mouse)
  data$ident <- factor(data$ident, levels = c("low yield", "human", "mouse", "mix"))
  median_hg38 <- median(data$count.hg38[data$ident == "human"])
  median_mm10 <- median(data$count.mm10[data$ident == "mouse"])
  drate <- round(100*sum(data$ident == "mix")/sum(data$ident != "low yield"), 3)
  data$Type <- data$ident
  data <- data %>% arrange(factor(data$ident, levels = c("mix", "low yield", "human", "mouse")))
  data$x <- data$count.hg38/1000
  data$y <- data$count.mm10/1000
  x_limit <- 30
  y_limit <- 50
  data_sum <- data.frame(table(data$Type))[2:4,]
  colnames(data_sum) <- c("Type", "Counts")
  data_sum$Type <- factor(data_sum$Type, levels = c("mouse", "human", "mix"))
  data_sum$x <- x_limit*0.7
  data_sum$y <- seq(from = y_limit*0.7-y_limit/10, to = y_limit*0.7+y_limit/10, by = y_limit/10)
  
  data$Type <- factor(data$Type, levels = c("low yield", "human", "mouse", "mix"))
  ggplot(data, aes(x = x, y = y)) +
    geom_point(size = 0.1, aes(colour = Type)) +
    scale_color_manual(values = c("#AAAAAA", "#4169E1", "#DC143C", "#c9c9c9")) +
    scale_x_continuous(limits = c(0, x_limit)) +
    scale_y_continuous(limits = c(0, y_limit)) +
    theme_classic() +
    # labs(title = paste0("RNA | min ", minFrags_human, " ", minFrags_mouse, " | ", sum(data$ident != "low yield"), " cells | ", drate, "%")) +
    labs(title = "RNA") +
    geom_text(data = data_sum, aes(label = paste0(Type, " (", Counts, ")"), x = x, y = y), size = 3) +
    xlab(bquote(paste('hg38 unique umis (', 10^3, ')', ' (', .(median_hg38), ')'))) +
    ylab(bquote(paste('mm10 unique umis (', 10^3, ')', ' (', .(median_mm10), ')'))) +
    xlab(bquote(paste('hg38 unique fragments (×', 10^3, ')'))) +
    ylab(bquote(paste('mm10 unique fragments (×', 10^3, ')'))) +
    theme(plot.title = element_text(hjust = 0.5)) +
    coord_fixed(ratio = x_limit/y_limit) +
    guides(color = guide_legend(override.aes = list(size = 2)))

}





# Fig S1E, NIH/3T3
{
  data <- read_tsv("source_fig_S1E_NIH3T3.tsv")
  
  model <- lm(y1~x1, data = data)
  model_summary <- summary(model)
  R_square <- model_summary$r.squared
  R_value <- round(sqrt(R_square), digits = 2)
  x_max <- 11
  y_max <- 11
  data_anno <- data.frame(x = 4+(x_max-4)/8, y = y_max-(y_max-4)/20, value = R_value)
  
  ggplot(data = data, aes(x = x, y = y)) +
    geom_point(size = 0.1, aes(color = density)) +
    geom_text(data = data_anno, aes(label = paste0("R = ", value), x = x, y = y), size = 5) +
    theme_classic() +
    scale_color_gradientn(colors = BuenColors::jdb_palette("solar_extra"), breaks = c(min(data$density), max(data$density)), labels = c("low", "high")) +
    guides(color = guide_colorbar(frame.colour = "black", frame.linewidth = 1, frame.linetype = 1, barwidth = 0.7, barheight = 2, ticks = FALSE)) +
    theme(axis.text.x = element_text(size = 13, color = "black"),
          axis.text.y = element_text(size = 13, color = "black"),
          axis.title.y = element_text(size = 15),
          axis.title.x = element_text(size = 15)) +
    scale_x_continuous(breaks = c(4, 5, 6, 7, 8, 9, 10, 11), limits = c(4, x_max)) +
    scale_y_continuous(breaks = c(4, 5, 6, 7, 8, 9, 10, 11), limits = c(4, y_max)) +
    xlab(bquote(Replicate~1~log[2](fragments))) +
    ylab(bquote(Replicate~2~log[2](fragments))) +
    coord_fixed(ratio = 1)
  
}


# Fig S1E, GM12878
{
  data <- read_tsv("source_fig_S1E_GM12878.tsv")
  model <- lm(y1~x1, data = data)
  model_summary <- summary(model)
  R_square <- model_summary$r.squared
  R_value <- round(sqrt(R_square), digits = 2)
  x_max <- 10
  y_max <- 10
  data_anno <- data.frame(x = 4+(x_max-4)/8, y = y_max-(y_max-4)/20, value = R_value)
  
  ggplot(data = data, aes(x = x, y = y)) +
    geom_point(size = 0.1, aes(color = density)) +
    geom_text(data = data_anno, aes(label = paste0("R = ", value), x = x, y = y), size = 5) +
    theme_classic() +
    scale_color_gradientn(colors = BuenColors::jdb_palette("solar_extra"), breaks = c(min(data$density), max(data$density)), labels = c("low", "high")) +
    guides(color = guide_colorbar(frame.colour = "black", frame.linewidth = 1, frame.linetype = 1, barwidth = 0.7, barheight = 2, ticks = FALSE)) +
    theme(axis.text.x = element_text(size = 13, color = "black"),
          axis.text.y = element_text(size = 13, color = "black"),
          axis.title.y = element_text(size = 15),
          axis.title.x = element_text(size = 15)) +
    scale_x_continuous(breaks = c(4, 5, 6, 7, 8, 9, 10, 11), limits = c(4, x_max)) +
    scale_y_continuous(breaks = c(4, 5, 6, 7, 8, 9, 10, 11), limits = c(4, y_max)) +
    xlab(bquote(Replicate~1~log[2](fragments))) +
    ylab(bquote(Replicate~2~log[2](fragments))) +
    coord_fixed(ratio = 1)
  
}

# Fig S1G, NIH3T3
{

  data <- read_tsv("source_fig_S1F_NIH3T3.tsv")
  model <- lm(y1~x1, data = data)
  model_summary <- summary(model)
  R_square <- model_summary$r.squared
  R_value <- round(sqrt(R_square), digits = 2)
  x_max <- 16
  y_max <- 16
  data_anno <- data.frame(x = 4+(x_max-4)/8, y = y_max-(y_max-4)/20, value = R_value)
  
  ggplot(data = data, aes(x = x, y = y)) +
    geom_point(size = 0.1, aes(color = density)) +
    geom_text(data = data_anno, aes(label = paste0("R = ", value), x = x, y = y), size = 5) +
    theme_classic() +
    scale_color_gradientn(colors = viridis::viridis(100), breaks = c(min(data$density), max(data$density)), labels = c("low", "high")) +
    guides(color = guide_colorbar(frame.colour = "black", frame.linewidth = 1, frame.linetype = 1, barwidth = 0.7, barheight = 2, ticks = FALSE)) +
    scale_x_continuous(breaks = c(4, 8, 12, 16), limits = c(4, x_max)) +
    scale_y_continuous(breaks = c(4, 8, 12, 16), limits = c(4, y_max)) +
    xlab(bquote(Replicate~1~log[2](umis))) +
    ylab(bquote(Replicate~2~log[2](umis))) +
    theme(axis.text.x = element_text(size = 13, color = "black"),
          axis.text.y = element_text(size = 13, color = "black"),
          axis.title.y = element_text(size = 15),
          axis.title.x = element_text(size = 15)) +
    coord_fixed(ratio = 14/14)

}

# Fig S1G, GM12878
{
  data <- read_tsv("source_fig_S1F_GM12878.tsv")
  model <- lm(y1~x1, data = data)
  model_summary <- summary(model)
  R_square <- model_summary$r.squared
  R_value <- round(sqrt(R_square), digits = 2)
  x_max <- 14
  y_max <- 14
  data_anno <- data.frame(x = 4+(x_max-4)/8, y = y_max-(y_max-4)/20, value = R_value)
  
  ggplot(data = data, aes(x = x, y = y)) +
    geom_point(size = 0.1, aes(color = density)) +
    geom_text(data = data_anno, aes(label = paste0("R = ", value), x = x, y = y), size = 5) +
    theme_classic() +
    scale_color_gradientn(colors = viridis::viridis(100), breaks = c(min(data$density), max(data$density)), labels = c("low", "high")) +
    guides(color = guide_colorbar(frame.colour = "black", frame.linewidth = 1, frame.linetype = 1, barwidth = 0.7, barheight = 2, ticks = FALSE)) +
    scale_x_continuous(breaks = c(4, 8, 12, 16), limits = c(4, x_max)) +
    scale_y_continuous(breaks = c(4, 8, 12, 16), limits = c(4, y_max)) +
    xlab(bquote(Replicate~1~log[2](umis))) +
    ylab(bquote(Replicate~2~log[2](umis))) +
    theme(axis.text.x = element_text(size = 13, color = "black"),
          axis.text.y = element_text(size = 13, color = "black"),
          axis.title.y = element_text(size = 15),
          axis.title.x = element_text(size = 15)) +
    coord_fixed(ratio = 14/14)
}


# Fig S1G
{
  data_combined <- read_tsv("source_fig_S1G.tsv")
  data_combined$sample <- factor(data_combined$sample, levels = c("GM12878_37", "GM12878_55", "3T3_37", "3T3_55"))
  data_combined$Type <- factor(data_combined$Type, levels = c("CDS Exons", "5\'UTR Exons", "3\'UTR Exons", "Introns", "Intergenic"))
  
  ggplot(data_combined, aes(x = sample, fill = Type, weight = proportion)) +
    theme_classic() +
    geom_bar(position = "stack") +
    ylab("Proportion %") +
    scale_fill_manual(values = c('#AF4304', '#F26D5B', '#FFBC42', '#D6ECFA', '#30A9DE')) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 13, color = c("black", "black")),
          axis.text.y = element_text(size = 13, color = "black"),
          axis.title.y = element_text(size = 15),
          axis.title.x = element_text(size = 15)) +
    coord_fixed(4/100)
}


# Fig S1H
{
  data <- read_tsv("source_fig_S1H.tsv")
  data$sample <- factor(data$sample, levels = c("GM12878_37", "GM12878_50", "3T3_37", "3T3_50"))
  
  ggplot(data = data, aes(x = reads, y = counts)) +
    theme_classic() +
    geom_smooth(aes(color = sample), formula = y~x+0, se = FALSE, fullrange = TRUE) +
    geom_point(aes(color = sample), size = 0.5, alpha = 1) +
    theme(axis.text.x = element_text(size = 13, color = c("black", "black")),
          axis.text.y = element_text(size = 13, color = "black"),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15)) +
    scale_color_manual(values = c('#0245a3', '#8fbaf3', '#a40a3c', '#ffc300')) +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab("Bowtie2 aligned reads per cell") +
    ylab("Unique fragments detected per cell") +
    ylim(c(0, 50000)) +
    xlim(c(0, 500000)) +
    guides(color = guide_legend(override.aes = list(size = 2))) +
    coord_fixed(ratio = 10) +
    facet_grid(.~sample)
}


# Fig S1J
{

  data <- read_tsv("source_fig_S1J.txt")
  data$reads <- factor(data$reads)
  data$dN_SMRT <- factor(data$dN_SMRT, levels =  c("No", "Yes"))
  
  out <- c()
  for(i in unique(data$reads)){
    data_used <- data %>% filter(reads == i)
    pout <- t.test(data_used$umis[data_used$dN_SMRT == "No"], data_used$umis[data_used$dN_SMRT == "Yes"])$p.value
    out <- c(out, pout)
  }
  
  color_JH8 <- '#d2ecf9'
  color_GJ <- '#fdb87d'
  
  ggplot(data, aes(x = reads, y = log10(umis))) +
    geom_boxplot(aes(fill = dN_SMRT), outlier.shape = NA, size = 0.5, alpha = 1) +
    coord_cartesian(ylim=c(2, 5)) +
    theme_classic() +
    scale_fill_manual(values = c(color_JH8, color_GJ)) +
    theme(axis.text.x = element_text(size = 13, color = c("black", "black")),
          axis.text.y = element_text(size = 13, color = "black"),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15)) +
    labs(title = "dN-SMRT") +
    xlab("Mean raw reads per cell") +
    ylab(expression(log[10](umis~per~cell))) +
    theme(plot.title = element_text(hjust = 0.5))
}

# Fig S1K
{
  
  data <- read_tsv("source_fig_S1K.txt")
  data$reads <- factor(data$reads)
  data$dN_SMRT <- factor(data$dN_SMRT, levels =  c("No", "Yes"))
  
  out <- c()
  for(i in unique(data$reads)){
    data_used <- data %>% filter(reads == i)
    pout <- t.test(data_used$genes[data_used$dN_SMRT == "No"], data_used$genes[data_used$dN_SMRT == "Yes"])$p.value
    out <- c(out, pout)
  }
  
  color_JH8 <- '#d2ecf9'
  color_GJ <- '#fdb87d'
  
  ggplot(data, aes(x = reads, y = log10(genes))) +
    geom_boxplot(aes(fill = dN_SMRT), outline = FALSE, outlier.shape = NA, size = 0.5, alpha = 1) +
    coord_cartesian(ylim=c(2, 4.2)) +
    theme_classic() +
    scale_fill_manual(values = c(color_JH8, color_GJ)) +
    theme(axis.text.x = element_text(size = 13, color = c("black", "black")),
          axis.text.y = element_text(size = 13, color = "black"),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15)) +
    labs(title = "dN-SMRT") +
    xlab("Mean raw reads per cell") +
    ylab(expression(log[10](genes~per~cell))) +
    theme(plot.title = element_text(hjust = 0.5))
}


# Fig S1L
{
  data <- read_tsv("source_fig_S1L.txt")
  
  ggplot(data = data, aes(x = reads, y = umis)) +
    theme_classic() +
    geom_point(aes(color = dN_SMRT), size = 1, alpha = 1, shape = 16) +
    geom_abline(slope = k_umis_JH8, intercept = 0, color = color_JH8, size = 1, linetype = 2) +
    geom_abline(slope = k_umis_GJ, intercept = 0, color = color_GJ, size = 1, linetype = 2) +
    scale_color_manual(values = c(color_JH8, color_GJ)) +
    theme(axis.text.x = element_text(size = 13, color = c("black", "black")),
          axis.text.y = element_text(size = 13, color = "black"),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15)) +
    labs(title = "dN-SMRT") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(limits = c(0, 50000)) +
    scale_y_continuous(limits = c(0, 30000)) +
    xlab("Raw reads per cell") +
    ylab("Umis detected per cell") +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    coord_fixed(ratio = 5/3)
}


# Fig S1M
{
  data <- read_tsv("source_fig_S1M.txt")
  
  ggplot(data = data, aes(x = reads, y = genes)) +
    theme_classic() +
    geom_point(aes(color = dN_SMRT), size = 1, alpha = 1, shape = 16) +
    # geom_abline(slope = k_umis_JH8, intercept = 0, color = color_JH8, size = 1, linetype = 2) +
    # geom_abline(slope = k_umis_GJ, intercept = 0, color = color_GJ, size = 1, linetype = 2) +
    scale_color_manual(values = c(color_JH8, color_GJ)) +
    theme(axis.text.x = element_text(size = 13, color = c("black", "black")),
          axis.text.y = element_text(size = 13, color = "black"),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15)) +
    labs(title = "dN-SMRT") +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab("Raw reads per cell") +
    ylab("Genes detected per cell") +
    ylim(c(0, 10000)) +
    xlim(c(0, 50000)) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    coord_fixed(ratio = 5)
}





