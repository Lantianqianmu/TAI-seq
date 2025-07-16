library(ggplot2)
library(viridis)
library(MASS)
# devtools::install_github("caleblareau/BuenColors")
library(BuenColors)
library(dplyr)
library(GenomicRanges)
library(SummarizedExperiment)
library(chromVAR)
library(parallel)
library(stringr)

count_matrix_from_fragments <- function(file, peaks, min_fragments = 1){
  # read fragment file
  fragments <- read.delim(file, header = FALSE, col.names = c("chr", "start", "end", "cell", "counts"), stringsAsFactors = FALSE)
  fragment_filtered <- fragments %>% group_by(cell) %>% filter(n()>=min_fragments)
  fragment_grouped <- fragment_filtered %>% group_split()
  
  results <- lapply(fragment_grouped, function(x, peaks) {
    fragment_ranges <- GRanges(x$chr, ranges = IRanges(x$start, x$end))
    return(list(counts = countOverlaps(peaks, fragment_ranges, type = "any", ignore.strand = TRUE), fragment_counts = length(fragment_ranges)))
  }, peaks = peaks)
  
  mat <- Matrix::Matrix(vapply(results, function(x) x[["counts"]], rep(0, length(peaks))))
  colnames(mat) <- (fragment_filtered %>% group_keys())$cell
  fragment_counts <- vapply(results, function(x) x[["fragment_counts"]], 0)
  
  colData <- data.frame(fragment_counts = fragment_counts)
  
  out <- SummarizedExperiment(assays = list(counts = mat), 
                              rowRanges = peaks, colData = colData)
  return(out)
}


read_peaks <- function(file, width = 250, non_overlapping = TRUE){
  cn <- c("chr", "start", "end", "name", "score", "strand", "fc", "pval", "qval", "summit")
  bed <- read.delim(file = file, header = FALSE, sep = "\t", stringsAsFactors = FALSE, col.names = cn)
  bed[, "summit"] <- bed[, "start"] + bed[, "summit"]
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


read_peaks_bed <- function(file, non_overlapping = TRUE){
  cn <- c("chr", "start", "end")
  bed <- read.delim(file = file, header = FALSE, sep = "\t", stringsAsFactors = FALSE, col.names = cn)
  bed <- as(bed, "DataFrame")
  
  bed <- makeGRangesFromDataFrame(bed, start.field = "start", end.field = "end", keep.extra.columns = TRUE)
  
  
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

solar_extra <- colorRampPalette(
  # blue - red
  c("#3361A5","#248AF3", "#14B3FF", "#88CEEF", "#C1D5DC", "#EAD397", "#FDB31A", "#E42A2A", "#A31D1D")
)

# remotes::install_local("/home/zeemeeuw/Setup/glossary-master.zip")
# remotes::install_local("/home/zeemeeuw/Setup/introdataviz-master.zip")

assignIdents <- function(x){
  human <- as.integer(x[2])
  mouse <- as.integer(x[3])
  if(human/mouse > 4){
    return("human")
  }else if(mouse/human > 4){
    return("mouse")
  }else{
    return("mix")
  }
}

{
# mix_assay <- function(file_hg38, file_mm39, min = 100){
#   hg38 <- read.delim(file_hg38, header = FALSE, stringsAsFactors = FALSE, col.names = c("barcode", "count.hg38"))
#   mm39 <- read.delim(file_mm39, header = FALSE, stringsAsFactors = FALSE, col.names = c("barcode", "count.mm39"))
#   data.all <- merge(hg38, mm39, by = 'barcode', all = TRUE)
#   data.all[is.na(data.all)] <- 0
#   data <- data.all[!(data.all$count.hg38 < min & data.all$count.mm39 < min),]
#   data$ident <- apply(data, 1, assignIdents)
#   data$ident <- factor(data$ident, levels = c("human", "mouse", "mix"))
#   colnames(data)[4] <- "cell type"
#   return(data)
# }
# 
# get_genes <- function(dir, min.cells = 3, min.features = 100){
#   library(Seurat)
#   mat <- Read10X(dir, gene.column = 2)
#   seu <- CreateSeuratObject(counts = mat, min.cells = min.cells, min.features = min.features)
#   n_genes <- as.numeric(apply(seu@assays$RNA@counts, 2, function(x) length(which(x != 0))))
#   return(n_genes)
# }
# 
# annotate_sublibrary <- function(seu, lib_AACCA = 'lib1', lib_CATAG = 'lib2', lib_GCACT = 'lib3'){
#   seu$barcode <- substr(colnames(seu), 1, 23)
#   seu$sub_library_barcode <- substr(colnames(seu), 25, 29)
#   seu$sub_library_annotation <- dplyr::recode(seu$sub_library_barcode, 'AACCA' = lib_AACCA, 'CATAG' = lib_CATAG, 'GCACT' = lib_GCACT)
#   seu$retain <- !duplicated(seu$barcode)
#   orig_length <- length(seu$retain)
#   new_length <- sum(seu$retain)
#   message(paste0("Found ", orig_length, " cells in the original Seurat object."))
#   message(paste0(orig_length-new_length, " cells dropped due to duplicated cell names."))
#   message("Final cell number: ", new_length, ".")
#   seu <- subset(seu, subset = retain)
#   seu@meta.data <- subset(seu@meta.data, select = -retain)
#   seu <- RenameCells(seu, new.names = seu$barcode)
#   return(seu)
# } 
}


runGenePeakcorr <- function(ATAC.se, # SummarizedExperiment object of scATAC data
                            RNAmat, # Paired normalized scRNA-seq data, with gene names as rownames
                            TSS, # Granges object
                            geneList = NULL, # 2 or more valid gene symbols (if only running on subset of genes)
                            windowPadSize = 50000, # base pairs padded on either side of gene TSS
                            normalizeATACmat = TRUE, # Whether or not to normalize scATAC counts (default is yes, assumes raw counts)
                            nCores = 4, # Number of cores if parallelization support
                            n_bg = 100, # Number of background peaks to use
                            genome = "mm10",
                            seed = 123
) {
  
  stopifnot(inherits(ATAC.se,"RangedSummarizedExperiment"))
  stopifnot(inherits(RNAmat,c("Matrix","matrix")))
  
  if(!all.equal(ncol(ATAC.se),ncol(RNAmat)))
    stop("Input ATAC and RNA objects must have same number of cells")
  
  message("Assuming paired scATAC/scRNA-seq data ..")
  
  peakRanges.OG <- granges(ATAC.se) # Peak ranges in reference input SE (pre-filtering)
  
  # Function needs rownames for both matrices or gives error
  rownames(ATAC.se) <- paste0("Peak",1:nrow(ATAC.se))
  ATACmat <- assay(ATAC.se) # Rownames preserved
  
  # Normalize peak counts
  if(normalizeATACmat)
    ATACmat <- centerCounts(ATACmat) # Rownames preserved
  
  if(is.null(rownames(RNAmat)))
    stop("RNA matrix must have gene names as rownames")
  
  # Check for peaks/genes with 0 accessibility/expression
  
  if(any(Matrix::rowSums(assay(ATAC.se))==0)){
    message("Peaks with 0 accessibility across cells exist ..")
    message("Removing these peaks prior to running correlations ..")
    peaksToKeep <- Matrix::rowSums(assay(ATAC.se))!=0
    ATAC.se <- ATAC.se[peaksToKeep,] # Subset ranges
    ATACmat <- ATACmat[peaksToKeep,]
    message("Important: peak indices in returned gene-peak maps are relative to original input SE")
  }
  
  
  peakRanges <- granges(ATAC.se) # Peak ranges
  
  if(any(Matrix::rowSums(RNAmat)==0)){
    message("Genes with 0 expression across cells exist ..")
    message("Removing these genes prior to running correlations ..")
    genesToKeep <- Matrix::rowSums(RNAmat)!=0
    RNAmat <- RNAmat[genesToKeep,]
  }
  
  cat("Number of peaks in ATAC data:",nrow(ATACmat),"\n")
  cat("Number of genes in RNA data:",nrow(RNAmat),"\n")
  
  # load TSS
  TSSg <- TSS
  
  # Keep genes that have annotation and are in RNA matrix
  names(TSSg) <- as.character(TSSg$gene_name)
  
  if(!is.null(geneList)){
    if(length(geneList)==1)
      stop("Please specify more than 1 valid gene symbol")
    
    if(any(!geneList %in% names(TSSg))){
      cat("One or more of the gene names supplied is not present in the TSS annotation specified: \n")
      cat(geneList[!geneList %in% names(TSSg)], sep = ", ")
      cat("\n")
      stop()
    }
    
    TSSg <- TSSg[geneList]
  }
  
  # Checking in case some genes in RNA don't overlap our TSS annotations
  genesToKeep <- intersect(names(TSSg), rownames(RNAmat))
  cat("\nNum genes overlapping TSS annotation and RNA matrix being considered: ",length(genesToKeep),"\n")
  
  # Match gene order in RNA matrix and TSS ranges
  RNAmat <- RNAmat[genesToKeep,]
  TSSg <- TSSg[genesToKeep]
  
  # Pad TSS by this much *either side*
  TSSflank <- GenomicRanges::flank(TSSg, width = windowPadSize, both = TRUE)
  
  # Get peak summit
  cat("\nTaking peak summits from peak windows ..\n")
  peakSummits <- resize(peakRanges, width = 1, fix = "center")
  
  # Find overlap of all peaks to all genes given window
  # Subject is Peaks, query is Gene
  cat("Finding overlapping peak-gene pairs ..\n")
  genePeakOv <- findOverlaps(query = TSSflank, subject = peakSummits)
  numPairs <- length(genePeakOv)
  
  cat("Found ", numPairs, "total gene-peak pairs for given TSS window ..\n")
  
  cat("Number of peak summits that overlap any gene TSS window: ", length(unique(subjectHits(genePeakOv))),"\n")
  cat("Number of gene TSS windows that overlap any peak summit: ", length(unique(queryHits(genePeakOv))),"\n\n")
  
  # For each gene, determine observed correlation of each overlapping peak to its associated gene (gene expression)
  
  # For each of those genes, also determine correlation based on background peaks (run in parallel) and save
  # Do this in parallel, and get data frame of gene-peak-pearson values
  # Fetch background peaks for each peak tested (i.e. that has overlap in window with gene)
  set.seed(123)
  cat("Determining background peaks ..\n")
  
  if(is.null(rowData(ATAC.se)$bias)){
    if(genome == "hg19")
      myGenome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    if(genome == "mm10")
      myGenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
    if(genome == "hg38")
      myGenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    
    ATAC.se <- chromVAR::addGCBias(ATAC.se, genome = myGenome) }
  
  cat("Using ", n_bg, " iterations ..\n\n")
  
  set.seed(seed)
  bg <- chromVAR::getBackgroundPeaks(ATAC.se, niterations = n_bg)
  
  cat("Computing gene-peak correlations ..\n")
  
  pairsPerChunk <- 1000
  
  # This defines the outer (larger chunks)
  largeChunkSize <- 10000
  
  startingPoint <- 1 # If for any reason the workers fail, resume from where it failed by specifying the starting point here
  chunkStarts <- seq(startingPoint, numPairs, largeChunkSize)
  chunkEnds <- chunkStarts + largeChunkSize -1
  chunkEnds[length(chunkEnds)] <- numPairs
  
  library(doParallel)
  
  dorcList <- list()
  for(i in 1:length(chunkStarts)){
    cat("Running pairs: ",chunkStarts[i], "to", chunkEnds[i], "\n")
    # This fill further chunk this up and run in parallel, saving the merged output ObsCor
    ObsCor <- PeakGeneCor(ATAC = ATACmat,
                          RNA = RNAmat,
                          OV = genePeakOv[chunkStarts[i]:chunkEnds[i]],
                          chunkSize = pairsPerChunk,
                          ncores = nCores,
                          bg = bg)
    gc()
    dorcList[[i]] <- ObsCor
  }
  
  cat("\nMerging results ..\n")
  dorcTab <- bind_rows(dorcList)
  
  cat("Performing Z-test for correlation significance ..\n")
  permCols <- 4:(ncol(bg)+3)
  
  
  # Swap gene number for gene symbol from TSS annotation lookup
  dorcTab$Gene <- as.character(TSSg$gene_name)[dorcTab$Gene]
  
  # Swap peak numbers to match reference input peak numbers
  # This only changes if some peaks had zero accessibility and were filtered out internally
  # Use rownames from reference matching
  dorcTab$Peak <- as.numeric(sapply(strsplit(as.character(rownames(ATACmat)[dorcTab$Peak]), "Peak", fixed=TRUE), "[[", 2))
  
  # # Z test pval
  dorcTab$rBgSD <- matrixStats::rowSds(as.matrix(dorcTab[,permCols]))
  dorcTab$rBgMean <- rowMeans(dorcTab[,permCols])
  dorcTab$pvalZ <- 1 - stats::pnorm(q = dorcTab$rObs, mean = dorcTab$rBgMean, sd = dorcTab$rBgSD)
  
  message("\nPeak-gene correlation calculation is finished!\n")

  
  # Add peak ranges to final data frame output
  dorcTab$PeakRanges <- paste(as.character(seqnames(peakRanges.OG[dorcTab$Peak])),paste(start(peakRanges.OG[dorcTab$Peak]),end(peakRanges.OG[dorcTab$Peak]),sep="-"),sep=":")
  
  return(as.data.frame(dorcTab[,c("Peak","PeakRanges","Gene","rObs","pvalZ")],stringsAsFactors=FALSE))
}



reduceMultiMappingPeaks <- function(dorcTab){
  message("Keeping max correlation for multi-mapping peaks ..\n")
  dorcTab <- dorcTab %>% dplyr::group_by(Peak) %>% dplyr::filter(rObs==max(rObs))
  return(dorcTab)
}


centerCounts <- function(obj,
                         doInChunks = TRUE,
                         chunkSize = 1000){
  # Check if SE or Matrix
  # Added to avoid error https://github.com/buenrostrolab/FigR/issues/16
  if(any(sapply(c("SummarizedExperiment", "RangedSummarizedExperiment"), function(x){inherits(obj,x)}))){
    cat("SummarizedExperiment object input detected .. Centering counts under assay")
    isSE <- TRUE
  } else {
    if(any(sapply(c("dgCMatrix", "dgeMatrix", "Matrix"),function(x){inherits(obj,x)}))){
      cat("Matrix object input detected")
      isSE <- FALSE
    } else {
      stop("Supplied object must be either of class SummarizedExperiment or Matrix ..\n")
    }
  }
  
  #if(!class(obj) %in% c("SummarizedExperiment","RangedSummarizedExperiment","dgCMatrix","dgeMatrix","Matrix"))
  
  
  if(ncol(obj) > 10000)
    doInChunks <- TRUE
  
  if(doInChunks){
    cat("Centering counts for cells sequentially in groups of size ",
        chunkSize, " ..\n\n")
    starts <- seq(1,ncol(obj),chunkSize)
  } else{
    starts <- 1
  }
  
  counts.l <- list()
  
  for(i in 1:length(starts)){
    beginning <- starts[i]
    if(i==length(starts)) # If it's the last one
    {
      ending <- ncol(obj)
    } else {
      ending <- starts[i]+chunkSize-1
    }
    
    cat("Computing centered counts for cells: ", beginning, " to ", ending, "..\n")
    
    #if(class(obj) == "RangedSummarizedExperiment" | class(obj)=="SummarizedExperiment"){
    if(isSE){
      m <- SummarizedExperiment::assay(obj[, beginning:ending])} else {
        m <- obj[,beginning:ending] # Assumes Matrix format
      }
    cellMeans <- Matrix::colMeans(m)
    cat("Computing centered counts per cell using mean reads in features ..\n\n")
    # Center cell counts based on its mean RIP count
    cCounts <- Matrix::t(Matrix::t(m)/cellMeans)
    
    counts.l[[i]] <- cCounts
    
    gc()
  }
  
  cat("Merging results..\n")
  centered.counts <- do.call("cbind",counts.l)
  cat("Done!\n")
  
  #if(class(obj) == "RangedSummarizedExperiment" | class(obj)=="SummarizedExperiment"){
  if(isSE){
    SummarizedExperiment::assay(obj) <- centered.counts
    return(obj)
  } else {
    return(centered.counts)
  }
}


PeakGeneCor <- function(ATAC, # Normalized reads in peaks counts (rownames should  be "Peak1","Peak2" etc.)
                        RNA, # Normalized gene expression counts
                        OV, # Gene TSS - Peak overlap pairs object (Genes: query, Peaks: subject)
                        ncores = 4,
                        chunkSize = 200,
                        metric = "spearman",
                        bg = NULL){
  
  stopifnot(ncol(ATAC)==ncol(RNA))
  
  if(chunkSize > 1000)
    stop("Do not specify very large chunk sizes. Please use chunkSize <= 1000")
  
  
  # Number of total gene-peak pairs to chunk up for parallelization
  n <- length(OV)
  starts <- seq(1, n, chunkSize)
  ends <- starts + chunkSize -1
  ends[length(ends)] <- n
  
  OVd <- OV %>% as.data.frame() %>% dplyr::rename("Gene"="queryHits","Peak"="subjectHits")
  
  chunkList <- mapply(c, starts, ends, SIMPLIFY = FALSE)
  
  time_elapsed <- Sys.time()
  
  cat("Running in parallel using ", ncores, "cores ..\n")
  
  cat("Computing observed correlations ..\n")
  
  corList <- parallel::mclapply(X = chunkList, function(x) {chunkCore(chunk=x,A=ATAC,R=RNA,O=OVd,met=metric)}, mc.cores = ncores)
  
  
  if(any(unlist(sapply(corList,is.null)))){
    message("One or more of the chunk processes failed unexpectedly (returned NULL) ..")
    message("Please check to see you have enough cores/memory allocated")
    message("Also make sure you have filtered down to non-zero peaks/genes")
  }
  
  OVd$rObs <- unlist(corList)
  
  
  time_elapsed <- Sys.time() - time_elapsed
  cat(paste("\nTime Elapsed: ",time_elapsed, units(time_elapsed)),"\n\n")
  
  if(!is.null(bg)){
    n_iter <- ncol(bg)
    cat("Computing background correlations ..\n")
    
    time_elapsed <- Sys.time()
    
    bgCor <- foreach(i= 1:n_iter, .combine = 'cbind', .export = c("chunkCore","t"), .packages = c("parallel","Matrix")) %do% {
                       OVdBg <- OVd[,1:2] # Initialize gene-peak pairing to observed
                       OVdBg$Peak <- bg[OVdBg$Peak,i] # Swap actual peaks with bg peaks for given iteration in pairing
                       bgCorList <- parallel::mclapply(X = chunkList, function(x) {chunkCore(chunk = x,A = ATAC,R = RNA,O = OVdBg,met = metric)}, mc.cores = ncores)
                       unlist(bgCorList) # Vector of background permuted correlation values for that set of background peaks
                     }
    
    
    if(sum(is.null(bgCor))!=0 | sum(is.na(bgCor))!=0)
      stop("One or more of the chunk processes failed unexpectedly (returned NULL) .. Please check to see you have enough cores/m
           emory allocated")
    
    time_elapsed <- Sys.time() - time_elapsed
    cat(paste("\nTime Elapsed: ",time_elapsed, units(time_elapsed)),"\n\n")
    
    colnames(bgCor) <- paste0("rBg",1:ncol(bgCor))
    OVd <- cbind(OVd,bgCor)
  }
  
  return(OVd)
}


chunkCore <- function(chunk,
                      A, # ATAC matrix
                      R, # RNA matrix
                      O, # Gene-Peak overlap pairing data.frame
                      met # Correlation method ("spearman" or "pearson")
){
  # Get indices of genes and peaks from overlap object for chunk
  # Assumes query hits are genes and subject hits are peaks in the overlap object
  geneIndices <- O$Gene[chunk[1]:chunk[2]]
  peakIndices <- O$Peak[chunk[1]:chunk[2]]
  
  pairnames <- cbind(rownames(A)[peakIndices],rownames(R)[geneIndices])
  
  uniquegenes <- unique(geneIndices)
  uniquepeaks <- unique(peakIndices)
  
  
  M1 <- as.matrix(Matrix::t(A[uniquepeaks,,drop=FALSE])) # In case only 1 match, keep matrix structure
  M2 <- as.matrix(Matrix::t(R[uniquegenes,,drop=FALSE])) # In case only 1 match, keep matrix structure
  
  # Peak x Gene correlation matrix, subset by peak-gene pair names to get corresponding correlation vector
  # NOTE: This correlation call fails if you have maps with just 1 gene / peak. This is unlikely for large chunk sizes
  cor(x = M1,y = M2, method = met)[pairnames]
  
}



dorcJPlot <- function(dorcTab,
                      cutoff = 7, # Eventually we should use a knee caller for this
                      labelTop = 25,
                      returnGeneList = FALSE,
                      cleanLabels = TRUE,
                      labelSize = 4,
                      ...){
  
  stopifnot(all(c("Peak", "Gene", "pvalZ") %in% colnames(dorcTab)))
  
  # Count the number of significant peak associations for each gene (without pre-filtering genes)
  numDorcs <- dorcTab  %>% group_by(Gene) %>% tally() %>% arrange(dplyr::desc(n))
  numDorcs$Index <- 1:nrow(numDorcs) # Add order index
  numDorcs %>% as.data.frame(stringsAsFactors=FALSE) -> numDorcs
  rownames(numDorcs) <- numDorcs$Gene
  
  dorcGenes <- numDorcs$Gene[numDorcs$n >= cutoff]
  
  numDorcs <- numDorcs %>%
    mutate(isDORC=ifelse(Gene %in% dorcGenes, "Yes", "No")) %>%
    mutate(Label=ifelse(Gene %in% dorcGenes[1:labelTop], Gene, ""))
  
  # Plot
  dorcG <- ggplot(numDorcs, aes(x = Index,y = n, color = isDORC, label = Label)) +
    geom_hline(linetype = "dotted", yintercept = cutoff)+
    geom_vline(linetype = "dotted", xintercept = max(numDorcs[numDorcs$Gene %in% dorcGenes, "Index"]))+
    geom_point(size = 0.8) +
    geom_line() +
    scale_color_manual(values = c("gray65", "firebrick"))+
    scale_y_continuous(breaks = scales::pretty_breaks())+
    theme_classic() +
    labs(y = "Number of correlated peaks", x = "Ranked genes", title = paste0("# DORCs: (n >= ", cutoff, ") = ", length(dorcGenes))) +
    theme(axis.text = element_text(color = "black"), legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    scale_x_reverse() # flip so we can add labels later, if needed, with more space
  
  if(labelTop != 0){
    if(cleanLabels){
      dorcG <- dorcG + ggrepel::geom_label_repel(size = labelSize, max.iter = 100, max.overlaps = Inf, fontface = "italic",...)
    } else {
      dorcG <- dorcG + ggplot2::geom_text(size = labelSize, fontface = "italic",...)
    }
  }

  
  print(dorcG)
  
  if(returnGeneList)
    return(dorcGenes)
  
}




getDORCScores <- function(ATAC.se,
                          dorcTab,
                          normalizeATACmat = TRUE,
                          geneList = NULL,
                          nCores = 4){
  if(!all(c("Peak","Gene") %in% colnames(dorcTab)))
    stop("The provided gene-peak table must have columns named Peak and Gene ..")
  
  if(any(dorcTab$Peak > nrow(ATAC.se)))
    stop("One or more peak indices in the gene-peak table are larger than the total number of peaks in the provided ATAC SE object ..\n Make sure the exact same SummarizedExperiment object is provided here as was for running the runGenePeakcorr function ..\n")
  
  if(!is.null(geneList)){
    if(!(all(geneList %in% as.character(dorcTab$Gene))))
      stop("One or more of the gene names supplied is not present in the gene-peak table provided..\n")
    
    if(length(geneList) > 50){
      message("Running DORC scoring for ", length(geneList), " genes: ", paste(geneList[1:20],collapse=", "),", ... , ... , ... (truncated display)")
    }else{
      message("Running DORC scoring for ", length(geneList), " genes: ", paste(geneList,collapse = "\n"))
    }
    
    cat("........\n")
    
    dorcTab <- dorcTab[dorcTab$Gene %in% geneList,] # Filter only to these genes
    dorcGenes <- sort(as.character(unique(dorcTab$Gene)))
  }else{
    dorcGenes <- sort(as.character(unique(dorcTab$Gene)))
    cat("Running DORC scoring for all genes in annotation! (n = ", length(dorcGenes), ")\n", sep="")
  }
  
  if(normalizeATACmat){
    # Normalize
    cat("Normalizing scATAC counts ..\n")
    ATAC.mat <- assay(centerCounts(ATAC.se, chunkSize = 5000))
    gc()
  }else{
    cat("Assuming provided scATAC counts are normalized ..\n")
    ATAC.mat <- assay(ATAC.se)
  }
  
  time_elapsed <- Sys.time()
  
  if(Sys.info()['sysname'] %in% "Windows"){
    message("Windows OS detected .. Cannot support parallilzation using mclapply for mc.cores > 1")
    message("Using 1 core instead ..\n")
    nCores <- 1
  }
  
  cat("Computing DORC scores ..\n")
  cat("Running in parallel using ", nCores, "cores ..\n")
  
  dorcMatL <- parallel::mclapply(X = dorcGenes, function(x) {
    dorcPeaks <- unique(dorcTab$Peak[dorcTab$Gene %in% x])
    if(length(dorcPeaks) > 1){
      dorcCounts <- Matrix::colSums(ATAC.mat[dorcPeaks,])
    } else if(length(dorcPeaks == 1)){
      dorcCounts <- ATAC.mat[dorcPeaks,]
    }
  }, mc.cores = nCores)
  
  dorcMat <- Matrix::Matrix(do.call('rbind', dorcMatL), sparse = TRUE)
  
  rownames(dorcMat) <- dorcGenes
  
  time_elapsed <- Sys.time() - time_elapsed
  cat(paste("\nTime Elapsed: ", time_elapsed, units(time_elapsed)), "\n\n")
  
  return(dorcMat)
  
}


choose_colorset <- function(colorset, n = 256){
  if(colorset == "solar_extra"){
    # colors <- ArchR::paletteContinuous(set = "solarExtra", n = 256, reverse = FALSE)
    colors <- colorRampPalette(colors = c("#3361A5", "#248AF3", "#14B3FF", "#88CEEF", "#C1D5DC", "#EAD397", "#FDB31A", "#E42A2A", "#A31D1D"))(n)
  }else if(colorset == "viridis"){
    colors <- viridis_pal(alpha = 1, begin = 0, end = 1, direction = 1, option = "D")(n)
  }else if(colorset == "greens"){
    colors <- colorRampPalette(colors = c("#f0f0f0", RColorBrewer::brewer.pal(9, "Greens")[9]))(n)
  }else if(colorset == "cosmic"){
    colors <- colorRampPalette(colors = c('#000000', '#010101', '#040305', '#0A060E', '#100916', '#160C1E', '#1C1027', '#211230', '#28153B', '#2E1745', '#341950', '#3A1A5B', '#3F1A67', '#461977', '#4C1885', '#521594', '#5710A4', '#5C09B5', '#6003C8', '#610AD5', '#6019DF', '#5D28E4', '#5936E7', '#5345E7', '#4E50E6', '#495AE5', '#4463E3', '#3E6DE1', '#3975DF', '#347CDD', '#3083DC', '#2C8ADA', '#2791D9', '#2498D8', '#209ED8', '#1DA4D7', '#1AAAD7', '#17B2D7', '#13B8D7', '#10BED7', '#0BC4D7', '#07CBD7', '#03D2D7', '#03D9D6', '#07DFD6', '#11E5D5', '#1DECD4', '#30F3D2', '#44F9D1', '#5FFDD0'))(n)
  }else if(colorset == "yp"){
    colors <- colorRampPalette(c("#60339B", "#805280", "#A07265", "#BF914A", "#DFB12F", "#FFD014"))(n)
  }else if(colorset == "imagine"){
    colors <- colorRampPalette(c("#1f6933", "#2e8c59", "#35a367", "#41b675", "#66c992", "#86d9ac",
                                 "#d8dcc1",
                                 "#FFC782", "#FFA12F", "#FF851D", "#ff6200", "#FF4E0A", "#CF3D00"))(n)
  }else if(colorset == "sunburst"){
    colors <- colorRampPalette(c('#FFFFFF', '#F1F1E6', '#E8E4CB', '#E2D5AD', '#E0C490', '#E0B275', '#E0A15E', '#DF8E47', '#DE7B32', '#DB661E', '#D8530F', '#D23C06', '#C92310', '#BB111E', '#AA0E29', '#961130', '#821434', '#6E1635', '#5C1633', '#49152F'))(n)
  }else if(colorset == "gates_of_heaven"){
    colors <- colorRampPalette(c("#85B0FF", "#98C4F2", "#ABD8E6", "#BEEBD9", "#D1FFCC", "#B5F5B3"))(n)
  }else if(colorset == "romance"){
    colors <- colorRampPalette(rev(c("#B22222", "#E25822", "#F1BC31", "#F6F052", "#C09ADB", "#A17BB9", "#825B97", "#623C74")))(n)
  }else if(colorset == "greyred"){
    colors <- colorRampPalette(c("#888888" ,"#989898", "#b8b8b8" ,"#d6d6d6", "#F6F052", "#F1BC31", "#E25822", "#B22222"))(n)
  }else if(colorset == "gothic"){
    colors <- colorRampPalette(c('#FFFFFF', '#EDF1FA', '#DDE4F8', '#CFD5F7', '#CAC3F2', '#C9B2E6', '#C8A2DB', '#C88FCD', '#C97BC1', '#CA64B7', '#C84DB8', '#BA3CC0', '#AA2EC5', '#9820C9', '#8811CB', '#7505C8', '#5F12B9', '#49219B', '#39247C', '#2B225D'))(n)
  }else if(colorset == "flamingo"){
    colors <- colorRampPalette(c('#FFFFFF', '#F5EEF4', '#EEDEEC', '#E9CCE3', '#E6BADA', '#E4A7D0', '#E395C6', '#E280BA', '#E16AAC', '#DF539A', '#DC3C87', '#D6226F', '#CB0755', '#BB023D', '#AA0E2D', '#951820', '#811D18', '#6D1E12', '#5B1D0F', '#491A0C'))(n)
  }else if(colorset == "tropical_ocean"){
    colors <- colorRampPalette(rev(c("#1831F5", "#3971FF", "#59C0CB", "#7ED6B7", "#EBE5E1", "#EFF4F3")))(n)
  }else if(colorset == "solar_flare"){
    colors <- colorRampPalette(c("#3361A5", "#2884E7", "#1BA7FF", "#76CEFF", "#FFFFFF", "#FFE060", "#FA8E24", "#DA2828", "#A31D1D"))(n)
  }else if(colorset == "mystic_sea"){
    colors <- colorRampPalette(rev(c("#12355C", "#1F5A79", "#2D83A1", "#61D2D1", "#92E0D0", "#F4FCFA")))(n)
  }else if(colorset == "hot"){
    colors <- colorRampPalette(c("#FFFFAA", "#FFEC77", "#FFB945", "#FE8902", "#C55A00", "#8E2C00"))(n)
  }else if(colorset == "ember"){
    colors <- colorRampPalette(c('#F1E13E', '#F3CE30', '#F4BD23', '#F3AB16', '#F2980A', '#EF8602', '#EC7601', '#E76406', '#E1510F', '#DA3F18', '#D22D21', '#C61A2C', '#B70C36', '#A60B3F', '#951043', '#821646', '#701945', '#5D1A42', '#4C1A3D', '#3B1835'))(n)
  }else if(colorset == "waterlily"){
    colors <- colorRampPalette(c('#E8F2F2', '#D1E5E7', '#B9DADD', '#A0CFD5', '#88C5D0', '#6FB9CD', '#5CADCA', '#539FC6', '#5090C0', '#5081B8', '#5073B0', '#4F64A8', '#4F55A0', '#4E4597', '#4D368E', '#4B2484', '#461076', '#3C045E', '#2C0347'))(n)
  }else if(colorset == "reds"){
    colors <- colorRampPalette(c("#deddda", "#d4c5c3", "#caadab", "#c09694", "#b67e7c", "#ad6665", "#a34e4d", "#993736", "#8f1f1e", "#850707"))(n)
  }else if(colorset == "oranges"){
    colors <- colorRampPalette(c("#FAF8EE", "#FBE6BF", "#FCD48F", "#FDC360", "#FEB130", "#FF9F01"))(n)
  }else if(colorset == "holy"){
    colors <- colorRampPalette(c("#f0f0f0", "#F9CC29", "#EE8300", "#EE8300", "#EE3B00", "#C80038"))(n)
  }else if(colorset == "evercode"){
    colors <- colorRampPalette(rev(c("#AC1926FF", "#B52027FF", "#BD2629FF", "#C83F01FF", "#E5570AFF", "#F67524FF", "#FD9346FF", "#F8B075FF", "#CACACAFF")))(n)
  }else if(colorset == "dream"){
    colors <- colorRampPalette(c("#4292c9", "#a0c9e5", "#35a153", "#afdd8b", "#f26a11", "#fe9376", "#817ab9", "#bcbddd"))(n)
  }else if(colorset == "bubble"){
    colors <- colorRampPalette(c("#0ddbf5", "#1d9bf7", "#8386fc", "#303cf9", "#fe5357", "#fd7c1a", "#ffbd15", "#fcff07"))(n)
  }else if(colorset == "RdBu"){
    colors <- colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(n)
  }else if(colorset == "colorful_1"){
    colors <- colorRampPalette(c("#DC050A", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", "#B17BA6", "#FF7F00", "#DFB462", "#E7298A", "#E78AC3", "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7"))(n)
  }else if(colorset == "solar_flare"){
    colors <-colorRampPalette(c("#3361A5", "#2884E7", "#1BA7FF", "#76CEFF", "#FFFFFF", "#FFE060", "#FA8E24", "#DA2828", "#A31D1D"))(n)
  }else if(colorset == "ocean_brick"){
    colors <-colorRampPalette(c("#0F3341", "#1563AA", "#0B99E6", "#3DCDFD", "#F7F7F7", "#EB9457", "#D1551F", "#B02F1B", "#8D1616"))(n)
  }else if(colorset == "calma_bosque"){
    colors <- colorRampPalette(c("#F2FEF3", "#B6EBBA", "#86D88B", "#56C35D", "#36A23D", "#1C7722", "#094D0E", "#032506"))(n)
  }else if(colorset == "calma_manudo"){
    colors <- colorRampPalette(c("#FFEEEE", "#F7B4B4", "#ED8080", "#DF4A4A", "#BE2A2A", "#8C1616", "#590707", "#290303"))(n)
  }else if(colorset == "calma_azules"){
    colors <- colorRampPalette(c("#F2FBFE", "#B6DDEB", "#86C2D8", "#56A6C3", "#3685A2", "#1C5F77", "#093B4D", "#031C25"))(n)
  }else if(colorset == "calma_morado"){
    colors <- colorRampPalette(c("#F2F4FE", "#B6BFEB", "#8694D8", "#5668C3", "#3648A2", "#1C2B77", "#09154D", "#030925"))(n)
  }else if(colorset == "calma_marino"){
    colors <-colorRampPalette(c("#F2FEF8", "#B6EBD2", "#86D8B2", "#56C390", "#36A26F", "#1C774D", "#094D2D", "#032515"))(n)
  }else if(colorset == "calma_musgos"){
    colors <- colorRampPalette(c("#FCFEF2", "#E4EBB6", "#CDD886", "#B4C356", "#93A236", "#6B771C", "#444D09", "#212503"))(n)
  }else if(colorset == "wolfgang_extra"){
    colors <- colorRampPalette(c("#FFFFFF", "#FCFED3", "#E3F4B1", "#ABDEB6", "#60C1BF", "#2A9EC1", "#206AAD", "#243996", "#081D58"))(n)
  }else if(colorset == "wolfgang_basic"){
    colors <- colorRampPalette(c("#FFFFD9", "#EDF8B1", "#C7E9B4", "#7FCDBB", "#41B6C4", "#1D91C0", "#225EA8", "#253494", "#081D58"))(n)
  }else if(colorset == "comet"){
    colors <- colorRampPalette(c("#E6E7E8", "#3A97FF", "#8816A7", "black"))(n)
  }else if(colorset == "white_purple"){
    colors <- colorRampPalette(c("#f7fcfd", "#e0ecf4", "#bfd3e6", "#9ebcda", "#8c96c6", "#8c6bb1", "#88419d", "#810f7c", "#4d004b"))(n)
  }else if(colorset == "white_blue"){
    colors <- colorRampPalette(c("#fff7fb", "#ece7f2", "#d0d1e6", "#a6bddb", "#74a9cf", "#3690c0", "#0570b0", "#045a8d", "#023858"))(n)
  }else if(colorset == "green_blue"){
    colors <- colorRampPalette(c("#e0f3db", "#ccebc5", "#a8ddb5", "#4eb3d3", "#2b8cbe", "#0868ac", "#084081"))(n)
  }else if(colorset == "purple_orange"){
    colors <- colorRampPalette(c("#581845", "#900C3F", "#C70039", "#FF5744", "#FFC30F"))(n)
  }else if(colorset == "horizon"){
    colors <- colorRampPalette(c("#000075", "#2E00FF", "#9408F7", "#C729D6", "#FA4AB5", "#FF6A95", "#FF8B74", "#FFAC53", "#FFCD32", "#FFFF60" ))(n)
  }else if(colorset == "horizon_extra"){
    colors <- colorRampPalette(c("#000436", "#021EA9", "#1632FB", "#6E34FC", "#C732D5", "#FD619D", "#FF9965", "#FFD32B", "#FFFC5A"))(n)
  }else if(colorset == "fireworks"){
    colors <- colorRampPalette(c("white", "#2488F0", "#7F3F98", "#E22929", "#FCB31A"))(n)
  }else if(colorset == "fireworks2"){
    colors <- colorRampPalette(c("black", "#2488F0", "#7F3F98", "#E22929", "#FCB31A"))(n)
  }else if(colorset == "blue_yellow"){
    colors <- colorRampPalette(c("#352A86", "#343DAE", "#0262E0", "#1389D2", "#2DB7A3", "#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D"))(n)
  }else if(colorset == "grey_magma"){
    colors <- colorRampPalette(c("grey", "#FB8861FF", "#B63679FF", "#51127CFF", "#000004FF" ))(n)
  }else if(colorset == "beach"){
    colors <- colorRampPalette(c("#87D2DB", "#5BB1CB", "#4F66AF", "#F15F30", "#F7962E", "#FCEE2B"))(n)
  }else if(colorset == "samba_night"){
    colors <- colorRampPalette(c("#1873CC", "#1798E5", "#00BFFF", "#4AC596", "#00CC00", "#A2E700", "#FFFF00", "#FFD200", "#FFA500"))(n)
  }else if(colorset == "nostelgia"){
    colors <- colorRampPalette(c("#4858A7", "#788FC8", "#D6DAE1", "#F49B7C", "#B51F29"))(n)
  }else if(colorset == "captain"){
    colors <- colorRampPalette(c("grey", "#A1CDE1", "#12477C", "#EC9274", "#67001E"))(n)
  }else if(colorset == "rushmore"){
    colors <- colorRampPalette(c("#E1BD6D", "#EABE94", "#0B775E", "#35274A", "#F2300F"))(n)
  }else if(colorset == "darjeeling"){
    colors <- colorRampPalette(c("#FF0000", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6"))(n)
  }else if(colorset == "stallion"){
    colors <- colorRampPalette(c("#D51F26", "#272E6A", "#208A42", "#89288F", "#F47D2B", "#FEE500", "#8A9FD1", "#C06CAB", "#E6C2DC", "#90D5E4",
                                 "#89C75F", "#F37B7D", "#9983BD", "#D24B27", "#3BBCA8", "#6E4B9E", "#0C727C", "#7E1416", "#D8A767", "#3D3D3D"))(n)
  }else if(colorset == "stallion2"){
    colors <- colorRampPalette(c("#D51F26", "#272E6A", "#208A42", "#89288F", "#F47D2B", "#FEE500", "#8A9FD1", "#C06CAB", "#E6C2DC", "#90D5E4",
                                 "#89C75F", "#F37B7D", "#9983BD", "#D24B27", "#3BBCA8", "#6E4B9E", "#0C727C", "#7E1416", "#D8A767"))(n)
  }else if(colorset == "calm"){
    colors <- colorRampPalette(c("#7DD06F", "#844081", "#688EC1", "#C17E73", "#484125", "#6CD3A7", "#597873", "#7B6FD0", "#CF4A31", "#D0CD47",
                                 "#722A2D", "#CBC594", "#D19EC4", "#5A7E36", "#D4477D", "#403552", "#76D73C", "#96CED5", "#CE54D1", "#C48736" ))(n)
  }else if(colorset == "kelly"){
    colors <- colorRampPalette(c("#FFB300", "#803E75", "#FF6800", "#A6BDD7", "#C10020", "#CEA262", "#817066", "#007D34", "#F6768E", "#00538A",
                                 "#FF7A5C", "#53377A", "#FF8E00", "#B32851", "#F4C800", "#7F180D", "#93AA00", "#593315", "#F13A13", "#232C16"))(n)
  }else if(colorset == "paired"){
    colors <- colorRampPalette(c("#A6CDE2", "#1E78B4", "#74C476", "#34A047", "#F59899", "#E11E26", "#FCBF6E", "#F47E1F", "#CAB2D6", "#6A3E98", "#FAF39B", "#B15928" ))(n)
  }else if(colorset == "grove"){
    colors <- colorRampPalette(c("#1a1334", "#01545a", "#017351", "#03c383", "#aad962", "#fbbf45", "#ef6a32", "#ed0345", "#a12a5e", "#710162", "#3B9AB2"))(n)
  }else if(colorset == "bear"){
    colors <- colorRampPalette(c("#faa818", "#41a30d", "#fbdf72", "#367d7d", "#d33502", "#6ebcbc", "#37526d", "#916848", "#f5b390", "#342739",
                                 "#bed678", "#a6d9ee", "#0d74b6", "#60824f", "#725ca5", "#e0598b"))(n)
  }else if(colorset == "ironman"){
    colors <- colorRampPalette(c("#371377", "#7700FF", "#9E0142", "#FF0080", "#DC494C", "#F88D51", "#FAD510", "#FFFF5F", "#88CFA4", "#238B45",
                                 "#02401B", "#0AD7D3", "#046C9A", "#A2A475", "grey35"))(n)
  }else if(colorset == "circus"){
    colors <- colorRampPalette(c("#D52126", "#88CCEE", "#FEE52C", "#117733", "#CC61B0", "#99C945", "#2F8AC4", "#332288", "#E68316", "#661101",
                                 "#F97B72", "#DDCC77", "#11A579", "#89288F", "#E73F74"))(n)
  }else if(colorset == "summer_night"){
    colors <- colorRampPalette(c("#2a7185", "#a64027", "#fbdf72", "#60824f", "#9cdff0", "#022336", "#725ca5" ))(n)
  }else if(colorset == "zissou"){
    colors <- colorRampPalette(c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00"))(n)
  }else if(colorset == "solar_dark"){
    colors <- colorRampPalette(c("#9C0824FF", "#A70C20FF", "#B3101BFF", "#BE1316FF", 
                                 "#C51517FF", "#CB1618FF", "#D21E1CFF", "#D73529FF", "#DC4636FF", 
                                 "#E25644FF", "#EC6857FF", "#F6796AFF", "#F78E80FF", "#E9A79DFF", 
                                 "#D6BFBBFF", "#BBC5CCFF", "#9CBBCFFF", "#78B1D3FF", "#5EA5CEFF", 
                                 "#4F98C4FF", "#3F8BBAFF", "#3482B6FF", "#2B7BB4FF", "#1F74B1FF", 
                                 "#1C6CAAFF", "#1C63A1FF", "#1C5A99FF", "#21538BFF", "#244C7CFF", 
                                 "#26456EFF"))(n)
  }else if(colorset == "solar_dim"){
    colors <- colorRampPalette(c( 
                                 "#C51517FF", "#CB1618FF", "#D21E1CFF", "#D73529FF", "#DC4636FF", 
                                 "#E25644FF", "#EC6857FF", "#F6796AFF", "#F78E80FF", "#E9A79DFF", 
                                 "#D6BFBBFF", 
                                 "#BBC5CCFF", "#9CBBCFFF", "#78B1D3FF", "#5EA5CEFF", "#4F98C4FF", 
                                 "#3F8BBAFF", "#3482B6FF", "#2B7BB4FF", "#1F74B1FF", "#1C6CAAFF", 
                                 "#1C63A1FF"))(n)
  }else if(colorset == "solar_blood"){
    colors <- colorRampPalette(c("#8B0000", "#FF0000",  "#FF6347",  "#FFCCCB",  "#F0F0F0",  "#CBEAFF",  "#87CEFA",  "#1E90FF",  "#00008B"))(n)
  }else if(colorset == "solar_blood2"){
    colors <- colorRampPalette(c("#8B0000", "#FF0000",  "#FF6347", "#F0F0F0", "#87CEFA",  "#1E90FF",  "#00008B"))(n)
  }else if(colorset == "solar_sky"){
    colors <- colorRampPalette(c("#8B0000", "#D55E5E",  "#F08080",  "#FFB6C1",  "#F5F5F5",  "#B3E0FF",  "#6495ED",  "#4169E1", "#00008B"))(n)
  }else if(colorset == "earth"){
    colors <- colorRampPalette(c("#417839FF", "#9BB655FF",  "#D3D5AEFF",  "#EDD6D5FF",  "#D89873FF",  "#803342FF"))(n)
  }else if(colorset == "laser"){
    colors <- colorRampPalette(c("#FCF5F2FF", "#FFEEE7FF", "#FFE4DBFF", "#FFD8D0FF", "#FFCAC6FF", 
                                 "#FFBCBDFF", "#FEADB6FF", "#FB9DB1FF", "#F68CADFF", "#F079AAFF", 
                                 "#EA66A9FF", "#E24FA8FF", "#DA32A9FF", "#C527A2FF", "#AF1F9AFF", 
                                 "#991690FF", "#840D84FF", "#6E0378FF", "#5A006CFF", "#490062FF"
    ))(n)
  }else if(colorset == "iridescent"){
    colors <- colorRampPalette(c("#F5F2D8FF", "#EBF0CCFF", "#DBEDC3FF", "#C9E9BDFF", "#B4E4B9FF", 
                                 "#9EDEB7FF", "#86D8B8FF", "#6DD0B9FF", "#55C7BCFF", "#3DBDBFFF", 
                                 "#2DB2C1FF", "#2DA6C2FF", "#3B98C2FF", "#4E89C0FF", "#6078BBFF", 
                                 "#7066B3FF", "#7C52A9FF", "#833E99FF", "#842A84FF", "#80146EFF"
    ))(n)
  }else{
    colors <- colorRampPalette(RColorBrewer::brewer.pal(11, name = colorset))(n)
  }
  return(colors)
}


     








Seurat_feature_plot <- function(seu, file = NULL, width = 4, height = 5, pt.size = 0.1, weight = NULL, slot = "umap", features = NULL, colorset = "greens", rasterize = FALSE, binarize = FALSE, split_by = NULL, ...){
  
  colors <- choose_colorset(colorset)
  
  if(!is.null(file)){
    pdf(file, width = width, height = height)
  }
  
  for(gene in features){
    # for rna space
    data <- data.frame(seu@reductions[[slot]]@cell.embeddings)
    colnames(data) <- c("x", "y")
    data <- cbind(data, seu@meta.data)
    x_ratio <- max(data$x) - min(data$x)
    y_ratio <- max(data$y) - min(data$y)
    # # for atac space
    # data <- proj2@embeddings$UMAP$df
    # colnames(data) <- c("x", "y")
    # rownames(data) <- substr(rownames(data), 14, 36)
    # data <- data[seu$barcode,]
    
    
    if(!gene %in% rownames(seu@assays$RNA@data)){
      message(paste0(gene, " not detected. Skip ", gene, "."))
      next
    }else{
      message(paste0("Plotting ", gene, "..."))
    }
    
    if(!is.null(weight)){
      weightList <- weight$Weights
      message("Using smoothed gene expression...")
      data$color <- calculated_smoothed_ge(proj = NULL, weight = weightList, mat = seu@assays$RNA@data[gene,,drop=FALSE], threads = 16)[rownames(data)]
    }else{
      data$color <- seu@assays$RNA@data[gene,]
    }
    
    if(binarize){
      data$color  <- ifelse(data$color > 0, "positive", "negative")
    }
    p <- ggplot(data, aes(x = x, y = y))
    if(rasterize){
      p <- p + geom_scattermore(aes(colour = color), pointsize = pt.size) + theme_classic()
    }else{
      p <- p + geom_point(aes(color = color), size = pt.size) + theme_classic()
    }
    if(binarize){
      p <- p + scale_color_manual(values = c("#deddda", "#c01c28"))
    }else{
      p <- p + scale_colour_gradientn(colors = colors, guide = guide_colourbar(direction = "horizontal", title.position = "top"))
    }
    p <- p +
      xlab(paste0(toupper(slot), " 1")) +
      ylab(paste0(toupper(slot), " 2")) +
      labs(color = "Counts") +
      ggtitle(gene) +
      theme(legend.position = c(0.5, -0.3),
            plot.margin = margin(c(5,5,80,5), unit = "pt"),
            plot.title = element_text(hjust = 0.5)) +
      coord_fixed(ratio = x_ratio/y_ratio)
    if(!is.null(split_by)){p <- p + facet_wrap(split_by)}
    if(!is.null(file)){
      print(p)
    }else{
      print(p)
      return(p)
    }
  }
  if(!is.null(file)){
    dev.off()
  }
}

FeaturePlotOb <- function(seu, file = NULL, width = 4, height = 5, pt.size = 0.1, weight = NULL, slot = "umap", assay = "SCT", features = NULL, colorset = "greens", rasterize = FALSE, binarize = FALSE, split_by = NULL, alpha = 1, ...){
  
  colors <- choose_colorset(colorset)
  
  if(!is.null(file)){
    pdf(file, width = width, height = height)
  }
  
  for(gene in features){
    # for rna space
    data <- data.frame(seu@reductions[[slot]]@cell.embeddings)
    colnames(data) <- c("x", "y")
    data <- cbind(data, seu@meta.data)
    x_ratio <- max(data$x) - min(data$x)
    y_ratio <- max(data$y) - min(data$y)
    # # for atac space
    # data <- proj2@embeddings$UMAP$df
    # colnames(data) <- c("x", "y")
    # rownames(data) <- substr(rownames(data), 14, 36)
    # data <- data[seu$barcode,]
    
    
    if(!gene %in% rownames(seu@assays$RNA@data)){
      message(paste0(gene, " not detected. Skip ", gene, "."))
      next
    }else{
      message(paste0("Plotting ", gene, "..."))
    }
    
    if(!is.null(weight)){
      weightList <- weight$Weights
      message("Using smoothed gene expression...")
      data$color <- calculated_smoothed_ge(proj = NULL, weight = weightList, mat = seu@assays$RNA@data[gene,,drop=FALSE], threads = 16)[rownames(data)]
    }else{
      if(!gene %in% rownames(seu@assays[[assay]]@data)){
        message(paste0(gene, " not detected. Skip ", gene, "."))
        next
      }else{
        message(paste0("Plotting ", gene, "..."))
      }
      data$color <- seu@assays[[assay]]@data[gene,]
    }
    data_bg <- data[data$color == 0,]
    data_fw <- data[data$color > 0,]
    data <- rbind(data_bg, data_fw)
    
    if(binarize){
      data$color  <- ifelse(data$color > 0, "positive", "negative")
    }
    p <- ggplot(data, aes(x = x, y = y))
    if(rasterize){
      p <- p + geom_scattermore(aes(colour = color), pointsize = pt.size) + theme_classic()
    }else{
      p <- p + geom_point(aes(color = color), size = pt.size, shape = 16, alpha = alpha) + theme_classic()
    }
    if(binarize){
      p <- p + scale_color_manual(values = c("#deddda", "#c01c28"))
    }else{
      p <- p + scale_colour_gradientn(colors = colors, guide = guide_colourbar(direction = "horizontal", title.position = "top"))
    }
    p <- p +
      xlab(paste0(toupper(slot), " 1")) +
      ylab(paste0(toupper(slot), " 2")) +
      labs(color = "Normalized counts") +
      ggtitle(gene) +
      theme(legend.position = c(0.5, -0.3),
            plot.margin = margin(c(5,5,80,5), unit = "pt"),
            plot.title = element_text(hjust = 0.5)) +
      coord_fixed(ratio = x_ratio/y_ratio)
    if(!is.null(split_by)){p <- p + facet_wrap(split_by)}
    if(!is.null(file)){
      print(p)
    }else{
      print(p)
      return(p)
    }
  }
  if(!is.null(file)){
    dev.off()
  }
}



calculated_smoothed_ge <- function (weightList = NULL, mat = NULL, threads = 16, hdf5 = TRUE){
  if(hdf5 == FALSE){
    smoothed_expression <- lapply(seq_along(weightList), function(x) {
      matx <- parallel::mclapply(seq_along(weightList[[x]]), function(y) {
        Matrix::t(as.matrix(weightList[[x]][[y]]) %*% Matrix::t(mat[, paste0(colnames(weightList[[x]][[y]])), drop = FALSE]))
      }, mc.cores = threads) %>% Reduce("cbind", .)
      return(matx[, colnames(mat)])
    }) %>% Reduce("+", .)
    smoothed_expression <- smoothed_expression/length(weightList)
    
  }else{
    smoothed_expression <- lapply(seq_along(weightList), function(x) {
      h5df <- h5ls(weightList[[x]])
      blocks <- gtools::mixedsort(grep("block", h5df$name, value = TRUE))
      matx <- mclapply(seq_along(blocks), function(y) {
        
        bn <- h5read(weightList[[x]], paste0(blocks[y], "/Names"))
        by <- h5read(weightList[[x]], paste0(blocks[y], "/Weights"))
        colnames(by) <- bn
        rownames(by) <- bn
        Matrix::t(by %*% Matrix::t(mat[, paste0(bn), drop = FALSE]))
      }, mc.cores = threads) %>% Reduce("cbind", .)
      matx[, colnames(mat)]
    }) %>% Reduce("+", .)
    smoothed_expression <- smoothed_expression/length(weightList)
  }

  return(smoothed_expression)  

}


ArchR_feature_plot <- function(proj, file = NULL, dimreduc = NULL, weight = NULL, width = 4, height = 5, pt.size = 0.1, matrix_to_use = "GeneScoreMatrix", features = NULL, colorset = "solar_extra", use_assay = "z"){
  
  gmat <- getMatrixFromProject(proj, useMatrix = matrix_to_use)
  
  colors <- choose_colorset(colorset)
  
 
  
  if(!is.null(file)){
    pdf(file, width = width, height = height)
  }
  for(gene in features){
    # ATAC space
    if(is.null(dimreduc)){
      data <- proj@embeddings$UMAP$df
    }else{
      data <- as.data.frame(dimreduc)
    }
    colnames(data) <- c("x", "y")
    x_ratio <- max(data$x) - min(data$x)
    y_ratio <- max(data$y) - min(data$y)
    
    if(is.null(weight)){
      weightList <- proj@imputeWeights$Weights
    }else{
      weightList <- weight$Weights
    }
    
    # # RNA space
    # if(!gene %in% rowData(gmat)$name){next}
    # data <- data.frame(seu@reductions$umap@cell.embeddings)
    # rownames(data) <- substr(rownames(data), 1, 23)
    # data <- data[substr(proj$cellNames, 14, 36),]
    # colnames(data) <- c("x", "y")
    
    message(paste0("Ploting ", gene, "..."))
    if(matrix_to_use == "MotifMatrix"){
      mat <- assays(gmat)[[use_assay]][which(rowData(gmat)$name == gene),,drop = FALSE]
    }else{
      mat <- assay(gmat[which(rowData(gmat)$name == gene),])
    }
    
    if(nrow(mat) == 0){
      message(paste0(gene, " not detected. Skip ", gene, "."))
      next
    }

    data$color <- calculated_smoothed_ge(proj = NULL, weight = weightList, mat = mat, threads = 16)[rownames(data)]
    
    p <- ggplot(data, aes(x = x, y = y)) +
      geom_point(aes(color = color), size = pt.size) +
      theme_classic() +
      scale_color_gradientn(colors = colors, guide = guide_colourbar(direction = "horizontal", title.position = "top")) +
      # guides(guide = guide_colourbar(direction = "horizontal", title.position = "top")) +
      # scale_color_viridis(guide = guide_colourbar(direction = "horizontal", title.position = "top")) +
      xlab("UMAP 1") +
      ylab("UMAP 2") +
      labs(color = "Log2(NormCounts +1)") +
      ggtitle(gene) +
      theme(legend.position = c(0.5, -0.3),
            plot.margin = margin(c(5,5,80,5), unit = "pt"),
            plot.title = element_text(hjust = 0.5)) +
      coord_fixed(ratio = x_ratio/y_ratio)
    if(!is.null(file)){
      print(p)
    }else{
      print(p)
      return(p)
    }
  }
  if(!is.null(file)){
    dev.off()
  }
}

DORC_feature_plot <- function(dorcmat, file = NULL, dimreduc = NULL, weight = NULL, width = 4, height = 5, pt.size = 0.1, features = NULL, colorset = "solar_extra", limits = NULL){
  dorcmat <- as.matrix(dorcmat)
  
  colors <- choose_colorset(colorset)
  
  if(!is.null(file)){
    pdf(file, width = width, height = height)
  }
  
  if(is.null(dimreduc)){
    stop("Dim reduction not provided!")
  }else{
    data <- as.data.frame(dimreduc)
  }
  
  for(gene in features){
    colnames(data) <- c("x", "y")
    x_ratio <- max(data$x) - min(data$x)
    y_ratio <- max(data$y) - min(data$y)
    
    if(!is.null(weight)){
      weightList <- weight$Weights
    }
    
    message(paste0("Ploting ", gene, "..."))
    
    
    if(!gene %in% rownames(dorcmat)){
      message(paste0(gene, " not detected. Skip ", gene, "."))
      next
    }else{
      mat <- dorcmat[gene,rownames(data),drop = FALSE]
    }
    
    if(!is.null(weight)){
      data$color <- calculated_smoothed_ge(weight = weightList, mat = mat, threads = 16)[rownames(data)]
    }else{
      data$color <- mat[gene,]
    }
    
    p <- ggplot(data, aes(x = x, y = y)) +
      geom_point(aes(color = color), size = pt.size) +
      theme_classic()
    if(is.null(limits)){
      p <- p + scale_color_gradientn(colors = colors, guide = guide_colourbar(direction = "horizontal", title.position = "top"))
    }else{
      p <-p + scale_color_gradientn(colors = colors, guide = guide_colourbar(direction = "horizontal", title.position = "top"), limits = limits, oob = scales::squish)
    }
      # guides(guide = guide_colourbar(direction = "horizontal", title.position = "top")) +
      # scale_color_viridis(guide = guide_colourbar(direction = "horizontal", title.position = "top")) +
    p <- p + xlab("UMAP 1") +
      ylab("UMAP 2") +
      labs(color = "Counts") +
      ggtitle(gene) +
      theme(legend.position = c(0.5, -0.3),
            plot.margin = margin(c(5,5,80,5), unit = "pt"),
            plot.title = element_text(hjust = 0.5)) +
      coord_fixed(ratio = x_ratio/y_ratio)
    if(!is.null(file)){
      print(p)
    }else{
      print(p)
      return(p)
    }
  }
  if(!is.null(file)){
    dev.off()
  }
}


DimPlotArchR <- function(proj, colors = NULL, group.by = NULL, pt.size = 0.1, slot = "UMAP"){
  data <- proj@embeddings[[slot]]$df
  colnames(data) <- c("x", "y")
  data[[group.by]] <- as.vector(proj@cellColData[[group.by]])
  x_ratio <- max(data$x) - min(data$x)
  y_ratio <- max(data$y) - min(data$y)
  
  p <- ggplot(data, aes(x = x, y = y)) +
    geom_point(aes_string(color = group.by), size = pt.size) +
    theme_classic() +
    scale_color_manual(values = colors) +
    # guides(guide = guide_colourbar(direction = "horizontal", title.position = "top")) +
    # scale_color_viridis(guide = guide_colourbar(direction = "horizontal", title.position = "top")) +
    xlab(paste0(slot, "_1")) +
    xlab(paste0(slot, "_2")) +
    theme(axis.text.x = element_text(size = 10, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.title.x = element_text(size = 13),
          axis.title.y = element_text(size = 13)) +
    guides(color = guide_legend(override.aes = list(size = 2))) +
    coord_fixed(ratio = x_ratio/y_ratio)
  print(p)
}

ArchR_to_Seurat <- function(proj, samples = NULL, fragment_files = NULL, pkm, annotation){
  seurat_list <- lapply(samples, function(cur_sample) {
    cur_fragments <- fragment_files[[which(samples == cur_sample)]]
    cur_pm <- pkm[, grepl(paste0(cur_sample, "#"), colnames(pkm))]
    
    cur_meta <- proj@cellColData %>% as.data.frame %>% subset(Sample == cur_sample)
    colnames(cur_pm) <- do.call(rbind, str_split(colnames(cur_pm), "#"))[, 2]
    rownames(cur_meta) <- do.call(rbind, str_split(rownames(cur_meta), "#"))[, 2]
    
    print(dim(cur_pm))
    if(dim(cur_pm)[2] == 0){return(NA)}
    cur_chromatin <- Signac::CreateChromatinAssay(counts = cur_pm, sep = c("-", "-"), fragments = cur_fragments, ranges = proj@peakSet, genome = unique(genome(annotation)), annotation = annotation)
    cur_atac <- Seurat::CreateSeuratObject(cur_chromatin, assay = "ATAC", meta.data = cur_meta)
  })
  # seurat_list <- seurat_list[-which(is.na(seurat_list))]
  print(seurat_list)
  seu_atac <- merge(x = seurat_list[[1]], y = seurat_list[2:length(seurat_list)])
  seu_atac <- RenameCells(seu_atac, new.names = paste0(seu_atac$Sample, "#", seu_atac$barcode))
  return(seu_atac)
}

seurat_reduction_from_ArchR <- function(seu, proj, reduc_name = "IterativeLSI"){
  if(reduc_name == "IterativeLSI"){
    LSI_matrix <- proj@reducedDims[[reduc_name]]$matSVD
    LSI_matrix <- LSI_matrix[colnames(seu),]
    rownames(LSI_matrix) <- colnames(seu)
    colnames(LSI_matrix) <- paste0("LSI_", 1:ncol(LSI_matrix))
    out <- Seurat::CreateDimReducObject(embeddings = LSI_matrix, assay = "ATAC")
  }else if(reduc_name == "LSI_ATAC"){
    LSI_matrix <- proj@reducedDims[[reduc_name]]$matSVD
    LSI_matrix <- LSI_matrix[colnames(seu),]
    rownames(LSI_matrix) <- colnames(seu)
    colnames(LSI_matrix) <- paste0("LSI_", 1:ncol(LSI_matrix))
    out <- Seurat::CreateDimReducObject(embeddings = LSI_matrix, assay = "ATAC")    
  }else{
    LSI_matrix <- proj@reducedDims[[reduc_name]]$matDR
    LSI_matrix <- LSI_matrix[colnames(seu),]
    rownames(LSI_matrix) <- colnames(seu)
    colnames(LSI_matrix) <- paste0("LSI_", 1:ncol(LSI_matrix))
    out <- Seurat::CreateDimReducObject(embeddings = LSI_matrix, assay = "ATAC")    
  }
  return(out)
}




add_magic_dimreduc <- function(matDR = NULL, dimsToUse = NULL, td = 3, ka = 4, sampleCells = 5000, nRep = 2, k = 15, epsilon = 1, threads = 24, seed = 1){
  set.seed(seed)
  if(is.null(rownames(matDR))){stop("Must provide rownames of matDR!")}
  N <- nrow(matDR)
  rn <- rownames(matDR)
  if (!is.null(sampleCells)){
    if (sampleCells > nrow(matDR)) {sampleCells <- NULL}
  }
  if (is.null(sampleCells)){
    binSize <- N
    nRep <- 1
  }else{
    cutoffs <- lapply(seq_len(1000), function(x) {
      N/x
    }) %>% unlist
    binSize <- min(cutoffs[order(abs(cutoffs - sampleCells))[1]] + 1, N)
  }
  
  weightList <- mclapply(seq_len(nRep), function(y) {
    if (!is.null(sampleCells)) {
      idx <- sample(seq_len(nrow(matDR)), nrow(matDR))
      blocks <- split(rownames(matDR)[idx], ceiling(seq_along(idx)/binSize))
    }else {
      blocks <- list(rownames(matDR))
    }
    
    blockList <- lapply(seq_along(blocks), function(x) {
      ix <- blocks[[x]]
      Nx <- length(ix)
      knnObj <- nabor::knn(data = matDR[ix,], query = matDR[ix,], k = k)
      knnIdx <- knnObj$nn.idx
      knnDist <- knnObj$nn.dists
      rm(knnObj)
      if (ka > 0) {
        knnDist <- knnDist/knnDist[, ka]
      }
      if (epsilon > 0) {
        W <- Matrix::sparseMatrix(rep(seq_len(Nx), k), c(knnIdx), x = c(knnDist), dims = c(Nx, Nx))
      }
      else {
        W <- Matrix::sparseMatrix(rep(seq_len(Nx), k), c(knnIdx), x = 1, dims = c(Nx, Nx))
      }
      W <- W + Matrix::t(W)
      if (epsilon > 0) {
        W@x <- exp(-(W@x/epsilon^2))
      }
      W <- W/Matrix::rowSums(W)
      Wt <- W
      for (i in seq_len(td)) {
        Wt <- Wt %*% W
      }
      rownames(Wt) <- ix
      colnames(Wt) <- ix
      rm(knnIdx)
      rm(knnDist)
      rm(W)
      gc()
      return(Wt)
    }) %>% SimpleList
    names(blockList) <- paste0("b", seq_along(blockList))
    return(blockList)
  }, mc.cores = threads) %>% SimpleList
  
  names(weightList) <- paste0("w", seq_along(weightList))
  SimpleList(Weights = weightList, Params = list(td = td, k = k, ka = ka, epsilon = epsilon))
}

add_magic_knn <- function(knn = NULL, td = 3, ka = 4, epsilon = 1, threads = 24, nRep = 1){
  k <- ncol(knn@nn.idx)
  if(is.null(knn@cell.names)){stop("Must provide cellnames in knn@cell.names!")}
  weightList <- mclapply(seq_len(nRep), function(y) {
    blockList <- lapply(seq_along(blocks), function(x) {
      Nx <- nrow(knn@nn.idx)
      knnIdx <- knn@nn.idx
      knnDist <- knn@nn.dist
      if (ka > 0) {
        knnDist <- knnDist/knnDist[, ka]
      }
      if (epsilon > 0) {
        W <- Matrix::sparseMatrix(rep(seq_len(Nx), k), c(knnIdx), x = c(knnDist), dims = c(Nx, Nx))
      }
      else {
        W <- Matrix::sparseMatrix(rep(seq_len(Nx), k), c(knnIdx), x = 1, dims = c(Nx, Nx))
      }
      W <- W + Matrix::t(W)
      if (epsilon > 0) {
        W@x <- exp(-(W@x/epsilon^2))
      }
      W <- W/Matrix::rowSums(W)
      Wt <- W
      for (i in seq_len(td)) {
        Wt <- Wt %*% W
      }
      rownames(Wt) <- knn@cell.names
      colnames(Wt) <- knn@cell.names
      rm(knnIdx)
      rm(knnDist)
      rm(W)
      gc()
      return(Wt)
    }) %>% SimpleList
    names(blockList) <- paste0("b", seq_along(blockList))
    return(blockList)
  }, mc.cores = threads) %>% SimpleList
  
  names(weightList) <- paste0("w", seq_along(weightList))
  SimpleList(Weights = weightList, Params = list(td = td, k = k, ka = ka, epsilon = epsilon))
}


smooth_mat_magic <- function (mat = NULL, weight = NULL, threads = 16){
  weightList <- weight
  smoothed_expression <- lapply(seq_along(weightList), function(x) {
    matx <- parallel::mclapply(seq_along(weightList[[x]]), function(y) {
      Matrix::t(as.matrix(weightList[[x]][[y]]) %*% Matrix::t(mat[, paste0(colnames(weightList[[x]][[y]])), drop = FALSE]))
    }, mc.cores = threads) %>% Reduce("cbind", .)
    return(matx[, colnames(mat)])
  }) %>% Reduce("+", .)
  smoothed_expression <- smoothed_expression/length(weightList)
  return(smoothed_expression)
}

smooth_mat_knn <- function (NNmat, mat, geneList = NULL, barcodesList = NULL, nCores = 16) {
  stopifnot(all.equal(nrow(NNmat), ncol(mat)))
  if (is.null(rownames(NNmat))) 
    stop("NN matrix has to have matching cell IDs as rownames\n")
  if (!all.equal(rownames(NNmat), colnames(mat))) 
    stop("Nearest-neighbor matrix and cell data matrix don't have matching cells barcodes ..\n")
  cat("Number of cells in supplied matrix: ", ncol(mat), "\n")
  cat("Number of genes in supplied matrix: ", nrow(mat), "\n")
  cat("Number of nearest neighbors being used per cell for smoothing: ", 
      ncol(NNmat), "\n")
  if (!is.null(geneList)) {
    if (!(all(geneList %in% rownames(mat)))) {
      cat("One or more of the gene names supplied is not present in the matrix provided: \n")
      cat(geneList[!geneList %in% rownames(mat)], sep = ", ")
      cat("\n")
      stop()
    }
    cat("Running smoothing only on genes:", geneList, sep = "\n")
    cat("........\n")
    mat <- mat[rownames(mat) %in% geneList, ]
  }
  else {
    if (nrow(mat) > 10000) {
      cat("Running smoothing for all genes in matrix! (n = ", 
          nrow(mat), ") This is bound to take more time than querying specific markers ..\n", 
          sep = "")
    }
  }

  if (!is.null(barcodesList)) {
    cat("Subsetting to ", length(barcodesList), " barcodes in dataset..\n")
    NNmat <- NNmat[barcodesList, ]
  }
  
  cat("Running in parallel using ", nCores, "cores ..\n")
  time_elapsed <- Sys.time()
  
  smoothedMat <- parallel::mclapply(1:nrow(NNmat), function(x){
    smoothedScore <- data.table(Matrix::rowMeans(mat[, NNmat[x, ], drop = FALSE]))
    rownames(smoothedScore) <- rownames(mat)
    colnames(smoothedScore) <- rownames(NNmat)[x]
    return(smoothedScore)
  }, mc.cores = nCores) %>% dplyr::bind_cols() %>% data.matrix() %>% Matrix::Matrix(sparse = TRUE)

  rownames(smoothedMat) <- rownames(mat)
  time_elapsed <- Sys.time() - time_elapsed
  cat(paste("\nTime Elapsed: ", time_elapsed, units(time_elapsed), "\n"))
  return(smoothedMat)
}







pseudotime_smoothed_exp <- function(expmat, pseudotime, span = 0.3, n = 100, normalize = TRUE, cores = 16){
  pred <- mclapply(rownames(expmat), function(x){
    fit <- loess(score ~ pseudotime, data = data.frame(score = expmat[x,], pseudotime = pseudotime), span = span)
    out <- predict(fit, data.frame(pseudotime = seq(0.5, 99.5, length.out = n)))
    # print(out)
    # out <- (out-min(out))/(max(out)-min(out))
    return(out)
  }, mc.cores = cores, mc.preschedule = TRUE) %>% Reduce("rbind", .)
  if(normalize){
    cat("min-max normalize smoothed matrix...")
    pred_out <- t(apply(pred, 1, function(x) (x-min(x))/(max(x)-min(x))))
  }else{
    cat("Output without normalization...")
    pred_out <- pred
  }
  rownames(pred_out) <- rownames(expmat)
  return(pred_out)
}
  

binarysort_heatmap <- function (m = NULL, scale = TRUE, cutOff = 1, clusterCols = TRUE, invert = TRUE) 
{
  if(invert){
    m <- -m
  }
  m <- as.matrix(m)
  if(scale){
    lmat <- sweep(m - rowMeans(m), 1, matrixStats::rowSds(m), `/`)
  }else{
    lmat <- m
  }
  lmat <- lmat >= cutOff
  if (clusterCols) {
    m <- t(m)
    lmat <- t(lmat)
    hc <- hclust(dist(m))
    colIdx <- hc$order
    m <- t(m[colIdx, ])
    lmat <- t(lmat[colIdx, ])
  }
  else {
    hc <- NULL
  }
  rowIdx <- do.call("order", c(as.data.frame(lmat)[seq_len(ncol(lmat))], list(decreasing = TRUE)))
  m <- m[rowIdx, ]
  if(invert){
    return(list(mat = -m, hclust = hc))
  }else{
    return(list(mat = m, hclust = hc))
  }
}









extract_TF_names <- function(motifIDs){
  if (all(grepl("_", motifIDs, fixed = TRUE))) {
    sapply(strsplit(sapply(strsplit(motifIDs, "_LINE.", fixed = FALSE), "[[", 2), "_", fixed = FALSE), "[[", 2)
  }else{
    motifIDs
  }
}

analyze_FigR_GRN <- function(ATAC.se, 
                             dorcK = 30, # peaks that are related to nearby k genes for each dorc are also calculated for motif enrichment analysis and TF-dorc correlation
                             dorcTab, # Filtered gene-peak correlation dataframe. With pvalZ > 0.05 and rObs > 0
                             n_bg = 50, # BG peak number for motif enrichment analysis
                             genome, # hg38, mm10
                             pwm = NULL, # motif PWM matrix
                             dorcMat, # Smoothed dorc * cell expression matrix. Mean normalized.
                             rnaMat, # Smoothed rna * cell expression matrix. Must cover all genes.
                             dorcGenes = NULL, # dorcGenes to use
                             nCores = 16,
                             seed = 123) 
{
  stopifnot(all.equal(ncol(dorcMat), ncol(rnaMat)))
  if (!all(c("Peak", "Gene") %in% colnames(dorcTab))) 
    stop("Expecting fields Peak and Gene in dorcTab data.frame .. see runGenePeakcorr function in BuenRTools")
  if (all(grepl("chr", dorcTab$Peak, ignore.case = TRUE))) {
    usePeakNames <- TRUE
    message("Detected peak region names in Peak field")
    if (!(all(grepl("chr", rownames(ATAC.se), ignore.case = TRUE)))) 
      stop("Peak regions provided in dorcTab data.frame but not found as rownames in input SE")
    if (!all(dorcTab$Peak %in% rownames(ATAC.se))) 
      stop("Found DORC peak region not present in input SE.. make sure DORC calling output corresponds to same input SE as the one provided here ..")
  }else{
    usePeakNames <- FALSE
    message("Assuming peak indices in Peak field")
    if (max(dorcTab$Peak) > nrow(ATAC.se)) 
      stop("Found DORC peak index outside range of input SE.. make sure DORC calling output corresponds to same input SE as the one provided here ..")
  }
  if (is.null(dorcGenes)) {
    dorcGenes <- rownames(dorcMat)
  }
  else {
    cat("Using specified list of dorc genes ..\n")
    if (!(all(dorcGenes %in% rownames(dorcMat)))) {
      cat("One or more of the gene names supplied is not present in the DORC matrix provided: \n")
      cat(dorcGenes[!dorcGenes %in% rownames(dorcMat)], sep = ", ")
      cat("\n")
      stop()
    }
  }
  DORC.knn <- FNN::get.knn(data = t(scale(Matrix::t(dorcMat))), k = dorcK)$nn.index
  rownames(DORC.knn) <- rownames(dorcMat)
  if (is.null(SummarizedExperiment::rowData(ATAC.se)$bias)) {
    if (genome %in% "hg19") 
      myGenome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    if (genome %in% "mm10") 
      myGenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
    if (genome %in% "hg38") 
      myGenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    ATAC.se <- chromVAR::addGCBias(ATAC.se, genome = myGenome)
  }
  packagePath <- find.package("FigR", lib.loc = NULL, quiet = TRUE)
  if (grepl("hg", genome)) {
    pwm <- readRDS(paste0(packagePath, "/data/cisBP_human_pfms_2021.rds"))
  }
  else {
    pwm <- readRDS(paste0(packagePath, "/data/cisBP_mouse_pfms_2021.rds"))
  }
  if (all(grepl("_", names(pwm), fixed = TRUE))) 
    names(pwm) <- extract_TF_names(names(pwm))
  message("Removing genes with 0 expression across cells ..\n")
  rnaMat <- rnaMat[Matrix::rowSums(rnaMat) != 0, ]
  myGeneNames <- gsub(x = rownames(rnaMat), pattern = "-", replacement = "")
  rownames(rnaMat) <- myGeneNames
  motifsToKeep <- intersect(names(pwm), myGeneNames)
  cat("Getting peak x motif matches ..\n")
  motif_ix <- motifmatchr::matchMotifs(subject = ATAC.se, pwms = pwm[motifsToKeep], genome = genome)
  motif_ix <- motif_ix[, Matrix::colSums(assay(motif_ix)) != 0]
  cat("Determining background peaks ..\n")
  cat("Using ", n_bg, " iterations ..\n\n")
  if (any(Matrix::rowSums(assay(ATAC.se)) == 0)) {
    ATAC.mat <- assay(ATAC.se)
    ATAC.mat <- cbind(ATAC.mat, 1)
    ATAC.se.new <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = ATAC.mat), rowRanges = granges(ATAC.se))
    set.seed(seed)
    bg <- chromVAR::getBackgroundPeaks(ATAC.se.new, niterations = n_bg)
  }
  else {
    set.seed(seed)
    bg <- chromVAR::getBackgroundPeaks(ATAC.se, niterations = n_bg)
  }
  cat("Testing ", length(motifsToKeep), " TFs\n")
  cat("Testing ", nrow(dorcMat), " DORCs\n")
  
  TFenrich.d <- mclapply(dorcGenes, function(x){
    # dorc peak motif enrichment
    DORCNNpeaks <- unique(dorcTab$Peak[dorcTab$Gene %in% c(x, rownames(dorcMat)[DORC.knn[x, ]])]) # peak-gene correlation of dorc and its knn nearby genes
    if (usePeakNames) DORCNNpeaks <- which(rownames(ATAC.se) %in% DORCNNpeaks)
    mZ <- motif_peak_ztest(peakSet = DORCNNpeaks, bgPeaks = bg, tfMat = assay(motif_ix))
    mZ <- mZ[, c("gene", "z_test")]
    colnames(mZ) <- c("Motif", "Enrichment.Z")
    mZ$Enrichment.P <- 2 * pnorm(abs(mZ$Enrichment.Z), lower.tail = FALSE)
    mZ$Enrichment.log10P <- sign(mZ$Enrichment.Z) * -log10(mZ$Enrichment.P)
    mZ <- cbind(DORC = x, mZ)
    # TF rna exp - dorc correlation
    corr.r <- cor(dorcMat[x, ], t(as.matrix(rnaMat[mZ$Motif, ])), method = "spearman")
    stopifnot(all.equal(colnames(corr.r), mZ$Motif))
    mZ$Corr <- corr.r[1, ]
    mZ$Corr.Z <- scale(mZ$Corr, center = TRUE, scale = TRUE)[,1]
    mZ$Corr.P <- 2 * pnorm(abs(mZ$Corr.Z), lower.tail = FALSE)
    mZ$Corr.log10P <- sign(mZ$Corr.Z) * -log10(mZ$Corr.P)
    return(mZ)
  }, mc.cores = nCores) %>% Reduce("rbind", .)
  cat("Finished!\n")
  TFenrich.d <- data.frame(TFenrich.d)
  rownames(TFenrich.d) <- NULL
  print(str(TFenrich.d))
  if(require(Rmpfr)){
    TFenrich.d <- TFenrich.d %>% dplyr::mutate(Score = sign(Corr) * as.numeric(-log10(1 - (1 - Rmpfr::mpfr(Enrichment.P, 100)) * (1 - Rmpfr::mpfr(Corr.P, 100)))))
  }else{
    message("Package Rmpfr is not installed. The scores calculated may lose precision.\n")
    TFenrich.d <- TFenrich.d %>% dplyr::mutate(Score = sign(Corr) * as.numeric(-log10(1 - (1 - Enrichment.P) * (1 - Corr.P))))
  }
  
  TFenrich.d$Score[TFenrich.d$Enrichment.Z < 0] <- 0
  return(TFenrich.d)
}

extract_TF_names <- function (motifIDs){
  if (all(grepl("_", motifIDs, fixed = TRUE))) {
    sapply(strsplit(sapply(strsplit(motifIDs, "_LINE.", fixed = FALSE), "[[", 2), "_", fixed = FALSE), "[[", 2)
  }else{
    motifIDs
  }
}


motif_peak_ztest <- function (peakSet, bgPeaks, tfMat){
  if (nrow(tfMat) != nrow(bgPeaks)) # tfMat is peak*motif 0-1 matrix
    stop("Reference peak set used for TF and background peaks matrix must match..\n")
  if (!all(peakSet %in% 1:nrow(bgPeaks))) 
    stop("One or more of the provided peak indices are out of the background peak set range ..\n")
  tfMat <- tfMat[, Matrix::colSums(tfMat) != 0]
  cat("Getting selected peak motif frequencies ..\n")
  p.tab <- Matrix::colMeans(tfMat[peakSet, ])
  cat("Getting background peak motif frequencies ..\n")
  bg.f <- as.matrix(bgPeaks[peakSet, ])
  bg.tab <- apply(bg.f[, c(1:ncol(bgPeaks))], 2, function(bg_iter) {
    b.i <- Matrix::colMeans(tfMat[bg_iter, ])
    return(b.i)
  })
  cat("Calculating empirical p values and z score p values ..\n")
  m.p <- dplyr::bind_rows(lapply(names(p.tab), function(motif) {
    s <- sd(bg.tab[motif, ])
    bg_freq <- mean(bg.tab[motif, ])
    z_score <- (p.tab[motif] - bg_freq)/s
    if (is.nan(z_score)) 
      z_score <- 0
    d <- data.frame(motifID = motif, gene = extract_TF_names(motif), 
                    motif_obs_freq = p.tab[motif], motif_bg_freq = mean(bg.tab[motif, 
                    ]), motif_counts = p.tab[motif] * length(peakSet), 
                    emp_pval = 1 - (sum(bg.tab[motif, ] < p.tab[motif])/ncol(bg.tab)), 
                    z_test = z_score, pval.z = 2 * pnorm(-abs(z_score)), 
                    signed.log10p = -log10(2 * pnorm(-abs(z_score))) * 
                      sign(z_score))
    return(d)
  }))
  return(m.p)
}



test_seurat_sparsity <- function(seu, threshold = 0.9, assay = "SCT"){
  mat <- seu@assays[[assay]]@data
  out <- rownames(mat)[Matrix::rowSums(mat > 0)/ncol(mat) > 1-threshold]
  return(out)
}



grn_TF_driver_plot <- function (grn, gene, score.cut = 1, label = TRUE, enrichment_threshold = 0.05, corr_threshold = 0.05, ...) 
{
  if(enrichment_threshold > 0){
    yline <- -log10(enrichment_threshold)
  }else{
    yline <- 0
  }
  if(corr_threshold > 0){
    xline <- -log10(corr_threshold)
  }else{
    xline <- 0
  }  

  if (!gene %in% grn$DORC) 
    stop("Marker specified is not a valid DORC symbol found in the data.frame")
  d <- grn %>% dplyr::filter(DORC %in% gene) %>% mutate(isSig = ifelse(abs(Score) >= score.cut, "Yes", "No"))
  d$Enrichment.P <= enrichment_threshold & d$Corr.P <= corr_threshold
  
  if (label){
    d$Label <- d$Motif
    d$Label[d$isSig %in% "No"] <- ""
  }
  else {
    d$Label <- ""
  }
  
  p <- ggplot(data = d, aes(x = Corr.log10P, y = Enrichment.log10P, color = isSig, label = Label)) + 
    geom_hline(yintercept = yline, color = "gray60", linetype = "dashed") + 
    geom_point(aes(size = isSig)) + 
    theme_classic() + 
    scale_color_manual(values = c("gray66", "firebrick3")) + 
    scale_size_manual(values = c(1, 3)) +
    scale_x_continuous(breaks = scales::pretty_breaks()) + 
    scale_y_continuous(breaks = scales::pretty_breaks()) + 
    labs(y = "Enrichment log10 P", x = "Correlation log10 P", title = gene) + 
    ylim(-ceiling(max(abs(d$Enrichment.log10P))), ceiling(max(abs(d$Enrichment.log10P)))) + 
    xlim(-ceiling(max(abs(d$Corr.log10P))), ceiling(max(abs(d$Corr.log10P)))) + 
    theme(legend.position = "none", 
          axis.text = element_text(color = "black"), 
          plot.title = element_text(hjust = 0.5, face = "italic"), 
          panel.background = element_rect(fill = NA)) + 
    coord_fixed(ceiling(max(abs(d$Corr.log10P))) / ceiling(max(abs(d$Enrichment.log10P)))) +
    geom_text(hjust = -0.2, fontface = "italic", color = "black", size = 5)
  if(corr_threshold > 0){
    p <- p + geom_vline(xintercept = xline, color = "gray60", linetype = "dashed") + geom_vline(xintercept = -xline, color = "gray60", linetype = "dashed")
  }else{
    p <- p + geom_vline(xintercept = xline, color = "gray60", linetype = "dashed")
  }
  
  print(p)
}




grn_rank_driver_plot <- function(grn, rankBy = c("meanScore", "nTargets"), myLabels = NULL, score.cut = NULL, interactive = FALSE, ...) 
{
  if (!rankBy %in% c("meanScore", "nTargets")) 
    stop("rankBy parameter has to be one of meanScore or nTargets to rank drivers using ..\n")
  if (rankBy %in% "meanScore") {
    message("Ranking TFs by mean regulation score across all DORCs ..\n")
    grn.summ <- grn %>% group_by(Motif) %>% 
      dplyr::summarise(Score = mean(Score)) %>% 
      arrange(dplyr::desc(Score)) %>% 
      mutate(Motif = factor(Motif, levels = as.character(Motif)))
    grn.summ$TF <- as.character(grn.summ$Motif)
    if (is.null(myLabels)) {
      grn.summ$TF[grn.summ$Score >= quantile(grn.summ$Score, 0.05) & grn.summ$Score <= quantile(grn.summ$Score, 0.95)] <- ""
    }else{
      grn.summ$TF[!grn.summ$TF %in% myLabels] <- ""
    }
    library(ggrepel)
    p <- ggplot(grn.summ, aes(x = Motif, y = Score, label = TF)) + 
      geom_bar(size = 0.1, stat = "identity", fill = "darkorange", color = NA) + 
      theme_classic() + 
      theme(axis.text.x = element_blank(), axis.text = element_text(color = "black")) + 
      ggrepel::geom_text_repel(size = 3, min.segment.length = 0.1, segment.size = 0.2, max.overlaps = 20) + 
      geom_hline(yintercept = 0) + 
      labs(x = "TF Motifs", y = "Regulation Score") +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1, size = 10), axis.text = element_text(color = "black"))
  }
  else {
    message("Ranking TFs by total number of associated DORCs ..\n")
    if (is.null(score.cut)) {
      message("Regulation score cut-off not specified ... Setting score.cut = 1...\n")
      score.cut <- 1
    }
    message("Using absolute score cut-off of: ", score.cut, " ..\n")
    grn.summ <- grn %>% dplyr::filter(abs(Score) >= score.cut) %>% 
      group_by(Motif) %>% 
      dplyr::select(-DORC) %>% 
      dplyr::summarize(numActivated = sum(Score > 0), numRepressed = sum(Score < 0)) %>% 
      dplyr::mutate(diff = numActivated - numRepressed) %>% 
      mutate(numActivatedY = ifelse(diff >= 0, numActivated, -numActivated), numRepressedY = ifelse(diff <= 0, -numRepressed, numRepressed)) %>% 
      dplyr::arrange(desc(diff)) %>% 
      mutate(Motif = factor(Motif, levels = as.character(Motif))) %>% 
      dplyr::select(-diff) %>% reshape2::melt(id.vars = c("Motif", "numActivated", "numRepressed"))
    p <- ggplot(grn.summ, aes(x = Motif, y = value, fill = variable)) + 
      geom_bar(stat = "identity", color = "lightgray", size = 0.1) + 
      theme_classic() + 
      geom_hline(yintercept = 0) + 
      scale_fill_manual(values = c("firebrick3", "steelblue4"), labels = c("Activated", "Repressed")) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6), axis.text = element_text(color = "black")) + 
      labs(x = "Ranked TF Motifs", y = paste0("# Associated genes \nabs(Score) >= ", score.cut), fill = "Class") + 
      scale_y_continuous(labels = abs) +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1, size = 10), axis.text = element_text(color = "black"))
  }
  if (!interactive) {
    p
  }
  else {
    if (rankBy %in% "meanScore") {
      plotly::ggplotly(p)
    }
    else {
      plotly::ggplotly(p + theme(legend.position = "none", axis.text.x = element_blank()), tooltip = c("Motif", "numActivated", "numRepressed"))
    }
  }
}



grn_heatmap <- function (figR.d, score.cut = 1, DORCs = NULL, TFs = NULL, transpose = TRUE, ...) {
  message("Using absolute score cut-off of: ", score.cut, " ..\n")
  DORCsToKeep <- figR.d %>% dplyr::filter(abs(Score) >= score.cut) %>% pull(DORC) %>% unique()
  TFsToKeep <- figR.d %>% dplyr::filter(abs(Score) >= score.cut) %>% pull(Motif) %>% unique()
  if (!is.null(DORCs)) {
    if (!all(DORCs %in% figR.d$DORC)) 
      stop("One or more DORCs specified is not a valid DORC symbol found in the data.frame")
    DORCsToKeep <- intersect(DORCsToKeep, DORCs)
    TFsToKeep <- figR.d %>% dplyr::filter(abs(Score) >= score.cut & DORC %in% DORCsToKeep) %>% pull(Motif) %>% unique()
  }
  if (!is.null(TFs)) {
    if (!all(TFs %in% figR.d$Motif)) 
      stop("One or more TFs specified is not a valid TF symbol found in the data.frame")
    TFsToKeep <- intersect(TFsToKeep, TFs)
    DORCsToKeep <- figR.d %>% dplyr::filter(abs(Score) >= score.cut & Motif %in% TFsToKeep) %>% pull(DORC) %>% unique()
  }
  net.d <- figR.d %>% dplyr::filter(DORC %in% DORCsToKeep & Motif %in% TFsToKeep) %>% reshape2::dcast(DORC ~ Motif, value.var = "Score") %>% tibble::column_to_rownames("DORC") %>% as.matrix()
  message("Plotting ", nrow(net.d), " DORCs x ", ncol(net.d), "TFs\n")
  if(transpose) net.d <- t(net.d)
  p <- pheatmap::pheatmap(net.d, 
                     border_color = NA,
                     show_colnames = TRUE,
                     show_rownames = TRUE,
                     cluster_rows = TRUE,
                     cluster_cols = TRUE,
                     legend_breaks = c(-2, 0, 2),
                     breaks = seq(-2, 2, length.out = 256),
                     scale = "none",
                     color = choose_colorset("solar_flare", 256),
                     ...)
  p
}




plot_GO <- function(geneid, OrgDb = org.Hs.eg.db, show = 40, ...){
  ontology <- "BP"
  ontology <- "MF"
  ontology <- "CC"
  go_bp <- enrichGO(geneid,
                    OrgDb = OrgDb,
                    ont = "BP",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    minGSSize = 1,
                    pAdjustMethod = "BH")
  go_mf <- enrichGO(geneid,
                    OrgDb = OrgDb,
                    ont = "MF",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    minGSSize = 1,
                    pAdjustMethod = "BH")
  go_cc <- enrichGO(geneid,
                    OrgDb = OrgDb,
                    ont = "CC",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    minGSSize = 1,
                    pAdjustMethod = "BH")
  # dotplot
  data_bp <- go_bp@result
  data_bp$ontology <- "BP"
  data_mf <- go_mf@result
  data_mf$ontology <- "MF"
  data_cc <- go_cc@result
  data_cc$ontology <- "CC"
  data_list <- list(BP = data_bp, MF = data_mf, CC = data_cc)
  
  show <- show
  pl <- list()
  for(data in data_list){
    data$log10p <- -log10(data$p.adjust)
    total_genes <- as.numeric(str_split(data$GeneRatio[1] ,"/")[[1]][2])
    data$Prop <- data$Count/total_genes
    data <- data[order(data$Prop, decreasing = TRUE),]
    p <- ggplot(data = data[1:show,], aes(x = Prop, y = Description)) +
      geom_point(aes(color = p.adjust, size = Count)) +
      scale_color_viridis() +
      scale_y_discrete(limits = rev(data[1:show,2])) +
      theme_bw() +
      theme(panel.grid.major = element_line(),
            panel.border = element_rect(),
            axis.text.y = element_text(size = 13, angle = 0, color = "black"),
            axis.text.x = element_text(size = 13, angle = 0, color = "black"),
            axis.title.y = element_text(size = 13),
            axis.title.x = element_text(size = 13)) +
      # coord_fixed(ratio = max(data$Prop)/show *1.5) +
      ylab("") +
      xlab("Gene ratio") +
      ggtitle(paste0("Top ", show, " of GO ", unique(data$ontology), " terms Enrichment"))
    pl <- c(pl, list(p))
  }
  return(pl)
}


get_corrs <- function(dmat, rmat, ncore){
  dlist <- t(scale(t(dmat))) %>% split(rownames(dmat))
  rlist <- t(scale(t(rmat))) %>% split(rownames(rmat))
  d <- mcmapply(cor, dlist, rlist, mc.cores = ncore)
  return(d)
}

filter_matrix <- function(dmat, rmat, scale_max_abs_val = 10, gene_corr_cutoff = 0, scaling = 'minmax', ncore = 12){
  
  if(scaling == 'minmax'){
    dmat_norm <- t(apply(dmat, 1, function(x) (x-min(x))/(max(x)-min(x))))
    rmat_norm <- t(apply(rmat, 1, function(x) (x-min(x))/(max(x)-min(x))))
  }else if(scaling == 'standard'){
    dmat_norm <- t(scale(t(dmat)))
    rmat_norm <- t(scale(t(rmat)))
    dmat_norm[dmat_norm > scale_max_abs_val] = scale_max_abs_val
    dmat_norm[dmat_norm < -scale_max_abs_val] = -scale_max_abs_val
    rmat_norm[rmat_norm > scale_max_abs_val] = scale_max_abs_val
    rmat_norm[rmat_norm < -scale_max_abs_val] = -scale_max_abs_val
  }
  
  corrs_genes <- get_corrs(dmat = dmat, rmat = rmat, ncore = ncore)
  idx <- which(corrs_genes > gene_corr_cutoff)
  message(paste0(round(sum(corrs_genes <= gene_corr_cutoff)/length(corrs_genes)*100, 2), "% genes filtered out (", sum(corrs_genes <= gene_corr_cutoff), " out of ", length(corrs_genes), " genes)."))
  
  # # For each cell, remove gene corr <= gene_corr_cutoff, and calculate pearson R between yp and yo for this given cell. Store pearson R of this cell in corrs_rows.
  # idx <- which(corrs > gene_corr_cutoff)
  # corrs_cells <- get_corrs(t(dmat_norm[idx,]), t(rmat_norm[idx,]), ncore = ncore)
  
  
  return(list(ndmat = as.matrix(dmat_norm[idx,]), nrmat = as.matrix(rmat_norm[idx,])))
}

calc_velocity <- function(ndmat, nrmat, embedding, max_per_cell = 10, metric = 'cosine', ncore = 12, cutoff = 0.95){
  # return arrow coordinates
  nn <- Seurat::FindNeighbors(t(ndmat), query = t(nrmat), k.param = max_per_cell, annoy.metric = metric, return.neighbor = TRUE)
  dists <- nn@nn.dist
  neighs <- nn@nn.idx
  cutoff <- cutoff
  
  # for each cell. calculate the coordinates of its neighbors (most compatible neighbors)
  umap_new_unsmooth <- parallel::mclapply(1:nrow(embedding), function(x){
    umap_unsmooth <- Matrix::colMeans(embedding[neighs[x,],])
    return(umap_unsmooth)
  }, mc.cores = ncore) %>% dplyr::bind_rows()
  
  arrow_length <- parallel::mclapply(1:nrow(embedding), function(x){
    arr_len <- sqrt(( (embedding[x,1] - umap_new_unsmooth[x,1]) **2 + (embedding[x,2] - umap_new_unsmooth[x,2]) **2))
    return(arr_len)
  }, mc.cores = ncore) %>% unlist()
  arrow_length <- arrow_length/max(arrow_length)
  
  umap_new_unsmooth[arrow_length > quantile(arrow_length, cutoff),] <- embedding[arrow_length > quantile(arrow_length, cutoff),]
  
  # umap neighbors
  nnumap <- Seurat::FindNeighbors(embedding, k.param = 15, annoy.metric = 'euclidean', return.neighbor = TRUE)
  umapneighs <- nnumap@nn.idx
  
  umap_new <- parallel::mclapply(1:nrow(embedding), function(x){
    umap_n <- Matrix::colMeans(umap_new_unsmooth[umapneighs[x,],])
    return(umap_n)
  }, mc.cores = ncore) %>% dplyr::bind_rows()
  
  # umap_new[arrow_length > quantile(arrow_length, cutoff),] <- embedding[arrow_length > quantile(arrow_length, cutoff),]
  
  dX = umap_new - embedding
  
  return(list(dX, arrow_length))
}

calc_velocity_nm <- function(ndmat, nrmat, embedding, max_per_cell = 10, metric = 'cosine', ncore = 12, cutoff = 0.95){
  # return arrow coordinates
  nn <- Seurat::FindNeighbors(t(ndmat), query = t(nrmat), k.param = max_per_cell, annoy.metric = metric, return.neighbor = TRUE)
  dists <- nn@nn.dist
  neighs <- nn@nn.idx
  cutoff <- cutoff
  
  # for each cell. calculate the coordinates of its neighbors (most compatible neighbors)
  umap_new_unsmooth <- parallel::mclapply(1:nrow(embedding), function(x){
    umap_unsmooth <- Matrix::colMeans(embedding[neighs[x,],])
    return(umap_unsmooth)
  }, mc.cores = ncore) %>% dplyr::bind_rows()
  
  arrow_length <- parallel::mclapply(1:nrow(embedding), function(x){
    arr_len <- sqrt(( (embedding[x,1] - umap_new_unsmooth[x,1]) **2 + (embedding[x,2] - umap_new_unsmooth[x,2]) **2))
    return(arr_len)
  }, mc.cores = ncore) %>% unlist()
  arrow_length <- arrow_length/max(arrow_length)
  
  umap_new_unsmooth[arrow_length > quantile(arrow_length, cutoff),] <- embedding[arrow_length > quantile(arrow_length, cutoff),]
  
  
  dX = umap_new_unsmooth - embedding
  
  return(list(dX, arrow_length))
}

# todo
calc_velocity_restrict <- function(ndmat, nrmat, pseudotime, embedding, max_per_cell = 10, metric = 'cosine', ncore = 12, cutoff = 0.95, time_range = 10){
  # return arrow coordinates
  nn <- Seurat::FindNeighbors(t(ndmat), query = t(nrmat), k.param = max_per_cell, annoy.metric = metric, return.neighbor = TRUE)
  dists <- nn@nn.dist
  neighs <- nn@nn.idx
  cutoff <- cutoff
  
  
  
  parallel::mclapply(1:nrow(t(ndmat)), function(x){
    q_cells <- (1:nrow(t(ndmat))) [(pseudotime < pseudotime[x] + 10) & (pseudotime > pseudotime[x] - 10)]
    nn <- Seurat::FindNeighbors(t(ndmat)[x,], query = t(nrmat)[q_cells,], k.param = max_per_cell, annoy.metric = metric, return.neighbor = TRUE)
    
  }, mc.cores = ncore) %>% dplyr::bind_rows()
  
  # for each cell. calculate the coordinates of its neighbors (most compatible neighbors)
  umap_new_unsmooth <- parallel::mclapply(1:nrow(embedding), function(x){
    umap_unsmooth <- Matrix::colMeans(embedding[neighs[x,],])
    return(umap_unsmooth)
  }, mc.cores = ncore) %>% dplyr::bind_rows()
  
  arrow_length <- parallel::mclapply(1:nrow(embedding), function(x){
    arr_len <- sqrt(( (embedding[x,1] - umap_new_unsmooth[x,1]) **2 + (embedding[x,2] - umap_new_unsmooth[x,2]) **2))
    return(arr_len)
  }, mc.cores = ncore) %>% unlist()
  arrow_length <- arrow_length/max(arrow_length)
  
  umap_new_unsmooth[arrow_length > quantile(arrow_length, cutoff),] <- embedding[arrow_length > quantile(arrow_length, cutoff),]
  
  # umap neighbors
  nnumap <- Seurat::FindNeighbors(embedding, k.param = 15, annoy.metric = 'euclidean', return.neighbor = TRUE)
  umapneighs <- nnumap@nn.idx
  
  umap_new <- parallel::mclapply(1:nrow(embedding), function(x){
    umap_n <- Matrix::colMeans(umap_new_unsmooth[umapneighs[x,],])
    return(umap_n)
  }, mc.cores = ncore) %>% dplyr::bind_rows()
  
  # umap_new[arrow_length > quantile(arrow_length, cutoff),] <- embedding[arrow_length > quantile(arrow_length, cutoff),]
  
  dX = umap_new - embedding
  
  return(list(dX, arrow_length))
}



smooth_arrows <- function(embedding, embedding_new, smooth_w = 50, min_count = 5, coef = 2, draw_all = FALSE){
  
  x <- embedding[,1]
  y <- embedding[,2]
  u <- embedding_new[,1]
  v <- embedding_new[,2]
  
  x_window <-  (max(x) - min(x)) / smooth_w
  y_window <-  (max(y) - min(y)) / smooth_w
  
  xwin <- seq(from = min(x), to = max(x), by = x_window)
  xwin <- xwin[1:length(xwin)]
  ywin <- seq(from = min(y), to = max(y), by = y_window)
  ywin <- ywin[1:length(ywin)]
  
  
  ll <- lapply(xwin, function(i){
    lapply(ywin, function(j){
      idx <- (x > i) & (x <= (i + x_window)) & (y > j) & (y <= (j + y_window))
      if (sum(idx) < min_count & draw_all == FALSE){
        return(data.frame(x = NULL, y = NULL, Ex = NULL, Ey = NULL))
      }else if (sum(idx) < min_count & draw_all){
        return(data.frame(x = x[idx], y = y[idx], Ex = u[idx], Ey = v[idx]))
      }else{
        x_m <- mean(x[idx])
        y_m <- mean(y[idx])
        Ex_m <- mean(u[idx])
        Ey_m <- mean(v[idx])
        return(data.frame(x = x_m, y = y_m, Ex = Ex_m, Ey = Ey_m))
      }
    })
  }) 
  l2 <- unlist(sapply(ll, function(x) unlist(x)))
  x_ms <- l2[which(grepl("^x", names(l2)))]
  y_ms <- l2[which(grepl("^y", names(l2)))]
  Ex_ms <- l2[which(grepl("^Ex", names(l2)))]
  Ey_ms <- l2[which(grepl("^Ey", names(l2)))]
  # x_ms <- l2[which(1:length(l2) %% 4 == 1)]
  # y_ms <- l2[which(1:length(l2) %% 4 == 2)]
  # Ex_ms <- l2[which(1:length(l2) %% 4 == 3)]
  # Ey_ms <- l2[which(1:length(l2) %% 4 == 0)]
  

  
  x_end <- x_ms +  Ex_ms/coef
  y_end <- y_ms +  Ey_ms/coef
  return(data.frame(x_start = x_ms, y_start = y_ms, x_end = x_end, y_end = y_end))
}





# ---------------------------------- get gene expression matrix from project for ArchR ---------------------------------- #
getGEXMatrixFromProject <- function (ArchRProj = NULL, useMatrix = "GeneExpressionMatrix", useSeqnames = NULL, 
                                     verbose = TRUE, binarize = FALSE, threads = getArchRThreads()) {
  
  ArrowFiles <- getArrowFiles(ArchRProj)
  cellNames <- ArchRProj$cellNames
  avMat <- getAvailableMatrices(ArchRProj)
  
  message("Extracting GEX matrix from arrow...")
  seL <- lapply(seq_along(ArrowFiles), function(x) {
    
    allCells <- availableCells(ArrowFile = ArrowFiles[x], subGroup = useMatrix)
    allCells <- allCells[allCells %in% cellNames]
    if (length(allCells) != 0) {
      o <- getMatrixFromArrow2(ArrowFile = ArrowFiles[x], 
                               useMatrix = useMatrix, useSeqnames = NULL, 
                               cellNames = allCells, ArchRProj = ArchRProj, 
                               verbose = FALSE, binarize = FALSE)
      
      o
    }
    else {
      NULL
    }
  })
  
  cD <- lapply(seq_along(seL), function(x) {
    colData(seL[[x]])
  }) %>% Reduce("rbind", .)
  
  rD1 <- rowData(seL[[1]])
  rD <- lapply(seq_along(seL), function(x) {
    identical(rowData(seL[[x]]), rD1)
  }) %>% unlist %>% all
  if (!rD) {
    stop("Error with rowData being equal for every sample!")
  }
  
  rR1 <- rowRanges(seL[[1]])
  rR <- lapply(seq_along(seL), function(x) {
    identical(rowRanges(seL[[x]]), rR1)
  }) %>% unlist %>% all
  if (!rR) {
    stop("Error with rowRanges being equal for every sample!")
  }
  nAssays <- names(assays(seL[[1]]))
  asy <- lapply(seq_along(nAssays), function(i) {
    m <- lapply(seq_along(seL), function(j) {
      assays(seL[[j]])[[nAssays[i]]]
    }) %>% Reduce("cbind", .)
    m
  }) %>% SimpleList()
  names(asy) <- nAssays
  if (!is.null(rR1)) {
    se <- SummarizedExperiment(assays = asy, colData = cD, 
                               rowRanges = rR1)
    se <- sort(se)
  }
  else {
    se <- SummarizedExperiment(assays = asy, colData = cD, 
                               rowData = rD1)
  }
  rm(seL)
  gc()
  se
}

sampleName <- function (ArrowFile = NULL) 
{
  o <- h5closeAll()
  sampleName <- h5read(ArrowFile, paste0("Metadata/Sample"))
  o <- h5closeAll()
  return(sampleName)
}


availableCells <- function (ArrowFile = NULL, subGroup = NULL, passQC = TRUE) 
{
  if (is.null(subGroup)) {
    o <- h5closeAll()
    cellNames <- h5read(ArrowFile, "Metadata/CellNames")
    if (passQC) {
      passQC <- tryCatch({
        h5read(ArrowFile, "Metadata/PassQC")
      }, error = function(x) {
        rep(1, length(cellNames))
      })
      cellNames <- cellNames[which(passQC == 1)]
    }
    sampleName <- h5read(ArrowFile, paste0("Metadata/Sample"))
    o <- h5closeAll()
  }
  else {
    o <- h5closeAll()
    cellNames <- h5read(ArrowFile, paste0(subGroup, "/Info/CellNames"))
    sampleName <- h5read(ArrowFile, paste0("Metadata/Sample"))
    o <- h5closeAll()
  }
  return(paste0(sampleName, "#", cellNames))
}

availableSeqnames <- function (ArrowFiles = NULL, subGroup = "Fragments") 
{
  
  o <- h5closeAll()
  seqList <- lapply(seq_along(ArrowFiles), function(x) {
    seqnames <- h5ls(ArrowFiles[x]) %>% {
      .[.$group == paste0("/", subGroup), ]$name
    }
    seqnames <- seqnames[!grepl("Info", seqnames)]
    seqnames
  })
  if (!all(unlist(lapply(seq_along(seqList), function(x) identical(seqList[[x]], seqList[[1]]))))) {
    stop("Not All Seqnames Identical!")
  }
  o <- h5closeAll()
  return(paste0(seqList[[1]]))
}

getMetaData <- function (ArrowFile = NULL) 
{
  o <- h5closeAll()
  sampleName <- h5read(ArrowFile, paste0("Metadata/Sample"))
  arrowMD <- summarizeArrowContent(ArrowFile)$Metadata
  arrowMD <- arrowMD[which(arrowMD$dim == arrowMD$dim[arrowMD$name == 
                                                        "CellNames"]), ]
  md <- lapply(seq_len(nrow(arrowMD)), function(x) {
    dfx <- DataFrame(h5read(ArrowFile, paste0(arrowMD$group[x], 
                                              "/", arrowMD$name[x])))
    colnames(dfx) <- arrowMD$name[x]
    dfx
  }) %>% Reduce("cbind", .)
  md$CellNames <- paste0(sampleName, "#", md$CellNames)
  md$Sample <- Rle(sampleName, nrow(md))
  rownames(md) <- md$CellNames
  md <- md[, -which(colnames(md) == "CellNames")]
  md <- md[, order(colnames(md))]
  o <- h5closeAll()
  return(md)
}

summarizeArrowContent <- function (ArrowFile = NULL) 
{
  o <- h5closeAll()
  h5DF <- h5ls(ArrowFile)
  h5DF <- h5DF[-which(h5DF$group == "/"), ]
  groups <- stringr::str_split(h5DF$group, pattern = "/", 
                               simplify = TRUE)[, 2]
  groupList <- split(h5DF, groups)
  groupList2 <- lapply(seq_along(groupList), function(x) {
    groupDFx <- groupList[[x]]
    groupx <- gsub(paste0("/", names(groupList)[x]), "", 
                   groupDFx$group)
    if (all(groupx == "")) {
      groupDFx
    }
    else {
      subDF <- groupDFx[-which(groupx == ""), ]
      split(subDF, stringr::str_split(subDF$group, pattern = "/", 
                                      simplify = TRUE)[, 3])
    }
  })
  names(groupList2) <- names(groupList)
  o <- h5closeAll()
  return(groupList2)
}

getFeatureDF <- function (ArrowFiles = NULL, subGroup = "GeneExpressionMatrix", threads = 12) 
{
  threads <- min(threads, length(ArrowFiles))
  .helpFeatureDF <- function(ArrowFile = NULL, subGroup = NULL) {
    o <- h5closeAll()
    featureDF <- DataFrame(h5read(ArrowFile, paste0(subGroup, "/Info/FeatureDF")))
    featureDF$seqnames <- Rle(as.character(featureDF$seqnames))
    o <- h5closeAll()
    return(featureDF)
  }
  fdf <- .helpFeatureDF(ArrowFiles[1], subGroup = subGroup)
  if (length(ArrowFiles) > 1) {
    ArrowFiles <- ArrowFiles[-1]
    checkIdentical <- lapply(seq_along(ArrowFiles), 
                             function(x) {
                               fdfx <- .helpFeatureDF(ArrowFiles[x], subGroup = subGroup)
                               identical(fdfx, fdf)
                             }) %>% unlist %>% all
    
  }
  
  newOrder <- split(seq_len(nrow(fdf)), fdf$seqnames) %>% 
    {
      lapply(seq_along(.), function(x) .[[x]])
    } %>% Reduce("c", .)
  fdf[newOrder, ]
}

getMatrixFromArrowInternal <- function (ArrowFile = NULL, featureDF = NULL, binarize = NULL, 
                                        cellNames = NULL, useMatrix = "GeneExpressionMatrix", useIndex = FALSE) {
  if (is.null(featureDF)) {
    featureDF <- getFeatureDF(ArrowFile, useMatrix)
  }
  
  matColNames <- paste0(sampleName(ArrowFile), "#", h5read(ArrowFile, paste0(useMatrix, "/Info/CellNames")))
  
  if (!is.null(cellNames)) {
    idxCols <- which(matColNames %in% cellNames)
  }
  else {
    idxCols <- seq_along(matColNames)
  }
  seqnames <- unique(featureDF$seqnames)
  
  mat <- lapply(seq_along(seqnames), function(x) {
    seqnamex <- seqnames[x]
    featureDFx <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin%  seqnamex), ]
    
    idxRows <- featureDFx$idx
    
    j <- Rle(values = h5read(ArrowFile, paste0(useMatrix, "/", seqnamex, "/jValues")), lengths = h5read(ArrowFile, paste0(useMatrix, "/", seqnamex, "/jLengths")))
    # print(j)
    
    matchJ <- S4Vectors::match(j, idxCols, nomatch = 0)
    idxJ <- BiocGenerics::which(matchJ > 0)
    if (useIndex) {
      i <- h5read(ArrowFile, paste0(useMatrix, "/", seqnamex, 
                                    "/i"), index = list(idxJ, 1))
    }
    else {
      i <- h5read(ArrowFile, paste0(useMatrix, "/", seqnamex, 
                                    "/i"))[idxJ]
    }
    
    j <- matchJ[idxJ]
    matchI <- match(i, idxRows, nomatch = 0)
    idxI <- which(matchI > 0)
    i <- i[idxI]
    j <- j[idxI]
    # i <- matchI[idxI]
    
    if (!binarize) {
      x <- h5read(ArrowFile, paste0(useMatrix, "/", seqnamex, "/x"))[idxJ][idxI]
    }
    else {
      x <- rep(1, length(j))
    }
    
    mat <- Matrix::sparseMatrix(i = as.vector(i), j = j, x = x, dims = c(length(idxRows), length(idxCols)))
    rownames(mat) <- rownames(featureDFx)
    rm(matchI, idxI, matchJ, idxJ, featureDFx, idxRows)
    return(mat)
  }) %>% Reduce("rbind", .)
  
  
  o <- h5closeAll()
  colnames(mat) <- matColNames[idxCols]
  
  # mat <- mat[rownames(featureDF), , drop = FALSE]
  rownames(mat) <- NULL
  if (!is.null(cellNames)) {
    mat <- mat[, cellNames, drop = FALSE]
  }
  return(mat)
}

getMatrixFromArrow2 <- function (ArrowFile = NULL, useMatrix = "GeneExpressionMatrix", useSeqnames = NULL, 
                                 cellNames = NULL, ArchRProj = NULL, verbose = TRUE, binarize = FALSE) 
{
  
  sampleName <- sampleName(ArrowFile)
  seqnames <- availableSeqnames(ArrowFile, subGroup = useMatrix)
  featureDF <- getFeatureDF(ArrowFile, subGroup = useMatrix)
  
  if (!is.null(useSeqnames)) {
    seqnames <- seqnames[seqnames %in% useSeqnames]
  }
  if (length(seqnames) == 0) {
    stop("No seqnames available!")
  }
  featureDF <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% seqnames), ]
  
  if (!is.null(cellNames)) {
    allCells <- availableCells(ArrowFile = ArrowFile, subGroup = useMatrix)
    if (!all(cellNames %in% allCells)) {
      stop("cellNames must all be within the ArrowFile!!!!")
    }
  }
  mat <- getMatrixFromArrowInternal(ArrowFile = ArrowFile, featureDF = featureDF, 
                                    cellNames = cellNames, useMatrix = useMatrix, binarize = binarize, 
                                    useIndex = FALSE)
  
  matrixClass <- h5read(ArrowFile, paste0(useMatrix, "/Info/Class"))
  if (matrixClass == "Sparse.Assays.Matrix") {
    rownames(mat) <- paste0(featureDF$name)
    splitIdx <- split(seq_len(nrow(mat)), featureDF$seqnames)
    mat <- lapply(seq_along(splitIdx), function(x) {
      mat[splitIdx[[x]], , drop = FALSE]
    }) %>% SimpleList
    names(mat) <- names(splitIdx)
    featureDF <- featureDF[!duplicated(paste0(featureDF$name)), 
                           , drop = FALSE]
    featureDF <- featureDF[, which(colnames(featureDF) %ni% 
                                     "seqnames"), drop = FALSE]
    rownames(featureDF) <- paste0(featureDF$name)
  }else{
    mat <- SimpleList(mat)
    names(mat) <- useMatrix
  }
  colData <- getMetaData(ArrowFile)
  colData <- colData[colnames(mat[[1]]), , drop = FALSE]
  if (!is.null(ArchRProj)) {
    projColData <- getCellColData(ArchRProj)[rownames(colData), 
    ]
    colData <- cbind(colData, projColData[, colnames(projColData) %ni% 
                                            colnames(colData)])
  }
  rowData <- tryCatch({
    makeGRangesFromDataFrame(featureDF, keep.extra.columns = TRUE)
  }, error = function(x) {
    featureDF
  })
  se <- SummarizedExperiment(assays = mat, rowData = rowData, colData = colData)
  se
}


# ---------------------------------- SCENT ---------------------------------- #

# Interpolate a p-value from quantiles that should be "null scaled"
# q: bootstrap quantiles, centered so that under the null, theta = 0
# return two-sided p-value
interp_pval <-  function(q) {
  R <- length(q)
  tstar <- sort(q)
  zero <- findInterval(0, tstar)
  if(zero == 0 || zero == R) return(2/R) # at/beyond extreme values
  pval <- 2 * min(zero/R, (R-zero)/R)
  pval
}

# Derive a p-value from a vector of bootstrap samples using the "basic" calculation
# obs: observed value of parameter (using actual data)
# boot: vector of bootstraps
# return p-value
basic_p <-  function(obs, boot, null = 0){
  interp_pval(2 * obs - boot - null)
}

# Perform poisson regression: exprs ~ peak + covariates

# data: contains expr values and associated peak and covariates for a gene.
# idx: rows of the data to use: argument for boot function (bootstrapping)
# formula: user defined formula based on initialization in CreateSCENTObj Constructor
# return vector: (coefficient of the peak effect on gene, variance of peak effect on gene)
assoc_poisson = function(data, idx = seq_len(nrow(data)), formula){
  gg <-  glm(formula, family = 'poisson', data = data[idx,,drop = FALSE])
  c(coef(gg)['atac'], diag(vcov(gg))['atac'])
}


# Perform negative binomial regression: exprs ~ peak + covariates
# data: contains expr values and associated peak and covariates for a gene.
# idx: rows of the data to use: argument for boot function (bootstrapping)
# formula: user defined formula based on initialization in CreateSCENTObj Constructor
# return vector: (coefficient of the peak effect on gene, variance of peak effect on gene)
assoc_negbin = function(data, idx = seq_len(nrow(data)), formula){
  gg = glm.nb(formula, data = data[idx,,drop = FALSE])
  c(coef(gg)['atac'], diag(vcov(gg))['atac'])
}


CreateSCENTObj <- setClass(
  Class = "SCENT",
  slots = c(
    rna = 'dgCMatrix',
    atac = 'dgCMatrix',
    meta.data = 'data.frame',
    peak.info = 'data.frame',  ###Must be gene (1st column) then peak (2nd column)
    peak.info.list = 'list',
    covariates = 'character',
    celltypes = 'character',
    SCENT.result = 'data.frame'
  )
)

SCENT_algorithm <- function(object, celltype, ncore_link = 10, ncore_boot = 4, regr = "poisson", bin = TRUE){
  print(paste0("Total ", nrow(object@peak.info), " gene-peak pairs."), quote = FALSE)
  res <- parallel::mclapply(1:nrow(object@peak.info), function(n){
    gene <- object@peak.info[n,1]
    this_peak <- object@peak.info[n,2]
    atac_target <- data.frame(cell = colnames(object@atac), atac = object@atac[this_peak,])
    
    if(n %% 1000 == 0) {
      system(sprintf('echo "\n%s\n"', paste0("Processed ", n, " pairs.")))
    }
    #binarize peaks:
    if(bin){
      atac_target[atac_target$atac>0,]$atac <- 1
    }
    
    mrna_target <- object@rna[gene,]
    df <- data.frame(cell = names(mrna_target), exprs = as.numeric(mrna_target))
    df <- merge(df, atac_target,by = "cell")
    df <- merge(df, object@meta.data,by = "cell")
    
    df2 <- df[df[[object@celltypes]] == celltype,]
    
    nonzero_m <- length( df2$exprs[df2$exprs > 0] ) / length( df2$exprs )
    nonzero_a <- length( df2$atac[df2$atac > 0] ) / length( df2$atac )
    if(nonzero_m > 0.05 & nonzero_a > 0.05){
      #Run Regression Once Before Bootstrapping:
      res_var <- "exprs"
      pred_var <- c("atac", object@covariates)
      formula <- as.formula(paste(res_var, paste(pred_var, collapse = "+"), sep = "~"))
      
      
      #Estimated Coefficients Obtained without Bootstrapping:
      if(regr == "poisson"){
        base <- glm(formula, family = 'poisson', data = df2)
        coefs <- summary(base)$coefficients["atac",]
        assoc <- assoc_poisson
      } else if (regr == "negbin"){
        base <- glm.nb(formula, data = df2)
        coefs <- summary(base)$coefficients["atac",]
        assoc <- assoc_negbin
      }
      
      ###Iterative Bootstrapping Procedure: Estimate the Beta coefficients and associate a 2-sided p-value.
      bs <- boot::boot(df2, assoc, R = 100, formula = formula, stype = 'i', parallel = "multicore", ncpus = ncore_boot)
      p0 <- basic_p(bs$t0[1], bs$t[,1])
      if(p0 < 0.1){
        bs <- boot::boot(df2,assoc, R = 500, formula = formula,  stype = 'i', parallel = "multicore", ncpus = ncore_boot)
        p0 <- basic_p(bs$t0[1], bs$t[,1])
      }
      if(p0 < 0.05){
        bs <- boot::boot(df2,assoc, R = 2500, formula = formula,  stype = 'i', parallel = "multicore", ncpus = ncore_boot)
        p0 <- basic_p(bs$t0[1], bs$t[,1])
      }
      if(p0 < 0.01){
        bs <- boot::boot(df2,assoc, R = 25000, formula = formula,  stype = 'i', parallel = "multicore", ncpus = ncore_boot)
        p0 <- basic_p(bs$t0[1], bs$t[,1])
      }
      if(p0 < 0.001){
        bs <- boot::boot(df2,assoc, R = 50000, formula = formula, stype = 'i', parallel = "multicore", ncpus = ncore_boot)
        p0 <-  basic_p(bs$t0[1], bs$t[,1])
      }
      out <- data.frame(gene = gene, 
                        peak = this_peak, 
                        beta = coefs[1],
                        se = coefs[2], 
                        z = coefs[3], 
                        p = coefs[4], 
                        boot_basic_p = p0)
    }
  }, mc.cores = ncore_link) %>% Reduce("rbind", .)
  
  
  #Update the SCENT.result field of the constructor in R:
  object@SCENT.result <- res
  return(object)
}


getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
get_trajectory_heatmap_cell_comp <- function(proj, trajectory, cell_type, breaks = 100){
  traj <- as.numeric(proj@cellColData[[trajectory]])
  group_df <- data.frame(type = unlist(lapply(seq_along(1:breaks), function(x) {
    getmode(as.character(proj@cellColData[[cell_type]])[which(traj <= x & traj > x-1)])
  }))) 
  group_df <- group_df %>% dplyr::mutate(x = 1:nrow(group_df))
}

minmax <- function(x){(x - min(x))/(max(x) - min(x))}


# cluster by group
cluster_cdr3 <- function(cdr3_seqs, similarity_threshold = 0.5, max_workers = 8) {
  # kmer
  create_kmer_features <- function(seqs, k = 3) {
    seqs %>%
      map(~strsplit(.x, "")[[1]] %>%
            embed(dimension = k) %>%
            apply(1, paste, collapse = "")) %>%
      map(table) %>%
      bind_rows() %>%
      dplyr::mutate(across(everything(), ~replace_na(.x, 0))) %>%
      as.matrix()
  }
  
  # parallel
  registerDoParallel(cores = max_workers)
  
  # group sequences by length
  grouped_seqs <- tibble(seq = cdr3_seqs, len = nchar(seq)) %>%
    dplyr::mutate(length_group = cut(len, breaks = seq(0, max(len)+5, by = 5), include.lowest = TRUE))
  
  # in group kmean clustering
  clustered <- grouped_seqs %>%
    group_by(length_group) %>%
    dplyr::mutate(subgroup = {
      # set the minimum number of sequences in each group for grouping
      if (n() > 10000) {
        features <- create_kmer_features(seq)
        k <- min(ceiling(n()/2000), 100)  
        as.integer(kmeans(features, centers = k, iter.max = 20)$cluster)
      } else {
        rep(1, n())
      }
    }) %>%
    ungroup()
  
  # mini-group clustering
  result <- clustered %>%
    group_by(length_group, subgroup) %>%
    partition() %>%  
    dplyr::mutate(cluster = {
      if (n() > 1) {
        # distance matrix
        dist_mat <- stringdistmatrix(seq, method = "lv")
        
        # siminarity
        similarity_mat <- 1 - dist_mat / outer(len, len, pmax)
        
        # clustering
        adj_matrix <- similarity_mat > similarity_threshold
        g <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
        as.integer(igraph::components(g)$membership)
      } else {
        1L  # 单一序列单独成组
      }
    }) %>%
    collect() %>%
    ungroup() %>%
    dplyr::mutate(final_cluster = group_indices(., length_group, subgroup, cluster))
  
  stopImplicitCluster()
  
  return(result)
}


generate_all_contingency_tables <- function(meta_df, slot_cluster = NULL, slot_group = NULL) {
  all_groups <- unique(meta_df[[slot_group]])
  all_clusters <- unique(meta_df[[slot_cluster]])
  
  contingency_tables <- lapply(all_clusters, function(cluster){
    meta_df$cluster_binary <- ifelse(meta_df[[slot_cluster]] == cluster, "Target", "Other")
    cluster_tables <- lapply(all_groups, function(group){
      meta_df$group_binary <- ifelse(meta_df[[slot_group]] == group, "Target", "Other")
      
      contingency_table <- table(
        Cluster = meta_df$cluster_binary,
        Group = meta_df$group_binary
      )
      contingency_table
    })
    names(cluster_tables) <- all_groups
    cluster_tables
  })
  names(contingency_tables) <- all_clusters
  return(contingency_tables)

}


or_from_contingency_tables <- function(contingency_tables){
  or_results <- lapply(names(contingency_tables), function(cluster) {
    list_cluster <- contingency_tables[[cluster]]
    or_subcluster <- lapply(names(list_cluster), function(group){
      tab <- list_cluster[[group]]
      
      a <- tab["Target", "Target"] 
      b <- tab["Target", "Other"]   
      c <- tab["Other", "Target"] 
      d <- tab["Other", "Other"]
      
      or_value <- (a * d) / (b * c)
  
      se_log_or <- sqrt(1/a + 1/b + 1/c + 1/d)
      lower_ci <- exp(log(or_value) - 1.96 * se_log_or)
      upper_ci <- exp(log(or_value) + 1.96 * se_log_or)
      fisher_test <- fisher.test(tab)
      
      data.frame(
        cluster = cluster,
        target_group = group,
        OR = or_value,
        lower_CI = lower_ci,
        upper_CI = upper_ci,
        p_value = fisher_test$p.value,
        cells_in_target = a,
        cells_in_other = b,
        total_target = a + c,
        total_other = b + d
      )
    }) %>% Reduce("rbind", .)
  }) %>% Reduce("rbind", .)
  
  or_results$fdr <- p.adjust(or_results$p_value, method = "BH")
  return(or_results)
}

ro_re_from_df <- function(meta_df, slot_cluster = NULL, slot_group = NULL){
  data <- as.data.frame.array(table(meta_df[[slot_cluster]], meta_df[[slot_group]]))
  chi_sq_test <- chisq.test(data)
  exp <- chi_sq_test$expected
  Roe <- data/exp
  return(Roe)
}









