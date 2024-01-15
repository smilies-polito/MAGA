#BiocManager::install("AUCell")
#BiocManager::install("qusage")
#BiocManager::install("JASPAR2020")
#BiocManager::install("TFBSTools")
#BiocManager::install("chromVAR")
#devtools::install_github('cole-trapnell-lab/monocle3')
#devtools::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")
#devtools::install_github("stuart-lab/signac", ref = "develop")
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
#BiocManager::install("motifmatchr")
########################################################################################
#
# INSTALL REQUIRED PACKAGES IF NEEDED
# 
# N.B. This may require the installation of local libraries. Please check the README
# file of the project for a list of required packages
#
# ######################################################################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.15")

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils',
                       'HDF5Array', 'terra', 'ggrastr'))



if (!requireNamespace("devtools", quietly = TRUE)) 
  install.packages("devtools", dependencies = c("Depends"))
devtools::install_github('cole-trapnell-lab/monocle3')

BiocManager::install(c("Gviz", "GenomicRanges", "rtracklayer"))
devtools::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")

if (!requireNamespace("tidyverse", quietly = TRUE)) 
  install.packages("tidyverse", dependencies = c("Depends"))

if (!requireNamespace("aricode", quietly = TRUE)) 
  install.packages("aricode", dependencies = c("Depends"))

if (!requireNamespace("Seurat", quietly = TRUE)) 
  install.packages("Seurat", dependencies = c("Depends"))

if (!requireNamespace("Seurat", quietly = TRUE)) 
  install.packages("R.utils", dependencies = c("Depends"))

if (!requireNamespace("SeuratWrappers", quietly = TRUE)) 
  remotes::install_github('satijalab/seurat-wrappers')

if (!requireNamespace("Signac", quietly = TRUE)) 
  devtools::install_github("stuart-lab/signac", ref = "develop")

if (!requireNamespace("Matrix", quietly = TRUE)) 
  install.packages("Matrix", dependencies = c("Depends"))

if (!requireNamespace("TFBSTools", quietly = TRUE)) 
  BiocManager::install("TFBSTools")

if (!requireNamespace("motifmatchr", quietly = TRUE)) 
  BiocManager::install("motifmatchr")

if (!requireNamespace("chromVAR", quietly = TRUE)) 
  BiocManager::install("chromVAR")

if (!requireNamespace("JASPAR2020", quietly = TRUE)) 
  BiocManager::install("JASPAR2020")

if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) 
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

if (!requireNamespace("tidyr", quietly = TRUE)) 
  install.packages("tidyr", dependencies = c("Depends"))

if (!requireNamespace("dplyr", quietly = TRUE)) 
  install.packages("dplyr", dependencies = c("Depends"))

if (!requireNamespace("ggplot2", quietly = TRUE)) 
  install.packages("ggplot2", dependencies = c("Depends"))

if (!requireNamespace("qusage", quietly = TRUE)) 
  BiocManager::install("qusage")

#######################################################################################
#
# LOAD REQUIRED PACKAGES
# 
#######################################################################################


library(cicero)
library(data.table)
library(Matrix)
library(proxy)
library(reshape2)
library(readsparse)
library(aricode)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Signac)
library(Seurat)
library(SeuratObject)
library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg38)
library(motifmatchr)

##########

DefaultAssay(motif_activity) <- "EnhdPeaks"
motif_activity
genome=BSgenome.Hsapiens.UCSC.hg38
motif.matrix = NULL
verbose = TRUE


 
  if (!requireNamespace("chromVAR", quietly = TRUE)) {
    stop("Please install chromVAR. https://greenleaflab.github.io/chromVAR/")
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Please install SummarizedExperiment")
  }
  motif.matrix <- GetMotifData(object = motif_activity, slot = "data")
  
  peak.matrix <- GetAssayData(object = motif_activity, slot = "counts")
  if (!(all(peak.matrix@x == floor(peak.matrix@x)))) {
    warning("Count matrix contains non-integer values.
            ChromVAR should only be run on integer counts.")
  }
  idx.keep <- rowSums(x = peak.matrix) > 0
  peak.matrix <- peak.matrix[idx.keep, , drop = FALSE]
  motif.matrix <- motif.matrix[idx.keep, , drop = FALSE]
  peak.ranges <- granges(x = motif_activity)
  peak.ranges <- peak.ranges[idx.keep]
  chromvar.obj <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = peak.matrix),
    rowRanges = peak.ranges
  )
  if (verbose) {
    message("Computing GC bias per region")
  }
  chromvar.obj <- chromVAR::addGCBias(
    object = chromvar.obj,
    genome = genome
  )
  # Remove NA values https://github.com/GreenleafLab/chromVAR/issues/26
  row.data <- data.frame(SummarizedExperiment::rowData(x = chromvar.obj))
  row.data[is.na(x = row.data)] <- 0
  SummarizedExperiment::rowData(x = chromvar.obj) <- row.data
  if (verbose) {
    message("Selecting background regions")
  }
  bg <- chromVAR::getBackgroundPeaks(
    object = chromvar.obj
    
  )
  if (verbose) {
    message("Computing deviations from background")
  }
  register(SerialParam())
  dev <- chromVAR::computeDeviations(
    object = chromvar.obj,
    annotations = motif.matrix,
    background_peaks = bg
  )
  
  
  variability_prom <- computeVariability(dev)
  variability_prom <- variability_prom[variability_prom$name %in% is_exprs$motif, ] %>%  mutate(name = is_exprs$V1)
  s1 <- plotVariability(variability_prom, use_plotly = FALSE) + geom_hline(yintercept = 1) + labs(y = "Promoter Motif Enrichment Variability")
  
  
  variability_enhd <- computeVariability(dev)
  variability_enhd <- variability_enhd[variability_enhd$name %in% is_exprs$motif, ] %>%  mutate(name = is_exprs$V1)
  s2 <- plotVariability(variability_enhd, use_plotly = FALSE) +geom_hline(yintercept = 1) + labs(y = "Enhacer Motif Enrichment Variability")
s1+s2
ggsave(path = paste0("TMPResults/IMAGES/", "MISC"), filename = "all_variability.pdf", width = 1080, height = 500, units= "px",scale = 3.5)

  
  
  foot_enhd_sds <-  data.frame(sds= rowSds(as.matrix(TF_footprint[,2:15])))
  rownames(foot_enhd_sds) <- rownames(TF_footprint)
  
  res_df <- foot_enhd_sds %>%  mutate(rank = rank(-1 * foot_enhd_sds$sds, 
                                                  ties.method = "random")) %>% mutate(annotation = rownames(foot_enhd_sds))
  
  ggplot(res_df, aes_string(x = "rank", y = "sds",
                            label = "annotation")) + geom_point() + 
    xlab("Sorted TFs") + ylab("variability") + chromVAR_theme() 
  
  

  TF_mat <- as.matrix(TF_footprint_prom[,2:15])
  TF_mat <- (TF_mat - colMeans(TF_mat))/colSds(TF_mat)
  foot_prom_sds <- data.frame(sds= rowSds(TF_mat))
  res_df_prom <- foot_prom_sds %>%  mutate(rank = rank(-1 * foot_prom_sds$sds, 
                                                  ties.method = "random")) %>% mutate(annotation = rownames(foot_enhd_sds))
  
  ggplot(res_df_prom, aes_string(x = "rank", y = "sds",
                            label = "annotation")) + geom_point() + 
    xlab("Sorted TFs") + ylab("variability") + chromVAR_theme() 
  
  
  TF_mat <- as.matrix(TF_footprint[,2:15])
  TF_mat <- (TF_mat - colMeans(TF_mat))/colSds(TF_mat)
  foot_enhd_sds <- data.frame(sds= rowSds(TF_mat))
  res_df <- foot_enhd_sds %>%  mutate(rank = rank(-1 * foot_enhd_sds$sds, 
                                                  ties.method = "random")) %>% mutate(annotation = rownames(foot_enhd_sds))
  
  ggplot(res_df, aes_string(x = "rank", y = "sds",
                            label = "annotation")) + geom_point() + 
    xlab("Sorted TFs") + ylab("variability") + chromVAR_theme() 
  