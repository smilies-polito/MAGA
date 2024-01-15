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
  install.packages("Signac", dependencies = c("Depends"))

if (!requireNamespace("readsparse", quietly = TRUE)) 
  install.packages("readsparse", dependencies = c("Depends"))

if (!requireNamespace("leiden", quietly = TRUE)) 
  install.packages("leiden", dependencies = c("Depends"))

if (!requireNamespace("hdf5r", quietly = TRUE)) 
  install.packages("hdf5r", dependencies = c("Depends"))

install.packages('BiocManager', dependencies = c("Depends"))
BiocManager::install()
BiocManager::install('multtest')
install.packages('Seurat')

if (!requireNamespace("Seurat", quietly = TRUE)) 
  install.packages("Seurat", dependencies = TRUE)

if (!requireNamespace("Seurat", quietly = TRUE)) 
  install.packages("R.utils", dependencies = c("Depends"))

if (!requireNamespace("SeuratWrappers", quietly = TRUE)) 
  remotes::install_github('satijalab/seurat-wrappers')

if (!requireNamespace("SeuratDisk", quietly = TRUE))
  remotes::install_github("mojaveazure/seurat-disk")

if (!requireNamespace("data.table", quietly = TRUE)) 
  install.packages("data.table", dependencies = c("Depends"))

if (!requireNamespace("cellranger", quietly = TRUE)) 
  install.packages("cellranger", dependencies = c("Depends"))

if (!requireNamespace("readxl", quietly = TRUE)) 
  install.packages("readxl", dependencies = c("Depends"))

if (!requireNamespace("sctransform", quietly = TRUE)) 
  install.packages("sctransform", dependencies = c("Depends"))

if (!requireNamespace("class", quietly = TRUE)) 
  install.packages("class", dependencies = c("Depends"))

if (!requireNamespace("GenomeInfoDb", quietly = TRUE)) 
  install.packages("GenomeInfoDb", dependencies = c("Depends"))


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
#library(SeuratDisk)
library(muscat)

#######################################################################################
#
# GAGAM COMPUTATION FOR MOUSE
# 
#######################################################################################
#
# NEEDED ELEMENTS FROM THE DATASETS PROCESSING SCRIPTS:
# 
# processed_ATAC_cds: THE SCATAC-SEQ DATA OBJECT
# labeled_peaks: LABELED PEAKS TABLE
# nmax
# PBMC: THE SEURAT OBJECT WITH THE ATAC DATA IN IT
# connection_table: THE CONNECTION OBJECT FROM THE CO-ACCESSIBILITY COMPUTATION
# 
#
#######################################################################################

if (!(dir.exists("../TMPResults"))){
  dir.create("../TMPResults")
}

if (!(dir.exists("../TMPResults/Tables"))){
  dir.create("../TMPResults/Tables")
}

if (!(dir.exists("../TMPResults/IMAGES"))){
  dir.create("../TMPResults/IMAGES")
}


#############
matrix <- readMM("../TMPDATA/10x_PBMC_Multiome_Controller/filtered_feature_bc_matrix/matrix.mtx")
cells <- read.table("../TMPDATA/10x_PBMC_Multiome_Controller/filtered_feature_bc_matrix/barcodes.tsv")
features <- read.delim("../TMPDATA/10x_PBMC_Multiome_Controller/filtered_feature_bc_matrix/features.tsv", header=FALSE)

#division of the fetures between genes and peaks
genes <- features[features$V3 == "Gene Expression",]
colnames(genes)[2] <- "gene_short_name"
peaks <- features[features$V3 == "Peaks",]


row.names(matrix) <- features$V2
colnames(matrix) <- cells$V1

#creation of the two matrices, accordingly to the features
RNA_matrix <- matrix[genes$gene_short_name,]
ATAC_matrix <- matrix[peaks$V2,]
ATAC_matrix@Dimnames[[1]] <- gsub("-","_",ATAC_matrix@Dimnames[[1]])
ATAC_matrix@Dimnames[[1]] <- gsub(":","_",ATAC_matrix@Dimnames[[1]])

peaks$V1 <- gsub("-","_",peaks$V1)
peaks$V1 <- gsub(":","_",peaks$V1)
peaks$V2 <- gsub("-","_",peaks$V2)
peaks$V2 <- gsub(":","_",peaks$V2)
#write.table(peaks, file='10x_multiome_peaks.csv', quote=FALSE, sep='\t', col.names = NA)

   ####### RNA ANALYSIS #######

CDS_RNA <- new_cell_data_set(RNA_matrix)
rowData(CDS_RNA)$gene_short_name <- genes$gene_short_name
CDS_RNA <- detect_genes(CDS_RNA)
CDS_RNA <- estimate_size_factors(CDS_RNA)
#preprocessing consisting in normalization, scaling, and dimansional reduction (both LSI and PCA)
CDS_RNA <- preprocess_cds(CDS_RNA, method = "LSI")
CDS_RNA <- preprocess_cds(CDS_RNA, method = "PCA")
CDS_RNA <- reduce_dimension(CDS_RNA, reduction_method = 'UMAP', 
                            preprocess_method = "LSI")
#clustering of the cells, based on the UMAP
CDS_RNA <- cluster_cells(CDS_RNA, resolution=0.4e-3)

ATAC_matrix@x[ATAC_matrix@x > 0] <- 1
#Creation of the CDS object for the ATAC data
CDS_ATAC <- new_cell_data_set(ATAC_matrix)
rowData(CDS_ATAC)$gene_short_name <- peaks$V2
#the process and function are totally analogous to  before
CDS_ATAC <- detect_genes(CDS_ATAC)
CDS_ATAC <- estimate_size_factors(CDS_ATAC)
CDS_ATAC <- preprocess_cds(CDS_ATAC, method = "LSI")
#CDS_ATAC <- preprocess_cds(CDS_ATAC, method = "PCA")
CDS_ATAC <- reduce_dimension(CDS_ATAC, reduction_method = 'UMAP', 
                             preprocess_method = "LSI")
CDS_ATAC <- cluster_cells(CDS_ATAC, resolution=0.5e-3)

####### SIGNAC #####
genome_ref = read.table("../DATA/genomes/hg38/hg38.p13.chrom.sizes.txt")
genome_ref <- genome_ref[1:24,]


hg38 <- genome_ref[1:24,]
hg38 <- Seqinfo(hg38$V1, seqlengths= hg38$V2)
hg38@genome[] <- "hg38"


PBMC_assay <- CreateChromatinAssay(
  counts = ATAC_matrix,
  sep = c("_", "_"),
  genome = hg38,
  min.cells = 1
)


PBMC <- CreateSeuratObject(
  counts = PBMC_assay,
  assay = 'peaks',
  project = 'ATAC'
)

 
labeled_peaks <- read.csv("../TMPDATA/10x_PBMC_Multiome_Controller/labeled peaks/encodeCcreCombined_hg38_ucscLabel_classifiedPeaks.csv", sep = "\t")
#labeled_peaks_multi <- read.csv("C:/Users/loren/IWBBIO/data/labeled peaks/classifiedPeaks_multiCols.csv", sep = "\t")

nmax <- max(stringr::str_count(labeled_peaks$encodeCcreCombined_hg38_ucscLabel, "\t")) + 1

labeled_peaks <- separate(labeled_peaks, col = encodeCcreCombined_hg38_ucscLabel, sep = "\t", into = paste0("RegFunc", seq_len(nmax)))

labeled_peaks$site_names <- paste0(labeled_peaks$X.chrom, "_", labeled_peaks$chromStart, "_", labeled_peaks$chromEnd)

labeled_peaks <- labeled_peaks[labeled_peaks$site_names %in% rownames(fData(CDS_ATAC)),]
labeled_peaks <- labeled_peaks[!duplicated(labeled_peaks),]

######

genome_ref = read.table("Data/genomes/hg38/hg38.p13.chrom.sizes.txt")
genome_ref <- genome_ref[1:24,]


#LOAD GENOME REFERENCE AND CHROMOSOMES INORMATION
refseq_anno <-  rtracklayer::readGFF("Data/genomes/hg38/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz")
chr2acc <- read.csv("Data/genomes/hg38/chr2acc.txt", sep = "\t")

#SUBSET TO ONLY GENES
refseq_gene_anno <- refseq_anno[refseq_anno$seqid %in% chr2acc$Accession.version,]
refseq_gene_anno$seqid <- as.factor(as.character(refseq_gene_anno$seqid))
levels(refseq_gene_anno$seqid) <- chr2acc$X.Chromosome
refseq_gene_anno$seqid <- paste0("chr", refseq_gene_anno$seqid)
refseq_gene_anno <- refseq_gene_anno[refseq_gene_anno$type == "gene",]
refseq_gene_anno <- refseq_gene_anno[refseq_gene_anno$gene_biotype %in% c("protein_coding"),]

#DEFINING TSS AND FINDING TASS PEAKS
pos <- subset(refseq_gene_anno, strand == "+")
pos <- pos[order(pos$start),] 
pos$end <- pos$start + 1 
neg <- subset(refseq_gene_anno, strand == "-")
neg <- neg[order(neg$start, decreasing = TRUE),] 
neg$start <- neg$end - 1

refseq_gene_annotation_sub <- rbind(pos, neg)
refseq_gene_annotation_sub <- refseq_gene_annotation_sub[,c("seqid", "start", "end", "gene_id")]

names(refseq_gene_annotation_sub)[4] <- "gene"
CDS_ATAC <- annotate_cds_by_site(CDS_ATAC, refseq_gene_annotation_sub)



#CONSTRUCT PEAKS INFORMATION FROM LABELED PEAKS
peaks_multi_info <- fData(CDS_ATAC)
peaks_multi_info$site_name <- rownames(fData(CDS_ATAC))
peaks_multi_info$is.prom <- FALSE
peaks_multi_info$is.enhD <- FALSE
peaks_multi_info <- cbind(peaks_multi_info,labeled_peaks[7:(nmax+6)])

#PROM PEAKS
ppp <- peaks_multi_info
ppp <- as.data.frame(peaks_multi_info)
ppp[ppp == "prom"] <- TRUE

peaks_multi_info$is.prom <- apply(ppp, 1, any)
peaks_multi_info[is.na(peaks_multi_info$is.prom),]$is.prom <- FALSE
peaks_prom <- peaks_multi_info[peaks_multi_info$is.prom == "TRUE",]
prom_list <- rownames(peaks_prom)
non_prom_list <- rownames(peaks_multi_info)
non_prom_list <- setdiff(non_prom_list, prom_list)
non_prom_peaks <- peaks_multi_info[non_prom_list,]
non_prom_peaks[!is.na(non_prom_peaks$gene),]
peaks_prom[!is.na(peaks_prom$gene),]


#ENHD PEAKS
ppp <- peaks_multi_info
ppp <- as.data.frame(peaks_multi_info)
ppp[ppp == "enhD"] <- TRUE

peaks_multi_info$is.enhD <- FALSE
peaks_multi_info$is.enhD <- apply(ppp,1, any)
peaks_multi_info[is.na(peaks_multi_info$is.enhD),]$is.enhD <- FALSE

enhD_peaks_list <- rownames(peaks_multi_info[peaks_multi_info$is.enhD == TRUE,])
enhD_peaks_list <- setdiff(enhD_peaks_list, prom_list)


peaks_multi_info$PT <- !(is.na(peaks_multi_info$gene)) | as.vector(peaks_multi_info$is.prom)


####### PROMOTER CONTRIBUTION ########

refseq_peaks_prom <- peaks_multi_info[peaks_multi_info$PT == "TRUE",]
refseq_prom_list <- rownames(refseq_peaks_prom)

hg38 <- genome_ref[1:24,]
hg38 <- Seqinfo(hg38$V1, seqlengths= hg38$V2)
hg38@genome[] <- "hg38"


refseq_GRanges <- makeGRangesFromDataFrame(refseq_gene_anno, seqinfo = hg38, seqnames.field = "seqid", keep.extra.columns = TRUE)

refseq_gene_closest_to_prom <- ClosestFeature(PBMC, gsub("_", "-", refseq_prom_list), annotation = refseq_GRanges)
refseq_gene_closest_to_prom <- refseq_gene_closest_to_prom[refseq_gene_closest_to_prom$distance <= 500,]

refseq_gene_closest_to_prom$query_region <- gsub("-","_",refseq_gene_closest_to_prom$query_region)

refseq_near_prom_list <- gsub("-","_",refseq_gene_closest_to_prom$query_region)

 
refseq_peaks_prom <-  refseq_peaks_prom[refseq_near_prom_list,]
refseq_peaks_prom$gene_anno <- refseq_gene_closest_to_prom$gene
refseq_peaks_prom$site_name <- rownames(refseq_peaks_prom)


refseq_promoter_peak_table <- refseq_peaks_prom[,c("gene","gene_anno","site_name")]
refseq_promoter_peak_table[is.na(refseq_promoter_peak_table$gene),]$gene <- refseq_promoter_peak_table[is.na(refseq_promoter_peak_table$gene),]$gene_anno
prom_gene_name <- levels(factor(refseq_promoter_peak_table$gene))

refseq_promoter_gene_mat <-
  Matrix::sparseMatrix(j=as.numeric(factor(refseq_promoter_peak_table$site_name)),
                       i=as.numeric(factor(refseq_promoter_peak_table$gene)),
                       x=1)

refseq_accessibility_mat <- exprs(CDS_ATAC)
refseq_accessibility_mat@x[refseq_accessibility_mat@x>0] <-1

refseq_promoter_activity_scores <- refseq_accessibility_mat[refseq_near_prom_list,, drop=FALSE]

colnames(refseq_promoter_gene_mat) = levels(factor(refseq_promoter_peak_table$site_name))
row.names(refseq_promoter_gene_mat) = levels(factor(refseq_promoter_peak_table$gene))
refseq_promoter_gene_mat <- refseq_promoter_gene_mat[,row.names(refseq_promoter_activity_scores)]

#refseq_first_gene_matrix2 <- refseq_promoter_gene_mat %*% refseq_promoter_activity_scores
refseq_first_gene_matrix <- refseq_promoter_gene_mat %*% refseq_promoter_activity_scores
refseq_first_gene_matrix2 <- refseq_first_gene_matrix
refseq_first_gene_matrix2@x[refseq_first_gene_matrix2@x > 0] <- 1


####### EXON CONTRIBUTION #########

prom_list <- rownames(refseq_promoter_peak_table)
non_prom_list <- rownames(peaks_multi_info)
non_prom_list <- setdiff(non_prom_list, prom_list)
non_prom_list <- setdiff(non_prom_list, enhD_peaks_list)


refseq_exon_anno <- refseq_anno[refseq_anno$seqid %in% chr2acc$Accession.version,]
refseq_exon_anno$seqid <- as.factor(as.character(refseq_exon_anno$seqid))
levels(refseq_exon_anno$seqid) <- chr2acc$X.Chromosome
refseq_exon_anno$seqid <- paste0("chr", refseq_exon_anno$seqid)
refseq_exon_anno <- refseq_exon_anno[refseq_exon_anno$type == "exon",]
refseq_GRanges_exon <- makeGRangesFromDataFrame(refseq_exon_anno, seqinfo = hg38, seqnames.field = "seqid", keep.extra.columns = TRUE)

intragenetic_non_prom_peaks <- ClosestFeature(PBMC, gsub("_", "-", non_prom_list), refseq_GRanges_exon)
intragenetic_non_prom_peaks <- intragenetic_non_prom_peaks[intragenetic_non_prom_peaks$distance == 0,]
intragenetic_non_prom_peaks_list <- gsub("-","_",intragenetic_non_prom_peaks$query_region)

intragenetic_peaks <-  peaks_multi_info[intragenetic_non_prom_peaks_list,]
intragenetic_peaks$gene_anno <- intragenetic_non_prom_peaks$gene
intragenetic_peaks$site_name <- rownames(intragenetic_peaks)

intragenetic_peak_table <- intragenetic_peaks[,c("gene_anno","site_name")]

refseq_gene_annotation_sub <- refseq_gene_annotation_sub[refseq_gene_annotation_sub$gene %in% intragenetic_peak_table$gene_anno,]
refseq_gene_annotation_sub$TSS  <- paste0(refseq_gene_annotation_sub$seqid,"_",refseq_gene_annotation_sub$start,"_", refseq_gene_annotation_sub$end)

intragenetic_peak_table <- intragenetic_peak_table[intragenetic_peak_table$gene_anno %in% prom_gene_name,]

prova <- intragenetic_peak_table
prova <- as.data.frame(prova)
prova <- prova %>% rowwise() %>% mutate(TSS = refseq_gene_annotation_sub[refseq_gene_annotation_sub$gene == gene_anno,]$TSS) 
intragenetic_peak_table$TSS <- prova$TSS

split_peak_names <- function(inp) {
  out <- stringr::str_split_fixed(stringi::stri_reverse(inp), 
                                  ":|-|_", 3)
  out[,1] <- stringi::stri_reverse(out[,1])
  out[,2] <- stringi::stri_reverse(out[,2])
  out[,3] <- stringi::stri_reverse(out[,3])
  out[,c(3,2,1), drop=FALSE]
}

if ("dist" %in% colnames(intragenetic_peak_table) == FALSE) {
  Peak1_cols <- split_peak_names(intragenetic_peak_table$site_name)
  Peak2_cols <- split_peak_names(intragenetic_peak_table$TSS)
  Peak1_bp <- round((as.integer(Peak1_cols[,3]) +
                       as.integer(Peak1_cols[,2])) / 2)
  Peak2_bp <- round((as.integer(Peak2_cols[,3]) +
                       as.integer(Peak2_cols[,2])) / 2)
  intragenetic_peak_table$dist <- abs(Peak2_bp - Peak1_bp)
}
rm(Peak1_cols)
rm(Peak2_cols)
rm(Peak1_bp)
rm(Peak2_bp)
intragenetic_peak_table$dist <- exp(-1*intragenetic_peak_table$dist/5000)

allgene_peaks_list <- rownames(intragenetic_peak_table)

refseq_allgene_gene_mat <-
  Matrix::sparseMatrix(j=as.numeric(factor(intragenetic_peak_table$site_name)),
                       i=as.numeric(factor(intragenetic_peak_table$gene_anno)),
                       x=intragenetic_peak_table$dist)

colnames(refseq_allgene_gene_mat) = levels(factor(intragenetic_peak_table$site_name))
row.names(refseq_allgene_gene_mat) = levels(factor(intragenetic_peak_table$gene_anno))
refseq_allgene_gene_mat <- refseq_allgene_gene_mat[row.names(refseq_allgene_gene_mat) %in% prom_gene_name,]



refseq_allgene_accessibility_mat <- exprs(CDS_ATAC)
refseq_allgene_accessibility_mat@x[refseq_allgene_accessibility_mat@x>0] <- 1
refseq_allgene_activity_scores <- refseq_allgene_accessibility_mat[allgene_peaks_list,, drop=FALSE]



refseq_allgene_gene_mat <- refseq_allgene_gene_mat[,row.names(refseq_allgene_activity_scores)]

refseq_allgene_gene_matrix <- refseq_allgene_gene_mat %*% refseq_allgene_activity_scores

refseq_allgene_gene_matrix_bin <-  refseq_allgene_gene_matrix * refseq_first_gene_matrix2[rownames(refseq_allgene_gene_matrix),]


####### FINAL CONSTRUCTION ########

make_sparse_matrix <- function(data,
                               i.name = "Peak1",
                               j.name = "Peak2",
                               x.name = "value") {
  if(!i.name %in% names(data) |
     !j.name %in% names(data) |
     !x.name %in% names(data)) {
    stop('i.name, j.name, and x.name must be columns in data')
  }
  
  data$i <- as.character(data[,i.name])
  data$j <- as.character(data[,j.name])
  data$x <- data[,x.name]
  
  if(!class(data$x) %in%  c("numeric", "integer"))
    stop('x.name column must be numeric')
  
  peaks <- data.frame(Peak = unique(c(data$i, data$j)),
                      index = seq_len(length(unique(c(data$i, data$j)))))
  
  data <- data[,c("i", "j", "x")]
  
  data <- rbind(data, data.frame(i=peaks$Peak, j = peaks$Peak, x = 0))
  data <- data[!duplicated(data[,c("i", "j", "x")]),]
  data <- data.table::as.data.table(data)
  peaks <- data.table::as.data.table(peaks)
  data.table::setkey(data, "i")
  data.table::setkey(peaks, "Peak")
  data <- data[peaks]
  data.table::setkey(data, "j")
  data <- data[peaks]
  data <- as.data.frame(data)
  
  data <- data[,c("index", "i.index", "x")]
  data2 <- data
  names(data2) <- c("i.index", "index", "x")
  
  data <- rbind(data, data2)
  
  data <- data[!duplicated(data[,c("index", "i.index")]),]
  data <- data[data$index >= data$i.index,]
  
  sp_mat <- Matrix::sparseMatrix(i=as.numeric(data$index),
                                 j=as.numeric(data$i.index),
                                 x=data$x,
                                 symmetric = TRUE)
  
  colnames(sp_mat) <- peaks[order(peaks$index),]$Peak
  row.names(sp_mat) <- peaks[order(peaks$index),]$Peak
  return(sp_mat)
}


connection_table <- readRDS("../TMPDATA/10x_PBMC_Multiome_Controller/conns/conns_10k_Multiome")

accessibility_mat <- exprs(CDS_ATAC)
accessibility_mat@x[accessibility_mat@x>0] <-1
#rownames(accessibility_mat) <- gsub(":", "_",rownames(accessibility_mat))
#rownames(accessibility_mat) <- gsub("-", "_",rownames(accessibility_mat))

con_val <- connection_table[connection_table$coaccess > 0,]
con_val <- con_val[!is.na(con_val$coaccess),]
coaccess <- signif(mean(con_val$coaccess), digits = 2)

if ("dist" %in% colnames(connection_table) == FALSE) {
  Peak1_cols <- split_peak_names(connection_table$Peak1)
  Peak2_cols <- split_peak_names(connection_table$Peak2)
  Peak1_bp <- round((as.integer(Peak1_cols[,3]) +
                       as.integer(Peak1_cols[,2])) / 2)
  Peak2_bp <- round((as.integer(Peak2_cols[,3]) +
                       as.integer(Peak2_cols[,2])) / 2)
  connection_table$dist <- abs(Peak2_bp - Peak1_bp)
}


rm(Peak1_cols)
rm(Peak2_cols)
rm(Peak1_bp)
rm(Peak2_bp)

site_weights = NULL
if (is.null(site_weights)) {
  site_weights <- Matrix::rowMeans(accessibility_mat) /
    Matrix::rowMeans(accessibility_mat)
  site_weights[names(site_weights)] <- 1
}
site_names <- names(site_weights)
site_weights <- as(Matrix::Diagonal(x=as.numeric(site_weights)),
                   "sparseMatrix")
row.names(site_weights) <- site_names
colnames(site_weights) <- site_names


promoter_peak_table <- refseq_promoter_peak_table[, c("site_name", "gene", "gene_anno")]



promoter_peak_table$site_name <- as.character(row.names(promoter_peak_table))
promoter_peak_table <- promoter_peak_table[!is.na(promoter_peak_table$gene),]
promoter_peak_table <- promoter_peak_table[,c("site_name", "gene")]

colnames(promoter_peak_table) <- c("peak", "gene")

dist_thresh=400000
coaccess
prom_enhD <- connection_table[(connection_table$Peak1 %in%
                                 promoter_peak_table$peak &
                                 connection_table$Peak2 %in%
                                 enhD_peaks_list) | (connection_table$Peak2 %in%
                                                       promoter_peak_table$peak &
                                                       connection_table$Peak1 %in%
                                                       enhD_peaks_list),]

#prom_enhD <- connection_table[(connection_table$Peak1 %in%
#                                 promoter_peak_table$peak &
#                                 connection_table$Peak2 %in%
#                                 enhD_peaks_list),]

prom_enhD <- prom_enhD[prom_enhD$coaccess >= coaccess & prom_enhD$dist <= dist_thresh,]
prom_enhD <- prom_enhD[!duplicated(prom_enhD),]

prom_enhD <- prom_enhD[,c("Peak1", "Peak2", "coaccess")]
prom_enhD <- prom_enhD[!duplicated(prom_enhD),]

prom_enhD$Peak1 <- as.character(prom_enhD$Peak1)
prom_enhD$Peak2 <- as.character(prom_enhD$Peak2)


prom_enhD <- rbind(prom_enhD,
                   data.frame(Peak1=unique(promoter_peak_table$peak),
                              Peak2=unique(promoter_peak_table$peak),
                              coaccess=0))

prom_enhD_connectivity <- make_sparse_matrix(prom_enhD, x.name = "coaccess")

promoter_conn_matrix <-
  prom_enhD_connectivity[unique(promoter_peak_table$peak),]


promoter_safe_sites <- intersect(rownames(promoter_conn_matrix),
                                 row.names(accessibility_mat))
distal_safe_sites <- intersect(colnames(promoter_conn_matrix),
                               row.names(accessibility_mat))
distal_safe_sites <- setdiff(distal_safe_sites, promoter_safe_sites)

promoter_access_mat_in_cicero_map <- accessibility_mat[promoter_safe_sites,, drop=FALSE]

distal_activity_scores <- accessibility_mat[distal_safe_sites,, drop=FALSE]

scaled_site_weights <- site_weights[distal_safe_sites,distal_safe_sites, drop=FALSE]
total_linked_site_weights <- promoter_conn_matrix[,distal_safe_sites, drop=FALSE] %*%
  scaled_site_weights
total_linked_site_weights <- 1/Matrix::rowSums(total_linked_site_weights,
                                               na.rm=TRUE)
total_linked_site_weights[is.finite(total_linked_site_weights) == FALSE] <- 0
total_linked_site_weights[is.na(total_linked_site_weights)] <- 0
total_linked_site_weights[is.nan(total_linked_site_weights)] <- 0
total_linked_site_weights <- Matrix::Diagonal(x=total_linked_site_weights)
scaled_site_weights <- total_linked_site_weights %*%
  promoter_conn_matrix[,distal_safe_sites, drop=FALSE] %*%
  scaled_site_weights
scaled_site_weights@x[scaled_site_weights@x > 1] <- 1

distal_activity_scores <- scaled_site_weights %*% distal_activity_scores

distal_activity_scores <- distal_activity_scores[row.names(promoter_access_mat_in_cicero_map),, drop=FALSE]

promoter_activity_scores <-  distal_activity_scores  + refseq_promoter_activity_scores 

promoter_gene_mat <-
  Matrix::sparseMatrix(j=as.numeric(factor(promoter_peak_table$peak)),
                       i=as.numeric(factor(promoter_peak_table$gene)),
                       x=1)
colnames(promoter_gene_mat) = levels(factor(promoter_peak_table$peak))
row.names(promoter_gene_mat) = levels(factor(promoter_peak_table$gene))
promoter_gene_mat <- promoter_gene_mat[,row.names(promoter_activity_scores)]
gene_activity_scores <- promoter_gene_mat %*% promoter_activity_scores


############### META-ANALYSIS ################### 
###### Proportions #####

enhd_contribution <- promoter_gene_mat %*% distal_activity_scores

mean(refseq_allgene_gene_matrix@x)
mean(refseq_first_gene_matrix@x)
mean(enhd_contribution@x)

info_table <- data.frame(gene = rownames(refseq_first_gene_matrix))
sum <- rowSums(refseq_promoter_gene_mat)
info_table$num_prom <-  as.data.frame(sum)$sum
rownames(info_table) <- info_table$gene

temp_bin <- refseq_first_gene_matrix
temp_bin@x[temp_bin@x>0] <- 1
info_table$active_prom_num <- rowSums(temp_bin)
info_table$active_prom_perc <- rowSums(temp_bin)/dim(refseq_first_gene_matrix)[2]


temp_bin <- refseq_allgene_gene_mat
temp_bin@x[temp_bin@x>0] <- 1
sum <- rowSums(temp_bin)
info_table$num_ex <- 0
info_table[rownames(temp_bin),]$num_ex <- info_table[rownames(temp_bin),]$num_ex + as.data.frame(sum)$sum

temp_bin <- refseq_allgene_gene_matrix
temp_bin@x[temp_bin@x>0] <- 1
info_table$active_ex_num <- 0
info_table[rownames(temp_bin),]$active_ex_num <- info_table[rownames(temp_bin),]$active_ex_num + rowSums(temp_bin)
info_table$active_ex_perc <- 0
info_table[rownames(temp_bin),]$active_ex_perc <- info_table[rownames(temp_bin),]$active_ex_perc + rowSums(temp_bin)/dim(refseq_allgene_gene_matrix)[2]



temp_bin <- promoter_conn_matrix[,distal_safe_sites, drop=FALSE]
temp_bin@x[temp_bin@x > 0] <- 1
sume <- rowSums(temp_bin)
sume <- as.data.frame(sume)
sume <- as.sparse(sume)
sum <- refseq_promoter_gene_mat %*%  sume
sum <- as.data.frame(sum)
info_table$num_enh <- sum$sume

temp_bin <- enhd_contribution
temp_bin@x[temp_bin@x > 0] <- 1
info_table$active_enh_num <- rowSums(temp_bin)
info_table$active_enh_perc <-  rowSums(temp_bin)/dim(temp_bin)[2]

info_table$active_enh <- NULL

active_prom_name <- rownames(info_table)
active_no_expr <- active_prom_name[!(rownames(info_table) %in% rownames(RNA_matrix))]
active_prom_name <- active_prom_name[!(active_prom_name %in% active_no_expr)]

temp_bin <- RNA_matrix[active_prom_name,]
temp_bin@x[temp_bin@x > 0] <- 1

info_table$expr_num <- 0
info_table[rownames(temp_bin),]$expr_num <- info_table[active_prom_name,]$expr_num + rowSums(temp_bin)
info_table$expr_perc <- 0
info_table[rownames(temp_bin),]$expr_perc <- info_table[active_prom_name,]$expr_perc + rowSums(temp_bin)/dim(temp_bin)[2]

write.table(info_table,file = "../TMPResults/Tables/Gene_Peak_info.csv")

######### Seurat object creation ############
DefaultAssay(SEU_RNA) <- "RNA"

SEU_RNA <-  CreateSeuratObject(RNA_matrix)
SEU_RNA <- NormalizeData(SEU_RNA)
SEU_RNA <- FindVariableFeatures(SEU_RNA)
SEU_RNA <- FindVariableFeatures(SEU_RNA, selection.method = "mvp")
#VariableFeaturePlot(SEU_RNA, log = FALSE, selection.method = "vst")
#VariableFeaturePlot(SEU_RNA, selection.method = "mvp", log = TRUE)
SEU_RNA <- ScaleData(SEU_RNA)
SEU_RNA <- RunPCA(SEU_RNA)
SEU_RNA <- FindNeighbors(SEU_RNA)
SEU_RNA <- FindClusters(SEU_RNA, resolution = 0.8)
SEU_RNA <- RunUMAP(SEU_RNA,dims = 1:10)
DimPlot(SEU_RNA)
AggExp <- AggregateExpression(SEU_RNA)

atac <- CreateAssayObject(refseq_first_gene_matrix)
SEU_RNA[["PROM"]] <- atac
DefaultAssay(SEU_RNA) <- "PROM"
SEU_RNA <- RunTFIDF(SEU_RNA)
SEU_RNA <- FindVariableFeatures(SEU_RNA)
SEU_RNA <- FindVariableFeatures(SEU_RNA, selection.method = "mvp")
VariableFeaturePlot(SEU_RNA)
SEU_RNA <- ScaleData(SEU_RNA)
SEU_RNA <- RunPCA(SEU_RNA)
SEU_RNA <- FindNeighbors(SEU_RNA)
SEU_RNA <- FindClusters(SEU_RNA) 
SEU_RNA <- RunUMAP(SEU_RNA,dims = 1:10)
AggProm <- AggregateExpression(SEU_RNA, group.by = "RNA_snn_res.0.8", features = active_prom_name)
AvgProm <- AverageExpression(SEU_RNA, group.by = "RNA_snn_res.0.8", features = active_prom_name)
DimPlot(SEU_RNA)
expr_markers <- FindAllMarkers(SEU_RNA)

atac <- CreateAssayObject(refseq_allgene_gene_matrix)
SEU_RNA[["EXON"]] <- atac
DefaultAssay(SEU_RNA) <- "EXON"
SEU_RNA <- NormalizeData(SEU_RNA)
SEU_RNA <- FindVariableFeatures(SEU_RNA)
SEU_RNA <- FindVariableFeatures(SEU_RNA, selection.method = "mvp")
SEU_RNA <- subset(SEU_RNA, features = active_prom_name)
Idents(SEU_RNA) <- "RNA_snn_res.0.8"

atac <- CreateAssayObject(enhd_contribution)
SEU_RNA[["ENHD"]] <- atac
DefaultAssay(SEU_RNA) <- "ENHD"
SEU_RNA <- NormalizeData(SEU_RNA)
SEU_RNA <- FindVariableFeatures(SEU_RNA)
SEU_RNA <- FindVariableFeatures(SEU_RNA, selection.method = "mvp")
SEU_RNA <- subset(SEU_RNA, features = active_prom_name)



###### Aggregated cells definition ########

variance_table <- data.frame(active_prom_name)
variance_ratio_table <- data.frame(active_prom_name)
mean_table <- data.frame(active_prom_name)
SEU_tot_var <- SEU_RNA@assays[["PROM"]]@meta.features[,c(2,4,7)]
DefaultAssay(SEU_RNA) <- "PROM"
for (i in levels(SEU_RNA@meta.data[["RNA_snn_res.0.8"]])){
    sub <- subset(SEU_RNA, idents = i)
    sub <- FindVariableFeatures(sub)
    sub <- FindVariableFeatures(sub, selection.method = "mvp")
    table_var <- sub@assays[["PROM"]]@meta.features[,c(2,4,7)]
    table_expr <- sub@assays[["PROM"]]@meta.features[,c(1,6)]
    table_ratio <- table_var/SEU_tot_var
    colnames(table_var) <- paste0("CL",i,"_", colnames(table_var))
    colnames(table_expr) <- paste0("CL",i,"_", colnames(table_expr))
    colnames(table_ratio) <- paste0("CL",i,"_", colnames(table_ratio))
    variance_table <- cbind(variance_table, table_var)
    mean_table <- cbind(mean_table, table_expr)
    variance_ratio_table <- cbind(variance_ratio_table, table_ratio)
}

variance_table_prom <- variance_table
mean_table_prom <- mean_table
variance_ratio_table_prom <- variance_ratio_table


variance_table <- data.frame(active_prom_name)
variance_ratio_table <- data.frame(active_prom_name)
mean_table <- data.frame(active_prom_name)
SEU_tot_var <- SEU_RNA@assays[["RNA"]]@meta.features[,c(2,4,7)]
DefaultAssay(SEU_RNA)<- "RNA"
for (i in levels(SEU_RNA@meta.data[["RNA_snn_res.0.8"]])){
  sub <- subset(SEU_RNA, idents = i)
  sub <- FindVariableFeatures(sub)
  sub <- FindVariableFeatures(sub, selection.method = "mvp")
  table_var <- sub@assays[["RNA"]]@meta.features[,c(2,4,7)]
  table_expr <- sub@assays[["RNA"]]@meta.features[,c(1,6)]
  table_ratio <- table_var/SEU_tot_var
  colnames(table_var) <- paste0("CL",i,"_", colnames(table_var))
  colnames(table_expr) <- paste0("CL",i,"_", colnames(table_expr))
  colnames(table_ratio) <- paste0("CL",i,"_", colnames(table_ratio))
  variance_table <- cbind(variance_table, table_var)
  mean_table <- cbind(mean_table, table_expr)
  variance_ratio_table <- cbind(variance_ratio_table, table_ratio)
}

variance_table_expr <- variance_table
mean_table_expr <- mean_table
variance_ratio_table_expr <- variance_ratio_table

exon_prom_name <- rownames(SEU_RNA@assays[["EXON"]])
variance_table <- data.frame(exon_prom_name)
variance_ratio_table <- data.frame(exon_prom_name)
mean_table <- data.frame(exon_prom_name)
SEU_tot_var <- SEU_RNA@assays[["EXON"]]@meta.features[,c(2,4,7)]
DefaultAssay(SEU_RNA)<- "EXON"
for (i in levels(SEU_RNA@meta.data[["RNA_snn_res.0.8"]])){
  sub <- subset(SEU_RNA, idents = i)
  sub <- FindVariableFeatures(sub)
  sub <- FindVariableFeatures(sub, selection.method = "mvp")
  table_var <- sub@assays[["EXON"]]@meta.features[,c(2,4,7)]
  table_expr <- sub@assays[["EXON"]]@meta.features[,c(1,6)]
  table_ratio <- table_var/SEU_tot_var
  colnames(table_var) <- paste0("CL",i,"_", colnames(table_var))
  colnames(table_expr) <- paste0("CL",i,"_", colnames(table_expr))
  colnames(table_ratio) <- paste0("CL",i,"_", colnames(table_ratio))
  variance_table <- cbind(variance_table, table_var)
  mean_table <- cbind(mean_table, table_expr)
  variance_ratio_table <- cbind(variance_ratio_table, table_ratio)
}

variance_table_exon <- variance_table
mean_table_exon <- mean_table
variance_ratio_table_exon <- variance_ratio_table


variance_table <- data.frame(active_prom_name)
variance_ratio_table <- data.frame(active_prom_name)
mean_table <- data.frame(active_prom_name)
SEU_tot_var <- SEU_RNA@assays[["ENHD"]]@meta.features[,c(2,4,7)]
DefaultAssay(SEU_RNA)<- "ENHD"
for (i in levels(SEU_RNA@meta.data[["RNA_snn_res.0.8"]])){
  sub <- subset(SEU_RNA, idents = i)
  sub <- FindVariableFeatures(sub)
  sub <- FindVariableFeatures(sub, selection.method = "mvp")
  table_var <- sub@assays[["ENHD"]]@meta.features[,c(2,4,7)]
  table_expr <- sub@assays[["ENHD"]]@meta.features[,c(1,6)]
  table_ratio <- table_var/SEU_tot_var
  colnames(table_var) <- paste0("CL",i,"_", colnames(table_var))
  colnames(table_expr) <- paste0("CL",i,"_", colnames(table_expr))
  colnames(table_ratio) <- paste0("CL",i,"_", colnames(table_ratio))
  variance_table <- cbind(variance_table, table_var)
  mean_table <- cbind(mean_table, table_expr)
  variance_ratio_table <- cbind(variance_ratio_table, table_ratio)
}

variance_table_enhd <- variance_table
mean_table_enhd <- mean_table
variance_ratio_table_enhd <- variance_ratio_table


######## Functions ########

cluster_scatter_plots(CL, name = paste0("CL", i-1),highlight_low_var = markers_genes_cl, log = TRUE, folder_name = "ENHANCER")

cluster_scatter_plots <- function (CL_0,name,  log = TRUE, division = "mean", highlight_low_var = NULL, folder_name = "PROM"){
  
  if (!(dir.exists(paste0("../TMPResults/IMAGES/",folder_name)))){
    dir.create(paste0("../TMPResults/IMAGES/",folder_name))
  }
  
  if(log){
    if (division == "mean"){
      #png(paste0("../TMPResults/IMAGES/",name,"_gene_scatter_plot.png"),width = 1080, height = 1080)
      plot <-ggplot(as.data.frame(log(CL_0)), aes(x=Activity, y = Expression)) + geom_point() + geom_hline(yintercept = colMeans(log(CL_0)[log(CL_0)[,2] > -Inf,])[2], linetype="dashed", color = "green", size=2) + 
        geom_vline(xintercept = colMeans(log(CL_0)[log(CL_0)[,1] > -Inf,])[1], linetype="dashed", color = "green", size=2)+ ggtitle(name) +
        theme(plot.title = element_text(hjust = 0.5, size = 30), axis.title = element_text(hjust = 0.5, size = 30))+
        geom_point(data=as.data.frame(log(CL_0))[highlight_low_var, ], aes(x=Activity, y = Expression), colour="red", size=2) #+ xlim(-150,3)
      #geom_point(data=as.data.frame(log(CL_0))[high_cor_genes, ], aes(x=Activity, y = Expression), colour="red", size=2)
      #dev.off()
      ggsave(path = paste0("../TMPResults/IMAGES/", folder_name,"/"), filename = paste0(name,"_gene_scatter_plot.pdf"), width = 1080, height = 1080, units= "px",scale = 3.5)
    }
    
  }
  if(!log){
    if(division == "mean"){
      #png(paste0("../TMPResults/IMAGES/",name,"_gene_scatter_plot.png"),width = 1080, height = 1080)
      ggplot(as.data.frame((CL_0)), aes(x=Activity, y = Expression)) + geom_point() + geom_hline(yintercept = colMeans((CL_0))[2], linetype="dashed", color = "green", size=2) + 
        geom_vline(xintercept = colMeans((CL_0))[1], linetype="dashed", color = "green", size=2)+ ggtitle(name) +
        theme(plot.title = element_text(hjust = 0.5, size = 30), axis.title = element_text(hjust = 0.5, size = 30))+
        geom_point(data=as.data.frame((CL_0))[highlight_low_var, ], aes(x=Activity, y = Expression), colour="red", size=2)
      #geom_point(data=as.data.frame((CL_0))[high_cor_genes, ], aes(x=Activity, y = Expression), colour="red", size=2)
      ggsave(path = paste0("../TMPResults/IMAGES/", folder_name,"/"), filename = paste0(name,"_gene_scatter_plot.pdf"), width = 1080, height = 1080, units= "px",scale = 3.5)
      #dev.off()
      
    }
  }
}
cluster_scatter_plots(CL, name = paste0("CL", i-1),highlight_low_var = markers_genes_cl, log = TRUE, folder_name = "ENHANCER")

number_correlated_genes <- function(CL_0,name, log = TRUE){
  
  table <- data.frame(H_H=0, H_L= 0, L_H = 0, L_L =0, row.names= name)
  if (log){
    table$H_H = sum(log(CL_0)[,1] >= colMeans(log(CL_0)[log(CL_0)[,1] > -Inf,])[1] & log(CL_0)[,2] >= colMeans(log(CL_0)[log(CL_0)[,2] > -Inf,])[2])
    table$H_L = sum(log(CL_0)[,1] >= colMeans(log(CL_0)[log(CL_0)[,1] > -Inf,])[1] & log(CL_0)[,2] < colMeans(log(CL_0)[log(CL_0)[,2] > -Inf,])[2])
    table$L_H = sum(log(CL_0)[,1] < colMeans(log(CL_0)[log(CL_0)[,1] > -Inf,])[1] & log(CL_0)[,2] >= colMeans(log(CL_0)[log(CL_0)[,2] > -Inf,])[2])
    table$L_L = sum(log(CL_0)[,1] < colMeans(log(CL_0)[log(CL_0)[,1] > -Inf,])[1] & log(CL_0)[,2] < colMeans(log(CL_0)[log(CL_0)[,2] > -Inf,])[2])
  }
  if(!log){
    table$H_H = sum((CL_0)[,1] >= colMeans((CL_0)[(CL_0)[,1] > -Inf,])[1] & (CL_0)[,2] >= colMeans((CL_0)[(CL_0)[,2] > -Inf,])[2])
    table$H_L = sum((CL_0)[,1] >= colMeans((CL_0)[(CL_0)[,1] > -Inf,])[1] & (CL_0)[,2] < colMeans((CL_0)[(CL_0)[,2] > -Inf,])[2])
    table$L_H = sum((CL_0)[,1] < colMeans((CL_0)[(CL_0)[,1] > -Inf,])[1] & (CL_0)[,2] >= colMeans((CL_0)[(CL_0)[,2] > -Inf,])[2])
    table$L_L = sum((CL_0)[,1] < colMeans((CL_0)[(CL_0)[,1] > -Inf,])[1] & (CL_0)[,2] < colMeans((CL_0)[(CL_0)[,2] > -Inf,])[2])
  }
  return(table)
}

number_correlated_genes_lines <- function(CL_0,name, limex, limacc ,log = TRUE){
  
  table <- data.frame(H_H=0, H_L= 0, L_H = 0, L_L =0, row.names= name)
  if (log){
    table$H_H = sum(log(CL_0)[,1] >= limacc & log(CL_0)[,2] >= limex)
    table$H_L = sum(log(CL_0)[,1] >= limacc & log(CL_0)[,2] < limex)
    table$L_H = sum(log(CL_0)[,1] < limacc & log(CL_0)[,2] >= limex)
    table$L_L = sum(log(CL_0)[,1] < limacc & log(CL_0)[,2] < limex)
  }
  if(!log){
    table$H_H = sum((CL_0)[,1] >= limacc & (CL_0)[,2] >= limex)
    table$H_L = sum((CL_0)[,1] >= limacc & (CL_0)[,2] < limex)
    table$L_H = sum((CL_0)[,1] < limacc& (CL_0)[,2] >= limex)
    table$L_L = sum((CL_0)[,1] < limacc & (CL_0)[,2] < limex)
  }
  return(table)
}




incoherence_factor <- function(CL_0,name ,log = FALSE){
  
  CL_0[CL_0 > 0] <-1
  table <- data.frame(H_H=0, H_L= 0, L_H = 0, L_L =0, row.names= name)
  if (log){
    table$H_H = sum(log(CL_0)[,1] > 0 & log(CL_0)[,2] > 0)
    table$H_L = sum(log(CL_0)[,1] > 0 & log(CL_0)[,2] <= 0)
    table$L_H = sum(log(CL_0)[,1] <= 0 & log(CL_0)[,2] > 0)
    table$L_L = sum(log(CL_0)[,1] <= 0 & log(CL_0)[,2] <= 0)
  }
  if(!log){
    table$H_H = sum((CL_0)[,1] > 0 & (CL_0)[,2] > 0)
    table$H_L = sum((CL_0)[,1] > 0 & (CL_0)[,2] <= 0)
    table$L_H = sum((CL_0)[,1] <= 0& (CL_0)[,2] > 0)
    table$L_L = sum((CL_0)[,1] <= 0 & (CL_0)[,2] <= 0)
  }
  return(table)
}



############ PROMOTER #######



key = "_mvp.mean"
key = "_vst.mean"
Sce_mean_Act = mean_table_prom[,paste0("CL", 0, key)]
Sce_mean_Expr = mean_table_expr[,paste0("CL", 0, key)]
for (i in (2:((length(colnames(mean_table_expr))-1)/2))){
  
  Sce_mean_Act = cbind(Sce_mean_Act, mean_table_prom[,paste0("CL", i-1, key)])
  Sce_mean_Expr = cbind (Sce_mean_Expr, mean_table_expr[,paste0("CL", i-1, key)])
  
}
rownames(Sce_mean_Act) <- mean_table_prom$active_prom_name
rownames(Sce_mean_Expr) <- mean_table_expr$active_prom_name


#### correlation ####
correlation_table_final <- info_table[,1:2]
correlation_table_final[,2] <- NULL

s <- sapply(correlation_table_final$gene, function(x) cor.test(Sce_mean_Act[x,], Sce_mean_Expr[x,])[["estimate"]][["cor"]] )
correlation_table_final$PcorAvg_m <- s
s <- sapply(correlation_table_final$gene, function(x) cor.test(Sce_mean_Act[x,], Sce_mean_Expr[x,])[["p.value"]] )
correlation_table_final$PcorAvg_m_pvalue <- s

write.table(correlation_table_final,file = "../TMPResults/Tables/correlation_table_promoter.csv")

SEU_tot_var <- SEU_RNA@assays[["RNA"]]@meta.features[,c(2,4,7)]
low_variance_genes <- rownames(top_n(SEU_tot_var, -1000, vst.variance.standardized))
low_variance_genes <- rownames(top_n(SEU_tot_var, -1000, mvp.dispersion))
prova <- correlation_table_final[correlation_table_final$gene %in% low_variance_genes,]
dim(prova[prova$PcorAvg_m >0.5,])[1]/dim(prova)[1]
dim(prova[correlation_table_final$PcorAvg_m >0.5,])[1]/dim(correlation_table_final)[1]


markers_genes_cl <- expr_markers %>% filter(cluster == i-1) %>% filter(avg_log2FC > 0)
markers_genes_cl <- markers_genes_cl$gene
prova <- correlation_table_final[correlation_table_final$gene %in% markers_genes_cl,]
dim(prova[prova$PcorAvg_m >0.5,])[1]/dim(prova)[1]
dim(prova[correlation_table_final$PcorAvg_m >0.5,])[1]/dim(correlation_table_final)[1]


#### Plots ####
CL <- cbind(mean_table_prom[,paste0("CL",0 , key)], mean_table_expr[,paste0("CL", 0, key)])
colnames(CL) <- c("Activity", "Expression")
row <- number_correlated_genes(CL, name = paste0("CL_0"))
number_genes_table_prom <- rbind(row)
for (i in (1:((length(colnames(mean_table_expr))-1)/2))){
  
  CL <- cbind(mean_table_prom[,paste0("CL", i-1, key)], mean_table_expr[,paste0("CL", i-1, key)])
  colnames(CL) <- c("Activity", "Expression")
  rownames(CL) <- rownames(mean_table_expr)
  
  low_variance_genes <- (top_n(variance_ratio_table_expr[rownames(CL[CL[,"Expression"] > 0,]),c("active_prom_name",paste0("CL", i-1, "_vst.variance.standardized"))], -1000)["active_prom_name"])$active_prom_name
  low_variance_genes_acc <- (top_n(variance_ratio_table_prom[rownames(CL[CL[,"Activity"] > 0,]),c("active_prom_name",paste0("CL", i-1, "_vst.variance.standardized"))], -1000)["active_prom_name"])$active_prom_name
  markers_genes_cl <- expr_markers %>% filter(cluster == i-1) %>% filter(avg_log2FC > 0)
  markers_genes_cl <- markers_genes_cl$gene
  cluster_scatter_plots(CL, name = paste0("CL", i-1),highlight_low_var = markers_genes_cl, log = TRUE, folder_name = "PROM")
  if(i>1){
    row <- number_correlated_genes(CL, name = paste0("CL_", i-1))
    number_genes_table_prom <- rbind(number_genes_table_prom,row)
  }
}

#### Coherence ####
CL <- cbind(mean_table_prom[,paste0("CL",0 , key)], mean_table_expr[,paste0("CL", 0, key)])
colnames(CL) <- c("Activity", "Expression")
row <- incoherence_factor(CL, name = paste0("CL_0"))
row <- (row/sum(row))*100
incoherence_factor_prom <- rbind(row)

for (i in (1:((length(colnames(mean_table_expr))-1)/2))){
  CL <- cbind(mean_table_prom[,paste0("CL", i-1, key)], mean_table_expr[,paste0("CL", i-1, key)])
  colnames(CL) <- c("Activity", "Expression")
  if(i>1){
    row <- incoherence_factor(CL, name = paste0("CL_", i-1), log = FALSE)
    row <- (row/sum(row))*100
    incoherence_factor_prom <- rbind(incoherence_factor_prom, row)
  }
  print(paste0("CL", i-1, key))
}

colMeans(incoherence_factor_prom)
incoherence_factor_prom <- rbind(incoherence_factor_prom, colMeans(incoherence_factor_prom))
rownames(incoherence_factor_prom)[length(rownames(incoherence_factor_prom))] <- "Average"
write.table(incoherence_factor_prom,file = "../TMPResults/Tables/incoherence_factor_prom.csv")

############## EXON #########

key = "_mvp.mean"
key = "_vst.mean"
Sce_mean_Act = mean_table_exon[,paste0("CL", 0, key)]
Sce_mean_Expr = mean_table_expr[exon_prom_name,paste0("CL", 0, key)]
for (i in (2:((length(colnames(mean_table_expr))-1)/2))){
  
  Sce_mean_Act = cbind(Sce_mean_Act, mean_table_exon[,paste0("CL", i-1, key)])
  Sce_mean_Expr = cbind (Sce_mean_Expr, mean_table_expr[exon_prom_name,paste0("CL", i-1, key)])
  
}
rownames(Sce_mean_Act) <- mean_table_exon$exon_prom_name
rownames(Sce_mean_Expr) <- exon_prom_name

#### correlation ####
correlation_table_final <- info_table[,1:2]
correlation_table_final[,2] <- NULL
correlation_table_final <- correlation_table[correlation_table$gene %in% exon_prom_name, ]


s <- sapply(correlation_table_final$gene, function(x) cor.test(Sce_mean_Act[x,], Sce_mean_Expr[x,])[["estimate"]][["cor"]] )
correlation_table_final$PcorAvg_m <- s
s <- sapply(correlation_table_final$gene, function(x) cor.test(Sce_mean_Act[x,], Sce_mean_Expr[x,])[["p.value"]] )
correlation_table_final$PcorAvg_m_pvalue <- s

write.table(correlation_table_final,file = "../TMPResults/Tables/correlation_table_exon.csv")

SEU_tot_var <- SEU_RNA@assays[["RNA"]]@meta.features[,c(2,4,7)]
low_variance_genes <- rownames(top_n(SEU_tot_var, -1000, vst.variance.standardized))
low_variance_genes <- rownames(top_n(SEU_tot_var, -1000, mvp.dispersion))
prova <- correlation_table_final[correlation_table_final$gene %in% low_variance_genes,]
dim(prova[prova$PcorAvg_m >0.5,])[1]/dim(prova)[1]
dim(prova[correlation_table_final$PcorAvg_m >0.5,])[1]/dim(correlation_table_final)[1]

markers_genes_cl <- expr_markers %>% filter(cluster == i-1) %>% filter(avg_log2FC > 0)
markers_genes_cl <- markers_genes_cl$gene
prova <- correlation_table_final[correlation_table_final$gene %in% markers_genes_cl,]
dim(prova[prova$PcorAvg_m >0.5,])[1]/dim(prova)[1]
dim(prova[correlation_table_final$PcorAvg_m >0.5,])[1]/dim(correlation_table_final)[1]


#### Plots ####
CL <- cbind(mean_table_exon[,paste0("CL",0 , key)], mean_table_expr[exon_prom_name,paste0("CL", 0, key)])
colnames(CL) <- c("Activity", "Expression")
row <- number_correlated_genes(CL, name = paste0("CL_0"))
number_genes_table_prom <- rbind(row)
for (i in (1:((length(colnames(mean_table_expr))-1)/2))){
  
  CL <- cbind(mean_table_exon[,paste0("CL", i-1, key)], mean_table_expr[exon_prom_name,paste0("CL", i-1, key)])
  colnames(CL) <- c("Activity", "Expression")
  rownames(CL) <- rownames(mean_table_exon)
  
  #low_variance_genes <- (top_n(variance_ratio_table_expr[rownames(CL[CL[,"Expression"] > 0,]),c("active_prom_name",paste0("CL", i-1, "_vst.variance.standardized"))], -1000)["active_prom_name"])$active_prom_name
  #low_variance_genes_acc <- (top_n(variance_ratio_table_prom[rownames(CL[CL[,"Activity"] > 0,]),c("active_prom_name",paste0("CL", i-1, "_vst.variance.standardized"))], -1000)["active_prom_name"])$active_prom_name
  markers_genes_cl <- expr_markers %>% filter(cluster == i-1) %>% filter(avg_log2FC > 0)
  markers_genes_cl <- markers_genes_cl$gene
  markers_genes_cl <- markers_genes_cl[markers_genes_cl %in% exon_prom_name]
  cluster_scatter_plots(CL, name = paste0("CL", i-1),highlight_low_var = markers_genes_cl, log = TRUE, folder_name = "EXON")
  if(i>1){
    row <- number_correlated_genes(CL, name = paste0("CL_", i-1))
    number_genes_table_prom <- rbind(number_genes_table_prom,row)
  }
}

#### Coherence ####
CL <- cbind(mean_table_exon[,paste0("CL",0 , key)], mean_table_expr[exon_prom_name,paste0("CL", 0, key)])
colnames(CL) <- c("Activity", "Expression")
row <- incoherence_factor(CL, name = paste0("CL_0"))
row <- (row/sum(row))*100
incoherence_factor_exon <- rbind(row)

for (i in (1:((length(colnames(mean_table_expr))-1)/2))){
  CL <- cbind(mean_table_exon[,paste0("CL", i-1, key)], mean_table_expr[exon_prom_name,paste0("CL", i-1, key)])
  colnames(CL) <- c("Activity", "Expression")
  if(i>1){
    row <- incoherence_factor(CL, name = paste0("CL_", i-1), log = FALSE)
    row <- (row/sum(row))*100
    incoherence_factor_exon <- rbind(incoherence_factor_exon, row)
  }
  print(paste0("CL", i-1, key))
}

colMeans(incoherence_factor_exon)

incoherence_factor_exon <- rbind(incoherence_factor_exon, colMeans(incoherence_factor_exon))
rownames(incoherence_factor_exon)[length(rownames(incoherence_factor_exon))] <- "Average"
write.table(incoherence_factor_exon,file = "../TMPResults/Tables/incoherence_factor_exon.csv")

############## ENHANCER #########

key = "_mvp.mean"
key = "_vst.mean"
Sce_mean_Act = mean_table_enhd[,paste0("CL", 0, key)]
Sce_mean_Expr = mean_table_expr[,paste0("CL", 0, key)]
for (i in (2:((length(colnames(mean_table_expr))-1)/2))){
  
  Sce_mean_Act = cbind(Sce_mean_Act, mean_table_enhd[,paste0("CL", i-1, key)])
  Sce_mean_Expr = cbind (Sce_mean_Expr, mean_table_expr[,paste0("CL", i-1, key)])
  
}
rownames(Sce_mean_Act) <- mean_table_enhd$active_prom_name
rownames(Sce_mean_Expr) <- mean_table_expr$active_prom_name

#### correlation ####
correlation_table_final <- info_table[,1:2]
correlation_table_final[,2] <- NULL

s <- sapply(correlation_table_final$gene, function(x) cor.test(Sce_mean_Act[x,], Sce_mean_Expr[x,])[["estimate"]][["cor"]] )
correlation_table_final$PcorAvg_m <- s
s <- sapply(correlation_table_final$gene, function(x) cor.test(Sce_mean_Act[x,], Sce_mean_Expr[x,])[["p.value"]])
correlation_table_final$PcorAvg_m_pvalue <- s

write.table(correlation_table_final,file = "../TMPResults/Tables/correlation_table_enhancer.csv")

SEU_tot_var <- SEU_RNA@assays[["RNA"]]@meta.features[,c(2,4,7)]
low_variance_genes <- rownames(top_n(SEU_tot_var, -1000, vst.variance.standardized))
low_variance_genes <- rownames(top_n(SEU_tot_var, -1000, mvp.dispersion))
prova <- correlation_table_final[correlation_table_final$gene %in% low_variance_genes,]
dim(prova[prova$PcorAvg_m >0.5,])[1]/dim(prova)[1]
dim(prova[correlation_table_final$PcorAvg_m >0.5,])[1]/dim(correlation_table_final)[1]

markers_genes_cl <- expr_markers %>% filter(cluster == i-1) %>% filter(avg_log2FC > 0)
markers_genes_cl <- markers_genes_cl$gene
prova <- correlation_table_final[correlation_table_final$gene %in% markers_genes_cl,]
dim(prova[prova$PcorAvg_m >0.5,])[1]/dim(prova)[1]
dim(prova[correlation_table_final$PcorAvg_m >0.5,])[1]/dim(correlation_table_final)[1]


#### Plots ####
CL <- cbind(mean_table_enhd[,paste0("CL",0 , key)], mean_table_expr[,paste0("CL", 0, key)])
colnames(CL) <- c("Activity", "Expression")
row <- number_correlated_genes(CL, name = paste0("CL_0"))
number_genes_table_prom <- rbind(row)
for (i in (1:((length(colnames(mean_table_expr))-1)/2))){
  
  CL <- cbind(mean_table_enhd[,paste0("CL", i-1, key)], mean_table_expr[,paste0("CL", i-1, key)])
  colnames(CL) <- c("Activity", "Expression")
  rownames(CL) <- rownames(mean_table_expr)
  
  low_variance_genes <- (top_n(variance_ratio_table_expr[rownames(CL[CL[,"Expression"] > 0,]),c("active_prom_name",paste0("CL", i-1, "_vst.variance.standardized"))], -1000)["active_prom_name"])$active_prom_name
  low_variance_genes_acc <- (top_n(variance_ratio_table_prom[rownames(CL[CL[,"Activity"] > 0,]),c("active_prom_name",paste0("CL", i-1, "_vst.variance.standardized"))], -1000)["active_prom_name"])$active_prom_name
  markers_genes_cl <- expr_markers %>% filter(cluster == i-1) %>% filter(avg_log2FC > 0)
  markers_genes_cl <- markers_genes_cl$gene
  cluster_scatter_plots(CL, name = paste0("CL", i-1),highlight_low_var = markers_genes_cl, log = TRUE, folder_name = "ENHANCER")
  if(i>1){
    row <- number_correlated_genes(CL, name = paste0("CL_", i-1))
    number_genes_table_prom <- rbind(number_genes_table_prom,row)
  }
}

#### Coherence ####
CL <- cbind(mean_table_enhd[,paste0("CL",0 , key)], mean_table_expr[,paste0("CL", 0, key)])
colnames(CL) <- c("Activity", "Expression")
row <- incoherence_factor(CL, name = paste0("CL_0"))
row <- (row/sum(row))*100
incoherence_factor_enhd <- rbind(row)

for (i in (1:((length(colnames(mean_table_expr))-1)/2))){
  CL <- cbind(mean_table_enhd[,paste0("CL", i-1, key)], mean_table_expr[,paste0("CL", i-1, key)])
  colnames(CL) <- c("Activity", "Expression")
  if(i>1){
    row <- incoherence_factor(CL, name = paste0("CL_", i-1), log = FALSE)
    row <- (row/sum(row))*100
    incoherence_factor_enhd <- rbind(incoherence_factor_enhd, row)
  }
  print(paste0("CL", i-1, key))
}

colMeans(incoherence_factor_enhd)

incoherence_factor_enhd <- rbind(incoherence_factor_enhd, colMeans(incoherence_factor_enhd))
rownames(incoherence_factor_enhd)[length(rownames(incoherence_factor_enhd))] <- "Average"
write.table(incoherence_factor_enhd,file = "../TMPResults/Tables/incoherence_factor_enhd.csv")



####################################





























temp_bin <- refseq_first_gene_matrix[active_prom_name,]
temp_bin@x[temp_bin@x>0] <- 1
X <- temp_bin

temp_bin <- RNA_matrix[active_prom_name,]
temp_bin@x[temp_bin@x > 0] <- 1
Y <- temp_bin

info_table$JPE <- 0
info_table_temp <- info_table
info_table_temp <- info_table_temp[active_prom_name,]

info_table_temp <- info_table_temp %>% rowwise() %>% mutate(JPE = jaccard(X[gene,],Y[gene,]))
write.table(info_table_temp, file = "info_table_temp.csv")
info_table_temp <-read.table("info_table_temp.csv")

s <- sapply(info_table_temp$gene, function(x) RI(X[x,],Y[x,]))
write.table(s, file = "RI.csv")

rownames(info_table_temp) <- info_table_temp$gene
info_table_temp$RPE <- 0
for (g in active_prom_name){
  info_table_temp[info_table_temp$gene == g,]$RPE <- RI(X[g,],Y[g,])
  print(g)
  
}

#####
low_variance_genes <- (top_n(variance_ratio_table_expr[rownames(CL[CL[,"Expression"] > 0,]),c("active_prom_name",paste0("CL", i-1, "_vst.variance.standardized"))], -1000)["active_prom_name"])$active_prom_name
prova <- correlation_table[correlation_table$gene %in% low_variance_genes,]
rm(prova)



SEU_RNA@assays[["PROM"]]@var.features

pie1 <- data.frame()

#cor <- cor.test(AggProm[["RNA"]][1,], AggProm[["PROM"]][1,])[["p.value"]]
################
s <- sapply(info_table_temp$gene, function(x) cor.test(AggProm[["RNA"]][x,], AggProm[["PROM"]][x,])[["estimate"]][["cor"]] )
info_table_temp$PcorAgg <- s
s <- sapply(info_table_temp$gene, function(x) cor.test(AggProm[["RNA"]][x,], AggProm[["PROM"]][x,])[["p.value"]] )
info_table_temp$PcorAgg_pvalue <- s

s <- sapply(info_table_temp$gene, function(x) cor.test(AvgProm[["RNA"]][x,], AvgProm[["PROM"]][x,])[["estimate"]][["cor"]] )
info_table_temp$PcorAvg <- s
s <- sapply(info_table_temp$gene, function(x) cor.test(AvgProm[["RNA"]][x,], AvgProm[["PROM"]][x,])[["p.value"]] )
info_table_temp$PcorAvg_pvalue <- s


SCE <- as.SingleCellExperiment(SEU_RNA)
SCE_MEAN <- aggregateData(SCE,  assay = "logcounts" , by = "RNA_snn_res.0.8", fun = "mean")

correlation_table <- info_table_temp[,1:2]
correlation_table[,2] <- NULL

AggProm <- AggregateExpression(SEU_RNA, group.by = "RNA_snn_res.0.8", features = active_prom_name)
AvgProm <- AverageExpression(SEU_RNA, group.by = "RNA_snn_res.0.8", features = active_prom_name)

s <- sapply(correlation_table$gene, function(x) cor.test(AggProm[["RNA"]][x,], AggProm[["PROM"]][x,])[["estimate"]][["cor"]] )
correlation_table$PcorAgg <- s
s <- sapply(correlation_table$gene, function(x) cor.test(AggProm[["RNA"]][x,], AggProm[["PROM"]][x,])[["p.value"]] )
correlation_table$PcorAgg_pvalue <- s

s <- sapply(correlation_table$gene, function(x) cor.test(AvgProm[["RNA"]][x,], AvgProm[["PROM"]][x,])[["estimate"]][["cor"]] )
correlation_table$PcorAvg <- s
s <- sapply(correlation_table$gene, function(x) cor.test(AvgProm[["RNA"]][x,], AvgProm[["PROM"]][x,])[["p.value"]] )
correlation_table$PcorAvg_pvalue <- s

###### AVG - COUNTS ######
DefaultAssay(SEU_RNA) <- "PROM"
SCE <- as.SingleCellExperiment(SEU_RNA)
Sce_mean_Act <- aggregateData(SCE,  assay = "counts" , by = "RNA_snn_res.0.8", fun = "mean")
Sce_mean_Act <- Sce_mean_Act@assays@data@listData[[1]]
Sce_mean_Act <- Sce_mean_Act[active_prom_name,]

DefaultAssay(SEU_RNA) <- "RNA"
SCE <- as.SingleCellExperiment(SEU_RNA)
Sce_mean_Expr <- aggregateData(SCE,  assay = "counts" , by = "RNA_snn_res.0.8", fun = "mean")
Sce_mean_Expr <- Sce_mean_Expr@assays@data@listData[[1]]
Sce_mean_Expr <- Sce_mean_Expr[active_prom_name,]



s <- sapply(correlation_table$gene, function(x) cor.test(Sce_mean_Act[x,], Sce_mean_Expr[x,])[["estimate"]][["cor"]] )
correlation_table$PcorAvg_m <- s
s <- sapply(correlation_table$gene, function(x) cor.test(Sce_mean_Act[x,], Sce_mean_Expr[x,])[["p.value"]] )
correlation_table$PcorAvg_m_pvalue <- s



CL <- cbind(Sce_mean_Act[,1], Sce_mean_Expr[,1])
colnames(CL) <- c("Activity", "Expression")
row <- number_correlated_genes(CL, name = paste0("CL_0"))
number_genes_table_prom <- rbind(row)
for (i in (1:length(colnames(Sce_mean_Expr)))){
  CL <- cbind(Sce_mean_Act[,i], Sce_mean_Expr[,i])
  colnames(CL) <- c("Activity", "Expression")
  cluster_scatter_plots(CL, name = paste0("CL_", i-1))
  if(i>1){
    row <- number_correlated_genes(CL, name = paste0("CL_", i-1))
    number_genes_table_prom <- rbind(number_genes_table_prom,row)
  }
}

CL <- cbind(Sce_mean_Act[,1], Sce_mean_Expr[,1])
colnames(CL) <- c("Activity", "Expression")
row <- incoherence_factor(CL, name = paste0("CL_0"))
row <- (row/sum(row))*100
incoherence_factor_prom <- rbind(row)
for (i in (1:length(colnames(Sce_mean_Expr)))){
  CL <- cbind(Sce_mean_Act[,i], Sce_mean_Expr[,i])
  colnames(CL) <- c("Activity", "Expression")
  if(i>1){
    row <- incoherence_factor(CL, name = paste0("CL_", i-1), log = FALSE)
    row <- (row/sum(row))*100
    incoherence_factor_prom <- rbind(incoherence_factor_prom, row)
  }
}


###### AVG - LOGCOUNTS ######
DefaultAssay(SEU_RNA) <- "PROM"
SCE <- as.SingleCellExperiment(SEU_RNA)
Sce_mean_Act <- aggregateData(SCE,  assay = "logcounts" , by = "RNA_snn_res.0.8", fun = "mean")
Sce_mean_Act <- Sce_mean_Act@assays@data@listData[[1]]
Sce_mean_Act <- Sce_mean_Act[active_prom_name,]

DefaultAssay(SEU_RNA) <- "RNA"
SCE <- as.SingleCellExperiment(SEU_RNA)
Sce_mean_Expr <- aggregateData(SCE,  assay = "logcounts" , by = "RNA_snn_res.0.8", fun = "mean")
Sce_mean_Expr <- Sce_mean_Expr@assays@data@listData[[1]]
Sce_mean_Expr <- Sce_mean_Expr[active_prom_name,]

s <- sapply(correlation_table$gene, function(x) cor.test(Sce_mean_Act[x,], Sce_mean_Expr[x,])[["estimate"]][["cor"]] )
correlation_table$PcorAvg_m_log <- s
s <- sapply(correlation_table$gene, function(x) cor.test(Sce_mean_Act[x,], Sce_mean_Expr[x,])[["p.value"]] )
correlation_table$PcorAvg_m_log_pvalue <- s



CL <- cbind(Sce_mean_Act[,1], Sce_mean_Expr[,1])
colnames(CL) <- c("Activity", "Expression")
row <- number_correlated_genes(CL, name = paste0("CL_0"))
number_genes_table_prom <- rbind(row)


for (i in (1:length(colnames(Sce_mean_Expr)))){
  CL <- cbind(Sce_mean_Act[,i], Sce_mean_Expr[,i])
  colnames(CL) <- c("Activity", "Expression")
   
  cluster_scatter_plots(CL, name = paste0("CL_", i-1))
  if(i>1){
    row <- number_correlated_genes(CL, name = paste0("CL_", i-1))
    number_genes_table_prom <- rbind(number_genes_table_prom,row)
  }
}


###### MEDIAN - COUNTS ######
DefaultAssay(SEU_RNA) <- "PROM"
SCE <- as.SingleCellExperiment(SEU_RNA)
Sce_mean_Act <- aggregateData(SCE,  assay = "counts" , by = "RNA_snn_res.0.8", fun = "median")
Sce_mean_Act <- Sce_mean_Act@assays@data@listData[[1]]
Sce_mean_Act <- Sce_mean_Act[active_prom_name,]

DefaultAssay(SEU_RNA) <- "RNA"
SCE <- as.SingleCellExperiment(SEU_RNA)
Sce_mean_Expr <- aggregateData(SCE,  assay = "counts" , by = "RNA_snn_res.0.8", fun = "median")
Sce_mean_Expr <- Sce_mean_Expr@assays@data@listData[[1]]
Sce_mean_Expr <- Sce_mean_Expr[active_prom_name,]

s <- sapply(correlation_table$gene, function(x) cor.test(Sce_mean_Act[x,], Sce_mean_Expr[x,])[["estimate"]][["cor"]] )
correlation_table$PcorMed_m <- s
s <- sapply(correlation_table$gene, function(x) cor.test(Sce_mean_Act[x,], Sce_mean_Expr[x,])[["p.value"]] )
correlation_table$PcorMed_m_pvalue <- s

###### MEDIAN - LOGCOUNTS ######
DefaultAssay(SEU_RNA) <- "PROM"
SCE <- as.SingleCellExperiment(SEU_RNA)
Sce_mean_Act <- aggregateData(SCE,  assay = "logcounts" , by = "RNA_snn_res.0.8", fun = "median")
Sce_mean_Act <- Sce_mean_Act@assays@data@listData[[1]]
Sce_mean_Act <- Sce_mean_Act[active_prom_name,]

DefaultAssay(SEU_RNA) <- "RNA"
SCE <- as.SingleCellExperiment(SEU_RNA)
Sce_mean_Expr <- aggregateData(SCE,  assay = "logcounts" , by = "RNA_snn_res.0.8", fun = "median")
Sce_mean_Expr <- Sce_mean_Expr@assays@data@listData[[1]]
Sce_mean_Expr <- Sce_mean_Expr[active_prom_name,]

s <- sapply(correlation_table$gene, function(x) cor.test(Sce_mean_Act[x,], Sce_mean_Expr[x,])[["estimate"]][["cor"]] )
correlation_table$PcorMed_m_log <- s
s <- sapply(correlation_table$gene, function(x) cor.test(Sce_mean_Act[x,], Sce_mean_Expr[x,])[["p.value"]] )
correlation_table$PcorMed_m_log_pvalue <- s
######################

##################### EXon contribution #############

DefaultAssay(SEU_RNA) <- "EXON"
SCE <- as.SingleCellExperiment(SEU_RNA)
Sce_mean_Act <- aggregateData(SCE,  assay = "counts" , by = "RNA_snn_res.0.8", fun = "mean")
Sce_mean_Act <- Sce_mean_Act@assays@data@listData[[1]]
axon_genes <- rownames(Sce_mean_Act)#[rownames(Sce_mean_Act) %in% active_prom_name]
Sce_mean_Act <- Sce_mean_Act[exon_prom_name,]

DefaultAssay(SEU_RNA) <- "RNA"
SCE <- as.SingleCellExperiment(SEU_RNA)
Sce_mean_Expr <- aggregateData(SCE,  assay = "counts" , by = "RNA_snn_res.0.8", fun = "mean")
Sce_mean_Expr <- Sce_mean_Expr@assays@data@listData[[1]]
Sce_mean_Expr <- Sce_mean_Expr[exon_prom_name,]

exon_prom_name <- rownames(Sce_mean_Act)
exon_no_expr <- exon_prom_name[!(rownames(info_table) %in% rownames(RNA_matrix))]
exon_prom_name <- exon_prom_name[!(exon_prom_name %in% active_no_expr)]

correlation_table_exon <- info_table_temp[,1:2]
correlation_table_exon <- correlation_table_exon[correlation_table_exon$gene %in% exon_prom_name,]
correlation_table_exon[,2] <- NULL


s <- sapply(correlation_table_exon$gene, function(x) cor.test(Sce_mean_Act[x,], Sce_mean_Expr[x,])[["estimate"]][["cor"]] )
correlation_table_exon$PcorAvg_m <- s
s <- sapply(correlation_table_exon$gene, function(x) cor.test(Sce_mean_Act[x,], Sce_mean_Expr[x,])[["p.value"]] )
correlation_table_exon$PcorAvg_m_pvalue <- s



CL <- cbind(Sce_mean_Act[,1], Sce_mean_Expr[,1])
colnames(CL) <- c("Activity", "Expression")
row <- number_correlated_genes(CL, name = paste0("CL_0"))
number_genes_table_exon <- rbind(row)
for (i in (1:length(colnames(Sce_mean_Expr)))){
  CL <- cbind(Sce_mean_Act[,i], Sce_mean_Expr[,i])
  colnames(CL) <- c("Activity", "Expression")
  cluster_scatter_plots(CL, name = paste0("CL_", i-1), log = TRUE)
  if(i>1){
    row <- number_correlated_genes(CL, name = paste0("CL_", i-1))
    number_genes_table_exon <- rbind(number_genes_table_exon, row)
  }
}


me <- mean (log(Sce_mean_Expr)[log(Sce_mean_Expr) > -Inf])
ma <- mean (log(Sce_mean_Act)[log(Sce_mean_Act) > -Inf])
CL <- cbind(Sce_mean_Act[,1], Sce_mean_Expr[,1])
colnames(CL) <- c("Activity", "Expression")
row <- number_correlated_genes_lines(CL,limex = me, limacc = ma, name = paste0("CL_0"), log = TRUE)
number_genes_table_exon <- rbind(row)
for (i in (1:length(colnames(Sce_mean_Expr)))){
  CL <- cbind(Sce_mean_Act[,i], Sce_mean_Expr[,i])
  colnames(CL) <- c("Activity", "Expression")
  cluster_scatter_plots(CL, name = paste0("CL_", i-1), log = TRUE)
  if(i>1){
    row <- number_correlated_genes_lines(CL,limex = me, limacc = ma, name = paste0("CL_", i-1), log = TRUE)
    number_genes_table_exon <- rbind(number_genes_table_exon, row)
    
  }
}



CL <- cbind(Sce_mean_Act[,1], Sce_mean_Expr[,1])
colnames(CL) <- c("Activity", "Expression")
row <- incoherence_factor(CL, name = paste0("CL_0"))
row <- (row/sum(row))*100
incoherence_factor_exon <- rbind(row)
for (i in (1:length(colnames(Sce_mean_Expr)))){
  CL <- cbind(Sce_mean_Act[,i], Sce_mean_Expr[,i])
  colnames(CL) <- c("Activity", "Expression")
  if(i>1){
    row <- incoherence_factor(CL, name = paste0("CL_", i-1), log = FALSE)
    row <- (row/sum(row))*100
    incoherence_factor_exon <- rbind(incoherence_factor_exon, row)
  }
}


##### LOG ######
DefaultAssay(SEU_RNA) <- "EXON"
SCE <- as.SingleCellExperiment(SEU_RNA)
Sce_mean_Act <- aggregateData(SCE,  assay = "logcounts" , by = "RNA_snn_res.0.8", fun = "mean")
Sce_mean_Act <- Sce_mean_Act@assays@data@listData[[1]]
axon_genes <- rownames(Sce_mean_Act)#[rownames(Sce_mean_Act) %in% active_prom_name]
Sce_mean_Act <- Sce_mean_Act[exon_prom_name,]

DefaultAssay(SEU_RNA) <- "RNA"
SCE <- as.SingleCellExperiment(SEU_RNA)
Sce_mean_Expr <- aggregateData(SCE,  assay = "logcounts" , by = "RNA_snn_res.0.8", fun = "mean")
Sce_mean_Expr <- Sce_mean_Expr@assays@data@listData[[1]]
Sce_mean_Expr <- Sce_mean_Expr[exon_prom_name,]

exon_prom_name <- rownames(Sce_mean_Act)
exon_no_expr <- exon_prom_name[!(rownames(info_table) %in% rownames(RNA_matrix))]
exon_prom_name <- exon_prom_name[!(exon_prom_name %in% active_no_expr)]

correlation_table_exon <- info_table_temp[,1:2]
correlation_table_exon <- correlation_table_exon[correlation_table_exon$gene %in% exon_prom_name,]
correlation_table_exon[,2] <- NULL


s <- sapply(correlation_table_exon$gene, function(x) cor.test(Sce_mean_Act[x,], Sce_mean_Expr[x,])[["estimate"]][["cor"]] )
correlation_table_exon$PcorAvg_m <- s
s <- sapply(correlation_table_exon$gene, function(x) cor.test(Sce_mean_Act[x,], Sce_mean_Expr[x,])[["p.value"]] )
correlation_table_exon$PcorAvg_m_pvalue <- s



CL <- cbind(Sce_mean_Act[,1], Sce_mean_Expr[,1])
colnames(CL) <- c("Activity", "Expression")
row <- number_correlated_genes(CL, name = paste0("CL_0"))
number_genes_table_exon <- rbind(row)
for (i in (1:length(colnames(Sce_mean_Expr)))){
  CL <- cbind(Sce_mean_Act[,i], Sce_mean_Expr[,i])
  colnames(CL) <- c("Activity", "Expression")
  cluster_scatter_plots(CL, name = paste0("CL_", i-1), log = TRUE)
  if(i>1){
    row <- number_correlated_genes(CL, name = paste0("CL_", i-1))
    number_genes_table_exon <- rbind(number_genes_table_exon, row)
  }
}


me <- mean (log(Sce_mean_Expr)[log(Sce_mean_Expr) > -Inf])
ma <- mean (log(Sce_mean_Act)[log(Sce_mean_Act) > -Inf])
CL <- cbind(Sce_mean_Act[,1], Sce_mean_Expr[,1])
colnames(CL) <- c("Activity", "Expression")
row <- number_correlated_genes_lines(CL,limex = me, limacc = ma, name = paste0("CL_0"), log = TRUE)
number_genes_table_exon <- rbind(row)
for (i in (1:length(colnames(Sce_mean_Expr)))){
  CL <- cbind(Sce_mean_Act[,i], Sce_mean_Expr[,i])
  colnames(CL) <- c("Activity", "Expression")
  cluster_scatter_plots(CL, name = paste0("CL_", i-1), log = TRUE)
  if(i>1){
    row <- number_correlated_genes_lines(CL,limex = me, limacc = ma, name = paste0("CL_", i-1), log = TRUE)
    number_genes_table_exon <- rbind(number_genes_table_exon, row)
  }
}



CL <- cbind(Sce_mean_Act[,1], Sce_mean_Expr[,1])
colnames(CL) <- c("Activity", "Expression")
row <- incoherence_factor(CL, name = paste0("CL_0"))
row <- (row/sum(row))*100
incoherence_factor_exon <- rbind(row)
for (i in (1:length(colnames(Sce_mean_Expr)))){
  CL <- cbind(Sce_mean_Act[,i], Sce_mean_Expr[,i])
  colnames(CL) <- c("Activity", "Expression")
  if(i>1){
    row <- incoherence_factor(CL, name = paste0("CL_", i-1), log = FALSE)
    row <- (row/sum(row))*100
    incoherence_factor_exon <- rbind(incoherence_factor_exon, row)
  }
}



mean(Sce_mean_Expr)
mean(Sce_mean_Act)

########### ENHD Contribution ##########


atac <- CreateAssayObject(enhd_contribution)
#SEU_RNA[["ENHD"]] <- atac
SEU_EHND <- CreateSeuratObject(atac, assay = "ENHD")
DefaultAssay(SEU_EHND) <- "ENHD"
SEU_EHND <- NormalizeData(SEU_EHND)
SEU_EHND <- AddMetaData(SEU_EHND, SEU_RNA@meta.data[["RNA_snn_res.0.8"]], col.name = "RNA_snn_res.0.8")

SCE <- as.SingleCellExperiment(SEU_EHND)
Sce_mean_Act <- aggregateData(SCE,  assay = "counts" , by = "RNA_snn_res.0.8", fun = "mean")
Sce_mean_Act <- Sce_mean_Act@assays@data@listData[[1]]
Sce_mean_Act <- Sce_mean_Act[active_prom_name,]

DefaultAssay(SEU_RNA) <- "RNA"
SCE <- as.SingleCellExperiment(SEU_RNA)
Sce_mean_Expr <- aggregateData(SCE,  assay = "counts" , by = "RNA_snn_res.0.8", fun = "mean")
Sce_mean_Expr <- Sce_mean_Expr@assays@data@listData[[1]]
Sce_mean_Expr <- Sce_mean_Expr[active_prom_name,]

correlation_table_enhd <- info_table_temp[,1:2]
correlation_table_enhd <- correlation_table_enhd[correlation_table_enhd$gene %in% active_prom_name,]
correlation_table_enhd[,2] <- NULL

s <- sapply(correlation_table_enhd$gene, function(x) cor.test(Sce_mean_Act[x,], Sce_mean_Expr[x,])[["estimate"]][["cor"]] )
correlation_table_enhd$PcorAvg_m <- s
s <- sapply(correlation_table_enhd$gene, function(x) cor.test(Sce_mean_Act[x,], Sce_mean_Expr[x,])[["p.value"]] )
correlation_table_enhd$PcorAvg_m_pvalue <- s



CL <- cbind(Sce_mean_Act[,1], Sce_mean_Expr[,1])
colnames(CL) <- c("Activity", "Expression")
row <- number_correlated_genes(CL, name = paste0("CL_0"))
number_genes_table_enhd <- rbind(row)
for (i in (1:length(colnames(Sce_mean_Expr)))){
  CL <- cbind(Sce_mean_Act[,i], Sce_mean_Expr[,i])
  colnames(CL) <- c("Activity", "Expression")
  cluster_scatter_plots(CL, name = paste0("CL_", i-1))
  if(i>1){
    row <- number_correlated_genes(CL, name = paste0("CL_", i-1))
    number_genes_table_enhd <- rbind(number_genes_table_enhd,row)
  }
}

CL <- cbind(Sce_mean_Act[,1], Sce_mean_Expr[,1])
colnames(CL) <- c("Activity", "Expression")
row <- incoherence_factor(CL, name = paste0("CL_0"))
row <- (row/sum(row))*100
incoherence_factor_enhd <- rbind(row)
for (i in (1:length(colnames(Sce_mean_Expr)))){
  CL <- cbind(Sce_mean_Act[,i], Sce_mean_Expr[,i])
  colnames(CL) <- c("Activity", "Expression")
  if(i>1){
    row <- incoherence_factor(CL, name = paste0("CL_", i-1), log = FALSE)
    row <- (row/sum(row))*100
    incoherence_factor_enhd <- rbind(incoherence_factor_enhd, row)
  }
}


#### LOG 

SCE <- as.SingleCellExperiment(SEU_EHND)
Sce_mean_Act <- aggregateData(SCE,  assay = "logcounts" , by = "RNA_snn_res.0.8", fun = "mean")
Sce_mean_Act <- Sce_mean_Act@assays@data@listData[[1]]
Sce_mean_Act <- Sce_mean_Act[active_prom_name,]

DefaultAssay(SEU_RNA) <- "RNA"
SCE <- as.SingleCellExperiment(SEU_RNA)
Sce_mean_Expr <- aggregateData(SCE,  assay = "logcounts" , by = "RNA_snn_res.0.8", fun = "mean")
Sce_mean_Expr <- Sce_mean_Expr@assays@data@listData[[1]]
Sce_mean_Expr <- Sce_mean_Expr[active_prom_name,]

correlation_table_enhd <- info_table_temp[,1:2]
correlation_table_enhd <- correlation_table_enhd[correlation_table_enhd$gene %in% active_prom_name,]
correlation_table_enhd[,2] <- NULL

s <- sapply(correlation_table_enhd$gene, function(x) cor.test(Sce_mean_Act[x,], Sce_mean_Expr[x,])[["estimate"]][["cor"]] )
correlation_table_enhd$PcorAvg_m <- s
s <- sapply(correlation_table_enhd$gene, function(x) cor.test(Sce_mean_Act[x,], Sce_mean_Expr[x,])[["p.value"]] )
correlation_table_enhd$PcorAvg_m_pvalue <- s



CL <- cbind(Sce_mean_Act[,1], Sce_mean_Expr[,1])
colnames(CL) <- c("Activity", "Expression")
row <- number_correlated_genes(CL, name = paste0("CL_0"))
number_genes_table_enhd <- rbind(row)
for (i in (1:length(colnames(Sce_mean_Expr)))){
  CL <- cbind(Sce_mean_Act[,i], Sce_mean_Expr[,i])
  colnames(CL) <- c("Activity", "Expression")
  cluster_scatter_plots(CL, name = paste0("CL_", i-1))
  if(i>1){
    row <- number_correlated_genes(CL, name = paste0("CL_", i-1))
    number_genes_table_enhd <- rbind(number_genes_table_enhd,row)
  }
}

CL <- cbind(Sce_mean_Act[,1], Sce_mean_Expr[,1])
colnames(CL) <- c("Activity", "Expression")
row <- incoherence_factor(CL, name = paste0("CL_0"))
row <- (row/sum(row))*100
incoherence_factor_enhd <- rbind(row)
for (i in (1:length(colnames(Sce_mean_Expr)))){
  CL <- cbind(Sce_mean_Act[,i], Sce_mean_Expr[,i])
  colnames(CL) <- c("Activity", "Expression")
  if(i>1){
    row <- incoherence_factor(CL, name = paste0("CL_", i-1), log = FALSE)
    row <- (row/sum(row))*100
    incoherence_factor_enhd <- rbind(incoherence_factor_enhd, row)
  }
}


###############

mean (log(Sce_mean_Expr)[log(Sce_mean_Expr) > -Inf])

mean ((Sce_mean_Expr)[(Sce_mean_Expr) > -Inf])







median(correlation_table$PcorAvg_m[1])


name = paste0("CL_", i-1)
png(paste0("../TMPResults/IMAGES/",name,"_gene_scatter_plot.png"), width = 1080, height = 1080)
ggplot(as.data.frame(log(CL_0)), aes(x=Activity, y = Expression)) + geom_point() + geom_hline(yintercept = colMeans(log(CL_0)[log(CL_0)[,2] > -Inf,])[2], linetype="dashed", color = "red", size=2) + 
  geom_vline(xintercept = colMeans(log(CL_0)[log(CL_0)[,1] > -Inf,])[1], linetype="dashed", color = "red", size=2)+ ggtitle(name)+ ggtitle("CL_0") +
  theme(plot.title = element_text(hjust = 0.5))#+
#geom_point(data=as.data.frame(log(CL_0))[high_cor_genes, ], aes(x=Activity, y = Expression), colour="red", size=2)
print(plot)
ggsave(path = "../TMPResults/IMAGES/", filename = paste0(name,"_gene_scatter_plot.png"), width = 1080, height = 1080, units= "px", scale = 3.5)

dev.off()







ggplot(as.data.frame(log(CL_0)), aes(x=Activity, y = Expression)) + geom_point() + geom_hline(yintercept = colMeans(log(CL_0)[log(CL_0)[,2] > -Inf,])[2], linetype="dashed", color = "red", size=2) + 
  geom_vline(xintercept = colMeans(log(CL_0)[log(CL_0)[,1] > -Inf,])[1], linetype="dashed", color = "red", size=2)#+
#geom_point(data=as.data.frame(log(CL_0))[high_cor_genes, ], aes(x=V1, y = V2), colour="red", size=2)


CL_0 <- cbind(Sce_mean_Act[,1], Sce_mean_Expr[,1])
colnames(CL_0) <- c("Activity", "Expression")
gene <- rbind(Sce_mean_Act["KLF4",], Sce_mean_Expr["KLF4",])
gene <- t(gene)

high_cor_genes <- correlation_table[correlation_table$PcorAvg_m > 0.8,]$gene

ggplot(as.data.frame((CL_0)), aes(x=Activity, y =Expression)) + geom_point() + geom_hline(yintercept = colMedians((CL_0))[2], linetype="dashed", color = "red", size=2) + geom_vline(xintercept = colMedians((CL_0))[1], linetype="dashed", color = "red", size=2)+
  geom_point(data=as.data.frame((CL_0))[high_cor_genes, ], aes(x=V1, y = V2), colour="red", size=2) +ylim(0,1)

ggplot(as.data.frame(log(CL_0)), aes(x=V1, y = V2)) + geom_point() + geom_hline(yintercept = colMeans(log(CL_0)[log(CL_0)[,2] > -Inf,])[2], linetype="dashed", color = "red", size=2) + geom_vline(xintercept = colMeans(log(CL_0)[log(CL_0)[,1] > -Inf,])[1], linetype="dashed", color = "red", size=2)+
  geom_point(data=as.data.frame(log(CL_0))[high_cor_genes, ], aes(x=V1, y = V2), colour="red", size=2)

ggplot(as.data.frame(gene), aes(x=V1, y = V2)) + geom_point()+ geom_text(label=rownames(as.data.frame(gene)))

colMeans(CL_0)
sum(log(CL_0)[,1] >= colMeans(log(CL_0)[log(CL_0)[,1] > -Inf,])[1] & log(CL_0)[,2] >= colMeans(log(CL_0)[log(CL_0)[,2] > -Inf,])[2])
sum(log(CL_0)[,1] >= colMeans(log(CL_0)[log(CL_0)[,1] > -Inf,])[1] & log(CL_0)[,2] < colMeans(log(CL_0)[log(CL_0)[,2] > -Inf,])[2])
sum(log(CL_0)[,1] < colMeans(log(CL_0)[log(CL_0)[,1] > -Inf,])[1] & log(CL_0)[,2] >= colMeans(log(CL_0)[log(CL_0)[,2] > -Inf,])[2])
sum(log(CL_0)[,1] < colMeans(log(CL_0)[log(CL_0)[,1] > -Inf,])[1] & log(CL_0)[,2] < colMeans(log(CL_0)[log(CL_0)[,2] > -Inf,])[2])

log(CL_0)[log(CL_0)[,1] > -Inf,]













temp_bin <- refseq_first_gene_matrix[active_prom_name,]
X <- temp_bin
temp_bin <- RNA_matrix[active_prom_name,]
Y <- temp_bin
s <- sapply(info_table_temp$gene, function(x) cor.test(X[x,], Y[x,])[["estimate"]][["cor"]] )
info_table_temp$Pcor_SC <- s

correlation_table <- info_table_temp[,1:2]
correlation_table[,2] <- NULL


filtered_info <- info_table_temp %>%
                  filter(PcorAgg_pvalue < 0.05) %>%
                  filter(PcorAvg_pvalue < 0.05)

ggscatter(filtered_info, x= "PcorAgg", y = "PcorAvg")
ggscatter(AggProm, x= "PcorAgg", y = "PcorAvg")

prova <- rbind(t(as.data.frame(AggProm[["RNA"]])),t(as.data.frame(AggProm[["RNA"]])))
prova <- t(prova)
pr <- cbind(as.data.frame(prova[,1]), as.data.frame(prova[,20]))
ggplot(pr) + geom_jitter(x= prova[,1], y = prova[,20])
ggscatter(prova, x= 1, y = 20 )





#########
s <- apply(refseq_promoter_gene_mat,1, function(x) names(which(x > 0)) )
info_table$prom_peaks <- s
s <- apply(refseq_allgene_gene_mat, 1, function(x) names(which(x > 0)))
info_table$ex_peaks <- NA
info_table[names(s),]$ ex_peaks <- s


s <- apply(refseq_allgene_gene_mat, 1, function(x) names(which(x > 0)))
write.table(s, file = "RI.csv")


for (i in rownames(refseq_promoter_gene_mat)){
  info_table[i,]$prom_peaks <- as.list(names(which(refseq_promoter_gene_mat[i,] > 0)))
}

s <- apply(refseq_promoter_gene_mat,1, function(x) names(which(x > 0)) )
info_table_temp$RPE <- s


RNA_ex <- CDS_RNA@assays@data@listData[["counts"]]["C3",]
temp_bin <- refseq_allgene_gene_matrix
temp_bin@x[temp_bin@x>0] <- 1
ACT_ex <- temp_bin["C3",]
sum(RNA_ex)
sum(ACT_ex)
open_ex_names <- names(which(ACT_ex > 0))
expr_ex_names <- names(which(RNA_ex > 0))
sum(expr_ex_names %in% open_ex_names)

RI(RNA_ex,ACT_ex)

CreateAssayObject(counts = refseq_allgene_gene_mat)

