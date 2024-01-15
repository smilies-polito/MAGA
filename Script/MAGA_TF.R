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

library(Matrix)
library(TFBSTools)
library(motifmatchr)
library(chromVAR)
library(JASPAR2020)
#library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BiocParallel)
library(Signac)
library(Seurat)
library(qusage)
library(tidyr)
library(dplyr)
library(cicero)
library(monocle3)
library(ggplot2)
#library(ggseqlogo)

matrix <- readMM("TMPDATA/10k_PBMC_Multiome_Controller/filtered_feature_bc_matrix/matrix.mtx")
cells <- read.table("TMPDATA/10k_PBMC_Multiome_Controller/filtered_feature_bc_matrix/barcodes.tsv")
features <- read.delim("TMPDATA/10k_PBMC_Multiome_Controller/filtered_feature_bc_matrix/features.tsv", header=FALSE)

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

genome_ref = read.table("Data/genomes/hg38/hg38.p13.chrom.sizes.txt")
genome_ref <- genome_ref[1:24,]
hg38 <- genome_ref[1:24,]
hg38 <- Seqinfo(hg38$V1, seqlengths= hg38$V2)
hg38@genome[] <- "hg38"

#creation of Suearta chromatin assay
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
rm(PBMC_assay)

#get human TF motif
human_motif <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', species = "Homo sapiens", all_versions = FALSE)
)

# peak labeling
labeled_peaks <- read.csv("Data/10x_PBMC_Multiome_Controller/labeled peaks/encodeCcreCombined_hg38_ucscLabel_classifiedPeaks.csv", sep = "\t")
#labeled_peaks_multi <- read.csv("C:/Users/loren/IWBBIO/data/labeled peaks/classifiedPeaks_multiCols.csv", sep = "\t")
nmax <- max(stringr::str_count(labeled_peaks$encodeCcreCombined_hg38_ucscLabel, "\t")) + 1
labeled_peaks <- separate(labeled_peaks, col = encodeCcreCombined_hg38_ucscLabel, sep = "\t", into = paste0("RegFunc", seq_len(nmax)))
labeled_peaks$site_names <- paste0(labeled_peaks$X.chrom, "_", labeled_peaks$chromStart, "_", labeled_peaks$chromEnd)
labeled_peaks <- labeled_peaks[labeled_peaks$site_names %in% rownames(fData(CDS_ATAC)),]
labeled_peaks <- labeled_peaks[!duplicated(labeled_peaks),]


###### PEAK LABELLING ###########
#labeled_peaks_multi <- read.csv("C:/Users/loren/IWBBIO/data/labeled peaks/classifiedPeaks_multiCols.csv", sep = "\t")
nmax <- max(stringr::str_count(labeled_peaks$encodeCcreCombined_hg38_ucscLabel, "\t")) + 1
labeled_peaks <- separate(labeled_peaks, col = encodeCcreCombined_hg38_ucscLabel, sep = "\t", into = paste0("RegFunc", seq_len(nmax)))
labeled_peaks$site_names <- paste0(labeled_peaks$X.chrom, "_", labeled_peaks$chromStart, "_", labeled_peaks$chromEnd)
labeled_peaks <- labeled_peaks[labeled_peaks$site_names %in% rownames(fData(CDS_ATAC)),]
labeled_peaks <- labeled_peaks[!duplicated(labeled_peaks),]

ATAC_matrix@x[ATAC_matrix@x > 0] <- 1
#Creation of the CDS object for the ATAC data
CDS_ATAC <- new_cell_data_set(ATAC_matrix)
rowData(CDS_ATAC)$gene_short_name <- peaks$V2
#the process and function are totally analogous to  before
CDS_ATAC <- detect_genes(CDS_ATAC)

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
refseq_peaks_prom <- peaks_multi_info[peaks_multi_info$PT == "TRUE",]
refseq_prom_list <- rownames(refseq_peaks_prom)
##########





prom_peaks_list <- data.frame(V1= refseq_near_prom_list) %>%  separate(col = V1, sep ="_", into = c("Chr","Start", "End")) %>%  mutate(Start = as.numeric(Start), End = as.numeric(End))
#prom_GRange <- GRanges(seqnames = prom_peaks_list$Chr, ranges =IRanges(start= prom_peaks_list$Start , end= prom_peaks_list$End))

#peak_motif <- matchMotifs(human_motif, prom_GRange, genome = BSgenome.Hsapiens.UCSC.hg38 , out = "scores")

mean(rowSums(peak_motif@assays@data@listData[["motifMatches"]]))
max(peak_motif@assays@data@listData[["motifScores"]])
max(rowSums(peak_motif@assays@data@listData[["motifMatches"]]))
rownames(CDS_RNA)


list <- data.frame(name = names(human_TF_Set))
list$short  <- sub("_.*", "", list$name)

sum(list$short %in% peak_motif@colData@listData[["name"]])
sum( peak_motif@colData@listData[["name"]] %in% list$short)


both_list <- list[list$short %in% peak_motif@colData@listData[["name"]],]
unique_TF_list <- unique(both_list$short)
human_TF_Set_short <- human_TF_Set[both_list$name]

is <- as.data.frame(t(as.data.frame(lapply(human_motif_short, function(x){x@name}))))
is$motif <- rownames(is)
#both_list$motif <- apply(both_list, MARGIN = 1, function(x){rownames(is[is$V1 == x[2],])})
#both_list$name_trattino <- gsub("_","-", both_list$name)  


#human_motif_short <- human_motif[lapply(human_motif, function(x){print(x@name)}) %in% list$short]

# motif enrichment calculation for promoters
PBMC <- PBMC[rownames(PBMC) %in% gsub("_","-",refseq_near_prom_list),]
PBMC <- AddMotifs(
  object = PBMC,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = human_motif,
  assay = "EnhdPeaks"
)

register(SerialParam())
motif_activity <- RunChromVAR(
  object = PBMC,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  new.assay.name = "chromvar_full"
)

# motif enrichment calculation for enhancers
motif_activity[["EnhdPeaks"]] <- CreateChromatinAssay(
  counts = ATAC_matrix[enhD_peaks_list,],
  sep = c("_", "_"),
  genome = hg38,
  min.cells = 1
)
motif_activity <- AddMotifs(
  object = motif_activity,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = human_motif,
  assay = "EnhdPeaks"
)
DefaultAssay(motif_activity) <- "EnhdPeaks"
register(SerialParam())
motif_activity <- RunChromVAR(
  object = motif_activity,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  new.assay.name = "chromvar_enhd",
  assay= "EnhdPeaks"
)

#motif_markers_enhd <- FindAllMarkers(motif_activity, assay = "chromvar_enhd")
#motif_markers_enhd$TF <- lapply(motif_markers_enhd$gene, function(x){is[is$motif == x,]$V1})
#motif_markers_enhd$TF <- as.character(motif_markers_enhd$TF)

########

# RNA data processing
motif_activity[["RNA"]] <- CreateAssayObject(
  counts = RNA_matrix,
  min.cells = 1,
)
DefaultAssay(motif_activity) <- "RNA"
motif_activity <- NormalizeData(motif_activity, normalization.method = "LogNormalize", scale.factor = 10000)
motif_activity <- FindVariableFeatures(motif_activity, selection.method = "vst")
motif_activity <- FindVariableFeatures(motif_activity, selection.method = "mvp")
motif_activity <- ScaleData(motif_activity)
motif_activity <- RunPCA(motif_activity, features = VariableFeatures(object = motif_activity))
motif_activity <- FindNeighbors(motif_activity, dims = 1:10)
motif_activity <- FindClusters(motif_activity, resolution = 0.3)
motif_activity <- RunUMAP(motif_activity, dims = 1:10)
#DimPlot(motif_activity, reduction = "umap")

DefaultAssay(motif_activity) <- "RNA"
pbmc_rna <- readRDS("TMPDATA/pbmc_10k_v3.rds")
transfer.anchors <- FindTransferAnchors(
  reference = pbmc_rna,
  query = motif_activity,
  reduction = 'cca'
)
predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = pbmc_rna$celltype,
  weight.reduction = motif_activity[['pca']],
  dims = 2:30
)

motif_activity <- AddMetaData(object = motif_activity, metadata = predicted.labels)
DimPlot(motif_activity,group.by = "predicted.id")
motif_activity@meta.data[["predicted.id"]] <- as.factor(motif_activity@meta.data[["predicted.id"]])

DefaultAssay(motif_activity) <- 'chromvar_full'
#motif_activity <- NormalizeData(motif_activity, normalization.method = "LogNormalize", scale.factor = 10000)
#motif_activity <- FindVariableFeatures(motif_activity, selection.method = "vst")
#motif_activity <- FindVariableFeatures(motif_activity, selection.method = "mvp")
#motif_markers_prom <- FindAllMarkers(motif_activity, assay = "chromvar_full")
#motif_markers_prom$TF <- lapply(motif_markers_prom$gene, function(x){is[is$motif == x,]$V1})
#motif_markers_prom$TF <- as.character(motif_markers_prom$TF)


####################### SHORT VERSION + ENRICH ############


unique_motif_list <- both_list[!duplicated(both_list$short),]$motif
variance_table <- data.frame(unique_motif_list)
variance_ratio_table <- data.frame(unique_motif_list)
mean_table <- data.frame(unique_motif_list)
SEU_tot_var <- motif_activity@assays[["chromvar"]]@meta.features[,c(2,3)][unique_motif_list,]
DefaultAssay(motif_activity) <- "chromvar"
for (i in levels(motif_activity@meta.data[["RNA_snn_res.0.3"]])){
  sub <- subset(motif_activity, idents = i)
  #sub <- FindVariableFeatures(sub)
  sub <- FindVariableFeatures(sub, selection.method = "mvp")
  table_var <- sub@assays[["chromvar"]]@meta.features[,c(2,3)][unique_motif_list,]
  table_expr <- sub@assays[["chromvar"]]@meta.features[,c(1,1)][unique_motif_list,]
  table_ratio <- table_var/SEU_tot_var[unique_motif_list,]
  colnames(table_var) <- paste0("CL",i,"_", colnames(table_var))
  colnames(table_expr) <- paste0("CL",i,"_", colnames(table_expr))
  colnames(table_ratio) <- paste0("CL",i,"_", colnames(table_ratio))
  variance_table <- cbind(variance_table, table_var)
  mean_table <- cbind(mean_table, table_expr)
  variance_ratio_table <- cbind(variance_ratio_table, table_ratio)
}

variance_table_motif <- variance_table
mean_table_motif <- mean_table
variance_ratio_table_motif <- variance_ratio_table


variance_table <- data.frame(unique_TF_list)
variance_ratio_table <- data.frame(unique_TF_list)
mean_table <- data.frame(unique_TF_list)
SEU_tot_var <- motif_activity@assays[["RNA"]]@meta.features[,c(2,4,7)]
DefaultAssay(motif_activity)<- "RNA"
for (i in levels(motif_activity@meta.data[["RNA_snn_res.0.3"]])){
  sub <- subset(motif_activity, idents = i)
  sub <- FindVariableFeatures(sub)
  sub <- FindVariableFeatures(sub, selection.method = "mvp")
  table_var <- sub@assays[["RNA"]]@meta.features[,c(2,4,7)][unique_TF_list,]
  table_expr <- sub@assays[["RNA"]]@meta.features[,c(1,6)][unique_TF_list,]
  table_ratio <- table_var/SEU_tot_var[unique_TF_list,]
  colnames(table_var) <- paste0("CL",i,"_", colnames(table_var))
  colnames(table_expr) <- paste0("CL",i,"_", colnames(table_expr))
  colnames(table_ratio) <- paste0("CL",i,"_", colnames(table_ratio))
  variance_table <- cbind(variance_table, table_var)
  mean_table <- cbind(mean_table, table_expr)
  variance_ratio_table <- cbind(variance_ratio_table, table_ratio)
}

variance_table_expr <- variance_table
variance_table_expr[is.na(variance_table_expr)] <- 0
rownames(variance_table_expr) <- variance_table_expr$unique_TF_list

mean_table_expr <- mean_table
mean_table_expr[is.na(mean_table_expr)] <- 0
rownames(mean_table_expr) <- mean_table_expr$unique_TF_list

variance_ratio_table_expr <- variance_ratio_table
variance_ratio_table_expr[is.na(variance_ratio_table_expr)] <- 0
rownames(variance_ratio_table_expr) <- variance_ratio_table_expr$unique_TF_list


unique_enrich_list <- both_list$name
variance_table <- data.frame(unique_enrich_list)
variance_ratio_table <- data.frame(unique_enrich_list)
mean_table <- data.frame(unique_enrich_list)
SEU_tot_var <- motif_activity@assays[["ENrich"]]@meta.features[,c(2,4,7)]
DefaultAssay(motif_activity)<- "ENrich"
for (i in levels(motif_activity@meta.data[["RNA_snn_res.0.3"]])){
  sub <- subset(motif_activity, idents = i)
  sub <- FindVariableFeatures(sub)
  sub <- FindVariableFeatures(sub, selection.method = "mvp")
  table_var <- sub@assays[["ENrich"]]@meta.features[,c(2,4,7)]
  table_expr <- sub@assays[["ENrich"]]@meta.features[,c(1,6)]
  table_ratio <- table_var/SEU_tot_var
  colnames(table_var) <- paste0("CL",i,"_", colnames(table_var))
  colnames(table_expr) <- paste0("CL",i,"_", colnames(table_expr))
  colnames(table_ratio) <- paste0("CL",i,"_", colnames(table_ratio))
  variance_table <- cbind(variance_table, table_var)
  mean_table <- cbind(mean_table, table_expr)
  variance_ratio_table <- cbind(variance_ratio_table, table_ratio)
}

variance_table_enrich <- variance_table
mean_table_enrich <- mean_table
variance_ratio_table_enrich <- variance_ratio_table




key = "_mvp.mean"
#key = "_vst.mean"
Sce_mean_motif = mean_table_motif[,paste0("CL", 0, key)]
Sce_mean_Expr = mean_table_expr[,paste0("CL", 0, key)]
for (i in (2:((length(colnames(mean_table_expr))-1)/2))){
  
  Sce_mean_motif = cbind(Sce_mean_motif, mean_table_motif[,paste0("CL", i-1, key)])
  Sce_mean_Expr = cbind (Sce_mean_Expr, mean_table_expr[,paste0("CL", i-1, key)])
  
}
rownames(Sce_mean_motif) <- mean_table_motif$active_prom_name
rownames(Sce_mean_Expr) <- mean_table_expr$unique_TF_list



#motif vs expression
CL <- cbind(mean_table_motif[,paste0("CL",0 , key)], mean_table_expr[,paste0("CL", 0, key)])

colnames(CL) <- c("Activity", "Expression")
row <- number_correlated_genes(CL, name = paste0("CL_0"), log = FALSE)
number_genes_table_prom <- rbind(row)
for (i in (1:((length(colnames(mean_table_expr))-1)/2))){
  
  CL <- cbind(mean_table_motif[,paste0("CL", i-1, key)], mean_table_expr[,paste0("CL", i-1, key)])
  colnames(CL) <- c("Activity", "Expression")
  rownames(CL) <- rownames(mean_table_expr)

  
  #low_variance_genes <- (top_n(variance_ratio_table_expr[rownames(CL[CL[,"Expression"] > 0,]),c("active_prom_name",paste0("CL", i-1, "_vst.variance.standardized"))], -1000)["active_prom_name"])$active_prom_name
  #low_variance_genes_acc <- (top_n(variance_ratio_table_motif[rownames(CL[CL[,"Activity"] > 0,]),c("active_prom_name",paste0("CL", i-1, "_vst.variance.standardized"))], -1000)["active_prom_name"])$active_prom_name
  #markers_genes_cl <- expr_markers %>% filter(cluster == i-1) %>% filter(avg_log2FC > 0)
  #markers_genes_cl <- markers_genes_cl$gene
  cluster_scatter_plots(CL, name = paste0("CL", i-1),highlight_low_var = C(), log = FALSE, folder_name = "MOTIF")
  if(i>1){
    row <- number_correlated_genes(CL, name = paste0("CL_", i-1), log = FALSE)
    number_genes_table_prom <- rbind(number_genes_table_prom,row)
  }
}


#enrich vs expression

CL <- select(both_list, c("name_trattino", "short")) %>% mutate(Activity = mean_table_enrich[,paste0("CL",0 , key)]) %>% rowwise() %>%  mutate(Expression = mean_table_expr[short,paste0("CL", 0, key)] ) %>% select(Activity, Expression)
row <- number_correlated_genes(CL, name = paste0("CL_0"), log = FALSE)
number_genes_table_prom <- rbind(row)

for (i in (1:((length(colnames(mean_table_expr))-1)/2))){
  
  CL <- select(both_list, c("name_trattino", "short")) %>% mutate(Activity = mean_table_enrich[,paste0("CL",i-1 , key)]) %>% rowwise() %>%  mutate(Expression = mean_table_expr[short,paste0("CL", i-1, key)] ) %>% select(Activity, Expression)
  
  rownames(CL) <- rownames(mean_table_enrich)
  
  
  #low_variance_genes <- (top_n(variance_ratio_table_expr[rownames(CL[CL[,"Expression"] > 0,]),c("active_prom_name",paste0("CL", i-1, "_vst.variance.standardized"))], -1000)["active_prom_name"])$active_prom_name
  #low_variance_genes_acc <- (top_n(variance_ratio_table_motif[rownames(CL[CL[,"Activity"] > 0,]),c("active_prom_name",paste0("CL", i-1, "_vst.variance.standardized"))], -1000)["active_prom_name"])$active_prom_name
  #markers_genes_cl <- expr_markers %>% filter(cluster == i-1) %>% filter(avg_log2FC > 0)
  #markers_genes_cl <- markers_genes_cl$gene
  cluster_scatter_plots(CL, name = paste0("CL", i-1),highlight_low_var = C(), log = FALSE, folder_name = "ENRICH")
  if(i>1){
    row <- number_correlated_genes(CL, name = paste0("CL_", i-1), log = FALSE)
    number_genes_table_prom <- rbind(number_genes_table_prom,row)
  }
}

#motif vs enrich

CL <- select(both_list, c("name_trattino", "motif")) %>% mutate(Activity = mean_table_enrich[,paste0("CL",0 , key)]) %>% rowwise() %>%  mutate(Expression = mean_table_motif[motif,paste0("CL", 0, key)] ) %>% select(Activity, Expression)
row <- number_correlated_genes(CL, name = paste0("CL_0"), log = FALSE)
number_genes_table_prom <- rbind(row)

for (i in (1:((length(colnames(mean_table_expr))-1)/2))){
  
  CL <- select(both_list, c("name_trattino", "motif")) %>% mutate(Activity = mean_table_enrich[,paste0("CL",i-1 , key)]) %>% rowwise() %>%  mutate(Expression = mean_table_motif[motif,paste0("CL", i-1, key)] ) %>% select(Activity, Expression)
  
  rownames(CL) <- rownames(mean_table_enrich)
  
  
  #low_variance_genes <- (top_n(variance_ratio_table_expr[rownames(CL[CL[,"Expression"] > 0,]),c("active_prom_name",paste0("CL", i-1, "_vst.variance.standardized"))], -1000)["active_prom_name"])$active_prom_name
  #low_variance_genes_acc <- (top_n(variance_ratio_table_motif[rownames(CL[CL[,"Activity"] > 0,]),c("active_prom_name",paste0("CL", i-1, "_vst.variance.standardized"))], -1000)["active_prom_name"])$active_prom_name
  #markers_genes_cl <- expr_markers %>% filter(cluster == i-1) %>% filter(avg_log2FC > 0)
  #markers_genes_cl <- markers_genes_cl$gene
  cluster_scatter_plots(CL, name = paste0("CL", i-1),highlight_low_var = C(), log = FALSE, folder_name = "ENRICHvsMOTIF")
  if(i>1){
    row <- number_correlated_genes(CL, name = paste0("CL_", i-1), log = FALSE)
    number_genes_table_prom <- rbind(number_genes_table_prom,row)
  }
}

####################### FULL MOTIF #######

is <- as.data.frame(t(as.data.frame(lapply(human_motif, function(x){x@name}))))
is$motif <- rownames(is)
is_exprs <- is[is$V1 %in% rownames(motif_activity[["RNA"]]),]

variance_table <- data.frame(motif = is$motif)
variance_ratio_table <- data.frame(motif = is$motif)
mean_table <- data.frame(motif = is$motif)
SEU_tot_var <- motif_activity@assays[["chromvar_full"]]@meta.features[,c(2,3)]#[unique_motif_list,]

DefaultAssay(motif_activity) <- "chromvar_full"

DefaultAssay(motif_activity) <- "chromvar_full"
for (i in levels(motif_activity@meta.data[["predicted.id"]])){
  sub <- subset(motif_activity, idents = i)
  #sub <- FindVariableFeatures(sub)
  sub <- FindVariableFeatures(sub, selection.method = "mvp")
  table_var <- sub@assays[["chromvar_full"]]@meta.features[,c(2,3)]#[unique_motif_list,]
  table_expr <- sub@assays[["chromvar_full"]]@meta.features[,c(1,1)]#[unique_motif_list,]
  table_expr[2]<- NULL 
  table_ratio <- table_var/SEU_tot_var#[unique_motif_list,]
  colnames(table_var) <- paste0("CL",i,"_", colnames(table_var))
  colnames(table_expr) <- paste0("CL",i,"_", colnames(table_expr))
  colnames(table_ratio) <- paste0("CL",i,"_", colnames(table_ratio))
  variance_table <- cbind(variance_table, table_var)
  mean_table <- cbind(mean_table, table_expr)
  variance_ratio_table <- cbind(variance_ratio_table, table_ratio)
}

variance_table_motif_prom <- variance_table
mean_table_motif_prom <- mean_table
variance_ratio_table_motif_prom <- variance_ratio_table

for (i in levels(motif_activity@meta.data[["RNA_snn_res.0.3"]])){
  sub <- subset(motif_activity, idents = i)
  #sub <- FindVariableFeatures(sub)
  sub <- FindVariableFeatures(sub, selection.method = "mvp")
  table_var <- sub@assays[["chromvar_full"]]@meta.features[,c(2,3)]#[unique_motif_list,]
  table_expr <- sub@assays[["chromvar_full"]]@meta.features[,c(1,1)]#[unique_motif_list,]
  table_ratio <- table_var/SEU_tot_var#[unique_motif_list,]
  colnames(table_var) <- paste0("CL",i,"_", colnames(table_var))
  colnames(table_expr) <- paste0("CL",i,"_", colnames(table_expr))
  colnames(table_ratio) <- paste0("CL",i,"_", colnames(table_ratio))
  variance_table <- cbind(variance_table, table_var)
  mean_table <- cbind(mean_table, table_expr)
  variance_ratio_table <- cbind(variance_ratio_table, table_ratio)
}

variance_table_motif <- variance_table
mean_table_motif <- mean_table
variance_ratio_table_motif_ <- variance_ratio_table


variance_table <- data.frame(TF = is$V1)
variance_ratio_table <- data.frame(TF = is$V1)
mean_table <- data.frame(TF = is$V1)
SEU_tot_var <- motif_activity@assays[["RNA"]]@meta.features[,c(2,4,7)]
DefaultAssay(motif_activity)<- "RNA"
for (i in levels(motif_activity@meta.data[["predicted.id"]])){
  sub <- subset(motif_activity, idents = i)
  sub <- FindVariableFeatures(sub)
  sub <- FindVariableFeatures(sub, selection.method = "mvp")
  table_var <- sub@assays[["RNA"]]@meta.features[,c(2,4,7)][is$V1,]
  table_expr <- sub@assays[["RNA"]]@meta.features[,c(1,6)][is$V1,]
  table_ratio <- table_var/SEU_tot_var[is$V1,]
  colnames(table_var) <- paste0("CL",i,"_", colnames(table_var))
  colnames(table_expr) <- paste0("CL",i,"_", colnames(table_expr))
  colnames(table_ratio) <- paste0("CL",i,"_", colnames(table_ratio))
  variance_table <- cbind(variance_table, table_var)
  mean_table <- cbind(mean_table, table_expr)
  variance_ratio_table <- cbind(variance_ratio_table, table_ratio)
}

variance_table_expr <- variance_table
variance_table_expr[is.na(variance_table_expr)] <- 0
rownames(variance_table_expr) <- variance_table_expr$TF

mean_table_expr <- mean_table
mean_table_expr[is.na(mean_table_expr)] <- 0
rownames(mean_table_expr) <- mean_table_expr$TF

variance_ratio_table_expr <- variance_ratio_table
variance_ratio_table_expr[is.na(variance_ratio_table_expr)] <- 0
rownames(variance_ratio_table_expr) <- variance_ratio_table_expr$TF


key = "_mvp.mean"
#key = "_vst.mean"
Sce_mean_motif = mean_table_motif[,paste0("CL", 0, key)]
Sce_mean_Expr = mean_table_expr[,paste0("CL", 0, key)]
for (i in (2:((length(colnames(mean_table_expr))-1)/2))){
  
  Sce_mean_motif = cbind(Sce_mean_motif, mean_table_motif[,paste0("CL", i-1, key)])
  Sce_mean_Expr = cbind (Sce_mean_Expr, mean_table_expr[,paste0("CL", i-1, key)])
  
}
rownames(Sce_mean_motif) <- mean_table_motif$active_prom_name
rownames(Sce_mean_Expr) <- mean_table_expr$unique_TF_list



#motif vs expression
CL <- cbind(mean_table_motif[,paste0("CL",0 , key)], mean_table_expr[,paste0("CL", 0, key)])
colnames(CL) <- c("Activity", "Expression")
row <- number_correlated_genes(CL, name = paste0("CL_0"), log = FALSE)
number_genes_table_prom <- rbind(row)

for (i in levels(motif_activity@meta.data[["predicted.id"]])){
  
  CL <- cbind(mean_table_motif_prom[,paste0("CL", i, key)], mean_table_expr[,paste0("CL", i, key)])
  colnames(CL) <- c("Activity", "Expression")
  rownames(CL) <- rownames(mean_table_expr)
  #low_variance_genes <- (top_n(variance_ratio_table_expr[rownames(CL[CL[,"Expression"] > 0,]),c("active_prom_name",paste0("CL", i-1, "_vst.variance.standardized"))], -1000)["active_prom_name"])$active_prom_name
  #low_variance_genes_acc <- (top_n(variance_ratio_table_motif[rownames(CL[CL[,"Activity"] > 0,]),c("active_prom_name",paste0("CL", i-1, "_vst.variance.standardized"))], -1000)["active_prom_name"])$active_prom_name
  #markers_genes_cl <- expr_markers %>% filter(cluster == i-1) %>% filter(avg_log2FC > 0)
  #markers_genes_cl <- markers_genes_cl$gene
  cluster_scatter_plots(CL, name = paste0(i),highlight_low_var = NULL, log = FALSE, folder_name = "MOTIF_FULL", quartile = 4, X= "Motif Enrichment (Promoter)", Y = "Expresssion")
}

if(i>1){
  row <- number_correlated_genes(CL, name = paste0("CL_", i-1), log = FALSE)
  number_genes_table_prom <- rbind(number_genes_table_prom,row)
}

levels <- levels(motif_activity@meta.data[["predicted.id"]])
key = "_mvp.mean"
#key = "_vst.mean"
Sce_mean_motif = mean_table_motif_prom[,paste0("CL", levels[1], key)]
Sce_mean_Expr = mean_table_expr[,paste0("CL", levels[1], key)]
for (i in levels[2:length(levels)]){
  
  Sce_mean_motif = cbind(Sce_mean_motif, mean_table_motif_prom[,paste0("CL", i, key)])
  Sce_mean_Expr = cbind (Sce_mean_Expr, mean_table_expr[,paste0("CL", i, key)])
  
}
rownames(Sce_mean_motif) <- is$V1
rownames(Sce_mean_Expr) <- is$V1

correlation_table_final_prom <- is_exprs[,1:2]
correlation_table_final_prom[,2] <- NULL


s <- sapply(correlation_table_final_prom$V1, function(x) cor.test(Sce_mean_motif[x,], Sce_mean_Expr[x,])[["estimate"]][["cor"]] )
correlation_table_final_prom$PcorAvg_m <- s
s <- sapply(correlation_table_final_prom$V1, function(x) cor.test(Sce_mean_motif[x,], Sce_mean_Expr[x,])[["p.value"]] )
correlation_table_final_prom$PcorAvg_m_pvalue <- s
dim(correlation_table_final_prom[correlation_table_final_prom$PcorAvg_m > 0.5,])[1]/dim(correlation_table_final_prom)[1]

top_variance_expr <- motif_activity@assays[["RNA"]]@meta.features[is_exprs$V1,c(2,4,7)] %>%  top_n(50, wt = vst.variance)
top_variance_expr <- rownames(top_variance_expr)
mean(correlation_table_final_prom[correlation_table_final_prom$V1 %in% top_variance_expr,]$PcorAvg_m)


mean(apply(mean_table_motif_prom[,2:ncol(mean_table_motif_prom)], 1 , sd))
mean(apply(mean_table_motif_prom[,2:ncol(mean_table_motif_prom)], 1 , mean))

prova <- mean_table_motif_prom[,2:ncol(mean_table_motif_prom)]
colnames(prova) <- levels
s1 <-gather(prova[-13], key = "CellType", value = "Expression")  %>% ggplot( aes(x = CellType, y = Expression, fill = CellType)) +
  geom_violin() +
  geom_jitter(color = "black", size = 0.1) + geom_hline(yintercept = 1) + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(title = "Promoter Motif Enrichment across Cell Types", 
       x = "", 
       y = "Motif Enrichment") + theme(axis.text.x = element_blank())


prova2 <- mean_table_motif_enhd[,2:ncol(mean_table_motif_enhd)]
colnames(prova2) <- levels
s2 <- gather(prova2[-13], key = "CellType", value = "Expression")  %>% ggplot( aes(x = CellType, y = Expression, fill = CellType)) +
  geom_violin() +
  geom_jitter(color = "black", size = 0.1) + geom_hline(yintercept = 1) + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(title = "Enhacer Motif Enrichment  across Cell Types", 
       x = "Cell Types", 
       y = "Motif Enrichment") + NoLegend()

s1 / s2
ggsave(path = paste0("TMPResults/IMAGES/", "MISC"), filename = "all_violin.pdf", width = 1080, height = 700, units= "px",scale = 3.5)


############# ENHANCERS ######


is <- as.data.frame(t(as.data.frame(lapply(human_motif, function(x){x@name}))))
is$motif <- rownames(is)

motif_activity <- FindVariableFeatures(motif_activity, selection.method = "mvp", assay = "chromvar_enhd")


variance_table <- data.frame(motif = is$motif)
variance_ratio_table <- data.frame(motif = is$motif)
mean_table <- data.frame(motif = is$motif)
SEU_tot_var <- motif_activity@assays[["chromvar_enhd"]]@meta.features[,c(2,3)]#[unique_motif_list,]
DefaultAssay(motif_activity) <- "chromvar_enhd"

for (i in levels(motif_activity@meta.data[["predicted.id"]])){
  sub <- subset(motif_activity, idents = i)
  #sub <- FindVariableFeatures(sub)
  sub <- FindVariableFeatures(sub, selection.method = "mvp")
  table_var <- sub@assays[["chromvar_enhd"]]@meta.features[,c(2,3)]#[unique_motif_list,]
  table_expr <- sub@assays[["chromvar_enhd"]]@meta.features[,c(1,1)]#[unique_motif_list,]
  table_expr[2]<- NULL
  table_ratio <- table_var/SEU_tot_var#[unique_motif_list,]
  colnames(table_var) <- paste0("CL",i,"_", colnames(table_var))
  colnames(table_expr) <- paste0("CL",i,"_", colnames(table_expr))
  colnames(table_ratio) <- paste0("CL",i,"_", colnames(table_ratio))
  variance_table <- cbind(variance_table, table_var)
  mean_table <- cbind(mean_table, table_expr)
  variance_ratio_table <- cbind(variance_ratio_table, table_ratio)
}

DefaultAssay(mid) <- "chromvar_enhd"
for (i in levels(mid@meta.data[["predicted.id"]])){
  sub <- subset(mid, idents = i)
  #sub <- FindVariableFeatures(sub)
  sub <- FindVariableFeatures(sub, selection.method = "mvp")
  table_var <- sub@assays[["chromvar_enhd"]]@meta.features[,c(2,3)]#[unique_motif_list,]
  table_expr <- sub@assays[["chromvar_enhd"]]@meta.features[,c(1,1)]#[unique_motif_list,]
  table_expr[2]<- NULL
  table_ratio <- table_var/SEU_tot_var#[unique_motif_list,]
  colnames(table_var) <- paste0("CL",i,"_", colnames(table_var))
  colnames(table_expr) <- paste0("CL",i,"_", colnames(table_expr))
  colnames(table_ratio) <- paste0("CL",i,"_", colnames(table_ratio))
  variance_table <- cbind(variance_table, table_var)
  mean_table <- cbind(mean_table, table_expr)
  variance_ratio_table <- cbind(variance_ratio_table, table_ratio)
}


variance_table_motif_enhd <- variance_table
mean_table_motif_enhd <- mean_table
variance_ratio_table_motif_enhd <- variance_ratio_table


variance_table <- data.frame(TF = is$V1)
variance_ratio_table <- data.frame(TF = is$V1)
mean_table <- data.frame(TF = is$V1)
SEU_tot_var <- motif_activity@assays[["RNA"]]@meta.features[,c(2,4,7)]
DefaultAssay(motif_activity)<- "RNA"

for (i in levels(motif_activity@meta.data[["RNA_snn_res.0.3"]])){
  sub <- subset(motif_activity, idents = i)
  sub <- FindVariableFeatures(sub)
  sub <- FindVariableFeatures(sub, selection.method = "mvp")
  table_var <- sub@assays[["RNA"]]@meta.features[,c(2,4,7)][is$V1,]
  table_expr <- sub@assays[["RNA"]]@meta.features[,c(1,6)][is$V1,]
  table_ratio <- table_var/SEU_tot_var[is$V1,]
  colnames(table_var) <- paste0("CL",i,"_", colnames(table_var))
  colnames(table_expr) <- paste0("CL",i,"_", colnames(table_expr))
  colnames(table_ratio) <- paste0("CL",i,"_", colnames(table_ratio))
  variance_table <- cbind(variance_table, table_var)
  mean_table <- cbind(mean_table, table_expr)
  variance_ratio_table <- cbind(variance_ratio_table, table_ratio)
}



Idents(motif_activity) <- "predicted.id"
mid <- motif_activity
mid[["peaks"]] <- NULL
mid[["EnhdPeaks"]] <- NULL
DefaultAssay(mid) <- "RNA"
for (i in levels(mid@meta.data[["predicted.id"]])){
  sub <- subset(mid, idents = i)
  if(length(colnames(sub)) > 1){
    sub <- FindVariableFeatures(sub)
  }
  sub <- FindVariableFeatures(sub, selection.method = "mvp")
  table_var <- sub@assays[["RNA"]]@meta.features[,c(2,4,7)][is$V1,]
  table_expr <- sub@assays[["RNA"]]@meta.features[,c(1,6)][is$V1,]
  table_ratio <- table_var/SEU_tot_var[is$V1,]
  colnames(table_var) <- paste0("CL",i,"_", colnames(table_var))
  colnames(table_expr) <- paste0("CL",i,"_", colnames(table_expr))
  colnames(table_ratio) <- paste0("CL",i,"_", colnames(table_ratio))
  variance_table <- cbind(variance_table, table_var)
  mean_table <- cbind(mean_table, table_expr)
  variance_ratio_table <- cbind(variance_ratio_table, table_ratio)
}

variance_table_expr <- variance_table
variance_table_expr[is.na(variance_table_expr)] <- 0
rownames(variance_table_expr) <- variance_table_expr$TF

mean_table_expr <- mean_table
mean_table_expr[is.na(mean_table_expr)] <- 0
rownames(mean_table_expr) <- mean_table_expr$TF

variance_ratio_table_expr <- variance_ratio_table
variance_ratio_table_expr[is.na(variance_ratio_table_expr)] <- 0
rownames(variance_ratio_table_expr) <- variance_ratio_table_expr$TF


tf_markers <- expr_markers[expr_markers$gene %in% is$V1,]


#motif vs expression
zero <- levels(motif_activity@meta.data[["predicted.id"]])[1]
CL <- cbind(mean_table_motif[,paste0("CL",zero , key)], mean_table_expr[,paste0("CL", zero, key)])
colnames(CL) <- c("Activity", "Expression")
row <- number_correlated_genes(CL, name = paste0("CL_", zero), log = FALSE)
number_genes_table_prom <- rbind(row)

for (i in levels(motif_activity@meta.data[["predicted.id"]])){
  
  CL <- cbind(mean_table_motif_enhd[,paste0("CL", i, key)], mean_table_expr[,paste0("CL", i, key)])
  colnames(CL) <- c("Activity", "Expression")
  rownames(CL) <- rownames(mean_table_expr)
  
  #low_variance_genes <- top_n(variance_ratio_table_expr[rownames(CL[CL[,"Expression"] > 0,]),c("TF",paste0("CL", i-1, "_vst.variance.standardized"))], -60)["TF"]$TF
  #low_variance_genes_acc <- top_n(variance_ratio_table_motif[,c("motif",paste0("CL", i-1, "_mvp.dispersion.scaled"))], -60)["motif"]$motif
  #low_variance_genes_acc <- top_n(variance_table_motif[,c("motif",paste0("CL", i-1, "_mvp.dispersion.scaled"))], -60)["motif"]$motif
  #low_variance_genes_acc <- is[is$motif %in% low_variance_genes_acc,]$V1
  #markers_genes_cl <- tf_markers %>% filter(cluster == i) %>% filter(avg_log2FC > 0)
  #markers_genes_cl <- markers_genes_cl$gene
  #markers_motif_cl <- motif_markers_enhd %>% filter(cluster == i) %>% filter(avg_log2FC > 0)
  #markers_motif_cl <- markers_motif_cl$TF
  #highlight <- markers_genes_cl[markers_genes_cl %in% markers_motif_cl]
  cluster_scatter_plots(CL, name = paste0( i),highlight_low_var = NULL, log = FALSE, folder_name = "MOTIF_ENHD", quartile = 4, X= "Motif Enrichment (Enhancer)", Y = "Expresssion")
}


for (i in (1:((length(colnames(mean_table_expr))-1)/2))){
  
  CL <- cbind(mean_table_motif[,paste0("CL", i-1, key)], mean_table_expr[,paste0("CL", i-1, key)])
  colnames(CL) <- c("Activity", "Expression")
  rownames(CL) <- rownames(mean_table_expr)
  
  #low_variance_genes <- top_n(variance_ratio_table_expr[rownames(CL[CL[,"Expression"] > 0,]),c("TF",paste0("CL", i-1, "_vst.variance.standardized"))], -60)["TF"]$TF
  #low_variance_genes_acc <- top_n(variance_ratio_table_motif[,c("motif",paste0("CL", i-1, "_mvp.dispersion.scaled"))], -60)["motif"]$motif
  #low_variance_genes_acc <- top_n(variance_table_motif[,c("motif",paste0("CL", i-1, "_mvp.dispersion.scaled"))], -60)["motif"]$motif
  #low_variance_genes_acc <- is[is$motif %in% low_variance_genes_acc,]$V1
  markers_genes_cl <- tf_markers %>% filter(cluster == i-1) %>% filter(avg_log2FC > 0)
  markers_genes_cl <- markers_genes_cl$gene
  markers_motif_cl <- motif_markers_enhd %>% filter(cluster == i-1) %>% filter(avg_log2FC > 0)
  markers_motif_cl <- markers_motif_cl$TF
  highlight <- markers_genes_cl[markers_genes_cl %in% markers_motif_cl]
  cluster_scatter_plots(CL, name = paste0("CL_", i-1),highlight_low_var = highlight, log = FALSE, folder_name = "MOTIF_ENHD", quartile = 4)
  if(i>1){
    row <- number_correlated_genes(CL, name = paste0("CL_", i-1), log = FALSE)
    number_genes_table_prom <- rbind(number_genes_table_prom,row)
  }
}


levels <- levels(motif_activity@meta.data[["predicted.id"]])
key = "_mvp.mean"
#key = "_vst.mean"
Sce_mean_motif = mean_table_motif_enhd[,paste0("CL", levels[1], key)]
Sce_mean_Expr = mean_table_expr[,paste0("CL", levels[1], key)]
for (i in levels[2:length(levels)]){
  
  Sce_mean_motif = cbind(Sce_mean_motif, mean_table_motif_enhd[,paste0("CL", i, key)])
  Sce_mean_Expr = cbind (Sce_mean_Expr, mean_table_expr[,paste0("CL", i, key)])
  
}
rownames(Sce_mean_motif) <- is$V1
rownames(Sce_mean_Expr) <- is$V1

correlation_table_final_enhd <- is_exprs[,1:2]
correlation_table_final_enhd[,2] <- NULL


s <- sapply(correlation_table_final_enhd$V1, function(x) cor.test(Sce_mean_motif[x,], Sce_mean_Expr[x,])[["estimate"]][["cor"]] )
correlation_table_final_enhd$PcorAvg_m <- s
s <- sapply(correlation_table_final_enhd$V1, function(x) cor.test(Sce_mean_motif[x,], Sce_mean_Expr[x,])[["p.value"]] )
correlation_table_final_enhd$PcorAvg_m_pvalue <- s
dim(correlation_table_final_enhd[correlation_table_final_enhd$PcorAvg_m > 0.5,])[1]/dim(correlation_table_final_enhd)[1]

top_variance_expr <- motif_activity@assays[["RNA"]]@meta.features[is_exprs$V1,c(2,4,7)] %>%  top_n(50, wt = vst.variance)
top_variance_expr <- rownames(top_variance_expr)
mean(correlation_table_final_enhd[correlation_table_final_enhd$V1 %in% top_variance_expr,]$PcorAvg_m)


mean(apply(mean_table_motif_enhd[,2:ncol(mean_table_motif_enhd)], 1 , sd))
mean(apply(mean_table_motif_enhd[,2:ncol(mean_table_motif_enhd)], 1 , mean))



s <- sapply(correlation_table_final$V1, function(x) cor.test(Sce_mean_motif[x,], Sce_mean_Expr[x,], method = "spearman")[["estimate"]][["rho"]] )
correlation_table_final$SrhoAvg_m <- s
s <- sapply(correlation_table_final$V1, function(x) cor.test(Sce_mean_motif[x,], Sce_mean_Expr[x,],method = "spearman")[["p.value"]] )
correlation_table_final$SrhoAvgAvg_m_pvalue <- s

correlation_table_final <- correlation_table_final[correlation_table_final$V1 %in% is_exprs$V1,]
dim(correlation_table_final[correlation_table_final$PcorAvg_m >0.5,])[1]/dim(correlation_table_final)[1]


ggplot(correlation_table_final, aes(x=PcorAvg_m, y = -log(PcorAvg_m_pvalue) )) + geom_point() + geom_point(data = correlation_table_final[correlation_table_final$PcorAvg_m < -0.5, ], aes(x=PcorAvg_m, y = -log(PcorAvg_m_pvalue)), colour = "blue")+
geom_point(data = correlation_table_final[correlation_table_final$PcorAvg_m > 0.5, ], aes(x=PcorAvg_m, y = -log(PcorAvg_m_pvalue)), colour = "red")

###############

gene <- "LEF1"
motif_activity@assays[["RNA"]]@data[gene,]
motif_activity@assays[["chromvar_enhd"]]@data[is[is$V1 == gene,]$motif,]

CC <- as.data.frame(t(rbind(motif_activity@assays[["chromvar_enhd"]]@data[is[is$V1 == gene,]$motif,],motif_activity@assays[["RNA"]]@data[gene,])))
colnames(CC) <- c("Activity", "Expression")
CC$type <- motif_activity@meta.data[["predicted.id"]]

CC <- CC %>%  filter(type = c("CD4 Naive", "CD4 Memory"))
ggplot(as.data.frame((CC)), aes(x=Activity, y = Expression)) + geom_point(aes(color= type)) +# geom_hline(yintercept = colMeans((CC))[2], linetype="dashed", color = "green", size=2) + 
  #geom_vline(xintercept = colMeans((CC))[1], linetype="dashed", color = "green", size=2) + 
  ggtitle(gene) +
  theme(plot.title = element_text(hjust = 0.5, size = 30), axis.title = element_text(hjust = 0.5, size = 30))
  #geom_point(data=as.data.frame((CC))[highlight_low_var, ], aes(x=Activity, y = Expression), colour="red", size=2)+
  #geom_text(aes(label = sub("-.*", "", rownames(CC))), vjust = -0.5)


CL <- cbind(mean_table_motif[,paste0("CL", i, key)], mean_table_expr[,paste0("CL", i, key)])
colnames(CL) <- c("Activity", "Expression")
rownames(CL) <- rownames(mean_table_expr)

cluster_scatter_plots(CL, name = paste0("CL_", i),highlight_low_var = NULL, log = FALSE, folder_name = "MOTIF_ENHD", quartile = 4)

#########

DefaultAssay(motif_activity) <- "peaks"
fragments <- CreateFragmentObject(
  path = "TMPDATA/10k_PBMC_Multiome_Controller/10k_PBMC_Multiome_nextgem_Chromium_Controller_atac_fragments.tsv.gz",
  cells = colnames(motif_activity), 
  validate.fragments = TRUE
)
Fragments(motif_activity) <- fragments


motif_activity <- Footprint(
  object = motif_activity,
  motif.name = c("CEBPD", "LEF1", "TCF7", "POU2F2", "FOS", "JUNB"),
  genome = BSgenome.Hsapiens.UCSC.hg38,
  in.peaks = 
)

PlotFootprint(motif_activity, features = c("CEBPD", "LEF1"))

DefaultAssay(motif_activity) <- "EnhdPeaks"

fragments <- CreateFragmentObject(
  path = "TMPDATA/10k_PBMC_Multiome_Controller/10k_PBMC_Multiome_nextgem_Chromium_Controller_atac_fragments.tsv.gz",
  cells = colnames(motif_activity), 
  validate.fragments = TRUE
)
Fragments(motif_activity) <- fragments


 
motif_activity <- Footprint(
  object = motif_activity,
  assay = "EnhdPeaks",
  motif.name = c("BATF"),
  genome = BSgenome.Hsapiens.UCSC.hg38,
  in.peaks = TRUE
)

motif_activity <- Footprint(
  object = motif_activity,
  assay = "peaks",
  motif.name = c("BATF"),
  genome = BSgenome.Hsapiens.UCSC.hg38,
  in.peaks = TRUE
)


PlotFootprint(motif_activity, features = c("JUNB", "FOS"),assay = "peaks", idents = c("CD4 Memory", "CD4 Naive")) / PlotFootprint(motif_activity, features = c("JUNB", "FOS"),assay = "EnhdPeaks", idents = c("CD4 Memory", "CD4 Naive"))

########## flanking footprint ######


plot.data <- GetFootprintData(
  object = motif_activity,
  features = c("JUNB", "FOS"),
  assay = "EnhdPeaks",
  group.by = "predicted.id"
)


motif.sizes <- GetMotifSize(
  object = motif_activity,
  features = c("JUNB", "FOS"),
  assay = "EnhdPeaks"
)
obs <- plot.data[plot.data$class == "Observed", ]
expect <- plot.data[plot.data$class == "Expected", ]

# flanks are motif edge to 50 bp each side
# add flank information (T/F)
base <- ceiling(motif.sizes / 2)
obs$flanks <- sapply(
  X = seq_len(length.out = nrow(x = obs)),
  FUN = function(x) {
    pos <- abs(obs[x, "position"])
    size <- base[[obs[x, "feature"]]]
    return((pos > size) & (pos < (size + 50)))
  })

flanks <- obs[obs$flanks,]
flanks <- group_by(.data = flanks, feature, group)
flankmeans <- summarize(.data = flanks, mn = mean(x = norm.value))
flankmeans <- flankmeans %>% pivot_wider(names_from = group, values_from = mn)
rownames(flankmeans) <- flankmeans$feature
# find top n groups for each feature
topmean <- top_n(x = flankmeans, n = label.top, wt = mn)

t(mean_table_expr %>%  filter(TF == "FOS")  %>%  select(contains("mvp")))

t(flankmeans["FOS",c(-1)])
CLTF <- as.data.frame(cbind(t(flankmeans["FOS",c(-1)]), t(mean_table_expr %>%  filter(TF == "FOS")  %>%  select(contains("mvp")))))
colnames(CLTF) <- c("flank_activity", "expression")
CLTF$type <- rownames(CLTF)

ggplot(as.data.frame((CLTF)), aes(x=flank_activity, y = expression)) + geom_point(aes(color = type)) +# geom_hline(yintercept = colMeans((CC))[2], linetype="dashed", color = "green", size=2) + 
  #geom_vline(xintercept = colMeans((CC))[1], linetype="dashed", color = "green", size=2) + 
  ggtitle(gene) + geom_text(aes(label = type), vjust = -0.5)
  theme(plot.title = element_text(hjust = 0.5, size = 30), axis.title = element_text(hjust = 0.5, size = 30))
#geom_point(data=as.data.frame((CC))[highlight_low_var, ], aes(x=Activity, y = Expression), colour="red", size=2)+
#geom_text(aes(label = sub("-.*", "", rownames(CC))), vjust = -0.5)


  
## Motif score for ENHD
mid <- motif_activity
DefaultAssay(mid) <- "EnhdPeaks"
mid[["RNA"]] <- NULL
mid[["chromvar_full"]] <- NULL
mid[["chromvar_enhd"]] <- NULL
mid[["peaks"]] <- NULL
TF_num <- length(is_exprs$V1)
register(SerialParam())
{
for (i in seq_along(is_exprs$V1)){
  
  TF <- is_exprs$V1[i]
mid <- Footprint(
  object = mid,
  motif.name = TF,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  in.peaks = TRUE
)
plot.data <- GetFootprintData(
  object = mid,
  features = TF,
  assay = "EnhdPeaks",
  group.by = "predicted.id"
)
motif.sizes <- GetMotifSize(
  object = mid,
  features = TF,
  assay = "EnhdPeaks"
)
obs <- plot.data[plot.data$class == "Observed", ]
base <- ceiling(motif.sizes / 2)
obs$flanks <- sapply(
  X = seq_len(length.out = nrow(x = obs)),
  FUN = function(x) {
    pos <- abs(obs[x, "position"])
    size <- base[[obs[x, "feature"]]]
    return((pos > size) & (pos < (size + 50)))
  })
flanks <- obs[obs$flanks,]
flanks <- group_by(.data = flanks, feature, group)
flankmeans <- summarize(.data = flanks, mn = mean(x = norm.value))
flankmeans <- flankmeans %>% pivot_wider(names_from = group, values_from = mn)
TF_footprint <-rbind(TF_footprint, flankmeans)
  
mid@assays[["EnhdPeaks"]]@positionEnrichment[[TF]]<- NULL
print(paste0(TF, ", ", i,"/",TF_num))
if (i %in% c(100,200,300,400)){
  write.csv(x = TF_footprint, file = "TMPResults/TF_footprint.csv")
}
}
write.csv(x = TF_footprint, file = "TMPResults/TF_footprint.csv")
}
rownames(TF_footprint) <- TF_footprint$feature


expr_foot <- mean_table_expr[is_exprs$V1,grep("mvp.mean$", names(mean_table_expr))]
colnames(expr_foot) <- colnames(TF_footprint)[-1]

for (i in levels(motif_activity@meta.data[["predicted.id"]])[-13]){
  
  CTF <- as.data.frame(cbind(TF_footprint[,i], expr_foot[,i]))
  colnames(CTF) <- c("Activity", "Expression")
  rownames(CTF) <- rownames(expr_foot)
  
  #low_variance_genes <- top_n(variance_ratio_table_expr[rownames(CL[CL[,"Expression"] > 0,]),c("TF",paste0("CL", i-1, "_vst.variance.standardized"))], -60)["TF"]$TF
  #low_variance_genes_acc <- top_n(variance_ratio_table_motif[,c("motif",paste0("CL", i-1, "_mvp.dispersion.scaled"))], -60)["motif"]$motif
  #low_variance_genes_acc <- top_n(variance_table_motif[,c("motif",paste0("CL", i-1, "_mvp.dispersion.scaled"))], -60)["motif"]$motif
  #low_variance_genes_acc <- is[is$motif %in% low_variance_genes_acc,]$V1
  #markers_genes_cl <- tf_markers %>% filter(cluster == i) %>% filter(avg_log2FC > 0)
  #markers_genes_cl <- markers_genes_cl$gene
  #markers_motif_cl <- motif_markers_enhd %>% filter(cluster == i) %>% filter(avg_log2FC > 0)
  #markers_motif_cl <- markers_motif_cl$TF
  #highlight <- markers_genes_cl[markers_genes_cl %in% markers_motif_cl]
  cluster_scatter_plots(CTF, name = paste0( i),highlight_low_var = NULL, log = FALSE, folder_name = "FOOT_ENHD", quartile = 3, X= "Footprint Score (Enhancer)", Y = "Expresssion")
}

foot_variance <- data.frame(variance = rowVars(as.matrix(TF_footprint[,-1])), row.names = is_exprs$V1)
variance_TF <- rownames(top_n(foot_variance, n = 50))

EXPR_variance <- data.frame(variance = rowVars(as.matrix(expr_foot[,])), row.names = is_exprs$V1)
variance_EXPR <- rownames(top_n(foot_variance, n = 50))

for (i in variance_TF){
  
  CTF <- as.data.frame(t(rbind(TF_footprint[i,c(-1,-14)], expr_foot[i,-13])))
  colnames(CTF) <- c("Activity", "Expression")
  
  #low_variance_genes <- top_n(variance_ratio_table_expr[rownames(CL[CL[,"Expression"] > 0,]),c("TF",paste0("CL", i-1, "_vst.variance.standardized"))], -60)["TF"]$TF
  #low_variance_genes_acc <- top_n(variance_ratio_table_motif[,c("motif",paste0("CL", i-1, "_mvp.dispersion.scaled"))], -60)["motif"]$motif
  #low_variance_genes_acc <- top_n(variance_table_motif[,c("motif",paste0("CL", i-1, "_mvp.dispersion.scaled"))], -60)["motif"]$motif
  #low_variance_genes_acc <- is[is$motif %in% low_variance_genes_acc,]$V1
  #markers_genes_cl <- tf_markers %>% filter(cluster == i) %>% filter(avg_log2FC > 0)
  #markers_genes_cl <- markers_genes_cl$gene
  #markers_motif_cl <- motif_markers_enhd %>% filter(cluster == i) %>% filter(avg_log2FC > 0)
  #markers_motif_cl <- markers_motif_cl$TF
  #highlight <- markers_genes_cl[markers_genes_cl %in% markers_motif_cl]
  cluster_scatter_plots(CTF, name = paste0(i),highlight_low_var = NULL, log = FALSE, folder_name = "FOOT_ENHD_GENE", quartile = FALSE, X= "Footprint Score (Enhancer)", Y = "Expresssion")
}


for (i in c("TCF7", "POU2F2", "LEF1")){
  
  CTF <- as.data.frame(t(rbind(TF_footprint[i,c(-1,-14)], expr_foot[i,-13])))
  colnames(CTF) <- c("Activity", "Expression")
  
  #low_variance_genes <- top_n(variance_ratio_table_expr[rownames(CL[CL[,"Expression"] > 0,]),c("TF",paste0("CL", i-1, "_vst.variance.standardized"))], -60)["TF"]$TF
  #low_variance_genes_acc <- top_n(variance_ratio_table_motif[,c("motif",paste0("CL", i-1, "_mvp.dispersion.scaled"))], -60)["motif"]$motif
  #low_variance_genes_acc <- top_n(variance_table_motif[,c("motif",paste0("CL", i-1, "_mvp.dispersion.scaled"))], -60)["motif"]$motif
  #low_variance_genes_acc <- is[is$motif %in% low_variance_genes_acc,]$V1
  #markers_genes_cl <- tf_markers %>% filter(cluster == i) %>% filter(avg_log2FC > 0)
  #markers_genes_cl <- markers_genes_cl$gene
  #markers_motif_cl <- motif_markers_enhd %>% filter(cluster == i) %>% filter(avg_log2FC > 0)
  #markers_motif_cl <- markers_motif_cl$TF
  #highlight <- markers_genes_cl[markers_genes_cl %in% markers_motif_cl]
  cluster_scatter_plots(CTF, name = paste0(i),highlight_low_var = NULL, log = FALSE, folder_name = "FOOT_ENHD_GENE", quartile = FALSE, X= "Footprint Score (Enhancer)", Y = "Expresssion")
}


correlation_table_FOOT <- is_exprs[,1:2]
correlation_table_FOOT[,2] <- NULL
rownames(correlation_table_FOOT) <- is_exprs$V1
s <- sapply(correlation_table_FOOT$V1, function(x) {cor.test(as.numeric(TF_footprint[x,c(-1,-14)]), as.numeric(expr_foot[x,-13]))[["estimate"]][["cor"]] })
correlation_table_FOOT$PcorAvg_m <- s
s <- sapply(correlation_table_FOOT$V1, function(x) cor.test(as.numeric(TF_footprint[x,c(-1,-14)]), as.numeric(expr_foot[x,-13]))[["p.value"]] )
correlation_table_FOOT$PcorAvg_m_pvalue <- s

dim(correlation_table_FOOT[correlation_table_FOOT$PcorAvg_m > 0.5,])[1]/dim(correlation_table_FOOT)[1]

top_variance_expr <- motif_activity@assays[["RNA"]]@meta.features[is_exprs$V1,c(2,4,7)] %>%  top_n(50, wt = mvp.dispersion)
top_variance_expr <- rownames(top_variance_expr)
mean(correlation_table_FOOT[correlation_table_FOOT$V1 %in% top_variance_expr,]$PcorAvg_m)


median <- quantile(TF_footprint[,2:ncol(TF_footprint)][,2])[4]
CL_0 <- CL_0[CL_0[,2] > median,]

mean(apply(TF_footprint[,2:ncol(TF_footprint)], 1 , sd))
mean(apply(TF_footprint[,2:ncol(TF_footprint)], 1 , mean))

TF_enhd_avg <- data.frame(row.names = is_exprs$V1, SD = apply(TF_footprint[,2:ncol(TF_footprint)], 1 , sd), mean = apply(TF_footprint[,2:ncol(TF_footprint)], 1 , mean))

mean(apply(TF_footprint[,2:ncol(TF_footprint)], 1 , sd))
mean(apply(TF_footprint[,2:ncol(TF_footprint)], 1 , mean))

q1 <- ggplot(TF_enhd_avg, aes(x=mean, y = SD)) + geom_point()



ggplot(correlation_table_FOOT, aes(x=PcorAvg_m, y = -log(PcorAvg_m_pvalue) )) + geom_point() + geom_point(data = correlation_table_FOOT[correlation_table_FOOT$PcorAvg_m < -0.5, ], aes(x=PcorAvg_m, y = -log(PcorAvg_m_pvalue)), colour = "blue")+
  geom_point(data = correlation_table_FOOT[correlation_table_FOOT$PcorAvg_m > 0.5, ], aes(x=PcorAvg_m, y = -log(PcorAvg_m_pvalue)), colour = "red") 

ggplot(correlation_table_final, aes(x=PcorAvg_m, y = -log(PcorAvg_m_pvalue) )) + geom_point() + geom_point(data = correlation_table_final[correlation_table_final$PcorAvg_m < -0.5, ], aes(x=PcorAvg_m, y = -log(PcorAvg_m_pvalue)), colour = "blue")+
  geom_point(data = correlation_table_final[correlation_table_final$PcorAvg_m > 0.5, ], aes(x=PcorAvg_m, y = -log(PcorAvg_m_pvalue)), colour = "red")


ggplot(correlation_table_FOOT, aes(x=PcorAvg_m, y = -log(PcorAvg_m_pvalue) )) + geom_point() +
  geom_point(data = correlation_table_FOOT[variance_EXPR, ], aes(x=PcorAvg_m, y = -log(PcorAvg_m_pvalue)), colour = "red") 



vln_foot_enhd <- TF_footprint[,2:ncol(TF_footprint)]
colnames(vln_foot_enhd) <- levels
gather(vln_foot_enhd[-13], key = "CellType", value = "Expression")  %>% ggplot( aes(x = CellType, y = Expression, fill = CellType)) +
  geom_violin() +
  geom_jitter(color = "black", size = 0.5, alpha = 0.5) + geom_hline(yintercept = 1) + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(title = "Gene Expression across Cell Types", 
       x = "Cell Types", 
       y = "Expression Level")


###

TF_footprint_prom <- NULL

mid <- motif_activity
DefaultAssay(mid) <- "peaks"
mid[["RNA"]] <- NULL
mid[["chromvar_full"]] <- NULL
mid[["chromvar_enhd"]] <- NULL
mid[["EnhdPeaks"]] <- NULL
TF_num <- length(is_exprs$V1)
register(SerialParam())

if (ftp == TRUE){
  for (i in seq_along(is_exprs$V1)){
    
    TF <- is_exprs$V1[i]
    mid <- Footprint(
      object = mid,
      motif.name = TF,
      genome = BSgenome.Hsapiens.UCSC.hg38,
      in.peaks = TRUE
    )
    plot.data <- GetFootprintData(
      object = mid,
      features = TF,
      assay = "peaks",
      group.by = "predicted.id"
    )
    motif.sizes <- GetMotifSize(
      object = mid,
      features = TF,
      assay = "peaks"
    )
    obs <- plot.data[plot.data$class == "Observed", ]
    base <- ceiling(motif.sizes / 2)
    obs$flanks <- sapply(
      X = seq_len(length.out = nrow(x = obs)),
      FUN = function(x) {
        pos <- abs(obs[x, "position"])
        size <- base[[obs[x, "feature"]]]
        return((pos > size) & (pos < (size + 50)))
      })
    flanks <- obs[obs$flanks,]
    flanks <- group_by(.data = flanks, feature, group)
    flankmeans <- summarize(.data = flanks, mn = mean(x = norm.value))
    flankmeans <- flankmeans %>% pivot_wider(names_from = group, values_from = mn)
    TF_footprint_prom <-rbind(TF_footprint_prom, flankmeans)
    
    mid@assays[["peaks"]]@positionEnrichment[[TF]]<- NULL
    print(paste0(TF, ", ", i,"/",TF_num))
    if (i %in% c(100,200,300,400)){
      write.csv(x = TF_footprint_prom, file = "TMPResults/TF_footprint_prom.csv")
    }
  }
  write.csv(x = TF_footprint_prom, file = "TMPResults/TF_footprint_prom.csv")
}
TF_footprint_prom <- read.csv(file = "TMPResults/TF_footprint_prom.csv", row.names = 1, )


rownames(TF_footprint_prom) <- TF_footprint_prom$feature
colnames(TF_footprint_prom) [2:15] <- levels

expr_foot <- mean_table_expr[is_exprs$V1,grep("mvp.mean$", names(mean_table_expr))]
colnames(expr_foot) <- colnames(TF_footprint_prom)[-1]


for (i in levels(motif_activity@meta.data[["predicted.id"]])[-13]){
  
  CTF <- as.data.frame(cbind(TF_footprint_prom[,i], expr_foot[,i]))
  colnames(CTF) <- c("Activity", "Expression")
  rownames(CTF) <- rownames(expr_foot)
  
  #low_variance_genes <- top_n(variance_ratio_table_expr[rownames(CL[CL[,"Expression"] > 0,]),c("TF",paste0("CL", i-1, "_vst.variance.standardized"))], -60)["TF"]$TF
  #low_variance_genes_acc <- top_n(variance_ratio_table_motif[,c("motif",paste0("CL", i-1, "_mvp.dispersion.scaled"))], -60)["motif"]$motif
  #low_variance_genes_acc <- top_n(variance_table_motif[,c("motif",paste0("CL", i-1, "_mvp.dispersion.scaled"))], -60)["motif"]$motif
  #low_variance_genes_acc <- is[is$motif %in% low_variance_genes_acc,]$V1
  #markers_genes_cl <- tf_markers %>% filter(cluster == i) %>% filter(avg_log2FC > 0)
  #markers_genes_cl <- markers_genes_cl$gene
  #markers_motif_cl <- motif_markers_enhd %>% filter(cluster == i) %>% filter(avg_log2FC > 0)
  #markers_motif_cl <- markers_motif_cl$TF
  #highlight <- markers_genes_cl[markers_genes_cl %in% markers_motif_cl]
  cluster_scatter_plots(CTF, name = paste0( i),highlight_low_var = NULL, log = FALSE, folder_name = "FOOT_PROM", quartile = FALSE, X= "Footprint Score (Promoter)", Y = "Expresssion")
}

foot_variance <- data.frame(variance = rowVars(as.matrix(TF_footprint_prom[,-1])), row.names = is_exprs$V1)
variance_TF <- rownames(top_n(foot_variance, n = 50))

EXPR_variance <- data.frame(variance = rowVars(as.matrix(expr_foot[,])), row.names = is_exprs$V1)
variance_EXPR <- rownames(top_n(foot_variance, n = 50))

for (i in rownames(correlation_table_FOOT_prom[correlation_table_FOOT_prom$PcorAvg_m > 0.5, ])){
  
  CTF <- as.data.frame(t(rbind(TF_footprint_prom[i,c(-1,-14)], expr_foot[i,-13])))
  colnames(CTF) <- c("Activity", "Expression")
  
  #low_variance_genes <- top_n(variance_ratio_table_expr[rownames(CL[CL[,"Expression"] > 0,]),c("TF",paste0("CL", i-1, "_vst.variance.standardized"))], -60)["TF"]$TF
  #low_variance_genes_acc <- top_n(variance_ratio_table_motif[,c("motif",paste0("CL", i-1, "_mvp.dispersion.scaled"))], -60)["motif"]$motif
  #low_variance_genes_acc <- top_n(variance_table_motif[,c("motif",paste0("CL", i-1, "_mvp.dispersion.scaled"))], -60)["motif"]$motif
  #low_variance_genes_acc <- is[is$motif %in% low_variance_genes_acc,]$V1
  #markers_genes_cl <- tf_markers %>% filter(cluster == i) %>% filter(avg_log2FC > 0)
  #markers_genes_cl <- markers_genes_cl$gene
  #markers_motif_cl <- motif_markers_enhd %>% filter(cluster == i) %>% filter(avg_log2FC > 0)
  #markers_motif_cl <- markers_motif_cl$TF
  #highlight <- markers_genes_cl[markers_genes_cl %in% markers_motif_cl]
  cluster_scatter_plots(CTF, name = paste0(i),highlight_low_var = NULL, log = FALSE, folder_name = "FOOT_PROM_GENE", quartile = FALSE, X= "Footprint Score (Promoter)", Y = "Expresssion")
}



for (i in c("TCF7", "POU2F2", "LEF1")){
  
  CTF <- as.data.frame(t(rbind(TF_footprint_prom[i,c(-1,-14)], expr_foot[i,-13])))
  colnames(CTF) <- c("Activity", "Expression")
  
  #low_variance_genes <- top_n(variance_ratio_table_expr[rownames(CL[CL[,"Expression"] > 0,]),c("TF",paste0("CL", i-1, "_vst.variance.standardized"))], -60)["TF"]$TF
  #low_variance_genes_acc <- top_n(variance_ratio_table_motif[,c("motif",paste0("CL", i-1, "_mvp.dispersion.scaled"))], -60)["motif"]$motif
  #low_variance_genes_acc <- top_n(variance_table_motif[,c("motif",paste0("CL", i-1, "_mvp.dispersion.scaled"))], -60)["motif"]$motif
  #low_variance_genes_acc <- is[is$motif %in% low_variance_genes_acc,]$V1
  #markers_genes_cl <- tf_markers %>% filter(cluster == i) %>% filter(avg_log2FC > 0)
  #markers_genes_cl <- markers_genes_cl$gene
  #markers_motif_cl <- motif_markers_enhd %>% filter(cluster == i) %>% filter(avg_log2FC > 0)
  #markers_motif_cl <- markers_motif_cl$TF
  #highlight <- markers_genes_cl[markers_genes_cl %in% markers_motif_cl]
  cluster_scatter_plots(CTF, name = paste0(i),highlight_low_var = NULL, log = FALSE, folder_name = "FOOT_PROM_GENE", quartile = FALSE, X= "Footprint Score (Promoter)", Y = "Expresssion")
}

correlation_table_FOOT_prom <- is_exprs[,1:2]
correlation_table_FOOT_prom[,2] <- NULL
rownames(correlation_table_FOOT_prom) <- is_exprs$V1
s <- sapply(correlation_table_FOOT_prom$V1, function(x) {cor.test(as.numeric(TF_footprint_prom[x,c(-1,-14)]), as.numeric(expr_foot[x,-13]))[["estimate"]][["cor"]] })
correlation_table_FOOT_prom$PcorAvg_m <- s
s <- sapply(correlation_table_FOOT_prom$V1, function(x) cor.test(as.numeric(TF_footprint_prom[x,c(-1,-14)]), as.numeric(expr_foot[x,-13]))[["p.value"]] )
correlation_table_FOOT_prom$PcorAvg_m_pvalue <- s

dim(correlation_table_FOOT_prom[correlation_table_FOOT_prom$PcorAvg_m > 0.5,])[1]/dim(correlation_table_FOOT_prom)[1]

top_variance_expr <- motif_activity@assays[["RNA"]]@meta.features[is_exprs$V1,c(2,4,7)] %>%  top_n(50, wt = mvp.dispersion)
top_variance_expr <- rownames(top_variance_expr)
mean(correlation_table_FOOT_prom[correlation_table_FOOT_prom$V1 %in% top_variance_expr,]$PcorAvg_m)

TF_prom_avg <- data.frame(row.names = is_exprs$V1, SD = apply(TF_footprint_prom[,2:ncol(TF_footprint_prom)], 1 , sd), mean = apply(TF_footprint_prom[,2:ncol(TF_footprint_prom)], 1 , mean))

mean(apply(TF_footprint_prom[,2:ncol(TF_footprint_prom)], 1 , sd))
mean(apply(TF_footprint_prom[,2:ncol(TF_footprint_prom)], 1 , mean))

q2 <- ggplot(TF_prom_avg, aes(x=mean, y = SD)) + geom_point()
plot_grid(q1, q2, ncol = 2, align = "hv")



ggplot(correlation_table_FOOT_prom, aes(x=PcorAvg_m, y = -log(PcorAvg_m_pvalue) )) + geom_point() + geom_point(data = correlation_table_FOOT_prom[correlation_table_FOOT_prom$PcorAvg_m < -0.5, ], aes(x=PcorAvg_m, y = -log(PcorAvg_m_pvalue)), colour = "blue")+
  geom_point(data = correlation_table_FOOT_prom[correlation_table_FOOT_prom$PcorAvg_m > 0.5, ], aes(x=PcorAvg_m, y = -log(PcorAvg_m_pvalue)), colour = "red") 

ggplot(correlation_table_final, aes(x=PcorAvg_m, y = -log(PcorAvg_m_pvalue) )) + geom_point() + geom_point(data = correlation_table_final[correlation_table_final$PcorAvg_m < -0.5, ], aes(x=PcorAvg_m, y = -log(PcorAvg_m_pvalue)), colour = "blue")+
  geom_point(data = correlation_table_final[correlation_table_final$PcorAvg_m > 0.5, ], aes(x=PcorAvg_m, y = -log(PcorAvg_m_pvalue)), colour = "red")


ggplot(correlation_table_FOOT_prom, aes(x=PcorAvg_m, y = -log(PcorAvg_m_pvalue) )) + geom_point() +
  geom_point(data = correlation_table_FOOT_prom[variance_TF, ], aes(x=PcorAvg_m, y = -log(PcorAvg_m_pvalue)), colour = "red") 

rownames(correlation_table_FOOT_prom[correlation_table_FOOT_prom$PcorAvg_m > 0.5, ])

correlation_table_FOOT_prom[correlation_table_FOOT_prom$PcorAvg_m < -0.5, ]



vln_foot_prom <- TF_footprint_prom[,2:ncol(TF_footprint_prom)]
colnames(vln_foot_prom) <- levels
gather(vln_foot_prom[-13], key = "CellType", value = "Expression")  %>% ggplot( aes(x = CellType, y = Expression, fill = CellType)) +
  geom_violin() +
  geom_jitter(color = "black", size = 0.5, alpha = 0.5) + geom_hline(yintercept = 1) + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(title = "Gene Expression across Cell Types", 
       x = "Cell Types", 
       y = "Expression Level")




###############

if (!(dir.exists(paste0("TMPResults/IMAGES/","MISC")))){
  dir.create(paste0("TMPResults/IMAGES/","MISC"))
}

gene<- "BATF"
p1 <- VlnPlot(motif_activity, features = is_exprs[is_exprs$V1 == gene,]$motif, assay = "chromvar_enhd", idents = levels[-13], pt.size = 0.5) + geom_boxplot()  + NoLegend()  +labs(title = "", y = "Motif Enrichment" )
#p3 <- VlnPlot(motif_activity, features = is_exprs[is_exprs$V1 == gene,]$motif, assay = "chromvar_full", idents = levels[-13]) + geom_boxplot() + geom_hline(yintercept = 0, linetype= "dashed", size =1) + NoLegend()

p2 <- VlnPlot(motif_activity, features = gene, assay = "RNA", idents = levels[-13]) + geom_boxplot() + theme(axis.text.x = element_blank()) +labs( x = "" )
p2 + p1
ggsave(path = paste0("TMPResults/IMAGES/", "MISC"), filename = paste0(gene,"_violin_plot.png"), width = 1080, height = 1080, units= "px",scale = 3.5)

PlotFootprint(motif_activity, features = gene ,assay = "EnhdPeaks", idents = c("CD4 Memory", "CD4 Naive")) / PlotFootprint(motif_activity, features = gene ,assay = "peaks", idents = c("CD4 Memory", "CD4 Naive"))
ggsave(path = paste0("TMPResults/IMAGES/", "MISC"), filename = paste0(gene,"_footprint.png"), width = 800, height = 1080, units= "px",scale = 3.5)

PlotFootprint(motif_activity, features = gene ,assay = "EnhdPeaks", idents = c("CD8 effector", "CD8 Naive")) / PlotFootprint(motif_activity, features = gene ,assay = "peaks", idents = c("CD8 effector", "CD8 Naive"))
PlotFootprint(motif_activity, features = gene ,assay = "EnhdPeaks", idents = levels[-13]) / PlotFootprint(motif_activity, features = gene ,assay = "peaks", idents = levels[-13])
PlotFootprint(motif_activity, features = gene ,assay = "EnhdPeaks", idents = c("B cell progenitor", "pre-B cell")) / PlotFootprint(motif_activity, features = gene ,assay = "peaks", idents = c("B cell progenitor", "pre-B cell"))
PlotFootprint(motif_activity, features = gene ,assay = "EnhdPeaks", idents = c("CD14+ Monocytes", "CD16+ Monocytes", "Dendritic cell")) / PlotFootprint(motif_activity, features = gene ,assay = "peaks", idents = c("CD14+ Monocytes", "CD16+ Monocytes", "Dendritic cell"))
PlotFootprint(motif_activity, features = gene ,assay = "EnhdPeaks", idents = c("CD4 Memory", "CD4 Naive","CD8 effector", "CD8 Naive")) / PlotFootprint(motif_activity, features = gene ,assay = "peaks", idents = c("CD4 Memory", "CD4 Naive","CD8 effector", "CD8 Naive"))


MotifPlot(
  object = motif_activity,
  motifs = "MA0476.1",
)

PlotFootprint(motif_activity, features = c("JUNB", "FOS"),assay = "peaks", idents = c("CD4 Memory", "CD4 Naive")) / PlotFootprint(motif_activity, features = c("JUNB", "FOS"),assay = "EnhdPeaks", idents = c("CD4 Memory", "CD4 Naive"))
PlotFootprint(motif_activity, features = c("JUNB", "FOS"),assay = "peaks", idents = levels[-13])  / PlotFootprint(motif_activity, features = c("JUNB", "FOS"),assay = "EnhdPeaks", idents = levels[-13])

######## functions #####

cluster_scatter_plots <- function (CL_0,name,  log = TRUE, division = "mean", highlight_low_var = NULL, folder_name = "PROM", quartile = FALSE, X="Activity", Y="Expression"){
  
  if (!(dir.exists(paste0("TMPResults/IMAGES/",folder_name)))){
    dir.create(paste0("TMPResults/IMAGES/",folder_name))
  }
  if (quartile != FALSE) {
    median <- quantile(CL_0[,2])[quartile]
    CL_0 <- CL_0[CL_0[,2] > median,]
  }
  if(log){
    if (division == "mean"){
      #png(paste0("TMPResults/IMAGES/",name,"_gene_scatter_plot.png"),width = 1080, height = 1080)
      plot <-ggplot(as.data.frame(CL_0), aes(x=1, y = log(2))) + geom_point() + geom_hline(yintercept = colMeans(log(CL_0)[log(CL_0)[,2] > -Inf,])[2], linetype="dashed", color = "green", size=2) + 
        geom_vline(xintercept = colMeans((CL_0))[1], linetype="dashed", color = "green", size=2)+ ggtitle(name) +
        theme(plot.title = element_text(hjust = 0.5, size = 30), axis.title = element_text(hjust = 0.5, size = 30))+
        geom_text(aes(label = sub("-.*", "", rownames(CL_0))), vjust = -0.5)+
        geom_point(data=as.data.frame(log(CL_0))[highlight_low_var, ], aes(x=Activity, y = Expression), colour="red", size=2) + theme_minimal() #+ xlim(-150,3)
      #geom_point(data=as.data.frame(log(CL_0))[high_cor_genes, ], aes(x=Activity, y = Expression), colour="red", size=2)
      #dev.off()
      ggsave(path = paste0("TMPResults/IMAGES/", folder_name,"/"), filename = paste0(name,"_gene_scatter_plot.pdf"), width = 1080, height = 1080, units= "px",scale = 3.5)
    }
    
  }
  if(!log){
    if(division == "mean"){
      #png(paste0("TMPResults/IMAGES/",name,"_gene_scatter_plot.png"),width = 1080, height = 1080)
      ggplot(as.data.frame((CL_0)), aes(x=Activity, y = Expression)) + geom_point() + geom_hline(yintercept = colMeans((CL_0))[2], linetype="dashed", color = "green", size=2) + 
        geom_vline(xintercept = colMeans((CL_0))[1], linetype="dashed", color = "green", size=2)+ ggtitle(name) +
        theme(plot.title = element_text(hjust = 0.5, size = 30), axis.title = element_text(hjust = 0.5, size = 30))+
        geom_point(data=as.data.frame((CL_0))[highlight_low_var, ], aes(x=Activity, y = Expression), colour="red", size=2)+
        geom_text(aes(label = sub("-.*", "", rownames(CL_0))), vjust = -0.5, size =40) + labs(x = X, y = Y) +  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "gray90"),
                                                                                                           panel.grid.minor = element_line(color = "gray90"), axis.line = element_line(color = "gray80", linewidth = 0.5))
      #geom_point(data=as.data.frame((CL_0))[high_cor_genes, ], aes(x=Activity, y = Expression), colour="red", size=2)
      ggsave(path = paste0("TMPResults/IMAGES/", folder_name,"/"), filename = paste0(name,"_gene_scatter_plot.pdf"), width = 1080, height = 1080, units= "px",scale = 3.5)
      #dev.off()
    }
  }
}


GetMotifSize <- function(
    object,
    features,
    assay = NULL
) {
  positionEnrichment <- GetAssayData(
    object = object,
    assay = assay,
    slot = "positionEnrichment"
  )
  sizes <- c()
  for (i in features) {
    motif <- positionEnrichment[[i]]["motif", ]
    sizes <- c(sizes, sum(motif))
  }
  names(x = sizes) <- features
  return(sizes)
}

