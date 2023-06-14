# Note: scDesign3 requre R version V.4.2.0
# R lib path: /home/gayan/R/x86_64-pc-linux-gnu-library/4.3
# Independent Rscript and python script since loading R packages are different
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("devtools", quietly = TRUE)) 
  install.packages("devtools")
if (!requireNamespace("Seurat", quietly = TRUE)) {
  install.packages('Seurat')
}
if (!requireNamespace("scDesign3", quietly = TRUE)) {
devtools::install_github("SONGDONGYUAN1994/scDesign3")
}
if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
BiocManager::install("SingleCellExperiment")
}

# Load libraries
cat(sprintf("[scReadSim] Loading R packages...\n"))
suppressPackageStartupMessages(library(scDesign3))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(SingleCellExperiment))

# Specify directories
samplename_RNA <- "10X_RNA_chr1_3073253_4526737.gene.countmatrix"
samplename_ATAC <- "10X_ATAC_chr1_4194444_4399104.peak.countmatrix"
directory <- "/home/gayan/Projects/scATAC_Simulator/results/Revision_test_scMultiOmics_20230611"
out_directory <- directory
# n_cell_new <- 10000
# n_cores <- 10
# scMultiOmics_runSyntheticCount(samplename_RNA, samplename_ATAC, directory, out_directory, n_cores=10)

# Functions
get_cluster_seurat <- function(data_mat, n_cluster,
                               res = 0.5){
  submat <- data_mat
  if(is.null(rownames(submat)))
    rownames(submat) <- 1:nrow(submat)
    count_seurat <- CreateSeuratObject(counts = submat, verbose = FALSE)
    ### normalization
    count_seurat <- NormalizeData(count_seurat,verbose = FALSE)
    ### select highly variable genes
    count_seurat <- FindVariableFeatures(count_seurat, selection.method = "vst", nfeatures = 4000,verbose = FALSE)
    ### scale the data
    count_seurat <- ScaleData(count_seurat,verbose = FALSE)
    ### PCA
    count_seurat <- RunPCA(count_seurat,
                            features = VariableFeatures(object = count_seurat),
                            verbose = F)
    ### clustering
    dims = 1:min(10, min(ncol(data_mat), nrow(data_mat))-1)
    count_seurat <- FindNeighbors(count_seurat, dims = dims)
  # If not specify target cluster number
    if (n_cluster == "default"){
        cat("No target cell cluster number detected.\n")        
        cat("Louvain clustering with default resolution 0.5...\n")        
        count_seurat <- FindClusters(count_seurat, resolution = res, verbose = FALSE)
        ### results
        cluster_predicted <- as.integer(Idents(count_seurat))
        colnames(submat) <- as.character(cluster_predicted)
        # rogue_scores <- check_cluster_quality(submat, unique(cluster_predicted), platform)
        cat(sprintf("%s cell clusters generated.\n", length(unique(cluster_predicted))))
    } else{
        n_cluster <- as.integer(n_cluster)
        cat(sprintf("Detected target cell cluster number %s\n", n_cluster))        
        cat("Generating target clusters with Louvain clustering through binary search...\n")        
        min_res <- 0
        max_res <- 3
        max_steps <- 20
        this_step <- 0
        while (this_step < max_steps) {
            cat('step', this_step,'\n')
            this_res <- min_res + ((max_res-min_res)/2)
            data_louvain <- Seurat::FindClusters(count_seurat, resolution = this_res, verbose = FALSE)
            this_clusters <- length(unique(Seurat::Idents(data_louvain)))
            if (this_clusters > n_cluster) {
            max_res <- this_res
            } else if (this_clusters < n_cluster) {
            min_res <- this_res
            } else {
            break
            }
            this_step <- this_step + 1    
        }
        Seurat_louvain <- Seurat::Idents(data_louvain)
        ### results
        cluster_predicted <- as.integer(Seurat_louvain)
        colnames(submat) <- as.character(cluster_predicted)
        cat(sprintf("%s cell clusters generated.\n", this_clusters))
    }
  return(list(clustering_result = cluster_predicted))
}

scMultiOmics_runSyntheticCount <- function(samplename_RNA, samplename_ATAC, directory, out_directory, n_cell_new="default", celllabel_file="default", n_cluster="default", n_cores=1){
    cat(sprintf("[scReadSim] Reading RNA count matrix %s.txt...\n", samplename_RNA))
    # Load in scRNA-seq UMI count matrix
    RNA_mat <- read.table(sprintf("%s/%s.txt", directory, samplename_RNA), sep="\t",header = FALSE,  row.names = 1) %>% as.matrix()
    cat(sprintf("[scReadSim] Reading ATAC count matrix %s.txt...\n", samplename_ATAC))
    # Load in scATAC-seq count matrix 
    ATAC_mat <- read.table(sprintf("%s/%s.txt", directory, samplename_ATAC), sep="\t",header = FALSE,  row.names = 1) %>% as.matrix()
    
    # Extract dimensionality
    n_RNA <- nrow(RNA_mat)
    n_ATAC <- nrow(ATAC_mat)

    # Clustering based on scRNA-seq
    matrix_num <- RNA_mat
    count_pergene_vec <- rowSums(matrix_num)
    matrix_num_nonzero <- matrix_num[count_pergene_vec>0,]
    if (celllabel_file == "default"){
        cat("[scReadSim] No cell label file detected. Louvain clustering before simulation...\n")
        clustering_result <- get_cluster_seurat(matrix_num_nonzero, n_cluster=n_cluster)
        colnames(matrix_num) <- clustering_result$clustering_result
        full_cell_label <- clustering_result$clustering_result
        cat("[scReadSim] Writing out Louvain clustering result...\n")
        cat("[scReadSim] Created:\n")
        cat(sprintf("[scReadSim] Real count matrix's Louvain clustering file: %s/%s.LouvainClusterResults.txt\n", out_directory, samplename_RNA))
        write.table(full_cell_label, sprintf("%s/%s.LouvainClusterResults.txt", out_directory, samplename_RNA), sep="\n", row.names = FALSE,col.names = FALSE)
    } else {
        cat(sprintf("[scReadSim] Loading cell label file %s...\n", celllabel_file))
        full_cell_label <- unlist(read.table(celllabel_file, header=FALSE))
        if (length(full_cell_label) == ncol(matrix_num)){
            colnames(matrix_num) <- full_cell_label
        } else {
            stop("[scReadSim] Number of cell labels differs from the cell number contained in the count matrix!\n ")
        }
    }

    # Generate sce object
    mat <- rbind(RNA_mat, ATAC_mat)
    sce <- SingleCellExperiment(list(counts = mat))
    colData(sce)$cell_type <- as.factor(full_cell_label)
    n_cell_old <- ncol(mat)
    if (n_cell_new == "default"){
        n_cell_new <- n_cell_old
    }

    # Simulate
    set.seed(2022)
    new_sce <- scdesign3(sce = sce,
                        assay_use = "counts",
                        celltype = "cell_type",
                        pseudotime = NULL,
                        spatial = NULL, 
                        other_covariates = NULL, 
                        ncell = n_cell_new, 
                        mu_formula = "cell_type", 
                        sigma_formula = "cell_type", 
                        family_use = "nb", 
                        n_cores = n_cores, 
                        usebam = FALSE, 
                        corr_formula = "cell_type", 
                        copula = "gaussian", 
                        fastmvn = FALSE, 
                        DT = TRUE, 
                        pseudo_obs = FALSE, 
                        family_set = "gauss", 
                        nonnegative = TRUE, 
                        nonzerovar = FALSE, 
                        return_model = FALSE, 
                        parallelization = "mclapply", 
                        trace = FALSE
                        )

    # Write out data files
    synthetic_RNA_mat <- new_sce$new_count[1:n_RNA, ]
    synthetic_ATAC_mat <- new_sce$new_count[(1+n_RNA):nrow(mat), ]
    if (n_cell_new == n_cell_old){
      synthetic_cell_label <- full_cell_label
    } else {
      synthetic_cell_label <- new_sce$new_covariate$cell_type
    }
    cat(sprintf("[scReadSim] Writing out synthetic cell labels to %s...\n", out_directory))
    write.table(synthetic_cell_label, sprintf("%s/scMultiOmics.scDesign3Simulated.CellTypeLabel.txt", out_directory), row.names = FALSE,col.names = FALSE, quote = FALSE)
    cat(sprintf("[scReadSim] Writing out synthetic count matrices to %s...\n", out_directory))
    write.table(synthetic_RNA_mat, sprintf("%s/%s.scMultiOmics.scDesign3Simulated.RNA.txt", out_directory, samplename_RNA), sep="\t", row.names = TRUE,col.names = FALSE)
    write.table(synthetic_ATAC_mat, sprintf("%s/%s.scMultiOmics.scDesign3Simulated.ATAC.txt", out_directory, samplename_ATAC), sep="\t", row.names = TRUE,col.names = FALSE)
    cat("[scReadSim] Done.\n")
}



