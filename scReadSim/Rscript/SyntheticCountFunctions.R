library(Matrix)
library(Rsubread)
library(pscl)
library(parallel)
library(MASS)
library(tidyverse)
library(ROGUE)
library(Seurat)

### Functions
fit_marginals_new <- function(x, marginal = c('auto_choose', 'zinb', 'nb', 'poisson'),
                          pval_cutoff = 0.05, epsilon = 1e-5,
                          jitter = TRUE, DT = TRUE){
  p <- nrow(x)
  n <- ncol(x)

  marginal <- match.arg(marginal)
  if(marginal == 'auto_choose'){
    params <- t(apply(x, 1, function(gene){
      m <- mean(gene)
      v <- var(gene)
      if(m >= v){
        mle_Poisson <- glm(gene ~ 1, family = poisson)
        tryCatch({
          mle_ZIP <- zeroinfl(gene ~ 1|1, dist = 'poisson')
          chisq_val <- 2 * (logLik(mle_ZIP) - logLik(mle_Poisson))
          pvalue <- as.numeric(1 - pchisq(chisq_val, 1))
          if(pvalue < pval_cutoff)
            c(plogis(mle_ZIP$coefficients$zero), Inf, exp(mle_ZIP$coefficients$count))
          else
            c(0.0, Inf, m)
        },
        error = function(cond){
          c(0.0, Inf, m)})
      }else{
        mle_NB <- glm.nb(gene ~ 1)
        if(min(gene) > 0)
          c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
        else
          tryCatch({
            mle_ZINB <- zeroinfl(gene ~ 1|1, dist = 'negbin')
            chisq_val <- 2 * (logLik(mle_ZINB) - logLik(mle_NB))
            pvalue <- as.numeric(1 - pchisq(chisq_val, 1))
            if(pvalue < pval_cutoff)
              c(plogis(mle_ZINB$coefficients$zero), mle_ZINB$theta, exp(mle_ZINB$coefficients$count))
            else
              c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
          },
          error = function(cond){
            c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
          })
      }
    }))
  }else if(marginal == 'zinb'){
    params <- t(apply(x, 1, function(gene){
      m <- mean(gene)
      v <- var(gene)
      if(m >= v)
      {
        mle_Poisson <- glm(gene ~ 1, family = poisson)
        tryCatch({
          mle_ZIP <- zeroinfl(gene ~ 1|1, dist = 'poisson')
          chisq_val <- 2 * (logLik(mle_ZIP) - logLik(mle_Poisson))
          pvalue <- as.numeric(1 - pchisq(chisq_val, 1))
          if(pvalue < pval_cutoff)
            c(plogis(mle_ZIP$coefficients$zero), Inf, exp(mle_ZIP$coefficients$count))
          else
            c(0.0, Inf, m)
        },
        error = function(cond){
          c(0.0, Inf, m)})
      }
      else
      {
        if(min(gene) > 0)
        {
          mle_NB <- glm.nb(gene ~ 1)
          c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
        }
        else
          tryCatch({
            mle_ZINB <- zeroinfl(gene ~ 1|1, dist = 'negbin')
            c(plogis(mle_ZINB$coefficients$zero), mle_ZINB$theta, exp(mle_ZINB$coefficients$count))
          },
          error = function(cond){
            mle_NB <- glm.nb(gene ~ 1)
            c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
          })
      }
    }))
  }else if(marginal == 'nb'){
    params <- t(apply(x, 1, function(gene){
      m <- mean(gene)
      v <- var(gene)
      if(m >= v){
        c(0.0, Inf, m)
      }else{
        mle_NB <- glm.nb(gene ~ 1)
        c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
      }
    }))
  }else if(marginal == 'poisson'){
    params <- t(apply(x, 1, function(gene){
      c(0.0, Inf, mean(gene))
    }))
  }

  if(DT){
    u <- t(sapply(1:p, function(iter){
      param <- params[iter, ]
      gene <- unlist(x[iter, ])
      prob0 <- param[1]
      u1 <- prob0 + (1 - prob0) * pnbinom(gene, size = param[2], mu = param[3])
      u2 <- (prob0 + (1 - prob0) * pnbinom(gene - 1, size = param[2], mu = param[3])) *
        as.integer(gene > 0)
      if(jitter)
        v <- runif(n)
      else
        v <- rep(0.5, n)
      r <- u1 * v + u2 * (1 - v)
      idx_adjust <- which(1-r < epsilon)
      r[idx_adjust] <- r[idx_adjust] - epsilon
      idx_adjust <- which(r < epsilon)
      r[idx_adjust] <- r[idx_adjust] + epsilon

      r
    }))
  }else{
    u <- NULL
  }

  return(list(params = params, u = u))
}

fit_Gaussian_copula_new <- function(x, marginal = c('auto_choose', 'zinb', 'nb', 'poisson'),
                                jitter = TRUE, zp_cutoff = 0.8,
                                min_nonzero_num = 2){
  marginal <- match.arg(marginal)
  n <- ncol(x)
  p <- nrow(x)

  gene_zero_prop <- apply(x, 1, function(y){
    sum(y < 1e-5) / n
  })

  gene_sel1 <- which(gene_zero_prop < zp_cutoff)
  gene_sel2 <- which(gene_zero_prop < 1.0 - min_nonzero_num/n &
                       gene_zero_prop >= zp_cutoff)
  gene_sel3 <- (1:p)[-c(gene_sel1, gene_sel2)]

  if(length(gene_sel1) > 0){
    marginal_result1 <- fit_marginals_new(x[gene_sel1, , drop = FALSE], marginal, jitter = jitter, DT = TRUE)
    quantile_normal <- qnorm(marginal_result1$u)
    cov_mat <- cor(t(quantile_normal))
  }else{
    cov_mat = NULL
    marginal_result1 = NULL
  }

  if(length(gene_sel2) > 0){
    marginal_result2 <- fit_marginals_new(x[gene_sel2, , drop = FALSE], marginal, DT = FALSE)
  }else{
    marginal_result2 = NULL
  }
  return(list(cov_mat = cov_mat, marginal_param1 = marginal_result1$params,
              marginal_param2 = marginal_result2$params,
              gene_sel1 = gene_sel1, gene_sel2 = gene_sel2, gene_sel3 = gene_sel3,
              zp_cutoff = zp_cutoff, min_nonzero_num = min_nonzero_num,
              sim_method = 'copula', n_cell = n, n_read = sum(x)))
}

fit_wo_copula_new <- function(x, marginal = c('auto_choose', 'zinb', 'nb', 'poisson'),
                          jitter = TRUE, min_nonzero_num = 2){
  marginal <- match.arg(marginal)
  n <- ncol(x)
  p <- nrow(x)

  gene_zero_prop <- apply(x, 1, function(y){
    sum(y < 1e-5) / n
  })

  gene_sel1 <- which(gene_zero_prop < 1.0 - min_nonzero_num/n)
  gene_sel2 <- (1:p)[-gene_sel1]

  if(length(gene_sel1) > 0){
    marginal_result1 <- fit_marginals_new(x[gene_sel1, ], marginal, jitter = jitter, DT = FALSE)
  }else{
    marginal_result1 = NULL
  }

  return(list(marginal_param1 = marginal_result1$params,
              gene_sel1 = gene_sel1, gene_sel2 = gene_sel2,
              min_nonzero_num = min_nonzero_num, sim_method = 'ind',
              n_cell = n, n_read = sum(x)))
}


fit_model_scDesign2_new <- function(data_mat, cell_type_sel, sim_method = c('copula', 'ind'),
                                marginal = c('auto_choose', 'zinb', 'nb', 'poisson'),
                                jitter = TRUE, zp_cutoff = 0.8,
                                min_nonzero_num = 2, ncores = 1){
  sim_method <- match.arg(sim_method)
  marginal <- match.arg(marginal)

  if(sum(abs(data_mat - round(data_mat))) > 1e-5){
    warning('The entries in the input matrix are not integers. Rounding is performed.')
    data_mat <- round(data_mat)
  }

  if(sim_method == 'copula'){
    param <- mclapply(1:length(cell_type_sel), function(iter){
      fit_Gaussian_copula_new(data_mat[, colnames(data_mat) == cell_type_sel[iter]], marginal,
                          jitter = jitter, zp_cutoff = zp_cutoff,
                          min_nonzero_num = min_nonzero_num)
    }, mc.cores = ncores)
  }else if(sim_method == 'ind'){
    param <- mclapply(1:length(cell_type_sel), function(iter){
      fit_wo_copula_new(data_mat[, colnames(data_mat) == cell_type_sel[iter]], marginal,
                    jitter = jitter,
                    min_nonzero_num = min_nonzero_num)
    }, mc.cores = ncores)
  }

  names(param) <- cell_type_sel
  return(param)
}

# function of getting a submatrix of selected cell types ----------------------------------
get_submat <- function(data_mat, cell_type_sel){
  if(is.null(colnames(data_mat))){
    data_mat
  }else{
    data_mat[, colnames(data_mat) %in% cell_type_sel]
  }
}
# function of using ROGUE scores to check cluster qualities -------------------------------
check_cluster_quality <- function(data_mat, cell_type_sel,
                                  platform = c("full-length", "UMI")){
  platform <- match.arg(platform)
  
  data_mat_sel <- get_submat(data_mat, cell_type_sel)
  
  expr <- data_mat_sel
  rogue.res <- rogue(expr, labels = colnames(expr),
                     samples = rep(1, ncol(expr)), platform = platform)
  unlist(rogue.res)
}
# function of using Seurat to cluster -----------------------------------------------------
get_cluster_seurat <- function(data_mat, platform = c("full-length", "UMI"),
                               res = 0.5){
  platform <- match.arg(platform)
  submat <- data_mat
  
  if(is.null(rownames(submat)))
    rownames(submat) <- 1:nrow(submat)
  count_seurat <- CreateSeuratObject(counts = submat)
  ### normalization
  count_seurat <- NormalizeData(count_seurat)
  ### select highly variable genes
  count_seurat <- FindVariableFeatures(count_seurat, selection.method = "vst", nfeatures = 4000)
  ### scale the data
  count_seurat <- ScaleData(count_seurat)
  ### PCA
  count_seurat <- RunPCA(count_seurat,
                         features = VariableFeatures(object = count_seurat),
                         verbose = F)
  ### clustering
  dims = 1:min(10, min(ncol(data_mat), nrow(data_mat))-1)
  count_seurat <- FindNeighbors(count_seurat, dims = dims)
  count_seurat <- FindClusters(count_seurat, resolution = res)
  ### results
  cluster_predicted <- as.integer(Idents(count_seurat))
  
  colnames(submat) <- as.character(cluster_predicted)
  # rogue_scores <- check_cluster_quality(submat, unique(cluster_predicted), platform)
  return(list(clustering_result = cluster_predicted))
}

# Calculate marginal count for each region
# args <- commandArgs(trailingOnly = TRUE) # Extract all arguments from Linux input
# samplename <- args[1]
# sample.format <- args[2]
# directory <- args[3]
# out.directory <- args[4]

# Test
# samplename <- "e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1.countmatrix"
# # samplename <- "e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1.COMPLE.countmatrix"
# sample.format <- "txt"
# directory <- "/home/gayan/Projects/scATAC_Simulator/package_development/package_results/20220116_e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1_NONINPUT_withCluster"
# out.directory <- "/home/gayan/Projects/scATAC_Simulator/package_development/package_results/20220116_e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1_NONINPUT_withCluster"
# dir.create(out.directory)

# samplename <- '10X_ATAC_chr1_4194444_4399104.countmatrix' 
# sample_format<-'txt'
# directory <- '/home/gayan/Projects/scATAC_Simulator/package_development/package_results/20220310_10X_scATACseq_NONINPUT'
# out_directory <- '/home/gayan/Projects/scATAC_Simulator/package_development/package_results/20220310_10X_scATACseq_NONINPUT'


## Read in count matrix
scATAC_runSyntheticCount <- function(samplename, directory, out_directory, n_cell_new="default", total_count_new="default", celllabel_file="default"){
  cat(sprintf("Reading count matrix %s.txt...\n", samplename))
  count_matrix <- read.table(sprintf("%s/%s.txt", directory, samplename), sep="\t",header = FALSE)
  matrix_num <- data.matrix(count_matrix[,2:ncol(count_matrix)])
  count_pergene_vec <- rowSums(ceiling(matrix_num/2))
  write.table(count_pergene_vec, sprintf("%s/%s.real.nPairsRegionmargional.txt",out_directory, samplename), row.names = FALSE,col.names = FALSE)
  matrix_num_nonzero <- matrix_num[count_pergene_vec>0,]
  ## Clustering
  if (celllabel_file == "default"){
    cat("No cell label file detected. Louvain clustering before simulation...\n")
    set.seed(2022)
    clustering_result <- get_cluster_seurat(matrix_num_nonzero, 'UMI')
    colnames(matrix_num) <- clustering_result$clustering_result
    write.table(clustering_result$clustering_result, sprintf("%s/%s.LouvainClusterResults.txt", out_directory, samplename), sep="\n", row.names = FALSE,col.names = FALSE)
  } else {
    cat(sprintf("Loading cell label file %s...\n", celllabel_file))
    clustering_result <- unlist(read.table(celllabel_file, header=FALSE))
    if (length(clustering_result) == ncol(count_matrix)){
    colnames(matrix_num) <- rep("Cluster1", ncol(matrix_num))
    } else {
      stop("Number of cell labels differs from the cell number contained in the count matrix!\n ")
    }
  }

## Use scDesign2 for training countmatrix
cat("Model fitting...\n")
n_cell_old <- ncol(matrix_num)
total_count_old <- sum(matrix_num)
if (n_cell_new == "default"){
  n_cell_new <- n_cell_old
}
if (total_count_new == "default"){
  total_count_new <- total_count_old
}
cell_type_sel <- unique(colnames(matrix_num))
cell_type_prop <- table(colnames(matrix_num))[cell_type_sel]
copula_result <- fit_model_scDesign2_new(matrix_num, cell_type_sel, sim_method = 'copula',
                                       ncores = length(cell_type_sel))
cat("Generating synthetic count matrix...\n")
cat(sprintf("Amount of synthetic cell: %s\n", n_cell_new))
cat(sprintf("Amount of (expected) sequencing depth: %s\n", total_count_new))

# simu_matrix <- scDesign2::simulate_count_scDesign2(copula_result, n_cell_new, sim_method = 'copula',
#                                                cell_type_prop = cell_type_prop)
simu_matrix <- scDesign2::simulate_count_scDesign2(copula_result, 
                           total_count_old = total_count_old,
                           n_cell_old = n_cell_old,
                           total_count_new = total_count_new,
                           n_cell_new = n_cell_new,
                           cell_type_prop = cell_type_prop,
                           reseq_method = 'mean_scale', cell_sample = TRUE)

write.table(colnames(simu_matrix), sprintf("%s/%s.scDesign2Simulated.CellTypeLabel.txt", out_directory, samplename), row.names = FALSE,col.names = FALSE)
rownames(simu_matrix) <- count_matrix[,1]
cat(sprintf("Writing out synthetic count matrix %s to %s...\n", out_directory, out_directory))
write.table(simu_matrix, sprintf("%s/%s.scDesign2Simulated.txt", out_directory, samplename), sep="\t", row.names = TRUE,col.names = TRUE)
simu_count_pergene_vec <- rowSums(ceiling(simu_matrix/2))
write.table(simu_count_pergene_vec, sprintf("%s/%s.scDesign2Simulated.nPairsRegionmargional.txt", out_directory, samplename), row.names = FALSE,col.names = FALSE)
cat("Done.\n")
}


# samplename <- '10X_RNA_chr1_3073253_4526737.countmatrix' 
# sample_format<-'txt'
# directory <- '/home/gayan/Projects/scATAC_Simulator/package_development/package_results/20220310_10X_scRNAseq_NONINPUT'
# out_directory <- '/home/gayan/Projects/scATAC_Simulator/package_development/package_results/20220310_10X_scRNAseq_NONINPUT'

scRNA_runSyntheticCount <- function(samplename, directory, out_directory, UMI_modeling=FALSE, UMIsamplename="UMI_countmat", n_cell_new="default", total_count_new="default", celllabel_file="default"){
  if (UMI_modeling == FALSE){
    scRNA_runSyntheticCount_Read (samplename, directory, out_directory, n_cell_new, total_count_new, celllabel_file)
  } else{
    scRNA_runSyntheticCount_UMIandRead(samplename, UMIsamplename, directory, out_directory, n_cell_new, total_count_new, celllabel_file)
  }
}


scRNA_runSyntheticCount_UMIandRead <- function(samplename, UMIsamplename, directory, out_directory, n_cell_new="default", total_count_new="default", celllabel_file="default"){
  cat(sprintf("Reading count matrix %s.txt...\n", samplename))
  count_matrix <- read.table(sprintf("%s/%s.txt", directory, samplename), sep="\t",header = FALSE)
  UMI_count_matrix <- read.table(sprintf("%s/%s.txt", directory, UMIsamplename), sep="\t",header = FALSE)
  matrix_num_UMI <- data.matrix(UMI_count_matrix[,-1])
  # Louvain clustering on read count matrix
  matrix_num <- data.matrix(count_matrix[,2:ncol(count_matrix)])
  count_pergene_vec <- rowSums(matrix_num)
  write.table(count_pergene_vec, sprintf("%s/%s.real.nReadRegionmargional.txt",out_directory, samplename), row.names = FALSE,col.names = FALSE)
  matrix_num_nonzero <- matrix_num[count_pergene_vec>0,]
  ## Clustering
  if (celllabel_file == "default"){
    cat("No cell label file detected. Louvain clustering before simulation...\n")
    set.seed(2022)
    clustering_result <- get_cluster_seurat(matrix_num_nonzero, 'UMI')
    colnames(matrix_num) <- clustering_result$clustering_result
    colnames(matrix_num_UMI) <- clustering_result$clustering_result

    write.table(clustering_result$clustering_result, sprintf("%s/%s.LouvainClusterResults.txt", out_directory, samplename), sep="\n", row.names = FALSE,col.names = FALSE)
  } else {
    cat(sprintf("Loading cell label file %s...\n", celllabel_file))
    clustering_result <- unlist(read.table(celllabel_file, header=FALSE))
    if (length(clustering_result) == ncol(count_matrix)){
    colnames(matrix_num) <- rep("Cluster1", ncol(matrix_num))
    colnames(matrix_num_UMI) <- rep("Cluster1", ncol(matrix_num_UMI))
    } else {
      stop("Number of cell labels differs from the cell number contained in the count matrix! \n")
    }
  }

## Use scDesign2 for training countmatrix
cat("Model fitting...\n")
n_cell_old <- ncol(matrix_num)
total_count_old <- sum(matrix_num)
total_count_UMI_old <- sum(matrix_num_UMI)
if (n_cell_new == "default"){
  n_cell_new <- n_cell_old
}
if (total_count_new == "default"){
  total_count_new <- total_count_old
  total_count_UMI_new <- total_count_UMI_old
} else {
  total_count_UMI_new <- round(total_count_new/total_count_old * total_count_UMI_old)
}
cell_type_sel <- unique(colnames(matrix_num))
cell_type_prop <- table(colnames(matrix_num))[cell_type_sel]
copula_result <- fit_model_scDesign2_new(matrix_num, cell_type_sel, sim_method = 'copula',
                                       ncores = length(cell_type_sel))
copula_result_UMI <- fit_model_scDesign2_new(matrix_num_UMI, cell_type_sel, sim_method = 'copula',
                                       ncores = length(cell_type_sel))
cat("Generating synthetic count matrix...\n")
cat(sprintf("Amount of synthetic cell: %s\n", n_cell_new))
cat(sprintf("Amount of (expected) sequencing depth: %s\n", total_count_new))
# simu_matrix <- scDesign2::simulate_count_scDesign2(copula_result, n_cell_new, sim_method = 'copula',
#                                                cell_type_prop = cell_type_prop)
simu_matrix <- scDesign2::simulate_count_scDesign2(copula_result, 
                           total_count_old = total_count_old,
                           n_cell_old = n_cell_old,
                           total_count_new = total_count_new,
                           n_cell_new = n_cell_new,
                           cell_type_prop = cell_type_prop,
                           reseq_method = 'mean_scale', cell_sample = TRUE)
simu_matrix_UMI <- scDesign2::simulate_count_scDesign2(copula_result_UMI, 
                           total_count_old = total_count_UMI_old,
                           n_cell_old = n_cell_old,
                           total_count_new = total_count_UMI_new,
                           n_cell_new = n_cell_new,
                           cell_type_prop = cell_type_prop,
                           reseq_method = 'mean_scale', cell_sample = TRUE)

rownames(simu_matrix) <- count_matrix[,1]
rownames(simu_matrix_UMI) <- count_matrix[,1]

write.table(colnames(simu_matrix), sprintf("%s/%s.scDesign2Simulated.CellTypeLabel.txt", out_directory, samplename), row.names = FALSE,col.names = FALSE)
cat(sprintf("Writing out synthetic count matrix %s to %s...\n", out_directory, out_directory))
write.table(simu_matrix, sprintf("%s/%s.scDesign2Simulated.txt", out_directory, samplename), sep="\t", row.names = TRUE,col.names = TRUE)
write.table(simu_matrix_UMI, sprintf("%s/%s.scDesign2Simulated.txt", out_directory, UMIsamplename), sep="\t", row.names = TRUE,col.names = TRUE)
simu_count_pergene_vec <- rowSums(simu_matrix)
write.table(simu_count_pergene_vec, sprintf("%s/%s.scDesign2Simulated.nReadRegionmargional.txt", out_directory, samplename), row.names = FALSE,col.names = FALSE)
cat("Done.\n")
}


scRNA_runSyntheticCount_Read <- function(samplename, directory, out_directory, n_cell_new="default", total_count_new="default", celllabel_file="default"){
  cat(sprintf("Reading count matrix %s.txt...\n", samplename))
  count_matrix <- read.table(sprintf("%s/%s.txt", directory, samplename), sep="\t",header = FALSE)
  matrix_num <- data.matrix(count_matrix[,2:ncol(count_matrix)])
count_pergene_vec <- rowSums(matrix_num)
  write.table(count_pergene_vec, sprintf("%s/%s.real.nReadRegionmargional.txt",out_directory, samplename), row.names = FALSE,col.names = FALSE)
  matrix_num_nonzero <- matrix_num[count_pergene_vec>0,]
  ## Clustering
  if (celllabel_file == "default"){
    cat("No cell label file detected. Louvain clustering before simulation...\n")
    set.seed(2022)
    clustering_result <- get_cluster_seurat(matrix_num_nonzero, 'UMI')
    colnames(matrix_num) <- clustering_result$clustering_result
    write.table(clustering_result$clustering_result, sprintf("%s/%s.LouvainClusterResults.txt", out_directory, samplename), sep="\n", row.names = FALSE,col.names = FALSE)
  } else {
    cat(sprintf("Loading cell label file %s...\n", celllabel_file))
    clustering_result <- unlist(read.table(celllabel_file, header=FALSE))
    if (length(clustering_result) == ncol(count_matrix)){
    colnames(matrix_num) <- rep("Cluster1", ncol(matrix_num))
    } else {
      stop("Number of cell labels differs from the cell number contained in the count matrix! \n")
    }
  }

## Use scDesign2 for training countmatrix
cat("Model fitting...\n")
n_cell_old <- ncol(matrix_num)
total_count_old <- sum(matrix_num)
if (n_cell_new == "default"){
  n_cell_new <- n_cell_old
}
if (total_count_new == "default"){
  total_count_new <- total_count_old
}
cell_type_sel <- unique(colnames(matrix_num))
cell_type_prop <- table(colnames(matrix_num))[cell_type_sel]
copula_result <- fit_model_scDesign2_new(matrix_num, cell_type_sel, sim_method = 'copula',
                                       ncores = length(cell_type_sel))
cat("Generating synthetic count matrix...\n")
cat(sprintf("Amount of synthetic cell: %s\n", n_cell_new))
cat(sprintf("Amount of (expected) sequencing depth: %s\n", total_count_new))
# simu_matrix <- scDesign2::simulate_count_scDesign2(copula_result, n_cell_new, sim_method = 'copula',
#                                                cell_type_prop = cell_type_prop)
simu_matrix <- scDesign2::simulate_count_scDesign2(copula_result, 
                           total_count_old = total_count_old,
                           n_cell_old = n_cell_old,
                           total_count_new = total_count_new,
                           n_cell_new = n_cell_new,
                           cell_type_prop = cell_type_prop,
                           reseq_method = 'mean_scale', cell_sample = TRUE)
rownames(simu_matrix) <- count_matrix[,1]
write.table(colnames(simu_matrix), sprintf("%s/%s.scDesign2Simulated.CellTypeLabel.txt", out_directory, samplename), row.names = FALSE,col.names = FALSE)
cat(sprintf("Writing out synthetic count matrix %s to %s...\n", out_directory, out_directory))
write.table(simu_matrix, sprintf("%s/%s.scDesign2Simulated.txt", out_directory, samplename), sep="\t", row.names = TRUE,col.names = TRUE)
simu_count_pergene_vec <- rowSums(simu_matrix)
write.table(simu_count_pergene_vec, sprintf("%s/%s.scDesign2Simulated.nReadRegionmargional.txt", out_directory, samplename), row.names = FALSE,col.names = FALSE)
cat("Done.\n")
}











