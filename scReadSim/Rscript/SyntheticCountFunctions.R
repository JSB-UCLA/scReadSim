library(Matrix)
library(Rsubread)
library(pscl)
library(parallel)
library(MASS)
library(tidyverse)
library(ROGUE)
library(Seurat)

### Functions
# fit_marginals_zinb <- function(x, mc = FALSE, ncores = 1, DT = TRUE, jitter_ind = TRUE, epsilon = 1e-5){
#   p <- nrow(x)
#   n <- ncol(x)
  
#   if (mc){
#     params_list <- mclapply(1:p, function(iter){
#       gene <- unlist(x[iter,])
#       if(min(gene) > 0){mle_NB <- glm.nb(gene ~ 1)
#       c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
#       }else{
#         mle_ZINB <- zeroinfl(gene ~ 1|1, dist = 'negbin')
#         c(plogis(mle_ZINB$coefficients$zero), mle_ZINB$theta, exp(mle_ZINB$coefficients$count))
#       }
#     }, mc.cores = ncores)
#     params <- matrix(unlist(params_list, use.names = FALSE), nrow = p, byrow = TRUE)
    
#   } else{
#     params <- t(apply(x, 1, function(gene){
#       if(min(gene) > 0)
#       {mle_NB <- glm.nb(gene ~ 1)
#       c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
#       }
#       else
#       {
#         mle_ZINB <- zeroinfl(gene ~ 1|1, dist = 'negbin')
#         c(plogis(mle_ZINB$coefficients$zero), mle_ZINB$theta, exp(mle_ZINB$coefficients$count))
#       }
#     }))
#   }
#   if(DT){
#     u_list <- mclapply(1:p, function(iter){
#       param <- params[iter, ]
#       gene <- unlist(x[iter, ])
#       prob0 <- param[1]
#       u1 <- prob0 + (1 - prob0) * pnbinom(gene, size = param[2], mu = param[3])
#       u2 <- (prob0 + (1 - prob0) * pnbinom(gene - 1, size = param[2], mu = param[3])) *
#         as.integer(gene > 0)
#       if(jitter_ind) {
#         v <- runif(n)
#       } else 
#         v <- rep(0.5, n)
#       r <- u1 * v + u2 * (1 - v)
#       idx_adjust <- which(1-r < epsilon)
#       r[idx_adjust] <- r[idx_adjust] - epsilon
#       idx_adjust <- which(r < epsilon)
#       r[idx_adjust] <- r[idx_adjust] + epsilon

#       r
#     }, mc.cores = ncores)
#     u <- matrix(unlist(u_list, use.names = FALSE), nrow = p, byrow = TRUE)
#     quantile_normal <- qnorm(u)
#     cov_mat <- cor(t(quantile_normal))
#   }else{
#     cov_mat <- NULL
#   }
#   return(list(params = params, cov_mat = cov_mat))

# }

# fit_Gaussian_copula <- function(x, mc=FALSE, ncores = 1, jitter = TRUE, zp_cutoff = 0.8,
#                                 min_nonzero_num = 2){
#   n <- ncol(x)
#   p <- nrow(x)

#   gene_zero_prop <- apply(x, 1, function(y){
#     sum(y < 1e-5) / n
#   })

#   gene_sel1 <- which(gene_zero_prop < zp_cutoff)
#   gene_sel2 <- which(gene_zero_prop < 1.0 - min_nonzero_num/n &
#                        gene_zero_prop >= zp_cutoff)
#   gene_sel3 <- (1:p)[-c(gene_sel1, gene_sel2)]

#   if(length(gene_sel1) > 0){
#     # print("Part 1 model fitting")
#     marginal_result1 <- fit_marginals_zinb(x[gene_sel1, , drop = FALSE], mc = mc, ncores = ncores, jitter = jitter, DT = TRUE)
#     # print("End Part 1 model fitting")
#     cov_mat <- marginal_result1$cov_mat
#   }else{
#     cov_mat = NULL
#     marginal_result1 = NULL
#   }

#   if(length(gene_sel2) > 0){
#     # print("Part 2 model fitting")
#     marginal_result2 <- fit_marginals_zinb(x[gene_sel2, , drop = FALSE], mc = mc, ncores = ncores, DT = FALSE)
#     # print("End Part 2 model fitting")
#   }else{
#     marginal_result2 = NULL
#   }
#   return(list(cov_mat = cov_mat, marginal_param1 = marginal_result1$params,
#               marginal_param2 = marginal_result2$params,
#               gene_sel1 = gene_sel1, gene_sel2 = gene_sel2, gene_sel3 = gene_sel3,
#               zp_cutoff = zp_cutoff, min_nonzero_num = min_nonzero_num,
#               sim_method = 'copula', n_cell = n, n_read = sum(x)))
# }

# fit_model_scDesign2 <- function(data_mat, cell_type_sel, jitter = TRUE, zp_cutoff = 0.8,
#                                 min_nonzero_num = 2, ncores = 1){
#   if(sum(abs(data_mat - round(data_mat))) > 1e-5){
#     warning('The entries in the input matrix are not integers. Rounding is performed.')
#     data_mat <- round(data_mat)
#   }
#     param <- mclapply(1:length(cell_type_sel), function(iter){
#       print(sprintf("CELL TYPE %s RUNING", iter))
#       return(fit_Gaussian_copula(data_mat[, colnames(data_mat) == cell_type_sel[iter]],
#                           jitter = jitter, zp_cutoff = zp_cutoff,
#                           min_nonzero_num = min_nonzero_num))
#       print(sprintf("CELL TYPE %s DONE.", iter))
#     }, mc.cores = ncores)
#   names(param) <- cell_type_sel
#   return(param)
# }

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


# simulate_count_copula <- function(copula_result, n = 100){
#   p1 <- length(copula_result$gene_sel1)
#   if(p1 > 0){
#     result1 <- mvrnorm(n = n, mu = rep(0.0, p1), Sigma = copula_result$cov_mat)
#     result1 <- matrix(result1, nrow = n)
#     result2 <- apply(result1, 2, pnorm)
#     result2 <- matrix(result2, nrow = n)
#   }
#   p2 <- length(copula_result$gene_sel2)
#   if(p1 > 0){
#     result31 <- t(sapply(1:p1, function(iter){
#       param <- copula_result$marginal_param1[iter, ]
#     qnbinom(pmax(0.0, result2[, iter] - param[1]) / (1-param[1]),
#               size = param[2], mu = param[3])
#     }))
#   }
#   if(p2 > 0){
#     result32 <- t(sapply(1:p2, function(iter){
#       param <- copula_result$marginal_param2[iter, ]
#       rbinom(n, 1, 1-param[1]) * rnbinom(n, size = param[2], mu = param[3])
#     }))
#   }

#   result <- matrix(0, nrow = p1 + p2 + length(copula_result$gene_sel3), ncol = n)
#   if(p1 > 0){
#     result[copula_result$gene_sel1, ] <- result31
#   }
#   if(p2 > 0){
#     result[copula_result$gene_sel2, ] <- result32
#   }
#   result
# }


# simulate_count_scDesign2 <- function(model_params, n_cell_new, cell_type_prop = 1,
#                                      total_count_new = NULL, total_count_old = NULL,
#                                      n_cell_old = NULL, cell_sample = FALSE){
#   n_cell_vec <- sapply(model_params, function(x) x$n_cell)
#   n_read_vec <- sapply(model_params, function(x) x$n_read)

#   # if(is.null(total_count_new)) total_count_new <- sum(n_read_vec)
#   # if(is.null(n_cell_new))      n_cell_new      <- sum(n_cell_vec)
#   # if(is.null(cell_type_prop))  cell_type_prop  <- n_cell_vec
#   if(is.null(total_count_old)) total_count_old <- sum(n_read_vec)
#   if(is.null(n_cell_old))      n_cell_old      <- sum(n_cell_vec)

#   if(is.null(total_count_new)) reseq_method <- 'mean_scale'

#   if(length(model_params)!=length(cell_type_prop)){
#     stop('Cell type proportion should have the same length as the number of models.')
#   }

#   n_cell_type <- length(cell_type_prop)
#   if(cell_sample == TRUE){
#     n_cell_each <- as.numeric(rmultinom(1, size = n_cell_new, prob = cell_type_prop))
#   }else{
#     cell_type_prop <- cell_type_prop / sum(cell_type_prop)
#     n_cell_each <- round(cell_type_prop * n_cell_new)
#     if(sum(n_cell_each) != n_cell_new){
#       idx <- sample(n_cell_type, size = 1)
#       n_cell_each[idx] <- n_cell_each[idx] + n_cell_new - sum(n_cell_each)
#     }
#   }

#   p <- length(model_params[[1]]$gene_sel1) + length(model_params[[1]]$gene_sel2) +
#     length(model_params[[1]]$gene_sel3)
#   new_count <- matrix(0, nrow = p, ncol = n_cell_new)
#     n_cell_each
#     if(is.null(total_count_new)){
#       r <- rep(1, n_cell_type)
#     }else if(length(total_count_new) == 1){
#       r <- rep(total_count_new / sum((total_count_old / n_cell_old) * n_cell_each),
#                n_cell_type)
#     }else{
#       r <- (total_count_new / n_cell_new) / (total_count_old / n_cell_old)
#     }
#     for(iter in 1:n_cell_type)
#       if(n_cell_each[iter] > 0){
#         ulim <- sum(n_cell_each[1:iter])
#         llim <- ulim - n_cell_each[iter] + 1
#         params_new <- model_params[[iter]]
#         params_new$marginal_param1[, 3] <- params_new$marginal_param1[, 3] * r[iter]
#         params_new$marginal_param2[, 3] <- params_new$marginal_param2[, 3] * r[iter]
#         new_count[, llim:ulim] <- simulate_count_copula(params_new, n = n_cell_each[iter])
#       }
#     if(is.null(names(model_params))){
#       colnames(new_count) <- unlist(lapply(1:n_cell_type, function(x){rep(x, n_cell_each[x])}))
#     }else{
#       colnames(new_count) <- unlist(lapply(1:n_cell_type, function(x){
#         rep(names(model_params)[x], n_cell_each[x])}))
#     }
#     return(new_count)
# }

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
scATAC_runSyntheticCount <- function(samplename, sample_format, directory, out_directory, cluster_prestep=TRUE){
  cat(sprintf("Reading count matrix %s.%s...\n", samplename, sample_format))
  count_matrix <- read.table(sprintf("%s/%s.%s", directory, samplename, sample_format), sep="\t",header = FALSE)
  matrix_num <- data.matrix(count_matrix[,2:ncol(count_matrix)])
  count_pergene_vec <- rowSums(ceiling(matrix_num/2))
  write.table(count_pergene_vec, sprintf("%s/%s.real.nPairsRegionmargional.txt",out_directory, samplename), row.names = FALSE,col.names = FALSE)
    matrix_num_nonzero <- matrix_num[count_pergene_vec>0,]
  ## Clustering
  if (cluster_prestep == TRUE){
    cat("Louvain Clustering Before Simulation...\n")
    set.seed(2022)
    clustering_result <- get_cluster_seurat(matrix_num_nonzero, 'UMI')
    # print(clustering_result$rogue_scores)
    # print(mean(clustering_result$rogue_scores))
    colnames(matrix_num) <- clustering_result$clustering_result
    write.table(clustering_result$clustering_result, sprintf("%s/%s.LouvainClusterResults.%s", out_directory, samplename, sample_format), sep="\n", row.names = FALSE,col.names = FALSE)
  } else {
    colnames(matrix_num) <- rep("Cluster1", ncol(matrix_num))
  }

## Use scDesign2 for training countmatrix
cat("Model fitting...\n")
n_cell_new <- ncol(matrix_num)
cell_type_sel <- unique(colnames(matrix_num))
cell_type_prop <- table(colnames(matrix_num))[cell_type_sel]
copula_result <- fit_model_scDesign2_new(matrix_num, cell_type_sel, sim_method = 'copula',
                                       ncores = length(cell_type_sel))
cat("Generating synthetic count matrix...\n")
simu_matrix <- scDesign2::simulate_count_scDesign2(copula_result, n_cell_new, sim_method = 'copula',
                                               cell_type_prop = cell_type_prop)
 rownames(simu_matrix) <- count_matrix[,1]
cat(sprintf("Writing out synthetic count matrix %s to %s...\n", out_directory, out_directory))
write.table(simu_matrix, sprintf("%s/%s.scDesign2Simulated.%s", out_directory, samplename, sample_format), sep="\t", row.names = TRUE,col.names = TRUE)
simu_count_pergene_vec <- rowSums(ceiling(simu_matrix/2))
write.table(simu_count_pergene_vec, sprintf("%s/%s.scDesign2Simulated.nPairsRegionmargional.txt", out_directory, samplename), row.names = FALSE,col.names = FALSE)
cat("Done.\n")
}


# samplename <- '10X_RNA_chr1_3073253_4526737.countmatrix' 
# sample_format<-'txt'
# directory <- '/home/gayan/Projects/scATAC_Simulator/package_development/package_results/20220310_10X_scRNAseq_NONINPUT'
# out_directory <- '/home/gayan/Projects/scATAC_Simulator/package_development/package_results/20220310_10X_scRNAseq_NONINPUT'

scRNA_runSyntheticCount <- function(samplename, sample_format, directory, out_directory, cluster_prestep=TRUE){
  cat(sprintf("Reading count matrix %s.%s...\n", samplename, sample_format))
  count_matrix <- read.table(sprintf("%s/%s.%s", directory, samplename, sample_format), sep="\t",header = FALSE)
  matrix_num <- data.matrix(count_matrix[,2:ncol(count_matrix)])
count_pergene_vec <- rowSums(matrix_num)
  write.table(count_pergene_vec, sprintf("%s/%s.real.nReadRegionmargional.txt",out_directory, samplename), row.names = FALSE,col.names = FALSE)
  matrix_num_nonzero <- matrix_num[count_pergene_vec>0,]
  ## Clustering
  if (cluster_prestep == TRUE){
    cat("Louvain Clustering Before Simulation...\n")
    set.seed(2022)
    clustering_result <- get_cluster_seurat(matrix_num_nonzero, 'UMI')
    colnames(matrix_num) <- clustering_result$clustering_result
    write.table(clustering_result$clustering_result, sprintf("%s/%s.LouvainClusterResults.%s", out_directory, samplename, sample_format), sep="\n", row.names = FALSE,col.names = FALSE)
  } else {
    colnames(matrix_num) <- rep("Cluster1", ncol(matrix_num))
  }

## Use scDesign2 for training countmatrix
cat("Model fitting...\n")
n_cell_new <- ncol(matrix_num)
cell_type_sel <- unique(colnames(matrix_num))
cell_type_prop <- table(colnames(matrix_num))[cell_type_sel]
copula_result <- fit_model_scDesign2_new(matrix_num, cell_type_sel, sim_method = 'copula',
                                       ncores = length(cell_type_sel))
cat("Generating synthetic count matrix...\n")
simu_matrix <- scDesign2::simulate_count_scDesign2(copula_result, n_cell_new, sim_method = 'copula',
                                               cell_type_prop = cell_type_prop)
rownames(simu_matrix) <- count_matrix[,1]
cat(sprintf("Writing out synthetic count matrix %s to %s...\n", out_directory, out_directory))
write.table(simu_matrix, sprintf("%s/%s.scDesign2Simulated.%s", out_directory, samplename, sample_format), sep="\t", row.names = TRUE,col.names = TRUE)
simu_count_pergene_vec <- rowSums(simu_matrix)
write.table(simu_count_pergene_vec, sprintf("%s/%s.scDesign2Simulated.nReadRegionmargional.txt", out_directory, samplename), row.names = FALSE,col.names = FALSE)
cat("Done.\n")
}










