# Note: scDblFinder installation is not straightforward
# Independent Rscript and python script since loading R packages are different
# Check if packages are installed
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("scDblFinder", quietly = TRUE)) {
BiocManager::install("plger/scDblFinder")
}

# Analysis
get_doublet_id <- function(samplename, directory, out_directory, omic_choice= c("RNA", "ATAC")){
    ## evaluate choices
    omic_choice <- match.arg(omic_choice)
    set.seed(123)
    # Load foreground feature data
    cat(sprintf("[scReadSim] Reading count matrix %s.txt...\n", samplename))
    count_matrix <- read.table(sprintf("%s/%s.txt", directory, samplename), sep="\t",header = FALSE)
    matrix_num <- data.matrix(count_matrix[,2:ncol(count_matrix)])
    # Identify doublet
    cat(sprintf("[scReadSim] Detecting doublets...\n"))
    if (omic_choice == "RNA"){
      suppressWarnings(sce_doub <- scDblFinder::scDblFinder(matrix_num))
    } else if (omic_choice == "ATAC") {
      # sce <- scDblFinder(sce, artificialDoublets=1, aggregateFeatures=TRUE, nfeatures=25, processing="normFeatures")
      suppressWarnings(sce_doub <- scDblFinder::scDblFinder(matrix_num))
    }
    # Output doublet cell index in Rdata for following synthetic data generation (load index from Rdata)
    cat(sprintf("[scReadSim] Saving doublets detection results...\n"))
    doub_classification_label <- sce_doub$scDblFinder.class
    save(doub_classification_label, file=sprintf("%s/doublet_classification.Rdata", out_directory))
    cat("[scReadSim] Done.\n")
}
