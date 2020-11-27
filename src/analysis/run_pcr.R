#Run Principal Component regression:

args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

set.seed(1000)

## Show arguments
print(data)
print(outputfile)
print(out_time)
print(meta)


#libraries
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scater)
  library(foreach)
  library(tictoc)
  library(dplyr)
  library(purrr)
  library(cluster)
  library(magrittr)
  library(tibble)
  library(PCAtools)
  library(scran)
})

# Read in data
sce <- readRDS(file = data)
meta <- readRDS(file = meta)

# vars
batch <- meta[["batch"]]
batch
celltype <- meta[["celltype"]]
celltype

### ------------ Principle component regression -----------------------###

# if (is.null(assays(sce)[["logcounts"]])) {
#   ### ------------ Standardize normalized counts-----------###
#   clusters <- quickCluster(sce, use.ranks=FALSE)
#   table(clusters)
#   sce <- computeSumFactors(sce, min.mean=0.1, cluster=clusters)
#   sce <-  logNormCounts(sce)
#   #assay_nam <- "counts"
# }
# # }else{
# #   assay_nam <- "logcounts"
# # }

assay_nam <- "logcounts"

#reduce sce to highly variable genes only to speed up
dec <- modelGeneVar(sce)
dec <- dec[order(dec$bio, decreasing=TRUE),] 
hvg <- getTopHVGs(dec, n=1000)
sce <- sce[hvg,]

tic.clearlog()
for( batch_var in batch ){
  tic(batch_var)
  print(batch_var)
  batch_d <- data.frame(batch = as.numeric(as.factor(colData(sce)[, batch_var]))) %>% 
    set_rownames(colnames(sce))
  p <- pca(assays(sce)[[assay_nam]], metadata = batch_d, rank = 100)
  data <- p$rotated
  metadata <- p$metadata
  components = getComponents(p, seq_len(100))
  xvals <- data.matrix(data[,which(colnames(data) %in% components)])
  yvals <- metadata$batch
  yvals <- data.matrix(yvals)
  
  # create correlation table
  corvals <- cor(xvals, yvals, use = 'pairwise.complete.obs', method = 'pearson')
  r_square <- corvals^2
  
  #get variance components for each PC
  var_p <- p$variance[1:length(r_square)]
  
  # sum of product of pcs variance component and regression coefficients
  sum_pcr <- sum(var_p * r_square)
  
  #scale to the overall variance reflected in PC1-PC100
  pcr_scaled <- sum_pcr/sum(var_p)
  
  colData(sce)[,paste0("pcr.", batch_var)] <- 
    rep(pcr_scaled, ncol(sce))
  toc(log = TRUE, quiet = TRUE)
}

#Get timings
log.txt <- tic.log(format = TRUE)
log.lst <- tic.log(format = FALSE)
tic.clearlog()
timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))
mean_t <- mean(timings)
names(mean_t) <- meta[["dataset_name"]]
writeLines(unlist(log.txt))


### -------------- save sce object ----------------------###
saveRDS(sce, file = outputfile)
saveRDS(mean_t, file = out_time)


