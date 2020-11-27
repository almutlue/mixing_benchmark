#Run kBet:

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
  library(scater)
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

### ------------ Average silhouette width -----------------------###

if (is.null(assays(sce)[["logcounts"]])) {
  assay_nam <- "counts"
}else{
  assay_nam <- "logcounts"
}

tic.clearlog()
for( batch_var in batch ){
  tic(batch_var)
  print(batch_var)
  #sce <- runPCA(sce, ncomponents = 10, exprs_values = assay_nam)
  dist_sce <- dist(reducedDim(sce, "PCA"))
  clust_sce <- as.numeric(as.factor(colData(sce)[, batch_var]))
  sil <- silhouette(clust_sce, dist = dist_sce)
  ct_mean <- sil[,c("cluster", "sil_width")] %>% as_tibble() %>% 
    group_by(cluster) %>% summarize(casw = 1 - abs(mean(sil_width)))
  ct_all <- data.frame(clust_sce) %>% set_colnames("cluster") %>%
    left_join(ct_mean)
  colData(sce)[,paste0("asw.", batch_var)] <- 
    rep(1 - abs(mean(sil[, "sil_width"])), ncol(sce))
  colData(sce)[,paste0("sw.", batch_var)] <- 1 - abs(sil[, "sil_width"])
  colData(sce)[,paste0("casw.", batch_var)] <- ct_all$casw
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


