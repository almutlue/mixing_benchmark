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
  library(igraph)
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

if( is.null(reducedDims(sce)[["PCA"]])){
  sce <- runPCA(sce, ncomponents = 50, exprs_values = assay_nam)
}

tic.clearlog()
for( batch_var in batch ){
  tic(batch_var)
  print(batch_var)
  knn <- buildKNNGraph(sce, use.dimred = "PCA", k = 5)
  comp <- components(knn)
  cellt_filtered <- names(table(colData(sce)[, celltype]))[table(colData(sce)[, celltype]) > 30]
  ct_g <- lapply(cellt_filtered, function(ct){
    ind <- which(colData(sce)[,celltype] %in% ct)
    # sub <- induced_subgraph(knn, ind, impl = "create_from_scratch")
    # sub_comp <- components(sub)
    lcc <- max(table(comp[[1]][ind]))
    lcc_all <- lcc/length(ind)
  }) %>% unlist() %>% sum()
 
  gcn <- ct_g/length(cellt_filtered)
  
  colData(sce)[,paste0("graph_connectivity.", batch_var)] <- 
    rep(gcn, ncol(sce))
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


