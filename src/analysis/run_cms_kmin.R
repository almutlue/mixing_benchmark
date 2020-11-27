#Run cellspecific mixing score using default settings:

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
  library(CellMixS)
  library(scater)
  library(scran)
  library(foreach)
  library(tictoc)
})

# Read in data
sce <- readRDS(file = data)
meta <- readRDS(file = meta)

# vars
batch <- meta[["batch"]]
celltype <- meta[["celltype"]]


### ------------ CellMixS-------------------------###
if (is.null(assays(sce)[["logcounts"]])) {
  assay_nam <- "counts"
}else{
  assay_nam <- "logcounts"
}

if( !is.null(reducedDims(sce)[["PCA"]]) && ncol(reducedDims(sce)[["PCA"]]) < 10 ){
  n_dim = ncol(reducedDims(sce)[["PCA"]])
}else{
  n_dim = 10
}

tic.clearlog()
foreach( batch_var = batch ) %dopar% {
  tic(batch_var)
  sce <- evalIntegration(metric = "cms",
                         sce = sce,
                         group = batch_var,
                         assay_name = assay_nam,
                         k = 200,
                         k_min = 80,
                         n_dim = n_dim,
                         res_name = paste0("kmin80_", batch_var))
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


