#Run inverse simpson index using non weighted settings:

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
  library(scran)
  library(CellMixS)
  library(foreach)
  library(tictoc)
})

# Read in data
sce <- readRDS(file = data)
meta <- readRDS(file = meta)

# vars
batch <- meta[["batch"]]
celltype <- meta[["celltype"]]


### ------------ Isi-------------------------###
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
for( batch_var in batch ){
  tic(batch_var)
  print(batch_var)
  sce <- evalIntegration(metric = "isi",
                         sce = sce,
                         group = batch_var,
                         k = 200,
                         assay_name = assay_nam,
                         weight = FALSE,
                         n_dim = n_dim,
                         res_name = paste0("isi_", batch_var))
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


