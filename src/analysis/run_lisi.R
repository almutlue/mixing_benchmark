#Run lisi using lisi package:

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
  library(lisi)
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
if( is.null(reducedDims(sce)[["PCA"]]) || ncol(reducedDims(sce)[["PCA"]]) < 10 ){
  if( !"logcounts" %in% names(assays(sce)) ){
    sce <- runPCA(sce, ncomponents = 10, exprs_values = "counts")
  }else{
    sce <- runPCA(sce, ncomponents = 10)
  }
}

tic.clearlog()
for( batch_var in batch ){
  tic(batch_var)
  print(batch_var)
  lisi_res <- lisi::compute_lisi(reducedDim(sce, "PCA")[,1:10], colData(sce), batch_var)
  colData(sce)[,paste0("lisi.", batch_var)] <- lisi_res
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


