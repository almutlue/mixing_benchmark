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
  library(kBET)
  library(foreach)
  library(tictoc)
  library(dplyr)
  library(purrr)
  library(magrittr)
})

# Read in data
sce <- readRDS(file = data)
meta <- readRDS(file = meta)

# vars
batch <- meta[["batch"]]
batch
celltype <- meta[["celltype"]]
celltype

### ------------ kBET-------------------------###

tic.clearlog()
for( batch_var in batch ){
  tic(batch_var)
  print(batch_var)
  #split by celltype
  cts <- as.factor(colData(sce)[,celltype])
  cts <- levels(droplevels(cts))
  kBet_list <- lapply(cts, function(ct){
    print(ct)
    print(batch_var)
    sce_sub <- sce[, colData(sce)[,celltype] %in% ct]
    names(assays(sce_sub))
    sub_data <- t(assays(sce_sub)[["counts"]])
    batch_vec <- colData(sce_sub)[, batch_var]
    batch.estimate <- kBET(sub_data, batch_vec, plot=FALSE)
    average <- ifelse(is.na(batch.estimate), NA, batch.estimate$summary$kBET.observed[1])[1]
  }) %>% setNames(cts)

  kbet_res <- unlist(kBet_list[as.character(colData(sce)[,celltype])])
  colData(sce)[,paste0("kbet.", batch_var)] <- kbet_res
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


