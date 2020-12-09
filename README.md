# Benchmark mixing metrics

## General

### Aim

This pipeline is build to evaluate metrics to detect batch effects and/or mixing bias in single cell transcriptome data in single cell RNAseq data sets. The following metrics are included by now:
+ [Cellspecific Mixing Score (cms)](https://bioconductor.org/packages/release/bioc/html/CellMixS.html)
+ [k-nearest neighbour batch effect test (kbet)](https://github.com/theislab/kBET)
+ [Seurat's mixingMetric (mm)](https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8)
+ [local inverse simpson index (lisi)](https://www.nature.com/articles/s41592-019-0619-0)
+ [Shannon's entrophy](https://www.biorxiv.org/content/10.1101/2020.05.22.111211v1.full)
+ [Graph connectivity (graph)](https://www.biorxiv.org/content/10.1101/2020.05.22.111161v2.full)
+ [Average silouette width (asw)](https://www.biorxiv.org/content/10.1101/2020.05.22.111161v2.full)
+ [Principal component regression (pcr)](https://www.biorxiv.org/content/10.1101/2020.05.22.111161v2.full)


### Structure

Metrics are evaluated by the following criteria. For each criteria different scenarios are constructed using batch effects in real data and/or simulated batch effects. 
These batch effects and their simulations have been characterized in detail before (see https://almutlue.github.io/batch_snakemake/).
+ Sensitivity
+ Scalability
+ Comparability
+ Flexibility
+ Practical meassures as Runtime ..


### Details evaluation criteria

#### Sensitivity 
+ Test each metrics performance on all batch datasets.
+ Simulation with realtive increase of the batch log fold changes

#### Scalabilty
+ Change simulated batch effect stepwise
+ Run metrics and stepwise randomly permuted batch label

#### Comparability
+ Use same lfc distributions for different datasets

#### Flexibility
+ Simulation unbalanced batch effect
+ Simulation with tuned cellspecificity


## Results
View results [here](https://almutlue.github.io/mixing_benchmark/index.html).

## Setup

### Preparations

To setup this pipeline follow these instructions (Step 1 -2 explain one possible way to setup and run snakemake):

1. Set up and activate an Anaconda enviroment with __Snakemake >= v.5.6.0__ (or sth. eqivalent)  
2. Make sure your path to __R__ is exported within snakemake  
  * e.g. adding `*export PATH="/your/prefered/R/bin:$PATH"*` in your `*~/.bashrc*`
3. Clone this repository 
  * Caution: If you don't want to get all analysis that came with this repo you need to clean the `docs` directory from all files except `_site.yaml`
4. . Install all required R packages using **renv**
5. Create `**log**` and `**out**` directories.
6. Run: `*snakemake dir_setup*` to set up the neccessary directory structure to make all rules work.
7. If you want to view or share your analysis as website, activate github pages within your corresponding repo and specify the `*/docs*` as source directory. 

### Run
To run the entire pipeline:
1. Copy your preprocessed `*SingleCellExperiment*` dataset into `*/src/datasets/*`
2. Generate a corresponding metadata file and save it at `*/src/meta_files/*`
3. Run **snakemake**
4. Push results to github and refresh it's web deployment.
