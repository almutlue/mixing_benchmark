# Main Workflow - Benchmark different mixing metrics to estimate batch effect strength 
#
# Contributors: @almutlue 

configfile: "config.yaml"

import glob
import os
import logging

datasets = glob_wildcards(config["src_data"] + "{dataset_name}.rds").dataset_name
logger.info("=====")
logger.info(datasets)
logger.info("=====")

metrics = glob_wildcards(config["src_analysis"] + "run_{metric_name}.R").metric_name
logger.info("=====")
logger.info(metrics)
logger.info("=====")


chars = glob_wildcards(config["out_cor"] + "cor_{chars_nam}.rds").chars_nam
logger.info("=====")
logger.info(chars)
logger.info("=====")


rand_batch = ["pbmc_roche", "csf_patient_random"]

real_data = [ele for ele in datasets if not ele.startswith("sim_")]
real_dataset = [ele for ele in real_data if not ele.endswith("random")]
real_datasets = [ele for ele in real_dataset if not ele.startswith("un_")]
logger.info("=====")
logger.info(real_datasets)
logger.info("=====")

sim_datasets = [ele for ele in datasets if ele.startswith("sim_")]
logger.info("=====")
logger.info(sim_datasets)
logger.info("=====")

un_datasets = [ele for ele in datasets if ele.startswith("un_")]
logger.info("=====")
logger.info(un_datasets)
logger.info("=====")


# --- Build --- #

## -------------------------------------------------------------------------- ##
## All
## -------------------------------------------------------------------------- ##
rule all:
	input:
	  expand(config["out_cms_def"] + "{datasets}_cms_default_sce.rds", datasets = datasets),
	  expand(config["out_cms_kmin"] + "{datasets}_cms_kmin_sce.rds", datasets = datasets),
	  expand(config["out_cms_bmin"] + "{datasets}_cms_bmin_sce.rds", datasets = datasets),
	  expand(config["out_mm"] + "{datasets}_mm_sce.rds", datasets = datasets),
	  expand(config["out_wisi"] + "{datasets}_wisi_sce.rds", datasets = datasets),
	  expand(config["out_isi"] + "{datasets}_isi_sce.rds", datasets = datasets),
	  expand(config["out_entr"] + "{datasets}_entropy_sce.rds", datasets = datasets),
	  expand(config["out_lisi"] + "{datasets}_lisi_sce.rds", datasets = datasets),
	  expand(config["out_kbet"] + "{datasets}_kbet_sce.rds", datasets = datasets),
	  expand(config["out_pcr"] + "{datasets}_pcr_sce.rds", datasets = datasets),
	  expand(config["out_gcon"] + "{datasets}_graph_connectivity_sce.rds", datasets = datasets),
	  expand(config["out_asw"] + "{datasets}_asw_sce.rds", datasets = datasets),
	  expand(config["docs"] + "metric_scaling_random_{rand_batch}.html", rand_batch = rand_batch),
	  config["docs"] + "batch_characteristics_metrics.html",
	  config["docs"] + "metric_scaling_simulations.html",
	  config["docs"] + "metric_scaling_unbalanced.html",
	  config["docs"] + "correlation_metrics.html",
	  config["docs"] + "index.html",
	  config["docs"] + "index_summary.html",
	  config["docs"] + "index_benchmark_tasks.html",
	  config["docs"] + "index_metrics.html",
	  config["docs"] + "index_characterization.html"


# --- Workflow --- #

## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## Calculate metric scores
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##



## -------------------------------------------------------------------------- ##
## Calculate cms scores
## -------------------------------------------------------------------------- ##

## ------------------ run cms default ---------------------------------- ##
rule cms_default:
    input: 
        script = config["src_analysis"] + "run_cms_default.R",
        data = config["src_data"] + "{datasets}.rds",
        meta = config["src_meta"] + "{datasets}_meta.rds"
    output:
        out = config["out_cms_def"] + "{datasets}_cms_default_sce.rds",
        out_time = config["out_time"] + "{datasets}_cms_default_time.rds"
    log:
        config["log_cms_def"] + "{datasets}_cms_default.Rout"
    shell:
        '''R CMD BATCH --no-restore --no-save "--args data='{input.data}' meta='{input.meta}' outputfile='{output.out}' out_time='{output.out_time}'" {input.script} {log}'''



## ------------------ run cms kmin ---------------------------------- ##
rule cms_kmin:
    input: 
        script = config["src_analysis"] + "run_cms_kmin.R",
        data = config["out_cms_def"] + "{datasets}_cms_default_sce.rds",
        meta = config["src_meta"] + "{datasets}_meta.rds"
    output:
        out = config["out_cms_kmin"] + "{datasets}_cms_kmin_sce.rds",
        out_time = config["out_time"] + "{datasets}_cms_kmin_time.rds"
    log:
        config["log_cms_kmin"] + "{datasets}_cms_kmin.Rout"
    shell:
        '''R CMD BATCH --no-restore --no-save "--args data='{input.data}' meta='{input.meta}' outputfile='{output.out}' out_time='{output.out_time}'" {input.script} {log}'''


    
## ------------------ run cms bmin ---------------------------------- ##
rule cms_bmin:
    input: 
        script = config["src_analysis"] + "run_cms_bmin.R",
        data = config["out_cms_kmin"] + "{datasets}_cms_kmin_sce.rds",
        meta = config["src_meta"] + "{datasets}_meta.rds"
    output:
        out = config["out_cms_bmin"] + "{datasets}_cms_bmin_sce.rds",
        out_time = config["out_time"] + "{datasets}_cms_bmin_time.rds"
    log:
        config["log_cms_bmin"] + "{datasets}_cms_bmin.Rout"
    shell:
        '''R CMD BATCH --no-restore --no-save "--args data='{input.data}' meta='{input.meta}' outputfile='{output.out}' out_time='{output.out_time}'" {input.script} {log}'''
        


## -------------------------------------------------------------------------- ##
## Calculate other scores from CellMixS package
## -------------------------------------------------------------------------- ##

## ------------------ run mm ---------------------------------- ##
rule mm:
    input: 
        script = config["src_analysis"] + "run_mm.R",
        data = config["out_cms_bmin"] + "{datasets}_cms_bmin_sce.rds",
        meta = config["src_meta"] + "{datasets}_meta.rds"
    output:
        out = config["out_mm"] + "{datasets}_mm_sce.rds",
        out_time = config["out_time"] + "{datasets}_mm_time.rds"
    log:
        config["log_mm"] + "{datasets}_mm.Rout"
    shell:
        '''R CMD BATCH --no-restore --no-save "--args data='{input.data}' meta='{input.meta}' outputfile='{output.out}' out_time='{output.out_time}'" {input.script} {log}'''
        
        
## ------------------ run wisi ---------------------------------- ##
rule wisi:
    input: 
        script = config["src_analysis"] + "run_wisi.R",
        data = config["out_mm"] + "{datasets}_mm_sce.rds",
        meta = config["src_meta"] + "{datasets}_meta.rds"
    output:
        out = config["out_wisi"] + "{datasets}_wisi_sce.rds",
        out_time = config["out_time"] + "{datasets}_wisi_time.rds"
    log:
        config["log_wisi"] + "{datasets}_wisi.Rout"
    shell:
        '''R CMD BATCH --no-restore --no-save "--args data='{input.data}' meta='{input.meta}' outputfile='{output.out}' out_time='{output.out_time}'" {input.script} {log}'''


## ------------------ run isi ---------------------------------- ##
rule isi:
    input: 
        script = config["src_analysis"] + "run_isi.R",
        data = config["out_wisi"] + "{datasets}_wisi_sce.rds",
        meta = config["src_meta"] + "{datasets}_meta.rds"
    output:
        out = config["out_isi"] + "{datasets}_isi_sce.rds",
        out_time = config["out_time"] + "{datasets}_isi_time.rds"
    log:
        config["log_isi"] + "{datasets}_isi.Rout"
    shell:
        '''R CMD BATCH --no-restore --no-save "--args data='{input.data}' meta='{input.meta}' outputfile='{output.out}' out_time='{output.out_time}'" {input.script} {log}'''
        
        
        
## ------------------ run entropy ---------------------------------- ##
rule entropy:
    input: 
        script = config["src_analysis"] + "run_entropy.R",
        data = config["out_isi"] + "{datasets}_isi_sce.rds",
        meta = config["src_meta"] + "{datasets}_meta.rds"
    output:
        out = config["out_entr"] + "{datasets}_entropy_sce.rds",
        out_time = config["out_time"] + "{datasets}_entropy_time.rds"
    log:
        config["log_entr"] + "{datasets}_entropy.Rout"
    shell:
        '''R CMD BATCH --no-restore --no-save "--args data='{input.data}' meta='{input.meta}' outputfile='{output.out}' out_time='{output.out_time}'" {input.script} {log}'''


## -------------------------------------------------------------------------- ##
## Calculate other scores NOT from CellMixS package
## -------------------------------------------------------------------------- ##

## ------------------ run lisi ---------------------------------- ##
rule lisi:
    input: 
        script = config["src_analysis"] + "run_lisi.R",
        data = config["out_entr"] + "{datasets}_entropy_sce.rds",
        meta = config["src_meta"] + "{datasets}_meta.rds"
    output:
        out = config["out_lisi"] + "{datasets}_lisi_sce.rds",
        out_time = config["out_time"] + "{datasets}_lisi_time.rds"
    log:
        config["log_lisi"] + "{datasets}_lisi.Rout"
    shell:
        '''R CMD BATCH --no-restore --no-save "--args data='{input.data}' meta='{input.meta}' outputfile='{output.out}' out_time='{output.out_time}'" {input.script} {log}'''
        

## ------------------ run kbet ---------------------------------- ##
rule kbet:
    input: 
        script = config["src_analysis"] + "run_kbet.R",
        data = config["out_lisi"] + "{datasets}_lisi_sce.rds",
        meta = config["src_meta"] + "{datasets}_meta.rds"
    output:
        out = config["out_kbet"] + "{datasets}_kbet_sce.rds",
        out_time = config["out_time"] + "{datasets}_kbet_time.rds"
    log:
        config["log_kbet"] + "{datasets}_kbet.Rout"
    shell:
        '''R CMD BATCH --no-restore --no-save "--args data='{input.data}' meta='{input.meta}' outputfile='{output.out}' out_time='{output.out_time}'" {input.script} {log}'''
        
        
## ------------------ run pcr ---------------------------------- ##        
rule pcr:
    input: 
        script = config["src_analysis"] + "run_pcr.R",
        data = config["out_kbet"] + "{datasets}_kbet_sce.rds",
        meta = config["src_meta"] + "{datasets}_meta.rds"
    output:
        out = config["out_pcr"] + "{datasets}_pcr_sce.rds",
        out_time = config["out_time"] + "{datasets}_pcr_time.rds"
    log:
        config["log_pcr"] + "{datasets}_pcr.Rout"
    shell:
        '''R CMD BATCH --no-restore --no-save "--args data='{input.data}' meta='{input.meta}' outputfile='{output.out}' out_time='{output.out_time}'" {input.script} {log}'''        
        

## ------------------ run graph_connectivity ---------------------------------- ##        
rule graph_connectivity:
    input: 
        script = config["src_analysis"] + "run_graph_connectivity.R",
        data = config["out_pcr"] + "{datasets}_pcr_sce.rds",
        meta = config["src_meta"] + "{datasets}_meta.rds"
    output:
        out = config["out_gcon"] + "{datasets}_graph_connectivity_sce.rds",
        out_time = config["out_gcon"] + "{datasets}_graph_connectivity_time.rds"
    log:
        config["log_gcon"] + "{datasets}_graph_connectivity.Rout"
    shell:
        '''R CMD BATCH --no-restore --no-save "--args data='{input.data}' meta='{input.meta}' outputfile='{output.out}' out_time='{output.out_time}'" {input.script} {log}'''

        
## ------------------ run asw ---------------------------------- ##
rule asw:
    input: 
        script = config["src_analysis"] + "run_asw.R",
        data = config["out_gcon"] + "{datasets}_graph_connectivity_sce.rds",
        meta = config["src_meta"] + "{datasets}_meta.rds"
    output:
        out = config["out_asw"] + "{datasets}_asw_sce.rds",
        out_time = config["out_time"] + "{datasets}_asw_time.rds"
    log:
        config["log_asw"] + "{datasets}_asw.Rout"
    shell:
        '''R CMD BATCH --no-restore --no-save "--args data='{input.data}' meta='{input.meta}' outputfile='{output.out}' out_time='{output.out_time}'" {input.script} {log}'''
        

## -------------------------------------------------------------------------- ##   
## -------------------------------------------------------------------------- ##
## Compare results
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##

## -------------------------------------------------------------------------- ##
## Scaling random mixing
## -------------------------------------------------------------------------- ##

rule scale_ran:
    input: 
        data = config["out_asw"] + "{rand_batch}_asw_sce.rds",
        meta = config["src_meta"] + "{rand_batch}_meta.rds",
        script = config["src_com"] + "metric_scaling_random.Rmd"
    params:
        outdir = config["docs"],
        metrics = ",".join(metrics)
    output:
        out = config["docs"] + "metric_scaling_random_{rand_batch}.html",
        out_cor = config["out_cor"] + "cor_random_{rand_batch}.rds"
    shell:
        '''R -e "rmarkdown::render(input ='{input.script}', output_file=basename('{output.out}'), output_dir='{params.outdir}', knit_root_dir=getwd(), params = list(data='{input.data}', meta ='{input.meta}', metrics = '{params.metrics}', out_cor ='{output.out_cor}'))"'''


## -------------------------------------------------------------------------- ##
## Comparison with batch characteristics
## -------------------------------------------------------------------------- ##

rule batch_chars:
    input: 
        check = expand(config["out_asw"] + "{datasets}_asw_sce.rds", datasets = real_datasets),
        script = config["src_com"] + "batch_characteristics_metrics.Rmd"
    params:
        outdir = config["docs"],
        sce = config["out_asw"],
        sce_name = ",".join(real_datasets),
        metrics = ",".join(metrics),
        meta = config["src_meta"],
        summary = config["src_summary"] + "summary_",
        de = config["src_de"] + "de_",
        last = "asw"
    output:
        out = config["docs"] + "batch_characteristics_metrics.html",
        out_cor = config["out_cor"] + "cor_batch_characteristics.rds"
    shell:
        '''R -e "rmarkdown::render(input ='{input.script}', output_file=basename('{output.out}'), output_dir='{params.outdir}', knit_root_dir=getwd(), params = list(sce='{params.sce}', sce_name ='{params.sce_name}', metrics='{params.metrics}', meta ='{params.meta}', summary = '{params.summary}', de = '{params.de}', last ='{params.last}', out_cor ='{output.out_cor}'))"'''


## -------------------------------------------------------------------------- ##
## Scaling simulations
## -------------------------------------------------------------------------- ##

rule sim_scaling:
    input: 
        check = expand(config["out_asw"] + "{datasets}_asw_sce.rds", datasets = sim_datasets),
        script = config["src_com"] + "metric_scaling_simulations.Rmd"
    params:
        outdir = config["docs"],
        sce = config["out_asw"],
        sce_name = ",".join(sim_datasets),
        metrics = ",".join(metrics),
        meta = config["src_meta"],
        last = "asw"
    output:
        out = config["docs"] + "metric_scaling_simulations.html",
        out_cor = config["out_cor"] + "cor_sim_scaling.rds",
        out_res = config["out_res"] + "res_sim_scaling.rds"
    shell:
        '''R -e "rmarkdown::render(input ='{input.script}', output_file=basename('{output.out}'), output_dir='{params.outdir}', knit_root_dir=getwd(), params = list(sce='{params.sce}', sce_name ='{params.sce_name}', metrics='{params.metrics}', meta ='{params.meta}', last ='{params.last}', out_cor ='{output.out_cor}', out_res ='{output.out_res}'))"'''


## -------------------------------------------------------------------------- ##
## Scaling unbalanced
## -------------------------------------------------------------------------- ##

rule sim_unbalanced:
    input: 
        check = expand(config["out_asw"] + "{datasets}_asw_sce.rds", datasets = un_datasets),
        script = config["src_com"] + "metric_scaling_unbalanced.Rmd"
    params:
        outdir = config["docs"],
        sce = config["out_asw"],
        sce_name = ",".join(un_datasets),
        metrics = ",".join(metrics),
        last = "asw"
    output:
        out = config["docs"] + "metric_scaling_unbalanced.html",
        out_cor = config["out_cor"] + "cor_unbalanced.rds",
        out_res = config["out_res"] + "res_unbalanced.rds"
    shell:
        '''R -e "rmarkdown::render(input ='{input.script}', output_file=basename('{output.out}'), output_dir='{params.outdir}', knit_root_dir=getwd(), params = list(sce='{params.sce}', sce_name ='{params.sce_name}', metrics='{params.metrics}', last ='{params.last}', out_cor ='{output.out_cor}', out_res ='{output.out_res}'))"'''


## -------------------------------------------------------------------------- ##
## Correlation metrics
## -------------------------------------------------------------------------- ##

rule correlation_metrics:
    input: 
        check = config["out_cor"] + "cor_unbalanced.rds",
        check1 = config["out_cor"] + "cor_sim_scaling.rds",
        check2 = config["out_cor"] + "cor_batch_characteristics.rds",
        check3 = expand(config["out_cor"] + "cor_random_{rand_batch}.rds", rand_batch = rand_batch),
        script = config["src_com"] + "correlation_metrics.Rmd"
    params:
        outdir = config["docs"],
        cor_path = config["out_cor"],
        chars = ",".join(chars),
        metrics = ",".join(metrics)
    output:
        out = config["docs"] + "correlation_metrics.html"
    shell:
        '''R -e "rmarkdown::render(input ='{input.script}', output_file=basename('{output.out}'), output_dir='{params.outdir}', knit_root_dir=getwd(), params = list(cor_path='{params.cor_path}', chars ='{params.chars}', metrics='{params.metrics}'))"'''


## -------------------------------------------------------------------------- ##   
## -------------------------------------------------------------------------- ##
## Summarize results
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##

rule summarize_ results:
    input: 
        check = config["out_cor"] + "cor_unbalanced.rds",
        check1 = config["out_cor"] + "cor_sim_scaling.rds",
        check2 = config["out_cor"] + "cor_batch_characteristics.rds",
        check3 = expand(config["out_cor"] + "cor_random_{rand_batch}.rds", rand_batch = rand_batch),
        script = config["src_com"] + "summary_results.Rmd"
    params:
        outdir = config["docs"],
        cor_path = config["out_cor"],
        chars = ",".join(chars),
        metrics = ",".join(metrics),
        unbalanced = config["out_res"] + "res_unbalanced.rds",
        sim_res = config["out_res"] + "res_sim_scaling.rds"
    output:
        out = config["docs"] + "summary_results.html"
    shell:
        '''R -e "rmarkdown::render(input ='{input.script}', output_file=basename('{output.out}'), output_dir='{params.outdir}', knit_root_dir=getwd(), params = list(cor_path='{params.cor_path}', chars ='{params.chars}', metrics='{params.metrics}', unbalanced='{params.unbalanced}', sim_res='{params.sim_res}'))"'''





## -------------------------------------------------------------------------- ##   
## -------------------------------------------------------------------------- ##
## Create index files
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##

rule index:
    input: 
        script = config["src_index"] + "index.Rmd"
    params:
        outdir = config["docs"]
    output:
        out = config["docs"] + "index.html"
    shell:
        '''R -e "rmarkdown::render(input ='{input.script}', output_file=basename('{output.out}'), output_dir='{params.outdir}', knit_root_dir=getwd())"'''

rule index_char:
    input: 
        script = config["src_index"] + "index_characterization.Rmd"
    params:
        outdir = config["docs"]
    output:
        out = config["docs"] + "index_characterization.html"
    shell:
        '''R -e "rmarkdown::render(input ='{input.script}', output_file=basename('{output.out}'), output_dir='{params.outdir}', knit_root_dir=getwd())"'''
        
rule index_bench:
    input: 
        script = config["src_index"] + "index_benchmark_tasks.Rmd"
    params:
        outdir = config["docs"]
    output:
        out = config["docs"] + "index_benchmark_tasks.html"
    shell:
        '''R -e "rmarkdown::render(input ='{input.script}', output_file=basename('{output.out}'), output_dir='{params.outdir}', knit_root_dir=getwd())"'''
        
rule index_summ:
    input: 
        script = config["src_index"] + "index_summary.Rmd"
    params:
        outdir = config["docs"]
    output:
        out = config["docs"] + "index_summary.html"
    shell:
        '''R -e "rmarkdown::render(input ='{input.script}', output_file=basename('{output.out}'), output_dir='{params.outdir}', knit_root_dir=getwd())"'''
        
rule index_met:
    input: 
        script = config["src_index"] + "index_metrics.Rmd"
    params:
        outdir = config["docs"]
    output:
        out = config["docs"] + "index_metrics.html"
    shell:
        '''R -e "rmarkdown::render(input ='{input.script}', output_file=basename('{output.out}'), output_dir='{params.outdir}', knit_root_dir=getwd())"'''



# --- Optional Rules  --- #

## -------------------------------------------------------------------------- ##
## Directory setup
## -------------------------------------------------------------------------- ##
rule dir_setup:
    input: 
        script = config["src_data_mgt"] + "dir_setup.R"
    params:
        dir_names = config["dir_names"].replace(" ", "")
    log:
        config["log_dir"] + "dir_setup.Rout"
    shell:
        '''R CMD BATCH --no-restore --no-save "--args dir_names='{params.dir_names}'" {input.script} {log}'''
        
        
        
        
        
