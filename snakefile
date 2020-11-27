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

time_data = glob_wildcards(config["src_time_dat"] + "{time_dat_nam}.rds").time_dat_nam
logger.info("=====")
logger.info(time_data)
logger.info("=====")

var = glob_wildcards(config["out_time_dat"] + "hca_{var}.rds").var
logger.info("=====")
logger.info(var)
logger.info("=====")

var_sub = [ele for ele in var if not ele.startswith("1")]
logger.info("=====")
logger.info(var_sub)
logger.info("=====")

rand_batch = ["pbmc_roche", "csf_patient_random"]
time_all = ["hca", "pbmc2_pat"]
example_dat = ["csf_patient"]

logger.info("=====")
logger.info(example_dat)
logger.info("=====")

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

# --- Subworkflows --- #

## -------------------------------------------------------------------------- ##
## Time and memory consumption
## -------------------------------------------------------------------------- ##

subworkflow time:
   workdir: config["ROOT"]
   snakefile:  config["src_time"] + "snakefile_time"


# --- Build --- #

## -------------------------------------------------------------------------- ##
## All
## -------------------------------------------------------------------------- ##
rule all:
	input:
	  expand(config["out_norm"] + "norm_{datasets}.rds", datasets = datasets),
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
	  time(expand(config["out_time"] + "{datasets}_{var}_cms_default_time.rds", datasets = time_all, var = var)),
	  time(expand(config["out_time"] + "{datasets}_{var}_cms_bmin_time.rds", datasets = time_all, var = var)),
	  time(expand(config["out_time"] + "{datasets}_{var}_cms_kmin_time.rds", datasets = time_all, var = var)),
	  time(expand(config["out_time"] + "{datasets}_{var}_lisi_time.rds", datasets = time_all, var = var)),
	  time(expand(config["out_time"] + "{datasets}_{var}_isi_time.rds", datasets = time_all, var = var)),
	  time(expand(config["out_time"] + "{datasets}_{var}_wisi_time.rds", datasets = time_all, var = var)),
	  time(expand(config["out_time"] + "{datasets}_{var}_mm_time.rds", datasets = time_all, var = var)),
	  time(expand(config["out_time"] + "{datasets}_{var}_entropy_time.rds", datasets = time_all, var = var)),
	  time(expand(config["out_time"] + "{datasets}_{var}_kbet_time.rds", datasets = time_all, var = var)),
	  time(expand(config["out_time"] + "{datasets}_{var}_graph_connectivity_time.rds", datasets = time_all, var = var)),
	  time(expand(config["out_time"] + "{datasets}_{var}_asw_time.rds", datasets = time_all, var = var_sub)),
	  time(expand(config["out_time"] + "{datasets}_{var}_pcr_time.rds", datasets = time_all, var = var)),
	  time(config["docs"] + "evaluate_time_data.html"),
	  config["docs"] + "vis_batches.html",
	  config["docs"] + "generate_time_data.html",
	  config["docs"] + "generate_plot_obj.html",
	  config["docs"] + "batch_characteristics_metrics.html",
	  config["docs"] + "metric_scaling_simulations.html",
	  config["docs"] + "metric_scaling_unbalanced.html",
	  config["docs"] + "correlation_metrics.html",
	  config["docs"] + "summary_results.html",
	  config["docs"] + "generate_benchmark_figures.html",
	  config["docs"] + "index.html",
	  config["docs"] + "index_summary.html",
	  config["docs"] + "index_benchmark_tasks.html",
	  config["docs"] + "index_metrics.html",
	  config["docs"] + "index_characterization.html",
	  config["docs"] + "metric_mm.html",
	  config["docs"] + "metric_cms.html",
	  config["docs"] + "metric_kbet.html",
	  config["docs"] + "metric_lisi.html",
	  config["docs"] + "metric_graph.html",
	  config["docs"] + "metric_asw.html",
	  config["docs"] + "metric_entropy.html"


# --- Workflow --- #

## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## Prepare datasets
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##

rule normalize:
    input: 
        data = config["src_data"] + "{datasets}.rds",
        param = config["src_meta"] + "{datasets}_meta.rds",
        script = config["src_norm"] + "normalization.R"
    output:
        out = config["out_norm"] + "norm_{datasets}.rds"
    log:
        config["log_norm"] + "normalization_{datasets}.Rout"
    shell:
        '''R CMD BATCH --no-restore --no-save "--args data='{input.data}' params='{input.param}' outputfile='{output.out}'" {input.script} {log}'''
        

rule vis_batch:
    input: 
        script = config["src_com"] + "vis_batches.Rmd",
        check = expand(config["out_norm"] + "norm_{real_datasets}.rds", real_datasets = real_datasets)
    params:
        sce_name = ",".join(real_datasets),
        meta = config["src_meta"],
        sce = config["out_norm"],
        out = config["out_fig"],
        outdir = config["docs"]
    output:
        out = config["docs"] + "vis_batches.html",
        out_fig = config["out_fig"] + "vis_all_batches.pdf"
    shell:
        '''R -e "rmarkdown::render(input ='{input.script}', output_file=basename('{output.out}'), output_dir='{params.outdir}', knit_root_dir=getwd(), params = list(sce ='{params.sce}', meta='{params.meta}', sce_name ='{params.sce_name}', out ='{params.out}'))"'''
        

        
        
rule time_data:
    input: 
        script = config["src_time"] + "generate_time_data.Rmd",
    params:
        data = ",".join(time_data),
        meta = config["src_meta"],
        data_path = config["src_time_dat"],
        out_path = config["out_time_dat"],
        outdir = config["docs"],
    output:
        out = config["docs"] + "generate_time_data.html",
        out_size = config["out_time_dat"] + "summary_size.rds",
    shell:
        '''R -e "rmarkdown::render(input ='{input.script}', output_file=basename('{output.out}'), output_dir='{params.outdir}', knit_root_dir=getwd(), params = list(data ='{params.data}', meta='{params.meta}', data_path ='{params.data_path}', out ='{params.out_path}'))"'''
        
        
        

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
        data = config["out_norm"] + "norm_{datasets}.rds",
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
        metrics = ",".join(metrics),
        out_figpath = config["out_fig_obj"] + "random/random_{rand_batch}"
    output:
        out = config["docs"] + "metric_scaling_random_{rand_batch}.html",
        out_cor = config["out_cor"] + "cor_random_{rand_batch}.rds",
        out_fig = config["out_fig_obj"] + "random/random_{rand_batch}_all.rds"
    shell:
        '''R -e "rmarkdown::render(input ='{input.script}', output_file=basename('{output.out}'), output_dir='{params.outdir}', knit_root_dir=getwd(), params = list(data='{input.data}', meta ='{input.meta}', metrics = '{params.metrics}', out_cor ='{output.out_cor}', fig_res ='{params.out_figpath}'))"'''


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
        last = "asw",
        out_figpath = config["out_fig_obj"] + "char/char"
    output:
        out = config["docs"] + "batch_characteristics_metrics.html",
        out_cor = config["out_cor"] + "cor_batch_characteristics.rds",
        out_fig = config["out_fig_obj"] + "char/char_summary_all.rds"
    shell:
        '''R -e "rmarkdown::render(input ='{input.script}', output_file=basename('{output.out}'), output_dir='{params.outdir}', knit_root_dir=getwd(), params = list(sce='{params.sce}', sce_name ='{params.sce_name}', metrics='{params.metrics}', meta ='{params.meta}', summary = '{params.summary}', de = '{params.de}', last ='{params.last}', out_cor ='{output.out_cor}', fig_res ='{params.out_figpath}'))"'''


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
        last = "asw",
        out_figpath = config["out_fig_obj"] + "scaling/scaling"
    output:
        out = config["docs"] + "metric_scaling_simulations.html",
        out_cor = config["out_cor"] + "cor_sim_scaling.rds",
        out_res = config["out_res"] + "res_sim_scaling.rds",
        out_fig = config["out_fig_obj"] + "scaling/scaling_cor.rds"
    shell:
        '''R -e "rmarkdown::render(input ='{input.script}', output_file=basename('{output.out}'), output_dir='{params.outdir}', knit_root_dir=getwd(), params = list(sce='{params.sce}', sce_name ='{params.sce_name}', metrics='{params.metrics}', meta ='{params.meta}', last ='{params.last}', out_cor ='{output.out_cor}', out_res ='{output.out_res}', fig_res ='{params.out_figpath}'))"'''


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
        last = "asw",
        out_figpath = config["out_fig_obj"] + "unbalanced/unbalanced"
    output:
        out = config["docs"] + "metric_scaling_unbalanced.html",
        out_cor = config["out_cor"] + "cor_unbalanced.rds",
        out_res = config["out_res"] + "res_unbalanced.rds",
        out_fig = config["out_fig_obj"] + "unbalanced/unbalanced_limits.rds"
    shell:
        '''R -e "rmarkdown::render(input ='{input.script}', output_file=basename('{output.out}'), output_dir='{params.outdir}', knit_root_dir=getwd(), params = list(sce='{params.sce}', sce_name ='{params.sce_name}', metrics='{params.metrics}', last ='{params.last}', out_cor ='{output.out_cor}', out_res ='{output.out_res}', fig_res ='{params.out_figpath}'))"'''


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
        metrics = ",".join(metrics),
        out_figpath = config["out_fig_obj"] + "correlation/cor"
    output:
        out = config["docs"] + "correlation_metrics.html",
        out_fig = config["out_fig_obj"] + "correlation/cor_all.rds"
    shell:
        '''R -e "rmarkdown::render(input ='{input.script}', output_file=basename('{output.out}'), output_dir='{params.outdir}', knit_root_dir=getwd(), params = list(cor_path='{params.cor_path}', chars ='{params.chars}', metrics='{params.metrics}', fig_res ='{params.out_figpath}'))"'''


## -------------------------------------------------------------------------- ##   
## -------------------------------------------------------------------------- ##
## Summarize results
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##

rule summarize_results:
    input: 
        check = config["out_cor"] + "cor_unbalanced.rds",
        check1 = config["out_cor"] + "cor_sim_scaling.rds",
        check2 = config["out_cor"] + "cor_batch_characteristics.rds",
        check3 = expand(config["out_cor"] + "cor_random_{rand_batch}.rds", rand_batch = rand_batch),
        script = config["src_com"] + "summary_results.Rmd",
        unbalanced = config["out_res"] + "res_unbalanced.rds",
        sim_res = config["out_res"] + "res_sim_scaling.rds",
        time_res = config["out_bench"] + "summary_time_mem.rds"
    params:
        outdir = config["docs"],
        cor_path = config["out_cor"],
        chars = ",".join(chars),
        metrics = ",".join(metrics),
        out_figpath = config["out_fig_obj"] + "summary/sum"
    output:
        out = config["docs"] + "summary_results.html",
        out_fig = config["out_fig_obj"] + "summary/sum_tab.rds"
    shell:
        '''R -e "rmarkdown::render(input ='{input.script}', output_file=basename('{output.out}'), output_dir='{params.outdir}', knit_root_dir=getwd(), params = list(cor_path='{params.cor_path}', chars ='{params.chars}', metrics='{params.metrics}', unbalanced='{input.unbalanced}', sim_res='{input.sim_res}', time_res='{input.time_res}', fig_res ='{params.out_figpath}'))"'''


## -------------------------------------------------------------------------- ##   
## -------------------------------------------------------------------------- ##
## Generate figures 
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##

rule plot_objects:
    input: 
        script = config["src_plot"] + "generate_plot_obj.Rmd",
    params:
        outdir = config["docs"],
        ex_data = ",".join(example_dat),
        out_path = config["out_fig_obj"] + "plot_obj/",
        data_path = config["src_data"],
        sce_name = ",".join(sim_datasets),
        un_name = ",".join(un_datasets),
        sce_type = config["src_fig"] + "type_pancreas_sce.rds"
    output:
        out = config["docs"] + "generate_plot_obj.html",
        out_sim_scal = expand(config["out_fig_obj"] + "plot_obj/sim_scaling_{example_dat}.rds", example_dat = example_dat)
    shell:
        '''R -e "rmarkdown::render(input ='{input.script}', output_file=basename('{output.out}'), output_dir='{params.outdir}', knit_root_dir=getwd(), params = list(data_path='{params.data_path}', out_path ='{params.out_path}', ex_data='{params.ex_data}', sce_name='{params.sce_name}', un_name='{params.un_name}', sce_type='{params.sce_type}'))"'''
        
        
rule plot_figures:
    input: 
        script = config["src_plot"] + "generate_benchmark_figures.Rmd",
        random_res_all = config["out_fig_obj"] + "random/random_csf_patient_random_all.rds",
        random_pbmc = config["out_fig_obj"] + "random/random_pbmc_roche_all.rds",
        random_cms = config["out_fig_obj"] + "random/random_pbmc_roche_cms.batch.rds",
        random_lisi = config["out_fig_obj"] + "random/random_pbmc_roche_lisi.rds",
        random_ent = config["out_fig_obj"] + "random/random_pbmc_roche_entropy.rds",
        random_mm = config["out_fig_obj"] + "random/random_pbmc_roche_mm.rds",
        random_pcr = config["out_fig_obj"] + "random/random_pbmc_roche_pcr.rds",
        random_kbet = config["out_fig_obj"] + "random/random_pbmc_roche_kbet.rds",
        sca_limit = config["out_fig_obj"] + "scaling/scaling_limits.rds",
        sca_cor = config["out_fig_obj"] + "scaling/scaling_cor.rds",
        sim_data = config["out_fig_obj"] + "plot_obj/sim_scaling_csf_patient.rds",
        sca_trade_off = config["out_fig_obj"] + "scaling/scaling_trade_off.rds",
        unb_moderate = config["out_fig_obj"] + "unbalanced/unbalanced_moderate_batch_effect.rds",
        unb_limit = config["out_fig_obj"] + "unbalanced/unbalanced_limits_lolli.rds",
        cor = config["out_fig_obj"] + "correlation/cor_all.rds",
        cor_legend = config["out_fig_obj"] + "correlation/cor_legend.rds",
        summary =  config["out_fig_obj"] + "summary/sum_tab.rds",
        cor_char =  config["out_fig_obj"] + "char/char_spear_cor.rds",
        tsne_type =  config["out_fig_obj"] + "plot_obj/batch_type_adjustment.rds",
        rel_cells =  config["out_fig_obj"] + "char/char_cor_rel_vars_n_cells.rds",
        rel_genes =  config["out_fig_obj"] + "char/char_cor_rel_vars_n_genes.rds",
        sim_scal_dist =  config["out_fig_obj"] + "scaling/scaling_pbmc2_media.rds",
        met_char = config["out_fig_obj"] + "char/char_task1.rds",
        time =  config["out_fig_obj"] + "time/rss_cpu_all.rds",
        mem =  config["out_fig_obj"] + "time/rss_cells.rds"
    params:
        outdir = config["docs"],
        out_path = config["out_fig"],
        h_size =  config["src_fig"] + "over_real_size_hm.rds",
        h_type =  config["src_fig"] + "over_real_batch_type.rds",
        h_ct_spec =  config["src_fig"] + "over_real_ct_spec.rds",
        lfc_dist =  config["src_fig"] + "over_real_pbmc2_pat_lfc_dist.rds",
        gi_legend =  config["src_fig"] + "over_real_general_info_short_lgd.rds",
        h_gi =  config["src_fig"] + "over_real_general_info_short.rds",
        tern =  config["src_fig"] + "ternary_ct.rds",
        tern_nopan =  config["src_fig"] + "ternary_ct_nopan.rds",
        den1 =  config["src_fig"] + "over_real_cellbench_lfc_dist.rds",
        den2 =  config["src_fig"] + "over_real_hca_lfc_dist.rds",
        den3 =  config["src_fig"] + "over_real_csf_media_lfc_dist.rds",
        den4 =  config["src_fig"] + "over_real_kang_lfc_dist.rds",
        den5 =  config["src_fig"] + "over_real_pbmc2_media_lfc_dist.rds",
        sim1 =  config["src_fig"] + "over_sim_csf_media_lfc_dist.rds",
        sim2 =  config["src_fig"] + "over_sim_kang_lfc_dist.rds",
        sim3 =  config["src_fig"] + "over_sim_pbmc2_media_lfc_dist.rds",
        random_tsne = config["src_fig"] + "tsne_csf_patient_random.rds",
        tern_sim =  config["src_fig"] + "ternary_ct_sim.rds",
    output:
        out = config["docs"] + "generate_benchmark_figures.html",
        out_fig = config["out_fig"] + "results_random.pdf"
    shell:
        '''R -e "rmarkdown::render(input ='{input.script}', output_file=basename('{output.out}'), output_dir='{params.outdir}', knit_root_dir=getwd(), params = list(out_path='{params.out_path}', random_tsne ='{params.random_tsne}', random_res_all='{input.random_res_all}', random_pbmc='{input.random_pbmc}', random_cms='{input.random_cms}', random_mm='{input.random_mm}', random_ent='{input.random_ent}', random_lisi='{input.random_lisi}', random_pcr='{input.random_pcr}', random_kbet='{input.random_kbet}', sca_limit='{input.sca_limit}', sca_cor='{input.sca_cor}', sim_data='{input.sim_data}', sca_trade_off='{input.sca_trade_off}', unb_moderate = '{input.unb_moderate}', unb_limit='{input.unb_limit}', cor='{input.cor}', cor_legend='{input.cor_legend}', summary='{input.summary}', cor_char='{input.cor_char}', tsne_type='{input.tsne_type}', h_size='{params.h_size}', h_type='{params.h_type}', h_ct_spec='{params.h_ct_spec}', lfc_dist='{params.lfc_dist}', rel_cells='{input.rel_cells}', rel_genes='{input.rel_genes}', sim_scal_dist='{input.sim_scal_dist}', gi_legend ='{params.gi_legend}', h_gi ='{params.h_gi}', tern ='{params.tern}', tern_nopan ='{params.tern_nopan}', den1 ='{params.den1}', den2 ='{params.den2}', met_char ='{input.met_char}', tern_sim ='{params.tern_sim}', den3 ='{params.den3}', den4 ='{params.den4}', den5 ='{params.den5}', sim1 ='{params.sim1}', sim2 ='{params.sim2}', sim3 ='{params.sim3}', mem ='{input.mem}', time ='{input.time}'))"'''


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


## -------------------------------------------------------------------------- ##   
## -------------------------------------------------------------------------- ##
## Create metric description files
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##

rule met_cms:
    input: 
        script = config["src_metric"] + "metric_cms.Rmd"
    params:
        outdir = config["docs"]
    output:
        out = config["docs"] + "metric_cms.html"
    shell:
        '''R -e "rmarkdown::render(input ='{input.script}', output_file=basename('{output.out}'), output_dir='{params.outdir}', knit_root_dir=getwd())"'''

rule met_mm:
    input: 
        script = config["src_metric"] + "metric_mm.Rmd"
    params:
        outdir = config["docs"]
    output:
        out = config["docs"] + "metric_mm.html"
    shell:
        '''R -e "rmarkdown::render(input ='{input.script}', output_file=basename('{output.out}'), output_dir='{params.outdir}', knit_root_dir=getwd())"'''
        
rule met_lisi:
    input: 
        script = config["src_metric"] + "metric_lisi.Rmd"
    params:
        outdir = config["docs"]
    output:
        out = config["docs"] + "metric_lisi.html"
    shell:
        '''R -e "rmarkdown::render(input ='{input.script}', output_file=basename('{output.out}'), output_dir='{params.outdir}', knit_root_dir=getwd())"'''
        
rule met_asw:
    input: 
        script = config["src_metric"] + "metric_asw.Rmd"
    params:
        outdir = config["docs"]
    output:
        out = config["docs"] + "metric_asw.html"
    shell:
        '''R -e "rmarkdown::render(input ='{input.script}', output_file=basename('{output.out}'), output_dir='{params.outdir}', knit_root_dir=getwd())"'''
        
        
rule met_kbet:
    input: 
        script = config["src_metric"] + "metric_kbet.Rmd"
    params:
        outdir = config["docs"]
    output:
        out = config["docs"] + "metric_kbet.html"
    shell:
        '''R -e "rmarkdown::render(input ='{input.script}', output_file=basename('{output.out}'), output_dir='{params.outdir}', knit_root_dir=getwd())"'''
        
rule met_entropy:
    input: 
        script = config["src_metric"] + "metric_entropy.Rmd"
    params:
        outdir = config["docs"]
    output:
        out = config["docs"] + "metric_entropy.html"
    shell:
        '''R -e "rmarkdown::render(input ='{input.script}', output_file=basename('{output.out}'), output_dir='{params.outdir}', knit_root_dir=getwd())"'''
        
rule met_graph:
    input: 
        script = config["src_metric"] + "metric_graph.Rmd"
    params:
        outdir = config["docs"]
    output:
        out = config["docs"] + "metric_graph.html"
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
        
        
        
        
        
