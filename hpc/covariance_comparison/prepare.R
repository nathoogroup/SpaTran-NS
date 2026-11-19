#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
if(length(args)<2 || length(args)>4) stop(
  "Usage: Rscript hpc/covariance_comparison/prepare.R data_dir NEW_run_dir [genes_per_task=100] [dataset_list]")
script_dir <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])))
source(file.path(script_dir,"workflow.R")); cc_load(script_dir)
config <- cc_config(
  slope_sd=as.numeric(Sys.getenv("CC_SLOPE_SD",as.character(sqrt(10)))),
  axis_ratio_upper=as.numeric(Sys.getenv("CC_AXIS_RATIO_UPPER","4")),
  axis_ratio_tail_probability=as.numeric(Sys.getenv("CC_AXIS_RATIO_TAIL_PROB","0.05")))
cc_prepare(args[1],args[2],if(length(args)>=3) args[3] else 100L,script_dir,
  if(length(args)>=4) args[4] else file.path(script_dir,"..","analysis_datasets.txt"),config)
