#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
if(length(args)<1 || length(args)>2) stop("Usage: Rscript hpc/covariance_comparison/run_task.R run_dir [task_id]")
script_dir <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])))
source(file.path(script_dir,"workflow.R")); cc_load(script_dir)
task <- if(length(args)==2) cc_integer(args[2],"task_id") else
  cc_integer(Sys.getenv("SLURM_ARRAY_TASK_ID"),"SLURM_ARRAY_TASK_ID") +
  cc_integer(Sys.getenv("CC_TASK_OFFSET","0"),"CC_TASK_OFFSET",0)
cc_run_task(args[1],task,script_dir,retry_failed=identical(Sys.getenv("CC_RETRY_FAILED","0"),"1"))
