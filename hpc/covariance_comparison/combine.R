#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
if(length(args)!=1) stop("Usage: Rscript hpc/covariance_comparison/combine.R run_dir")
script_dir <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])))
source(file.path(script_dir,"workflow.R")); cc_load(script_dir)
cc_combine(args[1],script_dir)
