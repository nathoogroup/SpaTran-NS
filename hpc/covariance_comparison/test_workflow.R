#!/usr/bin/env Rscript
# Fast numerical/integrity checks; --fit additionally runs four actual INLA fits.
script_dir <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])))
source(file.path(script_dir,"workflow.R")); cc_load(script_dir)
expect_error <- function(expr,pattern=NULL) {
  error <- tryCatch({force(expr); NULL},error=identity)
  stopifnot(inherits(error,"error"))
  if(!is.null(pattern)) stopifnot(grepl(pattern,conditionMessage(error),fixed=TRUE))
  invisible(error)
}
near <- function(a,b,tolerance=1e-8) stopifnot(max(abs(a-b))<tolerance)
coords <- expand.grid(x=seq(0,1,length.out=8),y=seq(0,1,length.out=8))
context <- cc_context(coords)

message("Testing SPDE precision, anisotropy geometry and normalized priors")
near(cc_anisotropy_matrix(c(0,0)),diag(2))
for(v in list(c(1,0),c(0,1),c(-.5,.8),c(.2,-.3))) {
  H <- cc_anisotropy_matrix(v); r <- sqrt(sum(v^2))
  near(det(H),1)
  near(eigen(H,symmetric=TRUE)$values,exp(c(r,-r)))
  direction <- c(cos(.5*atan2(v[2],v[1])),sin(.5*atan2(v[2],v[1])))
  near(as.vector(H %*% direction),exp(r)*direction)
  # Independent FEM implementation: H = exp(-r) [I + w w'].
  # gamma=1 also avoids version differences in fmesher's gamma convention.
  w <- sqrt(expm1(2*r))*c(direction,0)
  independent <- fmesher::fm_fem(context$mesh,order=2,
    aniso=list(gamma=rep(1,context$mesh$n),
               v=matrix(w,context$mesh$n,3,byrow=TRUE)))
  z <- context$fem
  G <- H[1,1]*z$gxx+H[1,2]*z$gxy+H[2,2]*z$gyy
  near(G,exp(-r)*independent$g1aniso)
}
rg <- cc_rgeneric(context)
for(delta in list(c(0,0),c(.3,-.4),c(-.4,.2))) {
  theta <- context$prior_mu+delta
  native <- INLA::inla.spde2.precision(context$S,theta)
  custom <- INLA::inla.rgeneric.q(rg,"Q",c(theta,0,0))
  near(Matrix::triu(native),Matrix::triu(custom))
  near(native,INLA::inla.spde2.precision(context$NS,c(theta,0,0)))
}
stopifnot(length(rg$f$rgeneric$definition("log.norm.const"))==0)
theta <- c(context$prior_mu,.3,-.5)
near(rg$f$rgeneric$definition("log.prior",theta),sum(dnorm(theta,
  c(context$prior_mu,0,0),c(rep(sqrt(10),2),rep(context$config$anisotropy_sd,2)),log=TRUE)))
near(exp(-log(4)^2/(2*context$config$anisotropy_sd^2)),.05)

mock_fit <- function(context,y,model) {
  ntheta <- if(model %in% c("S","S_rgeneric")) 2L else 4L
  hyper <- data.frame(mean=rep(.1,ntheta),sd=rep(.2,ntheta),
                      row.names=paste0("Theta",seq_len(ntheta)," for spatial"))
  list(mlik=matrix(switch(model,S=-10,NS=-8,A=-11,S_rgeneric=-10),1,1),
       summary.hyperpar=hyper,mode=list(mode.status=0))
}
message("Testing independent failure handling and Bayes-factor direction")
ok <- cc_compare_gene(seq_len(nrow(coords)),context,mock_fit)
stopifnot(ok$log_bf_NS_vs_S==2,ok$log_bf_A_vs_S==-1,
          ok$favors_nonstationarity,!ok$favors_anisotropy)
for(bad in c("S","NS","A")) {
  calls <- character()
  fit <- function(context,y,model) {
    calls <<- c(calls,model)
    if(model==bad) stop("intentional failed fit")
    mock_fit(context,y,model)
  }
  z <- cc_compare_gene(seq_len(nrow(coords)),context,fit)
  stopifnot(identical(calls,c("S","NS","A")),z[[paste0("status_",bad)]]=="failed")
  if(bad=="S") stopifnot(is.na(z$log_bf_NS_vs_S),is.na(z$log_bf_A_vs_S))
  if(bad=="NS") stopifnot(is.na(z$log_bf_NS_vs_S),z$log_bf_A_vs_S==-1)
  if(bad=="A") stopifnot(z$log_bf_NS_vs_S==2,is.na(z$log_bf_A_vs_S))
}
bad_status <- function(context,y,model) {
  z <- mock_fit(context,y,model); z$mode$mode.status <- 1; z
}
bad_lml <- function(context,y,model) {
  z <- mock_fit(context,y,model); z$mlik[1,1] <- Inf; z
}
bad_warning <- function(context,y,model) {
  warning("Newton optimizer did not converge"); mock_fit(context,y,model)
}
for(fun in list(bad_status,bad_lml,bad_warning))
  stopifnot(cc_safe_fit(context,seq_len(nrow(coords)),"S",fun)$status=="failed")

message("Testing in-tissue selection, serial-section separation and chunk integrity")
test_dir <- tempfile("covariance_comparison_test_"); dir.create(test_dir)
data_dir <- file.path(test_dir,"data"); dir.create(data_dir)
run_dir <- file.path(test_dir,"run")
set.seed(982)
counts <- matrix(rpois(4*2*nrow(coords),lambda=4),4)
counts[4,] <- 0 # QC removes one gene; duplicate display symbols remain valid.
dimnames(counts) <- list(paste0("ENSG",1:4),paste0("spot",seq_len(ncol(counts))))
spe <- SpatialExperiment::SpatialExperiment(
  assays=list(counts=Matrix::Matrix(counts,sparse=TRUE)),
  rowData=S4Vectors::DataFrame(symbol=c("same","same","third","excluded")),
  colData=S4Vectors::DataFrame(sample_id=rep(c("section1","section2"),each=nrow(coords)),
    in_tissue=rep(c(FALSE,rep(TRUE,nrow(coords)-1)),2)),
  spatialCoords=as.matrix(rbind(coords,coords)))
saveRDS(spe,file.path(data_dir,"Synthetic.rds"))
dataset_list <- file.path(test_dir,"datasets.txt"); writeLines("Synthetic",dataset_list)
manifest <- cc_prepare(data_dir,run_dir,2L,script_dir,dataset_list)
stopifnot(length(manifest$units)==2,nrow(manifest$tasks)==4,
          all(vapply(manifest$units,function(u) u$n_spots,numeric(1))==63),
          all(vapply(manifest$units,function(u) nrow(u$genes),integer(1))==3L))
expect_error(cc_prepare(data_dir,run_dir,2L,script_dir,dataset_list),"not empty")
expect_error(cc_combine(run_dir,script_dir),"Missing or running task")
first <- cc_run_task(run_dir,1L,script_dir,fit_fun=mock_fit)
for(i in 2:4) cc_run_task(run_dir,i,script_dir,fit_fun=mock_fit)
before <- cc_md5(file.path(run_dir,"chunks","task_000001.rds"))
cc_run_task(run_dir,1L,script_dir,fit_fun=function(...) stop("should be skipped"))
stopifnot(identical(before,cc_md5(file.path(run_dir,"chunks","task_000001.rds"))))
manifest <- cc_manifest(run_dir,script_dir)
broken <- first; broken$rows <- broken$rows[1,,drop=FALSE]
expect_error(cc_validate_chunk(broken,manifest,1L),"Truncated")
broken <- first; broken$rows$gene_id <- rev(broken$rows$gene_id)
expect_error(cc_validate_chunk(broken,manifest,1L),"gene identity")
broken <- first; broken$manifest_md5 <- "incorrect"
expect_error(cc_validate_chunk(broken,manifest,1L),"provenance")
broken <- first; broken$rows$log_bf_A_vs_S <- 999
expect_error(cc_validate_chunk(broken,manifest,1L),"Bayes factor")
broken <- first; broken$rows$status_NS[1] <- "failed"
expect_error(cc_validate_chunk(broken,manifest,1L),"fit/error")
broken <- first; broken$rows$axis_ratio_at_posterior_mean <- 0
expect_error(cc_validate_chunk(broken,manifest,1L),"anisotropy summary")
broken <- first; broken$software["R"] <- "different"
expect_error(cc_validate_chunk(broken,manifest,1L),"provenance")
rows <- cc_combine(run_dir,script_dir)
stopifnot(nrow(rows)==6,all(rows$log_bf_NS_vs_S==2),all(rows$log_bf_A_vs_S==-1))
# Check partial failures survive combination, and a retry replaces only failed chunks.
replacement <- first
replacement$rows[,names(ok)] <- cc_compare_gene(seq_len(nrow(coords)),context,
  function(context,y,model) {if(model=="NS") stop("test NS failure"); mock_fit(context,y,model)})
cc_validate_chunk(replacement,manifest,1L)
cc_atomic_rds(replacement,file.path(run_dir,"chunks","task_000001.rds"))
rows <- cc_combine(run_dir,script_dir)
stopifnot(sum(is.na(rows$log_bf_NS_vs_S))==2,all(rows$log_bf_A_vs_S==-1))
summary <- read.csv(file.path(run_dir,"comparison_summary.csv"))
stopifnot(summary$valid_comparisons[summary$unit_id=="unit_0001" & summary$comparison=="NS_vs_S"]==1)
cc_run_task(run_dir,1L,script_dir,retry_failed=TRUE,fit_fun=mock_fit)
stopifnot(all(cc_combine(run_dir,script_dir)$status_NS=="ok"))
lock <- file.path(run_dir,"chunks","task_000001.rds.lock"); dir.create(lock)
expect_error(cc_run_task(run_dir,1L,script_dir,fit_fun=mock_fit),"Task lock exists")
expect_error(cc_combine(run_dir,script_dir),"Missing or running task")
unlink(lock,recursive=TRUE)
# Provenance is checked even for an already complete chunk.
original <- readRDS(file.path(run_dir,"manifest.rds"))
modified <- original; modified$implementation_md5[1] <- "modified"
saveRDS(modified,file.path(run_dir,"manifest.rds"))
expect_error(cc_run_task(run_dir,1L,script_dir,fit_fun=mock_fit),"Manifest/code mismatch")
saveRDS(original,file.path(run_dir,"manifest.rds"))
expect_error(cc_preprocess(rep(0,nrow(coords)),coords),"Degenerate")
bad_spe <- spe; SummarizedExperiment::colData(bad_spe)$in_tissue[1] <- NA
stopifnot(ncol(select_in_tissue_spots(bad_spe))==126)
SummarizedExperiment::colData(bad_spe)$in_tissue <- rep("unknown",ncol(bad_spe))
expect_error(select_in_tissue_spots(bad_spe),"unrecognized")
message("Workflow fixtures: ",test_dir)

if("--fit" %in% commandArgs(TRUE)) {
  message("Running actual INLA fits: S, NS, A and isotropic rgeneric cross-check")
  set.seed(17); y <- sin(7*coords$x)+rnorm(nrow(coords),sd=.4)
  actual <- cc_compare_gene(y,context)
  print(actual[,c("log_ml_S","log_ml_NS","log_ml_A","log_bf_NS_vs_S","log_bf_A_vs_S")])
  stopifnot(actual$status_S=="ok",actual$status_NS=="ok",actual$status_A=="ok")
  isotropic_rgeneric <- cc_safe_fit(context,y,"S_rgeneric")
  stopifnot(isotropic_rgeneric$status=="ok")
  difference <- abs(isotropic_rgeneric$log_ml-actual$log_ml_S)
  message("Identical-model log marginal likelihood difference: ",signif(difference,4))
  stopifnot(difference<.05)
}
message("All covariance-comparison tests passed.")
