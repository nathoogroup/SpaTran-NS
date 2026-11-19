# Preparation, provenance, integrity checks and per-gene chunk execution.
# Only this workflow's run directory is written. Legacy results are untouched.

cc_load <- function(script_dir) {
  for(pkg in c("INLA","Matrix","SpatialExperiment","SummarizedExperiment","fmesher"))
    if(!requireNamespace(pkg,quietly=TRUE)) stop("Missing R package: ",pkg)
  source(file.path(script_dir,"models.R"),local=.GlobalEnv)
  source(file.path(script_dir,"..","spot_selection.R"),local=.GlobalEnv)
}

cc_integer <- function(x,name,minimum=1L) {
  z <- suppressWarnings(as.numeric(x))
  if(length(z)!=1 || !is.finite(z) || z!=floor(z) || z<minimum || z>.Machine$integer.max)
    stop(name," must be an integer >= ",minimum,call.=FALSE)
  as.integer(z)
}

cc_md5 <- function(path) {
  z <- unname(as.character(tools::md5sum(path)))
  if(anyNA(z)) stop("Cannot fingerprint file: ",paste(path,collapse=", "))
  z
}

cc_code_hashes <- function(script_dir) {
  files <- c("models.R","workflow.R","prepare.R","run_task.R","combine.R","../spot_selection.R")
  stats::setNames(cc_md5(file.path(script_dir,files)),files)
}

cc_versions <- function() {
  packages <- c("INLA","Matrix","fmesher","SpatialExperiment","SummarizedExperiment")
  c(R=as.character(getRversion()),stats::setNames(
    vapply(packages,function(p) as.character(utils::packageVersion(p)),character(1)),packages))
}

cc_atomic_rds <- function(object,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  tmp <- tempfile(".pending_",tmpdir=dirname(path))
  on.exit(unlink(tmp),add=TRUE)
  saveRDS(object,tmp,compress=FALSE)
  if(!file.rename(tmp,path)) stop("Atomic rename failed: ",path)
}

cc_atomic_csv <- function(object,path) {
  tmp <- tempfile(".pending_",tmpdir=dirname(path))
  on.exit(unlink(tmp),add=TRUE)
  utils::write.csv(object,tmp,row.names=FALSE,na="")
  if(!file.rename(tmp,path)) stop("Atomic rename failed: ",path)
}

cc_preprocess <- function(raw,coords) {
  if(any(!is.finite(raw) | raw<0)) stop("Non-finite or negative counts")
  z <- stats::qnorm((rank(log1p(raw),ties.method="average")-.5)/length(raw))
  x <- coords$x-mean(coords$x); y <- coords$y-mean(coords$y)
  X <- cbind(1,x,y,x*x,y*y,x*y)
  fit <- stats::lm.fit(X,z)
  e <- as.numeric(fit$residuals)
  if(fit$rank!=6 || any(!is.finite(e)) || stats::sd(e)<1e-12)
    stop("Degenerate quadratic detrending design or residual variance")
  e
}

cc_prepare <- function(data_dir,run_dir,genes_per_task,script_dir,dataset_list,
                       config=cc_config()) {
  genes_per_task <- cc_integer(genes_per_task,"genes_per_task")
  if(dir.exists(run_dir) && length(list.files(run_dir,all.files=TRUE,no..=TRUE)))
    stop("Run directory is not empty; choose a NEW directory: ",run_dir)
  datasets <- trimws(readLines(dataset_list,warn=FALSE))
  datasets <- datasets[nzchar(datasets) & !startsWith(datasets,"#")]
  if(!length(datasets) || anyDuplicated(datasets) ||
     any(!grepl("^[A-Za-z0-9_-]+$",datasets))) stop("Invalid dataset list")
  files <- file.path(data_dir,paste0(datasets,".rds"))
  if(any(!file.exists(files))) stop("Missing data: ",paste(files[!file.exists(files)],collapse=", "))
  dir.create(file.path(run_dir,"units"),recursive=TRUE,showWarnings=FALSE)
  dir.create(file.path(run_dir,"chunks"),showWarnings=FALSE)
  units <- list(); tasks <- list()
  code_hashes <- cc_code_hashes(script_dir); versions <- cc_versions()
  for(ds in datasets) {
    path <- file.path(data_dir,paste0(ds,".rds"))
    hash <- cc_md5(path)
    spe <- readRDS(path)
    n_supplied <- ncol(spe)
    spe <- select_in_tissue_spots(spe,ds)
    cd <- SummarizedExperiment::colData(spe)
    samples <- if("sample_id" %in% names(cd)) as.character(cd$sample_id) else rep(ds,ncol(spe))
    if(anyNA(samples) || any(!nzchar(samples))) stop("Missing sample_id in ",ds)
    for(sample in unique(samples)) {
      section <- spe[,samples==sample,drop=FALSE]
      coords <- as.data.frame(SpatialExperiment::spatialCoords(section)[,1:2,drop=FALSE])
      names(coords) <- c("x","y")
      if(nrow(coords)<7 || any(!is.finite(as.matrix(coords)))) stop("Invalid coordinates: ",ds,"/",sample)
      x <- coords$x-mean(coords$x); y <- coords$y-mean(coords$y)
      if(qr(cbind(1,x,y,x*x,y*y,x*y))$rank!=6) stop("Rank-deficient spatial design: ",ds,"/",sample)
      counts <- SummarizedExperiment::assay(section,"counts")
      if(!inherits(counts,"sparseMatrix")) counts <- Matrix::Matrix(as.matrix(counts),sparse=TRUE)
      counts <- methods::as(counts,"CsparseMatrix")
      if(!inherits(counts,"dMatrix")) counts <- methods::as(counts,"dMatrix")
      if(any(!is.finite(counts@x) | counts@x<0)) stop("Invalid count matrix")
      keep <- Matrix::rowSums(counts>0)>=ceiling(ncol(counts)*.005) & Matrix::rowSums(counts)>=3
      if(!any(keep)) stop("No eligible genes: ",ds,"/",sample)
      rd <- SummarizedExperiment::rowData(section)
      ids <- rownames(section)
      if(is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
        stop("Gene IDs must be nonempty and unique: ",ds)
      symbols <- if("symbol" %in% names(rd)) as.character(rd$symbol) else
        if("gene_name" %in% names(rd)) as.character(rd$gene_name) else ids
      symbols[is.na(symbols) | !nzchar(symbols)] <- ids[is.na(symbols) | !nzchar(symbols)]
      unit_id <- sprintf("unit_%04d",length(units)+1L)
      context <- cc_context(coords,config)
      genes <- data.frame(gene_index=seq_len(sum(keep)),source_gene_index=which(keep),
                          gene_id=ids[keep],gene_name=symbols[keep],stringsAsFactors=FALSE)
      unit <- list(unit_id=unit_id,dataset_name=ds,sample_id=sample,genes=genes,
        counts=counts[keep,,drop=FALSE],coords=coords,context=context,
        source_path=normalizePath(path),source_md5=hash,
        n_spots_supplied_dataset=n_supplied,n_spots_retained_dataset=ncol(spe),
        n_spots=ncol(section),n_genes_supplied=nrow(section))
      unit_file <- file.path("units",paste0(unit_id,".rds"))
      cc_atomic_rds(unit,file.path(run_dir,unit_file))
      units[[unit_id]] <- list(unit_id=unit_id,dataset_name=ds,sample_id=sample,
        path=unit_file,md5=cc_md5(file.path(run_dir,unit_file)),source_md5=hash,
        genes=genes,n_spots=unit$n_spots,prior_mu=context$prior_mu,
        n_spots_supplied_dataset=n_supplied,n_spots_retained_dataset=ncol(spe),
        tissue_annotation_available="in_tissue" %in% names(cd))
      for(start in seq.int(1,nrow(genes),by=genes_per_task)) {
        tasks[[length(tasks)+1L]] <- data.frame(task_id=length(tasks)+1L,unit_id=unit_id,
          first_gene=start,last_gene=min(start+genes_per_task-1L,nrow(genes)))
      }
      message(unit_id," ",ds," / ",sample,": ",ncol(section)," retained spots, ",sum(keep)," genes")
    }
    if(!identical(hash,cc_md5(path))) stop("Source dataset changed during preparation: ",ds)
  }
  if(!identical(code_hashes,cc_code_hashes(script_dir))) stop("Code changed during preparation")
  manifest <- list(schema=1L,config=config,implementation_md5=code_hashes,software=versions,
    dataset_list=datasets,units=units,tasks=do.call(rbind,tasks),genes_per_task=genes_per_task,
    prepared_at=format(Sys.time(),tz="UTC",usetz=TRUE))
  cc_atomic_rds(manifest,file.path(run_dir,"manifest.rds"))
  cc_atomic_csv(manifest$tasks,file.path(run_dir,"tasks.csv"))
  unit_table <- do.call(rbind,lapply(units,function(u) data.frame(
    unit_id=u$unit_id,dataset_name=u$dataset_name,sample_id=u$sample_id,
    n_spots=u$n_spots,n_genes=nrow(u$genes),
    n_spots_supplied_dataset=u$n_spots_supplied_dataset,
    n_spots_retained_dataset=u$n_spots_retained_dataset,
    tissue_annotation_available=u$tissue_annotation_available,
    log_tau_prior_mean=u$prior_mu[1],log_kappa_prior_mean=u$prior_mu[2],
    source_md5=u$source_md5,unit_md5=u$md5)))
  cc_atomic_csv(unit_table,file.path(run_dir,"units.csv"))
  # This completion marker is consumed by submit.sh only after manifest creation.
  writeLines(as.character(nrow(manifest$tasks)),file.path(run_dir,"task_count.txt"))
  message("Prepared ",length(units)," sections and ",nrow(manifest$tasks)," gene-chunk tasks.")
  invisible(manifest)
}

cc_manifest <- function(run_dir,script_dir) {
  path <- file.path(run_dir,"manifest.rds")
  x <- readRDS(path)
  if(!identical(x$schema,1L) || !identical(x$implementation_md5,cc_code_hashes(script_dir)))
    stop("Manifest/code mismatch; prepare a new run with this implementation")
  if(!identical(x$software,cc_versions())) stop("R/package versions differ from preparation")
  if(!is.data.frame(x$tasks) || !identical(x$tasks$task_id,seq_len(nrow(x$tasks))))
    stop("Invalid task manifest")
  attr(x,"manifest_md5") <- cc_md5(path)
  x
}

cc_expected <- function(manifest,task_id) {
  task_id <- cc_integer(task_id,"task_id")
  if(task_id>nrow(manifest$tasks)) stop("task_id exceeds manifest")
  task <- manifest$tasks[task_id,,drop=FALSE]
  unit <- manifest$units[[task$unit_id]]
  if(is.null(unit)) stop("Unknown unit in task manifest")
  list(task=task,unit=unit,genes=unit$genes[seq.int(task$first_gene,task$last_gene),,drop=FALSE])
}

cc_validate_chunk <- function(chunk,manifest,task_id) {
  expected <- cc_expected(manifest,task_id)
  if(!is.list(chunk) || !isTRUE(chunk$complete) ||
     !identical(chunk$manifest_md5,attr(manifest,"manifest_md5")) ||
     !identical(chunk$task_id,as.integer(task_id)) ||
     !identical(chunk$unit_md5,expected$unit$md5) ||
     !identical(chunk$software,manifest$software)) stop("Chunk provenance mismatch")
  rows <- chunk$rows
  if(!is.data.frame(rows) || nrow(rows)!=nrow(expected$genes)) stop("Truncated or invalid chunk")
  for(k in names(expected$genes))
    if(!identical(rows[[k]],expected$genes[[k]])) stop("Chunk gene identity/order mismatch: ",k)
  for(k in c("dataset_name","sample_id","unit_id"))
    if(!identical(unique(rows[[k]]),expected$unit[[k]])) stop("Chunk section identity mismatch")
  if(!identical(unique(rows$n_spots),expected$unit$n_spots)) stop("Chunk spot count mismatch")
  for(m in c("S","NS","A")) {
    fields <- paste0(c("log_ml","status","error","warnings","seconds"),"_",m)
    if(!all(fields %in% names(rows))) stop("Missing model result columns")
    ok <- rows[[paste0("status_",m)]]=="ok"
    if(anyNA(ok) || any(!rows[[paste0("status_",m)]] %in% c("ok","failed"))) stop("Invalid fit status")
    lml <- rows[[paste0("log_ml_",m)]]; err <- rows[[paste0("error_",m)]]
    if(any(!is.finite(lml[ok])) || any(!is.na(lml[!ok])) ||
       any(!is.na(err[ok])) || anyNA(err[!ok]) || any(!nzchar(err[!ok]))) stop("Inconsistent fit/error fields")
    elapsed <- rows[[paste0("seconds_",m)]]
    if(any(!is.finite(elapsed) | elapsed<0)) stop("Invalid fit timing")
    for(i in seq_len(if(m=="S") 2 else 4)) {
      z <- rows[[paste0("theta",i,"_",m)]]
      if(length(z)!=nrow(rows) || any(!is.finite(z[ok])) || any(!is.na(z[!ok]))) stop("Invalid theta output")
    }
  }
  for(m in c("NS","A")) {
    want <- rows[[paste0("log_ml_",m)]]-rows$log_ml_S
    got <- rows[[paste0("log_bf_",m,"_vs_S")]]
    if(!isTRUE(all.equal(got,want,tolerance=1e-12))) stop("Incorrect/missing Bayes factor")
    flag <- rows[[if(m=="NS") "favors_nonstationarity" else "favors_anisotropy"]]
    if(!identical(flag,want>0)) stop("Incorrect/missing model-preference flag")
  }
  r <- sqrt(rows$theta3_A^2+rows$theta4_A^2)
  angle <- rep(NA_real_,nrow(rows)); defined <- is.finite(r) & r>1e-8
  angle[defined] <- (.5*atan2(rows$theta4_A[defined],rows$theta3_A[defined])*180/pi) %% 180
  if(!isTRUE(all.equal(rows$axis_ratio_at_posterior_mean,exp(r),tolerance=1e-12)) ||
     !isTRUE(all.equal(rows$axis_angle_degrees_at_posterior_mean,angle,tolerance=1e-12)))
    stop("Incorrect/missing anisotropy summary")
  invisible(TRUE)
}

cc_run_task <- function(run_dir,task_id,script_dir,retry_failed=FALSE,fit_fun=cc_fit) {
  manifest <- cc_manifest(run_dir,script_dir)
  expected <- cc_expected(manifest,task_id); task <- expected$task
  path <- file.path(run_dir,"chunks",sprintf("task_%06d.rds",task$task_id))
  lock <- paste0(path,".lock")
  if(dir.exists(lock)) stop("Task lock exists: ",lock,"; check for a running task")
  if(file.exists(path)) {
    old <- readRDS(path); cc_validate_chunk(old,manifest,task$task_id)
    failed <- any(old$rows$status_S!="ok" | old$rows$status_NS!="ok" | old$rows$status_A!="ok")
    if(!retry_failed || !failed) { message("Validated existing chunk: ",path); return(invisible(old)) }
  }
  # Atomic directory lock prevents duplicate submissions from racing a write.
  if(!dir.create(lock,showWarnings=FALSE)) stop("Task lock exists: ",lock,"; check for a running task")
  on.exit(unlink(lock,recursive=TRUE),add=TRUE)
  unit_path <- file.path(run_dir,expected$unit$path)
  if(!identical(cc_md5(unit_path),expected$unit$md5)) stop("Prepared unit file changed")
  unit <- readRDS(unit_path)
  rows <- lapply(expected$genes$gene_index,function(i) {
    message("Task ",task$task_id," / ",unit$dataset_name," / ",unit$sample_id," / gene ",i)
    y <- tryCatch(cc_preprocess(as.numeric(unit$counts[i,]),unit$coords),error=identity)
    if(inherits(y,"error")) {
      problem <- paste0("Preprocessing: ",conditionMessage(y))
      result <- cc_compare_gene(numeric(),unit$context,fit_fun=function(...) stop(problem))
    } else result <- cc_compare_gene(y,unit$context,fit_fun)
    cbind(unit$genes[i,,drop=FALSE],data.frame(dataset_name=unit$dataset_name,
      sample_id=unit$sample_id,unit_id=unit$unit_id,n_spots=unit$n_spots),result)
  })
  chunk <- list(complete=TRUE,task_id=task$task_id,manifest_md5=attr(manifest,"manifest_md5"),
    unit_md5=expected$unit$md5,rows=do.call(rbind,rows),software=cc_versions(),
    completed_at=format(Sys.time(),tz="UTC",usetz=TRUE))
  cc_validate_chunk(chunk,manifest,task$task_id)
  if(!identical(manifest$implementation_md5,cc_code_hashes(script_dir))) stop("Code changed during fitting")
  cc_atomic_rds(chunk,path)
  invisible(chunk)
}

cc_combine <- function(run_dir,script_dir) {
  manifest <- cc_manifest(run_dir,script_dir)
  chunks <- lapply(manifest$tasks$task_id,function(i) {
    path <- file.path(run_dir,"chunks",sprintf("task_%06d.rds",i))
    if(!file.exists(path) || dir.exists(paste0(path,".lock"))) stop("Missing or running task ",i)
    x <- readRDS(path); cc_validate_chunk(x,manifest,i); x$rows
  })
  rows <- do.call(rbind,chunks); rownames(rows) <- NULL
  if(anyDuplicated(paste(rows$unit_id,rows$gene_index,sep="/"))) stop("Duplicate section/gene keys")
  summaries <- lapply(names(manifest$units),function(id) {
    d <- rows[rows$unit_id==id,,drop=FALSE]
    do.call(rbind,lapply(c("NS","A"),function(m) {
      bf <- d[[paste0("log_bf_",m,"_vs_S")]]; good <- is.finite(bf)
      data.frame(unit_id=id,dataset_name=d$dataset_name[1],sample_id=d$sample_id[1],
        comparison=paste0(m,"_vs_S"),genes_attempted=nrow(d),valid_comparisons=sum(good),
        unavailable_comparisons=sum(!good),baseline_failures=sum(d$status_S!="ok"),
        alternative_failures=sum(d[[paste0("status_",m)]]!="ok"),
        favor_alternative=sum(bf[good]>0),bf_above_3=sum(bf[good]>log(3)),
        bf_above_10=sum(bf[good]>log(10)),
        proportion_favor_alternative=if(any(good)) mean(bf[good]>0) else NA_real_)
    }))
  })
  summary <- do.call(rbind,summaries)
  cc_atomic_rds(list(rows=rows,summary=summary,manifest=manifest,
                    combined_at=format(Sys.time(),tz="UTC",usetz=TRUE)),
                file.path(run_dir,"comparisons.rds"))
  cc_atomic_csv(rows,file.path(run_dir,"comparisons.csv"))
  cc_atomic_csv(summary,file.path(run_dir,"comparison_summary.csv"))
  message("Combined ",nrow(rows)," gene-section rows. Failed fits remain explicit, not negative evidence.")
  invisible(rows)
}
