# Two separate alternatives against one stationary isotropic Matérn baseline.
# Stationary anisotropy: (kappa^2 - div H grad)(tau w) = W, det(H) = 1.
# Gaussian normalizing constants are retained for marginal-likelihood ratios.
# See README.md for parameterization, priors and references.

cc_config <- function(slope_sd = sqrt(10), axis_ratio_upper = 4,
                      axis_ratio_tail_probability = .05) {
  stopifnot(length(slope_sd) == 1, is.finite(slope_sd), slope_sd > 0,
            length(axis_ratio_upper) == 1, is.finite(axis_ratio_upper), axis_ratio_upper > 1,
            length(axis_ratio_tail_probability) == 1,
            is.finite(axis_ratio_tail_probability),
            axis_ratio_tail_probability > 0, axis_ratio_tail_probability < 1)
  list(schema = 1L, slope_sd = slope_sd,
       axis_ratio_upper = axis_ratio_upper,
       axis_ratio_tail_probability = axis_ratio_tail_probability,
       anisotropy_sd = log(axis_ratio_upper)/sqrt(-2*log(axis_ratio_tail_probability)),
       shared_theta_precision = .1, intercept_precision = .001,
       noise_precision_gamma = c(shape = 1, rate = 5e-5),
       spot_policy = "in_tissue_if_annotated_otherwise_supplied",
       sample_policy = "separate_sample_ids",
       normalization = "log1p_average_rank_normal_scores_quadratic_detrending",
       integration_strategy = "ccd")
}

cc_anisotropy_matrix <- function(v) {
  stopifnot(length(v) == 2, all(is.finite(v)))
  r <- sqrt(sum(v^2))
  a <- if (r < 1e-8) 1 + r^2/6 else sinh(r)/r
  matrix(c(cosh(r)+a*v[1], a*v[2], a*v[2], cosh(r)-a*v[1]), 2, 2)
}

# Assemble the three directional stiffness matrices once on the fixed mesh.
# G(H) = H11 Gxx + H12 Gxy + H22 Gyy, with Gxy containing both cross terms.
cc_directional_fem <- function(mesh) {
  tri <- mesh$graph$tv
  xy <- mesh$loc[,1:2,drop=FALSE]
  p1 <- xy[tri[,1],,drop=FALSE]; p2 <- xy[tri[,2],,drop=FALSE]
  p3 <- xy[tri[,3],,drop=FALSE]
  twice_area <- (p2[,1]-p1[,1])*(p3[,2]-p1[,2]) -
                (p3[,1]-p1[,1])*(p2[,2]-p1[,2])
  stopifnot(all(is.finite(twice_area)), all(abs(twice_area) > 0))
  area <- abs(twice_area)/2
  gx <- cbind(p2[,2]-p3[,2],p3[,2]-p1[,2],p1[,2]-p2[,2])/twice_area
  gy <- cbind(p3[,1]-p2[,1],p1[,1]-p3[,1],p2[,1]-p1[,1])/twice_area
  ij <- expand.grid(i=1:3,j=1:3)
  assemble <- function(cross = FALSE, direction = gx) {
    vals <- lapply(seq_len(nrow(ij)), function(k) {
      i <- ij$i[k]; j <- ij$j[k]
      area * if (cross) (gx[,i]*gy[,j]+gy[,i]*gx[,j]) else direction[,i]*direction[,j]
    })
    Matrix::sparseMatrix(i=as.vector(tri[,ij$i]), j=as.vector(tri[,ij$j]),
      x=unlist(vals,use.names=FALSE), dims=c(mesh$n,mesh$n))
  }
  fem <- fmesher::fm_fem(mesh,order=2)
  ans <- list(c0=fem$c0,gxx=assemble(),gyy=assemble(direction=gy),gxy=assemble(TRUE))
  discrepancy <- max(abs(ans$gxx+ans$gyy-fem$g1))
  if (discrepancy > 1e-8*max(1,max(abs(fem$g1))))
    stop("Directional FEM does not reproduce isotropic stiffness")
  mass <- Matrix::diag(ans$c0)
  stopifnot(all(is.finite(mass)),all(mass>0))
  ans$cinv <- Matrix::Diagonal(x=1/mass)
  adjacency <- abs(ans$gxx)+abs(ans$gyy)+abs(ans$gxy)
  adjacency@x[] <- 1
  ans$graph <- Matrix::triu(adjacency %*% adjacency + adjacency + Matrix::Diagonal(mesh$n))
  ans
}

cc_context <- function(coords, config = cc_config(), mesh = NULL) {
  coords <- as.data.frame(coords[,c("x","y"),drop=FALSE])
  stopifnot(nrow(coords)>6, all(is.finite(as.matrix(coords))))
  ranges <- vapply(coords,function(z) diff(range(z)),numeric(1))
  stopifnot(all(ranges>0))
  if (is.null(mesh)) mesh <- fmesher::fm_mesh_2d_inla(
    loc=as.matrix(coords), max.edge=c(max(ranges)/5,max(ranges)/3), cutoff=min(ranges)/20)
  # Use the existing INLA nominal variance/range centering on this mesh;
  # copy and store the SAME two prior means for all three fitted models.
  base <- INLA::inla.spde2.matern(mesh,alpha=2,
    B.tau=matrix(c(0,1,0),1,3),B.kappa=matrix(c(0,0,1),1,3),
    theta.prior.prec=config$shared_theta_precision)
  prior_mu <- as.numeric(base$param.inla$theta.mu)
  stopifnot(length(prior_mu)==2,all(is.finite(prior_mu)))
  Btau <- matrix(c(0,1,0,0,0),mesh$n,5,byrow=TRUE)
  Bkappa <- cbind(0,0,1,as.vector(scale(mesh$loc[,1])),as.vector(scale(mesh$loc[,2])))
  ns <- INLA::inla.spde2.matern(mesh,alpha=2,B.tau=Btau,B.kappa=Bkappa,
    theta.prior.mean=c(prior_mu,0,0),
    theta.prior.prec=diag(c(rep(config$shared_theta_precision,2),rep(1/config$slope_sd^2,2))))
  list(mesh=mesh,A=INLA::inla.spde.make.A(mesh,loc=as.matrix(coords)),
       S=base,NS=ns,fem=cc_directional_fem(mesh),prior_mu=prior_mu,config=config)
}

# Self-contained rgeneric callback: every non-package object is explicitly
# passed through inla.rgeneric.define (also works in INLA's external R process).
cc_anisotropic_rgeneric <- function(
    cmd=c("graph","Q","mu","initial","log.norm.const","log.prior","quit"),
    theta=NULL) {
  initial <- function() prior_mu
  if (is.null(theta)) theta <- initial()
  graph <- function() graph_matrix
  Q <- function() {
    if (anisotropic) {
      v <- theta[3:4]; r <- sqrt(sum(v*v))
      a <- if (r<1e-8) 1+r*r/6 else sinh(r)/r
      G <- (cosh(r)+a*v[1])*gxx + (cosh(r)-a*v[1])*gyy + a*v[2]*gxy
    } else G <- gxx+gyy
    k2 <- exp(2*theta[2])
    precision <- exp(2*theta[1])*(k2*k2*c0 + 2*k2*G + G %*% cinv %*% G)
    if (any(!is.finite(precision@x))) stop("Non-finite anisotropic precision")
    Matrix::triu(precision)
  }
  mu <- function() numeric(0)
  # INLA computes -n/2 log(2*pi) + log|Q|/2. Do NOT return zero here.
  log.norm.const <- function() numeric(0)
  log.prior <- function() sum(stats::dnorm(theta,prior_mu,prior_sd,log=TRUE))
  quit <- function() invisible(NULL)
  do.call(match.arg(cmd),list())
}

cc_rgeneric <- function(context, anisotropic = TRUE) {
  mu <- context$prior_mu
  sd <- rep(1/sqrt(context$config$shared_theta_precision),2)
  if (anisotropic) { mu <- c(mu,0,0); sd <- c(sd,rep(context$config$anisotropy_sd,2)) }
  z <- context$fem
  INLA::inla.rgeneric.define(cc_anisotropic_rgeneric,
    c0=z$c0,cinv=z$cinv,gxx=z$gxx,gyy=z$gyy,gxy=z$gxy,graph_matrix=z$graph,
    prior_mu=mu,prior_sd=sd,anisotropic=anisotropic)
}

cc_fit <- function(context, y, model=c("S","NS","A","S_rgeneric"), verbose=FALSE) {
  model <- match.arg(model)
  stopifnot(length(y)==nrow(context$A),all(is.finite(y)))
  spatial_model <- switch(model,S=context$S,NS=context$NS,
                          A=cc_rgeneric(context),S_rgeneric=cc_rgeneric(context,FALSE))
  stack <- INLA::inla.stack(data=list(y=as.numeric(y)),A=list(context$A,1),
    effects=list(list(spatial=seq_len(context$mesh$n)),list(intercept=rep(1,length(y)))))
  # All models have the same Gaussian observation and proper intercept priors.
  INLA::inla(y ~ -1 + intercept + f(spatial,model=spatial_model),
    data=INLA::inla.stack.data(stack),family="gaussian",
    control.predictor=list(A=INLA::inla.stack.A(stack),compute=FALSE),
    control.fixed=list(mean=0,prec=context$config$intercept_precision),
    control.family=list(hyper=list(prec=list(prior="loggamma",
      param=unname(context$config$noise_precision_gamma)))),
    control.inla=list(strategy="gaussian",int.strategy=context$config$integration_strategy),
    control.compute=list(config=FALSE,dic=FALSE,cpo=FALSE,mlik=TRUE),
    num.threads="1:1",verbose=verbose)
}

cc_safe_fit <- function(context,y,model,fit_fun=cc_fit) {
  warnings <- character()
  started <- proc.time()[[3]]
  result <- tryCatch(withCallingHandlers({
    fit <- fit_fun(context,y,model)
    log_ml <- as.numeric(fit$mlik[1,1])
    hyper <- fit$summary.hyperpar
    if (length(log_ml)!=1 || !is.finite(log_ml) || is.null(hyper) ||
        any(!is.finite(as.matrix(hyper)))) stop("Non-finite or missing model outputs")
    status <- fit$mode$mode.status
    if (length(status) && any(!is.finite(status) | status!=0))
      stop("INLA optimizer mode.status = ",paste(status,collapse=","))
    if (any(grepl("(fail|not|didn.t).*converg|Newton.*fail",warnings,ignore.case=TRUE)))
      stop("INLA reported convergence warning: ",paste(warnings,collapse=" | "))
    means <- vapply(seq_len(if(model %in% c("S","S_rgeneric")) 2L else 4L),function(i) {
      ix <- grep(paste0("^Theta",i," for spatial$"),rownames(hyper),ignore.case=TRUE)
      if(length(ix)!=1) stop("Missing spatial hyperparameter ",i)
      hyper$mean[ix]
    },numeric(1))
    list(log_ml=log_ml,theta=means,error=NA_character_,status="ok")
  },warning=function(w) {
    warnings <<- c(warnings,conditionMessage(w)); invokeRestart("muffleWarning")
  }),error=function(e) list(log_ml=NA_real_,theta=rep(NA_real_,if(model %in% c("S","S_rgeneric")) 2 else 4),
                           error=conditionMessage(e),status="failed"))
  result$warnings <- if(length(warnings)) paste(unique(warnings),collapse=" | ") else NA_character_
  result$seconds <- proc.time()[[3]]-started
  result
}

cc_compare_gene <- function(y,context,fit_fun=cc_fit) {
  fits <- lapply(c("S","NS","A"),function(m) cc_safe_fit(context,y,m,fit_fun))
  names(fits) <- c("S","NS","A")
  row <- list()
  for(m in names(fits)) {
    z <- fits[[m]]
    for(field in c("log_ml","status","error","warnings","seconds"))
      row[[paste0(field,"_",m)]] <- z[[field]]
    for(i in seq_along(z$theta)) row[[paste0("theta",i,"_",m)]] <- z$theta[i]
  }
  # Separate evidence axes, not a sequential gate or a four-way causal partition.
  row$log_bf_NS_vs_S <- fits$NS$log_ml-fits$S$log_ml
  row$log_bf_A_vs_S <- fits$A$log_ml-fits$S$log_ml
  row$favors_nonstationarity <- if(is.finite(row$log_bf_NS_vs_S)) row$log_bf_NS_vs_S>0 else NA
  row$favors_anisotropy <- if(is.finite(row$log_bf_A_vs_S)) row$log_bf_A_vs_S>0 else NA
  v <- fits$A$theta[3:4]
  r <- sqrt(sum(v*v))
  row$axis_ratio_at_posterior_mean <- exp(r)
  row$axis_angle_degrees_at_posterior_mean <- if(is.finite(r) && r>1e-8)
    (.5*atan2(v[2],v[1])*180/pi) %% 180 else NA_real_
  as.data.frame(row,stringsAsFactors=FALSE)
}
