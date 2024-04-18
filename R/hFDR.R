#' @import stats
NULL

#' The FDR estimator
#'
#' Document will be ready soon
#'
#' @export
hFDR <- function(X, y = NULL, model = c("gausslinear", "modelX", "gaussgraph"),
                 modelX = "auto.gaussian", select = c("lasso", "fs"), lambda,
                 nlambda = length(lambda), psi = "pval", se = T,
                 n_sample.hfdr = 20, n_sample.se = 10, n_cores = 1){
  tryCatch({
    model <- match.arg(model)
  }, error = function(msg){
    stop("Current version only support Gaussian linear models, model-X settings, and Gaussian graphical models.")
  })

  tryCatch({
    select <- match.arg(select)
  }, error = function(msg){
    if(!is.function(select)){
      stop('`select` must be a function that performs variable select or the built-in "lasso" or "fs".')
    }
  })
  if(!is.function(select) && select == "lasso" & !requireNamespace("glmnet", quietly=T)) {
    stop("R package 'glmnet' is required but is not installed.")
  }

  if(!is.function(psi) && psi != "pval"){
    stop('`psi` must be the built-in "pval" or a function that guess the noise variable (see paper).')
  }

  # if(is.function(select) & is.null(lambda)){
  #   stop('`lambda` must be specified if `select` is a user-defined function.')
  # }

  lambda <- lambda[round(seq(from = 1, to = length(lambda), length.out = nlambda))]

  n_sample.hfdr <- as.integer(n_sample.hfdr)
  n_sample.se <- as.integer(n_sample.se)
  n_cores <- as.integer(n_cores)
  if(n_sample.hfdr <= 0) stop("n_sample.hfdr must be a positive integer.")
  if(n_sample.se <= 0) stop("n_sample.se must be a positive integer.")
  if(n_cores <= 0) stop("n_cores must be a positive integer.")

  if(model == "gausslinear"){
    if(!is.function(psi) && psi == "pval") psi <- psi.guasslinear
    hFDR.val <- hFDR.gausslinear(X, y, select, lambda, psi, n_sample.hfdr, n_cores)
    hFDR.se <- if(se){
      se.gausslinear(X, y, select, lambda, psi, n_sample.hfdr, n_sample.se, n_cores)
    } else { NA }
  } else if(model == "modelX"){
    if(is.character(modelX) && modelX == "auto.gaussian"){
      # if(NROW(X) <= NCOL(X)) stop('`X` must have rows more than columns to use `auto.gaussian`.')
      X.mean <- colMeans(X)
      X.cov <- var(X)
      if(!is.function(psi) && psi == "pval") {
        psi <- function(X, y, variables){
          psi.modelX.gauss(X, y, variables, X.mean, X.cov, threshold = 0.1)
        }
      }
      sampler.modelX <- function(X, ind, sample_size, storage = NULL){
        sampler.modelX.gauss(X, ind, X.mean, X.cov, sample_size, storage)
      }
      predict <- NULL
    } else if(is.list(modelX)){
      if(!is.function(psi)) stop('if `modelX != "auto.gaussian"`, `psi` must be a self-defined function.')
      sampler.modelX <- modelX$sampler.modelX
      if(!is.function(sampler.modelX)){
        stop('if `modelX != "auto.gaussian"`, `modelX` must be a list with modelX$sampler.modelX(X, j, sample_size) being a function that generates `sample_size` samples of X_j conditional on X_{-j}.')
      }
      predict <- modelX$predict
    } else stop('if `modelX != "auto.gaussian"`, `modelX` must be a list with modelX$sampler.modelX(X, j, sample_size) being a function that generates `sample_size` samples of X_j conditional on X_{-j}.')

    if(!is.function(select)){
      if(select == "lasso"){
        select <- select.lasso
      } else if(select == "fs"){
        select <- select.fs
      }
    }
    hFDR.val <- hFDR.modelX(X, y, select, lambda, psi, sampler.modelX, n_sample.hfdr, n_cores)
    hFDR.se <- if(se){
      se.modelX(X, y, select, lambda, psi, sampler.modelX, predict, n_sample.hfdr, n_sample.se, n_cores)
    } else { NA }
  } else if(model == "gaussgraph"){
    if(!is.function(select)){
      stop('`select` must be a function, consider `select = select.glasso` for selection by graphical lasso.')
    }
    if(!is.function(psi) && psi == "pval") psi <- psi.gaussgraph
    hFDR.val <- hFDR.gaussgraph(X, select, lambda, psi, n_sample.hfdr, n_cores)
    hFDR.se <- if(se){
      se.gaussgraph(X, select, lambda, psi, n_sample.hfdr, n_sample.se, n_cores)
    } else { NA }
  }

  hFDR.obj <- structure(list(call = match.call(),
                             lambda = lambda,
                             hFDR = hFDR.val$hFDR,
                             hFDR.decompose = hFDR.val$hFDR.decompose,
                             hFDR.se = hFDR.se),
                        class = 'hFDR')

  return(hFDR.obj)
}
