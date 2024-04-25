#' @import stats
NULL

#' FDR estimator
#'
#' Estimates the FDR of a selection procedure at regularity parameter lambda.
#'
#' @param X n-by-p matrix of explanatory variables.
#' @param y response vector of length n. It is unused in the Gaussian graphical
#' model.
#' @param model specifies the model assumption. "gausslinear" for the Gaussian
#' linear model, "modelX" for the model-X setttings, and "gaussgraph" for the
#' Gaussian graphical model.
#' @param modelX only used in the model-X settings to specify the model
#' assumptions. When it equals "auto.gaussian" (default), the explanatory
#' variables X are considered multivariate Gaussian with the mean and covariance
#' estiamted from the observed data. If users want to specify the distribution
#' of X, then argument \code{modelX} should be a list of functions that contains
#' a conditional sampling function for X. See details.
#' @param select Variable selection procedure. Either a string "lasso" or "fs"
#' in the Gaussian linear model or a function that takes arguments (X, y, lambda)
#' and does variable selection. The return of the function must be a logical
#' matrix indicating being selected of size p by length(lambda).
#' If \code{select} is "lasso" or "fs", \code{model} must be "gausslinear". They
#' correspond to variable selection by lasso (hFDR::select.lasso) or by forward
#' stepwise regression (hFDR::select.fs). Built-in fast algorithms is used to
#' compute the FDR estimation.
#' If \code{select} is a function, Monte-Carlo with \code{n_sample.hfdr} samples
#' is used to compute the FDR estimation.
#' @param lambda regularity parameter sequence.
#' @param nlambda downsample the \code{lambda} sequence to be of length
#' \code{nlambda}.
#' @param psi the psi function to reduce the bias of the estimated FDR. Either a
#' string "pval" or a function that computes the psi value and its normalizer
#' (optional). If \code{psi="pval"}, default psi functions based on p-values is
#' used. For user-specific psi function, please see hFDR::psi.guasslinear for an
#' example of the format.
#' @param se logical. Estimate the standard error of FDR estimator or not
#' (default).
#' @param n_sample.hfdr number of Monte-Carlo samples used to estimate FDR when
#' \code{select} is a user-specified function.
#' @param n_sample.se number of Bootstrap samples used to estimate the standard
#' error of the FDR estimator when \code{se=True}.
#' @param n_cores number of cores to be used in computing the FDR estimator
#' parallelly.
#'
#' @return An object of class \code{hFDR}.
#'  \item{lambda}{the values of lambda used.}
#'  \item{hFDR}{the estimated FDR - a vector of length length(lambda).}
#'  \item{hFDR.decompose}{the estimated FDR contribution from each variable (or
#'  variable pair in the Gaussian graphical model). It is a matrix of size
#'  p-by-length(lambda) (or size (p(p-1)/2)-by-length(lambda) in the Gaussian
#'  graphical model). Summing over the columns of \code{hFDR.decompose} gives
#'  \code{hFDR}.}
#'  \item{hFDR.se}{the estimated standard error of \code{hFDR} - a vector of
#'  length length(lambda). Equals NULL if \code{se=False}.}
#'
#' @details
#'  In the model-X settings, if users want to specify the distribution of X,
#'  then argument \code{modelX} should be a list of functions that include items
#'  of the following.
#'
#'  \code{sampler.modelX}: a function that takes \code{(X, indices, sample_size)}
#'  and generates \code{sample_size} samples of explanatory variables indexed by
#'  \code{indices} conditioning on the others variables in the data matrix
#'  \code{X}.
#'
#'  \code{pred_fit}: Optional. a function that takes arguments (X.new, X, y,
#'  lambda). It train a model with (X, y) at regularity parameters lambda and
#'  predict the response at X.new. The return must be a matrix of size
#'  NROW(X.new) by length(lambda). s.e.(hFDR) is estimated via
#'  bootstrap sampling from a sparse model. If \code{pred_fit} is provided, the
#'  sparse model is obtained by cross-validation using \code{pred_fit}.
#'  Otherwise the sparse model consists of the variables with \code{psi!=0}.
#'
#' @examples
#' p <- 100; n <- 300; k <- 15
#' X <- matrix(rnorm(n*p), n)
#' nonzero <- sample(p, k)
#' beta <- 1 * (1:p %in% nonzero)
#' y <- X %*% beta + rnorm(n)
#'
#' n_lambda <- 40
#' lambda_max <- max(abs(t(scale(X, TRUE, FALSE)) %*% scale(y, TRUE, FALSE))) / n
#' sigma <- sqrt(sum(lm(y ~ X)$residuals^2) / (n-p-1))
#' lambda_min <- sigma / sqrt(n) / 10
#' lambda <- lambda_max * (lambda_min/lambda_max)^((0:n_lambda)/n_lambda)
#'
#' # simple example for Lasso selection. See more examples in the vignettes.
#' glmnet.cv <- glmnet::cv.glmnet(X, y, alpha = 1, nfolds = 10,
#'                        intercept = TRUE, standardize = TRUE, family = "gaussian")
#'
#' hFDR.res <- hFDR(X, y, model = "gausslinear", select = "lasso",
#'                  lambda = lambda, psi = "pval")
#'
#' plot(hFDR.res, glmnet.cv)
#'
#' @export
hFDR <- function(X, y = NULL, model = c("gausslinear", "modelX", "gaussgraph"),
                 modelX = "auto.gaussian", select = c("lasso", "fs"), lambda,
                 nlambda = length(lambda), psi = "pval", se = F,
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
    } else { NULL }
  } else if(model == "modelX"){
    if(is.character(modelX) && modelX == "auto.gaussian"){
      if(NROW(X) <= NCOL(X)) stop('`X` must have rows more than columns to use `auto.gaussian`.')
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
      pred_fit <- NULL
    } else if(is.list(modelX)){
      if(!is.function(psi)) stop('if `modelX != "auto.gaussian"`, `psi` must be a self-defined function.')
      sampler.modelX <- modelX$sampler.modelX
      if(!is.function(sampler.modelX)){
        stop('if `modelX != "auto.gaussian"`, `modelX` must be a list with modelX$sampler.modelX(X, j, sample_size) being a function that generates `sample_size` samples of X_j conditional on X_{-j}.')
      }
      pred_fit <- modelX$pred_fit
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
      se.modelX(X, y, select, lambda, psi, sampler.modelX, pred_fit, n_sample.hfdr, n_sample.se, n_cores)
    } else { NULL }
  } else if(model == "gaussgraph"){
    if(!is.function(select)){
      stop('`select` must be a function, consider `select = select.glasso` for selection by graphical lasso.')
    }
    if(!is.function(psi) && psi == "pval") psi <- psi.gaussgraph
    hFDR.val <- hFDR.gaussgraph(X, select, lambda, psi, n_sample.hfdr, n_cores)
    hFDR.se <- if(se){
      se.gaussgraph(X, select, lambda, psi, n_sample.hfdr, n_sample.se, n_cores)
    } else { NULL }
  }

  hFDR.obj <- structure(list(call = match.call(),
                             lambda = lambda,
                             hFDR = hFDR.val$hFDR,
                             hFDR.decompose = hFDR.val$hFDR.decompose,
                             hFDR.se = hFDR.se),
                        class = 'hFDR')

  return(hFDR.obj)
}
