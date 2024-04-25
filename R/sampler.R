sampler.gausslinear <- function(var, basis, sample_size){
  n <- length(var)

  var.proj <- lm(var ~ basis + 1)$fitted.values

  radius <- sqrt(sum(var^2) - sum(var.proj^2))

  var.sample <- matrix(rnorm(n * sample_size), nrow = n)

  var.sample.proj <- lm(var.sample ~ basis + 1)$fitted.values

  var.sample <- var.sample - var.sample.proj
  var.sample <- scale(var.sample, center = FALSE, scale = sqrt(colSums(var.sample^2)) / radius)
  var.sample <- var.proj + var.sample

  subspace_samples <- is.na(colSums(var.sample))
  if(sum(subspace_samples) > sample_size/2) stop()
  for(bad_i in which(subspace_samples)){
    var.sample[, bad_i] <- var.sample[, sample(which(!subspace_samples), 1)]
  }

  return(var.sample)
}

#' Model-X conditional sampler
#'
#' Generates samples of some variables conditioning on the others in the model-X
#' setting.
#'
#' @param X n-by-p matrix of explanatory variables.
#' @param indices non-empty vector of indices of variables that are to be generated
#' samples. If \code{indices=1:p}, it generates unconditional samples of all
#' variables.
#' @param X.mean length p mean vector of the explanatory variables.
#' @param X.cov p-by-p covariance matrix of explanatory variables.
#' @param sample_size the number of samples to be generated.
#' @param storage an R environment that holds intermediate results to save
#' computational time when the function is calling repeatively in
#' \code{hFDR:::se.modelX}.
#'
#' @return a n-by-length(indices)-by-sample_size array that holds the samples of
#' variables indexed by indices.
#'
#' @export
sampler.modelX.gauss <- function(X, indices, X.mean, X.cov, sample_size, storage = NULL){
  if(is.environment(storage) && !is.null(storage$R.cond)){
    trans.mat <- storage$trans.mat
    R.cond <- storage$R.cond
  } else{
    if(length(indices) < NCOL(X)){
      trans.mat <- X.cov[indices, -indices, drop = F] %*% solve(X.cov[-indices, -indices, drop = F])
      cov.cond <- X.cov[indices, indices, drop = F] - trans.mat %*% X.cov[-indices, indices, drop = F]
      R.cond <- chol(cov.cond)
    } else{
      trans.mat <- X.cov[indices, -indices, drop = F]
      R.cond <- chol(X.cov)
    }

    if(is.environment(storage)){
      storage$trans.mat <- trans.mat
      storage$R.cond <- R.cond
    }
  }

  n <- NROW(X)
  mean.cond <- t(X.mean[indices] + trans.mat %*% (t(X[, -indices, drop = F]) - X.mean[-indices]))
  X.indices.samples <- replicate(sample_size, mean.cond + matrix(rnorm(n * length(indices)), nrow = n) %*% R.cond)
}

sampler.gaussgraph <- function(X, i, j, sample_size){
  Xi.samples <- sampler.gausslinear(X[, i], X[, -c(i, j), drop = F], sample_size)
  Xj.samples <- sampler.gausslinear(X[, j], X[, -c(i, j), drop = F], sample_size)

  return(list(Xi = Xi.samples, Xj = Xj.samples))
}
