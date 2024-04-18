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

sampler.modelX.gauss <- function(X, ind, X.mean, X.cov, sample_size, storage = NULL){
  if(is.environment(storage) && !is.null(storage$R.cond)){
    trans.mat <- storage$trans.mat
    R.cond <- storage$R.cond
  } else{
    if(length(ind) < NCOL(X)){
      trans.mat <- X.cov[ind, -ind, drop = F] %*% solve(X.cov[-ind, -ind, drop = F])
      cov.cond <- X.cov[ind, ind, drop = F] - trans.mat %*% X.cov[-ind, ind, drop = F]
      R.cond <- chol(cov.cond)
    } else{
      trans.mat <- X.cov[ind, -ind, drop = F]
      R.cond <- chol(X.cov)
    }

    if(is.environment(storage)){
      storage$trans.mat <- trans.mat
      storage$R.cond <- R.cond
    }
  }

  n <- NROW(X)
  mean.cond <- t(X.mean[ind] + trans.mat %*% (t(X[, -ind, drop = F]) - X.mean[-ind]))
  X.ind.samples <- replicate(sample_size, mean.cond + matrix(rnorm(n * length(ind)), nrow = n) %*% R.cond)
}

sampler.gaussgraph <- function(X, i, j, sample_size){
  Xi.samples <- sampler.gausslinear(X[, i], X[, -c(i, j), drop = F], sample_size)
  Xj.samples <- sampler.gausslinear(X[, j], X[, -c(i, j), drop = F], sample_size)

  return(list(Xi = Xi.samples, Xj = Xj.samples))
}
