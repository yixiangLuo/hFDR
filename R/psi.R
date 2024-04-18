#' Document will be ready soon
#'
#' @export
psi.guasslinear <- function(X, y, variables, threshold = 0.1){
  X <- scale(X, center = T, scale = F)
  y <- scale(y, center = T, scale = F)

  Precision <- solve(t(X) %*% X)
  Xy <- t(X) %*% y
  df <- NROW(X) - NCOL(X) - 1

  zvals <- Precision %*% Xy
  sigmahat <- as.numeric(sqrt((sum(y^2) - t(Xy) %*% zvals) / df))
  tvals <- zvals / sqrt(diag(Precision)) / sigmahat

  pvals <- 2 * pt(abs(tvals), df = df, lower.tail = FALSE)
  pvals <- pvals[variables]

  psi <- (pvals > threshold)
  normalizer <- rep(1 - threshold, length(variables))

  return(list(psi = psi, normalizer = normalizer))
}

#' Document will be ready soon
#'
#' @export
psi.modelX.gauss <- function(X, y, variables, X.mean, X.cov, threshold = 0.1){
  y <- y - mean(y)
  pvals <- sapply(variables, function(j){
    trans.mat <- X.cov[j, -j, drop = F] %*% solve(X.cov[-j, -j, drop = F])
    mean.cond <- X.mean[j] + c(trans.mat %*% (t(X[, -j, drop = F]) - X.mean[-j]))
    cov.cond <- c(X.cov[j, j, drop = F] - trans.mat %*% X.cov[-j, j, drop = F])

    Xjy_std.cond <- sqrt(sum(y^2) * cov.cond)

    pval <- 2 * pnorm(-abs(sum((X[, j]-mean.cond) * y)), sd = Xjy_std.cond)
  })
  psi <- (pvals > threshold)
  normalizer <- rep(1 - threshold, length(variables))

  return(list(psi = psi, normalizer = normalizer))
}

#' Document will be ready soon
#'
#' @export
psi.gaussgraph <- function(X, pair_inds, threshold = 0.1){
  p <- NCOL(X)
  psi <- matrix(NA, p, p)
  normalizer <- matrix(NA, p, p)
  for(j in 2:p){
    if(any(pair_to_index(1:(j-1), j) %in% pair_inds)){
      res <- psi.guasslinear(X[, -j, drop = F], X[, j], 1:(j-1), threshold)
      psi[1:(j-1), j] <- res$psi
      normalizer[1:(j-1), j] <- res$normalizer
    }
  }
  psi <- psi[upper.tri(psi)]
  psi <- psi[pair_inds]
  normalizer <- normalizer[upper.tri(normalizer)]
  normalizer <- normalizer[pair_inds]

  return(list(psi = psi, normalizer = normalizer))
}
