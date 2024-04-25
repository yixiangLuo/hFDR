#' Default psi function in Gaussian linear model that thresholds variables by
#' their two-sided t-test p-values.
#'
#' @param X n-by-p matrix of explanatory variables.
#' @param y response vector of length n.
#' @param variables a sequence of indices of variables whose psi is to be
#' calculated.
#' @param threshold psi = 1 if the p-value of the variable is above threshold,
#' otherwise psi = 0.
#'
#' @return a list.
#'  \item{psi}{psi values of the variables}
#'  \item{normalizer}{the conditional expectation of psi given Sj and under the
#'  null Hj. This item is optional. If normalizer is not specified (NULL), it
#'  will be calculated approximately when called by hFDR::hFDR using Monte-Carlo.}
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

#' Default psi function in model-X settings with multivariate Gaussian X that
#' thresholds variables by their conditional randomization test p-values. The
#' test statistic is the absolute marginal covariance between the response y
#' and each explanatory variable.
#'
#' @param X n-by-p matrix of explanatory variables.
#' @param y response vector of length n.
#' @param variables a sequence of indices of variables whose psi is to be
#' calculated.
#' @param X.mean length p mean vector of the explanatory variables.
#' @param X.cov p-by-p covariance matrix of explanatory variables.
#' @param threshold psi = 1 if the p-value of the variable is above threshold,
#' otherwise psi = 0.
#'
#' @return a list.
#'  \item{psi}{psi values of the variables}
#'  \item{normalizer}{the conditional expectation of psi given Sj and under the
#'  null Hj. This item is optional. If normalizer is not specified (NULL), it
#'  will be calculated approximately when called by hFDR::hFDR using Monte-Carlo.}
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

#' Default psi function in the Gaussian graphical model that
#' thresholds variables by their two-sided t-test p-values.
#'
#' @param X n-by-p matrix of the variables.
#' @param pair_inds a sequence of indices of pairs of variables whose psi is to
#' be calculated. Converting a variable pair (i,j) to the pair index and vice
#' versa can be found in hFDR:::pair_to_index and hFDR:::index_to_pair.
#' @param threshold psi = 1 if the p-value of the variable is above threshold,
#' otherwise psi = 0.
#'
#' @return a list.
#'  \item{psi}{psi values of the variables}
#'  \item{normalizer}{the conditional expectation of psi given Sj and under the
#'  null Hj. This item is optional. If normalizer is not specified (NULL), it
#'  will be calculated approximately when called by hFDR::hFDR using Monte-Carlo.}
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
