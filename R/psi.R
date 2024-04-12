psi.guasslinear <- function(X, y, variables, threshold = 0.1){
  tvals <- lm_to_t(X, y)
  pvals <- pvals_t(tvals$tvals, tvals$df, side = "two")[variables]
  psi <- (pvals > threshold)
  normalizer <- rep(1 - threshold, length(variables))

  return(list(psi = psi, normalizer = normalizer))
}

psi.modelX.gauss <- function(X, y, variables, X.mean, X.cov, threshold = 0.1){
  pvals <- sapply(variables, function(j){
    trans.mat <- X.cov[j, -j, drop = F] %*% solve(X.cov[-j, -j, drop = F])
    mean.cond <- X.mean[j] + c(trans.mat %*% t(X[, -j, drop = F] - X.mean[-j]))
    cov.cond <- c(X.cov[j, j, drop = F] - trans.mat %*% X.cov[-j, j, drop = F])

    Xjy_mean.cond <- sum(y * mean.cond)
    Xjy_std.cond <- sqrt(sum(y^2) * cov.cond)

    pval <- 2 * pnorm(-abs(sum(X[, j] * y) - Xjy_mean.cond), sd = Xjy_std.cond)
  })
  psi <- (pvals > threshold)
  normalizer <- rep(1 - threshold, length(variables))

  return(list(psi = psi, normalizer = normalizer))
}
