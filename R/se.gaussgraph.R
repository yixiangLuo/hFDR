se.gaussgraph <- function(X, select, lambda, psi, n_sample.hfdr, n_sample.se, n_cores = 1){
  p <- NCOL(X)
  n <- NROW(X)

  precision.est <- function(X, lambda){
    precision.est.mle.sparse(X, select, lambda, precise = F)
  }

  cv.res <- cv.gaussgraph(X, lambda, precision.est, nfold = 10)

  hPrecision <- precision.est.mle.sparse(X, select, cv.res$lambda.min, precise = T)
  hSigma <- base::solve(hPrecision)
  R <- chol(hSigma)

  hFDR.samples <- matrix(NA, length(lambda), n_sample.se)

  for(boot_i in 1:n_sample.se){
    X.boot <- matrix(rnorm(n*p), n) %*% R

    hFDR.boot <- hFDR.gaussgraph(X.boot, select, lambda, psi, n_sample.hfdr, n_cores)$hFDR

    hFDR.samples[, boot_i] <- hFDR.boot
  }

  se <- sapply(1:NROW(hFDR.samples), function(lambda_i){
    sqrt(var(hFDR.samples[lambda_i, ]))
  })

  return(se)
}
