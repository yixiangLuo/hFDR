#' Document will be ready soon
#'
#' @export
cv.gaussgraph <- function(X, lambda, precision.est, select = NULL,
                          type.measure = c("-loglikelihood", "mse"), nfold = 10){
  type.measure <- match.arg(type.measure)
  n <- NROW(X)
  p <- NCOL(X)

  batch_size <- round(n / nfold)
  batches <- lapply(1:nfold, function(f_i){
    if(f_i < nfold) ((f_i-1) * batch_size + 1) : (f_i * batch_size)
    else ((f_i-1) * batch_size + 1) : n
  })

  cv.measure <- sapply(1:nfold, function(f_i){
    measure <- sapply(lambda, function(rho){
      test_batch <- batches[[f_i]]
      X.tr <- X[-test_batch, ]

      hPrecision <- precision.est(X.tr, rho)

      if(type.measure == "-loglikelihood"){
        measure <- -(log(det(hPrecision)) - sum(diag(var(X[test_batch, ]) %*% hPrecision)))
      } else{
        hSigma <- base::solve(hPrecision)
        X.predict <- sapply(1:p, function(j){
          hSigma[j, -j] %*% solve(hSigma[-j, -j]) %*% t(X[test_batch, -j])
        })
        measure <- mean(colMeans(X[test_batch, ] - X.predict)^2)
      }
      return(measure)
    })
  })

  measure.mean <- rowMeans(cv.measure)
  measure.std <- sapply(1:NROW(cv.measure), function(row){
    sqrt(var(cv.measure[row, ]) / nfold)
  })

  ind.min <- which.min(measure.mean)
  ind.1se <- min(which(measure.mean[1:ind.min] <= measure.mean[ind.min] + measure.std[ind.min]))
  index <- matrix(c(ind.min, ind.1se), 2, dimnames = list(c("min", "1se"), "Lambda"))

  nzero <- if(is.function(select)) colSums(select(X, lambda)) else NA

  name <- ifelse(type.measure == "-loglikelihood", "-Log Likelihood", "Mean-Squared Error")

  cv.obj <- structure(list(call = match.call(),
                              lambda = lambda,
                              cvm = measure.mean,
                              cvsd = measure.std,
                              cvup = measure.mean+measure.std,
                              cvlo = measure.std-measure.std,
                              nzero = nzero,
                              name = name,
                              lambda.min = lambda[ind.min],
                              lambda.1se = lambda[ind.1se],
                              index = index),
                         class = 'cv.obj')

  return(cv.obj)
}
