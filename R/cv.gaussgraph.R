#' Cross-validation in Gaussian graphical model
#'
#' Cross-validation for a estimating the precision matrix in the Gaussian
#' graphical model.
#'
#' @param X n-by-p matrix of variables.
#' @param lambda regularity parameter sequence.
#' @param precision.est a function that takes arguments (X, lambda). It
#' estimates the precision matrix at regularity parameters lambda based on the
#' observation X. The return is the estimated precision matrix of size p-by-p.
#' @param select a function that takes arguments (X, lambda) and does
#' variable pair selection. The return must be a logical matrix of size p(p-1)/2
#' by length(lambda). Converting a variable pair (i,j) to the column index and
#' vice versa can be found in hFDR:::pair_to_index and hFDR:::index_to_pair.
#' @param type.measure measurement of the error of the estimated precision
#' matrix based on the test data. If equals "-loglikelihood", it uses minus
#' log-likelihood, o.w. if equals "mse", it gets predicted values of \code{X[, j]}
#' based on \code{X[, -j]} from the test dataset and on the estimated precision
#' matrix from the training set, and use the averaged mean-squared prediction
#' error over p such prediction problems as the error measurement.
#' @param nfold number of folds in cross-validation.
#'
#' @return An object of class \code{cv.obj} that has similar structure as
#' \code{cv.glmnet}.
#'  \item{lambda}{the values of lambda used in the fits.}
#'  \item{cvm}{the mean cross-validated error - a vector of length
#'  length(lambda).}
#'  \item{cvsd}{estimate of standard error of cvm.}
#'  \item{cvup}{upper curve = cvm+cvsd.}
#'  \item{cvlo}{lower curve = cvm-cvsd.}
#'  \item{nzero}{number of selected variables at each lambda. NA if augument
#'  \code{select} is NULL.}
#'  \item{name}{a text string indicating type of measure (for plotting purposes).}
#'  \item{lambda.min}{value of lambda that gives minimum cvm.}
#'  \item{lambda.1se}{largest value of lambda such that error is within 1
#'  standard error of the minimum.}
#'  \item{index}{a one column matrix with the indices of lambda.min and
#'  lambda.1se.}
#'
#' @examples
#' # Please see vignettes "GaussianGraph"
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

  nzero <- if(is.function(select)) colSums(select(X, lambda)) else NULL

  name <- ifelse(type.measure == "-loglikelihood", "-Log Likelihood", "Mean-Squared Error")

  cv.obj <- structure(list(call = match.call(),
                           lambda = lambda,
                           cvm = measure.mean,
                           cvsd = measure.std,
                           cvup = measure.mean+measure.std,
                           cvlo = measure.mean-measure.std,
                           nzero = nzero,
                           name = name,
                           lambda.min = lambda[ind.min],
                           lambda.1se = lambda[ind.1se],
                           index = index),
                         class = 'cv.obj')

  return(cv.obj)
}
