#' Cross-validation for a generic predictive model.
#'
#' @param X n-by-p matrix of explanatory variables.
#' @param y response vector of length n.
#' @param lambda regularity parameter sequence.
#' @param pred_fit a function that takes arguments (X.new, X, y, lambda). It
#' train a model with (X, y) at regularity parameters lambda and predict the
#' response at X.new. The return must be a matrix of size NROW(X.new) by
#' length(lambda).
#' @param select a function that takes arguments (X, y, lambda) and does
#' variable selection. The return must be a logical matrix of size p by
#' length(lambda).
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
#' # Please see vignettes "GaussianLinear" and "ModelX"
#'
#' @export
cv.model <- function(X, y, lambda, pred_fit, select = NULL, nfold = 10){
  n <- NROW(X)
  p <- NCOL(X)

  batch_size <- round(n / nfold)
  batches <- lapply(1:nfold, function(f_i){
    if(f_i < nfold) ((f_i-1) * batch_size + 1) : (f_i * batch_size)
    else ((f_i-1) * batch_size + 1) : n
  })
  cv.measure <- sapply(1:nfold, function(f_i){
    test_batch <- batches[[f_i]]
    X.tr <- X[-test_batch, ]
    y.tr <- y[-test_batch]

    prediction <- pred_fit(X[test_batch, ], X.tr, y.tr, lambda)

    test_error <- sapply(1:NCOL(prediction), function(tune_i){
      mean((y[test_batch] - prediction[, tune_i])^2, na.rm = T)
    })

    return(test_error)
  })
  MSE.mean <- rowMeans(cv.measure)
  MSE.std <- sapply(1:NROW(cv.measure), function(row){
    sqrt(var(cv.measure[row, ]) / nfold)
  })

  ind.min <- which.min(MSE.mean)
  ind.1se <- min(which(MSE.mean[1:ind.min] <= MSE.mean[ind.min] + MSE.std[ind.min]))
  index <- matrix(c(ind.min, ind.1se), 2, dimnames = list(c("min", "1se"), "Lambda"))

  nzero <- if(is.function(select)) colSums(select(X, y, lambda)) else NULL

  cv.obj <- structure(list(call = match.call(),
                           lambda = lambda,
                           cvm = MSE.mean,
                           cvsd = MSE.std,
                           cvup = MSE.mean+MSE.std,
                           cvlo = MSE.mean-MSE.std,
                           nzero = nzero,
                           name = "Mean-Squared Error",
                           lambda.min = lambda[ind.min],
                           lambda.1se = lambda[ind.1se],
                           index = index),
                          class = 'cv.obj')

  return(cv.obj)
}
