cv.model <- function(X, y, select, lambda, predict, nfold = 10){
  n <- NROW(X)
  p <- NCOL(X)

  batch_size <- round(n / nfold)
  batches <- lapply(1:nfold, function(f_i){
    if(f_i < nfold) ((f_i-1) * batch_size + 1) : (f_i * batch_size)
    else ((f_i-1) * batch_size + 1) : n
  })
  cv.errors <- sapply(1:nfold, function(f_i){
    test_batch <- batches[[f_i]]
    X.tr <- X[-test_batch, ]
    y.tr <- y[-test_batch]

    prediction <- predict.model(X[test_batch, ], X.tr, y.tr, select, lambda, predict)

    test_error <- sapply(1:NCOL(prediction), function(tune_i){
      mean((y[test_batch] - prediction[, tune_i])^2, na.rm = T)
    })

    return(test_error)
  })
  MSE.mean <- rowMeans(cv.errors)
  MSE.std <- sapply(1:NROW(cv.errors), function(row){
    sqrt(var(cv.errors[row, ]) / nfold)
  })

  ind.min <- which.min(MSE.mean)
  ind.1se <- min(which(MSE.mean[1:ind.min] <= MSE.mean[ind.min] + MSE.std[ind.min]))
  index <- matrix(c(ind.min, ind.1se), 2, dimnames = list(c("min", "1se"), "Lambda"))

  selected <- select(X, y, lambda)

  cv.predict <- structure(list(call = match.call(),
                               lambda = lambda,
                               cvm = MSE.mean,
                               cvsd = MSE.std,
                               cvup = MSE.mean+MSE.std,
                               cvlo = MSE.std-MSE.std,
                               nzero = colSums(selected),
                               name = "Mean-Squared Error",
                               lambda.min = lambda[ind.min],
                               lambda.1se = lambda[ind.1se],
                               index = index),
                          class = 'cv.obj')

  return(cv.predict)
}
