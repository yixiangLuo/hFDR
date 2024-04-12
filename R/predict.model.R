predict.model <- function(X.new, X, y, select, lambda, predict.mle){
  selected <- select(X, y, lambda)
  prediction <- sapply(1:length(lambda), function(lambda_i){
    model <- selected[, lambda_i]
    predict.mle(X.new[, model, drop = F], X[, model, drop = F], y)
  })
  return(prediction)
}
