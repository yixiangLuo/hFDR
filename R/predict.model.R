#' Document will be ready soon
#'
#' @export
predict.model <- function(X.new, X, y, select, lambda, predict.mle){
  selected <- select(X, y, lambda)
  prediction <- sapply(1:length(lambda), function(lambda_i){
    model <- selected[, lambda_i]
    predict.mle(X.new[, model, drop = F], X[, model, drop = F], y)
  })
  return(prediction)
}

#' Document will be ready soon
#'
#' @export
predict.mle.gausslinear <- function(X.new, X, y){
  if(NCOL(X) > 0){
    lm.coefs <- stats::lm(y ~ X + 1)$coefficients
    coef <- lm.coefs[-1]
    a0 <- lm.coefs[1]

    if(any(is.na(coef))){
      coef[is.na(coef)] <- 0
      warning("singular X created by cross-validation.")
    }
    return(X.new %*% coef + a0)
  } else{
    return(rep(mean(y), NROW(X.new)))
  }
}
