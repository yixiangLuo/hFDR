#' Model fitting and predition
#'
#' Fit a sparse MLE model and predict the response value
#'
#' @param X.new matrix of new values for X at which predictions are to be made.
#' @param X matrix of explanatory variables for fitting the sparse MLE model.
#' @param y response vector for fitting the sparse MLE model.
#' @param select function that performs variable selection.
#' @param lambda regularity parameter sequence.
#' @param pred_fit.mle function that predicts the response value by the MLE model.
#'
#' @return A vector of predicted values of length NROW(X.new)
#'
#' @export
pred_fit.model <- function(X.new, X, y, select, lambda, pred_fit.mle){
  selected <- select(X, y, lambda)
  prediction <- sapply(1:length(lambda), function(lambda_i){
    model <- selected[, lambda_i]
    pred_fit.mle(X.new[, model, drop = F], X[, model, drop = F], y)
  })
  return(prediction)
}

#' Model fitting and predition
#'
#' Fit a MLE model and predict the response value in Gaussian linear model (OLS)
#'
#' @param X.new matrix of new values for X at which predictions are to be made.
#' @param X matrix of explanatory variables for fitting the MLE model.
#' @param y response vector for fitting the MLE model.
#'
#' @return a vector of predicted values of length NROW(X.new)
#'
#' @export
pred_fit.mle.gausslinear <- function(X.new, X, y){
  if(NCOL(X) > 0){
    lm.coefs <- lm(y ~ X + 1)$coefficients
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

#' Model fitting and predition
#'
#' Fit a MLE model and predict the response value by logistic regression
#'
#' @param X.new matrix of new values for X at which predictions are to be made.
#' @param X matrix of explanatory variables for fitting the MLE model.
#' @param y binary response vector for fitting the MLE model.
#'
#' @return a vector of predicted values of length NROW(X.new)
#'
#' @export
pred_fit.mle.logistic <- function(X.new, X, y){
  if(NCOL(X) > 1){
    mle <- glmnet::glmnet(X, y, lambda = 0,
                          intercept = T, standardize = T,
                          family = "binomial")
    return(c(predict(mle, newx = X.new, s = 0, type = "response")))
  } else{
    return(rep(mean(y), NROW(X.new)))
  }
}
