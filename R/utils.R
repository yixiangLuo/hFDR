lm_to_t <- function(X, y, Sigma = NULL){
  n <- NROW(X)
  p <- NCOL(X)

  if(is.null(Sigma)){
    Sigma <- solve(t(X) %*% X)
  }
  Xy <- t(X) %*% y
  df <- n - p

  zvals <- Sigma %*% Xy
  sigmahat <- as.numeric(sqrt((sum(y^2) - t(Xy) %*% zvals) / df))
  tvals <- zvals / sqrt(diag(Sigma)) / sigmahat

  return(list(tvals = tvals, df = df))
}

# compute the p-values of t-statistics
pvals_t <- function(tvals, df, side = "two"){
  if (side == "right"){
    pvals <- pt(tvals, df = df, lower.tail = FALSE)
  } else if (side == "left"){
    pvals <- pt(tvals, df = df, lower.tail = TRUE)
  } else if (side == "two"){
    pvals <- 2 * pt(abs(tvals), df = df, lower.tail = FALSE)
  }

  return(pvals)
}

Xjy_to_vjy <- function(Xjy, trans, j){
  (Xjy - trans$X_perpv_y[j]) / trans$Xv[j]
}

vjy_to_Xjy <- function(vjy, trans, j){
  trans$Xv[j] * vjy + trans$X_perpv_y[j]
}

# the distribution function of vjy conditional on Sj
vjy_CDF <- function(vjy, res_norm2, df){
  return(pt(vjy * sqrt(df) / sqrt(pmax(res_norm2 - vjy^2, 0)), df = df))
}

dir_trunc <- function(x, direction){
  if(length(x) == 0) return(NULL)

  # x <- sapply(x, function(e){
  #   if(direction * e > 0){e} else{direction * Inf}
  # })
  for(i in 1:length(x)){
    if(direction * x[i] <= 0){
      x[i] <- direction * Inf
    }
  }
  return(x)
}
