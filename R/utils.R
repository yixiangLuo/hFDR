pair_to_index <- function(i, j){
  (j-1)*(j-2)/2 + i
}

index_to_pair <- function(index, p){
  j <- sum(((2:p)-1)*((2:p)-2)/2 < index) + 1
  i <- index - (j-1)*(j-2)/2
  list(i = i, j = j)
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
