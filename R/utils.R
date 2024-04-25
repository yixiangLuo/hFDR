# convert a pair of variable (i,j) with i<j to an 1-d index
pair_to_index <- function(i, j){
  (j-1)*(j-2)/2 + i
}

# convert an 1-d index back to a pair of variable (i,j) with i<j
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

vjy_quantile <- function(prob, res_norm2, df){
  t_quantile <- qt(prob, df)
  if(abs(t_quantile) < Inf){
    vjy <- sign(t_quantile) * sqrt(t_quantile^2 * res_norm2 / (t_quantile^2 + df))
  } else{
    vjy <- sign(t_quantile) * sqrt(res_norm2)
  }
  return(vjy)
}

dir_trunc <- function(x, direction){
  if(length(x) == 0) return(NULL)
  for(i in 1:length(x)){
    if(direction * x[i] <= 0){
      x[i] <- direction * Inf
    }
  }
  return(x)
}
