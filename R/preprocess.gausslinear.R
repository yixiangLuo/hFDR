preprocess.gausslinear <- function(X, y, select, discard_tail_prob = 1e-4){
  n <- NROW(X)
  p <- NCOL(X)

  y <- scale(y, center = T, scale = F)
  X <- scale(X, center = T, scale = F)
  if(select == "lasso") X <- scale(X, center = F, scale = sqrt(colSums(X^2)/n))
  else if(select == "fs") X <- scale(X, center = F, scale = sqrt(colSums(X^2)))
  else stop("only lasso or fs supported.")

  QR <- qr(X)

  pivot_back <- sort.list(QR$pivot)

  Q_X <- qr.Q(QR, complete = F)
  R_X <- qr.R(QR, complete = F)[, pivot_back]

  # compute the a basic matrix needed by cknockoff
  # compute vj = unit(X_{j.-j}), X_{j.-j} = X_j orthogonal projected onto X_{-j}. time O(p^3)
  vj_mat <- sapply(1:p, function(j){
    # Q_X[, j] = X_j orthogonal projected onto X_{1:j-1}
    # X_{j.-j} = Q_X[, j] orthogonal projected onto S, S:=(X_{j+1:p} orthogonal projected onto X_{1:j-1})
    #          <=> find a vector in span(X_{j:p}) that is perpendicular to S
    # "coord" is the coordinate of such a vector under the basis Q_X[, j:p]
    coord <- forwardsolve(t(R_X[j:p,j:p]), c(1, rep(0, p-j)))
    vj <- Q_X[, j:p] %*% matrix(coord, nrow = p-j+1)
    vj <- vj / sqrt(sum(vj^2))
  })

  vjy_obs <- c(matrix(y, nrow=1) %*% vj_mat)
  RSS_X <- abs(sum(y^2) - sum((matrix(y, nrow=1) %*% Q_X)^2))
  RSS_Xnoj <- RSS_X + vjy_obs^2


  Xv <- colSums(X * vj_mat)
  X_perpv_y <- c(t(y) %*% X) - Xv * vjy_obs

  trans <- list(Xv = Xv, X_perpv_y = X_perpv_y)

  vy_bound <- pmax(abs(vjy_quantile(discard_tail_prob, RSS_Xnoj, df = n-p)), abs(vjy_obs))
  Xy_bound <- sapply(1:p, function(j){
    sort(c(vjy_to_Xjy(-vy_bound[j], trans, j),
           vjy_to_Xjy(vy_bound[j], trans, j)))
  })

  data.pack <- list(X = X, y = y, n = n, p = p, vj_mat = vj_mat, vjy_obs = vjy_obs,
                    RSS_X = RSS_X, RSS_Xnoj = RSS_Xnoj,
                    trans = trans, Xy_bound = Xy_bound)

  return(data.pack)
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
