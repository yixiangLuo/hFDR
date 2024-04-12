forward_stepwise <- function(X, y, lambda){
  X <- scale(X, center = T, scale = F)
  y <- scale(y, center = T, scale = F)
  X <- scale(X, center = F, scale = sqrt(colSums(X^2)))

  n <- NROW(X)
  p <- NCOL(X)
  n_sels <- lambda
  max_step <- max(n_sels)
  if(max_step > p) stop("Cannot select more than NCOL(X) variables")

  candidates <- 1:p
  sel_mat <- matrix(F, nrow = p, ncol = max_step)

  for(step in 1:max_step){
    sel <- which.max(abs(as.vector(matrix(y, nrow=1) %*% X[, candidates])))
    X_sel <- X[, candidates[sel]]
    sel_mat[candidates[sel], step:max_step] <- T

    candidates <- candidates[-sel]
    y <- y - sum(y * X_sel) * X_sel
    for(j in candidates){
      X[, j] <- X[, j] - sum(X[, j] * X_sel) * X_sel
      X[, j] <- X[, j] / sqrt(sum(X[, j]^2))
    }
  }
  sel_mat <- sel_mat[, n_sels]

  return(sel_mat)
}
