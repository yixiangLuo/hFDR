#' Document will be ready soon
#'
#' @export
precision.est.glasso <- function(X, lambda){
  glasso.res <- glasso::glasso(var(X), lambda)
  hPrecision <- glasso.res$wi
  return(hPrecision)
}

precision.est.mle.sparse <- function(X, select, lambda, precise = F){
  p <- NCOL(X)
  n <- NROW(X)

  model <- select(X, lambda)

  if(precise){
    # https://cvxr.rbind.io/cvxr_examples/cvxr_sparse_inverse_covariance_estimation/
    S <- var(X)
    hPrecision <- CVXR::Variable(p, p, PSD = TRUE)
    obj <- CVXR::Maximize(CVXR::log_det(hPrecision) - CVXR::matrix_trace(S %*% hPrecision))
    constraints <- lapply(which(!model), function(pair_i){
      pair <- index_to_pair(pair_i, p)
      i <- pair$i
      j <- pair$j
      hPrecision[i, j] == 0
    })
    problem <- CVXR::Problem(obj, constraints)
    result <- CVXR::psolve(problem)
    hPrecision <- result$getValue(hPrecision)
  } else{
    glasso.res <- glasso::glasso(var(X), rho = 0)
    hPrecision <- glasso.res$wi
    hPrecision <- (t(hPrecision) + hPrecision)/2
    for(pair_i in which(!model)){
      pair <- index_to_pair(pair_i, p)
      i <- pair$i
      j <- pair$j
      hPrecision[i, j] <- 0
      hPrecision[j, i] <- 0
    }
    eigvals <- eigen(hPrecision)$values
    min_eig <- min(eigvals)
    if(min_eig <= 0){
      hPrecision <- hPrecision + ((max(eigvals)-min_eig)*1e-5 - min_eig) * diag(p) # make condition number 1e5
    }
  }
  return(hPrecision)
}
