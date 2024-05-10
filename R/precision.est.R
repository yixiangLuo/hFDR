#' Estimate the precision matrix by graphical lasso
#'
#' @param X n-by-p matrix of multivariate Gaussian variables.
#' @param lambda regularity parameter sequence.
#'
#' @return the estimated precision matrix.
#'
#' @export
precision.est.glasso <- function(X, lambda){
  glasso.res <- glasso::glasso(var(X), lambda)
  hPrecision <- glasso.res$wi
  return(hPrecision)
}

# Estimate the MLE precision matrix under the constrained that Precision[i,j]=0
# if pair (i,j) is not selected by `select`
precision.est.mle.sparse <- function(X, select, lambda, precise = TRUE){
  p <- NCOL(X)
  n <- NROW(X)

  model <- select(X, lambda)

  if(precise){
    if(!requireNamespace("CVXR", quietly=T)) {
      warning("CVXR is not installed. Less accurate method will be used to find the MLE under the sparsity constraints. This would make the estimated standard error less reliable. ", call. = F, immediate. = T)
      precise <- FALSE
    } else{
      tryCatch({
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
      }, error = function(msg){
        warning(paste0("Error in CVXR optimization: ", msg, ". Less accurate method will be used to find the MLE under the sparsity constraints. This would make the estimated standard error less reliable. "), call. = F, immediate. = T)
        precise <- FALSE
      })
    }
  }

  if(!precise){
    Sigma <- var(X)
    glasso.res <- glasso::glasso(Sigma, rho = min(1, abs(Sigma))*1e-8)
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
