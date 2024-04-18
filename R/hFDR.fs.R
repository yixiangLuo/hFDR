hFDR.fs <- function(data.pack, lambda, psi){
  X <- data.pack$X
  y <- data.pack$y

  n <- NROW(X)
  p <- NCOL(X)

  n_sels <- lambda
  max_step <- max(n_sels)

  debias <- psi(X, y, 1:p)
  psi.val <- debias$psi
  normalizer <- debias$normalizer
  if(any(psi.val < 0)) stop('`psi` return a value less than zero.')
  if(is.null(normalizer)){
    stop('`psi` must return a non-null `normalizer` to use the built-in fast algorithm for lasso or fs. Otherwise please pass a user-specified function to the `selection` argument.')
  }

  candidates <- 1:p
  thresholds_up <- matrix(NA, nrow = p, ncol = max_step)
  thresholds_low <- matrix(NA, nrow = p, ncol = max_step)

  for(step in 1:max_step){
    res <- move_forward(X, y, candidates, move = step < max_step)
    selected <- candidates[res$cand_sel]
    for(cand in 1:length(candidates)){
      j <- candidates[cand]
      vjy_obs <- data.pack$vjy_obs[j]
      vj <- data.pack$vj_mat[, j]
      if(j != selected){
        y_proj_j <- res$y_proj[cand]
        vj_to_Xj <- sum(vj * X[, j])
        thresh_1 <- (res$y_proj[res$cand_sel] - (y_proj_j - vjy_obs*vj_to_Xj)) / vj_to_Xj
        thresh_2 <- (-res$y_proj[res$cand_sel] - (y_proj_j - vjy_obs*vj_to_Xj)) / vj_to_Xj
        thresholds_up[j, step] <- max(thresh_1, thresh_2)
        thresholds_low[j, step] <- min(thresh_1, thresh_2)
      } else if(psi.val[j] > 0){
        X_v <- X
        y_v <- y
        candidates_v <- candidates[-cand]
        for(step_v in step:max_step){
          if(step_v == p){
            thresholds_up[j, step_v] <- 0
            thresholds_low[j, step_v] <- 0
            break
          }
          res_v <- move_forward(X_v, y_v, candidates_v, aux = j, move = step_v < max_step)
          y_proj_j <- sum(X_v[, j] * y_v)
          vj_to_Xj <- sum(vj * X_v[, j])
          thresh_1 <- (res_v$y_proj[res_v$cand_sel] - (y_proj_j - vjy_obs*vj_to_Xj)) / vj_to_Xj
          thresh_2 <- (-res_v$y_proj[res_v$cand_sel] - (y_proj_j - vjy_obs*vj_to_Xj)) / vj_to_Xj
          thresholds_up[j, step_v] <- max(thresh_1, thresh_2)
          thresholds_low[j, step_v] <- min(thresh_1, thresh_2)

          X_v <- res_v$X
          y_v <- res_v$y
          candidates_v <- res_v$candidates
        }
      }
    }
    X <- res$X
    y <- res$y
    candidates <- res$candidates
  }

  hFDRj <- matrix(0, nrow = p, ncol = length(n_sels))
  for(j in 1:p){
    for(step_i in 1:length(n_sels)){
      step <- n_sels[step_i]
      if(psi.val[j] != 0){
        select_prob_low <- vjy_CDF(max(thresholds_low[j, 1:step]), res_norm2 = data.pack$RSS_Xnoj[j], df = n-p)
        select_prob_up <- 1 - vjy_CDF(min(thresholds_up[j, 1:step]), res_norm2 = data.pack$RSS_Xnoj[j], df = n-p)
        hFDRj[j, step_i] <- min(1, select_prob_low + select_prob_up) / step * psi.val[j] / normalizer[j]
      }
    }
  }
  hFDR <- colSums(hFDRj)
  return(list(hFDR = hFDR, hFDR.decompose = hFDRj))
}


move_forward <- function(X, y, candidates, aux = NULL, move = T){
  y_proj <- as.vector(matrix(y, nrow=1) %*% X[, candidates])
  sel <- which.max(abs(y_proj))
  sel_j <- candidates[sel]

  if(move){
    candidates <- candidates[-sel]
    X_sel <- X[, sel_j]
    y <- y - sum(y * X_sel) * X_sel

    updates <- c(candidates, aux)
    for(j in c(candidates, aux)){
      X[, j] <- X[, j] - sum(X[, j] * X_sel) * X_sel
      X[, j] <- X[, j] / sqrt(sum(X[, j]^2))
    }
  }

  return(list(X = X, y = y, candidates = candidates, cand_sel = sel, y_proj = y_proj))
}






