
hFDR.lasso <- function(data.pack, lambda, psi, n_cores){
  X <- data.pack$X
  y <- data.pack$y

  n <- NROW(X)
  p <- NCOL(X)

  homotopy.prepare <- lasso.homotopy.preprocess(X, y, lambda)
  tol <- homotopy.prepare$tol

  lambdas <- homotopy.prepare$lambdas
  nlambda <- length(lambdas)

  debias <- psi(X, y, 1:p)
  psi.val <- debias$psi
  normalizer <- debias$normalizer
  if(any(psi.val < 0)) stop('`psi` return a value less than zero.')
  if(is.null(normalizer)){
    stop('`psi` must return a non-null `normalizer` to use the built-in fast algorithm for lasso or fs. Otherwise please pass a user-specified function to the `selection` argument.')
  }

  parallel <- parallel.prepare(n_cores)
  forall <- parallel$iterator
  `%exec%` <- parallel$connector

  hFDRj.all <- sapply(1:nlambda, function(lambda_i){
    selected <- homotopy.prepare$selected.list[[lambda_i]]
    unselected <- which(!(1:p %in% selected))
    lasso.pack <- list(X = X, Xy = homotopy.prepare$Xy,
                       beta_lasso = homotopy.prepare$betas[, lambda_i],
                       selected = selected,
                       Sigma = homotopy.prepare$Sigma,
                       Sigma_selected_inv = homotopy.prepare$Sigma_selected_inv.list[[lambda_i]],
                       lambda = lambdas[lambda_i],
                       tol = tol)

    # get the range for Xj^T y condtional on Sj, with an adaptive early truncation
    truncate_prob <- 0.01 * max(1, length(selected)) / p
    vy_bound <- pmax(abs(vjy_quantile(truncate_prob, data.pack$RSS_Xnoj, df = n-p)), abs(data.pack$vjy_obs))
    Xy_range <- sapply(1:p, function(j){
      sort(c(vjy_to_Xjy(-vy_bound[j], data.pack$trans, j), vjy_to_Xjy(vy_bound[j], data.pack$trans, j)))
    })

    # compute the conditional probability of being selected given Sj for those
    # variables not selected in the realized data
    sel_probs <- rep(1, p)
    if(length(unselected) > 0){
      sel_probs[unselected] <- prob.selected(unselected, Xy_range, lasso.pack, data.pack)
    }
    # approximate hFDRj by the estimated PFER / #selection * phi
    hFDRj.approx <- sel_probs / max(1, length(selected)) * psi.val / normalizer

    # devide all varaibles into three categories:
    #   j.zero: hFDRj = 0
    #   j.exact: whose hFDRj roughly takes a large proportion of hFDR. Compute their hFDRj exactly.
    #   j.approx: whose hFDRj roughly takes a small proportion of hFDR. Compute their hFDRj approximately (but keeping conservative bias).
    j.zero <- which(hFDRj.approx == 0)
    j.exact <- which((1:p) %in% selected & psi.val > 0 | hFDRj.approx/sum(hFDRj.approx) > 0.1)
    j.approx <- setdiff(1:p, c(j.zero, j.exact))
    # for each j.approx, estimate hFDRj-hFDRj.approx by a random variable Dj = (hFDRj-hFDRj.approx)/j.approx.prob
    # with probability j.approx.prob and = 0 o.w.
    j.approx.prob <- pmax((hFDRj.approx/sum(hFDRj.approx))/0.1, 0.01)[j.approx]
    # index for those Dj != 0
    calc.index <- (runif(length(j.approx)) < j.approx.prob)

    hFDRj.all <- rep(0, p)
    # compute hFDRj exactly for some of the variables
    hFDRj.exact <- forall(j = c(j.exact, j.approx[calc.index]), .options.multicore = list(preschedule = T)) %exec% {
      Xjy_range <- Xy_range[, j]
      lasso.homopath <- lasso.homotopy(lasso.pack, j, Xjy_range)
      hFDRj <- integral.hFDRj_star(lasso.homopath, j, Xjy_range, data.pack, tol) * psi.val[j] / normalizer[j]
      return(hFDRj)
    }
    if(length(c(j.exact, calc.index)) > 0){
      hFDRj.exact <- do.call(c, hFDRj.exact)
      if(length(j.exact) > 0) hFDRj.all[j.exact] <- hFDRj.exact[1:length(j.exact)]
      if(length(calc.index) > 0){
        # sum of Dj over all j.approx
        approx_err.est <- sum((hFDRj.exact[(length(j.exact)+1):length(hFDRj.exact)] - hFDRj.approx[j.approx[calc.index]]) / j.approx.prob[calc.index])
        # compute hFDRj for j.approx approximately
        hFDRj.all[j.approx] <- pmax(0, hFDRj.approx[j.approx] + approx_err.est / length(j.approx))
      }
    }

    return(hFDRj.all)
  })

  hFDR <- colSums(hFDRj.all)

  return(list(hFDR = hFDR, hFDR.decompose = hFDRj.all))
}

# compute the conditional probability of being selected given Sj for those
# variables not selected in the realized data
prob.selected <- function(unselected, Xy_range, lasso.pack, data.pack){
  n <- data.pack$n
  p <- data.pack$p
  if(length(unselected) > 0){
    Sigmabeta <- lasso.pack$Sigma[unselected, , drop = F] %*% lasso.pack$beta_lasso
    Xy_at_select <- cbind(Sigmabeta - lasso.pack$lambda, Sigmabeta + lasso.pack$lambda)
    sapply(1:length(unselected), function(ind){
      j <- unselected[ind]
      Xjy_at_select <- Xy_at_select[ind, ]
      Xjy_range <- Xy_range[, j]
      if(max(Xjy_at_select) >= max(Xjy_range) && min(Xjy_at_select) <= min(Xjy_range)){
        0
      } else{
        Xjy_at_select[Xjy_at_select <= min(Xjy_range)] <- min(Xjy_range)
        Xjy_at_select[Xjy_at_select >= max(Xjy_range)] <- max(Xjy_range)
        int_nodes <- Xjy_to_vjy(Xjy_at_select, data.pack$trans, j)
        1 - abs(diff(vjy_CDF(int_nodes, data.pack$RSS_Xnoj[j], df = n-p)))
      }
    })
  } else {NULL}
}

integral.hFDRj_star <- function(lasso.homopath, j, Xjy_range, data.pack, tol = 1e-7){
  res_norm2 <- data.pack$RSS_Xnoj[j]
  df <- data.pack$n - data.pack$p
  trans <- data.pack$trans

  n_nodes <- length(lasso.homopath$Xjy_nodes)
  mid_beta <- (lasso.homopath$beta_at_nodes[, 1:(n_nodes-1)] +
    lasso.homopath$beta_at_nodes[, 2:n_nodes]) / 2
  mid_beta <- matrix(mid_beta, ncol = n_nodes-1)
  FDPj <- (abs(mid_beta[j, ]) > tol) / pmax(colSums(abs(mid_beta) > tol), 1)

  trunc.low <- sum(lasso.homopath$Xjy_nodes <= min(Xjy_range))
  trunc.up <- sum(lasso.homopath$Xjy_nodes >= max(Xjy_range))
  main_trunk <- if(trunc.low+1 <= n_nodes-trunc.up){
    lasso.homopath$Xjy_nodes[(trunc.low+1):(n_nodes-trunc.up)]
  } else { NULL }
  int_nodes <- main_trunk
  int_nodes <- Xjy_to_vjy(int_nodes, trans, j)
  FDPj <- FDPj[trunc.low:(n_nodes-trunc.up)]

  CDF <- c(trans$Xv[j] < 0, vjy_CDF(int_nodes, res_norm2, df), trans$Xv[j] > 0)
  hFDRj_star <- sum(abs(diff(CDF)) * FDPj)

  return(hFDRj_star)
}






# if(length(selected) < p){
#   unselected <- which(!(1:p %in% selected))
#   Gbeta <- G[unselected, , drop = F] %*% homotopy.prepare$betas[, lambda_i]
#   Xy_at_select <- cbind(Gbeta-lambdas[lambda_i], Gbeta+lambdas[lambda_i])
#   for(ind in 1:length(unselected)){
#     j <- unselected[ind]
#     Xjy_at_select <- Xy_at_select[ind, ]
#     Xjy_range <- Xy_range[, j]
#     if(max(Xjy_at_select) >= max(Xjy_range) && min(Xjy_at_select) <= min(Xjy_range)){
#       sel_probs[j] <- 0
#     } else{
#       Xjy_at_select[Xjy_at_select <= min(Xjy_range)] <- min(Xjy_range)
#       Xjy_at_select[Xjy_at_select >= max(Xjy_range)] <- max(Xjy_range)
#       int_nodes <- Xjy_to_vjy(Xjy_at_select, data.pack$trans, j)
#       sel_probs[j] <- 1 - abs(diff(vjy_CDF(int_nodes, data.pack$RSS_Xnoj[j], df = n-p)))
#     }
#   }
# }
# sel_probs[psi.val == 0] <- 0


