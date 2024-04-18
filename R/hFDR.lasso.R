
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

  hFDRj <- sapply(1:nlambda, function(lambda_i){
    selected <- homotopy.prepare$selected.list[[lambda_i]]
    lasso.pack <- list(X = X, Xy = homotopy.prepare$Xy,
                       beta_lasso = homotopy.prepare$betas[, lambda_i],
                       selected = selected,
                       Sigma = homotopy.prepare$Sigma,
                       Sigma_selected_inv = homotopy.prepare$Sigma_selected_inv.list[[lambda_i]],
                       lambda = lambdas[lambda_i],
                       tol = tol)

    hFDRj <- forall(j = 1:p, .options.multicore = list(preschedule = F)) %exec% {
      Xjy_range <- c(data.pack$Xy_bound[, j])
      if(psi.val[j] == 0){
        hFDRj <- 0
      } else{
        lasso.homopath <- lasso.homotopy(lasso.pack, j, Xjy_range)
        hFDRj <- integral.hFDRj_star(lasso.homopath, j, Xjy_range, data.pack, tol) * psi.val[j] / normalizer[j]
      }
      return(hFDRj)
    }
    hFDRj <- do.call(c, hFDRj)
  })

  hFDR <- colSums(hFDRj)

  return(list(hFDR = hFDR, hFDR.decompose = hFDRj))
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
  if(trunc.low < 1 || trunc.up < 1) browser()
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








