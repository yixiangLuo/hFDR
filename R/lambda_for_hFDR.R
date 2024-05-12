lambda_for_hFDR.lasso <- function(X, y, psi = "pval", nlambda = 40, trunc_range = 0.1){
  n <- NROW(X)
  p <- NCOL(X)

  if(!is.function(psi) && psi == "pval") psi <- psi.guasslinear
  debias <- psi(X, y, 1:p)
  psi.val <- debias$psi
  normalizer <- debias$normalizer
  if(any(psi.val < 0)) stop('`psi` return a value less than zero.')
  if(is.null(normalizer)){
    stop('`psi` must return a non-null `normalizer` to use the built-in fast algorithm for lasso or fs. Otherwise please pass a user-specified function to the `selection` argument.')
  }


  data.pack <- preprocess.gausslinear(X, y, "lasso")

  Sigma <- t(X) %*% X
  lasso.fit <- glmnet::glmnet(X, y,
                              intercept = T, standardize = T,
                              family = "gaussian")
  nnodes <- length(lasso.fit$lambda)

  hFDRj.all.approx <- matrix(NA, p, nnodes)
  n.precise_calc <- rep(NA, nnodes)

  for(lambda_i in 1:nnodes){
    selected <- sort(which(abs(lasso.fit$beta[, lambda_i]) > 1e-8))
    unselected <- which(!(1:p %in% selected))
    lasso.pack <- list(beta_lasso = lasso.fit$beta[, lambda_i],
                       Sigma = Sigma,
                       lambda = lasso.fit$lambda[lambda_i] * n)

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

    hFDRj.all.approx[, lambda_i] <- hFDRj.approx
    n.precise_calc[lambda_i] <- (length(j.exact) + sum(j.approx.prob)) / lasso.fit$lambda[lambda_i]
  }

  hFDR.approx <- colSums(hFDRj.all.approx)
  importance <- sapply(1:nnodes, function(i){abs(diff(range(pmin(1, hFDR.approx[i:nnodes]))))})
  # plot(importance)
  # lines(n.precise_calc / max(n.precise_calc) * max(importance))
  end <- max(which(importance >= trunc_range))

  lambda_max <- lasso.fit$lambda[1]
  lambda_min <- lasso.fit$lambda[end]
  lambda <- lambda_max * (lambda_min/lambda_max)^((0:nlambda)/nlambda)

  return(lambda)
}
