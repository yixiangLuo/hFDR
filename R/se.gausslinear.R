se.gausslinear <- function(X, y, select, lambda, psi, n_sample.hfdr, n_sample.se, n_cores = 1){
  p <- NCOL(X)
  n <- NROW(X)

  select.token <- select
  if(!is.function(select)){
    if(select == "lasso"){
      select <- select.lasso
    } else if(select == "fs"){
      select <- select.fs
    }
  }

  pred_fit <- function(X.new, X, y, lambda){
    pred_fit.model(X.new, X, y, select, lambda, pred_fit.mle.gausslinear)
  }

  cv.res <- cv.model(X, y, lambda, pred_fit, nfold = 10)
  model <- c(select(X, y, cv.res$lambda.min))
  beta.star <- rep(0, p)
  if(sum(model) > 0){
    lm.coefs <- lm(y ~ X[, model] + 1)$coefficients
    beta.star[model] <- lm.coefs[-1]
    a0 <- lm.coefs[1]
  } else{
    a0 <- mean(y)
  }

  RSS <- sum(lm(y ~ X + 1)$residuals^2)
  sigma.star <- sqrt(RSS / (n-p-1))

  hFDR.samples <- matrix(NA, length(lambda), n_sample.se)

  for(boot_i in 1:n_sample.se){
    X.boot <- X
    y.boot <- X %*% beta.star + a0 + rnorm(n, sd = sigma.star)

    hFDR.boot <- hFDR.gausslinear(X.boot, y.boot, select.token, lambda, psi, n_sample.hfdr, n_cores)$hFDR

    hFDR.samples[, boot_i] <- hFDR.boot
  }

  se <- sapply(1:NROW(hFDR.samples), function(lambda_i){
    sqrt(var(hFDR.samples[lambda_i, ]))
  })

  return(se)
}
