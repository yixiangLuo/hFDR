se.gausslinear <- function(X, y, select, lambda, psi, n_sample.hfdr, n_sample.se, n_cores = 1){
  p <- NCOL(X)
  n <- NROW(X)

  predict.mle <- function(X.new, X, y){
    if(NCOL(X) > 0){
      lm.coefs <- lm(y ~ X + 1)$coefficients
      coef <- lm.coefs[-1]
      a0 <- lm.coefs[1]

      if(any(is.na(coef))){
        coef[is.na(coef)] <- 0
        warning("singular X created by cross-validation.")
      }
      return(X.new %*% coef + a0)
    } else{
      return(rep(mean(y), NROW(X.new)))
    }
  }

  select.token <- select
  if(!is.function(select)){
    if(select == "lasso"){
      select <- function(X, y, lambda){
        res <- glmnet::glmnet(X, y, lambda = lambda,
                              intercept = T, standardize = T,
                              family = "gaussian")
        as.matrix(res$beta != 0)
      }
    } else if(select == "fs"){
      select <- forward_stepwise
    }
  }

  cv.res <- cv.model(X, y, select, lambda, predict.mle, nfold = 10)
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
