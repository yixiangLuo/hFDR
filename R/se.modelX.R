
se.modelX <- function(X, y, select, lambda, psi, sampler.modelX, predict = NULL, n_sample.hfdr, n_sample.se, n_cores = 1){
  p <- NCOL(X)
  n <- NROW(X)

  if(is.function(predict)){
    cv.res <- cv.model(X, y, lambda, predict, nfold = 10)
    model <- c(select(X, y, cv.res$lambda.min))
  } else{
    model <- (psi(X, y, 1:p)$psi == 0)
  }
  nulls.guess <- which(!model)

  storage <- new.env()

  hFDR.samples <- matrix(NA, length(lambda), n_sample.se)

  for(boot_i in 1:n_sample.se){
    bs.samples <- sample(1:n, n, replace = T)
    X.boot <- X[bs.samples, ]
    if(length(nulls.guess) > 0){
      X.boot[, nulls.guess] <- sampler.modelX(X.boot, nulls.guess, 1, storage)[, , 1]
    }
    y.boot <- y[bs.samples]

    hFDR.boot <- hFDR.modelX(X.boot, y.boot, select, lambda, psi, sampler.modelX, n_sample.hfdr, n_cores)$hFDR

    hFDR.samples[, boot_i] <- hFDR.boot
  }

  rm(list = ls(envir = storage), envir = storage)

  se <- sapply(1:NROW(hFDR.samples), function(lambda_i){
    sqrt(var(hFDR.samples[lambda_i, ]))
  })

  return(se)
}
