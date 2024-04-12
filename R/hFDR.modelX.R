hFDR.modelX <- function(X, y, select, lambda, psi, sampler.modelX, n_sample.hfdr, n_cores){
  nlambda <- length(lambda)
  p <- NCOL(X)

  debias <- psi(X, y, 1:p)
  psi.val <- debias$psi
  normalizer <- debias$normalizer
  if(any(psi.val < 0)) stop('`psi` return a value less than zero.')
  if(is.null(normalizer)){
    normalizer <- rep(NA, p)
  }

  parallel <- parallel.prepare(n_cores)
  forall <- parallel$iterator
  `%exec%` <- parallel$connector

  hFDRj <- matrix(0, nrow = p, ncol = nlambda)

  for(j in which(psi.val > 0)){
    Xj.samples <- sampler.modelX(X, j, n_sample.hfdr)

    mc.samples <- forall(mc_i = 1:n_sample.hfdr, .options.multicore = list(preschedule = F)) %exec% {
      X.sample <- X
      X.sample[, j] <- Xj.samples[, 1, mc_i]
      res.sample <- select(X.sample, y, lambda)
      FDPj.sample <- res.sample[j, ] / pmax(1, colSums(res.sample))

      if(is.na(normalizer[j])){
        psi_j.sample <- psi(X, y.sample, j)$psi
      } else{
        psi_j.sample <- NA
      }
      return(list(FDPj.sample = FDPj.sample, psi_j.sample = psi_j.sample))
    }
    FDPj.samples <- sapply(mc.samples, function(mc.sample){ mc.sample$FDPj.sample })
    psi_j.samples <- sapply(mc.samples, function(mc.sample){ mc.sample$psi_j.sample })

    if(is.na(normalizer[j])){
      normalizer[j] <- mean(psi_j.samples, na.rm = T)
      if(normalizer[j] == 0){
        warning('Monte-Carlo simulation with finite samples finds the `psi` normalizer is 0.
                  Please consider increase n_sample.hfdr until this warning disappear.')
        normalizer[j] <- NA
      }
    }
    hFDRj[j, ] <- rowMeans(FDPj.samples, na.rm = T) * psi.val[j] / normalizer[j]
  }
  hFDR <- colSums(hFDRj)

  hFDR.result <- list(hFDR = hFDR, hFDRj = hFDRj)

  return(hFDR.result)
}





