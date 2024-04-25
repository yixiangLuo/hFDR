hFDR.gaussgraph <- function(X, select, lambda, psi, n_sample.hfdr, n_cores){
  nlambda <- length(lambda)
  p <- NCOL(X)
  pair_num <- p*(p-1)/2

  debias <- psi(X, 1:pair_num)
  psi.val <- debias$psi
  normalizer <- debias$normalizer
  if(any(psi.val < 0)) stop('`psi` return a value less than zero.')
  if(is.null(normalizer)){
    normalizer <- rep(NA, p)
  }

  parallel <- parallel.prepare(n_cores)
  forall <- parallel$iterator
  `%exec%` <- parallel$connector

  hFDR_pair <- matrix(0, nrow = pair_num, ncol = nlambda)

  for(pair_i in which(psi.val > 0)){
    pair <- index_to_pair(pair_i, p)
    i <- pair$i
    j <- pair$j

    X.samples <- sampler.gaussgraph(X, i, j, n_sample.hfdr)

    mc.samples <- forall(mc_i = 1:n_sample.hfdr, .options.multicore = list(preschedule = T)) %exec% {
      X.sample <- X
      X.sample[, i] <- X.samples$Xi[, mc_i]
      X.sample[, j] <- X.samples$Xj[, mc_i]
      res.sample <- select(X.sample, lambda)
      FDPj.sample <- res.sample[pair_i, ] / pmax(1, colSums(res.sample))

      if(is.na(normalizer[pair_i])){
        psi_pair.sample <- psi(X.sample, pair_i)$psi
      } else{
        psi_pair.sample <- NA
      }
      return(list(FDPj.sample = FDPj.sample, psi_pair.sample = psi_pair.sample))
    }
    FDP_pair.samples <- sapply(mc.samples, function(mc.sample){ mc.sample$FDPj.sample })
    psi_pair.samples <- sapply(mc.samples, function(mc.sample){ mc.sample$psi_pair.sample })

    if(is.na(normalizer[pair_i])){
      normalizer[pair_i] <- mean(psi_pair.samples, na.rm = T)
      if(normalizer[pair_i] == 0){
        warning('Monte-Carlo simulation with finite samples finds the `psi` normalizer is 0.
                  Please consider increase n_sample.hfdr until this warning disappear.')
        normalizer[pair_i] <- NA
      }
    }
    hFDR_pair[pair_i, ] <- rowMeans(FDP_pair.samples, na.rm = T) * psi.val[pair_i] / normalizer[pair_i]
  }
  hFDR <- colSums(hFDR_pair)

  hFDR.result <- list(hFDR = hFDR, hFDR.decompose = hFDR_pair)

  return(hFDR.result)
}

