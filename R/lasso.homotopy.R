
lasso.homotopy <- function(lasso.pack, j, Xjy_range){

  Sigma <- lasso.pack$Sigma

  p <- length(lasso.pack$Xy)
  lasso.pack$beta_lasso <- as.vector(lasso.pack$beta_lasso)

  selected <- lasso.pack$selected
  not_selected <- setdiff(1:p, selected)
  if(length(selected) == 0){
    lasso.pack$Sigma_selected_inv <- NULL
  } else if(is.null(lasso.pack$Sigma_selected_inv)){
    lasso.pack$Sigma_selected_inv <- solve(Sigma[selected, selected])
  }

  Xjy_nodes <- c()
  beta_at_nodes <- c()
  selected.list <- list()
  Sigma_selected_inv.list <- list()

  for(direction in c(1, -1)){
    stop <- F
    lasso.pack.cur <- lasso.pack

    while (!stop) {
      if(direction < 0){
        selected.list <- c(selected.list, list(lasso.pack.cur$selected))
        Sigma_selected_inv.list <- c(Sigma_selected_inv.list, list(lasso.pack.cur$Sigma_selected_inv))
      }

      res_at_node <- lasso.homotopy.step(lasso.pack.cur, j, Xjy_range, direction)
      Xjy_nodes <- c(Xjy_nodes, res_at_node$Xjy)
      beta_at_nodes <- cbind(beta_at_nodes, res_at_node$beta_lasso)

      stop <- res_at_node$stop
      if(!stop){
        lasso.pack.cur$X <- lasso.pack$X
        lasso.pack.cur$Sigma_selected_inv <- update_selectedInv(Sigma,
                                                                lasso.pack.cur$Sigma_selected_inv,
                                                                lasso.pack.cur$selected,
                                                                res_at_node$selected)
        lasso.pack.cur$Xy[j] <- res_at_node$Xjy
        lasso.pack.cur$beta_lasso <- res_at_node$beta_lasso
        lasso.pack.cur$selected <- res_at_node$selected
      }

      if(direction > 0){
        selected.list <- c(selected.list, list(lasso.pack.cur$selected))
        Sigma_selected_inv.list <- c(Sigma_selected_inv.list, list(lasso.pack.cur$Sigma_selected_inv))
      }
    }
  }

  order_index <- order(Xjy_nodes)
  Xjy_nodes <- Xjy_nodes[order_index]
  beta_at_nodes <- beta_at_nodes[, order_index]
  selected.list <- selected.list[order_index]
  Sigma_selected_inv.list <- Sigma_selected_inv.list[order_index]

  lasso.homopath <- structure(list(Xjy_nodes = Xjy_nodes,
                                   beta_at_nodes = beta_at_nodes,
                                   selected_at_nodes_right = selected.list,
                                   SSI_at_nodes_right = Sigma_selected_inv.list),
                              class = 'lasso.homopath')

  return(lasso.homopath)
}

lasso.homotopy.preprocess <- function(X, y, lambda, tol = 1e-7){
  n <- NROW(X)
  p <- NCOL(X)

  Sigma <- t(X) %*% X
  Xy <- c(t(y) %*% X)

  lasso.fit <- glmnet::glmnet(X, y, lambda = lambda,
                              intercept = T, standardize = T,
                              family = "gaussian")

  lambdas <- lasso.fit$lambda * n
  nlambda <- length(lambdas)

  selected.list <- list()
  Sigma_selected_inv.list <- list()
  pre_selected <- NA
  for(lambda_i in 1:nlambda){
    selected <- sort(which(abs(lasso.fit$beta[, lambda_i]) > tol))

    if(!identical(pre_selected, selected)){
      Sigma_selected_inv <- if(length(selected) > 0){
        list(solve(Sigma[selected, selected]))
      } else{ list(NULL) }
    }

    Sigma_selected_inv.list[lambda_i] <- Sigma_selected_inv
    selected.list[[lambda_i]] <- selected

    pre_selected <- selected
  }

  return(list(Sigma = Sigma, Xy = Xy,
              lambdas = lambdas, betas = lasso.fit$beta,
              selected.list = selected.list,
              Sigma_selected_inv.list = Sigma_selected_inv.list,
              tol = tol))
}


lasso.homotopy.step <- function(lasso.pack, j, Xjy_range, direction){

  Xy <- lasso.pack$Xy
  beta_lasso <- lasso.pack$beta_lasso
  selected <- lasso.pack$selected
  Sigma <- lasso.pack$Sigma
  Sigma_selected_inv <- lasso.pack$Sigma_selected_inv
  lambda <- lasso.pack$lambda
  tol <- lasso.pack$tol

  beta_lasso <- as.vector(beta_lasso)
  if(!is.null(Sigma_selected_inv)){
    Sigma_selected_inv <- matrix(Sigma_selected_inv, nrow = NROW(Sigma_selected_inv))
  }
  direction <- sign(direction)
  p <- NCOL(Sigma)
  n_sel <- length(selected)
  not_selected <- which(!(1:p %in% selected))

  if(n_sel > 0){
    j_index <- which(selected == j)
    slope.beta <- if(length(j_index) > 0){
      c(Sigma_selected_inv[, j_index])
    } else{
      rep(0, n_sel)
    }
    intercept.beta <- beta_lasso[selected]
    PoC.beta <- (0 - intercept.beta) / slope.beta
    PoC.beta[is.na(PoC.beta)] <- 0
    PoC.beta <- dir_trunc(PoC.beta, direction)
  } else{
    PoC.beta <- NULL
  }

  if(n_sel < p){
    slope.lambda <- (not_selected == j) -
      ifelse(n_sel > 0, matrix(Sigma[not_selected, selected], ncol = n_sel)
             %*% slope.beta, 0)
    intercept.lambda <- Xy[not_selected] - Sigma[not_selected, ] %*% beta_lasso
    intercept.lambda <- sign(intercept.lambda) * pmin(abs(intercept.lambda), lambda*(1-tol))
    PoC.lambda <- (sign(slope.lambda) * direction * lambda - intercept.lambda) / slope.lambda
    PoC.lambda[is.na(PoC.lambda)] <- 0
    PoC.lambda <- dir_trunc(PoC.lambda, direction)
  } else{
    PoC.lambda <- NULL
  }

  PoC <- rep(NA, p)
  PoC[c(selected, not_selected)] <- c(PoC.beta, PoC.lambda)
  PoC_index <- which.min(abs(PoC))
  PoC <- PoC[PoC_index]

  if(abs(PoC) < Inf){
    exclude <- which(selected %in% PoC_index)
    selected.next <- if(length(exclude) > 0){
      selected[-exclude]
    } else{
      insert_index <- sum(selected < PoC_index)
      append(selected, PoC_index, after = insert_index)
    }
    next.Xjy <- Xy[j] + PoC
  } else{
    selected.next <- selected
    next.Xjy <- ifelse(PoC > 0, max(Xjy_range), min(Xjy_range))
    PoC <- next.Xjy - Xy[j]
  }
  next.beta_lasso <- beta_lasso
  if(n_sel > 0){
    next.beta_lasso[selected] <- beta_lasso[selected] + slope.beta * PoC
    next.beta_lasso[abs(next.beta_lasso) <= tol] <- 0
  }
  if(any(is.na(next.beta_lasso))) browser()

  stop <- if(direction > 0){
    next.Xjy >= max(Xjy_range)
  } else{
    next.Xjy <= min(Xjy_range)
  }

  X <- lasso.pack$X

  return(list(Xjy = next.Xjy, beta_lasso = next.beta_lasso,
              selected = selected.next, stop = stop))
}

update_selectedInv <- function(Sigma, Sigma_selected_inv, selected, selected.next){
  if(length(selected.next) == 0){
    return(NULL)
  }
  if(length(selected) == 0){
    return(solve(Sigma[selected.next, selected.next]))
  }

  include <- selected.next[!(selected.next %in% selected)]
  exclude <- selected[!(selected %in% selected.next)]

  if(length(include) + length(exclude) != 1){
    stop()
  }
  if(length(include) > 0){
    l <- length(selected)
    update.index <- order(c(selected, include))
    ## inverse of [A, b; b^T, d] by Woodbury formula
    org.inv <- cbind(rbind(Sigma_selected_inv, rep(0, l)), c(rep(0, l), 1/Sigma[include, include]))
    U <- matrix(c(rep(0, l), 1, Sigma[selected, include], 0), ncol = 2)
    V <- rbind(c(Sigma[selected, include], 0), c(rep(0, l), 1))
    org.inv_U <- org.inv %*% U
    update.inv <- org.inv - org.inv_U %*% solve_mat22(diag(2) + V %*% org.inv_U) %*% (V %*% org.inv)
    update.inv <- update.inv[update.index, update.index]
  } else if(length(exclude) > 0){
    l <- length(selected.next)
    update.index <- rank(c(selected.next, exclude))
    ## inverse of [A, 0; 0, d] by Woodbury formula
    org.inv <- Sigma_selected_inv[update.index, update.index]
    U <- matrix(c(rep(0, l), 1, -Sigma[selected.next, exclude], 0), ncol = 2)
    V <- rbind(c(-Sigma[selected.next, exclude], 0), c(rep(0, l), 1))
    org.inv_U <- org.inv %*% U
    update.inv <- org.inv - org.inv_U %*% solve_mat22(diag(2) + V %*% org.inv_U) %*% (V %*% org.inv)
    update.inv <- update.inv[1:l, 1:l]
  }

  return(update.inv)
}

solve_mat22 <- function(mat){
  a <- mat[1, 1]
  b <- mat[1, 2]
  c <- mat[2, 1]
  d <- mat[2, 2]
  det <- a*d-b*c
  inv <- matrix(c(d, -c, -b, a), ncol = 2) / det
  return(inv)
}

