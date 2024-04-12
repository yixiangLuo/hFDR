parallel.prepare <- function(n_cores){
  if(n_cores > 1 && !requireNamespace("doParallel", quietly=T)) {
    warning("doParallel is not installed. Using sequential computing instead.", call. = F, immediate. = T)
    n_cores <- 1
  }
  if(n_cores > 1 && !requireNamespace("parallel", quietly=T)) {
    warning("parallel is not installed. Using sequential computing instead.", call. = F, immediate. = T)
    n_cores <- 1
  }
  if(n_cores > 1 && !requireNamespace("foreach", quietly=T)) {
    warning("foreach is not installed. Using sequential computing instead.", call. = F, immediate. = T)
    n_cores <- 1
  }
  if(n_cores > 1){
    cores_avail <- parallel::detectCores(all.tests = TRUE, logical = TRUE)
    if(n_cores > cores_avail){
      warning(paste("The requested number of cores is not available. Using instead", cores_avail, "cores."), immediate. = T)
      n_cores <- cores_avail
    }
  }
  if(n_cores > 1){
    doParallel::registerDoParallel(n_cores)

    forall <- foreach::foreach
    `%exec%` <- foreach::`%dopar%`
  } else{
    forall <- iterate_seq
    `%exec%` <- `%do_seq%`
  }

  return(list(iterator = forall, connector = `%exec%`))
}

# sequence generator for a sequential iterator,
# need such format similar to the one used in foreach
iterate_seq <- function(...){
  arg <- list(...)[1]
  argname <- names(arg)
  arg <- arg[[1]]

  return(list(arg = arg, argname = argname))
}

# binary operator for a sequential iterator,
# need such format similar to the one used in foreach
`%do_seq%` <- function(obj, expr){
  result <- NULL
  if(length(obj$arg) > 0){
    for(item in obj$arg){
      envir <- parent.frame()
      envir[[obj$argname]] <- item
      result <- c(result, list(eval(substitute(expr), envir = envir)))
    }
  }
  return(result)
}
