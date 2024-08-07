#' @import graphics
NULL

#' Produces a plot of the estimated FDR
#'
#' @param hFDR.obj an \code{hFDR} object generated by \code{hFDR::hFDR}.
#' @param cv.obj a \code{cv.obj} object generated by \code{cv.model} or
#' \code{cv.gaussgraph}; or a \code{cv.glmnet} object generated by
#' \code{glmnet::cv.glmnet}.
#' @param sign.lambda either plot against lambda or its negative if
#' sign.lambda=-1 (default).
#' @param log_lambda either plot against log(lambda) (default) or lambda.
#' @param log_cv either plot cv.obj in log scale or not (default).
#' @param show_legend show legend if True.
#'
#' @examples
#' p <- 100; n <- 300; k <- 15
#' X <- matrix(rnorm(n*p), n)
#' nonzero <- sample(p, k)
#' beta <- 1 * (1:p %in% nonzero)
#' y <- X %*% beta + rnorm(n)
#'
#' n_lambda <- 40
#' lambda_max <- max(abs(t(scale(X, TRUE, FALSE)) %*% scale(y, TRUE, FALSE))) / n
#' sigma <- sqrt(sum(lm(y ~ X)$residuals^2) / (n-p-1))
#' lambda_min <- sigma / sqrt(n) / 10
#' lambda <- lambda_max * (lambda_min/lambda_max)^((0:n_lambda)/n_lambda)
#'
#' # simple example for Lasso selection. See more examples in the vignettes.
#' glmnet.cv <- glmnet::cv.glmnet(X, y, alpha = 1, nfolds = 10,
#'                        intercept = TRUE, standardize = TRUE, family = "gaussian")
#'
#' hFDR.res <- hFDR(X, y, model = "gausslinear", select = "lasso",
#'                  lambda = lambda, psi = "pval")
#'
#' plot(hFDR.res, glmnet.cv)
#'
#' @export
plot.hFDR <- function(hFDR.obj, cv.obj = NULL, sign.lambda = -1,
                      log_lambda = T, log_cv = F, show_legend = !is.null(cv.obj)){
  if(requireNamespace("latex2exp", quietly=T)) {
    hFDR_text <- latex2exp::TeX("$\\widehat{FDR}$")
    lambda_text <- latex2exp::TeX("$\\lambda$")
    loglambda_text <- latex2exp::TeX("$\\log(\\lambda)$")
  } else{
    hFDR_text <- "est FDR"
    lambda_text <- "lambda"
    loglambda_text <- "log(lambda)"
  }

  trans_x <- if(log_lambda) {function(x) log(x)} else identity

  hFDR_margin <- if(is.null(hFDR.obj$hFDR.se)) 0.05 else hFDR.obj$hFDR.se
  plot.range <- range(c(hFDR.obj$hFDR-hFDR_margin, hFDR.obj$hFDR+hFDR_margin))
  lambdas <- if(is.null(cv.obj)) hFDR.obj$lambda else c(hFDR.obj$lambda, cv.obj$lambda)
  x.range <- range(sign.lambda * trans_x(lambdas))
  plot.range[1] <- max(plot.range[1], 0)
  plot.range[2] <- min(plot.range[2], 1)
  plot.args <- list(x = sign.lambda * trans_x(hFDR.obj$lambda),
                    y = hFDR.obj$hFDR,
                    xlab = "", ylab = "",
                    xlim = x.range,
                    ylim = plot.range,
                    type = "n")

  par(mar = c(4.2,4.2,3.2,4.2), cex = 1, cex.axis = 1, cex.lab = 1)
  do.call("plot", plot.args)

  leg <- c()
  lty <- c()
  lwd <- c()
  pch <- c()
  col <- c()

  if(!is.null(cv.obj)){
    ep <- 1e-8
    cv_scale <- if(log_cv) function(x) log(pmax(x,ep)) else identity
    cv.mean <- cv_scale(cv.obj$cvm)
    cv.low <- cv_scale(cv.obj$cvlo)
    cv.up <- cv_scale(cv.obj$cvup)
    cv.range <- range(c(cv.low, cv.up))
    cv.transfer <- function(cv){
      (cv - min(cv.range)) / abs(diff(cv.range)) * abs(diff(plot.range)) + min(plot.range)
    }

    cv.color <- "#1874cd"

    lines(x = sign.lambda * trans_x(cv.obj$lambda),
          y = cv.transfer(cv.mean),
          col = cv.color)
    error.bars(sign.lambda * trans_x(cv.obj$lambda),
               cv.transfer(cv.up), cv.transfer(cv.low),
               width = 0.01, col = paste0(cv.color, as.hexmode(round(0.2*255))))

    abline(v = sign.lambda * trans_x(cv.obj$lambda.min), lty = 1, col = cv.color)
    abline(v = sign.lambda * trans_x(cv.obj$lambda.1se), lty = 2, col = cv.color)

    cv.leg <- if(cv.obj$name == "Mean-Squared Error") "CV MSE" else "CV"
    leg <- c(leg, cv.leg)
    lty <- c(lty, 1)
    lwd <- c(lwd, 1)
    pch <- c(pch, 26)
    col <- c(col, cv.color)
  }

  hFDR.color <- "#FF0000"

  if(!is.null(hFDR.obj$hFDR.se)){
    error.bars(sign.lambda * trans_x(hFDR.obj$lambda),
               hFDR.obj$hFDR+hFDR.obj$hFDR.se, hFDR.obj$hFDR-hFDR.obj$hFDR.se,
               width = 0.01, paste0(hFDR.color, as.hexmode(round(0.2*255))))
  }
  points(sign.lambda*trans_x(hFDR.obj$lambda), hFDR.obj$hFDR,
         pch = 20, col = hFDR.color)

  leg <- c(leg, hFDR_text)
  lty <- c(lty, 1)
  lwd <- c(lwd, 0)
  pch <- c(pch, 19)
  col <- c(col, hFDR.color)

  if(!is.null(cv.obj)){
    if(!is.null(cv.obj$nzero)){
      axis(side = 3, at = sign.lambda*trans_x(cv.obj$lambda),
           labels = paste(cv.obj$nzero), tick = FALSE, line = -0.5)
      mtext("# variables selected", side = 3, line = 2, cex = 1.3)
    }

    cv.ticks <- if(log_cv){
      seq(exp(cv.range[1]), exp(cv.range[2]), length.out = 6)
    } else { seq(cv.range[1], cv.range[2], length.out = 6) }
    axis(side = 4, at = cv.transfer(cv_scale(cv.ticks)),
         labels = formatC(cv.ticks, format = "g", digits = 2), tick = T, line = 0)
    mtext(cv.leg, side = 4, line = 2.5, cex = 1.3)  # "-log likelihood (cv)"

    if(show_legend) suppressWarnings({
      legend.pos <- legend("bottomright", inset = 0.05,
                           legend = leg, bg = "white",
                           lty = lty, lwd = lwd,
                           pch = pch, col = col, plot = F, y.intersp=1.2)[[1]]
      cap <- min((hFDR.obj$hFDR-hFDR_margin)[sign.lambda * trans_x(hFDR.obj$lambda) >= legend.pos$left])
      floor <- max(cv.transfer(cv.up[sign.lambda * trans_x(cv.obj$lambda) >= legend.pos$left]))
      top <- if(cap > floor) (min(cap, max(hFDR.obj$hFDR+hFDR_margin)) + legend.pos$top) / 2 else legend.pos$top
      legend(x = legend.pos$left, y = top, inset = 0.05,
             legend = leg, bg = "white",
             lty = lty, lwd = lwd,
             pch = pch, col = col, y.intersp = 1.2, bty = "o")
    })
  } else if(show_legend) {
    legend("bottomright", inset = 0.05,
           legend = leg, bg = "white",
           lty = lty, lwd = lwd,
           pch = pch, col = col, y.intersp = 1.2, bty = "o")
  }


  xlab <- if(log_lambda) loglambda_text else lambda_text
  if(sign.lambda < 0) xlab <- paste("-", xlab, sep = "")
  title(xlab = xlab, ylab = "False Discovery Rate", line = 2.5, cex.lab = 1.3)
  invisible()
}

error.bars <- function(x, upper, lower, width = 0.02, ...) {
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x, upper, x, lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
  range(upper, lower)
}
