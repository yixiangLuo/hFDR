#' @import graphics
NULL

#' The plot function
#'
#' Document will be ready soon
#'
#' @export
plot.hFDR <- function(hFDR.obj, cv.obj, sign.tune = -1, log_lambda = T, log_cv = F,
                      show_legend = T){
  trans_x <- if(log_lambda) {function(x) log(x)} else identity

  plot.range <- range(c(hFDR.obj$hFDR-hFDR.obj$hFDR.se, hFDR.obj$hFDR+hFDR.obj$hFDR.se))
  plot.range[1] <- max(plot.range[1], 0)
  plot.range[2] <- min(plot.range[2], 1)
  plot.args <- list(x = sign.tune * trans_x(hFDR.obj$lambda),
                    y = hFDR.obj$hFDR,
                    xlab = "", ylab = "",
                    xlim = range(sign.tune * trans_x(cv.obj$lambda)),
                    ylim = plot.range,
                    type = "n")

  # use log scale for MSE
  ep <- 1e-8
  cv_scale <- if(log_cv) function(x) log(pmax(x,ep)) else identity
  cv.mean <- cv_scale(cv.obj$cvm)
  cv.low <- cv_scale(cv.obj$cvlo)
  cv.up <- cv_scale(cv.obj$cvup)
  cv.range <- range(c(cv.low, cv.up))
  cv.transfer <- function(cv){
    (cv - min(cv.range)) / abs(diff(cv.range)) * abs(diff(plot.range)) + min(plot.range)
  }

  leg <- c()
  lty <- c()
  lwd <- c()
  pch <- c()
  col <- c()

  par(mar = c(4,4,3,4), cex = 1, cex.axis = 1, cex.lab = 1)
  do.call("plot", plot.args)

  cv.color <- "#1874cd"

  lines(x = sign.tune * trans_x(cv.obj$lambda),
        y = cv.transfer(cv.mean),
        col = cv.color)
  error.bars(sign.tune * trans_x(cv.obj$lambda),
             cv.transfer(cv.up), cv.transfer(cv.low),
             width = 0.01, col = paste0(cv.color, as.hexmode(round(0.2*255))))

  abline(v = sign.tune * trans_x(cv.obj$lambda[cv.obj$ind.min]), lty = 1, col = cv.color)
  abline(v = sign.tune * trans_x(cv.obj$lambda[cv.obj$ind.1se]), lty = 2, col = cv.color)

  cv.leg <- if(cv.obj$name == "Mean-Squared Error") "CV MSE" else "CV"
  leg <- c(leg, cv.leg)
  lty <- c(lty, 1)
  lwd <- c(lwd, 1)
  pch <- c(pch, 26)
  col <- c(col, cv.color)


  hFDR.color <- "#FF0000"

  error.bars(sign.tune * trans_x(hFDR.obj$lambda),
             hFDR.obj$hFDR+hFDR.obj$hFDR.se, hFDR.obj$hFDR-hFDR.obj$hFDR.se,
             width = 0.01, paste0(hFDR.color, as.hexmode(round(0.2*255))))
  points(sign.tune*trans_x(hFDR.obj$lambda), hFDR.obj$hFDR,
         pch = 20, col = hFDR.color)

  leg <- c(leg, "Est FDR")
  lty <- c(lty, 1)
  lwd <- c(lwd, 0)
  pch <- c(pch, 19)
  col <- c(col, hFDR.color)

  if(!is.na(cv.obj$nzero)){
    axis(side = 3, at = sign.tune*trans_x(cv.obj$lambda),
         labels = paste(cv.obj$nzero), tick = FALSE, line = -0.5)
    mtext("# selections", side = 3, line = 2)
  }

  cv.ticks <- if(log_cv){
    seq(exp(cv.range[1]), exp(cv.range[2]), length.out = 6)
  } else { seq(cv.range[1], cv.range[2], length.out = 6) }
  axis(side = 4, at = cv.transfer(cv_scale(cv.ticks)),
       labels = formatC(cv.ticks, format = "g", digits = 2), tick = T, line = 0)
  mtext(cv.leg, side = 4, line = 2.5)  # "-log likelihood (cv)"

  legend.pos <- legend("bottomright", inset = 0.05,
                       legend = leg, bg = "white",
                       lty = lty, lwd = lwd,
                       pch = pch, col = col, plot = F, y.intersp=1.2)[[1]]
  cap <- min((hFDR.obj$hFDR-hFDR.obj$hFDR.se)[sign.tune * trans_x(hFDR.obj$lambda) >= legend.pos$left])
  floor <- max(cv.transfer(cv.up[sign.tune * trans_x(cv.obj$lambda) >= legend.pos$left]))
  top <- if(cap > floor) (legend.pos$top + cap) / 2 else legend.pos$top
  if(show_legend) {
    legend(x = legend.pos$left, y = top, inset = 0.05,
           legend = leg, bg = "white",
           lty = lty, lwd = lwd,
           pch = pch, col = col, y.intersp = 1.2, bty = "o")
  }

  par(mar = c(4,4,3,4), cex = 1, cex.axis = 1, cex.lab = 1)
  xlab <- if(log_lambda) "log(lambda)" else "lambda"
  if(sign.tune < 0) xlab <- paste("-", xlab, sep = "")
  title(xlab = xlab, ylab = "False Discovery Rate", line = 2.5, cex.lab = 1)
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
