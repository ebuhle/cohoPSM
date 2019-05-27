vioplot2 <- function (x, ..., range = 1.5, h = NULL, quantiles = NULL, ylim = NULL, 
                      names = NULL,  horizontal = FALSE, col = "lightblue", border = "black", 
                      lty = 1, lwd = 1, lty.quantile = 1, lwd.quantile = 1,
                      rectCol = "black", colMed = "white", pchMed = 19, 
                      at, add = FALSE, wex = 1, drawRect = TRUE) 
{
  if(is.data.frame(x)) {
    datas <- x
  } else if(is.matrix(x)) {
    datas <- as.data.frame(x)
  } else datas <- list(x, ...)
  n <- length(datas)
  if (missing(at)) 
    at <- 1:n
  if(length(col) == 1) col <- rep(col,n)
  if(length(border) == 1) border <- rep(border,n)
  upper <- vector(mode = "numeric", length = n)
  lower <- vector(mode = "numeric", length = n)
  q1 <- vector(mode = "numeric", length = n)
  q3 <- vector(mode = "numeric", length = n)
  med <- vector(mode = "numeric", length = n)
  if(!is.null(quantiles)) {
    qq <- vector(mode = "list", length = n)
    qrank <- vector(mode = "list", length = n)
    if(length(lty.quantile)==1) lty.quantile <- rep(lty.quantile, length(quantile))
    if(length(lwd.quantile)==1) lwd.quantile <- rep(lwd.quantile, length(quantile))
  }
  base <- vector(mode = "list", length = n)
  height <- vector(mode = "list", length = n)
  baserange <- c(Inf, -Inf)
  args <- list(display = "none")
  if (!(is.null(h))) 
    args <- c(args, h = h)
  for (i in 1:n) {
    data <- datas[[i]]
    # data.min <- min(data)
    # data.max <- max(data)
    q1[i] <- quantile(data, 0.25)
    q3[i] <- quantile(data, 0.75)
    med[i] <- median(data)
    iqd <- q3[i] - q1[i]
    # upper[i] <- min(q3[i] + range * iqd, data.max)
    # lower[i] <- max(q1[i] - range * iqd, data.min)
    # est.xlim <- c(min(lower[i], data.min), max(upper[i], data.max))
    eval.points.test <- seq(min(data), max(data), length = 100)
    smout.test <- do.call("sm.density", 
                          c(list(data, xlim = range(data), eval.points = eval.points.test), args))
    height.test <- smout.test$estimate / max(smout.test$estimate)
    est.xlim <- range(quantile(data, c(0.001,0.999)), 
                      smout.test$eval.points[which(height.test >= 0.01)])
    eval.points <- seq(est.xlim[1], est.xlim[2], length  = 100)
    if(!is.null(quantiles)) {
      qq[[i]] <- quantile(data, quantiles)
      eval.points <- c(qq[[i]], eval.points)
      qrank[[i]] <- rank(eval.points, ties = "first")[1:length(quantiles)]
      eval.points <- sort(eval.points)
    }
    smout <- do.call("sm.density", 
                     c(list(data, xlim = est.xlim, eval.points = eval.points), args))
    hscale <- 0.4/max(smout$estimate) * wex
    base[[i]] <- smout$eval.points
    height[[i]] <- smout$estimate * hscale
    t <- range(base[[i]])
    baserange[1] <- min(baserange[1], t[1])
    baserange[2] <- max(baserange[2], t[2])
  }
  if (!add) {
    xlim <- if (n == 1) 
      at + c(-0.5, 0.5)
    else range(at) + min(diff(at))/2 * c(-1, 1)
    if (is.null(ylim)) {
      ylim <- baserange
    }
  }
  if (is.null(names)) {
    label <- 1:n
  }
  else {
    label <- names
  }
  boxwidth <- 0.05 * wex
  if (!add) 
    plot.new()
  if (!horizontal) {
    if (!add) {
      plot.window(xlim = xlim, ylim = ylim)
      axis(2)
      axis(1, at = at, label = label)
    }
    box()
    for (i in 1:n) {
      polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])), 
              c(base[[i]], rev(base[[i]])), 
              col = col[i], border = border[i], lty = lty, lwd = lwd)
      if (drawRect) {
        lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd, lty = lty)
        rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2, q3[i], col = rectCol)
        points(at[i], med[i], pch = pchMed, col = colMed)
      }
      if(!is.null(quantiles)) {
        segments(x0 = at[i] - height[[i]][qrank[[i]]],
                 x1 = at[i] + height[[i]][qrank[[i]]],
                 y0 = base[[i]][qrank[[i]]], col = border[i],
                 lty = lty.quantile, lwd = lwd.quantile)
      }
    }
  }
  else {
    if (!add) {
      plot.window(xlim = ylim, ylim = xlim)
      axis(1)
      axis(2, at = at, label = label)
    }
    box()
    for (i in 1:n) {
      polygon(c(base[[i]], rev(base[[i]])), 
              c(at[i] - height[[i]], rev(at[i] + height[[i]])), 
              col = col[i], border = border[i], lty = lty, lwd = lwd)
      if (drawRect) {
        lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd, lty = lty)
        rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] + boxwidth/2, col = rectCol)
        points(med[i], at[i], pch = pchMed, col = colMed)
      }
      if(!is.null(quantiles)) {
        segments(y0 = at[i] - height[[i]][qrank[[i]]],
                 y1 = at[i] + height[[i]][qrank[[i]]],
                 x0 = base[[i]][qrank[[i]]], col = border[i],
                 lty = lty.quantile, lwd = lwd.quantile)
      }
    }
  }
  invisible(list(upper = upper, lower = lower, median = med))
}
