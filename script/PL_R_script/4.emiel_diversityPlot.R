emiel.diversityPlot.family  <- function (data, samples = 1, iterations = 100, max.sample.size,
                                         step.size, col, show.legend = TRUE, legend.names, smooth.curves = TRUE, 
                                         print.diagonal = FALSE, col.diagonal = "grey", lwd.diagonal = 2,
                                         ylim, xlab, ylab, lwd.lines = 2, cex.lab = 1, cex.axis = 1, 
                                         lwd.axis = 1, cex.legend = 1.5, ...)
{
  if (missing(max.sample.size)) {
    max.sample.size <- max(apply(as.matrix(data$Data[, samples]),
                                 2, sum)) * 2
  }
  if (missing(step.size)) {
    step.size <- max.sample.size/100
  }
  if (missing(legend.names)) {
    legend.names = colnames(data$Data)[samples]
  }
  if (missing(xlab)) {
    xlab <- "Sampled reads"
  }
  if (missing(ylab)) {
    ylab <- "Number of unique gene families"
  }
  x.seq <- seq(1, max.sample.size, by = step.size)
  n.samples <- length(samples)
  if (missing(col)) {
    col = 1:n.samples
  }
  curves <- calculateRarefaction(data, samples = samples, iterations = iterations, 
                                 x.seq = x.seq)
  y.max <- max(unlist(curves))
  if (missing(ylim)) {
    ylim = c(0, y.max)
  }
  plot(NULL, xlim = c(0, max.sample.size), ylim = ylim, axes = FALSE, 
       xlab = xlab, ylab = ylab, cex.lab = cex.lab, ...)
  axis(1, cex.axis = cex.axis, lwd = lwd.axis)
  axis(2, cex.axis = cex.axis, lwd = lwd.axis)
  if (print.diagonal == TRUE) {
    abline(0, 1, col = col.diagonal, lwd = lwd.diagonal, 
           lty = 2)
  }
  for (s in 1:n.samples) {
    if (smooth.curves == TRUE) {
      lines(smooth.spline(x = x.seq, y = curves[[s]]), 
            lwd = lwd.lines, col = col[s], lty = 1)
    }
    else {
      lines(x = x.seq, y = curves[[s]], lwd = lwd.lines, 
            col = col[s], lty = 1)
    }
  }
  if (show.legend == TRUE) {
    legend(x = "bottomright", legend.names, col = col, lty = 1, 
           bty = "n", cex = cex.legend, lwd = lwd.lines)
  }
  return(curves)
}
