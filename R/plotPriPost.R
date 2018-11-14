#  plotPriPost.R
#' Plot 1 marginal prior and posterior pdf.
#' @param pri prior sample
#' @param pos posterior sample
#' @param tag name of the variable
#' @param xlim an interval
#' @param gPars a list of graphical parameters and colors
#' @return Produces a plot.
#' @author Pascal PERNOT
#' @export

plotPriPost      <- function(pri,pos,tag,xlim=range(c(pri,pos)), gPars) {
  # Extract graphical parameters
  for (n in names(gPars))
    assign(n,rlist::list.extract(gPars,n))
  
  # Plot overlapped densities for 2 samples
  d = density(pri)
  d$y = d$y/max(d$y)
  plot(d, type = 'l', col = cols[4],
       main = tag,
       xlab = '', xlim = xlim,
       ylab = 'Norm. density', ylim = c(0,1.1), yaxs = 'i')
  polygon(d$x,d$y,rev(d$w),0*d$y,col=col_tr2[4],border=NA)
  d = density(pos)
  d$y = d$y/max(d$y)
  lines(d$x,d$y,col=cols[6])
  polygon(d$x,d$y,rev(d$w),0*d$y,col=col_tr2[6],border=NA)
}
