#  plotNoise.R
#' Plot outputs from \code{estimaNoise}.
#' @param x a numeric vector
#' @param y a numeric vector of responses
#' @param uy a numeric vector of uncertainties
#' @param ySmooth a numeric vector of smoothed data
#' @param gPars a list of graphical parameters and colors
#' @return Produces a plot.
#' @author Pascal PERNOT
#' @export

plotNoise       <- function(x, y, uy, ySmooth, gPars, dataType) {
  # Extract graphical parameters
  for (n in names(gPars))
    assign(n,rlist::list.extract(gPars,n))
  
  par(mfrow=c(1,2),pty=pty,mar=mar,mgp=mgp,tcl=tcl,lwd=lwd, cex=cex)
  
  if (dataType==1){
    SNR = 20*log10(mean(y)/mean(abs(y-ySmooth)))
    ylabel = "mean amplitude (a.u.)"
  }
  if (dataType==2){
    SNR = 10*log10(mean(y)/mean(abs(y-ySmooth)))
    ylabel = "mean intensity (a.u.)"
  }
  
  # Smooth Fit
  plot(x,y,pch=20,cex=0.5,col=cols[6],
       main= plot_title,
       xlab= xlabel,
       ylab= ylabel)
  grid()
  lines(x,ySmooth,col=cols[7])
  
  legend('topright', bty='n',
         title = '', title.adj = 1,
         legend=c('data','smoother',as.expression(bquote("SNR" == .(formatC(SNR,digits=3))~"dB"))),
         pch=c(20,NA),lty=c(-1,1),
         col=c(cols[6],cols[7])
  )
  box()
  
  # Residuals
  resid   = y - ySmooth
  ylim=1.2*max(abs(resid))*c(-1,1)
  res = resid
  plot(x,res,type='n',
       ylim=ylim, main='Residuals',
       xlab= xlabel,
       ylab='residuals (a.u.)')
  grid()
  abline(h=0)
  polygon(c(x,rev(x)),c(-2*uy,rev(2*uy)),col=col_tr2[4],border = NA)
  points(x,res,pch=20,cex=0.75,col=cols[6])
  legend('topright', bty='n',
         title = '', title.adj = 1,
         legend=c('resid.','data 95% uncert.'),
         pch=c(20,NA),lty=c(-1,1),lwd=c(1,10),
         col=c(cols[6],col_tr2[4])
  )
  box()
}
