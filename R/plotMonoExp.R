#  plotMonoExp.R
#' Plot outputs from \code{fitMonoExp}.
#' @param x a numeric vector
#' @param y a numeric vector of responses/data
#' @param uy a numeric vector of uncertainties
#' @param ySmooth a numeric vector of smoothed data
#' @param mod a numeric vector of model values
#' @param resid a numeric vector of residuals
#' @param gPars a list of graphical parameters and colors
#' @return Produces a plot.
#' @author Pascal PERNOT
#' @export

plotMonoExp     <- function(x, y, uy, ySmooth, mod, resid, gPars, dataType) {
  # Extract graphical parameters
  for (n in names(gPars))
    assign(n,rlist::list.extract(gPars,n))
  
  par(mfrow=c(1,2),pty=pty,mar=mar,mgp=mgp,tcl=tcl,lwd=lwd, cex=cex)
  
  if (dataType==1){ylabel = "mean amplitude (a.u.)"}
  if (dataType==2){ylabel = "mean intensity (a.u.)"}
  
  # Fit
  plot(x,y,pch=20,cex=0.5,col=cols[6],
       main= plot_title,
       xlab= xlabel,
       ylab= ylabel)
  grid()
  lines(x,mod,col=cols[7])
  legend('topright', bty='n',
         title = '', title.adj = 1,
         legend=c('data','best fit'),
         pch=c(20,NA),lty=c(-1,1),
         col=c(cols[6],cols[7])
  )
  legend('topright', bty='n', legend=c('','','',as.expression(bquote("br    " == .(formatC(br,digits=3)))),
                                       as.expression(bquote("CI95" == .(paste0(signif(CI95,2),collapse='-'))))
                                      )
  )
  box()
  
  # Residus
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
  lines(x, ySmooth-mod, col=cols[7])
  legend('topright', bty='n',
         title = '', title.adj = 1,
         legend=c('mean resid.','data 95% uncert.','best fit - smooth'),
         pch=c(20,NA,NA),lty=c(-1,1,1),lwd=c(1,10,2),
         col=c(cols[6],col_tr2[4],cols[7])
  )
  box()
}
