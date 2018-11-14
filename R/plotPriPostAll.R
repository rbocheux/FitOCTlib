#  plotPriPostAll.R
#' Plot matrix of marginal prior and posterior pdfs.
#' @param pri a stanfit object with prior samples
#' @param pos a stanfit object with posterior samples
#' @param gPars a list of graphical parameters and colors
#' @return Produces a plot.
#' @author Pascal PERNOT
#' @export

plotPriPostAll   <- function(pri, pos, gPars) {
  # Extract graphical parameters
  for (n in names(gPars))
    assign(n,rlist::list.extract(gPars,n))
  
  theta_pri   = rstan::extract(pri,'theta')[[1]]
  yGP_pri     = rstan::extract(pri,'yGP')[[1]]
  lambda_pri  = rstan::extract(pri,'lambda')[[1]]
  sigma_pri   = rstan::extract(pri,'sigma')[[1]]
  
  theta_pos   = rstan::extract(pos,'theta')[[1]]
  yGP_pos     = rstan::extract(pos,'yGP')[[1]]
  lambda_pos  = rstan::extract(pos,'lambda')[[1]]
  sigma_pos   = rstan::extract(pos,'sigma')[[1]]
  
  nPar = ncol(theta_pri) + ncol(yGP_pri) + 2
  ncol = 4
  nrow = floor(nPar / ncol)
  if(ncol*nrow < nPar) nrow = nrow + 1
  
  par(mfrow=c(nrow,ncol),pty=pty,mar=mar,
      mgp=mgp,tcl=tcl,lwd=lwd, cex=cex)
  
  for(i in 1:ncol(theta_pri))
    plotPriPost(theta_pri[,i],theta_pos[,i],paste0('theta_',i),
               gPars=gPars)
  plotPriPost(lambda_pri,lambda_pos,'lambda',
             gPars=gPars)
  plotPriPost(sigma_pri,sigma_pos,'sigma',
             gPars=gPars)
  for(i in 1:ncol(yGP_pri))
    plotPriPost(yGP_pri[,i],yGP_pos[,i],paste0('yGP_',i),0.5*c(-1,1),
               gPars=gPars)
}