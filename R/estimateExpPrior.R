#  estimateExpPrior.R
#' Define/estimate normal multivariate prior pdf for exponential decay parameters.
#' @param x a numeric vector of depths
#' @param uy a numeric vector of uncertainties
#' @param dataType an numeric (1 or 2) defining the type of data 
#' @param priorType a string defining the type of prior
#' @param out ourput from \code{fitMonoExp}
#' @param ru_theta optional real defining the relative uncertainty on parameters
#' @return A list containing the center and covariance matrix of the prior pdf.  
#' @author Pascal PERNOT
#' @export

estimateExpPrior <- function(x, uy, dataType,
                             priorType='mono',
                             out, ru_theta = 0.05) {
  statsObs = function(x){
    # Stats for residuals are Q95 and 0.
    c(quantile(abs(x),probs=c(0.95)), 0.)
  }
  statsSim = function(x) {
    # Stats for residuals are PU95 and sd(PU)
    c( 1.96*sqrt(mean(x^2,na.rm=TRUE)),sd(x,na.rm=TRUE))
  }
  upModnoCorr    = function(p,parList) {
    # Use Propagation of variances to Estimate Prediction uncertainty
    # Covariances assumed to be null
    
    x  = parList$x
    p1 = parList$p1
    dt = parList$dataType
    
    ## Build covariance matrix from params
    sig = c(p[1],p[2],p[3])
    V = diag(sig^2)
    
    ## Linear UP
    up=rep(NA,length(x))
    for (j in 1:length(x)) {
      tmp = exp(- dt*x[j] / p1[3])
      J = c(1,
            tmp,
            dt*x[j]*p1[2]/p1[3]^2 * tmp)
      up[j] = sqrt(t(J) %*% V %*% J)
    }
    return(up)
  }
  modSM = function(par,parList){
    statsSim(
      upModnoCorr(par,parList)/parList$uy
    )
  }
  chi2SM = function (par, parList) {
    sum(
      (parList$Sobs-modSM(par,parList))^2
    )
  }
  opt2ndPass = function(x,uy,dataType,Sobs=NA, p1=NA,
                        lower=NULL,upper=NULL,start=NULL){
    
    # cl <- makeCluster(4); setDefaultCluster(cl=cl)
    # clusterExport(cl,
    #               c("residSM","wgt","modSM","statsSim",
    #                 "upModnoCorr"),
    #               envir = .GlobalEnv)
    # clusterExport(cl,
    #               c("x","y","Sobs","p1","weight"),
    #               envir = environment())
    
    # Bounds
    if(is.null(lower))
      lower = c(0,0,0)
    if(is.null(upper))
      upper = p1
    if(is.null(start))
      start = 0.5*(lower+upper)
    
    parList = list(
      x        = x,
      uy       = uy,
      dataType = dataType,
      p1       = p1,
      Sobs     = Sobs)
    
    best = rgenoud::genoud(
      fn              = chi2SM,
      parList         = parList,
      starting.values = start,
      nvars           = length(start),
      BFGS            = TRUE,
      BFGSburnin      = 1,
      print.level     = 0,
      max.generations = 30,
      wait.generations= 5,
      gradient.check  = FALSE,
      pop.size        = 100,
      Domains         = cbind(lower,upper),
      boundary.enforcement = 2,
      # cluster=cl,
      hessian         = FALSE
    )
    
    # setDefaultCluster(cl=NULL); stopCluster(cl)
    
    return(
      list(par     = best$par,
           hessian = best$hessian,
           parList = parList
      )
    )
  }
  fitUqOpt = function(x,uy,dataType,theta0,resid,sigFac=0.3) {
    
    # Estimate uncertainty by Statistics Matching model
    
    ## Data statistics
    S.obs = statsObs(resid/uy)
    
    ## Initial estimate
    start = sigFac/2 * theta0
    
    # Optimize
    bestu = opt2ndPass(
      x, uy, dataType,
      Sobs=S.obs, p1=theta0,
      upper=2*start,start=start
    )
    
    ## Posterior Prediction
    P.map = P.mean = c(bestu$par,0,0,0)
    
    ## Results summary
    S.map = modSM(bestu$par,bestu$parList)
    
    return(
      list(
        S.obs     = S.obs,
        S.map     = S.map,
        P.map     = P.map
      )
    )
  }
  
  theta0    = out$best.theta   # Used by next stage
  
  if(priorType == 'mono') {
    # Scale Monoexp covariance matrix by ru_theta
    cor_theta = out$cor.theta
    u_theta   = ru_theta * theta0
    
  } else {
    # ABC nocorr approx.
    cor_theta = diag(c(1,1,1))
    P = fitUqOpt(
      x, uy, dataType, theta0,
      resid=out$fit$par$resid
    )
    u_theta   = P$P.map[1:3]
  }
  
  Sigma0    = diag(u_theta) %*% cor_theta %*% diag(u_theta)
  
  return(
    list(
      theta0 = theta0,
      Sigma0 = Sigma0
    )
  )
}
