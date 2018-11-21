#  printBr.R
#' Print Birges's ratio and confidence interval.

#' @param fit an object (list or stanfit) issued from a stan code
#' @param prob a probability threshold
#' @param silent a logical controling standard output

#' @return Returns invisibly a list containing
#' \describe{
#'   \item{br}{the estimated value of the Birge's ratio} 
#'   \item{ndf}{the estimated number of degrees of freedom} 
#'   \item{CI95}{a 95\% interval for admissible values of br}
#'   \item{alert}{a possible warning message, or NULL}
#' }

#' @author Pascal PERNOT

#' @details Birges's ratio and a 95% interavl is estimated for 
#' fitMonoExp or fitExpGP models. For the latter, the number of
#' degrees of freedom is based on a count of active control points
#' in the GP model. Parameter \code{prob} defines the probability
#' above which a control point is considered active by differing 
#' from zero.

#' @export

printBr         <- function(fit, prob=0.90, silent = FALSE) {
  # Extract relevant infos from fit object
  
  ## Method and model
  if(class(fit) == 'stanfit') {
    model  = fit@model_name
    method = fit@stan_args[[1]]$method
  } else {
    method = 'optim'
    if(is.null(fit$par$yGP))
      model = 'modFitMonoExp'
    else
      model = 'modFitExpGP'
  }
  
  ## Birge's ratio and Nb of residuals
  if(method == 'optim') {
    br = fit$par$br
    N  = try(length(fit$par$resid),silent=TRUE)
  } else {
    br = try(mean(extract(fit,pars='br')[[1]]),silent=TRUE)
    N  = fit@par_dims$resid
  }
  if(is.null(br)              |
     class(br) == 'try-error' |
     is.null(N)               |
     class(N) == 'try-error'    )
    return(NULL)
  
  
  ## Nb of parameters (Np) and active control points (Nn)
  if(model == 'modFitExpGP') {
    Np   = 5
    if(method == 'optim')
      Nn0  = length(fit$par$yGP)
    else
      Nn0  = fit@par_dims$yGP
    nAct = nActCtrlPts(fit,prob)
    if(is.null(nAct))
      Nn = Nn0
    else
      Nn = nAct
  } else {
    Np   = 3
    Nn0  = Nn = 0
    nAct = NULL
  }
  
  # Degrees of freedom for residuals
  ndf0 = N - (Np + Nn0) # As computed in stan model
  ndf  = N - (Np + Nn)
  
  # Adjust br for the correct degrees of frredom
  br  = br * ndf0 / ndf
  
  # Confidence interval on br
  CI95 = c(qchisq(0.025,df=ndf),qchisq(0.975,df=ndf)) / ndf
  
  alert = NULL
  if(!silent) {
    if(!is.null(nAct))
      cat('Active pts.:',Nn,'/',Nn0,'\n')
    cat('ndf        :',ndf,'\n')
    cat('br         :',signif(br,2),'\n')
    cat('CI95(br)   :',paste0(signif(CI95,2),collapse='-'),'\n')
  }
  if(prod(CI95-br) >= 0)
    alert = '!!! WARNING: br out of interval !!!'
  
  if(model == 'modFitExpGP' & is.null(nAct)) {
    # Let the user decide by himself
    for(n in rev(0:Nn0)) {
      ndf1 = N - (Np + n)
      CI951 = c(qchisq(0.025,df=ndf1),qchisq(0.975,df=ndf1)) / ndf1
      br1  = br * ndf0 / ndf1
      if(prod(CI951-br1) < 0)
        break
    }
    alert = paste0(alert,'\n',
                   '--> OK if there are less than\n',
                   n+1,' active ctrl points')
  }

  if(!silent)
    cat(alert,'\n')
  
  invisible(
    list(
      br    = br,
      nsf   = ndf,
      CI95  = CI95,
      alert = alert
    )
  )
}
