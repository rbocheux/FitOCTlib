#  estimateExpPrior.R
#' Define/estimate normal multivariate prior pdf for exponential decay parameters.
#' @param x a numeric vector of depths
#' @param uy a numeric vector of uncertainties
#' @param dataType an numeric (1 or 2) defining the type of data 
#' @param priorType a string defining the type of prior
#' @param out ourput from \code{fitMonoExp}
#' @param ru_theta optional real defining the relative uncertainty on parameters
#' @param eps tolerance parameter for moments matching method
#' @param nb_chains  number of MCMC chains
#' @param nb_warmup number of warmup steps
#' @param nb_iter   number of steps
#'  
#' @return A list containing the center, covariance matrix 
#' of the prior pdf and the constraint and realized statistics.  
#' 
#' @author Pascal PERNOT
#' 
#' @details Provides two ways to buil a prior for the exponential decay 
#' parameters of the \code{ExpGP} model:
#' \describe{
#'   \item{priorType='mono'}{builds a covariance matrix from the correlation 
#'     matrix of the \code{fitMonoExp} model and a relative uncertainty 
#'     parameter \code{ru_theta}}
#'   \item{priorType='abc'}{estimates the parameters uncertainties by 
#'     a moments matching strategy, and assume no correlation}
#' } 
#' 
#' @export

estimateExpPrior <- function(x, uy, dataType, priorType = 'mono',
                             out, ru_theta = 0.05, eps = 1e-3,    
                             nb_chains = 4, nb_warmup = 800,
                             nb_iter = nb_warmup + 200) {

  statsObs = function(x){
    # Stats for residuals are Q95 and 0.
    c(quantile(abs(x),probs=c(0.95)), 0.)
  }
  
  theta0    = out$best.theta   # Used by next stage
  
  if(priorType == 'mono') {
    # Scale Monoexp covariance matrix by ru_theta
    cor_theta = out$cor.theta
    u_theta   = ru_theta * theta0
    rList = NULL
    
  } else {
    # ABC nocorr approx.
    resid = out$fit$par$resid
    Sobs  = statsObs(resid/uy)
    
    stanData = list(
      N        = length(x), 
      x        = x, 
      uy       = uy, 
      dataType = dataType,
      Sobs     = Sobs,
      Np       = 3,
      theta    = theta0,
      eps      = eps
      )
    init = list(
      u_theta = ru_theta * theta0
    )
    
    # Sample
    fit = rstan::sampling(
      stanmodels$modUQExp,
      data      = stanData,
      init      = function() {init},
      control   = list(adapt_delta=0.995),
      iter      = nb_iter,
      chains    = nb_chains,
      warmup    = nb_warmup,
      verbose   = FALSE
    )
    
    lp   = rstan::extract(fit,'lp__')[[1]]
    map  = which.max(lp)
    u_theta   = rstan::extract(fit,'u_theta')[[1]][map,]
    cor_theta = diag(c(1,1,1)) # Hyp. nocorr
    rList =     list(
      Sobs   = Sobs,
      Ssim   = rstan::extract(fit,'Ssim')[[1]][map,]
    )
    
  }
  
  # Cov. matrix
  Sigma0 = diag(u_theta) %*% cor_theta %*% diag(u_theta)
  
  return(
    list(
      theta0 = theta0,
      Sigma0 = Sigma0,
      rList  = rList
    )
  )
}
