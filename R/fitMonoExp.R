#  fitMonoExp.R
#' Monoexponential fit of OCT decay
#' @param x a numeric vector
#' @param y a numeric vector of responses
#' @param uy a numeric vector of uncertainty on 'y'
#' @param method choice of optimization method in c('optim','sample')
#' @return A list containing
#' \describe{
#'   \item{best.theta}{a vector of optimal parameters}
#'   \item{cor.theta}{correlation matrix of optimal parameters}
#'   \item{method}{choice of optimization method}
#'   \item{fit}{a \code{stanfit} object containg the results of the fit}
#' }
#' @author Pascal PERNOT
#' @details Bayesian inference of the parameters of an exponential
#' model assuming an uncorrelated normal noise
#' \code{y(x) ~ normal(m(x),uy(x));
#' m(x) = theta[1] + theta[2]*exp(-2*x/theta[3])}.
#' @export

fitMonoExp <- function(x, y, uy, method = 'optim',
                       nb_chains = 4, nb_warmup = 500,
                       nb_iter = nb_warmup + 500) {

  stanData = list(
    N =length(x),
    x=x, y=y, uy=uy,
    Np=3
  )

  init = function() {
    list(theta  = c(min(y),max(y)-min(y),mean(x)))
  }

  if(method == 'sample') {
    parOpt = c('theta')
    pars   = c(parOpt,'resid','m','br')

    fit = sampling(stanmodels$modFitExp,
                   data = stanData,
                   pars = pars,
                   init = init,
                   control = list(adapt_delta=0.99, max_treedepth=12),
                   iter = nb_iter, chains = nb_chains,
                   warmup = nb_warmup, verbose=FALSE)

    # Estimate decay params
    theta   = extract(fit,'theta')[[1]]
    theta0  = colMeans(theta)
    thetaCor= cor(theta)

  } else {

    fit = optimizing(stanmodels$modFitExp,
                     data = stanData,
                     init = init,
                     as_vector = FALSE,
                     verbose   = FALSE,
                     hessian   = TRUE )

    # Estimate decay params
    theta0  = fit$par$theta
    thetaCor= cov2cor(solve(-fit$hessian))

  }
  return(
    list(fit        = fit,
         best.theta = theta0,
         cor.theta  = thetaCor,
         method     = method
    )
  )
}
