#  fitMonoExp.R
#' Monoexponential decay fit
#' @param x a numeric vector
#' @param y a numeric vector of responses
#' @param uy a numeric vector of uncertainty on y
#' @param method choice of optimization method in c('optim','sample')
#' @return A list
#' @author Pascal PERNOT
#' @details Function which proceeds in two steps:
#' @details (1) get a set of residuals R using a smoothing splines model
#' @details (2) estimate the x-dependent standard deviation of the residuals
#' @details by bayesian inference: R(x) ~ normal(0,uy(x));
#' @details uy(x) = theta[1]*exp(-x/theta[2]) assuming a Poisson-type noise.
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
