#  estimateNoise.R
#' Estimate and model noise in signal
#' @param x a numeric vector
#' @param y a numeric vector of responses
#' @param df smoothing factor for `smooth.splines`
#' @return A list
#' @author Pascal PERNOT
#' @details Function which proceeds in two steps:
#' @details (1) get a set of residuals R using a smoothing splines model
#' @details (2) estimate the x-dependent standard deviation of the residuals
#' @details by bayesian inference: R(x) ~ normal(0,uy(x));
#' @details uy(x) = theta[1]*exp(-x/theta[2]) assuming a Poisson-type noise.
#' @export

estimateNoise <- function(x, y, df = 15) {

  # Smoothing
  ySpl   = smooth.spline(x,y,df=df)$y
  resSpl = y-ySpl

  # Fit smoothing residuals by exponential variance model

  stanData = list(N =length(x), x=x, y=resSpl)
  init = list(theta = c(max(resSpl),mean(x)))

  # Optimize
  fit = rstan::optimizing(stanmodels$modHetero,
                          data = stanData,
                          init = init,
                          as_vector = FALSE,
                          verbose   = FALSE)

  # Estimate data uncertainty
  theta = fit$par$theta
  sig = theta[1]*exp(-x/theta[2])

  return(list(theta = theta, uy = sig, ySmooth = ySpl))
}
