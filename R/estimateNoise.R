#  estimateNoise.R
#' Estimate and model noise in signal
#' @param x a numeric vector
#' @param y a numeric vector of responses
#' @param df smoothing factor for \code{smooth.splines}
#' @param maxRate max. value of rate parameter
#' @return A list containing
#' \describe{
#'   \item{fit}{a \code{stanfit} object containg the results of the fit}
#'   \item{method}{choice of optimization method}
#'   \item{theta}{a vector of optimal parameters}
#'   \item{uy}{a vector of estimated uncertainty values for \code{y}}
#'   \item{ySmooth}{ a vector of values for the smoother curve}
#' }
#' @author Pascal PERNOT
#' @details Function which proceeds in two steps:
#' \enumerate{
#'   \item get a set of residuals R using a smoothing splines model
#'   \item estimate the x-dependent standard deviation of the residuals
#' }
#' by bayesian inference: \code{R(x) ~ normal(0,uy(x));
#' uy(x) = theta[1]*exp(-x/theta[2])} assuming a Poisson-type noise.
#' @export

estimateNoise <- function(x, y, df = 15, maxRate = 10000) {

  # Smoothing
  ySpl   = stats::smooth.spline(x,y,df=df)$y
  resSpl = y-ySpl

  # Fit smoothing residuals by exponential variance model

  stanData = list(N =length(x), x=x, y=resSpl, maxRate = maxRate)
  init = list(theta = c(max(resSpl),mean(x)))

  # Optimize
  fit = rstan::optimizing(
    stanmodels$modHetero,
    data = stanData,
    init = init,
    as_vector = FALSE,
    verbose   = FALSE,
    hessian   = TRUE
  )
  
  # Estimate data uncertainty
  theta = fit$par$theta
  sig = theta[1]*exp(-x/theta[2])

  return(list(fit = fit, theta = theta, uy = sig, 
              ySmooth = ySpl, method='optim', data = stanData))
}
