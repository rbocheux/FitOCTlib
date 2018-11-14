#  fitExpGP.R
#' Decay fit with modulation of mean depth by Gaussian Process
#' @param x a numeric vector
#' @param y a numeric vector of responses
#' @param uy a numeric vector of uncertainty on 'y'
#' @param method choice of optimization method in c('optim','sample')
#' @param Nn number of control points
#' @param theta0    theta prior mean vbalues
#' @param cor_theta theta prior correlation matrix
#' @param ru_theta  theta prior uncertainty scale
#' @param lambda_rate scale of ctrl points prior
#' @param alpha_scale SD scale of GP
#' @param rho_scale relative correlation length of GP
#' @param gridType type of controle points grid (`internal` does not contain boundaries)
#' @param lasso flag to use lasso prior
#' @param iter max. number of iterations for 'optim'
#' @param prior_PD  flag to sample from prior pdf only
#' @param nb_chains  number of MCMC chains
#' @param nb_warmup number of warmup steps
#' @param nb_iter   number of steps
#' @return A list containing
#' \describe{
#'   \item{fit}{a \code{stanfit} object containg the results of the fit}
#'   \item{xGP}{a vector of coordinates for the control points}
#'   \item{method}{same as input}
#'   \item{prior_PD}{same as input}
#'   \item{lasso}{same as input}
#' }
#' @author Pascal PERNOT
#' @details Bayesian inference of the parameters of an exponential
#' model with modulation assuming an uncorrelated normal noise
#' \code{y(x) ~ normal(m(x),uy(x));
#' m(x) = theta[1] + theta[2]*exp(-2*x/theta[3]*(1+l(x)));}.
#' \code{l(x)} is defined by a GP with fixed positions 'xGP',
#' variance and correlation length.
#' @export

fitExpGP <- function(x, y, uy,
                     Nn        = 10,       # Nb of control points
                     theta0    = NULL,     # Theta prior
                     cor_theta = NULL,     # ...
                     ru_theta  = 0.1,      # ...
                     lambda_rate = 0.1,    # Scale of ctrl points prior
                     lasso     = FALSE,    # Use lasso prior ?
                     method    = 'sample', # One of 'sample','optim'
                     iter      = 50000,    # Max. iterations of optimizer
                     prior_PD  = 0,        # Flag to sample from prior only
                     alpha_scale = 0.1,    # SD scale of GP
                     rho_scale = 1/Nn,     # Relative correlation length of GP
                     gridType  = 'internal',
                     nb_chains = 4,
                     nb_warmup = 500,
                     nb_iter   = 1000,
                     verbose   = FALSE,
                     open_progress = TRUE) {

  # Grid of GP control points
  if(gridType == 'internal') {
    dx  = diff(range(x))/(Nn+1)
    xGP = seq(min(x)+dx/2,max(x)-dx/2,length.out = Nn)
  } else {
    xGP = seq(min(x),max(x),length.out = Nn)
  }

  # Initial monoexp params
  if(is.null(theta0))
    theta0 = c(min(y),max(y)-min(y),mean(x))

  if(is.null(cor_theta))
    cor_theta = diag(1,length(theta0))

  stanData = list(
    N =length(x),
    x=x, y=y, uy=uy,
    Np = 3,
    Nn = Nn,
    xGP = xGP,
    alpha_scale = alpha_scale,
    rho_scale   = rho_scale,
    theta0      = theta0,
    cor_theta   = cor_theta,
    ru_theta    = ru_theta,
    prior_PD    = prior_PD,
    lambda_rate = lambda_rate
  )

  init = list(
    theta  = theta0,
    yGP    = 0.01*rnorm(Nn),
    sigma  = 1.0
  )

  # Parameters to scatterplot
  parP = c('theta','sigma')
  if(!lasso) parP = c(parP,'lambda')

  # Parameters to report
  parOpt = c(parP,'yGP')

  # Parameters to save for plots
  pars   = parOpt
  if(prior_PD == 0)
    pars   = c(pars,'resid','br','m','dL')

  # Run Stan
  if(method == 'sample') {

    fit = rstan::sampling(
      stanmodels$modFitExpGP,
      data = stanData,
      pars = pars,
      init = function() {init},
      control = list(adapt_delta=0.99,max_treedepth=12),
      iter = nb_iter,
      chains = nb_chains,
      warmup = nb_warmup,
      verbose = verbose,
      open_progress = open_progress
    )

  } else {

    fit = rstan::optimizing(
      stanmodels$modFitExpGP,
      data = stanData,
      init = init,
      as_vector = FALSE,
      algorithm = 'LBFGS',
      hessian = TRUE,
      draws = nb_iter,
      iter = iter,
      refresh = 500,
      verbose = verbose
    )
    
    if(is.null(fit$theta_tilde))
      # No Hessian-based sampling done => try once to refine solution
      fit = rstan::optimizing(
        stanmodels$modFitExpGP,
        data = stanData,
        init = fit$par, # Restart
        as_vector = FALSE,
        algorithm = 'LBFGS',
        hessian = TRUE,
        draws = nb_iter,
        iter = iter,
        refresh = 500,
        verbose = verbose
      )
  }

  return(list(fit = fit, method = method, xGP = xGP,
              prior_PD = prior_PD, lasso = lasso))
}
