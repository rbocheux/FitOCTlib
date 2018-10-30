#  fitExpGP.R
#' Decay fit with modulation of mean depth
#' @param x a numeric vector
#' @param y a numeric vector of responses
#' @param uy a numeric vector of uncertainty on y
#' @param method choice of optimization method in c('optim','sample')
#' @return A list
#' @author Pascal PERNOT
#' @details Function which proceeds in two steps:
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
                     nb_chains = 4,
                     nb_warmup = 500,
                     nb_iter   = 1000) {

  # Grid of GP control points
  dx  = diff(range(x))/(Nn+1)
  xGP = seq(min(x)+dx/2,max(x)-dx/2,length.out = Nn)

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
    alpha_scale = 0.1,
    rho_scale   = 0.1,
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

  # model = stanmodels$fitExpGP
  # if(lasso) model = stanmodels$fitExpGPLasso

  # Run Stan
  if(method == 'sample') {

    fit = sampling(stanmodels$modFitExpGP,
                   data = stanData,
                   pars = pars,
                   init = function() {init},
                   control = list(adapt_delta=0.99,max_treedepth=12),
                   iter = nb_iter, chains = nb_chains,
                   warmup = nb_warmup, verbose=FALSE)

  } else {

    fit = optimizing(stanmodels$modFitExpGP,
                     data = stanData,
                     init = init,
                     as_vector = FALSE,
                     verbose   = FALSE,
                     algorithm = 'LBFGS',
                     hessian = TRUE,
                     draws = 500,
                     iter = iter,
                     refresh = 500)
    if(is.null(fit$theta_tilde))
      # No Hessian-based sampling done => try to refine solution
      fit = optimizing(stanmodels$modFitExpGP,
                       data = stanData,
                       init = fit$par, # Restart
                       as_vector = FALSE,
                       verbose   = FALSE,
                       algorithm = 'LBFGS',
                       hessian = TRUE,
                       draws = 500,
                       iter = iter,
                       refresh = 500)
  }

  return(list(fit = fit, method = method, xGP = xGP,
              prior_PD = prior_PD, lasso = lasso))
}
