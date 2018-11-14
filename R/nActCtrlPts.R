#  nActCtrlPts.R
#' Count active control points in fitExpGP output.

#' @param fit an object (list or stanfit) issued from a stan code
#' @param prob a probability threshold

#' @return The number of active points or NULL. 

#' @author Pascal PERNOT

#' @details Parameter \code{prob} defines the probability
#' above which a control point is considered active by differing 
#' from zero.

#' @export

nActCtrlPts     <- function(fit, prob=0.90){
  # Count active/non-zero ctrl points at prob% level
  # for modFitExpGP model
  
  ## Get a sample of ctrl points or return NULL
  if(class(fit) == 'stanfit') {
    yGP     = rstan::extract(fit,'yGP')[[1]]
  } else {
    if(!is.null(fit$theta_tilde)) {
      S = fit$theta_tilde
      c = which(grepl(pattern = 'yGP\\[', x=colnames(S)))
      yGP = S[,c]
    } else {
      return(NULL)
    }
  }
  ## Estimate p% CI on yGP
  Q = t(
    apply(
      yGP,
      2,
      function(x)
        quantile(x,probs = 0.5 + 0.5*prob*c(-1,1))
    )
  )
  ## Count zeroes out of p% CI
  return( sum( apply(Q,1,prod) > 0 ) )
}