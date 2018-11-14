#  expDecayModel.R
#' Reference OCT decay model.
#' @param x a numeric vector of depths
#' @param c a numeric vector of parameters
#' @return A numeric vector of model values 
#' @author Pascal PERNOT
#' @export

expDecayModel   <- function(x,c) {
  return( c[1] + c[2] * exp(-2*x/c[3]) )
}