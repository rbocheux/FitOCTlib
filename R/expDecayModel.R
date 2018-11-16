#  expDecayModel.R
#' Reference OCT decay model.
#' @param x a numeric vector of depths
#' @param c a numeric vector of parameters
#' @param dataType an numeric (1 or 2) defining the type of data 
#' @return A numeric vector of model values 
#' @author Pascal PERNOT
#' @export

expDecayModel   <- function(x,c,dataType = 2) {
  return( c[1] + c[2] * exp(-dataType*x / c[3]) )
}