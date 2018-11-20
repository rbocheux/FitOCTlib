#  selX.R
#' Subset OCT signal.

#' @param x a numeric vector
#' @param y a numeric vector of responses
#' @param depthSel range of \code{x} values: \code{c(xmin,xmax)}
#' @param subsample a numeric factor for regular subsampling

#' @return Returns a list containing new vectors \code{x} and \code{y}.

#' @author Pascal PERNOT

#' @details The default parameters retuen the input vectors unchanged.

#' @export


selX <- function(x,y,depthSel=NULL,subSample=1) {
  # Apply selectors to inputs
  
  if(!is.null(depthSel))
    xSel = which(x >= depthSel[1] &
                 x <= depthSel[2]  )
  else
    xSel = 1:length(x)
  
  x = x[xSel]; y = y[xSel]
  
  if(subSample != 1) {
    xSel = seq(1,length(x),by=subSample)
    x = x[xSel]; y = y[xSel]
  }
  return(list(x=x,y=y))
}