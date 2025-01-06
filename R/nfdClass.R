#' Define a Set of Multidimensional Non Functional Data Objects
#' 
#' @description
#' The `nfd` class represents a set of multidimensional vector data.
#' Vector data objects are constructed by specifying a matrix of data points.
#' 
#' @param data A numeric matrix of data points.
#'
#' @return 
#' @export
#'
#' @examples
nfd <- function(data) {
  if (!is.matrix(data) && !is.array(data)) {
    stop("Input must be a matrix or an array.")
  }
  class(data) <- c("nfd", class(data))
  return(data)
}





