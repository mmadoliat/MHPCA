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
nfd <- function(data,spars_par = NULL) {
  if (!is.matrix(data) && !is.array(data)) {
    stop("Input must be a matrix or an array.")
  }
  if (!is.null(spars_par)) {
    if (!is.numeric(spars_par) || any(spars_par != as.integer(spars_par))) {
      stop("`spars_par` must be an integer vector.", call. = FALSE)
    }
    spars_par <- as.integer(spars_par)
    
    ncol_data <- if (is.matrix(data)) ncol(data) else dim(data)[2]
    if (length(spars_par) >= ncol_data) {
      stop(
        sprintf("length(spars_par) = %d must be < number of columns (%d).",
                length(spars_par), ncol_data),
        call. = FALSE
      )
    }
  }
  class(data) <- c("nfd", class(data))
  attr(data, "spars_par") <- spars_par
  return(data)
}





