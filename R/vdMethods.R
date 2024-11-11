mean_vd <- function(vd_obj) {
  col_means <- t(as.matrix(colMeans(vd_obj$data)))
  return(Vd(data = col_means))
}

sd_vd <- function(vd_obj) {
  col_sd <- t(as.matrix(apply(vd_obj$data, 2, sd)))
  return(Vd(data = col_sd))
}

#' @title Compute the norm of an object of class `vd` 
#'
#' @param vd_obj An object of class `vd` 
#' @return The norm vector of the an object of class `vd` 
#' @export
norm_vd <- function(vd_obj) {
  return(sqrt(sum(vd_obj$data^2)))
}

