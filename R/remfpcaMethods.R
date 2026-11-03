plot_remfpca <- function(remfpca_obj, ...) {
  if (!(inherits(remfpca_obj, "remfpca"))) stop("Argument 'remfpca_obj' is not a remfpca object.")
  plot_components(functional_part(remfpca_obj$pc_mfd), functional_part(remfpca_obj$mean_mfd),
                  remfpca_obj$values, ...)
}


# Method for the internal generic behind plot()
#' @exportS3Method NULL
op_plot.remfpca <- function(x, ...) plot_remfpca(x, ...)
