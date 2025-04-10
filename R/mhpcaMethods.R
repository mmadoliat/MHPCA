plot_mhpca <- function(mhpca_obj, comp_index = NULL, var_index = NULL, ask = TRUE, expand = NULL, nx = 100, xlab = NULL, ylab = NULL, ...) {
  
  if (!(inherits(mhpca_obj, "mhpca"))) stop("Argument 'mhpca_obj' is not a mhpca object.")
  percentvar <- round(100 * (mhpca_obj$values) / sum(mhpca_obj$values), 1)
  mean_mfd <- mhpca_obj$mean_hd$mf
  pc_mfd <- mhpca_obj$pc_hd$mf
  n <- pc_mfd$nobs
  p <- pc_mfd$nvar
  if (is.null(var_index)) var_index <- 1:p
  flag <- FALSE
  if (is.null(comp_index)) {
    comp_index <- 1:n
    flag <- TRUE
  }
  if (is.null(ylab)) ylab <- paste("Variable", var_index)
  if (is.null(xlab)) xlab <- rep("time", length(var_index))
  if (is.null(expand)) expand <- 2 * sqrt(mhpca_obj$values[comp_index])
  old <- par()
  exclude_pars <- c("cin", "cra", "csi", "cxy", "din", "page")
  ind <- which(!(names(old) %in% exclude_pars))
  on.exit(par(old[ind]))
  if (flag) {
    par(mfrow = c(length(var_index), 1), ask = ask)
    for (ipc in comp_index) {
      for (j in var_index) {
        dimSupp <- pc_mfd[ipc, j]$basis$dimSupp
        supp <- pc_mfd[ipc, j]$basis$supp
        x_grids <- seq(supp[1, 1], supp[2, 1], len = nx)
        width <- expand[ipc] * pc_mfd[ipc, j]$eval(x_grids)
        mu <- mean_mfd[1, j]$eval(x_grids)
        ylim <- range(mu - width, mu + width)
        plot(x_grids, mu,
             type = "l", ylim = ylim, ylab = ylab[j], xlab = xlab[j],
             main = paste("FPC", ipc, "(", percentvar[ipc], "%)"), ...
        )
        points(x_grids, mu - width, pch = "-", col = 2, ...)
        points(x_grids, mu + width, pch = "+", col = 3, ...)
      }
    }
  } else {
    par(mfrow = c(length(var_index), length(comp_index)))
    for (j in var_index) {
      for (ipc in comp_index) {
        dimSupp <- pc_mfd[ipc, j]$basis$dimSupp
        supp <- pc_mfd[ipc, j]$basis$supp
        x_grids <- seq(supp[1, 1], supp[2, 1], len = nx)
        width <- expand[ipc] * pc_mfd[ipc, j]$eval(x_grids)
        mu <- mean_mfd[1, j]$eval(x_grids)
        ylim <- range(mu - width, mu + width)
        plot(x_grids, mu,
             type = "l", ylim = ylim, ylab = ylab[j], xlab = xlab[j],
             main = paste("FPC", ipc, "(", percentvar[ipc], "%)"), ...
        )
        points(x_grids, mu - width, pch = "-", col = 2, ...)
        points(x_grids, mu + width, pch = "+", col = 3, ...)
      }
    }
  }
}


