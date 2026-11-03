# Functional part of a data or component object as an mvmfd (NULL if none).
functional_part <- function(x) {
  if (is.hd(x)) x$mf else if (is.mvmfd(x)) x else if (is.mfd(x)) Mvmfd(x) else NULL
}

# Plot each functional principal component as the mean +/- a multiple of the
# component, one panel per variable. Used by plot() for mhpca and remfpca fits.
plot_components <- function(pc_mfd, mean_mfd, values, comp_index = NULL, var_index = NULL, ask = TRUE, expand = NULL, nx = 100, xlab = NULL, ylab = NULL, ...) {
  if (is.null(pc_mfd)) stop("This fit has no functional part to plot.", call. = FALSE)
  percentvar <- round(100 * (values) / sum(values), 1)
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
  # expand is indexed by component number below, so it covers all components
  if (is.null(expand)) expand <- 2 * sqrt(values)
  if (length(expand) == 1) expand <- rep(expand, length(values))
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

plot_mhpca <- function(mhpca_obj, ...) {
  if (!(inherits(mhpca_obj, "mhpca"))) stop("Argument 'mhpca_obj' is not a mhpca object.")
  plot_components(functional_part(mhpca_obj$pc_hd), functional_part(mhpca_obj$mean_hd),
                  mhpca_obj$values, ...)
}



#' Plot cross-validation trajectories for MHPCA
#'
#' Generates plots of cross-validation (CV) error scores against the tuning parameter \code{gamma}
#' for an \code{mhpca} R6 object.  Supports three modes:
#' \itemize{
#'   \item \code{"nfd"}: non-functional data (multiple variables per PC)
#'   \item \code{"fd"}: functional data (multiple variables per PC)
#'   \item \code{"u"}: univariate mode (one matrix per PC)
#' }
#'
#' @param mhpca_obj An R6 object of class \code{mhpca} containing slots:
#'   \code{CVs_nfd}, \code{CVs_fd}, or \code{CVs_u} and matching
#'   \code{sparse_tuning_nfd}, \code{sparse_tuning_fd}, or \code{sparse_tuning_u}.
#' @param type Character; one of \code{"nfd"} (default), \code{"fd"}, or \code{"u"}.
#'   Indicates which CV-list and tuning-list to use for plotting.
#' @param PC Integer or \code{NULL}; index of the principal component to plot.
#'   If \code{NULL} (default), plots all available components.
#' @param Var Integer or \code{NULL}; index of the variable (within a component) to plot.
#'   Ignored when \code{type = "u"}.  If \code{NULL} (default), plots all variables for each PC.
#'
#' @return Invisibly returns \code{NULL}.  The function is called for its side effect of
#'   producing base-R plots with CV trajectories, error bars, and selected-\code{gamma} lines.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Plot all non-functional CVs interactively:
#' plot_cv(res_hpca)
#'
#' # Plot only the 2nd variable of the 3rd principal component for functional data:
#' plot_cv(res_hpca, type = "fd", PC = 3, Var = 2)
#' }
plot_cv <- function(mhpca_obj,
                    type = c("nfd", "fd", "u"),
                    PC = NULL,
                    Var = NULL) {
  type <- match.arg(type)
  CVs <- mhpca_obj[[paste0("CVs_", type)]]
  tunes <- mhpca_obj[[paste0("sparse_tuning_", type)]]
  pc_names <- names(CVs)

  if (length(CVs) == 0) {
    if (identical(mhpca_obj$selection_method, "ebic")) {
      stop("plot_cv() plots cross-validation scores, but this fit selected sparsity by eBIC (see its eBICs_* fields).", call. = FALSE)
    }
    stop(sprintf("This fit has no cross-validation scores for type = \"%s\".", type), call. = FALSE)
  }
  # --- validate & convert PC index ---
  if (!is.null(PC)) {
    if (!is.numeric(PC) || length(PC) != 1 || PC != as.integer(PC)) {
      stop("`PC` must be a single integer index.")
    }
    if (PC < 1 || PC > length(CVs)) {
      stop("`PC` index out of range [1:", length(CVs), "].")
    }
    PC_name <- pc_names[PC]
  }

  # --- validate & convert Var index (only for type != "u") ---
  if (!is.null(Var)) {
    if (type == "u") {
      warning("`Var` is ignored when type = 'u'.")
      Var_name <- NULL
    } else {
      if (!is.numeric(Var) || length(Var) != 1 || Var != as.integer(Var)) {
        stop("`Var` must be a single integer index.")
      }
      # must specify PC when Var supplied
      if (is.null(PC)) {
        stop("You must supply `PC` if you want to plot a specific `Var`.")
      }
      mat_list <- CVs[[PC_name]]
      if (Var < 1 || Var > length(mat_list)) {
        stop("`Var` index out of range for ", PC_name, " [1:", length(mat_list), "].")
      }
      Var_name <- names(mat_list)[Var]
    }
  }

  # --- single-plot shortcut? ---
  single_plot <- !is.null(PC) && (type == "u" || !is.null(Var))
  if (single_plot) {
    # extract the one matrix + selected gam
    if (type == "u") {
      m <- CVs[[PC_name]]
      sel_val <- tunes[[PC_name]]
      ttl <- paste0(PC_name, " (u)")
    } else {
      m <- CVs[[PC_name]][[Var_name]]
      sel_vec <- tunes[[PC_name]]
      sel_val <- sel_vec[Var]
      desc <- if (type == "fd") "functional variable" else "non-functional variable"
      ttl <- paste0(ordinal(Var), " ", desc, " - ", PC_name)
    }

    # unpack data
    gam <- m[, "gamma"]
    cv <- m[, "CV"]
    se <- m[, "SE"]

    # y-limits to fit error bars
    y_lo <- min(cv - se)
    y_hi <- max(cv + se)
    pad <- 0.05 * (y_hi - y_lo)
    ylim <- c(y_lo - pad, y_hi + pad)

    # draw
    plot(gam, cv,
      type = "b", pch = 19, lwd = 1,
      xlab = expression(gamma), ylab = "CV",
      main = ttl, ylim = ylim
    )
    arrows(gam, cv - se, gam, cv + se, code = 3, angle = 90, length = 0.05)
    if (!sel_val %in% gam && sel_val %in% seq_along(gam)) {
      sel_val <- gam[sel_val]
    }
    abline(v = sel_val, lty = 2)

    return(invisible(NULL))
  }

  # --- otherwise: full interactive loop ---
  total_plots <- if (type == "u") {
    length(CVs)
  } else {
    sum(vapply(CVs, length, integer(1)))
  }
  plot_idx <- 0

  for (pc in pc_names) {
    if (type == "u") {
      mat_list <- list(Var1 = CVs[[pc]])
      sel_vec <- tunes[[pc]]
    } else {
      mat_list <- CVs[[pc]]
      sel_vec <- tunes[[pc]]
    }

    for (var in names(mat_list)) {
      plot_idx <- plot_idx + 1
      m <- mat_list[[var]]
      gam <- m[, "gamma"]
      cv <- m[, "CV"]
      se <- m[, "SE"]

      # dynamic ylim
      y_lo <- min(cv - se)
      y_hi <- max(cv + se)
      pad <- 0.05 * (y_hi - y_lo)
      ylim <- c(y_lo - pad, y_hi + pad)

      # title
      if (type == "u") {
        ttl <- paste0(pc, " (u)")
      } else {
        vnum <- as.integer(sub("Var", "", var))
        desc <- if (type == "fd") "functional variable" else "non-functional variable"
        ttl <- paste0(ordinal(vnum), " ", desc, " - ", pc)
      }

      # plot + bars + line
      plot(gam, cv,
        type = "b", pch = 19, lwd = 1,
        xlab = expression(gamma), ylab = "CV",
        main = ttl, ylim = ylim
      )
      arrows(gam, cv - se, gam, cv + se, code = 3, angle = 90, length = 0.05)
      sel_val <- if (type == "u") tunes[[pc]] else sel_vec[as.integer(sub("Var", "", var))]
      if (!sel_val %in% gam && sel_val %in% seq_along(gam)) {
        sel_val <- gam[sel_val]
      }
      abline(v = sel_val, lty = 2)

      # pause
      if (plot_idx < total_plots) {
        invisible(readline(prompt = "Hit <Return> to see next plot: "))
      }
    }
  }
}

# Method for the internal generic behind plot()
#' @exportS3Method NULL
op_plot.mhpca <- function(x, ...) plot_mhpca(x, ...)
