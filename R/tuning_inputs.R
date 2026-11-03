# Reshaping of user-supplied tuning parameters (smoothing / sparsity) into the
# per-variable lists the algorithms expect. Used by mhpca$initialize() and
# remfpca$initialize(); a new tuning input needs one call here instead of a
# copy of the reshaping logic.

# Reshape `x` into a list with one element per variable, named var1, var2, ...
#
# x            NULL, a vector, a matrix (if allow_matrix) or a list.
# arg          argument name, used in warnings.
# nvar         number of variables the result must cover.
# ncomp        number of components.
# grid         TRUE when `x` is a grid of candidates to search over (GCV / CV),
#              FALSE when it gives the values to use directly.
# mode         "per_component": each variable gets one value per component
#              (sequential power). "per_variable": each variable gets a single
#              value (joint power / eigen).
# vec_label, list_label
#              how the variable count is named in warnings about vector and
#              list inputs (e.g. "'hd_obj$mf$nvar'", "number of variables").
# allow_matrix whether a matrix input is reshaped; otherwise it is passed on
#              unchanged apart from the names.
#
# Inputs of the wrong shape are recycled with a warning.
normalize_tuning <- function(x, arg, nvar, ncomp, grid,
                             mode = c("per_component", "per_variable"),
                             vec_label, list_label, allow_matrix = FALSE) {
  mode <- match.arg(mode)
  per_var <- if (mode == "per_component") ncomp else 1
  warn <- function(msg) warning(msg, call. = FALSE)
  rows_to_list <- function(m) {
    if (mode == "per_component") split(m, row(m)) else as.list(m)
  }

  if (is.vector(x) & !is.list(x)) {
    if (grid == FALSE) {
      if (mode == "per_component") {
        if (length(x) != ncomp) {
          warn(sprintf("The length of '%s' did not match 'ncomp' and has been adjusted accordingly.", arg))
          x <- rep(x, length.out = ncomp)
        }
        x <- replicate(nvar, x, simplify = FALSE)
      } else {
        if (length(x) != nvar) {
          warn(sprintf("The length of '%s' did not match %s and has been adjusted accordingly.", arg, vec_label))
          x <- rep(x, length.out = nvar)
        }
        x <- lapply(1:nvar, function(i) x[i])
      }
    } else {
      # One grid for all variables: it is reused for each variable, which is
      # only worth a warning when there is more than one.
      if (nvar > 1) {
        warn(sprintf("The length of '%s' did not match %s and has been adjusted accordingly.", arg, vec_label))
      }
      x <- replicate(nvar, x, simplify = FALSE)
    }
  } else if (allow_matrix && is.matrix(x)) {
    dims_msg <- sprintf("The dimensions of '%s' did not match the expected size and have been adjusted accordingly.", arg)
    if (grid == FALSE) {
      if (!all(dim(x) == c(nvar, per_var))) {
        x <- x[rep(1:nrow(x), length.out = nvar), rep(1:ncol(x), length.out = per_var), drop = FALSE]
        x <- rows_to_list(x)
        warn(dims_msg)
      } else {
        x <- rows_to_list(x)
      }
    } else {
      if (dim(x)[1] != nvar) {
        x <- x[rep(1:nrow(x), length.out = nvar), , drop = FALSE][1:nvar, , drop = FALSE]
        x <- split(x, row(x))
        warn(dims_msg)
      } else {
        x <- split(x, row(x))
      }
    }
  } else if (is.list(x)) {
    if (grid == FALSE) {
      if (length(x) != nvar) {
        warn(sprintf("Adjusting '%s' to match %s.", arg, list_label))
        x <- rep(x, length.out = nvar)
      }
      x <- lapply(x, function(vec) {
        if (length(vec) != per_var) {
          warn(sprintf("Adjusting vector length in '%s' to match 'ncomp'.", arg))
          vec <- rep(vec, length.out = per_var)
        }
        vec
      })
    } else {
      if (length(x) != nvar) {
        warn(sprintf("Adjusting '%s' to match %s.", arg, list_label))
        x <- rep(x, length.out = nvar)
      }
    }
  }

  if (!is.null(x)) {
    names(x) <- paste0("var", 1:nvar)
  }
  x
}

# Sparsity on the scores (u) is a single vector over components: recycle it to
# ncomp when the values are given directly.
normalize_tuning_u <- function(x, arg, ncomp, grid) {
  if (grid == FALSE & length(x) != ncomp & !is.null(x)) {
    warning(sprintf("The length of '%s' did not match 'ncomp' and has been adjusted accordingly.", arg), call. = FALSE)
    x <- rep(x, length.out = ncomp)
  }
  x
}
