#' @title Define a Set of Multivariate Multidimensional Functional Data objects
#'
#' @description
#' The `mvmfd` class represents functional data ...

#' @field basis A `mvbasismfd` object
#' @field coefs a matrix of the coefficients.
#' @field nobs number of observation
#' @field nvar number of variables
#'
#' @examples
#' require(fda)
#' bs1 <- create.fourier.basis(c(0, 2 * pi), 5)
#' bs2 <- create.bspline.basis(c(0, 1), 7)
#' bs3 <- create.exponential.basis(c(0, 2), 3)
#' nobs <- 10
#' argval1 <- seq(0, 2 * pi, length.out = 12)
#' X1 <- outer(sin(argval1), seq(0.5, 1.5, length.out = nobs))
#' mdbs1 <- Basismfd(bs1)
#' mfd1 <- Mfd(argval1, X1, mdbs1)
#' mdbs2 <- Basismfd(bs1)
#' argval2 <- argval1
#' X2 <- outer(cos(argval2), seq(0.2, 1.5, length.out = nobs))
#' mfd2 <- Mfd(argval2, X2, mdbs1)
#' mvmfd1 <- Mvmfd(mfd1, mfd2)
#' mvmfd1[1]
#' mvmfd1[1, 1]
#' mvmfd1[1:5, 2]
#' mvmfd1[, 1]
#' mvmfd1[1:5, ]
#' evalarg <- list(argval1, argval2)
#' mvmfd1$eval(evalarg)
#' mvmfd1 + mvmfd1
#' mean(mvmfd1)
#' inprod_mvmfd(mvmfd1, mvmfd1)
#' norm_mvmfd(mvmfd1)
#' plot(mvmfd1)
#' bimfdplot(mvmfd1)
#'
#' @import R6
#' @importFrom fda is.basis eval.basis Data2fd
#'
#' @export
mvmfd <- R6::R6Class("mvmfd",
  public = list(
    #' @description
    #' Constructor for `mvmfd` objects (same as 'Mvmfd')
    #'
    #' @param ... A `mfd` objects which have separated by comma
    #'
    initialize = function(...) {
      mfd_list <- list(...)
      if (is.list(mfd_list[[1]])) mfd_list <- mfd_list[[1]]
      init_mfd_list_check(mfd_list)
      basis_list <- list()
      private$.nobs <- mfd_list[[1]]$nobs
      private$.nvar <- length(mfd_list)
      for (i in 1:private$.nvar) {
        mfd_list[[i]] <- mfd_list[[i]]$clone()
        basis_list[[i]] <- mfd_list[[i]]$basis
        private$.coefs[[i]] <- mfd_list[[i]]$coefs
      }
      private$.basis <- mvbasismfd$new(basis_list)
      sp_list <- lapply(mfd_list, function(comp) {
        v <- comp$spars_par
        if (length(v) == 0L) NULL else v
      })
      if (all(vapply(sp_list, is.null, logical(1)))) {
        private$.spars_par <- NULL
      } else {
        # replace NULLs with 0
        private$.spars_par <- lapply(sp_list, function(v) if (is.null(v)) 0L else v)
      }
      
      sm_list <- lapply(mfd_list, function(comp) {
        v <- comp$smooth_par
        if (length(v) == 0L) NULL else v
      })
      if (all(vapply(sp_list, is.null, logical(1)))) {
        private$.smooth_par <- NULL
      } else {
        # replace NULLs with 0
        private$.smooth_par <- lapply(sm_list, function(v) if (is.null(v)) 0L else v)
      }
    },

    #' @description
    #' Eval method for `mvmfd` objects
    #'
    #' @param evalarg A list of numeric vectors of argument values at which the `mvmfd` is to be evaluated.
    #' @return A list of evaluated values
    eval = function(evalarg) {
      Bmat <- private$.basis$eval(evalarg)
      Xhat <- list()
      for (i in 1:private$.nvar) {
        if (is.matrix(private$.coefs[[i]])) {
          Xhat[[i]] <- Bmat[[i]][[1]] %*% private$.coefs[[i]]
        } else {
          Xhat[[i]] <- (Bmat[[i]][[2]] %x% Bmat[[i]][[1]]) %*% apply(private$.coefs[[i]], 3, as.vector)
          Xhat[[i]] <- array(Xhat[[i]], dim = c(sapply(evalarg[[i]], length), private$.nobs))
        }
      }
      return(Xhat)
    },
    #' @description
    #' Print method for `mvmfd` objects
    #'
    #' @param ... Additional arguments to be passed to `print`
    #'
    print = function(...) {
      cat("A 'mvmfd' object with", private$.nvar, "variable(s):\n")
      for (i in 1:private$.nvar) {
        cat("\nVariable ", i, ":\n", sep = "")
        print(self[, i])
      }
      invisible(self)
    }
  ),
  active = list(
    # basis field
    basis = function(value) {
      if (missing(value)) {
        private$.basis
      } else {
        stop("`$basis` is read only", call. = FALSE)
      }
    },

    # coefs field
    coefs = function(value) {
      if (missing(value)) {
        private$.coefs
      } else {
        stop("`$coefs` is read only", call. = FALSE)
      }
    },

    # nvar field
    nvar = function(value) {
      if (missing(value)) {
        private$.nvar
      } else {
        stop("`$nvar` is read only", call. = FALSE)
      }
    },

    # nobs field
    nobs = function(value) {
      if (missing(value)) {
        private$.nobs
      } else {
        stop("`$nobs` is read only", call. = FALSE)
      }
    }, 
    spars_par = function(value) {
      if (missing(value)) {
        # getter
        return(private$.spars_par)
      }
      # setter
      if (!is.list(value)) {
        stop("`spars_par` must be a **list**.", call. = FALSE)
      }
      if (length(value) != private$.nvar) {
        stop(sprintf("`spars_par` must be length %d (nvar).", private$.nvar),
             call. = FALSE)
      }
      nbv <- private$.basis$nbasis
      # each element must be integer‐vector
      for (i in seq_along(value)) {
        elt <- value[[i]]
        if (!is.numeric(elt) || any(elt != as.integer(elt))) {
          stop(sprintf("`spars_par[[%d]]` must be an integer vector.", i),
               call. = FALSE)
        }
        value[[i]] <- as.integer(elt)
        if (length(elt) >= nbv[i]) {
          stop(sprintf("length(spars_par[[%d]]) = %d must be < %d (nbasis[%d]).",
                       i, length(vi), nbv[i], i),
               call. = FALSE)
        }
      }
      private$.spars_par <- value
    }, 
    smooth_par = function(value) {
      if (missing(value)) {
        # getter
        return(private$.smooth_par)
      }
      # setter
      if (!is.list(value)) {
        stop("`smooth_par` must be a **list**.", call. = FALSE)
      }
      if (length(value) != private$.nvar) {
        stop(sprintf("`smooth_par` must be length %d (nvar).", private$.nvar),
             call. = FALSE)
      }
      private$.smooth_par <- value
    }
  ),
  private = list(
    .basis = NULL,
    .coefs = list(), # we record vectorized of the coefs
    .nobs = NULL,
    .nvar = NULL,
    .spars_par = NULL,
    .smooth_par = NULL
  )
)
#' @rdname mvmfd
#' @seealso \code{\link{mvbasismfd}}, \code{\link{mfd}}

#' @title A Class of Multivariate Multidimensional Functional Data objects
#'
#' @description
#' Constructor for `mvmfd` objects (same as `Mvmfd`)
#'
#' @param ... A `mfd` objects which have separated by comma
#' @export
Mvmfd <- function(...) mvmfd$new(...)
