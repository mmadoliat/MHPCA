####################not completed yet

#' @title Define a Set of Multivariate Multidimensional Functional Data objects
#'
#' @description
#' The `hd` class represents functional data ...

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
#' hd1 <- Hd(mfd1, mfd2)
#' hd1[1]
#' hd1[1, 1]
#' hd1[1:5, 2]
#' hd1[, 1]
#' hd1[1:5, ]
#' evalarg <- list(argval1, argval2)
#' hd1$eval(evalarg)
#' hd1 + hd1
#' mean(hd1)
#' inprod_hd(hd1, hd1)
#' norm_hd(hd1)
#' plot(hd1)
#' bimfdplot(hd1)
#'
#' @import R6
#' @importFrom fda is.basis eval.basis Data2fd
#'
#' @export
hd <- R6::R6Class("hd",
  public = list(
    #' @description
    #' Constructor for `hd` objects (same as 'Hd')
    #'
    #' @param ... A `mfd` objects which have separated by comma
    #'
    initialize = function(...) {
      hd_list <- list(...)
      if (is.list(hd_list[[1]])) hd_list <- hd_list[[1]]
      init_hd_list_check(hd_list)
      mfd_list <- list()
      vd_list <- list()
      for (i in 1:length(hd_list)) {
        if(inherits(hd_list[[i]],"mfd")) {
          mfd_list[[length(mfd_list)+1]] <- hd_list[[i]]
        } else if(inherits(hd_list[[i]],"vd")) {
          vd_list[[length(vd_list)+1]] <- hd_list[[i]]
        }
      }
      #init_mfd_list_check(mfd_list)
      #init_vd_list_check(vd_list)
      basis_list <- list()
      private$.nobs <- hd_list[[1]]$nobs
      private$.nvar <- length(mfd_list)
      if (length(vd_list) > 0){
        private$.ncols <- sapply(vd_list,function(x) x$ncols)
      }
      if (private$.nvar > 0) {
        for (i in 1:private$.nvar) {
          mfd_list[[i]] <- mfd_list[[i]]$clone()
          basis_list[[i]] <- mfd_list[[i]]$basis
          private$.coefs[[i]] <- mfd_list[[i]]$coefs
        }
        private$.basis <- mvbasismfd$new(basis_list)
      }
      if (length(vd_list)) {
        private$.vdata <- Vd(do.call("cbind",lapply(vd_list,function(x) x$data)))
      }
    },

    #' @description
    #' Eval method for `hd` objects
    #'
    #' @param evalarg A list of numeric vectors of argument values at which the `hd` is to be evaluated.
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
    #' Print method for `hd` objects
    #'
    #' @param ... Additional arguments to be passed to `print`
    #'
    print = function(...) {
      cat("A 'hd' object with", private$.nvar, "variable(s):\n")
      for (i in 1:private$.nvar) {
        cat("\nVariable ", i, ":\n", sep = "")
        print(self[, i])
      }
      cat("\nVector Data", ":\n", sep = "")
      print(self$vdata)
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
    ncols = function(value){
      if (missing(value)) {
        private$.ncols
      } else {
        stop("`$ncols` is read only", call. = FALSE)
      }
    }, 
    vdata = function(value){
      if (missing(value)) {
        private$.vdata
      } else {
        stop("`$vdata` is read only", call. = FALSE)
      }
    }
    
  ),
  private = list(
    .basis = NULL,
    .coefs = NULL, # we record vectorized of the coefs
    .nobs = NULL,
    .nvar = NULL,
    .vdata = NULL,
    .ncols = NULL
  )
)
#' @rdname hd
#' @seealso \code{\link{mvbasismfd}}, \code{\link{mfd}}

#' @title A Class of Multivariate Multidimensional Functional Data objects
#'
#' @description
#' Constructor for `hd` objects (same as `Hd`)
#'
#' @param ... A `mfd` objects which have separated by comma
#' @export
Hd <- function(...) hd$new(...)
