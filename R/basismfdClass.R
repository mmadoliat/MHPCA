#' @title Define a Set of Multidimensional Functional Basis 
#'
#' @description
#' The `basismfd` class represents a set of  multidimensional basis functions. This class utilizes basis objects 
#' from the `fda` package, such as B-splines and Fourier bases.
#' @field basis A list of basis objects from the `fda` package.
#' @field dimSupp The dimension of the support domain of the `basismfd` object.
#' @field supp The matrix representing the ranges of the dimensions.
#' @field gram The Gram matrix.
#' @field nbasis A numeric vector containing the number of bases.
#' @examples
#' require(fda)
#' bs1 <- create.fourier.basis(c(0, 2 * pi), 5)
#' bs2 <- create.bspline.basis(c(0, 1), 7)
#' bs3 <- create.exponential.basis(c(0, 2), 3)
#' # 1-D Basis ######## (similar to the fd features)
#' mdbs1 <- Basismfd(bs1)
#' mdbs1$basis
#' mdbs1$dimSupp
#' mdbs1$nbasis
#' mdbs1$supp
#' mdbs1$gram
#' mdbs1$eval(1:7 / 10)
#' image(as.matrix(mdbs1$gram))
#' 
#' ####### 2-D Basis ######## (fd cannot handle this)
#' mdbs2 <- Basismfd(bs1, bs2)
#' mdbs2$basis
#' mdbs2$dimSupp
#' mdbs2$nbasis
#' mdbs2$supp
#' dim(mdbs2$gram)
#' arg_mdbs <- list(1:10, 1:9 / 10)
#' mdbs2$eval(arg_mdbs)
#' image(as.matrix(mdbs2$gram))

#' @import R6
#' @importFrom fda is.basis eval.basis inprod
#' @importFrom Matrix Matrix
#' @export
basismfd <- R6::R6Class("basismfd",
  public = list(
    #' @description
    #'  The constructor function for objects of the class `basismfd` (same as Basismfd(...) )
    #' @usage Basismfd(...)
    #' @param ... A list of `basisfd` objects
    initialize = function(...) {
      basis <- list(...)
      init_basismfd_check(basis)

      if (is.basis(basis)) {
        private$.gram <- Matrix(inprod(basis, basis))
        private$.basis <- list(basis)
        private$.dimSupp <- 1
        private$.supp <- matrix(basis$rangeval, nrow = 2, ncol = 1)
        private$.nbasis <- basis$nbasis
      } else {
        private$.basis <- basis
        private$.dimSupp <- length(basis)
        private$.supp <- matrix(0, nrow = 2, ncol = length(basis))
        private$.gram <- 1
        for (i in 1:private$.dimSupp) {
          private$.gram <- Matrix(inprod(basis[[i]], basis[[i]]) %x% private$.gram)
          private$.nbasis[i] <- basis[[i]]$nbasis
          private$.supp[, i] <- basis[[i]]$rangeval
        }
      }
    },
    #' @description
    #' Evaluate the `basismfd` object at given argument values
    #' @param evalarg A list of numeric vectors of argument values at which the `basismfd` is to be evaluated
    #' @return A list of evaluated values
    eval = function(evalarg) {
      eval_basismf_validity_check(evalarg, private$.dimSupp)

      if (is.numeric(evalarg)) {
        evalarg <- list(evalarg)
      }

      out <- list()
      for (i in 1:length(evalarg)) {
        out[[i]] <- eval.basis(evalarg[[i]], private$.basis[[i]])
      }
      return(out)
    },
    #' @description
    #' Print method for `basismfd` objects
    #'
    #' @param ... Additional arguments to be passed to `print`
    print = function(...) {
      for (i in 1:private$.dimSupp) {
        cat("basis ", i, ":\n", sep = "")
        cat("type:", private$.basis[[i]]$type)
        cat("\nnbasis:", private$.nbasis[i])
        cat("\nsupport:", private$.supp[, i], "\n")
        if (i != private$.dimSupp) cat("\n")
      }
      invisible(self)
    }
  ),
  active = list(
    #' Getter and setter for `basis` field
    basis = function(value) {
      if (missing(value)) {
        private$.basis
      } else {
        stop("`$basis` is read only", call. = FALSE)
      }
    },

    #' Getter and setter for `dimSupp` field
    dimSupp = function(value) {
      if (missing(value)) {
        private$.dimSupp
      } else {
        stop("`$dimSupp` is read only", call. = FALSE)
      }
    },

    #' Getter and setter for `nbasis` field
    nbasis = function(value) {
      if (missing(value)) {
        private$.nbasis
      } else {
        stop("`$nbasis` is read only", call. = FALSE)
      }
    },

    #' Getter and setter for `supp` field
    supp = function(value) {
      if (missing(value)) {
        private$.supp
      } else {
        stop("`$supp` is read only", call. = FALSE)
      }
    },

    #' Getter and setter for `gram` field
    gram = function(value) {
      if (missing(value)) {
        private$.gram
      } else {
        stop("`$gram` is read only", call. = FALSE)
      }
    }
  ),
  private = list(
    .basis = NULL,
    .dimSupp = NULL,
    .nbasis = NULL,
    .supp = NULL,
    .gram = NULL
  )
)
#' @rdname basismfd
#' 
#' @description
#'  Constructor for `basismfd` objects (same as Basismfd(...) )

#' @title A Class of Multidimensional Basis Functions
#' @param ... A list of `basisfd` objects
#' @export
Basismfd <- function(...) {
  basismfd$new(...)
}



# Operations that are not defined for basis objects: give a clear error instead
# of an obscure one from base R.
basis_op_error <- function(op, x) {
  stop(sprintf("`%s` is not defined for basis objects (%s); use it on the data objects (mfd, mvmfd, hd) instead.",
               op, class(x)[1]), call. = FALSE)
}

#' Operations on basis objects
#'
#' Arithmetic, comparison and subsetting are not defined for basis objects
#' (`basismfd`, and arithmetic for `mvbasismfd`); these methods only give an
#' informative error. Use them on the data objects (`mfd`, `mvmfd`, `hd`).
#' @param e1,e2 Operands.
#' @param x A `basismfd` object.
#' @param ... Ignored.
#' @return None; these methods always signal an error.
#' @export
Ops.basismfd <- function(e1, e2) {
  basis_op_error(.Generic, if (inherits(e1, c("basismfd", "mvbasismfd"))) e1 else e2)
}

#' @rdname Ops.basismfd
#' @export
Ops.mvbasismfd <- function(e1, e2) {
  basis_op_error(.Generic, if (inherits(e1, c("basismfd", "mvbasismfd"))) e1 else e2)
}

#' @rdname Ops.basismfd
#' @export
"[.basismfd" <- function(x, ...) basis_op_error("[", x)

# mean(), sd(), scale() and plot() of basis objects
#' @exportS3Method NULL
op_mean.basismfd <- function(x, ...) basis_op_error("mean", x)
#' @exportS3Method NULL
op_sd.basismfd <- function(x, ...) basis_op_error("sd", x)
#' @exportS3Method NULL
op_scale.basismfd <- function(x, ...) basis_op_error("scale", x)
#' @exportS3Method NULL
op_plot.basismfd <- function(x, ...) basis_op_error("plot", x)
#' @exportS3Method NULL
op_mean.mvbasismfd <- function(x, ...) basis_op_error("mean", x)
#' @exportS3Method NULL
op_sd.mvbasismfd <- function(x, ...) basis_op_error("sd", x)
#' @exportS3Method NULL
op_scale.mvbasismfd <- function(x, ...) basis_op_error("scale", x)
#' @exportS3Method NULL
op_plot.mvbasismfd <- function(x, ...) basis_op_error("plot", x)
