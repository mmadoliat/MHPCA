mean_mfd <- function(mfd_obj) {
  # average over observations (the last dimension; 2-D variables have 3)
  cof <- apply(mfd_obj$coefs, basis_margins(mfd_obj$coefs), mean)
  bs <- mfd_obj$basis
  return(Mfd(X = cof, mdbs = bs, method = "coefs"))
}


sd_mfd <- function(mfd_obj) {
  cof <- apply(mfd_obj$coefs, basis_margins(mfd_obj$coefs), sd)
  bs <- mfd_obj$basis
  return(Mfd(X = cof, mdbs = bs, method = "coefs"))
}


plot_mfd <- function(mfd_obj, obs = 1, xlab = "", ylab = "", main = "", type = "l", lty = 1, ...) {
  dimSupp <- mfd_obj$basis$dimSupp
  supp <- mfd_obj$basis$supp
  x_grids <- seq(supp[1, 1], supp[2, 1], len = 1000)
  if (dimSupp == 1) {
    X <- mfd_obj$eval(x_grids)
    matplot(x_grids, X, type = type, lty = lty, xlab = xlab, ylab = ylab, main = main, ...)
  } else {
    y_grids <- seq(from = supp[1, 2], to = supp[2, 2], length.out = 100)
    X <- mfd_obj$eval(list(x_grids, y_grids))[, , obs]
    image(X, xlab = xlab, ylab = ylab, axes = FALSE, main = paste(main, " Observation:", obs))
    axis(side = 1, at = seq(from = 0, to = 1, length.out = 10), labels = round(seq(supp[1, 1], supp[2, 1], len = 10), 1))
    axis(side = 2, at = seq(from = 0, to = 1, length.out = 10), labels = round(seq(supp[1, 2], supp[2, 2], len = 10), 1))
  }
}


#' @title Length of an object of classes `mfd`, `mvmfd`, `nfd` or `hd`.
#'
#' @description
#' Length of an object of an object of classes `mfd`, `mvmfd`,`nfd` or `hd`.
#'
#' @param x An object of classes `mfd`, `mvmfd`, `nfd` or `hd`.
#' @param ... all `length` function arguments.
#' @export
length <- function(x, ...) {
  op_length(x, ...)
}

#' @title plots an object of classes `mfd`, `mvmfd`, `hd`, `remfpca` or `mhpca`
#'
#' @description
#'  plot an object of classes `mfd`, `mvmfd`, `hd`, `remfpca` or `mhpca`
#' @param x An object of classes `mfd`, `mvmfd`, `hd`, `remfpca` or `mhpca`
#' @param ... all `plot` function arguments.
#' @export
plot <- function(x, ...) {
  op_plot(x, ...)
}
# mean should change to include vector data as well 
#' @title mean of an object of classes `mfd`, `mvmfd`, `nfd`, `mvnfd`, `mvmfd` or `hd`.
#'
#' @description
#' mean of an object of classes `mfd`, `mvmfd`, `nfd`, `mvnfd`, `mvmfd` or `hd`.
#' @param x An object of classes `mfd`, `mvmfd`, `nfd`, `mvnfd`, `mvmfd` or `hd`.
#' @param ... all `mean` function arguments.
#' @return An object of class `mfd`
#' @export
mean <- function(x, ...) {
  op_mean(x, ...)
}

#' @title Standard deviation of an object of class `mfd`, `mvmfd`, `nfd`, `mvnfd`, `mvmfd` or `hd`..
#'
#' @description
#' Standard deviation an object of class `mfd`, `mvmfd`, `nfd`, `mvnfd`, `mvmfd` or `hd`..
#' @param x An object of class `mfd`, `mvmfd`, `nfd`, `mvnfd`, `mvmfd` or `hd`.
#' @param ... all `sd` function arguments.
#' @return An object of class `mfd`, `mvmfd`, `nfd`, `mvnfd`, `mvmfd` or `hd`.
#' @export
sd <- function(x, ...) {
  op_sd(x, ...)
}

#' @title Compute the inner product between two objects of class `mfd`
#'
#' @param mfd_obj1 An `mfd` object
#' @param mfd_obj2 An `mfd` object
#' @return The inner products matrix between the two `mfd` objects
#' @seealso \code{\link{basismfd}}, \code{\link{mfd}}
#' @importFrom fda fd inprod
#' @export
inprod_mfd <- function(mfd_obj1, mfd_obj2) {
  bs1 <- mfd_obj1$basis$basis[[1]]
  bs2 <- mfd_obj2$basis$basis[[1]]
  cof1 <- mfd_obj1$coefs
  cof2 <- mfd_obj2$coefs
  fd1 <- fd(coef = cof1, basisobj = bs1)
  fd2 <- fd(coef = cof2, basisobj = bs2)
  inpr <- fda::inprod(fd1, fd2)
  return(inpr)
}

#' @title Compute the norm of an object of class `mfd` 
#'
#' @param mfd_obj An object of class `mfd` 
#' @seealso \code{\link{basismfd}}, \code{\link{mfd}}
#' @return The norm vector of the an object of class `mfd` 
#' @export
norm_mfd <- function(mfd_obj) {
  return(as.numeric(sqrt(diag(inprod_mfd(mfd_obj, mfd_obj)))))
}

#' Scale an `mfd` Object
#'
#' This function scales an `mfd` object by calculating a scaling factor based on the variance of its evaluations 
#' or using a provided weight. It returns a new scaled `mfd` object.
#'
#' @param mfd_obj An object of class `mfd`.
#' @param mfd_eval_length The number of evaluation points to use for scaling.
#' @param weight An optional numeric value to use as the scaling factor. If NULL, the scaling factor is calculated automatically.
#'
#' @return A scaled mfd object.
#'
#' @examples
#' # Example usage:
#' # Assuming `mfd_obj` is a valid mfd object:
#' # scaled_mfd <- scale_mfd(mfd_obj, mfd_eval_length = 100)
#' # scaled_mfd <- scale_mfd(mfd_obj, mfd_eval_length = 100, weight = 0.5)
scale_mfd <- function(mfd_obj,mfd_eval_length = 100,weight = NULL){
  if (!inherits(mfd_obj,"mfd")) stop("The object must be of the class 'mfd'.")
  v1 <- mfd_obj$eval(do.call(seq,c(as.list(mfd_obj$basis$supp),length.out = mfd_eval_length)))
  if (is.null(weight)){
    scaling_factor <- 1/sqrt(mean(diag(var(v1))))
  } else {
    scaling_factor <- weight
  }
  # Scaling is linear, so scale the coefficients and keep the original basis.
  return(mfd$new(X = scaling_factor * mfd_obj$coefs, mdbs = mfd_obj$basis, method = "coefs"))
}

#' @title Add two `mfd` objects
#'
#' @param obj1 An `mfd` object
#' @param obj2 An `mfd` object or a scalar
#' @return The sum of the two `mfd` objects
#' @seealso \code{\link{basismfd}}, \code{\link{mfd}}
#' @export
"+.mfd" <- function(obj1, obj2 = NULL) {
  if (is.null(obj2)) {
    return(obj1)
  }
  if (xor(is.mfd(obj1), is.mfd(obj2))) {
    if (!xor(is.numeric(obj1), is.numeric(obj2))) stop("At least one object must be a mfd, and the other one can be a scalar")
    if (is.numeric(obj1)) {
      temp <- obj1
      obj1 <- obj2
      obj2 <- temp
    }
    coef <- obj1$coefs + obj2
  } else {
    if (length(obj1) == length(obj2)) {
      if (obj1$basis$dimSupp != obj2$basis$dimSupp) stop("Two objects must have the same basis dimSupp.")
      if (obj1$basis$dimSupp == 1) {
        if (obj1$basis$basis[[1]]$type != obj2$basis$basis[[1]]$type) stop("Two objects must have the same basis types.")
      } else {
        for (i in 1:obj1$basis$dimSupp) {
          if (obj1$basis$basis[[i]]$type != obj2$basis$basis[[i]]$type) stop("Two objects must have the same basis types.")
        }
      }
      coef <- obj1$coefs + obj2$coefs
    } else {
      if (length(obj1) > 1 & length(obj2) > 1) stop("Two objects must have the same basis dimSupp.")
      if (length(obj1) == 1 & length(obj2) > 1) {
        temp <- obj1
        obj1 <- obj2
        obj2 <- temp
      }
      nobs <- obj1$nobs
      I_n <- matrix(1L, ncol = nobs, nrow = 1)
      coef <- obj1$coefs + as.matrix(obj2$coefs) %*% I_n
    }
  }
  return(mfd$new(X = coef, mdbs = obj1$basis$clone(), method = "coefs"))
}

#' @title Scalar multiplication of an `mfd` object
#' 
#' @description
#' Scalar multiplication of an `mfd` object. One object must be an `mfd`, and the other one a scalar
#' 
#' @param obj1 An `mfd` object or an scalar
#' @param obj2 An `mfd` object or an scalar
#' @return An `mfd` object
#' @seealso \code{\link{basismfd}}, \code{\link{mfd}}
#' @export
"*.mfd" <- function(obj1, obj2) {
  if (xor(is.mfd(obj1), is.mfd(obj2))) {
    if (xor(is.numeric(obj1), is.numeric(obj2))) {
      if (is.numeric(obj1)) {
        temp <- obj1
        obj1 <- obj2
        obj2 <- temp
      }
    }
    coef <- obj1$coefs * obj2
  } else {
    stop("One object must be an mfd, and the other one a scalar")
  }
  return(mfd$new(X = coef, mdbs = obj1$basis$clone(), method = "coefs"))
}

#' @title Subtract two `mfd` objects
#'
#' @param obj1 An `mfd` object
#' @param obj2 An `mfd` object or a scalar
#' @return The difference between the two `mfd` objects
#' @seealso \code{\link{basismfd}}, \code{\link{mfd}}
#' @export
"-.mfd" <- function(obj1, obj2 = NULL) {
  if (is.null(obj2)) {
    return((-1) * obj1)
  }
  if (xor(is.mfd(obj1), is.mfd(obj2))) {
    if (!xor(is.numeric(obj1), is.numeric(obj2))) stop("At least one object must be an mfd, and the other one can be a scalar")
    if (is.numeric(obj1)) {
      coef <- obj1 + (-1) * obj2$coefs
      obj1 <- obj2
    } else {
      coef <- obj1$coefs + (-obj2)
    }
  } else {
    return(obj1 + (-1) * obj2)
  }
  return(mfd$new(X = coef, mdbs = obj1$basis$clone(), method = "coefs"))
}

#' @title Extract subsets of an `mfd` object
#'
#' @param mfd_obj An `mfd` object
#' @param i An index or indices specifying the subsets to extract
#' @param j Variable index. An `mfd` holds a single variable, so `j` can only
#'   be left empty or 1: `mfd_obj[1:3, ]` and `mfd_obj[1:3, 1]` are the same as
#'   `mfd_obj[1:3]`.
#' @return An `mfd` object containing the specified subsets
#' @seealso \code{\link{basismfd}}, \code{\link{mfd}}
#' @export
"[.mfd" <- function(mfd_obj, i = "index", j) {
  if (!missing(j) && !identical(as.numeric(j), 1)) {
    stop("An `mfd` object holds a single variable, so its second index can only be 1 (or left empty).",
         call. = FALSE)
  }
  if (is.null(i) || identical(i, "index")) i <- 1:mfd_obj$nobs
  if (max(i) > mfd_obj$nobs | min(i) < 1) stop(" subscript i out of bounds")
  bs <- mfd_obj$basis$clone()
  if (mfd_obj$basis$dimSupp == 1) {
    coef <- mfd_obj$coefs[, i]
  } else {
    coef <- mfd_obj$coefs[, , i]
  }
  return(mfd$new(X = coef, mdbs = bs, method = "coefs"))
}

#' @title scale of an object of class `mfd`, `mvmfd`, `nfd`, `mvnfd` or `hd`..
#'
#' @description
#' Standard deviation an object of class `mfd`, `mvmfd`, `nfd`, `mvnfd` or `hd`..
#' @param x An object of class `mfd`, `mvmfd`, `nfd`, `mvnfd` or `hd`.
#' @param ... all `scale` function arguments.
#' @return An object of class `mfd`, `mvmfd`, `nfd`, `mvnfd` or `hd`.
#' @export
scale <- function(x, ...) {
  op_scale(x, ...)
}

# Internal generics behind the exported length(), plot(), mean(), sd() and
# scale(). To support a new class, define op_<name>.<class> next to that
# class's other methods; objects without a method fall back to base/stats.
# These are deliberately not methods of base's own generics: registering e.g.
# length.nfd would change what base R itself sees for these objects.
op_length <- function(x, ...) UseMethod("op_length")
op_plot <- function(x, ...) UseMethod("op_plot")
op_mean <- function(x, ...) UseMethod("op_mean")
op_sd <- function(x, ...) UseMethod("op_sd")
op_scale <- function(x, ...) UseMethod("op_scale")

#' @exportS3Method NULL
op_length.default <- function(x, ...) base::length(x, ...)
#' @exportS3Method NULL
op_plot.default <- function(x, ...) base::plot(x, ...)
#' @exportS3Method NULL
op_mean.default <- function(x, ...) base::mean(x, ...)
#' @exportS3Method NULL
op_sd.default <- function(x, ...) stats::sd(x, ...)
#' @exportS3Method NULL
op_scale.default <- function(x, ...) base::scale(x, ...)

#' @exportS3Method NULL
op_length.mfd <- function(x, ...) x$nobs
#' @exportS3Method NULL
op_plot.mfd <- function(x, ...) plot_mfd(x, ...)
#' @exportS3Method NULL
op_mean.mfd <- function(x, ...) mean_mfd(x, ...)
#' @exportS3Method NULL
op_sd.mfd <- function(x, ...) sd_mfd(x, ...)
#' @exportS3Method NULL
op_scale.mfd <- function(x, ...) scale_mfd(x, ...)

# Margins of a coefficient array that index basis functions: all dimensions
# except the last one, which indexes observations.
basis_margins <- function(coefs) seq_len(length(dim(coefs)) - 1)
