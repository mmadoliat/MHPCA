# Apply f to each part of an hd object that is present (mf, then nf) and
# rebuild an hd from the results.
map_hd <- function(hd_obj, f) {
  parts <- list(if (!is.null(hd_obj$mf)) f(hd_obj$mf), if (!is.null(hd_obj$nf)) f(hd_obj$nf))
  Hd(Filter(Negate(is.null), parts))
}

mean_hd <- function(hd_obj) map_hd(hd_obj, mean)

sd_hd <- function(hd_obj) map_hd(hd_obj, sd)

plot_hd <- function(hd_obj, xlab = NULL, ylab = NULL, ...) {
  if (is.null(hd_obj$mf)) stop("hd objects doesn't have functional part")
  plot(hd_obj$mf, xlab = xlab, ylab = ylab, ...)
}


#' Scale a `hd` Object
#'
#' This function scales a `hd` object, which may consist of multivariate functional 
#' data `mvmfd` and/or `mvnfd`. The function ensures proper scaling 
#' for all components based on their respective types and returns a scaled `hd` object.
#'
#' @param hd_obj A `hd` object.
#' @param method How to scale: "component" (default) scales each variable
#'   separately; "balanced" rescales the non-functional part so that both
#'   parts have the same total variation (sum of squared norms of the centred
#'   observations); "manual" multiplies the functional and
#'   non-functional parts by `f_weight` and `nf_weight` (rescaled to sum to 1;
#'   if only one is given, the other is 1 minus it).
#' @param f_weight Scaling factor(s) for the functional part: one per
#'   functional variable for "component", a single weight for "manual".
#' @param nf_weight Scaling factor(s) for the non-functional part: one per
#'   non-functional variable for "component", a single weight for "manual".
#' @param mf_eval_length A numeric vector giving the number of evaluation
#'   points for each functional variable ("component" only; default 100 each).
#' @return A scaled `hd` object.
#'
#' @examples
#' # Assuming `hd_obj` is a valid hd object:
#' # scaled_hd <- scale_hd(hd_obj)
#' # scaled_hd <- scale_hd(hd_obj, method = "manual", f_weight = 0.5)
scale_hd <- function(hd_obj,
                     method = c("component", "balanced","manual"),
                     f_weight = NULL,
                     nf_weight = NULL,
                     mf_eval_length = NULL) {
  method <- match.arg(method)
  mf <- hd_obj$mf
  nf <- hd_obj$nf
  if (method == "component") {
    if (is.null(mf_eval_length) && !is.null(mf)) mf_eval_length <- rep(100, mf$nvar)
    mvmfd_scaled <- if (!is.null(mf)) scale_mvmfd(mf, mf_eval_length, f_weight)
    mvnfd_scaled <- if (!is.null(nf)) scale_mvnfd(nf, nf_weight)
    
  } else if (method == "balanced") {
    if (is.null(mf) || is.null(nf)) stop("Balanced scaling requires both mf and nf")
    # total variation (sum of squared norms of the centred data) of each part
    w_top    <- sum(norm_mvmfd(mf - mean(mf))^2)
    nf_new <- nfd(do.call("cbind",lapply(seq_len(nf$nvar),function(i) nf[,i])))
    w_bottom <- sum(apply(sweep(nf_new, 2, colMeans(nf_new), "-"), 2, norm, "2")^2)
    w_half   <- sqrt(w_top / w_bottom)
    mvmfd_scaled <- mf
    mvnfd_scaled <- w_half * nf
    
  } else {
    if (is.null(f_weight) && is.null(nf_weight)) stop("Manual method requires f_weight or nf_weight")
    if (!is.null(f_weight) && is.null(nf_weight)) nf_weight <- 1 - f_weight
    if (is.null(f_weight) && !is.null(nf_weight)) f_weight <- 1 - nf_weight
    if ((f_weight + nf_weight) != 1) {
      total     <- f_weight + nf_weight
      f_weight  <- f_weight  / total
      nf_weight <- nf_weight / total
    }
    mvmfd_scaled <- if (!is.null(mf)) f_weight * mf
    mvnfd_scaled <- if (!is.null(nf)) nf_weight * nf
  }
  
  if (!is.null(mf) && !is.null(nf)) {
    Hd(mvmfd_scaled, mvnfd_scaled)
  } else if (!is.null(mf)) {
    Hd(mvmfd_scaled)
  } else if (!is.null(nf)) {
    Hd(mvnfd_scaled)
  }
}



#' @title Addition of two `hd` objects
#'
#' @param obj1 An `hd` object
#' @param obj2 An optional `hd` object
#' @return An `hd` object
#' @seealso \code{\link{hd}},\code{\link{mvbasismfd}}
#' @importFrom graphics image axis par points matplot
#' @export
"+.hd" <- function(obj1, obj2 = NULL) {
  if (is.null(obj2)) {
    return(obj1)
  }
  if (!is.null(obj1$mf) & !is.null(obj2$mf)){
    mf <- obj1$mf + obj2$mf
  } else {
    mf <- NULL
  }
  if (!is.null(obj1$nf) & !is.null(obj2$nf)){
    nf <- obj1$nf + obj2$nf
  } else {
    nf <- NULL
  }
  mvlist <- Filter(Negate(is.null), list(mf,nf))
  return(Hd(mvlist))
}

#' @title  Subtraction of two `hd` objects
#'
#' @param obj1 An `hd` object
#' @param obj2 An optional `hd` object
#' @return An `hd` object
#' @seealso \code{\link{hd}},\code{\link{mvbasismfd}}
#' @export
"-.hd" <- function(obj1, obj2 = NULL) {
  if (is.null(obj2)) {
    return((-1) * obj1)
  }
  return(obj1 + (-1) * obj2)
}

#'  Multiplication of an `hd` object with a scalar
#'
#'
#' @param obj1 An `hd` object or a scalar
#' @param obj2 An `hd` object or a scalar
#' @return An `hd` object
#' @seealso \code{\link{hd}},\code{\link{mvbasismfd}}
#' @export
"*.hd" <- function(obj1, obj2) {
  if (xor(is.hd(obj1), is.hd(obj2))) {
    if (xor(is.numeric(obj1), is.numeric(obj2))) {
      if (is.numeric(obj1)) {
        temp <- obj1
        obj1 <- obj2
        obj2 <- temp
      }
    }
    return(map_hd(obj1, function(part) obj2 * part))
  } else {
    stop("One object must be an hd, and the other one a scalar")
  }
}

#'  Extract subsets of an `hd` object
#'
#' @param hd_obj An `hd` object
#' @param i An index or indices specifying the subsets to extract for the first dimension
#' @param j An index or indices specifying the subsets to extract for the second dimension
#' @return An `hd` object containing the specified subsets
#' @seealso \code{\link{hd}},\code{\link{mvbasismfd}}
#' @export
"[.hd" <- function(hd_obj, i = NULL, j = NULL) {
  # j indexes the parts that are present: the functional part (if any) first,
  # then the non-functional part.
  parts <- Filter(Negate(is.null), list(hd_obj$mf, hd_obj$nf))
  nobs <- parts[[1]]$nobs
  if (is.null(i)) i <- 1:nobs
  if (is.null(j)) j <- seq_along(parts)
  if (max(i) > nobs | min(i) < 1) stop(" subscript i out of bounds")
  if (max(j) > length(parts) | min(j) < 1) stop(" subscript j out of bounds")
  mv_list <- lapply(parts[j], function(part) part[i, ])
  return(Hd(mv_list))
  
}


# bimfdplot() for hd objects is handled by bimfdplot() in mvmfdMethods.R.


#' @title Compute the inner product between two objects of class `hd`
#'
#' @param hd_obj1 An `hd` object
#' @param hd_obj2 An `hd` object
#' @return The inner products matrix between the two `hd` objects
#' @seealso \code{\link{hd}},\code{\link{mvbasismfd}}
#' @export
inprod_hd <- function(hd_obj1, hd_obj2) {
  inpr <- 0
  if (!is.null(hd_obj1$mf) & !is.null(hd_obj2$mf)){
    inpr <- inpr + inprod_mvmfd(hd_obj1$mf,hd_obj2$mf)
  }
  if (!is.null(hd_obj1$nf) & !is.null(hd_obj2$nf)){
    inpr <- inpr + inprod_mvnfd(hd_obj1$nf,hd_obj2$nf)
  }
  
  return(inpr)
}

#' @title Compute the norm of an object of class `hd` 
#'
#' @param hd_obj An `hd` object
#' @return The norm vector of the an object of class `hd` 
#' @seealso \code{\link{hd}},\code{\link{mvbasismfd}}
#' @export
norm_hd <- function(hd_obj) {
  inprod_hd(hd_obj,hd_obj)
}

center_hd <- function(hd_obj){
  if (!is.hd(hd_obj)) {
    stop("Input 'hd_obj' must be of class 'hd'.")
  }
  mhd <- mean(hd_obj)
  mfd_list <- NULL
  nfd_list <- NULL
  if (!is.null(hd_obj$mf)){
    mfd_list <- lapply(seq_len(mhd$mf$nvar),function(i){
      if (is.matrix(hd_obj$mf$coefs[[i]])){
        coefs <- sweep(hd_obj$mf$coefs[[i]],1,mhd$mf$coefs[[i]],"-")
        Mfd(X = coefs,mdbs = mhd$mf$basis[i],method = "coefs")
      } else if (is.array(hd_obj$mf$coefs[[i]])){
        # I'm not sure about this code check it later for array 
        flattened_current <- apply(hd_obj$mf$coefs[[i]], 3, as.vector)
        flattened_mean <- apply(mhd$mf$coefs[[i]], 3, as.vector)
        row_means <- rowMeans(flattened_mean)
        centered_flattened <- sweep(flattened_current, 1, row_means, "-")
        centered_coef <- centered_flattened
        Mfd(X = array(centered_coef, dim = c(mhd$mf$basis$nbasis[[i]], ncol(centered_coef))), mdbs = mhd$mf$basis[i], method = "coefs")
      }
      
    })
  }
  if (!is.null(hd_obj$nf)){
    nfd_list <- lapply(seq_len(mhd$nf$nvar),function(i){
      if (is.matrix(hd_obj$nf$data[[i]])){
        data <- sweep(hd_obj$nf$data[[i]],2,mhd$nf$data[[i]],"-")
        nfd(data)
      } else if (is.array(hd_obj$nf$data[[i]])){
        # I'm not sure about this code check it later for array 
        flattened_current <- apply(hd_obj$nf$data[[i]], 3, as.vector)
        flattened_mean <- apply(mhd$nf$data[[i]], 3, as.vector)
        row_means <- rowMeans(flattened_mean)
        centered_flattened <- sweep(flattened_current, 1, row_means, "-")
        centered_data <- centered_flattened
        nfd(centered_data)
      }
      
    })
  }
  return(Hd(c(mfd_list,nfd_list)))
}

# Methods for the internal generics behind length(), plot(), mean(), sd(), scale()
# length(): number of observations (hd itself has no `nobs` field)
#' @exportS3Method NULL
op_length.hd <- function(x, ...) if (!is.null(x$mf)) x$mf$nobs else x$nf$nobs
#' @exportS3Method NULL
op_plot.hd <- function(x, ...) plot_hd(x, ...)
#' @exportS3Method NULL
op_mean.hd <- function(x, ...) mean_hd(x, ...)
#' @exportS3Method NULL
op_sd.hd <- function(x, ...) sd_hd(x, ...)
#' @exportS3Method NULL
op_scale.hd <- function(x, ...) scale_hd(x, ...)
