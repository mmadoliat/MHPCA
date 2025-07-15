mean_hd <- function(hd_obj) {
  mf <- NULL
  nf <- NULL
  if (!is.null(hd_obj$mf)) mf <- mean(hd_obj$mf)
  if (!is.null(hd_obj$nf)) nf <- mean(hd_obj$nf)
  mvlist <- Filter(Negate(is.null),list(mf,nf))
  return(Hd(mvlist))
}

sd_hd <- function(hd_obj) {
  mf <- NULL
  nf <- NULL
  if (!is.null(hd_obj$mf)) mf <- sd(hd_obj$mf)
  if (!is.null(hd_obj$nf)) nf <- sd(hd_obj$nf)
  mvlist <- Filter(Negate(is.null),list(mf,nf))
  return(Hd(mvlist))
}

plot_hd <- function(hd_obj, xlab = NULL, ylab = NULL, ...) {
  if (is.null(hd_obj$mf)) stop("hd objects doesn't have functional part")
  plot(hd_obj$mf, ...)
}


#' Scale a `hd` Object
#'
#' This function scales a `hd` object, which may consist of multivariate functional 
#' data `mvmfd` and/or `mvnfd`. The function ensures proper scaling 
#' for all components based on their respective types and returns a scaled `hd` object.
#'
#' @param hd_obj A `hd` object, which can be an `mfd`, `mvmfd`, `nfd`, `mvnfd`.
#' @param mfd_eval_length A numeric vector specifying the number of evaluation points for functional variables.
#' @param f_weight A numeric vector of scaling factors for functional variables. The length must match number of functional variables.
#' @param nf_weight A numeric vector of scaling factors for non-functional variables. The length must match the number of non functional variables.

#' @return A scaled high-dimensional data object of class 'Hd'.
#'
#' @examples
#' # Example usage:
#' # Assuming `hd_obj` is a valid high-dimensional data object:
#' # scaled_hd <- scale_hd(hd_obj, mfd_eval_length = c(100, 200), weight = c(0.5, 0.8, 1.2))
scale_hd <- function(hd_obj,
                     method = c("component", "balanced","manual"),
                     f_weight = NULL,
                     nf_weight = NULL,
                     mf_eval_length = NULL) {
  method <- match.arg(method)
  mf <- hd_obj$mf
  nf <- hd_obj$nf
  if (method == "component") {
    mf_eval_length <- if (is.null(mf_eval_length)) rep(100, mf$nvar) else mf_eval_length
    mvmfd_scaled <- if (!is.null(mf)) scale_mvmfd(mf, mf_eval_length, f_weight)
    mvnfd_scaled <- if (!is.null(nf)) scale_mvnfd(nf, nf_weight)
    
  } else if (method == "balanced") {
    if (is.null(mf) || is.null(nf)) stop("Balanced scaling requires both mf and nf")
    w_top    <- sum(norm_mvmfd(mf - mean(mf)))
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
    if (xor(is.double(obj1), is.double(obj2))) {
      if (is.double(obj1)) {
        temp <- obj1
        obj1 <- obj2
        obj2 <- temp
      }
    }
    mf <- NULL
    nf <- NULL
    if (!is.null(obj1$mf)) mf <- obj2 * obj1$mf
    if (!is.null(obj1$nf)) nf <- obj2 * obj1$nf
    mv_list <- list(Filter(Negate(is.null), list(mf,nf)))
  } else {
    stop("One object must be an hd, and the other one a scalar")
  }
  return(Hd(mv_list))
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
  jj <- 0
  if (!is.null(hd_obj$mf)){
    nobs <- hd_obj$mf$nobs
    jj <- jj+1
  }
  if (!is.null(hd_obj$nf)){
    nobs <- hd_obj$nf$nobs
    jj <- jj+1
  }
  if (is.null(i)) i <- 1:nobs
  if (is.null(j)) j <- 1:jj
  if (max(i) > nobs | min(i) < 1) stop(" subscript i out of bounds")
  if (max(j) > jj | min(i) < 1) stop(" subscript j out of bounds")
  if (length(j) == 2){
    mf_obj_subset <- list(hd_obj$mf[i,])
    nf_obj_subset <- list(hd_obj$nf[i,])
  } else if (length(j) == 1){
    if (!is.null(hd_obj$mf)) {
      mf_obj_subset <- list(hd_obj$mf[i,])
      nf_obj_subset <- NULL
    } else {
        nf_obj_subset <- list(hd_obj$nf[i,])
        mf_obj_subset <- NULL
        }
  }
  mv_list <- Filter(Negate(is.null),c(mf_obj_subset,nf_obj_subset))
  return(Hd(mv_list))
  
}


#' @title   Bivariate plot for `hd` objects
#'
#' @param hd_obj An `hd` object
#' @param type Type of plot ('l' for lines, 'p' for points, etc.)
#' @param lty Line type
#' @param xlab Label for the x-axis
#' @param ylab Label for the y-axis
#' @param main Main title
#' @param ... Additional arguments for the matplot function
#' @seealso \code{\link{hd}}, \code{\link{mvbasismfd}}
#'
#' @export
bimfdplot <- function(hd_obj, type = "l", lty = 1, xlab = "", ylab = "", main = "", ...) {
  if (is.null(hd_obj$mf)) stop("There is no functional part in hf object")
  mf_obj <- hd_obj$mf
  nvar <- mf_obj$nvar
  stopifnot(nvar == 2)
  stopifnot(all(mf_obj$basis$supp[[1]] == mf_obj$basis$supp[[2]]))
  supp <- mf_obj$basis$supp[[1]]
  x_grids <- seq(supp[1, 1], supp[2, 1], len = 1000)
  X <- mf_obj[, 1]$eval(x_grids)
  Y <- mf_obj[, 2]$eval(x_grids)
  matplot(X, Y, type = type, lty = lty, xlab = xlab, ylab = ylab, main = main, ...)
}


#' @title Compute the inner product between two objects of class `hd`
#'
#' @param hd_obj1 An `hd` object
#' @param hd_obj2 An `hd` object
#' @return The inner products matrix between the two `hd` objects
#' @seealso \code{\link{hd}},\code{\link{mvbasismfd}}
#' @export
inprod_hd <- function(hd_obj1, hd_obj2) {
  inpr <- 0
  if (!is.null(hd_obj1$mf) & !is.null(hd_obj1$mf)){
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
        Mfd(X = center_coef,mdbs = mhd$mf$basis[i],method = "coefs")
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
