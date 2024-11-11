mean_hd <- function(hd_obj) {
  p <- hd_obj$nvar
  mvlist <- lapply(1:p, function(j) mean(hd_obj[, j]))
  return(Hd(mvlist))
}


plot_hd <- function(hd_obj, xlab = NULL, ylab = NULL, ...) {
  old <- par()
  exclude_pars <- c("cin", "cra", "csi", "cxy", "din", "page")
  ind <- which(!(names(old) %in% exclude_pars))
  on.exit(par(old[ind]))
  p <- hd_obj$nvar
  par(mfrow = c(p, 1))
  if (is.null(ylab)) ylab <- paste("Variable ", 1:p)
  if (is.null(xlab)) xlab <- rep("time", p)
  for (i in 1:p) {
    plot(hd_obj[, i], ylab = ylab[i], xlab = xlab[i], ...)
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
  p <- obj1$nvar
  mvlist <- list()
  for (j in 1:p) {
    mvlist[[j]] <- obj1[, j] + obj2[, j]
  }
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
  p <- obj1$nvar
  mvlist <- list()
  for (j in 1:p) {
    mvlist[[j]] <- obj1[, j] + (-1) * obj2[, j]
  }
  return(Hd(mvlist))
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
    p <- obj1$nvar
    mvlist <- list()
    for (j in 1:p) {
      mvlist[[j]] <- obj1[, j] * obj2
    }
  } else {
    stop("One object must be an hd, and the other one a scalar")
  }
  return(Hd(mvlist))
}

#'  Extract subsets of an `hd` object
#'
#' @param hd_obj An `hd` object
#' @param i An index or indices specifying the subsets to extract for the first dimension
#' @param j An index or indices specifying the subsets to extract for the second dimension
#' @return An `hd` object containing the specified subsets
#' @seealso \code{\link{hd}},\code{\link{mvbasismfd}}
#' @export
"[.hd" <- function(hd_obj, i = "index", j = "index") {
  if (i[1] == "index") i <- 1:hd_obj$nobs
  if (j[1] == "index") j <- 1:hd_obj$nvar
  if (max(i) > hd_obj$nobs | min(i) < 1) stop(" subscript i out of bounds")
  if (max(j) > hd_obj$nvar | min(i) < 1) stop(" subscript j out of bounds")
  dimSupp <- hd_obj$basis$dimSupp
  bs <- hd_obj$basis[j]
  if (length(j) == 1) {
    if (dimSupp[j] == 1) {
      coef <- hd_obj$coefs[[j]][, i]
    } else {
      coef <- hd_obj$coefs[[j]][, , i]
    }
    return(mfd$new(X = coef, mdbs = bs, method = "coefs"))
  } else {
    hd_list <- list()
    for (k in 1:length(j)) {
      if (dimSupp[k] == 1) {
        coef <- hd_obj$coefs[[k]][, i]
      } else {
        coef <- hd_obj$coefs[[k]][, , i]
      }
      hd_list[[k]] <- mfd$new(X = coef, mdbs = bs[k], method = "coefs")
    }
    return(Hd(hd_list))
  }
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
  nvar <- hd_obj$nvar
  stopifnot(nvar == 2)
  stopifnot(all(hd_obj$basis$supp[[1]] == hd_obj$basis$supp[[2]]))
  supp <- hd_obj$basis$supp[[1]]
  x_grids <- seq(supp[1, 1], supp[2, 1], len = 1000)
  X <- hd_obj[, 1]$eval(x_grids)
  Y <- hd_obj[, 2]$eval(x_grids)
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
  p <- hd_obj1$nvar
  if (p != hd_obj2$nvar) stop("The number of variables must be equal.")
  m <- hd_obj1$nobs
  n <- hd_obj2$nobs
  inpr <- matrix(0, nrow = m, ncol = n)
  for (j in 1:p) {
    inpr <- inpr + inprod_mfd(hd_obj1[, j], hd_obj2[, j])
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
  return(as.numeric(sqrt(diag(inprod_hd(hd_obj, hd_obj)))))
}
