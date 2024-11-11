#' Check if an object is of class 'basismfd'
#'
#' @param fdobj The object to check.
#' @return TRUE if the object is of class 'basismfd', FALSE otherwise.
#' @seealso \code{\link{is.mvbasismfd}}, \code{\link{is.mfd}}, \code{\link{is.hd}}
#' @export
is.basismfd <- function(fdobj) {
  inherits(fdobj, "basismfd")
}

#' Check if an object is of class 'mvbasismfd'
#'
#' @param fdobj The object to check.
#' @return TRUE if the object is of class 'mvbasismfd', FALSE otherwise.
#' @seealso \code{\link{is.basismfd}}, \code{\link{is.mfd}}, \code{\link{is.hd}}
#' @export
is.mvbasismfd <- function(fdobj) {
  inherits(fdobj, "mvbasismfd")
}

#' Check if an object is of class 'mfd'
#'
#' @param fdobj The object to check.
#' @return TRUE if the object is of class 'mfd', FALSE otherwise.
#' @seealso \code{\link{is.mvbasismfd}}, \code{\link{is.basismfd}}, \code{\link{is.hd}}
#' @export
is.mfd <- function(fdobj) {
  inherits(fdobj, "mfd")
}

#' Check if an object is of class 'vd'
#'
#' @param obj The object to check.
#' @return TRUE if the object is of class 'vd', FALSE otherwise.
#' @seealso \code{\link{is.vd}}, \code{\link{is.mfd}}, \code{\link{is.basismfd}}
#' @export
is.vd <- function(obj) {
  inherits(obj, "vd")
}

#' Check if an object is of class 'hd'
#'
#' @param obj The object to check.
#' @return TRUE if the object is of class 'hd', FALSE otherwise.
#' @seealso \code{\link{is.hd}}, \code{\link{is.mfd}}, \code{\link{is.basismfd}}
#' @export
is.hd <- function(obj) {
  inherits(obj, "hd")
}
