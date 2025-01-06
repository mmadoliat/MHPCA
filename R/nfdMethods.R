mean_nfd <- function(nfd_obj, ...){
  nfd(t(colMeans(nfd_obj,...)))
}

sd_nfd <- function(nfd_obj, ...){
  nfd(t(apply(nfd_obj,2,sd,...)))
}

#' Print Method for `nfd` Objects
#'
#' This function provides a custom print method for objects of class `nfd`.
#' It displays a header indicating that the object is an `nfd` object
#' and then proceeds to print the underlying matrix data dimension.
#'
#' @param x An object of class `nfd`.
#' @param ... Additional arguments passed to the default print method.
#'
#' @return Invisibly returns the input object `x` after printing.
#'
#' @examples
#' # Create a standard matrix
#' mat <- matrix(1:9, nrow = 3, ncol = 3)
#'
#' # Convert the matrix to an 'nfd' object
#' nfd_obj <- nfd(mat)
#'
#' # Print the 'nfd' object
#' print(nfd_obj)
#'
#' @export
print.nfd = function(x,...) {
  cat("A ",ncol(x),"-Dimensional 'nfd' object:", sep="")
  cat("\nnobs:", nrow(x), "\n")
  #NextMethod("print")
  invisible(x)
}

#' Accessor for nfd Class
#'
#' @param x An object of class `nfd`.
#' @param name The name of the element to access.
#'
#' @return The requested element or `NULL` if it doesn't exist.
#' @export
#'
#' @examples
#' mat <- matrix(1:9, nrow = 3)
#' nfd_obj <- nfd(mat)
#' nfd_obj$nobs # Returns 3
`$.nfd` <- function(x, name) {
  if (name == "nobs") {
    return(nrow(x))
  } else if (name == "nfeatures") {
    return(ncol(x))
  } else if (name == "features") {
    if (!is.null(colnames(x))) {
      return(colnames(x))
    } else {
      return(paste0("V"), seq_len(ncol(x)))
    }
  } else if (name == "data"){
    return(unclass(x))
  } else {
    stop(paste("Unknown field:", name))
  }
}

#' Subset Method for 'nfd' Objects
#'
#' This function allows subsetting of 'nfd' objects using the `[ ]` operator.
#' The result of the subsetting operation is returned as an 'nfd' object.
#'
#' @param x An object of class `nfd`.
#' @param i (optional) Numeric, logical, or character vector specifying the rows to extract.
#' @param j (optional) Numeric, logical, or character vector specifying the columns to extract.
#' @param drop Logical indicating whether to drop dimensions (default is `FALSE`).
#' @param ... Further arguments passed to or from other methods.
#'
#' @return An `nfd` object containing the subset of the original data.
#' @export
#'
#' @examples
#' # Define the 'nfd' class constructor
#' nfd <- function(data) {
#'   if (!is.matrix(data) && !is.array(data)) {
#'     stop("Input must be a matrix or an array.")
#'   }
#'   class(data) <- c("nfd", class(data))
#'   return(data)
#' }
#'
#' # Create an 'nfd' object
#' mat <- matrix(1:12, nrow = 3, ncol = 4)
#' nfd_obj <- nfd(mat)
#'
#' # Subset rows and columns
#' subset_nfd <- nfd_obj[1:2, 3:4]
#' print(subset_nfd)
#' print(class(subset_nfd))  # Should include 'nfd'
`[.nfd` <- function(x, ...) {
  # Perform the subsetting using the underlying matrix or array
  subset_data <- NextMethod("[")
  
  # Check if the subsetted data is still a matrix or array
  if (is.matrix(subset_data) || is.array(subset_data)) {
    # If so, convert it back to an 'nfd' object
    return(nfd(subset_data))
  } else {
    # If not (i.e., it's a vector), return it as-is
    return(subset_data)
  }
}

#' @title Compute the inner product between two objects of class `nfd`
#'
#' @param nfd_obj1 An `nfd` object
#' @param nfd_obj2 An `nfd` object
#' @return The inner products matrix between the two `nfd` objects
#' @seealso \code{\link{nfd}}
#' @export
inprod_nfd <- function(nfd_obj1, nfd_obj2) {
  if (nfd_obj1$nfeatures != nfd_obj2$nfeatures) stop("number of features must be same")
  return(nfd_obj1%*%t(nfd_obj2))
}

#' @title Compute the norm objects of object `nfd`
#'
#' @param nfd_obj An `nfd` object
#' @return The norm of `nfd` objects
#' @seealso \code{\link{nfd}}
#' @export
norm_nfd <- function(nfd_obj1) {
  return(sqrt(diag(inprod_nfd(nfd_obj1,nfd_obj1))))
}