
#' Print Method for `img` Objects
#'
#' @param x a matrix
#' @param ... other args
#'
#' @return Invisibly returns `x`.
#' @export
#'
#' @examples
#'  x <- matrix(rnorm(100),ncol = 10)
#' img_obj <- img(x)
#' print(img_obj)
print.img <- function(x, ...) {
  cat("An object of class 'img':\n")
  
  orig_dim <- attr(x, "original_dimensions")
  if (!is.null(orig_dim)) {
    cat("Original Dimensions: ", paste(orig_dim, collapse = " x "), "\n", sep = "")
  }
  
  # Print dimensions of the vectorized data
  cat("Vectorized Dimensions: ", paste(dim(x), collapse = " x "), "\n", sep = "")
  
  # Print class hierarchy
  cat("Classes: ", paste(class(x), collapse = ", "), "\n", sep = "")
  
  # Optionally, print first few rows and columns
  num_rows <- nrow(x)
  num_cols <- ncol(x)
  max_print_rows <- 5
  max_print_cols <- 5
  
  cat("\nFirst few entries:\n")
  
  if (num_rows > max_print_rows) {
    rows_to_print <- 1:max_print_rows
    cat("Showing first", max_print_rows, "images out of", num_rows, "\n")
  } else {
    rows_to_print <- 1:num_rows
  }
  
  if (num_cols > max_print_cols) {
    cols_to_print <- 1:max_print_cols
    cat("Showing first", max_print_cols, "pixels out of", num_cols, "\n")
  } else {
    cols_to_print <- 1:num_cols
  }
  
  print(x$data[rows_to_print, cols_to_print, drop = FALSE])
  
  if (num_rows > max_print_rows || num_cols > max_print_cols) {
    cat("... (truncated)\n")
  }
  
  invisible(x)
}

#' Check if an object is of class 'img'
#'
#' @param x The object to check.
#' @return TRUE if the object is of class 'img', FALSE otherwise.
#' @export
is.img <- function(x) {
  inherits(x, "img")
}

#' Original dimensions of the images in an `img` object
#'
#' @param x An object of class `img`.
#' @return The dimensions of a single image before vectorization.
#' @examples
#' img_obj <- img(matrix(1:6, 2), matrix(7:12, 2))
#' dim_original(img_obj)
#' @export
dim_original <- function(x) {
  if (!is.img(x)) {
    stop("Object is not of class 'img'.")
  }
  attr(x, "original_dimensions")
}
