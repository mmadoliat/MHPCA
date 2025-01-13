#' @title Define a Set of Multivariate Image Objects
#'
#' @description
#' The `mvimg` class represents a set of multivariate image data objects.
#' Each `img` object within `mvimg` corresponds to a different variable (e.g., different time points, modalities).
#'
#' @field data A list of `img` objects.
#' @field nobs Number of observations (pixels).
#' @field nvar Number of variables (images).
#'
#' @examples
#' \dontrun{
#' # Assuming 'img_obj1' and 'img_obj2' are previously created 'img' objects
#' mvimg_obj <- Mvimg(img_obj1, img_obj2)
#' print(mvimg_obj)
#' mvimg_obj$data
#' mvimg_obj$nobs
#' mvimg_obj$nvar
#' }
#'
#' @import R6
#' @export
mvimg <- R6::R6Class(
  "mvimg",
  public = list(
    #' @description
    #' Constructor for `mvimg` objects (same as 'Mvimg')
    #'
    #' @param ... Multiple `img` objects to include in the multivariate image set.
    #'
    initialize = function(...) {
      img_list <- list(...)
      
      # If the first argument is a list, assume it's a list of img objects
      if (length(img_list) == 1 && is.list(img_list[[1]])) {
        img_list <- img_list[[1]]
      }
      
      # Validate that all elements are 'img' objects
      if (!all(sapply(img_list, function(x) inherits(x, "img")))) {
        stop("All inputs must be of class 'img'.")
      }
      
      # Ensure all 'img' objects have the same number of observations (pixels)
      nobs_values <- sapply(img_list, function(x) ncol(x))
      if (length(unique(nobs_values)) != 1) {
        stop("All 'img' objects must have the same number of observations (pixels).")
      }
      
      # Assign nobs and nvar
      private$.nobs <- nobs_values[1]
      private$.nvar <- length(img_list)
      
      # Store the data
      private$.data <- img_list
    },
    
    #' @description
    #' Print method for `mvimg` objects
    #'
    #' @param ... Additional arguments to be passed to `print`
    #'
    print = function(...) {
      cat("A 'mvimg' object with", private$.nvar, "variable(s):\n")
      cat("Number of observations (pixels):", private$.nobs, "\n")
      for (i in seq_along(private$.data)) {
        cat("\nVariable ", i, ":\n", sep = "")
        print(private$.data[[i]])
      }
      invisible(self)
    }
  ),
  active = list(
    # Getter for `data` field
    data = function(value) {
      if (missing(value)) {
        private$.data
      } else {
        stop("`$data` is read-only", call. = FALSE)
      }
    },
    
    # Getter for `nvar` field
    nvar = function(value) {
      if (missing(value)) {
        private$.nvar
      } else {
        stop("`$nvar` is read-only", call. = FALSE)
      }
    },
    
    # Getter for `nobs` field
    nobs = function(value) {
      if (missing(value)) {
        private$.nobs
      } else {
        stop("`$nobs` is read-only", call. = FALSE)
      }
    }
  ),
  private = list(
    .data = list(),    # List to store 'img' objects
    .nobs = NULL,      # Number of observations (pixels)
    .nvar = NULL       # Number of variables (images)
  )
)

#' @title A Class of Multivariate Image Data Objects
#'
#' @description
#' Constructor for `mvimg` objects (same as Mvimg(...))
#'
#' @param ... Multiple `img` objects to include in the multivariate image set.
#'
#' @export
Mvimg <- function(...) mvimg$new(...)
