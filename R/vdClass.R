#' @title  Define a Set of Multidimensional Vector Data Objects
#'
#' @description
#' The `vd` class represents a set of multidimensional vector data.
#' Vector data objects are constructed by specifying a matrix of data points.
#' 
#' @field data A numeric matrix of data points.
#' @field nobs Number of observations.
#'
#' @examples
#' # Create some vector data
#' data_matrix <- matrix(rnorm(100), nrow = 10, ncol = 10)
#' vd_obj <- Vd(data = data_matrix)
#' print(vd_obj)
#' vd_obj$data
#' vd_obj$nobs
#' 
#' @import R6
#' @export
vd <- R6::R6Class("vd",
                  public = list(
                    #' @description
                    #' Constructor for `vd` objects (same as Vd(...))
                    #'
                    #' @param data A numeric matrix of data points.
                    initialize = function(data) {
                      if (!is.matrix(data)) {
                        stop("`data` must be a matrix", call. = FALSE)
                      }
                      private$.data <- data
                      private$.nobs <- nrow(data)
                      private$.ncols <- ncol(data)
                      
                    },
                    #' @description
                    #' Print method for `vd` objects.
                    #'
                    #' @param ... Additional arguments to be passed to `print`.
                    print = function(...) {
                      cat("A 'vd' object with", private$.nobs, "observations\n")
                      cat("Data dimensions:", dim(private$.data), "\n")
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
                    # Getter for `nobs` field
                    nobs = function(value) {
                      if (missing(value)) {
                        private$.nobs
                      } else {
                        stop("`$nobs` is read-only", call. = FALSE)
                      }
                    },
                    ncols = function(value){
                      if (missing(value)){
                        private$.ncols
                      } else {
                        stop("`$ncols` is read-only", call. = FALSE)
                      }
                    }
                  ),
                  private = list(
                    .data = NULL,
                    .nobs = NULL,
                    .ncols = NULL
                  )
)

#' @title A Class of Multidimensional Vector Data Objects
#'
#' @description
#' Constructor for `vd` objects (same as Vd(...))
#'
#' @rdname vd
#' @param data A numeric matrix of data points.
#' @export
Vd <- function(data) vd$new(data)
