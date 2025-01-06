#' @title Define a Set of Multivariate Multidimensional Vector Data objects
#'
#' @description
#' The `mvnfd` class represents a set of multivariate multidimensional vector data.

#' @field data A numeric matrix of data points.
#' @field nobs Number of observations.
#' @field nvar number of variables

#' @examples
#' # Create some vector data
#' data_matrix1 <- matrix(rnorm(100), nrow = 10, ncol = 10)
#' data_matrix2 <- matrix(rnorm(100), nrow = 10, ncol = 10)

#' nfd_obj1 <- Nfd(data = data_matrix1)
#' nfd_obj2 <- Nfd(data = data_matrix2)
#' mvnfd_obj <- Mvnfd(nfd_obj1,nfd_obj2)
#' print(mvnfd_obj)
#' nfd_obj$data
#' nfd_obj$nobs
#'
#' @import R6
#'
#' @export
mvnfd <- R6::R6Class("mvnfd",
                     public = list(
                       #' @description
                       #' Constructor for `mvnfd` objects (same as 'Mvnfd')
                       #'
                       #' @param ... A `nfd` objects which have separated by comma
                       #'
                       initialize = function(...) {
                         nfd_list <- list(...)
                         if (is.list(nfd_list[[1]])) nfd_list <- nfd_list[[1]]
                         init_nfd_list_check(nfd_list)
                         private$.nobs <- nfd_list[[1]]$nobs
                         private$.nvar <- length(nfd_list)
                         for (i in 1:private$.nvar) {
                           private$.data[[i]] <- unclass(nfd_list[[i]])
                           private$.features[[i]] <- nfd_list[[i]]$features
                         }
                       },
                       
                       #' @description
                       #' Print method for `mvnfd` objects
                       #'
                       #' @param ... Additional arguments to be passed to `print`
                       #'
                       print = function(...) {
                         cat("A 'mvnfd' object with", private$.nvar, "variable(s):\n")
                         for (i in 1:private$.nvar) {
                           cat("\nVariable ", i, ":\n", sep = "")
                           print(self[, i])
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
                       
                       features = function(value){
                         if (missing(value)){
                           private$.features
                         } else {
                           stop("`$features` is read-only", call. = FALSE)
                         }
                       }, 
                       # nvar field
                       nvar = function(value) {
                         if (missing(value)) {
                           private$.nvar
                         } else {
                           stop("`$nvar` is read only", call. = FALSE)
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
                       .features = list(),
                       .data = list(), # we record data
                       .nobs = NULL,
                       .nvar = NULL
                     )
)
#' @title A Class of Multivariate Multidimensional Vector Data Objects
#'
#' @description
#' Constructor for `mvnfd` objects (same as Mvnfd(...))
#'
#' @rdname mvnfd
#' @param ... A `nfd` objects which have separated by comma
#' @export
Mvnfd <- function(...) mvnfd$new(...)