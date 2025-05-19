#' @title Define a Set of Hybrid Data objects
#'
#' @description
#' The `hd` class represents functional data ...

#' @field mf a `mvmfd` object
#' @field nf  a `mvnfd` object
#'
#' @examples
#' require(fda)
#' bs1 <- create.fourier.basis(c(0, 2 * pi), 5)
#' bs2 <- create.bspline.basis(c(0, 1), 7)
#' bs3 <- create.exponential.basis(c(0, 2), 3)
#' nobs <- 10
#' argval1 <- seq(0, 2 * pi, length.out = 12)
#' X1 <- outer(sin(argval1), seq(0.5, 1.5, length.out = nobs))
#' mdbs1 <- Basismfd(bs1)
#' mfd1 <- Mfd(argval1, X1, mdbs1)
#' mdbs2 <- Basismfd(bs1)
#' argval2 <- argval1
#' X2 <- outer(cos(argval2), seq(0.2, 1.5, length.out = nobs))
#' mfd2 <- Mfd(argval2, X2, mdbs1)
#' hd1 <- Hd(mfd1, mfd2)
#' hd1[1]
#' hd1[1, 1]
#' hd1[1:5, 2]
#' hd1[, 1]
#' hd1[1:5, ]
#' evalarg <- list(argval1, argval2)
#' hd1$eval(evalarg)
#' hd1 + hd1
#' mean(hd1)
#' inprod_hd(hd1, hd1)
#' norm_hd(hd1)
#' plot(hd1)
#' bimfdplot(hd1)
#'
#' @import R6
#' @importFrom fda is.basis eval.basis Data2fd
#'
#' @export
hd <- R6::R6Class("hd",
                  public = list(
                    #' @description
                    #' Constructor for `hd` objects (same as 'Hd')
                    #'
                    #' @param ... A `mfd` or `nfd` objects which have separated by comma
                    #' @param spars_par assign sparsity parameters of observations to hd object
                    #'
                    initialize = function(..., spars_par = NULL) {
                      hd_list <- list(...)
                      if (is.list(hd_list[[1]])) hd_list <- hd_list[[1]]
                      init_hd_list_check(hd_list)
                      mfd_list <- list()
                      nfd_list <- list()
                      mvmfd_list <- list()
                      mvnfd_list <- list()
                      for (i in seq_along(hd_list)) {
                        if(inherits(hd_list[[i]],"mfd")) {
                          mfd_list[[length(mfd_list)+1]] <- hd_list[[i]]
                        } else if (inherits(hd_list[[i]],"mvmfd")){
                          mvmfd_list[[length(mvmfd_list)+1]] <- hd_list[[i]]
                        } else if(inherits(hd_list[[i]],"nfd")) {
                          nfd_list[[length(nfd_list)+1]] <- hd_list[[i]]
                        } else if(inherits(hd_list[[i]],"mvnfd")) {
                          mvnfd_list[[length(nfd_list)+1]] <- hd_list[[i]]
                        }
                      }
                      if (length(mfd_list) > 0 & length(mvmfd_list) > 0){
                        stop("Only one object type can be provided: either `mfd` or `mvmfd`, not both.")
                      } 
                      if (length(nfd_list) > 0 & length(mvnfd_list) > 0){
                        stop("Only one object type can be provided: either `nfd` or `mvnfd`, not both.")
                      }
                      if (length(mvmfd_list) > 1){
                        stop("Only one `mvmfd' object type can be provided")
                      }
                      if (length(mvnfd_list) > 1){
                        stop("Only one `mvnfd' object type can be provided")
                      }
                      if (length(nfd_list) > 0) mvnfd_list <- list(Mvnfd(nfd_list))
                      if (length(mfd_list) > 0) mvmfd_list <- list(Mvmfd(mfd_list))
                      if (length(mvnfd_list) > 0) private$.nf <- mvnfd_list[[1]] 
                      if (length(mvmfd_list) > 0) private$.mf <- mvmfd_list[[1]]  
                      if (!is.null(spars_par)) self$spars_par <- spars_par
                    },
                    
                    #' @description
                    #' Eval method for `hd` objects
                    #'
                    #' @param evalarg A list of numeric vectors of argument values at which the `hd` is to be evaluated.
                    #' @return A list of evaluated values
                    eval = function(evalarg) {
                      if (!is.null(private$.mf)) {
                        fdata <- private$.mf$eval(evalarg)
                        names(fdata) <- paste0("fd",seq_along(fdata))
                      } else {
                          fdata <- NULL
                        }
                      if (!is.null(private$.nf)) {
                        nfdata <- private$.nf$data
                        names(nfdata) <- paste0("nfd",seq_along(nfdata))
                      } else {
                          nfdata <- NULL
                        }
                      return(c(fdata,nfdata))
                    },
                    #' @description
                    #' Print method for `hd` objects
                    #'
                    #' @param ... Additional arguments to be passed to `print`
                    #'
                    print = function(...) {
                      cat("An 'hd' object contains:\n")
                      if (!is.null(private$.mf)) print(private$.mf)
                      cat("\n")
                      if (!is.null(private$.nf)) print(private$.nf)
                      invisible(self)
                    }),
                  active = list(
                    # basis field
                    nf = function(value) {
                      if (missing(value)) {
                        private$.nf
                      } else {
                        stop("`$nf` is read only", call. = FALSE)
                      }
                    },
                    mf = function(value) {
                      if (missing(value)) {
                        private$.mf
                      } else {
                        stop("`$mf` is read only", call. = FALSE)
                      }
                    },
                    spars_par = function(value) {
                      if (missing(value)) {
                        # getter
                        return(private$.spars_par)
                      }
                      # setter: 1) must be integer
                      if (!is.numeric(value) || any(value != as.integer(value))) {
                        stop("`spars_par` must be an integer vector.", call. = FALSE)
                      }
                      value <- as.integer(value)
                      
                      # 2) max index must be less than total nbasis
                      if (!is.null(self$mf)) {
                        total_nb <- self$mf$nobs
                      } else if (!is.null(self$nf)) {
                        total_nb <- self$nf$nobs
                      } else {
                        stop("`hd` has neither `mf` nor `nf`.", call. = FALSE)
                      }
                      if (length(value) && max(value) >= total_nb) {
                        stop(
                          sprintf("all(spars_par) must be < %d (total number of basis functions)", 
                                  total_nb),
                          call. = FALSE
                        )
                      }
                      
                      # assign into private storage
                      private$.spars_par <- value
                    }
                  ),
                  private = list(
                    .mf = NULL,
                    .nf = NULL,
                    .spars_par = NULL
                  )
)
#' @rdname hd
#' @seealso \code{\link{mvbasismfd}}, \code{\link{mfd}}

#' @title A Class of Multivariate Multidimensional Functional Data objects
#'
#' @description
#' Constructor for `hd` objects (same as `Hd`)
#'
#' @param ... A `mfd` objects which have separated by comma
#' @param spars_par assign sparsity parameters of observations to hd object
#' @export
Hd <- function(..., spars_par = NULL) hd$new(..., spars_par = spars_par)
