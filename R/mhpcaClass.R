#' @title A Class for `MHPCA` objects
#'
#' @description
#' The `mhpca` class represents regularized functional principal components components.
#'

#' @field pc_mfd an object of class `hd` where the first indices (fields)
#' represents harmonics and  second indices represents variables
#' @field lsv = Left singular values vectors
#' @field values = the set of eigenvalues
#' @field alpha = The vector of penalties parameters
#' @field GCVs = generalized cross validations
#' @field mean_mfd a multivariate functional data object giving the mean function

#' @examples
#' require(fda)
#' # Brownian Bridge simulation on [0,1]
#' M <- 110 # number of components
#' N <- 20 # number of instances
#' n <- 100 # number of grides
#' t0 <- seq(0, 1, len = n)
#' j <- 1:M
#' alpha1 <- list(a1 = 2^seq(0, 1, length.out = 3), a2 = 2^seq(0, 1, length.out = 3))
#' psi_1 <- function(t, m) sin(m * pi * t) # eigenfunction of BB
#' psi_2 <- function(t, m) sin((2 * m - 1) * pi / 2 * t) # eigenfunction of BM
#' PC_1 <- outer(t0, j, FUN = psi_1) # n by M matrix
#' PC_2 <- outer(t0, j, FUN = psi_2) # n by M matrix
#' Z <- matrix(rnorm(N * M), nr = M)
#' lambda <- matrix(2 / (pi * (2 * j - 1)), nr = M, nc = N)
#' X_1t <- PC_1 %*% (lambda * Z)
#' X_2t <- PC_2 %*% (lambda * Z)
#' noise <- rnorm(n * N, 0, 0.1)
#' X_1 <- X_1t + noise
#' X_2 <- X_2t + noise
#' bs <- create.bspline.basis(c(0, 1), 51)
#' mdbs <- Basismfd(bs)
#' mfd1 <- Mfd(X = X_1, mdbs = mdbs) #123
#' mfd2 <- Mfd(X = X_2, mdbs = mdbs)
#' hd_obj <- Hd(mfd1, mfd2)
#' k <- 2
#' Re0 <- Mhpca(hd_obj, ncomp = k, alpha = c(0, 0))
#' fpc0 <- Re0$pc_mfd
#' scores0 <- inprod_hd(hd_obj, fpc0)
#' dim(scores0)
#' Re0$alpha
#' Re1 <- Mhpca(hd_obj, ncomp = k, alpha = alpha1)
#' Re1$alpha
#' Re3 <- Mhpca(mfd1, ncomp = k, alpha = alpha1$a1)
#' Re3$alpha

#' @import R6
#' @importFrom fda getbasispenalty
#' @importFrom Matrix bdiag
#' @importFrom expm sqrtm
#' @export
mhpca <- R6::R6Class("mhpca",
                             public = list(
                               # initialize = function(hd_obj, method = "power", ncomp, smooth_tuning = NULL, sparse_tuning = 0, centerfns = TRUE, alpha_orth = FALSE, smoothing_type = "coefpen", sparse_type = "soft", K_fold = 30, sparsity_CV = "marginal") {
                               initialize = function(hd_obj, method = "power", ncomp, smooth_tuning = NULL, sparse_tuning = 0, centerfns = TRUE, alpha_orth = FALSE, smoothing_type = "coefpen", sparse_type = "soft", K_fold = 30, sparse_CV, smooth_GCV) {
                                 # if (is.numeric(smooth_tuning)) smooth_tuning <- as.list(smooth_tuning)
                                 # if (is.vector(smooth_tuning)) smooth_tuning <- as.list(smooth_tuning)
                                 # if (is.vector(smooth_tuning)&& !is.list(smooth_tuning)) smooth_tuning <- list(smooth_tuning)
                                 if (is.mfd(hd_obj)) hd_obj <- hd$new(hd_obj)
                                 
                                 if (method == "power" & alpha_orth == "FALSE") {
                                   # Adjust the vector length to match the required dimensions if they are incorrect
                                   if (is.vector(smooth_tuning)& !is.list(smooth_tuning)) {
                                     if (smooth_GCV == FALSE) {
                                     if (length(smooth_tuning) != ncomp) {
                                       warning("The length of 'smooth_tuning' did not match 'ncomp' and has been adjusted accordingly.", call. = FALSE)
                                       smooth_tuning <- rep(smooth_tuning, length.out = ncomp)
                                     }
                                     smooth_tuning <- replicate(hd_obj$nvar, smooth_tuning, simplify = FALSE)
                                     }
                                     else{
                                       warning("The length of 'smooth_tuning' did not match 'hd_obj$nvar' and has been adjusted accordingly.", call. = FALSE)
                                       smooth_tuning <- replicate(hd_obj$nvar, smooth_tuning, simplify = FALSE)
                                     }
                                   }
                                   
                                   # Adjust the matrix to match the required dimensions if they are incorrect
                                   else if (is.matrix(smooth_tuning)) {
                                     if (smooth_GCV == FALSE) {
                                     if (!all(dim(smooth_tuning) == c(hd_obj$nvar, ncomp))) {
                                       smooth_tuning <- smooth_tuning[rep(1:nrow(smooth_tuning), length.out = hd_obj$nvar), rep(1:ncol(smooth_tuning), length.out = ncomp)]
                                       # print(smooth_tuning)
                                       smooth_tuning <- split(smooth_tuning, row(smooth_tuning))
                                       # print(smooth_tuning)
                                       warning("The dimensions of 'smooth_tuning' did not match the expected size and have been adjusted accordingly.", call. = FALSE)
                                     } else{
                                       smooth_tuning <- split(smooth_tuning, row(smooth_tuning))
                                     }
                                     }
                                     else{
                                       if (dim(smooth_tuning)[1] != hd_obj$nvar) {
                                       smooth_tuning <- smooth_tuning[rep(1:nrow(smooth_tuning), length.out = hd_obj$nvar), , drop = FALSE][1:hd_obj$nvar, , drop = FALSE]
                                       smooth_tuning <- split(smooth_tuning, row(smooth_tuning))
                                       warning("The dimensions of 'smooth_tuning' did not match the expected size and have been adjusted accordingly.", call. = FALSE)
                                       }
                                       else{
                                         smooth_tuning <- split(smooth_tuning, row(smooth_tuning))
                                       }
                                     }
                                   }
                                   
                                   # Adjust the list length and element sizes to match the required dimensions if they are incorrect
                                   else if (is.list(smooth_tuning)) {
                                     if (smooth_GCV == FALSE) {
                                     if (length(smooth_tuning) != hd_obj$nvar) {
                                       warning("Adjusting 'smooth_tuning' to match 'hd_obj$nvar'.", call. = FALSE)
                                       smooth_tuning <- rep(smooth_tuning, length.out = hd_obj$nvar)
                                     }
                                     smooth_tuning <- lapply(smooth_tuning, function(vec) {
                                       if (length(vec) != ncomp) {
                                         warning("Adjusting vector length in 'smooth_tuning' to match 'ncomp'.", call. = FALSE)
                                         vec <- rep(vec, length.out = ncomp)
                                       }
                                       vec
                                     })
                                     }
                                     else{
                                       if (length(smooth_tuning) != hd_obj$nvar) {
                                         warning("Adjusting 'smooth_tuning' to match 'hd_obj$nvar'.", call. = FALSE)
                                         smooth_tuning <- rep(smooth_tuning, length.out = hd_obj$nvar)
                                       }
                                       
                                       # smooth_tuning <- rep(smooth_tuning, length.out = hd_obj$nvar)[1:hd_obj$nvar]
                                     }
                                   }
                                   
                                   names(smooth_tuning) <- paste0("var", 1:hd_obj$nvar)
                                   
                                   # Adjust the list length and element sizes to match the required dimensions if they are incorrect
                                   if (sparse_CV == FALSE & length(sparse_tuning) != ncomp) {
                                     warning("The length of 'sparse_tuning' did not match 'ncomp' and has been adjusted accordingly.", call. = FALSE)
                                     sparse_tuning <- rep(sparse_tuning, length.out = ncomp)
                                   }
                                   
                                   result <- sequential_power(hd_obj = hd_obj, n = ncomp, smooth_tuning = smooth_tuning, sparse_tuning=sparse_tuning, centerfns = centerfns, alpha_orth = alpha_orth, smooth_tuning_type = smoothing_type, sparse_tuning_type = sparse_type, K_fold = K_fold, sparse_CV, smooth_GCV)
                                 } 
                                 
                                 else if (method == "power" & alpha_orth == "TRUE") {
                                   # Adjust the vector to match the required lengths if they are incorrect
                                   if (is.vector(smooth_tuning) & !is.list(smooth_tuning)) {
                                     if (smooth_GCV == FALSE) {
                                     if (length(smooth_tuning) != hd_obj$nvar) {
                                       warning("The length of 'smooth_tuning' did not match number of variables and has been adjusted accordingly.", call. = FALSE)
                                       smooth_tuning <- rep(smooth_tuning, length.out = hd_obj$nvar)
                                     }
                                     smooth_tuning <- lapply(1:hd_obj$nvar, function(i) smooth_tuning[i])
                                     }
                                     else{
                                       warning("The length of 'smooth_tuning' did not match number of variables and has been adjusted accordingly.", call. = FALSE)
                                       smooth_tuning <- replicate(hd_obj$nvar, smooth_tuning, simplify = FALSE)
                                     }
                                   }
                                   
                                   # Adjust the matrix to match the required if they are incorrect
                                   else if (is.matrix(smooth_tuning)) {
                                     if (smooth_GCV == FALSE) {
                                     if (!all(dim(smooth_tuning) == c(hd_obj$nvar, 1))) {
                                       smooth_tuning <- smooth_tuning[rep(1:nrow(smooth_tuning), length.out = hd_obj$nvar), rep(1:ncol(smooth_tuning), length.out = 1)]
                                       smooth_tuning <- as.list(smooth_tuning)
                                       warning("The dimensions of 'smooth_tuning' did not match the expected size and have been adjusted accordingly.", call. = FALSE)
                                     } else{
                                       smooth_tuning <- as.list(smooth_tuning)
                                     }
                                     } 
                                     else{
                                       if (dim(smooth_tuning)[1] != hd_obj$nvar) {
                                       smooth_tuning <- smooth_tuning[rep(1:nrow(smooth_tuning), length.out = hd_obj$nvar), , drop = FALSE][1:hd_obj$nvar, , drop = FALSE]
                                       smooth_tuning <- split(smooth_tuning, row(smooth_tuning))
                                       warning("The dimensions of 'smooth_tuning' did not match the expected size and have been adjusted accordingly.", call. = FALSE)
                                       }
                                       else{
                                         smooth_tuning <- split(smooth_tuning, row(smooth_tuning))
                                       }
                                     }
                                   }
                                   
                                   # Adjust the list length and element sizes to match the required dimensions if they are incorrect
                                   else if (is.list(smooth_tuning)) {
                                     if (smooth_GCV == FALSE) {
                                     if (length(smooth_tuning) != hd_obj$nvar) {
                                       warning("Adjusting 'smooth_tuning' to match 'hd_obj$nvar'.", call. = FALSE)
                                       smooth_tuning <- rep(smooth_tuning, length.out = hd_obj$nvar)
                                     }
                                     smooth_tuning <- lapply(smooth_tuning, function(vec) {
                                       if (length(vec) != 1) {
                                         warning("Adjusting vector length in 'smooth_tuning' to match 'ncomp'.", call. = FALSE)
                                         vec <- rep(vec, length.out = 1)
                                       }
                                       vec
                                     })
                                     }
                                     else{
                                       if (length(smooth_tuning) != hd_obj$nvar) {
                                        warning("Adjusting 'smooth_tuning' to match 'hd_obj$nvar'.", call. = FALSE)
                                        smooth_tuning <- rep(smooth_tuning, length.out = hd_obj$nvar)[1:hd_obj$nvar]
                                       }
                                     }
                                   }
                                   names(smooth_tuning) <- paste0("var", 1:hd_obj$nvar)
                                   
                                   result <- joint_power(hd_obj = hd_obj, n = ncomp, smooth_tuning = smooth_tuning, centerfns = centerfns, alpha_orth = alpha_orth, smooth_tuning_type = smoothing_type)
                                 } else if (method == "eigen" ) {
                                   result <- eigen_approach(hd_obj = hd_obj, n = ncomp, alpha = smooth_tuning, centerfns = centerfns, penalty_type = smoothing_type)
                                 }
                                 coef <- result[[1]]
                                 pcmfd <- list()
                                 for (i in 1:hd_obj$nvar) {
                                   if (hd_obj$basis$dimSupp[i] == 1) {
                                     pcmfd[[i]] <- Mfd(X = coef[[i]], mdbs = hd_obj$basis$basis[[i]], method = "coefs")
                                   } else {
                                     coef_new <- array(coef[[i]], dim = c(hd_obj$basis$basis[[i]]$nbasis, ncol(coef[[i]])))
                                     pcmfd[[i]] <- Mfd(X = coef_new, mdbs = hd_obj$basis$basis[[i]], method = "coefs")
                                   }
                                 }
                                 out <- Hd(pcmfd)
                                 if(hd_obj$nvar==1) {
                                   private$.pc_mfd <- pcmfd[[1]]
                                 } else {
                                   private$.pc_mfd <- Hd(pcmfd) 
                                 }
                                 private$.lsv <- result[[2]]
                                 private$.values <- result[[3]]
                                 private$.smooth_tuning <- result[[4]]
                                 if (alpha_orth == "FALSE") {
                                   private$.sparse_tuning <- result[[5]]
                                   private$.CVs <- result[[6]]
                                   private$.GCVs <- result[[7]]
                                 } else{
                                   private$.GCVs <- result[[5]]
                                 }
                                 # private$.sparse_tuning <- result[[5]]
                                 private$.mean_mfd <- mean(hd_obj)
                               }
                             ),
                             active = list(
                               pc_mfd = function(value) {
                                 if (missing(value)) {
                                   private$.pc_mfd
                                 } else {
                                   stop("`$pc_mfd` is read only", call. = FALSE)
                                 }
                               },
                               lsv = function(value) {
                                 if (missing(value)) {
                                   private$.lsv
                                 } else {
                                   stop("`$lsv` is read only", call. = FALSE)
                                 }
                               },
                               values = function(value) {
                                 if (missing(value)) {
                                   private$.values
                                 } else {
                                   stop("`$sigma` is read only", call. = FALSE)
                                 }
                               },
                               smooth_tuning = function(value) {
                                 if (missing(value)) {
                                   private$.smooth_tuning
                                 } else {
                                   stop("`$smooth_tuning` is read only", call. = FALSE)
                                 }
                               },
                               sparse_tuning = function(value) {
                                 if (missing(value)) {
                                   private$.sparse_tuning
                                 } else {
                                   stop("`$sparse_tuning` is read only", call. = FALSE)
                                 }
                               },
                               GCVs = function(value) {
                                 if (missing(value)) {
                                   private$.GCVs
                                 } else {
                                   stop("`$GCVs` is read only", call. = FALSE)
                                 }
                               },
                               CVs = function(value) {
                                 if (missing(value)) {
                                   private$.CVs
                                 } else {
                                   stop("`$CVs` is read only", call. = FALSE)
                                 }
                               },
                               mean_mfd = function(value) {
                                 if (missing(value)) {
                                   private$.mean_mfd
                                 } else {
                                   stop("`$mean_mfd` is read only", call. = FALSE)
                                 }
                               }
                             ),
                             private = list(
                               .pc_mfd = NULL,
                               .lsv = NULL,
                               .values = NULL,
                               .smooth_tuning = NULL,
                               .sparse_tuning = NULL,
                               .GCVs = NULL,
                               .CVs = NULL,
                               .mean_mfd = NULL
                             )
)

#' @rdname mhpca
#' @seealso \code{\link{hd}}

#' @title A Class for 'mhpca' objects
#'
#' @description
#' The `mhpca` class represents regularized functional principal components ('ReMFPCs') components.
#'
#' @param hd_obj An `hd` object representing the multivariate functional data.
#' @param ncomp The number of functional principal components to retain.
#' @param alpha A list or vector specifying the regularization parameter(s) for each variable.
#'              If NULL, the regularization parameter is estimated internally.
#' @param centerfns Logical indicating whether to center the functional data before analysis.
#' @param alpha_orth Logical indicating whether to perform orthogonalization of the regularization parameters.
#' @param penalty_type The type of penalty to be applied on the coefficients. The types "coefpen" and "basispen" is supported. Default is "coefpen".
#' @export
Mhpca <- function(hd_obj, method = "power", ncomp, smooth_tuning = NULL, sparse_tuning, centerfns = TRUE, alpha_orth = FALSE, smoothing_type = "coefpen", sparse_type = "soft", K_fold=30, sparse_CV = TRUE, smooth_GCV = TRUE) {
  mhpca$new(hd_obj, method, ncomp, smooth_tuning, sparse_tuning, centerfns, alpha_orth, smoothing_type, sparse_type, K_fold, sparse_CV, smooth_GCV)
}
#' @rdname mhpca
