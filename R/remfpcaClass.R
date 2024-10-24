#' @title A Class for `ReMFPCA` objects
#'
#' @description
#' The `remfpca` class represents regularized functional principal components components.
#'

#' @field pc_mfd an object of class `mvmfd` where the first indices (fields)
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
#' mfd1 <- Mfd(X = X_1, mdbs = mdbs)
#' mfd2 <- Mfd(X = X_2, mdbs = mdbs)
#' mvmfd_obj <- Mvmfd(mfd1, mfd2)
#' k <- 2
#' Re0 <- Remfpca(mvmfd_obj, ncomp = k, alpha = c(0, 0))
#' fpc0 <- Re0$pc_mfd
#' scores0 <- inprod_mvmfd(mvmfd_obj, fpc0)
#' dim(scores0)
#' Re0$alpha
#' Re1 <- Remfpca(mvmfd_obj, ncomp = k, alpha = alpha1)
#' Re1$alpha
#' Re3 <- Remfpca(mfd1, ncomp = k, alpha = alpha1$a1)
#' Re3$alpha

#' @import R6
#' @importFrom fda getbasispenalty
#' @importFrom Matrix bdiag
#' @importFrom expm sqrtm
#' @export
remfpca <- R6::R6Class("remfpca",
                             public = list(
                               initialize = function(mvmfd_obj, 
                                                     method = "power", 
                                                     ncomp, 
                                                     smooth_tuning = NULL,
                                                     sparse_tuning = 0, 
                                                     centerfns = TRUE, 
                                                     alpha_orth = FALSE,
                                                     smooth_tuning_type = "coefpen",
                                                     sparse_tuning_type = "SCAD",
                                                     K_fold = 30,
                                                     cv_type = "marginal",
                                                     vdata = NULL) {
                                 if (is.numeric(smooth_tuning)) smooth_tuning <- as.list(smooth_tuning)
                                 if (is.mfd(mvmfd_obj)) mvmfd_obj <- mvmfd$new(mvmfd_obj)
                                 if (method == "power") {
                                   result <- s_power_algorithm(mvmfd_obj = mvmfd_obj, 
                                                               n = ncomp, 
                                                               smooth_tuning = smooth_tuning,
                                                               sparse_tuning=sparse_tuning, 
                                                               centerfns = centerfns, 
                                                               alpha_orth = alpha_orth,
                                                               smooth_tuning_type = smooth_tuning_type,
                                                               sparse_tuning_type = sparse_tuning_type,
                                                               K_fold = K_fold,
                                                               cv_type,
                                                               vdata = vdata)
                                 }
                                 coef <- result[[1]]
                                 pcmfd <- list()
                                 for (i in 1:mvmfd_obj$nvar) {
                                   if (mvmfd_obj$basis$dimSupp[i] == 1) {
                                     pcmfd[[i]] <- Mfd(X = coef[[i]], mdbs = mvmfd_obj$basis$basis[[i]], method = "coefs")
                                   } else {
                                     coef_new <- array(coef[[i]], dim = c(mvmfd_obj$basis$basis[[i]]$nbasis, ncol(coef[[i]])))
                                     pcmfd[[i]] <- Mfd(X = coef_new, mdbs = mvmfd_obj$basis$basis[[i]], method = "coefs")
                                   }
                                 }
                                 out <- Mvmfd(pcmfd)
                                 if(mvmfd_obj$nvar==1) {
                                   private$.pc_mfd <- pcmfd[[1]]
                                 } else {
                                   private$.pc_mfd <- Mvmfd(pcmfd) 
                                 }
                                 private$.lsv <- result[[2]]
                                 private$.values <- result[[3]]
                                 private$.smooth_tuning <- result[[4]]
                                 private$.sparse_tuning <- result[[5]]
                                 private$.mean_mfd <- mean(mvmfd_obj)
                                 private$.vector_pc <- result[[6]]
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
                                   stop("`$GCVs` is read only", call. = FALSE)
                                 }
                               },
                               mean_mfd = function(value) {
                                 if (missing(value)) {
                                   private$.mean_mfd
                                 } else {
                                   stop("`$mean_mfd` is read only", call. = FALSE)
                                 }
                               } , 
                               vector_pc = function(values){
                                 if (missing(values)) {
                                   private$.vector_pc
                                 } else {
                                   stop("`$vector_pc` is read only", call. = FALSE)
                                 }
                               }
                             ),
                             private = list(
                               .pc_mfd = NULL,
                               .lsv = NULL,
                               .values = NULL,
                               .smooth_tuning = NULL,
                               .sparse_tuning = NULL,
                               .mean_mfd = NULL,
                               .vector_pc = NULL 
                             )
)

#' @rdname remfpca
#' @seealso \code{\link{mvmfd}}

#' @title A Class for 'remfpca' objects
#'
#' @description
#' The `remfpca` class represents regularized functional principal components ('ReMFPCs') components.
#'
#' @param mvmfd_obj An `mvmfd` object representing the multivariate functional data.
#' @param ncomp The number of functional principal components to retain.
#' @param alpha A list or vector specifying the regularization parameter(s) for each variable.
#'              If NULL, the regularization parameter is estimated internally.
#' @param centerfns Logical indicating whether to center the functional data before analysis.
#' @param alpha_orth Logical indicating whether to perform orthogonalization of the regularization parameters.
#' @param penalty_type The type of penalty to be applied on the coefficients. The types "coefpen" and "basispen" is supported. Default is "coefpen".
#' @export
Remfpca <- function(mvmfd_obj, method = "power", ncomp, smooth_tuning = NULL,sparse_tuning, centerfns = TRUE, alpha_orth = FALSE,smooth_tuning_type = "coefpen",sparse_tuning_type = "soft",K_fold=30,cv_type = "marginal",vdata=NULL) {
  remfpca$new(mvmfd_obj, method, ncomp, smooth_tuning,sparse_tuning, centerfns, alpha_orth,smooth_tuning_type,sparse_tuning_type,K_fold,cv_type,vdata = vdata)
}
#' @rdname remfpca
