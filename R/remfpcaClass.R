#' @title A Class for `ReMFPCA` objects
#'
#' @description
#' The `remfpca` class represents regularized functional principal components.
#'

#' @field pc_mfd An object of class `mvmfd` where the first indices (fields)
#' represents harmonics and  second indices represents variables
#' @field lsv = Left singular values vectors
#' @field values = The set of eigenvalues
#' @field smooth_tuning = The list of smoothing penalties parameters
#' @field sparse_tuning = The list of sparse penalties parameters
#' @field GCVs = Generalized cross validations scores of smoothing penalties parameters.
#'               If both smoothing and sparse tuning penalties are used in the ReMFPCA method, 
#'               this represents the conditional generalized cross-validation scores, which 
#'               means it is computed based on the optimal sparse tuning parameter selected via cross validation.
#' @field CVs = Cross validations scores of sparse penalties parameters
#' @field mean_mfd A multivariate functional data object giving the mean function

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
#' # Smoothing parameters selected by GCV from a grid for each variable
#' Re0 <- Remfpca(mvmfd_obj, ncomp = k, smooth_tuning = alpha1)
#' fpc0 <- Re0$pc_mfd
#' scores0 <- inprod_mvmfd(mvmfd_obj, fpc0)
#' dim(scores0)
#' Re0$smooth_tuning
#' # Fixed smoothing parameters, one per component
#' Re1 <- Remfpca(mvmfd_obj, ncomp = k, smooth_tuning = c(0, 0), smooth_GCV = FALSE)
#' Re1$values
#' # A single functional variable
#' Re2 <- Remfpca(mfd1, ncomp = k, smooth_tuning = list(alpha1$a1))
#' Re2$smooth_tuning

#' @import R6
#' @importFrom fda getbasispenalty
#' @importFrom Matrix bdiag
#' @importFrom expm sqrtm
#' @export
remfpca <- R6::R6Class("remfpca",
                             public = list(
                               # initialize = function(mvmfd_obj, method = "power", ncomp, smooth_tuning = NULL, sparse_tuning = 0, centerfns = TRUE, alpha_orth = FALSE, smoothing_type = "coefpen", sparse_type = "soft", K_fold = 30, sparsity_CV = "marginal") {
                               initialize = function(mvmfd_obj, method = "power", ncomp, smooth_tuning = NULL, sparse_tuning = NULL, centerfns = TRUE, alpha_orth = FALSE, smoothing_type = "coefpen", sparse_type = "soft", K_fold = 30, sparse_CV, smooth_GCV) {
      
                                 # if (is.numeric(smooth_tuning)) smooth_tuning <- as.list(smooth_tuning)
                                 # if (is.vector(smooth_tuning)) smooth_tuning <- as.list(smooth_tuning)
                                 # if (is.vector(smooth_tuning)&& !is.list(smooth_tuning)) smooth_tuning <- list(smooth_tuning)
                                 if (is.mfd(mvmfd_obj)) mvmfd_obj <- mvmfd$new(mvmfd_obj)
                                 method <- check_choice(method, c("power", "eigen"), "method")
                                 smoothing_type <- check_choice(smoothing_type, c("coefpen", "basispen"), "smoothing_type")
                                 sparse_type <- check_choice(sparse_type, sparse_penalty_types, "sparse_type")
                                 alpha_orth <- check_flag(alpha_orth, "alpha_orth")
                                 sequential <- method == "power" && !alpha_orth
                                 
                                 if (sequential) {
                                   smooth_tuning <- normalize_tuning(
                                     smooth_tuning, "smooth_tuning", nvar = mvmfd_obj$nvar, ncomp = ncomp,
                                     grid = smooth_GCV, mode = "per_component",
                                     vec_label = "'mvmfd_obj$nvar'", list_label = "'mvmfd_obj$nvar'",
                                     allow_matrix = TRUE
                                   )
                                   sparse_tuning <- normalize_tuning_u(sparse_tuning, "sparse_tuning", ncomp, sparse_CV)
                                   
                                   result <- sequential_power(mvmfd_obj = mvmfd_obj, n = ncomp, smooth_tuning = smooth_tuning, sparse_tuning=sparse_tuning, centerfns = centerfns, alpha_orth = alpha_orth, smooth_tuning_type = smoothing_type, sparse_tuning_type = sparse_type, K_fold = K_fold, sparse_CV, smooth_GCV)
                                 } 
                                 
                                 else {
                                   smooth_tuning <- normalize_tuning(
                                     smooth_tuning, "smooth_tuning", nvar = mvmfd_obj$nvar, ncomp = ncomp,
                                     grid = smooth_GCV, mode = "per_variable",
                                     vec_label = "number of variables", list_label = "'mvmfd_obj$nvar'",
                                     allow_matrix = TRUE
                                   )
                                   if (method == "power") {
                                     result <- joint_power(mvmfd_obj = mvmfd_obj, n = ncomp, smooth_tuning = smooth_tuning, centerfns = centerfns, alpha_orth = alpha_orth, smooth_tuning_type = smoothing_type)
                                   } else{
                                     result <- eigen_approach(mvmfd_obj = mvmfd_obj, n = ncomp, alpha = smooth_tuning, centerfns = centerfns, penalty_type = smoothing_type)
                                   }
                                   
                                   # result <- joint_power(mvmfd_obj = mvmfd_obj, n = ncomp, smooth_tuning = smooth_tuning, centerfns = centerfns, alpha_orth = alpha_orth, smooth_tuning_type = smoothing_type)
                                 } 
                                 # else if (method == "eigen") {
                                 #   result <- eigen_approach(mvmfd_obj = mvmfd_obj, n = ncomp, alpha = smooth_tuning, centerfns = centerfns, penalty_type = smoothing_type)
                                 # }
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
                                 if (sequential) {
                                   private$.sparse_tuning <- result[[5]]
                                   private$.CVs <- result[[6]]
                                   private$.GCVs <- result[[7]]
                                 } else{
                                   private$.GCVs <- result[[5]]
                                 }
                                 # private$.sparse_tuning <- result[[5]]
                                 private$.mean_mfd <- mean(mvmfd_obj)
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

#' @rdname remfpca
#' @seealso \code{\link{mvmfd}}

#' @title A Class for 'remfpca' objects
#'
#' @description
#' The `remfpca` class represents regularized functional principal components ('ReMFPCs') components.
#'
#' @param mvmfd_obj An `mvmfd` object representing the multivariate functional data.
#' @param method A character string specifying the approach to be used for MFPCA computation. 
#'               Options are "power" (the default), which uses the power algorithm, or "eigen", 
#'               which uses the eigen decomposition approach.
#' @param ncomp The number of functional principal components to retain.
#' @param smooth_tuning A list or vector specifying the smoothing regularization parameter(s) for each variable. 
#'                      If NULL, non-smoothing MFPCA is estimated.
#' @param sparse_tuning A list or vector specifying the sparsity regularization parameter(s) for each variable. 
#'                      If NULL, non-sparse MFPCA is estimated.
#' @param centerfns Logical indicating whether to center the functional data before analysis. Default is TRUE.
#' @param alpha_orth Logical indicating whether to perform orthogonalization of the regularization parameters. 
#'                   If `method` is "power", setting `alpha_orth = FALSE` (default) uses the sequential power approach, 
#'                   while setting `alpha_orth = TRUE` uses the joint power approach.
#' @param smoothing_type The type of smoothing penalty to be applied on the coefficients. The types "coefpen" and "basispen" is supported. Default is "coefpen".
#' @param sparse_type The type of sparse penalty to be applied on the coefficients. The types "soft", "hard", "SCAD" and "adaptive" are supported. Default is "soft".
#' @param K_fold  An integer specifying the number of folds in the sparse cross-validation process. Default is 30. 
#' @param sparse_CV @param sparse_CV Logical indicating whether cross-validation should be applied to select the optimal sparse tuning parameter in sequential power approach. 
#'                                        If `sparse_CV = TRUE`, a series of tuning parameters should be provided as a vector with positive number with max equals to number of subjects. 
#'                                        If `sparse_CV = FALSE`, specific tuning parameters are given directly to each principal components. Tuning parameters should be provided as a vector with length equal to `ncomp`.
#'                                        If the dimensions of input tuning parameters are incorrect, it will be converted to a list internally, and a warning will be issued. 
#' @param smooth_GCV @param smooth_GCV Logical indicating whether generalized cross-validation should be applied to select the optimal smooth tuning parameter. 
#'                                        If `smooth_GCV = TRUE`, a series of tuning parameters should be provided as a list with length equal to the number of variables. 
#'                                        If a list with incorrect dimensions is provided, it will be converted to a correct list internally, and a warning will be issued. 
#'                                        If `smooth_GCV = FALSE`, specific tuning parameters are given directly. If `method` is "power" and `alpha_orth = FALSE` (sequential power), 
#'                                        tuning parameters should be provided as a list with length equal to the number of variables, where each element is a vector of length `ncomp`. 
#'                                        If `method` is "power" and `alpha_orth = TRUE` (joint power), tuning parameters should be provided as a vector with length equal to the number of variables.
#'                                        If the dimensions of input tuning parameters are incorrect, it will be converted to a list internally, and a warning will be issued. 
#' @export
Remfpca <- function(mvmfd_obj, method = "power", ncomp, smooth_tuning = NULL, sparse_tuning = NULL, centerfns = TRUE, alpha_orth = FALSE, smoothing_type = "basispen", sparse_type = "soft", K_fold=30, sparse_CV = TRUE, smooth_GCV = TRUE) {
  remfpca$new(mvmfd_obj, method, ncomp, smooth_tuning, sparse_tuning, centerfns, alpha_orth, smoothing_type, sparse_type, K_fold, sparse_CV, smooth_GCV)
}
#' @rdname remfpca
