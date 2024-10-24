#' @import R6
#' @importFrom fda getbasispenalty
#' @importFrom Matrix bdiag
#' @importFrom expm sqrtm
#' @export
s_power_class <- R6::R6Class("s_power_class",
                             public = list(
                               initialize = function(mvmfd_obj, 
                                                     method = "power", 
                                                     ncomp , 
                                                     smooth_tuning = NULL,
                                                     sparse_tuning = NULL, 
                                                     centerfns = TRUE, 
                                                     alpha_orth = TRUE,
                                                     smooth_tuning_type = "coefpen",
                                                     sparse_tuning_type = "soft",
                                                     K_fold,
                                                     cv_type = "joint") {
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
                                                               cv_type)
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
                               }
                             ),
                             private = list(
                               .pc_mfd = NULL,
                               .lsv = NULL,
                               .values = NULL,
                               .smooth_tuning = NULL,
                               .sparse_tuning = NULL,
                               .mean_mfd = NULL
                             )
)


#' @title A Class for 'ReMFPCA' objects
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
Remfpca_fsvd <- function(mvmfd_obj, method = "power", ncomp, smooth_tuning = NULL,sparse_tuning = NULL, centerfns = TRUE, alpha_orth = TRUE,smooth_tuning_type = "coefpen",sparse_tuning_type = "soft",K_fold,cv_type = "joint") {
  s_power_class$new(mvmfd_obj, method = "power", ncomp, smooth_tuning = NULL,sparse_tuning = NULL, centerfns = TRUE, alpha_orth = TRUE,smooth_tuning_type = "coefpen",sparse_tuning_type = "soft",K_fold,cv_type = "joint")
# class
  }