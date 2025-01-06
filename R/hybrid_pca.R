#' Factory function to create either mhpca or remfpca objects
#'
#' @param hd_obj 
#' @param method 
#' @param ncomp 
#' @param smooth_tuning 
#' @param sparse_tuning 
#' @param centerfns 
#' @param alpha_orth 
#' @param smoothing_type 
#' @param sparse_type 
#' @param K_fold 
#' @param sparse_CV 
#' @param smooth_GCV 
#' @param penalize_nfd 
#' @param penalize_u 
#'
#' @returns
#' @export
#'
#' @examples
hybrid_pca <- function(hd_obj, 
                       method = "power", 
                       ncomp = 3, 
                       smooth_tuning = NULL, 
                       sparse_tuning = NULL, 
                       centerfns = TRUE, 
                       alpha_orth = FALSE, 
                       smoothing_type = "coefpen", 
                       sparse_type = "soft", 
                       K_fold = 30, 
                       sparse_CV = TRUE, 
                       smooth_GCV = TRUE, 
                       penalize_nfd = FALSE, 
                       penalize_u = FALSE) {
  
  # Define your criteria here
  if (is.mfd(hd_obj) || is.mvmfd(hd_obj)) {
    if (is.mfd(hd_obj)) {
      hd_obj <- Mvmfd(hd_obj)
    }
    # Instantiate and return a Remfpca object
    return(Remfpca(
      mvmfd_obj = hd_obj,
      method = method,
      ncomp = ncomp,
      smooth_tuning = smooth_tuning,
      sparse_tuning = sparse_tuning,
      centerfns = centerfns,
      alpha_orth = alpha_orth,
      smoothing_type = smoothing_type,
      sparse_type = sparse_type,
      K_fold = K_fold,
      sparse_CV = sparse_CV, 
      smooth_GCV = smooth_GCV
    ))
  } else {
    # Instantiate and return an mhpca object
    return(Mhpca(
      hd_obj = hd_obj, 
      method = method, 
      ncomp = ncomp, 
      smooth_tuning = smooth_tuning, 
      sparse_tuning = sparse_tuning, 
      centerfns = centerfns, 
      alpha_orth = alpha_orth, 
      smoothing_type = smoothing_type, 
      sparse_type = sparse_type, 
      K_fold = K_fold, 
      sparse_CV = sparse_CV, 
      smooth_GCV = smooth_GCV,
      penalize_nfd = penalize_nfd, 
      penalize_u = penalize_u
    ))
  }
}

