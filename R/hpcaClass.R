#' @title A Class for `MHPCA` objects
#'
#' @description
#' The `mhpca` class represents hybrid principal components components.
#'

#' @field pc_mfd An object of class `mvmfd` where the first indices (fields)
#' represents harmonics and  second indices represents variables
#' @field pc_nfd An object of class `mvnfd` where the first indices (fields)
#' represents harmonics and  second indices represents variables
#' @field lsv = Left singular values vectors
#' @field values = The set of eigenvalues
#' @field smooth_tuning = The list of smoothing penalties parameters
#' @field sparse_tuning_u = The list of sparse penalties parameters
#' @field GCVs = Generalized cross validations scores of smoothing penalties parameters.
#'               If both smoothing and sparse tuning penalties are used in the MHPCA method,
#'               this represents the conditional generalized cross-validation scores, which
#'               means it is computed based on the optimal sparse tuning parameter selected via cross validation.
#' @field CVs_u = Cross validations scores of sparse penalties on u parameters
#' @field CVs_nfd = Cross validations scores of sparse penalties on nfd parameters
#' @field mean_mfd A multivariate functional data object giving the mean function
#' @field mean_nfd A data object giving the mean of non functional objects

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
#' hd_obj <- hd(mfd1, mfd2)
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
    initialize = function(hd_obj,
                          method = "power",
                          ncomp = 3,
                          smooth_tuning = NULL,
                          sparse_tuning_u = NULL,
                          sparse_tuning_nfd = NULL,
                          sparse_tuning_fd = NULL,
                          centerfns = TRUE,
                          alpha_orth = FALSE,
                          smoothing_type = "coefpen",
                          sparse_type_u = "soft",
                          sparse_type_nfd = "soft",
                          sparse_type_fd = "soft",
                          K_fold_u = 30,
                          K_fold_nfd = 30,
                          K_fold_fd = 30,
                          sparse_CV,
                          smooth_GCV,
                          penalize_nfd = FALSE,
                          penalize_fd = FALSE,
                          penalize_u = FALSE) {
      
      if (is.null(K_fold_nfd) && !is.null(K_fold_fd) && penalize_nfd == TRUE) {
        K_fold_nfd <- K_fold_fd
        warning("`K_fold_nfd` was not provided. It has been set to match `K_fold_fd`.")
      } else if (!is.null(K_fold_nfd) && is.null(K_fold_fd) && penalize_fd == TRUE) {
        K_fold_fd <- K_fold_nfd
        warning("`K_fold_fd` was not provided. It has been set to match `K_fold_nfd`.")
      } else if (!is.null(K_fold_nfd) && !is.null(K_fold_fd) && K_fold_fd != K_fold_nfd) {
        K_fold_fd <- K_fold_nfd
        warning("`K_fold_fd` and `K_fold_nfd` were not equal. `K_fold_fd` has been updated to match `K_fold_nfd`.")
      }
      
      if (is.mfd(hd_obj) || is.mvmfd(hd_obj) || is.nfd(hd_obj) || is.mvnfd(hd_obj)) {
        hd_obj <- Hd(hd_obj)
      }
      if (method == "power" & alpha_orth == "FALSE") {
        if (!is.null(hd_obj$mf)) {
          # Adjust the vector length to match the required dimensions if they are incorrect
          if (is.vector(smooth_tuning) & !is.list(smooth_tuning)) {
            if (smooth_GCV == FALSE) {
              if (length(smooth_tuning) != ncomp) {
                warning("The length of 'smooth_tuning' did not match 'ncomp' and has been adjusted accordingly.", call. = FALSE)
                smooth_tuning <- rep(smooth_tuning, length.out = ncomp)
              }
              smooth_tuning <- replicate(hd_obj$mf$nvar, smooth_tuning, simplify = FALSE)
            } else {
              warning("The length of 'smooth_tuning' did not match 'hd_obj$mf$nvar' and has been adjusted accordingly.", call. = FALSE)
              smooth_tuning <- replicate(hd_obj$mf$nvar, smooth_tuning, simplify = FALSE)
            }
          }

          # Adjust the matrix to match the required dimensions if they are incorrect
          else if (is.matrix(smooth_tuning)) {
            if (smooth_GCV == FALSE) {
              if (!all(dim(smooth_tuning) == c(hd_obj$mf$nvar, ncomp))) {
                smooth_tuning <- smooth_tuning[rep(1:nrow(smooth_tuning), length.out = hd_obj$mf$nvar), rep(1:ncol(smooth_tuning), length.out = ncomp)]
                # print(smooth_tuning)
                smooth_tuning <- split(smooth_tuning, row(smooth_tuning))
                # print(smooth_tuning)
                warning("The dimensions of 'smooth_tuning' did not match the expected size and have been adjusted accordingly.", call. = FALSE)
              } else {
                smooth_tuning <- split(smooth_tuning, row(smooth_tuning))
              }
            } else {
              if (dim(smooth_tuning)[1] != hd_obj$mf$nvar) {
                smooth_tuning <- smooth_tuning[rep(1:nrow(smooth_tuning), length.out = hd_obj$mf$nvar), , drop = FALSE][1:hd_obj$mf$nvar, , drop = FALSE]
                smooth_tuning <- split(smooth_tuning, row(smooth_tuning))
                warning("The dimensions of 'smooth_tuning' did not match the expected size and have been adjusted accordingly.", call. = FALSE)
              } else {
                smooth_tuning <- split(smooth_tuning, row(smooth_tuning))
              }
            }
          }

          # Adjust the list length and element sizes to match the required dimensions if they are incorrect
          else if (is.list(smooth_tuning)) {
            if (smooth_GCV == FALSE) {
              if (length(smooth_tuning) != hd_obj$mf$nvar) {
                warning("Adjusting 'smooth_tuning' to match 'hd_obj$mf$nvar'.", call. = FALSE)
                smooth_tuning <- rep(smooth_tuning, length.out = hd_obj$mf$nvar)
              }
              smooth_tuning <- lapply(smooth_tuning, function(vec) {
                if (length(vec) != ncomp) {
                  warning("Adjusting vector length in 'smooth_tuning' to match 'ncomp'.", call. = FALSE)
                  vec <- rep(vec, length.out = ncomp)
                }
                vec
              })
            } else {
              if (length(smooth_tuning) != hd_obj$mf$nvar) {
                warning("Adjusting 'smooth_tuning' to match 'hd_obj$mf$nvar'.", call. = FALSE)
                smooth_tuning <- rep(smooth_tuning, length.out = hd_obj$mf$nvar)
              }

              # smooth_tuning <- rep(smooth_tuning, length.out = hd_obj$mf$nvar)[1:hd_obj$mf$nvar]
            }
          }

          if (!is.null(smooth_tuning)) {
            names(smooth_tuning) <- paste0("var", 1:hd_obj$mf$nvar)
          }
          
          if (is.vector(sparse_tuning_fd) & !is.list(sparse_tuning_fd)) {
            if (sparse_CV == FALSE) {
              if (length(sparse_tuning_fd) != ncomp) {
                warning("The length of 'sparse_tuning_fd' did not match 'ncomp' and has been adjusted accordingly.", call. = FALSE)
                sparse_tuning_fd <- rep(sparse_tuning_fd, length.out = ncomp)
              }
              sparse_tuning_fd <- replicate(hd_obj$mf$nvar, sparse_tuning_fd, simplify = FALSE)
            } else {
              
              warning("The length of 'sparse_tuning_fd' did not match 'hd_obj$mf$nvar' and has been adjusted accordingly.", call. = FALSE)
              sparse_tuning_fd <- replicate(hd_obj$mf$nvar, sparse_tuning_fd, simplify = FALSE)
            }
          } 
          else if (is.list(sparse_tuning_fd)) {
            if (sparse_CV == FALSE) {
              if (length(sparse_tuning_fd) != hd_obj$mf$nvar) {
                warning("Adjusting 'sparse_tuning_fd' to match 'hd_obj$mf$nvar'.", call. = FALSE)
                sparse_tuning_fd <- rep(sparse_tuning_fd, length.out = hd_obj$mf$nvar)
              }
              sparse_tuning_fd <- lapply(sparse_tuning_fd, function(vec) {
                if (length(vec) != ncomp) {
                  warning("Adjusting vector length in 'sparse_tuning_fd' to match 'ncomp'.", call. = FALSE)
                  vec <- rep(vec, length.out = ncomp)
                }
                vec
              })
            } else {
              if (length(sparse_tuning_fd) != hd_obj$mf$nvar) {
                warning("Adjusting 'sparse_tuning_fd' to match 'hd_obj$mf$nvar'.", call. = FALSE)
                sparse_tuning_fd <- rep(sparse_tuning_fd, length.out = hd_obj$mf$nvar)
              }
              
              
            }
          }
          if (!is.null(sparse_tuning_fd)) {
            names(sparse_tuning_fd) <- paste0("var", 1:hd_obj$mf$nvar)
          }
        }


        # Adjust the list length and element sizes to match the required dimensions if they are incorrect
        if (sparse_CV == FALSE & length(sparse_tuning_u) != ncomp & !is.null(sparse_tuning_u)) {
          warning("The length of 'sparse_tuning_u' did not match 'ncomp' and has been adjusted accordingly.", call. = FALSE)
          sparse_tuning_u <- rep(sparse_tuning_u, length.out = ncomp)
        }
        if (sparse_CV == FALSE & length(sparse_tuning_nfd) != ncomp & !is.null(sparse_tuning_nfd)) {
          warning("The length of 'sparse_tuning_nfd' did not match 'ncomp' and has been adjusted accordingly.", call. = FALSE)
          sparse_tuning_nfd <- rep(sparse_tuning_nfd, length.out = ncomp)
        }
        # if (sparse_CV == FALSE & length(sparse_tuning_fd) != ncomp & !is.null(sparse_tuning_fd)) {
        #   warning("The length of 'sparse_tuning_nfd' did not match 'ncomp' and has been adjusted accordingly.", call. = FALSE)
        #   sparse_tuning_fd <- rep(sparse_tuning_fd, length.out = ncomp)
        # }
        #browser()
        result <- sequential_power_hybrid(
          hd_obj = hd_obj,
          n = ncomp,
          smooth_tuning = smooth_tuning,
          sparse_tuning_u = sparse_tuning_u,
          sparse_tuning_nfd = sparse_tuning_nfd,
          sparse_tuning_fd = sparse_tuning_fd,
          centerfns = centerfns,
          alpha_orth = alpha_orth,
          smooth_tuning_type = smoothing_type,
          sparse_tuning_type_u = sparse_type_u,
          sparse_tuning_type_nfd = sparse_type_nfd,
          sparse_tuning_type_fd = sparse_type_fd,
          K_fold_u = K_fold_u,
          K_fold_nfd = K_fold_nfd,
          K_fold_fd = K_fold_fd,
          sparse_CV = sparse_CV,
          smooth_GCV = smooth_GCV,
          penalize_nfd = penalize_nfd,
          penalize_fd = penalize_fd,
          penalize_u = penalize_u
        )
      } else if (method == "eigen" || alpha_orth == "TRUE") {
        if (!is.null(hd_obj$mf)) {
          # Adjust the vector to match the required lengths if they are incorrect
          if (is.vector(smooth_tuning) & !is.list(smooth_tuning)) {
            if (smooth_GCV == FALSE) {
              if (length(smooth_tuning) != hd_obj$mf$nvar) {
                warning("The length of 'smooth_tuning' did not match number of variables and has been adjusted accordingly.", call. = FALSE)
                smooth_tuning <- rep(smooth_tuning, length.out = hd_obj$mf$nvar)
              }
              smooth_tuning <- lapply(1:hd_obj$mf$nvar, function(i) smooth_tuning[i])
            } else {
              warning("The length of 'smooth_tuning' did not match number of variables and has been adjusted accordingly.", call. = FALSE)
              smooth_tuning <- replicate(hd_obj$mf$nvar, smooth_tuning, simplify = FALSE)
            }
          } else if (is.matrix(smooth_tuning)) { # Adjust the matrix to match the required if they are incorrect
            if (smooth_GCV == FALSE) {
              if (!all(dim(smooth_tuning) == c(hd_obj$mf$nvar, 1))) {
                smooth_tuning <- smooth_tuning[rep(1:nrow(smooth_tuning), length.out = hd_obj$mf$nvar), rep(1:ncol(smooth_tuning), length.out = 1)]
                smooth_tuning <- as.list(smooth_tuning)
                warning("The dimensions of 'smooth_tuning' did not match the expected size and have been adjusted accordingly.", call. = FALSE)
              } else {
                smooth_tuning <- as.list(smooth_tuning)
              }
            } else {
              if (dim(smooth_tuning)[1] != hd_obj$mf$nvar) {
                smooth_tuning <- smooth_tuning[rep(1:nrow(smooth_tuning), length.out = hd_obj$mf$nvar), , drop = FALSE][1:hd_obj$mf$nvar, , drop = FALSE]
                smooth_tuning <- split(smooth_tuning, row(smooth_tuning))
                warning("The dimensions of 'smooth_tuning' did not match the expected size and have been adjusted accordingly.", call. = FALSE)
              } else {
                smooth_tuning <- split(smooth_tuning, row(smooth_tuning))
              }
            }
          } else if (is.list(smooth_tuning)) { # Adjust the list length and element sizes to match the required dimensions if they are incorrect
            if (smooth_GCV == FALSE) {
              if (length(smooth_tuning) != hd_obj$mf$nvar) {
                warning("Adjusting 'smooth_tuning' to match 'hd_obj$mf$nvar'.", call. = FALSE)
                smooth_tuning <- rep(smooth_tuning, length.out = hd_obj$mf$nvar)
              }
              smooth_tuning <- lapply(smooth_tuning, function(vec) {
                if (length(vec) != 1) {
                  warning("Adjusting vector length in 'smooth_tuning' to match 'ncomp'.", call. = FALSE)
                  vec <- rep(vec, length.out = 1)
                }
                vec
              })
            } else {
              if (length(smooth_tuning) != hd_obj$mf$nvar) {
                warning("Adjusting 'smooth_tuning' to match 'hd_obj$mf$nvar'.", call. = FALSE)
                smooth_tuning <- rep(smooth_tuning, length.out = hd_obj$mf$nvar)[1:hd_obj$mf$nvar]
              }
            }
          }

          if (!is.null(smooth_tuning)) {
            names(smooth_tuning) <- paste0("var", 1:hd_obj$mf$nvar)
          }
        }

        if (method == "power") {
          result <- joint_power_hybrid(hd_obj = hd_obj, n = ncomp, smooth_tuning = smooth_tuning, centerfns = centerfns, alpha_orth = alpha_orth, smooth_tuning_type = smoothing_type)
        } else {
          result <- eigen_approach_hybrid(hd_obj = hd_obj, n = ncomp, alpha = smooth_tuning, centerfns = centerfns, penalty_type = smoothing_type)
        }

        # result <- joint_power(hd_obj = hd_obj, n = ncomp, smooth_tuning = smooth_tuning, centerfns = centerfns, alpha_orth = alpha_orth, smooth_tuning_type = smoothing_type)
      }
      # else if (method == "eigen") {
      #   result <- eigen_approach(hd_obj = hd_obj, n = ncomp, alpha = smooth_tuning, centerfns = centerfns, penalty_type = smoothing_type)
      # }

      coef <- result$pc_fd
      pcmfd <- list()
      if (!is.null(hd_obj$mf)) {
        for (i in 1:hd_obj$mf$nvar) {
          if (hd_obj$mf$basis$dimSupp[i] == 1) {
            pcmfd[[i]] <- Mfd(X = coef[[i]], mdbs = hd_obj$mf$basis$basis[[i]], method = "coefs")
          } else {
            coef_new <- array(coef[[i]], dim = c(hd_obj$mf$basis$basis[[i]]$nbasis, ncol(coef[[i]])))
            pcmfd[[i]] <- Mfd(X = coef_new, mdbs = hd_obj$mf$basis$basis[[i]], method = "coefs")
          }
        }
      } else {
        pcmfd <- NULL
      }

      if (!is.null(hd_obj$nf)) {
        pc_nfd <- result$pc_nfd
        pc_nfd_len <- sapply(hd_obj$nf$features, length)
        cumulative_indices <- c(0, cumsum(pc_nfd_len))
        pcnfd <- lapply(seq_along(pc_nfd_len), function(i) {
          start_row <- cumulative_indices[i] + 1
          end_row <- cumulative_indices[i + 1]
          t(nfd(pc_nfd[start_row:end_row, 1:ncomp, drop = FALSE]))
        })
      } else {
        pcnfd <- NULL
      }

      if (!is.null(hd_obj$mf) && !is.null(hd_obj$nf)) {
        mfd_obj <- Mvmfd(pcmfd)
        nfd_obj <- Mvnfd(pcnfd)
        hd_obj <- Hd(mfd_obj, nfd_obj)
      } else if (!is.null(hd_obj$mf) && is.null(hd_obj$nf)) {
        mfd_obj <- Mvmfd(pcmfd)
        nfd_obj <- NULL
        hd_obj <- mfd_obj # Hd(mfd_obj)
      } else if (is.null(hd_obj$mf) && !is.null(hd_obj$nf)) {
        mfd_obj <- NULL
        nfd_obj <- Mvnfd(pcnfd)
        hd_obj <- nfd_obj # Hd(nfd_obj)
      }

      private$.pc_hd <- hd_obj
      # private$.pc_mfd <- mfd_obj
      # private$.pc_nfd <- pcnfd #nfd_obj

      private$.lsv <- result$lsv
      private$.values <- result$variance
      private$.smooth_tuning <- result$smooth_tuning
      if (alpha_orth == "FALSE" && method == "power") {
        private$.sparse_tuning_u <- result$sparse_tuning_result_u
        private$.sparse_tuning_nfd <- result$sparse_tuning_result_nfd
        private$.sparse_tuning_fd <- result$sparse_tuning_result_fd
        private$.CVs_u <- result$CV_score_u
        private$.CVs_nfd <- result$CV_score_nfd
        private$.CVs_fd <- result$CV_score_fd
        private$.GCVs <- result$GCV_score
      } else {
        private$.GCVs <- result$GCV_score
      }
      # private$.sparse_tuning <- result[[5]]
      # private$.mean_mfd <- mean(mfd_obj)
      # private$.mean_nfd <- mean(mfd_obj) # change it later
      private$.mean_hd <- mean(hd_obj)
    }
  ),
  active = list(
    pc_hd = function(value) {
      if (missing(value)) {
        private$.pc_hd
      } else {
        stop("`$pc_hd` is read only", call. = FALSE)
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
    sparse_tuning_u = function(value) {
      if (missing(value)) {
        private$.sparse_tuning_u
      } else {
        stop("`$sparse_tuning_u` is read only", call. = FALSE)
      }
    },
    sparse_tuning_nfd = function(value) {
      if (missing(value)) {
        private$.sparse_tuning_nfd
      } else {
        stop("`$sparse_tuning_nfd` is read only", call. = FALSE)
      }
    },
    sparse_tuning_fd = function(value) {
      if (missing(value)) {
        private$.sparse_tuning_fd
      } else {
        stop("`$sparse_tuning_fd` is read only", call. = FALSE)
      }
    },
    GCVs = function(value) {
      if (missing(value)) {
        private$.GCVs
      } else {
        stop("`$GCVs` is read only", call. = FALSE)
      }
    },
    CVs_u = function(value) {
      if (missing(value)) {
        private$.CVs_u
      } else {
        stop("`$CVs_u` is read only", call. = FALSE)
      }
    },
    CVs_nfd = function(value) {
      if (missing(value)) {
        private$.CVs_nfd
      } else {
        stop("`$CVs_nfd` is read only", call. = FALSE)
      }
    },
    CVs_fd = function(value) {
      if (missing(value)) {
        private$.CVs_fd
      } else {
        stop("`$CVs_fd` is read only", call. = FALSE)
      }
    },
    mean_hd = function(value) {
      if (missing(value)) {
        private$.mean_hd
      } else {
        stop("`$mean_hd` is read only", call. = FALSE)
      }
    }
  ),
  private = list(
    .pc_hd = NULL,
    .lsv = NULL,
    .values = NULL,
    .smooth_tuning = NULL,
    .sparse_tuning_u = NULL,
    .sparse_tuning_nfd = NULL,
    .sparse_tuning_fd = NULL,
    .GCVs = NULL,
    .CVs_u = NULL,
    .CVs_nfd = NULL,
    .CVs_fd = NULL,
    .mean_hd = NULL
  )
)

#' @rdname mhpca
#' @seealso \code{\link{hd}}

#' @title A Class for 'mhpca' objects
#'
#' @description
#' The `mhpca` class represents regularized functional principal components ('MHFPCs') components.
#'
#' @param hd_obj An `hd` object representing the multivariate functional data.
#' @param method A character string specifying the approach to be used for MFPCA computation.
#'               Options are "power" (the default), which uses the power algorithm, or "eigen",
#'               which uses the eigen decomposition approach.
#' @param ncomp The number of functional principal components to retain.
#' @param smooth_tuning A list or vector specifying the smoothing regularization parameter(s) for each variable.
#'                      If NULL, non-smoothing MFPCA is estimated.
#' @param sparse_tuning_u A list or vector specifying the sparsity regularization parameter(s) for each variable.
#'                      If NULL, non-sparse MHPCA is estimated.
#' @param sparse_tuning_nfd A list or vector specifying the sparsity regularization parameter(s) for each non functional variable.
#'                      If NULL, non-sparse MHPCA is estimated.
#' @param sparse_tuning_fd A list or vector specifying the sparsity regularization parameter(s) for each functional variable.
#'                      If NULL, non-sparse MHPCA is estimated.
#'
#' @param centerfns Logical indicating whether to center the functional data before analysis. Default is TRUE.
#' @param alpha_orth Logical indicating whether to perform orthogonalization of the regularization parameters.
#'                   If `method` is "power", setting `alpha_orth = FALSE` (default) uses the sequential power approach,
#'                   while setting `alpha_orth = TRUE` uses the joint power approach.
#' @param smoothing_type The type of smoothing penalty to be applied on the coefficients. The types "coefpen" and "basispen" is supported. Default is "coefpen".
#' @param sparse_type_u The type of sparse penalty to be applied on the coefficients. The types "soft", "hard" and "SCAD" is supported. Default is "soft".
#' @param sparse_type_nfd The type of sparse penalty to be applied on the nfd right singular vectors. The types "soft", "hard" and "SCAD" is supported. Default is "soft".
#' @param sparse_type_fd The type of sparse penalty to be applied on the fd right singular vectors. The types "soft", "hard" and "SCAD" is supported. Default is "soft".
#' @param K_fold_u  An integer specifying the number of folds in the sparse cross-validation process for u. Default is 30.
#' @param K_fold_nfd  An integer specifying the number of folds in the sparse cross-validation process for nfd. Default is 30.
#' @param K_fold_fd  An integer specifying the number of folds in the sparse cross-validation process for fd. Default is 30.
#' @param sparse_CV Logical indicating whether cross-validation should be applied to select the optimal sparse tuning parameter in sequential power approach.
#'                                        If `sparse_CV = TRUE`, a series of tuning parameters should be provided as a vector with positive number with max equals to number of subjects.
#'                                        If `sparse_CV = FALSE`, specific tuning parameters are given directly to each principal components. Tuning parameters should be provided as a vector with length equal to `ncomp`.
#'                                        If the dimensions of input tuning parameters are incorrect, it will be converted to a list internally, and a warning will be issued.
#' @param smooth_GCV Logical indicating whether generalized cross-validation should be applied to select the optimal smooth tuning parameter.
#'                                        If `smooth_GCV = TRUE`, a series of tuning parameters should be provided as a list with length equal to the number of variables.
#'                                        If a list with incorrect dimensions is provided, it will be converted to a correct list internally, and a warning will be issued.
#'                                        If `smooth_GCV = FALSE`, specific tuning parameters are given directly. If `method` is "power" and `alpha_orth = FALSE` (sequential power),
#'                                        tuning parameters should be provided as a list with length equal to the number of variables, where each element is a vector of length `ncomp`.
#'                                        If `method` is "power" and `alpha_orth = TRUE` (joint power), tuning parameters should be provided as a vector with length equal to the number of variables.
#'                                        If the dimensions of input tuning parameters are incorrect, it will be converted to a list internally, and a warning will be issued.
#' @param penalize_nfd Logical indicating whether sparsity penalty in sequential power approach should be applied on nfd right singular vector.
#' @param penalize_fd Logical indicating whether sparsity penalty in sequential power approach should be applied on fd right singular vector.
#' @param penalize_u Logical indicating whether penalize non functional object or left singular vector or not.
#' @export
Mhpca <- function(hd_obj,
                  method = "power",
                  ncomp = 3,
                  smooth_tuning = NULL,
                  sparse_tuning_u = NULL,
                  sparse_tuning_nfd = NULL,
                  sparse_tuning_fd = NULL,
                  centerfns = TRUE,
                  alpha_orth = FALSE,
                  smoothing_type = "basispen",
                  sparse_type_u = "soft",
                  sparse_type_nfd = "soft",
                  sparse_type_fd = "soft",
                  K_fold_u = 30,
                  K_fold_nfd = 30,
                  K_fold_fd = 30,
                  sparse_CV = TRUE,
                  smooth_GCV = TRUE,
                  penalize_nfd = FALSE,
                  penalize_fd = FALSE,
                  penalize_u = FALSE) {
  mhpca$new(
    hd_obj = hd_obj,
    method = method,
    ncomp = ncomp,
    smooth_tuning = smooth_tuning,
    sparse_tuning_u = sparse_tuning_u,
    sparse_tuning_nfd = sparse_tuning_nfd,
    sparse_tuning_fd = sparse_tuning_fd,
    centerfns = centerfns,
    alpha_orth = alpha_orth,
    smoothing_type = smoothing_type,
    sparse_type_u = sparse_type_u,
    sparse_type_nfd = sparse_type_nfd,
    sparse_type_fd = sparse_type_fd,
    K_fold_u = K_fold_u,
    K_fold_nfd = K_fold_nfd,
    K_fold_fd = K_fold_fd,
    sparse_CV = sparse_CV,
    smooth_GCV = smooth_GCV,
    penalize_nfd = penalize_nfd,
    penalize_fd = penalize_fd,
    penalize_u = penalize_u
  )
}
#' @rdname mhpca
