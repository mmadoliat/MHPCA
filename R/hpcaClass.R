# Output fields of an mhpca object. Each is stored in private$.<name> and
# exposed as a read-only active binding; to add an output, add its name here
# and assign private$.<name> in initialize().
mhpca_fields <- c(
  "pc_hd", "lsv", "values", "smooth_tuning", "sparse_tuning_u",
  "sparse_tuning_nfd", "sparse_tuning_fd", "GCVs", "CVs_u", "CVs_nfd",
  "CVs_fd", "eBICs_u", "eBICs_nfd", "eBICs_fd", "selection_method", "mean_hd"
)

# Read-only active bindings for R6 fields stored as private$.<name>.
# `labels` overrides the field name shown in the error message.
# The bodies are built with the field name written in literally because R6
# re-parents active-binding functions, so they cannot rely on closures.
read_only_fields <- function(fields, labels = character()) {
  bindings <- lapply(fields, function(field) {
    label <- if (field %in% names(labels)) labels[[field]] else field
    f <- function(value) NULL
    body(f) <- substitute(
      if (missing(value)) private[[PRIVATE]] else stop(MESSAGE, call. = FALSE),
      list(PRIVATE = paste0(".", field), MESSAGE = sprintf("`$%s` is read only", label))
    )
    f
  })
  names(bindings) <- fields
  bindings
}

# Private storage (all NULL initially) for the fields above.
private_fields <- function(fields) {
  stats::setNames(rep(list(NULL), length(fields)), paste0(".", fields))
}

# Call `target` with every argument of the calling function, passed by name
# and still unevaluated: forward_args(f) inside g(a, b) is the same as writing
# f(a = a, b = b). Used so a wrapper such as Mhpca() does not have to repeat
# its whole argument list.
forward_args <- function(target) {
  caller <- parent.frame()
  arg_names <- names(formals(sys.function(sys.parent())))
  args <- stats::setNames(lapply(arg_names, as.name), arg_names)
  eval(as.call(c(substitute(target), args)), caller)
}

#' @title A Class for `MHPCA` objects
#'
#' @description
#' The `mhpca` class represents multivariate hybrid principal components.
#'
#' @field pc_hd The principal components: an `hd` object for hybrid data, or
#'   an `mvmfd` / `mvnfd` object when the data has only a functional or only a
#'   non-functional part. The first index selects components.
#' @field lsv Left singular vectors (scores direction), one column per component.
#' @field values The variance explained by each component.
#' @field smooth_tuning The selected smoothing parameters.
#' @field sparse_tuning_u The selected sparsity parameters for the scores (u).
#' @field sparse_tuning_nfd The selected sparsity parameters for the non-functional loadings.
#' @field sparse_tuning_fd The selected sparsity parameters for the functional loadings.
#' @field GCVs Generalized cross-validation scores of the smoothing parameters.
#'   If sparsity is also tuned, these are computed at the selected sparsity.
#' @field CVs_u Cross-validation scores of the sparsity parameters for u.
#' @field CVs_nfd Cross-validation scores of the sparsity parameters for the non-functional loadings.
#' @field CVs_fd Cross-validation scores of the sparsity parameters for the functional loadings.
#' @field eBICs_u Extended BIC scores of the sparsity parameters for u.
#' @field eBICs_nfd Extended BIC scores of the sparsity parameters for the non-functional loadings.
#' @field eBICs_fd Extended BIC scores of the sparsity parameters for the functional loadings.
#' @field selection_method How sparsity parameters were selected ("cv" or "ebic").
#' @field mean_hd The mean of the data, with the same structure as `pc_hd`.
#'
#' @examples
#' require(fda)
#' # Brownian Bridge simulation on [0,1]
#' M <- 110 # number of components
#' N <- 20 # number of instances
#' n <- 100 # number of grid points
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
#' hd_obj <- Hd(mfd1, mfd2)
#' k <- 2
#' # Smoothing parameters selected by GCV from a grid for each variable
#' Re0 <- Mhpca(hd_obj, ncomp = k, smooth_tuning = alpha1)
#' Re0$smooth_tuning
#' scores0 <- inprod_mvmfd(hd_obj$mf, Re0$pc_hd)
#' dim(scores0)
#' # Fixed smoothing parameters, one per component
#' Re1 <- Mhpca(hd_obj, ncomp = k, smooth_tuning = c(0, 0), smooth_GCV = FALSE)
#' Re1$values
#' # A single functional variable
#' Re2 <- Mhpca(mfd1, ncomp = k, smooth_tuning = list(alpha1$a1))
#' Re2$smooth_tuning

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
                          nfold_u = 30,
                          nfold_nfd = 30,
                          nfold_fd = 30,
                          parallel = FALSE,#ncor = 1,
                          sparse_iter = 1,
                          tol = 1e-4, 
                          max_iter = 1000,
                          sparse_CV,
                          cv.pick = "1se",
                          smooth_GCV,
                          pen_nfd = FALSE,
                          pen_fd = FALSE,
                          pen_u = FALSE,
                          sparse_selection_method = "cv",
                          ebic_xi = 0.5,
                          ebic_N_type = "product",
                          ebic_RSS_type = "direct") {
      
      # Validate sparse_selection_method
      sparse_selection_method <- match.arg(sparse_selection_method, c("cv", "ebic"))
      method <- check_choice(method, c("power", "eigen"), "method")
      smoothing_type <- check_choice(smoothing_type, c("coefpen", "basispen"), "smoothing_type")
      sparse_type_u <- check_choice(sparse_type_u, sparse_penalty_types, "sparse_type_u")
      sparse_type_nfd <- check_choice(sparse_type_nfd, sparse_penalty_types, "sparse_type_nfd")
      sparse_type_fd <- check_choice(sparse_type_fd, sparse_penalty_types, "sparse_type_fd")
      cv.pick <- check_choice(cv.pick, c("1se", "min"), "cv.pick")
      ebic_RSS_type <- check_choice(ebic_RSS_type, c("direct", "refitted"), "ebic_RSS_type")
      alpha_orth <- check_flag(alpha_orth, "alpha_orth")
      # sequential power (one component at a time, with sparsity) vs joint power / eigen
      sequential <- method == "power" && !alpha_orth
      
      # Validate ebic_xi
      if (sparse_selection_method == "ebic") {
        if (!is.numeric(ebic_xi) || length(ebic_xi) != 1 || ebic_xi < 0 || ebic_xi > 1) {
          stop("ebic_xi must be a single numeric value between 0 and 1", call. = FALSE)
        }
      }
      
      # Validate ebic_N_type
      if (sparse_selection_method == "ebic") {
        if (!is.numeric(ebic_N_type)) {
          ebic_N_type <- match.arg(ebic_N_type, c("product", "sum"))
        }
      }
      
      if (is.mfd(hd_obj) || is.mvmfd(hd_obj) || is.nfd(hd_obj) || is.mvnfd(hd_obj)) {
        hd_obj <- Hd(hd_obj)
      }
      
      if (is.null(smooth_tuning) && !is.null(hd_obj$mf)) {
        smooth_tuning <- hd_obj$mf$smooth_par
      }
      
      if (is.null(sparse_tuning_fd) && !is.null(hd_obj$mf)) {
        sparse_tuning_fd <- hd_obj$mf$spars_par
      }
      
      if (is.null(sparse_tuning_nfd) && !is.null(hd_obj$nf)) {
        sparse_tuning_nfd <- hd_obj$nf$spars_par
      }
      if (is.null(sparse_tuning_u)) {
        sparse_tuning_u <- hd_obj$spars_par
      }
      # Sparsity grids for a data part that is not present are ignored (tuning
      # them would run on empty data).
      if (is.null(hd_obj$mf)) sparse_tuning_fd <- NULL
      if (is.null(hd_obj$nf)) sparse_tuning_nfd <- NULL
      # A block is only penalized when it has sparsity parameters.
      pen_u <- isTRUE(pen_u) && !is.null(sparse_tuning_u)
      pen_nfd <- isTRUE(pen_nfd) && !is.null(sparse_tuning_nfd)
      pen_fd <- isTRUE(pen_fd) && !is.null(sparse_tuning_fd)
      if (isTRUE(parallel)) {
        ## use workers that the user already registered
        ncor <- foreach::getDoParWorkers()
        if (ncor < 2L)                         # user forgot to register?
          ncor <- max(1L, parallel::detectCores() - 1L)
      } else {                                    # sequential path
        ncor <- 1L
      }
      
      if (sequential) {
        if (!is.null(hd_obj$mf)) {
          smooth_tuning <- normalize_tuning(
            smooth_tuning, "smooth_tuning", nvar = hd_obj$mf$nvar, ncomp = ncomp,
            grid = smooth_GCV, mode = "per_component",
            vec_label = "'hd_obj$mf$nvar'", list_label = "'hd_obj$mf$nvar'",
            allow_matrix = TRUE
          )
          sparse_tuning_fd <- normalize_tuning(
            sparse_tuning_fd, "sparse_tuning_fd", nvar = hd_obj$mf$nvar, ncomp = ncomp,
            grid = sparse_CV, mode = "per_component",
            vec_label = "'hd_obj$mf$nvar'", list_label = "'hd_obj$mf$nvar'"
          )
        }

        if (!is.null(hd_obj$nf)) {
          sparse_tuning_nfd <- normalize_tuning(
            sparse_tuning_nfd, "sparse_tuning_nfd", nvar = hd_obj$nf$nvar, ncomp = ncomp,
            grid = sparse_CV, mode = "per_component",
            vec_label = "'hd_obj$nf$nvar'", list_label = "'hd_obj$nf$nvar'"
          )
        }

        sparse_tuning_u <- normalize_tuning_u(sparse_tuning_u, "sparse_tuning_u", ncomp, sparse_CV)

        settings <- mhpca_settings(
          smoothing_type = smoothing_type, centerfns = centerfns,
          smooth_GCV = smooth_GCV, sparse_CV = sparse_CV,
          sparse_type_u = sparse_type_u, sparse_type_nfd = sparse_type_nfd,
          sparse_type_fd = sparse_type_fd,
          nfold_u = nfold_u, nfold_nfd = nfold_nfd, nfold_fd = nfold_fd,
          pen_u = pen_u, pen_nfd = pen_nfd, pen_fd = pen_fd,
          sparse_iter = sparse_iter, cv.pick = cv.pick,
          sparse_selection_method = sparse_selection_method,
          ebic_xi = ebic_xi, ebic_N_type = ebic_N_type, ebic_RSS_type = ebic_RSS_type,
          tol = tol, max_iter = max_iter, ncor = ncor
        )
        result <- sequential_power_hybrid(
          hd_obj = hd_obj,
          n = ncomp,
          smooth_tuning = smooth_tuning,
          sparse_tuning_u = sparse_tuning_u,
          sparse_tuning_nfd = sparse_tuning_nfd,
          sparse_tuning_fd = sparse_tuning_fd,
          settings = settings
        )
      } else {
        if (!is.null(hd_obj$mf)) {
          smooth_tuning <- normalize_tuning(
            smooth_tuning, "smooth_tuning", nvar = hd_obj$mf$nvar, ncomp = ncomp,
            grid = smooth_GCV, mode = "per_variable",
            vec_label = "number of variables", list_label = "'hd_obj$mf$nvar'",
            allow_matrix = TRUE
          )
        }

        if (method == "power") {
          settings <- mhpca_settings(smoothing_type = smoothing_type, centerfns = centerfns,
                                     tol = tol, max_iter = max_iter)
          result <- joint_power_hybrid(hd_obj = hd_obj, n = ncomp, smooth_tuning = smooth_tuning,
                                       settings = settings)
        } else {
          result <- eigen_approach_hybrid(hd_obj = hd_obj, n = ncomp, alpha = smooth_tuning, centerfns = centerfns, penalty_type = smoothing_type)
        }

      }

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

      # Mean of the input data (hd_obj is replaced by the components below).
      data_mean <- mean(hd_obj)

      if (!is.null(hd_obj$mf) && !is.null(hd_obj$nf)) {
        mfd_obj <- Mvmfd(pcmfd)
        nfd_obj <- Mvnfd(pcnfd)
        hd_obj <- Hd(mfd_obj, nfd_obj)
      } else if (!is.null(hd_obj$mf) && is.null(hd_obj$nf)) {
        mfd_obj <- Mvmfd(pcmfd)
        nfd_obj <- NULL
        hd_obj <- mfd_obj 
      } else if (is.null(hd_obj$mf) && !is.null(hd_obj$nf)) {
        mfd_obj <- NULL
        nfd_obj <- Mvnfd(pcnfd)
        hd_obj <- nfd_obj 
      }

      private$.pc_hd <- hd_obj

      private$.lsv <- result$lsv
      private$.values <- result$variance
      private$.smooth_tuning <- result[["smooth_tuning_result"]]
      if (sequential) {
        
        private$.sparse_tuning_u <- result$sparse_tuning_result_u
        private$.sparse_tuning_nfd <- result$sparse_tuning_result_nfd
        private$.sparse_tuning_fd <- result$sparse_tuning_result_fd
        private$.selection_method <- sparse_selection_method
        private$.GCVs <- result[["GCV_score"]]
        # Store CV or eBIC scores based on selection method
        if (sparse_selection_method == "cv") {
          private$.CVs_u <- result$CV_score_u
          private$.CVs_nfd <- result$CV_score_nfd
          private$.CVs_fd <- result$CV_score_fd
          private$.eBICs_u <- NULL
          private$.eBICs_nfd <- NULL
          private$.eBICs_fd <- NULL
        } else {
          private$.CVs_u <- NULL
          private$.CVs_nfd <- NULL
          private$.CVs_fd <- NULL
          private$.eBICs_u <- result$eBIC_score_u
          private$.eBICs_nfd <- result$eBIC_score_nfd
          private$.eBICs_fd <- result$eBIC_score_fd
        }
        
        # private$.CV_se_u <- result$cv_se_u
        # private$.CV_se_nfd <- result$cv_se_nfd
        # private$.CV_se_fd <- result$cv_se_fd
        
      } else {
        private$.GCVs <- result[["GCV_score"]]
      }
      # Same form as pc_hd: hd for hybrid data, otherwise the single part.
      private$.mean_hd <- if (is.hd(hd_obj)) data_mean else if (is.mvmfd(hd_obj)) data_mean$mf else data_mean$nf
    }
  ),
  active = read_only_fields(mhpca_fields),
  private = private_fields(mhpca_fields)
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
#' @param sparse_type_u The type of sparse penalty to be applied on the coefficients. The types "soft", "hard", "SCAD" and "adaptive" are supported. Default is "soft".
#' @param sparse_type_nfd The type of sparse penalty to be applied on the nfd right singular vectors. The types "soft", "hard", "SCAD" and "adaptive" are supported. Default is "soft".
#' @param sparse_type_fd The type of sparse penalty to be applied on the fd right singular vectors. The types "soft", "hard", "SCAD" and "adaptive" are supported. Default is "soft".
#' @param nfold_u  An integer specifying the number of folds in the sparse cross-validation process for u. Default is 30.
#' @param nfold_nfd  An integer specifying the number of folds in the sparse cross-validation process for nfd. Default is 30.
#' @param nfold_fd  An integer specifying the number of folds in the sparse cross-validation process for fd. Default is 30.
#' @param parallel `TRUE`/`FALSE` parallel computing of Cross Validation. 
#' @param sparse_iter number of iteration for sparse parameter selection 
#' @param tol power algorithm tolerance 
#' @param max_iter power algorithm max iteration 
#' @param sparse_CV Logical indicating whether cross-validation should be applied to select the optimal sparse tuning parameter in sequential power approach.
#'                                        If `sparse_CV = TRUE`, a series of tuning parameters should be provided as a vector with positive number with max equals to number of subjects.
#'                                        If `sparse_CV = FALSE`, specific tuning parameters are given directly to each principal components. Tuning parameters should be provided as a vector with length equal to `ncomp`.
#'                                        If the dimensions of input tuning parameters are incorrect, it will be converted to a list internally, and a warning will be issued.
#' @param cv.pick CV selection either based on "1se" or "min"
#' @param smooth_GCV Logical indicating whether generalized cross-validation should be applied to select the optimal smooth tuning parameter.
#'                                        If `smooth_GCV = TRUE`, a series of tuning parameters should be provided as a list with length equal to the number of variables.
#'                                        If a list with incorrect dimensions is provided, it will be converted to a correct list internally, and a warning will be issued.
#'                                        If `smooth_GCV = FALSE`, specific tuning parameters are given directly. If `method` is "power" and `alpha_orth = FALSE` (sequential power),
#'                                        tuning parameters should be provided as a list with length equal to the number of variables, where each element is a vector of length `ncomp`.
#'                                        If `method` is "power" and `alpha_orth = TRUE` (joint power), tuning parameters should be provided as a vector with length equal to the number of variables.
#'                                        If the dimensions of input tuning parameters are incorrect, it will be converted to a list internally, and a warning will be issued.
#' @param pen_nfd Logical indicating whether sparsity penalty in sequential power approach should be applied on nfd right singular vector.
#' @param pen_fd Logical indicating whether sparsity penalty in sequential power approach should be applied on fd right singular vector.
#' @param pen_u Logical indicating whether penalize non functional object or left singular vector or not.
#' @param sparse_selection_method Character: "cv" for cross-validation or "ebic" for extended BIC. Default is "cv".
#' @param ebic_xi eBIC tuning parameter in [0,1]. 0 = standard BIC, 1 = most conservative. Default is 0.5.
#' @param ebic_N_type How to compute effective sample size for eBIC: "product", "sum", or a numeric value.
#' @param ebic_RSS_type Residual sum of squares used by eBIC: "direct" (default) or "refitted".
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
                  nfold_u = 30,
                  nfold_nfd = 30,
                  nfold_fd = 30,
                  parallel = FALSE,
                  sparse_iter = 1,
                  tol = 1e-4, 
                  max_iter = 1000,
                  sparse_CV = TRUE,
                  cv.pick = "1se",
                  smooth_GCV = TRUE,
                  pen_nfd = FALSE,
                  pen_fd = FALSE,
                  pen_u = FALSE,
                  sparse_selection_method = "cv",
                  ebic_xi = 0.5,
                  ebic_N_type = "product",
                  ebic_RSS_type = "direct") {
  forward_args(mhpca$new)
}
#' @rdname mhpca
