# Small synthetic data and a catalogue of Mhpca() calls covering the option
# combinations the yeast tests do not reach (eBIC, fixed tuning, cv.pick,
# sparse_iter, penalty types, single-part data, ...). Each fit runs from its
# own seed, so the cases are independent of each other.

paths_data <- function() {
  if (!is.null(yeast_cache$paths_data)) return(yeast_cache$paths_data)
  set.seed(42)
  n <- 30
  grid <- seq(0, 1, length.out = 25)
  scores <- matrix(stats::rnorm(n * 2), n, 2)
  curves <- outer(sin(2 * pi * grid), scores[, 1]) + outer(cos(2 * pi * grid), scores[, 2]) +
    matrix(stats::rnorm(length(grid) * n, sd = 0.1), length(grid), n)
  vectors <- cbind(scores %*% matrix(stats::rnorm(8), 2, 4), matrix(stats::rnorm(n * 2), n, 2))
  colnames(vectors) <- paste0("f", 1:6)
  mf <- Mfd(grid, curves, Basismfd(fda::create.bspline.basis(c(0, 1), 8)))
  nf <- nfd(vectors)
  yeast_cache$paths_data <- list(hd = Hd(mf, nf), mf_only = Hd(mf), nf_only = Hd(nf))
  yeast_cache$paths_data
}

# Arguments shared by the sparse cases: J = 8 basis functions and m = 6
# features are both divisible by nfold_u = 2; n = 30 rows by nfold_nfd/fd = 5.
sparse_args <- list(
  ncomp = 2, smooth_tuning = list(2^c(-8, -2, 2)),
  sparse_tuning_u = c(0L, 5L, 10L), sparse_tuning_nfd = c(0L, 2L, 4L),
  sparse_tuning_fd = c(0L, 2L, 4L), nfold_u = 2, nfold_nfd = 5, nfold_fd = 5,
  pen_u = TRUE, pen_nfd = TRUE, pen_fd = TRUE
)

paths_cases <- function() {
  d <- paths_data()
  with_args <- function(...) utils::modifyList(sparse_args, list(...))
  list(
    cv_1se = c(list(hd_obj = d$hd), sparse_args),
    cv_min = c(list(hd_obj = d$hd), with_args(cv.pick = "min")),
    cv_iter2 = c(list(hd_obj = d$hd), with_args(sparse_iter = 2)),
    cv_types = c(list(hd_obj = d$hd), with_args(sparse_type_u = "hard",
                                                 sparse_type_nfd = "SCAD",
                                                 sparse_type_fd = "adaptive")),
    cv_coefpen_uncentered = c(list(hd_obj = d$hd), with_args(smoothing_type = "coefpen",
                                                            centerfns = FALSE)),
    ebic_direct = c(list(hd_obj = d$hd), with_args(sparse_selection_method = "ebic")),
    ebic_refitted = c(list(hd_obj = d$hd), with_args(sparse_selection_method = "ebic",
                                                     ebic_RSS_type = "refitted",
                                                     ebic_xi = 0, ebic_N_type = "sum")),
    ebic_numeric_N = c(list(hd_obj = d$hd), with_args(sparse_selection_method = "ebic",
                                                      ebic_N_type = 100)),
    # Only some blocks penalized (grids for the others are ignored).
    cv_pen_u_only = c(list(hd_obj = d$hd), with_args(pen_nfd = FALSE, pen_fd = FALSE)),
    cv_pen_nfd_only = c(list(hd_obj = d$hd), with_args(pen_u = FALSE, pen_fd = FALSE)),
    cv_pen_fd_only = c(list(hd_obj = d$hd), with_args(pen_u = FALSE, pen_nfd = FALSE)),
    ebic_pen_fd_only = c(list(hd_obj = d$hd), with_args(sparse_selection_method = "ebic",
                                                        pen_u = FALSE, pen_nfd = FALSE)),
    ebic_pen_u_nfd = c(list(hd_obj = d$hd), with_args(sparse_selection_method = "ebic",
                                                      pen_fd = FALSE)),
    # smooth_GCV = FALSE (smoothing given per PC) combined with other options.
    fixed_smooth_cv = c(list(hd_obj = d$hd), with_args(smooth_GCV = FALSE,
                                                       smooth_tuning = c(0.01, 0.1))),
    fixed_smooth_ebic = c(list(hd_obj = d$hd), with_args(smooth_GCV = FALSE,
                                                         smooth_tuning = c(0.01, 0.1),
                                                         sparse_selection_method = "ebic")),
    fixed_smooth_mf_only = list(hd_obj = d$mf_only, ncomp = 2, smooth_GCV = FALSE,
                                smooth_tuning = c(0.01, 0.1), sparse_tuning_u = c(0L, 5L),
                                sparse_tuning_fd = c(0L, 2L), nfold_u = 2, nfold_fd = 5,
                                pen_u = TRUE, pen_fd = TRUE),
    fixed_smooth_nf_only = list(hd_obj = d$nf_only, ncomp = 2, smooth_GCV = FALSE,
                                sparse_tuning_u = c(0L, 5L), sparse_tuning_nfd = c(0L, 2L),
                                nfold_u = 2, nfold_nfd = 5, pen_u = TRUE, pen_nfd = TRUE),
    fixed_smooth_null = list(hd_obj = d$hd, ncomp = 2, smooth_GCV = FALSE),
    grid_smooth_null = list(hd_obj = d$hd, ncomp = 2, smooth_GCV = TRUE),
    fixed_tuning = list(hd_obj = d$hd, ncomp = 2, sparse_CV = FALSE, smooth_GCV = FALSE,
                        smooth_tuning = c(0.01, 0.1), sparse_tuning_u = c(5L, 3L),
                        sparse_tuning_nfd = c(2L, 1L), sparse_tuning_fd = c(1L, 2L),
                        pen_u = TRUE, pen_nfd = TRUE, pen_fd = TRUE),
    mf_only_cv = list(hd_obj = d$mf_only, ncomp = 2, smooth_tuning = list(2^c(-8, 2)),
                      sparse_tuning_u = c(0L, 5L), sparse_tuning_fd = c(0L, 2L),
                      nfold_u = 2, nfold_fd = 5, pen_u = TRUE, pen_fd = TRUE),
    nf_only_cv = list(hd_obj = d$nf_only, ncomp = 2, sparse_tuning_u = c(0L, 5L),
                      sparse_tuning_nfd = c(0L, 2L), nfold_u = 2, nfold_nfd = 5,
                      pen_u = TRUE, pen_nfd = TRUE),
    joint = list(hd_obj = d$hd, ncomp = 2, alpha_orth = TRUE, smooth_tuning = list(2^c(-8, 2))),
    joint_mf_only = list(hd_obj = d$mf_only, ncomp = 2, alpha_orth = TRUE),
    joint_nf_only = list(hd_obj = d$nf_only, ncomp = 2, alpha_orth = TRUE),
    eigen_nf_only = list(hd_obj = d$nf_only, ncomp = 2, method = "eigen"),
    eigen_fixed = list(hd_obj = d$hd, ncomp = 2, method = "eigen", smooth_GCV = FALSE,
                       smooth_tuning = 0.1)
  )
}

# Fit one case from a fixed seed; returns the fit (or error) and warnings.
run_path_case <- function(args, seed = 1) {
  set.seed(seed)
  out <- tryCatch(quiet_fit(do.call(Mhpca, args)),
                  error = function(e) list(value = NULL, warnings = character(),
                                           error = conditionMessage(e)))
  out
}

# The numeric content of a fit, for snapshots.
fit_summary <- function(out) {
  fit <- out$value
  if (is.null(fit)) return(list(error = out$error))
  pc <- fit$pc_hd
  list(
    values = fit$values,
    lsv = fit$lsv,
    pc_fd = if (is.hd(pc)) pc$mf$coefs else if (is.mvmfd(pc)) pc$coefs,
    pc_nfd = if (is.hd(pc)) pc$nf$data else if (is.mvnfd(pc)) pc$data,
    smooth_tuning = fit$smooth_tuning,
    sparse_tuning_u = fit$sparse_tuning_u,
    sparse_tuning_nfd = fit$sparse_tuning_nfd,
    sparse_tuning_fd = fit$sparse_tuning_fd,
    GCVs = fit$GCVs,
    CVs_u = fit$CVs_u, CVs_nfd = fit$CVs_nfd, CVs_fd = fit$CVs_fd,
    eBICs_u = fit$eBICs_u, eBICs_nfd = fit$eBICs_nfd, eBICs_fd = fit$eBICs_fd,
    selection_method = fit$selection_method,
    warnings = unique(out$warnings)
  )
}
