# Fixtures reproducing dev/real_data.Rmd (yeast cell-cycle analysis).
#
# The Rmd sets the seed once and then fits all models in sequence, so every fit
# depends on the RNG state left behind by the previous one (CV fold shuffling
# calls sample() even when no sparse tuning is requested). yeast_fits() runs
# the same calls in the same order, once per test session, and caches them.
#
# Only deviation from the Rmd: res_power_sparse uses parallel = FALSE instead of
# a 9-worker cluster. Fold assignment happens on the main process and the
# workers use no RNG, so the results are identical (checked when these tests
# were written).

skip_if_no_yeast <- function() {
  skip_on_cran()
  skip_if_not_installed("spls")
}

yeast_cache <- new.env(parent = emptyenv())

# Run expr silently (the package prints progress with cat()) and record the
# warnings it raises instead of letting them escape.
quiet_fit <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(
    {
      utils::capture.output(res <- suppressMessages(expr))
      res
    },
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = warnings)
}

yeast_data <- function() {
  if (!is.null(yeast_cache$data)) return(yeast_cache$data)

  yeast <- NULL
  utils::data("yeast", package = "spls", envir = environment())
  x <- yeast$x
  colnames(x) <- sub("_YPD$", "", colnames(x))
  y <- yeast$y

  v1 <- t(y)
  v1 <- 1 / sqrt(mean(diag(stats::var(v1)))) * v1
  t0_f <- 0
  t1_f <- 119
  t <- seq(t0_f, t1_f, length.out = ncol(y))
  bspline_basis <- fda::create.bspline.basis(c(t0_f, t1_f), 16)
  mfd_basis <- basismfd$new(bspline_basis)
  mfd_v1 <- mfd$new(t, v1, mfd_basis)
  hd_obj1f <- Hd(list(mfd_v1))

  scaling_factor_vector <- sqrt(mean(diag(stats::var(x))))
  vector_data <- 1 / scaling_factor_vector * x
  nfd_obj1 <- nfd(vector_data)

  hd_obj <- Hd(mfd_v1, nfd_obj1)

  yeast_cache$data <- list(
    x = x, y = y, t0_f = t0_f, t1_f = t1_f,
    mfd_v1 = mfd_v1, nfd_obj1 = nfd_obj1, hd_obj1f = hd_obj1f, hd_obj = hd_obj
  )
  yeast_cache$data
}

yeast_fits <- function() {
  if (!is.null(yeast_cache$fits)) return(yeast_cache$fits)
  d <- yeast_data()
  hd_obj <- d$hd_obj

  old_seed <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (!is.null(old_seed)) assign(".Random.seed", old_seed, globalenv()))

  set.seed(123)
  alpha1 <- list(a1 = 2^seq(-30, 0, length.out = 40))
  sparse_tuning_u <- as.integer(seq(0, 80, length.out = 3))
  sparse_tuning_nfd <- as.integer(seq(0, 105, length.out = 3))
  sparse_tuning_fd <- as.integer(seq(0, 15, length.out = 3))

  fits <- list()
  fits$power_smooth <- quiet_fit(Mhpca(
    hd_obj = hd_obj, method = "power", ncomp = 4, smooth_tuning = alpha1,
    centerfns = TRUE, alpha_orth = FALSE, smoothing_type = "basispen"
  ))
  fits$power_sparse <- quiet_fit(Mhpca(
    hd_obj = hd_obj, method = "power", ncomp = 4, parallel = FALSE,
    sparse_iter = 1, centerfns = TRUE, alpha_orth = FALSE,
    smooth_tuning = alpha1, smooth_GCV = TRUE, smoothing_type = "basispen",
    sparse_tuning_nfd = sparse_tuning_nfd, sparse_tuning_fd = sparse_tuning_fd,
    sparse_tuning_u = sparse_tuning_u, nfold_nfd = 10, nfold_u = 10,
    nfold_fd = 10, sparse_CV = TRUE, pen_u = TRUE, pen_nfd = TRUE, pen_fd = TRUE
  ))
  fits$joint <- quiet_fit(Mhpca(
    hd_obj = hd_obj, method = "power", ncomp = 4, smoothing_type = "basispen",
    smooth_tuning = alpha1, alpha_orth = TRUE, centerfns = TRUE,
    smooth_GCV = TRUE
  ))
  fits$eigen <- quiet_fit(Mhpca(
    hd_obj = hd_obj, method = "eigen", ncomp = 4, smooth_tuning = alpha1,
    smoothing_type = "basispen", smooth_GCV = TRUE, centerfns = TRUE
  ))
  # The Rmd assigns both of these to res_remfpca; the first one is kept here
  # too because it is run (and consumes RNG) before the second.
  fits$remfpca_sparse <- quiet_fit(Mhpca(
    hd_obj = hd_obj$mf, method = "power", ncomp = 4, centerfns = TRUE,
    alpha_orth = FALSE, smooth_tuning = alpha1, smooth_GCV = TRUE,
    smoothing_type = "basispen", sparse_tuning_fd = sparse_tuning_fd,
    sparse_tuning_u = sparse_tuning_u, nfold_u = 10, nfold_fd = 10,
    sparse_CV = TRUE, pen_u = TRUE, pen_nfd = FALSE, pen_fd = TRUE
  ))
  fits$remfpca <- quiet_fit(Mhpca(hd_obj = hd_obj$mf, method = "power", ncomp = 4))

  yeast_cache$fits <- fits
  fits
}

yeast_fit <- function(name) yeast_fits()[[name]]$value

# Score matrix exactly as computed in the Rmd: t(inprod_hd(pc_hd, hd_obj)).
yeast_scores <- function(name) t(inprod_hd(yeast_fit(name)$pc_hd, yeast_data()$hd_obj))
