# "Model" chunk of dev/real_data.Rmd: the six Mhpca() fits.
# These pin down the CURRENT outputs of the package (regression tests), not
# necessarily correct ones; see the notes on individual tests.
testthat::local_edition(3)

# Fields of an mhpca object that hold numeric results.
fit_numbers <- function(fit) {
  pc <- fit$pc_hd
  list(
    values = fit$values,
    pc_fd_coefs = if (is.hd(pc)) pc$mf$coefs else if (is.mvmfd(pc)) pc$coefs,
    pc_nfd = if (is.hd(pc)) pc$nf$data
  )
}

test_that("all fits return mhpca objects with 4 components", {
  skip_if_no_yeast()
  for (name in names(yeast_fits())) {
    fit <- yeast_fit(name)
    expect_s3_class(fit, "mhpca")
    expect_length(fit$values, 4)
    expect_equal(dim(fit$lsv), c(542L, 4L), info = name)
  }
})

test_that("pc_hd type depends on the input", {
  skip_if_no_yeast()
  expect_true(is.hd(yeast_fit("power_smooth")$pc_hd))
  expect_true(is.hd(yeast_fit("power_sparse")$pc_hd))
  expect_true(is.hd(yeast_fit("joint")$pc_hd))
  expect_true(is.hd(yeast_fit("eigen")$pc_hd))
  # functional-only input returns an mvmfd, not an hd
  expect_true(is.mvmfd(yeast_fit("remfpca")$pc_hd))
  expect_true(is.mvmfd(yeast_fit("remfpca_sparse")$pc_hd))
})

test_that("sequential power, smoothing only", {
  skip_if_no_yeast()
  fit <- yeast_fit("power_smooth")
  expect_snapshot_small(fit_numbers(fit))
  expect_snapshot_large(fit$lsv)
  expect_snapshot_small(fit$smooth_tuning)
  expect_snapshot_small(fit$GCVs)
  expect_snapshot_small(fit$sparse_tuning_u)
})

test_that("sequential power, smoothing + sparsity on u, nfd and fd", {
  skip_if_no_yeast()
  fit <- yeast_fit("power_sparse")

  # Selected sparsity (number of entries set to zero) per PC. With the current
  # code CV selects no sparsity for nfd or fd, so all 106 TFs stay non-zero;
  # the comments in dev/real_data.Rmd (6 / 8 / 12 TFs) came from older code.
  expect_equal(unname(unlist(fit$sparse_tuning_u)), c(40L, 80L, 80L, 80L))
  expect_equal(unname(unlist(fit$sparse_tuning_nfd)), c(0, 0, 0, 0))
  expect_equal(unname(unlist(fit$sparse_tuning_fd)), c(0, 0, 0, 0))
  expect_equal(colSums(fit$lsv == 0), c(40, 80, 80, 80), ignore_attr = TRUE)

  expect_equal(fit$selection_method, "cv")
  expect_snapshot_small(fit_numbers(fit))
  expect_snapshot_large(fit$lsv)
  expect_snapshot_small(fit$smooth_tuning)
  expect_snapshot_small(fit$GCVs)
  expect_snapshot_small(fit$CVs_u)
  expect_snapshot_small(fit$CVs_nfd)
  expect_snapshot_small(fit$CVs_fd)
  expect_null(fit$eBICs_u)
})

test_that("joint power, smoothing only", {
  skip_if_no_yeast()
  fit <- yeast_fit("joint")
  expect_snapshot_small(fit_numbers(fit))
  expect_snapshot_large(fit$lsv)
  expect_snapshot_small(fit$smooth_tuning)
  expect_snapshot_small(fit$GCVs)
})

test_that("eigen approach, smoothing only", {
  skip_if_no_yeast()
  fit <- yeast_fit("eigen")
  expect_snapshot_small(fit_numbers(fit))
  expect_snapshot_large(fit$lsv)
  # The eigen path reports the selected smoothing parameter and GCV scores.
  expect_length(fit$GCVs, 40)
  expect_snapshot_small(fit$smooth_tuning)
  expect_snapshot_small(fit$GCVs)
})

test_that("functional-only fits (res_remfpca in the Rmd)", {
  skip_if_no_yeast()
  sparse <- yeast_fit("remfpca_sparse")
  expect_snapshot_small(fit_numbers(sparse))
  expect_snapshot_large(sparse$lsv)
  expect_snapshot_small(sparse$sparse_tuning_u)
  expect_snapshot_small(sparse$sparse_tuning_fd)
  expect_snapshot_small(sparse$CVs_u)
  expect_snapshot_small(sparse$CVs_fd)

  plain <- yeast_fit("remfpca")
  expect_snapshot_small(fit_numbers(plain))
  expect_snapshot_large(plain$lsv)
})

test_that("warnings raised while fitting", {
  skip_if_no_yeast()
  warnings <- lapply(yeast_fits(), function(f) unique(f$warnings))
  expect_snapshot_small(warnings)
})
