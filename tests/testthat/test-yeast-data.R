# Data preparation chunk of dev/real_data.Rmd.
testthat::local_edition(3)

test_that("yeast inputs have the expected shape", {
  skip_if_no_yeast()
  d <- yeast_data()

  expect_equal(dim(d$y), c(542L, 18L))
  expect_equal(dim(d$x), c(542L, 106L))
  expect_false(any(grepl("_YPD$", colnames(d$x))))
})

test_that("functional part (mfd) is built as in the Rmd", {
  skip_if_no_yeast()
  d <- yeast_data()

  expect_true(is.mfd(d$mfd_v1))
  expect_equal(d$mfd_v1$nobs, 542L)
  expect_equal(d$mfd_v1$basis$nbasis, 16)
  expect_equal(as.vector(d$mfd_v1$basis$supp), c(0, 119))
  expect_snapshot_large(d$mfd_v1$coefs)
})

test_that("non-functional part (nfd) is the scaled ChIP-chip matrix", {
  skip_if_no_yeast()
  d <- yeast_data()

  expect_true(is.nfd(d$nfd_obj1))
  expect_equal(d$nfd_obj1$nobs, 542L)
  expect_equal(d$nfd_obj1$nfeatures, 106L)
  expect_equal(d$nfd_obj1$features, colnames(d$x))
  expect_equal(
    unclass(d$nfd_obj1$data),
    d$x / sqrt(mean(diag(stats::var(d$x)))),
    ignore_attr = TRUE
  )
})

test_that("hd objects combine the parts", {
  skip_if_no_yeast()
  d <- yeast_data()

  expect_true(is.hd(d$hd_obj))
  expect_equal(d$hd_obj$mf$nvar, 1L)
  expect_equal(d$hd_obj$nf$nvar, 1L)
  expect_equal(d$hd_obj$mf$nobs, 542L)
  expect_equal(d$hd_obj$nf$nobs, 542L)

  # Hd(list(mfd)) gives a functional-only hd
  expect_true(is.hd(d$hd_obj1f))
  expect_null(d$hd_obj1f$nf)
  expect_equal(d$hd_obj1f$mf$coefs, d$hd_obj$mf$coefs)
})
