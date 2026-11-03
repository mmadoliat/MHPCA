# Characterization tests for the exported dispatchers length(), mean(), sd(),
# plot() and scale() (R/mfdMethods.R). They pin which helper each class is
# routed to, that `...` is passed through, the fall-back to base/stats for
# other objects, and the values produced on small objects.
testthat::local_edition(3)

# A stand-in helper that records which route was taken and what it received.
route <- function(name) function(...) list(route = name, args = list(...))

# Bare objects carrying only a class attribute are enough to test routing.
fake <- function(...) structure(list(), class = c(...))

small_objects <- function() {
  bs <- fda::create.fourier.basis(c(0, 2 * pi), 5)
  argval <- seq(0, 2 * pi, length.out = 12)
  X1 <- outer(sin(argval), seq(0.5, 1.5, length.out = 6))
  X2 <- outer(cos(argval), seq(0.2, 1.5, length.out = 6))
  m1 <- Mfd(argval, X1, Basismfd(bs))
  m2 <- Mfd(argval, X2, Basismfd(bs))
  v1 <- nfd(matrix(seq_len(18) / 7, nrow = 6))
  v2 <- nfd(matrix(cos(seq_len(12)), nrow = 6))
  list(
    mfd = m1, mvmfd = Mvmfd(m1, m2), nfd = v1, mvnfd = Mvnfd(v1, v2),
    hd = Hd(m1, v1), img = img(matrix(1:4, 2), matrix(5:8, 2))
  )
}

test_that("mean() routes each class to its helper and passes ... through", {
  local_mocked_bindings(
    mean_mfd = route("mean_mfd"), mean_mvmfd = route("mean_mvmfd"),
    mean_hd = route("mean_hd"), mean_nfd = route("mean_nfd"),
    mean_mvnfd = route("mean_mvnfd")
  )
  expect_equal(mean(fake("mfd", "R6"))$route, "mean_mfd")
  expect_equal(mean(fake("mvmfd", "R6"))$route, "mean_mvmfd")
  expect_equal(mean(fake("hd", "R6"))$route, "mean_hd")
  expect_equal(mean(fake("nfd", "matrix"))$route, "mean_nfd")
  expect_equal(mean(fake("img", "nfd", "matrix"))$route, "mean_nfd")
  expect_equal(mean(fake("mvnfd", "R6"))$route, "mean_mvnfd")
  expect_equal(mean(fake("mfd", "R6"), extra = 1)$args$extra, 1)
})

test_that("sd() routes each class to its helper and passes ... through", {
  local_mocked_bindings(
    sd_mfd = route("sd_mfd"), sd_mvmfd = route("sd_mvmfd"),
    sd_hd = route("sd_hd"), sd_nfd = route("sd_nfd"),
    sd_mvnfd = route("sd_mvnfd")
  )
  expect_equal(sd(fake("mfd", "R6"))$route, "sd_mfd")
  expect_equal(sd(fake("mvmfd", "R6"))$route, "sd_mvmfd")
  expect_equal(sd(fake("hd", "R6"))$route, "sd_hd")
  expect_equal(sd(fake("nfd", "matrix"))$route, "sd_nfd")
  expect_equal(sd(fake("img", "nfd", "matrix"))$route, "sd_nfd")
  expect_equal(sd(fake("mvnfd", "R6"))$route, "sd_mvnfd")
  expect_equal(sd(fake("nfd", "matrix"), na.rm = TRUE)$args$na.rm, TRUE)
})

test_that("plot() routes each class to its helper and passes ... through", {
  local_mocked_bindings(
    plot_mfd = route("plot_mfd"), plot_mvmfd = route("plot_mvmfd"),
    plot_hd = route("plot_hd"), plot_mhpca = route("plot_mhpca"),
    plot_remfpca = route("plot_remfpca")
  )
  expect_equal(plot(fake("mfd", "R6"))$route, "plot_mfd")
  expect_equal(plot(fake("mvmfd", "R6"))$route, "plot_mvmfd")
  expect_equal(plot(fake("hd", "R6"))$route, "plot_hd")
  expect_equal(plot(fake("mhpca", "R6"))$route, "plot_mhpca")
  expect_equal(plot(fake("remfpca", "R6"))$route, "plot_remfpca")
  expect_equal(plot(fake("mhpca", "R6"), ask = FALSE)$args$ask, FALSE)
})

test_that("scale() routes each class to its helper and passes ... through", {
  local_mocked_bindings(
    scale_mfd = route("scale_mfd"), scale_mvmfd = route("scale_mvmfd"),
    scale_hd = route("scale_hd"), scale_nfd = route("scale_nfd"),
    scale_mvnfd = route("scale_mvnfd")
  )
  expect_equal(scale(fake("mfd", "R6"))$route, "scale_mfd")
  expect_equal(scale(fake("mvmfd", "R6"))$route, "scale_mvmfd")
  expect_equal(scale(fake("hd", "R6"))$route, "scale_hd")
  expect_equal(scale(fake("nfd", "matrix"))$route, "scale_nfd")
  expect_equal(scale(fake("img", "nfd", "matrix"))$route, "scale_nfd")
  expect_equal(scale(fake("mvnfd", "R6"))$route, "scale_mvnfd")
  expect_equal(scale(fake("hd", "R6"), method = "balanced")$args$method, "balanced")
})

test_that("other objects fall back to base / stats", {
  expect_equal(length(1:5), 5L)
  expect_equal(length(list(1, 2)), 2L)
  expect_equal(mean(c(1, 2, 6)), 3)
  expect_equal(mean(c(1, NA, 3), na.rm = TRUE), 2)
  expect_equal(sd(c(1, 2, 3)), 1)
  m <- matrix(c(1, 2, 3, 4, 6, 8), 3)
  expect_equal(scale(m), base::scale(m))
  expect_equal(scale(m, center = FALSE), base::scale(m, center = FALSE))
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(plot(1:3, 3:1))
  expect_no_error(plot(nfd(matrix(1:6, 3))))
})

test_that("length() on the package classes", {
  obj <- small_objects()
  expect_equal(length(obj$mfd), 6L)
  expect_equal(length(obj$mvmfd), 6L)
  expect_equal(length(obj$nfd), 6L)
  expect_equal(length(obj$img), 2L)
  expect_equal(length(obj$hd), 6L)
  expect_equal(length(Hd(obj$nfd)), 6L)
  # mvnfd and basismfd are not handled specially: base::length() of the R6
  # environment (its number of bindings) is returned.
  expect_equal(length(obj$mvnfd), base::length(obj$mvnfd))
  expect_equal(length(obj$mfd$basis), base::length(obj$mfd$basis))
})

test_that("mean(), sd() and scale() values on small objects", {
  obj <- small_objects()

  expect_snapshot_small(list(
    mean_mfd = mean(obj$mfd)$coefs,
    mean_mvmfd = mean(obj$mvmfd)$coefs,
    mean_nfd = unclass(mean(obj$nfd)),
    mean_mvnfd = mean(obj$mvnfd)$data,
    mean_hd_mf = mean(obj$hd)$mf$coefs,
    mean_hd_nf = mean(obj$hd)$nf$data
  ))
  expect_snapshot_small(list(
    sd_mfd = sd(obj$mfd)$coefs,
    sd_mvmfd = sd(obj$mvmfd)$coefs,
    sd_nfd = unclass(sd(obj$nfd)),
    sd_mvnfd = sd(obj$mvnfd)$data,
    sd_hd_mf = sd(obj$hd)$mf$coefs,
    sd_hd_nf = sd(obj$hd)$nf$data
  ))
  expect_snapshot_small(list(
    scale_mfd = scale(obj$mfd)$coefs,
    scale_mvmfd = scale(obj$mvmfd)$coefs,
    scale_nfd = unclass(scale(obj$nfd)),
    scale_mvnfd = scale(obj$mvnfd)$data,
    scale_hd_mf = scale(obj$hd)$mf$coefs,
    scale_hd_nf = scale(obj$hd)$nf$data
  ))
})
