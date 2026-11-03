# Regression tests for fixed bugs in the data classes (mfd, mvmfd, nfd, mvnfd,
# hd). Each test states the correct behaviour.
testthat::local_edition(3)

class_objects <- function() {
  argval <- seq(0, 2 * pi, length.out = 12)
  fourier <- Basismfd(fda::create.fourier.basis(c(0, 2 * pi), 5))
  m1 <- Mfd(argval, outer(sin(argval), seq(0.5, 1.5, length.out = 6)), fourier)
  m2 <- Mfd(argval, outer(cos(argval), seq(0.2, 1.5, length.out = 6)), fourier)
  grid <- seq(0, 1, length.out = 20)
  bspline <- Basismfd(fda::create.bspline.basis(c(0, 1), 7))
  mb <- Mfd(grid, outer(grid^2, 1:6) + outer(sin(3 * grid), 6:1), bspline)
  n1 <- nfd(matrix(seq_len(18) / 7, nrow = 6))
  n2 <- nfd(matrix(cos(seq_len(12)), nrow = 6))
  n3 <- nfd(matrix(sin(seq_len(6)), nrow = 6))
  list(m1 = m1, m2 = m2, mb = mb, n1 = n1, n2 = n2, n3 = n3)
}

test_that("scalar multiplication, negation and subtraction of hd", {
  o <- class_objects()
  h <- Hd(o$m1, o$n1)
  expect_equal((2 * h)$mf$coefs, lapply(h$mf$coefs, function(x) 2 * x))
  expect_equal((2 * h)$nf$data, lapply(h$nf$data, function(x) 2 * x), ignore_attr = TRUE)
  expect_equal((h * 2)$mf$coefs, (2 * h)$mf$coefs)
  expect_equal((-h)$mf$coefs, lapply(h$mf$coefs, function(x) -x))
  diff <- h - h
  expect_equal(max(abs(unlist(diff$mf$coefs))), 0)
  expect_equal(max(abs(unlist(diff$nf$data))), 0)
  # one-part hd objects too
  expect_equal((3 * Hd(o$n1))$nf$data[[1]], 3 * unclass(o$n1), ignore_attr = TRUE)
})

test_that("integer scalars work like doubles", {
  o <- class_objects()
  expect_equal((2L * o$m1)$coefs, (2 * o$m1)$coefs)
  expect_equal((o$m1 + 1L)$coefs, (o$m1 + 1)$coefs)
  expect_equal((2L * Mvmfd(o$m1, o$m2))$coefs, (2 * Mvmfd(o$m1, o$m2))$coefs)
  expect_equal((2L * Hd(o$m1))$mf$coefs, (2 * Hd(o$m1))$mf$coefs)
})

test_that("mvnfd addition adds matching variables", {
  o <- class_objects()
  mv <- Mvnfd(o$n1, o$n2)
  s <- mv + mv
  expect_equal(s$data[[1]], 2 * unclass(o$n1), ignore_attr = TRUE)
  expect_equal(s$data[[2]], 2 * unclass(o$n2), ignore_attr = TRUE)
  one <- Mvnfd(list(o$n1))
  expect_equal((one + one)$data[[1]], 2 * unclass(o$n1), ignore_attr = TRUE)
  expect_equal(max(abs(unlist((mv - mv)$data))), 0)
  expect_error(Mvnfd(o$n1) + Mvnfd(o$n2), "same dimension")
})

test_that("inner products of mvnfd and hd", {
  o <- class_objects()
  mv3 <- Mvnfd(o$n1, o$n2, o$n3)
  expected <- inprod_nfd(o$n1, o$n1) + inprod_nfd(o$n2, o$n2) + inprod_nfd(o$n3, o$n3)
  expect_equal(inprod_mvnfd(mv3, mv3), expected)
  # second object without a functional part: only the shared part counts
  expect_equal(inprod_hd(Hd(o$m1, o$n1), Hd(o$n1)), inprod_nfd(o$n1, o$n1))
})

test_that("subsetting hd selects parts by position", {
  o <- class_objects()
  h <- Hd(o$m1, o$n1)
  expect_true(is.null(h[, 1]$nf) && !is.null(h[, 1]$mf))
  expect_true(is.null(h[, 2]$mf) && !is.null(h[, 2]$nf))
  expect_equal(h[, 2]$nf$data[[1]], unclass(o$n1), ignore_attr = TRUE)
  expect_equal(length(h[1:3, ]), 3L)
  expect_equal(h[2:4, 2]$nf$data[[1]], unclass(o$n1)[2:4, ], ignore_attr = TRUE)
  expect_error(Hd(o$m1)[, 2], "subscript j out of bounds")
  expect_error(h[, 0], "subscript j out of bounds")
  expect_error(Mvmfd(o$m1, o$m2)[, 0], "subscript j out of bounds")
})

test_that("length() and subsetting of mfd", {
  o <- class_objects()
  h <- Hd(o$m1, o$n1)
  expect_equal(length(h), 6L)
  expect_equal(length(Hd(o$n1)), 6L)
  expect_equal(o$m1[]$coefs, o$m1$coefs)
})

test_that("scaling keeps the basis and scales the coefficients", {
  o <- class_objects()
  s <- scale(o$m1)
  expect_identical(s$basis$basis[[1]]$type, "fourier")
  factor <- s$coefs[2, 1] / o$m1$coefs[2, 1]
  expect_equal(s$coefs, factor * o$m1$coefs)
  # For a B-spline basis this equals the previous evaluate-and-refit method.
  grid <- seq(0, 1, length.out = 100)
  v <- o$mb$eval(grid)
  f <- 1 / sqrt(mean(diag(stats::var(v))))
  refit <- Mfd(grid, f * v, Basismfd(fda::create.bspline.basis(c(0, 1), 7)))
  expect_equal(scale(o$mb)$coefs, refit$coefs, tolerance = 1e-8)
  # vector-only hd can be scaled
  expect_equal(scale(Hd(o$n1))$nf$data[[1]], unclass(scale_nfd(o$n1)), ignore_attr = TRUE)
})

test_that("balanced scaling equalizes the total variation of both parts", {
  o <- class_objects()
  total_variation <- function(h) {
    nf_centred <- lapply(h$nf$data, function(x) sweep(unclass(x), 2, colMeans(x)))
    c(mf = sum(norm_mvmfd(h$mf - mean(h$mf))^2),
      nf = sum(vapply(nf_centred, function(x) sum(x^2), numeric(1))))
  }
  for (h in list(Hd(o$m1, o$n1), Hd(o$m1, o$m2, o$n1, o$n2))) {
    s <- scale(h, method = "balanced")
    tv <- total_variation(s)
    expect_equal(tv[["nf"]], tv[["mf"]])
    # the functional part is left unchanged
    expect_equal(s$mf$coefs, h$mf$coefs)
  }
})

test_that("operations that basis objects do not support give a clear error", {
  o <- class_objects()
  b <- o$m1$basis
  mb <- Mvmfd(o$m1, o$m2)$basis
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  msg <- "is not defined for basis objects"
  expect_error(b + b, paste0("`\\+` ", msg))
  expect_error(b - b, msg)
  expect_error(2 * b, paste0("`\\*` ", msg))
  expect_error(-b, msg)
  expect_error(b == b, msg)
  expect_error(b[1], paste0("`\\[` ", msg))
  expect_error(mean(b), paste0("`mean` ", msg))
  expect_error(sd(b), paste0("`sd` ", msg))
  expect_error(scale(b), paste0("`scale` ", msg))
  expect_error(plot(b), paste0("`plot` ", msg))
  expect_error(mb + mb, msg)
  expect_error(mean(mb), msg)
  # supported things still work
  expect_equal(b$nbasis, 5)
  expect_true(is.basismfd(mb[1]))
  expect_no_error(utils::capture.output(print(b)))
  expect_no_error(b$clone())
})

test_that("an mfd can be subset with an empty or unit second index", {
  o <- class_objects()
  expect_equal(o$m1[1:3, ]$coefs, o$m1[1:3]$coefs)
  expect_equal(o$m1[1:3, 1]$coefs, o$m1[1:3]$coefs)
  expect_equal(o$m1[, 1]$coefs, o$m1$coefs)
  expect_error(o$m1[, 2], "holds a single variable")
  expect_error(o$m1[1:2, 1:2], "holds a single variable")
})

test_that("bimfdplot accepts hd objects", {
  o <- class_objects()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(bimfdplot(Hd(o$m1, o$m2)))
  expect_no_error(bimfdplot(Mvmfd(o$m1, o$m2)))
  expect_error(bimfdplot(Hd(o$n1)), "no functional part")
})

test_that("2-D functional data: mean, centring and all PCA paths agree", {
  set.seed(3)
  b1 <- fda::create.bspline.basis(c(0, 1), 5)
  b2 <- fda::create.bspline.basis(c(0, 1), 4)
  m2d <- Mfd(list(seq(0, 1, length.out = 10), seq(0, 1, length.out = 8)),
             array(stats::rnorm(10 * 8 * 12), c(10, 8, 12)), Basismfd(b1, b2))
  expect_equal(as.vector(mean(m2d)$coefs), as.vector(apply(m2d$coefs, 1:2, mean)))
  centred <- center_hd(Hd(m2d))
  expect_equal(dim(centred$mf$coefs[[1]]), c(5L, 4L, 12L))
  expect_lt(max(abs(apply(centred$mf$coefs[[1]], 1:2, mean))), 1e-12)

  fit <- function(...) {
    set.seed(1)
    utils::capture.output(f <- suppressMessages(suppressWarnings(Mhpca(Hd(m2d), ncomp = 2, ...))))
    f
  }
  sequential <- fit(smooth_tuning = list(0), smooth_GCV = FALSE)
  joint <- fit(alpha_orth = TRUE, smooth_tuning = list(0), smooth_GCV = FALSE)
  eigen <- fit(method = "eigen", smooth_tuning = list(0), smooth_GCV = FALSE)
  expect_equal(sequential$values, eigen$values, tolerance = 1e-8)
  expect_equal(joint$values, eigen$values, tolerance = 1e-8)
  cosines <- abs(diag(crossprod(sequential$lsv, eigen$lsv))) /
    sqrt(colSums(sequential$lsv^2) * colSums(eigen$lsv^2))
  expect_equal(cosines, c(1, 1), tolerance = 1e-8)
})

test_that("constructor and setter errors are informative", {
  o <- class_objects()
  expect_error(Hd(Mvnfd(o$n1), Mvnfd(o$n2)), "Only one `mvnfd' object")
  mv <- Mvmfd(o$m1, o$m2)
  expect_error(mv$spars_par <- list(1:5, 1L), "must be < 5")
})
