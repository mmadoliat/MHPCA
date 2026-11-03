# Plotting of fitted objects and CV scores.
testthat::local_edition(3)

quiet_mhpca <- function(...) {
  set.seed(1)
  utils::capture.output(fit <- suppressMessages(suppressWarnings(Mhpca(...))))
  fit
}

test_that("plot() works for every fit with a functional part", {
  d <- paths_data()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(plot(quiet_mhpca(d$hd, ncomp = 2), ask = FALSE))
  expect_no_error(plot(quiet_mhpca(d$mf_only, ncomp = 2), ask = FALSE))
  expect_no_error(plot(quiet_mhpca(d$hd, ncomp = 2), comp_index = 2, ask = FALSE))
  utils::capture.output(r1 <- suppressMessages(Remfpca(d$hd$mf, ncomp = 2)))
  expect_no_error(plot(r1, ask = FALSE))
  utils::capture.output(r2 <- suppressMessages(suppressWarnings(Remfpca(Mvmfd(d$hd$mf[, 1], d$hd$mf[, 1]), ncomp = 2))))
  expect_no_error(plot(r2, ask = FALSE))
})

test_that("plot() of a vector-only fit explains why it cannot plot", {
  d <- paths_data()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_error(plot(quiet_mhpca(d$nf_only, ncomp = 2), ask = FALSE),
               "no functional part to plot")
})

test_that("plot_cv() plots CV scores and explains when there are none", {
  d <- paths_data()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  cv <- quiet_mhpca(d$hd, ncomp = 2, sparse_tuning_u = c(0L, 5L), nfold_u = 2, pen_u = TRUE)
  expect_no_error(plot_cv(cv, type = "u", PC = 1))
  ebic <- quiet_mhpca(d$hd, ncomp = 2, sparse_tuning_u = c(0L, 5L), pen_u = TRUE,
                      sparse_selection_method = "ebic")
  expect_error(plot_cv(ebic, type = "u", PC = 1), "selected sparsity by eBIC")
  expect_error(plot_cv(cv, type = "nfd"), "no cross-validation scores")
})

test_that("ordinal numbers", {
  expect_equal(vapply(c(1, 2, 3, 4, 10, 11, 12, 13, 21, 22, 23, 101, 111, 112), ordinal, ""),
               c("1st", "2nd", "3rd", "4th", "10th", "11th", "12th", "13th", "21st", "22nd",
                 "23rd", "101st", "111th", "112th"))
})
