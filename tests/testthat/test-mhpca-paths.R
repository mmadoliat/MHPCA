# Characterization tests for Mhpca() option combinations not reached by the
# yeast tests (see helper-paths.R). They pin current results, including errors.
testthat::local_edition(3)

for (case in names(paths_cases())) {
  local({
    case <- case
    test_that(paste("Mhpca path:", case), {
      expect_snapshot_small(fit_summary(run_path_case(paths_cases()[[case]])))
    })
  })
}

test_that("selection outputs depend on the path", {
  cases <- paths_cases()
  cv <- run_path_case(cases$cv_1se)$value
  expect_equal(cv$selection_method, "cv")
  expect_false(is.null(cv$CVs_u))
  expect_null(cv$eBICs_u)

  ebic <- run_path_case(cases$ebic_direct)$value
  expect_equal(ebic$selection_method, "ebic")
  expect_null(ebic$CVs_u)
  expect_false(is.null(ebic$eBICs_u))

  # CVs_u is filled whenever u is penalized (also without pen_fd)
  expect_false(is.null(run_path_case(cases$nf_only_cv)$value$CVs_u))
  expect_false(is.null(run_path_case(cases$cv_pen_u_only)$value$CVs_u))

  # Blocks that are not tuned report no selection
  only_u <- run_path_case(cases$cv_pen_u_only)$value
  expect_length(only_u$sparse_tuning_fd, 0)
  expect_length(only_u$sparse_tuning_nfd, 0)

  # Single-part data works on every path, with or without smoothing grids
  for (case in c("joint_mf_only", "joint_nf_only", "fixed_smooth_nf_only", "cv_pen_fd_only")) {
    expect_null(fit_summary(run_path_case(cases[[case]]))$error, label = case)
  }
})

test_that("eBIC results do not depend on the order of the smoothing grid", {
  for (case in c("ebic_direct", "ebic_pen_fd_only")) {
    args <- paths_cases()[[case]]
    reversed <- args
    reversed$smooth_tuning <- list(rev(args$smooth_tuning[[1]]))
    a <- fit_summary(run_path_case(args))
    b <- fit_summary(run_path_case(reversed))
    expect_equal(b$values, a$values, info = case)
    expect_equal(b$sparse_tuning_fd, a$sparse_tuning_fd, info = case)
    expect_equal(b$sparse_tuning_nfd, a$sparse_tuning_nfd, info = case)
    expect_equal(b$sparse_tuning_u, a$sparse_tuning_u, info = case)
  }
})

test_that("mhpca fields are read-only with the current error messages", {
  fit <- run_path_case(paths_cases()$cv_1se)$value
  fields <- c("pc_hd", "lsv", "values", "smooth_tuning", "sparse_tuning_u",
              "sparse_tuning_nfd", "sparse_tuning_fd", "GCVs", "CVs_u", "CVs_nfd",
              "CVs_fd", "eBICs_u", "eBICs_nfd", "eBICs_fd", "selection_method", "mean_hd")
  # Every field can be read...
  for (f in fields) expect_no_error(fit[[f]])
  # ...and assigning to any of them fails with a fixed message.
  messages <- vapply(fields, function(f) {
    tryCatch({ fit[[f]] <- 1; "<no error>" }, error = conditionMessage)
  }, character(1))
  expect_snapshot_small(messages)
  # Public (non-private) members of the object.
  expect_snapshot_small(sort(setdiff(ls(fit), ".__enclos_env__")))
})

test_that("Mhpca() and mhpca$new() pass arguments the same way", {
  # fixed_tuning sets sparse_CV and smooth_GCV, which mhpca$new() has no
  # defaults for; smoothing_type is the one other default that differs.
  args <- paths_cases()$fixed_tuning
  via_wrapper <- run_path_case(args)
  set.seed(1)
  via_new <- quiet_fit(do.call(mhpca$new, c(args, list(smoothing_type = "basispen"))))
  expect_equal(fit_summary(via_new), fit_summary(via_wrapper))

  # Both signatures are public API: pin names, order and defaults.
  expect_equal(formals(Mhpca)$smoothing_type, "basispen")
  expect_equal(formals(mhpca$public_methods$initialize)$smoothing_type, "coefpen")
  expect_snapshot_small(lapply(formals(Mhpca), deparse))
  expect_snapshot_small(lapply(formals(mhpca$public_methods$initialize), deparse))
})

test_that("parallel = TRUE gives the same results as parallel = FALSE", {
  # Opt-in: workers load the *installed* MHPCA, so this installs the current
  # source into a temporary library and points the workers at it.
  skip_if_not(identical(Sys.getenv("MHPCA_TEST_PARALLEL"), "true"),
              "set MHPCA_TEST_PARALLEL=true to run")
  lib <- withr::local_tempdir()
  utils::install.packages(testthat::test_path("..", ".."), lib = lib, repos = NULL,
                          type = "source", quiet = TRUE)
  withr::local_envvar(R_LIBS = paste(c(lib, .libPaths()), collapse = .Platform$path.sep))
  foreach::registerDoSEQ()

  for (case in c("cv_1se", "cv_iter2", "mf_only_cv", "nf_only_cv")) {
    args <- paths_cases()[[case]]
    sequential <- run_path_case(c(args, list(parallel = FALSE)))
    parallel <- run_path_case(c(args, list(parallel = TRUE)))
    expect_equal(fit_summary(parallel), fit_summary(sequential), info = case)
  }
})

test_that("invalid choice arguments fail early with a clear message", {
  d <- paths_data()
  expect_error(Mhpca(d$hd, method = "Power"), "`method` must be one of \"power\", \"eigen\"")
  expect_error(Mhpca(d$hd, smoothing_type = "basis"), "`smoothing_type` must be one of")
  expect_error(Mhpca(d$hd, sparse_type_u = "lasso"), "`sparse_type_u` must be one of")
  expect_error(Mhpca(d$hd, sparse_type_fd = "soft-threshold"), "`sparse_type_fd` must be one of")
  expect_error(Mhpca(d$hd, cv.pick = "max"), "`cv.pick` must be one of")
  expect_error(Mhpca(d$hd, alpha_orth = "yes"), "`alpha_orth` must be TRUE or FALSE")
  expect_error(Remfpca(d$hd$mf, ncomp = 2, method = "eig"), "`method` must be one of")
  expect_error(Remfpca(d$hd$mf, ncomp = 2, sparse_type = "soft-threshold"), "`sparse_type` must be one of")
})

test_that("alpha_orth accepts logical and \"TRUE\"/\"FALSE\" strings alike", {
  d <- paths_data()
  joint <- fit_summary(run_path_case(list(hd_obj = d$hd, ncomp = 2, alpha_orth = TRUE)))
  joint_chr <- fit_summary(run_path_case(list(hd_obj = d$hd, ncomp = 2, alpha_orth = "TRUE")))
  expect_equal(joint_chr, joint)
})

test_that("progress text is a message that can be suppressed", {
  d <- paths_data()
  messages <- testthat::capture_messages(utils::capture.output(Mhpca(d$hd, ncomp = 1)))
  expect_true(any(grepl("Computing the 1st PC", messages)))
  printed <- utils::capture.output(suppressMessages(Mhpca(d$hd, ncomp = 1)))
  expect_false(any(grepl("Computing", printed)))
})
