# Characterization tests for how Mhpca() / Remfpca() reshape tuning inputs
# (vector / matrix / list, recycling, warnings, var1.. names) before calling
# the algorithms. The algorithm entry points are mocked to capture the
# reshaped arguments and stop, so nothing is fitted.
testthat::local_edition(3)

tuning_data <- function() {
  bs <- fda::create.fourier.basis(c(0, 2 * pi), 5)
  argval <- seq(0, 2 * pi, length.out = 12)
  m1 <- Mfd(argval, outer(sin(argval), seq(0.5, 1.5, length.out = 6)), Basismfd(bs))
  m2 <- Mfd(argval, outer(cos(argval), seq(0.2, 1.5, length.out = 6)), Basismfd(bs))
  n1 <- nfd(matrix(seq_len(12) / 5, nrow = 6))
  n2 <- nfd(matrix(sin(seq_len(18)), nrow = 6))
  n3 <- nfd(matrix(cos(seq_len(6)), nrow = 6))
  # 2 functional and 3 non-functional variables, so that mf$nvar and nf$nvar
  # differ and mix-ups between them are visible.
  list(hd = Hd(m1, m2, n1, n2, n3), mvmfd = Mvmfd(m1, m2))
}

# Run `expr` with the algorithm functions mocked. Returns which function was
# reached and the requested arguments, the warnings raised, or the error.
capture_dispatch <- function(expr, fns, keep) {
  captured <- NULL
  mocks <- lapply(fns, function(fn) {
    force(fn)
    function(...) {
      args <- list(...)
      captured <<- list(fn = fn, args = args[intersect(keep, names(args))])
      stop("<captured>")
    }
  })
  names(mocks) <- fns
  do.call(local_mocked_bindings, c(mocks, list(.env = environment())))

  warnings <- character()
  error <- NULL
  withCallingHandlers(
    tryCatch(utils::capture.output(suppressMessages(expr)), error = function(e) {
      if (!identical(conditionMessage(e), "<captured>")) error <<- conditionMessage(e)
    }),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(called = captured$fn, args = captured$args, warnings = warnings, error = error)
}

mhpca_fns <- c("sequential_power_hybrid", "joint_power_hybrid", "eigen_approach_hybrid")
mhpca_keep <- c("smooth_tuning", "sparse_tuning_u", "sparse_tuning_nfd", "sparse_tuning_fd", "alpha")
remfpca_fns <- c("sequential_power", "joint_power", "eigen_approach")
remfpca_keep <- c("smooth_tuning", "sparse_tuning", "alpha")

smooth_inputs <- list(
  null = NULL,
  scalar = 0.5,
  vec_ncomp = c(0.1, 0.2),
  vec_3 = c(0.1, 0.2, 0.3),
  mat_nvar_x_ncomp = matrix(c(0.1, 0.2, 0.3, 0.4), 2, 2),
  mat_1x3 = matrix(c(0.1, 0.2, 0.3), 1, 3),
  mat_nvar_x_1 = matrix(c(0.1, 0.2), 2, 1),
  list_ok = list(c(0.1, 0.2), c(0.3, 0.4)),
  list_short = list(c(0.1, 0.2)),
  list_elem_3 = list(c(0.1, 0.2, 0.3), 0.4)
)

sparse_inputs <- list(
  null = NULL,
  scalar = 1L,
  vec_ncomp = c(1L, 2L),
  vec_3 = c(1L, 2L, 3L),
  mat = matrix(1:4, 2, 2),
  list_1 = list(c(1L, 2L)),
  list_2 = list(c(1L, 2L), c(0L, 1L)),
  list_3 = list(1L, c(1L, 2L), c(0L, 1L, 2L))
)

run_cases <- function(inputs, call) {
  lapply(inputs, function(x) call(x))
}

test_that("Mhpca sequential path: smooth_tuning", {
  hd <- tuning_data()$hd
  for (grid in c(TRUE, FALSE)) {
    res <- run_cases(smooth_inputs, function(x) capture_dispatch(
      Mhpca(hd, ncomp = 2, smooth_tuning = x, smooth_GCV = grid),
      mhpca_fns, mhpca_keep
    ))
    expect_true(all(vapply(res, function(r) identical(r$called, "sequential_power_hybrid") || !is.null(r$error), logical(1))))
    expect_snapshot_small(res)
  }
})

test_that("Mhpca sequential path: sparse_tuning_fd / _nfd / _u", {
  hd <- tuning_data()$hd
  for (arg in c("sparse_tuning_fd", "sparse_tuning_nfd", "sparse_tuning_u")) {
    for (grid in c(TRUE, FALSE)) {
      res <- run_cases(sparse_inputs, function(x) {
        call_args <- list(hd, ncomp = 2, sparse_CV = grid)
        call_args[arg] <- list(x)
        capture_dispatch(do.call(Mhpca, call_args), mhpca_fns, mhpca_keep)
      })
      expect_snapshot_small(list(arg = arg, sparse_CV = grid, cases = res))
    }
  }
})

test_that("Mhpca joint and eigen paths: smooth_tuning", {
  hd <- tuning_data()$hd
  for (path in c("joint", "eigen")) {
    for (grid in c(TRUE, FALSE)) {
      res <- run_cases(smooth_inputs, function(x) capture_dispatch(
        if (path == "joint") {
          Mhpca(hd, ncomp = 2, alpha_orth = TRUE, smooth_tuning = x, smooth_GCV = grid)
        } else {
          Mhpca(hd, ncomp = 2, method = "eigen", smooth_tuning = x, smooth_GCV = grid)
        },
        mhpca_fns, mhpca_keep
      ))
      expect_snapshot_small(list(path = path, smooth_GCV = grid, cases = res))
    }
  }
})

test_that("Mhpca with only one data part", {
  d <- tuning_data()
  mf_only <- capture_dispatch(
    Mhpca(d$hd$mf, ncomp = 2, smooth_tuning = c(0.1, 0.2), sparse_tuning_fd = 1L),
    mhpca_fns, mhpca_keep
  )
  nf_only <- capture_dispatch(
    Mhpca(d$hd$nf, ncomp = 2, smooth_tuning = c(0.1, 0.2), sparse_tuning_nfd = 1L,
          sparse_tuning_fd = 1L),
    mhpca_fns, mhpca_keep
  )
  expect_snapshot_small(list(mf_only = mf_only, nf_only = nf_only))
})

test_that("Remfpca: smooth_tuning on all paths", {
  mv <- tuning_data()$mvmfd
  for (path in c("sequential", "joint", "eigen")) {
    for (grid in c(TRUE, FALSE)) {
      res <- run_cases(smooth_inputs, function(x) capture_dispatch(
        switch(path,
          sequential = Remfpca(mv, ncomp = 2, smooth_tuning = x, smooth_GCV = grid),
          joint = Remfpca(mv, ncomp = 2, alpha_orth = TRUE, smooth_tuning = x, smooth_GCV = grid),
          eigen = Remfpca(mv, method = "eigen", ncomp = 2, smooth_tuning = x, smooth_GCV = grid)
        ),
        remfpca_fns, remfpca_keep
      ))
      expect_snapshot_small(list(path = path, smooth_GCV = grid, cases = res))
    }
  }
})

test_that("Remfpca: sparse_tuning", {
  mv <- tuning_data()$mvmfd
  for (grid in c(TRUE, FALSE)) {
    res <- run_cases(sparse_inputs, function(x) capture_dispatch(
      Remfpca(mv, ncomp = 2, sparse_tuning = x, sparse_CV = grid),
      remfpca_fns, remfpca_keep
    ))
    expect_snapshot_small(list(sparse_CV = grid, cases = res))
  }
})
