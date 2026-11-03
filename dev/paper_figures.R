# =============================================================================
#  dev/paper_figures.R
#
#  Figures for "Hybrid Principal Component Analysis of Functional and Vector
#  Data via a Penalized Hybrid Singular Value Decomposition".
#
#  1. Runs a reduced version of the simulation study of Section 6 with MHPCA:
#     HPCA (no penalties) and ReHPCA (roughness + sparsity penalties, tuned by
#     GCV and CV), plus "separate" baselines that analyse the curves alone
#     (FPCA / ReFPCA) or the vectors alone (PCA / sparse PCA).
#  2. Builds the figures as objects p1-p9: a schematic of the method, the
#     simulation design, estimation results, tuning and explained variance.
#  3. Displays the figures. Nothing is written to disk: saving is a separate
#     step at the end of the file and only runs when `save_figures <- TRUE`.
#
#  Usage, from the package root (the folder that contains DESCRIPTION):
#     source("dev/paper_figures.R")
#  With the default settings the simulation takes about 10 minutes on 7
#  worker processes (8-core laptop). Lower `n_rep` for a quick look, or set
#  `sim_cache` to reuse the simulation results while you work on the figures;
#  after one run you can also re-run Section 6 alone.
# =============================================================================


## 0. Settings ------------------------------------------------------------------

seed         <- 20261012   # master seed; every data set derives its own seed
n_rep        <- 20         # replicates per setting
N            <- 500        # units per data set (paper: 500, 1000, 2000)
W            <- 30         # vector dimension m (paper: 30, 60, 90); multiple of 3
sigma_levels <- c(0.1, 0.3, 0.5)              # noise SD (paper)
scenarios    <- list(I   = c(0.05, 0.60, 0.35), # group sizes N1:N2:N3 (paper)
                     II  = c(0.15, 0.50, 0.35),
                     III = c(0.25, 0.40, 0.35))
kappa        <- c(70, 65, 60)  # scale factors of the three components (paper)
n_grid       <- 101        # time points per curve on [0, 1]
n_basis      <- 20         # cubic B-spline basis functions for the curves (d)

# Tuning grids. alpha (roughness, GCV): the paper's yeast grid 2^-10, ..., 2^1
# is for a [0, 119] time axis; the penalty-to-Gram ratio scales with the
# domain length^-4, so on [0, 1] that grid is about 2^-38, ..., 2^-27. The grid
# below covers it and continues to much heavier smoothing.
alpha_grid   <- 2^seq(-38, -10, by = 2)
grid_u       <- as.integer(round(seq(0, 0.70, by = 0.05) * N)) # zero scores
grid_fd      <- 0:(n_basis - 2)                                 # zero B-spline coefficients
grid_nfd     <- seq(0L, W - 2L, by = 2L)                        # zero vector loadings
nfold        <- 10         # CV folds for all blocks (must divide n_basis, W and N)

n_cores      <- max(1L, parallel::detectCores() - 1L) # worker processes
sigma_show   <- 0.3        # noise level used by figures that show one setting
sim_cache    <- NULL       # e.g. "dev/paper_figures_sim.rds" to cache the simulation
use_load_all <- TRUE       # TRUE: load the package source with pkgload::load_all()
preview_at_size <- FALSE   # TRUE: open each figure in its own window at print size
save_figures <- FALSE      # nothing is written to disk unless TRUE (Section 9)


## 1. Setup ---------------------------------------------------------------------

needed <- c("ggplot2", "patchwork", "scales", "fda", if (use_load_all) "pkgload")
missing_pkgs <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop("Please install: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}
library(ggplot2)
library(patchwork)

find_pkg_root <- function() {
  for (dir in c(".", "..")) {
    desc <- file.path(dir, "DESCRIPTION")
    if (file.exists(desc) && identical(unname(read.dcf(desc, "Package")[1, 1]), "MHPCA")) {
      return(normalizePath(dir))
    }
  }
  stop("Run this script from the MHPCA package root (or its dev/ folder).", call. = FALSE)
}

# Also called on the worker processes of the simulation.
load_mhpca <- function(root, use_load_all) {
  if (use_load_all) {
    suppressMessages(pkgload::load_all(root, quiet = TRUE))
  } else {
    suppressPackageStartupMessages(library(MHPCA))
  }
  invisible(NULL)
}

pkg_root <- find_pkg_root()
load_mhpca(pkg_root, use_load_all)
stopifnot(W %% 3 == 0, n_basis %% nfold == 0, W %% nfold == 0, N %% nfold == 0,
          sigma_show %in% sigma_levels)


## 2. Figure style --------------------------------------------------------------

# Design size of each figure in inches (6.5 in = full text width). The fonts
# are sized for these dimensions, so preview and export at these sizes.
fig_sizes <- list(
  p1 = c(6.5, 3.3), p2 = c(6.5, 6.3), p3 = c(6.5, 4.9), p4 = c(6.5, 5.0),
  p5 = c(6.5, 4.9), p6 = c(6.5, 4.6), p7 = c(6.5, 6.4), p8 = c(6.5, 4.8),
  p9 = c(6.5, 4.6)
)
base_pt <- 8  # base font size (pt) at the design size

theme_paper <- function(base_size = base_pt) {
  theme_bw(base_size = base_size) +
    theme(
      text              = element_text(colour = "grey10"),
      panel.grid.major  = element_line(colour = "grey92", linewidth = 0.25),
      panel.grid.minor  = element_blank(),
      panel.border      = element_rect(fill = NA, colour = "grey35", linewidth = 0.3),
      axis.ticks        = element_line(colour = "grey35", linewidth = 0.3),
      axis.ticks.length = unit(2, "pt"),
      axis.text         = element_text(size = rel(0.85), colour = "grey20"),
      axis.title        = element_text(size = rel(0.95)),
      strip.background  = element_rect(fill = "grey94", colour = "grey35", linewidth = 0.3),
      strip.text        = element_text(size = rel(0.9), colour = "grey10",
                                       margin = margin(2, 2, 2, 2)),
      legend.position   = "bottom",
      legend.key.size   = unit(9, "pt"),
      legend.text       = element_text(size = rel(0.85)),
      legend.title      = element_text(size = rel(0.9)),
      legend.margin     = margin(0, 0, 0, 0),
      legend.box.spacing = unit(3, "pt"),
      plot.title        = element_text(size = rel(0.95), face = "bold", hjust = 0,
                                       margin = margin(b = 3)),
      plot.tag          = element_text(size = rel(1.2), face = "bold"),
      plot.margin       = margin(4, 5, 4, 4)
    )
}

# Colour-blind-safe palettes (Okabe-Ito). Regularized methods are cool,
# unregularized warm; hybrid methods blue/orange, separate ones green/vermillion.
pal_method <- c(ReHPCA = "#0072B2", HPCA = "#E69F00", ReFPCA = "#009E73",
                FPCA = "#D55E00", sPCA = "#009E73", PCA = "#D55E00")
lab_method <- c(ReHPCA = "ReHPCA (proposed)", HPCA = "HPCA (no penalties)",
                ReFPCA = "ReFPCA (curves only)", FPCA = "FPCA (curves only)",
                sPCA = "Sparse PCA (vectors only)", PCA = "PCA (vectors only)")
# Component m and group S_m (the units in which component m is absent).
pal_comp <- c(`1` = "#CC79A7", `2` = "#56B4E9", `3` = "#7F7F00")
col_truth <- "grey5"
pt <- function(x) x / ggplot2::.pt   # font size in points -> ggplot text size

# Sequential (one hue, light to dark) and diverging (blue / grey / red) fills.
scale_fill_seq <- function(...) {
  scale_fill_gradientn(colours = c("#FFFFFF", "#CDE2FB", "#86B6EF", "#3987E5",
                                   "#1C5CAB", "#0D366B"), ...)
}
scale_fill_div <- function(lim, ...) {
  scale_fill_gradient2(low = "#2166AC", mid = "#F0EFEC", high = "#B2182B",
                       midpoint = 0, limits = c(-lim, lim), ...)
}

# Facet labels are plotmath; plain text is quoted.
f_scen  <- function(x, long = FALSE) {
  lab <- if (long) {
    sprintf("'Scenario %s:'~~abs(S[1])*' : '*abs(S[2])*' : '*abs(S[3]) == '%s'",
            names(scenarios),
            vapply(scenarios, function(p) paste0(100 * p, "%", collapse = " : "), ""))
  } else {
    sprintf("'Scenario %s'", names(scenarios))
  }
  factor(x, levels = names(scenarios), labels = lab)
}
f_sigma <- function(x) factor(x, levels = sigma_levels, labels = paste("sigma ==", sigma_levels))
f_comp  <- function(x, what = c("psi", "v", "u", "text")) {
  lab <- switch(match.arg(what),
    psi  = sprintf("'Component %d ('*psi[%d]*')'", 1:3, 1:3),
    v    = sprintf("'Component %d ('*bold(v)[%d]*')'", 1:3, 1:3),
    u    = sprintf("'Component %d (absent in '*S[%d]*')'", 1:3, 1:3),
    text = sprintf("'Component %d'", 1:3))
  factor(x, levels = 1:3, labels = lab)
}
s_labels <- parse(text = sprintf("S[%d]", 1:3))


## 3. Data-generating process (Section 6) ----------------------------------------

t_grid <- seq(0, 1, length.out = n_grid)   # observation grid of the curves
t_eval <- seq(0, 1, length.out = 1001)     # grid for errors and plotted curves
keep_eval <- seq(1, length(t_eval), by = 5) # points of t_eval stored for figures

# Trapezoidal inner product of two functions sampled on x.
trapz_ip <- function(f, g, x) {
  h <- f * g
  sum(diff(x) * (h[-1] + h[-length(h)]) / 2)
}

# Functional loadings on I1 = [0, 1/3], I2 = (1/3, 2/3), I3 = [2/3, 1].
psi_raw <- function(t) {
  I1 <- t <= 1 / 3
  I3 <- t >= 2 / 3
  I2 <- !I1 & !I3
  cbind(ifelse(I1, sin(2 * pi * t), ifelse(I3, -sin(pi * t), 0)),
        ifelse(I2, sin(2 * pi * t), 0),
        ifelse(I3, sin(2 * pi * t), 0))
}
interval_of <- function(t) ifelse(t <= 1 / 3, 1L, ifelse(t < 2 / 3, 2L, 3L))
t_fine <- seq(0, 1, length.out = 30001)
psi_norms <- apply(psi_raw(t_fine), 2, function(f) sqrt(trapz_ip(f, f, t_fine)))
psi_true <- function(t) sweep(psi_raw(t), 2, psi_norms, "/")  # unit L2 norm
psi_obs  <- psi_true(t_grid)
psi_eval <- apply(psi_true(t_eval), 2, function(f) f / sqrt(trapz_ip(f, f, t_eval)))

# One hybrid data set: z_i = sum_m kappa_m rho_im v_m + e_i (Section 6).
simulate_hybrid <- function(N, W, sigma, props, kappa, psi_obs) {
  # score variances lambda_m^(k): rows = group S_k, columns = component m
  lam <- rbind(c(0, 0.5, 0.9), c(0.9, 0, 0.5), c(0.5, 0.9, 0))
  sizes <- round(N * props)
  sizes[3] <- N - sum(sizes[1:2])
  group <- rep(1:3, sizes)
  b <- W / 3
  V <- matrix(0, W, 3)
  for (m in 1:3) {
    V[(m - 1) * b + seq_len(b), m] <- stats::rnorm(b, c(5, 3, 1)[m], 0.1)
  }
  V <- sweep(V, 2, sqrt(colSums(V^2)), "/")
  rho <- matrix(stats::rnorm(N * 3), N, 3) * sqrt(lam[group, ])
  scores <- sweep(rho, 2, kappa, "*")
  n_t <- nrow(psi_obs)
  Y <- psi_obs %*% t(scores) + matrix(stats::rnorm(n_t * N, sd = sigma), n_t, N)
  X <- scores %*% t(V) + matrix(stats::rnorm(N * W, sd = sigma), N, W)
  list(Y = Y, X = X, V = V, group = group, scores = scores)
}

# Number of B-spline coefficients of psi_m that should be zero: basis functions
# whose support lies inside the region where psi_m vanishes.
true_zero_coefs <- function(basis) {
  knots <- c(rep(basis$rangeval[1], 4), basis$params, rep(basis$rangeval[2], 4))
  k <- seq_len(basis$nbasis)
  inside <- t(vapply(k, function(i) {
    tt <- seq(knots[i], knots[i + 4], length.out = 202)[-c(1, 202)]
    colSums(abs(psi_raw(tt)) > 0) == 0
  }, logical(3)))
  colSums(inside)
}


# Smallest Err(psi_m) that any estimate in the span of the basis can reach
# (the psi_m jump at 1/3 and 2/3): the unit-norm L2 projection of psi_m.
basis_floor <- function(basis, t = t_fine) {
  w <- (c(diff(t), 0) + c(0, diff(t))) / 2
  B <- fda::eval.basis(t, basis)
  P <- psi_true(t)
  proj <- B %*% solve(crossprod(B * w, B), crossprod(B * w, P))
  proj <- sweep(proj, 2, sqrt(colSums(proj^2 * w)), "/")
  sqrt(colSums((proj - P)^2 * w))
}


## 4. Fitting and evaluation ----------------------------------------------------

# Run an Mhpca() call without its progress output; keep warnings and errors.
run_quietly <- function(expr) {
  warns <- character()
  value <- tryCatch(
    withCallingHandlers(
      {
        utils::capture.output(out <- suppressMessages(expr))
        out
      },
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  ok <- !inherits(value, "error")
  list(value = if (ok) value, error = if (ok) NA_character_ else conditionMessage(value),
       warnings = unique(warns))
}

# The six fits. Penalized blocks are tuned as in Section 4.2: sparsity by
# 10-fold CV (1-SE rule) on the degree of sparsity, then smoothing by GCV.
method_specs <- function(mf, nf, cfg) {
  smooth <- list(smooth_tuning = list(cfg$alpha_grid))
  sp_u   <- list(sparse_tuning_u = cfg$grid_u, nfold_u = cfg$nfold, pen_u = TRUE)
  sp_fd  <- list(sparse_tuning_fd = cfg$grid_fd, nfold_fd = cfg$nfold, pen_fd = TRUE)
  sp_nfd <- list(sparse_tuning_nfd = cfg$grid_nfd, nfold_nfd = cfg$nfold, pen_nfd = TRUE)
  list(
    HPCA   = list(hd_obj = Hd(mf, nf)),
    ReHPCA = c(list(hd_obj = Hd(mf, nf)), smooth, sp_u, sp_fd, sp_nfd),
    FPCA   = list(hd_obj = Hd(mf)),
    ReFPCA = c(list(hd_obj = Hd(mf)), smooth, sp_u, sp_fd),
    PCA    = list(hd_obj = Hd(nf)),
    sPCA   = c(list(hd_obj = Hd(nf)), sp_u, sp_nfd)
  )
}

# Functional loadings (evaluated and as coefficients), vector loadings and
# scores of a fit; parts that a fit does not have are NULL.
loadings_of <- function(fit, t_eval) {
  pc <- fit$pc_hd
  fpart <- if (inherits(pc, "hd")) pc$mf else if (inherits(pc, "mvmfd")) pc else NULL
  vpart <- if (inherits(pc, "hd")) pc$nf else if (inherits(pc, "mvnfd")) pc else NULL
  list(psi  = if (!is.null(fpart)) fpart$eval(t_eval)[[1]],
       coef = if (!is.null(fpart)) as.matrix(fpart$coefs[[1]]),
       v    = if (!is.null(vpart)) t(vpart$data[[1]]),
       u    = fit$lsv)
}

# Estimated components are matched to the true ones by the permutation that
# maximizes the total absolute cosine similarity in the hybrid space (only the
# parts a method estimates are used). Components 2 and 3 have similar
# variances, so the sequential algorithm can extract them in either order.
perms3 <- list(c(1, 2, 3), c(1, 3, 2), c(2, 1, 3), c(2, 3, 1), c(3, 1, 2), c(3, 2, 1))
match_components <- function(ld, psi_ref, V_ref, t_eval) {
  sim <- matrix(0, 3, 3)
  for (k in 1:3) for (m in 1:3) {
    num <- den_est <- den_ref <- 0
    if (!is.null(ld$psi)) {
      num     <- num + trapz_ip(ld$psi[, k], psi_ref[, m], t_eval)
      den_est <- den_est + trapz_ip(ld$psi[, k], ld$psi[, k], t_eval)
      den_ref <- den_ref + trapz_ip(psi_ref[, m], psi_ref[, m], t_eval)
    }
    if (!is.null(ld$v)) {
      num     <- num + sum(ld$v[, k] * V_ref[, m])
      den_est <- den_est + sum(ld$v[, k]^2)
      den_ref <- den_ref + sum(V_ref[, m]^2)
    }
    sim[k, m] <- abs(num) / sqrt(den_est * den_ref)
  }
  perms3[[which.max(vapply(perms3, function(p) sum(sim[cbind(1:3, p)]), numeric(1)))]]
}

# Hybrid sign of estimated component k relative to true component m.
hybrid_sign <- function(ld, k, m, psi_ref, V_ref, t_eval) {
  s <- 0
  if (!is.null(ld$psi)) s <- s + trapz_ip(ld$psi[, k], psi_ref[, m], t_eval)
  if (!is.null(ld$v)) s <- s + sum(ld$v[, k] * V_ref[, m])
  if (s < 0) -1 else 1
}

# Scale an estimate to unit norm and align its sign with the reference.
aligned_unit <- function(est, ref, ip) {
  est <- est / sqrt(ip(est, est))
  if (ip(est, ref) < 0) -est else est
}

# Error measures of one fit, per true component m (Section 6): Err(psi_m) and
# Err(v_m) after unit scaling and sign alignment, correctly identified zeros
# (CIZ), misclassification rate (MCR) of the vector loadings, the share of
# zero scores within each group, and the selected tuning parameters.
evaluate_method <- function(name, fit, dat, cfg) {
  ld <- loadings_of(fit, cfg$t_eval)
  match <- match_components(ld, cfg$psi_eval, dat$V, cfg$t_eval)
  ip_f <- function(a, b) trapz_ip(a, b, cfg$t_eval)
  ip_v <- function(a, b) sum(a * b)
  pick <- function(x, k) {
    if (length(x) < k || is.null(x[[k]])) return(NA_real_)
    as.numeric(unlist(x[[k]]))[1]
  }
  metrics <- curves <- vsel <- NULL
  for (k in 1:3) {
    m <- match[k]
    err_psi <- n_zero_coef <- NA_real_
    if (!is.null(ld$psi)) {
      f <- aligned_unit(ld$psi[, k], cfg$psi_eval[, m], ip_f)
      err_psi <- sqrt(ip_f(f - cfg$psi_eval[, m], f - cfg$psi_eval[, m]))
      cf <- ld$coef[, k]
      n_zero_coef <- sum(abs(cf) <= 1e-8 * max(abs(cf)))
      curves <- rbind(curves, data.frame(method = name, comp = m, t = cfg$t_eval[cfg$keep],
                                         value = f[cfg$keep]))
    }
    err_v <- ciz <- mcr <- n_zero_v <- NA_real_
    if (!is.null(ld$v)) {
      g <- aligned_unit(ld$v[, k], dat$V[, m], ip_v)
      err_v <- sqrt(sum((g - dat$V[, m])^2))
      est_zero <- ld$v[, k] == 0
      true_zero <- dat$V[, m] == 0
      n_zero_v <- sum(est_zero)
      ciz <- sum(est_zero & true_zero)
      mcr <- mean(est_zero != true_zero)
      vsel <- rbind(vsel, data.frame(method = name, comp = m, j = seq_along(est_zero),
                                     nonzero = as.integer(!est_zero)))
    }
    zero_by_group <- tapply(ld$u[, k] == 0, dat$group, mean)
    alpha <- pick(fit$smooth_tuning, k)
    metrics <- rbind(metrics, data.frame(
      method = name, k = k, comp = m, err_psi = err_psi, err_v = err_v,
      ciz = ciz, n_true_zero_v = sum(dat$V[, m] == 0), mcr = mcr,
      n_zero_v = n_zero_v, n_zero_coef = n_zero_coef,
      zero_S1 = zero_by_group[[1]], zero_S2 = zero_by_group[[2]], zero_S3 = zero_by_group[[3]],
      sel_u = pick(fit$sparse_tuning_u, k), sel_fd = pick(fit$sparse_tuning_fd, k),
      sel_nfd = pick(fit$sparse_tuning_nfd, k),
      log2_alpha = if (is.na(alpha) || alpha <= 0) NA_real_ else log2(alpha)
    ))
  }
  list(metrics = metrics, curves = curves, vsel = vsel, loadings = ld, match = match)
}

mat_sqrt <- function(G) {
  e <- eigen(G, symmetric = TRUE)
  e$vectors %*% (sqrt(pmax(e$values, 0)) * t(e$vectors))
}

# The centred data in the coordinates of Section 4, Z = [X, C G^{1/2}], in
# which the hybrid inner product is the Euclidean one.
hybrid_matrix <- function(dat, mf, basis) {
  C <- t(matrix(mf$coefs, nrow = basis$nbasis))
  C <- sweep(C, 2, colMeans(C))
  X <- sweep(dat$X, 2, colMeans(dat$X))
  G_half <- mat_sqrt(fda::inprod(basis, basis))
  list(Z = cbind(X, C %*% G_half), G_half = G_half)
}

# Cumulative percentage of explained variance (Section 5) of the first r
# components: "adjusted" projects the data on the span of the loadings
# (||Z^(r)||^2 / ||Z||^2); "naive" adds up the variances of the individual
# projections, which can overstate it for non-orthogonal loadings.
cpev_curve <- function(ld, Z, G_half) {
  Vs <- rbind(ld$v, G_half %*% ld$coef)
  total <- sum(Z^2)
  adjusted <- vapply(1:3, function(r) {
    Vr <- Vs[, seq_len(r), drop = FALSE]
    sum((Z %*% Vr %*% solve(crossprod(Vr), t(Vr)))^2) / total
  }, numeric(1))
  Vn <- sweep(Vs, 2, sqrt(colSums(Vs^2)), "/")
  data.frame(r = 1:3, adjusted = adjusted, naive = cumsum(colSums((Z %*% Vn)^2)) / total)
}

# Simulate the data set of one task and set up the fits.
setup_task <- function(task, cfg) {
  set.seed(task$seed)
  dat <- simulate_hybrid(cfg$N, cfg$W, task$sigma, cfg$scenarios[[task$scenario]],
                         cfg$kappa, cfg$psi_obs)
  basis <- fda::create.bspline.basis(c(0, 1), cfg$n_basis)
  mf <- Mfd(argval = cfg$t_grid, X = dat$Y, mdbs = Basismfd(basis))
  list(dat = dat, basis = basis, mf = mf, specs = method_specs(mf, nfd(dat$X), cfg))
}

# Fit one method. Each method gets its own seed for its CV folds, so a fit
# does not depend on the other fits and can be reproduced on its own.
fit_method <- function(setup, name, task) {
  i <- match(name, names(setup$specs))
  set.seed(task$seed * 10 + i)
  run_quietly(do.call(Mhpca, c(setup$specs[[i]], list(ncomp = 3))))
}

# One data set: simulate, fit all methods, evaluate.
run_task <- function(task, cfg) {
  s <- setup_task(task, cfg)
  zinfo <- hybrid_matrix(s$dat, s$mf, s$basis)
  id <- data.frame(scenario = task$scenario, sigma = task$sigma, rep = task$rep)
  out <- list(metrics = NULL, curves = NULL, vsel = NULL, cpev = NULL, errors = NULL)
  for (name in names(s$specs)) {
    res <- fit_method(s, name, task)
    if (is.null(res$value)) {
      out$errors <- rbind(out$errors, cbind(id, method = name, error = res$error))
      next
    }
    ev <- evaluate_method(name, res$value, s$dat, cfg)
    out$metrics <- rbind(out$metrics, cbind(id, ev$metrics))
    if (!is.null(ev$curves)) out$curves <- rbind(out$curves, cbind(id, ev$curves))
    if (name %in% c("ReHPCA", "sPCA")) out$vsel <- rbind(out$vsel, cbind(id, ev$vsel))
    if (name %in% c("HPCA", "ReHPCA")) {
      out$cpev <- rbind(out$cpev, cbind(id, method = name,
                                        cpev_curve(ev$loadings, zinfo$Z, zinfo$G_half)))
    }
  }
  out
}

sim_helpers <- c("trapz_ip", "psi_raw", "simulate_hybrid", "run_quietly", "method_specs",
                 "loadings_of", "perms3", "match_components", "aligned_unit",
                 "evaluate_method", "mat_sqrt", "hybrid_matrix", "cpev_curve",
                 "setup_task", "fit_method", "run_task")


## 5. Run the simulation ----------------------------------------------------------

design <- expand.grid(rep = seq_len(n_rep), sigma = sigma_levels,
                      scenario = names(scenarios), stringsAsFactors = FALSE)
# Seeds depend only on (scenario, sigma, replicate), so adding replicates or
# changing the number of workers leaves existing results unchanged.
design$seed <- seed + 1e5 * match(design$scenario, names(scenarios)) +
  1e3 * match(design$sigma, sigma_levels) + design$rep

sim_cfg <- list(N = N, W = W, kappa = kappa, scenarios = scenarios, n_basis = n_basis,
                t_grid = t_grid, psi_obs = psi_obs, t_eval = t_eval, psi_eval = psi_eval,
                keep = keep_eval, alpha_grid = alpha_grid, grid_u = grid_u,
                grid_fd = grid_fd, grid_nfd = grid_nfd, nfold = nfold)

run_simulation <- function(design, cfg, n_cores) {
  tasks <- split(design, seq_len(nrow(design)))
  n_cores <- min(n_cores, length(tasks))
  message(sprintf("Simulation: %d data sets x 6 fits on %d worker process(es)...",
                  length(tasks), n_cores))
  start <- Sys.time()
  if (n_cores > 1) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, load_mhpca, pkg_root, use_load_all)
    parallel::clusterExport(cl, sim_helpers, envir = globalenv())
    out <- parallel::parLapplyLB(cl, tasks, run_task, cfg = cfg)
  } else {
    out <- lapply(tasks, run_task, cfg = cfg)
  }
  minutes <- as.numeric(difftime(Sys.time(), start, units = "mins"))
  message(sprintf("Simulation finished in %.1f minutes.", minutes))
  stack <- function(name) do.call(rbind, lapply(out, `[[`, name))
  list(metrics = stack("metrics"), curves = stack("curves"), vsel = stack("vsel"),
       cpev = stack("cpev"), errors = stack("errors"), minutes = minutes)
}

sim <- NULL
cache_key <- list(sim_cfg, design)
if (!is.null(sim_cache) && file.exists(sim_cache)) {
  cached <- readRDS(sim_cache)
  if (identical(cached$key, cache_key)) {
    sim <- cached$sim
    message("Using cached simulation results from ", sim_cache)
  } else {
    message("Cached results do not match the current settings; running the simulation.")
  }
}
if (is.null(sim)) {
  sim <- run_simulation(design, sim_cfg, n_cores)
  if (!is.null(sim_cache)) saveRDS(list(key = cache_key, sim = sim), sim_cache)
}
if (!is.null(sim$errors) && nrow(sim$errors)) {
  warning(sprintf("%d fit(s) failed and are left out of the figures, e.g.: %s",
                  nrow(sim$errors), sim$errors$error[1]), call. = FALSE)
}


## 6. Figures ---------------------------------------------------------------------

M <- sim$metrics
M$scenario <- factor(M$scenario, levels = names(scenarios))

# Summaries by group: `fun` maps a data frame to a one-row data frame.
summarise_by <- function(df, by, fun) {
  pieces <- split(df, interaction(df[by], drop = TRUE, lex.order = TRUE))
  out <- do.call(rbind, lapply(pieces, function(d) cbind(d[1, by, drop = FALSE], fun(d))))
  rownames(out) <- NULL
  out
}

# Example data set for Figures 1, 2, 6 and 7: the replicate of scenario II at
# sigma_show whose ReHPCA error (mean Err(psi_m) over the components) is
# closest to the median. It is simulated and fitted again with the seeds of
# the simulation, so it is exactly that replicate's fit.
example_fit <- function(task, cfg) {
  s <- setup_task(task, cfg)
  fit <- fit_method(s, "ReHPCA", task)$value
  ev <- evaluate_method("ReHPCA", fit, s$dat, cfg)
  list(rep = task$rep, data = s$dat, match = ev$match, loadings = ev$loadings,
       cv = list(u = fit$CVs_u, fd = fit$CVs_fd, nfd = fit$CVs_nfd), gcv = fit$GCVs,
       selected = list(u = fit$sparse_tuning_u, fd = fit$sparse_tuning_fd,
                       nfd = fit$sparse_tuning_nfd, alpha = fit$smooth_tuning))
}
re_ii <- subset(M, method == "ReHPCA" & scenario == "II" & sigma == sigma_show)
rep_err <- tapply(re_ii$err_psi, re_ii$rep, mean)
ex_rep <- as.integer(names(rep_err)[which.min(abs(rep_err - stats::median(rep_err)))])
message(sprintf("Example data set: scenario II, sigma = %s, replicate %d (refitting ReHPCA)...",
                sigma_show, ex_rep))
ex <- example_fit(design[design$scenario == "II" & design$sigma == sigma_show &
                           design$rep == ex_rep, ], sim_cfg)
ex_basis <- fda::create.bspline.basis(c(0, 1), n_basis)
fd_true_zeros <- true_zero_coefs(ex_basis)          # per component
err_floor <- basis_floor(ex_basis)                  # per component
b_size <- W / 3
group_sizes <- function(s) {
  sizes <- round(N * scenarios[[s]])
  sizes[3] <- N - sum(sizes[1:2])
  sizes
}

# Aligned (unit norm, true sign) loadings and scores of the example ReHPCA fit.
ex_aligned <- lapply(1:3, function(k) {
  m <- ex$match[k]
  ld <- ex$loadings
  s <- hybrid_sign(ld, k, m, psi_eval, ex$data$V, t_eval)
  list(k = k, comp = m,
       psi = s * ld$psi[, k] / sqrt(trapz_ip(ld$psi[, k], ld$psi[, k], t_eval)),
       v = s * ld$v[, k] / sqrt(sum(ld$v[, k]^2)),
       u = s * ld$u[, k])
})


### Figure 1: schematic of the method -------------------------------------------

# Small plots that sit inside the schematic (real data and estimates from the
# example data set).
ex_group <- ex$data$group
show_units <- unlist(lapply(1:3, function(g) head(which(ex_group == g), 6)))
mini_curves <- data.frame(t = rep(t_grid, length(show_units)),
                          value = c(ex$data$Y[, show_units]),
                          unit = rep(show_units, each = n_grid),
                          group = factor(rep(ex_group[show_units], each = n_grid)))
ins_curves <- ggplot(mini_curves, aes(t, value, group = unit, colour = group)) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.2) +
  geom_line(linewidth = 0.25) +
  scale_colour_manual(values = pal_comp, guide = "none") +
  theme_void()

ord <- order(ex_group)
rows_shown <- ord[seq(1, N, by = 4)]
mini_X <- data.frame(j = rep(seq_len(W), each = length(rows_shown)),
                     row = rep(seq_along(rows_shown), W),
                     value = c(ex$data$X[rows_shown, ]))
lim_X <- stats::quantile(abs(mini_X$value), 0.98)
ins_X <- ggplot(mini_X, aes(j, row, fill = pmax(pmin(value, lim_X), -lim_X))) +
  geom_raster() +
  scale_fill_div(lim_X, guide = "none") +
  scale_y_reverse(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(panel.border = element_rect(fill = NA, colour = "grey40", linewidth = 0.3))

k2 <- which(ex$match == 2)  # estimated component matched to component 2
a2 <- ex_aligned[[k2]]
ins_psi <- ggplot(data.frame(t = t_eval, value = a2$psi), aes(t, value)) +
  annotate("rect", xmin = 1 / 3, xmax = 2 / 3, ymin = -Inf, ymax = Inf, fill = "grey93") +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.25) +
  geom_line(colour = pal_method[["ReHPCA"]], linewidth = 0.5) +
  theme_void()
ins_v <- ggplot(data.frame(j = seq_len(W), value = a2$v), aes(j, value)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.25) +
  geom_segment(aes(xend = j, y = 0, yend = value), colour = pal_method[["ReHPCA"]],
               linewidth = 0.4) +
  geom_point(data = function(d) d[d$value != 0, ], colour = pal_method[["ReHPCA"]],
             size = 0.5) +
  geom_point(data = function(d) d[d$value == 0, ], colour = "grey55", size = 0.35) +
  theme_void()
ins_u <- ggplot(data.frame(i = seq_len(N), value = a2$u[ord]), aes(i, value)) +
  geom_segment(aes(xend = i, y = 0, yend = value), colour = pal_method[["ReHPCA"]],
               linewidth = 0.15) +
  geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.25) +
  theme_void()

# Canvas in "schematic units": x in [0, 100], y in [0, 50].
xr <- c(0, 100); yr <- c(0, 50)
stage <- data.frame(x0 = c(0.4, 23.2, 48.6, 80.8), x1 = c(19.8, 45.2, 77.4, 99.6))
stage$mid <- (stage$x0 + stage$x1) / 2
stage$title <- c("Hybrid data", "Basis and half-smoothing", "Penalized power algorithm",
                 "Regularized PCs")
zx <- stage$x0[2] + 4.4   # left edge of the sketch of Z = [X  C~] in stage 2
top <- 49.4; bottom <- 0.6; head_h <- 4.6
math <- function(x, y, label, size = 6.8, hjust = 0.5, colour = "grey10") {
  annotate("text", x = x, y = y, label = label, parse = TRUE, size = pt(size),
           hjust = hjust, colour = colour)
}
note <- function(x, y, label, size = 6, hjust = 0.5, colour = "grey35", fontface = "italic") {
  annotate("text", x = x, y = y, label = label, size = pt(size), hjust = hjust,
           colour = colour, fontface = fontface)
}
arrow_spec <- arrow(length = unit(3.5, "pt"), type = "closed", angle = 25)
dbar <- intToUtf8(0x2016)  # double vertical bar for norms (plotmath's '||' draws one bar)

canvas <- ggplot() +
  coord_cartesian(xlim = xr, ylim = yr, expand = FALSE) +
  # stage boxes and header bands
  annotate("rect", xmin = stage$x0, xmax = stage$x1, ymin = bottom, ymax = top,
           fill = "white", colour = "grey55", linewidth = 0.35) +
  annotate("rect", xmin = stage$x0, xmax = stage$x1, ymin = top - head_h, ymax = top,
           fill = "#E7EEF6", colour = "grey55", linewidth = 0.35) +
  annotate("point", x = stage$x0 + 2.1, y = top - head_h / 2, shape = 21, size = 3.6,
           fill = pal_method[["ReHPCA"]], colour = pal_method[["ReHPCA"]], stroke = 0.3) +
  annotate("text", x = stage$x0 + 2.1, y = top - head_h / 2, label = 1:4, size = pt(6.5),
           colour = "white", fontface = "bold") +
  annotate("text", x = stage$x0 + 4.1, y = top - head_h / 2, label = stage$title,
           size = pt(7.2), hjust = 0, fontface = "bold", colour = "grey10") +
  # arrows between stages
  annotate("segment", x = stage$x1[1:3] + 0.5, xend = stage$x0[2:4] - 0.5,
           y = 25, yend = 25, colour = "grey30", linewidth = 0.45, arrow = arrow_spec) +

  # 1. data
  math(stage$mid[1], 41.2, "bold(z)[i] == group('(', list(bold(x)[i], y[i]), ')')*','~~i == 1*','*ldots*','*n") +
  math(stage$mid[1], 24.3, "'functional part'~~y[i]*(t)", size = 6.3) +
  math(stage$mid[1], 4.0, "'vector part'~~bold(x)[i] %in% R^m", size = 6.3) +

  # 2. basis expansion, Z = [X  C~], half-smoothing
  math(stage$x0[2] + 1.2, 41.2, "y[ij]*(t) == nu[j]*(t)^T*bold(c)[ij]", hjust = 0) +
  math(stage$x0[2] + 1.2, 37.6, "tilde(bold(C)) == bold(C)*bold(G)^{1/2}*','~~tilde(bold(b)) == bold(G)^{1/2}*bold(b)", hjust = 0) +
  math(zx - 1.6, 27.0, "bold(Z) == ''", hjust = 0.5) +
  annotate("rect", xmin = zx, xmax = zx + 6, ymin = 21.0, ymax = 33.0,
           fill = "#EDEDED", colour = "grey40", linewidth = 0.3) +
  annotate("rect", xmin = zx + 6, xmax = zx + 15.6, ymin = 21.0, ymax = 33.0,
           fill = "#DDEAF6", colour = "grey40", linewidth = 0.3) +
  math(zx + 3, 28.4, "bold(X)", size = 7.5) +
  math(zx + 3, 25.0, "n %*% m", size = 6) +
  math(zx + 10.8, 28.4, "tilde(bold(C))", size = 7.5) +
  math(zx + 10.8, 25.0, "n %*% d", size = 6) +
  note(zx + 10.8, 19.6, "smoothed", size = 5.8) +
  note(zx + 3, 19.6, "not smoothed", size = 5.8) +
  math(stage$x0[2] + 1.2, 14.6, "tilde(bold(S))[alpha] == bold(G)^{1/2}*(bold(G) + bold(D)[alpha])^{-1}*bold(G)^{1/2}", hjust = 0, size = 6.5) +
  math(stage$x0[2] + 1.2, 9.6, sprintf("min~~'%s'*bold(Z) - bold(u)*tilde(bold(v))^T*'%s'[F]^2 + bold(u)^T*bold(u)~tilde(bold(b))^T*Omega[alpha]*tilde(bold(b))", dbar, dbar), hjust = 0, size = 6.3) +
  note(stage$x0[2] + 1.2, 5.6, "+ sparsity penalties on u, b and v  (3.2)", hjust = 0, size = 5.8) +

  # 3. algorithm
  annotate("rect", xmin = stage$x0[3] + 1.0, xmax = stage$x1[3] - 2.6, ymin = 19.2, ymax = 42.2,
           fill = "#FAFAFA", colour = "grey60", linewidth = 0.3) +
  note(stage$x0[3] + 2.0, 40.4, "repeat until convergence", hjust = 0) +
  math(stage$x0[3] + 2.0, 36.4, "(a)~~bold(u) %<-% h[gamma[1]]*group('(', tilde(bold(C))*tilde(bold(b)) + bold(X)*bold(v), ')')", hjust = 0, size = 6.5) +
  math(stage$x0[3] + 2.0, 31.6, "(b)~~tilde(bold(b)) %<-% bold(G)^{1/2}*h[gamma[2]]*group('(', bold(G)^{-1/2}*tilde(bold(S))[alpha]*tilde(bold(C))^T*bold(u), ')')", hjust = 0, size = 6.5) +
  math(stage$x0[3] + 2.0, 26.8, "(c)~~bold(v) %<-% h[gamma[3]]*group('(', bold(X)^T*bold(u), ')')", hjust = 0, size = 6.5) +
  math(stage$x0[3] + 2.0, 22.4, "(d)~~'normalize'~~bold(u)*','~~(bold(v)*','~bold(G)^{-1/2}*tilde(bold(b)))", hjust = 0, size = 6.5) +
  annotate("curve", x = stage$x1[3] - 2.4, y = 22.4, xend = stage$x1[3] - 2.4, yend = 36.4,
           curvature = 0.55, colour = "grey30", linewidth = 0.4, arrow = arrow_spec) +
  math(stage$x0[3] + 1.4, 15.6, "'deflate:'~~tilde(bold(C)) %<-% tilde(bold(C)) - bold(u)*tilde(bold(b))^T*','~~bold(X) %<-% bold(X) - bold(u)*bold(v)^T", hjust = 0, size = 6.3) +
  note(stage$x0[3] + 1.4, 12.6, "then the next component (separate tuning per PC)", hjust = 0, size = 5.8) +
  annotate("rect", xmin = stage$x0[3] + 1.0, xmax = stage$x1[3] - 1.0, ymin = 1.6, ymax = 10.0,
           fill = "#F3F7FB", colour = "grey70", linewidth = 0.3) +
  annotate("text", x = stage$x0[3] + 2.0, y = 8.2, label = "Tuning (Section 4.2)",
           size = pt(6.3), hjust = 0, fontface = "bold", colour = "grey10") +
  math(stage$x0[3] + 2.0, 5.5, "'1.'~~gamma[1]*','~gamma[2]*','~gamma[3]*':'~~'K-fold CV on the degree of sparsity'", hjust = 0, size = 6) +
  math(stage$x0[3] + 2.0, 3.0, "'2.'~~alpha*':'~~'closed-form GCV (4.5)'", hjust = 0, size = 6) +

  # 4. output
  math(stage$mid[4], 31.0, "'smooth, localized'~~hat(psi)[m]*(t)", size = 6.3) +
  math(stage$mid[4], 17.6, "'sparse vector loading'~~hat(bold(v))[m]", size = 6.3) +
  math(stage$mid[4], 4.0, "'sparse scores'~~hat(bold(u))[m]", size = 6.3) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

# Position of an inset given in schematic units.
inset <- function(p, x0, x1, y0, y1) {
  inset_element(p, left = (x0 - xr[1]) / diff(xr), right = (x1 - xr[1]) / diff(xr),
                bottom = (y0 - yr[1]) / diff(yr), top = (y1 - yr[1]) / diff(yr),
                align_to = "full", clip = TRUE)
}
p1 <- canvas +
  inset(ins_curves, stage$x0[1] + 1.2, stage$x1[1] - 1.2, 26.2, 38.6) +
  inset(ins_X, stage$x0[1] + 2.4, stage$x1[1] - 2.4, 6.4, 21.6) +
  inset(ins_psi, stage$x0[4] + 1.6, stage$x1[4] - 1.6, 33.0, 42.6) +
  inset(ins_v, stage$x0[4] + 1.6, stage$x1[4] - 1.6, 19.8, 28.8) +
  inset(ins_u, stage$x0[4] + 1.6, stage$x1[4] - 1.6, 6.4, 15.6)


### Figure 2: simulation design ---------------------------------------------------

tA <- seq(0, 1, length.out = 1201)
dA <- data.frame(t = rep(tA, 3), comp = factor(rep(1:3, each = length(tA))),
                 value = c(psi_true(tA)), seg = rep(interval_of(tA), 3))
# draw each psi_m only on the subintervals where it is nonzero
support <- list(`1` = c(1L, 3L), `2` = 2L, `3` = 3L)
dA <- dA[mapply(function(m, s) s %in% support[[m]], as.character(dA$comp), dA$seg), ]
labA <- data.frame(t = c(0.25, 0.43, 0.84), y = c(2.18, 2.45, -2.45),
                   lab = c("psi[1]", "psi[2]", "psi[3]"))
pA <- ggplot(dA, aes(t, value)) +
  annotate("rect", xmin = 1 / 3, xmax = 2 / 3, ymin = -Inf, ymax = Inf, fill = "grey95") +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.25) +
  geom_line(aes(colour = comp, group = interaction(comp, seg)), linewidth = 0.6) +
  annotate("text", x = c(1, 3, 5) / 6, y = 3.05, size = pt(6.8), colour = "grey30",
           label = c("italic(I)[1]", "italic(I)[2]", "italic(I)[3]"), parse = TRUE) +
  geom_text(data = labA, aes(t, y, label = lab), parse = TRUE, size = pt(7.5),
            colour = "grey10") +
  scale_colour_manual(values = pal_comp, guide = "none") +
  scale_x_continuous(breaks = c(0, 1 / 3, 2 / 3, 1), labels = c("0", "1/3", "2/3", "1"),
                     expand = expansion(mult = 0.01)) +
  scale_y_continuous(limits = c(-3.05, 3.35)) +
  labs(x = expression(t), y = expression(psi[m]*(t)), title = "True functional loadings") +
  theme_paper()

dB <- data.frame(j = rep(seq_len(W), 3), comp = factor(rep(1:3, each = W)),
                 value = c(ex$data$V))
dB <- dB[dB$value != 0, ]
pB <- ggplot(dB, aes(j, value, colour = comp)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_segment(aes(xend = j, y = 0, yend = value), linewidth = 0.45) +
  geom_point(size = 0.9) +
  annotate("text", x = (0:2) * b_size + (b_size + 1) / 2, y = max(dB$value) + 0.07,
           label = sprintf("bold(v)[%d]", 1:3), parse = TRUE, size = pt(7.5),
           colour = "grey10") +
  scale_colour_manual(values = pal_comp, guide = "none") +
  scale_x_continuous(breaks = c(1, seq(5, W, 5)), expand = expansion(add = 0.8)) +
  scale_y_continuous(limits = c(0, max(dB$value) + 0.11), expand = expansion(mult = c(0, 0.02))) +
  labs(x = expression("Vector variable"~italic(j)), y = expression(v[mj]),
       title = "True vector loadings (zero outside their block)") +
  theme_paper()

lam <- rbind(c(0, 0.5, 0.9), c(0.9, 0, 0.5), c(0.5, 0.9, 0))
dC <- do.call(rbind, lapply(names(scenarios), function(s) {
  edges <- c(0, cumsum(scenarios[[s]]))
  g <- expand.grid(group = 1:3, comp = 1:3)
  data.frame(scenario = s, group = g$group, comp = g$comp, xmin = edges[g$group],
             xmax = edges[g$group + 1], lambda = lam[cbind(g$group, g$comp)])
}))
dC$scen_lab <- f_scen(dC$scenario, long = TRUE)
dC$label <- ifelse(dC$lambda == 0, "0", format(dC$lambda))
dC$y <- 4 - dC$comp   # component 1 on top
labC <- unique(dC[, c("scenario", "scen_lab", "group", "xmin", "xmax")])
pC <- ggplot(dC) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = y - 0.42, ymax = y + 0.42,
                fill = factor(comp), alpha = lambda / 0.9),
            colour = "grey60", linewidth = 0.25) +
  geom_text(aes(x = (xmin + xmax) / 2, y = y, label = label), size = pt(6.2),
            colour = "grey10") +
  geom_text(data = labC, aes(x = (xmin + xmax) / 2, y = 3.74, label = sprintf("S[%d]", group)),
            parse = TRUE, size = pt(6.5), colour = "grey10") +
  facet_wrap(~scen_lab, ncol = 1, labeller = label_parsed) +
  scale_fill_manual(values = pal_comp, guide = "none") +
  scale_alpha_identity() +
  scale_y_continuous(breaks = 3:1, labels = parse(text = sprintf("m == %d", 1:3)),
                     limits = c(0.52, 4.02), expand = c(0, 0)) +
  scale_x_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(x = "Units, ordered by group (cumulative proportion)", y = NULL,
       title = expression(bold("Score variances")~lambda[m]^{(k)}~bold("of component")~italic(m)~bold("in group")~S[k]~bold("(0 = component absent)"))) +
  theme_paper() +
  theme(panel.grid.major = element_blank(), panel.spacing = unit(3, "pt"))

d_units <- unlist(lapply(1:3, function(g) head(which(ex_group == g), 8)))
dD <- data.frame(t = rep(t_grid, length(d_units)), value = c(ex$data$Y[, d_units]),
                 unit = rep(d_units, each = n_grid),
                 group = factor(rep(ex_group[d_units], each = n_grid)))
pD <- ggplot(dD, aes(t, value, group = unit, colour = group)) +
  annotate("rect", xmin = 1 / 3, xmax = 2 / 3, ymin = -Inf, ymax = Inf, fill = "grey95") +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.25) +
  geom_line(linewidth = 0.3) +
  scale_colour_manual(values = pal_comp, labels = s_labels, name = NULL) +
  scale_x_continuous(breaks = c(0, 1 / 3, 2 / 3, 1), labels = c("0", "1/3", "2/3", "1"),
                     expand = expansion(mult = 0.01)) +
  labs(x = expression(t), y = expression(y[i]*(t)),
       title = bquote(bold("Curves, 8 units per group (scenario II,")~sigma == .(sigma_show)*bold(")"))) +
  theme_paper() +
  theme(legend.key.width = unit(12, "pt"))

ordE <- order(ex_group)
dE <- data.frame(row = rep(seq_len(N), W), j = rep(seq_len(W), each = N),
                 value = c(ex$data$X[ordE, ]))
lim_E <- stats::quantile(abs(dE$value), 0.99)
cuts_E <- cumsum(table(ex_group))
pE <- ggplot(dE, aes(j, row, fill = pmax(pmin(value, lim_E), -lim_E))) +
  geom_raster() +
  geom_hline(yintercept = cuts_E[1:2] + 0.5, colour = "grey10", linewidth = 0.3) +
  geom_vline(xintercept = c(b_size, 2 * b_size) + 0.5, colour = "grey10", linewidth = 0.3) +
  scale_y_reverse(breaks = (c(0, cuts_E[1:2]) + cuts_E) / 2, labels = s_labels,
                  expand = c(0, 0)) +
  scale_x_continuous(breaks = c(1, seq(5, W, 5)), expand = c(0, 0)) +
  scale_fill_div(lim_E, name = expression(x[ij]),
                 breaks = c(-1, 0, 1) * signif(0.75 * lim_E, 1)) +
  labs(x = expression("Vector variable"~italic(j)), y = NULL,
       title = "Vector part, all units") +
  theme_paper() +
  theme(legend.position = "right", legend.key.width = unit(6, "pt"),
        legend.key.height = unit(16, "pt"), panel.grid.major = element_blank())

p2 <- ((pA | pB) + plot_layout(widths = c(1.25, 1))) /
  pC /
  ((pD | pE) + plot_layout(widths = c(1.25, 1))) +
  plot_layout(heights = c(1, 1.3, 1.05)) +
  plot_annotation(tag_levels = "A")


### Figure 3: estimated vs true functional loadings ---------------------------------

cur3 <- subset(sim$curves, sigma == sigma_show & method %in% c("ReHPCA", "HPCA"))
band3 <- summarise_by(cur3, c("method", "scenario", "comp", "t"), function(d) {
  q <- stats::quantile(d$value, c(0.05, 0.5, 0.95), names = FALSE)
  data.frame(lo = q[1], mid = q[2], hi = q[3])
})
band3$method <- factor(band3$method, levels = c("ReHPCA", "HPCA"))
band3$scen_lab <- f_scen(band3$scenario)
band3$comp_lab <- factor(band3$comp, levels = 1:3, labels = sprintf("psi[%d]", 1:3))
truth3 <- data.frame(t = rep(tA, 3), value = c(psi_true(tA)),
                     seg = rep(interval_of(tA), 3), comp = rep(1:3, each = length(tA)))
truth3$comp_lab <- factor(truth3$comp, levels = 1:3, labels = sprintf("psi[%d]", 1:3))
p3 <- ggplot() +
  geom_vline(xintercept = c(1 / 3, 2 / 3), colour = "grey85", linewidth = 0.25) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.25) +
  geom_ribbon(data = band3, aes(t, ymin = lo, ymax = hi, fill = method), alpha = 0.25) +
  geom_line(data = band3, aes(t, mid, colour = method), linewidth = 0.5) +
  geom_line(data = truth3, aes(t, value, group = seg, linetype = "truth"),
            colour = col_truth, linewidth = 0.4) +
  facet_grid(scen_lab ~ comp_lab, labeller = label_parsed) +
  scale_colour_manual(NULL, values = pal_method, labels = lab_method[c("ReHPCA", "HPCA")],
                      breaks = c("ReHPCA", "HPCA"), aesthetics = c("colour", "fill")) +
  scale_linetype_manual(NULL, values = c(truth = "22"),
                        labels = expression("True"~psi[m])) +
  scale_x_continuous(breaks = c(0, 1 / 3, 2 / 3, 1), labels = c("0", "1/3", "2/3", "1")) +
  labs(x = expression(t), y = expression(hat(psi)[m]*(t))) +
  theme_paper() +
  theme(legend.key.width = unit(16, "pt"))


### Figure 4: estimation error of the functional loadings ---------------------------

d4 <- subset(M, method %in% c("ReHPCA", "HPCA", "ReFPCA", "FPCA") & !is.na(err_psi))
d4$method <- factor(d4$method, levels = c("ReHPCA", "HPCA", "ReFPCA", "FPCA"))
d4$sigma_lab <- f_sigma(d4$sigma)
d4$comp_lab <- f_comp(d4$comp, "psi")
log_breaks <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2)
breaks4 <- c(0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 1)
floor4 <- data.frame(comp = 1:3, floor = err_floor)
floor4$comp_lab <- f_comp(floor4$comp, "psi")
p4 <- ggplot(d4, aes(scenario, err_psi, fill = method)) +
  geom_hline(data = floor4, aes(yintercept = floor, linetype = "floor"), colour = "grey25",
             linewidth = 0.35) +
  geom_boxplot(width = 0.74, position = position_dodge(width = 0.84), linewidth = 0.25,
               colour = "grey20", outlier.size = 0.45, outlier.stroke = 0.2) +
  facet_grid(sigma_lab ~ comp_lab, labeller = label_parsed) +
  scale_y_log10(breaks = breaks4, labels = format(breaks4, drop0trailing = TRUE, trim = TRUE)) +
  scale_fill_manual(NULL, values = pal_method, labels = lab_method) +
  scale_linetype_manual(NULL, values = c(floor = "22"),
                        labels = sprintf("smallest error possible with %d B-splines", n_basis)) +
  labs(x = "Score-sparsity scenario", y = expression(Err(hat(psi)[m]))) +
  theme_paper() +
  theme(legend.box = "vertical", legend.spacing.y = unit(1, "pt")) +
  guides(fill = guide_legend(nrow = 1, order = 1), linetype = guide_legend(order = 2))


### Figure 5: vector loadings -----------------------------------------------------

freq5 <- aggregate(nonzero ~ scenario + comp + j,
                   data = subset(sim$vsel, method == "ReHPCA" & sigma == sigma_show),
                   FUN = mean)
freq5$scen_lab <- f_scen(freq5$scenario)
blocks5 <- data.frame(xmin = (0:2) * b_size + 0.5, xmax = (1:3) * b_size + 0.5,
                      ymin = (1:3) - 0.5, ymax = (1:3) + 0.5)
p5a <- ggplot(freq5, aes(j, comp, fill = nonzero)) +
  geom_tile(width = 0.86, height = 0.8) +
  geom_rect(data = blocks5, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = NA, colour = "grey10", linewidth = 0.35) +
  facet_wrap(~scen_lab, nrow = 1, labeller = label_parsed) +
  scale_y_reverse(breaks = 1:3, labels = parse(text = sprintf("bold(v)[%d]", 1:3)),
                  expand = c(0, 0)) +
  scale_x_continuous(breaks = c(1, seq(10, W, 10)), expand = c(0, 0)) +
  scale_fill_seq(name = "Selection\nfrequency", limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(x = expression("Vector variable"~italic(j)), y = NULL,
       title = "ReHPCA: how often each loading is nonzero (black frames: true support)") +
  theme_paper() +
  theme(panel.grid.major = element_blank(), legend.position = "right",
        legend.key.width = unit(6, "pt"), legend.key.height = unit(12, "pt"))

d5 <- subset(M, method %in% c("ReHPCA", "HPCA", "sPCA", "PCA") & !is.na(err_v) &
               sigma == sigma_show)
d5$method <- factor(d5$method, levels = c("ReHPCA", "HPCA", "sPCA", "PCA"))
d5$comp_lab <- f_comp(d5$comp, "v")
p5b <- ggplot(d5, aes(scenario, err_v, fill = method)) +
  geom_boxplot(width = 0.74, position = position_dodge(width = 0.84), linewidth = 0.25,
               colour = "grey20", outlier.size = 0.45, outlier.stroke = 0.2) +
  facet_grid(. ~ comp_lab, labeller = label_parsed) +
  scale_y_log10(breaks = log_breaks, labels = format(log_breaks, drop0trailing = TRUE,
                                                      scientific = FALSE, trim = TRUE)) +
  scale_fill_manual(NULL, values = pal_method, labels = lab_method) +
  labs(x = "Score-sparsity scenario", y = expression(Err(hat(bold(v))[m])),
       title = "Estimation error of the vector loadings") +
  theme_paper() +
  guides(fill = guide_legend(nrow = 1))

p5 <- p5a / p5b + plot_layout(heights = c(1, 1.35)) + plot_annotation(tag_levels = "A")


### Figure 6: sparse scores --------------------------------------------------------

d6a <- do.call(rbind, lapply(ex_aligned, function(a) {
  data.frame(comp = a$comp, group = factor(ex_group), score = sqrt(N) * a$u, zero = a$u == 0)
}))
d6a$colour_key <- ifelse(d6a$zero, "zero", as.character(d6a$group))
d6a$comp_lab <- f_comp(d6a$comp, "u")
lab6a <- aggregate(zero ~ comp_lab + group, data = d6a, FUN = mean)
lab6a$label <- sprintf("%.0f%% zero", 100 * lab6a$zero)
p6a <- ggplot(d6a, aes(group, score)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.25) +
  geom_point(data = d6a[!d6a$zero, ], aes(colour = colour_key),
             position = position_jitter(width = 0.33, height = 0, seed = 1),
             size = 0.55, stroke = 0, alpha = 0.85) +
  geom_point(data = d6a[d6a$zero, ], aes(colour = colour_key),  # zeros drawn on top
             position = position_jitter(width = 0.33, height = 0, seed = 1),
             size = 0.6, stroke = 0) +
  geom_text(data = lab6a, aes(group, Inf, label = label), vjust = 1.6, size = pt(6.2),
            colour = "grey20") +
  facet_wrap(~comp_lab, nrow = 1, labeller = label_parsed) +
  scale_colour_manual(NULL, values = c(pal_comp, zero = "grey25"),
                      breaks = c("1", "2", "3", "zero"),
                      labels = c(s_labels, expression("score set to 0"))) +
  scale_x_discrete(labels = s_labels) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  labs(x = "Group", y = expression(sqrt(n)~hat(u)[im]),
       title = bquote(bold("ReHPCA scores in the example data set (scenario II,")~sigma == .(sigma_show)*bold(")"))) +
  theme_paper() +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 1.6, alpha = 1)))

d6 <- subset(M, method == "ReHPCA" & sigma == sigma_show)
d6l <- do.call(rbind, lapply(1:3, function(g) {
  data.frame(scenario = d6$scenario, comp = d6$comp, group = g, rate = d6[[paste0("zero_S", g)]])
}))
agg6 <- summarise_by(d6l, c("scenario", "comp", "group"), function(d) {
  data.frame(mean = mean(d$rate), sd = stats::sd(d$rate))
})
agg6$comp_lab <- f_comp(agg6$comp, "u")
agg6$group <- factor(agg6$group)
p6b <- ggplot(agg6, aes(scenario, 100 * mean, fill = group)) +
  geom_col(position = position_dodge(width = 0.82), width = 0.74) +
  geom_errorbar(aes(ymin = 100 * pmax(mean - sd, 0), ymax = 100 * pmin(mean + sd, 1)),
                position = position_dodge(width = 0.82), width = 0.25, linewidth = 0.3,
                colour = "grey20") +
  facet_wrap(~comp_lab, nrow = 1, labeller = label_parsed) +
  scale_fill_manual(values = pal_comp, labels = s_labels, name = "Units in") +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.03))) +
  labs(x = "Score-sparsity scenario", y = expression("% of units with"~hat(u)[im] == 0),
       title = sprintf("Across %d replicates (mean and SD)", n_rep)) +
  theme_paper()

p6 <- p6a / p6b + plot_layout(heights = c(1.1, 1)) + plot_annotation(tag_levels = "A")


### Figure 7: tuning on the example data set ---------------------------------------

cv_table <- function(x) if (is.matrix(x)) x else x[[1]]
pc_lab <- function(k) {
  factor(k, levels = 1:3, labels = sprintf("'PC %d (component %d)'", 1:3, ex$match))
}
true_count <- list(
  u = round(N * scenarios[["II"]])[ex$match],
  fd = fd_true_zeros[ex$match],
  nfd = rep(W - b_size, 3)
)
cv_row <- function(block, xlab, title) {
  tabs <- ex$cv[[block]]
  d <- do.call(rbind, lapply(1:3, function(k) {
    m <- cv_table(tabs[[k]])
    data.frame(k = k, gamma = m[, "gamma"], cv = m[, "CV"] / min(m[, "CV"]),
               se = m[, "SE"] / min(m[, "CV"]))
  }))
  sel <- vapply(1:3, function(k) as.numeric(unlist(ex$selected[[block]][[k]]))[1], numeric(1))
  marks <- do.call(rbind, lapply(1:3, function(k) {
    dk <- d[d$k == k, ]
    i <- which.min(dk$cv)
    data.frame(k = k, cutoff = dk$cv[i] + dk$se[i], selected = sel[k],
               cv_sel = dk$cv[match(sel[k], dk$gamma)], truth = true_count[[block]][k])
  }))
  d$pc <- pc_lab(d$k)
  marks$pc <- pc_lab(marks$k)
  ggplot(d, aes(gamma, cv)) +
    geom_hline(data = marks, aes(yintercept = cutoff), colour = "grey45", linewidth = 0.3,
               linetype = "dotted") +
    geom_vline(data = marks, aes(xintercept = truth, linetype = "truth"), colour = col_truth,
               linewidth = 0.35) +
    geom_errorbar(aes(ymin = cv - se, ymax = cv + se), width = 0, colour = "grey65",
                  linewidth = 0.3) +
    geom_line(colour = "grey35", linewidth = 0.35) +
    geom_point(colour = "grey35", size = 0.6) +
    geom_point(data = marks, aes(selected, cv_sel, shape = "selected"),
               colour = pal_method[["ReHPCA"]], size = 1.9) +
    facet_wrap(~pc, nrow = 1, scales = "free_y", labeller = label_parsed) +
    scale_linetype_manual(NULL, values = c(truth = "22"), labels = "true number of zeros") +
    scale_shape_manual(NULL, values = c(selected = 16), labels = "selected value") +
    labs(x = xlab, y = "CV error / min", title = title) +
    theme_paper()
}
p7u <- cv_row("u", "Number of zero scores", "Scores: sparsity of u by CV over the columns of Z")
p7f <- cv_row("fd", "Number of zero B-spline coefficients",
              "Functional loading: sparsity of b by CV over the rows")
p7v <- cv_row("nfd", "Number of zero vector loadings",
              "Vector loading: sparsity of v by CV over the rows")
d7g <- do.call(rbind, lapply(1:3, function(k) {
  g <- ex$gcv[[k]]
  data.frame(k = k, log2a = log2(alpha_grid), gcv = g / min(g))
}))
d7g$pc <- pc_lab(d7g$k)
sel7g <- data.frame(k = 1:3, log2a = log2(vapply(1:3, function(k) {
  as.numeric(unlist(ex$selected$alpha[[k]]))[1]
}, numeric(1))))
sel7g$gcv <- vapply(1:3, function(k) d7g$gcv[d7g$k == k][match(sel7g$log2a[k], d7g$log2a[d7g$k == k])], numeric(1))
sel7g$pc <- pc_lab(sel7g$k)
p7g <- ggplot(d7g, aes(log2a, gcv)) +
  geom_line(colour = "grey35", linewidth = 0.35) +
  geom_point(colour = "grey35", size = 0.6) +
  geom_point(data = sel7g, aes(shape = "selected"), colour = pal_method[["ReHPCA"]], size = 1.9) +
  facet_wrap(~pc, nrow = 1, scales = "free_y", labeller = label_parsed) +
  scale_shape_manual(NULL, values = c(selected = 16), labels = "selected value") +
  scale_y_log10() +
  labs(x = expression(log[2]~alpha), y = "GCV / min",
       title = "Smoothing: roughness penalty alpha by GCV (computed without the sparsity penalties)") +
  theme_paper()
p7 <- (p7u / p7f / p7v / p7g) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")


### Figure 8: selected sparsity across replicates -----------------------------------

d8 <- subset(M, method == "ReHPCA" & sigma == sigma_show)
block_lab <- factor(c("u", "fd", "nfd"), levels = c("u", "fd", "nfd"), labels = c(
  sprintf("atop('Zero scores', '(of'~n == %d*')')", N),
  sprintf("atop('Zero B-spline coef.', '(of'~d == %d*')')", n_basis),
  sprintf("atop('Zero vector loadings', '(of'~m == %d*')')", W)))
d8l <- do.call(rbind, lapply(c("u", "fd", "nfd"), function(b) {
  data.frame(block = b, scenario = d8$scenario, comp = d8$comp, value = d8[[paste0("sel_", b)]])
}))
truth8 <- do.call(rbind, lapply(names(scenarios), function(s) {
  rbind(data.frame(block = "u", scenario = s, comp = 1:3, truth = group_sizes(s)),
        data.frame(block = "fd", scenario = s, comp = 1:3, truth = fd_true_zeros),
        data.frame(block = "nfd", scenario = s, comp = 1:3, truth = W - b_size))
}))
for (nm in c("d8l", "truth8")) {
  d <- get(nm)
  d$block_lab <- block_lab[match(d$block, c("u", "fd", "nfd"))]
  d$comp_lab <- f_comp(d$comp, "text")
  d$scenario <- factor(d$scenario, levels = names(scenarios))
  assign(nm, d)
}
p8 <- ggplot(d8l, aes(scenario, value)) +
  geom_boxplot(width = 0.55, fill = "#CDE2FB", colour = "grey20", linewidth = 0.25,
               outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.12, height = 0, seed = 2), size = 0.5,
             colour = pal_method[["ReHPCA"]], alpha = 0.7, stroke = 0) +
  geom_errorbar(data = truth8, aes(scenario, ymin = truth, ymax = truth, linetype = "truth"),
                inherit.aes = FALSE, width = 0.75, colour = col_truth, linewidth = 0.45) +
  facet_grid(block_lab ~ comp_lab, scales = "free_y", labeller = label_parsed) +
  scale_linetype_manual(NULL, values = c(truth = "22"), labels = "true number of zeros") +
  labs(x = "Score-sparsity scenario", y = "Selected degree of sparsity (number of zeros)") +
  theme_paper()


### Figure 9: cumulative percentage of explained variance ----------------------------

cp <- subset(sim$cpev, sigma == sigma_show)
cp_long <- rbind(
  data.frame(cp[cp$method == "HPCA", c("scenario", "rep", "r")], measure = "HPCA",
             value = cp$adjusted[cp$method == "HPCA"]),
  data.frame(cp[cp$method == "ReHPCA", c("scenario", "rep", "r")], measure = "adjusted",
             value = cp$adjusted[cp$method == "ReHPCA"]),
  data.frame(cp[cp$method == "ReHPCA", c("scenario", "rep", "r")], measure = "naive",
             value = cp$naive[cp$method == "ReHPCA"])
)
agg9 <- summarise_by(cp_long, c("scenario", "measure", "r"), function(d) {
  q <- stats::quantile(d$value, c(0.05, 0.95), names = FALSE)
  data.frame(mean = mean(d$value), lo = q[1], hi = q[2])
})
agg9$measure <- factor(agg9$measure, levels = c("adjusted", "naive", "HPCA"))
agg9$scen_lab <- f_scen(agg9$scenario)
lab9 <- c(adjusted = "ReHPCA, adjusted CPEV (Section 5)",
          naive = "ReHPCA, sum of individual variances",
          HPCA = "HPCA (both measures agree)")
pd9 <- position_dodge(width = 0.3)
p9a <- ggplot(agg9, aes(r, 100 * mean, colour = measure, linetype = measure, shape = measure)) +
  geom_hline(yintercept = 100, colour = "grey55", linewidth = 0.3) +
  geom_line(linewidth = 0.45, position = pd9) +
  geom_errorbar(aes(ymin = 100 * lo, ymax = 100 * hi), width = 0, linewidth = 0.35,
                position = pd9, linetype = "solid") +
  geom_point(size = 1.3, position = pd9, fill = "white") +
  facet_wrap(~scen_lab, nrow = 1, labeller = label_parsed) +
  scale_colour_manual(NULL, values = c(adjusted = pal_method[["ReHPCA"]],
                                       naive = pal_method[["ReHPCA"]],
                                       HPCA = pal_method[["HPCA"]]), labels = lab9) +
  scale_linetype_manual(NULL, values = c(adjusted = "solid", naive = "22", HPCA = "solid"),
                        labels = lab9) +
  scale_shape_manual(NULL, values = c(adjusted = 16, naive = 21, HPCA = 17), labels = lab9) +
  scale_x_continuous(breaks = 1:3, expand = expansion(add = 0.25)) +
  labs(x = expression("Number of components"~italic(r)),
       y = "Cumulative % of variance",
       title = sprintf("Mean over %d replicates, with 5%%-95%% range", n_rep)) +
  theme_paper() +
  theme(legend.key.width = unit(18, "pt")) +
  guides(colour = guide_legend(nrow = 1))

# How much the naive sum overstates the adjusted CPEV for ReHPCA.
re9 <- cp[cp$method == "ReHPCA", ]
re9$excess <- 100 * (re9$naive - re9$adjusted)
re9$scen_lab <- f_scen(re9$scenario)
p9b <- ggplot(re9, aes(factor(r), excess)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_boxplot(width = 0.5, fill = "#CDE2FB", colour = "grey20", linewidth = 0.25,
               outlier.size = 0.45, outlier.stroke = 0.2) +
  facet_wrap(~scen_lab, nrow = 1, labeller = label_parsed) +
  labs(x = expression("Number of components"~italic(r)),
       y = "Naive minus adjusted\n(percentage points)",
       title = "ReHPCA: overstatement by the naive sum") +
  theme_paper()

p9 <- p9a / p9b + plot_layout(heights = c(1.25, 1)) + plot_annotation(tag_levels = "A")


## 7. Display ---------------------------------------------------------------------

figures <- list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5, p6 = p6, p7 = p7,
                p8 = p8, p9 = p9)

# Draw a figure at its design size. In a non-interactive session (Rscript)
# nothing is drawn, since R would otherwise write the plots to Rplots.pdf.
show_figure <- function(p, size) {
  if (!interactive()) return(invisible(NULL))
  if (preview_at_size) {
    grDevices::dev.new(width = size[1], height = size[2], noRStudioGD = TRUE)
    print(p)
  } else if (all(grDevices::dev.size("in") >= size)) {
    grid::grid.newpage()
    print(p, vp = grid::viewport(width = grid::unit(size[1], "in"),
                                 height = grid::unit(size[2], "in")))
  } else {
    print(p)  # device smaller than the design size: the figure is squeezed to fit
  }
  invisible(NULL)
}
for (nm in names(figures)) show_figure(figures[[nm]], fig_sizes[[nm]])
if (interactive() && !preview_at_size &&
    any(vapply(fig_sizes, function(s) any(grDevices::dev.size("in") < s), logical(1)))) {
  message("Some figures are larger than the plot window and were squeezed to fit; ",
          "enlarge the window or set preview_at_size <- TRUE to see them at print size.")
}


## 8. Numbers for the text and captions -------------------------------------------

show <- subset(M, sigma == sigma_show)
med <- function(x) if (all(is.na(x))) NA_real_ else stats::median(x, na.rm = TRUE)
summary_errors <- aggregate(cbind(err_psi, err_v) ~ method + scenario + comp, data = show,
                            FUN = med, na.action = stats::na.pass)
summary_errors[, c("err_psi", "err_v")] <- round(summary_errors[, c("err_psi", "err_v")], 3)
summary_sparsity <- aggregate(cbind(CIZ_pct = 100 * ciz / n_true_zero_v, MCR_pct = 100 * mcr) ~
                                scenario + sigma + comp,
                              data = subset(M, method == "ReHPCA"), FUN = mean)
summary_sparsity[, c("CIZ_pct", "MCR_pct")] <- round(summary_sparsity[, c("CIZ_pct", "MCR_pct")], 1)
gcv_at_min <- mean(M$log2_alpha[M$method == "ReHPCA"] == min(log2(alpha_grid)), na.rm = TRUE)
message(sprintf("Median errors at sigma = %s (Figures 4-5):", sigma_show))
print(summary_errors[order(summary_errors$method, summary_errors$comp, summary_errors$scenario), ],
      row.names = FALSE)
message("ReHPCA vector loadings, mean CIZ % and MCR % (cf. Table 1):")
print(summary_sparsity, row.names = FALSE)
message(sprintf("GCV chose the smallest alpha in the grid for %.0f%% of ReHPCA components.",
                100 * gcv_at_min))
message(sprintf("Smallest Err(psi_m) attainable with %d B-splines (m = 1, 2, 3): %s.",
                n_basis, paste(round(err_floor, 3), collapse = ", ")))


## 9. Saving (off by default) --------------------------------------------------------

# Nothing below runs unless save_figures is TRUE. Choose the figures and the
# formats first; each figure is written at its design size (fig_sizes).
if (save_figures) {
  save_which   <- names(figures)        # e.g. c("p1", "p4")
  save_formats <- c("pdf")              # e.g. c("pdf", "png", "tiff")
  out_dir <- file.path(pkg_root, "dev", "figures")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  for (nm in save_which) {
    for (fmt in save_formats) {
      file <- file.path(out_dir, paste0("figure_", nm, ".", fmt))
      device <- if (fmt == "pdf") grDevices::cairo_pdf else fmt
      ggplot2::ggsave(file, figures[[nm]], width = fig_sizes[[nm]][1],
                      height = fig_sizes[[nm]][2], units = "in", device = device,
                      dpi = 600, bg = "white")
    }
  }
  message("Figures written to ", out_dir)
}
