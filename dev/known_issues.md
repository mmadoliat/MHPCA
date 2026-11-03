# Known issues

## Open: need a decision from the authors

These are not clear-cut bugs; changing them needs a methodological or
packaging decision.

| # | Issue | Where |
|---|---|---|
| A | The DOI in `DESCRIPTION` is a placeholder (`10.48550/arXiv.xxxx.yyyyy`). CRAN's incoming check reports it as invalid. | `DESCRIPTION` |
| B | The TF counts in the comments of `dev/real_data.Rmd` (6 / 8 / 12 non-zero TFs) cannot come from the Rmd's own call: its grid `seq(0, 105, length.out = 3)` only allows 106, 54 or 1 non-zero TFs. With that grid CV keeps all 106 (in every commit since 2025-05-20 that runs the call); with a finer grid (`seq(0, 105, by = 5)`) it keeps 71-81 (1se rule). The comments must come from older code or a different setup. | `dev/real_data.Rmd` |
| C | With smoothing, components are normalized in the penalized norm (Silverman-type regularized FPCA), so `values` are not ordered and can be tiny when GCV picks a large smoothing parameter (seen on synthetic data where GCV chose the top of the grid). Check whether GCV should choose the grid edge in such cases. | `R/power_approach_hybrid.R` (`gcv_local_hybrid`) |
| F | The 1se CV rule picks the largest candidate index within one SE of the minimum, so grids must be ordered by increasing sparsity. | `R/sparse_blocks.R` (`pick_cv_index`) |
| G | Dataset documentation (`motion_sense_data`, `electrical_power_data`) describes the structure only; the source / citation should be added. | `R/MotionSenseDataset.R` |
| H | `method = "eigen"` with `smooth_tuning = NULL` is not the unsmoothed fit that the `smooth_tuning` documentation promises ("If NULL, non-smoothing MFPCA is estimated"): `eigen_approach_hybrid()` then searches the default grid `2^seq(-20, 20, length.out = 5)` by GCV. On a [0, 1] domain with 15 B-splines it picked 2^-20, which smooths noticeably; the power path with `NULL` does not smooth. | `R/eigen_approach_hybrid.R` |
| I | Sparse functional loadings have no exact zeros: the thresholded coefficients go through `G_half` and back (`G_half_inverse %*% G_half`), so they come back as ~1e-18. `coef == 0` finds none (the vector loadings and scores are exact zeros). Counting zeros needs a tolerance (as in `dev/paper_figures.R`). | `R/power_approach_hybrid.R` (`init_sequential_hybrid`) |
| J | The smoothing GCV ignores the selected sparsity: `cv_gcv_sequential_hybrid()` calls `handle_smooth_tuning_hybrid()` without `pen_*`, so each GCV fit has no sparsity penalties (the code comment says this is deliberate). The `GCVs` field documentation ("computed at the selected sparsity") and Section 4.2 of the paper (alpha chosen "given the sparsity pattern") say otherwise. | `R/power_approach_hybrid.R` (`cv_gcv_sequential_hybrid`), `R/hpcaClass.R` (`GCVs` docs) |

## Fixed

All fixes have regression tests (`tests/testthat/test-data-classes.R`,
`test-mhpca-paths.R`, `test-generics.R`, `test-yeast-*.R`).

| # | Issue (now fixed) |
|---|---|
| 1 | `length()` of an `hd` object returned `NULL`. |
| 2 | Spurious "did not match ... adjusted" warning for a grid given for a single variable. |
| 3 | The `sparse_tuning_nfd` warning named `hd_obj$mf$nvar`. |
| 4 | A wrong-length `sparse_tuning_nfd` list was recycled to `hd_obj$mf$nvar` (crash when the counts differ). |
| 5 | `CVs_u` was gated on `pen_fd` instead of `pen_u` (missing table, or crash). |
| 6 | Joint path without `smooth_tuning` crashed on single-part data. |
| 7 | Assigning to `fit$values` reported "`$sigma` is read only". |
| 8 | Inconsistent selection reported for untuned sparsity blocks. |
| 9 | Selections had no `PC`/`Var` names with `smooth_GCV = FALSE`. |
| 10 | `smooth_GCV = FALSE` crashed on vector-only data (`C_temp` not found). |
| 11 | `fit$mean_hd` was the mean of the principal components, not of the data. |
| 12 | The eigen path returned `smooth_tuning` and `GCVs` as `NULL`; results were read by `$` partial matching. |
| 13 | eBIC tuned sparsity on untransformed coefficients instead of `C %*% G_half`. |
| 14 | For PCs after the first, tuning used the non-deflated non-functional data. |
| 15 | Sparsity grids and `pen_*` flags for absent data parts (or without a grid) were tuned on empty data, crashing about 30% of random configurations. |
| 16 | eBIC indexed an empty smoothing-matrix list; a one-variable matrix `smooth_tuning` collapsed to a vector. |
| 17 | `2 * hd`, `-hd`, `hd - hd`, integer scalars (`2L * mfd`), `mfd[]` errored. |
| 18 | `mvnfd + mvnfd` added the wrong variables; `inprod_mvnfd` failed with 3+ variables; `inprod_hd` checked the wrong object. |
| 19 | `hd[, 2]` returned the functional part; j bounds checks tested `i`. |
| 20 | `scale()` of an `mfd` silently converted it to a B-spline basis; `scale()` of a vector-only `hd` errored. |
| 21 | 2-D functional data: `mean()`/`sd()` averaged over the wrong dimensions, `center_hd()` failed, and the sequential and joint paths failed. All three PCA paths now agree on 2-D data. |
| 22 | Duplicate `bimfdplot()` definitions (the `hd` version was unreachable). |
| 23 | Broken help examples, undocumented objects and datasets, Rd argument mismatches, non-ASCII code, `:::` self-calls, unescaped `%` in Rd. |
| 24 | `scale(hd, method = "balanced")` compared a sum of norms with a sum of squared norms; both parts now get the same total variation. |
| 25 | eBIC sparsity tuning used the smoothing matrix of the first point of the GCV grid, so results depended on the grid order. With a GCV grid it now tunes without smoothing (as CV does); with fixed smoothing it uses that smoothing. |
