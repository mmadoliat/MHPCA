# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

MHPCA is an R package for multivariate hybrid PCA. It handles data that mixes functional variables (curves or surfaces, possibly on domains of different dimension) with ordinary vector variables. Components can be smoothed (roughness penalty, chosen by GCV) and made sparse (chosen by CV or eBIC). It is built on R6 and uses `fda` bases. It is not yet on CRAN; install from GitHub with `remotes::install_github("mmadoliat/MHPCA")`.

## Commands

R is not on PATH on this machine. Use `"C:\Program Files\R\R-4.3.3\bin\Rscript.exe"` (or `R.exe --vanilla --slave` with a heredoc from Bash; PowerShell 5.1 strips embedded quotes from `-e` strings).

```r
devtools::document()        # regenerate NAMESPACE + man/ from roxygen (never edit these by hand)
devtools::load_all()        # load the package in memory for interactive testing
devtools::check(args = "--as-cran")   # full R CMD check
tools::undoc(dir = "."); tools::checkDocFiles(dir = ".")   # fast doc-consistency checks on source
```

```r
devtools::test()                          # full suite (~45 s, mostly fitting the six yeast models)
devtools::test(filter = "yeast-models")   # one file: tests/testthat/test-yeast-models.R
testthat::snapshot_review()               # inspect snapshot diffs after an intentional change
```

The suite in `tests/testthat/` is a regression harness built from `dev/real_data.Rmd` (yeast data from `spls`):
- `helper-yeast.R` rebuilds the Rmd's data and runs its six `Mhpca()` fits once per session, in the Rmd's order from `set.seed(123)`. Order matters: each fit consumes RNG (CV fold shuffling) that the next one depends on, so don't reorder or insert fits.
- Expected values are mostly snapshots in `tests/testthat/_snaps/` (tolerance 1e-6; `helper-snapshots.R`). They pin current results, so any change that alters results changes specific snapshots. Review each diff (`testthat::snapshot_review()`) and never accept snapshot changes blindly.
- Tests use `skip_on_cran()`, so a plain `R CMD check` skips them unless `NOT_CRAN=true`; `devtools::test()` sets it.
- When patching functions in memory under `load_all()`, patch both `asNamespace("MHPCA")` and `as.environment("package:MHPCA")`. Tests see the attached copies.

## Architecture

Three layers:

1. **Data classes (R6).** Each class has a lowercase R6 generator plus a capitalized wrapper (`mfd`/`Mfd`, `hd`/`Hd`, `mhpca`/`Mhpca`, ...). `nfd()` and `img()` are exceptions: they are S3-classed matrices, not R6.
   - `basismfd` / `mvbasismfd` hold fda bases and the Gram matrix. `mfd` / `mvmfd` hold functional data as basis coefficients.
   - `nfd` / `mvnfd` hold vector data. `img` is a subclass of `nfd`.
   - `hd` is the hybrid container: `$mf` is an `mvmfd` and `$nf` is an `mvnfd`, and either may be NULL. Everything downstream branches on which parts are present.
2. **Methods.** Each class has a `*Methods.R` file with S3 operators (`+ - * [`) and helpers named `<op>_<class>` (`mean_mfd`, `inprod_hd`, `norm_mvnfd`, `scale_hd`, `center_hd`, ...).
   - `length`, `mean`, `sd`, `plot` and `scale` are exported functions in `R/mfdMethods.R` that mask base and forward to internal generics `op_length`, `op_mean`, ... (with `.default` falling back to base/stats). A new class adds `op_<name>.<class>` methods in its own `*Methods.R` file, tagged `@exportS3Method NULL`. Do not register methods on base's own generics (e.g. `length.nfd`); that would change what base R sees for these objects.
3. **Models and algorithms.**
   - `Mhpca()` (`R/hpcaClass.R`) is the main entry point. Its `initialize` first normalizes the tuning inputs (vector/matrix/list → per-variable lists), then dispatches:
     - `method = "power"`, `alpha_orth = FALSE` → `sequential_power_hybrid()` (one PC at a time with deflation; the sparse/smooth tuning path).
     - `alpha_orth = TRUE` → `joint_power_hybrid()`.
     - `method = "eigen"` → `eigen_approach_hybrid()`.
   - All of these live in `R/power_approach_hybrid.R` and `R/eigen_approach_hybrid.R`.
   - Sparse tuning goes `cv_gcv_sequential_hybrid` → `handle_sparse_tuning_hybrid` (CV, via `cv_local_hybrid`) or `eBIC_sequential_hybrid` (`R/eBIC.R`). Then `handle_smooth_tuning_hybrid` (GCV) runs, and `init_sequential_hybrid` does the actual power iteration.
   - `Remfpca()` (`R/remfpcaClass.R` + `R/power_approach.R` + `R/eigen_approach.R`) is an older, functional-only duplicate of this stack. `Mhpca` also accepts mfd-only input.

## Things that bite

- **Three penalty blocks.** Sparsity has three parallel blocks: `u` (scores), `nfd` (vector loadings) and `fd` (functional loadings). Each needs its own `sparse_tuning_*`, `sparse_type_*`, `nfold_*` and `pen_*`, and a block is only penalized when its `pen_* = TRUE`, it has a grid, and its data part exists (grids for a missing part are dropped in `mhpca$initialize`). A sparsity value is an integer count of entries set to zero (see `sparse_pen_fun`), not a lambda. Block-specific behaviour is described once in `sparse_blocks` (`R/sparse_blocks.R`); CV and eBIC tuning both run through `tune_sparse_blocks()`, each supplying its own `evaluate` function, and CV score tables come from `format_cv_scores()`.
- **Adding a model option or output.** Fit-wide options travel as one `settings` list built by `mhpca_settings()` (`R/mhpca_settings.R`). A new option goes into the `Mhpca()` and `mhpca$initialize()` signatures, `mhpca_settings()`, and is read as `settings$<name>` where it is used. `Mhpca()` forwards its arguments automatically (`forward_args()`). A new output goes into `mhpca_fields` (top of `R/hpcaClass.R`), which generates its private slot and read-only binding. Tuning inputs are reshaped by `normalize_tuning()` (`R/tuning_inputs.R`). Leaf functions (`init_sequential_hybrid`, `init_joint_hybrid`) and the exported eBIC functions still take explicit arguments.
- **Defaults differ.** The `Mhpca()` wrapper defaults are not the same as the `mhpca$new()` defaults (e.g. `smoothing_type`). `sequential_power_hybrid` has a single PC loop for both `smooth_GCV` settings; the places where they differ are explicit `smooth_grid` / `concat_scores` conditions.
- **Return types vary.** `pc_hd` and `mean_hd` are an `hd`, `mvmfd` or `mvnfd` depending on the input. All three algorithms return result lists with the same names, which `mhpca$initialize` reads with exact `[[ ]]` indexing.
- **Parallel mode** (`parallel = TRUE`) only reads `foreach::getDoParWorkers()`. It then builds its own PSOCK cluster with `library(MHPCA)` on the workers, so workers run the *installed* package, not `load_all()` code.
- **R files are collated alphabetically.** A function defined in two files silently resolves to the later one; define each function once.
- **2-D functional variables** store coefficients as an array (nbasis1 x nbasis2 x nobs). Average over observations with `basis_margins()` and flatten with `coef_matrix()`; plain `apply(coefs, 1, ...)` or `rbind()` silently does the wrong thing.
- **Conventions for new code.** Validate choice / flag arguments up front with `check_choice()` / `check_flag()` (`R/utils.R`). Progress text goes through `message()`, never `cat()`. Build block-diagonal matrices with `block_diag()`; apply something to each part of an `hd` with `map_hd()`; get the functional part of any data or component object with `functional_part()`. Component plots for both `mhpca` and `remfpca` go through `plot_components()`.
- **Roxygen markdown is off.** Escape `%` as `\%` in roxygen text; R code must be ASCII-only.
- **`dev/` is excluded from the build** (`.Rbuildignore`). `dev/real_data.Rmd` is an analysis script on the yeast cell-cycle data (`spls::yeast`) that depends on external files and many extra packages. It is not a vignette.

## Rules for this project

- Never run git commands that change the repository (no commit, push,
  merge, rebase, reset, or checkout of other branches). I handle git myself.
- The public API must not change: exported function names, argument names,
  argument order, defaults, and the structure of returned objects stay the
  same. If an API change seems worthwhile, propose it to me instead of
  making it.
- The testthat suite in tests/testthat/ (built from dev/real_data.Rmd) is
  the acceptance test. It must pass after every change. If a change is meant
  to alter results, say which snapshots change and why before updating them.
- Work in small steps. Run the tests after each step and show me what
  changed before moving on.
- Known bugs go in dev/known_issues.md, not in this file.
