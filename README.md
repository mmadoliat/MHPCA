# MHPCA: Multivariate Hybrid Principal Component Analysis

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: GPL (>= 2)](https://img.shields.io/badge/license-GPL%20(%3E%3D%202)-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

**MHPCA** is an R package for principal component analysis of *hybrid* data:
observations that combine functional variables (curves, or surfaces on
two-dimensional domains) with ordinary vector variables. All variables are
analyzed jointly, so each principal component has a functional part and a
vector part that describe the same mode of variation.

Components can be **smoothed**, with a roughness penalty for each functional
variable chosen by generalized cross-validation, and made **sparse** in the
scores, the vector loadings and the functional loadings, chosen by
cross-validation or the extended BIC.

## Features

- **Hybrid data.** Functional and vector variables in a single analysis;
  either part may be absent, so the package also covers multivariate
  functional PCA and ordinary (sparse) PCA.
- **Mixed domains.** Functional variables may live on domains of different
  dimension: 1-D curves alongside 2-D surfaces built from tensor-product
  bases. Grayscale images can enter as vector variables through `img()`.
- **Smooth components.** A separate roughness-penalty parameter for each
  functional variable, selected by GCV over a user-supplied grid.
- **Sparse components.** Independent sparsity penalties (soft, hard, SCAD or
  adaptive thresholding) on the scores, the vector loadings and the
  functional loadings, tuned by k-fold cross-validation (minimum or one-SE
  rule) or the extended BIC.
- **Three algorithms.** A sequential power method with deflation (the
  default, and the one that supports sparsity), a joint power method, and an
  eigendecomposition approach.
- **Object-oriented design.** R6 classes built on `fda` bases, with
  arithmetic, subsetting, inner products, norms, scaling and plotting methods.
- **Parallel tuning.** Cross-validation can run on parallel workers.

## Installation

MHPCA is not yet on CRAN. Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("mmadoliat/MHPCA")
```

MHPCA requires R (>= 4.0) and depends on `R6`, `fda`, `Matrix`, `expm` and
`foreach`.

## Quick start

The example below simulates 60 observations, each made up of one curve and a
vector of 12 features, both driven by the same two latent factors.

```r
library(MHPCA)
library(fda)

set.seed(1)
n <- 60                                   # observations
t_grid <- seq(0, 1, length.out = 101)     # grid for the functional variable
z <- matrix(rnorm(2 * n), nrow = 2)       # two latent factors

# Functional variable: one curve per observation (grid points x observations)
curves <- cbind(sin(2 * pi * t_grid), cos(2 * pi * t_grid)) %*% z +
  matrix(rnorm(101 * n, sd = 0.2), nrow = 101)

# Vector variable: 12 features per observation (observations x features)
features <- t(z) %*% matrix(rnorm(2 * 12), nrow = 2) +
  matrix(rnorm(n * 12, sd = 0.2), nrow = n)

# Build the hybrid data object
basis  <- Basismfd(create.bspline.basis(c(0, 1), nbasis = 20))
x_fd   <- Mfd(argval = t_grid, X = curves, mdbs = basis)
x_nfd  <- nfd(features)
hd_obj <- Hd(x_fd, x_nfd)
```

### Smooth components

```r
fit <- Mhpca(hd_obj, ncomp = 2,
             smooth_tuning = list(2^seq(-20, 0, length.out = 10)))

fit$values          # variance explained by each component
fit$smooth_tuning   # smoothing parameter selected by GCV

scores <- t(inprod_hd(fit$pc_hd, hd_obj))   # 60 x 2 matrix of PC scores
plot(fit)                                   # mean curve +/- each component
```

### Sparse components

A sparsity parameter is the **number of entries set to zero**, not a penalty
weight. Give each grid in increasing order.

```r
# Sparse vector loadings, tuned by 5-fold cross-validation
fit_cv <- Mhpca(hd_obj, ncomp = 2,
                sparse_tuning_nfd = 0:10, pen_nfd = TRUE, nfold_nfd = 5)
fit_cv$sparse_tuning_nfd    # features set to zero, per component
plot_cv(fit_cv, type = "nfd")

# The same model, tuned by the extended BIC
fit_ebic <- Mhpca(hd_obj, ncomp = 2,
                  sparse_tuning_nfd = 0:10, pen_nfd = TRUE,
                  sparse_selection_method = "ebic")
```

Sparsity on the scores and on the functional loadings works the same way,
through `sparse_tuning_u` / `pen_u` and `sparse_tuning_fd` / `pen_fd`. The
smoothing and sparsity grids can be combined in a single call.

## Main functions

### Building data objects

| Function | Description |
|---|---|
| `Basismfd()` | Basis for a functional variable, from one or more `fda` bases. Several bases give a tensor-product basis on a multidimensional domain. |
| `Mfd()` | A functional variable, stored as basis coefficients. Takes raw values on a grid (`method = "data"`) or coefficients (`method = "coefs"`). |
| `nfd()` | A vector variable: a matrix with one row per observation. |
| `img()` | Grayscale images, vectorized into an `nfd` object. |
| `Mvmfd()`, `Mvnfd()` | Collections of several functional or vector variables. |
| `Hd()` | Hybrid data: any combination of `mfd` and `nfd` objects. |

Data objects support `+`, `-`, `*`, subsetting with `[`, `mean()`, `sd()`,
`scale()`, `plot()`, inner products (`inprod_hd()`, `inprod_mfd()`, ...) and
norms (`norm_hd()`, `norm_mfd()`, ...). `bimfdplot()` plots two functional
variables against each other.

### Fitting a model

`Mhpca()` is the main entry point. Its most important arguments are:

| Argument | Description |
|---|---|
| `ncomp` | Number of components to estimate. |
| `method`, `alpha_orth` | `"power"` with `alpha_orth = FALSE` (default) estimates one component at a time with deflation; `alpha_orth = TRUE` estimates all components jointly; `"eigen"` uses an eigendecomposition. |
| `smooth_tuning`, `smooth_GCV` | Smoothing parameters for each functional variable. With `smooth_GCV = TRUE` (default) each entry is a grid searched by GCV; otherwise the values are used as given. |
| `sparse_tuning_u`, `sparse_tuning_nfd`, `sparse_tuning_fd` | Sparsity grids for the scores, the vector loadings and the functional loadings. |
| `pen_u`, `pen_nfd`, `pen_fd` | Switch each sparsity penalty on (all `FALSE` by default). |
| `sparse_type_u`, `sparse_type_nfd`, `sparse_type_fd` | Thresholding rule: `"soft"` (default), `"hard"`, `"SCAD"` or `"adaptive"`. |
| `sparse_selection_method` | `"cv"` (default) or `"ebic"`. |
| `nfold_u`, `nfold_nfd`, `nfold_fd`, `cv.pick` | Number of CV folds for each block, and the selection rule (`"1se"` (default) or `"min"`). |
| `parallel` | Run cross-validation on the number of workers registered with `foreach` (for example with `doParallel::registerDoParallel()`). |

See `?mhpca` for the full list.

### Working with the fit

`Mhpca()` returns an `mhpca` object with read-only fields:

| Field | Description |
|---|---|
| `pc_hd` | The principal components, with the same structure as the input (`hd`, `mvmfd` or `mvnfd`). |
| `values` | Variance explained by each component. |
| `lsv` | Left singular vectors, one column per component. |
| `mean_hd` | The mean of the data. |
| `smooth_tuning`, `sparse_tuning_u`, `sparse_tuning_nfd`, `sparse_tuning_fd` | The selected tuning parameters. |
| `GCVs`, `CVs_u`, `CVs_nfd`, `CVs_fd`, `eBICs_u`, `eBICs_nfd`, `eBICs_fd` | The tuning criteria over each grid. |

`plot()` shows each functional component as perturbations of the mean, and
`plot_cv()` shows the cross-validation curves used to select sparsity.

## Datasets

| Dataset | Description |
|---|---|
| `motion_sense_data` | Smartphone accelerometer and gyroscope signals (user acceleration and pitch attitude) from 24 participants performing four activities. |
| `electrical_power_data` | Daily curves of household electric power and voltage at 5-minute resolution, starting 17 December 2006. |

## Development

The test suite in `tests/testthat/` covers the data classes, their methods
and every model-fitting path. It includes a regression harness built on the
yeast cell-cycle data from the `spls` package, which fits each algorithm and
compares the results with stored snapshots. To run it:

```r
devtools::test()
```

Bug reports and feature requests are welcome on the
[issue tracker](https://github.com/mmadoliat/MHPCA/issues).

## Citation

The methodology is described in Mahmoudiandehkordi et al. (2025). A full
reference will be added here once the paper is available.

## Authors

- Mehdi Maadooliat ([ORCID](https://orcid.org/0000-0002-5408-2676))
- Hossein Haghbin, maintainer ([ORCID](https://orcid.org/0000-0001-8416-2354))
- Yue Zhao ([ORCID](https://orcid.org/0009-0000-4561-9163))
- Soroush Mahmoudiandehkordi

## License

GPL (>= 2)
