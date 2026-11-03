# Analysis chunks of dev/real_data.Rmd that use the fitted models: scores,
# TF selection, evaluated loadings / mean, clustering of scores, sparsity vs
# cluster enrichment, and the TF-time map.
#
# Not covered (they need files outside this repo or only draw plots):
# the happ comparison (source("happ_hpca_file.R")), the SupSFPCA / SupPCA
# .mat results, and the ggplot / factoextra / GGally figures.
testthat::local_edition(3)

val_tf <- c("ACE2", "FKH1", "FKH2", "GAT3", "MBP1", "MCM1", "NDD1", "REB1",
            "STB1", "STE12", "SWI4", "SWI5", "SWI6")
non_val_tf <- c("ARG80", "CIN5", "DOT6", "HAA1", "HIR1", "HIR2", "HSF1", "JNO4",
                "MAL13", "MET4", "MSN2", "PHD1", "PHO2", "PUT3", "RIM101",
                "SFL1", "SIP4", "SOK2", "YAP5")

kmeans_rows <- function(scores, k) {
  set.seed(123)
  km <- stats::kmeans(scores, centers = k)
  split(seq_len(nrow(scores)), km$cluster)
}

test_that("scores t(inprod_hd(pc_hd, hd_obj)) for each hybrid fit", {
  skip_if_no_yeast()
  for (name in c("power_smooth", "power_sparse", "joint", "eigen")) {
    scores <- yeast_scores(name)
    expect_equal(dim(scores), c(542L, 4L), info = name)
    expect_snapshot_large(scores)
  }
})

test_that("non-zero TF loadings per PC (sparse fit)", {
  skip_if_no_yeast()
  nfpc_sparse <- t(yeast_fit("power_sparse")$pc_hd$nf$data[[1]])
  expect_equal(dim(nfpc_sparse), c(106L, 4L))

  non_sparse_tf <- lapply(seq_len(ncol(nfpc_sparse)), function(i) {
    rownames(nfpc_sparse)[nfpc_sparse[, i] != 0]
  })
  # Current code: no TF loading is set to zero (see test-yeast-models.R).
  expect_equal(lengths(non_sparse_tf), c(106L, 106L, 106L, 106L))

  all_tf <- unique(unlist(non_sparse_tf))
  expect_equal(sum(all_tf %in% val_tf), 13L)
  expect_equal(sum(all_tf %in% non_val_tf), 18L)
  expect_equal(sum(!all_tf %in% c(val_tf, non_val_tf)), 75L)
  expect_snapshot_small(nfpc_sparse)
})

test_that("loadings and mean evaluated on the time grid", {
  skip_if_no_yeast()
  d <- yeast_data()
  fit <- yeast_fit("power_sparse")
  tt <- seq(d$t0_f, d$t1_f, length.out = 119)

  pc_loadings <- fit$pc_hd$eval(tt)
  expect_named(pc_loadings, c("fd1", "nfd1"))
  expect_equal(dim(pc_loadings[[1]]), c(119L, 4L))
  expect_snapshot_large(pc_loadings)

  # mean_hd is the mean of the data
  ave <- fit$mean_hd$eval(tt)
  expect_equal(dim(ave[[1]]), c(119L, 1L))
  expect_equal(as.vector(ave[[1]]), rowMeans(d$hd_obj$mf$eval(tt)[[1]]))
  expect_equal(as.vector(ave[[2]]), colMeans(d$hd_obj$nf$data[[1]]), ignore_attr = TRUE)
  expect_snapshot_large(ave)
})

test_that("normalised functional loadings (method comparison chunk)", {
  skip_if_no_yeast()
  d <- yeast_data()
  sparse_mf <- yeast_fit("power_sparse")$pc_hd$mf[, 1]
  tt <- seq(d$t0_f, d$t1_f, length.out = 120)

  norms <- norm_mfd(sparse_mf)
  expect_length(norms, 4)
  expect_snapshot_small(norms)

  scaled <- Mfd(X = sparse_mf$coefs %*% diag(1 / norms), mdbs = sparse_mf$basis,
                method = "coefs")
  expect_equal(norm_mfd(scaled), rep(1, 4))
  expect_snapshot_large(scaled$eval(tt))

  remfpca_eval <- yeast_fit("remfpca")$pc_hd$eval(tt)
  expect_equal(dim(remfpca_eval[[1]]), c(120L, 4L))
  expect_snapshot_large(remfpca_eval)
})

test_that("k-means clusters of the sparse-fit scores", {
  skip_if_no_yeast()
  scores <- yeast_scores("power_sparse")
  for (k in c(2, 5, 6, 10)) {
    rows <- kmeans_rows(scores, k)
    expect_length(rows, k)
    expect_equal(sum(lengths(rows)), 542L)
    expect_snapshot_small(lengths(rows))
  }
  # The Rmd's downstream chunks use the k = 5 clustering.
  expect_equal(unname(lengths(kmeans_rows(scores, 5))), c(29L, 77L, 128L, 211L, 97L))
})

test_that("sparsity in u vs cluster membership (overlap, Fisher tests)", {
  skip_if_no_yeast()
  fit <- yeast_fit("power_sparse")
  cluster_rows <- kmeans_rows(yeast_scores("power_sparse"), 5)

  sparse_us <- lapply(seq_len(ncol(fit$lsv)), function(i) unname(which(fit$lsv[, i] == 0)))
  overlap_matrix <- sapply(cluster_rows, function(b) {
    sapply(sparse_us, function(a) length(intersect(a, b)))
  })
  sparse_counts <- lengths(sparse_us)
  cluster_counts <- lengths(cluster_rows)
  expect_equal(sparse_counts, c(40L, 80L, 80L, 80L))
  expect_snapshot_small(overlap_matrix)

  expected_matrix <- outer(sparse_counts, cluster_counts,
                           FUN = function(r, c) (r * c) / sum(overlap_matrix))
  expect_snapshot_small(overlap_matrix / expected_matrix)

  total_items <- sum(cluster_counts)
  p_val_matrix <- matrix(NA_real_, nrow = length(sparse_counts), ncol = length(cluster_counts))
  for (i in seq_along(sparse_counts)) {
    for (j in seq_along(cluster_counts)) {
      A <- overlap_matrix[i, j]
      B <- sparse_counts[i] - A
      C <- cluster_counts[j] - A
      D <- total_items - (A + B + C)
      p_val_matrix[i, j] <- stats::fisher.test(
        matrix(c(A, B, C, D), nrow = 2, byrow = TRUE), alternative = "greater"
      )$p.value
    }
  }
  expect_snapshot_small(p_val_matrix)
})

test_that("TF-time map from reconstructed data", {
  skip_if_no_yeast()
  d <- yeast_data()
  fit <- yeast_fit("power_sparse")
  scores <- yeast_scores("power_sparse")

  y_reconstruction <- scores %*% t(fit$pc_hd$mf$eval(seq(d$t0_f, d$t1_f, length = 120))[[1]])
  x_reconstruction <- scores %*% fit$pc_hd$nf$data[[1]]
  xy_map <- t(x_reconstruction) %*% y_reconstruction
  expect_equal(dim(xy_map), c(106L, 120L))
  expect_snapshot_large(xy_map)
})
