# Literal ports of the Hotspot (Python) loops, used as references for the
# vectorized implementation.
ref_non_redundant <- function(nbr, w) {
  w2 <- w
  for (i in seq_len(nrow(nbr))) {
    for (k in seq_len(ncol(nbr))) {
      j <- nbr[i, k]
      if (j < i) next
      for (k2 in seq_len(ncol(nbr))) {
        if (nbr[j, k2] == i) {
          w_ji <- w2[j, k2]
          w2[j, k2] <- 0
          w2[i, k] <- w2[i, k] + w_ji
        }
      }
    }
  }
  w2
}

ref_local_cov <- function(x, nbr, w) {
  out <- 0
  for (i in seq_along(x)) for (k in seq_len(ncol(nbr))) {
    out <- out + x[i] * x[nbr[i, k]] * w[i, k]
  }
  out
}

ref_cond_eg2 <- function(x, nbr, w) {
  t1 <- numeric(length(x))
  for (i in seq_along(x)) for (k in seq_len(ncol(nbr))) {
    j <- nbr[i, k]
    t1[i] <- t1[i] + w[i, k] * x[j]
    t1[j] <- t1[j] + w[i, k] * x[i]
  }
  sum(t1^2)
}

ref_pair_lc <- function(x, y, nbr, w) {
  out <- 0
  for (i in seq_along(x)) for (k in seq_len(ncol(nbr))) {
    j <- nbr[i, k]
    out <- out + w[i, k] * (x[i] * y[j] + y[i] * x[j]) / 2
  }
  2 * out
}

test_that("kNN graph excludes self and matches non-redundant weights", {
  set.seed(1)
  X <- scale(matrix(rnorm(40 * 6), 40, 6))
  for (weighted in c(FALSE, TRUE)) {
    g <- transmap:::.hotspot_knn_graph(X, n_neighbors = 5,
                                       weighted_graph = weighted)
    expect_false(any(g$neighbors == seq_len(40)))
    w_nr <- ref_non_redundant(g$neighbors, g$weights)
    expect_equal(g$Wtot2, sum(w_nr^2))
    D_ref <- numeric(40)
    for (i in 1:40) for (k in 1:5) {
      D_ref[i] <- D_ref[i] + w_nr[i, k]
      D_ref[g$neighbors[i, k]] <- D_ref[g$neighbors[i, k]] + w_nr[i, k]
    }
    expect_equal(unname(g$D), D_ref)
  }
})

test_that("autocorrelation and pairwise statistics match the reference loops", {
  set.seed(2)
  X <- scale(matrix(rnorm(50 * 4), 50, 4))
  colnames(X) <- paste0("g", 1:4)
  g <- transmap:::.hotspot_knn_graph(X, n_neighbors = 6, weighted_graph = TRUE)
  w_nr <- ref_non_redundant(g$neighbors, g$weights)
  auto <- transmap:::.hotspot_autocorrelation(X, g)
  G_ref <- apply(X, 2, ref_local_cov, nbr = g$neighbors, w = w_nr)
  Z_ref <- G_ref / sqrt(sum(w_nr^2))
  expect_equal(auto$Z[match(colnames(X), auto$gene)], unname(Z_ref))

  lc <- transmap:::.hotspot_local_correlation(X, g)
  eg2 <- unname(apply(X, 2, ref_cond_eg2, nbr = g$neighbors, w = w_nr))
  for (a in 1:3) for (b in (a + 1):4) {
    lcv <- ref_pair_lc(X[, a], X[, b], g$neighbors, w_nr)
    zab <- lcv / sqrt(eg2[a])
    zba <- lcv / sqrt(eg2[b])
    expect_equal(lc$Z[a, b], if (abs(zab) < abs(zba)) zab else zba)
  }
  expect_equal(unname(diag(lc$Z)), rep(0, 4))
  expect_true(isSymmetric(lc$Z))
})

ref_smooth <- function(x, nbr, w, lambda = 0.9) {
  out <- numeric(length(x))
  den <- numeric(length(x))
  for (i in seq_along(x)) for (k in seq_len(ncol(nbr))) {
    j <- nbr[i, k]
    out[i] <- out[i] + x[j] * w[i, k]
    out[j] <- out[j] + x[i] * w[i, k]
    den[i] <- den[i] + w[i, k]
    den[j] <- den[j] + w[i, k]
  }
  lambda * out / den + (1 - lambda) * x
}

test_that(".hotspot_module_score is the signed PC1 of smoothed expression", {
  set.seed(6)
  X <- scale(matrix(rnorm(40 * 5), 40, 5) + rnorm(40))
  g <- transmap:::.hotspot_knn_graph(X, n_neighbors = 5, weighted_graph = TRUE)
  w_nr <- ref_non_redundant(g$neighbors, g$weights)
  S <- apply(X, 2, ref_smooth, nbr = g$neighbors, w = w_nr)
  pc <- prcomp(S, center = TRUE, scale. = FALSE)
  ref <- pc$x[, 1] * if (mean(pc$rotation[, 1]) < 0) -1 else 1
  score <- transmap:::.hotspot_module_score(X, g)
  expect_equal(unname(score), unname(ref))
})

test_that(".hotspot_hclust_to_linkage follows the scipy layout", {
  hc <- hclust(dist(c(0, 1, 5, 6.5)), method = "average")
  L <- transmap:::.hotspot_hclust_to_linkage(hc)
  expect_equal(L, cbind(child_1 = c(0, 2, 4), child_2 = c(1, 3, 5),
                        height = c(1, 1.5, 5.25), size = c(2, 2, 4)))
})

block_z <- function() {
  genes <- paste0("g", 1:13)
  Z <- matrix(0, 13, 13, dimnames = list(genes, genes))
  Z[1:6, 1:6] <- 10
  Z[7:12, 7:12] <- 10
  Z[13, 1:6] <- Z[1:6, 13] <- 1
  diag(Z) <- 0
  Z
}

test_that(".hotspot_modules splits two correlated blocks", {
  for (core in c(FALSE, TRUE)) {
    m <- transmap:::.hotspot_modules(block_z(), min_gene = 3, fdr = 0.05,
                                     core_only = core)
    a <- m$assignment
    expect_equal(names(a), paste0("g", 1:13))
    expect_length(unique(a[1:6]), 1)
    expect_length(unique(a[7:12]), 1)
    expect_setequal(c(a[1], a[7]), 1:2)
    expect_equal(unname(a[13]), -1L)
  }
  expect_equal(m$offset, 10)
  expect_equal(m$z_threshold, 10)
  expect_equal(nrow(m$linkage), 12)
})

test_that(".hotspot_modules leaves genes unassigned below min_gene", {
  m <- transmap:::.hotspot_modules(block_z(), min_gene = 7)
  expect_true(all(m$assignment == -1L))
})

test_that("r_hotspot uses half the samples as neighbours by default", {
  hs <- r_hotspot(transmap_example$Y, min_gene = 5)
  expect_equal(ncol(hs$graph$neighbors), floor(nrow(transmap_example$Y) / 2))
})

test_that("r_hotspot accepts a latent space and unscaled input", {
  Y <- transmap_example$Y
  X <- scale(Y)
  a <- r_hotspot(Y, n_neighbors = 20, min_gene = 5)
  b <- r_hotspot(Y, latent = X[rev(rownames(X)), ], n_neighbors = 20,
                 min_gene = 5)
  expect_equal(a$autocorrelation, b$autocorrelation)
  c <- r_hotspot(X, n_neighbors = 20, min_gene = 5, scale = FALSE)
  expect_equal(a$modules, c$modules)
  bad <- X[-1, ]
  expect_error(r_hotspot(Y, latent = bad), "same samples")
})

test_that("r_hotspot returns no modules without autocorrelated genes", {
  set.seed(8)
  Y <- matrix(rnorm(60 * 8), 60, 8,
              dimnames = list(paste0("s", 1:60), paste0("g", 1:8)))
  hs <- r_hotspot(Y, fdr = 1e-12)
  expect_length(hs$modules, 0)
  expect_null(hs$module_scores)
})

test_that("r_hotspot recovers the simulated modules", {
  hs <- r_hotspot(transmap_example$Y, min_gene = 5)
  expect_s3_class(hs, "r_hotspot")
  truth <- transmap_example$modules
  for (m in truth) {
    hit <- vapply(hs$modules, function(x) length(intersect(x, m)), numeric(1))
    expect_equal(max(hit), length(m))
  }
  expect_equal(dim(hs$module_scores),
               c(nrow(transmap_example$Y), length(hs$modules)))
  expect_true(all(unlist(hs$modules) %in% colnames(transmap_example$Y)))
  expect_equal(nrow(hs$linkage), nrow(hs$local_correlation_z) - 1)
})

test_that("r_hotspot output feeds transmap_local", {
  hs <- r_hotspot(transmap_example$Y, min_gene = 5)
  res <- transmap_local(transmap_example$Y, hs, method = "PC1")
  expect_equal(colnames(res$scores), names(hs$modules))
})

test_that("r_hotspot validates input", {
  Y <- transmap_example$Y
  Y[, 1] <- 1
  expect_error(r_hotspot(Y), "zero variance")
  expect_error(r_hotspot(transmap_example$Y, n_neighbors = 1000),
               "n_neighbors")
  expect_error(r_hotspot(transmap_example$Y, min_gene = 1), "min_gene")
})

test_that(".hotspot_modules gives a merged node the label of the looser block", {
  # two blocks (Z = 10 and Z = 8) linked by Z = 5; g13 joins both at Z = 4,
  # after the blocks merge, so it inherits the label chosen for that merge
  mk <- function(z1, z2) {
    g <- paste0("g", 1:13)
    Z <- matrix(4, 13, 13, dimnames = list(g, g))
    Z[1:6, 1:6] <- z1
    Z[7:12, 7:12] <- z2
    Z[1:6, 7:12] <- Z[7:12, 1:6] <- 5
    diag(Z) <- 0
    Z
  }
  a <- transmap:::.hotspot_modules(mk(10, 8), min_gene = 3)$assignment
  expect_equal(unname(a[13]), unname(a[7]))
  b <- transmap:::.hotspot_modules(mk(8, 10), min_gene = 3)$assignment
  expect_equal(unname(b[13]), unname(b[1]))
  expect_false(a[1] == a[7])
  # looser block formed first (first child): g1-g3 at Z = 8 complete before
  # g4-g9, whose mean merge height is lower (Z = 10 plus one gene at Z = 7)
  g <- paste0("g", 1:10)
  Z <- matrix(4, 10, 10, dimnames = list(g, g))
  Z[1:3, 1:3] <- 8
  Z[4:8, 4:8] <- 10
  Z[9, 4:8] <- Z[4:8, 9] <- 7
  Z[1:3, 4:9] <- Z[4:9, 1:3] <- 5
  diag(Z) <- 0
  m <- transmap:::.hotspot_modules(Z, min_gene = 3)
  expect_equal(unname(m$linkage[8, c("child_1", "child_2")]), c(15, 16))
  expect_equal(unname(m$assignment[10]), unname(m$assignment[1]))
  expect_false(m$assignment[1] == m$assignment[4])
})

test_that(".hotspot_modules sets the Z threshold above all Z without signal", {
  Z <- matrix(0.1, 5, 5, dimnames = list(letters[1:5], letters[1:5]))
  diag(Z) <- 0
  m <- transmap:::.hotspot_modules(Z, min_gene = 2)
  expect_equal(m$z_threshold, 1.1)
  expect_true(all(m$assignment == -1L))
})

test_that("r_hotspot reports progress when verbose", {
  msg <- testthat::capture_messages(
    r_hotspot(transmap_example$Y, min_gene = 5, verbose = TRUE))
  expect_match(msg[1], "60 of 60 genes are autocorrelated at FDR < 0.05")
  expect_match(msg[2], "Found [0-9]+ module")
  set.seed(8)
  Yn <- matrix(rnorm(60 * 8), 60, 8,
               dimnames = list(paste0("s", 1:60), paste0("g", 1:8)))
  msg <- testthat::capture_messages(
    r_hotspot(Yn, fdr = 1e-12, verbose = TRUE))
  expect_match(msg[2], "Fewer than two")
})
