Y <- transmap_example$Y
null_z <- transmap_example$null_z
E1 <- Y[, transmap_example$modules$module_1]
E12 <- Y[, c(transmap_example$modules$module_1, transmap_example$modules$module_2)]
pinv_root <- transmap:::.transmap_pseudo_inverse_root

# .local_single_gene ------------------------------------------------------

test_that(".local_single_gene returns the gene itself", {
  fit <- transmap:::.local_single_gene(Y[, "GENE50", drop = FALSE], "m", "_WPC1")
  expect_equal(colnames(fit$scores), "m_WPC1")
  expect_equal(unname(fit$scores[, 1]), unname(Y[, "GENE50"]))
  expect_equal(fit$loadings,
               matrix(1, 1, 1, dimnames = list("GENE50", "m_WPC1")))
  expect_equal(fit$k, 1L)
  expect_true(is.na(fit$eigenvalues))
  expect_true(fit$info$single_gene)
})

# .local_pc1 ----------------------------------------------------------------

test_that(".local_pc1 equals the first principal component", {
  fit <- transmap:::.local_pc1(E1, "m")
  pc <- prcomp(E1, center = TRUE, scale. = TRUE)
  expect_equal(colnames(fit$scores), "m")
  expect_equal(rownames(fit$scores), rownames(E1))
  expect_equal(unname(fit$scores[, 1]), unname(pc$x[, 1]))
  expect_equal(unname(fit$loadings[, 1]), unname(pc$rotation[, 1]))
  expect_equal(unname(fit$eigenvalues), pc$sdev[1]^2)
  expect_equal(fit$info$var_explained, pc$sdev[1]^2 / sum(pc$sdev^2))
})

test_that(".local_pc1 handles a single gene", {
  fit <- transmap:::.local_pc1(E1[, 1, drop = FALSE], "m")
  expect_equal(colnames(fit$scores), "m")
  expect_equal(unname(fit$scores[, 1]), unname(E1[, 1]))
})

# .local_apc ----------------------------------------------------------------

test_that(".local_apc keeps one component per latent factor", {
  fit <- transmap:::.local_apc(E12, "m")
  expect_equal(fit$k, 2L)
  expect_equal(colnames(fit$scores), c("m_PC1", "m_PC2"))
  expect_equal(fit$info$mp_upper, (1 + sqrt(24 / 300))^2 / 24)
  pc <- prcomp(E12, center = TRUE, scale. = TRUE)
  expect_equal(unname(fit$scores), unname(pc$x[, 1:2]))
})

test_that(".local_apc stops at max_cumvar", {
  fit <- transmap:::.local_apc(E12, "m", max_cumvar = 0.01)
  expect_equal(fit$k, 1L)
  fit1 <- transmap:::.local_apc(E1[, 1, drop = FALSE], "m")
  expect_equal(colnames(fit1$scores), "m_PC1")
})

# .local_wpc ----------------------------------------------------------------

ref_wpc <- function(E, Sigma) {
  pc <- prcomp(scale(E) %*% pinv_root(Sigma), center = FALSE, scale. = FALSE)
  k <- sum(pc$sdev^2 > 1)
  list(scores = pc$x[, seq_len(k), drop = FALSE], ev = pc$sdev[seq_len(k)]^2)
}

test_that(".local_wpc whitens with the null correlation (|z| < 2)", {
  fit <- transmap:::.local_wpc(E1, "m", null_z)
  z <- null_z[, colnames(E1)]
  idx <- which(apply(abs(z), 1, max) < 2)
  Sigma <- as.matrix(Matrix::nearPD(cor(z[idx, ]), corr = TRUE)$mat)
  ref <- ref_wpc(E1, Sigma)
  expect_equal(fit$info$null_threshold, 2)
  expect_equal(fit$info$null_variants, length(idx))
  expect_equal(unname(fit$scores), unname(ref$scores))
  expect_equal(unname(fit$eigenvalues), ref$ev)
  expect_true(all(fit$eigenvalues > 1))
  expect_equal(colnames(fit$scores), paste0("m_WPC", seq_len(fit$k)))
  expect_equal(rownames(fit$loadings), colnames(E1))
})

test_that(".local_wpc relaxes the null threshold to |z| < 3", {
  set.seed(4)
  p <- ncol(E1)
  small <- matrix(runif(5 * p, -1.9, 1.9), 5, p)
  mid <- matrix(runif(30 * p, 2.1, 2.9) * sample(c(-1, 1), 30 * p, TRUE), 30, p)
  big <- matrix(5, 10, p)
  z <- rbind(small, mid, big)
  colnames(z) <- colnames(E1)
  fit <- transmap:::.local_wpc(E1, "m", z)
  expect_equal(fit$info$null_threshold, 3)
  expect_equal(fit$info$null_variants, 35)
})

test_that(".local_wpc falls back to the identity without null variants", {
  z <- matrix(5, 20, ncol(E1), dimnames = list(NULL, colnames(E1)))
  fit <- transmap:::.local_wpc(E1, "m", z)
  expect_true(is.na(fit$info$null_threshold))
  expect_equal(fit$info$null_variants, 0)
  ref <- ref_wpc(E1, diag(ncol(E1)))
  expect_equal(unname(fit$scores), unname(ref$scores))
})

test_that(".local_wpc returns no component when none has eigenvalue > 1", {
  # two anti-correlated genes; the null correlation is ~1, so the whitening
  # keeps only the direction (1, 1), which carries no expression variance
  set.seed(5)
  x <- rnorm(100)
  E <- cbind(a = x, b = -x)
  r <- runif(200, -1.5, 1.5)
  z <- cbind(a = r, b = r + rnorm(200, sd = 0.03))
  fit <- transmap:::.local_wpc(E, "m", z)
  expect_equal(fit$k, 0L)
  expect_null(fit$scores)
})

test_that(".local_wpc errors on a non-finite null correlation", {
  z <- null_z[, colnames(E1)]
  z[, 1] <- 0
  expect_error(transmap:::.local_wpc(E1, "m", z), "not finite")
})

test_that(".local_wpc handles a single gene", {
  fit <- transmap:::.local_wpc(E1[, 1, drop = FALSE], "m", null_z)
  expect_equal(colnames(fit$scores), "m_WPC1")
})

# .local_aucell -------------------------------------------------------------

test_that(".local_aucell matches a direct AUCell run", {
  skip_if_not_installed("AUCell")
  mods <- transmap_example$modules
  set.seed(9)
  got <- transmap:::.local_aucell(Y, mods, auc_max_rank = 0.2)
  set.seed(9)
  rk <- AUCell::AUCell_buildRankings(t(Y), plotStats = FALSE, verbose = FALSE)
  auc <- AUCell::AUCell_calcAUC(mods, rk, aucMaxRank = nrow(rk) * 0.2,
                                verbose = FALSE)
  ref <- t(AUCell::getAUC(auc))
  expect_equal(rownames(got), rownames(Y))
  expect_equal(got, ref[rownames(Y), ], ignore_attr = TRUE)
})
