Y <- transmap_example$Y
modules <- transmap_example$modules
null_z <- transmap_example$null_z

test_that("PC1 equals the first principal component of each module", {
  res <- transmap_local(Y, modules, method = "PC1")
  expect_s3_class(res, "transmap_local")
  expect_equal(colnames(res$scores), names(modules))
  expect_equal(rownames(res$scores), rownames(Y))
  pc <- prcomp(Y[, modules$module_2], center = TRUE, scale. = TRUE)
  expect_equal(unname(res$scores[, "module_2"]), unname(pc$x[, 1]))
})

test_that("aPC keeps more components for a two-factor module", {
  two <- list(both = c(modules$module_1, modules$module_2))
  res <- transmap_local(Y, two, method = "aPC")
  expect_equal(res$summary$k, 2L)
  expect_equal(colnames(res$scores), c("both_PC1", "both_PC2"))
})

test_that("WPC whitens with the null correlation and keeps eigenvalues > 1", {
  res <- transmap_local(Y, modules, method = "WPC", null_z = null_z)
  expect_true(all(res$eigenvalues > 1))
  expect_true(all(res$summary$null_threshold == 2))
  expect_equal(sum(res$summary$k), ncol(res$scores))
  # list(z = ...) input, as stored in null-z RDS files, gives the same result
  res2 <- transmap_local(Y, modules, method = "WPC", null_z = list(z = null_z))
  expect_identical(res$scores, res2$scores)
})

test_that("WPC requires null z-scores covering the module genes", {
  expect_error(transmap_local(Y, modules), "null_z")
  expect_error(transmap_local(Y, modules, null_z = null_z[, -1]),
               "missing 1 module gene")
})

test_that("single-gene and empty modules are handled", {
  mods <- list(one = "GENE50", none = c("NOT_A_GENE"), m1 = modules$module_1)
  for (m in c("WPC", "PC1", "aPC")) {
    res <- transmap_local(Y, mods, method = m, null_z = null_z)
    expect_equal(res$summary$k, c(1L, 0L, res$summary$k[3]))
    first <- colnames(res$scores)[1]
    expect_equal(unname(res$scores[, first]), unname(Y[, "GENE50"]))
  }
  expect_error(transmap_local(Y, "NOT_A_GENE", method = "PC1"), "No module genes")
})

test_that("WPC skips modules without a component of eigenvalue > 1", {
  set.seed(5)
  x <- rnorm(nrow(Y))
  Y2 <- cbind(Y, ANTI1 = x, ANTI2 = -x)
  r <- runif(nrow(null_z), -1.5, 1.5)
  z2 <- cbind(null_z, ANTI1 = r, ANTI2 = r + rnorm(nrow(null_z), sd = 0.03))
  mods <- list(anti = c("ANTI1", "ANTI2"), m1 = modules$module_1)
  res <- transmap_local(Y2, mods, null_z = z2)
  expect_equal(res$summary$k[1], 0L)
  expect_false(any(grepl("^anti", colnames(res$scores))))
  expect_null(res$loadings$anti)
  expect_error(transmap_local(Y2, mods["anti"], null_z = z2),
               "No program phenotype")
})

test_that("module input forms are accepted", {
  a <- transmap_local(Y, modules$module_1, method = "PC1")
  expect_equal(colnames(a$scores), "module")
  b <- transmap_local(Y, unname(modules), method = "PC1")
  expect_equal(colnames(b$scores), paste0("module_", 1:3))
})

test_that("covariates are regressed out before building programs", {
  cov <- transmap_example$covariates
  adj <- qr.resid(qr(cbind(1, cov)), Y)
  a <- transmap_local(adj, modules, method = "WPC", null_z = null_z)
  b <- transmap_local(Y, modules, method = "WPC", null_z = null_z,
                      covariates = cov)
  expect_equal(a$scores, b$scores)
  # samples missing from covariates are dropped
  c <- transmap_local(Y, modules, method = "PC1", covariates = cov[1:250, ])
  expect_equal(nrow(c$scores), 250)
})

test_that("AUCell scores are returned per module", {
  skip_if_not_installed("AUCell")
  set.seed(1)
  res <- transmap_local(Y, modules, method = "AUCell")
  expect_equal(dim(res$scores), c(nrow(Y), length(modules)))
  expect_true(all(res$scores >= 0 & res$scores <= 1))
})

test_that("transmap_local reports progress when verbose", {
  mods <- list(none = "NOT_A_GENE", m1 = modules$module_1)
  msg <- testthat::capture_messages(
    transmap_local(Y, mods, method = "PC1", verbose = TRUE))
  expect_match(msg, "Skipping 1 module", all = FALSE)
  expect_match(msg, "Constructed 1 program phenotype", all = FALSE)
  set.seed(5)
  x <- rnorm(nrow(Y))
  Y2 <- cbind(Y, ANTI1 = x, ANTI2 = -x)
  r <- runif(nrow(null_z), -1.5, 1.5)
  z2 <- cbind(null_z, ANTI1 = r, ANTI2 = r + rnorm(nrow(null_z), sd = 0.03))
  mods2 <- list(anti = c("ANTI1", "ANTI2"), m1 = modules$module_1)
  msg <- testthat::capture_messages(
    transmap_local(Y2, mods2, null_z = z2, verbose = TRUE))
  expect_match(msg, "anti: no component kept", all = FALSE)
})

test_that("missing null z genes are listed, truncated after ten", {
  expect_error(transmap_local(Y, modules["module_1"], null_z = null_z[, -(1:11)]),
               "missing 11 module gene\\(s\\): GENE01, .*GENE10, \\.\\.\\.")
})

test_that("AUCell scores are residualized on covariates", {
  skip_if_not_installed("AUCell")
  cov <- transmap_example$covariates
  adj <- qr.resid(qr(cbind(1, cov)), Y)
  a <- transmap_local(adj, modules, method = "AUCell")
  b <- transmap_local(Y, modules, method = "AUCell", covariates = cov)
  expect_equal(b$scores, qr.resid(qr(cbind(1, cov)), a$scores),
               ignore_attr = TRUE)
  expect_equal(b$summary$k, rep(1L, 3))
})
