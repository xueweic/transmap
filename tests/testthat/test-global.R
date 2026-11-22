z <- transmap_example$z
n <- transmap_example$n
Y <- transmap_example$Y
dl <- list(Z = z, N = n)

test_that("r_transcca returns named sparse programs", {
  prog <- r_transcca(dl, K = 3, sumabs = 0.4)
  expect_s3_class(prog, "r_transcca")
  expect_equal(dim(prog$V), c(ncol(z), 3))
  expect_equal(dim(prog$U), c(nrow(z), 3))
  expect_equal(rownames(prog$V), colnames(z))
  expect_equal(rownames(prog$U), rownames(z))
  expect_equal(colnames(prog$V), paste0("Program_", 1:3))
  expect_true(any(prog$V == 0))
  expect_equal(prog$sparsity, 0.4)
  # the strongest program is driven by the simulated variant block
  top <- order(abs(prog$U[, 1]), decreasing = TRUE)[1:10]
  expect_true(all(top <= 20))
})

test_that("r_transcca agrees with a direct PMA::PMD fit", {
  prog <- r_transcca(dl, K = 2, sumabs = 0.5)
  W <- z / sqrt(n - 2 + z^2)
  utils::capture.output(fit <- suppressWarnings(
    PMA::PMD(W, K = 2, sumabs = 0.5, type = "standard", center = FALSE)))
  expect_equal(unname(prog$V), fit$v)
  expect_equal(unname(prog$d), fit$d)
  expect_equal(unname(prog$d2_score[1]), fit$d[1]^2)
})

test_that("different_sparsity tunes each program", {
  set.seed(3)
  prog <- r_transcca(dl, K = 2, different_sparsity = TRUE)
  expect_length(prog$sparsity, 2)
  expect_equal(dim(prog$V), c(ncol(z), 2))
})

test_that("r_transcca validates input", {
  expect_error(r_transcca(list(Z = z, N = 1)), "`N`")
  expect_error(r_transcca(dl, K = 0), "`K`")
  expect_null(r_transcca(dl, K = 1, null_thres = 10)$V)
})

test_that("transmap_global scores are Y %*% V on shared genes", {
  prog <- r_transcca(dl, K = 2, sumabs = 0.4)
  res <- transmap_global(Y, prog)
  expect_s3_class(res, "transmap_global")
  expect_equal(res$scores, Y %*% prog$V, ignore_attr = TRUE)
  expect_equal(colnames(res$scores), colnames(prog$V))

  V <- prog$V[c(1:10, 1:2), ]
  rownames(V)[11:12] <- c("absent_1", "absent_2")
  res2 <- transmap_global(Y, V)
  expect_equal(res2$genes_missing, c("absent_1", "absent_2"))
  expect_equal(res2$genes_used, rownames(prog$V)[1:10])
})

test_that("transmap_global accepts a plain loading matrix and covariates", {
  V <- matrix(c(1, 1, 0, 0, 0, 1), ncol = 2,
              dimnames = list(c("GENE01", "GENE02", "GENE03"), NULL))
  res <- transmap_global(Y, V)
  expect_equal(colnames(res$scores), c("Program_1", "Program_2"))
  expect_equal(unname(res$scores[, 1]), unname(Y[, "GENE01"] + Y[, "GENE02"]))
  cov <- transmap_example$covariates
  res_c <- transmap_global(Y, V, covariates = cov)
  adj <- qr.resid(qr(cbind(1, cov)), Y)
  expect_equal(unname(res_c$scores[, 2]), unname(adj[, "GENE03"]))
})

test_that("d2_score follows the chosen normalization", {
  W <- z / sqrt(n - 2 + z^2)
  base <- r_transcca(dl, K = 2, sumabs = 0.5)
  fn <- r_transcca(dl, K = 2, sumabs = 0.5, normalization = "F_norm")
  mx <- r_transcca(dl, K = 2, sumabs = 0.5, normalization = "max")
  expect_equal(fn$V, base$V)
  expect_equal(fn$d2_score, base$d2_score / sum(W^2)^2)
  expect_equal(mx$d2_score, base$d2_score / max(abs(W))^2)
})

test_that("r_transcca chooses sumabs by cross-validation when not given", {
  W <- z / sqrt(n - 2 + z^2)
  set.seed(31)
  prog <- r_transcca(dl, K = 2)
  set.seed(31)
  expect_equal(prog$sparsity, transmap:::.transcca_cv_sumabs(W))
})

test_that("different_sparsity selects K from the singular values", {
  W <- z / sqrt(n - 2 + z^2)
  set.seed(32)
  prog <- r_transcca(dl, different_sparsity = TRUE)
  expect_equal(ncol(prog$V), transmap:::.transcca_select_K(svd(W)$d, "max"))
})

test_that("r_transcca reports the fitting mode when verbose", {
  utils::capture.output(msg <- testthat::capture_messages(suppressWarnings(
    r_transcca(dl, K = 1, sumabs = 0.5, verbose = TRUE))))
  expect_match(msg, "one sparsity", all = FALSE)
  set.seed(33)
  msg <- testthat::capture_messages(
    r_transcca(dl, K = 1, different_sparsity = TRUE, verbose = TRUE))
  expect_match(msg, "separate sparsity", all = FALSE)
})

test_that("transmap_global scale = TRUE scores standardized genes", {
  prog <- r_transcca(dl, K = 2, sumabs = 0.4)
  res <- transmap_global(Y, prog, scale = TRUE)
  expect_equal(res$scores, scale(Y) %*% prog$V, ignore_attr = TRUE)
  Yc <- cbind(Y, CONST = 1)
  V <- rbind(prog$V, CONST = 1)
  expect_equal(transmap_global(Yc, V, scale = TRUE)$genes_missing, "CONST")
})

test_that("transmap_global validates programs and reports missing genes", {
  empty <- r_transcca(dl, K = 1, null_thres = 10)
  expect_error(transmap_global(Y, empty), "contains no programs")
  V <- matrix(1, 2, 1, dimnames = list(c("absent_1", "absent_2"), "P"))
  expect_error(transmap_global(Y, V), "No program genes")
  V2 <- rbind(V, GENE01 = 1)
  expect_message(transmap_global(Y, V2, verbose = TRUE),
                 "2 program gene\\(s\\) are not in `Y`")
})
