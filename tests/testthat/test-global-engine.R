# Unit tests for the internal transCCA engine in R/transmap_global.R.

set.seed(10)
W_tall <- matrix(rnorm(12 * 5), 12, 5)
W_wide <- t(W_tall)

test_that(".transcca_singular_values: svd and eigen styles agree", {
  for (W in list(W_tall, W_wide)) {
    d <- transmap:::.transcca_singular_values(W, "svd")
    expect_equal(d, svd(W)$d)
    expect_equal(transmap:::.transcca_singular_values(W, "eigen")[1:5], d^2)
  }
})

test_that(".transcca_select_K picks the largest gap or a variance cut-off", {
  expect_equal(transmap:::.transcca_select_K(c(10, 9, 3, 2.5), "max"), 2)
  expect_equal(transmap:::.transcca_select_K(c(5, 3, 1, 1), "var_explain",
                                             var_exp = 2 / 3), 2)
  expect_equal(transmap:::.transcca_select_K(c(5, 3, 1, 1), "var_explain",
                                             var_exp = 0.85), 3)
})

test_that(".transcca_min_abs is the smallest admissible sparsity", {
  expect_equal(transmap:::.transcca_min_abs(matrix(0, 4, 9)), 0.51)
  expect_equal(transmap:::.transcca_min_abs(matrix(0, 9, 4)), 0.51)
})

test_that(".transcca_cv_sumabs returns a sparsity in the search grid", {
  set.seed(11)
  s <- transmap:::.transcca_cv_sumabs(W_tall)
  expect_gte(s, transmap:::.transcca_min_abs(W_tall))
  expect_lte(s, 0.8)
})

test_that(".transcca_soft soft-thresholds", {
  expect_equal(transmap:::.transcca_soft(c(-3, -0.5, 0, 2), 1), c(-2, 0, 0, 1))
})

test_that(".transcca_l2n is the L2 norm, with 0.05 for a zero vector", {
  expect_equal(transmap:::.transcca_l2n(c(3, 4)), 5)
  expect_equal(transmap:::.transcca_l2n(c(0, 0)), 0.05)
})

test_that(".transcca_binary_search hits the requested L1/L2 ratio", {
  expect_equal(transmap:::.transcca_binary_search(c(1, 0, 0), 2), 0)
  argu <- as.numeric(1:10)
  lam <- transmap:::.transcca_binary_search(argu, 2)
  su <- transmap:::.transcca_soft(argu, lam)
  expect_gt(lam, 0)
  expect_equal(sum(abs(su)) / sqrt(sum(su^2)), 2, tolerance = 1e-3)
})

test_that(".transcca_safe_svd equals svd", {
  expect_equal(transmap:::.transcca_safe_svd(W_tall), svd(W_tall))
})

test_that(".transcca_init_v spans the leading right singular vectors", {
  for (W in list(W_tall, W_wide)) {
    v <- transmap:::.transcca_init_v(W, 2)
    expect_equal(dim(v), c(ncol(W), 2))
    expect_equal(sqrt(colSums(v^2)), c(1, 1))
    expect_equal(abs(diag(crossprod(v, svd(W)$v[, 1:2]))), c(1, 1),
                 tolerance = 1e-6)
  }
})

test_that(".transcca_same_sparsity matches PMA::PMD and handles the null", {
  res <- transmap:::.transcca_same_sparsity(W_tall, K = 2, sumabs = 0.6)
  utils::capture.output(ref <- suppressWarnings(
    PMA::PMD(W_tall, K = 2, sumabs = 0.6, type = "standard", center = FALSE)))
  expect_equal(res$v, ref$v)
  expect_equal(res$u, ref$u)
  expect_equal(res$d, ref$d)
  expect_equal(res$sparsity, 0.6)
  expect_equal(res$s_values, svd(W_tall)$d)
  null <- transmap:::.transcca_same_sparsity(W_tall, K = 2, sumabs = 0.6,
                                             null_thres = 100)
  expect_null(null$v)
  # K = NULL selects K from the singular values
  auto <- transmap:::.transcca_same_sparsity(W_tall, sumabs = 0.6)
  expect_equal(ncol(auto$v),
               transmap:::.transcca_select_K(svd(W_tall)$d, "max"))
})

test_that(".transcca_diff_sparsity deflates one program at a time", {
  set.seed(12)
  res <- transmap:::.transcca_diff_sparsity(W_tall, K = 2)
  expect_equal(dim(res$u), c(12, 2))
  expect_equal(dim(res$v), c(5, 2))
  expect_length(res$sparsity, 2)
  W_use <- W_tall
  for (k in 1:2) {
    u <- res$u[, k]
    v <- res$v[, k]
    expect_true(sum(v^2) == 0 || abs(sum(v^2) - 1) < 1e-8)
    expect_equal(res$d[k], as.numeric(t(u) %*% W_use %*% v))
    W_use <- W_use - res$d[k] * u %*% t(v)
  }
})

test_that(".transcca_same_sparsity prints PMD progress when verbose", {
  quiet <- transmap:::.transcca_same_sparsity(W_tall, K = 2, sumabs = 0.6)
  utils::capture.output(loud <- suppressWarnings(
    transmap:::.transcca_same_sparsity(W_tall, K = 2, sumabs = 0.6,
                                       verbose = TRUE)))
  expect_equal(loud, quiet)
})

test_that(".transcca_diff_sparsity reports every tenth program", {
  set.seed(34)
  W <- matrix(rnorm(15 * 12), 15, 12)
  expect_message(res <- transmap:::.transcca_diff_sparsity(W, K = 10,
                                                           verbose = TRUE),
                 "Working on program K = 10")
  expect_equal(dim(res$v), c(12, 10))
})

test_that(".transcca_binary_search warns when the interval does not shrink", {
  expect_warning(lam <- transmap:::.transcca_binary_search(c(1e45, 1e45, 1), 1.2),
                 "converge")
  expect_true(is.finite(lam))
})

test_that(".transcca_safe_svd falls back to a random matrix when svd fails", {
  bad <- W_tall
  bad[1, 1] <- Inf
  set.seed(35)
  s <- transmap:::.transcca_safe_svd(bad)
  expect_equal(dim(s$u), dim(W_tall))
  expect_true(all(is.finite(s$d)))
})

test_that(".transcca_init_v handles non-finite wide matrices as PMA does", {
  # NA: both SVDs fail, the random fallback gives finite unit vectors
  bad <- W_wide
  bad[1, 1] <- NA
  set.seed(36)
  v <- transmap:::.transcca_init_v(bad, 2)
  expect_equal(dim(v), c(ncol(bad), 2))
  expect_false(anyNA(v))
  expect_equal(sqrt(colSums(v^2)), c(1, 1))
  # Inf survives into v and normalizing gives NaN, which is an error
  bad[1, 1] <- Inf
  expect_error(transmap:::.transcca_init_v(bad, 2), "some are NA")
})
