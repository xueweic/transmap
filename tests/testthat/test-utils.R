check_matrix <- transmap:::.transmap_check_matrix
residualize <- transmap:::.transmap_residualize
drop_constant <- transmap:::.transmap_drop_constant
pinv_root <- transmap:::.transmap_pseudo_inverse_root
as_modules <- transmap:::.transmap_as_module_list

# .transmap_check_matrix --------------------------------------------------

test_that(".transmap_check_matrix coerces data frames", {
  df <- data.frame(a = 1:3, b = c(2, 4, 6), row.names = c("s1", "s2", "s3"))
  m <- check_matrix(df, "x")
  expect_true(is.matrix(m))
  expect_equal(dimnames(m), list(c("s1", "s2", "s3"), c("a", "b")))
})

test_that(".transmap_check_matrix rejects invalid input", {
  expect_error(check_matrix(matrix(letters[1:4], 2), "x"), "numeric matrix")
  expect_error(check_matrix(1:4, "x"), "numeric matrix")
  x <- matrix(1:4, 2)
  expect_error(check_matrix(x, "x"), "row names")
  rownames(x) <- c("a", "b")
  expect_error(check_matrix(x, "x"), "column names")
  expect_silent(check_matrix(x, "x", need_colnames = FALSE))
  colnames(x) <- c("g", "g")
  expect_error(check_matrix(x, "x"), "duplicated column names")
  y <- matrix(1:4, 2, dimnames = list(c("a", "a"), c("g1", "g2")))
  expect_error(check_matrix(y, "x"), "duplicated row names")
  y <- matrix(c(1, NA, 3, 4), 2, dimnames = list(c("a", "b"), c("g1", "g2")))
  expect_error(check_matrix(y, "x"), "missing values")
})

# .transmap_residualize ----------------------------------------------------

set.seed(1)
Y_res <- matrix(rnorm(50 * 3), 50, 3,
                dimnames = list(paste0("s", 1:50), paste0("g", 1:3)))
C_res <- matrix(rnorm(50 * 2), 50, 2, dimnames = list(paste0("s", 1:50), NULL))

test_that(".transmap_residualize matches lm residuals", {
  expect_identical(residualize(Y_res, NULL), Y_res)
  r <- residualize(Y_res, C_res)
  expect_equal(r, unname(residuals(lm(Y_res ~ C_res))), ignore_attr = TRUE)
  expect_equal(dimnames(r), dimnames(Y_res))
  expect_equal(crossprod(cbind(1, C_res), r), matrix(0, 3, 3),
               tolerance = 1e-10, ignore_attr = TRUE)
})

test_that(".transmap_residualize aligns and drops samples", {
  r <- residualize(Y_res, C_res[50:11, ])
  expect_equal(rownames(r), rownames(Y_res)[11:50])
  expect_equal(r, residualize(Y_res[11:50, ], C_res[11:50, ]))
  expect_message(residualize(Y_res, C_res[1:40, ], verbose = TRUE),
                 "keeps 40 of 50")
  expect_error(residualize(Y_res, C_res[1:3, ]), "Too few samples")
})

# .transmap_drop_constant --------------------------------------------------

test_that(".transmap_drop_constant removes zero-variance columns", {
  Y <- cbind(a = c(1, 2, 3, 4, 5), b = rep(2, 5), c = c(1, 3, 2, 5, 4))
  expect_equal(colnames(drop_constant(Y)), c("a", "c"))
  expect_message(drop_constant(Y, verbose = TRUE), "Dropping 1 genes")
  expect_silent(drop_constant(Y[, c("a", "c")], verbose = TRUE))
})

# .transmap_pseudo_inverse_root -------------------------------------------

test_that(".transmap_pseudo_inverse_root whitens a full-rank matrix", {
  set.seed(2)
  A <- matrix(rnorm(25), 5)
  S <- crossprod(A) + diag(5)
  P <- pinv_root(S, tol = 1)
  expect_true(isSymmetric(P, tol = 1e-10))
  expect_equal(P %*% S %*% P, diag(5), tolerance = 1e-8)
})

test_that(".transmap_pseudo_inverse_root truncates small eigenvalues", {
  # the first eigenvalue explains > 99.9%, so only it is kept (L = 1)
  P <- pinv_root(diag(c(1000, 1e-3)))
  expect_equal(P, diag(c(1 / sqrt(1000), 0)))
  expect_error(pinv_root(matrix(1:6, 2)), "square")
})

# .transmap_as_module_list -------------------------------------------------

test_that(".transmap_as_module_list normalises module input", {
  expect_equal(as_modules(c("a", "b")), list(module = c("a", "b")))
  expect_equal(names(as_modules(list("a", "b"))), c("module_1", "module_2"))
  expect_equal(names(as_modules(list(x = "a", "b"))), c("x", "module_2"))
  expect_identical(as_modules(list(m = c("a", "a", NA, "b")))$m, c("a", "b"))
  hs <- structure(list(modules = list(m1 = c("a", "b"))), class = "r_hotspot")
  expect_equal(as_modules(hs), list(m1 = c("a", "b")))
})

test_that(".transmap_as_module_list rejects invalid input", {
  expect_error(as_modules(list()), "`module` must be")
  expect_error(as_modules(1:3), "`module` must be")
  expect_error(as_modules(list(1:3)), "character vector")
  expect_error(as_modules(list(a = "x", a = "y")), "duplicated names")
})

test_that(".transmap_pseudo_inverse_root keeps every eigenvalue when tol is unreachable", {
  set.seed(21)
  A <- crossprod(matrix(rnorm(30), 10, 3))
  P <- transmap:::.transmap_pseudo_inverse_root(A, tol = 2)
  e <- eigen(A, symmetric = TRUE)
  expect_equal(P, e$vectors %*% diag(1 / sqrt(e$values)) %*% t(e$vectors))
})

test_that(".transmap_pseudo_inverse_root strict = TRUE needs the ratio above tol", {
  # eigenvalues 3, 1: the first explains exactly 0.75
  A <- diag(c(3, 1))
  expect_equal(transmap:::.transmap_pseudo_inverse_root(A, tol = 0.75),
               diag(c(1 / sqrt(3), 0)))
  expect_equal(transmap:::.transmap_pseudo_inverse_root(A, tol = 0.75,
                                                        strict = TRUE),
               diag(c(1 / sqrt(3), 1)))
})
