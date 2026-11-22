# Unit tests for w_construction() and its helpers: building W from a
# transBoost-style data_list (the W_calculation step).

set.seed(40)
p <- 25
g <- 8
N <- 150
nm <- list(paste0("v", 1:p), paste0("g", 1:g))
Gi <- scale(matrix(rbinom(N * p, 2, 0.3), N, p,
                   dimnames = list(paste0("s", 1:N), nm[[1]])))
Ei <- matrix(rnorm(N * g), N, g, dimnames = list(paste0("s", 1:N), nm[[2]]))
Ei[, 1:3] <- Ei[, 1:3] + Gi[, 1]
Ei <- scale(Ei)
Z <- crossprod(Gi, Ei) / sqrt(N)
dimnames(Z) <- nm
maf <- runif(p, 0.05, 0.5)
sdy <- runif(g, 0.5, 2)
beta <- Z / sqrt(N) * 0.7
se <- matrix(1 / sqrt(N), p, g, dimnames = nm)
pv <- pchisq((beta / se)^2, 1, lower.tail = FALSE)
Wz <- Z / sqrt(N - 2 + Z^2)
W_of <- function(...) unclass(w_construction(list(...)))
strip <- function(W) {
  attr(W, "input") <- NULL
  W
}

test_that(".transcca_input_type picks the input type and checks companions", {
  it <- transmap:::.transcca_input_type
  expect_equal(it(c("Z", "N")), "Z_N")
  expect_equal(it(c("Z", "N", "MAF", "sdY")), "Z_N")
  expect_equal(it("Beta"), "Beta")
  expect_equal(it(c("Beta", "N", "Pvalue")), "Beta_N")
  expect_equal(it(c("Beta", "N", "Z", "sdY")), "Beta_N")
  expect_equal(it(c("G", "E")), "G_E")
  expect_equal(it(c("log2FC", "N1", "N0")), "log2FC")
  expect_equal(it("Sigma_GE"), "Sigma_GE")
  expect_error(it(character()), "must hold association statistics")
  expect_error(it("Z"), "also needs `N`")
  expect_error(it("G"), "also needs `E`")
  expect_error(it(c("Sigma_GE", "Z")), "`Z` cannot be combined")
  expect_error(it(c("Beta", "N")), "one of `Pvalue`, `Z` or `seBeta`")
  expect_error(it(c("Beta", "N", "seBeta", "MAF")), "`MAF` cannot be")
  expect_error(it(c("Z", "N", "sdY")), "`sdY` is used with `MAF`")
})

test_that(".transcca_check_like, _check_vector and _check_n validate shapes", {
  expect_equal(transmap:::.transcca_check_like(5, Z, "N", allow_scalar = TRUE),
               5)
  expect_equal(transmap:::.transcca_check_like(as.data.frame(se), Z, "x"),
               as.matrix(as.data.frame(se)))
  expect_error(transmap:::.transcca_check_like(se[-1, ], Z, "seBeta"),
               "`seBeta` must be a numeric matrix")
  expect_error(transmap:::.transcca_check_like(c(1, 2), Z, "N", TRUE),
               "single number or")
  expect_equal(transmap:::.transcca_check_vector(1:3, 3, "x"), c(1, 2, 3))
  expect_error(transmap:::.transcca_check_vector(1:2, 3, "MAF"),
               "length 1 or 3")
  expect_error(transmap:::.transcca_check_vector(c(NA, 1, 1), 3, "MAF"),
               "`MAF`")
  expect_error(transmap:::.transcca_check_n(2, Z), "greater than 2")
  nmat <- matrix(N, p, g)
  expect_equal(transmap:::.transcca_check_n(nmat, Z), nmat)
  expect_error(transmap:::.transcca_check_n(matrix(NA_real_, p, g), Z),
               "`N` contains missing values")
})

test_that(".transcca_rescale scales rows by genotype sd and columns by sdY", {
  ref <- diag(sqrt(2 * maf * (1 - maf))) %*% Wz %*% diag(1 / sdy)
  expect_equal(transmap:::.transcca_rescale(Wz, maf, sdy), ref,
               ignore_attr = TRUE)
  expect_equal(transmap:::.transcca_rescale(Wz), Wz)
  expect_equal(transmap:::.transcca_rescale(Wz, maf = 0.5), Wz * sqrt(0.5))
})

test_that(".transcca_beta_scale averages genotype scales over signal genes", {
  sig <- which(colSums(pv < 1e-4) > 0)
  expect_gt(length(sig), 0)
  v <- sqrt(1 / (beta^2 + (N - 2) * se^2))
  expect_equal(transmap:::.transcca_beta_scale(beta, N, se_beta = se),
               rowMeans(v[, sig, drop = FALSE]))
  expect_equal(transmap:::.transcca_beta_scale(beta, N, pvalue = pv),
               transmap:::.transcca_beta_scale(beta, N, se_beta = se))
  expect_equal(transmap:::.transcca_beta_scale(beta, N, z = beta / se),
               transmap:::.transcca_beta_scale(beta, N, se_beta = se))
  vy <- sqrt(rep(sdy^2, each = p) / (beta^2 + (N - 2) * se^2))
  expect_equal(transmap:::.transcca_beta_scale(beta, N, se_beta = se,
                                               sd_y = sdy),
               rowMeans(vy[, sig, drop = FALSE]))
  # without signal genes, the 20% of genes with the smallest p-values are used
  top <- order(apply(pv, 2, min))[seq_len(ceiling(0.2 * g))]
  expect_equal(transmap:::.transcca_beta_scale(beta, N, se_beta = se,
                                               sig_cut = 0),
               rowMeans(v[, top, drop = FALSE]))
  # 0 / 0 scales (beta = 0 with p = 1) are set to 0
  b0 <- beta
  b0[1, ] <- 0
  expect_equal(unname(transmap:::.transcca_beta_scale(
    b0, N, pvalue = pv * 0 + 1, sig_cut = 0)[1]), 0)
})

test_that(".transcca_sigma returns given or in-sample covariance", {
  S <- cor(Gi)
  expect_identical(transmap:::.transcca_sigma(S, NULL, "Sigma_GG", p), S)
  expect_equal(transmap:::.transcca_sigma("insample", Gi, "Sigma_GG", p),
               cov(Gi))
  expect_equal(transmap:::.transcca_sigma(as.data.frame(S), NULL, "x", p), S)
  expect_error(transmap:::.transcca_sigma("insample", NULL, "Sigma_GG", p),
               "only with `G` and `E`")
  expect_error(transmap:::.transcca_sigma(S[-1, -1], NULL, "Sigma_GG", p),
               "25 x 25")
})

test_that(".transcca_W validates the data_list", {
  expect_error(transmap:::.transcca_W(Z), "named list")
  expect_error(transmap:::.transcca_W(list(Z, N)), "named list")
  expect_error(transmap:::.transcca_W(list(z = Z, N = N)),
               "Unknown element\\(s\\) in `data_list`: `z`")
  # NULL elements are ignored, as in transBoost
  expect_equal(transmap:::.transcca_W(list(Z = Z, N = N, MAF = NULL))$W, Wz)
})

test_that("w_construction follows the W_calculation formula of each input", {
  W <- w_construction(list(Z = Z, N = N))
  expect_equal(attr(W, "input"), "Z_N")
  expect_equal(strip(W), Wz)
  expect_equal(W_of(Sigma_GE = Wz), Wz, ignore_attr = TRUE)
  expect_equal(W_of(Z = Z, N = matrix(N, p, g)), Wz, ignore_attr = TRUE)
  expect_equal(W_of(Beta = beta), beta, ignore_attr = TRUE)
  expect_equal(W_of(Z = Z, N = N, MAF = maf, sdY = sdy),
               diag(sqrt(2 * maf * (1 - maf))) %*% Wz %*% diag(1 / sdy),
               ignore_attr = TRUE)
  expect_equal(W_of(G = Gi, E = Ei), crossprod(Gi, Ei) / N,
               ignore_attr = TRUE)
  n1 <- sample(20:40, p, TRUE)
  n0 <- sample(20:40, p, TRUE)
  expect_equal(W_of(log2FC = Z, N1 = n1, N0 = n0),
               diag(sqrt(n1 * n0) / (n1 + n0)) %*% Z, ignore_attr = TRUE)
  expect_equal(W_of(log2FC = Z, N1 = 30, N0 = 30), Z * sqrt(900) / 60,
               ignore_attr = TRUE)
  d_g <- transmap:::.transcca_beta_scale(beta, N, se_beta = se, sd_y = sdy)
  expect_equal(W_of(Beta = beta, N = N, seBeta = se, sdY = sdy),
               sweep(beta * d_g, 2, sdy, "/"), ignore_attr = TRUE)
  expect_equal(dimnames(W_of(Beta = beta, N = N, Pvalue = pv)), nm)
})

test_that("missing values in companion statistics are rejected", {
  se_na <- se
  se_na[1, 1] <- NA
  expect_error(W_of(Beta = beta, N = N, seBeta = se_na),
               "`seBeta` contains missing values")
  expect_error(W_of(Beta = beta, N = N, Pvalue = se_na),
               "`Pvalue` contains missing values")
})

test_that("exclude masks W before whitening", {
  ex <- matrix(FALSE, p, g)
  ex[1:5, 1:2] <- TRUE
  masked <- Wz
  masked[ex] <- 0
  expect_equal(strip(w_construction(list(Z = Z, N = N), exclude = ex)),
               masked)
  expect_equal(strip(w_construction(list(Z = Z, N = N), exclude = ex * 1)),
               masked)
  expect_equal(strip(w_construction(list(Z = Z, N = N),
                                    exclude = as.data.frame(ex))), masked)
  S <- cor(Ei)
  P <- transmap:::.transmap_pseudo_inverse_root(S, strict = TRUE)
  expect_equal(strip(w_construction(list(Z = Z, N = N, Sigma_EE = S),
                                    exclude = ex)),
               masked %*% P, ignore_attr = TRUE)
  expect_error(w_construction(list(Z = Z, N = N), exclude = ex[-1, ]),
               "`exclude` must be")
  expect_error(w_construction(list(Z = Z, N = N),
                              exclude = matrix("a", p, g)),
               "`exclude` must be")
})

test_that("Sigma_GG and Sigma_EE whiten rows and columns", {
  SG <- cor(Gi)
  SE <- cor(Ei)
  PG <- transmap:::.transmap_pseudo_inverse_root(SG, strict = TRUE)
  PE <- transmap:::.transmap_pseudo_inverse_root(SE, strict = TRUE)
  expect_equal(W_of(Z = Z, N = N, Sigma_GG = SG), PG %*% Wz,
               ignore_attr = TRUE)
  out <- W_of(Z = Z, N = N, Sigma_GG = SG, Sigma_EE = SE)
  expect_equal(out, PG %*% Wz %*% PE, ignore_attr = TRUE)
  expect_equal(dimnames(out), nm)
  # in-sample whitening of individual-level data equals whitening G and E
  Gt <- Gi %*% transmap:::.transmap_pseudo_inverse_root(cov(Gi), strict = TRUE)
  Et <- Ei %*% transmap:::.transmap_pseudo_inverse_root(cov(Ei), strict = TRUE)
  expect_equal(W_of(G = Gi, E = Ei, Sigma_GG = "insample",
                    Sigma_EE = "insample"),
               crossprod(Gt, Et) / N, ignore_attr = TRUE)
  expect_error(W_of(Z = Z, N = N, Sigma_EE = "insample"), "only with `G`")
})

test_that("individual-level input checks samples", {
  expect_error(W_of(G = Gi, E = Ei[-1, ]), "same samples \\(rows\\)")
  E2 <- Ei
  rownames(E2) <- rev(rownames(Ei))
  expect_error(W_of(G = Gi, E = E2), "same order")
})

test_that("r_transcca takes a data_list of any input type", {
  a <- r_transcca(list(Z = Z, N = N), K = 2, sumabs = 0.5)
  b <- r_transcca(list(Sigma_GE = Wz), K = 2, sumabs = 0.5)
  expect_equal(a$V, b$V)
  expect_equal(a$input, "Z_N")
  expect_equal(b$input, "Sigma_GE")
  c <- r_transcca(list(G = Gi, E = Ei), K = 1, sumabs = 0.5)
  expect_equal(rownames(c$V), nm[[2]])
  expect_equal(rownames(c$U), nm[[1]])
  expect_equal(c$input, "G_E")
  ex <- matrix(FALSE, p, g)
  ex[, 1:3] <- TRUE
  e <- r_transcca(list(Z = Z, N = N), K = 1, sumabs = 0.5, exclude = ex)
  expect_true(all(e$V[1:3, ] == 0))
  expect_error(r_transcca(Z, K = 1), "named list")
  utils::capture.output(msg <- testthat::capture_messages(suppressWarnings(
    r_transcca(list(Z = Z, N = N), K = 1, sumabs = 0.5, verbose = TRUE))))
  expect_match(msg[1], "Built W \\(25 variants x 8 genes\\) from input 'Z_N'")
})
