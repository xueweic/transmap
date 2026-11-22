# Unit tests for the transCCA null simulation, p-values and plot.

# Literal port of transBoost's Generate_NullData_one (transBoost_utils.R).
ref_sample_maf <- function(d, cdf) {
  u <- runif(1)
  v <- d$x[which.min(abs(cdf - u))]
  if (v <= 0.01 || v >= 0.5) return(ref_sample_maf(d, cdf))
  v
}
ref_null_one <- function(n, p, MAF, Sigma_EE) {
  MAF <- as.numeric(MAF[!is.na(MAF)])
  d <- density(MAF)
  cdf <- cumsum(d$y) / sum(d$y)
  G <- lapply(1:p, function(i) rbinom(n, 2, ref_sample_maf(d, cdf)))
  G <- do.call(cbind, G)
  list(G_scaled = apply(G, 2, scale),
       E_trans_scaled = MASS::mvrnorm(n = n, mu = rep(0, ncol(Sigma_EE)),
                                      Sigma = Sigma_EE))
}

set.seed(50)
maf_obs <- c(runif(300, 0.02, 0.5), NA)
S_ee <- cor(transmap_example$Y[, 1:6])

test_that("generate_null_data reproduces transBoost's null simulation", {
  set.seed(51)
  ref <- ref_null_one(80, 15, maf_obs, S_ee)
  set.seed(51)
  out <- generate_null_data(80, 15, maf_obs, S_ee)
  expect_length(out, 1)
  expect_equal(unname(out[[1]]$G), ref$G_scaled)
  expect_equal(unname(out[[1]]$E), unname(ref$E_trans_scaled))
  expect_equal(colnames(out[[1]]$E), colnames(S_ee))
  expect_equal(colnames(out[[1]]$G), paste0("SNP", 1:15))
  expect_identical(rownames(out[[1]]$G), rownames(out[[1]]$E))
})

test_that("generate_null_data simulates several views", {
  set.seed(52)
  out <- generate_null_data(n = c(40, 60), p = 10, maf = maf_obs,
                            sigma_ee = list(S_ee, diag(3)), R = 2)
  expect_length(out, 2)
  expect_equal(dim(out[[1]]$G), c(40, 10))
  expect_equal(dim(out[[2]]$E), c(60, 3))
  expect_equal(colnames(out[[2]]$E), paste0("Gene", 1:3))
  expect_error(generate_null_data(c(40, 50, 60), 10, maf_obs, S_ee, R = 2),
               "one per view")
  expect_error(generate_null_data(40, 10, maf_obs, S_ee, R = 0), "`R`")
  expect_error(generate_null_data(40, 10, c(0.2, NA), S_ee), "two non-missing")
})

test_that(".null_sample_maf rejects frequencies outside (0.01, 0.5)", {
  d <- density(runif(500, 0.02, 0.5))
  cdf <- cumsum(d$y) / sum(d$y)
  set.seed(53)
  f <- replicate(200, transmap:::.null_sample_maf(d, cdf))
  expect_true(all(f > 0.01 & f < 0.5))
  bad <- density(c(0.6, 0.7, 0.8))
  expect_error(transmap:::.null_sample_maf(bad, cumsum(bad$y) / sum(bad$y),
                                           max_tries = 50),
               "Could not sample")
})

test_that(".null_data_one standardizes genotypes and zeroes monomorphic ones", {
  d <- density(rep(0.011, 5) + c(0, 1e-4, 2e-4, 3e-4, 4e-4), bw = 1e-4)
  cdf <- cumsum(d$y) / sum(d$y)
  set.seed(54)
  out <- transmap:::.null_data_one(30, 40, d, cdf, diag(2))
  sds <- apply(out$G, 2, sd)
  expect_true(all(sds == 0 | abs(sds - 1) < 1e-12))
  expect_true(any(sds == 0))
  expect_false(anyNA(out$G))
  expect_null(attr(out$G, "scaled:center"))
})

test_that("transcca_null equals refitting r_transcca on simulated data", {
  set.seed(55)
  null <- transcca_null(n_null = 3, n = 80, p = 15, maf = maf_obs,
                        sigma_ee = S_ee, K = 2, sumabs = 0.6)
  expect_equal(dim(null), c(3, 2))
  expect_equal(colnames(null), c("Program_1", "Program_2"))
  set.seed(55)
  ref <- t(vapply(1:3, function(i) {
    dl <- generate_null_data(80, 15, maf_obs, S_ee)[[1]]
    r_transcca(dl, K = 2, sumabs = 0.6)$d2_score
  }, numeric(2)))
  expect_equal(null, ref, ignore_attr = TRUE)
  # "insample" whitens the expression of each replicate
  set.seed(56)
  ins <- transcca_null(1, 80, 15, maf_obs, S_ee, K = 1, sumabs = 0.6,
                       method = "insample")
  set.seed(56)
  dl <- generate_null_data(80, 15, maf_obs, S_ee)[[1]]
  dl$Sigma_EE <- "insample"
  expect_equal(unname(ins[1, 1]),
               unname(r_transcca(dl, K = 1, sumabs = 0.6)$d2_score))
  expect_message(transcca_null(1, 80, 15, maf_obs, S_ee, K = 1, sumabs = 0.6,
                               verbose = TRUE), "replicate 1 of 1")
  expect_error(transcca_null(0, 80, 15, maf_obs, S_ee, K = 1), "`n_null`")
  # replicates without a fit (below null_thres) are NA
  na <- transcca_null(1, 80, 15, maf_obs, S_ee, K = 1, sumabs = 0.6,
                      null_thres = 10)
  expect_true(is.na(na[1, 1]))
})

test_that("transcca_pvalue is the fraction of null scores above observed", {
  null <- cbind(c(0.1, 0.2, 0.3, 0.4), c(0.5, 0.6, 0.7, 0.8))
  expect_equal(transcca_pvalue(c(0.25, 0.9), null),
               c(Program_1 = 0.5, Program_2 = 0))
  prog <- r_transcca(list(Z = transmap_example$z, N = transmap_example$n),
                     K = 2, sumabs = 0.4)
  expect_equal(transcca_pvalue(prog, null),
               transcca_pvalue(prog$d2_score, null))
  expect_error(transcca_pvalue("a", null), "`observed`")
  expect_error(transcca_pvalue(c(1, 2, 3), null), "one column per program")
})

test_that("transcca_plot draws bar and box plots and returns the summary", {
  set.seed(57)
  null <- cbind(rnorm(50, 0.2, 0.05), rnorm(50, 0.3, 0.05),
                rnorm(50, 0.1, 0.05))
  obs <- c(0.6, 0.25, 0.4)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  bar <- transcca_plot(obs, null, alpha = 0.05)
  expect_equal(bar$program, paste0("Program_", 1:3))
  expect_equal(bar$observed, obs)
  expect_equal(bar$null_mean, unname(colMeans(null)))
  expect_equal(bar$significant, c(TRUE, FALSE, TRUE))
  box <- transcca_plot(obs, null, alpha = 0.05, type = "box")
  expect_equal(box$program, paste0("Program_", c(2, 1, 3)))
  expect_equal(box$p_value, unname(transcca_pvalue(obs, null)[c(2, 1, 3)]))
  # nothing significant: no stars, still returns the table
  none <- transcca_plot(c(0, 0, 0), null, type = "box")
  expect_false(any(none$significant))
  none_bar <- transcca_plot(c(0, 0, 0), null)
  expect_false(any(none_bar$significant))
})
