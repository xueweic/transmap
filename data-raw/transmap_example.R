# Simulate the small example data shipped with transmap.
# Run from the package root: Rscript data-raw/transmap_example.R
set.seed(2026)

n_sample <- 300
n_gene <- 60
genes <- sprintf("GENE%02d", seq_len(n_gene))
samples <- sprintf("S%03d", seq_len(n_sample))

# Three co-expression modules of 12 genes each; the remaining genes are noise.
modules <- list(
  module_1 = genes[1:12],
  module_2 = genes[13:24],
  module_3 = genes[25:36]
)
factors <- matrix(stats::rnorm(n_sample * 3), n_sample, 3)
covariates <- cbind(cov1 = stats::rnorm(n_sample), cov2 = stats::rnorm(n_sample))
rownames(covariates) <- samples

Y <- matrix(stats::rnorm(n_sample * n_gene), n_sample, n_gene,
            dimnames = list(samples, genes))
for (m in seq_along(modules)) {
  Y[, modules[[m]]] <- Y[, modules[[m]]] + 1.2 * factors[, m]
}
Y <- Y + covariates %*% matrix(stats::rnorm(2 * n_gene, sd = 0.3), 2, n_gene)

# Null trans z-scores: correlated across genes as the module structure implies.
R <- diag(n_gene)
for (m in modules) {
  idx <- match(m, genes)
  R[idx, idx] <- 0.5
}
diag(R) <- 1
n_null <- 400
null_z <- matrix(stats::rnorm(n_null * n_gene), n_null, n_gene) %*% chol(R)
dimnames(null_z) <- list(sprintf("chr1:%d:A:G", seq_len(n_null) * 1000), genes)

# Marginal trans z-scores with two global programs.
n_var <- 150
z <- matrix(stats::rnorm(n_var * n_gene), n_var, n_gene)
z[1:10, 1:12] <- z[1:10, 1:12] + 4
z[11:20, c(25:36, 40:45)] <- z[11:20, c(25:36, 40:45)] - 4
dimnames(z) <- list(sprintf("chr2:%d:C:T", seq_len(n_var) * 1000), genes)

transmap_example <- list(
  Y = round(Y, 3),
  covariates = round(covariates, 3),
  modules = modules,
  null_z = round(null_z, 3),
  z = round(z, 3),
  n = n_sample
)
save(transmap_example, file = "data/transmap_example.rda", compress = "xz")
