#' Construct the transCCA input matrix W
#'
#' @description `w_construction()` turns a list of variant-gene association
#'   statistics into the variant-by-gene matrix of standardized effects
#'   \eqn{W} analysed by [r_transcca()]. It is the `W_calculation()` step of
#'   transBoost and takes the same `data_list` format.
#'
#' @details `data_list` is a named list holding exactly one of the following
#' sets of elements (variants in rows, genes in columns; matrices need row and
#' column names):
#' \describe{
#'   \item{`Z`, `N`}{Marginal z-scores and sample size:
#'     \eqn{W = Z / \sqrt{N - 2 + Z^2}}. With `MAF` (and optionally `sdY`),
#'     \eqn{W} is rescaled to the effect scale,
#'     \eqn{W_{ij} \sqrt{2 f_i (1 - f_i)} / s_j}.}
#'   \item{`Beta`}{Standardized effects used as \eqn{W} directly, or rescaled
#'     with `MAF` and `sdY` as above.}
#'   \item{`Beta`, `N` and one of `Pvalue`, `Z`, `seBeta`}{Unstandardized
#'     effects. Each variant's genotype scale is estimated as the mean over
#'     genes of \eqn{s_j / \sqrt{\beta_{ij}^2 + (N - 2) se_{ij}^2}}, using
#'     genes with at least one p-value below `sig_cut` (or, if none, the 20\%
#'     of genes with the smallest p-values); \eqn{W} is `Beta` scaled by it
#'     and divided by `sdY`.}
#'   \item{`log2FC`, `N1`, `N0`}{Differential expression between two groups
#'     of sizes `N1` and `N0`: \eqn{W = \sqrt{N_1 N_0} / (N_1 + N_0) \cdot
#'     log2FC}.}
#'   \item{`G`, `E`}{Individual-level genotypes and expression with the same
#'     samples in rows: \eqn{W = G^T E / N}.}
#'   \item{`Sigma_GE`}{A precomputed variant-gene matrix used as \eqn{W}.}
#' }
#' `N` is a single number or a matrix like `Z`; `MAF` has one value per
#' variant, `sdY` one per gene, and `N1`, `N0` one value or one per variant.
#'
#' The list may also hold `Sigma_GG` (variant correlation) and/or `Sigma_EE`
#' (gene correlation). \eqn{W} is then whitened as
#' \eqn{\Sigma_{GG}^{-1/2} W \Sigma_{EE}^{-1/2}}, using pseudo inverse square
#' roots that keep the leading eigenvalues explaining more than `tol` of the
#' total. With `G` and `E`, `Sigma_GG = "insample"` or `Sigma_EE = "insample"`
#' uses the sample covariance of `G` or `E`. `exclude` is applied before
#' whitening.
#'
#' @param data_list Named list of association statistics; see Details.
#' @param exclude Optional logical matrix with the dimension of `W`; `TRUE`
#'   entries (e.g. cis variant-gene pairs) are set to zero.
#' @param tol Fraction of the eigenvalue sum kept when whitening.
#' @param sig_cut P-value cut-off used to estimate genotype scales from
#'   `Beta`, `N` and `Pvalue`/`Z`/`seBeta`.
#'
#' @return The variants by genes matrix \eqn{W}, with attribute `"input"`
#'   naming the input type used (`"Z_N"`, `"Beta"`, `"Beta_N"`, `"log2FC"`,
#'   `"G_E"` or `"Sigma_GE"`).
#'
#' @examples
#' data(transmap_example)
#' W <- w_construction(list(Z = transmap_example$z, N = transmap_example$n))
#' dim(W)
#' attr(W, "input")
#'
#' @seealso [r_transcca()]
#' @export
w_construction <- function(data_list, exclude = NULL, tol = 0.999,
                           sig_cut = 1e-4) {
  built <- .transcca_W(data_list, exclude = exclude, tol = tol,
                       sig_cut = sig_cut)
  structure(built$W, input = built$input)
}

#' Identify global trans programs with transCCA
#'
#' @description `r_transcca()` finds global trans programs from variant-gene
#'   association statistics. The statistics in `data_list` are converted to
#'   the variant-by-gene matrix \eqn{W} by [w_construction()], and a sparse
#'   rank-K decomposition of \eqn{W} is fitted by penalized matrix
#'   decomposition (PMD, Witten et al. 2009). Each component is a program: a
#'   sparse set of genes (`V`) jointly associated with a sparse set of
#'   variants (`U`). Use [transmap_global()] to turn the gene loadings into
#'   program phenotypes.
#'
#' @details With `different_sparsity = FALSE` (default) one sparsity level,
#'   chosen by cross-validation with [PMA::PMD.cv()] unless `sumabs` is given,
#'   is used for all K components. With `different_sparsity = TRUE`, the
#'   sparsity is re-tuned by cross-validation for every component on the
#'   deflated matrix, and components whose tuned sparsity falls to the lower
#'   bound are set to zero. Cross-validation is random; set a seed for
#'   reproducible results.
#'
#' @param data_list Named list of association statistics, in the format of
#'   transBoost's `data.list`, e.g. `list(Z = z, N = n)`,
#'   `list(Sigma_GE = W)` or `list(G = G, E = E, Sigma_EE = "insample")`. See
#'   [w_construction()] for all accepted inputs.
#' @param K Number of programs. If `NULL`, chosen from the singular values of
#'   `W` with `select_method`.
#' @param different_sparsity Logical; tune the sparsity separately for each
#'   program.
#' @param sumabs Sparsity level passed to [PMA::PMD()] when
#'   `different_sparsity = FALSE`. If `NULL`, chosen by cross-validation.
#' @param select_style How singular values are computed when `K` is `NULL`:
#'   `"svd"` or `"eigen"`.
#' @param select_method How `K` is chosen when `NULL`: `"max"` (largest gap
#'   between consecutive singular values) or `"var_explain"`.
#' @param var_exp Cumulative fraction of singular values used with
#'   `select_method = "var_explain"`.
#' @param null_thres If the largest absolute entry of `W` is below this value
#'   no program is fitted.
#' @param normalization Scale used for the program scores `d2_score`:
#'   `"original"`, `"F_norm"` or `"max"`.
#' @inheritParams w_construction
#' @param verbose Logical; print progress messages.
#' @param ... Further arguments passed to [PMA::PMD()] when
#'   `different_sparsity = FALSE`.
#'
#' @return An object of class `"r_transcca"`, a list with
#' \describe{
#'   \item{`V`}{Genes by programs matrix of gene loadings.}
#'   \item{`U`}{Variants by programs matrix of variant loadings.}
#'   \item{`d`}{Program strengths (PMD factor values).}
#'   \item{`d2_score`}{Squared program strengths after `normalization`.}
#'   \item{`sparsity`}{Sparsity level(s) used.}
#'   \item{`singular_values`}{Singular values of `W` (used to choose `K`).}
#'   \item{`input`}{The input type used to build `W`.}
#' }
#' `V`, `U`, `d` and `d2_score` are `NULL` when `W` falls below `null_thres`.
#'
#' @references Witten DM, Tibshirani R, Hastie T (2009). A penalized matrix
#'   decomposition, with applications to sparse principal components and
#'   canonical correlation analysis. Biostatistics 10(3):515-534.
#'
#' @examples
#' data(transmap_example)
#' data_list <- list(Z = transmap_example$z, N = transmap_example$n)
#' prog <- r_transcca(data_list, K = 3, sumabs = 0.4)
#' dim(prog$V)
#'
#' @seealso [w_construction()], [transmap_global()], [transcca_null()]
#' @export
r_transcca <- function(data_list, K = NULL,
                       different_sparsity = FALSE,
                       sumabs = NULL,
                       select_style = c("svd", "eigen"),
                       select_method = c("max", "var_explain"),
                       var_exp = 2 / 3,
                       null_thres = 0,
                       normalization = c("original", "F_norm", "max"),
                       exclude = NULL, tol = 0.999, sig_cut = 1e-4,
                       verbose = FALSE, ...) {
  select_style <- match.arg(select_style)
  select_method <- match.arg(select_method)
  normalization <- match.arg(normalization)
  built <- .transcca_W(data_list, exclude = exclude, tol = tol,
                       sig_cut = sig_cut)
  W <- built$W
  if (verbose) {
    message("Built W (", nrow(W), " variants x ", ncol(W),
            " genes) from input '", built$input, "'.")
  }
  if (!is.null(K) && (length(K) != 1 || K < 1 || K > min(dim(W)))) {
    stop("`K` must be between 1 and min(dim(W)).", call. = FALSE)
  }

  res <- if (isTRUE(different_sparsity)) {
    if (verbose) message("Fitting PMD with a separate sparsity per program.")
    .transcca_diff_sparsity(W, K = K, select_style = select_style,
                            select_method = select_method, var_exp = var_exp,
                            verbose = verbose)
  } else {
    if (verbose) message("Fitting PMD with one sparsity for all programs.")
    .transcca_same_sparsity(W, K = K, select_style = select_style,
                            select_method = select_method, var_exp = var_exp,
                            null_thres = null_thres, sumabs = sumabs,
                            verbose = verbose, ...)
  }

  out <- list(V = NULL, U = NULL, d = NULL, d2_score = NULL,
              sparsity = res$sparsity, singular_values = res$s_values,
              input = built$input)
  class(out) <- "r_transcca"
  if (is.null(res$v)) return(out)

  U <- as.matrix(res$u)
  V <- as.matrix(res$v)
  d <- as.numeric(res$d)
  K_used <- ncol(V)
  d_norm <- d
  W_k <- W
  for (ii in seq_len(K_used)) {
    u_ii <- U[, ii, drop = FALSE]
    v_ii <- V[, ii, drop = FALSE]
    W_scaled <- switch(normalization,
      F_norm = W_k / sum(W^2),
      max = W_k / max(abs(W)),
      original = W_k
    )
    d_norm[ii] <- as.numeric(t(u_ii) %*% W_scaled %*% v_ii)
    W_k <- W_k - d[ii] * u_ii %*% t(v_ii)
  }
  programs <- paste0("Program_", seq_len(K_used))
  dimnames(V) <- list(colnames(W), programs)
  dimnames(U) <- list(rownames(W), programs)
  out$V <- V
  out$U <- U
  out$d <- stats::setNames(d, programs)
  out$d2_score <- stats::setNames(d_norm^2, programs)
  out
}

#' Construct global program phenotypes
#'
#' @description `transmap_global()` turns global programs into program
#'   phenotypes: for each program, the expression of its genes weighted by the
#'   program gene loadings, \eqn{Y V}. Programs usually come from
#'   [r_transcca()].
#'
#' @param Y Numeric matrix of expression, samples in rows and genes in columns,
#'   with sample row names and gene column names.
#' @param program Gene loadings: the result of [r_transcca()], or a numeric
#'   matrix with genes in rows (gene row names) and programs in columns.
#' @param covariates Optional numeric matrix of covariates with sample row
#'   names, regressed out of `Y` before scoring. Only samples present in both
#'   inputs are kept.
#' @param scale Logical; standardize each gene of `Y` before scoring.
#' @param verbose Logical; print progress messages.
#'
#' @return An object of class `"transmap_global"`, a list with
#' \describe{
#'   \item{`scores`}{Samples by programs matrix of program phenotypes.}
#'   \item{`loadings`}{Gene loadings restricted to the genes used.}
#'   \item{`genes_used`}{Genes present in both `Y` and `program`.}
#'   \item{`genes_missing`}{Program genes absent from `Y`.}
#' }
#'
#' @examples
#' data(transmap_example)
#' prog <- r_transcca(list(Z = transmap_example$z, N = transmap_example$n),
#'                     K = 3, sumabs = 0.4)
#' res <- transmap_global(transmap_example$Y, prog)
#' dim(res$scores)
#'
#' @seealso [r_transcca()], [transmap_local()]
#' @export
transmap_global <- function(Y, program, covariates = NULL, scale = FALSE,
                            verbose = FALSE) {
  Y <- .transmap_check_matrix(Y, "Y")
  V <- if (inherits(program, "r_transcca")) program$V else program
  if (is.null(V)) stop("`program` contains no programs.", call. = FALSE)
  V <- .transmap_check_matrix(V, "program", need_colnames = FALSE)
  if (is.null(colnames(V))) colnames(V) <- paste0("Program_", seq_len(ncol(V)))
  Y <- .transmap_residualize(Y, covariates, verbose = verbose)
  if (isTRUE(scale)) {
    Y <- .transmap_drop_constant(Y, verbose = verbose)
    Y <- scale(Y, center = TRUE, scale = TRUE)
  }
  genes <- intersect(colnames(Y), rownames(V))
  if (length(genes) == 0) {
    stop("No program genes are present in `Y`.", call. = FALSE)
  }
  missing <- setdiff(rownames(V), genes)
  if (length(missing) > 0 && verbose) {
    message(length(missing), " program gene(s) are not in `Y` and are ignored.")
  }
  loadings <- V[genes, , drop = FALSE]
  scores <- Y[, genes, drop = FALSE] %*% loadings
  dimnames(scores) <- list(rownames(Y), colnames(V))
  out <- list(scores = scores, loadings = loadings, genes_used = genes,
              genes_missing = missing)
  class(out) <- "transmap_global"
  out
}

# ---------------------------------------------------------------------------
# transCCA engine (internal). Ported from the transBoost implementation.
# ---------------------------------------------------------------------------

.transcca_singular_values <- function(W, select_style = "svd") {
  if (select_style == "svd") return(svd(W, nu = 0, nv = 0)$d)
  WW <- if (nrow(W) >= ncol(W)) crossprod(W) else tcrossprod(W)
  eigen(WW, symmetric = TRUE, only.values = TRUE)$values
}

.transcca_select_K <- function(s_value, select_method = "max", var_exp = 2 / 3) {
  if (select_method == "max") {
    which.max(-diff(s_value))
  } else {
    which(cumsum(s_value) / sum(s_value) > var_exp)[1]
  }
}

.transcca_min_abs <- function(W) {
  max(1 / sqrt(nrow(W)), 1 / sqrt(ncol(W))) + 0.01
}

.transcca_cv_sumabs <- function(W) {
  min_abs <- .transcca_min_abs(W)
  cv <- PMA::PMD.cv(W, type = "standard", center = FALSE,
                    sumabss = seq(min_abs, 0.8, len = 20), trace = FALSE)
  cv$bestsumabs
}

.transcca_same_sparsity <- function(W, K = NULL, select_style = "svd",
                                    select_method = "max", var_exp = 2 / 3,
                                    null_thres = 0, sumabs = NULL,
                                    verbose = FALSE, ...) {
  s_value <- .transcca_singular_values(W, select_style)
  if (is.null(K)) K <- .transcca_select_K(s_value, select_method, var_exp)
  if (max(abs(W)) < null_thres) {
    return(list(u = NULL, v = NULL, d = NULL, s_values = s_value,
                sparsity = sumabs))
  }
  if (is.null(sumabs)) sumabs <- .transcca_cv_sumabs(W)
  fit <- function() {
    PMA::PMD(x = W, K = K, sumabs = sumabs, type = "standard",
             center = FALSE, ...)
  }
  if (verbose) {
    res <- fit()
  } else {
    utils::capture.output(res <- suppressWarnings(suppressMessages(fit())),
                          type = "output")
  }
  list(u = res$u, v = res$v, d = res$d, s_values = s_value,
       sparsity = sumabs)
}

.transcca_diff_sparsity <- function(W, K = NULL, select_style = "svd",
                                    select_method = "max", var_exp = 2 / 3,
                                    niter = 100, verbose = FALSE) {
  s_value <- .transcca_singular_values(W, select_style)
  if (is.null(K)) K <- .transcca_select_K(s_value, select_method, var_exp)
  min_abs <- .transcca_min_abs(W)
  V_init <- .transcca_init_v(W, K)
  maxW <- max(abs(W))
  W_use <- W
  Ds <- numeric(K)
  Us <- matrix(0, nrow = nrow(W), ncol = K)
  Vs <- matrix(0, nrow = ncol(W), ncol = K)
  sparsity <- numeric(K)
  for (k in seq_len(K)) {
    if (verbose && k %% 10 == 0) message("Working on program K = ", k)
    sumabs <- .transcca_cv_sumabs(W_use) * max(abs(W_use)) / maxW
    if (sumabs <= min_abs) {
      sparsity[k] <- 0
      u <- matrix(0, nrow = nrow(W), ncol = 1)
      V_k <- matrix(0, nrow = ncol(W), ncol = 1)
    } else {
      sparsity[k] <- sumabs
      sumabsu <- sqrt(nrow(W_use)) * sumabs
      sumabsv <- sqrt(ncol(W_use)) * sumabs
      old_V <- stats::rnorm(ncol(W))
      V_k <- V_init[, k, drop = FALSE]
      diff_V <- sum(abs(old_V - V_k))
      u <- matrix(0, nrow = nrow(W), ncol = 1)
      for (iter in seq_len(niter)) {
        if (diff_V <= 1e-7) break
        old_V <- V_k
        argu <- W_use %*% V_k
        su <- .transcca_soft(argu, .transcca_binary_search(argu, sumabsu))
        u <- matrix(su / .transcca_l2n(su), ncol = 1)
        argv <- t(u) %*% W_use
        sv <- .transcca_soft(argv, .transcca_binary_search(argv, sumabsv))
        V_k <- matrix(sv / .transcca_l2n(sv), ncol = 1)
        diff_V <- sum(abs(old_V - V_k))
      }
    }
    Ds[k] <- as.numeric(t(u) %*% (W_use %*% V_k))
    Us[, k] <- u
    Vs[, k] <- V_k
    W_use <- W_use - Ds[k] * u %*% t(V_k)
  }
  list(u = Us, v = Vs, d = Ds, sparsity = sparsity, s_values = s_value)
}

.transcca_soft <- function(x, d) sign(x) * pmax(0, abs(x) - d)

.transcca_l2n <- function(vec) {
  a <- sqrt(sum(vec^2))
  if (a == 0) a <- 0.05
  a
}

# Soft-threshold level giving ||u / ||u||_2||_1 = sumabs.
.transcca_binary_search <- function(argu, sumabs) {
  if (.transcca_l2n(argu) == 0 ||
      sum(abs(argu / .transcca_l2n(argu))) <= sumabs) return(0)
  lam1 <- 0
  lam2 <- max(abs(argu)) - 1e-5
  for (iter in seq_len(149)) {
    su <- .transcca_soft(argu, (lam1 + lam2) / 2)
    if (sum(abs(su / .transcca_l2n(su))) < sumabs) {
      lam2 <- (lam1 + lam2) / 2
    } else {
      lam1 <- (lam1 + lam2) / 2
    }
    if ((lam2 - lam1) < 1e-6) return((lam1 + lam2) / 2)
  }
  warning("Didn't quite converge")
  (lam1 + lam2) / 2
}

.transcca_safe_svd <- function(x) {
  for (i in seq_len(10)) {
    out <- try(svd(x), silent = TRUE)
    if (!inherits(out, "try-error")) return(out)
  }
  svd(matrix(stats::rnorm(nrow(x) * ncol(x)), ncol = ncol(x)))
}

# Initial right singular vectors, as in PMA's CheckPMDV.
.transcca_init_v <- function(W, K) {
  if (ncol(W) > nrow(W)) {
    v <- matrix(t(W) %*% (.transcca_safe_svd(W %*% t(W))$v[, seq_len(K)]),
                ncol = K)
    if (anyNA(v)) v <- matrix(.transcca_safe_svd(W)$v[, seq_len(K)], ncol = K)
    v <- sweep(v, 2, apply(v, 2, .transcca_l2n), "/")
    if (anyNA(v)) stop("some are NA")
  } else {
    v <- matrix(.transcca_safe_svd(t(W) %*% W)$v[, seq_len(K)], ncol = K)
  }
  v
}

#' Simulate null data for transCCA
#'
#' @description `generate_null_data()` simulates genotypes and expression with
#'   no variant-gene association, to build the null distribution of transCCA
#'   program strengths (see [transcca_null()]). Each variant is drawn
#'   independently as Binomial(2, f), with its allele frequency f sampled from
#'   the density of the observed `maf` (restricted to 0.01 < f < 0.5), and
#'   standardized. Expression is drawn from \eqn{N(0, \Sigma_{EE})},
#'   independently of the genotypes.
#'
#' @param n Sample size: one number, or one per view.
#' @param p Number of variants: one number, or one per view.
#' @param maf Observed minor allele frequencies whose density is sampled.
#' @param sigma_ee Gene covariance matrix, or a list of one matrix per view.
#'   Its column names name the simulated genes.
#' @param R Number of independent views (data sets) to simulate.
#'
#' @return A list of `R` data lists, each `list(G = , E = )` in the
#'   `data_list` format of [r_transcca()]: `G` is the `n` by `p` matrix of
#'   standardized genotypes (monomorphic variants are set to 0) and `E` the
#'   `n` by genes expression matrix.
#'
#' @examples
#' data(transmap_example)
#' set.seed(1)
#' null <- generate_null_data(n = 100, p = 50, maf = runif(200, 0.05, 0.5),
#'                            sigma_ee = cor(transmap_example$Y[, 1:20]))
#' str(null[[1]])
#'
#' @seealso [transcca_null()]
#' @export
generate_null_data <- function(n, p, maf, sigma_ee, R = 1) {
  if (!is.numeric(R) || length(R) != 1 || R < 1) {
    stop("`R` must be a positive integer.", call. = FALSE)
  }
  per_view <- function(x, name) {
    if (!(length(x) %in% c(1, R))) {
      stop("Provide one `", name, "` for all views or one per view (R = ",
           R, ").", call. = FALSE)
    }
    rep_len(x, R)
  }
  n <- per_view(n, "n")
  p <- per_view(p, "p")
  if (is.matrix(sigma_ee) || is.data.frame(sigma_ee)) {
    sigma_ee <- list(as.matrix(sigma_ee))
  }
  sigma_ee <- per_view(sigma_ee, "sigma_ee")
  maf <- as.numeric(maf[!is.na(maf)])
  if (length(maf) < 2) {
    stop("`maf` needs at least two non-missing values.", call. = FALSE)
  }
  dens <- stats::density(maf)
  cdf <- cumsum(dens$y) / sum(dens$y)
  lapply(seq_len(R), function(r) {
    .null_data_one(n[r], p[r], dens, cdf, sigma_ee[[r]])
  })
}

# Draw one allele frequency from the density of the observed MAF, rejecting
# values outside (0.01, 0.5).
.null_sample_maf <- function(dens, cdf, max_tries = 10000) {
  for (i in seq_len(max_tries)) {
    f <- dens$x[which.min(abs(cdf - stats::runif(1)))]
    if (f > 0.01 && f < 0.5) return(f)
  }
  stop("Could not sample an allele frequency in (0.01, 0.5) from `maf`.",
       call. = FALSE)
}

.null_data_one <- function(n, p, dens, cdf, sigma_ee) {
  G <- vapply(seq_len(p), function(i) {
    as.numeric(stats::rbinom(n, 2, .null_sample_maf(dens, cdf)))
  }, numeric(n))
  G <- matrix(G, nrow = n)
  G <- scale(G)
  G[is.nan(G)] <- 0
  attr(G, "scaled:center") <- NULL
  attr(G, "scaled:scale") <- NULL
  samples <- paste0("sample", seq_len(n))
  dimnames(G) <- list(samples, paste0("SNP", seq_len(p)))
  E <- MASS::mvrnorm(n = n, mu = rep(0, ncol(sigma_ee)), Sigma = sigma_ee)
  E <- matrix(E, nrow = n)
  genes <- colnames(sigma_ee)
  if (is.null(genes)) genes <- paste0("Gene", seq_len(ncol(sigma_ee)))
  dimnames(E) <- list(samples, genes)
  list(G = G, E = E)
}

#' Null distribution of transCCA program strengths
#'
#' @description `transcca_null()` repeatedly simulates null data with
#'   [generate_null_data()] and fits [r_transcca()] with the same `K` and
#'   `sumabs` as the observed analysis, returning the null `d2_score` of each
#'   program. Use [transcca_pvalue()] to compare observed and null scores, and
#'   [transcca_plot()] to display them.
#'
#' @param n_null Number of null replicates.
#' @inheritParams generate_null_data
#' @param K Number of programs fitted in each replicate.
#' @param sumabs Sparsity level for each replicate; `NULL` re-tunes it by
#'   cross-validation every time (slow).
#' @param method `"diag"` fits the null data as is; `"insample"` also whitens
#'   the expression by its sample covariance (`Sigma_EE = "insample"`).
#' @param verbose Logical; report progress.
#' @param ... Further arguments passed to [r_transcca()], e.g.
#'   `normalization`.
#'
#' @return An `n_null` by `K` matrix of null `d2_score` values.
#'
#' @examples
#' data(transmap_example)
#' set.seed(1)
#' null <- transcca_null(n_null = 5, n = 300, p = 150,
#'                       maf = runif(500, 0.05, 0.5),
#'                       sigma_ee = cor(transmap_example$Y), K = 3,
#'                       sumabs = 0.4)
#' dim(null)
#'
#' @seealso [transcca_pvalue()], [transcca_plot()]
#' @export
transcca_null <- function(n_null, n, p, maf, sigma_ee, K, sumabs = NULL,
                          method = c("diag", "insample"), verbose = FALSE,
                          ...) {
  method <- match.arg(method)
  if (!is.numeric(n_null) || length(n_null) != 1 || n_null < 1) {
    stop("`n_null` must be a positive integer.", call. = FALSE)
  }
  out <- matrix(NA_real_, n_null, K,
                dimnames = list(NULL, paste0("Program_", seq_len(K))))
  for (i in seq_len(n_null)) {
    data_list <- generate_null_data(n, p, maf, sigma_ee)[[1]]
    if (method == "insample") data_list$Sigma_EE <- "insample"
    fit <- r_transcca(data_list, K = K, sumabs = sumabs, ...)
    if (!is.null(fit$d2_score)) out[i, seq_along(fit$d2_score)] <- fit$d2_score
    if (verbose) message("Finished null replicate ", i, " of ", n_null, ".")
  }
  out
}

#' Empirical p-values of transCCA programs
#'
#' @description The p-value of program k is the fraction of null replicates
#'   whose `d2_score` for program k exceeds the observed one.
#'
#' @param observed An [r_transcca()] result, or a numeric vector of observed
#'   `d2_score` values.
#' @param null Null replicates by programs matrix, e.g. from
#'   [transcca_null()].
#'
#' @return Named numeric vector of p-values, one per program.
#'
#' @examples
#' transcca_pvalue(c(0.9, 0.1), cbind(runif(100, 0, 0.5), runif(100, 0, 0.5)))
#'
#' @seealso [transcca_null()], [transcca_plot()]
#' @export
transcca_pvalue <- function(observed, null) {
  observed <- .transcca_observed(observed)
  null <- .transcca_check_null(null, observed)
  stats::setNames(colMeans(sweep(null, 2, observed, ">"), na.rm = TRUE),
                  names(observed))
}

.transcca_observed <- function(observed) {
  if (inherits(observed, "r_transcca")) observed <- observed$d2_score
  if (!is.numeric(observed) || length(observed) == 0) {
    stop("`observed` must be an r_transcca result or numeric d2 scores.",
         call. = FALSE)
  }
  if (is.null(names(observed))) {
    names(observed) <- paste0("Program_", seq_along(observed))
  }
  observed
}

.transcca_check_null <- function(null, observed) {
  null <- as.matrix(null)
  if (!is.numeric(null) || ncol(null) != length(observed)) {
    stop("`null` must be a numeric matrix with one column per program (",
         length(observed), ").", call. = FALSE)
  }
  null
}

# ---------------------------------------------------------------------------
# Internal: building W from data_list (the W_calculation step of transBoost)
# ---------------------------------------------------------------------------

.transcca_list_names <- c("Sigma_GE", "Z", "N", "MAF", "sdY", "Beta",
                          "seBeta", "Pvalue", "G", "E", "log2FC", "N1", "N0",
                          "Sigma_GG", "Sigma_EE")

# Elements accepted by each input type (the required ones first).
.transcca_input_types <- list(
  Sigma_GE = list(required = "Sigma_GE", optional = character()),
  G_E      = list(required = c("G", "E"), optional = character()),
  log2FC   = list(required = c("log2FC", "N1", "N0"),
                  optional = character()),
  Beta_N   = list(required = c("Beta", "N"),
                  optional = c("Z", "Pvalue", "seBeta", "sdY")),
  Z_N      = list(required = c("Z", "N"), optional = c("MAF", "sdY")),
  Beta     = list(required = "Beta", optional = c("MAF", "sdY"))
)

# Decide the input type from the names of the supplied statistics.
.transcca_input_type <- function(supplied) {
  has <- function(x) x %in% supplied
  type <- if (has("Sigma_GE")) {
    "Sigma_GE"
  } else if (has("G") || has("E")) {
    "G_E"
  } else if (has("log2FC")) {
    "log2FC"
  } else if (has("Beta") && has("N")) {
    "Beta_N"
  } else if (has("Z")) {
    "Z_N"
  } else if (has("Beta")) {
    "Beta"
  } else {
    stop("`data_list` must hold association statistics: `Z` and `N`, ",
         "`Beta`, `Beta` and `N` with `Pvalue`/`Z`/`seBeta`, `log2FC` with ",
         "`N1` and `N0`, `G` and `E`, or `Sigma_GE`.", call. = FALSE)
  }
  spec <- .transcca_input_types[[type]]
  missing <- setdiff(spec$required, supplied)
  if (length(missing) > 0) {
    stop("Input type '", type, "' also needs `",
         paste(missing, collapse = "`, `"), "` in `data_list`.",
         call. = FALSE)
  }
  extra <- setdiff(supplied, c(spec$required, spec$optional))
  if (length(extra) > 0) {
    stop("`", paste(extra, collapse = "`, `"), "` cannot be combined with ",
         "input type '", type, "'.", call. = FALSE)
  }
  if (type == "Beta_N" && !any(c("Z", "Pvalue", "seBeta") %in% supplied)) {
    stop("`Beta` with `N` also needs one of `Pvalue`, `Z` or `seBeta`.",
         call. = FALSE)
  }
  if (type %in% c("Z_N", "Beta") && has("sdY") && !has("MAF")) {
    stop("`sdY` is used with `MAF` for input type '", type, "'.",
         call. = FALSE)
  }
  type
}

# Check a matrix of the same dimension as W (or a single number for `N`).
.transcca_check_like <- function(x, ref, name, allow_scalar = FALSE) {
  if (allow_scalar && is.numeric(x) && length(x) == 1) return(x)
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.numeric(x) || !identical(dim(x), dim(ref))) {
    stop("`", name, "` must be ",
         if (allow_scalar) "a single number or " else "",
         "a numeric matrix with the dimension of the statistics.",
         call. = FALSE)
  }
  if (anyNA(x)) stop("`", name, "` contains missing values.", call. = FALSE)
  x
}

.transcca_check_vector <- function(x, len, name) {
  if (!is.numeric(x) || anyNA(x) || !(length(x) %in% c(1, len))) {
    stop("`", name, "` must be a numeric vector of length 1 or ", len, ".",
         call. = FALSE)
  }
  as.numeric(x)
}

.transcca_check_n <- function(n, ref) {
  n <- .transcca_check_like(n, ref, "N", allow_scalar = TRUE)
  if (anyNA(n) || any(n <= 2)) {
    stop("`N` must be greater than 2.", call. = FALSE)
  }
  n
}

# Rescale W from the correlation to the effect scale: rows by the genotype
# standard deviation sqrt(2 MAF (1 - MAF)), columns by 1 / sdY.
.transcca_rescale <- function(W, maf = NULL, sd_y = NULL) {
  if (!is.null(maf)) {
    maf <- .transcca_check_vector(maf, nrow(W), "MAF")
    W <- W * rep_len(sqrt(2 * maf * (1 - maf)), nrow(W))
  }
  if (!is.null(sd_y)) {
    sd_y <- .transcca_check_vector(sd_y, ncol(W), "sdY")
    W <- sweep(W, 2, rep_len(sd_y, ncol(W)), "/")
  }
  W
}

# Genotype scale of each variant from unstandardized summary statistics.
.transcca_beta_scale <- function(beta, n, se_beta = NULL, pvalue = NULL,
                                 z = NULL, sd_y = NULL, sig_cut = 1e-4) {
  if (!is.null(pvalue)) {
    pvalue <- .transcca_check_like(pvalue, beta, "Pvalue")
    z <- sign(beta) * sqrt(stats::qchisq(pvalue, 1, lower.tail = FALSE))
    se_beta <- beta / z
  } else if (!is.null(z)) {
    z <- .transcca_check_like(z, beta, "Z")
    se_beta <- beta / z
    pvalue <- stats::pchisq(z^2, 1, lower.tail = FALSE)
  } else {
    se_beta <- .transcca_check_like(se_beta, beta, "seBeta")
    pvalue <- stats::pchisq((beta / se_beta)^2, 1, lower.tail = FALSE)
  }
  var_y <- if (is.null(sd_y)) 1 else {
    rep(rep_len(.transcca_check_vector(sd_y, ncol(beta), "sdY"),
                ncol(beta))^2, each = nrow(beta))
  }
  scale_ij <- sqrt(var_y / (beta^2 + (n - 2) * se_beta^2))
  scale_ij[is.nan(scale_ij)] <- 0
  sig <- which(colSums(pvalue < sig_cut) != 0)
  if (length(sig) == 0) {
    min_p <- apply(pvalue, 2, min)
    sig <- order(min_p)[seq_len(ceiling(0.2 * length(min_p)))]
  }
  rowMeans(scale_ij[, sig, drop = FALSE])
}

.transcca_sigma <- function(sigma, X, name, size) {
  if (is.character(sigma)) {
    if (!identical(sigma, "insample") || is.null(X)) {
      stop("`", name, "` can be \"insample\" only with `G` and `E`.",
           call. = FALSE)
    }
    return(stats::cov(X))
  }
  if (is.data.frame(sigma)) sigma <- as.matrix(sigma)
  if (!is.numeric(sigma) || !is.matrix(sigma) ||
      any(dim(sigma) != size)) {
    stop("`", name, "` must be a ", size, " x ", size, " numeric matrix.",
         call. = FALSE)
  }
  sigma
}

# Build W from data_list; then mask and whiten.
.transcca_W <- function(data_list, exclude = NULL, tol = 0.999,
                        sig_cut = 1e-4) {
  if (!is.list(data_list) || is.data.frame(data_list) ||
      is.null(names(data_list)) || any(names(data_list) == "")) {
    stop("`data_list` must be a named list, e.g. list(Z = z, N = n).",
         call. = FALSE)
  }
  unknown <- setdiff(names(data_list), .transcca_list_names)
  if (length(unknown) > 0) {
    stop("Unknown element(s) in `data_list`: `",
         paste(unknown, collapse = "`, `"), "`. Allowed: ",
         paste(.transcca_list_names, collapse = ", "), ".", call. = FALSE)
  }
  d <- Filter(Negate(is.null), data_list)
  type <- .transcca_input_type(setdiff(names(d), c("Sigma_GG", "Sigma_EE")))

  W <- switch(type,
    Sigma_GE = .transmap_check_matrix(d$Sigma_GE, "Sigma_GE"),
    G_E = {
      G <- .transmap_check_matrix(d$G, "G", need_rownames = FALSE)
      E <- .transmap_check_matrix(d$E, "E", need_rownames = FALSE)
      if (nrow(G) != nrow(E)) {
        stop("`G` and `E` must have the same samples (rows).", call. = FALSE)
      }
      if (!is.null(rownames(G)) && !is.null(rownames(E)) &&
          !identical(rownames(G), rownames(E))) {
        stop("`G` and `E` must have the same samples in the same order.",
             call. = FALSE)
      }
      crossprod(G, E) / nrow(G)
    },
    log2FC = {
      lfc <- .transmap_check_matrix(d$log2FC, "log2FC")
      n1 <- .transcca_check_vector(d$N1, nrow(lfc), "N1")
      n0 <- .transcca_check_vector(d$N0, nrow(lfc), "N0")
      lfc * rep_len(sqrt(n1 * n0) / (n1 + n0), nrow(lfc))
    },
    Beta_N = {
      beta <- .transmap_check_matrix(d$Beta, "Beta")
      n <- .transcca_check_n(d$N, beta)
      d_g <- .transcca_beta_scale(beta, n, se_beta = d$seBeta,
                                  pvalue = d$Pvalue, z = d$Z, sd_y = d$sdY,
                                  sig_cut = sig_cut)
      .transcca_rescale(beta * d_g, sd_y = d$sdY)
    },
    Z_N = {
      z <- .transmap_check_matrix(d$Z, "Z")
      n <- .transcca_check_n(d$N, z)
      .transcca_rescale(z / sqrt(n - 2 + z^2), maf = d$MAF, sd_y = d$sdY)
    },
    Beta = {
      beta <- .transmap_check_matrix(d$Beta, "Beta")
      .transcca_rescale(beta, maf = d$MAF, sd_y = d$sdY)
    }
  )
  W <- as.matrix(W)

  if (!is.null(exclude)) {
    if (is.data.frame(exclude)) exclude <- as.matrix(exclude)
    if (!identical(dim(exclude), dim(W)) || anyNA(exclude) ||
        !(is.logical(exclude) || all(exclude %in% c(0, 1)))) {
      stop("`exclude` must be a logical matrix with the dimension of `W`.",
           call. = FALSE)
    }
    W[exclude == 1] <- 0
  }

  dn <- dimnames(W)
  if (!is.null(d$Sigma_GG)) {
    S <- .transcca_sigma(d$Sigma_GG, if (type == "G_E") G, "Sigma_GG",
                         nrow(W))
    W <- .transmap_pseudo_inverse_root(S, tol = tol, strict = TRUE) %*% W
  }
  if (!is.null(d$Sigma_EE)) {
    S <- .transcca_sigma(d$Sigma_EE, if (type == "G_E") E, "Sigma_EE",
                         ncol(W))
    W <- W %*% .transmap_pseudo_inverse_root(S, tol = tol, strict = TRUE)
  }
  dimnames(W) <- dn
  list(W = W, input = type)
}
