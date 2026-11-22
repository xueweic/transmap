# Internal program constructors used by transmap_local().
#
# Each constructor takes the expression of one module (samples x genes, genes
# already restricted to the module and to non-constant columns) and returns a
# list with:
#   scores      samples x k matrix of program phenotypes
#   loadings    genes x k matrix of gene weights (NULL when not defined)
#   eigenvalues named numeric vector, one value per score column (NA if none)
#   k           number of program phenotypes kept
#   info        list of method-specific diagnostics

# Single-gene module: the program phenotype is the gene itself.
.local_single_gene <- function(E, name, suffix = "") {
  score <- E[, 1, drop = FALSE]
  colnames(score) <- paste0(name, suffix)
  loadings <- matrix(1, 1, 1, dimnames = list(colnames(E), colnames(score)))
  list(scores = score, loadings = loadings,
       eigenvalues = stats::setNames(NA_real_, colnames(score)),
       k = 1L, info = list(single_gene = TRUE))
}

# PC1: first principal component of the standardized module expression.
.local_pc1 <- function(E, name) {
  if (ncol(E) == 1) return(.local_single_gene(E, name))
  pc <- stats::prcomp(E, center = TRUE, scale. = TRUE)
  var_each <- pc$sdev^2 / sum(pc$sdev^2)
  score <- pc$x[, 1, drop = FALSE]
  rownames(score) <- rownames(E)
  colnames(score) <- name
  loadings <- pc$rotation[, 1, drop = FALSE]
  colnames(loadings) <- name
  list(scores = score, loadings = loadings,
       eigenvalues = stats::setNames(pc$sdev[1]^2, name), k = 1L,
       info = list(var_explained = var_each[1]))
}

# aPC: adaptive number of principal components. Components are added while the
# cumulative variance of the previous components is below `max_cumvar` and the
# next component explains more variance than the Marchenko-Pastur upper bound
# (1 + sqrt(p / n))^2 / p.
.local_apc <- function(E, name, max_cumvar = 0.8) {
  if (ncol(E) == 1) return(.local_single_gene(E, name, "_PC1"))
  pc <- stats::prcomp(E, center = TRUE, scale. = TRUE)
  var_each <- pc$sdev^2 / sum(pc$sdev^2)
  cumvar <- cumsum(var_each)
  n <- nrow(E)
  p <- ncol(E)
  mp_upper <- (1 + sqrt(p / n))^2 / p
  k <- 1L
  max_possible <- length(var_each)
  if (max_possible > 1) {
    for (j in 2:max_possible) {
      if (cumvar[j - 1] >= max_cumvar || var_each[j] < mp_upper) break
      k <- j
    }
  }
  cols <- paste0(name, "_PC", seq_len(k))
  score <- pc$x[, seq_len(k), drop = FALSE]
  rownames(score) <- rownames(E)
  colnames(score) <- cols
  loadings <- pc$rotation[, seq_len(k), drop = FALSE]
  colnames(loadings) <- cols
  list(scores = score, loadings = loadings,
       eigenvalues = stats::setNames(pc$sdev[seq_len(k)]^2, cols), k = k,
       info = list(var_explained = var_each[seq_len(k)],
                   mp_upper = mp_upper))
}

# WPC: whitened principal components. The gene-gene correlation expected under
# the null is estimated from null trans z-scores (variants with max |z| < 2,
# relaxed to < 3, else identity). Standardized expression is whitened by the
# pseudo inverse square root of that correlation, and components with
# eigenvalue > 1 are kept.
.local_wpc <- function(E, name, null_z) {
  if (ncol(E) == 1) return(.local_single_gene(E, name, "_WPC1"))
  z <- null_z[, colnames(E), drop = FALSE]
  max_abs <- apply(abs(z), 1, max)
  null_idx <- which(max_abs < 2)
  threshold <- 2
  if (length(null_idx) < ncol(z) + 1) {
    null_idx <- which(max_abs < 3)
    threshold <- 3
  }
  if (length(null_idx) < ncol(z) + 1) {
    Sigma <- diag(ncol(z))
    null_idx <- integer()
    threshold <- NA_real_
  } else {
    Sigma <- suppressWarnings(stats::cor(z[null_idx, , drop = FALSE]))
    if (any(!is.finite(Sigma))) {
      stop("Null z-score correlation is not finite for module ", name,
           "; check `null_z` for constant columns.", call. = FALSE)
    }
    Sigma <- as.matrix(Matrix::nearPD(Sigma, corr = TRUE)$mat)
  }
  whitened <- scale(E, center = TRUE, scale = TRUE) %*%
    .transmap_pseudo_inverse_root(Sigma)
  pca_w <- stats::prcomp(whitened, center = FALSE, scale. = FALSE)
  eigenvalues <- pca_w$sdev^2
  k <- sum(eigenvalues > 1)
  info <- list(null_threshold = threshold, null_variants = length(null_idx))
  if (k == 0) {
    return(list(scores = NULL, loadings = NULL, eigenvalues = NULL, k = 0L,
                info = info))
  }
  cols <- paste0(name, "_WPC", seq_len(k))
  score <- pca_w$x[, seq_len(k), drop = FALSE]
  rownames(score) <- rownames(E)
  colnames(score) <- cols
  loadings <- pca_w$rotation[, seq_len(k), drop = FALSE]
  rownames(loadings) <- colnames(E)
  colnames(loadings) <- cols
  list(scores = score, loadings = loadings,
       eigenvalues = stats::setNames(eigenvalues[seq_len(k)], cols), k = k,
       info = info)
}

# AUCell: area under the recovery curve of each module, computed on gene
# rankings of the full expression matrix. Requires the Bioconductor package
# AUCell.
.local_aucell <- function(Y, modules, auc_max_rank = 0.2, verbose = FALSE) {
  if (!requireNamespace("AUCell", quietly = TRUE)) {
    stop("method = \"AUCell\" requires the Bioconductor package AUCell.",
         call. = FALSE)
  }
  rankings <- AUCell::AUCell_buildRankings(t(Y), plotStats = FALSE,
                                           verbose = verbose)
  auc <- AUCell::AUCell_calcAUC(modules, rankings,
                                aucMaxRank = nrow(rankings) * auc_max_rank,
                                verbose = verbose)
  scores <- t(AUCell::getAUC(auc))
  scores <- scores[rownames(Y), , drop = FALSE]
  storage.mode(scores) <- "double"
  scores
}
