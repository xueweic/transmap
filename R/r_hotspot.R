#' Detect gene modules with Hotspot in R
#'
#' @description `r_hotspot()` is an R implementation of the Hotspot method
#'   (DeTomaso and Yosef, 2021) for pre-standardized data (Hotspot
#'   `model = "none"`). It finds genes whose expression is autocorrelated over
#'   a k-nearest-neighbour graph of samples, groups them into modules from
#'   their pairwise local correlations, and scores each module. The modules
#'   can be passed directly to [transmap_local()].
#'
#' @details The steps follow the Python package `hotspot-sc`:
#' \enumerate{
#'   \item Build an exact k-nearest-neighbour graph of samples from `latent`
#'     (Euclidean distance, each sample excluded from its own neighbours).
#'     With `weighted_graph = TRUE`, edges get Gaussian kernel weights whose
#'     width is the distance to the `ceiling(n_neighbors /
#'     neighborhood_factor)`-th neighbour; otherwise all edges have weight 1.
#'   \item Compute the local autocorrelation \eqn{G}, its Z-score, one-sided
#'     p-value and Benjamini-Hochberg FDR for every gene.
#'   \item For genes with FDR below `fdr`, compute pairwise local
#'     correlation Z-scores.
#'   \item Cluster genes by average linkage on these Z-scores and cut the
#'     tree into modules of at least `min_gene` genes, stopping joins whose
#'     correlation Z-score is not significant at `fdr` (BH-adjusted).
#'   \item Score each module as the first principal component of its
#'     neighbour-smoothed expression.
#' }
#' By default `Y` is standardized (`scale = TRUE`) and the standardized
#' expression is also used as `latent`, as in the TransMap module detection
#' pipeline.
#'
#' @param Y Numeric matrix of expression, samples in rows and genes in columns,
#'   with sample row names and gene column names. Genes must not be constant.
#' @param latent Optional numeric matrix (samples x dimensions) defining the
#'   sample neighbourhoods, with the same samples as `Y`. Defaults to the
#'   (standardized) expression matrix.
#' @param n_neighbors Number of nearest neighbours per sample. Defaults to
#'   half the number of samples, `floor(nrow(Y) / 2)`, as in the TransMap
#'   whole-genome module detection.
#' @param weighted_graph Logical; use Gaussian kernel edge weights.
#' @param neighborhood_factor Kernel width factor used when
#'   `weighted_graph = TRUE`.
#' @param fdr FDR threshold for autocorrelated genes and for joining genes
#'   into modules.
#' @param min_gene Minimum number of genes in a module (at least 2).
#' @param core_only Logical; if `TRUE`, only genes in the module cores are
#'   assigned and ambiguous genes are left unassigned.
#' @param scale Logical; standardize each gene of `Y` before analysis.
#' @param module_prefix Prefix for module names, giving names such as
#'   `"hotspot_module_1"`.
#' @param verbose Logical; print progress messages.
#'
#' @return An object of class `"r_hotspot"`, a list with
#' \describe{
#'   \item{`modules`}{Named list of character vectors, one per module. Pass
#'     this, or the whole object, as `module` to [transmap_local()].}
#'   \item{`assignment`}{Named integer vector mapping each autocorrelated
#'     gene to its module (-1 when unassigned).}
#'   \item{`autocorrelation`}{Data frame with columns `gene`, `C`, `Z`,
#'     `Pval` and `FDR`, sorted by decreasing `Z`.}
#'   \item{`local_correlation_z`}{Pairwise local correlation Z-scores of the
#'     autocorrelated genes.}
#'   \item{`local_correlation_c`}{Scaled pairwise local correlations.}
#'   \item{`module_scores`}{Samples by modules matrix of module scores.}
#'   \item{`linkage`}{Average-linkage merge matrix in the layout of
#'     `scipy.cluster.hierarchy.linkage`, with columns `child_1`, `child_2`
#'     (0-based ids), `height` and `size`.}
#'   \item{`graph`}{The neighbour indices and weights used.}
#' }
#'
#' @references DeTomaso D, Yosef N (2021). Hotspot identifies informative gene
#'   modules across modalities of single-cell genomics. Cell Systems
#'   12(5):446-456.
#'
#' @examples
#' data(transmap_example)
#' hs <- r_hotspot(transmap_example$Y, min_gene = 5)
#' lengths(hs$modules)
#' head(hs$autocorrelation)
#'
#' @seealso [transmap_local()]
#' @export
r_hotspot <- function(Y,
                      latent = NULL,
                      n_neighbors = floor(nrow(Y) / 2),
                      weighted_graph = FALSE,
                      neighborhood_factor = 3,
                      fdr = 0.05,
                      min_gene = 10,
                      core_only = FALSE,
                      scale = TRUE,
                      module_prefix = "hotspot_module",
                      verbose = FALSE) {
  Y <- .transmap_check_matrix(Y, "Y")
  gene_var <- apply(Y, 2, stats::var)
  if (any(!is.finite(gene_var) | gene_var == 0)) {
    stop("`Y` has ", sum(!is.finite(gene_var) | gene_var == 0),
         " gene(s) with zero variance; remove them first.", call. = FALSE)
  }
  if (min_gene < 2) stop("`min_gene` must be at least 2.", call. = FALSE)
  X <- if (isTRUE(scale)) scale(Y, center = TRUE, scale = TRUE) else Y
  attr(X, "scaled:center") <- NULL
  attr(X, "scaled:scale") <- NULL
  if (is.null(latent)) {
    latent <- X
  } else {
    latent <- .transmap_check_matrix(latent, "latent", need_colnames = FALSE)
    if (!identical(rownames(latent), rownames(Y))) {
      if (!setequal(rownames(latent), rownames(Y))) {
        stop("`latent` must have the same samples as `Y`.", call. = FALSE)
      }
      latent <- latent[rownames(Y), , drop = FALSE]
    }
  }

  graph <- .hotspot_knn_graph(latent, n_neighbors = n_neighbors,
                              weighted_graph = weighted_graph,
                              neighborhood_factor = neighborhood_factor)
  auto <- .hotspot_autocorrelation(X, graph)
  sig <- auto$gene[auto$FDR < fdr]
  if (verbose) {
    message(length(sig), " of ", ncol(X), " genes are autocorrelated at FDR < ",
            fdr, ".")
  }

  out <- list(modules = list(), assignment = integer(),
              autocorrelation = auto, local_correlation_z = NULL,
              local_correlation_c = NULL, module_scores = NULL,
              linkage = NULL,
              graph = graph[c("neighbors", "weights")])
  class(out) <- "r_hotspot"
  if (length(sig) < 2) {
    if (verbose) message("Fewer than two autocorrelated genes; no modules.")
    return(out)
  }

  lc <- .hotspot_local_correlation(X[, sig, drop = FALSE], graph)
  mods <- .hotspot_modules(lc$Z, min_gene = min_gene, fdr = fdr,
                           core_only = core_only)
  ids <- sort(unique(mods$assignment[mods$assignment > 0]))
  modules <- lapply(ids, function(id) names(mods$assignment)[mods$assignment == id])
  names(modules) <- paste0(module_prefix, "_", ids)
  if (verbose) message("Found ", length(modules), " module(s).")

  module_scores <- NULL
  if (length(modules) > 0) {
    module_scores <- vapply(modules, function(g) {
      .hotspot_module_score(X[, g, drop = FALSE], graph)
    }, numeric(nrow(X)))
    module_scores <- matrix(module_scores, nrow = nrow(X),
                            dimnames = list(rownames(X), names(modules)))
  }

  out$modules <- modules
  out$assignment <- mods$assignment
  out$local_correlation_z <- lc$Z
  out$local_correlation_c <- lc$C
  out$module_scores <- module_scores
  out$linkage <- mods$linkage
  out
}

# Exact kNN graph excluding self. Returns the directed weight matrix W (row i
# holds the edges from sample i to its neighbours), its symmetrization
# W + t(W), node degrees and the total squared weight after merging reciprocal
# edges (Hotspot's non-redundant weights).
.hotspot_knn_graph <- function(latent, n_neighbors = 30, weighted_graph = FALSE,
                               neighborhood_factor = 3) {
  N <- nrow(latent)
  K <- as.integer(n_neighbors)
  if (K < 1 || K >= N) {
    stop("`n_neighbors` must be between 1 and the number of samples - 1.",
         call. = FALSE)
  }
  nn <- RANN::nn2(latent, latent, k = K + 1)
  ind <- matrix(0L, N, K)
  dist <- matrix(0, N, K)
  for (i in seq_len(N)) {
    keep <- which(nn$nn.idx[i, ] != i)[seq_len(K)]
    ind[i, ] <- nn$nn.idx[i, keep]
    dist[i, ] <- nn$nn.dists[i, keep]
  }
  if (isTRUE(weighted_graph)) {
    radius <- ceiling(K / neighborhood_factor)
    sigma <- dist[, radius]
    sigma[sigma == 0] <- 1
    w <- exp(-dist^2 / sigma^2)
    wnorm <- rowSums(w)
    wnorm[wnorm == 0] <- 1
    w <- w / wnorm
  } else {
    w <- matrix(1, N, K)
  }
  W <- Matrix::sparseMatrix(i = rep(seq_len(N), times = K),
                            j = as.vector(ind), x = as.vector(w),
                            dims = c(N, N))
  Wsym <- W + Matrix::t(W)
  list(neighbors = ind, weights = w, W = W, Wsym = Wsym,
       D = Matrix::rowSums(Wsym),
       Wtot2 = sum(Wsym^2) / 2)
}

# Local autocorrelation of every gene (Hotspot compute_autocorrelations with
# centered = TRUE and the "none" model).
.hotspot_autocorrelation <- function(X, graph) {
  WX <- as.matrix(graph$W %*% X)
  G <- colSums(X * WX)
  Z <- G / sqrt(graph$Wtot2)
  G_max <- colSums(graph$D * X^2) / 2
  pval <- stats::pnorm(Z, lower.tail = FALSE)
  out <- data.frame(gene = colnames(X), C = G / G_max, Z = Z, Pval = pval,
                    FDR = stats::p.adjust(pval, method = "BH"),
                    stringsAsFactors = FALSE)
  out <- out[order(out$Z, decreasing = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Pairwise local correlations (Hotspot compute_local_correlations). For each
# pair the Z-score uses the larger of the two conditional variances.
.hotspot_local_correlation <- function(X, graph) {
  SX <- as.matrix(graph$Wsym %*% X)
  LC <- crossprod(X, SX)
  LC <- (LC + t(LC)) / 2
  eg2 <- colSums(SX^2)
  Z <- LC / sqrt(outer(eg2, eg2, pmax))
  diag(Z) <- 0
  G_max <- colSums(graph$D * X^2) / 2
  C <- LC / (outer(G_max, G_max, "+") / 2)
  diag(C) <- 0
  dimnames(Z) <- dimnames(C) <- list(colnames(X), colnames(X))
  list(Z = Z, C = C)
}

# Module score: first principal component of neighbour-smoothed expression,
# signed so that the mean loading is non-negative.
.hotspot_module_score <- function(Xm, graph, lambda = 0.9) {
  smooth <- as.matrix(graph$Wsym %*% Xm) / graph$D
  smooth <- lambda * smooth + (1 - lambda) * Xm
  pc <- stats::prcomp(smooth, center = TRUE, scale. = FALSE, rank. = 1)
  score <- pc$x[, 1]
  if (mean(pc$rotation[, 1]) < 0) score <- -score
  score
}

# Convert an hclust tree to a scipy-style linkage matrix: rows are merges,
# columns are the two child ids (0-based; leaves 0..N-1, merge r is N + r - 1),
# the merge height and the cluster size.
.hotspot_hclust_to_linkage <- function(hc) {
  N <- length(hc$order)
  M <- N - 1L
  ids <- ifelse(hc$merge < 0, -hc$merge - 1, hc$merge + N - 1)
  ids <- matrix(ids, ncol = 2)
  size <- numeric(M)
  for (r in seq_len(M)) {
    s <- 0
    for (c in ids[r, ]) s <- s + if (c < N) 1 else size[c - N + 1]
    size[r] <- s
  }
  out <- cbind(pmin(ids[, 1], ids[, 2]), pmax(ids[, 1], ids[, 2]),
               hc$height, size)
  dimnames(out) <- list(NULL, c("child_1", "child_2", "height", "size"))
  out
}

# Hotspot compute_modules: average-linkage clustering on (offset - Z) and
# module assignment by walking the merges.
.hotspot_modules <- function(Zmat, min_gene = 10, fdr = 0.05,
                             core_only = FALSE) {
  genes <- rownames(Zmat)
  N <- length(genes)
  Zs <- (Zmat + t(Zmat)) / 2
  allZ <- sort(Zs[upper.tri(Zs)])
  allP <- stats::p.adjust(stats::pnorm(allZ, lower.tail = FALSE), method = "BH")
  ii <- which(allP < fdr)
  z_threshold <- if (length(ii) > 0) allZ[ii[1]] else allZ[length(allZ)] + 1

  dd <- Zmat
  diag(dd) <- 0
  offset <- max(dd[upper.tri(dd)])
  hc <- stats::hclust(stats::as.dist(offset - dd), method = "average")
  L <- .hotspot_hclust_to_linkage(hc)
  M <- nrow(L)
  child <- function(id) if (id < N) NA_integer_ else as.integer(id - N + 1)

  # mean merge height of the sub-tree under each merge
  mean_dists <- numeric(M)
  n_merges <- numeric(M)
  for (r in seq_len(M)) {
    a <- child(L[r, 1])
    b <- child(L[r, 2])
    am <- if (is.na(a)) 0 else n_merges[a]
    bm <- if (is.na(b)) 0 else n_merges[b]
    aa <- if (is.na(a)) 0 else mean_dists[a]
    ba <- if (is.na(b)) 0 else mean_dists[b]
    n_merges[r] <- am + bm + 1
    mean_dists[r] <- (aa * am + ba * bm + L[r, 3]) / n_merges[r]
  }

  labels <- rep(-1, M)
  clust_i <- 0
  for (r in seq_len(M)) {
    a <- child(L[r, 1])
    b <- child(L[r, 2])
    na <- if (is.na(a)) 1 else L[a, 4]
    nb <- if (is.na(b)) 1 else L[b, 4]
    la <- if (is.na(a)) -1 else labels[a]
    lb <- if (is.na(b)) -1 else labels[b]
    too_far <- L[r, 3] > offset - z_threshold
    both_big <- na >= min_gene && nb >= min_gene
    if (!core_only && too_far) {
      new <- -1
    } else if (both_big) {
      if (core_only) {
        new <- -1
      } else {
        new <- if (mean_dists[a] >= mean_dists[b]) la else lb
      }
    } else if (core_only && too_far) {
      new <- -1
    } else if (na >= min_gene) {
      new <- la
    } else if (nb >= min_gene) {
      new <- lb
    } else if (na + nb >= min_gene) {
      new <- clust_i
      clust_i <- clust_i + 1
    } else {
      new <- -1
    }
    labels[r] <- new
  }

  # propagate labels from the root down to the leaves
  parent <- integer(M)
  for (r in seq_len(M)) {
    for (c in L[r, 1:2]) {
      ci <- child(c)
      if (!is.na(ci)) parent[ci] <- r
    }
  }
  eff <- numeric(M)
  eff[M] <- labels[M]
  if (M > 1) {
    for (r in (M - 1):1) {
      p <- eff[parent[r]]
      if (core_only) {
        eff[r] <- if (p != -1) p else labels[r]
      } else {
        eff[r] <- if (labels[r] != -1) labels[r] else p
      }
    }
  }
  leaf_label <- numeric(N)
  for (r in seq_len(M)) {
    for (c in L[r, 1:2]) if (c < N) leaf_label[c + 1] <- eff[r]
  }
  ids <- sort(unique(leaf_label[leaf_label != -1]))
  assignment <- ifelse(leaf_label == -1, -1L, match(leaf_label, ids))
  assignment <- stats::setNames(as.integer(assignment), genes)
  list(assignment = assignment, linkage = L, z_threshold = z_threshold,
       offset = offset)
}
