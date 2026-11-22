#' Construct local program phenotypes from gene modules
#'
#' @description `transmap_local()` summarizes the expression of each gene
#'   module into one or more program phenotypes that can be used as outcomes in
#'   trans-QTL mapping. Modules are typically co-expression modules from
#'   [r_hotspot()], but any named list of gene sets can be used.
#'
#' @details Four constructors are available:
#' \describe{
#'   \item{`"WPC"` (default)}{Whitened principal components. The gene-gene
#'     correlation expected under the null is estimated from null trans
#'     z-scores in `null_z`, using variants whose maximum absolute z-score
#'     across the module genes is below 2 (relaxed to 3, otherwise the
#'     identity matrix is used). Standardized expression is whitened by the
#'     pseudo inverse square root of this correlation, and all components with
#'     eigenvalue greater than 1 are kept. Modules with no such component are
#'     skipped.}
#'   \item{`"PC1"`}{First principal component of the standardized module
#'     expression.}
#'   \item{`"aPC"`}{Adaptive principal components. Components are added while
#'     the cumulative variance explained by the previous components is below
#'     `max_cumvar` and the next component explains more variance than the
#'     Marchenko-Pastur bound \eqn{(1 + \sqrt{p/n})^2 / p}.}
#'   \item{`"AUCell"`}{AUCell enrichment score of each module, computed from
#'     gene rankings of the full expression matrix `Y`. Requires the
#'     Bioconductor package AUCell. AUCell breaks ties in the rankings at
#'     random, so set a seed for reproducible scores.}
#' }
#' A module with a single available gene uses that gene's expression as its
#' program phenotype for WPC, PC1 and aPC.
#'
#' @param Y Numeric matrix of expression, samples in rows and genes in columns,
#'   with sample row names and gene column names.
#' @param module Gene modules: a character vector of genes (one module), a
#'   named list of character vectors, or the result of [r_hotspot()]. Genes not
#'   in `Y` are ignored.
#' @param method Program constructor, one of `"WPC"`, `"PC1"`, `"aPC"` or
#'   `"AUCell"`.
#' @param null_z Null trans z-scores used by WPC: a numeric matrix with
#'   variants in rows and genes in columns, or a list with element `z` holding
#'   that matrix. Required for WPC when a module has more than one gene, and
#'   must contain every module gene present in `Y`.
#' @param covariates Optional numeric matrix of covariates with sample row
#'   names. When supplied, covariates are regressed out of `Y` before programs
#'   are built (and, for AUCell, also out of the resulting scores). Only
#'   samples present in both `Y` and `covariates` are kept.
#' @param max_cumvar Cumulative variance cut-off for aPC.
#' @param auc_max_rank Fraction of the ranked genes used as `aucMaxRank` in
#'   AUCell.
#' @param verbose Logical; print progress messages.
#'
#' @return An object of class `"transmap_local"`, a list with
#' \describe{
#'   \item{`scores`}{Samples by programs matrix of program phenotypes.}
#'   \item{`loadings`}{Named list of gene by component loading matrices, one
#'     per module (`NULL` for AUCell).}
#'   \item{`eigenvalues`}{Named numeric vector, one value per column of
#'     `scores` (`NA` for single-gene modules and AUCell).}
#'   \item{`summary`}{Data frame with one row per module: number of input
#'     genes, genes used, programs kept and method diagnostics.}
#'   \item{`method`}{The method used.}
#' }
#'
#' @examples
#' data(transmap_example)
#' Y <- transmap_example$Y
#' modules <- transmap_example$modules
#'
#' # Default WPC, which needs null trans z-scores
#' res <- transmap_local(Y, modules, null_z = transmap_example$null_z)
#' dim(res$scores)
#' res$summary
#'
#' # First principal component per module
#' res_pc1 <- transmap_local(Y, modules, method = "PC1")
#'
#' @seealso [r_hotspot()] to detect gene modules, [transmap_global()] for
#'   global programs.
#' @export
transmap_local <- function(Y, module,
                           method = c("WPC", "PC1", "aPC", "AUCell"),
                           null_z = NULL,
                           covariates = NULL,
                           max_cumvar = 0.8,
                           auc_max_rank = 0.2,
                           verbose = FALSE) {
  method <- match.arg(method)
  Y <- .transmap_check_matrix(Y, "Y")
  modules <- .transmap_as_module_list(module)
  Y <- .transmap_residualize(Y, covariates, verbose = verbose)
  Y <- .transmap_drop_constant(Y, verbose = verbose)

  genes_in <- lapply(modules, function(g) intersect(g, colnames(Y)))
  n_used <- lengths(genes_in)
  if (all(n_used == 0)) {
    stop("No module genes are present in `Y`.", call. = FALSE)
  }
  if (any(n_used == 0) && verbose) {
    message("Skipping ", sum(n_used == 0),
            " module(s) with no genes in `Y`.")
  }

  if (method == "WPC" && any(n_used > 1)) {
    if (is.list(null_z) && !is.data.frame(null_z)) null_z <- null_z$z
    if (is.null(null_z)) {
      stop("`null_z` is required for method = \"WPC\".", call. = FALSE)
    }
    null_z <- .transmap_check_matrix(null_z, "null_z", need_rownames = FALSE)
    missing_z <- setdiff(unique(unlist(genes_in[n_used > 1])), colnames(null_z))
    if (length(missing_z) > 0) {
      stop("`null_z` is missing ", length(missing_z), " module gene(s): ",
           paste(utils::head(missing_z, 10), collapse = ", "),
           if (length(missing_z) > 10) ", ..." else "", ".", call. = FALSE)
    }
  }

  summary <- data.frame(module = names(modules),
                        n_genes = lengths(modules),
                        n_used = n_used,
                        k = 0L,
                        stringsAsFactors = FALSE)
  rownames(summary) <- NULL

  if (method == "AUCell") {
    keep <- n_used > 0
    scores <- .local_aucell(Y, genes_in[keep], auc_max_rank = auc_max_rank,
                            verbose = verbose)
    if (!is.null(covariates)) scores <- .transmap_residualize(scores, covariates)
    summary$k[summary$module %in% colnames(scores)] <- 1L
    out <- list(scores = scores, loadings = NULL,
                eigenvalues = stats::setNames(rep(NA_real_, ncol(scores)),
                                              colnames(scores)),
                summary = summary, method = method)
    class(out) <- "transmap_local"
    return(out)
  }

  scores <- list()
  loadings <- list()
  eigenvalues <- c()
  diag_info <- vector("list", length(modules))
  for (i in seq_along(modules)) {
    if (n_used[i] == 0) next
    name <- names(modules)[i]
    E <- Y[, genes_in[[i]], drop = FALSE]
    fit <- switch(method,
      WPC = .local_wpc(E, name, null_z),
      PC1 = .local_pc1(E, name),
      aPC = .local_apc(E, name, max_cumvar = max_cumvar)
    )
    summary$k[i] <- fit$k
    diag_info[[i]] <- fit$info
    if (fit$k == 0) {
      if (verbose) message("Module ", name, ": no component kept; skipped.")
      next
    }
    scores[[name]] <- fit$scores
    loadings[[name]] <- fit$loadings
    eigenvalues <- c(eigenvalues, fit$eigenvalues)
  }
  if (length(scores) == 0) {
    stop("No program phenotype could be constructed.", call. = FALSE)
  }
  if (method == "WPC") {
    summary$null_threshold <- vapply(diag_info, function(x) {
      if (is.null(x$null_threshold)) NA_real_ else x$null_threshold
    }, numeric(1))
    summary$null_variants <- vapply(diag_info, function(x) {
      if (is.null(x$null_variants)) NA_integer_ else as.integer(x$null_variants)
    }, integer(1))
  }
  scores <- do.call(cbind, unname(scores))
  rownames(scores) <- rownames(Y)
  if (verbose) {
    message("Constructed ", ncol(scores), " program phenotype(s) from ",
            sum(summary$k > 0), " module(s) with method ", method, ".")
  }
  out <- list(scores = scores, loadings = loadings, eigenvalues = eigenvalues,
              summary = summary, method = method)
  class(out) <- "transmap_local"
  out
}
