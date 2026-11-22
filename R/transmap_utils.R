# Shared internal helpers for transmap.

# Coerce and validate a numeric matrix with sample row names and feature
# column names.
.transmap_check_matrix <- function(x, name, need_rownames = TRUE,
                                   need_colnames = TRUE) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("`", name, "` must be a numeric matrix.", call. = FALSE)
  }
  if (need_rownames && is.null(rownames(x))) {
    stop("`", name, "` must have row names.", call. = FALSE)
  }
  if (need_colnames && is.null(colnames(x))) {
    stop("`", name, "` must have column names.", call. = FALSE)
  }
  if (need_rownames && anyDuplicated(rownames(x))) {
    stop("`", name, "` has duplicated row names.", call. = FALSE)
  }
  if (need_colnames && anyDuplicated(colnames(x))) {
    stop("`", name, "` has duplicated column names.", call. = FALSE)
  }
  if (anyNA(x)) {
    stop("`", name, "` contains missing values.", call. = FALSE)
  }
  x
}

# Regress covariates (plus an intercept) out of every column of Y. Samples are
# matched by row name; samples missing from either input are dropped.
.transmap_residualize <- function(Y, covariates, verbose = FALSE) {
  if (is.null(covariates)) return(Y)
  covariates <- .transmap_check_matrix(covariates, "covariates",
                                       need_colnames = FALSE)
  common <- intersect(rownames(Y), rownames(covariates))
  if (length(common) < ncol(covariates) + 2) {
    stop("Too few samples shared by `Y` and `covariates` (", length(common),
         ").", call. = FALSE)
  }
  if (length(common) < nrow(Y) && verbose) {
    message("Residualization keeps ", length(common), " of ", nrow(Y),
            " samples shared with `covariates`.")
  }
  Y <- Y[common, , drop = FALSE]
  design <- cbind(1, covariates[common, , drop = FALSE])
  out <- qr.resid(qr(design), Y)
  dimnames(out) <- dimnames(Y)
  out
}

# Drop columns with zero or non-finite variance.
.transmap_drop_constant <- function(Y, verbose = FALSE, what = "genes") {
  v <- apply(Y, 2, stats::var)
  keep <- is.finite(v) & v > 0
  if (!all(keep) && verbose) {
    message("Dropping ", sum(!keep), " ", what, " with zero variance.")
  }
  Y[, keep, drop = FALSE]
}

# Pseudo inverse square root based on the leading eigenvalues that explain at
# least `tol` of the total. Used to whiten genes in WPC.
# `strict = TRUE` uses the `> tol` cut-off of transBoost's W calculation;
# the default `>= tol` is the one used by WPC.
.transmap_pseudo_inverse_root <- function(Sigma, tol = 0.999, strict = FALSE) {
  if (nrow(Sigma) != ncol(Sigma)) stop("Please provide a square matrix!")
  eig <- eigen(Sigma, symmetric = TRUE)
  ratio <- cumsum(eig$values) / sum(eig$values)
  L <- which(if (strict) ratio > tol else ratio >= tol)[1]
  if (is.na(L)) L <- length(eig$values)
  vec <- eig$vectors[, seq_len(L), drop = FALSE]
  vec %*% diag(1 / sqrt(eig$values[seq_len(L)]), nrow = L) %*% t(vec)
}

# Normalise `module` input: a character vector, a list of character vectors,
# or an `r_hotspot` object.
.transmap_as_module_list <- function(module) {
  if (inherits(module, "r_hotspot")) module <- module$modules
  if (is.character(module)) module <- list(module = module)
  if (!is.list(module) || length(module) == 0) {
    stop("`module` must be a character vector of genes, a list of character ",
         "vectors, or an `r_hotspot` result.", call. = FALSE)
  }
  ok <- vapply(module, function(x) is.character(x) || is.null(x), logical(1))
  if (!all(ok)) {
    stop("Every element of `module` must be a character vector of gene names.",
         call. = FALSE)
  }
  if (is.null(names(module))) names(module) <- paste0("module_", seq_along(module))
  empty <- is.na(names(module)) | !nzchar(names(module))
  names(module)[empty] <- paste0("module_", which(empty))
  if (anyDuplicated(names(module))) {
    stop("`module` has duplicated names.", call. = FALSE)
  }
  lapply(module, function(x) unique(as.character(x[!is.na(x)])))
}
