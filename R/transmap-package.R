#' transmap: Multiscale trans-QTL analysis
#'
#' transmap builds program phenotypes for trans-QTL mapping at two scales.
#'
#' Local programs:
#' * [r_hotspot()] detects co-expression gene modules (R implementation of
#'   Hotspot).
#' * [transmap_local()] turns gene modules into program phenotypes with WPC
#'   (default), PC1, aPC or AUCell.
#'
#' Global programs:
#' * [w_construction()] builds the variant-gene matrix W from a
#'   transBoost-style `data_list`.
#' * [r_transcca()] identifies sparse global programs from W.
#' * [transmap_global()] turns program gene loadings into program phenotypes.
#' * [generate_null_data()], [transcca_null()], [transcca_pvalue()] and
#'   [transcca_plot()] compare programs with a simulated null distribution.
#'
#' All functions take in-memory matrices; reading files, covariate selection
#' and chromosome-specific filtering are left to the calling pipeline.
#'
#' @keywords internal
"_PACKAGE"
