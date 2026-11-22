#' Simulated example data for transmap
#'
#' A small simulated data set used in the examples and vignettes.
#'
#' @format A list with elements
#' \describe{
#'   \item{`Y`}{300 x 60 expression matrix (samples x genes). Genes
#'     `GENE01`-`GENE12`, `GENE13`-`GENE24` and `GENE25`-`GENE36` form three
#'     co-expression modules; the other genes are noise.}
#'   \item{`covariates`}{300 x 2 covariate matrix with sample row names.}
#'   \item{`modules`}{Named list of the three simulated gene modules.}
#'   \item{`null_z`}{400 x 60 matrix of null trans z-scores (variants x
#'     genes), used by WPC.}
#'   \item{`z`}{150 x 60 matrix of marginal trans z-scores (variants x genes)
#'     with two global programs, used by [r_transcca()].}
#'   \item{`n`}{Sample size of the z-scores.}
#' }
#' @source Simulated; see `data-raw/transmap_example.R` in the package source.
"transmap_example"
