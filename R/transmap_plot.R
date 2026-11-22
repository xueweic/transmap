#' Plot observed against null transCCA program strengths
#'
#' @description `transcca_plot()` compares the observed `d2_score` of each
#'   program with its null distribution (from [transcca_null()]) and marks
#'   programs with an empirical p-value at or below `alpha`.
#'
#' @details `type = "bar"` draws the null mean of each program as a bar with
#'   a \eqn{\pm 1.96} SD error bar, the observed score as a red triangle and
#'   a blue star above significant programs. `type = "box"` draws the null
#'   distribution of each program as a box plot, ordered by decreasing null
#'   mean, with the observed score as a red square and a star above
#'   significant programs.
#'
#' @inheritParams transcca_pvalue
#' @param alpha Significance level for the stars.
#' @param type `"bar"` or `"box"`.
#' @param ylab Y-axis label.
#' @param ... Further graphical arguments passed to [graphics::barplot()] or
#'   [graphics::boxplot()].
#'
#' @return Invisibly, a data frame with one row per program (in plotting
#'   order): `program`, `observed`, `null_mean`, `null_sd`, `p_value` and
#'   `significant`.
#'
#' @examples
#' set.seed(1)
#' null <- cbind(Program_1 = rnorm(50, 0.2, 0.05),
#'               Program_2 = rnorm(50, 0.15, 0.05),
#'               Program_3 = rnorm(50, 0.1, 0.05))
#' transcca_plot(c(0.6, 0.18, 0.4), null, alpha = 0.05)
#' transcca_plot(c(0.6, 0.18, 0.4), null, alpha = 0.05, type = "box")
#'
#' @seealso [transcca_null()], [transcca_pvalue()]
#' @export
transcca_plot <- function(observed, null, alpha = 0.005,
                          type = c("bar", "box"),
                          ylab = "Program strength (d2 score)", ...) {
  type <- match.arg(type)
  observed <- .transcca_observed(observed)
  null <- .transcca_check_null(null, observed)
  colnames(null) <- names(observed)
  ord <- seq_along(observed)
  if (type == "box") {
    ord <- order(colMeans(null, na.rm = TRUE), decreasing = TRUE)
  }
  tab <- data.frame(
    program = names(observed)[ord],
    observed = unname(observed[ord]),
    null_mean = unname(colMeans(null, na.rm = TRUE)[ord]),
    null_sd = unname(apply(null, 2, stats::sd, na.rm = TRUE)[ord]),
    p_value = unname(transcca_pvalue(observed, null)[ord]),
    stringsAsFactors = FALSE
  )
  tab$significant <- tab$p_value <= alpha
  sig <- which(tab$significant)

  if (type == "bar") {
    lower <- tab$null_mean - 1.96 * tab$null_sd
    upper <- tab$null_mean + 1.96 * tab$null_sd
    top <- max(upper, tab$observed, na.rm = TRUE)
    ylim <- c(min(0, lower, na.rm = TRUE), top * 1.15)
    mids <- graphics::barplot(tab$null_mean, names.arg = tab$program,
                              col = "grey", border = "grey", ylim = ylim,
                              ylab = ylab, ...)
    graphics::arrows(mids, lower, mids, upper, angle = 90, code = 3,
                     length = 0.05)
    graphics::points(mids, tab$observed, pch = 17, col = "red", cex = 1.5)
    if (length(sig) > 0) {
      star_y <- pmax(upper[sig], tab$observed[sig]) + 0.06 * top
      graphics::text(mids[sig], star_y, "*", col = "blue", cex = 2)
    }
    graphics::legend("topright", legend = c("Observed", "Null mean +/- 1.96 SD"),
                     pch = c(17, 15), col = c("red", "grey"), bty = "n")
  } else {
    null_ord <- null[, ord, drop = FALSE]
    colmax <- apply(rbind(null_ord, tab$observed), 2, max, na.rm = TRUE)
    top <- max(colmax)
    graphics::boxplot(null_ord, col = "grey", pch = 19,
                      ylim = c(min(0, null_ord, na.rm = TRUE), top * 1.15),
                      ylab = ylab, ...)
    graphics::points(seq_along(ord), tab$observed, pch = 15, col = "red",
                     cex = 1.5)
    if (length(sig) > 0) {
      graphics::points(sig, colmax[sig] + 0.06 * top, pch = 8, cex = 1.5)
    }
    graphics::legend("topright", legend = "Observed", pch = 15, col = "red",
                     bty = "n")
  }
  invisible(tab)
}
