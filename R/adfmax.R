#' Limiting distribution of the ADFmax test of Leybourne (1995)
#'
#' `adfmax_cdf()` and `adfmax_pdf()` give the limit of the distribution function and the density of the ADFmax
#' statistic of Leybourne (1995), the larger of the Dickey-Fuller t-ratios of the forward regression (\eqn{y_t} on
#' \eqn{y_{t-1}}) and of the reverse regression (\eqn{y_t} on \eqn{y_{t+1}}), under
#' \eqn{\phi = 1 + c/n}{phi = 1 + c/n}, no drift. In the limit both regressions have the denominator \eqn{V} and the
#' numerators \eqn{U} and \eqn{-U - 1}, so that ADFmax converges to \eqn{(|R| - 1/2)/\sqrt{V}}{(|R| - 1/2)/sqrt(V)},
#' \eqn{R = U + 1/2}, whose distribution is computed by the conditional integral (Section 6 of the paper).
#' `adfmax_icdf()` gives its quantiles.
#'
#' The cases are `"DU"` and `"DT"` (a constant, or a constant and a trend, zero start) and `"SU"` and `"ST"` (the
#' same with a stationary start, `c < 0`). The values are accurate to about 1e-11.
#'
#' @inheritParams ltratio_cdf
#' @param cas The case, one of `"DU"`, `"DT"`, `"SU"`, `"ST"`.
#' @param info If `TRUE`, return a list with the values `value`, the error estimates `err` and the statuses
#'   `status` (nonzero: the error estimate is not met) of the conditional integral.
#' @return A numeric vector with the attributes of `x` (or the list described under `info`). A warning is given if
#'   the conditional integral reports, at some points, that its error estimate is not met.
#' @references Leybourne, S. J. (1995), "Testing for Unit Roots Using Forward and Reverse Dickey-Fuller
#'   Regressions", Oxford Bulletin of Economics and Statistics 57, 559-571.
#' @seealso [ltratio_cdf()].
#' @examples
#' adfmax_cdf(-2.43542, "DU")       # 0.05 (Table I of the paper)
#' adfmax_pdf(c(-3, -2), "DT", -5)
#' @export
adfmax_cdf <- function(x, cas, c = 0, info = FALSE) {
  dfl_adfmax(2L, x, cas, c, info)
}

#' @rdname adfmax_cdf
#' @export
adfmax_pdf <- function(x, cas, c = 0, info = FALSE) {
  dfl_adfmax(1L, x, cas, c, info)
}

dfl_adfmax <- function(kind, x, cas, c, info) {
  k <- dfl_case(cas, c("DU", "DT", "SU", "ST"))
  r <- .Call(C_dfl_r_adfmax, kind, dfl_real(x, "x"), k$model, dfl_num(c, "c"), k$stat)
  names(r) <- c("value", "err", "status")
  if (any(r$status != 0L)) {
    warning(sprintf("the conditional integral did not meet its error estimate at %d of %d points; the result may be inaccurate there",
                    sum(r$status != 0L), length(x)), call. = FALSE)
  }
  r <- lapply(r, dfl_shape, x = x)
  if (isTRUE(info)) r else r$value
}

#' @rdname adfmax_cdf
#' @param p Numeric vector or array of probabilities in (0, 1).
#' @param tol The tolerance on the quantile.
#' @param x0 The starting point of the search.
#' @details `adfmax_icdf()` has no entry point of its own in the C engine: as `adfmax_icdf` of the MATLAB and Julia
#'   programs, it finds each quantile by Newton's method on `adfmax_cdf()` with `adfmax_pdf()` as the derivative,
#'   safeguarded by bisection on a bracket grown from `[x0 - 1, x0 + 1]`.
#' @examples
#' adfmax_icdf(0.05, "DT")          # -3.11891 (Table I of the paper)
#' @export
adfmax_icdf <- function(p, cas, c = 0, tol = 1e-10, x0 = -2) {
  p0 <- p
  p <- dfl_real(p, "p")
  if (!all(!is.na(p) & p > 0 & p < 1)) stop("p must lie in (0, 1)", call. = FALSE)
  dfl_case(cas, c("DU", "DT", "SU", "ST"))
  c <- dfl_num(c, "c")
  tol <- dfl_num(tol, "tol")
  x0 <- dfl_num(x0, "x0")
  cdf <- function(y) adfmax_cdf(y, cas, c)
  q <- vapply(p, function(pk) {
    lo <- x0 - 1
    hi <- x0 + 1
    while (cdf(lo) > pk) lo <- lo - 1
    while (cdf(hi) < pk) hi <- hi + 1
    x <- (lo + hi) / 2
    for (it in 1:100) {
      F <- cdf(x)
      f <- adfmax_pdf(x, cas, c)
      if (F > pk) hi <- x else lo <- x
      xn <- x - (F - pk) / f
      if (!(is.finite(xn) && lo < xn && xn < hi)) xn <- (lo + hi) / 2
      if (abs(xn - x) < tol) return(xn)
      x <- xn
    }
    warning(sprintf("no convergence at p = %g", pk), call. = FALSE)
    x
  }, 0)
  dfl_shape(q, p0)
}
