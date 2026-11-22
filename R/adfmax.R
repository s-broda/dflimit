#' Limiting distribution of the ADFmax test of Leybourne (1995)
#'
#' `adfmax_cdf()` and `adfmax_pdf()` give the limit of the distribution function and the density of the ADFmax
#' statistic of Leybourne (1995), the larger of the Dickey-Fuller t-ratios of the forward regression (\eqn{y_t} on
#' \eqn{y_{t-1}}) and of the reverse regression (\eqn{y_t} on \eqn{y_{t+1}}), under
#' \eqn{\phi = 1 + c/n}{phi = 1 + c/n}, no drift. In the limit both regressions have the denominator \eqn{V} and the
#' numerators \eqn{U} and \eqn{-U - 1}, so that ADFmax converges to \eqn{(|R| - 1/2)/\sqrt{V}}{(|R| - 1/2)/sqrt(V)},
#' \eqn{R = U + 1/2}, whose distribution is computed by the conditional integral (the ADFmax test of the paper).
#' `adfmax_icdf()` gives its quantiles.
#'
#' The cases are `"DU"` and `"DT"` (a constant, or a constant and a trend, zero start) and `"SU"` and `"ST"` (the
#' same with a stationary start, `c <= -1e-4`). The values are accurate to about 1e-11.
#'
#' In place of the case a design of the general engine ([dfl_design]) may be given, with `start = "stationary"` for a
#' stationary start; `engine` chooses between the two engines as in [ltratio_cdf()] (its section Engines): with a case
#' code and no `engine` given the classical engine computes as before; a design of one of the four cases takes the
#' classical engine (with `engine = "auto"` the general engine with a stationary start for `c <= -0.1`, with the values
#' below 1e-3 from the classical engine), every other design the general engine.
#'
#' @inheritParams ltratio_cdf
#' @param cas The case, one of `"DU"`, `"DT"`, `"SU"`, `"ST"`, or a design ([dfl_design]).
#' @param start With a design: `"zero"` (the default, `NULL`) or `"stationary"`.
#' @param info If `TRUE`, return a list with the values `value`, the error estimates `err` and the statuses
#'   `status` (nonzero: the error estimate is not met) of the conditional integral (`NA` where the router computes).
#' @return A numeric vector with the attributes of `x` (or the list described under `info`). A warning is given if
#'   the conditional integral reports, at some points, that its error estimate is not met.
#' @references Leybourne, S. J. (1995), "Testing for Unit Roots Using Forward and Reverse Dickey-Fuller
#'   Regressions", Oxford Bulletin of Economics and Statistics 57, 559-571.
#' @seealso [ltratio_cdf()].
#' @examples
#' adfmax_cdf(-2.43542, "DU")       # 0.05 (the quantiles of the ADFmax statistic in the paper)
#' adfmax_pdf(c(-3, -2), "DT", -5)
#' adfmax_cdf(-3, dfl_perron("A", 0.5))
#' @export
adfmax_cdf <- function(x, cas, c = 0, info = FALSE, engine = c("auto", "classical", "general"), start = NULL) {
  dfl_adfmax(2L, x, cas, c, info, dfl_engine(engine), !missing(engine), start)
}

#' @rdname adfmax_cdf
#' @export
adfmax_pdf <- function(x, cas, c = 0, info = FALSE, engine = c("auto", "classical", "general"), start = NULL) {
  dfl_adfmax(1L, x, cas, c, info, dfl_engine(engine), !missing(engine), start)
}

dfl_adfmax_cases <- c("DU", "DT", "SU", "ST")
# the path of ADFmax: "legacy" (the classical entry point), or the router with the spec (for a design: the design,
# which the C router gives to the classical engine where it is one of the four cases)
dfl_adfmax_plan <- function(cas, engine, given, start) {
  tg <- dfl_target(cas, 0, start, NULL, dfl_adfmax_cases)
  cls <- !is.null(tg$case) && tg$case %in% dfl_adfmax_cases
  if (is.null(tg$D)) {
    if (!given || engine == "classical") return(list(path = "legacy", case = tg$case))
    return(list(path = "router", spec = dfl_spec(tg), engine = dfl_engine_code(engine)))
  }
  if (engine == "classical") {
    if (!cls) dfl_no_classical("ADFmax of the classical engine has the cases DU, DT, SU, ST only")
    return(list(path = "legacy", case = tg$case))
  }
  list(path = "router", spec = dfl_spec(tg, FALSE), engine = dfl_engine_code(engine), tg = tg)
}

dfl_adfmax <- function(kind, x, cas, c, info, engine = "auto", given = FALSE, start = NULL) {
  pl <- dfl_adfmax_plan(cas, engine, given, start)
  if (pl$path == "router") {
    v <- dfl_route_call(if (kind == 2L) "adf" else "adfpdf", dfl_real(x, "x"), pl, dfl_num(c, "c"))
    v <- dfl_shape(v, x)
    if (!isTRUE(info)) return(v)
    na <- dfl_shape(rep(NA_real_, length(v)), x)
    return(list(value = v, err = na, status = dfl_shape(rep(NA_integer_, length(v)), x)))
  }
  k <- dfl_case(pl$case, dfl_adfmax_cases)
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
#'   safeguarded by bisection on a bracket grown from `[x0 - 1, x0 + 1]`; where the general engine computes (a design
#'   other than the four cases, or `engine = "general"`), by its own search from `x0` ([adfmax_quantile()]).
#' @examples
#' adfmax_icdf(0.05, "DT")          # -3.11891 (the quantiles of the ADFmax statistic in the paper)
#' @export
adfmax_icdf <- function(p, cas, c = 0, tol = 1e-10, x0 = -2, engine = c("auto", "classical", "general"), start = NULL) {
  p0 <- p
  p <- dfl_real(p, "p")
  if (!all(!is.na(p) & p > 0 & p < 1)) stop("p must lie in (0, 1)", call. = FALSE)
  given <- !missing(engine)
  engine <- dfl_engine(engine)
  pl <- dfl_adfmax_plan(cas, engine, given, start)
  c <- dfl_num(c, "c")
  tol <- dfl_num(tol, "tol")
  x0 <- dfl_num(x0, "x0")
  if (pl$path == "router") {
    tg <- dfl_target(cas, 0, start, NULL, dfl_adfmax_cases)
    if (engine == "general" || is.null(tg$case) || !tg$case %in% dfl_adfmax_cases) {
      # (the general engine's own search, from x0)
      G <- dfl_general_target(tg)
      return(dfl_shape(dfl_det_call("adfquantile", G$D, p, c(c, x0, tol), dfl_opts(0, G$start, 0, NULL))[[1]], p0))
    }
  }
  cdf <- function(y) adfmax_cdf(y, cas, c, engine = engine, start = start)
  pdf <- function(y) adfmax_pdf(y, cas, c, engine = engine, start = start)
  if (pl$path == "legacy") {
    cas <- pl$case
    cdf <- function(y) adfmax_cdf(y, cas, c)
    pdf <- function(y) adfmax_pdf(y, cas, c)
  }
  q <- vapply(p, function(pk) {
    lo <- x0 - 1
    hi <- x0 + 1
    while (cdf(lo) > pk) lo <- lo - 1
    while (cdf(hi) < pk) hi <- hi + 1
    x <- (lo + hi) / 2
    for (it in 1:100) {
      F <- cdf(x)
      f <- pdf(x)
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
