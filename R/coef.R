#' Limiting distribution of the Dickey-Fuller coefficient statistic
#'
#' `coef_cdf()` gives the limit of \eqn{P[n(\hat\phi - \phi_0) \le y]}{P[n(phi^ - phi0) <= y]} for the least squares
#' estimate \eqn{\hat\phi}{phi^} of \eqn{\phi}{phi} in the regression of [ltratio_cdf()], under
#' \eqn{\phi = 1 + c/n}{phi = 1 + c/n}, with \eqn{\phi_0 = 1 + c_0/n}{phi0 = 1 + c0/n} (`c0`, default 0: the
#' coefficient test of Dickey and Fuller, \eqn{n(\hat\phi - 1)}{n(phi^ - 1)}). `coef_icdf()` gives its quantiles.
#'
#' The limit is \eqn{(U - c_0 V)/V}{(U - c0 V)/V}, so that \eqn{P[(U - c_0 V)/V \le y] = P[U - (c_0 + y)V \le 0]}{P[(U
#' - c0 V)/V <= y] = P[U - (c0 + y)V <= 0]}: the distribution function of the t-ratio at 0 with the local null value
#' \eqn{c_0 + y}{c0 + y}. The C engine has no entry point of its own for the coefficient statistic: as `coef_cdf` and
#' `coef_icdf` of the MATLAB programs, `coef_cdf()` evaluates `ltratio_cdf(0, cas, c, c0 = c0 + y, g0, method)` at
#' each `y`, in every case and with every start and option of [ltratio_cdf()], and `coef_icdf()` finds each quantile
#' by the Illinois variant of regula falsi on `coef_cdf()`, on a bracket grown from `[-40, 10]` until its width is
#' below `tol` (about ten evaluations of `coef_cdf()`). With a design, or with `engine` given, the router of the C
#' engine computes (the same construction in C, with the engine of each point; the general engine's own functions
#' where it takes the general engine).
#'
#' @inheritParams ltratio_cdf
#' @param y Numeric vector or array of points. `-Inf`, `Inf` and `NaN` give 0, 1 and `NaN`.
#' @param c0 The local null value \eqn{c_0}{c0} of the statistic \eqn{n(\hat\phi - \phi_0)}{n(phi^ - phi0)}.
#' @return `coef_cdf()`: the distribution function, with the attributes of `y`; `coef_icdf()`: the quantiles, with
#'   the attributes of `p`.
#' @inheritSection ltratio_cdf Cases
#' @inheritSection ltratio_cdf Engines
#' @seealso [ltratio_cdf()], [coef_pdf()], [coef_quantile()].
#' @examples
#' coef_cdf(-14.0936, "DU", -10)    # power of the 5% coefficient test of rho_mu at c = -10
#' coef_icdf(0.05, "DK")            # -8.039, the 5% critical value of rho
#' coef_cdf(-20, dfl_perron("A", 0.5))
#' @export
coef_cdf <- function(y, cas, c = 0, c0 = 0, g0 = 0, method = "auto", verified = FALSE,
                     engine = c("auto", "classical", "general"), start = NULL, gls = NULL) {
  engine <- dfl_engine(engine)
  tg <- dfl_target(cas, g0, start, gls)
  me <- dfl_method(method)
  ve <- dfl_flag(verified, "verified")
  c <- dfl_num(c, "c")
  c0 <- dfl_num(c0, "c0")
  pl <- dfl_plan(tg, engine, c, c0, me, ve)
  if (pl$path == "legacy") return(dfl_once(coef_cdf_(y, pl$case, c, c0, pl$g0, method, verified)))
  if (me != 0L) stop("method names a representation of the classical engine; it does not apply with the general engine", call. = FALSE)
  if (ve == 1L) {
    # (the checked computation of the general engine: the router has none for the coefficient statistic)
    G <- dfl_general_target(tg)
    return(dfl_det_finite("coefcdf", G$D, y, c, dfl_opts(c0, G$start, G$g0, G$gls, TRUE)))
  }
  dfl_once(dfl_shape(dfl_route_call("coef", dfl_real(y, "y"), pl, c(c, c0)), y))
}
coef_cdf_ <- function(y, cas, c, c0, g0, method, verified) {
  dfl_case(cas)
  c <- dfl_num(c, "c")
  c0 <- dfl_num(c0, "c0")
  g0 <- dfl_num(g0, "g0")
  dfl_method(method)
  yv <- dfl_real(y, "y")
  F <- vapply(yv, function(yk) {
    if (is.nan(yk) || is.na(yk)) NaN
    else if (yk == Inf) 1
    else if (yk == -Inf) 0
    else ltratio_cdf(0, cas, c, c0 = c0 + yk, g0 = g0, method = method, verified = verified, engine = "classical")
  }, 0)
  dfl_shape(F, y)
}

#' @rdname coef_cdf
#' @param p Numeric vector or array of probabilities in (0, 1).
#' @param tol The width of the final bracket.
#' @export
coef_icdf <- function(p, cas, c = 0, c0 = 0, g0 = 0, method = "auto", tol = 1e-9, verified = FALSE,
                      engine = c("auto", "classical", "general"), start = NULL, gls = NULL) {
  engine <- dfl_engine(engine)
  tg <- dfl_target(cas, g0, start, gls)
  me <- dfl_method(method)
  ve <- dfl_flag(verified, "verified")
  c <- dfl_num(c, "c")
  c0 <- dfl_num(c0, "c0")
  pl <- dfl_plan(tg, engine, c, c0, me, ve)
  if (pl$path == "legacy") return(dfl_once(coef_icdf_(p, pl$case, c, c0, pl$g0, method, tol, verified)))
  if (me != 0L) stop("method names a representation of the classical engine; it does not apply with the general engine", call. = FALSE)
  pv <- dfl_real(p, "p")
  if (!all(!is.na(pv) & pv > 0 & pv < 1)) stop("p must lie in (0, 1)", call. = FALSE)
  dfl_once(dfl_shape(dfl_route_call("coefq", pv, pl, c(c, c0, dfl_num(tol, "tol"))), p))
}
coef_icdf_ <- function(p, cas, c, c0, g0, method, tol, verified) {
  p0 <- p
  p <- dfl_real(p, "p")
  if (!all(!is.na(p) & p > 0 & p < 1)) stop("p must lie in (0, 1)", call. = FALSE)
  tol <- dfl_num(tol, "tol")
  dfl_case(cas)
  c <- dfl_num(c, "c")
  c0 <- dfl_num(c0, "c0")
  g0 <- dfl_num(g0, "g0")
  dfl_method(method)
  y <- vapply(p, function(pk) {
    f <- function(x) coef_cdf_(x, cas, c, c0, g0, method, verified) - pk
    lo <- -40
    hi <- 10
    flo <- f(lo)
    fhi <- f(hi)
    while (flo > 0) {
      hi <- lo; fhi <- flo; lo <- 2 * lo - 10; flo <- f(lo)
    }
    while (fhi < 0) {
      lo <- hi; flo <- fhi; hi <- 2 * hi + 10; fhi <- f(hi)
    }
    side <- 0
    x <- (lo + hi) / 2
    for (it in 1:200) {
      x <- (lo * fhi - hi * flo) / (fhi - flo)
      fx <- f(x)
      if (fx == 0) break
      if (fx < 0) {
        lo <- x; flo <- fx
        if (side == -1) fhi <- fhi / 2
        side <- -1
      } else {
        hi <- x; fhi <- fx
        if (side == 1) flo <- flo / 2
        side <- 1
      }
      if (hi - lo < tol) {
        x <- (lo * fhi - hi * flo) / (fhi - flo)
        break
      }
    }
    x
  }, 0)
  dfl_shape(y, p0)
}
