#' Limiting distribution of the Dickey-Fuller coefficient statistic
#'
#' `coef_cdf()` gives the limit of \eqn{P[n(\hat\phi - \phi_0) \le y]}{P[n(phi^ - phi0) <= y]} for the least squares
#' estimate \eqn{\hat\phi}{phi^} of \eqn{\phi}{phi} in the regression of [ltratio_cdf()], under
#' \eqn{\phi = 1 + c/n}{phi = 1 + c/n}, with \eqn{\phi_0 = 1 + c_0/n}{phi0 = 1 + c0/n} (`c0`, default 0: the
#' coefficient test of Dickey and Fuller, \eqn{n(\hat\phi - 1)}{n(phi^ - 1)}). `coef_icdf()` gives its quantiles and
#' `coef_mom()` its moments; [coef_pdf()] its density.
#'
#' The limit is \eqn{(U - c_0 V)/V}{(U - c0 V)/V}, so that \eqn{P[(U - c_0 V)/V \le y] = P[U - (c_0 + y)V \le 0]}{P[(U
#' - c0 V)/V <= y] = P[U - (c0 + y)V <= 0]}.
#'
#' In the cases DK, DU, DT (zero or fixed start `g0`) and SK, SU, ST the default (`method = "auto"` with `engine =
#' "auto"`, up to `c = 20`) is the single inversion of the joint characteristic function along the ray of the half-plane,
#' the method of Yong Bao and Raymond Kan, "Exact Finite-Sample Distribution and Moments of the Least Squares Estimator in
#' First-Order Autoregressions: O(1) Algorithms": adaptive quadrature over decades, the tail completed by integration by
#' parts with the derivatives from jets, the closed-form transform with the branch fixed at each point, one evaluation
#' with an error estimate per value (the C engine's `dfl_engine_coef_x`, as `coef_cdf` of DFLimit.jl and the MATLAB
#' programs). `tol` is its absolute tolerance, `err = TRUE` returns the estimates as well, and a warning names the points
#' at which they exceed it. `method = "inversion"` forces it (also above `c = 20`); `method = "engine"` is the computation
#' of before, the distribution function of the t-ratio at 0 with the local null value \eqn{c_0 + y}{c0 + y}
#' (`ltratio_cdf(0, cas, c, c0 = c0 + y, g0)`, with its choice of engine and representation), which reproduces the
#' values of earlier versions; the t-ratio's representations (`"single"` or `"line"`, `"double"` or `"conditional"`,
#' `"contour"`), `engine = "classical"` or `"general"` and case DG take that computation too. `verified = TRUE` keeps the
#' values of the single inversion whose estimate meets `tol` and recomputes the others by the engines' checked
#' computation. `coef_icdf()` finds each quantile by the Illinois variant of regula falsi to the width `tol`: by default
#' on the single inversion's distribution function (tolerance 1e-12) from the bracket \eqn{E[U - cV]/E[V]}{E[U - cV]/E[V]}
#' -/+ twice the scale of the ratio (`err = TRUE`: the estimate of the distribution function at the quantile divided by
#' the density there), otherwise on the engines' distribution function from a bracket grown from `[-40, 10]`.
#'
#' @inheritParams ltratio_cdf
#' @param y Numeric vector or array of points. `-Inf`, `Inf` and `NaN` give 0, 1 and `NaN`.
#' @param c0 The local null value \eqn{c_0}{c0} of the statistic \eqn{n(\hat\phi - \phi_0)}{n(phi^ - phi0)}.
#' @param method `"auto"` (the default), `"inversion"` (the single inversion of Bao and Kan), `"engine"` (the engines'
#'   computation of before), or a representation of the t-ratio (`"single"` or `"line"`, `"double"` or
#'   `"conditional"`, `"contour"`), which takes the engines.
#' @param tol `coef_cdf()`: the absolute tolerance of each value of the single inversion (1e-12); `coef_icdf()`: the width
#'   of the final bracket; `coef_mom()`: the relative tolerance of the quadrature (1e-10).
#' @param err `TRUE`: a list of the values (`value`) and their error estimates (`err`; `NaN` where the engines computed).
#' @return `coef_cdf()`: the distribution function, with the attributes of `y`; `coef_icdf()`: the quantiles, with the
#'   attributes of `p`; `coef_mom()`: the moments \eqn{E[\{(U - c_0 V)/V\}^j]}{E[{(U - c0 V)/V}^j]}, \eqn{j = 1..p};
#'   with `err = TRUE` a list (`value`, `err`).
#' @inheritSection ltratio_cdf Cases
#' @inheritSection ltratio_cdf Engines
#' @seealso [ltratio_cdf()], [coef_pdf()], [coef_quantile()].
#' @examples
#' coef_cdf(-14.0936, "DU", -10)    # power of the 5% coefficient test of rho_mu at c = -10
#' coef_icdf(0.05, "DK")            # -8.039, the 5% critical value of rho
#' coef_mom(1, "DK")                # -1.7814: the mean of n(phi^ - 1) at the unit root
#' coef_cdf(-20, dfl_perron("A", 0.5))
#' @export
coef_cdf <- function(y, cas, c = 0, c0 = 0, g0 = 0, method = "auto", verified = FALSE,
                     engine = c("auto", "classical", "general"), start = NULL, gls = NULL, tol = 1e-12, err = FALSE) {
  engine <- dfl_engine(engine)
  tg <- dfl_target(cas, g0, start, gls)
  im <- dfl_imethod(method, "coef_cdf", TRUE)
  ve <- dfl_flag(verified, "verified")
  c <- dfl_num(c, "c")
  c0 <- dfl_num(c0, "c0")
  if (im %in% 1:4 || (engine == "classical" && im != 5L)) {
    # (a representation named, or the classical engine: the computation of before)
    me <- if (im %in% 1:4) im else 0L
    pl <- dfl_plan(tg, engine, c, c0, me, ve)
    if (pl$path == "legacy") {
      v <- dfl_once(coef_cdf_(y, pl$case, c, c0, pl$g0, if (im %in% 1:4) method else "auto", verified))
      return(dfl_with_err(v, rep(NaN, length(v)), y, err))
    }
    if (me != 0L) stop("method names a representation of the classical engine; it does not apply with the general engine", call. = FALSE)
  }
  pl <- dfl_plan(tg, engine, c, c0)
  if (pl$path == "legacy") pl <- list(path = "router", spec = dfl_spec(tg, TRUE), engine = dfl_engine_code(engine))
  # (method = "inversion": the engine as given, so that the router names the cases of the single inversion)
  if (im == 5L) pl$engine <- dfl_engine_code(engine)
  r <- dfl_once(dfl_route_call2("coef", dfl_real(y, "y"), pl, c(c, c0, if (im %in% 1:4) 0L else im, dfl_tol(tol), ve)))
  dfl_with_err(r$v, r$e, y, err)
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
#' @param p `coef_icdf()`: numeric vector or array of probabilities in (0, 1); `coef_mom()`: the number of moments.
#' @export
coef_icdf <- function(p, cas, c = 0, c0 = 0, g0 = 0, method = "auto", tol = 1e-9, verified = FALSE,
                      engine = c("auto", "classical", "general"), start = NULL, gls = NULL, err = FALSE) {
  engine <- dfl_engine(engine)
  tg <- dfl_target(cas, g0, start, gls)
  im <- dfl_imethod(method, "coef_icdf", TRUE)
  ve <- dfl_flag(verified, "verified")
  c <- dfl_num(c, "c")
  c0 <- dfl_num(c0, "c0")
  if (im %in% 1:4 || (engine == "classical" && im != 5L)) {
    me <- if (im %in% 1:4) im else 0L
    pl <- dfl_plan(tg, engine, c, c0, me, ve)
    if (pl$path == "legacy") {
      v <- dfl_once(coef_icdf_(p, pl$case, c, c0, pl$g0, if (im %in% 1:4) method else "auto", tol, verified))
      return(dfl_with_err(v, rep(NaN, length(v)), p, err))
    }
    if (me != 0L) stop("method names a representation of the classical engine; it does not apply with the general engine", call. = FALSE)
  }
  pl <- dfl_plan(tg, engine, c, c0)
  if (pl$path == "legacy") pl <- list(path = "router", spec = dfl_spec(tg, TRUE), engine = dfl_engine_code(engine))
  # (method = "inversion": the engine as given, so that the router names the cases of the single inversion)
  if (im == 5L) pl$engine <- dfl_engine_code(engine)
  pv <- dfl_real(p, "p")
  if (!all(!is.na(pv) & pv > 0 & pv < 1)) stop("p must lie in (0, 1)", call. = FALSE)
  r <- dfl_once(dfl_route_call2("coefq", pv, pl, c(c, c0, dfl_num(tol, "tol"), if (im %in% 1:4) 0L else im, ve)))
  dfl_with_err(r$v, r$e, p, err)
}

#' @rdname coef_cdf
#' @export
coef_mom <- function(p, cas, c = 0, c0 = 0, g0 = 0, tol = 1e-10, err = FALSE, engine = c("auto", "classical", "general"),
                     method = "auto", start = NULL, gls = NULL) {
  if (!is.numeric(p) || length(p) != 1L || is.na(p) || p != round(p) || p < 1) stop("p must be a positive integer", call. = FALSE)
  engine <- dfl_engine(engine)
  tg <- dfl_target(cas, g0, start, gls)
  im <- dfl_imethod(method, "coef_mom")
  if (im == 6L) stop("coef_mom: method must be \"auto\" (the jets of Bao and Kan; the engines have no moments of the coefficient statistic)", call. = FALSE)
  c <- dfl_num(c, "c")
  c0 <- dfl_num(c0, "c0")
  pl <- list(path = "router", spec = dfl_spec(tg, TRUE), engine = dfl_engine_code(engine))
  r <- dfl_route_call2("coefmom", 0, pl, c(c, c0, p, im, dfl_tol(tol)))
  if (isTRUE(err)) list(value = r$v, err = r$e) else r$v
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
