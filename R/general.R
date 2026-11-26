#' The t-ratio by the general engine
#'
#' The distribution function, density, quantiles and moments of the limit of the Dickey-Fuller t-ratio for any design
#' of deterministic terms ([dfl_design]), computed by the general engine: the joint moment generating function of
#' \eqn{(U, V)} as a small determinant and the inversion formula on paths around its singularities in the first
#' argument (the inversion on paths around the roots), the same representation for every design; and, for the designs
#' whose numerator has rank one or two, the line integral with the explicit kernel (see `method`).
#'
#' The statistic is the t-ratio of the hypothesis \eqn{\phi = \phi_0 = 1 + c_0/n}{phi = phi0 = 1 + c0/n} in the
#' regression of \eqn{y_t} on the regressors \eqn{d(t/n)} of the design and \eqn{y_{t-1}}, under
#' \eqn{\phi = 1 + c/n}{phi = 1 + c/n}; its limit is \eqn{(U - c_0 V)/\sqrt{V}}{(U - c0 V)/sqrt(V)} with
#' \eqn{U = \int_0^1 \tilde J\,dJ}{U = int J~ dJ} and \eqn{V = \int_0^1 \tilde J^2\,dr}{V = int J~^2 dr},
#' \eqn{J} the Ornstein-Uhlenbeck process \eqn{dJ = cJ\,dr + dW}{dJ = cJ dr + dW} and \eqn{\tilde J}{J~} its residual
#' on the design (by least squares, or with `gls` by GLS detrending at \eqn{\bar\alpha = 1 + \bar c/n}{abar = 1 +
#' cbar/n}, as Elliott, Rothenberg and Stock (1996) detrend; with [dfl_detrend_first()] the statistic that detrends
#' first). The start is \eqn{J(0) = 0} (`start = "zero"`), \eqn{J(0) = \gamma_0}{J(0) = gamma0} (`"fixed"`, the
#' initial observation \eqn{y_0 = \gamma_0 \sigma \sqrt{n}}{y_0 = gamma0 sigma sqrt(n)}) or a draw from the
#' stationary law \eqn{N(0, -1/(2c))} (`"stationary"`, `c < 0`). [ltratio_cdf()] with a design is the same
#' computation where the classical engine does not take it (see its section Engines); the functions here always use
#' the general engine.
#'
#' By default (`verified = FALSE`) every value is one evaluation; `verified = TRUE` evaluates it a second time with
#' other rules on a steeper path, and a third time with doubled rules where the two differ by more than 1e-9, with a
#' warning where the last two still differ. One evaluation of the inversion on paths takes the rule of `rule`: by
#' default (`"auto"`) the outer hyperbola with 1.5 times the step and the paths around the roots with 16 nodes per panel
#' and 32 on the tails of their rays (24 and 40 before November 2026; 32 and 48 for 48 and 64 with a fixed start) where
#' that holds, about half the time, and by regions of the start, the deterministic terms, `c` and `c0` the rule of before
#' or a finer one (near the unit root with `c0 < 0`; the rules and regions of DFLimit.jl, within 7.6e-10 of the rule of
#' before in the distribution function and 4e-9 of the peak in the density on a dense grid); `verified = TRUE` takes the
#' rule of before. The points of one call share the paths of the inversion, so a vector of
#' points is much cheaper in one call than point by point. The engine is checked for `c` from -480 to 10 (above the
#' unit root the coefficient statistic, the moments and the t-ratio centred at `c0 = c` up to 16 with deterministic terms,
#' 10 with them and a fixed start, and 20 without them; beyond those bounds every function stops with an error); designs
#' whose roots run far out on the vertical line (above the unit root, and some unusual designs) can need several
#' gigabytes of memory; where the paths around the roots of one node of the outer contour would need more than a million
#' nodes the engine stops with an error (for instance the design \eqn{(1 + r, 1\{r \le 0.4\})}{(1 + r, 1{r <= 0.4})},
#' without a proper intercept, at every `c`). The engine warns where its estimate of the rounding error of a value is
#' large.
#'
#' `dist_quantile()` finds each quantile by Newton's method from `x0` (the safeguarded search where it does not
#' converge); [ltratio_icdf()] with a design starts it at the normal approximation from the moments. `dist_moments()`
#' gives \eqn{E[\tau^j]}{E[tau^j]}, \eqn{j = 1, \dots, p}{j = 1..p}, by the trapezoidal rule in the logarithm of the
#' argument of an integral of the moment generating function (not for the statistic that detrends first).
#'
#' @param x Numeric vector or array of points. `-Inf`, `Inf` and `NaN` give 0, 1 (density 0) and `NaN`.
#' @param p `dist_quantile()`: numeric vector or array of probabilities in (0, 1). `dist_moments()`: the number of
#'   moments.
#' @param D A design ([dfl_design]).
#' @param c The local-to-unity parameter, \eqn{\phi = 1 + c/n}{phi = 1 + c/n}.
#' @param c0 The local null value, \eqn{\phi_0 = 1 + c_0/n}{phi0 = 1 + c0/n}.
#' @param start The start: `"zero"`, `"fixed"` (with `gamma0`) or `"stationary"` (`c < 0`); `NULL` (the default) is
#'   `"fixed"` where `gamma0` is not 0 and `"zero"` otherwise.
#' @param gamma0 The fixed start \eqn{J(0) = \gamma_0}{J(0) = gamma0}.
#' @param gls `NULL` (the default) for least squares detrending, or \eqn{\bar c < 0}{cbar < 0} for GLS detrending at
#'   \eqn{1 + \bar c/n}{1 + cbar/n}.
#' @param verified `FALSE` (the default): one evaluation; `TRUE`: the checked computation (see Details).
#' @param method `"auto"` (the default), `"line"` or `"paths"`: the line integral with the explicit kernel (the
#'   classical engine's line integral with \eqn{P(t)} and the two roots supplied by the general engine's transform) or
#'   the inversion on paths around the roots. The line integral covers the designs whose numerator has rank one or
#'   two: no deterministic terms, continuous designs that do not vanish at the origin and whose span is closed under
#'   differentiation (polynomial trends of any order, an intercept with sines and cosines), and [dfl_detrend_first()]
#'   of a continuous design that does not vanish at the origin, with a zero or a stationary start and least squares
#'   detrending. `"auto"` takes it at `c0 = 0` for `c >= -20`, where it is 1.1 to 1.7 times faster than the paths per
#'   value for `-20 <= c <= 0` (0.9 to 1.3 times per quantile; at `c = 2` the paths are 1.35 to 2.1 times faster for designs
#'   with deterministic terms; C engine, one thread, November 2026), and the paths otherwise; `"line"` stops with the reason where it does not
#'   cover the design, the start or the options (GLS detrending, a fixed start, a jump).
#' @param rule The rule of one evaluation of the inversion on paths: `"auto"` (the default, by regions), `"present"`
#'   (the rule before November 2026, also taken by `verified = TRUE`), `"thin"` (1.5 times the outer step, 16 and 32
#'   nodes), `"inner"` (the thin nodes on the present step), `"outer"` (the present nodes on the longer step), `"fine"`
#'   (the thin nodes on half the present step) or `"finestep"` (a fixed start without deterministic terms: the present
#'   nodes on the step min(1/2, 2.5/(gamma0 |c0|)) times the present one where the outer hyperbola is flatter than 70
#'   degrees). The line integral (see `method`) is not affected.
#' @param x0 The starting point of the search.
#' @param tol The tolerance of the search.
#' @param info `dist_moments()`: if `TRUE`, a list with the moments `m` and the differences `err` to the rule of
#'   twice the step.
#' @return `dist_cdf()`, `dist_pdf()`, `dist_quantile()`: a numeric vector with the attributes of `x` or `p`;
#'   `dist_moments()`: a numeric vector of length `p` (or the list described under `info`).
#' @references Elliott, G., Rothenberg, T. J. and Stock, J. H. (1996), "Efficient Tests for an Autoregressive Unit
#'   Root", Econometrica 64, 813-836.
#' @seealso [dfl_design], [ltratio_cdf()], [coef_pdf()], [lin_cdf()], [local_power()], [pt_cdf()].
#' @examples
#' D <- dfl_perron("C", 0.5)               # a break in level and slope at mid-sample
#' dist_cdf(c(-5, -4.5, -4), D)            # the left tail at the unit root
#' dist_quantile(0.05, D, x0 = -4.5)       # the 5% critical value
#' dist_cdf(-4, D, -10, start = "stationary")
#' dist_moments(2, dfl_trend(1), -5)
#' @export
dist_cdf <- function(x, D, c = 0, c0 = 0, start = NULL, gamma0 = 0, gls = NULL, verified = FALSE,
                     method = c("auto", "line", "paths"), rule = "auto") {
  dfl_det_finite("cdf", D, x, dfl_num(c, "c"), dfl_opts(c0, start, gamma0, gls, verified, method[1], rule))
}

#' @rdname dist_cdf
#' @export
dist_pdf <- function(x, D, c = 0, c0 = 0, start = NULL, gamma0 = 0, gls = NULL, verified = FALSE,
                     method = c("auto", "line", "paths"), rule = "auto") {
  dfl_det_finite("pdf", D, x, dfl_num(c, "c"), dfl_opts(c0, start, gamma0, gls, verified, method[1], rule), hi = 0)
}

#' @rdname dist_cdf
#' @export
dist_quantile <- function(p, D, c = 0, c0 = 0, start = NULL, gamma0 = 0, gls = NULL, x0 = -2, tol = 1e-9,
                          method = c("auto", "line", "paths"), rule = "auto") {
  pv <- dfl_prob(p)
  dfl_shape(dfl_det_call("quantile", D, pv, c(dfl_num(c, "c"), dfl_num(x0, "x0"), dfl_num(tol, "tol")),
                         dfl_opts(c0, start, gamma0, gls, method = method[1], rule = rule))[[1]], p)
}

#' @rdname dist_cdf
#' @export
dist_moments <- function(p, D, c = 0, c0 = 0, start = NULL, gamma0 = 0, gls = NULL, info = FALSE) {
  p <- dfl_int(p, "p")
  if (p < 1L) stop("p must be a positive integer", call. = FALSE)
  r <- dfl_det_call("moments", D, 0, c(dfl_num(c, "c"), p), dfl_opts(c0, start, gamma0, gls))
  if (isTRUE(info)) list(m = r[[1]], err = r[[2]]) else r[[1]]
}

#' The coefficient statistic by the general engine
#'
#' The limit \eqn{(U - c_0 V)/V}{(U - c0 V)/V} of the coefficient statistic \eqn{n(\hat\phi - \phi_0)}{n(phi^ -
#' phi0)} for any design ([dfl_design]), by the general engine (the notation of [dist_cdf()]):
#' `coef_pdf()` its density, `coef_quantile()` its quantiles (Newton's method from `y0`), `coef_cdf_gp()` its
#' distribution function by the inversion formula of Gil-Pelaez, a check of the contour integral. The distribution
#' function by the inversion on paths around the roots is [coef_cdf()] with a design (or with `engine = "general"`).
#'
#' `coef_pdf()` also takes a case code (`"DK"`, `"DU"`, `"DT"` with the fixed start `gamma0`, `"SK"`, `"SU"`, `"ST"`,
#' `"DG"`) in place of `D`. In those six classical cases, and for a design that is one of them, the default (`method =
#' "auto"` with `engine = "auto"`, up to `c = 20`) is the single inversion of Yong Bao and Raymond Kan, "Exact
#' Finite-Sample Distribution and Moments of the Least Squares Estimator in First-Order Autoregressions: O(1)
#' Algorithms": the density of the ratio by the formula of Geary, one integral of the derivative of the joint
#' characteristic function along the ray, from jets, the tail completed by integration by parts. `tol` is its tolerance
#' relative to the scale \eqn{\max(1, E[V]/\mathrm{sd}(U))}{max(1, E[V]/sd(U))} of the density, `err = TRUE` returns the
#' estimates as well, and a warning names the points at which they exceed it (`verified = TRUE` recomputes those by the
#' general engine's checked computation). `method = "inversion"` forces it, `method = "engine"` (and `engine =
#' "general"`, case DG, the other designs) is the general engine's density, the kernel of its contour integral. Near
#' the point \eqn{y = -c_0}{y = -c0}, where the integrand decays slowly, one value of the single inversion can take a
#' second, with a warning.
#'
#' `coef_cdf_gp()` is by default (a zero or stationary start, up to `c = 5`) the single inversion of Bao and Kan for the
#' general engine's transform, with the branch of the square root of its determinant fixed at each point from the
#' pivots of the precision matrix and the tail completed by integration by parts; `tol` its absolute tolerance, `err =
#' TRUE` the error estimates. `method = "grid"` (and `N` or `umax` given, a fixed start, `c` above 5) is the trapezoidal
#' rule of before in \eqn{u = \sinh v}{u = sinh(v)} with `N` steps (60000) up to `umax` (4000), refused where the grid
#' does not resolve the singularity of the transform nearest the real line.
#'
#' @inheritParams dist_cdf
#' @param D A design ([dfl_design]); for `coef_pdf()` also a case code.
#' @param y Numeric vector or array of points. `-Inf`, `Inf` and `NaN` give 0, 1 (density 0) and `NaN`.
#' @param p Numeric vector or array of probabilities in (0, 1).
#' @param y0 The starting point of the search.
#' @param N,umax The number of steps and the end of the fixed grid of the inversion formula of Gil-Pelaez (`NULL`:
#'   60000 and 4000 where the grid computes; given, they choose the grid).
#' @param method `coef_pdf()`: `"auto"`, `"inversion"` (the single inversion of Bao and Kan) or `"engine"` (the general
#'   engine); `coef_cdf_gp()`: `"auto"`, `"inversion"` or `"grid"`.
#' @param engine `coef_pdf()`: `"auto"` (the default), `"classical"` (the six classical cases, by the single inversion)
#'   or `"general"`.
#' @param tol The tolerance of the single inversion (1e-12).
#' @param err `TRUE`: a list of the values (`value`) and their error estimates (`err`; `NaN` where the general engine
#'   computed).
#' @return A numeric vector with the attributes of `y` or `p` (with `err = TRUE` a list, `value` and `err`).
#' @seealso [coef_cdf()], [dist_cdf()].
#' @examples
#' D <- dfl_perron("A", 0.3)
#' coef_pdf(c(-20, -10), D)
#' coef_quantile(0.05, D)
#' coef_cdf_gp(-12, dfl_trend(0))           # = coef_cdf(-12, "DU")
#' @export
coef_pdf <- function(y, D, c = 0, c0 = 0, start = NULL, gamma0 = 0, gls = NULL, verified = FALSE, method = "auto",
                     engine = c("auto", "classical", "general"), tol = 1e-12, err = FALSE) {
  engine <- dfl_engine(engine)
  tg <- dfl_target(D, gamma0, start, gls)
  im <- dfl_imethod(method, "coef_pdf")
  ve <- dfl_flag(verified, "verified")
  c <- dfl_num(c, "c")
  c0 <- dfl_num(c0, "c0")
  pl <- dfl_plan(tg, engine, c, c0)
  if (pl$path == "legacy") pl <- list(path = "router", spec = dfl_spec(tg, TRUE), engine = dfl_engine_code(engine))
  # (the engine as given: the router decides between the single inversion and the general engine)
  pl$engine <- dfl_engine_code(engine)
  yv <- dfl_real(y, "y")
  out <- rep(NaN, length(yv))
  out[!is.na(yv) & is.infinite(yv)] <- 0
  e <- ifelse(is.finite(yv), NaN, 0)
  ok <- is.finite(yv)
  if (any(ok)) {
    r <- dfl_route_call2("coefpdf", yv[ok], pl, c(c, c0, im, dfl_tol(tol), ve))
    out[ok] <- r$v
    e[ok] <- r$e
  }
  dfl_with_err(out, e, y, err)
}

#' @rdname coef_pdf
#' @export
coef_quantile <- function(p, D, c = 0, c0 = 0, start = NULL, gamma0 = 0, gls = NULL, y0 = -10, tol = 1e-9) {
  pv <- dfl_prob(p)
  dfl_shape(dfl_det_call("coefquantile", D, pv, c(dfl_num(c, "c"), dfl_num(y0, "y0"), dfl_num(tol, "tol")),
                         dfl_opts(c0, start, gamma0, gls))[[1]], p)
}

#' @rdname coef_pdf
#' @export
coef_cdf_gp <- function(y, D, c = 0, c0 = 0, start = NULL, gamma0 = 0, gls = NULL, N = NULL, umax = NULL, method = "auto",
                        tol = 1e-12, err = FALSE) {
  m <- dfl_ptmethod(method)
  r <- dfl_det_finite2("coefgp", D, y, c(dfl_num(c, "c"), dfl_ptgrid(N, umax), m, dfl_tol(tol)), dfl_opts(c0, start, gamma0, gls))
  if (isTRUE(err)) r else r$value
}

#' Local power by the general engine
#'
#' The asymptotic local power at each `c` of the level-`level` test of \eqn{\phi = 1 + c_0/n}{phi = 1 + c0/n} that
#' rejects for small values of the t-ratio (`statistic = "tratio"`) or of the coefficient statistic (`"coef"`), for
#' any design ([dfl_design]): the critical value is the `level`-quantile of the limit at `c = c0`
#' ([dist_quantile()], [coef_quantile()]), the power its distribution function at `c` ([dist_cdf()], [coef_cdf()]).
#'
#' @inheritParams dist_cdf
#' @param c Numeric vector or array of the local-to-unity parameters at which the power is computed.
#' @param level The level of the test.
#' @param statistic `"tratio"` (the default) or `"coef"`.
#' @return The power, with the attributes of `c`.
#' @seealso [dist_cdf()], [pt_envelope()], [power_envelope()].
#' @examples
#' local_power(c(-5, -10), dfl_perron("C", 0.5))
#' local_power(-10, dfl_trend(1), statistic = "coef")
#' @export
local_power <- function(c, D, level = 0.05, statistic = c("tratio", "coef"), c0 = 0, start = NULL, gamma0 = 0,
                        gls = NULL, verified = FALSE) {
  statistic <- match.arg(statistic)
  cv <- dfl_real(c, "c")
  if (!all(is.finite(cv))) stop("c must be finite", call. = FALSE)
  dfl_shape(dfl_det_call("power", D, cv, c(dfl_num(level, "level"), if (statistic == "coef") 1 else 0),
                         dfl_opts(c0, start, gamma0, gls, verified))[[1]], c)
}

#' The linear form U - aV by the general engine
#'
#' `lin_cdf()` and `lin_pdf()` give the limit of \eqn{P(U - aV \le b)}{P(U - aV <= b)} and its density for any design
#' ([dfl_design]; \eqn{U} and \eqn{V} as in [dist_cdf()]), `lin_quantile()` its quantiles (Newton's method from
#' `b0`), and `v_cdf()` and `v_pdf()` the law of \eqn{V} (`a = NaN` in `lin_cdf()`), by the inversion on paths around
#' the roots of the general engine, which takes \eqn{U - aV \le b}{U - aV <= b} as the event below a parabola in
#' \eqn{\sqrt{V}}{sqrt(V)}. [ltratio_lincdf()] with a design is the same where the classical engine does not take it.
#' The law of \eqn{V} is 0 at `v <= 0` (\eqn{V > 0}).
#'
#' @inheritParams dist_cdf
#' @param b Numeric vector or array of points. `-Inf`, `Inf` and `NaN` give 0, 1 (density 0) and `NaN`.
#' @param a The coefficient of \eqn{V}, a single number; `NaN` for the law of \eqn{V}.
#' @param p Numeric vector or array of probabilities in (0, 1).
#' @param b0 The starting point of the search.
#' @param v Numeric vector or array of points of \eqn{V}.
#' @return A numeric vector with the attributes of `b`, `p` or `v`.
#' @seealso [ltratio_lincdf()], [pt_cdf()], [dist_cdf()].
#' @examples
#' D <- dfl_trend(1)
#' lin_cdf(c(-0.6, 0, 0.5), -7.3, D, gls = -13.5)   # the trend case GLS-detrended at -13.5
#' lin_quantile(0.05, -3.5, dfl_trend(-1))         # -0.2639: the point-optimal test against -7
#' v_cdf(c(0.02, 0.2), dfl_perron("B", 0.5), -5)
#' @export
lin_cdf <- function(b, a, D, c = 0, start = NULL, gamma0 = 0, gls = NULL, verified = FALSE) {
  dfl_det_finite("lincdf", D, b, c(dfl_a(a), dfl_num(c, "c")), dfl_opts(0, start, gamma0, gls, verified))
}

#' @rdname lin_cdf
#' @export
lin_pdf <- function(b, a, D, c = 0, start = NULL, gamma0 = 0, gls = NULL, verified = FALSE) {
  dfl_det_finite("linpdf", D, b, c(dfl_a(a), dfl_num(c, "c")), dfl_opts(0, start, gamma0, gls, verified), hi = 0)
}

#' @rdname lin_cdf
#' @export
lin_quantile <- function(p, a, D, c = 0, start = NULL, gamma0 = 0, gls = NULL, b0 = 0, tol = 1e-10) {
  pv <- dfl_prob(p)
  a <- dfl_num(a, "a")
  dfl_shape(dfl_det_call("linquantile", D, pv, c(a, dfl_num(c, "c"), dfl_num(b0, "b0"), dfl_num(tol, "tol")),
                         dfl_opts(0, start, gamma0, gls))[[1]], p)
}

#' @rdname lin_cdf
#' @export
v_cdf <- function(v, D, c = 0, start = NULL, gamma0 = 0, gls = NULL) {
  lin_cdf(v, NaN, D, c, start, gamma0, gls)
}

#' @rdname lin_cdf
#' @export
v_pdf <- function(v, D, c = 0, start = NULL, gamma0 = 0, gls = NULL) {
  lin_pdf(v, NaN, D, c, start, gamma0, gls)
}

#' ADFmax and regions by the general engine
#'
#' `adfmax_quantile()` gives the quantiles of the limit of the ADFmax statistic of Leybourne (1995) for any design
#' ([dfl_design]), by Newton's method from `x0` on the distribution function of the general engine (which
#' [adfmax_cdf()] and [adfmax_pdf()] give with a design or with `engine = "general"`). `region_prob()` gives
#' \eqn{P\{U \le g(V)\}}{P{U <= g(V)}} for the boundaries that the general engine has: the line
#' \eqn{g(v) = b + a v} (`bnd = list(type = "line", a = , b = )`) and the parabola of the t-ratio
#' \eqn{g(v) = x \sqrt{v} + c_0 v}{g(v) = x sqrt(v) + c0 v} (`list(type = "parabola", x = , c0 = 0)`).
#'
#' @inheritParams dist_cdf
#' @param p Numeric vector or array of probabilities in (0, 1).
#' @param bnd The boundary, a list (see Details).
#' @return `adfmax_quantile()`: the quantiles, with the attributes of `p`; `region_prob()`: the probability.
#' @references Leybourne, S. J. (1995), "Testing for Unit Roots Using Forward and Reverse Dickey-Fuller
#'   Regressions", Oxford Bulletin of Economics and Statistics 57, 559-571.
#' @seealso [adfmax_cdf()], [ltratio_region()].
#' @examples
#' adfmax_quantile(0.05, dfl_perron("A", 0.5))
#' region_prob(list(type = "line", a = -3.5, b = -0.2639), dfl_trend(-1))   # 0.05
#' @export
adfmax_quantile <- function(p, D, c = 0, start = NULL, gamma0 = 0, x0 = -2, tol = 1e-10) {
  pv <- dfl_prob(p)
  dfl_shape(dfl_det_call("adfquantile", D, pv, c(dfl_num(c, "c"), dfl_num(x0, "x0"), dfl_num(tol, "tol")),
                         dfl_opts(0, start, gamma0, NULL))[[1]], p)
}

#' @rdname adfmax_quantile
#' @export
region_prob <- function(bnd, D, c = 0, start = NULL, gamma0 = 0, gls = NULL, verified = FALSE) {
  if (!is.list(bnd) || is.null(bnd$type)) stop("bnd must be a list with the element type", call. = FALSE)
  type <- match(tolower(bnd$type), c("line", "parabola"))
  if (is.na(type)) {
    stop("the boundary type must be \"line\" or \"parabola\" (the kernels of the general engine are in closed form for boundaries quadratic in sqrt(V) only)",
         call. = FALSE)
  }
  el <- function(name, default = NULL) {
    v <- bnd[[name]]
    if (is.null(v)) {
      if (is.null(default)) stop(sprintf("a boundary of type %s needs the element %s", bnd$type, name), call. = FALSE)
      v <- default
    }
    dfl_num(v, paste0("bnd$", name))
  }
  par <- if (type == 1L) c(el("a"), el("b")) else c(el("x"), el("c0", 0))
  dfl_det_call("region", D, 0, c(type, par, dfl_num(c, "c")), dfl_opts(0, start, gamma0, gls, verified))[[1]]
}

# probabilities in (0, 1), as a double vector
dfl_prob <- function(p) {
  pv <- dfl_real(p, "p")
  if (!all(!is.na(pv) & pv > 0 & pv < 1)) stop("p must lie in (0, 1)", call. = FALSE)
  pv
}
# the coefficient a of the linear form (NaN: the law of V)
dfl_a <- function(a) {
  if (!is.numeric(a) || length(a) != 1L || (is.na(a) && !is.nan(a))) {
    stop("a must be a single number (NaN: the distribution of V)", call. = FALSE)
  }
  as.double(a)
}
