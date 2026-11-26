#' Limiting distribution of the Dickey-Fuller t-ratio
#'
#' `ltratio_cdf()` and `ltratio_pdf()` give the distribution function and the density of the limit, as the sample
#' size grows, of the least squares t-ratio of the hypothesis \eqn{\phi = \phi_0}{phi = phi0} in the Gaussian
#' first-order autoregression \eqn{y_t = \alpha + \delta t + \phi y_{t-1} + \sigma e_t}{y_t = alpha + delta t +
#' phi y_{t-1} + sigma e_t} under the local alternatives \eqn{\phi = 1 + c/n}{phi = 1 + c/n},
#' \eqn{\phi_0 = 1 + c_0/n}{phi0 = 1 + c0/n}, with no drift (\eqn{\alpha = \delta = 0}{alpha = delta = 0}).
#' With `c0 = 0` (the default) these are the Dickey-Fuller tests of a unit root; `c = 0` is the null hypothesis.
#'
#' The limit is \eqn{\tau = (U - c_0 V)/\sqrt{V}}{tau = (U - c0 V)/sqrt(V)}, with \eqn{U = \int_0^1 \tilde J\,dJ}{U
#' = int J~ dJ} and \eqn{V = \int_0^1 \tilde J^2\,dr}{V = int J~^2 dr}, \eqn{J} the Ornstein-Uhlenbeck process
#' \eqn{dJ = cJ\,dr + dW}{dJ = cJ dr + dW} and \eqn{\tilde J}{J~} its residual on the deterministic terms of the
#' regression.
#'
#' @section Cases:
#' `cas` is a two-letter code, the initial observation (`D`: deterministic, i.e. fixed; `S`: drawn from the
#' stationary distribution, which requires `c <= -1e-4`: closer to zero the variance \eqn{-1/(2c)} of the start is
#' too large for the quadrature, and the functions raise an error) followed by the deterministic terms of the regression
#' (`K`: none, the intercept is known to be zero; `U`: an unknown intercept; `T`: an intercept and a trend):
#' \describe{
#'   \item{`"DK"`, `"DU"`, `"DT"`}{the statistics \eqn{\tau}{tau}, \eqn{\tau_\mu}{tau_mu},
#'     \eqn{\tau_\tau}{tau_tau}, with the start \eqn{y_0 = g_0 \sigma \sqrt{n}}{y_0 = g0 sigma sqrt(n)}
#'     (`g0`, default 0);}
#'   \item{`"SK"`, `"SU"`, `"ST"`}{the same regressions with a stationary start;}
#'   \item{`"DG"`}{the GLS-detrended trend case of Elliott, Rothenberg and Stock (1996), detrended at
#'     \eqn{\bar c = -13.5}{cbar = -13.5} (the DF-GLS tests): zero start, the line integral only, no moments.}
#' }
#'
#' @section Methods and ranges:
#' Three exact representations are implemented: the line integral representation of the companion paper
#' (`method = "single"`, or `"line"`), its conditional integral representation (`"double"`, or `"conditional"`) and,
#' for a fixed start `g0 != 0` in cases DU and DT only, the contour integral representation (`"contour"`): a double
#' integral from the inversion formula on paths of the lead paper, around the singularities of the moment generating
#' function, with the integral over
#' v in closed form and the paths in s and t deformed, an explicit integrand where the line integral needs a noncentral
#' kernel that is itself an integral.
#'
#' By default (`verified = FALSE`) every value is one evaluation by one method: no second evaluation, no refinement
#' driven by an error estimate, no fallback to another representation and no recomputation. `"auto"` (the default)
#' takes at each point the fastest representation that one evaluation gives to well below 1e-5: the contour integral
#' with a fixed start in cases DU, DT (the conditional integral there for `c > 1` with `c0 < 0` and, in case DU, for
#' `c > 1.5` with `c0 > 3`, where the contour integral fails), the line integral with a fixed start in case DK, at
#' `c0 = 0`, with a zero start and `c0 > 0` (but the conditional integral at `x < 0` for `c >= 0` in cases DU, DT), in
#' case DK with `c0 < 0` for `c > -5`, with a stationary start in cases
#' SU, ST and in case SK for `c > -0.5`, and far from the unit root (`c < -100` with a zero start, `c < -50` with a
#' stationary start), and the conditional integral otherwise. Where the parameters lie in a regime in which one evaluation is
#' known to be brittle a warning names it (once per call). `verified = TRUE` gives the checked computation, at 1 to
#' 4.5 times the cost of the default (about 1.8 with the contour integral; more in the tails below 1e-3, which it
#' recomputes): the conditional integral with the agreement of two paths and the refinement of their rules,
#' the contour integral as two evaluations on different paths refined until they agree, a warning where the error
#' estimate of `"double"` or `"contour"` is not met, and with `"auto"` the rule before (the conditional integral with
#' `c0 != 0` and a zero start except for `c0 > 0` at `x >= 0`), the flagged values given another representation's
#' value and the values below 1e-3 recomputed by the conditional integral. A computed density below zero is returned
#' as zero.
#'
#' `c` must be at least -550 (the range in which the programs are checked); the rules are set for `c` from about
#' -480 to 10. With a fixed start `g0 != 0` the functions are checked, and an error is raised outside, in case DK
#' for `c <= 2`, `|g0| <= 2` (`|g0| <= 1` for `c > 0`), in case DU for `-100 <= c <= 2`, `|g0| <= 1.5` and
#' `(c g0)^2/2 <= 4050`, in case DT for `-100 <= c <= 2`, `|g0| <= 3` and `(c g0)^2/2 <= 7200`.
#'
#' In case DK with `c0 > 0` the density is infinite at `x = -sqrt(2 c0)`, and `Inf` is returned there.
#'
#' In two regions the line integral would need more quadrature nodes than the engine allows: with `c0 < 0` near zero
#' and `x > 0`, where a root of the parabola \eqn{X(v) = c_0 v^2 + x v + 1/2}{X(v) = c0 v^2 + x v + 1/2} lies far out
#' (near `x/|c0|`), more than two million nodes on one path of its inner integrals (from `c0 = -1e-4` at `x = 0.05`,
#' at `c0 = -1e-3` near `x = 1.5`, at `c0 = -0.01` in places, at `c0 = -0.1` only far above the unit root, as at
#' `c = 5`, `x = 14`), and above the unit root
#' at large `x`, a Gauss-Legendre rule of more than 4096 nodes (at `c = 7` from `x` of about 660 to 780). The line
#' integral refuses such a point before evaluating it: `"auto"` takes the conditional integral there, with its checks,
#' and warns where such a value near `c0 = 0` is flagged; `method = "single"`, and case DG, which has the line integral
#' only, raise an error that names the arguments.
#'
#' The computation runs on the engine's worker threads; see [ltratio_threads()].
#'
#' @section Engines:
#' The package holds two engines. The classical engine computes the cases above by the line, conditional and contour
#' integral representations. The general engine computes the same limits for any piecewise smooth deterministic design
#' (polynomial trends of any order, breaks in level and slope at several dates, regressors given by a function or by the
#' regressor matrix of a sample, GLS detrending, the statistic that detrends first; see [dfl_design]) by the inversion
#' formula on paths around the singularities of the moment generating function in its first argument, the inversion on
#' paths around the roots. In place of the case code `cas` a design may be given, with the start by `g0` (a fixed start)
#' or `start = "stationary"`, and GLS detrending by `gls`; a design that is a classical case is that case:
#' `dfl_trend(0)` is case DU, with `start = "stationary"` case SU, and `dfl_trend(1)` with `gls = -13.5` is case DG.
#'
#' `engine = "auto"` takes, for a computation that the classical engine accepts (a case code, or a design that is a
#' classical case within the classical engine's checks of the arguments), the faster of the two engines at the same
#' accuracy, or the more accurate where they differ, by a rule measured on the two engines (for the t-ratio point by
#' point: the general engine for instance with a stationary start and `c0 != 0`, with a zero start and `c0 > 0` near the
#' unit root, and far below it; the values below 1e-3 of the general engine go back to the classical conditional
#' integral, whose error is relative), and the general engine for every other design. `"classical"` forces the classical
#' engine (an error for a design it does not have) and `"general"` the general one (without the classical engine's
#' checks; `method` then names the general engine's own: `"auto"`, `"line"`, the line integral with the explicit kernel,
#' or `"paths"`, the inversion on paths around the roots, see [dist_cdf()]). With `engine = "auto"`, a `method` other than
#' `"auto"` and `verified = TRUE` keep the classical engine for a classical case; for a design the classical engine does
#' not have, `method = "line"` or `"paths"` forces the general engine's line integral or its paths, and for a design that
#' is a classical case `"line"` is the classical line integral and `"paths"` the general engine's paths. `"auto"` is the default, for a case code as for a design: where it takes
#' the general engine for a case code the two engines agree to about 1e-9 (the values move in the tenth digit or below
#' against the earlier versions of the package), and `engine = "classical"` gives the values of the classical engine, those
#' of the earlier versions. The functions of the general engine itself, with all its options, are [dist_cdf()] and the others
#' listed there.
#'
#' @param x Numeric vector or array of points. `-Inf`, `Inf` and `NaN` give 0, 1 (density 0) and `NaN`.
#' @param cas The case, one of `"DK"`, `"DU"`, `"DT"`, `"SK"`, `"SU"`, `"ST"`, `"DG"` (see the section Cases), or a
#'   design of the general engine ([dfl_design]; see the section Engines).
#' @param c The local-to-unity parameter \eqn{c}, \eqn{\phi = 1 + c/n}{phi = 1 + c/n}.
#' @param c0 The local null value \eqn{c_0}{c0}, \eqn{\phi_0 = 1 + c_0/n}{phi0 = 1 + c0/n}.
#' @param g0 The fixed start, \eqn{y_0 = g_0 \sigma \sqrt{n}}{y_0 = g0 sigma sqrt(n)} (cases DK, DU, DT).
#' @param method `"auto"`, `"single"` (or `"line"`), `"double"` (or `"conditional"`) or `"contour"` (a fixed start
#'   in cases DU, DT only); see the section Methods and ranges. With a design, or `engine = "general"`, [ltratio_cdf()],
#'   [ltratio_pdf()] and [ltratio_icdf()] also take `"paths"` (the general engine's inversion on paths; see the section
#'   Engines).
#' @param verified `FALSE` (the default): one evaluation by one method; `TRUE`: the checked computation (see the
#'   section Methods and ranges; with the general engine, a second evaluation by other rules, see [dist_cdf()]).
#' @param engine `"auto"` (the default), `"classical"` or `"general"` (see the section Engines).
#' @param start With a design: the start, `"zero"`, `"fixed"` (with `g0`) or `"stationary"`; `NULL` (the default) is
#'   `"fixed"` where `g0` is not 0 and `"zero"` otherwise. With a case code it must be `NULL` (the code names the
#'   start).
#' @param gls With a design: `NULL` (the default) for least squares detrending, or \eqn{\bar c < 0}{cbar < 0} for GLS
#'   detrending at \eqn{1 + \bar c/n}{1 + cbar/n}. With a case code it must be `NULL`.
#' @param vol `NULL` (the default) or a volatility path ([dfl_volpath()]): volatility and persistence regimes (`c` may
#'   then be a vector of one value per regime), computed by the general engine for the case or the design (`engine =
#'   "classical"` refuses it; `method` `"auto"` or `"paths"`). See [dfl_volpath()].
#' @param rho A normal component (1, the default: none): the limit \eqn{\rho \tau + \sqrt{1 - \rho^2} Z}{rho tau +
#'   sqrt(1 - rho^2) Z} of the covariate-augmented t-ratio ([cadf_c0()]) and of the predictive regression
#'   ([predictive_c0()]), computed by the general engine (`method` `"auto"`, `"lines"` or `"mixture"`). See
#'   [dfl_volpath()].
#' @return A numeric vector with the attributes (`dim`, `dimnames`, `names`) of `x`. A warning is given in a regime
#'   in which one evaluation is known to be brittle, where points that the line integral refuses get flagged values
#'   from the conditional integral (see the section Methods and ranges), and, with `verified = TRUE` and
#'   `method = "double"` or `"contour"`, if the integral reports at some points that its error estimate is not met.
#' @references Bao, Y., Broda, S. A. and Kan, R., "The Exact Limiting Distribution of the Dickey-Fuller t-Ratio
#'   with General Deterministic Terms", working paper, 2026, \url{https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7745921}.
#'
#'   Bao, Y., Broda, S. A. and Kan, R., "Low-Rank Representations of the Classical Dickey-Fuller Limits, with an
#'   Extension to Predictive Regressions", working paper, 2026,
#'   \url{https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7725619}.
#'
#'   Elliott, G., Rothenberg, T. J. and Stock, J. H. (1996), "Efficient Tests for an Autoregressive Unit Root",
#'   Econometrica 64, 813-836.
#' @seealso [ltratio_icdf()], [ltratio_mom()], [coef_cdf()], [adfmax_cdf()], [ltratio_region()], [ltratio_threads()].
#' @examples
#' ltratio_cdf(-2.86157, "DU", -5)            # limiting power of the 5% test of tau_mu at c = -5
#' ltratio_cdf(-2.86157, "SU", -5)            # the same with a stationary start
#' ltratio_cdf(0, "DK", 2, c0 = 2)            # t-ratio centred at the true value
#' ltratio_pdf(c(-3, -2, -1), "DT", -5)       # density
#' ltratio_cdf(-2.86157, "DU", -5, engine = "auto")     # the engine of the faster representation
#' ltratio_cdf(c(-4.5, -4), dfl_perron("C", 0.5))       # a break in level and slope at mid-sample
#' @export
ltratio_cdf <- function(x, cas, c = 0, c0 = 0, g0 = 0, method = "auto", verified = FALSE,
                        engine = c("auto", "classical", "general"), start = NULL, gls = NULL, vol = NULL, rho = 1) {
  dfl_dist(2L, x, cas, c, c0, g0, method, verified, dfl_engine(engine), start, gls, vol, rho)
}

#' @rdname ltratio_cdf
#' @export
ltratio_pdf <- function(x, cas, c = 0, c0 = 0, g0 = 0, method = "auto", verified = FALSE,
                        engine = c("auto", "classical", "general"), start = NULL, gls = NULL, vol = NULL, rho = 1) {
  dfl_dist(1L, x, cas, c, c0, g0, method, verified, dfl_engine(engine), start, gls, vol, rho)
}

dfl_dist <- function(kind, x, cas, c, c0, g0, method, verified, engine = "auto", start = NULL, gls = NULL, vol = NULL, rho = 1) {
  tg <- dfl_target(cas, g0, start, gls)
  e <- dfl_ext(vol, rho, c)
  if (!is.null(e$v)) {
    # (the extensions of the limit: the general engine for the case or the design, through the router; engine =
    # "classical" refuses them)
    pl <- list(spec = dfl_spec(tg, is.null(tg$D)), engine = dfl_engine_code(engine))
    a <- c(e$c, dfl_num(c0, "c0"), dfl_ext_method(method), dfl_flag(verified, "verified"), e$v)
    return(dfl_shape(dfl_route_call(if (kind == 2L) "cdfx" else "pdfx", dfl_real(x, "x"), pl, a), x))
  }
  me <- dfl_method(method)
  ve <- dfl_flag(verified, "verified")
  c <- dfl_num(c, "c")
  c0 <- dfl_num(c0, "c0")
  pl <- dfl_plan(tg, engine, c, c0, me, ve)
  if (pl$path == "router") {
    return(dfl_shape(dfl_route_call(if (kind == 2L) "cdf" else "pdf", dfl_real(x, "x"), pl, c(c, c0, me, ve)), x))
  }
  k <- dfl_case(pl$case)
  r <- .Call(C_dfl_r_dist, kind, dfl_real(x, "x"), k$model, c, c0, as.double(pl$g0), k$stat, me, ve)
  if (ve == 1L && me >= 2L && r[[2]] > 0L) {
    warning(sprintf("the %s integral did not meet its error estimate at %d of %d points",
                    if (me == 2L) "conditional" else "contour", r[[2]], length(x)), call. = FALSE)
  }
  if (nzchar(r[[3]])) warning(r[[3]], call. = FALSE)
  if (nzchar(r[[4]])) warning(r[[4]], call. = FALSE)
  dfl_shape(r[[1]], x)
}

#' Quantiles of the limiting Dickey-Fuller t-ratio
#'
#' The p-quantiles of the limit of [ltratio_cdf()], by Newton's method with [ltratio_pdf()] as the derivative,
#' safeguarded by bisection, started at the normal approximation from the first two moments ([ltratio_mom()]).
#'
#' @inheritParams ltratio_cdf
#' @param p Numeric vector or array of probabilities in (0, 1).
#' @param tol The tolerance on the quantile.
#' @param info If `TRUE`, return a list with the quantiles `q`, the distribution function `F` at `q`, the number
#'   of iterations `iter` and the convergence flags `conv`.
#' @return The quantiles, with the attributes of `p` (or the list described under `info`; where the router computes,
#'   `F` is the distribution function at `q` by the same engine and `iter` and `conv` are `NA`). A warning is given if
#'   Newton's method did not converge at some points.
#' @inheritSection ltratio_cdf Cases
#' @inheritSection ltratio_cdf Engines
#' @seealso [ltratio_cdf()], [dist_quantile()].
#' @examples
#' ltratio_icdf(c(0.01, 0.05, 0.10), "DU")   # asymptotic critical values of tau_mu
#' ltratio_icdf(0.05, dfl_trend(2))          # the quadratic trend (the general engine)
#' @export
ltratio_icdf <- function(p, cas, c = 0, c0 = 0, g0 = 0, method = "auto", tol = 1e-10, info = FALSE, verified = FALSE,
                         engine = c("auto", "classical", "general"), start = NULL, gls = NULL, vol = NULL, rho = 1) {
  engine <- dfl_engine(engine)
  tg <- dfl_target(cas, g0, start, gls)
  e <- dfl_ext(vol, rho, c)
  if (!is.null(e$v)) {
    pv <- dfl_real(p, "p")
    if (!all(!is.na(pv) & pv > 0 & pv < 1)) stop("p must lie in (0, 1)", call. = FALSE)
    pl <- list(spec = dfl_spec(tg, is.null(tg$D)), engine = dfl_engine_code(engine))
    me <- dfl_ext_method(method)
    q <- dfl_route_call("icdfx", pv, pl, c(e$c, dfl_num(c0, "c0"), dfl_num(tol, "tol"), me, e$v))
    if (!isTRUE(info)) return(dfl_shape(q, p))
    F <- dfl_route_call("cdfx", q, pl, c(e$c, dfl_num(c0, "c0"), me, 0, e$v))
    na <- rep(NA_integer_, length(q))
    return(lapply(list(q = q, F = F, iter = na, conv = as.logical(na)), dfl_shape, x = p))
  }
  me <- dfl_method(method)
  ve <- dfl_flag(verified, "verified")
  c <- dfl_num(c, "c")
  c0 <- dfl_num(c0, "c0")
  tol <- dfl_num(tol, "tol")
  pl <- dfl_plan(tg, engine, c, c0, me, ve)
  if (pl$path == "router") {
    pv <- dfl_real(p, "p")
    if (!all(!is.na(pv) & pv > 0 & pv < 1)) stop("p must lie in (0, 1)", call. = FALSE)
    q <- dfl_route_call("icdf", pv, pl, c(c, c0, tol, me))
    if (!isTRUE(info)) return(dfl_shape(q, p))
    F <- dfl_route_call("cdf", q, pl, c(c, c0, me, 0))
    na <- rep(NA_integer_, length(q))
    return(lapply(list(q = q, F = F, iter = na, conv = as.logical(na)), dfl_shape, x = p))
  }
  k <- dfl_case(pl$case)
  r <- .Call(C_dfl_r_icdf, dfl_real(p, "p"), k$model, c, c0, as.double(pl$g0),
             k$stat, me, tol, ve)
  if (nzchar(r[[5]])) warning(r[[5]], call. = FALSE)
  if (nzchar(r[[6]])) warning(r[[6]], call. = FALSE)
  r <- r[1:4]
  names(r) <- c("q", "F", "iter", "conv")
  if (!all(r$conv)) {
    warning(sprintf("Newton's method did not converge at %d of %d points", sum(!r$conv), length(r$conv)),
            call. = FALSE)
  }
  r <- lapply(r, dfl_shape, x = p)
  if (isTRUE(info)) r else r$q
}

#' Moments of the limiting Dickey-Fuller t-ratio
#'
#' The moments \eqn{E[\tau^j]}{E[tau^j]}, \eqn{j = 1, \dots, p}{j = 1..p}, of the limit
#' \eqn{\tau = (U - c_0 V)/\sqrt{V}}{tau = (U - c0 V)/sqrt(V)} of [ltratio_cdf()], as one-dimensional integrals of
#' the joint moment generating function of \eqn{(U, V)} (the moments of the limit in the Supplementary Material of the
#' paper): \eqn{E\tau^j = \Gamma(j/2)^{-1} \int_0^\infty s^{j/2-1} m_j(s)\,ds}{E tau^j = Gamma(j/2)^(-1) int_0^Inf
#' s^(j/2-1) m_j(s) ds}, \eqn{m_j(s)}{m_j(s)} the j-th derivative of the moment generating function along its affine path.
#' Not available for case DG in the classical engine (a design of case DG takes the general engine).
#'
#' By default (`method = "auto"` with `engine = "auto"`, the cases DK, DU, DT with a zero or fixed start and SK, SU, ST,
#' up to `c = 20` and 31 moments) the derivative comes from one Taylor-mode jet of order j at each node of an adaptive
#' quadrature in \eqn{\log s}{log s}, the method of Yong Bao and Raymond Kan, "Exact Finite-Sample Distribution and
#' Moments of the Least Squares Estimator in First-Order Autoregressions: O(1) Algorithms", with the closed-form
#' transform of the six cases; `tol` is the relative tolerance of its quadrature, and a warning names the moments whose
#' error estimate exceeds `tol` times their size (`verified = TRUE` recomputes them by the engines). `method = "jets"`
#' forces it (also above `c = 20`); `method = "engine"` is the computation of before (the derivatives as Cauchy
#' integrals on a circle, by the classical engine or the general engine), which reproduces those values; `engine =
#' "classical"` or `"general"` takes that computation too.
#'
#' @inheritParams ltratio_cdf
#' @param p The number of moments, a positive integer.
#' @param info If `TRUE`, return a list with the moments `m` and the error estimates `err` of the integrals (`NA`
#'   where the router computes by the engines).
#' @param method `"auto"` (the default), `"jets"` (the moments from jets of Bao and Kan) or `"engine"` (the engines'
#'   circles).
#' @param tol The relative tolerance of the quadrature of the moments from jets (1e-10).
#' @return A numeric vector of length `p` (or the list described under `info`).
#' @inheritSection ltratio_cdf Cases
#' @inheritSection ltratio_cdf Engines
#' @seealso [ltratio_cdf()], [dist_moments()].
#' @examples
#' ltratio_mom(4, "DT")      # the first four moments of tau_tau
#' ltratio_mom(2, "DU")[2]   # 3.0559953450 (Nabeya, 1999)
#' ltratio_mom(2, dfl_perron("A", 0.5))
#' @export
ltratio_mom <- function(p, cas, c = 0, c0 = 0, g0 = 0, info = FALSE, engine = c("auto", "classical", "general"),
                        start = NULL, gls = NULL, method = "auto", tol = 1e-10, verified = FALSE) {
  if (!is.numeric(p) || length(p) != 1L || is.na(p) || p != round(p) || p < 1) {
    stop("p must be a positive integer", call. = FALSE)
  }
  engine <- dfl_engine(engine)
  tg <- dfl_target(cas, g0, start, gls)
  c <- dfl_num(c, "c")
  c0 <- dfl_num(c0, "c0")
  im <- dfl_imethod(method, "ltratio_mom")
  ve <- dfl_flag(verified, "verified")
  pl <- dfl_plan(tg, engine, c, c0)
  # (the moments from jets with engine = "classical": the router, which takes them for method = "jets")
  if (pl$path == "legacy" && im == 5L) pl <- list(path = "router", spec = dfl_spec(tg, TRUE), engine = dfl_engine_code(engine))
  # (a design of case DG, which has no moments in the classical engine: the general engine)
  if (pl$path == "router" && identical(tg$case, "DG") && !is.null(tg$D)) pl <- list(path = "router", spec = dfl_spec(tg, FALSE), engine = 2L)
  if (pl$path == "router") {
    if (im == 5L) pl$engine <- dfl_engine_code(engine)
    r <- dfl_route_call2("mom", 0, pl, c(c, c0, p, im, dfl_tol(tol), ve))
    return(if (isTRUE(info)) list(m = r$v, err = ifelse(is.nan(r$e), NA_real_, r$e)) else r$v)
  }
  k <- dfl_case(pl$case)
  r <- .Call(C_dfl_r_moments, as.integer(p), k$model, c, c0, as.double(pl$g0), k$stat)
  names(r) <- c("m", "err")
  if (isTRUE(info)) r else r$m
}
