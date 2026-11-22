#' The point-optimal statistic and the power envelope for any design
#'
#' The limit of the point-optimal invariant statistic \eqn{P_T}{P_T} of Elliott, Rothenberg and Stock (1996) against
#' \eqn{\bar c < 0}{cbar < 0}, for any design ([dfl_design]), with a zero start: \eqn{P_T}{P_T} rejects the unit root
#' for small values, and its limit is a quadratic form \eqn{Q} in the Ornstein-Uhlenbeck process (with the data
#' GLS-detrended at \eqn{\bar c}{cbar}), whose law is computed by one inversion of Gil-Pelaez (the trapezoidal rule in
#' \eqn{u = \sinh v}{u = sinh(v)} with `N` steps up to `umax`). `pt_cdf()` and `pt_pdf()` give its distribution
#' function and density at `q` under \eqn{\phi = 1 + c/n}{phi = 1 + c/n}, `pt_quantile()` its quantiles at the unit
#' root, `pt_power()` the power at `c` of the level-`level` test against `cbar`, and `pt_envelope()` the Gaussian
#' power envelope of the design: the power at `c` of the point-optimal test against `c` itself (`level` for
#' `c >= 0`, `NaN` at `NaN`). For the intercept and the linear trend it is the envelope of Elliott, Rothenberg and
#' Stock ([power_envelope()]).
#'
#' The inversion is refused with an error where its grid does not resolve the transform: far above the unit root, where
#' the singularity of the transform of \eqn{Q} nearest to the real line comes closer than the step of the grid allows
#' (with the default grid from about `c = 3.7` with a linear trend at `cbar = -13.5`, `c = 3.1` with an intercept at
#' `cbar = -7`; the message names the `N` that would do), and at \eqn{|q|} beyond \eqn{0.15 N/\mathrm{asinh}(u_{max})}{0.15
#' N/asinh(umax)} (1000 with the default grid). `-Inf`, `Inf` and `NaN` give 0, 1 (density 0) and `NaN`; `p` and `level`
#' must lie in (0, 1).
#'
#' @param q Numeric vector or array of points of the statistic.
#' @param p Numeric vector or array of probabilities in (0, 1).
#' @param D A design ([dfl_design]).
#' @param cbar The alternative \eqn{\bar c < 0}{cbar < 0} of the point-optimal test.
#' @param c `pt_cdf()`, `pt_pdf()`: the local-to-unity parameter. `pt_power()`, `pt_envelope()`: numeric vector or
#'   array of them.
#' @param level The level of the test.
#' @param N,umax The number of steps and the end of the grid of the inversion formula of Gil-Pelaez.
#' @return A numeric vector with the attributes of `q`, `p` or `c`; that of `pt_power()` has the attribute
#'   `critical`, the critical value of the test.
#' @references Elliott, G., Rothenberg, T. J. and Stock, J. H. (1996), "Efficient Tests for an Autoregressive Unit
#'   Root", Econometrica 64, 813-836.
#' @seealso [power_envelope()], [lm_cdf()], [mza_cdf()], [dist_cdf()].
#' @examples
#' D <- dfl_trend(1)
#' pt_quantile(0.05, D, -13.5)          # 5.7045: the 5% critical value of the limit of P_T, trend case
#' pt_cdf(c(4, 7), dfl_perron("C", 0.5), -10, -10)
#' pt_envelope(-10, dfl_perron("A", 0.5))
#' @export
pt_cdf <- function(q, D, cbar, c = 0, N = 60000, umax = 4000) {
  dfl_det_finite("ptcdf", D, q, c(dfl_cbar(cbar), dfl_num(c, "c"), dfl_int(N, "N"), dfl_num(umax, "umax")))
}

#' @rdname pt_cdf
#' @export
pt_pdf <- function(q, D, cbar, c = 0, N = 60000, umax = 4000) {
  dfl_det_finite("ptpdf", D, q, c(dfl_cbar(cbar), dfl_num(c, "c"), dfl_int(N, "N"), dfl_num(umax, "umax")), hi = 0)
}

#' @rdname pt_cdf
#' @export
pt_quantile <- function(p, D, cbar, N = 60000, umax = 4000) {
  pv <- dfl_prob(p)
  dfl_shape(dfl_det_call("ptquantile", D, pv, c(dfl_cbar(cbar), dfl_int(N, "N"), dfl_num(umax, "umax")))[[1]], p)
}

#' @rdname pt_cdf
#' @export
pt_power <- function(c, D, cbar, level = 0.05, N = 60000, umax = 4000) {
  cv <- dfl_finite(c, "c")
  r <- dfl_det_call("ptpower", D, cv, c(dfl_cbar(cbar), dfl_num(level, "level"), dfl_int(N, "N"), dfl_num(umax, "umax")))
  v <- dfl_shape(r[[1]], c)
  attr(v, "critical") <- if (length(r[[2]])) r[[2]][1] else NA_real_
  v
}

#' @rdname pt_cdf
#' @export
pt_envelope <- function(c, D, level = 0.05, N = 60000, umax = 4000) {
  cv <- dfl_real(c, "c")
  out <- rep(NaN, length(cv))
  ok <- !is.na(cv)
  if (!all(is.finite(cv[ok]))) stop("c must be finite (or NaN)", call. = FALSE)
  if (any(ok)) out[ok] <- dfl_det_call("ptenvelope", D, cv[ok], c(dfl_num(level, "level"), dfl_int(N, "N"), dfl_num(umax, "umax")))[[1]]
  dfl_shape(out, c)
}

#' The LM statistic for any design
#'
#' The limit of the LM (score) statistic of Schmidt and Phillips (1992) and Amsler and Lee (1995), with the breaks of
#' the design, under \eqn{\phi = 1 + c/n}{phi = 1 + c/n}, zero start: its limit is \eqn{-1/(2\sqrt{V})}{-1/(2
#' sqrt(V))} with \eqn{V} the integral of the square of the detrended partial sums, whose law is one inversion of the
#' transform of \eqn{V} on a path of its own. The design ([dfl_design]) must contain a linear trend. `lm_cdf()` gives the
#' distribution function, `lm_quantile()` the quantiles and `lm_power()` the power at `c` of the level-`level` test (its
#' critical value at the unit root in the attribute `critical`). `c` must be at most 10 (above, the inversion is not
#' accurate: off by 2e-7 at `c = 12` and by 3e-5 at `c = 15`), and an error is raised beyond.
#'
#' @param x Numeric vector or array of points. `-Inf`, `Inf` and `NaN` give 0, 1 and `NaN`.
#' @param p Numeric vector or array of probabilities in (0, 1).
#' @param D A design ([dfl_design]) with a linear trend among its regressors.
#' @param c `lm_cdf()`, `lm_quantile()`: the local-to-unity parameter. `lm_power()`: numeric vector or array of them.
#' @param level The level of the test.
#' @return A numeric vector with the attributes of `x`, `p` or `c`.
#' @references Schmidt, P. and Phillips, P. C. B. (1992), "LM Tests for a Unit Root in the Presence of Deterministic
#'   Trends", Oxford Bulletin of Economics and Statistics 54, 257-287.
#'
#'   Amsler, C. and Lee, J. (1995), "An LM Test for a Unit Root in the Presence of a Structural Change", Econometric
#'   Theory 11, 359-368.
#' @seealso [pt_cdf()], [dist_cdf()].
#' @examples
#' lm_quantile(0.05, dfl_trend(1))       # -3.0230: the 5% critical value, linear trend
#' lm_power(c(-5, -10), dfl_perron("B", 0.5))
#' @export
lm_cdf <- function(x, D, c = 0) {
  dfl_det_finite("lmcdf", D, x, dfl_num(c, "c"))
}

#' @rdname lm_cdf
#' @export
lm_quantile <- function(p, D, c = 0) {
  pv <- dfl_prob(p)
  dfl_shape(dfl_det_call("lmquantile", D, pv, dfl_num(c, "c"))[[1]], p)
}

#' @rdname lm_cdf
#' @export
lm_power <- function(c, D, level = 0.05) {
  cv <- dfl_finite(c, "c")
  r <- dfl_det_call("lmpower", D, cv, dfl_num(level, "level"))
  v <- dfl_shape(r[[1]], c)
  attr(v, "critical") <- r[[2]][1]
  v
}

#' MZ-alpha and MSB for any design
#'
#' The companions of the GLS-detrended t-ratio among the modified tests of Ng and Perron (2001), for any design
#' ([dfl_design]) with the data GLS-detrended at \eqn{1 + \bar c/n}{1 + cbar/n} (`gls`, by default -13.5), zero start,
#' under \eqn{\phi = 1 + c/n}{phi = 1 + c/n}: \eqn{MZ_\alpha}{MZ_alpha}, whose limit is that of the coefficient
#' statistic with GLS detrending ([coef_cdf()]), and MSB, whose limit is \eqn{\sqrt{V}}{sqrt(V)} with GLS detrending.
#' `mza_cdf()` and `msb_cdf()` give the distribution functions, `mza_quantile()` (Newton's method from `y0`) and
#' `msb_quantile()` the quantiles.
#'
#' @param y,x Numeric vectors or arrays of points. `-Inf`, `Inf` and `NaN` give 0, 1 and `NaN`.
#' @param p Numeric vector or array of probabilities in (0, 1).
#' @param D A design ([dfl_design]).
#' @param c The local-to-unity parameter.
#' @param gls The GLS detrending parameter \eqn{\bar c < 0}{cbar < 0}.
#' @param y0 The starting point of the search.
#' @param tol The tolerance of the search.
#' @return A numeric vector with the attributes of `y`, `x` or `p`.
#' @references Ng, S. and Perron, P. (2001), "Lag Length Selection and the Construction of Unit Root Tests with Good
#'   Size and Power", Econometrica 69, 1519-1554.
#' @seealso [pt_cdf()], [coef_cdf()], [dist_cdf()].
#' @examples
#' D <- dfl_trend(1)
#' mza_quantile(0.05, D)                 # -16.591: the 5% critical value of MZ_alpha, trend case
#' msb_cdf(c(0.15, 0.25), dfl_perron("A", 0.5))
#' @export
mza_cdf <- function(y, D, c = 0, gls = -13.5) {
  dfl_det_finite("mzacdf", D, y, c(dfl_num(c, "c"), dfl_cbar(gls, "gls")))
}

#' @rdname mza_cdf
#' @export
mza_quantile <- function(p, D, c = 0, gls = -13.5, y0 = -20, tol = 1e-9) {
  pv <- dfl_prob(p)
  dfl_shape(dfl_det_call("mzaquantile", D, pv, c(dfl_num(c, "c"), dfl_cbar(gls, "gls"), dfl_num(y0, "y0"),
                                                  dfl_num(tol, "tol")))[[1]], p)
}

#' @rdname mza_cdf
#' @export
msb_cdf <- function(x, D, c = 0, gls = -13.5) {
  dfl_det_finite("msbcdf", D, x, c(dfl_num(c, "c"), dfl_cbar(gls, "gls")))
}

#' @rdname mza_cdf
#' @export
msb_quantile <- function(p, D, c = 0, gls = -13.5) {
  pv <- dfl_prob(p)
  dfl_shape(dfl_det_call("msbquantile", D, pv, c(dfl_num(c, "c"), dfl_cbar(gls, "gls")))[[1]], p)
}

# a negative detrending parameter or alternative
dfl_cbar <- function(v, name = "cbar") {
  v <- dfl_num(v, name)
  if (!(v < 0)) stop(name, " must be negative", call. = FALSE)
  v
}
# finite numbers, as a double vector
dfl_finite <- function(v, name) {
  v <- dfl_real(v, name)
  if (!all(is.finite(v))) stop(name, " must be finite", call. = FALSE)
  v
}
