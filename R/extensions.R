#' The extensions of the limit: volatility regimes, a normal component, the first seasonal t-ratio
#'
#' Three extensions of the limit of the Dickey-Fuller t-ratio, computed by the general engine (the functions of
#' [dist_cdf()] and [ltratio_cdf()] take the arguments `vol` and `rho`):
#'
#' \describe{
#'   \item{Volatility and persistence regimes}{`dfl_volpath(breaks, omega)`: the standard deviation of the innovations
#'     piecewise constant between known dates (the fractions `breaks`, strictly increasing in (0, 1)), `omega` on the
#'     `length(breaks) + 1` regimes (any positive scale: only the ratios matter). With `vol = dfl_volpath(...)` the
#'     t-ratio's limit is that of the process \eqn{dJ = c(r) J\,dr + \omega(r)\,dW}{dJ = c(r) J dr + omega(r) dW}, and
#'     `c` may be a vector of one value per regime (changing persistence). Every `c` in \[-480, 0\]; GLS detrending
#'     (`gls`), the statistic that detrends first ([dfl_detrend_first()]) and the stationary start as without `vol`,
#'     the stationary start being that of the first regime, \eqn{J(0) \sim N(0, \omega_1^2/(-2 c_1))}{J(0) ~ N(0,
#'     omega_1^2/(-2 c_1))} (in the scale of `omega`), which needs `c < 0` on the first regime; GLS detrending and the
#'     statistic that detrends first take the zero start only, as without `vol`. A fixed start (`gamma0 != 0`) is
#'     refused with more than one regime. `method` `"auto"` or `"paths"` (the inversion on paths of the volatility
#'     transform). One regime with a scalar `c` is the engine itself.}
#'   \item{A normal component}{`rho`: the limit \eqn{\rho \tau_{c_0} + \sqrt{1 - \rho^2} Z}{rho tau_c0 + sqrt(1 - rho^2) Z}
#'     with \eqn{Z} standard normal and independent of the t-ratio's limit \eqn{\tau_{c_0}}{tau_c0}: the limit of the
#'     covariate-augmented t-ratio of Hansen (1995) at `c0 = cadf_c0(c, rho, R)` and of the t-ratio of the slope in a
#'     predictive regression with a regressor local to unity at `c0 = predictive_c0(c, rho, b)`. `rho` in \[-1, 1\]
#'     (1, the default, is the t-ratio itself, 0 the standard normal law, a negative `rho` by the reflection); `method`
#'     `"auto"`, `"lines"` (the inversion on two lines with the normal component in the kernel: a zero or stationary
#'     start, `c <= 0`) or `"mixture"` (the normal average of the t-ratio's distribution function); `"auto"` takes the
#'     lines for \eqn{0.1 \le \rho^2 \le 0.99}{0.1 <= rho^2 <= 0.99} with a zero or stationary start at `c <= 0` and the
#'     mixture otherwise. Least squares detrending and the one-step regression; not together with `vol`.}
#'   \item{The first seasonal t-ratio}{`hegy_t3_cdf()`, `hegy_t3_pdf()`, `hegy_t3_quantile()`: the t-ratio
#'     \eqn{t_3 = (R - c_0 V)/\sqrt{V}}{t3 = (R - c0 V)/sqrt(V)} of the regression of Hylleberg, Engle, Granger and Yoo
#'     at a harmonic frequency (\eqn{\pi/2}{pi/2} for quarterly data, the coefficient of \eqn{y_{3,t-2}}{y_{3,t-2}}; each
#'     pair of harmonic frequencies of monthly data), the one of Dickey-Fuller shape: \eqn{R = Re \int \tilde Z^* dZ}{R =
#'     Re int Z~* dZ}, \eqn{V = \int |\tilde Z|^2}{V = int |Z~|^2} of the complex Ornstein-Uhlenbeck process and its
#'     residual on the deterministic terms of the design `D` (`dfl_trend(-1)` no seasonal terms, `dfl_trend(0)` seasonal
#'     dummies, `dfl_trend(1)` seasonal dummies and seasonal trends, any design of the general engine), by the sum of
#'     the residues at the roots of the engine's determinant along its outer contour (where two roots of a group come
#'     within 1e-4 of each other, the inversion on paths of the doubled transform); zero start, `c <= 0` (in the units
#'     of the roots' modulus \eqn{1 + c/n}{1 + c/n}, \eqn{n} the number of years), with `c0` and `vol`. `verified =
#'     TRUE` checks the residues against the doubled transform.}
#' }
#'
#' @param breaks The interior regime boundaries, strictly increasing in (0, 1) (`numeric(0)`: one regime).
#' @param omega The standard deviations of the innovations on the regimes, finite and positive (`NULL`: all equal).
#' @param c The local-to-unity parameter; with `vol` a number or a vector of one value per regime.
#' @param rho The correlation parameter: Hansen's long-run correlation (`cadf_c0()`), the correlation of the innovations
#'   of the predictive regression (`predictive_c0()`); in \[-1, 1\] and not 0.
#' @param R The ratio of the long-run standard deviations of the series and the error of the augmented regression
#'   (Hansen's; `abs(rho)` for the well-specified regression).
#' @param b The local alternative of the slope, \eqn{\beta = b \sigma_u/(\sigma_v n)}{beta = b sigma_u/(sigma_v n)}.
#' @param x Numeric vector or array of points (-Inf, Inf and NaN give 0, 1 (density 0) and NaN).
#' @param p Numeric vector or array of probabilities in (0, 1).
#' @param D A design ([dfl_design]).
#' @param c0 The local null value (the units of `c`).
#' @param vol `NULL` or a volatility path (`dfl_volpath()`).
#' @param verified `FALSE`: one evaluation; `TRUE`: the residues checked against the inversion on paths of the doubled
#'   transform (the latter returned, with a warning, where they differ by more than 1e-9).
#' @param x0 The starting point of the quantile search.
#' @param tol The tolerance of the quantile search.
#' @return `dfl_volpath()`: a volatility path (class `"dfl_volpath"`: the boundaries `lam` and the scaled standard
#'   deviations `om`). `cadf_c0()`, `predictive_c0()`: the local null value. `hegy_t3_cdf()`, `hegy_t3_pdf()`,
#'   `hegy_t3_quantile()`: numeric vectors with the attributes of `x` or `p`.
#' @references Hansen, B. E. (1995), "Rethinking the Univariate Approach to Unit Root Testing: Using Covariates to
#'   Increase Power", Econometric Theory 11, 1148-1171.
#'
#'   Hylleberg, S., Engle, R. F., Granger, C. W. J. and Yoo, B. S. (1990), "Seasonal Integration and Cointegration",
#'   Journal of Econometrics 44, 215-238.
#' @seealso [dist_cdf()], [ltratio_cdf()].
#' @examples
#' V <- dfl_volpath(0.3, c(1, 3))         # the standard deviation triples at 30%
#' dist_cdf(-2.86, dfl_trend(0), 0, vol = V)
#' dist_cdf(-2.86, dfl_trend(0), c(0, -20), vol = dfl_volpath(0.5))   # changing persistence
#' dist_quantile(0.05, dfl_trend(1), 0, gls = -13.5, vol = V)   # DF-GLS: -2.9967, not -2.8456
#' # the stationary start of the first regime (c = -10 up to mid-sample, -5 after it)
#' dist_cdf(-2.86, dfl_trend(0), c(-10, -5), start = "stationary", vol = dfl_volpath(0.5))
#' cadf_c0(-5, sqrt(0.5))                 # 5
#' ltratio_cdf(1.645, "DU", 0, rho = -0.95, c0 = predictive_c0(0, -0.95))   # 0.5788
#' hegy_t3_cdf(-3.3442, dfl_trend(0))     # 0.05, seasonal dummies
#' hegy_t3_quantile(0.05, dfl_trend(0))
#' @name dfl_volpath
#' @export
dfl_volpath <- function(breaks = numeric(0), omega = NULL) {
  if (!is.numeric(breaks) || anyNA(breaks) || !all(is.finite(breaks) & breaks > 0 & breaks < 1)) {
    stop("volpath: the break dates must lie in (0, 1)", call. = FALSE)
  }
  br <- as.double(breaks)
  if (length(br) > 1L && any(diff(br) <= 0)) stop("volpath: the break dates must be strictly increasing", call. = FALSE)
  if (is.null(omega)) omega <- rep(1, length(br) + 1L)
  if (length(omega) != length(br) + 1L) {
    stop(sprintf("volpath: one standard deviation omega per regime: %d for %d break dates, %d given",
                 length(br) + 1L, length(br), length(omega)), call. = FALSE)
  }
  if (!is.numeric(omega) || anyNA(omega) || !all(is.finite(omega) & omega > 0)) {
    stop("volpath: the standard deviations omega must be finite and positive", call. = FALSE)
  }
  om <- as.double(omega)
  lam <- c(0, br, 1)
  sc <- om / sqrt(sum(om^2 * diff(lam)))
  structure(list(lam = lam, om = sc, breaks = br, omega = om), class = "dfl_volpath")
}

#' @export
print.dfl_volpath <- function(x, ...) {
  cat("volatility path: ", length(x$om), " regime", if (length(x$om) > 1L) "s", " (breaks ",
      if (length(x$breaks)) paste(format(x$breaks), collapse = ", ") else "none", "; omega ",
      paste(format(signif(x$om, 6)), collapse = ", "), ")\n", sep = "")
  invisible(x)
}

#' @rdname dfl_volpath
#' @export
cadf_c0 <- function(c, rho, R = abs(rho)) {
  if (!is.numeric(rho) || length(rho) != 1L || !is.finite(rho) || rho < -1 || rho > 1 || rho == 0) {
    stop("cadf_c0: rho must lie in [-1, 1] and not be 0", call. = FALSE)
  }
  if (!is.numeric(R) || length(R) != 1L || !is.finite(R) || R <= 0) stop("cadf_c0: R must be finite and positive", call. = FALSE)
  if (!is.numeric(c) || length(c) != 1L || !is.finite(c)) stop("cadf_c0: c must be finite", call. = FALSE)
  if (c == 0) 0 else as.double(c) * (1 - 1 / (as.double(rho) * as.double(R)))
}

#' @rdname dfl_volpath
#' @export
predictive_c0 <- function(c, rho, b = 0) {
  if (!is.numeric(rho) || length(rho) != 1L || !is.finite(rho) || rho < -1 || rho > 1 || rho == 0) {
    stop("predictive_c0: rho must lie in [-1, 1] and not be 0", call. = FALSE)
  }
  if (!is.numeric(c) || length(c) != 1L || !is.finite(c) || !is.numeric(b) || length(b) != 1L || !is.finite(b)) {
    stop("predictive_c0: c and b must be finite", call. = FALSE)
  }
  as.double(c) - as.double(b) / as.double(rho)
}

#' @rdname dfl_volpath
#' @export
hegy_t3_cdf <- function(x, D, c = 0, c0 = 0, vol = NULL, verified = FALSE) {
  e <- dfl_ext(vol, 1, c)
  dfl_det_finite("t3cdf", D, x, e$c, dfl_opts(c0, NULL, 0, NULL, verified, ext = e$v))
}

#' @rdname dfl_volpath
#' @export
hegy_t3_pdf <- function(x, D, c = 0, c0 = 0, vol = NULL, verified = FALSE) {
  e <- dfl_ext(vol, 1, c)
  dfl_det_finite("t3pdf", D, x, e$c, dfl_opts(c0, NULL, 0, NULL, verified, ext = e$v), hi = 0)
}

#' @rdname dfl_volpath
#' @export
hegy_t3_quantile <- function(p, D, c = 0, c0 = 0, vol = NULL, x0 = -2, tol = 1e-10) {
  pv <- dfl_prob(p)
  e <- dfl_ext(vol, 1, c)
  dfl_shape(dfl_det_call("t3quantile", D, pv, c(e$c, dfl_num(x0, "x0"), dfl_num(tol, "tol")),
                         dfl_opts(c0, NULL, 0, NULL, ext = e$v))[[1]], p)
}

# the extensions of a call for the engine: list(c = the scalar c of the call, v = c(rho, nreg, breaks, omega, c of the
# regimes) or NULL without them); c a vector needs vol (one value per regime)
dfl_ext <- function(vol, rho, c) {
  rho <- dfl_num(rho, "rho")
  if (is.numeric(c) && length(c) > 1L && is.null(vol)) {
    stop("c is a vector (one value per regime): that needs a volatility path, vol = dfl_volpath(breaks, omega)", call. = FALSE)
  }
  if (length(c) == 1L || is.null(vol)) c <- dfl_num(c, "c")
  else if (!is.numeric(c) || anyNA(c)) stop("c must be a number or a vector of numbers (one value per regime)", call. = FALSE)
  if (is.null(vol)) return(list(c = as.double(c), v = if (rho == 1) NULL else c(rho, 0)))
  if (!inherits(vol, "dfl_volpath")) stop("vol must be a volatility path, dfl_volpath(breaks, omega)", call. = FALSE)
  k <- length(vol$om)
  if (length(c) != 1L && length(c) != k) {
    stop(sprintf("c must be a number or a vector with one value per regime of the volatility path (%d), not %d", k, length(c)), call. = FALSE)
  }
  list(c = as.double(c[1L]), v = c(rho, k, vol$breaks, vol$omega, if (length(c) > 1L) as.double(c)))
}
# the router's method with the extensions: 0 "auto", 1 "single" ("line"), 2 "double" ("conditional"), 3 "contour",
# 4 "paths", 8 "lines", 9 "mixture"
dfl_ext_method <- function(method) {
  msg <- "with vol or rho the method must be \"auto\", \"paths\" (vol), \"lines\" or \"mixture\" (rho)"
  if (!is.character(method) || length(method) != 1L) stop(msg, call. = FALSE)
  m <- match(tolower(method), c("auto", "single", "line", "double", "conditional", "contour", "paths", "lines", "mixture"))
  if (is.na(m)) stop(msg, call. = FALSE)
  c(0L, 1L, 1L, 2L, 2L, 3L, 4L, 8L, 9L)[m]
}
