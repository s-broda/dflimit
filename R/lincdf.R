#' Limiting distribution of the linear form U - aV
#'
#' The limit of \eqn{P(U - aV \le b)}{P(U - aV <= b)} under \eqn{\phi = 1 + c/n}{phi = 1 + c/n}, zero start, no
#' drift, with \eqn{U} and \eqn{V} as in [ltratio_cdf()]: case `"DK"` (no deterministic terms, and the constant case
#' with GLS demeaning, whose limit is the same) or `"DG"` (the GLS-detrended trend case of Elliott, Rothenberg and
#' Stock (1996), detrended at `cbar`). With \eqn{a = \bar c/2}{a = cbar/2}, \eqn{\bar c < 0}{cbar < 0}, the
#' point-optimal test of Elliott, Rothenberg and Stock against \eqn{\bar c}{cbar} rejects for \eqn{U - aV < b} (its
#' log likelihood ratio is \eqn{\bar c U - \bar c^2 V/2}{cbar U - cbar^2 V/2}); see [power_envelope()]. `a = NaN`
#' gives \eqn{P(V \le b)}{P(V <= b)} (case DK).
#'
#' The probability is computed by the inversion formula of Gil-Pelaez,
#' \eqn{P(X < b) = 1/2 - (1/\pi) \int_0^\infty \mathrm{Im}[\varphi(u) e^{-iub}]/u \, du}{P(X < b) = 1/2 - (1/pi)
#' int_0^Inf Im[phi(u) exp(-iub)]/u du}, from the joint moment generating function of \eqn{(U, V)}, by the
#' trapezoidal rule in \eqn{u = \sinh v}{u = sinh(v)} with 60000 steps up to \eqn{u = 4000}, as `ltratio_lincdf` of
#' the MATLAB and Julia programs.
#'
#' @param b Numeric vector or array of points. `-Inf`, `Inf` and `NaN` give 0, 1 and `NaN`.
#' @param a The coefficient of \eqn{V}, a single number; `NaN` for the distribution of \eqn{V} (case DK).
#' @param cas The case, `"DK"` or `"DG"`.
#' @param c The local-to-unity parameter, finite.
#' @param cbar The detrending parameter \eqn{\bar c}{cbar} of case DG (-13.5: the DF-GLS tests); not used in case DK.
#' @return The probabilities, with the attributes of `b`.
#' @references Elliott, G., Rothenberg, T. J. and Stock, J. H. (1996), "Efficient Tests for an Autoregressive Unit
#'   Root", Econometrica 64, 813-836.
#' @seealso [power_envelope()], [ltratio_region()], [ltratio_cdf()].
#' @examples
#' ltratio_lincdf(-0.2639, -3.5, "DK", 0)    # 0.05: the size of the point-optimal test at cbar = -7
#' ltratio_lincdf(-0.2639, -3.5, "DK", -7)   # its power at c = -7
#' ltratio_lincdf(c(0.1, 0.5, 2), NaN, "DK") # P(V <= b) at the unit root
#' @export
ltratio_lincdf <- function(b, a, cas, c = 0, cbar = -13.5) {
  k <- dfl_case(cas, c("DK", "DG"))
  if (!is.numeric(a) || length(a) != 1L || (is.na(a) && !is.nan(a))) {
    stop("a must be a single number (NaN: the distribution of V)", call. = FALSE)
  }
  P <- .Call(C_dfl_r_lincdf, dfl_real(b, "b"), as.double(a), k$model, dfl_num(c, "c"), dfl_num(cbar, "cbar"))
  dfl_shape(P, b)
}

#' Asymptotic power envelope of Elliott, Rothenberg and Stock (1996)
#'
#' The power at `c` of the point-optimal level-`alpha` test against `c` itself, the envelope of the powers of all
#' tests of the unit root, zero start, no drift: case `"DK"` (no deterministic terms; also the constant case,
#' `"DU"`, where the point-optimal test uses the GLS-demeaned data and the envelope is the same curve) and case
#' `"DG"` (the trend case, the point-optimal test against `c` with the data GLS-detrended at `c`).
#'
#' The test against \eqn{\bar c}{cbar} rejects for small \eqn{U - aV}, \eqn{a = \bar c/2}{a = cbar/2}
#' ([ltratio_lincdf()]): its critical value is the `alpha`-quantile at `c = 0`, found by bisection on
#' \eqn{[-50, 50]} to the width 1e-9, and the envelope is its power at \eqn{\bar c = c}{cbar = c}; `alpha` for
#' `c >= 0`. As `power_envelope` of the MATLAB and Julia programs.
#'
#' @param c Numeric vector or array of the local-to-unity parameters. `NaN` gives `NaN`.
#' @param cas The case, `"DK"`, `"DU"` (the same as `"DK"`) or `"DG"`.
#' @param alpha The level of the tests.
#' @return The envelope, with the attributes of `c`.
#' @references Elliott, G., Rothenberg, T. J. and Stock, J. H. (1996), "Efficient Tests for an Autoregressive Unit
#'   Root", Econometrica 64, 813-836.
#' @seealso [ltratio_lincdf()], [ltratio_cdf()].
#' @examples
#' power_envelope(-7, "DK")                   # 0.4994: at c = -7 the envelope has power one half
#' power_envelope(c(-5, -13.5), "DG")
#' @export
power_envelope <- function(c, cas = "DK", alpha = 0.05) {
  k <- dfl_case(cas, c("DK", "DU", "DG"))
  model <- if (k$model == 4L) 4L else 1L
  P <- .Call(C_dfl_r_envelope, dfl_real(c, "c"), model, dfl_num(alpha, "alpha"))
  dfl_shape(P, c)
}
