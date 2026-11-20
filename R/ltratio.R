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
#' Three exact representations are implemented: the line integral of Theorem S1 of the Supplementary Material of the
#' paper (`method = "single"`, or `"line"`), the conditional integral of Theorem S2 of the Supplementary Material
#' (`"double"`, or `"conditional"`) and, for a fixed start `g0 != 0` in cases DU and DT only, the contour integral
#' (`"contour"`): the double integral of Corollary 1 of the paper (the second form of Theorem 1), with the integral over
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
#' The computation runs on the engine's worker threads; see [ltratio_threads()].
#'
#' @param x Numeric vector or array of points. `-Inf`, `Inf` and `NaN` give 0, 1 (density 0) and `NaN`.
#' @param cas The case, one of `"DK"`, `"DU"`, `"DT"`, `"SK"`, `"SU"`, `"ST"`, `"DG"` (see the section Cases).
#' @param c The local-to-unity parameter \eqn{c}, \eqn{\phi = 1 + c/n}{phi = 1 + c/n}.
#' @param c0 The local null value \eqn{c_0}{c0}, \eqn{\phi_0 = 1 + c_0/n}{phi0 = 1 + c0/n}.
#' @param g0 The fixed start, \eqn{y_0 = g_0 \sigma \sqrt{n}}{y_0 = g0 sigma sqrt(n)} (cases DK, DU, DT).
#' @param method `"auto"`, `"single"` (or `"line"`), `"double"` (or `"conditional"`) or `"contour"` (a fixed start
#'   in cases DU, DT only); see the section Methods and ranges.
#' @param verified `FALSE` (the default): one evaluation by one method; `TRUE`: the checked computation (see the
#'   section Methods and ranges).
#' @return A numeric vector with the attributes (`dim`, `dimnames`, `names`) of `x`. A warning is given in a regime
#'   in which one evaluation is known to be brittle, and, with `verified = TRUE` and `method = "double"` or
#'   `"contour"`, if the integral reports at some points that its error estimate is not met.
#' @references Bao, Y., Broda, S. A. and Kan, R., "The Exact Limiting Distribution of the Dickey-Fuller t-Ratio
#'   and Some Related Statistics", working paper.
#'
#'   Elliott, G., Rothenberg, T. J. and Stock, J. H. (1996), "Efficient Tests for an Autoregressive Unit Root",
#'   Econometrica 64, 813-836.
#' @seealso [ltratio_icdf()], [ltratio_mom()], [coef_cdf()], [adfmax_cdf()], [ltratio_region()], [ltratio_threads()].
#' @examples
#' ltratio_cdf(-2.86157, "DU", -5)            # limiting power of the 5% test of tau_mu at c = -5
#' ltratio_cdf(-2.86157, "SU", -5)            # the same with a stationary start
#' ltratio_cdf(0, "DK", 2, c0 = 2)            # t-ratio centred at the true value
#' ltratio_pdf(c(-3, -2, -1), "DT", -5)       # density
#' @export
ltratio_cdf <- function(x, cas, c = 0, c0 = 0, g0 = 0, method = "auto", verified = FALSE) {
  dfl_dist(2L, x, cas, c, c0, g0, method, verified)
}

#' @rdname ltratio_cdf
#' @export
ltratio_pdf <- function(x, cas, c = 0, c0 = 0, g0 = 0, method = "auto", verified = FALSE) {
  dfl_dist(1L, x, cas, c, c0, g0, method, verified)
}

dfl_dist <- function(kind, x, cas, c, c0, g0, method, verified) {
  k <- dfl_case(cas)
  me <- dfl_method(method)
  ve <- dfl_flag(verified, "verified")
  r <- .Call(C_dfl_r_dist, kind, dfl_real(x, "x"), k$model, dfl_num(c, "c"), dfl_num(c0, "c0"),
             dfl_num(g0, "g0"), k$stat, me, ve)
  if (ve == 1L && me >= 2L && r[[2]] > 0L) {
    warning(sprintf("the %s integral did not meet its error estimate at %d of %d points",
                    if (me == 2L) "conditional" else "contour", r[[2]], length(x)), call. = FALSE)
  }
  if (nzchar(r[[3]])) warning(r[[3]], call. = FALSE)
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
#' @return The quantiles, with the attributes of `p` (or the list described under `info`). A warning is given if
#'   Newton's method did not converge at some points.
#' @inheritSection ltratio_cdf Cases
#' @seealso [ltratio_cdf()].
#' @examples
#' ltratio_icdf(c(0.01, 0.05, 0.10), "DU")   # asymptotic critical values of tau_mu
#' @export
ltratio_icdf <- function(p, cas, c = 0, c0 = 0, g0 = 0, method = "auto", tol = 1e-10, info = FALSE, verified = FALSE) {
  k <- dfl_case(cas)
  r <- .Call(C_dfl_r_icdf, dfl_real(p, "p"), k$model, dfl_num(c, "c"), dfl_num(c0, "c0"), dfl_num(g0, "g0"),
             k$stat, dfl_method(method), dfl_num(tol, "tol"), dfl_flag(verified, "verified"))
  if (nzchar(r[[5]])) warning(r[[5]], call. = FALSE)
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
#' the joint moment generating function of \eqn{(U, V)} (Proposition S1 of the Supplementary Material of the paper).
#' Not available for case DG.
#'
#' @inheritParams ltratio_cdf
#' @param p The number of moments, a positive integer.
#' @param info If `TRUE`, return a list with the moments `m` and the error estimates `err` of the integrals.
#' @return A numeric vector of length `p` (or the list described under `info`).
#' @inheritSection ltratio_cdf Cases
#' @seealso [ltratio_cdf()].
#' @examples
#' ltratio_mom(4, "DT")      # the first four moments of tau_tau
#' ltratio_mom(2, "DU")[2]   # 3.0559953450 (Nabeya, 1999)
#' @export
ltratio_mom <- function(p, cas, c = 0, c0 = 0, g0 = 0, info = FALSE) {
  if (!is.numeric(p) || length(p) != 1L || is.na(p) || p != round(p) || p < 1) {
    stop("p must be a positive integer", call. = FALSE)
  }
  k <- dfl_case(cas)
  r <- .Call(C_dfl_r_moments, as.integer(p), k$model, dfl_num(c, "c"), dfl_num(c0, "c0"), dfl_num(g0, "g0"), k$stat)
  names(r) <- c("m", "err")
  if (isTRUE(info)) r else r$m
}
