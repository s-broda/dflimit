#' Limiting distribution of the linear form U - aV
#'
#' The limit of \eqn{P(U - aV \le b)}{P(U - aV <= b)} under \eqn{\phi = 1 + c/n}{phi = 1 + c/n}, zero start, no
#' drift, with \eqn{U} and \eqn{V} as in [ltratio_cdf()]: case `"DK"` (no deterministic terms, and the constant case
#' with GLS demeaning, whose limit is the same) or `"DG"` (the GLS-detrended trend case of Elliott, Rothenberg and
#' Stock (1996), detrended at `cbar`). The point-optimal test of Elliott, Rothenberg and Stock against
#' \eqn{\bar c < 0}{cbar < 0} rejects for \eqn{U - aV < b}: in case DK with \eqn{a = \bar c/2}{a = cbar/2} (its log
#' likelihood ratio is \eqn{\bar c U - \bar c^2 V/2}{cbar U - cbar^2 V/2}), in case DG, with the data detrended at
#' \eqn{\bar c}{cbar}, with \eqn{a = \bar c^2/(2(\bar c - 1))}{a = cbar^2/(2 (cbar - 1))} (it rejects for small
#' \eqn{\bar c^2 V + (1 - \bar c) J(1)^2}{cbar^2 V + (1 - cbar) J(1)^2}, \eqn{J(1)^2 = 2U + 1}); see
#' [power_envelope()]. `a = NaN` gives \eqn{P(V \le b)}{P(V <= b)} (case DK).
#'
#' By default (`method = "auto"` with `engine = "auto"`, a finite \eqn{a \ne 0}{a != 0}) the probability is the single
#' inversion of the joint characteristic function along the ray of the half-plane, the method of Yong Bao and Raymond
#' Kan, "Exact Finite-Sample Distribution and Moments of the Least Squares Estimator in First-Order Autoregressions: O(1)
#' Algorithms", with the tail beyond \eqn{S = 400/|w|}{S = 400/|w|} (\eqn{w = 1/2 + b}) completed by integration by
#' parts: in case DK for \eqn{-100 \le c \le 4}{-100 <= c <= 4}, and for a design that is one of the six classical
#' cases (`dfl_trend(-1)`, `dfl_trend(0)`, `dfl_trend(1)` with any start) for \eqn{-100 \le c \le 2}{-100 <= c <= 2}
#' (as `ltratio_lincdf` of DFLimit.jl and the MATLAB programs). `tol` is its absolute tolerance, `err = TRUE` returns
#' the error estimates as well, and a warning names the points at which they exceed it (`verified = TRUE` recomputes
#' those by the engines). `method = "inversion"` forces it (any `c`); `method = "engine"`, case DG, `a = NaN`, `engine =
#' "classical"` or `"general"` and the points outside those ranges take the computation of before, which reproduces its
#' values and is described next.
#'
#' The engines compute the probability by the inversion formula of Gil-Pelaez,
#' \eqn{P(X < b) = 1/2 - (1/\pi) \int_0^\infty \mathrm{Im}[\varphi(u) e^{-iub}]/u \, du}{P(X < b) = 1/2 - (1/pi)
#' int_0^Inf Im[phi(u) exp(-iub)]/u du}, from the joint moment generating function of \eqn{(U, V)}, by the
#' trapezoidal rule in \eqn{u = \sinh v}{u = sinh(v)} with 60000 steps up to \eqn{u = 4000}, as `ltratio_lincdf` of
#' the MATLAB and Julia programs.
#'
#' Where the transform has not fallen below 5e-9 by the end of that grid, or `c > 0` (\eqn{|a|} below about 0.3 to 3,
#' depending on `c`, which includes `a = 0`, the law of \eqn{U} alone), that inversion is not accurate, and the
#' probability comes from a conditional integral instead: in case DK that of the region
#' \eqn{U \le aV + b}{U <= aV + b} ([ltratio_region()]), in case DG the conditional integral given \eqn{U + 1/2}, whose
#' law is a gamma law of shape 1/2 under every tilt in \eqn{V}; \eqn{P(V \le b)}{P(V <= b)} takes a longer grid for
#' `c <= -50` and a finer one for `c > 2`. The values are accurate to about 1e-10 and lie in [0, 1]. `c` must be at least -550, in case DG at most
#' about 3.5 (at `cbar = -13.5`; less for `cbar` nearer 0), and for `a = NaN` at most about 7.
#'
#' In place of the case a design of the general engine ([dfl_design]) may be given, with its start (`g0`, `start`)
#' and GLS detrending (`gls`); `engine` chooses between the two engines as in [ltratio_cdf()] (its section Engines):
#' `dfl_trend(-1)` (case DK) and `dfl_trend(1)` with `gls = -13.5` (case DG), with a zero start, are the classical cases;
#' by default (`engine = "auto"`) the rule of the C router chooses for them (the general engine for the law of \eqn{V} in
#' case DK and in case DG with `a != 0` above the unit root for `cbar >= -5`, and beyond `c = 2`, the classical engine
#' elsewhere), and `engine = "classical"` gives the classical values throughout; every other design takes the general
#' engine ([lin_cdf()]).
#'
#' @param b Numeric vector or array of points. `-Inf`, `Inf` and `NaN` give 0, 1 and `NaN`.
#' @param a The coefficient of \eqn{V}, a single number; `NaN` for the distribution of \eqn{V} (case DK).
#' @param cas The case, `"DK"` or `"DG"`, or a design ([dfl_design]).
#' @param engine `"auto"` (the default), `"classical"` or `"general"`; see the section Engines of [ltratio_cdf()].
#' @param start,g0,gls With a design: its start (`"zero"`, `"fixed"` with the fixed start `g0`, `"stationary"`) and
#'   GLS detrending (`NULL`: least squares); see [ltratio_cdf()].
#' @param c The local-to-unity parameter, finite.
#' @param cbar The detrending parameter \eqn{\bar c}{cbar} of case DG (-13.5: the DF-GLS tests); not used in case DK.
#' @param method `"auto"` (the default), `"inversion"` (the single inversion of Bao and Kan) or `"engine"` (the
#'   engines' computation of before).
#' @param tol The absolute tolerance of the single inversion.
#' @param verified `TRUE`: the values whose error estimate exceeds `tol` are recomputed by the engines.
#' @param err `TRUE`: a list of the values (`value`) and their error estimates (`err`; `NaN` where the engines computed).
#' @return The probabilities, with the attributes of `b` (with `err = TRUE` a list, `value` and `err`).
#' @references Elliott, G., Rothenberg, T. J. and Stock, J. H. (1996), "Efficient Tests for an Autoregressive Unit
#'   Root", Econometrica 64, 813-836.
#' @seealso [power_envelope()], [ltratio_region()], [ltratio_cdf()].
#' @examples
#' ltratio_lincdf(-0.2639, -3.5, "DK", 0)    # 0.05: the size of the point-optimal test at cbar = -7
#' ltratio_lincdf(-0.2639, -3.5, "DK", -7)   # its power at c = -7
#' ltratio_lincdf(c(0.1, 0.5, 2), NaN, "DK") # P(V <= b) at the unit root
#' ltratio_lincdf(0, -3.5, dfl_perron("A", 0.5), -5)
#' @export
ltratio_lincdf <- function(b, a, cas, c = 0, cbar = -13.5, engine = c("auto", "classical", "general"), start = NULL,
                           g0 = 0, gls = NULL, method = "auto", tol = 1e-12, verified = FALSE, err = FALSE) {
  engine <- dfl_engine(engine)
  if (!is.numeric(a) || length(a) != 1L || (is.na(a) && !is.nan(a))) {
    stop("a must be a single number (NaN: the distribution of V)", call. = FALSE)
  }
  im <- dfl_imethod(method, "ltratio_lincdf")
  ve <- dfl_flag(verified, "verified")
  bv <- dfl_real(b, "b")
  if (dfl_is_design(cas)) {
    tg <- dfl_target(cas, g0, start, gls)
    # (the router: the single inversion of Bao and Kan for a design that is one of the six classical cases, the classical
    # linear form for cases DK and DG with a zero start, the general engine otherwise; the general engine's points at
    # -Inf, Inf, NaN: 0, 1, NaN)
    if (engine == "classical" && is.null(tg$case)) dfl_no_classical("the linear form of the classical engine has the cases DK and DG only")
    pl <- list(spec = dfl_spec(tg, FALSE), engine = dfl_engine_code(engine))
    out <- rep(NaN, length(bv))
    out[!is.na(bv) & bv > 0] <- 1
    out[!is.na(bv) & bv < 0] <- 0
    e <- ifelse(is.finite(bv), NaN, 0)
    ok <- is.finite(bv)
    if (any(ok)) {
      r <- dfl_route_call2("lind", bv[ok], pl, c(as.double(a), dfl_num(c, "c"), im, dfl_tol(tol), ve))
      out[ok] <- r$v
      e[ok] <- r$e
    }
    return(dfl_with_err(out, e, b, err))
  } else if (dfl_num(g0, "g0") != 0 || !is.null(start) || !is.null(gls)) {
    stop("start, g0 and gls apply to a design; the linear form of the cases DK and DG has a zero start", call. = FALSE)
  }
  k <- dfl_case(cas, c("DK", "DG"))
  if (engine != "classical" || im == 5L) {
    pl <- list(spec = list(NULL, 1L, 0L, 0, NaN), engine = dfl_engine_code(engine))
    r <- dfl_route_call2("lin", bv, pl, c(as.double(a), k$model, dfl_num(c, "c"), dfl_num(cbar, "cbar"), im, dfl_tol(tol), ve))
    return(dfl_with_err(r$v, r$e, b, err))
  }
  P <- .Call(C_dfl_r_lincdf, bv, as.double(a), k$model, dfl_num(c, "c"), dfl_num(cbar, "cbar"))
  dfl_with_err(P, rep(NaN, length(P)), b, err)
}

#' Asymptotic power envelope of Elliott, Rothenberg and Stock (1996)
#'
#' The power at `c` of the point-optimal level-`alpha` test against `c` itself, the envelope of the powers of all
#' tests of the unit root, zero start, no drift: case `"DK"` (no deterministic terms; also the constant case,
#' `"DU"`, where the point-optimal test uses the GLS-demeaned data and the envelope is the same curve) and case
#' `"DG"` (the trend case, the point-optimal test against `c` with the data GLS-detrended at `c`).
#'
#' The test against \eqn{\bar c}{cbar} rejects for small \eqn{U - aV} ([ltratio_lincdf()]), with
#' \eqn{a = \bar c/2}{a = cbar/2} in case DK and \eqn{a = \bar c^2/(2(\bar c - 1))}{a = cbar^2/(2 (cbar - 1))} in case
#' DG, where the limit of the invariant statistic of Elliott, Rothenberg and Stock is
#' \eqn{\bar c^2 V + (1 - \bar c) J(1)^2}{cbar^2 V + (1 - cbar) J(1)^2}, \eqn{J(1)^2 = 2U + 1}: its critical value is the `alpha`-quantile at `c = 0`, found by bisection on
#' \eqn{[-50, 50]} to the width 1e-9, and the envelope is its power at \eqn{\bar c = c}{cbar = c}; `alpha` for
#' `c >= 0`. Near `c = 0`, where `a` is small and that inversion is not accurate, the critical value and the power come
#' from the conditional integral of [ltratio_lincdf()], the critical value by bisection in \eqn{b + 1/2} to
#' 1e-12 of its scale. As `power_envelope` of the MATLAB and Julia programs.
#'
#' In place of the case a design of the general engine ([dfl_design]) may be given: `dfl_trend(-1)` and
#' `dfl_trend(0)` are case DK, `dfl_trend(1)` case DG, and every other design has the envelope of the point-optimal
#' statistic of the design ([pt_envelope()]). `engine` chooses between the two engines as in [ltratio_cdf()]; by default
#' (`engine = "auto"`) the rule of the C router keeps the classical engine for the cases DK and DG (1.6 to 11 times
#' faster; its critical value is bisected to 1e-9).
#'
#' In case DK (and DU) the default (`method = "auto"` with `engine = "auto"`, `c <= -1e-4`) computes the critical value
#' and the power by the single inversion of Bao and Kan (see [ltratio_lincdf()]), the critical value by regula falsi
#' to 1e-12; `tol` is the tolerance of each inversion, `err = TRUE` returns the error estimates (the sum of the two
#' inversions' estimates), and a warning names the points at which they exceed `tol` (`verified = TRUE` recomputes those
#' by the engines). `method = "inversion"` forces it; `method = "engine"`, case DG, `engine = "classical"` or `"general"`
#' and \eqn{-10^{-4} < c < 0}{-1e-4 < c < 0} take the computation described above, which reproduces the values of
#' earlier versions.
#'
#' @param c Numeric vector or array of the local-to-unity parameters. `NaN` gives `NaN`.
#' @param cas The case, `"DK"`, `"DU"` (the same as `"DK"`) or `"DG"`, or a design ([dfl_design]).
#' @param engine `"auto"` (the default), `"classical"` or `"general"`; see the section Engines of [ltratio_cdf()].
#' @param alpha The level of the tests.
#' @param method `"auto"` (the default), `"inversion"` (the single inversion of Bao and Kan) or `"engine"` (the
#'   engines' computation of before).
#' @param tol The absolute tolerance of each single inversion.
#' @param verified `TRUE`: the values whose error estimate exceeds `tol` are recomputed by the engines.
#' @param err `TRUE`: a list of the values (`value`) and their error estimates (`err`; `NaN` where the engines computed).
#' @return The envelope, with the attributes of `c` (with `err = TRUE` a list, `value` and `err`).
#' @references Elliott, G., Rothenberg, T. J. and Stock, J. H. (1996), "Efficient Tests for an Autoregressive Unit
#'   Root", Econometrica 64, 813-836.
#' @seealso [ltratio_lincdf()], [ltratio_cdf()].
#' @examples
#' power_envelope(-7, "DK")                   # 0.4994: at c = -7 the envelope has power one half
#' power_envelope(c(-5, -13.5), "DG")
#' power_envelope(-10, dfl_perron("C", 0.5))
#' @export
power_envelope <- function(c, cas = "DK", alpha = 0.05, engine = c("auto", "classical", "general"), method = "auto",
                           tol = 1e-12, verified = FALSE, err = FALSE) {
  engine <- dfl_engine(engine)
  im <- dfl_imethod(method, "power_envelope")
  ve <- dfl_flag(verified, "verified")
  if (dfl_is_design(cas)) {
    # (dfl_trend(-1) and dfl_trend(0): case DK; dfl_trend(1): case DG; the other designs: the envelope of the
    # point-optimal statistic of the design, pt_envelope())
    cas <- if (dfl_is_trend(cas, -1L) || dfl_is_trend(cas, 0L)) "DK" else if (dfl_is_trend(cas, 1L)) "DG" else cas
    if (dfl_is_design(cas)) {
      if (engine == "classical") dfl_no_classical("the envelope of the classical engine has the cases DK (DU) and DG only")
      if (im != 0L) stop("power_envelope: method applies to the cases DK (DU) and DG; the envelope of a design is pt_envelope()", call. = FALSE)
      v <- pt_envelope(c, cas, alpha)
      return(dfl_with_err(v, rep(NaN, length(v)), c, err))
    }
  }
  k <- dfl_case(cas, c("DK", "DU", "DG"))
  model <- if (k$model == 4L) 4L else 1L
  if (engine != "classical" || im == 5L) {
    pl <- list(spec = list(NULL, 1L, 0L, 0, NaN), engine = dfl_engine_code(engine))
    r <- dfl_route_call2("env", dfl_real(c, "c"), pl, c(model, dfl_num(alpha, "alpha"), im, dfl_tol(tol), ve))
    return(dfl_with_err(r$v, r$e, c, err))
  }
  P <- .Call(C_dfl_r_envelope, dfl_real(c, "c"), model, dfl_num(alpha, "alpha"))
  dfl_with_err(P, rep(NaN, length(P)), c, err)
}
