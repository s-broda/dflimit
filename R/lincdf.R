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
#' The probability is computed by the inversion formula of Gil-Pelaez,
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
#' with a case code and no `engine` given the classical engine computes as before; `dfl_trend(-1)` (case DK) and
#' `dfl_trend(1)` with `gls = -13.5` (case DG) are classical cases (with `engine = "auto"` the general engine for the law
#' of \eqn{V} in case DK and in case DG above the unit root), every other design takes the general engine
#' ([lin_cdf()]).
#'
#' @param b Numeric vector or array of points. `-Inf`, `Inf` and `NaN` give 0, 1 and `NaN`.
#' @param a The coefficient of \eqn{V}, a single number; `NaN` for the distribution of \eqn{V} (case DK).
#' @param cas The case, `"DK"` or `"DG"`, or a design ([dfl_design]).
#' @param engine `"auto"`, `"classical"` or `"general"`. Not given: the classical engine for a case code, `"auto"` for
#'   a design.
#' @param start,g0,gls With a design: its start (`"zero"`, `"fixed"` with the fixed start `g0`, `"stationary"`) and
#'   GLS detrending (`NULL`: least squares); see [ltratio_cdf()].
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
#' ltratio_lincdf(0, -3.5, dfl_perron("A", 0.5), -5)
#' @export
ltratio_lincdf <- function(b, a, cas, c = 0, cbar = -13.5, engine = c("auto", "classical", "general"), start = NULL,
                           g0 = 0, gls = NULL) {
  given <- !missing(engine)
  engine <- dfl_engine(engine)
  if (!is.numeric(a) || length(a) != 1L || (is.na(a) && !is.nan(a))) {
    stop("a must be a single number (NaN: the distribution of V)", call. = FALSE)
  }
  if (dfl_is_design(cas)) {
    tg <- dfl_target(cas, g0, start, gls)
    # (the classical cases of the linear form: DK, and DG detrended at -13.5, zero start)
    cls <- !is.null(tg$case) && tg$case %in% c("DK", "DG")
    if (engine == "classical" && !cls) dfl_no_classical("the linear form of the classical engine has the cases DK and DG only")
    if (engine == "general" || !cls) {
      G <- dfl_general_target(tg)
      return(dfl_det_finite("lincdf", G$D, b, c(as.double(a), dfl_num(c, "c")), dfl_opts(0, G$start, G$g0, G$gls)))
    }
    cas <- tg$case
    cbar <- -13.5
    given <- TRUE
  } else if (dfl_num(g0, "g0") != 0 || !is.null(start) || !is.null(gls)) {
    stop("start, g0 and gls apply to a design; the linear form of the cases DK and DG has a zero start", call. = FALSE)
  }
  k <- dfl_case(cas, c("DK", "DG"))
  if (given && engine != "classical") {
    pl <- list(spec = list(NULL, 1L, 0L, 0, NaN), engine = dfl_engine_code(engine))
    return(dfl_shape(dfl_route_call("lin", dfl_real(b, "b"), pl, c(as.double(a), k$model, dfl_num(c, "c"), dfl_num(cbar, "cbar"))), b))
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
#' statistic of the design ([pt_envelope()]). `engine` chooses between the two engines as in [ltratio_cdf()]: with a
#' case code and no `engine` given the classical engine computes as before.
#'
#' @param c Numeric vector or array of the local-to-unity parameters. `NaN` gives `NaN`.
#' @param cas The case, `"DK"`, `"DU"` (the same as `"DK"`) or `"DG"`, or a design ([dfl_design]).
#' @param engine `"auto"`, `"classical"` or `"general"`. Not given: the classical engine for a case code, `"auto"` for
#'   a design.
#' @param alpha The level of the tests.
#' @return The envelope, with the attributes of `c`.
#' @references Elliott, G., Rothenberg, T. J. and Stock, J. H. (1996), "Efficient Tests for an Autoregressive Unit
#'   Root", Econometrica 64, 813-836.
#' @seealso [ltratio_lincdf()], [ltratio_cdf()].
#' @examples
#' power_envelope(-7, "DK")                   # 0.4994: at c = -7 the envelope has power one half
#' power_envelope(c(-5, -13.5), "DG")
#' power_envelope(-10, dfl_perron("C", 0.5))
#' @export
power_envelope <- function(c, cas = "DK", alpha = 0.05, engine = c("auto", "classical", "general")) {
  given <- !missing(engine)
  engine <- dfl_engine(engine)
  if (dfl_is_design(cas)) {
    # (dfl_trend(-1) and dfl_trend(0): case DK; dfl_trend(1): case DG; the other designs: the envelope of the
    # point-optimal statistic of the design, pt_envelope())
    cas <- if (dfl_is_trend(cas, -1L) || dfl_is_trend(cas, 0L)) "DK" else if (dfl_is_trend(cas, 1L)) "DG" else cas
    if (dfl_is_design(cas)) {
      if (engine == "classical") dfl_no_classical("the envelope of the classical engine has the cases DK (DU) and DG only")
      return(pt_envelope(c, cas, alpha))
    }
    given <- TRUE
  }
  k <- dfl_case(cas, c("DK", "DU", "DG"))
  model <- if (k$model == 4L) 4L else 1L
  if (given && engine != "classical") {
    pl <- list(spec = list(NULL, 1L, 0L, 0, NaN), engine = dfl_engine_code(engine))
    return(dfl_shape(dfl_route_call("env", dfl_real(c, "c"), pl, c(model, dfl_num(alpha, "alpha"))), c))
  }
  P <- .Call(C_dfl_r_envelope, dfl_real(c, "c"), model, dfl_num(alpha, "alpha"))
  dfl_shape(P, c)
}
