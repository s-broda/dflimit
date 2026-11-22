#' Limiting probability of a region of the plane
#'
#' The limit of the probability \eqn{P\{U \le g(V)\}}{P{U <= g(V)}} under \eqn{\phi = 1 + c/n}{phi = 1 + c/n}, no
#' drift, zero start, case DK (no deterministic terms; also the constant case with GLS demeaning), with \eqn{U} and
#' \eqn{V} as in [ltratio_cdf()]. Every test of the unit root in the limit experiment of case DK rejects on such a
#' region (the paper's section on the weighted-average-power test), and the first part of the conditional integral
#' representation of the paper, the probability of a region of the plane, gives its probability by one outer integral
#' over \eqn{V}.
#' Accurate to about 1e-11.
#'
#' `bnd` is a list with the element `type` and the parameters of the boundary:
#' \describe{
#'   \item{`"line"`}{\eqn{g(v) = b + a v}, elements `a` and `b` (the point-optimal test, \eqn{U - aV \le b}{U - aV
#'     <= b});}
#'   \item{`"parabola"`}{\eqn{g(v) = x \sqrt{v} + c_0 v}{g(v) = x sqrt(v) + c0 v}, elements `x` and `c0` (default 0):
#'     the t-ratio, the same probability as `ltratio_cdf(x, "DK", c, c0 = c0)`;}
#'   \item{`"wap"`}{the weighted-average-power test with the half-normal weight of scale \eqn{\omega}{omega} on
#'     \eqn{c < 0}, \eqn{(1 + \omega^2 V)^{-1/2} \mathrm{erfcx}(U/\sqrt{2 (V + 1/\omega^2)}) \ge k}{(1 + omega^2
#'     V)^(-1/2) erfcx(U/sqrt(2 (V + 1/omega^2))) >= k}, elements `k` and `omega`;}
#' }
#' and optionally `v0` (\eqn{g = +\infty}{g = +Inf} on \eqn{(0, v_0)}{(0, v0)}, whose probability
#' \eqn{P(V < v_0)}{P(V < v0)} is added) and `brks` (kinks of \eqn{g}, breakpoints of the outer rule). The boundary
#' type `"fun"` of the MATLAB programs (any function \eqn{g}) is not available: the C engine has the three families
#' above only.
#'
#' In place of the case a design of the general engine ([dfl_design]) may be given, with its start (`g0`, `start`)
#' and GLS detrending (`gls`): `dfl_trend(-1)` with a zero start is case DK, and every other design, a fixed start
#' `g0 != 0` (also with the case code `"DK"`), or `engine = "general"`, takes the general engine, for the types `"line"`
#' and `"parabola"` without `v0` and `brks` ([region_prob()]); the classical engine has the zero start only.
#'
#' @param bnd The boundary, a list as described above.
#' @param cas The case: `"DK"` only (with `g0` a fixed start), or a design ([dfl_design]).
#' @param engine `"auto"` (the default; the classical engine for case DK), `"classical"` or `"general"`.
#' @param start,g0,gls With a design: its start (`"zero"`, `"fixed"` with the fixed start `g0`, `"stationary"`) and
#'   GLS detrending (`NULL`: least squares); see [ltratio_cdf()].
#' @param c The local-to-unity parameter.
#' @param info If `TRUE`, return a list with the probability `P`, the error estimate `err` and the status `status`
#'   (nonzero: the error estimate is not met).
#' @return The probability (or the list described under `info`). A warning is given if the error estimate is not
#'   met.
#' @seealso [ltratio_cdf()].
#' @examples
#' # the weighted-average-power test of the paper: size 0.05
#' ltratio_region(list(type = "wap", k = 3.6658, omega = 10), "DK", 0)
#' # the point-optimal test against cbar = -7: size 0.05
#' ltratio_region(list(type = "line", a = -3.5, b = -0.2639), "DK", 0)
#' ltratio_region(list(type = "line", a = -3.5, b = -0.2639), dfl_perron("A", 0.5), 0)
#' @export
ltratio_region <- function(bnd, cas = "DK", c = 0, info = FALSE, engine = c("auto", "classical", "general"),
                           start = NULL, g0 = 0, gls = NULL) {
  engine <- dfl_engine(engine)
  if (dfl_is_design(cas)) {
    tg <- dfl_target(cas, g0, start, gls)
    # (the classical regions: case DK with a zero start)
    if (identical(tg$case, "DK") && tg$start == 0L && engine != "general") return(ltratio_region(bnd, "DK", c, info))
    if (engine == "classical") dfl_no_classical("the regions of the classical engine have case DK only, with a zero start")
    return(dfl_region_general(bnd, list(path = "router", spec = dfl_spec(tg, FALSE), engine = 2L), c, info))
  }
  if (!is.character(cas) || length(cas) != 1L || toupper(cas) != "DK") {
    stop("the case must be DK", call. = FALSE)
  }
  if (!is.null(start) || !is.null(gls)) {
    stop("start and gls apply to a design; the case code DK takes a fixed start by g0", call. = FALSE)
  }
  g0 <- dfl_num(g0, "g0")
  if (g0 != 0) {
    # (a fixed start: the classical regions have the zero start only; the router takes the general engine, and with
    # engine = "classical" it refuses)
    return(dfl_region_general(bnd, list(path = "router", spec = list(NULL, 1L, 1L, g0, NaN), engine = dfl_engine_code(engine)),
                              c, info))
  }
  if (engine == "general") {
    return(dfl_region_general(bnd, list(path = "router", spec = list(NULL, 1L, 0L, 0, NaN), engine = 2L), c, info))
  }
  if (!is.list(bnd) || is.null(bnd$type)) stop("bnd must be a list with the element type", call. = FALSE)
  el <- function(name, default = NULL) {
    v <- bnd[[name]]
    if (is.null(v)) {
      if (is.null(default)) stop(sprintf("a boundary of type %s needs the element %s", bnd$type, name), call. = FALSE)
      v <- default
    }
    dfl_num(v, paste0("bnd$", name))
  }
  type <- match(tolower(bnd$type), c("line", "parabola", "wap"))
  if (is.na(type)) {
    stop("the boundary type must be \"line\", \"parabola\" or \"wap\" (the type \"fun\" of the MATLAB programs is not available in the C engine)",
         call. = FALSE)
  }
  par <- switch(type, c(el("a"), el("b")), c(el("x"), el("c0", 0)), c(el("k"), el("omega")))
  v0 <- if (is.null(bnd$v0)) 0 else el("v0")
  brks <- if (is.null(bnd$brks)) double() else dfl_real(bnd$brks, "bnd$brks")
  if (!all(is.finite(brks) & brks > 0)) stop("bnd$brks must be positive and finite", call. = FALSE)
  if (!is.finite(v0) || v0 < 0) stop("bnd$v0 must be finite and nonnegative", call. = FALSE)
  r <- .Call(C_dfl_r_region, type, par, v0, brks, dfl_num(c, "c"))
  names(r) <- c("P", "err", "status")
  if (r$status != 0L) {
    warning(sprintf("the error estimate %.1e is above the tolerance; the result may be inaccurate", r$err), call. = FALSE)
  }
  if (isTRUE(info)) r else r$P
}

# the region with the general engine: the line and the parabola only, without v0 and brks
dfl_region_general <- function(bnd, pl, c, info) {
  if (!is.list(bnd) || is.null(bnd$type)) stop("bnd must be a list with the element type", call. = FALSE)
  type <- match(tolower(bnd$type), c("line", "parabola"))
  if (is.na(type)) {
    stop("with the general engine the boundary type must be \"line\" or \"parabola\" (its kernels are in closed form for boundaries quadratic in sqrt(V) only)",
         call. = FALSE)
  }
  if (!is.null(bnd$v0) || !is.null(bnd$brks)) stop("bnd$v0 and bnd$brks apply to the classical engine only", call. = FALSE)
  el <- function(name, default = NULL) {
    v <- bnd[[name]]
    if (is.null(v)) {
      if (is.null(default)) stop(sprintf("a boundary of type %s needs the element %s", bnd$type, name), call. = FALSE)
      v <- default
    }
    dfl_num(v, paste0("bnd$", name))
  }
  par <- if (type == 1L) c(el("a"), el("b")) else c(el("x"), el("c0", 0))
  P <- dfl_route_call("region", 0, pl, c(type, par, dfl_num(c, "c")))
  if (isTRUE(info)) list(P = P, err = NA_real_, status = NA_integer_) else P
}
