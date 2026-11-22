# argument conversion shared by the user-level functions (the checks of the ranges are the engine's, dfl_check_args)

# the case code: model 1, 2, 3 (no deterministic terms, intercept, intercept and trend), 4 (the GLS-detrended trend
# case DG), and the start (stat = 1: stationary)
dfl_case <- function(cas, allowed = c("DK", "DU", "DT", "SK", "SU", "ST", "DG")) {
  if (!is.character(cas) || length(cas) != 1L || is.na(cas)) {
    stop("the case must be one of ", paste(allowed, collapse = ", "), call. = FALSE)
  }
  cas <- toupper(cas)
  if (!cas %in% allowed) stop("the case must be one of ", paste(allowed, collapse = ", "), call. = FALSE)
  if (cas == "DG") return(list(case = cas, model = 4L, stat = 0L))
  list(case = cas, model = match(substr(cas, 2, 2), c("K", "U", "T")), stat = as.integer(substr(cas, 1, 1) == "S"))
}

# a logical flag (TRUE or FALSE), as an integer for the engine
dfl_flag <- function(v, name) {
  if (!is.logical(v) || length(v) != 1L || is.na(v)) stop(name, " must be TRUE or FALSE", call. = FALSE)
  as.integer(v)
}

# a real scalar (finiteness and range are checked by the engine)
dfl_num <- function(v, name) {
  if (!is.numeric(v) || length(v) != 1L || is.na(v)) stop(name, " must be a single number", call. = FALSE)
  as.double(v)
}

# the method: 0 automatic, 1 the line (single) integral, 2 the conditional (double) integral, 3 the contour integral
dfl_method <- function(method) {
  msg <- paste0("method must be \"auto\", \"single\" (\"line\"), \"double\" (\"conditional\") or \"contour\" ",
                "(\"paths\", the general engine's inversion on paths around the roots, with a design or engine = \"general\")")
  if (!is.character(method) || length(method) != 1L) stop(msg, call. = FALSE)
  m <- match(tolower(method), c("auto", "single", "line", "double", "conditional", "contour", "paths"))
  if (is.na(m)) stop(msg, call. = FALSE)
  c(0L, 1L, 1L, 2L, 2L, 3L, 4L)[m]
}

dfl_real <- function(x, name) {
  if (!(is.numeric(x) || is.logical(x) && all(is.na(x)))) stop(name, " must be numeric", call. = FALSE)
  as.double(x)
}

# the attributes of the input (dim, dimnames, names) on the output
dfl_shape <- function(out, x) {
  if (!is.null(dim(x))) {
    dim(out) <- dim(x)
    dimnames(out) <- dimnames(x)
  } else {
    names(out) <- names(x)
  }
  out
}

# evaluate expr giving each distinct warning once (the calls of ltratio_cdf inside coef_cdf and coef_icdf would repeat
# the warning of a brittle regime at every point)
dfl_once <- function(expr) {
  seen <- character(0)
  withCallingHandlers(expr, warning = function(w) {
    m <- conditionMessage(w)
    if (!(m %in% seen)) {
      seen <<- c(seen, m)
      warning(m, call. = FALSE)
    }
    invokeRestart("muffleWarning")
  })
}

# ---- the general engine and the router over the two engines ----

# warnings of the engines (a message may hold several, one per line)
dfl_warn <- function(msgs) {
  for (m in msgs) {
    if (!nzchar(m)) next
    for (l in strsplit(m, "\n", fixed = TRUE)[[1L]]) warning(l, call. = FALSE)
  }
}

# the engine: "auto" (the default), "classical" or "general", as its code 0, 1, 2
dfl_engine <- function(engine) {
  if (identical(engine, c("auto", "classical", "general"))) return("auto")
  if (!is.character(engine) || length(engine) != 1L || is.na(engine) ||
      !tolower(engine) %in% c("auto", "classical", "general")) {
    stop("engine must be \"auto\", \"classical\" or \"general\"", call. = FALSE)
  }
  tolower(engine)
}
dfl_engine_code <- function(engine) match(engine, c("auto", "classical", "general")) - 1L

# the start of a design: NULL (fixed where the fixed start g0 is not 0, zero otherwise), "zero", "fixed", "stationary";
# the code 0, 1, 2 (a fixed start with g0 = 0 is the zero start)
dfl_start <- function(start, g0, g0name = "gamma0") {
  if (is.null(start)) return(if (g0 != 0) 1L else 0L)
  if (!is.character(start) || length(start) != 1L || is.na(start) ||
      is.na(k <- match(tolower(start), c("zero", "fixed", "stationary")))) {
    stop("start must be \"zero\", \"fixed\" or \"stationary\"", call. = FALSE)
  }
  k <- k - 1L
  if (k == 0L && g0 != 0) stop(sprintf("start = \"zero\" and %s != 0 (a fixed start) exclude each other", g0name), call. = FALSE)
  if (k == 2L && g0 != 0) stop(sprintf("start = \"stationary\" and %s != 0 (a fixed start) exclude each other", g0name), call. = FALSE)
  if (k == 1L && g0 == 0) k <- 0L
  k
}

# GLS detrending at abar = 1 + gls/n (gls < 0), or least squares (NULL; NaN for the engine)
dfl_gls <- function(gls, name = "gls") {
  if (is.null(gls)) return(NaN)
  v <- dfl_num(gls, name)
  if (!(v < 0)) stop(name, " must be negative (GLS detrending at 1 + ", name, "/n), or NULL for least squares", call. = FALSE)
  v
}

# the options of the general engine's functions: c(c0, start, gamma0, gls, verified)
dfl_opts <- function(c0 = 0, start = NULL, gamma0 = 0, gls = NULL, verified = FALSE, method = NULL) {
  c0 <- dfl_num(c0, "c0")
  gamma0 <- dfl_num(gamma0, "gamma0")
  o <- c(c0, dfl_start(start, gamma0), gamma0, dfl_gls(gls), dfl_flag(verified, "verified"))
  if (!is.null(method)) o <- c(o, dfl_gmethod(method))
  o
}
# the method of the general engine's t-ratio (dist_cdf, dist_pdf, dist_quantile): 0 "auto", 1 "line" (the line integral
# with the explicit kernel), 2 "paths" (the inversion on paths around the roots)
dfl_gmethod <- function(method) {
  msg <- "method must be \"auto\", \"line\" (the line integral with the explicit kernel) or \"paths\" (the inversion on paths around the roots)"
  if (!is.character(method) || length(method) != 1L) stop(msg, call. = FALSE)
  m <- match(tolower(method), c("auto", "line", "paths"))
  if (is.na(m)) stop(msg, call. = FALSE)
  m - 1L
}

# a function of the general engine (dfl_r_det): the value, with the warnings of the engine given
dfl_det_call <- function(what, D, x, args, opts = c(0, 0, 0, NaN, 0)) {
  dfl_check_design(D)
  r <- .Call(C_dfl_r_det, what, D, as.double(x), as.double(args), as.double(opts))
  dfl_warn(r[[3]])
  r
}
# the same at the finite points of x only (the others get lo, hi and NaN at -Inf, Inf, NaN)
dfl_det_finite <- function(what, D, x, args, opts = c(0, 0, 0, NaN, 0), lo = 0, hi = 1) {
  dfl_check_design(D)
  xv <- dfl_real(x, "x")
  out <- rep(NaN, length(xv))
  out[!is.na(xv) & xv > 0] <- hi
  out[!is.na(xv) & xv < 0] <- lo
  ok <- is.finite(xv)
  if (any(ok)) out[ok] <- dfl_det_call(what, D, xv[ok], args, opts)[[1]]
  dfl_shape(out, x)
}

# the classical case of a design with its start (0, 1, 2) and GLS detrending (NaN: least squares): "DK", ..., "ST",
# "DG" (dfl_trend(1) GLS-detrended at -13.5 with a zero start) or NULL (as the router of the C engine decides it)
dfl_is_trend <- function(D, p) {
  K <- p + 1L
  deg <- max(p, 0L)
  !D$first && !D$loc && D$K == K && length(D$br) == 2L && D$br[1] == 0 && D$br[2] == 1 && nrow(D$C) == deg + 1L &&
    (K == 0L || isTRUE(all(D$C == diag(1, deg + 1L, K))))
}
dfl_design_case <- function(D, start, gls) {
  for (p in -1:1) {
    if (!dfl_is_trend(D, p)) next
    if (is.nan(gls)) return(paste0(if (start == 2L) "S" else "D", c("K", "U", "T")[p + 2L]))
    if (p == 1L && start == 0L && gls == -13.5) return("DG")
    return(NULL)
  }
  NULL
}

# the computation of a function of the classical engine: a case code, or a design with its start (g0 a fixed start,
# start, gls); case: the classical case (NULL where the design is none), D: the design (NULL for a case code)
dfl_target <- function(cas, g0, start, gls, allowed = c("DK", "DU", "DT", "SK", "SU", "ST", "DG")) {
  g0 <- dfl_num(g0, "g0")
  if (dfl_is_design(cas)) {
    st <- dfl_start(start, g0, "g0")
    gl <- dfl_gls(gls)
    return(list(case = dfl_design_case(cas, st, gl), D = cas, start = st, g0 = if (st == 1L) g0 else 0, gls = gl))
  }
  if (!is.null(start) || !is.null(gls)) {
    stop("start and gls apply to a design; with a case code the code names the start (D: fixed, S: stationary)",
         call. = FALSE)
  }
  k <- dfl_case(cas, allowed)
  list(case = k$case, D = NULL, start = if (k$stat == 1L) 2L else if (g0 != 0) 1L else 0L, g0 = g0,
       gls = if (k$model == 4L) -13.5 else NaN)
}

# the spec of the router: list(design or NULL, model, start, g0, gls); the case code where the computation is a
# classical case and case is TRUE
dfl_spec <- function(tg, case = TRUE) {
  if (case && !is.null(tg$case)) {
    k <- dfl_case(tg$case)
    return(list(NULL, k$model, tg$start, as.double(tg$g0), NaN))
  }
  list(tg$D, 0L, tg$start, as.double(tg$g0), as.double(tg$gls))
}

# whether the classical engine accepts the arguments of a classical case (its checks of c, c0, g0, the start, and
# of the method), as the router of the C engine decides it for a design
dfl_classical_valid <- function(tg, c, c0, method = 0L) {
  if (is.null(tg$case)) return(FALSE)
  k <- dfl_case(tg$case)
  if (nzchar(.Call(C_dfl_r_check, k$model, as.double(c), as.double(c0), as.double(tg$g0), k$stat))) return(FALSE)
  if (method == 4L) return(FALSE)        # (the general engine's inversion on paths)
  if (k$model == 4L && method == 2L) return(FALSE)
  if (method == 3L && !(k$model %in% 2:3 && tg$g0 != 0 && k$stat == 0L)) return(FALSE)
  TRUE
}

# the path of a computation: "legacy" (the classical engine's own entry points, as before the router: with a case code
# and no engine named, engine = "classical", or a method or the checked computation with engine = "auto") or "router"
# (dfl_engine_* of the C engine, with the spec and the engine code)
dfl_plan <- function(tg, engine, c, c0, method = 0L, verified = 0L) {
  legacy <- function() list(path = "legacy", case = tg$case, g0 = tg$g0)
  router <- function(case, e = engine) list(path = "router", spec = dfl_spec(tg, case), engine = dfl_engine_code(e))
  if (is.null(tg$D)) {
    # (a case code: the router's rule, as for a design; engine = "classical" the classical entry points, with the values
    # of the earlier versions; a method other than "auto" and verified = TRUE keep the classical engine, as in Julia)
    if (engine == "classical") return(legacy())
    if (engine == "auto" && (method != 0L || verified != 0L)) return(legacy())
    return(router(TRUE))
  }
  if (engine == "general") return(router(FALSE))
  if (is.null(tg$case)) {
    if (engine == "classical") dfl_no_classical()
    return(router(FALSE, "general"))
  }
  if (engine == "classical") return(legacy())
  if (!dfl_classical_valid(tg, c, c0, method)) return(router(FALSE, "general"))
  if (method != 0L || verified != 0L) return(legacy())
  router(TRUE)
}
dfl_no_classical <- function(what = "the classical engine has the cases DK, DU, DT, SK, SU, ST and DG only") {
  stop("engine = \"classical\": ", what, " (dfl_trend(-1), dfl_trend(0), dfl_trend(1), the last GLS-detrended at ",
       "gls = -13.5 for case DG)", call. = FALSE)
}

# a function of the router (dfl_r_route): the value, with the warnings of the engines given
dfl_route_call <- function(what, x, plan, args) {
  r <- .Call(C_dfl_r_route, what, as.double(x), plan$spec, as.double(args), plan$engine)
  dfl_warn(r[[2]])
  r[[1]]
}

# the design, start, fixed start and detrending of a computation for the general engine's own functions (a case code:
# its trend; case DG: the linear trend GLS-detrended at -13.5)
dfl_general_target <- function(tg) {
  if (!is.null(tg$D)) return(list(D = tg$D, start = c("zero", "fixed", "stationary")[tg$start + 1L], g0 = tg$g0,
                                  gls = if (is.nan(tg$gls)) NULL else tg$gls))
  k <- dfl_case(tg$case)
  list(D = dfl_trend(if (k$model == 4L) 1L else k$model - 2L), start = c("zero", "fixed", "stationary")[tg$start + 1L],
       g0 = tg$g0, gls = if (k$model == 4L) -13.5 else NULL)
}
