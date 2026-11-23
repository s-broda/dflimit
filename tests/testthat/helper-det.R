# the reference values of the general engine and of the router (ref/det_ref.txt, the lines design, val, route, err and
# rerr of c/test/det_ref.txt of dflimit-matlab, copied by tools/det-ref-subset.R; the values of the Julia engine
# DFLimit.jl)

det_ref_lines <- function() {
  l <- readLines(test_path("ref", "det_ref.txt"))
  l[!startsWith(l, "#")]
}

# the designs of the lines 'design id spec...', built with the constructors of the package (the design 'raw np deg K
# loc br... C...' as a list filled by hand, a design given by its coefficients)
det_designs <- function(lines = det_ref_lines()) {
  fd <- list(sincos = function(r) c(sin(2 * pi * r), cos(2 * pi * r)),
             perronC04 = function(r) c(1, r, if (r > 0.4) 1 else 0, max(r - 0.4, 0)),
             expr = function(r) exp(r),
             trig = function(r) c(1, sin(2 * pi * r), cos(2 * pi * r)),
             rlin = function(r) r)
  fb <- list(sincos = numeric(0), perronC04 = 0.4, expr = numeric(0), trig = numeric(0), rlin = numeric(0))
  D <- list()
  for (l in lines[startsWith(lines, "design ")]) {
    w <- strsplit(l, " ", fixed = TRUE)[[1]]
    id <- as.integer(w[2])
    D[[id]] <- switch(w[3],
      trend = dfl_trend(as.integer(w[4])),
      ltrend = dfl_ltrend(as.integer(w[4])),
      perron = dfl_perron(w[4], as.numeric(w[5 + seq_len(as.integer(w[5]))])),
      fdesign = dfl_fdesign(fd[[w[4]]], fb[[w[4]]]),
      matrix = dfl_design_from_matrix(dfl_regmat(D[[as.integer(w[4])]], as.integer(w[5]))),
      first = dfl_detrend_first(D[[as.integer(w[4])]]),
      raw = {
        v <- as.numeric(w[-(1:3)])
        np <- as.integer(v[1]); deg <- as.integer(v[2]); K <- as.integer(v[3])
        structure(list(br = v[5 + 0:np], C = matrix(v[5 + np + seq_len(np * (deg + 1) * K)], deg + 1, K * np),
                       K = K, loc = v[4] != 0, first = FALSE, name = "a design given by its coefficients"),
                  class = "dfl_design")
      })
  }
  D
}

# a line 'val func id cbar start g0 c c0 verified par n x... nv v... tol [brittle]'
det_parse_val <- function(l) {
  w <- strsplit(l, " ", fixed = TRUE)[[1]]
  num <- function(i) as.numeric(w[i])
  n <- as.integer(w[11])
  nv <- as.integer(w[12 + n])
  list(fn = w[2], id = as.integer(w[3]), cb = num(4), st = as.integer(w[5]), g0 = num(6), c = num(7), c0 = num(8),
       ver = w[9] == "1", par = num(10), x = num(11 + seq_len(n)), want = num(12 + n + seq_len(nv)),
       tol = num(13 + n + nv), brittle = length(w) > 13 + n + nv)
}
# a line 'err func id cbar start g0 c c0 verified par n x... | text': the arguments of 'val', the call must fail with
# an error whose message contains text
det_parse_err <- function(l) {
  parts <- strsplit(l, " | ", fixed = TRUE)[[1]]
  w <- strsplit(parts[1], " ", fixed = TRUE)[[1]]
  num <- function(i) as.numeric(w[i])
  n <- as.integer(w[11])
  list(fn = w[2], id = as.integer(w[3]), cb = num(4), st = as.integer(w[5]), g0 = num(6), c = num(7), c0 = num(8),
       ver = w[9] == "1", par = num(10), x = num(11 + seq_len(n)), text = parts[2])
}
det_start <- function(st) c("zero", "fixed", "stationary")[st + 1L]
det_gls <- function(cb) if (is.nan(cb)) NULL else cb

# the value of a line 'val' by the functions of the general engine ('cdf@line', 'pdf@paths', 'quant@auto': the method of
# the t-ratio)
det_eval_val <- function(v, D) {
  d <- D[[v$id]]
  st <- det_start(v$st)
  gl <- det_gls(v$cb)
  x <- v$x
  # ('cdf#thin': the rule of one evaluation of the paths)
  fr <- strsplit(v$fn, "#", fixed = TRUE)[[1]]
  ru <- if (length(fr) > 1L) fr[2] else "auto"
  fm <- strsplit(fr[1], "@", fixed = TRUE)[[1]]
  me <- if (length(fm) > 1L) fm[2] else "auto"
  switch(fm[1],
    cdf = dist_cdf(x, d, v$c, v$c0, st, v$g0, gl, v$ver, method = me, rule = ru),
    pdf = dist_pdf(x, d, v$c, v$c0, st, v$g0, gl, v$ver, method = me, rule = ru),
    quant = dist_quantile(x, d, v$c, v$c0, st, v$g0, gl, x0 = v$par, tol = 1e-9, method = me, rule = ru),
    coefcdf = coef_cdf(x, d, v$c, v$c0, v$g0, verified = v$ver, engine = "general", start = st, gls = gl),
    coefpdf = coef_pdf(x, d, v$c, v$c0, st, v$g0, gl, v$ver),
    coefq = coef_quantile(x, d, v$c, v$c0, st, v$g0, gl, y0 = v$par, tol = 1e-9),
    coefgp = coef_cdf_gp(x, d, v$c, v$c0, st, v$g0, gl),
    mom = dist_moments(length(x), d, v$c, v$c0, st, v$g0, gl),
    lpow = local_power(x, d, statistic = if (v$par == 1) "coef" else "tratio"),
    lincdf = lin_cdf(x, v$par, d, v$c, st, v$g0, gl),
    linpdf = lin_pdf(x, v$par, d, v$c, st, v$g0, gl),
    linq = lin_quantile(x, v$par, d, v$c, st, v$g0, gl, b0 = 0, tol = 1e-10),
    adfcdf = adfmax_cdf(x, d, v$c, engine = "general", start = st, g0 = v$g0),
    adfpdf = adfmax_pdf(x, d, v$c, engine = "general", start = st, g0 = v$g0),
    adfq = adfmax_quantile(x, d, v$c, st, v$g0, x0 = v$par, tol = 1e-10),
    env = power_envelope(x, d, alpha = v$par, engine = "general"),
    ptcdf = pt_cdf(x, d, v$par, v$c),
    ptpdf = pt_pdf(x, d, v$par, v$c),
    ptcdfN = pt_cdf(x, d, v$par, v$c, N = v$c0),
    ptcdfG = pt_cdf(x, d, v$par, v$c, method = "grid"),
    ptq = pt_quantile(x, d, v$par),
    ptpow = { p <- pt_power(v$c, d, v$par, level = x); c(attr(p, "critical"), p) },
    ptenv = pt_envelope(x, d, level = v$par),
    lmcdf = lm_cdf(x, d, v$c),
    lmq = lm_quantile(x, d, v$c),
    lmpow = { p <- lm_power(x, d, level = v$par); c(attr(p, "critical"), p) },
    mzacdf = mza_cdf(x, d, v$c, gls = v$par),
    mzaq = mza_quantile(x, d, v$c, gls = v$par),
    msbcdf = msb_cdf(x, d, v$c, gls = v$par),
    msbq = msb_quantile(x, d, v$c, gls = v$par),
    stop("unknown function ", v$fn))
}

# a line 'route func engine spec a b g0 gls c c0 par n x... nv v... tol [brittle]' (spec 'case': a the model, b the
# start; 'design': a the design, b the start)
det_parse_route <- function(l) {
  w <- strsplit(l, " ", fixed = TRUE)[[1]]
  num <- function(i) as.numeric(w[i])
  n <- as.integer(w[12])
  nv <- as.integer(w[13 + n])
  list(fn = w[2], engine = c("auto", "classical", "general")[as.integer(w[3]) + 1L], spec = w[4], a = as.integer(w[5]),
       st = as.integer(w[6]), g0 = num(7), gls = num(8), c = num(9), c0 = num(10), par = num(11),
       x = num(12 + seq_len(n)), want = num(13 + n + seq_len(nv)), tol = num(14 + n + nv),
       brittle = length(w) > 14 + n + nv)
}
# a line 'rerr func engine spec a b g0 gls c c0 par n x... | text': the arguments of 'route', the call must fail with
# an error whose message contains text
det_parse_rerr <- function(l) {
  parts <- strsplit(l, " | ", fixed = TRUE)[[1]]
  w <- strsplit(parts[1], " ", fixed = TRUE)[[1]]
  num <- function(i) as.numeric(w[i])
  n <- as.integer(w[12])
  list(fn = w[2], engine = c("auto", "classical", "general")[as.integer(w[3]) + 1L], spec = w[4], a = as.integer(w[5]),
       st = as.integer(w[6]), g0 = num(7), gls = num(8), c = num(9), c0 = num(10), par = num(11),
       x = num(12 + seq_len(n)), text = parts[2])
}
det_case <- function(m, st) if (m == 4L) "DG" else paste0(if (st == 2L) "S" else "D", c("K", "U", "T")[m])

# the value of a line 'route' by the functions of the classical engine with the argument engine
det_eval_route <- function(r, D) {
  e <- r$engine
  if (r$spec == "design") {
    cas <- D[[r$a]]
    st <- det_start(r$st)
    gl <- det_gls(r$gls)
  } else {
    # (a case code with gls: the argument gls with the case code, which the R functions refuse)
    cas <- det_case(r$a, r$st)
    st <- NULL
    gl <- if (r$a != 4L && !is.nan(r$gls)) r$gls else NULL
  }
  x <- r$x
  # ('cdf@3': the method of the t-ratio, the codes of the C router: 1 line, 2 conditional, 3 contour, 4 paths)
  fm <- strsplit(r$fn, "@", fixed = TRUE)[[1]]
  me <- if (length(fm) > 1L) c("line", "conditional", "contour", "paths")[as.integer(fm[2])] else "auto"
  switch(fm[1],
    cdf = ltratio_cdf(x, cas, r$c, r$c0, r$g0, method = me, engine = e, start = st, gls = gl),
    pdf = ltratio_pdf(x, cas, r$c, r$c0, r$g0, method = me, engine = e, start = st, gls = gl),
    icdf = ltratio_icdf(x, cas, r$c, r$c0, r$g0, engine = e),
    mom = ltratio_mom(as.integer(x), cas, r$c, r$c0, r$g0, engine = e),
    coef = coef_cdf(x, cas, r$c, r$c0, r$g0, engine = e),
    coefq = coef_icdf(x, cas, r$c, r$c0, r$g0, engine = e),
    lin = ltratio_lincdf(x, r$par, cas, r$c, cbar = if (r$a == 4L) r$gls else -13.5, engine = e),
    env = power_envelope(x, cas, alpha = r$par, engine = e),
    adf = adfmax_cdf(x, cas, r$c, engine = e, start = st, g0 = r$g0),
    adfpdf = adfmax_pdf(x, cas, r$c, engine = e, start = st, g0 = r$g0),
    region1 = ltratio_region(list(type = "line", a = r$par, b = x[1]), cas, r$c, engine = e, start = st, g0 = r$g0),
    region2 = ltratio_region(list(type = "parabola", x = x[1]), cas, r$c, engine = e, start = st, g0 = r$g0),
    stop("unknown function ", r$fn))
}

# the deviations from the reference values of a line (the moments relative to max(1, |value|), as in the C test)
det_dev <- function(got, want, fn) {
  got <- as.numeric(got)
  d <- abs(got - want)
  d[is.nan(want) & is.nan(got)] <- 0
  d[is.infinite(want) & got == want] <- 0
  if (fn == "mom") d <- d / pmax(1, abs(want))
  d
}
