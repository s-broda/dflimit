# The extensions of the limit (volatility and persistence regimes, with their variants GLS detrending, the statistic
# that detrends first and the stationary start of the first regime; a normal component; the first seasonal t-ratio)
# against the values of the Julia engine DFLimit.jl (ref/ext_ref.txt, written by dev/ext_ref.jl of DFLimit.jl; the
# second seasonal t-ratio and the F statistic are not part of the release 0.5.0), the argument errors,
# and the limiting cases. On CRAN the lines of ext_cran (about five seconds on two threads); otherwise every line.
# The tags vdval, vroots, vtstar, rdval, rroots and rtstar are layers of the engine (the transform, its roots, t_*)
# that the package does not reach: they are skipped.

ext_cran <- c("vval cdf 2 1 2 0 0", "vval pdf 2 1 2 0 0", "vquant 2 1 2 0 0 0 0.05", "cval cdf 2 0 0 0 0 0.7071",
              "cquant 2 0 0 0 0 0.7071", "sval t3cdf 2 0 0 0", "sval t3pdf 2 0 0 0", "squant t3 2 0 0.05",
              # (the regime variants: GLS detrending, the statistic that detrends first, the stationary start)
              "rval cdf 2 5 2 0 0 0 -7 0", "rval cdf 8 9 2 0 0 0 NaN 1", "rcase cdf SU 1 2 -4 -4")

ext_lines <- function(cran) {
  l <- readLines(test_path("ref", "ext_ref.txt"))
  l <- l[!startsWith(l, "#")]
  if (cran) {
    keep <- startsWith(l, "design ") | startsWith(l, "vol ") |
      vapply(l, function(x) any(startsWith(gsub("(\\.[0-9]{4})[0-9]*", "\\1", x), ext_cran)), logical(1))
    l <- l[keep]
  }
  l
}
ext_check <- function(cran) {
  L <- ext_lines(cran)
  D <- list(); V <- list(); n <- 0L
  num <- function(w) as.numeric(w)
  # a regime variant from w[k] on, <var> = did vid nc c_1..c_nc start gls first (start 0 zero, 2 stationary; gls NaN:
  # least squares; first 1: the statistic that detrends first): list(D, v, c, start, gls, k the position after it)
  rvar <- function(w, k) {
    nc <- as.integer(w[k + 2]); j <- k + 3 + nc
    d <- D[[as.integer(w[k])]]; g <- num(w[j + 1])
    list(D = if (w[j + 2] == "1") dfl_detrend_first(d) else d, v = V[[as.integer(w[k + 1])]], c = num(w[k + 2 + seq_len(nc)]),
         start = if (w[j] == "2") "stationary" else "zero", gls = if (is.nan(g)) NULL else g, k = j + 3)
  }
  for (l in L) {
    w <- strsplit(l, " ", fixed = TRUE)[[1]]
    tag <- w[1]
    if (tag == "design") {
      id <- as.integer(w[2])
      D[[id]] <- if (w[3] == "trend") dfl_trend(as.integer(w[4])) else dfl_perron(w[4], num(w[5 + seq_len(as.integer(w[5]))]))
    } else if (tag == "vol") {
      nb <- as.integer(w[3])
      V[[as.integer(w[2])]] <- dfl_volpath(num(w[3 + seq_len(nb)]), num(w[3 + nb + seq_len(nb + 1)]))
    } else if (tag %in% c("vval", "vquant")) {
      k0 <- if (tag == "vval") 3L else 2L
      d <- D[[as.integer(w[k0])]]; v <- V[[as.integer(w[k0 + 1])]]; nc <- as.integer(w[k0 + 2])
      cc <- num(w[k0 + 2 + seq_len(nc)]); r <- num(w[(k0 + 3 + nc):length(w)])
      if (tag == "vval") {
        nx <- r[2]; x <- r[2 + seq_len(nx)]; ref <- r[3 + nx + seq_len(nx)]; tol <- r[length(r)]
        got <- if (w[2] == "cdf") dist_cdf(x, d, cc, c0 = r[1], vol = v) else dist_pdf(x, d, cc, c0 = r[1], vol = v)
      } else {
        got <- dist_quantile(r[2], d, cc, c0 = r[1], vol = v); ref <- r[3]; tol <- r[4]
      }
      expect_lt(max(abs(got - ref)), tol, label = l)
      n <- n + 1L
    } else if (tag == "cval") {
      d <- D[[as.integer(w[3])]]; st <- c("zero", "fixed", "stationary")[as.integer(w[4]) + 1L]
      g0 <- num(w[5]); cc <- num(w[6]); c0 <- num(w[7]); rho <- num(w[8]); m <- if (rho == 0) "auto" else w[9]
      r <- num(w[10:length(w)]); nx <- r[1]; x <- r[1 + seq_len(nx)]; ref <- r[2 + nx + seq_len(nx)]; tol <- r[length(r)]
      f <- if (w[2] == "cdf") dist_cdf else dist_pdf
      got <- f(x, d, cc, c0 = c0, start = st, gamma0 = g0, rho = rho, method = m)
      expect_lt(max(abs(got - ref)), tol, label = l)
      n <- n + 1L
    } else if (tag == "cquant") {
      d <- D[[as.integer(w[2])]]; r <- num(w[5:10])      # (c, c0, rho, p, x, tol)
      expect_lt(abs(dist_quantile(r[4], d, r[1], c0 = r[2], rho = r[3]) - r[5]), r[6], label = l)
      n <- n + 1L
    } else if (tag == "sval" && startsWith(w[2], "t3")) {
      d <- D[[as.integer(w[3])]]; cc <- num(w[4]); c0 <- num(w[5]); vid <- as.integer(w[6])
      r <- num(w[7:length(w)]); nx <- r[1]; x <- r[1 + seq_len(nx)]; ref <- r[2 + nx + seq_len(nx)]; tol <- r[length(r)]
      v <- if (vid > 0L) V[[vid]] else NULL; c2 <- if (vid > 0L) c(cc, cc) else cc
      got <- if (w[2] == "t3cdf") hegy_t3_cdf(x, d, c2, c0 = c0, vol = v) else hegy_t3_pdf(x, d, c2, c0 = c0, vol = v)
      expect_lt(max(abs(got - ref)), tol, label = l)
      n <- n + 1L
    } else if (tag == "squant" && w[2] == "t3") {
      d <- D[[as.integer(w[3])]]; r <- num(w[4:7])
      expect_lt(abs(hegy_t3_quantile(r[2], d, r[1]) - r[3]), r[4], label = l)
      n <- n + 1L
    } else if (tag == "rval") {
      a <- rvar(w, 3L); r <- num(w[a$k:length(w)])     # (c0, n, x.., n, v.., tol)
      nx <- r[2]; x <- r[2 + seq_len(nx)]; ref <- r[3 + nx + seq_len(nx)]; tol <- r[length(r)]
      f <- if (w[2] == "cdf") dist_cdf else dist_pdf
      got <- f(x, a$D, a$c, c0 = r[1], start = a$start, gls = a$gls, vol = a$v)
      expect_lt(max(abs(got - ref)), tol, label = l)
      n <- n + 1L
    } else if (tag == "rquant") {
      a <- rvar(w, 2L); r <- num(w[a$k:length(w)])     # (c0, p, x, tol)
      got <- dist_quantile(r[2], a$D, a$c, c0 = r[1], start = a$start, gls = a$gls, vol = a$v)
      expect_lt(abs(got - r[3]), r[4], label = l)
      n <- n + 1L
    } else if (tag == "rcase") {
      v <- V[[as.integer(w[4])]]; nc <- as.integer(w[5]); cc <- num(w[5 + seq_len(nc)])
      r <- num(w[(6 + nc):length(w)]); nx <- r[2]; x <- r[2 + seq_len(nx)]; ref <- r[3 + nx + seq_len(nx)]; tol <- r[length(r)]
      got <- switch(w[2], cdf = ltratio_cdf(x, w[3], cc, c0 = r[1], vol = v), pdf = ltratio_pdf(x, w[3], cc, c0 = r[1], vol = v),
                    icdf = ltratio_icdf(x, w[3], cc, c0 = r[1], vol = v))
      expect_lt(max(abs(got - ref)), tol, label = l)
      n <- n + 1L
    } else if (tag == "rwarn") {
      p <- strsplit(l, " | ", fixed = TRUE)[[1]]; w <- strsplit(p[1], " ", fixed = TRUE)[[1]]
      a <- rvar(w, 2L); r <- num(w[a$k:length(w)])     # (c0, x, v, tol)
      expect_warning(got <- dist_cdf(r[2], a$D, a$c, c0 = r[1], start = a$start, gls = a$gls, vol = a$v), p[2], fixed = TRUE)
      expect_lt(abs(got - r[3]), r[4], label = l)
      n <- n + 1L
    }
  }
  n
}

test_that("the extensions against the Julia engine (the CRAN subset)", {
  expect_equal(ext_check(TRUE), length(ext_cran))
})

test_that("the extensions against the Julia engine (every line)", {
  skip_on_cran()
  expect_gt(ext_check(FALSE), 120)
})

test_that("the argument errors of the extensions", {
  # (the err lines of ref/ext_ref.txt, in its order, then the R functions' own errors)
  expect_equal(sum(startsWith(ext_lines(FALSE), "err ")), 18)
  V2 <- dfl_volpath(0.3, c(1, 3))
  T0 <- dfl_trend(0); T1 <- dfl_trend(1); F1 <- dfl_detrend_first(T1)
  expect_error(dfl_volpath(1.2), "must lie in \\(0, 1\\)")
  expect_error(dfl_volpath(c(0.5, 0.3)), "strictly increasing")
  expect_error(dfl_volpath(0.5, 1), "one standard deviation omega per regime")
  expect_error(dfl_volpath(0.5, c(1, -1)), "finite and positive")
  expect_error(dist_cdf(-2, T0, 2, vol = V2), "must be at most 0")
  expect_error(dist_cdf(-2, T0, -5, vol = V2, start = "fixed", gamma0 = 1), "fixed start (gamma0 != 0) is not computed", fixed = TRUE)
  expect_error(dist_cdf(-2, T0, c(0, -5), vol = V2, start = "stationary"), "needs c < 0 on the first regime", fixed = TRUE)
  expect_error(dist_cdf(-2, F1, -5, vol = V2, gls = -13.5), "detrends first already")
  expect_error(dist_cdf(-2, T1, -5, vol = V2, gls = -13.5, start = "stationary"), "zero start only")
  expect_error(dist_cdf(-2, F1, -5, vol = V2, start = "stationary"), "zero start only")
  expect_error(dist_cdf(-2, T0, c(0, -1, 0), vol = V2), "one value per regime")
  expect_error(dist_cdf(-2, T0, c(0, -1)), "needs a volatility path")
  expect_error(dist_cdf(-2, T0, -5, vol = V2, rho = 0.5), "rho and vol together")
  expect_error(dist_cdf(-2, T0, -5, rho = 1.5), "rho must lie in \\[-1, 1\\]")
  expect_error(dist_cdf(-2, T1, -5, rho = 0.5, gls = -13.5), "GLS detrending")
  expect_error(dist_cdf(-2, T0, -5, rho = 0.5, gamma0 = 1, method = "lines"), "zero or a stationary start")
  expect_error(dist_cdf(-2, T0, -5, vol = V2, method = "line"), "method must be 'auto' or 'paths'")
  # (the methods of rho without it: the refusal and the message of DFLimit.jl, November 26, 2026)
  lmsg <- "method must be 'auto', 'line' (the line integral with the explicit kernel) or 'paths' (the inversion on paths around the roots)"
  for (m in c("lines", "mixture")) {
    expect_error(dist_cdf(-2, T0, -5, method = m), lmsg, fixed = TRUE)
    expect_error(dist_pdf(-2, T0, -5, method = m), lmsg, fixed = TRUE)
    expect_error(dist_quantile(0.05, T0, -5, method = m), lmsg, fixed = TRUE)
  }
  expect_error(dist_cdf(-2, T0, -5, vol = list(1)), "must be a volatility path")
  expect_error(cadf_c0(-5, 0), "not be 0")
  expect_error(predictive_c0(-5, 1.5), "not be 0")
  expect_error(hegy_t3_cdf(1, T0, 1), "must be at most 0")
  expect_error(hegy_t3_quantile(1.5, T0), "p must lie in")
  expect_error(ltratio_cdf(-2, "DU", 0, vol = V2, engine = "classical"), "engine classical")
  expect_error(ltratio_cdf(-2, "DU", c(0, -1)), "needs a volatility path")
  # (the refusals of the regime variants through the case codes and the quantiles)
  expect_error(ltratio_cdf(-2, "DU", -5, g0 = 1, vol = V2), "fixed start (gamma0 != 0) is not computed", fixed = TRUE)
  expect_error(ltratio_icdf(0.05, "SU", c(0, -5), vol = V2), "needs c < 0 on the first regime", fixed = TRUE)
  expect_error(dist_quantile(0.05, F1, -5, vol = V2, gls = -13.5), "detrends first already")
})

test_that("the limiting cases of the extensions", {
  xs <- c(-4, -2.86, -1, 0.5)
  T0 <- dfl_trend(0)
  # one regime is the engine itself; a dummy regime boundary at a constant variance is the engine's paths to rounding
  expect_identical(dist_cdf(xs, T0, -10, vol = dfl_volpath()), dist_cdf(xs, T0, -10))
  expect_identical(dist_quantile(0.05, dfl_trend(1), 0, vol = dfl_volpath()), dist_quantile(0.05, dfl_trend(1), 0))
  for (cc in c(0, -10)) {
    expect_lt(max(abs(dist_cdf(xs, T0, cc, vol = dfl_volpath(0.4)) - dist_cdf(xs, T0, cc, method = "paths"))), 1e-12)
    expect_lt(max(abs(dist_pdf(xs, T0, cc, vol = dfl_volpath(0.4)) - dist_pdf(xs, T0, cc, method = "paths"))), 1e-12)
  }
  # the regime variants at a dummy boundary, constant variance and persistence: the engine with the same detrending,
  # design and start (GLS detrending, the statistic that detrends first, the stationary start; about 1e-15)
  V1 <- dfl_volpath(0.4, c(1, 1)); T1 <- dfl_trend(1); F1 <- dfl_detrend_first(T1)
  for (cc in c(0, -10)) {
    expect_lt(max(abs(dist_cdf(xs, T1, c(cc, cc), gls = -13.5, vol = V1) - dist_cdf(xs, T1, cc, gls = -13.5, method = "paths"))), 1e-10)
    expect_lt(max(abs(dist_pdf(xs, T1, c(cc, cc), gls = -13.5, vol = V1) - dist_pdf(xs, T1, cc, gls = -13.5, method = "paths"))), 1e-10)
    expect_lt(max(abs(dist_cdf(xs, F1, c(cc, cc), vol = V1) - dist_cdf(xs, F1, cc, method = "paths"))), 1e-10)
    expect_lt(max(abs(dist_pdf(xs, F1, c(cc, cc), vol = V1) - dist_pdf(xs, F1, cc, method = "paths"))), 1e-10)
  }
  for (cc in c(-5, -10)) {
    expect_lt(max(abs(dist_cdf(xs, T0, c(cc, cc), start = "stationary", vol = V1) - dist_cdf(xs, T0, cc, start = "stationary", method = "paths"))), 1e-10)
    expect_lt(max(abs(dist_pdf(xs, T0, c(cc, cc), start = "stationary", vol = V1) - dist_pdf(xs, T0, cc, start = "stationary", method = "paths"))), 1e-10)
  }
  # the endpoints and the reflection of rho
  x4 <- c(-3, -1, 0.5, 2)
  expect_lt(max(abs(dist_cdf(x4, T0, -5, rho = 0) - pnorm(x4))), 1e-15)
  expect_lt(max(abs(dist_cdf(x4, T0, -5, rho = -1) - (1 - dist_cdf(-x4, T0, -5)))), 1e-15)
  expect_lt(max(abs(dist_cdf(x4, T0, -5, rho = -0.7, method = "lines") - (1 - dist_cdf(-x4, T0, -5, rho = 0.7, method = "lines")))), 1e-15)
  expect_equal(cadf_c0(-5, sqrt(0.5)), 5)
  expect_equal(predictive_c0(-5, -0.95, 1), -5 + 1 / 0.95)
  # the classical functions take vol and rho (the general engine for the case's design)
  V2 <- dfl_volpath(0.3, c(1, 3))
  expect_equal(ltratio_cdf(c(-3, -2), "DU", -5, vol = V2), dist_cdf(c(-3, -2), T0, -5, vol = V2), tolerance = 1e-15)
  # (case SU: the intercept with a stationary start; case DG: the linear trend GLS-detrended at -13.5)
  x3 <- c(-3.5, -2.5, -1)
  expect_equal(ltratio_cdf(x3, "SU", c(-4, -4), vol = V2), dist_cdf(x3, T0, c(-4, -4), start = "stationary", vol = V2), tolerance = 1e-15)
  expect_equal(ltratio_pdf(x3, "DG", c(0, -5), vol = V2), dist_pdf(x3, T1, c(0, -5), gls = -13.5, vol = V2), tolerance = 1e-15)
  q <- ltratio_icdf(0.05, "DU", 0, vol = V2)
  expect_lt(abs(ltratio_cdf(q, "DU", 0, vol = V2) - 0.05), 1e-9)
  expect_lt(abs(ltratio_cdf(1.645, "DU", 0, c0 = predictive_c0(0, -0.95), rho = -0.95) - 0.5788), 1e-4)
  # non-finite points
  expect_equal(dist_cdf(c(-Inf, Inf), T0, -5, vol = V2), c(0, 1))
  expect_true(is.nan(hegy_t3_cdf(NaN, T0)))
})

test_that("the outer contour under regimes: a root near the imaginary axis far out", {
  # (the testset of test/extensions/regimes.jl of DFLimit.jl: cases in which the thin step on the engine's 45-degree
  # hyperbola lost 1.4e-6 to 3.8e-4 at every x, and in the far upper tail the range check of dist_cdf() refused the
  # value: an early regime of high variance with a break in the deterministic terms or a persistence that changes; the
  # hyperbola under regimes starts at 70 degrees since the C engine of dflimit-matlab ffeb48e)
  skip_on_cran()
  xs <- c(-4, -2, -1, 0, 0.5, 1, 1.5, 2, 3, 5)
  V3 <- dfl_volpath(c(0.2, 0.7), c(2, 1, 0.5))
  cases <- list(list(D = dfl_perron("B", 0.5), v = V3, c = c(-5, 0, -10), start = NULL, c0 = 0),
                list(D = dfl_perron("C", 0.4), v = dfl_volpath(0.3, c(3, 1)), c = c(-5, -5), start = "stationary", c0 = 0),
                list(D = dfl_perron("B", 0.3), v = dfl_volpath(0.2, c(5, 1)), c = c(-5, -20), start = NULL, c0 = 0),
                list(D = dfl_perron("B", 0.5), v = V3, c = c(-20, 0, -10), start = NULL, c0 = 1))
  for (a in cases) {
    F <- dist_cdf(xs, a$D, a$c, c0 = a$c0, start = a$start, vol = a$v)
    expect_lt(max(abs(F - dist_cdf(xs, a$D, a$c, c0 = a$c0, start = a$start, verified = TRUE, vol = a$v))), 1e-10)
    expect_true(all(diff(F) > -1e-11) && 1 - F[length(F)] < 1e-11)
    expect_gt(min(dist_pdf(xs, a$D, a$c, c0 = a$c0, start = a$start, vol = a$v)), -1e-11)
  }
})
