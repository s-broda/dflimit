# The extensions of the limit (volatility and persistence regimes, a normal component, the first seasonal t-ratio)
# against the values of the Julia engine DFLimit.jl (ref/ext_ref.txt, written by dev/ext_ref.jl of DFLimit.jl; the
# second seasonal t-ratio and the F statistic, the lines t4cdf, t4pdf, fcdf, sdval, are not ported), the argument errors,
# and the limiting cases. On CRAN the lines of ext_cran (about five seconds on two threads); otherwise every line.

ext_cran <- c("vval cdf 2 1 2 0 0", "vval pdf 2 1 2 0 0", "vquant 2 1 2 0 0 0 0.05", "cval cdf 2 0 0 0 0 0.7071",
              "cquant 2 0 0 0 0 0.7071", "sval t3cdf 2 0 0 0", "sval t3pdf 2 0 0 0", "squant t3 2 0 0.05")

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
    }
  }
  n
}

test_that("the extensions against the Julia engine (the CRAN subset)", {
  expect_equal(ext_check(TRUE), length(ext_cran))
})

test_that("the extensions against the Julia engine (every line)", {
  skip_on_cran()
  expect_gt(ext_check(FALSE), 50)
})

test_that("the argument errors of the extensions", {
  V2 <- dfl_volpath(0.3, c(1, 3))
  T0 <- dfl_trend(0); T1 <- dfl_trend(1)
  expect_error(dfl_volpath(1.2), "must lie in \\(0, 1\\)")
  expect_error(dfl_volpath(c(0.5, 0.3)), "strictly increasing")
  expect_error(dfl_volpath(0.5, 1), "one standard deviation omega per regime")
  expect_error(dfl_volpath(0.5, c(1, -1)), "finite and positive")
  expect_error(dist_cdf(-2, T0, 2, vol = V2), "must be at most 0")
  expect_error(dist_cdf(-2, T0, -5, vol = V2, start = "stationary"), "zero start only")
  expect_error(dist_cdf(-2, T1, -5, vol = V2, gls = -13.5), "GLS detrending")
  expect_error(dist_cdf(-2, T0, c(0, -1, 0), vol = V2), "one value per regime")
  expect_error(dist_cdf(-2, T0, c(0, -1)), "needs a volatility path")
  expect_error(dist_cdf(-2, T0, -5, vol = V2, rho = 0.5), "rho and vol together")
  expect_error(dist_cdf(-2, T0, -5, rho = 1.5), "rho must lie in \\[-1, 1\\]")
  expect_error(dist_cdf(-2, T1, -5, rho = 0.5, gls = -13.5), "GLS detrending")
  expect_error(dist_cdf(-2, T0, -5, rho = 0.5, gamma0 = 1, method = "lines"), "zero or a stationary start")
  expect_error(dist_cdf(-2, T0, -5, vol = V2, method = "line"), "method must be 'auto' or 'paths'")
  expect_error(dist_cdf(-2, T0, -5, vol = list(1)), "must be a volatility path")
  expect_error(cadf_c0(-5, 0), "not be 0")
  expect_error(predictive_c0(-5, 1.5), "not be 0")
  expect_error(hegy_t3_cdf(1, T0, 1), "must be at most 0")
  expect_error(hegy_t3_quantile(1.5, T0), "p must lie in")
  expect_error(ltratio_cdf(-2, "DU", 0, vol = V2, engine = "classical"), "engine classical")
  expect_error(ltratio_cdf(-2, "DU", c(0, -1)), "needs a volatility path")
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
  q <- ltratio_icdf(0.05, "DU", 0, vol = V2)
  expect_lt(abs(ltratio_cdf(q, "DU", 0, vol = V2) - 0.05), 1e-9)
  expect_lt(abs(ltratio_cdf(1.645, "DU", 0, c0 = predictive_c0(0, -0.95), rho = -0.95) - 0.5788), 1e-4)
  # non-finite points
  expect_equal(dist_cdf(c(-Inf, Inf), T0, -5, vol = V2), c(0, 1))
  expect_true(is.nan(hegy_t3_cdf(NaN, T0)))
})
