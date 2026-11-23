# The devices of Bao and Kan (the single inversion of the half-plane statistics, the moments from jets, the single inversion
# of the general engine's point-optimal statistic and coef_cdf_gp) against the values of the Julia engine DFLimit.jl
# (ref/kan_ref.txt, the user-level lines of c/test/kan_ref.txt of dflimit-matlab), with their options, error estimates,
# warnings and refusals. The method is that of Yong Bao and Raymond Kan, "Exact Finite-Sample Distribution and Moments of
# the Least Squares Estimator in First-Order Autoregressions: O(1) Algorithms".

kan_lines <- function() {
  l <- readLines(test_path("ref", "kan_ref.txt"))
  l[!startsWith(l, "#") & nzchar(l)]
}
kan_design <- function(v) {
  # (v: kind np par..., then the arguments; kind 0 the trend of order par, 1, 2, 3 Perron's models A, B, C)
  kind <- v[1]; np <- v[2]; par <- v[2 + seq_len(np)]
  D <- if (kind == 0) dfl_trend(as.integer(par[1])) else dfl_perron(c("A", "B", "C")[kind], par)
  list(D = D, rest = v[-seq_len(2 + np)])
}
kan_start <- function(k) c("zero", "fixed", "stationary")[k + 1]
# the tolerances of the C engine's test (dfl_test_kan.c): absolute for probabilities, relative to max(1, |value|) for
# densities and quantiles, relative to the size of the moments (max(|m_j|, |m_2|^(j/2)))
kan_tol <- c(rcoef = 1e-11, rcoefpdf = 1e-11, rcoefq = 1e-8, rcoefmom = 1e-10, rmom = 1e-10, rlin = 1e-11, renv = 1e-11,
             ptcdf = 1e-11, ptpdf = 1e-11, ptq = 1e-10, ptpow = 1e-10, ptenv = 1e-11, gp = 1e-11)

kan_dev <- function(tag, w) {
  v <- suppressWarnings(as.numeric(w[-1]))
  if (tag %in% c("ptcdf", "ptpdf", "ptq", "ptpow", "ptenv", "gp")) {
    d <- kan_design(v[-1])
    D <- d$D; a <- d$rest
  }
  rel <- function(g, r) abs(g - r) / pmax(1, abs(r))
  same <- function(g, r) ifelse(is.nan(g) & is.nan(r), 0, abs(g - r))
  switch(tag,
    rcoef = , rcoefpdf = , rcoefq = , rcoefmom = , rmom = {
      cas <- case_code(v[1], v[2] == 2)
      g0 <- if (v[2] == 1) v[3] else 0
      c <- v[4]; c0 <- v[5]; n <- v[6]
      if (tag %in% c("rcoefmom", "rmom")) {
        r <- v[6 + seq_len(n)]
        g <- if (tag == "rmom") ltratio_mom(n, cas, c, c0, g0) else coef_mom(n, cas, c, c0, g0)
        return(abs(g - r) / pmax(abs(r), abs(r[min(2, n)])^(seq_len(n) / 2)))
      }
      x <- v[6 + seq_len(n)]; r <- v[6 + n + seq_len(n)]
      switch(tag,
        rcoef = same(coef_cdf(x, cas, c, c0, g0), r),
        rcoefpdf = rel(coef_pdf(x, cas, c, c0, gamma0 = g0), r),
        rcoefq = rel(coef_icdf(x, cas, c, c0, g0, tol = 1e-9), r))
    },
    rlin = {
      n <- v[6]; b <- v[6 + seq_len(n)]; r <- v[6 + n + seq_len(n)]
      g <- ltratio_lincdf(b, v[5], dfl_trend(as.integer(v[1])), v[4], start = kan_start(v[2]), g0 = if (v[2] == 1) v[3] else 0)
      abs(g - r)
    },
    renv = {
      n <- v[2]; c <- v[2 + seq_len(n)]; r <- v[2 + n + seq_len(n)]
      same(power_envelope(c, "DK", alpha = v[1]), r)
    },
    ptcdf = , ptpdf = {
      n <- a[3]; q <- a[3 + seq_len(n)]; r <- a[3 + n + seq_len(n)]
      if (tag == "ptcdf") abs(pt_cdf(q, D, a[1], a[2]) - r) else rel(pt_pdf(q, D, a[1], a[2]), r)
    },
    ptq = rel(pt_quantile(a[2], D, a[1]), a[3]),
    ptpow = {
      P <- pt_power(a[2], D, a[1], level = a[3])
      max(rel(attr(P, "critical"), a[4]), abs(as.numeric(P) - a[5]))
    },
    ptenv = {
      n <- a[2]; c <- a[2 + seq_len(n)]; r <- a[2 + n + seq_len(n)]
      same(pt_envelope(c, D, level = a[1]), r)
    },
    gp = {
      n <- a[5]; y <- a[5 + seq_len(n)]; r <- a[5 + n + seq_len(n)]
      same(coef_cdf_gp(y, D, a[3], a[4], start = kan_start(a[1]), gls = if (is.nan(a[2])) NULL else a[2]), r)
    })
}

test_that("the devices of Bao and Kan reproduce the Julia values", {
  skip_on_cran()
  L <- kan_lines()
  n <- 0L
  for (l in L) {
    w <- strsplit(l, " ", fixed = TRUE)[[1]]
    d <- suppressWarnings(kan_dev(w[1], w))
    expect(all(d <= kan_tol[[w[1]]]), sprintf("%s: deviation %.3g above the tolerance %.0g", substr(l, 1, 100), max(d), kan_tol[[w[1]]]))
    n <- n + 1L
  }
  expect_equal(n, length(L))
})

test_that("the defaults of the half-plane statistics (a few values for CRAN)", {
  # (ref/kan_ref.txt: coef_cdf at c = -5 in case DK, ltratio_lincdf of the point-optimal test, the envelope at c = -7)
  expect_equal(coef_cdf(c(-9, -3), "DK", -5), c(0.25024110084038587, 0.85251049438730919), tolerance = 1e-11)
  expect_equal(ltratio_lincdf(c(-0.2639, 0.5), -3.5, "DK"), c(0.049996781931337997, 0.41434224345768855), tolerance = 1e-11)
  expect_equal(power_envelope(-7), 0.49940198281556841, tolerance = 1e-11)
  expect_equal(ltratio_mom(4, "DK"), c(-0.42309564509731556, 1.1415950751360362, -1.0617121847180102, 3.5238029768231982),
               tolerance = 1e-10)
  expect_equal(coef_mom(1, "DK"), -1.7814, tolerance = 1e-4)
})

test_that("method = \"engine\" reproduces the engines' values; the defaults agree with them", {
  y <- c(-12, -5, 0.5)
  Fe <- coef_cdf(y, "DU", -5, method = "engine")
  expect_identical(Fe, vapply(y, function(x) ltratio_cdf(0, "DU", -5, c0 = x), 0))
  expect_lt(max(abs(coef_cdf(y, "DU", -5) - Fe)), 1e-9)
  expect_identical(coef_cdf(y, "DU", -5, engine = "classical"), coef_cdf(y, "DU", -5, method = "engine", engine = "classical"))
  b <- c(-0.2639, 0.3)
  Pe <- ltratio_lincdf(b, -3.5, "DK", -7, method = "engine")
  expect_identical(Pe, ltratio_lincdf(b, -3.5, "DK", -7, engine = "classical"))
  expect_lt(max(abs(ltratio_lincdf(b, -3.5, "DK", -7) - Pe)), 1e-9)
  Ee <- power_envelope(c(-5, -13.5), method = "engine")
  expect_identical(Ee, power_envelope(c(-5, -13.5), engine = "classical"))
  expect_lt(max(abs(power_envelope(c(-5, -13.5)) - Ee)), 1e-9)
  me <- ltratio_mom(3, "DT", -5, method = "engine")
  expect_equal(me, ltratio_mom(3, "DT", -5, engine = "classical"), tolerance = 1e-12)
  expect_lt(max(abs(ltratio_mom(3, "DT", -5) - me) / abs(me)), 1e-9)
  fe <- coef_pdf(c(-10, -3), "DU", -5, method = "engine")
  expect_identical(fe, coef_pdf(c(-10, -3), dfl_trend(0), -5, engine = "general"))
  expect_lt(max(abs(coef_pdf(c(-10, -3), "DU", -5) - fe)), 1e-9)
  # (the point-optimal statistic and coef_cdf_gp: the fixed grid of before)
  D <- dfl_trend(1)
  expect_lt(max(abs(pt_cdf(c(4, 7), D, -13.5, -5) - pt_cdf(c(4, 7), D, -13.5, -5, method = "grid"))), 1e-8)
  expect_identical(pt_cdf(c(4, 7), D, -13.5, -5, method = "grid"), pt_cdf(c(4, 7), D, -13.5, -5, N = 60000, umax = 4000))
  expect_lt(abs(coef_cdf_gp(-12, dfl_trend(0), -5) - coef_cdf_gp(-12, dfl_trend(0), -5, method = "grid")), 1e-9)
  expect_lt(abs(coef_cdf_gp(-12, dfl_trend(0), -5) - coef_cdf(-12, "DU", -5)), 1e-11)
})

test_that("the error estimates, the warnings and verified", {
  r <- coef_cdf(c(-9, -3), "DK", -5, err = TRUE)
  expect_named(r, c("value", "err"))
  expect_true(all(r$err >= 0 & r$err <= 1e-12))
  expect_true(all(is.nan(coef_cdf(c(-9, -3), "DK", -5, method = "engine", err = TRUE)$err)))
  expect_true(all(is.nan(coef_cdf(-9, "DG", -5, err = TRUE)$err)))
  # (a tolerance that cannot be met: the warning names the points; verified recomputes them by the engines)
  expect_warning(coef_cdf(-9, "DK", -5, tol = 1e-30), "exceeds")
  expect_warning(ltratio_lincdf(-0.2639, -3.5, "DK", -7, tol = 1e-30), "exceeds")
  expect_identical(suppressWarnings(coef_cdf(-9, "DK", -5, tol = 1e-30, verified = TRUE)),
                   coef_cdf(-9, "DK", -5, method = "engine", verified = TRUE))
  expect_identical(ltratio_lincdf(-0.2639, -3.5, "DK", -7, tol = 1e-30, verified = TRUE),
                   ltratio_lincdf(-0.2639, -3.5, "DK", -7, method = "engine"))
  e <- power_envelope(-7, err = TRUE)
  expect_lt(e$err, 1e-11)
  m <- ltratio_mom(2, "DU", info = TRUE)
  expect_true(all(m$err >= 0 & m$err < 1e-9))
  expect_true(all(is.na(ltratio_mom(2, "DU", info = TRUE, engine = "general")$err)))
  p <- pt_cdf(c(4, 7), dfl_trend(1), -13.5, -5, err = TRUE)
  expect_true(all(p$err <= 1e-12))
  expect_identical(pt_cdf(c(-Inf, Inf, NaN), dfl_trend(1), -13.5, err = TRUE)$err, c(0, 0, 0))
  g <- coef_cdf_gp(c(-12, Inf), dfl_trend(0), -5, err = TRUE)
  expect_true(g$err[1] <= 1e-12 && g$err[2] == 0)
  expect_true(is.nan(coef_cdf_gp(-12, dfl_trend(0), -5, method = "grid", err = TRUE)$err))
})

test_that("the rules of the defaults: the regions and the cases of the single inversion", {
  # (above c = 20 the engines; method = "inversion" forces the single inversion)
  expect_identical(coef_cdf(-2, "DU", 25), coef_cdf(-2, "DU", 25, method = "engine"))
  expect_lt(abs(coef_cdf(-2, "DU", 25, method = "inversion") - coef_cdf(-2, "DU", 25)), 1e-8)
  # (the linear form beyond c = 4 in case DK, the envelope at -1e-4 < c < 0: the engines)
  expect_identical(ltratio_lincdf(1, -1, "DK", 5), ltratio_lincdf(1, -1, "DK", 5, method = "engine"))
  expect_identical(power_envelope(-1e-5), power_envelope(-1e-5, method = "engine"))
  # (the point-optimal statistic above c = 4, coef_cdf_gp above c = 5 and with a fixed start: the fixed grid)
  # (the grid's refusal above c = 3.1 with an intercept at cbar = -7)
  D <- dfl_trend(0)
  vm <- function(e) tryCatch(e, error = conditionMessage)
  expect_identical(vm(pt_cdf(30, D, -7, 4.2)), vm(pt_cdf(30, D, -7, 4.2, method = "grid")))
  expect_match(vm(pt_cdf(30, D, -7, 4.2)), "does not resolve")
  expect_identical(vm(coef_cdf_gp(-2, D, 5.5)), vm(coef_cdf_gp(-2, D, 5.5, method = "grid")))
  expect_identical(coef_cdf_gp(-12, D, -5, start = "fixed", gamma0 = 1), coef_cdf_gp(-12, D, -5, start = "fixed", gamma0 = 1, method = "grid"))
  # (a design that is a classical case takes the single inversion; others the general engine)
  expect_identical(coef_pdf(-8, dfl_trend(-1)), coef_pdf(-8, "DK"))
  expect_identical(ltratio_lincdf(-0.26, -3.5, dfl_trend(0), -5), ltratio_lincdf(-0.26, -3.5, dfl_trend(0), -5, method = "inversion"))
  expect_identical(coef_pdf(-8, dfl_perron("A", 0.5)), coef_pdf(-8, dfl_perron("A", 0.5), method = "engine"))
})

test_that("the refusals of the single inversion", {
  expect_error(coef_cdf(-9, "DK", method = "inversion", engine = "general"), "does not apply with engine")
  expect_error(coef_cdf(-9, "DG", method = "inversion"), "applies to the cases DK, DU, DT, SK, SU, ST")
  expect_error(coef_cdf(-9, "DK", method = "fast"), "method must be")
  expect_error(coef_cdf(-9, "DK", tol = -1), "tol must be")
  expect_error(coef_pdf(-8, dfl_perron("A", 0.5), method = "inversion"), "applies to the cases")
  expect_error(coef_mom(2, "DG"), "available for the cases")
  expect_error(coef_mom(2, "DK", engine = "general"), "no moments of the coefficient statistic")
  expect_error(coef_mom(2, "DK", 25), "at most 20")
  expect_error(coef_mom(2, "DK", method = "engine"), "method must be")
  expect_error(ltratio_mom(2, "DK", method = "jets", engine = "general"), "does not apply with engine")
  expect_error(ltratio_lincdf(0, -3.5, "DG", method = "inversion"), "not to case DG")
  expect_error(ltratio_lincdf(0, NaN, "DK", method = "inversion"), "finite a")
  expect_error(power_envelope(-5, "DG", method = "inversion"), "not to case DG")
  expect_error(power_envelope(-5, method = "inversion", engine = "general"), "does not apply with engine")
  expect_error(pt_quantile(0.05, dfl_trend(1), -13.5, method = "inversion", umax = 10), "N and umax belong")
  expect_error(coef_cdf_gp(-12, dfl_trend(0), -5, start = "fixed", gamma0 = 1, method = "inversion"), "zero or stationary start")
})
