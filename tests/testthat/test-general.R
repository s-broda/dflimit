# The general engine (any deterministic design) and the router over the two engines, against the values of the Julia
# engine DFLimit.jl (ref/det_ref.txt; see helper-det.R). On CRAN a subset of the lines (named by their computation,
# the first tokens of each line, below) is checked, about four seconds on two threads; otherwise every line (about
# fifty seconds).

det_cran <- c(
  # every constructor (perron A, two and three breaks, ltrend, fdesign, design_from_matrix, detrend_first), a cubic
  "val cdf 7 NaN 0 0 0 0 0 0 5 -5", "val pdf 7 NaN 0 0 0 0 0 0 5 -5", "val cdf 13 NaN 0 0 0 0 0 0 5 -5",
  "val cdf 14 NaN 0 0 0 0 0 0 5 -5", "val cdf 16 NaN 0 0 0 0 0 0 5 -5", "val cdf 18 NaN 0 0 0 0 0 0 5 -5",
  "val cdf 21 NaN 0 0 0 0 0 0 5 -5", "val cdf 10 NaN 0 0 0 0 0 0 5 -5", "val cdf 4 NaN 0 0 0 0 0 0 6 -4",
  # the three starts, GLS detrending, a local null value (and the tangency of case DK), the checked computation
  "val cdf 1 NaN 1 1.5 -5 0 0 0 4 -4", "val pdf 1 NaN 1 1.5 -5 0 0 0 4 -4", "val cdf 3 NaN 2 0 -10 0 0 0 4 -4",
  "val cdf 13 NaN 2 0 -5 0 0 0 4 -4", "val cdf 2 -7 0 0 -5 0 0 0 5 -3.5", "val pdf 2 -7 0 0 -5 0 0 0 5 -3.5",
  "val cdf 3 NaN 0 0 -5 3 0 0 5 -4.5", "val cdf 1 NaN 0 0 -10 2 0 0 5 -4.5", "val cdf 1 NaN 0 0 -5 2 0 0 1 -2",
  "val cdf 3 NaN 0 0 -5 0 1 0 2 -3",
  # each function of the general engine and each companion
  "val quant 1 NaN 0 0 0 0 0 -2 1 0.050000000000000003", "val coefcdf 2 NaN 0 0 0 0 0 0 4 -25",
  "val coefpdf 2 NaN 0 0 0 0 0 0 4 -25", "val coefq 2 NaN 0 0 0 0 0 -10 1 0.050000000000000003",
  "val mom 1 NaN 0 0 -50 0 0 4 4 1", "val lpow 3 NaN 0 0 0 0 0 0 2 -5",
  "val lincdf 1 NaN 0 0 -5 0 0 NaN 4 0.02", "val lincdf 1 NaN 0 0 0 0 0 0 4 -0.59999999999999998",
  "val linpdf 1 NaN 0 0 -5 0 0 NaN 4 0.02", "val linpdf 1 NaN 0 0 0 0 0 0 4 -0.59999999999999998",
  "val linq 1 NaN 0 0 0 0 0 -3.5 1 0.050000000000000003", "val adfcdf 2 NaN 0 0 0 0 0 0 3 -3",
  "val adfpdf 2 NaN 0 0 0 0 0 0 3 -3", "val adfq 2 NaN 0 0 0 0 0 -2 1 0.050000000000000003",
  "val env 1 NaN 0 0 0 0 0 0.050000000000000003 2 -5", "val ptcdf 2 -7 0 0 0 0 0 -7 4 2",
  "val ptpdf 2 -7 0 0 0 0 0 -7 4 2", "val lmcdf 3 NaN 0 0 -10 0 0 0 5 -3", "val lmq 3 NaN 0 0 0 0 0 0 1 0.050000000000000003",
  "val lmpow 8 NaN 0 0 0 0 0 0.050000000000000003 2 -5", "val mzacdf 3 -13.5 0 0 0 0 0 -13.5 4 -25",
  "val msbcdf 3 -13.5 0 0 0 0 0 -13.5 4 0.080000000000000002",
  "val msbq 3 -13.5 0 0 0 0 0 -13.5 1 0.050000000000000003",
  # the router: each function with engine auto (0), classical (1) and general (2); designs
  "route cdf 0 case 4 0 0 NaN -5 0 0 5 -3.5", "route cdf 1 case 4 0 0 NaN -5 0 0 5 -3.5",
  "route cdf 2 case 1 0 0 NaN -5 2 0 5 -3.5", "route pdf 0 case 4 0 0 NaN -5 0 0 5 -3.5",
  "route pdf 1 case 4 0 0 NaN -5 0 0 5 -3.5", "route pdf 2 case 1 0 0 NaN -5 2 0 5 -3.5",
  "route cdf 0 design 3 0 0 -13.5 -5 0 0 5 -3.5", "route cdf 1 design 3 0 0 -13.5 -5 0 0 5 -3.5",
  "route cdf 2 design 3 0 0 NaN -5 0 0 5 -3.5", "route cdf 0 design 9 0 0 NaN -5 0 0 5 -3.5",
  "route cdf 0 design 2 2 0 NaN -5 -5 0 5 -3.5",
  "route icdf 0 case 1 0 0 NaN 0 0 0 2 0.050000000000000003", "route icdf 2 case 1 0 0 NaN 0 0 0 2 0.050000000000000003",
  "route mom 0 case 2 0 0 NaN -5 0 0 1 4", "route mom 1 case 2 0 0 NaN -5 0 0 1 4", "route mom 2 case 2 0 0 NaN -5 0 0 1 4",
  "route coef 0 case 2 0 0 NaN 0 0 0 3 -20", "route coef 1 case 2 0 0 NaN 0 0 0 3 -20",
  "route coef 2 case 2 0 0 NaN 0 0 0 3 -20", "route coefq 2 case 2 0 0 NaN 0 0 0 1 0.050000000000000003",
  "route lin 0 case 1 0 0 NaN -5 0 NaN 3 0.02", "route lin 0 case 1 0 0 NaN 0 0 -3.5 3 -0.59999999999999998",
  "route lin 1 case 1 0 0 NaN 0 0 -3.5 3 -0.59999999999999998", "route lin 1 case 1 0 0 NaN -5 0 NaN 3 0.02",
  "route lin 2 case 1 0 0 NaN -5 0 NaN 3 0.02", "route lin 2 case 1 0 0 NaN 0 0 -3.5 3 -0.59999999999999998",
  "route lin 0 case 4 0 0 -13.5 0 0 -7.2999999999999998 3 -0.90000000000000002",
  "route lin 1 case 4 0 0 -13.5 0 0 -7.2999999999999998 3 -0.90000000000000002",
  "route lin 2 case 4 0 0 -13.5 0 0 -7.2999999999999998 3 -0.90000000000000002",
  "route env 0 case 1 0 0 NaN 0 0 0.050000000000000003 2 -5", "route env 1 case 4 0 0 NaN 0 0 0.050000000000000003 1 -5",
  "route env 2 case 1 0 0 NaN 0 0 0.050000000000000003 2 -5",
  "route adf 0 case 2 0 0 NaN 0 0 0 3 -3", "route adf 1 case 2 0 0 NaN 0 0 0 3 -3", "route adf 2 case 2 0 0 NaN 0 0 0 3 -3",
  "route adfpdf 0 case 2 0 0 NaN 0 0 0 3 -3", "route adfpdf 1 case 2 0 0 NaN 0 0 0 3 -3",
  "route adfpdf 2 case 2 0 0 NaN 0 0 0 3 -3",
  "route region1 0 case 1 0 0 NaN -5 0 -3.5 1 -0.26390000000000002",
  "route region1 1 case 1 0 0 NaN -5 0 -3.5 1 -0.26390000000000002",
  "route region1 2 case 1 0 0 NaN -5 0 -3.5 1 -0.26390000000000002",
  "route region2 0 case 1 0 0 NaN -7 0 0 1 -1.9408700000000001",
  "route region2 1 case 1 0 0 NaN -7 0 0 1 -1.9408700000000001",
  "route region2 2 case 1 0 0 NaN -7 0 0 1 -1.9408700000000001")

det_key <- function(l) {
  w <- strsplit(l, " ", fixed = TRUE)[[1]]
  paste(w[seq_len(if (w[1] == "val") 12L else 13L)], collapse = " ")
}

# the lines (all, or those of det_cran) against their tolerances; the warnings of the engines are not checked here
# (the lines that the Julia engine flags as brittle have the tolerance 1e-6)
det_check_lines <- function(cran) {
  L <- det_ref_lines()
  D <- det_designs(L)
  L <- L[startsWith(L, "val ") | startsWith(L, "route ")]
  if (cran) {
    k <- vapply(L, det_key, "", USE.NAMES = FALSE)
    expect_true(all(det_cran %in% k), label = "every line of the CRAN subset is in ref/det_ref.txt")
    L <- L[k %in% det_cran]
  }
  for (l in L) {
    isval <- startsWith(l, "val ")
    v <- if (isval) det_parse_val(l) else det_parse_route(l)
    got <- suppressWarnings(if (isval) det_eval_val(v, D) else det_eval_route(v, D))
    dfl_empty_caches()
    d <- det_dev(got, v$want, v$fn)
    expect(length(got) == length(v$want) && all(d <= v$tol),
           sprintf("%s: deviation %.3g above the tolerance %.1g", substr(l, 1, 90), max(d), v$tol))
  }
  length(L)
}

test_that("the general engine and the router reproduce the Julia values (a subset)", {
  expect_equal(det_check_lines(TRUE), length(det_cran))
})

test_that("the general engine and the router reproduce the Julia values (every line)", {
  skip_on_cran()
  expect_gt(det_check_lines(FALSE), 450)
})

test_that("the constructors make the designs of the Julia engine", {
  D <- dfl_perron("C", c(0.7, 0.3))
  expect_s3_class(D, "dfl_design")
  expect_identical(D$br, c(0, 0.3, 0.7, 1))
  expect_identical(c(D$K, nrow(D$C), ncol(D$C)), c(6L, 2L, 18L))
  expect_output(print(D), "Perron's model C, lambda = 0.3, 0.7")
  expect_output(print(D), "6 regressors, 3 pieces")
  expect_identical(dfl_trend(-1)$K, 0L)
  expect_identical(dfl_trend(1)$C, diag(2))
  expect_equal(dfl_regmat(dfl_trend(1), 4), cbind(1, (1:4) / 4))
  X <- dfl_regmat(dfl_perron("A", 0.4), 50)
  expect_equal(X[, 3], as.numeric(1:50 > 20))
  # the break detected and the break given, and a dependent column dropped
  D1 <- dfl_design_from_matrix(X)
  D2 <- dfl_design_from_matrix(X, breaks = 20)
  expect_equal(D1$br, c(0, 0.4, 1))
  expect_equal(D1$C, D2$C, tolerance = 1e-12)
  expect_warning(D3 <- dfl_design_from_matrix(cbind(X, 2 * X[, 2])), "linearly dependent")
  expect_identical(D3$K, 3L)
  expect_true(dfl_detrend_first(dfl_trend(1))$first)
  expect_output(print(dfl_detrend_first(dfl_trend(1))), "detrends first")
  # the same span in the Legendre basis and from a function: the same statistic
  F1 <- dist_cdf(c(-3.5, -2), dfl_trend(2))
  expect_equal(dist_cdf(c(-3.5, -2), dfl_ltrend(2)), F1, tolerance = 1e-10)
  expect_equal(dist_cdf(c(-3.5, -2), dfl_fdesign(function(r) c(1, r, r^2))), F1, tolerance = 1e-10)
  expect_identical(dfl_fdesign(function(r) c(1, r), name = "line")$name, "line")
})

test_that("the constructors check their arguments; errors of the engine and of the function become R errors", {
  expect_error(dfl_trend(-2), "p >= -1")
  expect_error(dfl_trend(0.5), "single integer")
  expect_error(dfl_perron("D", 0.5), "model must be")
  expect_error(dfl_perron("A", 1.2), "0 < lambda_1 < ... < lambda_q < 1", fixed = TRUE)
  expect_error(dfl_ltrend(-1), "p >= 0")
  expect_error(dfl_fdesign(function(r) stop("no value here"), K = 2), "no value here")
  expect_error(dfl_fdesign(function(r) if (r > 0.7) c(NA, r) else c(1, r)), "not finite")
  expect_error(dfl_fdesign(function(r) if (r > 0.7) 1 else c(1, r)), "does not return 2 numbers")
  expect_error(dfl_fdesign(function(r) c(1, r), breaks = 1.5), "strictly inside")
  expect_error(dfl_fdesign(function(r) c(1, abs(r - 0.5))), "do not fall below")
  pulse <- cbind(1, as.numeric(1:100 == 50))
  expect_error(dfl_design_from_matrix(pulse))
  expect_error(dfl_design_from_matrix(matrix(1:4, 2), breaks = 3), "breaks must be")
  expect_error(dfl_detrend_first(1), "must be a design")
  expect_error(dist_cdf(0, "DU"), "must be a design")
  bad <- dfl_trend(1)
  bad$C <- bad$C[, 1, drop = FALSE]
  expect_error(dist_cdf(0, bad), "coefficient matrix")
})

test_that("the starts, GLS detrending and the options are checked", {
  D <- dfl_trend(0)
  expect_error(dist_cdf(0, D, 1, start = "stationary"), "stationary start needs c < 0")
  expect_error(dist_cdf(0, D, start = "random"), "start must be")
  expect_error(dist_cdf(0, D, start = "zero", gamma0 = 1), "exclude each other")
  expect_error(dist_cdf(0, D, -5, start = "stationary", gamma0 = 1), "exclude each other")
  expect_error(dist_cdf(0, D, gls = 1), "gls must be negative")
  expect_error(dist_quantile(1, D), "(0, 1)", fixed = TRUE)
  expect_error(dist_moments(2, dfl_detrend_first(D)), "not implemented")
  expect_error(lm_cdf(-2, D), "no trend function that vanishes at the origin")
  expect_error(pt_cdf(2, D, cbar = 1), "cbar must be negative")
  expect_error(region_prob(list(type = "wap", k = 1, omega = 1), D), "line\" or \"parabola")
  expect_error(ltratio_cdf(0, "DU", start = "stationary"), "apply to a design")
  expect_error(ltratio_cdf(0, "DU", engine = "fast"), "engine must be")
  expect_error(ltratio_cdf(0, dfl_perron("A", 0.5), engine = "classical"), "engine = \"classical\"")
  expect_error(ltratio_cdf(0, "DU", method = "single", engine = "general"), "representation of the classical engine")
  expect_error(adfmax_cdf(0, dfl_perron("A", 0.5), engine = "classical"), "ADFmax of the classical engine")
  expect_error(ltratio_lincdf(0, 1, dfl_perron("A", 0.5), engine = "classical"), "cases DK and DG only")
  expect_error(ltratio_region(list(type = "line", a = 1, b = 0), dfl_trend(0), engine = "classical"), "case DK only")
  expect_error(ltratio_region(list(type = "line", a = 1, b = 0, v0 = 0.1), "DK", engine = "general"), "v0")
})

test_that("with a case code and no engine the classical engine computes as before", {
  # (values computed by separate calls are compared to rounding: with threads the order of summation may differ by one
  # ulp)
  expect_identical <- function(a, b) expect_equal(a, b, tolerance = 1e-14)
  x <- c(-3, -2, -1)
  F <- ltratio_cdf(x, "DU", -5)
  expect_identical(ltratio_cdf(x, "DU", -5, engine = "classical"), F)
  expect_identical(ltratio_cdf(x, dfl_trend(0), -5, engine = "classical"), F)
  expect_identical(ltratio_pdf(x, "ST", -5, engine = "classical"), ltratio_pdf(x, "ST", -5))
  expect_identical(ltratio_pdf(x, dfl_trend(1), -5, start = "stationary", engine = "classical"), ltratio_pdf(x, "ST", -5))
  expect_identical(ltratio_cdf(x, dfl_trend(1), -5, gls = -13.5, engine = "classical"), ltratio_cdf(x, "DG", -5))
  expect_identical(ltratio_mom(2, dfl_trend(-1), -5, g0 = 1, engine = "classical"), ltratio_mom(2, "DK", -5, g0 = 1))
  expect_identical(coef_cdf(-10, dfl_trend(0), engine = "classical"), coef_cdf(-10, "DU"))
  expect_identical(ltratio_lincdf(-0.2639, -3.5, dfl_trend(-1), -7, engine = "classical"),
                   ltratio_lincdf(-0.2639, -3.5, "DK", -7))
  expect_identical(power_envelope(-7, dfl_trend(0), engine = "classical"), power_envelope(-7, "DK"))
  # a design that is a classical case takes the router's rule, as the case code with engine = "auto"
  expect_identical(ltratio_cdf(x, dfl_trend(0), -5), ltratio_cdf(x, "DU", -5, engine = "auto"))
  expect_identical(ltratio_cdf(x, dfl_trend(0), -5, c0 = 2), ltratio_cdf(x, "DU", -5, c0 = 2, engine = "auto"))
  # the general engine, through the router and on its own
  expect_identical(ltratio_cdf(x, "DT", -5, engine = "general"), dist_cdf(x, dfl_trend(1), -5))
  expect_equal(ltratio_cdf(x, "DT", -5, engine = "general"), ltratio_cdf(x, "DT", -5), tolerance = 1e-9)
  # the moments of a design of case DG: the general engine (the classical engine has none)
  expect_error(ltratio_mom(2, "DG"), "not available")
  expect_equal(ltratio_mom(2, dfl_trend(1), gls = -13.5), dist_moments(2, dfl_trend(1), gls = -13.5))
})

test_that("the companions agree with each other", {
  D <- dfl_trend(0)
  q <- pt_quantile(0.05, D, -7)
  p <- pt_power(-7, D, -7)
  expect_equal(attr(p, "critical"), q, tolerance = 1e-8)
  expect_equal(as.numeric(p), pt_cdf(q, D, -7, -7), tolerance = 1e-8)
  expect_equal(pt_envelope(c(-7, 0, NaN), D), c(as.numeric(p), 0.05, NaN), tolerance = 1e-8)
  # (power_envelope with a design other than a trend of order -1, 0, 1: the envelope of its point-optimal statistic)
  expect_identical(power_envelope(c(0, NaN), dfl_perron("A", 0.5)), c(0.05, NaN))
  expect_equal(pt_cdf(c(-Inf, Inf, NaN), D, -7), c(0, 1, NaN))
  expect_equal(msb_cdf(c(-Inf, 0, Inf, NaN), dfl_trend(1)), c(0, 0, 1, NaN))
  expect_equal(lm_cdf(c(-Inf, 1, NaN), dfl_trend(1)), c(0, 1, NaN))
  # the law of V below its support, by the general engine and through the router
  expect_identical(v_cdf(c(-1, -Inf), D), c(0, 0))
  expect_identical(v_pdf(-1, D), 0)
  expect_identical(ltratio_lincdf(-1, NaN, "DK", engine = "general"), 0)
  P <- lm_power(-5, dfl_trend(1))
  expect_equal(attr(P, "critical"), lm_quantile(0.05, dfl_trend(1)), tolerance = 1e-12)
  expect_equal(as.numeric(P), lm_cdf(attr(P, "critical"), dfl_trend(1), -5), tolerance = 1e-12)
})

test_that("the shapes of the inputs are kept and non-finite points have their limits", {
  x <- matrix(c(-Inf, -2, NaN, Inf), 2, dimnames = list(c("a", "b"), NULL))
  D <- dfl_perron("A", 0.5)
  F <- dist_cdf(x, D)
  expect_identical(dim(F), dim(x))
  expect_identical(dimnames(F), dimnames(x))
  expect_identical(F[c(1, 3, 4)], c(0, NaN, 1))
  expect_identical(dist_pdf(c(-Inf, Inf), D), c(0, 0))
  expect_named(ltratio_cdf(c(a = -3, b = -2), D), c("a", "b"))
  expect_identical(dim(ltratio_icdf(matrix(c(0.05, 0.1), 1), D)), c(1L, 2L))
  expect_identical(coef_pdf(numeric(0), D), numeric(0))
  r <- ltratio_icdf(0.05, D, info = TRUE)
  expect_named(r, c("q", "F", "iter", "conv"))
  expect_equal(r$F, 0.05, tolerance = 1e-9)
})

test_that("warnings of the engines become R warnings, once, and are not given again by a later call", {
  # (the general engine on the vertical line with GLS detrending at small |cbar|: a warning of the engine)
  expect_warning(dist_cdf(-2.5, dfl_trend(1), -5, gls = -2), "vertical line")
  expect_silent(ltratio_cdf(-2.86157, "DU", -5, engine = "general"))
  # (a brittle regime of the classical engine's one evaluation, through the router)
  expect_warning(ltratio_cdf(-2, "DU", -20, c0 = 1, g0 = 1, engine = "auto"), "c0 > 0 far from the unit root")
})

test_that("a classical warning is not given again by a call of the general engine (not on CRAN)", {
  skip_on_cran()
  expect_warning(ltratio_cdf(0.05, "SK", -1e-3, c0 = -1e-3, engine = "auto"), "flags some of them", fixed = TRUE)
  expect_silent(ltratio_cdf(-2, "SK", -5, engine = "general"))
})
