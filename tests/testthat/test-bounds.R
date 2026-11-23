# The repairs of 0.5.0: the points where the line integral exceeds the bounds of the engine and the fixed start whose path
# climbed.  The values: the Julia engine DFLimit.jl 0.5.0, which has the same repairs, and the contour integral of a newer
# engine (the values of the Julia tests).

test_that("large x above the unit root goes to the conditional integral (the engine crashed here)", {
  expect_equal(ltratio_cdf(751, "DT", 7), 0.990823749760809, tolerance = 1e-12)
  expect_equal(ltratio_pdf(751, "DT", 7), 9.2924954e-5, tolerance = 1e-6)
  expect_error(ltratio_cdf(751, "DT", 7, method = "single"), "at large x above the unit root", fixed = TRUE)
})

test_that("c0 < 0 near zero with x > 0 goes to the conditional integral (more than 2 GB before)", {
  # (the classical engine: the default takes the general engine with a stationary start and c0 != 0 for c <= -3e-3)
  expect_warning(F <- ltratio_cdf(0.05, "SK", -1e-3, c0 = -1e-3, engine = "classical"), "flags some of them", fixed = TRUE)
  expect_equal(F, 0.54749158678338294, tolerance = 1e-10)
  expect_warning(F <- ltratio_cdf(1, "SK", -1e-2, c0 = -1e-2, engine = "classical"), "flags some of them", fixed = TRUE)
  expect_equal(F, 0.87075930276798486, tolerance = 1e-10)
  expect_warning(F <- ltratio_cdf(c(-1, 0.05, 0.5, 1, 1.5), "SK", -1e-4, c0 = -1e-4), "exceed the bounds of the engine",
                 fixed = TRUE)
  expect_equal(F, c(0.1666135587511895, 0.531233538152, 0.700742367308, 0.847290600876, 0.936188673511),
               tolerance = 1e-6)
})

test_that("the line integral requested at such a point is an error that names the arguments", {
  e <- tryCatch(ltratio_cdf(1, "SK", -1e-4, c0 = -1e-4, method = "single"), error = conditionMessage)
  expect_match(e, "the line integral needs more nodes than the engine allows", fixed = TRUE)
  expect_match(e, "c0 = -0.0001, x = 1", fixed = TRUE)
  expect_error(ltratio_cdf(1, "DG", -1, c0 = -1e-4), "case DG has the line integral representation only", fixed = TRUE)
})

test_that("a fixed start with c0 near zero at x > 0 (6e74 before)", {
  expect_equal(ltratio_cdf(1.5, "DK", 0, c0 = -0.01, g0 = 1), 0.960991862318, tolerance = 1e-10)
  expect_equal(ltratio_pdf(1.5, "DK", 0, c0 = -0.01, g0 = 1), 0.079242719297, tolerance = 1e-10)
  expect_equal(ltratio_cdf(1.5, "DK", 2, c0 = -0.01, g0 = 1), 0.086764388544, tolerance = 1e-10)
})

# The upper tail above the unit root (the C engine of dflimit-matlab, November 23, 2026): case DK with c0 != c beyond the
# line integral's bound, for c > 4 at x > 300 and with c0 > c at x < -20, the conditional integral with its checks (before,
# its one evaluation where the checked value was flagged: 0.8415, -1.3674, 1.0479 here); above c = 14 an error; the general
# engine's line cut to the rate with c0; the final check of a distribution function. The values: an independent
# computation (the law of V given J(1) by Girsanov's theorem, inverted and integrated over J(1), DFLimit.jl).
test_that("the upper tail above the unit root, case DK", {
  expect_equal(ltratio_cdf(1e4, "DK", 10, c0 = 2), 0.743623945123, tolerance = 1e-10)
  expect_equal(ltratio_cdf(3e4, "DK", 10, c0 = 2), 0.999338346257, tolerance = 1e-10)
  expect_equal(ltratio_cdf(3e3, "DK", 10, c0 = 9), 0.993550212928, tolerance = 1e-10)
  expect_equal(ltratio_cdf(-200, "DK", 10, c0 = 10.5), 0.716459029660, tolerance = 1e-10)
  expect_equal(ltratio_cdf(c(1e3, 3e3), "DK", 10, c0 = 9, engine = "general"), c(0.636120922790, 0.993550212928),
               tolerance = 1e-10)
  expect_error(ltratio_cdf(1e7, "DK", 16, c0 = 2), "is computed for c <= 14 only", fixed = TRUE)
})
