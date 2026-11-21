test_that("published critical values and moments are reproduced", {
  expect_equal(ltratio_icdf(0.05, "DU"), -2.86157, tolerance = 1e-5)       # 5% critical value of tau_mu
  expect_equal(ltratio_icdf(0.05, "DK"), -1.94087, tolerance = 1e-5)       # 5% critical value of tau
  expect_equal(ltratio_mom(2, "DU")[2], 3.0559953450, tolerance = 1e-10)   # Nabeya (1999)
  expect_equal(adfmax_cdf(-2.43542, "DU"), 0.05, tolerance = 1e-5)         # the quantiles of the ADFmax statistic in the paper
  expect_equal(ltratio_region(list(type = "wap", k = 3.6658, omega = 10)), 0.05, tolerance = 1e-4)   # the weighted-average-power test of the paper
  expect_equal(ltratio_lincdf(-0.2639, -3.5, "DK"), 0.05, tolerance = 1e-4)   # the point-optimal test at cbar = -7
  expect_equal(power_envelope(-7), 0.5, tolerance = 1e-2)                     # power one half at c = -7
})

test_that("the power at c = -5 is the distribution function at the critical value", {
  q <- ltratio_icdf(0.05, "DU")
  expect_equal(ltratio_cdf(q, "DU"), 0.05, tolerance = 1e-10)
  expect_equal(ltratio_cdf(-2.86157, "DU", -5), 0.116513562093, tolerance = 1e-9)
})

test_that("the quantiles invert the distribution function", {
  r <- ltratio_icdf(c(0.1, 0.9), "ST", -10, info = TRUE)
  expect_true(all(r$conv))
  expect_equal(r$F, c(0.1, 0.9), tolerance = 1e-9)
  expect_equal(ltratio_cdf(r$q, "ST", -10), c(0.1, 0.9), tolerance = 1e-9)
  expect_equal(adfmax_cdf(adfmax_icdf(0.05, "DU"), "DU"), 0.05, tolerance = 1e-9)
})

test_that("the two representations agree", {
  x <- c(-3, -1.5, 0.5)
  expect_equal(ltratio_cdf(x, "DU", -2, method = "line"), ltratio_cdf(x, "DU", -2, method = "conditional"),
               tolerance = 1e-10)
  expect_equal(ltratio_cdf(-1, "DK", 2, c0 = 2), ltratio_region(list(type = "parabola", x = -1, c0 = 2), "DK", 2),
               tolerance = 1e-10)
  # the line of the point-optimal test by the conditional integral and by the inversion of U - aV
  expect_equal(ltratio_region(list(type = "line", a = -3.5, b = -0.2639), "DK", -5),
               ltratio_lincdf(-0.2639, -3.5, "DK", -5), tolerance = 1e-10)
})

test_that("the linear form near a = 0 and outside its support", {
  # a = 0: P(U <= b) = P(chi^2_1 <= 2 kappa (b + 1/2)), kappa = 2c/(e^{2c} - 1) (1 at c = 0), zero below the support
  for (c in c(-50, -5, 0, 2)) {
    kap <- if (c == 0) 1 else 2 * c / expm1(2 * c)
    b <- -0.5 + c(-0.1, 0, 1e-6, 1e-3, 0.05, 0.5, 2) / kap
    expect_equal(ltratio_lincdf(b, 0, "DK", c), pchisq(pmax(2 * kap * (b + 0.5), 0), 1), tolerance = 1e-13)
  }
  expect_identical(ltratio_lincdf(-0.6, 0, "DK", -5), 0)
  # small |a| and c > 0, against the Julia values (the conditional integral given U in case DG)
  expect_equal(ltratio_lincdf(-0.51, 0.3, "DG", -5), 0.291315235195, tolerance = 1e-11)
  expect_equal(ltratio_lincdf(-0.51, 0.3, "DG", 0), 0.287795093100, tolerance = 1e-11)
  P <- ltratio_lincdf(seq(-0.7, 1, by = 0.1), 1e-3, "DG", 2)
  expect_true(all(P >= 0 & P <= 1) && all(diff(P) >= 0))
  # the power envelope near c = 0: alpha + O(c^2)
  expect_equal(power_envelope(-0.01), 0.050249919491, tolerance = 1e-10)
  e <- power_envelope(-1e-3, "DG")
  expect_true(e > 0.05 && e < 0.05 + 1e-8)
  # the ends of the range
  expect_error(ltratio_lincdf(0, 0, "DK", -600), "at least -550")
  expect_error(ltratio_lincdf(0, 0, "DG", 6), "case DG is available for c up to about 3.5")
})
