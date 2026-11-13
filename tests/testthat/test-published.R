test_that("published critical values and moments are reproduced", {
  expect_equal(ltratio_icdf(0.05, "DU"), -2.86157, tolerance = 1e-5)       # 5% critical value of tau_mu
  expect_equal(ltratio_icdf(0.05, "DK"), -1.94087, tolerance = 1e-5)       # 5% critical value of tau
  expect_equal(ltratio_mom(2, "DU")[2], 3.0559953450, tolerance = 1e-10)   # Nabeya (1999)
  expect_equal(adfmax_cdf(-2.43542, "DU"), 0.05, tolerance = 1e-5)         # Table I of the paper
  expect_equal(ltratio_region(list(type = "wap", k = 3.6658, omega = 10)), 0.05, tolerance = 1e-4)   # Table V
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
})
