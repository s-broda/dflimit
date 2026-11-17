# independent references: the Julia values stored in public/ltratio_test.m of DF-Limit (from julia/DFExact.jl and
# julia/application) and the 5% critical values printed in the paper (public/ltratio_paper.m, public/Contents.m)
test_that("coef_cdf and coef_icdf reproduce the stored Julia values", {
  expect_lt(abs(coef_cdf(-14, "DU", -5) - 0.1938495145682877), 1e-10)
  expect_lt(abs(coef_cdf(-14, "SU", -10) - 0.45401150754561903), 1e-10)
  expect_lt(abs(coef_cdf(-16.59, "DG", -13.5) - 0.51132404681639121), 1e-10)
  expect_lt(abs(coef_icdf(0.05, "DU") + 14.093590831267115), 1e-7)
})

test_that("coef_icdf reproduces the published 5% critical values", {
  expect_equal(round(coef_icdf(0.05, "DU"), 4), -14.0936)                       # public/Contents.m
  expect_equal(round(coef_icdf(0.05, "DK"), 3), -8.039)                      # the paper: -8.039 -14.094 -21.711
  expect_equal(round(coef_icdf(0.05, "DT"), 3), -21.711)
  expect_equal(round(coef_icdf(0.05, "DG"), 2), -16.59)                         # the paper: -16.59
})

# self-consistency only (no stored value): the size at the critical value, round trips, and the definition
test_that("coef_icdf inverts coef_cdf", {
  expect_equal(coef_cdf(-14.0936, "DU"), 0.05, tolerance = 1e-5)
  expect_equal(coef_cdf(coef_icdf(0.1, "DK"), "DK"), 0.1, tolerance = 1e-8)
  expect_equal(coef_cdf(coef_icdf(0.9, "DU"), "DU"), 0.9, tolerance = 1e-8)
  expect_equal(coef_cdf(coef_icdf(0.05, "DT", -5), "DT", -5), 0.05, tolerance = 1e-8)
  expect_equal(coef_cdf(coef_icdf(0.1, "SU", -5), "SU", -5), 0.1, tolerance = 1e-8)
})

test_that("coef_cdf is the distribution function of the t-ratio at 0 with the shifted null value", {
  # equal to rounding, not identical: with threads the order of summation in the engine can differ by one ulp
  expect_equal(coef_cdf(-5, "DT", -2, c0 = 1), ltratio_cdf(0, "DT", -2, c0 = -4), tolerance = 1e-13)
  expect_equal(coef_cdf(-3, "DK", 0, g0 = 1), ltratio_cdf(0, "DK", 0, c0 = -3, g0 = 1), tolerance = 1e-13)
})

test_that("coef_cdf and coef_icdf keep the shape and check their arguments", {
  y <- matrix(c(-Inf, -5, NaN, Inf), 2, dimnames = list(c("a", "b"), NULL))
  F <- coef_cdf(y, "DK")
  expect_identical(dim(F), dim(y))
  expect_identical(dimnames(F), dimnames(y))
  expect_identical(F[c(1, 3, 4)], c(0, NaN, 1))
  expect_named(coef_cdf(c(a = -1), "DU"), "a")
  expect_identical(coef_cdf(numeric(0), "DU"), numeric(0))
  expect_error(coef_cdf(0, "XX"), "the case must be")
  expect_error(coef_cdf(0, "DU", method = "fast"), "method must be")
  expect_error(coef_cdf(-1, "SU", 1), "requires c <= -1e-4", fixed = TRUE)
  expect_error(coef_icdf(0, "DU"), "(0, 1)", fixed = TRUE)
  expect_error(coef_icdf(0.5, "DU", c = NA), "single number")
})
