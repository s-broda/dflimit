# The repairs of 0.5.0: the points where the line integral exceeds the bounds of the engine and the fixed start whose path
# climbed.  The values: the Julia engine DFLimit.jl 0.5.0, which has the same repairs, and the contour integral of a newer
# engine (the values of the Julia tests).

test_that("large x above the unit root goes to the conditional integral (the engine crashed here)", {
  expect_equal(ltratio_cdf(751, "DT", 7), 0.990823749760809, tolerance = 1e-12)
  expect_equal(ltratio_pdf(751, "DT", 7), 9.2924954e-5, tolerance = 1e-6)
  expect_error(ltratio_cdf(751, "DT", 7, method = "single"), "at large x above the unit root", fixed = TRUE)
})

test_that("c0 < 0 near zero with x > 0 goes to the conditional integral (more than 2 GB before)", {
  expect_warning(F <- ltratio_cdf(0.05, "SK", -1e-3, c0 = -1e-3), "flags some of them", fixed = TRUE)
  expect_equal(F, 0.54749158678338294, tolerance = 1e-10)
  expect_warning(F <- ltratio_cdf(1, "SK", -1e-2, c0 = -1e-2), "flags some of them", fixed = TRUE)
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
