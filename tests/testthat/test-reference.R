test_that("the line integral at c0 = 0 reproduces the Julia values", {
  d <- read_ref("single_ref.txt", c("m", "c", "st", "x", "F", "f"))
  for (i in seq_len(nrow(d))) {
    r <- d[i, ]
    cas <- case_code(r$m, r$st)
    expect_close_Ff(ltratio_cdf(r$x, cas, r$c, method = "single"), ltratio_pdf(r$x, cas, r$c, method = "single"),
                    r$F, r$f, 1e-11, 1e-9, 1e-11, sprintf("%s c = %g x = %g", cas, r$c, r$x))
  }
})

test_that("the line integral with c0 != 0 reproduces the Julia values", {
  d <- read_ref("gen_ref.txt", c("m", "c", "c0", "st", "x", "F", "f"))
  for (i in seq_len(nrow(d))) {
    r <- d[i, ]
    cas <- case_code(r$m, r$st)
    expect_close_Ff(ltratio_cdf(r$x, cas, r$c, r$c0, method = "single"),
                    ltratio_pdf(r$x, cas, r$c, r$c0, method = "single"),
                    r$F, r$f, 1e-11, 1e-9, 1e-11, sprintf("%s c = %g c0 = %g x = %g", cas, r$c, r$c0, r$x))
  }
})

test_that("the line integral with a fixed start reproduces the Julia values", {
  d <- read_ref("fixed_ref.txt", c("m", "c", "c0", "g0", "st", "x", "F", "f"))
  for (i in seq_len(nrow(d))) {
    r <- d[i, ]
    cas <- case_code(r$m, r$st)
    expect_close_Ff(ltratio_cdf(r$x, cas, r$c, r$c0, r$g0, method = "single"),
                    ltratio_pdf(r$x, cas, r$c, r$c0, r$g0, method = "single"),
                    r$F, r$f, 1e-11, 1e-9, 1e-11, sprintf("%s c = %g c0 = %g g0 = %g x = %g", cas, r$c, r$c0, r$g0, r$x))
  }
})

test_that("the conditional integral reproduces the Julia values", {
  d <- read_ref("dbl_ref.txt", c("m", "c", "c0", "g0", "st", "x", "F", "f"))
  for (i in seq_len(nrow(d))) {
    r <- d[i, ]
    cas <- case_code(r$m, r$st)
    lab <- sprintf("%s c = %g c0 = %g g0 = %g x = %g", cas, r$c, r$c0, r$g0, r$x)
    # verified: the checked computation at the tolerances of the C engine's tests, without a warning
    expect_no_warning(Fc <- ltratio_cdf(r$x, cas, r$c, r$c0, r$g0, method = "double", verified = TRUE))
    expect_no_warning(fc <- ltratio_pdf(r$x, cas, r$c, r$c0, r$g0, method = "double", verified = TRUE))
    expect_close_Ff(Fc, fc, r$F, r$f, 1e-10, 1e-8, 1e-10, lab)
    # the default, one evaluation: at the tolerance it meets (5e-7 in F, 1e-6 relative or 1e-5 absolute in f), with
    # the warning of a brittle regime where it lies in one
    w <- NULL
    F1 <- withCallingHandlers(ltratio_cdf(r$x, cas, r$c, r$c0, r$g0, method = "double"),
                              warning = function(e) { w <<- conditionMessage(e); invokeRestart("muffleWarning") })
    f1 <- suppressWarnings(ltratio_pdf(r$x, cas, r$c, r$c0, r$g0, method = "double"))
    if (is.null(w)) expect_close_Ff(F1, f1, r$F, r$f, 5e-7, 1e-6, 1e-5, paste(lab, "(one evaluation)"))
  }
})

test_that("the moments reproduce the Julia values", {
  d <- read_ref("mom_ref.txt", c("m", "c", "c0", "g0", "st", "m1", "m2", "m3", "m4"))
  for (i in seq_len(nrow(d))) {
    r <- d[i, ]
    mr <- c(r$m1, r$m2, r$m3, r$m4)
    expect_lt(max(abs(ltratio_mom(4, case_code(r$m, r$st), r$c, r$c0, r$g0) / mr - 1)), 1e-10)
  }
})

test_that("ADFmax, the regions, the linear form and the power envelope reproduce the Julia values", {
  rows <- strsplit(grep("^#", readLines(test_path("ref", "app_ref.txt")), invert = TRUE, value = TRUE), " ")
  # the pairs (point, value) of a line from position 'from' on
  pairs <- function(v, from) matrix(v[from:length(v)], 2)
  for (v in rows) {
    tag <- v[1]
    v <- as.numeric(v[-1])
    if (tag == "lincdf") {
      p <- pairs(v, 6)
      expect_lt(max(abs(ltratio_lincdf(p[1, ], v[2], c("DK", "", "", "DG")[v[1]], v[3], v[4]) - p[2, ])), 1e-10)
    } else if (tag == "vcdf") {
      p <- pairs(v, 3)
      expect_lt(max(abs(ltratio_lincdf(p[1, ], NaN, "DK", v[1]) - p[2, ])), 1e-10)
    } else if (tag == "envelope") {
      p <- pairs(v, 4)
      expect_lt(max(abs(power_envelope(p[1, ], c("DK", "", "", "DG")[v[1]], v[2]) - p[2, ])), 1e-10)
    } else if (tag == "adfmax") {
      cas <- paste0(if (v[2] == 1) "S" else "D", c("K", "U", "T")[v[1]])
      expect_close_Ff(adfmax_cdf(v[4], cas, v[3]), adfmax_pdf(v[4], cas, v[3]), v[5], v[6], 1e-10, 0, 1e-10,
                      sprintf("ADFmax %s c = %g x = %g", cas, v[3], v[4]))
    } else {
      bnd <- switch(v[1], list(type = "line", a = v[2], b = v[3]), list(type = "parabola", x = v[2], c0 = v[3]),
                    list(type = "wap", k = v[2], omega = v[3]))
      bnd$v0 <- v[4]
      expect_lt(abs(ltratio_region(bnd, "DK", v[5]) - v[6]), 1e-10)
    }
  }
})

test_that("the contour integral with a fixed start in cases DU, DT reproduces the Julia values", {
  d <- read_ref("contour_ref.txt", c("m", "c", "c0", "g0", "x", "F", "f", "eF", "ef"))
  for (i in seq_len(nrow(d))) {
    r <- d[i, ]
    cas <- case_code(r$m, FALSE)
    # verified: the stored values, without a warning; method "auto" takes the contour integral there
    expect_no_warning(Fc <- ltratio_cdf(r$x, cas, r$c, r$c0, r$g0, method = "contour", verified = TRUE))
    expect_no_warning(fc <- ltratio_pdf(r$x, cas, r$c, r$c0, r$g0, method = "contour", verified = TRUE))
    lab <- sprintf("%s c = %g c0 = %g g0 = %g x = %g", cas, r$c, r$c0, r$g0, r$x)
    expect_close_Ff(Fc, fc, r$F, r$f, max(1e-12, r$eF), 1e-9, max(1e-12, r$ef), lab)
    # (equal to rounding, not identical: with threads the order of summation in the engine can differ by one ulp)
    if (Fc >= 1e-3 && fc >= 1e-3) expect_equal(ltratio_cdf(r$x, cas, r$c, r$c0, r$g0, verified = TRUE), Fc, tolerance = 1e-13, label = lab)
    # the default, one evaluation (with the warning of a brittle regime where it lies in one, c0 > 0 at c < -15)
    F1 <- suppressWarnings(ltratio_cdf(r$x, cas, r$c, r$c0, r$g0, method = "contour"))
    f1 <- suppressWarnings(ltratio_pdf(r$x, cas, r$c, r$c0, r$g0, method = "contour"))
    expect_equal(suppressWarnings(ltratio_cdf(r$x, cas, r$c, r$c0, r$g0)), F1, tolerance = 1e-13, label = lab)
    expect_close_Ff(F1, f1, r$F, r$f, 5e-9, 1e-8, 5e-9, paste(lab, "(one evaluation)"))
  }
  expect_error(ltratio_cdf(-2, "DK", -5, g0 = 1, method = "contour"), "applies to a fixed start")
  expect_error(ltratio_cdf(-2, "DU", -5, method = "contour"), "applies to a fixed start")
})

test_that("the warnings of the brittle regimes fire where they should and not elsewhere", {
  expect_warning(ltratio_cdf(-3, "DU", -5, c0 = 2, method = "double"), "one evaluation")
  expect_warning(ltratio_cdf(c(-3, 1), "DK", 2, c0 = -5, method = "double"), "one evaluation")
  expect_warning(ltratio_cdf(-5, "DU", -20, c0 = 2, g0 = 0.5), "one evaluation")
  expect_warning(ltratio_cdf(-2, "SK", -0.1, c0 = 1, method = "double"), "one evaluation")
  expect_no_warning(ltratio_cdf(-3, "DU", -5, c0 = 2))
  expect_no_warning(ltratio_cdf(-2, "DU", -5))
  expect_no_warning(ltratio_cdf(-2, "DU", -5, g0 = 1.5))
  expect_no_warning(ltratio_cdf(-3, "DU", -5, c0 = 2, method = "double", verified = TRUE))
  expect_warning(ltratio_cdf(-8.4, "DU", 2, c0 = 5, g0 = 1.5, method = "contour"), "contour integral by one evaluation fails")
  # case DU with a fixed start, c0 > 3 at c > 1.5: the contour integral fails there, "auto" takes the conditional integral
  expect_no_warning(Fa <- ltratio_cdf(-8.4075, "DU", 2, c0 = 5, g0 = 1.5))
  expect_lt(abs(Fa - ltratio_cdf(-8.4075, "DU", 2, c0 = 5, g0 = 1.5, method = "double", verified = TRUE)), 1e-9)
  # once per call of the functions that evaluate the distribution many times
  w <- character(0)
  withCallingHandlers(coef_cdf(c(5, 10), "DU", -10, method = "double"),
                      warning = function(e) { w <<- c(w, conditionMessage(e)); invokeRestart("muffleWarning") })
  expect_equal(length(w), 1L)
  # (the quantiles below take about 3 s on two threads, a fifth of the suite: not on CRAN)
  skip_on_cran()
  w <- character(0)
  withCallingHandlers(ltratio_icdf(c(0.05, 0.5), "DU", -20, c0 = 2, g0 = 0.5),
                      warning = function(e) { w <<- c(w, conditionMessage(e)); invokeRestart("muffleWarning") })
  expect_equal(length(w), 1L)
})

test_that("verified = TRUE reproduces the checked computation", {
  # the values of the engine before the default became one evaluation (November 2026; quickcheck.jl, section 21)
  vr <- list(list("DU", -5, 2, 0, -3, 0.1976948237625359, 0.4648778393052269),
             list("DK", -10, -5, 0, -2, 0.16726497409337868, 0.3125724838563969),
             list("SU", -5, -5, 0, -1, 0.44272919976045844, 0.4264723317698583),
             list("DU", -5, 0, 0, -4.5, 0.0005234905959262681, 0.0022074999304135296))
  for (v in vr) {
    expect_lt(abs(ltratio_cdf(v[[5]], v[[1]], v[[2]], v[[3]], v[[4]], verified = TRUE) - v[[6]]), 1e-12)
    expect_lt(abs(ltratio_pdf(v[[5]], v[[1]], v[[2]], v[[3]], v[[4]], verified = TRUE) - v[[7]]), 1e-12)
    expect_lt(abs(ltratio_cdf(v[[5]], v[[1]], v[[2]], v[[3]], v[[4]]) - v[[6]]), 1e-9)
  }
})
