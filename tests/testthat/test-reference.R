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
    expect_no_warning(Fc <- ltratio_cdf(r$x, cas, r$c, r$c0, r$g0, method = "double"))
    expect_no_warning(fc <- ltratio_pdf(r$x, cas, r$c, r$c0, r$g0, method = "double"))
    expect_close_Ff(Fc, fc, r$F, r$f, 1e-10, 1e-8, 1e-10,
                    sprintf("%s c = %g c0 = %g g0 = %g x = %g", cas, r$c, r$c0, r$g0, r$x))
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
