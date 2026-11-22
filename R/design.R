#' Designs of the general engine
#'
#' The general engine computes the limits of the Dickey-Fuller statistics for any piecewise smooth deterministic
#' design: the regression of \eqn{y_t} on \eqn{d(t/n)} and \eqn{y_{t-1}}, with \eqn{d(r)} a vector of \eqn{K}
#' regressors that are polynomials or smooth functions of \eqn{r = t/n} between known break dates. A design is made by
#' one of these constructors and passed in place of the case code to [ltratio_cdf()] and the other functions of the
#' classical engine (which then choose between the two engines, see the section Engines of [ltratio_cdf()]) or to the
#' functions of the general engine ([dist_cdf()] and the others).
#'
#' \describe{
#'   \item{`dfl_trend(p)`}{the polynomial trend \eqn{d(r) = (1, r, \dots, r^p)}{d(r) = (1, r, ..., r^p)}; `p = -1`
#'     is the regression without deterministic terms (case DK), `p = 0` the intercept (DU), `p = 1` the linear trend
#'     (DT), and any higher order is allowed;}
#'   \item{`dfl_perron(model, lambda)`}{Perron's models with breaks at the fractions `lambda` (one or several, inside
#'     (0, 1)): model `"A"` a break in the level, \eqn{(1, r, 1\{r > \lambda_i\})}{(1, r, 1{r > lambda_i})}, `"B"` a
#'     break in the slope, \eqn{(1, r, (r - \lambda_i)_+)}{(1, r, (r - lambda_i)+)}, `"C"` both at every break;}
#'   \item{`dfl_fdesign(f, breaks)`}{regressors given by a function `f(r)` that returns the `K` values at `r`, smooth
#'     between the `breaks`: each piece is interpolated by Chebyshev polynomials of the smallest degree
#'     \eqn{2^m - 1 \ge 15}{2^m - 1 >= 15} whose coefficients fall below `tol` times the largest, up to `maxdeg`;}
#'   \item{`dfl_ltrend(p)`}{the polynomial trend of order `p` in the orthonormal shifted Legendre basis, the same span
#'     as `dfl_trend(p)` and hence the same statistic, better conditioned for large `p`;}
#'   \item{`dfl_design_from_matrix(X, breaks)`}{the design of the limit from the `n` x `K` matrix `X` of the
#'     deterministic regressors of a regression (column `k` sampled at `t = 1, ..., n`, any scaling): the known
#'     breaks are the indices \eqn{T_B}{T_B} after which the regressors change (\eqn{\lambda = T_B/n}{lambda =
#'     T_B/n}), and without them the breaks are detected from the differences of the columns; each piece is then
#'     fitted by Chebyshev polynomials of increasing degree up to `maxdeg` until the largest residual is below `tol`.
#'     Columns linearly dependent on the others are dropped with a warning; columns that are not piecewise smooth at
#'     the sample's resolution (a piece shorter than `minlen` observations, as a pulse dummy makes, or seasonal
#'     dummies) are refused with an error: they do not affect the limit and should be left out;}
#'   \item{`dfl_detrend_first(D)`}{the statistic that detrends first: the data are detrended by least squares on the
#'     design `D`, and the t-ratio is that of the regression of the detrended data on their lag without deterministic
#'     terms.}
#' }
#'
#' `dfl_regmat(D, n)` gives the `n` x `K` matrix of the regressors \eqn{d(t/n)}, \eqn{t = 1, \dots, n}{t = 1..n}, of
#' a design (a regressor is left-continuous at a break), and `dfl_empty_caches()` frees the prepared designs and the
#' paths of the inversion that the general engine keeps between calls (up to 64 of each, reused when the same design
#' is evaluated again).
#'
#' A design is a list of class `"dfl_design"` with the break points `br` (from 0 to 1), the coefficient matrix `C`
#' (the matrices of the pieces side by side, one column per regressor and piece, one row per power of \eqn{r}, or of
#' \eqn{r} less the left end of the piece when `loc` is `TRUE`), the number of regressors `K`, the flags `loc` and
#' `first` (the statistic that detrends first) and a `name`.
#'
#' @param p The order of the polynomial trend (`-1` or more in `dfl_trend()`, `0` or more in `dfl_ltrend()`).
#' @param model Perron's model, `"A"`, `"B"` or `"C"`.
#' @param lambda The break fractions, one or more numbers in (0, 1).
#' @param f A function of a single number `r` in \[0, 1\] that returns the `K` regressors at `r`.
#' @param breaks `dfl_fdesign()`: the break fractions inside (0, 1) between which `f` is smooth (none by default).
#'   `dfl_design_from_matrix()`: the indices \eqn{T_B}{T_B} of the last observations before the breaks, or `NULL`
#'   (the default) to detect the breaks.
#' @param K The number of regressors (by default the length of `f(0.5)`).
#' @param tol `dfl_fdesign()`: the relative size below which the Chebyshev coefficients are dropped.
#'   `dfl_design_from_matrix()`: the largest residual of the fit of each column on each piece.
#' @param maxdeg The largest degree of the Chebyshev interpolants or fits.
#' @param name A name for the design, used by its print method.
#' @param X A numeric matrix of the deterministic regressors, one row per observation.
#' @param minlen The shortest piece, in observations, between two detected or given breaks.
#' @param D A design.
#' @param n The number of observations.
#' @param x A design (for `print()`).
#' @param ... Not used.
#' @return The constructors return a design (a list of class `"dfl_design"`); `dfl_regmat()` a numeric matrix;
#'   `dfl_empty_caches()` `NULL`, invisibly.
#' @seealso [ltratio_cdf()] (the router over the two engines), [dist_cdf()] (the functions of the general engine).
#' @examples
#' D <- dfl_perron("C", c(0.3, 0.7))          # breaks in level and slope at 0.3 and 0.7
#' D
#' dfl_trend(2)                               # a quadratic trend
#' dfl_fdesign(function(r) c(1, sin(2 * pi * r), cos(2 * pi * r)))
#' X <- dfl_regmat(dfl_perron("A", 0.4), 200) # the regressors of a sample of 200
#' dfl_design_from_matrix(X)                  # the break is found at t = 80
#' @name dfl_design
NULL

# a design from the list of the C gateway (br, C, K, loc, first, warning)
dfl_new_design <- function(r, name) {
  if (nzchar(r$warning)) warning(r$warning, call. = FALSE)
  r$warning <- NULL
  r$name <- name
  structure(r, class = "dfl_design")
}

dfl_int <- function(v, name) {
  if (!is.numeric(v) || length(v) != 1L || is.na(v) || v != round(v) || abs(v) > .Machine$integer.max) {
    stop(name, " must be a single integer", call. = FALSE)
  }
  as.integer(v)
}

#' @rdname dfl_design
#' @export
dfl_trend <- function(p) {
  p <- dfl_int(p, "p")
  nm <- switch(as.character(p), "-1" = "no deterministic terms (case DK)", "0" = "intercept (case DU)",
               "1" = "linear trend (case DT)", "2" = "quadratic trend", "3" = "cubic trend",
               sprintf("polynomial trend of order %d", p))
  dfl_new_design(.Call(C_dfl_r_design, "trend", p, NULL, NULL, NULL, NULL), nm)
}

#' @rdname dfl_design
#' @export
dfl_perron <- function(model, lambda) {
  if (!is.character(model) || length(model) != 1L || !toupper(model) %in% c("A", "B", "C")) {
    stop("model must be \"A\", \"B\" or \"C\"", call. = FALSE)
  }
  lambda <- dfl_real(lambda, "lambda")
  if (length(lambda) < 1L) stop("lambda must have at least one break fraction", call. = FALSE)
  model <- toupper(model)
  r <- .Call(C_dfl_r_design, "perron", model, lambda, NULL, NULL, NULL)
  dfl_new_design(r, sprintf("Perron's model %s, lambda = %s", model, paste(format(sort(lambda)), collapse = ", ")))
}

#' @rdname dfl_design
#' @export
dfl_fdesign <- function(f, breaks = numeric(0), K = NULL, tol = 1e-14, maxdeg = 96, name = "function design") {
  if (!is.function(f)) stop("f must be a function", call. = FALSE)
  breaks <- dfl_real(breaks, "breaks")
  if (is.null(K)) K <- length(f(0.5))
  K <- dfl_int(K, "K")
  if (K < 1L) stop("K must be at least 1", call. = FALSE)
  msg <- new.env()
  g <- function(r) {
    tryCatch({
      v <- f(r)
      if (!is.numeric(v) || length(v) != K) stop(sprintf("f(%g) does not return %d numbers", r, K), call. = FALSE)
      v <- as.double(v)
      if (!all(is.finite(v))) stop(sprintf("f(%g) is not finite", r), call. = FALSE)
      v
    }, error = function(e) {
      if (is.null(msg$text)) msg$text <- conditionMessage(e)
      NULL
    })
  }
  opt <- c(dfl_num(tol, "tol"), dfl_int(maxdeg, "maxdeg"), 0)
  r <- .Call(C_dfl_r_design, "fdesign", g, environment(), K, breaks, opt)
  if (is.null(r)) {
    stop("dfl_fdesign: ", if (is.null(msg$text)) "the evaluation of f failed" else msg$text, call. = FALSE)
  }
  dfl_new_design(r, dfl_name(name))
}

#' @rdname dfl_design
#' @export
dfl_ltrend <- function(p) {
  p <- dfl_int(p, "p")
  dfl_new_design(.Call(C_dfl_r_design, "ltrend", p, NULL, NULL, NULL, NULL),
                 sprintf("polynomial trend of order %d (Legendre basis)", p))
}

#' @rdname dfl_design
#' @export
dfl_design_from_matrix <- function(X, breaks = NULL, tol = 1e-8, maxdeg = 40, minlen = 10,
                                   name = "design from a matrix") {
  if (is.data.frame(X)) X <- as.matrix(X)
  if (!is.numeric(X) || !is.matrix(X)) stop("X must be a numeric matrix", call. = FALSE)
  storage.mode(X) <- "double"
  if (!all(is.finite(X))) stop("X must be finite", call. = FALSE)
  if (!is.null(breaks)) {
    if (!is.numeric(breaks) || anyNA(breaks) || any(breaks != round(breaks)) || any(breaks < 1) || any(breaks >= nrow(X))) {
      stop("breaks must be indices T_B with 1 <= T_B < nrow(X)", call. = FALSE)
    }
    breaks <- as.integer(sort(breaks))
  }
  opt <- c(dfl_num(tol, "tol"), dfl_int(maxdeg, "maxdeg"), dfl_int(minlen, "minlen"))
  dfl_new_design(.Call(C_dfl_r_design, "matrix", X, breaks, opt, NULL, NULL), dfl_name(name))
}

#' @rdname dfl_design
#' @export
dfl_detrend_first <- function(D) {
  dfl_check_design(D)
  dfl_new_design(.Call(C_dfl_r_design, "first", D, NULL, NULL, NULL, NULL), paste0(D$name, " (detrended first)"))
}

#' @rdname dfl_design
#' @export
dfl_regmat <- function(D, n) {
  dfl_check_design(D)
  .Call(C_dfl_r_regmat, D, dfl_int(n, "n"))
}

#' @rdname dfl_design
#' @export
dfl_empty_caches <- function() {
  .Call(C_dfl_r_empty_caches)
  invisible(NULL)
}

#' @rdname dfl_design
#' @export
print.dfl_design <- function(x, ...) {
  np <- length(x$br) - 1L
  cat("dflimit design: ", x$name, "\n", sep = "")
  cat(sprintf("  %d regressor%s, %d piece%s", x$K, if (x$K == 1L) "" else "s", np, if (np == 1L) "" else "s"))
  if (np > 1L) cat(" (", if (x$loc) "boundaries" else "breaks", " at ", paste(format(x$br[2:np], digits = 4), collapse = ", "),
                   ")", sep = "")
  cat(sprintf(", degree %d%s\n", nrow(x$C) - 1L, if (x$loc) " (local coefficients)" else ""))
  if (x$first) cat("  the statistic that detrends first\n")
  invisible(x)
}

dfl_name <- function(name) {
  if (!is.character(name) || length(name) != 1L || is.na(name)) stop("name must be a string", call. = FALSE)
  name
}

dfl_is_design <- function(D) inherits(D, "dfl_design")

dfl_check_design <- function(D) {
  if (!dfl_is_design(D)) stop("D must be a design (dfl_trend(), dfl_perron(), dfl_fdesign(), ...)", call. = FALSE)
  invisible(D)
}
