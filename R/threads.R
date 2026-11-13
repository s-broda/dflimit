#' Number of worker threads of the engine
#'
#' The C engine computes the points of its quadrature rules on a pool of worker threads (POSIX threads, or Windows
#' threads on Windows). `ltratio_threads()` returns the number of workers; `ltratio_threads(n)` sets it.
#'
#' When the package is loaded the number is taken, in this order, from the option `dflimit.threads`, from the
#' environment variable `DFL_NTHREADS`, as 2 when the environment variable `_R_CHECK_LIMIT_CORES_` is set (as by
#' `R CMD check --as-cran`), and otherwise it is the number of processors, capped by the environment variable
#' `OMP_THREAD_LIMIT` if that is set.
#'
#' @param n `NULL` (the default) to query, or the number of worker threads, a positive integer (1: no threads).
#' @return The number of worker threads in effect, invisibly when `n` is given.
#' @examples
#' old <- ltratio_threads()
#' ltratio_threads(2)
#' ltratio_threads(old)
#' @export
ltratio_threads <- function(n = NULL) {
  if (is.null(n)) return(.Call(C_dfl_r_threads, NULL))
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n != round(n) || n < 1) {
    stop("n must be a positive integer", call. = FALSE)
  }
  invisible(.Call(C_dfl_r_threads, as.integer(n)))
}

dfl_default_threads <- function() {
  opt <- getOption("dflimit.threads")
  if (!is.null(opt)) return(opt)
  if (nzchar(Sys.getenv("DFL_NTHREADS"))) return(NULL)                 # (the engine reads it)
  lim <- Sys.getenv("_R_CHECK_LIMIT_CORES_")
  if (nzchar(lim) && !tolower(lim) %in% c("false", "no", "0")) return(2L)
  omp <- suppressWarnings(as.integer(Sys.getenv("OMP_THREAD_LIMIT")))
  if (!is.na(omp) && omp >= 1L) return(min(omp, .Call(C_dfl_r_threads, NULL)))
  NULL
}
