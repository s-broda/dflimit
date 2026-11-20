#' dflimit: Exact limiting distribution of the Dickey-Fuller t-ratio
#'
#' The distribution function, density, quantiles and moments of the limits of the Dickey-Fuller t-ratios under
#' local alternatives, by the exact representations of Bao, Broda and Kan, computed by their C engine (a port of
#' the Julia reference implementation of the paper; the same functions, with the same names and arguments, exist in
#' MATLAB and Julia).
#'
#' \describe{
#'   \item{[ltratio_cdf()], [ltratio_pdf()], [ltratio_icdf()], [ltratio_mom()]}{the t-ratio, cases DK, DU, DT
#'     (zero or fixed start), SK, SU, ST (stationary start) and DG (GLS-detrended trend case);}
#'   \item{[coef_cdf()], [coef_icdf()]}{the coefficient statistic \eqn{n(\hat\phi - \phi_0)}{n(phi^ - phi0)};}
#'   \item{[adfmax_cdf()], [adfmax_pdf()], [adfmax_icdf()]}{the ADFmax test of Leybourne (1995);}
#'   \item{[ltratio_region()]}{the probability of a region \eqn{\{U \le g(V)\}}{{U <= g(V)}}, case DK;}
#'   \item{[ltratio_lincdf()]}{the distribution of the linear form \eqn{U - aV}, cases DK and DG;}
#'   \item{[power_envelope()]}{the asymptotic power envelope of Elliott, Rothenberg and Stock (1996);}
#'   \item{[ltratio_threads()]}{the number of worker threads.}
#' }
#'
#' @references Bao, Y., Broda, S. A. and Kan, R., "The Exact Limiting Distribution of the Dickey-Fuller t-Ratio
#'   and Some Related Statistics", working paper, 2026, \url{https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7725619}.
#' @keywords internal
#' @useDynLib dflimit, .registration = TRUE, .fixes = "C_"
"_PACKAGE"

.onLoad <- function(libname, pkgname) {
  n <- dfl_default_threads()
  if (!is.null(n)) ltratio_threads(n)
}
