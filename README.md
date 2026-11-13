# dflimit

The exact limiting distribution of the Dickey-Fuller t-ratio in R.

`dflimit` computes the distribution function, density, quantiles and moments of the limits of the Dickey-Fuller
t-ratios under local alternatives, phi = 1 + c/n, by the exact representations of

> Bao, Y., Broda, S. A. and Kan, R., "The Exact Limiting Distribution of the Dickey-Fuller t-Ratio and Some Related
> Statistics", working paper.

for the regression without deterministic terms, with an intercept and with an intercept and a trend; with a zero,
a fixed or a stationary initial observation; for the t-ratio centred at the unit root or at a local null value c0;
and for the GLS-detrended trend case of Elliott, Rothenberg and Stock (1996). It also has the ADFmax test of
Leybourne (1995) and the probabilities of regions {U <= g(V)} of the limit experiment (point-optimal and
weighted-average-power tests). The numbers come from the multithreaded C engine of the paper (a port of its Julia
reference implementation, checked against it to 1e-10 or better); the same functions, with the same names and
arguments, exist in MATLAB and Julia.

| function | computes |
|---|---|
| `ltratio_cdf(x, cas, c, c0, g0, method)` | distribution function P[t <= x] of the limit |
| `ltratio_pdf(x, cas, c, c0, g0, method)` | its density |
| `ltratio_icdf(p, cas, c, c0, g0, method, tol)` | quantiles (asymptotic critical values) |
| `ltratio_mom(p, cas, c, c0, g0)` | moments E[t^j], j = 1..p |
| `coef_cdf(y, cas, c, c0, g0, method)`, `coef_icdf(p, ...)` | the coefficient statistic n(phi^ - phi0) |
| `adfmax_cdf`, `adfmax_pdf`, `adfmax_icdf` | the ADFmax statistic (cases DU, DT, SU, ST) |
| `ltratio_region(bnd, "DK", c)` | P{U <= g(V)} for a line, the parabola of the t-ratio, the WAP boundary |
| `ltratio_threads(n)` | the number of worker threads of the engine |

The case `cas` is a two-letter code: `D` (fixed start, `g0`, default 0) or `S` (stationary start, c < 0),
followed by `K` (no deterministic terms, tau), `U` (intercept, tau_mu) or `T` (intercept and trend, tau_tau);
`"DG"` is the GLS-detrended trend case (DF-GLS, cbar = -13.5). `method` is `"auto"` (the faster representation at
each point), `"single"` (the line integral, Theorem 2) or `"double"` (the conditional integral, Theorem S1).

All numbers come from the C engine. `coef_cdf`, `coef_icdf` and `adfmax_icdf` have no entry point of their own
there: as in the MATLAB programs, they are short R-level wrappers (the distribution function of the t-ratio at 0
with a shifted local null value, and root searches on the engine's distribution functions). Not in the package,
because the C engine does not have them: `ltratio_lincdf` (the linear form U - aV) and `power_envelope` of the
MATLAB programs, and the boundary type `"fun"` (any function g) of `ltratio_region`.

## Example

```r
library(dflimit)
ltratio_icdf(c(0.01, 0.05, 0.10), "DU")        # asymptotic critical values of tau_mu
#> [1] -3.430353 -2.861569 -2.566791
ltratio_cdf(-2.86157, "DU", -5)                # limiting power of the 5% test at c = -5
#> [1] 0.1165136
ltratio_cdf(-2.86157, "SU", -5)                # the same with a stationary start
#> [1] 0.1249641
ltratio_cdf(0, "DK", 2, c0 = 2)                # t-ratio centred at the true value
#> [1] 0.6504356
ltratio_cdf(-2, "DU", -5, g0 = 1.5)            # fixed start
#> [1] 0.8649772
ltratio_pdf(-2, "DT", -5)                      # density
#> [1] 0.5088008
ltratio_mom(4, "DT")                           # first four moments of tau_tau
#> [1]  -2.181356   5.320676 -14.035800  39.554430
coef_cdf(-14.0936, "DU", -10)                  # power of the 5% coefficient test at c = -10
#> [1] 0.4592526
```

The engine uses all processors by default; `ltratio_threads(n)` or `options(dflimit.threads = n)` (before the
package is loaded) changes that.

## Installation

`dflimit` is an ordinary R source package with C code (C11 and the C++ of the Faddeeva package), so the usual
routes apply.

**From source** (needs a compiler: [Rtools](https://cran.r-project.org/bin/windows/Rtools/) on Windows matching your
R version, the Xcode command line tools on macOS (`xcode-select --install`), the compilers of the distribution on
Linux, e.g. `build-essential`):

```r
# install.packages("remotes")
remotes::install_github("s-broda/dflimit")       # while the repository is private: set GITHUB_PAT first
```

or, from a downloaded source tarball, `install.packages("dflimit_0.1.0.tar.gz", repos = NULL, type = "source")`
(on the command line `R CMD INSTALL dflimit_0.1.0.tar.gz`).

**From a release binary** (no compiler needed): the releases of the GitHub repository (the Releases page of
`s-broda/dflimit`) carry binary packages for Windows (`.zip`) and macOS
(`.tgz`, Apple silicon and Intel), built for the R release named in the file name (a binary package works only
with that minor version of R, e.g. 4.5.x). Download the one for your platform and

```r
install.packages("dflimit_0.1.0_R4.5_windows-x86_64.zip", repos = NULL)   # Windows
install.packages("dflimit_0.1.0_R4.5_macos-arm64.tgz", repos = NULL)      # macOS, Apple silicon
```

Linux users install from source (R has no portable Linux binary format; Posit Public Package Manager and
r-universe build such binaries for packages that are public).

## Licence and citation

MIT (see `LICENSE`); the Faddeeva package in `src/Faddeeva.cc` is by Steven G. Johnson, MIT licence (see
`inst/COPYRIGHTS`). Please cite the paper when using the package.

The C engine in `src/` is copied from the paper's repository by `tools/sync-engine.sh` (the upstream commit is in
`inst/UPSTREAM`); do not edit those files here.
