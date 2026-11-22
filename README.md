# dflimit

The exact limiting distribution of the Dickey-Fuller t-ratio in R.

`dflimit` computes the distribution function, density, quantiles and moments of the limits of the Dickey-Fuller
t-ratios under local alternatives, phi = 1 + c/n, by the exact representations of

> Bao, Y., Broda, S. A. and Kan, R., "The Exact Limiting Distribution of the Dickey-Fuller t-Ratio and Some Related
> Statistics", working paper, 2026, available at <https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7725619>.

for the regression without deterministic terms, with an intercept and with an intercept and a trend; with a zero,
a fixed or a stationary initial observation; for the t-ratio centred at the unit root or at a local null value c0;
and for the GLS-detrended trend case of Elliott, Rothenberg and Stock (1996). It also has the ADFmax test of
Leybourne (1995), the probabilities of regions {U <= g(V)} of the limit experiment (point-optimal and
weighted-average-power tests), the distribution of the linear form U - aV and the asymptotic power envelope of
Elliott, Rothenberg and Stock, and, by a second engine, the same limits for any deterministic terms (see "The
general engine" below). The numbers come from the multithreaded C engine of the paper (a port of its Julia
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
| `ltratio_lincdf(b, a, cas, c, cbar)` | P(U - aV <= b), cases DK and DG (a = NaN: P(V <= b)) |
| `power_envelope(c, cas, alpha)` | the asymptotic power envelope of Elliott, Rothenberg and Stock (cases DK, DU, DG) |
| `ltratio_threads(n)` | the number of worker threads of the engine |

The case `cas` is a two-letter code: `D` (fixed start, `g0`, default 0) or `S` (stationary start, c <= -1e-4),
followed by `K` (no deterministic terms, tau), `U` (intercept, tau_mu) or `T` (intercept and trend, tau_tau);
`"DG"` is the GLS-detrended trend case (DF-GLS, cbar = -13.5). `method` is `"auto"` (at each point the fastest
representation that one evaluation gives to well below 1e-5), `"single"` (the line integral representation of the
paper), `"double"` (the conditional integral representation) or `"contour"`
(the contour integral representation, a fixed start in cases DU, DT only). Every
value is one evaluation by one method, with a warning where the parameters lie in a regime in which one evaluation is
known to be brittle; `verified = TRUE` gives the checked computation (error estimates and fallbacks).

All numbers come from the C engine. `coef_cdf`, `coef_icdf` and `adfmax_icdf` have no entry point of their own
there: as in the MATLAB programs, they are short R-level wrappers (the distribution function of the t-ratio at 0
with a shifted local null value, and root searches on the engine's distribution functions). Not in the package,
because the C engine does not have it: the boundary type `"fun"` (any function g) of `ltratio_region`.

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
ltratio_lincdf(-0.2639, -3.5, "DK", -7)        # power of the point-optimal test against cbar = -7 at c = -7
#> [1] 0.4993816
power_envelope(c(-5, -10), "DG")               # the power envelope of the trend case
#> [1] 0.1148460 0.3100516
```

The engine uses all processors by default; `ltratio_threads(n)` or `options(dflimit.threads = n)` (before the
package is loaded) changes that.

## The general engine: any deterministic terms

A second engine computes the same limits for any piecewise smooth deterministic design: polynomial trends of any order,
breaks in the level and the slope at several known dates (Perron's models), regressors given by a function of r = t/n
or by the regressor matrix of a sample, GLS detrending, and the statistic that detrends first. It writes the joint
moment generating function of the numerator and the denominator as a small determinant and inverts it on paths around
its singularities in the first argument, with a zero, a fixed or a stationary start and a local null value c0. A
design is made by a constructor and passed in place of the case code:

| design | regressors d(r) |
|---|---|
| `dfl_trend(p)` | 1, r, ..., r^p (p = -1: none, case DK; 0: case DU; 1: case DT) |
| `dfl_perron(model, lambda)` | Perron's model `"A"` (level), `"B"` (slope) or `"C"` (both), breaks at the fractions `lambda` |
| `dfl_fdesign(f, breaks)` | a function `f(r)` returning the regressors, smooth between the breaks |
| `dfl_ltrend(p)` | the polynomial trend of order p in the orthonormal Legendre basis |
| `dfl_design_from_matrix(X, breaks)` | the regressor matrix of a sample (breaks given as indices, or detected) |
| `dfl_detrend_first(D)` | the statistic that detrends first, for the design D |

`ltratio_cdf`, `ltratio_pdf`, `ltratio_icdf`, `ltratio_mom`, `coef_cdf`, `coef_icdf`, `adfmax_cdf`, `adfmax_pdf`,
`adfmax_icdf`, `ltratio_lincdf`, `power_envelope` and `ltratio_region` take a design where they take the case code,
with the start given by `g0` (a fixed start) or `start = "stationary"` and GLS detrending by `gls = cbar`; a design
that is a classical case is that case (`dfl_trend(0)` is case DU, with `start = "stationary"` case SU, and
`dfl_trend(1)` with `gls = -13.5` case DG). Their argument `engine` chooses the engine: `"auto"` takes, for what the
classical engine accepts, the faster of the two engines at the same accuracy by a rule measured on both (for the
t-ratio point by point) and the general engine for every other design; `"classical"` and `"general"` force one of
them. **`"auto"` is the default, for a case code as for a design**: where it takes the general engine for a case code
the two engines agree to about 1e-9, so that values can move in the tenth digit against the versions before 0.5.0;
`engine = "classical"` gives the values of the classical engine, those of the earlier versions.

For the designs whose numerator has rank one or two (polynomial trends of any order, an intercept with sines and
cosines, `dfl_detrend_first()` of a continuous design that does not vanish at the origin; a zero or a stationary start
and least squares detrending) the general engine also has the line integral with the explicit kernel: the classical line
integral with P(t) and the two roots supplied by the general engine's transform. `dist_cdf()`, `dist_pdf()` and
`dist_quantile()` choose by `method`: `"auto"` (the default) takes it at `c0 = 0` for `c >= -20`, where it is 2.2 to 3.5
times faster than the inversion on paths per value and 1.8 to 2.7 times per quantile, and the paths otherwise; `"line"`
forces it (an error with the reason where it does not cover the design) and `"paths"` the paths. `ltratio_cdf()`,
`ltratio_pdf()` and `ltratio_icdf()` with a design pass `method = "line"` or `"paths"` on (for a design that is a
classical case, `"line"` is the classical line integral). The classical engine's line integral takes its inner
integrals on fixed node layouts at `c <= 0`, `c0 = 0` (8 to 20 times faster than the adaptive rules in case DK, 2 to 3
times in the other cases).

The general engine's own functions, with all its options (`c`, `c0`, `start`, `gamma0`, `gls`, `verified`, `method`), are
`dist_cdf`, `dist_pdf`, `dist_quantile`, `dist_moments` (the t-ratio), `coef_pdf`, `coef_quantile`, `coef_cdf_gp`
(the coefficient statistic), `local_power`, `lin_cdf`, `lin_pdf`, `lin_quantile`, `v_cdf`, `v_pdf` (the linear form
U - aV and the law of V), `adfmax_quantile` and `region_prob`. The companions of the t-ratio for any design are the
point-optimal statistic of Elliott, Rothenberg and Stock and the Gaussian power envelope (`pt_cdf`, `pt_pdf`,
`pt_quantile`, `pt_power`, `pt_envelope`), the LM statistic of Schmidt and Phillips with breaks (`lm_cdf`,
`lm_quantile`, `lm_power`) and MZ_alpha and MSB with GLS detrending (`mza_cdf`, `mza_quantile`, `msb_cdf`,
`msb_quantile`).

```r
D <- dfl_perron("C", c(0.3, 0.7))              # breaks in level and slope at 0.3 and 0.7
ltratio_icdf(c(0.01, 0.05, 0.10), D)           # its critical values
#> [1] -5.287359 -4.750570 -4.468876
ltratio_cdf(-4, D, -10, start = "stationary")  # power at c = -10 with a stationary start
#> [1] 0.3873863
ltratio_cdf(-2.86157, "DU", -5, engine = "general")   # a classical case by the general engine
#> [1] 0.1165136
local_power(c(-5, -10, -20), dfl_perron("C", 0.5))    # the power of the 5% test
#> [1] 0.06689364 0.11488561 0.34595669
pt_envelope(c(-5, -10, -20), dfl_perron("C", 0.5))    # the power envelope of the design
#> [1] 0.08062918 0.17888854 0.60615345
lm_quantile(0.05, dfl_perron("B", 0.5))        # the LM test with a slope break
#> [1] -3.647478
```

By default every value is one evaluation (`verified = TRUE` evaluates it again with other rules on a steeper path). The
points of one call share the paths of the inversion, so a vector of points costs much less in one call than point by
point, and the engine keeps up to 64 prepared designs and paths between calls (`dfl_empty_caches()` frees them). The
general engine is checked for c from -480 to 10; designs whose roots run far out on the vertical line (above the unit
root, and some unusual designs) can need several gigabytes of memory, and where one point of the outer contour would
need more than a million nodes it stops with an error. It refuses, with an error, what it does not compute accurately
(the LM statistic above c = 10, the point-optimal statistic far above the unit root, where the grid of its inversion
does not resolve the transform). Its values are checked against those of the
Julia engine (DFLimit.jl) to 1e-9 or better (`tests/testthat/test-general.R`).

## Installation

`dflimit` is an ordinary R source package with C code (C11 and the C++ of the Faddeeva package), so the usual
routes apply.

**From source** (needs a compiler: [Rtools](https://cran.r-project.org/bin/windows/Rtools/) on Windows matching your
R version, the Xcode command line tools on macOS (`xcode-select --install`), the compilers of the distribution on
Linux, e.g. `build-essential`):

```r
# install.packages("remotes")
remotes::install_github("s-broda/dflimit")
```

or, from a downloaded source tarball, `install.packages("dflimit_0.5.0.tar.gz", repos = NULL, type = "source")`
(on the command line `R CMD INSTALL dflimit_0.5.0.tar.gz`).

**From a release binary** (no compiler needed): the releases of the GitHub repository (the Releases page of
`s-broda/dflimit`) carry binary packages for Windows (`.zip`) and macOS
(`.tgz`, Apple silicon and Intel), built for the R release named in the release notes and in the names of the
macOS files (a binary package works only with that minor version of R, e.g. 4.6.x). Download the one for your
platform and

```r
install.packages("dflimit_0.5.0.zip", repos = NULL)                       # Windows (do not rename the file)
install.packages("dflimit_0.5.0_R4.6_macos-arm64.tgz", repos = NULL)      # macOS, Apple silicon
```

Linux users install from source (R has no portable Linux binary format; Posit Public Package Manager and
r-universe build such binaries for packages that are public).

## Licence and citation

MIT (see `LICENSE`); the Faddeeva package in `src/Faddeeva.cc` is by Steven G. Johnson, MIT licence (see
`inst/COPYRIGHTS`). Please cite the paper when using the package.

The C engine in `src/` is copied from dflimit-matlab, the repository of the MATLAB package, which holds the C engine
(its directory `c/`), by `tools/sync-engine.sh` (the upstream commit is in `inst/UPSTREAM`); do not edit those files
here. The engine was developed in the paper's repository up to commit b484d4f.
