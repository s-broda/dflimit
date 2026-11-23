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
| `coef_cdf(y, cas, c, c0, g0, method)`, `coef_icdf(p, ...)`, `coef_pdf(y, cas, ...)`, `coef_mom(p, cas, ...)` | the coefficient statistic n(phi^ - phi0) |
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

All numbers come from the C engine. The half-plane statistics of the six classical cases (the coefficient statistic,
the linear form, the power envelope) and the moments of the t-ratio come by default from the devices of Bao and Kan
(see "The devices of Bao and Kan" below); `adfmax_icdf` and, with `method = "engine"` or a representation of the
t-ratio, `coef_cdf` and `coef_icdf` are short R-level wrappers (the distribution function of the t-ratio at 0 with a
shifted local null value, and root searches on the engine's distribution functions). Not in the package, because the C
engine does not have it: the boundary type `"fun"` (any function g) of `ltratio_region`.

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

By default every value is one evaluation (`verified = TRUE` evaluates it again with other rules on a steeper path), by
the thinner rules of one evaluation where they hold (`rule = "auto"` of `dist_cdf()`, `dist_pdf()`, `dist_quantile()`: the
outer step 1.5 times as long and 16 and 32 nodes on the paths around the roots, by regions of the start, the deterministic
terms, `c` and `c0`, as in DFLimit.jl; `rule = "present"` the rule before November 2026). The
points of one call share the paths of the inversion, so a vector of points costs much less in one call than point by
point, and the engine keeps up to 64 prepared designs and paths between calls (`dfl_empty_caches()` frees them). The
general engine is checked for c from -480 to 10, and above the unit root, for the coefficient statistic, the moments and
the t-ratio centred at the true value (c0 = c), up to c = 16 with deterministic terms (10 with them and a fixed start)
and c = 20 without them; above those bounds every function of the general engine stops with an error that names the
statistic (the roots of its transform are no longer resolved in double precision there), and `engine = "auto"` takes the
classical engine for the classical cases. Designs whose roots run far out on the vertical line (above the unit root, and
some unusual designs) can need several gigabytes of memory, and where one point of the outer contour would need more
than a million nodes it stops with an error. It refuses, with an error, what it does not compute accurately (the LM
statistic above c = 10, the point-optimal statistic on its fixed grid far above the unit root, where the grid does not
resolve the transform, and `coef_cdf_gp()` likewise, from about c = 10 with the default grid; the message names the grid
size `N` that serves). Its values are checked against those of the Julia engine (DFLimit.jl) to 1e-9 or better
(`tests/testthat/test-general.R`).

Above the unit root the t-ratio with `c0 != c` has its upper tail at x of the order (c - c0) e^c/(2c); the classical
engine computes there by the conditional integral (within 1.3e-8 of an independent computation in case DK and 5e-7 of the
large-c law in cases DU, DT up to c = 14) and stops with an error that names the case, c, c0, x and the region computed
above c = 14. Every distribution
function of both engines ends with a check of its values: outside [0, 1] by more than 1e-8 an error that names the call.

## The devices of Bao and Kan

The half-plane statistics of the six classical cases and the moments of the t-ratio are computed by default with the
numerical devices of Yong Bao and Raymond Kan, "Exact Finite-Sample Distribution and Moments of the Least Squares
Estimator in First-Order Autoregressions: O(1) Algorithms", as in the Julia package (its module `DFKan`) and the MATLAB
programs (the C engine's `dfl_kan.c` and `dfl_det_kan.c`). The joint characteristic function of the limit pair comes in
closed forms that are free of cancellation on both sides of the unit root, with the branch of the square root fixed at
each point, and in Taylor-mode arithmetic (jets), so that its derivatives come out of the same assembly. A statistic that
rejects on a half-plane is one inversion along a ray, by adaptive quadrature over decades with the tail completed by
repeated integration by parts; the moments of the t-ratio are one-dimensional integrals of the j-th derivative of the
moment generating function, which one jet of order j gives at each node.

| function | the default | warm time in C against `method = "engine"` |
|---|---|---|
| `coef_cdf`, `coef_icdf`, `coef_pdf`, `coef_mom` | the cases DK, DU, DT (zero or fixed start), SK, SU, ST, c <= 20 | 0.005 to 0.18 (quantiles 0.02; the density 0.03 against the general engine's) |
| `ltratio_lincdf` | case DK, a != 0, -100 <= c <= 4; the six cases through their designs, -100 <= c <= 2 | 0.07 to 0.24 |
| `power_envelope` | case DK (DU), c <= -1e-4 | 0.13 to 0.16 |
| `ltratio_mom` | the six cases, any c0 and fixed start, c <= 20, up to 31 moments | 0.013 to 0.23 |
| `pt_cdf`, `pt_pdf`, `pt_quantile`, `pt_power`, `pt_envelope` | any design, c <= 4 | 0.08 to 0.15 (quantiles 0.3 to 0.4, envelopes 0.2) |
| `coef_cdf_gp` | zero or stationary start, c <= 5 | 0.05 to 0.1 |

`method = "inversion"` (`"jets"` for the moments) forces them, `method = "engine"` is the computation of before, which
reproduces those values (`method = "grid"` the fixed grid of the point-optimal functions and of `coef_cdf_gp`, also taken
where `N` or `umax` is given); `tol` is the tolerance (1e-12; 1e-10 relative for the moments), `verified = TRUE` recomputes
by the engines the values whose error estimate exceeds it, `err = TRUE` returns the estimates with the values (`info =
TRUE` for the moments), and a warning names the values whose estimate exceeds the tolerance. The regions, the tolerances
and the warnings are those of the Julia package, whose README lists the comparisons behind them and the cases where the
engines' one evaluation was wrong (the coefficient statistic of case DU with a fixed start far from the unit root, the
linear form and the moments with a fixed start above the unit root). The values agree with DFLimit.jl to 6e-15 in the
single inversions and the moments (`tests/testthat/test-kan.R`). The same devices fix the branches inside the engines: P(t),
M_0(t) and the labels of the roots of the classical line and contour integrals and the argument of log P(t) in the general
engine's line integral are determined at each point instead of being continued along the grids, to 4e-13 the values of
before, at the same cost.

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
