## Submission

This is dflimit 0.5.0. Version 0.4.0 was submitted to CRAN on 20 November 2026 as a first submission; the incoming
check still reports a new submission, and 0.5.0 replaces it. Changes since 0.4.0 (NEWS.md has the details):

* the general engine of the lead paper (any piecewise smooth deterministic terms: polynomial trends of any order,
  breaks in level and slope, GLS detrending, the statistic that detrends first), with designs (`dfl_trend()`,
  `dfl_perron()`, `dfl_fdesign()`, ...), its own functions (`dist_cdf()`, `region_prob()`, `pt_cdf()`, `lm_cdf()`,
  ...), and the argument `engine` of the functions of the classical engine, whose default takes the faster or the
  more accurate engine;
* three extensions of the limit: volatility and persistence regimes (`dfl_volpath()`, the argument `vol`), a normal
  component (`rho`, `cadf_c0()`, `predictive_c0()`) and the first seasonal t-ratio (`hegy_t3_cdf()`, ...);
* the devices of Bao and Kan for the half-plane statistics and the moments, and the thinner rules of one evaluation;
* repairs of the C engine (wrong values with a fixed start, unbounded memory near a zero local null value and at large x
  above the unit root, the linear form for small |a|, NaN or a hang at some finite points, the line integral of the
  general engine far above the unit root at large x); the engine is that of dflimit-matlab 1914e58.

## Test environments

* local macOS Tahoe 26.5 (aarch64-apple-darwin20), R 4.5.2 (CRAN binary), Apple clang 17.0.0: `R CMD check
  --as-cran` on the tarball of 0.5.0 built by `R CMD build` (27 November 2026); the tests took 94 s of CPU time in 54 s
  elapsed.
* The checks on win-builder and GitHub Actions and the builds with strict compiler warnings and with the sanitizers
  were run for 0.4.0 (20 November 2026) and have not been repeated for 0.5.0.

## R CMD check results

0 errors | 0 warnings | 3 notes

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Simon A. Broda <simon.broda@hslu.ch>'
  New submission
  Found the following (possibly) invalid URLs: the two SSRN addresses (status 403; see below).

* checking for future file timestamps ... NOTE
  unable to verify current time (a property of the local machine's network at the time of the check)

* checking HTML version of manual ... NOTE
  Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy (the local installation).

## Notes for the reviewers

* Threads: the C engine evaluates its quadrature rules on a pool of worker threads (POSIX threads; Windows threads on
  Windows). Under `R CMD check` the package uses two threads: the thread count is set to 2 when the package is
  loaded with `_R_CHECK_LIMIT_CORES_` or `_R_CHECK_PACKAGE_NAME_` set, and the tests set it to 2 explicitly. In the
  check above (0.5.0, 27 Nov 2026) the tests took 94 s of CPU time in 54 s elapsed.
* Run time: no example takes more than 0.8 s (user time; the longest, `mza_cdf`, 0.79 s user and 0.52 s elapsed). Four test blocks are skipped on CRAN.
* The method is that of two working papers by Bao, Broda and Kan (2026), "The Exact Limiting Distribution of the
  Dickey-Fuller t-Ratio with General Deterministic Terms", on SSRN at
  <https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7745921>, and its companion, "Low-Rank Representations of
  the Classical Dickey-Fuller Limits, with an Extension to Predictive Regressions", on SSRN at
  <https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7725619>; the DESCRIPTION and two help pages give these
  addresses. SSRN answers automated requests with status 403, so the
  URL check may list them as possibly invalid; they open in a browser. Their DOIs were not yet registered when this
  was submitted.
* `src/Faddeeva.cc` and `src/Faddeeva.h` are the Faddeeva package by Steven G. Johnson (MIT licence, copyright
  Massachusetts Institute of Technology); Authors@R lists him as a contributor and MIT as a copyright holder, and
  `inst/COPYRIGHTS` reproduces the licence.
* `src/tables.c` (6 MB of source) holds precomputed coefficients that the engine needs; the tarball of 0.5.0 is
  3.4 MB and the installed package 3.5 MB.
