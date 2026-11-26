## Submission

This is the first submission of dflimit to CRAN.

## Test environments

* local macOS Tahoe 26.5 (aarch64-apple-darwin20), R 4.5.2 (CRAN binary), Apple clang 17.0.0:
  `R CMD check --as-cran` on the tarball built by `R CMD build`.
* the same machine, R 4.5.3 (conda-forge build), clang 23.1.2: `R CMD check --as-cran` (the same result, and a
  second NOTE of that installation only: "Skipping checking math rendering: package 'V8' unavailable").
* win-builder, R-devel (2026-11-17 r90926 ucrt) and R-release (4.6.2), x86_64-w64-mingw32: 1 NOTE each (New
  submission, with the spelling flags listed below); tests 77 s and 74 s.
* GitHub Actions, R-release on windows-latest, macos-latest and ubuntu-latest, `R CMD check --as-cran --no-manual`:
  Status OK on all three.
* the same machine, R 4.5.2: the package installed with strict compiler warnings (`-Wall -pedantic
  -Wstrict-prototypes` for C, `-Wall -pedantic` for C++); the C sources also compiled with GCC 16.2 at `-O2 -Wall
  -pedantic -Wstrict-prototypes`. Neither reports a warning that `R CMD check` counts as significant.
* the same machine, R 4.5.2: the shared library built with clang 23.1.2 and `-fsanitize=address,undefined`; the
  full test suite and all examples ran without a sanitizer report.

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Simon A. Broda <simon.broda@hslu.ch>'
  New submission

  The spelling check of the DESCRIPTION flags Bao and Broda (authors of the method), Leybourne and Rothenberg
  (authors cited), ADFmax (the name of Leybourne's statistic), GLS (generalized least squares, spelled out in the
  Description) and multithreaded; they are spelled correctly.

## Notes for the reviewers

* Threads: the C engine evaluates its quadrature rules on a pool of worker threads (POSIX threads; Windows threads on
  Windows). Under `R CMD check` the package uses two threads: the thread count is set to 2 when the package is
  loaded with `_R_CHECK_LIMIT_CORES_` or `_R_CHECK_PACKAGE_NAME_` set, and the tests set it to 2 explicitly. In the
  check above (0.5.0, 26 Nov 2026) the tests took 87 s of CPU time in 48 s elapsed, and no example had a ratio of CPU
  to elapsed time above 2.
* Run time: the whole `R CMD check --as-cran` takes about 130 s on the machine above (0.5.0); no example takes more
  than 0.8 s. Four test blocks are skipped on CRAN.
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
* `src/tables.c` (6 MB of source) holds precomputed coefficients that the engine needs; the tarball is 3.0 MB and
  the installed package 2.6 MB.
