# tools/det-ref-subset.R: copy the reference values of the general engine and of the router, the subset of
# c/test/det_ref.txt of dflimit-matlab that the tests of this package use, into tests/testthat/ref/det_ref.txt.
#
#   Rscript tools/det-ref-subset.R /path/to/dflimit-matlab
#
# c/test/det_ref.txt holds the values of the Julia engine (DFLimit.jl, written by c/test/det_ref.jl) that the C test
# c/dfl_test_det.c checks. Kept here: the lines 'design' (the designs, which the tests build with the constructors of
# the package), 'val' (the user-level functions of the general engine) and 'route' (the router over the two engines,
# with engine auto, classical and general). Left out: the layers of the engine that the R interface does not reach
# (the prepared designs 'prep' and 'model', the determinant 'dval', the roots 'roots', 'tstar', the outer paths 'path')
# and the rule of the router 'auto' (its choice of engine, which the values of 'route' cover). The lines are copied
# unchanged, in their order; the first line of the output names the source and its commit. Only reads the checkout.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: Rscript tools/det-ref-subset.R /path/to/dflimit-matlab")
up <- normalizePath(args[1])
src <- file.path(up, "c", "test", "det_ref.txt")
if (!file.exists(src)) stop(up, " does not look like a dflimit-matlab checkout (no c/test/det_ref.txt)")
if (!file.exists("DESCRIPTION")) stop("run from the package root")
l <- readLines(src)
tag <- sub(" .*", "", l)
keep <- tag %in% c("design", "val", "route")
commit <- tryCatch(system2("git", c("--no-optional-locks", "-C", shQuote(up), "log", "-1", "--format=%h", "--",
                                    "c/test/det_ref.txt"), stdout = TRUE), error = function(e) "unknown")
head <- sprintf(paste("# the lines design, val and route (%d of %d lines) of c/test/det_ref.txt of dflimit-matlab (%s; Julia",
                      "reference values of c/test/det_ref.jl), copied by tools/det-ref-subset.R; the format is described",
                      "in c/test/det_ref.jl"), sum(keep), length(l), commit)
out <- file.path("tests", "testthat", "ref", "det_ref.txt")
writeLines(c(head, l[keep]), out)
cat(sprintf("wrote %s: %d lines (design %d, val %d, route %d) from %s at %s\n", out, sum(keep), sum(tag == "design"),
            sum(tag == "val"), sum(tag == "route"), src, commit))
