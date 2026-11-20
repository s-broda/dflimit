# the stored reference values of the Julia engine (subsets of c/test/*_ref.txt of dflimit-matlab; see the first line of
# each file in ref/)
read_ref <- function(name, cols) {
  utils::read.table(test_path("ref", name), col.names = cols, comment.char = "#")
}
case_code <- function(m, st) if (m == 4) "DG" else paste0(if (st) "S" else "D", c("K", "U", "T")[m])

# the criterion of the C engine's own tests: F to tolF absolute; f to relf relative or absf absolute (Inf matching Inf)
expect_close_Ff <- function(F, f, Fr, fr, tolF, relf, absf, label) {
  okF <- abs(F - Fr) < tolF
  okf <- abs(f - fr) / pmax(abs(fr), 1e-300) < relf || abs(f - fr) < absf || (is.infinite(f) && is.infinite(fr))
  expect(isTRUE(okF), sprintf("%s: F = %.15g, reference %.15g", label, F, Fr))
  expect(isTRUE(okf), sprintf("%s: f = %.15g, reference %.15g", label, f, fr))
}
