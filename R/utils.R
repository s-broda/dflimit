# argument conversion shared by the user-level functions (the checks of the ranges are the engine's, dfl_check_args)

# the case code: model 1, 2, 3 (no deterministic terms, intercept, intercept and trend), 4 (the GLS-detrended trend
# case DG), and the start (stat = 1: stationary)
dfl_case <- function(cas, allowed = c("DK", "DU", "DT", "SK", "SU", "ST", "DG")) {
  if (!is.character(cas) || length(cas) != 1L || is.na(cas)) {
    stop("the case must be one of ", paste(allowed, collapse = ", "), call. = FALSE)
  }
  cas <- toupper(cas)
  if (!cas %in% allowed) stop("the case must be one of ", paste(allowed, collapse = ", "), call. = FALSE)
  if (cas == "DG") return(list(case = cas, model = 4L, stat = 0L))
  list(case = cas, model = match(substr(cas, 2, 2), c("K", "U", "T")), stat = as.integer(substr(cas, 1, 1) == "S"))
}

# a real scalar (finiteness and range are checked by the engine)
dfl_num <- function(v, name) {
  if (!is.numeric(v) || length(v) != 1L || is.na(v)) stop(name, " must be a single number", call. = FALSE)
  as.double(v)
}

# the method: 0 automatic, 1 the line (single) integral, 2 the conditional (double) integral
dfl_method <- function(method) {
  if (!is.character(method) || length(method) != 1L) {
    stop("method must be \"auto\", \"single\" (\"line\") or \"double\" (\"conditional\")", call. = FALSE)
  }
  m <- match(tolower(method), c("auto", "single", "line", "double", "conditional"))
  if (is.na(m)) stop("method must be \"auto\", \"single\" (\"line\") or \"double\" (\"conditional\")", call. = FALSE)
  c(0L, 1L, 1L, 2L, 2L)[m]
}

dfl_real <- function(x, name) {
  if (!(is.numeric(x) || is.logical(x) && all(is.na(x)))) stop(name, " must be numeric", call. = FALSE)
  as.double(x)
}

# the attributes of the input (dim, dimnames, names) on the output
dfl_shape <- function(out, x) {
  if (!is.null(dim(x))) {
    dim(out) <- dim(x)
    dimnames(out) <- dimnames(x)
  } else {
    names(out) <- names(x)
  }
  out
}
