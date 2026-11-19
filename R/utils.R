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

# a logical flag (TRUE or FALSE), as an integer for the engine
dfl_flag <- function(v, name) {
  if (!is.logical(v) || length(v) != 1L || is.na(v)) stop(name, " must be TRUE or FALSE", call. = FALSE)
  as.integer(v)
}

# a real scalar (finiteness and range are checked by the engine)
dfl_num <- function(v, name) {
  if (!is.numeric(v) || length(v) != 1L || is.na(v)) stop(name, " must be a single number", call. = FALSE)
  as.double(v)
}

# the method: 0 automatic, 1 the line (single) integral, 2 the conditional (double) integral, 3 the contour integral
dfl_method <- function(method) {
  msg <- "method must be \"auto\", \"single\" (\"line\"), \"double\" (\"conditional\") or \"contour\""
  if (!is.character(method) || length(method) != 1L) stop(msg, call. = FALSE)
  m <- match(tolower(method), c("auto", "single", "line", "double", "conditional", "contour"))
  if (is.na(m)) stop(msg, call. = FALSE)
  c(0L, 1L, 1L, 2L, 2L, 3L)[m]
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

# evaluate expr giving each distinct warning once (the calls of ltratio_cdf inside coef_cdf and coef_icdf would repeat
# the warning of a brittle regime at every point)
dfl_once <- function(expr) {
  seen <- character(0)
  withCallingHandlers(expr, warning = function(w) {
    m <- conditionMessage(w)
    if (!(m %in% seen)) {
      seen <<- c(seen, m)
      warning(m, call. = FALSE)
    }
    invokeRestart("muffleWarning")
  })
}
