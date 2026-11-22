/* dflimit_r.h: the helpers shared by the two files of the R gateway, dflimit_r.c (the classical engine) and
 * dflimit_r_det.c (the general engine and the router over the two engines).  Hand-written, as the gateway files; the
 * engine's files in src/ are copied from c/ of dflimit-matlab by tools/sync-engine.sh. */
#ifndef DFLIMIT_R_H
#define DFLIMIT_R_H
#include "dflimit.h"
#include <string.h>
#include <limits.h>
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

#define ERRLEN 512

static inline double num1(SEXP a, const char *name)
{
    if (!Rf_isReal(a) || XLENGTH(a) != 1) Rf_error("%s must be a single number", name);
    return REAL(a)[0];
}
static inline int int1(SEXP a, const char *name)
{
    if (!Rf_isInteger(a) || XLENGTH(a) != 1 || INTEGER(a)[0] == NA_INTEGER) Rf_error("%s must be a single integer", name);
    return INTEGER(a)[0];
}
static inline int len(SEXP a, const char *name)
{
    if (!Rf_isReal(a)) Rf_error("%s must be a double vector", name);
    if (XLENGTH(a) > INT_MAX / 2) Rf_error("%s is too long", name);
    return (int) XLENGTH(a);
}
/* (no call in the message: the R functions check their arguments and report with call. = FALSE as well) */
static inline void engine_error(int rc, const char *err) { if (rc) Rf_errorcall(R_NilValue, "%s", err[0] ? err : "the engine failed"); }

/* the entry points of dflimit_r_det.c, registered in R_init_dflimit (dflimit_r.c) */
SEXP dfl_r_design(SEXP kind, SEXP a1, SEXP a2, SEXP a3, SEXP a4, SEXP a5);
SEXP dfl_r_regmat(SEXP D, SEXP n);
SEXP dfl_r_det(SEXP what, SEXP D, SEXP x, SEXP args, SEXP opts);
SEXP dfl_r_route(SEXP what, SEXP x, SEXP spec, SEXP args, SEXP engine);
SEXP dfl_r_check(SEXP model, SEXP c, SEXP c0, SEXP g0, SEXP stat);
SEXP dfl_r_empty_caches(void);
#endif
