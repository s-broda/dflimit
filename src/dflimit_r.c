/* dflimit_r.c: the R gateway to the C engine (the other files of src/ are copied from DF-Limit's c/ by
 * tools/sync-engine.sh and are not edited here).  Plain .Call entry points, registered in R_init_dflimit; the R
 * functions in R/ check and convert the arguments (the case code to the model and the start, the method name to its
 * number), this file checks their types again and calls the engine.
 *
 * The engine runs its own worker threads (dfl_threads.c); no R API is called while it runs: the inputs are copied
 * and the outputs allocated before the call, and an error of the engine (a return code with a message in err[])
 * becomes an R error on the main thread after the call has returned.  Scratch space is taken with R_alloc, which R
 * releases also when an error is raised. */
#include "dflimit.h"
#include <string.h>
#include <limits.h>
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#define ERRLEN 512

static double num1(SEXP a, const char *name)
{
    if (!Rf_isReal(a) || XLENGTH(a) != 1) Rf_error("%s must be a single number", name);
    return REAL(a)[0];
}
static int int1(SEXP a, const char *name)
{
    if (!Rf_isInteger(a) || XLENGTH(a) != 1 || INTEGER(a)[0] == NA_INTEGER) Rf_error("%s must be a single integer", name);
    return INTEGER(a)[0];
}
static int len(SEXP a, const char *name)
{
    if (!Rf_isReal(a)) Rf_error("%s must be a double vector", name);
    if (XLENGTH(a) > INT_MAX / 2) Rf_error("%s is too long", name);
    return (int) XLENGTH(a);
}
/* (no call in the message: the R functions check their arguments and report with call. = FALSE as well) */
static void engine_error(int rc, const char *err) { if (rc) Rf_errorcall(R_NilValue, "%s", err[0] ? err : "the engine failed"); }

/* .Call("dfl_r_dist", kind, x, model, c, c0, g0, stat, method, verified): list(value, nstatus, the brittle-regime message
 * of the one-evaluation default, "" if none) */
static SEXP dfl_r_dist(SEXP kind, SEXP x, SEXP model, SEXP c, SEXP c0, SEXP g0, SEXP stat, SEXP method, SEXP verified)
{
    int n = len(x, "x"), ns = 0, k = int1(kind, "kind"), m = int1(model, "model"), st = int1(stat, "stat"), me = int1(method, "method"), ve = int1(verified, "verified");
    double cc = num1(c, "c"), cc0 = num1(c0, "c0"), gg0 = num1(g0, "g0");
    SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
    char err[ERRLEN] = "";
    dfl_double_opts dv; dfl_double_defaults(&dv); dv.verified = ve != 0;
    int rc = dfl_dist(k, n, REAL(x), m, cc, cc0, gg0, st != 0, me, NULL, &dv, REAL(out), &ns, err, sizeof err);
    engine_error(rc, err);
    const char *bm = ve ? NULL : dfl_brittle(k, m, cc, cc0, gg0, st != 0, me);
    SEXP res = PROTECT(Rf_allocVector(VECSXP, 3));
    SET_VECTOR_ELT(res, 0, out);
    SET_VECTOR_ELT(res, 1, Rf_ScalarInteger(ns));
    SET_VECTOR_ELT(res, 2, Rf_mkString(bm ? bm : ""));
    UNPROTECT(2);
    return res;
}

/* .Call("dfl_r_icdf", p, model, c, c0, g0, stat, method, tol, verified): list(q, F, iter, conv, message) */
static SEXP dfl_r_icdf(SEXP p, SEXP model, SEXP c, SEXP c0, SEXP g0, SEXP stat, SEXP method, SEXP tol, SEXP verified)
{
    int n = len(p, "p"), m = int1(model, "model"), st = int1(stat, "stat"), me = int1(method, "method"), ve = int1(verified, "verified");
    double cc = num1(c, "c"), cc0 = num1(c0, "c0"), gg0 = num1(g0, "g0"), tl = num1(tol, "tol");
    SEXP q = PROTECT(Rf_allocVector(REALSXP, n)), F = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP it = PROTECT(Rf_allocVector(INTSXP, n)), cv = PROTECT(Rf_allocVector(LGLSXP, n));
    char err[ERRLEN] = "";
    dfl_double_opts dv; dfl_double_defaults(&dv); dv.verified = ve != 0;
    int rc = dfl_icdf(n, REAL(p), m, cc, cc0, gg0, st != 0, me, tl, NULL, &dv, REAL(q), REAL(F), INTEGER(it), LOGICAL(cv), err, sizeof err);
    engine_error(rc, err);
    const char *bm = ve ? NULL : dfl_brittle(2, m, cc, cc0, gg0, st != 0, me);
    SEXP res = PROTECT(Rf_allocVector(VECSXP, 5));
    SET_VECTOR_ELT(res, 0, q); SET_VECTOR_ELT(res, 1, F); SET_VECTOR_ELT(res, 2, it); SET_VECTOR_ELT(res, 3, cv);
    SET_VECTOR_ELT(res, 4, Rf_mkString(bm ? bm : ""));
    UNPROTECT(5);
    return res;
}

/* .Call("dfl_r_moments", pmax, model, c, c0, g0, stat): list(m, err) */
static SEXP dfl_r_moments(SEXP pmax, SEXP model, SEXP c, SEXP c0, SEXP g0, SEXP stat)
{
    int p = int1(pmax, "p"), m = int1(model, "model"), st = int1(stat, "stat");
    if (p < 1) Rf_error("p must be a positive integer");
    double cc = num1(c, "c"), cc0 = num1(c0, "c0"), gg0 = num1(g0, "g0");
    SEXP mo = PROTECT(Rf_allocVector(REALSXP, p)), e = PROTECT(Rf_allocVector(REALSXP, p));
    char err[ERRLEN] = "";
    int rc = dfl_moments(m, cc, cc0, gg0, st != 0, p, REAL(mo), REAL(e), err, sizeof err);
    engine_error(rc, err);
    SEXP res = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(res, 0, mo); SET_VECTOR_ELT(res, 1, e);
    UNPROTECT(3);
    return res;
}

/* .Call("dfl_r_adfmax", kind, x, model, c, stat): list(value, err, status) */
static SEXP dfl_r_adfmax(SEXP kind, SEXP x, SEXP model, SEXP c, SEXP stat)
{
    int n = len(x, "x"), k = int1(kind, "kind"), m = int1(model, "model"), st = int1(stat, "stat");
    double cc = num1(c, "c");
    SEXP out = PROTECT(Rf_allocVector(REALSXP, n)), er = PROTECT(Rf_allocVector(REALSXP, n)), su = PROTECT(Rf_allocVector(INTSXP, n));
    dfl_double_info *info = (dfl_double_info *) R_alloc((size_t) n + 1, sizeof(dfl_double_info));
    char err[ERRLEN] = "";
    int rc = dfl_adfmax(k, n, REAL(x), m, cc, st != 0, NULL, REAL(out), info, err, sizeof err);
    engine_error(rc, err);
    for (int j = 0; j < n; j++) { REAL(er)[j] = info[j].err; INTEGER(su)[j] = info[j].status; }
    SEXP res = PROTECT(Rf_allocVector(VECSXP, 3));
    SET_VECTOR_ELT(res, 0, out); SET_VECTOR_ELT(res, 1, er); SET_VECTOR_ELT(res, 2, su);
    UNPROTECT(4);
    return res;
}

/* .Call("dfl_r_region", type, par, v0, brks, c): list(P, err, status) */
static SEXP dfl_r_region(SEXP type, SEXP par, SEXP v0, SEXP brks, SEXP c)
{
    int t = int1(type, "type"), nb = len(brks, "brks");
    if (len(par, "par") != 2) Rf_error("par must have two elements");
    double vv0 = num1(v0, "v0"), cc = num1(c, "c"), P = NA_REAL;
    dfl_double_info info;
    char err[ERRLEN] = "";
    int rc = dfl_region(t, REAL(par), vv0, nb, REAL(brks), cc, NULL, &P, &info, err, sizeof err);
    engine_error(rc, err);
    SEXP res = PROTECT(Rf_allocVector(VECSXP, 3));
    SET_VECTOR_ELT(res, 0, Rf_ScalarReal(P)); SET_VECTOR_ELT(res, 1, Rf_ScalarReal(info.err)); SET_VECTOR_ELT(res, 2, Rf_ScalarInteger(info.status));
    UNPROTECT(1);
    return res;
}

/* .Call("dfl_r_lincdf", b, a, model, c, cbar): P(U - aV <= b) (a = NaN: P(V <= b)) */
static SEXP dfl_r_lincdf(SEXP b, SEXP a, SEXP model, SEXP c, SEXP cbar)
{
    int n = len(b, "b"), m = int1(model, "model");
    double aa = num1(a, "a"), cc = num1(c, "c"), cb = num1(cbar, "cbar");
    SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
    char err[ERRLEN] = "";
    int rc = dfl_lincdf(n, REAL(b), aa, m, cc, cb, NULL, REAL(out), err, sizeof err);
    engine_error(rc, err);
    UNPROTECT(1);
    return out;
}

/* .Call("dfl_r_envelope", c, model, alpha): the power envelope at the points c */
static SEXP dfl_r_envelope(SEXP c, SEXP model, SEXP alpha)
{
    int n = len(c, "c"), m = int1(model, "model");
    double al = num1(alpha, "alpha");
    SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
    char err[ERRLEN] = "";
    int rc = dfl_power_envelope(n, REAL(c), m, al, NULL, REAL(out), err, sizeof err);
    engine_error(rc, err);
    UNPROTECT(1);
    return out;
}

/* .Call("dfl_r_threads", n): the number of worker threads, after setting it to n (NULL: unchanged) */
static SEXP dfl_r_threads(SEXP n)
{
    if (!Rf_isNull(n)) { int k = int1(n, "n"); if (k < 1) Rf_error("n must be a positive integer"); dfl_set_nworkers(k); }
    return Rf_ScalarInteger(dfl_nworkers());
}

static const R_CallMethodDef call_methods[] = {
    {"dfl_r_dist", (DL_FUNC) &dfl_r_dist, 9},
    {"dfl_r_icdf", (DL_FUNC) &dfl_r_icdf, 9},
    {"dfl_r_moments", (DL_FUNC) &dfl_r_moments, 6},
    {"dfl_r_adfmax", (DL_FUNC) &dfl_r_adfmax, 5},
    {"dfl_r_region", (DL_FUNC) &dfl_r_region, 5},
    {"dfl_r_lincdf", (DL_FUNC) &dfl_r_lincdf, 5},
    {"dfl_r_envelope", (DL_FUNC) &dfl_r_envelope, 3},
    {"dfl_r_threads", (DL_FUNC) &dfl_r_threads, 1},
    {NULL, NULL, 0}
};

void R_init_dflimit(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
