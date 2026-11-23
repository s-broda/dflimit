/* dflimit_r_det.c: the R gateway to the general engine (dfl_det_*.c: any piecewise smooth deterministic design) and to
 * the router over the two engines (dfl_engine_*, dfl_det_route.c).  Hand-written, as dflimit_r.c, whose R_init_dflimit
 * registers these entry points too.
 *
 * A design crosses the gateway as a list of numeric arrays, the fields of dfl_design (br: the np + 1 break points;
 * C: the (deg + 1) x (K np) matrix of the coefficient matrices of the pieces side by side; K; loc; first), which R owns:
 * the engine only reads a design given to it (it copies what it caches).  A design made by a constructor of the engine
 * is copied into such a list and freed at once; until then it is held by an external pointer whose finalizer frees it,
 * so that an R error raised in between (an allocation that fails) does not leak it.
 *
 * As in dflimit_r.c, no R API is called while the engine computes on its workers; the one exception is the function
 * of a design given by an R function (dfl_design_fdesign), which the engine calls on the calling thread only: it is
 * evaluated by R_tryEval, so that an R error in it does not jump over the engine's own frames (the R side catches
 * errors in the function and passes back NULL, and the design is then freed and the error raised).
 *
 * Warnings: the general engine's warning of a call (dfl_det_warning) and, for the router, the classical engine's
 * warning of the points that the line integral refused (dfl_warning) are thread-local; they are cleared here before
 * each call (dfl_det_warning_clear, dfl_warning_clear), and the R functions give them as R warnings. */
#include "dflimit_r.h"
#include <stdlib.h>
#include <math.h>

/* both warnings of the calling thread cleared before a call whose warnings are read (dfl_det_warning_clear also drops a
 * hold that an earlier call may have left) */
static void clear_warnings(void)
{
    dfl_det_warning_clear();
    dfl_warning_clear();
}
static SEXP mkstr(const char *s) { return Rf_mkString(s ? s : ""); }

/* ---- a design from R (read only) ---- */
static SEXP getel(SEXP list, const char *name)
{
    SEXP nm = Rf_getAttrib(list, R_NamesSymbol);
    if (!Rf_isString(nm)) return R_NilValue;
    for (R_xlen_t i = 0; i < XLENGTH(list); i++) if (!strcmp(CHAR(STRING_ELT(nm, i)), name)) return VECTOR_ELT(list, i);
    return R_NilValue;
}
static int flag1(SEXP a, const char *name)
{
    if (!Rf_isLogical(a) || XLENGTH(a) != 1 || LOGICAL(a)[0] == NA_LOGICAL) Rf_error("the design's %s must be TRUE or FALSE", name);
    return LOGICAL(a)[0] != 0;
}
static double ZERO1[1] = {0};
static void design_in(SEXP L, dfl_design *D)
{
    if (!Rf_isNewList(L)) Rf_error("not a design of the general engine");
    SEXP br = getel(L, "br"), C = getel(L, "C"), K = getel(L, "K");
    if (!Rf_isReal(br) || XLENGTH(br) < 2 || XLENGTH(br) > INT_MAX / 2) Rf_error("the design's break points br must be a double vector of length at least 2");
    int np = (int) XLENGTH(br) - 1;
    const double *b = REAL(br);
    if (b[0] != 0 || b[np] != 1) Rf_error("the design's break points must run from 0 to 1");
    for (int j = 0; j < np; j++) if (!(b[j + 1] > b[j])) Rf_error("the design's break points must increase");
    SEXP dm = Rf_getAttrib(C, R_DimSymbol);
    if (!Rf_isReal(C) || !Rf_isInteger(dm) || XLENGTH(dm) != 2) Rf_error("the design's coefficients C must be a double matrix");
    int rows = INTEGER(dm)[0], cols = INTEGER(dm)[1], k = int1(K, "the design's K");
    if (rows < 1 || k < 0 || cols != k * np) Rf_error("the design's coefficient matrix must have deg + 1 rows and K np columns");
    for (R_xlen_t i = 0; i < XLENGTH(C); i++) if (!R_FINITE(REAL(C)[i])) Rf_error("the design's coefficients must be finite");
    D->K = k; D->np = np; D->deg = rows - 1; D->br = REAL(br);
    D->C = XLENGTH(C) > 0 ? REAL(C) : ZERO1;
    D->loc = flag1(getel(L, "loc"), "loc"); D->first = flag1(getel(L, "first"), "first");
}
/* ---- a design made by the engine ---- */
static void design_fin(SEXP p)
{
    dfl_design *D = (dfl_design *) R_ExternalPtrAddr(p);
    if (D) { dfl_design_free(D); free(D); R_ClearExternalPtr(p); }
}
/* an empty design held by the external pointer *h (protected by the caller) */
static dfl_design *design_new(SEXP *h)
{
    *h = PROTECT(R_MakeExternalPtr(NULL, R_NilValue, R_NilValue));
    R_RegisterCFinalizerEx(*h, design_fin, TRUE);
    dfl_design *D = calloc(1, sizeof(dfl_design));
    if (!D) Rf_error("out of memory");
    R_SetExternalPtrAddr(*h, D);
    UNPROTECT(1);
    return D;
}
/* the list of a design and the warning of its constructor ("" if none; the caller frees the design afterwards) */
static SEXP design_out(const dfl_design *D, const char *warning)
{
    SEXP res = PROTECT(Rf_allocVector(VECSXP, 6)), nm = PROTECT(Rf_allocVector(STRSXP, 6));
    SEXP br = Rf_allocVector(REALSXP, D->np + 1);
    SET_VECTOR_ELT(res, 0, br);
    memcpy(REAL(br), D->br, (size_t) (D->np + 1) * sizeof(double));
    SEXP C = Rf_allocMatrix(REALSXP, D->deg + 1, D->K * D->np);
    SET_VECTOR_ELT(res, 1, C);
    if (D->K > 0) memcpy(REAL(C), D->C, (size_t) (D->deg + 1) * D->K * D->np * sizeof(double));
    SET_VECTOR_ELT(res, 2, Rf_ScalarInteger(D->K));
    SET_VECTOR_ELT(res, 3, Rf_ScalarLogical(D->loc != 0));
    SET_VECTOR_ELT(res, 4, Rf_ScalarLogical(D->first != 0));
    SET_VECTOR_ELT(res, 5, mkstr(warning));
    const char *names[6] = {"br", "C", "K", "loc", "first", "warning"};
    for (int i = 0; i < 6; i++) SET_STRING_ELT(nm, i, Rf_mkChar(names[i]));
    Rf_setAttrib(res, R_NamesSymbol, nm);
    UNPROTECT(2);
    return res;
}

/* the regressors of a design given by an R function: g(r) returns the K values, or NULL where the function failed (the
 * R side records the message) */
typedef struct { SEXP fun, rho; int K, bad; } rfun_ctx;
static void rfun(double r, double *out, void *v)
{
    rfun_ctx *m = (rfun_ctx *) v;
    for (int i = 0; i < m->K; i++) out[i] = 0;
    if (m->bad) return;
    SEXP a = PROTECT(Rf_ScalarReal(r)), call = PROTECT(Rf_lang2(m->fun, a));
    int er = 0;
    SEXP res = R_tryEval(call, m->rho, &er);
    if (er || !Rf_isReal(res) || XLENGTH(res) != m->K) m->bad = 1;
    else memcpy(out, REAL(res), (size_t) m->K * sizeof(double));
    UNPROTECT(2);
}

/* .Call("dfl_r_design", kind, a1, ..., a5): a design of the general engine as a list (see the top of this file) with a
 * sixth element, the constructor's warning ("" if none: dependent columns dropped by design_from_matrix):
 *   "trend" p | "ltrend" p | "perron" model ("A", "B", "C") lambda | "first" D |
 *   "fdesign" g rho K breaks c(tol, maxdeg, sub) (NULL where g failed) |
 *   "matrix" X TB (integer, or NULL: the breaks detected) c(tol, maxdeg, minlen) */
SEXP dfl_r_design(SEXP kind, SEXP a1, SEXP a2, SEXP a3, SEXP a4, SEXP a5)
{
    if (!Rf_isString(kind) || XLENGTH(kind) != 1) Rf_error("kind must be a string");
    const char *k = CHAR(STRING_ELT(kind, 0));
    char err[ERRLEN] = "";
    SEXP h;
    dfl_design *D = design_new(&h);
    PROTECT(h);
    int rc = 1;
    clear_warnings();
    if (!strcmp(k, "trend")) rc = dfl_design_trend(int1(a1, "p"), D, err, sizeof err);
    else if (!strcmp(k, "ltrend")) rc = dfl_design_ltrend(int1(a1, "p"), D, err, sizeof err);
    else if (!strcmp(k, "perron")) {
        if (!Rf_isString(a1) || XLENGTH(a1) != 1) Rf_error("model must be a string");
        int q = len(a2, "lambda");
        rc = dfl_design_perron(CHAR(STRING_ELT(a1, 0))[0], q, REAL(a2), D, err, sizeof err);
    } else if (!strcmp(k, "first")) {
        dfl_design B; design_in(a1, &B);
        rc = dfl_design_detrend_first(&B, D, err, sizeof err);
    } else if (!strcmp(k, "fdesign")) {
        rfun_ctx m = {a1, a2, int1(a3, "K"), 0};
        int nb = len(a4, "breaks");
        if (len(a5, "the options") != 3) Rf_error("the options of fdesign are tol, maxdeg and sub");
        const double *o = REAL(a5);
        rc = dfl_design_fdesign(rfun, &m, m.K, nb, REAL(a4), o[0], (int) o[1], (int) o[2], D, err, sizeof err);
        if (m.bad) { design_fin(h); UNPROTECT(1); return R_NilValue; }
    } else if (!strcmp(k, "matrix")) {
        SEXP dm = Rf_getAttrib(a1, R_DimSymbol);
        if (!Rf_isReal(a1) || !Rf_isInteger(dm) || XLENGTH(dm) != 2) Rf_error("X must be a double matrix");
        int n = INTEGER(dm)[0], K = INTEGER(dm)[1], nb = -1;
        const int *TB = NULL;
        if (!Rf_isNull(a2)) {
            if (!Rf_isInteger(a2)) Rf_error("the breaks must be an integer vector");
            nb = (int) XLENGTH(a2); TB = INTEGER(a2);
        }
        if (len(a3, "the options") != 3) Rf_error("the options of design_from_matrix are tol, maxdeg and minlen");
        const double *o = REAL(a3);
        rc = dfl_design_from_matrix(n, K, REAL(a1), nb, TB, o[0], (int) o[1], (int) o[2], D, err, sizeof err);
    } else Rf_error("unknown kind of design %s", k);
    if (rc) { design_fin(h); engine_error(rc, err); }
    SEXP res = PROTECT(design_out(D, dfl_det_warning()));
    design_fin(h);
    UNPROTECT(2);
    return res;
}

/* .Call("dfl_r_regmat", D, n): the n x K matrix of the regressors d(t/n), t = 1..n */
SEXP dfl_r_regmat(SEXP Ds, SEXP n)
{
    dfl_design D; design_in(Ds, &D);
    int nn = int1(n, "n");
    if (nn < 1) Rf_error("n must be a positive integer");
    SEXP X = PROTECT(Rf_allocMatrix(REALSXP, nn, D.K));
    if (D.K > 0) dfl_design_regmat(&D, nn, REAL(X));
    UNPROTECT(1);
    return X;
}

/* the options of the general engine from c(c0, start, gamma0, gls, verified), c(c0, start, gamma0, gls, verified, method)
 * or c(c0, start, gamma0, gls, verified, method, rule) (gls NaN: least squares; method 0 auto, 1 the line integral with
 * the explicit kernel, 2 the inversion on paths; rule the rule of one evaluation of the paths, DFL_RULE_*) */
static dfl_det_opts opts_in(SEXP o)
{
    int no = len(o, "the options");
    if (no < 5 || no > 7) Rf_error("the options are c0, start, gamma0, gls, verified and (optionally) method and rule");
    const double *v = REAL(o);
    dfl_det_opts d; dfl_det_defaults(&d);
    d.c0 = v[0]; d.start = (int) v[1]; d.gamma0 = v[2]; d.gls = v[3]; d.verified = v[4] != 0;
    if (no >= 6) d.method = (int) v[5];
    if (no >= 7) d.rule = (int) v[6];
    return d;
}

/* .Call("dfl_r_det", what, D, x, args, opts): a function of the general engine; list(value, extra, warning):
 *   "cdf" "pdf" x | args c                       "quantile" p | c x0 tol
 *   "coefcdf" "coefpdf" y | c                    "coefquantile" p | c y0 tol         "coefgp" y | c N umax method tol
 *   (extra: the error estimates of the single inversion of Bao and Kan, NaN on the fixed grid)
 *   "moments" (x unused) | c pmax (extra: the error estimates)
 *   "power" cs | level statistic (0 the t-ratio, 1 the coefficient statistic)
 *   "lincdf" "linpdf" b | a c (a NaN: the law of V)   "linquantile" p | a c b0 tol
 *   "adfcdf" "adfpdf" x | c                      "adfquantile" p | c x0 tol
 *   "region" (x unused) | type par1 par2 c
 *   "ptcdf" "ptpdf" q | cbar c N umax method tol (extra: the error estimates)   "ptquantile" p | cbar N umax method tol
 *   "ptpower" cs | cbar level N umax method tol (extra: the critical values)   "ptenvelope" cs | level N umax method tol
 * (method of the point-optimal functions and coef_cdf_gp: 0 automatic, the single inversion of Bao and Kan up to c = 4, 5
 * for coef_cdf_gp, and the fixed grid beyond or where N or umax is given (> 0); 5 the single inversion; 7 the fixed grid;
 * tol NaN for 1e-12)
 *   "lmcdf" x | c      "lmquantile" p | c      "lmpower" cs | level (extra: the critical value)
 *   "mzacdf" y | c gls      "mzaquantile" p | c gls y0 tol      "msbcdf" x | c gls      "msbquantile" p | c gls
 * opts = c(c0, start, gamma0, gls, verified) where the function takes the options; the probabilities are checked by the
 * R functions. */
SEXP dfl_r_det(SEXP what, SEXP Ds, SEXP x, SEXP args, SEXP opts)
{
    if (!Rf_isString(what) || XLENGTH(what) != 1) Rf_error("what must be a string");
    const char *w = CHAR(STRING_ELT(what, 0));
    dfl_design D; design_in(Ds, &D);
    int n = len(x, "x"), na = len(args, "args");
    const double *a = REAL(args), *xv = REAL(x);
    dfl_det_opts o = opts_in(opts);
#define NEED(k) do { if (na < (k)) Rf_error("%s needs %d arguments", w, (k)); } while (0)
    char err[ERRLEN] = "";
    int rc = 1, nout = n, nextra = 0;
    if (!strcmp(w, "moments")) { NEED(2); nout = (int) a[1]; nextra = nout; if (nout < 1) Rf_error("p must be a positive integer"); }
    if (!strcmp(w, "ptpower") || !strcmp(w, "ptcdf") || !strcmp(w, "ptpdf") || !strcmp(w, "ptquantile") || !strcmp(w, "coefgp")) nextra = n;
    if (!strcmp(w, "lmpower")) nextra = 1;
    SEXP out = PROTECT(Rf_allocVector(REALSXP, nout)), ex = PROTECT(Rf_allocVector(REALSXP, nextra));
    double *v = REAL(out), *e = REAL(ex);
    clear_warnings();
    if (!strcmp(w, "cdf") || !strcmp(w, "pdf")) { NEED(1); rc = dfl_det_dist(w[0] == 'c' ? 2 : 1, &D, n, xv, a[0], &o, v, err, sizeof err); }
    else if (!strcmp(w, "quantile")) { NEED(3); rc = dfl_det_quantile(&D, n, xv, a[0], a[1], a[2], &o, v, err, sizeof err); }
    else if (!strcmp(w, "coefcdf") || !strcmp(w, "coefpdf")) { NEED(1); rc = dfl_det_coef(w[4] == 'c' ? 2 : 1, &D, n, xv, a[0], &o, v, err, sizeof err); }
    else if (!strcmp(w, "coefquantile")) { NEED(3); rc = dfl_det_coef_quantile(&D, n, xv, a[0], a[1], a[2], &o, v, err, sizeof err); }
    else if (!strcmp(w, "coefgp")) { NEED(5); rc = dfl_det_coef_cdf_gp_x(&D, n, xv, a[0], &o, (int) a[3], (int) a[1], a[2], a[4], v, e, err, sizeof err); }
    else if (!strcmp(w, "moments")) rc = dfl_det_moments(&D, a[0], nout, &o, v, e, err, sizeof err);
    else if (!strcmp(w, "power")) { NEED(2); rc = dfl_det_local_power(&D, n, xv, a[0], (int) a[1], &o, v, err, sizeof err); }
    else if (!strcmp(w, "lincdf") || !strcmp(w, "linpdf")) { NEED(2); rc = dfl_det_lin(w[3] == 'c' ? 2 : 1, &D, n, xv, a[0], a[1], &o, v, err, sizeof err); }
    else if (!strcmp(w, "linquantile")) {
        NEED(4); rc = 0;
        for (int k = 0; k < n && !rc; k++) rc = dfl_det_lin_quantile(&D, xv[k], a[0], a[1], a[2], a[3], &o, &v[k], err, sizeof err);
    } else if (!strcmp(w, "adfcdf") || !strcmp(w, "adfpdf")) { NEED(1); rc = dfl_det_adfmax(w[3] == 'c' ? 2 : 1, &D, n, xv, a[0], &o, v, err, sizeof err); }
    else if (!strcmp(w, "adfquantile")) {
        NEED(3); rc = 0;
        for (int k = 0; k < n && !rc; k++) rc = dfl_det_adfmax_quantile(&D, xv[k], a[0], a[1], a[2], &o, &v[k], err, sizeof err);
    } else if (!strcmp(w, "region")) {
        NEED(4); double par[2] = {a[1], a[2]};
        if (nout != 1) Rf_error("region: one value");
        rc = dfl_det_region(&D, (int) a[0], par, a[3], &o, v, err, sizeof err);
    } else if (!strcmp(w, "ptcdf") || !strcmp(w, "ptpdf")) { NEED(6); rc = dfl_det_pt_x(w[2] == 'c' ? 2 : 1, &D, n, xv, a[0], a[1], (int) a[4], (int) a[2], a[3], a[5], v, e, err, sizeof err); }
    else if (!strcmp(w, "ptquantile")) {
        NEED(5); rc = 0;
        for (int k = 0; k < n && !rc; k++) rc = dfl_det_pt_quantile_x(&D, xv[k], a[0], (int) a[3], (int) a[1], a[2], a[4], &v[k], &e[k], err, sizeof err);
    } else if (!strcmp(w, "ptpower")) {
        NEED(6); rc = 0;
        for (int k = 0; k < n && !rc; k++) rc = dfl_det_pt_power_x(&D, a[0], xv[k], a[1], (int) a[4], (int) a[2], a[3], a[5], &e[k], &v[k], err, sizeof err);
    } else if (!strcmp(w, "ptenvelope")) { NEED(5); rc = dfl_det_pt_envelope_x(&D, n, xv, a[0], (int) a[3], (int) a[1], a[2], a[4], v, err, sizeof err); }
    else if (!strcmp(w, "lmcdf")) { NEED(1); rc = dfl_det_lm_cdf(&D, n, xv, a[0], v, err, sizeof err); }
    else if (!strcmp(w, "lmquantile")) {
        NEED(1); rc = 0;
        for (int k = 0; k < n && !rc; k++) rc = dfl_det_lm_quantile(&D, xv[k], a[0], &v[k], err, sizeof err);
    } else if (!strcmp(w, "lmpower")) { NEED(1); rc = dfl_det_lm_power(&D, n, xv, a[0], e, v, err, sizeof err); }
    else if (!strcmp(w, "mzacdf")) { NEED(2); rc = dfl_det_mza_cdf(&D, n, xv, a[0], a[1], &o, v, err, sizeof err); }
    else if (!strcmp(w, "mzaquantile")) {
        NEED(4); rc = 0;
        for (int k = 0; k < n && !rc; k++) rc = dfl_det_mza_quantile(&D, xv[k], a[0], a[1], a[2], a[3], &o, &v[k], err, sizeof err);
    } else if (!strcmp(w, "msbcdf")) { NEED(2); rc = dfl_det_msb_cdf(&D, n, xv, a[0], a[1], v, err, sizeof err); }
    else if (!strcmp(w, "msbquantile")) {
        NEED(2); rc = 0;
        for (int k = 0; k < n && !rc; k++) rc = dfl_det_msb_quantile(&D, xv[k], a[0], a[1], &v[k], err, sizeof err);
    } else Rf_error("unknown function %s of the general engine", w);
#undef NEED
    engine_error(rc, err);
    SEXP res = PROTECT(Rf_allocVector(VECSXP, 3));
    SET_VECTOR_ELT(res, 0, out); SET_VECTOR_ELT(res, 1, ex); SET_VECTOR_ELT(res, 2, mkstr(dfl_det_warning()));
    UNPROTECT(3);
    return res;
}

/* the computation of the router: spec = list(D, model, start, g0, gls), D a design list or NULL (a classical case: the
 * model 1 to 4 with its start) */
static void spec_in(SEXP s, dfl_spec *S, dfl_design *D)
{
    if (!Rf_isNewList(s) || XLENGTH(s) != 5) Rf_error("spec must be a list of five elements");
    SEXP d = VECTOR_ELT(s, 0);
    if (Rf_isNull(d)) S->D = NULL; else { design_in(d, D); S->D = D; }
    S->model = int1(VECTOR_ELT(s, 1), "model"); S->start = int1(VECTOR_ELT(s, 2), "start");
    S->g0 = num1(VECTOR_ELT(s, 3), "g0"); S->gls = num1(VECTOR_ELT(s, 4), "gls");
}

/* .Call("dfl_r_route", what, x, spec, args, engine): the router (engine 0 auto, 1 classical, 2 general; with a classical
 * case and engine auto the R functions take the classical engine's own entry points where a method other than 0 or the
 * checked computation is asked for); list(value, the warnings, a character vector):
 *   "cdf" "pdf" x | args c c0 method verified          "icdf" p | c c0 tol method
 *   "mom" (x unused) | c c0 pmax [imethod tol verified]           "coefmom" (x unused) | c c0 pmax imethod tol
 *   "coef" y | c c0 [imethod tol verified]       "coefq" p | c c0 tol [imethod verified]     "coefpdf" y | c c0 imethod tol verified
 *   "lin" b | a model c cbar [imethod tol verified] (spec unused)   "lind" b | a c imethod tol verified (the design of spec)
 *   "env" cs | model alpha [imethod tol verified] (spec unused)
 * imethod (dfl_inv_opts): 0 automatic (the single inversion of Bao and Kan, the moments from jets, in their cases and
 * regions), 5 forces them, 6 the engines' computation of before, 1-4 the representations of the coefficient statistic; tol
 * NaN for the default; the error estimates of these functions are the third element of the list (NaN where the engines
 * computed)
 *   "adf" "adfpdf" x | c               "region" (x unused) | type par1 par2 c
 * The warnings: where a computation of a classical case with engine auto takes the classical engine's one evaluation
 * in a regime in which it is known to be brittle, the message of dfl_brittle (as the classical entry points give it);
 * dfl_warning; dfl_det_warning. */
SEXP dfl_r_route(SEXP what, SEXP x, SEXP spec, SEXP args, SEXP engine)
{
    if (!Rf_isString(what) || XLENGTH(what) != 1) Rf_error("what must be a string");
    const char *w = CHAR(STRING_ELT(what, 0));
    dfl_spec S; dfl_design D; spec_in(spec, &S, &D);
    int n = len(x, "x"), na = len(args, "args"), eng = int1(engine, "engine");
    const double *a = REAL(args), *xv = REAL(x);
#define NEED(k) do { if (na < (k)) Rf_error("%s needs %d arguments", w, (k)); } while (0)
    char err[ERRLEN] = "";
    int rc = 1, nout = n;
    if (!strcmp(w, "mom") || !strcmp(w, "coefmom")) { NEED(3); nout = (int) a[2]; if (nout < 1) Rf_error("p must be a positive integer"); }
    if (!strcmp(w, "region")) nout = 1;
    SEXP out = PROTECT(Rf_allocVector(REALSXP, nout)), eo = PROTECT(Rf_allocVector(REALSXP, nout));
    double *v = REAL(out), *ev = REAL(eo);
    for (int k = 0; k < nout; k++) ev[k] = NAN;
    /* the options of the single inversion from the trailing arguments from position k0 (dfl_inv_opts) */
#define IOPTS(k0) do { dfl_inv_defaults(&io); if (na > (k0)) io.method = (int) a[k0]; if (na > (k0) + 1) io.tol = a[(k0) + 1]; if (na > (k0) + 2) io.verified = a[(k0) + 2] != 0; } while (0)
    dfl_inv_opts io;
    /* the classical case of the computation (spec->D NULL), for the warnings of its brittle regimes */
    int cls = S.D == NULL && eng == DFL_ENGINE_AUTO && S.model >= 1 && S.model <= 4, stat = S.start == 2;
    double g0 = stat ? 0.0 : S.g0, gls = S.model == 4 ? DFL_CBAR_GLS : NAN;
    const char *bm[64]; int nb = 0;
    clear_warnings();
    if (!strcmp(w, "cdf") || !strcmp(w, "pdf")) {
        NEED(4); int kind = w[0] == 'c' ? 2 : 1, me = (int) a[2], ve = a[3] != 0;
        rc = dfl_engine_dist(kind, n, xv, &S, a[0], a[1], eng, me, ve, v, err, sizeof err);
        if (!rc && cls && me == 0 && !ve) {
            /* (the points of the classical engine, as merged_dist of dfl_det_route.c chooses them) */
            int any = 0;
            for (int k = 0; k < n && !any; k++) {
                if (!isfinite(xv[k])) continue;
                if (dfl_autoengine(DFL_F_DIST, S.model, stat, g0, gls, a[0], a[1], xv[k]) == DFL_ENGINE_CLASSICAL) any = 1;
                else if (v[k] < 1e-3 && dfl_automethod(S.model, a[0], a[1], g0, stat, xv[k]) == 2) any = 1;
            }
            if (any) bm[nb] = dfl_brittle(kind, S.model, a[0], a[1], g0, stat, 0), nb += bm[nb] != NULL;
        }
    } else if (!strcmp(w, "icdf")) {
        NEED(4);
        rc = dfl_engine_icdf(n, xv, &S, a[0], a[1], eng, (int) a[3], a[2], v, err, sizeof err);
        if (!rc && cls && a[3] == 0 && dfl_autoengine(DFL_F_ICDF, S.model, stat, g0, gls, a[0], a[1], 0.0) == DFL_ENGINE_CLASSICAL) {
            bm[nb] = dfl_brittle(2, S.model, a[0], a[1], g0, stat, 0); nb += bm[nb] != NULL;
        }
    } else if (!strcmp(w, "mom")) { IOPTS(3); rc = dfl_engine_mom_x(nout, &S, a[0], a[1], eng, &io, v, ev, err, sizeof err); }
    else if (!strcmp(w, "coefmom")) { IOPTS(3); rc = dfl_engine_coef_mom(nout, &S, a[0], a[1], eng, &io, v, ev, err, sizeof err); }
    else if (!strcmp(w, "coefpdf")) { NEED(2); IOPTS(2); rc = dfl_engine_coef_pdf(n, xv, &S, a[0], a[1], eng, &io, v, ev, err, sizeof err); }
    else if (!strcmp(w, "coef") || !strcmp(w, "coefq")) {
        NEED(2);
        if (w[4] == 'q') {
            NEED(3); dfl_inv_defaults(&io); if (na > 3) io.method = (int) a[3]; if (na > 4) io.verified = a[4] != 0;
            rc = dfl_engine_coef_icdf_x(n, xv, &S, a[0], a[1], eng, a[2], &io, v, ev, err, sizeof err);
        } else { IOPTS(2); rc = dfl_engine_coef_x(n, xv, &S, a[0], a[1], eng, &io, v, ev, err, sizeof err); }
        /* (where the engines computed: the t-ratio at 0 with the null value c0 + y) */
        if (!rc && cls && !io.verified && !dfl_inv_takes(&S, a[0], a[1], eng, io.method)) {
            /* (the t-ratio at 0 with the null value c0 + y: its engine and its brittle regime at each y, or at each
             * quantile) */
            for (int k = 0; k < n && nb < 64; k++) {
                double y = w[4] == 'q' ? v[k] : xv[k], c0y = a[1] + y;
                if (!isfinite(y) || dfl_autoengine(DFL_F_DIST, S.model, stat, g0, gls, a[0], c0y, 0.0) != DFL_ENGINE_CLASSICAL) continue;
                const char *m = dfl_brittle(2, S.model, a[0], c0y, g0, stat, 0);
                if (!m) continue;
                int seen = 0; for (int j = 0; j < nb; j++) if (!strcmp(bm[j], m)) seen = 1;
                /* (dfl_brittle writes into one buffer of its own: the message is copied below before the next call) */
                if (!seen) { char *cp = R_alloc(strlen(m) + 1, 1); strcpy(cp, m); bm[nb++] = cp; }
            }
        }
    } else if (!strcmp(w, "lin")) {
        NEED(4); IOPTS(4); dfl_spec S1 = {NULL, (int) a[1], 0, 0.0, NAN};
        rc = dfl_engine_lincdf_x(n, xv, a[0], &S1, a[2], a[3], eng, &io, v, ev, err, sizeof err);
    } else if (!strcmp(w, "lind")) { NEED(2); IOPTS(2); rc = dfl_engine_lincdf_x(n, xv, a[0], &S, a[1], NAN, eng, &io, v, ev, err, sizeof err); }
    else if (!strcmp(w, "env")) { NEED(2); IOPTS(2); rc = dfl_engine_envelope_x(n, xv, (int) a[0], a[1], eng, &io, v, ev, err, sizeof err); }
    else if (!strcmp(w, "adf") || !strcmp(w, "adfpdf")) { NEED(1); rc = dfl_engine_adfmax(w[3] == 'p' ? 1 : 2, n, xv, &S, a[0], eng, v, err, sizeof err); }
    else if (!strcmp(w, "region")) { NEED(4); double par[2] = {a[1], a[2]}; rc = dfl_engine_region((int) a[0], par, &S, a[3], eng, v, err, sizeof err); }
    else Rf_error("unknown function %s of the router", w);
#undef NEED
#undef IOPTS
    engine_error(rc, err);
    const char *w1 = dfl_warning(), *w2 = dfl_det_warning();
    int nw = nb + (w1 != NULL) + (w2 != NULL);
    SEXP res = PROTECT(Rf_allocVector(VECSXP, 3)), ws = PROTECT(Rf_allocVector(STRSXP, nw));
    int j = 0;
    for (int i = 0; i < nb; i++) SET_STRING_ELT(ws, j++, Rf_mkChar(bm[i]));
    if (w1) SET_STRING_ELT(ws, j++, Rf_mkChar(w1));
    if (w2) SET_STRING_ELT(ws, j++, Rf_mkChar(w2));
    SET_VECTOR_ELT(res, 0, out); SET_VECTOR_ELT(res, 1, ws); SET_VECTOR_ELT(res, 2, eo);
    UNPROTECT(4);
    return res;
}

/* .Call("dfl_r_check", model, c, c0, g0, stat): the message of the classical engine's checks of the arguments
 * (dfl_check_args), "" where it accepts them; for the R side's choice of the engine of a design that is a classical
 * case */
SEXP dfl_r_check(SEXP model, SEXP c, SEXP c0, SEXP g0, SEXP stat)
{
    char err[ERRLEN] = "";
    int rc = dfl_check_args(int1(model, "model"), num1(c, "c"), num1(c0, "c0"), num1(g0, "g0"), int1(stat, "stat") != 0, err, sizeof err);
    return mkstr(rc ? (err[0] ? err : "the classical engine refuses the arguments") : "");
}

/* .Call("dfl_r_empty_caches"): frees the cached prepared designs, paths and refined lines of the general engine */
SEXP dfl_r_empty_caches(void)
{
    dfl_det_empty_caches();
    return R_NilValue;
}
