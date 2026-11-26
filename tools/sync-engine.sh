#!/bin/sh
# tools/sync-engine.sh: copy the C engine from a checkout of dflimit-matlab, the repository of the MATLAB package that
# holds the C engine, into src/ of this package.
#
#   sh tools/sync-engine.sh /path/to/dflimit-matlab
#
# The script only reads the dflimit-matlab checkout.  It
#   1. copies the engine sources of c/ (and c/faddeeva/, flattened into src/) into src/,
#   2. converts them to LF line endings (a Windows checkout with core.autocrlf=true has CRLF, which R CMD check
#      reports),
#   3. applies the patches in tools/patches/*.patch, if there are any, in name order (a patch is a change that R
#      requires and that should go upstream into c/ of dflimit-matlab; once it is upstream, delete the patch file; there
#      are none at present),
#   4. records the upstream commit in inst/UPSTREAM.
# The engine has three parts: the classical engine (dfl_api.c to tables.c), the general engine for any deterministic
# terms (dfl_det_*.c with the internal header dfl_det_int.h; dfl_det_tables.c holds its generated tables) and the
# devices of Bao and Kan (dfl_kan.c with dfl_kan_int.h and the templates dfl_kan_jet.h, dfl_kan_logcf.h, which
# dfl_kan.c includes; dfl_det_kan.c for the general engine's transform).
# Left out: the test programs (dfl_test_*.c), dfl_probe.c, gen_tables.jl, gen_det_tables.jl, the MEX gateway (c/mex/)
# and c/test/.
set -e

if [ $# -ne 1 ]; then echo "usage: sh tools/sync-engine.sh /path/to/dflimit-matlab" >&2; exit 2; fi
UP=$(cd "$1" && pwd)
PKG=$(cd "$(dirname "$0")/.." && pwd)
[ -f "$UP/c/dflimit.h" ] || { echo "$UP does not look like a dflimit-matlab checkout (no c/dflimit.h)" >&2; exit 1; }
[ -f "$PKG/DESCRIPTION" ] || { echo "$PKG is not the package root" >&2; exit 1; }

ENGINE="dflimit.h dfl_lock.h dfl_vmath.h dfl_api.c dfl_double.c dfl_gen.c dfl_gilpelaez.c dfl_mgf.c dfl_moments.c
        dfl_single.c dfl_batch.c dfl_special.c dfl_threads.c tables.c
        dfl_det_int.h dfl_det_api.c dfl_det_design.c dfl_det_inner.c dfl_det_la.c dfl_det_mgf.c dfl_det_paper2.c
        dfl_det_path.c dfl_det_route.c dfl_det_line.c dfl_det_tables.c dfl_det_ext.c
        dfl_kan_int.h dfl_kan.c dfl_kan_jet.h dfl_kan_logcf.h dfl_det_kan.c"
FADDEEVA="Faddeeva.cc Faddeeva.h faddeeva_shim.cpp"

lf() { tr -d '\r' < "$1" > "$2"; }
for f in $ENGINE; do lf "$UP/c/$f" "$PKG/src/$f"; done
for f in $FADDEEVA; do lf "$UP/c/faddeeva/$f" "$PKG/src/$f"; done

# the patches (unified diffs relative to the package root, paths src/...)
PATCHES=""
for p in "$PKG"/tools/patches/*.patch; do
    [ -f "$p" ] || continue
    patch -p1 -s --forward --no-backup-if-mismatch -d "$PKG" < "$p"
    PATCHES="$PATCHES $(basename "$p")"
done

# the upstream commit (read-only git commands: --no-optional-locks keeps git from refreshing the index of the checkout)
G="git --no-optional-locks -C $UP"
COMMIT=$($G rev-parse HEAD 2>/dev/null || echo unknown)
CDATE=$($G log -1 --format=%cI 2>/dev/null || echo unknown)
if $G diff --quiet HEAD -- c 2>/dev/null; then DIRTY=no; else DIRTY=yes; fi
mkdir -p "$PKG/inst"
{
    echo "# The C engine in src/ was copied from dflimit-matlab (its directory c/) by tools/sync-engine.sh"
    echo "upstream_commit: $COMMIT"
    echo "upstream_commit_date: $CDATE"
    echo "upstream_c_modified: $DIRTY"
    echo "patches:$PATCHES"
} > "$PKG/inst/UPSTREAM"
echo "synced the engine from $UP at $COMMIT (c/ modified: $DIRTY); patches:${PATCHES:- none}"
