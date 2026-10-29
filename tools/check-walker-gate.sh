#!/bin/sh
#
# Prove that tools/check-verb-walkers.sh fails on what it exists to catch.
#
# The real script runs against a small tree built here: one library with a source file whose
# functions are audited, and a test that feeds one of them an arc.  The control must pass; each
# plant has one defect and the gate must exit non-zero and name it:
#
#   dropped-ARC          a site audited as "updated" stops naming the arc verb
#   unaudited caller     a new file calls gshgeom_path_element; and a new function in an audited file
#   deleted audit line   the line that audits a file is gone; and a line names a test that does not exist
#   empty                no library at all: the gate must not report success over nothing
#
# Usage: check-walker-gate.sh (run by `make check-verb-walkers`)

set -u

HERE="$(cd "$(dirname "$0")" && pwd)"
GATE="$HERE/check-verb-walkers.sh"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT INT TERM
failures=0
checks=0

fixture() {
  # fixture <dir>: a library named like this one, so that the gate scans it
  root="$1/$(basename "$(dirname "$HERE")")"
  mkdir -p "$root/src" "$root/tests/unit" "$root/tools"
  printf 'LIBS := shape-geometry\n' > "$root/Makefile"
  cat > "$root/src/walk.c" <<'EOC'
#include <stddef.h>
static int handles(const void * path) {
  int verb = GSHGEOM_VERB_ARC;
  return gshgeom_path_verb(path, 0, &verb);
}
static int refuses(const void * path) {
  GSHGEOM_Element el;
  return gshgeom_path_element(path, 0, &el);
}
EOC
  cat > "$root/tests/unit/test_walk.cpp" <<'EOC'
TEST(Walk, FeedsTheSitesAnArc) { int v = GSHGEOM_VERB_ARC; }
EOC
  cat > "$1/audit" <<EOA
$(basename "$root")|src/walk.c|handles|updated|handles an arc|$(basename "$root"):tests/unit/test_walk.cpp:Walk.FeedsTheSitesAnArc
$(basename "$root")|src/walk.c|refuses|confirmed|passes ERR_UNSUPPORTED on|$(basename "$root"):tests/unit/test_walk.cpp:Walk.FeedsTheSitesAnArc
$(basename "$root")|tests/unit/test_walk.cpp|*|support|names no verb of its own|-
EOA
}

# expect <label> <needle> <expected exit: pass|fail> <libs> <audit>
run_gate() {
  VW_LIBS="$1" VW_AUDIT="$2" "$GATE" 2>&1
}

expect_pass() {
  label="$1"; shift
  checks=$((checks + 1))
  if out="$(run_gate "$@")"; then
    printf '  ok   %s passes\n' "$label"
  else
    printf 'check-walker-gate: FAIL: %s should pass but failed:\n%s\n' "$label" "$out" >&2
    failures=$((failures + 1))
  fi
}

expect_fail() {
  label="$1"; needle="$2"; shift 2
  checks=$((checks + 1))
  out="$(run_gate "$@")"; rc=$?
  if [ "$rc" -eq 0 ]; then
    printf 'check-walker-gate: FAIL: %s should fail but exited 0:\n%s\n' "$label" "$out" >&2
    failures=$((failures + 1))
  elif ! printf '%s' "$out" | grep -qF -- "$needle"; then
    printf 'check-walker-gate: FAIL: %s failed without naming "%s":\n%s\n' "$label" "$needle" "$out" >&2
    failures=$((failures + 1))
  else
    printf '  ok   %s fails, naming %s\n' "$label" "$needle"
  fi
}

mkdir -p "$work/control"
fixture "$work/control"
lib="$(basename "$(dirname "$HERE")")"
expect_pass 'walkers/control' "$work/control" "$work/control/audit"

# dropped-ARC: the site stays audited as "updated" but no longer names the arc verb or reader
mkdir -p "$work/drop"; fixture "$work/drop"
sed -i 's/GSHGEOM_VERB_ARC/GSHGEOM_VERB_LINE/; s/gshgeom_path_verb(path, 0, &verb)/0/' "$work/drop/$lib/src/walk.c"
expect_fail 'walkers/planted (a site drops ARC)' 'audited as updated for arcs' "$work/drop" "$work/drop/audit"

# an unaudited caller in a new file
mkdir -p "$work/new"; fixture "$work/new"
printf 'int late(const void * p) { GSHGEOM_Element e; return gshgeom_path_element(p, 0, &e); }\n' > "$work/new/$lib/src/late.c"
expect_fail 'walkers/planted (a caller outside the audit)' "src/late.c: the function late" "$work/new" "$work/new/audit"

# an unaudited function in an audited file
mkdir -p "$work/fn"; fixture "$work/fn"
printf 'static int extra(const void * p) { GSHGEOM_Element e; return gshgeom_path_element(p, 1, &e); }\n' >> "$work/fn/$lib/src/walk.c"
expect_fail 'walkers/planted (a new function in an audited file)' "the function extra" "$work/fn" "$work/fn/audit"

# a deleted audit line
mkdir -p "$work/del"; fixture "$work/del"
grep -v '|refuses|' "$work/del/audit" > "$work/del/audit.new"
expect_fail 'walkers/planted (an audit line deleted)' "the function refuses" "$work/del" "$work/del/audit.new"

# an audit line whose test does not exist
mkdir -p "$work/test"; fixture "$work/test"
sed -i 's/Walk.FeedsTheSitesAnArc/Walk.NoSuchTest/' "$work/test/audit"
expect_fail 'walkers/planted (a test that does not exist)' 'Walk.NoSuchTest is not in' "$work/test" "$work/test/audit"

# a file with no audit line at all (tests, tools and examples need one too)
mkdir -p "$work/file"; fixture "$work/file"
grep -v 'test_walk.cpp|\*|' "$work/file/audit" > "$work/file/audit.new"
expect_fail 'walkers/planted (a test file outside the audit)' 'test_walk.cpp: names a path verb' "$work/file" "$work/file/audit.new"

# nothing at all
mkdir -p "$work/empty"
expect_fail 'walkers/empty (no library)' 'found no walker site' "$work/empty" "$work/control/audit"

if [ "$failures" -ne 0 ]; then
  exit 1
fi
printf 'check-walker-gate: all %d checks passed; the gate fails on each planted defect\n' "$checks"
