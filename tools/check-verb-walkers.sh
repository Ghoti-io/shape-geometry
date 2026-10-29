#!/bin/sh
#
# The walker gate: every site that reads a path by verb is in tools/verb-walkers.audit.
#
# A path can hold an arc (arc.h).  A function that reads a path element by element, or switches
# on the verb, either handles the arc verb or refuses it loudly, and the audit says which, with
# a reason and a test that feeds the site an arc.  This gate finds the sites again, by the same
# search each time, so that a caller added later cannot land unaudited:
#
#   * it scans every library next to this one that links shape-geometry (and shape-geometry
#     itself): source, tests, tools and examples, but not the public headers of shape-geometry,
#     which only declare;
#   * a site is a line that names gshgeom_path_element, gshgeom_path_element_count,
#     gshgeom_path_vertices, gshgeom_path_vertex_count, a GSHGEOM_VERB_ constant or
#     GSHGEOM_Element;
#   * every file with a site needs a line in the audit; for a file under src/ every function
#     that holds a site needs to be named in one (the symbol column; "<file>" for a site outside
#     any function);
#   * every audit line needs a known disposition and a reason; a file and symbols that exist; a
#     test that exists (library:path:Suite.Name, or "-" for a support file, or "A2" for a site
#     that story edits); a line that says "updated" needs the file to name the arc verb, the
#     arc reader or the lowering, so that a site which drops its arc handling fails.
#
# Dispositions:
#   updated      handles an arc (the reason cites the rule)
#   confirmed    receives ::GSHGEOM_ERR_UNSUPPORTED from gshgeom_path_element (or passes it on) and
#                fails or declines loudly, never drawing the chord
#   unreachable  an arc cannot get there; the reason says why
#   support      a test, tool or example that builds or reads arc-free paths
#   A2           a shape-scene site the named-primitives story edits (pending)
#   owner        a site in a library this one does not own; it fails loudly meanwhile, and the
#                owner updates it
#
# Usage: check-verb-walkers.sh [--list]
#   --list       print the sites found, as library|file|function, and exit
# Environment (used by check-gates.sh to plant defects):
#   VW_LIBS      the directory that holds the libraries (default: the parent of this one)
#   VW_AUDIT     the audit file (default: tools/verb-walkers.audit)

set -u

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(dirname "$HERE")"
LIBS="${VW_LIBS:-$(dirname "$ROOT")}"
AUDIT="${VW_AUDIT:-$HERE/verb-walkers.audit}"
SELF="$(basename "$ROOT")"

PATTERN='gshgeom_path_element([^A-Za-z0-9_]|_count|$)|gshgeom_path_vertices|gshgeom_path_vertex_count|GSHGEOM_VERB_[A-Z_]+|GSHGEOM_Element([^A-Za-z0-9_]|$)'

work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT INT TERM

# The libraries that link shape-geometry, and shape-geometry.
libs=""
for d in "$LIBS"/*/; do
  name="$(basename "$d")"
  if [ "$name" = "$SELF" ] || { [ -f "$d/Makefile" ] && grep -q 'shape-geometry' "$d/Makefile"; }; then
    libs="$libs $name"
  fi
done

# sites.tsv: library|file|function, one per (file, function) holding a site.
: > "$work/sites"
for lib in $libs; do
  dir="$LIBS/$lib"
  for sub in src tests tools examples; do
    [ -d "$dir/$sub" ] || continue
    find "$dir/$sub" -type f \( -name '*.c' -o -name '*.h' -o -name '*.cpp' -o -name '*.py' -o -name '*.json' -o -name '*.sh' \) \
        ! -path '*/build/*' ! -path '*/corpus/*' ! -name 'verb-walkers.audit' ! -name 'check-verb-walkers.sh' ! -name 'check-walker-gate.sh' \
        ! -name 'arc_vectors.h' ! -path '*/include/ghoti.io/shape-geometry/*' | sort > "$work/files"
    while IFS= read -r f; do
      rel="${f#"$dir"/}"
      awk -v lib="$lib" -v rel="$rel" -v pat="$PATTERN" '
        BEGIN { cur = "<file>" }
        /^[A-Za-z_]/ && /\(/ && !/^(typedef|struct|enum|union|#)/ { l = $0; sub(/\(.*/, "", l); n = split(l, a, /[ \t*]+/); cur = a[n] }
        $0 ~ pat { seen[cur] = 1 }
        /^}/ { cur = "<file>" }
        END { for (s in seen) print lib "|" rel "|" s }
      ' "$f"
    done < "$work/files" >> "$work/sites"
  done
done
sort -u "$work/sites" -o "$work/sites"

if [ "${1:-}" = "--list" ]; then
  cat "$work/sites"
  exit 0
fi

fail=0
say() {
  printf 'check-verb-walkers: %s\n' "$*" >&2
  fail=1
}

if [ ! -s "$work/sites" ]; then
  say "found no walker site at all in $LIBS (libraries:$libs); the search itself is broken"
fi
if [ ! -f "$AUDIT" ]; then
  say "no audit file at $AUDIT"
  exit 1
fi

# Audit lines: library|file|symbols|disposition|reason|test
: > "$work/audited"
n=0
skipped=0
while IFS= read -r line; do
  n=$((n + 1))
  case "$line" in ''|'#'*) continue ;; esac
  oldifs="$IFS"; IFS='|'
  set -f
  set -- $line
  set +f
  IFS="$oldifs"
  if [ "$#" -ne 6 ]; then
    say "$AUDIT:$n: a line has six fields (library|file|symbols|disposition|reason|test): $line"
    continue
  fi
  lib="$1"; file="$2"; syms="$3"; disp="$4"; reason="$5"; test="$6"
  path="$LIBS/$lib/$file"
  if [ ! -d "$LIBS/$lib" ]; then
    skipped=$((skipped + 1)) # a library that is not checked out beside this one is not audited here
    continue
  fi
  case "$disp" in
    updated|confirmed|unreachable|support|A2|owner) ;;
    *) say "$lib/$file: unknown disposition '$disp'" ;;
  esac
  [ -n "$reason" ] || say "$lib/$file: no reason"
  if [ ! -f "$path" ]; then
    say "$lib/$file: the audited file does not exist"
    continue
  fi
  # symbols
  oldifs="$IFS"; IFS=','
  set -f
  for s in $syms; do
    IFS="$oldifs"
    printf '%s|%s|%s\n' "$lib" "$file" "$s" >> "$work/audited"
    case "$s" in
      '<file>'|'*') ;;
      *) grep -q "$s" "$path" || say "$lib/$file: the symbol $s is not in the file" ;;
    esac
    IFS=','
  done
  IFS="$oldifs"
  set +f
  # the test
  case "$disp" in
    support) [ "$test" = "-" ] || say "$lib/$file: a support file names no test (-)" ;;
    owner) [ "$test" = "-" ] || [ -f "$LIBS/${test%%:*}/$(echo "$test" | cut -d: -f2)" ] || say "$lib/$file: the test file '$test' does not exist" ;;
    A2) [ "$test" = "A2" ] || say "$lib/$file: an A2 site's test is the word A2" ;;
    *)
      tlib="${test%%:*}"; rest="${test#*:}"; tfile="${rest%%:*}"; tname="${rest#*:}"
      suite="${tname%%.*}"; case_="${tname#*.}"
      if [ "$tlib" = "$test" ] || [ ! -f "$LIBS/$tlib/$tfile" ]; then
        say "$lib/$file: the test file '$test' does not exist"
      elif ! grep -qE "^TEST(_F)?\($suite, $case_\)" "$LIBS/$tlib/$tfile"; then
        say "$lib/$file: the test $tname is not in $tlib/$tfile"
      fi
      ;;
  esac
  if [ "$disp" = "updated" ] && ! grep -qE 'GSHGEOM_VERB_ARC|gshgeom_path_lower_arcs|gshgeom_path_arc\(|gshgeom_path_verb\(' "$path"; then
    say "$lib/$file: audited as updated for arcs but the file never names GSHGEOM_VERB_ARC, gshgeom_path_arc, gshgeom_path_verb or gshgeom_path_lower_arcs"
  fi
  if [ "$disp" = "support" ]; then
    case "$file" in tests/*|tools/*|examples/*) ;; *) say "$lib/$file: 'support' is for tests, tools and examples" ;; esac
  fi
done < "$AUDIT"
sort -u "$work/audited" -o "$work/audited"

# Every site is audited; a source file's every function holding a site is named.
while IFS='|' read -r lib file func; do
  case "$file" in
    src/*)
      if ! grep -qxF "$lib|$file|$func" "$work/audited"; then
        say "$lib/$file: the function $func reads the path by verb and is not in the audit"
      fi
      ;;
    *)
      if ! grep -q "^$lib|$file|" "$work/audited"; then
        say "$lib/$file: names a path verb or element and is not in the audit"
      fi
      ;;
  esac
done < "$work/sites"

# A line that audits a file with no site is stale, unless it says the site is unreachable by design.
cut -d'|' -f1,2 "$work/audited" | sort -u | while IFS='|' read -r lib file; do
  if ! grep -q "^$lib|$file|" "$work/sites" && ! grep -qE "^$lib\|$file\|[^|]*\|(unreachable|owner)\|" "$AUDIT"; then
    printf 'check-verb-walkers: %s/%s: audited, but the file names no walker site\n' "$lib" "$file" >&2
    exit 1
  fi
done || fail=1

if [ "$fail" -ne 0 ]; then
  exit 1
fi
printf 'check-verb-walkers: %d sites in %d libraries are all audited (%d lines for libraries not checked out)\n' \
  "$(wc -l < "$work/sites")" "$(echo $libs | wc -w)" "$skipped"
