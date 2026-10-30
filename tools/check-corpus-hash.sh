#!/bin/sh
# SPDX-License-Identifier: LGPL-3.0-only
#
# Fail unless the hash a probe prints is the committed one (CAP-13).
#
# Usage: check-corpus-hash.sh PROBE HASH_FILE
#
# PROBE is tools/corpus-hash.c built for this host (or a planted variant);
# HASH_FILE is tests/corpus/geometry.hash, whose `hash` and `count` lines are
# the committed values, and whose `predicate-hash` and `predicate-count` lines are the predicate section's (a hash of its
# own, so the first pair never moves) and whose `vertex-hash` and `vertex-count` lines are the vertex section's, a third; the arc,
# arc-vertex and conic pairs follow in the same way.
# No value may be empty, so a probe that crashed, or a
# hash file that was emptied, fails rather than comparing nothing with nothing.
probe="${1:?usage: check-corpus-hash.sh PROBE HASH_FILE}"
file="${2:?usage: check-corpus-hash.sh PROBE HASH_FILE}"

got="$("$probe")" || { printf 'check-corpus-hash: the probe failed\n' >&2; exit 1; }
got_hash=$(printf '%s\n' "$got" | awk '$1=="shape-geometry.hash"{print $2}')
got_count=$(printf '%s\n' "$got" | awk '$1=="shape-geometry.count"{print $2}')
want_hash=$(awk '$1=="hash"{print $2}' "$file" 2>/dev/null)
want_count=$(awk '$1=="count"{print $2}' "$file" 2>/dev/null)
got_phash=$(printf '%s\n' "$got" | awk '$1=="shape-geometry.predicate-hash"{print $2}')
got_pcount=$(printf '%s\n' "$got" | awk '$1=="shape-geometry.predicate-count"{print $2}')
want_phash=$(awk '$1=="predicate-hash"{print $2}' "$file" 2>/dev/null)
want_pcount=$(awk '$1=="predicate-count"{print $2}' "$file" 2>/dev/null)
got_vhash=$(printf '%s\n' "$got" | awk '$1=="shape-geometry.vertex-hash"{print $2}')
got_vcount=$(printf '%s\n' "$got" | awk '$1=="shape-geometry.vertex-count"{print $2}')
want_vhash=$(awk '$1=="vertex-hash"{print $2}' "$file" 2>/dev/null)
want_vcount=$(awk '$1=="vertex-count"{print $2}' "$file" 2>/dev/null)
got_ahash=$(printf '%s\n' "$got" | awk '$1=="shape-geometry.arc-hash"{print $2}')
got_acount=$(printf '%s\n' "$got" | awk '$1=="shape-geometry.arc-count"{print $2}')
want_ahash=$(awk '$1=="arc-hash"{print $2}' "$file" 2>/dev/null)
want_acount=$(awk '$1=="arc-count"{print $2}' "$file" 2>/dev/null)
got_avhash=$(printf '%s\n' "$got" | awk '$1=="shape-geometry.arc-vertex-hash"{print $2}')
got_avcount=$(printf '%s\n' "$got" | awk '$1=="shape-geometry.arc-vertex-count"{print $2}')
want_avhash=$(awk '$1=="arc-vertex-hash"{print $2}' "$file" 2>/dev/null)
want_avcount=$(awk '$1=="arc-vertex-count"{print $2}' "$file" 2>/dev/null)
got_chash=$(printf '%s\n' "$got" | awk '$1=="shape-geometry.conic-hash"{print $2}')
got_ccount=$(printf '%s\n' "$got" | awk '$1=="shape-geometry.conic-count"{print $2}')
want_chash=$(awk '$1=="conic-hash"{print $2}' "$file" 2>/dev/null)
want_ccount=$(awk '$1=="conic-count"{print $2}' "$file" 2>/dev/null)
got_shash=$(printf '%s\n' "$got" | awk '$1=="shape-geometry.shapes-hash"{print $2}')
got_scount=$(printf '%s\n' "$got" | awk '$1=="shape-geometry.shapes-count"{print $2}')
want_shash=$(awk '$1=="shapes-hash"{print $2}' "$file" 2>/dev/null)
want_scount=$(awk '$1=="shapes-count"{print $2}' "$file" 2>/dev/null)
if [ -z "$want_hash" ] || [ -z "$want_count" ] || [ -z "$got_hash" ] || [ -z "$got_count" ] ||
   [ -z "$want_phash" ] || [ -z "$want_pcount" ] || [ -z "$got_phash" ] || [ -z "$got_pcount" ] ||
   [ -z "$want_vhash" ] || [ -z "$want_vcount" ] || [ -z "$got_vhash" ] || [ -z "$got_vcount" ] ||
   [ -z "$want_ahash" ] || [ -z "$want_acount" ] || [ -z "$got_ahash" ] || [ -z "$got_acount" ] ||
   [ -z "$want_avhash" ] || [ -z "$want_avcount" ] || [ -z "$got_avhash" ] || [ -z "$got_avcount" ] ||
   [ -z "$want_chash" ] || [ -z "$want_ccount" ] || [ -z "$got_chash" ] || [ -z "$got_ccount" ] ||
   [ -z "$want_shash" ] || [ -z "$want_scount" ] || [ -z "$got_shash" ] || [ -z "$got_scount" ]; then
  printf 'check-corpus-hash: nothing to compare (probe %s/%s, predicate %s/%s and vertex %s/%s, committed %s/%s, %s/%s and %s/%s in %s)\n' \
    "$got_hash" "$got_count" "$got_phash" "$got_pcount" "$got_vhash" "$got_vcount" "$want_hash" "$want_count" "$want_phash" "$want_pcount" \
    "$want_vhash" "$want_vcount" "$file" >&2
  exit 1
fi
if [ "$got_hash" != "$want_hash" ] || [ "$got_count" != "$want_count" ]; then
  printf 'check-corpus-hash: FAIL: %s over %s results, committed %s over %s\n' \
    "$got_hash" "$got_count" "$want_hash" "$want_count" >&2
  exit 1
fi
if [ "$got_phash" != "$want_phash" ] || [ "$got_pcount" != "$want_pcount" ]; then
  printf 'check-corpus-hash: FAIL: predicate %s over %s results, committed %s over %s\n' \
    "$got_phash" "$got_pcount" "$want_phash" "$want_pcount" >&2
  exit 1
fi
if [ "$got_vhash" != "$want_vhash" ] || [ "$got_vcount" != "$want_vcount" ]; then
  printf 'check-corpus-hash: FAIL: vertex %s over %s results, committed %s over %s\n' \
    "$got_vhash" "$got_vcount" "$want_vhash" "$want_vcount" >&2
  exit 1
fi
if [ "$got_ahash" != "$want_ahash" ] || [ "$got_acount" != "$want_acount" ]; then
  printf 'check-corpus-hash: FAIL: arc %s over %s results, committed %s over %s\n' \
    "$got_ahash" "$got_acount" "$want_ahash" "$want_acount" >&2
  exit 1
fi
if [ "$got_avhash" != "$want_avhash" ] || [ "$got_avcount" != "$want_avcount" ]; then
  printf 'check-corpus-hash: FAIL: arc-vertex %s over %s results, committed %s over %s\n' \
    "$got_avhash" "$got_avcount" "$got_chash" "$got_ccount" "$want_avhash" "$want_avcount" >&2
  exit 1
fi
if [ "$got_chash" != "$want_chash" ] || [ "$got_ccount" != "$want_ccount" ]; then
  printf 'check-corpus-hash: FAIL: conic %s over %s results, committed %s over %s\n' \
    "$got_chash" "$got_ccount" "$want_chash" "$want_ccount" >&2
  exit 1
fi
if [ "$got_shash" != "$want_shash" ] || [ "$got_scount" != "$want_scount" ]; then
  printf 'check-corpus-hash: FAIL: shapes %s over %s results, committed %s over %s\n' \
    "$got_shash" "$got_scount" "$want_shash" "$want_scount" >&2
  exit 1
fi
printf 'check-corpus-hash: ok, %s over %s results; predicate %s over %s; vertex %s over %s; arc %s over %s; arc-vertex %s over %s; conic %s over %s; shapes %s over %s\n' "$got_hash" "$got_count" "$got_phash" "$got_pcount" "$got_vhash" "$got_vcount" "$got_ahash" "$got_acount" "$got_avhash" "$got_avcount" "$got_chash" "$got_ccount" "$got_shash" "$got_scount"
