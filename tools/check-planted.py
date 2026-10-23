#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-only
#
# Copyright (C) 2026 Corey Pennycuff
#
# This file is part of Ghoti.io Shape Geometry.
#
# Ghoti.io Shape Geometry is free software: you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License version 3
# as published by the Free Software Foundation.
"""Show that the tests fail on planted defects (CAP-2, CAP-13).

A test that has never been seen to fail may measure nothing.  For each defect
below this builds the library's sources with one literal replacement applied
(the text must match exactly once and the replacement must change it), links
the tests that are meant to catch it against the result, and requires one of
them to fail.  The unplanted build is the control: every test must pass on it,
which is what makes a failure mean the defect and not the build.

A plant is a replacement in a source file (``src``), or in a public header
(``hdr``: the whole library is rebuilt with the changed header ahead of the
real one on the include path; the tests keep the real one, which is how a
changed constant disagrees with its pinned value), or a compile flag
(``flags``) caught by the cross-target probe's hash.  Families:

  hash     the hash spelling: the basis and prime, xor before multiply, byte
           order and width of each field, the string prefix, the bits of a
           double; and the path's fields, their order, the verb codes
  refuse   the validity rules: non-finite, the magnitude and weight bounds,
           the order in which a NaN and a limit are reported, the ignored z
  eval     de Casteljau, the homogeneous coordinates, the ends, the range of t
  bounds   the candidate parameters, the Bernstein form of N'D - ND', the root
           isolation and the bisection, the min/max
  path     the builder's rules, the limit, the start of a subpath, readback,
           the box of a path, allocation failure and release
  probe    a bit flipped in the cross-target probe, and a fused multiply-add

Story 6 adds (ids start with x-, o-, m-, mh-, md- and route-):

  x-       the exact sum of products: the cut of a double into mantissa and
           exponent, the product of mantissas, the shifted two's complement
           accumulator and its carries and borrows, the sign, the determinant
           expansion
  o-       the predicates: the filter's bounds and guards, the expanded
           determinants term by term, the refusals, the normals, the area
  m-       the maps: the constants and enumerations of map.h, the angle
           reduction and the exact quarters, the sense of every rotation, the
           classifier branch by branch, the composer's operands, order of
           operations, tag union and collapse, application, curves
  route-   the route predicate: each class the tag lets through or not, each
           value test, the grid it returns
  mh-, md- the map's hash fields and dump text

Story 7 adds (ids start with r-, dp-, b-, h-, fm-, rat-, fl- and fh-):

  r-       the rounding rule: the bounds of the range and whether each edge is
           inclusive, the half-way comparison, the double expression `floor(t +
           0.5)` that is wrong at 0.5 - 2^-54, round-half-away and truncation,
           the factor 64, a write on failure, the refusal of NaN
  dp-      the point function: the dimension, which coordinate is rounded
  b-       the device path builder: the verb codes, the limit and the all-or-
           nothing move, what a close stores, the verb after a close, readback,
           release and allocation failure
  h-       the device path's hash fields and dump text
  fm-      a model path under a map: the dimensions, the order of elements, the
           mapped points of each verb, the leak on failure, the range rule
  rat-     the rational subdivision: the tolerance, the flatness test, the
           clamp to the chord's segment, the de Casteljau halves, the weights
  fl-      the flattener: the midpoint's floor, the deviation sums term by term,
           the threshold, the order and the propagation of errors, the limit,
           the subpath bookkeeping
  fh-      the flattened path's readback, hash fields and dump text

A plant of kind ``srch`` changes an internal header under src/ and rebuilds
what includes it (device.c and flatten.c).

Story 8 adds (ids start with al-, alh- and ald-):

  al-      along-path parametrization: the Gauss-Legendre constants (one unit in
           the last place each, which only the pinned tables and the probe can
           see), the homogeneous form, de Casteljau and the differences, the
           derivative of E = N'D - ND' by Leibniz's rule, the speed, the rule's
           scaling, the polygon length that scales the tolerance, the table's
           rows and the exact row limit, the acceptance test of a panel and the
           order and error propagation of its halves, the tolerance argument and
           its bounds, the single-curve length, building (subpaths, the first
           point, release on failure), the accessors, the query (range, clamping,
           the boundary rule, the end of a row, the frame), the inverse of the
           arc length (Newton, the bracket, the stop) and the tangent (the
           threshold, the orders of derivative, the sign at the end)
  alh-     the measured path's hash fields, their order and width
  ald-     its dump text

Story 14a adds (ids start with fx- and sc-):

  fx-      the fixed-point builder: the bound on the fraction bits (each verb's
           check, the constant), the exact division, the coordinate order, the
           z, the order of a curve's points
  sc-      the integer spelling of the rule: the flip and its argument check, the
           half, the floor of a negative number, the int32 range of the result

Strokes add (ids start with sk-): the parameters, their order of refusal and
their bounds (SG1 to SG3), the determinant (SG15), the statuses and the parts
(SG12, SG14, SG16), the effective tangent (SG11), every join and its reversal
limit and SG9's exact decision (SG8, SG9), the 2018 CR's arcs join (its
fallbacks, adjustment, meeting point, clip and limit), the body's bands and fans
(SG6), circular pieces, the emission's inflation, orientation and rounding, the
outline's hash and dump (SG18), and the wide exact sum.  A test list may name
``probe``: the plant is then caught by the cross-target probe's hash as well.

Run by ``make check-planted`` from the library directory, with the build's
flags in the PLANT_* environment variables.

What is NOT planted, and why.  A defect that no input can reach, or that
changes no result, cannot be seen by any test:

  - the scale of the derivative in ``roots01`` (``m *`` in ``d[k]``): a positive
    multiple of a polynomial has the same roots
  - a root found exactly at a critical point (``fa == 0`` at an interior
    breakpoint): that is a touching zero of x'(t), which is not an extremum of x
    (the end points are candidates anyway)
  - ``fb != 0.0`` before the sign test in ``roots01``: with ``fb`` zero the
    bisection converges to ``b``, which the next interval finds again
  - the sign of the whole Bernstein polynomial (``x_i - x_{j+1}``): a negated
    polynomial has the same roots
  - the explicit ``!w`` check for a rational kind in ``make_curve``: without it
    a NULL weight array leaves every weight zero, which is refused as well
    (the explicit check is kept so that the refusal does not rest on that)
  - ``choose`` row 3: never read (``n - 1`` is at most 2 and ``m`` is 2 or 4)
  - a scale on the whole of the Bernstein polynomial (the factor ``n`` that the
    comment on ``axis_extrema`` leaves out, or any other positive constant):
    a positive multiple has the same roots
  - the result of ``bisect`` after the loop (``a`` for the midpoint): the loop
    ends when the interval is one unit in the last place wide
  - a root exactly at t = 1 (``fa == 0`` at the end of ``roots01``): the end
    point is a candidate already
  - ``s[6]`` in ``accumulate`` (the word the last mantissa word spills into): a
    product of at most three 53-bit mantissas is 159 bits and fills five words,
    so the sixth is zero and so is its spill
  - the coordinate limit of the predicates' filter (``FILTER_COORD_LIMIT`` and
    ``small_enough``): the filter is valid far beyond it, and where a product
    overflows the result is an infinity of the right sign or a NaN that fails
    every comparison and goes to the exact path; the limit is a margin
  - the shortcut ``rest == 0.0`` in ``angle_cs``: ``gmath_cos(0)`` is 1 and
    ``gmath_sin(0)`` is 0 exactly, so the shortcut changes nothing; it is kept
    so that "a quarter turn is exact" does not rest on ``math``'s special values
  - the explicit NULL checks on the two maps in ``gshgeom_map_compose``: the
    validity check that follows refuses a NULL map as well, with the same
    result
  - the check of the first write in ``gshgeom_map_dump`` as a check: every
    later write is checked too, and a stream that failed once fails again, so
    the header's check only changes which line reports the error
  - ``RATIONAL_MAX_DEPTH`` raised, and the depth test of ``rational_lines`` made
    ``< 0``, or ``depth - 1`` left as ``depth``: the flat test is met long before
    the cap, so the cap is reached by no input within the int32 range; lowering
    it is planted
  - the flattener's depth cap (story 7a), as the cap and not as a number: the
    proof (flatten.c, (g)) is that no int32 input reaches it, so removing the
    ``depth <= 0`` test, testing ``depth < 0``, or passing ``depth + 1`` to the
    right half changes no output, and the compiler may not be asked to show
    otherwise.  What is planted is a cap that is too low: 3, 4 levels early (both
    degrees), 8, and 16, which the tolerance of 2 (threshold 1,024 widened)
    reaches on a curve spanning the int32 range, as the 20 is for; the cap
    clamp of ``flattener_begin`` at its upper end (a ``levels`` past 20 is
    brought to 20) is the same equivalence, as is ``GSHGEOM_FLATTEN_MAX_DEPTH``
    above 20
  - the shift spelling of the widening (``<< FLAT_EXTRA`` for ``* FLAT_SCALE``):
    the same number
  - a threshold of ``3 * t`` widened units for ``4 * (t - 1)`` is planted: it is
    the same number at the default tolerance (3,072) and differs at every other,
    which the tolerance tests pin (when the tolerance was fixed, this was an
    equivalent mutant)
  - the ``_Static_assert`` that the minimum tolerance is at least 2 and that the
    default is not below the minimum: a plant that does not compile cannot be run
    here (a ``GSHGEOM_FLATTEN_MIN_TOLERANCE`` of 1 is refused by it, which is the
    gate working); the minimum above 2 and the runtime refusal are planted
  - the guard that the first flattened element is a move: the builder
    guarantees it
  - the guard `path->el_n > path->max_elements` in `append` (the element count
    never exceeds the limit, so the subtraction after it cannot wrap), and the
    `grown > SIZE_MAX / size` test of `gshgeom_grow` (an array of that size is
    refused by the allocator long before; the doubling test above it is the one
    a caller can reach)
  - the check of one write in a dump other than the first (the newline after a
    device path element, the numbers of a device path element, a flattened path's
    subpath line or point line): a
    glibc stream that fails a write keeps the data it could not write, so the
    final `fflush` fails too, and the dump returns ``ERR_IO`` from there; only a
    libc whose streams discard failed data would show the difference. The first
    write's check and the flush are planted, and the first write is the one a
    read-only stream fails
  - the `!map` in `gshgeom_device_path_from_model`'s first test (the map check
    that follows refuses NULL with the same result), the result of
    `gshgeom_map_apply_point` in `map_scaled` (a model point and a map that were
    both validated cannot give a non-finite point: 2^100 times 2^100 is far
    below DBL_MAX), the result of `gshgeom_device_path_element` in the
    flattener (the index is in range), and the `!out` test of
    `gshgeom_path_flatten` (the flattener it calls refuses NULL the same way,
    after wasting the conversion)
  - ``fp->sub_n--`` after a refused first point of a subpath: the whole result
    is freed on any failure
  - the explicit ``len2`` guard's comparison (``> 0`` against ``!= 0``) and the
    dropped ``map_check`` in ``gshgeom_device_point``: the callee repeats them
  - the ``p[0]`` against ``p[1]`` of a move in ``convert_element``: for a move
    the two are the same point
  - ``grow.h``'s starting capacity, doubling and the ``*cap`` update: they
    change only how often the allocator is called (the pointer update and the
    refusal are planted)
  - ``add_pt`` without the stored move: the array has slack up to its capacity
    of 16 and a plain run cannot see a store past the need (``make test-asan``
    is the instrument for it, and the sweep in ``DevicePathAllocation`` grows
    the arrays through it)
  - story 8: the order test ``k > m`` in ``deriv`` (``k > m + 1``): a derivative
    beyond the degree is asked for only at a point where every lower derivative
    of E vanishes too, which a curve with extent does not have; ``z`` of E in 2D
    (``mag`` ignores it); the polynomial branch of ``speed`` (the denominator of a
    polynomial curve is a de Casteljau blend of ones, which is exactly one for
    these parameters); the result of ``measure`` in the single-curve length and
    of ``gshgeom_path_element`` and ``add_segment`` in ``create`` (a failure of
    the first is reached only by more than 4 Mi rows, the next two stop the loop
    on their own condition); the initial ``end`` of a subpath (overwritten by the
    segment that follows or kept by none); a line solved by Newton instead of the
    closed form (the first guess is exact, so the answer is the same bits); the
    result of ``gshgeom_curve_eval`` in a query and the zeroing of the frame
    (its inputs are validated); the direction of Newton's step and its first
    guess (the stop test and the bracket fix the answer, only the iteration count
    differs); the threshold of "vanishes" at zero and the terms of ``e_scale``
    (the cancellations in the tests are exact or far from the threshold); the
    binomial weights of ``binom`` and their use in ``e_derivative`` (the rational
    fallback reaches E'' only at a triple coincident point, where every term of
    it points along P3 - P0, so the weights scale the direction and never turn
    it)
  - the 200-step cap of ``bisect`` as a cap: it is reached only for roots
    closer than 2^-200 to an end, where 2^-200 is not a double of interest
"""
import concurrent.futures
import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def env(name, default=None):
    value = os.environ.get(name, default)
    if value is None:
        sys.exit("check-planted: %s is not set (run it with `make check-planted`)" % name)
    return value


CC = shlex.split(env("PLANT_CC"))
CFLAGS = shlex.split(env("PLANT_CFLAGS"))
PROBE_CFLAGS = shlex.split(env("PLANT_PROBE_CFLAGS"))
PROBE_LIBS = env("PLANT_PROBE_LIBS")
CXX = env("PLANT_CXX")
CXX_INCLUDE = env("PLANT_CXX_INCLUDE")
LDFLAGS = env("PLANT_LDFLAGS")
TESTS = env("PLANT_TESTS").split()
OBJ_DIR = Path(env("PLANT_OBJ_DIR"))
LD_PATH = env("PLANT_LD_PATH")
WORK = Path(env("PLANT_WORK"))
JOBS = int(env("PLANT_JOBS", "4"))
# PLANT_ONLY=<regex> runs only the plants whose id matches (and not the probe's
# own plants): for re-running the few that a change to a test touched.  The
# control always runs.  Unset, every plant runs, which is what `make test` wants.
ONLY = os.environ.get("PLANT_ONLY")
# PLANT_PROBE_ONLY=1: the control and the probe's own plants, no source plant
# (for re-running the probe's gate after a change to the probe alone).
PROBE_ONLY = bool(os.environ.get("PLANT_PROBE_ONLY"))


def test_object(test):
    """testCurve -> build/.../tests/test_curve.o"""
    snake = re.sub(r"(?<!^)(?=[A-Z])", "_", test[len("test"):]).lower()
    return OBJ_DIR / "tests" / ("test_%s.o" % snake)


def sources():
    return sorted(str(p.relative_to(ROOT)) for p in (ROOT / "src").rglob("*.c"))


# --- the plants -------------------------------------------------------------
# (id, kind, file, old, new, tests that must catch it)
# kind: "src", "hdr", or "flags" (file is then a flag string, old/new unused).

H = "src/hash/hash.c"
C = "src/curve/curve.c"
P = "src/path/path.c"
K = "src/core/core.c"
HH = "include/ghoti.io/shape-geometry/hash.h"
CH = "include/ghoti.io/shape-geometry/curve.h"
PH = "include/ghoti.io/shape-geometry/path.h"
KH = "include/ghoti.io/shape-geometry/core.h"
FILS = "src/fill/fill.c"
FILH = "include/ghoti.io/shape-geometry/fill.h"
FILT = ["testFill"]

HT = ["testHash"]
PT = ["testPath"]
CT = ["testCurve"]

PLANTS = [
    # --- hash spelling ------------------------------------------------------
    ("basis", "hdr", HH, "UINT64_C(0xcbf29ce484222325)", "UINT64_C(0xcbf29ce484222326)", HT),
    ("prime", "hdr", HH, "UINT64_C(0x00000100000001b3)", "UINT64_C(0x00000100000001b5)", HT),
    ("fnv1-not-1a", "src", H, "hash->state = (hash->state ^ byte) * GSHGEOM_FNV_PRIME;",
     "hash->state = (hash->state * GSHGEOM_FNV_PRIME) ^ byte;", HT),
    ("big-endian-bytes", "src", H, "mix_byte(hash, (uint8_t)((value >> (8 * i)) & 0xffu));",
     "mix_byte(hash, (uint8_t)((value >> (8 * (bytes - 1 - i))) & 0xffu));", HT),
    ("u16-wide", "src", H, "void gshgeom_hash_u16(GSHGEOM_Hash * hash, uint16_t value) {\n  if (hash) {\n    mix_le(hash, value, 2);",
     "void gshgeom_hash_u16(GSHGEOM_Hash * hash, uint16_t value) {\n  if (hash) {\n    mix_le(hash, value, 4);", HT),
    ("u32-narrow", "src", H, "void gshgeom_hash_u32(GSHGEOM_Hash * hash, uint32_t value) {\n  if (hash) {\n    mix_le(hash, value, 4);",
     "void gshgeom_hash_u32(GSHGEOM_Hash * hash, uint32_t value) {\n  if (hash) {\n    mix_le(hash, value, 3);", HT),
    ("u64-narrow", "src", H, "void gshgeom_hash_u64(GSHGEOM_Hash * hash, uint64_t value) {\n  if (hash) {\n    mix_le(hash, value, 8);",
     "void gshgeom_hash_u64(GSHGEOM_Hash * hash, uint64_t value) {\n  if (hash) {\n    mix_le(hash, value, 7);", HT),
    ("i32-clamped", "src", H, "gshgeom_hash_u32(hash, (uint32_t)value);",
     "gshgeom_hash_u32(hash, (uint32_t)(value < 0 ? 0 : value));", HT),
    ("i64-clamped", "src", H, "gshgeom_hash_u64(hash, (uint64_t)value);",
     "gshgeom_hash_u64(hash, (uint64_t)(value < 0 ? 0 : value));", HT),
    ("string-no-length", "src", H, "gshgeom_hash_u32(hash, (uint32_t)(length & 0xffffffffu));", "(void)length;", HT),
    ("f64-numeric-not-bits", "src", H, "gshgeom_hash_u64(hash, gshgeom_f64_bits(value));",
     "gshgeom_hash_u64(hash, (uint64_t)(int64_t)value);", HT),
    ("bytes-drops-last", "src", H, "for (size_t i = 0; i < size; i++) {\n    mix_byte(hash, p[i]);",
     "for (size_t i = 0; i + 1 < size; i++) {\n    mix_byte(hash, p[i]);", HT),
    ("value-of-null", "src", H, "return hash ? hash->state : GSHGEOM_FNV_OFFSET_BASIS;",
     "return hash ? hash->state : 0;", HT),
    ("init-zero", "src", H, "hash->state = GSHGEOM_FNV_OFFSET_BASIS;\n  }\n}\n\nvoid gshgeom_hash_u8",
     "hash->state = 0;\n  }\n}\n\nvoid gshgeom_hash_u8", HT),
    ("f64-bits-halves", "src", H, "memcpy(&bits, &value, sizeof bits);",
     "memcpy(&bits, &value, sizeof bits);\n  bits = (bits << 32) | (bits >> 32);", HT + PT),

    # --- the hash of a path -------------------------------------------------
    ("path-version", "src", P, "#define PATH_FORMAT_VERSION 1u", "#define PATH_FORMAT_VERSION 2u", PT),
    ("path-no-dim", "src", P, "  gshgeom_hash_u32(&h, (uint32_t)path->dim);\n", "", PT),
    ("path-count-u32", "src", P, "gshgeom_hash_u64(&h, (uint64_t)path->el_n);", "gshgeom_hash_u32(&h, (uint32_t)path->el_n);", PT),
    ("path-no-verb", "src", P, "    gshgeom_hash_u8(&h, e->verb);\n", "", PT),
    ("path-close-hashes-a-point", "src", P, "    if (e->verb == GSHGEOM_VERB_CLOSE) {\n      continue;\n    }\n    for (size_t k = 0; k < verb_points(e->verb); k++) {\n      hash_point",
     "    for (size_t k = 0; k < verb_points(e->verb); k++) {\n      hash_point", PT),
    ("path-no-z", "src", P, "  if (dim == 3) {\n    gshgeom_hash_f64(h, p->z);", "  if (dim == 4) {\n    gshgeom_hash_f64(h, p->z);", PT),
    ("path-no-weights", "src", P, "      gshgeom_hash_f64(&h, path->wt[e->weight + k]);\n", "", PT),
    ("path-ends-only", "src", P, "    for (size_t k = 0; k < verb_points(e->verb); k++) {\n      hash_point(&h, &path->pt[e->point + k], path->dim);",
     "    for (size_t k = verb_points(e->verb) - 1; k < verb_points(e->verb); k++) {\n      hash_point(&h, &path->pt[e->point + k], path->dim);", PT),
    ("path-weights-first-only", "src", P, "    for (size_t k = 0; k < verb_weights(e->verb); k++) {\n      gshgeom_hash_f64",
     "    for (size_t k = 0; k < (verb_weights(e->verb) ? 1u : 0u); k++) {\n      gshgeom_hash_f64", PT),
    ("verb-code-quadratic", "hdr", PH, "GSHGEOM_VERB_QUADRATIC = 2,", "GSHGEOM_VERB_QUADRATIC = 9,", PT),
    ("verb-code-close", "hdr", PH, "GSHGEOM_VERB_CLOSE = 6,", "GSHGEOM_VERB_CLOSE = 8,", PT),
    ("kind-code-cubic", "hdr", CH, "GSHGEOM_CURVE_CUBIC = 2,", "GSHGEOM_CURVE_CUBIC = 6,", CT + PT),

    # --- the dump -----------------------------------------------------------
    ("dump-hash-width", "src", P, "hash=%016llx", "hash=%llx", PT),
    ("dump-number-width", "src", P, 'return fprintf(out, " %016llx"', 'return fprintf(out, " %llx"', PT),
    ("dump-no-z", "src", P, "if (n >= 0 && dim == 3) n = put_f64(out, p->z);", "if (n >= 0 && dim == 4) n = put_f64(out, p->z);", PT),
    ("dump-name-order", "src", P, '"M", "L", "Q", "C", "RQ", "RC", "Z"}', '"M", "L", "C", "Q", "RQ", "RC", "Z"}', PT),
    ("dump-no-weights", "src", P, "      for (size_t k = 0; k < verb_weights(e->verb); k++) {\n        if (put_f64",
     "      for (size_t k = 0; k < verb_weights(e->verb) / 5; k++) {\n        if (put_f64", PT),
    ("dump-ignores-flush", "src", P, "return (fflush(out) == 0 && !ferror(out)) ? GSHGEOM_OK : GSHGEOM_ERR_IO;", "return GSHGEOM_OK;", PT),
    ("dump-ignores-ferror", "src", P, "(fflush(out) == 0 && !ferror(out))", "(fflush(out) == 0)", PT),
    ("dump-ignores-write", "src", P, "(unsigned long long)gshgeom_path_hash(path)) < 0) {\n    return GSHGEOM_ERR_IO;",
     "(unsigned long long)gshgeom_path_hash(path)) < 0) {\n    return GSHGEOM_OK;", PT),

    # --- refusal ------------------------------------------------------------
    ("infinity-is-finite", "src", C, "return (v - v) == 0.0;", "return v == v;", CT + PT),
    ("limit-edge", "src", C, "if (v > GSHGEOM_COORD_LIMIT || v < -GSHGEOM_COORD_LIMIT) {", "if (v >= GSHGEOM_COORD_LIMIT || v < -GSHGEOM_COORD_LIMIT) {", CT),
    ("limit-negative-side", "src", C, "if (v > GSHGEOM_COORD_LIMIT || v < -GSHGEOM_COORD_LIMIT) {", "if (v > GSHGEOM_COORD_LIMIT) {", CT),
    ("weight-zero-ok", "src", C, "if (!is_finite(w) || !(w > 0.0)) {", "if (!is_finite(w) || !(w >= 0.0)) {", CT),
    ("weight-nonfinite-ok", "src", C, "if (!is_finite(w) || !(w > 0.0)) {", "if (!(w > 0.0)) {", CT),
    ("weight-max-unchecked", "src", C, "if (w < GSHGEOM_WEIGHT_MIN || w > GSHGEOM_WEIGHT_MAX) {", "if (w < GSHGEOM_WEIGHT_MIN) {", CT),
    ("weight-min-unchecked", "src", C, "if (w < GSHGEOM_WEIGHT_MIN || w > GSHGEOM_WEIGHT_MAX) {", "if (w > GSHGEOM_WEIGHT_MAX) {", CT),
    ("z-checked-in-2d", "src", C, "for (int a = 0; a < curve->dim; a++) {\n      double v = axis_of", "for (int a = 0; a < 3; a++) {\n      double v = axis_of", CT + PT),
    ("limit-hides-nan", "src", C, "if (v > GSHGEOM_COORD_LIMIT || v < -GSHGEOM_COORD_LIMIT) {\n        too_big = 1;",
     "if (v > GSHGEOM_COORD_LIMIT || v < -GSHGEOM_COORD_LIMIT) {\n        return GSHGEOM_ERR_LIMIT;", CT),
    ("z-kept-in-2d", "src", C, "    if (dim == 2) {\n      curve.p[i].z = 0.0;\n    }\n", "", CT + PT),
    ("dim-4-ok", "src", C, "(curve->dim != 2 && curve->dim != 3)", "(curve->dim < 2)", CT),
    ("dim-3-only", "src", C, "(curve->dim != 2 && curve->dim != 3)", "(curve->dim != 3)", CT),
    ("kind-count-ok", "src", C, "case GSHGEOM_CURVE_KIND_COUNT:\n    break;\n  }\n  return 0;", "case GSHGEOM_CURVE_KIND_COUNT:\n    return 4;\n  }\n  return 0;", CT),
    ("rational-flag-quadratic", "src", C, "return kind == GSHGEOM_CURVE_RATIONAL_QUADRATIC ||\n      kind == GSHGEOM_CURVE_RATIONAL_CUBIC;",
     "return kind == GSHGEOM_CURVE_RATIONAL_CUBIC;", CT),
    ("coord-limit-const", "hdr", CH, "#define GSHGEOM_COORD_LIMIT 0x1p100", "#define GSHGEOM_COORD_LIMIT 0x1p101", CT + PT),
    ("weight-max-const", "hdr", CH, "#define GSHGEOM_WEIGHT_MAX 0x1p100", "#define GSHGEOM_WEIGHT_MAX 0x1p99", CT),
    ("weight-min-const", "hdr", CH, "#define GSHGEOM_WEIGHT_MIN 0x1p-100", "#define GSHGEOM_WEIGHT_MIN 0x1p-99", CT),
    # --- evaluation ---------------------------------------------------------
    ("t-nan-ok", "src", C, "if (!out || !(t >= 0.0 && t <= 1.0)) {", "if (!out || t < 0.0 || t > 1.0) {", CT),
    ("t-one-refused", "src", C, "!(t >= 0.0 && t <= 1.0)", "!(t >= 0.0 && t < 1.0)", CT),
    ("t-zero-refused", "src", C, "!(t >= 0.0 && t <= 1.0)", "!(t > 0.0 && t <= 1.0)", CT),
    ("ends-swapped", "src", C, "int end = (t == 0.0) ? 0 : n - 1;", "int end = (t == 0.0) ? n - 1 : 0;", CT),
    ("no-exact-ends", "src", C, "if (t == 0.0 || t == 1.0) {", "if (t == 2.0) {", CT),
    ("lerp-wrong-operand", "src", C, "h[i][a] = s * h[i][a] + t * h[i + 1][a];", "h[i][a] = s * h[i][a] + t * h[i][a];", CT),
    ("lerp-s-t-swapped", "src", C, "h[i][a] = s * h[i][a] + t * h[i + 1][a];", "h[i][a] = t * h[i][a] + s * h[i + 1][a];", CT),
    ("weights-ignored-in-numerator", "src", C, "h[i][a] = work->rational ? work->w[i] * work->c[i][a] : work->c[i][a];",
     "h[i][a] = work->c[i][a];", CT),
    ("weight-not-blended", "src", C, "hw[i] = s * hw[i] + t * hw[i + 1];", "hw[i] = hw[i];", CT),
    ("no-division", "src", C, "out[a] = work->rational ? h[0][a] / hw[0] : h[0][a];", "out[a] = h[0][a];", CT),
    ("one-level-short", "src", C, "for (int level = n - 1; level > 0; level--) {", "for (int level = n - 1; level > 1; level--) {", CT),
    ("load-ignores-weights", "src", C, "work->w[i] = work->rational ? curve->w[i] : 1.0;", "work->w[i] = 1.0;", CT),

    # --- bounds -------------------------------------------------------------
    ("quadratic-no-extrema", "src", C, "if (work.n > 2) {", "if (work.n > 3) {", CT),
    ("cubic-no-extrema", "src", C, "if (work.n > 2) {", "if (work.n > 2 && work.n != 4) {", CT),
    ("no-end-candidate", "src", C, "  t[count++] = 1.0;\n  if (work.n > 2) {", "  if (work.n > 2) {", CT),
    ("max-flipped", "src", C, "if (v[a] > hi[a]) {", "if (v[a] < hi[a]) {", CT),
    ("min-flipped", "src", C, "if (v[a] < lo[a]) {", "if (v[a] > lo[a]) {", CT),
    ("only-x-extrema", "src", C, "for (int a = 0; a < work.dim; a++) {\n      count += axis_extrema", "for (int a = 0; a < 1; a++) {\n      count += axis_extrema", CT + PT),
    ("no-z-extrema", "src", C, "for (int a = 0; a < work.dim; a++) {\n      count += axis_extrema", "for (int a = 0; a < 2; a++) {\n      count += axis_extrema", CT),
    ("bisect-wrong-side", "src", C, "if ((fm < 0.0) == (fa < 0.0)) {\n      a = mid;", "if ((fm < 0.0) != (fa < 0.0)) {\n      a = mid;", CT),
    ("bisect-two-steps", "src", C, "for (int step = 0; step < 200; step++) {", "for (int step = 0; step < 2; step++) {", CT),
    ("linear-root", "src", C, "r[0] = c[0] / (c[0] - c[1]);", "r[0] = c[1] / (c[0] - c[1]);", CT),
    ("linear-needs-sign", "src", C, "if ((c[0] < 0.0 && c[1] > 0.0) || (c[0] > 0.0 && c[1] < 0.0)) {", "if (c[0] < 0.0 && c[1] > 0.0) {", CT),
    ("derivative-sum", "src", C, "d[k] = (double)m * (c[k + 1] - c[k]);", "d[k] = (double)m * (c[k + 1] + c[k]);", CT),
    ("critical-points-dropped", "src", C, "if (crit[i] > 0.0 && crit[i] < 1.0) {", "if (crit[i] > 2.0) {", CT),
    ("first-interval-only", "src", C, "for (int i = 0; i + 1 < k; i++) {\n    double a = pts[i];", "for (int i = 0; i + 1 < k && i < 1; i++) {\n    double a = pts[i];", CT),
    ("sign-change-test", "src", C, "} else if (fb != 0.0 && ((fa < 0.0) != (fb < 0.0))) {", "} else if (fb != 0.0 && ((fa < 0.0) == (fb < 0.0))) {", CT),
    ("choose-no-division", "src", C, "choose(n - 1, j) / choose(m, i + j))", "choose(n - 1, j))", CT),
    ("choose-row-4", "src", C, "{1, 4, 6, 4, 1}}", "{1, 4, 5, 4, 1}}", CT),
    ("choose-row-2", "src", C, "{1, 2, 1, 0, 0}", "{1, 3, 1, 0, 0}", CT),
    ("choose-row-1", "src", C, "{1, 1, 0, 0, 0},", "{1, 2, 0, 0, 0},", CT),
    ("term-difference-dropped", "src", C, "(work->c[j + 1][axis] - work->c[i][axis]);", "work->c[j + 1][axis];", CT),
    ("term-weight-dropped", "src", C, "double term = work->w[i] * work->w[j + 1] *", "double term = work->w[i] *", CT),
    ("term-weights-i-only", "src", C, "double term = work->w[i] * work->w[j + 1] *", "double term = work->w[i] * work->w[i] *", CT),
    ("loop-i-short", "src", C, "for (int i = 0; i < n; i++) {\n    for (int j = 0; j < n; j++) {\n      double term", "for (int i = 0; i < n - 1; i++) {\n    for (int j = 0; j < n; j++) {\n      double term", CT),
    ("loop-j-short", "src", C, "for (int j = 0; j < n; j++) {\n      double term", "for (int j = 0; j < n - 1; j++) {\n      double term", CT),
    ("degree-off", "src", C, "int m = 2 * n - 2;", "int m = 2 * n - 1;", CT),
    ("index-i-only", "src", C, "g[i + j] += (choose", "g[i] += (choose", CT),
    ("constructor-null-out", "src", C, "if (!out || (gshgeom_curve_is_rational(kind) && !w)) {", "if (gshgeom_curve_is_rational(kind) && !w) {", CT),
    ("candidates-not-evaluated", "src", C, "eval_work(&work, t[i], v);", "eval_work(&work, t[i] * 0.5, v);", CT),
    ("bounds-null-out", "src", C, "  if (!out) {\n    return GSHGEOM_ERR_INVALID;\n  }\n  Work work;\n  load(curve, &work);\n\n  double t[",
     "  Work work;\n  load(curve, &work);\n\n  double t[", CT),

    # --- the path -----------------------------------------------------------
    ("limit-off-by-one", "src", P, "if (path->el_n >= path->max_elements) {", "if (path->el_n > path->max_elements) {", PT),
    ("limit-off-by-one-low", "src", P, "if (path->el_n >= path->max_elements) {", "if (path->el_n + 1 >= path->max_elements) {", PT),
    ("zero-limit-is-zero", "src", P, "path->max_elements = (limits && limits->max_elements)\n      ? limits->max_elements\n      : GSHGEOM_LIMITS_DEFAULT_MAX_ELEMENTS;",
     "path->max_elements = limits ? limits->max_elements\n      : GSHGEOM_LIMITS_DEFAULT_MAX_ELEMENTS;", PT),
    ("default-limit-const", "hdr", KH, "#define GSHGEOM_LIMITS_DEFAULT_MAX_ELEMENTS ((size_t)1 << 22)", "#define GSHGEOM_LIMITS_DEFAULT_MAX_ELEMENTS ((size_t)1 << 21)", PT),
    ("close-without-move", "src", P, "GSHGEOM_Result gshgeom_path_close(GSHGEOM_Path * path) {\n  if (!path || !path->has_current) {",
     "GSHGEOM_Result gshgeom_path_close(GSHGEOM_Path * path) {\n  if (!path) {", PT),
    ("line-without-move", "src", P, "GSHGEOM_Result gshgeom_path_line_to(GSHGEOM_Path * path, GSHGEOM_Point p) {\n  if (!path || !path->has_current) {",
     "GSHGEOM_Result gshgeom_path_line_to(GSHGEOM_Path * path, GSHGEOM_Point p) {\n  if (!path) {", PT),
    ("cubic-without-move", "src", P, "GSHGEOM_Point c2, GSHGEOM_Point p) {\n  if (!path || !path->has_current) {", "GSHGEOM_Point c2, GSHGEOM_Point p) {\n  if (!path) {", PT),
    ("quadratic-without-move", "src", P, "GSHGEOM_Point cp, GSHGEOM_Point p) {\n  if (!path || !path->has_current) {", "GSHGEOM_Point cp, GSHGEOM_Point p) {\n  if (!path) {", PT),
    ("rational-quadratic-without-move", "src", P, "GSHGEOM_Point cp, GSHGEOM_Point p, const double w[3]) {\n  if (!path || !path->has_current) {",
     "GSHGEOM_Point cp, GSHGEOM_Point p, const double w[3]) {\n  if (!path) {", PT),
    ("rational-cubic-without-move", "src", P, "GSHGEOM_Point c1, GSHGEOM_Point c2, GSHGEOM_Point p, const double w[4]) {\n  if (!path || !path->has_current) {",
     "GSHGEOM_Point c1, GSHGEOM_Point c2, GSHGEOM_Point p, const double w[4]) {\n  if (!path) {", PT),
    ("move-forgets-start", "src", P, "    path->start = c.p[0];\n", "", PT),
    ("close-to-current", "src", P, "GSHGEOM_Point start = path->start;", "GSHGEOM_Point start = path->current;", PT),
    ("move-unchecked", "src", P, "  GSHGEOM_Result r = gshgeom_curve_line(path->dim, p, p, &c);\n  if (r != GSHGEOM_OK) {\n    return r;\n  }\n  r = append(path, GSHGEOM_VERB_MOVE",
     "  GSHGEOM_Result r = GSHGEOM_OK;\n  c.p[0] = p;\n  r = append(path, GSHGEOM_VERB_MOVE", PT),
    ("start-is-first-point", "src", P, "from = path->pt[prev->point + verb_points(prev->verb) - 1];", "from = path->pt[prev->point];", PT),
    ("start-off-by-one", "src", P, "from = path->pt[prev->point + verb_points(prev->verb) - 1];", "from = path->pt[prev->point + verb_points(prev->verb)];", PT),
    ("weights-from-first", "src", P, "c->w[i] = path->wt[e->weight + i];", "c->w[i] = path->wt[e->weight];", PT),
    ("close-kind", "src", P, "  case GSHGEOM_VERB_LINE:\n  case GSHGEOM_VERB_CLOSE:\n    c->kind = GSHGEOM_CURVE_LINE;",
     "  case GSHGEOM_VERB_LINE:\n    c->kind = GSHGEOM_CURVE_LINE;\n    break;\n  case GSHGEOM_VERB_CLOSE:\n    c->kind = GSHGEOM_CURVE_QUADRATIC;", PT),
    ("element-index-unchecked", "src", P, "if (!path || !out || index >= path->el_n) {", "if (!path || !out) {", PT),
    ("bounds-includes-move", "src", P, "if (verb == GSHGEOM_VERB_MOVE || verb == GSHGEOM_VERB_CLOSE) {", "if (verb == GSHGEOM_VERB_CLOSE) {", PT),
    ("bounds-includes-close", "src", P, "if (verb == GSHGEOM_VERB_MOVE || verb == GSHGEOM_VERB_CLOSE) {", "if (verb == GSHGEOM_VERB_MOVE) {", PT),
    ("bounds-min-x", "src", P, "if (b.min.x < box.min.x) box.min.x = b.min.x;", "if (b.min.x > box.min.x) box.min.x = b.min.x;", PT),
    ("bounds-min-z", "src", P, "if (b.min.z < box.min.z) box.min.z = b.min.z;", "if (b.min.z > box.min.z) box.min.z = b.min.z;", PT),
    ("bounds-max-y", "src", P, "if (b.max.y > box.max.y) box.max.y = b.max.y;", "if (b.max.y < box.max.y) box.max.y = b.max.y;", PT),
    ("bounds-max-x", "src", P, "if (b.max.x > box.max.x) box.max.x = b.max.x;", "if (b.max.x < box.max.x) box.max.x = b.max.x;", PT),
    ("bounds-first-only", "src", P, "    if (box.empty) {\n      box = b;\n    } else {", "    if (box.empty || 1) {\n      box = b;\n    } else {", PT),
    ("bounds-never-empty", "src", P, "box.empty = 1;\n  for", "box.empty = 0;\n  for", PT),
    ("create-dim-4", "src", P, "if (!out || (dim != 2 && dim != 3)) {", "if (!out || dim < 2) {", PT),
    ("create-no-out-check", "src", P, "if (!out || (dim != 2 && dim != 3)) {", "if (dim != 2 && dim != 3) {", PT),
    ("destroy-leaks-weights", "src", P, "  gcu_allocator_free(&alloc, path->wt);\n", "", PT),
    ("destroy-leaks-points", "src", P, "  gcu_allocator_free(&alloc, path->pt);\n", "", PT),
    ("destroy-leaks-elements", "src", P, "  gcu_allocator_free(&alloc, path->el);\n", "", PT),
    ("destroy-leaks-path", "src", P, "  gcu_allocator_free(&alloc, path);\n", "", PT),
    ("oom-as-internal", "src", P, "  void * bigger = gcu_allocator_realloc(alloc, array, grown * size);\n  if (!bigger) {\n    return GSHGEOM_ERR_OOM;",
     "  void * bigger = gcu_allocator_realloc(alloc, array, grown * size);\n  if (!bigger) {\n    return GSHGEOM_ERR_INTERNAL;", PT),
    ("create-oom-as-invalid", "src", P, "if (!path) {\n    return GSHGEOM_ERR_OOM;", "if (!path) {\n    return GSHGEOM_ERR_INVALID;", PT),
    ("verb-points-quadratic", "src", P, "  case GSHGEOM_VERB_QUADRATIC:\n  case GSHGEOM_VERB_RATIONAL_QUADRATIC:\n    return 2;",
     "  case GSHGEOM_VERB_QUADRATIC:\n  case GSHGEOM_VERB_RATIONAL_QUADRATIC:\n    return 3;", PT),
    ("verb-points-cubic", "src", P, "  case GSHGEOM_VERB_CUBIC:\n  case GSHGEOM_VERB_RATIONAL_CUBIC:\n    return 3;",
     "  case GSHGEOM_VERB_CUBIC:\n  case GSHGEOM_VERB_RATIONAL_CUBIC:\n    return 2;", PT),
    ("verb-weights-quadratic", "src", P, "return verb == GSHGEOM_VERB_RATIONAL_QUADRATIC ? 3", "return verb == GSHGEOM_VERB_RATIONAL_QUADRATIC ? 2", PT),
    ("verb-weights-cubic", "src", P, "verb == GSHGEOM_VERB_RATIONAL_CUBIC     ? 4", "verb == GSHGEOM_VERB_RATIONAL_CUBIC     ? 3", PT),
    ("hash-null-path", "src", P, "  if (!path) {\n    return gshgeom_hash_value(&h);\n  }", "  if (!path) {\n    return 0;\n  }", PT),

    # --- core ---------------------------------------------------------------
    ("string-duplicate", "src", K, 'return "invalid argument";', 'return "out of memory";', PT),
    ("string-empty", "src", K, 'return "a limit was exceeded";', 'return "";', PT),
    ("string-unknown-in-range", "src", K, 'case GSHGEOM_RESULT_COUNT:\n    break;', 'case GSHGEOM_RESULT_COUNT:\n    return "ok";', PT),
    ("limits-default-zero", "src", K, "limits->max_elements = GSHGEOM_LIMITS_DEFAULT_MAX_ELEMENTS;", "limits->max_elements = 0;", PT),
    ("limits-default-keeps-reserved", "src", K, "  *limits = (GSHGEOM_Limits){0};\n", "", PT),
    ("allocator-null", "src", "src/core/allocator.c", "return gcu_allocator_default();", "return NULL;", PT),
]


# --- story 6: the exact sum, the predicates, the maps -----------------------
X = "src/exact/exact.c"
O = "src/orient/orient.c"
M = "src/map/map.c"
MH = "include/ghoti.io/shape-geometry/map.h"
XT = ["testOrient", "testMap"]
OT = ["testOrient"]
MT = ["testMap"]


def kind_sets(allowed_old, extra):
    """The route predicate's allowed set with one more class let through."""
    return ("route-allows-" + extra.lower(), "src", M, allowed_old,
            allowed_old.replace("GSHGEOM_MAP_FLIPS_Y;", "GSHGEOM_MAP_FLIPS_Y | GSHGEOM_MAP_" + extra + ";"), MT)


ALLOWED = "GSHGEOM_MAP_TRANSLATES | GSHGEOM_MAP_SCALES | GSHGEOM_MAP_FLIPS_Y;"

PLANTS += [
    # --- exact sum: decomposition, product, accumulator, sign ---------------
    ("x-subnormal-exponent", "src", X, "*e = -1074;", "*e = -1075;", XT),
    ("x-no-implicit-bit", "src", X, "*m = frac | (UINT64_C(1) << 52);", "*m = frac;", XT),
    ("x-exponent-bias", "src", X, "*e = (int)ex - 1075;", "*e = (int)ex - 1074;", XT),
    ("x-factor-sign-ignored", "src", X, "p->neg ^= neg;", "(void)neg;", XT),
    ("x-term-sign-ignored", "src", X, "p->neg = t->neg ? 1 : 0;", "p->neg = 0;", XT),
    ("x-exponent-not-summed", "src", X, "p->e += e;", "p->e = e;", XT),
    ("x-first-high-word", "src", X, "v[1] = (uint32_t)(m >> 32);", "v[1] = 0;", XT),
    ("x-product-low-word-only", "src", X, "for (unsigned j = 0; j < 2; j++) {\n      uint64_t mj", "for (unsigned j = 0; j < 1; j++) {\n      uint64_t mj", XT),
    ("x-product-carry-dropped", "src", X, "out[len + j] = (uint32_t)carry;", "out[len + j] = 0;", XT),
    ("x-product-length", "src", X, "    len += 2;\n", "    len += 1;\n", XT),
    ("x-shift-no-spill", "src", X, "s[i] = (p->v[i] << r) | (p->v[i - 1] >> (32u - r));", "s[i] = (p->v[i] << r);", XT),
    ("x-shift-zero-lost", "src", X, "      s[i] = p->v[i];\n    }\n    s[6] = 0;", "      s[i] = 0;\n    }\n    s[6] = 0;", XT),
    ("x-add-carry-in-window", "src", X, "acc[q + i] = (uint32_t)(cur & 0xffffffffu);\n      carry = cur >> 32;", "acc[q + i] = (uint32_t)(cur & 0xffffffffu);\n      carry = 0;", XT),
    ("x-add-carry-above-window", "src", X, "carry && i < ACC_WORDS; i++) {\n      uint64_t cur = (uint64_t)acc[i] + carry;", "0 && i < ACC_WORDS; i++) {\n      uint64_t cur = (uint64_t)acc[i] + carry;", XT),
    ("x-sub-borrow-in-window", "src", X, "acc[q + i] = (uint32_t)(cur & 0xffffffffu);\n      borrow = (cur >> 63) & 1u;", "acc[q + i] = (uint32_t)(cur & 0xffffffffu);\n      borrow = 0;", XT),
    ("x-sub-borrow-above-window", "src", X, "acc[i] = (uint32_t)(cur & 0xffffffffu);\n      borrow = (cur >> 63) & 1u;", "acc[i] = (uint32_t)(cur & 0xffffffffu);\n      borrow = 0;", XT),
    ("x-minimum-exponent-is-maximum", "src", X, "(!have || ps[i].e < emin)", "(!have || ps[i].e > emin)", XT),
    ("x-negative-is-positive", "src", X, "    return -1;\n  }\n  for (unsigned i = 0; i < ACC_WORDS; i++) {", "    return 1;\n  }\n  for (unsigned i = 0; i < ACC_WORDS; i++) {", XT),
    ("x-word-zero-not-read", "src", X, "  for (unsigned i = 0; i < ACC_WORDS; i++) {\n    if (acc[i]) {", "  for (unsigned i = 1; i < ACC_WORDS; i++) {\n    if (acc[i]) {", XT),
    ("x-term-limit", "src", X, "count > GSHGEOM_EXACT_MAX_TERMS ||", "count >= GSHGEOM_EXACT_MAX_TERMS ||", XT),
    ("x-det3-sign-2", "src", X, "gshgeom_exact_put(t, n, !neg, 3, r1[0], r2[2], r3[1]);", "gshgeom_exact_put(t, n, neg, 3, r1[0], r2[2], r3[1]);", XT),
    ("x-det3-sign-4", "src", X, "gshgeom_exact_put(t, n, neg, 3, r1[1], r2[2], r3[0]);", "gshgeom_exact_put(t, n, !neg, 3, r1[1], r2[2], r3[0]);", XT),
    ("x-det3-sign-6", "src", X, "gshgeom_exact_put(t, n, !neg, 3, r1[2], r2[1], r3[0]);", "gshgeom_exact_put(t, n, neg, 3, r1[2], r2[1], r3[0]);", XT),
    ("x-det3-index", "src", X, "gshgeom_exact_put(t, n, neg, 3, r1[2], r2[0], r3[1]);", "gshgeom_exact_put(t, n, neg, 3, r1[2], r2[0], r3[0]);", XT),
    ("x-put-factors", "src", X, "t[*n].n = factors;", "(void)factors;\n  t[*n].n = 3;", XT),

    # --- predicates -----------------------------------------------------------
    ("o-ccw-bound-zero", "src", O, "#define CCW_ERR_BOUND ((3.0 + 16.0 * EPS) * EPS)", "#define CCW_ERR_BOUND 0.0", OT),
    ("o-ccw-bound-small", "src", O, "#define CCW_ERR_BOUND ((3.0 + 16.0 * EPS) * EPS)", "#define CCW_ERR_BOUND (0.25 * EPS)", OT),
    ("o-3d-bound-zero", "src", O, "#define O3D_ERR_BOUND ((7.0 + 56.0 * EPS) * EPS)", "#define O3D_ERR_BOUND 0.0", OT),
    ("o-3d-bound-small", "src", O, "#define O3D_ERR_BOUND ((7.0 + 56.0 * EPS) * EPS)", "#define O3D_ERR_BOUND (0.25 * EPS)", OT),
    ("o-filter-minimum", "src", O, "#define FILTER_MIN 0x1p-900", "#define FILTER_MIN 0.0", OT),
    ("o-filter-always-accepts-2d", "src", O, "(det >= errbound || -det >= errbound)", "(det >= errbound * 0.0 || -det >= errbound * 0.0)", OT),
    ("o-filter-always-accepts-3d", "src", O, "(det > errbound || -det > errbound)", "(det > errbound * 0.0 || -det > errbound * 0.0)", OT),
    ("o-filter-sign-2d", "src", O, "*sign = det > 0.0 ? 1 : -1;", "*sign = det > 0.0 ? -1 : 1;", OT),
    ("o-filter-sign-3d", "src", O, "*sign = det > 0.0 ? -1 : 1;", "*sign = det > 0.0 ? 1 : -1;", OT),
    ("o-detsum-difference", "src", O, "double detsum = absolute(detleft) + absolute(detright);", "double detsum = absolute(detleft) - absolute(detright);", OT),
    ("o-permanent-third-term", "src", O, "(absolute(adxbdy) + absolute(bdxady)) * absolute(cdz)", "0.0", OT),
    ("o-permanent-first-term", "src", O, "(absolute(bdxcdy) + absolute(cdxbdy)) * absolute(adz)", "0.0", OT),
    ("o-3d-determinant", "src", O, "cdz * (adxbdy - bdxady)", "cdz * (adxbdy + bdxady)", OT),
    ("o-2d-term-1", "src", O, "put(t, &n, 1, a.x, c.y);", "put(t, &n, 0, a.x, c.y);", OT),
    ("o-2d-term-5", "src", O, "put(t, &n, 0, b.x, c.y);", "put(t, &n, 1, b.x, c.y);", OT),
    ("o-2d-term-6", "src", O, "put(t, &n, 1, b.y, c.x);", "put(t, &n, 1, b.y, c.y);", OT),
    ("o-3d-sign-bcd", "src", O, "put_det3(t, &n, 0, b, c, d);", "put_det3(t, &n, 1, b, c, d);", OT),
    ("o-3d-sign-acd", "src", O, "put_det3(t, &n, 1, a, c, d);", "put_det3(t, &n, 0, a, c, d);", OT),
    ("o-3d-sign-abd", "src", O, "put_det3(t, &n, 0, a, b, d);", "put_det3(t, &n, 1, a, b, d);", OT),
    ("o-3d-sign-abc", "src", O, "put_det3(t, &n, 1, a, b, c);", "put_det3(t, &n, 0, a, b, c);", OT),
    ("o-3d-rows", "src", O, "const double c[3] = {r3.x, r3.y, r3.z};", "const double c[3] = {r3.x, r3.y, r3.y};", OT),
    ("o-2d-no-finite-check", "src", O, "if (!sign || !finite2(a) || !finite2(b) || !finite2(c)) {\n    return GSHGEOM_ERR_INVALID;\n  }\n  int s;", "if (!sign) {\n    return GSHGEOM_ERR_INVALID;\n  }\n  int s;", OT),
    ("o-2d-exact-no-finite-check", "src", O, "if (!sign || !finite2(a) || !finite2(b) || !finite2(c)) {\n    return GSHGEOM_ERR_INVALID;\n  }\n  *sign = exact2", "if (!sign) {\n    return GSHGEOM_ERR_INVALID;\n  }\n  *sign = exact2", OT),
    ("o-3d-no-finite-check", "src", O, "if (!sign || !finite3(a) || !finite3(b) || !finite3(c) || !finite3(d)) {\n    return GSHGEOM_ERR_INVALID;\n  }\n  int s;", "if (!sign) {\n    return GSHGEOM_ERR_INVALID;\n  }\n  int s;", OT),
    ("o-3d-exact-no-finite-check", "src", O, "if (!sign || !finite3(a) || !finite3(b) || !finite3(c) || !finite3(d)) {\n    return GSHGEOM_ERR_INVALID;\n  }\n  *sign = exact3", "if (!sign) {\n    return GSHGEOM_ERR_INVALID;\n  }\n  *sign = exact3", OT),
    ("o-finite3-no-z", "src", O, "return is_finite(p.x) && is_finite(p.y) && is_finite(p.z);", "return is_finite(p.x) && is_finite(p.y);", OT),
    ("o-finite-2d-reads-nan", "src", O, "return (v - v) == 0.0;", "return v == v;", OT),
    ("o-normal-left", "src", O, "GSHGEOM_Point n = {tangent.y, 0.0 - tangent.x, 0.0};", "GSHGEOM_Point n = {0.0 - tangent.y, tangent.x, 0.0};", OT),
    ("o-normal-right", "src", O, "GSHGEOM_Point n = {0.0 - tangent.y, tangent.x, 0.0};", "GSHGEOM_Point n = {tangent.y, 0.0 - tangent.x, 0.0};", OT),
    ("o-normal-negative-zero", "src", O, "GSHGEOM_Point n = {tangent.y, 0.0 - tangent.x, 0.0};", "GSHGEOM_Point n = {tangent.y, -tangent.x, 0.0};", OT),
    ("o-normal-keeps-z", "src", O, "GSHGEOM_Point n = {0.0 - tangent.y, tangent.x, 0.0};", "GSHGEOM_Point n = {0.0 - tangent.y, tangent.x, tangent.z};", OT),
    ("o-area-sign", "src", O, "sum += ux * vy - uy * vx;", "sum += uy * vx - ux * vy;", OT),
    ("o-area-half", "src", O, "double area = sum * 0.5;", "double area = sum;", OT),
    ("o-area-fan-start", "src", O, "for (size_t i = 1; i + 1 < count; i++) {", "for (size_t i = 2; i + 1 < count; i++) {", OT),
    ("o-area-fan-end", "src", O, "for (size_t i = 1; i + 1 < count; i++) {", "for (size_t i = 1; i + 2 < count; i++) {", OT),
    ("o-area-from-origin", "src", O, "double ux = points[i].x - points[0].x;", "double ux = points[i].x;", OT),
    ("o-area-no-limit", "src", O, "if (!is_finite(area)) {\n    return GSHGEOM_ERR_LIMIT;", "if (0) {\n    return GSHGEOM_ERR_LIMIT;", OT),
    ("o-area-no-finite-check", "src", O, "    if (!finite2(points[i])) {\n      return GSHGEOM_ERR_INVALID;", "    if (0) {\n      return GSHGEOM_ERR_INVALID;", OT),
    ("o-area-no-null-check", "src", O, "if (!out || (count && !points)) {", "if (!out) {", OT),

    # --- maps: constants and enumerations -------------------------------------
    ("m-limit-bound", "hdr", MH, "#define GSHGEOM_MAP_LIMIT 0x1p100", "#define GSHGEOM_MAP_LIMIT 0x1p101", MT),
    ("m-quarters-bound-high", "hdr", MH, "#define GSHGEOM_MAP_MAX_QUARTERS 0x1p40", "#define GSHGEOM_MAP_MAX_QUARTERS 0x1p41", MT),
    ("m-quarters-bound-low", "hdr", MH, "#define GSHGEOM_MAP_MAX_QUARTERS 0x1p40", "#define GSHGEOM_MAP_MAX_QUARTERS 0x1p39", MT),
    ("m-unit-gradians", "hdr", MH, "GSHGEOM_ANGLE_GRADIANS = 2,", "GSHGEOM_ANGLE_GRADIANS = 3,", MT),
    ("m-unit-degrees", "hdr", MH, "GSHGEOM_ANGLE_DEGREES = 1,", "GSHGEOM_ANGLE_DEGREES = 2,", MT),
    ("m-axis-z", "hdr", MH, "GSHGEOM_AXIS_Z = 2,", "GSHGEOM_AXIS_Z = 1,", MT),
    ("m-kind-translates", "hdr", MH, "GSHGEOM_MAP_TRANSLATES = 1u << 0,", "GSHGEOM_MAP_TRANSLATES = 1u << 9,", MT),
    ("m-kind-flips-y", "hdr", MH, "GSHGEOM_MAP_FLIPS_Y = 1u << 3,", "GSHGEOM_MAP_FLIPS_Y = 1u << 4,", MT),
    ("m-kind-quarter", "hdr", MH, "GSHGEOM_MAP_QUARTER_TURNS = 1u << 5,", "GSHGEOM_MAP_QUARTER_TURNS = 1u << 6,", MT),
    ("m-kind-all", "hdr", MH, "GSHGEOM_MAP_KIND_ALL = (1u << 10) - 1u", "GSHGEOM_MAP_KIND_ALL = (1u << 9) - 1u", MT),
    ("m-rad-per-degree", "src", M, "#define RAD_PER_DEGREE 0x1.1df46a2529d39p-6", "#define RAD_PER_DEGREE 0x1.1df46a2529d00p-6", MT),
    ("m-rad-per-gradian", "src", M, "#define RAD_PER_GRADIAN 0x1.015bf9217271ap-6", "#define RAD_PER_GRADIAN 0x1.015bf92172700p-6", MT),
    ("m-rad-per-turn", "src", M, "#define RAD_PER_TURN 0x1.921fb54442d18p+2", "#define RAD_PER_TURN 0x1.921fb54442d00p+2", MT),
    ("m-format-version", "src", M, "#define MAP_FORMAT_VERSION 1u", "#define MAP_FORMAT_VERSION 2u", MT),

    # --- maps: angles ---------------------------------------------------------
    ("m-degree-quarter", "src", M, "quarter = 90.0;", "quarter = 89.0;", MT),
    ("m-gradian-quarter", "src", M, "quarter = 100.0;", "quarter = 99.0;", MT),
    ("m-turn-quarter", "src", M, "u = angle * 4.0;", "u = angle * 2.0;", MT),
    ("m-turn-quarter-size", "src", M, "quarter = 0.25;", "quarter = 0.5;", MT),
    ("m-nearest-quarter-floor", "src", M, "double nr = gmath_round(u);", "double nr = gmath_floor(u);", MT),
    ("m-no-reduction", "src", M, "double rest = angle - nr * quarter;", "double rest = angle;", MT),
    ("m-quarter-mod", "src", M, "double kq = nr - 4.0 * gmath_floor(nr * 0.25);", "double kq = nr - 4.0 * gmath_trunc(nr * 0.25);", MT),
    ("m-quarter-mod-period", "src", M, "double kq = nr - 4.0 * gmath_floor(nr * 0.25);", "double kq = nr - 2.0 * gmath_floor(nr * 0.5);", MT),
    ("m-quarter-limit", "src", M, "if (nr > GSHGEOM_MAP_MAX_QUARTERS || nr < -GSHGEOM_MAP_MAX_QUARTERS) {", "if (nr > GSHGEOM_MAP_MAX_QUARTERS) {", MT),
    ("m-quarter-1", "src", M, "*c = 0.0 - ss;\n    *s = cc;", "*c = ss;\n    *s = cc;", MT),
    ("m-quarter-2", "src", M, "*c = 0.0 - cc;\n    *s = 0.0 - ss;", "*c = 0.0 - cc;\n    *s = ss;", MT),
    ("m-quarter-3", "src", M, "*c = ss;\n    *s = 0.0 - cc;", "*c = ss;\n    *s = cc;", MT),
    ("m-rest-zero-exact", "src", M, "    cc = 1.0;\n    ss = 0.0;\n  } else {\n    double theta", "    cc = 1.0;\n    ss = 1.0;\n  } else {\n    double theta", MT),
    ("m-angle-not-finite", "src", M, "|| !is_finite(angle)) {", ") {", MT),
    ("m-unit-range", "src", M, "(unsigned)unit >= (unsigned)GSHGEOM_ANGLE_UNIT_COUNT", "(unsigned)unit > (unsigned)GSHGEOM_ANGLE_UNIT_COUNT", MT),
    ("m-radians-sine-cosine", "src", M, "*c = gmath_cos(angle);\n      *s = gmath_sin(angle);", "*c = gmath_sin(angle);\n      *s = gmath_cos(angle);", MT),
    ("m-rotation-sense", "src", M, "double m[3][3] = {{c, 0.0 - s, 0}, {s, c, 0}, {0, 0, 1}};", "double m[3][3] = {{c, s, 0}, {0.0 - s, c, 0}, {0, 0, 1}};", MT),
    ("m-axis-range", "src", M, "(unsigned)axis >= (unsigned)GSHGEOM_AXIS_COUNT", "(unsigned)axis > (unsigned)GSHGEOM_AXIS_COUNT", MT),
    ("m-axis-x-sense", "src", M, "m[1][2] = 0.0 - s;\n    m[2][1] = s;", "m[1][2] = s;\n    m[2][1] = 0.0 - s;", MT),
    ("m-axis-y-sense", "src", M, "m[2][0] = 0.0 - s;\n    m[0][2] = s;", "m[2][0] = s;\n    m[0][2] = 0.0 - s;", MT),
    ("m-axis-z-sense", "src", M, "m[0][1] = 0.0 - s;\n    m[1][0] = s;", "m[0][1] = s;\n    m[1][0] = 0.0 - s;", MT),
    ("m-axis-branches", "src", M, "  } else if (axis == GSHGEOM_AXIS_X) {", "  } else if (axis == GSHGEOM_AXIS_Y) {", MT),

    # --- maps: construction and the classifier -------------------------------
    ("m-zero-sign", "src", M, "return v == 0.0 ? 0.0 : v;", "return v;", MT),
    ("m-dim-range", "src", M, "return dim == 2 || dim == 3;", "return dim >= 2 && dim <= 4;", MT),
    ("m-number-limit", "src", M, "if (v[i] > GSHGEOM_MAP_LIMIT || v[i] < -GSHGEOM_MAP_LIMIT) {", "if (v[i] >= GSHGEOM_MAP_LIMIT || v[i] <= -GSHGEOM_MAP_LIMIT) {", MT),
    ("m-number-lower-limit", "src", M, "if (v[i] > GSHGEOM_MAP_LIMIT || v[i] < -GSHGEOM_MAP_LIMIT) {", "if (v[i] > GSHGEOM_MAP_LIMIT) {", MT),
    ("m-fill-unused-diagonal", "src", M, "(i == j ? 1.0 : 0.0);\n    }\n    r.t[i]", "(i == j ? 0.0 : 0.0);\n    }\n    r.t[i]", MT),
    ("m-fill-unused-translation", "src", M, "r.t[i] = i < dim ? nz(t[i]) : 0.0;", "r.t[i] = i < dim ? nz(t[i]) : t[i];", MT),
    ("m-from-matrix-transposed", "src", M, "mm[i][j] = m[i * dim + j];", "mm[i][j] = m[j * dim + i];", MT),
    ("m-from-matrix-translation", "src", M, "tt[i] = t[i];", "tt[i] = t[0];", MT),
    ("m-translates-bit", "src", M, "if (t[i] != 0.0) {\n      kind |= GSHGEOM_MAP_TRANSLATES;", "if (t[i] == 0.0) {\n      kind |= GSHGEOM_MAP_TRANSLATES;", MT),
    ("m-singular-sign", "src", M, "  if (sign == 0) {\n    return kind | GSHGEOM_MAP_SINGULAR;", "  if (sign == 1) {\n    return kind | GSHGEOM_MAP_SINGULAR;", MT),
    ("m-det2-second-term", "src", M, "gshgeom_exact_put(t, &n, 1, 2, m[0][1], m[1][0], 0.0);", "gshgeom_exact_put(t, &n, 0, 2, m[0][1], m[1][0], 0.0);", MT),
    ("m-mono-row-count", "src", M, "    if (nonzero != 1) {\n      mono = 0;", "    if (nonzero > 2) {\n      mono = 0;", MT),
    ("m-mono-column", "src", M, "      if (col[i] != i) {\n        identity_perm = 0;", "      if (col[i] == i) {\n        identity_perm = 0;", MT),
    ("m-negative-count", "src", M, "      if (v < 0.0) {\n        negatives++;", "      if (v > 0.0) {\n        negatives++;", MT),
    ("m-y-row", "src", M, "        if (i == 1) {\n          y_negative = 1;", "        if (i == 0) {\n          y_negative = 1;", MT),
    ("m-flip-or-mirror", "src", M, "kind |= y_negative ? GSHGEOM_MAP_FLIPS_Y : GSHGEOM_MAP_MIRRORS;", "kind |= y_negative ? GSHGEOM_MAP_MIRRORS : GSHGEOM_MAP_FLIPS_Y;", MT),
    ("m-half-turn", "src", M, "kind |= GSHGEOM_MAP_QUARTER_TURNS; /* a half turn */", "kind |= GSHGEOM_MAP_MIRRORS; /* a half turn */", MT),
    ("m-point-reflection", "src", M, "} else if (negatives == 3) {", "} else if (negatives == 4) {", MT),
    ("m-axis-swap-mirror", "src", M, "      if (sign < 0) {\n        kind |= GSHGEOM_MAP_MIRRORS;", "      if (sign > 0) {\n        kind |= GSHGEOM_MAP_MIRRORS;", MT),
    ("m-axis-swap-quarter", "src", M, "    } else {\n      kind |= GSHGEOM_MAP_QUARTER_TURNS;\n      if (sign < 0) {", "    } else {\n      kind |= GSHGEOM_MAP_SCALES;\n      if (sign < 0) {", MT),
    ("m-magnitude-sign", "src", M, "mag[i] = absolute(v);", "mag[i] = v;", MT),
    ("m-uniform-test", "src", M, "if (mag[i] != mag[0]) {", "if (mag[i] == mag[0]) {", MT),
    ("m-unit-scale", "src", M, "} else if (mag[0] != 1.0) {", "} else if (mag[0] != 2.0) {", MT),
    ("m-nonuniform-bit", "src", M, "kind |= GSHGEOM_MAP_SCALES_NONUNIFORM;", "kind |= GSHGEOM_MAP_SCALES;", MT),
    ("m-rotation-test", "src", M, "if (a == d && b == -c) {", "if (a == d && b == c) {", MT),
    ("m-mirrored-rotation-test", "src", M, "} else if (a == -d && b == c) {", "} else if (a == -d && b == -c) {", MT),
    ("m-skew-bit", "src", M, "      kind |= GSHGEOM_MAP_SKEWS;\n    }", "      kind |= GSHGEOM_MAP_ROTATES;\n    }", MT),
    ("m-general-3d", "src", M, "return kind | GSHGEOM_MAP_GENERAL;", "return kind;", MT),

    # --- maps: check, composition, application ------------------------------
    ("m-check-kind-bits", "src", M, "(map->kind & ~(uint32_t)GSHGEOM_MAP_KIND_ALL)", "0", MT),
    ("m-check-translation", "src", M, "  for (int i = 0; i < dim; i++) {\n    v[n++] = map->t[i];\n  }", "  for (int i = 0; i < 0; i++) {\n    v[n++] = map->t[i];\n  }", MT),
    ("m-check-matrix-rows", "src", M, "  for (int i = 0; i < dim; i++) {\n    for (int j = 0; j < dim; j++) {\n      v[n++] = map->m[i][j];", "  for (int i = 0; i < dim; i++) {\n    for (int j = 0; j < 1; j++) {\n      v[n++] = map->m[i][j];", MT),
    ("m-compose-dimension", "src", M, "if (outer->dim != inner->dim) {", "if (0) {", MT),
    ("m-compose-operands", "src", M, "double acc = outer->m[i][0] * inner->m[0][j];", "double acc = inner->m[i][0] * outer->m[0][j];", MT),
    ("m-compose-association", "src", M, "acc = acc + outer->m[i][2] * inner->m[2][j];", "acc = outer->m[i][0] * inner->m[0][j] + (outer->m[i][1] * inner->m[1][j] + outer->m[i][2] * inner->m[2][j]);", MT),
    ("m-compose-no-third-term", "src", M, "      if (dim == 3) {\n        acc = acc + outer->m[i][2] * inner->m[2][j];", "      if (dim == 4) {\n        acc = acc + outer->m[i][2] * inner->m[2][j];", MT),
    ("m-compose-translation-order", "src", M, "tv = tv + outer->m[i][1] * inner->t[1];", "tv = outer->m[i][1] * inner->t[1];", MT),
    ("m-compose-translation-add", "src", M, "res.t[i] = nz(tv + outer->t[i]);", "res.t[i] = nz(tv);", MT),
    ("m-compose-translation-sign", "src", M, "res.t[i] = nz(tv + outer->t[i]);", "res.t[i] = tv + outer->t[i];", MT),
    ("m-compose-entry-sign", "src", M, "res.m[i][j] = nz(acc);", "res.m[i][j] = acc;", MT),
    ("m-compose-outer-kind-only", "src", M, "res.kind = outer->kind | inner->kind;", "res.kind = outer->kind;", MT),
    ("m-compose-inner-kind-only", "src", M, "res.kind = outer->kind | inner->kind;", "res.kind = inner->kind;", MT),
    ("m-compose-kind-intersection", "src", M, "res.kind = outer->kind | inner->kind;", "res.kind = outer->kind & inner->kind;", MT),
    ("m-compose-collapse", "src", M, "if (det_sign(dim, res.m) == 0) {\n    res.kind |= GSHGEOM_MAP_SINGULAR;", "if (0) {\n    res.kind |= GSHGEOM_MAP_SINGULAR;", MT),
    ("m-compose-result-limit", "src", M, "  res.kind = outer->kind | inner->kind;\n  r = gshgeom_map_check(&res);", "  res.kind = outer->kind | inner->kind;\n  r = GSHGEOM_OK;", MT),
    ("m-chain-leaf-first", "src", M, "r = compose_one(&acc, &chain[i], &acc);", "r = compose_one(&chain[i], &acc, &acc);", MT),
    ("m-chain-right-to-left", "src", M, "    for (size_t i = 1; i < count; i++) {\n      r = compose_one(&acc, &chain[i], &acc);", "    acc = chain[count - 1];\n    for (size_t i = count - 1; i-- > 0;) {\n      r = compose_one(&chain[i], &acc, &acc);", MT),
    ("m-chain-target-first", "src", M, "r = compose_one(target, &acc, &acc);", "r = compose_one(&acc, target, &acc);", MT),
    ("m-chain-target-ignored", "src", M, "    if (target) {\n      r = compose_one(target, &acc, &acc);", "    if (0) {\n      r = compose_one(target, &acc, &acc);", MT),
    ("m-chain-empty", "src", M, "(!count && !target)", "(!count && target)", MT),
    ("m-design-scale-sign", "src", M, "r = gshgeom_map_scale(2, scale, scale, 1.0, &chain[1]);", "r = gshgeom_map_scale(2, scale, -scale, 1.0, &chain[1]);", MT),
    ("m-design-scale-positive", "src", M, "|| !(scale > 0.0)) {", "|| !(scale >= 0.0)) {", MT),
    ("m-design-order", "src", M, "r = gshgeom_map_translate(2, tx, ty, 0.0, &chain[0]);", "r = gshgeom_map_translate(2, ty, tx, 0.0, &chain[0]);", MT),
    ("m-apply-translation", "src", M, "    if (translate) {\n      v = v + map->t[i];", "    if (!translate) {\n      v = v + map->t[i];", MT),
    ("m-apply-association", "src", M, "    v = v + map->m[i][1] * in[1];\n    if (map->dim == 3) {\n      v = v + map->m[i][2] * in[2];\n    }", "    if (map->dim == 3) {\n      v = v + (map->m[i][1] * in[1] + map->m[i][2] * in[2]);\n    } else {\n      v = v + map->m[i][1] * in[1];\n    }", MT),
    ("m-apply-no-limit", "src", M, "    if (!is_finite(v)) {\n      return GSHGEOM_ERR_LIMIT;", "    if (0) {\n      return GSHGEOM_ERR_LIMIT;", MT),
    ("m-apply-reads-z-in-2d", "src", M, "(map->dim == 3 && !is_finite(p.z))", "!is_finite(p.z)", MT),
    ("m-apply-no-output-check", "src", M, "if (!out || !is_finite(p.x)", "if (!is_finite(p.x)", MT),
    ("m-curve-point-count", "src", M, "int n = gshgeom_curve_point_count(curve->kind);", "int n = gshgeom_curve_point_count(curve->kind) - 1;", MT),
    ("m-curve-weights", "src", M, "  GSHGEOM_Curve res = *curve;\n", "  GSHGEOM_Curve res = *curve;\n  res.w[1] = 1.0;\n", MT),
    ("m-curve-dimension", "src", M, "if (curve->dim != map->dim) {", "if (0) {", MT),
    ("m-curve-result-check", "src", M, "r = gshgeom_curve_check(&res);", "r = GSHGEOM_OK;", MT),

    # --- maps: the route predicate -------------------------------------------
    kind_sets("      " + ALLOWED, "QUARTER_TURNS"),
    kind_sets("      " + ALLOWED, "SCALES_NONUNIFORM"),
    kind_sets("      " + ALLOWED, "MIRRORS"),
    kind_sets("      " + ALLOWED, "ROTATES"),
    kind_sets("      " + ALLOWED, "SKEWS"),
    kind_sets("      " + ALLOWED, "SINGULAR"),
    kind_sets("      " + ALLOWED, "GENERAL"),
    ("route-excludes-flip", "src", M, "GSHGEOM_MAP_TRANSLATES | GSHGEOM_MAP_SCALES | GSHGEOM_MAP_FLIPS_Y;", "GSHGEOM_MAP_TRANSLATES | GSHGEOM_MAP_SCALES;", MT),
    ("route-excludes-scale", "src", M, "GSHGEOM_MAP_TRANSLATES | GSHGEOM_MAP_SCALES | GSHGEOM_MAP_FLIPS_Y;", "GSHGEOM_MAP_TRANSLATES | GSHGEOM_MAP_FLIPS_Y;", MT),
    ("route-excludes-translation", "src", M, "GSHGEOM_MAP_TRANSLATES | GSHGEOM_MAP_SCALES | GSHGEOM_MAP_FLIPS_Y;", "GSHGEOM_MAP_SCALES | GSHGEOM_MAP_FLIPS_Y;", MT),
    ("route-dimension", "src", M, "  if (map->dim != 2) {\n    return GSHGEOM_OK;", "  if (map->dim != 3) {\n    return GSHGEOM_OK;", MT),
    ("route-ignores-tag", "src", M, "  if (map->kind & ~allowed) {\n    return GSHGEOM_OK;", "  if (map->kind & ~allowed & 0u) {\n    return GSHGEOM_OK;", MT),
    ("route-positive-scale", "src", M, "!(a > 0.0)) {", "!(a >= 0.0)) {", MT),
    ("route-negative-scale", "src", M, "!(a > 0.0)) {", "a == 0.0) {", MT),
    ("route-skew-upper", "src", M, "if (map->m[0][1] != 0.0 || map->m[1][0] != 0.0 ||", "if (map->m[1][0] != 0.0 ||", MT),
    ("route-skew-lower", "src", M, "if (map->m[0][1] != 0.0 || map->m[1][0] != 0.0 ||", "if (map->m[0][1] != 0.0 ||", MT),
    ("route-no-flip", "src", M, "  if (d == a) {\n    flip = 0;", "  if (d == a) {\n    flip = 1;", MT),
    ("route-flip-flag", "src", M, "  } else if (d == -a) {\n    flip = 1;", "  } else if (d == -a) {\n    flip = 0;", MT),
    ("route-flip-magnitude", "src", M, "  } else if (d == -a) {", "  } else if (d * d == a * a) {", MT),
    ("route-aligned-value", "src", M, "  *aligned = 1;\n  if (grid) {", "  *aligned = 2;\n  if (grid) {", MT),
    ("route-not-aligned-value", "src", M, "  *aligned = 0;\n  if (map->dim != 2) {", "  if (map->dim != 2) {", MT),
    ("route-grid-scale", "src", M, "grid->scale = a;", "grid->scale = d;", MT),
    ("route-grid-tx", "src", M, "grid->tx = map->t[0];", "grid->tx = map->t[1];", MT),
    ("route-grid-ty", "src", M, "grid->ty = map->t[1];", "grid->ty = map->t[0];", MT),
    ("route-null-aligned", "src", M, "  if (!aligned) {\n    return GSHGEOM_ERR_INVALID;\n  }\n  GSHGEOM_Result r = gshgeom_map_check(map);", "  GSHGEOM_Result r = gshgeom_map_check(map);", MT),

    # --- maps: hash and dump ---------------------------------------------------
    ("mh-no-dimension", "src", M, "  gshgeom_hash_u32(&h, (uint32_t)map->dim);\n", "", MT),
    ("mh-no-kind", "src", M, "  gshgeom_hash_u32(&h, map->kind);\n", "", MT),
    ("mh-kind-u64", "src", M, "gshgeom_hash_u32(&h, map->kind);", "gshgeom_hash_u64(&h, map->kind);", MT),
    ("mh-column-major", "src", M, "gshgeom_hash_f64(&h, map->m[i][j]);", "gshgeom_hash_f64(&h, map->m[j][i]);", MT),
    ("mh-two-columns", "src", M, "    for (int j = 0; j < map->dim; j++) {\n      gshgeom_hash_f64(&h", "    for (int j = 0; j < 2; j++) {\n      gshgeom_hash_f64(&h", MT),
    ("mh-translation", "src", M, "gshgeom_hash_f64(&h, map->t[i]);", "gshgeom_hash_f64(&h, map->t[0]);", MT),
    ("mh-invalid-map", "src", M, "  if (gshgeom_map_check(map) != GSHGEOM_OK) {\n    return gshgeom_hash_value(&h);", "  if (0) {\n    return gshgeom_hash_value(&h);", MT),
    ("md-kind-width", "src", M, "kind=%08llx", "kind=%llx", MT),
    ("md-hash-width", "src", M, "hash=%016llx", "hash=%llx", MT),
    ("md-row-letter", "src", M, "if (fputc('m', out) == EOF) {", "if (fputc('M', out) == EOF) {", MT),
    ("md-translation-letter", "src", M, "if (fputc('t', out) == EOF) {", "if (fputc('T', out) == EOF) {", MT),
    ("md-entry-width", "src", M, 'fprintf(out, " %016llx", (unsigned long long)gshgeom_f64_bits(map->m[i][j]))', 'fprintf(out, " %llx", (unsigned long long)gshgeom_f64_bits(map->m[i][j]))', MT),
    ("md-translation-width", "src", M, 'fprintf(out, " %016llx", (unsigned long long)gshgeom_f64_bits(map->t[i]))', 'fprintf(out, " %llx", (unsigned long long)gshgeom_f64_bits(map->t[i]))', MT),
    ("md-ignores-flush", "src", M, "return (fflush(out) == 0 && !ferror(out)) ? GSHGEOM_OK : GSHGEOM_ERR_IO;", "return GSHGEOM_OK;", MT),
    ("md-ignores-ferror", "src", M, "(fflush(out) == 0 && !ferror(out))", "(fflush(out) == 0)", MT),
    ("md-invalid-map", "src", M, "if (!out || gshgeom_map_check(map) != GSHGEOM_OK) {", "if (!out) {", MT),
]


# Plants caught by the probe's hash alone.
FLAG_PLANTS = [
    ("fused-multiply-add", "-mfma -ffp-contract=fast"),
]


# --- story 7: the rounding, the device path, the flattener ---------------------
DV = "src/device/device.c"
FL = "src/flatten/flatten.c"
GR = "src/core/grow.h"
DH = "include/ghoti.io/shape-geometry/device.h"
FH = "include/ghoti.io/shape-geometry/flatten.h"
DT = ["testDevice"]
FT = ["testFlatten", "testTolerance"]
DFT = ["testDevice", "testFlatten"]

PLANTS += [
    # --- the rounding rule ----------------------------------------------------
    ("r-max-bound", "src", DV, "#define ROUND_T_MAX 2147483647.5", "#define ROUND_T_MAX 2147483648.5", DT),
    ("r-max-bound-low", "src", DV, "#define ROUND_T_MAX 2147483647.5", "#define ROUND_T_MAX 2147483646.5", DT),
    ("r-min-bound", "src", DV, "#define ROUND_T_MIN (-2147483648.5)", "#define ROUND_T_MIN (-2147483649.5)", DT),
    ("r-min-bound-high", "src", DV, "#define ROUND_T_MIN (-2147483648.5)", "#define ROUND_T_MIN (-2147483647.5)", DT),
    ("r-max-inclusive", "src", DV, "if (!(t < ROUND_T_MAX) || t < ROUND_T_MIN) {", "if (!(t <= ROUND_T_MAX) || t < ROUND_T_MIN) {", DT),
    ("r-min-inclusive", "src", DV, "if (!(t < ROUND_T_MAX) || t < ROUND_T_MIN) {", "if (!(t < ROUND_T_MAX) || t <= ROUND_T_MIN) {", DT),
    ("r-max-unchecked", "src", DV, "if (!(t < ROUND_T_MAX) || t < ROUND_T_MIN) {", "if (t < ROUND_T_MIN) {", DT),
    ("r-min-unchecked", "src", DV, "if (!(t < ROUND_T_MAX) || t < ROUND_T_MIN) {", "if (!(t < ROUND_T_MAX)) {", DT),
    ("r-half-goes-down", "src", DV, "double r = (t - f) >= 0.5 ? f + 1.0 : f;", "double r = (t - f) > 0.5 ? f + 1.0 : f;", DT),
    ("r-half-always-up", "src", DV, "double r = (t - f) >= 0.5 ? f + 1.0 : f;", "double r = (t - f) >= 0.25 ? f + 1.0 : f;", DT),
    ("r-double-expression", "src", DV, "double r = (t - f) >= 0.5 ? f + 1.0 : f;", "double r = f * 0.0 + gmath_floor(t + 0.5);", DT),
    ("r-half-away-from-zero", "src", DV, "double r = (t - f) >= 0.5 ? f + 1.0 : f;", "double r = f * 0.0 + gmath_round(t);", DT),
    ("r-truncated", "src", DV, "double f = gmath_floor(t);", "double f = (double)(int64_t)t;", DT),
    ("r-ceiling", "src", DV, "double r = (t - f) >= 0.5 ? f + 1.0 : f;", "double r = (t - f) > 0.0 ? f + 1.0 : f;", DT),
    ("r-scale", "src", DV, "return round_scaled(x * SCALE, out);", "return round_scaled(x * 32.0, out);", DT),
    ("r-writes-on-limit", "src", DV, "    return GSHGEOM_ERR_LIMIT;\n  }\n  double f = gmath_floor(t);",
     "    *out = 0;\n    return GSHGEOM_ERR_LIMIT;\n  }\n  double f = gmath_floor(t);", DT),
    ("r-nan-not-refused", "src", DV, "if (!out || !((x - x) == 0.0)) {", "if (!out) {", DT),
    ("r-out-unchecked", "src", DV, "if (!out || !((x - x) == 0.0)) {", "if (!((x - x) == 0.0)) {", DT),
    ("r-infinity-as-finite", "src", DV, "if (!out || !((x - x) == 0.0)) {", "if (!out || x != x) {", DT),
    ("r-device-one", "hdr", DH, "#define GSHGEOM_DEVICE_ONE 64", "#define GSHGEOM_DEVICE_ONE 32", DT),

    # --- the point function ---------------------------------------------------
    ("dp-no-2d-check", "src", DV, "  if (map->dim != 2) {\n    return GSHGEOM_ERR_INVALID;\n  }\n  GSHGEOM_Point m;",
     "  GSHGEOM_Point m;", DT),
    ("dp-y-from-x", "src", DV, "r = gshgeom_round_26_6(m.y, &d.y);", "r = gshgeom_round_26_6(m.x, &d.y);", DT),
    ("dp-x-from-y", "src", DV, "r = gshgeom_round_26_6(m.x, &d.x);", "r = gshgeom_round_26_6(m.y, &d.x);", DT),
    ("dp-writes-on-failure", "src", DV, "  if (r == GSHGEOM_OK) {\n    *out = d;\n  }\n  return r;\n}", "  *out = d;\n  return r;\n}", DT),
    ("dp-null-out", "src", DV, "GSHGEOM_Device_Point * out) {\n  if (!map || !out) {", "GSHGEOM_Device_Point * out) {\n  if (!map) {", DT),

    # --- the device path builder ----------------------------------------------
    ("b-verb-move", "hdr", DH, "GSHGEOM_DEVICE_MOVE = 0,", "GSHGEOM_DEVICE_MOVE = 8,", DT),
    ("b-verb-line", "hdr", DH, "GSHGEOM_DEVICE_LINE = 1,", "GSHGEOM_DEVICE_LINE = 8,", DT),
    ("b-verb-quadratic", "hdr", DH, "GSHGEOM_DEVICE_QUADRATIC = 2,", "GSHGEOM_DEVICE_QUADRATIC = 7,", DT),
    ("b-verb-cubic", "hdr", DH, "GSHGEOM_DEVICE_CUBIC = 3,", "GSHGEOM_DEVICE_CUBIC = 6,", DT),
    ("b-verb-close", "hdr", DH, "GSHGEOM_DEVICE_CLOSE = 4,", "GSHGEOM_DEVICE_CLOSE = 9,", DT),
    ("b-limit-allows-one-more", "src", DV, "add_el > path->max_elements - path->el_n) {", "add_el > path->max_elements - path->el_n + 1) {", DT),
    ("b-limit-one-short", "src", DV, "add_el > path->max_elements - path->el_n) {", "add_el >= path->max_elements - path->el_n) {", DT),
        ("b-move-counts-one", "src", DV, "size_t add_el = with_move ? 2 : 1;", "size_t add_el = 1;", DT),
    ("b-zero-limit-is-zero", "src", DV, "path->max_elements = (limits && limits->max_elements)\n      ? limits->max_elements\n      : GSHGEOM_LIMITS_DEFAULT_MAX_ELEMENTS;",
     "path->max_elements = limits ? limits->max_elements\n      : GSHGEOM_LIMITS_DEFAULT_MAX_ELEMENTS;", DT),
    ("b-close-flag-never", "src", DV, "path->closed = verb == GSHGEOM_DEVICE_CLOSE;", "path->closed = 0;", DT),
    ("b-close-flag-sticky", "src", DV, "path->closed = verb == GSHGEOM_DEVICE_CLOSE;", "path->closed = path->closed || verb == GSHGEOM_DEVICE_CLOSE;", DT),
    ("b-close-flag-not-move", "src", DV, "path->closed = verb == GSHGEOM_DEVICE_CLOSE;", "path->closed = verb != GSHGEOM_DEVICE_MOVE;", DT),
    ("b-move-forgets-start", "src", DV, "    path->has_current = 1;\n    path->start = p;", "    path->has_current = 1;", DT),
    ("b-close-to-current", "src", DV, "GSHGEOM_Device_Point start = path->start;\n  return append(path, GSHGEOM_DEVICE_CLOSE", "GSHGEOM_Device_Point start = path->pt[path->pt_n - 1];\n  return append(path, GSHGEOM_DEVICE_CLOSE", DT),
    ("b-stored-move-is-origin", "src", DV, "    path->pt[path->pt_n] = path->start;", "    path->pt[path->pt_n] = pts[0];", DT),
    ("b-stored-move-verb", "src", DV, "    path->el[path->el_n].verb = GSHGEOM_DEVICE_MOVE;\n    path->el[path->el_n].point = path->pt_n;\n    path->pt[path->pt_n] = path->start;",
     "    path->el[path->el_n].verb = GSHGEOM_DEVICE_LINE;\n    path->el[path->el_n].point = path->pt_n;\n    path->pt[path->pt_n] = path->start;", DT),
    ("b-line-without-move", "src", DV, "GSHGEOM_Result gshgeom_device_path_line_to(\n    GSHGEOM_Device_Path * path, GSHGEOM_Device_Point p) {\n  if (!path || !path->has_current) {",
     "GSHGEOM_Result gshgeom_device_path_line_to(\n    GSHGEOM_Device_Path * path, GSHGEOM_Device_Point p) {\n  if (!path) {", DT),
    ("b-quadratic-without-move", "src", DV, "GSHGEOM_Device_Point c, GSHGEOM_Device_Point p) {\n  if (!path || !path->has_current) {", "GSHGEOM_Device_Point c, GSHGEOM_Device_Point p) {\n  if (!path) {", DT),
    ("b-cubic-without-move", "src", DV, "GSHGEOM_Device_Point c2, GSHGEOM_Device_Point p) {\n  if (!path || !path->has_current) {", "GSHGEOM_Device_Point c2, GSHGEOM_Device_Point p) {\n  if (!path) {", DT),
    ("b-close-without-move", "src", DV, "GSHGEOM_Result gshgeom_device_path_close(GSHGEOM_Device_Path * path) {\n  if (!path || !path->has_current) {",
     "GSHGEOM_Result gshgeom_device_path_close(GSHGEOM_Device_Path * path) {\n  if (!path) {", DT),
    ("b-line-no-stored-move", "src", DV, "return append(path, GSHGEOM_DEVICE_LINE, &p, path->closed);", "return append(path, GSHGEOM_DEVICE_LINE, &p, 0);", DT),
    ("b-quadratic-no-stored-move", "src", DV, "return append(path, GSHGEOM_DEVICE_QUADRATIC, pts, path->closed);", "return append(path, GSHGEOM_DEVICE_QUADRATIC, pts, 0);", DT),
    ("b-cubic-no-stored-move", "src", DV, "return append(path, GSHGEOM_DEVICE_CUBIC, pts, path->closed);", "return append(path, GSHGEOM_DEVICE_CUBIC, pts, 0);", DT),
    ("b-close-stored-move", "src", DV, "return append(path, GSHGEOM_DEVICE_CLOSE, &start, 0);", "return append(path, GSHGEOM_DEVICE_CLOSE, &start, path->closed);", DT),
    ("b-move-stored-move", "src", DV, "GSHGEOM_Result r = append(path, GSHGEOM_DEVICE_MOVE, &p, 0);", "GSHGEOM_Result r = append(path, GSHGEOM_DEVICE_MOVE, &p, path->closed);", DT),
    ("b-quadratic-order", "src", DV, "GSHGEOM_Device_Point pts[2] = {c, p};", "GSHGEOM_Device_Point pts[2] = {p, c};", DT),
    ("b-cubic-order", "src", DV, "GSHGEOM_Device_Point pts[3] = {c1, c2, p};", "GSHGEOM_Device_Point pts[3] = {c2, c1, p};", DT),
    ("b-points-line", "src", DV, "  case GSHGEOM_DEVICE_MOVE:\n  case GSHGEOM_DEVICE_LINE:\n  case GSHGEOM_DEVICE_CLOSE:\n    return 1;",
     "  case GSHGEOM_DEVICE_MOVE:\n  case GSHGEOM_DEVICE_LINE:\n  case GSHGEOM_DEVICE_CLOSE:\n    return 2;", DT),
    ("b-points-quadratic", "src", DV, "  case GSHGEOM_DEVICE_QUADRATIC:\n    return 2;", "  case GSHGEOM_DEVICE_QUADRATIC:\n    return 3;", DT),
    ("b-points-cubic", "src", DV, "  case GSHGEOM_DEVICE_CUBIC:\n    return 3;", "  case GSHGEOM_DEVICE_CUBIC:\n    return 2;", DT),
    ("b-readback-points", "src", DV, "    el.p[i] = path->pt[e->point + i];", "    el.p[i] = path->pt[e->point];", DT),
    ("b-readback-count", "src", DV, "  el.count = verb_points(e->verb);", "  el.count = 1;", DT),
    ("b-readback-verb", "src", DV, "  el.verb = (GSHGEOM_Device_Verb)e->verb;", "  el.verb = GSHGEOM_DEVICE_LINE;", DT),
    ("b-readback-unzeroed", "src", DV, "  GSHGEOM_Device_Element el;\n  memset(&el, 0, sizeof el);", "  GSHGEOM_Device_Element el;\n  memset(&el, 0xff, sizeof el);", DT),
    ("b-element-index", "src", DV, "if (!path || !out || index >= path->el_n) {", "if (!path || !out) {", DT),
    ("b-element-index-off-by-one", "src", DV, "if (!path || !out || index >= path->el_n) {", "if (!path || !out || index > path->el_n) {", DT),
    ("b-count-null", "src", DV, "return path ? path->el_n : 0;", "return path->el_n;", DT),
    ("b-destroy-leaks-elements", "src", DV, "  gcu_allocator_free(&alloc, path->el);\n  gcu_allocator_free(&alloc, path->pt);\n  gcu_allocator_free(&alloc, path);",
     "  gcu_allocator_free(&alloc, path->pt);\n  gcu_allocator_free(&alloc, path);", DT),
    ("b-destroy-leaks-points", "src", DV, "  gcu_allocator_free(&alloc, path->el);\n  gcu_allocator_free(&alloc, path->pt);\n  gcu_allocator_free(&alloc, path);",
     "  gcu_allocator_free(&alloc, path->el);\n  gcu_allocator_free(&alloc, path);", DT),
    ("b-destroy-leaks-path", "src", DV, "  gcu_allocator_free(&alloc, path->el);\n  gcu_allocator_free(&alloc, path->pt);\n  gcu_allocator_free(&alloc, path);",
     "  gcu_allocator_free(&alloc, path->el);\n  gcu_allocator_free(&alloc, path->pt);", DT),
    ("b-create-oom-code", "src", DV, "  if (!path) {\n    return GSHGEOM_ERR_OOM;\n  }\n  path->alloc = *alloc;", "  if (!path) {\n    return GSHGEOM_ERR_INVALID;\n  }\n  path->alloc = *alloc;", DT),
    ("b-create-null-out", "src", DV, "GSHGEOM_Device_Path ** out) {\n  if (!out) {\n    return GSHGEOM_ERR_INVALID;\n  }\n  if (!alloc) {", "GSHGEOM_Device_Path ** out) {\n  if (!alloc) {", DT),
    ("b-elements-grow-oom-code", "src", DV, "  GSHGEOM_Result r = gshgeom_grow(&path->alloc, &grown, &path->el_cap,\n      path->el_n + add_el, sizeof *path->el);\n  if (r != GSHGEOM_OK) {\n    return r;\n  }",
     "  GSHGEOM_Result r = gshgeom_grow(&path->alloc, &grown, &path->el_cap,\n      path->el_n + add_el, sizeof *path->el);\n  if (r != GSHGEOM_OK) {\n    return GSHGEOM_ERR_INTERNAL;\n  }", DT),
    ("b-points-grow-unchecked", "src", DV, "  r = gshgeom_grow(&path->alloc, &grown, &path->pt_cap, path->pt_n + add_pt,\n      sizeof *path->pt);\n  if (r != GSHGEOM_OK) {\n    return r;\n  }",
     "  r = gshgeom_grow(&path->alloc, &grown, &path->pt_cap, path->pt_n + add_pt,\n      sizeof *path->pt);", DT),
    ("b-grow-keeps-old-pointer", "srch", GR, "  *array = bigger;\n", "", DFT),
    ("b-grow-oom-ok", "srch", GR, "  if (!bigger) {\n    return GSHGEOM_ERR_OOM;\n  }\n  *array = bigger;", "  if (!bigger) {\n    return GSHGEOM_OK;\n  }\n  *array = bigger;", DFT),
        # --- the device path's hash and dump --------------------------------------
    ("h-version", "src", DV, "#define DEVICE_FORMAT_VERSION 1u", "#define DEVICE_FORMAT_VERSION 2u", DT),
    ("h-kind", "src", DV, "#define DEVICE_KIND 0x50564544u", "#define DEVICE_KIND 0x50564545u", DT),
    ("h-no-kind", "src", DV, "  gshgeom_hash_u32(&h, DEVICE_KIND);\n", "", DT),
    ("h-count-u32", "src", DV, "gshgeom_hash_u64(&h, (uint64_t)path->el_n);", "gshgeom_hash_u32(&h, (uint32_t)path->el_n);", DT),
    ("h-no-verb", "src", DV, "    gshgeom_hash_u8(&h, e->verb);\n    if (e->verb == GSHGEOM_DEVICE_CLOSE) {", "    if (e->verb == GSHGEOM_DEVICE_CLOSE) {", DT),
    ("h-close-hashes-a-point", "src", DV, "    if (e->verb == GSHGEOM_DEVICE_CLOSE) {\n      continue;\n    }\n    for (unsigned k = 0; k < verb_points(e->verb); k++) {\n      gshgeom_hash_i32",
     "    for (unsigned k = 0; k < verb_points(e->verb); k++) {\n      gshgeom_hash_i32", DT),
    ("h-y-for-x", "src", DV, "      gshgeom_hash_i32(&h, path->pt[e->point + k].x);", "      gshgeom_hash_i32(&h, path->pt[e->point + k].y);", DT),
    ("h-no-y", "src", DV, "      gshgeom_hash_i32(&h, path->pt[e->point + k].y);\n    }\n  }\n  return gshgeom_hash_value(&h);\n}\n\nstatic int put_i32", "    }\n  }\n  return gshgeom_hash_value(&h);\n}\n\nstatic int put_i32", DT),
    ("h-last-point-only", "src", DV, "    for (unsigned k = 0; k < verb_points(e->verb); k++) {\n      gshgeom_hash_i32", "    for (unsigned k = verb_points(e->verb) - 1; k < verb_points(e->verb); k++) {\n      gshgeom_hash_i32", DT),
    ("h-null-path", "src", DV, "  if (!path) {\n    return gshgeom_hash_value(&h);\n  }\n  gshgeom_hash_u32(&h, DEVICE_FORMAT_VERSION);", "  if (!path) {\n    return 0;\n  }\n  gshgeom_hash_u32(&h, DEVICE_FORMAT_VERSION);", DT),
    ("h-dump-header", "src", DV, '"gshgeom.device-path %u elements=%llu hash=%016llx\\n"', '"gshgeom.device_path %u elements=%llu hash=%016llx\\n"', DT),
    ("h-dump-hash-width", "src", DV, "elements=%llu hash=%016llx", "elements=%llu hash=%llx", DT),
    ("h-dump-letters", "src", DV, "{'M', 'L', 'Q', 'C', 'Z'}", "{'M', 'L', 'C', 'Q', 'Z'}", DT),
    ("h-dump-number-width", "src", DV, 'return fprintf(out, " %08x", (unsigned)(uint32_t)v);', 'return fprintf(out, " %x", (unsigned)(uint32_t)v);', DT),
    ("h-dump-signed", "src", DV, 'return fprintf(out, " %08x", (unsigned)(uint32_t)v);', 'return fprintf(out, " %08d", (int)v);', DT),
    ("h-dump-close-points", "src", DV, "    if (e->verb != GSHGEOM_DEVICE_CLOSE) {\n      for (unsigned k = 0; k < verb_points(e->verb); k++) {\n        if (put_i32",
     "    if (1) {\n      for (unsigned k = 0; k < verb_points(e->verb); k++) {\n        if (put_i32", DT),
    ("h-dump-ignores-flush", "src", DV, "  return (fflush(out) == 0 && !ferror(out)) ? GSHGEOM_OK : GSHGEOM_ERR_IO;\n}\n", "  return GSHGEOM_OK;\n}\n", DT + ["testFlatten"][:0]),
    ("h-dump-ignores-ferror", "src", DV, "(fflush(out) == 0 && !ferror(out))", "(fflush(out) == 0)", DT),
    ("h-dump-ignores-header-error", "src", DV, "(unsigned long long)gshgeom_device_path_hash(path)) < 0) {\n    return GSHGEOM_ERR_IO;", "(unsigned long long)gshgeom_device_path_hash(path)) < 0) {\n    return GSHGEOM_OK;", DT),
        # --- a model path under a map ---------------------------------------------
    ("fm-3d-path-ok", "src", DV, "if (!path || !map || !out || gshgeom_path_dimension(path) != 2) {", "if (!path || !map || !out) {", DT),
    ("fm-null-out", "src", DV, "if (!path || !map || !out || gshgeom_path_dimension(path) != 2) {", "if (!path || !map || gshgeom_path_dimension(path) != 2) {", DT),
        ("fm-3d-map-ok", "src", DV, "  if (map->dim != 2) {\n    return GSHGEOM_ERR_INVALID;\n  }\n  GSHGEOM_Device_Path * dp;", "  GSHGEOM_Device_Path * dp;", DT),
    ("fm-skips-first-element", "src", DV, "for (size_t i = 0; i < n && r == GSHGEOM_OK; i++) {\n    GSHGEOM_Element el;", "for (size_t i = 1; i < n && r == GSHGEOM_OK; i++) {\n    GSHGEOM_Element el;", DT),
    ("fm-skips-last-element", "src", DV, "for (size_t i = 0; i < n && r == GSHGEOM_OK; i++) {\n    GSHGEOM_Element el;", "for (size_t i = 0; i + 1 < n && r == GSHGEOM_OK; i++) {\n    GSHGEOM_Element el;", DT),
    ("fm-leaks-on-failure", "src", DV, "  if (r != GSHGEOM_OK) {\n    gshgeom_device_path_destroy(dp);\n    return r;\n  }\n  *out = dp;", "  if (r != GSHGEOM_OK) {\n    return r;\n  }\n  *out = dp;", DT),
    ("fm-line-from-start", "src", DV, "    r = map_rounded(map, c->p[1], &q[0]);\n    return r != GSHGEOM_OK ? r : gshgeom_device_path_line_to(dp, q[0]);",
     "    r = map_rounded(map, c->p[0], &q[0]);\n    return r != GSHGEOM_OK ? r : gshgeom_device_path_line_to(dp, q[0]);", DT),
    ("fm-quadratic-swapped", "src", DV, "gshgeom_device_path_quadratic_to(dp, q[0], q[1]);", "gshgeom_device_path_quadratic_to(dp, q[1], q[0]);", DT),
    ("fm-quadratic-end-twice", "src", DV, "      r = map_rounded(map, c->p[2], &q[1]);\n    }\n    return r != GSHGEOM_OK ? r : gshgeom_device_path_quadratic_to",
     "      r = map_rounded(map, c->p[1], &q[1]);\n    }\n    return r != GSHGEOM_OK ? r : gshgeom_device_path_quadratic_to", DT),
    ("fm-cubic-swapped", "src", DV, "gshgeom_device_path_cubic_to(dp, q[0], q[1], q[2]);", "gshgeom_device_path_cubic_to(dp, q[0], q[2], q[1]);", DT),
    ("fm-cubic-end", "src", DV, "      r = map_rounded(map, c->p[3], &q[2]);", "      r = map_rounded(map, c->p[2], &q[2]);", DT),
    ("fm-close-as-line", "src", DV, "  case GSHGEOM_VERB_CLOSE:\n    return gshgeom_device_path_close(dp);", "  case GSHGEOM_VERB_CLOSE:\n    return GSHGEOM_OK;", DT),
    ("fm-scale", "src", DV, "double sx = m.x * SCALE;", "double sx = m.x * 32.0;", DT),
    ("fm-scale-y", "src", DV, "double sy = m.y * SCALE;", "double sy = m.y * 32.0;", DT),
    ("fm-y-range-unchecked", "src", DV, "  r = round_scaled(sx, &probe);\n  if (r == GSHGEOM_OK) {\n    r = round_scaled(sy, &probe);\n  }", "  r = round_scaled(sx, &probe);", DT),
    ("fm-x-range-unchecked", "src", DV, "  r = round_scaled(sx, &probe);\n  if (r == GSHGEOM_OK) {\n    r = round_scaled(sy, &probe);\n  }", "  r = round_scaled(sy, &probe);", DT),
        # --- the rational subdivision ---------------------------------------------
    ("rat-tolerance-loose", "src", DV, "#define RATIONAL_TOLERANCE_SQ 4.0", "#define RATIONAL_TOLERANCE_SQ 400.0", DT),
    ("rat-tolerance-tight", "src", DV, "#define RATIONAL_TOLERANCE_SQ 4.0", "#define RATIONAL_TOLERANCE_SQ 0.0001", DT),
    ("rat-depth-shallow", "src", DV, "#define RATIONAL_MAX_DEPTH 16", "#define RATIONAL_MAX_DEPTH 2", DT),
    ("rat-always-flat", "src", DV, "  if (depth <= 0 || flat) {\n    GSHGEOM_Device_Point d;", "  (void)flat;\n  if (1) {\n    GSHGEOM_Device_Point d;", DT),
    ("rat-never-flat", "src", DV, "  if (depth <= 0 || flat) {\n    GSHGEOM_Device_Point d;", "  (void)flat;\n  if (depth <= 0) {\n    GSHGEOM_Device_Point d;", DT),
    ("rat-flat-inverted", "src", DV, "if (dist2_to_segment(px, py, sx, sy, ex, ey) > RATIONAL_TOLERANCE_SQ) {", "if (dist2_to_segment(px, py, sx, sy, ex, ey) < RATIONAL_TOLERANCE_SQ) {", DT),
    ("rat-first-control-only", "src", DV, "  for (int i = 1; i < n - 1; i++) {\n    double px", "  for (int i = 1; i < 2; i++) {\n    double px", DT),
    ("rat-last-control-only", "src", DV, "  for (int i = 1; i < n - 1; i++) {\n    double px", "  for (int i = n - 2; i < n - 1; i++) {\n    double px", DT),
    ("rat-controls-unprojected", "src", DV, "    double px = h[i].x / h[i].w;\n    double py = h[i].y / h[i].w;", "    double px = h[i].x;\n    double py = h[i].y;", DT),
    ("rat-y-unprojected", "src", DV, "    double py = h[i].y / h[i].w;", "    double py = h[i].y;", DT),
    ("rat-distance-to-line", "src", DV, "  if (u < 0.0) {\n    u = 0.0;\n  } else if (u > 1.0) {\n    u = 1.0;\n  }\n", "", DT),
    ("rat-clamp-low-only", "src", DV, "  } else if (u > 1.0) {\n    u = 1.0;\n  }\n", "  }\n", DT),
    ("rat-clamp-high-only", "src", DV, "  if (u < 0.0) {\n    u = 0.0;\n  } else if (u > 1.0) {", "  if (u > 1.0) {", DT),
    ("rat-degenerate-chord", "src", DV, "  if (!(len2 > 0.0)) {\n    return ax * ax + ay * ay;\n  }", "  if (len2 < 0.0) {\n    return ax * ax + ay * ay;\n  }", DT),
    ("rat-degenerate-chord-flat", "src", DV, "  if (!(len2 > 0.0)) {\n    return ax * ax + ay * ay;\n  }", "  if (!(len2 > 0.0)) {\n    return 0.0;\n  }", DT),
    ("rat-dist-no-x", "src", DV, "  double qx = ax - u * dx;\n  double qy = ay - u * dy;\n  return qx * qx + qy * qy;", "  double qx = 0.0;\n  double qy = ay - u * dy;\n  return qx * qx + qy * qy;", DT),
    ("rat-mid-weight", "src", DV, "  m.w = 0.5 * (a.w + b.w);", "  m.w = 0.5 * (a.w + a.w);", DT),
    ("rat-mid-x", "src", DV, "  m.x = 0.5 * (a.x + b.x);", "  m.x = 0.5 * (a.x + a.x);", DT),
    ("rat-mid-y", "src", DV, "  m.y = 0.5 * (a.y + b.y);", "  m.y = 0.5 * (b.y + b.y);", DT),
    ("rat-quadratic-left", "src", DV, "    left[0] = h[0]; left[1] = a; left[2] = m;", "    left[0] = h[0]; left[1] = b; left[2] = m;", DT),
    ("rat-quadratic-right", "src", DV, "    right[0] = m; right[1] = b; right[2] = h[2];", "    right[0] = m; right[1] = a; right[2] = h[2];", DT),
    ("rat-quadratic-mid-point", "src", DV, "    double mx = m.x / m.w;\n    double my = m.y / m.w;\n    GSHGEOM_Result r = rational_lines(dp, 3,", "    double mx = m.x;\n    double my = m.y / m.w;\n    GSHGEOM_Result r = rational_lines(dp, 3,", DT),
    ("rat-quadratic-order", "src", DV, "    return r != GSHGEOM_OK\n        ? r\n        : rational_lines(dp, 3, right, mx, my, ex, ey, depth - 1);",
     "    return r != GSHGEOM_OK\n        ? r\n        : rational_lines(dp, 3, right, sx, sy, ex, ey, depth - 1);", DT),
    ("rat-cubic-left", "src", DV, "left[0] = h[0]; left[1] = a; left[2] = d; left[3] = m;", "left[0] = h[0]; left[1] = a; left[2] = e; left[3] = m;", DT),
    ("rat-cubic-right", "src", DV, "right[0] = m; right[1] = e; right[2] = c; right[3] = h[3];", "right[0] = m; right[1] = d; right[2] = c; right[3] = h[3];", DT),
    ("rat-cubic-mid", "src", DV, "  Hp m = hp_mid(d, e);", "  Hp m = hp_mid(a, c);", DT),
    ("rat-cubic-d", "src", DV, "  Hp d = hp_mid(a, b);\n  Hp e = hp_mid(b, c);", "  Hp d = hp_mid(a, c);\n  Hp e = hp_mid(b, c);", DT),
    ("rat-cubic-order", "src", DV, "      : rational_lines(dp, 4, right, mx, my, ex, ey, depth - 1);\n}", "      : rational_lines(dp, 4, right, sx, sy, ex, ey, depth - 1);\n}", DT),
    ("rat-cubic-error-lost", "src", DV, "  GSHGEOM_Result r = rational_lines(dp, 4, left, sx, sy, mx, my, depth - 1);\n  return r != GSHGEOM_OK\n      ? r\n      : rational_lines(dp, 4, right",
     "  rational_lines(dp, 4, left, sx, sy, mx, my, depth - 1);\n  GSHGEOM_Result r = GSHGEOM_OK;\n  return r != GSHGEOM_OK\n      ? r\n      : rational_lines(dp, 4, right", DT),
    ("rat-quadratic-error-lost", "src", DV, "    GSHGEOM_Result r = rational_lines(dp, 3, left, sx, sy, mx, my, depth - 1);\n    return r != GSHGEOM_OK\n        ? r\n        : rational_lines(dp, 3, right",
     "    rational_lines(dp, 3, left, sx, sy, mx, my, depth - 1);\n    GSHGEOM_Result r = GSHGEOM_OK;\n    return r != GSHGEOM_OK\n        ? r\n        : rational_lines(dp, 3, right", DT),
    ("rat-end-y", "src", DV, "r = round_scaled(ey, &d.y);", "r = round_scaled(ex, &d.y);", DT),
    ("rat-end-x", "src", DV, "    GSHGEOM_Result r = round_scaled(ex, &d.x);", "    GSHGEOM_Result r = round_scaled(ey, &d.x);", DT),
    ("rat-weights-ignored-x", "src", DV, "      h[i].x = c->w[i] * sx[i];", "      h[i].x = sx[i];", DT),
    ("rat-weights-ignored-y", "src", DV, "      h[i].y = c->w[i] * sy[i];", "      h[i].y = sy[i];", DT),
    ("rat-weight-one", "src", DV, "      h[i].w = c->w[i];", "      h[i].w = 1.0;", DT),
    ("rat-kind-swapped", "src", DV, "int n = el->verb == GSHGEOM_VERB_RATIONAL_QUADRATIC ? 3 : 4;", "int n = el->verb == GSHGEOM_VERB_RATIONAL_QUADRATIC ? 4 : 3;", DT),
    ("rat-last-control-end", "src", DV, "sx[0], sy[0], sx[n - 1], sy[n - 1],", "sx[0], sy[0], sx[n - 2], sy[n - 2],", DT),
    ("rat-start-end", "src", DV, "sx[0], sy[0], sx[n - 1], sy[n - 1],", "sx[1], sy[1], sx[n - 1], sy[n - 1],", DT),
    ("rat-depth-start", "src", DV, "        RATIONAL_MAX_DEPTH);\n  }", "        0);\n  }", DT),
    ("rat-map-error-lost", "src", DV, "      r = map_scaled(map, c->p[i], &sx[i], &sy[i]);\n", "      map_scaled(map, c->p[i], &sx[i], &sy[i]);\n", DT),
    ("rat-skips-last-point", "src", DV, "for (int i = 0; i < n && r == GSHGEOM_OK; i++) {\n      r = map_scaled", "for (int i = 0; i < n - 1 && r == GSHGEOM_OK; i++) {\n      r = map_scaled", DT),

    # --- the flattener --------------------------------------------------------
    # The floor division and the midpoint (revision 2: widened, half up).
    ("fl-floordiv-truncates", "src", FL, "  if (a % b != 0 && a < 0) {\n    q -= 1;\n  }", "", FT),
    ("fl-floordiv-floors-positive", "src", FL, "  if (a % b != 0 && a < 0) {", "  if (a % b != 0 && a > 0) {", FT),
    ("fl-mid-floors-x", "src", FL, "  m.x = floor_div(a.x + b.x + 1, 2);", "  m.x = floor_div(a.x + b.x, 2);", FT),
    ("fl-mid-floors-y", "src", FL, "  m.y = floor_div(a.y + b.y + 1, 2);", "  m.y = floor_div(a.y + b.y, 2);", FT),
    ("fl-mid-rounds-up-x", "src", FL, "  m.x = floor_div(a.x + b.x + 1, 2);", "  m.x = floor_div(a.x + b.x + 2, 2);", FT),
    ("fl-mid-truncates-x", "src", FL, "  m.x = floor_div(a.x + b.x + 1, 2);", "  m.x = (a.x + b.x + 1) / 2;", FT),
    ("fl-mid-truncates-y", "src", FL, "  m.y = floor_div(a.y + b.y + 1, 2);", "  m.y = (a.y + b.y + 1) / 2;", FT),
    ("fl-mid-of-x", "src", FL, "  m.y = floor_div(a.y + b.y + 1, 2);", "  m.y = floor_div(a.x + b.y + 1, 2);", FT),
    ("fl-mid-of-y", "src", FL, "  m.y = floor_div(a.y + b.y + 1, 2);", "  m.y = floor_div(a.y + a.y + 1, 2);", FT),
    # The widening and the one narrowing.
    ("fl-extra-bits-7", "src", FL, "#define FLAT_EXTRA 8", "#define FLAT_EXTRA 7", FT),
    ("fl-extra-bits-9", "src", FL, "#define FLAT_EXTRA 8", "#define FLAT_EXTRA 9", FT),
    ("fl-extra-bits-3", "src", FL, "#define FLAT_EXTRA 8", "#define FLAT_EXTRA 3", FT),
    ("fl-extra-bits-0", "src", FL, "#define FLAT_EXTRA 8", "#define FLAT_EXTRA 0", FT),
    ("fl-widen-x-unscaled", "src", FL, "  w.x = (int64_t)p.x * FLAT_SCALE;", "  w.x = (int64_t)p.x;", FT),
    ("fl-widen-y-unscaled", "src", FL, "  w.y = (int64_t)p.y * FLAT_SCALE;", "  w.y = (int64_t)p.y;", FT),
    ("fl-widen-fields-32bit", "src", FL, "typedef struct WPt {\n  int64_t x;\n  int64_t y;\n} WPt;", "typedef struct WPt {\n  int32_t x;\n  int32_t y;\n} WPt;", FT),
    ("fl-narrow-floors-x", "src", FL, "  p.x = (int32_t)floor_div(w.x + FLAT_SCALE / 2, FLAT_SCALE);", "  p.x = (int32_t)floor_div(w.x, FLAT_SCALE);", FT),
    ("fl-narrow-floors-y", "src", FL, "  p.y = (int32_t)floor_div(w.y + FLAT_SCALE / 2, FLAT_SCALE);", "  p.y = (int32_t)floor_div(w.y, FLAT_SCALE);", FT),
    ("fl-narrow-ties-down-x", "src", FL, "  p.x = (int32_t)floor_div(w.x + FLAT_SCALE / 2, FLAT_SCALE);", "  p.x = (int32_t)floor_div(w.x + FLAT_SCALE / 2 - 1, FLAT_SCALE);", FT),
    ("fl-narrow-ties-down-y", "src", FL, "  p.y = (int32_t)floor_div(w.y + FLAT_SCALE / 2, FLAT_SCALE);", "  p.y = (int32_t)floor_div(w.y + FLAT_SCALE / 2 - 1, FLAT_SCALE);", FT),
    ("fl-narrow-truncates-x", "src", FL, "  p.x = (int32_t)floor_div(w.x + FLAT_SCALE / 2, FLAT_SCALE);", "  p.x = (int32_t)((w.x + FLAT_SCALE / 2) / FLAT_SCALE);", FT),
    ("fl-narrow-truncates-y", "src", FL, "  p.y = (int32_t)floor_div(w.y + FLAT_SCALE / 2, FLAT_SCALE);", "  p.y = (int32_t)((w.y + FLAT_SCALE / 2) / FLAT_SCALE);", FT),
    ("fl-narrow-ceils", "src", FL, "  p.x = (int32_t)floor_div(w.x + FLAT_SCALE / 2, FLAT_SCALE);", "  p.x = (int32_t)floor_div(w.x + FLAT_SCALE - 1, FLAT_SCALE);", FT),
    ("fl-narrow-x-of-y", "src", FL, "  p.y = (int32_t)floor_div(w.y + FLAT_SCALE / 2, FLAT_SCALE);", "  p.y = (int32_t)floor_div(w.x + FLAT_SCALE / 2, FLAT_SCALE);", FT),
    # The flatness threshold: its factor, its "- 1", the comparison.
    ("fl-threshold-t", "src", FL, "  fl->threshold = (int64_t)4 * FLAT_SCALE * ((int64_t)tolerance - 1);", "  fl->threshold = (int64_t)4 * FLAT_SCALE * (int64_t)tolerance;", FT),
    ("fl-threshold-t-minus-2", "src", FL, "  fl->threshold = (int64_t)4 * FLAT_SCALE * ((int64_t)tolerance - 1);", "  fl->threshold = (int64_t)4 * FLAT_SCALE * ((int64_t)tolerance - 2);", FT),
    ("fl-threshold-factor-8", "src", FL, "  fl->threshold = (int64_t)4 * FLAT_SCALE * ((int64_t)tolerance - 1);", "  fl->threshold = (int64_t)8 * FLAT_SCALE * ((int64_t)tolerance - 1);", FT),
    ("fl-threshold-factor-2", "src", FL, "  fl->threshold = (int64_t)4 * FLAT_SCALE * ((int64_t)tolerance - 1);", "  fl->threshold = (int64_t)2 * FLAT_SCALE * ((int64_t)tolerance - 1);", FT),
    ("fl-threshold-not-widened", "src", FL, "  fl->threshold = (int64_t)4 * FLAT_SCALE * ((int64_t)tolerance - 1);", "  fl->threshold = (int64_t)4 * ((int64_t)tolerance - 1);", FT),
    ("fl-threshold-plus-one", "src", FL, "  fl->threshold = (int64_t)4 * FLAT_SCALE * ((int64_t)tolerance - 1);", "  fl->threshold = (int64_t)4 * FLAT_SCALE * ((int64_t)tolerance - 1) + 1;", FT),
    ("fl-threshold-minus-one", "src", FL, "  fl->threshold = (int64_t)4 * FLAT_SCALE * ((int64_t)tolerance - 1);", "  fl->threshold = (int64_t)4 * FLAT_SCALE * ((int64_t)tolerance - 1) - 1;", FT),
    ("fl-quad-threshold-strict", "src", FL, "  if (sum <= fl->threshold) {\n    return emit_vertex(fl, p2, index, level);", "  if (sum < fl->threshold) {\n    return emit_vertex(fl, p2, index, level);", FT),
    ("fl-quad-x-only", "src", FL, "  int64_t sum = abs64(2 * p1.x - p0.x - p2.x) + abs64(2 * p1.y - p0.y - p2.y);", "  int64_t sum = abs64(2 * p1.x - p0.x - p2.x);", FT),
    ("fl-quad-y-only", "src", FL, "  int64_t sum = abs64(2 * p1.x - p0.x - p2.x) + abs64(2 * p1.y - p0.y - p2.y);", "  int64_t sum = abs64(2 * p1.y - p0.y - p2.y);", FT),
    ("fl-quad-sum-32bit", "src", FL, "  int64_t sum = abs64(2 * p1.x - p0.x - p2.x)", "  int32_t sum = abs64(2 * p1.x - p0.x - p2.x)", FT),
    ("fl-quad-dx-coefficient", "src", FL, "abs64(2 * p1.x - p0.x - p2.x)", "abs64(p1.x - p0.x - p2.x)", FT),
    ("fl-quad-dx-no-abs", "src", FL, "abs64(2 * p1.x - p0.x - p2.x)", "(2 * p1.x - p0.x - p2.x)", FT),
    ("fl-quad-dx-ends", "src", FL, "abs64(2 * p1.x - p0.x - p2.x)", "abs64(2 * p1.x - p0.x - p0.x)", FT),
    ("fl-quad-dy-ends", "src", FL, "abs64(2 * p1.y - p0.y - p2.y);", "abs64(2 * p1.y - p2.y - p2.y);", FT),
    ("fl-quad-dy-coefficient", "src", FL, "abs64(2 * p1.y - p0.y - p2.y);", "abs64(3 * p1.y - p0.y - p2.y);", FT),
    ("fl-quad-order", "src", FL, "  GSHGEOM_Result r = quad(fl, p0, a, m, index * 2u, level + 1);\n  return r != GSHGEOM_OK ? r : quad(fl, m, b, p2, index * 2u + 1u, level + 1);",
     "  GSHGEOM_Result r = quad(fl, m, b, p2, index * 2u + 1u, level + 1);\n  return r != GSHGEOM_OK ? r : quad(fl, p0, a, m, index * 2u, level + 1);", FT),
    ("fl-quad-error-lost", "src", FL, "  GSHGEOM_Result r = quad(fl, p0, a, m, index * 2u, level + 1);\n  return r != GSHGEOM_OK ? r : quad(fl, m, b, p2, index * 2u + 1u, level + 1);",
     "  quad(fl, p0, a, m, index * 2u, level + 1);\n  return quad(fl, m, b, p2, index * 2u + 1u, level + 1);", FT),
    ("fl-quad-halves", "src", FL, "  WPt b = mid_pt(p1, p2);\n  WPt m = mid_pt(a, b);\n  GSHGEOM_Result r = quad(", "  WPt b = mid_pt(p0, p2);\n  WPt m = mid_pt(a, b);\n  GSHGEOM_Result r = quad(", FT),
    ("fl-quad-emits-start", "src", FL, "    return emit_vertex(fl, p2, index, level);\n  }\n  WPt a = mid_pt(p0, p1);", "    return emit_vertex(fl, p0, index, level);\n  }\n  WPt a = mid_pt(p0, p1);", FT),
    ("fl-quad-emits-control", "src", FL, "    return emit_vertex(fl, p2, index, level);\n  }\n  WPt a = mid_pt(p0, p1);", "    return emit_vertex(fl, p1, index, level);\n  }\n  WPt a = mid_pt(p0, p1);", FT),
    ("fl-cubic-d1-coefficient", "src", FL, "  int64_t d1 = abs64(3 * p1.x - 2 * p0.x - p3.x) +", "  int64_t d1 = abs64(3 * p1.x - p0.x - p3.x) +", FT),
    ("fl-cubic-d1-end", "src", FL, "  int64_t d1 = abs64(3 * p1.x - 2 * p0.x - p3.x) +", "  int64_t d1 = abs64(3 * p1.x - 2 * p0.x - 2 * p3.x) +", FT),
    ("fl-cubic-d1-y", "src", FL, "      abs64(3 * p1.y - 2 * p0.y - p3.y);", "      abs64(3 * p1.y - p0.y - 2 * p3.y);", FT),
    ("fl-cubic-d1-no-y", "src", FL, "      abs64(3 * p1.y - 2 * p0.y - p3.y);", "      0;", FT),
    ("fl-cubic-d1-32bit", "src", FL, "  int64_t d1 = abs64(", "  int32_t d1 = abs64(", FT),
    ("fl-cubic-d2-x", "src", FL, "  int64_t d2 = abs64(3 * p2.x - p0.x - 2 * p3.x) +", "  int64_t d2 = abs64(3 * p2.x - 2 * p0.x - p3.x) +", FT),
    ("fl-cubic-d2-y", "src", FL, "      abs64(3 * p2.y - p0.y - 2 * p3.y);", "      abs64(3 * p2.y - 2 * p0.y - p3.y);", FT),
    ("fl-cubic-d2-no-y", "src", FL, "      abs64(3 * p2.y - p0.y - 2 * p3.y);", "      0;", FT),
    ("fl-cubic-d2-32bit", "src", FL, "  int64_t d2 = abs64(", "  int32_t d2 = abs64(", FT),
    ("fl-cubic-d2-unused", "src", FL, "  if (d1 <= fl->threshold && d2 <= fl->threshold) {", "  (void)d2;\n  if (d1 <= fl->threshold) {", FT),
    ("fl-cubic-d1-unused", "src", FL, "  if (d1 <= fl->threshold && d2 <= fl->threshold) {", "  (void)d1;\n  if (d2 <= fl->threshold) {", FT),
    ("fl-cubic-or", "src", FL, "  if (d1 <= fl->threshold && d2 <= fl->threshold) {", "  if (d1 <= fl->threshold || d2 <= fl->threshold) {", FT),
    ("fl-cubic-strict", "src", FL, "d2 <= fl->threshold) {", "d2 < fl->threshold) {", FT),
    ("fl-cubic-strict-1", "src", FL, "  if (d1 <= fl->threshold && d2", "  if (d1 < fl->threshold && d2", FT),
    ("fl-cubic-halves-d", "src", FL, "  WPt d = mid_pt(a, b);\n  WPt e = mid_pt(b, c);\n  WPt m = mid_pt(d, e);", "  WPt d = mid_pt(a, c);\n  WPt e = mid_pt(b, c);\n  WPt m = mid_pt(d, e);", FT),
    ("fl-cubic-halves-e", "src", FL, "  WPt e = mid_pt(b, c);\n  WPt m = mid_pt(d, e);", "  WPt e = mid_pt(a, c);\n  WPt m = mid_pt(d, e);", FT),
    ("fl-cubic-halves-m", "src", FL, "  WPt m = mid_pt(d, e);\n  GSHGEOM_Result r = cubic(", "  WPt m = mid_pt(a, e);\n  GSHGEOM_Result r = cubic(", FT),
    ("fl-cubic-halves-c", "src", FL, "  WPt c = mid_pt(p2, p3);", "  WPt c = mid_pt(p1, p3);", FT),
    ("fl-cubic-left", "src", FL, "  GSHGEOM_Result r = cubic(fl, p0, a, d, m, index * 2u, level + 1);", "  GSHGEOM_Result r = cubic(fl, p0, a, e, m, index * 2u, level + 1);", FT),
    ("fl-cubic-right", "src", FL, "                         : cubic(fl, m, e, c, p3, index * 2u + 1u, level + 1);", "                         : cubic(fl, m, d, c, p3, index * 2u + 1u, level + 1);", FT),
    ("fl-cubic-order", "src", FL, "  GSHGEOM_Result r = cubic(fl, p0, a, d, m, index * 2u, level + 1);\n  return r != GSHGEOM_OK ? r\n                         : cubic(fl, m, e, c, p3, index * 2u + 1u, level + 1);",
     "  GSHGEOM_Result r = cubic(fl, m, e, c, p3, index * 2u + 1u, level + 1);\n  return r != GSHGEOM_OK ? r\n                         : cubic(fl, p0, a, d, m, index * 2u, level + 1);", FT),
    ("fl-cubic-error-lost", "src", FL, "  GSHGEOM_Result r = cubic(fl, p0, a, d, m, index * 2u, level + 1);\n  return r != GSHGEOM_OK ? r\n                         : cubic(fl, m, e, c, p3, index * 2u + 1u, level + 1);",
     "  cubic(fl, p0, a, d, m, index * 2u, level + 1);\n  return cubic(fl, m, e, c, p3, index * 2u + 1u, level + 1);", FT),
    ("fl-cubic-emits-start", "src", FL, "    return emit_vertex(fl, p3, index, level);\n  }\n  WPt a = mid_pt(p0, p1);", "    return emit_vertex(fl, p0, index, level);\n  }\n  WPt a = mid_pt(p0, p1);", FT),
    ("fl-cubic-emits-control", "src", FL, "    return emit_vertex(fl, p3, index, level);\n  }\n  WPt a = mid_pt(p0, p1);", "    return emit_vertex(fl, p2, index, level);\n  }\n  WPt a = mid_pt(p0, p1);", FT),
    # The depth cap and the tolerance.
    ("fl-tolerance", "hdr", FH, "#define GSHGEOM_FLATTEN_TOLERANCE 4", "#define GSHGEOM_FLATTEN_TOLERANCE 8", FT),
    ("fl-tolerance-tight", "hdr", FH, "#define GSHGEOM_FLATTEN_TOLERANCE 4", "#define GSHGEOM_FLATTEN_TOLERANCE 2", FT),
    ("fl-depth-shallow", "hdr", FH, "#define GSHGEOM_FLATTEN_MAX_DEPTH 20", "#define GSHGEOM_FLATTEN_MAX_DEPTH 3", FT),
    ("fl-depth-8", "hdr", FH, "#define GSHGEOM_FLATTEN_MAX_DEPTH 20", "#define GSHGEOM_FLATTEN_MAX_DEPTH 8", FT),
    ("fl-depth-16", "hdr", FH, "#define GSHGEOM_FLATTEN_MAX_DEPTH 20", "#define GSHGEOM_FLATTEN_MAX_DEPTH 16", FT),
    ("fl-min-tolerance-3", "hdr", FH, "#define GSHGEOM_FLATTEN_MIN_TOLERANCE 2", "#define GSHGEOM_FLATTEN_MIN_TOLERANCE 3", FT),
    # The tolerance as a parameter: the threshold derived from it, its refusal.
    ("fl-threshold-3t", "src", FL, "  fl->threshold = (int64_t)4 * FLAT_SCALE * ((int64_t)tolerance - 1);", "  fl->threshold = (int64_t)3 * FLAT_SCALE * (int64_t)tolerance;", FT),
    ("fl-threshold-32bit", "src", FL, "  fl->threshold = (int64_t)4 * FLAT_SCALE * ((int64_t)tolerance - 1);", "  fl->threshold = (int32_t)((int64_t)4 * FLAT_SCALE * ((int64_t)tolerance - 1));", FT),
    ("fl-threshold-constant", "src", FL, "  fl->threshold = (int64_t)4 * FLAT_SCALE * ((int64_t)tolerance - 1);", "  fl->threshold = (int64_t)4 * FLAT_SCALE * (GSHGEOM_FLATTEN_TOLERANCE - 1);", FT),
    ("fl-tolerance-refuses-minimum", "src", FL, "  if (tolerance < GSHGEOM_FLATTEN_MIN_TOLERANCE || !emit) {", "  if (tolerance <= GSHGEOM_FLATTEN_MIN_TOLERANCE || !emit) {", FT),
    ("fl-tolerance-accepts-one", "src", FL, "  if (tolerance < GSHGEOM_FLATTEN_MIN_TOLERANCE || !emit) {", "  if (tolerance < GSHGEOM_FLATTEN_MIN_TOLERANCE - 1 || !emit) {", FT),
    ("fl-tolerance-unchecked", "src", FL, "  if (tolerance < GSHGEOM_FLATTEN_MIN_TOLERANCE || !emit) {", "  if (!emit) {", FT),
    ("fl-tolerance-clamped-not-refused", "src", FL, "  if (tolerance < GSHGEOM_FLATTEN_MIN_TOLERANCE || !emit) {\n    return GSHGEOM_ERR_INVALID;\n  }\n  fl->threshold = (int64_t)4 * FLAT_SCALE * ((int64_t)tolerance - 1);",
     "  if (!emit) {\n    return GSHGEOM_ERR_INVALID;\n  }\n  if (tolerance < GSHGEOM_FLATTEN_MIN_TOLERANCE) {\n    tolerance = GSHGEOM_FLATTEN_MIN_TOLERANCE;\n  }\n  fl->threshold = (int64_t)4 * FLAT_SCALE * ((int64_t)tolerance - 1);", FT),
    ("fl-emit-unchecked", "src", FL, "  if (tolerance < GSHGEOM_FLATTEN_MIN_TOLERANCE || !emit) {", "  if (tolerance < GSHGEOM_FLATTEN_MIN_TOLERANCE) {", FT),
    ("fl-cap-clamp-low", "src", FL, "  fl->levels = levels < 1 ? 1\n", "  fl->levels = levels < 0 ? 0\n", FT),
    ("fl-cap-unclamped-low", "src", FL, "  fl->levels = levels < 1 ? 1\n      : (levels > GSHGEOM_FLATTEN_MAX_DEPTH ? GSHGEOM_FLATTEN_MAX_DEPTH : levels);", "  fl->levels = levels > GSHGEOM_FLATTEN_MAX_DEPTH ? GSHGEOM_FLATTEN_MAX_DEPTH : levels;", FT),
    ("fl-default-tolerance-wrapper", "src", FL, "      path, GSHGEOM_FLATTEN_TOLERANCE, alloc, limits, out);", "      path, GSHGEOM_FLATTEN_TOLERANCE + 1, alloc, limits, out);", FT),
    ("fl-path-tolerance-ignored-quad", "src", FL, "tolerance, current, el.p[0], el.p[1], append_vertex, fp);", "GSHGEOM_FLATTEN_TOLERANCE, current, el.p[0], el.p[1], append_vertex, fp);", FT),
    ("fl-path-tolerance-ignored-cubic", "src", FL, "tolerance, current, el.p[0], el.p[1],\n          el.p[2], append_vertex, fp);", "GSHGEOM_FLATTEN_TOLERANCE, current, el.p[0], el.p[1],\n          el.p[2], append_vertex, fp);", FT),
    # The pieces' places and the count of forced ones.
    ("fl-quad-forced-not-counted", "src", FL, "  if (level >= fl->levels) {\n    fl->forced += 1;\n    return emit_vertex(fl, p2, index, level);", "  if (level >= fl->levels) {\n    return emit_vertex(fl, p2, index, level);", FT),
    ("fl-cubic-forced-not-counted", "src", FL, "  if (level >= fl->levels) {\n    fl->forced += 1;\n    return emit_vertex(fl, p3, index, level);", "  if (level >= fl->levels) {\n    return emit_vertex(fl, p3, index, level);", FT),
    ("fl-quad-forced-counted-always", "src", FL, "  if (sum <= fl->threshold) {\n    return emit_vertex(fl, p2, index, level);", "  if (sum <= fl->threshold) {\n    fl->forced += 1;\n    return emit_vertex(fl, p2, index, level);", FT),
    ("fl-cubic-forced-counted-always", "src", FL, "  if (d1 <= fl->threshold && d2 <= fl->threshold) {\n    return emit_vertex(fl, p3, index, level);", "  if (d1 <= fl->threshold && d2 <= fl->threshold) {\n    fl->forced += 1;\n    return emit_vertex(fl, p3, index, level);", FT),
    ("fl-quad-forced-not-reported", "src", FL, "  r = quad(&fl, widen(p0), widen(p1), widen(p2), 0u, 0);\n  if (forced) {\n    *forced = fl.forced;", "  r = quad(&fl, widen(p0), widen(p1), widen(p2), 0u, 0);\n  if (forced) {\n    *forced = 0;", FT),
    ("fl-cubic-forced-not-reported", "src", FL, "  r = cubic(&fl, widen(p0), widen(p1), widen(p2), widen(p3), 0u, 0);\n  if (forced) {\n    *forced = fl.forced;", "  r = cubic(&fl, widen(p0), widen(p1), widen(p2), widen(p3), 0u, 0);\n  if (forced) {\n    *forced = 0;", FT),
    ("fl-emit-level-lost", "src", FL, "  return fl->emit(fl->user, narrow(to), index, level);", "  (void)level;\n  return fl->emit(fl->user, narrow(to), index, 0);", FT),
    ("fl-emit-index-lost", "src", FL, "  return fl->emit(fl->user, narrow(to), index, level);", "  (void)index;\n  return fl->emit(fl->user, narrow(to), 0u, level);", FT),
    ("fl-emit-user-lost", "src", FL, "  return fl->emit(fl->user, narrow(to), index, level);", "  return fl->emit(NULL, narrow(to), index, level);", FT),
    ("fl-emit-result-lost", "src", FL, "  return fl->emit(fl->user, narrow(to), index, level);", "  fl->emit(fl->user, narrow(to), index, level);\n  return GSHGEOM_OK;", FT),
    ("fl-quad-left-index", "src", FL, "quad(fl, p0, a, m, index * 2u, level + 1);\n  return r != GSHGEOM_OK ? r : quad(", "quad(fl, p0, a, m, index, level + 1);\n  return r != GSHGEOM_OK ? r : quad(", FT),
    ("fl-quad-right-index", "src", FL, "quad(fl, m, b, p2, index * 2u + 1u, level + 1);", "quad(fl, m, b, p2, index * 2u, level + 1);", FT),
    ("fl-cubic-left-index", "src", FL, "cubic(fl, p0, a, d, m, index * 2u, level + 1);\n  return r != GSHGEOM_OK ? r\n", "cubic(fl, p0, a, d, m, index, level + 1);\n  return r != GSHGEOM_OK ? r\n", FT),
    ("fl-cubic-right-index", "src", FL, "cubic(fl, m, e, c, p3, index * 2u + 1u, level + 1);", "cubic(fl, m, e, c, p3, index * 2u, level + 1);", FT),
    ("fl-append-vertex-drops", "src", FL, "  return add_point(user, to);", "  (void)user;\n  (void)to;\n  return GSHGEOM_OK;", FT),
    ("fl-append-vertex-error-lost", "src", FL, "  return add_point(user, to);", "  add_point(user, to);\n  return GSHGEOM_OK;", FT),
    # The revision.
    ("fl-revision-constant", "hdr", FH, "#define GSHGEOM_FLATTEN_REVISION 2", "#define GSHGEOM_FLATTEN_REVISION 1", FT),
    ("fl-revision-constant-3", "hdr", FH, "#define GSHGEOM_FLATTEN_REVISION 2", "#define GSHGEOM_FLATTEN_REVISION 3", FT),
    ("fl-revision-function", "src", FL, "  return GSHGEOM_FLATTEN_REVISION;", "  return 1;", FT),
    ("fl-revision-in-hash", "src", FL, "  gshgeom_hash_u32(&h, FLAT_KIND);", "  gshgeom_hash_u32(&h, FLAT_KIND);\n  gshgeom_hash_u32(&h, (uint32_t)GSHGEOM_FLATTEN_REVISION);", FT),
    ("fl-revision-in-dump", "src", FL, "  if (fprintf(out, \"gshgeom.flat-path %u subpaths=%llu points=%llu hash=%016llx\\n\",\n          FLAT_FORMAT_VERSION,", "  if (fprintf(out, \"gshgeom.flat-path %u subpaths=%llu points=%llu hash=%016llx\\n\",\n          FLAT_FORMAT_VERSION + (unsigned)GSHGEOM_FLATTEN_REVISION - 1u,", FT),
    ("fl-limit-off-by-one", "src", FL, "if (fp->pt_n >= fp->max_points) {", "if (fp->pt_n > fp->max_points) {", FT),
    ("fl-limit-one-short", "src", FL, "if (fp->pt_n >= fp->max_points) {", "if (fp->pt_n + 1 >= fp->max_points) {", FT),
    ("fl-zero-limit-is-zero", "src", FL, "  fp->max_points = (limits && limits->max_points)\n      ? limits->max_points\n      : GSHGEOM_LIMITS_DEFAULT_MAX_POINTS;",
     "  fp->max_points = limits ? limits->max_points\n      : GSHGEOM_LIMITS_DEFAULT_MAX_POINTS;", FT),
    ("fl-limit-ignored", "src", FL, "  fp->max_points = (limits && limits->max_points)\n      ? limits->max_points\n      : GSHGEOM_LIMITS_DEFAULT_MAX_POINTS;",
     "  (void)limits;\n  fp->max_points = GSHGEOM_LIMITS_DEFAULT_MAX_POINTS;", FT),
    ("fl-point-count-lost", "src", FL, "  fp->sub[fp->sub_n - 1].count++;\n", "", FT),
    ("fl-subpath-first", "src", FL, "  fp->sub[fp->sub_n].first = fp->pt_n;", "  fp->sub[fp->sub_n].first = 0;", FT),
    ("fl-subpath-not-reset", "src", FL, "  fp->sub[fp->sub_n].count = 0;\n", "", FT + ["testDevice"][:0]),
    ("fl-subpath-open", "src", FL, "  fp->sub[fp->sub_n].closed = 0;\n", "  fp->sub[fp->sub_n].closed = 1;\n", FT),
    ("fl-close-lost", "src", FL, "      fp->sub[fp->sub_n - 1].closed = 1;", "      fp->sub[fp->sub_n - 1].closed = 0;", FT),
    ("fl-close-marks-first", "src", FL, "      fp->sub[fp->sub_n - 1].closed = 1;", "      fp->sub[0].closed = 1;", FT),
    ("fl-move-forgets-current", "src", FL, "      r = begin_subpath(fp, el.p[0]);\n      current = el.p[0];", "      r = begin_subpath(fp, el.p[0]);", FT),
    ("fl-line-forgets-current", "src", FL, "      r = add_point(fp, el.p[0]);\n      current = el.p[0];", "      r = add_point(fp, el.p[0]);", FT),
    ("fl-quadratic-forgets-current", "src", FL, "      current = el.p[1];\n      break;\n    case GSHGEOM_DEVICE_CUBIC:", "      break;\n    case GSHGEOM_DEVICE_CUBIC:", FT),
    ("fl-cubic-forgets-current", "src", FL, "      current = el.p[2];\n      break;\n    case GSHGEOM_DEVICE_CLOSE:", "      break;\n    case GSHGEOM_DEVICE_CLOSE:", FT),
    ("fl-quadratic-args", "src", FL, "tolerance, current, el.p[0], el.p[1], append_vertex, fp);", "tolerance, current, el.p[1], el.p[0], append_vertex, fp);", FT),
    ("fl-cubic-args", "src", FL, "tolerance, current, el.p[0], el.p[1],\n          el.p[2], append_vertex, fp);", "tolerance, current, el.p[1], el.p[0],\n          el.p[2], append_vertex, fp);", FT),
    ("fl-line-adds-start", "src", FL, "      r = add_point(fp, el.p[0]);\n      current = el.p[0];", "      r = add_point(fp, current);\n      current = el.p[0];", FT),
        ("fl-skips-first", "src", FL, "for (size_t i = 0; i < n && r == GSHGEOM_OK; i++) {\n    GSHGEOM_Device_Element el;", "for (size_t i = 1; i < n && r == GSHGEOM_OK; i++) {\n    GSHGEOM_Device_Element el;", FT),
    ("fl-begin-subpath-error", "src", FL, "  r = add_point(fp, p);\n  if (r != GSHGEOM_OK) {\n    fp->sub_n--;\n  }\n  return r;", "  r = add_point(fp, p);\n  return GSHGEOM_OK;", FT),
    ("fl-subpath-grow-error", "src", FL, "  fp->sub = grown;\n  fp->sub[fp->sub_n].first", "  if (r != GSHGEOM_OK) { fp->sub = grown; }\n  fp->sub[fp->sub_n].first", FT),
    ("fl-add-point-grow", "src", FL, "  fp->pt = grown;\n  fp->pt[fp->pt_n++] = p;", "  fp->pt = grown;\n  fp->pt[fp->pt_n] = p;\n  fp->pt_n += 1;\n  fp->pt[0].x ^= (fp->pt_n > 40) ? 1 : 0;", FT),
    ("fl-leaks-points", "src", FL, "  gcu_allocator_free(&alloc, fp->pt);\n  gcu_allocator_free(&alloc, fp->sub);\n  gcu_allocator_free(&alloc, fp);", "  gcu_allocator_free(&alloc, fp->sub);\n  gcu_allocator_free(&alloc, fp);", FT),
    ("fl-leaks-subpaths", "src", FL, "  gcu_allocator_free(&alloc, fp->pt);\n  gcu_allocator_free(&alloc, fp->sub);\n  gcu_allocator_free(&alloc, fp);", "  gcu_allocator_free(&alloc, fp->pt);\n  gcu_allocator_free(&alloc, fp);", FT),
    ("fl-leaks-path", "src", FL, "  gcu_allocator_free(&alloc, fp->pt);\n  gcu_allocator_free(&alloc, fp->sub);\n  gcu_allocator_free(&alloc, fp);", "  gcu_allocator_free(&alloc, fp->pt);\n  gcu_allocator_free(&alloc, fp->sub);", FT),
    ("fl-failure-leaks", "src", FL, "  if (r != GSHGEOM_OK) {\n    flat_destroy_parts(fp);\n    return r;\n  }\n  *out = fp;", "  if (r != GSHGEOM_OK) {\n    return r;\n  }\n  *out = fp;", FT),
    ("fl-create-oom-code", "src", FL, "  if (!fp) {\n    return GSHGEOM_ERR_OOM;\n  }", "  if (!fp) {\n    return GSHGEOM_ERR_INVALID;\n  }", FT),
    ("fl-null-out", "src", FL, "  if (!path || !out || tolerance < GSHGEOM_FLATTEN_MIN_TOLERANCE) {\n    return GSHGEOM_ERR_INVALID;\n  }\n  if (!alloc) {", "  if (!path || tolerance < GSHGEOM_FLATTEN_MIN_TOLERANCE) {\n    return GSHGEOM_ERR_INVALID;\n  }\n  if (!alloc) {", FT),
        ("fl-one-call-leaks-device", "src", FL, "  r = gshgeom_device_path_flatten(dp, alloc, limits, out);\n  gshgeom_device_path_destroy(dp);", "  r = gshgeom_device_path_flatten(dp, alloc, limits, out);", FT),
    ("fl-one-call-error-lost", "src", FL, "  GSHGEOM_Result r = gshgeom_device_path_from_model(path, map, alloc, limits, &dp);\n  if (r != GSHGEOM_OK) {\n    return r;\n  }",
     "  GSHGEOM_Result r = gshgeom_device_path_from_model(path, map, alloc, limits, &dp);\n  if (r != GSHGEOM_OK) {\n    return GSHGEOM_ERR_INTERNAL;\n  }", FT),

    # --- the flattened path's readback, hash and dump ---------------------------
    ("fh-point-count", "src", FL, "return path ? path->pt_n : 0;", "return path ? path->pt_n + 1 : 0;", FT),
    ("fh-subpath-count", "src", FL, "return path ? path->sub_n : 0;", "return path ? path->pt_n : 0;", FT),
    ("fh-points-null", "src", FL, "return (path && path->pt_n) ? path->pt : NULL;", "(void)path;\n  return NULL;", FT),
    ("fh-subpath-first", "src", FL, "  out->first = path->sub[index].first;", "  out->first = path->sub[index].count;", FT),
    ("fh-subpath-count-field", "src", FL, "  out->count = path->sub[index].count;", "  out->count = path->sub[index].first;", FT),
    ("fh-subpath-closed", "src", FL, "  out->closed = path->sub[index].closed;", "  out->closed = 0;", FT),
    ("fh-subpath-index", "src", FL, "if (!path || !out || index >= path->sub_n) {", "if (!path || !out) {", FT),
    ("fh-subpath-index-off-by-one", "src", FL, "if (!path || !out || index >= path->sub_n) {", "if (!path || !out || index > path->sub_n) {", FT),
    ("fh-version", "src", FL, "#define FLAT_FORMAT_VERSION 1u", "#define FLAT_FORMAT_VERSION 2u", FT),
    ("fh-kind", "src", FL, "#define FLAT_KIND 0x54414c46u", "#define FLAT_KIND 0x54414c47u", FT),
    ("fh-no-kind", "src", FL, "  gshgeom_hash_u32(&h, FLAT_KIND);\n", "", FT),
    ("fh-count-u32", "src", FL, "gshgeom_hash_u64(&h, (uint64_t)path->sub_n);", "gshgeom_hash_u32(&h, (uint32_t)path->sub_n);", FT),
    ("fh-no-closed", "src", FL, "    gshgeom_hash_u8(&h, s->closed);\n", "", FT),
    ("fh-closed-wide", "src", FL, "    gshgeom_hash_u8(&h, s->closed);", "    gshgeom_hash_u32(&h, s->closed);", FT),
    ("fh-no-point-count", "src", FL, "    gshgeom_hash_u64(&h, (uint64_t)s->count);\n", "", FT),
    ("fh-point-count-u32", "src", FL, "    gshgeom_hash_u64(&h, (uint64_t)s->count);", "    gshgeom_hash_u32(&h, (uint32_t)s->count);", FT),
    ("fh-y-for-x", "src", FL, "      gshgeom_hash_i32(&h, path->pt[s->first + k].x);", "      gshgeom_hash_i32(&h, path->pt[s->first + k].y);", FT),
    ("fh-no-y", "src", FL, "      gshgeom_hash_i32(&h, path->pt[s->first + k].y);\n", "", FT),
    ("fh-first-point-only", "src", FL, "    for (size_t k = 0; k < s->count; k++) {\n      gshgeom_hash_i32", "    for (size_t k = 0; k < (s->count ? 1u : 0u); k++) {\n      gshgeom_hash_i32", FT),
    ("fh-null-path", "src", FL, "  if (!path) {\n    return gshgeom_hash_value(&h);\n  }\n  gshgeom_hash_u32(&h, FLAT_FORMAT_VERSION);", "  if (!path) {\n    return 0;\n  }\n  gshgeom_hash_u32(&h, FLAT_FORMAT_VERSION);", FT),
    ("fh-dump-header", "src", FL, '"gshgeom.flat-path %u subpaths=%llu points=%llu hash=%016llx\\n"', '"gshgeom.flat_path %u subpaths=%llu points=%llu hash=%016llx\\n"', FT),
    ("fh-dump-hash-width", "src", FL, "points=%llu hash=%016llx", "points=%llu hash=%llx", FT),
    ("fh-dump-closed", "src", FL, '"S closed=%d points=%llu\\n", s->closed ? 1 : 0,', '"S closed=%d points=%llu\\n", s->closed ? 0 : 1,', FT),
    ("fh-dump-point-letter", "src", FL, '"P %08x %08x\\n"', '"p %08x %08x\\n"', FT),
    ("fh-dump-width", "src", FL, '"P %08x %08x\\n"', '"P %x %x\\n"', FT),
    ("fh-dump-swapped", "src", FL, "(unsigned)(uint32_t)path->pt[s->first + k].x,\n              (unsigned)(uint32_t)path->pt[s->first + k].y) < 0) {", "(unsigned)(uint32_t)path->pt[s->first + k].y,\n              (unsigned)(uint32_t)path->pt[s->first + k].x) < 0) {", FT),
    ("fh-dump-ignores-flush", "src", FL, "  return (fflush(out) == 0 && !ferror(out)) ? GSHGEOM_OK : GSHGEOM_ERR_IO;\n}\n", "  return GSHGEOM_OK;\n}\n", FT),
    ("fh-dump-ignores-ferror", "src", FL, "(fflush(out) == 0 && !ferror(out))", "(fflush(out) == 0)", FT),
    ("fh-dump-ignores-header-error", "src", FL, "(unsigned long long)gshgeom_flat_path_hash(path)) < 0) {\n    return GSHGEOM_ERR_IO;", "(unsigned long long)gshgeom_flat_path_hash(path)) < 0) {\n    return GSHGEOM_OK;", FT),
            ("fh-dump-null", "src", FL, "GSHGEOM_Result gshgeom_flat_path_dump(\n    const GSHGEOM_Flat_Path * path, FILE * out) {\n  if (!path || !out) {", "GSHGEOM_Result gshgeom_flat_path_dump(\n    const GSHGEOM_Flat_Path * path, FILE * out) {\n  if (!path) {", FT),

    # --- the core -------------------------------------------------------------
    ("max-points-default-zero", "src", K, "  limits->max_points = GSHGEOM_LIMITS_DEFAULT_MAX_POINTS;\n", "", PT),
    ("max-points-default-const", "hdr", KH, "#define GSHGEOM_LIMITS_DEFAULT_MAX_POINTS ((size_t)1 << 22)", "#define GSHGEOM_LIMITS_DEFAULT_MAX_POINTS ((size_t)1 << 21)", PT + FT),
]


# --- story 8: along-path parametrization --------------------------------------
A = "src/along/along.c"
AH = "include/ghoti.io/shape-geometry/along.h"
AT = ["testAlong"]

PLANTS += [
    # --- the Gauss-Legendre constants (one unit in the last place each; the
    # pinned tables and the probe see what the lengths cannot) -----------------
    ('al-gl-node-0', "src", A, '0x1.77ac94f3c7345p-3', '0x1.77ac95f3c7345p-3', AT),
    ("al-gl-node-1", "src", A, "0x1.0d129583284b4p-1", "0x1.0d129583284b5p-1", AT),
    ('al-gl-node-2', "src", A, '0x1.97e4ab249f41ep-1', '0x1.97e4ac249f41ep-1', AT),
    ('al-gl-node-3', "src", A, '0x1.ebab1cb0acc67p-1', '0x1.ebab1db0acc67p-1', AT),
    ("al-gl-weight-0", "src", A, "0x1.736360b199343p-2", "0x1.736360b199344p-2", AT),
    ("al-gl-weight-1", "src", A, "0x1.413c50a255615p-2", "0x1.413c50a255616p-2", AT),
    ("al-gl-weight-2", "src", A, "0x1.c76fb531d2b96p-3", "0x1.c76fb531d2b97p-3", AT),
    ("al-gl-weight-3", "src", A, "0x1.9ea1d04ca0374p-4", "0x1.9ea1d04ca0375p-4", AT),
    ("al-gl-weight-big", "src", A, "0x1.9ea1d04ca0374p-4", "0x1.9ea1d04ca0374p-3", AT),
    ('al-gl-pairs', "src", A, 'for (int i = 0; i < GL_HALF; i++) {\n    double off', 'for (int i = 0; i < GL_HALF - 1; i++) {\n    double off', AT),

    # --- the curve in homogeneous form, de Casteljau, differences --------------
    ("al-load-numerator-weights", "src", A, "h->c[a][i] = h->rational ? curve->w[i] * v : v;", "h->c[a][i] = v;", AT),
    ("al-load-denominator", "src", A, "h->d[i] = h->rational ? curve->w[i] : 1.0;", "h->d[i] = 1.0;", AT),
    ("al-load-dimension", "src", A, "h->dim = curve->dim;", "h->dim = 2;", AT),
    ("al-load-degree", "src", A, "h->n = gshgeom_curve_point_count(curve->kind) - 1;", "h->n = gshgeom_curve_point_count(curve->kind) > 3 ? 3 : 2;", AT),
    ("al-load-rational-flag", "src", A, "h->rational = gshgeom_curve_is_rational(curve->kind);", "h->rational = 0;", AT),
    ('al-mag-no-z', "src", A, 'return dim == 3 ? gmath_hypot(h, v[2]) : h;', 'return (void)dim, h;', AT),
    ("al-bez-s", "src", A, "  double s = 1.0 - t;\n  for (int level = m;", "  double s = t;\n  for (int level = m;", AT),
    ("al-bez-levels", "src", A, "for (int level = m; level > 0; level--) {", "for (int level = m; level > 1; level--) {", AT),
    ("al-bez-lerp", "src", A, "h[i] = s * h[i] + t * h[i + 1];", "h[i] = s * h[i + 1] + t * h[i];", AT),
    ("al-deriv-difference", "src", A, "h[i] = h[i + 1] - h[i];", "h[i] = h[i + 1] + h[i];", AT),
    ("al-deriv-factor", "src", A, "factor *= (double)(m - j);", "factor *= (double)(m - j - 1 + 1 - (j == 0));", AT),
    ("al-deriv-degree", "src", A, "return factor * bez(h, m - k, t);", "return factor * bez(h, m - k + (k > 1), t);", AT),
    ("al-e-polynomial-order", "src", A, "e[a] = deriv(h->c[a], h->n, k + 1, t);", "e[a] = deriv(h->c[a], h->n, k, t);", AT),
    ("al-e-rational-sign", "src", A, "deriv(h->c[a], h->n, j + 1, t) * deriv(h->d, h->n, k - j, t) -", "deriv(h->c[a], h->n, j + 1, t) * deriv(h->d, h->n, k - j, t) +", AT),
    ("al-e-rational-second", "src", A, "deriv(h->c[a], h->n, j, t) * deriv(h->d, h->n, k - j + 1, t);", "deriv(h->c[a], h->n, j, t) * deriv(h->d, h->n, k - j, t);", AT),
        ("al-speed-denominator", "src", A, "return m / (d * d);", "return m / d;", AT),
        ("al-gl-half-width", "src", A, "double half = (b - a) * 0.5;\n  double sum = 0.0;", "double half = (b - a);\n  double sum = 0.0;", AT),
    ("al-gl-centre", "src", A, "double c = (a + b) * 0.5;\n  double half", "double c = a;\n  double half", AT),
    ("al-gl-one-side", "src", A, "sum += gl_weight[i] * (speed(h, c - off) + speed(h, c + off));", "sum += gl_weight[i] * (speed(h, c + off) + speed(h, c + off));", AT),
    ("al-gl-no-scale", "src", A, "return half * sum;", "return sum;", AT),
    ("al-poly-dropped-edge", "src", A, "  int n = gshgeom_curve_point_count(curve->kind) - 1;\n  double total = 0.0;", "  int n = gshgeom_curve_point_count(curve->kind) - 2;\n  double total = 0.0;", AT),
    ("al-poly-last-edge", "src", A, "total += mag(curve->dim, d);", "total = mag(curve->dim, d);", AT),
    ("al-poly-no-z", "src", A, "curve->p[i + 1].y - curve->p[i].y, curve->p[i + 1].z - curve->p[i].z};", "curve->p[i + 1].y - curve->p[i].y, 0.0};", AT),

    # --- the table ---------------------------------------------------------------
    ("al-emit-limit-plus-one", "src", A, "if (tab->n >= tab->max) {", "if (tab->n > tab->max) {", AT),
    ("al-emit-limit-minus-one", "src", A, "if (tab->n >= tab->max) {", "if (tab->n + 1 >= tab->max) {", AT),
    ("al-emit-always-stores", "src", A, "  if (tab->store) {\n    void * grown = tab->row;", "  if (1) {\n    void * grown = tab->row;", AT),
    ('al-emit-grow-error', "src", A, '    if (r != GSHGEOM_OK) {\n      return r;\n    }\n    tab->row = grown;', '    (void)r;\n    tab->row = grown;', AT),
    ("al-emit-segment", "src", A, "row->seg = tab->cur;", "row->seg = 0;", AT),
    ("al-emit-parameters-swapped", "src", A, "row->t0 = t0;\n    row->t1 = t1;", "row->t0 = t1;\n    row->t1 = t0;", AT),
    ("al-emit-start", "src", A, "row->s0 = tab->sum;", "row->s0 = 0.0;", AT),
    ("al-emit-length", "src", A, "row->len = len;", "row->len = len * 0.5;", AT),
    ("al-emit-count", "src", A, "  tab->n++;\n  tab->sum += len;", "  tab->sum += len;", AT),
    ("al-emit-sum", "src", A, "  tab->n++;\n  tab->sum += len;", "  tab->n++;\n  tab->sum = len;", AT),
    ("al-refine-depth", "src", A, "if (depth >= GSHGEOM_MEASURE_MAX_DEPTH || diff <= scaled * (b - a)) {", "if (depth >= 3 || diff <= scaled * (b - a)) {", AT),
    ("al-refine-depth-const", "hdr", AH, "#define GSHGEOM_MEASURE_MAX_DEPTH 52", "#define GSHGEOM_MEASURE_MAX_DEPTH 4", AT),
    ("al-refine-width", "src", A, "diff <= scaled * (b - a)) {", "diff <= scaled) {", AT),
    ("al-refine-tighter", "src", A, "diff <= scaled * (b - a)) {", "diff <= scaled * (b - a) * 0.5) {", AT),
    ("al-refine-looser", "src", A, "diff <= scaled * (b - a)) {", "diff <= scaled * (b - a) * 2.0) {", AT),
    ("al-refine-diff-left-only", "src", A, "double diff = dabs(whole - (left + right));", "double diff = dabs(whole - left);", AT),
    ("al-refine-diff-sign", "src", A, "double diff = dabs(whole - (left + right));", "double diff = whole - (left + right);", AT),
    ("al-refine-midpoint", "src", A, "double m = (a + b) * 0.5;\n  double left = gl(h, a, m);", "double m = a + (b - a) * 0.25;\n  double left = gl(h, a, m);", AT),
    ("al-refine-right-half", "src", A, "double right = gl(h, m, b);\n  double diff", "double right = gl(h, m, b) * 1.0000000000001;\n  double diff", AT),
    ("al-refine-order", "src", A, "  GSHGEOM_Result r = refine(h, tab, a, m, left, scaled, depth + 1);\n  return r != GSHGEOM_OK ? r : refine(h, tab, m, b, right, scaled, depth + 1);",
     "  GSHGEOM_Result r = refine(h, tab, m, b, right, scaled, depth + 1);\n  return r != GSHGEOM_OK ? r : refine(h, tab, a, m, left, scaled, depth + 1);", AT),
    ('al-refine-left-error', "src", A, 'return r != GSHGEOM_OK ? r : refine(h, tab, m, b, right, scaled, depth + 1);', '(void)r;\n  return refine(h, tab, m, b, right, scaled, depth + 1);', AT),
    ('al-refine-emit-error', "src", A, 'return r != GSHGEOM_OK ? r : emit(tab, m, b, right);', '(void)r;\n    return emit(tab, m, b, right);', AT),
    ("al-refine-whole-passed", "src", A, "GSHGEOM_Result r = refine(h, tab, a, m, left, scaled, depth + 1);", "GSHGEOM_Result r = refine(h, tab, a, m, left + right, scaled, depth + 1);", AT),
    ("al-measure-line-closed-form", "src", A, "if (curve->kind == GSHGEOM_CURVE_LINE) {\n    return emit", "if (0) {\n    return emit", AT),
    ("al-measure-line-parameters", "src", A, "return emit(tab, 0.0, 1.0, poly);", "return emit(tab, 0.0, 0.5, poly);", AT),
    ("al-measure-whole", "src", A, "double whole = gl(&h, 0.0, 1.0);", "double whole = gl(&h, 0.0, 0.5);", AT),
    ("al-measure-scale", "src", A, "tolerance * poly, 0);", "tolerance, 0);", AT),
    ("al-measure-start-depth", "src", A, "tolerance * poly, 0);", "tolerance * poly, GSHGEOM_MEASURE_MAX_DEPTH - 2);", AT),

    # --- the tolerance argument --------------------------------------------------------
    ("al-tol-nan", "src", A, "if (!is_finite(tolerance) || tolerance < 0.0) {", "if (tolerance < 0.0) {", AT),
    ("al-tol-negative", "src", A, "if (!is_finite(tolerance) || tolerance < 0.0) {", "if (!is_finite(tolerance) || tolerance < -1.0) {", AT),
    ("al-tol-negative-limit", "src", A, "if (!is_finite(tolerance) || tolerance < 0.0) {\n    return GSHGEOM_ERR_INVALID;", "if (!is_finite(tolerance) || tolerance < 0.0) {\n    return GSHGEOM_ERR_LIMIT;", AT),
    ("al-tol-zero", "src", A, "if (tolerance == 0.0) {\n    *out = GSHGEOM_MEASURE_TOLERANCE_DEFAULT;", "if (tolerance == 1.0) {\n    *out = GSHGEOM_MEASURE_TOLERANCE_DEFAULT;", AT),
    ("al-tol-min-exclusive", "src", A, "if (tolerance < GSHGEOM_MEASURE_TOLERANCE_MIN ||", "if (tolerance <= GSHGEOM_MEASURE_TOLERANCE_MIN ||", AT),
    ("al-tol-max-exclusive", "src", A, "tolerance > GSHGEOM_MEASURE_TOLERANCE_MAX) {", "tolerance >= GSHGEOM_MEASURE_TOLERANCE_MAX) {", AT),
    ("al-tol-min-unchecked", "src", A, "if (tolerance < GSHGEOM_MEASURE_TOLERANCE_MIN ||\n      tolerance > GSHGEOM_MEASURE_TOLERANCE_MAX) {", "if (tolerance > GSHGEOM_MEASURE_TOLERANCE_MAX) {", AT),
    ("al-tol-max-unchecked", "src", A, "if (tolerance < GSHGEOM_MEASURE_TOLERANCE_MIN ||\n      tolerance > GSHGEOM_MEASURE_TOLERANCE_MAX) {", "if (tolerance < GSHGEOM_MEASURE_TOLERANCE_MIN) {", AT),
    ("al-tol-limit-code", "src", A, "tolerance > GSHGEOM_MEASURE_TOLERANCE_MAX) {\n    return GSHGEOM_ERR_LIMIT;", "tolerance > GSHGEOM_MEASURE_TOLERANCE_MAX) {\n    return GSHGEOM_ERR_INVALID;", AT),
    ("al-tol-min-const", "hdr", AH, "#define GSHGEOM_MEASURE_TOLERANCE_MIN 0x1p-46", "#define GSHGEOM_MEASURE_TOLERANCE_MIN 0x1p-45", AT),
    ("al-tol-max-const", "hdr", AH, "#define GSHGEOM_MEASURE_TOLERANCE_MAX 0x1p-4", "#define GSHGEOM_MEASURE_TOLERANCE_MAX 0x1p-5", AT),
    ("al-tol-default-const", "hdr", AH, "#define GSHGEOM_MEASURE_TOLERANCE_DEFAULT 0x1p-40", "#define GSHGEOM_MEASURE_TOLERANCE_DEFAULT 0x1p-30", AT),

    # --- the single-curve length ---------------------------------------------------
    ("al-cl-null-out", "src", A, "  if (!curve || !out) {\n    return GSHGEOM_ERR_INVALID;\n  }\n  GSHGEOM_Result r = gshgeom_curve_check(curve);", "  if (!curve) {\n    return GSHGEOM_ERR_INVALID;\n  }\n  GSHGEOM_Result r = gshgeom_curve_check(curve);", AT),
    ("al-cl-unchecked-curve", "src", A, "  GSHGEOM_Result r = gshgeom_curve_check(curve);\n  if (r != GSHGEOM_OK) {\n    return r;\n  }\n  double tol = 0.0;", "  GSHGEOM_Result r = GSHGEOM_OK;\n  double tol = 0.0;", AT),
    ("al-cl-unchecked-tolerance", "src", A, "  r = check_tolerance(tolerance, &tol);\n  if (r != GSHGEOM_OK) {\n    return r;\n  }\n  double poly = polygon_length(curve);", "  tol = tolerance == 0.0 ? GSHGEOM_MEASURE_TOLERANCE_DEFAULT : tolerance;\n  double poly = polygon_length(curve);", AT),
    ("al-cl-max-rows", "src", A, "tab.max = GSHGEOM_LIMITS_DEFAULT_MAX_POINTS;\n  r = measure", "tab.max = 3;\n  r = measure", AT),
    ("al-cl-result", "src", A, "*out = tab.sum;\n  return GSHGEOM_OK;\n}\n\n/* ---", "*out = tab.sum * 1.000001;\n  return GSHGEOM_OK;\n}\n\n/* ---", AT),
    
    # --- building a measured path ----------------------------------------------------
    ("al-create-null-out", "src", A, "  if (!path || !out) {\n    return GSHGEOM_ERR_INVALID;\n  }\n  double tol = 0.0;", "  if (!path) {\n    return GSHGEOM_ERR_INVALID;\n  }\n  double tol = 0.0;", AT),
    ("al-create-null-path", "src", A, "  if (!path || !out) {\n    return GSHGEOM_ERR_INVALID;\n  }\n  double tol = 0.0;", "  if (!out) {\n    return GSHGEOM_ERR_INVALID;\n  }\n  double tol = 0.0;", AT),
    ("al-create-tolerance-unchecked", "src", A, "  GSHGEOM_Result r = check_tolerance(tolerance, &tol);\n  if (r != GSHGEOM_OK) {\n    return r;\n  }\n  if (!alloc) {", "  GSHGEOM_Result r = GSHGEOM_OK;\n  tol = tolerance == 0.0 ? GSHGEOM_MEASURE_TOLERANCE_DEFAULT : tolerance;\n  if (!alloc) {", AT),
    ("al-create-oom-code", "src", A, "  if (!m) {\n    return GSHGEOM_ERR_OOM;\n  }\n  m->alloc = *alloc;", "  if (!m) {\n    return GSHGEOM_ERR_INVALID;\n  }\n  m->alloc = *alloc;", AT),
    ("al-create-dimension", "src", A, "m->dim = gshgeom_path_dimension(path);", "m->dim = 2;", AT),
    ("al-create-zero-limit", "src", A, "m->tab.max = (limits && limits->max_points) ? limits->max_points\n                                               : GSHGEOM_LIMITS_DEFAULT_MAX_POINTS;", "m->tab.max = limits ? limits->max_points : GSHGEOM_LIMITS_DEFAULT_MAX_POINTS;", AT),
    ('al-create-limit-ignored', "src", A, 'm->tab.max = (limits && limits->max_points) ? limits->max_points\n                                               : GSHGEOM_LIMITS_DEFAULT_MAX_POINTS;', '(void)limits;\n  m->tab.max = GSHGEOM_LIMITS_DEFAULT_MAX_POINTS;', AT),
    ("al-create-origin-first", "src", A, "      if (!m->has_origin) {\n        m->origin = el.curve.p[0];\n        m->has_origin = 1;\n      }", "      m->origin = el.curve.p[0];\n      m->has_origin = 1;", AT),
    ("al-create-move-no-subpath", "src", A, "      need_sub = 1;\n      continue;", "      continue;", AT),
    ("al-create-close-no-subpath", "src", A, "      m->sub[m->sub_n - 1].closed = 1;\n      need_sub = 1;", "      m->sub[m->sub_n - 1].closed = 1;", AT),
    ("al-create-close-flag", "src", A, "      m->sub[m->sub_n - 1].closed = 1;\n      need_sub = 1;", "      m->sub[m->sub_n - 1].closed = 0;\n      need_sub = 1;", AT),
    ("al-create-subpath-end", "src", A, "    m->sub[m->sub_n - 1].end = m->tab.sum;\n    if (el.verb", "    if (el.verb", AT),
    ("al-create-subpath-open-error", "src", A, "      r = open_subpath(m);\n      need_sub = 0;\n      if (r != GSHGEOM_OK) {\n        break;\n      }", "      r = open_subpath(m);\n      need_sub = 0;", AT),
            ("al-create-leaks-on-failure", "src", A, "  if (r != GSHGEOM_OK) {\n    gshgeom_measured_destroy(m);\n    return r;\n  }\n  *out = m;", "  if (r != GSHGEOM_OK) {\n    return r;\n  }\n  *out = m;", AT),
    ("al-subpath-begin", "src", A, "m->sub[m->sub_n].begin = m->tab.sum;", "m->sub[m->sub_n].begin = 0.0;", AT),
        ("al-subpath-initial-closed", "src", A, "m->sub[m->sub_n].closed = 0;\n  m->sub_n++;", "m->sub[m->sub_n].closed = 1;\n  m->sub_n++;", AT),
    ('al-subpath-grow-error', "src", A, '  if (r != GSHGEOM_OK) {\n    return r;\n  }\n  m->sub = grown;', '  (void)r;\n  m->sub = grown;', AT),
    ("al-segment-extent", "src", A, "  if (poly == 0.0) {\n    return GSHGEOM_OK; /* no extent: nothing selects it */\n  }", "  if (poly < 0.0) {\n    return GSHGEOM_OK;\n  }", AT),
    ('al-segment-grow-error', "src", A, '  if (r != GSHGEOM_OK) {\n    return r;\n  }\n  m->seg = grown;', '  (void)r;\n  m->seg = grown;', AT),
    ('al-segment-element', "src", A, 'm->seg[m->seg_n].element = element;', '(void)element;\n  m->seg[m->seg_n].element = m->seg_n;', AT),
    ("al-segment-curve", "src", A, "m->seg[m->seg_n].curve = *curve;", "m->seg[m->seg_n].curve = m->seg_n ? m->seg[0].curve : *curve;", AT),
    ("al-segment-current", "src", A, "m->tab.cur = m->seg_n;\n  m->seg_n++;", "m->tab.cur = 0;\n  m->seg_n++;", AT),
    ("al-destroy-leaks-segments", "src", A, "  gcu_allocator_free(&alloc, m->seg);\n", "", AT),
    ("al-destroy-leaks-rows", "src", A, "  gcu_allocator_free(&alloc, m->tab.row);\n", "", AT),
    ("al-destroy-leaks-subpaths", "src", A, "  gcu_allocator_free(&alloc, m->sub);\n", "", AT),
    ("al-destroy-leaks-object", "src", A, "  gcu_allocator_free(&alloc, m->sub);\n  gcu_allocator_free(&alloc, m);", "  gcu_allocator_free(&alloc, m->sub);", AT),

    # --- accessors ----------------------------------------------------------------------
    ("al-dimension", "src", A, "return m ? m->dim : 0;", "return m ? 2 : 0;", AT),
    ("al-length", "src", A, "return m ? m->tab.sum : 0.0;", "return m ? m->tab.sum + 1.0 : 0.0;", AT),
    ("al-segment-count", "src", A, "return m ? m->seg_n : 0;", "return m ? m->tab.n : 0;", AT),
    ("al-row-count", "src", A, "return m ? m->tab.n : 0;", "return m ? m->seg_n : 0;", AT),
    ("al-subpath-count", "src", A, "return m ? m->sub_n : 0;", "return m ? m->seg_n : 0;", AT),
    ("al-subpath-index", "src", A, "if (!m || !out || index >= m->sub_n) {", "if (!m || !out || index > m->sub_n) {", AT),
    ("al-subpath-null-out", "src", A, "if (!m || !out || index >= m->sub_n) {", "if (!m || index >= m->sub_n) {", AT),
    ("al-subpath-value", "src", A, "*out = m->sub[index];", "*out = m->sub[0];", AT),

    # --- queries -----------------------------------------------------------------------------
    ("al-query-nan", "src", A, "if (!m || !out || !is_finite(s) || !m->has_origin) {", "if (!m || !out || !m->has_origin) {", AT),
    ("al-query-empty", "src", A, "if (!m || !out || !is_finite(s) || !m->has_origin) {", "if (!m || !out || !is_finite(s)) {", AT),
    ("al-query-null-out", "src", A, "if (!m || !out || !is_finite(s) || !m->has_origin) {", "if (!m || !is_finite(s) || !m->has_origin) {", AT),
    ("al-query-end-excluded", "src", A, "if (s < 0.0 || s > length) {", "if (s < 0.0 || s >= length) {", AT),
    ("al-query-start-excluded", "src", A, "if (s < 0.0 || s > length) {", "if (s <= 0.0 || s > length) {", AT),
    ("al-query-high-unchecked", "src", A, "if (s < 0.0 || s > length) {", "if (s < 0.0) {", AT),
    ("al-query-low-unchecked", "src", A, "if (s < 0.0 || s > length) {", "if (s > length) {", AT),
    ("al-query-clamp-inverted", "src", A, "    if (!clamp) {\n      return GSHGEOM_ERR_INVALID;", "    if (clamp) {\n      return GSHGEOM_ERR_INVALID;", AT),
    ('al-query-clamp-never', "src", A, '    if (!clamp) {\n      return GSHGEOM_ERR_INVALID;', '    (void)clamp;\n    if (1) {\n      return GSHGEOM_ERR_INVALID;', AT),
    ("al-query-clamp-swapped", "src", A, "s = s < 0.0 ? 0.0 : length;", "s = s < 0.0 ? length : 0.0;", AT),
    ("al-query-no-rows-point", "src", A, "    f.point = m->origin;\n    *out = f;", "    *out = f;", AT),
    ("al-query-no-rows-test", "src", A, "  if (m->tab.n == 0) {\n    f.point", "  if (m->tab.n == 1) {\n    f.point", AT),
    ("al-query-boundary-earlier", "src", A, "if (m->tab.row[mid].s0 <= s) {", "if (m->tab.row[mid].s0 < s) {", AT),
    ("al-query-search-wrong-way", "src", A, "if (m->tab.row[mid].s0 <= s) {\n      lo = mid + 1;\n    } else {\n      hi = mid;", "if (m->tab.row[mid].s0 <= s) {\n      hi = mid;\n    } else {\n      lo = mid + 1;", AT),
    ("al-query-target", "src", A, "double target = s - row->s0;", "double target = row->s0 - s;", AT),
    ("al-query-start-parameter", "src", A, "  if (target <= 0.0) {\n    u = row->t0;", "  if (target <= 0.0) {\n    u = row->t1;", AT),
    ("al-query-end-test", "src", A, "} else if (s >= row->s0 + row->len) {", "} else if (s > row->s0 + row->len) {", AT),
    ("al-query-end-uses-length", "src", A, "} else if (s >= row->s0 + row->len) {", "} else if (target >= row->len) {", AT),
    ("al-query-end-parameter", "src", A, "    u = row->t1;\n  } else if (seg->curve.kind == GSHGEOM_CURVE_LINE) {", "    u = row->t0;\n  } else if (seg->curve.kind == GSHGEOM_CURVE_LINE) {", AT),
        ("al-query-line-parameter", "src", A, "    u = target / row->len;", "    u = row->t0 + target / row->len * 0.999;", AT),
    ('al-query-tangent-skipped', "src", A, '  if (tangent(&h, u, tan)) {\n    f.has_tangent = 1;', '  if (tangent(&h, u, tan) && 0) {\n    f.has_tangent = 1;', AT),
    ("al-query-tangent-flag", "src", A, "    f.has_tangent = 1;", "    f.has_tangent = 0;", AT),
    ("al-query-tangent-z", "src", A, "f.tangent.z = h.dim == 3 ? tan[2] : 0.0;", "f.tangent.z = 0.0;", AT),
    ("al-query-normal-3d", "src", A, "    if (h.dim == 2) {\n      f.normal = gshgeom_normal_left(f.tangent);", "    if (h.dim == 3) {\n      f.normal = gshgeom_normal_left(f.tangent);", AT),
    ("al-query-normal-right", "src", A, "f.normal = gshgeom_normal_left(f.tangent);", "f.normal = gshgeom_normal_right(f.tangent);", AT),
    ("al-query-normal-flag", "src", A, "      f.has_normal = 1;", "      f.has_normal = 0;", AT),
    ("al-query-element", "src", A, "  f.element = seg->element;", "  f.element = row->seg;", AT),
    ("al-query-parameter", "src", A, "  f.t = u;\n  *out = f;", "  f.t = row->t0;\n  *out = f;", AT),
        
    # --- the inverse of the arc length -----------------------------------------------------------
        ("al-solve-bracket", "src", A, "    if (f > 0.0) {\n      hi = u;\n    } else {\n      lo = u;\n    }", "    if (f > 0.0) {\n      lo = u;\n    } else {\n      hi = u;\n    }", AT),
    ("al-solve-stop", "src", A, "double stop = row->len * 0x1p-52;", "double stop = row->len * 0x1p-10;", AT),
    ("al-solve-target", "src", A, "double f = gl(h, row->t0, u) - target;", "double f = gl(h, row->t0, u) - target * 1.0001;", AT),
    ("al-solve-from-start", "src", A, "double f = gl(h, row->t0, u) - target;", "double f = gl(h, 0.0, u) - target;", AT),
    ("al-solve-one-step", "src", A, "for (int iteration = 0; iteration < 64; iteration++) {", "for (int iteration = 0; iteration < 1; iteration++) {", AT),
    ("al-solve-never-bisects", "src", A, "    if (!(next > lo && next < hi)) {\n      next = lo + (hi - lo) * 0.5;\n    }\n    if (next == u) {", "    if (next < lo || next > hi) {\n      next = lo;\n    }\n    if (next == u) {", AT),
    ("al-solve-no-bracket-update", "src", A, "    if (f > 0.0) {\n      hi = u;", "    if (f > 1.0e300) {\n      hi = u;", AT),
    
    # --- the tangent ----------------------------------------------------------------------------------
    ("al-tangent-threshold-large", "src", A, "double thr = e_scale(h) * 0x1p-48;", "double thr = e_scale(h) * 0x1p-1;", AT),
        ("al-tangent-orders", "src", A, "for (int k = 0; k <= 4; k++) {", "for (int k = 0; k <= 1; k++) {", AT),
    ("al-tangent-no-fallback", "src", A, "for (int k = 0; k <= 4; k++) {", "for (int k = 0; k <= 0; k++) {", AT),
    ("al-tangent-first-nonzero", "src", A, "    if (m > thr) {\n      /* The right-hand", "    if (m > thr * 0.0 + 1e300) {\n      /* The right-hand", AT),
    ("al-tangent-end-sign", "src", A, "double sign = (u >= 1.0 && (k & 1)) ? -1.0 : 1.0;", "double sign = 1.0;", AT),
    ("al-tangent-end-test", "src", A, "double sign = (u >= 1.0 && (k & 1)) ? -1.0 : 1.0;", "double sign = (u > 1.0 && (k & 1)) ? -1.0 : 1.0;", AT),
    ("al-tangent-end-parity", "src", A, "double sign = (u >= 1.0 && (k & 1)) ? -1.0 : 1.0;", "double sign = (u >= 1.0 && (k & 2)) ? -1.0 : 1.0;", AT),
    ("al-tangent-end-always-flipped", "src", A, "double sign = (u >= 1.0 && (k & 1)) ? -1.0 : 1.0;", "double sign = (u >= 1.0 || (k & 1)) ? -1.0 : 1.0;", AT),
    ("al-tangent-normalised", "src", A, "out[a] = (sign * e[a]) / m;", "out[a] = (sign * e[a]) / (m * 1.0000001);", AT),
    ("al-tangent-undefined-flag", "src", A, "      return 1;\n    }\n  }\n  return 0;\n}", "      return 1;\n    }\n  }\n  return 1;\n}", AT),
            
    # --- hash ------------------------------------------------------------------------------------------
    ("alh-version", "src", A, "#define MEASURED_FORMAT_VERSION 1u", "#define MEASURED_FORMAT_VERSION 2u", AT),
    ("alh-kind", "src", A, "#define MEASURED_KIND 0x5341454du", "#define MEASURED_KIND 0x5341454eu", AT),
    ("alh-no-kind", "src", A, "  gshgeom_hash_u32(&h, MEASURED_KIND);\n", "", AT),
    ("alh-no-dimension", "src", A, "  gshgeom_hash_u32(&h, (uint32_t)m->dim);\n", "", AT),
    ("alh-segments-u32", "src", A, "gshgeom_hash_u64(&h, (uint64_t)m->seg_n);", "gshgeom_hash_u32(&h, (uint32_t)m->seg_n);", AT),
    ("alh-rows-u32", "src", A, "gshgeom_hash_u64(&h, (uint64_t)m->tab.n);", "gshgeom_hash_u32(&h, (uint32_t)m->tab.n);", AT),
    ("alh-subpaths-u32", "src", A, "gshgeom_hash_u64(&h, (uint64_t)m->sub_n);", "gshgeom_hash_u32(&h, (uint32_t)m->sub_n);", AT),
    ("alh-no-origin", "src", A, "  hash_point(&h, &m->origin, m->dim);\n", "", AT),
    ("alh-no-length", "src", A, "  gshgeom_hash_f64(&h, m->tab.sum);\n  for (size_t i = 0; i < m->seg_n", "  for (size_t i = 0; i < m->seg_n", AT),
    ("alh-no-z", "src", A, "  if (dim == 3) {\n    gshgeom_hash_f64(h, p->z);", "  if (dim == 4) {\n    gshgeom_hash_f64(h, p->z);", AT),
    ("alh-no-element", "src", A, "    gshgeom_hash_u64(&h, (uint64_t)s->element);\n", "", AT),
    ("alh-no-curve-kind", "src", A, "    gshgeom_hash_u8(&h, (uint8_t)s->curve.kind);\n", "", AT),
    ("alh-kind-wide", "src", A, "gshgeom_hash_u8(&h, (uint8_t)s->curve.kind);", "gshgeom_hash_u32(&h, (uint32_t)s->curve.kind);", AT),
    ("alh-first-point-only", "src", A, "    for (int k = 0; k < n; k++) {\n      hash_point(&h, &s->curve.p[k], m->dim);", "    for (int k = 0; k < 1; k++) {\n      hash_point(&h, &s->curve.p[k], m->dim);", AT),
    ("alh-no-weights", "src", A, "    if (gshgeom_curve_is_rational(s->curve.kind)) {\n      for (int k = 0; k < n; k++) {\n        gshgeom_hash_f64(&h, s->curve.w[k]);", "    if (gshgeom_curve_is_rational(s->curve.kind)) {\n      for (int k = 0; k < 0; k++) {\n        gshgeom_hash_f64(&h, s->curve.w[k]);", AT),
    ("alh-no-row-segment", "src", A, "    gshgeom_hash_u64(&h, (uint64_t)r->seg);\n", "", AT),
    ("alh-no-row-t0", "src", A, "    gshgeom_hash_f64(&h, r->t0);\n", "", AT),
    ("alh-no-row-t1", "src", A, "    gshgeom_hash_f64(&h, r->t1);\n", "", AT),
    ("alh-no-row-start", "src", A, "    gshgeom_hash_f64(&h, r->s0);\n", "", AT),
    ("alh-no-row-length", "src", A, "    gshgeom_hash_f64(&h, r->len);\n", "", AT),
    ("alh-no-subpath-begin", "src", A, "    gshgeom_hash_f64(&h, m->sub[i].begin);\n", "", AT),
    ("alh-no-subpath-end", "src", A, "    gshgeom_hash_f64(&h, m->sub[i].end);\n", "", AT),
    ("alh-no-subpath-closed", "src", A, "    gshgeom_hash_u8(&h, m->sub[i].closed ? 1 : 0);\n", "", AT),
    ("alh-closed-inverted", "src", A, "gshgeom_hash_u8(&h, m->sub[i].closed ? 1 : 0);", "gshgeom_hash_u8(&h, m->sub[i].closed ? 0 : 1);", AT),
    ("alh-null", "src", A, "  if (!m) {\n    return gshgeom_hash_value(&h);\n  }\n  gshgeom_hash_u32(&h, MEASURED_FORMAT_VERSION);", "  if (!m) {\n    return 0;\n  }\n  gshgeom_hash_u32(&h, MEASURED_FORMAT_VERSION);", AT),

    # --- dump -------------------------------------------------------------------------------------------
    ("ald-header", "src", A, '"gshgeom.measured %u dim=%d', '"gshgeom_measured %u dim=%d', AT),
    ("ald-length-width", "src", A, '"length=%016llx hash=%016llx\\n"', '"length=%llx hash=%016llx\\n"', AT),
    ("ald-hash-width", "src", A, '"length=%016llx hash=%016llx\\n"', '"length=%016llx hash=%llx\\n"', AT),
    ("ald-number-width", "src", A, 'return fprintf(out, " %016llx", (unsigned long long)gshgeom_f64_bits(v));', 'return fprintf(out, " %llx", (unsigned long long)gshgeom_f64_bits(v));', AT),
    ("ald-no-z", "src", A, "if (n >= 0 && dim == 3) n = put_f64(out, p->z);", "if (n >= 0 && dim == 4) n = put_f64(out, p->z);", AT),
    ("ald-origin-letter", "src", A, "fputc('O', out)", "fputc('o', out)", AT),
    ("ald-segment-letter", "src", A, '"G %llu %d"', '"g %llu %d"', AT),
    ("ald-segment-kind", "src", A, "(int)s->curve.kind) < 0) {", "(int)s->curve.kind + 1) < 0) {", AT),
    ("ald-no-weights", "src", A, "    if (gshgeom_curve_is_rational(s->curve.kind)) {\n      for (int k = 0; k < n; k++) {\n        if (put_f64(out, s->curve.w[k]) < 0) {", "    if (gshgeom_curve_is_rational(s->curve.kind)) {\n      for (int k = 0; k < 0; k++) {\n        if (put_f64(out, s->curve.w[k]) < 0) {", AT),
    ("ald-row-letter", "src", A, '"R %llu"', '"r %llu"', AT),
    ("ald-row-swapped", "src", A, "put_f64(out, r->t0) < 0 || put_f64(out, r->t1) < 0 ||", "put_f64(out, r->t1) < 0 || put_f64(out, r->t0) < 0 ||", AT),
    ("ald-subpath-letter", "src", A, "fputc('S', out)", "fputc('s', out)", AT),
    ("ald-subpath-closed", "src", A, 'fprintf(out, " %d\\n", m->sub[i].closed ? 1 : 0)', 'fprintf(out, " %d\\n", m->sub[i].closed ? 0 : 1)', AT),
    ("ald-ignores-flush", "src", A, "  return (fflush(out) == 0 && !ferror(out)) ? GSHGEOM_OK : GSHGEOM_ERR_IO;\n}\n", "  return GSHGEOM_OK;\n}\n", AT),
    ("ald-ignores-ferror", "src", A, "(fflush(out) == 0 && !ferror(out))", "(fflush(out) == 0)", AT),
    ("ald-ignores-header-error", "src", A, "(unsigned long long)gshgeom_measured_hash(m)) < 0) {\n    return GSHGEOM_ERR_IO;", "(unsigned long long)gshgeom_measured_hash(m)) < 0) {\n    return GSHGEOM_OK;", AT),
    ("ald-null-stream", "src", A, "GSHGEOM_Result gshgeom_measured_dump(const GSHGEOM_Measured * m, FILE * out) {\n  if (!m || !out) {", "GSHGEOM_Result gshgeom_measured_dump(const GSHGEOM_Measured * m, FILE * out) {\n  if (!m) {", AT),

    # --- story 9: the fill rule (AD-4) --------------------------------------
    ("fill-value-even-odd", "hdr", FILH, "GSHGEOM_FILL_EVEN_ODD = 1,", "GSHGEOM_FILL_EVEN_ODD = 2,", FILT),
    ("fill-value-nonzero", "hdr", FILH, "GSHGEOM_FILL_NONZERO = 0,", "GSHGEOM_FILL_NONZERO = 3,", FILT),
    ("fill-name-nonzero", "src", FILS, '"non-zero"', '"nonzero"', FILT),
    ("fill-name-even-odd", "src", FILS, '"even-odd"', '"even_odd"', FILT),
    ("fill-name-unknown", "src", FILS, '"unknown fill rule"', '"unknown"', FILT),
    ("fill-names-swapped", "src", FILS, "case GSHGEOM_FILL_NONZERO:\n    return \"non-zero\";\n  case GSHGEOM_FILL_EVEN_ODD:\n    return \"even-odd\";",
     "case GSHGEOM_FILL_NONZERO:\n    return \"even-odd\";\n  case GSHGEOM_FILL_EVEN_ODD:\n    return \"non-zero\";", FILT),
    ("fill-count-names-itself", "src", FILS, "case GSHGEOM_FILL_RULE_COUNT:\n    break;", "case GSHGEOM_FILL_RULE_COUNT:\n    return \"even-odd\";", FILT),
]

# --- story 14a: the fixed-point builder and the integer spelling of the rule ----
PLANTS += [
    ("fx-move-bits-inclusive", "src", P, "gshgeom_path_move_to_fixed(\n    GSHGEOM_Path * path, int32_t x, int32_t y, unsigned bits) {\n  if (bits > GSHGEOM_FIXED_MAX_BITS) {",
     "gshgeom_path_move_to_fixed(\n    GSHGEOM_Path * path, int32_t x, int32_t y, unsigned bits) {\n  if (bits >= GSHGEOM_FIXED_MAX_BITS) {", PT),
    ("fx-move-bits-unchecked", "src", P, "gshgeom_path_move_to_fixed(\n    GSHGEOM_Path * path, int32_t x, int32_t y, unsigned bits) {\n  if (bits > GSHGEOM_FIXED_MAX_BITS) {",
     "gshgeom_path_move_to_fixed(\n    GSHGEOM_Path * path, int32_t x, int32_t y, unsigned bits) {\n  if (bits > 31u) {", PT),
    ("fx-line-bits-unchecked", "src", P, "gshgeom_path_line_to_fixed(\n    GSHGEOM_Path * path, int32_t x, int32_t y, unsigned bits) {\n  if (bits > GSHGEOM_FIXED_MAX_BITS) {",
     "gshgeom_path_line_to_fixed(\n    GSHGEOM_Path * path, int32_t x, int32_t y, unsigned bits) {\n  if (bits > 31u) {", PT),
    ("fx-quadratic-bits-unchecked", "src", P, "int32_t cx, int32_t cy, int32_t x, int32_t y, unsigned bits) {\n  if (bits > GSHGEOM_FIXED_MAX_BITS) {",
     "int32_t cx, int32_t cy, int32_t x, int32_t y, unsigned bits) {\n  if (bits > 31u) {", PT),
    ("fx-cubic-bits-unchecked", "src", P, "    unsigned bits) {\n  if (bits > GSHGEOM_FIXED_MAX_BITS) {\n    return GSHGEOM_ERR_INVALID;\n  }\n  return gshgeom_path_cubic_to(",
     "    unsigned bits) {\n  if (bits > 31u) {\n    return GSHGEOM_ERR_INVALID;\n  }\n  return gshgeom_path_cubic_to(", PT),
    ("fx-max-bits-const", "hdr", PH, "#define GSHGEOM_FIXED_MAX_BITS 16u", "#define GSHGEOM_FIXED_MAX_BITS 15u", PT),
    ("fx-unit-off", "src", P, "double unit = (double)((int32_t)1 << bits);", "double unit = (double)((int32_t)1 << bits) + 1.0;", PT),
    ("fx-unit-shift", "src", P, "double unit = (double)((int32_t)1 << bits);", "double unit = (double)((int32_t)2 << bits);", PT),
    ("fx-swapped", "src", P, "GSHGEOM_Point p = {(double)x / unit, (double)y / unit, 0.0};", "GSHGEOM_Point p = {(double)y / unit, (double)x / unit, 0.0};", PT),
    ("fx-z", "src", P, "GSHGEOM_Point p = {(double)x / unit, (double)y / unit, 0.0};", "GSHGEOM_Point p = {(double)x / unit, (double)y / unit, 1.0};", PT),
    ("fx-quadratic-order", "src", P, "      path, fixed_point(cx, cy, bits), fixed_point(x, y, bits));", "      path, fixed_point(x, y, bits), fixed_point(cx, cy, bits));", PT),
    ("fx-cubic-order", "src", P, "  return gshgeom_path_cubic_to(path, fixed_point(c1x, c1y, bits),\n      fixed_point(c2x, c2y, bits), fixed_point(x, y, bits));",
     "  return gshgeom_path_cubic_to(path, fixed_point(c2x, c2y, bits),\n      fixed_point(c1x, c1y, bits), fixed_point(x, y, bits));", PT),
    ("sc-null-out", "src", DV, "if (!out || (flip != 0 && flip != 1)) {", "if (flip != 0 && flip != 1) {", DT),
    ("sc-flip-two", "src", DV, "if (!out || (flip != 0 && flip != 1)) {", "if (!out || (flip != 0 && flip != 1 && flip != 2)) {", DT),
    ("sc-flip-ignored", "src", DV, "  if (flip) {\n    v = -v;\n  }", "  (void)flip;", DT),
    ("sc-flip-inverted", "src", DV, "  if (flip) {\n    v = -v;\n  }", "  if (!flip) {\n    v = -v;\n  }", DT),
    ("sc-half-low", "src", DV, "int64_t n = v + 32768;", "int64_t n = v + 32767;", DT),
    ("sc-half-high", "src", DV, "int64_t n = v + 32768;", "int64_t n = v + 32769;", DT),
    ("sc-negative-truncates", "src", DV, ": -(int64_t)(((uint64_t)(-n) + 65535u) >> 16);", ": -(int64_t)((uint64_t)(-n) >> 16);", DT),
    ("sc-negative-ceiling-off", "src", DV, ": -(int64_t)(((uint64_t)(-n) + 65535u) >> 16);", ": -(int64_t)(((uint64_t)(-n) + 65536u) >> 16);", DT),
    ("sc-shift-15", "src", DV, "int64_t q = n >= 0 ? (int64_t)((uint64_t)n >> 16)", "int64_t q = n >= 0 ? (int64_t)((uint64_t)n >> 15)", DT),
    ("sc-max-wide", "src", DV, "if (q > INT32_MAX || q < INT32_MIN) {", "if (q > (int64_t)INT32_MAX + 1 || q < INT32_MIN) {", DT),
    ("sc-min-narrow", "src", DV, "if (q > INT32_MAX || q < INT32_MIN) {", "if (q > INT32_MAX || q <= INT32_MIN) {", DT),
    ("sc-min-wide", "src", DV, "if (q > INT32_MAX || q < INT32_MIN) {", "if (q > INT32_MAX || q < (int64_t)INT32_MIN - 128) {", DT),
]


# --- strokes (stroke.h, SG1 to SG18): ids start with sk- --------------------------------------------------------
# A test list may name "probe": the plant is then caught too when the cross-target probe no longer prints the
# committed hash (a vertex moved by less than the oracle can see moves the probe's hash).
SKP = "src/stroke/params.c"
SKE = "src/stroke/expand.c"
SKJ = "src/stroke/join.c"
SKB = "src/stroke/body.c"
SKA = "src/stroke/arc.c"
SKM = "src/stroke/emit.c"
SKO = "src/stroke/outline.c"
SKS = "src/exact/sum.c"
SKD = "src/stroke/dash.c"
SKX = "src/stroke/dxf.c"
STD = ["testDash"]
STK = ["testStroke"]
STM = ["testStroke_matrix"]
STO = ["testStroke", "testStroke_matrix"]
STROKE_PLANTS = [
    # SG1 to SG3: the parameters
    ("sk-init-limit", "src", SKP, "  stroke->miter_limit = GSHGEOM_STROKE_MITER_LIMIT_DEFAULT;", "  stroke->miter_limit = 10.0;", STK),
    ("sk-init-cap", "src", SKP, "  stroke->cap = GSHGEOM_STROKE_CAP_BUTT;", "  stroke->cap = GSHGEOM_STROKE_CAP_ROUND;", STK),
    ("sk-check-null-dashes", "src", SKP, "  if (s->dash_count && !s->dashes) {\n    return GSHGEOM_ERR_INVALID;\n  }\n  if (!is_finite_value(s->width)",
     "  if (!is_finite_value(s->width)", STK),
    ("sk-check-profile-first", "src", SKP, "  if (s->dash_count && !s->dashes) {\n    return GSHGEOM_ERR_INVALID;\n  }\n  if (!is_finite_value(s->width)",
     "  if (s->width_profile) {\n    return GSHGEOM_ERR_UNSUPPORTED;\n  }\n  if (s->dash_count && !s->dashes) {\n    return GSHGEOM_ERR_INVALID;\n  }\n  if (!is_finite_value(s->width)", STK),
    ("sk-check-enum-range", "src", SKP, "s->join >= (uint32_t)GSHGEOM_STROKE_JOIN_COUNT", "s->join > (uint32_t)GSHGEOM_STROKE_JOIN_COUNT", STK),
    ("sk-check-negative-width", "src", SKP, "  if (s->width < 0.0) {\n    return GSHGEOM_ERR_INVALID;", "  if (s->width < -1.0) {\n    return GSHGEOM_ERR_INVALID;", STK),
    ("sk-check-reserved", "src", SKP, "    if (s->_reserved[i]) {", "    if (s->_reserved[i] && i < 4) {", STK),
    ("sk-value-min-open", "src", SKP, "(m >= GSHGEOM_STROKE_VALUE_MIN && m <= GSHGEOM_STROKE_VALUE_MAX)", "(m > GSHGEOM_STROKE_VALUE_MIN && m <= GSHGEOM_STROKE_VALUE_MAX)", STK),
    ("sk-value-max-open", "src", SKP, "(m >= GSHGEOM_STROKE_VALUE_MIN && m <= GSHGEOM_STROKE_VALUE_MAX)", "(m >= GSHGEOM_STROKE_VALUE_MIN && m < GSHGEOM_STROKE_VALUE_MAX)", STK),
    ("sk-miter-limit-low", "src", SKP, "!(s->miter_limit >= 1.0 &&", "!(s->miter_limit > 1.0 &&", STK),
    ("sk-dashes-max", "src", SKP, "  if (canonical_dash_count(s) > GSHGEOM_STROKE_DASHES_MAX) {", "  if (canonical_dash_count(s) > 2 * GSHGEOM_STROKE_DASHES_MAX) {", STK),
    ("sk-odd-not-doubled", "src", SKP, "    return s->dash_count > SIZE_MAX / 2u ? SIZE_MAX : s->dash_count * 2u;", "    return s->dash_count;", STK),
    ("sk-rule-reaches-output", "src", SKP, "  out->alignment_rule = in->alignment == (uint32_t)GSHGEOM_STROKE_ALIGN_CENTER ? GSHGEOM_FILL_NONZERO\n                                                                                : in->alignment_rule;",
     "  out->alignment_rule = in->alignment_rule;", STK),
    ("sk-degenerate-reaches-butt", "src", SKP, "  out->degenerate = in->cap == (uint32_t)GSHGEOM_STROKE_CAP_BUTT ? (uint32_t)GSHGEOM_STROKE_DEGENERATE_DRAW : in->degenerate;", "  out->degenerate = in->degenerate;", STK),
    ("sk-negative-zero-width", "src", SKP, "  out->width = in->width == 0.0 ? 0.0 : in->width;", "  out->width = in->width;", STK),
    # SG15: the determinant
    ("sk-det-one-product", "src", SKP, "  if (gshgeom_orient2d(o, c0, c1, &sign) != GSHGEOM_OK) {\n    return GSHGEOM_ERR_INTERNAL;\n  }",
     "  sign = map->m[0][0] * map->m[1][1] != 0.0;\n  (void)o;\n  (void)c0;\n  (void)c1;", STK),
    # SG14 to SG16: statuses and parts
    ("sk-zero-width-draws", "src", SKE, "  if (c.st.width == 0.0) {", "  if (0) {", STK),
    ("sk-status-order", "src", SKE, "  if (!invertible) {\n    o->status", "  if (!invertible && c.st.width != 0.0) {\n    o->status", STK),
    ("sk-aligned-at-w", "src", SKE, "    c.h = to_aligned ? 2.0 * h : h;", "    c.h = h;", STM),
    ("sk-open-subpath-clipped", "src", SKE, "    int to_aligned = aligned && rd.subs[i].closed && rd.subs[i].seg_count;", "    int to_aligned = aligned && rd.subs[i].seg_count;", STM),
    ("sk-complement-inverted", "src", SKE, "  o->complement = c.st.alignment == GSHGEOM_STROKE_ALIGN_OUTSIDE;", "  o->complement = c.st.alignment == GSHGEOM_STROKE_ALIGN_INSIDE;", STK),
    ("sk-clip-rule-ignored", "src", SKE, "  o->rule = c.st.alignment_rule;", "  o->rule = GSHGEOM_FILL_NONZERO;", STK),
    ("sk-clip-open-subpaths", "src", SKE, "    if (!s->closed || s->seg_count == 0) {\n      continue;\n    }\n    r = gshgeom_path_move_to(mp", "    if (s->seg_count == 0) {\n      continue;\n    }\n    r = gshgeom_path_move_to(mp", STK),
    ("sk-local-built-in-device", "src", SKE, "    if (!c->local) {\n      r = gshgeom_map_apply_curve", "    if (1) {\n      r = gshgeom_map_apply_curve", STM),
    ("sk-draw-omit-swapped", "src", SKE, "    if (!s->draws || c->st.degenerate == GSHGEOM_STROKE_DEGENERATE_OMIT) {", "    if (!s->draws || c->st.degenerate == GSHGEOM_STROKE_DEGENERATE_DRAW) {", STM),
    ("sk-lone-move-draws", "src", SKE, "    if (!s->draws || c->st.degenerate", "    if (c->st.degenerate", STM),
    ("sk-start-cap-forward", "src", SKE, "    r = gshgeom_sk_cap(c, seg[0].p[0], gshgeom_sk_scale(-1.0, d));", "    r = gshgeom_sk_cap(c, seg[0].p[0], d);", STM),
    ("sk-closed-start-join-dropped", "src", SKE, "    } else if (r == GSHGEOM_OK && s->closed) {\n      r = gshgeom_sk_join(c, &seg[k], &seg[0]);", "    } else if (r == GSHGEOM_OK && s->closed && 0) {\n      r = gshgeom_sk_join(c, &seg[k], &seg[0]);", STM),
    ("sk-tangent-wrong-neighbour", "src", SKE, "  for (size_t j = i; j-- > 0;) {", "  for (size_t j = 0; j > 0;) {", STM),
    ("sk-tangent-always-x", "src", SKE, "static GSHGEOM_SK_Vec effective_tangent(const GSHGEOM_SK_Ctx * c, const Reader * rd, size_t i) {\n",
     "static GSHGEOM_SK_Vec effective_tangent(const GSHGEOM_SK_Ctx * c, const Reader * rd, size_t i) {\n  if (rd) {\n    return gshgeom_sk_v(1.0, 0.0);\n  }\n", STM),
    ("sk-tangent-device-under-local", "src", SKE, "  return c->local ? gshgeom_sk_v(1.0, 0.0) : gshgeom_sk_v(c->map.m[0][0], c->map.m[1][0]);",
     "  return gshgeom_sk_v(c->map.m[0][0], c->map.m[1][0]);", STM),
    ("sk-magnitude-rule-removed", "src", SKE, "  if (r == GSHGEOM_OK && !magnitudes_ok(&c, &rd,", "  if (r == GSHGEOM_OK && 0 && !magnitudes_ok(&c, &rd,", STK),
    # SG8, SG9: joins
    ("sk-turn-side-inverted", "src", SKJ, "  int outer_left = sigma > 0;", "  int outer_left = sigma < 0;", STM),
    ("sk-smooth-gets-a-join", "src", SKJ, "  if (sigma == 0 && ds > 0) {\n    return GSHGEOM_OK; /* a smooth vertex: no join */\n  }", "  (void)ds;", STM),
    ("sk-sg9-strict", "src", SKJ, "    return D <= 0; /* q >= 0: |P||a||b| <= L^2 q */", "    return D < 0; /* q >= 0: |P||a||b| <= L^2 q */", STO),
    ("sk-sg9-filter-loose", "src", SKJ, "    double bound = abs * 0x1p-40;", "    double bound = abs * 0x1p-80;", STO),
    ("sk-sg9-dot-sign", "src", SKJ, "  int qs = dot_sign(q, v, s);", "  int qs = -dot_sign(q, v, s);", STO),
    ("sk-miter-in-device", "src", SKJ, "    if (gshgeom_sk_miter_within(q, v, s, L)) {\n      return miter_clip(c, anchor, ah, bh, oa, ob, 1);",
     "    if (gshgeom_sk_miter_within(gshgeom_sk_lin(c, q), gshgeom_sk_lin(c, v), gshgeom_sk_lin(c, s), L)) {\n      return miter_clip(c, anchor, ah, bh, oa, ob, 1);", STK),
    ("sk-miter-clip-at-lw", "src", SKJ, "  double la = (L * h - gshgeom_sk_dot(A, u))", "  double la = (2.0 * L * h - gshgeom_sk_dot(A, u))", STM),
    ("sk-miter-clip-never-clips", "src", SKJ, "    return miter_clip(c, anchor, ah, bh, oa, ob, gshgeom_sk_miter_within(q, v, s, L));\n  case GSHGEOM_STROKE_JOIN_ROUND:",
     "    return miter_clip(c, anchor, ah, bh, oa, ob, 1);\n  case GSHGEOM_STROKE_JOIN_ROUND:", STM),
    ("sk-round-join-full-disc", "src", SKJ, "    return gshgeom_sk_sector(c, anchor, oa, ob, 0.0, 0.0);", "    return gshgeom_sk_disc(c, anchor, 0.0, 0.0);", STM),
    ("sk-round-sector-inner", "src", SKJ, "    return gshgeom_sk_sector(c, anchor, oa, ob, 0.0, 0.0);",
     "    return gshgeom_sk_sector(c, anchor, gshgeom_sk_scale(-1.0, oa), gshgeom_sk_scale(-1.0, ob), 0.0, 0.0);", STM),
    ("sk-bevel-dropped", "src", SKJ, "  case GSHGEOM_STROKE_JOIN_BEVEL:\n    return poly(c, anchor, bevel, 3, 0.0);", "  case GSHGEOM_STROKE_JOIN_BEVEL:\n    return GSHGEOM_OK;", STM),
    ("sk-reversal-round-dropped", "src", SKJ, "    case GSHGEOM_STROKE_JOIN_ROUND:\n      return half_disc(c, anchor, n, ah);", "    case GSHGEOM_STROKE_JOIN_ROUND:\n      return GSHGEOM_OK;", STM),
    ("sk-reversal-rect-behind", "src", SKJ, "      GSHGEOM_SK_Vec hn = gshgeom_sk_scale(h, n), f = gshgeom_sk_scale(L * h, ah);", "      GSHGEOM_SK_Vec hn = gshgeom_sk_scale(h, n), f = gshgeom_sk_scale(-L * h, ah);", STM),
    ("sk-square-cap-half", "src", SKJ, "  GSHGEOM_SK_Vec hm = gshgeom_sk_scale(h, m), hd = gshgeom_sk_scale(h, dh);\n  GSHGEOM_SK_Vec rect[4]", "  GSHGEOM_SK_Vec hm = gshgeom_sk_scale(h, m), hd = gshgeom_sk_scale(h / 2, dh);\n  GSHGEOM_SK_Vec rect[4]", STM),
    ("sk-butt-degenerate-dot", "src", SKJ, "  if (c->st.cap == GSHGEOM_STROKE_CAP_BUTT) {\n    return GSHGEOM_OK;\n  }\n  GSHGEOM_SK_Vec anchor;\n  GSHGEOM_Result r = gshgeom_sk_anchor(c, p, &anchor);",
     "  GSHGEOM_SK_Vec anchor;\n  GSHGEOM_Result r = gshgeom_sk_anchor(c, p, &anchor);", STM),
    ("sk-degenerate-square-small", "src", SKJ, "  GSHGEOM_SK_Vec hm = gshgeom_sk_scale(c->h, m), hd = gshgeom_sk_scale(c->h, dh);", "  GSHGEOM_SK_Vec hm = gshgeom_sk_scale(c->h, m), hd = gshgeom_sk_scale(c->h / 2, dh);", STM),
    ("sk-miter-tip-error-ignored", "src", SKJ, "  *err = 0x1p-49 * ratio * ratio * len * c->norm * (1.0 + 0x1p-40);", "  *err = 0.0 * ratio * len;", ["testStroke_matrix", "testStroke"]),
    # the arcs join (SG8, the 2018 CR)
    ("sk-arcs-round-at-ge", "src", SKA, "  return gshgeom_exact_sum_sign(&s) > 0;", "  return gshgeom_exact_sum_sign(&s) >= 0;", STM),
    ("sk-arcs-round-threshold-halved", "src", SKA, "  double nm1 = (double)(n - 1), nn = (double)n;", "  double nm1 = (double)(n - 1) * 2.0, nn = (double)n;", STM),
    ("sk-arcs-miter-clip-fallback-dropped", "src", SKJ, "    if (r == GSHGEOM_OK && mc) {", "    if (0) {", STM),
    ("sk-arcs-enclosing-as-apart", "src", SKA, "      int enclosing = dz < diff;", "      int enclosing = 0 * (int)dz;", STM),
    ("sk-arcs-clip-at-lw", "src", SKA, "  double dn = gshgeom_sk_norm(X), Lh = L * h;", "  double dn = gshgeom_sk_norm(X), Lh = 2.0 * L * h;", STM),
    ("sk-arcs-limit-along-chord", "src", SKA, "      if (2.0 * beta * rho > Lh) {", "      if (dn > Lh + 0.0 * beta * rho) {", STM),
    ("sk-arcs-reversal-one-side", "src", SKJ, "        r = gshgeom_sk_arcs_join(c, in, out, 0, 1, &mc);", "        r = GSHGEOM_OK;", STM),
    ("sk-arcs-sign-of-curvature", "src", SKA, "  e.kappa = inner_positive ? k : -k;", "  e.kappa = inner_positive ? -k : k;", STM),
    ("sk-arcs-grow-direction", "src", SKA, "    ea.nu = rr > 0.0 ? gshgeom_sk_scale(-1.0, oa) : oa;", "    ea.nu = rr > 0.0 ? oa : gshgeom_sk_scale(-1.0, oa);", STM),
    ("sk-arcs-fan-from-vertex", "src", SKA, "  size_t n = b->n;\n  GSHGEOM_Result res = GSHGEOM_OK;\n  double * ys = NULL;",
     "  size_t n = b->n;\n  GSHGEOM_Result res = GSHGEOM_OK;\n  for (size_t fi = 0; fi < n && res == GSHGEOM_OK; fi++) {\n    GSHGEOM_SK_Vec tri[3] = {gshgeom_sk_v(0.0, 0.0), b->e[fi].p, b->e[fi].q};\n    res = gshgeom_sk_emit(c, anchor, tri, 3, dd, dd);\n  }\n  if (n) {\n    return res;\n  }\n  double * ys = NULL;", STM),
    ("sk-arcs-crossings-not-heights", "src", SKA, "      if (t > 0.0 && t < 1.0 && u > 0.0 && u < 1.0) {", "      if (0 && t < u) {", STM),
    ("sk-arcs-ray-ignored", "src", SKA, "    int ok0 = gshgeom_sk_dot(gshgeom_sk_sub(cand[0], le->a), le->dir) >= 0.0;\n    int ok1 = gshgeom_sk_dot(gshgeom_sk_sub(cand[1], le->a), le->dir) >= 0.0;", "    int ok0 = 1;\n    int ok1 = 1;\n    (void)le;", STM),
    ("sk-arcs-grow-behind", "src", SKA, "            if (gshgeom_sk_dot(gshgeom_sk_sub(foot, le->a), le->dir) >= 0.0) {", "            if (1 || foot.x == 0.0) {", STM),
    ("sk-arcs-positive-winding-only", "src", SKA, "      if (before == 0 && w != 0) {\n        start = i;\n      } else if (before != 0 && w == 0) {", "      if (before <= 0 && w > 0) {\n        start = i;\n      } else if (before > 0 && w <= 0) {", STM),
    ("sk-arcs-closer-the-other-way", "src", SKA, "      pick = side > 0.0 ? 1 : 0;", "      pick = side > 0.0 ? 0 : 1;", STM),
    ("sk-arcs-reversal-tie-last", "src", SKA, "      if (along_key(&ea, cand[1]) < along_key(&ea, cand[0])) {", "      if (along_key(&ea, cand[1]) > along_key(&ea, cand[0])) {", STM),
    ("sk-arcs-arc-not-clipped", "src", SKA, "          if (clip_f(b, m) > 0.0) {", "          if (clip_f(b, m) > 1e300) {", STM),
    ("sk-arcs-chords-uncancelled", "src", SKA, "    for (int rep = 0; rep < (net < 0 ? -net : net); rep++) {", "    net = net != 0 ? (net > 0 ? 1 : -1) : 0;\n    for (int rep = 0; rep < (net < 0 ? -net : net); rep++) {", STM),
    ("sk-arcs-polyline-unbounded", "src", SKA, "  bd.allowed = room * 0.999 - eps;", "  bd.allowed = room * 4.0;", STM),
    ("sk-body-coincident-end-ignored", "src", SKB, "  p.zn = s->p[s->n].x == s->p[s->n - 1].x && s->p[s->n].y == s->p[s->n - 1].y && s->w[s->n] == s->w[s->n - 1];", "  p.zn = 0;", STM),
    ("sk-body-coincident-start-ignored", "src", SKB, "  p.z0 = s->p[1].x == s->p[0].x && s->p[1].y == s->p[0].y && s->w[1] == s->w[0];", "  p.z0 = 0;", STM + ["probe"]),
    # SG6: the body
    ("sk-band-fold-test", "src", SKB, "  if (!(kappa * c->h < cosb)) {", "  if (!(kappa * c->h < 4.0 * cosb)) {", STM + ["probe"]),
    ("sk-band-delta-halved", "src", SKB, "  double delta = 0.5 * lambda * sind * UP;", "  double delta = 0.01 * lambda * sind * UP;", STM),
    ("sk-cusp-disc-dropped", "src", SKB, "  if (gshgeom_sk_dot(t0, t1) < 0.0) {", "  if (0) {", STM),
    ("sk-fan-overshoot-ignored", "src", SKB, "        double extra = (dm + c->h_dev * over) * UP;", "        double extra = (dm + 0.0 * over) * UP;", STM),
    ("sk-fan-sector-one-side", "src", SKB, "          *r = gshgeom_sk_sector(c, centre, gshgeom_sk_scale(-1.0, na), gshgeom_sk_scale(-1.0, nb), d_in, extra);", "          *r = GSHGEOM_OK;", STM),
    ("sk-depth-cap-accepts", "src", SKB, "    return GSHGEOM_ERR_LIMIT; /* a forced piece is refused, never drawn coarser */", "    return GSHGEOM_OK; /* a forced piece is refused, never drawn coarser */", STK),
    ("sk-line-normal-right", "src", SKB, "    GSHGEOM_SK_Vec hn = gshgeom_sk_scale(c->h, gshgeom_sk_left(gshgeom_sk_unit(a)));", "    GSHGEOM_SK_Vec hn = gshgeom_sk_scale(c->h / 2, gshgeom_sk_left(gshgeom_sk_unit(a)));", STM),
    ("sk-split-error-dropped", "src", SKB, "      en[i] = SPLIT_C * en[i] + SPLIT * en[i + 1] + 3.0 * U * (SPLIT_C * ma + SPLIT * mb);", "      en[i] = 0.0 * (ma + mb);", ["testStroke_matrix", "probe"]),
    # arcs (circular pieces) and emission
    ("sk-arc-norm-ignored", "src", SKA, "    double d = (cpq > -0.5) ? r * c->norm * (gmath_sqrt(2.0 / (1.0 + cpq)) - 1.0) * UP : -1.0;", "    double d = (cpq > -0.5) ? r * (gmath_sqrt(2.0 / (1.0 + cpq)) - 1.0) * UP : -1.0;", STM),
    ("sk-arc-inscribed", "src", SKA, "      double k = r / (1.0 + cpq);", "      double k = r / 2.0;", STM + ["probe"]),
    ("sk-norm-unscaled", "src", SKM, "  double x = ax / m, y = ay / m;\n  return m * gmath_sqrt(x * x + y * y);", "  return gmath_sqrt(a.x * a.x + a.y * a.y);", STK),
    ("sk-build-min-dropped", "src", SKE, "        if (dm != 0.0 && dm < GSHGEOM_SK_BUILD_MIN) {\n          return 0;\n        }", "        (void)dm;", STK),
    ("sk-inflation-dropped", "src", SKM, "  double rc = (d_in + c->rho) * 1.0824;", "  double rc = d_in * 1.0824;", ["testStroke"]),
    ("sk-ring-wound-backwards", "src", SKM, "    if (whole) {\n      r = write_ring(c->part, v, h);",
     "    if (whole) {\n      for (size_t i = 0; i < h / 2; i++) {\n        GSHGEOM_Device_Point tmp = v[i];\n        v[i] = v[h - 1 - i];\n        v[h - 1 - i] = tmp;\n      }\n      r = write_ring(c->part, v, h);", STM),
    ("sk-truncate-not-round", "src", SKM, "    r = gshgeom_round_26_6(hull[i].x, &v[i].x);\n    if (r == GSHGEOM_OK) {\n      r = gshgeom_round_26_6(hull[i].y, &v[i].y);\n    }",
     "    v[i].x = (int32_t)(int64_t)(hull[i].x * 64.0);\n    v[i].y = (int32_t)(int64_t)(hull[i].y * 64.0);", ["testStroke", "probe"]),
    ("sk-fits-ignores-octagon", "src", SKM, "  return d_out + GSHGEOM_SK_OCT * (d_in + c->rho) + c->rho <= c->budget;", "  return d_out <= c->budget * 4.0;", STM),
    # SG18: the hash
    ("sk-hash-kind", "src", SKO, "#define OUTLINE_KIND 0x4b525453u", "#define OUTLINE_KIND 0x4b525454u", STK),
    ("sk-hash-parts-swapped", "src", SKO, "  gshgeom_hash_u64(&h, gshgeom_flat_path_hash(o->body));\n  gshgeom_hash_u64(&h, gshgeom_flat_path_hash(o->aligned));",
     "  gshgeom_hash_u64(&h, gshgeom_flat_path_hash(o->aligned));\n  gshgeom_hash_u64(&h, gshgeom_flat_path_hash(o->body));", STK),
    ("sk-hash-status-dropped", "src", SKO, "  gshgeom_hash_u8(&h, (uint8_t)o->status);\n", "", STK),
    ("sk-hash-complement-dropped", "src", SKO, "  gshgeom_hash_u8(&h, (uint8_t)(o->complement ? 1 : 0));\n", "", STK),
    ("sk-dump-absent-part", "src", SKO, '    return fprintf(out, "part %s -\\n", name) < 0 ? GSHGEOM_ERR_IO : GSHGEOM_OK;', '    return fprintf(out, "part %s\\n", name) < 0 ? GSHGEOM_ERR_IO : GSHGEOM_OK;', STK),
    ("sk-dump-flush", "src", SKO, "  return (fflush(out) == 0 && !ferror(out)) ? GSHGEOM_OK : GSHGEOM_ERR_IO;\n}", "  (void)fflush(out);\n  return GSHGEOM_OK;\n}", STK),
    # the wide exact sum
    ("sk-sum-base", "src", SKS, "#define SUM_BASE 12888L", "#define SUM_BASE 12000L", STK),
    ("sk-sum-borrow", "src", SKS, "      borrow = (cur >> 63) & 1u;\n    }\n    for (i = q + len + 1;", "      borrow = 0;\n    }\n    for (i = q + len + 1;", STK),
    # SG19 to SG30: dashes and DXF linetypes
    ('sk-dash-refused-again', "src", SKP, '  size_t dash_count = canonical_dash_count(in);\n', '  size_t dash_count = canonical_dash_count(in);\n  if (dash_count) {\n    return GSHGEOM_ERR_UNSUPPORTED;\n  }\n', STK),
    ('sk-dash-ignored', "src", SKP, '  out->dash_count = dash_count;\n', '  out->dash_count = 0;\n', STO),
    ('sk-dash-offset-ignored', "src", SKP, '  out->dash_offset = dash_count && in->dash_offset != 0.0 ? in->dash_offset : 0.0;\n', '  out->dash_offset = 0.0;\n', STO),
    ('sk-dash-scope-ignored', "src", SKP, '  out->dash_scope = dash_count ? in->dash_scope : (uint32_t)GSHGEOM_STROKE_DASH_PER_SUBPATH;\n', '  out->dash_scope = (uint32_t)GSHGEOM_STROKE_DASH_PER_SUBPATH;\n', STO),
    ('sk-dash-negative-offset', "src", SKD, '    r = q == 0.0 ? 0.0 : s - q;\n', '    r = q;\n', STD),
    ('sk-dash-phase-strict', "src", SKD, '    if (s_in >= o) {\n', '    if (s_in > o) {\n', STD),
    ('sk-dash-stop-inclusive', "src", SKD, '      if (!(start < length)) {\n', '      if (!(start <= length)) {\n', STD),
    ('sk-dash-last-not-clipped', "src", SKD, '      if (end > length) {\n        end = length;\n      }\n', '', STD),
    ('sk-dash-period-step', "src", SKD, '      q += 1.0;\n', '      q += 2.0;\n', STD),
    ('sk-dash-odd-not-repeated', "src", SKD, '  return (count % 2u) ? count * 2u : count;\n', '  return count;\n', STD),
    ('sk-dash-negative-kept', "src", SKD, '    if (dashes[i] < 0.0) {\n      return 0;\n    }\n', '', STD),
    ('sk-dash-zero-kept', "src", SKD, '  if (!any) {\n    return 0;\n  }\n  return (count', '  return (count', STD),
    ('sk-dash-first-piece-unclipped', "src", SKD, '    double e = remaining < length ? remaining : length;\n', '    double e = remaining;\n', STD),
    ('sk-dash-steps-uncounted', "src", SKD, '    if (++*steps > max_steps) {\n', '    if (*steps > max_steps) {\n', STD),
    ('sk-dash-bound-offset', "src", SKD, '  if (!in_band(offset)) {\n    return GSHGEOM_ERR_LIMIT;\n  }\n', '', STD),
    ('sk-dash-scope-swapped', "src", SKD, '  D->per_path = c->st.dash_scope == (uint32_t)GSHGEOM_STROKE_DASH_PATH;\n', '  D->per_path = c->st.dash_scope != (uint32_t)GSHGEOM_STROKE_DASH_PATH;\n', STM),
    ('sk-dash-path-phase', "src", SKD, 'gshgeom_sk_dash_phase(D->offset + ms.begin, D->period)', 'gshgeom_sk_dash_phase(D->offset, D->period)', STM),
    ('sk-dash-device-in-model-units', "src", SKD, '  if (!c->local) {\n    /* the mapped centerline', '  if (0) {\n    /* the mapped centerline', STM),
    ('sk-dash-join-dropped', "src", SKD, '    if (r == GSHGEOM_OK && k < k1) {\n', '    if (0) {\n', STM),
    ('sk-dash-start-cap-forward', "src", SKD, '  r = gshgeom_sk_cap(c, first.p[0], gshgeom_sk_scale(-1.0, d0));\n', '  r = gshgeom_sk_cap(c, first.p[0], d0);\n', STM),
    ('sk-dash-end-cap-backward', "src", SKD, '    r = gshgeom_sk_cap(c, lastp.p[lastp.n], seg_dir(&segs[k1], t1));\n', '    r = gshgeom_sk_cap(c, lastp.p[lastp.n], gshgeom_sk_scale(-1.0, seg_dir(&segs[k1], t1)));\n', STM),
    ('sk-dash-dot-dropped-under-omit', "src", SKD, '  return gshgeom_sk_degenerate(cb->c, seg_point(s, t), seg_dir(s, t));\n', '  if (cb->c->st.degenerate == (uint32_t)GSHGEOM_STROKE_DEGENERATE_OMIT) {\n    return GSHGEOM_OK;\n  }\n  return gshgeom_sk_degenerate(cb->c, seg_point(s, t), seg_dir(s, t));\n', STM),
    ('sk-dash-degenerate-ignores-omit', "src", SKD, '  if (sub->seg_count == 0 && c->st.degenerate == (uint32_t)GSHGEOM_STROKE_DEGENERATE_OMIT) {\n    return GSHGEOM_OK;\n  }\n', '', STM),
    ('sk-dash-dir-end-at-start', "src", SKD, '  if (t <= 0.0) {\n    (void)gshgeom_sk_start_dir(s, &d);\n    return d;\n  }\n  if (t >= 1.0) {\n    (void)gshgeom_sk_end_dir(s, &d);\n    return d;\n  }\n', '  if (t <= 0.0) {\n    (void)gshgeom_sk_end_dir(s, &d);\n    return d;\n  }\n  if (t >= 1.0) {\n    (void)gshgeom_sk_end_dir(s, &d);\n    return d;\n  }\n', STM),
    ('sk-dash-split-weights', "src", SKD, '        W[i] = u * W[i] + t * W[i + 1];\n', '        W[i] = W[i];\n', STM),
    ('sk-dash-cut-at-half', "src", SKD, '    double u = a > 0.0 ? (b - a) / (1.0 - a) : b;\n', '    double u = b;\n', STM),
    ('sk-dash-gate-removed', "src", SKD, '  if (!(eps * c->norm <= 0x1p-20)) {\n', '  if (0) {\n', STM),
    ('sk-dash-gate-measurer-term', "src", SKD, 'double eps = (0x1p-46 * poly', 'double eps = (0x1p-60 * poly', STM),
    ('sk-dash-measure-loose', "src", SKD, 'GSHGEOM_MEASURE_TOLERANCE_MIN, &D->m', 'GSHGEOM_MEASURE_TOLERANCE_MAX, &D->m', STM),
    ('sk-dxf-no-cyclic-merge', "src", SKX, '  } else if (k > 1 && first_down == last_down) {\n', '  } else if (0) {\n', STD),
    ('sk-dxf-offset-lost', "src", SKX, '  double off = lead;\n', '  double off = 0.0;\n', STD),
    ('sk-dxf-rotation-offset', "src", SKX, '    off = lead - first;\n', '    off = lead;\n', STD),
    ('sk-dxf-dot-merged', "src", SKX, '  if (first_down && dashes[0] == 0.0 && last_down) {\n', '  if (0) {\n', STD),
    ('sk-dxf-dot-not-pen-down', "src", SKX, '    if (v >= 0.0) {\n      any_down = 1;\n', '    if (v > 0.0) {\n      any_down = 1;\n', STD),
    ('sk-dxf-scale-ignored', "src", SKX, '    double v = elements[i] * scale;\n    if (!finite(v)', '    double v = elements[i];\n    if (!finite(v)', STD),
    ('sk-dxf-solid-kept-dashed', "src", SKX, '  if (!any_up) {\n', '  if (0) {\n', STD),
    ('sk-dxf-capacity-unchecked', "src", SKX, '  if (need > capacity || need > GSHGEOM_STROKE_DASHES_MAX) {\n', '  if (need > GSHGEOM_STROKE_DASHES_MAX) {\n', STD),
    ('sk-dxf-rotation-skipped', "src", SKX, '  if (!first_down) {\n    /* rotate', '  if (0) {\n    /* rotate', STD),
    ('sk-dash-spent-dash-kept', "src", SKD, '(remaining > 0.0 || d_index == 0.0)', '(1)', STD),
    ('sk-dash-dot-at-phase-dropped', "src", SKD, '(remaining > 0.0 || d_index == 0.0)', '(remaining > 0.0)', STD),
    ('sk-dash-wrap-never', "src", SKD, '(k.first_flags & GSHGEOM_SK_POS_BEGINS_BEFORE) && (k.last_flags & GSHGEOM_SK_POS_ENDS_AFTER)', '0', STM),
    ('sk-dash-wrap-needs-only-last', "src", SKD, '(k.first_flags & GSHGEOM_SK_POS_BEGINS_BEFORE) && (k.last_flags & GSHGEOM_SK_POS_ENDS_AFTER)', '(k.last_flags & GSHGEOM_SK_POS_ENDS_AFTER)', STO),
    ('sk-dash-wrap-needs-only-first', "src", SKD, '(k.first_flags & GSHGEOM_SK_POS_BEGINS_BEFORE) && (k.last_flags & GSHGEOM_SK_POS_ENDS_AFTER)', '(k.first_flags & GSHGEOM_SK_POS_BEGINS_BEFORE)', STO),
    ('sk-dash-wrap-whole-not-solid', "src", SKD, '      if (k.n == 1 || !(k.first_end < k.last_start)) {\n', '      if (0) {\n', STO),
    ('sk-dash-wrap-flag-at-phase', "src", SKD, 'int flags = (s_ex < o ?', 'int flags = (s_ex <= o ?', STO),
    ('sk-dash-wrap-end-at-length', "src", SKD, 'flags = GSHGEOM_SK_POS_ENDS_AFTER;\n      }\n      GSHGEOM_Result r = fn(', 'flags = GSHGEOM_SK_POS_ENDS_AFTER;\n      }\n      flags = GSHGEOM_SK_POS_ENDS_AFTER;\n      GSHGEOM_Result r = fn(', STO),
    ('sk-dash-wrap-no-join', "src", SKD, '      r = gshgeom_sk_join(c, &segs[k], &segs[f + (r0 + j0 + jj + 1) % N]);', '      r = jj + 1 == m && m > 0 && wrap && jj + 1 > 0 && k == f + N - 1 ? GSHGEOM_OK : gshgeom_sk_join(c, &segs[k], &segs[f + (r0 + j0 + jj + 1) % N]);', STO),
    ('sk-dash-wrap-held-first-dropped', "src", SKD, '    cb->held_end = b; /* the first piece goes with the last dash */\n    return GSHGEOM_OK;\n', '    cb->held_end = 0.0;\n    return GSHGEOM_OK;\n', STO),
    ('sk-dxf-no-nudge', "src", SKX, '        off = next_above(at);', '        off = at;', STD),
    ('sk-dxf-nudge-down', "src", SKX, '  bits++;\n', '  bits--;\n', STD),
    ('sk-dxf-nudge-sum-short', "src", SKX, '  for (size_t i = 0; i + 2 < k; i++) {\n        at += dashes[i];', '  for (size_t i = 0; i + 3 < k; i++) {\n        at += dashes[i];', STD),
    ('sk-dxf-nudge-always', "src", SKX, '    int dot_last = lead == 0.0 && dashes[k - 1] == 0.0;', '    int dot_last = lead == 0.0;', STD),
    ('sk-dxf-lone-dot-no-offset', "src", SKX, '        off = period * 0x1p-52; /* the only dot', '        off = period * 0x1p-60; /* the only dot', STD),
    # SG31: the canonical map (strokes story S4)
    ('sg31-floor-ceiling', 'src', SKP, '    double f = gmath_floor(map->t[i]);', '    double f = -gmath_floor(-map->t[i]);', ["testStroke_canonical"]),
    ('sg31-fraction-unrounded', 'src', SKP, '    if (d >= 1.0) {\n      f += 1.0;\n      d = 0.0;\n    }', '', ["testStroke_canonical"]),
    ('sg31-fraction-not-reset', 'src', SKP, '      d = 0.0;\n', '      d = d;\n', ["testStroke_canonical"]),
    ('sg31-negative-zero', 'src', SKP, '    frac[i] = d + 0.0; /* -0 is +0 */', '    frac[i] = d;', ["testStroke_canonical"]),
    ('sg31-one-axis', 'src', SKP, '  out.t[1] = frac[1];\n', '  out.t[1] = map->t[1];\n', ["testStroke_canonical"]),
    ('sg31-shift-x-sign', 'src', SKP, '  *shift_x = (int64_t)whole[0];', '  *shift_x = -(int64_t)whole[0];', ["testStroke_canonical"]),
    ('sg31-shift-y-swapped', 'src', SKP, '  *shift_y = (int64_t)whole[1];', '  *shift_y = (int64_t)whole[0];', ["testStroke_canonical"]),
    ('sg31-limit-open', 'src', SKP, '  const double limit = 0x1p40;', '  const double limit = 0x1p80;', ["testStroke_canonical"]),
    ('sg31-limit-y-unchecked', 'src', SKP, ' && map->t[1] > -limit && map->t[1] < limit)) {', ')) {', ["testStroke_canonical"]),
    ('sg31-linear-part-dropped', 'src', SKP, '  GSHGEOM_Map out = *map;', '  GSHGEOM_Map out = *map;\n  out.m[0][1] = 0.0;', ["testStroke_canonical"]),
    ('sg31-3d-accepted', 'src', SKP, 'gshgeom_map_check(map) != GSHGEOM_OK || map->dim != 2) {\n    return GSHGEOM_ERR_INVALID;\n  }\n  const double limit', 'gshgeom_map_check(map) != GSHGEOM_OK) {\n    return GSHGEOM_ERR_INVALID;\n  }\n  const double limit', ["testStroke_canonical"]),
    ('sg31-kind-stale', 'src', SKP, '  out.kind = frac[0] != 0.0 || frac[1] != 0.0 ?', '  out.kind = 1 ?', ["testStroke_canonical"]),
    ('sk-join-class-ignored', "src", SKJ, '    if (out->vtx_class == 1) {\n', '    if (0) {\n', STM),
    ('sk-join-class-on-mapped-points', "src", SKE, '    pm.p[k] = pm.model_p[k];\n', '    pm.p[k] = pm.p[k];\n', STM),
    ('sk-join-class-next-mapped', "src", SKE, '    nm.p[k] = nm.model_p[k];\n', '    nm.p[k] = nm.p[k];\n', STM),
    ('sk-join-class-wrap-dropped', "src", SKE, '    if (s->closed && s->seg_count > 0) {\n      classify(', '    if (0) {\n      classify(', STM),
]

PLANTS += STROKE_PLANTS

# Strokes story S6: the curve outline for writers (stroke_curves.h, SG32 to SG43).  The ids start with cv-.
SCL = "src/stroke_curves/lower.c"
TCL = ["testStroke_curves_lower"]
SCC = "src/stroke_curves/contour.c"
SCO = "src/stroke_curves/output.c"
TCC = ["testStroke_curves"]
TCO = ["testStroke_curves_oracle"]
TCB = ["testStroke_curves", "testStroke_curves_oracle"]
CURVE_PLANTS = [
    # SG40: the lowering
    ('cv-lower-handle-short', "src", SCL, '4.0 * w / (3.0 * (1.0 + w))', '2.0 * w / (3.0 * (1.0 + w))', TCL),
    ('cv-lower-never-bisects', "src", SCL, '  if (bound <= tol) {\n    if (w < 1.0) {', '  if (1) {\n    if (w < 1.0) {', TCL),
    ('cv-lower-half-weight', "src", SCL, '  double wh = gmath_sqrt((1.0 + w) / 2.0);', '  double wh = w;', TCL),
    ('cv-lower-quadratic-inexact', "src", SCL, '    return cubic(out, lerp(p0, p1, 2.0 / 3.0), lerp(p2, p1, 2.0 / 3.0), p2);\n  }\n  double bound;', '    return cubic(out, lerp(p0, p1, 0.5), lerp(p2, p1, 0.5), p2);\n  }\n  double bound;', TCL),
    ('cv-lower-weight-not-standardized', "src", SCL, '      double ws = w[1] / gmath_sqrt(w[0] * w[2]);', '      double ws = w[1];', TCL),
    ('cv-lower-error-too-loose', "src", SCL, '    double e = (1.0 - w) * (1.0 - w) * (1.0 - w) / (54.0 * (1.0 + w));', '    double e = (1.0 - w) * (1.0 - w) * (1.0 - w) / (5400.0 * (1.0 + w));', TCL),
    ('cv-lower-ellipse-sag-dropped', "src", SCL, 'sag * sag * w * w / (s2 * s2)', '0.0', TCL),
    ('cv-lower-general-spread-dropped', "src", SCL, '  double bound = (wmax - wmin) / wmin * d * UPL;', '  double bound = (wmax - wmin) / wmax * d * UPL * 0.001;', TCL),
    # SG37: the exact pieces
    ('cv-arc-as-chord', "src", SCC, '  return gshgeom_skc_conic(c, p0, p1, p2, gmath_sqrt((1.0 + cpq) / 2.0));', '  return gshgeom_skc_line(c, p0, p2);', TCB),
    ('cv-arc-weight-wrong', "src", SCC, 'gmath_sqrt((1.0 + cpq) / 2.0)', 'gmath_sqrt((1.0 + cpq) / 3.0)', TCB),
    ('cv-arc-control-point-at-chord', "src", SCC, 'dev(c, centre, h / (1.0 + cpq), gshgeom_sk_add(p, q))', 'dev(c, centre, h / 2.0, gshgeom_sk_add(p, q))', TCB),
    ('cv-arc-no-quarter-split', "src", SCC, '  if (cpq < 0.0) {\n    GSHGEOM_SK_Vec m;', '  if (0) {\n    GSHGEOM_SK_Vec m;', TCB),
    ('cv-conic-weight-scaled', "src", SCC, '  e.w[1] = w;\n  e.w[2] = 1.0;', '  e.w[1] = w * 0.5;\n  e.w[2] = 1.0;', TCB),
    ('cv-sector-inner-side', "src", SCC, '  if (cr < 0.0) {\n    GSHGEOM_SK_Vec t = p;\n    p = q;\n    q = t;\n  }\n  gshgeom_skc_begin(c);', '  if (cr > 0.0) {\n    GSHGEOM_SK_Vec t = p;\n    p = q;\n    q = t;\n  }\n  gshgeom_skc_begin(c);', TCB),
    ('cv-half-disc-wrong-side', "src", SCC, '  if (gshgeom_sk_cross(m, d) < 0.0) {\n    m = gshgeom_sk_scale(-1.0, m);\n  }\n  GSHGEOM_SK_Vec nm', '  m = gshgeom_sk_scale(-1.0, m);\n  GSHGEOM_SK_Vec nm', TCB),
    ('cv-contour-wound-backwards', "src", SCC, '  if (!cv->flip) {\n    r = gshgeom_path_move_to(cv->part, gp(cv->el[0].p[0]));', '  if (cv->flip) {\n    r = gshgeom_path_move_to(cv->part, gp(cv->el[0].p[0]));', TCB),
    ('cv-contour-never-flipped', "src", SKE, '      cvs.flip = sg < 0;', '      cvs.flip = 0;', TCB),
    ('cv-conic-weights-not-reversed', "src", SCC, '    w[0] = e->w[reversed ? 2 : 0];', '    w[0] = e->w[0];', TCB),
    ('cv-line-polygon-dropped', "src", SCC, '  GSHGEOM_Result r = gshgeom_path_move_to(path, gp(pts[0]));\n  for (size_t i = 1; i < n', '  GSHGEOM_Result r = gshgeom_path_move_to(path, gp(pts[1]));\n  for (size_t i = 1; i < n', TCB),
    ('cv-fan-not-inflated', "src", SKB, '          *r = gshgeom_skc_fan_sector(c, centre, na, nb, d_in * UP);', '          *r = gshgeom_skc_fan_sector(c, centre, na, nb, 0.0);', TCO),
    ('cv-cusp-disc-dropped', "src", SKB, '        *r = gshgeom_skc_disc(c, centre, c->h + d_in * UP);', '        *r = GSHGEOM_OK;', TCO),
    ('cv-fan-back-line-dropped', "src", SCC, '  double la = (-r - gshgeom_sk_dot(oa, u)) / gshgeom_sk_dot(p, u);', '  double la = (-100.0 * r - gshgeom_sk_dot(oa, u)) / gshgeom_sk_dot(p, u);', TCO),
    # SG37, SG38: the curved offsets
    ('cv-rails-swapped', "src", SKB, '    double sgn = side == 0 ? 1.0 : -1.0;\n    GSHGEOM_SK_Vec b0', '    double sgn = side == 0 ? -1.0 : 1.0;\n    GSHGEOM_SK_Vec b0', TCB),
    ('cv-fold-band-accepted', "src", SKB, '  int fold = !(kappa * c->h < bi->cosb);', '  int fold = 0;', TCO),
    ('cv-fold-uses-the-smooth-band', "src", SKB, '    if (st == 2 && p->has0 && p->has1 && chord_curve(sg, p, &bi, &r)) {\n      return r;\n    }', '    if (st == 2 && p->has0 && p->has1 && band_curve(sg, p, &bi, &r)) {\n      return r;\n    }', TCO),
    ('cv-budget-not-divided-by-the-norm', "src", SKB, '  double d = eps * c->norm + corner;\n  if (!(eps <= 0.25 * h)', '  double d = eps + corner;\n  if (!(eps <= 0.25 * h)', TCC),
    ('cv-fit-error-bound-too-small', "src", SKB, '  double eps = o4 / 384.0 * UP;', '  double eps = o4 / 384.0e6 * UP;', TCO),
    ('cv-fit-error-fourth-derivative-dropped', "src", SKB, '  return (G[4] + h * F[4]) * UP;', '  return (h * F[4]) * UP * 1e-6;', TCO),
    ('cv-chord-bound-too-small', "src", SKB, '  double d2 = second_bound(G, m, h) / 8.0 * UP;\n', '  double d2 = second_bound(G, m, h) / 8.0e6 * UP;\n', TCO),
    ('cv-chord-lobes-dropped', "src", SKB, '  if (s1 * s2 < 0 && s3 * s4 < 0) {', '  if (0) {', TCO),
    ('cv-chord-wrong-diagonal', "src", SKB, '  if (a1 * a2 > 0 || (a1 == 0 && a2 == 0)) {', '  if (a1 * a2 <= 0) {', TCO),
    ('cv-rail-start-moved', "src", SKB, '    GSHGEOM_SK_Vec b0 = gshgeom_sk_add(P0, gshgeom_sk_scale(sgn * h, gshgeom_sk_left(t0)));', '    GSHGEOM_SK_Vec b0 = gshgeom_sk_add(P0, gshgeom_sk_scale(sgn * h * 1.05, gshgeom_sk_left(t0)));', TCO),
    ('cv-rail-handle-sign', "src", SKB, '    rail[side][1] = gshgeom_sk_add(b0, gshgeom_sk_scale(1.0 / 3.0, o0));', '    rail[side][1] = gshgeom_sk_sub(b0, gshgeom_sk_scale(1.0 / 3.0, o0));', TCO),
    ('cv-tangent-turning-dropped', "src", SKB, '    GSHGEOM_SK_Vec o0 = gshgeom_sk_add(c1s, gshgeom_sk_scale(sgn * h, gshgeom_sk_left(tp0)));', '    GSHGEOM_SK_Vec o0 = c1s;', TCO),
    ('cv-shared-tangent-from-one-side', "src", SKB, '      l->dir1 = r->dir0 = dv;', '      l->dir1 = dv;\n      r->dir0 = gshgeom_sk_v(dv.x, dv.y + 1e-3 * dv.x);', TCO),
    ('cv-arcs-join-monotone-unchecked', "src", SCC, '  if (!side_monotone(&ea, sigma, A, B) || !side_monotone(&eb, sigma, A, B)) {', '  if (0) {', TCO),
    ('cv-arcs-join-orientation-ignored', "src", SCC, '  if (sigma < 0.0) {\n    reverse_elements(c->cv);\n  }', '', TCO),
    ('cv-arcs-join-arc-as-chord', "src", SCC, '    r = gshgeom_skc_conic(c, dev(c, anchor, 1.0, p0), dev(c, anchor, 1.0, p1), dev(c, anchor, 1.0, p2),\n        gmath_sqrt((1.0 + cpq) / 2.0));', '    (void)p1;\n    r = gshgeom_skc_line(c, dev(c, anchor, 1.0, p0), dev(c, anchor, 1.0, p2));', TCC),
    ('cv-hulls-not-apart-ignored', "src", SKB, '  if (!hulls_apart(rail[0], rail[1], 0x1p-30 * scale + 4.0 * eps)) {', '  if (0) {', TCO),
    ('cv-band-run-never-flushed', "src", SCC, '  size_t n = cv->run_n;\n  cv->run_n = 0;\n  if (n == 0) {', '  size_t n = cv->run_n;\n  cv->run_n = 0;\n  if (1) {', TCO),
    ('cv-depth-cap-truncated', "src", SKB, '  if (depth >= sg->c->depth_cap) {\n    return GSHGEOM_ERR_LIMIT;', '  if (depth >= sg->c->depth_cap) {\n    return GSHGEOM_OK;', TCB),
    # SG32 to SG36, SG41: the parts, the forms, the hash
    ('cv-clip-flattened', "src", SKE, '      case GSHGEOM_VERB_CUBIC:\n        r = gshgeom_path_cubic_to(mp, p[1], p[2], p[3]);\n        break;\n      case GSHGEOM_VERB_RATIONAL_QUADRATIC:\n        r = gshgeom_path_rational_quadratic_to(mp, p[1], p[2], cv.w);', '      case GSHGEOM_VERB_CUBIC:\n        r = gshgeom_path_line_to(mp, p[3]);\n        break;\n      case GSHGEOM_VERB_RATIONAL_QUADRATIC:\n        r = gshgeom_path_rational_quadratic_to(mp, p[1], p[2], cv.w);', TCC),
    ('cv-clip-not-mapped', "src", SKE, '      r = gshgeom_map_apply_curve(&c->map, &el.curve, &cv);\n      if (r != GSHGEOM_OK) {\n        break;\n      }\n      const GSHGEOM_Point * p = cv.p;\n      switch (el.verb) {\n      case GSHGEOM_VERB_LINE:\n        r = gshgeom_path_line_to(mp, p[1]);\n        break;\n      case GSHGEOM_VERB_QUADRATIC:\n        r = gshgeom_path_quadratic_to(mp, p[1], p[2]);\n        break;\n      case GSHGEOM_VERB_CUBIC:\n        r = gshgeom_path_cubic_to(mp, p[1], p[2], p[3]);\n        break;\n      case GSHGEOM_VERB_RATIONAL_QUADRATIC:\n        r = gshgeom_path_rational_quadratic_to(mp, p[1], p[2], cv.w);', '      cv = el.curve;\n      const GSHGEOM_Point * p = cv.p;\n      switch (el.verb) {\n      case GSHGEOM_VERB_LINE:\n        r = gshgeom_path_line_to(mp, p[1]);\n        break;\n      case GSHGEOM_VERB_QUADRATIC:\n        r = gshgeom_path_quadratic_to(mp, p[1], p[2]);\n        break;\n      case GSHGEOM_VERB_CUBIC:\n        r = gshgeom_path_cubic_to(mp, p[1], p[2], p[3]);\n        break;\n      case GSHGEOM_VERB_RATIONAL_QUADRATIC:\n        r = gshgeom_path_rational_quadratic_to(mp, p[1], p[2], cv.w);', TCC),
    ('cv-clip-rule-ignored', "src", SKE, '    cu->rule = c.st.alignment_rule;', '    cu->rule = GSHGEOM_FILL_NONZERO;', TCC),
    ('cv-clip-complement-ignored', "src", SKE, '    cu->complement = c.st.alignment == GSHGEOM_STROKE_ALIGN_OUTSIDE;', '    cu->complement = 0;', TCC),
    ('cv-cubic-form-at-the-wrong-budget', "src", SKE, 'gshgeom_stroke_curves_lower_cubic(*parts[k], cr->t / 4.0, &c.alloc, &c.limits, &low);', 'gshgeom_stroke_curves_lower_cubic(*parts[k], cr->t / 2.0, &c.alloc, &c.limits, &low);', TCB),
    ('cv-cubic-form-not-lowered', "src", SKE, '  if (r == GSHGEOM_OK && cr && cr->form == GSHGEOM_STROKE_CURVES_CUBIC) {', '  if (r == GSHGEOM_OK && cr && cr->form == (GSHGEOM_Stroke_Curves_Form)99) {', TCC),
    ('cv-magnitude-bound-removed', "src", SKE, '    c.mag_limit = lim < GSHGEOM_SK_MAG_LIMIT ? lim : GSHGEOM_SK_MAG_LIMIT;', '    c.mag_limit = GSHGEOM_SK_MAG_LIMIT;', TCC),
    ('cv-budget-widened', "src", SKE, '    cvs.budget = cr->t / 2.0;', '    cvs.budget = cr->t * 40.0;', TCO),
    ('cv-tolerance-range-open-below', "src", SKE, '  if (!(tolerance >= GSHGEOM_STROKE_CURVES_MIN_TOLERANCE) ||', '  if (!(tolerance >= 0.0) ||', TCC),
    ('cv-tolerance-range-open-above', "src", SKE, '!(tolerance <= GSHGEOM_STROKE_CURVES_MAX_TOLERANCE) ||', '!(tolerance <= 1e30) ||', TCC),
    ('cv-form-range-unchecked', "src", SKE, '      (unsigned)form >= (unsigned)GSHGEOM_STROKE_CURVES_FORM_COUNT) {', '      0) {', TCC),
    ('cv-hash-parts-swapped', "src", SCO, '  gshgeom_hash_u64(&h, gshgeom_path_hash(o->body));\n  gshgeom_hash_u64(&h, gshgeom_path_hash(o->aligned));', '  gshgeom_hash_u64(&h, gshgeom_path_hash(o->aligned));\n  gshgeom_hash_u64(&h, gshgeom_path_hash(o->body));', TCC),
    ('cv-hash-kind-wrong', "src", SCO, '#define CURVES_KIND 0x43525453u', '#define CURVES_KIND 0x4b525453u', TCC),
    ('cv-dump-form-missing', "src", SCO, 'status=%u form=%u align=%u', 'status=%u form=%u form=%u', TCC),
    ('cv-dump-flush-unchecked', "src", SCO, '  return (fflush(out) == 0 && !ferror(out)) ? GSHGEOM_OK : GSHGEOM_ERR_IO;', '  return GSHGEOM_OK;', TCC),
    ('cv-hash-tolerance-included', "src", SCO, '  gshgeom_hash_u8(&h, (uint8_t)(o->complement ? 1 : 0));', '  gshgeom_hash_u8(&h, (uint8_t)(o->complement ? 1 : 0));\n  gshgeom_hash_u8(&h, (uint8_t)o->form);', TCC),
]
PLANTS += CURVE_PLANTS

def run(cmd, shell=False, env_extra=None, cwd=None):
    environment = dict(os.environ)
    if env_extra:
        environment.update(env_extra)
    return subprocess.run(cmd, shell=shell, capture_output=True, text=True,
                          env=environment, cwd=cwd or ROOT)


def compile_all(destdir, srcs, extra=(), pre=()):
    """Compile every source to destdir/<name>.o; return the object list."""
    destdir.mkdir(parents=True, exist_ok=True)
    objs = []
    for src in srcs:
        obj = destdir / (src.replace("/", "_") + ".o")
        r = run(CC + list(pre) + CFLAGS + list(extra) + ["-c", src, "-o", str(obj)])
        if r.returncode != 0:
            raise RuntimeError("cannot compile %s:\n%s" % (src, r.stderr))
        objs.append(str(obj))
    return objs


def link_test(test, objs, out):
    cmd = "%s -o %s %s %s %s" % (CXX, out, test_object(test), " ".join(objs), LDFLAGS)
    r = run(cmd, shell=True)
    if r.returncode != 0:
        raise RuntimeError("cannot link %s:\n%s" % (test, r.stderr))


def _limits():
    """A planted defect may loop or allocate without end; bound both."""
    import resource
    resource.setrlimit(resource.RLIMIT_AS, (3 << 30, 3 << 30))
    resource.setrlimit(resource.RLIMIT_CPU, (120, 120))


def run_test(exe):
    environment = dict(os.environ)
    environment["LD_LIBRARY_PATH"] = LD_PATH
    try:
        r = subprocess.run([str(exe), "--gtest_brief=1"], capture_output=True, text=True,
                           env=environment, cwd=ROOT, timeout=180, preexec_fn=_limits)
    except subprocess.TimeoutExpired:
        return 1, ["[  FAILED  ] (timed out)"]
    failing = [l for l in (r.stdout + r.stderr).splitlines() if l.startswith("[  FAILED  ]")]
    return r.returncode, failing


def probe_hash(objs, d, extra=()):
    exe = d / "corpus-hash"
    cmd = "%s -o %s tools/corpus-hash.c %s %s" % (
        " ".join(CC + PROBE_CFLAGS + list(extra)), exe, " ".join(objs), PROBE_LIBS)
    r = run(cmd, shell=True)
    if r.returncode != 0:
        raise RuntimeError("cannot build the probe:\n%s" % r.stderr)
    c = run(["tools/check-corpus-hash.sh", str(exe), "tests/corpus/geometry.hash"],
            env_extra={"LD_LIBRARY_PATH": LD_PATH})
    return c.returncode, (c.stdout + c.stderr).strip()


def apply(text, old, new, label):
    n = text.count(old)
    if n != 1:
        raise RuntimeError("plant %s: the text to replace occurs %d times, not once:\n%s" % (label, n, old))
    if old == new:
        raise RuntimeError("plant %s: the replacement changes nothing" % label)
    return text.replace(old, new)


# The sources that include src/core/grow.h.
INTERNAL_HEADER_USERS = ["src/device/device.c", "src/flatten/flatten.c"]


def do_plant(plant, control_objs):
    pid, kind, fname, old, new, tests = plant
    d = WORK / pid
    shutil.rmtree(d, ignore_errors=True)
    d.mkdir(parents=True)
    srcs = sources()
    if kind == "src":
        text = apply((ROOT / fname).read_text(), old, new, pid)
        planted = d / "planted" / fname
        planted.parent.mkdir(parents=True, exist_ok=True)
        # A source may include a neighbour by a relative path ("../exact/exact.h"),
        # so the planted file sits in a copy of the whole tree.
        shutil.copytree(ROOT / "src", d / "planted" / "src", dirs_exist_ok=True)
        planted.write_text(text)
        obj = d / "planted.o"
        r = run(CC + CFLAGS + ["-c", str(planted), "-o", str(obj)])
        if r.returncode != 0:
            raise RuntimeError("plant %s does not compile:\n%s" % (pid, r.stderr))
        base = str(ROOT / fname)
        objs = [str(obj) if o.endswith(fname.replace("/", "_") + ".o") else o for o in control_objs]
        assert str(obj) in objs, pid
    elif kind == "srch":
        # An internal header: change it in a copy of src/ and rebuild what
        # includes it.
        shutil.copytree(ROOT / "src", d / "planted" / "src", dirs_exist_ok=True)
        planted = d / "planted" / fname
        planted.write_text(apply((ROOT / fname).read_text(), old, new, pid))
        objs = list(control_objs)
        for dep in INTERNAL_HEADER_USERS:
            obj = d / (dep.replace("/", "_") + ".o")
            r = run(CC + CFLAGS + ["-c", str(d / "planted" / dep), "-o", str(obj)])
            if r.returncode != 0:
                raise RuntimeError("plant %s does not compile:\n%s" % (pid, r.stderr))
            objs = [str(obj) if o.endswith(dep.replace("/", "_") + ".o") else o for o in objs]
    else:
        text = apply((ROOT / fname).read_text(), old, new, pid)
        hdr = d / "include" / Path(fname).relative_to("include")
        shutil.copytree(ROOT / "include", d / "include")
        hdr.write_text(text)
        objs = compile_all(d / "obj", srcs, pre=["-I" + str(d / "include")])
    caught_by = []
    for test in tests:
        if test == "probe":
            rc, msg = probe_hash(objs, d)
            if rc != 0:
                caught_by.append(("probe", [msg.splitlines()[-1] if msg else "hash differs"]))
            continue
        exe = d / test
        link_test(test, objs, exe)
        rc, failing = run_test(exe)
        if rc != 0:
            caught_by.append((test, failing[:1]))
    return pid, caught_by


def main():
    shutil.rmtree(WORK, ignore_errors=True)
    WORK.mkdir(parents=True)
    srcs = sources()
    failures = 0

    # The control: everything passes.
    control = compile_all(WORK / "control", srcs)
    for test in TESTS:
        exe = WORK / ("control-" + test)
        link_test(test, control, exe)
        rc, failing = run_test(exe)
        if rc != 0:
            print("FAIL  control: %s fails on the unplanted build: %s" % (test, failing))
            return 1
    rc, msg = probe_hash(control, WORK)
    if rc != 0:
        print("FAIL  control: the probe does not print the committed hash: %s" % msg)
        return 1
    print("ok    control: %d tests and the probe pass on the unplanted build" % len(TESTS))

    ids = [p[0] for p in PLANTS]
    if len(set(ids)) != len(ids):
        print("FAIL  duplicate plant ids")
        return 1
    plants = [p for p in PLANTS if not ONLY or re.search(ONLY, p[0])]
    if PROBE_ONLY:
        plants = []
    if not plants and not PROBE_ONLY:
        print("FAIL  PLANT_ONLY matches no plant")
        return 1

    with concurrent.futures.ThreadPoolExecutor(max_workers=JOBS) as pool:
        futures = {pool.submit(do_plant, p, control): p for p in plants}
        results = {}
        for fut in concurrent.futures.as_completed(futures):
            p = futures[fut]
            try:
                results[p[0]] = fut.result()[1]
            except Exception as e:  # a plant that cannot be built is a defect of the harness
                results[p[0]] = e
    for p in plants:
        r = results[p[0]]
        if isinstance(r, Exception):
            print("FAIL  %-34s %s" % (p[0], r))
            failures += 1
        elif not r:
            print("FAIL  %-34s not caught by %s" % (p[0], " ".join(p[5])))
            failures += 1
        else:
            print("ok    %-34s caught by %s" % (p[0], ", ".join("%s (%s)" % (t, (f[0][14:70] if f else "crash")) for t, f in r)))

    if ONLY and not PROBE_ONLY:
        print("%d plants run (PLANT_ONLY=%s), %d not caught" % (len(plants), ONLY, failures))
        return 1 if failures else 0

    # The probe's own plants.
    (WORK / "bit").mkdir(parents=True, exist_ok=True)
    for name, flag in (("probe-planted-bit", "-DSHAPE_GEOM_PLANT_BIT"),
                       ("probe-planted-int-bit", "-DSHAPE_GEOM_PLANT_BIT_INT")):
        (WORK / name).mkdir(parents=True, exist_ok=True)
        rc, msg = probe_hash(control, WORK / name, extra=[flag])
        if rc == 0:
            print("FAIL  %s: the probe printed the committed hash with a bit flipped" % name)
            failures += 1
        else:
            print("ok    %s: %s" % (name, msg.splitlines()[-1]))
    for pid, flags in FLAG_PLANTS:
        d = WORK / pid
        d.mkdir(parents=True, exist_ok=True)
        objs = compile_all(d / "obj", srcs, extra=shlex.split(flags))
        rc, msg = probe_hash(objs, d, extra=shlex.split(flags))
        if rc == 0:
            print("FAIL  %-34s the probe did not see %s" % (pid, flags))
            failures += 1
        else:
            print("ok    %-34s the probe's hash differs: %s" % (pid, msg.splitlines()[-1]))

    print("%d plants, %d not caught" % (len(PLANTS) + 2 + len(FLAG_PLANTS), failures))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
