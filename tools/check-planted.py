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
    ("verb-code-quadratic", "hdr", PH, "GSHGEOM_VERB_QUADRATIC = 2,", "GSHGEOM_VERB_QUADRATIC = 7,", PT),
    ("verb-code-close", "hdr", PH, "GSHGEOM_VERB_CLOSE = 6,", "GSHGEOM_VERB_CLOSE = 8,", PT),
    ("kind-code-cubic", "hdr", CH, "GSHGEOM_CURVE_CUBIC = 2,", "GSHGEOM_CURVE_CUBIC = 6,", CT + PT),

    # --- the dump -----------------------------------------------------------
    ("dump-hash-width", "src", P, "hash=%016llx", "hash=%llx", PT),
    ("dump-number-width", "src", P, 'return fprintf(out, " %016llx"', 'return fprintf(out, " %llx"', PT),
    ("dump-no-z", "src", P, "if (n >= 0 && dim == 3) n = put_f64(out, p->z);", "if (n >= 0 && dim == 4) n = put_f64(out, p->z);", PT),
    ("dump-name-order", "src", P, '"M", "L", "Q", "C", "RQ", "RC", "Z"}', '"M", "L", "C", "Q", "RQ", "RC", "Z"}', PT),
    ("dump-no-weights", "src", P, "      for (size_t k = 0; k < verb_weights(e->verb); k++) {\n        if (put_f64",
     "      for (size_t k = 0; k < verb_weights(e->verb) / 5; k++) {\n        if (put_f64", PT),
    ("dump-ignores-flush", "src", P, "return fflush(out) == 0 ? GSHGEOM_OK : GSHGEOM_ERR_IO;", "return GSHGEOM_OK;", PT),
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
    ("default-limit-const", "hdr", KH, "((size_t)1 << 22)", "((size_t)1 << 21)", PT),
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


# Plants caught by the probe's hash alone.
FLAG_PLANTS = [
    ("fused-multiply-add", "-mfma -ffp-contract=fast"),
]


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
        planted.write_text(text)
        obj = d / "planted.o"
        r = run(CC + CFLAGS + ["-c", str(planted), "-o", str(obj)])
        if r.returncode != 0:
            raise RuntimeError("plant %s does not compile:\n%s" % (pid, r.stderr))
        base = str(ROOT / fname)
        objs = [str(obj) if o.endswith(fname.replace("/", "_") + ".o") else o for o in control_objs]
        assert str(obj) in objs, pid
    else:
        text = apply((ROOT / fname).read_text(), old, new, pid)
        hdr = d / "include" / Path(fname).relative_to("include")
        shutil.copytree(ROOT / "include", d / "include")
        hdr.write_text(text)
        objs = compile_all(d / "obj", srcs, pre=["-I" + str(d / "include")])
    caught_by = []
    for test in tests:
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

    with concurrent.futures.ThreadPoolExecutor(max_workers=JOBS) as pool:
        futures = {pool.submit(do_plant, p, control): p for p in PLANTS}
        results = {}
        for fut in concurrent.futures.as_completed(futures):
            p = futures[fut]
            try:
                results[p[0]] = fut.result()[1]
            except Exception as e:  # a plant that cannot be built is a defect of the harness
                results[p[0]] = e
    for p in PLANTS:
        r = results[p[0]]
        if isinstance(r, Exception):
            print("FAIL  %-34s %s" % (p[0], r))
            failures += 1
        elif not r:
            print("FAIL  %-34s not caught by %s" % (p[0], " ".join(p[5])))
            failures += 1
        else:
            print("ok    %-34s caught by %s" % (p[0], ", ".join("%s (%s)" % (t, (f[0][14:70] if f else "crash")) for t, f in r)))

    # The probe's own plants.
    (WORK / "bit").mkdir(parents=True, exist_ok=True)
    rc, msg = probe_hash(control, WORK / "bit", extra=["-DSHAPE_GEOM_PLANT_BIT"])
    if rc == 0:
        print("FAIL  probe-planted-bit: the probe printed the committed hash with a bit flipped")
        failures += 1
    else:
        print("ok    probe-planted-bit: %s" % msg.splitlines()[-1])
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

    print("%d plants, %d not caught" % (len(PLANTS) + 1 + len(FLAG_PLANTS), failures))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
