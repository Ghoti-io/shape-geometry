/*
 * SPDX-License-Identifier: LGPL-3.0-only
 *
 * Copyright (C) 2026 Corey Pennycuff
 *
 * This file is part of Ghoti.io Shape Geometry.
 *
 * Ghoti.io Shape Geometry is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3 as
 * published by the Free Software Foundation.
 *
 * Ghoti.io Shape Geometry is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * @file
 *
 * Path building, mapping, flattening and measuring driven by bytes (story 16 of
 * the shape family's first milestone; CAP-15, AD-19): one entry point for
 * libFuzzer and for the replay gate (`make fuzz-replay`).
 *
 * Two header bytes choose the dimension, the limits (the defaults, or caps
 * small enough to be met) and the allocation that fails.  The rest is a stream
 * of ops, each an op byte and argument bytes: every verb of the model path in
 * `double` and in fixed point, a map built from translations, scales, rotations
 * in every unit and raw matrices, and the passes that read the path: bounds,
 * hash and dump, the 26.6 device path, the integer flattener at the default
 * tolerance and at another, the measured path and points along it, and the
 * orientation predicates.  Arcs (arc.h) are built by all three constructors (the endpoint form, the centre and
 * conjugate forms, the tangent form and the bulge) and by appending a raw record whose ends agree with its centre
 * form to a chosen bound or disagree past it; one pass maps the last arc, lowers every arc to its pieces, copies
 * the path, and reads every arc back (the element is refused, the pieces chain end to end and meet the arc's own
 * ends, the vertices are one per element); the device path and the measured path take the arcs the ops built.  A
 * coordinate is drawn from a table of the degenerate
 * values (zeros of both signs, subnormals, the 26.6 range's edges, the map
 * limit, values past it, NaN and the infinities) or from the bytes.
 *
 * After each pass the library's stated properties are checked, and a violation
 * is an `abort`:
 *
 * - a call returns only the results its header lists, and a refusal leaves
 *   nothing behind (the counting allocator holds nothing live at the end);
 * - a path never holds more elements than `max_elements`, a flattened path
 *   more points than `max_points`, and a hash is the same when asked twice;
 * - every flattened point lies within the control box of the device path,
 *   widened by one unit (the flattener emits points of the curve, rounded once);
 * - the adaptive orientation predicates agree with the exact ones in sign;
 * - a point along a measured path is finite, and the clamped query at the
 *   length is the end of the path;
 * - an arc is refused as an element and read back whole, its pieces (at most four) chain bit for bit from its start
 *   to its end, the lowered path holds no arc and has the arc path's box, a copy is equal, a mapped arc keeps the
 *   circular flag only under a similarity, and the endpoint form's result kind says what it added.
 *
 * The replay fails if the corpus never reaches an outcome (each refusal kind
 * and each success): tools/gen-fuzz-corpus.py writes it.
 */

#include <ghoti.io/shape-geometry/shape-geometry.h>

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(cond) \
  do { \
    if (!(cond)) { \
      fprintf(stderr, "fuzz_path: invariant failed: %s (line %d)\n", #cond, __LINE__); \
      abort(); \
    } \
  } while (0)

/* --- the outcomes the corpus must reach ----------------------------------------------------------- */

enum {
  OUT_PATH_INVALID = 0, OUT_PATH_LIMIT, OUT_DEVICE_OK, OUT_DEVICE_RANGE, OUT_FLAT_OK, OUT_FLAT_LIMIT, OUT_FLAT_TOLERANCE,
  OUT_MEASURED_OK, OUT_MEASURED_REFUSED, OUT_OOM, OUT_GRID_ALIGNED, OUT_MAP_INVALID, OUT_ORIENT_ZERO, OUT_CURVED,
  OUT_ARC, OUT_ARC_LINE, OUT_ARC_NOTHING, OUT_ARC_CIRCULAR, OUT_ARC_DISAGREE, OUT_ARC_MAPPED, OUT_ARC_MAP_REFUSED,
  OUT_ARC_LOWERED, OUT_ARC_VERTICES, OUT_ARC_COPIED, OUT_ARC_FULL_TURN, OUT_ARC_DEVICE,
  OUT_SHAPE_RECT, OUT_SHAPE_ROUNDED, OUT_SHAPE_ELLIPSE, OUT_SHAPE_FLIPPED, OUT_SHAPE_OVERLAP, OUT_SHAPE_BAD_FLIP,
  OUT_COUNT
};

static const char * const outcome_names[OUT_COUNT] = {"a non-finite point refused", "a path's element cap",
    "a device path", "a device coordinate out of the 26.6 range", "a flattened path", "a flattened path's point cap",
    "a tolerance below the minimum", "a measured path", "a measured path refused", "a failed allocation",
    "a grid-aligned map", "an invalid map refused", "a degenerate orientation (sign 0)",
    "a device path holding a curve", "an arc added", "an arc constructor's line", "an endpoint form that adds nothing",
    "a circular arc", "a record whose ends disagree with its centre form, refused", "a mapped arc", "a mapped arc refused",
    "a path with arcs lowered", "the vertices of a path with an arc", "a path with an arc copied", "a full-turn arc",
    "a device path of an arc", "a rect added", "a rounded rect with arcs added", "an ellipse added", "a flipped primitive added",
    "a rounded rect refused for an overlap", "a flip that is not 0 or 1, refused"};

static unsigned long outcomes[OUT_COUNT];

int gshgeom_fuzz_unreached(void);

int gshgeom_fuzz_unreached(void) {
  int missing = 0;
  for (int i = 0; i < OUT_COUNT; i++) {
    if (outcomes[i] == 0) {
      fprintf(stderr, "fuzz_path: the corpus never reached: %s\n", outcome_names[i]);
      missing++;
    }
  }
  return missing;
}

/* --- a counting allocator that can fail one call -------------------------------------------------- */

typedef struct Counter {
  long live;
  long calls;
  long fail_at; /* -1 never */
  int failed;
} Counter;

static int fails(Counter * c) {
  if (c->fail_at >= 0 && c->calls++ == c->fail_at) {
    c->failed = 1;
    return 1;
  }
  return 0;
}

static void * c_malloc(void * ctx, size_t n) {
  Counter * c = ctx;
  void * p = fails(c) ? NULL : malloc(n ? n : 1);
  c->live += p != NULL;
  return p;
}

static void * c_calloc(void * ctx, size_t k, size_t n) {
  Counter * c = ctx;
  void * p = fails(c) ? NULL : calloc(k ? k : 1, n ? n : 1);
  c->live += p != NULL;
  return p;
}

static void * c_realloc(void * ctx, void * old, size_t n) {
  Counter * c = ctx;
  void * p = fails(c) ? NULL : realloc(old, n ? n : 1);
  c->live += p != NULL && old == NULL;
  return p;
}

static void c_free(void * ctx, void * p) {
  ((Counter *)ctx)->live -= p != NULL;
  free(p);
}

/* --- the input ---------------------------------------------------------------------------------------- */

typedef struct Input {
  const uint8_t * data;
  size_t size;
  size_t at;
} Input;

static uint8_t next(Input * in) {
  return in->at < in->size ? in->data[in->at++] : 0;
}

/** A coordinate: a degenerate value from the table, or a small number in 1/16 from two bytes. */
static double coordinate(Input * in) {
  static const double table[] = {0.0, -0.0, 1.0, -1.0, 0.5, 1.0 / 64.0, 4.9406564584124654e-324, 2.2250738585072014e-308,
      33554431.0, 33554432.0, -33554432.0, -33554433.0, 0x1p100, -0x1p100, 0x1p101, 1e300, -1e300,
      0x1p53, 0x1p53 + 2.0, 1e-9};
  const uint8_t pick = next(in);
  if (pick < 24) {
    const size_t n = sizeof table / sizeof table[0];
    if (pick < n) {
      return table[pick];
    }
    /* NaN and the infinities, made without <math.h>: the library must refuse them. */
    volatile double zero = 0.0;
    return pick == 21 ? zero / zero : (pick == 22 ? 1.0 / zero : -1.0 / zero);
  }
  const unsigned hi = next(in);
  const unsigned lo = next(in);
  return (double)(int16_t)(uint16_t)((hi << 8) | lo) / 16.0;
}

/* Every draw is a statement of its own: C leaves the order of a call's arguments and of an initializer's
 * expressions unspecified, and compilers differ (story 14c's probe met it), so one input must mean one thing to
 * the replay's compiler and to the fuzzer's. */
static void coordinates(Input * in, double * out, int n) {
  for (int i = 0; i < n; i++) {
    out[i] = coordinate(in);
  }
}

static GSHGEOM_Point point(Input * in, int dim) {
  GSHGEOM_Point p = {0.0, 0.0, 0.0};
  p.x = coordinate(in);
  p.y = coordinate(in);
  if (dim == 3) {
    p.z = coordinate(in);
  }
  return p;
}

static void weights(Input * in, double * w, int n) {
  static const double table[] = {1.0, 0.70710678118654752, 2.0, 0.0, -1.0, 1e-300, 1e300, 0.5};
  for (int i = 0; i < n; i++) {
    w[i] = table[next(in) % (sizeof table / sizeof table[0])];
  }
}

static double absd(double x) {
  return x < 0.0 ? -x : x;
}

/** A radius: a table of the degenerate and ordinary values, or a number in 1/16 from two bytes. */
static double radius(Input * in) {
  static const double table[] = {0.0, -0.0, 1.0, 2.0, 0.5, 100.0, -3.0, 1e-9, 1e9, 1e300, 0x1p100, 0x1p-1000, 4.0, 8.0};
  const uint8_t pick = next(in);
  if (pick < 16) {
    const size_t n = sizeof table / sizeof table[0];
    volatile double zero = 0.0;
    return pick < n ? table[pick] : (pick == 14 ? zero / zero : 1.0 / zero);
  }
  const unsigned hi = next(in);
  const unsigned lo = next(in);
  return (double)(int16_t)(uint16_t)((hi << 8) | lo) / 16.0;
}

/** An angle (in whatever unit the op names): the quarter turns, odd angles, past the limit, and bytes. */
static double angle(Input * in) {
  static const double table[] = {0.0, 90.0, 180.0, 270.0, 360.0, -90.0, 720.0, 45.5, -30.0, 1.5707963267948966, 3.141592653589793,
      6.283185307179586, 0.25, 0.5, 1.0, 1e9, 1e300, 0x1p40, 4.0, 100.0, 200.0, 400.0};
  const uint8_t pick = next(in);
  if (pick < 24) {
    const size_t n = sizeof table / sizeof table[0];
    volatile double zero = 0.0;
    return pick < n ? table[pick] : (pick == 22 ? zero / zero : 1.0 / zero);
  }
  const unsigned hi = next(in);
  const unsigned lo = next(in);
  return (double)(int16_t)(uint16_t)((hi << 8) | lo) / 16.0;
}

static int32_t fixed(Input * in) {
  static const int32_t table[] = {0, 1, -1, INT32_MAX, INT32_MIN, 1 << 24, -(1 << 24), 64};
  const uint8_t pick = next(in);
  if (pick < 8) {
    return table[pick];
  }
  const unsigned hi = next(in);
  const unsigned lo = next(in);
  return (int32_t)((hi << 8) | lo) - 32768;
}

/* --- the map ---------------------------------------------------------------------------------------- */

static GSHGEOM_Result build_map(Input * in, GSHGEOM_Map * map) {
  GSHGEOM_Map step;
  GSHGEOM_Result r = gshgeom_map_identity(2, map);
  const int steps = next(in) % 4;
  for (int i = 0; i < steps && r == GSHGEOM_OK; i++) {
    double v[6];
    switch (next(in) % 6) {
      case 0:
        coordinates(in, v, 2);
        r = gshgeom_map_translate(2, v[0], v[1], 0.0, &step);
        break;
      case 1:
        coordinates(in, v, 2);
        r = gshgeom_map_scale(2, v[0], v[1], 1.0, &step);
        break;
      case 2:
      {
        const GSHGEOM_Angle_Unit unit = (GSHGEOM_Angle_Unit)(next(in) % 5);
        r = gshgeom_map_rotate(unit, coordinate(in), &step);
      }
        break;
      case 3: {
        coordinates(in, v, 6);
        r = gshgeom_map_from_matrix(2, v, v + 4, &step);
        break;
      }
      case 4:
        coordinates(in, v, 3);
        r = gshgeom_map_design_to_model(v[0], v[1], v[2], &step);
        break;
      default:
        r = gshgeom_map_flip_y(2, &step);
        break;
    }
    if (r == GSHGEOM_OK) {
      GSHGEOM_Map outer = *map;
      r = gshgeom_map_compose(&step, &outer, map);
    }
    CHECK(r == GSHGEOM_OK || r == GSHGEOM_ERR_INVALID || r == GSHGEOM_ERR_LIMIT);
  }
  if (r != GSHGEOM_OK) {
    outcomes[OUT_MAP_INVALID]++;
    return r;
  }
  int aligned = 0;
  CHECK(gshgeom_map_grid_aligned(map, &aligned, NULL) == GSHGEOM_OK);
  outcomes[OUT_GRID_ALIGNED] += aligned != 0;
  return GSHGEOM_OK;
}

/* --- the passes over a finished path ---------------------------------------------------------------- */

static void check_flat(const GSHGEOM_Device_Path * device, const GSHGEOM_Flat_Path * flat, const GSHGEOM_Limits * limits) {
  int64_t x0 = INT64_MAX, y0 = INT64_MAX, x1 = INT64_MIN, y1 = INT64_MIN;
  const size_t elements = gshgeom_device_path_element_count(device);
  for (size_t i = 0; i < elements; i++) {
    GSHGEOM_Device_Element e;
    CHECK(gshgeom_device_path_element(device, i, &e) == GSHGEOM_OK);
    for (unsigned k = 0; k < e.count; k++) {
      x0 = e.p[k].x < x0 ? e.p[k].x : x0;
      y0 = e.p[k].y < y0 ? e.p[k].y : y0;
      x1 = e.p[k].x > x1 ? e.p[k].x : x1;
      y1 = e.p[k].y > y1 ? e.p[k].y : y1;
    }
  }
  const size_t count = gshgeom_flat_path_point_count(flat);
  CHECK(count <= limits->max_points);
  const GSHGEOM_Device_Point * pts = gshgeom_flat_path_points(flat);
#ifdef GSHGEOM_FUZZ_PLANT
  const int planted = 1; /* the replay gate's plant: the first point moved past the box, which this must see */
#else
  const int planted = 0;
#endif
  for (size_t i = 0; i < count; i++) {
    const int64_t x = pts[i].x + (planted && i == 0 ? x1 - x0 + 2 : 0);
    CHECK(x >= x0 - 1 && x <= x1 + 1 && pts[i].y >= y0 - 1 && pts[i].y <= y1 + 1);
  }
  size_t seen = 0;
  for (size_t s = 0; s < gshgeom_flat_path_subpath_count(flat); s++) {
    GSHGEOM_Flat_Subpath sub;
    CHECK(gshgeom_flat_path_subpath(flat, s, &sub) == GSHGEOM_OK);
    CHECK(sub.count >= 1 && sub.first == seen && sub.first + sub.count <= count);
    seen += sub.count;
  }
  CHECK(seen == count);
  CHECK(gshgeom_flat_path_hash(flat) == gshgeom_flat_path_hash(flat));
}

static void device_passes(Input * in, const GSHGEOM_Path * path, const GSHGEOM_Allocator * alloc,
    const GSHGEOM_Limits * limits, FILE * sink) {
  GSHGEOM_Map map;
  if (build_map(in, &map) != GSHGEOM_OK) {
    return;
  }
  GSHGEOM_Device_Path * device = NULL;
  GSHGEOM_Result r = gshgeom_device_path_from_model(path, &map, alloc, limits, &device);
  CHECK(r == GSHGEOM_OK || r == GSHGEOM_ERR_LIMIT || r == GSHGEOM_ERR_INVALID || r == GSHGEOM_ERR_OOM);
  CHECK((r == GSHGEOM_OK) == (device != NULL));
  outcomes[OUT_OOM] += r == GSHGEOM_ERR_OOM;
  if (r == GSHGEOM_ERR_LIMIT) {
    outcomes[OUT_DEVICE_RANGE]++;
  }
  if (!device) {
    return;
  }
  outcomes[OUT_DEVICE_OK]++;
  {
    int has_arc = 0;
    for (size_t i = 0; i < gshgeom_path_element_count(path); i++) {
      GSHGEOM_Verb verb;
      has_arc |= gshgeom_path_verb(path, i, &verb) == GSHGEOM_OK && verb == GSHGEOM_VERB_ARC;
    }
    outcomes[OUT_ARC_DEVICE] += has_arc;
  }
  CHECK(gshgeom_device_path_element_count(device) <= limits->max_elements);
  for (size_t i = 0; i < gshgeom_device_path_element_count(device); i++) {
    GSHGEOM_Device_Element e;
    CHECK(gshgeom_device_path_element(device, i, &e) == GSHGEOM_OK);
    outcomes[OUT_CURVED] += e.verb == GSHGEOM_DEVICE_QUADRATIC || e.verb == GSHGEOM_DEVICE_CUBIC;
  }
  CHECK(gshgeom_device_path_hash(device) == gshgeom_device_path_hash(device));
  CHECK(gshgeom_device_path_dump(device, sink) == GSHGEOM_OK);

  /* The default tolerance, then the caller's (sometimes below the minimum, which is refused before anything). */
  GSHGEOM_Flat_Path * flat = NULL;
  const uint8_t which = next(in);
  const int32_t tolerance = which < 128 ? GSHGEOM_FLATTEN_TOLERANCE : (int32_t)(which & 63) - 2;
  r = which < 128 ? gshgeom_device_path_flatten(device, alloc, limits, &flat)
                  : gshgeom_device_path_flatten_with_tolerance(device, tolerance, alloc, limits, &flat);
  CHECK(r == GSHGEOM_OK || r == GSHGEOM_ERR_LIMIT || r == GSHGEOM_ERR_INVALID || r == GSHGEOM_ERR_OOM);
  CHECK(r != GSHGEOM_ERR_INVALID || tolerance < GSHGEOM_FLATTEN_MIN_TOLERANCE);
  CHECK((r == GSHGEOM_OK) == (flat != NULL));
  outcomes[OUT_FLAT_LIMIT] += r == GSHGEOM_ERR_LIMIT;
  outcomes[OUT_FLAT_TOLERANCE] += r == GSHGEOM_ERR_INVALID;
  outcomes[OUT_OOM] += r == GSHGEOM_ERR_OOM;
  if (flat) {
    outcomes[OUT_FLAT_OK]++;
    check_flat(device, flat, limits);
    CHECK(gshgeom_flat_path_dump(flat, sink) == GSHGEOM_OK);
    gshgeom_flat_path_destroy(flat);
  }
  gshgeom_device_path_destroy(device);
}

static void measure(Input * in, const GSHGEOM_Path * path, const GSHGEOM_Allocator * alloc,
    const GSHGEOM_Limits * limits, FILE * sink) {
  static const double tolerances[] = {0.0, 1e-3, 0x1p-4, 0x1p-46, 0x1p-47, 0.5, -1.0};
  const double tolerance = tolerances[next(in) % (sizeof tolerances / sizeof tolerances[0])];
  GSHGEOM_Measured * measured = NULL;
  GSHGEOM_Result r = gshgeom_measured_create(alloc, limits, path, tolerance, &measured);
  CHECK(r == GSHGEOM_OK || r == GSHGEOM_ERR_LIMIT || r == GSHGEOM_ERR_INVALID || r == GSHGEOM_ERR_OOM);
  CHECK((r == GSHGEOM_OK) == (measured != NULL));
  outcomes[OUT_OOM] += r == GSHGEOM_ERR_OOM;
  outcomes[OUT_MEASURED_REFUSED] += r == GSHGEOM_ERR_LIMIT || r == GSHGEOM_ERR_INVALID;
  if (!measured) {
    return;
  }
  outcomes[OUT_MEASURED_OK]++;
  const double length = gshgeom_measured_length(measured);
  CHECK(length >= 0.0 && length - length == 0.0);
  for (int i = 0; i < 4; i++) {
    const double s = length * (double)next(in) / 255.0;
    GSHGEOM_Measured_Frame frame;
    r = gshgeom_measured_at(measured, s, &frame);
    CHECK(r == GSHGEOM_OK || r == GSHGEOM_ERR_INVALID);
    if (r == GSHGEOM_OK) {
      CHECK(frame.point.x - frame.point.x == 0.0 && frame.point.y - frame.point.y == 0.0);
      CHECK(frame.element < gshgeom_path_element_count(path));
    }
  }
  GSHGEOM_Measured_Frame end;
  r = gshgeom_measured_at_clamped(measured, coordinate(in), &end);
  CHECK(r == GSHGEOM_OK || r == GSHGEOM_ERR_INVALID);
  CHECK(gshgeom_measured_hash(measured) == gshgeom_measured_hash(measured));
  CHECK(gshgeom_measured_dump(measured, sink) == GSHGEOM_OK);
  gshgeom_measured_destroy(measured);
}

static void orient(Input * in, int dim) {
  GSHGEOM_Point p[4];
  for (int i = 0; i < 4; i++) {
    p[i] = point(in, dim);
    if (next(in) & 1) {
      /* a near-collinear point: the hard case for a naive determinant */
      p[i].x = p[0].x + (p[1].x - p[0].x) * 0.5;
      p[i].y = p[0].y + (p[1].y - p[0].y) * 0.5;
    }
  }
  int fast = 2, exact = 2;
  GSHGEOM_Result a = gshgeom_orient2d(p[0], p[1], p[2], &fast);
  GSHGEOM_Result b = gshgeom_orient2d_exact(p[0], p[1], p[2], &exact);
  CHECK(a == b);
  if (a == GSHGEOM_OK) {
    CHECK(fast == exact && (fast == -1 || fast == 0 || fast == 1));
    outcomes[OUT_ORIENT_ZERO] += exact == 0;
  }
  a = gshgeom_orient3d(p[0], p[1], p[2], p[3], &fast);
  b = gshgeom_orient3d_exact(p[0], p[1], p[2], p[3], &exact);
  CHECK(a == b);
  CHECK(a != GSHGEOM_OK || fast == exact);
}

/* --- arcs ------------------------------------------------------------------------------------------- */

static int same_bits(double a, double b) {
  return memcmp(&a, &b, sizeof a) == 0;
}

static int same_point(GSHGEOM_Point a, GSHGEOM_Point b) {
  return same_bits(a.x, b.x) && same_bits(a.y, b.y);
}

static int finite(double v) {
  return v - v == 0.0;
}

/** The pieces of an arc chain bit for bit from its start to its end (AR8), at most four, with weights in (0, 1]. */
static void check_pieces(const GSHGEOM_Arc * a) {
  GSHGEOM_Curve pieces[GSHGEOM_ARC_MAX_PIECES];
  int n = 0;
  const GSHGEOM_Result pr = gshgeom_arc_pieces(a, pieces, &n);
  if (pr != GSHGEOM_OK) {
    fprintf(stderr, "arc pieces %d: start %a %a end %a %a c %a %a u %a %a v %a %a unit %d t0 %a sweep %a flags %u\n", (int)pr, a->start.x,
        a->start.y, a->end.x, a->end.y, a->center.x, a->center.y, a->u.x, a->u.y, a->v.x, a->v.y, (int)a->unit, a->t0, a->sweep, a->flags);
  }
  CHECK(pr == GSHGEOM_OK);
  CHECK(n >= 1 && n <= GSHGEOM_ARC_MAX_PIECES);
  CHECK(same_point(pieces[0].p[0], a->start) && same_point(pieces[n - 1].p[2], a->end));
  for (int k = 0; k < n; k++) {
    CHECK(pieces[k].kind == GSHGEOM_CURVE_RATIONAL_QUADRATIC);
    CHECK(pieces[k].w[0] == 1.0 && pieces[k].w[2] == 1.0 && pieces[k].w[1] > 0.0 && pieces[k].w[1] <= 1.0);
    if (k + 1 < n) {
      CHECK(same_point(pieces[k].p[2], pieces[k + 1].p[0]));
    }
  }
}

/** Element @p i is an arc: refused as an element, read back whole and chained. */
static void check_arc_at(const GSHGEOM_Path * path, size_t i) {
  GSHGEOM_Verb verb;
  GSHGEOM_Element el;
  GSHGEOM_Arc a;
  CHECK(gshgeom_path_verb(path, i, &verb) == GSHGEOM_OK && verb == GSHGEOM_VERB_ARC);
  CHECK(gshgeom_path_element(path, i, &el) == GSHGEOM_ERR_UNSUPPORTED);
  CHECK(gshgeom_path_arc(path, i, &a) == GSHGEOM_OK);
  CHECK(finite(a.t0) && finite(a.sweep) && a.sweep != 0.0);
  CHECK((a.flags & ~GSHGEOM_ARC_CIRCULAR) == 0u);
  check_pieces(&a);
  if (i > 0) {
    /* the start is the end of the element before it, bit for bit */
    GSHGEOM_Verb pv;
    CHECK(gshgeom_path_verb(path, i - 1, &pv) == GSHGEOM_OK);
    if (pv == GSHGEOM_VERB_ARC) {
      GSHGEOM_Arc pa;
      CHECK(gshgeom_path_arc(path, i - 1, &pa) == GSHGEOM_OK && same_point(pa.end, a.start));
    } else {
      GSHGEOM_Element pe;
      CHECK(gshgeom_path_element(path, i - 1, &pe) == GSHGEOM_OK);
      CHECK(same_point(pe.curve.p[gshgeom_curve_point_count(pe.curve.kind) - 1], a.start));
    }
  }
}

typedef struct VertexCount {
  size_t n;
  size_t last_element;
  int ordered;
} VertexCount;

static GSHGEOM_Result count_vertex(const GSHGEOM_Vertex * v, void * ctx) {
  VertexCount * c = ctx;
  c->ordered = c->ordered && (c->n == 0 || v->element == c->last_element + 1) && v->element == c->n;
  c->last_element = v->element;
  c->n++;
  CHECK(finite(v->point.x) && finite(v->point.y));
  CHECK(finite(v->direction.x) && finite(v->direction.y));
  const double len2 = v->direction.x * v->direction.x + v->direction.y * v->direction.y;
  CHECK(len2 > 0.99999 && len2 < 1.00001);
  return GSHGEOM_OK;
}

/** Every reader that takes an arc, on the path as it is: the readback, the vertices, the copy, the lowering. */
static void arc_passes(const GSHGEOM_Path * path, const GSHGEOM_Allocator * alloc, const GSHGEOM_Limits * limits) {
  const size_t n = gshgeom_path_element_count(path);
  size_t arcs = 0;
  for (size_t i = 0; i < n; i++) {
    GSHGEOM_Verb verb;
    CHECK(gshgeom_path_verb(path, i, &verb) == GSHGEOM_OK);
    if (verb == GSHGEOM_VERB_ARC) {
      check_arc_at(path, i);
      arcs++;
    } else {
      GSHGEOM_Arc a;
      CHECK(gshgeom_path_arc(path, i, &a) == GSHGEOM_ERR_INVALID);
    }
  }
  CHECK(gshgeom_path_verb(path, n, &(GSHGEOM_Verb){GSHGEOM_VERB_MOVE}) == GSHGEOM_ERR_INVALID);
  if (gshgeom_path_dimension(path) == 2) {
    VertexCount vc = {0, 0, 1};
    size_t want = 0;
    CHECK(gshgeom_path_vertex_count(path, &want) == GSHGEOM_OK && want == n);
    CHECK(gshgeom_path_vertices(path, count_vertex, &vc) == GSHGEOM_OK);
    CHECK(vc.n == n && vc.ordered);
    outcomes[OUT_ARC_VERTICES] += arcs > 0;
  }
  GSHGEOM_Path * copy = NULL;
  GSHGEOM_Result r = gshgeom_path_copy(path, alloc, limits, &copy);
  CHECK(r == GSHGEOM_OK || r == GSHGEOM_ERR_OOM);
  if (copy) {
    CHECK(gshgeom_path_equal(path, copy) && gshgeom_path_hash(path) == gshgeom_path_hash(copy));
    gshgeom_path_destroy(copy);
    outcomes[OUT_ARC_COPIED] += arcs > 0;
  }
  GSHGEOM_Path * low = NULL;
  r = gshgeom_path_lower_arcs(path, alloc, limits, &low);
  CHECK(r == GSHGEOM_OK || r == GSHGEOM_ERR_LIMIT || r == GSHGEOM_ERR_OOM);
  CHECK((r == GSHGEOM_OK) == (low != NULL));
  if (low) {
    for (size_t i = 0; i < gshgeom_path_element_count(low); i++) {
      GSHGEOM_Verb verb;
      CHECK(gshgeom_path_verb(low, i, &verb) == GSHGEOM_OK && verb != GSHGEOM_VERB_ARC);
    }
    CHECK(arcs > 0 || gshgeom_path_equal(low, path));
    GSHGEOM_Bounds a, b;
    CHECK(gshgeom_path_bounds(path, &a) == GSHGEOM_OK && gshgeom_path_bounds(low, &b) == GSHGEOM_OK);
    CHECK(a.empty == b.empty);
    CHECK(a.empty || (same_bits(a.min.x, b.min.x) && same_bits(a.min.y, b.min.y) && same_bits(a.max.x, b.max.x) &&
                      same_bits(a.max.y, b.max.y)));
    gshgeom_path_destroy(low);
    outcomes[OUT_ARC_LOWERED] += arcs > 0;
  }
}

/** Map the last arc of the path by a map built from the input (AR11). */
static void arc_map(Input * in, const GSHGEOM_Path * path) {
  GSHGEOM_Map map;
  const size_t n = gshgeom_path_element_count(path);
  GSHGEOM_Arc a;
  size_t at = n;
  for (size_t i = n; i > 0 && at == n; i--) {
    GSHGEOM_Verb verb;
    if (gshgeom_path_verb(path, i - 1, &verb) == GSHGEOM_OK && verb == GSHGEOM_VERB_ARC) {
      at = i - 1;
    }
  }
  if (build_map(in, &map) != GSHGEOM_OK || at == n || gshgeom_path_arc(path, at, &a) != GSHGEOM_OK) {
    return;
  }
  GSHGEOM_Arc out;
  GSHGEOM_Result r = gshgeom_arc_map(&map, &a, &out);
  CHECK(r == GSHGEOM_OK || r == GSHGEOM_ERR_INVALID || r == GSHGEOM_ERR_LIMIT);
  if (r != GSHGEOM_OK) {
    outcomes[OUT_ARC_MAP_REFUSED]++;
    return;
  }
  outcomes[OUT_ARC_MAPPED]++;
  /* the circular flag is kept only by a similarity, and only if it was there */
  const unsigned not_similar = GSHGEOM_MAP_SCALES_NONUNIFORM | GSHGEOM_MAP_SKEWS | GSHGEOM_MAP_SINGULAR | GSHGEOM_MAP_GENERAL;
  CHECK((out.flags & GSHGEOM_ARC_CIRCULAR) == 0u || ((a.flags & GSHGEOM_ARC_CIRCULAR) && (map.kind & not_similar) == 0u));
  CHECK((a.flags & GSHGEOM_ARC_CIRCULAR) == 0u || (map.kind & not_similar) != 0u || (out.flags & GSHGEOM_ARC_CIRCULAR));
  CHECK(out.unit == a.unit && same_bits(out.t0, a.t0) && same_bits(out.sweep, a.sweep));
  GSHGEOM_Point centre;
  double rad;
  int ccw;
  r = gshgeom_arc_circle(&out, &centre, &rad, &ccw);
  CHECK((r == GSHGEOM_OK) == ((out.flags & GSHGEOM_ARC_CIRCULAR) != 0u));
  CHECK(r != GSHGEOM_OK || (rad >= 0.0 && finite(rad)));
  /* a map is a map: mapping the mapped arc's start and end again by the inverse is not asked; its pieces chain */
  check_pieces(&out);
}

/** The result of an arc-adding op: a refusal changes nothing; an addition is an arc whose ends and flag are as asked. */
static void arc_added(const GSHGEOM_Path * path, size_t before) {
  const size_t n = gshgeom_path_element_count(path);
  CHECK(n >= before);
  if (n > before) {
    GSHGEOM_Verb verb;
    CHECK(gshgeom_path_verb(path, n - 1, &verb) == GSHGEOM_OK);
    if (verb == GSHGEOM_VERB_ARC) {
      GSHGEOM_Arc a;
      check_arc_at(path, n - 1);
      CHECK(gshgeom_path_arc(path, n - 1, &a) == GSHGEOM_OK);
      const double turn = a.unit == GSHGEOM_ANGLE_RADIANS ? 6.283185307179586
          : (a.unit == GSHGEOM_ANGLE_DEGREES ? 360.0 : (a.unit == GSHGEOM_ANGLE_GRADIANS ? 400.0 : 1.0));
      outcomes[OUT_ARC]++;
      outcomes[OUT_ARC_CIRCULAR] += (a.flags & GSHGEOM_ARC_CIRCULAR) != 0u;
      outcomes[OUT_ARC_FULL_TURN] += absd(a.sweep) >= turn;
      /* a full turn closes to the bit */
      CHECK(absd(a.sweep) < turn || same_point(a.start, a.end));
    } else {
      CHECK(verb == GSHGEOM_VERB_LINE);
      outcomes[OUT_ARC_LINE]++;
    }
  }
}


/** A primitive of shapes.h from the input: the result is one of the four the path verbs have, a refusal changes nothing (the caller
 *  checks), and an added primitive is a closed subpath of the elements the rules name, its arcs valid, positive or negative as
 *  the parity of the flips says, and its copy equal. */
static GSHGEOM_Result shape_op(Input * in, GSHGEOM_Path * path, int dim, size_t before) {
  const unsigned kind = next(in) % 3u;
  const uint8_t fl = next(in);
  int fx = fl & 1;
  int fy = (fl >> 1) & 1;
  const int bad_flip = (fl & 0x10) != 0;
  if (bad_flip) {
    fx = (fl & 0x20) ? -1 : 2;
  }
  const double x = coordinate(in);
  const double y = coordinate(in);
  const double w = radius(in);
  const double h = radius(in);
  GSHGEOM_Corner c[4];
  GSHGEOM_Result r;
  int overlap = 0;
  if (kind == 0u) {
    r = gshgeom_path_add_rect(path, x, y, w, h, fx, fy);
  } else if (kind == 1u) {
    for (int i = 0; i < 4; i++) {
      c[i].rx = radius(in);
      c[i].ry = radius(in);
    }
    r = gshgeom_path_add_rounded_rect(path, x, y, w, h, c);
    overlap = r == GSHGEOM_ERR_INVALID && finite(w) && finite(h) && w >= 0.0 && h >= 0.0;
    for (int i = 0; i < 4; i++) {
      overlap = overlap && finite(c[i].rx) && finite(c[i].ry) && c[i].rx >= 0.0 && c[i].ry >= 0.0;
    }
  } else {
    r = gshgeom_path_add_ellipse(path, x, y, w, h, fx, fy);
  }
  CHECK(dim == 2 || r == GSHGEOM_ERR_INVALID);
  outcomes[OUT_SHAPE_OVERLAP] += overlap;
  outcomes[OUT_SHAPE_BAD_FLIP] += bad_flip && kind != 1u && r == GSHGEOM_ERR_INVALID;
  if (r != GSHGEOM_OK) {
    return r;
  }
  const size_t n = gshgeom_path_element_count(path);
  const size_t added = n - before;
  CHECK(added >= 5u && added <= 10u);
  GSHGEOM_Verb verb;
  CHECK(gshgeom_path_verb(path, before, &verb) == GSHGEOM_OK && verb == GSHGEOM_VERB_MOVE);
  CHECK(gshgeom_path_verb(path, n - 1, &verb) == GSHGEOM_OK && verb == GSHGEOM_VERB_CLOSE);
  size_t arcs = 0;
  for (size_t i = before; i < n; i++) {
    CHECK(gshgeom_path_verb(path, i, &verb) == GSHGEOM_OK);
    if (verb == GSHGEOM_VERB_ARC) {
      GSHGEOM_Arc a;
      arcs++;
      check_arc_at(path, i);
      CHECK(gshgeom_path_arc(path, i, &a) == GSHGEOM_OK);
      CHECK(a.unit == GSHGEOM_ANGLE_DEGREES && a.sweep == 90.0);
      const double cross = a.u.x * a.v.y - a.u.y * a.v.x;
      CHECK(kind == 1u ? cross > 0.0 : (cross < 0.0) == (fx != fy));
      GSHGEOM_Curve pieces[GSHGEOM_ARC_MAX_PIECES];
      int np = 0;
      CHECK(gshgeom_arc_pieces(&a, pieces, &np) == GSHGEOM_OK && np == 1);
    }
  }
  CHECK(kind != 0u || (added == 5u && arcs == 0u));
  CHECK(kind != 2u || (added == 6u && arcs == 4u));
  CHECK(kind != 1u || (arcs == 0u ? added == 5u : added == 6u + arcs));
  outcomes[OUT_SHAPE_RECT] += kind == 0u;
  outcomes[OUT_SHAPE_ROUNDED] += kind == 1u && arcs > 0u;
  outcomes[OUT_SHAPE_ELLIPSE] += kind == 2u;
  outcomes[OUT_SHAPE_FLIPPED] += kind != 1u && (fx || fy);
  return r;
}

/* --- the entry point ---------------------------------------------------------------------------------- */

enum { OP_MOVE, OP_LINE, OP_QUAD, OP_CUBIC, OP_RQUAD, OP_RCUBIC, OP_CLOSE, OP_MOVE_FIXED, OP_LINE_FIXED,
       OP_QUAD_FIXED, OP_CUBIC_FIXED, OP_DEVICE, OP_MEASURE, OP_ORIENT, OP_BOUNDS, OP_ARC_ENDPOINT, OP_ARC_CENTER, OP_ARC_CONJUGATE,
       OP_ARC_TANGENT, OP_ARC_BULGE, OP_ARC_RAW, OP_ARC_PASS, OP_ARC_MAP, OP_SHAPE, OP_COUNT };

int LLVMFuzzerTestOneInput(const uint8_t * data, size_t size);

int LLVMFuzzerTestOneInput(const uint8_t * data, size_t size) {
  Input in = {data, size, 0};
  const uint8_t options = next(&in);
  const uint8_t fail_pick = next(&in);
  const int dim = (options & 1) ? 3 : 2;
  GSHGEOM_Limits limits;
  gshgeom_limits_default(&limits);
  if (options & 2) {
    limits.max_elements = 1 + ((options >> 2) & 15);
  }
  if (options & 64) {
    limits.max_points = 2 + ((options >> 2) & 15);
  }
  Counter counter = {0, 0, (options & 128) ? (long)fail_pick : -1, 0};
  const GSHGEOM_Allocator alloc = {&counter, c_malloc, c_calloc, c_realloc, c_free};
  FILE * sink = fopen("/dev/null", "wb");
  if (!sink) {
    return 0;
  }

  GSHGEOM_Path * path = NULL;
  GSHGEOM_Result r = gshgeom_path_create(&alloc, &limits, dim, &path);
  CHECK(r == GSHGEOM_OK || r == GSHGEOM_ERR_OOM);
  outcomes[OUT_OOM] += r == GSHGEOM_ERR_OOM;
  for (int ops = 0; path && in.at < in.size && ops < 512; ops++) {
    const size_t before = gshgeom_path_element_count(path);
    const uint64_t hash_before = gshgeom_path_hash(path);
    double w[4];
    r = GSHGEOM_OK;
    switch (next(&in) % OP_COUNT) {
      case OP_MOVE:
        r = gshgeom_path_move_to(path, point(&in, dim));
        break;
      case OP_LINE:
        r = gshgeom_path_line_to(path, point(&in, dim));
        break;
      case OP_QUAD: {
        GSHGEOM_Point c = point(&in, dim);
        r = gshgeom_path_quadratic_to(path, c, point(&in, dim));
        break;
      }
      case OP_CUBIC: {
        GSHGEOM_Point c1 = point(&in, dim);
        GSHGEOM_Point c2 = point(&in, dim);
        r = gshgeom_path_cubic_to(path, c1, c2, point(&in, dim));
        break;
      }
      case OP_RQUAD: {
        GSHGEOM_Point c = point(&in, dim);
        GSHGEOM_Point p = point(&in, dim);
        weights(&in, w, 3);
        r = gshgeom_path_rational_quadratic_to(path, c, p, w);
        break;
      }
      case OP_RCUBIC: {
        GSHGEOM_Point c1 = point(&in, dim);
        GSHGEOM_Point c2 = point(&in, dim);
        GSHGEOM_Point p = point(&in, dim);
        weights(&in, w, 4);
        r = gshgeom_path_rational_cubic_to(path, c1, c2, p, w);
        break;
      }
      case OP_CLOSE:
        r = gshgeom_path_close(path);
        break;
      case OP_MOVE_FIXED:
      {
        const int32_t x = fixed(&in);
        const int32_t y = fixed(&in);
        r = gshgeom_path_move_to_fixed(path, x, y, next(&in) % 40);
      }
        break;
      case OP_LINE_FIXED:
      {
        const int32_t x = fixed(&in);
        const int32_t y = fixed(&in);
        r = gshgeom_path_line_to_fixed(path, x, y, next(&in) % 40);
      }
        break;
      case OP_QUAD_FIXED: {
        int32_t v[4];
        for (int k = 0; k < 4; k++) {
          v[k] = fixed(&in);
        }
        r = gshgeom_path_quadratic_to_fixed(path, v[0], v[1], v[2], v[3], next(&in) % 40);
        break;
      }
      case OP_CUBIC_FIXED: {
        int32_t v[6];
        for (int k = 0; k < 6; k++) {
          v[k] = fixed(&in);
        }
        r = gshgeom_path_cubic_to_fixed(path, v[0], v[1], v[2], v[3], v[4], v[5], next(&in) % 40);
        break;
      }
      case OP_ARC_ENDPOINT: {
        GSHGEOM_Point e = point(&in, dim);
        const double rx = radius(&in);
        const double ry = radius(&in);
        const double phi = angle(&in);
        const GSHGEOM_Angle_Unit unit = (GSHGEOM_Angle_Unit)(next(&in) % 5);
        const uint8_t fl = next(&in);
        /* a flag that is neither 0 nor 1 is refused: one input in sixteen asks for it */
        const int large = (fl & 1) | ((fl & 0xF0) == 0xF0 ? 2 : 0), sweep = ((fl >> 1) & 1) | ((fl & 0xF0) == 0xE0 ? 2 : 0);
        GSHGEOM_Arc_Kind kind = GSHGEOM_ARC_KIND_ARC;
        r = gshgeom_path_arc_endpoint_to(path, e, rx, ry, phi, unit, large, sweep, &kind);
        if (r == GSHGEOM_OK) {
          CHECK(dim == 2);
          const size_t now = gshgeom_path_element_count(path);
          outcomes[OUT_ARC_NOTHING] += kind == GSHGEOM_ARC_KIND_NOTHING;
          CHECK(kind == GSHGEOM_ARC_KIND_NOTHING ? (now == before && gshgeom_path_hash(path) == hash_before) : now == before + 1);
          if (kind == GSHGEOM_ARC_KIND_ARC) {
            GSHGEOM_Arc a;
            CHECK(gshgeom_path_arc(path, now - 1, &a) == GSHGEOM_OK);
            /* the end is what the caller gave, to the bit; the circular flag says the radii were equal */
            CHECK(same_bits(a.end.x, e.x) && same_bits(a.end.y, e.y));
            CHECK(((a.flags & GSHGEOM_ARC_CIRCULAR) != 0u) == (absd(rx) == absd(ry)));
            CHECK(a.unit == GSHGEOM_ANGLE_RADIANS && absd(a.sweep) <= 6.2831853071795873);
          }
        }
        if (r == GSHGEOM_OK && kind != GSHGEOM_ARC_KIND_NOTHING) {
          arc_added(path, before);
        }
        break;
      }
      case OP_ARC_CENTER: {
        GSHGEOM_Point c = point(&in, dim);
        const double rx = radius(&in);
        const double ry = radius(&in);
        const double phi = angle(&in);
        const GSHGEOM_Angle_Unit unit = (GSHGEOM_Angle_Unit)(next(&in) % 5);
        const double t0 = angle(&in);
        const double sweep = angle(&in);
        const int has_end = next(&in) & 1;
        GSHGEOM_Point e = point(&in, dim);
        r = gshgeom_path_arc_center_to(path, c, rx, ry, phi, unit, t0, sweep, has_end ? &e : NULL);
        if (r == GSHGEOM_OK) {
          arc_added(path, before);
        }
        break;
      }
      case OP_ARC_CONJUGATE: {
        GSHGEOM_Arc_Form f;
        f.center = point(&in, dim);
        f.u = point(&in, dim);
        f.v = point(&in, dim);
        f.unit = (GSHGEOM_Angle_Unit)(next(&in) % 5);
        f.t0 = angle(&in);
        f.sweep = angle(&in);
        f.flags = next(&in) % 3; /* 2 is no flag the record knows */
        const int has_end = next(&in) & 1;
        GSHGEOM_Point e = point(&in, dim);
        r = gshgeom_path_arc_conjugate_to(path, &f, has_end ? &e : NULL);
        if (r == GSHGEOM_OK) {
          arc_added(path, before);
        }
        break;
      }
      case OP_ARC_TANGENT: {
        GSHGEOM_Point p1 = point(&in, dim);
        GSHGEOM_Point p2 = point(&in, dim);
        r = gshgeom_path_arc_tangent_to(path, p1, p2, radius(&in));
        if (r == GSHGEOM_OK) {
          arc_added(path, before);
        }
        break;
      }
      case OP_ARC_BULGE: {
        static const double bulges[] = {0.0, 1.0, -1.0, 0.5, 1e-9, 1e9, 0.41421356237309503, -2.0, 1e300};
        GSHGEOM_Point q = point(&in, dim);
        const uint8_t pick = next(&in);
        volatile double zero = 0.0;
        const double bulge = pick < 9 ? bulges[pick] : (pick == 9 ? zero / zero : (pick == 10 ? 1.0 / zero : 0.25 * (double)(pick - 128)));
        r = gshgeom_path_arc_bulge_to(path, q, bulge);
        if (r == GSHGEOM_OK) {
          arc_added(path, before);
        }
        break;
      }
      case OP_ARC_RAW: {
        /* A record built from a centre form, begun by a move to its start: its end is the form's to the bit, a nudge
         * inside the agreement bound, or a nudge past it (which must be refused). */
        GSHGEOM_Arc_Form f;
        f.center = point(&in, dim);
        f.u = point(&in, dim);
        f.v = point(&in, dim);
        f.unit = (GSHGEOM_Angle_Unit)(next(&in) % 5);
        f.t0 = angle(&in);
        f.sweep = angle(&in);
        f.flags = next(&in) % 3;
        const uint8_t nudge = next(&in) % 3;
        GSHGEOM_Arc a;
        GSHGEOM_Point s, e;
        if (gshgeom_arc_point(&f, f.t0, &s) != GSHGEOM_OK || gshgeom_arc_point(&f, f.t0 + f.sweep, &e) != GSHGEOM_OK ||
            !finite(s.x) || !finite(s.y) || !finite(e.x) || !finite(e.y)) {
          continue;
        }
        const double scale = absd(f.u.x) + absd(f.u.y) + absd(f.v.x) + absd(f.v.y) + absd(f.center.x) + absd(f.center.y) +
            absd(s.x) + absd(s.y) + absd(e.x) + absd(e.y);
        if (nudge == 1) {
          e.x += scale * 0x1p-30; /* far past the agreement bound, 2^-40 of the scale */
        } else if (nudge == 2) {
          e.x += scale * 0x1p-48;
        }
        r = gshgeom_path_move_to(path, s);
        if (r != GSHGEOM_OK) {
          break;
        }
        const size_t b2 = gshgeom_path_element_count(path);
        const uint64_t h2 = gshgeom_path_hash(path);
        a.start = s;
        a.end = e;
        a.center = f.center;
        a.u = f.u;
        a.v = f.v;
        a.unit = f.unit;
        a.t0 = f.t0;
        a.sweep = f.sweep;
        a.flags = f.flags;
        r = gshgeom_path_append_arc(path, &a);
        CHECK(r == GSHGEOM_OK || r == GSHGEOM_ERR_INVALID || r == GSHGEOM_ERR_LIMIT || r == GSHGEOM_ERR_OOM);
        if (r != GSHGEOM_OK) {
          CHECK(gshgeom_path_element_count(path) == b2 && gshgeom_path_hash(path) == h2);
          outcomes[OUT_ARC_DISAGREE] += nudge == 1 && r == GSHGEOM_ERR_INVALID;
        } else {
          CHECK(nudge != 1 && dim == 2);
          arc_added(path, b2);
        }
        outcomes[OUT_PATH_INVALID] += r == GSHGEOM_ERR_INVALID;
        outcomes[OUT_PATH_LIMIT] += r == GSHGEOM_ERR_LIMIT;
        outcomes[OUT_OOM] += r == GSHGEOM_ERR_OOM;
        CHECK(gshgeom_path_element_count(path) <= limits.max_elements);
        continue;
      }
      case OP_ARC_PASS:
        arc_passes(path, &alloc, &limits);
        continue;
      case OP_ARC_MAP:
        if (dim == 2) {
          arc_map(&in, path);
        }
        continue;
      case OP_SHAPE:
        r = shape_op(&in, path, dim, before);
        break;
      case OP_DEVICE:
        if (dim == 2) {
          device_passes(&in, path, &alloc, &limits, sink);
        }
        continue;
      case OP_MEASURE:
        measure(&in, path, &alloc, &limits, sink);
        continue;
      case OP_ORIENT:
        orient(&in, dim);
        continue;
      default: {
        GSHGEOM_Bounds box;
        CHECK(gshgeom_path_bounds(path, &box) == GSHGEOM_OK);
        CHECK(box.empty || (box.min.x <= box.max.x && box.min.y <= box.max.y && box.min.z <= box.max.z));
        CHECK(gshgeom_path_dump(path, sink) == GSHGEOM_OK);
        continue;
      }
    }
    /* A verb either adds what it says or is refused and changes nothing. */
    CHECK(r == GSHGEOM_OK || r == GSHGEOM_ERR_INVALID || r == GSHGEOM_ERR_LIMIT || r == GSHGEOM_ERR_OOM);
    outcomes[OUT_PATH_INVALID] += r == GSHGEOM_ERR_INVALID;
    outcomes[OUT_PATH_LIMIT] += r == GSHGEOM_ERR_LIMIT;
    outcomes[OUT_OOM] += r == GSHGEOM_ERR_OOM;
    if (r != GSHGEOM_OK) {
      CHECK(gshgeom_path_element_count(path) == before && gshgeom_path_hash(path) == hash_before);
    }
    CHECK(gshgeom_path_element_count(path) <= limits.max_elements);
  }
  gshgeom_path_destroy(path);
  fclose(sink);
  CHECK(counter.live == 0);
  return 0;
}
