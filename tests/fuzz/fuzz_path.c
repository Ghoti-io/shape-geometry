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
 * orientation predicates.  A coordinate is drawn from a table of the degenerate
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
 *   length is the end of the path.
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
  OUT_COUNT
};

static const char * const outcome_names[OUT_COUNT] = {"a non-finite point refused", "a path's element cap",
    "a device path", "a device coordinate out of the 26.6 range", "a flattened path", "a flattened path's point cap",
    "a tolerance below the minimum", "a measured path", "a measured path refused", "a failed allocation",
    "a grid-aligned map", "an invalid map refused", "a degenerate orientation (sign 0)",
    "a device path holding a curve"};

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
  const int16_t raw = (int16_t)(uint16_t)(((unsigned)next(in) << 8) | next(in));
  return (double)raw / 16.0;
}

static GSHGEOM_Point point(Input * in, int dim) {
  GSHGEOM_Point p = {coordinate(in), coordinate(in), 0.0};
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

static int32_t fixed(Input * in) {
  static const int32_t table[] = {0, 1, -1, INT32_MAX, INT32_MIN, 1 << 24, -(1 << 24), 64};
  const uint8_t pick = next(in);
  if (pick < 8) {
    return table[pick];
  }
  return (int32_t)(uint32_t)(((unsigned)next(in) << 8) | next(in)) - 32768;
}

/* --- the map ---------------------------------------------------------------------------------------- */

static GSHGEOM_Result build_map(Input * in, GSHGEOM_Map * map) {
  GSHGEOM_Map step;
  GSHGEOM_Result r = gshgeom_map_identity(2, map);
  const int steps = next(in) % 4;
  for (int i = 0; i < steps && r == GSHGEOM_OK; i++) {
    switch (next(in) % 6) {
      case 0:
        r = gshgeom_map_translate(2, coordinate(in), coordinate(in), 0.0, &step);
        break;
      case 1:
        r = gshgeom_map_scale(2, coordinate(in), coordinate(in), 1.0, &step);
        break;
      case 2:
        r = gshgeom_map_rotate((GSHGEOM_Angle_Unit)(next(in) % 5), coordinate(in), &step);
        break;
      case 3: {
        double m[4] = {coordinate(in), coordinate(in), coordinate(in), coordinate(in)};
        double t[2] = {coordinate(in), coordinate(in)};
        r = gshgeom_map_from_matrix(2, m, t, &step);
        break;
      }
      case 4:
        r = gshgeom_map_design_to_model(coordinate(in), coordinate(in), coordinate(in), &step);
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

/* --- the entry point ---------------------------------------------------------------------------------- */

enum { OP_MOVE, OP_LINE, OP_QUAD, OP_CUBIC, OP_RQUAD, OP_RCUBIC, OP_CLOSE, OP_MOVE_FIXED, OP_LINE_FIXED,
       OP_QUAD_FIXED, OP_CUBIC_FIXED, OP_DEVICE, OP_MEASURE, OP_ORIENT, OP_BOUNDS, OP_COUNT };

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
        r = gshgeom_path_move_to_fixed(path, fixed(&in), fixed(&in), next(&in) % 40);
        break;
      case OP_LINE_FIXED:
        r = gshgeom_path_line_to_fixed(path, fixed(&in), fixed(&in), next(&in) % 40);
        break;
      case OP_QUAD_FIXED: {
        int32_t v[4] = {fixed(&in), fixed(&in), fixed(&in), fixed(&in)};
        r = gshgeom_path_quadratic_to_fixed(path, v[0], v[1], v[2], v[3], next(&in) % 40);
        break;
      }
      case OP_CUBIC_FIXED: {
        int32_t v[6] = {fixed(&in), fixed(&in), fixed(&in), fixed(&in), fixed(&in), fixed(&in)};
        r = gshgeom_path_cubic_to_fixed(path, v[0], v[1], v[2], v[3], v[4], v[5], next(&in) % 40);
        break;
      }
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
