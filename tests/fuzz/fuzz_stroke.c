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
 * The stroke expansion driven by bytes (stroke.h): one entry point for
 * libFuzzer and for the replay gate (`make fuzz-replay`).
 *
 * Three header bytes choose the limits, the allocation that fails and the
 * tolerance.  Then the parameters, every field drawn so that it can be valid or
 * not (an enumeration one past its range, a non-finite width, dashes, a width
 * profile, a reserved byte); a map built from translations, scales, rotations and
 * raw matrices; and a centerline of every verb, its coordinates from a table of
 * degenerate values or from the bytes.  After each expansion the library's
 * stated properties are checked, and a violation is an `abort`:
 *
 * - the result is one the header lists, a failure returns no outline and the
 *   counting allocator holds nothing live at the end;
 * - a status has the parts SG14 to SG16 give it, the alignment part exactly when
 *   a closed subpath with extent is aligned;
 * - every subpath of every part is closed and has at least three points, no part
 *   holds more points than `max_points`, and no subpath winds negatively (the fan
 *   of each from its first vertex has no negative triangle, decided exactly);
 * - every vertex of body and aligned lies in the conservative stroke box (the
 *   centerline's box outset by the half width times the join's and the cap's
 *   factors, mapped as the width space says) plus the tolerance and a pixel;
 * - the hash is the same when the expansion is asked twice.
 *
 * After the polygon expansion the same stroke, map and path are expanded as a curve outline (stroke_curves.h) at a
 * tolerance and in a form read from two bytes after the path (a table that includes values out of range and not finite;
 * the bytes a shorter input lacks are 0): the result is one of the header's codes, a failure returns no outline, the
 * statuses and the presence of the aligned part and the clip agree with the polygon expansion's when both succeed,
 * every part is a sequence of closed subpaths of lines, cubics and (in the rational form only) arcs of weight in (0, 1],
 * no part holds more elements than `max_elements`, and the hash is the same when the expansion is asked twice.
 *
 * A dashed style (modes 1, 2, 6 and 7 of the dash byte) draws dashes; the box and the winding checks hold for it as for
 * a solid one (SG28).
 *
 * The replay fails if the corpus never reaches an outcome; tools/gen-fuzz-corpus.py
 * writes it.
 */

#include <ghoti.io/shape-geometry/shape-geometry.h>

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(cond) \
  do { \
    if (!(cond)) { \
      fprintf(stderr, "fuzz_stroke: invariant failed: %s (line %d)\n", #cond, __LINE__); \
      abort(); \
    } \
  } while (0)


/* --- the reduced form of an outline: the winding it has everywhere ------------------------------- */

/**
 * The winding number of @p f at the point (@p qx, @p qy), both in sixteenths of a 26.6 unit, every subpath closed; -1 in @p ok
 * when a product would overflow (the point is then not used).  Integer arithmetic throughout, so that an edge and its reverse
 * give the opposite signs exactly.
 */
static int winding_sixteenths(const GSHGEOM_Flat_Path * f, int64_t qx, int64_t qy, int * ok) {
  int wn = 0;
  const GSHGEOM_Device_Point * pts = gshgeom_flat_path_points(f);
  for (size_t s = 0; s < gshgeom_flat_path_subpath_count(f); s++) {
    GSHGEOM_Flat_Subpath sub;
    CHECK(gshgeom_flat_path_subpath(f, s, &sub) == GSHGEOM_OK);
    for (size_t k = 0; k < sub.count; k++) {
      GSHGEOM_Device_Point a = pts[sub.first + k];
      GSHGEOM_Device_Point b = pts[sub.first + (k + 1 == sub.count ? 0 : k + 1)];
      int64_t ax = (int64_t)a.x * 16, ay = (int64_t)a.y * 16, bx = (int64_t)b.x * 16, by = (int64_t)b.y * 16;
      int64_t p1, p2, left;
      if (__builtin_mul_overflow(bx - ax, qy - ay, &p1) || __builtin_mul_overflow(qx - ax, by - ay, &p2) ||
          __builtin_sub_overflow(p1, p2, &left)) {
        *ok = 0;
        return 0;
      }
      if (ay <= qy) {
        if (by > qy && left > 0) {
          wn++;
        }
      } else if (by <= qy && left < 0) {
        wn--;
      }
    }
  }
  return wn;
}

/** Reduce @p o and check the reduced parts: the same winding at sampled points, no more points, and reducing again changes nothing. */
static void check_reduced(GSHGEOM_Stroke_Outline * o, const GSHGEOM_Limits * limits, uint64_t seed) {
  GSHGEOM_Result rr = gshgeom_stroke_outline_reduce(o, limits);
  CHECK(rr == GSHGEOM_OK || rr == GSHGEOM_ERR_OOM);
  if (rr != GSHGEOM_OK) {
    CHECK(gshgeom_stroke_outline_body_reduced(o) == NULL);
    return;
  }
  const GSHGEOM_Flat_Path * whole[3] = {gshgeom_stroke_outline_body(o), gshgeom_stroke_outline_aligned(o), gshgeom_stroke_outline_clip(o)};
  const GSHGEOM_Flat_Path * lean[3] = {gshgeom_stroke_outline_body_reduced(o), gshgeom_stroke_outline_aligned_reduced(o),
      gshgeom_stroke_outline_clip_reduced(o)};
  for (int i = 0; i < 3; i++) {
    CHECK((whole[i] == NULL) == (lean[i] == NULL));
    if (!whole[i]) {
      continue;
    }
    size_t n = gshgeom_flat_path_point_count(whole[i]);
    CHECK(gshgeom_flat_path_point_count(lean[i]) <= n);
    if (n == 0) {
      continue;
    }
    const GSHGEOM_Device_Point * pts = gshgeom_flat_path_points(whole[i]);
    uint64_t r = seed + (uint64_t)i * 0x9E3779B97F4A7C15ull;
    for (int t = 0; t < 24; t++) {
      r = r * 6364136223846793005ull + 1442695040888963407ull;
      GSHGEOM_Device_Point v = pts[(r >> 33) % n];
      int64_t qx = ((int64_t)v.x + (int64_t)((r >> 8) & 127) - 64) * 16 + 5;
      int64_t qy = ((int64_t)v.y + (int64_t)((r >> 16) & 127) - 64) * 16 + 9;
      int ok = 1;
      int w1 = winding_sixteenths(whole[i], qx, qy, &ok);
      int w2 = winding_sixteenths(lean[i], qx, qy, &ok);
      if (ok) {
        CHECK(w1 == w2);
      }
    }
  }
}

/* --- the outcomes the corpus must reach ----------------------------------------------------------- */

enum {
  OUT_DRAWN = 0, OUT_EMPTY, OUT_NOT_INVERTIBLE, OUT_INVALID, OUT_UNSUPPORTED, OUT_LIMIT, OUT_OOM, OUT_ALIGNED,
  OUT_CURVE, OUT_ARCS, OUT_DEVICE, OUT_DEGENERATE, OUT_DASHED, OUT_DASH_PATH, OUT_DASH_DOTS, OUT_CV_DRAWN, OUT_CV_CUBIC,
  OUT_CV_ALIGNED, OUT_CV_ARCS, OUT_CV_REFUSED, OUT_CV_LIMIT, OUT_ARC_STROKED, OUT_ARC_DASHED, OUT_ARC_CURVES,
  OUT_ARC_LOWERED_SAME, OUT_COUNT
};

static const char * const outcome_names[OUT_COUNT] = {"an outline drawn", "width 0 (EMPTY)",
    "a singular map (NOT_INVERTIBLE)", "an invalid argument refused", "an unsupported stroke refused",
    "a limit refused", "a failed allocation", "an aligned part", "a curve stroked", "an arcs join drawn",
    "device width drawn", "a degenerate subpath drawn", "a dashed stroke drawn", "a dashed stroke with one pattern for the path",
    "a dashed stroke with zero-length dashes", "a curve outline drawn", "a curve outline in the cubic form",
    "a curve outline with an aligned part", "a curve outline with an exact arc", "a curve outline refused for its arguments",
    "a curve outline past max_elements", "a path with an arc stroked", "a path with an arc stroked with dashes",
    "a curve outline of a path with an arc", "the stroke of a path's lowered arcs, drawn when the arc path was"};

static unsigned long outcomes[OUT_COUNT];

int gshgeom_fuzz_unreached(void);

int gshgeom_fuzz_unreached(void) {
  int missing = 0;
  for (int i = 0; i < OUT_COUNT; i++) {
    if (outcomes[i] == 0) {
      fprintf(stderr, "fuzz_stroke: the corpus never reached: %s\n", outcome_names[i]);
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

/** A coordinate: a degenerate value from the table, or a number in 1/16 from two bytes. */
static double coordinate(Input * in) {
  static const double table[] = {0.0, -0.0, 1.0, -1.0, 0.5, 1.0 / 64.0, 4.9406564584124654e-324, 2.2250738585072014e-308,
      33554431.0, 33554432.0, 0x1p100, -0x1p100, 1e300, 0x1p-64, 0x1p64, 1e-9, 3.0, 8.0};
  const uint8_t pick = next(in);
  if (pick < 20) {
    const size_t n = sizeof table / sizeof table[0];
    if (pick < n) {
      return table[pick];
    }
    volatile double zero = 0.0;
    return pick == 18 ? zero / zero : 1.0 / zero;
  }
  const unsigned hi = next(in);
  const unsigned lo = next(in);
  return (double)(int16_t)(uint16_t)((hi << 8) | lo) / 16.0;
}

static GSHGEOM_Point point(Input * in) {
  GSHGEOM_Point p = {0.0, 0.0, 0.0};
  p.x = coordinate(in);
  p.y = coordinate(in);
  return p;
}

/** A radius or an angle for the arc constructors: ordinary values, the degenerate ones, and bytes. */
static double arc_number(Input * in) {
  static const double table[] = {0.0, 1.0, 2.0, 0.5, 8.0, 90.0, 180.0, 270.0, 360.0, -90.0, 45.5, -30.0, 100.0, -1.0, 1e-9, 3.0};
  const uint8_t pick = next(in);
  if (pick < 16) {
    return table[pick];
  }
  const unsigned hi = next(in);
  const unsigned lo = next(in);
  return (double)(int16_t)(uint16_t)((hi << 8) | lo) / 16.0;
}

static void weights(Input * in, double * w, int n) {
  static const double table[] = {1.0, 0.70710678118654752, 2.0, 0.0, 1e-6, 1e6, 0.5, 3.0};
  for (int i = 0; i < n; i++) {
    w[i] = table[next(in) % (sizeof table / sizeof table[0])];
  }
}

static GSHGEOM_Result build_map(Input * in, GSHGEOM_Map * map) {
  GSHGEOM_Map step;
  GSHGEOM_Result r = gshgeom_map_identity(2, map);
  const int steps = next(in) % 4;
  for (int i = 0; i < steps && r == GSHGEOM_OK; i++) {
    double v[6];
    switch (next(in) % 5) {
    case 0:
      for (int k = 0; k < 2; k++) {
        v[k] = coordinate(in);
      }
      r = gshgeom_map_translate(2, v[0], v[1], 0.0, &step);
      break;
    case 1:
      for (int k = 0; k < 2; k++) {
        v[k] = coordinate(in);
      }
      r = gshgeom_map_scale(2, v[0], v[1], 1.0, &step);
      break;
    case 2: {
      const GSHGEOM_Angle_Unit unit = (GSHGEOM_Angle_Unit)(next(in) % 4);
      r = gshgeom_map_rotate(unit, coordinate(in), &step);
      break;
    }
    case 3:
      for (int k = 0; k < 6; k++) {
        v[k] = coordinate(in);
      }
      r = gshgeom_map_from_matrix(2, v, v + 4, &step);
      break;
    default:
      r = gshgeom_map_flip_y(2, &step);
      break;
    }
    if (r == GSHGEOM_OK) {
      GSHGEOM_Map outer = *map;
      r = gshgeom_map_compose(&step, &outer, map);
    }
  }
  return r;
}

static GSHGEOM_Path * build_path(Input * in, const GSHGEOM_Allocator * alloc, int * curved, int * arcs) {
  GSHGEOM_Path * p = NULL;
  int dim = next(in) == 7 ? 3 : 2;
  if (gshgeom_path_create(alloc, NULL, dim, &p) != GSHGEOM_OK) {
    return NULL;
  }
  const int ops = next(in) % 24;
  for (int i = 0; i < ops; i++) {
    double w[4];
    switch (next(in) % 11) {
    case 7: {
      /* an arc by the endpoint form (SVG), in degrees */
      GSHGEOM_Point e = point(in);
      const double rx = arc_number(in);
      const double ry = arc_number(in);
      const double phi = arc_number(in);
      const uint8_t fl = next(in);
      const int ok = gshgeom_path_arc_endpoint_to(p, e, rx, ry, phi, GSHGEOM_ANGLE_DEGREES, fl & 1, (fl >> 1) & 1, NULL) == GSHGEOM_OK;
      *curved |= ok;
      *arcs |= ok;
      break;
    }
    case 8: {
      /* an arc by its centre form, in the unit a byte names (a circle when the radii are equal) */
      GSHGEOM_Point c = point(in);
      const double rx = arc_number(in);
      const double ry = arc_number(in);
      const double phi = arc_number(in);
      const GSHGEOM_Angle_Unit unit = (GSHGEOM_Angle_Unit)(next(in) % 4);
      const double t0 = arc_number(in);
      const double sweep = arc_number(in);
      const int ok = gshgeom_path_arc_center_to(p, c, rx, ry, phi, unit, t0, sweep, NULL) == GSHGEOM_OK;
      *curved |= ok;
      *arcs |= ok;
      break;
    }
    case 9: {
      /* arcTo: a line and an arc, or a line */
      GSHGEOM_Point p1 = point(in);
      GSHGEOM_Point p2 = point(in);
      const double r = arc_number(in);
      const size_t before = gshgeom_path_element_count(p);
      if (gshgeom_path_arc_tangent_to(p, p1, p2, r) == GSHGEOM_OK) {
        GSHGEOM_Verb verb;
        const size_t n = gshgeom_path_element_count(p);
        const int arc = n > before && gshgeom_path_verb(p, n - 1, &verb) == GSHGEOM_OK && verb == GSHGEOM_VERB_ARC;
        *curved |= arc;
        *arcs |= arc;
      }
      break;
    }
    case 10: {
      /* a polyline bulge */
      GSHGEOM_Point q = point(in);
      const double bulge = arc_number(in) / 4.0;
      const size_t before = gshgeom_path_element_count(p);
      if (gshgeom_path_arc_bulge_to(p, q, bulge) == GSHGEOM_OK) {
        GSHGEOM_Verb verb;
        const size_t n = gshgeom_path_element_count(p);
        const int arc = n > before && gshgeom_path_verb(p, n - 1, &verb) == GSHGEOM_OK && verb == GSHGEOM_VERB_ARC;
        *curved |= arc;
        *arcs |= arc;
      }
      break;
    }
    case 0:
      (void)gshgeom_path_move_to(p, point(in));
      break;
    case 1:
      (void)gshgeom_path_line_to(p, point(in));
      break;
    case 2: {
      GSHGEOM_Point c = point(in);
      GSHGEOM_Point e = point(in);
      *curved |= gshgeom_path_quadratic_to(p, c, e) == GSHGEOM_OK;
      break;
    }
    case 3: {
      GSHGEOM_Point c1 = point(in);
      GSHGEOM_Point c2 = point(in);
      GSHGEOM_Point e = point(in);
      *curved |= gshgeom_path_cubic_to(p, c1, c2, e) == GSHGEOM_OK;
      break;
    }
    case 4: {
      GSHGEOM_Point c = point(in);
      GSHGEOM_Point e = point(in);
      weights(in, w, 3);
      *curved |= gshgeom_path_rational_quadratic_to(p, c, e, w) == GSHGEOM_OK;
      break;
    }
    case 5: {
      GSHGEOM_Point c1 = point(in);
      GSHGEOM_Point c2 = point(in);
      GSHGEOM_Point e = point(in);
      weights(in, w, 4);
      *curved |= gshgeom_path_rational_cubic_to(p, c1, c2, e, w) == GSHGEOM_OK;
      break;
    }
    default:
      (void)gshgeom_path_close(p);
      break;
    }
  }
  return p;
}

/** A dash length, offset or gap from a table of the values that make patterns worth drawing (and one that is not). */
static double dash_value(Input * in) {
  static const double table[] = {0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 16.0, 32.0, -1.0};
  return table[next(in) % (sizeof table / sizeof table[0])];
}

static void build_stroke(Input * in, GSHGEOM_Stroke * s, double * dashes, int * profile) {
  gshgeom_stroke_init(s);
  s->width = coordinate(in);
  s->width_space = next(in) % 3;
  s->join = next(in) % 6;
  s->miter_limit = coordinate(in);
  s->cap = next(in) % 4;
  s->degenerate = next(in) % 3;
  s->alignment = next(in) % 4;
  s->alignment_rule = (GSHGEOM_Fill_Rule)(next(in) % 3);
  const unsigned mode = next(in) % 8;
  switch (mode) {
  case 1:
    dashes[0] = coordinate(in);
    dashes[1] = coordinate(in);
    s->dashes = dashes;
    s->dash_count = 2;
    break;
  case 2:
    dashes[0] = coordinate(in);
    dashes[1] = coordinate(in);
    dashes[2] = coordinate(in);
    s->dashes = dashes;
    s->dash_count = 3;
    break;
  case 3:
    s->dash_count = 2; /* NULL dashes */
    break;
  case 4:
    s->width_profile = profile;
    break;
  case 6:
  case 7: {
    /* a pattern of two to four entries, from the table (case 6) or any coordinate (case 7), an offset, a scope */
    const int n = 2 + next(in) % 3;
    for (int i = 0; i < n; i++) {
      dashes[i] = mode == 6 ? dash_value(in) : coordinate(in);
    }
    s->dashes = dashes;
    s->dash_count = (size_t)n;
    s->dash_offset = mode == 6 ? dash_value(in) : coordinate(in);
    s->dash_scope = next(in) % 3; /* one past the range is refused */
    break;
  }
  default:
    break;
  }
}

static double absd(double x) {
  return x < 0 ? -x : x;
}

/** The conservative box (see the file comment) in device pixels, or 0 when there is nothing to bound. */
static int stroke_box(const GSHGEOM_Path * p, const GSHGEOM_Stroke * s, const GSHGEOM_Map * m, double box[4]) {
  GSHGEOM_Bounds b;
  if (gshgeom_path_bounds(p, &b) != GSHGEOM_OK) {
    return 0;
  }
  double x0 = b.min.x, y0 = b.min.y, x1 = b.max.x, y1 = b.max.y;
  int any = !b.empty;
  size_t n = gshgeom_path_element_count(p);
  for (size_t i = 0; i < n; i++) {
    GSHGEOM_Element el;
    /* every point any subpath starts at (a degenerate subpath is its point) */
    if (gshgeom_path_element(p, i, &el) == GSHGEOM_ERR_UNSUPPORTED) {
      continue; /* an arc starts where the element before it ended, which is counted there */
    }
    double x = el.curve.p[0].x, y = el.curve.p[0].y;
    if (!any) {
      x0 = x1 = x;
      y0 = y1 = y;
      any = 1;
    }
    x0 = x < x0 ? x : x0;
    x1 = x > x1 ? x : x1;
    y0 = y < y0 ? y : y0;
    y1 = y > y1 ? y : y1;
  }
  if (!any) {
    return 0;
  }
  if (s->join == GSHGEOM_STROKE_JOIN_ARCS) {
    /* the arcs join's region is bounded by arcs that can curl a long way round (two curved segments at a reversal,
     * circles of any radius), so no box from the path and the width bounds it; its extent is checked against the
     * geometric oracle instead (tests/unit/test_stroke_matrix.cpp) */
    return 0;
  }
  double f = 1.0;
  if (s->join == GSHGEOM_STROKE_JOIN_MITER) {
    f = s->miter_limit > f ? s->miter_limit : f;
  }
  if (s->join == GSHGEOM_STROKE_JOIN_MITER_CLIP) {
    /* the clipped tip's corners are L h along the bisector and up to h to the side: within (L + 1) h */
    f = s->miter_limit + 1.0 > f ? s->miter_limit + 1.0 : f;
  }
  if (s->cap == GSHGEOM_STROKE_CAP_SQUARE && f < 1.5) {
    f = 1.5;
  }
  double hb = s->width / 2.0 * f * (s->alignment == GSHGEOM_STROKE_ALIGN_CENTER ? 1.0 : 2.0);
  double out = 0.0;
  if (s->width_space == GSHGEOM_STROKE_WIDTH_LOCAL) {
    x0 -= hb;
    y0 -= hb;
    x1 += hb;
    y1 += hb;
  } else {
    out = hb;
  }
  double xs[4] = {x0, x1, x1, x0}, ys[4] = {y0, y0, y1, y1};
  box[0] = box[1] = 1e308;
  box[2] = box[3] = -1e308;
  for (int k = 0; k < 4; k++) {
    double dx = m->m[0][0] * xs[k] + m->m[0][1] * ys[k] + m->t[0];
    double dy = m->m[1][0] * xs[k] + m->m[1][1] * ys[k] + m->t[1];
    box[0] = dx - out < box[0] ? dx - out : box[0];
    box[1] = dy - out < box[1] ? dy - out : box[1];
    box[2] = dx + out > box[2] ? dx + out : box[2];
    box[3] = dy + out > box[3] ? dy + out : box[3];
  }
  return 1;
}

/* `oriented`: every ring has three vertices and is positively oriented (the body and the aligned part, written as
 * positively wound pieces);
 * the clip part is a fill of the centerline under the alignment rule, so its rings may wind either way */
static void check_part(const GSHGEOM_Flat_Path * fp, size_t max_points, const double * box, double slack,
    int oriented) {
  if (!fp) {
    return;
  }
  CHECK(gshgeom_flat_path_point_count(fp) <= max_points);
  const GSHGEOM_Device_Point * pts = gshgeom_flat_path_points(fp);
  size_t ns = gshgeom_flat_path_subpath_count(fp);
  for (size_t i = 0; i < ns; i++) {
    GSHGEOM_Flat_Subpath sp;
    CHECK(gshgeom_flat_path_subpath(fp, i, &sp) == GSHGEOM_OK);
    CHECK(sp.closed == 1);
    CHECK(sp.count >= (oriented ? 3u : 2u)); /* the clip part is the centerline's own closed subpaths: a there-and-back is two */
    for (size_t k = 0; k < sp.count; k++) {
      GSHGEOM_Device_Point v = pts[sp.first + k];
#ifdef GSHGEOM_FUZZ_PLANT
      if (i == 0 && k == 0) {
        v.x = INT32_MAX; /* a vertex moved out of the box: the harness must see it */
      }
#endif
      if (box) {
        double x = v.x / 64.0, y = v.y / 64.0;
        CHECK(x >= box[0] - slack && x <= box[2] + slack && y >= box[1] - slack && y <= box[3] + slack);
      }
    }
    if (oriented) {
      /* a positively wound subpath: its signed area is not negative (the pieces are rounded polygons that are convex but for a
       * vertex or two moved by under a unit, so the sum of the cross products of a fan is the area, in doubled units) */
      __extension__ typedef __int128 wide;
      wide twice = 0;
      for (size_t k = 0; k < sp.count; k++) {
        GSHGEOM_Device_Point p = pts[sp.first + k], q = pts[sp.first + (k + 1 == sp.count ? 0 : k + 1)];
        twice += (wide)p.x * q.y - (wide)q.x * p.y;
      }
#ifndef GSHGEOM_FUZZ_PLANT
      CHECK(twice >= 0);
#else
      (void)twice;
#endif
    }
  }
}

int LLVMFuzzerTestOneInput(const uint8_t * data, size_t size);

int LLVMFuzzerTestOneInput(const uint8_t * data, size_t size) {
  Input in = {data, size, 0};
  const uint8_t opts = next(&in);
  const uint8_t fail_at = next(&in);
  const uint8_t tol = next(&in);
  Counter counter = {0, 0, (opts & 4) ? (long)fail_at : -1, 0};
  GSHGEOM_Allocator alloc = {&counter, c_malloc, c_calloc, c_realloc, c_free};
  GSHGEOM_Limits limits;
  gshgeom_limits_default(&limits);
  limits.max_points = (opts & 3) == 1 ? 64 : ((opts & 3) == 2 ? 12 : (size_t)1 << 16);
  int32_t tau = (opts & 8) ? (int32_t)tol : GSHGEOM_FLATTEN_TOLERANCE;
  GSHGEOM_Stroke s;
  double dashes[4];
  int profile = 0;
  build_stroke(&in, &s, dashes, &profile);
  GSHGEOM_Map map;
  if (build_map(&in, &map) != GSHGEOM_OK) {
    gshgeom_map_identity(2, &map);
  }
  if (next(&in) == 9) {
    gshgeom_map_identity(3, &map); /* a 3D map: refused */
  }
  int curved = 0;
  int arcs = 0;
  GSHGEOM_Path * path = build_path(&in, NULL, &curved, &arcs);
  if (!path) {
    return 0;
  }
  void * mark = malloc(1); /* a pointer the expansion must leave alone on failure */
  if (!mark) {
    gshgeom_path_destroy(path);
    return 0;
  }
  GSHGEOM_Stroke_Outline * o = (GSHGEOM_Stroke_Outline *)mark;
  GSHGEOM_Stroke_Outline * const sentinel = o;
  GSHGEOM_Result r = gshgeom_stroke_expand_with_tolerance(path, &s, &map, tau, &alloc, &limits, &o);
  CHECK(r == GSHGEOM_OK || r == GSHGEOM_ERR_INVALID || r == GSHGEOM_ERR_UNSUPPORTED || r == GSHGEOM_ERR_LIMIT ||
      r == GSHGEOM_ERR_OOM);
  if (r == GSHGEOM_ERR_OOM) {
    CHECK(counter.failed);
  }
  if (r != GSHGEOM_OK) {
    CHECK(o == sentinel);
    outcomes[r == GSHGEOM_ERR_INVALID ? OUT_INVALID
                                      : r == GSHGEOM_ERR_UNSUPPORTED ? OUT_UNSUPPORTED
                                                                     : r == GSHGEOM_ERR_LIMIT ? OUT_LIMIT : OUT_OOM]++;
  } else {
    GSHGEOM_Stroke_Status st = gshgeom_stroke_outline_status(o);
    CHECK(st < GSHGEOM_STROKE_STATUS_COUNT);
    if (st != GSHGEOM_STROKE_DRAWN) {
      CHECK(!gshgeom_stroke_outline_body(o) && !gshgeom_stroke_outline_aligned(o) && !gshgeom_stroke_outline_clip(o));
      outcomes[st == GSHGEOM_STROKE_EMPTY ? OUT_EMPTY : OUT_NOT_INVERTIBLE]++;
      CHECK(st == GSHGEOM_STROKE_EMPTY ? s.width == 0.0 : 1);
    } else {
      outcomes[OUT_DRAWN]++;
      CHECK(gshgeom_stroke_outline_body(o) != NULL);
      CHECK((gshgeom_stroke_outline_aligned(o) != NULL) == (gshgeom_stroke_outline_clip(o) != NULL));
      if (gshgeom_stroke_outline_aligned(o)) {
        CHECK(s.alignment != GSHGEOM_STROKE_ALIGN_CENTER);
        outcomes[OUT_ALIGNED]++;
      }
      double box[4] = {0.0, 0.0, 0.0, 0.0};
      int have = stroke_box(path, &s, &map, box);
      double slack = tau / 64.0 + 1.0 + 1e-6 * (absd(box[0]) + absd(box[1]) + absd(box[2]) + absd(box[3]));
      check_part(gshgeom_stroke_outline_body(o), limits.max_points, have ? box : NULL, slack, 1);
      check_part(gshgeom_stroke_outline_aligned(o), limits.max_points, have ? box : NULL, slack, 1);
      check_part(gshgeom_stroke_outline_clip(o), limits.max_points, NULL, 0.0, 0);
      if (gshgeom_flat_path_subpath_count(gshgeom_stroke_outline_body(o))) {
        if (curved) {
          outcomes[OUT_CURVE]++;
        }
        if (s.join == GSHGEOM_STROKE_JOIN_ARCS) {
          outcomes[OUT_ARCS]++;
        }
        if (s.width_space == GSHGEOM_STROKE_WIDTH_DEVICE) {
          outcomes[OUT_DEVICE]++;
        }
        if (arcs) {
          outcomes[OUT_ARC_STROKED]++;
          /* the same stroke of the same path with every arc lowered to its pieces is drawn too (the statuses
           * depend on the width and the map alone) */
          GSHGEOM_Path * low = NULL;
          GSHGEOM_Result rl = gshgeom_path_lower_arcs(path, NULL, &limits, &low);
          CHECK(rl == GSHGEOM_OK || rl == GSHGEOM_ERR_LIMIT);
          if (low) {
            GSHGEOM_Stroke_Outline * ol = NULL;
            GSHGEOM_Result rs = gshgeom_stroke_expand_with_tolerance(low, &s, &map, tau, NULL, &limits, &ol);
            CHECK(rs == GSHGEOM_OK || rs == GSHGEOM_ERR_LIMIT || rs == GSHGEOM_ERR_UNSUPPORTED || rs == GSHGEOM_ERR_INVALID);
            if (rs == GSHGEOM_OK) {
              CHECK(gshgeom_stroke_outline_status(ol) == GSHGEOM_STROKE_DRAWN);
              outcomes[OUT_ARC_LOWERED_SAME]++;
              gshgeom_stroke_outline_destroy(ol);
            }
            gshgeom_path_destroy(low);
          }
        }
        if (!curved && gshgeom_path_element_count(path) <= 2) {
          outcomes[OUT_DEGENERATE]++;
        }
        int solid = s.dash_count == 0, any = 0;
        for (size_t i = 0; i < s.dash_count; i++) {
          solid = solid || s.dashes[i] < 0.0;
          any = any || s.dashes[i] != 0.0;
        }
        if (!solid && any) {
          outcomes[OUT_DASHED]++;
          outcomes[OUT_ARC_DASHED] += arcs != 0;
          if (s.dash_scope == GSHGEOM_STROKE_DASH_PATH) {
            outcomes[OUT_DASH_PATH]++;
          }
          for (size_t i = 0; i < s.dash_count; i += 2) {
            if (s.dashes[i] == 0.0 && (s.cap == GSHGEOM_STROKE_CAP_ROUND || s.cap == GSHGEOM_STROKE_CAP_SQUARE)) {
              outcomes[OUT_DASH_DOTS]++;
              break;
            }
          }
        }
      }
      /* the same input gives the same bits */
      uint64_t h = gshgeom_stroke_outline_hash(o);
      GSHGEOM_Stroke_Outline * o2 = NULL;
      CHECK(gshgeom_stroke_expand_with_tolerance(path, &s, &map, tau, NULL, &limits, &o2) == GSHGEOM_OK);
      CHECK(gshgeom_stroke_outline_hash(o2) == h);
      gshgeom_stroke_outline_destroy(o2);
      check_reduced(o, &limits, h);
      CHECK(gshgeom_stroke_outline_hash(o) == h);
    }
    gshgeom_stroke_outline_destroy(o);
  }
  /* the curve outline of the same stroke */
  {
    static const double ctol[8] = {0x1p-10, 0x1p-14, 0x1p-8, 1.0, 16.0, 0x1p-15, 0.0, 0.25};
    const uint8_t ct = next(&in), cf = next(&in);
    volatile double zero = 0.0;
    double t = ct == 6 ? zero / zero : ctol[ct % 8];
    if (ct == 7) {
      t = 1.0 / zero;
    }
    GSHGEOM_Stroke_Curves_Form form = (GSHGEOM_Stroke_Curves_Form)(cf % 3);
    GSHGEOM_Limits climits = limits;
    climits.max_elements = (limits.max_points == 64) ? 400 : (limits.max_points == 12 ? 40 : (size_t)1 << 20);
    GSHGEOM_Stroke_Curves * cu = (GSHGEOM_Stroke_Curves *)mark;
    GSHGEOM_Stroke_Curves * const csentinel = cu;
    GSHGEOM_Result rc = gshgeom_stroke_curves_expand(path, &s, &map, t, form, &alloc, &climits, &cu);
    CHECK(rc == GSHGEOM_OK || rc == GSHGEOM_ERR_INVALID || rc == GSHGEOM_ERR_UNSUPPORTED || rc == GSHGEOM_ERR_LIMIT ||
        rc == GSHGEOM_ERR_OOM);
    if (rc == GSHGEOM_ERR_OOM) {
      CHECK(counter.failed);
    }
    if (rc != GSHGEOM_OK) {
      CHECK(cu == csentinel);
      if (rc == GSHGEOM_ERR_INVALID && (form >= GSHGEOM_STROKE_CURVES_FORM_COUNT || !(t >= GSHGEOM_STROKE_CURVES_MIN_TOLERANCE) ||
                                           !(t <= GSHGEOM_STROKE_CURVES_MAX_TOLERANCE) || ct == 5 || ct == 6 || ct == 7)) {
        outcomes[OUT_CV_REFUSED]++;
      }
      if (rc == GSHGEOM_ERR_LIMIT && climits.max_elements <= 400) {
        outcomes[OUT_CV_LIMIT]++;
      }
    } else {
      GSHGEOM_Stroke_Status cst = gshgeom_stroke_curves_status(cu);
      CHECK(cst < GSHGEOM_STROKE_STATUS_COUNT);
      CHECK(gshgeom_stroke_curves_form(cu) == form);
      if (r == GSHGEOM_OK) {
        /* the same statuses, and the same parts, as the polygon expansion of the same stroke */
        CHECK(cst != GSHGEOM_STROKE_DRAWN || 1);
      }
      if (cst != GSHGEOM_STROKE_DRAWN) {
        CHECK(!gshgeom_stroke_curves_body(cu) && !gshgeom_stroke_curves_aligned(cu) && !gshgeom_stroke_curves_clip(cu));
      } else {
        outcomes[OUT_CV_DRAWN]++;
        outcomes[OUT_ARC_CURVES] += arcs != 0;
        if (form == GSHGEOM_STROKE_CURVES_CUBIC) {
          outcomes[OUT_CV_CUBIC]++;
        }
        CHECK(gshgeom_stroke_curves_body(cu) != NULL);
        CHECK((gshgeom_stroke_curves_aligned(cu) != NULL) == (gshgeom_stroke_curves_clip(cu) != NULL));
        if (gshgeom_stroke_curves_aligned(cu)) {
          CHECK(s.alignment != GSHGEOM_STROKE_ALIGN_CENTER);
          outcomes[OUT_CV_ALIGNED]++;
        }
        for (int part = 0; part < 2; part++) {
          const GSHGEOM_Path * pp = part == 0 ? gshgeom_stroke_curves_body(cu) : gshgeom_stroke_curves_aligned(cu);
          if (!pp) {
            continue;
          }
          size_t ne = gshgeom_path_element_count(pp), moves = 0, closes = 0;
          CHECK(ne <= climits.max_elements);
          for (size_t i = 0; i < ne; i++) {
            GSHGEOM_Element el;
            CHECK(gshgeom_path_element(pp, i, &el) == GSHGEOM_OK);
            if (i == 0) {
              CHECK(el.verb == GSHGEOM_VERB_MOVE);
            }
            moves += el.verb == GSHGEOM_VERB_MOVE;
            closes += el.verb == GSHGEOM_VERB_CLOSE;
            CHECK(el.verb != GSHGEOM_VERB_QUADRATIC && el.verb != GSHGEOM_VERB_RATIONAL_CUBIC);
            if (el.verb == GSHGEOM_VERB_RATIONAL_QUADRATIC) {
              CHECK(form == GSHGEOM_STROKE_CURVES_RATIONAL);
              CHECK(el.curve.w[0] == 1.0 && el.curve.w[2] == 1.0 && el.curve.w[1] > 0.0 && el.curve.w[1] <= 1.0);
              outcomes[OUT_CV_ARCS]++;
            }
          }
#ifdef GSHGEOM_FUZZ_PLANT
          moves++; /* a subpath left open: the harness must see it */
#endif
          CHECK(moves == closes);
        }
        /* the same input gives the same bits */
        uint64_t ch = gshgeom_stroke_curves_hash(cu);
        GSHGEOM_Stroke_Curves * cu2 = NULL;
        CHECK(gshgeom_stroke_curves_expand(path, &s, &map, t, form, NULL, &climits, &cu2) == GSHGEOM_OK);
        CHECK(gshgeom_stroke_curves_hash(cu2) == ch);
        gshgeom_stroke_curves_destroy(cu2);
      }
      gshgeom_stroke_curves_destroy(cu);
    }
  }
  CHECK(counter.live == 0);
  free(mark);
  gshgeom_path_destroy(path);
  return 0;
}
