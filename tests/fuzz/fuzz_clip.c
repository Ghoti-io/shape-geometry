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
 * The clip of a flattened path to a convex quadrilateral (clip.h, CK1 to CK6) driven by bytes: one entry point for libFuzzer
 * and for the replay gate (`make fuzz-replay`).
 *
 * Two header bytes choose the range of the coordinates (small, where the winding can be checked in exact 64-bit integers,
 * or the whole int32 range), the limit (the default or a small `max_points`) and the allocation that fails.  Then come the
 * quadrilateral (a rectangle or a sheared one from a centre and half sizes, reversed, with a repeated point, a bowtie, or
 * four raw points that are mostly refused) and one to four contours of up to twelve points from a table of degenerate
 * values (zeros, the int32 limits) or from the bytes, with a closed flag.
 *
 * After each clip the library's stated properties are checked, and a violation is an `abort`:
 *
 * - the result code is one the header lists, and it is ::GSHGEOM_ERR_INVALID exactly when the four turns of the
 *   quadrilateral (::gshgeom_orient2d_exact) are not all of one sign or zero, or are all zero, or its area is zero;
 * - a refusal leaves nothing behind (the counting allocator holds nothing live at the end) and writes nothing;
 * - the input is untouched (its hash is the same after), and the same call twice gives the same hash (CK6);
 * - the result has no more subpaths than the input, at most `n + c` points for `n` input points and `c` crossings and no
 *   more than the limit (CK5, and `c` is at most nine times `n`: the proved bound is 8.125), every point lies in the box of the input's
 *   points, and a contour wholly inside the closed quadrilateral is in the result bit for bit, in order (CK2);
 * - in the small range, the winding number of the result equals the input's inside the interior of the quadrilateral and
 *   is zero outside it, at the probes at least three units from every edge of the input and of the quadrilateral, by an
 *   independent count of ray crossings with the sign of each crossing from ::gshgeom_orient2d_exact (CK2).
 *
 * The replay fails if the corpus never reaches an outcome (each refusal kind and each success): tools/gen-fuzz-corpus.py
 * writes it.
 */

#include <ghoti.io/shape-geometry/shape-geometry.h>

#include "../../src/clip/clip_internal.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(cond) \
  do { \
    if (!(cond)) { \
      fprintf(stderr, "fuzz_clip: invariant failed: %s (line %d)\n", #cond, __LINE__); \
      abort(); \
    } \
  } while (0)

/* --- the outcomes the corpus must reach ----------------------------------------------------------------------------- */

enum {
  OUT_INVALID_QUAD = 0, OUT_LIMIT, OUT_OOM, OUT_EMPTY, OUT_CUT, OUT_COPIED, OUT_COPIED_OPEN, OUT_WIDE, OUT_CLOCKWISE,
  OUT_COUNTERCLOCKWISE, OUT_REPEATED_POINT, OUT_ON_THE_LINE, OUT_WINDING_CHECKED, OUT_WINDING_NONZERO, OUT_COUNT
};

static const char * const outcome_names[OUT_COUNT] = {"a degenerate or non-convex quadrilateral refused", "the point limit",
    "a failed allocation", "an empty result", "a contour cut", "a contour copied bit for bit", "an open contour copied",
    "the whole int32 range", "a clockwise quadrilateral", "a counterclockwise quadrilateral", "a quadrilateral with a repeated point",
    "a point exactly on a side of the quadrilateral", "a winding compared", "a nonzero winding compared"};

static unsigned long outcomes[OUT_COUNT];

int gshgeom_fuzz_unreached(void);

int gshgeom_fuzz_unreached(void) {
  int missing = 0;
  for (int i = 0; i < OUT_COUNT; i++) {
    if (outcomes[i] == 0) {
      fprintf(stderr, "fuzz_clip: the corpus never reached: %s\n", outcome_names[i]);
      missing++;
    }
  }
  return missing;
}

/* --- a counting allocator that can fail one call -------------------------------------------------------------------- */

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

/* --- the input ------------------------------------------------------------------------------------------------------- */

typedef struct Input {
  const uint8_t * data;
  size_t size;
  size_t at;
} Input;

static uint8_t next(Input * in) {
  return in->at < in->size ? in->data[in->at++] : 0;
}

/* A coordinate: a degenerate value from the table, or two bytes (small: a 16-bit value over 4; wide: that shifted up 15). */
static int32_t coordinate(Input * in, int wide) {
  static const int32_t small_table[] = {0, 1, -1, 63, 64, 65, -64, 127, 128, 8191, -8192, 4096};
  static const int32_t wide_table[] = {INT32_MAX, INT32_MIN, INT32_MAX - 1, INT32_MIN + 1, 1 << 30, -(1 << 30), 0, 1, -1};
  const uint8_t pick = next(in);
  if (pick < 12) {
    return small_table[pick];
  }
  if (wide && pick < 21) {
    return wide_table[pick - 12];
  }
  const unsigned hi = next(in);
  const unsigned lo = next(in);
  const int32_t v = (int32_t)(int16_t)(uint16_t)((hi << 8) | lo);
  return wide ? (int32_t)((int64_t)v * 32768) : v / 4;
}

static GSHGEOM_Device_Point point(Input * in, int wide) {
  GSHGEOM_Device_Point p;
  p.x = coordinate(in, wide);
  p.y = coordinate(in, wide);
  return p;
}

static GSHGEOM_Point as_point(GSHGEOM_Device_Point p) {
  GSHGEOM_Point q = {(double)p.x, (double)p.y, 0.0};
  return q;
}

/* The sign of the side of p of the line a -> b, from the public exact predicate. */
static int orient(GSHGEOM_Device_Point a, GSHGEOM_Device_Point b, GSHGEOM_Device_Point p) {
  int s = 0;
  CHECK(gshgeom_orient2d_exact(as_point(a), as_point(b), as_point(p), &s) == GSHGEOM_OK);
  return s;
}

/* Whether the quadrilateral is valid under CK1, and its orientation (+1 or -1) when it is. */
static int quad_valid(const GSHGEOM_Device_Point q[4], int * sigma_out) {
  int turn[4];
  int sigma = 0;
  for (int i = 0; i < 4; i++) {
    turn[i] = orient(q[(i + 3) % 4], q[i], q[(i + 1) % 4]);
    if (turn[i]) {
      if (sigma && sigma != turn[i]) {
        return 0;
      }
      sigma = turn[i];
    }
  }
  if (sigma == 0 || (turn[1] == 0 && turn[3] == 0)) {
    return 0;
  }
  *sigma_out = sigma;
  return 1;
}

static int path_equal(const GSHGEOM_Flat_Path * a, const GSHGEOM_Flat_Path * b) {
  return gshgeom_flat_path_hash(a) == gshgeom_flat_path_hash(b);
}

/* --- the independent winding ---------------------------------------------------------------------------------------- */

/* The winding number of the closed subpaths of a flat path about (px, py): the crossings of the ray to the right, a downward
 * edge counting +1, the side of each from the exact predicate. */
static int winding(const GSHGEOM_Flat_Path * path, GSHGEOM_Device_Point p) {
  int w = 0;
  const GSHGEOM_Device_Point * pts = gshgeom_flat_path_points(path);
  for (size_t s = 0; s < gshgeom_flat_path_subpath_count(path); s++) {
    GSHGEOM_Flat_Subpath sub;
    CHECK(gshgeom_flat_path_subpath(path, s, &sub) == GSHGEOM_OK);
    for (size_t k = 0; k < sub.count; k++) {
      const GSHGEOM_Device_Point a = pts[sub.first + k];
      const GSHGEOM_Device_Point b = pts[sub.first + (k + 1 == sub.count ? 0 : k + 1)];
      if ((a.y <= p.y) == (b.y <= p.y)) {
        continue;
      }
      const int o = orient(a, b, p);
      /* o = dy * (x_edge - px): the edge is to the right of p when o has the sign of dy */
      if (b.y > a.y ? o > 0 : o < 0) {
        w += b.y > a.y ? 1 : -1;
      }
    }
  }
  return w;
}

/* Whether p is within 3 units of the segment a-b (int64 exact: the coordinates are within 2^14). */
static int near_segment(GSHGEOM_Device_Point a, GSHGEOM_Device_Point b, GSHGEOM_Device_Point p) {
  const int64_t ex = (int64_t)b.x - a.x, ey = (int64_t)b.y - a.y;
  const int64_t fx = (int64_t)p.x - a.x, fy = (int64_t)p.y - a.y;
  const int64_t len2 = ex * ex + ey * ey;
  const int64_t dot = fx * ex + fy * ey;
  if (len2 == 0 || dot <= 0) {
    return fx * fx + fy * fy <= 9;
  }
  if (dot >= len2) {
    const int64_t gx = (int64_t)p.x - b.x, gy = (int64_t)p.y - b.y;
    return gx * gx + gy * gy <= 9;
  }
  const int64_t cross = ex * fy - ey * fx;
  return cross * cross <= 9 * len2;
}

static int near_path(const GSHGEOM_Flat_Path * path, GSHGEOM_Device_Point p) {
  const GSHGEOM_Device_Point * pts = gshgeom_flat_path_points(path);
  for (size_t s = 0; s < gshgeom_flat_path_subpath_count(path); s++) {
    GSHGEOM_Flat_Subpath sub;
    CHECK(gshgeom_flat_path_subpath(path, s, &sub) == GSHGEOM_OK);
    for (size_t k = 0; k < sub.count; k++) {
      if (near_segment(pts[sub.first + k], pts[sub.first + (k + 1 == sub.count ? 0 : k + 1)], p)) {
        return 1;
      }
    }
  }
  return 0;
}

/* --- the entry point --------------------------------------------------------------------------------------------------- */

#define MAX_CONTOURS 4
#define MAX_POINTS 12

int LLVMFuzzerTestOneInput(const uint8_t * data, size_t size);

int LLVMFuzzerTestOneInput(const uint8_t * data, size_t size) {
  Input in = {data, size, 0};
  const uint8_t opts = next(&in);
  const uint8_t fail = next(&in);
  const int wide = (opts & 1) != 0;
  const int small_limit = (opts & 2) != 0;
  const int fail_alloc = (opts & 4) != 0;

  /* The quadrilateral. */
  GSHGEOM_Device_Point quad[4];
  const uint8_t qkind = next(&in) % 6u;
  if (qkind == 5) {
    for (int i = 0; i < 4; i++) {
      quad[i] = point(&in, wide);
    }
  } else {
    const int64_t cx = coordinate(&in, wide) / (wide ? 4 : 1);
    const int64_t cy = coordinate(&in, wide) / (wide ? 4 : 1);
    int64_t hx = coordinate(&in, wide) / (wide ? 4 : 1);
    int64_t hy = coordinate(&in, wide) / (wide ? 4 : 1);
    hx = (hx < 0 ? -hx : hx) / 2 + 1;
    hy = (hy < 0 ? -hy : hy) / 2 + 1;
    const int64_t shear = (int64_t)(int8_t)next(&in) * (wide ? 4096 : 1) / 8;
    const int64_t px[4] = {cx - hx, cx + hx, cx + hx + shear, cx - hx};
    const int64_t py[4] = {cy - hy, cy - hy + shear, cy + hy, cy + hy - shear};
    for (int i = 0; i < 4; i++) {
      quad[i].x = (int32_t)px[i];
      quad[i].y = (int32_t)py[i];
    }
    if (qkind == 1) {
      GSHGEOM_Device_Point t = quad[0];
      quad[0] = quad[3];
      quad[3] = t;
      t = quad[1];
      quad[1] = quad[2];
      quad[2] = t;
    } else if (qkind == 2) {
      quad[3] = quad[2];
    } else if (qkind == 3) {
      GSHGEOM_Device_Point t = quad[1];
      quad[1] = quad[2];
      quad[2] = t;
    }
  }

  /* The contours. */
  GSHGEOM_Device_Point pts[MAX_CONTOURS * MAX_POINTS];
  GSHGEOM_Flat_Subpath subs[MAX_CONTOURS];
  size_t np = 0;
  const size_t nc = 1u + next(&in) % MAX_CONTOURS;
  for (size_t c = 0; c < nc; c++) {
    const size_t n = 1u + next(&in) % MAX_POINTS;
    subs[c].first = np;
    subs[c].count = n;
    for (size_t k = 0; k < n; k++) {
      const uint8_t kind = next(&in) % 8u;
      if (kind == 0) {
        pts[np++] = quad[next(&in) % 4u];
      } else if (kind == 1) {
        const GSHGEOM_Device_Point a = quad[next(&in) % 4u];
        const GSHGEOM_Device_Point b = quad[next(&in) % 4u];
        const int64_t f = next(&in) % 9u;
        pts[np].x = (int32_t)(a.x + ((int64_t)b.x - a.x) * f / 8);
        pts[np].y = (int32_t)(a.y + ((int64_t)b.y - a.y) * f / 8);
        np++;
      } else {
        pts[np++] = point(&in, wide);
      }
    }
    subs[c].closed = next(&in) & 1;
  }

  Counter counter = {0, 0, fail_alloc ? (long)(fail % 24u) : -1, 0};
  GSHGEOM_Allocator alloc = {&counter, c_malloc, c_calloc, c_realloc, c_free};
  GSHGEOM_Limits limits;
  gshgeom_limits_default(&limits);
  if (small_limit) {
    limits.max_points = 1u + fail % 20u;
  }

  /* The input is built with the default allocator: the failing allocation belongs to the clip. */
  GSHGEOM_Flat_Path * flat = NULL;
  CHECK(gshgeom_flat_path_from_points(pts, np, subs, nc, NULL, NULL, &flat) == GSHGEOM_OK);
  const uint64_t input_hash = gshgeom_flat_path_hash(flat);

  int sigma = 0;
  const int valid = quad_valid(quad, &sigma);
  if (valid) {
    outcomes[sigma > 0 ? OUT_COUNTERCLOCKWISE : OUT_CLOCKWISE]++;
  }
  if (valid && (quad[0].x == quad[1].x && quad[0].y == quad[1].y)) {
    outcomes[OUT_REPEATED_POINT]++;
  }
  if (valid && quad[2].x == quad[3].x && quad[2].y == quad[3].y) {
    outcomes[OUT_REPEATED_POINT]++;
  }
  if (wide) {
    outcomes[OUT_WIDE]++;
  }

  GSHGEOM_Flat_Path * out = NULL;
  GSHGEOM_Clip_Stats stats = {0, 0, 0, 0};
  const GSHGEOM_Result r = gshgeom_clip_convex_counted(flat, quad, &alloc, &limits, &out, &stats);
  CHECK(r == GSHGEOM_OK || r == GSHGEOM_ERR_INVALID || r == GSHGEOM_ERR_LIMIT || r == GSHGEOM_ERR_OOM);
  CHECK((r == GSHGEOM_ERR_INVALID) == !valid);
  CHECK(gshgeom_flat_path_hash(flat) == input_hash);
  if (r != GSHGEOM_OK) {
    CHECK(out == NULL);
    CHECK(counter.live == 0);
    if (r == GSHGEOM_ERR_INVALID) {
      outcomes[OUT_INVALID_QUAD]++;
    } else if (r == GSHGEOM_ERR_LIMIT) {
      CHECK(small_limit);
      outcomes[OUT_LIMIT]++;
    } else {
      CHECK(fail_alloc && counter.failed);
      outcomes[OUT_OOM]++;
    }
    gshgeom_flat_path_destroy(flat);
    return 0;
  }
  CHECK(out != NULL);

  /* The numbers of CK5. */
  const size_t n_out = gshgeom_flat_path_point_count(out);
  CHECK(n_out <= limits.max_points);
  CHECK(n_out <= np + stats.crossings);
  CHECK(stats.crossings <= 9 * np);
  CHECK(gshgeom_flat_path_subpath_count(out) <= nc);
  if (gshgeom_flat_path_subpath_count(out) == 0) {
    outcomes[OUT_EMPTY]++;
  }
  if (stats.contours_in) {
    outcomes[OUT_CUT]++;
  }

  /* Every point lies in the box of the input's points. */
  {
    int32_t lo_x = pts[0].x, hi_x = pts[0].x, lo_y = pts[0].y, hi_y = pts[0].y;
    for (size_t i = 1; i < np; i++) {
      lo_x = pts[i].x < lo_x ? pts[i].x : lo_x;
      hi_x = pts[i].x > hi_x ? pts[i].x : hi_x;
      lo_y = pts[i].y < lo_y ? pts[i].y : lo_y;
      hi_y = pts[i].y > hi_y ? pts[i].y : hi_y;
    }
    const GSHGEOM_Device_Point * op = gshgeom_flat_path_points(out);
    for (size_t i = 0; i < n_out; i++) {
      CHECK(op[i].x >= lo_x && op[i].x <= hi_x && op[i].y >= lo_y && op[i].y <= hi_y);
    }
  }

  /* Determinism, and a clip with the default allocator and limits agrees when nothing was limited. */
  if (!small_limit) {
    GSHGEOM_Flat_Path * again = NULL;
    CHECK(gshgeom_flat_path_clip_convex(flat, quad, NULL, NULL, &again) == GSHGEOM_OK);
    CHECK(path_equal(out, again));
    gshgeom_flat_path_destroy(again);
  }

  /* A contour wholly inside the closed quadrilateral is in the result bit for bit, in order. */
  {
    const GSHGEOM_Device_Point * op = gshgeom_flat_path_points(out);
    size_t next_out = 0;
    for (size_t c = 0; c < nc; c++) {
      int inside = 1;
      for (size_t k = 0; k < subs[c].count && inside; k++) {
        for (int e = 0; e < 4 && inside; e++) {
          const GSHGEOM_Device_Point a = quad[e];
          const GSHGEOM_Device_Point b = quad[(e + 1) % 4];
          if ((a.x == b.x && a.y == b.y) || orient(a, b, pts[subs[c].first + k]) * sigma >= 0) {
            continue;
          }
          inside = 0;
        }
        for (int e = 0; e < 4; e++) {
          const GSHGEOM_Device_Point a = quad[e];
          const GSHGEOM_Device_Point b = quad[(e + 1) % 4];
          if (!(a.x == b.x && a.y == b.y) && orient(a, b, pts[subs[c].first + k]) == 0) {
            outcomes[OUT_ON_THE_LINE]++;
          }
        }
      }
      if (!inside) {
        continue;
      }
      int found = 0;
      for (; next_out < gshgeom_flat_path_subpath_count(out) && !found; next_out++) {
        GSHGEOM_Flat_Subpath so;
        CHECK(gshgeom_flat_path_subpath(out, next_out, &so) == GSHGEOM_OK);
        found = so.count == subs[c].count && so.closed == subs[c].closed &&
            memcmp(op + so.first, pts + subs[c].first, so.count * sizeof *op) == 0;
      }
      CHECK(found);
      outcomes[OUT_COPIED]++;
      if (!subs[c].closed) {
        outcomes[OUT_COPIED_OPEN]++;
      }
    }
  }

  /* The winding, in the small range. */
  if (!wide) {
    int32_t lo_x = quad[0].x, hi_x = quad[0].x, lo_y = quad[0].y, hi_y = quad[0].y;
    for (int i = 1; i < 4; i++) {
      lo_x = quad[i].x < lo_x ? quad[i].x : lo_x;
      hi_x = quad[i].x > hi_x ? quad[i].x : hi_x;
      lo_y = quad[i].y < lo_y ? quad[i].y : lo_y;
      hi_y = quad[i].y > hi_y ? quad[i].y : hi_y;
    }
    for (size_t i = 0; i < np; i++) {
      lo_x = pts[i].x < lo_x ? pts[i].x : lo_x;
      hi_x = pts[i].x > hi_x ? pts[i].x : hi_x;
      lo_y = pts[i].y < lo_y ? pts[i].y : lo_y;
      hi_y = pts[i].y > hi_y ? pts[i].y : hi_y;
    }
    GSHGEOM_Flat_Path * quad_path = NULL;
    GSHGEOM_Flat_Subpath qs = {0, 4, 1};
    CHECK(gshgeom_flat_path_from_points(quad, 4, &qs, 1, NULL, NULL, &quad_path) == GSHGEOM_OK);
    for (int i = 0; i < 9; i++) {
      for (int j = 0; j < 9; j++) {
        GSHGEOM_Device_Point p;
        p.x = (int32_t)(lo_x - 4 + ((int64_t)hi_x - lo_x + 8) * (2 * i + 1) / 18);
        p.y = (int32_t)(lo_y - 4 + ((int64_t)hi_y - lo_y + 8) * (2 * j + 1) / 18);
        if (near_path(flat, p) || near_path(quad_path, p)) {
          continue;
        }
        int in_quad = 1;
        for (int e = 0; e < 4 && in_quad; e++) {
          const GSHGEOM_Device_Point a = quad[e];
          const GSHGEOM_Device_Point b = quad[(e + 1) % 4];
          if (!(a.x == b.x && a.y == b.y) && orient(a, b, p) * sigma <= 0) {
            in_quad = 0;
          }
        }
        const int want = in_quad ? winding(flat, p) : 0;
        int got = winding(out, p);
#ifdef GSHGEOM_FUZZ_PLANT
        got += 1; /* the replay gate's plant: a result whose winding is off by one everywhere, which this must see */
#endif
        CHECK(got == want);
        outcomes[OUT_WINDING_CHECKED]++;
        if (want != 0) {
          outcomes[OUT_WINDING_NONZERO]++;
        }
      }
    }
    gshgeom_flat_path_destroy(quad_path);
  }

  gshgeom_flat_path_destroy(out);
  gshgeom_flat_path_destroy(flat);
  /* The failing allocation was never reached: a clean run holds nothing. */
  CHECK(counter.live == 0);
  return 0;
}
