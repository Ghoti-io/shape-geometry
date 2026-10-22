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
 * The integer flattener (flatten.h), revision 2, the proof that its polyline
 * stays within the tolerance of the curve, and the flattened path.
 *
 * This was `font`'s private flattener (`src/raster/flatten.c` at `font` commit
 * 69d1860, whose recursion it keeps exactly: the tests compare its output with
 * the vertices `font` produced for the same input, bit for bit).  `font` now
 * calls this one and keeps no copy, so there is one flattener in the family.
 * The recursion hands each vertex to an emit callback; the flattened path is
 * the callback that appends it.  The tolerance is a parameter of the curve
 * functions and of ::gshgeom_device_path_flatten_with_tolerance, not less than
 * ::GSHGEOM_FLATTEN_MIN_TOLERANCE; the path functions without one use
 * ::GSHGEOM_FLATTEN_TOLERANCE.
 *
 * A failed call frees what it made and writes nothing to its output.
 *
 * ## What is done
 *
 * A curve is halved at its midpoint until its control polygon is flat, and the
 * end point of each flat piece is a vertex.  Three things make the bound below
 * true, and each was a defect of revision 1 (story 7's port of `font`):
 *
 * 1. **The subdivision is done in widened integers**: every point is multiplied
 *    by 2^FLAT_EXTRA (256) on entry and all halving and both flatness tests
 *    run on those.  Revision 1 rounded every midpoint to 26.6, so each level
 *    moved the control points of every later level and the errors added up
 *    with the depth: measured against the exact curve it reached 13 units on a
 *    large curve, where 4 was promised.
 * 2. **A vertex is rounded to 26.6 once**, when it is emitted, to nearest with
 *    ties upward.  That is the only rounding of a coordinate that is not
 *    2^-9 of a unit.
 * 3. **The flatness threshold leaves room for that one rounding**: a piece is
 *    flat when its L1 deviation sum is at most `4 * 256 * (t - 1)` widened
 *    units, i.e. `4 * (t - 1)` units, for the tolerance `t` (at the default 4:
 *    3,072 widened, 12 units).  Revision 1 spent the whole tolerance on the flatness
 *    test (`4 * t`) and had none to spare for rounding.
 *
 * No float, no libm, no shift of a negative value (a widening is a
 * multiplication, a narrowing a floor division), no overflow: the proof below
 * gives the ranges.  The result is the same on every target.
 *
 * ## The theorem
 *
 * Let the curve have control points with int32 coordinates (any value of the
 * type) and let `t >= 2` be the tolerance in 26.6 units (the default is
 * ::GSHGEOM_FLATTEN_TOLERANCE, 4).  Then every point of the curve is within
 * `t - 0.12` units (so, at the default, 3.88, and within `t`) of the emitted
 * polyline, and the depth cap never forces a piece.  The proof is for every
 * `t >= 2` up to INT32_MAX because the threshold is derived from `t`; a smaller
 * tolerance is refused with ::GSHGEOM_ERR_INVALID, before anything is emitted.
 *
 * ## The proof
 *
 * Units are 26.6 units unless called "widened" (1 widened = 1/256 unit).  Write
 * W = 256, T = 4W(t-1), E = 30.
 *
 * **(a) Ranges.** Every coordinate of every point the recursion makes lies in
 * the interval spanned by the corresponding input coordinates, widened, so
 * within +-2^39 widened (and the threshold `4W(t-1)`, with `t - 1 < 2^31`, is
 * below 2^41): the inputs are there, and the half-up midpoint
 * `floor((u+v+1)/2)` of two integers in an interval `[lo, hi]` is in it (the
 * sum is at most `2 hi + 1`, whose half floors to `hi`).  Hence 3p, 2p and sums
 * of up to six such terms are below 2^42, far inside int64.  Narrowing
 * `floor((v + 128) / 256)` of a widened value in that interval returns a value
 * in [-2^31, 2^31-1] (monotonic, and exact at the end points), so the cast to
 * int32 does not truncate.  No statement depends on the input being a glyph.
 *
 * **(b) The midpoints are close to the true ones.** Take the control points Q
 * of a piece reached by `k` halvings and the control points P of the same
 * parameter interval of the exact curve (the exact de Casteljau points, which
 * are rational).  If every point of the parent is within `e` of its exact
 * counterpart on each axis, each new point is: it is the average of two parent
 * points (error at most `e`) plus the half-up rounding, at most 1/2 widened.  A
 * quadratic produces `a, b` (error +1/2) and `m = mid(a, b)` (+1); a cubic
 * `a, b, c` (+1/2), `d, e` (+1) and `m` (+3/2).  The root's points are exact.
 * So after `k <= GSHGEOM_FLATTEN_MAX_DEPTH = 20` halvings every point is within
 * `1.5 * 20 = 30 = E` widened per axis, for either degree.
 *
 * **(c) The flatness test measures the piece it has.** The test is made on Q.
 * For a quadratic piece, `B(s) - L(s) = s(1-s)(2Q1 - Q0 - Q2)`, where L is the
 * chord, because the Bernstein form minus linear interpolation leaves only the
 * middle coefficient.  So the piece is within `|2Q1-Q0-Q2|_2 / 4 <= (dx+dy) / 4`
 * of the chord, and the test accepts at `dx+dy <= T`: within `T/4 = W(t-1)`
 * widened, which is `t - 1` units.  For a cubic, with `e1 = Q1 - (2Q0+Q3)/3`
 * and `e2 = Q2 - (Q0+2Q3)/3`, the same identity gives
 * `B(s) - L(s) = 3s(1-s)((1-s) e1 + s e2)`, at most `(3/4) max(|e1|, |e2|)`;
 * and `|e1|_2 <= |e1|_1 = d1/3`, `|e2|_2 <= d2/3` for the sums d1, d2 the code
 * computes.  So the piece is within `max(d1, d2) / 4 <= T/4`: `t - 1` units.
 *
 * **(d) The piece is near the exact curve.** The piece's Bezier curve, built on
 * Q, and the exact curve over the same interval, built on P, are the same
 * convex combination of Q and of P at every parameter, so they are within
 * `max |Q_i - P_i|_2 <= sqrt(2) E / 256 <= 0.166` units of each other at the
 * same parameter, by (b).
 *
 * **(e) The vertex rounding.** The segment emitted joins `narrow(Q0)` and
 * `narrow(Qk)`: the previous vertex is the narrowing of the same point (the
 * shared `m`; for the first piece the input point, which narrows to itself).
 * Each narrowing moves a point by at most 1/2 unit per axis, so at most
 * `sqrt(2)/2 = 0.7072` units, and a point of the chord is a convex combination
 * of the two ends, so it moves by at most that.
 *
 * **(f) Together**, for any point B(s) of the curve and the piece containing s:
 * `0.166 (d) + (t - 1) (c) + 0.7072 (e) = t - 0.1268 < t - 0.12`.  QED, given
 * that the cap does not bind:
 *
 * **(g) The depth cap does not bind** for `t >= 2`, over the whole int32 range.
 * A piece is forced only after 20 halvings; its control polygon then has the
 * exact deviation of the sub-curve plus the error of (b).  The exact part: a
 * quadratic's `2Q1-Q0-Q2` for a piece of parameter length `h` is
 * `h^2 (P0 - 2P1 + P2)` (per axis at most `2^33 h^2` units), so with
 * `h = 2^-20` it is at most `2^34 / 2^40` units for both axes together: 4
 * widened.  A cubic's `3Q1-2Q0-Q3` is the negative of a Taylor remainder of the
 * piece, `<= (h^2/2) max |B''|` per axis, and `|B''| <= 6 * 2^33` per axis for
 * int32 control points, so `d1` and `d2` are each at most
 * `2 * 3 * 2^33 * 2^-40` units = 12 widened.  The error part: coefficients
 * (2,1,1) and (3,2,1) times `E` per axis, both axes: at most `8 * 20 = 160`
 * widened for a quadratic and `12 * 30 = 360` for a cubic.  The most is 372,
 * and the smallest threshold, at `t = 2`, is `T = 1024` (at the default `t = 4`
 * it is 3072).  So every piece is accepted at or before depth 20.
 * (This is why the depth is 20 and was 16: at 16 a curve spanning the int32
 * range is still 4 (quadratic) or 12 (cubic) units of deviation sum from flat.)
 *
 * ## What checks the proof
 *
 * tests/unit/test_tolerance.cpp evaluates the curve exactly and measures the
 * distance to the polyline, on every size class up to the whole int32 range and
 * at tolerances from 2 to 1000, and asserts the bound and a zero cap count; it
 * pins vertex lists from an independent big-integer model.
 * tests/unit/test_flatten.cpp does the same at the default tolerance through
 * the path functions and pins the vectors `font` recorded.  That they find a
 * violation was shown by planting (tools/check-planted.py).
 */

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/shape-geometry/flatten.h>
#include <ghoti.io/shape-geometry/hash.h>

#include "../core/grow.h"
#include "flatten_internal.h"

#include <string.h>

typedef struct Sub {
  size_t first;
  size_t count;
  uint8_t closed;
} Sub;

struct GSHGEOM_Flat_Path {
  GSHGEOM_Allocator alloc;
  size_t max_points;
  GSHGEOM_Device_Point * pt;
  size_t pt_n, pt_cap;
  Sub * sub;
  size_t sub_n, sub_cap;
};

typedef GSHGEOM_Device_Point Pt;

static GSHGEOM_Result add_point(GSHGEOM_Flat_Path * fp, Pt p) {
  if (fp->pt_n >= fp->max_points) {
    return GSHGEOM_ERR_LIMIT;
  }
  void * grown = fp->pt;
  GSHGEOM_Result r = gshgeom_grow(
      &fp->alloc, &grown, &fp->pt_cap, fp->pt_n + 1, sizeof *fp->pt);
  if (r != GSHGEOM_OK) {
    return r;
  }
  fp->pt = grown;
  fp->pt[fp->pt_n++] = p;
  fp->sub[fp->sub_n - 1].count++;
  return GSHGEOM_OK;
}

static GSHGEOM_Result begin_subpath(GSHGEOM_Flat_Path * fp, Pt p) {
  void * grown = fp->sub;
  GSHGEOM_Result r = gshgeom_grow(
      &fp->alloc, &grown, &fp->sub_cap, fp->sub_n + 1, sizeof *fp->sub);
  if (r != GSHGEOM_OK) {
    return r;
  }
  fp->sub = grown;
  fp->sub[fp->sub_n].first = fp->pt_n;
  fp->sub[fp->sub_n].count = 0;
  fp->sub[fp->sub_n].closed = 0;
  fp->sub_n++;
  r = add_point(fp, p);
  if (r != GSHGEOM_OK) {
    fp->sub_n--;
  }
  return r;
}

/** Extra fractional bits the subdivision carries beyond 26.6. */
#define FLAT_EXTRA 8

/** 2^FLAT_EXTRA, the widening factor. */
#define FLAT_SCALE ((int64_t)1 << FLAT_EXTRA)

/* The proof (a theorem for t >= 2) needs the tolerance to be at least 2, and
 * the default is one the proof covers. */
_Static_assert(GSHGEOM_FLATTEN_MIN_TOLERANCE >= 2,
    "the flattener's guarantee is proved for tolerances of 2 units and up");
_Static_assert(GSHGEOM_FLATTEN_TOLERANCE >= GSHGEOM_FLATTEN_MIN_TOLERANCE,
    "the default tolerance is below the smallest one accepted");

/** A point in widened units: 26.6 times 2^FLAT_EXTRA. */
typedef struct WPt {
  int64_t x;
  int64_t y;
} WPt;

/** One curve's recursion: the threshold derived from the tolerance, and where
 * the vertices go. */
typedef struct Flattener {
  int64_t threshold;          /* The L1 flatness limit, widened: 4 * 256 * (t - 1). */
  GSHGEOM_Flatten_Emit emit;
  void * user;
  int levels;                 /* The depth cap: 1 to GSHGEOM_FLATTEN_MAX_DEPTH. */
  int forced;                 /* Pieces the depth cap ended although not flat. */
} Flattener;

/* floor(a / b) for b > 0, without shifting a negative value. */
static int64_t floor_div(int64_t a, int64_t b) {
  int64_t q = a / b;
  if (a % b != 0 && a < 0) {
    q -= 1;
  }
  return q;
}

/* A multiplication, not a shift: the operand may be negative. */
static WPt widen(Pt p) {
  WPt w;
  w.x = (int64_t)p.x * FLAT_SCALE;
  w.y = (int64_t)p.y * FLAT_SCALE;
  return w;
}

/* A widened point back to 26.6: to nearest, ties upward. */
static Pt narrow(WPt w) {
  Pt p;
  p.x = (int32_t)floor_div(w.x + FLAT_SCALE / 2, FLAT_SCALE);
  p.y = (int32_t)floor_div(w.y + FLAT_SCALE / 2, FLAT_SCALE);
  return p;
}

/* The midpoint of two widened points: to nearest, ties upward. */
static WPt mid_pt(WPt a, WPt b) {
  WPt m;
  m.x = floor_div(a.x + b.x + 1, 2);
  m.y = floor_div(a.y + b.y + 1, 2);
  return m;
}

static int64_t abs64(int64_t v) {
  return v < 0 ? -v : v;
}

/* The vertex a piece ends at, rounded to 26.6 once, handed to the sink. */
static GSHGEOM_Result emit_vertex(
    Flattener * fl, WPt to, uint32_t index, int level) {
  return fl->emit(fl->user, narrow(to), index, level);
}

/* A quadratic piece is flat when |2 p1 - p0 - p2| summed over both axes is at
 * most the threshold (its deviation from its chord is a quarter of that).  The
 * flatness test comes before the cap, as in the proof. */
static GSHGEOM_Result quad(
    Flattener * fl, WPt p0, WPt p1, WPt p2, uint32_t index, int level) {
  int64_t sum = abs64(2 * p1.x - p0.x - p2.x) + abs64(2 * p1.y - p0.y - p2.y);
  if (sum <= fl->threshold) {
    return emit_vertex(fl, p2, index, level);
  }
  if (level >= fl->levels) {
    fl->forced += 1;
    return emit_vertex(fl, p2, index, level);
  }
  WPt a = mid_pt(p0, p1);
  WPt b = mid_pt(p1, p2);
  WPt m = mid_pt(a, b);
  GSHGEOM_Result r = quad(fl, p0, a, m, index * 2u, level + 1);
  return r != GSHGEOM_OK ? r : quad(fl, m, b, p2, index * 2u + 1u, level + 1);
}

/* A cubic piece is flat when both control-point deviations d1 and d2, each
 * summed over both axes, are at most the threshold. */
static GSHGEOM_Result cubic(Flattener * fl, WPt p0, WPt p1, WPt p2, WPt p3,
    uint32_t index, int level) {
  int64_t d1 = abs64(3 * p1.x - 2 * p0.x - p3.x) +
      abs64(3 * p1.y - 2 * p0.y - p3.y);
  int64_t d2 = abs64(3 * p2.x - p0.x - 2 * p3.x) +
      abs64(3 * p2.y - p0.y - 2 * p3.y);
  if (d1 <= fl->threshold && d2 <= fl->threshold) {
    return emit_vertex(fl, p3, index, level);
  }
  if (level >= fl->levels) {
    fl->forced += 1;
    return emit_vertex(fl, p3, index, level);
  }
  WPt a = mid_pt(p0, p1);
  WPt b = mid_pt(p1, p2);
  WPt c = mid_pt(p2, p3);
  WPt d = mid_pt(a, b);
  WPt e = mid_pt(b, c);
  WPt m = mid_pt(d, e);
  GSHGEOM_Result r = cubic(fl, p0, a, d, m, index * 2u, level + 1);
  return r != GSHGEOM_OK ? r
                         : cubic(fl, m, e, c, p3, index * 2u + 1u, level + 1);
}

/* The state for one curve, or ::GSHGEOM_ERR_INVALID for a tolerance the proof
 * does not cover (or no sink).  The threshold is derived from the tolerance:
 * 4 * 256 * (t - 1) is below 2^41 for any int32 t. */
static GSHGEOM_Result flattener_begin(Flattener * fl, int32_t tolerance,
    int levels, GSHGEOM_Flatten_Emit emit, void * user) {
  if (tolerance < GSHGEOM_FLATTEN_MIN_TOLERANCE || !emit) {
    return GSHGEOM_ERR_INVALID;
  }
  fl->threshold = (int64_t)4 * FLAT_SCALE * ((int64_t)tolerance - 1);
  fl->emit = emit;
  fl->user = user;
  fl->levels = levels < 1 ? 1
      : (levels > GSHGEOM_FLATTEN_MAX_DEPTH ? GSHGEOM_FLATTEN_MAX_DEPTH : levels);
  fl->forced = 0;
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_flatten_quadratic_capped(int32_t tolerance, int levels,
    Pt p0, Pt p1, Pt p2, GSHGEOM_Flatten_Emit emit, void * user, int * forced) {
  Flattener fl;
  if (forced) {
    *forced = 0;
  }
  GSHGEOM_Result r = flattener_begin(&fl, tolerance, levels, emit, user);
  if (r != GSHGEOM_OK) {
    return r;
  }
  r = quad(&fl, widen(p0), widen(p1), widen(p2), 0u, 0);
  if (forced) {
    *forced = fl.forced;
  }
  return r;
}

GSHGEOM_Result gshgeom_flatten_cubic_capped(int32_t tolerance, int levels,
    Pt p0, Pt p1, Pt p2, Pt p3, GSHGEOM_Flatten_Emit emit, void * user,
    int * forced) {
  Flattener fl;
  if (forced) {
    *forced = 0;
  }
  GSHGEOM_Result r = flattener_begin(&fl, tolerance, levels, emit, user);
  if (r != GSHGEOM_OK) {
    return r;
  }
  r = cubic(&fl, widen(p0), widen(p1), widen(p2), widen(p3), 0u, 0);
  if (forced) {
    *forced = fl.forced;
  }
  return r;
}

GSHGEOM_Result gshgeom_flatten_quadratic(int32_t tolerance, Pt p0, Pt p1,
    Pt p2, GSHGEOM_Flatten_Emit emit, void * user) {
  return gshgeom_flatten_quadratic_capped(
      tolerance, GSHGEOM_FLATTEN_MAX_DEPTH, p0, p1, p2, emit, user, NULL);
}

GSHGEOM_Result gshgeom_flatten_cubic(int32_t tolerance, Pt p0, Pt p1, Pt p2,
    Pt p3, GSHGEOM_Flatten_Emit emit, void * user) {
  return gshgeom_flatten_cubic_capped(tolerance, GSHGEOM_FLATTEN_MAX_DEPTH, p0,
      p1, p2, p3, emit, user, NULL);
}

int gshgeom_flatten_revision(void) {
  return GSHGEOM_FLATTEN_REVISION;
}

static void flat_destroy_parts(GSHGEOM_Flat_Path * fp) {
  GSHGEOM_Allocator alloc = fp->alloc;
  gcu_allocator_free(&alloc, fp->pt);
  gcu_allocator_free(&alloc, fp->sub);
  gcu_allocator_free(&alloc, fp);
}

void gshgeom_flat_path_destroy(GSHGEOM_Flat_Path * path) {
  if (path) {
    flat_destroy_parts(path);
  }
}

/* The sink the flattened path gives the curve functions: the vertex is a point. */
static GSHGEOM_Result append_vertex(
    void * user, GSHGEOM_Device_Point to, uint32_t index, int level) {
  (void)index;
  (void)level;
  return add_point(user, to);
}

GSHGEOM_Result gshgeom_device_path_flatten(const GSHGEOM_Device_Path * path,
    const GSHGEOM_Allocator * alloc, const GSHGEOM_Limits * limits,
    GSHGEOM_Flat_Path ** out) {
  return gshgeom_device_path_flatten_with_tolerance(
      path, GSHGEOM_FLATTEN_TOLERANCE, alloc, limits, out);
}

GSHGEOM_Result gshgeom_device_path_flatten_with_tolerance(
    const GSHGEOM_Device_Path * path, int32_t tolerance,
    const GSHGEOM_Allocator * alloc, const GSHGEOM_Limits * limits,
    GSHGEOM_Flat_Path ** out) {
  if (!path || !out || tolerance < GSHGEOM_FLATTEN_MIN_TOLERANCE) {
    return GSHGEOM_ERR_INVALID;
  }
  if (!alloc) {
    alloc = gshgeom_allocator_default();
  }
  GSHGEOM_Flat_Path * fp = gcu_allocator_calloc(alloc, 1, sizeof *fp);
  if (!fp) {
    return GSHGEOM_ERR_OOM;
  }
  fp->alloc = *alloc;
  fp->max_points = (limits && limits->max_points)
      ? limits->max_points
      : GSHGEOM_LIMITS_DEFAULT_MAX_POINTS;
  GSHGEOM_Result r = GSHGEOM_OK;
  Pt current = {0, 0};
  size_t n = gshgeom_device_path_element_count(path);
  for (size_t i = 0; i < n && r == GSHGEOM_OK; i++) {
    GSHGEOM_Device_Element el;
    if (gshgeom_device_path_element(path, i, &el) != GSHGEOM_OK) {
      r = GSHGEOM_ERR_INTERNAL;
      break;
    }
    /* The builder guarantees a move first and after every close. */
    if (el.verb != GSHGEOM_DEVICE_MOVE && fp->sub_n == 0) {
      r = GSHGEOM_ERR_INTERNAL;
      break;
    }
    switch (el.verb) {
    case GSHGEOM_DEVICE_MOVE:
      r = begin_subpath(fp, el.p[0]);
      current = el.p[0];
      break;
    case GSHGEOM_DEVICE_LINE:
      r = add_point(fp, el.p[0]);
      current = el.p[0];
      break;
    case GSHGEOM_DEVICE_QUADRATIC:
      r = gshgeom_flatten_quadratic(
          tolerance, current, el.p[0], el.p[1], append_vertex, fp);
      current = el.p[1];
      break;
    case GSHGEOM_DEVICE_CUBIC:
      r = gshgeom_flatten_cubic(tolerance, current, el.p[0], el.p[1],
          el.p[2], append_vertex, fp);
      current = el.p[2];
      break;
    case GSHGEOM_DEVICE_CLOSE:
      fp->sub[fp->sub_n - 1].closed = 1;
      break;
    case GSHGEOM_DEVICE_VERB_COUNT:
      r = GSHGEOM_ERR_INTERNAL;
      break;
    }
  }
  if (r != GSHGEOM_OK) {
    flat_destroy_parts(fp);
    return r;
  }
  *out = fp;
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_path_flatten(const GSHGEOM_Path * path,
    const GSHGEOM_Map * map, const GSHGEOM_Allocator * alloc,
    const GSHGEOM_Limits * limits, GSHGEOM_Flat_Path ** out) {
  if (!out) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Device_Path * dp;
  GSHGEOM_Result r = gshgeom_device_path_from_model(path, map, alloc, limits, &dp);
  if (r != GSHGEOM_OK) {
    return r;
  }
  r = gshgeom_device_path_flatten(dp, alloc, limits, out);
  gshgeom_device_path_destroy(dp);
  return r;
}

/* ---------------------------------------------------------------------- */
/* The builder the stroke expansion uses (flatten_internal.h)              */
/* ---------------------------------------------------------------------- */

GSHGEOM_Result gshgeom_flat_path_new(const GSHGEOM_Allocator * alloc,
    const GSHGEOM_Limits * limits, GSHGEOM_Flat_Path ** out) {
  if (!alloc) {
    alloc = gshgeom_allocator_default();
  }
  GSHGEOM_Flat_Path * fp = gcu_allocator_calloc(alloc, 1, sizeof *fp);
  if (!fp) {
    return GSHGEOM_ERR_OOM;
  }
  fp->alloc = *alloc;
  fp->max_points = (limits && limits->max_points) ? limits->max_points : GSHGEOM_LIMITS_DEFAULT_MAX_POINTS;
  *out = fp;
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_flat_path_begin(GSHGEOM_Flat_Path * path, GSHGEOM_Device_Point p) {
  return begin_subpath(path, p);
}

GSHGEOM_Result gshgeom_flat_path_add(GSHGEOM_Flat_Path * path, GSHGEOM_Device_Point p) {
  return add_point(path, p);
}

void gshgeom_flat_path_close_last(GSHGEOM_Flat_Path * path) {
  if (path->sub_n) {
    path->sub[path->sub_n - 1].closed = 1;
  }
}

/* ---------------------------------------------------------------------- */
/* Readback                                                                */
/* ---------------------------------------------------------------------- */

size_t gshgeom_flat_path_point_count(const GSHGEOM_Flat_Path * path) {
  return path ? path->pt_n : 0;
}

size_t gshgeom_flat_path_subpath_count(const GSHGEOM_Flat_Path * path) {
  return path ? path->sub_n : 0;
}

const GSHGEOM_Device_Point * gshgeom_flat_path_points(
    const GSHGEOM_Flat_Path * path) {
  return (path && path->pt_n) ? path->pt : NULL;
}

GSHGEOM_Result gshgeom_flat_path_subpath(
    const GSHGEOM_Flat_Path * path, size_t index, GSHGEOM_Flat_Subpath * out) {
  if (!path || !out || index >= path->sub_n) {
    return GSHGEOM_ERR_INVALID;
  }
  out->first = path->sub[index].first;
  out->count = path->sub[index].count;
  out->closed = path->sub[index].closed;
  return GSHGEOM_OK;
}

/* ---------------------------------------------------------------------- */
/* Hash and dump                                                           */
/* ---------------------------------------------------------------------- */

/** The format version of the hash and the dump. */
#define FLAT_FORMAT_VERSION 1u
/** The kind field: "FLAT" as little-endian bytes. */
#define FLAT_KIND 0x54414c46u

uint64_t gshgeom_flat_path_hash(const GSHGEOM_Flat_Path * path) {
  GSHGEOM_Hash h;
  gshgeom_hash_init(&h);
  if (!path) {
    return gshgeom_hash_value(&h);
  }
  gshgeom_hash_u32(&h, FLAT_FORMAT_VERSION);
  gshgeom_hash_u32(&h, FLAT_KIND);
  gshgeom_hash_u64(&h, (uint64_t)path->sub_n);
  for (size_t i = 0; i < path->sub_n; i++) {
    const Sub * s = &path->sub[i];
    gshgeom_hash_u8(&h, s->closed);
    gshgeom_hash_u64(&h, (uint64_t)s->count);
    for (size_t k = 0; k < s->count; k++) {
      gshgeom_hash_i32(&h, path->pt[s->first + k].x);
      gshgeom_hash_i32(&h, path->pt[s->first + k].y);
    }
  }
  return gshgeom_hash_value(&h);
}

GSHGEOM_Result gshgeom_flat_path_dump(
    const GSHGEOM_Flat_Path * path, FILE * out) {
  if (!path || !out) {
    return GSHGEOM_ERR_INVALID;
  }
  if (fprintf(out, "gshgeom.flat-path %u subpaths=%llu points=%llu hash=%016llx\n",
          FLAT_FORMAT_VERSION, (unsigned long long)path->sub_n,
          (unsigned long long)path->pt_n,
          (unsigned long long)gshgeom_flat_path_hash(path)) < 0) {
    return GSHGEOM_ERR_IO;
  }
  for (size_t i = 0; i < path->sub_n; i++) {
    const Sub * s = &path->sub[i];
    if (fprintf(out, "S closed=%d points=%llu\n", s->closed ? 1 : 0,
            (unsigned long long)s->count) < 0) {
      return GSHGEOM_ERR_IO;
    }
    for (size_t k = 0; k < s->count; k++) {
      if (fprintf(out, "P %08x %08x\n",
              (unsigned)(uint32_t)path->pt[s->first + k].x,
              (unsigned)(uint32_t)path->pt[s->first + k].y) < 0) {
        return GSHGEOM_ERR_IO;
      }
    }
  }
  /* A buffered stream reports a write error at the flush, and a stream's error flag is sticky:
   * a failed write whose return value was lost (an unbuffered stream's fprintf) still sets it. */
  return (fflush(out) == 0 && !ferror(out)) ? GSHGEOM_OK : GSHGEOM_ERR_IO;
}
