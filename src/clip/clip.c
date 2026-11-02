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
 * Clipping a flattened path to a convex quadrilateral (clip.h, CK1 to CK6).
 *
 * **The cut.**  A contour is cut by one half-plane of the quadrilateral at a time.  Walking the closed polygon edge by
 * edge from `S` to `E`, with each point classified against the half-plane as strictly inside (+1), on the line (0) or
 * strictly outside (-1): a crossing is emitted when the two have opposite strict signs (`S` is +1 and `E` is -1, or the
 * other way round), and then `E` is emitted if it is not strictly outside.  A point on the line is inside, so it is kept
 * and needs no crossing.  The emitted points are joined in order, and the join from the last emitted point back to the
 * first is implicit, so wherever the contour left the half-plane and came back, the straight run along the clip line
 * stands in for the part that was outside.
 *
 * **Why the winding survives.**  Take a maximal run of the contour that is strictly outside the half-plane: it leaves
 * at a point `X` of the closed inside (a crossing, or a vertex on the line) and returns at a point `Y` of it.  The cut
 * replaces the run, `X ... Y`, by the segment `X Y`.  The run and that segment form a closed polygon, and it lies in
 * the closed outside half-plane (the segment's ends are on the clip line and the half-plane is convex), so its winding
 * number about every point strictly inside the half-plane is zero.  The cut contour is the contour plus that closed
 * polygon taken once with the sign that cancels the run, and so winds as the contour does about every point strictly
 * inside the half-plane; about a point strictly outside it winds zero, because every point emitted is in the closed
 * inside and a polygon in a closed half-plane winds zero about a point strictly outside it.  Nothing is unioned or
 * merged, so the argument holds for any closed polylines, whatever they do to one another.  Applied for each of the
 * four half-planes in turn, the result winds as the input inside the interior of the quadrilateral and zero outside it.
 *
 * **Exactness.**  The side of a point is the sign of a 2 by 2 determinant of 33-bit differences: two products of at
 * most 64 bits, compared as sign and magnitude, with no overflow.  The side value itself is twice the area of a triangle
 * inside a square of side 2^32 - 1, so it is below 2^64 in magnitude.  The crossing is `P = A + (B - A) t`,
 * `t = |fA| / (|fA| + |fB|)` for the side values `fA`, `fB` of the end points (opposite signs), so each coordinate is
 * `A + round(d |fA| / (|fA| + |fB|))` with `d` the coordinate's difference (33 bits): the numerator is below 2^97 and the
 * denominator below 2^65, and the rounding `floor(v + 1/2)` is `floor((2 d |fA| + den) / (2 den))` in integers below 2^99.  They
 * are held in a 128-bit unsigned type built from two `uint64_t` (the library avoids the compiler's 128-bit integer, which
 * a 32-bit target lacks), with a signed wrapper that keeps a sign and a magnitude, so no shift of a negative value
 * occurs anywhere.  The quotient is below 2^34.  The result is an integer function of the integers, so it is the same on
 * every target (CK6).  After the lexicographic ordering `d` is not negative for the x coordinate (the end points are
 * ordered by x first); the y coordinate's can be, and the sign is carried for both.
 */

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/shape-geometry/clip.h>
#include <ghoti.io/shape-geometry/orient.h>

#include "../core/grow.h"
#include "../flatten/flatten_internal.h"
#include "clip_internal.h"

#include <stdint.h>
#include <string.h>

typedef GSHGEOM_Device_Point Pt;

/* --- 128-bit unsigned arithmetic -------------------------------------------------------------------------------- */

typedef struct U128 {
  uint64_t hi;
  uint64_t lo;
} U128;

#define M32 UINT64_C(0xFFFFFFFF)

/** The 128-bit product of two 64-bit values. */
static U128 mul64(uint64_t a, uint64_t b) {
  const uint64_t a0 = a & M32;
  const uint64_t a1 = a >> 32;
  const uint64_t b0 = b & M32;
  const uint64_t b1 = b >> 32;
  const uint64_t p00 = a0 * b0;
  const uint64_t p01 = a0 * b1;
  const uint64_t p10 = a1 * b0;
  const uint64_t p11 = a1 * b1;
  const uint64_t mid = (p00 >> 32) + (p01 & M32) + (p10 & M32);
  U128 r;
  r.lo = (p00 & M32) | (mid << 32);
  r.hi = p11 + (p01 >> 32) + (p10 >> 32) + (mid >> 32);
  return r;
}

/** The sum; the callers' operands are below 2^127, so it does not wrap. */
static U128 add128(U128 a, U128 b) {
  U128 r;
  r.lo = a.lo + b.lo;
  r.hi = a.hi + b.hi + (r.lo < a.lo ? 1u : 0u);
  return r;
}

/** The difference `a - b` for `a >= b`. */
static U128 sub128(U128 a, U128 b) {
  U128 r;
  r.lo = a.lo - b.lo;
  r.hi = a.hi - b.hi - (a.lo < b.lo ? 1u : 0u);
  return r;
}

static int cmp128(U128 a, U128 b) {
  if (a.hi != b.hi) {
    return a.hi < b.hi ? -1 : 1;
  }
  if (a.lo != b.lo) {
    return a.lo < b.lo ? -1 : 1;
  }
  return 0;
}

static int is_zero128(U128 a) {
  return a.hi == 0 && a.lo == 0;
}

/** `a * b` for a product below 2^128 (the callers' are below 2^100). */
static U128 mul128_64(U128 a, uint64_t b) {
  U128 r = mul64(a.lo, b);
  r.hi += a.hi * b;
  return r;
}

/** The bit `i` (0 is the least significant). */
static unsigned bit128(U128 a, unsigned i) {
  return (unsigned)(i >= 64 ? (a.hi >> (i - 64)) & 1u : (a.lo >> i) & 1u);
}

/** The number of significant bits (0 for 0). */
static unsigned bits128(U128 a) {
  unsigned n = 0;
  uint64_t v = a.hi ? a.hi : a.lo;
  if (a.hi) {
    n = 64;
  }
  while (v) {
    n++;
    v >>= 1;
  }
  return n;
}

/** Quotient and remainder of `n / d` for `d != 0` and `d < 2^127`, by restoring division. */
static void divmod128(U128 n, U128 d, U128 * q, U128 * r) {
  if (n.hi == 0 && d.hi == 0) {
    q->hi = 0;
    q->lo = n.lo / d.lo;
    r->hi = 0;
    r->lo = n.lo % d.lo;
    return;
  }
  U128 quo = {0, 0};
  U128 rem = {0, 0};
  for (unsigned i = bits128(n); i-- > 0;) {
    rem.hi = (rem.hi << 1) | (rem.lo >> 63);
    rem.lo = (rem.lo << 1) | bit128(n, i);
    if (cmp128(rem, d) >= 0) {
      rem = sub128(rem, d);
      if (i >= 64) {
        quo.hi |= UINT64_C(1) << (i - 64);
      } else {
        quo.lo |= UINT64_C(1) << i;
      }
    }
  }
  *q = quo;
  *r = rem;
}

/* --- signed values as a sign and a magnitude -------------------------------------------------------------------- */

typedef struct S128 {
  int neg; /* 1 for a negative value; 0 for zero and the positive ones */
  U128 m;
} S128;

static S128 s_make(int neg, U128 m) {
  S128 r;
  r.neg = is_zero128(m) ? 0 : neg;
  r.m = m;
  return r;
}

static S128 s_add(S128 a, S128 b) {
  if (a.neg == b.neg) {
    return s_make(a.neg, add128(a.m, b.m));
  }
  if (cmp128(a.m, b.m) >= 0) {
    return s_make(a.neg, sub128(a.m, b.m));
  }
  return s_make(b.neg, sub128(b.m, a.m));
}

static S128 s_neg(S128 a) {
  return s_make(!a.neg, a.m);
}

/** The magnitude of an int64 as an unsigned value (no negation of INT64_MIN: the differences here are 33 bits). */
static uint64_t abs64(int64_t v) {
  return v < 0 ? (uint64_t)0 - (uint64_t)v : (uint64_t)v;
}

/** A product of two differences: a sign and a magnitude below 2^64. */
typedef struct Prod {
  int neg;
  uint64_t m;
} Prod;

static Prod prod(int64_t a, int64_t b) {
  Prod p;
  p.m = abs64(a) * abs64(b);
  p.neg = p.m != 0 && ((a < 0) != (b < 0));
  return p;
}

/** The side value `(l1.x - l0.x) (p.y - l0.y) - (l1.y - l0.y) (p.x - l0.x)` as the pair of products. */
static void side_products(Pt l0, Pt l1, Pt p, Prod * p1, Prod * p2) {
  *p1 = prod((int64_t)l1.x - l0.x, (int64_t)p.y - l0.y);
  *p2 = prod((int64_t)l1.y - l0.y, (int64_t)p.x - l0.x);
}

int gshgeom_clip_side(Pt l0, Pt l1, Pt p) {
  Prod p1;
  Prod p2;
  side_products(l0, l1, p, &p1, &p2);
  /* The sign of p1 - p2. */
  if (p1.neg != p2.neg) {
    return p1.neg ? -1 : 1;
  }
  if (p1.m == p2.m) {
    return 0;
  }
  const int larger = p1.m > p2.m ? 1 : -1; /* the sign of m1 - m2 */
  return p1.neg ? -larger : larger;
}

/** The side value itself (below 2^64 in magnitude). */
static S128 side_value(Pt l0, Pt l1, Pt p) {
  Prod p1;
  Prod p2;
  side_products(l0, l1, p, &p1, &p2);
  return s_add(s_make(p1.neg, (U128){0, p1.m}), s_neg(s_make(p2.neg, (U128){0, p2.m})));
}

/** Whether point (ax, ay) comes before (bx, by) in the order of the coordinates: x, then y. */
static int before(Pt a, Pt b) {
  return a.x < b.x || (a.x == b.x && a.y < b.y);
}

/**
 * `floor(v + 1/2)` for `v = (neg ? -1 : 1) num / den`, `den > 0`: the nearest integer with ties upward, exactly.  It is
 * `floor(N / D)` for `N = den + (neg ? -2 num : 2 num)` and `D = 2 den`.  The result is below 2^63 in magnitude (the
 * callers' quotients are below 2^34).
 */
static int64_t round_half_up(int neg, U128 num, U128 den) {
  const U128 twice = add128(num, num);
  const U128 d2 = add128(den, den);
  S128 n = s_add(s_make(0, den), s_make(neg, twice));
  U128 q;
  U128 r;
  divmod128(n.m, d2, &q, &r);
  if (n.neg) {
    /* floor of a negative quotient: minus the ceiling of its magnitude. */
    if (!is_zero128(r)) {
      q = add128(q, (U128){0, 1});
    }
    return -(int64_t)q.lo;
  }
  return (int64_t)q.lo;
}

int gshgeom_clip_crossing(Pt a, Pt b, Pt l0, Pt l1, Pt * out) {
  const int sa = gshgeom_clip_side(l0, l1, a);
  const int sb = gshgeom_clip_side(l0, l1, b);
  if (sa * sb >= 0) {
    return 0;
  }
  /* The end points of the edge and the points of the line in lexicographic order, so the arithmetic below is the same
   * for the same edge and line however they were named (CK3). */
  if (before(b, a)) {
    const Pt t = a;
    a = b;
    b = t;
  }
  if (before(l1, l0)) {
    const Pt t = l0;
    l0 = l1;
    l1 = t;
  }
  const S128 fa = side_value(l0, l1, a);
  const S128 fb = side_value(l0, l1, b);
  /* The side values have opposite strict signs, so |fa - fb| = |fa| + |fb| and t = |fa| / (|fa| + |fb|). */
  const U128 den = add128(fa.m, fb.m);
  const int64_t dx = (int64_t)b.x - a.x;
  const int64_t dy = (int64_t)b.y - a.y;
  const U128 nx = mul128_64(fa.m, abs64(dx));
  const U128 ny = mul128_64(fa.m, abs64(dy));
  const int64_t x = (int64_t)a.x + round_half_up(dx < 0, nx, den);
  const int64_t y = (int64_t)a.y + round_half_up(dy < 0, ny, den);
  out->x = (int32_t)x;
  out->y = (int32_t)y;
  return 1;
}

/* --- the quadrilateral ------------------------------------------------------------------------------------------- */

typedef struct Side {
  Pt a;
  Pt b;
  int sigma; /* +1 when the interior is on the positive side of a -> b, else -1 */
} Side;

static GSHGEOM_Point as_point(Pt p) {
  GSHGEOM_Point q = {(double)p.x, (double)p.y, 0.0};
  return q;
}

/**
 * The sides of the quadrilateral that cut (a side of zero length is skipped), or ::GSHGEOM_ERR_INVALID.  The turn at
 * each corner is ::gshgeom_orient2d_exact of the corner and its neighbours (CK1).
 */
static GSHGEOM_Result sides_of(const Pt quad[4], Side sides[4], int * count) {
  int turn[4];
  int sigma = 0;
  for (int i = 0; i < 4; i++) {
    const GSHGEOM_Result r = gshgeom_orient2d_exact(as_point(quad[(i + 3) % 4]), as_point(quad[i]), as_point(quad[(i + 1) % 4]), &turn[i]);
    if (r != GSHGEOM_OK) {
      return GSHGEOM_ERR_INTERNAL;
    }
    if (turn[i] != 0) {
      if (sigma != 0 && turn[i] != sigma) {
        return GSHGEOM_ERR_INVALID; /* turns of both signs: not convex */
      }
      sigma = turn[i];
    }
  }
  /* Not all collinear, and an area: the fan from corner 0 has the triangles at corners 1 and 3 (the turns there). */
  if (sigma == 0 || (turn[1] == 0 && turn[3] == 0)) {
    return GSHGEOM_ERR_INVALID;
  }
  int n = 0;
  for (int i = 0; i < 4; i++) {
    const Pt a = quad[i];
    const Pt b = quad[(i + 1) % 4];
    if (a.x == b.x && a.y == b.y) {
      continue;
    }
    sides[n].a = a;
    sides[n].b = b;
    sides[n].sigma = sigma;
    n++;
  }
  *count = n;
  return GSHGEOM_OK;
}

/** -1 strictly outside, 0 on the side's line, +1 strictly inside. */
static int classify(const Side * s, Pt p) {
  return gshgeom_clip_side(s->a, s->b, p) * s->sigma;
}

/* --- growable point buffers ------------------------------------------------------------------------------------- */

typedef struct Buf {
  Pt * p;
  size_t n;
  size_t cap;
} Buf;

static GSHGEOM_Result push(const GSHGEOM_Allocator * alloc, Buf * b, Pt v) {
  void * grown = b->p;
  const GSHGEOM_Result r = gshgeom_grow(alloc, &grown, &b->cap, b->n + 1, sizeof *b->p);
  if (r != GSHGEOM_OK) {
    return r;
  }
  b->p = grown;
  b->p[b->n++] = v;
  return GSHGEOM_OK;
}

typedef enum Cut { CUT_SAME, CUT_EMPTY, CUT_MADE } Cut;

/**
 * Cut the closed polygon (in, n) by one side.  CUT_SAME when every point is inside or on the line (the polygon is
 * unchanged and @p out is not touched), CUT_EMPTY when none is, else CUT_MADE with the cut polygon in @p out.
 */
static GSHGEOM_Result cut(const GSHGEOM_Allocator * alloc, const Side * s, const Pt * in, size_t n, Buf * out, Cut * what,
    size_t * crossings) {
  int any_in = 0;
  int any_out = 0;
  for (size_t k = 0; k < n && !(any_in && any_out); k++) {
    if (classify(s, in[k]) < 0) {
      any_out = 1;
    } else {
      any_in = 1;
    }
  }
  if (!any_out) {
    *what = CUT_SAME;
    return GSHGEOM_OK;
  }
  if (!any_in) {
    *what = CUT_EMPTY;
    return GSHGEOM_OK;
  }
  out->n = 0;
  int prev = classify(s, in[n - 1]);
  for (size_t k = 0; k < n; k++) {
    const Pt e = in[k];
    const int cur = classify(s, e);
    if (prev * cur < 0) {
      Pt x;
      const Pt from = in[k == 0 ? n - 1 : k - 1];
      if (!gshgeom_clip_crossing(from, e, s->a, s->b, &x)) {
        return GSHGEOM_ERR_INTERNAL;
      }
      (*crossings)++;
      const GSHGEOM_Result r = push(alloc, out, x);
      if (r != GSHGEOM_OK) {
        return r;
      }
    }
    if (cur >= 0) {
      const GSHGEOM_Result r = push(alloc, out, e);
      if (r != GSHGEOM_OK) {
        return r;
      }
    }
    prev = cur;
  }
  *what = CUT_MADE;
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_clip_convex_counted(const GSHGEOM_Flat_Path * path, const GSHGEOM_Device_Point quad[4],
    const GSHGEOM_Allocator * alloc, const GSHGEOM_Limits * limits, GSHGEOM_Flat_Path ** out,
    GSHGEOM_Clip_Stats * stats) {
  if (!path || !quad || !out) {
    return GSHGEOM_ERR_INVALID;
  }
  Side sides[4];
  int side_count = 0;
  GSHGEOM_Result status = sides_of(quad, sides, &side_count);
  if (status != GSHGEOM_OK) {
    return status;
  }
  if (!alloc) {
    alloc = gshgeom_allocator_default();
  }
  const size_t max_points = (limits && limits->max_points) ? limits->max_points : GSHGEOM_LIMITS_DEFAULT_MAX_POINTS;
  GSHGEOM_Clip_Stats local = {0, 0, 0};
  const Pt * points = gshgeom_flat_path_points(path);
  const size_t subpaths = gshgeom_flat_path_subpath_count(path);
  GSHGEOM_Flat_Path * result = NULL;
  Buf bufs[2] = {{NULL, 0, 0}, {NULL, 0, 0}};
  status = gshgeom_flat_path_new(alloc, limits, &result);
  for (size_t sp = 0; status == GSHGEOM_OK && sp < subpaths; sp++) {
    GSHGEOM_Flat_Subpath sub;
    if (gshgeom_flat_path_subpath(path, sp, &sub) != GSHGEOM_OK) {
      status = GSHGEOM_ERR_INTERNAL;
      break;
    }
    const Pt * cur = points + sub.first;
    size_t count = sub.count;
    int which = 0; /* the buffer the next cut writes */
    int changed = 0;
    int empty = 0;
    for (int i = 0; i < side_count && !empty && status == GSHGEOM_OK; i++) {
      Cut what = CUT_SAME;
      status = cut(alloc, &sides[i], cur, count, &bufs[which], &what, &local.crossings);
      if (status != GSHGEOM_OK) {
        break;
      }
      if (what == CUT_EMPTY) {
        empty = 1;
      } else if (what == CUT_MADE) {
        cur = bufs[which].p;
        count = bufs[which].n;
        which = 1 - which;
        changed = 1;
      }
    }
    if (status != GSHGEOM_OK || empty) {
      continue;
    }
    if (changed) {
      local.contours_in++;
      if (count < 3) {
        continue; /* a cut that leaves fewer than three points encloses nothing */
      }
    } else {
      local.contours_copied++;
    }
    if (count > max_points - gshgeom_flat_path_point_count(result)) {
      status = GSHGEOM_ERR_LIMIT;
      break;
    }
    status = gshgeom_flat_path_begin(result, cur[0]);
    for (size_t k = 1; status == GSHGEOM_OK && k < count; k++) {
      status = gshgeom_flat_path_add(result, cur[k]);
    }
    if (status == GSHGEOM_OK && (changed || sub.closed)) {
      gshgeom_flat_path_close_last(result);
    }
  }
  gcu_allocator_free(alloc, bufs[0].p);
  gcu_allocator_free(alloc, bufs[1].p);
  if (status != GSHGEOM_OK) {
    gshgeom_flat_path_destroy(result);
    return status;
  }
  if (stats) {
    *stats = local;
  }
  *out = result;
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_flat_path_clip_convex(const GSHGEOM_Flat_Path * path, const GSHGEOM_Device_Point quad[4],
    const GSHGEOM_Allocator * alloc, const GSHGEOM_Limits * limits, GSHGEOM_Flat_Path ** out) {
  return gshgeom_clip_convex_counted(path, quad, alloc, limits, out, NULL);
}
