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
 * The exact orientation predicates, the normals and the signed area (see
 * orient.h for the rules).
 */

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/shape-geometry/orient.h>

#include "../exact/exact.h"

/** A finite double: v - v is zero for every finite v. */
static int is_finite(double v) {
  return (v - v) == 0.0;
}

static double absolute(double v) {
  return v < 0.0 ? -v : v;
}

/* The filter's static error bounds are Shewchuk's published ones, for round to
 * nearest with eps = 2^-53 (half an ulp of 1):
 *   orient2d  (3 + 16 eps) eps  times |detleft| + |detright|
 *   orient3d  (7 + 56 eps) eps  times the permanent of the 3x3 of differences
 * They hold while nothing overflows or underflows.  The filter is therefore
 * used only for coordinates up to 2^100 (differences and products are far
 * below the overflow threshold) and when the quantity the bound scales is at
 * least 2^-900 (an underflowed product then adds an absolute error of 2^-1074,
 * which is 2^-174 times smaller than the bound), and everything else goes to
 * the exact path. */
#define EPS 0x1p-53
#define CCW_ERR_BOUND ((3.0 + 16.0 * EPS) * EPS)
#define O3D_ERR_BOUND ((7.0 + 56.0 * EPS) * EPS)
#define FILTER_COORD_LIMIT 0x1p100
#define FILTER_MIN 0x1p-900

static int small_enough(double v) {
  return v <= FILTER_COORD_LIMIT && v >= -FILTER_COORD_LIMIT;
}

/** Whether every coordinate is within the filter's limit (which no NaN and no infinity is). */
static int all_small2(GSHGEOM_Point a, GSHGEOM_Point b, GSHGEOM_Point c) {
  return small_enough(a.x) & small_enough(a.y) & small_enough(b.x) & small_enough(b.y) &
      small_enough(c.x) & small_enough(c.y);
}

/** 1 and *sign set if the filter decides, for coordinates that are all within the limit. */
static int filter2_small(GSHGEOM_Point a, GSHGEOM_Point b, GSHGEOM_Point c, int * sign) {
  double detleft = (a.x - c.x) * (b.y - c.y);
  double detright = (a.y - c.y) * (b.x - c.x);
  double det = detleft - detright;
  double detsum = absolute(detleft) + absolute(detright);
  double errbound = CCW_ERR_BOUND * detsum;
  if (detsum >= FILTER_MIN && (det >= errbound || -det >= errbound)) {
    *sign = det > 0.0 ? 1 : -1;
    return 1;
  }
  return 0;
}

static int filter3(GSHGEOM_Point a, GSHGEOM_Point b, GSHGEOM_Point c,
    GSHGEOM_Point d, int * sign) {
  if (!(small_enough(a.x) && small_enough(a.y) && small_enough(a.z) &&
          small_enough(b.x) && small_enough(b.y) && small_enough(b.z) &&
          small_enough(c.x) && small_enough(c.y) && small_enough(c.z) &&
          small_enough(d.x) && small_enough(d.y) && small_enough(d.z))) {
    return 0;
  }
  double adx = a.x - d.x;
  double bdx = b.x - d.x;
  double cdx = c.x - d.x;
  double ady = a.y - d.y;
  double bdy = b.y - d.y;
  double cdy = c.y - d.y;
  double adz = a.z - d.z;
  double bdz = b.z - d.z;
  double cdz = c.z - d.z;
  double bdxcdy = bdx * cdy;
  double cdxbdy = cdx * bdy;
  double cdxady = cdx * ady;
  double adxcdy = adx * cdy;
  double adxbdy = adx * bdy;
  double bdxady = bdx * ady;
  double det = adz * (bdxcdy - cdxbdy) + bdz * (cdxady - adxcdy) +
      cdz * (adxbdy - bdxady);
  double permanent = (absolute(bdxcdy) + absolute(cdxbdy)) * absolute(adz) +
      (absolute(cdxady) + absolute(adxcdy)) * absolute(bdz) +
      (absolute(adxbdy) + absolute(bdxady)) * absolute(cdz);
  double errbound = O3D_ERR_BOUND * permanent;
  if (permanent >= FILTER_MIN && (det > errbound || -det > errbound)) {
    /* This det is the determinant of a-d, b-d, c-d, the negative of ours. */
    *sign = det > 0.0 ? -1 : 1;
    return 1;
  }
  return 0;
}

/* The exact forms, in the raw coordinates (so that no difference is rounded).
 *
 * orient2d is the determinant of the rows (a, 1), (b, 1), (c, 1):
 *   ax by - ax cy - ay bx + ay cx + bx cy - by cx.
 *
 * orient3d is the determinant of b-a, c-a, d-a.  The determinant of the four
 * rows (a, 1) .. (d, 1) expanded along its last column is
 *   -D(b,c,d) + D(a,c,d) - D(a,b,d) + D(a,b,c)
 * with D the 3x3 determinant of the coordinates, and it is the determinant of
 * a-d, b-d, c-d; subtracting the first row instead gives it as minus the
 * determinant of b-a, c-a, d-a.  So ours is
 *   D(b,c,d) - D(a,c,d) + D(a,b,d) - D(a,b,c). */
static void put(GSHGEOM_ExactTerm * t, size_t * n, int neg, double f0, double f1) {
  gshgeom_exact_put(t, n, neg, 2, f0, f1, 0.0);
}

/** Add `(-1)^neg D(r1, r2, r3)`, the 3x3 determinant of the coordinates. */
static void put_det3(GSHGEOM_ExactTerm * t, size_t * n, int neg,
    GSHGEOM_Point r1, GSHGEOM_Point r2, GSHGEOM_Point r3) {
  const double a[3] = {r1.x, r1.y, r1.z};
  const double b[3] = {r2.x, r2.y, r2.z};
  const double c[3] = {r3.x, r3.y, r3.z};
  gshgeom_exact_put_det3(t, n, neg, a, b, c);
}

static int exact2(GSHGEOM_Point a, GSHGEOM_Point b, GSHGEOM_Point c) {
  GSHGEOM_ExactTerm t[6];
  size_t n = 0;
  put(t, &n, 0, a.x, b.y);
  put(t, &n, 1, a.x, c.y);
  put(t, &n, 1, a.y, b.x);
  put(t, &n, 0, a.y, c.x);
  put(t, &n, 0, b.x, c.y);
  put(t, &n, 1, b.y, c.x);
  return gshgeom_exact_sign(t, n);
}

/* The exact sign of orient2d by floating-point expansions (Shewchuk's two-product and grow-expansion), for the cases the filter cannot
 * certify, which are mostly points that are exactly collinear (the filter can never certify a zero).
 *
 * The sum is the one exact2 forms, ax by - ax cy - ay bx + ay cx + bx cy - by cx, each product split into a rounded product and the
 * exact rounding error (Dekker's split, which needs no fused multiply-add and so no particular target), and the twelve doubles added
 * into a nonoverlapping expansion whose largest component carries the sign.  Every step is exact while no product or error term
 * overflows or underflows, which holds when each coordinate is zero or has a magnitude in [2^-400, 2^400]: a product is then within
 * [2^-800, 2^800], its error term within 2^-53 of that and above 2^-1000 (a split half is at least 2^-453 and the product of the two
 * low halves is the smallest piece, above 2^-906 and so above the smallest normal double).  A coordinate outside that range, which
 * no drawing makes, goes to the multiword sum, so that the result is the exact sign everywhere (the tests compare the two). */
#define EXPANSION_LO 0x1p-400
#define EXPANSION_HI 0x1p400

static int expansion_ok(double v) {
  const double m = absolute(v);
  return m == 0.0 || (m >= EXPANSION_LO && m <= EXPANSION_HI);
}

static void split_double(double a, double * hi, double * lo) {
  const double c = 134217729.0 * a; /* 2^27 + 1 */
  const double big = c - a;
  *hi = c - big;
  *lo = a - *hi;
}

static void two_product(double a, double b, double * x, double * y) {
  double ahi, alo, bhi, blo;
  split_double(a, &ahi, &alo);
  split_double(b, &bhi, &blo);
  const double p = a * b;
  const double e1 = p - ahi * bhi;
  const double e2 = e1 - alo * bhi;
  const double e3 = e2 - ahi * blo;
  *x = p;
  *y = alo * blo - e3;
}

/** Add @p v to the expansion e[0..*n) (increasing magnitude, nonoverlapping, no zeros), dropping zeros: Shewchuk's grow-expansion. */
static void grow(double * e, size_t * n, double v) {
  double q = v;
  size_t k = 0;
  for (size_t i = 0; i < *n; i++) {
    const double sum = q + e[i];
    const double bv = sum - q;
    const double av = sum - bv;
    const double h = (q - av) + (e[i] - bv);
    q = sum;
    if (h != 0.0) {
      e[k++] = h;
    }
  }
  if (q != 0.0 || k == 0) {
    e[k++] = q;
  }
  *n = k;
}

/** 1 and *sign set when every coordinate is in the range the expansions are exact for, else 0. */
static int exact2_expansion(GSHGEOM_Point a, GSHGEOM_Point b, GSHGEOM_Point c, int * sign) {
  if (!(expansion_ok(a.x) && expansion_ok(a.y) && expansion_ok(b.x) && expansion_ok(b.y) && expansion_ok(c.x) && expansion_ok(c.y))) {
    return 0;
  }
  const double fa[6] = {a.x, a.x, a.y, a.y, b.x, b.y};
  const double fb[6] = {b.y, c.y, b.x, c.x, c.y, c.x};
  const int neg[6] = {0, 1, 1, 0, 0, 1};
  double e[12];
  size_t n = 0;
  for (int i = 0; i < 6; i++) {
    double hi, lo;
    two_product(fa[i], fb[i], &hi, &lo);
    if (lo != 0.0) {
      grow(e, &n, neg[i] ? 0.0 - lo : lo);
    }
    if (hi != 0.0) {
      grow(e, &n, neg[i] ? 0.0 - hi : hi);
    }
  }
  *sign = n == 0 || e[n - 1] == 0.0 ? 0 : e[n - 1] > 0.0 ? 1 : -1;
  return 1;
}

/** The exact sign of orient2d: expansions where they are exact, else the multiword sum. */
static int exact2_fast(GSHGEOM_Point a, GSHGEOM_Point b, GSHGEOM_Point c) {
  int s;
  if (exact2_expansion(a, b, c, &s)) {
    return s;
  }
  return exact2(a, b, c);
}

static int exact3(GSHGEOM_Point a, GSHGEOM_Point b, GSHGEOM_Point c,
    GSHGEOM_Point d) {
  GSHGEOM_ExactTerm t[24];
  size_t n = 0;
  put_det3(t, &n, 0, b, c, d);
  put_det3(t, &n, 1, a, c, d);
  put_det3(t, &n, 0, a, b, d);
  put_det3(t, &n, 1, a, b, c);
  return gshgeom_exact_sign(t, n);
}

static int finite2(GSHGEOM_Point p) {
  return is_finite(p.x) && is_finite(p.y);
}

static int finite3(GSHGEOM_Point p) {
  return is_finite(p.x) && is_finite(p.y) && is_finite(p.z);
}

GSHGEOM_Result gshgeom_orient2d_exact(
    GSHGEOM_Point a, GSHGEOM_Point b, GSHGEOM_Point c, int * sign) {
  if (!sign || !finite2(a) || !finite2(b) || !finite2(c)) {
    return GSHGEOM_ERR_INVALID;
  }
  *sign = exact2(a, b, c);
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_orient2d(
    GSHGEOM_Point a, GSHGEOM_Point b, GSHGEOM_Point c, int * sign) {
  if (!sign) {
    return GSHGEOM_ERR_INVALID;
  }
  if (all_small2(a, b, c)) {
    int s;
    /* Within the limit is finite, so the checks below would pass: the filter is all that is left of them. */
    if (!filter2_small(a, b, c, &s)) {
      s = exact2_fast(a, b, c);
    }
    *sign = s;
    return GSHGEOM_OK;
  }
  if (!finite2(a) || !finite2(b) || !finite2(c)) {
    return GSHGEOM_ERR_INVALID;
  }
  /* A coordinate past the limit: the filter's bound does not hold there, and the exact form decides. */
  *sign = exact2_fast(a, b, c);
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_orient3d_exact(GSHGEOM_Point a, GSHGEOM_Point b,
    GSHGEOM_Point c, GSHGEOM_Point d, int * sign) {
  if (!sign || !finite3(a) || !finite3(b) || !finite3(c) || !finite3(d)) {
    return GSHGEOM_ERR_INVALID;
  }
  *sign = exact3(a, b, c, d);
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_orient3d(GSHGEOM_Point a, GSHGEOM_Point b,
    GSHGEOM_Point c, GSHGEOM_Point d, int * sign) {
  if (!sign || !finite3(a) || !finite3(b) || !finite3(c) || !finite3(d)) {
    return GSHGEOM_ERR_INVALID;
  }
  int s;
  if (!filter3(a, b, c, d, &s)) {
    s = exact3(a, b, c, d);
  }
  *sign = s;
  return GSHGEOM_OK;
}

GSHGEOM_Point gshgeom_normal_left(GSHGEOM_Point tangent) {
  /* 0.0 - v rather than -v: the sign of a zero is +0 for either zero. */
  GSHGEOM_Point n = {tangent.y, 0.0 - tangent.x, 0.0};
  return n;
}

GSHGEOM_Point gshgeom_normal_right(GSHGEOM_Point tangent) {
  GSHGEOM_Point n = {0.0 - tangent.y, tangent.x, 0.0};
  return n;
}

GSHGEOM_Result gshgeom_polygon_signed_area(
    const GSHGEOM_Point * points, size_t count, double * out) {
  if (!out || (count && !points)) {
    return GSHGEOM_ERR_INVALID;
  }
  for (size_t i = 0; i < count; i++) {
    if (!finite2(points[i])) {
      return GSHGEOM_ERR_INVALID;
    }
  }
  double sum = 0.0;
  for (size_t i = 1; i + 1 < count; i++) {
    double ux = points[i].x - points[0].x;
    double uy = points[i].y - points[0].y;
    double vx = points[i + 1].x - points[0].x;
    double vy = points[i + 1].y - points[0].y;
    sum += ux * vy - uy * vx;
  }
  double area = sum * 0.5;
  if (!is_finite(area)) {
    return GSHGEOM_ERR_LIMIT;
  }
  *out = area;
  return GSHGEOM_OK;
}
