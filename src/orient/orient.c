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

/** 1 and *sign set if the filter decides. */
static int filter2(GSHGEOM_Point a, GSHGEOM_Point b, GSHGEOM_Point c, int * sign) {
  if (!(small_enough(a.x) && small_enough(a.y) && small_enough(b.x) &&
          small_enough(b.y) && small_enough(c.x) && small_enough(c.y))) {
    return 0;
  }
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
  if (!sign || !finite2(a) || !finite2(b) || !finite2(c)) {
    return GSHGEOM_ERR_INVALID;
  }
  int s;
  if (!filter2(a, b, c, &s)) {
    s = exact2(a, b, c);
  }
  *sign = s;
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
