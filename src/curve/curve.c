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
 * Curves: construction, validation, evaluation and bounds (curve.h).
 *
 * Nothing here calls the C library's math: evaluation is de Casteljau's
 * algorithm with the operations in a fixed order, and the bounds find the
 * extrema of each coordinate as the roots of a polynomial by recursion on its
 * derivative and bisection.  No operation depends on the host, and
 * -ffp-contract=off keeps the compiler from fusing any of them.
 */

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/shape-geometry/curve.h>

/** A finite double: v - v is zero for every finite v, NaN for NaN and infinity. */
static int is_finite(double v) {
  return (v - v) == 0.0;
}

static double axis_of(const GSHGEOM_Point * p, int axis) {
  return axis == 0 ? p->x : (axis == 1 ? p->y : p->z);
}

int gshgeom_curve_point_count(GSHGEOM_Curve_Kind kind) {
  switch (kind) {
  case GSHGEOM_CURVE_LINE:
    return 2;
  case GSHGEOM_CURVE_QUADRATIC:
  case GSHGEOM_CURVE_RATIONAL_QUADRATIC:
    return 3;
  case GSHGEOM_CURVE_CUBIC:
  case GSHGEOM_CURVE_RATIONAL_CUBIC:
    return 4;
  case GSHGEOM_CURVE_KIND_COUNT:
    break;
  }
  return 0;
}

int gshgeom_curve_is_rational(GSHGEOM_Curve_Kind kind) {
  return kind == GSHGEOM_CURVE_RATIONAL_QUADRATIC ||
      kind == GSHGEOM_CURVE_RATIONAL_CUBIC;
}

GSHGEOM_Result gshgeom_curve_check(const GSHGEOM_Curve * curve) {
  if (!curve) {
    return GSHGEOM_ERR_INVALID;
  }
  /* A kind outside the enum has no points, so a structure filled in by hand
   * with a wrong one is refused here. */
  int n = gshgeom_curve_point_count(curve->kind);
  if (n == 0 || (curve->dim != 2 && curve->dim != 3)) {
    return GSHGEOM_ERR_INVALID;
  }
  int rational = gshgeom_curve_is_rational(curve->kind);

  /* First everything that is wrong whatever its size: a NaN, an infinity, a
   * weight that is not positive.  Then the size.  So a NaN is never reported
   * as a limit because an earlier number was large. */
  int too_big = 0;
  for (int i = 0; i < n; i++) {
    for (int a = 0; a < curve->dim; a++) {
      double v = axis_of(&curve->p[i], a);
      if (!is_finite(v)) {
        return GSHGEOM_ERR_INVALID;
      }
      if (v > GSHGEOM_COORD_LIMIT || v < -GSHGEOM_COORD_LIMIT) {
        too_big = 1;
      }
    }
    if (rational) {
      double w = curve->w[i];
      if (!is_finite(w) || !(w > 0.0)) {
        return GSHGEOM_ERR_INVALID;
      }
      if (w < GSHGEOM_WEIGHT_MIN || w > GSHGEOM_WEIGHT_MAX) {
        too_big = 1;
      }
    }
  }
  return too_big ? GSHGEOM_ERR_LIMIT : GSHGEOM_OK;
}

/* Fill a curve with the used slots only (the unused are zero, and z is zero
 * in 2D), check it, and hand it over on success. */
static GSHGEOM_Result make_curve(GSHGEOM_Curve_Kind kind, int dim,
    const GSHGEOM_Point * pts, const double * w, GSHGEOM_Curve * out) {
  if (!out || (gshgeom_curve_is_rational(kind) && !w)) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Curve curve = {0};
  int n = gshgeom_curve_point_count(kind);
  curve.kind = kind;
  curve.dim = dim;
  for (int i = 0; i < n; i++) {
    curve.p[i] = pts[i];
    if (dim == 2) {
      curve.p[i].z = 0.0;
    }
    if (w) {
      curve.w[i] = w[i];
    }
  }
  GSHGEOM_Result r = gshgeom_curve_check(&curve);
  if (r == GSHGEOM_OK) {
    *out = curve;
  }
  return r;
}

GSHGEOM_Result gshgeom_curve_line(
    int dim, GSHGEOM_Point a, GSHGEOM_Point b, GSHGEOM_Curve * out) {
  GSHGEOM_Point pts[2] = {a, b};
  return make_curve(GSHGEOM_CURVE_LINE, dim, pts, NULL, out);
}

GSHGEOM_Result gshgeom_curve_quadratic(int dim, GSHGEOM_Point a,
    GSHGEOM_Point c, GSHGEOM_Point b, GSHGEOM_Curve * out) {
  GSHGEOM_Point pts[3] = {a, c, b};
  return make_curve(GSHGEOM_CURVE_QUADRATIC, dim, pts, NULL, out);
}

GSHGEOM_Result gshgeom_curve_cubic(int dim, GSHGEOM_Point a, GSHGEOM_Point c1,
    GSHGEOM_Point c2, GSHGEOM_Point b, GSHGEOM_Curve * out) {
  GSHGEOM_Point pts[4] = {a, c1, c2, b};
  return make_curve(GSHGEOM_CURVE_CUBIC, dim, pts, NULL, out);
}

GSHGEOM_Result gshgeom_curve_rational_quadratic(int dim, GSHGEOM_Point a,
    GSHGEOM_Point c, GSHGEOM_Point b, const double w[3], GSHGEOM_Curve * out) {
  GSHGEOM_Point pts[3] = {a, c, b};
  return make_curve(GSHGEOM_CURVE_RATIONAL_QUADRATIC, dim, pts, w, out);
}

GSHGEOM_Result gshgeom_curve_rational_cubic(int dim, GSHGEOM_Point a,
    GSHGEOM_Point c1, GSHGEOM_Point c2, GSHGEOM_Point b, const double w[4],
    GSHGEOM_Curve * out) {
  GSHGEOM_Point pts[4] = {a, c1, c2, b};
  return make_curve(GSHGEOM_CURVE_RATIONAL_CUBIC, dim, pts, w, out);
}

/* ---------------------------------------------------------------------- */
/* Evaluation                                                              */
/* ---------------------------------------------------------------------- */

/** A checked curve in arrays, so that an axis is an index. */
typedef struct Work {
  int n;        /* control points */
  int dim;      /* axes in use */
  int rational; /* weights in use */
  double c[4][3];
  double w[4];  /* 1 for a polynomial curve */
} Work;

static void load(const GSHGEOM_Curve * curve, Work * work) {
  work->n = gshgeom_curve_point_count(curve->kind);
  work->dim = curve->dim;
  work->rational = gshgeom_curve_is_rational(curve->kind);
  for (int i = 0; i < work->n; i++) {
    for (int a = 0; a < 3; a++) {
      work->c[i][a] = (a < work->dim) ? axis_of(&curve->p[i], a) : 0.0;
    }
    work->w[i] = work->rational ? curve->w[i] : 1.0;
  }
}

/* The point at t, by de Casteljau; the ends exactly. */
static void eval_work(const Work * work, double t, double out[3]) {
  int n = work->n;
  if (t == 0.0 || t == 1.0) {
    int end = (t == 0.0) ? 0 : n - 1;
    for (int a = 0; a < 3; a++) {
      out[a] = work->c[end][a];
    }
    return;
  }
  double s = 1.0 - t;
  double h[4][3];
  double hw[4];
  for (int i = 0; i < n; i++) {
    hw[i] = work->w[i];
    for (int a = 0; a < work->dim; a++) {
      /* Homogeneous coordinates; for a polynomial curve the weight is 1 and
       * the product is exact. */
      h[i][a] = work->rational ? work->w[i] * work->c[i][a] : work->c[i][a];
    }
  }
  for (int level = n - 1; level > 0; level--) {
    for (int i = 0; i < level; i++) {
      for (int a = 0; a < work->dim; a++) {
        h[i][a] = s * h[i][a] + t * h[i + 1][a];
      }
      if (work->rational) {
        hw[i] = s * hw[i] + t * hw[i + 1];
      }
    }
  }
  for (int a = 0; a < 3; a++) {
    if (a >= work->dim) {
      out[a] = 0.0;
    } else {
      out[a] = work->rational ? h[0][a] / hw[0] : h[0][a];
    }
  }
}

GSHGEOM_Result gshgeom_curve_eval(
    const GSHGEOM_Curve * curve, double t, GSHGEOM_Point * out) {
  GSHGEOM_Result r = gshgeom_curve_check(curve);
  if (r != GSHGEOM_OK) {
    return r;
  }
  /* Written so that NaN fails, and infinity too. */
  if (!out || !(t >= 0.0 && t <= 1.0)) {
    return GSHGEOM_ERR_INVALID;
  }
  Work work;
  load(curve, &work);
  double v[3];
  eval_work(&work, t, v);
  out->x = v[0];
  out->y = v[1];
  out->z = v[2];
  return GSHGEOM_OK;
}

/* ---------------------------------------------------------------------- */
/* Bounds                                                                  */
/* ---------------------------------------------------------------------- */

/* The Bernstein polynomial c[0..m] (c[k] is the coefficient of B_k^m) at t, by
 * de Casteljau's algorithm. */
static double bernstein_eval(const double * c, int m, double t) {
  double s = 1.0 - t;
  double v[8];
  for (int i = 0; i <= m; i++) {
    v[i] = c[i];
  }
  for (int level = m; level > 0; level--) {
    for (int i = 0; i < level; i++) {
      v[i] = s * v[i] + t * v[i + 1];
    }
  }
  return v[0];
}

/* A root of the polynomial in (a, b), given that it changes sign there: the
 * value at a is fa.  Bisection until the interval cannot be halved in double
 * precision (at most 200 steps, which is 2^-200 of the unit interval). */
static double bisect(const double * c, int m, double a, double b, double fa) {
  for (int step = 0; step < 200; step++) {
    double mid = a + (b - a) * 0.5;
    if (!(mid > a && mid < b)) {
      break;
    }
    double fm = bernstein_eval(c, m, mid);
    if (fm == 0.0) {
      return mid;
    }
    if ((fm < 0.0) == (fa < 0.0)) {
      a = mid;
      fa = fm;
    } else {
      b = mid;
    }
  }
  return a + (b - a) * 0.5;
}

/* The roots in [0, 1] of the Bernstein polynomial c[0..m] (m >= 1) at which it
 * changes sign or is zero, ascending, into r (room for 8).  The roots of the
 * derivative cut [0, 1] into intervals on which the polynomial is monotone, so
 * each holds at most one.  The derivative's coefficients are m times the
 * differences of the polynomial's. */
static int roots01(const double * c, int m, double * r) {
  if (m == 1) {
    if ((c[0] < 0.0 && c[1] > 0.0) || (c[0] > 0.0 && c[1] < 0.0)) {
      r[0] = c[0] / (c[0] - c[1]);
      return 1;
    }
    return 0;
  }
  double d[8] = {0};
  for (int k = 0; k < m; k++) {
    d[k] = (double)m * (c[k + 1] - c[k]);
  }
  double crit[8] = {0};
  int nc = roots01(d, m - 1, crit);

  double pts[10];
  int k = 0;
  pts[k++] = 0.0;
  for (int i = 0; i < nc; i++) {
    if (crit[i] > 0.0 && crit[i] < 1.0) {
      pts[k++] = crit[i];
    }
  }
  pts[k++] = 1.0;

  int count = 0;
  double fa = bernstein_eval(c, m, pts[0]);
  for (int i = 0; i + 1 < k; i++) {
    double a = pts[i];
    double b = pts[i + 1];
    double fb = bernstein_eval(c, m, b);
    if (fa == 0.0) {
      r[count++] = a;
    } else if (fb != 0.0 && ((fa < 0.0) != (fb < 0.0))) {
      r[count++] = bisect(c, m, a, b, fa);
    }
    fa = fb;
  }
  if (fa == 0.0) {
    r[count++] = 1.0;
  }
  return count;
}

/* Binomial coefficients C(n, k) for n <= 4. */
static double choose(int n, int k) {
  static const double table[5][5] = {{1, 0, 0, 0, 0}, {1, 1, 0, 0, 0},
      {1, 2, 1, 0, 0}, {1, 3, 3, 1, 0}, {1, 4, 6, 4, 1}};
  return table[n][k];
}

/* The parameters in [0, 1] where coordinate `axis` may have an extremum: the
 * roots of N'D - ND', where the coordinate is N/D (N = sum w_i x_i B_i,
 * D = sum w_i B_i, degree n; D = 1 for a polynomial curve).
 *
 * By de Casteljau's recurrence N = (1-t) N0 + t N1 and D = (1-t) A0 + t A1,
 * where N0, A0 are the sums over the first n control points and N1, A1 over
 * the last n, each of degree n - 1.  N' = n (N1 - N0) and D' = n (A1 - A0), and
 * the terms in t cancel, leaving
 *
 *   N'D - ND' = n (A0 N1 - A1 N0),
 *
 * a polynomial of degree m = 2n - 2.  In the Bernstein basis, with
 * B_i^{n-1} B_j^{n-1} = C(n-1,i) C(n-1,j) / C(m,i+j) B_{i+j}^m, its coefficient
 * of B_k^m is the sum over i + j = k (i, j < n) of
 *
 *   C(n-1,i) C(n-1,j) / C(m,k) w_i w_{j+1} (x_{j+1} - x_i),
 *
 * weights times a coordinate difference: nothing large is subtracted from
 * anything large, as it would be in the power basis when the weights are 2^60
 * apart.  The positive factor n does not move a root and is left out. */
static int axis_extrema(const Work * work, int axis, double * roots) {
  int n = work->n - 1;
  int m = 2 * n - 2;
  double g[8] = {0};
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      double term = work->w[i] * work->w[j + 1] *
          (work->c[j + 1][axis] - work->c[i][axis]);
      g[i + j] += (choose(n - 1, i) * choose(n - 1, j) / choose(m, i + j)) * term;
    }
  }
  return roots01(g, m, roots);
}

GSHGEOM_Result gshgeom_curve_bounds(
    const GSHGEOM_Curve * curve, GSHGEOM_Bounds * out) {
  GSHGEOM_Result r = gshgeom_curve_check(curve);
  if (r != GSHGEOM_OK) {
    return r;
  }
  if (!out) {
    return GSHGEOM_ERR_INVALID;
  }
  Work work;
  load(curve, &work);

  double t[2 + 3 * 8];
  int count = 0;
  t[count++] = 0.0;
  t[count++] = 1.0;
  if (work.n > 2) {
    for (int a = 0; a < work.dim; a++) {
      count += axis_extrema(&work, a, &t[count]);
    }
  }

  double lo[3];
  double hi[3];
  eval_work(&work, 0.0, lo);
  eval_work(&work, 0.0, hi);
  for (int i = 1; i < count; i++) {
    double v[3];
    eval_work(&work, t[i], v);
    for (int a = 0; a < 3; a++) {
      if (v[a] < lo[a]) {
        lo[a] = v[a];
      }
      if (v[a] > hi[a]) {
        hi[a] = v[a];
      }
    }
  }
  out->empty = 0;
  out->min.x = lo[0];
  out->min.y = lo[1];
  out->min.z = lo[2];
  out->max.x = hi[0];
  out->max.y = hi[1];
  out->max.z = hi[2];
  return GSHGEOM_OK;
}
