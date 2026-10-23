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
 * The lowering of a path to lines and cubics (stroke_curves.h, SG40).
 *
 * ## A conic arc
 *
 * A rational quadratic with control points `P0 P1 P2` and weights `(w0, w1, w2)` is
 * the same curve, reparametrized, as the one with weights `(1, w, 1)`, `w = w1 /
 * sqrt(w0 w2)`.  For `w < 1` it is an ellipse arc: the affine image `A` of the
 * circular arc of radius 1 and angle `psi = 2 acos w` whose control triangle `A`
 * maps onto `P0 P1 P2`.  Write `c = cos(psi / 2) = w`, `s = sin(psi / 2)`,
 * `T = tan(psi / 4)`; the cubic with the same ends and end tangents whose handles
 * are `k = (4/3) T` of the radius (the fraction `lambda = 4 w / (3 (1 + w))` of
 * the control legs `|P0 P1|`, `|P2 P1|`, since a leg is `tan(psi / 2)` of the
 * radius) is, on the unit circle with the arc centred on the x axis,
 *
 *     Q(t) = (1-t)^3 (c, -s) + 3 (1-t)^2 t (c + k s, -s + k c)
 *          + 3 (1-t) t^2 (c + k s, s - k c) + t^3 (c, s).
 *
 * `|Q(t)|^2 - 1` is a polynomial of degree 6 with a double root at `t = 0`, at
 * `t = 1` (the ends are on the circle and the tangent there is the circle's) and at
 * `t = 1/2` (`Q(1/2)` is at distance `c + (3/4) k s = 1` by the choice of `k`, and the
 * polynomial is symmetric), so it is `kappa t^2 (1-t)^2 (1-2t)^2`; its `t^6`
 * coefficient `|a3|^2` with `a3 = -P0 + 3 Q1 - 3 Q2 + P3 = (0, 6 k c - 4 s)` gives
 * `kappa = (3 k c - 2 s)^2`, and `3 k c - 2 s = -4 T^3 / (1 + T^2)`.  The maximum
 * of `t^2 (1-t)^2 (1-2t)^2` is `1/108` (at `t (1-t) = 1/6`), so the cubic stays
 * outside the circle and within
 *
 *     e = sqrt(1 + 16 T^6 / (108 (1 + T^2)^2)) - 1 = sqrt(1 + (1-w)^3 / (27 (1+w))) - 1
 *       <= (1-w)^3 / (54 (1+w))
 *
 * of it in the radial direction (`T^2 = (1-w) / (1+w)`), and, the cubic's angle
 * being monotone, the arc is within the same of the cubic: the Hausdorff distance
 * is at most `e` (a quarter turn: 2.7e-4; 60 degrees: 2.4e-5; it falls as
 * `psi^6`).  Under the affine map `A` it is at most `|A| e`, and `|A|` is
 * bounded from the control triangle: `A` takes the unit along the arc's
 * bisector to `(P1 - M) c / s^2` and the unit along the chord to `(P2 - P0) / (2 s)`
 * (`M` the chord's midpoint), so `|A|_F^2 = |P2 - P0|^2 / (4 s^2) + |P1 - M|^2 c^2 / s^4`
 * (and `s^2 = (1 - w)(1 + w)`).  Halving the angle (exact bisection of the conic,
 * a half having the weight `sqrt((1 + w) / 2)`) cuts `e` by about 64.
 *
 * ## A rational curve in general
 *
 * `R(u) = sum w_i P_i B_i(u) / W(u)` and `Q(u) = sum P_i B_i(u)` differ by
 * `R - Q = sum B_i (w_i / W - 1)(P_i - P_0)`, and `W` lies between the least and
 * greatest weight, so `|R - Q| <= rho D` with `rho = (w_max - w_min) / w_min` and
 * `D = max |P_i - P_0|`.  Bisection in homogeneous coordinates halves `D` and
 * the weights' spread together, so the bound falls by 4 each time.  This bound
 * serves a rational cubic and a rational quadratic whose middle weight is above 1.
 *
 * Every computed bound is raised by 2^-20 for the rounding of its own
 * computation; the control points the lowering writes carry the rounding of a
 * handful of operations (relative 2^-52), which stroke_curves.h counts in its
 * reserve.
 */

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/math/basic.h>
#include <ghoti.io/shape-geometry/stroke_curves.h>

#include <string.h>

/** The most halvings of one curve. */
#define LOWER_DEPTH 20

/** How much a computed bound is raised for the rounding of its own computation. */
#define UPL (1.0 + 0x1p-20)

typedef struct Pt {
  double x, y;
} Pt;

static Pt pt(double x, double y) {
  Pt p;
  p.x = x;
  p.y = y;
  return p;
}

static GSHGEOM_Point gp(Pt p) {
  GSHGEOM_Point g = {p.x, p.y, 0.0};
  return g;
}

static double fabs_(double x) {
  return x < 0 ? -x : x;
}

static double len(Pt a, Pt b) {
  double dx = fabs_(a.x - b.x), dy = fabs_(a.y - b.y);
  double m = dx > dy ? dx : dy;
  if (m == 0.0) {
    return 0.0;
  }
  double x = dx / m, y = dy / m;
  return m * gmath_sqrt(x * x + y * y);
}

static Pt lerp(Pt a, Pt b, double t) {
  return pt(a.x + t * (b.x - a.x), a.y + t * (b.y - a.y));
}

static GSHGEOM_Result cubic(GSHGEOM_Path * out, Pt p1, Pt p2, Pt p3) {
  return gshgeom_path_cubic_to(out, gp(p1), gp(p2), gp(p3));
}

/** The conic with ends `a`, `c` and middle point `b`, standard weight `w` (ends 1), as lines of cubics. */
static GSHGEOM_Result conic(GSHGEOM_Path * out, Pt p0, Pt p1, Pt p2, double w, double tol, int depth) {
  if (w == 1.0) {
    /* an ordinary quadratic: the degree elevation is exact */
    return cubic(out, lerp(p0, p1, 2.0 / 3.0), lerp(p2, p1, 2.0 / 3.0), p2);
  }
  double bound;
  if (w < 1.0) {
    double s2 = (1.0 - w) * (1.0 + w);
    Pt m = pt(0.5 * (p0.x + p2.x), 0.5 * (p0.y + p2.y));
    double chord = len(p0, p2), sag = len(p1, m);
    double fro = gmath_sqrt(chord * chord / (4.0 * s2) + sag * sag * w * w / (s2 * s2));
    double e = (1.0 - w) * (1.0 - w) * (1.0 - w) / (54.0 * (1.0 + w));
    bound = fro * e * UPL;
  } else {
    double d = len(p0, p1) > len(p0, p2) ? len(p0, p1) : len(p0, p2);
    bound = (w - 1.0) * d * UPL;
  }
  if (bound <= tol) {
    if (w < 1.0) {
      double lam = 4.0 * w / (3.0 * (1.0 + w));
      return cubic(out, lerp(p0, p1, lam), lerp(p2, p1, lam), p2);
    }
    return cubic(out, lerp(p0, p1, 2.0 / 3.0), lerp(p2, p1, 2.0 / 3.0), p2);
  }
  if (depth >= LOWER_DEPTH) {
    return GSHGEOM_ERR_LIMIT;
  }
  double a = 1.0 + w;
  Pt q1l = pt((p0.x + w * p1.x) / a, (p0.y + w * p1.y) / a);
  Pt q1r = pt((w * p1.x + p2.x) / a, (w * p1.y + p2.y) / a);
  Pt mid = pt((p0.x + 2.0 * w * p1.x + p2.x) / (2.0 * a), (p0.y + 2.0 * w * p1.y + p2.y) / (2.0 * a));
  double wh = gmath_sqrt((1.0 + w) / 2.0);
  GSHGEOM_Result r = conic(out, p0, q1l, mid, wh, tol, depth + 1);
  if (r == GSHGEOM_OK) {
    r = conic(out, mid, q1r, p2, wh, tol, depth + 1);
  }
  return r;
}

/** A homogeneous Bezier piece of degree 2 or 3. */
typedef struct Hom {
  int n;
  double x[4], y[4], w[4]; ///< `w P`, `w`
} Hom;

static void hom_point(const Hom * h, int i, Pt * p) {
  *p = pt(h->x[i] / h->w[i], h->y[i] / h->w[i]);
}

static void hom_split(const Hom * h, Hom * l, Hom * r) {
  double x[4], y[4], w[4];
  memcpy(x, h->x, sizeof x);
  memcpy(y, h->y, sizeof y);
  memcpy(w, h->w, sizeof w);
  l->n = r->n = h->n;
  l->x[0] = x[0];
  l->y[0] = y[0];
  l->w[0] = w[0];
  r->x[h->n] = x[h->n];
  r->y[h->n] = y[h->n];
  r->w[h->n] = w[h->n];
  for (int lev = 1; lev <= h->n; lev++) {
    for (int i = 0; i + lev <= h->n; i++) {
      x[i] = 0.5 * (x[i] + x[i + 1]);
      y[i] = 0.5 * (y[i] + y[i + 1]);
      w[i] = 0.5 * (w[i] + w[i + 1]);
    }
    l->x[lev] = x[0];
    l->y[lev] = y[0];
    l->w[lev] = w[0];
    r->x[h->n - lev] = x[h->n - lev];
    r->y[h->n - lev] = y[h->n - lev];
    r->w[h->n - lev] = w[h->n - lev];
  }
}

static GSHGEOM_Result general(GSHGEOM_Path * out, const Hom * h, double tol, int depth) {
  Pt p[4];
  double wmin = h->w[0], wmax = h->w[0];
  for (int i = 0; i <= h->n; i++) {
    hom_point(h, i, &p[i]);
    wmin = h->w[i] < wmin ? h->w[i] : wmin;
    wmax = h->w[i] > wmax ? h->w[i] : wmax;
  }
  double d = 0.0;
  for (int i = 1; i <= h->n; i++) {
    double l = len(p[0], p[i]);
    d = l > d ? l : d;
  }
  double bound = (wmax - wmin) / wmin * d * UPL;
  if (bound <= tol) {
    if (h->n == 2) {
      return cubic(out, lerp(p[0], p[1], 2.0 / 3.0), lerp(p[2], p[1], 2.0 / 3.0), p[2]);
    }
    return cubic(out, p[1], p[2], p[3]);
  }
  if (depth >= LOWER_DEPTH) {
    return GSHGEOM_ERR_LIMIT;
  }
  Hom l, r;
  hom_split(h, &l, &r);
  GSHGEOM_Result res = general(out, &l, tol, depth + 1);
  if (res == GSHGEOM_OK) {
    res = general(out, &r, tol, depth + 1);
  }
  return res;
}

GSHGEOM_Result gshgeom_stroke_curves_lower_cubic(const GSHGEOM_Path * in, double tolerance,
    const GSHGEOM_Allocator * alloc, const GSHGEOM_Limits * limits, GSHGEOM_Path ** out) {
  if (!in || !out || !(tolerance > 0.0) || !(tolerance - tolerance == 0.0) || gshgeom_path_dimension(in) != 2) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Path * dst = NULL;
  GSHGEOM_Result r = gshgeom_path_create(alloc, limits, 2, &dst);
  size_t n = gshgeom_path_element_count(in);
  for (size_t i = 0; i < n && r == GSHGEOM_OK; i++) {
    GSHGEOM_Element el;
    r = gshgeom_path_element(in, i, &el);
    if (r != GSHGEOM_OK) {
      break;
    }
    const GSHGEOM_Point * p = el.curve.p;
    const double * w = el.curve.w;
    switch (el.verb) {
    case GSHGEOM_VERB_MOVE:
      r = gshgeom_path_move_to(dst, p[0]);
      break;
    case GSHGEOM_VERB_LINE:
      r = gshgeom_path_line_to(dst, p[1]);
      break;
    case GSHGEOM_VERB_CLOSE:
      r = gshgeom_path_close(dst);
      break;
    case GSHGEOM_VERB_CUBIC:
      r = gshgeom_path_cubic_to(dst, p[1], p[2], p[3]);
      break;
    case GSHGEOM_VERB_QUADRATIC:
      r = cubic(dst, lerp(pt(p[0].x, p[0].y), pt(p[1].x, p[1].y), 2.0 / 3.0),
          lerp(pt(p[2].x, p[2].y), pt(p[1].x, p[1].y), 2.0 / 3.0), pt(p[2].x, p[2].y));
      break;
    case GSHGEOM_VERB_RATIONAL_QUADRATIC: {
      double ws = w[1] / gmath_sqrt(w[0] * w[2]);
      r = conic(dst, pt(p[0].x, p[0].y), pt(p[1].x, p[1].y), pt(p[2].x, p[2].y), ws, tolerance, 0);
      break;
    }
    case GSHGEOM_VERB_RATIONAL_CUBIC: {
      Hom h;
      h.n = 3;
      for (int k = 0; k < 4; k++) {
        h.x[k] = w[k] * p[k].x;
        h.y[k] = w[k] * p[k].y;
        h.w[k] = w[k];
      }
      r = general(dst, &h, tolerance, 0);
      break;
    }
    default:
      r = GSHGEOM_ERR_INTERNAL;
      break;
    }
  }
  if (r != GSHGEOM_OK) {
    gshgeom_path_destroy(dst);
    return r;
  }
  *out = dst;
  return GSHGEOM_OK;
}
