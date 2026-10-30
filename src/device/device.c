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
 * The rounding rule, the device path (builder, readback, hash and dump) and the
 * conversion of a model path under a map into it (device.h).
 *
 * Elements store their points in one array, in readback order: control points
 * then the end point, and a close stores the subpath's start.  The builders
 * reserve everything an operation needs before they store anything, so a
 * failure changes capacity and nothing else.
 */

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/math/basic.h>
#include <ghoti.io/shape-geometry/arc.h>
#include <ghoti.io/shape-geometry/device.h>
#include <ghoti.io/shape-geometry/flatten.h>
#include <ghoti.io/shape-geometry/hash.h>

#include "../core/grow.h"
#include "conic_internal.h"

#include <string.h>

/* One pixel in device units: the factor of the rule (an exact power of two). */
#define SCALE ((double)GSHGEOM_DEVICE_ONE)

/* The bounds of the rule, on t = x * 64 (exact): floor(t + 0.5) fits an int32
 * exactly when -2^31 - 0.5 <= t < 2^31 - 0.5. */
#define ROUND_T_MAX 2147483647.5
#define ROUND_T_MIN (-2147483648.5)

/**
 * floor(t + 0.5) for a finite, already scaled t, without ever adding 0.5: the
 * floor and the fraction are both exact, and the half-way point is a
 * comparison.  Infinite t (a finite x so large that x * 64 overflowed) falls to
 * the range check.
 */
static GSHGEOM_Result round_scaled(double t, int32_t * out) {
  if (!(t < ROUND_T_MAX) || t < ROUND_T_MIN) {
    return GSHGEOM_ERR_LIMIT;
  }
  double f = gmath_floor(t);
  double r = (t - f) >= 0.5 ? f + 1.0 : f;
  *out = (int32_t)r;
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_round_26_6(double x, int32_t * out) {
  if (!out || !((x - x) == 0.0)) {
    return GSHGEOM_ERR_INVALID;
  }
  return round_scaled(x * SCALE, out);
}

GSHGEOM_Result gshgeom_scale_26_6(
    int32_t p, int32_t scale, int flip, int32_t * out) {
  if (!out || (flip != 0 && flip != 1)) {
    return GSHGEOM_ERR_INVALID;
  }
  /* |p * scale| <= 2^62: the product and its negation fit an int64_t. */
  int64_t v = (int64_t)p * (int64_t)scale;
  if (flip) {
    v = -v;
  }
  /* floor((v + 2^15) / 2^16) by shifts of a magnitude, never of a negative
   * number (whose right shift C leaves to the implementation). */
  int64_t n = v + 32768;
  int64_t q = n >= 0 ? (int64_t)((uint64_t)n >> 16)
                     : -(int64_t)(((uint64_t)(-n) + 65535u) >> 16);
  if (q > INT32_MAX || q < INT32_MIN) {
    return GSHGEOM_ERR_LIMIT;
  }
  *out = (int32_t)q;
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_device_point(const GSHGEOM_Map * map, GSHGEOM_Point p,
    GSHGEOM_Device_Point * out) {
  if (!map || !out) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Result r = gshgeom_map_check(map);
  if (r != GSHGEOM_OK) {
    return r;
  }
  if (map->dim != 2) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Point m;
  r = gshgeom_map_apply_point(map, p, &m);
  if (r != GSHGEOM_OK) {
    return r;
  }
  GSHGEOM_Device_Point d;
  r = gshgeom_round_26_6(m.x, &d.x);
  if (r == GSHGEOM_OK) {
    r = gshgeom_round_26_6(m.y, &d.y);
  }
  if (r == GSHGEOM_OK) {
    *out = d;
  }
  return r;
}

/* ---------------------------------------------------------------------- */
/* The device path                                                         */
/* ---------------------------------------------------------------------- */

typedef struct Elem {
  uint8_t verb;
  size_t point; /* index of its first stored point */
} Elem;

struct GSHGEOM_Device_Path {
  GSHGEOM_Allocator alloc;
  size_t max_elements;
  Elem * el;
  size_t el_n, el_cap;
  GSHGEOM_Device_Point * pt;
  size_t pt_n, pt_cap;
  int has_current; /* a move has been made */
  int closed;      /* the last element is a close */
  GSHGEOM_Device_Point start;
};

static unsigned verb_points(int verb) {
  switch (verb) {
  case GSHGEOM_DEVICE_MOVE:
  case GSHGEOM_DEVICE_LINE:
  case GSHGEOM_DEVICE_CLOSE:
    return 1;
  case GSHGEOM_DEVICE_QUADRATIC:
    return 2;
  case GSHGEOM_DEVICE_CUBIC:
    return 3;
  default:
    return 0;
  }
}

GSHGEOM_Result gshgeom_device_path_create(const GSHGEOM_Allocator * alloc,
    const GSHGEOM_Limits * limits, GSHGEOM_Device_Path ** out) {
  if (!out) {
    return GSHGEOM_ERR_INVALID;
  }
  if (!alloc) {
    alloc = gshgeom_allocator_default();
  }
  GSHGEOM_Device_Path * path = gcu_allocator_calloc(alloc, 1, sizeof *path);
  if (!path) {
    return GSHGEOM_ERR_OOM;
  }
  path->alloc = *alloc;
  path->max_elements = (limits && limits->max_elements)
      ? limits->max_elements
      : GSHGEOM_LIMITS_DEFAULT_MAX_ELEMENTS;
  *out = path;
  return GSHGEOM_OK;
}

void gshgeom_device_path_destroy(GSHGEOM_Device_Path * path) {
  if (!path) {
    return;
  }
  GSHGEOM_Allocator alloc = path->alloc;
  gcu_allocator_free(&alloc, path->el);
  gcu_allocator_free(&alloc, path->pt);
  gcu_allocator_free(&alloc, path);
}

size_t gshgeom_device_path_element_count(const GSHGEOM_Device_Path * path) {
  return path ? path->el_n : 0;
}

/* Add `verb` with its points, preceded by a move to the closed subpath's start
 * when `with_move` (a drawing verb after a close).  All or nothing. */
static GSHGEOM_Result append(GSHGEOM_Device_Path * path, int verb,
    const GSHGEOM_Device_Point * pts, int with_move) {
  size_t np = verb_points(verb);
  size_t add_el = with_move ? 2 : 1;
  size_t add_pt = np + (with_move ? 1 : 0);
  if (path->el_n > path->max_elements ||
      add_el > path->max_elements - path->el_n) {
    return GSHGEOM_ERR_LIMIT;
  }
  void * grown = path->el;
  GSHGEOM_Result r = gshgeom_grow(&path->alloc, &grown, &path->el_cap,
      path->el_n + add_el, sizeof *path->el);
  if (r != GSHGEOM_OK) {
    return r;
  }
  path->el = grown;
  grown = path->pt;
  r = gshgeom_grow(&path->alloc, &grown, &path->pt_cap, path->pt_n + add_pt,
      sizeof *path->pt);
  if (r != GSHGEOM_OK) {
    return r;
  }
  path->pt = grown;
  if (with_move) {
    path->el[path->el_n].verb = GSHGEOM_DEVICE_MOVE;
    path->el[path->el_n].point = path->pt_n;
    path->pt[path->pt_n] = path->start;
    path->el_n++;
    path->pt_n++;
  }
  path->el[path->el_n].verb = (uint8_t)verb;
  path->el[path->el_n].point = path->pt_n;
  for (size_t i = 0; i < np; i++) {
    path->pt[path->pt_n + i] = pts[i];
  }
  path->el_n++;
  path->pt_n += np;
  path->closed = verb == GSHGEOM_DEVICE_CLOSE;
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_device_path_move_to(
    GSHGEOM_Device_Path * path, GSHGEOM_Device_Point p) {
  if (!path) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Result r = append(path, GSHGEOM_DEVICE_MOVE, &p, 0);
  if (r == GSHGEOM_OK) {
    path->has_current = 1;
    path->start = p;
  }
  return r;
}

GSHGEOM_Result gshgeom_device_path_line_to(
    GSHGEOM_Device_Path * path, GSHGEOM_Device_Point p) {
  if (!path || !path->has_current) {
    return GSHGEOM_ERR_INVALID;
  }
  return append(path, GSHGEOM_DEVICE_LINE, &p, path->closed);
}

GSHGEOM_Result gshgeom_device_path_quadratic_to(GSHGEOM_Device_Path * path,
    GSHGEOM_Device_Point c, GSHGEOM_Device_Point p) {
  if (!path || !path->has_current) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Device_Point pts[2] = {c, p};
  return append(path, GSHGEOM_DEVICE_QUADRATIC, pts, path->closed);
}

GSHGEOM_Result gshgeom_device_path_cubic_to(GSHGEOM_Device_Path * path,
    GSHGEOM_Device_Point c1, GSHGEOM_Device_Point c2, GSHGEOM_Device_Point p) {
  if (!path || !path->has_current) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Device_Point pts[3] = {c1, c2, p};
  return append(path, GSHGEOM_DEVICE_CUBIC, pts, path->closed);
}

GSHGEOM_Result gshgeom_device_path_close(GSHGEOM_Device_Path * path) {
  if (!path || !path->has_current) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Device_Point start = path->start;
  return append(path, GSHGEOM_DEVICE_CLOSE, &start, 0);
}

GSHGEOM_Result gshgeom_device_path_element(const GSHGEOM_Device_Path * path,
    size_t index, GSHGEOM_Device_Element * out) {
  if (!path || !out || index >= path->el_n) {
    return GSHGEOM_ERR_INVALID;
  }
  const Elem * e = &path->el[index];
  GSHGEOM_Device_Element el;
  memset(&el, 0, sizeof el);
  el.verb = (GSHGEOM_Device_Verb)e->verb;
  el.count = verb_points(e->verb);
  for (unsigned i = 0; i < el.count; i++) {
    el.p[i] = path->pt[e->point + i];
  }
  *out = el;
  return GSHGEOM_OK;
}

/* ---------------------------------------------------------------------- */
/* A model path under a map                                                */
/* ---------------------------------------------------------------------- */

/*
 * ## The conic route: the definition (flatten.h, revision 3)
 *
 * Every rational curve of a model path is treated by this one route: a
 * rational quadratic, a rational cubic and each piece of an arc.  The control
 * points are mapped in `double` by ::gshgeom_map_apply_point and scaled by 64
 * (exact), so they are the curve's definition in 26.6 units (X_i); the weights
 * W_i are unchanged (an affine map keeps them).  The homogeneous curve
 * (W_i X_i, W_i) is split at the parameter 1/2 by de Casteljau in `double`
 * until the piece is flat, which is when every interior control point of the
 * projected piece is within `delta = t - 1` units of the piece's chord, a
 * segment; a flat piece becomes a line to its end, rounded once to 26.6 (to
 * nearest, ties upward).  The ends of the whole curve are the mapped control
 * points themselves, rounded once, so the device ends of an arc are its
 * authoritative ends and a closed arc closes in device space.  A split point is
 * the quotient of its homogeneous point, clamped into the box of the control
 * points (which holds the exact point): the same value ends the left piece and
 * starts the right one.  The pieces are visited left to right, and a piece
 * deeper than the cap is ended as it stands.
 *
 * The tolerance `t` is the flattener's (flatten.h), at least
 * ::GSHGEOM_FLATTEN_MIN_TOLERANCE and the default ::GSHGEOM_FLATTEN_TOLERANCE
 * for the calls without one; `t` is the whole budget and `delta` is `t - 1`
 * because the vertex rounding takes 0.7072 of the unit that remains.  The cap
 * is ::GSHGEOM_CONIC_MAX_DEPTH, 240.  Everything is one `double` operation in
 * the written order: no `long double`, no fused multiply-add, no library
 * function (`-ffp-contract=off`), so the same bits on every target.
 *
 * ## The theorem
 *
 * Let the curve be a rational Bezier of degree 2 or 3 whose control points X_i
 * (26.6 units, after the map) each round into `int32` (so |coordinate| <=
 * 2^31 + 1/2) and whose weights W_i are in [2^-100, 2^100] (the range the path
 * accepts), and let `t >= 2`.  Then every point of the exact curve (the
 * rational Bezier of the X_i and W_i in real arithmetic) is within `t - 0.29`
 * units of the emitted polyline (so within the `t - 0.12` of the integer
 * flattener's headline), and the depth cap never ends a piece that is not flat.
 *
 * ## The proof
 *
 * Units are 26.6 units.  u = 2^-53.  Xmax = 2^31 + 1 bounds every coordinate
 * of every control point, and hence of every point of the hull.  "Exact" means
 * real arithmetic on the doubles X_i and W_i.
 *
 * **(a) The hull.**  A rational Bezier curve with positive weights is the
 * convex combination sum_i lambda_i(s) X_i with lambda_i = b_i W_i / sum_j
 * b_j W_j >= 0.  The piece of it over a parameter interval [a, b] is again a
 * rational Bezier curve, on the exact control points P_i that de Casteljau (or
 * the blossom) gives, whose weights are positive (convex combinations of the
 * W_i).  So the exact piece lies in the convex hull of its exact projected
 * control points P_i, and, by the same token, so do the P_i in the box of the
 * X_i.
 *
 * **(b) The computed pieces are near the exact ones.**  Let k <= 240 be the
 * number of halvings.  Every new homogeneous point is the average of two
 * parents, at most three averages deep (a cubic's middle point; a quadratic's
 * is two deep), each computed as 0.5 * (a + b), the sum rounded and the
 * halving exact.  Write e_w and e_x for the error of a computed weight and
 * of a computed `W X` against the exact ones.  The root's weights are exact,
 * and its `W X` is rounded once: |e_x| <= u W Xmax.  An average of two points
 * averages their errors, and a rounded sum of terms of one sign adds a relative
 * error of at most u (and, for `W X`, of at most u (|W_a X_a| + |W_b X_b|) / 2
 * <= u W Xmax in absolute terms, W being the average's weight).  So after
 * one halving, with the relative-to-weight errors c_w = |e_w| / (u W) and
 * c_x = |e_x| / (u W Xmax), c_w grows by at most 3 and c_x by at most 3, from 0
 * and 1: c_w <= 3k and c_x <= 1 + 3k, each times 1.01 for the second-order
 * terms.  The projection is (W X + e_x) / (W + e_w) against W X / W, so each
 * coordinate is off by at most (c_x + c_w) u Xmax (1 + 10^-9) plus the
 * division's rounding u Xmax: at k = 240, (1 + 720 + 720) 1.01 + 1 <= 1457
 * u Xmax = 3.47 * 10^-4 per axis, and E <= 2^-10 units (0.00098) as a
 * distance.  The same bound holds for a clamped split point, which is nearer
 * to the exact one than the unclamped.  An underflow adds at most 2^-900
 * (the products W X cannot overflow: |W X| < 2^132).  So every computed
 * projected control point Q_i, and every computed vertex, is within E = 2^-10
 * of the exact one at every depth up to the cap.
 *
 * **(c) The test measures what it measures.**  The distance of Q_i to the
 * chord is computed from differences of values below 2^32 in magnitude, a
 * parameter clamped to [0, 1] and a squared comparison with `delta^2`; each
 * operation is one rounding, and the result is within 7 u |Q_i - S| <= 2^-16
 * units of the true distance (an operand pair with a chord shorter than 2^-500
 * adds at most that length, since the parameter is clamped).  So a piece the test
 * calls flat has every Q_i within delta + 2^-16 of its chord.
 *
 * **(d) The vertex rounding** moves each end of a chord by at most 1/2 unit per
 * axis, so any point of the chord by at most sqrt(2)/2 <= 0.7072 units.
 *
 * **(e) Together.**  For a point of the exact curve, take the piece that
 * contains its parameter.  By (a) the point is in the hull of that piece's exact
 * P_i.  The P_i are within E of the computed Q_i by (b), and the Q_i are within
 * delta + 2^-16 of the computed chord S-S' (c) (and the ends P_0 and P_n are
 * within E of S and S', by (b)), so the hull is within delta + 2^-16 + E of the
 * segment S-S' (the neighbourhood of a convex set is convex), and the emitted
 * segment, between the rounded S and S' (d), is 0.7072 further:
 * `delta + 2^-16 + 2^-10 + 0.7072 = t - 1 + 0.7072 + 0.00099 = t - 0.2918`.
 * QED, given that the cap never ends a piece that is not flat:
 *
 * **(f) The cap never binds.**  Take a piece of exact control points P_i from
 * a parameter interval of length h = 2^-k, with the translation X_0 taken as
 * the origin in the proof (nothing in the code depends on it): every
 * |X_j - X_0|, and every projected point of the hull, is at most B = 2^32 + 1.
 * The Bernstein points of the piece's homogeneous polynomial F(s) = (W X, W)
 * on [a, a + h] are values of its blossom f, and two neighbours differ by
 * f(a.., b, ..) - f(a.., a, ..) = h g(..), where g is the blossom of F' / m
 * (m the degree), at arguments in [0, 1]: a convex combination of the
 * Bernstein coefficients of F' / m, which are the differences of the original
 * homogeneous points.  So neighbouring homogeneous points of the piece differ,
 * in each of the three components, by at most h times the largest difference
 * of two neighbours of the original: |Delta (W X)| <= h (W_max B + W_max B) =
 * h 2^101 B and |Delta W| <= h 2^100.  The projected points then differ by
 * Delta Q = (Delta (W X) - Q Delta W) / W' with W' >= 2^-100, so by at most
 * (2^101 B + B 2^100) h 2^100 = 3 B 2^200 h per axis.  An interior control
 * point of the piece is at most two neighbours from the start of the chord,
 * so it is within 2 sqrt(2) 3 B 2^200 h <= 2^235.1 h of S.  At the cap, k =
 * 240, that is 2^-4.9 = 0.034 units; the computed Q_i and the computed S differ from
 * the exact ones by E each, so the test sees at most 0.034 + 2 E + 2^-16 <
 * 0.04 units, which is less than delta >= 1: the piece is flat.  (236 levels
 * would do, at 0.54; the cap has margin.)  The same holds for any `t`, since
 * delta only grows.  So the cap is never what ends a piece, and the theorem
 * holds with the cap in place.  (A weight ratio of 2^200 is what makes the
 * proof need 236 levels rather than 20: a uniform split of such a curve spends
 * about a hundred levels before it divides the curve's length at all, and each of
 * them emits one line, which is the price of an input no one writes.)
 *
 * ## What checks the proof
 *
 * tests/unit/test_conic.cpp evaluates conics and arcs densely and exactly (the
 * Bernstein form in `long double` at dyadic parameters, which are exact in the
 * type) at thirty-three parameters inside every piece the route made, and, where
 * the weights are far apart, at 2^-k and 1 - 2^-k for every depth the proof
 * allows, on every size class up to the whole `int32` range and at tolerances
 * from 2 up; it asserts the bound, the piece count and a zero cap count.  A
 * variant of the route with a smaller cap (gshgeom_conic_route_capped) does
 * bind on a curve of extreme weights, the count says so and the bound is lost,
 * which shows that the check sees a violation.
 */

/* The cap (flatten.h) bounds the explicit stack: a pending right half exists at
 * each depth below the piece being worked on. */
#define CONIC_STACK (GSHGEOM_CONIC_MAX_DEPTH + 1)

/* A point in homogeneous coordinates: (w X, w Y, w). */
typedef struct Hp {
  double x, y, w;
} Hp;

static Hp hp_mid(Hp a, Hp b) {
  Hp m;
  m.x = 0.5 * (a.x + b.x);
  m.y = 0.5 * (a.y + b.y);
  m.w = 0.5 * (a.w + b.w);
  return m;
}

/* The squared distance from (px, py) to the segment (sx, sy) - (ex, ey). */
static double dist2_to_segment(
    double px, double py, double sx, double sy, double ex, double ey) {
  double dx = ex - sx;
  double dy = ey - sy;
  double ax = px - sx;
  double ay = py - sy;
  double len2 = dx * dx + dy * dy;
  if (!(len2 > 0.0)) {
    return ax * ax + ay * ay;
  }
  double u = (ax * dx + ay * dy) / len2;
  if (u < 0.0) {
    u = 0.0;
  } else if (u > 1.0) {
    u = 1.0;
  }
  double qx = ax - u * dx;
  double qy = ay - u * dy;
  return qx * qx + qy * qy;
}

/* One piece of the curve, with its chord (the computed, clamped ends) and where it is
 * on the parameter interval: [leaf / 2^depth, (leaf + 1) / 2^depth], with leaf
 * exact only while depth <= 32. */
typedef struct Piece {
  Hp h[4];
  double sx, sy, ex, ey;
  int depth;
  uint32_t leaf;
} Piece;

static double clamp_to(double v, double lo, double hi) {
  return v < lo ? lo : (v > hi ? hi : v);
}

GSHGEOM_Result gshgeom_conic_route_capped(int n, const double * sx, const double * sy,
    const double * w, int32_t tolerance, int levels, GSHGEOM_Flatten_Emit emit,
    void * user, int * forced, int * deepest) {
  if ((n != 3 && n != 4) || !sx || !sy || !w || !emit || tolerance < GSHGEOM_FLATTEN_MIN_TOLERANCE ||
      levels < 0 || levels > GSHGEOM_CONIC_MAX_DEPTH) {
    return GSHGEOM_ERR_INVALID;
  }
  if (forced) {
    *forced = 0;
  }
  if (deepest) {
    *deepest = 0;
  }
  const double delta = (double)(tolerance - 1);
  const double delta_sq = delta * delta;
  double xlo = sx[0], xhi = sx[0], ylo = sy[0], yhi = sy[0];
  for (int i = 1; i < n; i++) {
    xlo = sx[i] < xlo ? sx[i] : xlo;
    xhi = sx[i] > xhi ? sx[i] : xhi;
    ylo = sy[i] < ylo ? sy[i] : ylo;
    yhi = sy[i] > yhi ? sy[i] : yhi;
  }
  Piece stack[CONIC_STACK];
  size_t sp = 0;
  Piece cur;
  for (int i = 0; i < n; i++) {
    cur.h[i].x = w[i] * sx[i];
    cur.h[i].y = w[i] * sy[i];
    cur.h[i].w = w[i];
  }
  cur.sx = sx[0];
  cur.sy = sy[0];
  cur.ex = sx[n - 1];
  cur.ey = sy[n - 1];
  cur.depth = 0;
  cur.leaf = 0;
  for (;;) {
    int flat = 1;
    for (int i = 1; i < n - 1; i++) {
      double px = cur.h[i].x / cur.h[i].w;
      double py = cur.h[i].y / cur.h[i].w;
      if (dist2_to_segment(px, py, cur.sx, cur.sy, cur.ex, cur.ey) > delta_sq) {
        flat = 0;
        break;
      }
    }
    if (flat || cur.depth >= levels) {
      if (!flat && forced) {
        (*forced)++;
      }
      if (deepest && cur.depth > *deepest) {
        *deepest = cur.depth;
      }
      GSHGEOM_Device_Point d;
      GSHGEOM_Result r = round_scaled(cur.ex, &d.x);
      if (r == GSHGEOM_OK) {
        r = round_scaled(cur.ey, &d.y);
      }
      if (r == GSHGEOM_OK) {
        r = emit(user, d, cur.leaf, cur.depth);
      }
      if (r != GSHGEOM_OK) {
        return r;
      }
      if (sp == 0) {
        return GSHGEOM_OK;
      }
      cur = stack[--sp];
      continue;
    }
    Piece left, right;
    Hp m;
    if (n == 3) {
      Hp a = hp_mid(cur.h[0], cur.h[1]);
      Hp b = hp_mid(cur.h[1], cur.h[2]);
      m = hp_mid(a, b);
      left.h[0] = cur.h[0]; left.h[1] = a; left.h[2] = m;
      right.h[0] = m; right.h[1] = b; right.h[2] = cur.h[2];
    } else {
      Hp a = hp_mid(cur.h[0], cur.h[1]);
      Hp b = hp_mid(cur.h[1], cur.h[2]);
      Hp c = hp_mid(cur.h[2], cur.h[3]);
      Hp d = hp_mid(a, b);
      Hp e = hp_mid(b, c);
      m = hp_mid(d, e);
      left.h[0] = cur.h[0]; left.h[1] = a; left.h[2] = d; left.h[3] = m;
      right.h[0] = m; right.h[1] = e; right.h[2] = c; right.h[3] = cur.h[3];
    }
    double mx = clamp_to(m.x / m.w, xlo, xhi);
    double my = clamp_to(m.y / m.w, ylo, yhi);
    left.sx = cur.sx; left.sy = cur.sy; left.ex = mx; left.ey = my;
    right.sx = mx; right.sy = my; right.ex = cur.ex; right.ey = cur.ey;
    left.depth = right.depth = cur.depth + 1;
    left.leaf = cur.leaf * 2u;
    right.leaf = cur.leaf * 2u + 1u;
    stack[sp++] = right;
    cur = left;
  }
}

GSHGEOM_Result gshgeom_conic_route(int n, const double * sx, const double * sy,
    const double * w, int32_t tolerance, GSHGEOM_Flatten_Emit emit, void * user) {
  int forced = 0;
  GSHGEOM_Result r = gshgeom_conic_route_capped(n, sx, sy, w, tolerance,
      GSHGEOM_CONIC_MAX_DEPTH, emit, user, &forced, NULL);
  /* The proof is that the cap never ends a piece that is not flat. */
  return r == GSHGEOM_OK && forced ? GSHGEOM_ERR_INTERNAL : r;
}

static GSHGEOM_Result emit_line(
    void * user, GSHGEOM_Device_Point to, uint32_t leaf, int level) {
  (void)leaf;
  (void)level;
  return gshgeom_device_path_line_to((GSHGEOM_Device_Path *)user, to);
}


/* Map a model point to device units, scaled but not rounded (x * 64 is exact),
 * and check that it rounds into range. */
static GSHGEOM_Result map_scaled(const GSHGEOM_Map * map, GSHGEOM_Point p,
    double * x, double * y) {
  GSHGEOM_Point m;
  GSHGEOM_Result r = gshgeom_map_apply_point(map, p, &m);
  if (r != GSHGEOM_OK) {
    return r;
  }
  double sx = m.x * SCALE;
  double sy = m.y * SCALE;
  int32_t probe;
  r = round_scaled(sx, &probe);
  if (r == GSHGEOM_OK) {
    r = round_scaled(sy, &probe);
  }
  if (r == GSHGEOM_OK) {
    *x = sx;
    *y = sy;
  }
  return r;
}

static GSHGEOM_Result map_rounded(const GSHGEOM_Map * map, GSHGEOM_Point p,
    GSHGEOM_Device_Point * out) {
  double x, y;
  GSHGEOM_Result r = map_scaled(map, p, &x, &y);
  if (r == GSHGEOM_OK) {
    r = round_scaled(x, &out->x);
  }
  if (r == GSHGEOM_OK) {
    r = round_scaled(y, &out->y);
  }
  return r;
}

/* The conic route (above) for one rational curve of @p n control points. */
static GSHGEOM_Result convert_conic(GSHGEOM_Device_Path * dp, const GSHGEOM_Map * map,
    int n, const GSHGEOM_Point * p, const double * w, int32_t tolerance) {
  double sx[4], sy[4];
  for (int i = 0; i < n; i++) {
    GSHGEOM_Result r = map_scaled(map, p[i], &sx[i], &sy[i]);
    if (r != GSHGEOM_OK) {
      return r;
    }
  }
  return gshgeom_conic_route(n, sx, sy, w, tolerance, emit_line, dp);
}

/* An arc is its pieces, each a rational quadratic: the same route, so that the pieces
 * meet at the same device point (a shared control point is mapped to the same bits), the
 * ends are the arc's own ends rounded once and a closed arc closes. */
static GSHGEOM_Result convert_arc(GSHGEOM_Device_Path * dp, const GSHGEOM_Map * map,
    const GSHGEOM_Arc * arc, int32_t tolerance) {
  GSHGEOM_Curve pieces[GSHGEOM_ARC_MAX_PIECES];
  int count = 0;
  GSHGEOM_Result r = gshgeom_arc_pieces(arc, pieces, &count);
  for (int k = 0; k < count && r == GSHGEOM_OK; k++) {
    r = convert_conic(dp, map, 3, pieces[k].p, pieces[k].w, tolerance);
  }
  return r;
}

static GSHGEOM_Result convert_element(GSHGEOM_Device_Path * dp,
    const GSHGEOM_Map * map, const GSHGEOM_Element * el, int32_t tolerance) {
  const GSHGEOM_Curve * c = &el->curve;
  GSHGEOM_Device_Point q[3];
  GSHGEOM_Result r = GSHGEOM_OK;
  switch (el->verb) {
  case GSHGEOM_VERB_MOVE:
    r = map_rounded(map, c->p[1], &q[0]);
    return r != GSHGEOM_OK ? r : gshgeom_device_path_move_to(dp, q[0]);
  case GSHGEOM_VERB_LINE:
    r = map_rounded(map, c->p[1], &q[0]);
    return r != GSHGEOM_OK ? r : gshgeom_device_path_line_to(dp, q[0]);
  case GSHGEOM_VERB_QUADRATIC:
    r = map_rounded(map, c->p[1], &q[0]);
    if (r == GSHGEOM_OK) {
      r = map_rounded(map, c->p[2], &q[1]);
    }
    return r != GSHGEOM_OK ? r : gshgeom_device_path_quadratic_to(dp, q[0], q[1]);
  case GSHGEOM_VERB_CUBIC:
    r = map_rounded(map, c->p[1], &q[0]);
    if (r == GSHGEOM_OK) {
      r = map_rounded(map, c->p[2], &q[1]);
    }
    if (r == GSHGEOM_OK) {
      r = map_rounded(map, c->p[3], &q[2]);
    }
    return r != GSHGEOM_OK
        ? r
        : gshgeom_device_path_cubic_to(dp, q[0], q[1], q[2]);
  case GSHGEOM_VERB_RATIONAL_QUADRATIC:
    return convert_conic(dp, map, 3, c->p, c->w, tolerance);
  case GSHGEOM_VERB_RATIONAL_CUBIC:
    return convert_conic(dp, map, 4, c->p, c->w, tolerance);
  case GSHGEOM_VERB_CLOSE:
    return gshgeom_device_path_close(dp);
  case GSHGEOM_VERB_ARC: /* never reached: an arc is not an element (gshgeom_path_element), see the loop below */
  case GSHGEOM_VERB_COUNT:
    break;
  }
  return GSHGEOM_ERR_INTERNAL;
}

GSHGEOM_Result gshgeom_device_path_from_model(const GSHGEOM_Path * path,
    const GSHGEOM_Map * map, const GSHGEOM_Allocator * alloc,
    const GSHGEOM_Limits * limits, GSHGEOM_Device_Path ** out) {
  return gshgeom_device_path_from_model_with_tolerance(
      path, map, GSHGEOM_FLATTEN_TOLERANCE, alloc, limits, out);
}

GSHGEOM_Result gshgeom_device_path_from_model_with_tolerance(const GSHGEOM_Path * path,
    const GSHGEOM_Map * map, int32_t tolerance, const GSHGEOM_Allocator * alloc,
    const GSHGEOM_Limits * limits, GSHGEOM_Device_Path ** out) {
  if (!path || !map || !out || gshgeom_path_dimension(path) != 2 ||
      tolerance < GSHGEOM_FLATTEN_MIN_TOLERANCE) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Result r = gshgeom_map_check(map);
  if (r != GSHGEOM_OK) {
    return r;
  }
  if (map->dim != 2) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Device_Path * dp;
  r = gshgeom_device_path_create(alloc, limits, &dp);
  if (r != GSHGEOM_OK) {
    return r;
  }
  size_t n = gshgeom_path_element_count(path);
  for (size_t i = 0; i < n && r == GSHGEOM_OK; i++) {
    GSHGEOM_Verb verb;
    r = gshgeom_path_verb(path, i, &verb);
    if (r != GSHGEOM_OK) {
      r = GSHGEOM_ERR_INTERNAL;
      break;
    }
    if (verb == GSHGEOM_VERB_ARC) {
      GSHGEOM_Arc arc;
      r = gshgeom_path_arc(path, i, &arc);
      r = r != GSHGEOM_OK ? GSHGEOM_ERR_INTERNAL : convert_arc(dp, map, &arc, tolerance);
      continue;
    }
    GSHGEOM_Element el;
    r = gshgeom_path_element(path, i, &el);
    if (r != GSHGEOM_OK) {
      r = GSHGEOM_ERR_INTERNAL;
      break;
    }
    r = convert_element(dp, map, &el, tolerance);
  }
  if (r != GSHGEOM_OK) {
    gshgeom_device_path_destroy(dp);
    return r;
  }
  *out = dp;
  return GSHGEOM_OK;
}

/* ---------------------------------------------------------------------- */
/* Hash and dump                                                           */
/* ---------------------------------------------------------------------- */

/** The format version of the hash and the dump. */
#define DEVICE_FORMAT_VERSION 1u
/** The kind field: "DEVP" as little-endian bytes. */
#define DEVICE_KIND 0x50564544u

uint64_t gshgeom_device_path_hash(const GSHGEOM_Device_Path * path) {
  GSHGEOM_Hash h;
  gshgeom_hash_init(&h);
  if (!path) {
    return gshgeom_hash_value(&h);
  }
  gshgeom_hash_u32(&h, DEVICE_FORMAT_VERSION);
  gshgeom_hash_u32(&h, DEVICE_KIND);
  gshgeom_hash_u64(&h, (uint64_t)path->el_n);
  for (size_t i = 0; i < path->el_n; i++) {
    const Elem * e = &path->el[i];
    gshgeom_hash_u8(&h, e->verb);
    if (e->verb == GSHGEOM_DEVICE_CLOSE) {
      continue;
    }
    for (unsigned k = 0; k < verb_points(e->verb); k++) {
      gshgeom_hash_i32(&h, path->pt[e->point + k].x);
      gshgeom_hash_i32(&h, path->pt[e->point + k].y);
    }
  }
  return gshgeom_hash_value(&h);
}

static int put_i32(FILE * out, int32_t v) {
  return fprintf(out, " %08x", (unsigned)(uint32_t)v);
}

GSHGEOM_Result gshgeom_device_path_dump(
    const GSHGEOM_Device_Path * path, FILE * out) {
  static const char names[GSHGEOM_DEVICE_VERB_COUNT] = {'M', 'L', 'Q', 'C', 'Z'};
  if (!path || !out) {
    return GSHGEOM_ERR_INVALID;
  }
  /* The count goes through an unsigned long long: the width of size_t is a
   * property of the host. */
  if (fprintf(out, "gshgeom.device-path %u elements=%llu hash=%016llx\n",
          DEVICE_FORMAT_VERSION, (unsigned long long)path->el_n,
          (unsigned long long)gshgeom_device_path_hash(path)) < 0) {
    return GSHGEOM_ERR_IO;
  }
  for (size_t i = 0; i < path->el_n; i++) {
    const Elem * e = &path->el[i];
    if (fputc(names[e->verb], out) == EOF) {
      return GSHGEOM_ERR_IO;
    }
    if (e->verb != GSHGEOM_DEVICE_CLOSE) {
      for (unsigned k = 0; k < verb_points(e->verb); k++) {
        if (put_i32(out, path->pt[e->point + k].x) < 0 ||
            put_i32(out, path->pt[e->point + k].y) < 0) {
          return GSHGEOM_ERR_IO;
        }
      }
    }
    if (fputc('\n', out) == EOF) {
      return GSHGEOM_ERR_IO;
    }
  }
  /* A buffered stream reports a write error at the flush, and a stream's error flag is sticky:
   * a failed write whose return value was lost (an unbuffered stream's fprintf) still sets it. */
  return (fflush(out) == 0 && !ferror(out)) ? GSHGEOM_OK : GSHGEOM_ERR_IO;
}
