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
 * How a piece of a stroke becomes positively wound 26.6 polygons, and the proof
 * that the outline keeps SG17 (construction revision 1).
 *
 * ## The construction
 *
 * The region `R` of a part (stroke.h, SG6 to SG12) is the union of finitely
 * many **pieces** `P_j`: the body of a stretch of one segment, a join, a cap, a
 * degenerate shape.  For each the expansion forms a **core polygon** `K0_j`, the
 * convex hull of a few points it computes, with two bounds in device pixels that
 * the file building the piece proves:
 *
 *     P_j  is within d_in of K0_j   (P_j <= K0_j + D(d_in)),
 *     K0_j is within d_out of P_j   (K0_j <= P_j + D(d_out)).
 *
 * It then emits `K_j`, the convex hull of the core points each replaced by the
 * eight vertices of a regular octagon of inradius `r = d_in + rho`, where
 * `rho = (ROUND + ETA) / 64` pixel; rounds every vertex of `K_j` to 26.6
 * (::gshgeom_round_26_6); and writes the rounded hull as one closed subpath when
 * every triangle of its fan from the first vertex has a non-negative exact
 * orientation, and otherwise each triangle of the fan as a closed subpath of its
 * own, turned to positive orientation, a triangle of zero area being dropped.
 * The outline of a part is every subpath so written, in SG18's order.
 *
 * ## The theorem
 *
 * For every input the expansion accepts, with `tau` the tolerance in 26.6 units
 * (at least 2): every point of `R` has a positive winding number in the
 * outline; every point with a nonzero winding number is within
 * `(tau - 1/4) / 64` pixel of `R`; and no subpath has a negative winding number
 * anywhere.  So `R <= fill(O) <= R + D(tau - 1/4)` (stronger than SG17), and the
 * outline's nonzero fill is a union (SG13: overlaps never cancel).
 *
 * ## The proof
 *
 * **(a) The arithmetic.**  Every device coordinate the expansion computes is a
 * fixed composition of at most a few hundred `+ - * /` and `gmath_sqrt` (each
 * correctly rounded, unit roundoff `u = 2^-53`), the de Casteljau subdivisions
 * included (at most ::GSHGEOM_STROKE_MAX_DEPTH levels of degree at most 3), on
 * quantities whose magnitudes the expansion checks against
 * ::GSHGEOM_SK_MAG_LIMIT (2^30 pixels) before it builds anything: anchors, the
 * mapped control points' offsets, `h` in device pixels, the miter limit times
 * that.  A bound of `2^10 u` relative to `2^31` pixels is `2^-12` pixel, which is
 * ETA (1/64 unit).  The quantities a piece's bounds depend on and whose
 * relative error is not small (a direction from a short hodograph coefficient, a
 * miter tip from nearly parallel lines) carry their own error bound into that
 * piece's `d_in` and `d_out` instead, and a piece whose bounds do not fit is
 * split or refused, never emitted.  So each computed core point is within
 * ETA of the point the piece's proof is about.
 *
 * **(b) The inflation contains the piece with room for the rounding.**  An
 * octagon of inradius `r` contains the disc `D(r)`, so `K_j` contains
 * `K0_j + D(r)` (the hull of the octagons around the core points is the
 * Minkowski sum of the hull and the octagon).  With (a), the computed core points
 * are within ETA of the proof's, so `K_j` contains `K0_j + D(d_in + ROUND)` and
 * hence `P_j + D(ROUND)` (the first bound above).
 *
 * **(c) The inflation stays near the piece.**  The octagon's axis vertices are at
 * `1.0824 r` (1.0824 >= 1 / cos(pi / 8), so its inradius is at least `r`) and
 * its diagonal ones, written with coordinates `0.7072 >= sqrt(2) / 2`, at
 * `1.0824 * 0.7072 * sqrt(2) r < OCT r` (OCT = 1.0826), so `K_j` is
 * within `OCT r` of `K0_j`, within `d_out + OCT r` of `P_j` (the second bound),
 * and after the rounding (each vertex moves at most `ROUND` = sqrt(2) / 2 unit,
 * and a point of the hull of moved points is within that of the hull) within
 * `d_out + OCT (d_in + rho) + rho` of `P_j`.  ::gshgeom_sk_fits accepts a piece
 * only when that is at most `budget = (tau - MARGIN) / 64` pixel.
 *
 * **(d) Coverage.**  `K_j` is convex and its exact hull (exact orientation
 * predicates on the computed vertices) winds once, positively, around every
 * point inside it.  Rounding moves every vertex at most ROUND, and a point
 * of a polygon edge at most that (a convex combination of its ends), so for every
 * point `x` farther than ROUND from `K_j`'s boundary the rounded hull winds as
 * the exact one does (the straight-line homotopy between the two boundaries
 * avoids `x`).  By (b) every point of `P_j` is such a point inside `K_j`, so the
 * rounded hull winds once around it.  The fan of the rounded hull from its first
 * vertex has signed windings that sum to the hull's (the inner edges cancel); so
 * some fan triangle winds positively around `x`.  Whether the hull is written
 * whole (every fan triangle non-negative) or triangle by triangle (each turned
 * positive, which keeps the positive ones as they are), `x` has a positive
 * winding in what is written.
 *
 * **(e) Nothing outside.**  Every fan triangle lies in the convex hull of the
 * rounded vertices, which by (c) is within the budget of `P_j`, hence of `R`.  A
 * subpath written whole has as its winding the sum of its fan triangles', each 0
 * or 1, and a triangle written alone is positive: no subpath winds negatively,
 * and a nonzero winding number of the outline means some subpath winds
 * positively, so the point lies in some piece's fan, within the budget of `R`.
 *
 * **(f) The margin.**  A piece is emitted only within `budget`, so the outline is
 * within `tau - 1/4` units of `R` and contains `R`: the proof closes with a
 * margin of 1/4 unit, the rounding counted once on each side of the
 * inflation.  For the least tolerance, 2 units, a piece's own bounds have
 * `budget - (1 + OCT) rho = 1.75 - 1.5052 = 0.2448` units between them.
 *
 * The bounds of each kind of piece are proved where it is built: body.c (the
 * band of a stretch of curve, the sector of a short stretch), join.c (joins,
 * caps, degenerate shapes), arc.c (circular arcs and the arcs join).  What
 * checks the theorem is tests/unit/test_stroke_oracle.cpp: an oracle that
 * decides membership from the definition alone.
 */

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/math/basic.h>
#include <ghoti.io/shape-geometry/orient.h>

#include "../core/grow.h"
#include "../flatten/flatten_internal.h"
#include "stroke_internal.h"

#include <string.h>

/** The octagon's vertices for unit circumradius scale: the diagonals written at 0.7072 >= sqrt(2) / 2. */
static const double OCT_DIR[8][2] = {{1.0, 0.0}, {0.7072, 0.7072}, {0.0, 1.0}, {-0.7072, 0.7072}, {-1.0, 0.0},
    {-0.7072, -0.7072}, {0.0, -1.0}, {0.7072, -0.7072}};

double gshgeom_sk_norm(GSHGEOM_SK_Vec a) {
  /* scaled by the larger coordinate, so that neither the squares of a huge vector overflow nor those of a tiny one
   * underflow to zero (a direction of a segment of 1e-200 is as good as any other) */
  double ax = a.x < 0.0 ? -a.x : a.x, ay = a.y < 0.0 ? -a.y : a.y;
  double m = ax > ay ? ax : ay;
  if (m == 0.0 || !(m <= 1.7976931348623157e308)) {
    return m; /* zero, an infinity, or a NaN stays itself */
  }
  double x = ax / m, y = ay / m;
  return m * gmath_sqrt(x * x + y * y);
}

GSHGEOM_SK_Vec gshgeom_sk_unit(GSHGEOM_SK_Vec a) {
  double n = gshgeom_sk_norm(a);
  return gshgeom_sk_v(a.x / n, a.y / n);
}

GSHGEOM_SK_Vec gshgeom_sk_lin(const GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec d) {
  double x = c->lin[0][0] * d.x;
  x = x + c->lin[0][1] * d.y;
  double y = c->lin[1][0] * d.x;
  y = y + c->lin[1][1] * d.y;
  return gshgeom_sk_v(x, y);
}

GSHGEOM_Result gshgeom_sk_anchor(const GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec o, GSHGEOM_SK_Vec * out) {
  if (!c->local) {
    *out = o;
    return GSHGEOM_OK;
  }
  GSHGEOM_Point p = {o.x, o.y, 0.0};
  GSHGEOM_Point q;
  GSHGEOM_Result r = gshgeom_map_apply_point(&c->map, p, &q);
  if (r != GSHGEOM_OK) {
    return r;
  }
  *out = gshgeom_sk_v(q.x, q.y);
  return GSHGEOM_OK;
}

int gshgeom_sk_fits(const GSHGEOM_SK_Ctx * c, double d_in, double d_out) {
  if (!(d_in >= 0.0) || !(d_out >= 0.0)) {
    return 0; /* NaN or a negative bound: a computation that went wrong is never accepted */
  }
  if (c->cv) {
    /* curve mode (SG35): nothing is inflated or rounded, the piece is written as it is, so each bound is its own */
    return d_in <= c->cv->budget && d_out <= c->cv->budget;
  }
  return d_out + GSHGEOM_SK_OCT * (d_in + c->rho) + c->rho <= c->budget;
}

double gshgeom_sk_room(const GSHGEOM_SK_Ctx * c, double d_in) {
  if (c->cv) {
    return d_in <= c->cv->budget ? c->cv->budget : -1.0;
  }
  return c->budget - GSHGEOM_SK_OCT * (d_in + c->rho) - c->rho;
}

/** Lexicographic order for the hull. */
static int before(GSHGEOM_SK_Vec a, GSHGEOM_SK_Vec b) {
  return a.x < b.x || (a.x == b.x && a.y < b.y);
}

/** Sort by (x, y): a merge sort, so the order (and the hull) does not depend on the C library. */
static void sort_points(GSHGEOM_SK_Vec * a, GSHGEOM_SK_Vec * tmp, size_t n) {
  for (size_t width = 1; width < n; width *= 2) {
    for (size_t lo = 0; lo < n; lo += 2 * width) {
      size_t mid = lo + width < n ? lo + width : n;
      size_t hi = lo + 2 * width < n ? lo + 2 * width : n;
      size_t i = lo, j = mid, k = lo;
      while (i < mid && j < hi) {
        tmp[k++] = before(a[j], a[i]) ? a[j++] : a[i++];
      }
      while (i < mid) {
        tmp[k++] = a[i++];
      }
      while (j < hi) {
        tmp[k++] = a[j++];
      }
    }
    memcpy(a, tmp, n * sizeof *a);
  }
}

static int orient(GSHGEOM_SK_Vec a, GSHGEOM_SK_Vec b, GSHGEOM_SK_Vec c) {
  GSHGEOM_Point p = {a.x, a.y, 0.0};
  GSHGEOM_Point q = {b.x, b.y, 0.0};
  GSHGEOM_Point r = {c.x, c.y, 0.0};
  int s = 0;
  (void)gshgeom_orient2d(p, q, r, &s); /* finite by (a); the predicate only refuses a non-finite input */
  return s;
}

/**
 * The convex hull of pts[0..n) by the monotone chain with exact orientations,
 * in positive orientation, collinear points dropped.  Writes it to hull (room
 * for n + 1) and returns its size.
 */
static size_t convex_hull(GSHGEOM_SK_Vec * pts, GSHGEOM_SK_Vec * tmp, size_t n, GSHGEOM_SK_Vec * hull) {
  sort_points(pts, tmp, n);
  size_t k = 0;
  for (size_t i = 0; i < n; i++) {
    while (k >= 2 && orient(hull[k - 2], hull[k - 1], pts[i]) <= 0) {
      k--;
    }
    hull[k++] = pts[i];
  }
  size_t lower = k + 1;
  for (size_t i = n - 1; i-- > 0;) {
    while (k >= lower && orient(hull[k - 2], hull[k - 1], pts[i]) <= 0) {
      k--;
    }
    hull[k++] = pts[i];
  }
  return k > 1 ? k - 1 : k;
}

static int orient_dev(GSHGEOM_Device_Point a, GSHGEOM_Device_Point b, GSHGEOM_Device_Point c) {
  return orient(gshgeom_sk_v((double)a.x, (double)a.y), gshgeom_sk_v((double)b.x, (double)b.y),
      gshgeom_sk_v((double)c.x, (double)c.y));
}

static GSHGEOM_Result write_ring(GSHGEOM_Flat_Path * fp, const GSHGEOM_Device_Point * v, size_t n) {
  GSHGEOM_Result r = gshgeom_flat_path_begin(fp, v[0]);
  for (size_t i = 1; i < n && r == GSHGEOM_OK; i++) {
    r = gshgeom_flat_path_add(fp, v[i]);
  }
  if (r == GSHGEOM_OK) {
    gshgeom_flat_path_close_last(fp);
  }
  return r;
}

/** Make room for @p need vectors in the context's scratch. */
static GSHGEOM_Result room_for(GSHGEOM_SK_Ctx * c, size_t need) {
  void * grown = c->scratch;
  GSHGEOM_Result r = gshgeom_grow(&c->alloc, &grown, &c->scratch_cap, need, sizeof *c->scratch);
  if (r == GSHGEOM_OK) {
    c->scratch = grown;
  }
  return r;
}

GSHGEOM_Result gshgeom_sk_emit(GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec anchor, const GSHGEOM_SK_Vec * d, size_t n,
    double d_in, double d_out) {
  if (n == 0) {
    return GSHGEOM_OK;
  }
  if (!gshgeom_sk_fits(c, d_in, d_out)) {
    return GSHGEOM_ERR_INTERNAL;
  }
  if (n > SIZE_MAX / 32) {
    return GSHGEOM_ERR_OOM;
  }
  /* scratch: the core points, a sort buffer and the core hull (3n + 1), then the candidates, their sort buffer and
   * their hull (each at most 8 per core hull vertex: 24n + 3) */
  GSHGEOM_Result r = room_for(c, 27 * n + 8);
  if (r != GSHGEOM_OK) {
    return r;
  }
  GSHGEOM_SK_Vec * core = c->scratch;
  GSHGEOM_SK_Vec * tmp = core + n;
  GSHGEOM_SK_Vec * hull1 = tmp + n;
  for (size_t i = 0; i < n; i++) {
    core[i] = gshgeom_sk_add(anchor, gshgeom_sk_lin(c, d[i]));
  }
  size_t h1 = convex_hull(core, tmp, n, hull1);
  if (c->cv) {
    /* curve mode: the core polygon is the piece, as lines, in target units; the hull is positive in target units */
    return gshgeom_skc_ring(c, hull1, h1);
  }
  /* An octagon of inradius at least d_in + rho: its axis vertices at 1.0824 >= 1 / cos(pi / 8) times that, its
   * diagonal ones (written 0.7072) at 1.0824 * 0.7072 * sqrt(2) < 1.0826 = OCT times it. */
  double rc = (d_in + c->rho) * 1.0824;
  GSHGEOM_SK_Vec * cand = hull1 + h1 + 1;
  size_t m = 0;
  for (size_t i = 0; i < h1; i++) {
    GSHGEOM_SK_Vec nin = gshgeom_sk_v(0.0, 0.0), nout = nin;
    if (h1 >= 3) {
      GSHGEOM_SK_Vec ein = gshgeom_sk_sub(hull1[i], hull1[(i + h1 - 1) % h1]);
      GSHGEOM_SK_Vec eout = gshgeom_sk_sub(hull1[(i + 1) % h1], hull1[i]);
      nin = gshgeom_sk_left(ein);
      nout = gshgeom_sk_left(eout);
    }
    for (int k = 0; k < 8; k++) {
      GSHGEOM_SK_Vec o = gshgeom_sk_v(OCT_DIR[k][0], OCT_DIR[k][1]);
      /* A vertex of the sum is a hull vertex plus an octagon vertex whose direction lies between the outward
       * normals of the hull's two edges there (an angle below pi): such a direction is within a quarter turn of
       * one of them.  Every other candidate is interior, and keeping a few is harmless. */
      if (h1 < 3 || gshgeom_sk_dot(o, nin) >= 0.0 || gshgeom_sk_dot(o, nout) >= 0.0) {
        cand[m++] = gshgeom_sk_v(hull1[i].x + rc * o.x, hull1[i].y + rc * o.y);
      }
    }
  }
  GSHGEOM_SK_Vec * ctmp = cand + m;
  GSHGEOM_SK_Vec * hull = ctmp + m;
  size_t h = convex_hull(cand, ctmp, m, hull);
  if (h < 3) {
    return GSHGEOM_ERR_INTERNAL; /* an octagon has area: so does the hull of a translate of it */
  }
  GSHGEOM_Device_Point stack[64];
  GSHGEOM_Device_Point * v = stack;
  if (h > 64) {
    v = gcu_allocator_malloc(&c->alloc, h * sizeof *v);
    if (!v) {
      return GSHGEOM_ERR_OOM;
    }
  }
  for (size_t i = 0; i < h && r == GSHGEOM_OK; i++) {
    r = gshgeom_round_26_6(hull[i].x, &v[i].x);
    if (r == GSHGEOM_OK) {
      r = gshgeom_round_26_6(hull[i].y, &v[i].y);
    }
  }
  if (r == GSHGEOM_OK) {
    int whole = 1;
    for (size_t i = 1; i + 1 < h; i++) {
      if (orient_dev(v[0], v[i], v[i + 1]) < 0) {
        whole = 0;
        break;
      }
    }
    if (whole) {
      r = write_ring(c->part, v, h);
    } else {
      for (size_t i = 1; i + 1 < h && r == GSHGEOM_OK; i++) {
        int s = orient_dev(v[0], v[i], v[i + 1]);
        if (s == 0) {
          continue;
        }
        GSHGEOM_Device_Point t[3] = {v[0], s > 0 ? v[i] : v[i + 1], s > 0 ? v[i + 1] : v[i]};
        r = write_ring(c->part, t, 3);
      }
    }
  }
  if (v != stack) {
    gcu_allocator_free(&c->alloc, v);
  }
  return r;
}
