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
  if (c->conform && !c->cv) {
    c->inflated = 1; /* the pieces beside this one are conforming and do not overlap it: the pass is void */
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


/* --- conforming pieces: shared lattice points, no inflation ---------------------------------------------------------------- */

/*
 * A conforming piece is the hull of its core points with no inflation, each vertex put on the 26.6 lattice through a table
 * that gives two computed points that are the same point one lattice point.  The pieces that meet along a seam (the bodies
 * and bands of a segment, the joins and caps beside them) compute the seam's two ends the same way, so they write the seam
 * as two edges that are exact reverses, which ::gshgeom_flat_path_cancel_seams removes.  Nothing overlaps, so the error is
 * two-sided, `err + rho <= budget`: the polygon is within `err` of the piece and the lattice within `rho` of the polygon.
 *
 * A piece that only approximates its region (a fan, a disc, an arcs join) has no edge in common with its neighbours, and
 * relies on its inflation to cover the seam.  The covering is of the exact seam, and a conforming neighbour's edge is the
 * exact one moved by up to a lattice rounding, so a hairline gap opens between them.  A stroke that emits one such piece
 * while conforming (`inflated`) is therefore expanded again with every piece inflated, as it was before this existed.
 */

int gshgeom_sk_fits_conform(const GSHGEOM_SK_Ctx * c, double d_in, double err) {
  if (!(err >= 0.0)) {
    return 0; /* NaN or a negative bound: a computation that went wrong is never accepted */
  }
  if (c->cv) {
    return err <= c->cv->budget;
  }
  if (!c->conform) {
    return gshgeom_sk_fits(c, d_in, err);
  }
  return err + c->rho <= c->budget;
}

double gshgeom_sk_room_conform(const GSHGEOM_SK_Ctx * c, double d_in) {
  if (c->cv) {
    return c->cv->budget;
  }
  if (!c->conform) {
    return gshgeom_sk_room(c, d_in);
  }
  return c->budget - c->rho;
}

/** Bins of 2^-18 pixel, and the distance within which two computed points are the same point: far above the rounding of the
 * arithmetic that computes them (2^-40 pixel at a 2^30 pixel extent and below), far below the lattice (2^-6). */
#define SNAP_BIN 0x1p18
#define SNAP_TOL 0x1p-19

static uint64_t snap_hash(int64_t qx, int64_t qy) {
  uint64_t h = ((uint64_t)qx * UINT64_C(0x9E3779B97F4A7C15)) ^ ((uint64_t)qy * UINT64_C(0xC2B2AE3D27D4EB4F));
  return h ^ (h >> 29);
}

static GSHGEOM_Result snap_grow(GSHGEOM_SK_Ctx * c) {
  size_t cap = c->snap_cap ? c->snap_cap * 2 : 256;
  GSHGEOM_SK_Snap * t = gcu_allocator_calloc(&c->alloc, cap, sizeof *t);
  if (!t) {
    return GSHGEOM_ERR_OOM;
  }
  for (size_t i = 0; i < c->snap_cap; i++) {
    if (c->snap[i].used) {
      size_t k = snap_hash(c->snap[i].qx, c->snap[i].qy) & (cap - 1);
      while (t[k].used) {
        k = (k + 1) & (cap - 1);
      }
      t[k] = c->snap[i];
    }
  }
  gcu_allocator_free(&c->alloc, c->snap);
  c->snap = t;
  c->snap_cap = cap;
  return GSHGEOM_OK;
}

/*
 * The bin of a coordinate in range (at most 2^31 sixty-fourths, so a bin of at most 2^49), as the floor of its scaled value.  It
 * is taken in two int32 steps, since the conversion of a double to an int64_t is an x87 instruction on i686, which the build
 * keeps out of this library.
 */
static int64_t snap_bin(double x) {
  double t = x * SNAP_BIN;
  double hi = gmath_floor(t / 1048576.0);
  return (int64_t)(int32_t)hi * 1048576 + (int32_t)(t - hi * 1048576.0);
}

/**
 * The lattice point of the pixel point (@p x, @p y): the one an earlier point within SNAP_TOL of it was given, else the point
 * rounded to nearest, which is remembered.  The pieces of an outline that meet along an edge compute its ends by different
 * expressions that agree to a few ulps; a rounding that is right for both (an end within an ulp of a half unit would round
 * two ways) is what makes the edge theirs together.
 */
static GSHGEOM_Result snap_round(GSHGEOM_SK_Ctx * c, double x, double y, GSHGEOM_Device_Point * out) {
  GSHGEOM_Device_Point p;
  GSHGEOM_Result r = gshgeom_round_26_6(x, &p.x); /* first: it refuses what is not finite or is out of range */
  if (r == GSHGEOM_OK) {
    r = gshgeom_round_26_6(y, &p.y);
  }
  if (r != GSHGEOM_OK) {
    return r;
  }
  int64_t qx = snap_bin(x), qy = snap_bin(y);
  if (c->snap_cap) {
    for (int dx = -1; dx <= 1; dx++) {
      for (int dy = -1; dy <= 1; dy++) {
        size_t k = snap_hash(qx + dx, qy + dy) & (c->snap_cap - 1);
        while (c->snap[k].used) {
          const GSHGEOM_SK_Snap * e = &c->snap[k];
          if (e->qx == qx + dx && e->qy == qy + dy) {
            double ex = e->x - x, ey = e->y - y;
            if (ex < SNAP_TOL && ex > -SNAP_TOL && ey < SNAP_TOL && ey > -SNAP_TOL) {
              out->x = e->rx;
              out->y = e->ry;
              return GSHGEOM_OK;
            }
          }
          k = (k + 1) & (c->snap_cap - 1);
        }
      }
    }
  }
  *out = p;
  if ((c->snap_n + 1) * 2 > c->snap_cap) {
    r = snap_grow(c);
    if (r != GSHGEOM_OK) {
      return r;
    }
  }
  size_t k = snap_hash(qx, qy) & (c->snap_cap - 1);
  while (c->snap[k].used) {
    k = (k + 1) & (c->snap_cap - 1);
  }
  GSHGEOM_SK_Snap * e = &c->snap[k];
  e->qx = qx;
  e->qy = qy;
  e->x = x;
  e->y = y;
  e->rx = out->x;
  e->ry = out->y;
  e->used = 1;
  c->snap_n++;
  return GSHGEOM_OK;
}

/** Insert @p p into the ring @p v of @p *m vertices on the edge it lies on (within 2 units, between its ends), unless it is one. */
static void ring_insert(GSHGEOM_Device_Point * v, size_t * m, GSHGEOM_Device_Point p) {
  for (size_t i = 0; i < *m; i++) {
    if (v[i].x == p.x && v[i].y == p.y) {
      return;
    }
  }
  for (size_t i = 0; i < *m; i++) {
    GSHGEOM_Device_Point a = v[i], b = v[(i + 1) % *m];
    double ex = (double)b.x - a.x, ey = (double)b.y - a.y, px = (double)p.x - a.x, py = (double)p.y - a.y;
    double len2 = ex * ex + ey * ey, along = px * ex + py * ey, across = ex * py - ey * px;
    /* an edge of under 8 units is too short to take a vertex off its line by a unit and keep the ring simple */
    if (len2 >= 64.0 && along > 0.0 && along < len2 && across * across <= 4.0 * len2) {
      for (size_t j = *m; j > i + 1; j--) {
        v[j] = v[j - 1];
      }
      v[i + 1] = p;
      (*m)++;
      return;
    }
  }
}

GSHGEOM_Result gshgeom_sk_emit_conform(GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec anchor, const GSHGEOM_SK_Vec * d, size_t n,
    double err, double d_in, const GSHGEOM_SK_Vec * piv, size_t npiv) {
  if (c->cv) {
    return gshgeom_sk_emit(c, anchor, d, n, err, err);
  }
  if (!c->conform) {
    return gshgeom_sk_emit(c, anchor, d, n, d_in, err);
  }
  if (n == 0) {
    return GSHGEOM_OK;
  }
  if (!gshgeom_sk_fits_conform(c, d_in, err)) {
    return GSHGEOM_ERR_INTERNAL;
  }
  if (n > SIZE_MAX / 32) {
    return GSHGEOM_ERR_OOM;
  }
  GSHGEOM_Result r = room_for(c, 3 * n + 1);
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
  GSHGEOM_Device_Point stack[72];
  GSHGEOM_Device_Point * v = stack;
  size_t room = h1 + npiv + 2;
  if (room > 72) {
    v = gcu_allocator_malloc(&c->alloc, room * sizeof *v);
    if (!v) {
      return GSHGEOM_ERR_OOM;
    }
  }
  size_t m = 0;
  for (size_t i = 0; i < h1 && r == GSHGEOM_OK; i++) {
    GSHGEOM_Device_Point p;
    r = snap_round(c, hull1[i].x, hull1[i].y, &p);
    if (r == GSHGEOM_OK && (m == 0 || p.x != v[m - 1].x || p.y != v[m - 1].y)) {
      v[m++] = p;
    }
  }
  while (r == GSHGEOM_OK && m > 1 && v[0].x == v[m - 1].x && v[0].y == v[m - 1].y) {
    m--;
  }
  /* fewer than three distinct lattice points: a piece of no area at this resolution, which draws nothing */
  if (r == GSHGEOM_OK && m >= 3) {
    int whole = 1;
    for (size_t i = 1; i + 1 < m; i++) {
      if (orient_dev(v[0], v[i], v[i + 1]) < 0) {
        whole = 0;
        break;
      }
    }
    if (whole) {
      for (size_t k = 0; k < npiv && r == GSHGEOM_OK; k++) {
        GSHGEOM_SK_Vec pv = gshgeom_sk_add(anchor, gshgeom_sk_lin(c, piv[k]));
        GSHGEOM_Device_Point p;
        r = snap_round(c, pv.x, pv.y, &p);
        if (r == GSHGEOM_OK) {
          ring_insert(v, &m, p);
        }
      }
      if (r == GSHGEOM_OK) {
        r = write_ring(c->part, v, m);
      }
    } else {
      /* rounding turned the hull's fan from its first vertex negative somewhere: each triangle of it, turned positive */
      for (size_t i = 1; i + 1 < m && r == GSHGEOM_OK; i++) {
        int sg = orient_dev(v[0], v[i], v[i + 1]);
        if (sg == 0) {
          continue;
        }
        GSHGEOM_Device_Point t[3] = {v[0], sg > 0 ? v[i] : v[i + 1], sg > 0 ? v[i + 1] : v[i]};
        r = write_ring(c->part, t, 3);
      }
    }
  }
  if (v != stack) {
    gcu_allocator_free(&c->alloc, v);
  }
  return r;
}
