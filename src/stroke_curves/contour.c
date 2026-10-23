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
 * The pieces of the curve outline that are exact (stroke_curves.h, SG37): contours of lines, cubics and rational
 * quadratic arcs, written to the part in target units.
 *
 * ## The arcs
 *
 * The arc of radius `h` about a centre from unit `p` to unit `q` (an angle `psi` of at most a quarter turn, so
 * `p . q >= 0`) is the rational quadratic with control points `h p`, `h (p + q) / (1 + p . q)` (where the tangents
 * at the ends meet) and `h q`, weights `(1, cos(psi / 2), 1)`, `cos(psi / 2) = sqrt((1 + p . q) / 2)`: the one
 * square root of the construction, by `gmath_sqrt`.  An angle between a quarter turn and a half is cut at the
 * bisector `p + q` (or, near a half turn, where the sum loses its direction, at a quarter turn from `p` in the arc's
 * own sense).  The control points are offsets in the build space, and reach target units as `A + lin d`, `A` the
 * centre's position: an affine map takes a rational quadratic to a rational quadratic with the same weights, so
 * under a non-uniform map the arc is an ellipse arc and costs no approximation.
 *
 * ## The orientation
 *
 * Every contour is built in the build space's positive orientation (stroke.h: the boundary goes from +x toward +y
 * around the inside), so that its winding is +1 inside.  A map with a negative determinant reverses that in target
 * units; `flip` then writes the contour backwards (each element with its control points reversed), which restores
 * it.  A sector is built from the unit vector whose rotation to the other is positive.
 */

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/math/basic.h>

#include "../core/grow.h"
#include "../stroke/stroke_internal.h"

static GSHGEOM_Point gp(GSHGEOM_SK_Vec v) {
  GSHGEOM_Point p = {v.x, v.y, 0.0};
  return p;
}

void gshgeom_skc_begin(GSHGEOM_SK_Ctx * c) {
  c->cv->el_n = 0;
}

static GSHGEOM_Result push(GSHGEOM_SK_Ctx * c, const GSHGEOM_SK_CElem * e) {
  GSHGEOM_SK_Curves * cv = c->cv;
  void * g = cv->el;
  GSHGEOM_Result r = gshgeom_grow(&c->alloc, &g, &cv->el_cap, cv->el_n + 1, sizeof *cv->el);
  if (r != GSHGEOM_OK) {
    return r;
  }
  cv->el = g;
  cv->el[cv->el_n++] = *e;
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_skc_line(GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec a, GSHGEOM_SK_Vec b) {
  if (a.x == b.x && a.y == b.y) {
    return GSHGEOM_OK;
  }
  GSHGEOM_SK_CElem e;
  e.deg = 1;
  e.p[0] = a;
  e.p[1] = b;
  e.p[2] = e.p[3] = b;
  e.w[0] = e.w[1] = e.w[2] = 1.0;
  return push(c, &e);
}

GSHGEOM_Result gshgeom_skc_cubic(GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec p0, GSHGEOM_SK_Vec p1, GSHGEOM_SK_Vec p2,
    GSHGEOM_SK_Vec p3) {
  GSHGEOM_SK_CElem e;
  e.deg = 3;
  e.p[0] = p0;
  e.p[1] = p1;
  e.p[2] = p2;
  e.p[3] = p3;
  e.w[0] = e.w[1] = e.w[2] = 1.0;
  return push(c, &e);
}

GSHGEOM_Result gshgeom_skc_conic(GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec p0, GSHGEOM_SK_Vec p1, GSHGEOM_SK_Vec p2, double w) {
  GSHGEOM_SK_CElem e;
  e.deg = 2;
  e.p[0] = p0;
  e.p[1] = p1;
  e.p[2] = p2;
  e.p[3] = p2;
  e.w[0] = 1.0;
  e.w[1] = w;
  e.w[2] = 1.0;
  return push(c, &e);
}

static GSHGEOM_Result write_el(GSHGEOM_Path * path, const GSHGEOM_SK_CElem * e, int reversed) {
  GSHGEOM_SK_Vec q[4];
  double w[3];
  for (int i = 0; i <= e->deg; i++) {
    q[i] = reversed ? e->p[e->deg - i] : e->p[i];
  }
  switch (e->deg) {
  case 1:
    return gshgeom_path_line_to(path, gp(q[1]));
  case 2:
    w[0] = e->w[reversed ? 2 : 0];
    w[1] = e->w[1];
    w[2] = e->w[reversed ? 0 : 2];
    return gshgeom_path_rational_quadratic_to(path, gp(q[1]), gp(q[2]), w);
  default:
    return gshgeom_path_cubic_to(path, gp(q[1]), gp(q[2]), gp(q[3]));
  }
}

GSHGEOM_Result gshgeom_skc_finish(GSHGEOM_SK_Ctx * c) {
  GSHGEOM_SK_Curves * cv = c->cv;
  size_t n = cv->el_n;
  cv->el_n = 0;
  if (n < 2) {
    return GSHGEOM_OK;
  }
  GSHGEOM_Result r;
  if (!cv->flip) {
    r = gshgeom_path_move_to(cv->part, gp(cv->el[0].p[0]));
    for (size_t i = 0; i < n && r == GSHGEOM_OK; i++) {
      r = write_el(cv->part, &cv->el[i], 0);
    }
  } else {
    r = gshgeom_path_move_to(cv->part, gp(cv->el[n - 1].p[cv->el[n - 1].deg]));
    for (size_t i = n; i-- > 0 && r == GSHGEOM_OK;) {
      r = write_el(cv->part, &cv->el[i], 1);
    }
  }
  if (r == GSHGEOM_OK) {
    r = gshgeom_path_close(cv->part);
  }
  return r;
}

GSHGEOM_Result gshgeom_skc_ring(GSHGEOM_SK_Ctx * c, const GSHGEOM_SK_Vec * pts, size_t n) {
  if (n < 3) {
    return GSHGEOM_OK;
  }
  GSHGEOM_Path * path = c->cv->part;
  GSHGEOM_Result r = gshgeom_path_move_to(path, gp(pts[0]));
  for (size_t i = 1; i < n && r == GSHGEOM_OK; i++) {
    r = gshgeom_path_line_to(path, gp(pts[i]));
  }
  if (r == GSHGEOM_OK) {
    r = gshgeom_path_close(path);
  }
  return r;
}

/** The target point of the build-space offset `h v` from the centre. */
static GSHGEOM_SK_Vec dev(const GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec centre, double k, GSHGEOM_SK_Vec v) {
  return gshgeom_sk_add(centre, gshgeom_sk_lin(c, gshgeom_sk_scale(k, v)));
}

/** Append the arcs from unit @p p to unit @p q, the short way, about @p centre at radius @p rad. */
static GSHGEOM_Result arcs(GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec centre, double rad, GSHGEOM_SK_Vec p, GSHGEOM_SK_Vec q) {
  double cpq = gshgeom_sk_dot(p, q);
  if (cpq < 0.0) {
    GSHGEOM_SK_Vec m;
    if (cpq > -0.5) {
      m = gshgeom_sk_unit(gshgeom_sk_add(p, q));
    } else {
      m = gshgeom_sk_cross(p, q) >= 0.0 ? gshgeom_sk_right(p) : gshgeom_sk_left(p);
    }
    GSHGEOM_Result r = arcs(c, centre, rad, p, m);
    return r != GSHGEOM_OK ? r : arcs(c, centre, rad, m, q);
  }
  double h = rad;
  GSHGEOM_SK_Vec p0 = dev(c, centre, h, p), p2 = dev(c, centre, h, q);
  GSHGEOM_SK_Vec p1 = dev(c, centre, h / (1.0 + cpq), gshgeom_sk_add(p, q));
  return gshgeom_skc_conic(c, p0, p1, p2, gmath_sqrt((1.0 + cpq) / 2.0));
}

GSHGEOM_Result gshgeom_skc_sector(GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec centre, GSHGEOM_SK_Vec p, GSHGEOM_SK_Vec q) {
  double cr = gshgeom_sk_cross(p, q);
  if (cr == 0.0 && gshgeom_sk_dot(p, q) > 0.0) {
    return GSHGEOM_OK; /* no angle, no area */
  }
  if (cr < 0.0) {
    GSHGEOM_SK_Vec t = p;
    p = q;
    q = t;
  }
  gshgeom_skc_begin(c);
  GSHGEOM_Result r = gshgeom_skc_line(c, centre, dev(c, centre, c->h, p));
  if (r == GSHGEOM_OK) {
    r = arcs(c, centre, c->h, p, q);
  }
  if (r == GSHGEOM_OK) {
    r = gshgeom_skc_line(c, dev(c, centre, c->h, q), centre);
  }
  return r != GSHGEOM_OK ? r : gshgeom_skc_finish(c);
}

GSHGEOM_Result gshgeom_skc_half_disc(GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec centre, GSHGEOM_SK_Vec m, GSHGEOM_SK_Vec d) {
  if (gshgeom_sk_cross(m, d) < 0.0) {
    m = gshgeom_sk_scale(-1.0, m);
  }
  GSHGEOM_SK_Vec nm = gshgeom_sk_scale(-1.0, m);
  gshgeom_skc_begin(c);
  GSHGEOM_Result r = arcs(c, centre, c->h, m, d);
  if (r == GSHGEOM_OK) {
    r = arcs(c, centre, c->h, d, nm);
  }
  if (r == GSHGEOM_OK) {
    r = gshgeom_skc_line(c, dev(c, centre, c->h, nm), dev(c, centre, c->h, m));
  }
  return r != GSHGEOM_OK ? r : gshgeom_skc_finish(c);
}

GSHGEOM_Result gshgeom_skc_fan_sector(GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec centre, GSHGEOM_SK_Vec p, GSHGEOM_SK_Vec q,
    double r) {
  double cr = gshgeom_sk_cross(p, q);
  if (cr < 0.0) {
    GSHGEOM_SK_Vec t = p;
    p = q;
    q = t;
  }
  double cpq = gshgeom_sk_dot(p, q);
  if (cpq < 0.0) {
    /* more than a quarter turn: the union of the inflations of the two halves contains the inflation of the whole */
    GSHGEOM_SK_Vec m = cpq > -0.5 ? gshgeom_sk_unit(gshgeom_sk_add(p, q)) : gshgeom_sk_right(p);
    GSHGEOM_Result res = gshgeom_skc_fan_sector(c, centre, p, m, r);
    return res != GSHGEOM_OK ? res : gshgeom_skc_fan_sector(c, centre, m, q, r);
  }
  double h = c->h;
  GSHGEOM_SK_Vec u = gshgeom_sk_unit(gshgeom_sk_add(p, q));
  GSHGEOM_SK_Vec na = gshgeom_sk_left(p), nb = gshgeom_sk_right(q);
  GSHGEOM_SK_Vec oa = gshgeom_sk_scale(r, na), ob = gshgeom_sk_scale(r, nb);
  /* the offset side lines meet the back line {x : x . u = -r} */
  double la = (-r - gshgeom_sk_dot(oa, u)) / gshgeom_sk_dot(p, u);
  double lb = (-r - gshgeom_sk_dot(ob, u)) / gshgeom_sk_dot(q, u);
  GSHGEOM_SK_Vec a0 = gshgeom_sk_add(oa, gshgeom_sk_scale(la, p)), b0 = gshgeom_sk_add(ob, gshgeom_sk_scale(lb, q));
  /* and the circle of radius h + r about the apex, where |h p + r| ... = sqrt(h^2 + 2 h r) along the side */
  double sl = gmath_sqrt(h * h + 2.0 * h * r);
  GSHGEOM_SK_Vec a1 = gshgeom_sk_add(oa, gshgeom_sk_scale(sl, p)), b1 = gshgeom_sk_add(ob, gshgeom_sk_scale(sl, q));
  double rad = h + r;
  GSHGEOM_SK_Vec ua = gshgeom_sk_scale(1.0 / rad, a1), ub = gshgeom_sk_scale(1.0 / rad, b1);
  gshgeom_skc_begin(c);
  GSHGEOM_Result res = gshgeom_skc_line(c, dev(c, centre, 1.0, a0), dev(c, centre, 1.0, a1));
  if (res == GSHGEOM_OK) {
    res = arcs(c, centre, rad, ua, ub);
  }
  if (res == GSHGEOM_OK) {
    res = gshgeom_skc_line(c, dev(c, centre, 1.0, b1), dev(c, centre, 1.0, b0));
  }
  if (res == GSHGEOM_OK) {
    res = gshgeom_skc_line(c, dev(c, centre, 1.0, b0), dev(c, centre, 1.0, a0));
  }
  return res != GSHGEOM_OK ? res : gshgeom_skc_finish(c);
}

GSHGEOM_Result gshgeom_skc_disc(GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec centre, double radius) {
  static const double ax[5][2] = {{1, 0}, {0, 1}, {-1, 0}, {0, -1}, {1, 0}};
  gshgeom_skc_begin(c);
  GSHGEOM_Result r = GSHGEOM_OK;
  for (int k = 0; k < 4 && r == GSHGEOM_OK; k++) {
    r = arcs(c, centre, radius, gshgeom_sk_v(ax[k][0], ax[k][1]), gshgeom_sk_v(ax[k + 1][0], ax[k + 1][1]));
  }
  return r != GSHGEOM_OK ? r : gshgeom_skc_finish(c);
}

GSHGEOM_Result gshgeom_skc_add_band(GSHGEOM_SK_Ctx * c, const GSHGEOM_SK_Band * b) {
  GSHGEOM_SK_Curves * cv = c->cv;
  void * g = cv->run;
  GSHGEOM_Result r = gshgeom_grow(&c->alloc, &g, &cv->run_cap, cv->run_n + 1, sizeof *cv->run);
  if (r != GSHGEOM_OK) {
    return r;
  }
  cv->run = g;
  cv->run[cv->run_n] = *b;
  if (cv->run_n > 0) {
    /* a point two pieces share is one point (SG37): this piece starts where the one before ended */
    cv->run[cv->run_n].l[0] = cv->run[cv->run_n - 1].l[3];
    cv->run[cv->run_n].r[0] = cv->run[cv->run_n - 1].r[3];
  }
  cv->run_n++;
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_skc_flush_run(GSHGEOM_SK_Ctx * c) {
  GSHGEOM_SK_Curves * cv = c->cv;
  size_t n = cv->run_n;
  cv->run_n = 0;
  if (n == 0) {
    return GSHGEOM_OK;
  }
  const GSHGEOM_SK_Band * b = cv->run;
  gshgeom_skc_begin(c);
  GSHGEOM_Result r = GSHGEOM_OK;
  for (size_t i = 0; i < n && r == GSHGEOM_OK; i++) {
    r = gshgeom_skc_cubic(c, b[i].l[0], b[i].l[1], b[i].l[2], b[i].l[3]);
  }
  if (r == GSHGEOM_OK) {
    r = gshgeom_skc_line(c, b[n - 1].l[3], b[n - 1].r[3]);
  }
  for (size_t i = n; i-- > 0 && r == GSHGEOM_OK;) {
    r = gshgeom_skc_cubic(c, b[i].r[3], b[i].r[2], b[i].r[1], b[i].r[0]);
  }
  if (r == GSHGEOM_OK) {
    r = gshgeom_skc_line(c, b[0].r[0], b[0].l[0]);
  }
  return r != GSHGEOM_OK ? r : gshgeom_skc_finish(c);
}

/* --- the arcs join, exactly where its boundary is a polar graph --------------------------------------------- */

/** One side of the arcs join's boundary (build space, offsets from the vertex): a line, or a circle (z, r) traversed in sense s. */
typedef struct Side {
  int line;
  GSHGEOM_SK_Vec from, to;
  GSHGEOM_SK_Vec z;
  double r, s;
} Side;

/** The unit vector from the circle's centre to a point. */
static GSHGEOM_SK_Vec toward(const Side * e, GSHGEOM_SK_Vec p) {
  return gshgeom_sk_unit(gshgeom_sk_sub(p, e->z));
}

/** The quarter-turn boundaries between two units in the sense s: ends[0] = p, ..., ends[n-1] = q (n at most 6). */
static int quarter_ends(GSHGEOM_SK_Vec p, GSHGEOM_SK_Vec q, double s, GSHGEOM_SK_Vec * ends) {
  int n = 0;
  ends[n++] = p;
  GSHGEOM_SK_Vec cur = p;
  for (int k = 0; k < 4; k++) {
    if (s * gshgeom_sk_cross(cur, q) >= 0.0 && gshgeom_sk_dot(cur, q) >= 0.0) {
      break;
    }
    cur = s > 0 ? gshgeom_sk_right(cur) : gshgeom_sk_left(cur);
    ends[n++] = cur;
  }
  ends[n++] = q;
  return n;
}

/** Whether the unit d lies on the arc from p to q (at most a quarter turn) in the sense s. */
static int inside_arc(GSHGEOM_SK_Vec p, GSHGEOM_SK_Vec q, double s, GSHGEOM_SK_Vec d) {
  return s * gshgeom_sk_cross(p, d) >= 0.0 && s * gshgeom_sk_cross(d, q) >= 0.0;
}

/**
 * Whether the loop's polar angle about the vertex runs one way, strictly, along one side, in the sense @p sigma (the
 * sign of cross(A, B)), and every point stays in the wedge between A and B.  For a line cross(P, dP) is constant; for
 * a circle (z, r) in the sense s it is s r (w . u + r) at the point z + r u (w = z), and w . u + r is a cosine plus a
 * constant, extreme at the ends of a quarter arc or at the direction of w.
 */
static int side_monotone(const Side * e, double sigma, GSHGEOM_SK_Vec A, GSHGEOM_SK_Vec B) {
  GSHGEOM_SK_Vec pts[8];
  int np = 0;
  if (e->line) {
    double cr = gshgeom_sk_cross(e->from, gshgeom_sk_sub(e->to, e->from));
    if (!(sigma * cr > 0.0)) {
      return 0;
    }
    pts[np++] = e->from;
    pts[np++] = e->to;
  } else {
    GSHGEOM_SK_Vec ends[8];
    int ne = quarter_ends(toward(e, e->from), toward(e, e->to), e->s, ends);
    double wl = gshgeom_sk_norm(e->z);
    GSHGEOM_SK_Vec wh = wl > 0.0 ? gshgeom_sk_scale(1.0 / wl, e->z) : gshgeom_sk_v(0.0, 0.0);
    for (int k = 0; k + 1 < ne; k++) {
      GSHGEOM_SK_Vec p = ends[k], q = ends[k + 1];
      double fp = gshgeom_sk_dot(e->z, p) + e->r, fq = gshgeom_sk_dot(e->z, q) + e->r;
      double lo = fp < fq ? fp : fq, hi = fp < fq ? fq : fp;
      if (wl > 0.0 && inside_arc(p, q, e->s, gshgeom_sk_scale(-1.0, wh))) {
        lo = e->r - wl < lo ? e->r - wl : lo;
      }
      if (wl > 0.0 && inside_arc(p, q, e->s, wh)) {
        hi = e->r + wl > hi ? e->r + wl : hi;
      }
      double margin = 0x1p-30 * e->r;
      if (e->s * sigma > 0.0 ? !(lo > margin) : !(hi < -margin)) {
        return 0;
      }
      pts[np++] = gshgeom_sk_add(e->z, gshgeom_sk_scale(e->r, p));
    }
    pts[np++] = e->to;
  }
  for (int k = 0; k < np; k++) {
    if (!(sigma * gshgeom_sk_cross(A, pts[k]) >= 0.0) || !(sigma * gshgeom_sk_cross(pts[k], B) >= 0.0) ||
        !(gshgeom_sk_dot(pts[k], pts[k]) > 0.0)) {
      return 0;
    }
  }
  return 1;
}

static void reverse_elements(GSHGEOM_SK_Curves * cv) {
  size_t n = cv->el_n;
  for (size_t i = 0; i < n / 2; i++) {
    GSHGEOM_SK_CElem t = cv->el[i];
    cv->el[i] = cv->el[n - 1 - i];
    cv->el[n - 1 - i] = t;
  }
  for (size_t i = 0; i < n; i++) {
    GSHGEOM_SK_CElem * e = &cv->el[i];
    for (int k = 0; k < (e->deg + 1) / 2; k++) {
      GSHGEOM_SK_Vec t = e->p[k];
      e->p[k] = e->p[e->deg - k];
      e->p[e->deg - k] = t;
    }
  }
}

static GSHGEOM_Result side_elements(GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec anchor, const Side * e) {
  if (e->line) {
    return gshgeom_skc_line(c, dev(c, anchor, 1.0, e->from), dev(c, anchor, 1.0, e->to));
  }
  GSHGEOM_SK_Vec ends[8];
  int ne = quarter_ends(toward(e, e->from), toward(e, e->to), e->s, ends);
  GSHGEOM_Result r = GSHGEOM_OK;
  for (int k = 0; k + 1 < ne && r == GSHGEOM_OK; k++) {
    GSHGEOM_SK_Vec p = ends[k], q = ends[k + 1];
    double cpq = gshgeom_sk_dot(p, q);
    GSHGEOM_SK_Vec p0 = k == 0 ? e->from : gshgeom_sk_add(e->z, gshgeom_sk_scale(e->r, p));
    GSHGEOM_SK_Vec p2 = k + 2 == ne ? e->to : gshgeom_sk_add(e->z, gshgeom_sk_scale(e->r, q));
    GSHGEOM_SK_Vec p1 = gshgeom_sk_add(e->z, gshgeom_sk_scale(e->r / (1.0 + cpq), gshgeom_sk_add(p, q)));
    r = gshgeom_skc_conic(c, dev(c, anchor, 1.0, p0), dev(c, anchor, 1.0, p1), dev(c, anchor, 1.0, p2),
        gmath_sqrt((1.0 + cpq) / 2.0));
  }
  return r;
}

GSHGEOM_Result gshgeom_skc_arcs_join(GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec anchor, GSHGEOM_SK_Vec A, GSHGEOM_SK_Vec B,
    GSHGEOM_SK_Vec X, int line_a, GSHGEOM_SK_Vec za, double ra, double sa, int line_b, GSHGEOM_SK_Vec zb, double rb,
    double sb, int * done) {
  *done = 0;
  double sigma = gshgeom_sk_cross(A, B) > 0.0 ? 1.0 : (gshgeom_sk_cross(A, B) < 0.0 ? -1.0 : 0.0);
  if (sigma == 0.0) {
    return GSHGEOM_OK;
  }
  GSHGEOM_SK_Vec zero = gshgeom_sk_v(0.0, 0.0);
  Side ea = {line_a, A, X, za, ra, sa}, eb = {line_b, X, B, zb, rb, sb};
  Side va = {1, zero, A, zero, 0.0, 0.0}, vb = {1, B, zero, zero, 0.0, 0.0};
  (void)va;
  (void)vb;
  if (!side_monotone(&ea, sigma, A, B) || !side_monotone(&eb, sigma, A, B)) {
    return GSHGEOM_OK;
  }
  gshgeom_skc_begin(c);
  GSHGEOM_Result r = gshgeom_skc_line(c, dev(c, anchor, 1.0, zero), dev(c, anchor, 1.0, A));
  if (r == GSHGEOM_OK) {
    r = side_elements(c, anchor, &ea);
  }
  if (r == GSHGEOM_OK) {
    r = side_elements(c, anchor, &eb);
  }
  if (r == GSHGEOM_OK) {
    r = gshgeom_skc_line(c, dev(c, anchor, 1.0, B), dev(c, anchor, 1.0, zero));
  }
  if (r != GSHGEOM_OK) {
    return r;
  }
  if (sigma < 0.0) {
    reverse_elements(c->cv);
  }
  r = gshgeom_skc_finish(c);
  if (r == GSHGEOM_OK) {
    *done = 1;
  }
  return r;
}
