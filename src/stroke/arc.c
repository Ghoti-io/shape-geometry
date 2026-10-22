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
 * Circular arcs as pieces (round joins and caps, discs, sectors), and the arcs
 * join of SVG 2's Candidate Recommendation of 2018-10-04 (stroke.h, SG8).
 *
 * ## An arc
 *
 * The arc of radius `r` from unit `p` to unit `q` (angle `psi` below pi) is the
 * rational quadratic with control points `r p`, `r (p + q) / (1 + p . q)` (the
 * tangents' intersection) and `r q` and middle weight `cos(psi / 2)`; it lies in
 * its control triangle (positive weights), and the triangle's apex is
 * `r (1 / cos(psi / 2) - 1) = r (sqrt(2 / (1 + p . q)) - 1)` outside it.  An arc
 * is halved (de Casteljau at 1/2 of that rational quadratic: the new end is
 * `(p + q) / |p + q|`) until that distance, in device pixels (times the map's
 * norm under local width), is within what the piece may spend.  The polygon of
 * the centre and the control points of the halves contains the sector (the
 * control triangles contain the arcs, and the fan from the centre contains the
 * rest) and is within that distance of it: a piece with `d_out` that distance.
 * No angle is formed and no trigonometric function is called.
 *
 * ## The arcs join
 *
 * The join of SG8 is built as the CR states, each gap settled as SG8 says:
 * curvatures from the end control points (exact zero and `|kappa| > 1 / h`
 * decided exactly), the circles concentric with the centres of curvature, the
 * meeting point closest to the vertex, the radii adjusted by the least amount that
 * makes them touch, the CR's rectangle when nothing can meet.  Its region is the
 * bevel together with what the closed boundary V, A, A's extension to X, B's
 * extension back to B, V encloses by nonzero winding.  The boundary is a set of
 * directed edges: the straight ones, and for each extension the circumscribed
 * polylines of its sub-arcs of at most a quarter turn (gshgeom_sk_arc), each
 * within `dev` of its arc.  The clip (`s > L h`) keeps `(x - Y) . T <= 0`; it is
 * applied to the boundary, not to the pieces: the parts of the edges in the half
 * plane are kept, and the gaps between them are closed along the clip line.  A
 * closed curve and its clipping wind the same about every point of the half
 * plane and the clipping winds zero about every other, so the clipped boundary
 * encloses exactly the region cut by the half plane; and a stretch of the clip
 * line two chords cover in opposite directions carries nothing, so those cancel
 * before anything is written (a perturbation of them would leave a sliver far
 * from the region).  The polygon of the edges winds like the exact boundary about
 * every point farther than `dev` from it (the straight-line homotopy between
 * them never reaches the point), so the region it encloses contains the exact
 * region without the points within `dev` of the boundary and lies within `dev`
 * of it; each piece is therefore emitted with `d_in = d_out = dev` plus the
 * arithmetic's reach.  The region is cut into slabs at every height where an
 * edge ends or two edges cross; within a slab the edges are ordered by x, the
 * winding between neighbours is the running sum of their directions, and each
 * stretch of nonzero winding is one trapezoid, so every piece is convex and
 * the pieces are positively wound (emit.c).  The clip's point `Y` and the arc
 * length `s` are the one place math's `gmath_atan2`, `gmath_sin` and `gmath_cos`
 * are called (SG18); they are constructions, and the region is continuous in
 * them.  A crossing of an arc with the clip line is found from the line's normal
 * and placed along its sub-arc by the signs of cross products: no angle.
 */

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/math/basic.h>
#include <ghoti.io/math/elementary.h>
#include <ghoti.io/shape-geometry/orient.h>

#include "../core/grow.h"
#include "../exact/exact.h"
#include "stroke_internal.h"

#include <string.h>

static double fabs_d(double x) {
  return x < 0.0 ? -x : x;
}

/** How much a computed bound is raised for the rounding of its own computation. */
#define UP (1.0 + 0x1p-40)

/** A sub-arc on the work stack: unit ends. */
typedef struct Sub {
  GSHGEOM_SK_Vec p, q;
} Sub;

/** Append one point to the arc buffer. */
static GSHGEOM_Result put(GSHGEOM_SK_Ctx * c, size_t * n, GSHGEOM_SK_Vec v) {
  size_t most = c->limits.max_points;
  if (c->arc_limit && c->arc_limit < most) {
    most = c->arc_limit;
  }
  if (*n >= most) {
    return GSHGEOM_ERR_LIMIT; /* the part could not hold the piece */
  }
  void * g = c->arcbuf;
  GSHGEOM_Result r = gshgeom_grow(&c->alloc, &g, &c->arcbuf_cap, *n + 1, sizeof *c->arcbuf);
  if (r != GSHGEOM_OK) {
    return r;
  }
  c->arcbuf = g;
  c->arcbuf[(*n)++] = v;
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_sk_arc(GSHGEOM_SK_Ctx * c, double r, GSHGEOM_SK_Vec p, GSHGEOM_SK_Vec q, double allowed,
    size_t * n, double * dev) {
  Sub stack[2 * 64 + 4];
  int depth[2 * 64 + 4];
  size_t top = 0;
  double worst = 0.0;
  stack[top].p = p;
  stack[top].q = q;
  depth[top] = 0;
  top++;
  GSHGEOM_Result res = put(c, n, gshgeom_sk_scale(r, p));
  while (top > 0 && res == GSHGEOM_OK) {
    top--;
    Sub s = stack[top];
    int dp = depth[top];
    double cpq = gshgeom_sk_dot(s.p, s.q);
    /* the apex's distance outside the arc, in device pixels */
    double d = (cpq > -0.5) ? r * c->norm * (gmath_sqrt(2.0 / (1.0 + cpq)) - 1.0) * UP : -1.0;
    if (d >= 0.0 && d <= allowed) {
      double k = r / (1.0 + cpq);
      res = put(c, n, gshgeom_sk_scale(k, gshgeom_sk_add(s.p, s.q)));
      if (res == GSHGEOM_OK) {
        res = put(c, n, gshgeom_sk_scale(r, s.q));
      }
      worst = d > worst ? d : worst;
      continue;
    }
    if (dp >= 64) {
      return GSHGEOM_ERR_LIMIT; /* halved 64 times and still too far: no double can say more */
    }
    /* halve: the new end is the bisector; push the second half first so the first is taken first */
    GSHGEOM_SK_Vec m = gshgeom_sk_add(s.p, s.q);
    if (cpq <= -0.5) {
      /* near a half turn the sum loses its direction: take the quarter turn in the arc's own sense */
      m = gshgeom_sk_cross(s.p, s.q) >= 0.0 ? gshgeom_sk_right(s.p) : gshgeom_sk_left(s.p);
      /* and then halve toward q from there if that quarter overshoots */
    }
    m = gshgeom_sk_unit(m);
    stack[top].p = m;
    stack[top].q = s.q;
    depth[top] = dp + 1;
    top++;
    stack[top].p = s.p;
    stack[top].q = m;
    depth[top] = dp + 1;
    top++;
  }
  *dev = worst;
  return res;
}

GSHGEOM_Result gshgeom_sk_sector(GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec anchor_dev, GSHGEOM_SK_Vec p, GSHGEOM_SK_Vec q,
    double d_in, double extra_out) {
  double allowed = gshgeom_sk_room(c, d_in) - extra_out;
  if (!(allowed > 0.0)) {
    return GSHGEOM_ERR_INTERNAL;
  }
  size_t n = 0;
  double dev = 0.0;
  GSHGEOM_Result r = put(c, &n, gshgeom_sk_v(0.0, 0.0));
  if (r == GSHGEOM_OK) {
    r = gshgeom_sk_arc(c, c->h, p, q, allowed, &n, &dev);
  }
  if (r == GSHGEOM_OK) {
    r = gshgeom_sk_emit(c, anchor_dev, c->arcbuf, n, d_in, dev + extra_out);
  }
  return r;
}

GSHGEOM_Result gshgeom_sk_disc(GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec anchor_dev, double d_in, double extra_out) {
  double allowed = gshgeom_sk_room(c, d_in) - extra_out;
  if (!(allowed > 0.0)) {
    return GSHGEOM_ERR_INTERNAL;
  }
  static const double axes[5][2] = {{1, 0}, {0, 1}, {-1, 0}, {0, -1}, {1, 0}};
  size_t n = 0;
  double worst = 0.0;
  GSHGEOM_Result r = GSHGEOM_OK;
  for (int k = 0; k < 4 && r == GSHGEOM_OK; k++) {
    double dev = 0.0;
    if (k > 0) {
      n--; /* each quarter starts where the last ended */
    }
    r = gshgeom_sk_arc(c, c->h, gshgeom_sk_v(axes[k][0], axes[k][1]), gshgeom_sk_v(axes[k + 1][0], axes[k + 1][1]),
        allowed, &n, &dev);
    worst = dev > worst ? dev : worst;
  }
  if (r == GSHGEOM_OK) {
    r = gshgeom_sk_emit(c, anchor_dev, c->arcbuf, n, d_in, worst + extra_out);
  }
  return r;
}

/* --- the arcs join ----------------------------------------------------------------- */

/**
 * `|kappa| > 1 / h` at an end, decided exactly (SG8):
 * `(n-1)^2 h^2 w0^2 w2^2 X^2 > n^2 w1^4 |d1|^6` with `X = d1 x d2`, every product
 * expanded into the raw coordinates.  @p q holds the three control points from
 * the end, @p w their weights.
 */
static int exceeds_inverse_h(int n, const GSHGEOM_SK_Vec q[3], const double w[3], double h) {
  /* X = (q1 - q0) x (q2 - q1) for a start; for an end the caller passes the points from the end and the sign does
   * not matter (it is squared).  Expanded: q1x q2y - q1x q1y - q0x q2y + q0x q1y - q1y q2x + q1y q1x + q0y q2x - q0y
   * q1x = q1x q2y - q0x q2y + q0x q1y - q1y q2x + q0y q2x - q0y q1x. */
  double xs[6][2] = {{q[1].x, q[2].y}, {q[0].x, q[2].y}, {q[0].x, q[1].y}, {q[1].y, q[2].x}, {q[0].y, q[2].x},
      {q[0].y, q[1].x}};
  int xneg[6] = {0, 1, 0, 1, 0, 1};
  /* |d1|^2 = q1x^2 - 2 q1x q0x + q0x^2 + (same in y): terms with a constant factor */
  double ds[6][2] = {{q[1].x, q[1].x}, {q[1].x, q[0].x}, {q[0].x, q[0].x}, {q[1].y, q[1].y}, {q[1].y, q[0].y},
      {q[0].y, q[0].y}};
  double dk[6] = {1, -2, 1, 1, -2, 1};
  GSHGEOM_ExactSum s;
  gshgeom_exact_sum_init(&s);
  double nm1 = (double)(n - 1), nn = (double)n;
  /* + (n-1)^2 h^2 w0^2 w2^2 X^2: 36 terms */
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      double f[11] = {nm1 * nm1, h, h, w[0], w[0], w[2], w[2], xs[i][0], xs[i][1], xs[j][0], xs[j][1]};
      gshgeom_exact_sum_add(&s, xneg[i] ^ xneg[j], 11, f);
    }
  }
  /* - n^2 w1^4 |d1|^6: 216 terms */
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      for (int k = 0; k < 6; k++) {
        double kc = dk[i] * dk[j] * dk[k];
        double f[11] = {nn * nn * (kc < 0 ? -kc : kc), w[1], w[1], w[1], w[1], ds[i][0], ds[i][1], ds[j][0], ds[j][1],
            ds[k][0], ds[k][1]};
        gshgeom_exact_sum_add(&s, kc < 0 ? 0 : 1, 11, f);
      }
    }
  }
  return gshgeom_exact_sum_sign(&s) > 0;
}

/** The curvature data of one end of a segment at a vertex. */
typedef struct EndCurv {
  int line;      ///< The extension is a line (a line segment, a vanishing first difference, or a zero curvature).
  int over;      ///< |kappa| > 1/h, decided exactly.
  double kappa;  ///< Signed: positive when the centre of curvature is on the inner side.
} EndCurv;

static EndCurv end_curvature(const GSHGEOM_SK_Seg * s, int at_end, int outer_left, double h) {
  EndCurv e = {1, 0, 0.0};
  if (s->n == 1) {
    return e;
  }
  int i0 = at_end ? s->n : 0, i1 = at_end ? s->n - 1 : 1, i2 = at_end ? s->n - 2 : 2;
  GSHGEOM_SK_Vec q[3] = {s->p[i0], s->p[i1], s->p[i2]};
  double w[3] = {s->w[i0], s->w[i1], s->w[i2]};
  /* a vanishing first difference (q1 = q0): the CR says the curvature is infinite and a line is used; the three
   * points are then collinear, so the zero turn below returns the line */
  /* the turning sense along travel, exactly: at a start q0 -> q1 -> q2; at an end q2 -> q1 -> q0 */
  GSHGEOM_Point a = {at_end ? q[2].x : q[0].x, at_end ? q[2].y : q[0].y, 0.0};
  GSHGEOM_Point b = {q[1].x, q[1].y, 0.0};
  GSHGEOM_Point cc = {at_end ? q[0].x : q[2].x, at_end ? q[0].y : q[2].y, 0.0};
  int turn = 0;
  (void)gshgeom_orient2d(a, b, cc, &turn);
  if (turn == 0) {
    return e; /* a zero curvature: a line */
  }
  e.line = 0;
  e.over = exceeds_inverse_h(s->n, q, w, h);
  GSHGEOM_SK_Vec d1 = gshgeom_sk_sub(q[1], q[0]);
  GSHGEOM_SK_Vec d2 = gshgeom_sk_sub(q[2], q[1]);
  double x = gshgeom_sk_cross(d1, d2);
  double nd = gshgeom_sk_norm(d1);
  double k = ((double)(s->n - 1) / (double)s->n) * (w[0] * w[2] / (w[1] * w[1])) * (x < 0 ? -x : x) / (nd * nd * nd);
  /* a positive turn puts the centre on the right of travel; the inner side is the right when the outer is the left */
  int inner_positive = (turn > 0) == (outer_left != 0);
  e.kappa = inner_positive ? k : -k;
  return e;
}

/** An extension: a line through `a` along unit `dir`, or the circle (z, r) through `a` with unit tangent `dir` there. */
typedef struct Ext {
  int line;
  GSHGEOM_SK_Vec a, dir;
  GSHGEOM_SK_Vec z;
  double r;
  GSHGEOM_SK_Vec nu; ///< unit from `a` toward the centre (for the adjustment)
} Ext;

static int circle_circle(const Ext * e1, const Ext * e2, GSHGEOM_SK_Vec out[2]) {
  GSHGEOM_SK_Vec dz = gshgeom_sk_sub(e2->z, e1->z);
  double d = gshgeom_sk_norm(dz);
  double r1 = e1->r, r2 = e2->r;
  if (d == 0.0 || d > r1 + r2 || d < (r1 > r2 ? r1 - r2 : r2 - r1)) {
    return 0;
  }
  double a = (r1 * r1 - r2 * r2 + d * d) / (2.0 * d);
  double hh = r1 * r1 - a * a;
  double hl = hh > 0.0 ? gmath_sqrt(hh) : 0.0;
  GSHGEOM_SK_Vec e = gshgeom_sk_scale(1.0 / d, dz);
  GSHGEOM_SK_Vec m = gshgeom_sk_add(e1->z, gshgeom_sk_scale(a, e));
  GSHGEOM_SK_Vec perp = gshgeom_sk_v(-e.y, e.x);
  out[0] = gshgeom_sk_add(m, gshgeom_sk_scale(hl, perp));
  out[1] = gshgeom_sk_sub(m, gshgeom_sk_scale(hl, perp));
  return 2;
}

static int circle_line(const Ext * circ, const Ext * line, GSHGEOM_SK_Vec out[2]) {
  GSHGEOM_SK_Vec m = gshgeom_sk_sub(line->a, circ->z);
  double bq = gshgeom_sk_dot(m, line->dir), cq = gshgeom_sk_dot(m, m) - circ->r * circ->r;
  double disc = bq * bq - cq;
  if (disc < 0.0) {
    return 0;
  }
  double sq = gmath_sqrt(disc);
  out[0] = gshgeom_sk_add(line->a, gshgeom_sk_scale(-bq - sq, line->dir));
  out[1] = gshgeom_sk_add(line->a, gshgeom_sk_scale(-bq + sq, line->dir));
  return 2;
}

/** The least positive root of q2 x^2 + q1 x + q0; -1 for none. */
static double least_positive(double q2, double q1, double q0) {
  double best = -1.0;
  double roots[2];
  int nr = 0;
  if (q2 == 0.0) {
    if (q1 != 0.0) {
      roots[nr++] = -q0 / q1;
    }
  } else {
    double disc = q1 * q1 - 4.0 * q2 * q0;
    if (disc < 0.0) {
      disc = 0.0; /* a double root the rounding pushed below zero: the circles graze */
    }
    double sq = gmath_sqrt(disc);
    /* the stable pair of roots */
    double t = q1 >= 0.0 ? -0.5 * (q1 + sq) : -0.5 * (q1 - sq);
    if (t != 0.0) {
      roots[nr++] = q0 / t;
    }
    roots[nr++] = t / q2;
  }
  for (int i = 0; i < nr; i++) {
    if (roots[i] > 0.0 && (best < 0.0 || roots[i] < best)) {
      best = roots[i];
    }
  }
  return best;
}

/** The sense of travel around an extension's circle from its start: +1 positive (orient.h), -1 negative. */
static double ext_sense(const Ext * e) {
  return gshgeom_sk_cross(gshgeom_sk_sub(e->a, e->z), e->dir) >= 0.0 ? 1.0 : -1.0;
}

/** The angle-free order of a point along a circle from the extension's start, as a monotone key in [0, 4). */
static double along_key(const Ext * e, GSHGEOM_SK_Vec p) {
  if (e->line) {
    return gshgeom_sk_dot(gshgeom_sk_sub(p, e->a), e->dir);
  }
  double s = ext_sense(e);
  GSHGEOM_SK_Vec u = gshgeom_sk_sub(e->a, e->z), v = gshgeom_sk_sub(p, e->z);
  double cr = s * gshgeom_sk_cross(u, v), dt = gshgeom_sk_dot(u, v);
  double nu = gshgeom_sk_norm(u) * gshgeom_sk_norm(v);
  if (nu == 0.0) {
    return 0.0;
  }
  double cs = dt / nu;
  /* pseudo-angle: 0 at a, increasing with the swept angle up to 4 at a full turn */
  return cr >= 0.0 ? 1.0 - cs : 3.0 + cs;
}

/** The clip of the arcs join: keep `(x - y) . t <= 0` (build space, offsets from V; `t` a unit vector). */
typedef struct Clip {
  int on;
  GSHGEOM_SK_Vec y, t;
} Clip;

/** Clip a convex polygon (build-space offsets from V) by the half plane; returns the new count. */
static size_t clip_poly(const GSHGEOM_SK_Vec * in, size_t n, GSHGEOM_SK_Vec * out, const Clip * cl) {
  if (!cl->on) {
    memcpy(out, in, n * sizeof *in);
    return n;
  }
  size_t m = 0;
  for (size_t i = 0; i < n; i++) {
    GSHGEOM_SK_Vec a = in[i], b = in[(i + 1) % n];
    double fa = gshgeom_sk_dot(gshgeom_sk_sub(a, cl->y), cl->t);
    double fb = gshgeom_sk_dot(gshgeom_sk_sub(b, cl->y), cl->t);
    if (fa <= 0.0) {
      out[m++] = a;
    }
    if ((fa < 0.0 && fb > 0.0) || (fa > 0.0 && fb < 0.0)) {
      double s = fa / (fa - fb);
      out[m++] = gshgeom_sk_add(a, gshgeom_sk_scale(s, gshgeom_sk_sub(b, a)));
    }
  }
  return m;
}

/* --- the region of the arcs join: a boundary, then slabs ------------------------------ */

/** The most directed edges one arcs join's boundary may have (ERR_LIMIT past it). */
#define BOUNDARY_EDGES 2048
/** The most gaps the clip line may cut in the boundary. */
#define BOUNDARY_CHORDS GSHGEOM_SK_CHORDS_MAX
/** The most slab boundaries one join may have. */
#define BOUNDARY_EVENTS (1u << 17)

typedef struct Edge {
  GSHGEOM_SK_Vec p, q;
} Edge;

typedef struct Span {
  double xm, x0, x1;
  int dir;
} Span;

typedef struct Boundary {
  GSHGEOM_SK_Ctx * c;
  Edge * e;
  size_t n, cap;
  GSHGEOM_SK_Vec * pts;
  size_t pts_cap;
  int have;
  GSHGEOM_SK_Vec first, last;
  GSHGEOM_SK_Vec chord[BOUNDARY_CHORDS][2];
  size_t nchords;
  Clip cl;
  double allowed; ///< device pixels: how far a polygon may lie outside the arc it replaces
  double worst;   ///< the largest bound the arcs reported
} Boundary;

static int same(GSHGEOM_SK_Vec a, GSHGEOM_SK_Vec b) {
  return a.x == b.x && a.y == b.y;
}

static double clip_f(const Boundary * b, GSHGEOM_SK_Vec p) {
  return gshgeom_sk_dot(gshgeom_sk_sub(p, b->cl.y), b->cl.t);
}

static GSHGEOM_Result add_edge(Boundary * b, GSHGEOM_SK_Vec p, GSHGEOM_SK_Vec q) {
  if (same(p, q)) {
    return GSHGEOM_OK;
  }
  if (b->n >= BOUNDARY_EDGES) {
    return GSHGEOM_ERR_LIMIT;
  }
  void * g = b->e;
  GSHGEOM_Result r = gshgeom_grow(&b->c->alloc, &g, &b->cap, b->n + 1, sizeof *b->e);
  if (r != GSHGEOM_OK) {
    return r;
  }
  b->e = g;
  b->e[b->n].p = p;
  b->e[b->n].q = q;
  b->n++;
  return GSHGEOM_OK;
}

static GSHGEOM_Result add_chord(Boundary * b, GSHGEOM_SK_Vec from, GSHGEOM_SK_Vec to) {
  if (b->nchords >= BOUNDARY_CHORDS) {
    return GSHGEOM_ERR_LIMIT;
  }
  b->chord[b->nchords][0] = from;
  b->chord[b->nchords][1] = to;
  b->nchords++;
  return GSHGEOM_OK;
}

/** Append one kept piece of the curve (a polyline); a gap from the last kept piece is a chord of the clip line. */
static GSHGEOM_Result add_piece(Boundary * b, const GSHGEOM_SK_Vec * pts, size_t np) {
  if (np < 2) {
    return GSHGEOM_OK;
  }
  GSHGEOM_Result r = GSHGEOM_OK;
  if (!b->have) {
    b->have = 1;
    b->first = pts[0];
  } else if (!same(b->last, pts[0])) {
    r = add_chord(b, b->last, pts[0]);
  }
  for (size_t i = 0; i + 1 < np && r == GSHGEOM_OK; i++) {
    r = add_edge(b, pts[i], pts[i + 1]);
  }
  b->last = pts[np - 1];
  return r;
}

static GSHGEOM_Result add_line(Boundary * b, GSHGEOM_SK_Vec p, GSHGEOM_SK_Vec q) {
  GSHGEOM_SK_Vec pq[2] = {p, q};
  if (!b->cl.on) {
    return add_piece(b, pq, 2);
  }
  double fp = clip_f(b, p), fq = clip_f(b, q);
  if (fp <= 0.0 && fq <= 0.0) {
    return add_piece(b, pq, 2);
  }
  if (fp > 0.0 && fq > 0.0) {
    return GSHGEOM_OK;
  }
  double t = fp / (fp - fq);
  GSHGEOM_SK_Vec x = gshgeom_sk_add(p, gshgeom_sk_scale(t, gshgeom_sk_sub(q, p)));
  if (fp <= 0.0) {
    pq[1] = x;
  } else {
    pq[0] = x;
  }
  return add_piece(b, pq, 2);
}

/**
 * The arc of the circle (z, r) from p0 to p1 in the travel sense s (+1 positive), as the circumscribed polylines of
 * its quarter-turn sub-arcs, those parts of it in the clip's half plane.  No angle is formed: a crossing of the clip
 * line is a unit vector found from the line's normal, and its place along a sub-arc (at most a quarter turn) is the
 * sign and size of cross products.
 */
static GSHGEOM_Result add_arc(Boundary * b, GSHGEOM_SK_Vec z, double r, GSHGEOM_SK_Vec p0, GSHGEOM_SK_Vec p1,
    double s) {
  GSHGEOM_SK_Ctx * c = b->c;
  GSHGEOM_SK_Vec u0 = gshgeom_sk_unit(gshgeom_sk_sub(p0, z));
  GSHGEOM_SK_Vec u1 = gshgeom_sk_unit(gshgeom_sk_sub(p1, z));
  GSHGEOM_SK_Vec ends[7];
  int ne = 0;
  ends[ne++] = u0;
  GSHGEOM_SK_Vec cur = u0;
  for (int k = 0; k < 4; k++) {
    double cr = s * gshgeom_sk_cross(cur, u1), dt = gshgeom_sk_dot(cur, u1);
    if (cr >= 0.0 && dt >= 0.0) {
      break;
    }
    cur = s > 0 ? gshgeom_sk_right(cur) : gshgeom_sk_left(cur);
    ends[ne++] = cur;
  }
  ends[ne++] = u1;
  GSHGEOM_Result res = GSHGEOM_OK;
  for (int k = 0; k + 1 < ne && res == GSHGEOM_OK; k++) {
    GSHGEOM_SK_Vec a = ends[k], q = ends[k + 1];
    GSHGEOM_SK_Vec bound[4];
    int nb = 0;
    bound[nb++] = a;
    if (b->cl.on) {
      double cc = -gshgeom_sk_dot(gshgeom_sk_sub(z, b->cl.y), b->cl.t) / r;
      if (cc > -1.0 && cc < 1.0) {
        double s2 = gmath_sqrt(1.0 - cc * cc);
        GSHGEOM_SK_Vec perp = gshgeom_sk_right(b->cl.t);
        GSHGEOM_SK_Vec cu[2];
        double key[2];
        int nc = 0;
        for (int sg = -1; sg <= 1; sg += 2) {
          GSHGEOM_SK_Vec u = gshgeom_sk_add(gshgeom_sk_scale(cc, b->cl.t), gshgeom_sk_scale((double)sg * s2, perp));
          if (s * gshgeom_sk_cross(a, u) >= 0.0 && s * gshgeom_sk_cross(u, q) >= 0.0) {
            cu[nc] = u;
            key[nc] = s * gshgeom_sk_cross(a, u);
            nc++;
          }
        }
        if (nc == 2 && key[1] < key[0]) {
          GSHGEOM_SK_Vec t = cu[0];
          cu[0] = cu[1];
          cu[1] = t;
        }
        for (int i = 0; i < nc; i++) {
          bound[nb++] = cu[i];
        }
      }
    }
    bound[nb++] = q;
    for (int j = 0; j + 1 < nb && res == GSHGEOM_OK; j++) {
      GSHGEOM_SK_Vec pa = bound[j], pb = bound[j + 1];
      if (same(pa, pb)) {
        continue;
      }
      if (b->cl.on) {
        GSHGEOM_SK_Vec mid = gshgeom_sk_add(pa, pb);
        double mn = gshgeom_sk_norm(mid);
        if (mn > 0.0) {
          GSHGEOM_SK_Vec m = gshgeom_sk_add(z, gshgeom_sk_scale(r / mn, mid));
          if (clip_f(b, m) > 0.0) {
            continue; /* this piece is on the far side of the clip line */
          }
        }
      }
      size_t n = 0;
      double dev = 0.0;
      c->arc_limit = BOUNDARY_EDGES - b->n + 1; /* a polyline the boundary could not hold is refused as it grows */
      res = gshgeom_sk_arc(c, r, pa, pb, b->allowed, &n, &dev);
      c->arc_limit = 0;
      if (res != GSHGEOM_OK) {
        break;
      }
      b->worst = dev > b->worst ? dev : b->worst;
      void * g = b->pts;
      res = gshgeom_grow(&c->alloc, &g, &b->pts_cap, n, sizeof *b->pts);
      if (res != GSHGEOM_OK) {
        break;
      }
      b->pts = g;
      for (size_t i = 0; i < n; i++) {
        b->pts[i] = gshgeom_sk_add(z, c->arcbuf[i]);
      }
      if (k == 0 && j == 0) {
        b->pts[0] = p0; /* the edge's own ends, not their rounding from the circle */
      }
      if (k + 2 == ne && j + 2 == nb) {
        b->pts[n - 1] = p1;
      }
      res = add_piece(b, b->pts, n);
    }
  }
  return res;
}

/** Sort doubles ascending: a merge sort, so the result does not depend on the C library. */
static void sort_doubles(double * a, double * tmp, size_t n) {
  for (size_t w = 1; w < n; w *= 2) {
    for (size_t lo = 0; lo < n; lo += 2 * w) {
      size_t mid = lo + w < n ? lo + w : n, hi = lo + 2 * w < n ? lo + 2 * w : n;
      size_t i = lo, j = mid, k = lo;
      while (i < mid && j < hi) {
        tmp[k++] = a[j] < a[i] ? a[j++] : a[i++];
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

static void sort_spans(Span * a, Span * tmp, size_t n) {
  for (size_t w = 1; w < n; w *= 2) {
    for (size_t lo = 0; lo < n; lo += 2 * w) {
      size_t mid = lo + w < n ? lo + w : n, hi = lo + 2 * w < n ? lo + 2 * w : n;
      size_t i = lo, j = mid, k = lo;
      while (i < mid && j < hi) {
        int take_j = a[j].xm < a[i].xm || (a[j].xm == a[i].xm && a[j].x0 + a[j].x1 < a[i].x0 + a[i].x1);
        tmp[k++] = take_j ? a[j++] : a[i++];
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

/** The x of the edge at height y: its own end where y is one, else interpolated and kept between its ends. */
static double edge_x(const Edge * e, double y) {
  if (y == e->p.y) {
    return e->p.x;
  }
  if (y == e->q.y) {
    return e->q.x;
  }
  double x = e->p.x + (y - e->p.y) * (e->q.x - e->p.x) / (e->q.y - e->p.y);
  double lo = e->p.x < e->q.x ? e->p.x : e->q.x, hi = e->p.x < e->q.x ? e->q.x : e->p.x;
  return x < lo ? lo : (x > hi ? hi : x);
}

GSHGEOM_Result gshgeom_sk_net_chords(const GSHGEOM_SK_Vec (*chord)[2], size_t k, GSHGEOM_SK_Vec y,
    GSHGEOM_SK_Vec t, GSHGEOM_SK_EdgeFn fn, void * ctx) {
  if (k == 0) {
    return GSHGEOM_OK;
  }
  if (k > GSHGEOM_SK_CHORDS_MAX) {
    return GSHGEOM_ERR_LIMIT;
  }
  GSHGEOM_SK_Vec perp = gshgeom_sk_right(t);
  GSHGEOM_SK_Vec pt[2 * GSHGEOM_SK_CHORDS_MAX];
  double tau[2 * GSHGEOM_SK_CHORDS_MAX];
  size_t m = 0;
  for (size_t i = 0; i < k; i++) {
    for (int e = 0; e < 2; e++) {
      GSHGEOM_SK_Vec p = chord[i][e];
      double v = gshgeom_sk_dot(gshgeom_sk_sub(p, y), perp);
      size_t j = 0;
      while (j < m && tau[j] != v) {
        j++;
      }
      if (j == m) {
        pt[m] = p;
        tau[m] = v;
        m++;
      }
    }
  }
  /* order the distinct ends along the line */
  for (size_t i = 1; i < m; i++) {
    GSHGEOM_SK_Vec vp = pt[i];
    double vt = tau[i];
    size_t j = i;
    while (j > 0 && tau[j - 1] > vt) {
      pt[j] = pt[j - 1];
      tau[j] = tau[j - 1];
      j--;
    }
    pt[j] = vp;
    tau[j] = vt;
  }
  for (size_t i = 0; i + 1 < m; i++) {
    int net = 0;
    for (size_t j = 0; j < k; j++) {
      double ta = gshgeom_sk_dot(gshgeom_sk_sub(chord[j][0], y), perp);
      double tb = gshgeom_sk_dot(gshgeom_sk_sub(chord[j][1], y), perp);
      if (ta < tb && ta <= tau[i] && tb >= tau[i + 1]) {
        net++;
      } else if (tb < ta && tb <= tau[i] && ta >= tau[i + 1]) {
        net--;
      }
    }
    for (int rep = 0; rep < (net < 0 ? -net : net); rep++) {
      GSHGEOM_Result r = net > 0 ? fn(ctx, pt[i], pt[i + 1]) : fn(ctx, pt[i + 1], pt[i]);
      if (r != GSHGEOM_OK) {
        return r;
      }
    }
  }
  return GSHGEOM_OK;
}

static GSHGEOM_Result net_edge(void * ctx, GSHGEOM_SK_Vec from, GSHGEOM_SK_Vec to) {
  return add_edge((Boundary *)ctx, from, to);
}

/** Close the boundary: the chords the clip line cut into it, as the edges that remain once the cancelled stretches are gone. */
static GSHGEOM_Result add_net_chords(Boundary * b) {
  if (b->have && !same(b->last, b->first)) {
    GSHGEOM_Result r = add_chord(b, b->last, b->first);
    if (r != GSHGEOM_OK) {
      return r;
    }
  }
  return gshgeom_sk_net_chords((const GSHGEOM_SK_Vec (*)[2])b->chord, b->nchords, b->cl.y, b->cl.t, net_edge, b);
}

/**
 * The region the boundary's edges enclose by nonzero winding, as positively wound convex pieces: the plane is cut
 * at every height where an edge ends or two edges cross, so that within a slab the edges do not cross and are
 * ordered by x; the winding between two neighbours is the running sum of the edges' directions, and each stretch of
 * nonzero winding is one trapezoid of its outer two edges.  @p dd bounds, in device pixels, how far the polygon of
 * edges is from the region and the region from it.
 */
static GSHGEOM_Result fill_boundary(Boundary * b, GSHGEOM_SK_Vec anchor, double dd) {
  GSHGEOM_SK_Ctx * c = b->c;
  size_t n = b->n;
  GSHGEOM_Result res = GSHGEOM_OK;
  double * ys = NULL;
  double * ytmp = NULL;
  Span * sp = NULL;
  Span * stmp = NULL;
  size_t ny = 0, ycap = 0;
  void * g = NULL;
  /* the heights */
  for (size_t i = 0; i < n && res == GSHGEOM_OK; i++) {
    for (int e = 0; e < 2 && res == GSHGEOM_OK; e++) {
      g = ys;
      res = gshgeom_grow(&c->alloc, &g, &ycap, ny + 1, sizeof *ys);
      if (res == GSHGEOM_OK) {
        ys = g;
        ys[ny++] = e ? b->e[i].q.y : b->e[i].p.y;
      }
    }
  }
  for (size_t i = 0; i < n && res == GSHGEOM_OK; i++) {
    GSHGEOM_SK_Vec d1 = gshgeom_sk_sub(b->e[i].q, b->e[i].p);
    if (d1.y == 0.0) {
      continue;
    }
    for (size_t j = i + 1; j < n && res == GSHGEOM_OK; j++) {
      GSHGEOM_SK_Vec d2 = gshgeom_sk_sub(b->e[j].q, b->e[j].p);
      if (d2.y == 0.0) {
        continue;
      }
      double den = gshgeom_sk_cross(d1, d2);
      if (den == 0.0) {
        continue;
      }
      GSHGEOM_SK_Vec w = gshgeom_sk_sub(b->e[j].p, b->e[i].p);
      double t = gshgeom_sk_cross(w, d2) / den, u = gshgeom_sk_cross(w, d1) / den;
      if (t > 0.0 && t < 1.0 && u > 0.0 && u < 1.0) {
        if (ny >= BOUNDARY_EVENTS) {
          res = GSHGEOM_ERR_LIMIT;
          break;
        }
        g = ys;
        res = gshgeom_grow(&c->alloc, &g, &ycap, ny + 1, sizeof *ys);
        if (res == GSHGEOM_OK) {
          ys = g;
          ys[ny++] = b->e[i].p.y + t * d1.y;
        }
      }
    }
  }
  if (res == GSHGEOM_OK) {
    ytmp = gcu_allocator_malloc(&c->alloc, (ny ? ny : 1) * sizeof *ytmp);
    sp = gcu_allocator_malloc(&c->alloc, (n ? n : 1) * sizeof *sp);
    stmp = gcu_allocator_malloc(&c->alloc, (n ? n : 1) * sizeof *stmp);
    if (!ytmp || !sp || !stmp) {
      res = GSHGEOM_ERR_OOM;
    }
  }
  if (res == GSHGEOM_OK) {
    sort_doubles(ys, ytmp, ny);
  }
  for (size_t si = 0; si + 1 < ny && res == GSHGEOM_OK; si++) {
    double y0 = ys[si], y1 = ys[si + 1];
    if (!(y1 > y0)) {
      continue;
    }
    double ym = y0 + (y1 - y0) * 0.5;
    size_t k = 0;
    for (size_t i = 0; i < n; i++) {
      const Edge * e = &b->e[i];
      double lo = e->p.y < e->q.y ? e->p.y : e->q.y, hi = e->p.y < e->q.y ? e->q.y : e->p.y;
      if (lo <= y0 && hi >= y1 && lo < hi) {
        sp[k].xm = edge_x(e, ym);
        sp[k].x0 = edge_x(e, y0);
        sp[k].x1 = edge_x(e, y1);
        sp[k].dir = e->q.y > e->p.y ? 1 : -1;
        k++;
      }
    }
    sort_spans(sp, stmp, k);
    int w = 0;
    size_t start = 0;
    for (size_t i = 0; i < k && res == GSHGEOM_OK; i++) {
      int before = w;
      w += sp[i].dir;
      if (before == 0 && w != 0) {
        start = i;
      } else if (before != 0 && w == 0) {
        double l0 = sp[start].x0, l1 = sp[start].x1, r0 = sp[i].x0, r1 = sp[i].x1;
        if (r0 == l0 && r1 == l1) {
          continue; /* no area */
        }
        GSHGEOM_SK_Vec quad[4] = {gshgeom_sk_v(l0, y0), gshgeom_sk_v(r0, y0), gshgeom_sk_v(r1, y1),
            gshgeom_sk_v(l1, y1)};
        res = gshgeom_sk_emit(c, anchor, quad, 4, dd, dd);
      }
    }
  }
  gcu_allocator_free(&c->alloc, ys);
  gcu_allocator_free(&c->alloc, ytmp);
  gcu_allocator_free(&c->alloc, sp);
  gcu_allocator_free(&c->alloc, stmp);
  return res;
}

GSHGEOM_Result gshgeom_sk_arcs_join(GSHGEOM_SK_Ctx * c, const GSHGEOM_SK_Seg * in, const GSHGEOM_SK_Seg * out,
    int outer_left, int reversal, int * miter_clip) {
  *miter_clip = 0;
  double h = c->h, L = c->st.miter_limit;
  GSHGEOM_SK_Vec a, b;
  (void)gshgeom_sk_end_dir(in, &a);
  (void)gshgeom_sk_start_dir(out, &b);
  GSHGEOM_SK_Vec v = in->p[in->n];
  GSHGEOM_SK_Vec ah = gshgeom_sk_unit(a), bh = gshgeom_sk_unit(b);
  GSHGEOM_SK_Vec oa = outer_left ? gshgeom_sk_left(ah) : gshgeom_sk_right(ah);
  GSHGEOM_SK_Vec ob = outer_left ? gshgeom_sk_left(bh) : gshgeom_sk_right(bh);
  GSHGEOM_SK_Vec A = gshgeom_sk_scale(h, oa), B = gshgeom_sk_scale(h, ob); /* offsets from V */
  GSHGEOM_SK_Vec anchor;
  GSHGEOM_Result r = gshgeom_sk_anchor(c, v, &anchor);
  if (r != GSHGEOM_OK) {
    return r;
  }
  EndCurv ka = end_curvature(in, 1, outer_left, h);
  EndCurv kb = end_curvature(out, 0, outer_left, h);
  if (ka.line && kb.line) {
    /* both curvatures zero: miter-clip's join (the caller's rule; at a reversal, its rectangle) */
    if (reversal) {
      GSHGEOM_SK_Vec rect[4] = {A, gshgeom_sk_add(A, gshgeom_sk_scale(L * h, ah)),
          gshgeom_sk_add(gshgeom_sk_scale(-1.0, A), gshgeom_sk_scale(L * h, ah)), gshgeom_sk_scale(-1.0, A)};
      return gshgeom_sk_emit(c, anchor, rect, 4, 0.0, 0.0);
    }
    *miter_clip = 1; /* both curvatures zero: the caller draws miter-clip's join */
    return GSHGEOM_OK;
  }
  if (ka.over || kb.over) {
    if (reversal) {
      GSHGEOM_Result rr = gshgeom_sk_sector(c, anchor, oa, ah, 0.0, 0.0);
      if (rr == GSHGEOM_OK) {
        rr = gshgeom_sk_sector(c, anchor, ah, gshgeom_sk_scale(-1.0, oa), 0.0, 0.0);
      }
      return rr;
    }
    return gshgeom_sk_sector(c, anchor, oa, ob, 0.0, 0.0);
  }
  /* the extensions, in build space relative to V */
  Ext ea, eb;
  memset(&ea, 0, sizeof ea);
  memset(&eb, 0, sizeof eb);
  ea.a = A;
  ea.dir = ah;
  ea.line = ka.line;
  eb.a = B;
  eb.dir = gshgeom_sk_scale(-1.0, bh);
  eb.line = kb.line;
  if (!ea.line) {
    ea.z = gshgeom_sk_scale(-1.0 / ka.kappa, oa);
    double rr = 1.0 / ka.kappa + h;
    ea.r = rr < 0 ? -rr : rr;
    ea.nu = rr > 0.0 ? gshgeom_sk_scale(-1.0, oa) : oa;
  }
  if (!eb.line) {
    eb.z = gshgeom_sk_scale(-1.0 / kb.kappa, ob);
    double rr = 1.0 / kb.kappa + h;
    eb.r = rr < 0 ? -rr : rr;
    eb.nu = rr > 0.0 ? gshgeom_sk_scale(-1.0, ob) : ob;
  }
  GSHGEOM_SK_Vec cand[2];
  int nc = 0;
  if (ea.line && eb.line) {
    nc = 0;
  } else if (ea.line) {
    nc = circle_line(&eb, &ea, cand);
  } else if (eb.line) {
    nc = circle_line(&ea, &eb, cand);
  } else {
    nc = circle_circle(&ea, &eb, cand);
  }
  if (nc == 2 && ea.line != eb.line) {
    /* a line's extension is a ray: forward from A along `a`, backward from B along `-b` (SG8); a meeting point behind
     * its start is not on it */
    const Ext * le = ea.line ? &ea : &eb;
    int ok0 = gshgeom_sk_dot(gshgeom_sk_sub(cand[0], le->a), le->dir) >= 0.0;
    int ok1 = gshgeom_sk_dot(gshgeom_sk_sub(cand[1], le->a), le->dir) >= 0.0;
    if (!ok0 && !ok1) {
      nc = 0;
    } else if (!ok0) {
      cand[0] = cand[1];
      nc = 1;
    } else if (!ok1) {
      nc = 1;
    }
  }
  if (nc == 0) {
    /* adjust the radii (the CR's 13.5.9) */
    if (!ea.line && !eb.line) {
      GSHGEOM_SK_Vec w0 = gshgeom_sk_sub(ea.z, eb.z);
      double dz = gshgeom_sk_norm(w0);
      double diff = ea.r > eb.r ? ea.r - eb.r : eb.r - ea.r;
      int enclosing = dz < diff;
      int a_big = ea.r > eb.r;
      double sa = 1.0, sb = 1.0, c0, k;
      if (enclosing) {
        sa = a_big ? -1.0 : 1.0;
        sb = -sa;
        c0 = diff;
        k = -2.0;
      } else {
        c0 = ea.r + eb.r;
        k = 2.0;
      }
      GSHGEOM_SK_Vec g = gshgeom_sk_sub(gshgeom_sk_scale(sa, ea.nu), gshgeom_sk_scale(sb, eb.nu));
      double dl = least_positive(gshgeom_sk_dot(g, g) - k * k, 2.0 * gshgeom_sk_dot(w0, g) - 2.0 * c0 * k,
          gshgeom_sk_dot(w0, w0) - c0 * c0);
      if (dl > 0.0) {
        ea.z = gshgeom_sk_add(ea.z, gshgeom_sk_scale(sa * dl, ea.nu));
        ea.r += sa * dl;
        eb.z = gshgeom_sk_add(eb.z, gshgeom_sk_scale(sb * dl, eb.nu));
        eb.r += sb * dl;
        if (enclosing && !a_big) {
          cand[0] = gshgeom_sk_add(eb.z, gshgeom_sk_scale(eb.r, gshgeom_sk_unit(gshgeom_sk_sub(ea.z, eb.z))));
        } else {
          cand[0] = gshgeom_sk_add(ea.z, gshgeom_sk_scale(ea.r, gshgeom_sk_unit(gshgeom_sk_sub(eb.z, ea.z))));
        }
        nc = 1;
      }
    } else if (ea.line != eb.line) {
      Ext * ce = ea.line ? &eb : &ea;
      const Ext * le = ea.line ? &ea : &eb;
      GSHGEOM_SK_Vec m = gshgeom_sk_left(le->dir);
      double d0 = gshgeom_sk_dot(gshgeom_sk_sub(ce->z, le->a), m), g = gshgeom_sk_dot(ce->nu, m);
      /* the circle grows until it touches the line, at a point on the ray: the least such growth */
      double best = -1.0;
      for (int sg = -1; sg <= 1; sg += 2) {
        double den = g - (double)sg;
        if (den != 0.0) {
          double dl = ((double)sg * ce->r - d0) / den;
          if (dl > 0.0 && (best < 0.0 || dl < best)) {
            GSHGEOM_SK_Vec zz = gshgeom_sk_add(ce->z, gshgeom_sk_scale(dl, ce->nu));
            GSHGEOM_SK_Vec foot = gshgeom_sk_sub(zz, gshgeom_sk_scale(gshgeom_sk_dot(gshgeom_sk_sub(zz, le->a), m), m));
            if (gshgeom_sk_dot(gshgeom_sk_sub(foot, le->a), le->dir) >= 0.0) {
              best = dl;
            }
          }
        }
      }
      if (best > 0.0) {
        ce->z = gshgeom_sk_add(ce->z, gshgeom_sk_scale(best, ce->nu));
        ce->r += best;
        cand[0] = gshgeom_sk_sub(ce->z, gshgeom_sk_scale(gshgeom_sk_dot(gshgeom_sk_sub(ce->z, le->a), m), m));
        nc = 1;
      }
    }
    if (nc == 0) {
      /* the CR's rectangle: parallel tangents that never meet */
      GSHGEOM_SK_Vec rect[4] = {A, gshgeom_sk_add(A, gshgeom_sk_scale(L * h, ah)),
          gshgeom_sk_add(gshgeom_sk_scale(-1.0, A), gshgeom_sk_scale(L * h, ah)), gshgeom_sk_scale(-1.0, A)};
      return gshgeom_sk_emit(c, anchor, rect, 4, 0.0, 0.0);
    }
  }
  int pick = 0;
  if (nc == 2) {
    /* The candidate closer to V is the one on V's side of the line of centres (two circles), or the one nearer V's
     * foot along the line (a circle and a line): the two are mirror images across that line, so the difference of
     * their squared distances is a multiple of V's signed distance from it.  They tie exactly when V is on it,
     * which is exactly a reversal (the offsets are then opposite, so both centres and V are collinear); the tie
     * goes to the one reached first along A's extension. */
    double side;
    if (!ea.line && !eb.line) {
      side = -gshgeom_sk_cross(gshgeom_sk_sub(eb.z, ea.z), gshgeom_sk_scale(-1.0, ea.z));
    } else {
      const Ext * ce = ea.line ? &eb : &ea;
      const Ext * le = ea.line ? &ea : &eb;
      GSHGEOM_SK_Vec m = gshgeom_sk_sub(le->a, ce->z);
      double bq = gshgeom_sk_dot(m, le->dir);
      side = bq - gshgeom_sk_dot(le->a, le->dir); /* V (the origin) minus the circle's foot, along the line */
    }
    if (reversal || side == 0.0) {
      if (along_key(&ea, cand[1]) < along_key(&ea, cand[0])) {
        pick = 1;
      }
    } else {
      pick = side > 0.0 ? 1 : 0;
    }
  }
  GSHGEOM_SK_Vec X = cand[pick];
  /* Gamma, the miter length along it, and the clip (offsets from V) */
  GSHGEOM_SK_Vec u = reversal ? ah : gshgeom_sk_unit(gshgeom_sk_add(oa, ob));
  Clip cl;
  memset(&cl, 0, sizeof cl);
  double dn = gshgeom_sk_norm(X), Lh = L * h;
  if (dn > 0.0) {
    double cr = gshgeom_sk_cross(u, X), dt = gshgeom_sk_dot(u, X);
    if (cr == 0.0) {
      if (dt <= 0.0 || dn > Lh) {
        cl.on = 1;
        cl.y = gshgeom_sk_scale(Lh, u);
        cl.t = u;
      }
    } else {
      GSHGEOM_SK_Vec w = cr > 0.0 ? gshgeom_sk_right(u) : gshgeom_sk_left(u);
      if (gshgeom_sk_dot(w, X) < 0.0) {
        w = gshgeom_sk_scale(-1.0, w);
      }
      double rho = dn * dn / (2.0 * gshgeom_sk_dot(X, w));
      double beta = gmath_atan2(cr < 0 ? -cr : cr, dt);
      if (2.0 * beta * rho > Lh) {
        double psi = Lh / rho;
        double cp = gmath_cos(psi), sp = gmath_sin(psi);
        GSHGEOM_SK_Vec O = gshgeom_sk_scale(rho, w);
        cl.on = 1;
        cl.y = gshgeom_sk_add(O, gshgeom_sk_scale(rho, gshgeom_sk_add(gshgeom_sk_scale(-cp, w), gshgeom_sk_scale(sp, u))));
        cl.t = gshgeom_sk_unit(gshgeom_sk_add(gshgeom_sk_scale(sp, w), gshgeom_sk_scale(cp, u)));
      }
    }
  }
  /* The region (SG8): the bevel, together with what the closed boundary V, A, A's extension to X, B's extension back
   * to B, V encloses by nonzero winding, the whole cut by the clip's half plane.  The clip is applied to the
   * boundary itself (the parts of it in the half plane, closed along the clip line), which cuts the winding region
   * exactly: a curve and its clipping wind the same about every point of the half plane. */
  Boundary bd;
  memset(&bd, 0, sizeof bd);
  bd.c = c;
  bd.cl = cl;
  /* the largest coordinate the construction touches; anything not finite, or past 2^150 (a circle of curvature so
   * slight that its square or its cube is out of a double's range), is refused: a stated cap */
  double s0 = 0.0;
  int bad = 0;
  {
    double vals[14] = {A.x, A.y, B.x, B.y, X.x, X.y, ea.z.x, ea.z.y, eb.z.x, eb.z.y, ea.r, eb.r, cl.y.x, cl.y.y};
    for (int i = 0; i < 14; i++) {
      double m = fabs_d(vals[i]);
      if (!(m <= 0x1p150)) {
        bad = 1;
      } else if (m > s0) {
        s0 = m;
      }
    }
    if (!ea.line) {
      s0 = s0 + ea.r;
    }
    if (!eb.line) {
      s0 = s0 + eb.r;
    }
  }
  if (bad) {
    return GSHGEOM_ERR_LIMIT;
  }
  /* what the double arithmetic of the construction may move a point by: 2^-30 of the largest coordinate, in pixels */
  double eps = (s0 * 0x1p-30) * c->norm;
  double room = gshgeom_sk_room(c, 0.0) / (1.0 + GSHGEOM_SK_OCT);
  bd.allowed = room * 0.999 - eps;
  if (!(bd.allowed > 0.0) || !(s0 * c->norm < GSHGEOM_SK_MAG_LIMIT)) {
    return GSHGEOM_ERR_LIMIT;
  }
  GSHGEOM_SK_Vec zero = gshgeom_sk_v(0.0, 0.0);
  /* the bevel, cut exactly */
  if (!reversal) {
    GSHGEOM_SK_Vec tri[3] = {zero, A, B}, cut[8];
    size_t m = clip_poly(tri, 3, cut, &cl);
    if (m >= 3) {
      r = gshgeom_sk_emit(c, anchor, cut, m, eps, eps);
    }
  }
  if (r == GSHGEOM_OK) {
    r = add_line(&bd, zero, A);
  }
  if (r == GSHGEOM_OK) {
    r = ea.line ? add_line(&bd, A, X) : add_arc(&bd, ea.z, ea.r, A, X, ext_sense(&ea));
  }
  if (r == GSHGEOM_OK) {
    r = eb.line ? add_line(&bd, X, B) : add_arc(&bd, eb.z, eb.r, X, B, -ext_sense(&eb));
  }
  if (r == GSHGEOM_OK) {
    r = add_line(&bd, B, zero);
  }
  if (r == GSHGEOM_OK) {
    r = add_net_chords(&bd);
  }
  if (r == GSHGEOM_OK) {
    double dd = bd.worst + eps;
    if (!gshgeom_sk_fits(c, dd, dd)) {
      r = GSHGEOM_ERR_LIMIT;
    } else {
      r = fill_boundary(&bd, anchor, dd);
    }
  }
  gcu_allocator_free(&c->alloc, bd.e);
  gcu_allocator_free(&c->alloc, bd.pts);
  return r;
}
