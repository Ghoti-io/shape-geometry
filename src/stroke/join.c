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
 * Joins (SG8), the miter limit decided exactly (SG9), caps (SG7) and degenerate
 * shapes (SG11), as pieces (emit.c).
 *
 * Every decision is an exact predicate on the build-space doubles: the turn's
 * sign `sigma` (orient2d), the sign of `a . b` (eight products), and SG9's
 * comparison (below).  Every shape here is a polygon of SG8's own corners (the
 * bevel, the miter quadrilateral, the pentagon, the caps' rectangles, the
 * degenerate square) or a circular sector (arc.c), so a polygon piece has
 * `d_in = d_out = 0` beyond the arithmetic of emit.c (a), except the miter tip,
 * whose computation from nearly parallel edges carries the bound below.
 *
 * ## SG9
 *
 * The miter is within the limit when `(L^2 - 2) |a||b| >= -L^2 (a . b)`.  With
 * `P = L^2 - 2` and `q = a . b` (both signs exact): `P >= 0, q >= 0` is within;
 * `P <= 0, q < 0` is not; otherwise both sides have one sign and the test is the
 * sign of `D = P^2 |a|^2 |b|^2 - L^4 q^2`, written with
 * `|a|^2 |b|^2 = q^2 + c^2` (`c = a x b`) as `D = 4 (1 - L^2) q^2 + P^2 c^2`:
 * 236 products of at most eight of `L` and the raw coordinates of `Q`, `V`, `S`
 * (`a = V - Q`, `b = S - V` are never formed).  A filter takes the sign from the
 * products in double when it exceeds `2^-40` of their absolute sum (each product
 * of at most eight factors is within `7u` of itself and the sum of 236 within
 * `243 u < 2^-44` of their absolute sum, while nothing underflows: the absolute
 * sum is at least `2^-900`); otherwise src/exact's wide sum decides.
 */

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/shape-geometry/orient.h>

#include "../exact/exact.h"
#include "stroke_internal.h"

#include <string.h>

int gshgeom_sk_start_dir(const GSHGEOM_SK_Seg * s, GSHGEOM_SK_Vec * dir) {
  for (int k = 1; k <= s->n; k++) {
    if (s->p[k].x != s->p[0].x || s->p[k].y != s->p[0].y) {
      *dir = gshgeom_sk_sub(s->p[k], s->p[0]);
      return k;
    }
  }
  return -1;
}

int gshgeom_sk_end_dir(const GSHGEOM_SK_Seg * s, GSHGEOM_SK_Vec * dir) {
  for (int k = s->n - 1; k >= 0; k--) {
    if (s->p[k].x != s->p[s->n].x || s->p[k].y != s->p[s->n].y) {
      *dir = gshgeom_sk_sub(s->p[s->n], s->p[k]);
      return k;
    }
  }
  return -1;
}

/** The sign of a . b, a = V - Q, b = S - V, exactly: eight products of raw coordinates. */
static int dot_sign(GSHGEOM_SK_Vec q, GSHGEOM_SK_Vec v, GSHGEOM_SK_Vec s) {
  GSHGEOM_ExactTerm t[8];
  size_t n = 0;
  /* (Vx - Qx)(Sx - Vx) = VxSx - Vx^2 - QxSx + QxVx, and the same in y */
  gshgeom_exact_put(t, &n, 0, 2, v.x, s.x, 0.0);
  gshgeom_exact_put(t, &n, 1, 2, v.x, v.x, 0.0);
  gshgeom_exact_put(t, &n, 1, 2, q.x, s.x, 0.0);
  gshgeom_exact_put(t, &n, 0, 2, q.x, v.x, 0.0);
  gshgeom_exact_put(t, &n, 0, 2, v.y, s.y, 0.0);
  gshgeom_exact_put(t, &n, 1, 2, v.y, v.y, 0.0);
  gshgeom_exact_put(t, &n, 1, 2, q.y, s.y, 0.0);
  gshgeom_exact_put(t, &n, 0, 2, q.y, v.y, 0.0);
  return gshgeom_exact_sign(t, n);
}

/** The terms of q = a . b (8) and c = a x b (6) in raw coordinates: a factor pair and a sign each. */
typedef struct RawTerm {
  double f0, f1;
  int neg;
} RawTerm;

static void raw_terms(GSHGEOM_SK_Vec q, GSHGEOM_SK_Vec v, GSHGEOM_SK_Vec s, RawTerm qd[8], RawTerm cr[6]) {
  RawTerm d[8] = {{v.x, s.x, 0}, {v.x, v.x, 1}, {q.x, s.x, 1}, {q.x, v.x, 0}, {v.y, s.y, 0}, {v.y, v.y, 1},
      {q.y, s.y, 1}, {q.y, v.y, 0}};
  /* (Vx - Qx)(Sy - Vy) - (Vy - Qy)(Sx - Vx) = VxSy - QxSy + QxVy - VySx + QySx - QyVx */
  RawTerm c[6] = {{v.x, s.y, 0}, {q.x, s.y, 1}, {q.x, v.y, 0}, {v.y, s.x, 1}, {q.y, s.x, 0}, {q.y, v.x, 1}};
  memcpy(qd, d, sizeof d);
  memcpy(cr, c, sizeof c);
}

/** One product of D's expansion, handed to the filter and to the exact sum. */
typedef struct Prod {
  double f[8];
  int n;
  int neg;
} Prod;

/** Enumerate D = 4 (1 - L^2) q^2 + (L^2 - 2)^2 c^2 = 4q^2 - 4L^2q^2 + L^4c^2 - 4L^2c^2 + 4c^2. */
static size_t d_products(GSHGEOM_SK_Vec q, GSHGEOM_SK_Vec v, GSHGEOM_SK_Vec s, double L, Prod * out) {
  RawTerm qd[8], cr[6];
  raw_terms(q, v, s, qd, cr);
  size_t n = 0;
  for (int i = 0; i < 8; i++) {
    for (int j = 0; j < 8; j++) {
      int sg = qd[i].neg ^ qd[j].neg;
      Prod a = {{4.0, qd[i].f0, qd[i].f1, qd[j].f0, qd[j].f1, 0, 0, 0}, 5, sg};
      out[n++] = a;
      Prod b = {{4.0, L, L, qd[i].f0, qd[i].f1, qd[j].f0, qd[j].f1, 0}, 7, !sg};
      out[n++] = b;
    }
  }
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      int sg = cr[i].neg ^ cr[j].neg;
      Prod a = {{L, L, L, L, cr[i].f0, cr[i].f1, cr[j].f0, cr[j].f1}, 8, sg};
      out[n++] = a;
      Prod b = {{4.0, L, L, cr[i].f0, cr[i].f1, cr[j].f0, cr[j].f1, 0}, 7, !sg};
      out[n++] = b;
      Prod c = {{4.0, cr[i].f0, cr[i].f1, cr[j].f0, cr[j].f1, 0, 0, 0}, 5, sg};
      out[n++] = c;
    }
  }
  return n; /* 128 + 108 = 236 */
}

static int d_sign(GSHGEOM_SK_Vec q, GSHGEOM_SK_Vec v, GSHGEOM_SK_Vec s, double L) {
  Prod p[236];
  size_t n = d_products(q, v, s, L, p);
  double sum = 0.0, abs = 0.0;
  for (size_t i = 0; i < n; i++) {
    double t = 1.0;
    for (int k = 0; k < p[i].n; k++) {
      t = t * p[i].f[k];
    }
    sum = p[i].neg ? sum - t : sum + t;
    abs = abs + (t < 0 ? -t : t);
  }
  if (abs >= 0x1p-900 && abs - abs == 0.0) {
    double bound = abs * 0x1p-40;
    if (sum > bound) {
      return 1;
    }
    if (sum < -bound) {
      return -1;
    }
  }
  GSHGEOM_ExactSum es;
  gshgeom_exact_sum_init(&es);
  for (size_t i = 0; i < n; i++) {
    gshgeom_exact_sum_add(&es, p[i].neg, p[i].n, p[i].f);
  }
  return gshgeom_exact_sum_sign(&es);
}

int gshgeom_sk_miter_within(GSHGEOM_SK_Vec q, GSHGEOM_SK_Vec v, GSHGEOM_SK_Vec s, double L) {
  GSHGEOM_ExactTerm t[2];
  size_t n = 0;
  gshgeom_exact_put(t, &n, 0, 2, L, L, 0.0);
  gshgeom_exact_put(t, &n, 1, 1, 2.0, 0.0, 0.0);
  int P = gshgeom_exact_sign(t, n);
  int qs = dot_sign(q, v, s);
  if (P >= 0 && qs >= 0) {
    return 1;
  }
  if (P <= 0 && qs < 0) {
    return 0;
  }
  int D = d_sign(q, v, s, L);
  if (P < 0) {
    return D <= 0; /* q >= 0: |P||a||b| <= L^2 q */
  }
  return D >= 0; /* P > 0, q < 0: P|a||b| >= L^2 |q| */
}

/** The miter quadrilateral's tip offset from V, and a bound of its arithmetic error in device pixels. */
static GSHGEOM_SK_Vec miter_tip(const GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec oa, GSHGEOM_SK_Vec ob, double * err) {
  double k = c->h / (1.0 + gshgeom_sk_dot(oa, ob));
  GSHGEOM_SK_Vec t = gshgeom_sk_scale(k, gshgeom_sk_add(oa, ob));
  /* 1 + oa . ob carries an absolute error of a few u (the unit vectors and the dot), so a relative one of a few u
   * over itself; the tip's length is h / cos(phi/2) = |T| and 1 + cos(phi) = 2 (h / |T|)^2.  So the tip moves by at
   * most 8 u (|T| / h)^2 |T| in build space: 2^-49 (|T|/h)^2 |T|, times the norm in device pixels. */
  double len = gshgeom_sk_norm(t);
  double ratio = len / c->h;
  *err = 0x1p-49 * ratio * ratio * len * c->norm * (1.0 + 0x1p-40);
  return t;
}

/** Emit a polygon of offsets from V with d_in = d_out = e. */
static GSHGEOM_Result poly(GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec anchor, const GSHGEOM_SK_Vec * p, size_t n, double e) {
  if (!gshgeom_sk_fits(c, e, e)) {
    return GSHGEOM_ERR_LIMIT; /* a miter tip the arithmetic cannot place within the tolerance */
  }
  return gshgeom_sk_emit(c, anchor, p, n, e, e);
}

/** The miter-clip join, the outer side given (also the arcs join's fallback). */
static GSHGEOM_Result miter_clip(GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec anchor, GSHGEOM_SK_Vec ah, GSHGEOM_SK_Vec bh,
    GSHGEOM_SK_Vec oa, GSHGEOM_SK_Vec ob, int within) {
  double h = c->h, L = c->st.miter_limit;
  GSHGEOM_SK_Vec A = gshgeom_sk_scale(h, oa), B = gshgeom_sk_scale(h, ob);
  if (within) {
    double e = 0.0;
    GSHGEOM_SK_Vec T = miter_tip(c, oa, ob, &e);
    GSHGEOM_SK_Vec q[4] = {gshgeom_sk_v(0.0, 0.0), A, T, B};
    return poly(c, anchor, q, 4, e);
  }
  GSHGEOM_SK_Vec u = gshgeom_sk_unit(gshgeom_sk_add(oa, ob));
  double la = (L * h - gshgeom_sk_dot(A, u)) / gshgeom_sk_dot(ah, u);
  double lb = (L * h - gshgeom_sk_dot(B, u)) / -gshgeom_sk_dot(bh, u);
  GSHGEOM_SK_Vec q[5] = {gshgeom_sk_v(0.0, 0.0), A, gshgeom_sk_add(A, gshgeom_sk_scale(la, ah)),
      gshgeom_sk_sub(B, gshgeom_sk_scale(lb, bh)), B};
  return poly(c, anchor, q, 5, 0.0);
}

/** The half disc {E + v : |v| <= h, v . d >= 0} as one convex piece (two quarter arcs). */
static GSHGEOM_Result half_disc(GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec anchor, GSHGEOM_SK_Vec m, GSHGEOM_SK_Vec d) {
  if (c->cv) {
    return gshgeom_skc_half_disc(c, anchor, m, d); /* two quarter arcs, exact (SG37) */
  }
  double allowed = gshgeom_sk_room(c, 0.0);
  size_t n = 1;
  double d1 = 0.0, d2 = 0.0;
  GSHGEOM_Result r = gshgeom_sk_arc(c, c->h, m, d, allowed, &n, &d1);
  if (r == GSHGEOM_OK) {
    n--; /* the second arc starts where the first ended */
    r = gshgeom_sk_arc(c, c->h, d, gshgeom_sk_scale(-1.0, m), allowed, &n, &d2);
  }
  if (r == GSHGEOM_OK) {
    c->arcbuf[0] = gshgeom_sk_v(0.0, 0.0); /* the centre; slot 0 was kept for it */
    r = gshgeom_sk_emit(c, anchor, c->arcbuf, n, 0.0, d1 > d2 ? d1 : d2);
  }
  return r;
}

int gshgeom_sk_vertex_class(GSHGEOM_SK_Vec q, GSHGEOM_SK_Vec v, GSHGEOM_SK_Vec s) {
  GSHGEOM_Point pq = {q.x, q.y, 0.0}, pv = {v.x, v.y, 0.0}, ps = {s.x, s.y, 0.0};
  int sigma = 0;
  if (gshgeom_orient2d(pq, pv, ps, &sigma) != GSHGEOM_OK || sigma != 0) {
    return 0;
  }
  return dot_sign(q, v, s) > 0 ? 1 : 2;
}

GSHGEOM_Result gshgeom_sk_join(GSHGEOM_SK_Ctx * c, const GSHGEOM_SK_Seg * in, const GSHGEOM_SK_Seg * out) {
  GSHGEOM_SK_Vec a, b;
  int qi = gshgeom_sk_end_dir(in, &a);
  int si = gshgeom_sk_start_dir(out, &b);
  if (qi < 0 || si < 0) {
    return GSHGEOM_ERR_INTERNAL; /* segments without extent were dropped */
  }
  GSHGEOM_SK_Vec v = in->p[in->n], q = in->p[qi], s = out->p[si];
  GSHGEOM_Point pq = {q.x, q.y, 0.0}, pv = {v.x, v.y, 0.0}, ps = {s.x, s.y, 0.0};
  int sigma = 0;
  if (gshgeom_orient2d(pq, pv, ps, &sigma) != GSHGEOM_OK) {
    return GSHGEOM_ERR_INTERNAL;
  }
  int ds = dot_sign(q, v, s);
  if (out->vtx_valid) {
    /* the class decided before the map (see gshgeom_sk_vertex_class): the mapped doubles can show a turn of a few ulps
     * at a vertex that is smooth, or miss one that is not */
    if (out->vtx_class == 1) {
      return GSHGEOM_OK;
    }
    if (out->vtx_class == 2) {
      sigma = 0;
      ds = -1;
    }
  }
  if (sigma == 0 && ds > 0) {
    return GSHGEOM_OK; /* a smooth vertex: no join */
  }
  GSHGEOM_SK_Vec anchor;
  GSHGEOM_Result r = gshgeom_sk_anchor(c, v, &anchor);
  if (r != GSHGEOM_OK) {
    return r;
  }
  double h = c->h, L = c->st.miter_limit;
  GSHGEOM_SK_Vec ah = gshgeom_sk_unit(a), bh = gshgeom_sk_unit(b);
  if (sigma == 0) {
    /* a reversal: each join's limit */
    GSHGEOM_SK_Vec n = gshgeom_sk_left(ah);
    switch (c->st.join) {
    case GSHGEOM_STROKE_JOIN_ROUND:
      return half_disc(c, anchor, n, ah);
    case GSHGEOM_STROKE_JOIN_MITER_CLIP: {
      GSHGEOM_SK_Vec hn = gshgeom_sk_scale(h, n), f = gshgeom_sk_scale(L * h, ah);
      GSHGEOM_SK_Vec rect[4] = {hn, gshgeom_sk_add(hn, f), gshgeom_sk_sub(f, hn), gshgeom_sk_scale(-1.0, hn)};
      return poly(c, anchor, rect, 4, 0.0);
    }
    case GSHGEOM_STROKE_JOIN_ARCS: {
      int mc = 0;
      r = gshgeom_sk_arcs_join(c, in, out, 1, 1, &mc);
      if (r == GSHGEOM_OK) {
        r = gshgeom_sk_arcs_join(c, in, out, 0, 1, &mc);
      }
      return r;
    }
    default:
      return GSHGEOM_OK; /* bevel and miter add nothing of area */
    }
  }
  int outer_left = sigma > 0;
  GSHGEOM_SK_Vec oa = outer_left ? gshgeom_sk_left(ah) : gshgeom_sk_right(ah);
  GSHGEOM_SK_Vec ob = outer_left ? gshgeom_sk_left(bh) : gshgeom_sk_right(bh);
  GSHGEOM_SK_Vec A = gshgeom_sk_scale(h, oa), B = gshgeom_sk_scale(h, ob);
  GSHGEOM_SK_Vec bevel[3] = {gshgeom_sk_v(0.0, 0.0), A, B};
  switch (c->st.join) {
  case GSHGEOM_STROKE_JOIN_BEVEL:
    return poly(c, anchor, bevel, 3, 0.0);
  case GSHGEOM_STROKE_JOIN_MITER:
    if (gshgeom_sk_miter_within(q, v, s, L)) {
      return miter_clip(c, anchor, ah, bh, oa, ob, 1);
    }
    return poly(c, anchor, bevel, 3, 0.0);
  case GSHGEOM_STROKE_JOIN_MITER_CLIP:
    return miter_clip(c, anchor, ah, bh, oa, ob, gshgeom_sk_miter_within(q, v, s, L));
  case GSHGEOM_STROKE_JOIN_ROUND:
    return gshgeom_sk_sector(c, anchor, oa, ob, 0.0, 0.0);
  case GSHGEOM_STROKE_JOIN_ARCS: {
    int mc = 0;
    r = gshgeom_sk_arcs_join(c, in, out, outer_left, 0, &mc);
    if (r == GSHGEOM_OK && mc) {
      r = miter_clip(c, anchor, ah, bh, oa, ob, gshgeom_sk_miter_within(q, v, s, L));
    }
    return r;
  }
  default:
    return GSHGEOM_ERR_INTERNAL;
  }
}

GSHGEOM_Result gshgeom_sk_cap(GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec e, GSHGEOM_SK_Vec d) {
  if (c->st.cap == GSHGEOM_STROKE_CAP_BUTT) {
    return GSHGEOM_OK;
  }
  GSHGEOM_SK_Vec anchor;
  GSHGEOM_Result r = gshgeom_sk_anchor(c, e, &anchor);
  if (r != GSHGEOM_OK) {
    return r;
  }
  GSHGEOM_SK_Vec dh = gshgeom_sk_unit(d), m = gshgeom_sk_left(dh);
  double h = c->h;
  if (c->st.cap == GSHGEOM_STROKE_CAP_ROUND) {
    return half_disc(c, anchor, m, dh);
  }
  GSHGEOM_SK_Vec hm = gshgeom_sk_scale(h, m), hd = gshgeom_sk_scale(h, dh);
  GSHGEOM_SK_Vec rect[4] = {hm, gshgeom_sk_add(hm, hd), gshgeom_sk_sub(hd, hm), gshgeom_sk_scale(-1.0, hm)};
  return poly(c, anchor, rect, 4, 0.0);
}

GSHGEOM_Result gshgeom_sk_degenerate(GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec p, GSHGEOM_SK_Vec d) {
  if (c->st.cap == GSHGEOM_STROKE_CAP_BUTT) {
    return GSHGEOM_OK;
  }
  GSHGEOM_SK_Vec anchor;
  GSHGEOM_Result r = gshgeom_sk_anchor(c, p, &anchor);
  if (r != GSHGEOM_OK) {
    return r;
  }
  if (c->st.cap == GSHGEOM_STROKE_CAP_ROUND) {
    return gshgeom_sk_disc(c, anchor, 0.0, 0.0);
  }
  GSHGEOM_SK_Vec dh = gshgeom_sk_unit(d), m = gshgeom_sk_left(dh);
  GSHGEOM_SK_Vec hm = gshgeom_sk_scale(c->h, m), hd = gshgeom_sk_scale(c->h, dh);
  GSHGEOM_SK_Vec sq[4] = {gshgeom_sk_add(hm, hd), gshgeom_sk_sub(hd, hm), gshgeom_sk_scale(-1.0, gshgeom_sk_add(hm, hd)),
      gshgeom_sk_sub(hm, hd)};
  return poly(c, anchor, sq, 4, 0.0);
}
