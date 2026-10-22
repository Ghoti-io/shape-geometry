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
 * The dash positions (stroke.h, SG21 and SG29) on a raw dash array.
 *
 * The positions are SVG 2's ("dash positions", 13.5.7 of the 2018 Candidate
 * Recommendation), read literally (stroke.h, SG21).  The text
 * keeps a running position, adding each dash length to it, so the error of a
 * boundary grows with the number of dashes before it.  Here a boundary is
 * `(q P + S) - o`: the period count `q` (an integer, exact in a double) times
 * the period `P`, plus a prefix sum `S` formed the same way in every period, minus
 * the phase `o`; its error is that of three roundings of the magnitudes involved
 * and does not depend on the dashes before it.  The phase is reduced to `[0, P)` by
 * an exact long division: the largest `P 2^k` not above the remainder is
 * subtracted, then `P 2^(k-1)`, and so on; each subtraction is of a number not
 * above the remainder and below twice it, so it is exact (Sterbenz), and the
 * result is the remainder of the exact division, not an approximation of it.
 */

#include <ghoti.io/shape-geometry/macros.h>

#include "stroke_internal.h"

#include <ghoti.io/shape-geometry/along.h>

/** The remainder of x >= 0 modulo s > 0, exactly (+ - * only; no libm). */
static double mod_pos(double x, double s) {
  if (x < s) {
    return x;
  }
  double t = s;
  int k = 0;
  while (t * 2.0 <= x) {
    t *= 2.0;
    k++;
  }
  double r = x;
  for (;;) {
    if (r >= t) {
      r -= t;
    }
    if (k == 0) {
      break;
    }
    t *= 0.5;
    k--;
  }
  return r;
}

/** The phase reduced to [0, s): the floored modulo (SG21). */
static double reduce(double x, double s) {
  double r;
  if (x >= 0.0) {
    r = mod_pos(x, s);
  } else {
    double q = mod_pos(-x, s);
    r = q == 0.0 ? 0.0 : s - q;
  }
  return r >= s ? 0.0 : r;
}

/** The canonical entry i of a raw array (SG3: an odd array is repeated once). */
static double entry(const double * dashes, size_t count, size_t i) {
  return dashes[i % count];
}

/** SG3 on a raw array: the canonical count, 0 for a solid array. */
static size_t canonical_count(const double * dashes, size_t count) {
  int any = 0;
  for (size_t i = 0; i < count; i++) {
    if (dashes[i] < 0.0) {
      return 0;
    }
    any = any || dashes[i] != 0.0;
  }
  if (!any) {
    return 0;
  }
  return (count % 2u) ? count * 2u : count;
}

static int finite(double x) {
  return x - x == 0.0;
}

static int in_band(double x) {
  double m = x < 0.0 ? -x : x;
  return m == 0.0 || (m >= GSHGEOM_STROKE_VALUE_MIN && m <= GSHGEOM_STROKE_VALUE_MAX);
}

double gshgeom_sk_dash_period(const double * dashes, size_t count, size_t n) {
  double period = 0.0;
  for (size_t i = 0; i < n; i++) {
    period += entry(dashes, count, i);
  }
  return period;
}

double gshgeom_sk_dash_phase(double offset, double period) {
  return reduce(offset, period);
}

GSHGEOM_Result gshgeom_sk_dash_positions(const double * dashes, size_t count, size_t n, double period, double o,
    double length, size_t * steps, size_t max_steps, GSHGEOM_Dash_Fn fn, void * ctx) {
  /* the least index whose prefix sum reaches the phase */
  double s_ex = 0.0, s_in = 0.0;
  size_t index = n - 1;
  for (size_t i = 0; i < n; i++) {
    s_in = s_ex + entry(dashes, count, i);
    if (s_in >= o) {
      index = i;
      break;
    }
    s_ex = s_in;
  }
  double remaining = s_in - o;
  /* the first piece: a dash when the index is even (SG21, the text's rule: a phase exactly at the end of a dash
   * leaves a dash of no length at the start) */
  if ((index % 2u) == 0) {
    double e = remaining < length ? remaining : length;
    GSHGEOM_Result r = fn(ctx, 0.0, e);
    if (r != GSHGEOM_OK) {
      return r;
    }
  }
  double q = 0.0;
  size_t i = index + 1;
  s_ex = s_in;
  for (;;) {
    if (i == n) {
      i = 0;
      q += 1.0;
      s_ex = 0.0;
    }
    double d = entry(dashes, count, i);
    double s_end = s_ex + d;
    if ((i % 2u) == 0) {
      double start = (q * period + s_ex) - o;
      if (!(start < length)) {
        break;
      }
      double end = (q * period + s_end) - o;
      if (end > length) {
        end = length;
      }
      GSHGEOM_Result r = fn(ctx, start, end);
      if (r != GSHGEOM_OK) {
        return r;
      }
    }
    s_ex = s_end;
    i++;
    if (++*steps > max_steps) {
      return GSHGEOM_ERR_LIMIT;
    }
  }
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_stroke_dash_each(const double * dashes, size_t count, double offset, double length,
    size_t max_steps, GSHGEOM_Dash_Fn fn, void * ctx) {
  if (!fn || (count && !dashes) || !finite(offset) || !finite(length) || length < 0.0) {
    return GSHGEOM_ERR_INVALID;
  }
  for (size_t i = 0; i < count; i++) {
    if (!finite(dashes[i])) {
      return GSHGEOM_ERR_INVALID;
    }
  }
  if (!in_band(offset)) {
    return GSHGEOM_ERR_LIMIT;
  }
  for (size_t i = 0; i < count; i++) {
    if (!in_band(dashes[i])) {
      return GSHGEOM_ERR_LIMIT;
    }
  }
  if (count > GSHGEOM_STROKE_DASHES_MAX * 2u) {
    return GSHGEOM_ERR_LIMIT;
  }
  size_t n = canonical_count(dashes, count);
  if (n > GSHGEOM_STROKE_DASHES_MAX) {
    return GSHGEOM_ERR_LIMIT;
  }
  if (n == 0) {
    return fn(ctx, 0.0, length); /* solid: the single pair */
  }
  double period = gshgeom_sk_dash_period(dashes, count, n);
  size_t steps = 0;
  return gshgeom_sk_dash_positions(dashes, count, n, period, reduce(offset, period), length, &steps, max_steps, fn,
      ctx);
}

/* ---------------------------------------------------------------------- */
/* The expansion of a dashed subpath (SG19 to SG29)                        */
/* ---------------------------------------------------------------------- */

typedef GSHGEOM_SK_Seg Seg;
typedef GSHGEOM_SK_Vec Vec;
typedef GSHGEOM_SK_Sub Sub;

static double fabs_(double x) {
  return x < 0.0 ? -x : x;
}

/** De Casteljau at t on m values. */
static double dc(const double * v, int m, double t) {
  double a[4];
  for (int i = 0; i < m; i++) {
    a[i] = v[i];
  }
  double u = 1.0 - t;
  for (int lev = 1; lev < m; lev++) {
    for (int i = 0; i + lev < m; i++) {
      a[i] = u * a[i] + t * a[i + 1];
    }
  }
  return a[0];
}

/**
 * Split a segment at t in (0, 1) into the piece [0, t] and the piece [t, 1], in homogeneous coordinates for a rational
 * kind (the weights stay positive).  The outer end points are the segment's own, bit for bit, so a piece that keeps an
 * end of the segment keeps that vertex exactly.
 */
static void seg_split(const Seg * s, double t, Seg * l, Seg * r) {
  int n = s->n;
  double u = 1.0 - t;
  double X[4], Y[4], W[4];
  for (int i = 0; i <= n; i++) {
    W[i] = s->rational ? s->w[i] : 1.0;
    X[i] = s->rational ? s->w[i] * s->p[i].x : s->p[i].x;
    Y[i] = s->rational ? s->w[i] * s->p[i].y : s->p[i].y;
  }
  *l = *s;
  *r = *s;
  for (int lev = 1; lev <= n; lev++) {
    for (int i = 0; i + lev <= n; i++) {
      X[i] = u * X[i] + t * X[i + 1];
      Y[i] = u * Y[i] + t * Y[i + 1];
      if (s->rational) {
        W[i] = u * W[i] + t * W[i + 1];
      }
    }
    double w0 = s->rational ? W[0] : 1.0;
    double wr = s->rational ? W[n - lev] : 1.0;
    {
      l->p[lev] = s->rational ? gshgeom_sk_v(X[0] / w0, Y[0] / w0) : gshgeom_sk_v(X[0], Y[0]);
      l->w[lev] = w0;
      r->p[n - lev] = s->rational ? gshgeom_sk_v(X[n - lev] / wr, Y[n - lev] / wr) : gshgeom_sk_v(X[n - lev], Y[n - lev]);
      r->w[n - lev] = wr;
    }
  }
  l->p[0] = s->p[0];
  l->w[0] = s->w[0];
  r->p[n] = s->p[n];
  r->w[n] = s->w[n];
}

/** Whether a piece has extent (SG5) and no nonzero difference of neighbours below the arithmetic's floor. */
static int piece_ok(const Seg * s) {
  int any = 0;
  for (int k = 1; k <= s->n; k++) {
    if (s->p[k].x != s->p[0].x || s->p[k].y != s->p[0].y) {
      any = 1;
    }
    double dx = fabs_(s->p[k].x - s->p[k - 1].x), dy = fabs_(s->p[k].y - s->p[k - 1].y);
    double dm = dx > dy ? dx : dy;
    if (dm != 0.0 && dm < GSHGEOM_SK_BUILD_MIN) {
      return 0;
    }
  }
  return any;
}

/** The sub-curve of a segment between the parameters a < b; 0 when it has no extent after the cut. */
static int piece(const Seg * s, double a, double b, Seg * out) {
  if (a <= 0.0 && b >= 1.0) {
    *out = *s;
    return piece_ok(s);
  }
  Seg cur = *s;
  if (a > 0.0) {
    Seg l, r;
    seg_split(&cur, a, &l, &r);
    cur = r;
  }
  if (b < 1.0) {
    double u = a > 0.0 ? (b - a) / (1.0 - a) : b;
    if (!(u > 0.0)) {
      return 0;
    }
    if (u < 1.0) {
      Seg l, r;
      seg_split(&cur, u, &l, &r);
      cur = l;
    }
  }
  *out = cur;
  return piece_ok(&cur);
}

/** The point at t. */
static Vec seg_point(const Seg * s, double t) {
  if (t <= 0.0) {
    return s->p[0];
  }
  if (t >= 1.0) {
    return s->p[s->n];
  }
  double x[4], y[4], w[4];
  for (int i = 0; i <= s->n; i++) {
    w[i] = s->rational ? s->w[i] : 1.0;
    x[i] = s->rational ? s->w[i] * s->p[i].x : s->p[i].x;
    y[i] = s->rational ? s->w[i] * s->p[i].y : s->p[i].y;
  }
  double X = dc(x, s->n + 1, t), Y = dc(y, s->n + 1, t);
  if (!s->rational) {
    return gshgeom_sk_v(X, Y);
  }
  double W = dc(w, s->n + 1, t);
  return gshgeom_sk_v(X / W, Y / W);
}

/**
 * The direction of the centerline at t: SG5's start direction at 0 and end direction at 1, so the exact input
 * differences where a cut is a vertex; inside, the hodograph (for a rational kind `N'D - ND'`, which has the direction
 * of `B'`), from the segment's own control points and not from a cut piece's; where that vanishes, the start direction
 * of the piece after t.
 */
static Vec seg_dir(const Seg * s, double t) {
  Vec d;
  if (t <= 0.0) {
    (void)gshgeom_sk_start_dir(s, &d);
    return d;
  }
  if (t >= 1.0) {
    (void)gshgeom_sk_end_dir(s, &d);
    return d;
  }
  int n = s->n;
  double dx[4], dy[4], dw[4], x[4], y[4], w[4];
  for (int i = 0; i <= n; i++) {
    w[i] = s->rational ? s->w[i] : 1.0;
    x[i] = s->rational ? s->w[i] * s->p[i].x : s->p[i].x;
    y[i] = s->rational ? s->w[i] * s->p[i].y : s->p[i].y;
  }
  for (int i = 0; i < n; i++) {
    dx[i] = x[i + 1] - x[i];
    dy[i] = y[i + 1] - y[i];
    dw[i] = w[i + 1] - w[i];
  }
  if (!s->rational) {
    d = gshgeom_sk_v(dc(dx, n, t), dc(dy, n, t));
  } else {
    double N0 = dc(x, n + 1, t), N1 = dc(y, n + 1, t), D = dc(w, n + 1, t);
    double Nd0 = dc(dx, n, t), Nd1 = dc(dy, n, t), Dd = dc(dw, n, t);
    d = gshgeom_sk_v(Nd0 * D - N0 * Dd, Nd1 * D - N1 * Dd);
  }
  if ((d.x != 0.0 || d.y != 0.0) && d.x - d.x == 0.0 && d.y - d.y == 0.0) {
    return d;
  }
  Seg l, r;
  seg_split(s, t, &l, &r);
  if (gshgeom_sk_start_dir(&r, &d) >= 0) {
    return d;
  }
  if (gshgeom_sk_end_dir(&l, &d) >= 0) {
    return d;
  }
  (void)gshgeom_sk_start_dir(s, &d);
  return d;
}

/** What one subpath's dashes need. */
typedef struct Cb {
  GSHGEOM_SK_Ctx * c;
  const GSHGEOM_SK_Reader * rd;
  size_t sub;
  double begin, end; ///< the subpath's distances on the path
} Cb;

/** The segment and parameter at a distance (SG20): a boundary belongs to the later segment, the end to the last. */
static GSHGEOM_Result locate(const Cb * cb, double s, size_t * k, double * t) {
  const Sub * sub = &cb->rd->subs[cb->sub];
  size_t f = sub->first_seg, last = f + sub->seg_count - 1;
  if (s >= cb->end) {
    *k = last;
    *t = 1.0;
    return GSHGEOM_OK;
  }
  if (s <= cb->begin) {
    *k = f;
    *t = 0.0;
    return GSHGEOM_OK;
  }
  GSHGEOM_Measured_Frame fr;
  if (gshgeom_measured_at(cb->c->dash.m, s, &fr) != GSHGEOM_OK) {
    return GSHGEOM_ERR_INTERNAL;
  }
  size_t lo = f, hi = last + 1;
  while (lo < hi) {
    size_t mid = lo + (hi - lo) / 2;
    if (cb->rd->segs[mid].elem < fr.element) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }
  if (lo > last || cb->rd->segs[lo].elem != fr.element) {
    return GSHGEOM_ERR_INTERNAL;
  }
  *k = lo;
  *t = fr.t;
  return GSHGEOM_OK;
}

/** A zero-length dash at (k, t) (SG24). */
static GSHGEOM_Result dot(const Cb * cb, size_t k, double t) {
  const Seg * s = &cb->rd->segs[k];
  return gshgeom_sk_degenerate(cb->c, seg_point(s, t), seg_dir(s, t));
}

/** One dash <a, b> of the subpath (SG22 to SG25). */
static GSHGEOM_Result on_dash(void * ctx, double a, double b) {
  Cb * cb = ctx;
  GSHGEOM_SK_Ctx * c = cb->c;
  const Sub * sub = &cb->rd->subs[cb->sub];
  if (++c->dash.pieces > c->work_cap) {
    return GSHGEOM_ERR_LIMIT;
  }
  if (sub->seg_count == 0) {
    /* a degenerate subpath: SG11's shape, the pattern being on at its point */
    return gshgeom_sk_degenerate(c, sub->start, gshgeom_sk_effective_tangent(c, cb->rd, cb->sub));
  }
  const Seg * segs = cb->rd->segs;
  size_t f = sub->first_seg, last = f + sub->seg_count - 1;
  double s0 = cb->begin + a, s1 = cb->begin + b;
  if (s1 > cb->end) {
    s1 = cb->end;
  }
  size_t k0, k1;
  double t0, t1;
  GSHGEOM_Result r = locate(cb, s0, &k0, &t0);
  if (r != GSHGEOM_OK) {
    return r;
  }
  if (t0 >= 1.0 && k0 < last) {
    k0++;
    t0 = 0.0;
  }
  if (!(s1 > s0)) {
    return dot(cb, k0, t0);
  }
  r = locate(cb, s1, &k1, &t1);
  if (r != GSHGEOM_OK) {
    return r;
  }
  if (t1 <= 0.0 && k1 > f) {
    k1--;
    t1 = 1.0;
  }
  if (k1 < k0 || (k1 == k0 && !(t1 > t0))) {
    return dot(cb, k0, t0);
  }
  /* the first and last spans may collapse under the cut: the dash then starts or ends at a vertex */
  Seg first, lastp;
  for (;;) {
    if (piece(&segs[k0], t0, k0 == k1 ? t1 : 1.0, &first)) {
      break;
    }
    if (k0 == k1) {
      return dot(cb, k0, t0); /* nothing of the dash has extent: a dot at its start */
    }
    k0++;
    t0 = 0.0;
  }
  if (k1 == k0) {
    lastp = first;
  } else {
    for (;;) {
      if (piece(&segs[k1], 0.0, t1, &lastp)) {
        break;
      }
      k1--;
      t1 = 1.0;
      if (k1 == k0) {
        lastp = first; /* computed to the end of its segment */
        break;
      }
    }
  }
  /* the start cap, the spans with the joins between, the end cap (SG25) */
  Vec d0 = seg_dir(&segs[k0], t0);
  r = gshgeom_sk_cap(c, first.p[0], gshgeom_sk_scale(-1.0, d0));
  for (size_t k = k0; k <= k1 && r == GSHGEOM_OK; k++) {
    if (k == k0) {
      r = gshgeom_sk_body(c, &first);
    } else if (k == k1) {
      r = gshgeom_sk_body(c, &lastp);
    } else {
      r = gshgeom_sk_body(c, &segs[k]);
    }
    if (r == GSHGEOM_OK && k < k1) {
      r = gshgeom_sk_join(c, &segs[k], &segs[k + 1]);
    }
  }
  if (r == GSHGEOM_OK) {
    r = gshgeom_sk_cap(c, lastp.p[lastp.n], seg_dir(&segs[k1], t1));
  }
  return r;
}

GSHGEOM_Result gshgeom_sk_dash_begin(GSHGEOM_SK_Ctx * c, const GSHGEOM_Path * centerline,
    const GSHGEOM_SK_Reader * rd, const GSHGEOM_Stroke * stroke) {
  GSHGEOM_SK_Dash * D = &c->dash;
  size_t n = c->st.dash_count;
  D->dashes = gcu_allocator_calloc(&c->alloc, n, sizeof *D->dashes);
  if (!D->dashes) {
    return GSHGEOM_ERR_OOM;
  }
  for (size_t i = 0; i < n; i++) {
    double v = stroke->dashes[i % stroke->dash_count];
    D->dashes[i] = v == 0.0 ? 0.0 : v;
  }
  D->n = n;
  D->period = gshgeom_sk_dash_period(D->dashes, n, n);
  D->offset = c->st.dash_offset;
  D->per_path = c->st.dash_scope == (uint32_t)GSHGEOM_STROKE_DASH_PATH;
  D->phase0 = gshgeom_sk_dash_phase(D->offset, D->period);
  const GSHGEOM_Path * mp = centerline;
  GSHGEOM_Result r = GSHGEOM_OK;
  if (!c->local) {
    /* the mapped centerline, element for element: its arc length is the one the pattern is measured in */
    r = gshgeom_path_create(&c->alloc, &c->limits, 2, &D->bpath);
    size_t ne = gshgeom_path_element_count(centerline);
    for (size_t i = 0; i < ne && r == GSHGEOM_OK; i++) {
      GSHGEOM_Element el;
      r = gshgeom_path_element(centerline, i, &el);
      if (r != GSHGEOM_OK) {
        break;
      }
      GSHGEOM_Curve cv;
      r = gshgeom_map_apply_curve(&c->map, &el.curve, &cv);
      if (r != GSHGEOM_OK) {
        break;
      }
      const GSHGEOM_Point * p = cv.p;
      switch (el.verb) {
      case GSHGEOM_VERB_MOVE:
        r = gshgeom_path_move_to(D->bpath, p[0]);
        break;
      case GSHGEOM_VERB_LINE:
        r = gshgeom_path_line_to(D->bpath, p[1]);
        break;
      case GSHGEOM_VERB_QUADRATIC:
        r = gshgeom_path_quadratic_to(D->bpath, p[1], p[2]);
        break;
      case GSHGEOM_VERB_CUBIC:
        r = gshgeom_path_cubic_to(D->bpath, p[1], p[2], p[3]);
        break;
      case GSHGEOM_VERB_RATIONAL_QUADRATIC:
        r = gshgeom_path_rational_quadratic_to(D->bpath, p[1], p[2], cv.w);
        break;
      case GSHGEOM_VERB_RATIONAL_CUBIC:
        r = gshgeom_path_rational_cubic_to(D->bpath, p[1], p[2], p[3], cv.w);
        break;
      case GSHGEOM_VERB_CLOSE:
        r = gshgeom_path_close(D->bpath);
        break;
      default:
        r = GSHGEOM_ERR_INTERNAL;
        break;
      }
    }
    mp = D->bpath;
  }
  if (r == GSHGEOM_OK) {
    r = gshgeom_measured_create(&c->alloc, &c->limits, mp, GSHGEOM_MEASURE_TOLERANCE_MIN, &D->m);
  }
  if (r != GSHGEOM_OK) {
    return r;
  }
  size_t draws = 0;
  for (size_t i = 0; i < rd->sub_n; i++) {
    draws += rd->subs[i].draws ? 1u : 0u;
  }
  if (draws != gshgeom_measured_subpath_count(D->m)) {
    return GSHGEOM_ERR_INTERNAL;
  }
  /* SG26: the cut position's error bound, and its gate */
  double poly = 0.0;
  for (size_t i = 0; i < rd->seg_n; i++) {
    const Seg * s = &rd->segs[i];
    for (int k = 0; k < s->n; k++) {
      poly += gshgeom_sk_norm(gshgeom_sk_sub(s->p[k + 1], s->p[k]));
    }
  }
  double L = gshgeom_measured_length(D->m);
  double rows = (double)gshgeom_measured_row_count(D->m);
  double eps = (0x1p-46 * poly + (rows + 4.0) * 0x1p-53 * L + 0x1p-52 * L +
                   ((double)n + 8.0) * 0x1p-53 * (L + 2.0 * D->period)) *
      (1.0 + 0x1p-30);
  if (!(eps * c->norm <= 0x1p-20)) {
    return GSHGEOM_ERR_LIMIT;
  }
  D->active = 1;
  return GSHGEOM_OK;
}

void gshgeom_sk_dash_end(GSHGEOM_SK_Ctx * c) {
  GSHGEOM_SK_Dash * D = &c->dash;
  gshgeom_measured_destroy(D->m);
  gshgeom_path_destroy(D->bpath);
  gcu_allocator_free(&c->alloc, D->dashes);
  D->m = NULL;
  D->bpath = NULL;
  D->dashes = NULL;
  D->active = 0;
}

GSHGEOM_Result gshgeom_sk_dash_sub(GSHGEOM_SK_Ctx * c, const GSHGEOM_SK_Reader * rd, size_t i) {
  const Sub * sub = &rd->subs[i];
  GSHGEOM_SK_Dash * D = &c->dash;
  if (!sub->draws) {
    return GSHGEOM_OK; /* a lone move */
  }
  GSHGEOM_Measured_Subpath ms;
  if (gshgeom_measured_subpath(D->m, D->msub++, &ms) != GSHGEOM_OK) {
    return GSHGEOM_ERR_INTERNAL;
  }
  if (sub->seg_count == 0 && c->st.degenerate == (uint32_t)GSHGEOM_STROKE_DEGENERATE_OMIT) {
    return GSHGEOM_OK;
  }
  Cb cb = {c, rd, i, ms.begin, ms.end};
  double o = D->per_path ? gshgeom_sk_dash_phase(D->offset + ms.begin, D->period) : D->phase0;
  return gshgeom_sk_dash_positions(D->dashes, D->n, D->n, D->period, o, ms.end - ms.begin, &D->steps, c->work_cap,
      on_dash, &cb);
}
