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
 * The copies of a repeated cell (lattice.h, LT1 to LT5).
 *
 * **The set.**  A copy `(i, j)` is a copy of the use when the convex hull `H_ij = conv(Q_k + w_ij)` of its four exact corners
 * meets the closed rectangle `R`, with `w_ij = tau + i A + j B` formed exactly.  `H_ij` meets `R` exactly when `w_ij` lies in the
 * convex polygon `K = R - H_00` (the Minkowski difference), and `K` is the intersection of the half-planes `u . w <= h_K(u)`
 * over the edge normals `u` of `R` and of `H_00` (the support function of a Minkowski sum is the sum of the support
 * functions).  The edge normals of `H_00` are among the normals of the six segments between the four points of `Q`, both
 * signs, so the sixteen directions used here (the two axes, both signs, and those twelve) give `K` and any more would
 * give the same set.  Each half-plane is `c i + e j + g <= 0`, with `c`, `e` and `g` exact sums of products of doubles.
 *
 * **The search.**  `S(i, j) = c i + e j + g` is evaluated exactly, by a filter that decides the clear cases from a
 * rounded evaluation and its error bound and by the exact sum for the rest.  For a fixed row `S <= 0` is monotone in
 * `i`, so the last (or first) feasible `i` is found from an estimate in doubles by galloping (steps 1, 2, 4, ...) to a
 * change of sign and bisecting: the estimate only decides how many steps it takes, never what is found.
 *
 * Nothing here uses libm, `long double`, an integer to double cast of a 64-bit value (which x87 hosts make with x87
 * instructions) or a fused multiply-add.
 */

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/shape-geometry/hash.h>
#include <ghoti.io/shape-geometry/lattice.h>
#include <ghoti.io/shape-geometry/map.h>
#include <ghoti.io/shape-geometry/stroke.h>

#include "../exact/exact.h"
#include "lattice_internal.h"

#include <stddef.h>
#include <stdint.h>

#define INDEX_MAX GSHGEOM_LATTICE_INDEX_MAX
#define TWO53 ((int64_t)1 << 53)

/* --- integers and doubles without a 64-bit conversion ------------------------------------------------------------------ */

/** The double of an integer of magnitude at most 2^53 (exact). */
static double i2d(int64_t v) {
  const int32_t hi = (int32_t)(v >> 32);
  const uint32_t lo = (uint32_t)((uint64_t)v & 0xffffffffu);
  return (double)hi * 4294967296.0 + (double)lo;
}

/** @p x toward zero as an integer; @p x is finite with magnitude below 2^53 (the callers clamp first). */
static int64_t d2i_trunc(double x) {
  const uint64_t bits = gshgeom_f64_bits(x);
  const int neg = (int)(bits >> 63);
  const int ex = (int)((bits >> 52) & 0x7ffu);
  if (ex < 1023) {
    return 0;
  }
  const int e = ex - 1075; /* x = m * 2^e with m the 53-bit mantissa */
  const uint64_t m = (bits & ((UINT64_C(1) << 52) - 1u)) | (UINT64_C(1) << 52);
  const uint64_t v = e >= 0 ? m << e : m >> (-e);
  return neg ? -(int64_t)v : (int64_t)v;
}

/** floor(@p x) clamped to [-INDEX_MAX, INDEX_MAX]; a NaN is 0. */
static int64_t d2i_floor_clamped(double x) {
  if (!(x == x)) {
    return 0;
  }
  if (x >= i2d(INDEX_MAX)) {
    return INDEX_MAX;
  }
  if (x <= -i2d(INDEX_MAX)) {
    return -INDEX_MAX;
  }
  int64_t t = d2i_trunc(x);
  if (i2d(t) > x) {
    t--;
  }
  return t < -INDEX_MAX ? -INDEX_MAX : t;
}

static int64_t floor_mod(int64_t a, int64_t p) {
  int64_t m = a % p;
  return m < 0 ? m + p : m;
}

static int64_t floor_div(int64_t a, int64_t p) {
  int64_t q = a / p;
  return (a % p != 0 && a < 0) ? q - 1 : q;
}

/* --- the lattice and the map ------------------------------------------------------------------------------------------- */

static int is_fin(double v) {
  return v - v == 0.0;
}

static int lattice_valid(const GSHGEOM_Lattice * l) {
  if (!l || !is_fin(l->x) || !is_fin(l->y) || !is_fin(l->w) || !is_fin(l->h) || !is_fin(l->x_step) || !is_fin(l->y_step) ||
      !(l->w > 0.0) || !(l->h > 0.0) || l->x_step < 0.0 || l->y_step < 0.0 || l->repeat_x > 1u || l->repeat_y > 1u) {
    return 0;
  }
  return (l->repeat_x == 0u || l->x_step > 0.0) && (l->repeat_y == 0u || l->y_step > 0.0);
}

GSHGEOM_Result gshgeom_lattice_cell(const GSHGEOM_Lattice * l, int64_t i, int64_t j, double out[4]) {
  if (!l || !out || i < -INDEX_MAX || i > INDEX_MAX || j < -INDEX_MAX || j > INDEX_MAX || !lattice_valid(l)) {
    return GSHGEOM_ERR_INVALID;
  }
  const double left = l->x + i2d(i) * l->x_step;
  const double top = l->y + i2d(j) * l->y_step;
  const double right = l->x_step == l->w ? l->x + i2d(i + 1) * l->x_step : left + l->w;
  const double bottom = l->y_step == l->h ? l->y + i2d(j + 1) * l->y_step : top + l->h;
  if (!is_fin(left) || !is_fin(top) || !is_fin(right) || !is_fin(bottom)) {
    return GSHGEOM_ERR_INVALID;
  }
  out[0] = left;
  out[1] = top;
  out[2] = right;
  out[3] = bottom;
  return GSHGEOM_OK;
}

static int map_ok(const GSHGEOM_Map * m) {
  return m && m->dim == 2 && gshgeom_map_check(m) == GSHGEOM_OK;
}

/** The step vectors of LT2: a = (fl(m00 x_step), fl(m10 x_step)), b = (fl(m01 y_step), fl(m11 y_step)). */
static void steps_of(const GSHGEOM_Map * m, const GSHGEOM_Lattice * l, double a[2], double b[2]) {
  a[0] = m->m[0][0] * l->x_step;
  a[1] = m->m[1][0] * l->x_step;
  b[0] = m->m[0][1] * l->y_step;
  b[1] = m->m[1][1] * l->y_step;
}

GSHGEOM_Result gshgeom_lattice_translation(const GSHGEOM_Map * m, const GSHGEOM_Lattice * l, int64_t i, int64_t j, double out[2]) {
  if (!out || i < -INDEX_MAX || i > INDEX_MAX || j < -INDEX_MAX || j > INDEX_MAX || !map_ok(m) || !lattice_valid(l)) {
    return GSHGEOM_ERR_INVALID;
  }
  double a[2];
  double b[2];
  steps_of(m, l, a, b);
  const double di = i2d(i);
  const double dj = i2d(j);
  const double x = (m->t[0] + di * a[0]) + dj * b[0];
  const double y = (m->t[1] + di * a[1]) + dj * b[1];
  if (!is_fin(x) || !is_fin(y)) {
    return GSHGEOM_ERR_INVALID;
  }
  out[0] = x;
  out[1] = y;
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_lattice_quad(const GSHGEOM_Map * m, const GSHGEOM_Lattice * l, double out[8]) {
  double c[4];
  if (!out || !map_ok(m) || gshgeom_lattice_cell(l, 0, 0, c) != GSHGEOM_OK) {
    return GSHGEOM_ERR_INVALID;
  }
  const double px[4] = {c[0], c[2], c[2], c[0]};
  const double py[4] = {c[1], c[1], c[3], c[3]};
  double q[8];
  for (int k = 0; k < 4; k++) {
    q[2 * k] = m->m[0][0] * px[k] + m->m[0][1] * py[k];
    q[2 * k + 1] = m->m[1][0] * px[k] + m->m[1][1] * py[k];
    if (!is_fin(q[2 * k]) || !is_fin(q[2 * k + 1])) {
      return GSHGEOM_ERR_INVALID;
    }
  }
  for (int k = 0; k < 8; k++) {
    out[k] = q[k];
  }
  return GSHGEOM_OK;
}

/* --- exact constraints ------------------------------------------------------------------------------------------------- */

/** A direction `u = (ux1 - ux0, uy1 - uy0)`, kept as the differences' operands so that it is exact. */
typedef struct Dir {
  double ux1, ux0, uy1, uy0;
} Dir;

/** What the exact evaluation reads. */
typedef struct Ctx {
  double tau[2];
  double A[2];
  double B[2];
  double Q[4][2];
  double R[4][2]; /* the corners of the rectangle, in the order (x0,y0), (x1,y0), (x1,y1), (x0,y1) */
} Ctx;

/** One half-plane `c i + e j + g <= 0`. */
typedef struct Con {
  Dir u;
  double r[2];      /* the corner of R that maximizes u . r */
  double q[2];      /* the point of Q that minimizes u . Q */
  double tau[2];
  double A[2];
  double B[2];
  double c, e, g;   /* the coefficients, each rounded once from the exact sum */
  int cs, es;       /* the exact signs of c and of e */
  int filter;       /* c, e and g are in the range where the rounded evaluation has a proved error bound */
} Con;

/** Append the terms of `(-1)^neg (u . p) m` (m omitted when @p has_m is 0). */
static void dot_terms(GSHGEOM_ExactTerm * t, size_t * n, int neg, const Dir * u, double px, double py, int has_m, double m) {
  const int f = has_m ? 3 : 2;
  gshgeom_exact_put(t, n, neg, f, u->ux1, px, m);
  gshgeom_exact_put(t, n, !neg, f, u->ux0, px, m);
  gshgeom_exact_put(t, n, neg, f, u->uy1, py, m);
  gshgeom_exact_put(t, n, !neg, f, u->uy0, py, m);
}

static void dot_sum(GSHGEOM_ExactSum * s, int neg, const Dir * u, double px, double py) {
  const double f1[2] = {u->ux1, px};
  const double f2[2] = {u->ux0, px};
  const double f3[2] = {u->uy1, py};
  const double f4[2] = {u->uy0, py};
  gshgeom_exact_sum_add(s, neg, 2, f1);
  gshgeom_exact_sum_add(s, !neg, 2, f2);
  gshgeom_exact_sum_add(s, neg, 2, f3);
  gshgeom_exact_sum_add(s, !neg, 2, f4);
}

static int in_filter_range(double v) {
  if (!is_fin(v)) {
    return 0;
  }
  const double a = v < 0.0 ? -v : v;
  return a == 0.0 || (a >= 0x1p-900 && a <= 0x1p900);
}

static void con_build(const Ctx * cx, Dir u, Con * k) {
  k->u = u;
  /* h(u): the corner of R with the greatest u . r, the point of Q with the least u . Q, each decided exactly. */
  int best = 0;
  for (int c = 1; c < 4; c++) {
    GSHGEOM_ExactTerm t[8];
    size_t n = 0;
    dot_terms(t, &n, 0, &u, cx->R[c][0], cx->R[c][1], 0, 0.0);
    dot_terms(t, &n, 1, &u, cx->R[best][0], cx->R[best][1], 0, 0.0);
    if (gshgeom_exact_sign(t, n) > 0) {
      best = c;
    }
  }
  k->r[0] = cx->R[best][0];
  k->r[1] = cx->R[best][1];
  best = 0;
  for (int c = 1; c < 4; c++) {
    GSHGEOM_ExactTerm t[8];
    size_t n = 0;
    dot_terms(t, &n, 0, &u, cx->Q[c][0], cx->Q[c][1], 0, 0.0);
    dot_terms(t, &n, 1, &u, cx->Q[best][0], cx->Q[best][1], 0, 0.0);
    if (gshgeom_exact_sign(t, n) < 0) {
      best = c;
    }
  }
  k->q[0] = cx->Q[best][0];
  k->q[1] = cx->Q[best][1];
  for (int d = 0; d < 2; d++) {
    k->tau[d] = cx->tau[d];
    k->A[d] = cx->A[d];
    k->B[d] = cx->B[d];
  }
  GSHGEOM_ExactSum sum;
  gshgeom_exact_sum_init(&sum);
  dot_sum(&sum, 0, &u, cx->A[0], cx->A[1]);
  k->c = gshgeom_exact_sum_to_double(&sum);
  gshgeom_exact_sum_init(&sum);
  dot_sum(&sum, 0, &u, cx->B[0], cx->B[1]);
  k->e = gshgeom_exact_sum_to_double(&sum);
  gshgeom_exact_sum_init(&sum);
  dot_sum(&sum, 0, &u, cx->tau[0], cx->tau[1]);
  dot_sum(&sum, 1, &u, k->r[0], k->r[1]);
  dot_sum(&sum, 0, &u, k->q[0], k->q[1]);
  k->g = gshgeom_exact_sum_to_double(&sum);
  GSHGEOM_ExactTerm t[4];
  size_t n = 0;
  dot_terms(t, &n, 0, &u, cx->A[0], cx->A[1], 0, 0.0);
  k->cs = gshgeom_exact_sign(t, n);
  n = 0;
  dot_terms(t, &n, 0, &u, cx->B[0], cx->B[1], 0, 0.0);
  k->es = gshgeom_exact_sign(t, n);
  /* A coefficient whose exact sign is not zero but whose rounded value is must not be trusted (underflow). */
  k->filter = in_filter_range(k->c) && in_filter_range(k->e) && in_filter_range(k->g) && (k->c != 0.0 || k->cs == 0) &&
      (k->e != 0.0 || k->es == 0);
}

/** The sign of `S(i, j) = u . tau + i (u . A) + j (u . B) - u . r + u . q`, exactly. */
static int con_sign(const Con * k, int64_t i, int64_t j) {
  const double di = i2d(i);
  const double dj = i2d(j);
  if (k->filter) {
    const double ci = k->c * di;
    const double ej = k->e * dj;
    const double s = (ci + ej) + k->g;
    const double abs_sum = (ci < 0.0 ? -ci : ci) + (ej < 0.0 ? -ej : ej) + (k->g < 0.0 ? -k->g : k->g);
    const double bound = abs_sum * 0x1p-48; /* five roundings of at most 2^-53 each, a wide margin */
    if (s > bound) {
      return 1;
    }
    if (s < -bound) {
      return -1;
    }
  }
  GSHGEOM_ExactTerm t[20];
  size_t n = 0;
  dot_terms(t, &n, 0, &k->u, k->tau[0], k->tau[1], 0, 0.0);
  dot_terms(t, &n, 0, &k->u, k->A[0], k->A[1], 1, di);
  dot_terms(t, &n, 0, &k->u, k->B[0], k->B[1], 1, dj);
  dot_terms(t, &n, 1, &k->u, k->r[0], k->r[1], 0, 0.0);
  dot_terms(t, &n, 0, &k->u, k->q[0], k->q[1], 0, 0.0);
  return gshgeom_exact_sign(t, n);
}

/* --- the interval of one variable -------------------------------------------------------------------------------------- */

typedef struct Probe {
  const Con * k;
  int var;        /* 0: the variable is i; 1: it is j */
  int64_t other;  /* the other index */
  int dir;        /* +1: feasible for small values; -1: for large ones (t = -v) */
  uint64_t * work; /* counts the evaluations, or NULL */
} Probe;

static int feasible(const Probe * p, int64_t t) {
  if (p->work) {
    ++*p->work;
  }
  const int64_t v = p->dir > 0 ? t : -t;
  return con_sign(p->k, p->var == 0 ? v : p->other, p->var == 0 ? p->other : v) <= 0;
}

/**
 * The greatest @p t in [-INDEX_MAX, INDEX_MAX] with `feasible(t)` (true for small t, false for large: monotone), or
 * -INDEX_MAX - 1 when there is none; @p est is the guess it starts from.
 */
static int64_t last_feasible(const Probe * p, int64_t est) {
  int64_t good;
  int64_t bad;
  if (feasible(p, est)) {
    good = est;
    bad = INDEX_MAX + 1;
    int64_t step = 1;
    while (good < INDEX_MAX) {
      const int64_t cand = good + step > INDEX_MAX ? INDEX_MAX : good + step;
      if (feasible(p, cand)) {
        good = cand;
        step *= 2;
      } else {
        bad = cand;
        break;
      }
    }
    if (bad == INDEX_MAX + 1) {
      return INDEX_MAX;
    }
  } else {
    bad = est;
    good = -INDEX_MAX - 1;
    int64_t step = 1;
    while (bad > -INDEX_MAX) {
      const int64_t cand = bad - step < -INDEX_MAX ? -INDEX_MAX : bad - step;
      if (!feasible(p, cand)) {
        bad = cand;
        step *= 2;
      } else {
        good = cand;
        break;
      }
    }
    if (good == -INDEX_MAX - 1) {
      return good;
    }
  }
  while (bad - good > 1) {
    const int64_t mid = good + (bad - good) / 2;
    if (feasible(p, mid)) {
      good = mid;
    } else {
      bad = mid;
    }
  }
  return good;
}

/**
 * The integers `v` of the variable @p var (0: i, 1: j), with the other index @p other, that satisfy every constraint: an
 * interval, clamped to the index range.  Returns 0 when there are none.
 */
static int interval(const Con * ks, size_t n, int var, int64_t other, int64_t * lo, int64_t * hi, uint64_t * work) {
  int64_t l = -INDEX_MAX;
  int64_t h = INDEX_MAX;
  for (size_t m = 0; m < n; m++) {
    const Con * k = &ks[m];
    const double cv = var == 0 ? k->c : k->e;
    const double co = var == 0 ? k->e : k->c;
    const int cvs = var == 0 ? k->cs : k->es;
    if (cvs == 0) {
      /* The variable does not appear: the row is feasible or not as a whole. */
      if (work) {
        ++*work;
      }
      if (con_sign(k, var == 0 ? 0 : other, var == 0 ? other : 0) > 0) {
        return 0;
      }
      continue;
    }
    Probe p;
    p.k = k;
    p.var = var;
    p.other = other;
    p.dir = cvs > 0 ? 1 : -1;
    p.work = work;
    const double vstar = -(co * i2d(other) + k->g) / cv;
    const int64_t est = d2i_floor_clamped(p.dir > 0 ? vstar : -vstar);
    /* `t` is -INDEX_MAX - 1 when no index is feasible: then `h` falls below `l` (or `l` rises past `h`) and the test below ends the search. */
    const int64_t t = last_feasible(&p, est);
    if (p.dir > 0) {
      h = t < h ? t : h;
    } else {
      l = -t > l ? -t : l;
    }
    /* An interval that is empty is not searched further: each of the constraints that follow costs a gallop from its estimate. */
    if (l > h) {
      return 0;
    }
  }
  *lo = l;
  *hi = h;
  return l <= h;
}

/* --- the enumeration --------------------------------------------------------------------------------------------------- */

typedef struct Emit {
  GSHGEOM_Lattice_Row_Fn fn;
  void * ctx;
  uint64_t max;
  uint64_t count;
} Emit;

enum { EMIT_GO = 0, EMIT_STOP = 1, EMIT_LIMIT = 2 };

/** Report the row `(j, lo .. hi)`; the count of the copies reported so far is kept. */
static int emit_row(Emit * e, int64_t j, int64_t lo, int64_t hi) {
  const uint64_t len = (uint64_t)(hi - lo) + 1u;
  const uint64_t room = e->max - e->count;
  if (len > room) {
    if (room > 0 && e->fn) {
      (void)e->fn(e->ctx, j, lo, lo + (int64_t)room - 1);
    }
    e->count = e->max + 1u;
    return EMIT_LIMIT;
  }
  e->count += len;
  if (e->fn && e->fn(e->ctx, j, lo, hi)) {
    return EMIT_STOP;
  }
  return EMIT_GO;
}

static void dir_set(Dir * d, double ux1, double ux0, double uy1, double uy0) {
  d->ux1 = ux1;
  d->ux0 = ux0;
  d->uy1 = uy1;
  d->uy0 = uy0;
}

GSHGEOM_Result gshgeom_lattice_each_counted(const GSHGEOM_Map * map, const GSHGEOM_Lattice * l, const int64_t rect[4], uint64_t max_count,
    GSHGEOM_Lattice_Row_Fn row_fn, void * ctx, uint64_t * count, GSHGEOM_Lattice_Stats * stats) {
  uint64_t * work = stats ? &stats->evaluations : NULL;
  if (stats) {
    stats->evaluations = 0;
  }
  if (!count || !rect || !map_ok(map) || !lattice_valid(l) || max_count > (uint64_t)INDEX_MAX) {
    return GSHGEOM_ERR_INVALID;
  }
  for (int k = 0; k < 4; k++) {
    if (rect[k] > TWO53 || rect[k] < -TWO53) {
      return GSHGEOM_ERR_INVALID;
    }
  }
  if (rect[0] > rect[2] || rect[1] > rect[3]) {
    return GSHGEOM_ERR_INVALID;
  }
  {
    /* The linear part must be nonsingular, exactly. */
    GSHGEOM_ExactTerm t[2];
    size_t n = 0;
    gshgeom_exact_put(t, &n, 0, 2, map->m[0][0], map->m[1][1], 0.0);
    gshgeom_exact_put(t, &n, 1, 2, map->m[0][1], map->m[1][0], 0.0);
    if (gshgeom_exact_sign(t, n) == 0) {
      return GSHGEOM_ERR_INVALID;
    }
  }
  Ctx cx;
  double quad[8];
  if (gshgeom_lattice_quad(map, l, quad) != GSHGEOM_OK) {
    return GSHGEOM_ERR_INVALID;
  }
  for (int k = 0; k < 4; k++) {
    cx.Q[k][0] = quad[2 * k];
    cx.Q[k][1] = quad[2 * k + 1];
  }
  steps_of(map, l, cx.A, cx.B);
  cx.tau[0] = map->t[0];
  cx.tau[1] = map->t[1];
  if (!is_fin(cx.A[0]) || !is_fin(cx.A[1]) || !is_fin(cx.B[0]) || !is_fin(cx.B[1])) {
    return GSHGEOM_ERR_INVALID;
  }
  const double x0 = i2d(rect[0]);
  const double y0 = i2d(rect[1]);
  const double x1 = i2d(rect[2]);
  const double y1 = i2d(rect[3]);
  const double rc[4][2] = {{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}};
  for (int k = 0; k < 4; k++) {
    cx.R[k][0] = rc[k][0];
    cx.R[k][1] = rc[k][1];
  }
  const int rx = l->repeat_x != 0u;
  const int ry = l->repeat_y != 0u;
  if (rx && ry) {
    GSHGEOM_ExactTerm t[2];
    size_t n = 0;
    gshgeom_exact_put(t, &n, 0, 2, cx.A[0], cx.B[1], 0.0);
    gshgeom_exact_put(t, &n, 1, 2, cx.A[1], cx.B[0], 0.0);
    if (gshgeom_exact_sign(t, n) == 0) {
      return GSHGEOM_ERR_INVALID; /* parallel step vectors: infinitely many copies on a line */
    }
  }
  /* The sixteen directions: the axes and the normals of the six segments of Q, both signs. */
  Con ks[16];
  size_t nk = 0;
  {
    Dir d;
    dir_set(&d, 1.0, 0.0, 0.0, 0.0);
    con_build(&cx, d, &ks[nk++]);
    dir_set(&d, 0.0, 1.0, 0.0, 0.0);
    con_build(&cx, d, &ks[nk++]);
    dir_set(&d, 0.0, 0.0, 1.0, 0.0);
    con_build(&cx, d, &ks[nk++]);
    dir_set(&d, 0.0, 0.0, 0.0, 1.0);
    con_build(&cx, d, &ks[nk++]);
    for (int a = 0; a < 4; a++) {
      for (int b = a + 1; b < 4; b++) {
        if (cx.Q[a][0] == cx.Q[b][0] && cx.Q[a][1] == cx.Q[b][1]) {
          continue; /* the direction zero: its half-plane `0 <= 0` bounds nothing, and each one is a test of every row */
        }
        /* n = (-(by - ay), bx - ax) and -n */
        dir_set(&d, cx.Q[a][1], cx.Q[b][1], cx.Q[b][0], cx.Q[a][0]);
        con_build(&cx, d, &ks[nk++]);
        dir_set(&d, cx.Q[b][1], cx.Q[a][1], cx.Q[a][0], cx.Q[b][0]);
        con_build(&cx, d, &ks[nk++]);
      }
    }
  }
  Emit em = {row_fn, ctx, max_count, 0};
  GSHGEOM_Result result = GSHGEOM_OK;
  if (!rx && !ry) {
    int64_t lo;
    int64_t hi;
    if (interval(ks, nk, 0, 0, &lo, &hi, work) && lo <= 0 && 0 <= hi) {
      (void)emit_row(&em, 0, 0, 0);
    }
  } else if (rx && !ry) {
    int64_t lo;
    int64_t hi;
    if (interval(ks, nk, 0, 0, &lo, &hi, work)) {
      if (emit_row(&em, 0, lo, hi) == EMIT_LIMIT) {
        result = GSHGEOM_ERR_LIMIT;
      }
    }
  } else if (!rx && ry) {
    int64_t lo;
    int64_t hi;
    if (interval(ks, nk, 1, 0, &lo, &hi, work)) {
      const uint64_t len = (uint64_t)(hi - lo) + 1u;
      if (!row_fn) {
        if (len > max_count) {
          em.count = max_count + 1u;
          result = GSHGEOM_ERR_LIMIT;
        } else {
          em.count = len;
        }
      } else {
        for (int64_t j = lo; j <= hi; j++) {
          const int s = emit_row(&em, j, 0, 0);
          if (s == EMIT_LIMIT) {
            result = GSHGEOM_ERR_LIMIT;
            break;
          }
          if (s == EMIT_STOP) {
            break;
          }
        }
      }
    }
  } else {
    /* The rows: the integers j between the two supporting lines of K parallel to A, from the direction (-A_y, A_x). */
    Con vk[2];
    Dir d;
    dir_set(&d, 0.0, cx.A[1], cx.A[0], 0.0);
    con_build(&cx, d, &vk[0]);
    dir_set(&d, cx.A[1], 0.0, 0.0, cx.A[0]);
    con_build(&cx, d, &vk[1]);
    int64_t jlo;
    int64_t jhi;
    if (interval(vk, 2, 1, 0, &jlo, &jhi, work)) {
      const uint64_t budget = 2u * max_count + GSHGEOM_LATTICE_ROW_WORK;
      uint64_t visited = 0;
      for (int64_t j = jlo; j <= jhi; j++) {
        if (++visited > budget) {
          em.count = max_count + 1u;
          result = GSHGEOM_ERR_LIMIT;
          break;
        }
        int64_t lo;
        int64_t hi;
        if (!interval(ks, nk, 0, j, &lo, &hi, work)) {
          continue;
        }
        const int s = emit_row(&em, j, lo, hi);
        if (s == EMIT_LIMIT) {
          result = GSHGEOM_ERR_LIMIT;
          break;
        }
        if (s == EMIT_STOP) {
          break;
        }
      }
    }
  }
  *count = em.count;
  return result;
}

GSHGEOM_Result gshgeom_lattice_each(const GSHGEOM_Map * map, const GSHGEOM_Lattice * l, const int64_t rect[4], uint64_t max_count,
    GSHGEOM_Lattice_Row_Fn row_fn, void * ctx, uint64_t * count) {
  return gshgeom_lattice_each_counted(map, l, rect, max_count, row_fn, ctx, count, NULL);
}

/* --- the integer lattice and the canonical position of a copy (LT5) --------------------------------------------------- */

/** Whether @p v is an integer of magnitude below 2^31; then @p out is it. */
static int exact_int31(double v, int64_t * out) {
  if (!(v > -2147483648.0 && v < 2147483648.0)) {
    return 0;
  }
  const int64_t t = d2i_trunc(v);
  if (i2d(t) != v) {
    return 0;
  }
  *out = t;
  return 1;
}

/** Both device step vectors of @p m and @p l as integers, when every component is an exact integer of magnitude below 2^31. */
static int integral_steps(const GSHGEOM_Map * m, const GSHGEOM_Lattice * l, int64_t a[2], int64_t b[2]) {
  double da[2];
  double db[2];
  steps_of(m, l, da, db);
  return exact_int31(da[0], &a[0]) && exact_int31(da[1], &a[1]) && exact_int31(db[0], &b[0]) && exact_int31(db[1], &b[1]);
}

GSHGEOM_Result gshgeom_lattice_steps(const GSHGEOM_Map * m, const GSHGEOM_Lattice * l, double a[2], double b[2]) {
  if (!a || !b || !map_ok(m) || !lattice_valid(l)) {
    return GSHGEOM_ERR_INVALID;
  }
  steps_of(m, l, a, b);
  if (!is_fin(a[0]) || !is_fin(a[1]) || !is_fin(b[0]) || !is_fin(b[1])) {
    return GSHGEOM_ERR_INVALID;
  }
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_lattice_integer_steps(const GSHGEOM_Map * m, const GSHGEOM_Lattice * l, int64_t a[2], int64_t b[2], int * integer) {
  if (!a || !b || !integer || !map_ok(m) || !lattice_valid(l)) {
    return GSHGEOM_ERR_INVALID;
  }
  *integer = 0;
  int64_t ia[2];
  int64_t ib[2];
  if (l->repeat_x != 1u || l->repeat_y != 1u || !integral_steps(m, l, ia, ib)) {
    return GSHGEOM_OK;
  }
  const int64_t det = ia[0] * ib[1] - ia[1] * ib[0]; /* each product is below 2^62 */
  const int64_t adet = det < 0 ? -det : det;
  if (adet == 0 || adet > GSHGEOM_LATTICE_PERIOD_MAX) {
    return GSHGEOM_OK;
  }
  a[0] = ia[0];
  a[1] = ia[1];
  b[0] = ib[0];
  b[1] = ib[1];
  *integer = 1;
  return GSHGEOM_OK;
}

/** @p i times @p a, or 0 when its magnitude passes 2^40 (a copy whose whole pixels are past what any sink draws). */
static int mul_bounded(int64_t i, int64_t a, int64_t * out) {
  const int64_t limit = (int64_t)1 << 40;
  if (a == 0 || i == 0) {
    *out = 0;
    return 1;
  }
  const int64_t ai = a < 0 ? -a : a;
  const int64_t ii = i < 0 ? -i : i;
  if (ii > limit / ai) {
    return 0;
  }
  *out = i * a;
  return 1;
}

GSHGEOM_Result gshgeom_lattice_canonical(const GSHGEOM_Map * m, const GSHGEOM_Lattice * l, int64_t i, int64_t j, GSHGEOM_Map * canonical,
    int64_t * shift_x, int64_t * shift_y) {
  if (!canonical || !shift_x || !shift_y || i < -INDEX_MAX || i > INDEX_MAX || j < -INDEX_MAX || j > INDEX_MAX || !map_ok(m) ||
      !lattice_valid(l)) {
    return GSHGEOM_ERR_INVALID;
  }
  int64_t ia[2];
  int64_t ib[2];
  if (!integral_steps(m, l, ia, ib)) {
    double w[2];
    GSHGEOM_Map full = *m;
    const GSHGEOM_Result r = gshgeom_lattice_translation(m, l, i, j, w);
    if (r != GSHGEOM_OK) {
      return r;
    }
    full.t[0] = w[0];
    full.t[1] = w[1];
    return gshgeom_stroke_canonical_map(&full, canonical, shift_x, shift_y);
  }
  /* Integer steps: the copy's whole pixels are floor(tau) + i a + j b in integers and its fraction is the fraction of tau, so that
   * no sum is rounded and every copy of the lattice is the same bits shifted whatever the fraction is. */
  int64_t sx = 0;
  int64_t sy = 0;
  GSHGEOM_Map base = *m;
  const GSHGEOM_Result r = gshgeom_stroke_canonical_map(m, &base, &sx, &sy);
  if (r != GSHGEOM_OK) {
    return r;
  }
  int64_t t[4];
  if (!mul_bounded(i, ia[0], &t[0]) || !mul_bounded(j, ib[0], &t[1]) || !mul_bounded(i, ia[1], &t[2]) || !mul_bounded(j, ib[1], &t[3])) {
    return GSHGEOM_ERR_LIMIT;
  }
  const int64_t limit = (int64_t)1 << 40;
  const int64_t x = sx + t[0] + t[1];
  const int64_t y = sy + t[2] + t[3];
  if (x <= -limit || x >= limit || y <= -limit || y >= limit) {
    return GSHGEOM_ERR_LIMIT;
  }
  *canonical = base;
  *shift_x = x;
  *shift_y = y;
  return GSHGEOM_OK;
}

/* --- the period (LT3) -------------------------------------------------------------------------------------------------- */

GSHGEOM_Result gshgeom_lattice_period(const int64_t a[2], const int64_t b[2], GSHGEOM_Lattice_Period * out) {
  const int64_t lim = (int64_t)1 << 31;
  if (!a || !b || !out) {
    return GSHGEOM_ERR_INVALID;
  }
  for (int k = 0; k < 2; k++) {
    if (a[k] <= -lim || a[k] >= lim || b[k] <= -lim || b[k] >= lim) {
      return GSHGEOM_ERR_INVALID;
    }
  }
  const int64_t det = a[0] * b[1] - a[1] * b[0]; /* each product is below 2^62 */
  if (det == 0) {
    return GSHGEOM_ERR_INVALID;
  }
  const int64_t adet = det < 0 ? -det : det;
  if (adet > GSHGEOM_LATTICE_PERIOD_MAX) {
    return GSHGEOM_ERR_LIMIT;
  }
  /* The extended Euclid algorithm on |a1| and |b1|: u |a1| + v |b1| = g. */
  int64_t r0 = a[1] < 0 ? -a[1] : a[1];
  int64_t r1 = b[1] < 0 ? -b[1] : b[1];
  int64_t s0 = 1;
  int64_t s1 = 0;
  int64_t t0 = 0;
  int64_t t1 = 1;
  while (r1 != 0) {
    const int64_t q = r0 / r1;
    int64_t tmp = r0 - q * r1;
    r0 = r1;
    r1 = tmp;
    tmp = s0 - q * s1;
    s0 = s1;
    s1 = tmp;
    tmp = t0 - q * t1;
    t0 = t1;
    t1 = tmp;
  }
  const int64_t g = r0; /* gcd(a1, b1) >= 1, since det != 0 */
  const int64_t u = a[1] < 0 ? -s0 : s0;
  const int64_t v = b[1] < 0 ? -t0 : t0;
  const int64_t rr = u * a[0] + v * b[0];
  out->s = g;
  out->p = adet / g;
  out->r = floor_mod(rr, out->p);
  return GSHGEOM_OK;
}

void gshgeom_lattice_reduce(const GSHGEOM_Lattice_Period * pd, int64_t x, int64_t y, int64_t out[2]) {
  const int64_t k = floor_div(y, pd->s);
  const int64_t km = floor_mod(k, pd->p);
  out[0] = floor_mod(x - km * pd->r, pd->p);
  out[1] = y - k * pd->s;
}
