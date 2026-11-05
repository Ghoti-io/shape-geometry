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
 * Hatching (hatch.h, HX4 to HX12).
 *
 * **The method.**  For each line of a family the device line is the one through the double point `A_k` with the double
 * direction `w`; the side of every contour point is decided exactly (a sign filter, then an exact sum of products).  The
 * region along a line is read twice, along the line shifted infinitesimally to its left and to its right: a point on the
 * line counts as right of the left-shifted line and left of the right-shifted one, so an edge that lies along the line is
 * crossed by neither, a vertex on it is crossed consistently, and no case is special.  Each shifted line meets the
 * contour edges at isolated crossings; the winding number, and each contour's own parity, change by one at each, and the
 * region's membership is read between consecutive distinct crossing positions.  The trim intervals are the union of the
 * intervals of the two readings, which is `closure(G) cap L` (a segment of `L` in the closure of `G` is on the boundary of a
 * cell of `G`, so an open set of `G` lies beside it on one side at least).  Positions are computed in `double` from the
 * vertices' projections, so a vertex shared by two edges is one position; two crossings within 2^-26 pixel are one.
 *
 * Nothing here uses libm, `long double`, an integer to double conversion of a 64-bit value or a fused multiply-add.
 */

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/shape-geometry/hash.h>
#include <ghoti.io/shape-geometry/hatch.h>
#include <ghoti.io/shape-geometry/orient.h>

#include <ghoti.io/math/basic.h>

#include "../exact/exact.h"
#include "../stroke/stroke_internal.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* --- integers and doubles without a 64-bit conversion ------------------------------------------------------------------ */

static double i2d(int64_t v) {
  const int32_t hi = (int32_t)(v >> 32);
  const uint32_t lo = (uint32_t)((uint64_t)v & 0xffffffffu);
  return (double)hi * 4294967296.0 + (double)lo;
}

/** The floor of a finite double of magnitude below 2^52, as an integer. */
static int64_t d2i_floor(double x) {
  const uint64_t bits = gshgeom_f64_bits(x);
  const int neg = (int)(bits >> 63);
  const int ex = (int)((bits >> 52) & 0x7ffu);
  if (ex < 1023) { /* |x| < 1 */
    return neg && x != 0.0 ? -1 : 0;
  }
  const int e = ex - 1075;
  const uint64_t m = (bits & ((UINT64_C(1) << 52) - 1u)) | (UINT64_C(1) << 52);
  uint64_t v;
  int frac = 0;
  if (e >= 0) {
    v = m << e;
  } else {
    v = m >> (-e);
    frac = (m & ((UINT64_C(1) << (-e)) - 1u)) != 0;
  }
  if (neg) {
    return -(int64_t)v - (frac ? 1 : 0);
  }
  return (int64_t)v;
}

static double dabs(double x) {
  return x < 0.0 ? -x : x;
}

static int is_fin(double x) {
  return x - x == 0.0;
}

/* --- the pieces object ------------------------------------------------------------------------------------------------- */

typedef struct Piece {
  double x0, y0, x1, y1; /* local space */
  double dx, dy;         /* the mapped direction of the line (used by a dot) */
  int dot;
} Piece;

struct GSHGEOM_Hatch_Pieces {
  GSHGEOM_Allocator alloc;
  Piece * v;
  size_t n, cap;
  uint64_t lines;
  int h_singular; /* the hatch map's linear part is exactly singular: nothing is drawn (HX8) */
};

void gshgeom_hatch_pieces_destroy(GSHGEOM_Hatch_Pieces * p) {
  if (p) {
    GSHGEOM_Allocator a = p->alloc;
    gcu_allocator_free(&a, p->v);
    gcu_allocator_free(&a, p);
  }
}

size_t gshgeom_hatch_pieces_count(const GSHGEOM_Hatch_Pieces * p) {
  return p ? p->n : 0;
}

uint64_t gshgeom_hatch_pieces_lines(const GSHGEOM_Hatch_Pieces * p) {
  return p ? p->lines : 0;
}

GSHGEOM_Result gshgeom_hatch_pieces_get(const GSHGEOM_Hatch_Pieces * p, size_t index, double xy[4], int * dot) {
  if (!p || !xy || !dot || index >= p->n) {
    return GSHGEOM_ERR_INVALID;
  }
  xy[0] = p->v[index].x0;
  xy[1] = p->v[index].y0;
  xy[2] = p->v[index].x1;
  xy[3] = p->v[index].y1;
  *dot = p->v[index].dot;
  return GSHGEOM_OK;
}

uint64_t gshgeom_hatch_pieces_hash(const GSHGEOM_Hatch_Pieces * p) {
  GSHGEOM_Hash h;
  gshgeom_hash_init(&h);
  gshgeom_hash_u32(&h, 1);
  gshgeom_hash_u32(&h, UINT32_C(0x43435048));
  gshgeom_hash_u64(&h, p ? (uint64_t)p->n : 0);
  for (size_t i = 0; p && i < p->n; i++) {
    gshgeom_hash_u8(&h, (uint8_t)(p->v[i].dot ? 1 : 0));
    gshgeom_hash_f64(&h, p->v[i].x0);
    gshgeom_hash_f64(&h, p->v[i].y0);
    gshgeom_hash_f64(&h, p->v[i].x1);
    gshgeom_hash_f64(&h, p->v[i].y1);
  }
  return gshgeom_hash_value(&h);
}

GSHGEOM_Result gshgeom_hatch_pieces_dump(const GSHGEOM_Hatch_Pieces * p, FILE * out) {
  if (!p || !out) {
    return GSHGEOM_ERR_INVALID;
  }
  int bad = fprintf(out, "gshgeom.hatch-pieces 1 pieces=%llu lines=%llu hash=%016llx\n", (unsigned long long)p->n,
                    (unsigned long long)p->lines, (unsigned long long)gshgeom_hatch_pieces_hash(p)) < 0;
  for (size_t i = 0; i < p->n && !bad; i++) {
    bad = fprintf(out, "piece %s %016llx %016llx %016llx %016llx\n", p->v[i].dot ? "dot" : "seg",
                  (unsigned long long)gshgeom_f64_bits(p->v[i].x0), (unsigned long long)gshgeom_f64_bits(p->v[i].y0),
                  (unsigned long long)gshgeom_f64_bits(p->v[i].x1), (unsigned long long)gshgeom_f64_bits(p->v[i].y1)) < 0;
  }
  return bad || fflush(out) != 0 ? GSHGEOM_ERR_IO : GSHGEOM_OK;
}

int gshgeom_hatch_revision(void) {
  return GSHGEOM_HATCH_REVISION;
}

/* --- HX2: the geometric checks ----------------------------------------------------------------------------------------- */

static int in_band(double x) {
  const double m = dabs(x);
  return m == 0.0 || (m >= GSHGEOM_HATCH_VALUE_MIN && m <= GSHGEOM_HATCH_VALUE_MAX);
}

GSHGEOM_Result gshgeom_hatch_check(const GSHGEOM_Hatch * h) {
  if (!h) {
    return GSHGEOM_ERR_INVALID;
  }
  if (h->angle_unit >= (uint32_t)GSHGEOM_ANGLE_UNIT_COUNT || h->nesting >= (uint32_t)GSHGEOM_HATCH_NESTING_COUNT ||
      h->family_count < 1 || !h->families) {
    return GSHGEOM_ERR_INVALID;
  }
  for (size_t i = 0; i < h->family_count; i++) {
    const GSHGEOM_Hatch_Family * f = &h->families[i];
    if (!is_fin(f->angle) || !is_fin(f->origin_x) || !is_fin(f->origin_y) || !is_fin(f->offset_along) ||
        !is_fin(f->offset_across) || (f->dash_count && !f->dashes)) {
      return GSHGEOM_ERR_INVALID;
    }
    for (size_t k = 0; k < f->dash_count; k++) {
      if (!is_fin(f->dashes[k])) {
        return GSHGEOM_ERR_INVALID;
      }
    }
  }
  GSHGEOM_Result r = gshgeom_map_check(&h->map);
  if (r != GSHGEOM_OK) {
    return GSHGEOM_ERR_INVALID;
  }
  if (h->map.dim != 2) {
    return GSHGEOM_ERR_UNSUPPORTED;
  }
  /* a family whose elements have no pen-down element or no length draws nothing and has no period: invalid, as
   * ::gshgeom_stroke_dashes_from_dxf says (the bands are the limit stage's) */
  for (size_t i = 0; i < h->family_count; i++) {
    const GSHGEOM_Hatch_Family * f = &h->families[i];
    if (f->dash_count) {
      int down = 0;
      double total = 0.0;
      for (size_t k = 0; k < f->dash_count; k++) {
        down |= f->dashes[k] >= 0.0;
        total += dabs(f->dashes[k]);
      }
      if (!down || !(total > 0.0)) {
        return GSHGEOM_ERR_INVALID;
      }
    }
  }
  if (h->family_count > GSHGEOM_HATCH_FAMILIES_MAX) {
    return GSHGEOM_ERR_LIMIT;
  }
  for (size_t i = 0; i < h->family_count; i++) {
    const GSHGEOM_Hatch_Family * f = &h->families[i];
    if (f->dash_count > GSHGEOM_HATCH_DASHES_MAX || !in_band(f->origin_x) || !in_band(f->origin_y) ||
        !in_band(f->offset_along) || !in_band(f->offset_across) || dabs(f->angle) > GSHGEOM_HATCH_ANGLE_MAX) {
      return GSHGEOM_ERR_LIMIT;
    }
    for (size_t k = 0; k < f->dash_count; k++) {
      if (!in_band(f->dashes[k])) {
        return GSHGEOM_ERR_LIMIT;
      }
    }
    double c, s;
    if (gshgeom_sincos((GSHGEOM_Angle_Unit)h->angle_unit, f->angle, &c, &s) != GSHGEOM_OK) {
      return GSHGEOM_ERR_LIMIT;
    }
  }
  for (int a = 0; a < 2; a++) {
    for (int b = 0; b < 2; b++) {
      if (!in_band(h->map.m[a][b])) {
        return GSHGEOM_ERR_LIMIT;
      }
    }
    if (!in_band(h->map.t[a])) {
      return GSHGEOM_ERR_LIMIT;
    }
  }
  return GSHGEOM_OK;
}

/* --- the exact side of a point ----------------------------------------------------------------------------------------- */

/** The sign of `w x (P - A)`, exact, and an approximate value of it (never 0 for a nonzero sign). */
static int side_of(double wx, double wy, double ax, double ay, double px, double py, double * value) {
  const double t1 = wx * (py - ay);
  const double t2 = wy * (px - ax);
  const double c = t1 - t2;
  const double bound = (dabs(t1) + dabs(t2)) * 0x1p-49;
  int sg;
  if (c > bound) {
    sg = 1;
  } else if (c < -bound) {
    sg = -1;
  } else {
    GSHGEOM_ExactTerm t[4];
    size_t n = 0;
    gshgeom_exact_put(t, &n, 0, 2, wx, py, 0.0);
    gshgeom_exact_put(t, &n, 1, 2, wx, ay, 0.0);
    gshgeom_exact_put(t, &n, 1, 2, wy, px, 0.0);
    gshgeom_exact_put(t, &n, 0, 2, wy, ax, 0.0);
    sg = gshgeom_exact_sign(t, n);
  }
  if (sg != 0 && (sg > 0 ? !(c > 0.0) : !(c < 0.0))) {
    *value = sg > 0 ? 0x1p-1000 : -0x1p-1000;
  } else {
    *value = c;
  }
  return sg;
}

/* --- one region: the contours ------------------------------------------------------------------------------------------ */

typedef struct Region {
  size_t np, nc;
  const GSHGEOM_Device_Point * pt;
  GSHGEOM_Flat_Subpath * sub; /* nc subpaths */
  double * px;                /* pixels, exact */
  double * py;
  int32_t lo_x, lo_y, hi_x, hi_y; /* the closed bounding box, 26.6 */
} Region;

typedef struct Ev {
  double s;
  uint32_t contour;
  int8_t dw;
  uint32_t seq;
} Ev;

typedef struct Iv {
  double t0, t1;
} Iv;

typedef struct Work {
  GSHGEOM_Allocator * alloc;
  int8_t * sg;     /* np */
  double * sp;     /* np: the parameter of the point's projection */
  double * cv;     /* np: approximate cross value */
  Ev * ev;         /* np (one crossing per edge at most) */
  Iv * iv;         /* 2 * np / 2 + 2 */
  Iv * out;        /* the merged intervals */
  uint8_t * par;   /* nc */
  size_t iv_cap;
  size_t out_n;
} Work;

static int ev_cmp(const void * a, const void * b) {
  const Ev * x = (const Ev *)a;
  const Ev * y = (const Ev *)b;
  if (x->s < y->s) {
    return -1;
  }
  if (x->s > y->s) {
    return 1;
  }
  return x->seq < y->seq ? -1 : x->seq > y->seq;
}

static int iv_cmp(const void * a, const void * b) {
  const Iv * x = (const Iv *)a;
  const Iv * y = (const Iv *)b;
  return x->t0 < y->t0 ? -1 : x->t0 > y->t0;
}

static int member(uint32_t nesting, GSHGEOM_Fill_Rule rule, int64_t wind, int64_t odd_count) {
  switch (nesting) {
    case GSHGEOM_HATCH_NESTING_OUTER:
      return odd_count == 1;
    case GSHGEOM_HATCH_NESTING_IGNORE:
      return odd_count >= 1;
    default:
      return rule == GSHGEOM_FILL_NONZERO ? wind != 0 : (wind % 2) != 0;
  }
}

/** Read the region along one shifted line (@p up: left-shifted) and append its intervals to w->iv at *n. */
static void read_line(const Region * R, Work * w, uint32_t nesting, GSHGEOM_Fill_Rule rule, int up, double band_s, size_t * n) {
  size_t ne = 0;
  uint32_t seq = 0;
  for (size_t c = 0; c < R->nc; c++) {
    const GSHGEOM_Flat_Subpath sp = R->sub[c];
    for (size_t i = 0; i < sp.count; i++) {
      const size_t a = sp.first + i;
      const size_t b = sp.first + (i + 1 == sp.count ? 0 : i + 1);
      const int sa = w->sg[a];
      const int sb = w->sg[b];
      const int ka = up ? (sa > 0 ? 1 : -1) : (sa < 0 ? -1 : 1);
      const int kb = up ? (sb > 0 ? 1 : -1) : (sb < 0 ? -1 : 1);
      if (ka == kb) {
        continue;
      }
      double s;
      if (sa == 0) {
        s = w->sp[a];
      } else if (sb == 0) {
        s = w->sp[b];
      } else {
        const double ca = w->cv[a];
        const double cb = w->cv[b];
        s = w->sp[a] + (ca / (ca - cb)) * (w->sp[b] - w->sp[a]);
      }
      w->ev[ne].s = s;
      w->ev[ne].contour = (uint32_t)c;
      w->ev[ne].dw = (int8_t)(ka < kb ? -1 : 1);
      w->ev[ne].seq = seq++;
      ne++;
    }
  }
  if (ne == 0) {
    return;
  }
  qsort(w->ev, ne, sizeof w->ev[0], ev_cmp);
  /* w->par is all zero here: a closed contour crosses the shifted line an even number of times (the signs around a cycle change
   * an even number of times), so every parity flag is back at zero when a reading ends. */
  int64_t wind = 0, odd = 0;
  size_t i = 0;
  int in = 0;
  double start = 0.0;
  while (i < ne) {
    size_t j = i;
    /* a group: crossings within the band of their predecessor are one position */
    for (;;) {
      wind += w->ev[j].dw;
      const uint32_t c = w->ev[j].contour;
      w->par[c] ^= 1u;
      odd += w->par[c] ? 1 : -1;
      if (j + 1 < ne && w->ev[j + 1].s - w->ev[j].s <= band_s) {
        j++;
      } else {
        break;
      }
    }
    const double first = w->ev[i].s;
    const double last = w->ev[j].s;
    const int m = member(nesting, rule, wind, odd);
    if (in && !m) {
      /* The interval ends at the first position of the group that leaves the region.  It is longer than the band without a test:
       * a group is closed under "within the band of the last crossing", so this group's first position is more than a band past
       * the last position of the group that entered, where the interval began. */
      w->iv[(*n)].t0 = start;
      w->iv[(*n)].t1 = first;
      (*n)++;
      in = 0;
    } else if (!in && m) {
      /* it begins at the last position of the group that enters the region; a group that leaves the state as it was changes nothing,
       * so a run of crossings each within the band of the last never cuts the interval it lies in */
      in = 1;
      start = last;
    }
    i = j + 1;
  }
}

/** Merge the n intervals of w->iv into w->out (sorted, touching or overlapping ones joined). */
static void merge_intervals(Work * w, size_t n, double band_s) {
  w->out_n = 0;
  if (n == 0) {
    return;
  }
  qsort(w->iv, n, sizeof w->iv[0], iv_cmp);
  Iv cur = w->iv[0];
  for (size_t i = 1; i < n; i++) {
    if (w->iv[i].t0 <= cur.t1 + band_s) {
      if (w->iv[i].t1 > cur.t1) {
        cur.t1 = w->iv[i].t1;
      }
    } else {
      w->out[w->out_n++] = cur;
      cur = w->iv[i];
    }
  }
  w->out[w->out_n++] = cur;
}

/* --- pieces -------------------------------------------------------------------------------------------------------------- */

typedef struct Emit {
  GSHGEOM_Hatch_Pieces * P;
  size_t max_pieces;
  const GSHGEOM_Map * H;
  double ox, oy, ux, uy; /* the line's base point and direction, in hatch space */
  double dirx, diry;     /* H_lin u */
  const Iv * iv;
  size_t ivn, j;
  double base, slack, band; /* band: the shortest length that counts, in s (HX4) */
  GSHGEOM_Result err;
} Emit;

static GSHGEOM_Result push(Emit * e, double s0, double s1, int dot) {
  GSHGEOM_Hatch_Pieces * P = e->P;
  if (P->n >= e->max_pieces) {
    return GSHGEOM_ERR_LIMIT;
  }
  if (P->n == P->cap) {
    size_t cap = P->cap ? 2 * P->cap : 64;
    if (cap > e->max_pieces) {
      cap = e->max_pieces;
    }
    Piece * bigger = gcu_allocator_realloc(&P->alloc, P->v, cap * sizeof *bigger);
    if (!bigger) {
      return GSHGEOM_ERR_OOM;
    }
    P->v = bigger;
    P->cap = cap;
  }
  GSHGEOM_Point a = {e->ox + s0 * e->ux, e->oy + s0 * e->uy, 0.0};
  GSHGEOM_Point b = {e->ox + s1 * e->ux, e->oy + s1 * e->uy, 0.0};
  GSHGEOM_Point la, lb;
  if (gshgeom_map_apply_point(e->H, a, &la) != GSHGEOM_OK || gshgeom_map_apply_point(e->H, b, &lb) != GSHGEOM_OK ||
      !is_fin(la.x) || !is_fin(la.y) || !is_fin(lb.x) || !is_fin(lb.y)) {
    return GSHGEOM_ERR_LIMIT;
  }
  Piece * p = &P->v[P->n++];
  p->x0 = la.x;
  p->y0 = la.y;
  p->x1 = lb.x; /* a dot is pushed with s0 == s1, so lb == la */
  p->y1 = lb.y;
  p->dx = e->dirx;
  p->dy = e->diry;
  p->dot = dot;
  return GSHGEOM_OK;
}

static GSHGEOM_Result dash_cb(void * ctx, double a, double b) {
  Emit * e = (Emit *)ctx;
  const double s0 = e->base + a;
  const double s1 = e->base + b;
  /* intervals wholly before this dash can never meet a later one */
  while (e->j < e->ivn && e->iv[e->j].t1 + e->slack < s0) {
    e->j++;
  }
  if (b == a) {
    for (size_t j = e->j; j < e->ivn && e->iv[j].t0 - e->slack <= s0; j++) {
      if (s0 <= e->iv[j].t1 + e->slack) {
        double d = s0 < e->iv[j].t0 ? e->iv[j].t0 : (s0 > e->iv[j].t1 ? e->iv[j].t1 : s0);
        return push(e, d, d, 1);
      }
    }
    return GSHGEOM_OK;
  }
  for (size_t j = e->j; j < e->ivn && e->iv[j].t0 < s1; j++) {
    const double lo = s0 > e->iv[j].t0 ? s0 : e->iv[j].t0;
    const double hi = s1 < e->iv[j].t1 ? s1 : e->iv[j].t1;
    if (hi - lo > e->band) {
      GSHGEOM_Result r = push(e, lo, hi, 0);
      if (r != GSHGEOM_OK) {
        return r;
      }
    }
  }
  return GSHGEOM_OK;
}

/* --- the driver ---------------------------------------------------------------------------------------------------------- */

static void work_free(GSHGEOM_Allocator * a, Work * w, Region * R) {
  gcu_allocator_free(a, w->sg);
  gcu_allocator_free(a, w->sp);
  gcu_allocator_free(a, w->cv);
  gcu_allocator_free(a, w->ev);
  gcu_allocator_free(a, w->iv);
  gcu_allocator_free(a, w->out);
  gcu_allocator_free(a, w->par);
  gcu_allocator_free(a, R->px);
  gcu_allocator_free(a, R->py);
  gcu_allocator_free(a, R->sub);
}

static GSHGEOM_Result pieces_core(const GSHGEOM_Flat_Path * flat, GSHGEOM_Fill_Rule rule, const GSHGEOM_Hatch * hatch,
    const GSHGEOM_Map * chain, uint64_t max_lines, const GSHGEOM_Allocator * alloc, const GSHGEOM_Limits * limits_in,
    GSHGEOM_Hatch_Pieces ** out) {
  if (!flat || !hatch || !chain || !out || (unsigned)rule >= (unsigned)GSHGEOM_FILL_RULE_COUNT) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Result r = gshgeom_hatch_check(hatch);
  if (r != GSHGEOM_OK) {
    return r;
  }
  r = gshgeom_map_check(chain);
  if (r != GSHGEOM_OK || chain->dim != 2) {
    return r == GSHGEOM_OK ? GSHGEOM_ERR_UNSUPPORTED : GSHGEOM_ERR_INVALID;
  }
  if (!alloc) {
    alloc = gshgeom_allocator_default();
  }
  GSHGEOM_Limits limits;
  gshgeom_limits_default(&limits);
  if (limits_in) {
    if (limits_in->max_elements) {
      limits.max_elements = limits_in->max_elements;
    }
    if (limits_in->max_points) {
      limits.max_points = limits_in->max_points;
    }
    limits.max_dash_steps = limits_in->max_dash_steps;
  }
  if (max_lines == 0) {
    max_lines = GSHGEOM_HATCH_LINES_DEFAULT;
  }
  GSHGEOM_Hatch_Pieces * P = gcu_allocator_calloc(alloc, 1, sizeof *P);
  if (!P) {
    return GSHGEOM_ERR_OOM;
  }
  P->alloc = *alloc;
  int inv_chain = 0, inv_h = 0;
  r = gshgeom_stroke_map_invertible(chain, &inv_chain);
  if (r == GSHGEOM_OK) {
    r = gshgeom_stroke_map_invertible(&hatch->map, &inv_h);
  }
  if (r != GSHGEOM_OK) {
    gshgeom_hatch_pieces_destroy(P);
    return r;
  }
  P->h_singular = !inv_h;
  GSHGEOM_Map C;
  r = gshgeom_map_compose(chain, &hatch->map, &C);
  const size_t np = gshgeom_flat_path_point_count(flat);
  const size_t nc = gshgeom_flat_path_subpath_count(flat);
  if (r != GSHGEOM_OK || !inv_chain || !inv_h || np == 0 || nc == 0) {
    if (r != GSHGEOM_OK && inv_chain && inv_h) {
      gshgeom_hatch_pieces_destroy(P);
      return GSHGEOM_ERR_LIMIT; /* the composition overflowed */
    }
    *out = P; /* nothing is drawn, no lines */
    return GSHGEOM_OK;
  }
  Region R;
  memset(&R, 0, sizeof R);
  Work w;
  memset(&w, 0, sizeof w);
  GSHGEOM_Allocator a = *alloc;
  w.alloc = &a;
  R.np = np;
  R.nc = nc;
  R.pt = gshgeom_flat_path_points(flat);
  R.sub = gcu_allocator_calloc(&a, nc, sizeof *R.sub);
  R.px = gcu_allocator_calloc(&a, np, sizeof(double));
  R.py = gcu_allocator_calloc(&a, np, sizeof(double));
  w.sg = gcu_allocator_calloc(&a, np, sizeof(int8_t));
  w.sp = gcu_allocator_calloc(&a, np, sizeof(double));
  w.cv = gcu_allocator_calloc(&a, np, sizeof(double));
  w.ev = gcu_allocator_calloc(&a, np, sizeof(Ev));
  w.iv_cap = np + 2;
  w.iv = gcu_allocator_calloc(&a, w.iv_cap, sizeof(Iv));
  w.out = gcu_allocator_calloc(&a, w.iv_cap, sizeof(Iv));
  w.par = gcu_allocator_calloc(&a, nc, 1);
  double * dashes = NULL;
  if (!R.sub || !R.px || !R.py || !w.sg || !w.sp || !w.cv || !w.ev || !w.iv || !w.out || !w.par) {
    r = GSHGEOM_ERR_OOM;
    goto done;
  }
  for (size_t c = 0; c < nc; c++) {
    (void)gshgeom_flat_path_subpath(flat, c, &R.sub[c]);
  }
  for (size_t i = 0; i < np; i++) {
    R.px[i] = (double)R.pt[i].x / 64.0;
    R.py[i] = (double)R.pt[i].y / 64.0;
    if (i == 0 || R.pt[i].x < R.lo_x) {
      R.lo_x = R.pt[i].x;
    }
    if (i == 0 || R.pt[i].y < R.lo_y) {
      R.lo_y = R.pt[i].y;
    }
    if (i == 0 || R.pt[i].x > R.hi_x) {
      R.hi_x = R.pt[i].x;
    }
    if (i == 0 || R.pt[i].y > R.hi_y) {
      R.hi_y = R.pt[i].y;
    }
  }
  const double cx[4] = {(double)R.lo_x / 64.0, (double)R.hi_x / 64.0, (double)R.hi_x / 64.0, (double)R.lo_x / 64.0};
  const double cy[4] = {(double)R.lo_y / 64.0, (double)R.lo_y / 64.0, (double)R.hi_y / 64.0, (double)R.hi_y / 64.0};
  const size_t max_pieces = limits.max_elements / 2 ? limits.max_elements / 2 : 1;
  const size_t dash_steps = limits.max_dash_steps && limits.max_dash_steps < (uint64_t)SIZE_MAX ? (size_t)limits.max_dash_steps
                                                                                                  : 4 * limits.max_points + 4096;
  Emit e;
  memset(&e, 0, sizeof e);
  e.P = P;
  e.max_pieces = max_pieces;
  e.H = &hatch->map;
  for (size_t fi = 0; fi < hatch->family_count && r == GSHGEOM_OK; fi++) {
    const GSHGEOM_Hatch_Family * f = &hatch->families[fi];
    if (f->offset_across == 0.0) {
      continue; /* HX9: infinitely many coincident lines draw nothing */
    }
    double ux, uy;
    {
      double cs, sn;
      if (gshgeom_sincos((GSHGEOM_Angle_Unit)hatch->angle_unit, f->angle, &cs, &sn) != GSHGEOM_OK) {
        r = GSHGEOM_ERR_LIMIT;
        break;
      }
      ux = cs;
      uy = sn;
    }
    const double nx = -uy, ny = ux;
    /* the lowered dash array of the family */
    gcu_allocator_free(&a, dashes);
    dashes = NULL;
    size_t nd = 0;
    double off = 0.0, period = 0.0;
    if (f->dash_count) {
      dashes = gcu_allocator_calloc(&a, f->dash_count + 2, sizeof(double));
      if (!dashes) {
        r = GSHGEOM_ERR_OOM;
        break;
      }
      r = gshgeom_stroke_dashes_from_dxf(f->dashes, f->dash_count, 1.0, dashes, f->dash_count + 2, &nd, &off);
      if (r != GSHGEOM_OK) {
        break;
      }
      for (size_t k = 0; k < nd; k++) {
        period += dashes[k];
      }
    }
    GSHGEOM_Point wv, wn, dv;
    GSHGEOM_Point uv = {ux, uy, 0.0}, nv = {nx, ny, 0.0};
    if (gshgeom_map_apply_vector(&C, uv, &wv) != GSHGEOM_OK || gshgeom_map_apply_vector(&C, nv, &wn) != GSHGEOM_OK ||
        gshgeom_map_apply_vector(&hatch->map, uv, &dv) != GSHGEOM_OK) {
      r = GSHGEOM_ERR_LIMIT;
      break;
    }
    const double wx = wv.x, wy = wv.y;
    const double ww = wx * wx + wy * wy;
    if (!is_fin(ww) || !(ww > 0.0)) {
      r = GSHGEOM_ERR_LIMIT;
      break;
    }
    const double band_s = GSHGEOM_HATCH_BAND / gmath_sqrt(ww);
    e.ux = ux;
    e.uy = uy;
    e.dirx = dv.x;
    e.diry = dv.y;
    e.slack = 0.0;
    e.band = band_s;
    /* the transversal measure: line k lies at k * kappa from line 0 */
    GSHGEOM_Point o0 = {f->origin_x, f->origin_y, 0.0}, a0;
    if (gshgeom_map_apply_point(&C, o0, &a0) != GSHGEOM_OK) {
      r = GSHGEOM_ERR_LIMIT;
      break;
    }
    const double kappa = f->offset_across * (wx * wn.y - wy * wn.x);
    if (!is_fin(kappa) || kappa == 0.0) {
      r = GSHGEOM_ERR_LIMIT;
      break;
    }
    double klo = 0.0, khi = 0.0;
    for (int k = 0; k < 4; k++) {
      const double tau = wx * (cy[k] - a0.y) - wy * (cx[k] - a0.x);
      const double kk = tau / kappa;
      if (!is_fin(kk) || dabs(kk) > 0x1p50) {
        r = GSHGEOM_ERR_LIMIT;
        break;
      }
      if (k == 0 || kk < klo) {
        klo = kk;
      }
      if (k == 0 || kk > khi) {
        khi = kk;
      }
    }
    if (r != GSHGEOM_OK) {
      break;
    }
    const int64_t k_lo = d2i_floor(klo) - 2;
    const int64_t k_hi = d2i_floor(khi) + 2;
    if ((uint64_t)(k_hi - k_lo + 1) > (max_lines - (P->lines < max_lines ? P->lines : max_lines)) + 6u) {
      r = GSHGEOM_ERR_LIMIT;
      break;
    }
    for (int64_t k = k_lo; k <= k_hi && r == GSHGEOM_OK; k++) {
      const double kd = i2d(k);
      const double kdx = kd * f->offset_along;
      const double kdy = kd * f->offset_across;
      const double ox = (f->origin_x + kdx * ux) + kdy * nx;
      const double oy = (f->origin_y + kdx * uy) + kdy * ny;
      GSHGEOM_Point ok = {ox, oy, 0.0}, ak;
      if (gshgeom_map_apply_point(&C, ok, &ak) != GSHGEOM_OK || !is_fin(ak.x) || !is_fin(ak.y)) {
        r = GSHGEOM_ERR_LIMIT;
        break;
      }
      /* does the line meet the closed bounding box? */
      int pos = 0, neg = 0;
      for (int q = 0; q < 4; q++) {
        double v;
        const int sg = side_of(wx, wy, ak.x, ak.y, cx[q], cy[q], &v);
        pos += sg > 0;
        neg += sg < 0;
      }
      if (pos == 4 || neg == 4) {
        continue;
      }
      if (P->lines >= max_lines) {
        r = GSHGEOM_ERR_LIMIT;
        break;
      }
      P->lines++;
      for (size_t i = 0; i < np; i++) {
        double v;
        w.sg[i] = (int8_t)side_of(wx, wy, ak.x, ak.y, R.px[i], R.py[i], &v);
        w.cv[i] = v;
        w.sp[i] = ((R.px[i] - ak.x) * wx + (R.py[i] - ak.y) * wy) / ww;
      }
      size_t n = 0;
      read_line(&R, &w, hatch->nesting, rule, 1, band_s, &n);
      read_line(&R, &w, hatch->nesting, rule, 0, band_s, &n);
      merge_intervals(&w, n, band_s);
      if (w.out_n == 0) {
        continue;
      }
      e.ox = ox;
      e.oy = oy;
      e.iv = w.out;
      e.ivn = w.out_n;
      e.j = 0;
      if (nd == 0) {
        e.slack = 0.0;
        for (size_t q = 0; q < w.out_n && r == GSHGEOM_OK; q++) {
          r = push(&e, w.out[q].t0, w.out[q].t1, 0);
        }
      } else {
        const double t_first = w.out[0].t0;
        const double t_last = w.out[w.out_n - 1].t1;
        e.base = t_first - period;
        {
          const double sc = dabs(t_first) > dabs(t_last) ? dabs(t_first) : dabs(t_last);
          e.slack = 0x1p-46 * (period > sc ? period : (sc > 1.0 ? sc : 1.0));
        }
        r = gshgeom_stroke_dash_each(dashes, nd, off + e.base, (t_last - t_first) + 2.0 * period, dash_steps, dash_cb, &e);
      }
    }
  }
done:
  gcu_allocator_free(&a, dashes);
  work_free(&a, &w, &R);
  if (r != GSHGEOM_OK) {
    gshgeom_hatch_pieces_destroy(P);
    return r;
  }
  *out = P;
  return GSHGEOM_OK;
}

/* --- the public entries -------------------------------------------------------------------------------------------------- */

GSHGEOM_Result gshgeom_hatch_pieces_from_flat(const GSHGEOM_Flat_Path * region, GSHGEOM_Fill_Rule rule, const GSHGEOM_Hatch * hatch,
    const GSHGEOM_Map * chain, uint64_t max_lines, const GSHGEOM_Allocator * alloc, const GSHGEOM_Limits * limits,
    GSHGEOM_Hatch_Pieces ** out) {
  return pieces_core(region, rule, hatch, chain, max_lines, alloc, limits, out);
}

GSHGEOM_Result gshgeom_hatch_pieces_with_tolerance(const GSHGEOM_Path * region, GSHGEOM_Fill_Rule rule, const GSHGEOM_Hatch * hatch,
    const GSHGEOM_Map * chain, int32_t tolerance, uint64_t max_lines, const GSHGEOM_Allocator * alloc,
    const GSHGEOM_Limits * limits, GSHGEOM_Hatch_Pieces ** out) {
  if (!region || !chain || !out || tolerance < GSHGEOM_FLATTEN_MIN_TOLERANCE) {
    return GSHGEOM_ERR_INVALID;
  }
  if (gshgeom_path_dimension(region) != 2) {
    return GSHGEOM_ERR_UNSUPPORTED;
  }
  GSHGEOM_Result r = gshgeom_map_check(chain);
  if (r != GSHGEOM_OK || chain->dim != 2) {
    return r == GSHGEOM_OK ? GSHGEOM_ERR_UNSUPPORTED : GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Device_Path * dp = NULL;
  GSHGEOM_Flat_Path * flat = NULL;
  r = gshgeom_device_path_from_model_with_tolerance(region, chain, tolerance, alloc, limits, &dp);
  if (r == GSHGEOM_OK) {
    r = gshgeom_device_path_flatten_with_tolerance(dp, tolerance, alloc, limits, &flat);
  }
  gshgeom_device_path_destroy(dp);
  if (r == GSHGEOM_OK) {
    r = pieces_core(flat, rule, hatch, chain, max_lines, alloc, limits, out);
  }
  gshgeom_flat_path_destroy(flat);
  return r;
}

GSHGEOM_Result gshgeom_hatch_pieces(const GSHGEOM_Path * region, GSHGEOM_Fill_Rule rule, const GSHGEOM_Hatch * hatch,
    const GSHGEOM_Map * chain, uint64_t max_lines, const GSHGEOM_Allocator * alloc, const GSHGEOM_Limits * limits,
    GSHGEOM_Hatch_Pieces ** out) {
  return gshgeom_hatch_pieces_with_tolerance(region, rule, hatch, chain, GSHGEOM_FLATTEN_TOLERANCE, max_lines, alloc, limits, out);
}

/* --- the stroke of the pieces (HX7) ------------------------------------------------------------------------------------- */

GSHGEOM_Result gshgeom_hatch_expand_pieces(const GSHGEOM_Hatch_Pieces * pieces, const GSHGEOM_Stroke * stroke,
    const GSHGEOM_Map * chain, int32_t tolerance, const GSHGEOM_Allocator * alloc, const GSHGEOM_Limits * limits,
    GSHGEOM_Stroke_Outline ** out) {
  if (!pieces || !stroke || !chain || !out) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Result r = gshgeom_stroke_check(stroke);
  if (r != GSHGEOM_OK) {
    return r;
  }
  if (stroke->dash_count != 0 || stroke->alignment != (uint32_t)GSHGEOM_STROKE_ALIGN_CENTER) {
    return GSHGEOM_ERR_INVALID;
  }
  if (!alloc) {
    alloc = gshgeom_allocator_default();
  }
  GSHGEOM_Stroke st = *stroke;
  GSHGEOM_Map map = *chain;
  if (pieces->h_singular) {
    /* the hatch map draws nothing: SG15's status, by a map that is singular */
    r = gshgeom_map_scale(2, 0.0, 0.0, 0.0, &map);
    if (r != GSHGEOM_OK) {
      return r;
    }
  }
  GSHGEOM_Path * path = NULL;
  r = gshgeom_path_create(alloc, limits, 2, &path);
  if (r != GSHGEOM_OK) {
    return r;
  }
  size_t dots = 0;
  for (size_t i = 0; i < pieces->n; i++) {
    dots += pieces->v[i].dot || (pieces->v[i].x0 == pieces->v[i].x1 && pieces->v[i].y0 == pieces->v[i].y1);
  }
  double * dirs = dots ? gcu_allocator_calloc(alloc, dots, 2 * sizeof(double)) : NULL;
  if (dots && !dirs) {
    gshgeom_path_destroy(path);
    return GSHGEOM_ERR_OOM;
  }
  size_t di = 0;
  for (size_t i = 0; i < pieces->n && r == GSHGEOM_OK; i++) {
    const Piece * p = &pieces->v[i];
    const int dot = p->dot || (p->x0 == p->x1 && p->y0 == p->y1);
    GSHGEOM_Point from = {p->x0, p->y0, 0.0};
    GSHGEOM_Point to = {dot ? p->x0 : p->x1, dot ? p->y0 : p->y1, 0.0};
    r = gshgeom_path_move_to(path, from);
    if (r == GSHGEOM_OK) {
      r = gshgeom_path_line_to(path, to);
    }
    if (dot) {
      dirs[2 * di] = p->dx;
      dirs[2 * di + 1] = p->dy;
      di++;
    }
  }
  if (r == GSHGEOM_OK) {
    r = gshgeom_stroke_expand_dots(path, &st, &map, tolerance, alloc, limits, dirs, dots, out);
  }
  gcu_allocator_free(alloc, dirs);
  gshgeom_path_destroy(path);
  return r;
}

static GSHGEOM_Result expand_check(const GSHGEOM_Stroke * stroke) {
  if (!stroke) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Result r = gshgeom_stroke_check(stroke);
  if (r != GSHGEOM_OK) {
    return r;
  }
  return stroke->dash_count != 0 || stroke->alignment != (uint32_t)GSHGEOM_STROKE_ALIGN_CENTER ? GSHGEOM_ERR_INVALID : GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_hatch_expand_from_flat(const GSHGEOM_Flat_Path * region, GSHGEOM_Fill_Rule rule, const GSHGEOM_Hatch * hatch,
    const GSHGEOM_Stroke * stroke, const GSHGEOM_Map * chain, uint64_t max_lines, const GSHGEOM_Allocator * alloc,
    const GSHGEOM_Limits * limits, GSHGEOM_Stroke_Outline ** out) {
  GSHGEOM_Result r = expand_check(stroke);
  if (r != GSHGEOM_OK) {
    return r;
  }
  GSHGEOM_Hatch_Pieces * pieces = NULL;
  if (stroke->width == 0.0) {
    /* SG14: width 0 is the status EMPTY, decided before any line is made (and after the arguments are checked) */
    r = gshgeom_hatch_check(hatch);
    if (r != GSHGEOM_OK) {
      return r;
    }
    pieces = gcu_allocator_calloc(alloc ? alloc : gshgeom_allocator_default(), 1, sizeof *pieces);
    if (!pieces) {
      return GSHGEOM_ERR_OOM;
    }
    pieces->alloc = *(alloc ? alloc : gshgeom_allocator_default());
  } else {
    r = pieces_core(region, rule, hatch, chain, max_lines, alloc, limits, &pieces);
    if (r != GSHGEOM_OK) {
      return r;
    }
  }
  r = gshgeom_hatch_expand_pieces(pieces, stroke, chain, GSHGEOM_FLATTEN_TOLERANCE, alloc, limits, out);
  gshgeom_hatch_pieces_destroy(pieces);
  return r;
}

GSHGEOM_Result gshgeom_hatch_expand_with_tolerance(const GSHGEOM_Path * region, GSHGEOM_Fill_Rule rule, const GSHGEOM_Hatch * hatch,
    const GSHGEOM_Stroke * stroke, const GSHGEOM_Map * chain, int32_t tolerance, uint64_t max_lines,
    const GSHGEOM_Allocator * alloc, const GSHGEOM_Limits * limits, GSHGEOM_Stroke_Outline ** out) {
  GSHGEOM_Result r = expand_check(stroke);
  if (r != GSHGEOM_OK) {
    return r;
  }
  GSHGEOM_Hatch_Pieces * pieces = NULL;
  if (stroke->width == 0.0) {
    r = gshgeom_hatch_check(hatch);
    if (r != GSHGEOM_OK) {
      return r;
    }
    pieces = gcu_allocator_calloc(alloc ? alloc : gshgeom_allocator_default(), 1, sizeof *pieces);
    if (!pieces) {
      return GSHGEOM_ERR_OOM;
    }
    pieces->alloc = *(alloc ? alloc : gshgeom_allocator_default());
  } else {
    r = gshgeom_hatch_pieces_with_tolerance(region, rule, hatch, chain, tolerance, max_lines, alloc, limits, &pieces);
    if (r != GSHGEOM_OK) {
      return r;
    }
  }
  r = gshgeom_hatch_expand_pieces(pieces, stroke, chain, tolerance, alloc, limits, out);
  gshgeom_hatch_pieces_destroy(pieces);
  return r;
}

GSHGEOM_Result gshgeom_hatch_expand(const GSHGEOM_Path * region, GSHGEOM_Fill_Rule rule, const GSHGEOM_Hatch * hatch,
    const GSHGEOM_Stroke * stroke, const GSHGEOM_Map * chain, uint64_t max_lines, const GSHGEOM_Allocator * alloc,
    const GSHGEOM_Limits * limits, GSHGEOM_Stroke_Outline ** out) {
  return gshgeom_hatch_expand_with_tolerance(region, rule, hatch, stroke, chain, GSHGEOM_FLATTEN_TOLERANCE, max_lines, alloc,
                                             limits, out);
}
