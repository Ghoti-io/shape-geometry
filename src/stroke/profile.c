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
 * A stroke with a width profile (stroke.h, SW1 to SW9): the centerline is flattened under the map by the integer flattener, each
 * of its pieces is cut where the profile has a stop, and the region of a piece is the convex hull of the pen at each end, the
 * pen being the disc of the piece's radius there (an ellipse under the map's linear part).  The union of the hulls, wound the one
 * way so that nonzero never cancels them, is the body.
 *
 * Positions are measured in model space: a device vector `d` of the flattened centerline is `M^-1 d` there, `M` the linear part
 * of the map (SW3).  A jump of the profile (two stops at one position) gives the pieces on each side their own limit and, at the
 * position itself, a zero-length piece with both, so the region has both discs.
 */

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/math/basic.h>
#include <ghoti.io/math/elementary.h>
#include <ghoti.io/shape-geometry/device.h>

#include "../core/grow.h"
#include "../flatten/flatten_internal.h"
#include "stroke_internal.h"

#include <stdlib.h>
#include <string.h>

/** Side counts a pen polygon may have, in order: the first that is within the error is used. */
static const unsigned kSides[] = {8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512};

/** A point of a hull's construction. */
typedef struct Pt {
  int64_t x, y;
} Pt;

/** One pen: the offsets of its polygon's vertices, in 26.6, from its center. */
typedef struct Pen {
  int32_t * x;
  int32_t * y;
  unsigned n;
} Pen;

/** The state of one expansion. */
typedef struct Prof {
  GSHGEOM_Allocator alloc;
  GSHGEOM_Limits limits;
  const GSHGEOM_Width_Profile * wp;
  double width;       /* the stroke's width */
  int local;          /* 1: the pen is the image of a disc under the linear part; 0: a disc of device pixels */
  double lin[2][2];   /* the map's linear part, pixels per model unit */
  double inv[2][2];   /* its inverse */
  double err;         /* the most a pen polygon may lie inside its disc, in pixels */
  double cs[512], sn[512]; /* cos and sin of the vertex angles of the largest polygon, and the others are every k-th */
  unsigned sides;     /* the side count in use: fixed for the stroke, by its largest radius */
  Pen a, b;           /* the pens at the two ends of the piece being made */
  int32_t * hx;       /* the hull's points */
  int32_t * hy;
  Pt * pts;           /* the points of the hull construction, and its chain */
  Pt * chain;
  GSHGEOM_Flat_Path * body;
  uint64_t effort;
} Prof;

static double factor_at(const GSHGEOM_Width_Profile * wp, double s, int side) {
  const size_t n = wp->count;
  const GSHGEOM_Width_Stop * st = wp->stops;
  if (s < st[0].position || (s == st[0].position && side < 0)) {
    return st[0].factor;
  }
  if (s > st[n - 1].position || (s == st[n - 1].position && side > 0)) {
    return st[n - 1].factor;
  }
  size_t lo = 0;
  size_t hi = n - 1;
  if (side > 0) {
    /* the last stop with position <= s; then the next one is past s, so it is not equal in position */
    while (lo < hi) {
      size_t mid = lo + (hi - lo + 1) / 2;
      if (st[mid].position <= s) {
        lo = mid;
      } else {
        hi = mid - 1;
      }
    }
    if (lo + 1 >= n) {
      return st[lo].factor;
    }
    const double a = st[lo].position, b = st[lo + 1].position;
    return st[lo].factor + (st[lo + 1].factor - st[lo].factor) * ((s - a) / (b - a));
  }
  /* the first stop with position >= s; the one before it is before s */
  while (lo < hi) {
    size_t mid = lo + (hi - lo) / 2;
    if (st[mid].position >= s) {
      hi = mid;
    } else {
      lo = mid + 1;
    }
  }
  if (lo == 0) {
    return st[0].factor;
  }
  const double a = st[lo - 1].position, b = st[lo].position;
  return st[lo - 1].factor + (st[lo].factor - st[lo - 1].factor) * ((s - a) / (b - a));
}

/** Fill @p pen with the polygon of the disc of @p radius (model units, or pixels under device width), at the stroke's side count. */
static GSHGEOM_Result make_pen(const Prof * p, Pen * pen, double radius) {
  const unsigned n = p->sides;
  const unsigned step = 512u / n;
  pen->n = n;
  for (unsigned k = 0; k < n; k++) {
    const double c = p->cs[k * step], s = p->sn[k * step];
    double dx, dy;
    if (p->local) {
      dx = (p->lin[0][0] * c + p->lin[0][1] * s) * radius;
      dy = (p->lin[1][0] * c + p->lin[1][1] * s) * radius;
    } else {
      dx = c * radius;
      dy = s * radius;
    }
    if (gshgeom_round_26_6(dx, &pen->x[k]) != GSHGEOM_OK || gshgeom_round_26_6(dy, &pen->y[k]) != GSHGEOM_OK) {
      return GSHGEOM_ERR_LIMIT;
    }
  }
  return GSHGEOM_OK;
}

static int64_t cross(int64_t ax, int64_t ay, int64_t bx, int64_t by, int64_t cx, int64_t cy) {
  return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
}

static int pt_less(const void * a, const void * b) {
  const Pt * p = a;
  const Pt * q = b;
  return p->x < q->x ? -1 : p->x > q->x ? 1 : p->y < q->y ? -1 : p->y > q->y ? 1 : 0;
}

/** The convex hull of the @p n points in `pts`, strictly convex and wound one way, into hx, hy; returns the count (0 when it has no area). */
static unsigned hull(Prof * p, unsigned n) {
  Pt * P = p->pts;
  qsort(P, n, sizeof *P, pt_less);
  unsigned m = 0;
  for (unsigned i = 0; i < n; i++) { /* drop duplicates */
    if (m == 0 || P[i].x != P[m - 1].x || P[i].y != P[m - 1].y) {
      P[m++] = P[i];
    }
  }
  if (m < 3) {
    return 0;
  }
  Pt * H = p->chain; /* Andrew's monotone chain: the lower hull, then the upper */
  unsigned k = 0;
  for (unsigned i = 0; i < m; i++) {
    while (k >= 2 && cross(H[k - 2].x, H[k - 2].y, H[k - 1].x, H[k - 1].y, P[i].x, P[i].y) <= 0) {
      k--;
    }
    H[k++] = P[i];
  }
  const unsigned lower = k + 1;
  for (unsigned i = m - 1; i-- > 0;) {
    while (k >= lower && cross(H[k - 2].x, H[k - 2].y, H[k - 1].x, H[k - 1].y, P[i].x, P[i].y) <= 0) {
      k--;
    }
    H[k++] = P[i];
  }
  k--; /* the last point is the first */
  if (k < 3) {
    return 0;
  }
  for (unsigned i = 0; i < k; i++) {
    p->hx[i] = (int32_t)H[i].x;
    p->hy[i] = (int32_t)H[i].y;
  }
  return k;
}

/** Emit the hull of pen @p a at (ax, ay) and pen @p b at (bx, by), centers in 26.6.  ERR_LIMIT past the points or the work or the range. */
static GSHGEOM_Result emit_hull(Prof * p, int64_t ax, int64_t ay, int64_t bx, int64_t by) {
  const unsigned n = p->a.n;
  if (p->effort > p->limits.max_work || 2u * n > p->limits.max_work - p->effort) {
    return GSHGEOM_ERR_LIMIT;
  }
  p->effort += 2u * n;
  for (unsigned k = 0; k < n; k++) {
    p->pts[k].x = ax + p->a.x[k];
    p->pts[k].y = ay + p->a.y[k];
    p->pts[n + k].x = bx + p->b.x[k];
    p->pts[n + k].y = by + p->b.y[k];
  }
  const int64_t lim = INT32_MAX;
  for (unsigned k = 0; k < 2u * n; k++) {
    if (p->pts[k].x > lim || p->pts[k].x < -lim || p->pts[k].y > lim || p->pts[k].y < -lim) {
      return GSHGEOM_ERR_LIMIT;
    }
  }
  const unsigned m = hull(p, 2u * n);
  if (m == 0) {
    return GSHGEOM_OK; /* a hull with no area draws nothing */
  }
  GSHGEOM_Device_Point q = {p->hx[0], p->hy[0]};
  GSHGEOM_Result r = gshgeom_flat_path_begin(p->body, q);
  for (unsigned i = 1; i < m && r == GSHGEOM_OK; i++) {
    q.x = p->hx[i];
    q.y = p->hy[i];
    r = gshgeom_flat_path_add(p->body, q);
  }
  if (r == GSHGEOM_OK) {
    gshgeom_flat_path_close_last(p->body);
  }
  return r;
}

/** The model length of the device vector (dx, dy) in 26.6. */
static double model_length(const Prof * p, double dx, double dy) {
  const double ux = (p->inv[0][0] * dx + p->inv[0][1] * dy) / 64.0;
  const double uy = (p->inv[1][0] * dx + p->inv[1][1] * dy) / 64.0;
  return gmath_sqrt(ux * ux + uy * uy);
}

/** The piece from device point (x0, y0) at position s0 to (x1, y1) at s1, cut at the profile's interior stops. */
static GSHGEOM_Result piece(Prof * p, int64_t x0, int64_t y0, double s0, int64_t x1, int64_t y1, double s1) {
  const GSHGEOM_Width_Profile * wp = p->wp;
  /* the stops strictly inside (s0, s1), in order */
  size_t lo = 0, hi = wp->count;
  while (lo < hi) { /* the first stop past s0 */
    size_t mid = lo + (hi - lo) / 2;
    if (wp->stops[mid].position > s0) {
      hi = mid;
    } else {
      lo = mid + 1;
    }
  }
  int64_t cx = x0, cy = y0;
  double cs = s0;
  int side = 1;
  GSHGEOM_Result r = GSHGEOM_OK;
  for (size_t i = lo; r == GSHGEOM_OK; i++) {
    const int last = i >= wp->count || wp->stops[i].position >= s1;
    const double es = last ? s1 : wp->stops[i].position;
    int64_t ex, ey;
    if (last) {
      ex = x1;
      ey = y1;
    } else {
      const double t = s1 > s0 ? (es - s0) / (s1 - s0) : 0.0;
      int32_t ox, oy;
      if (gshgeom_round_26_6((double)(x1 - x0) * t / 64.0, &ox) != GSHGEOM_OK || gshgeom_round_26_6((double)(y1 - y0) * t / 64.0, &oy) != GSHGEOM_OK) {
        return GSHGEOM_ERR_LIMIT;
      }
      ex = x0 + ox;
      ey = y0 + oy;
    }
    const double ra = p->width * factor_at(wp, cs, side) / 2.0;
    const double rb = p->width * factor_at(wp, es, -1) / 2.0;
    r = make_pen(p, &p->a, ra);
    if (r == GSHGEOM_OK) {
      r = make_pen(p, &p->b, rb);
    }
    if (r == GSHGEOM_OK) {
      r = emit_hull(p, cx, cy, ex, ey);
    }
    if (last) {
      break;
    }
    cx = ex;
    cy = ey;
    cs = es;
    side = 1;
  }
  return r;
}

/** One subpath of the flattened centerline. */
static GSHGEOM_Result subpath(Prof * p, const GSHGEOM_Device_Point * pt, size_t n, int closed, uint32_t degenerate) {
  if (n == 0) {
    return GSHGEOM_OK;
  }
  double s = 0.0;
  const size_t segs = closed ? n : n - 1;
  if (n == 1) {
    /* SW8: a point */
    if (degenerate == (uint32_t)GSHGEOM_STROKE_DEGENERATE_OMIT) {
      return GSHGEOM_OK;
    }
    return piece(p, pt[0].x, pt[0].y, 0.0, pt[0].x, pt[0].y, 0.0);
  }
  /* a jump at the very start has a limit that no piece begins with */
  if (factor_at(p->wp, 0.0, -1) != factor_at(p->wp, 0.0, 1)) {
    GSHGEOM_Result r0 = piece(p, pt[0].x, pt[0].y, 0.0, pt[0].x, pt[0].y, 0.0);
    if (r0 != GSHGEOM_OK) {
      return r0;
    }
  }
  for (size_t i = 0; i < segs; i++) {
    const GSHGEOM_Device_Point a = pt[i];
    const GSHGEOM_Device_Point b = pt[(i + 1) % n];
    const double len = model_length(p, (double)b.x - (double)a.x, (double)b.y - (double)a.y);
    GSHGEOM_Result r = piece(p, a.x, a.y, s, b.x, b.y, s + len);
    if (r != GSHGEOM_OK) {
      return r;
    }
    s += len;
  }
  /* and a jump at the very end has one that no piece ends with */
  if (factor_at(p->wp, s, -1) != factor_at(p->wp, s, 1)) {
    const GSHGEOM_Device_Point e = pt[segs % n];
    return piece(p, e.x, e.y, s, e.x, e.y, s);
  }
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_sk_profile_expand(const GSHGEOM_Path * centerline, const GSHGEOM_Stroke * stroke, const GSHGEOM_Map * map,
    int32_t tolerance, const GSHGEOM_Allocator * alloc, const GSHGEOM_Limits * limits, GSHGEOM_Stroke_Outline ** out,
    GSHGEOM_Stroke_Work * work) {
  GSHGEOM_Result r = gshgeom_stroke_check(stroke);
  if (r != GSHGEOM_OK) {
    return r;
  }
  r = gshgeom_map_check(map);
  if (r != GSHGEOM_OK || map->dim != 2) {
    return GSHGEOM_ERR_INVALID;
  }
  if (gshgeom_path_dimension(centerline) != 2) {
    return GSHGEOM_ERR_UNSUPPORTED;
  }
  if (!alloc) {
    alloc = gshgeom_allocator_default();
  }
  Prof p;
  memset(&p, 0, sizeof p);
  p.alloc = *alloc;
  gshgeom_limits_default(&p.limits);
  if (limits) {
    if (limits->max_elements) {
      p.limits.max_elements = limits->max_elements;
    }
    if (limits->max_points) {
      p.limits.max_points = limits->max_points;
    }
    if (limits->max_work) {
      p.limits.max_work = limits->max_work;
    }
  }
  p.wp = stroke->width_profile;
  p.width = stroke->width;
  p.local = stroke->width_space == (uint32_t)GSHGEOM_STROKE_WIDTH_LOCAL;
  GSHGEOM_Stroke_Outline * o = gcu_allocator_calloc(&p.alloc, 1, sizeof *o);
  if (!o) {
    return GSHGEOM_ERR_OOM;
  }
  o->alloc = p.alloc;
  o->alignment = (uint32_t)GSHGEOM_STROKE_ALIGN_CENTER;
  o->rule = GSHGEOM_FILL_NONZERO;
  int invertible = 0;
  r = gshgeom_stroke_map_invertible(map, &invertible);
  if (r != GSHGEOM_OK) {
    gshgeom_stroke_outline_destroy(o);
    return r;
  }
  if (!invertible) {
    o->status = GSHGEOM_STROKE_NOT_INVERTIBLE; /* SG15, before SG14 */
    *out = o;
    return GSHGEOM_OK;
  }
  if (stroke->width == 0.0) {
    o->status = GSHGEOM_STROKE_EMPTY;
    *out = o;
    return GSHGEOM_OK;
  }
  p.lin[0][0] = map->m[0][0];
  p.lin[0][1] = map->m[0][1];
  p.lin[1][0] = map->m[1][0];
  p.lin[1][1] = map->m[1][1];
  const double det = p.lin[0][0] * p.lin[1][1] - p.lin[0][1] * p.lin[1][0];
  p.inv[0][0] = p.lin[1][1] / det;
  p.inv[0][1] = -p.lin[0][1] / det;
  p.inv[1][0] = -p.lin[1][0] / det;
  p.inv[1][1] = p.lin[0][0] / det;
  /* The largest radius in pixels bounds the side count (SW6): the polygon lies inside its disc by at most R (1 - cos(pi / n)). */
  double fmax = 0.0;
  for (size_t i = 0; i < p.wp->count; i++) {
    fmax = p.wp->stops[i].factor > fmax ? p.wp->stops[i].factor : fmax;
  }
  const double norm = p.local ? gmath_sqrt(p.lin[0][0] * p.lin[0][0] + p.lin[0][1] * p.lin[0][1] + p.lin[1][0] * p.lin[1][0] + p.lin[1][1] * p.lin[1][1]) : 1.0;
  const double rmax = stroke->width * fmax / 2.0 * norm;
  p.err = (double)tolerance / 128.0; /* half the tolerance, in pixels */
  p.sides = kSides[sizeof kSides / sizeof kSides[0] - 1];
  for (size_t i = 0; i < sizeof kSides / sizeof kSides[0]; i++) {
    const double a = 3.14159265358979323846 / (double)kSides[i];
    if (rmax * (1.0 - gmath_cos(a)) <= p.err) {
      p.sides = kSides[i];
      break;
    }
  }
  const unsigned step = 512u / p.sides;
  (void)step;
  for (unsigned k = 0; k < 512; k++) {
    const double a = 6.283185307179586476925 * (double)k / 512.0;
    p.cs[k] = gmath_cos(a);
    p.sn[k] = gmath_sin(a);
  }
  const size_t cap = 1024;
  p.a.x = gcu_allocator_calloc(&p.alloc, cap, sizeof *p.a.x);
  p.a.y = gcu_allocator_calloc(&p.alloc, cap, sizeof *p.a.y);
  p.b.x = gcu_allocator_calloc(&p.alloc, cap, sizeof *p.b.x);
  p.b.y = gcu_allocator_calloc(&p.alloc, cap, sizeof *p.b.y);
  p.hx = gcu_allocator_calloc(&p.alloc, 2048, sizeof *p.hx);
  p.hy = gcu_allocator_calloc(&p.alloc, 2048, sizeof *p.hy);
  p.pts = gcu_allocator_calloc(&p.alloc, 2048, sizeof *p.pts);
  p.chain = gcu_allocator_calloc(&p.alloc, 2048, sizeof *p.chain);
  r = p.a.x && p.a.y && p.b.x && p.b.y && p.hx && p.hy && p.pts && p.chain ? GSHGEOM_OK : GSHGEOM_ERR_OOM;
  GSHGEOM_Flat_Path * flat = NULL;
  if (r == GSHGEOM_OK) {
    r = gshgeom_flat_path_new(&p.alloc, &p.limits, &o->body);
    p.body = o->body;
  }
  if (r == GSHGEOM_OK) {
    GSHGEOM_Device_Path * device = NULL;
    r = gshgeom_device_path_from_model_with_tolerance(centerline, map, tolerance, &p.alloc, &p.limits, &device);
    if (r == GSHGEOM_OK) {
      r = gshgeom_device_path_flatten_with_tolerance(device, tolerance, &p.alloc, &p.limits, &flat);
    }
    gshgeom_device_path_destroy(device);
  }
  if (r == GSHGEOM_OK) {
    const GSHGEOM_Device_Point * pts = gshgeom_flat_path_points(flat);
    const size_t subs = gshgeom_flat_path_subpath_count(flat);
    for (size_t i = 0; i < subs && r == GSHGEOM_OK; i++) {
      GSHGEOM_Flat_Subpath sp;
      r = gshgeom_flat_path_subpath(flat, i, &sp);
      if (r == GSHGEOM_OK) {
        r = subpath(&p, pts + sp.first, sp.count, sp.closed, stroke->degenerate);
      }
    }
  }
  gshgeom_flat_path_destroy(flat);
  gcu_allocator_free(&p.alloc, p.a.x);
  gcu_allocator_free(&p.alloc, p.a.y);
  gcu_allocator_free(&p.alloc, p.b.x);
  gcu_allocator_free(&p.alloc, p.b.y);
  gcu_allocator_free(&p.alloc, p.hx);
  gcu_allocator_free(&p.alloc, p.hy);
  gcu_allocator_free(&p.alloc, p.pts);
  gcu_allocator_free(&p.alloc, p.chain);
  if (work) {
    work->spent = p.effort;
  }
  if (r != GSHGEOM_OK) {
    gshgeom_stroke_outline_destroy(o);
    return r;
  }
  o->status = GSHGEOM_STROKE_DRAWN;
  *out = o;
  return GSHGEOM_OK;
}
