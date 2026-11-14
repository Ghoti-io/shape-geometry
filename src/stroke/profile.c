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
#include <ghoti.io/shape-geometry/orient.h>

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

/** One segment of the flattened centerline: its two end points in 26.6 device units, the profile positions at them, and the dash measure there. */
typedef struct Seg {
  int64_t x0, y0, x1, y1;       /* the part drawn: the segment, or what a dash leaves of it */
  int64_t ax0, ay0, ax1, ay1;   /* the whole segment, which is what a direction and a join are decided from (SG22) */
  double s0, s1; /* the position (SW3) at the start and the end */
  double d0, d1; /* the distance the pattern is measured in (SG19) */
  double mx0, my0, mx1, my1; /* the segment's two points before the map, when the centerline gives them (a line): what its direction is read from */
  int has_m;
} Seg;

/** The vertices of one subpath of the centerline before the map, when it is made of lines alone. */
typedef struct MSub {
  double * xy;
  size_t n, cap;
  int ok;
} MSub;

/** The first and the last step of a subpath, in device units: what a degenerate subpath's tangent is read from. */
typedef struct Ends {
  int has;
  double fx, fy, lx, ly;
} Ends;

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
  GSHGEOM_Flat_Path * body;  /* the part being drawn: the outline's body, or its aligned part */
  uint64_t effort;
  int hull_mode;      /* 1: round joins and round caps, the swept disc of SW4; 0: the trapezoids of SW10 */
  uint32_t join, cap; /* the canonical join and cap */
  double miter;       /* the miter limit */
  double detsign;     /* the sign of the linear part's determinant, 1 under device width */
  double sdet;        /* the sign of the linear part's determinant */
  int64_t ox, oy;     /* the origin of the build space, in 26.6 device units */
  double (*bq)[2];    /* the points of a shape under construction, in the build space */
  Seg * segs;         /* the subpath's segments */
  size_t seg_n, seg_cap;
  unsigned char * dotflag; /* with `dirs`: whether the subpath is a dot (a zero-length piece of the centerline), from the model path itself */
  double * piece_dir;      /* and the model direction of each piece that is not */
  const double * dirs; /* a hatch's dots (HX7): the model-space direction of each degenerate subpath, in order; NULL for a stroke */
  size_t dir_n, dir_i;
  double dtx, dty;    /* the effective tangent of the degenerate subpath being drawn (SG11), a unit vector of the build space */
  Seg * st;           /* the stretch being drawn: a dash, or the whole subpath */
  size_t st_n, st_cap;
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
  pen->n = n;
  for (unsigned k = 0; k < n; k++) {
    const double c = p->cs[k], s = p->sn[k];
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

/* ---------------------------------------------------------------------- */
/* The build space and the shapes of SW10                                  */
/* ---------------------------------------------------------------------- */

/** The point (x, y) of the device in the build space: model units relative to the origin under local width, pixels under device width. */
static void to_build(const Prof * p, int64_t x, int64_t y, double * u, double * v) {
  const double dx = (double)(x - p->ox) / 64.0, dy = (double)(y - p->oy) / 64.0;
  if (p->local) {
    *u = p->inv[0][0] * dx + p->inv[0][1] * dy;
    *v = p->inv[1][0] * dx + p->inv[1][1] * dy;
  } else {
    *u = dx;
    *v = dy;
  }
}

static GSHGEOM_Result from_build(const Prof * p, double u, double v, int64_t * x, int64_t * y) {
  double dx = u, dy = v;
  if (p->local) {
    dx = p->lin[0][0] * u + p->lin[0][1] * v;
    dy = p->lin[1][0] * u + p->lin[1][1] * v;
  }
  int32_t ix, iy;
  if (gshgeom_round_26_6(dx, &ix) != GSHGEOM_OK || gshgeom_round_26_6(dy, &iy) != GSHGEOM_OK) {
    return GSHGEOM_ERR_LIMIT;
  }
  *x = p->ox + ix;
  *y = p->oy + iy;
  return GSHGEOM_OK;
}

/** Emit the polygon of the @p m build points in `bq`, wound the one way whatever the map does; one with no area draws nothing. */
static GSHGEOM_Result emit_poly(Prof * p, unsigned m) {
  if (p->effort > p->limits.max_work || m > p->limits.max_work - p->effort) {
    return GSHGEOM_ERR_LIMIT;
  }
  p->effort += m;
  unsigned k = 0;
  for (unsigned i = 0; i < m; i++) {
    int64_t x, y;
    GSHGEOM_Result r = from_build(p, p->bq[i][0], p->bq[i][1], &x, &y);
    if (r != GSHGEOM_OK) {
      return r;
    }
    if (x > INT32_MAX || x < -(int64_t)INT32_MAX || y > INT32_MAX || y < -(int64_t)INT32_MAX) {
      return GSHGEOM_ERR_LIMIT;
    }
    if (k == 0 || p->pts[k - 1].x != x || p->pts[k - 1].y != y) {
      p->pts[k].x = x;
      p->pts[k].y = y;
      k++;
    }
  }
  while (k > 1 && p->pts[k - 1].x == p->pts[0].x && p->pts[k - 1].y == p->pts[0].y) {
    k--;
  }
  if (k < 3) {
    return GSHGEOM_OK;
  }
  double area = 0.0; /* twice the signed area, relative to the first point */
  for (unsigned i = 1; i + 1 < k; i++) {
    area += (double)(p->pts[i].x - p->pts[0].x) * (double)(p->pts[i + 1].y - p->pts[0].y) -
        (double)(p->pts[i + 1].x - p->pts[0].x) * (double)(p->pts[i].y - p->pts[0].y);
  }
  if (area == 0.0) {
    return GSHGEOM_OK;
  }
  GSHGEOM_Result r = GSHGEOM_OK;
  for (unsigned j = 0; j < k && r == GSHGEOM_OK; j++) {
    const unsigned i = area > 0.0 ? j : k - 1 - j;
    GSHGEOM_Device_Point q = {(int32_t)p->pts[i].x, (int32_t)p->pts[i].y};
    r = j == 0 ? gshgeom_flat_path_begin(p->body, q) : gshgeom_flat_path_add(p->body, q);
  }
  if (r == GSHGEOM_OK) {
    gshgeom_flat_path_close_last(p->body);
  }
  return r;
}

#define BQ(i, a, b) (p->bq[(i)][0] = (a), p->bq[(i)][1] = (b))
#define PI_ 3.14159265358979323846

/** The unit direction of a segment in the build space; 0 when it has none. */
static int seg_dir(const Prof * p, const Seg * s, double * dx, double * dy) {
  double u0, v0, u1, v1;
  if (s->has_m) {
    /* a line of the centerline: its direction is the one of its own points, not of what the lattice made of them */
    u0 = 0.0;
    v0 = 0.0;
    const double mx = s->mx1 - s->mx0, my = s->my1 - s->my0;
    u1 = p->local ? mx : p->lin[0][0] * mx + p->lin[0][1] * my;
    v1 = p->local ? my : p->lin[1][0] * mx + p->lin[1][1] * my;
  } else {
    to_build(p, s->ax0, s->ay0, &u0, &v0);
    to_build(p, s->ax1, s->ay1, &u1, &v1);
  }
  const double ex = u1 - u0, ey = v1 - v0;
  const double len = gmath_sqrt(ex * ex + ey * ey);
  if (!(len > 0.0)) {
    return 0;
  }
  *dx = ex / len;
  *dy = ey / len;
  return 1;
}

/** The radius at position @p s from the @p side, for a stroke whose width is `p->width`. */
static double radius_at(const Prof * p, double s, int side) {
  return p->width * factor_at(p->wp, s, side) / 2.0;
}

/** A trapezoid for the part of a segment between two stops: the normals at the two ends are the segment's, the radii are the profile's. */
static GSHGEOM_Result trap_seg(Prof * p, const Seg * sg) {
  double u0, v0, u1, v1;
  to_build(p, sg->x0, sg->y0, &u0, &v0);
  to_build(p, sg->x1, sg->y1, &u1, &v1);
  const double ex = u1 - u0, ey = v1 - v0;
  if (!(ex != 0.0 || ey != 0.0)) {
    return GSHGEOM_OK;
  }
  double ddx, ddy;
  if (!seg_dir(p, sg, &ddx, &ddy)) {
    return GSHGEOM_OK;
  }
  const double nx = -ddy, ny = ddx;
  const GSHGEOM_Width_Profile * wp = p->wp;
  size_t lo = 0, hi = wp->count;
  while (lo < hi) { /* the first stop past s0 */
    size_t mid = lo + (hi - lo) / 2;
    if (wp->stops[mid].position > sg->s0) {
      hi = mid;
    } else {
      lo = mid + 1;
    }
  }
  double cu = u0, cv = v0, cs = sg->s0;
  int side = 1;
  for (size_t i = lo;; i++) {
    const int last = i >= wp->count || wp->stops[i].position >= sg->s1;
    const double es = last ? sg->s1 : wp->stops[i].position;
    double eu, ev;
    if (last) {
      eu = u1;
      ev = v1;
    } else {
      const double t = sg->s1 > sg->s0 ? (es - sg->s0) / (sg->s1 - sg->s0) : 0.0;
      eu = u0 + ex * t;
      ev = v0 + ey * t;
    }
    const double ra = radius_at(p, cs, side), rb = radius_at(p, es, -1);
    BQ(0, cu + ra * nx, cv + ra * ny);
    BQ(1, eu + rb * nx, ev + rb * ny);
    BQ(2, eu - rb * nx, ev - rb * ny);
    BQ(3, cu - ra * nx, cv - ra * ny);
    GSHGEOM_Result r = emit_poly(p, 4);
    if (r != GSHGEOM_OK || last) {
      return r;
    }
    cu = eu;
    cv = ev;
    cs = es;
    side = 1;
  }
}

/** The join of SG8 at the vertex where segment @p a ends and @p b begins, for the radius @p r. */
static GSHGEOM_Result join_shape(Prof * p, const Seg * a, const Seg * b, double r) {
  if (!(r > 0.0)) {
    return GSHGEOM_OK;
  }
  double ax, ay, bx, by;
  if (!seg_dir(p, a, &ax, &ay) || !seg_dir(p, b, &bx, &by)) {
    return GSHGEOM_OK;
  }
  double vu, vv;
  to_build(p, a->x1, a->y1, &vu, &vv);
  int sg = 0;
  double dot_dev, turn;
  if (a->has_m && b->has_m) {
    GSHGEOM_Point qa = {a->mx0, a->my0, 0.0}, qv = {a->mx1, a->my1, 0.0}, qb = {b->mx1, b->my1, 0.0};
    if (gshgeom_orient2d(qa, qv, qb, &sg) != GSHGEOM_OK) {
      return GSHGEOM_ERR_INTERNAL;
    }
    dot_dev = (a->mx1 - a->mx0) * (b->mx1 - b->mx0) + (a->my1 - a->my0) * (b->my1 - b->my0);
    turn = (double)sg * (p->local ? 1.0 : p->sdet); /* the model's turn, in the build space */
  } else {
    GSHGEOM_Point pa = {(double)a->ax0, (double)a->ay0, 0.0}, pv = {(double)a->ax1, (double)a->ay1, 0.0}, pb = {(double)b->ax1, (double)b->ay1, 0.0};
    if (gshgeom_orient2d(pa, pv, pb, &sg) != GSHGEOM_OK) {
      return GSHGEOM_ERR_INTERNAL;
    }
    dot_dev = (double)(a->ax1 - a->ax0) * (double)(b->ax1 - b->ax0) + (double)(a->ay1 - a->ay0) * (double)(b->ay1 - b->ay0);
    turn = (double)sg * p->detsign; /* the lattice's turn, in the build space */
  }
  if (sg == 0 && dot_dev > 0.0) {
    return GSHGEOM_OK; /* a smooth vertex has no join */
  }
  const double dt = ax * bx + ay * by;
  if (sg == 0) {
    /* a reversal: the limit of each join as the turn tends to pi */
    const double nx = -ay, ny = ax;
    if (p->join == (uint32_t)GSHGEOM_STROKE_JOIN_ROUND) {
      const unsigned M = (p->sides + 1u) / 2u;
      for (unsigned k = 0; k <= M; k++) {
        const double al = PI_ / 2.0 - PI_ * (double)k / (double)M;
        double c = gmath_cos(al), s2 = gmath_sin(al);
        if (k == 0) {
          c = 0.0;
          s2 = 1.0;
        } else if (k == M) {
          c = 0.0;
          s2 = -1.0;
        }
        BQ(k, vu + r * (s2 * nx + c * ax), vv + r * (s2 * ny + c * ay));
      }
      return emit_poly(p, M + 1u);
    }
    if (p->join == (uint32_t)GSHGEOM_STROKE_JOIN_MITER_CLIP || p->join == (uint32_t)GSHGEOM_STROKE_JOIN_ARCS) {
      const double L = p->miter * r;
      BQ(0, vu + r * nx, vv + r * ny);
      BQ(1, vu + r * nx + L * ax, vv + r * ny + L * ay);
      BQ(2, vu - r * nx + L * ax, vv - r * ny + L * ay);
      BQ(3, vu - r * nx, vv - r * ny);
      return emit_poly(p, 4);
    }
    return GSHGEOM_OK;
  }
  /* the outer normals: the right normal of a left turn, the left normal of a right turn */
  const double oax = turn > 0.0 ? ay : -ay, oay = turn > 0.0 ? -ax : ax;
  const double obx = turn > 0.0 ? by : -by, oby = turn > 0.0 ? -bx : bx;
  const double Au = vu + r * oax, Av = vv + r * oay, Bu = vu + r * obx, Bv = vv + r * oby;
  if (p->join == (uint32_t)GSHGEOM_STROKE_JOIN_ROUND) {
    double cr = ax * by - ay * bx;
    double phi = gmath_atan2(cr < 0.0 ? -cr : cr, dt);
    unsigned m = (unsigned)gmath_ceil(phi * (double)p->sides / (2.0 * PI_));
    m = m < 1u ? 1u : m;
    BQ(0, vu, vv);
    BQ(1, Au, Av);
    const double ca = turn > 0.0 ? 1.0 : -1.0;
    for (unsigned k = 1; k < m; k++) {
      const double th = ca * phi * (double)k / (double)m;
      const double c = gmath_cos(th), s2 = gmath_sin(th);
      BQ(1u + k, vu + r * (c * oax - s2 * oay), vv + r * (s2 * oax + c * oay));
    }
    BQ(1u + m, Bu, Bv);
    return emit_poly(p, 2u + m);
  }
  const int within = (1.0 + dt) / 2.0 >= 1.0 / (p->miter * p->miter);
  const int is_miter = p->join == (uint32_t)GSHGEOM_STROKE_JOIN_MITER;
  if (p->join == (uint32_t)GSHGEOM_STROKE_JOIN_BEVEL || (is_miter && !within)) {
    BQ(0, vu, vv);
    BQ(1, Au, Av);
    BQ(2, Bu, Bv);
    return emit_poly(p, 3);
  }
  const double den = 1.0 + dt;
  if (within) {
    BQ(0, vu, vv);
    BQ(1, Au, Av);
    BQ(2, vu + r * (oax + obx) / den, vv + r * (oay + oby) / den);
    BQ(3, Bu, Bv);
    return emit_poly(p, 4);
  }
  /* miter-clip and arcs past the limit (the polyline's curvature is zero): the pentagon V A A' B' B */
  double ux = oax + obx, uy = oay + oby;
  const double ul = gmath_sqrt(ux * ux + uy * uy);
  ux /= ul;
  uy /= ul;
  const double c2 = oax * ux + oay * uy;                    /* cos(phi / 2) */
  const double au = ax * ux + ay * uy, bu = -(bx * ux + by * uy); /* a . u and -b . u, both positive */
  const double ta = (p->miter * r - r * c2) / au, tb = (p->miter * r - r * c2) / bu;
  BQ(0, vu, vv);
  BQ(1, Au, Av);
  BQ(2, Au + ta * ax, Av + ta * ay);
  BQ(3, Bu - tb * bx, Bv - tb * by);
  BQ(4, Bu, Bv);
  return emit_poly(p, 5);
}

/** The cap of SG7 at the end (@p sx, @p sy) with the outward unit tangent (dx, dy) and radius @p r. */
static GSHGEOM_Result cap_shape(Prof * p, int64_t sx, int64_t sy, double dx, double dy, double r) {
  if (!(r > 0.0) || p->cap == (uint32_t)GSHGEOM_STROKE_CAP_BUTT) {
    return GSHGEOM_OK;
  }
  double eu, ev;
  to_build(p, sx, sy, &eu, &ev);
  const double mx = -dy, my = dx;
  if (p->cap == (uint32_t)GSHGEOM_STROKE_CAP_SQUARE) {
    BQ(0, eu + r * mx, ev + r * my);
    BQ(1, eu + r * mx + r * dx, ev + r * my + r * dy);
    BQ(2, eu - r * mx + r * dx, ev - r * my + r * dy);
    BQ(3, eu - r * mx, ev - r * my);
    return emit_poly(p, 4);
  }
  const unsigned M = (p->sides + 1u) / 2u;
  for (unsigned k = 0; k <= M; k++) {
    double c, s2;
    if (k == 0) {
      c = 0.0;
      s2 = 1.0;
    } else if (k == M) {
      c = 0.0;
      s2 = -1.0;
    } else {
      const double al = PI_ / 2.0 - PI_ * (double)k / (double)M;
      c = gmath_cos(al);
      s2 = gmath_sin(al);
    }
    BQ(k, eu + r * (s2 * mx + c * dx), ev + r * (s2 * my + c * dy));
  }
  return emit_poly(p, M + 1u);
}

/** The dot of SG24 (and SG11) at a point with the tangent (dx, dy): the disc, the square, or nothing, for the radius @p r. */
static GSHGEOM_Result dot_shape(Prof * p, int64_t x, int64_t y, double dx, double dy, double r) {
  if (!(r > 0.0) || p->cap == (uint32_t)GSHGEOM_STROKE_CAP_BUTT) {
    return GSHGEOM_OK;
  }
  double eu, ev;
  to_build(p, x, y, &eu, &ev);
  if (p->cap == (uint32_t)GSHGEOM_STROKE_CAP_SQUARE) {
    const double mx = -dy, my = dx;
    BQ(0, eu + r * (dx + mx), ev + r * (dy + my));
    BQ(1, eu + r * (-dx + mx), ev + r * (-dy + my));
    BQ(2, eu + r * (-dx - mx), ev + r * (-dy - my));
    BQ(3, eu + r * (dx - mx), ev + r * (dy - my));
    return emit_poly(p, 4);
  }
  for (unsigned k = 0; k < p->sides; k++) {
    BQ(k, eu + r * p->cs[k], ev + r * p->sn[k]);
  }
  return emit_poly(p, p->sides);
}

/* ---------------------------------------------------------------------- */
/* A run of segments: a dash, or a subpath                                 */
/* ---------------------------------------------------------------------- */

/** The disc (a hull of the pens at the one point) at a position, for both limits: SW8. */
static GSHGEOM_Result point_pen(Prof * p, int64_t x, int64_t y, double s) {
  return piece(p, x, y, s, x, y, s);
}

/** A dot at the point (x, y), position @p s, tangent (dx, dy), of both limits of the profile there. */
static GSHGEOM_Result dot(Prof * p, int64_t x, int64_t y, double s, double dx, double dy) {
  if (p->hull_mode) {
    return point_pen(p, x, y, s);
  }
  GSHGEOM_Result r = dot_shape(p, x, y, dx, dy, radius_at(p, s, -1));
  const double r1 = radius_at(p, s, 1);
  if (r == GSHGEOM_OK && r1 != radius_at(p, s, -1)) {
    r = dot_shape(p, x, y, dx, dy, r1);
  }
  return r;
}

/**
 * Draw the @p n segments of @p s: a closed loop (a join at every vertex, the last to the first, and no caps) or an open run (a cap at each
 * end), the join at a vertex between two segments being the one the profile's limits there ask for.
 */
static GSHGEOM_Result run(Prof * p, const Seg * s, size_t n, int closed) {
  GSHGEOM_Result r = GSHGEOM_OK;
  if (p->hull_mode) {
    /* a jump at the first or the last position of the run has a limit that no piece begins or ends with: its disc is the point's */
    if (factor_at(p->wp, s[0].s0, -1) != factor_at(p->wp, s[0].s0, 1)) {
      r = point_pen(p, s[0].x0, s[0].y0, s[0].s0);
    }
    for (size_t i = 0; i < n && r == GSHGEOM_OK; i++) {
      r = piece(p, s[i].x0, s[i].y0, s[i].s0, s[i].x1, s[i].y1, s[i].s1);
    }
    if (r == GSHGEOM_OK && factor_at(p->wp, s[n - 1].s1, -1) != factor_at(p->wp, s[n - 1].s1, 1)) {
      r = point_pen(p, s[n - 1].x1, s[n - 1].y1, s[n - 1].s1);
    }
    return r;
  }
  for (size_t i = 0; i < n && r == GSHGEOM_OK; i++) {
    r = trap_seg(p, &s[i]);
  }
  const size_t joins = closed ? n : n - 1;
  for (size_t i = 0; i < joins && r == GSHGEOM_OK; i++) {
    const Seg * a = &s[i];
    const Seg * b = &s[(i + 1) % n];
    const double r0 = radius_at(p, a->s1, -1), r1 = radius_at(p, b->s0, 1);
    r = join_shape(p, a, b, r0);
    if (r == GSHGEOM_OK && r1 != r0) {
      r = join_shape(p, a, b, r1);
    }
  }
  if (!closed && r == GSHGEOM_OK) {
    double dx, dy;
    if (seg_dir(p, &s[0], &dx, &dy)) {
      for (int side = -1; side <= 1 && r == GSHGEOM_OK; side += 2) {
        const double rr = radius_at(p, s[0].s0, side);
        if (side == 1 && rr == radius_at(p, s[0].s0, -1)) {
          break;
        }
        r = cap_shape(p, s[0].x0, s[0].y0, -dx, -dy, rr);
      }
    }
    if (r == GSHGEOM_OK && seg_dir(p, &s[n - 1], &dx, &dy)) {
      for (int side = -1; side <= 1 && r == GSHGEOM_OK; side += 2) {
        const double rr = radius_at(p, s[n - 1].s1, side);
        if (side == 1 && rr == radius_at(p, s[n - 1].s1, -1)) {
          break;
        }
        r = cap_shape(p, s[n - 1].x1, s[n - 1].y1, dx, dy, rr);
      }
    }
  }
  return r;
}

/** Append the part of segment @p g between its parameters @p t0 and @p t1 to the stretch. */
static GSHGEOM_Result push_part(Prof * p, const Seg * g, double t0, double t1) {
  void * mem = p->st;
  GSHGEOM_Result r = gshgeom_grow(&p->alloc, &mem, &p->st_cap, p->st_n + 1, sizeof *p->st);
  if (r != GSHGEOM_OK) {
    return r;
  }
  p->st = mem;
  Seg o = *g;
  const double dx = (double)(g->x1 - g->x0), dy = (double)(g->y1 - g->y0);
  if (t0 > 0.0) {
    o.x0 = g->x0 + (int64_t)gmath_floor(dx * t0 + 0.5);
    o.y0 = g->y0 + (int64_t)gmath_floor(dy * t0 + 0.5);
    o.s0 = g->s0 + (g->s1 - g->s0) * t0;
    o.d0 = g->d0 + (g->d1 - g->d0) * t0;
  }
  if (t1 < 1.0) {
    o.x1 = g->x0 + (int64_t)gmath_floor(dx * t1 + 0.5);
    o.y1 = g->y0 + (int64_t)gmath_floor(dy * t1 + 0.5);
    o.s1 = g->s0 + (g->s1 - g->s0) * t1;
    o.d1 = g->d0 + (g->d1 - g->d0) * t1;
  }
  p->st[p->st_n++] = o;
  return GSHGEOM_OK;
}

/** The segments of the subpath between the distances @p a and @p b (a < b), appended to the stretch. */
static GSHGEOM_Result take(Prof * p, double a, double b) {
  const Seg * g = p->segs;
  const size_t n = p->seg_n;
  size_t lo = 0, hi = n; /* the first segment that ends past a */
  while (lo < hi) {
    size_t mid = lo + (hi - lo) / 2;
    if (g[mid].d1 > a) {
      hi = mid;
    } else {
      lo = mid + 1;
    }
  }
  for (size_t k = lo; k < n && g[k].d0 < b; k++) {
    const double w = g[k].d1 - g[k].d0;
    double t0 = a > g[k].d0 ? (a - g[k].d0) / w : 0.0;
    double t1 = b < g[k].d1 ? (b - g[k].d0) / w : 1.0;
    if (t1 > t0) {
      GSHGEOM_Result r = push_part(p, &g[k], t0, t1);
      if (r != GSHGEOM_OK) {
        return r;
      }
    }
  }
  return GSHGEOM_OK;
}

/** The dot of a dash at the distance @p a: the position, the point and the tangent of the segment that begins there. */
static GSHGEOM_Result dash_dot(Prof * p, double a) {
  const Seg * g = p->segs;
  size_t lo = 0, hi = p->seg_n;
  while (lo < hi) {
    size_t mid = lo + (hi - lo) / 2;
    if (g[mid].d1 > a) {
      hi = mid;
    } else {
      lo = mid + 1;
    }
  }
  const int at_end = lo >= p->seg_n;
  const Seg * k = &g[at_end ? p->seg_n - 1 : lo];
  const double t = at_end ? 1.0 : (a - k->d0) / (k->d1 - k->d0);
  int64_t x = k->x1, y = k->y1;
  double s = k->s1;
  if (!at_end && t > 0.0) {
    x = k->x0 + (int64_t)gmath_floor((double)(k->x1 - k->x0) * t + 0.5);
    y = k->y0 + (int64_t)gmath_floor((double)(k->y1 - k->y0) * t + 0.5);
    s = k->s0 + (k->s1 - k->s0) * t;
  } else if (!at_end) {
    x = k->x0;
    y = k->y0;
    s = k->s0;
  }
  double dx = 1.0, dy = 0.0;
  (void)seg_dir(p, k, &dx, &dy);
  return dot(p, x, y, s, dx, dy);
}

/** What the dashes of one subpath need. */
typedef struct DashCtx {
  Prof * p;
  int closed;
  int wrap;      /* the first and the last dash are one (SG23) */
  size_t index;
  size_t total;
  double held;   /* the end of the first dash, kept for the wrap */
} DashCtx;

static GSHGEOM_Result census_fn(void * ctx, double a, double b, int flags) {
  double * k = ctx; /* n, first_end, last_start, first_flags, last_flags */
  if (k[0] == 0.0) {
    k[1] = b;
    k[3] = (double)flags;
  }
  k[2] = a;
  k[4] = (double)flags;
  k[0] += 1.0;
  return GSHGEOM_OK;
}

static GSHGEOM_Result on_dash(void * ctx, double a, double b, int flags) {
  (void)flags;
  DashCtx * d = ctx;
  Prof * p = d->p;
  const size_t idx = d->index++;
  const double total = p->segs[p->seg_n - 1].d1;
  if (d->wrap && idx == 0) {
    d->held = b;
    return GSHGEOM_OK;
  }
  p->st_n = 0;
  GSHGEOM_Result r;
  if (d->wrap && idx + 1 == d->total) {
    r = take(p, a, total);
    if (r == GSHGEOM_OK) {
      r = take(p, 0.0, d->held);
    }
  } else {
    if (b > total) {
      b = total;
    }
    if (!(b > a)) {
      return dash_dot(p, a);
    }
    r = take(p, a, b);
  }
  if (r != GSHGEOM_OK) {
    return r;
  }
  if (p->st_n == 0) {
    return dash_dot(p, a);
  }
  return run(p, p->st, p->st_n, 0);
}

/** The dashes of the subpath held in `segs`, `begin` being the distance along the path at its start. */
static GSHGEOM_Result dashed(Prof * p, const GSHGEOM_Stroke * canon, const GSHGEOM_Stroke * raw, int closed, double begin, size_t * steps) {
  const size_t n = canon->dash_count;
  double dashes_stack[64];
  double * dashes = dashes_stack;
  if (n > 64) {
    dashes = gcu_allocator_calloc(&p->alloc, n, sizeof *dashes);
    if (!dashes) {
      return GSHGEOM_ERR_OOM;
    }
  }
  for (size_t i = 0; i < n; i++) {
    const double v = raw->dashes[i % raw->dash_count];
    dashes[i] = v == 0.0 ? 0.0 : v;
  }
  const double period = gshgeom_sk_dash_period(dashes, n, n);
  const double L = p->segs[p->seg_n - 1].d1;
  const double o = gshgeom_sk_dash_phase(canon->dash_scope == (uint32_t)GSHGEOM_STROKE_DASH_PATH ? canon->dash_offset + begin : canon->dash_offset, period);
  const size_t cap = p->limits.max_dash_steps && p->limits.max_dash_steps < (uint64_t)SIZE_MAX ? (size_t)p->limits.max_dash_steps : 4 * p->limits.max_points + 4096;
  DashCtx d = {p, closed, 0, 0, 0, 0.0};
  GSHGEOM_Result r = GSHGEOM_OK;
  if (closed) {
    double k[5] = {0, 0, 0, 0, 0};
    size_t st = 0;
    r = gshgeom_sk_dash_positions(dashes, n, n, period, o, L, &st, cap, census_fn, k);
    if (r == GSHGEOM_OK && k[0] > 0.0 && ((int)k[3] & GSHGEOM_SK_POS_BEGINS_BEFORE) && ((int)k[4] & GSHGEOM_SK_POS_ENDS_AFTER)) {
      if (k[0] == 1.0 || !(k[1] < k[2])) {
        r = run(p, p->segs, p->seg_n, 1); /* the dash covers the whole subpath: its stroke, closed */
        goto done;
      }
      d.wrap = 1;
      d.total = (size_t)k[0];
    }
  }
  if (r == GSHGEOM_OK) {
    r = gshgeom_sk_dash_positions(dashes, n, n, period, o, L, steps, cap, on_dash, &d);
  }
done:
  if (dashes != dashes_stack) {
    gcu_allocator_free(&p->alloc, dashes);
  }
  return r;
}

typedef struct DegenCtx {
  Prof * p;
  int64_t x, y;
} DegenCtx;

static GSHGEOM_Result degen_fn(void * ctx, double a, double b, int flags) {
  (void)a;
  (void)b;
  (void)flags;
  DegenCtx * d = ctx;
  return dot(d->p, d->x, d->y, 0.0, d->p->dtx, d->p->dty);
}

/** One subpath of the flattened centerline, drawn into the part @p p->body at the width in `p->width`.  @p begin: the distance at its start. */
static GSHGEOM_Result subpath(Prof * p, const GSHGEOM_Stroke * canon, const GSHGEOM_Stroke * raw, const GSHGEOM_Device_Point * pt,
    size_t n, int closed, const MSub * ms, double * begin, size_t * steps) {
  const int use_m = ms && ms->ok && ms->n == n;
  size_t prev_i = 0;
  /* the distinct points: a repeated one makes a segment with no extent */
  p->seg_n = 0;
  p->ox = n ? pt[0].x : 0;
  p->oy = n ? pt[0].y : 0;
  double s = 0.0, dd = 0.0;
  GSHGEOM_Device_Point prev = n ? pt[0] : (GSHGEOM_Device_Point){0, 0};
  const size_t total = closed ? n + 1 : n;
  for (size_t i = 1; i < total; i++) {
    const GSHGEOM_Device_Point b = pt[i % n];
    if (b.x == prev.x && b.y == prev.y) {
      continue;
    }
    void * mem = p->segs;
    GSHGEOM_Result r = gshgeom_grow(&p->alloc, &mem, &p->seg_cap, p->seg_n + 1, sizeof *p->segs);
    if (r != GSHGEOM_OK) {
      return r;
    }
    p->segs = mem;
    const double dx = (double)b.x - (double)prev.x, dy = (double)b.y - (double)prev.y;
    const double len = model_length(p, dx, dy);
    const double dl = p->local ? len : gmath_sqrt(dx * dx + dy * dy) / 64.0;
    Seg g = {prev.x, prev.y, b.x, b.y, prev.x, prev.y, b.x, b.y, s, s + len, dd, dd + dl, 0.0, 0.0, 0.0, 0.0, 0};
    if (use_m) {
      g.has_m = 1;
      g.mx0 = ms->xy[2 * prev_i];
      g.my0 = ms->xy[2 * prev_i + 1];
      g.mx1 = ms->xy[2 * (i % n)];
      g.my1 = ms->xy[2 * (i % n) + 1];
    }
    prev_i = i % n;
    p->segs[p->seg_n++] = g;
    s += len;
    dd += dl;
    prev = b;
  }
  if (p->seg_n == 0) {
    /* SW8 and SG11: a subpath with no extent is a point */
    if (n == 0) {
      return GSHGEOM_OK;
    }
    if (canon->degenerate == (uint32_t)GSHGEOM_STROKE_DEGENERATE_OMIT && !p->dirs) {
      return GSHGEOM_OK;
    }
    /* with dashes the pattern is on at its point or it is not (SG21 with a length of 0) */
    if (canon->dash_count) {
      const size_t nd = canon->dash_count;
      double ds[64];
      if (nd > 64) {
        return GSHGEOM_ERR_LIMIT;
      }
      for (size_t i = 0; i < nd; i++) {
        const double v = raw->dashes[i % raw->dash_count];
        ds[i] = v == 0.0 ? 0.0 : v;
      }
      const double period = gshgeom_sk_dash_period(ds, nd, nd);
      const double o = gshgeom_sk_dash_phase(canon->dash_scope == (uint32_t)GSHGEOM_STROKE_DASH_PATH ? canon->dash_offset + *begin : canon->dash_offset, period);
      DegenCtx dc = {p, pt[0].x, pt[0].y};
      size_t st = 0;
      return gshgeom_sk_dash_positions(ds, nd, nd, period, o, 0.0, &st, 4, degen_fn, &dc);
    }
    return dot(p, pt[0].x, pt[0].y, 0.0, p->dtx, p->dty);
  }
  /* an open subpath that closed back on its own point is a closed loop of the segments it has; a closed one whose last point repeats the first */
  const int loop = closed;
  GSHGEOM_Result r;
  if (canon->dash_count) {
    r = dashed(p, canon, raw, loop, *begin, steps);
  } else {
    r = run(p, p->segs, p->seg_n, loop);
  }
  *begin += dd;
  return r;
}

/** The vertices of each subpath of @p path that has lines alone, before the map; @p subs must be how many subpaths the flattener made of it. */
static GSHGEOM_Result read_model(Prof * p, const GSHGEOM_Path * path, size_t subs, MSub ** out) {
  *out = NULL;
  const size_t ne = gshgeom_path_element_count(path);
  MSub * v = gcu_allocator_calloc(&p->alloc, subs ? subs : 1, sizeof *v);
  if (!v) {
    return GSHGEOM_ERR_OOM;
  }
  size_t k = 0;
  int open = 0; /* a subpath is being read */
  GSHGEOM_Result r = GSHGEOM_OK;
  double sx = 0.0, sy = 0.0; /* the start of the subpath being read, for a drawing verb after a close */
  for (size_t e = 0; e < ne && r == GSHGEOM_OK; e++) {
    GSHGEOM_Verb verb;
    (void)gshgeom_path_verb(path, e, &verb);
    GSHGEOM_Element el;
    const int have = verb != GSHGEOM_VERB_ARC && gshgeom_path_element(path, e, &el) == GSHGEOM_OK;
    if (verb == GSHGEOM_VERB_MOVE || (!open && verb != GSHGEOM_VERB_CLOSE)) {
      if (k >= subs) {
        k = subs + 1;
        break;
      }
      MSub * m = &v[k++];
      m->ok = 1;
      open = 1;
      void * mem = m->xy;
      r = gshgeom_grow(&p->alloc, &mem, &m->cap, 2, sizeof(double) * 2);
      m->xy = mem;
      if (r != GSHGEOM_OK) {
        break;
      }
      if (verb == GSHGEOM_VERB_MOVE && have) {
        sx = el.curve.p[0].x;
        sy = el.curve.p[0].y;
      }
      m->xy[0] = sx;
      m->xy[1] = sy;
      m->n = 1;
      if (verb == GSHGEOM_VERB_MOVE) {
        continue;
      }
    }
    MSub * m = &v[k - 1];
    if (verb == GSHGEOM_VERB_CLOSE) {
      open = 0;
      continue;
    }
    if (verb != GSHGEOM_VERB_LINE || !have) {
      m->ok = 0;
      continue;
    }
    void * mem = m->xy;
    r = gshgeom_grow(&p->alloc, &mem, &m->cap, m->n + 1, sizeof(double) * 2);
    m->xy = mem;
    if (r != GSHGEOM_OK) {
      break;
    }
    m->xy[2 * m->n] = el.curve.p[1].x;
    m->xy[2 * m->n + 1] = el.curve.p[1].y;
    m->n++;
  }
  if (r == GSHGEOM_OK && k != subs) {
    for (size_t i = 0; i < subs; i++) {
      v[i].ok = 0; /* the subpaths of the path and of the flat path do not line up: read no direction from the path */
    }
  }
  if (r != GSHGEOM_OK) {
    for (size_t i = 0; i < subs; i++) {
      gcu_allocator_free(&p->alloc, v[i].xy);
    }
    gcu_allocator_free(&p->alloc, v);
    return r;
  }
  *out = v;
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_sk_profile_expand(const GSHGEOM_Path * centerline, const GSHGEOM_Stroke * stroke, const GSHGEOM_Map * map,
    int32_t tolerance, const GSHGEOM_Allocator * alloc, const GSHGEOM_Limits * limits, const double * dirs, size_t dir_count,
    GSHGEOM_Stroke_Outline ** out, GSHGEOM_Stroke_Work * work) {
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
  GSHGEOM_Stroke canon;
  r = gshgeom_stroke_canonical(stroke, &canon);
  if (r != GSHGEOM_OK) {
    return r;
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
    p.limits.max_dash_steps = limits->max_dash_steps;
  }
  p.wp = stroke->width_profile;
  p.dirs = dirs;
  p.dir_n = dir_count;
  p.local = canon.width_space == (uint32_t)GSHGEOM_STROKE_WIDTH_LOCAL;
  p.join = canon.join;
  p.cap = canon.cap;
  p.miter = canon.miter_limit;
  p.hull_mode = canon.join == (uint32_t)GSHGEOM_STROKE_JOIN_ROUND && canon.cap == (uint32_t)GSHGEOM_STROKE_CAP_ROUND;
  const int aligned = canon.alignment != (uint32_t)GSHGEOM_STROKE_ALIGN_CENTER;
  GSHGEOM_Stroke_Outline * o = gcu_allocator_calloc(&p.alloc, 1, sizeof *o);
  if (!o) {
    return GSHGEOM_ERR_OOM;
  }
  o->alloc = p.alloc;
  o->alignment = canon.alignment;
  o->rule = canon.alignment_rule;
  o->complement = canon.alignment == (uint32_t)GSHGEOM_STROKE_ALIGN_OUTSIDE;
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
  p.detsign = !p.local || det > 0.0 ? 1.0 : -1.0;
  p.sdet = det > 0.0 ? 1.0 : -1.0;
  /* The largest radius in pixels bounds the side count (SW6): the polygon lies inside its disc by at most R (1 - cos(pi / n)). */
  double fmax = 0.0;
  for (size_t i = 0; i < p.wp->count; i++) {
    fmax = p.wp->stops[i].factor > fmax ? p.wp->stops[i].factor : fmax;
  }
  const double norm = p.local ? gmath_sqrt(p.lin[0][0] * p.lin[0][0] + p.lin[0][1] * p.lin[0][1] + p.lin[1][0] * p.lin[1][0] + p.lin[1][1] * p.lin[1][1]) : 1.0;
  const double rmax = stroke->width * (aligned ? 2.0 : 1.0) * fmax / 2.0 * norm;
  p.err = (double)tolerance / 128.0; /* half the tolerance, in pixels */
  p.sides = kSides[sizeof kSides / sizeof kSides[0] - 1];
  for (size_t i = 0; i < sizeof kSides / sizeof kSides[0]; i++) {
    const double a = 3.14159265358979323846 / (double)kSides[i];
    if (rmax * (1.0 - gmath_cos(a)) <= p.err) {
      p.sides = kSides[i];
      break;
    }
  }
  for (unsigned k = 0; k < p.sides; k++) {
    const double a = 6.283185307179586476925 * (double)k / (double)p.sides;
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
  p.bq = gcu_allocator_calloc(&p.alloc, 1100, sizeof *p.bq);
  r = p.a.x && p.a.y && p.b.x && p.b.y && p.hx && p.hy && p.pts && p.chain && p.bq ? GSHGEOM_OK : GSHGEOM_ERR_OOM;
  GSHGEOM_Flat_Path * flat = NULL;
  Ends * ends = NULL;
  MSub * msubs = NULL;
  size_t msub_n = 0;
  if (r == GSHGEOM_OK) {
    r = gshgeom_flat_path_new(&p.alloc, &p.limits, &o->body);
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
    msub_n = subs;
    /* SG12: the closed subpaths with extent are R_c and the clip F; a distinct point count of 2 or more is extent */
    if (aligned) {
      int any = 0;
      for (size_t i = 0; i < subs && r == GSHGEOM_OK; i++) {
        GSHGEOM_Flat_Subpath sp;
        r = gshgeom_flat_path_subpath(flat, i, &sp);
        if (r != GSHGEOM_OK || !sp.closed) {
          continue;
        }
        int ext = 0;
        for (size_t k = 1; k < sp.count; k++) {
          ext |= pts[sp.first + k].x != pts[sp.first].x || pts[sp.first + k].y != pts[sp.first].y;
        }
        if (!ext) {
          continue;
        }
        if (!any) {
          any = 1;
          r = gshgeom_flat_path_new(&p.alloc, &p.limits, &o->aligned);
          if (r == GSHGEOM_OK) {
            r = gshgeom_flat_path_new(&p.alloc, &p.limits, &o->clip);
          }
        }
        for (size_t k = 0; k < sp.count && r == GSHGEOM_OK; k++) {
          r = k == 0 ? gshgeom_flat_path_begin(o->clip, pts[sp.first]) : gshgeom_flat_path_add(o->clip, pts[sp.first + k]);
        }
        if (r == GSHGEOM_OK) {
          gshgeom_flat_path_close_last(o->clip);
        }
      }
    }
    /* SG11's effective tangent of a degenerate subpath: the end direction of the last segment with extent before it, else the start
     * direction of the first after it, else the map's image of +x; found from the first and last steps of each subpath */
    ends = gcu_allocator_calloc(&p.alloc, subs ? subs : 1, sizeof *ends);
    if (!ends) {
      r = GSHGEOM_ERR_OOM;
    }
    for (size_t i = 0; i < subs && r == GSHGEOM_OK; i++) {
      GSHGEOM_Flat_Subpath sp;
      r = gshgeom_flat_path_subpath(flat, i, &sp);
      if (r != GSHGEOM_OK) {
        break;
      }
      const GSHGEOM_Device_Point * q = pts + sp.first;
      GSHGEOM_Device_Point prev = q[0];
      const size_t total = sp.closed ? sp.count + 1 : sp.count;
      for (size_t k = 1; k < total && sp.count; k++) {
        const GSHGEOM_Device_Point b = q[k % sp.count];
        if (b.x == prev.x && b.y == prev.y) {
          continue;
        }
        if (!ends[i].has) {
          ends[i].has = 1;
          ends[i].fx = (double)b.x - (double)prev.x;
          ends[i].fy = (double)b.y - (double)prev.y;
        }
        ends[i].lx = (double)b.x - (double)prev.x;
        ends[i].ly = (double)b.y - (double)prev.y;
        prev = b;
      }
    }
    if (p.dirs && r == GSHGEOM_OK) {
      /* the hatch's centerline is a move and a line for each piece: a dot is the one whose two points are the same */
      const size_t ne = gshgeom_path_element_count(centerline);
      if (ne == 2 * subs && subs) {
        p.dotflag = gcu_allocator_calloc(&p.alloc, subs, sizeof *p.dotflag);
        p.piece_dir = gcu_allocator_calloc(&p.alloc, subs, 2 * sizeof *p.piece_dir);
        if (!p.dotflag || !p.piece_dir) {
          r = GSHGEOM_ERR_OOM;
        }
        for (size_t i = 0; i < subs && r == GSHGEOM_OK; i++) {
          GSHGEOM_Element e0, e1;
          if (gshgeom_path_element(centerline, 2 * i, &e0) != GSHGEOM_OK || gshgeom_path_element(centerline, 2 * i + 1, &e1) != GSHGEOM_OK) {
            r = GSHGEOM_ERR_INTERNAL;
            break;
          }
          const double dx = e1.curve.p[1].x - e1.curve.p[0].x, dy = e1.curve.p[1].y - e1.curve.p[0].y;
          p.dotflag[i] = dx == 0.0 && dy == 0.0;
          p.piece_dir[2 * i] = dx;
          p.piece_dir[2 * i + 1] = dy;
        }
      }
    }
    if (r == GSHGEOM_OK) {
      r = read_model(&p, centerline, subs, &msubs);
    }
    double begin = 0.0;
    size_t steps = 0;
    for (size_t i = 0; i < subs && r == GSHGEOM_OK; i++) {
      GSHGEOM_Flat_Subpath sp;
      r = gshgeom_flat_path_subpath(flat, i, &sp);
      if (r != GSHGEOM_OK) {
        break;
      }
      int ext = 0;
      for (size_t k = 1; k < sp.count; k++) {
        ext |= pts[sp.first + k].x != pts[sp.first].x || pts[sp.first + k].y != pts[sp.first].y;
      }
      const int to_aligned = aligned && sp.closed && ext;
      p.body = to_aligned ? o->aligned : o->body;
      p.width = stroke->width * (to_aligned ? 2.0 : 1.0);
      if (!ext && p.dirs) {
        /* a hatch's dot lies along its own direction (HX7), whatever the style's `degenerate` and caps say of the others; a piece that
         * only collapsed on the lattice lies along its own */
        double ux, uy;
        if (p.dotflag && !p.dotflag[i] && p.piece_dir) {
          ux = p.piece_dir[2 * i];
          uy = p.piece_dir[2 * i + 1];
        } else {
          ux = p.dir_i < p.dir_n ? p.dirs[2 * p.dir_i] : 1.0;
          uy = p.dir_i < p.dir_n ? p.dirs[2 * p.dir_i + 1] : 0.0;
          p.dir_i++;
        }
        if (!p.local) {
          const double vx = p.lin[0][0] * ux + p.lin[0][1] * uy, vy = p.lin[1][0] * ux + p.lin[1][1] * uy;
          ux = vx;
          uy = vy;
        }
        const double l = gmath_sqrt(ux * ux + uy * uy);
        p.dtx = l > 0.0 ? ux / l : 1.0;
        p.dty = l > 0.0 ? uy / l : 0.0;
      } else if (!ext) {
        double tx = 1.0, ty = 0.0; /* the map's image of +x under device width, +x itself under local */
        double dx = p.local ? 0.0 : p.lin[0][0], dy = p.local ? 0.0 : p.lin[1][0];
        int found = 0;
        for (size_t j = i; j-- > 0 && !found;) {
          if (ends[j].has) {
            dx = ends[j].lx;
            dy = ends[j].ly;
            found = 1;
          }
        }
        for (size_t j = i + 1; j < subs && !found; j++) {
          if (ends[j].has) {
            dx = ends[j].fx;
            dy = ends[j].fy;
            found = 1;
          }
        }
        if (found) {
          double ux = dx, uy = dy;
          if (p.local) {
            ux = p.inv[0][0] * dx + p.inv[0][1] * dy;
            uy = p.inv[1][0] * dx + p.inv[1][1] * dy;
          }
          const double l = gmath_sqrt(ux * ux + uy * uy);
          if (l > 0.0) {
            tx = ux / l;
            ty = uy / l;
          }
        } else if (!p.local) {
          const double l = gmath_sqrt(dx * dx + dy * dy);
          if (l > 0.0) {
            tx = dx / l;
            ty = dy / l;
          }
        }
        p.dtx = tx;
        p.dty = ty;
      }
      r = subpath(&p, &canon, stroke, pts + sp.first, sp.count, sp.closed, msubs ? &msubs[i] : NULL, &begin, &steps);
    }
  }
  gshgeom_flat_path_destroy(flat);
  for (size_t i = 0; msubs && i < msub_n; i++) {
    gcu_allocator_free(&p.alloc, msubs[i].xy);
  }
  gcu_allocator_free(&p.alloc, msubs);
  gcu_allocator_free(&p.alloc, ends);
  gcu_allocator_free(&p.alloc, p.dotflag);
  gcu_allocator_free(&p.alloc, p.piece_dir);
  gcu_allocator_free(&p.alloc, p.a.x);
  gcu_allocator_free(&p.alloc, p.a.y);
  gcu_allocator_free(&p.alloc, p.b.x);
  gcu_allocator_free(&p.alloc, p.b.y);
  gcu_allocator_free(&p.alloc, p.hx);
  gcu_allocator_free(&p.alloc, p.hy);
  gcu_allocator_free(&p.alloc, p.pts);
  gcu_allocator_free(&p.alloc, p.chain);
  gcu_allocator_free(&p.alloc, p.bq);
  gcu_allocator_free(&p.alloc, p.segs);
  gcu_allocator_free(&p.alloc, p.st);
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
