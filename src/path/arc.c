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
 * Exact arcs (arc.h): the record's validation, readback, the exact decomposition into
 * rational quadratics, lowering, and the three constructors.
 *
 * Every number a rule produces is one `double` operation in the written order (no
 * contraction, no `long double`, `sin`, `cos`, `tan`, `atan2` and `sqrt` from `math`).  The
 * decisions that go wrong on a rounded sum (is the radius too small; is the ellipse flat)
 * are made on exact sums of products (exact.h), and the radicand of the endpoint form's
 * centre is evaluated exactly and rounded once, which is what keeps the centre accurate
 * as the two ends approach a half ellipse.
 */

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/shape-geometry/arc.h>
#include <ghoti.io/shape-geometry/hash.h>
#include <ghoti.io/shape-geometry/orient.h>
#include <ghoti.io/shape-geometry/path.h>

#include <ghoti.io/math/basic.h>
#include <ghoti.io/math/elementary.h>

#include "../exact/exact.h"
#include "path_internal.h"

#include <string.h>

/** Radians per degree, per gradian and per turn, correctly rounded (as map.c's). */
#define RAD_PER_DEGREE 0x1.1df46a2529d39p-6
#define RAD_PER_GRADIAN 0x1.015bf9217271ap-6
#define RAD_PER_TURN 0x1.921fb54442d18p+2
/** The double nearest pi / 2 (a quarter turn in radians), pi and 2 pi. */
#define QUARTER_RADIANS 0x1.921fb54442d18p+0
#define PI_DOUBLE 0x1.921fb54442d18p+1
#define TWO_PI_DOUBLE 0x1.921fb54442d18p+2
/** The double just below 2 pi. */
#define TWO_PI_BELOW 0x1.921fb54442d17p+2

static int is_finite(double v) {
  return (v - v) == 0.0;
}

static double absf(double v) {
  return v < 0.0 ? 0.0 - v : v;
}

static double maxf(double a, double b) {
  return a > b ? a : b;
}

static int same_bits(double a, double b) {
  return gshgeom_f64_bits(a) == gshgeom_f64_bits(b);
}

static int same_xy(GSHGEOM_Point a, GSHGEOM_Point b) {
  return same_bits(a.x, b.x) && same_bits(a.y, b.y);
}

static GSHGEOM_Point pt(double x, double y) {
  GSHGEOM_Point p = {x, y, 0.0};
  return p;
}

/** A quarter turn in @p unit, and the radians per unit. */
static double quarter_of(GSHGEOM_Angle_Unit unit) {
  switch (unit) {
  case GSHGEOM_ANGLE_DEGREES:
    return 90.0;
  case GSHGEOM_ANGLE_GRADIANS:
    return 100.0;
  case GSHGEOM_ANGLE_TURNS:
    return 0.25;
  default:
    return QUARTER_RADIANS;
  }
}

static double per_unit_of(GSHGEOM_Angle_Unit unit) {
  switch (unit) {
  case GSHGEOM_ANGLE_DEGREES:
    return RAD_PER_DEGREE;
  case GSHGEOM_ANGLE_GRADIANS:
    return RAD_PER_GRADIAN;
  case GSHGEOM_ANGLE_TURNS:
    return RAD_PER_TURN;
  default:
    return 1.0;
  }
}

/** One full turn in the unit (four quarters: exact in every unit). */
static double turn_of(GSHGEOM_Angle_Unit unit) {
  return 4.0 * quarter_of(unit);
}

GSHGEOM_Result gshgeom_arc_point(const GSHGEOM_Arc_Form * form, double t, GSHGEOM_Point * out) {
  if (!form || !out || (unsigned)form->unit >= (unsigned)GSHGEOM_ANGLE_UNIT_COUNT) {
    return GSHGEOM_ERR_INVALID;
  }
  double c;
  double s;
  GSHGEOM_Result r = gshgeom_sincos(form->unit, t, &c, &s);
  if (r != GSHGEOM_OK) {
    return r;
  }
  double ix = (form->u.x * c) + (form->v.x * s);
  double iy = (form->u.y * c) + (form->v.y * s);
  out->x = form->center.x + ix;
  out->y = form->center.y + iy;
  out->z = 0.0;
  return GSHGEOM_OK;
}

static GSHGEOM_Arc_Form form_of(const GSHGEOM_Arc * a) {
  GSHGEOM_Arc_Form f;
  f.center = a->center;
  f.u = a->u;
  f.v = a->v;
  f.unit = a->unit;
  f.t0 = a->t0;
  f.sweep = a->sweep;
  f.flags = a->flags;
  return f;
}

/** The arc's scale: the largest coordinate magnitude, and `|u|_inf + |v|_inf`. */
static double scale_of(const GSHGEOM_Arc * a) {
  double s = maxf(maxf(absf(a->start.x), absf(a->start.y)), maxf(absf(a->end.x), absf(a->end.y)));
  s = maxf(s, maxf(absf(a->center.x), absf(a->center.y)));
  double uv = maxf(absf(a->u.x), absf(a->u.y)) + maxf(absf(a->v.x), absf(a->v.y));
  return maxf(s, uv);
}

/**
 * The decomposition of an accepted record must be representable: a control point of a piece (the meeting of two tangents) can lie past
 * the coordinate limit when the centre form is close to it, and an arc whose pieces cannot be made could not be read by anything that
 * walks a path by its pieces.  Such a record is a limit, as a coordinate past the limit is.
 */
static GSHGEOM_Result pieces_within_limit(const GSHGEOM_Arc * a) {
  GSHGEOM_Curve pieces[GSHGEOM_ARC_MAX_PIECES];
  int n = 0;
  GSHGEOM_Result r = gshgeom_arc_pieces(a, pieces, &n);
  return r == GSHGEOM_ERR_LIMIT ? r : GSHGEOM_OK;
}

/**
 * AR2's checks, in order, for an arc to be appended where the current point is @p current.
 * @p max_elements_room is nonzero if the path can take another element.
 */
static GSHGEOM_Result validate(const GSHGEOM_Arc * a, GSHGEOM_Point current, int room) {
  const double nums[] = {a->start.x, a->start.y, a->end.x, a->end.y, a->center.x, a->center.y, a->u.x,
      a->u.y, a->v.x, a->v.y, a->t0, a->sweep};
  for (size_t i = 0; i < sizeof nums / sizeof nums[0]; i++) {
    if (!is_finite(nums[i])) {
      return GSHGEOM_ERR_INVALID;
    }
  }
  if ((unsigned)a->unit >= (unsigned)GSHGEOM_ANGLE_UNIT_COUNT || (a->flags & ~GSHGEOM_ARC_CIRCULAR) != 0u) {
    return GSHGEOM_ERR_INVALID;
  }
  {
    /* the ellipse must have area: decided exactly (the sign of cross(u, v)) */
    GSHGEOM_Point z = pt(0.0, 0.0);
    int sign = 0;
    if (gshgeom_orient2d(z, pt(a->u.x, a->u.y), pt(a->v.x, a->v.y), &sign) != GSHGEOM_OK || sign == 0) {
      return GSHGEOM_ERR_INVALID;
    }
  }
  double turn = turn_of(a->unit);
  if (a->sweep == 0.0 || absf(a->sweep) > turn) {
    return GSHGEOM_ERR_INVALID;
  }
  for (size_t i = 0; i < 10; i++) {
    if (absf(nums[i]) > GSHGEOM_COORD_LIMIT) {
      return GSHGEOM_ERR_LIMIT;
    }
  }
  double qt = absf(a->t0) / quarter_of(a->unit);
  if (qt > GSHGEOM_MAP_MAX_QUARTERS) {
    return GSHGEOM_ERR_LIMIT;
  }
  if (!room) {
    return GSHGEOM_ERR_LIMIT;
  }
  /* AR4: the start is the current point, and the centre form reproduces both ends */
  if (!same_xy(a->start, current)) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Arc_Form f = form_of(a);
  double tol = GSHGEOM_ARC_AGREE * scale_of(a);
  GSHGEOM_Point p;
  GSHGEOM_Result r = gshgeom_arc_point(&f, a->t0, &p);
  if (r != GSHGEOM_OK) {
    return r;
  }
  if (absf(p.x - a->start.x) > tol || absf(p.y - a->start.y) > tol) {
    return GSHGEOM_ERR_INVALID;
  }
  if (absf(a->sweep) == turn) {
    return same_xy(a->start, a->end) ? pieces_within_limit(a) : GSHGEOM_ERR_INVALID;
  }
  r = gshgeom_arc_point(&f, a->t0 + a->sweep, &p);
  if (r != GSHGEOM_OK) {
    return r;
  }
  if (absf(p.x - a->end.x) > tol || absf(p.y - a->end.y) > tol) {
    return GSHGEOM_ERR_INVALID;
  }
  return pieces_within_limit(a);
}

static GSHGEOM_Result append_checked(GSHGEOM_Path * path, const GSHGEOM_Arc * a) {
  if (path->dim != 2 || !path->has_current) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Result r = validate(a, path->current, path->el_n < path->max_elements);
  if (r != GSHGEOM_OK) {
    return r;
  }
  GSHGEOM_Point pts[4] = {pt(a->center.x, a->center.y), pt(a->u.x, a->u.y), pt(a->v.x, a->v.y),
      pt(a->end.x, a->end.y)};
  const double w[2] = {a->t0, a->sweep};
  return gshgeom_path_append_raw(path, GSHGEOM_VERB_ARC, (unsigned)a->unit, a->flags, pts, w);
}

GSHGEOM_Result gshgeom_path_append_arc(GSHGEOM_Path * path, const GSHGEOM_Arc * arc) {
  if (!path || !arc) {
    return GSHGEOM_ERR_INVALID;
  }
  return append_checked(path, arc);
}

GSHGEOM_Result gshgeom_path_arc(const GSHGEOM_Path * path, size_t index, GSHGEOM_Arc * out) {
  if (!path || !out || index >= path->el_n || path->el[index].verb != GSHGEOM_VERB_ARC) {
    return GSHGEOM_ERR_INVALID;
  }
  const Elem * e = &path->el[index];
  const GSHGEOM_Point * q = &path->pt[e->point];
  GSHGEOM_Arc a;
  a.start = gshgeom_path_element_start(path, index);
  a.end = q[3];
  a.center = q[0];
  a.u = q[1];
  a.v = q[2];
  a.unit = (GSHGEOM_Angle_Unit)e->unit;
  a.t0 = path->wt[e->weight];
  a.sweep = path->wt[e->weight + 1];
  a.flags = e->flags;
  *out = a;
  return GSHGEOM_OK;
}

/* ---------------------------------------------------------------------- */
/* The exact decomposition (AR8)                                           */
/* ---------------------------------------------------------------------- */

/** t0 + (k sweep) / n, the product rounded once, then the quotient. */
static double param_at(const GSHGEOM_Arc * a, int k, int n) {
  if (k == 0) {
    return a->t0;
  }
  return a->t0 + ((double)k * a->sweep) / (double)n;
}

GSHGEOM_Result gshgeom_arc_pieces(
    const GSHGEOM_Arc * arc, GSHGEOM_Curve out[GSHGEOM_ARC_MAX_PIECES], int * np) {
  if (!arc || !out || !np || (unsigned)arc->unit >= (unsigned)GSHGEOM_ANGLE_UNIT_COUNT || !is_finite(arc->sweep)
      || !is_finite(arc->t0) || arc->sweep == 0.0) {
    return GSHGEOM_ERR_INVALID;
  }
  const double q = quarter_of(arc->unit);
  const double a = absf(arc->sweep);
  /* how many quarter turns the sweep spans, compared exactly: a - 2q is exact for a in [2q, 4q] */
  int n;
  if (a <= q) {
    n = 1;
  } else if (a <= 2.0 * q) {
    n = 2;
  } else if ((a - 2.0 * q) <= q) {
    n = 3;
  } else {
    n = 4;
  }
  const double sgn = arc->sweep < 0.0 ? -1.0 : 1.0;
  const double hu = a / (2.0 * (double)n);
  const int quarter_piece = (hu == 0.5 * q);
  const double h = hu * per_unit_of(arc->unit);
  const double w = quarter_piece ? gmath_sqrt(0.5) : gmath_cos(h);
  const double th = quarter_piece ? 1.0 : gmath_tan(h);
  GSHGEOM_Arc_Form f = form_of(arc);
  double cs_a = 0.0;
  double sn_a = 0.0;
  GSHGEOM_Result r = gshgeom_sincos(arc->unit, param_at(arc, 0, n), &cs_a, &sn_a);
  if (r != GSHGEOM_OK) {
    return r;
  }
  GSHGEOM_Point pa = arc->start;
  for (int k = 0; k < n; k++) {
    double cs_b;
    double sn_b;
    r = gshgeom_sincos(arc->unit, param_at(arc, k + 1, n), &cs_b, &sn_b);
    if (r != GSHGEOM_OK) {
      return r;
    }
    GSHGEOM_Point pb;
    if (k == n - 1) {
      pb = arc->end;
    } else {
      pb.x = f.center.x + ((f.u.x * cs_b) + (f.v.x * sn_b));
      pb.y = f.center.y + ((f.u.y * cs_b) + (f.v.y * sn_b));
      pb.z = 0.0;
    }
    /* the tangent direction at each end: sign (-u sin t + v cos t) */
    double tax = sgn * ((f.v.x * cs_a) - (f.u.x * sn_a));
    double tay = sgn * ((f.v.y * cs_a) - (f.u.y * sn_a));
    double tbx = sgn * ((f.v.x * cs_b) - (f.u.x * sn_b));
    double tby = sgn * ((f.v.y * cs_b) - (f.u.y * sn_b));
    GSHGEOM_Point ctl;
    ctl.z = 0.0;
    if (tay == 0.0 && tbx == 0.0) {
      ctl.x = pb.x;
      ctl.y = pa.y;
    } else if (tax == 0.0 && tby == 0.0) {
      ctl.x = pa.x;
      ctl.y = pb.y;
    } else {
      ctl.x = pa.x + (th * tax);
      ctl.y = pa.y + (th * tay);
    }
    const double wts[3] = {1.0, w, 1.0};
    r = gshgeom_curve_rational_quadratic(2, pa, ctl, pb, wts, &out[k]);
    if (r != GSHGEOM_OK) {
      return r;
    }
    pa = pb;
    cs_a = cs_b;
    sn_a = sn_b;
  }
  *np = n;
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_path_lower_arcs(const GSHGEOM_Path * path, const GSHGEOM_Allocator * alloc,
    const GSHGEOM_Limits * limits, GSHGEOM_Path ** out) {
  if (!path || !out) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Path * dst;
  GSHGEOM_Result r = gshgeom_path_create(alloc, limits, path->dim, &dst);
  if (r != GSHGEOM_OK) {
    return r;
  }
  for (size_t i = 0; i < path->el_n && r == GSHGEOM_OK; i++) {
    const Elem * e = &path->el[i];
    if (e->verb == GSHGEOM_VERB_ARC) {
      GSHGEOM_Arc arc;
      GSHGEOM_Curve pieces[GSHGEOM_ARC_MAX_PIECES];
      int n = 0;
      r = gshgeom_path_arc(path, i, &arc);
      if (r == GSHGEOM_OK) {
        r = gshgeom_arc_pieces(&arc, pieces, &n);
      }
      for (int k = 0; k < n && r == GSHGEOM_OK; k++) {
        r = gshgeom_path_rational_quadratic_to(dst, pieces[k].p[1], pieces[k].p[2], pieces[k].w);
      }
    } else {
      r = gshgeom_path_append_raw(dst, e->verb, e->unit, e->flags, gshgeom_path_elem_points(path, e),
          gshgeom_path_elem_weights(path, e));
      if (r == GSHGEOM_OK && e->verb == GSHGEOM_VERB_MOVE) {
        dst->has_current = 1;
        dst->start = path->pt[e->point];
      }
    }
  }
  if (r != GSHGEOM_OK) {
    gshgeom_path_destroy(dst);
    return r;
  }
  *out = dst;
  return GSHGEOM_OK;
}

/* ---------------------------------------------------------------------- */
/* Maps and the circular flag (AR11)                                       */
/* ---------------------------------------------------------------------- */

GSHGEOM_Result gshgeom_arc_map(const GSHGEOM_Map * map, const GSHGEOM_Arc * arc, GSHGEOM_Arc * out) {
  if (!map || !arc || !out || map->dim != 2) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Result r = gshgeom_map_check(map);
  if (r != GSHGEOM_OK) {
    return r;
  }
  if ((map->kind & GSHGEOM_MAP_SINGULAR) != 0u) {
    return GSHGEOM_ERR_INVALID; /* the image is a segment: a consumer that needs it maps the pieces */
  }
  GSHGEOM_Arc a = *arc;
  r = gshgeom_map_apply_point(map, arc->start, &a.start);
  if (r == GSHGEOM_OK) {
    r = gshgeom_map_apply_point(map, arc->end, &a.end);
  }
  if (r == GSHGEOM_OK) {
    r = gshgeom_map_apply_point(map, arc->center, &a.center);
  }
  if (r == GSHGEOM_OK) {
    r = gshgeom_map_apply_vector(map, arc->u, &a.u);
  }
  if (r == GSHGEOM_OK) {
    r = gshgeom_map_apply_vector(map, arc->v, &a.v);
  }
  if (r != GSHGEOM_OK) {
    return r;
  }
  const unsigned not_similar = GSHGEOM_MAP_SCALES_NONUNIFORM | GSHGEOM_MAP_SKEWS | GSHGEOM_MAP_SINGULAR | GSHGEOM_MAP_GENERAL;
  if ((map->kind & not_similar) != 0u) {
    a.flags &= ~GSHGEOM_ARC_CIRCULAR;
  }
  r = validate(&a, a.start, 1);
  if (r != GSHGEOM_OK) {
    return r;
  }
  *out = a;
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_arc_circle(const GSHGEOM_Arc * arc, GSHGEOM_Point * centre, double * radius, int * ccw) {
  if (!arc || !centre || !radius || (arc->flags & GSHGEOM_ARC_CIRCULAR) == 0u) {
    return GSHGEOM_ERR_INVALID;
  }
  int sign = 0;
  GSHGEOM_Point z = pt(0.0, 0.0);
  if (gshgeom_orient2d(z, pt(arc->u.x, arc->u.y), pt(arc->v.x, arc->v.y), &sign) != GSHGEOM_OK || sign == 0) {
    return GSHGEOM_ERR_INVALID;
  }
  *centre = pt(arc->center.x, arc->center.y);
  *radius = gmath_hypot(arc->u.x, arc->u.y);
  if (ccw) {
    *ccw = ((arc->sweep < 0.0) ? -1 : 1) * sign > 0;
  }
  return GSHGEOM_OK;
}

/* ---------------------------------------------------------------------- */
/* Constructors                                                            */
/* ---------------------------------------------------------------------- */

static GSHGEOM_Result begin_check(const GSHGEOM_Path * path) {
  if (!path || path->dim != 2) {
    return GSHGEOM_ERR_INVALID;
  }
  return GSHGEOM_OK;
}

/** A linear form `sum s_i p_i q_i` of up to four products, for the exact squares. */
typedef struct Lin {
  int n;
  int neg[4];
  double p[4];
  double q[4];
} Lin;

static void lin_add(Lin * l, int neg, double p, double q) {
  l->neg[l->n] = neg;
  l->p[l->n] = p;
  l->q[l->n] = q;
  l->n++;
}

/** Add `(-1)^neg * k * L^2` to the sum (k a power of two, 1 or 4). */
static void add_square(GSHGEOM_ExactSum * sum, int neg, double k, const Lin * l) {
  for (int i = 0; i < l->n; i++) {
    for (int j = 0; j < l->n; j++) {
      const double f[5] = {l->p[i], l->q[i], l->p[j], l->q[j], k};
      gshgeom_exact_sum_add(sum, neg ^ l->neg[i] ^ l->neg[j], 5, f);
    }
  }
}

static double lin_value(const Lin * l) {
  GSHGEOM_ExactSum s;
  gshgeom_exact_sum_init(&s);
  for (int i = 0; i < l->n; i++) {
    const double f[2] = {l->p[i], l->q[i]};
    gshgeom_exact_sum_add(&s, l->neg[i], 2, f);
  }
  return gshgeom_exact_sum_to_double(&s);
}

GSHGEOM_Result gshgeom_path_arc_endpoint_to(GSHGEOM_Path * path, GSHGEOM_Point end, double rx, double ry,
    double phi, GSHGEOM_Angle_Unit unit, int large, int sweep, GSHGEOM_Arc_Kind * kind) {
  GSHGEOM_Result r = begin_check(path);
  if (r != GSHGEOM_OK || !path->has_current || (large != 0 && large != 1) || (sweep != 0 && sweep != 1)
      || (unsigned)unit >= (unsigned)GSHGEOM_ANGLE_UNIT_COUNT || !is_finite(end.x) || !is_finite(end.y)
      || !is_finite(rx) || !is_finite(ry) || !is_finite(phi)) {
    return GSHGEOM_ERR_INVALID;
  }
  const GSHGEOM_Point S = path->current;
  const GSHGEOM_Point E = pt(end.x, end.y);
  if (same_xy(S, E)) {
    if (kind) {
      *kind = GSHGEOM_ARC_KIND_NOTHING;
    }
    return GSHGEOM_OK;
  }
  if (rx == 0.0 || ry == 0.0) {
    r = gshgeom_path_line_to(path, E);
    if (r == GSHGEOM_OK && kind) {
      *kind = GSHGEOM_ARC_KIND_LINE;
    }
    return r;
  }
  rx = absf(rx);
  ry = absf(ry);
  double cp;
  double sp;
  r = gshgeom_sincos(unit, phi, &cp, &sp);
  if (r != GSHGEOM_OK) {
    return r;
  }
  GSHGEOM_Point u = pt(rx * cp, rx * sp);
  GSHGEOM_Point v = pt(ry * (0.0 - sp), ry * cp);
  if (absf(u.x) > GSHGEOM_COORD_LIMIT || absf(u.y) > GSHGEOM_COORD_LIMIT || absf(v.x) > GSHGEOM_COORD_LIMIT
      || absf(v.y) > GSHGEOM_COORD_LIMIT || !is_finite(u.x) || !is_finite(u.y) || !is_finite(v.x) || !is_finite(v.y)
      || absf(E.x) > GSHGEOM_COORD_LIMIT || absf(E.y) > GSHGEOM_COORD_LIMIT) {
    return GSHGEOM_ERR_LIMIT;
  }
  /* d = E - S stays an expansion: a = v.y dx - v.x dy, b = u.x dy - u.y dx, D = u.x v.y - u.y v.x */
  Lin a = {0}, b = {0}, d = {0};
  lin_add(&a, 0, v.y, E.x);
  lin_add(&a, 1, v.y, S.x);
  lin_add(&a, 1, v.x, E.y);
  lin_add(&a, 0, v.x, S.y);
  lin_add(&b, 0, u.x, E.y);
  lin_add(&b, 1, u.x, S.y);
  lin_add(&b, 1, u.y, E.x);
  lin_add(&b, 0, u.y, S.x);
  lin_add(&d, 0, u.x, v.y);
  lin_add(&d, 1, u.y, v.x);
  if (lin_value(&d) <= 0.0) {
    GSHGEOM_Point z = pt(0.0, 0.0);
    int sign = 0;
    if (gshgeom_orient2d(z, u, v, &sign) != GSHGEOM_OK || sign <= 0) {
      return GSHGEOM_ERR_INVALID;
    }
  }
  /* N = 4 D^2 - a^2 - b^2 (positive iff the radii fit), AB = a^2 + b^2, DEN = 4 D^2 */
  GSHGEOM_ExactSum n_sum, ab_sum, den_sum;
  gshgeom_exact_sum_init(&n_sum);
  gshgeom_exact_sum_init(&ab_sum);
  gshgeom_exact_sum_init(&den_sum);
  add_square(&n_sum, 0, 4.0, &d);
  add_square(&n_sum, 1, 1.0, &a);
  add_square(&n_sum, 1, 1.0, &b);
  add_square(&ab_sum, 0, 1.0, &a);
  add_square(&ab_sum, 0, 1.0, &b);
  add_square(&den_sum, 0, 4.0, &d);
  const int fits = gshgeom_exact_sum_sign(&n_sum) > 0;
  const double ad = lin_value(&a);
  const double bd = lin_value(&b);
  const double abd = gshgeom_exact_sum_to_double(&ab_sum);
  const double dend = gshgeom_exact_sum_to_double(&den_sum);
  GSHGEOM_Arc arc;
  memset(&arc, 0, sizeof arc);
  arc.start = pt(S.x, S.y);
  arc.end = E;
  arc.unit = GSHGEOM_ANGLE_RADIANS;
  arc.flags = rx == ry ? GSHGEOM_ARC_CIRCULAR : 0u;
  const GSHGEOM_Point mid = pt((S.x * 0.5) + (E.x * 0.5), (S.y * 0.5) + (E.y * 0.5));
  if (!fits) {
    /* F.6.6: scale the radii by sqrt(Lambda); the centre is the midpoint, the sweep half a turn */
    const double lambda = abd / dend;
    const double sc = gmath_sqrt(lambda);
    arc.u = pt(u.x * sc, u.y * sc);
    arc.v = pt(v.x * sc, v.y * sc);
    arc.center = mid;
    arc.t0 = gmath_atan2(0.0 - bd, 0.0 - ad);
    arc.sweep = sweep ? PI_DOUBLE : 0.0 - PI_DOUBLE;
  } else {
    const double nd = gshgeom_exact_sum_to_double(&n_sum);
    const double dd = lin_value(&d);
    double kk = gmath_sqrt(nd / abd);
    if (large == sweep) {
      kk = 0.0 - kk; /* F.6.5: the plus sign is chosen if the flags differ */
    }
    const double g = kk / (2.0 * dd);
    const double cx = g * (0.0 - bd); /* the centre, in the unit circle's frame, from the midpoint */
    const double cy = g * ad;
    arc.u = u;
    arc.v = v;
    arc.center = pt(mid.x + ((u.x * cx) + (v.x * cy)), mid.y + ((u.y * cx) + (v.y * cy)));
    const double sx = (kk * bd) - ad;
    const double sy = (0.0 - bd) - (kk * ad);
    const double t0 = gmath_atan2(sy, sx);
    /* The angle between the two ends seen from the centre, in the unit circle's frame: the chord is 2 sqrt(Lambda), so
     * the half angle is atan2(sqrt(a^2 + b^2), sqrt(N)) from the two exact sums, accurate however small the chord is
     * (a difference of two atan2 results cannot show a sweep below the spacing of t0).  F.6.5's large-arc flag
     * chooses the angle or its complement, the sweep flag the sign. */
    const double half = gmath_atan2(gmath_sqrt(abd), gmath_sqrt(nd));
    double mag = large ? TWO_PI_DOUBLE - (2.0 * half) : 2.0 * half;
    if (mag >= TWO_PI_DOUBLE) {
      /* The ends differ, so this is not the full turn that closes: keep the sweep one step below it (the
       * chord this leaves open is a few 1e-16 of the radius, far inside the agreement bound). */
      mag = TWO_PI_BELOW;
    }
    const double dth = sweep ? mag : 0.0 - mag;
    arc.t0 = t0;
    arc.sweep = dth;
  }
  r = append_checked(path, &arc);
  if (r == GSHGEOM_OK && kind) {
    *kind = GSHGEOM_ARC_KIND_ARC;
  }
  return r;
}

/** Add the arc for @p form from the current point (or beginning at its start), atomically. */
static GSHGEOM_Result add_form(GSHGEOM_Path * path, const GSHGEOM_Arc_Form * form, const GSHGEOM_Point * end) {
  GSHGEOM_Result r = begin_check(path);
  if (r != GSHGEOM_OK || !form) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Point p0;
  r = gshgeom_arc_point(form, form->t0, &p0);
  if (r != GSHGEOM_OK) {
    return r;
  }
  GSHGEOM_Arc a;
  memset(&a, 0, sizeof a);
  a.start = path->has_current ? pt(path->current.x, path->current.y) : p0;
  a.center = form->center;
  a.center.z = 0.0;
  a.u = form->u;
  a.u.z = 0.0;
  a.v = form->v;
  a.v.z = 0.0;
  a.unit = form->unit;
  a.t0 = form->t0;
  a.sweep = form->sweep;
  a.flags = form->flags;
  if (end) {
    a.end = pt(end->x, end->y);
  } else if (is_finite(form->sweep) && absf(form->sweep) == turn_of(form->unit)) {
    a.end = a.start;
  } else {
    r = gshgeom_arc_point(form, form->t0 + form->sweep, &a.end);
    if (r != GSHGEOM_OK) {
      return r;
    }
  }
  PathMark mark;
  path_mark(path, &mark);
  if (!path->has_current) {
    /* a subpath begins at the start; validate first so a bad arc leaves no move behind */
    r = validate(&a, a.start, path->el_n + 1 < path->max_elements);
    if (r != GSHGEOM_OK) {
      return r;
    }
    r = gshgeom_path_move_to(path, a.start);
    if (r != GSHGEOM_OK) {
      return r;
    }
  }
  r = append_checked(path, &a);
  if (r != GSHGEOM_OK) {
    path_rollback(path, &mark);
  }
  return r;
}

GSHGEOM_Result gshgeom_path_arc_conjugate_to(
    GSHGEOM_Path * path, const GSHGEOM_Arc_Form * form, const GSHGEOM_Point * end) {
  return add_form(path, form, end);
}

GSHGEOM_Result gshgeom_path_arc_center_to(GSHGEOM_Path * path, GSHGEOM_Point center, double rx, double ry,
    double phi, GSHGEOM_Angle_Unit unit, double t0, double sweep, const GSHGEOM_Point * end) {
  GSHGEOM_Result r = begin_check(path);
  if (r != GSHGEOM_OK || !is_finite(rx) || !is_finite(ry) || !(rx > 0.0) || !(ry > 0.0) || !is_finite(phi)
      || (unsigned)unit >= (unsigned)GSHGEOM_ANGLE_UNIT_COUNT) {
    return GSHGEOM_ERR_INVALID;
  }
  double cp;
  double sp;
  r = gshgeom_sincos(unit, phi, &cp, &sp);
  if (r != GSHGEOM_OK) {
    return r;
  }
  GSHGEOM_Arc_Form f;
  f.center = center;
  f.u = pt(rx * cp, rx * sp);
  f.v = pt(ry * (0.0 - sp), ry * cp);
  f.unit = unit;
  f.t0 = t0;
  f.sweep = sweep;
  f.flags = rx == ry ? GSHGEOM_ARC_CIRCULAR : 0u;
  return add_form(path, &f, end);
}

GSHGEOM_Result gshgeom_path_arc_tangent_to(GSHGEOM_Path * path, GSHGEOM_Point p1, GSHGEOM_Point p2, double r) {
  GSHGEOM_Result res = begin_check(path);
  if (res != GSHGEOM_OK || !path->has_current || !is_finite(r) || r < 0.0 || !is_finite(p1.x) || !is_finite(p1.y)
      || !is_finite(p2.x) || !is_finite(p2.y)) {
    return GSHGEOM_ERR_INVALID;
  }
  const GSHGEOM_Point p0 = pt(path->current.x, path->current.y);
  p1.z = 0.0;
  p2.z = 0.0;
  int sign = 0;
  res = gshgeom_orient2d(p0, p1, p2, &sign);
  if (res != GSHGEOM_OK) {
    return res;
  }
  if (same_xy(p0, p1) || same_xy(p1, p2) || r == 0.0 || sign == 0) {
    return gshgeom_path_line_to(path, p1);
  }
  const double l0 = gmath_hypot(p0.x - p1.x, p0.y - p1.y);
  const double l2 = gmath_hypot(p2.x - p1.x, p2.y - p1.y);
  const GSHGEOM_Point e0 = pt((p0.x - p1.x) / l0, (p0.y - p1.y) / l0);
  const GSHGEOM_Point e2 = pt((p2.x - p1.x) / l2, (p2.y - p1.y) / l2);
  const double cosq = (e0.x * e2.x) + (e0.y * e2.y);
  const double crs = (e0.x * e2.y) - (e0.y * e2.x);
  const double sinq = crs < 0.0 ? 0.0 - crs : crs;
  const double dist = (r * (1.0 + cosq)) / sinq;
  const GSHGEOM_Point t1 = pt(p1.x + (dist * e0.x), p1.y + (dist * e0.y));
  const GSHGEOM_Point t2 = pt(p1.x + (dist * e2.x), p1.y + (dist * e2.y));
  const double sg = sign > 0 ? 1.0 : -1.0;
  const GSHGEOM_Point c = pt(t1.x + ((r * sg) * e0.y), t1.y + ((r * sg) * (0.0 - e0.x)));
  const GSHGEOM_Point ra = pt(t1.x - c.x, t1.y - c.y);
  const GSHGEOM_Point rb = pt(t2.x - c.x, t2.y - c.y);
  GSHGEOM_Arc a;
  memset(&a, 0, sizeof a);
  a.start = t1;
  a.end = t2;
  a.center = c;
  a.u = pt(r, 0.0);
  a.v = pt(0.0, r);
  a.unit = GSHGEOM_ANGLE_RADIANS;
  a.t0 = gmath_atan2(ra.y, ra.x);
  a.sweep = gmath_atan2((ra.x * rb.y) - (ra.y * rb.x), (ra.x * rb.x) + (ra.y * rb.y));
  a.flags = GSHGEOM_ARC_CIRCULAR;
  PathMark mark;
  path_mark(path, &mark);
  if (!same_xy(t1, p0)) {
    res = gshgeom_path_line_to(path, t1);
    if (res != GSHGEOM_OK) {
      return res;
    }
  }
  res = append_checked(path, &a);
  if (res != GSHGEOM_OK) {
    path_rollback(path, &mark);
  }
  return res;
}

GSHGEOM_Result gshgeom_path_arc_bulge_to(GSHGEOM_Path * path, GSHGEOM_Point q, double bulge) {
  GSHGEOM_Result res = begin_check(path);
  if (res != GSHGEOM_OK || !path->has_current || !is_finite(bulge) || !is_finite(q.x) || !is_finite(q.y)) {
    return GSHGEOM_ERR_INVALID;
  }
  const GSHGEOM_Point s = pt(path->current.x, path->current.y);
  q.z = 0.0;
  if (same_xy(s, q)) {
    return GSHGEOM_OK;
  }
  if (bulge == 0.0) {
    return gshgeom_path_line_to(path, q);
  }
  const double sweep = 4.0 * gmath_atan2(bulge, 1.0);
  const double half = sweep * 0.5;
  const double dx = q.x - s.x;
  const double dy = q.y - s.y;
  const double len = gmath_hypot(dx, dy);
  const double radius = (len * 0.5) / absf(gmath_sin(half));
  const double off = (len * 0.5) / gmath_tan(half);
  const GSHGEOM_Point nrm = pt((0.0 - dy) / len, dx / len); /* the left normal of the chord */
  const GSHGEOM_Point c = pt(((s.x * 0.5) + (q.x * 0.5)) + (off * nrm.x), ((s.y * 0.5) + (q.y * 0.5)) + (off * nrm.y));
  GSHGEOM_Arc a;
  memset(&a, 0, sizeof a);
  a.start = s;
  a.end = q;
  a.center = c;
  a.u = pt(radius, 0.0);
  a.v = pt(0.0, radius);
  a.unit = GSHGEOM_ANGLE_RADIANS;
  a.t0 = gmath_atan2(s.y - c.y, s.x - c.x);
  a.sweep = sweep;
  a.flags = GSHGEOM_ARC_CIRCULAR;
  return append_checked(path, &a);
}
