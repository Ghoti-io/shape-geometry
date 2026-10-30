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
 * The named primitives (shapes.h): the rectangle, the rounded rectangle and the ellipse, with the exact overlap
 * predicate.
 *
 * Every coordinate is one exact sum rounded once (exact.h), so it is the correctly rounded value of the real expression named.
 * The mirror of shapes.h is one exchange: a point is a distance from the left or the right edge, and from the top or the bottom
 * edge, and a flip exchanges the edge.
 */

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/shape-geometry/shapes.h>

#include "../exact/exact.h"
#include "path_internal.h"

#include <stddef.h>

static double absd(double v) {
  return v < 0.0 ? -v : v;
}

static int is_finite_num(double v) {
  return v - v == 0.0;
}

/** The correctly rounded `a + b`, `+0` when it is exactly zero. */
static double sum2(double a, double b) {
  GSHGEOM_ExactSum s;
  gshgeom_exact_sum_init(&s);
  gshgeom_exact_sum_add(&s, 0, 1, &a);
  gshgeom_exact_sum_add(&s, 0, 1, &b);
  return gshgeom_exact_sum_to_double(&s);
}

/** The correctly rounded `a + b - c`. */
static double sum3(double a, double b, double c) {
  GSHGEOM_ExactSum s;
  gshgeom_exact_sum_init(&s);
  gshgeom_exact_sum_add(&s, 0, 1, &a);
  gshgeom_exact_sum_add(&s, 0, 1, &b);
  gshgeom_exact_sum_add(&s, 1, 1, &c);
  return gshgeom_exact_sum_to_double(&s);
}

/** The sign of `a + b - c` in the reals. */
static int over(double a, double b, double c) {
  GSHGEOM_ExactSum s;
  gshgeom_exact_sum_init(&s);
  gshgeom_exact_sum_add(&s, 0, 1, &a);
  gshgeom_exact_sum_add(&s, 0, 1, &b);
  gshgeom_exact_sum_add(&s, 1, 1, &c);
  return gshgeom_exact_sum_sign(&s) > 0;
}

/** The distance @p a from the left edge (`x + a`; a zero distance is `x`, bit for bit). */
static double from_near(double x, double a) {
  return a == 0.0 ? x : sum2(x, a);
}

/** The distance @p a from the far edge (`x + w - a`). */
static double from_far(double x, double w, double a) {
  return sum3(x, w, a);
}

/** The coordinate a distance @p a from the near or the far edge of the box, the edge exchanged when @p flip is set. */
static double edge(double x, double w, int far_edge, int flip, double a) {
  return (far_edge ^ flip) ? from_far(x, w, a) : from_near(x, a);
}

static GSHGEOM_Point pt(double x, double y) {
  GSHGEOM_Point p = {x, y, 0.0};
  return p;
}

static int out_of_range(double v) {
  return !(absd(v) <= GSHGEOM_COORD_LIMIT);
}

static GSHGEOM_Corner canonical(GSHGEOM_Corner c) {
  if (c.rx == 0.0 || c.ry == 0.0) {
    c.rx = 0.0;
    c.ry = 0.0;
  }
  return c;
}

GSHGEOM_Result gshgeom_rounded_rect_check(double w, double h, const GSHGEOM_Corner r[4]) {
  if (!r || !is_finite_num(w) || !is_finite_num(h)) {
    return GSHGEOM_ERR_INVALID;
  }
  for (unsigned i = 0; i < 4; i++) {
    if (!is_finite_num(r[i].rx) || !is_finite_num(r[i].ry)) {
      return GSHGEOM_ERR_INVALID;
    }
  }
  if (w < 0.0 || h < 0.0) {
    return GSHGEOM_ERR_INVALID;
  }
  for (unsigned i = 0; i < 4; i++) {
    if (r[i].rx < 0.0 || r[i].ry < 0.0) {
      return GSHGEOM_ERR_INVALID;
    }
  }
  GSHGEOM_Corner c[4];
  for (unsigned i = 0; i < 4; i++) {
    c[i] = canonical(r[i]);
  }
  if (over(c[0].rx, c[1].rx, w) || over(c[2].rx, c[3].rx, w) || over(c[3].ry, c[0].ry, h) || over(c[1].ry, c[2].ry, h)) {
    return GSHGEOM_ERR_INVALID;
  }
  return GSHGEOM_OK;
}

static GSHGEOM_Result begin(GSHGEOM_Path * path) {
  return path && gshgeom_path_dimension(path) == 2 ? GSHGEOM_OK : GSHGEOM_ERR_INVALID;
}

static GSHGEOM_Result rect_points(double x, double y, double w, double h, int fx, int fy, GSHGEOM_Point p[4]) {
  static const int far_x[4] = {0, 1, 1, 0};
  static const int far_y[4] = {0, 0, 1, 1};
  for (unsigned i = 0; i < 4; i++) {
    p[i] = pt(edge(x, w, far_x[i], fx, 0.0), edge(y, h, far_y[i], fy, 0.0));
    if (out_of_range(p[i].x) || out_of_range(p[i].y)) {
      return GSHGEOM_ERR_LIMIT;
    }
  }
  return GSHGEOM_OK;
}

static GSHGEOM_Result emit_rect(GSHGEOM_Path * path, const GSHGEOM_Point p[4]) {
  GSHGEOM_Result r = gshgeom_path_move_to(path, p[0]);
  for (unsigned i = 1; r == GSHGEOM_OK && i < 4; i++) {
    r = gshgeom_path_line_to(path, p[i]);
  }
  return r == GSHGEOM_OK ? gshgeom_path_close(path) : r;
}

static GSHGEOM_Result add_rect(
    GSHGEOM_Path * path, double x, double y, double w, double h, int flip_x, int flip_y) {
  GSHGEOM_Result r = begin(path);
  if (r != GSHGEOM_OK || (flip_x != 0 && flip_x != 1) || (flip_y != 0 && flip_y != 1)) {
    return GSHGEOM_ERR_INVALID;
  }
  if (!is_finite_num(x) || !is_finite_num(y) || !is_finite_num(w) || !is_finite_num(h) || w < 0.0 || h < 0.0) {
    return GSHGEOM_ERR_INVALID;
  }
  if (out_of_range(w) || out_of_range(h)) {
    return GSHGEOM_ERR_LIMIT;
  }
  GSHGEOM_Point p[4];
  r = rect_points(x, y, w, h, flip_x, flip_y, p);
  return r == GSHGEOM_OK ? emit_rect(path, p) : r;
}

static GSHGEOM_Result quarter(GSHGEOM_Path * path, GSHGEOM_Point start, GSHGEOM_Point end, GSHGEOM_Point centre,
    GSHGEOM_Point u, GSHGEOM_Point v, double t0) {
  GSHGEOM_Arc a;
  a.start = start;
  a.end = end;
  a.center = centre;
  a.u = u;
  a.v = v;
  a.unit = GSHGEOM_ANGLE_DEGREES;
  a.t0 = t0;
  a.sweep = 90.0;
  a.flags = absd(u.x) == absd(v.y) ? GSHGEOM_ARC_CIRCULAR : 0u;
  return gshgeom_path_append_arc(path, &a);
}

static GSHGEOM_Result add_rounded_rect(
    GSHGEOM_Path * path, double x, double y, double w, double h, const GSHGEOM_Corner r[4]) {
  if (begin(path) != GSHGEOM_OK || gshgeom_rounded_rect_check(w, h, r) != GSHGEOM_OK || !is_finite_num(x) || !is_finite_num(y)) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Corner c[4];
  int round[4];
  int any = 0;
  for (unsigned i = 0; i < 4; i++) {
    c[i] = canonical(r[i]);
    round[i] = c[i].rx != 0.0;
    any |= round[i];
  }
  if (out_of_range(w) || out_of_range(h)) {
    return GSHGEOM_ERR_LIMIT;
  }
  if (!any) {
    GSHGEOM_Point p[4];
    GSHGEOM_Result e = rect_points(x, y, w, h, 0, 0, p);
    return e == GSHGEOM_OK ? emit_rect(path, p) : e;
  }
  /* The corners in the order top-left, top-right, bottom-right, bottom-left.  A radius is at most the side it lies on (the overlap
   * rule), so the limit on the width and the height bounds them. */
  const double rtl_x = c[0].rx, rtl_y = c[0].ry, rtr_x = c[1].rx, rtr_y = c[1].ry;
  const double rbr_x = c[2].rx, rbr_y = c[2].ry, rbl_x = c[3].rx, rbl_y = c[3].ry;
  const GSHGEOM_Point start = pt(from_near(x, rtl_x), y);
  const GSHGEOM_Point top_end = pt(from_far(x, w, rtr_x), y);
  const GSHGEOM_Point tr_end = pt(from_far(x, w, 0.0), from_near(y, rtr_y));
  const GSHGEOM_Point tr_c = pt(top_end.x, tr_end.y);
  const GSHGEOM_Point right_end = pt(tr_end.x, from_far(y, h, rbr_y));
  const GSHGEOM_Point br_end = pt(from_far(x, w, rbr_x), from_far(y, h, 0.0));
  const GSHGEOM_Point br_c = pt(br_end.x, right_end.y);
  const GSHGEOM_Point bottom_end = pt(from_near(x, rbl_x), br_end.y);
  const GSHGEOM_Point bl_end = pt(x, from_far(y, h, rbl_y));
  const GSHGEOM_Point bl_c = pt(bottom_end.x, bl_end.y);
  const GSHGEOM_Point left_end = pt(x, from_near(y, rtl_y));
  const GSHGEOM_Point tl_c = pt(start.x, left_end.y);
  const GSHGEOM_Point all[] = {start, top_end, tr_end, tr_c, right_end, br_end, br_c, bottom_end, bl_end, bl_c, left_end, tl_c};
  for (size_t i = 0; i < sizeof all / sizeof all[0]; i++) {
    if (out_of_range(all[i].x) || out_of_range(all[i].y)) {
      return GSHGEOM_ERR_LIMIT;
    }
  }
  GSHGEOM_Result e = gshgeom_path_move_to(path, start);
  if (e == GSHGEOM_OK) {
    e = gshgeom_path_line_to(path, top_end);
  }
  if (e == GSHGEOM_OK && round[1]) {
    e = quarter(path, top_end, tr_end, tr_c, pt(rtr_x, 0.0), pt(0.0, rtr_y), -90.0);
  }
  if (e == GSHGEOM_OK) {
    e = gshgeom_path_line_to(path, right_end);
  }
  if (e == GSHGEOM_OK && round[2]) {
    e = quarter(path, right_end, br_end, br_c, pt(rbr_x, 0.0), pt(0.0, rbr_y), 0.0);
  }
  if (e == GSHGEOM_OK) {
    e = gshgeom_path_line_to(path, bottom_end);
  }
  if (e == GSHGEOM_OK && round[3]) {
    e = quarter(path, bottom_end, bl_end, bl_c, pt(rbl_x, 0.0), pt(0.0, rbl_y), 90.0);
  }
  if (e == GSHGEOM_OK) {
    e = gshgeom_path_line_to(path, left_end);
  }
  if (e == GSHGEOM_OK && round[0]) {
    e = quarter(path, left_end, start, tl_c, pt(rtl_x, 0.0), pt(0.0, rtl_y), 180.0);
  }
  return e == GSHGEOM_OK ? gshgeom_path_close(path) : e;
}

static GSHGEOM_Result add_ellipse(
    GSHGEOM_Path * path, double cx, double cy, double rx, double ry, int flip_x, int flip_y) {
  if (begin(path) != GSHGEOM_OK || (flip_x != 0 && flip_x != 1) || (flip_y != 0 && flip_y != 1)) {
    return GSHGEOM_ERR_INVALID;
  }
  if (!is_finite_num(cx) || !is_finite_num(cy) || !is_finite_num(rx) || !is_finite_num(ry) || !(rx > 0.0) || !(ry > 0.0)) {
    return GSHGEOM_ERR_INVALID;
  }
  const double sx = flip_x ? -rx : rx;
  const double sy = flip_y ? -ry : ry;
  const GSHGEOM_Point c = pt(cx, cy);
  const GSHGEOM_Point u = pt(sx, 0.0);
  const GSHGEOM_Point v = pt(0.0, sy);
  /* The ends `cx +- sx` and `cy +- sy` are checked, and a centre or a radius past the limit puts one of the pair past it. */
  const GSHGEOM_Point e[4] = {
      pt(sum2(cx, sx), cy), pt(cx, sum2(cy, sy)), pt(sum2(cx, -sx), cy), pt(cx, sum2(cy, -sy))};
  for (unsigned i = 0; i < 4; i++) {
    if (out_of_range(e[i].x) || out_of_range(e[i].y)) {
      return GSHGEOM_ERR_LIMIT;
    }
  }
  GSHGEOM_Result r = gshgeom_path_move_to(path, e[0]);
  for (unsigned k = 0; r == GSHGEOM_OK && k < 4; k++) {
    const GSHGEOM_Point from = e[k];
    const GSHGEOM_Point to = e[(k + 1) % 4];
    r = quarter(path, from, to, c, u, v, 90.0 * (double)k);
  }
  return r == GSHGEOM_OK ? gshgeom_path_close(path) : r;
}

/* Each public builder is all or nothing: a failure part way (the element limit, an allocation) takes back what it added. */

GSHGEOM_Result gshgeom_path_add_rect(
    GSHGEOM_Path * path, double x, double y, double w, double h, int flip_x, int flip_y) {
  PathMark mark;
  if (begin(path) != GSHGEOM_OK) {
    return GSHGEOM_ERR_INVALID;
  }
  path_mark(path, &mark);
  GSHGEOM_Result r = add_rect(path, x, y, w, h, flip_x, flip_y);
  if (r != GSHGEOM_OK) {
    path_rollback(path, &mark);
  }
  return r;
}

GSHGEOM_Result gshgeom_path_add_rounded_rect(
    GSHGEOM_Path * path, double x, double y, double w, double h, const GSHGEOM_Corner r[4]) {
  PathMark mark;
  if (begin(path) != GSHGEOM_OK) {
    return GSHGEOM_ERR_INVALID;
  }
  path_mark(path, &mark);
  GSHGEOM_Result e = add_rounded_rect(path, x, y, w, h, r);
  if (e != GSHGEOM_OK) {
    path_rollback(path, &mark);
  }
  return e;
}

GSHGEOM_Result gshgeom_path_add_ellipse(
    GSHGEOM_Path * path, double cx, double cy, double rx, double ry, int flip_x, int flip_y) {
  PathMark mark;
  if (begin(path) != GSHGEOM_OK) {
    return GSHGEOM_ERR_INVALID;
  }
  path_mark(path, &mark);
  GSHGEOM_Result r = add_ellipse(path, cx, cy, rx, ry, flip_x, flip_y);
  if (r != GSHGEOM_OK) {
    path_rollback(path, &mark);
  }
  return r;
}
