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
 * Maps: constructors, the one composer, the exact classifier, application, the
 * route predicate, the hash and the dump (see map.h).
 */

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/shape-geometry/hash.h>
#include <ghoti.io/shape-geometry/map.h>

#include <ghoti.io/math/basic.h>
#include <ghoti.io/math/elementary.h>

#include "../exact/exact.h"

#include <stdint.h>

/** The hash and dump format of a map. */
#define MAP_FORMAT_VERSION 1u

/** Radians per degree, per gradian and per turn, correctly rounded. */
#define RAD_PER_DEGREE 0x1.1df46a2529d39p-6
#define RAD_PER_GRADIAN 0x1.015bf9217271ap-6
#define RAD_PER_TURN 0x1.921fb54442d18p+2

static int is_finite(double v) {
  return (v - v) == 0.0;
}

static double absolute(double v) {
  return v < 0.0 ? -v : v;
}

/** Make a negative zero a positive one (and change nothing else). */
static double nz(double v) {
  return v == 0.0 ? 0.0 : v;
}

/** NaN or infinity anywhere is INVALID (first), a finite number past the bound LIMIT. */
static GSHGEOM_Result check_numbers(const double * v, int count) {
  for (int i = 0; i < count; i++) {
    if (!is_finite(v[i])) {
      return GSHGEOM_ERR_INVALID;
    }
  }
  for (int i = 0; i < count; i++) {
    if (v[i] > GSHGEOM_MAP_LIMIT || v[i] < -GSHGEOM_MAP_LIMIT) {
      return GSHGEOM_ERR_LIMIT;
    }
  }
  return GSHGEOM_OK;
}

static int dim_ok(int dim) {
  return dim == 2 || dim == 3;
}

/* ---------------------------------------------------------------------- */
/* The classifier                                                          */
/* ---------------------------------------------------------------------- */

/** The exact sign of the determinant of the linear part. */
static int det_sign(int dim, double m[3][3]) {
  GSHGEOM_ExactTerm t[6];
  size_t n = 0;
  if (dim == 2) {
    gshgeom_exact_put(t, &n, 0, 2, m[0][0], m[1][1], 0.0);
    gshgeom_exact_put(t, &n, 1, 2, m[0][1], m[1][0], 0.0);
  } else {
    gshgeom_exact_put_det3(t, &n, 0, m[0], m[1], m[2]);
  }
  return gshgeom_exact_sign(t, n);
}

/**
 * The tag of a matrix and translation, from the values.
 *
 * A row with exactly one nonzero entry in every row (and a nonzero
 * determinant, so one in every column) is axis aligned: diagonal, or a
 * permutation of the axes (a quarter turn, with a mirror if the determinant is
 * negative).  A diagonal one is classified by the signs and magnitudes of its
 * entries.  Otherwise a 2D matrix `[[a, -b], [b, a]]` is a rotation and
 * `[[a, b], [b, -a]]` a mirrored one (both exact equalities of entries), and
 * anything else a skew; a 3D one of no axis-aligned structure is general.
 */
static uint32_t classify(int dim, double m[3][3], const double t[3]) {
  uint32_t kind = 0;
  for (int i = 0; i < dim; i++) {
    if (t[i] != 0.0) {
      kind |= GSHGEOM_MAP_TRANSLATES;
    }
  }
  int sign = det_sign(dim, m);
  if (sign == 0) {
    return kind | GSHGEOM_MAP_SINGULAR;
  }
  int col[3] = {0, 0, 0};
  int mono = 1;
  for (int i = 0; i < dim; i++) {
    int nonzero = 0;
    for (int j = 0; j < dim; j++) {
      if (m[i][j] != 0.0) {
        nonzero++;
        col[i] = j;
      }
    }
    if (nonzero != 1) {
      mono = 0;
    }
  }
  if (mono) {
    int identity_perm = 1;
    int negatives = 0;
    int y_negative = 0;
    double mag[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < dim; i++) {
      double v = m[i][col[i]];
      if (col[i] != i) {
        identity_perm = 0;
      }
      if (v < 0.0) {
        negatives++;
        if (i == 1) {
          y_negative = 1;
        }
      }
      mag[i] = absolute(v);
    }
    if (identity_perm) {
      if (negatives == 1) {
        kind |= y_negative ? GSHGEOM_MAP_FLIPS_Y : GSHGEOM_MAP_MIRRORS;
      } else if (negatives == 2) {
        kind |= GSHGEOM_MAP_QUARTER_TURNS; /* a half turn */
      } else if (negatives == 3) {
        kind |= GSHGEOM_MAP_MIRRORS; /* a point reflection */
      }
    } else {
      kind |= GSHGEOM_MAP_QUARTER_TURNS;
      if (sign < 0) {
        kind |= GSHGEOM_MAP_MIRRORS;
      }
    }
    int uniform = 1;
    for (int i = 1; i < dim; i++) {
      if (mag[i] != mag[0]) {
        uniform = 0;
      }
    }
    if (!uniform) {
      kind |= GSHGEOM_MAP_SCALES_NONUNIFORM;
    } else if (mag[0] != 1.0) {
      kind |= GSHGEOM_MAP_SCALES;
    }
    return kind;
  }
  if (dim == 2) {
    double a = m[0][0];
    double b = m[0][1];
    double c = m[1][0];
    double d = m[1][1];
    if (a == d && b == -c) {
      kind |= GSHGEOM_MAP_ROTATES;
    } else if (a == -d && b == c) {
      kind |= GSHGEOM_MAP_ROTATES | GSHGEOM_MAP_MIRRORS;
    } else {
      kind |= GSHGEOM_MAP_SKEWS;
    }
    return kind;
  }
  return kind | GSHGEOM_MAP_GENERAL;
}

/* ---------------------------------------------------------------------- */
/* Construction                                                            */
/* ---------------------------------------------------------------------- */

/** Fill a map from the dim-by-dim part of `m` and `t`, with the tag classified. */
static void fill(GSHGEOM_Map * out, int dim, double m[3][3], const double t[3]) {
  GSHGEOM_Map r;
  r.dim = dim;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      r.m[i][j] = (i < dim && j < dim) ? nz(m[i][j]) : (i == j ? 1.0 : 0.0);
    }
    r.t[i] = i < dim ? nz(t[i]) : 0.0;
  }
  r.kind = classify(dim, r.m, r.t);
  *out = r;
}

GSHGEOM_Result gshgeom_map_identity(int dim, GSHGEOM_Map * out) {
  if (!out || !dim_ok(dim)) {
    return GSHGEOM_ERR_INVALID;
  }
  double m[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  const double t[3] = {0, 0, 0};
  fill(out, dim, m, t);
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_map_translate(
    int dim, double tx, double ty, double tz, GSHGEOM_Map * out) {
  if (!out || !dim_ok(dim)) {
    return GSHGEOM_ERR_INVALID;
  }
  const double t[3] = {tx, ty, tz};
  GSHGEOM_Result r = check_numbers(t, dim);
  if (r != GSHGEOM_OK) {
    return r;
  }
  double m[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  fill(out, dim, m, t);
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_map_scale(
    int dim, double sx, double sy, double sz, GSHGEOM_Map * out) {
  if (!out || !dim_ok(dim)) {
    return GSHGEOM_ERR_INVALID;
  }
  const double s[3] = {sx, sy, sz};
  GSHGEOM_Result r = check_numbers(s, dim);
  if (r != GSHGEOM_OK) {
    return r;
  }
  double m[3][3] = {{sx, 0, 0}, {0, sy, 0}, {0, 0, sz}};
  const double t[3] = {0, 0, 0};
  fill(out, dim, m, t);
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_map_flip_y(int dim, GSHGEOM_Map * out) {
  return gshgeom_map_scale(dim, 1.0, -1.0, 1.0, out);
}

/**
 * The cosine and sine of an angle, with the quarter turns exact.
 *
 * In degrees, gradians and turns: the nearest whole number of quarter turns
 * `n` is found, the remainder `angle - n * quarter` is formed (exact: `n *
 * quarter` is an integer below 2^47 or a power-of-two multiple of one, and the
 * difference of two doubles within a factor of two is exact), the remainder
 * goes through `gmath_cos` and `gmath_sin` after one multiplication by the
 * radians per unit, and the quarter turns are applied by swapping and negating
 * `(c, s)`, which is exact.
 */
static GSHGEOM_Result angle_cs(GSHGEOM_Angle_Unit unit, double angle, double * c, double * s) {
  if ((unsigned)unit >= (unsigned)GSHGEOM_ANGLE_UNIT_COUNT || !is_finite(angle)) {
    return GSHGEOM_ERR_INVALID;
  }
  if (unit == GSHGEOM_ANGLE_RADIANS) {
    if (angle == 0.0) {
      *c = 1.0;
      *s = 0.0;
    } else {
      *c = gmath_cos(angle);
      *s = gmath_sin(angle);
    }
    return GSHGEOM_OK;
  }
  double quarter;
  double per_unit;
  double u;
  if (unit == GSHGEOM_ANGLE_DEGREES) {
    quarter = 90.0;
    per_unit = RAD_PER_DEGREE;
    u = angle / quarter;
  } else if (unit == GSHGEOM_ANGLE_GRADIANS) {
    quarter = 100.0;
    per_unit = RAD_PER_GRADIAN;
    u = angle / quarter;
  } else {
    quarter = 0.25;
    per_unit = RAD_PER_TURN;
    u = angle * 4.0;
  }
  double nr = gmath_round(u);
  if (nr > GSHGEOM_MAP_MAX_QUARTERS || nr < -GSHGEOM_MAP_MAX_QUARTERS) {
    return GSHGEOM_ERR_LIMIT;
  }
  double rest = angle - nr * quarter;
  /* The number of quarter turns modulo 4, without a conversion of a double to
   * a 64-bit integer (on i686 that is an x87 instruction, which the
   * determinism gate refuses): nr is a whole number below 2^40, nr / 4 is
   * exact, and so are the floor and the difference. */
  double kq = nr - 4.0 * gmath_floor(nr * 0.25);
  int k = (int)kq;
  double cc;
  double ss;
  if (rest == 0.0) {
    cc = 1.0;
    ss = 0.0;
  } else {
    double theta = rest * per_unit;
    cc = gmath_cos(theta);
    ss = gmath_sin(theta);
  }
  switch (k) {
  case 1:
    *c = 0.0 - ss;
    *s = cc;
    break;
  case 2:
    *c = 0.0 - cc;
    *s = 0.0 - ss;
    break;
  case 3:
    *c = ss;
    *s = 0.0 - cc;
    break;
  default:
    *c = cc;
    *s = ss;
    break;
  }
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_map_rotate(
    GSHGEOM_Angle_Unit unit, double angle, GSHGEOM_Map * out) {
  if (!out) {
    return GSHGEOM_ERR_INVALID;
  }
  double c;
  double s;
  GSHGEOM_Result r = angle_cs(unit, angle, &c, &s);
  if (r != GSHGEOM_OK) {
    return r;
  }
  double m[3][3] = {{c, 0.0 - s, 0}, {s, c, 0}, {0, 0, 1}};
  const double t[3] = {0, 0, 0};
  fill(out, 2, m, t);
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_map_rotate_axis(
    GSHGEOM_Axis axis, GSHGEOM_Angle_Unit unit, double angle, GSHGEOM_Map * out) {
  if (!out || (unsigned)axis >= (unsigned)GSHGEOM_AXIS_COUNT) {
    return GSHGEOM_ERR_INVALID;
  }
  double c;
  double s;
  GSHGEOM_Result r = angle_cs(unit, angle, &c, &s);
  if (r != GSHGEOM_OK) {
    return r;
  }
  double m[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  /* A positive angle turns the next axis toward the one after. */
  if (axis == GSHGEOM_AXIS_Z) {
    m[0][0] = c;
    m[0][1] = 0.0 - s;
    m[1][0] = s;
    m[1][1] = c;
  } else if (axis == GSHGEOM_AXIS_X) {
    m[1][1] = c;
    m[1][2] = 0.0 - s;
    m[2][1] = s;
    m[2][2] = c;
  } else {
    m[2][2] = c;
    m[2][0] = 0.0 - s;
    m[0][2] = s;
    m[0][0] = c;
  }
  const double t[3] = {0, 0, 0};
  fill(out, 3, m, t);
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_map_from_matrix(
    int dim, const double * m, const double * t, GSHGEOM_Map * out) {
  if (!out || !m || !t || !dim_ok(dim)) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Result r = check_numbers(m, dim * dim);
  if (r == GSHGEOM_OK) {
    r = check_numbers(t, dim);
  }
  if (r != GSHGEOM_OK) {
    /* A NaN in the translation beats a large number in the matrix. */
    for (int i = 0; i < dim; i++) {
      if (!is_finite(t[i])) {
        return GSHGEOM_ERR_INVALID;
      }
    }
    return r;
  }
  double mm[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  double tt[3] = {0, 0, 0};
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      mm[i][j] = m[i * dim + j];
    }
    tt[i] = t[i];
  }
  fill(out, dim, mm, tt);
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_map_check(const GSHGEOM_Map * map) {
  if (!map || !dim_ok(map->dim) || (map->kind & ~(uint32_t)GSHGEOM_MAP_KIND_ALL)) {
    return GSHGEOM_ERR_INVALID;
  }
  int dim = map->dim;
  double v[12];
  int n = 0;
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      v[n++] = map->m[i][j];
    }
  }
  for (int i = 0; i < dim; i++) {
    v[n++] = map->t[i];
  }
  return check_numbers(v, n);
}

/* ---------------------------------------------------------------------- */
/* Composition: the one place maps are multiplied                          */
/* ---------------------------------------------------------------------- */

/** `outer . inner`, in the stated order of operations. */
static GSHGEOM_Result compose_one(
    const GSHGEOM_Map * outer, const GSHGEOM_Map * inner, GSHGEOM_Map * out) {
  GSHGEOM_Result r = gshgeom_map_check(outer);
  if (r == GSHGEOM_OK) {
    r = gshgeom_map_check(inner);
  }
  if (r != GSHGEOM_OK) {
    return r;
  }
  if (outer->dim != inner->dim) {
    return GSHGEOM_ERR_INVALID;
  }
  int dim = outer->dim;
  GSHGEOM_Map res;
  res.dim = dim;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      res.m[i][j] = i == j ? 1.0 : 0.0;
    }
    res.t[i] = 0.0;
  }
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      double acc = outer->m[i][0] * inner->m[0][j];
      acc = acc + outer->m[i][1] * inner->m[1][j];
      if (dim == 3) {
        acc = acc + outer->m[i][2] * inner->m[2][j];
      }
      res.m[i][j] = nz(acc);
    }
    double tv = outer->m[i][0] * inner->t[0];
    tv = tv + outer->m[i][1] * inner->t[1];
    if (dim == 3) {
      tv = tv + outer->m[i][2] * inner->t[2];
    }
    res.t[i] = nz(tv + outer->t[i]);
  }
  /* The union of the factors' tags, and singular if the product has
   * collapsed (two small scales multiplied to zero). */
  res.kind = outer->kind | inner->kind;
  r = gshgeom_map_check(&res);
  if (r != GSHGEOM_OK) {
    return r;
  }
  if (det_sign(dim, res.m) == 0) {
    res.kind |= GSHGEOM_MAP_SINGULAR;
  }
  *out = res;
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_map_compose(
    const GSHGEOM_Map * outer, const GSHGEOM_Map * inner, GSHGEOM_Map * out) {
  if (!outer || !inner || !out) {
    return GSHGEOM_ERR_INVALID;
  }
  return compose_one(outer, inner, out);
}

GSHGEOM_Result gshgeom_map_compose_chain(const GSHGEOM_Map * chain,
    size_t count, const GSHGEOM_Map * target, GSHGEOM_Map * out) {
  if (!out || (count && !chain) || (!count && !target)) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Map acc;
  if (count) {
    acc = chain[0];
    GSHGEOM_Result r = gshgeom_map_check(&acc);
    if (r != GSHGEOM_OK) {
      return r;
    }
    for (size_t i = 1; i < count; i++) {
      r = compose_one(&acc, &chain[i], &acc);
      if (r != GSHGEOM_OK) {
        return r;
      }
    }
    if (target) {
      r = compose_one(target, &acc, &acc);
      if (r != GSHGEOM_OK) {
        return r;
      }
    }
  } else {
    acc = *target;
    GSHGEOM_Result r = gshgeom_map_check(&acc);
    if (r != GSHGEOM_OK) {
      return r;
    }
  }
  *out = acc;
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_map_design_to_model(
    double scale, double tx, double ty, GSHGEOM_Map * out) {
  if (!out || !is_finite(scale) || !(scale > 0.0)) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Map chain[3];
  GSHGEOM_Result r = gshgeom_map_translate(2, tx, ty, 0.0, &chain[0]);
  if (r == GSHGEOM_OK) {
    r = gshgeom_map_scale(2, scale, scale, 1.0, &chain[1]);
  }
  if (r == GSHGEOM_OK) {
    r = gshgeom_map_flip_y(2, &chain[2]);
  }
  if (r != GSHGEOM_OK) {
    return r;
  }
  return gshgeom_map_compose_chain(chain, 3, NULL, out);
}

/* ---------------------------------------------------------------------- */
/* Application                                                             */
/* ---------------------------------------------------------------------- */

static GSHGEOM_Result apply(
    const GSHGEOM_Map * map, GSHGEOM_Point p, int translate, GSHGEOM_Point * out) {
  GSHGEOM_Result r = gshgeom_map_check(map);
  if (r != GSHGEOM_OK) {
    return r;
  }
  if (!out || !is_finite(p.x) || !is_finite(p.y) || (map->dim == 3 && !is_finite(p.z))) {
    return GSHGEOM_ERR_INVALID;
  }
  double in[3] = {p.x, p.y, map->dim == 3 ? p.z : 0.0};
  double res[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < map->dim; i++) {
    double v = map->m[i][0] * in[0];
    v = v + map->m[i][1] * in[1];
    if (map->dim == 3) {
      v = v + map->m[i][2] * in[2];
    }
    if (translate) {
      v = v + map->t[i];
    }
    if (!is_finite(v)) {
      return GSHGEOM_ERR_LIMIT;
    }
    res[i] = v;
  }
  out->x = res[0];
  out->y = res[1];
  out->z = res[2];
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_map_apply_point(
    const GSHGEOM_Map * map, GSHGEOM_Point p, GSHGEOM_Point * out) {
  return apply(map, p, 1, out);
}

GSHGEOM_Result gshgeom_map_apply_vector(
    const GSHGEOM_Map * map, GSHGEOM_Point v, GSHGEOM_Point * out) {
  return apply(map, v, 0, out);
}

GSHGEOM_Result gshgeom_map_apply_curve(
    const GSHGEOM_Map * map, const GSHGEOM_Curve * curve, GSHGEOM_Curve * out) {
  if (!map || !curve || !out) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Result r = gshgeom_map_check(map);
  if (r == GSHGEOM_OK) {
    r = gshgeom_curve_check(curve);
  }
  if (r != GSHGEOM_OK) {
    return r;
  }
  if (curve->dim != map->dim) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Curve res = *curve;
  int n = gshgeom_curve_point_count(curve->kind);
  for (int i = 0; i < n; i++) {
    r = apply(map, curve->p[i], 1, &res.p[i]);
    if (r != GSHGEOM_OK) {
      return r;
    }
  }
  r = gshgeom_curve_check(&res);
  if (r != GSHGEOM_OK) {
    return r;
  }
  *out = res;
  return GSHGEOM_OK;
}

/* ---------------------------------------------------------------------- */
/* The route predicate                                                     */
/* ---------------------------------------------------------------------- */

GSHGEOM_Result gshgeom_map_grid_aligned(
    const GSHGEOM_Map * map, int * aligned, GSHGEOM_Grid * grid) {
  if (!aligned) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Result r = gshgeom_map_check(map);
  if (r != GSHGEOM_OK) {
    return r;
  }
  *aligned = 0;
  if (map->dim != 2) {
    return GSHGEOM_OK;
  }
  const uint32_t allowed =
      GSHGEOM_MAP_TRANSLATES | GSHGEOM_MAP_SCALES | GSHGEOM_MAP_FLIPS_Y;
  if (map->kind & ~allowed) {
    return GSHGEOM_OK;
  }
  double a = map->m[0][0];
  double d = map->m[1][1];
  if (map->m[0][1] != 0.0 || map->m[1][0] != 0.0 || !(a > 0.0)) {
    return GSHGEOM_OK;
  }
  int flip;
  if (d == a) {
    flip = 0;
  } else if (d == -a) {
    flip = 1;
  } else {
    return GSHGEOM_OK;
  }
  *aligned = 1;
  if (grid) {
    grid->scale = a;
    grid->tx = map->t[0];
    grid->ty = map->t[1];
    grid->flip_y = flip;
  }
  return GSHGEOM_OK;
}

/* ---------------------------------------------------------------------- */
/* Hash and dump                                                           */
/* ---------------------------------------------------------------------- */

uint64_t gshgeom_map_hash(const GSHGEOM_Map * map) {
  GSHGEOM_Hash h;
  gshgeom_hash_init(&h);
  if (gshgeom_map_check(map) != GSHGEOM_OK) {
    return gshgeom_hash_value(&h);
  }
  gshgeom_hash_u32(&h, MAP_FORMAT_VERSION);
  gshgeom_hash_u32(&h, (uint32_t)map->dim);
  gshgeom_hash_u32(&h, map->kind);
  for (int i = 0; i < map->dim; i++) {
    for (int j = 0; j < map->dim; j++) {
      gshgeom_hash_f64(&h, map->m[i][j]);
    }
  }
  for (int i = 0; i < map->dim; i++) {
    gshgeom_hash_f64(&h, map->t[i]);
  }
  return gshgeom_hash_value(&h);
}

GSHGEOM_Result gshgeom_map_dump(const GSHGEOM_Map * map, FILE * out) {
  if (!out || gshgeom_map_check(map) != GSHGEOM_OK) {
    return GSHGEOM_ERR_INVALID;
  }
  if (fprintf(out, "gshgeom.map %u dim=%d kind=%08llx hash=%016llx\n",
          MAP_FORMAT_VERSION, map->dim, (unsigned long long)map->kind,
          (unsigned long long)gshgeom_map_hash(map)) < 0) {
    return GSHGEOM_ERR_IO;
  }
  for (int i = 0; i < map->dim; i++) {
    if (fputc('m', out) == EOF) {
      return GSHGEOM_ERR_IO;
    }
    for (int j = 0; j < map->dim; j++) {
      if (fprintf(out, " %016llx", (unsigned long long)gshgeom_f64_bits(map->m[i][j])) < 0) {
        return GSHGEOM_ERR_IO;
      }
    }
    if (fputc('\n', out) == EOF) {
      return GSHGEOM_ERR_IO;
    }
  }
  if (fputc('t', out) == EOF) {
    return GSHGEOM_ERR_IO;
  }
  for (int i = 0; i < map->dim; i++) {
    if (fprintf(out, " %016llx", (unsigned long long)gshgeom_f64_bits(map->t[i])) < 0) {
      return GSHGEOM_ERR_IO;
    }
  }
  if (fputc('\n', out) == EOF) {
    return GSHGEOM_ERR_IO;
  }
  return fflush(out) == 0 ? GSHGEOM_OK : GSHGEOM_ERR_IO;
}
