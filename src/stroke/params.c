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
 * A stroke's parameters: the defaults (SG1), validity and bounds (SG2), the
 * canonical meaning (SG3), and the not-invertible predicate (SG15).
 */

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/shape-geometry/orient.h>
#include <ghoti.io/shape-geometry/stroke.h>

#include "stroke_internal.h"

#include <string.h>

/** A finite double: neither NaN nor an infinity (no libm: x - x is 0 exactly when x is finite). */
static int is_finite_value(double x) {
  return x - x == 0.0;
}

/** SG2: zero, or a magnitude in [VALUE_MIN, VALUE_MAX] (@p x finite). */
static int within(double x) {
  double m = x < 0.0 ? -x : x;
  return m == 0.0 || (m >= GSHGEOM_STROKE_VALUE_MIN && m <= GSHGEOM_STROKE_VALUE_MAX);
}

/** Whether a join reads the miter limit (SG3). */
static int join_reads_miter(uint32_t join) {
  return join == (uint32_t)GSHGEOM_STROKE_JOIN_MITER || join == (uint32_t)GSHGEOM_STROKE_JOIN_MITER_CLIP ||
      join == (uint32_t)GSHGEOM_STROKE_JOIN_ARCS;
}

void gshgeom_stroke_init(GSHGEOM_Stroke * stroke) {
  if (!stroke) {
    return;
  }
  memset(stroke, 0, sizeof *stroke);
  stroke->width = 1.0;
  stroke->width_space = GSHGEOM_STROKE_WIDTH_LOCAL;
  stroke->join = GSHGEOM_STROKE_JOIN_MITER;
  stroke->miter_limit = GSHGEOM_STROKE_MITER_LIMIT_DEFAULT;
  stroke->cap = GSHGEOM_STROKE_CAP_BUTT;
  stroke->degenerate = GSHGEOM_STROKE_DEGENERATE_DRAW;
  stroke->dash_scope = GSHGEOM_STROKE_DASH_PER_SUBPATH;
  stroke->alignment = GSHGEOM_STROKE_ALIGN_CENTER;
  stroke->alignment_rule = GSHGEOM_FILL_NONZERO;
}

/** SG3 (a): the canonical dash count: 0 for a solid stroke, the count doubled when it is odd. */
static size_t canonical_dash_count(const GSHGEOM_Stroke * s) {
  int any = 0;
  for (size_t i = 0; i < s->dash_count; i++) {
    if (s->dashes[i] < 0.0) {
      return 0; /* a negative length: solid */
    }
    any = any || s->dashes[i] != 0.0;
  }
  if (!any) {
    return 0; /* none, or every length zero: solid */
  }
  if (s->dash_count % 2u) {
    return s->dash_count > SIZE_MAX / 2u ? SIZE_MAX : s->dash_count * 2u;
  }
  return s->dash_count;
}

GSHGEOM_Result gshgeom_stroke_check(const GSHGEOM_Stroke * s) {
  if (!s) {
    return GSHGEOM_ERR_INVALID;
  }
  if (s->dash_count && !s->dashes) {
    return GSHGEOM_ERR_INVALID;
  }
  if (!is_finite_value(s->width) || !is_finite_value(s->miter_limit) || !is_finite_value(s->dash_offset)) {
    return GSHGEOM_ERR_INVALID;
  }
  for (size_t i = 0; i < s->dash_count; i++) {
    if (!is_finite_value(s->dashes[i])) {
      return GSHGEOM_ERR_INVALID;
    }
  }
  if (s->width_space >= (uint32_t)GSHGEOM_STROKE_WIDTH_SPACE_COUNT || s->join >= (uint32_t)GSHGEOM_STROKE_JOIN_COUNT ||
      s->cap >= (uint32_t)GSHGEOM_STROKE_CAP_COUNT || s->degenerate >= (uint32_t)GSHGEOM_STROKE_DEGENERATE_COUNT ||
      s->dash_scope >= (uint32_t)GSHGEOM_STROKE_DASH_SCOPE_COUNT || s->alignment >= (uint32_t)GSHGEOM_STROKE_ALIGN_COUNT ||
      (unsigned)s->alignment_rule >= (unsigned)GSHGEOM_FILL_RULE_COUNT) {
    return GSHGEOM_ERR_INVALID;
  }
  if (s->width < 0.0) {
    return GSHGEOM_ERR_INVALID;
  }
  for (size_t i = 0; i < sizeof s->_reserved; i++) {
    if (s->_reserved[i]) {
      return GSHGEOM_ERR_INVALID;
    }
  }
  if (s->width_profile) {
    return GSHGEOM_ERR_UNSUPPORTED; /* reserved until a profile is defined */
  }
  if (!within(s->width) || !within(s->dash_offset) ||
      !(s->miter_limit >= 1.0 && s->miter_limit <= GSHGEOM_STROKE_MITER_LIMIT_MAX)) {
    return GSHGEOM_ERR_LIMIT;
  }
  for (size_t i = 0; i < s->dash_count; i++) {
    if (!within(s->dashes[i])) {
      return GSHGEOM_ERR_LIMIT;
    }
  }
  if (canonical_dash_count(s) > GSHGEOM_STROKE_DASHES_MAX) {
    return GSHGEOM_ERR_LIMIT;
  }
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_stroke_canonical(const GSHGEOM_Stroke * in, GSHGEOM_Stroke * out) {
  GSHGEOM_Result r = gshgeom_stroke_check(in);
  if (r != GSHGEOM_OK) {
    return r;
  }
  size_t dash_count = canonical_dash_count(in);
  memset(out, 0, sizeof *out);
  out->width = in->width == 0.0 ? 0.0 : in->width;
  out->width_space = in->width_space;
  out->join = in->join;
  out->miter_limit = join_reads_miter(in->join) ? in->miter_limit : GSHGEOM_STROKE_MITER_LIMIT_DEFAULT;
  out->cap = in->cap;
  /* a butt cap draws nothing for a degenerate subpath whether it is drawn or omitted */
  out->degenerate = in->cap == (uint32_t)GSHGEOM_STROKE_CAP_BUTT ? (uint32_t)GSHGEOM_STROKE_DEGENERATE_DRAW : in->degenerate;
  /* the canonical count (an odd array repeated once, SG3); the entries are the caller's, entry i being
   * in->dashes[i % in->dash_count] */
  out->dash_count = dash_count;
  out->dashes = dash_count ? in->dashes : NULL;
  out->dash_offset = dash_count && in->dash_offset != 0.0 ? in->dash_offset : 0.0;
  out->dash_scope = dash_count ? in->dash_scope : (uint32_t)GSHGEOM_STROKE_DASH_PER_SUBPATH;
  out->alignment = in->alignment;
  out->alignment_rule = in->alignment == (uint32_t)GSHGEOM_STROKE_ALIGN_CENTER ? GSHGEOM_FILL_NONZERO
                                                                                : in->alignment_rule;
  out->width_profile = NULL;
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_stroke_map_invertible(const GSHGEOM_Map * map, int * invertible) {
  if (!map || !invertible || gshgeom_map_check(map) != GSHGEOM_OK || map->dim != 2) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Point o = {0.0, 0.0, 0.0};
  GSHGEOM_Point c0 = {map->m[0][0], map->m[1][0], 0.0};
  GSHGEOM_Point c1 = {map->m[0][1], map->m[1][1], 0.0};
  int sign = 0;
  if (gshgeom_orient2d(o, c0, c1, &sign) != GSHGEOM_OK) {
    return GSHGEOM_ERR_INTERNAL;
  }
  *invertible = sign != 0;
  return GSHGEOM_OK;
}

int gshgeom_stroke_revision(void) {
  return GSHGEOM_STROKE_REVISION;
}
