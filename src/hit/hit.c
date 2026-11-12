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
 * Hit testing against a flattened path (hit.h, HT1 to HT6).  One pass over the edges of every subpath, the closing edge
 * included: each edge is first asked whether the point lies on it, and otherwise whether it crosses the ray to the left of
 * the point.  Every decision is the exact sign of the side test of the clip (a product of two 33-bit differences, held in
 * 64 bits, and their difference compared without forming it), so no coordinate overflows.
 */

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/shape-geometry/hit.h>

#include "../clip/clip_internal.h"

typedef GSHGEOM_Device_Point Pt;

/** Whether @p p is on the closed segment @p a to @p b, given that they differ (HT3): collinear, and inside the box of the two. */
static int on_segment(Pt a, Pt b, Pt p) {
  const int32_t x0 = a.x < b.x ? a.x : b.x;
  const int32_t x1 = a.x < b.x ? b.x : a.x;
  const int32_t y0 = a.y < b.y ? a.y : b.y;
  const int32_t y1 = a.y < b.y ? b.y : a.y;
  return p.x >= x0 && p.x <= x1 && p.y >= y0 && p.y <= y1 && gshgeom_clip_side(a, b, p) == 0;
}

/**
 * What the edge @p a to @p b adds to the winding number of @p p (HT4), given that @p p is not on it: +1 for an edge directed
 * down that has @p p to its right at the height of @p p, -1 for one directed up, 0 for any other.  The edge spans the height
 * of @p p when `min <= p.y < max`, so of two edges that meet at a vertex on the ray exactly one counts when the contour
 * crosses the ray there and neither or both when it only touches it.
 */
static int crossing(Pt a, Pt b, Pt p) {
  if (a.y <= p.y && p.y < b.y) {          /* down */
    return gshgeom_clip_side(a, b, p) < 0 ? 1 : 0;
  }
  if (b.y <= p.y && p.y < a.y) {          /* up: the side test wants the upper end first */
    return gshgeom_clip_side(b, a, p) < 0 ? -1 : 0;
  }
  return 0;
}

GSHGEOM_Result gshgeom_flat_path_winding(const GSHGEOM_Flat_Path * path, GSHGEOM_Device_Point p, int64_t * winding, int * boundary) {
  if (!winding) {
    return GSHGEOM_ERR_INVALID;
  }
  *winding = 0;
  if (boundary) {
    *boundary = 0;
  }
  const size_t subpaths = gshgeom_flat_path_subpath_count(path);
  const Pt * pts = gshgeom_flat_path_points(path);
  int64_t w = 0;
  for (size_t s = 0; s < subpaths; s++) {
    GSHGEOM_Flat_Subpath sp;
    if (gshgeom_flat_path_subpath(path, s, &sp) != GSHGEOM_OK) {
      return GSHGEOM_ERR_INVALID;
    }
    for (size_t i = 0; i < sp.count; i++) {
      const Pt a = pts[sp.first + i];
      const Pt b = pts[sp.first + (i + 1 < sp.count ? i + 1 : 0)];
      if (a.x == b.x && a.y == b.y) {
        continue; /* no edge (HT2); a point on it is on a neighbour's edge if the subpath has one */
      }
      if (on_segment(a, b, p)) {
        if (boundary) {
          *boundary = 1;
        }
        return GSHGEOM_OK;
      }
      w += crossing(a, b, p);
    }
  }
  *winding = w;
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_flat_path_hit(const GSHGEOM_Flat_Path * path, GSHGEOM_Device_Point p, GSHGEOM_Fill_Rule rule,
    GSHGEOM_Hit * hit) {
  if (!hit || (rule != GSHGEOM_FILL_NONZERO && rule != GSHGEOM_FILL_EVEN_ODD)) {
    return GSHGEOM_ERR_INVALID;
  }
  int64_t w;
  int on;
  const GSHGEOM_Result r = gshgeom_flat_path_winding(path, p, &w, &on);
  if (r != GSHGEOM_OK) {
    return r;
  }
  if (on) {
    *hit = GSHGEOM_HIT_BOUNDARY;
  } else if (rule == GSHGEOM_FILL_NONZERO) {
    *hit = w != 0 ? GSHGEOM_HIT_INSIDE : GSHGEOM_HIT_OUTSIDE;
  } else {
    *hit = (w & 1) != 0 ? GSHGEOM_HIT_INSIDE : GSHGEOM_HIT_OUTSIDE;
  }
  return GSHGEOM_OK;
}

const char * gshgeom_hit_string(GSHGEOM_Hit hit) {
  switch (hit) {
  case GSHGEOM_HIT_OUTSIDE:
    return "outside";
  case GSHGEOM_HIT_INSIDE:
    return "inside";
  case GSHGEOM_HIT_BOUNDARY:
    return "boundary";
  case GSHGEOM_HIT_COUNT:
    break;
  }
  return "unknown hit";
}
