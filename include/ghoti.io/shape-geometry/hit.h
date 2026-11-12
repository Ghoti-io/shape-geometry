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
 * @file hit.h
 *
 * Hit testing: where a point of the device lattice lies against a flattened path (HT1 to HT6).
 *
 * - **HT1 The point and the path are both on the lattice.**  A point is a 26.6 device point (an `int32_t` count of 1/64 pixel
 *   in x and y, y down), and a flattened path is the polylines `flatten.h` makes.  A caller that has a position in pixels
 *   makes the point by the library's one rounding rule (`device.h`), and one that has a model path flattens it first, under
 *   the map and the tolerance it draws with, so that a hit is a statement about the pixels that are drawn.
 * - **HT2 Every subpath is closed.**  A fill closes an open subpath with the line from its last point to its first, so a hit
 *   test does too, whatever the subpath's `closed` flag says.  A subpath of one point or of points that are all the same has
 *   no edge and takes no part.
 * - **HT3 The boundary is decided exactly.**  A point is *on the boundary* when it lies on an edge of a subpath, an end point
 *   included, decided by the exact sign of the side test in integers (no rounding, no tolerance, and no overflow at any
 *   `int32_t` coordinate).  A boundary point is neither in nor out: its state is ::GSHGEOM_HIT_BOUNDARY under both rules, so
 *   a caller that wants the closed region takes it as inside and one that wants the open region as outside, and neither
 *   depends on which way a contour runs.
 * - **HT4 The winding number.**  For a point that is not on the boundary, the winding number `w` is the signed count of the
 *   edges that cross the horizontal ray to the left of the point, an edge directed down (y increasing) counting +1 and one
 *   directed up -1, the same sign `shape-cpu`'s scan converter gives a contour.  An edge with its end points at the heights
 *   `y0` and `y1` crosses the ray at the height `y` when `min(y0, y1) <= y < max(y0, y1)`, so a vertex on the ray, and an edge
 *   along it, are counted once and consistently: no pair of edges that share a vertex is counted twice or not at all.
 * - **HT5 The rules.**  The point is inside under the nonzero rule when `w != 0` and under the even-odd rule when `w` is odd
 *   (`fill.h`), outside otherwise, and these are the pixels a fill of the same path under the same rule covers wholly.
 * - **HT6 The count is exact and bounded.**  `w` is an `int64_t` that holds any count a path within `max_points` can make (every
 *   edge adds one at most), so it does not overflow; the work is one pass over the edges.
 */

#ifndef GHOTI_IO_GSHGEOM_HIT_H
#define GHOTI_IO_GSHGEOM_HIT_H

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/shape-geometry/core.h>
#include <ghoti.io/shape-geometry/device.h>
#include <ghoti.io/shape-geometry/fill.h>
#include <ghoti.io/shape-geometry/flatten.h>

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Where a point is against a path (the values are pinned). */
typedef enum GSHGEOM_Hit {
  GSHGEOM_HIT_OUTSIDE = 0,  ///< Not on the boundary and not enclosed under the rule.
  GSHGEOM_HIT_INSIDE = 1,   ///< Not on the boundary and enclosed under the rule.
  GSHGEOM_HIT_BOUNDARY = 2, ///< On an edge of a subpath (HT3), the end points included.
  GSHGEOM_HIT_COUNT         ///< Closes the enum.
} GSHGEOM_Hit;

/**
 * @brief The winding number of a point (HT2, HT4), and whether it is on the boundary (HT3).
 *
 * @param path     The flattened path; NULL is the empty path (winding 0, not on the boundary).
 * @param p        The point, in 26.6.
 * @param winding  Receives the winding number when the point is not on the boundary and 0 when it is.  Must not be NULL.
 * @param boundary Receives 1 when the point is on an edge and 0 when it is not.  May be NULL.
 * @return ::GSHGEOM_OK, or ::GSHGEOM_ERR_INVALID for a NULL @p winding (nothing is written).
 */
GSHGEOM_API GSHGEOM_Result gshgeom_flat_path_winding(const GSHGEOM_Flat_Path * path, GSHGEOM_Device_Point p, int64_t * winding,
    int * boundary);

/**
 * @brief Where a point is against a flattened path under a fill rule (HT2 to HT5).
 *
 * @param path The flattened path; NULL is the empty path.
 * @param p    The point, in 26.6.
 * @param rule The fill rule.
 * @param hit  Receives the answer.  Must not be NULL.
 * @return ::GSHGEOM_OK, or ::GSHGEOM_ERR_INVALID for a NULL @p hit or a @p rule outside ::GSHGEOM_Fill_Rule (nothing is
 *         written).
 */
GSHGEOM_API GSHGEOM_Result gshgeom_flat_path_hit(const GSHGEOM_Flat_Path * path, GSHGEOM_Device_Point p, GSHGEOM_Fill_Rule rule,
    GSHGEOM_Hit * hit);

/** @brief Name a hit; "unknown hit" for a value outside the enumeration, never NULL. */
GSHGEOM_API const char * gshgeom_hit_string(GSHGEOM_Hit hit);

#ifdef __cplusplus
}
#endif

#endif // GHOTI_IO_GSHGEOM_HIT_H
