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
 * The clip's seams for the tests (internal): the exact side test, the crossing, and the clip with the number of crossings
 * it computed.  Not installed and not part of the interface; the tests link the archive.
 */

#ifndef GHOTI_IO_GSHGEOM_CLIP_INTERNAL_H
#define GHOTI_IO_GSHGEOM_CLIP_INTERNAL_H

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/shape-geometry/clip.h>

#ifdef __cplusplus
extern "C" {
#endif

/** What a clip did, for the size rule (CK5). */
typedef struct GSHGEOM_Clip_Stats {
  size_t crossings;       /**< Crossing points computed, over all four cuts and every contour. */
  size_t contours_in;     /**< Contours cut (partly inside and partly outside some side). */
  size_t contours_copied; /**< Contours copied bit for bit (wholly inside). */
  size_t cuts;            /**< Passes over a contour's points by one side of the quadrilateral (the work of the cut). */
} GSHGEOM_Clip_Stats;

/** ::gshgeom_flat_path_clip_convex, also filling @p stats (NULL to leave it). */
GSHGEOM_Result gshgeom_clip_convex_counted(const GSHGEOM_Flat_Path * path, const GSHGEOM_Device_Point quad[4],
    const GSHGEOM_Allocator * alloc, const GSHGEOM_Limits * limits, GSHGEOM_Flat_Path ** out,
    GSHGEOM_Clip_Stats * stats);

/**
 * The sign (-1, 0, +1) of `(l1.x - l0.x) (p.y - l0.y) - (l1.y - l0.y) (p.x - l0.x)`, exactly: the side of `p` of the line
 * from `l0` to `l1`, by the convention of ::gshgeom_orient2d.
 */
int gshgeom_clip_side(GSHGEOM_Device_Point l0, GSHGEOM_Device_Point l1, GSHGEOM_Device_Point p);

/**
 * The crossing of the edge `a`-`b` with the line through `l0` and `l1` (CK3): 1 and `*out` set when `a` and `b` are
 * strictly on opposite sides of the line, else 0 and `*out` untouched.  `l0` and `l1` must differ.
 */
int gshgeom_clip_crossing(GSHGEOM_Device_Point a, GSHGEOM_Device_Point b, GSHGEOM_Device_Point l0,
    GSHGEOM_Device_Point l1, GSHGEOM_Device_Point * out);

#ifdef __cplusplus
}
#endif

#endif
