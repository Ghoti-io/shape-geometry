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
 * The flattener's depth-cap seam (internal).
 *
 * The proof (flatten.c, (g)) is that no input reaches the depth cap, so
 * nothing the public functions do can show that a piece the cap forces is
 * counted, emitted and tiles the curve.  These variants take a smaller cap and
 * report how many pieces it forced, so a test can use a cap a large curve does
 * reach.  They are not exported from the shared library; the tests link the
 * archive.
 */

#ifndef GHOTI_IO_GSHGEOM_FLATTEN_INTERNAL_H
#define GHOTI_IO_GSHGEOM_FLATTEN_INTERNAL_H

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/shape-geometry/flatten.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * ::gshgeom_flatten_quadratic with the depth cap at @p levels (clamped to 1
 * through ::GSHGEOM_FLATTEN_MAX_DEPTH).
 *
 * @param forced If not NULL, receives how many pieces the cap ended although
 *   they were not flat: zero for any input at the real cap, and zero when the
 *   call is refused.
 */
GSHGEOM_Result gshgeom_flatten_quadratic_capped(int32_t tolerance, int levels,
    GSHGEOM_Device_Point p0, GSHGEOM_Device_Point p1, GSHGEOM_Device_Point p2,
    GSHGEOM_Flatten_Emit emit, void * user, int * forced);

/** ::gshgeom_flatten_cubic with a smaller depth cap; see the quadratic. */
GSHGEOM_Result gshgeom_flatten_cubic_capped(int32_t tolerance, int levels,
    GSHGEOM_Device_Point p0, GSHGEOM_Device_Point p1, GSHGEOM_Device_Point p2,
    GSHGEOM_Device_Point p3, GSHGEOM_Flatten_Emit emit, void * user,
    int * forced);

/**
 * An empty flattened path, built point by point by the stroke expansion
 * (src/stroke/).  `limits->max_points` bounds it as it bounds the flattener's.
 *
 * @return ::GSHGEOM_OK, ::GSHGEOM_ERR_OOM.
 */
GSHGEOM_Result gshgeom_flat_path_new(const GSHGEOM_Allocator * alloc,
    const GSHGEOM_Limits * limits, GSHGEOM_Flat_Path ** out);

/** Start a subpath at @p p.  ::GSHGEOM_ERR_LIMIT past `max_points`; ::GSHGEOM_ERR_OOM. */
GSHGEOM_Result gshgeom_flat_path_begin(GSHGEOM_Flat_Path * path, GSHGEOM_Device_Point p);

/** Add @p p to the last subpath.  ::GSHGEOM_ERR_LIMIT past `max_points`; ::GSHGEOM_ERR_OOM. */
GSHGEOM_Result gshgeom_flat_path_add(GSHGEOM_Flat_Path * path, GSHGEOM_Device_Point p);

/** Mark the last subpath closed. */
void gshgeom_flat_path_close_last(GSHGEOM_Flat_Path * path);

#ifdef __cplusplus
}
#endif

#endif
