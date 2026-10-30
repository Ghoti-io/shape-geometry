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
 * The conic route's seam (internal): the route that makes the device path of a rational curve
 * (device.c, "The conic route"), with the pieces' parameter intervals handed out and a depth cap
 * a test can lower.
 *
 * The proof is that no input reaches the cap with a piece that is not flat, so nothing the public
 * functions do can show that the count and the bound are looked at.  The capped variant takes a
 * smaller cap and reports how many pieces it forced, so a test can use a cap a curve does reach.
 * These are not exported from the shared library; the tests link the archive.
 */

#ifndef GHOTI_IO_GSHGEOM_CONIC_INTERNAL_H
#define GHOTI_IO_GSHGEOM_CONIC_INTERNAL_H

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/shape-geometry/flatten.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * The route for a rational curve of @p n (3 or 4) control points, given in device units (26.6, each
 * rounding into `int32`) with weights @p w, to the tolerance @p tolerance (at least
 * ::GSHGEOM_FLATTEN_MIN_TOLERANCE).  Each flat piece is handed to @p emit as its end vertex, rounded
 * once, with the piece's place in the parameter interval as in ::GSHGEOM_Flatten_Emit (@p leaf is
 * exact while @p level is at most 32).  The first point is not emitted.
 *
 * @return ::GSHGEOM_OK; ::GSHGEOM_ERR_INVALID for a bad argument; ::GSHGEOM_ERR_LIMIT for a vertex that
 *   does not round into `int32`; whatever @p emit returns; ::GSHGEOM_ERR_INTERNAL if the depth cap ended a
 *   piece that was not flat (which the proof says is never so).
 */
GSHGEOM_Result gshgeom_conic_route(int n, const double * sx, const double * sy,
    const double * w, int32_t tolerance, GSHGEOM_Flatten_Emit emit, void * user);

/**
 * ::gshgeom_conic_route with the depth cap at @p levels (0 through ::GSHGEOM_CONIC_MAX_DEPTH).
 *
 * @param forced If not NULL, receives how many pieces the cap ended although they were not flat:
 *   zero at the real cap.  A call that forces pieces is not an error here.
 * @param deepest If not NULL, receives the greatest depth of any piece.
 */
GSHGEOM_Result gshgeom_conic_route_capped(int n, const double * sx, const double * sy,
    const double * w, int32_t tolerance, int levels, GSHGEOM_Flatten_Emit emit,
    void * user, int * forced, int * deepest);

#ifdef __cplusplus
}
#endif

#endif /* GHOTI_IO_GSHGEOM_CONIC_INTERNAL_H */
