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
 * The measured path's seam for the tests (internal): ::gshgeom_measured_create, also reporting the units of the work budget
 * (the limits' `max_work`) that the refinement spent.  Not installed and not part of the interface; the tests link the archive.
 */

#ifndef GHOTI_IO_GSHGEOM_ALONG_INTERNAL_H
#define GHOTI_IO_GSHGEOM_ALONG_INTERNAL_H

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/shape-geometry/along.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * ::gshgeom_measured_create, also writing the work units spent to @p work (NULL to leave it): on success the units the
 * refinement charged, on ::GSHGEOM_ERR_LIMIT the units charged before the one that was refused (or before the row cap was met).
 */
GSHGEOM_Result gshgeom_measured_create_counted(const GSHGEOM_Allocator * alloc, const GSHGEOM_Limits * limits,
    const GSHGEOM_Path * path, double tolerance, GSHGEOM_Measured ** out, uint64_t * work);

#ifdef __cplusplus
}
#endif

#endif /* GHOTI_IO_GSHGEOM_ALONG_INTERNAL_H */
