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
 * @file core.h
 *
 * The version this build reports.
 */

#ifndef GHOTI_IO_GSHGEOM_CORE_H
#define GHOTI_IO_GSHGEOM_CORE_H

#include <ghoti.io/shape-geometry/macros.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief This build's version, as the string the Makefile generated.
 *
 * @return A static string, never NULL.  "0.0.0", or "0.0.0-dev" when BRANCH
 *   was overridden.
 */
GSHGEOM_API const char * gshgeom_version_string(void);

/**
 * @brief This build's version, packed as ::GSHGEOM_MAKE_VERSION packs it.
 *
 * @return `(major << 16) | (minor << 8) | patch`.
 */
GSHGEOM_API unsigned gshgeom_version_number(void);

#ifdef __cplusplus
}
#endif

#endif /* GHOTI_IO_GSHGEOM_CORE_H */
