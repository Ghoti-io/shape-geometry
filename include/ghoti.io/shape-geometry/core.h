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
 * Result codes, limits and the version this build reports.
 */

#ifndef GHOTI_IO_GSHGEOM_CORE_H
#define GHOTI_IO_GSHGEOM_CORE_H

#include <ghoti.io/shape-geometry/macros.h>

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Result of an operation.
 *
 * Zero is success.  ::GSHGEOM_RESULT_COUNT closes the enum so that a test can
 * check the string table is complete.  This is the suite's fixed vocabulary
 * with nothing added; some of it is never returned by this library (it reads
 * no bytes), and is here so that the one enum is the same on every library of
 * the family.
 */
typedef enum {
  GSHGEOM_OK = 0,          ///< The operation succeeded.
  GSHGEOM_ERR_IO,          ///< A write to a caller's stream failed.
  GSHGEOM_ERR_FORMAT,      ///< Not returned by this library.
  GSHGEOM_ERR_UNSUPPORTED, ///< A stroke this library does not define yet (stroke.h: a 3D centerline, a width profile, a DXF pattern that cannot be lowered).
  GSHGEOM_ERR_LIMIT,       ///< A stated cap was exceeded; nothing was truncated.
  GSHGEOM_ERR_CORRUPT,     ///< Not returned by this library.
  GSHGEOM_ERR_OOM,         ///< The allocator returned NULL.
  GSHGEOM_ERR_INVALID,     ///< A caller-supplied argument is wrong.
  GSHGEOM_ERR_INTERNAL,    ///< The library's own invariant failed.
  GSHGEOM_RESULT_COUNT
} GSHGEOM_Result;

/**
 * @brief Static description of a result code.
 *
 * The string is never NULL, never allocated, and never contains caller data.
 *
 * @param result The result code, including values outside the enum.
 * @return A static string.
 */
GSHGEOM_API const char * gshgeom_result_string(GSHGEOM_Result result);

/** Default ::GSHGEOM_Limits.max_elements: 4 Mi path elements. */
#define GSHGEOM_LIMITS_DEFAULT_MAX_ELEMENTS ((size_t)1 << 22)

/** Default ::GSHGEOM_Limits.max_points: 4 Mi points in a flattened path. */
#define GSHGEOM_LIMITS_DEFAULT_MAX_POINTS ((size_t)1 << 22)

/**
 * @brief Caps for the quantities a caller can make unbounded.
 *
 * A field of 0 means the default above.  NULL where a `const GSHGEOM_Limits *`
 * is taken means every default.  Exceeding a cap is ::GSHGEOM_ERR_LIMIT, never
 * a silent truncation, and the object is left as it was.
 */
typedef struct GSHGEOM_Limits {
  size_t max_elements;  ///< Most elements (verbs) one model or device path may hold.
  size_t max_points;    ///< Most points one flattened path may hold, and most rows (quadrature panels) one measured path may hold.
  uint8_t _reserved[8]; ///< Reserved; zero.
} GSHGEOM_Limits;

/**
 * @brief Fill @p limits with the built-in defaults.
 *
 * @param limits The structure to fill; NULL is ignored.
 */
GSHGEOM_API void gshgeom_limits_default(GSHGEOM_Limits * limits);

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
