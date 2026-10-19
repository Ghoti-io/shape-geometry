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
 * @file fill.h
 *
 * The fill rules of the shape family (AD-4): which points a path encloses.
 * This header holds an enumeration and a string function and nothing else, so
 * that a library with no use for the rest of this one (`font`, whose own
 * headers carry no floating-point type, AD-5) can include it alone.
 *
 * The values are pinned: 0 is nonzero (the default, and what every font format
 * means), 1 is even-odd.  A further rule is appended after them; the count
 * moves, the values do not.  `font`'s `GFNT_FillRule` is this type under its
 * own name (story 9).
 */

#ifndef GHOTI_IO_GSHGEOM_FILL_H
#define GHOTI_IO_GSHGEOM_FILL_H

#include <ghoti.io/shape-geometry/macros.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Which points a path encloses.
 */
typedef enum GSHGEOM_Fill_Rule {
  /** A point is inside where the winding number is not zero.  The default. */
  GSHGEOM_FILL_NONZERO = 0,
  /** A point is inside where the winding number is odd. */
  GSHGEOM_FILL_EVEN_ODD = 1,
  GSHGEOM_FILL_RULE_COUNT ///< Closes the enum.
} GSHGEOM_Fill_Rule;

/**
 * @brief Name a fill rule.
 *
 * @param rule The rule, including a value outside the enum.
 * @return A static string, never NULL: "non-zero", "even-odd", or "unknown
 *   fill rule".
 */
GSHGEOM_API const char * gshgeom_fill_rule_string(GSHGEOM_Fill_Rule rule);

#ifdef __cplusplus
}
#endif

#endif /* GHOTI_IO_GSHGEOM_FILL_H */
