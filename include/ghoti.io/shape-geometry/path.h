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
 * @file path.h
 *
 * The model path: an opaque sequence of subpaths made of lines and quadratic,
 * cubic and rational curves, in 2D or 3D (AD-4).
 *
 * A path is built by verbs, as in SVG and canvas: ::gshgeom_path_move_to
 * starts a subpath, a drawing verb adds a segment from the current point to a
 * new one, ::gshgeom_path_close ends the subpath with a straight line back to
 * its start.  After a close the current point is the subpath's start, so a
 * drawing verb then begins a new subpath there.  A drawing verb with no
 * current point, and a close with no subpath, are ::GSHGEOM_ERR_INVALID.
 *
 * Every point is checked as ::gshgeom_curve_check checks it (non-finite:
 * ::GSHGEOM_ERR_INVALID; out of range: ::GSHGEOM_ERR_LIMIT), and a failed call
 * leaves the path as it was.  An element count past the limit is
 * ::GSHGEOM_ERR_LIMIT.  A path is used from one thread at a time.
 *
 * ### The hash of a path (version 1)
 *
 * In the spelling of hash.h: `u32` 1 (this format), `u32` dimension, `u64`
 * element count; then for each element a `u8` verb and, in order, the
 * coordinates of its points as `f64` (`x`, `y`, and `z` in 3D): for
 * ::GSHGEOM_VERB_MOVE the point, for a drawing verb its control points after
 * the start and then the end, and for a rational verb the weights (3 or 4,
 * the start's first) as `f64` after the points; a close adds nothing.  The
 * limits and the allocator are not part of it.
 *
 * ### The fixed-point builder (story 14a, AD-5)
 *
 * ::gshgeom_path_move_to_fixed and the three drawing verbs after it take each
 * coordinate as an `int32_t` `v` and a number of fraction bits `b` from 0 to
 * ::GSHGEOM_FIXED_MAX_BITS, and add the point `v / 2^b`, which is exact in a
 * `double` (an `int32_t` has 31 bits of magnitude, a `double` 53, and the
 * division is by a power of two), with `z = 0`.  Nothing rounds: the builder
 * is how a library that names no floating type (`font`, whose design units are
 * 26.6 integers, `b` = 6) hands its coordinates to a model path.  A fraction
 * count past the bound is ::GSHGEOM_ERR_INVALID; otherwise each is exactly the
 * verb it names, with that verb's checks and its atomicity.
 */

#ifndef GHOTI_IO_GSHGEOM_PATH_H
#define GHOTI_IO_GSHGEOM_PATH_H

#include <ghoti.io/shape-geometry/allocator.h>
#include <ghoti.io/shape-geometry/core.h>
#include <ghoti.io/shape-geometry/curve.h>
#include <ghoti.io/shape-geometry/macros.h>

#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief A path under construction or in use.  Opaque. */
typedef struct GSHGEOM_Path GSHGEOM_Path;

/** @brief The verb of a path element.  The values are pinned (they are hashed). */
typedef enum GSHGEOM_Verb {
  GSHGEOM_VERB_MOVE = 0,              ///< Start a subpath.
  GSHGEOM_VERB_LINE = 1,              ///< A line to a point.
  GSHGEOM_VERB_QUADRATIC = 2,         ///< A quadratic curve.
  GSHGEOM_VERB_CUBIC = 3,             ///< A cubic curve.
  GSHGEOM_VERB_RATIONAL_QUADRATIC = 4, ///< A rational quadratic curve.
  GSHGEOM_VERB_RATIONAL_CUBIC = 5,    ///< A rational cubic curve.
  GSHGEOM_VERB_CLOSE = 6,             ///< A line back to the subpath's start.
  GSHGEOM_VERB_COUNT                  ///< Closes the enum.
} GSHGEOM_Verb;

/**
 * @brief One element as read back.
 *
 * For a drawing verb, `curve` is the segment, its first control point the
 * start.  For ::GSHGEOM_VERB_CLOSE it is the line from the current point to the
 * subpath's start (possibly of length zero).  For ::GSHGEOM_VERB_MOVE it is a
 * line of length zero at the new point.
 */
typedef struct GSHGEOM_Element {
  GSHGEOM_Verb verb;   ///< The verb.
  GSHGEOM_Curve curve; ///< The segment, as above.
} GSHGEOM_Element;

/**
 * @brief Make an empty path.
 *
 * @param alloc The allocator, or NULL for the default; copied.
 * @param limits The limits, or NULL for the defaults; copied.
 * @param dim 2 or 3.
 * @param out The path; written only on success.  Free it with
 *   ::gshgeom_path_destroy.
 * @return ::GSHGEOM_OK, ::GSHGEOM_ERR_INVALID (dim, NULL @p out) or
 *   ::GSHGEOM_ERR_OOM.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_path_create(const GSHGEOM_Allocator * alloc,
    const GSHGEOM_Limits * limits, int dim, GSHGEOM_Path ** out);

/** @brief Free a path.  NULL is ignored. */
GSHGEOM_API void gshgeom_path_destroy(GSHGEOM_Path * path);

/** @brief The dimension the path was made with, 2 or 3; 0 for NULL. */
GSHGEOM_API int gshgeom_path_dimension(const GSHGEOM_Path * path);

/** @brief How many elements the path holds; 0 for NULL. */
GSHGEOM_API size_t gshgeom_path_element_count(const GSHGEOM_Path * path);

/**
 * @brief Start a subpath at @p p.
 *
 * @return ::GSHGEOM_OK, ::GSHGEOM_ERR_INVALID, ::GSHGEOM_ERR_LIMIT or
 *   ::GSHGEOM_ERR_OOM.  Likewise for every verb below.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_path_move_to(
    GSHGEOM_Path * path, GSHGEOM_Point p);

/** @brief A line from the current point to @p p. */
GSHGEOM_API GSHGEOM_Result gshgeom_path_line_to(
    GSHGEOM_Path * path, GSHGEOM_Point p);

/** @brief A quadratic curve from the current point, control point @p c, to @p p. */
GSHGEOM_API GSHGEOM_Result gshgeom_path_quadratic_to(
    GSHGEOM_Path * path, GSHGEOM_Point c, GSHGEOM_Point p);

/** @brief A cubic curve from the current point, control points @p c1, @p c2, to @p p. */
GSHGEOM_API GSHGEOM_Result gshgeom_path_cubic_to(GSHGEOM_Path * path,
    GSHGEOM_Point c1, GSHGEOM_Point c2, GSHGEOM_Point p);

/**
 * @brief A rational quadratic curve from the current point, control point
 * @p c, to @p p.
 *
 * @param w Three weights, for the start, @p c and @p p (positive, in range).
 *   The start's weight belongs to this curve alone.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_path_rational_quadratic_to(
    GSHGEOM_Path * path, GSHGEOM_Point c, GSHGEOM_Point p, const double w[3]);

/**
 * @brief A rational cubic curve from the current point, control points
 * @p c1 and @p c2, to @p p.
 *
 * @param w Four weights, in the order of the points.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_path_rational_cubic_to(GSHGEOM_Path * path,
    GSHGEOM_Point c1, GSHGEOM_Point c2, GSHGEOM_Point p, const double w[4]);

/** @brief End the subpath with a line back to its start. */
GSHGEOM_API GSHGEOM_Result gshgeom_path_close(GSHGEOM_Path * path);

/** The most fraction bits a fixed-point coordinate may have (see the file comment). */
#define GSHGEOM_FIXED_MAX_BITS 16u

/**
 * @brief ::gshgeom_path_move_to at `(x / 2^bits, y / 2^bits)`, exactly.
 *
 * @param bits Fraction bits, 0 to ::GSHGEOM_FIXED_MAX_BITS.
 * @return As ::gshgeom_path_move_to; ::GSHGEOM_ERR_INVALID for @p bits past the bound.  Likewise for the three below.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_path_move_to_fixed(
    GSHGEOM_Path * path, int32_t x, int32_t y, unsigned bits);

/** @brief ::gshgeom_path_line_to at `(x, y) / 2^bits`, exactly. */
GSHGEOM_API GSHGEOM_Result gshgeom_path_line_to_fixed(
    GSHGEOM_Path * path, int32_t x, int32_t y, unsigned bits);

/** @brief ::gshgeom_path_quadratic_to with the control point `(cx, cy)` and the end `(x, y)`, each `/ 2^bits`. */
GSHGEOM_API GSHGEOM_Result gshgeom_path_quadratic_to_fixed(GSHGEOM_Path * path,
    int32_t cx, int32_t cy, int32_t x, int32_t y, unsigned bits);

/** @brief ::gshgeom_path_cubic_to with the control points `(c1x, c1y)`, `(c2x, c2y)` and the end `(x, y)`, each `/ 2^bits`. */
GSHGEOM_API GSHGEOM_Result gshgeom_path_cubic_to_fixed(GSHGEOM_Path * path,
    int32_t c1x, int32_t c1y, int32_t c2x, int32_t c2y, int32_t x, int32_t y,
    unsigned bits);

/**
 * @brief Read element @p index back.
 *
 * @param path The path.
 * @param index From 0 to the element count minus one.
 * @param out The element; written only on success.
 * @return ::GSHGEOM_OK or ::GSHGEOM_ERR_INVALID (NULL, index out of range).
 */
GSHGEOM_API GSHGEOM_Result gshgeom_path_element(
    const GSHGEOM_Path * path, size_t index, GSHGEOM_Element * out);

/**
 * @brief The box of the path: the union of the boxes of its segments.
 *
 * A move that starts no segment and a closing line add nothing (a closing line
 * lies between points already in the box).  A path with no segment gives
 * `empty`.
 *
 * @param path The path.
 * @param out The box; written only on success.
 * @return ::GSHGEOM_OK or ::GSHGEOM_ERR_INVALID.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_path_bounds(
    const GSHGEOM_Path * path, GSHGEOM_Bounds * out);

/**
 * @brief The hash of the path, in the spelling in this file's comment.
 *
 * @param path The path, or NULL.
 * @return The hash; for NULL, the offset basis.
 */
GSHGEOM_API uint64_t gshgeom_path_hash(const GSHGEOM_Path * path);

/**
 * @brief Write the path as text and end with its hash.
 *
 * The first line is `gshgeom.path 1 dim=D elements=N hash=H`; then one line
 * per element: `M`, `L`, `Q`, `C`, `RQ`, `RC` or `Z`, then the numbers.  A
 * rational element lists its weights first.  Every number is the 16
 * hexadecimal digits of the bit pattern of a `double`, so the text does not
 * depend on the locale or on the C library's number formatting.
 *
 * @param path The path.
 * @param out The stream.
 * @return ::GSHGEOM_OK, ::GSHGEOM_ERR_INVALID (NULL) or ::GSHGEOM_ERR_IO.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_path_dump(
    const GSHGEOM_Path * path, FILE * out);

#ifdef __cplusplus
}
#endif

#endif /* GHOTI_IO_GSHGEOM_PATH_H */
