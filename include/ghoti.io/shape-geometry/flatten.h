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
 * @file flatten.h
 *
 * The integer flattener and the flattened path (CAP-3, CAP-13; AD-11, AD-12).
 *
 * Flattening is a call into this library after the map, with error in device
 * pixels.  The flattener here is the definition's: it was `font`'s private one,
 * and a raster sink uses it so that no two sinks flatten differently.  Its
 * tolerance is part of the definition (AD-12) and is not a parameter.
 *
 * ## The guarantee (flattener revision 2)
 *
 * For every curve whose control points are any `int32` values, every point of
 * the curve is within the tolerance less 0.12 units of the emitted polyline: at
 * most 3.88 units (a little over 1/17 pixel) for the tolerance of 4 units
 * (::GSHGEOM_FLATTEN_TOLERANCE, 1/16 pixel), and the depth cap is never what
 * ends a piece.  The proof is stated above the code (src/flatten/flatten.c) and
 * in documentation/design.md, and tests/unit/test_flatten.cpp checks the bound
 * against dense exact evaluation on every size class up to the whole `int32`
 * range.
 *
 * The proof holds for any tolerance of 2 units (1/32 pixel) and up, with the
 * flatness threshold derived from it.  The tolerance is not a parameter of
 * this library, so there is no value below 2 to refuse here; the compile-time
 * check in flatten.c refuses a ::GSHGEOM_FLATTEN_TOLERANCE under 2.  `font`'s
 * rasteriser, which does take a tolerance, returns `GFNT_ERR_INVALID` for 1 or
 * a negative value; this library's contract is the same one, held at the only
 * value it has.
 *
 * ## The algorithm
 *
 * Points are widened by eight fractional bits and a curve is split at the
 * midpoint of its widened control points (each midpoint rounded half up, as
 * `floor((u + v + 1) / 2)`) until a quadratic's `|2 p1 - p0 - p2|`, or both of
 * a cubic's `|3 p1 - 2 p0 - p3|` and `|3 p2 - p0 - 2 p3|` (each summed over
 * both axes), are at most `4 * (t - 1)` units, i.e. 12 units at `t = 4`; to at
 * most ::GSHGEOM_FLATTEN_MAX_DEPTH levels.  The end point of every piece is
 * rounded once to 26.6, to nearest with ties upward (`floor((v + 128) / 256)`
 * of the widened value), and emitted; duplicates are included, in order.
 * Arithmetic is 64-bit and uses no shift of a negative value.
 *
 * ## The revision
 *
 * The flattener is versioned, and its polyline is what a rebaseline of every
 * sink's goldens follows.  ::GSHGEOM_FLATTEN_REVISION (and
 * ::gshgeom_flatten_revision) is the number; it is recorded where a rebaseline
 * lives (a golden's header, a vector file) and is in no hash: the flattened
 * path's spelling (::gshgeom_flat_path_hash) is unchanged, and two revisions
 * that gave the same polyline hash the same.  There is no way to ask for an
 * older revision.
 *
 * | Revision | Meaning |
 * |---|---|
 * | 1 | `font`'s flattener as ported in story 7: every midpoint floored to 26.6, flat at `4 * t` units, depth 16.  Its error against the exact curve reached 13 units on a large curve (13.1 over the `int32` range), 2 to 4 on a tiny one. |
 * | 2 | The above: widened by eight bits, a vertex rounded once, flat at `4 * (t - 1)`, depth 20, proved within `t - 0.12` units. |
 *
 */

#ifndef GHOTI_IO_GSHGEOM_FLATTEN_H
#define GHOTI_IO_GSHGEOM_FLATTEN_H

#include <ghoti.io/shape-geometry/allocator.h>
#include <ghoti.io/shape-geometry/core.h>
#include <ghoti.io/shape-geometry/device.h>
#include <ghoti.io/shape-geometry/macros.h>
#include <ghoti.io/shape-geometry/map.h>
#include <ghoti.io/shape-geometry/path.h>

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/** The flattening tolerance: 4 units of 1/64 pixel, i.e. 1/16 pixel. */
#define GSHGEOM_FLATTEN_TOLERANCE 4

/**
 * How deep a curve may be subdivided before its chord is used as it is.  The
 * proof shows the cap is never reached; it is 20, not 16, because at 16 a
 * curve spanning the `int32` range would still be over the threshold.
 */
#define GSHGEOM_FLATTEN_MAX_DEPTH 20

/** The revision of the flattener's definition; see the table above. */
#define GSHGEOM_FLATTEN_REVISION 2

/** The opaque flattened path: polylines in 26.6, open or closed. */
typedef struct GSHGEOM_Flat_Path GSHGEOM_Flat_Path;

/**
 * One polyline of a flattened path: `count` points starting at `first` in
 * ::gshgeom_flat_path_points.  The first is the subpath's start.  A closed
 * subpath's closing line, from the last point back to the first, is implicit.
 */
typedef struct GSHGEOM_Flat_Subpath {
  size_t first;  ///< Index of the first point.
  size_t count;  ///< Number of points (at least 1).
  int closed;    ///< 1 if the subpath was closed.
} GSHGEOM_Flat_Subpath;

/**
 * @brief The flattener revision this library implements
 *        (::GSHGEOM_FLATTEN_REVISION), for a caller that records or checks it
 *        across a dynamic-link boundary.
 */
GSHGEOM_API int gshgeom_flatten_revision(void);

/**
 * @brief Flatten a device path with the integer flattener.
 *
 * A move starts a subpath (a lone move is a subpath of one point); a line adds
 * its end point; a quadratic or cubic adds the end points of the lines it
 * flattens to; a close marks the subpath closed.
 *
 * @param path   The device path.
 * @param alloc  The allocator; NULL for the default.
 * @param limits Caps (`max_points`); NULL for the defaults.
 * @param out    Receives the flattened path, only on success.
 * @return ::GSHGEOM_OK; ::GSHGEOM_ERR_INVALID for NULL;
 *         ::GSHGEOM_ERR_LIMIT if more than `max_points` points would result;
 *         ::GSHGEOM_ERR_OOM.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_device_path_flatten(
    const GSHGEOM_Device_Path * path, const GSHGEOM_Allocator * alloc,
    const GSHGEOM_Limits * limits, GSHGEOM_Flat_Path ** out);

/**
 * @brief The whole of CAP-3 in one call: map a model path, round it to 26.6
 *        and flatten it, after the map.
 *
 * `limits->max_elements` bounds the intermediate device path and `max_points`
 * the result.  Errors are those of ::gshgeom_device_path_from_model and
 * ::gshgeom_device_path_flatten.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_path_flatten(const GSHGEOM_Path * path,
    const GSHGEOM_Map * map, const GSHGEOM_Allocator * alloc,
    const GSHGEOM_Limits * limits, GSHGEOM_Flat_Path ** out);

/** @brief Destroy a flattened path; NULL is ignored. */
GSHGEOM_API void gshgeom_flat_path_destroy(GSHGEOM_Flat_Path * path);

/** @brief The number of points (0 for NULL). */
GSHGEOM_API size_t gshgeom_flat_path_point_count(const GSHGEOM_Flat_Path * path);

/** @brief The number of subpaths (0 for NULL). */
GSHGEOM_API size_t gshgeom_flat_path_subpath_count(
    const GSHGEOM_Flat_Path * path);

/**
 * @brief The points, all subpaths in order; valid until the path is destroyed.
 *
 * @return NULL for a NULL path or one with no points.
 */
GSHGEOM_API const GSHGEOM_Device_Point * gshgeom_flat_path_points(
    const GSHGEOM_Flat_Path * path);

/**
 * @brief Read subpath @p index back.
 *
 * @return ::GSHGEOM_OK, ::GSHGEOM_ERR_INVALID for NULL or an index past the
 *         end.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_flat_path_subpath(
    const GSHGEOM_Flat_Path * path, size_t index, GSHGEOM_Flat_Subpath * out);

/**
 * @brief The FNV-1a 64 hash of a flattened path, in the spelling of hash.h
 *        (AD-21).
 *
 * `u32` format version (1), `u32` kind (`0x54414c46`, "FLAT" little-endian),
 * `u64` subpath count, and for each subpath `u8` closed, `u64` point count and
 * the points, each `i32` x then `i32` y.  NULL hashes as the empty stream.
 */
GSHGEOM_API uint64_t gshgeom_flat_path_hash(const GSHGEOM_Flat_Path * path);

/**
 * @brief Write the dump: a header line `gshgeom.flat-path 1 subpaths=S
 *        points=P hash=<16 hex digits>`, then for each subpath a line `S
 *        closed=<0|1> points=<n>` and one line per point, `P x y`, each number
 *        the two's complement of the `int32_t` in 8 hexadecimal digits.
 *
 * @return ::GSHGEOM_OK, ::GSHGEOM_ERR_INVALID for NULL, ::GSHGEOM_ERR_IO when
 *         the stream fails (including at the flush).
 */
GSHGEOM_API GSHGEOM_Result gshgeom_flat_path_dump(
    const GSHGEOM_Flat_Path * path, FILE * out);

#ifdef __cplusplus
}
#endif

#endif /* GHOTI_IO_GSHGEOM_FLATTEN_H */
