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
 * tolerance is part of the definition (AD-12): ::GSHGEOM_FLATTEN_TOLERANCE is
 * the one a path is flattened to unless the caller states another, which a
 * sink with its own tolerance option (`font`'s) does, through the curve
 * functions or ::gshgeom_device_path_flatten_with_tolerance.
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
 * The proof holds for every tolerance from ::GSHGEOM_FLATTEN_MIN_TOLERANCE, 2
 * units (1/32 pixel), up to INT32_MAX, with the flatness threshold derived from
 * the tolerance and the same depth cap.  A smaller tolerance is refused with
 * ::GSHGEOM_ERR_INVALID before anything is emitted: a tolerance asked for is a
 * promise, so one the proof cannot keep is not rounded up to another.  The
 * result is the same polyline on every target.
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

/** The default flattening tolerance: 4 units of 1/64 pixel, i.e. 1/16 pixel. */
#define GSHGEOM_FLATTEN_TOLERANCE 4

/**
 * The smallest tolerance accepted, in 26.6 units (two units, 1/32 pixel): the
 * proof of the guarantee is closed from here up (flatten.c, "Tolerances").
 */
#define GSHGEOM_FLATTEN_MIN_TOLERANCE 2

/**
 * How deep a curve may be subdivided before its chord is used as it is.  The
 * proof shows the cap is never reached; it is 20, not 16, because at 16 a
 * curve spanning the `int32` range would still be over the threshold.
 */
#define GSHGEOM_FLATTEN_MAX_DEPTH 20

/** The revision of the flattener's definition; see the table above. */
#define GSHGEOM_FLATTEN_REVISION 2

/**
 * Receives one vertex of a flattened curve, in order.
 *
 * @param user  The pointer given to the flattener.
 * @param to    The vertex, 26.6.  The polyline starts at the curve's first
 *              point, which is not emitted; the last vertex is its last point.
 * @param leaf  Which piece of the curve this vertex ends: the piece covers the
 *              parameter interval `[leaf / 2^level, (leaf + 1) / 2^level]`, so
 *              a test, or a caller that interpolates an attribute along the
 *              curve, knows where on the curve the segment replaces it.
 * @param level How many times the curve was halved to reach that piece.
 * @return ::GSHGEOM_OK to go on; anything else stops the flattening and is
 *         returned by it (a sink that is out of room says so this way).
 */
typedef GSHGEOM_Result (*GSHGEOM_Flatten_Emit)(
    void * user, GSHGEOM_Device_Point to, uint32_t leaf, int level);

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
 * @brief Flatten one quadratic Bezier with the integer flattener.
 *
 * The polyline's vertices are handed to @p emit in order, each rounded to 26.6
 * once; the guarantee above holds for any @p p0, @p p1 and @p p2 and any
 * accepted @p tolerance, and the depth cap is never what ends a piece.  The
 * caller flips, maps or places the points before the call and not after: the
 * half-up midpoint is not symmetric under negation, so flattening and then
 * flipping gives other vertices than flipping and then flattening.
 *
 * @param tolerance How far the polyline may stray, in 26.6 units; at least
 *                  ::GSHGEOM_FLATTEN_MIN_TOLERANCE.
 * @param p0   First point (any `int32` coordinates).
 * @param p1   The control point.
 * @param p2   Last point; the last vertex emitted.
 * @param emit Called once per vertex.
 * @param user Handed to @p emit.
 * @return ::GSHGEOM_OK; ::GSHGEOM_ERR_INVALID for a tolerance below the
 *         minimum or a NULL @p emit (nothing is emitted); whatever @p emit
 *         returned that was not ::GSHGEOM_OK, at which the flattening stopped.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_flatten_quadratic(int32_t tolerance,
    GSHGEOM_Device_Point p0, GSHGEOM_Device_Point p1, GSHGEOM_Device_Point p2,
    GSHGEOM_Flatten_Emit emit, void * user);

/** @brief Flatten one cubic Bezier; see ::gshgeom_flatten_quadratic. */
GSHGEOM_API GSHGEOM_Result gshgeom_flatten_cubic(int32_t tolerance,
    GSHGEOM_Device_Point p0, GSHGEOM_Device_Point p1, GSHGEOM_Device_Point p2,
    GSHGEOM_Device_Point p3, GSHGEOM_Flatten_Emit emit, void * user);

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
 * @brief ::gshgeom_device_path_flatten at a tolerance of the caller's.
 *
 * Exactly ::gshgeom_device_path_flatten when @p tolerance is
 * ::GSHGEOM_FLATTEN_TOLERANCE.
 *
 * @param tolerance As for ::gshgeom_flatten_quadratic.
 * @return As for ::gshgeom_device_path_flatten, and ::GSHGEOM_ERR_INVALID for a
 *         tolerance below ::GSHGEOM_FLATTEN_MIN_TOLERANCE, before anything is
 *         allocated.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_device_path_flatten_with_tolerance(
    const GSHGEOM_Device_Path * path, int32_t tolerance,
    const GSHGEOM_Allocator * alloc, const GSHGEOM_Limits * limits,
    GSHGEOM_Flat_Path ** out);

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

/**
 * @brief The same winding everywhere, without the edges that cancel.
 *
 * Every subpath is read as closed, so a path is a multiset of directed edges.  An edge `a -> b` and an edge `b -> a`
 * cross every ray from a point together or not at all, with opposite signs, so a pair of them changes the winding
 * number nowhere; an edge of zero length changes nothing either.  The result is @p path without such pairs (as many of
 * a repeated edge as have no twin remain, with their multiplicity) and without zero-length edges, regrouped into closed
 * subpaths: the edges left have as many leaving each point as entering it, so they split into closed walks.  The winding
 * number at every point not on an edge is the same as @p path's, and so is the nonzero and the even-odd fill, which are
 * functions of it.  The stroke expansion builds its outline of positively wound pieces whose neighbours share their sides
 * in opposite directions, so about two fifths of an outline's edges cancel, and they are the long ones.
 *
 * The result depends on the multiset of edges and not on the order of the subpaths or where each starts.  Its points are
 * no more than @p path's.  `limits->max_points` bounds it as it bounds the flattener's.
 *
 * @param path   The flattened path; its subpaths are read as closed whatever their flag.
 * @param alloc  The allocator for the result and the scratch, or NULL for the default.
 * @param limits The limits, or NULL for the defaults.
 * @param out    The result, to be destroyed with ::gshgeom_flat_path_destroy; written only on success.  It has no subpaths
 *               when everything cancels.
 * @return ::GSHGEOM_OK, ::GSHGEOM_ERR_INVALID (NULL @p path or @p out), ::GSHGEOM_ERR_LIMIT, ::GSHGEOM_ERR_OOM, or
 *         ::GSHGEOM_ERR_INTERNAL for a path whose edges are not a union of closed subpaths (none that this library builds).
 */
GSHGEOM_API GSHGEOM_Result gshgeom_flat_path_cancel_seams(const GSHGEOM_Flat_Path * path,
    const GSHGEOM_Allocator * alloc, const GSHGEOM_Limits * limits, GSHGEOM_Flat_Path ** out);

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
