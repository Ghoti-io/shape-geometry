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
 * @file device.h
 *
 * The 26.6 device path (CAP-3, CAP-13; AD-4, AD-6, AD-12): the one
 * fixed-point path any library exchanges, and the one rounding rule that
 * reaches it.
 *
 * Device space has its origin at the top-left corner of the top-left pixel,
 * y down, and a coordinate is an `int32_t` count of 1/64 pixel.  A device path
 * is never interned as a resource or stored in a scene: a sink makes one from
 * a model path and a map at render time, or a glyph path provider returns one.
 *
 * ## The rounding rule
 *
 * `floor(x * 64 + 0.5)`: round half up.  Half-way points go up, so
 * -0.5/64 gives 0 and -1.5/64 gives -1.  It is the mathematical function.
 * Computed as the `double` expression `floor(t + 0.5)` with `t = x * 64` it is
 * wrong for the one `t` just below one half, `0.5 - 2^-54`, where the sum
 * rounds to 1; here the floor of `t` and the fraction `t - floor(t)` (both
 * exact) are used instead, so no addition rounds.  No FMA is involved, and none
 * could change the result: the product by 64 is exact.
 *
 * ## The integer spelling (story 14a, AD-5)
 *
 * ::gshgeom_scale_26_6 is the same rule for a point a fixed-point source
 * scales itself: a 26.6 design coordinate `p` (1/64 unit) times a 16.16 scale
 * `S` (pixels per unit times 2^16), optionally negated (the y flip of a y-up
 * design space into y-down device space), is `floor((s p S + 2^15) / 2^16)`
 * 1/64 pixel with `s` the sign: half up, in y down, computed in 64-bit
 * integers.  It equals ::gshgeom_round_26_6 of the `double` route
 * `(s (p / 64)) (S / 65536)` wherever it answers: an answer fits an `int32_t`,
 * so `|p S| < 2^47 + 2^15`, below the `2^53` under which every step of that
 * route is exact.  A result that does not fit is ::GSHGEOM_ERR_LIMIT, as in the
 * `double` rule.
 * It is what a glyph path provider's fast tier rounds an unhinted point with,
 * so the fast tier with hinting off is the general route under the snapped
 * map bit for bit (shape-scene's glyph_provider.h).
 *
 * ## Out of range
 *
 * A value whose rounding does not fit an `int32_t` is ::GSHGEOM_ERR_LIMIT:
 * it is refused, never wrapped and never clipped.  Clipping in `double` before
 * asking for a device path is a sink's decision, and a sink that does not do it
 * reports the node refused.
 */

#ifndef GHOTI_IO_GSHGEOM_DEVICE_H
#define GHOTI_IO_GSHGEOM_DEVICE_H

#include <ghoti.io/shape-geometry/allocator.h>
#include <ghoti.io/shape-geometry/core.h>
#include <ghoti.io/shape-geometry/curve.h>
#include <ghoti.io/shape-geometry/macros.h>
#include <ghoti.io/shape-geometry/map.h>
#include <ghoti.io/shape-geometry/path.h>

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/** 1/64 pixel units in one pixel. */
#define GSHGEOM_DEVICE_ONE 64

/** A device point: 26.6 fixed point, y down. */
typedef struct GSHGEOM_Device_Point {
  int32_t x; ///< Horizontal, in 1/64 pixel.
  int32_t y; ///< Vertical, in 1/64 pixel, down.
} GSHGEOM_Device_Point;

/** The verbs of a device path.  The values are fixed. */
typedef enum GSHGEOM_Device_Verb {
  GSHGEOM_DEVICE_MOVE = 0,      ///< Start a subpath.
  GSHGEOM_DEVICE_LINE = 1,      ///< A line to a point.
  GSHGEOM_DEVICE_QUADRATIC = 2, ///< A quadratic: control point, end point.
  GSHGEOM_DEVICE_CUBIC = 3,     ///< A cubic: two control points, end point.
  GSHGEOM_DEVICE_CLOSE = 4,     ///< A line back to the subpath's start.
  GSHGEOM_DEVICE_VERB_COUNT     ///< Closes the enum.
} GSHGEOM_Device_Verb;

/**
 * One element read back.  `p` holds the control points and then the end point:
 * a move or a line has `p[0]`, a quadratic `p[0]` (control) and `p[1]` (end), a
 * cubic `p[0]`, `p[1]` (controls) and `p[2]` (end).  A close has `p[0]`, the
 * subpath's start.  `count` is the number of points in `p`.  The rest of `p`
 * is zero.
 */
typedef struct GSHGEOM_Device_Element {
  GSHGEOM_Device_Verb verb;      ///< The verb.
  unsigned count;                ///< How many entries of `p` are set.
  GSHGEOM_Device_Point p[3];     ///< Control points, then the end point.
} GSHGEOM_Device_Element;

/** The opaque device path. */
typedef struct GSHGEOM_Device_Path GSHGEOM_Device_Path;

/**
 * @brief The one conversion from `double` to 26.6: `floor(x * 64 + 0.5)`.
 *
 * @param x   The coordinate.
 * @param out Receives the result, only on success.
 * @return ::GSHGEOM_OK; ::GSHGEOM_ERR_INVALID for NULL, NaN or infinity;
 *         ::GSHGEOM_ERR_LIMIT if the result does not fit an `int32_t`.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_round_26_6(double x, int32_t * out);

/**
 * @brief The integer spelling of the rule (see the file comment):
 *        `floor((s * p * scale + 2^15) / 2^16)`, `s` = -1 when @p flip is 1.
 *
 * @param p     A coordinate in 26.6 of design units.
 * @param scale A 16.16 scale, pixels per design unit times 2^16.
 * @param flip  0, or 1 to negate (y up to y down).
 * @param out   Receives the 26.6 result, only on success.
 * @return ::GSHGEOM_OK; ::GSHGEOM_ERR_INVALID for NULL or a @p flip other than
 *         0 and 1; ::GSHGEOM_ERR_LIMIT when the result does not fit an
 *         `int32_t`.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_scale_26_6(
    int32_t p, int32_t scale, int flip, int32_t * out);

/**
 * @brief Apply a 2D map to a model point and round the result (the point of a
 *        device path).
 *
 * @return ::GSHGEOM_ERR_INVALID for NULL, an invalid or 3D map, or a
 *         non-finite point; ::GSHGEOM_ERR_LIMIT if either coordinate does not
 *         fit.  @p out is written only on success.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_device_point(const GSHGEOM_Map * map,
    GSHGEOM_Point p, GSHGEOM_Device_Point * out);

/**
 * @brief Create an empty device path.
 *
 * @param alloc  The allocator; NULL for cutil's default.
 * @param limits Caps (`max_elements`); NULL for the defaults.
 * @param out    Receives the path.
 * @return ::GSHGEOM_OK, ::GSHGEOM_ERR_INVALID for a NULL @p out,
 *         ::GSHGEOM_ERR_OOM.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_device_path_create(
    const GSHGEOM_Allocator * alloc, const GSHGEOM_Limits * limits,
    GSHGEOM_Device_Path ** out);

/** @brief Destroy a device path; NULL is ignored. */
GSHGEOM_API void gshgeom_device_path_destroy(GSHGEOM_Device_Path * path);

/** @brief The number of elements stored (0 for NULL). */
GSHGEOM_API size_t gshgeom_device_path_element_count(
    const GSHGEOM_Device_Path * path);

/**
 * @brief Start a subpath.
 *
 * The builders below are atomic: on any failure the path is as it was.
 *
 * @return ::GSHGEOM_OK, ::GSHGEOM_ERR_INVALID for NULL,
 *         ::GSHGEOM_ERR_LIMIT at `max_elements`, ::GSHGEOM_ERR_OOM.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_device_path_move_to(
    GSHGEOM_Device_Path * path, GSHGEOM_Device_Point p);

/**
 * @brief A line to a point.  Fails with ::GSHGEOM_ERR_INVALID before any move.
 *
 * After a close, a drawing verb starts a new subpath at the closed one's start:
 * the path stores a move there first, so every stored subpath begins with a
 * move, and the two elements are added together or not at all (a limit that
 * leaves room for one is ::GSHGEOM_ERR_LIMIT).
 */
GSHGEOM_API GSHGEOM_Result gshgeom_device_path_line_to(
    GSHGEOM_Device_Path * path, GSHGEOM_Device_Point p);

/** @brief A quadratic to @p p with control point @p c; see line_to. */
GSHGEOM_API GSHGEOM_Result gshgeom_device_path_quadratic_to(
    GSHGEOM_Device_Path * path, GSHGEOM_Device_Point c, GSHGEOM_Device_Point p);

/** @brief A cubic to @p p with control points @p c1 and @p c2; see line_to. */
GSHGEOM_API GSHGEOM_Result gshgeom_device_path_cubic_to(
    GSHGEOM_Device_Path * path, GSHGEOM_Device_Point c1,
    GSHGEOM_Device_Point c2, GSHGEOM_Device_Point p);

/**
 * @brief Close the subpath: a line back to its start.  Fails with
 *        ::GSHGEOM_ERR_INVALID before any move.  A second close in a row is a
 *        second close element.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_device_path_close(
    GSHGEOM_Device_Path * path);

/**
 * @brief Read element @p index back.
 *
 * @return ::GSHGEOM_OK, ::GSHGEOM_ERR_INVALID for NULL or an index past the
 *         end.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_device_path_element(
    const GSHGEOM_Device_Path * path, size_t index,
    GSHGEOM_Device_Element * out);

/**
 * @brief Make the device path of a model path under a map (CAP-3).
 *
 * Every point of every segment is mapped in `double` by ::gshgeom_map_apply_point and then
 * rounded by ::gshgeom_round_26_6.  A line, quadratic or cubic keeps its kind.  A rational curve
 * has no 26.6 form, and neither has an arc, which is read as its pieces (::gshgeom_arc_pieces,
 * rational quadratics): each is subdivided in `double` after the map, in homogeneous
 * coordinates at the parameter 1/2, until every interior control point of the projected piece
 * is within `t - 1` units of the piece's chord (`t` being ::GSHGEOM_FLATTEN_TOLERANCE, 4, so 3
 * units, 3/64 pixel), at most ::GSHGEOM_CONIC_MAX_DEPTH levels deep, and becomes a line to the
 * end of each piece (the ends of the whole curve are the mapped control points themselves, rounded
 * once, so an arc's device ends are its own ends and a closed arc closes).  The route is part of the
 * flattener's definition and its proof (flatten.h): every point of the exact curve
 * is within `t - 0.29` units of the polyline.  Every control point of every segment must round
 * into range.
 *
 * @param path   A 2D model path.
 * @param map    A 2D map (model to device pixels).
 * @param alloc  The allocator; NULL for the default.
 * @param limits Caps (`max_elements`); NULL for the defaults.
 * @param out    Receives the path, only on success.
 * @return ::GSHGEOM_OK; ::GSHGEOM_ERR_INVALID for NULL, a 3D path or map, or an
 *         invalid map; ::GSHGEOM_ERR_LIMIT for a coordinate out of range or
 *         too many elements; ::GSHGEOM_ERR_OOM.  Nothing is returned on
 *         failure.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_device_path_from_model(
    const GSHGEOM_Path * path, const GSHGEOM_Map * map,
    const GSHGEOM_Allocator * alloc, const GSHGEOM_Limits * limits,
    GSHGEOM_Device_Path ** out);

/**
 * @brief The FNV-1a 64 hash of a device path, in the spelling of hash.h (AD-21).
 *
 * `u32` format version (1), `u32` kind (`0x50564544`, "DEVP" little-endian),
 * `u64` element count, and for each element `u8` verb followed, for every verb
 * but a close, by its points in readback order (control points, then the end),
 * each as `i32` x then `i32` y.  NULL hashes as the empty stream.
 */
GSHGEOM_API uint64_t gshgeom_device_path_hash(const GSHGEOM_Device_Path * path);

/**
 * @brief Write the dump: a header line `gshgeom.device-path 1 elements=N
 *        hash=<16 hex digits>`, then one line per element, a letter (`M`,
 *        `L`, `Q`, `C`, `Z`) and its points as `x y` pairs, each number the
 *        two's complement of the `int32_t` in 8 hexadecimal digits.
 *
 * @return ::GSHGEOM_OK, ::GSHGEOM_ERR_INVALID for NULL, ::GSHGEOM_ERR_IO when
 *         the stream fails (including at the flush).
 */
GSHGEOM_API GSHGEOM_Result gshgeom_device_path_dump(
    const GSHGEOM_Device_Path * path, FILE * out);

#ifdef __cplusplus
}
#endif

#endif /* GHOTI_IO_GSHGEOM_DEVICE_H */
