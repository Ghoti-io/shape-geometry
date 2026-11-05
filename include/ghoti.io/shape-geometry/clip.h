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
 * @file clip.h
 *
 * Clipping a flattened path to a convex quadrilateral, keeping the winding.
 *
 * ::gshgeom_flat_path_clip_convex cuts the contours of a ::GSHGEOM_Flat_Path (26.6 points) to a convex quadrilateral
 * `Q` of 26.6 points and returns them as a new flat path.  It is built for a rasterizer or a pattern that draws one
 * shape many times and needs each copy cut at the edge of its cell: the pieces of two neighbours must add to the whole,
 * so the clip changes the geometry, and it must not change what the geometry means.  That is the first rule (CK2).
 *
 * ## The rules
 *
 * - **CK1, the function and what it refuses.**  `Q` is four points in order.  It must be convex and have an area: the
 *   four turns (the sign of ::gshgeom_orient2d_exact at each corner) are all of one sign or zero, they are not all zero,
 *   and the area is not zero; otherwise the call is ::GSHGEOM_ERR_INVALID.  Either orientation is accepted, and points
 *   may repeat (a triangle given as four points).  When nothing of the path is inside `Q` the call succeeds and the
 *   result has no subpath: that is the status "empty", not an error, and there is no separate "empty" result code
 *   (::gshgeom_flat_path_cancel_seams does the same).
 *   ::GSHGEOM_ERR_LIMIT is returned when the result would hold more than `max_points` points.
 * - **CK2, winding is preserved.**  Every subpath is read as a closed polygon, as a rasterizer reads it.  For every point
 *   `q` that is not within the rounding bound of CK4 of an edge of the input or of `Q`, the winding number of the result
 *   at `q` equals the input's if `q` is in the interior of `Q`, and is 0 outside it.  This holds for any closed
 *   polylines: self-intersecting, nested and overlapping contours included.  Each contour is cut by one half-plane of
 *   `Q` at a time, and the part of the contour outside the half-plane is replaced by the straight run along the clip line
 *   between the point where the contour left and the point where it came back; the replaced part and the run form a closed
 *   curve outside the half-plane, which winds zero times about every point inside it.  No contour is dropped (except one
 *   that is wholly outside, or whose cut leaves fewer than three points and so encloses nothing), merged with another,
 *   re-oriented or unioned.  A contour wholly inside `Q` is copied bit for bit, its `closed` flag included; a cut
 *   contour is a closed polygon and is flagged closed.
 * - **CK3, exact decisions, one rounding.**  "Inside" is the closed quadrilateral: a point exactly on a side of `Q` is
 *   inside, so a contour that lies along the boundary is kept as it is.  Which side of a clip line a point lies on is
 *   decided exactly on the integer coordinates (a determinant of two 64-bit products, no `double`).  An edge of the contour yields a
 *   crossing only when its end points lie strictly on opposite sides of a clip line; an end point on the line is kept
 *   as it is and is its own crossing.  The crossing is the exact rational intersection of the edge with the line
 *   through two points of `Q`, rounded once to 26.6 to nearest with ties upward (`floor(v + 1/2)` in each coordinate), by
 *   exact integer arithmetic of 128 bits (no `double`, no `long double`, no libm, no fused multiply-add).  It is
 *   computed with the edge's end points in lexicographic order and the line's points in lexicographic order, so the
 *   same edge and line give the same point whichever way the edge is traversed or the line is named.  The result is the
 *   same bits on every target.
 * - **CK4, the seam bound.**  Each crossing is within half a unit (1/128 pixel) of the exact one in each coordinate.
 *   Later cuts decide on the rounded points, so a point of the result can lie outside `Q` by the rounding, and every edge of
 *   the result lies within three units of an edge of the input or of `Q`; CK2 excepts that distance.  The proved bound on
 *   the area of the symmetric difference between the result and the exact clip of the input (where the windings differ,
 *   weighted by the difference) is `c (sqrt(2)/2) D` square units for `c` crossings and a quadrilateral of diameter `D`:
 *   a crossing moves by at most `sqrt(2)/2` units, and the two result edges that meet there are at most `D` long.  There
 *   is no bound proportional to the length `P` of the clip lines crossed: an edge that is long against a line moves its
 *   end by the rounding over its whole length, so a thin cell crossed at a shallow angle by a long edge exceeds
 *   `P/128` pixels squared by a large factor, and over the tests' random sweep the area reached 1.40 times that
 *   figure.  `P/128` is therefore not a guarantee; the bound above is.  (documentation/design.md, "Clipping a
 *   flattened path to a convex quadrilateral", gives the proof and the measured cases.)
 * - **CK5, size.**  The result has at most `n + c/2` points, so at most `n + c`, for `n` input points and `c` crossings
 *   computed (each cut removes at least one point for every two crossings it makes), and at most `(3/2)^4 n`, a little
 *   over `5n`.  The crossings are at most `65 n / 8 = 8.125 n` in all (a cut makes at most one crossing per point of
 *   the polygon it cuts).  That `c <= 4 n` is **not proved**: a search over polygons reached `c = 4 n` exactly, for
 *   every even `n` it tried, and never more, which is evidence and not a theorem.  Each cut is `O(n)` in the points of
 *   the polygon it cuts.
 * - **CK6, determinism.**  The result, and so ::gshgeom_flat_path_hash of it, is a function of the input and `Q` alone:
 *   the same on every target, in the cross-target probe.
 *
 * The cuts are made in the order of `Q`'s sides (the side from `quad[0]` to `quad[1]` first), and a side whose two
 * points are the same is skipped.
 */

#ifndef GHOTI_IO_GSHGEOM_CLIP_H
#define GHOTI_IO_GSHGEOM_CLIP_H

#include <ghoti.io/shape-geometry/allocator.h>
#include <ghoti.io/shape-geometry/core.h>
#include <ghoti.io/shape-geometry/device.h>
#include <ghoti.io/shape-geometry/flatten.h>
#include <ghoti.io/shape-geometry/macros.h>

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * The revision of the clip's construction: it changes when the points a clip produces for the same input change, so that a
 * sink that keeps what a clip made (a cached picture of a tile) can tell the bits are no longer the ones it holds.
 */
#define GSHGEOM_CLIP_REVISION 1

/**
 * @brief Clip a flattened path to a convex quadrilateral, preserving the winding (CK1 to CK6 above).
 *
 * The signature is that of ::gshgeom_flat_path_cancel_seams: the path, the quadrilateral, an allocator, the limits and
 * the result.  The result is `OK` with no subpath when nothing is inside, so a caller tests the subpath count.
 *
 * @param path   The flattened path; every subpath is read as closed whatever its flag.  Not changed.
 * @param quad   The four corners of `Q`, in order, in 26.6 (either orientation; repeats allowed).
 * @param alloc  The allocator for the result and the scratch, or NULL for the default.
 * @param limits The limits (`max_points` bounds the result), or NULL for the defaults.
 * @param out    The result, to be destroyed with ::gshgeom_flat_path_destroy; written only on success.  It has no subpath
 *               when nothing of @p path is inside `Q`.
 * @return ::GSHGEOM_OK (an empty result included); ::GSHGEOM_ERR_INVALID for a NULL @p path, @p quad or @p out, or a
 *         quadrilateral that is not convex or has no area; ::GSHGEOM_ERR_LIMIT if the result would hold more than
 *         `max_points` points; ::GSHGEOM_ERR_OOM.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_flat_path_clip_convex(const GSHGEOM_Flat_Path * path,
    const GSHGEOM_Device_Point quad[4], const GSHGEOM_Allocator * alloc, const GSHGEOM_Limits * limits,
    GSHGEOM_Flat_Path ** out);

#ifdef __cplusplus
}
#endif

#endif /* GHOTI_IO_GSHGEOM_CLIP_H */
