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
 * @file vertices.h
 *
 * The vertices of a path and the directions at them (CAP-5, CAP-13; AD-6): where a
 * marker goes and which way it points.
 *
 * A *vertex* of a 2D path is the end point of one of its elements: a move gives
 * the point it moves to, a line or a curve its last control point, and a close
 * the start of the subpath it closes.  They come in element order.  A drawing
 * verb that follows a close starts a subpath at the closed one's start without a
 * move (stroke.h, SG5) and so adds no vertex of its own.
 *
 * ### Directions
 *
 * A segment of a subpath is a line, a curve or a close.  It **has extent** when
 * some control point differs, exactly, from its first (SG5).  Its **start
 * direction** is `P_k - P_0` for the least `k` with `P_k != P_0` and its **end
 * direction** `P_n - P_k` for the greatest `k` with `P_k != P_n`, the two
 * directions the stroke expansion joins and caps with (the same code reads
 * both).  A *unit* direction is that vector divided by its length, with the
 * length computed after scaling by the larger coordinate, so that neither a huge
 * nor a tiny vector overflows or underflows.
 *
 * The **incoming** direction of a vertex is the end direction of the nearest
 * segment with extent at or before the vertex in its subpath, and the
 * **outgoing** direction the start direction of the nearest one after it: a
 * segment without extent has no direction of its own and the vertex looks past
 * it, to either side, so a run of zero-length segments is invisible to the
 * directions of every vertex around and inside it.  A vertex with no such
 * segment on a side has no direction on that side.
 *
 * In a *closed* subpath (one that ends in a close), the vertex that starts it and
 * the vertex of the close, which are the same point, each take the end
 * direction of the subpath's last segment with extent as incoming and the start
 * direction of its first as outgoing: the corner at the start of a closed
 * subpath is a corner like any other.
 *
 * ### The marker direction
 *
 * `direction` is the direction an `orient="auto"` marker takes (SVG 2, markers):
 * the outgoing direction at a vertex with no incoming one, the incoming at one
 * with no outgoing, the bisector at one with both, and `(1, 0)` (angle zero) at
 * one with neither.  The bisector is the unit vector of the sum of the unit
 * incoming and outgoing directions; where they are exactly opposite, so that the
 * sum is the zero vector (a reversal), it is the incoming direction turned by a
 * quarter turn, ::gshgeom_normal_right of it (SVG leaves this case open).  The angle is
 * left to the caller: `gmath_atan2` of the components is bit-exact on every target.
 *
 * ### Exactness
 *
 * Every number is a product of IEEE double operations the library defines (a
 * difference, a square root, a quotient, a scale by a power of two), computed
 * without contraction, so a vertex has the same bits on every target the library
 * is built for (the committed corpus hash reads them).
 */

#ifndef GHOTI_IO_GSHGEOM_VERTICES_H
#define GHOTI_IO_GSHGEOM_VERTICES_H

#include <ghoti.io/shape-geometry/core.h>
#include <ghoti.io/shape-geometry/curve.h>
#include <ghoti.io/shape-geometry/macros.h>
#include <ghoti.io/shape-geometry/path.h>

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief What is known of a vertex.  Bit flags of ::GSHGEOM_Vertex.flags. */
typedef enum GSHGEOM_Vertex_Flag {
  GSHGEOM_VERTEX_HAS_IN = 1u << 0,        ///< `in` is set: some segment with extent lies at or before it in its subpath.
  GSHGEOM_VERTEX_HAS_OUT = 1u << 1,       ///< `out` is set: some segment with extent lies after it in its subpath.
  GSHGEOM_VERTEX_PATH_START = 1u << 2,    ///< The first vertex of the path (where a start marker goes).
  GSHGEOM_VERTEX_PATH_END = 1u << 3,      ///< The last vertex of the path (where an end marker goes).
  GSHGEOM_VERTEX_SUBPATH_START = 1u << 4, ///< The move that starts a subpath.
  GSHGEOM_VERTEX_SUBPATH_END = 1u << 5,   ///< The last vertex of its subpath.
  GSHGEOM_VERTEX_CLOSED = 1u << 6         ///< Its subpath ends in a close.
} GSHGEOM_Vertex_Flag;

/** @brief One vertex of a path. */
typedef struct GSHGEOM_Vertex {
  size_t element;        ///< The index of the path element that ends at this vertex.
  size_t subpath;        ///< The index of its subpath, counted from 0 in element order (a drawing verb after a close starts the next).
  uint32_t flags;        ///< ::GSHGEOM_Vertex_Flag bits.
  GSHGEOM_Point point;   ///< The position.
  GSHGEOM_Point in;      ///< The unit incoming direction; zero unless ::GSHGEOM_VERTEX_HAS_IN.
  GSHGEOM_Point out;     ///< The unit outgoing direction; zero unless ::GSHGEOM_VERTEX_HAS_OUT.
  GSHGEOM_Point direction; ///< The direction of an `orient="auto"` marker: a unit vector, always set.
} GSHGEOM_Vertex;

/**
 * @brief The callback of ::gshgeom_path_vertices.
 *
 * @param vertex The vertex, valid for the call only.
 * @param ctx The caller's pointer.
 * @return ::GSHGEOM_OK to go on; anything else stops the walk and is returned from it.
 */
typedef GSHGEOM_Result (*GSHGEOM_Vertex_Fn)(const GSHGEOM_Vertex * vertex, void * ctx);

/**
 * @brief The number of vertices of a path: its elements that end at a point.
 *
 * @param path The path (2D).
 * @param count The number; written only on success.
 * @return ::GSHGEOM_OK or ::GSHGEOM_ERR_INVALID (NULL, or a 3D path).
 */
GSHGEOM_API GSHGEOM_Result gshgeom_path_vertex_count(const GSHGEOM_Path * path, size_t * count);

/**
 * @brief Give each vertex of a path to @p fn, in element order.
 *
 * Nothing is allocated, and the work is linear in the number of elements (a run of
 * segments without extent is passed over once, not once per vertex in it).  A
 * path with no element gives no call.
 *
 * @param path The path (2D).
 * @param fn The callback.
 * @param ctx Passed to @p fn.
 * @return ::GSHGEOM_OK, ::GSHGEOM_ERR_INVALID (NULL, or a 3D path), or what @p fn returned when it stopped the walk.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_path_vertices(const GSHGEOM_Path * path, GSHGEOM_Vertex_Fn fn, void * ctx);

#ifdef __cplusplus
}
#endif

#endif /* GHOTI_IO_GSHGEOM_VERTICES_H */
