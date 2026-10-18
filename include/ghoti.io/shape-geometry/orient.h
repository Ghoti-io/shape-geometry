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
 * @file orient.h
 *
 * The orientation rules of the shape family, stated once, and the exact
 * predicates that decide with them (AD-6, CAP-2, CAP-4).
 *
 * **The rules.**  Model space is `double`, x to the right and y down.  Every
 * rule below is a statement about that space, and no other library restates
 * one.  Their values follow SVG and canvas:
 *
 * | Rule | Value |
 * |------|-------|
 * | Positive angle | turns +x toward +y (clockwise on a screen) |
 * | Positive winding, positive signed area | the same sense: the boundary goes from +x toward +y around the inside |
 * | `orient2d(a, b, c)` | positive when `a -> b -> c` turns that way, equivalently when `c` is on the positive side of the line `a -> b` |
 * | Arc sweep | a positive sweep goes that way (SVG `sweep-flag` 1) |
 * | Left of a tangent | the tangent turned by a negative angle: for a tangent (1, 0), left is (0, -1), up on the screen |
 * | Right of a tangent | the tangent turned by a positive angle |
 * | Positive offset | grows a positively wound region: the boundary moves along its left normal (the outward side of a positively wound boundary) |
 * | Sweep gradient | advances in the positive angle direction |
 * | `orient3d(a, b, c, d)` | the sign of the determinant of `b-a, c-a, d-a`, positive when they are a positively oriented basis of (x, y, z); with `d - a = (0, 0, 1)` it is `orient2d` of the xy projections |
 *
 * Arcs, offsets and gradients are not built yet; their rule is stated here so
 * that the first library to build one does not choose again, and the chiral
 * fixture in the tests pins the rules that have code (angles, the signed area,
 * `orient2d`, `orient3d`, the two normals).
 *
 * **The predicates** return the exactly correct sign for every finite input.
 * Where plain `double` arithmetic would get the sign wrong (a nearly
 * collinear triple, a nearly coplanar quadruple), the answer is still exact.
 * A filter decides the clear cases at the cost of a few multiplications; the
 * rest are decided in integer arithmetic on the mantissas and exponents.
 * `..._exact` always takes the exact path (for tests and audits).
 *
 * Copyright 2026 by Corey Pennycuff
 */

#ifndef GHOTI_IO_GSHGEOM_ORIENT_H
#define GHOTI_IO_GSHGEOM_ORIENT_H

#include <ghoti.io/shape-geometry/core.h>
#include <ghoti.io/shape-geometry/curve.h>
#include <ghoti.io/shape-geometry/macros.h>

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * The orientation of three points in the plane: the sign of
 * `(b.x-a.x)(c.y-a.y) - (b.y-a.y)(c.x-a.x)`, computed exactly.
 *
 * @param a,b,c The points; `z` is ignored.
 * @param sign  Receives -1, 0 or +1: negative, collinear, positive (see the
 *              rules above).
 * @return ::GSHGEOM_OK, or ::GSHGEOM_ERR_INVALID for a NULL `sign` or a
 *         non-finite coordinate (`sign` is then untouched).
 */
GSHGEOM_API GSHGEOM_Result gshgeom_orient2d(
    GSHGEOM_Point a, GSHGEOM_Point b, GSHGEOM_Point c, int * sign);

/** As ::gshgeom_orient2d, always by the exact path. */
GSHGEOM_API GSHGEOM_Result gshgeom_orient2d_exact(
    GSHGEOM_Point a, GSHGEOM_Point b, GSHGEOM_Point c, int * sign);

/**
 * The orientation of four points in space: the sign of the determinant of the
 * rows `b-a`, `c-a`, `d-a`, computed exactly.  It is positive when `d` lies on
 * the side of the plane through `a`, `b`, `c` that the right-hand rule gives
 * for `a -> b -> c` taken as `orient2d` reads it, in the axes (x, y, z); for
 * points with equal `z` and `d` above them (+z) it equals `orient2d` of the xy
 * projections.  (It is the negative of the determinant of `a-d`, `b-d`,
 * `c-d`, which is Shewchuk's sign.)
 *
 * @return As ::gshgeom_orient2d.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_orient3d(GSHGEOM_Point a, GSHGEOM_Point b,
    GSHGEOM_Point c, GSHGEOM_Point d, int * sign);

/** As ::gshgeom_orient3d, always by the exact path. */
GSHGEOM_API GSHGEOM_Result gshgeom_orient3d_exact(GSHGEOM_Point a,
    GSHGEOM_Point b, GSHGEOM_Point c, GSHGEOM_Point d, int * sign);

/**
 * The left-hand normal of a tangent: the tangent turned by exactly a negative
 * quarter turn, `(t.y, -t.x)`, not normalized.  `z` is dropped.
 */
GSHGEOM_API GSHGEOM_Point gshgeom_normal_left(GSHGEOM_Point tangent);

/** The right-hand normal: `(-t.y, t.x)`, the opposite of the left. */
GSHGEOM_API GSHGEOM_Point gshgeom_normal_right(GSHGEOM_Point tangent);

/**
 * The signed area of the closed polygon through `points` (the last joined to
 * the first), positive when its boundary winds in the positive sense.  An
 * ordinary floating-point sum in a fixed order (a fan from the first point),
 * not a predicate: use ::gshgeom_orient2d to decide.
 *
 * @param points The vertices; `z` is ignored.
 * @param count  The number of vertices; fewer than 3 give 0.
 * @param out    Receives the area.
 * @return ::GSHGEOM_OK; ::GSHGEOM_ERR_INVALID for NULL or a non-finite
 *         coordinate; ::GSHGEOM_ERR_LIMIT if the area is not finite.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_polygon_signed_area(
    const GSHGEOM_Point * points, size_t count, double * out);

#ifdef __cplusplus
}
#endif

#endif /* GHOTI_IO_GSHGEOM_ORIENT_H */
