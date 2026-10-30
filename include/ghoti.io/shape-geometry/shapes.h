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
 * @file shapes.h
 *
 * The named primitives of a 2D path (rect, rounded rect, ellipse), defined exactly (PR2 to PR5, PR13).
 *
 * A front-end that knows it is drawing a rectangle, a rounded rectangle or an ellipse says so by calling one of these builders instead of
 * deriving control points itself, and so every front-end gets the same elements: lines and quarter-turn arcs (arc.h), the start point of
 * SVG 2's equivalent paths, and the positive direction of AD-6 (clockwise on a y-down screen).  A builder begins a new subpath at the start
 * point and closes it.  The builders add geometry only: the primitive tag, which says "this is what the front-end meant", lives on a
 * path resource in the scene library, never on a path.
 *
 * ### The box
 *
 * `(x, y)` is the top-left corner of the box, `w` and `h` its width and height (not negative; zero is allowed: a rectangle of no
 * width is a degenerate but legitimate fact, and Canvas strokes one).  Every coordinate a builder derives is the correctly rounded
 * value of the real expression named (`x + w`, `x + w - rx`, `cx - rx`): a sum is rounded once, from the exact sum.  A sum with a
 * zero term is the other term, bit for bit (`x + 0` is `x`, so a point at `-0` stays), and a sum that is exactly zero is `+0`.
 *
 * ### The flips
 *
 * `flip_x` and `flip_y` (each 0 or 1) mirror the shape in the box's vertical and horizontal axis, the start corner and the turn
 * following the mirror.  They are Canvas's four sign cases of `rect(x, y, w, h)`: no flip is `w >= 0, h >= 0` (top-left, clockwise);
 * `flip_x` is `w < 0` (top-right, counter-clockwise); `flip_y` is `h < 0` (bottom-left, counter-clockwise); both is `w < 0, h < 0`
 * (bottom-right, clockwise).  An ellipse follows: no flip starts at angle 0 and turns the positive way; `flip_y` alone starts at 0
 * and turns the other way; `flip_x` alone starts at pi and turns the other way; both start at pi and turn the positive way.  The
 * direction is the parity of the flips.  SVG uses no flip.  A rounded rectangle takes flips too, with a different path: see
 * ::gshgeom_path_add_rounded_rect.
 *
 * ### Failure
 *
 * Every argument is checked before the path is touched, and a failure part way (the path's element limit, an allocation) takes back
 * what the call added: a failed call leaves the path as it was.
 */

#ifndef GHOTI_IO_GSHGEOM_SHAPES_H
#define GHOTI_IO_GSHGEOM_SHAPES_H

#include <ghoti.io/shape-geometry/allocator.h>
#include <ghoti.io/shape-geometry/arc.h>
#include <ghoti.io/shape-geometry/core.h>
#include <ghoti.io/shape-geometry/macros.h>
#include <ghoti.io/shape-geometry/path.h>

#ifdef __cplusplus
extern "C" {
#endif

/** The radii of one corner of a rounded rectangle: the x and y semi-axes of its quarter ellipse. */
typedef struct GSHGEOM_Corner {
  double rx; ///< The horizontal radius.
  double ry; ///< The vertical radius.
} GSHGEOM_Corner;

/** The most elements a rect adds (a move, three lines, a close). */
#define GSHGEOM_SHAPE_RECT_ELEMENTS 5u
/** The most elements a rounded rectangle adds (a move, four lines, four arcs, a close). */
#define GSHGEOM_SHAPE_ROUNDED_RECT_ELEMENTS 10u
/** The elements an ellipse adds (a move, four arcs, a close). */
#define GSHGEOM_SHAPE_ELLIPSE_ELEMENTS 6u

/**
 * @brief Whether @p r are radii a rounded rectangle of @p w by @p h can have (PR2).
 *
 * @p r holds the corners in the order top-left, top-right, bottom-right, bottom-left.  A corner with a zero component is square,
 * whatever its other component, and counts as `{0, 0}`.  The curves of two corners must not overlap on any side: the top
 * (`r[0].rx + r[1].rx > w`), the right (`r[1].ry + r[2].ry > h`), the bottom (`r[2].rx + r[3].rx > w`) or the left
 * (`r[3].ry + r[0].ry > h`).  Each is decided **exactly**: the sum of two doubles is compared with a double in the reals, never on a
 * rounded sum, so a pair that sums to the side exactly is valid (a stadium, a circle) and one an ulp over is not, however the
 * sum rounds.  The waist does not clamp or scale radii (D5).
 *
 * @return ::GSHGEOM_OK, or ::GSHGEOM_ERR_INVALID for a NULL, a non-finite number, a negative width, height or radius component, or an
 *   overlap.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_rounded_rect_check(double w, double h, const GSHGEOM_Corner r[4]);

/**
 * @brief A rectangle (PR3, PR13): `M`, three `L`, `Z`.
 *
 * No flip: `(x, y)`, `(x + w, y)`, `(x + w, y + h)`, `(x, y + h)`.  `flip_x` exchanges `x` with `x + w`, `flip_y` exchanges `y` with
 * `y + h`.
 *
 * @return ::GSHGEOM_OK; ::GSHGEOM_ERR_INVALID for a NULL or 3D path, a flip that is not 0 or 1, a non-finite number or a negative
 *   `w` or `h`; ::GSHGEOM_ERR_LIMIT for a coordinate past ::GSHGEOM_COORD_LIMIT or the path's limit; ::GSHGEOM_ERR_OOM.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_path_add_rect(
    GSHGEOM_Path * path, double x, double y, double w, double h, int flip_x, int flip_y);

/**
 * @brief A rounded rectangle (PR4): a move, four lines and up to four quarter-turn arcs, closed.
 *
 * @p r are the radii in the order the sign case visits the corners (below), valid by ::gshgeom_rounded_rect_check.  With no flip that is
 * top-left, top-right, bottom-right, bottom-left: the path starts after the top-left radius, `(x + r[0].rx, y)`, and turns the positive
 * way.  The four lines are **always** emitted, even when a side's radii sum to the side and a line has no length (the structure is a
 * function of which corners are round, nothing else); an arc is emitted for a round corner only.  Each arc is a quarter turn in degrees
 * (so its control polygon is exact), its ends the named points to the bit, so the subpath closes to the bit; it is circular iff `rx` and
 * `ry` are bit equal.  With no flip and no round corner the path is the rectangle's.
 *
 * ### The flips: Canvas `roundRect` with a negative width or height, as Chromium and Firefox draw it
 *
 * `flip_x` is a negative width and `flip_y` a negative height; `(x, y, w, h)` is the normalized box.  A flip does two things, and **not**
 * what it does to `rect` (the path does not start at the mirrored corner):
 *  - the radii go to the corners the mirror gives them: `r[0]` to the corner `rect` would start at.  `r` is top-left, top-right,
 *    bottom-right, bottom-left for no flip; top-right, top-left, bottom-left, bottom-right for `flip_x`; bottom-left, bottom-right,
 *    top-right, top-left for `flip_y`; bottom-right, bottom-left, top-left, top-right for both;
 *  - the path is the clockwise path of the box with those radii when the flips are equal (it starts after the box's top-left radius), and
 *    the **reversal** of it when they differ: the same start `(x + tl.rx, y)`, the top-left arc and the left edge first, then the bottom,
 *    the right edge and the top-right arc, the arcs running the other way (sweep -90).
 * A flipped rounded rectangle with no round corner still has the four lines and is a rounded rectangle, not a rect.  Firefox agrees with
 * Chromium for an even number of flips and starts an odd one at `(x + w - tr.rx, y)` instead, which this does not follow.
 *
 * @return As ::gshgeom_path_add_rect, and ::GSHGEOM_ERR_INVALID for radii ::gshgeom_rounded_rect_check refuses.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_path_add_rounded_rect(GSHGEOM_Path * path, double x, double y, double w, double h,
    const GSHGEOM_Corner r[4], int flip_x, int flip_y);

/**
 * @brief An ellipse (PR5, PR13): a move and four quarter-turn arcs about `(cx, cy)`, closed.
 *
 * No flip: from `(cx + rx, cy)` ("3 o'clock"), arcs through `(cx, cy + ry)`, `(cx - rx, cy)`, `(cx, cy - ry)` and back, each a positive
 * quarter turn in degrees starting at 0, 90, 180 and 270, with the half-diameters `u = (rx, 0)` and `v = (0, ry)`.  `flip_x` negates
 * `u`, `flip_y` negates `v`; the start and the other ends follow.  Circular iff `rx` and `ry` are bit equal.
 *
 * @return As ::gshgeom_path_add_rect, with ::GSHGEOM_ERR_INVALID for an `rx` or `ry` that is not positive.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_path_add_ellipse(
    GSHGEOM_Path * path, double cx, double cy, double rx, double ry, int flip_x, int flip_y);

#ifdef __cplusplus
}
#endif

#endif /* GHOTI_IO_GSHGEOM_SHAPES_H */
