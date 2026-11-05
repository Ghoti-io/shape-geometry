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
 * @file lattice.h
 *
 * The copies of a repeated cell: which copies of a lattice meet a rectangle, where each is translated to, and the
 * integer period of a lattice of whole-pixel steps.  A pattern, a hatch or any picture drawn on a lattice of copies uses
 * these so that every renderer draws the same set of copies and forms each copy's position with the same bits.
 *
 * A **lattice** is a cell `[x, x + w] x [y, y + h]` in a space `S` and two steps: the copy `(i, j)` is the cell
 * translated by `(i x_step, j y_step)`.  A map `M` takes `S` to the device.  `repeat_x` (`repeat_y`) of 0 keeps only
 * index 0 on that axis.
 *
 * ## The rules
 *
 * - **LT1, the copies of a use.**  ::gshgeom_lattice_each reports, in **lattice order** (increasing `j`, then
 *   increasing `i`), the copies `(i, j)` whose device quadrilateral meets the **closed** rectangle of the use (touching
 *   counts).  The quadrilateral of the copy is the convex hull of the four points `Q_k + tau + i A + j B`, taken as exact
 *   real numbers, where `Q_k` are the corners of the cell of copy `(0, 0)` under the linear part of `M`
 *   (::gshgeom_lattice_quad), `tau` is the translation of `M` and `A`, `B` are the device step vectors of
 *   ::gshgeom_lattice_translation, all of them doubles.  The decision is exact: no comparison is made in `double`.  The
 *   copy's drawn position is the rounded one of LT2 (::gshgeom_lattice_translation); it differs from the exact one by a
 *   few units in the last place of its magnitude, so a copy that touches the rectangle at zero area may be drawn at a
 *   position that no longer touches it, and one is never reported for that reason alone.  The set, and so its count, is a
 *   function of the doubles and the rectangle only.
 *   - *Rows are intervals.*  Fix `j`.  The copies of the row are the `i` with `w(i) = tau + j B + i A` in
 *     `K = R - hull(Q)` (the Minkowski difference of the rectangle `R` and the hull of the quadrilateral), a convex
 *     polygon, since a translate of the hull meets `R` exactly when its translation lies in `K`.  The points `w(i)` lie
 *     on a line at equal spacing, so the `i` for which they lie in a convex set are consecutive integers.  A row is
 *     reported as `i_lo .. i_hi`; the rows with copies are *not* always consecutive in `j` (a thin region can hold lattice
 *     points in rows `0` and `2` and none in row `1`), and a row with none is not reported.
 *   - *How the interval is found.*  `K` is the intersection of the half-planes `u . w <= h(u)` for `u` in the sixteen
 *     directions made of the two axes (both signs) and the normals (both signs) of the six segments between the four
 *     points of `Q` (a hull edge is among them), where `h(u) = max over the corners of R of u . r  -  min over k of u .
 *     Q_k`.  Each is `c i + e j + g <= 0` with `c = u . A`, `e = u . B` and `g = u . tau - u . r + u . Q_k` exact sums of
 *     products of doubles.  For one row each constraint bounds `i` on one side; the bound is estimated in `double`, and
 *     then **fixed exactly**: the predicate `c i + e j + g <= 0` is evaluated exactly (a sign filter decides the clear
 *     cases, the exact sum the rest) and is monotone in `i`, so galloping from the estimate in steps of 1, 2, 4 and then
 *     bisecting finds the exact integer bound whatever the error of the estimate (an error of the estimate costs steps and
 *     never a wrong answer).  The rows are the integers `j` of the strip between the two supporting lines of `K` parallel
 *     to `A`, found the same way from the direction `(-A_y, A_x)`, so the first and last row are exact.
 *   - *Cost.*  A row costs a constant number of constraint evaluations and logarithmically many more only when an estimate
 *     is far off.  The count is the sum of the row lengths, found without visiting the copies, so a lattice of `2^-64`
 *     steps under a full-target rectangle reports ::GSHGEOM_ERR_LIMIT on its first row.  The rows visited, empty ones
 *     included, are bounded by `2 max_count + GSHGEOM_LATTICE_ROW_WORK`: a lattice so sparse that its rows are mostly
 *     empty is ::GSHGEOM_ERR_LIMIT as well, decided the same way on every target (the count in that case is `max_count + 1`).
 *   - *Index range.*  Indices are in `[-GSHGEOM_LATTICE_INDEX_MAX, GSHGEOM_LATTICE_INDEX_MAX]`, the range in which an index is
 *     an exact double; a copy outside it is not a copy of the lattice.
 * - **LT2, the position of a copy.**  ::gshgeom_lattice_translation is `((tau + i a) + j b)` per component, `a =
 *   (fl(m00 x_step), fl(m10 x_step))` and `b = (fl(m01 y_step), fl(m11 y_step))` formed once, `i` and `j` converted to
 *   double exactly, each operation rounded once and evaluated left to right, no fused multiply-add: nothing accumulates, so
 *   the position of copy `(i, j)` carries no error from the copies before it.  ::gshgeom_lattice_cell is the cell of
 *   copy `(i, j)` in the source space: `left = x + i x_step`, `right = x + (i + 1) x_step` when `x_step` equals `w` (so the
 *   edge shared by two abutting cells is one double) and `left + w` otherwise, likewise `top` and `bottom`.
 * - **LT3, the period.**  ::gshgeom_lattice_period turns two integer step vectors into the Hermite normal form `(p, r, s)`
 *   of the lattice `L = Z a + Z b`: `L` is generated by `(p, 0)` and `(r, s)`, `p s = |D|`, `0 <= r < p`.  The reduction
 *   `rep(x, y) = ((x - k r) mod p, y - k s)`, `k = floor(y / s)`, takes every integer point to the one point of the window
 *   `[0, p) x [0, s)` that is congruent to it modulo `L`: `rep(q + lambda) = rep(q)` for every `lambda` in `L`, and the window
 *   holds exactly one point of each of the `|D|` classes.  Integer arithmetic only.
 * - **LT4, determinism.**  Nothing is decided in `double` but the formulas LT2 states; there is no libm, no `long double`, no
 *   fused multiply-add and no x87 instruction, so the results are the same bits on every target.
 * - **LT5, the canonical position of a copy.**  A sink that draws the copies translates a copy by whole pixels and draws it at
 *   the fraction, so that a copy is a function of the linear part, the fraction and the shift (the split of
 *   ::gshgeom_stroke_canonical_map).  ::gshgeom_lattice_canonical returns that split for the copy `(i, j)`.  When the device step
 *   vectors `a` and `b` are **integer vectors** (every component an exact integer of magnitude below `2^31`; see
 *   ::gshgeom_lattice_integer_steps) the whole pixels are `floor(tau) + i a + j b` in integers and the fraction is the fraction of
 *   `tau`, for every copy: no sum is rounded, so all the copies are one set of bits shifted, which is what a periodic lattice
 *   needs (LT2's `(tau + i a) + j b` in `double` rounds `tau`'s fraction once the sum grows a bit, so two copies would differ in
 *   their last bits).  Otherwise the position is LT2's, split by ::gshgeom_stroke_canonical_map.  A whole-pixel coordinate of
 *   magnitude `2^40` or more is ::GSHGEOM_ERR_LIMIT.
 */

#ifndef GHOTI_IO_GSHGEOM_LATTICE_H
#define GHOTI_IO_GSHGEOM_LATTICE_H

#include <ghoti.io/shape-geometry/core.h>
#include <ghoti.io/shape-geometry/macros.h>
#include <ghoti.io/shape-geometry/map.h>

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/** The largest index magnitude: `2^53 - 1`, so that an index and its successor are exact doubles. */
#define GSHGEOM_LATTICE_INDEX_MAX ((((int64_t)1) << 53) - 1)

/** The rows ::gshgeom_lattice_each may visit beyond `2 max_count` before it gives up with ::GSHGEOM_ERR_LIMIT (LT1, cost). */
#define GSHGEOM_LATTICE_ROW_WORK ((uint64_t)65536)

/** The largest determinant magnitude ::gshgeom_lattice_period accepts: `2^31`. */
#define GSHGEOM_LATTICE_PERIOD_MAX (((int64_t)1) << 31)

/** The revision of the lattice's positions (LT2, LT5): it changes when the position or the split of a copy changes for the same input. */
#define GSHGEOM_LATTICE_REVISION 1

/** A lattice: the cell of copy `(0, 0)` and the steps. */
typedef struct GSHGEOM_Lattice {
  double x, y, w, h;              /**< The cell `[x, x + w] x [y, y + h]`; finite, `w` and `h` greater than zero. */
  double x_step, y_step;          /**< The steps; finite, not negative, and greater than zero on an axis that repeats. */
  uint32_t repeat_x, repeat_y;    /**< 1 repeats along the axis; 0 keeps index 0 on it.  Other values are invalid. */
} GSHGEOM_Lattice;

/** Receives one row of copies: the copies `(i_lo .. i_hi, j)`, `i_lo <= i_hi`.  A nonzero result stops the enumeration. */
typedef int (*GSHGEOM_Lattice_Row_Fn)(void * ctx, int64_t j, int64_t i_lo, int64_t i_hi);

/** The Hermite normal form of an integer lattice (LT3). */
typedef struct GSHGEOM_Lattice_Period {
  int64_t p, r, s;
} GSHGEOM_Lattice_Period;

/**
 * @brief The cell of copy `(i, j)` (LT2): `out = {left, top, right, bottom}`.
 *
 * @return ::GSHGEOM_OK; ::GSHGEOM_ERR_INVALID for NULL, an index past ::GSHGEOM_LATTICE_INDEX_MAX or a lattice that is not
 *   valid (a non-finite value, a size that is not positive, a step that is negative, or a `repeat` other than 0 and 1).
 *   @p out is written only on success.  An axis that does not repeat is not restricted here.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_lattice_cell(const GSHGEOM_Lattice * lattice, int64_t i, int64_t j, double out[4]);

/**
 * @brief The device translation of copy `(i, j)` (LT2): `((tau + i a) + j b)` per component.
 *
 * @param map A valid 2D map from the lattice's space to the device; `tau` is its translation.
 * @return ::GSHGEOM_OK; ::GSHGEOM_ERR_INVALID for NULL, an invalid map or lattice, an index past the range, or a result
 *   that is not finite.  @p out is written only on success.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_lattice_translation(const GSHGEOM_Map * map, const GSHGEOM_Lattice * lattice, int64_t i,
    int64_t j, double out[2]);

/**
 * @brief The quadrilateral of copy `(0, 0)` under the linear part of @p map alone: `out = x0, y0, ..., x3, y3` for the
 *   corners `(left, top)`, `(right, top)`, `(right, bottom)`, `(left, bottom)`, each coordinate one expression
 *   `m00 px + m01 py` (or `m10 px + m11 py`) with the product and the sum rounded and no fused multiply-add.
 *
 * The device corners of copy `(i, j)` are these plus ::gshgeom_lattice_translation, one addition each.
 *
 * @return As ::gshgeom_lattice_translation.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_lattice_quad(const GSHGEOM_Map * map, const GSHGEOM_Lattice * lattice, double out[8]);

/**
 * @brief The copies of a use (LT1), in lattice order, as rows.
 *
 * @param map The map from the lattice's space to the device: valid, 2D, and with a linear part of nonzero determinant
 *   (otherwise ::GSHGEOM_ERR_INVALID; a caller handles a map that draws nothing first).
 * @param lattice The lattice.  With both axes repeating, the device step vectors must not be parallel (otherwise
 *   ::GSHGEOM_ERR_INVALID).
 * @param rect `{x0, y0, x1, y1}`: the closed rectangle, whole numbers of magnitude at most `2^53`, `x0 <= x1`, `y0 <= y1`.
 * @param max_count The most copies to report, at most `2^53 - 1`.
 * @param row_fn Called once for each row that has copies, in increasing `j`; may be NULL (count only).  A nonzero result
 *   stops the enumeration, and the call returns ::GSHGEOM_OK with @p count the copies reported so far, that row included.
 * @param ctx Passed to @p row_fn.
 * @param count Receives the number of copies reported (and, for ::GSHGEOM_ERR_LIMIT, `max_count + 1`).  Written whenever
 *   the arguments are valid.
 * @return ::GSHGEOM_OK; ::GSHGEOM_ERR_LIMIT when the copies of the use are more than @p max_count (or the rows to visit are
 *   more than LT1's bound): the rows reported before are whole rows and the last is cut so that exactly @p max_count
 *   copies were reported (when the row that crosses the bound would be cut to nothing, it is not reported);
 *   ::GSHGEOM_ERR_INVALID.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_lattice_each(const GSHGEOM_Map * map, const GSHGEOM_Lattice * lattice,
    const int64_t rect[4], uint64_t max_count, GSHGEOM_Lattice_Row_Fn row_fn, void * ctx, uint64_t * count);

/**
 * @brief The device step vectors of LT2 (@p a for the first step, @p b for the second), as the doubles they are.
 *
 * @return As ::gshgeom_lattice_translation.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_lattice_steps(const GSHGEOM_Map * map, const GSHGEOM_Lattice * lattice, double a[2], double b[2]);

/**
 * @brief Whether the lattice has an integer period (LT3, LT5): both axes repeat, every component of the device step vectors is an
 *   exact integer of magnitude below `2^31`, and their determinant is not zero and has magnitude at most
 *   ::GSHGEOM_LATTICE_PERIOD_MAX.
 *
 * @param a,b Receive the vectors when @p integer is 1 (and are not written otherwise).
 * @param integer Receives 1 or 0.
 * @return ::GSHGEOM_OK; ::GSHGEOM_ERR_INVALID for NULL, an invalid map or lattice.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_lattice_integer_steps(const GSHGEOM_Map * map, const GSHGEOM_Lattice * lattice, int64_t a[2],
    int64_t b[2], int * integer);

/**
 * @brief LT5: the canonical map and the whole-pixel shift of the copy `(i, j)`.
 *
 * @p canonical is @p map with its translation replaced by the copy's fraction in `[0, 1)` and `(shift_x, shift_y)` are the whole
 * pixels, so that the copy is the cell under @p canonical moved by the shift.
 *
 * @return ::GSHGEOM_OK; ::GSHGEOM_ERR_INVALID for NULL, an invalid map or lattice, an index past the range or a position that is
 *   not finite; ::GSHGEOM_ERR_LIMIT for a whole-pixel coordinate of magnitude `2^40` or more.  Nothing is written on failure.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_lattice_canonical(const GSHGEOM_Map * map, const GSHGEOM_Lattice * lattice, int64_t i, int64_t j,
    GSHGEOM_Map * canonical, int64_t * shift_x, int64_t * shift_y);

/**
 * @brief The Hermite normal form of the lattice `Z a + Z b` (LT3).
 *
 * @param a,b Integer vectors with components of magnitude less than `2^31` and a determinant `D = a0 b1 - a1 b0` of
 *   magnitude at most ::GSHGEOM_LATTICE_PERIOD_MAX.
 * @return ::GSHGEOM_OK; ::GSHGEOM_ERR_INVALID for NULL, a component out of range or `D = 0`; ::GSHGEOM_ERR_LIMIT for `|D|`
 *   past the bound.  @p out is written only on success.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_lattice_period(const int64_t a[2], const int64_t b[2], GSHGEOM_Lattice_Period * out);

/**
 * @brief `rep(x, y)` of LT3: the point of the window `[0, p) x [0, s)` congruent to `(x, y)`.  @p x and @p y have magnitude at most
 *   `2^62`; @p period must be one ::gshgeom_lattice_period returned.
 */
GSHGEOM_API void gshgeom_lattice_reduce(const GSHGEOM_Lattice_Period * period, int64_t x, int64_t y, int64_t out[2]);

#ifdef __cplusplus
}
#endif

#endif /* GHOTI_IO_GSHGEOM_LATTICE_H */
