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
 * @file curve.h
 *
 * Curves: lines and quadratic, cubic and rational curves in 2D or 3D, with
 * evaluation and bounds.
 *
 * A curve is a Bezier curve of degree 1, 2 or 3 given by its control points
 * (the first and last are the ends), or a rational one of degree 2 or 3,
 * which also has a positive weight at each control point and is the curve
 * `sum(w_i P_i B_i(t)) / sum(w_i B_i(t))`.  A rational quadratic with the
 * right weights is exactly a conic section, a circle's arc among them.
 *
 * Model space is `double`, y down (AD-6).  The dimension (2 or 3) is a
 * property of the curve; in 2D the `z` of a point is ignored on input and is
 * zero on output.
 *
 * ### Validity
 *
 * Every function that takes a ::GSHGEOM_Curve checks it, so a structure filled
 * in by hand is as safe as one a constructor made.  A curve is valid when its
 * kind and dimension are in range, every coordinate it uses is finite and at
 * most ::GSHGEOM_COORD_LIMIT in magnitude, and, for a rational curve, every
 * weight is finite, positive and within [::GSHGEOM_WEIGHT_MIN,
 * ::GSHGEOM_WEIGHT_MAX].  A NaN or infinity anywhere, a weight that is zero or
 * negative, a kind or dimension out of range: ::GSHGEOM_ERR_INVALID.  A finite
 * number outside its range: ::GSHGEOM_ERR_LIMIT (never clamped).  The range
 * is what guarantees that every intermediate of an evaluation or a bound is
 * finite, so no function here returns a non-finite number.
 */

#ifndef GHOTI_IO_GSHGEOM_CURVE_H
#define GHOTI_IO_GSHGEOM_CURVE_H

#include <ghoti.io/shape-geometry/core.h>
#include <ghoti.io/shape-geometry/macros.h>

#ifdef __cplusplus
extern "C" {
#endif

/** The largest coordinate magnitude accepted: 2^100. */
#define GSHGEOM_COORD_LIMIT 0x1p100

/** The smallest rational weight accepted: 2^-100. */
#define GSHGEOM_WEIGHT_MIN 0x1p-100

/** The largest rational weight accepted: 2^100. */
#define GSHGEOM_WEIGHT_MAX 0x1p100

/**
 * @brief The kind of a curve.  The values are pinned; any further kind is
 * appended before ::GSHGEOM_CURVE_KIND_COUNT.
 */
typedef enum GSHGEOM_Curve_Kind {
  GSHGEOM_CURVE_LINE = 0,               ///< Two control points.
  GSHGEOM_CURVE_QUADRATIC = 1,          ///< Three control points.
  GSHGEOM_CURVE_CUBIC = 2,              ///< Four control points.
  GSHGEOM_CURVE_RATIONAL_QUADRATIC = 3, ///< Three control points and weights.
  GSHGEOM_CURVE_RATIONAL_CUBIC = 4,     ///< Four control points and weights.
  GSHGEOM_CURVE_KIND_COUNT              ///< Closes the enum.
} GSHGEOM_Curve_Kind;

/** @brief A point: `z` is ignored in 2D and zero on output. */
typedef struct GSHGEOM_Point {
  double x; ///< Across, positive to the right.
  double y; ///< Down, positive toward the bottom (AD-6).
  double z; ///< Depth; zero in 2D.
} GSHGEOM_Point;

/**
 * @brief A curve: a value, copied freely.
 *
 * The constructors set every unused slot to zero.  Only the first
 * ::gshgeom_curve_point_count() points are used, and `w` only for a rational
 * kind.
 */
typedef struct GSHGEOM_Curve {
  GSHGEOM_Curve_Kind kind; ///< The kind.
  int dim;                 ///< 2 or 3.
  GSHGEOM_Point p[4];      ///< Control points; the first and last are the ends.
  double w[4];             ///< Weights of the control points (rational only).
} GSHGEOM_Curve;

/**
 * @brief An axis-aligned box.
 *
 * `empty` is set when there is nothing to bound (a path with no segment), and
 * then `min` and `max` are zero.  In 2D the `z` of both is zero.
 */
typedef struct GSHGEOM_Bounds {
  int empty;          ///< Nonzero if there is nothing in the box.
  GSHGEOM_Point min;  ///< The smallest coordinate on each axis.
  GSHGEOM_Point max;  ///< The largest coordinate on each axis.
} GSHGEOM_Bounds;

/**
 * @brief How many control points a kind has.
 *
 * @param kind The kind.
 * @return 2, 3 or 4; 0 for a value that is not a kind.
 */
GSHGEOM_API int gshgeom_curve_point_count(GSHGEOM_Curve_Kind kind);

/**
 * @brief Whether a kind carries weights.
 *
 * @param kind The kind.
 * @return 1 for the rational kinds, 0 otherwise (and for a non-kind).
 */
GSHGEOM_API int gshgeom_curve_is_rational(GSHGEOM_Curve_Kind kind);

/**
 * @brief Make a line.
 *
 * @param dim 2 or 3.
 * @param a Start.
 * @param b End.
 * @param out The curve; written only on success.
 * @return ::GSHGEOM_OK, ::GSHGEOM_ERR_INVALID or ::GSHGEOM_ERR_LIMIT.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_curve_line(
    int dim, GSHGEOM_Point a, GSHGEOM_Point b, GSHGEOM_Curve * out);

/**
 * @brief Make a quadratic Bezier curve.
 *
 * @param dim 2 or 3.
 * @param a Start.
 * @param c Control point.
 * @param b End.
 * @param out The curve; written only on success.
 * @return ::GSHGEOM_OK, ::GSHGEOM_ERR_INVALID or ::GSHGEOM_ERR_LIMIT.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_curve_quadratic(int dim, GSHGEOM_Point a,
    GSHGEOM_Point c, GSHGEOM_Point b, GSHGEOM_Curve * out);

/**
 * @brief Make a cubic Bezier curve.
 *
 * @param dim 2 or 3.
 * @param a Start.
 * @param c1 First control point.
 * @param c2 Second control point.
 * @param b End.
 * @param out The curve; written only on success.
 * @return ::GSHGEOM_OK, ::GSHGEOM_ERR_INVALID or ::GSHGEOM_ERR_LIMIT.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_curve_cubic(int dim, GSHGEOM_Point a,
    GSHGEOM_Point c1, GSHGEOM_Point c2, GSHGEOM_Point b, GSHGEOM_Curve * out);

/**
 * @brief Make a rational quadratic Bezier curve.
 *
 * @param dim 2 or 3.
 * @param a Start.
 * @param c Control point.
 * @param b End.
 * @param w Three weights, for `a`, `c` and `b`; positive and in range.
 * @param out The curve; written only on success.
 * @return ::GSHGEOM_OK, ::GSHGEOM_ERR_INVALID or ::GSHGEOM_ERR_LIMIT.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_curve_rational_quadratic(int dim,
    GSHGEOM_Point a, GSHGEOM_Point c, GSHGEOM_Point b, const double w[3],
    GSHGEOM_Curve * out);

/**
 * @brief Make a rational cubic Bezier curve.
 *
 * @param dim 2 or 3.
 * @param a Start.
 * @param c1 First control point.
 * @param c2 Second control point.
 * @param b End.
 * @param w Four weights, in the order of the points; positive and in range.
 * @param out The curve; written only on success.
 * @return ::GSHGEOM_OK, ::GSHGEOM_ERR_INVALID or ::GSHGEOM_ERR_LIMIT.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_curve_rational_cubic(int dim,
    GSHGEOM_Point a, GSHGEOM_Point c1, GSHGEOM_Point c2, GSHGEOM_Point b,
    const double w[4], GSHGEOM_Curve * out);

/**
 * @brief Check a curve against the validity rules in the file comment.
 *
 * @param curve The curve.
 * @return ::GSHGEOM_OK, ::GSHGEOM_ERR_INVALID (NULL, kind or dimension out of
 *   range, a non-finite number, a weight that is not positive) or
 *   ::GSHGEOM_ERR_LIMIT (a finite number out of range).
 */
GSHGEOM_API GSHGEOM_Result gshgeom_curve_check(const GSHGEOM_Curve * curve);

/**
 * @brief The point of a curve at a parameter.
 *
 * By de Casteljau's algorithm (in homogeneous coordinates for a rational
 * curve).  `t = 0` and `t = 1` give the end control points exactly.  The
 * result is a function of the arguments alone, bit for bit, on every target.
 *
 * @param curve The curve.
 * @param t The parameter, in [0, 1].
 * @param out The point; written only on success.  `z` is zero in 2D.
 * @return ::GSHGEOM_OK, or ::GSHGEOM_ERR_INVALID (an invalid curve, or `t`
 *   NaN, infinite or outside [0, 1]), or ::GSHGEOM_ERR_LIMIT.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_curve_eval(
    const GSHGEOM_Curve * curve, double t, GSHGEOM_Point * out);

/**
 * @brief The box of the curve itself, not of its control polygon.
 *
 * The box of the curve evaluated at both ends and at every parameter in (0, 1)
 * where a coordinate has an extremum.  Those parameters are the roots of
 * `N'D - ND'` for that coordinate (`N/D` the coordinate as a rational function
 * of `t`; `D = 1` for a polynomial curve), found by recursion on the
 * derivative and bisection, to about 2^-60 in `t`.  The coordinate's value at
 * an extremum is second order in that error, so the box is the curve's own to
 * the rounding of the evaluation: the end points are inside it exactly, and a
 * point of the curve is inside it to a few units in the last place.  It is not
 * a conservative box.
 *
 * @param curve The curve.
 * @param out The box (never empty); written only on success.
 * @return ::GSHGEOM_OK, ::GSHGEOM_ERR_INVALID or ::GSHGEOM_ERR_LIMIT.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_curve_bounds(
    const GSHGEOM_Curve * curve, GSHGEOM_Bounds * out);

#ifdef __cplusplus
}
#endif

#endif /* GHOTI_IO_GSHGEOM_CURVE_H */
