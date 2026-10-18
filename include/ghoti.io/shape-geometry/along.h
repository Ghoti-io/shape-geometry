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
 * @file along.h
 *
 * Along-path parametrization (CAP-5, CAP-13; AD-1, AD-6): the arc length of a
 * model path and the point, tangent and normal at a distance along it.
 *
 * A *measured path* is built once from a model path and queried many times.
 * It is opaque and immutable; it keeps a copy of every segment that has
 * extent and a table of the quadrature's panels, so it does not refer to the
 * path it was made from.  A measured path is used from one thread at a time
 * for building and destroying, and may be queried from several (a query
 * writes nothing).
 *
 * ### Distance
 *
 * The distance along a path is the sum of the arc lengths of its segments, in
 * element order.  A move contributes nothing: where a move separates two
 * segments the point jumps and the distance does not, so a distance runs across
 * the subpaths in order.  A close is the straight segment back to the start (of
 * length zero if the path is there already).  A segment of length exactly zero
 * has no extent and no distance selects it.  The valid distances are `[0, L]`,
 * `L` being ::gshgeom_measured_length.  A distance that is exactly a boundary
 * between two segments (or two subpaths) belongs to the later one; `L` is the
 * end of the last segment with extent.
 *
 * ### Arc length
 *
 * A line's length is the hypotenuse of its difference.  A curve's is the
 * integral of its speed `|B'(t)|` over `[0, 1]` (for a rational curve
 * `B' = (N'D - ND') / D^2`, `N` and `D` the homogeneous numerator and
 * denominator), by adaptive eight-point Gauss-Legendre: a panel `[a, b]` is
 * accepted when its value differs from the sum of its two halves by at most
 * `tol * S * (b - a)`, `S` being the length of the segment's control polygon
 * (its projected control points, for a rational curve); the halves are then the
 * panels of the table.  At most ::GSHGEOM_MEASURE_MAX_DEPTH levels are used; a
 * panel at the last level is accepted whatever it estimates, which only a cusp
 * or an extreme weight ratio reaches, and then to a negligible error.
 *
 * The estimate is the difference of two evaluations of a smooth integrand, so it
 * bounds the error of a segment in practice and not by proof.  The contract is:
 * *the estimated error of a segment is at most `tol * S`*.  The tests measure
 * the real error against an independent high-precision reference.  The running
 * sum of the table adds at most `rows * 2^-53` of itself.
 *
 * Everything is `+ - * /`, ::gmath_hypot and a fixed order of operations, so
 * the table, the queries and the hash are the same bits on every target.
 *
 * ### Tangent and normal
 *
 * The tangent is the unit vector `B'(u) / |B'(u)|`.  Where `B'` vanishes (a
 * cusp, or a control point repeated at an end) it is the direction of the first
 * nonzero derivative of `E = N'D - ND'`, which has the direction of `B'`, taken
 * on the later side: the right-hand limit and, at `u = 1`, the left-hand limit,
 * which is the direction the curve arrives in.  "Vanishes" means `|E|` is at most
 * `2^-48` of the largest term of `E` in that segment.  If none of the first five
 * derivatives is nonzero the segment has no direction and `has_tangent` is 0.
 *
 * In 2D the normal is ::gshgeom_normal_left of the tangent (orient.h; AD-6): the
 * tangent turned by a negative angle, `(t.y, -t.x)`, so for a tangent along +x
 * it is `(0, -1)`, up on the screen, and on a positively wound circle it points
 * out.  It is unit and exactly orthogonal to the tangent.  The right normal is
 * its negation (::gshgeom_normal_right), applied by the caller.  There is no
 * normal in 3D, where "left" has no meaning: `has_normal` is 0 and it is zero.
 *
 * ### The hash of a measured path (version 1)
 *
 * In the spelling of hash.h: `u32` 1 (this format), `u32` kind (`MEAS`, so that
 * a measured path never shares a stream with another intermediate), `u32`
 * dimension, `u64` segment count, `u64` row count, `u64` subpath count, the first point of the
 * path (`f64` `x`, `y`, and `z` in 3D; zeros for a path with no element), the
 * total length `f64`; then for each segment the `u64` element index, the `u8`
 * curve kind, the control points `f64` and, for a rational kind, the weights
 * `f64`; for each row the `u64` segment index and the `f64` values `t0`, `t1`,
 * `s0` and `len`; for each subpath `f64` begin, `f64` end and `u8` closed.  The
 * tolerance, the limits and the allocator are not part of it.
 */

#ifndef GHOTI_IO_GSHGEOM_ALONG_H
#define GHOTI_IO_GSHGEOM_ALONG_H

#include <ghoti.io/shape-geometry/allocator.h>
#include <ghoti.io/shape-geometry/core.h>
#include <ghoti.io/shape-geometry/curve.h>
#include <ghoti.io/shape-geometry/macros.h>
#include <ghoti.io/shape-geometry/orient.h>
#include <ghoti.io/shape-geometry/path.h>

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/** The smallest tolerance accepted: 2^-46 (about 1.4e-14). */
#define GSHGEOM_MEASURE_TOLERANCE_MIN 0x1p-46

/** The largest tolerance accepted: 2^-4. */
#define GSHGEOM_MEASURE_TOLERANCE_MAX 0x1p-4

/** The tolerance a request of 0 means: 2^-40 (about 9.1e-13). */
#define GSHGEOM_MEASURE_TOLERANCE_DEFAULT 0x1p-40

/** The deepest level of subdivision of one segment. */
#define GSHGEOM_MEASURE_MAX_DEPTH 52

/** @brief A measured path.  Opaque. */
typedef struct GSHGEOM_Measured GSHGEOM_Measured;

/**
 * @brief The answer to a query by distance.
 *
 * `tangent` and `normal` are unit vectors when `has_tangent` and `has_normal`
 * are set and zero otherwise.  `z` is zero in 2D.
 */
typedef struct GSHGEOM_Measured_Frame {
  GSHGEOM_Point point;   ///< The point on the path.
  GSHGEOM_Point tangent; ///< The unit tangent, in the direction of travel.
  GSHGEOM_Point normal;  ///< The unit left normal (2D only).
  size_t element;        ///< The index of the path element the point is on (0 for a path with no extent).
  double t;              ///< The parameter of the point on that element, in [0, 1].
  int has_tangent;       ///< Nonzero if the point has a tangent.
  int has_normal;        ///< Nonzero in 2D when it has a tangent.
} GSHGEOM_Measured_Frame;

/** @brief The extent of one subpath in distance. */
typedef struct GSHGEOM_Measured_Subpath {
  double begin; ///< The distance at which the subpath starts.
  double end;   ///< The distance at which it ends (equal to `begin` for zero length).
  int closed;   ///< Nonzero if the subpath ended in a close.
} GSHGEOM_Measured_Subpath;

/**
 * @brief Measure a path.
 *
 * @param alloc The allocator, or NULL for the default; copied.
 * @param limits The limits, or NULL for the defaults; `max_points` bounds the
 *   number of rows of the table, exactly (a path whose rows number `max_points`
 *   is accepted, one more is ::GSHGEOM_ERR_LIMIT).
 * @param path The path, 2D or 3D, any number of subpaths.
 * @param tolerance The relative tolerance of the arc lengths, or 0 for
 *   ::GSHGEOM_MEASURE_TOLERANCE_DEFAULT.
 * @param out The measured path; written only on success.  Free it with
 *   ::gshgeom_measured_destroy.
 * @return ::GSHGEOM_OK; ::GSHGEOM_ERR_INVALID (NULL, a tolerance that is
 *   negative, NaN or infinite); ::GSHGEOM_ERR_LIMIT (a tolerance outside
 *   [::GSHGEOM_MEASURE_TOLERANCE_MIN, ::GSHGEOM_MEASURE_TOLERANCE_MAX], or too
 *   many rows); ::GSHGEOM_ERR_OOM.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_measured_create(
    const GSHGEOM_Allocator * alloc, const GSHGEOM_Limits * limits,
    const GSHGEOM_Path * path, double tolerance, GSHGEOM_Measured ** out);

/** @brief Free a measured path.  NULL is ignored. */
GSHGEOM_API void gshgeom_measured_destroy(GSHGEOM_Measured * measured);

/** @brief The dimension of the path, 2 or 3; 0 for NULL. */
GSHGEOM_API int gshgeom_measured_dimension(const GSHGEOM_Measured * measured);

/** @brief The total length; 0 for NULL and for a path with no extent. */
GSHGEOM_API double gshgeom_measured_length(const GSHGEOM_Measured * measured);

/** @brief How many segments with extent the path holds; 0 for NULL. */
GSHGEOM_API size_t gshgeom_measured_segment_count(const GSHGEOM_Measured * measured);

/** @brief How many rows (quadrature panels) the table holds; 0 for NULL. */
GSHGEOM_API size_t gshgeom_measured_row_count(const GSHGEOM_Measured * measured);

/** @brief How many subpaths the path has (see along.h's description); 0 for NULL. */
GSHGEOM_API size_t gshgeom_measured_subpath_count(const GSHGEOM_Measured * measured);

/**
 * @brief The extent of subpath @p index in distance.
 *
 * @return ::GSHGEOM_OK or ::GSHGEOM_ERR_INVALID (NULL, index out of range).
 */
GSHGEOM_API GSHGEOM_Result gshgeom_measured_subpath(const GSHGEOM_Measured * measured,
    size_t index, GSHGEOM_Measured_Subpath * out);

/**
 * @brief The point, tangent and normal at distance @p s.
 *
 * @param measured The measured path.
 * @param s The distance, in [0, length].
 * @param out The answer; written only on success.
 * @return ::GSHGEOM_OK; ::GSHGEOM_ERR_INVALID for NULL, NaN, infinity, a
 *   distance outside [0, length], or a path with no element (it has no point).
 */
GSHGEOM_API GSHGEOM_Result gshgeom_measured_at(const GSHGEOM_Measured * measured,
    double s, GSHGEOM_Measured_Frame * out);

/**
 * @brief As ::gshgeom_measured_at, with a finite distance clamped into
 * [0, length].
 *
 * @return As ::gshgeom_measured_at, except that only NaN and infinity (and
 *   NULL and an empty path) are refused.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_measured_at_clamped(
    const GSHGEOM_Measured * measured, double s, GSHGEOM_Measured_Frame * out);

/**
 * @brief The arc length of one curve, by the rule of this file.
 *
 * Allocates nothing.
 *
 * @param curve The curve (valid as ::gshgeom_curve_check requires).
 * @param tolerance As for ::gshgeom_measured_create.
 * @param out The length; written only on success.
 * @return ::GSHGEOM_OK; ::GSHGEOM_ERR_INVALID; ::GSHGEOM_ERR_LIMIT (a
 *   tolerance out of range, an out-of-range number in the curve, or more rows
 *   than ::GSHGEOM_LIMITS_DEFAULT_MAX_POINTS).
 */
GSHGEOM_API GSHGEOM_Result gshgeom_measured_curve_length(
    const GSHGEOM_Curve * curve, double tolerance, double * out);

/**
 * @brief The hash of the measured path, in the spelling in this file's comment.
 *
 * @param measured The measured path, or NULL.
 * @return The hash; for NULL, the offset basis.
 */
GSHGEOM_API uint64_t gshgeom_measured_hash(const GSHGEOM_Measured * measured);

/**
 * @brief Write the measured path as text and end with its hash.
 *
 * The first line is `gshgeom.measured 1 dim=D segments=N rows=P subpaths=S
 * length=L hash=H`, then one line per segment (`G element kind` and the
 * control points, then the weights of a rational kind), per row (`R segment t0 t1
 * s0 len`) and per subpath (`S begin end closed`).  Integers are decimal and
 * every `double` is the 16 hexadecimal digits of its bit pattern.
 *
 * @return ::GSHGEOM_OK, ::GSHGEOM_ERR_INVALID (NULL) or ::GSHGEOM_ERR_IO.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_measured_dump(
    const GSHGEOM_Measured * measured, FILE * out);

#ifdef __cplusplus
}
#endif

#endif /* GHOTI_IO_GSHGEOM_ALONG_H */
