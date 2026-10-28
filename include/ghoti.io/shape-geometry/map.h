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
 * @file map.h
 *
 * Maps: the affine transforms that place model paths, composed in one place
 * only and carrying a kind tag that survives composition (AD-24, CAP-2,
 * CAP-4).
 *
 * A map sends a point `p` to `M p + t`, `M` a row-major 3 by 3 and `t` a
 * vector, in 2D or 3D.  A 2D map uses the upper left 2 by 2 and `t[0..1]`;
 * every other entry is ignored on input and is the identity (zero off the
 * diagonal) on output, and a 2D map takes a 2D point's `x` and `y` and gives
 * `z = 0`.  Every entry is finite with magnitude at most ::GSHGEOM_MAP_LIMIT,
 * so that applying a map to a point of story 5's range never makes a
 * non-finite number; a finite number outside the bound is
 * ::GSHGEOM_ERR_LIMIT and is never clamped.
 *
 * **One composer.**  ::gshgeom_map_compose is the only place that multiplies
 * two maps, and ::gshgeom_map_compose_chain, which a scene uses, calls it.
 * The order is the spine's: root to leaf, then the target map last.  With
 * `chain[0]` the root and `chain[n-1]` the leaf, the result sends a leaf-local
 * point `p` to `target(root(...(leaf(p))))`.  Each entry of a product is
 * formed as `((a0 b0 + a1 b1) + a2 b2)` with no fused multiply-add, so the
 * composed map has the same bits on every target and whoever composed it.
 *
 * **The kind tag** is a set of the transform classes a map was built from.
 * A constructor sets it from the matrix it builds (exactly: the sign of a
 * determinant is decided by ::gshgeom_orient3d's exact arithmetic); a
 * composition is the union of its factors' tags, so a tag survives
 * composition and never loses a class.  A tag may claim more than the net
 * matrix has: a quarter turn followed by its inverse is the identity and is
 * still tagged as a quarter turn, which sends a glyph to the general route.
 * The predicates read the tag (and, so that a hand-filled map cannot lie to
 * them, the values as well).
 *
 * Rotations take the source's angle unit.  In degrees, gradians and turns a
 * quarter turn is exact; the angle in radians is not snapped (the double
 * nearest pi/2 is not a quarter turn): ask in a unit that has it.
 *
 * Copyright 2026 by Corey Pennycuff
 */

#ifndef GHOTI_IO_GSHGEOM_MAP_H
#define GHOTI_IO_GSHGEOM_MAP_H

#include <ghoti.io/shape-geometry/core.h>
#include <ghoti.io/shape-geometry/curve.h>
#include <ghoti.io/shape-geometry/macros.h>

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/** The largest magnitude of a map entry: `2^100`, story 5's coordinate limit. */
#define GSHGEOM_MAP_LIMIT 0x1p100

/** The most quarter turns a rotation in degrees, gradians or turns takes (`2^40`). */
#define GSHGEOM_MAP_MAX_QUARTERS 0x1p40

/** The unit of an angle: what the source states, never converted by the caller. */
typedef enum GSHGEOM_Angle_Unit {
  GSHGEOM_ANGLE_RADIANS = 0, ///< Radians.
  GSHGEOM_ANGLE_DEGREES = 1, ///< Degrees: a quarter turn is 90.
  GSHGEOM_ANGLE_GRADIANS = 2, ///< Gradians (gon): a quarter turn is 100.
  GSHGEOM_ANGLE_TURNS = 3,   ///< Whole turns: a quarter turn is 0.25.
  GSHGEOM_ANGLE_UNIT_COUNT   ///< Closes the enum.
} GSHGEOM_Angle_Unit;

/** The axis of a 3D rotation.  A positive angle turns the next axis toward the one after, cyclically: x to y about z, y to z about x, z to x about y. */
typedef enum GSHGEOM_Axis {
  GSHGEOM_AXIS_X = 0, ///< About x: y turns toward z.
  GSHGEOM_AXIS_Y = 1, ///< About y: z turns toward x.
  GSHGEOM_AXIS_Z = 2, ///< About z: x turns toward y (the 2D rotation).
  GSHGEOM_AXIS_COUNT  ///< Closes the enum.
} GSHGEOM_Axis;

/**
 * The classes a map was built from, as a bit set (the kind tag).  Zero is the
 * identity.
 */
typedef enum GSHGEOM_Map_Kind {
  GSHGEOM_MAP_TRANSLATES = 1u << 0,        ///< A translation other than zero.
  GSHGEOM_MAP_SCALES = 1u << 1,            ///< A uniform scale other than 1.
  GSHGEOM_MAP_SCALES_NONUNIFORM = 1u << 2, ///< Scales that differ between axes.
  GSHGEOM_MAP_FLIPS_Y = 1u << 3,           ///< The canonical y flip (design y-up to model y-down).
  GSHGEOM_MAP_MIRRORS = 1u << 4,           ///< Any other mirror.
  GSHGEOM_MAP_QUARTER_TURNS = 1u << 5,     ///< An exact quarter or half turn, or an axis swap.
  GSHGEOM_MAP_ROTATES = 1u << 6,           ///< A rotation by any other angle.
  GSHGEOM_MAP_SKEWS = 1u << 7,             ///< A skew, or any other non-similar linear part.
  GSHGEOM_MAP_SINGULAR = 1u << 8,          ///< Not invertible (the determinant is exactly zero).
  GSHGEOM_MAP_GENERAL = 1u << 9,           ///< A 3D linear part of no known structure.
  GSHGEOM_MAP_KIND_ALL = (1u << 10) - 1u   ///< Every defined bit.
} GSHGEOM_Map_Kind;

/** A map.  Plain data, validated by every function that takes one. */
typedef struct GSHGEOM_Map {
  uint32_t kind;  ///< The kind tag: a union of ::GSHGEOM_Map_Kind bits.
  int dim;        ///< 2 or 3.
  double m[3][3]; ///< The linear part, row major: `out[i] = sum_j m[i][j] in[j] + t[i]`.
  double t[3];    ///< The translation.
} GSHGEOM_Map;

/** What the route predicate returns for an accepted map. */
typedef struct GSHGEOM_Grid {
  double scale; ///< The uniform scale `s > 0`.
  double tx;    ///< The translation in x: `x' = s x + tx`.
  double ty;    ///< The translation in y: `y' = (flip_y ? -s : s) y + ty`.
  int flip_y;   ///< Nonzero if the map carries the canonical y flip.
} GSHGEOM_Grid;

/** The identity of dimension `dim` (2 or 3). */
GSHGEOM_API GSHGEOM_Result gshgeom_map_identity(int dim, GSHGEOM_Map * out);

/** A translation.  `tz` is ignored in 2D. */
GSHGEOM_API GSHGEOM_Result gshgeom_map_translate(
    int dim, double tx, double ty, double tz, GSHGEOM_Map * out);

/**
 * A scale about the origin.  `sz` is ignored in 2D.  A zero factor is
 * allowed and makes a singular map; a negative one is a mirror, or a half turn
 * if two are negative.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_map_scale(
    int dim, double sx, double sy, double sz, GSHGEOM_Map * out);

/** The canonical y flip, `y' = -y`: the one a y-up source takes into model space (AD-6). */
GSHGEOM_API GSHGEOM_Result gshgeom_map_flip_y(int dim, GSHGEOM_Map * out);

/**
 * The rotation in the plane (a 2D map), positive from +x toward +y.
 *
 * In degrees, gradians and turns the angle is reduced by an exact subtraction
 * to within half a quarter turn, so the quarter turns are exact (entries 0 and
 * +-1, never -0) and only the remainder goes through `gmath_sin` and
 * `gmath_cos`.  In radians the angle goes through them whole.
 *
 * @return ::GSHGEOM_ERR_INVALID for a NaN or an infinite angle or an unknown
 *         unit; ::GSHGEOM_ERR_LIMIT for more than ::GSHGEOM_MAP_MAX_QUARTERS
 *         quarter turns in a unit other than radians.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_map_rotate(
    GSHGEOM_Angle_Unit unit, double angle, GSHGEOM_Map * out);

/** The rotation about an axis (a 3D map), with the sense ::GSHGEOM_Axis states. */
GSHGEOM_API GSHGEOM_Result gshgeom_map_rotate_axis(
    GSHGEOM_Axis axis, GSHGEOM_Angle_Unit unit, double angle, GSHGEOM_Map * out);

/**
 * A map from a matrix and a translation, with its tag classified from the
 * values.
 *
 * @param dim 2 or 3.
 * @param m   `dim * dim` entries, row major.
 * @param t   `dim` entries.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_map_from_matrix(
    int dim, const double * m, const double * t, GSHGEOM_Map * out);

/**
 * The map a y-up design space takes into model space at one size: `x' = s x +
 * tx`, `y' = -s y + ty` (a 2D map, the canonical flip composed with the
 * scale and the translation).  It is the map ::gshgeom_map_grid_aligned
 * accepts.
 *
 * @param scale Finite and greater than zero.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_map_design_to_model(
    double scale, double tx, double ty, GSHGEOM_Map * out);

/** Whether the map is valid: dimension, kind bits, every entry finite and at most ::GSHGEOM_MAP_LIMIT. */
GSHGEOM_API GSHGEOM_Result gshgeom_map_check(const GSHGEOM_Map * map);

/**
 * The composition `outer . inner`: the map that applies `inner` first and
 * `outer` second.  The one composer.  `out` may be either input.
 *
 * @return ::GSHGEOM_ERR_INVALID for an invalid map or maps of different
 *         dimensions; ::GSHGEOM_ERR_LIMIT if an entry of the result passes the
 *         bound.  `out` is untouched on an error.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_map_compose(
    const GSHGEOM_Map * outer, const GSHGEOM_Map * inner, GSHGEOM_Map * out);

/**
 * The composed map of a scene path: `chain[0]` is the root and `chain[count-1]`
 * the leaf; the result is `target . chain[0] . ... . chain[count-1]`, formed by
 * composing the chain from the root down and the target last.  A NULL `target`
 * is the identity; `count` may be 0.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_map_compose_chain(const GSHGEOM_Map * chain,
    size_t count, const GSHGEOM_Map * target, GSHGEOM_Map * out);

/**
 * Apply the map to a point.  The point's dimension is the map's: a 2D map
 * ignores `z` and gives `z = 0`.
 *
 * @return ::GSHGEOM_ERR_INVALID for a NULL, an invalid map or a non-finite
 *         point; ::GSHGEOM_ERR_LIMIT if the result is not finite.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_map_apply_point(
    const GSHGEOM_Map * map, GSHGEOM_Point p, GSHGEOM_Point * out);

/** Apply only the linear part (a vector is not translated). */
GSHGEOM_API GSHGEOM_Result gshgeom_map_apply_vector(
    const GSHGEOM_Map * map, GSHGEOM_Point v, GSHGEOM_Point * out);

/**
 * Apply the map to a curve: the control points move and the weights stay (an
 * affine map preserves a rational curve's form).  The curve's dimension must
 * be the map's.  The result must be a valid curve: a point past
 * ::GSHGEOM_COORD_LIMIT is ::GSHGEOM_ERR_LIMIT.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_map_apply_curve(
    const GSHGEOM_Map * map, const GSHGEOM_Curve * curve, GSHGEOM_Curve * out);

/**
 * The route predicate (AD-5): whether a glyph under this map takes the
 * grid-aligned route.
 *
 * It holds when the map is 2D, its tag is a subset of translation, uniform
 * scale and the y flip, and its values are exactly `x' = s x + tx`,
 * `y' = -s y + ty` (the canonical design-to-model map; ::GSHGEOM_Grid
 * `flip_y` is set) or `y' = s y + ty` (a map with no flip), with `s > 0`.  It
 * rejects every rotation (a quarter turn
 * included), every skew, every other mirror, every non-uniform scale, a
 * singular map and a 3D map.  Both the tag and the values must pass: a map
 * whose tag says rotation but whose values are the identity is rejected, and
 * so is one whose tag is zero but whose values rotate.  The size `s` at which
 * a hinter takes the route is for the caller to judge; it is returned.
 *
 * @param aligned Receives 1 or 0.
 * @param grid    Receives the scale and translation when aligned; may be NULL.
 * @return ::GSHGEOM_OK (not aligned is not an error); ::GSHGEOM_ERR_INVALID
 *         for NULL or an invalid map.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_map_grid_aligned(
    const GSHGEOM_Map * map, int * aligned, GSHGEOM_Grid * grid);

/**
 * Whether the map stretches every direction of the plane, exactly: the smallest singular value of its 2x2 linear part is
 * greater than 1.
 *
 * With the columns `a = (m00, m10)` and `b = (m01, m11)` and the exact reals `A = a.a`, `B = a.b`, `C = b.b`, that is
 * `A > 1`, `C > 1` and `(A - 1)(C - 1) > B^2`, decided on the entries' own doubles by exact arithmetic (no rounding in the
 * products, no square root).  An exactly unit scale does not stretch (the smallest singular value is 1), and neither does a
 * map that stretches one direction and shrinks another: a determinant, which is an area ratio, would call a 7 by 0.3
 * map a stretch.  A 3D map is judged by the plane's linear part (`m00, m01, m10, m11`), as a 2D map is.  The translation does
 * not matter.
 *
 * @param map The map; a zero-filled map (dimension 0) is the identity.
 * @return 1 if every direction is stretched; 0 otherwise, and for NULL, an invalid map or a non-finite entry.
 */
GSHGEOM_API int gshgeom_map_stretches_every_direction(const GSHGEOM_Map * map);

/**
 * The map's hash, in the spelling of hash.h: the format version, the
 * dimension and the kind tag as `u32`, the entries row by row (`dim * dim`
 * doubles) and then the translation.  An invalid map hashes as the empty
 * stream.
 */
GSHGEOM_API uint64_t gshgeom_map_hash(const GSHGEOM_Map * map);

/**
 * Write the map as text: a header line with the version, dimension, kind tag
 * and hash, then one line per row (`m` and the row) and `t`, every double as 16
 * hexadecimal digits of its bit pattern.
 *
 * @return ::GSHGEOM_ERR_IO if the sink fails; ::GSHGEOM_ERR_INVALID for NULL or
 *         an invalid map.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_map_dump(const GSHGEOM_Map * map, FILE * out);

#ifdef __cplusplus
}
#endif

#endif /* GHOTI_IO_GSHGEOM_MAP_H */
