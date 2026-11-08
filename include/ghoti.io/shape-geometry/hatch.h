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
 * @file hatch.h
 *
 * Hatching: families of parallel lines clipped to a closed region, cut into dashes and stroked, stated once (HX4 to
 * HX12, HX16).  A hatch is geometry, not a raster pattern: a plotter follows its lines, a CAD writer emits them as a
 * hatch entity, and a pixel sink draws the outline this file expands them into.
 *
 * A **hatch** is a set of **families**.  A family is a pattern definition line of the DXF HATCH entity (the Autodesk DXF
 * Reference, Pattern Data: groups 53, 43 and 44, 45 and 46, 79 and 49; the sense of each value is the one of the AutoCAD
 * hatch pattern definition `angle, x-origin, y-origin, delta-x, delta-y, dash-1, dash-2, ...`: a line is of infinite length, `delta-x`
 * is the displacement between members of the family in the direction of the line and matters to dashed lines only, `delta-y`
 * the spacing between members, measured perpendicular to the lines, and a dash length is positive for a pen-down dash,
 * negative for a space and zero for a dot).  The values are in **hatch space**, which a map takes to the local space of the
 * item that is hatched.
 *
 * ## The rules
 *
 * **HX4 The trim region.**  `P` is the region as closed contours of 26.6 points (an open subpath closes), made by the
 * caller or from a model path by the chain and the flattener at the tolerance.  For a point `p`, `c(p)` is the number of
 * contours whose own even-odd interior contains `p` and `w(p)` the sum of their winding numbers.  The region `G` is, for
 * ::GSHGEOM_HATCH_NESTING_RULE, `{w != 0}` under the nonzero rule or `{w odd}` under the even-odd rule; for
 * ::GSHGEOM_HATCH_NESTING_OUTER `{c = 1}`; for ::GSHGEOM_HATCH_NESTING_IGNORE `{c >= 1}` (the item's rule is not used for
 * the last two).  For a family line `L` the **trim set** is `closure(G) cap L` and the **trim intervals** are its
 * connected components of positive length: a point where a line only touches the boundary is no interval, and a line
 * that lies along a boundary edge is in the set exactly where the edge has an interior side, once (the decision is that of
 * the line shifted infinitesimally to either side: every side test is exact, only a crossing's position is computed in
 * `double`).  Two crossings of one line whose positions differ by at most 2^-26 pixel are one position and an interval or
 * a gap that short is no interval and no gap: the positions are computed, not exact, and a vertex shared by two edges is
 * one position.
 *
 * **HX5 The lines of a family and the exact decisions.**  In hatch space `u = (cos a, sin a)` and `n = (-sin a, cos a)`
 * with `a` the family's angle in the hatch's unit, by ::gshgeom_sincos (exact at a quarter turn: 0 and +-1, never -0).
 * Line `k` (every integer) is `L_k(s) = O_k + s u` with `O_k = O + (k dx) u + (k dy) n`, `s` real, so the stagger `dx`
 * moves each line's dashes along it and the spacing `dy` moves the lines apart; `O`, the origin, is applied in hatch space
 * **before** the map.  `O_k` is formed in `double`, left to right, each product and sum rounded once, no fused
 * multiply-add, from the integer `k` converted exactly; nothing accumulates.  The line in device pixels is the line through the
 * double point `A_k = C O_k` with the double direction `w = C_lin u` (`C` the chain composed with the hatch map, formed once
 * by ::gshgeom_map_compose), and `s` is its parameter, so a distance on the line is in hatch units.  Which side of that line a
 * contour point lies on, and whether it lies on it, is decided exactly (the sign of `w x (P - A_k)`, an exact sum of
 * products).  The lines of the use are those whose line meets the closed bounding box of `P` (an exact test of the four
 * corners), for the families in order and `k` increasing; a family with `delta-y = 0` has none.  Their number is the **lines of the
 * use**.
 *
 * **HX6 Dash positions: the continuous pattern.**  A family's elements are lowered to a dash array `A` of period `P_A` and
 * an offset `o` by ::gshgeom_stroke_dashes_from_dxf with scale 1; the pattern at line distance `s` is at `(o + s) mod P_A`
 * (floored), running without a break along the whole line and in both directions, and it does **not** restart at a trim
 * interval.  A family whose elements are all dashes (no space), or that has none, is solid: its pieces are its trim
 * intervals.  For a dashed one, per line, with `t_first` and `t_last` the ends of the first and last trim interval,
 * ::gshgeom_stroke_dash_each gives the dashes `<a, b>` of the window that starts one period before `t_first` and ends one
 * period after `t_last` (offset `o + (t_first - P_A)`, length `(t_last - t_first) + 2 P_A`), at `s = (t_first - P_A) + a` and `+ b`, so
 * that the first interval sees the unbroken pattern, SG21's half-open tie at the start of a window falls outside every interval,
 * and a dot that lies exactly at the end of the last one is made (SG21 makes no entry that begins at the end of its subpath).
 * Each dash is intersected with the trim intervals in order: a dash cut by an interval to a length of at most 2^-26 pixel (HX4's band) is dropped (a
 * boundary that only touches a dash does not make a dot); a **dot** (a dash of no length) inside an interval, its ends included, is kept, and
 * is kept at an end to within 2^-46 of the largest of the period, 1 and the magnitude of the line's first and last trim position (and then placed on the interval), because a pattern that begins with a
 * space and ends with a dot has its line start one unit in the last place after the dot (SG30), so the dot that belongs at
 * a multiple of the period lies that much before it.  The **pieces** of a line are the resulting `<s0, s1>` in increasing
 * `s`, a piece with `s0 = s1` a **dot**; they are emitted family by family in the given order, each family's lines by
 * increasing `k`.
 *
 * **HX7 The stroke of a piece.**  A piece with `s0 < s1` is the straight open centerline from `H(O_k + s0 u)` to
 * `H(O_k + s1 u)` in the item's local space (`H` the hatch map), stroked by stroke.h's SG6 and SG7 with the style's width,
 * width space and cap at both ends; there is no join, a piece being one segment.  A dot is SG24's shape along the mapped
 * line direction `H_lin u`: the disc of radius `h` for a round cap, the square of side `2h` for a square cap, nothing for
 * a butt cap, drawn whatever the style's `degenerate` says.  All the pieces of all the families are one union (SG13), so
 * translucent ink does not darken where two lines cross.  The expansion is **the code of ::gshgeom_stroke_expand_with_tolerance**
 * run on the pieces as an explicit path, one open two-point subpath each (a dot a degenerate subpath carrying its direction),
 * in the order HX6 gives; there is no second construction of a body or a cap.
 *
 * **HX8 The functions.**  ::gshgeom_hatch_pieces and ::gshgeom_hatch_expand take the region as a model path and the
 * chain's map from local space to device pixels; the `_from_flat` forms take the already flattened contours (a sink with a
 * glyph run has them).  ::gshgeom_hatch_expand returns a ::GSHGEOM_Stroke_Outline whose body is the pieces' and which has
 * no aligned or clip part; the hash and dump are stroke.h's.  The status is ::GSHGEOM_STROKE_NOT_INVERTIBLE when the chain or the
 * hatch map has a linear part of determinant exactly zero (each tested alone, SG15's predicate) and
 * ::GSHGEOM_STROKE_EMPTY for width 0 (SG14).
 *
 * **HX10 Tolerance.**  A piece end lies within the flattening tolerance of the true boundary crossing (the boundary is
 * the flattener's polygon), and the stroke of a piece is within SG17.
 *
 * **HX11 Hatch space.**  Dash lengths, origins and offsets are in hatch space and are mapped by the hatch map and the chain whichever the
 * style's width space; the width is the style's.
 *
 * **HX12 Limits.**  The structural bounds (::gshgeom_hatch_check), `max_lines` (the lines of the use, over all families;
 * more is ::GSHGEOM_ERR_LIMIT, decided before a line past it is generated), the pieces (at most the limits' `max_elements` / 2)
 * and the dash steps (the limits' `max_dash_steps`, as the stroke's) and SG26's gate for a dash against its period.  The limits'
 * `max_work` bounds the work the count of lines multiplies: each line reads every point of the region (a unit for each 16
 * points, after the first 256 units of the line) and walks its dash pattern across its length (3 units a step of the walk, after
 * the first 256 units of the line), and the stroke of the pieces spends the stroke's (stroke.h).  A call that makes the pieces and
 * strokes them (::gshgeom_hatch_expand and its variants) spends one budget over both stages; the pieces' and the stroke's entries
 * called one after the other are two calls and have a budget each.  More is ::GSHGEOM_ERR_LIMIT (documentation/design.md,
 * "Work budget").
 *
 * **HX16 Revision.**  ::GSHGEOM_HATCH_REVISION is 1 and is in no hash.
 *
 * ## The pieces' hash and dump (format 1)
 *
 * The hash, in the spelling of hash.h: `u32` 1, `u32` kind `0x43435048` ("HPCC" little-endian), `u64` the piece count,
 * then per piece `u8` the dot flag and `f64` the four coordinates (the second point is the first for a dot).  The dump: a
 * line `gshgeom.hatch-pieces 1 pieces=N lines=N hash=H` and a line `piece dot|seg x0 y0 x1 y1` per piece, each double as 16 hex digits.
 */
#ifndef GHOTI_IO_GSHGEOM_HATCH_H
#define GHOTI_IO_GSHGEOM_HATCH_H

#include <ghoti.io/shape-geometry/allocator.h>
#include <ghoti.io/shape-geometry/core.h>
#include <ghoti.io/shape-geometry/fill.h>
#include <ghoti.io/shape-geometry/flatten.h>
#include <ghoti.io/shape-geometry/macros.h>
#include <ghoti.io/shape-geometry/map.h>
#include <ghoti.io/shape-geometry/path.h>
#include <ghoti.io/shape-geometry/stroke.h>

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/** The most families of a hatch (HX2): 256. */
#define GSHGEOM_HATCH_FAMILIES_MAX ((size_t)256)
/** The most dash elements of a family (HX2): 65535. */
#define GSHGEOM_HATCH_DASHES_MAX ((size_t)65535)
/** The largest angle magnitude, in its unit (HX2): 2^40. */
#define GSHGEOM_HATCH_ANGLE_MAX 0x1p40
/** The smallest magnitude of a nonzero origin, offset or dash (HX2): 2^-64; the largest is 2^64. */
#define GSHGEOM_HATCH_VALUE_MIN 0x1p-64
/** The largest magnitude of an origin, offset or dash (HX2): 2^64. */
#define GSHGEOM_HATCH_VALUE_MAX 0x1p64
/** The lines of a use a caller that states no bound allows: 2^16. */
#define GSHGEOM_HATCH_LINES_DEFAULT ((uint64_t)1 << 16)
/** The shortest crossing gap or interval that counts (HX4), in device pixels: 2^-26. */
#define GSHGEOM_HATCH_BAND 0x1p-26
/** The revision of the construction (HX16). */
#define GSHGEOM_HATCH_REVISION 1

/** How loops nest (HX4).  The values are pinned (they are hashed by the scene). */
typedef enum GSHGEOM_Hatch_Nesting {
  GSHGEOM_HATCH_NESTING_RULE = 0,   ///< The fill rule given.
  GSHGEOM_HATCH_NESTING_OUTER = 1,  ///< Only the outermost ring (DXF group 75, style 1).
  GSHGEOM_HATCH_NESTING_IGNORE = 2, ///< Everything inside the outer loops (DXF group 75, style 2).
  GSHGEOM_HATCH_NESTING_COUNT       ///< Closes the enum.
} GSHGEOM_Hatch_Nesting;

/** One family (HX1): a pattern definition line. */
typedef struct GSHGEOM_Hatch_Family {
  double angle;           ///< The line's angle, in the hatch's unit, from +x toward +y of hatch space.
  double origin_x;        ///< A point of line 0, in hatch space.
  double origin_y;
  double offset_along;    ///< `delta-x`: the shift of each line's dashes along it, line to line.
  double offset_across;   ///< `delta-y`: the spacing between lines; 0 draws nothing.
  size_t dash_count;      ///< The elements (DXF group 49); 0 is solid.
  const double * dashes;  ///< Positive a dash, negative a space, zero a dot.  Borrowed; NULL only when @p dash_count is 0.
} GSHGEOM_Hatch_Family;

/** A hatch (HX1), plain data. */
typedef struct GSHGEOM_Hatch {
  uint32_t angle_unit;                   ///< A ::GSHGEOM_Angle_Unit.
  uint32_t nesting;                      ///< A ::GSHGEOM_Hatch_Nesting.
  GSHGEOM_Map map;                       ///< From hatch space to the local space of the hatched item; 2D.
  size_t family_count;                   ///< At least 1.
  const GSHGEOM_Hatch_Family * families; ///< Borrowed.
} GSHGEOM_Hatch;

/** The pieces of a hatch (HX6), in local space.  Opaque. */
typedef struct GSHGEOM_Hatch_Pieces GSHGEOM_Hatch_Pieces;

/**
 * @brief HX2's geometric checks of a hatch, in order: ::GSHGEOM_ERR_INVALID (NULL, a non-finite double, an enumeration
 *   out of range, no family, dashes counted and not given, an invalid map, a family whose elements have no dash or no period);
 *   ::GSHGEOM_ERR_UNSUPPORTED (a 3D map); then ::GSHGEOM_ERR_LIMIT (the families, the dashes, the bands of every value, the
 *   map's entries, the angle).
 */
GSHGEOM_API GSHGEOM_Result gshgeom_hatch_check(const GSHGEOM_Hatch * hatch);

/**
 * @brief The pieces of a hatch over a region given as a model path (HX4 to HX6, HX12), at the default tolerance.
 *
 * @param region A 2D model path in the item's local space; open subpaths close.
 * @param rule   The fill rule, read by ::GSHGEOM_HATCH_NESTING_RULE only.
 * @param chain  A 2D map from local space to device pixels (already composed).
 * @param max_lines The most lines of the use; 0 is ::GSHGEOM_HATCH_LINES_DEFAULT.
 * @param out    Receives the pieces, only on success.  A chain or a hatch map that is not invertible gives no pieces and
 *   no lines (the caller reports it).
 * @return ::GSHGEOM_OK; ::GSHGEOM_ERR_INVALID (NULL, an invalid hatch, map or rule); ::GSHGEOM_ERR_UNSUPPORTED (a 3D path or map);
 *   ::GSHGEOM_ERR_LIMIT (more lines than @p max_lines, more pieces than the limits, the dash steps, a cut error past SG26's
 *   gate); ::GSHGEOM_ERR_OOM.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_hatch_pieces(const GSHGEOM_Path * region, GSHGEOM_Fill_Rule rule, const GSHGEOM_Hatch * hatch,
    const GSHGEOM_Map * chain, uint64_t max_lines, const GSHGEOM_Allocator * alloc, const GSHGEOM_Limits * limits,
    GSHGEOM_Hatch_Pieces ** out);

/** @brief As ::gshgeom_hatch_pieces at a tolerance, in 26.6 units, of at least ::GSHGEOM_FLATTEN_MIN_TOLERANCE. */
GSHGEOM_API GSHGEOM_Result gshgeom_hatch_pieces_with_tolerance(const GSHGEOM_Path * region, GSHGEOM_Fill_Rule rule,
    const GSHGEOM_Hatch * hatch, const GSHGEOM_Map * chain, int32_t tolerance, uint64_t max_lines, const GSHGEOM_Allocator * alloc,
    const GSHGEOM_Limits * limits, GSHGEOM_Hatch_Pieces ** out);

/** @brief As ::gshgeom_hatch_pieces with the region as flattened device contours (HX4's `P`), which the chain maps back to local space. */
GSHGEOM_API GSHGEOM_Result gshgeom_hatch_pieces_from_flat(const GSHGEOM_Flat_Path * region, GSHGEOM_Fill_Rule rule,
    const GSHGEOM_Hatch * hatch, const GSHGEOM_Map * chain, uint64_t max_lines, const GSHGEOM_Allocator * alloc,
    const GSHGEOM_Limits * limits, GSHGEOM_Hatch_Pieces ** out);

/** @brief Free the pieces; NULL is ignored. */
GSHGEOM_API void gshgeom_hatch_pieces_destroy(GSHGEOM_Hatch_Pieces * pieces);
/** @brief How many pieces (0 for NULL). */
GSHGEOM_API size_t gshgeom_hatch_pieces_count(const GSHGEOM_Hatch_Pieces * pieces);
/** @brief The lines of the use (HX5), the count HX12 bounds (0 for NULL). */
GSHGEOM_API uint64_t gshgeom_hatch_pieces_lines(const GSHGEOM_Hatch_Pieces * pieces);
/**
 * @brief Piece @p index: its two points in local space as `xy = {x0, y0, x1, y1}`, and @p dot nonzero for a dot (then both
 *   points are the dot's).
 * @return ::GSHGEOM_OK; ::GSHGEOM_ERR_INVALID for NULL or an index past the end.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_hatch_pieces_get(const GSHGEOM_Hatch_Pieces * pieces, size_t index, double xy[4], int * dot);
/** @brief The pieces' hash (format 1, this file's comment); the offset basis for NULL. */
GSHGEOM_API uint64_t gshgeom_hatch_pieces_hash(const GSHGEOM_Hatch_Pieces * pieces);
/** @brief Write the pieces' dump (format 1); ::GSHGEOM_ERR_INVALID for NULL, ::GSHGEOM_ERR_IO when the stream fails. */
GSHGEOM_API GSHGEOM_Result gshgeom_hatch_pieces_dump(const GSHGEOM_Hatch_Pieces * pieces, FILE * out);

/**
 * @brief The expansion of a hatch at the default tolerance (HX7, HX8): its pieces stroked by @p stroke, as one outline.
 *
 * @p stroke must have no dashes and center alignment (otherwise ::GSHGEOM_ERR_INVALID) and is otherwise checked as
 * ::gshgeom_stroke_check does.  Its join, miter limit and `degenerate` are carried and unused.  Errors as
 * ::gshgeom_hatch_pieces and ::gshgeom_stroke_expand.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_hatch_expand(const GSHGEOM_Path * region, GSHGEOM_Fill_Rule rule, const GSHGEOM_Hatch * hatch,
    const GSHGEOM_Stroke * stroke, const GSHGEOM_Map * chain, uint64_t max_lines, const GSHGEOM_Allocator * alloc,
    const GSHGEOM_Limits * limits, GSHGEOM_Stroke_Outline ** out);

/** @brief As ::gshgeom_hatch_expand at a tolerance in 26.6 units. */
GSHGEOM_API GSHGEOM_Result gshgeom_hatch_expand_with_tolerance(const GSHGEOM_Path * region, GSHGEOM_Fill_Rule rule,
    const GSHGEOM_Hatch * hatch, const GSHGEOM_Stroke * stroke, const GSHGEOM_Map * chain, int32_t tolerance, uint64_t max_lines,
    const GSHGEOM_Allocator * alloc, const GSHGEOM_Limits * limits, GSHGEOM_Stroke_Outline ** out);

/** @brief As ::gshgeom_hatch_expand with the region as flattened device contours (the tolerance is the contours' own). */
GSHGEOM_API GSHGEOM_Result gshgeom_hatch_expand_from_flat(const GSHGEOM_Flat_Path * region, GSHGEOM_Fill_Rule rule,
    const GSHGEOM_Hatch * hatch, const GSHGEOM_Stroke * stroke, const GSHGEOM_Map * chain, uint64_t max_lines,
    const GSHGEOM_Allocator * alloc, const GSHGEOM_Limits * limits, GSHGEOM_Stroke_Outline ** out);

/**
 * @brief Stroke already-made pieces (the second half of ::gshgeom_hatch_expand): for a caller that wants the pieces as well
 *   (a test sink records their hash).  @p chain and @p tolerance are those the pieces were made with.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_hatch_expand_pieces(const GSHGEOM_Hatch_Pieces * pieces, const GSHGEOM_Stroke * stroke,
    const GSHGEOM_Map * chain, int32_t tolerance, const GSHGEOM_Allocator * alloc, const GSHGEOM_Limits * limits,
    GSHGEOM_Stroke_Outline ** out);

/** @brief The construction's revision (::GSHGEOM_HATCH_REVISION). */
GSHGEOM_API int gshgeom_hatch_revision(void);

#ifdef __cplusplus
}
#endif

#endif /* GHOTI_IO_GSHGEOM_HATCH_H */
