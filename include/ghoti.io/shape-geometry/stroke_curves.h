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
 * @file stroke_curves.h
 *
 * The region a stroke covers (stroke.h, SG6 to SG13) as an outline of curves:
 * lines, cubics and rational quadratics in target units, with a tolerance stated
 * in those units (SG32 to SG43).
 *
 * stroke.h's expansion returns 26.6 polygons accurate to 1/16 of the pixel of
 * the map it was given, which is what a rasterizer wants.  A writer (SVG, PDF)
 * that meets a style field its format cannot express must expand the stroke and
 * write the region as a fill; written as that polygon the file is correct at the
 * size it was expanded for and faceted when a viewer zooms.  This outline is the
 * same region as a path of curves: joins, caps and discs are lines and exact
 * rational-quadratic arcs, the pen of a non-uniform map costs no approximation
 * (an affine map takes each of these curves to a curve of its kind), and only
 * the offsets of curved segments are approximated, by cubics with a certified
 * bound.  `h` is half the width, `t` the caller's tolerance in target units,
 * `M` the map; the region `R` and the statuses are stroke.h's.
 *
 * ## The rules
 *
 * **SG32 What is returned.**  ::gshgeom_stroke_curves_expand takes the inputs of
 * ::gshgeom_stroke_expand plus a tolerance `t` (a `double`, target units) and a
 * form (::GSHGEOM_Stroke_Curves_Form), and returns an opaque
 * ::GSHGEOM_Stroke_Curves: a status (stroke.h's ::GSHGEOM_Stroke_Status, with SG14
 * and SG15's conditions; ::gshgeom_stroke_map_invertible is the one predicate);
 * a **body**, an **aligned** part and a **clip** part, each a 2D ::GSHGEOM_Path
 * in target units, and the clip's rule and complement flag (SG12, SG16).  Every
 * subpath of body and aligned is closed and holds only moves, closes, lines,
 * cubics and, in form ::GSHGEOM_STROKE_CURVES_RATIONAL, rational quadratics
 * (weights positive and at most 1).  **The clip part is `M` applied to the closed
 * subpaths of the centerline that have extent**, element for element, so it may
 * hold any element kind of the input: it is the fill region `F` of SG12 as the
 * centerline states it.  The parts and the status are stroke.h's, in the same
 * cases, with the same meaning; nothing is intersected here either.
 *
 * **SG33 The tolerance.**  With `D(r)` the closed disc of radius `r`, `R` the
 * region of SG6 to SG13 (SG12's parts for alignment) and `O` a part's path:
 * `R (-) D(t) <= fill(O) <= R (+) D(t)`, `fill` the nonzero fill, **the arcs of
 * `O` taken as they are, not flattened**, each part against its own region.  The
 * winding number of body and aligned is at least 0 at every point off them, which
 * is what makes nonzero their union (SG13).  Where `R_c`'s boundary crosses `F`'s
 * at an angle `alpha` the composed error may reach `t / sin(alpha / 2)`.
 *
 * **SG34 The tolerance parameter.**  `t` is finite and in
 * [::GSHGEOM_STROKE_CURVES_MIN_TOLERANCE, ::GSHGEOM_STROKE_CURVES_MAX_TOLERANCE]
 * (2^-14 to 2^4 target units); outside it is ::GSHGEOM_ERR_INVALID.
 * ::GSHGEOM_STROKE_CURVES_TOLERANCE (2^-10) is the named default.  The minimum is
 * the finest the dash cut's gate (SG26) guarantees: it refuses a cut whose
 * position error times the map's norm passes 2^-20 target units, which is `t / 64`
 * at the minimum.
 *
 * **SG35 The budget.**  The part's deviation from `R` is at most `t`, spent as:
 * approximation at most `3t / 4`, a reserve of `t / 4` for the double arithmetic
 * (at most `t / 8`, checked against the magnitudes of the input), the dash cut
 * (at most `t / 64`) and the rest.  The approximation is `t / 2` for offsets,
 * fans and the arcs join's polylines and `t / 4` for the lowering of arcs to
 * cubics (SG40).  The offsets get `t / 2` in **both** forms, so that form
 * ::GSHGEOM_STROKE_CURVES_CUBIC is exactly form ::GSHGEOM_STROKE_CURVES_RATIONAL
 * lowered by ::gshgeom_stroke_curves_lower_cubic at `t / 4`, bit for bit.
 *
 * **SG36 Width spaces.**  Local: the pieces are built in the centerline's model
 * space at the budget divided by an upper bound of `M`'s spectral norm (its
 * Frobenius norm rounded up), then `M` is applied to every control point and the
 * weights are kept: an affine map takes a line to a line, a cubic to a cubic and
 * a rational quadratic to a rational quadratic exactly, so the elliptical pen of
 * a non-uniform `M` costs no approximation beyond the build's.  Device: the
 * pieces are built around the mapped centerline (each segment mapped by
 * ::gshgeom_map_apply_curve) in target units with `h` in target units.  Joins,
 * the miter limit and the degenerate square's axis are decided where stroke.h
 * decides them (SG9, SG10, SG11).
 *
 * **SG37 The pieces.**  The outline is a sum of closed pieces whose winding is 0
 * or +1 at every point, so its nonzero fill is their union and an overlap never
 * cancels.  Exact (no approximation): the body of a line segment, a bevel, a
 * miter, a miter-clip pentagon, a square cap and the reversal's rectangle are
 * polygons of lines; a round join (the sector), a round cap (the half disc), a
 * degenerate disc and the reversal's half disc are lines and rational-quadratic
 * arcs of at most a quarter turn each (middle weight `cos(psi / 2)` from the unit
 * vectors and `gmath_sqrt`).  Approximated: the body of a **curved** segment.  A
 * piece of the segment's parameter whose band (stroke.h, SG6) is certified
 * simple is bounded by two offset cubics (the cubic Hermite interpolants of the
 * offsets `c(u) +- h n(u)` at the piece's ends) and the two end normals, and
 * consecutive pieces of a segment are written as one contour whose shared
 * normals do not appear.  A piece the certificate refuses is split; a zone where
 * the band is not simple (the offset folds, or a cusp) is covered, down to the
 * budget's extent, by the fans of stroke.h (a double sector, or a disc at a
 * cusp), whose arcs are exact.  The arcs join (stroke.h, SG8) is built as the
 * boundary of its region decomposed into convex trapezoids, its sub-arcs
 * circumscribed polylines within the budget.
 *
 * **SG38 The certified bound for an offset fit.**  The cubic Hermite
 * interpolant `C` of the offset `O(u) = c(u) + s h n(u)` (`s = +1` or `-1`, `n`
 * the unit left normal) over a parameter interval of length `d` satisfies, by the
 * interpolation remainder, `|O - C| <= (d^4 / 384) sup |O''''|`.  The supremum
 * is certified, not measured: from Bernstein bounds of the derivatives of `c` up
 * to order five (for a rational segment, through `D c = N`), a certified lower
 * bound `m` of `|c'|` on the piece (the hodograph's cone, stroke.h's) and the
 * recurrences `|s| T = c'`, `s^2 = c' . c'` that bound the derivatives of the unit
 * tangent `T`; documentation/design.md states each constant.  A piece is
 * accepted when the bound is within the offset budget, a piece whose `m` the cone
 * cannot certify away from zero is split, and a forced depth past
 * ::GSHGEOM_STROKE_MAX_DEPTH is ::GSHGEOM_ERR_LIMIT, never a coarser piece.
 *
 * **SG39 Dashes, alignment, zero width, degenerate subpaths.**  Unchanged in
 * meaning (stroke.h, SG12, SG14, SG15, SG19 to SG30).  Dash pieces are sub-curves
 * of the original segments cut at the positions stroke.h's cut finds, each
 * stroked by SG37 with caps; a join inside a dash is decided from the original
 * segments.  Alignment returns body (`R_o`), aligned (`R_c` at `2w`) and the
 * exact clip.
 *
 * **SG40 The lowering to cubics.**  ::gshgeom_stroke_curves_lower_cubic turns a
 * path of lines, quadratics, cubics and rational curves into lines and cubics.
 * Form ::GSHGEOM_STROKE_CURVES_CUBIC is the expansion followed by exactly this
 * function, at `t / 4`, on body and aligned.  A quadratic is raised to a cubic
 * exactly.  A rational quadratic with ends of weight 1 and middle weight `w < 1`
 * is a conic arc, the affine image of a circular arc of angle `psi = 2 acos w`;
 * it is split by exact bisection of the conic (a half has middle weight
 * `sqrt((1 + w) / 2)`) until a certified bound is within the tolerance, and each
 * piece becomes the cubic with the same ends and end tangents whose handles are
 * the fraction `4 w / (3 (1 + w))` of the control legs.  For a circular arc of
 * radius `r` the cubic bulges outward by at most `r e` with
 * `e = sqrt(1 + (1 - w)^3 / (27 (1 + w))) - 1 <= (1 - w)^3 / (54 (1 + w))`
 * (exact: `|P|^2 / r^2 - 1 = (3 k c - 2 s)^2 t^2 (1 - t)^2 (1 - 2 t)^2`, maximal
 * at `t(1 - t) = 1/6`), which falls as `psi^6`; for the image of the circle under
 * an affine map `A` the bound is `|A| r` times that, with `|A| r` taken from
 * the control triangle.  Any other rational segment is bisected, in homogeneous
 * coordinates, until the variation of its weights makes the polynomial curve on
 * the same control points within the tolerance.
 *
 * **SG41 Determinism, revision, hash, dump.**  The outline is a function of the
 * centerline, the canonical parameters, the map, `t` and the form, bit for bit on
 * every target: `+ - * /`, `gmath_sqrt`, the exact predicates and, in the arcs
 * join's clip only, `math`'s trigonometry, in a fixed order.  Pieces are emitted
 * in stroke.h's order (SG18, SG25).  The construction is numbered
 * (::GSHGEOM_STROKE_CURVES_REVISION, ::gshgeom_stroke_curves_revision), recorded
 * where a rebaseline lives and in no hash.  **Hash** (format 1): `u32` 1, `u32`
 * kind `0x43525453` ("STRC" little-endian), `u8` status, `u8` alignment, `u8`
 * clip rule, `u8` complement, then `u64` the ::gshgeom_path_hash of body,
 * aligned and clip (the offset basis for an absent part).  **Dump:** a line
 * `gshgeom.stroke-curves 1 status=N form=N align=N rule=N complement=N hash=H`,
 * then for each part a line `part body`, `part aligned` or `part clip` followed by
 * its path's dump, or `part <name> -` when absent.  The tolerance and the form are
 * inputs and are in no hash.  A coordinate outside ::GSHGEOM_COORD_LIMIT or an
 * element count past `max_elements` is ::GSHGEOM_ERR_LIMIT, never clipped.
 *
 * **SG42 The agreement with stroke.h's polyline.**  With `A` the fill of the body
 * of form ::GSHGEOM_STROKE_CURVES_RATIONAL flattened by
 * ::gshgeom_device_path_flatten at stroke.h's tolerance `tau` and `B` the fill of
 * stroke.h's body, both in target units under the same map: `A (-) D(d) <= B <=
 * A (+) D(d)` with `d = 2 tau + t` (`tau` is `tau / 64` target units), by the
 * triangle inequality: `B` is within `tau` of `R`, `R` within `t` of `fill(O)`,
 * and `fill(O)` within `tau` of `A`.
 *
 * **SG43 The zoom.**  The error is a constant `t` of the file's own unit, so a
 * viewer at zoom `z` sees at most `t z`, where a polygon expanded at 1x shows
 * its facets in proportion to `z` (tests/unit/test_stroke_curves_zoom.cpp and
 * tools/resvg-stroke/zoom.py render both against an external renderer at 1x,
 * 16x and 64x).
 */

#ifndef GHOTI_IO_GSHGEOM_STROKE_CURVES_H
#define GHOTI_IO_GSHGEOM_STROKE_CURVES_H

#include <ghoti.io/shape-geometry/allocator.h>
#include <ghoti.io/shape-geometry/core.h>
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

/** The finest tolerance accepted (SG34): 2^-14 target units. */
#define GSHGEOM_STROKE_CURVES_MIN_TOLERANCE 0x1p-14
/** The coarsest tolerance accepted (SG34): 2^4 target units. */
#define GSHGEOM_STROKE_CURVES_MAX_TOLERANCE 0x1p4
/** The default tolerance (SG34): 2^-10 target units, about a thousandth of a pixel. */
#define GSHGEOM_STROKE_CURVES_TOLERANCE 0x1p-10
/** The revision of the construction (SG41). */
#define GSHGEOM_STROKE_CURVES_REVISION 1

/** Which curves the outline may hold (SG32, SG40).  The values are pinned (they are dumped). */
typedef enum GSHGEOM_Stroke_Curves_Form {
  GSHGEOM_STROKE_CURVES_RATIONAL = 0, ///< Lines, cubics and rational quadratics (arcs stay exact).
  GSHGEOM_STROKE_CURVES_CUBIC = 1,    ///< Lines and cubics only: the rational form lowered by SG40.
  GSHGEOM_STROKE_CURVES_FORM_COUNT    ///< Closes the enum.
} GSHGEOM_Stroke_Curves_Form;

/** The result of ::gshgeom_stroke_curves_expand (SG32).  Opaque. */
typedef struct GSHGEOM_Stroke_Curves GSHGEOM_Stroke_Curves;

/**
 * @brief Expand a stroke into a curve outline in target units (SG32, SG33).
 *
 * @param centerline A 2D model path.
 * @param stroke     The parameters (stroke.h, SG2 checks them).
 * @param map        A 2D map from model space to target units (a stream target composes its scale into it).
 * @param tolerance  `t`, in target units, in [::GSHGEOM_STROKE_CURVES_MIN_TOLERANCE,
 *   ::GSHGEOM_STROKE_CURVES_MAX_TOLERANCE].
 * @param form       ::GSHGEOM_Stroke_Curves_Form.
 * @param alloc      The allocator; NULL for the default.
 * @param limits     `max_elements` bounds each part; NULL for the defaults.
 * @param out        Receives the outline, only on success; free it with ::gshgeom_stroke_curves_destroy.
 * @return ::GSHGEOM_OK (a not-invertible map and width 0 included); ::GSHGEOM_ERR_INVALID for NULL, an invalid
 *   map, a 3D map, a tolerance outside its range or not finite, a form out of range, and as stroke.h's SG2;
 *   ::GSHGEOM_ERR_UNSUPPORTED for a 3D centerline or a width profile; ::GSHGEOM_ERR_LIMIT as stroke.h states
 *   for the polygon expansion, and for a magnitude past what the tolerance can be certified for (the double
 *   arithmetic's reach must be within `t / 8`), a coordinate past ::GSHGEOM_COORD_LIMIT, a part past
 *   `max_elements`, a piece the depth cap would force; ::GSHGEOM_ERR_OOM.  Nothing is returned on failure.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_stroke_curves_expand(const GSHGEOM_Path * centerline,
    const GSHGEOM_Stroke * stroke, const GSHGEOM_Map * map, double tolerance, GSHGEOM_Stroke_Curves_Form form,
    const GSHGEOM_Allocator * alloc, const GSHGEOM_Limits * limits, GSHGEOM_Stroke_Curves ** out);

/** @brief Free an outline; NULL is ignored. */
GSHGEOM_API void gshgeom_stroke_curves_destroy(GSHGEOM_Stroke_Curves * curves);

/** @brief The status; ::GSHGEOM_STROKE_STATUS_COUNT for NULL. */
GSHGEOM_API GSHGEOM_Stroke_Status gshgeom_stroke_curves_status(const GSHGEOM_Stroke_Curves * curves);

/** @brief The form the outline was made in; ::GSHGEOM_STROKE_CURVES_FORM_COUNT for NULL. */
GSHGEOM_API GSHGEOM_Stroke_Curves_Form gshgeom_stroke_curves_form(const GSHGEOM_Stroke_Curves * curves);

/** @brief The body part (SG32), or NULL when absent; owned by the outline. */
GSHGEOM_API const GSHGEOM_Path * gshgeom_stroke_curves_body(const GSHGEOM_Stroke_Curves * curves);

/** @brief The aligned part (SG32), or NULL when absent; owned by the outline. */
GSHGEOM_API const GSHGEOM_Path * gshgeom_stroke_curves_aligned(const GSHGEOM_Stroke_Curves * curves);

/** @brief The clip part (SG32), or NULL when absent; owned by the outline. */
GSHGEOM_API const GSHGEOM_Path * gshgeom_stroke_curves_clip(const GSHGEOM_Stroke_Curves * curves);

/** @brief The clip part's fill rule (the canonical alignment rule). */
GSHGEOM_API GSHGEOM_Fill_Rule gshgeom_stroke_curves_clip_rule(const GSHGEOM_Stroke_Curves * curves);

/** @brief 1 when the aligned part is kept outside the clip (outside alignment), else 0. */
GSHGEOM_API int gshgeom_stroke_curves_clip_complement(const GSHGEOM_Stroke_Curves * curves);

/** @brief The hash (format 1, SG41); the offset basis for NULL. */
GSHGEOM_API uint64_t gshgeom_stroke_curves_hash(const GSHGEOM_Stroke_Curves * curves);

/**
 * @brief Write the dump (format 1, SG41).
 *
 * @return ::GSHGEOM_OK, ::GSHGEOM_ERR_INVALID for NULL, ::GSHGEOM_ERR_IO when the stream fails (including at the flush).
 */
GSHGEOM_API GSHGEOM_Result gshgeom_stroke_curves_dump(const GSHGEOM_Stroke_Curves * curves, FILE * out);

/** @brief The construction's revision (::GSHGEOM_STROKE_CURVES_REVISION). */
GSHGEOM_API int gshgeom_stroke_curves_revision(void);

/**
 * @brief SG40: lower a 2D path to lines and cubics, each curve within @p tolerance of its original.
 *
 * Moves, lines, cubics and closes are copied; a quadratic is raised exactly; a rational curve is bisected until
 * SG40's certified bound is within @p tolerance and then written as cubics.  The expansion's form
 * ::GSHGEOM_STROKE_CURVES_CUBIC is this function at a quarter of the caller's tolerance, so a writer that lowers the
 * clip part or a path of its own gets the same code.
 *
 * @param in        A 2D path.
 * @param tolerance Finite and positive, in the path's units.
 * @param alloc     The allocator; NULL for the default.
 * @param limits    `max_elements` bounds the result; NULL for the defaults.
 * @param out       Receives the lowered path, only on success; free it with ::gshgeom_path_destroy.
 * @return ::GSHGEOM_OK; ::GSHGEOM_ERR_INVALID for NULL, a 3D path or a tolerance that is not finite and positive;
 *   ::GSHGEOM_ERR_LIMIT when a curve would need more than 2^20 pieces (a depth past 20 halvings) or the result
 *   passes `max_elements`; ::GSHGEOM_ERR_OOM.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_stroke_curves_lower_cubic(const GSHGEOM_Path * in, double tolerance,
    const GSHGEOM_Allocator * alloc, const GSHGEOM_Limits * limits, GSHGEOM_Path ** out);

#ifdef __cplusplus
}
#endif

#endif /* GHOTI_IO_GSHGEOM_STROKE_CURVES_H */
