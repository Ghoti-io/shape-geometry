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
 * @file stroke.h
 *
 * What a stroke covers, stated once (SG1 to SG30), and the expansion of a
 * stroke into an outline whose nonzero fill is that region within the device
 * tolerance.
 *
 * A stroke is a centerline (a 2D model path) and a set of parameters.  The
 * region it covers is defined below in real arithmetic, by rules numbered SG1 to
 * SG18 (the numbering is stable: tests, planted defects, the generators and
 * documentation/design.md cite it).  The rules follow SVG 2 (the W3C Candidate
 * Recommendation of 2018-10-04, painting and paths chapters, which is the one
 * text that defines the miter-clip and arcs joins) and the WHATWG HTML
 * standard's canvas line styles; documentation/design.md records each reading
 * of a point those texts leave open.  `h` is half the width throughout.
 * Orientation (the left side, a positive turn) is orient.h's, cited and not
 * restated.
 *
 * ## The rules
 *
 * **SG1 The parameters.**  ::GSHGEOM_Stroke is plain data, checked by every
 * function that takes it.  Its enumerations carry fixed values (the scene's
 * stroke style uses the same ones, so a conversion is the identity), and
 * ::gshgeom_stroke_init writes the defaults: width 1, local width, miter join,
 * miter limit 4, butt caps, degenerate subpaths drawn, no dashes, offset 0, dashes
 * restart per subpath, centered, nonzero alignment rule, no width profile.  A
 * zero-filled structure is not valid.
 *
 * **SG2 Validity and bounds.**  ::gshgeom_stroke_check decides, in this order:
 * ::GSHGEOM_ERR_INVALID for NULL, `dash_count > 0` with NULL `dashes`, a
 * non-finite double, an enumeration out of range, `width < 0`, or a nonzero
 * `_reserved` byte; ::GSHGEOM_ERR_UNSUPPORTED for a non-NULL `width_profile`;
 * then ::GSHGEOM_ERR_LIMIT: `width`, each dash (by magnitude) and `dash_offset`
 * must be 0 or of magnitude in [::GSHGEOM_STROKE_VALUE_MIN,
 * ::GSHGEOM_STROKE_VALUE_MAX], `miter_limit` in [1,
 * ::GSHGEOM_STROKE_MITER_LIMIT_MAX], and the canonical dash count (SG3) at most
 * ::GSHGEOM_STROKE_DASHES_MAX.
 *
 * **SG3 Canonical meaning.**  The outline is a function of the canonical
 * parameters alone: a dash array with a negative length, or whose lengths are
 * all zero, is solid; an odd-length array is repeated once; a field that does
 * not apply has no effect (`miter_limit` under a join other than miter,
 * miter-clip and arcs; the dash fields with no dashes; `alignment_rule` under
 * center; `degenerate` under a butt cap, which draws nothing either way); `-0`
 * is `+0`.
 *
 * **SG4 Inputs and spaces.**  The centerline is a 2D ::GSHGEOM_Path (a 3D path
 * is ::GSHGEOM_ERR_UNSUPPORTED).  The map is a 2D ::GSHGEOM_Map from the
 * centerline's model space to device space, already composed (map.h); device
 * space is measured in pixels.  A stream target composes its unit-to-pixel scale
 * into the map.  The tolerance `tau` is in 26.6 units: ::GSHGEOM_FLATTEN_TOLERANCE
 * (4, 1/16 pixel) by default, at least ::GSHGEOM_FLATTEN_MIN_TOLERANCE (2).
 *
 * **SG5 Subpaths, segments, vertices, ends.**  A subpath runs from a move to the
 * next move; a close adds its closing line (the line back to the subpath's
 * start), and a drawing verb after a close starts a subpath at the closed one's
 * start.  A segment **has extent** when some control point differs, exactly,
 * from its first; a segment without extent is dropped before anything else.  A
 * **vertex** is the point between two consecutive segments with extent, and, for
 * a closed subpath with extent, its start (between the last and the first).  The
 * **ends** are the first and last points of an open subpath with extent.  The
 * **start direction** of a segment is `P_k - P_0` for the least `k` with
 * `P_k != P_0`; its **end direction** `P_n - P_k` for the greatest `k` with
 * `P_k != P_n`.  The left normal is ::gshgeom_normal_left.
 *
 * **SG6 The body.**  A segment `c` with extent covers the closure of the union,
 * over every `t` with `c'(t) != 0`, of the normal segment
 * `{c(t) + s n(t) : |s| <= h}` (`n` the unit left normal), and, at each `t*` in
 * (0, 1) with `c'(t*) = 0` (a cusp), the closed disc of radius `h` at `c(t*)`:
 * the limit of the normals of every nearby curve, which sweep a half turn there.
 *
 * **SG7 Caps.**  At an end `E` with outward unit tangent `d` (the end direction
 * at a subpath's end, the negated start direction at its start) and left normal
 * `m` of `d`: butt adds nothing; square the rectangle with corners `E + h m`,
 * `E + h m + h d`, `E - h m + h d`, `E - h m`; round the half disc
 * `{E + v : |v| <= h, v . d >= 0}`.
 *
 * **SG8 Joins.**  At a vertex `V` with incoming direction `a` (the end direction
 * of the segment before) and outgoing `b` (the start direction of the one
 * after), `Q` and `S` the control points that define them, and `sigma` the exact
 * sign of `orient2d(Q, V, S)`: `sigma = 0` with `a . b > 0` is a smooth vertex
 * and has no join; `sigma = 0` with `a . b < 0` is a reversal (below).
 * Otherwise the outer side is the left when `sigma > 0` and the right when
 * `sigma < 0`; `o_a`, `o_b` are the unit outer normals of `a` and `b`,
 * `A = V + h o_a`, `B = V + h o_b`, `phi` in (0, pi) the turn angle,
 * `u` the unit bisector of `o_a` and `o_b`, `L` the miter limit.  Every join
 * contains the triangle `V A B`.
 *  - bevel: the triangle `V A B`.
 *  - miter: within the limit (SG9), the quadrilateral `V A T B`, `T` the meeting
 *    point of the outer edges' tangent lines (`T = V + (h / cos(phi/2)) u`);
 *    otherwise the bevel.
 *  - miter-clip: within the limit, the miter; otherwise the pentagon
 *    `V A A' B' B`, `A'` and `B'` those tangent lines' intersections with the
 *    line `{x : (x - V) . u = L h}`.
 *  - round: the circular sector of radius `h` at `V` from `o_a` to `o_b` on the
 *    outer side (not a full disc).
 *  - arcs (SVG 2 CR 2018-10-04, 13.5.4 to 13.5.9): the curvatures `kappa_a` (the
 *    incoming segment's at its end) and `kappa_b` (the outgoing one's at its
 *    start) are `((n - 1) / n) (w_0 w_2 / w_1^2) |d_1 x d_2| / |d_1|^3` over the
 *    three control points nearest that end (`d_1` the first difference from the
 *    end, `d_2` the next; weights 1 for a polynomial kind; a line has 0), signed
 *    positive when the centre of curvature lies on the join's inner side; where
 *    the end's first difference vanishes the extension is a line.  If both are
 *    zero the join is miter-clip's; if either `|kappa| > 1 / h` (decided exactly)
 *    it is round's.  Otherwise each outer edge is extended beyond `A` (forward,
 *    along `a`) and `B` (backward, along `-b`) by the circle of curvature
 *    `kappa / (1 + h kappa)` tangent to it there, centre `V - o / kappa` (a line
 *    for a zero curvature).  `X` is their meeting point closest to `V`; where
 *    they do not meet, the radii are adjusted as the CR's 13.5.9 states (if one
 *    circle encloses the other the larger shrinks and the smaller grows by the
 *    same amount, if they are apart both grow, a circle and a line: the circle
 *    grows; each circle stays tangent at its own point) by the least amount that
 *    makes them touch, and `X` is the point of touching; where the offset
 *    tangents are parallel and nothing can meet, the join is the rectangle of
 *    width `2h` and length `L h` ahead of `V`.  The region is the bevel together
 *    with what the closed boundary `V`, `A`, the first extension to `X`, the
 *    second extension back to `B`, `V` encloses by nonzero winding (the CR's
 *    lines and arcs; at a reversal the boundary can curl past `V`).  `Gamma` is the circle tangent to `u` at `V` through `X`; the
 *    miter length `s` is `Gamma`'s arc length from `V` to `X`; when
 *    `s > L h` the region is cut by the line perpendicular to `Gamma` at the
 *    point `Y` at arc length `L h`, keeping `V`'s side.
 *  - reversal (`phi = pi`): the limit of each join as `phi` tends to `pi`: bevel
 *    and miter add nothing; round the half disc `{V + v : |v| <= h, v . a >= 0}`;
 *    miter-clip the rectangle from `V +- h n_a` to `V +- h n_a + L h a^`; arcs the
 *    union of its joins with the outer side taken on the left and on the right.
 *
 * **SG9 The miter limit, decided exactly.**  The miter is within the limit when
 * `1 / sin(theta / 2) <= L` (`theta = pi - phi` the angle between the
 * segments), that is `(L^2 - 2) |a||b| >= -L^2 (a . b)` for the unnormalized
 * `a`, `b`, decided by the signs of its sides and the comparison of their
 * squares, polynomials in `L` and the raw coordinates of `Q`, `V`, `S` whose
 * signs are taken exactly.  No angle is computed.
 *
 * **SG10 Width spaces.**  Local: the region is built by SG6 to SG9 and SG11 in the
 * centerline's model space and mapped by the map (a non-uniform map makes the
 * elliptical pen; joins, the miter limit and the degenerate shape's axis are
 * decided in model space).  Device: the region is built by the same rules around
 * the centerline mapped segment by segment (::gshgeom_map_apply_curve), with `h`
 * in pixels.
 *
 * **SG11 Degenerate subpaths.**  A subpath with a drawing verb or a close and no
 * segment with extent is degenerate at its point `P`.  Drawn: a round cap gives
 * the disc of radius `h` at `P`; a square cap the square of side `2h` centred at
 * `P` with two sides parallel to the effective tangent `d`: the end direction of
 * the last segment with extent before the subpath in the path's element order,
 * else the start direction of the first segment with extent after it, else
 * (no segment with extent anywhere in the path) local `+x`; decided on the raw
 * coordinates of model space under local width and of the mapped curves under
 * device width (whose `+x` fallback is the map's image of `+x`).  Butt draws
 * nothing.  Omitted: nothing.  A lone move (no drawing verb, no close) is no
 * subpath and draws nothing.
 *
 * **SG12 Alignment.**  Center: the region of SG6 to SG11 at width `w`.  Inside and
 * outside apply to the closed subpaths with extent: with `R_c` their centered
 * region at width `2w`, `R_o` the centered region of every other subpath at `w`,
 * and `F` the fill region, under `alignment_rule`, of the closed subpaths with
 * extent mapped by the map, inside is `R_o u (R_c n F)` and outside
 * `R_o u (R_c \ F)`.  The expansion returns the three separately and intersects
 * nothing: a sink composes them.
 *
 * **SG13 One union.**  The region is a set: a point is covered once however many
 * pieces, segments or subpaths contain it.
 *
 * **SG14 Zero width.**  Width 0 gives the status ::GSHGEOM_STROKE_EMPTY and no
 * outline, degenerate subpaths included.
 *
 * **SG15 Not invertible.**  For both width spaces, when the 2 by 2 linear part of
 * the map has determinant exactly zero (the exact sign of
 * `orient2d((0, 0), (m00, m10), (m01, m11))`), the status is
 * ::GSHGEOM_STROKE_NOT_INVERTIBLE, every part is absent and the call succeeds.
 * ::gshgeom_stroke_map_invertible is the predicate.  It is decided before
 * SG14.
 *
 * **SG16 The outline.**  An opaque ::GSHGEOM_Stroke_Outline: a status; **body**,
 * a flattened path (26.6, every subpath closed) whose nonzero fill is the
 * region for center and `R_o` otherwise; **aligned**, a flattened path for `R_c`
 * (absent when no closed subpath with extent is aligned); **clip**, the
 * flattened path of `F` made by ::gshgeom_device_path_from_model and the
 * flattener at the same tolerance, with its rule and a complement flag (1 for
 * outside).  The winding number of body and of aligned is at least 0 at every
 * point, which is what makes nonzero their union.  Every vertex is rounded by
 * ::gshgeom_round_26_6.
 *
 * **SG17 Tolerance.**  In device space, with `D(r)` the closed disc of radius `r`:
 * `R - D(tau) <= fill(O) <= R + D(tau)` (erosion and dilation) for body and
 * for aligned, each against its own region.  The expansion keeps the stronger
 * `R <= fill(O) <= R + D(tau - 1/4)`, proved above src/stroke/emit.c and in
 * documentation/design.md.  Where the aligned part's boundary crosses `F`'s at
 * an angle `alpha`, the composed error may reach `tau / sin(alpha / 2)`.
 *
 * **SG18 Determinism, revision and hash.**  The outline is a function of the
 * centerline, the canonical parameters, the map and `tau`, the same bits on
 * every target: `+ - * /`, `gmath_sqrt`, the exact predicates and, in the arcs
 * join's clip only, `math`'s trigonometry, in a fixed order.  Pieces are emitted
 * in a fixed order: subpaths in order; in a subpath the start cap, then for each
 * segment its body pieces and the join at its end vertex (a closed subpath's
 * start join after its last segment), then the end cap; a degenerate shape in
 * its subpath's place.  The construction is numbered
 * (::GSHGEOM_STROKE_REVISION, ::gshgeom_stroke_revision), recorded where a
 * rebaseline lives and in no hash.
 *
 * **SG19 The pattern.**  The canonical dash array (SG3) has an even count and a
 * positive sum `P`, the period; an entry of an even index is a dash, of an odd
 * index a gap.  Lengths are measured along the centerline in the width space's
 * unit: under local width the arc length of the model path (before the map), under
 * device width the arc length in pixels of the centerline mapped by the map.  A
 * solid style is SG6 to SG13 unchanged.
 *
 * **SG20 The distance.**  The distance along a subpath is along.h's: the sum of
 * the arc lengths of its segments with extent in element order, the closing
 * segment of a closed subpath included, arc length at the finest tolerance
 * along.h accepts (2^-46).  The distance along the path is the sum over its
 * subpaths in order; a move adds nothing.  A subpath's `begin` is the path
 * distance at its first point; its length is `L` (0 for a degenerate subpath).
 *
 * **SG21 Dash positions** (SVG 2 Candidate Recommendation of 2018-10-04, 13.5.7,
 * "dash positions").  For a subpath of length `L` the phase `o` is `dash_offset`
 * under ::GSHGEOM_STROKE_DASH_PER_SUBPATH (the pattern restarts at each subpath)
 * and `dash_offset + begin` under ::GSHGEOM_STROKE_DASH_PATH (one pattern for the
 * path, cut at the subpath's ends), reduced to `[0, P)` by the floored modulo.
 * With the prefix sums `S_i = d_0 + ... + d_i`, `index` is the least `i` with
 * `S_i >= o`.  The first piece is `<0, min(S_index - o, L)>`, a dash when `index`
 * is even, except that a positive dash the phase has used up entirely (`S_index =
 * o`) is over and leaves no dash of no length at the start: the one tie where this
 * differs from the text's loop, which would make one (a dot under a round or square
 * cap; every shipping renderer draws nothing, and a DXF linetype that begins with a
 * space could not be lowered).  A dash of no length exactly at the phase is kept.
 * Each following entry in turn, wrapping to entry 0 after the last, starts at
 * `position` and has length `min(d_i, L - position)`, its pair
 * `<position, position + length>` a dash when `i` is even, and the
 * sequence stops at the first entry that starts at or after `L`: a dash that
 * would begin exactly at the end of the subpath is not made.  For `L = 0` the
 * first piece only.
 *
 * **SG22 The shape of one dash.**  The dash `<a, b>` is the stretch of the
 * centerline between the distances `a` and `b`, stroked as an open subpath by SG5
 * to SG10: the body of the stretch, a cap at each end in the direction the
 * centerline has there, and the join at every vertex strictly inside the stretch,
 * decided from the original segments that meet at the vertex (never from the cut
 * pieces).  A cut exactly at a vertex is an end: a cap, not a join (SG23 adds the
 * join of a dash across the start of a closed subpath).
 *
 * **SG23 Closed subpaths.**  The closing line is a segment of the subpath: it is
 * measured and dashes run along it.  The distance ends at the start point.  A dash
 * that crosses the start point is one dash: when the first piece began before the
 * path (the phase is inside its dash) and the last dash would go on past the end of
 * the subpath, they are joined, the last dash running on through the closing line into
 * the first segment with a join (SG8) at the start vertex and a cap only at its two
 * ends; it is emitted where the last dash is.  A dash that ends exactly at the start
 * point, or a first piece that begins exactly there, stays separate, each with its
 * cap.  A dash that would cover the whole subpath is the stroke of the closed subpath
 * (joined at the start vertex, no caps).
 *
 * **SG24 Zero-length dashes.**  A dash with equal ends is a dot: round gives the
 * disc of radius `h`, square the square of side `2h` with two sides parallel to
 * the tangent there, butt nothing.  The tangent is the start direction (SG5) of
 * the segment that begins at the point (at a vertex the later segment's), the end
 * direction of the last segment at the end of the subpath, SG11's effective
 * tangent at a degenerate subpath.  A dot is drawn whatever `degenerate` is; a
 * degenerate subpath is itself omitted when `degenerate` says omit.
 *
 * **SG25 Overlap and order.**  Dashes are one union (SG13).  Pieces are emitted
 * in this order: subpaths in order; in a subpath its dashes in increasing
 * position; in a dash the start cap, then for each span (the part of a segment
 * inside it) the body pieces and the join at the vertex that ends the span, then
 * the end cap; a dot in its place.
 *
 * **SG26 The cut position's error.**  A dash end at distance `p` is placed within
 * `eps` of the point at that distance, and the call is ::GSHGEOM_ERR_LIMIT unless
 * `eps * norm <= 2^-20` pixel, `norm` the map's Frobenius norm under local width
 * and 1 under device width, `eps = 2^-46 S + (rows + 4) 2^-53 L_total + 2^-52
 * L_total + (n + 8) 2^-53 (L_total + 2P)`, `S` the total length of the control
 * polygons, `L_total` the path's length, `rows` along.h's table rows and `n` the
 * canonical count.  The first three terms are along.h's contract (an estimate:
 * the tests measure the real error against an independent quadrature), the last
 * SG29's arithmetic.  Moving a dash end by `delta` along the centerline moves
 * every boundary point of the piece by at most `(1 + h kappa) delta`, `kappa` the
 * curvature there: the outline error from the cut is at most `(1 + h kappa)
 * 2^-20` pixel, within the reserve SG17 keeps (2^-8 pixel) where `h kappa <=
 * 2^12`.
 *
 * **SG27 Alignment.**  The aligned region is the dashed region of the closed
 * subpaths at width `2w`, the body the dashed region of the others at `w`; the
 * clip region is of the undashed closed subpaths.
 *
 * **SG28 Bounds.**  A dash piece lies within the box of the undashed stroke.
 *
 * **SG29 Limits and arithmetic.**  A dash boundary is `(q P + S) - o` for the
 * period count `q`, a prefix sum `S` and the phase `o`, reduced by an exact long
 * division, so that its error does not grow with the dashes before it.  The
 * entries visited by one call are at most the work cap (`4 max_points + 4096`);
 * more is ::GSHGEOM_ERR_LIMIT.  A dash shorter than the spacing of doubles at its
 * position is a dot, never dropped.
 *
 * **SG30 A DXF simple linetype.**  A pattern of elements (DXF group 49: positive,
 * a dash with the pen down; negative, a space; zero, a dot), repeated from the
 * start of the line and scaled, lowers to a dash array and an offset by
 * ::gshgeom_stroke_dashes_from_dxf: consecutive pen-down elements are one run
 * (a dot adding zero), consecutive spaces one run; the first and last runs are
 * one run when they are of one kind (the pattern repeats, and the offset places
 * the line's start inside it); the runs are rotated to start with a pen-down
 * run.  A dot first with a dash last keeps the dot apart (a zero-length space
 * between).  No pen-down element draws nothing (::GSHGEOM_ERR_INVALID), no space
 * is solid, and a pattern that begins with a space which no trailing space
 * continues and ends with a dot cannot be expressed (::GSHGEOM_ERR_UNSUPPORTED:
 * the line would start exactly at that dot, which SG21 keeps).
 *
 * ## The hash and the dump (format 1)
 *
 * The hash, in the spelling of hash.h: `u32` 1, `u32` kind `0x4b525453` ("STRK"
 * little-endian), `u8` status, `u8` alignment, `u8` clip rule, `u8` complement,
 * then `u64` the ::gshgeom_flat_path_hash of body, aligned and clip (the offset
 * basis for an absent part).  The dump: a line
 * `gshgeom.stroke-outline 1 status=N align=N rule=N complement=N hash=H`, then
 * for each part a line `part body`, `part aligned` or `part clip` followed by
 * its flattened path's dump, or the line `part <name> -` when it is absent.
 */

#ifndef GHOTI_IO_GSHGEOM_STROKE_H
#define GHOTI_IO_GSHGEOM_STROKE_H

#include <ghoti.io/shape-geometry/allocator.h>
#include <ghoti.io/shape-geometry/core.h>
#include <ghoti.io/shape-geometry/fill.h>
#include <ghoti.io/shape-geometry/flatten.h>
#include <ghoti.io/shape-geometry/macros.h>
#include <ghoti.io/shape-geometry/map.h>
#include <ghoti.io/shape-geometry/path.h>

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/** The width space (SG10).  The values are pinned. */
typedef enum GSHGEOM_Stroke_Width_Space {
  GSHGEOM_STROKE_WIDTH_LOCAL = 0,  ///< The map applies to the width (an elliptical pen under a non-uniform map).
  GSHGEOM_STROKE_WIDTH_DEVICE = 1, ///< The width is in device pixels around the mapped centerline.
  GSHGEOM_STROKE_WIDTH_SPACE_COUNT ///< Closes the enum.
} GSHGEOM_Stroke_Width_Space;

/** The join (SG8).  The values are pinned (PDF's `j` for the first three). */
typedef enum GSHGEOM_Stroke_Join {
  GSHGEOM_STROKE_JOIN_MITER = 0,      ///< Miter, a bevel past the limit.
  GSHGEOM_STROKE_JOIN_ROUND = 1,      ///< A circular sector.
  GSHGEOM_STROKE_JOIN_BEVEL = 2,      ///< The triangle.
  GSHGEOM_STROKE_JOIN_MITER_CLIP = 3, ///< Miter, clipped at the limit.
  GSHGEOM_STROKE_JOIN_ARCS = 4,       ///< SVG 2's arcs (the 2018 Candidate Recommendation).
  GSHGEOM_STROKE_JOIN_COUNT           ///< Closes the enum.
} GSHGEOM_Stroke_Join;

/** The cap (SG7).  The values are pinned (PDF's `J`). */
typedef enum GSHGEOM_Stroke_Cap {
  GSHGEOM_STROKE_CAP_BUTT = 0,   ///< Nothing.
  GSHGEOM_STROKE_CAP_ROUND = 1,  ///< A half disc.
  GSHGEOM_STROKE_CAP_SQUARE = 2, ///< A half square.
  GSHGEOM_STROKE_CAP_COUNT       ///< Closes the enum.
} GSHGEOM_Stroke_Cap;

/** What a degenerate subpath draws (SG11).  The values are pinned. */
typedef enum GSHGEOM_Stroke_Degenerate {
  GSHGEOM_STROKE_DEGENERATE_DRAW = 0, ///< The cap's shape, along the effective tangent.
  GSHGEOM_STROKE_DEGENERATE_OMIT = 1, ///< Nothing.
  GSHGEOM_STROKE_DEGENERATE_COUNT     ///< Closes the enum.
} GSHGEOM_Stroke_Degenerate;

/** Where a dash pattern restarts.  The values are pinned. */
typedef enum GSHGEOM_Stroke_Dash_Scope {
  GSHGEOM_STROKE_DASH_PER_SUBPATH = 0, ///< At each subpath (SVG 2).
  GSHGEOM_STROKE_DASH_PATH = 1,        ///< Once for the path.
  GSHGEOM_STROKE_DASH_SCOPE_COUNT      ///< Closes the enum.
} GSHGEOM_Stroke_Dash_Scope;

/** The alignment (SG12).  The values are pinned. */
typedef enum GSHGEOM_Stroke_Alignment {
  GSHGEOM_STROKE_ALIGN_CENTER = 0,  ///< Centered on the centerline.
  GSHGEOM_STROKE_ALIGN_INSIDE = 1,  ///< Twice the width, inside the fill region.
  GSHGEOM_STROKE_ALIGN_OUTSIDE = 2, ///< Twice the width, outside the fill region.
  GSHGEOM_STROKE_ALIGN_COUNT        ///< Closes the enum.
} GSHGEOM_Stroke_Alignment;

/** The smallest magnitude of a nonzero width, dash or dash offset (SG2): 2^-64. */
#define GSHGEOM_STROKE_VALUE_MIN 0x1p-64
/** The largest magnitude of a width, dash or dash offset (SG2): 2^64. */
#define GSHGEOM_STROKE_VALUE_MAX 0x1p64
/** The largest miter limit (SG2): 2^16.  The smallest is 1. */
#define GSHGEOM_STROKE_MITER_LIMIT_MAX 0x1p16
/** The default miter limit (SG1), and SG3's value where the limit does not apply. */
#define GSHGEOM_STROKE_MITER_LIMIT_DEFAULT 4.0
/** The most dashes after an odd array is repeated (SG2). */
#define GSHGEOM_STROKE_DASHES_MAX ((size_t)65535)
/** The revision of the expansion's construction (SG18). */
#define GSHGEOM_STROKE_REVISION 1
/** The most dash steps one call visits when the caller states no bound (::gshgeom_stroke_dash_each). */
#define GSHGEOM_STROKE_DASH_STEPS_DEFAULT ((size_t)1 << 24)
/** How many times a curve segment may be halved before a piece is refused (::GSHGEOM_ERR_LIMIT). */
#define GSHGEOM_STROKE_MAX_DEPTH 64

/** A stroke's parameters (SG1).  Plain data; start from ::gshgeom_stroke_init. */
typedef struct GSHGEOM_Stroke {
  double width;                     ///< 0 or of magnitude in [2^-64, 2^64]; 0 draws nothing (SG14).
  uint32_t width_space;             ///< A ::GSHGEOM_Stroke_Width_Space.
  uint32_t join;                    ///< A ::GSHGEOM_Stroke_Join.
  double miter_limit;               ///< In [1, 2^16]; read by the miter, miter-clip and arcs joins.
  uint32_t cap;                     ///< A ::GSHGEOM_Stroke_Cap.
  uint32_t degenerate;              ///< A ::GSHGEOM_Stroke_Degenerate.
  size_t dash_count;                ///< How many dash lengths; 0 is solid.
  const double * dashes;            ///< The dash lengths; NULL only when @p dash_count is 0.  Borrowed.
  double dash_offset;               ///< Read only with dashes.
  uint32_t dash_scope;              ///< A ::GSHGEOM_Stroke_Dash_Scope; read only with dashes.
  uint32_t alignment;               ///< A ::GSHGEOM_Stroke_Alignment.
  GSHGEOM_Fill_Rule alignment_rule; ///< The rule of `F` (SG12); read only off center.
  const void * width_profile;       ///< Reserved: NULL; anything else is ::GSHGEOM_ERR_UNSUPPORTED.
  uint8_t _reserved[8];             ///< Zero.
} GSHGEOM_Stroke;

/** The status of an outline (SG14 to SG16).  The values are pinned (they are hashed). */
typedef enum GSHGEOM_Stroke_Status {
  GSHGEOM_STROKE_DRAWN = 0,          ///< The parts are as SG16 states.
  GSHGEOM_STROKE_EMPTY = 1,          ///< Width 0: nothing is covered (SG14).
  GSHGEOM_STROKE_NOT_INVERTIBLE = 2, ///< The map's linear part is singular: nothing is drawn (SG15).
  GSHGEOM_STROKE_STATUS_COUNT        ///< Closes the enum.
} GSHGEOM_Stroke_Status;

/** The expansion's result (SG16).  Opaque. */
typedef struct GSHGEOM_Stroke_Outline GSHGEOM_Stroke_Outline;

/** @brief Write SG1's defaults into @p stroke; NULL is ignored. */
GSHGEOM_API void gshgeom_stroke_init(GSHGEOM_Stroke * stroke);

/**
 * @brief Check the parameters (SG2).
 *
 * @return ::GSHGEOM_OK; ::GSHGEOM_ERR_INVALID, ::GSHGEOM_ERR_UNSUPPORTED or
 *   ::GSHGEOM_ERR_LIMIT, decided in SG2's order.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_stroke_check(const GSHGEOM_Stroke * stroke);

/**
 * @brief SG15's predicate: whether a 2D map's linear part has a nonzero
 *        determinant, decided exactly.
 *
 * @param invertible Receives 1 or 0.
 * @return ::GSHGEOM_OK; ::GSHGEOM_ERR_INVALID for NULL, an invalid map or a 3D map.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_stroke_map_invertible(
    const GSHGEOM_Map * map, int * invertible);

/**
 * @brief Expand a stroke at the default tolerance (SG16).
 *
 * Exactly ::gshgeom_stroke_expand_with_tolerance at ::GSHGEOM_FLATTEN_TOLERANCE.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_stroke_expand(const GSHGEOM_Path * centerline,
    const GSHGEOM_Stroke * stroke, const GSHGEOM_Map * map,
    const GSHGEOM_Allocator * alloc, const GSHGEOM_Limits * limits,
    GSHGEOM_Stroke_Outline ** out);

/**
 * @brief Expand a stroke into its outline at a tolerance (SG16, SG17).
 *
 * @param centerline A 2D model path.
 * @param stroke     The parameters (SG2 checks them).
 * @param map        A 2D map from model space to device pixels.
 * @param tolerance  `tau`, in 26.6 units; at least ::GSHGEOM_FLATTEN_MIN_TOLERANCE.
 * @param alloc      The allocator; NULL for the default.
 * @param limits     `max_points` bounds each part; NULL for the defaults.
 * @param out        Receives the outline, only on success; free it with
 *   ::gshgeom_stroke_outline_destroy.
 * @return ::GSHGEOM_OK (a not-invertible map and width 0 included);
 *   ::GSHGEOM_ERR_INVALID for NULL, an invalid map, a 3D map or a tolerance
 *   below the minimum, and as SG2; ::GSHGEOM_ERR_UNSUPPORTED for a 3D
 *   centerline or a width profile; ::GSHGEOM_ERR_LIMIT
 *   as SG2, for an outline vertex outside the 26.6 range, a part past
 *   `max_points`, a magnitude past the expansion's arithmetic bound
 *   (documentation/design.md, "The expansion's arithmetic"), a piece the
 *   depth cap would force, more dash steps than the work cap (`4 max_points +
 *   4096`) or a cut position's error past SG26's gate; ::GSHGEOM_ERR_OOM.  Nothing is returned on failure.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_stroke_expand_with_tolerance(
    const GSHGEOM_Path * centerline, const GSHGEOM_Stroke * stroke,
    const GSHGEOM_Map * map, int32_t tolerance, const GSHGEOM_Allocator * alloc,
    const GSHGEOM_Limits * limits, GSHGEOM_Stroke_Outline ** out);

/**
 * @brief Receives one dash of ::gshgeom_stroke_dash_each: the distances of its two ends.
 *
 * @return ::GSHGEOM_OK to go on; anything else stops the enumeration and is returned by it.
 */
typedef GSHGEOM_Result (*GSHGEOM_Dash_Fn)(void * ctx, double start, double end);

/**
 * @brief SG21: the dash positions of a subpath of @p length, in increasing order.
 *
 * SG3 is applied to the raw array: a negative length or all zeros is solid and gives the one pair `<0, length>`,
 * and an odd count is repeated once.  The phase is @p offset reduced by the floored modulo.
 *
 * @param max_steps At most this many entries of the pattern are visited after the first piece
 *   (::GSHGEOM_STROKE_DASH_STEPS_DEFAULT is a bound that cannot be reached in a useful time).
 * @return ::GSHGEOM_OK; ::GSHGEOM_ERR_INVALID for NULL, a non-finite value or a negative length;
 *   ::GSHGEOM_ERR_LIMIT for SG2's bounds, more entries than @p max_steps, or what @p fn returns.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_stroke_dash_each(const double * dashes, size_t count, double offset, double length,
    size_t max_steps, GSHGEOM_Dash_Fn fn, void * ctx);

/**
 * @brief SG30: lower a DXF simple linetype to a dash array and an offset.
 *
 * @param elements  DXF group 49 values: positive a dash, negative a space, zero a dot.
 * @param scale     A positive factor applied to every element.
 * @param dashes    Receives the array (at most @p capacity entries; a pattern of `m` runs needs `m + 1`).
 * @param dash_count Receives the count; 0 means solid.
 * @param offset    Receives the offset (SG21's phase at the start of the line).
 * @return ::GSHGEOM_OK; ::GSHGEOM_ERR_INVALID (NULL, a non-finite value, no pen-down element, a zero total);
 *   ::GSHGEOM_ERR_LIMIT (SG2's bounds, the capacity); ::GSHGEOM_ERR_UNSUPPORTED (a space first and a dot last).
 */
GSHGEOM_API GSHGEOM_Result gshgeom_stroke_dashes_from_dxf(const double * elements, size_t count, double scale,
    double * dashes, size_t capacity, size_t * dash_count, double * offset);

/** @brief Free an outline; NULL is ignored. */
GSHGEOM_API void gshgeom_stroke_outline_destroy(GSHGEOM_Stroke_Outline * outline);

/** @brief The outline's status; ::GSHGEOM_STROKE_STATUS_COUNT for NULL. */
GSHGEOM_API GSHGEOM_Stroke_Status gshgeom_stroke_outline_status(
    const GSHGEOM_Stroke_Outline * outline);

/** @brief The body part (SG16), or NULL when absent; owned by the outline. */
GSHGEOM_API const GSHGEOM_Flat_Path * gshgeom_stroke_outline_body(
    const GSHGEOM_Stroke_Outline * outline);

/** @brief The aligned part (SG16), or NULL when absent; owned by the outline. */
GSHGEOM_API const GSHGEOM_Flat_Path * gshgeom_stroke_outline_aligned(
    const GSHGEOM_Stroke_Outline * outline);

/** @brief The clip part (SG16), or NULL when absent; owned by the outline. */
GSHGEOM_API const GSHGEOM_Flat_Path * gshgeom_stroke_outline_clip(
    const GSHGEOM_Stroke_Outline * outline);

/** @brief The clip part's fill rule (the canonical alignment rule). */
GSHGEOM_API GSHGEOM_Fill_Rule gshgeom_stroke_outline_clip_rule(
    const GSHGEOM_Stroke_Outline * outline);

/** @brief 1 when the aligned part is kept outside the clip (outside alignment), else 0. */
GSHGEOM_API int gshgeom_stroke_outline_clip_complement(
    const GSHGEOM_Stroke_Outline * outline);

/** @brief The outline's hash (format 1, this file's comment); the offset basis for NULL. */
GSHGEOM_API uint64_t gshgeom_stroke_outline_hash(const GSHGEOM_Stroke_Outline * outline);

/**
 * @brief Write the outline's dump (format 1, this file's comment).
 *
 * @return ::GSHGEOM_OK, ::GSHGEOM_ERR_INVALID for NULL, ::GSHGEOM_ERR_IO when
 *   the stream fails (including at the flush).
 */
GSHGEOM_API GSHGEOM_Result gshgeom_stroke_outline_dump(
    const GSHGEOM_Stroke_Outline * outline, FILE * out);

/** @brief The construction's revision (::GSHGEOM_STROKE_REVISION). */
GSHGEOM_API int gshgeom_stroke_revision(void);

#ifdef __cplusplus
}
#endif

#endif /* GHOTI_IO_GSHGEOM_STROKE_H */
