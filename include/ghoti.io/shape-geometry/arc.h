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
 * @file arc.h
 *
 * Exact arcs: an elliptical arc as one element of a 2D model path (::GSHGEOM_VERB_ARC).
 *
 * ### The element
 *
 * An arc keeps its start and end exactly as the caller gave them, and its centre form: a
 * centre `c`, two conjugate half-diameters `u` and `v`, a start parameter `t0`, a signed
 * `sweep` and the unit the two angles are in.  The curve is
 * `P(t) = c + u cos t + v sin t` for `t` from `t0` to `t0 + sweep`.  Conjugate half-diameters
 * are closed under every affine map (the map is applied to the three vectors), which a centre,
 * two radii and a rotation is not.  The radii and rotation of the SVG and Canvas forms are
 * derived, never stored.
 *
 * The ends are authoritative: nothing recomputes them.  The centre form must reproduce them
 * to within ::GSHGEOM_ARC_AGREE times the arc's scale `s`, the largest magnitude of any
 * coordinate of `start`, `end` and `c` and of `|u|_inf + |v|_inf`.  A sweep of exactly one
 * turn requires `end` equal to `start` bit for bit.
 *
 * A path that holds an arc is read back with ::gshgeom_path_verb and ::gshgeom_path_arc.
 * ::gshgeom_path_element on an arc is ::GSHGEOM_ERR_UNSUPPORTED, so a walker that does not
 * know arcs fails instead of drawing the chord.  ::gshgeom_path_lower_arcs gives the path
 * with every arc replaced by the rational quadratics of ::gshgeom_arc_pieces.
 *
 * ### Orientation
 *
 * `t` increases from `t0`.  Every constructor builds `cross(u, v) > 0`, so a positive sweep
 * turns from +x toward +y (SVG's sweep-flag 1).  `sin` and `cos` of a parameter come from
 * ::gshgeom_sincos, so a quarter turn in degrees, gradians or turns is exact.
 *
 * ### The decomposition
 *
 * ::gshgeom_arc_pieces gives at most four rational quadratics of at most a quarter turn each,
 * weights `(1, w, 1)` with `w = cos h` for the half sweep `h` of the piece, exactly
 * `gmath_sqrt(0.5)` for a quarter turn.  The control point is the meeting of the tangents at
 * the piece's ends; where both tangents are parallel to an axis (an exact zero component) it
 * takes each coordinate from the end point on that tangent, so a rounded rectangle's corner
 * has exact control points.  The first piece starts at `start` and the last ends at `end`.  The
 * pieces lie on the ellipse of the centre form to within ::GSHGEOM_ARC_DECOMP times the scale.
 * A boundary between two pieces is not a vertex.
 *
 * ### The hash (version 1) and the dump
 *
 * After the `u8` verb 7: `u8` flags, `u8` unit, then `f64` for `start.x`, `start.y`, `end.x`,
 * `end.y`, `center.x`, `center.y`, `u.x`, `u.y`, `v.x`, `v.y`, `t0`, `sweep`.  The dump line is
 * `A <flags> <unit>` and the twelve numbers, each the 16 hexadecimal digits of the bit pattern.
 */

#ifndef GHOTI_IO_GSHGEOM_ARC_H
#define GHOTI_IO_GSHGEOM_ARC_H

#include <ghoti.io/shape-geometry/allocator.h>
#include <ghoti.io/shape-geometry/core.h>
#include <ghoti.io/shape-geometry/curve.h>
#include <ghoti.io/shape-geometry/macros.h>
#include <ghoti.io/shape-geometry/map.h>
#include <ghoti.io/shape-geometry/path.h>

#ifdef __cplusplus
extern "C" {
#endif

/** The arc's radii are equal and nothing since has stretched them unevenly (set by a constructor). */
#define GSHGEOM_ARC_CIRCULAR 1u

/**
 * An upper limit on the distance, relative to the arc's scale, by which the centre form may miss
 * the stored start and end: `2^-40`.
 */
#define GSHGEOM_ARC_AGREE 0x1p-40

/** An upper limit on the distance, relative to the arc's scale, of the pieces from the ellipse: `2^-44`. */
#define GSHGEOM_ARC_DECOMP 0x1p-44

/** The most pieces an arc decomposes into. */
#define GSHGEOM_ARC_MAX_PIECES 4

/** The centre form of an arc: everything but its ends. */
typedef struct GSHGEOM_Arc_Form {
  GSHGEOM_Point center;    ///< The centre `c`.
  GSHGEOM_Point u;         ///< A half-diameter.
  GSHGEOM_Point v;         ///< The conjugate half-diameter; `cross(u, v)` must not be zero.
  GSHGEOM_Angle_Unit unit; ///< The unit of `t0` and `sweep`.
  double t0;               ///< The start parameter.
  double sweep;            ///< The signed sweep: not zero, at most one turn in magnitude.
  unsigned flags;          ///< ::GSHGEOM_ARC_CIRCULAR or 0.
} GSHGEOM_Arc_Form;

/** One arc element: its ends as given and its centre form. */
typedef struct GSHGEOM_Arc {
  GSHGEOM_Point start;     ///< The start, equal to the current point bit for bit; `z` is 0.
  GSHGEOM_Point end;      ///< The end.
  GSHGEOM_Point center;   ///< The centre `c`.
  GSHGEOM_Point u;        ///< A half-diameter.
  GSHGEOM_Point v;        ///< The conjugate half-diameter.
  GSHGEOM_Angle_Unit unit; ///< The unit of `t0` and `sweep`.
  double t0;              ///< The start parameter.
  double sweep;           ///< The signed sweep.
  unsigned flags;         ///< ::GSHGEOM_ARC_CIRCULAR or 0.
} GSHGEOM_Arc;

/** What an endpoint-form call added. */
typedef enum GSHGEOM_Arc_Kind {
  GSHGEOM_ARC_KIND_NOTHING = 0, ///< The end was the current point: nothing was added.
  GSHGEOM_ARC_KIND_LINE = 1,    ///< A radius was zero: a line was added.
  GSHGEOM_ARC_KIND_ARC = 2      ///< An arc was added.
} GSHGEOM_Arc_Kind;

/**
 * @brief The point `P(t)` of a centre form, as the constructors compute it.
 *
 * `c + (u cos t + v sin t)`, the inner sum first, each operation rounded, `cos` and `sin` from
 * ::gshgeom_sincos in the form's unit.
 *
 * @return ::GSHGEOM_OK, or ::GSHGEOM_ERR_INVALID (NULL, non-finite, unit) or ::GSHGEOM_ERR_LIMIT.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_arc_point(
    const GSHGEOM_Arc_Form * form, double t, GSHGEOM_Point * out);

/**
 * @brief Append an arc record raw.
 *
 * Checked in this order, a failed call leaving the path as it was: ::GSHGEOM_ERR_INVALID for a NULL or 3D
 * path, no current point, a non-finite number, an unknown unit or flag bit, `cross(u, v)` zero (decided
 * exactly), a sweep of zero or past one turn; then ::GSHGEOM_ERR_LIMIT for a coordinate or a component of
 * `u` or `v` past ::GSHGEOM_COORD_LIMIT, `|t0|` past ::GSHGEOM_MAP_MAX_QUARTERS quarter turns, or the
 * element limit; then ::GSHGEOM_ERR_INVALID if `start` is not the current point bit for bit or the centre
 * form misses an end by more than ::GSHGEOM_ARC_AGREE (or a full turn does not close); then
 * ::GSHGEOM_ERR_OOM.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_path_append_arc(GSHGEOM_Path * path, const GSHGEOM_Arc * arc);

/** @brief Read arc element @p index back; ::GSHGEOM_ERR_INVALID if it is not an arc or out of range. */
GSHGEOM_API GSHGEOM_Result gshgeom_path_arc(const GSHGEOM_Path * path, size_t index, GSHGEOM_Arc * out);

/**
 * @brief The rational quadratics an arc decomposes into.
 *
 * @param arc The arc (not checked beyond what the pieces need: finite, unit valid).
 * @param out Receives @p n curves, dimension 2.
 * @param n Receives the count, 1 to ::GSHGEOM_ARC_MAX_PIECES: 1 to 4 for a sweep up to one, two, three, four quarter turns.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_arc_pieces(
    const GSHGEOM_Arc * arc, GSHGEOM_Curve out[GSHGEOM_ARC_MAX_PIECES], int * n);

/**
 * @brief The path with every arc replaced by its pieces, every other element as it was.
 *
 * @param path The path.
 * @param alloc The allocator for the result, or NULL for the default.
 * @param limits The limits for the result, or NULL for the defaults.
 * @param out The new path; written only on success.
 * @return ::GSHGEOM_OK, ::GSHGEOM_ERR_INVALID, ::GSHGEOM_ERR_LIMIT or ::GSHGEOM_ERR_OOM.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_path_lower_arcs(const GSHGEOM_Path * path,
    const GSHGEOM_Allocator * alloc, const GSHGEOM_Limits * limits, GSHGEOM_Path ** out);

/**
 * @brief Map an arc: an affine map acts linearly on `P(t)`, so the ends and the centre move as points, the
 * half-diameters as vectors, and the parameters stay.
 *
 * `start`, `end` and `center` go through ::gshgeom_map_apply_point and `u`, `v` through
 * ::gshgeom_map_apply_vector, and nothing else is multiplied.  `unit`, `t0` and `sweep` are kept; a mirror
 * flips the sign of `cross(u, v)` and with it the arc's direction.  The result is checked as
 * ::gshgeom_path_append_arc checks a record (the ends are still the map's images of the ends, and the
 * centre form still reproduces them).  ::GSHGEOM_ARC_CIRCULAR is kept if and only if the map's kind tag has
 * none of ::GSHGEOM_MAP_SCALES_NONUNIFORM, ::GSHGEOM_MAP_SKEWS, ::GSHGEOM_MAP_SINGULAR and ::GSHGEOM_MAP_GENERAL
 * (the similarities), and is cleared otherwise; a rotation by an arbitrary angle makes `u . v` a few ulps from
 * zero, so the tag, not a test of the numbers, is the predicate.
 *
 * @param map A 2D map.
 * @param arc The arc.
 * @param out The mapped arc; written only on success.
 * @return ::GSHGEOM_OK; ::GSHGEOM_ERR_INVALID for a NULL, a map that is not 2D or invalid, a singular map (its
 *   image is a segment; a consumer that needs one maps the arc's pieces, ::gshgeom_arc_pieces, instead) or an
 *   arc that is not valid; ::GSHGEOM_ERR_LIMIT if a result is past ::GSHGEOM_COORD_LIMIT.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_arc_map(const GSHGEOM_Map * map, const GSHGEOM_Arc * arc, GSHGEOM_Arc * out);

/**
 * @brief The centre, radius and sense of a circular arc, for a toolpath or a writer.
 *
 * Refuses (::GSHGEOM_ERR_INVALID) an arc that does not carry ::GSHGEOM_ARC_CIRCULAR, never inferring a circle
 * from the numbers.  The radius is `|u|` (::gmath_hypot).
 *
 * @param centre Receives the centre.
 * @param radius Receives the radius.
 * @param ccw Receives nonzero when the arc runs in the positive sense (`sign(sweep) * sign(cross(u, v))` is +1,
 *   from +x toward +y), zero when it runs the other way.  May be NULL.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_arc_circle(
    const GSHGEOM_Arc * arc, GSHGEOM_Point * centre, double * radius, int * ccw);

/**
 * @brief Endpoint form (SVG 2's endpoint-to-centre conversion and radius correction).
 *
 * From the current point `S` to @p end: the same point adds nothing (`NOTHING`); a zero radius adds a
 * line (`LINE`); otherwise the radii are made absolute, scaled up if too small (the half ellipse,
 * decided exactly), and the arc the two flags name is added (`ARC`).  The rotation is applied to the
 * radii as `u = rx (cos phi, sin phi)`, `v = ry (-sin phi, cos phi)`.  The `CIRCULAR` flag is set iff
 * the radii are bit equal after the correction.  `S` and @p end are kept as given.
 *
 * @param large 0 or 1 (else ::GSHGEOM_ERR_INVALID).
 * @param sweep 0 or 1: 1 is a positive sweep.
 * @param kind Receives what was added; may be NULL.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_path_arc_endpoint_to(GSHGEOM_Path * path, GSHGEOM_Point end,
    double rx, double ry, double phi, GSHGEOM_Angle_Unit unit, int large, int sweep,
    GSHGEOM_Arc_Kind * kind);

/**
 * @brief Centre form from its fields (DXF `ARC`, `ELLIPSE`; Canvas `ellipse` after resolution).
 *
 * The start is the current point if there is one, else a subpath begins at `P(t0)`.  @p end overrides
 * `P(t0 + sweep)` when not NULL; a sweep of exactly one turn closes (`end = start`).  The call is
 * atomic.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_path_arc_conjugate_to(
    GSHGEOM_Path * path, const GSHGEOM_Arc_Form * form, const GSHGEOM_Point * end);

/**
 * @brief Centre form from radii and rotation (Canvas `arc`, DXF `CIRCLE`): `u` and `v` as for the endpoint form.
 *
 * @param rx,ry Positive radii; the `CIRCULAR` flag is set iff bit equal.
 * @param phi The rotation, in @p unit as are @p t0 and @p sweep.
 * @param end An explicit end, or NULL.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_path_arc_center_to(GSHGEOM_Path * path, GSHGEOM_Point center,
    double rx, double ry, double phi, GSHGEOM_Angle_Unit unit, double t0, double sweep,
    const GSHGEOM_Point * end);

/**
 * @brief Canvas `arcTo`: a line to the first tangent point and a circular arc to the second.
 *
 * The result is a line to @p p1 when the current point equals @p p1 or @p p1 equals @p p2 (bit
 * equality), @p r is zero, or the three points are collinear (exactly).  A negative or non-finite
 * @p r is ::GSHGEOM_ERR_INVALID.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_path_arc_tangent_to(
    GSHGEOM_Path * path, GSHGEOM_Point p1, GSHGEOM_Point p2, double r);

/**
 * @brief A DXF polyline bulge from the current point to @p q: `bulge = tan(sweep / 4)`.
 *
 * Zero is a line; @p q equal to the current point adds nothing; NaN or infinity is ::GSHGEOM_ERR_INVALID.
 * The arc is circular, its sweep `4 atan(bulge)` in radians.
 */
GSHGEOM_API GSHGEOM_Result gshgeom_path_arc_bulge_to(
    GSHGEOM_Path * path, GSHGEOM_Point q, double bulge);

#ifdef __cplusplus
}
#endif

#endif /* GHOTI_IO_GSHGEOM_ARC_H */
