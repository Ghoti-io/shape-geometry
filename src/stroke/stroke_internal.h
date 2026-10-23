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
 * @file
 *
 * What the files of the stroke expansion share (internal; stroke.h states the
 * rules, emit.c the guarantee and its proof).
 *
 * Geometry is built in the **build space**: model space under local width,
 * device pixels under device width (SG10).  A point of a piece is given as an
 * offset `d` from a build-space anchor `o`, and reaches device space as
 * `A + lin d`, `A` the anchor's device position (for local width
 * ::gshgeom_map_apply_point of `o`, as a fill's point is formed; for device
 * width `o` itself) and `lin` the map's linear part (the identity under device
 * width).  Every distance a piece's bound is stated in is device pixels.
 */

#ifndef GHOTI_IO_GSHGEOM_STROKE_INTERNAL_H
#define GHOTI_IO_GSHGEOM_STROKE_INTERNAL_H

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/shape-geometry/device.h>
#include <ghoti.io/shape-geometry/along.h>
#include <ghoti.io/shape-geometry/flatten.h>
#include <ghoti.io/shape-geometry/map.h>
#include <ghoti.io/shape-geometry/path.h>
#include <ghoti.io/shape-geometry/stroke.h>
#include <ghoti.io/shape-geometry/stroke_curves.h>

#include <stddef.h>
#include <stdint.h>

/** A point or vector of the plane, in build space or in device pixels. */
typedef struct GSHGEOM_SK_Vec {
  double x;
  double y;
} GSHGEOM_SK_Vec;

/** One segment as the joins and caps read it: its build-space control points and weights, its degree. */
typedef struct GSHGEOM_SK_Seg {
  int n;                 ///< Degree, 1 to 3.
  int rational;          ///< Nonzero for the rational kinds.
  GSHGEOM_SK_Vec p[4];   ///< Control points, build space.
  double w[4];           ///< Weights (1 for a polynomial kind).
  size_t elem;           ///< The index of the path element it was read from.
  GSHGEOM_SK_Vec model_p[4]; ///< The control points before the map (what the vertex class is decided on).
  int vtx_valid;         ///< Nonzero when `vtx_class` was decided from the model points.
  int vtx_class;         ///< The vertex at this segment's start, against the segment before it: 0 a turn, 1 smooth, 2 a reversal.
} GSHGEOM_SK_Seg;

/** One subpath as SG5 reads it (expand.c). */
typedef struct GSHGEOM_SK_Sub {
  size_t first_seg;          ///< Its first segment with extent in the segment array.
  size_t seg_count;          ///< How many it has.
  int closed;                ///< Ended by a close.
  int draws;                 ///< Has a drawing verb or a close.
  GSHGEOM_SK_Vec start;      ///< Its start, build space.
  GSHGEOM_Point start_model; ///< Its start, model space (for the clip's path).
  size_t first_el;           ///< Its first drawing element in the path.
  size_t end_el;             ///< One past its last element.
} GSHGEOM_SK_Sub;

/** The centerline as SG5 reads it: subpaths and, in element order, their segments with extent (build space). */
typedef struct GSHGEOM_SK_Reader {
  GSHGEOM_SK_Sub * subs;
  size_t sub_n, sub_cap;
  GSHGEOM_SK_Seg * segs;
  size_t seg_n, seg_cap;
} GSHGEOM_SK_Reader;

/** The 26.6 rounding's reach: sqrt(2)/2 of a unit, rounded up (in 26.6 units). */
#define GSHGEOM_SK_ROUND 0.70711
/** The double arithmetic's reach on any device coordinate, under SK_MAG_LIMIT (in 26.6 units; emit.c). */
#define GSHGEOM_SK_ETA 0.015625
/** An upper bound of the emitted octagon's circumradius over its inradius (emit.c, (c)). */
#define GSHGEOM_SK_OCT 1.0826
/** What the proof keeps back of the tolerance (in 26.6 units). */
#define GSHGEOM_SK_MARGIN 0.25
/** The largest magnitude, in device pixels, the expansion's arithmetic bound is proved for (emit.c, (a)). */
#define GSHGEOM_SK_MAG_LIMIT 0x1p30
/** The smallest nonzero difference between neighbouring control points the expansion's build-space arithmetic (cubes
 * and fourth powers of a length) stays inside a double's range for. */
#define GSHGEOM_SK_BUILD_MIN 0x1p-200

/** One element of a contour of the curve outline, in target units (curve mode, src/stroke_curves/contour.c). */
typedef struct GSHGEOM_SK_CElem {
  int deg;                ///< 1 a line, 2 a rational quadratic, 3 a cubic.
  GSHGEOM_SK_Vec p[4];    ///< Control points; the first is the start, the last the end.
  double w[3];            ///< The weights of a rational quadratic.
} GSHGEOM_SK_CElem;

/** One accepted band piece of a run (offset.c): the control points of its two offset cubics, in target units. */
typedef struct GSHGEOM_SK_Band {
  GSHGEOM_SK_Vec l[4]; ///< The offset on the left of the centerline.
  GSHGEOM_SK_Vec r[4]; ///< The offset on the right.
} GSHGEOM_SK_Band;

/**
 * The state of the curve outline's expansion (stroke_curves.h; curve mode).  When ::GSHGEOM_SK_Ctx::cv is not NULL
 * every piece is written as lines, cubics and rational quadratics to `part`, in target units, instead of as 26.6
 * polygons.  `budget` is what a piece may spend on approximation (SG35), `flip` is set when the build space's
 * positive orientation is the device space's negative (a mirroring map), so contours are written reversed.
 */
typedef struct GSHGEOM_SK_Curves {
  GSHGEOM_Path * part;         ///< The part being written (body or aligned).
  double t;                    ///< The caller's tolerance.
  double budget;               ///< t / 2: the approximation each piece may spend.
  int flip;                    ///< The map reverses orientation.
  GSHGEOM_SK_CElem * el;       ///< The contour being built.
  size_t el_n, el_cap;
  GSHGEOM_SK_Band * run;       ///< The band pieces of the segment being built, in order.
  size_t run_n, run_cap;
} GSHGEOM_SK_Curves;

/** The state of a dashed expansion (dash.c). */
typedef struct GSHGEOM_SK_Dash {
  int active;                 ///< Nonzero when the style has a non-solid canonical array.
  double * dashes;            ///< The canonical array (an odd array repeated), owned.
  size_t n;                   ///< Its count (even, at most ::GSHGEOM_STROKE_DASHES_MAX).
  double period;              ///< Its sum.
  double offset;              ///< The style's offset.
  int per_path;               ///< ::GSHGEOM_STROKE_DASH_PATH.
  GSHGEOM_Measured * m;///< The centerline measured in build space.
  GSHGEOM_Path * bpath;       ///< The mapped centerline (device width), owned; NULL under local width.
  size_t msub;                ///< The next measured subpath.
  size_t steps;               ///< Pattern entries visited so far.
  size_t pieces;              ///< Dashes handed out so far.
  double phase0;              ///< The reduced offset (per subpath).
} GSHGEOM_SK_Dash;

/** The state of one expansion. */
typedef struct GSHGEOM_SK_Ctx {
  GSHGEOM_Allocator alloc;       ///< The caller's allocator.
  GSHGEOM_Limits limits;         ///< The caller's limits (defaults filled in).
  int32_t tau;                   ///< The tolerance, 26.6 units.
  double budget;                 ///< (tau - margin) in device pixels: how far outside its region a piece may reach.
  double rho;                    ///< (ROUND + ETA) in device pixels.
  GSHGEOM_Stroke st;             ///< The canonical parameters.
  int local;                     ///< 1 under local width.
  GSHGEOM_Map map;               ///< The map (2D).
  double lin[2][2];              ///< The linear part from build space to device pixels.
  double norm;                   ///< An upper bound of lin's spectral norm (its Frobenius norm, rounded up).
  double h;                      ///< Half the width in build space, of the part being built.
  double h_dev;                  ///< An upper bound of h's reach in device pixels: h * norm.
  GSHGEOM_Flat_Path * part;      ///< The flattened path the part is written to.
  GSHGEOM_SK_Vec * scratch;      ///< Hull scratch.
  size_t scratch_cap;            ///< Its capacity.
  GSHGEOM_SK_Vec * arcbuf;       ///< The arc points of the piece being built (arc.c).
  size_t arcbuf_cap;             ///< Its capacity.
  size_t arc_limit;              ///< When nonzero, the most points one gshgeom_sk_arc call may produce (arc.c).
  size_t work;                   ///< Curve pieces visited so far.
  size_t work_cap;               ///< The most curve pieces one expansion visits (ERR_LIMIT past it).
  int depth_cap;                 ///< How deep a curve may be split (::GSHGEOM_STROKE_MAX_DEPTH but in tests).
  GSHGEOM_SK_Dash dash;          ///< Dashes (dash.c).
  double mag_limit;              ///< The largest magnitude, in device units, the arithmetic bound is proved for.
  GSHGEOM_SK_Curves * cv;        ///< Non-NULL in curve mode.
} GSHGEOM_SK_Ctx;

/** The outline (stroke.h, SG16). */
struct GSHGEOM_Stroke_Outline {
  GSHGEOM_Allocator alloc;      ///< What it was made with.
  GSHGEOM_Stroke_Status status; ///< SG14 to SG16.
  uint32_t alignment;           ///< The canonical alignment.
  GSHGEOM_Fill_Rule rule;       ///< The clip's rule (the canonical alignment rule).
  int complement;               ///< 1 for outside alignment.
  GSHGEOM_Flat_Path * body;     ///< NULL when absent.
  GSHGEOM_Flat_Path * aligned;  ///< NULL when absent.
  GSHGEOM_Flat_Path * clip;     ///< NULL when absent.
};

/** The curve outline (stroke_curves.h, SG32). */
struct GSHGEOM_Stroke_Curves {
  GSHGEOM_Allocator alloc;           ///< What it was made with.
  GSHGEOM_Stroke_Status status;      ///< SG14, SG15.
  GSHGEOM_Stroke_Curves_Form form;   ///< The form it was asked for.
  uint32_t alignment;                ///< The canonical alignment.
  GSHGEOM_Fill_Rule rule;            ///< The clip's rule.
  int complement;                    ///< 1 for outside alignment.
  GSHGEOM_Path * body;               ///< NULL when absent.
  GSHGEOM_Path * aligned;            ///< NULL when absent.
  GSHGEOM_Path * clip;               ///< NULL when absent.
};

/**
 * ::gshgeom_stroke_expand_with_tolerance with the depth cap at @p depth_cap (1 to
 * ::GSHGEOM_STROKE_MAX_DEPTH): the seam a test forces the cap through, since no
 * ordinary input reaches the real one.  Not exported; the tests link the archive.
 */
GSHGEOM_Result gshgeom_stroke_expand_capped(const GSHGEOM_Path * centerline, const GSHGEOM_Stroke * stroke,
    const GSHGEOM_Map * map, int32_t tolerance, const GSHGEOM_Allocator * alloc, const GSHGEOM_Limits * limits,
    int depth_cap, GSHGEOM_Stroke_Outline ** out);

/** ::gshgeom_stroke_curves_expand with the depth cap at @p depth_cap (1 to ::GSHGEOM_STROKE_MAX_DEPTH): the test seam. */
GSHGEOM_Result gshgeom_stroke_curves_expand_capped(const GSHGEOM_Path * centerline, const GSHGEOM_Stroke * stroke,
    const GSHGEOM_Map * map, double tolerance, GSHGEOM_Stroke_Curves_Form form, const GSHGEOM_Allocator * alloc,
    const GSHGEOM_Limits * limits, int depth_cap, GSHGEOM_Stroke_Curves ** out);

/* --- params.c --------------------------------------------------------------- */

/** SG2 then SG3: the canonical parameters (a non-solid array keeps the caller's entries; see params.c). */
GSHGEOM_Result gshgeom_stroke_canonical(const GSHGEOM_Stroke * in, GSHGEOM_Stroke * out);

/* --- dash.c ----------------------------------------------------------------- */

/** The sum of the first @p n canonical entries of a raw array (count entries, repeated), in order. */
double gshgeom_sk_dash_period(const double * dashes, size_t count, size_t n);

/** The phase reduced to [0, period) exactly (SG21, SG29). */
double gshgeom_sk_dash_phase(double offset, double period);

/** What a dash says of itself beyond its ends (SG21, SG23). */
#define GSHGEOM_SK_POS_BEGINS_BEFORE 1 /**< the first piece began before the path: the phase is inside its dash */
#define GSHGEOM_SK_POS_ENDS_AFTER 2    /**< the dash would have gone on past the end of the subpath */

/** Receives one dash of ::gshgeom_sk_dash_positions with its flags. */
typedef GSHGEOM_Result (*GSHGEOM_SK_PosFn)(void * ctx, double start, double end, int flags);

/** SG21 with a reduced phase @p o: the pairs of a subpath of @p length, handed to @p fn. */
GSHGEOM_Result gshgeom_sk_dash_positions(const double * dashes, size_t count, size_t n, double period, double o,
    double length, size_t * steps, size_t max_steps, GSHGEOM_SK_PosFn fn, void * ctx);

/** SG26's gate and the measured path of the dashes; the canonical array comes from @p stroke. */
GSHGEOM_Result gshgeom_sk_dash_begin(GSHGEOM_SK_Ctx * c, const GSHGEOM_Path * centerline,
    const GSHGEOM_SK_Reader * rd, const GSHGEOM_Stroke * stroke);

/** Free what dash_begin made. */
void gshgeom_sk_dash_end(GSHGEOM_SK_Ctx * c);

/** The pieces of subpath @p i of the reader, dashed (SG21 to SG25), in the context's part. */
GSHGEOM_Result gshgeom_sk_dash_sub(GSHGEOM_SK_Ctx * c, const GSHGEOM_SK_Reader * rd, size_t i);

/** The solid stroke of subpath @p i (expand.c): what a dash that covers the whole of a closed subpath is. */
GSHGEOM_Result gshgeom_sk_solid_sub(GSHGEOM_SK_Ctx * c, const GSHGEOM_SK_Reader * rd, size_t i);

/** SG11's effective tangent of the degenerate subpath @p i (expand.c). */
GSHGEOM_SK_Vec gshgeom_sk_effective_tangent(const GSHGEOM_SK_Ctx * c, const GSHGEOM_SK_Reader * rd, size_t i);

/* --- emit.c ----------------------------------------------------------------- */

/** Map a build-space offset to a device offset (lin d). */
GSHGEOM_SK_Vec gshgeom_sk_lin(const GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec d);

/** The device position of a build-space anchor (the map applied as to a fill's point, or the point itself). */
GSHGEOM_Result gshgeom_sk_anchor(const GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec o, GSHGEOM_SK_Vec * out);

/**
 * Whether a piece whose core polygon is within @p d_in of containing its region
 * and within @p d_out of it (device pixels) fits the tolerance after the
 * inflation and the rounding (emit.c, (c)).
 */
int gshgeom_sk_fits(const GSHGEOM_SK_Ctx * c, double d_in, double d_out);

/** The most a piece with @p d_in may reach outside its region and still fit: what is left for d_out. */
double gshgeom_sk_room(const GSHGEOM_SK_Ctx * c, double d_in);

/**
 * Emit one piece: the convex hull of the core points `A + lin d[i]`, inflated by
 * an octagon of inradius `d_in + rho`, rounded to 26.6 and written to the part as
 * positively wound closed subpaths (emit.c).  The caller has checked
 * ::gshgeom_sk_fits; a piece that does not fit is ::GSHGEOM_ERR_INTERNAL.
 *
 * @return ::GSHGEOM_OK, ::GSHGEOM_ERR_LIMIT (the 26.6 range, `max_points`),
 *   ::GSHGEOM_ERR_OOM, ::GSHGEOM_ERR_INTERNAL.
 */
GSHGEOM_Result gshgeom_sk_emit(GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec anchor, const GSHGEOM_SK_Vec * d, size_t n,
    double d_in, double d_out);

/* --- arc.c ------------------------------------------------------------------ */

/**
 * Append to the context's arc buffer, from index `*n`, the circumscribed polygon
 * of the arc of radius @p r about the build-space origin from unit @p p to unit
 * @p q, the short way (angle below pi): `r p`, then for each sub-arc the tangent
 * intersection and the sub-arc's end, the last being `r q`.  `*dev` receives an
 * upper bound of its distance outside the arc in device pixels.  Sub-arcs are
 * halved until that bound is at most @p allowed.
 *
 * @return ::GSHGEOM_OK; ::GSHGEOM_ERR_LIMIT past `max_points` or when 64 halvings
 *   do not reach @p allowed; ::GSHGEOM_ERR_OOM.
 */
GSHGEOM_Result gshgeom_sk_arc(GSHGEOM_SK_Ctx * c, double r, GSHGEOM_SK_Vec p, GSHGEOM_SK_Vec q, double allowed,
    size_t * n, double * dev);

/** A sector of radius h about the anchor from unit p to unit q (the short way), as one convex piece. */
GSHGEOM_Result gshgeom_sk_sector(GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec anchor_dev, GSHGEOM_SK_Vec p, GSHGEOM_SK_Vec q,
    double d_in, double extra_out);

/** The full disc of radius h about the anchor, as one convex piece. */
GSHGEOM_Result gshgeom_sk_disc(GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec anchor_dev, double d_in, double extra_out);

/** The most chords of the clip line gshgeom_sk_net_chords takes. */
#define GSHGEOM_SK_CHORDS_MAX 64

/** Receives one directed edge. */
typedef GSHGEOM_Result (*GSHGEOM_SK_EdgeFn)(void * ctx, GSHGEOM_SK_Vec from, GSHGEOM_SK_Vec to);

/**
 * The chords a clip line cut into a boundary (each a directed segment along the line through @p y with unit
 * direction @p t), reduced to the edges that remain once the stretches two chords cover in opposite directions have
 * cancelled: the ends of the chords cut the line into stretches, each stretch carries the sum of the chords over it,
 * and a stretch with a sum of n is passed n times to @p fn in that direction (a perturbation of two chords that
 * cancel would otherwise leave a sliver along their whole overlap).
 *
 * @return ::GSHGEOM_OK, what @p fn returns when it is not that, ::GSHGEOM_ERR_LIMIT for more than
 *   ::GSHGEOM_SK_CHORDS_MAX chords.
 */
GSHGEOM_Result gshgeom_sk_net_chords(const GSHGEOM_SK_Vec (*chord)[2], size_t k, GSHGEOM_SK_Vec y,
    GSHGEOM_SK_Vec t, GSHGEOM_SK_EdgeFn fn, void * ctx);

/* --- join.c ----------------------------------------------------------------- */



/** SG5's start direction and the index of the control point that defines it; -1 for no extent. */
int gshgeom_sk_start_dir(const GSHGEOM_SK_Seg * s, GSHGEOM_SK_Vec * dir);

/** SG5's end direction and the index of the control point that defines it; -1 for no extent. */
int gshgeom_sk_end_dir(const GSHGEOM_SK_Seg * s, GSHGEOM_SK_Vec * dir);

/**
 * Whether the vertex of three points is a turn (0), smooth (1: `orient2d` zero and `a . b > 0`) or a reversal (2),
 * decided exactly on the given points.  An invertible linear map keeps this class, so under device width it is decided
 * on the points before the map: the mapped doubles are rounded and a vertex smooth in the model can be a turn of a few
 * ulps there (SG8, SG10).
 */
int gshgeom_sk_vertex_class(GSHGEOM_SK_Vec q, GSHGEOM_SK_Vec v, GSHGEOM_SK_Vec s);

/** The join at the vertex between @p in and @p out (SG8, SG9). */
GSHGEOM_Result gshgeom_sk_join(GSHGEOM_SK_Ctx * c, const GSHGEOM_SK_Seg * in, const GSHGEOM_SK_Seg * out);

/** The cap at an end @p e with the raw outward direction @p d (SG7). */
GSHGEOM_Result gshgeom_sk_cap(GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec e, GSHGEOM_SK_Vec d);

/** The degenerate shape at @p p with the raw effective tangent @p d (SG11). */
GSHGEOM_Result gshgeom_sk_degenerate(GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec p, GSHGEOM_SK_Vec d);

/**
 * SG9: whether the miter at `V` between `a = V - Q` and `b = S - V` is within
 * the limit @p L, decided exactly.
 */
int gshgeom_sk_miter_within(GSHGEOM_SK_Vec q, GSHGEOM_SK_Vec v, GSHGEOM_SK_Vec s, double L);

/**
 * The arcs join with the outer side given (arc.c, SG8).  `*miter_clip` is set,
 * and nothing drawn, when both curvatures are zero away from a reversal: the
 * caller then draws miter-clip's join.
 */
GSHGEOM_Result gshgeom_sk_arcs_join(GSHGEOM_SK_Ctx * c, const GSHGEOM_SK_Seg * in, const GSHGEOM_SK_Seg * out,
    int outer_left, int reversal, int * miter_clip);

/* --- contour.c (curve mode) ------------------------------------------------- */

/** Start a contour. */
void gshgeom_skc_begin(GSHGEOM_SK_Ctx * c);
/** Append a line (a zero-length line is dropped). */
GSHGEOM_Result gshgeom_skc_line(GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec a, GSHGEOM_SK_Vec b);
/** Append a cubic. */
GSHGEOM_Result gshgeom_skc_cubic(GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec p0, GSHGEOM_SK_Vec p1, GSHGEOM_SK_Vec p2,
    GSHGEOM_SK_Vec p3);
/** Append a rational quadratic of weights (1, w, 1). */
GSHGEOM_Result gshgeom_skc_conic(GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec p0, GSHGEOM_SK_Vec p1, GSHGEOM_SK_Vec p2, double w);
/**
 * Write the contour as one closed subpath of the part.  The contour was built in the build space's positive
 * orientation; it is reversed when the map mirrors.  A contour of fewer than two elements is dropped.
 */
GSHGEOM_Result gshgeom_skc_finish(GSHGEOM_SK_Ctx * c);
/** Write a convex polygon already in positive device orientation as lines (the hull of a polygon piece). */
GSHGEOM_Result gshgeom_skc_ring(GSHGEOM_SK_Ctx * c, const GSHGEOM_SK_Vec * pts, size_t n);
/** The sector of radius h about the device point @p centre from unit @p p to unit @p q (the short way). */
GSHGEOM_Result gshgeom_skc_sector(GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec centre, GSHGEOM_SK_Vec p, GSHGEOM_SK_Vec q);
/** The half disc `{v : |v| <= h, v . d >= 0}` about @p centre, @p m the left normal of the unit @p d. */
GSHGEOM_Result gshgeom_skc_half_disc(GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec centre, GSHGEOM_SK_Vec m, GSHGEOM_SK_Vec d);
/** The disc of radius @p radius (build space) about @p centre. */
GSHGEOM_Result gshgeom_skc_disc(GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec centre, double radius);
/**
 * The inflated sector of a fan: the sector of radius h about @p centre from unit @p p to unit @p q (the short way),
 * grown by @p r (build space): the half planes of its sides moved out by @p r, the circle of radius `h + r`, and the
 * line `r` behind the apex.  It contains the sector's Minkowski sum with a disc of radius @p r and lies within
 * `sqrt(2) r` of the sector (SG37).
 */
GSHGEOM_Result gshgeom_skc_fan_sector(GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec centre, GSHGEOM_SK_Vec p, GSHGEOM_SK_Vec q,
    double r);
/** Write the run of band pieces as one contour and empty the run. */
GSHGEOM_Result gshgeom_skc_flush_run(GSHGEOM_SK_Ctx * c);
/** Add a band piece to the run. */
GSHGEOM_Result gshgeom_skc_add_band(GSHGEOM_SK_Ctx * c, const GSHGEOM_SK_Band * b);

/* --- body.c ----------------------------------------------------------------- */

/** The body of one segment with extent (SG6). */
GSHGEOM_Result gshgeom_sk_body(GSHGEOM_SK_Ctx * c, const GSHGEOM_SK_Seg * s);

/* --- vector helpers --------------------------------------------------------- */

static inline GSHGEOM_SK_Vec gshgeom_sk_v(double x, double y) {
  GSHGEOM_SK_Vec v;
  v.x = x;
  v.y = y;
  return v;
}

static inline GSHGEOM_SK_Vec gshgeom_sk_add(GSHGEOM_SK_Vec a, GSHGEOM_SK_Vec b) {
  return gshgeom_sk_v(a.x + b.x, a.y + b.y);
}

static inline GSHGEOM_SK_Vec gshgeom_sk_sub(GSHGEOM_SK_Vec a, GSHGEOM_SK_Vec b) {
  return gshgeom_sk_v(a.x - b.x, a.y - b.y);
}

static inline GSHGEOM_SK_Vec gshgeom_sk_scale(double s, GSHGEOM_SK_Vec a) {
  return gshgeom_sk_v(s * a.x, s * a.y);
}

static inline double gshgeom_sk_dot(GSHGEOM_SK_Vec a, GSHGEOM_SK_Vec b) {
  return a.x * b.x + a.y * b.y;
}

static inline double gshgeom_sk_cross(GSHGEOM_SK_Vec a, GSHGEOM_SK_Vec b) {
  return a.x * b.y - a.y * b.x;
}

/** orient.h's left normal (the tangent turned by a negative quarter turn). */
static inline GSHGEOM_SK_Vec gshgeom_sk_left(GSHGEOM_SK_Vec t) {
  return gshgeom_sk_v(t.y, -t.x);
}

static inline GSHGEOM_SK_Vec gshgeom_sk_right(GSHGEOM_SK_Vec t) {
  return gshgeom_sk_v(-t.y, t.x);
}

/** The length, by gmath_sqrt (emit.c). */
double gshgeom_sk_norm(GSHGEOM_SK_Vec a);

/** The unit vector of a nonzero vector. */
GSHGEOM_SK_Vec gshgeom_sk_unit(GSHGEOM_SK_Vec a);

#endif /* GHOTI_IO_GSHGEOM_STROKE_INTERNAL_H */
