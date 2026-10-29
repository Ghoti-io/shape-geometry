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
 * The expansion (stroke.h, SG4 to SG18): the checks, the not-invertible and
 * zero-width statuses, the centerline read into subpaths of segments with
 * extent (SG5) in build space, the magnitude rule of the arithmetic bound
 * (emit.c (a)), the parts of SG12, and the pieces in SG18's order.
 */

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/math/basic.h>
#include <ghoti.io/shape-geometry/arc.h>
#include <ghoti.io/shape-geometry/device.h>

#include "../core/grow.h"
#include "../flatten/flatten_internal.h"
#include "stroke_internal.h"

#include <string.h>

typedef GSHGEOM_SK_Sub Sub;
typedef GSHGEOM_SK_Reader Reader;

static GSHGEOM_Result push_sub(GSHGEOM_SK_Ctx * c, Reader * rd, GSHGEOM_SK_Vec start, GSHGEOM_Point model,
    size_t first_el) {
  void * g = rd->subs;
  GSHGEOM_Result r = gshgeom_grow(&c->alloc, &g, &rd->sub_cap, rd->sub_n + 1, sizeof *rd->subs);
  if (r != GSHGEOM_OK) {
    return r;
  }
  rd->subs = g;
  Sub * s = &rd->subs[rd->sub_n++];
  memset(s, 0, sizeof *s);
  s->first_seg = rd->seg_n;
  s->start = start;
  s->start_model = model;
  s->first_el = first_el;
  s->end_el = first_el;
  return GSHGEOM_OK;
}

static int has_extent(const GSHGEOM_SK_Seg * s) {
  for (int k = 1; k <= s->n; k++) {
    if (s->p[k].x != s->p[0].x || s->p[k].y != s->p[0].y) {
      return 1;
    }
  }
  return 0;
}

/** Add the segment @p s (with extent) to subpath @p sb. */
static GSHGEOM_Result push_seg(GSHGEOM_SK_Ctx * c, Reader * rd, Sub * sb, const GSHGEOM_SK_Seg * s) {
  void * g = rd->segs;
  GSHGEOM_Result r = gshgeom_grow(&c->alloc, &g, &rd->seg_cap, rd->seg_n + 1, sizeof *rd->segs);
  if (r != GSHGEOM_OK) {
    return r;
  }
  rd->segs = g;
  rd->segs[rd->seg_n++] = *s;
  sb->seg_count++;
  return GSHGEOM_OK;
}

/** The segment of element @p i, piece @p piece, from the curve @p el (model space) and its image @p cv (build space). */
static GSHGEOM_SK_Seg make_seg(size_t i, size_t piece, const GSHGEOM_Curve * el, const GSHGEOM_Curve * cv) {
  GSHGEOM_SK_Seg s;
  memset(&s, 0, sizeof s);
  s.elem = i;
  s.piece = piece;
  s.n = gshgeom_curve_point_count(cv->kind) - 1;
  s.rational = gshgeom_curve_is_rational(cv->kind);
  for (int k = 0; k <= s.n; k++) {
    s.model_p[k] = gshgeom_sk_v(el->p[k].x, el->p[k].y);
    s.p[k] = gshgeom_sk_v(cv->p[k].x, cv->p[k].y);
    s.w[k] = s.rational ? cv->w[k] : 1.0;
  }
  return s;
}

/** Read the centerline: subpaths and segments with extent, in build space (SG5, SG10). */
static GSHGEOM_Result read_centerline(GSHGEOM_SK_Ctx * c, const GSHGEOM_Path * path, Reader * rd) {
  size_t n = gshgeom_path_element_count(path);
  int open = 0;
  GSHGEOM_SK_Vec start = {0.0, 0.0};
  GSHGEOM_Point start_model = {0.0, 0.0, 0.0};
  for (size_t i = 0; i < n; i++) {
    GSHGEOM_Verb verb;
    GSHGEOM_Result r = gshgeom_path_verb(path, i, &verb);
    if (r != GSHGEOM_OK) {
      return GSHGEOM_ERR_INTERNAL;
    }
    /* an element is one curve, or an arc's pieces: the arc is one segment of its subpath (SG5) whose body is the
     * union of its pieces' bodies (AR13) */
    GSHGEOM_Curve model[GSHGEOM_ARC_MAX_PIECES];
    int nc = 1;
    if (verb == GSHGEOM_VERB_ARC) {
      GSHGEOM_Arc arc;
      r = gshgeom_path_arc(path, i, &arc);
      if (r == GSHGEOM_OK) {
        r = gshgeom_arc_pieces(&arc, model, &nc);
      }
      if (r != GSHGEOM_OK) {
        return GSHGEOM_ERR_INTERNAL;
      }
    } else {
      GSHGEOM_Element el;
      r = gshgeom_path_element(path, i, &el);
      if (r != GSHGEOM_OK) {
        return GSHGEOM_ERR_INTERNAL;
      }
      model[0] = el.curve;
    }
    if (verb == GSHGEOM_VERB_MOVE) {
      GSHGEOM_Curve cv = model[0];
      if (!c->local) {
        r = gshgeom_map_apply_curve(&c->map, &model[0], &cv);
        if (r != GSHGEOM_OK) {
          return r;
        }
      }
      start = gshgeom_sk_v(cv.p[0].x, cv.p[0].y);
      start_model = model[0].p[0];
      r = push_sub(c, rd, start, start_model, i + 1);
      if (r != GSHGEOM_OK) {
        return r;
      }
      open = 1;
      continue;
    }
    if (!open) {
      /* a drawing verb after a close starts a subpath at the closed one's start */
      r = push_sub(c, rd, start, start_model, i);
      if (r != GSHGEOM_OK) {
        return r;
      }
      open = 1;
    }
    Sub * sb = &rd->subs[rd->sub_n - 1];
    sb->draws = 1;
    sb->end_el = i + 1;
    int kept = 0; /* a piece of this arc has been kept: the next one continues it, and no join is made between them */
    for (int k = 0; k < nc; k++) {
      GSHGEOM_Curve cv = model[k];
      if (!c->local) {
        r = gshgeom_map_apply_curve(&c->map, &model[k], &cv);
        if (r != GSHGEOM_OK) {
          return r;
        }
      }
      GSHGEOM_SK_Seg s = make_seg(i, (size_t)k, &model[k], &cv);
      if (has_extent(&s)) {
        s.seam = kept;
        kept = 1;
        r = push_seg(c, rd, sb, &s);
        if (r != GSHGEOM_OK) {
          return r;
        }
      }
    }
    if (verb == GSHGEOM_VERB_CLOSE) {
      sb->closed = 1;
      open = 0;
    }
  }
  return GSHGEOM_OK;
}

/** The class of the vertex at the start of segment @p next, against @p prev, on the points before the map (SG8, SG10). */
static void classify(const GSHGEOM_SK_Seg * prev, GSHGEOM_SK_Seg * next) {
  GSHGEOM_SK_Seg pm = *prev, nm = *next;
  for (int k = 0; k <= pm.n; k++) {
    pm.p[k] = pm.model_p[k];
  }
  for (int k = 0; k <= nm.n; k++) {
    nm.p[k] = nm.model_p[k];
  }
  GSHGEOM_SK_Vec a, b;
  int qi = gshgeom_sk_end_dir(&pm, &a), si = gshgeom_sk_start_dir(&nm, &b);
  if (qi < 0 || si < 0) {
    return;
  }
  next->vtx_class = gshgeom_sk_vertex_class(pm.p[qi], pm.p[pm.n], nm.p[si]);
  next->vtx_valid = 1;
}

static void classify_all(Reader * rd) {
  for (size_t i = 0; i < rd->sub_n; i++) {
    const Sub * s = &rd->subs[i];
    for (size_t k = 1; k < s->seg_count; k++) {
      classify(&rd->segs[s->first_seg + k - 1], &rd->segs[s->first_seg + k]);
    }
    if (s->closed && s->seg_count > 0) {
      classify(&rd->segs[s->first_seg + s->seg_count - 1], &rd->segs[s->first_seg]);
    }
  }
}

static double fabs_(double x) {
  return x < 0 ? -x : x;
}

/** The magnitude rule (emit.c (a)): every quantity the arithmetic bound is proved for is within 2^30 pixels. */
static int point_ok(const GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec p) {
  if (c->local) {
    const GSHGEOM_Map * m = &c->map;
    double r0 = fabs_(m->m[0][0] * p.x) + fabs_(m->m[0][1] * p.y) + fabs_(m->t[0]);
    double r1 = fabs_(m->m[1][0] * p.x) + fabs_(m->m[1][1] * p.y) + fabs_(m->t[1]);
    return r0 <= c->mag_limit && r1 <= c->mag_limit;
  }
  return fabs_(p.x) <= c->mag_limit && fabs_(p.y) <= c->mag_limit;
}

static int magnitudes_ok(const GSHGEOM_SK_Ctx * c, const Reader * rd, double h_dev_widest) {
  double L = c->st.miter_limit > 2.0 ? c->st.miter_limit : 2.0;
  if (!(h_dev_widest * L * 2.0 <= c->mag_limit)) {
    return 0;
  }
  for (size_t i = 0; i < rd->sub_n; i++) {
    if (!point_ok(c, rd->subs[i].start)) {
      return 0;
    }
  }
  for (size_t i = 0; i < rd->seg_n; i++) {
    const GSHGEOM_SK_Seg * s = &rd->segs[i];
    for (int k = 0; k <= s->n; k++) {
      if (!point_ok(c, s->p[k])) {
        return 0;
      }
      /* the build-space arithmetic (cubes of a difference, squares of a cross product) stays inside the range of a
       * double: a nonzero difference of neighbours is at least 2^-200 (the other end of the range is bounded by the
       * device magnitude and the map's smallest scale) */
      if (k > 0) {
        double dx = fabs_(s->p[k].x - s->p[k - 1].x), dy = fabs_(s->p[k].y - s->p[k - 1].y);
        double dm = dx > dy ? dx : dy;
        if (dm != 0.0 && dm < GSHGEOM_SK_BUILD_MIN) {
          return 0;
        }
      }
      GSHGEOM_SK_Vec d = gshgeom_sk_sub(s->p[k], s->p[0]);
      if (!(c->norm * gshgeom_sk_norm(d) <= c->mag_limit)) {
        return 0;
      }
    }
  }
  return 1;
}

/** Append the pieces of arc element @p e of @p path to @p mp as rational quadratics, mapped by @p map when it is not NULL. */
static GSHGEOM_Result clip_arc(GSHGEOM_Path * mp, const GSHGEOM_Path * path, size_t e, const GSHGEOM_Map * map) {
  GSHGEOM_Arc arc;
  GSHGEOM_Curve pieces[GSHGEOM_ARC_MAX_PIECES];
  int np = 0;
  GSHGEOM_Result r = gshgeom_path_arc(path, e, &arc);
  if (r == GSHGEOM_OK) {
    r = gshgeom_arc_pieces(&arc, pieces, &np);
  }
  for (int k = 0; k < np && r == GSHGEOM_OK; k++) {
    GSHGEOM_Curve cv = pieces[k];
    if (map) {
      r = gshgeom_map_apply_curve(map, &pieces[k], &cv);
    }
    if (r == GSHGEOM_OK) {
      r = gshgeom_path_rational_quadratic_to(mp, cv.p[1], cv.p[2], cv.w);
    }
  }
  return r;
}

/** The clip's region F (SG12): the closed subpaths with extent, flattened after the map at the tolerance. */
static GSHGEOM_Result build_clip(GSHGEOM_SK_Ctx * c, const GSHGEOM_Path * path, const Reader * rd,
    GSHGEOM_Flat_Path ** out) {
  GSHGEOM_Path * mp = NULL;
  GSHGEOM_Result r = gshgeom_path_create(&c->alloc, &c->limits, 2, &mp);
  for (size_t i = 0; i < rd->sub_n && r == GSHGEOM_OK; i++) {
    const Sub * s = &rd->subs[i];
    if (!s->closed || s->seg_count == 0) {
      continue;
    }
    r = gshgeom_path_move_to(mp, s->start_model);
    for (size_t e = s->first_el; e < s->end_el && r == GSHGEOM_OK; e++) {
      GSHGEOM_Element el;
      r = gshgeom_path_element(path, e, &el);
      if (r == GSHGEOM_ERR_UNSUPPORTED) {
        r = clip_arc(mp, path, e, NULL);
        continue;
      }
      if (r != GSHGEOM_OK) {
        break;
      }
      const GSHGEOM_Point * p = el.curve.p;
      switch (el.verb) {
      case GSHGEOM_VERB_LINE:
        r = gshgeom_path_line_to(mp, p[1]);
        break;
      case GSHGEOM_VERB_QUADRATIC:
        r = gshgeom_path_quadratic_to(mp, p[1], p[2]);
        break;
      case GSHGEOM_VERB_CUBIC:
        r = gshgeom_path_cubic_to(mp, p[1], p[2], p[3]);
        break;
      case GSHGEOM_VERB_RATIONAL_QUADRATIC:
        r = gshgeom_path_rational_quadratic_to(mp, p[1], p[2], el.curve.w);
        break;
      case GSHGEOM_VERB_RATIONAL_CUBIC:
        r = gshgeom_path_rational_cubic_to(mp, p[1], p[2], p[3], el.curve.w);
        break;
      case GSHGEOM_VERB_CLOSE:
        r = gshgeom_path_close(mp);
        break;
      default:
        r = GSHGEOM_ERR_INTERNAL;
        break;
      }
    }
  }
  GSHGEOM_Device_Path * dp = NULL;
  if (r == GSHGEOM_OK) {
    r = gshgeom_device_path_from_model(mp, &c->map, &c->alloc, &c->limits, &dp);
  }
  if (r == GSHGEOM_OK) {
    r = gshgeom_device_path_flatten_with_tolerance(dp, c->tau, &c->alloc, &c->limits, out);
  }
  gshgeom_device_path_destroy(dp);
  gshgeom_path_destroy(mp);
  return r;
}

/** The curve outline's clip (SG32): the closed subpaths with extent, `M` applied to every element. */
static GSHGEOM_Result build_clip_curves(GSHGEOM_SK_Ctx * c, const GSHGEOM_Path * path, const Reader * rd,
    GSHGEOM_Path ** out) {
  GSHGEOM_Path * mp = NULL;
  GSHGEOM_Result r = gshgeom_path_create(&c->alloc, &c->limits, 2, &mp);
  for (size_t i = 0; i < rd->sub_n && r == GSHGEOM_OK; i++) {
    const Sub * s = &rd->subs[i];
    if (!s->closed || s->seg_count == 0) {
      continue;
    }
    GSHGEOM_Point sp;
    r = gshgeom_map_apply_point(&c->map, s->start_model, &sp);
    if (r == GSHGEOM_OK) {
      r = gshgeom_path_move_to(mp, sp);
    }
    for (size_t e = s->first_el; e < s->end_el && r == GSHGEOM_OK; e++) {
      GSHGEOM_Element el;
      r = gshgeom_path_element(path, e, &el);
      if (r == GSHGEOM_ERR_UNSUPPORTED) {
        r = clip_arc(mp, path, e, &c->map);
        continue;
      }
      if (r != GSHGEOM_OK) {
        break;
      }
      GSHGEOM_Curve cv;
      r = gshgeom_map_apply_curve(&c->map, &el.curve, &cv);
      if (r != GSHGEOM_OK) {
        break;
      }
      const GSHGEOM_Point * p = cv.p;
      switch (el.verb) {
      case GSHGEOM_VERB_LINE:
        r = gshgeom_path_line_to(mp, p[1]);
        break;
      case GSHGEOM_VERB_QUADRATIC:
        r = gshgeom_path_quadratic_to(mp, p[1], p[2]);
        break;
      case GSHGEOM_VERB_CUBIC:
        r = gshgeom_path_cubic_to(mp, p[1], p[2], p[3]);
        break;
      case GSHGEOM_VERB_RATIONAL_QUADRATIC:
        r = gshgeom_path_rational_quadratic_to(mp, p[1], p[2], cv.w);
        break;
      case GSHGEOM_VERB_RATIONAL_CUBIC:
        r = gshgeom_path_rational_cubic_to(mp, p[1], p[2], p[3], cv.w);
        break;
      case GSHGEOM_VERB_CLOSE:
        r = gshgeom_path_close(mp);
        break;
      default:
        r = GSHGEOM_ERR_INTERNAL;
        break;
      }
    }
  }
  if (r != GSHGEOM_OK) {
    gshgeom_path_destroy(mp);
    return r;
  }
  *out = mp;
  return GSHGEOM_OK;
}

/** SG11's effective tangent of the degenerate subpath i (raw, build space). */
GSHGEOM_SK_Vec gshgeom_sk_effective_tangent(const GSHGEOM_SK_Ctx * c, const GSHGEOM_SK_Reader * rd, size_t i) {
  for (size_t j = i; j-- > 0;) {
    if (rd->subs[j].seg_count) {
      GSHGEOM_SK_Vec d;
      (void)gshgeom_sk_end_dir(&rd->segs[rd->subs[j].first_seg + rd->subs[j].seg_count - 1], &d);
      return d;
    }
  }
  for (size_t j = i + 1; j < rd->sub_n; j++) {
    if (rd->subs[j].seg_count) {
      GSHGEOM_SK_Vec d;
      (void)gshgeom_sk_start_dir(&rd->segs[rd->subs[j].first_seg], &d);
      return d;
    }
  }
  /* no segment with extent anywhere: local +x (under device width, the map's image of it) */
  return c->local ? gshgeom_sk_v(1.0, 0.0) : gshgeom_sk_v(c->map.m[0][0], c->map.m[1][0]);
}

/** One subpath's pieces in SG18's order. */

static GSHGEOM_Result stroke_sub(GSHGEOM_SK_Ctx * c, const Reader * rd, size_t i) {
  if (c->dash.active) {
    return gshgeom_sk_dash_sub(c, rd, i);
  }
  return gshgeom_sk_solid_sub(c, rd, i);
}

GSHGEOM_Result gshgeom_sk_solid_sub(GSHGEOM_SK_Ctx * c, const GSHGEOM_SK_Reader * rd, size_t i) {
  const Sub * s = &rd->subs[i];
  if (s->seg_count == 0) {
    if (!s->draws || c->st.degenerate == GSHGEOM_STROKE_DEGENERATE_OMIT) {
      return GSHGEOM_OK;
    }
    return gshgeom_sk_degenerate(c, s->start, gshgeom_sk_effective_tangent(c, rd, i));
  }
  const GSHGEOM_SK_Seg * seg = &rd->segs[s->first_seg];
  size_t n = s->seg_count;
  GSHGEOM_Result r = GSHGEOM_OK;
  if (!s->closed) {
    GSHGEOM_SK_Vec d;
    (void)gshgeom_sk_start_dir(&seg[0], &d);
    r = gshgeom_sk_cap(c, seg[0].p[0], gshgeom_sk_scale(-1.0, d));
  }
  for (size_t k = 0; k < n && r == GSHGEOM_OK; k++) {
    r = gshgeom_sk_body(c, &seg[k]);
    if (r == GSHGEOM_OK && k + 1 < n) {
      if (!seg[k + 1].seam) { /* a boundary between two pieces of one arc is not a vertex (AR13) */
        r = gshgeom_sk_join(c, &seg[k], &seg[k + 1]);
      }
    } else if (r == GSHGEOM_OK && s->closed) {
      r = gshgeom_sk_join(c, &seg[k], &seg[0]);
    }
  }
  if (r == GSHGEOM_OK && !s->closed) {
    GSHGEOM_SK_Vec d;
    (void)gshgeom_sk_end_dir(&seg[n - 1], &d);
    r = gshgeom_sk_cap(c, seg[n - 1].p[seg[n - 1].n], d);
  }
  return r;
}

GSHGEOM_Result gshgeom_stroke_expand(const GSHGEOM_Path * centerline, const GSHGEOM_Stroke * stroke,
    const GSHGEOM_Map * map, const GSHGEOM_Allocator * alloc, const GSHGEOM_Limits * limits,
    GSHGEOM_Stroke_Outline ** out) {
  return gshgeom_stroke_expand_with_tolerance(
      centerline, stroke, map, GSHGEOM_FLATTEN_TOLERANCE, alloc, limits, out);
}

GSHGEOM_Result gshgeom_stroke_expand_with_tolerance(const GSHGEOM_Path * centerline, const GSHGEOM_Stroke * stroke,
    const GSHGEOM_Map * map, int32_t tolerance, const GSHGEOM_Allocator * alloc, const GSHGEOM_Limits * limits,
    GSHGEOM_Stroke_Outline ** out) {
  return gshgeom_stroke_expand_capped(centerline, stroke, map, tolerance, alloc, limits, GSHGEOM_STROKE_MAX_DEPTH, out);
}

/** What a curve outline asks of the shared expansion (stroke_curves.h). */
typedef struct CurveReq {
  double t;
  GSHGEOM_Stroke_Curves_Form form;
} CurveReq;

/** The expansion of stroke.h (@p cr NULL: 26.6 polygons in @p po) or of stroke_curves.h (curves in @p pc). */
static GSHGEOM_Result expand_core(const GSHGEOM_Path * centerline, const GSHGEOM_Stroke * stroke,
    const GSHGEOM_Map * map, int32_t tolerance, const CurveReq * cr, const GSHGEOM_Allocator * alloc,
    const GSHGEOM_Limits * limits, int depth_cap, int conform, int * redo, GSHGEOM_Stroke_Outline ** po,
    GSHGEOM_Stroke_Curves ** pc) {
  GSHGEOM_Result r = gshgeom_map_check(map);
  if (r != GSHGEOM_OK) {
    return r;
  }
  if (map->dim != 2) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_SK_Ctx c;
  memset(&c, 0, sizeof c);
  r = gshgeom_stroke_canonical(stroke, &c.st);
  if (r != GSHGEOM_OK) {
    return r;
  }
  if (gshgeom_path_dimension(centerline) != 2) {
    return GSHGEOM_ERR_UNSUPPORTED; /* the definition is planar */
  }
  if (!alloc) {
    alloc = gshgeom_allocator_default();
  }
  c.alloc = *alloc;
  gshgeom_limits_default(&c.limits);
  if (limits) {
    if (limits->max_elements) {
      c.limits.max_elements = limits->max_elements;
    }
    if (limits->max_points) {
      c.limits.max_points = limits->max_points;
    }
    c.limits.max_dash_steps = limits->max_dash_steps;
  }
  GSHGEOM_SK_Curves cvs;
  memset(&cvs, 0, sizeof cvs);
  c.conform = conform && !cr;
  c.tau = tolerance;
  c.budget = ((double)tolerance - GSHGEOM_SK_MARGIN) / 64.0;
  c.rho = (GSHGEOM_SK_ROUND + GSHGEOM_SK_ETA) / 64.0;
  c.mag_limit = GSHGEOM_SK_MAG_LIMIT;
  if (cr) {
    /* SG35: t / 2 for the approximation of pieces; the arithmetic's reach, 2^-43 of the largest magnitude, must be
     * within t / 8, which bounds the magnitudes the construction is certified for by t 2^40 (and by 2^30) */
    cvs.t = cr->t;
    cvs.budget = cr->t / 2.0;
    c.cv = &cvs;
    double lim = cr->t * 0x1p40;
    c.mag_limit = lim < GSHGEOM_SK_MAG_LIMIT ? lim : GSHGEOM_SK_MAG_LIMIT;
  }
  c.map = *map;
  c.local = c.st.width_space == GSHGEOM_STROKE_WIDTH_LOCAL;
  c.work_cap = 4 * c.limits.max_points + 4096;
  c.dash_cap = c.limits.max_dash_steps && c.limits.max_dash_steps < (uint64_t)SIZE_MAX ? (size_t)c.limits.max_dash_steps : c.work_cap;
  c.depth_cap = depth_cap < 1 ? 1 : (depth_cap > GSHGEOM_STROKE_MAX_DEPTH ? GSHGEOM_STROKE_MAX_DEPTH : depth_cap);
  GSHGEOM_Stroke_Outline * o = NULL;
  GSHGEOM_Stroke_Curves * cu = NULL;
  if (cr) {
    cu = gcu_allocator_calloc(&c.alloc, 1, sizeof *cu);
    if (!cu) {
      return GSHGEOM_ERR_OOM;
    }
    cu->alloc = c.alloc;
    cu->form = cr->form;
    cu->alignment = c.st.alignment;
    cu->rule = c.st.alignment_rule;
    cu->complement = c.st.alignment == GSHGEOM_STROKE_ALIGN_OUTSIDE;
  } else {
    o = gcu_allocator_calloc(&c.alloc, 1, sizeof *o);
    if (!o) {
      return GSHGEOM_ERR_OOM;
    }
    o->alloc = c.alloc;
    o->alignment = c.st.alignment;
    o->rule = c.st.alignment_rule;
    o->complement = c.st.alignment == GSHGEOM_STROKE_ALIGN_OUTSIDE;
  }
#define EXPAND_FAIL(res)                                                                                               \
  do {                                                                                                                 \
    gshgeom_stroke_outline_destroy(o);                                                                                 \
    gshgeom_stroke_curves_destroy(cu);                                                                                 \
    return (res);                                                                                                      \
  } while (0)
#define EXPAND_DONE(st)                                                                                                \
  do {                                                                                                                 \
    if (cr) {                                                                                                          \
      cu->status = (st);                                                                                               \
      *pc = cu;                                                                                                        \
    } else {                                                                                                           \
      o->status = (st);                                                                                                \
      *po = o;                                                                                                         \
    }                                                                                                                  \
    return GSHGEOM_OK;                                                                                                 \
  } while (0)
  int invertible = 0;
  r = gshgeom_stroke_map_invertible(map, &invertible);
  if (r != GSHGEOM_OK) {
    EXPAND_FAIL(r);
  }
  if (!invertible) {
    EXPAND_DONE(GSHGEOM_STROKE_NOT_INVERTIBLE); /* SG15, decided before SG14 */
  }
  if (c.st.width == 0.0) {
    EXPAND_DONE(GSHGEOM_STROKE_EMPTY); /* SG14 */
  }
  if (c.local) {
    c.lin[0][0] = map->m[0][0];
    c.lin[0][1] = map->m[0][1];
    c.lin[1][0] = map->m[1][0];
    c.lin[1][1] = map->m[1][1];
    double f = map->m[0][0] * map->m[0][0] + map->m[0][1] * map->m[0][1] + map->m[1][0] * map->m[1][0] +
        map->m[1][1] * map->m[1][1];
    c.norm = gmath_sqrt(f) * (1.0 + 0x1p-40);
    if (cr) {
      /* the determinant's sign, exactly: a mirroring map reverses the orientation of every contour */
      GSHGEOM_Point o0 = {0.0, 0.0, 0.0}, ex = {map->m[0][0], map->m[1][0], 0.0}, ey = {map->m[0][1], map->m[1][1], 0.0};
      int sg = 0;
      (void)gshgeom_orient2d(o0, ex, ey, &sg);
      cvs.flip = sg < 0;
    }
  } else {
    c.lin[0][0] = c.lin[1][1] = 1.0;
    c.norm = 1.0;
  }
  Reader rd;
  memset(&rd, 0, sizeof rd);
  r = read_centerline(&c, centerline, &rd);
  if (r == GSHGEOM_OK) {
    classify_all(&rd);
  }
  int aligned = c.st.alignment != GSHGEOM_STROKE_ALIGN_CENTER;
  int any_closed = 0;
  for (size_t i = 0; i < rd.sub_n; i++) {
    if (rd.subs[i].closed && rd.subs[i].seg_count) {
      any_closed = 1;
    }
  }
  double h = c.st.width / 2.0;
  if (r == GSHGEOM_OK && !magnitudes_ok(&c, &rd, (aligned ? 2.0 : 1.0) * h * c.norm)) {
    r = GSHGEOM_ERR_LIMIT;
  }
  if (r == GSHGEOM_OK && c.st.dash_count) {
    r = gshgeom_sk_dash_begin(&c, centerline, &rd, stroke);
  }
  if (r == GSHGEOM_OK) {
    if (cr) {
      r = gshgeom_path_create(&c.alloc, &c.limits, 2, &cu->body);
    } else {
      r = gshgeom_flat_path_new(&c.alloc, &c.limits, &o->body);
    }
  }
  if (r == GSHGEOM_OK && aligned && any_closed) {
    if (cr) {
      r = gshgeom_path_create(&c.alloc, &c.limits, 2, &cu->aligned);
      if (r == GSHGEOM_OK) {
        r = build_clip_curves(&c, centerline, &rd, &cu->clip);
      }
    } else {
      r = gshgeom_flat_path_new(&c.alloc, &c.limits, &o->aligned);
      if (r == GSHGEOM_OK) {
        r = build_clip(&c, centerline, &rd, &o->clip);
      }
    }
  }
  for (size_t i = 0; i < rd.sub_n && r == GSHGEOM_OK; i++) {
    int to_aligned = aligned && rd.subs[i].closed && rd.subs[i].seg_count;
    c.h = to_aligned ? 2.0 * h : h;
    c.h_dev = c.h * c.norm;
    if (cr) {
      cvs.part = to_aligned ? cu->aligned : cu->body;
    } else {
      c.part = to_aligned ? o->aligned : o->body;
    }
    r = stroke_sub(&c, &rd, i);
  }
  gshgeom_sk_dash_end(&c);
  gcu_allocator_free(&c.alloc, rd.subs);
  gcu_allocator_free(&c.alloc, rd.segs);
  gcu_allocator_free(&c.alloc, c.scratch);
  gcu_allocator_free(&c.alloc, c.arcbuf);
  gcu_allocator_free(&c.alloc, c.snap);
  gcu_allocator_free(&c.alloc, cvs.el);
  gcu_allocator_free(&c.alloc, cvs.run);
  if (r == GSHGEOM_OK && cr && cr->form == GSHGEOM_STROKE_CURVES_CUBIC) {
    /* SG40: the cubic form is the rational form lowered at a quarter of the tolerance, by the public function */
    GSHGEOM_Path ** parts[2] = {&cu->body, &cu->aligned};
    for (int k = 0; k < 2 && r == GSHGEOM_OK; k++) {
      if (*parts[k]) {
        GSHGEOM_Path * low = NULL;
        r = gshgeom_stroke_curves_lower_cubic(*parts[k], cr->t / 4.0, &c.alloc, &c.limits, &low);
        if (r == GSHGEOM_OK) {
          gshgeom_path_destroy(*parts[k]);
          *parts[k] = low;
        }
      }
    }
  }
  if (r != GSHGEOM_OK) {
    EXPAND_FAIL(r);
  }
  if (c.inflated) {
    /* an inflated piece among conforming ones leaves hairline gaps along the seams it shares with them: draw it all again
     * with every piece inflated */
    *redo = 1;
    EXPAND_FAIL(GSHGEOM_OK);
  }
  EXPAND_DONE(GSHGEOM_STROKE_DRAWN);
#undef EXPAND_FAIL
#undef EXPAND_DONE
}

GSHGEOM_Result gshgeom_stroke_expand_capped(const GSHGEOM_Path * centerline, const GSHGEOM_Stroke * stroke,
    const GSHGEOM_Map * map, int32_t tolerance, const GSHGEOM_Allocator * alloc, const GSHGEOM_Limits * limits,
    int depth_cap, GSHGEOM_Stroke_Outline ** out) {
  if (!centerline || !stroke || !map || !out || tolerance < GSHGEOM_FLATTEN_MIN_TOLERANCE) {
    return GSHGEOM_ERR_INVALID;
  }
  int redo = 0;
  GSHGEOM_Result r = expand_core(centerline, stroke, map, tolerance, NULL, alloc, limits, depth_cap, 1, &redo, out, NULL);
  if (r == GSHGEOM_OK && redo) {
    r = expand_core(centerline, stroke, map, tolerance, NULL, alloc, limits, depth_cap, 0, &redo, out, NULL);
  }
  return r;
}

GSHGEOM_Result gshgeom_stroke_curves_expand_capped(const GSHGEOM_Path * centerline, const GSHGEOM_Stroke * stroke,
    const GSHGEOM_Map * map, double tolerance, GSHGEOM_Stroke_Curves_Form form, const GSHGEOM_Allocator * alloc,
    const GSHGEOM_Limits * limits, int depth_cap, GSHGEOM_Stroke_Curves ** out) {
  if (!centerline || !stroke || !map || !out) {
    return GSHGEOM_ERR_INVALID;
  }
  if (!(tolerance >= GSHGEOM_STROKE_CURVES_MIN_TOLERANCE) || !(tolerance <= GSHGEOM_STROKE_CURVES_MAX_TOLERANCE) ||
      (unsigned)form >= (unsigned)GSHGEOM_STROKE_CURVES_FORM_COUNT) {
    return GSHGEOM_ERR_INVALID;
  }
  CurveReq cr = {tolerance, form};
  int redo = 0;
  return expand_core(centerline, stroke, map, GSHGEOM_FLATTEN_TOLERANCE, &cr, alloc, limits, depth_cap, 0, &redo, NULL, out);
}

GSHGEOM_Result gshgeom_stroke_curves_expand(const GSHGEOM_Path * centerline, const GSHGEOM_Stroke * stroke,
    const GSHGEOM_Map * map, double tolerance, GSHGEOM_Stroke_Curves_Form form, const GSHGEOM_Allocator * alloc,
    const GSHGEOM_Limits * limits, GSHGEOM_Stroke_Curves ** out) {
  return gshgeom_stroke_curves_expand_capped(
      centerline, stroke, map, tolerance, form, alloc, limits, GSHGEOM_STROKE_MAX_DEPTH, out);
}
