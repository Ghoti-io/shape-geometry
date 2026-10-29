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
 * The vertices of a path and the directions at them (vertices.h).  The directions are SG5's, read by the code
 * the stroke expansion reads them with (gshgeom_sk_start_dir, gshgeom_sk_end_dir, gshgeom_sk_unit): this file
 * only decides which segment a vertex looks to on each side, and adds the bisector.
 */

#include <ghoti.io/shape-geometry/vertices.h>

#include "../stroke/stroke_internal.h"

#include <string.h>

/** One segment as the walk reads it: whether it has extent, and its start and end directions (raw vectors). */
typedef struct Seg {
  int extent;
  GSHGEOM_SK_Vec start;
  GSHGEOM_SK_Vec end;
} Seg;

/** The segment of element @p el: a line, a curve or a close has one; a move has none. */
static void read_seg(const GSHGEOM_Element * el, Seg * s) {
  memset(s, 0, sizeof *s);
  if (el->verb == GSHGEOM_VERB_MOVE) {
    return;
  }
  GSHGEOM_SK_Seg k;
  memset(&k, 0, sizeof k);
  k.n = gshgeom_curve_point_count(el->curve.kind) - 1;
  for (int i = 0; i <= k.n; i++) {
    k.p[i] = gshgeom_sk_v(el->curve.p[i].x, el->curve.p[i].y);
  }
  s->extent = gshgeom_sk_start_dir(&k, &s->start) >= 0 && gshgeom_sk_end_dir(&k, &s->end) >= 0;
}

/** Element @p i, or ::GSHGEOM_ERR_INTERNAL. */
static GSHGEOM_Result element(const GSHGEOM_Path * path, size_t i, GSHGEOM_Element * el) {
  return gshgeom_path_element(path, i, el) == GSHGEOM_OK ? GSHGEOM_OK : GSHGEOM_ERR_INTERNAL;
}

static GSHGEOM_Point to_point(GSHGEOM_SK_Vec v) {
  GSHGEOM_Point p = {v.x, v.y, 0.0};
  return p;
}

/** The marker direction: the bisector of two unit directions, the first turned a quarter turn at a reversal. */
static GSHGEOM_Point bisector(GSHGEOM_Point a, GSHGEOM_Point b) {
  GSHGEOM_SK_Vec sum = gshgeom_sk_add(gshgeom_sk_v(a.x, a.y), gshgeom_sk_v(b.x, b.y));
  if (sum.x == 0.0 && sum.y == 0.0) {
    return gshgeom_normal_right(a);
  }
  return to_point(gshgeom_sk_unit(sum));
}

GSHGEOM_Result gshgeom_path_vertex_count(const GSHGEOM_Path * path, size_t * count) {
  if (!path || !count || gshgeom_path_dimension(path) != 2) {
    return GSHGEOM_ERR_INVALID;
  }
  *count = gshgeom_path_element_count(path);
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_path_vertices(const GSHGEOM_Path * path, GSHGEOM_Vertex_Fn fn, void * ctx) {
  if (!path || !fn || gshgeom_path_dimension(path) != 2) {
    return GSHGEOM_ERR_INVALID;
  }
  const size_t n = gshgeom_path_element_count(path);
  size_t sa = 0; /* the first element of the subpath */
  size_t sub = 0;
  while (sa < n) {
    /* The subpath: from sa to the element before the next move, or through the close that ends it. */
    GSHGEOM_Element el;
    GSHGEOM_Result r = element(path, sa, &el);
    if (r != GSHGEOM_OK) {
      return r;
    }
    size_t sb = sa + 1;
    int closed = el.verb == GSHGEOM_VERB_CLOSE;
    while (!closed && sb < n) {
      r = element(path, sb, &el);
      if (r != GSHGEOM_OK) {
        return r;
      }
      if (el.verb == GSHGEOM_VERB_MOVE) {
        break;
      }
      sb++;
      closed = el.verb == GSHGEOM_VERB_CLOSE;
    }
    /* A closed subpath is read round: what lies before its first segment with extent is its last, and after its last is its first. */
    Seg first, last;
    int have_first = 0, have_last = 0;
    memset(&first, 0, sizeof first);
    memset(&last, 0, sizeof last);
    if (closed) {
      for (size_t i = sa; i < sb && !have_first; i++) {
        Seg s;
        if ((r = element(path, i, &el)) != GSHGEOM_OK) {
          return r;
        }
        read_seg(&el, &s);
        if (s.extent) {
          first = s;
          have_first = 1;
        }
      }
      for (size_t i = sb; i > sa && !have_last; i--) {
        Seg s;
        if ((r = element(path, i - 1, &el)) != GSHGEOM_OK) {
          return r;
        }
        read_seg(&el, &s);
        if (s.extent) {
          last = s;
          have_last = 1;
        }
      }
    }
    Seg cur_in;           /* the nearest segment with extent at or before the vertex */
    int have_in = 0;
    Seg next;             /* the nearest one after it, at element `nj` (`sb` when there is none) */
    size_t nj = sa;
    memset(&cur_in, 0, sizeof cur_in);
    memset(&next, 0, sizeof next);
    for (size_t i = sa; i < sb; i++) {
      if ((r = element(path, i, &el)) != GSHGEOM_OK) {
        return r;
      }
      const GSHGEOM_Verb verb = el.verb;
      GSHGEOM_Vertex v;
      memset(&v, 0, sizeof v);
      v.element = i;
      v.subpath = sub;
      v.point = el.curve.p[gshgeom_curve_point_count(el.curve.kind) - 1];
      v.point.z = 0.0;
      Seg here;
      read_seg(&el, &here);
      if (here.extent) {
        cur_in = here;
        have_in = 1;
      }
      /* The scan for the next segment with extent only moves forward, so a run of segments without extent is
       * passed over once. */
      if (nj <= i) {
        nj = sb;
        for (size_t j = i + 1; j < sb; j++) {
          Seg s;
          if ((r = element(path, j, &el)) != GSHGEOM_OK) {
            return r;
          }
          read_seg(&el, &s);
          if (s.extent) {
            next = s;
            nj = j;
            break;
          }
        }
      }
      const int have_next = nj < sb;
      const Seg * in_seg = have_in ? &cur_in : (have_last ? &last : NULL);
      const Seg * out_seg = have_next ? &next : (have_first ? &first : NULL);
      GSHGEOM_Point in = {0.0, 0.0, 0.0}, out = {0.0, 0.0, 0.0};
      if (in_seg) {
        in = to_point(gshgeom_sk_unit(in_seg->end));
        v.flags |= GSHGEOM_VERTEX_HAS_IN;
      }
      if (out_seg) {
        out = to_point(gshgeom_sk_unit(out_seg->start));
        v.flags |= GSHGEOM_VERTEX_HAS_OUT;
      }
      v.in = in;
      v.out = out;
      if (in_seg && out_seg) {
        v.direction = bisector(in, out);
      } else if (in_seg) {
        v.direction = in;
      } else if (out_seg) {
        v.direction = out;
      } else {
        v.direction.x = 1.0;
      }
      if (i == 0) {
        v.flags |= GSHGEOM_VERTEX_PATH_START;
      }
      if (i + 1 == n) {
        v.flags |= GSHGEOM_VERTEX_PATH_END;
      }
      if (verb == GSHGEOM_VERB_MOVE) {
        v.flags |= GSHGEOM_VERTEX_SUBPATH_START;
      }
      if (i + 1 == sb) {
        v.flags |= GSHGEOM_VERTEX_SUBPATH_END;
      }
      if (closed) {
        v.flags |= GSHGEOM_VERTEX_CLOSED;
      }
      r = fn(&v, ctx);
      if (r != GSHGEOM_OK) {
        return r;
      }
    }
    sa = sb;
    sub++;
  }
  return GSHGEOM_OK;
}
