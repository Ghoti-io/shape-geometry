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
 * The seams of a flattened path (flatten.h, ::gshgeom_flat_path_cancel_seams).
 *
 * **What is cancelled.**  Every subpath is read as closed, as a rasterizer reads it, so a path is a multiset of directed
 * edges `a -> b` (the closing line of each subpath included).  The winding number at a point not on an edge is the signed
 * number of edges that cross a ray from it, an edge going down counting +1 and one going up -1, so the edge `a -> b` and
 * the edge `b -> a` cross every ray together or not at all and contribute `+1` and `-1`: a pair of them changes no winding
 * number anywhere.  The stroke expansion builds its outline from positively wound pieces, one for each segment, each join
 * and each cap (SG18), and the pieces of neighbours share their sides in opposite directions, so most edges of an outline
 * have such a twin.  An edge of zero length crosses no ray and is dropped too.
 *
 * **What is left is closed again.**  A closed subpath leaves and enters each of its points as often as it has edges there,
 * so every point of the multiset has as many edges leaving as entering, and removing a pair `a -> b`, `b -> a` takes one
 * of each from `a` and from `b`: the multiset that remains has that balance too.  It therefore splits into closed
 * walks, which is what the result holds, one subpath each.  A walk starts at the first unused edge in the order of
 * (tail, head) and follows, at each point, the first unused edge that leaves it, so the result is a function of the edges
 * alone and not of the order of the subpaths.  The walk can only end at its start (the balance), and an input for
 * which it does not is ::GSHGEOM_ERR_INTERNAL.
 *
 * **Cost.**  Two sorts and a binary search for each step of a walk: `O(E log E)` for `E` edges, and 16 bytes for each edge
 * twice over besides the result.
 */

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/shape-geometry/flatten.h>

#include "flatten_internal.h"

#include <stdlib.h>

/** An edge between two points in the order of their coordinates (`lo` first), and whether the path goes that way. */
typedef struct Pair {
  int32_t x0, y0, x1, y1;
  int32_t forward;
} Pair;

/** A directed edge. */
typedef struct Arc {
  int32_t x0, y0, x1, y1;
} Arc;

static int cmp_pair(const void * p, const void * q) {
  const Pair * a = p;
  const Pair * b = q;
  if (a->x0 != b->x0) {
    return a->x0 < b->x0 ? -1 : 1;
  }
  if (a->y0 != b->y0) {
    return a->y0 < b->y0 ? -1 : 1;
  }
  if (a->x1 != b->x1) {
    return a->x1 < b->x1 ? -1 : 1;
  }
  if (a->y1 != b->y1) {
    return a->y1 < b->y1 ? -1 : 1;
  }
  return 0;
}

static int cmp_arc(const void * p, const void * q) {
  const Arc * a = p;
  const Arc * b = q;
  if (a->x0 != b->x0) {
    return a->x0 < b->x0 ? -1 : 1;
  }
  if (a->y0 != b->y0) {
    return a->y0 < b->y0 ? -1 : 1;
  }
  if (a->x1 != b->x1) {
    return a->x1 < b->x1 ? -1 : 1;
  }
  if (a->y1 != b->y1) {
    return a->y1 < b->y1 ? -1 : 1;
  }
  return 0;
}

/** Whether the point (x0, y0) comes before (x1, y1). */
static int before(int32_t x0, int32_t y0, int32_t x1, int32_t y1) {
  return x0 < x1 || (x0 == x1 && y0 < y1);
}

/** The index of the first arc whose tail is not before (x, y). */
static size_t lower_bound(const Arc * arcs, size_t n, int32_t x, int32_t y) {
  size_t lo = 0;
  size_t hi = n;
  while (lo < hi) {
    size_t mid = lo + (hi - lo) / 2;
    if (before(arcs[mid].x0, arcs[mid].y0, x, y)) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }
  return lo;
}

GSHGEOM_Result gshgeom_flat_path_cancel_seams(const GSHGEOM_Flat_Path * path, const GSHGEOM_Allocator * alloc,
    const GSHGEOM_Limits * limits, GSHGEOM_Flat_Path ** out) {
  if (!path || !out) {
    return GSHGEOM_ERR_INVALID;
  }
  if (!alloc) {
    alloc = gshgeom_allocator_default();
  }
  const GSHGEOM_Device_Point * points = gshgeom_flat_path_points(path);
  size_t subpaths = gshgeom_flat_path_subpath_count(path);
  size_t edges = gshgeom_flat_path_point_count(path); /* one edge for each point: the next, or the closing line */
  if (edges > SIZE_MAX / sizeof(Pair)) {
    return GSHGEOM_ERR_OOM;
  }
  Pair * pairs = gcu_allocator_malloc(alloc, (edges ? edges : 1) * sizeof *pairs);
  Arc * arcs = gcu_allocator_malloc(alloc, (edges ? edges : 1) * sizeof *arcs);
  unsigned char * used = gcu_allocator_calloc(alloc, edges ? edges : 1, 1);
  GSHGEOM_Flat_Path * result = NULL;
  GSHGEOM_Result status = GSHGEOM_OK;
  if (!pairs || !arcs || !used) {
    status = GSHGEOM_ERR_OOM;
    goto done;
  }
  size_t n = 0;
  for (size_t s = 0; s < subpaths; s++) {
    GSHGEOM_Flat_Subpath sub;
    if (gshgeom_flat_path_subpath(path, s, &sub) != GSHGEOM_OK) {
      status = GSHGEOM_ERR_INTERNAL;
      goto done;
    }
    for (size_t k = 0; k < sub.count; k++) {
      GSHGEOM_Device_Point a = points[sub.first + k];
      GSHGEOM_Device_Point b = points[sub.first + (k + 1 == sub.count ? 0 : k + 1)];
      if (a.x == b.x && a.y == b.y) {
        continue;
      }
      Pair * e = &pairs[n++];
      e->forward = before(a.x, a.y, b.x, b.y);
      e->x0 = e->forward ? a.x : b.x;
      e->y0 = e->forward ? a.y : b.y;
      e->x1 = e->forward ? b.x : a.x;
      e->y1 = e->forward ? b.y : a.y;
    }
  }
  qsort(pairs, n, sizeof *pairs, cmp_pair);
  /* Each run of one edge between two points: as many go one way as the other cancel, and the rest stay as they were. */
  size_t m = 0;
  for (size_t i = 0; i < n;) {
    size_t j = i;
    long net = 0;
    while (j < n && cmp_pair(&pairs[i], &pairs[j]) == 0) {
      net += pairs[j].forward ? 1 : -1;
      j++;
    }
    for (long r = net < 0 ? -net : net; r > 0; r--) {
      Arc * a = &arcs[m++];
      a->x0 = net > 0 ? pairs[i].x0 : pairs[i].x1;
      a->y0 = net > 0 ? pairs[i].y0 : pairs[i].y1;
      a->x1 = net > 0 ? pairs[i].x1 : pairs[i].x0;
      a->y1 = net > 0 ? pairs[i].y1 : pairs[i].y0;
    }
    i = j;
  }
  qsort(arcs, m, sizeof *arcs, cmp_arc);
  status = gshgeom_flat_path_new(alloc, limits, &result);
  if (status != GSHGEOM_OK) {
    goto done;
  }
  for (size_t i = 0; i < m; i++) {
    if (used[i]) {
      continue;
    }
    GSHGEOM_Device_Point start = {arcs[i].x0, arcs[i].y0};
    status = gshgeom_flat_path_begin(result, start);
    if (status != GSHGEOM_OK) {
      goto done;
    }
    used[i] = 1;
    int32_t hx = arcs[i].x1;
    int32_t hy = arcs[i].y1;
    while (hx != start.x || hy != start.y) {
      GSHGEOM_Device_Point at = {hx, hy};
      status = gshgeom_flat_path_add(result, at);
      if (status != GSHGEOM_OK) {
        goto done;
      }
      size_t j = lower_bound(arcs, m, hx, hy);
      while (j < m && arcs[j].x0 == hx && arcs[j].y0 == hy && used[j]) {
        j++;
      }
      if (j >= m || arcs[j].x0 != hx || arcs[j].y0 != hy) {
        status = GSHGEOM_ERR_INTERNAL; /* a point with more edges entering than leaving: not a union of closed subpaths */
        goto done;
      }
      used[j] = 1;
      hx = arcs[j].x1;
      hy = arcs[j].y1;
    }
    gshgeom_flat_path_close_last(result);
  }
done:
  gcu_allocator_free(alloc, pairs);
  gcu_allocator_free(alloc, arcs);
  gcu_allocator_free(alloc, used);
  if (status != GSHGEOM_OK) {
    gshgeom_flat_path_destroy(result);
    return status;
  }
  *out = result;
  return GSHGEOM_OK;
}
