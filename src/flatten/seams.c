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
 * **Cost.**  Two sorts and a binary search for each step of a walk: `O(E log E)` for `E` edges, and 49 bytes for each edge
 * (the edges and the sorts' room, 24 bytes each, and a flag) besides the result and the walk's links.
 */

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/shape-geometry/flatten.h>

#include "flatten_internal.h"

#include <stdint.h>

/*
 * A point is one 64-bit key: x in the high half and y in the low, each with its sign bit flipped, so that the unsigned order of
 * keys is the order of points by x and then y.  The sorts below compare keys inline; the C library's qsort, through a comparison
 * function and a copy of each element by its size, was a fifth of a dashed stroke's whole render.  Two records the sort cannot
 * tell apart are the same edge, so the order the sort leaves them in does not reach the result: any correct sort gives the
 * same bytes.
 */
static uint64_t point_key(int32_t x, int32_t y) {
  return ((uint64_t)((uint32_t)x ^ 0x80000000u) << 32) | (uint64_t)((uint32_t)y ^ 0x80000000u);
}

static GSHGEOM_Device_Point key_point(uint64_t k) {
  GSHGEOM_Device_Point p;
  p.x = (int32_t)((uint32_t)(k >> 32) ^ 0x80000000u);
  p.y = (int32_t)((uint32_t)k ^ 0x80000000u);
  return p;
}

/** A directed edge, tail and head as keys.  The pass that cancels holds each edge with its ends in key order (`t <= h`)
 * and the direction the path went in @c forward; the walk holds the edges that remain, tail first. */
typedef struct Arc {
  uint64_t t, h;
  uint64_t forward;
} Arc;

static int arc_before(const Arc * a, const Arc * b) {
  return a->t < b->t || (a->t == b->t && a->h < b->h);
}

/**
 * Sort @p a[0..n) by (t, h): runs of 16 by insertion, then merged in passes between @p a and @p tmp (room for n), so no pass
 * copies back.  Returns the buffer that holds the result, @p a or @p tmp.
 */
static Arc * sort_arcs(Arc * a, Arc * tmp, size_t n) {
  enum { RUN = 16 };
  for (size_t lo = 0; lo < n; lo += RUN) {
    size_t hi = n - lo < RUN ? n : lo + RUN;
    for (size_t i = lo + 1; i < hi; i++) {
      Arc v = a[i];
      size_t j = i;
      while (j > lo && arc_before(&v, &a[j - 1])) {
        a[j] = a[j - 1];
        j--;
      }
      a[j] = v;
    }
  }
  Arc * src = a;
  Arc * dst = tmp;
  for (size_t width = RUN; width < n; width *= 2) {
    for (size_t lo = 0; lo < n; lo += 2 * width) {
      size_t mid = n - lo < width ? n : lo + width;
      size_t hi = n - lo < 2 * width ? n : lo + 2 * width;
      size_t i = lo, j = mid, k = lo;
      while (i < mid && j < hi) {
        dst[k++] = arc_before(&src[j], &src[i]) ? src[j++] : src[i++];
      }
      while (i < mid) {
        dst[k++] = src[i++];
      }
      while (j < hi) {
        dst[k++] = src[j++];
      }
    }
    Arc * t = src;
    src = dst;
    dst = t;
  }
  return src;
}

/**
 * The least index at or past @p j whose arc is not used, by the links of @p skip (each used arc links to the next index; a link
 * is shortened when it is followed).  The scan past the used arcs at one point is what made a point that many closed subpaths pass
 * through cost the square of their number: 160 000 loops through one vertex took 7 seconds.
 */
static size_t first_unused(size_t * skip, size_t j) {
  while (skip[j] != j) {
    skip[j] = skip[skip[j]];
    j = skip[j];
  }
  return j;
}

/** The index of the first arc whose tail is not before @p key. */
static size_t lower_bound(const Arc * arcs, size_t n, uint64_t key) {
  size_t lo = 0;
  size_t hi = n;
  while (lo < hi) {
    size_t mid = lo + (hi - lo) / 2;
    if (arcs[mid].t < key) {
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
  if (edges > SIZE_MAX / (2 * sizeof(Arc))) {
    return GSHGEOM_ERR_OOM;
  }
  Arc * buf = gcu_allocator_malloc(alloc, 2 * (edges ? edges : 1) * sizeof *buf); /* the edges, and the sorts' room */
  unsigned char * used = gcu_allocator_calloc(alloc, edges ? edges : 1, 1);
  size_t * skip = NULL;
  GSHGEOM_Flat_Path * result = NULL;
  GSHGEOM_Result status = GSHGEOM_OK;
  if (!buf || !used) {
    status = GSHGEOM_ERR_OOM;
    goto done;
  }
  Arc * pairs = buf;
  Arc * room = buf + (edges ? edges : 1);
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
      uint64_t ka = point_key(a.x, a.y), kb = point_key(b.x, b.y);
      if (ka == kb) {
        continue;
      }
      Arc * e = &pairs[n++];
      e->forward = ka < kb;
      e->t = e->forward ? ka : kb;
      e->h = e->forward ? kb : ka;
    }
  }
  pairs = sort_arcs(pairs, pairs == buf ? room : buf, n);
  Arc * arcs = pairs == buf ? room : buf; /* the other half: the remaining edges are no more than the pairs */
  /* Each run of one edge between two points: as many go one way as the other cancel, and the rest stay as they were. */
  size_t m = 0;
  for (size_t i = 0; i < n;) {
    size_t j = i;
    long net = 0;
    while (j < n && pairs[j].t == pairs[i].t && pairs[j].h == pairs[i].h) {
      net += pairs[j].forward ? 1 : -1;
      j++;
    }
    for (long r = net < 0 ? -net : net; r > 0; r--) {
      Arc * a = &arcs[m++];
      a->t = net > 0 ? pairs[i].t : pairs[i].h;
      a->h = net > 0 ? pairs[i].h : pairs[i].t;
      a->forward = 1;
    }
    i = j;
  }
  arcs = sort_arcs(arcs, pairs, m);
  skip = gcu_allocator_malloc(alloc, (m + 1) * sizeof *skip);
  if (!skip) {
    status = GSHGEOM_ERR_OOM;
    goto done;
  }
  for (size_t i = 0; i <= m; i++) {
    skip[i] = i;
  }
  status = gshgeom_flat_path_new(alloc, limits, &result);
  if (status != GSHGEOM_OK) {
    goto done;
  }
  for (size_t i = 0; i < m; i++) {
    if (used[i]) {
      continue;
    }
    const uint64_t start = arcs[i].t;
    status = gshgeom_flat_path_begin(result, key_point(start));
    if (status != GSHGEOM_OK) {
      goto done;
    }
    used[i] = 1;
    skip[i] = i + 1;
    uint64_t head = arcs[i].h;
    while (head != start) {
      status = gshgeom_flat_path_add(result, key_point(head));
      if (status != GSHGEOM_OK) {
        goto done;
      }
      size_t j = first_unused(skip, lower_bound(arcs, m, head));
      if (j >= m || arcs[j].t != head) {
        status = GSHGEOM_ERR_INTERNAL; /* a point with more edges entering than leaving: not a union of closed subpaths */
        goto done;
      }
      used[j] = 1;
      skip[j] = j + 1;
      head = arcs[j].h;
    }
    gshgeom_flat_path_close_last(result);
  }
done:
  gcu_allocator_free(alloc, buf);
  gcu_allocator_free(alloc, used);
  gcu_allocator_free(alloc, skip);
  if (status != GSHGEOM_OK) {
    gshgeom_flat_path_destroy(result);
    return status;
  }
  *out = result;
  return GSHGEOM_OK;
}
