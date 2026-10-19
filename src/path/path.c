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
 * The model path: builder, readback, bounds, hash and dump (path.h).
 *
 * Three parallel arrays: the elements (a verb and where its points and
 * weights start), the points (always three doubles, z zero in 2D) and the
 * weights.  Every element stores its end point, and a close stores the
 * subpath's start, so the start of element i is the end of element i - 1.
 * Every segment is validated as a ::GSHGEOM_Curve before it is stored, which
 * is what lets readback and bounds trust what they find.
 */

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/shape-geometry/hash.h>
#include <ghoti.io/shape-geometry/path.h>

#include <stdlib.h>
#include <string.h>

typedef struct Elem {
  uint8_t verb;
  size_t point;  /* index of its first stored point */
  size_t weight; /* index of its first weight (rational verbs) */
} Elem;

struct GSHGEOM_Path {
  GSHGEOM_Allocator alloc;
  size_t max_elements;
  int dim;
  Elem * el;
  size_t el_n, el_cap;
  GSHGEOM_Point * pt;
  size_t pt_n, pt_cap;
  double * wt;
  size_t wt_n, wt_cap;
  int has_current; /* a move has been made */
  GSHGEOM_Point current;
  GSHGEOM_Point start;
};

/* Points stored with a verb (the control points after the start, and the end),
 * and weights. */
static size_t verb_points(int verb) {
  switch (verb) {
  case GSHGEOM_VERB_MOVE:
  case GSHGEOM_VERB_LINE:
  case GSHGEOM_VERB_CLOSE:
    return 1;
  case GSHGEOM_VERB_QUADRATIC:
  case GSHGEOM_VERB_RATIONAL_QUADRATIC:
    return 2;
  case GSHGEOM_VERB_CUBIC:
  case GSHGEOM_VERB_RATIONAL_CUBIC:
    return 3;
  default:
    return 0;
  }
}

static size_t verb_weights(int verb) {
  return verb == GSHGEOM_VERB_RATIONAL_QUADRATIC ? 3
      : verb == GSHGEOM_VERB_RATIONAL_CUBIC     ? 4
                                                : 0;
}

GSHGEOM_Result gshgeom_path_create(const GSHGEOM_Allocator * alloc,
    const GSHGEOM_Limits * limits, int dim, GSHGEOM_Path ** out) {
  if (!out || (dim != 2 && dim != 3)) {
    return GSHGEOM_ERR_INVALID;
  }
  if (!alloc) {
    alloc = gshgeom_allocator_default();
  }
  GSHGEOM_Path * path = gcu_allocator_calloc(alloc, 1, sizeof *path);
  if (!path) {
    return GSHGEOM_ERR_OOM;
  }
  path->alloc = *alloc;
  path->dim = dim;
  path->max_elements = (limits && limits->max_elements)
      ? limits->max_elements
      : GSHGEOM_LIMITS_DEFAULT_MAX_ELEMENTS;
  *out = path;
  return GSHGEOM_OK;
}

void gshgeom_path_destroy(GSHGEOM_Path * path) {
  if (!path) {
    return;
  }
  GSHGEOM_Allocator alloc = path->alloc;
  gcu_allocator_free(&alloc, path->el);
  gcu_allocator_free(&alloc, path->pt);
  gcu_allocator_free(&alloc, path->wt);
  gcu_allocator_free(&alloc, path);
}

int gshgeom_path_dimension(const GSHGEOM_Path * path) {
  return path ? path->dim : 0;
}

size_t gshgeom_path_element_count(const GSHGEOM_Path * path) {
  return path ? path->el_n : 0;
}

/* Make room for `need` items in an array.  Only the capacity changes. */
static GSHGEOM_Result reserve(const GSHGEOM_Allocator * alloc, void * array,
    size_t * cap, size_t need, size_t size, void ** grown_out) {
  *grown_out = array;
  if (need <= *cap) {
    return GSHGEOM_OK;
  }
  size_t grown = *cap ? *cap : 16;
  while (grown < need) {
    if (grown > SIZE_MAX / 2) {
      return GSHGEOM_ERR_OOM;
    }
    grown *= 2;
  }
  if (grown > SIZE_MAX / size) {
    return GSHGEOM_ERR_OOM;
  }
  void * bigger = gcu_allocator_realloc(alloc, array, grown * size);
  if (!bigger) {
    return GSHGEOM_ERR_OOM;
  }
  *grown_out = bigger;
  *cap = grown;
  return GSHGEOM_OK;
}

/* Add one element: the limit first, then every reservation, and only then the
 * stores, so a failure changes nothing but capacity. */
static GSHGEOM_Result append(GSHGEOM_Path * path, int verb,
    const GSHGEOM_Point * pts, const double * weights) {
  size_t np = verb_points(verb);
  size_t nw = verb_weights(verb);
  if (path->el_n >= path->max_elements) {
    return GSHGEOM_ERR_LIMIT;
  }
  void * grown;
  GSHGEOM_Result r = reserve(&path->alloc, path->el, &path->el_cap,
      path->el_n + 1, sizeof *path->el, &grown);
  if (r != GSHGEOM_OK) {
    return r;
  }
  path->el = grown;
  r = reserve(&path->alloc, path->pt, &path->pt_cap, path->pt_n + np,
      sizeof *path->pt, &grown);
  if (r != GSHGEOM_OK) {
    return r;
  }
  path->pt = grown;
  r = reserve(&path->alloc, path->wt, &path->wt_cap, path->wt_n + nw,
      sizeof *path->wt, &grown);
  if (r != GSHGEOM_OK) {
    return r;
  }
  path->wt = grown;
  Elem * e = &path->el[path->el_n];
  e->verb = (uint8_t)verb;
  e->point = path->pt_n;
  e->weight = path->wt_n;
  for (size_t i = 0; i < np; i++) {
    path->pt[path->pt_n + i] = pts[i];
  }
  for (size_t i = 0; i < nw; i++) {
    path->wt[path->wt_n + i] = weights[i];
  }
  path->el_n++;
  path->pt_n += np;
  path->wt_n += nw;
  path->current = pts[np - 1];
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_path_move_to(GSHGEOM_Path * path, GSHGEOM_Point p) {
  if (!path) {
    return GSHGEOM_ERR_INVALID;
  }
  /* A move is a point; checking it as a degenerate line checks it as the
   * curves are checked, in one place. */
  GSHGEOM_Curve c;
  GSHGEOM_Result r = gshgeom_curve_line(path->dim, p, p, &c);
  if (r != GSHGEOM_OK) {
    return r;
  }
  r = append(path, GSHGEOM_VERB_MOVE, &c.p[0], NULL);
  if (r == GSHGEOM_OK) {
    path->has_current = 1;
    path->start = c.p[0];
  }
  return r;
}

GSHGEOM_Result gshgeom_path_line_to(GSHGEOM_Path * path, GSHGEOM_Point p) {
  if (!path || !path->has_current) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Curve c;
  GSHGEOM_Result r = gshgeom_curve_line(path->dim, path->current, p, &c);
  return r != GSHGEOM_OK ? r : append(path, GSHGEOM_VERB_LINE, &c.p[1], NULL);
}

GSHGEOM_Result gshgeom_path_quadratic_to(
    GSHGEOM_Path * path, GSHGEOM_Point cp, GSHGEOM_Point p) {
  if (!path || !path->has_current) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Curve c;
  GSHGEOM_Result r = gshgeom_curve_quadratic(path->dim, path->current, cp, p, &c);
  return r != GSHGEOM_OK ? r
                         : append(path, GSHGEOM_VERB_QUADRATIC, &c.p[1], NULL);
}

GSHGEOM_Result gshgeom_path_cubic_to(GSHGEOM_Path * path, GSHGEOM_Point c1,
    GSHGEOM_Point c2, GSHGEOM_Point p) {
  if (!path || !path->has_current) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Curve c;
  GSHGEOM_Result r =
      gshgeom_curve_cubic(path->dim, path->current, c1, c2, p, &c);
  return r != GSHGEOM_OK ? r : append(path, GSHGEOM_VERB_CUBIC, &c.p[1], NULL);
}

GSHGEOM_Result gshgeom_path_rational_quadratic_to(GSHGEOM_Path * path,
    GSHGEOM_Point cp, GSHGEOM_Point p, const double w[3]) {
  if (!path || !path->has_current) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Curve c;
  GSHGEOM_Result r = gshgeom_curve_rational_quadratic(
      path->dim, path->current, cp, p, w, &c);
  return r != GSHGEOM_OK
      ? r
      : append(path, GSHGEOM_VERB_RATIONAL_QUADRATIC, &c.p[1], c.w);
}

GSHGEOM_Result gshgeom_path_rational_cubic_to(GSHGEOM_Path * path,
    GSHGEOM_Point c1, GSHGEOM_Point c2, GSHGEOM_Point p, const double w[4]) {
  if (!path || !path->has_current) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Curve c;
  GSHGEOM_Result r = gshgeom_curve_rational_cubic(
      path->dim, path->current, c1, c2, p, w, &c);
  return r != GSHGEOM_OK
      ? r
      : append(path, GSHGEOM_VERB_RATIONAL_CUBIC, &c.p[1], c.w);
}

GSHGEOM_Result gshgeom_path_close(GSHGEOM_Path * path) {
  if (!path || !path->has_current) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Point start = path->start;
  return append(path, GSHGEOM_VERB_CLOSE, &start, NULL);
}

GSHGEOM_Result gshgeom_path_element(
    const GSHGEOM_Path * path, size_t index, GSHGEOM_Element * out) {
  if (!path || !out || index >= path->el_n) {
    return GSHGEOM_ERR_INVALID;
  }
  const Elem * e = &path->el[index];
  size_t np = verb_points(e->verb);
  GSHGEOM_Element el = {0};
  el.verb = (GSHGEOM_Verb)e->verb;
  el.curve.dim = path->dim;

  /* The start is the end of the element before. */
  GSHGEOM_Point from = path->pt[e->point];
  if (index > 0) {
    const Elem * prev = &path->el[index - 1];
    from = path->pt[prev->point + verb_points(prev->verb) - 1];
  }
  GSHGEOM_Curve * c = &el.curve;
  c->p[0] = from;
  for (size_t i = 0; i < np; i++) {
    c->p[1 + i] = path->pt[e->point + i];
  }
  switch (e->verb) {
  case GSHGEOM_VERB_MOVE:
    c->kind = GSHGEOM_CURVE_LINE;
    c->p[0] = path->pt[e->point]; /* p[1] is the same point, from the loop */
    break;
  case GSHGEOM_VERB_LINE:
  case GSHGEOM_VERB_CLOSE:
    c->kind = GSHGEOM_CURVE_LINE;
    break;
  case GSHGEOM_VERB_QUADRATIC:
    c->kind = GSHGEOM_CURVE_QUADRATIC;
    break;
  case GSHGEOM_VERB_CUBIC:
    c->kind = GSHGEOM_CURVE_CUBIC;
    break;
  case GSHGEOM_VERB_RATIONAL_QUADRATIC:
    c->kind = GSHGEOM_CURVE_RATIONAL_QUADRATIC;
    break;
  case GSHGEOM_VERB_RATIONAL_CUBIC:
    c->kind = GSHGEOM_CURVE_RATIONAL_CUBIC;
    break;
  default:
    return GSHGEOM_ERR_INTERNAL;
  }
  for (size_t i = 0; i < verb_weights(e->verb); i++) {
    c->w[i] = path->wt[e->weight + i];
  }
  *out = el;
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_path_bounds(
    const GSHGEOM_Path * path, GSHGEOM_Bounds * out) {
  if (!path || !out) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Bounds box = {0};
  box.empty = 1;
  for (size_t i = 0; i < path->el_n; i++) {
    int verb = path->el[i].verb;
    if (verb == GSHGEOM_VERB_MOVE || verb == GSHGEOM_VERB_CLOSE) {
      continue;
    }
    GSHGEOM_Element el;
    GSHGEOM_Bounds b;
    GSHGEOM_Result r = gshgeom_path_element(path, i, &el);
    if (r == GSHGEOM_OK) {
      r = gshgeom_curve_bounds(&el.curve, &b);
    }
    if (r != GSHGEOM_OK) {
      /* Every stored segment was validated, so this is our own invariant. */
      return GSHGEOM_ERR_INTERNAL;
    }
    if (box.empty) {
      box = b;
    } else {
      if (b.min.x < box.min.x) box.min.x = b.min.x;
      if (b.min.y < box.min.y) box.min.y = b.min.y;
      if (b.min.z < box.min.z) box.min.z = b.min.z;
      if (b.max.x > box.max.x) box.max.x = b.max.x;
      if (b.max.y > box.max.y) box.max.y = b.max.y;
      if (b.max.z > box.max.z) box.max.z = b.max.z;
    }
  }
  *out = box;
  return GSHGEOM_OK;
}

/* ---------------------------------------------------------------------- */
/* Hash and dump                                                           */
/* ---------------------------------------------------------------------- */

/** The format version of the hash and the dump. */
#define PATH_FORMAT_VERSION 1u

static void hash_point(GSHGEOM_Hash * h, const GSHGEOM_Point * p, int dim) {
  gshgeom_hash_f64(h, p->x);
  gshgeom_hash_f64(h, p->y);
  if (dim == 3) {
    gshgeom_hash_f64(h, p->z);
  }
}

uint64_t gshgeom_path_hash(const GSHGEOM_Path * path) {
  GSHGEOM_Hash h;
  gshgeom_hash_init(&h);
  if (!path) {
    return gshgeom_hash_value(&h);
  }
  gshgeom_hash_u32(&h, PATH_FORMAT_VERSION);
  gshgeom_hash_u32(&h, (uint32_t)path->dim);
  gshgeom_hash_u64(&h, (uint64_t)path->el_n);
  for (size_t i = 0; i < path->el_n; i++) {
    const Elem * e = &path->el[i];
    gshgeom_hash_u8(&h, e->verb);
    if (e->verb == GSHGEOM_VERB_CLOSE) {
      continue;
    }
    for (size_t k = 0; k < verb_points(e->verb); k++) {
      hash_point(&h, &path->pt[e->point + k], path->dim);
    }
    for (size_t k = 0; k < verb_weights(e->verb); k++) {
      gshgeom_hash_f64(&h, path->wt[e->weight + k]);
    }
  }
  return gshgeom_hash_value(&h);
}

static int put_f64(FILE * out, double v) {
  return fprintf(out, " %016llx", (unsigned long long)gshgeom_f64_bits(v));
}

static int put_point(FILE * out, const GSHGEOM_Point * p, int dim) {
  int n = put_f64(out, p->x);
  if (n >= 0) n = put_f64(out, p->y);
  if (n >= 0 && dim == 3) n = put_f64(out, p->z);
  return n;
}

GSHGEOM_Result gshgeom_path_dump(const GSHGEOM_Path * path, FILE * out) {
  static const char * const names[GSHGEOM_VERB_COUNT] = {
      "M", "L", "Q", "C", "RQ", "RC", "Z"};
  if (!path || !out) {
    return GSHGEOM_ERR_INVALID;
  }
  /* The count is printed through an unsigned long long: the width of size_t
   * is a property of the host. */
  if (fprintf(out, "gshgeom.path %u dim=%d elements=%llu hash=%016llx\n",
          PATH_FORMAT_VERSION, path->dim, (unsigned long long)path->el_n,
          (unsigned long long)gshgeom_path_hash(path)) < 0) {
    return GSHGEOM_ERR_IO;
  }
  for (size_t i = 0; i < path->el_n; i++) {
    const Elem * e = &path->el[i];
    if (fputs(names[e->verb], out) < 0) {
      return GSHGEOM_ERR_IO;
    }
    if (e->verb != GSHGEOM_VERB_CLOSE) {
      for (size_t k = 0; k < verb_weights(e->verb); k++) {
        if (put_f64(out, path->wt[e->weight + k]) < 0) {
          return GSHGEOM_ERR_IO;
        }
      }
      for (size_t k = 0; k < verb_points(e->verb); k++) {
        if (put_point(out, &path->pt[e->point + k], path->dim) < 0) {
          return GSHGEOM_ERR_IO;
        }
      }
    }
    if (fputc('\n', out) == EOF) {
      return GSHGEOM_ERR_IO;
    }
  }
  /* A buffered stream reports a write error at the flush, and a stream's error flag is sticky:
   * a failed write whose return value was lost (an unbuffered stream's fprintf) still sets it. */
  return (fflush(out) == 0 && !ferror(out)) ? GSHGEOM_OK : GSHGEOM_ERR_IO;
}
