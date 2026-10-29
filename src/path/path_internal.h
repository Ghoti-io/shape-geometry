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
 * The model path's representation, shared by path.c and arc.c (internal).
 *
 * Three parallel arrays: the elements (a verb and where its points and weights start), the
 * points (always three doubles, z zero in 2D) and the weights.  Every element stores its end
 * point as its last point, and a close stores the subpath's start, so the start of element i
 * is the end of element i - 1.  An arc stores `center`, `u`, `v`, `end` (four points) and
 * `t0`, `sweep` (two weights); its unit and flags live in the element.
 */

#ifndef GHOTI_IO_GSHGEOM_PATH_INTERNAL_H
#define GHOTI_IO_GSHGEOM_PATH_INTERNAL_H

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/shape-geometry/arc.h>
#include <ghoti.io/shape-geometry/path.h>

typedef struct Elem {
  uint8_t verb;
  uint8_t unit;  /* an arc's angle unit */
  uint8_t flags; /* an arc's flags */
  size_t point;  /* index of its first stored point */
  size_t weight; /* index of its first weight (rational verbs, and an arc's t0 and sweep) */
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

/** Points stored with a verb (the control points after the start, and the end), and weights. */
size_t gshgeom_path_verb_points(int verb);
size_t gshgeom_path_verb_weights(int verb);

/**
 * Append one element whose points and weights are already valid (no checking beyond the element limit and
 * allocation); a failure changes nothing but capacity.  @p unit and @p flags matter for an arc only.
 */
GSHGEOM_Result gshgeom_path_append_raw(GSHGEOM_Path * path, int verb, unsigned unit,
    unsigned flags, const GSHGEOM_Point * pts, const double * weights);

/** The state to return to when a constructor that adds several elements fails part way. */
typedef struct PathMark {
  size_t el_n, pt_n, wt_n;
  int has_current;
  GSHGEOM_Point current, start;
} PathMark;

static inline void path_mark(const GSHGEOM_Path * p, PathMark * m) {
  m->el_n = p->el_n;
  m->pt_n = p->pt_n;
  m->wt_n = p->wt_n;
  m->has_current = p->has_current;
  m->current = p->current;
  m->start = p->start;
}

static inline void path_rollback(GSHGEOM_Path * p, const PathMark * m) {
  p->el_n = m->el_n;
  p->pt_n = m->pt_n;
  p->wt_n = m->wt_n;
  p->has_current = m->has_current;
  p->current = m->current;
  p->start = m->start;
}

/** The start of element @p i (the end of the element before it; the point itself for a move). */
GSHGEOM_Point gshgeom_path_element_start(const GSHGEOM_Path * path, size_t i);

/** The current point (valid when ::GSHGEOM_Path has one). */
static inline GSHGEOM_Point path_current(const GSHGEOM_Path * p) {
  return p->current;
}

#endif /* GHOTI_IO_GSHGEOM_PATH_INTERNAL_H */
