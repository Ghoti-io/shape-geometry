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
 * Named primitives: a rounded rectangle, an ellipse and the four sign cases of a Canvas `rect` as two flips, and an overlap of the
 * corner curves refused on the exact sum.
 *
 *     $ ./shapes
 *     rounded rectangle: 10 elements, starts at (10, 0), 4 arcs
 *     ...
 */

#include <ghoti.io/shape-geometry/shape-geometry.h>

#include <math.h>
#include <stdio.h>

static int check(GSHGEOM_Result r, const char * what) {
  if (r != GSHGEOM_OK) {
    fprintf(stderr, "shapes: %s: %s\n", what, gshgeom_result_string(r));
    return 1;
  }
  return 0;
}

int main(void) {
  int failed = 0;

  /* A 100 by 60 rectangle with different corners: it starts after the top-left radius and turns clockwise on a y-down screen. */
  GSHGEOM_Path * rr = NULL;
  const GSHGEOM_Corner corners[4] = {{10, 10}, {20, 10}, {5, 5}, {0, 0}};
  failed |= check(gshgeom_path_create(NULL, NULL, 2, &rr), "create");
  failed |= check(gshgeom_path_add_rounded_rect(rr, 0, 0, 100, 60, corners, 0, 0), "rounded rectangle");
  GSHGEOM_Element first;
  failed |= check(gshgeom_path_element(rr, 0, &first), "first element");
  size_t arcs = 0;
  for (size_t i = 0; i < gshgeom_path_element_count(rr); i++) {
    GSHGEOM_Verb verb;
    failed |= check(gshgeom_path_verb(rr, i, &verb), "verb");
    arcs += verb == GSHGEOM_VERB_ARC;
  }
  printf("rounded rectangle: %zu elements, starts at (%g, %g), %zu arcs (the square corner has none; its lines stay)\n",
      gshgeom_path_element_count(rr), first.curve.p[0].x, first.curve.p[0].y, arcs);
  failed |= gshgeom_path_element_count(rr) != 9 || arcs != 3 || first.curve.p[0].x != 10.0;

  /* The four sign cases of Canvas rect(x, y, w, h) as two flips: where the path starts. */
  for (int fy = 0; fy < 2; fy++) {
    for (int fx = 0; fx < 2; fx++) {
      GSHGEOM_Path * r = NULL;
      failed |= check(gshgeom_path_create(NULL, NULL, 2, &r), "create");
      failed |= check(gshgeom_path_add_rect(r, 10, 20, 100, 50, fx, fy), "rect");
      GSHGEOM_Element m;
      failed |= check(gshgeom_path_element(r, 0, &m), "move");
      printf("rect flip_x=%d flip_y=%d: starts at (%g, %g), %s\n", fx, fy, m.curve.p[0].x, m.curve.p[0].y,
          fx == fy ? "clockwise" : "counter-clockwise");
      failed |= m.curve.p[0].x != (fx ? 110.0 : 10.0) || m.curve.p[0].y != (fy ? 70.0 : 20.0);
      gshgeom_path_destroy(r);
    }
  }

  /* A circle is an ellipse with equal radii; its arcs are circular, so a toolpath or a writer can say so. */
  GSHGEOM_Path * e = NULL;
  GSHGEOM_Arc a;
  failed |= check(gshgeom_path_create(NULL, NULL, 2, &e), "create");
  failed |= check(gshgeom_path_add_ellipse(e, 50, 30, 20, 20, 0, 0), "circle");
  failed |= check(gshgeom_path_arc(e, 1, &a), "arc");
  printf("circle: %zu elements, four quarter arcs, circular %d, first arc ends at (%g, %g)\n", gshgeom_path_element_count(e),
      (a.flags & GSHGEOM_ARC_CIRCULAR) != 0u, a.end.x, a.end.y);
  failed |= gshgeom_path_element_count(e) != 6 || (a.flags & GSHGEOM_ARC_CIRCULAR) == 0u || a.end.x != 50.0 || a.end.y != 50.0;

  /* Radii that overlap are refused, decided on the exact sum: 0.1 + 0.9 is 1.0 in double and a little over one in the reals. */
  const GSHGEOM_Corner tight[4] = {{0.1, 0.1}, {0.9, 0.1}, {0, 0}, {0, 0}};
  const GSHGEOM_Result narrow = gshgeom_rounded_rect_check(1.0, 10.0, tight);
  const GSHGEOM_Result wider = gshgeom_rounded_rect_check(nextafter(1.0, 2.0), 10.0, tight);
  printf("radii 0.1 and 0.9 on a width of 1: %s; one ulp wider: %s\n", gshgeom_result_string(narrow), gshgeom_result_string(wider));
  failed |= narrow != GSHGEOM_ERR_INVALID || wider != GSHGEOM_OK;

  gshgeom_path_destroy(e);
  gshgeom_path_destroy(rr);
  if (failed) {
    fprintf(stderr, "shapes: FAILED\n");
  }
  return failed;
}
