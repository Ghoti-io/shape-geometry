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
 * Expand a stroke into its outline (stroke.h): a polyline and a curve, round
 * joins, square caps, a width of 3 pixels under a map that doubles everything.
 * The outline is 26.6 polygons whose nonzero fill is the region the stroke
 * covers; a raster sink fills it like any path.  It then dashes a straight line
 * (a dash of 20 and a gap of 10) and lowers a DXF linetype.  The example checks
 * what it shows: the outline is drawn, every polygon is closed, and a point on the
 * centerline is inside while a point far from it is not.  It exits non-zero
 * otherwise.
 */

#include <ghoti.io/shape-geometry/shape-geometry.h>

#include <stdio.h>

/** The winding number of the outline around (x, y) pixels: the sum over its closed polygons. */
static int winding(const GSHGEOM_Flat_Path * fp, double x, double y) {
  const GSHGEOM_Device_Point * pts = gshgeom_flat_path_points(fp);
  double px = x * 64.0, py = y * 64.0;
  int w = 0;
  for (size_t s = 0; s < gshgeom_flat_path_subpath_count(fp); s++) {
    GSHGEOM_Flat_Subpath sp;
    gshgeom_flat_path_subpath(fp, s, &sp);
    for (size_t i = 0; i < sp.count; i++) {
      GSHGEOM_Device_Point a = pts[sp.first + i], b = pts[sp.first + (i + 1) % sp.count];
      double c = (double)(b.x - a.x) * (py - a.y) - (px - a.x) * (double)(b.y - a.y);
      if (a.y <= py && b.y > py && c > 0) {
        w++;
      } else if (a.y > py && b.y <= py && c < 0) {
        w--;
      }
    }
  }
  return w;
}

int main(void) {
  GSHGEOM_Path * path = NULL;
  GSHGEOM_Point a = {10, 10, 0}, b = {60, 10, 0}, c = {80, 10, 0}, d = {80, 40, 0}, e = {40, 50, 0};
  if (gshgeom_path_create(NULL, NULL, 2, &path) != GSHGEOM_OK || gshgeom_path_move_to(path, a) != GSHGEOM_OK ||
      gshgeom_path_line_to(path, b) != GSHGEOM_OK || gshgeom_path_quadratic_to(path, c, d) != GSHGEOM_OK ||
      gshgeom_path_line_to(path, e) != GSHGEOM_OK) {
    fprintf(stderr, "stroke: could not build the path\n");
    return 1;
  }
  GSHGEOM_Stroke stroke;
  gshgeom_stroke_init(&stroke);
  stroke.width = 1.5; /* local width: 3 pixels after the map */
  stroke.join = GSHGEOM_STROKE_JOIN_ROUND;
  stroke.cap = GSHGEOM_STROKE_CAP_SQUARE;
  GSHGEOM_Map map;
  gshgeom_map_scale(2, 2.0, 2.0, 1.0, &map);
  GSHGEOM_Stroke_Outline * outline = NULL;
  GSHGEOM_Result r = gshgeom_stroke_expand(path, &stroke, &map, NULL, NULL, &outline);
  if (r != GSHGEOM_OK || gshgeom_stroke_outline_status(outline) != GSHGEOM_STROKE_DRAWN) {
    fprintf(stderr, "stroke: the expansion failed: %s\n", gshgeom_result_string(r));
    return 1;
  }
  const GSHGEOM_Flat_Path * body = gshgeom_stroke_outline_body(outline);
  int ok = gshgeom_flat_path_subpath_count(body) > 0;
  for (size_t s = 0; s < gshgeom_flat_path_subpath_count(body); s++) {
    GSHGEOM_Flat_Subpath sp;
    gshgeom_flat_path_subpath(body, s, &sp);
    ok = ok && sp.closed && sp.count >= 3;
  }
  /* (35, 10) model is (70, 20) pixels: on the centerline; (100, 100) pixels is far from it */
  ok = ok && winding(body, 70.0, 20.0) != 0 && winding(body, 100.0, 100.0) == 0;
  printf("stroke: %zu polygons, %zu points, hash %016llx\n", gshgeom_flat_path_subpath_count(body),
      gshgeom_flat_path_point_count(body), (unsigned long long)gshgeom_stroke_outline_hash(outline));
  gshgeom_stroke_outline_destroy(outline);
  gshgeom_path_destroy(path);

  /* a dashed line: dashes of 20 model units and gaps of 10, from (0, 0) to (100, 0); with the identity map the dash
   * [30, 50] covers the pixel (40, 0) and the gap [50, 60] does not cover (55, 0) */
  GSHGEOM_Path * line = NULL;
  GSHGEOM_Point o = {0, 0, 0}, far_end = {100, 0, 0};
  if (gshgeom_path_create(NULL, NULL, 2, &line) != GSHGEOM_OK || gshgeom_path_move_to(line, o) != GSHGEOM_OK ||
      gshgeom_path_line_to(line, far_end) != GSHGEOM_OK) {
    fprintf(stderr, "stroke: could not build the line\n");
    return 1;
  }
  double dashes[2] = {20.0, 10.0};
  gshgeom_stroke_init(&stroke);
  stroke.width = 4.0;
  stroke.dashes = dashes;
  stroke.dash_count = 2;
  GSHGEOM_Map identity;
  gshgeom_map_identity(2, &identity);
  outline = NULL;
  r = gshgeom_stroke_expand(line, &stroke, &identity, NULL, NULL, &outline);
  if (r != GSHGEOM_OK) {
    fprintf(stderr, "stroke: the dashed expansion failed: %s\n", gshgeom_result_string(r));
    return 1;
  }
  body = gshgeom_stroke_outline_body(outline);
  int dashed_ok = winding(body, 40.0, 0.0) != 0 && winding(body, 55.0, 0.0) == 0 && winding(body, 65.0, 0.0) != 0;
  printf("stroke: dashed, %zu polygons, hash %016llx\n", gshgeom_flat_path_subpath_count(body),
      (unsigned long long)gshgeom_stroke_outline_hash(outline));
  gshgeom_stroke_outline_destroy(outline);
  gshgeom_path_destroy(line);

  /* a DXF linetype of a dash 0.5, a space 0.25, a dot and a space 0.25, scaled by 8: the dash array and offset */
  const double dxf[4] = {0.5, -0.25, 0.0, -0.25};
  double lowered[8];
  size_t lowered_count = 0;
  double lowered_offset = 0.0;
  r = gshgeom_stroke_dashes_from_dxf(dxf, 4, 8.0, lowered, 8, &lowered_count, &lowered_offset);
  dashed_ok = dashed_ok && r == GSHGEOM_OK && lowered_count == 4 && lowered[0] == 4.0 && lowered[1] == 2.0 &&
      lowered[2] == 0.0 && lowered[3] == 2.0 && lowered_offset == 0.0;
  if (!ok) {
    fprintf(stderr, "stroke: the outline does not cover what the stroke covers\n");
    return 1;
  }
  if (!dashed_ok) {
    fprintf(stderr, "stroke: the dashed outline or the lowering is wrong\n");
    return 1;
  }
  return 0;
}
