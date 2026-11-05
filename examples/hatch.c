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
 * Hatch a region (hatch.h): a crosshatch of two families over a rectangle with a hole, one of them dashed, made as pieces (the
 * pen-down segments a plotter follows) and then expanded to the outline a raster sink fills.  The example checks what it
 * shows: the pieces are made and counted, a point on a line between two dashes is outside the outline and a point on a dash is
 * inside, and a point in the hole (with the outer-ring style) is not hatched.  It exits non-zero otherwise.
 */

#include <ghoti.io/shape-geometry/shape-geometry.h>

#include <stdio.h>

/** The winding number of the outline around (x, y) pixels. */
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
  GSHGEOM_Path * region = NULL;
  if (gshgeom_path_create(NULL, NULL, 2, &region) != GSHGEOM_OK || gshgeom_path_add_rect(region, 0, 0, 40, 20, 0, 0) != GSHGEOM_OK ||
      gshgeom_path_add_rect(region, 10, 5, 20, 10, 1, 0) != GSHGEOM_OK) {
    fprintf(stderr, "hatch: could not build the region\n");
    return 1;
  }
  /* A horizontal family, spacing 4, dashes 6 on and 2 off (a DXF pattern line); and a vertical solid one, spacing 8. */
  const double elements[2] = {6.0, -2.0};
  GSHGEOM_Hatch_Family families[2];
  families[0] = (GSHGEOM_Hatch_Family){0.0, 0.0, 0.0, 0.0, 4.0, 2, elements};
  families[1] = (GSHGEOM_Hatch_Family){90.0, 0.0, 0.0, 0.0, 8.0, 0, NULL};
  GSHGEOM_Hatch hatch;
  hatch.angle_unit = GSHGEOM_ANGLE_DEGREES;
  hatch.nesting = GSHGEOM_HATCH_NESTING_OUTER; /* DXF's "outer" style: the hole stays empty */
  gshgeom_map_identity(2, &hatch.map);
  hatch.family_count = 2;
  hatch.families = families;
  GSHGEOM_Map chain;
  gshgeom_map_scale(2, 3.0, 3.0, 1.0, &chain);

  GSHGEOM_Hatch_Pieces * pieces = NULL;
  GSHGEOM_Result r = gshgeom_hatch_pieces(region, GSHGEOM_FILL_NONZERO, &hatch, &chain, 0, NULL, NULL, &pieces);
  if (r != GSHGEOM_OK || gshgeom_hatch_pieces_count(pieces) == 0) {
    fprintf(stderr, "hatch: no pieces: %s\n", gshgeom_result_string(r));
    return 1;
  }
  printf("hatch: %llu lines of the use, %zu pieces\n", (unsigned long long)gshgeom_hatch_pieces_lines(pieces),
         gshgeom_hatch_pieces_count(pieces));

  GSHGEOM_Stroke stroke;
  gshgeom_stroke_init(&stroke);
  stroke.width = 0.5; /* local width: 1.5 pixels after the map */
  GSHGEOM_Stroke_Outline * outline = NULL;
  r = gshgeom_hatch_expand_pieces(pieces, &stroke, &chain, GSHGEOM_FLATTEN_TOLERANCE, NULL, NULL, &outline);
  if (r != GSHGEOM_OK || gshgeom_stroke_outline_status(outline) != GSHGEOM_STROKE_DRAWN) {
    fprintf(stderr, "hatch: the expansion failed: %s\n", gshgeom_result_string(r));
    return 1;
  }
  const GSHGEOM_Flat_Path * body = gshgeom_stroke_outline_body(outline);
  /* the line y = 4 is dashed: on from x = 0 to 6, off from 6 to 8: in pixels (3x) a dash at (9, 12), a gap at (21, 12) */
  int on_dash = winding(body, 9.0, 12.0) != 0;
  int in_gap = winding(body, 21.0, 12.0) == 0;
  /* the vertical line x = 16 (48 pixels) crosses the hole: it is drawn above and below it, and not through it */
  int above_hole = winding(body, 48.0, 6.0) != 0;
  int in_hole = winding(body, 48.0, 30.0) == 0;
  gshgeom_stroke_outline_destroy(outline);
  gshgeom_hatch_pieces_destroy(pieces);
  gshgeom_path_destroy(region);
  if (!on_dash || !in_gap || !above_hole || !in_hole) {
    fprintf(stderr, "hatch: the outline is not what the definition says (dash %d, gap %d, line %d, hole %d)\n", on_dash, in_gap, above_hole, in_hole);
    return 1;
  }
  printf("hatch: a dash is drawn, a gap is not, and the hole is empty\n");
  return 0;
}
