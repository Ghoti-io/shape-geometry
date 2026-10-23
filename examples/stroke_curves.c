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
 * Write a stroke as an SVG fill (stroke_curves.h): the case of a writer that meets a style field its format cannot
 * express (here inside alignment) and expands the stroke instead.  The map carries the scale from the model to the
 * document's pixels; the outline comes back in those pixels as lines and cubics, which is path data as it stands, and the
 * clip part (the closed subpaths the stroke is aligned against) as a clip path.  A stroke written this way stays exact to a
 * thousandth of a pixel at any zoom, where the polygon a rasterizer wants would show its facets.  The example checks
 * what it shows: the path data is the outline's, every `C` is a cubic, and the document is a fill.  It exits non-zero
 * otherwise.
 */

#include <ghoti.io/shape-geometry/shape-geometry.h>

#include <stdio.h>
#include <string.h>

/** Append the path data of @p path (moves, lines, cubics and closes) to @p out; the number of cubics is added to *cubics. */
static int path_data(const GSHGEOM_Path * path, char * out, size_t cap, int * cubics) {
  size_t len = 0;
  for (size_t i = 0; i < gshgeom_path_element_count(path); i++) {
    GSHGEOM_Element e;
    gshgeom_path_element(path, i, &e);
    const GSHGEOM_Point * p = e.curve.p;
    int n = 0;
    switch (e.verb) {
    case GSHGEOM_VERB_MOVE:
      n = snprintf(out + len, cap - len, "M%.9g %.9g", p[0].x, p[0].y);
      break;
    case GSHGEOM_VERB_LINE:
      n = snprintf(out + len, cap - len, "L%.9g %.9g", p[1].x, p[1].y);
      break;
    case GSHGEOM_VERB_CUBIC:
      n = snprintf(out + len, cap - len, "C%.9g %.9g %.9g %.9g %.9g %.9g", p[1].x, p[1].y, p[2].x, p[2].y, p[3].x, p[3].y);
      (*cubics)++;
      break;
    case GSHGEOM_VERB_CLOSE:
      n = snprintf(out + len, cap - len, "Z");
      break;
    default:
      return 0; /* form CUBIC holds nothing else */
    }
    if (n < 0 || (size_t)n >= cap - len) {
      return 0;
    }
    len += (size_t)n;
  }
  return 1;
}

int main(void) {
  GSHGEOM_Path * path = NULL;
  GSHGEOM_Point a = {20, 20, 0}, b = {100, 20, 0}, c = {100, 70, 0}, d = {20, 70, 0};
  if (gshgeom_path_create(NULL, NULL, 2, &path) != GSHGEOM_OK || gshgeom_path_move_to(path, a) != GSHGEOM_OK ||
      gshgeom_path_cubic_to(path, b, c, d) != GSHGEOM_OK || gshgeom_path_close(path) != GSHGEOM_OK) {
    fprintf(stderr, "stroke_curves: could not build the path\n");
    return 1;
  }
  GSHGEOM_Stroke stroke;
  gshgeom_stroke_init(&stroke);
  stroke.width = 6.0;
  stroke.join = GSHGEOM_STROKE_JOIN_ROUND;
  stroke.alignment = GSHGEOM_STROKE_ALIGN_INSIDE; /* SVG 1.1 has no inside stroke: the writer expands it */
  GSHGEOM_Map map;
  gshgeom_map_scale(2, 1.5, 1.5, 1.0, &map); /* the document's pixels are 1.5 model units */
  GSHGEOM_Stroke_Curves * curves = NULL;
  GSHGEOM_Result r = gshgeom_stroke_curves_expand(path, &stroke, &map, GSHGEOM_STROKE_CURVES_TOLERANCE,
      GSHGEOM_STROKE_CURVES_CUBIC, NULL, NULL, &curves);
  if (r != GSHGEOM_OK || gshgeom_stroke_curves_status(curves) != GSHGEOM_STROKE_DRAWN) {
    fprintf(stderr, "stroke_curves: the expansion failed: %s\n", gshgeom_result_string(r));
    return 1;
  }
  static char aligned[1 << 18], clip[1 << 12];
  int cubics = 0, clip_cubics = 0;
  const GSHGEOM_Path * pa = gshgeom_stroke_curves_aligned(curves);
  const GSHGEOM_Path * pc = gshgeom_stroke_curves_clip(curves);
  /* the clip part keeps the input's own elements (here a cubic and a close), which are path data as well */
  /* every subpath here is closed, so with inside alignment the body is empty and the aligned part is the stroke */
  int ok = pa && pc && path_data(pa, aligned, sizeof aligned, &cubics) && path_data(pc, clip, sizeof clip, &clip_cubics);
  ok = ok && cubics > 4 && clip_cubics == 1 && aligned[0] == 'M' && gshgeom_stroke_curves_clip_complement(curves) == 0 &&
      gshgeom_path_element_count(gshgeom_stroke_curves_body(curves)) == 0;
  if (ok) {
    /* the document: the aligned part clipped to the shape (inside alignment), under the nonzero rule */
    printf("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"200\" height=\"140\">\n"
           "<clipPath id=\"c\"><path clip-rule=\"nonzero\" d=\"%s\"/></clipPath>\n"
           "<path clip-path=\"url(#c)\" fill-rule=\"nonzero\" d=\"%s\"/>\n</svg>\n",
        clip, aligned);
  }
  printf("stroke_curves: %d cubics, hash %016llx\n", cubics, (unsigned long long)gshgeom_stroke_curves_hash(curves));
  gshgeom_stroke_curves_destroy(curves);
  gshgeom_path_destroy(path);
  if (!ok) {
    fprintf(stderr, "stroke_curves: the outline is not path data of lines and cubics\n");
    return 1;
  }
  return 0;
}
