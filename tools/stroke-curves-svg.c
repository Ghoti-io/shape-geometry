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
 * The documents of the zoom gate (stroke_curves.h, SG43; tools/resvg-stroke/zoom.py).  For each case it writes three
 * SVG files into a directory: NAME.stroke.svg, the stroke as a stroke (join, cap, dashes, a transform for the elliptical
 * pen) for the external renderer to stroke itself; NAME.fill.svg, the same stroke expanded by gshgeom_stroke_curves_expand
 * in form CUBIC at the default tolerance and written as the fill of its body; and NAME.poly.svg, the control: the polygon
 * of gshgeom_stroke_expand at the default tolerance of the polyline, written as a fill.  The files are not scaled: the
 * renderer zooms them.  The path data of an outline is the example's writer (examples/stroke_curves.c).
 *
 *   usage: stroke-curves-svg DIR
 */

#include <ghoti.io/shape-geometry/shape-geometry.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BOX 12.0

typedef struct Case {
  const char * name;
  const char * d;           /* the path as SVG path data, in model units */
  double width;
  int join, cap;
  double limit;
  const double * dashes;
  size_t dash_count;
  double dash_offset;
  double map[4];            /* a b c d of x' = a x + b y, y' = c x + d y */
  int curved;               /* the control must fail at 16x and 64x */
  const char * reference;   /* when not NULL, the SVG element that is the region exactly, drawn in place of the renderer's stroke */
} Case;

static const double dash_a[2] = {1.6, 1.1};
static const double dash_b[2] = {0.0, 1.3};

static const Case cases[] = {
    {"scurve_round", "M 1 6 C 4 0 8 12 11 6", 1.6, GSHGEOM_STROKE_JOIN_ROUND, GSHGEOM_STROKE_CAP_ROUND, 4, NULL, 0, 0,
        {1, 0, 0, 1}, 1, NULL},
    {"polyline_miter", "M 1 10 L 6 1.5 L 11 10", 1.4, GSHGEOM_STROKE_JOIN_MITER, GSHGEOM_STROKE_CAP_BUTT, 4, NULL, 0, 0,
        {1, 0, 0, 1}, 0, NULL},
    {"polyline_bevel", "M 1 10 L 6 1.5 L 11 10", 1.4, GSHGEOM_STROKE_JOIN_BEVEL, GSHGEOM_STROKE_CAP_SQUARE, 4, NULL, 0, 0,
        {1, 0, 0, 1}, 0, NULL},
    {"circle_fold", "M 7.2 6 A 1.2 1.2 0 0 1 6 7.2 A 1.2 1.2 0 0 1 4.8 6 A 1.2 1.2 0 0 1 6 4.8 A 1.2 1.2 0 0 1 7.2 6 Z", 5.0,
        GSHGEOM_STROKE_JOIN_ROUND, GSHGEOM_STROKE_CAP_BUTT, 4, NULL, 0, 0, {1, 0, 0, 1}, 1,
        "<circle cx=\"6\" cy=\"6\" r=\"3.7\" fill=\"#000000\"/>"},
    {"dashed_cubic", "M 1 10 C 3 0 9 0 11 10", 1.2, GSHGEOM_STROKE_JOIN_ROUND, GSHGEOM_STROKE_CAP_ROUND, 4, dash_a, 2, 0.5,
        {1, 0, 0, 1}, 1, NULL},
    {"dotted_cubic", "M 1 10 C 3 0 9 0 11 10", 1.0, GSHGEOM_STROKE_JOIN_ROUND, GSHGEOM_STROKE_CAP_ROUND, 4, dash_b, 2, 0,
        {1, 0, 0, 1}, 1, NULL},
    {"elliptic_pen", "M 1 1.5 C 4 -0.5 8 3.5 11 1.5", 1.0, GSHGEOM_STROKE_JOIN_ROUND, GSHGEOM_STROKE_CAP_ROUND, 4, NULL, 0, 0,
        {1, 0, 0, 3}, 1, NULL},
    {"cubic_loop", "M 3 9 C 14 0 -2 0 9 9", 1.5, GSHGEOM_STROKE_JOIN_ROUND, GSHGEOM_STROKE_CAP_BUTT, 4, NULL, 0, 0,
        {1, 0, 0, 1}, 1, NULL},
    {"curve_miter", "M 1 2 L 6 2 Q 11 2 11 7 L 6 11", 1.2, GSHGEOM_STROKE_JOIN_MITER, GSHGEOM_STROKE_CAP_BUTT, 4, NULL, 0, 0,
        {1, 0, 0, 1}, 1, NULL},
};

/** Parse the subset of path data the cases use (M L Q C A Z) into a model path; A is written as four rational arcs. */
static int parse(const char * d, GSHGEOM_Path ** out) {
  GSHGEOM_Path * p = NULL;
  if (gshgeom_path_create(NULL, NULL, 2, &p) != GSHGEOM_OK) {
    return 0;
  }
  const char * s = d;
  char verb = 0;
  double cx = 0, cy = 0;
  while (*s) {
    while (*s == ' ') {
      s++;
    }
    if (!*s) {
      break;
    }
    if ((*s >= 'A' && *s <= 'Z')) {
      verb = *s++;
    }
    double v[7];
    int need = verb == 'M' || verb == 'L' ? 2 : verb == 'Q' ? 4 : verb == 'C' ? 6 : verb == 'A' ? 7 : 0;
    for (int i = 0; i < need; i++) {
      char * end;
      v[i] = strtod(s, &end);
      s = end;
    }
    GSHGEOM_Point a = {v[0], v[1], 0}, b = {v[2], v[3], 0}, c = {v[4], v[5], 0};
    GSHGEOM_Result r = GSHGEOM_OK;
    switch (verb) {
    case 'M':
      r = gshgeom_path_move_to(p, a);
      cx = v[0];
      cy = v[1];
      break;
    case 'L':
      r = gshgeom_path_line_to(p, a);
      cx = v[0];
      cy = v[1];
      break;
    case 'Q':
      r = gshgeom_path_quadratic_to(p, a, b);
      cx = v[2];
      cy = v[3];
      break;
    case 'C':
      r = gshgeom_path_cubic_to(p, a, b, c);
      cx = v[4];
      cy = v[5];
      break;
    case 'A': {
      /* a quarter circle of radius v[0] from the current point to (v[5], v[6]), sweep positive: the control point is the
       * corner of the two tangents */
      double ex = v[5], ey = v[6];
      double mx = cx + (ex - cx) / 2, my = cy + (ey - cy) / 2;
      /* the centre is on the perpendicular bisector at distance rho from the midpoint; for a quarter arc of radius r the
       * chord is r sqrt 2 and the corner lies r from the chord's midpoint away from the centre: the centre is the
       * other corner of the square on the chord */
      double dx = ex - cx, dy = ey - cy;
      double zx = mx - dy / 2, zy = my + dx / 2; /* the centre on the left of the chord direction (sweep 1) */
      double qx = mx + dy / 2, qy = my - dx / 2; /* the tangents' corner, on the other side */
      (void)zx;
      (void)zy;
      double w[3] = {1, 0.70710678118654752440, 1};
      GSHGEOM_Point cp = {qx, qy, 0}, ep = {ex, ey, 0};
      r = gshgeom_path_rational_quadratic_to(p, cp, ep, w);
      cx = ex;
      cy = ey;
      break;
    }
    case 'Z':
      r = gshgeom_path_close(p);
      break;
    default:
      return 0;
    }
    if (r != GSHGEOM_OK) {
      return 0;
    }
  }
  *out = p;
  return 1;
}

/** SVG path data of a path of lines and cubics, 10 significant digits. */
static void write_path(FILE * f, const GSHGEOM_Path * p) {
  for (size_t i = 0; i < gshgeom_path_element_count(p); i++) {
    GSHGEOM_Element e;
    gshgeom_path_element(p, i, &e);
    const GSHGEOM_Point * q = e.curve.p;
    switch (e.verb) {
    case GSHGEOM_VERB_MOVE:
      fprintf(f, "M%.10g %.10g", q[0].x, q[0].y);
      break;
    case GSHGEOM_VERB_LINE:
      fprintf(f, "L%.10g %.10g", q[1].x, q[1].y);
      break;
    case GSHGEOM_VERB_CUBIC:
      fprintf(f, "C%.10g %.10g %.10g %.10g %.10g %.10g", q[1].x, q[1].y, q[2].x, q[2].y, q[3].x, q[3].y);
      break;
    case GSHGEOM_VERB_CLOSE:
      fputc('Z', f);
      break;
    default:
      fprintf(stderr, "stroke-curves-svg: an element the writer does not write\n");
      exit(1);
    }
  }
}

static FILE * open_doc(const char * dir, const char * name, const char * kind) {
  char path[512];
  snprintf(path, sizeof path, "%s/%s.%s.svg", dir, name, kind);
  FILE * f = fopen(path, "w");
  if (!f) {
    fprintf(stderr, "stroke-curves-svg: cannot write %s\n", path);
    exit(1);
  }
  fprintf(f, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%g\" height=\"%g\" viewBox=\"0 0 %g %g\" "
             "shape-rendering=\"crispEdges\">\n",
      BOX, BOX, BOX, BOX);
  return f;
}

static const char * const JOIN[3] = {"miter", "round", "bevel"};
static const char * const CAP[3] = {"butt", "round", "square"};

int main(int argc, char ** argv) {
  if (argc != 2) {
    fprintf(stderr, "usage: stroke-curves-svg DIR\n");
    return 2;
  }
  for (size_t ci = 0; ci < sizeof cases / sizeof cases[0]; ci++) {
    const Case * c = &cases[ci];
    GSHGEOM_Path * path = NULL;
    if (!parse(c->d, &path)) {
      fprintf(stderr, "stroke-curves-svg: cannot parse %s\n", c->name);
      return 1;
    }
    /* the stroke, for the renderer to stroke */
    FILE * f = open_doc(argv[1], c->name, "stroke");
    if (c->reference) {
      /* a stroke wider than twice its circle's radius is the disc of radius r + h, which the renderer's stroker does not
       * paint (it leaves the middle out): the exact region stands in for its stroke */
      fprintf(f, "%s\n</svg>\n", c->reference);
    } else {
      /* a stroke join that is miter-clip or arcs would need svg2: the cases use none */
      fprintf(f, "<path d=\"%s\" transform=\"matrix(%g %g %g %g 0 0)\" fill=\"none\" stroke=\"#000000\" "
                 "stroke-width=\"%g\" stroke-linejoin=\"%s\" stroke-linecap=\"%s\" stroke-miterlimit=\"%g\"",
          c->d, c->map[0], c->map[2], c->map[1], c->map[3], c->width, JOIN[c->join], CAP[c->cap], c->limit);
      if (c->dash_count) {
        fprintf(f, " stroke-dasharray=\"");
        for (size_t i = 0; i < c->dash_count; i++) {
          fprintf(f, i ? " %g" : "%g", c->dashes[i]);
        }
        fprintf(f, "\" stroke-dashoffset=\"%g\"", c->dash_offset);
      }
      fprintf(f, "/>\n</svg>\n");
    }
    fclose(f);
    /* the expansions */
    GSHGEOM_Stroke st;
    gshgeom_stroke_init(&st);
    st.width = c->width;
    st.join = (uint32_t)c->join;
    st.cap = (uint32_t)c->cap;
    st.miter_limit = c->limit;
    st.dashes = c->dashes;
    st.dash_count = c->dash_count;
    st.dash_offset = c->dash_offset;
    double mm[4] = {c->map[0], c->map[1], c->map[2], c->map[3]}, tt[2] = {0, 0};
    GSHGEOM_Map map;
    if (gshgeom_map_from_matrix(2, mm, tt, &map) != GSHGEOM_OK) {
      return 1;
    }
    GSHGEOM_Stroke_Curves * cu = NULL;
    GSHGEOM_Result r = gshgeom_stroke_curves_expand(path, &st, &map, GSHGEOM_STROKE_CURVES_TOLERANCE,
        GSHGEOM_STROKE_CURVES_CUBIC, NULL, NULL, &cu);
    if (r != GSHGEOM_OK) {
      fprintf(stderr, "stroke-curves-svg: %s: the curve expansion failed: %s\n", c->name, gshgeom_result_string(r));
      return 1;
    }
    f = open_doc(argv[1], c->name, "fill");
    fprintf(f, "<path d=\"");
    write_path(f, gshgeom_stroke_curves_body(cu));
    fprintf(f, "\" fill=\"#000000\" fill-rule=\"nonzero\"/>\n</svg>\n");
    fclose(f);
    gshgeom_stroke_curves_destroy(cu);
    GSHGEOM_Stroke_Outline * so = NULL;
    r = gshgeom_stroke_expand(path, &st, &map, NULL, NULL, &so);
    if (r != GSHGEOM_OK) {
      fprintf(stderr, "stroke-curves-svg: %s: the polygon expansion failed: %s\n", c->name, gshgeom_result_string(r));
      return 1;
    }
    f = open_doc(argv[1], c->name, "poly");
    fprintf(f, "<path d=\"");
    const GSHGEOM_Flat_Path * fp = gshgeom_stroke_outline_body(so);
    const GSHGEOM_Device_Point * pts = gshgeom_flat_path_points(fp);
    for (size_t i = 0; i < gshgeom_flat_path_subpath_count(fp); i++) {
      GSHGEOM_Flat_Subpath sp;
      gshgeom_flat_path_subpath(fp, i, &sp);
      for (size_t k = 0; k < sp.count; k++) {
        fprintf(f, "%c%.10g %.10g", k ? 'L' : 'M', (double)pts[sp.first + k].x / 64.0, (double)pts[sp.first + k].y / 64.0);
      }
      fputc('Z', f);
    }
    fprintf(f, "\" fill=\"#000000\" fill-rule=\"nonzero\"/>\n</svg>\n");
    fclose(f);
    gshgeom_stroke_outline_destroy(so);
    gshgeom_path_destroy(path);
    printf("%s %d\n", c->name, c->curved);
  }
  return 0;
}
