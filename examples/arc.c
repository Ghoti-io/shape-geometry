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
 * Exact arcs: a rounded rectangle built from corner arcs in degrees (each a quarter turn, so exact), an SVG `A` (the endpoint form), a circle built the way Canvas
 * `arc(x, y, r, 0, 2 pi)` is (the centre form), what a reader can ask of them (the verb, the record, the circle a
 * toolpath or a writer wants, the vertices a marker sits on), the exact pieces a rounded corner decomposes into, the
 * path placed by a map and flattened by the conic route, and a reader that does not know arcs being refused instead of
 * drawing the chord.
 *
 *     $ ./arc
 *     rounded rectangle: 8 elements, 4 arcs, bounds (0, 0)..(100, 60), hash ...
 *     ...
 */

#include <ghoti.io/shape-geometry/shape-geometry.h>

#include <inttypes.h>
#include <stdio.h>

static int check(GSHGEOM_Result r, const char * what) {
  if (r != GSHGEOM_OK) {
    fprintf(stderr, "arc: %s: %s\n", what, gshgeom_result_string(r));
    return 1;
  }
  return 0;
}

static GSHGEOM_Point pt(double x, double y) {
  GSHGEOM_Point p = {x, y, 0.0};
  return p;
}

/** One corner of radius @p r about (@p cx, @p cy), from the angle @p t0 degrees through a quarter turn: a quarter turn in degrees is
 *  exact, so the arc's control point is the corner of the box, not a rounding of it. */
static int corner(GSHGEOM_Path * path, double cx, double cy, double r, double t0) {
  GSHGEOM_Arc_Form f = {pt(cx, cy), pt(r, 0), pt(0, r), GSHGEOM_ANGLE_DEGREES, t0, 90.0, GSHGEOM_ARC_CIRCULAR};
  return check(gshgeom_path_arc_conjugate_to(path, &f, NULL), "corner");
}

typedef struct Marks {
  size_t vertices;
  size_t arc_ends;
  int exact_corner_directions;
} Marks;

static GSHGEOM_Result mark(const GSHGEOM_Vertex * v, void * ctx) {
  Marks * m = ctx;
  m->vertices++;
  if ((v->flags & GSHGEOM_VERTEX_HAS_IN) && (v->flags & GSHGEOM_VERTEX_HAS_OUT)) {
    /* at the end of a corner arc the incoming direction is an axis vector, and so is the line that follows */
    m->arc_ends += 1;
    m->exact_corner_directions += (v->in.x == 0.0 || v->in.y == 0.0) && (v->out.x == 0.0 || v->out.y == 0.0);
  }
  return GSHGEOM_OK;
}

int main(void) {
  int failed = 0;
  GSHGEOM_Path * path = NULL;
  if (check(gshgeom_path_create(NULL, NULL, 2, &path), "create")) {
    return 1;
  }

  /* A rounded rectangle 100 by 60 with corners of radius 10, clockwise on a y-down screen: a line, then an arc whose
   * sweep is positive (from +x toward +y).  The ends of each arc are exactly the numbers written here. */
  failed |= check(gshgeom_path_move_to(path, pt(10, 0)), "move");
  failed |= check(gshgeom_path_line_to(path, pt(90, 0)), "line");
  failed |= corner(path, 90, 10, 10, -90);
  failed |= check(gshgeom_path_line_to(path, pt(100, 50)), "line");
  failed |= corner(path, 90, 50, 10, 0);
  failed |= check(gshgeom_path_line_to(path, pt(10, 60)), "line");
  failed |= corner(path, 10, 50, 10, 90);
  failed |= check(gshgeom_path_line_to(path, pt(0, 10)), "line");
  failed |= corner(path, 10, 10, 10, 180);
  failed |= check(gshgeom_path_close(path), "close");

  GSHGEOM_Bounds box;
  failed |= check(gshgeom_path_bounds(path, &box), "bounds");
  size_t arcs = 0;
  for (size_t i = 0; i < gshgeom_path_element_count(path); i++) {
    GSHGEOM_Verb verb;
    failed |= check(gshgeom_path_verb(path, i, &verb), "verb");
    arcs += verb == GSHGEOM_VERB_ARC;
  }
  printf("rounded rectangle: %zu elements, %zu arcs, bounds (%g, %g)..(%g, %g), hash %016" PRIx64 "\n",
      gshgeom_path_element_count(path), arcs, box.min.x, box.min.y, box.max.x, box.max.y, gshgeom_path_hash(path));
  failed |= arcs != 4 || box.min.x != 0.0 || box.min.y != 0.0 || box.max.x != 100.0 || box.max.y != 60.0;

  /* A corner is a quarter turn, so it is one piece, a rational quadratic whose middle weight is exactly sqrt(1/2) and
   * whose control point is the corner of the box: (100, 0), exactly. */
  GSHGEOM_Arc a;
  GSHGEOM_Curve pieces[GSHGEOM_ARC_MAX_PIECES];
  int n = 0;
  failed |= check(gshgeom_path_arc(path, 2, &a), "read the arc back");
  failed |= check(gshgeom_arc_pieces(&a, pieces, &n), "pieces");
  printf("corner: %d piece, control point (%g, %g), weight %.17g, circular %d\n", n, pieces[0].p[1].x, pieces[0].p[1].y,
      pieces[0].w[1], (a.flags & GSHGEOM_ARC_CIRCULAR) != 0u);
  failed |= n != 1 || pieces[0].p[1].x != 100.0 || pieces[0].p[1].y != 0.0 || (a.flags & GSHGEOM_ARC_CIRCULAR) == 0u;

  /* A reader that does not know arcs is refused, never given the chord. */
  GSHGEOM_Element element;
  GSHGEOM_Result refused = gshgeom_path_element(path, 2, &element);
  printf("reading the arc as an ordinary element: %s\n", gshgeom_result_string(refused));
  failed |= refused != GSHGEOM_ERR_UNSUPPORTED;

  /* The vertices a marker sits on: one per element, the arc's end among them, with exact axis directions. */
  Marks marks = {0, 0, 0};
  failed |= check(gshgeom_path_vertices(path, mark, &marks), "vertices");
  printf("vertices: %zu (one per element); %zu are joins, %d of them with axis-parallel directions\n", marks.vertices,
      marks.arc_ends, marks.exact_corner_directions);
  failed |= marks.vertices != gshgeom_path_element_count(path) || marks.exact_corner_directions < 8;

  /* SVG's `A 10 10 0 0 1 20 0` from (0, 0): the endpoint form keeps both ends as given, whatever the rounding of its angles. */
  GSHGEOM_Path * svg = NULL;
  GSHGEOM_Arc_Kind kind = GSHGEOM_ARC_KIND_NOTHING;
  GSHGEOM_Arc sa;
  failed |= check(gshgeom_path_create(NULL, NULL, 2, &svg), "create");
  failed |= check(gshgeom_path_move_to(svg, pt(0, 0)), "move");
  failed |= check(gshgeom_path_arc_endpoint_to(svg, pt(20, 0), 10.0, 10.0, 0.0, GSHGEOM_ANGLE_DEGREES, 0, 1, &kind), "A");
  failed |= check(gshgeom_path_arc(svg, 1, &sa), "read the A back");
  printf("SVG A: kind %d (2 is an arc), ends (%g, %g)..(%g, %g), sweep %.17g radians\n", (int)kind, sa.start.x, sa.start.y, sa.end.x,
      sa.end.y, sa.sweep);
  failed |= kind != GSHGEOM_ARC_KIND_ARC || sa.end.x != 20.0 || sa.end.y != 0.0 || sa.sweep <= 0.0;
  gshgeom_path_destroy(svg);

  /* A circle the way Canvas builds one: the centre form with a sweep of exactly one turn closes to the bit. */
  GSHGEOM_Path * circle = NULL;
  failed |= check(gshgeom_path_create(NULL, NULL, 2, &circle), "create");
  failed |= check(gshgeom_path_arc_center_to(circle, pt(50, 30), 20.0, 20.0, 0.0, GSHGEOM_ANGLE_DEGREES, 0.0, 360.0, NULL),
      "circle");
  GSHGEOM_Arc c;
  GSHGEOM_Point centre;
  double radius = 0.0;
  int ccw = 0;
  failed |= check(gshgeom_path_arc(circle, 1, &c), "read the circle back");
  failed |= check(gshgeom_arc_circle(&c, &centre, &radius, &ccw), "circle");
  printf("circle: centre (%g, %g), radius %g, positive sense %d, closed to the bit %d, %s\n", centre.x, centre.y, radius, ccw,
      c.start.x == c.end.x && c.start.y == c.end.y, "a toolpath can say G2 or G3 and a writer <circle>");
  failed |= centre.x != 50.0 || centre.y != 30.0 || radius != 20.0 || !ccw || c.start.x != c.end.x || c.start.y != c.end.y;

  /* Under a map: a rotation keeps it a circle, a non-uniform scale does not (the flag follows the map's kind tag). */
  GSHGEOM_Map turn, stretch;
  GSHGEOM_Arc mapped;
  failed |= check(gshgeom_map_rotate(GSHGEOM_ANGLE_DEGREES, 30.0, &turn), "rotate");
  failed |= check(gshgeom_map_scale(2, 2.0, 1.0, 1.0, &stretch), "scale");
  failed |= check(gshgeom_arc_map(&turn, &c, &mapped), "map by a rotation");
  int kept = (mapped.flags & GSHGEOM_ARC_CIRCULAR) != 0u;
  failed |= check(gshgeom_arc_map(&stretch, &c, &mapped), "map by a stretch");
  int stretched = (mapped.flags & GSHGEOM_ARC_CIRCULAR) != 0u;
  printf("circular flag: kept under a rotation %d, kept under a stretch %d\n", kept, stretched);
  failed |= !kept || stretched;

  /* The device path under a rotation: the conic route flattens each piece to the tolerance, ends rounded once. */
  GSHGEOM_Device_Path * device = NULL;
  GSHGEOM_Flat_Path * flat = NULL;
  failed |= check(gshgeom_device_path_from_model(path, &turn, NULL, NULL, &device), "device path");
  failed |= check(gshgeom_device_path_flatten(device, NULL, NULL, &flat), "flatten");
  printf("rounded rectangle under a 30 degree turn: %zu device elements, flattened (revision %d) to %zu points, hash %016" PRIx64
         "\n",
      gshgeom_device_path_element_count(device), gshgeom_flatten_revision(), gshgeom_flat_path_point_count(flat),
      gshgeom_flat_path_hash(flat));
  failed |= gshgeom_flat_path_point_count(flat) < 8;

  gshgeom_flat_path_destroy(flat);
  gshgeom_device_path_destroy(device);
  gshgeom_path_destroy(circle);
  gshgeom_path_destroy(path);
  if (failed) {
    fprintf(stderr, "arc: FAILED\n");
  }
  return failed;
}
