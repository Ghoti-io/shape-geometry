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
 * A path from building to pixels' geometry: a model path made of a line, a
 * cubic and an exact quarter circle (a rational quadratic), placed with a map,
 * turned into the one 26.6 device path and flattened by the definition's integer
 * flattener; then measured, to find the point and the direction halfway along
 * the arc, as an animation or text on a path would.
 *
 * Model space is y down, like a screen (AD-6).  Nothing here touches a pixel:
 * a sink (shape-cpu, a GPU) takes the device path from here.
 *
 *     $ ./path
 *     model path: 4 elements, bounds (0, 0)..(100, 100), hash ...
 *     ...
 */

#include <ghoti.io/shape-geometry/shape-geometry.h>

#include <inttypes.h>
#include <stdio.h>

#define HALF_SQRT2 0.70710678118654752440 /* the weight that makes a rational quadratic a circular arc */
#define PI 3.14159265358979323846

static int check(GSHGEOM_Result r, const char * what) {
  if (r != GSHGEOM_OK) {
    fprintf(stderr, "path: %s: %s\n", what, gshgeom_result_string(r));
    return 1;
  }
  return 0;
}

static GSHGEOM_Point pt(double x, double y) {
  GSHGEOM_Point p = {x, y, 0.0};
  return p;
}

int main(void) {
  int failed = 0;
  GSHGEOM_Path * path = NULL;
  GSHGEOM_Device_Path * device = NULL;
  GSHGEOM_Flat_Path * flat = NULL;
  GSHGEOM_Measured * measured = NULL;

  /* A 2D path (NULL allocator and limits: the defaults).  Every call is checked; a NaN or an infinity is refused,
   * never stored. */
  if (check(gshgeom_path_create(NULL, NULL, 2, &path), "create")) {
    return 1;
  }
  const double w[3] = {1.0, HALF_SQRT2, 1.0};
  failed |= check(gshgeom_path_move_to(path, pt(0, 100)), "move");
  failed |= check(gshgeom_path_line_to(path, pt(0, 50)), "line");
  failed |= check(gshgeom_path_cubic_to(path, pt(0, 20), pt(20, 0), pt(50, 0)), "cubic");
  /* An exact quarter circle of radius 50 about (50, 50): no approximation by cubics. */
  failed |= check(gshgeom_path_rational_quadratic_to(path, pt(100, 0), pt(100, 50), w), "arc");
  failed |= check(gshgeom_path_close(path), "close");

  GSHGEOM_Bounds box;
  failed |= check(gshgeom_path_bounds(path, &box), "bounds");
  printf("model path: %zu elements, bounds (%g, %g)..(%g, %g), hash %016" PRIx64 "\n", gshgeom_path_element_count(path),
      box.min.x, box.min.y, box.max.x, box.max.y, gshgeom_path_hash(path));

  /* Place it: scale by a half, turn it by 30 degrees, move it by (20, 4).  Maps are composed only here
   * (outer after inner), and each carries a tag saying what kind it is. */
  GSHGEOM_Map scale, turn, move, tmp, place;
  failed |= check(gshgeom_map_scale(2, 0.5, 0.5, 1.0, &scale), "scale");
  failed |= check(gshgeom_map_rotate(GSHGEOM_ANGLE_DEGREES, 30.0, &turn), "rotate");
  failed |= check(gshgeom_map_translate(2, 20.0, 4.0, 0.0, &move), "translate");
  failed |= check(gshgeom_map_compose(&turn, &scale, &tmp), "compose");
  failed |= check(gshgeom_map_compose(&move, &tmp, &place), "compose");
  int aligned = 1;
  failed |= check(gshgeom_map_grid_aligned(&place, &aligned, NULL), "grid aligned");
  printf("placement: grid aligned %d (a rotated map: a glyph under it takes the general route)\n", aligned);

  /* The device path: every point rounded once, floor(x * 64 + 0.5), into 26.6 fixed point.  Curves stay curves. */
  failed |= check(gshgeom_device_path_from_model(path, &place, NULL, NULL, &device), "device path");
  printf("device path: %zu elements, hash %016" PRIx64 "\n", gshgeom_device_path_element_count(device),
      gshgeom_device_path_hash(device));

  /* Flattened by the integer flattener: within 1/16 pixel of the curves, and the same points on every target. */
  failed |= check(gshgeom_device_path_flatten(device, NULL, NULL, &flat), "flatten");
  printf("flattened (revision %d): %zu points in %zu subpath, hash %016" PRIx64 "\n", gshgeom_flatten_revision(),
      gshgeom_flat_path_point_count(flat), gshgeom_flat_path_subpath_count(flat), gshgeom_flat_path_hash(flat));
  const GSHGEOM_Device_Point * points = gshgeom_flat_path_points(flat);
  printf("  first point (%d, %d) in 1/64 px = (%.4f, %.4f) px\n", points[0].x, points[0].y, points[0].x / 64.0,
      points[0].y / 64.0);

  /* Distance along the path, in model units.  Element 3 (after the move, the line and the cubic) is the arc. */
  failed |= check(gshgeom_measured_create(NULL, NULL, path, 0.0, &measured), "measure");
  GSHGEOM_Element element;
  double arc = 0.0;
  failed |= check(gshgeom_path_element(path, 3, &element), "element");
  failed |= check(gshgeom_measured_curve_length(&element.curve, 0.0, &arc), "arc length");
  const double exact = PI * 50.0 / 2.0;
  printf("arc length %.12f, closed form %.12f\n", arc, exact);
  failed |= !(arc - exact < 1e-9 && exact - arc < 1e-9);

  /* Halfway along the arc: past the line (50) and the cubic.  The point is on the circle at 45 degrees, and the
   * tangent is the direction of travel. */
  GSHGEOM_Element cubic;
  double before = 0.0;
  failed |= check(gshgeom_path_element(path, 2, &cubic), "element");
  failed |= check(gshgeom_measured_curve_length(&cubic.curve, 0.0, &before), "cubic length");
  printf("whole path %.6f long (the closing line included)\n", gshgeom_measured_length(measured));
  GSHGEOM_Measured_Frame frame;
  failed |= check(gshgeom_measured_at(measured, 50.0 + before + arc / 2.0, &frame), "at");
  printf("halfway along the arc: point (%.6f, %.6f), tangent (%.6f, %.6f), element %zu\n", frame.point.x, frame.point.y,
      frame.tangent.x, frame.tangent.y, frame.element);
  const double on_circle = 50.0 + 50.0 * HALF_SQRT2;
  failed |= frame.element != 3 || !frame.has_tangent;
  failed |= !(frame.point.x - on_circle < 1e-9 && on_circle - frame.point.x < 1e-9);

  gshgeom_measured_destroy(measured);
  gshgeom_flat_path_destroy(flat);
  gshgeom_device_path_destroy(device);
  gshgeom_path_destroy(path);
  if (failed) {
    fprintf(stderr, "path: FAILED\n");
  }
  return failed;
}
