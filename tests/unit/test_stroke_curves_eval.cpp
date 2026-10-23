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
 * The curve-winding evaluator's own check (tests/curve_eval.h), first, where the answer is exact: circles and
 * ellipses written as rational quadratics, and the cubic approximations of them (stroke_curves.h, SG40).  Its winding
 * agrees with the analytic answer at every point farther than `m` from the boundary, and a defect planted in the
 * evaluator (a flipped sign of a rational term) fails the same check.
 */

#include "../test_helpers.h"

#include "../curve_eval.h"

#include <random>

namespace ce = curve_eval;

namespace {

/** A circle (or, under (a, d), an ellipse with semi-axes r a and r d) as four rational quadratic quarter arcs. */
GSHGEOM_Path * circle_path(double cx, double cy, double r, double a, double d) {
  GSHGEOM_Path * p = nullptr;
  EXPECT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &p), GSHGEOM_OK);
  const double w[3] = {1, 0.70710678118654752440, 1};
  auto pt = [&](double x, double y) { return GSHGEOM_Point{cx + a * x * r, cy + d * y * r, 0.0}; };
  gshgeom_path_move_to(p, pt(1, 0));
  gshgeom_path_rational_quadratic_to(p, pt(1, 1), pt(0, 1), w);
  gshgeom_path_rational_quadratic_to(p, pt(-1, 1), pt(-1, 0), w);
  gshgeom_path_rational_quadratic_to(p, pt(-1, -1), pt(0, -1), w);
  gshgeom_path_rational_quadratic_to(p, pt(1, -1), pt(1, 0), w);
  gshgeom_path_close(p);
  return p;
}

/** Points that must wind as the ellipse says, farther than m from the boundary; returns the number that do not. */
int check(const GSHGEOM_Path * path, double cx, double cy, double rx, double ry, double m) {
  ce::Winding wd;
  wd.build(path, 1e-7L);
  std::mt19937_64 rng(7);
  std::uniform_real_distribution<double> u(-1.5, 1.5);
  int bad = 0;
  for (int i = 0; i < 4000; i++) {
    double px = cx + u(rng) * rx, py = cy + u(rng) * ry;
    /* the normalized radius: below 1 inside; the distance to the boundary is at least |rho - 1| min(rx, ry) */
    double rho = std::sqrt((px - cx) * (px - cx) / (rx * rx) + (py - cy) * (py - cy) / (ry * ry));
    if (std::fabs(rho - 1.0) * std::min(rx, ry) <= m) {
      continue;
    }
    int w = wd.at({px, py});
    if (w != (rho < 1.0 ? 1 : 0)) {
      bad++;
    }
  }
  return bad;
}

} // namespace

TEST(CurveEval, CirclesAndEllipsesWindAsTheyShould) {
  for (double a : {1.0, 2.0, 8.0}) {
    GSHGEOM_Path * p = circle_path(3, -2, 5, a, 1);
    EXPECT_EQ(check(p, 3, -2, 5 * a, 5, 1e-5), 0) << "scale " << a;
    gshgeom_path_destroy(p);
  }
}

TEST(CurveEval, TheCubicApproximationsWindAsTheyShould) {
  for (double a : {1.0, 3.0}) {
    GSHGEOM_Path * p = circle_path(0, 0, 4, a, 1);
    GSHGEOM_Path * low = nullptr;
    ASSERT_EQ(gshgeom_stroke_curves_lower_cubic(p, 1e-4, nullptr, nullptr, &low), GSHGEOM_OK);
    /* the cubics are within 1e-4 of the ellipse */
    EXPECT_EQ(check(low, 0, 0, 4 * a, 4, 2e-4), 0);
    gshgeom_path_destroy(low);
    gshgeom_path_destroy(p);
  }
}

TEST(CurveEval, AWindingNeverNegativeForAPositivelyWoundRing) {
  GSHGEOM_Path * p = circle_path(0, 0, 4, 1, 1);
  ce::Winding wd;
  wd.build(p, 1e-6L);
  for (double x = -5; x <= 5; x += 0.37) {
    for (double y = -5; y <= 5; y += 0.41) {
      EXPECT_GE(wd.least_ring({x, y}), 0);
    }
  }
  gshgeom_path_destroy(p);
}

TEST(CurveEval, ADefectPlantedInTheEvaluatorFailsTheCheck) {
  GSHGEOM_Path * p = circle_path(0, 0, 5, 1, 1);
  EXPECT_EQ(check(p, 0, 0, 5, 5, 1e-5), 0);
  ce::plant = 1;
  EXPECT_GT(check(p, 0, 0, 5, 5, 1e-5), 0);
  ce::plant = 0;
  gshgeom_path_destroy(p);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
