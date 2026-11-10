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
 * A quadratic or cubic that is a straight line in all but name strokes as that line (stroke.h, SG10): its control points lie
 * within `m * sqrt(1e-5)` of the chord, `m` the larger coordinate difference of its ends, and run along it from the first end
 * to the last.  Without this the direction at an end is the direction of a control point a hair from it, and one that is a
 * millionth off swings the cap or join through any angle (resvg's line-as-curve and control-points-clamping cases).
 */

#include "../test_helpers.h"

#include <cstdint>

namespace {

GSHGEOM_Point P(double x, double y) {
  return GSHGEOM_Point{x, y, 0.0};
}

GSHGEOM_Stroke style(int join, int cap) {
  GSHGEOM_Stroke st;
  gshgeom_stroke_init(&st);
  st.width = 20.0;
  st.join = (GSHGEOM_Stroke_Join)join;
  st.cap = (GSHGEOM_Stroke_Cap)cap;
  return st;
}

uint64_t hash_of(GSHGEOM_Path * p, const GSHGEOM_Stroke & st) {
  GSHGEOM_Map m;
  EXPECT_EQ(gshgeom_map_identity(2, &m), GSHGEOM_OK);
  GSHGEOM_Stroke_Outline * o = nullptr;
  EXPECT_EQ(gshgeom_stroke_expand(p, &st, &m, nullptr, nullptr, &o), GSHGEOM_OK);
  const uint64_t h = o ? gshgeom_stroke_outline_hash(o) : 0;
  gshgeom_stroke_outline_destroy(o);
  return h;
}

struct Pt {
  double x, y;
};

/** M a, then one cubic (c1, c2, end), then an optional line to @p tail. */
GSHGEOM_Path * cubic(Pt a, Pt c1, Pt c2, Pt e, const Pt * tail = nullptr) {
  GSHGEOM_Path * p = nullptr;
  EXPECT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &p), GSHGEOM_OK);
  gshgeom_path_move_to(p, P(a.x, a.y));
  gshgeom_path_cubic_to(p, P(c1.x, c1.y), P(c2.x, c2.y), P(e.x, e.y));
  if (tail) {
    gshgeom_path_line_to(p, P(tail->x, tail->y));
  }
  return p;
}

GSHGEOM_Path * quad(Pt a, Pt c, Pt e) {
  GSHGEOM_Path * p = nullptr;
  EXPECT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &p), GSHGEOM_OK);
  gshgeom_path_move_to(p, P(a.x, a.y));
  gshgeom_path_quadratic_to(p, P(c.x, c.y), P(e.x, e.y));
  return p;
}

GSHGEOM_Path * line(Pt a, Pt e, const Pt * tail = nullptr) {
  GSHGEOM_Path * p = nullptr;
  EXPECT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &p), GSHGEOM_OK);
  gshgeom_path_move_to(p, P(a.x, a.y));
  gshgeom_path_line_to(p, P(e.x, e.y));
  if (tail) {
    gshgeom_path_line_to(p, P(tail->x, tail->y));
  }
  return p;
}

}  // namespace

TEST(StrokeNearLine, ACubicWithControlPointsAHairFromItsEndsStrokesAsTheLineBetweenThem) {
  // line-as-curve-1 and -2 of the resvg suite: control points 0.0001 and 0.1 off the end points of a 198-unit diagonal.
  for (double e : {0.0001, 0.1}) {
    GSHGEOM_Path * c = cubic({20, 160}, {20 + e, 160}, {160, 20}, {160 + e, 20});
    GSHGEOM_Path * l = line({20, 160}, {160 + e, 20});
    for (int cap : {GSHGEOM_STROKE_CAP_BUTT, GSHGEOM_STROKE_CAP_ROUND, GSHGEOM_STROKE_CAP_SQUARE}) {
      const GSHGEOM_Stroke st = style(GSHGEOM_STROKE_JOIN_MITER, cap);
      EXPECT_EQ(hash_of(c, st), hash_of(l, st)) << "offset " << e << ", cap " << cap;
    }
    gshgeom_path_destroy(c);
    gshgeom_path_destroy(l);
  }
}

TEST(StrokeNearLine, TheJoinAfterAControlPointThatOvershootsItsEndIsTheLinesJoin) {
  // control-points-clamping-1: the second control point is 0.1 beside the end, so the end direction was the horizontal blip
  // and the bevel was never made.
  const Pt tail{15, 175};
  GSHGEOM_Path * c = cubic({175, 15}, {175, 30}, {175.1, 175}, {175, 175}, &tail);
  GSHGEOM_Path * l = line({175, 15}, {175, 175}, &tail);
  const GSHGEOM_Stroke st = style(GSHGEOM_STROKE_JOIN_BEVEL, GSHGEOM_STROKE_CAP_BUTT);
  EXPECT_EQ(hash_of(c, st), hash_of(l, st));
  gshgeom_path_destroy(c);
  gshgeom_path_destroy(l);
}

TEST(StrokeNearLine, AQuadraticIsReducedTheSameWay) {
  GSHGEOM_Path * q = quad({0, 0}, {50, 0.01}, {100, 0});
  GSHGEOM_Path * l = line({0, 0}, {100, 0});
  const GSHGEOM_Stroke st = style(GSHGEOM_STROKE_JOIN_ROUND, GSHGEOM_STROKE_CAP_SQUARE);
  EXPECT_EQ(hash_of(q, st), hash_of(l, st));
  gshgeom_path_destroy(q);
  gshgeom_path_destroy(l);
}

TEST(StrokeNearLine, ACurveOffItsChordByMoreThanTheSlopOrThatDoublesBackIsNotALine) {
  const GSHGEOM_Stroke st = style(GSHGEOM_STROKE_JOIN_ROUND, GSHGEOM_STROKE_CAP_BUTT);
  GSHGEOM_Path * l = line({0, 0}, {100, 0});
  const uint64_t line_hash = hash_of(l, st);
  // a control point 1 unit off a 100-unit chord: 1 > 100 * sqrt(1e-5) = 0.316
  GSHGEOM_Path * off = quad({0, 0}, {50, 1.0}, {100, 0});
  EXPECT_NE(hash_of(off, st), line_hash);
  // the first control point behind the start: the curve leaves the way it came
  GSHGEOM_Path * back = cubic({0, 0}, {-0.001, 0}, {60, 0}, {100, 0});
  EXPECT_NE(hash_of(back, st), line_hash);
  // a control point past the far end: it overshoots and doubles back to it
  GSHGEOM_Path * over = cubic({0, 0}, {40, 0}, {100.001, 0}, {100, 0});
  EXPECT_NE(hash_of(over, st), line_hash);
  // control points in the wrong order along the chord: it goes out, comes back and goes out again
  GSHGEOM_Path * order = cubic({0, 0}, {70, 0}, {30, 0}, {100, 0});
  EXPECT_NE(hash_of(order, st), line_hash);
  gshgeom_path_destroy(order);
  gshgeom_path_destroy(l);
  gshgeom_path_destroy(off);
  gshgeom_path_destroy(back);
  gshgeom_path_destroy(over);
}

TEST(StrokeNearLine, TheDecisionIsMadeBeforeTheMapSoAScaleDoesNotChangeIt) {
  // the same curve under a scale of 100 in x only: the reduction is on the points as the path holds them, so the outline is
  // the line's under the same map
  GSHGEOM_Path * c = cubic({20, 160}, {20.0001, 160}, {160, 20}, {160.0001, 20});
  GSHGEOM_Path * l = line({20, 160}, {160.0001, 20});
  GSHGEOM_Map m;
  ASSERT_EQ(gshgeom_map_scale(2, 100.0, 1.0, 1.0, &m), GSHGEOM_OK);
  const GSHGEOM_Stroke st = style(GSHGEOM_STROKE_JOIN_MITER, GSHGEOM_STROKE_CAP_BUTT);
  GSHGEOM_Stroke_Outline * a = nullptr;
  GSHGEOM_Stroke_Outline * b = nullptr;
  ASSERT_EQ(gshgeom_stroke_expand(c, &st, &m, nullptr, nullptr, &a), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_stroke_expand(l, &st, &m, nullptr, nullptr, &b), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_stroke_outline_hash(a), gshgeom_stroke_outline_hash(b));
  gshgeom_stroke_outline_destroy(a);
  gshgeom_stroke_outline_destroy(b);
  gshgeom_path_destroy(c);
  gshgeom_path_destroy(l);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
