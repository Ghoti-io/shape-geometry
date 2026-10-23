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
 * The lowering of a path to lines and cubics (stroke_curves.h, SG40) against the exact model of
 * tools/gen-curves-vectors.py: the handle fraction, the weight of a bisected half, the radial error (the closed form
 * against the model's own scan of the cubic), the number of halvings, the bound over circles and ellipses, a
 * quadratic raised exactly, a rational cubic, and the limits.
 */

#include "../test_helpers.h"

#include "../curve_eval.h"
#include "../curves_vectors.h"

#include <vector>

namespace ce = curve_eval;
namespace cv = curves_vectors;

namespace {

GSHGEOM_Point pt(double x, double y) {
  return GSHGEOM_Point{x, y, 0.0};
}

/** The arc of radius r with middle weight w (ends weight 1), centred on the origin and symmetric about +x. */
GSHGEOM_Path * arc_path(double w, double r, double sx = 1.0, double sy = 1.0, const double * weights = nullptr) {
  GSHGEOM_Path * p = nullptr;
  EXPECT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &p), GSHGEOM_OK);
  double s = std::sqrt(1.0 - w * w);
  double ww[3] = {1, w, 1};
  if (weights) {
    ww[0] = weights[0];
    ww[1] = weights[1];
    ww[2] = weights[2];
  }
  EXPECT_EQ(gshgeom_path_move_to(p, pt(sx * r * w, -sy * r * s)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_rational_quadratic_to(p, pt(sx * r / w, 0), pt(sx * r * w, sy * r * s), ww), GSHGEOM_OK);
  return p;
}

std::vector<GSHGEOM_Element> elements(const GSHGEOM_Path * p) {
  std::vector<GSHGEOM_Element> v;
  for (size_t i = 0; i < gshgeom_path_element_count(p); i++) {
    GSHGEOM_Element e;
    gshgeom_path_element(p, i, &e);
    v.push_back(e);
  }
  return v;
}

size_t cubic_count(const GSHGEOM_Path * p) {
  size_t n = 0;
  for (auto & e : elements(p)) {
    n += e.verb == GSHGEOM_VERB_CUBIC;
  }
  return n;
}

/** The farthest any point of `a` is from the polyline of `b`, both sampled densely. */
ce::Real one_sided(const std::vector<ce::P2> & a, const std::vector<ce::P2> & b) {
  ce::Real worst = 0;
  for (auto & p : a) {
    ce::Real best = 1e300L;
    for (size_t i = 0; i + 1 < b.size(); i++) {
      best = std::min(best, ce::dist_to_segment(p, b[i], b[i + 1]));
    }
    worst = std::max(worst, best);
  }
  return worst;
}

std::vector<ce::P2> dense(const GSHGEOM_Path * p, int per) {
  return ce::sample_path(p, per);
}

} // namespace

TEST(CurvesLower, TheModelsHandleWeightAndError) {
  for (size_t i = 0; i < cv::arc_count; i++) {
    const cv::Arc & a = cv::arcs[i];
    SCOPED_TRACE(a.angle_deg);
    GSHGEOM_Path * p = arc_path(a.w, 1.0);
    GSHGEOM_Path * low = nullptr;
    /* a tolerance no arc needs a halving for: one cubic, the model's */
    ASSERT_EQ(gshgeom_stroke_curves_lower_cubic(p, 1.0, nullptr, nullptr, &low), GSHGEOM_OK);
    auto el = elements(low);
    ASSERT_EQ(el.size(), 2u);
    ASSERT_EQ(el[1].verb, GSHGEOM_VERB_CUBIC);
    EXPECT_NEAR(el[1].curve.p[1].x, a.q1x, 4e-15 * std::max(1.0, std::fabs(a.q1x)));
    EXPECT_NEAR(el[1].curve.p[1].y, a.q1y, 4e-15 * std::max(1.0, std::fabs(a.q1y)));
    /* the handle fraction of the leg P0 -> P1 */
    double lam = (el[1].curve.p[1].x - el[0].curve.p[0].x) / (1.0 / a.w - el[0].curve.p[0].x);
    EXPECT_NEAR(lam, a.lambda, 5e-13);
    /* the closed form is the model's scan: e = sqrt(1 + (1-w)^3 / (27 (1+w))) - 1 */
    double x = (1 - a.w) * (1 - a.w) * (1 - a.w) / (27 * (1 + a.w));
    double e = x / (std::sqrt(1 + x) + 1);
    EXPECT_NEAR(e, a.error, 1e-9 * a.error + 1e-17);
    /* and the cubic's measured radial error never passes it */
    ce::Real worst = 0;
    ce::Elem c = ce::elem_of(el[1]);
    c.p[0] = {el[0].curve.p[0].x, el[0].curve.p[0].y};
    for (int k = 0; k <= 2000; k++) {
      ce::P2 q = ce::eval(c, (ce::Real)k / 2000);
      worst = std::max(worst, std::sqrt(q.x * q.x + q.y * q.y) - 1);
    }
    EXPECT_LE((double)worst, e * (1 + 1e-9) + 1e-15);
    EXPECT_GE((double)worst, e * 0.999 - 1e-15);
    /* the bisected half's weight is the model's: the first half runs from P0 to the arc's midpoint (1, 0) */
    GSHGEOM_Path * halved = nullptr;
    ASSERT_EQ(gshgeom_stroke_curves_lower_cubic(p, a.error * 0.9, nullptr, nullptr, &halved), GSHGEOM_OK);
    if (a.error > 1e-14) {
      auto h = elements(halved);
      ASSERT_GE(h.size(), 3u);
      EXPECT_NEAR(h[1].curve.p[3].x, 1.0, 1e-14); /* the exact bisection of the conic lands on the circle */
      EXPECT_NEAR(h[1].curve.p[3].y, 0.0, 1e-14);
      /* one half is the arc of weight w_half: its handle fraction is the model's formula at w_half */
      double wh = a.w_half;
      double lh = 4 * wh / (3 * (1 + wh));
      GSHGEOM_Point p0 = h[0].curve.p[0], q1 = h[1].curve.p[1];
      double s = std::sqrt(1 - a.w * a.w);
      /* the first half's control point P1 is (P0 + w P1) / (1 + w) */
      double p1x = (p0.x + a.w * (1.0 / a.w)) / (1 + a.w), p1y = (p0.y + a.w * 0.0) / (1 + a.w);
      (void)s;
      EXPECT_NEAR(q1.x, p0.x + lh * (p1x - p0.x), 1e-14);
      EXPECT_NEAR(q1.y, p0.y + lh * (p1y - p0.y), 1e-14);
    }
    gshgeom_path_destroy(halved);
    gshgeom_path_destroy(low);
    gshgeom_path_destroy(p);
  }
}

TEST(CurvesLower, TheNumberOfHalvingsIsTheModels) {
  for (size_t i = 0; i < cv::halvings_count; i++) {
    const cv::Halvings & h = cv::halvings_table[i];
    GSHGEOM_Path * p = arc_path(h.w, h.radius);
    GSHGEOM_Path * low = nullptr;
    ASSERT_EQ(gshgeom_stroke_curves_lower_cubic(p, h.tolerance, nullptr, nullptr, &low), GSHGEOM_OK);
    EXPECT_EQ(cubic_count(low), (size_t)1 << h.k) << "w " << h.w << " r " << h.radius << " tol " << h.tolerance;
    /* and it is within the tolerance: the cubics against the exact arc */
    /* and it is within the tolerance of the exact circle: every sample of the cubics is within it radially */
    auto got = dense(low, 400);
    ce::Real worst = 0;
    for (auto & q : got) {
      worst = std::max<ce::Real>(worst, std::fabs(std::sqrt(q.x * q.x + q.y * q.y) - h.radius));
    }
    EXPECT_LE((double)worst, h.tolerance * 1.001) << "w " << h.w << " r " << h.radius;
    gshgeom_path_destroy(low);
    gshgeom_path_destroy(p);
  }
}

TEST(CurvesLower, AnEllipseArcIsWithinTheBoundTimesItsLargestSemiAxis) {
  for (double scale : {1.0, 3.0, 8.0}) {
    for (double w : {0.70710678118654752440, 0.5, 0.2}) {
      double tol = 1e-3;
      GSHGEOM_Path * p = arc_path(w, 10.0, scale, 1.0);
      GSHGEOM_Path * low = nullptr;
      ASSERT_EQ(gshgeom_stroke_curves_lower_cubic(p, tol, nullptr, nullptr, &low), GSHGEOM_OK);
      auto got = dense(low, 600);
      std::vector<ce::P2> exact;
      double s = std::sqrt(1 - w * w);
      double a0 = std::atan2(-s, w), a1 = std::atan2(s, w);
      for (int k = 0; k <= 3000; k++) {
        double th = a0 + (a1 - a0) * k / 3000.0;
        exact.push_back({scale * 10.0 * std::cos(th), 10.0 * std::sin(th)});
      }
      EXPECT_LE((double)one_sided(got, exact), tol * 1.001) << "scale " << scale << " w " << w;
      EXPECT_LE((double)one_sided(exact, got), tol * 1.001) << "scale " << scale << " w " << w;
      gshgeom_path_destroy(low);
      gshgeom_path_destroy(p);
    }
  }
}

TEST(CurvesLower, WeightsThatAreNotInStandardFormGiveTheSameCurve) {
  double w = 0.6;
  double scaled[3] = {2, 2 * w, 2}; /* weights (a, b, c) = (2, 1.2, 2): w = b / sqrt(a c) = 0.6 */
  double repar[3] = {1, 2 * w, 4};  /* reparametrized: b / sqrt(a c) = 1.2 / 2 = 0.6 */
  GSHGEOM_Path * p0 = arc_path(w, 5.0);
  GSHGEOM_Path * p1 = arc_path(w, 5.0, 1, 1, scaled);
  GSHGEOM_Path * p2 = arc_path(w, 5.0, 1, 1, repar);
  GSHGEOM_Path *l0 = nullptr, *l1 = nullptr, *l2 = nullptr;
  ASSERT_EQ(gshgeom_stroke_curves_lower_cubic(p0, 1e-3, nullptr, nullptr, &l0), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_stroke_curves_lower_cubic(p1, 1e-3, nullptr, nullptr, &l1), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_stroke_curves_lower_cubic(p2, 1e-3, nullptr, nullptr, &l2), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_hash(l0), gshgeom_path_hash(l1));
  EXPECT_EQ(gshgeom_path_hash(l0), gshgeom_path_hash(l2));
  for (GSHGEOM_Path * q : {p0, p1, p2, l0, l1, l2}) {
    gshgeom_path_destroy(q);
  }
}

TEST(CurvesLower, AQuadraticIsRaisedExactlyAndOtherElementsAreCopied) {
  GSHGEOM_Path * p = nullptr;
  ASSERT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &p), GSHGEOM_OK);
  gshgeom_path_move_to(p, pt(0, 0));
  gshgeom_path_line_to(p, pt(3, 0));
  gshgeom_path_quadratic_to(p, pt(6, 9), pt(9, 0));
  gshgeom_path_cubic_to(p, pt(10, 1), pt(11, 2), pt(12, 3));
  double w1[3] = {1, 1, 1};
  gshgeom_path_rational_quadratic_to(p, pt(14, 5), pt(16, 3), w1); /* an ordinary quadratic by its weights */
  gshgeom_path_close(p);
  GSHGEOM_Path * low = nullptr;
  ASSERT_EQ(gshgeom_stroke_curves_lower_cubic(p, 1e-6, nullptr, nullptr, &low), GSHGEOM_OK);
  auto el = elements(low);
  ASSERT_EQ(el.size(), 6u);
  EXPECT_EQ(el[1].verb, GSHGEOM_VERB_LINE);
  EXPECT_EQ(el[2].verb, GSHGEOM_VERB_CUBIC);
  EXPECT_EQ(el[3].verb, GSHGEOM_VERB_CUBIC);
  EXPECT_EQ(el[4].verb, GSHGEOM_VERB_CUBIC);
  EXPECT_EQ(el[5].verb, GSHGEOM_VERB_CLOSE);
  /* the raised quadratic is the same curve at every parameter */
  GSHGEOM_Element q;
  gshgeom_path_element(p, 2, &q);
  ce::Elem a = ce::elem_of(q), b = ce::elem_of(el[2]);
  for (int k = 0; k <= 20; k++) {
    ce::P2 u = ce::eval(a, k / 20.0L), v = ce::eval(b, k / 20.0L);
    EXPECT_NEAR((double)u.x, (double)v.x, 1e-14);
    EXPECT_NEAR((double)u.y, (double)v.y, 1e-14);
  }
  /* and so is the one given by weights of 1: the same curve at every parameter, not a nearby one */
  GSHGEOM_Element rq;
  gshgeom_path_element(p, 4, &rq);
  ce::Elem ra = ce::elem_of(rq), rb = ce::elem_of(el[4]);
  for (int k = 0; k <= 20; k++) {
    ce::P2 u = ce::eval(ra, k / 20.0L), v = ce::eval(rb, k / 20.0L);
    EXPECT_NEAR((double)u.x, (double)v.x, 1e-14) << k;
    EXPECT_NEAR((double)u.y, (double)v.y, 1e-14) << k;
  }
  /* the cubic and the line are copied unchanged */
  EXPECT_EQ(el[3].curve.p[1].x, 10.0);
  EXPECT_EQ(el[3].curve.p[3].y, 3.0);
  gshgeom_path_destroy(low);
  gshgeom_path_destroy(p);
}

TEST(CurvesLower, ARationalCubicIsWithinTheToleranceAndTheWeightsMakeItConverge) {
  GSHGEOM_Path * p = nullptr;
  ASSERT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &p), GSHGEOM_OK);
  gshgeom_path_move_to(p, pt(0, 0));
  double w[4] = {1, 2, 0.5, 1.5};
  gshgeom_path_rational_cubic_to(p, pt(3, 6), pt(7, 6), pt(10, 0), w);
  for (double tol : {1e-2, 1e-4}) {
    GSHGEOM_Path * low = nullptr;
    ASSERT_EQ(gshgeom_stroke_curves_lower_cubic(p, tol, nullptr, nullptr, &low), GSHGEOM_OK);
    auto got = dense(low, 40);
    auto exact = dense(p, 1500);
    EXPECT_LE((double)one_sided(got, exact), tol) << tol;
    EXPECT_LE((double)one_sided(exact, got), tol) << tol;
    gshgeom_path_destroy(low);
  }
  /* constant weights: the polynomial curve itself, in one piece */
  GSHGEOM_Path * q = nullptr;
  ASSERT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &q), GSHGEOM_OK);
  gshgeom_path_move_to(q, pt(0, 0));
  double wc[4] = {2, 2, 2, 2};
  gshgeom_path_rational_cubic_to(q, pt(3, 6), pt(7, 6), pt(10, 0), wc);
  GSHGEOM_Path * low = nullptr;
  ASSERT_EQ(gshgeom_stroke_curves_lower_cubic(q, 1e-9, nullptr, nullptr, &low), GSHGEOM_OK);
  EXPECT_EQ(cubic_count(low), 1u);
  gshgeom_path_destroy(low);
  gshgeom_path_destroy(q);
  gshgeom_path_destroy(p);
}

TEST(CurvesLower, ExtremeAndHyperbolicConicsLower) {
  for (double w : {1e-6, 1e-3, 1.0000001, 3.0, 40.0}) {
    GSHGEOM_Path * p = nullptr;
    ASSERT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &p), GSHGEOM_OK);
    gshgeom_path_move_to(p, pt(0, 0));
    double ww[3] = {1, w, 1};
    gshgeom_path_rational_quadratic_to(p, pt(5, 8), pt(10, 0), ww);
    GSHGEOM_Path * low = nullptr;
    ASSERT_EQ(gshgeom_stroke_curves_lower_cubic(p, 1e-3, nullptr, nullptr, &low), GSHGEOM_OK) << w;
    auto got = dense(low, 100);
    auto exact = dense(p, 4000);
    EXPECT_LE((double)one_sided(got, exact), 1.001e-3) << w;
    EXPECT_LE((double)one_sided(exact, got), 1.001e-3) << w;
    gshgeom_path_destroy(low);
    gshgeom_path_destroy(p);
  }
}

TEST(CurvesLower, ErrorsAndLimits) {
  GSHGEOM_Path * p = arc_path(0.7, 1.0);
  GSHGEOM_Path * low = nullptr;
  EXPECT_EQ(gshgeom_stroke_curves_lower_cubic(nullptr, 1e-3, nullptr, nullptr, &low), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_stroke_curves_lower_cubic(p, 1e-3, nullptr, nullptr, nullptr), GSHGEOM_ERR_INVALID);
  for (double t : {0.0, -1.0, std::nan(""), HUGE_VAL}) {
    EXPECT_EQ(gshgeom_stroke_curves_lower_cubic(p, t, nullptr, nullptr, &low), GSHGEOM_ERR_INVALID);
  }
  EXPECT_EQ(low, nullptr);
  GSHGEOM_Path * p3 = nullptr;
  ASSERT_EQ(gshgeom_path_create(nullptr, nullptr, 3, &p3), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_stroke_curves_lower_cubic(p3, 1e-3, nullptr, nullptr, &low), GSHGEOM_ERR_INVALID);
  /* a result past max_elements */
  GSHGEOM_Limits lim;
  gshgeom_limits_default(&lim);
  lim.max_elements = 3;
  EXPECT_EQ(gshgeom_stroke_curves_lower_cubic(p, 1e-9, nullptr, &lim, &low), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(low, nullptr);
  /* a depth no tolerance can reach: a tolerance far below the arithmetic */
  GSHGEOM_Path * big = arc_path(0.7, 1e12);
  EXPECT_EQ(gshgeom_stroke_curves_lower_cubic(big, 1e-30, nullptr, nullptr, &low), GSHGEOM_ERR_LIMIT);
  gshgeom_path_destroy(big);
  gshgeom_path_destroy(p3);
  gshgeom_path_destroy(p);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
