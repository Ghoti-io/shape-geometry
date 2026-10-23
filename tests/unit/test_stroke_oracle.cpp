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
 * The geometric oracle of a stroke (tests/stroke_oracle.h), checked first on the
 * case it can be checked against exactly, and against the exact verdicts of
 * tools/gen-stroke-vectors.py.
 *
 * With round joins and round caps a stroke covers exactly the points within half
 * the width of its centerline, so the oracle's union of bodies, cusp discs,
 * sectors and half discs must agree with a direct distance computation wherever
 * the distance is not within `m` of `h`: on lines, every curve kind, closed and
 * open subpaths, degenerate ones, under a non-uniform map in both width spaces.
 * The oracle with a sector put on the inner side of a join (its own planted
 * defect) must disagree somewhere, every run.
 */

#include "../test_helpers.h"

#include "../stroke_matrix.h"
#include "../stroke_oracle.h"
#include "../stroke_vectors.h"

#include <cstdio>
#include <random>
#include <string>
#include <vector>

namespace so = stroke_oracle;

namespace {

GSHGEOM_Point pt(double x, double y) {
  GSHGEOM_Point p = {x, y, 0.0};
  return p;
}

/** A path built from a compact description: M x y, L x y, Q cx cy x y, C ..., R (rational quadratic) cx cy x y w0 w1
 * w2, K (rational cubic) c1 c2 x y w0..w3, Z. */
GSHGEOM_Path * build(const std::vector<std::pair<char, std::vector<double>>> & ops) {
  GSHGEOM_Path * p = nullptr;
  EXPECT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &p), GSHGEOM_OK);
  for (const auto & op : ops) {
    const std::vector<double> & a = op.second;
    GSHGEOM_Result r = GSHGEOM_OK;
    switch (op.first) {
    case 'M':
      r = gshgeom_path_move_to(p, pt(a[0], a[1]));
      break;
    case 'L':
      r = gshgeom_path_line_to(p, pt(a[0], a[1]));
      break;
    case 'Q':
      r = gshgeom_path_quadratic_to(p, pt(a[0], a[1]), pt(a[2], a[3]));
      break;
    case 'C':
      r = gshgeom_path_cubic_to(p, pt(a[0], a[1]), pt(a[2], a[3]), pt(a[4], a[5]));
      break;
    case 'R': {
      double w[3] = {a[4], a[5], a[6]};
      r = gshgeom_path_rational_quadratic_to(p, pt(a[0], a[1]), pt(a[2], a[3]), w);
      break;
    }
    case 'K': {
      double w[4] = {a[6], a[7], a[8], a[9]};
      r = gshgeom_path_rational_cubic_to(p, pt(a[0], a[1]), pt(a[2], a[3]), pt(a[4], a[5]), w);
      break;
    }
    case 'Z':
      r = gshgeom_path_close(p);
      break;
    }
    EXPECT_EQ(r, GSHGEOM_OK) << op.first;
  }
  return p;
}

typedef std::vector<std::pair<char, std::vector<double>>> Ops;

struct SelfCase {
  const char * name;
  Ops ops;
  double width;
};

/** The circle of radius r about (cx, cy) as four rational quadratics (weight sqrt(2)/2 in the middle). */
Ops circle(double cx, double cy, double r) {
  double w = 0.70710678118654752440;
  return Ops{{'M', {cx + r, cy}}, {'R', {cx + r, cy + r, cx, cy + r, 1, w, 1}}, {'R', {cx - r, cy + r, cx - r, cy, 1, w, 1}},
      {'R', {cx - r, cy - r, cx, cy - r, 1, w, 1}}, {'R', {cx + r, cy - r, cx + r, cy, 1, w, 1}}, {'Z', {}}};
}

std::vector<SelfCase> self_cases() {
  std::vector<SelfCase> c;
  c.push_back({"line", {{'M', {0, 0}}, {'L', {10, 3}}}, 4});
  c.push_back({"polyline_acute", {{'M', {0, 0}}, {'L', {10, 0}}, {'L', {1, 2}}, {'L', {6, 9}}}, 3});
  c.push_back({"polyline_short_segment", {{'M', {0, 0}}, {'L', {10, 0}}, {'L', {10.25, 0.5}}, {'L', {0, 6}}}, 5});
  c.push_back({"reversal", {{'M', {0, 0}}, {'L', {10, 0}}, {'L', {3, 0}}}, 4});
  c.push_back({"closed_square", {{'M', {0, 0}}, {'L', {8, 0}}, {'L', {8, 8}}, {'L', {0, 8}}, {'Z', {}}}, 3});
  c.push_back({"quadratic", {{'M', {0, 0}}, {'Q', {5, 9, 10, 0}}}, 3});
  c.push_back({"cubic_inflection", {{'M', {0, 0}}, {'C', {4, 8, 6, -8, 10, 0}}}, 2.5});
  c.push_back({"cubic_loop", {{'M', {0, 0}}, {'C', {12, 8, -2, 8, 10, 0}}}, 2});
  c.push_back({"cubic_cusp", {{'M', {0, 0}}, {'C', {8, 4, 4, 4, 4, 0}}}, 3});
  c.push_back({"cubic_near_cusp", {{'M', {0, 0}}, {'C', {8, 4, 4.01, 4, 4, 0}}}, 3});
  c.push_back({"circle_h_below_r", circle(0, 0, 5), 4});
  c.push_back({"circle_h_is_r", circle(0, 0, 3), 6});
  c.push_back({"circle_h_above_r", circle(0, 0, 2), 7});
  c.push_back({"rational_cubic", {{'M', {0, 0}}, {'K', {3, 6, 7, 6, 10, 0, 1, 2, 0.5, 1.5}}}, 3});
  c.push_back({"rational_extreme_weights", {{'M', {0, 0}}, {'R', {5, 8, 10, 0, 1, 40, 1}}}, 2});
  c.push_back({"degenerate_and_lone_move", {{'M', {0, 0}}, {'L', {6, 0}}, {'M', {3, 5}}, {'L', {3, 5}}, {'M', {9, 9}}}, 3});
  c.push_back({"overlapping_subpaths", {{'M', {0, 0}}, {'L', {10, 10}}, {'M', {0, 10}}, {'Q', {5, -6, 10, 0}}}, 3});
  return c;
}

/** The distance from a model point to the drawn centerline: every segment with extent and every degenerate point. */
so::Real centerline_distance(const so::Centerline & cl, so::V2 x) {
  so::Real best = 1e300L;
  for (const so::Sub & s : cl.subs) {
    if (s.segs.empty()) {
      if (s.draws) {
        best = std::min(best, so::norm(x - s.start));
      }
      continue;
    }
    for (const so::Seg & g : s.segs) {
      best = std::min(best, so::distance_to(so::make_poly(g), x));
    }
  }
  return best;
}

struct Agreement {
  long decided = 0;
  long disagree = 0;
  std::string first;
};

/**
 * Compare the oracle's verdict with dist <= h over a grid and random points about the centerline's box, in device
 * space, for both width spaces under the given map.
 */
Agreement self_check(const SelfCase & sc, const so::Affine & map, int space, int plant) {
  GSHGEOM_Path * path = build(sc.ops);
  so::Centerline cl = so::read_centerline(path);
  so::Style st;
  st.width = sc.width;
  st.space = space;
  st.join = so::JOIN_ROUND;
  st.cap = so::CAP_ROUND;
  st.plant = plant;
  so::Stroked sk = so::build_stroke(cl, st, map, 0);
  so::Centerline mapped = cl;
  if (space == 1) {
    for (so::Sub & s : mapped.subs) {
      s.start = map.apply(s.start);
      for (so::Seg & g : s.segs) {
        g = so::map_seg(g, map);
      }
    }
  }
  GSHGEOM_Bounds b;
  EXPECT_EQ(gshgeom_path_bounds(path, &b), GSHGEOM_OK);
  if (b.empty) {
    b.min = b.max = pt(0, 0);
    for (const so::Sub & s : cl.subs) {
      b.min.x = std::min(b.min.x, (double)s.start.x);
      b.min.y = std::min(b.min.y, (double)s.start.y);
      b.max.x = std::max(b.max.x, (double)s.start.x);
      b.max.y = std::max(b.max.y, (double)s.start.y);
    }
  }
  for (const so::Sub & s : cl.subs) {
    b.min.x = std::min(b.min.x, (double)s.start.x);
    b.min.y = std::min(b.min.y, (double)s.start.y);
    b.max.x = std::max(b.max.x, (double)s.start.x);
    b.max.y = std::max(b.max.y, (double)s.start.y);
  }
  so::Real h = sc.width / 2;
  so::Real reach = space == 0 ? h + 1 : h / std::max((so::Real)1e-9L, std::fabs(map.det()) / map.frobenius()) + 1;
  /* device-space box of the model box outset by reach */
  so::V2 corners[4] = {{b.min.x - reach, b.min.y - reach}, {b.max.x + reach, b.min.y - reach},
      {b.max.x + reach, b.max.y + reach}, {b.min.x - reach, b.max.y + reach}};
  so::Real x0 = 1e300L, y0 = 1e300L, x1 = -1e300L, y1 = -1e300L;
  for (so::V2 c : corners) {
    so::V2 d = map.apply(c);
    x0 = std::min(x0, d.x);
    y0 = std::min(y0, d.y);
    x1 = std::max(x1, d.x);
    y1 = std::max(y1, d.y);
  }
  const so::Real m = 1e-6L;
  Agreement ag;
  auto check = [&](so::V2 xd) {
    so::Real dist;
    if (space == 0) {
      dist = centerline_distance(cl, map.inverse(xd));
    } else {
      dist = centerline_distance(mapped, xd);
    }
    if (std::fabs(dist - h) <= m) {
      return;
    }
    so::Verdict v = so::region_contains(sk.body, xd);
    if (v.ambiguous) {
      return;
    }
    ag.decided++;
    if (v.inside != (dist <= h)) {
      if (ag.disagree++ == 0) {
        char buf[200];
        std::snprintf(buf, sizeof buf, "%s at (%.9Lg, %.9Lg): oracle %d, distance %.9Lg (h %.9Lg)", sc.name, xd.x, xd.y,
            (int)v.inside, dist, h);
        ag.first = buf;
      }
    }
  };
  const int N = 60;
  for (int i = 0; i <= N; i++) {
    for (int j = 0; j <= N; j++) {
      check({x0 + (x1 - x0) * i / N, y0 + (y1 - y0) * j / N});
    }
  }
  std::mt19937_64 rng(0x5eedULL);
  std::uniform_real_distribution<double> ux((double)x0, (double)x1), uy((double)y0, (double)y1);
  for (int k = 0; k < 1000; k++) {
    check({ux(rng), uy(rng)});
  }
  gshgeom_path_destroy(path);
  return ag;
}

so::Affine skew_map() {
  so::Affine f;
  f.a = 1.5;
  f.b = 0.75;
  f.c = -0.25;
  f.d = 0.5;
  f.tx = 3;
  f.ty = -2;
  return f;
}

} // namespace

TEST(OracleSelfCheck, ADashThatStartsOrEndsAtACornerHasNoPieceOnTheOtherSide) {
  /* L-shaped path of two segments of 12; dashes cut at exactly the corner.  Segment lengths are summed in long double,
   * so on a wide one the corner is at 12 - 1e-33 and bisecting that gives a parameter of 1e-25: a piece of that length
   * is kept in 128 bits and collapses in 80, which put a join (or a dash piece) at the corner on some targets only. */
  GSHGEOM_Path * path = stroke_matrix::build({{'M', {0, 0}}, {'L', {12, 0}}, {'L', {12, 12}}});
  ASSERT_NE(path, nullptr);
  so::Centerline cl = so::read_centerline(path);
  for (std::vector<so::Real> dashes : {std::vector<so::Real>{12, 4}, std::vector<so::Real>{8, 4}}) {
    so::Style st;
    st.width = 3;
    st.dashes = dashes;
    std::vector<so::V2> tangents(cl.subs.size(), so::V2{1, 0});
    std::vector<so::Sub> subs = so::dash_subs(cl.subs, st, tangents);
    ASSERT_FALSE(subs.empty());
    for (const so::Sub & d : subs) {
      EXPECT_FALSE(d.dot) << dashes[0];
      EXPECT_EQ(d.segs.size(), 1u) << "a dash with a piece on each side of the corner is not one the corner starts or ends";
      EXPECT_TRUE(d.joins.empty()) << dashes[0];
    }
  }
  gshgeom_path_destroy(path);
}

TEST(OracleSelfCheck, ABernsteinRootIsFoundWhenTheIntervalStartsAtAnotherRoot) {
  /* p(t) = t^2 - 0.7 t: a root at 0 and one at 0.7.  The front coefficient is exactly zero (on a wide long double a
   * split at a root leaves such a coefficient, and on a narrow one its noise), so the bisection must take its sign from
   * the first coefficient beyond the noise. */
  stroke_oracle::Bern p = {0.0L, -0.35L, 0.3L};
  std::vector<stroke_oracle::Real> roots;
  bool ambiguous = false;
  stroke_oracle::bern_roots_rec(p, 0, 1, 0, roots, &ambiguous);
  EXPECT_FALSE(ambiguous);
  bool found = false;
  for (stroke_oracle::Real r : roots) {
    found = found || std::fabs(r - 0.7L) < 1e-15L;
  }
  EXPECT_TRUE(found) << "roots found: " << roots.size();
  /* and the same polynomial with noise of either sign in the front coefficient */
  for (stroke_oracle::Real noise : {1e-36L, -1e-36L}) {
    stroke_oracle::Bern q = {noise, -0.35L, 0.3L};
    roots.clear();
    stroke_oracle::bern_roots_rec(q, 0, 1, 0, roots, &ambiguous);
    found = false;
    for (stroke_oracle::Real r : roots) {
      found = found || std::fabs(r - 0.7L) < 1e-15L;
    }
    EXPECT_TRUE(found) << noise;
  }
}

TEST(OracleSelfCheck, RoundJoinsAndCapsAreTheDistanceToTheCenterline) {
  GSHGEOM_SKIP_UNDER_VALGRIND();
  for (const SelfCase & sc : self_cases()) {
    for (int space = 0; space <= 1; space++) {
      for (const so::Affine & map : {so::Affine{}, skew_map()}) {
        Agreement ag = self_check(sc, map, space, 0);
        EXPECT_GT(ag.decided, 1000) << sc.name;
        EXPECT_EQ(ag.disagree, 0) << ag.first;
      }
    }
  }
}

TEST(OracleSelfCheck, AnInnerSideSectorIsSeenToFail) {
  /* The oracle's own planted defect: the round join's sector on the inner side.  The self-check must see it on a
   * polyline whose inner corner is not covered by the bodies. */
  long seen = 0;
  for (const SelfCase & sc : self_cases()) {
    if (std::string(sc.name).rfind("polyline", 0) != 0) {
      continue; /* the polylines have the inner corners a sector on the wrong side shows */
    }
    Agreement ag = self_check(sc, so::Affine{}, 0, 1);
    seen += ag.disagree;
  }
  EXPECT_GT(seen, 0) << "the self-check cannot see a sector on the inner side";
}

/* --- the exact vectors --------------------------------------------------------------- */

namespace {

so::Style vector_style(const stroke_vectors::Case & c) {
  so::Style st;
  st.width = c.width;
  st.space = c.space;
  st.join = c.join;
  st.limit = c.limit;
  st.cap = c.cap;
  st.degenerate = c.degenerate;
  return st;
}

so::Affine vector_map(const stroke_vectors::Case & c) {
  so::Affine f;
  f.a = c.map[0];
  f.b = c.map[1];
  f.c = c.map[2];
  f.d = c.map[3];
  f.tx = c.map[4];
  f.ty = c.map[5];
  return f;
}

GSHGEOM_Path * vector_path(const stroke_vectors::Case & c) {
  Ops ops;
  for (size_t i = 0; i < c.op_count; i++) {
    const stroke_vectors::Op & o = c.ops[i];
    std::vector<double> a(o.v, o.v + o.n);
    ops.push_back({o.verb, a});
  }
  return build(ops);
}

} // namespace

TEST(OracleExactVectors, EveryVerdictMatches) {
  size_t checked = 0;
  for (size_t i = 0; i < stroke_vectors::case_count; i++) {
    const stroke_vectors::Case & c = stroke_vectors::cases[i];
    GSHGEOM_Path * path = vector_path(c);
    so::Centerline cl = so::read_centerline(path);
    so::Stroked sk = so::build_stroke(cl, vector_style(c), vector_map(c), c.alignment);
    for (size_t k = 0; k < c.point_count; k++) {
      const stroke_vectors::Point & p = c.points[k];
      const so::Region & r = p.part == 0 ? sk.body : sk.aligned;
      so::Verdict v = so::region_contains(r, {p.x, p.y});
      EXPECT_FALSE(v.ambiguous) << c.name << " point " << k;
      EXPECT_EQ((int)v.inside, p.inside) << c.name << " point " << k << " (" << p.x << ", " << p.y << ")";
      checked++;
    }
    gshgeom_path_destroy(path);
  }
  EXPECT_GT(checked, (size_t)100);
}

TEST(OracleExactVectors, TheArcsJoinsConstructionsMatch) {
  /* X, the clip point Y and its tangent, for each arcs case, against the exact model. */
  size_t checked = 0;
  for (size_t i = 0; i < stroke_vectors::arcs_count; i++) {
    const stroke_vectors::Arcs & a = stroke_vectors::arcs[i];
    const stroke_vectors::Case & c = stroke_vectors::cases[a.case_index];
    GSHGEOM_Path * path = vector_path(c);
    so::Centerline cl = so::read_centerline(path);
    so::Stroked sk = so::build_stroke(cl, vector_style(c), vector_map(c), c.alignment);
    const so::Shape * found = nullptr;
    size_t seen = 0;
    for (const so::Shape & s : sk.body.shapes) {
      if (s.kind == so::Shape::ARCS) {
        if (seen++ == a.join_index) {
          found = &s;
        }
      }
    }
    ASSERT_NE(found, nullptr) << c.name << " (" << a.branch << ")";
    SCOPED_TRACE(std::string(c.name) + " join " + std::to_string(a.join_index) + " " + a.branch);
    const so::Real tol = 1e-9L;
    EXPECT_NEAR((double)found->X.x, a.x[0], tol * (1 + std::fabs(a.x[0]))) << c.name;
    EXPECT_NEAR((double)found->X.y, a.x[1], tol * (1 + std::fabs(a.x[1]))) << c.name;
    EXPECT_EQ((int)found->clip, a.clipped) << c.name;
    if (a.clipped && found->clip) {
      EXPECT_NEAR((double)found->clip_y.x, a.y[0], tol * (1 + std::fabs(a.y[0]))) << c.name;
      EXPECT_NEAR((double)found->clip_y.y, a.y[1], tol * (1 + std::fabs(a.y[1]))) << c.name;
      EXPECT_NEAR((double)found->clip_t.x, a.t[0], tol) << c.name;
      EXPECT_NEAR((double)found->clip_t.y, a.t[1], tol) << c.name;
    }
    checked++;
    gshgeom_path_destroy(path);
  }
  EXPECT_GT(checked, (size_t)8);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
