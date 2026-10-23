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
 * The expansion against the geometric oracle over the stroke's case matrix
 * (stroke.h, SG17): lines and every join, the miter boundary, joins against caps,
 * every curve kind (cusps, loops, circles with h below, at and above the
 * radius), zero extent, degenerate subpaths, lone moves, both width spaces under
 * every kind of map, alignment, overlaps, and every tolerance.  For each the
 * outline's nonzero fill agrees with the region at every decided point, no
 * subpath winds negatively, every vertex is inside the conservative box plus the
 * tolerance, and the deviation profile's maximum is recorded as a test property
 * and is at most the tolerance; and the exact vectors of tools/gen-stroke-vectors.py
 * are classified as the exact model classifies them.
 */

#include "../test_helpers.h"

#include "../stroke_cases.h"
#include "../stroke_matrix.h"
#include "../stroke_vectors.h"

#include <string>
#include <vector>

namespace sm = stroke_matrix;
namespace so = stroke_oracle;
namespace sc = stroke_cases;

namespace {

using namespace stroke_cases;

/** The worst deviation seen over every case, in 26.6 units. */
double g_worst = 0.0;
double g_worst_margin_tau = 0.0;

void run(const Case & c) {
  SCOPED_TRACE(c.name);
  GSHGEOM_Stroke st = c.st;
  if (!c.dashes.empty()) {
    st.dashes = c.dashes.data();
    st.dash_count = c.dashes.size();
  }
  GSHGEOM_Path * path = sm::build(c.ops);
  ASSERT_NE(path, nullptr);
  GSHGEOM_Stroke_Outline * o = nullptr;
  ASSERT_EQ(gshgeom_stroke_expand_with_tolerance(path, &st, &c.map, c.tau, nullptr, nullptr, &o), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_stroke_outline_status(o), GSHGEOM_STROKE_DRAWN);
  so::Centerline cl = so::read_centerline(path);
  so::Stroked sk = so::build_stroke(cl, sm::oracle_style(st), so::affine_of(c.map), (int)st.alignment);
  sm::Box box = sm::stroke_box(path, st, c.map);
  struct Part {
    const GSHGEOM_Flat_Path * fp;
    const so::Region * r;
    const char * name;
  } parts[2] = {{gshgeom_stroke_outline_body(o), &sk.body, "body"},
      {gshgeom_stroke_outline_aligned(o), &sk.aligned, "aligned"}};
  EXPECT_EQ(gshgeom_stroke_outline_aligned(o) != nullptr, sk.has_aligned);
  for (const Part & p : parts) {
    if (!p.fp) {
      continue;
    }
    sm::PartResult pr = sm::compare_part(p.fp, *p.r, box, c.tau, 2500, 0x51ULL);
    EXPECT_GT(pr.decided, 0) << p.name;
    EXPECT_EQ(pr.wrong, 0) << p.name << ": " << pr.first;
    EXPECT_EQ(pr.negative, 0) << p.name;
    EXPECT_EQ(pr.outside_box, 0) << p.name;
    double dev_units = (double)(pr.worst_outward * 64);
    EXPECT_LE(dev_units, (double)c.tau) << p.name;
    if (dev_units > g_worst) {
      g_worst = dev_units;
    }
    if (dev_units / c.tau > g_worst_margin_tau) {
      g_worst_margin_tau = dev_units / c.tau;
    }
  }
  gshgeom_stroke_outline_destroy(o);
  gshgeom_path_destroy(path);
}

} // namespace

TEST(StrokeMatrix, LinesEveryJoinAndCap) {
  for (const Case & c : sc::linesEveryJoinAndCap()) {
    run(c);
  }
}

TEST(StrokeMatrix, TheMiterBoundaryIsDecidedExactly) {
  for (const Case & c : sc::theMiterBoundaryIsDecidedExactly()) {
    run(c);
  }
}

TEST(StrokeMatrix, EveryCurveKind) {
  for (const Case & c : sc::everyCurveKind()) {
    run(c);
  }
}

TEST(StrokeMatrix, ZeroExtentAndDegenerateSubpaths) {
  std::vector<Ops> paths = {
      {{'M', {2, 2}}, {'L', {2, 2}}},
      {{'M', {2, 2}}, {'Z', {}}},
      {{'M', {2, 2}}, {'C', {2, 2, 2, 2, 2, 2}}},
      {{'M', {0, 0}}, {'L', {6, 0}}, {'L', {6, 0}}, {'L', {6, 6}}},
      {{'M', {0, 0}}, {'L', {3, 4}}, {'M', {10, 10}}, {'Z', {}}},
      {{'M', {10, 10}}, {'L', {10, 10}}, {'M', {0, 0}}, {'L', {4, -3}}},
      {{'M', {0, 0}}, {'L', {5, 5}}, {'M', {3, 3}}, {'M', {10, 10}}, {'L', {10, 10}}},
      {{'M', {0, 0}}, {'L', {6, 0}}, {'L', {6, 6}}, {'Z', {}}, {'M', {20, 20}}, {'Z', {}}},
      {{'M', {4, 4}}},
  };
  std::vector<GSHGEOM_Map> maps = {map_id(), map_of(0.6, -0.8, 0.8, 0.6, 1, 2), map_of(1, 1, 0, 1, 0, 0)};
  for (const Ops & p : paths) {
    for (int cap : {BUTT, RCAP, SQUARE}) {
      for (uint32_t deg : {0u, 1u}) {
        for (uint32_t space : {0u, 1u}) {
          for (const GSHGEOM_Map & m : maps) {
            if ((cap != SQUARE || deg != 0) && (&m != &maps[0] || space != 0)) {
              continue; /* every cap and both settings under the identity; the square under every map and space */
            }
            GSHGEOM_Stroke s = style(3, MITER, 4, cap);
            s.degenerate = deg;
            s.width_space = space;
            GSHGEOM_Path * pp = sm::build(p);
            so::Centerline cl = so::read_centerline(pp);
            gshgeom_path_destroy(pp);
            bool draws = false;
            for (const so::Sub & sb : cl.subs) {
              draws = draws || (sb.draws && (!sb.segs.empty() || (deg == 0 && cap != BUTT)));
            }
            if (!draws) {
              /* nothing is drawn: the body is empty */
              GSHGEOM_Path * q = sm::build(p);
              GSHGEOM_Stroke_Outline * o = nullptr;
              ASSERT_EQ(gshgeom_stroke_expand(q, &s, &m, nullptr, nullptr, &o), GSHGEOM_OK);
              EXPECT_EQ(gshgeom_flat_path_subpath_count(gshgeom_stroke_outline_body(o)), 0u);
              gshgeom_stroke_outline_destroy(o);
              gshgeom_path_destroy(q);
              continue;
            }
            run({"degenerate", p, s, m, GSHGEOM_FLATTEN_TOLERANCE});
          }
        }
      }
    }
  }
}

TEST(StrokeMatrix, BothWidthSpacesUnderEveryMap) {
  for (const Case & c : sc::bothWidthSpacesUnderEveryMap()) {
    run(c);
  }
}

TEST(StrokeMatrix, Alignment) {
  for (const Case & c : sc::alignment()) {
    run(c);
  }
}

TEST(StrokeMatrix, OverlapsAreOneUnion) {
  for (const Case & c : sc::overlapsAreOneUnion()) {
    run(c);
  }
}

TEST(StrokeMatrix, EveryTolerance) {
  for (const Case & c : sc::everyTolerance()) {
    run(c);
  }
}

TEST(StrokeMatrix, ArcsJoins) {
  /* the exact vectors' arcs cases, through the expansion, at every tolerance the matrix uses */
  for (size_t i = 0; i < stroke_vectors::case_count; i++) {
    const stroke_vectors::Case & vc = stroke_vectors::cases[i];
    if (vc.join != ARCS) {
      continue;
    }
    Ops ops;
    for (size_t k = 0; k < vc.op_count; k++) {
      ops.push_back({vc.ops[k].verb, std::vector<double>(vc.ops[k].v, vc.ops[k].v + vc.ops[k].n)});
    }
    GSHGEOM_Stroke s = style(vc.width, vc.join, vc.limit, vc.cap);
    s.width_space = (uint32_t)vc.space;
    for (int32_t tau : {2, 4}) {
      run({vc.name, ops, s, map_of(vc.map[0], vc.map[1], vc.map[2], vc.map[3], vc.map[4], vc.map[5]), tau});
    }
  }
}

TEST(StrokeMatrix, TheExactVectorsThroughTheExpansion) {
  /* every exact verdict, through the outline: the points are 0.2 px from every boundary the exact model drew, beyond
   * the default tolerance's reach */
  size_t checked = 0;
  for (size_t i = 0; i < stroke_vectors::case_count; i++) {
    const stroke_vectors::Case & vc = stroke_vectors::cases[i];
    Ops ops;
    for (size_t k = 0; k < vc.op_count; k++) {
      ops.push_back({vc.ops[k].verb, std::vector<double>(vc.ops[k].v, vc.ops[k].v + vc.ops[k].n)});
    }
    GSHGEOM_Stroke s = style(vc.width, vc.join, vc.limit, vc.cap);
    s.width_space = (uint32_t)vc.space;
    s.degenerate = (uint32_t)vc.degenerate;
    s.alignment = (uint32_t)vc.alignment;
    GSHGEOM_Map m = map_of(vc.map[0], vc.map[1], vc.map[2], vc.map[3], vc.map[4], vc.map[5]);
    GSHGEOM_Path * path = sm::build(ops);
    ASSERT_NE(path, nullptr);
    GSHGEOM_Stroke_Outline * o = nullptr;
    ASSERT_EQ(gshgeom_stroke_expand(path, &s, &m, nullptr, nullptr, &o), GSHGEOM_OK) << vc.name;
    so::Winding body, aligned;
    body.build(so::polygons_of(gshgeom_stroke_outline_body(o)));
    aligned.build(so::polygons_of(gshgeom_stroke_outline_aligned(o)));
    for (size_t k = 0; k < vc.point_count; k++) {
      const stroke_vectors::Point & pt = vc.points[k];
      int w = (pt.part == 0 ? body : aligned).at({pt.x, pt.y});
      EXPECT_EQ(w != 0, pt.inside != 0) << vc.name << " point " << k << " (" << pt.x << ", " << pt.y << ")";
      checked++;
    }
    gshgeom_stroke_outline_destroy(o);
    gshgeom_path_destroy(path);
  }
  EXPECT_GT(checked, (size_t)1000);
}


TEST(StrokeMatrix, DashesOnLines) {
  for (const Case & c : sc::dashesOnLines()) {
    run(c);
  }
}

TEST(StrokeMatrix, DashesAtVerticesAndDots) {
  for (const Case & c : sc::dashesAtVerticesAndDots()) {
    run(c);
  }
}

TEST(StrokeMatrix, DashesAcrossTheStartOfAClosedSubpath) {
  for (const Case & c : sc::dashesAcrossTheStartOfAClosedSubpath()) {
    run(c);
  }
}

TEST(StrokeMatrix, ASmoothClosedCurveUnderADeviceMapHasNoArcsJoinAtItsSmoothVertices) {
  /* A lens whose two cubics meet smoothly at both ends (the tangents are (-8, 0) and (-8, 0) at (12, 16), (8, 0) and
   * (8, 0) at (12, 4)).  Under a rotation by 0.9/0.4 the mapped control points no longer line up exactly, the vertex is
   * a turn of a few ulps on the doubles, and under device width the arcs join was built from it (a filled disc where the
   * region is a ring).  SG8 gives a smooth vertex no join: the region must not depend on rounding of the mapped doubles. */
  Ops lens = {{'M', {12, 4}}, {'C', {20, 4, 20, 16, 12, 16}}, {'C', {4, 16, 4, 4, 12, 4}}, {'Z', {}}};
  for (const GSHGEOM_Map & m : {map_of(0.9, -0.4, 0.4, 0.9, 1.3125, 0.6875), map_of(0.9, 0.4, -0.4, 0.9, 1.3125, 0.6875),
           map_of(0.95, 0.2, -0.2, 0.95, 0.375, 0.5), map_id()}) {
    for (int join : {ARCS, MITER, ROUND, BEVEL, MCLIP}) {
      for (uint32_t space : {1u, 0u}) {
        GSHGEOM_Stroke st = style(3, join, 4, BUTT);
        st.width_space = space;
        run({"lens, a smooth closed curve under a map", lens, st, m, GSHGEOM_FLATTEN_TOLERANCE});
        GSHGEOM_Path * path = sm::build(lens);
        GSHGEOM_Stroke_Outline * o = nullptr;
        ASSERT_EQ(gshgeom_stroke_expand(path, &st, &m, nullptr, nullptr, &o), GSHGEOM_OK);
        /* every join of the same case is about 3000 points; the defect made 7647 */
        EXPECT_LT(gshgeom_flat_path_point_count(gshgeom_stroke_outline_body(o)), (size_t)4500) << "join " << join << " space " << space;
        gshgeom_stroke_outline_destroy(o);
        gshgeom_path_destroy(path);
      }
    }
  }
}

TEST(StrokeMatrix, DashesOnCurves) {
  for (const Case & c : sc::dashesOnCurves()) {
    run(c);
  }
}

TEST(StrokeMatrix, DashesUnderBothWidthSpaces) {
  for (const Case & c : sc::dashesUnderBothWidthSpaces()) {
    run(c);
  }
}

TEST(StrokeMatrix, DashedAlignment) {
  for (const Case & c : sc::dashedAlignment()) {
    run(c);
  }
}

TEST(StrokeMatrix, DashEdgesOnACircleAreAtTheExactArcLength) {
  /* a circle of radius 20 stroked at width 2 with butt caps and the pattern [10, 5]: the dash k covers the arc
   * [15 k, 15 k + 10], so the centerline point at the arc length s is inside exactly when s mod 15 < 10; probed 0.2
   * beyond every edge (past the tolerance) at the centerline's radius */
  Ops c = circle(0, 0, 20);
  GSHGEOM_Path * path = sm::build(c);
  GSHGEOM_Stroke st = style(2, MITER, 4, BUTT);
  double d[2] = {10, 5};
  st.dashes = d;
  st.dash_count = 2;
  GSHGEOM_Map m = map_id();
  GSHGEOM_Stroke_Outline * o = nullptr;
  ASSERT_EQ(gshgeom_stroke_expand(path, &st, &m, nullptr, nullptr, &o), GSHGEOM_OK);
  so::Winding w;
  w.build(so::polygons_of(gshgeom_stroke_outline_body(o)));
  double total = 2 * 3.14159265358979323846 * 20;
  int checked = 0;
  for (double s0 = 0; s0 < total - 0.5; s0 += 5.0) {
    for (double delta : {-0.2, 0.2}) {
      double edge = s0;
      double s = edge + delta;
      if (s < 0 || s > total - 0.3) {
        continue;
      }
      double mod = std::fmod(s, 15.0);
      /* only points 0.2 from an edge: edges are at multiples of 5 */
      bool inside = mod < 10.0;
      double th = s / 20.0;
      so::V2 pt = {20 * std::cos(th), 20 * std::sin(th)};
      EXPECT_EQ(w.at(pt) != 0, inside) << "arc length " << s;
      checked++;
    }
  }
  EXPECT_GT(checked, 20);
  gshgeom_stroke_outline_destroy(o);
  gshgeom_path_destroy(path);
}

TEST(StrokeMatrix, DashesAreRefusedPastTheErrorGate) {
  /* SG26: a cut error past the gate is ERR_LIMIT, never a coarser cut: 2^-46 of a polygon 2^29 long is 2^-17 */
  GSHGEOM_Path * path = sm::build(Ops{{'M', {0, 0}}, {'L', {0x1p29, 0}}});
  GSHGEOM_Stroke st = style(2, MITER, 4, BUTT);
  double d[2] = {0x1p22, 0x1p22}; /* a few dozen dashes: only the cut position's error bound refuses this */
  st.dashes = d;
  st.dash_count = 2;
  GSHGEOM_Map m = map_id();
  GSHGEOM_Stroke_Outline * o = nullptr;
  EXPECT_EQ(gshgeom_stroke_expand(path, &st, &m, nullptr, nullptr, &o), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(o, nullptr);
  gshgeom_path_destroy(path);
}

TEST(StrokeMatrix, TheDeviationProfileIsRecorded) {
  /* runs last in this file (gtest runs tests in order of definition): the worst outward deviation of any outline
   * vertex from its region over every case above, as a property and against the proof's margin */
  RecordProperty("worst_outward_deviation_units", std::to_string(g_worst));
  RecordProperty("worst_deviation_over_tolerance", std::to_string(g_worst_margin_tau));
  std::printf("stroke matrix: worst outward deviation %.4f units (%.4f of the tolerance)\n", g_worst,
      g_worst_margin_tau);
  EXPECT_GT(g_worst, 0.0) << "nothing was measured";
  EXPECT_LE(g_worst_margin_tau, 1.0);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
