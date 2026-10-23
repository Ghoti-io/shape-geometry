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
 * The curve outline (stroke_curves.h) against the geometric oracle over the stroke case matrix of
 * tests/stroke_cases.h (SG33): at every decided point (farther than `t + m` from the boundary) the nonzero fill of the
 * curves agrees with the region, no ring winds negatively, the outline lies inside the conservative box plus `t`, and the
 * deviation profile's maximum is recorded as a test property and is at most `t`.  The form `CUBIC` is the form
 * `RATIONAL` lowered by ::gshgeom_stroke_curves_lower_cubic, bit for bit (SG40).  And the agreement with stroke.h's
 * polyline (SG42): flattened at the polyline's tolerance `tau`, the curves and the polyline agree within `2 tau + t`, and
 * the worst measured disagreement is recorded and is at most `2 tau`.
 */

#include "../test_helpers.h"

#include "../stroke_cases.h"
#include "../stroke_curves_matrix.h"
#include "../stroke_vectors.h"

#include <chrono>
#include <random>
#include <string>
#include <vector>

namespace sm = stroke_matrix;
namespace so = stroke_oracle;
namespace scm = stroke_curves_matrix;
namespace sc = stroke_cases;
namespace ce = curve_eval;

namespace {

/** The worst deviation over every case, as a fraction of the tolerance. */
double g_worst_dev = 0.0;
/** The same, over the cases run with boundary samples t / 8 apart: the deviation at t's own resolution. */
double g_worst_fine = 0.0;
/** The worst disagreement with the polyline, in units of the polyline's tolerance (SG42). */
double g_worst_poly = 0.0;
double g_worst_poly_over_bound = 0.0;
double g_tm[5];
std::vector<std::pair<double, std::string>> g_slow;
double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

bool has_rational(const GSHGEOM_Path * p) {
  size_t n = gshgeom_path_element_count(p);
  for (size_t i = 0; i < n; i++) {
    GSHGEOM_Element e;
    gshgeom_path_element(p, i, &e);
    if (e.verb == GSHGEOM_VERB_RATIONAL_QUADRATIC || e.verb == GSHGEOM_VERB_RATIONAL_CUBIC ||
        e.verb == GSHGEOM_VERB_QUADRATIC) {
      return true;
    }
  }
  return false;
}

/** The distance from x to the boundary of a set of rings, in pixels. */
double dist_to_rings(const so::Polygons & polys, so::V2 x) {
  so::Real best = 1e300L;
  for (const auto & ring : polys.rings) {
    size_t n = ring.size();
    for (size_t i = 0; i < n; i++) {
      ce::P2 a = {(ce::Real)ring[i].first / 64, (ce::Real)ring[i].second / 64};
      ce::P2 b = {(ce::Real)ring[(i + 1) % n].first / 64, (ce::Real)ring[(i + 1) % n].second / 64};
      best = std::min(best, ce::dist_to_segment({x.x, x.y}, a, b));
    }
  }
  return (double)best;
}

/** The curve outline's body, flattened at tau (SG42): the polyline `A`. */
GSHGEOM_Flat_Path * flatten_curves(const GSHGEOM_Path * body, int32_t tau) {
  GSHGEOM_Map id;
  gshgeom_map_identity(2, &id);
  GSHGEOM_Device_Path * dp = nullptr;
  EXPECT_EQ(gshgeom_device_path_from_model(body, &id, nullptr, nullptr, &dp), GSHGEOM_OK);
  GSHGEOM_Flat_Path * fp = nullptr;
  EXPECT_EQ(gshgeom_device_path_flatten_with_tolerance(dp, tau, nullptr, nullptr, &fp), GSHGEOM_OK);
  gshgeom_device_path_destroy(dp);
  return fp;
}

void check_polyline_agreement(const sc::Case & c, const GSHGEOM_Stroke & st, const GSHGEOM_Path * path,
    const GSHGEOM_Path * curves_body, double t) {
  int32_t tau = GSHGEOM_FLATTEN_TOLERANCE;
  GSHGEOM_Stroke_Outline * so_out = nullptr;
  ASSERT_EQ(gshgeom_stroke_expand_with_tolerance(path, &st, &c.map, tau, nullptr, nullptr, &so_out), GSHGEOM_OK);
  GSHGEOM_Flat_Path * A = flatten_curves(curves_body, tau);
  so::Polygons pa = so::polygons_of(A), pb = so::polygons_of(gshgeom_stroke_outline_body(so_out));
  if (pa.rings.empty() && pb.rings.empty()) {
    gshgeom_flat_path_destroy(A);
    gshgeom_stroke_outline_destroy(so_out);
    return;
  }
  so::Winding wa, wb;
  wa.build(pa);
  wb.build(pb);
  double bound = 2.0 * tau / 64.0 + t; /* SG42's proved value, in pixels */
  double brief = 2.0 * tau / 64.0;     /* the brief's number, a recorded property */
  /* sample at the vertices of both polygons, displaced, and on a coarse grid of their box */
  std::mt19937_64 rng(11);
  std::vector<so::V2> pts;
  for (const so::Polygons * pp : {&pa, &pb}) {
    size_t stride = 1;
    size_t total = 0;
    for (auto & r : pp->rings) {
      total += r.size();
    }
    stride = total / 1500 + 1;
    size_t k = 0;
    for (auto & r : pp->rings) {
      for (auto & q : r) {
        if (k++ % stride) {
          continue;
        }
        std::uniform_real_distribution<double> u(-0.3, 0.3);
        pts.push_back({(so::Real)q.first / 64 + u(rng), (so::Real)q.second / 64 + u(rng)});
      }
    }
  }
  for (const so::V2 & x : pts) {
    bool ina = wa.at(x) != 0, inb = wb.at(x) != 0;
    if (ina == inb) {
      continue;
    }
    double da = dist_to_rings(pa, x);
    EXPECT_LE(da, bound) << c.name << " at (" << (double)x.x << ", " << (double)x.y << ")";
    g_worst_poly = std::max(g_worst_poly, da * 64.0 / tau);
    g_worst_poly_over_bound = std::max(g_worst_poly_over_bound, da / bound);
    EXPECT_LE(da, brief) << c.name << ": the brief's twice the tolerance";
  }
  gshgeom_flat_path_destroy(A);
  gshgeom_stroke_outline_destroy(so_out);
}

void run(const sc::Case & c, double t, bool polyline, double sigma_over_t = 1.0) {
  SCOPED_TRACE(c.name);
  GSHGEOM_Stroke st = c.st;
  if (!c.dashes.empty()) {
    st.dashes = c.dashes.data();
    st.dash_count = c.dashes.size();
  }
  GSHGEOM_Path * path = sm::build(c.ops);
  ASSERT_NE(path, nullptr);
  double t0 = now();
  GSHGEOM_Stroke_Curves * o = nullptr;
  ASSERT_EQ(gshgeom_stroke_curves_expand(path, &st, &c.map, t, GSHGEOM_STROKE_CURVES_RATIONAL, nullptr, nullptr, &o),
      GSHGEOM_OK);
  ASSERT_EQ(gshgeom_stroke_curves_status(o), GSHGEOM_STROKE_DRAWN);
  g_tm[0] += now() - t0;
  t0 = now();
  so::Centerline cl = so::read_centerline(path);
  so::Stroked sk = so::build_stroke(cl, sm::oracle_style(st), so::affine_of(c.map), (int)st.alignment);
  sm::Box box = sm::stroke_box(path, st, c.map);
  if (st.join == GSHGEOM_STROKE_JOIN_MITER_CLIP || st.join == GSHGEOM_STROKE_JOIN_ARCS) {
    /* the scene's conservative box outsets the path by L h; a miter-clip corner lies at distance up to h sqrt(L^2 + 1)
     * from the vertex (the clip line is L h along the bisector and the corner up to h across it), so it can pass the
     * box by up to h (sqrt(L^2 + 1) - L): a finding recorded in the story's change log */
    double L = st.miter_limit;
    double slack = 0.5 * st.width * (st.alignment == GSHGEOM_STROKE_ALIGN_CENTER ? 1.0 : 2.0) *
        (std::sqrt(L * L + 1.0) - L);
    if (st.width_space == GSHGEOM_STROKE_WIDTH_LOCAL) {
      slack *= (double)so::affine_of(c.map).frobenius();
    }
    box.x0 -= slack;
    box.y0 -= slack;
    box.x1 += slack;
    box.y1 += slack;
  }
  struct Part {
    const GSHGEOM_Path * p;
    const so::Region * r;
    const char * name;
  } parts[2] = {{gshgeom_stroke_curves_body(o), &sk.body, "body"}, {gshgeom_stroke_curves_aligned(o), &sk.aligned, "aligned"}};
  EXPECT_EQ(gshgeom_stroke_curves_aligned(o) != nullptr, sk.has_aligned);
  for (const Part & p : parts) {
    if (!p.p) {
      continue;
    }
    scm::PartResult pr = scm::compare_curves(p.p, *p.r, box, t, 2500, 0x51ULL, sigma_over_t);
    EXPECT_GT(pr.decided, 0) << p.name;
    EXPECT_EQ(pr.wrong, 0) << p.name << ": " << pr.first;
    EXPECT_EQ(pr.negative, 0) << p.name;
    EXPECT_EQ(pr.outside_box, 0) << p.name;
    double dev = (double)pr.worst_outward / t;
    EXPECT_LE(dev, 1.0 + sigma_over_t + 0.01) << p.name; /* the oracle's boundary samples are sigma_over_t t apart */
    (sigma_over_t < 0.5 ? g_worst_fine : g_worst_dev) = std::max(sigma_over_t < 0.5 ? g_worst_fine : g_worst_dev, dev);
  }
  g_tm[1] += now() - t0;
  g_slow.push_back({now() - t0, c.name + " t=" + std::to_string(t) + " elements=" + std::to_string(gshgeom_path_element_count(gshgeom_stroke_curves_body(o)))});
  t0 = now();
  /* SG40: the cubic form is the rational form lowered at a quarter of the tolerance */
  GSHGEOM_Stroke_Curves * oc = nullptr;
  ASSERT_EQ(gshgeom_stroke_curves_expand(path, &st, &c.map, t, GSHGEOM_STROKE_CURVES_CUBIC, nullptr, nullptr, &oc),
      GSHGEOM_OK);
  for (int k = 0; k < 2; k++) {
    const GSHGEOM_Path * rat = k == 0 ? gshgeom_stroke_curves_body(o) : gshgeom_stroke_curves_aligned(o);
    const GSHGEOM_Path * cub = k == 0 ? gshgeom_stroke_curves_body(oc) : gshgeom_stroke_curves_aligned(oc);
    ASSERT_EQ(rat != nullptr, cub != nullptr);
    if (!rat) {
      continue;
    }
    GSHGEOM_Path * low = nullptr;
    ASSERT_EQ(gshgeom_stroke_curves_lower_cubic(rat, t / 4, nullptr, nullptr, &low), GSHGEOM_OK);
    EXPECT_EQ(gshgeom_path_hash(low), gshgeom_path_hash(cub));
    EXPECT_FALSE(has_rational(cub));
    gshgeom_path_destroy(low);
  }
  g_tm[2] += now() - t0;
  t0 = now();
  if (polyline && st.alignment == GSHGEOM_STROKE_ALIGN_CENTER) {
    check_polyline_agreement(c, st, path, gshgeom_stroke_curves_body(o), t);
  }
  g_tm[3] += now() - t0;
  gshgeom_stroke_curves_destroy(oc);
  gshgeom_stroke_curves_destroy(o);
  gshgeom_path_destroy(path);
}

} // namespace

TEST(StrokeCurvesOracle, TheCaseMatrixAtTheDefaultTolerance) {
  GSHGEOM_SKIP_UNDER_VALGRIND();
  for (const sc::Case & c : sc::everything()) {
    run(c, GSHGEOM_STROKE_CURVES_TOLERANCE, true);
  }
}

TEST(StrokeCurvesOracle, TheArcsJoinsOfTheExactVectors) {
  /* the exact vectors' arcs cases, through the curve outline, at the default and a coarse tolerance */
  for (size_t i = 0; i < stroke_vectors::case_count; i++) {
    const stroke_vectors::Case & vc = stroke_vectors::cases[i];
    if (vc.join != sc::ARCS) {
      continue;
    }
    sc::Ops ops;
    for (size_t k = 0; k < vc.op_count; k++) {
      ops.push_back({vc.ops[k].verb, std::vector<double>(vc.ops[k].v, vc.ops[k].v + vc.ops[k].n)});
    }
    GSHGEOM_Stroke st = sc::style(vc.width, vc.join, vc.limit, vc.cap);
    st.width_space = (uint32_t)vc.space;
    for (double t : {GSHGEOM_STROKE_CURVES_TOLERANCE, 0.0625}) {
      run({vc.name, ops, st, sc::map_of(vc.map[0], vc.map[1], vc.map[2], vc.map[3], vc.map[4], vc.map[5]), 4}, t, false);
    }
  }
}

TEST(StrokeCurvesOracle, TheExactVectorsThroughTheCurveOutline) {
  /* every exact verdict, through the outline's nonzero fill */
  size_t checked = 0;
  for (size_t i = 0; i < stroke_vectors::case_count; i++) {
    const stroke_vectors::Case & vc = stroke_vectors::cases[i];
    sc::Ops ops;
    for (size_t k = 0; k < vc.op_count; k++) {
      ops.push_back({vc.ops[k].verb, std::vector<double>(vc.ops[k].v, vc.ops[k].v + vc.ops[k].n)});
    }
    GSHGEOM_Stroke st = sc::style(vc.width, vc.join, vc.limit, vc.cap);
    st.width_space = (uint32_t)vc.space;
    st.degenerate = (uint32_t)vc.degenerate;
    st.alignment = (uint32_t)vc.alignment;
    GSHGEOM_Map m = sc::map_of(vc.map[0], vc.map[1], vc.map[2], vc.map[3], vc.map[4], vc.map[5]);
    GSHGEOM_Path * path = sm::build(ops);
    ASSERT_NE(path, nullptr);
    GSHGEOM_Stroke_Curves * o = nullptr;
    ASSERT_EQ(gshgeom_stroke_curves_expand(path, &st, &m, GSHGEOM_STROKE_CURVES_TOLERANCE, GSHGEOM_STROKE_CURVES_RATIONAL,
                  nullptr, nullptr, &o), GSHGEOM_OK) << vc.name;
    ce::Winding body, aligned;
    body.build(gshgeom_stroke_curves_body(o), 1e-5L);
    aligned.build(gshgeom_stroke_curves_aligned(o), 1e-5L);
    for (size_t k = 0; k < vc.point_count; k++) {
      const stroke_vectors::Point & pt = vc.points[k];
      int w = (pt.part == 0 ? body : aligned).at({pt.x, pt.y});
      EXPECT_EQ(w != 0, pt.inside != 0) << vc.name << " point " << k << " (" << pt.x << ", " << pt.y << ")";
      checked++;
    }
    gshgeom_stroke_curves_destroy(o);
    gshgeom_path_destroy(path);
  }
  EXPECT_GT(checked, (size_t)1000);
}

TEST(StrokeCurvesOracle, EveryToleranceOfTheRange) {
  sc::Ops p = {{'M', {0, 0}}, {'C', {8, 12, 14, -6, 20, 4}}, {'L', {12, 10}}};
  for (double t : {GSHGEOM_STROKE_CURVES_MIN_TOLERANCE, 0x1p-12, 0x1p-8, 0.125, 1.0, GSHGEOM_STROKE_CURVES_MAX_TOLERANCE}) {
    for (int join : {sc::MITER, sc::ROUND, sc::ARCS}) {
      sc::Case c{"tolerance " + std::to_string(t), p, sc::style(5, join, 4, sc::RCAP), sc::map_id(), 4};
      run(c, t, false);
    }
  }
}

TEST(StrokeCurvesOracle, TheDeviationAtTheResolutionOfTheTolerance) {
  /* boundary samples t / 8 apart: the deviation profile at t's own resolution, over a representative set */
  sc::Ops cubic = {{'M', {0, 0}}, {'C', {4, 8, 6, -8, 10, 0}}};
  sc::Ops loop = {{'M', {0, 0}}, {'C', {12, 8, -2, 8, 10, 0}}};
  sc::Ops poly = {{'M', {0, 0}}, {'L', {10, 0}}, {'L', {10, 10}}, {'L', {2, 3}}};
  for (double t : {GSHGEOM_STROKE_CURVES_TOLERANCE, 0x1p-8}) {
    for (int join : {sc::MITER, sc::ROUND, sc::BEVEL}) {
      for (const sc::Ops * ops : {&cubic, &loop, &poly}) {
        sc::Case c{"fine", *ops, sc::style(3, join, 4, sc::RCAP), sc::map_id(), 4};
        run(c, t, false, 0.125);
      }
    }
  }
  sc::Case c{"circle h below r", sc::circle(0, 0, 6), sc::style(3, sc::ROUND, 4, sc::BUTT), sc::map_id(), 4};
  run(c, GSHGEOM_STROKE_CURVES_TOLERANCE, false, 0.125);
}

TEST(StrokeCurvesOracle, TheDeviationAndTheAgreementAreRecorded) {
  RecordProperty("worst_deviation_over_tolerance", std::to_string(g_worst_dev));
  RecordProperty("worst_deviation_over_tolerance_fine", std::to_string(g_worst_fine));
  RecordProperty("worst_polyline_disagreement_over_tau", std::to_string(g_worst_poly));
  std::printf("stroke curves: worst deviation %.4f (%.4f at t / 8 resolution) of the tolerance; worst disagreement with the polyline %.4f tau "
              "(%.4f of the proved bound)\n",
      g_worst_dev, g_worst_fine, g_worst_poly, g_worst_poly_over_bound);
  std::sort(g_slow.begin(), g_slow.end());
  for (size_t i = g_slow.size() > 12 ? g_slow.size() - 12 : 0; i < g_slow.size(); i++) {
    std::printf("slow %.2f %s\n", g_slow[i].first, g_slow[i].second.c_str());
  }
  std::printf("time: expand %.1f compare %.1f cubic %.1f polyline %.1f\n", g_tm[0], g_tm[1], g_tm[2], g_tm[3]);
  EXPECT_GT(g_worst_dev, 0.0) << "nothing was measured";
  EXPECT_LE(g_worst_dev, 1.26);
  EXPECT_LE(g_worst_fine, 1.0 + 0.125 + 0.01);
  EXPECT_LE(g_worst_poly, 2.0);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
