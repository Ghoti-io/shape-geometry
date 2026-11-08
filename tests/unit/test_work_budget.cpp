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
 * The work budget (`GSHGEOM_Limits.max_work`, design.md "Work budget"): the work a call may spend on what its output count does
 * not bound.  The units a call spends are counted through the internal seams (along_internal.h, stroke_internal.h) and the
 * pieces' accessor, so the tests say exactly what the budget is made of:
 *
 * - the work of fixed inputs is pinned, so a charge that is taken out, doubled or moved shows as a changed number;
 * - the boundary is exact: a call that spends W units succeeds with a budget of W and is refused with W - 1, and the answer at W
 *   is the answer without a budget;
 * - the refusal depends on nothing but the call: not on a call made before it, a refusal before it, or the objects it was
 *   called with;
 * - the four inputs that a fuzzer found to take seconds (pass B of the shape family) are refused at the default budget, by the
 *   budget (not by the count that ended them after the work was done), and in a fraction of the time they took.
 *
 * Copyright 2026 by Corey Pennycuff
 */

#include "../test_helpers.h"

extern "C" {
#include "../../src/along/along_internal.h"
#include "../../src/stroke/stroke_internal.h"
}
#include "../alloc_tracker.h"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

double seconds_since(Clock::time_point t0) {
  return std::chrono::duration<double>(Clock::now() - t0).count();
}

GSHGEOM_Point P(double x, double y) {
  return GSHGEOM_Point{x, y, 0.0};
}

GSHGEOM_Map identity() {
  GSHGEOM_Map m;
  EXPECT_EQ(gshgeom_map_identity(2, &m), GSHGEOM_OK);
  return m;
}

GSHGEOM_Limits with_work(uint64_t work) {
  GSHGEOM_Limits l;
  gshgeom_limits_default(&l);
  l.max_work = work;
  return l;
}

/** The units of the largest charge any call here makes at once, plus slack: a refusal for the budget leaves the budget this close. */
constexpr uint64_t kLargestCharge = 4u * 8192u;

// ---------------------------------------------------------------------------
// The measured path
// ---------------------------------------------------------------------------

struct Measure {
  GSHGEOM_Result result;
  uint64_t work;
  uint64_t hash;
  size_t rows;
};

Measure measure(const GSHGEOM_Path * path, const GSHGEOM_Limits * limits, double tolerance = 0.0,
    const GSHGEOM_Allocator * alloc = nullptr) {
  Measure out = {GSHGEOM_OK, 12345, 0, 0};
  GSHGEOM_Measured * m = reinterpret_cast<GSHGEOM_Measured *>(0x1);
  out.result = gshgeom_measured_create_counted(alloc, limits, path, tolerance, &m, &out.work);
  if (out.result == GSHGEOM_OK) {
    out.hash = gshgeom_measured_hash(m);
    out.rows = gshgeom_measured_row_count(m);
    gshgeom_measured_destroy(m);
  } else {
    EXPECT_EQ(m, reinterpret_cast<GSHGEOM_Measured *>(0x1)) << "a refusal leaves the output as it was";
  }
  return out;
}

GSHGEOM_Path * path_of_cubic(GSHGEOM_Point a, GSHGEOM_Point b, GSHGEOM_Point c, GSHGEOM_Point d) {
  GSHGEOM_Path * p = nullptr;
  EXPECT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &p), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_move_to(p, a), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_cubic_to(p, b, c, d), GSHGEOM_OK);
  return p;
}

/** A closed path of @p n copies of one cubic (a loop with a tight turn), each its own curve. */
GSHGEOM_Path * path_of_many(size_t n) {
  GSHGEOM_Path * p = nullptr;
  EXPECT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &p), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_move_to(p, P(0, 0)), GSHGEOM_OK);
  for (size_t i = 0; i < n; i++) {
    EXPECT_EQ(gshgeom_path_cubic_to(p, P(900, 300), P(-700, 300), P(0, 0)), GSHGEOM_OK);
  }
  return p;
}

struct Pinned {
  const char * name;
  GSHGEOM_Path * (*make)();
  double tolerance;
  uint64_t work;
};

GSHGEOM_Path * make_s_curve() {
  return path_of_cubic(P(0, 0), P(100, 300), P(200, -300), P(300, 0));
}
GSHGEOM_Path * make_loop() {
  return path_of_cubic(P(0, 0), P(900, 300), P(-700, 300), P(0, 0));
}
GSHGEOM_Path * make_cusp() {
  return path_of_cubic(P(0, 0), P(100, 100), P(0, 100), P(100, 0));
}
GSHGEOM_Path * make_flat_cubic() {
  return path_of_cubic(P(0, 0), P(1, 0), P(2, 0), P(3, 0));
}
GSHGEOM_Path * make_quarter_circle() {
  GSHGEOM_Path * p = nullptr;
  EXPECT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &p), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_move_to(p, P(100, 0)), GSHGEOM_OK);
  const double w[3] = {1.0, std::sqrt(0.5), 1.0};
  EXPECT_EQ(gshgeom_path_rational_quadratic_to(p, P(100, 100), P(0, 100), w), GSHGEOM_OK);
  return p;
}

/**
 * The arc of the slow input of fuzz_path (pass B, 66 seconds): a half turn of radius 87 whose centre is 10^30 away, reached by a
 * line from (704, 2^-6).  At that magnitude the doubles cannot tell the curve's points from one another, so the quadrature's two
 * estimates never agree to the tolerance, whatever the refinement, and the table grew to the row limit.
 */
GSHGEOM_Path * make_slow_arc() {
  GSHGEOM_Path * p = nullptr;
  EXPECT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &p), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_move_to(p, P(0x1.60cp+9, 0x1p-6)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_line_to(p, P(0x1.f7406effcc5f6p+6, -0x1.8fa12267a88e5p+99)), GSHGEOM_OK);
  GSHGEOM_Arc a;
  std::memset(&a, 0, sizeof a);
  a.start = P(0x1.f7406effcc5f6p+6, -0x1.8fa12267a88e5p+99);
  a.end = P(-0x1.827f220067415p+5, -0x1.8fa12267a88e5p+99);
  a.center = P(0x1.3600ddff98becp+5, -0x1.8fa12267a88e5p+99);
  a.u = P(0x1.5c4p+6, 0.0);
  a.v = P(0.0, 0x1.5c4p+6);
  a.unit = GSHGEOM_ANGLE_RADIANS;
  a.t0 = 0.0;
  a.sweep = 0x1.921fb54442d18p+1;
  a.flags = GSHGEOM_ARC_CIRCULAR;
  EXPECT_EQ(gshgeom_path_append_arc(p, &a), GSHGEOM_OK);
  return p;
}

}  // namespace

// ---------------------------------------------------------------------------
// Core: the field and its default
// ---------------------------------------------------------------------------

TEST(WorkBudgetLimits, TheDefaultIsFilledInAndZeroIsTheDefaultNotUnlimited) {
  GSHGEOM_Limits l;
  std::memset(&l, 0xff, sizeof l);
  gshgeom_limits_default(&l);
  EXPECT_EQ(l.max_work, GSHGEOM_LIMITS_DEFAULT_MAX_WORK);
  EXPECT_EQ(GSHGEOM_LIMITS_DEFAULT_MAX_WORK, uint64_t(1) << 22);
  // 0 is the default whatever else the structure says, and NULL limits are the defaults: the same work, the same answer
  GSHGEOM_Path * p = make_loop();
  GSHGEOM_Limits zero;
  std::memset(&zero, 0, sizeof zero);
  GSHGEOM_Limits dflt;
  gshgeom_limits_default(&dflt);
  const Measure a = measure(p, nullptr);
  const Measure b = measure(p, &zero);
  const Measure c = measure(p, &dflt);
  ASSERT_EQ(a.result, GSHGEOM_OK);
  EXPECT_EQ(b.result, GSHGEOM_OK);
  EXPECT_EQ(c.result, GSHGEOM_OK);
  EXPECT_EQ(a.hash, b.hash);
  EXPECT_EQ(a.hash, c.hash);
  EXPECT_EQ(a.work, b.work);
  gshgeom_path_destroy(p);
}

TEST(WorkBudgetLimits, UnlimitedIsTheLargestValueAndNothingElse) {
  GSHGEOM_Path * p = make_loop();
  const GSHGEOM_Limits max = with_work(UINT64_MAX);
  const Measure m = measure(p, &max);
  EXPECT_EQ(m.result, GSHGEOM_OK);
  EXPECT_GT(m.work, 0u);
  GSHGEOM_Limits one = with_work(1);
  EXPECT_EQ(measure(p, &one).result, GSHGEOM_ERR_LIMIT) << "1 is a budget of one unit, not a flag";
  gshgeom_path_destroy(p);
}

// ---------------------------------------------------------------------------
// The measured path
// ---------------------------------------------------------------------------

// The work of fixed inputs, pinned: a subdivision costs SPLIT_WORK (256) units, the first rule and the rules of the first two halves
// cost nothing, and a curve that meets the tolerance at once costs none.
TEST(WorkBudgetMeasured, TheWorkOfFixedCurvesIsPinned) {
  const Pinned cases[] = {
      {"a cubic that meets the tolerance at once", make_flat_cubic, 0.0, 0},
      {"an S curve", make_s_curve, 0.0, 2304},
      {"a loop", make_loop, 0.0, 2560},
      {"a cusp", make_cusp, 0.0, 256},
      {"a quarter circle", make_quarter_circle, 0.0, 256},
      {"an S curve at the tightest tolerance", make_s_curve, GSHGEOM_MEASURE_TOLERANCE_MIN, 2816},
      {"a loop at the tightest tolerance", make_loop, GSHGEOM_MEASURE_TOLERANCE_MIN, 3584},
      {"a loop at the loosest tolerance", make_loop, GSHGEOM_MEASURE_TOLERANCE_MAX, 0},
  };
  const GSHGEOM_Limits max = with_work(UINT64_MAX);
  for (const Pinned & c : cases) {
    GSHGEOM_Path * p = c.make();
    const Measure m = measure(p, &max, c.tolerance);
    ASSERT_EQ(m.result, GSHGEOM_OK) << c.name;
    EXPECT_EQ(m.work % 256u, 0u) << c.name << ": whole subdivisions";
    EXPECT_EQ(m.work, c.work) << c.name;
    gshgeom_path_destroy(p);
  }
}

TEST(WorkBudgetMeasured, TheBoundaryIsExactAndTheAnswerAtItIsTheAnswerWithoutABudget) {
  const GSHGEOM_Limits max = with_work(UINT64_MAX);
  GSHGEOM_Path * (*const makers[])() = {make_s_curve, make_loop, make_cusp, make_quarter_circle};
  for (auto make : makers) {
    for (double tol : {0.0, GSHGEOM_MEASURE_TOLERANCE_MIN}) {
      GSHGEOM_Path * p = make();
      const Measure free_run = measure(p, &max, tol);
      ASSERT_EQ(free_run.result, GSHGEOM_OK);
      ASSERT_GT(free_run.work, 0u) << "the curve must need a subdivision for the boundary to mean anything";
      const GSHGEOM_Limits at = with_work(free_run.work);
      const Measure ok = measure(p, &at, tol);
      EXPECT_EQ(ok.result, GSHGEOM_OK);
      EXPECT_EQ(ok.work, free_run.work);
      EXPECT_EQ(ok.hash, free_run.hash) << "the table at the bound is the table without one";
      const GSHGEOM_Limits below = with_work(free_run.work - 1);
      const Measure refused = measure(p, &below, tol);
      EXPECT_EQ(refused.result, GSHGEOM_ERR_LIMIT);
      EXPECT_EQ(refused.work, free_run.work - 256u) << "the refused charge is the last subdivision, and what was spent before it is reported";
      // the charge is made before the work it pays for: nothing is done once the budget is gone
      EXPECT_LE(refused.work, free_run.work - 1);
      gshgeom_path_destroy(p);
    }
  }
}

TEST(WorkBudgetMeasured, AClosedPathOfManyCurvesSpendsTheSumOfTheirWork) {
  GSHGEOM_Path * one = path_of_many(1);
  const GSHGEOM_Limits max = with_work(UINT64_MAX);
  const Measure m1 = measure(one, &max);
  ASSERT_EQ(m1.result, GSHGEOM_OK);
  ASSERT_GT(m1.work, 0u);
  for (size_t n : {size_t(2), size_t(7)}) {
    GSHGEOM_Path * many = path_of_many(n);
    const Measure mn = measure(many, &max);
    EXPECT_EQ(mn.work, n * m1.work) << "the budget is the call's, over every curve";
    const GSHGEOM_Limits at = with_work(n * m1.work);
    EXPECT_EQ(measure(many, &at).result, GSHGEOM_OK);
    const GSHGEOM_Limits below = with_work(n * m1.work - 1);
    EXPECT_EQ(measure(many, &below).result, GSHGEOM_ERR_LIMIT);
    // and no one curve is allowed the whole budget: the second curve is what is refused
    const GSHGEOM_Limits first_only = with_work(m1.work);
    EXPECT_EQ(measure(many, &first_only).result, GSHGEOM_ERR_LIMIT);
    gshgeom_path_destroy(many);
  }
  gshgeom_path_destroy(one);
}

TEST(WorkBudgetMeasured, TheRefusalDependsOnTheCallAloneAndLeavesNothingBehind) {
  tracker::Tracker t;
  const GSHGEOM_Allocator alloc = tracker::allocator_for(&t);
  GSHGEOM_Path * p = make_loop();
  const GSHGEOM_Limits max = with_work(UINT64_MAX);
  const Measure ref = measure(p, &max, 0.0, &alloc);
  ASSERT_EQ(ref.result, GSHGEOM_OK);
  const GSHGEOM_Limits below = with_work(ref.work - 1);
  const GSHGEOM_Limits at = with_work(ref.work);
  // success, refusal, success, refusal in any order: the same answers, from a path and limits that are never changed
  for (int round = 0; round < 3; round++) {
    const Measure r1 = measure(p, &below, 0.0, &alloc);
    EXPECT_EQ(r1.result, GSHGEOM_ERR_LIMIT);
    EXPECT_EQ(r1.work, ref.work - 256u);
    EXPECT_TRUE(t.live.empty()) << "a refusal frees everything it allocated";
    const Measure r2 = measure(p, &at, 0.0, &alloc);
    EXPECT_EQ(r2.result, GSHGEOM_OK);
    EXPECT_EQ(r2.hash, ref.hash);
    EXPECT_EQ(r2.work, ref.work);
    EXPECT_TRUE(t.live.empty());
  }
  // a measured path already made does not charge a later one
  GSHGEOM_Measured * held = nullptr;
  ASSERT_EQ(gshgeom_measured_create(&alloc, &at, p, 0.0, &held), GSHGEOM_OK);
  EXPECT_EQ(measure(p, &below, 0.0, &alloc).result, GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(measure(p, &at, 0.0, &alloc).result, GSHGEOM_OK);
  gshgeom_measured_destroy(held);
  EXPECT_TRUE(t.live.empty());
  gshgeom_path_destroy(p);
}

TEST(WorkBudgetMeasured, ItComposesWithTheRowLimitEitherOneRefuses) {
  GSHGEOM_Path * p = make_loop();
  const GSHGEOM_Limits max = with_work(UINT64_MAX);
  const Measure ref = measure(p, &max);
  ASSERT_EQ(ref.result, GSHGEOM_OK);
  ASSERT_GT(ref.rows, 8u);
  GSHGEOM_Limits both;
  gshgeom_limits_default(&both);
  both.max_points = ref.rows;
  both.max_work = ref.work;
  EXPECT_EQ(measure(p, &both).result, GSHGEOM_OK) << "both met exactly";
  both.max_work = ref.work - 1;
  EXPECT_EQ(measure(p, &both).result, GSHGEOM_ERR_LIMIT) << "the work is what is short";
  both.max_work = ref.work;
  both.max_points = ref.rows - 1;
  EXPECT_EQ(measure(p, &both).result, GSHGEOM_ERR_LIMIT) << "the rows are what is short";
  both.max_points = 0;  // the defaults
  both.max_work = 0;
  EXPECT_EQ(measure(p, &both).result, GSHGEOM_OK);
  gshgeom_path_destroy(p);
}

TEST(WorkBudgetMeasured, ASingleCurveLengthHasTheDefaultBudgetToo) {
  // gshgeom_measured_curve_length takes no limits: its budget is the default.  A curve far below it is measured, and it is the
  // length of the same curve as a path.
  GSHGEOM_Curve curve;
  std::memset(&curve, 0, sizeof curve);
  curve.kind = GSHGEOM_CURVE_CUBIC;
  curve.dim = 2;
  curve.p[0] = P(0, 0);
  curve.p[1] = P(900, 300);
  curve.p[2] = P(-700, 300);
  curve.p[3] = P(0, 0);
  double length = 0.0;
  ASSERT_EQ(gshgeom_measured_curve_length(&curve, 0.0, &length), GSHGEOM_OK);
  GSHGEOM_Path * p = make_loop();
  GSHGEOM_Measured * m = nullptr;
  ASSERT_EQ(gshgeom_measured_create(nullptr, nullptr, p, 0.0, &m), GSHGEOM_OK);
  EXPECT_EQ(length, gshgeom_measured_length(m));
  gshgeom_measured_destroy(m);
  gshgeom_path_destroy(p);
}

TEST(WorkBudgetMeasured, ASingleCurveThatNeverMeetsTheToleranceIsRefusedAtTheDefaultToo) {
  GSHGEOM_SKIP_UNDER_VALGRIND();
  GSHGEOM_Path * p = make_slow_arc();
  GSHGEOM_Arc a;
  ASSERT_EQ(gshgeom_path_arc(p, 2, &a), GSHGEOM_OK);
  GSHGEOM_Curve pieces[GSHGEOM_ARC_MAX_PIECES];
  int n = 0;
  ASSERT_EQ(gshgeom_arc_pieces(&a, pieces, &n), GSHGEOM_OK);
  size_t refused = 0;
  for (int i = 0; i < n; i++) {
    double length = -1.0;
    const GSHGEOM_Result r = gshgeom_measured_curve_length(&pieces[i], GSHGEOM_MEASURE_TOLERANCE_MIN, &length);
    EXPECT_TRUE(r == GSHGEOM_OK || r == GSHGEOM_ERR_LIMIT);
    refused += r == GSHGEOM_ERR_LIMIT;
    if (r == GSHGEOM_ERR_LIMIT) {
      EXPECT_EQ(length, -1.0) << "a refusal writes nothing";
    }
  }
  EXPECT_GT(refused, 0u) << "at least one of the pieces of the slow arc is past the default budget";
  gshgeom_path_destroy(p);
}

// The slow unit of pass B (fuzz_path, 66 s before the budget): an arc whose refinement never meets the tolerance, emitted panels
// until the row limit of 4 Mi.  The budget refuses it at the default after its
// own 2^22 units, and in a fraction of a second.
TEST(WorkBudgetMeasured, TheSlowArcOfThePassBFuzzerIsRefusedByTheBudgetAtOnce) {
  GSHGEOM_SKIP_UNDER_VALGRIND();
  GSHGEOM_Path * p = make_slow_arc();
  const auto t0 = Clock::now();
  const Measure m = measure(p, nullptr, GSHGEOM_MEASURE_TOLERANCE_MIN);
  const double seconds = seconds_since(t0);
  EXPECT_EQ(m.result, GSHGEOM_ERR_LIMIT);
  EXPECT_LE(m.work, GSHGEOM_LIMITS_DEFAULT_MAX_WORK);
  EXPECT_GT(m.work + 256u, GSHGEOM_LIMITS_DEFAULT_MAX_WORK) << "it was the budget that ended it, not the row limit";
  // a guard, not the requirement: 66 s before, a few tenths of a second with the budget; sanitizers are several times slower
  EXPECT_LT(seconds, 20.0);
  GSHGEOM_Limits zero;
  std::memset(&zero, 0, sizeof zero);
  EXPECT_EQ(measure(p, &zero, GSHGEOM_MEASURE_TOLERANCE_MIN).result, GSHGEOM_ERR_LIMIT) << "0 is the default, not unlimited";
  gshgeom_path_destroy(p);
}

// A smaller budget refuses the same arc sooner, after exactly its own units, whatever the machine: the count is the work.
TEST(WorkBudgetMeasured, ASmallerBudgetRefusesTheSameArcAfterExactlyItsOwnUnits) {
  GSHGEOM_Path * p = make_slow_arc();
  for (uint64_t budget : {uint64_t(1) << 10, uint64_t(1) << 14, (uint64_t(1) << 16) + 130}) {
    const GSHGEOM_Limits l = with_work(budget);
    const Measure m = measure(p, &l, GSHGEOM_MEASURE_TOLERANCE_MIN);
    EXPECT_EQ(m.result, GSHGEOM_ERR_LIMIT) << budget;
    EXPECT_EQ(m.work, budget / 256u * 256u) << "whole subdivisions fill the budget to within one";
  }
  gshgeom_path_destroy(p);
}

// ---------------------------------------------------------------------------
// The stroke
// ---------------------------------------------------------------------------

namespace {

struct Expand {
  GSHGEOM_Result result;
  uint64_t work;
  uint64_t hash;
};

GSHGEOM_Stroke make_stroke(double width, uint32_t join, uint32_t cap) {
  GSHGEOM_Stroke s;
  gshgeom_stroke_init(&s);
  s.width = width;
  s.join = join;
  s.cap = cap;
  return s;
}

Expand expand(const GSHGEOM_Path * path, const GSHGEOM_Stroke & stroke, const GSHGEOM_Limits * limits,
    const GSHGEOM_Allocator * alloc = nullptr) {
  GSHGEOM_Map m = identity();
  GSHGEOM_Stroke_Outline * o = nullptr;
  GSHGEOM_Stroke_Work w;
  w.spent = 12345;
  Expand out = {GSHGEOM_OK, 0, 0};
  out.result = gshgeom_stroke_expand_counted(path, &stroke, &m, GSHGEOM_FLATTEN_TOLERANCE, alloc, limits,
      GSHGEOM_STROKE_MAX_DEPTH, &o, &w);
  out.work = w.spent;
  if (out.result == GSHGEOM_OK) {
    EXPECT_NE(o, nullptr);
    out.hash = gshgeom_stroke_outline_hash(o);
    gshgeom_stroke_outline_destroy(o);
  } else {
    EXPECT_EQ(o, nullptr) << "nothing is returned on a refusal";
  }
  return out;
}

GSHGEOM_Path * new_path() {
  GSHGEOM_Path * p = nullptr;
  EXPECT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &p), GSHGEOM_OK);
  return p;
}

GSHGEOM_Path * line_path(double length) {
  GSHGEOM_Path * p = new_path();
  EXPECT_EQ(gshgeom_path_move_to(p, P(10, 10)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_line_to(p, P(10 + length, 10)), GSHGEOM_OK);
  return p;
}

/** A polyline of @p n segments that turns every time, from a fixed recurrence. */
GSHGEOM_Path * zigzag_path(size_t n) {
  GSHGEOM_Path * p = new_path();
  EXPECT_EQ(gshgeom_path_move_to(p, P(10, 10)), GSHGEOM_OK);
  uint32_t state = 12345u;
  double x = 10.0, y = 10.0;
  for (size_t i = 0; i < n; i++) {
    state = state * 1664525u + 1013904223u;
    const double dx = 4.0 + double(state >> 24) / 16.0;
    state = state * 1664525u + 1013904223u;
    const double dy = (double(state >> 24) - 128.0) / 8.0;
    x += dx;
    y += dy;
    EXPECT_EQ(gshgeom_path_line_to(p, P(x, y)), GSHGEOM_OK);
  }
  return p;
}

GSHGEOM_Path * curve_path() {
  GSHGEOM_Path * p = new_path();
  EXPECT_EQ(gshgeom_path_move_to(p, P(0, 0)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_cubic_to(p, P(100, 300), P(200, -300), P(300, 0)), GSHGEOM_OK);
  return p;
}

struct StrokeCase {
  const char * name;
  GSHGEOM_Path * path;
  GSHGEOM_Stroke stroke;
  std::vector<double> dashes;
  uint64_t work;
};

std::vector<StrokeCase> stroke_cases() {
  std::vector<StrokeCase> v;
  v.push_back({"one segment, butt and miter", line_path(100), make_stroke(4, GSHGEOM_STROKE_JOIN_MITER, GSHGEOM_STROKE_CAP_BUTT), {}, 16});
  v.push_back({"one segment with round caps", line_path(100), make_stroke(40, GSHGEOM_STROKE_JOIN_ROUND, GSHGEOM_STROKE_CAP_ROUND), {}, 676});
  v.push_back({"a wide one segment with round caps", line_path(100), make_stroke(400, GSHGEOM_STROKE_JOIN_ROUND, GSHGEOM_STROKE_CAP_ROUND), {}, 2596});
  v.push_back({"a polyline of round joins", zigzag_path(60), make_stroke(6, GSHGEOM_STROKE_JOIN_ROUND, GSHGEOM_STROKE_CAP_ROUND), {}, 4230});
  v.push_back({"a polyline of miter joins", zigzag_path(60), make_stroke(6, GSHGEOM_STROKE_JOIN_MITER, GSHGEOM_STROKE_CAP_SQUARE), {}, 1936});
  v.push_back({"a cubic", curve_path(), make_stroke(10, GSHGEOM_STROKE_JOIN_ROUND, GSHGEOM_STROKE_CAP_ROUND), {}, 6540});
  v.push_back({"round dots along a line", line_path(300), make_stroke(4, GSHGEOM_STROKE_JOIN_ROUND, GSHGEOM_STROKE_CAP_ROUND), {0.0, 3.0}, 103756});
  return v;
}

}  // namespace

// The work of fixed strokes, pinned: a regression of what each charge is and where it is made.
TEST(WorkBudgetStroke, TheWorkOfFixedStrokesIsPinned) {
  std::vector<StrokeCase> cases = stroke_cases();
  const GSHGEOM_Limits max = with_work(UINT64_MAX);
  for (StrokeCase & c : cases) {
    if (!c.dashes.empty()) {
      c.stroke.dash_count = c.dashes.size();
      c.stroke.dashes = c.dashes.data();
    }
    const Expand e = expand(c.path, c.stroke, &max);
    ASSERT_EQ(e.result, GSHGEOM_OK) << c.name;
    EXPECT_EQ(e.work, c.work) << c.name;
    gshgeom_path_destroy(c.path);
  }
}

TEST(WorkBudgetStroke, TheBoundaryIsExactAndTheAnswerAtItIsTheAnswerWithoutABudget) {
  std::vector<StrokeCase> cases = stroke_cases();
  const GSHGEOM_Limits max = with_work(UINT64_MAX);
  for (StrokeCase & c : cases) {
    if (!c.dashes.empty()) {
      c.stroke.dash_count = c.dashes.size();
      c.stroke.dashes = c.dashes.data();
    }
    const Expand free_run = expand(c.path, c.stroke, &max);
    ASSERT_EQ(free_run.result, GSHGEOM_OK) << c.name;
    ASSERT_GT(free_run.work, 0u) << c.name << ": it must be charged for the boundary to mean anything";
    const GSHGEOM_Limits at = with_work(free_run.work);
    const Expand ok = expand(c.path, c.stroke, &at);
    EXPECT_EQ(ok.result, GSHGEOM_OK) << c.name;
    EXPECT_EQ(ok.work, free_run.work) << c.name;
    EXPECT_EQ(ok.hash, free_run.hash) << c.name << ": the outline at the bound is the outline without one";
    const GSHGEOM_Limits below = with_work(free_run.work - 1);
    const Expand refused = expand(c.path, c.stroke, &below);
    EXPECT_EQ(refused.result, GSHGEOM_ERR_LIMIT) << c.name;
    EXPECT_LT(refused.work, free_run.work) << c.name;
    EXPECT_LE(refused.work, free_run.work - 1) << c.name;
    EXPECT_GE(refused.work + kLargestCharge, free_run.work - 1) << c.name << ": it stopped at the charge that passed the budget";
    gshgeom_path_destroy(c.path);
  }
}

// A round cap is an approximating piece, so the first, conforming pass finds an inflated piece and the stroke is drawn again with
// every piece inflated: the budget is the call's, over both passes.
TEST(WorkBudgetStroke, ADrawnTwiceStrokeSpendsOneBudgetOverBothPasses) {
  GSHGEOM_Path * p = line_path(100);
  const GSHGEOM_Stroke round = make_stroke(40, GSHGEOM_STROKE_JOIN_ROUND, GSHGEOM_STROKE_CAP_ROUND);
  const GSHGEOM_Limits max = with_work(UINT64_MAX);
  const Expand both = expand(p, round, &max);
  ASSERT_EQ(both.result, GSHGEOM_OK);
  // the same stroke with butt caps is one conforming pass: its work is the segment's alone, and the round caps add the rest
  const GSHGEOM_Stroke butt = make_stroke(40, GSHGEOM_STROKE_JOIN_ROUND, GSHGEOM_STROKE_CAP_BUTT);
  const Expand once = expand(p, butt, &max);
  ASSERT_EQ(once.result, GSHGEOM_OK);
  EXPECT_GT(both.work, 2u * once.work);
  const GSHGEOM_Limits below = with_work(both.work - 1);
  EXPECT_EQ(expand(p, round, &below).result, GSHGEOM_ERR_LIMIT);
  // a budget that holds the first pass and not the second is a refusal: the second pass does not start from zero
  const GSHGEOM_Limits first_pass_only = with_work(both.work / 2);
  EXPECT_EQ(expand(p, round, &first_pass_only).result, GSHGEOM_ERR_LIMIT);
  gshgeom_path_destroy(p);
}

TEST(WorkBudgetStroke, TheRefusalDependsOnTheCallAloneAndLeavesNothingBehind) {
  tracker::Tracker t;
  const GSHGEOM_Allocator alloc = tracker::allocator_for(&t);
  GSHGEOM_Path * p = zigzag_path(60);
  const GSHGEOM_Stroke s = make_stroke(6, GSHGEOM_STROKE_JOIN_ROUND, GSHGEOM_STROKE_CAP_ROUND);
  const GSHGEOM_Limits max = with_work(UINT64_MAX);
  const Expand ref = expand(p, s, &max, &alloc);
  ASSERT_EQ(ref.result, GSHGEOM_OK);
  ASSERT_GT(ref.work, 1u);
  const GSHGEOM_Limits below = with_work(ref.work - 1);
  const GSHGEOM_Limits at = with_work(ref.work);
  for (int round = 0; round < 3; round++) {
    const Expand r1 = expand(p, s, &below, &alloc);
    EXPECT_EQ(r1.result, GSHGEOM_ERR_LIMIT);
    EXPECT_TRUE(t.live.empty()) << "a refusal frees everything it allocated";
    const Expand r2 = expand(p, s, &at, &alloc);
    EXPECT_EQ(r2.result, GSHGEOM_OK);
    EXPECT_EQ(r2.hash, ref.hash);
    EXPECT_EQ(r2.work, ref.work);
    EXPECT_TRUE(t.live.empty());
  }
  // the same stroke under another map is another stroke: the work is the ring's and the map's, never a leftover of the last call
  GSHGEOM_Map scaled;
  ASSERT_EQ(gshgeom_map_scale(2, 2.0, 2.0, 1.0, &scaled), GSHGEOM_OK);
  GSHGEOM_Stroke_Outline * o = nullptr;
  GSHGEOM_Stroke_Work w;
  ASSERT_EQ(gshgeom_stroke_expand_counted(p, &s, &scaled, GSHGEOM_FLATTEN_TOLERANCE, &alloc, &max, GSHGEOM_STROKE_MAX_DEPTH, &o, &w), GSHGEOM_OK);
  gshgeom_stroke_outline_destroy(o);
  EXPECT_EQ(expand(p, s, &at, &alloc).work, ref.work) << "unchanged by the call between";
  EXPECT_TRUE(t.live.empty());
  gshgeom_path_destroy(p);
}

// The expansion of a hatch's dots reports its work from zero, whatever the structure held, and refuses at the same boundary.
TEST(WorkBudgetStroke, TheDotsExpansionCountsFromZeroAndItsBoundaryIsExact) {
  GSHGEOM_Path * p = new_path();
  ASSERT_EQ(gshgeom_path_move_to(p, P(20, 20)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(p, P(20, 20)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_move_to(p, P(40, 20)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(p, P(40, 20)), GSHGEOM_OK);
  const double dirs[4] = {1.0, 0.0, 0.0, 1.0};
  const GSHGEOM_Stroke s = make_stroke(6, GSHGEOM_STROKE_JOIN_ROUND, GSHGEOM_STROKE_CAP_ROUND);
  GSHGEOM_Map m = identity();
  auto run = [&](uint64_t budget, uint64_t start, uint64_t * spent) {
    const GSHGEOM_Limits l = with_work(budget);
    GSHGEOM_Stroke_Outline * o = nullptr;
    GSHGEOM_Stroke_Work w;
    w.spent = start;
    const GSHGEOM_Result r = gshgeom_stroke_expand_dots_counted(p, &s, &m, GSHGEOM_FLATTEN_TOLERANCE, nullptr, &l, dirs, 2, &o, &w);
    gshgeom_stroke_outline_destroy(o);
    *spent = w.spent;
    return r;
  };
  uint64_t a = 0, b = 0;
  ASSERT_EQ(run(UINT64_MAX, 0, &a), GSHGEOM_OK);
  ASSERT_GT(a, 0u);
  ASSERT_EQ(run(UINT64_MAX, 987654321, &b), GSHGEOM_OK);
  EXPECT_EQ(a, b) << "what the structure held before is not spent";
  EXPECT_EQ(run(a, 5, &b), GSHGEOM_OK);
  EXPECT_EQ(run(a - 1, 5, &b), GSHGEOM_ERR_LIMIT);
  gshgeom_path_destroy(p);
}

TEST(WorkBudgetStroke, ItComposesWithThePointLimitEitherOneRefuses) {
  GSHGEOM_Path * p = zigzag_path(60);
  const GSHGEOM_Stroke s = make_stroke(6, GSHGEOM_STROKE_JOIN_ROUND, GSHGEOM_STROKE_CAP_ROUND);
  const GSHGEOM_Limits max = with_work(UINT64_MAX);
  const Expand ref = expand(p, s, &max);
  ASSERT_EQ(ref.result, GSHGEOM_OK);
  GSHGEOM_Limits l;
  gshgeom_limits_default(&l);
  l.max_work = ref.work;
  EXPECT_EQ(expand(p, s, &l).result, GSHGEOM_OK);
  l.max_points = 8;
  EXPECT_EQ(expand(p, s, &l).result, GSHGEOM_ERR_LIMIT) << "the points are what is short";
  l.max_points = 0;
  l.max_work = ref.work - 1;
  EXPECT_EQ(expand(p, s, &l).result, GSHGEOM_ERR_LIMIT) << "the work is what is short";
  l.max_work = 0;
  EXPECT_EQ(expand(p, s, &l).result, GSHGEOM_OK) << "0 is the default";
  gshgeom_path_destroy(p);
}

TEST(WorkBudgetStroke, TheCurveOutlinesAreChargedToo) {
  GSHGEOM_Path * p = line_path(100);
  const GSHGEOM_Stroke s = make_stroke(40, GSHGEOM_STROKE_JOIN_ROUND, GSHGEOM_STROKE_CAP_ROUND);
  GSHGEOM_Map m = identity();
  GSHGEOM_Stroke_Curves * c = nullptr;
  const GSHGEOM_Limits max = with_work(UINT64_MAX);
  ASSERT_EQ(gshgeom_stroke_curves_expand(p, &s, &m, 0.25, GSHGEOM_STROKE_CURVES_RATIONAL, nullptr, &max, &c), GSHGEOM_OK);
  gshgeom_stroke_curves_destroy(c);
  c = nullptr;
  const GSHGEOM_Limits one = with_work(1);
  EXPECT_EQ(gshgeom_stroke_curves_expand(p, &s, &m, 0.25, GSHGEOM_STROKE_CURVES_RATIONAL, nullptr, &one, &c), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(c, nullptr);
  gshgeom_path_destroy(p);
}

// ---------------------------------------------------------------------------
// The hatch
// ---------------------------------------------------------------------------

#include "../../src/hatch/hatch_internal.h"

namespace {

GSHGEOM_Hatch one_family(GSHGEOM_Hatch_Family & f, double angle, double across, const std::vector<double> & dashes) {
  std::memset(&f, 0, sizeof f);
  f.angle = angle;
  f.offset_across = across;
  f.dash_count = dashes.size();
  f.dashes = dashes.empty() ? nullptr : dashes.data();
  GSHGEOM_Hatch h;
  std::memset(&h, 0, sizeof h);
  h.map = identity();
  h.angle_unit = GSHGEOM_ANGLE_DEGREES;
  h.family_count = 1;
  h.families = &f;
  return h;
}

/** A regular polygon of @p n points about (500, 500), closed: the region whose every point every line of a hatch reads. */
GSHGEOM_Path * ring_path(size_t n) {
  GSHGEOM_Path * p = new_path();
  for (size_t i = 0; i < n; i++) {
    const double a = 6.283185307179586 * double(i) / double(n);
    const GSHGEOM_Point q = P(500.0 + 400.0 * std::cos(a), 500.0 + 400.0 * std::sin(a));
    EXPECT_EQ(i == 0 ? gshgeom_path_move_to(p, q) : gshgeom_path_line_to(p, q), GSHGEOM_OK);
  }
  EXPECT_EQ(gshgeom_path_close(p), GSHGEOM_OK);
  return p;
}

struct Pieces {
  GSHGEOM_Result result;
  uint64_t work;
  uint64_t hash;
  uint64_t lines;
};

Pieces pieces_of(const GSHGEOM_Path * region, const GSHGEOM_Hatch & h, const GSHGEOM_Limits * limits,
    const GSHGEOM_Allocator * alloc = nullptr) {
  GSHGEOM_Map chain = identity();
  GSHGEOM_Hatch_Pieces * out = nullptr;
  Pieces r = {GSHGEOM_OK, 0, 0, 0};
  r.result = gshgeom_hatch_pieces(region, GSHGEOM_FILL_NONZERO, &h, &chain, 0, alloc, limits, &out);
  if (r.result == GSHGEOM_OK) {
    r.work = gshgeom_hatch_pieces_work(out);
    r.hash = gshgeom_hatch_pieces_hash(out);
    r.lines = gshgeom_hatch_pieces_lines(out);
    gshgeom_hatch_pieces_destroy(out);
  } else {
    EXPECT_EQ(out, nullptr);
  }
  return r;
}

/** The least budget under which @p fn succeeds, by doubling and bisection (the budget is monotone: more never refuses what less passes). */
template <class Fn>
uint64_t least_budget(Fn fn) {
  uint64_t hi = 1;
  while (fn(hi) != GSHGEOM_OK) {
    EXPECT_LT(hi, uint64_t(1) << 40) << "no budget passes";
    if (hi >= (uint64_t(1) << 40)) {
      return 0;
    }
    hi *= 2;
  }
  uint64_t lo = hi / 2 + 1;
  if (hi == 1) {
    return 1;
  }
  while (lo < hi) {
    const uint64_t mid = lo + (hi - lo) / 2;
    if (fn(mid) == GSHGEOM_OK) {
      hi = mid;
    } else {
      lo = mid + 1;
    }
  }
  return lo;
}

}  // namespace

// Each line reads every point of the region: a unit for each sixteen points after the first 256 units, for each line.
TEST(WorkBudgetHatch, TheLinesOfARegionOfManyPointsSpendWhatTheyReadAndTheBoundaryIsExact) {
  GSHGEOM_Hatch_Family f;
  const std::vector<double> none;
  const GSHGEOM_Hatch h = one_family(f, 0.0, 100.0, none);
  for (size_t points : {size_t(2000), size_t(20000)}) {
    GSHGEOM_Path * region = ring_path(points);
    const GSHGEOM_Limits max = with_work(UINT64_MAX);
    const Pieces ref = pieces_of(region, h, &max);
    ASSERT_EQ(ref.result, GSHGEOM_OK);
    ASSERT_GT(ref.lines, 5u);
    EXPECT_EQ(ref.work, points >= 4352 ? ref.lines * (points / 16 - 256) : 0u) << points << " points, " << ref.lines << " lines";
    if (ref.work > 0) {
      const GSHGEOM_Limits at = with_work(ref.work);
      const Pieces ok = pieces_of(region, h, &at);
      EXPECT_EQ(ok.result, GSHGEOM_OK);
      EXPECT_EQ(ok.hash, ref.hash);
      const GSHGEOM_Limits below = with_work(ref.work - 1);
      EXPECT_EQ(pieces_of(region, h, &below).result, GSHGEOM_ERR_LIMIT);
    }
    gshgeom_path_destroy(region);
  }
}

// The first 256 units of a line are free: a region of 4096 points costs a line nothing (256 units), one of 4112 costs it one unit.
TEST(WorkBudgetHatch, TheFirstTwoHundredAndFiftySixUnitsOfALineAreFree) {
  GSHGEOM_Hatch_Family f;
  const std::vector<double> none;
  const GSHGEOM_Hatch h = one_family(f, 0.0, 100.0, none);
  const GSHGEOM_Limits max = with_work(UINT64_MAX);
  for (size_t points : {size_t(4096), size_t(4111), size_t(4112), size_t(4127), size_t(4128)}) {
    GSHGEOM_Path * region = ring_path(points);
    const Pieces ref = pieces_of(region, h, &max);
    ASSERT_EQ(ref.result, GSHGEOM_OK);
    EXPECT_EQ(ref.work, ref.lines * (points / 16 > 256 ? points / 16 - 256 : 0)) << points;
    gshgeom_path_destroy(region);
  }
}

// The steps of a dash pattern across a line, from the lengths, three units each after the first 256 units, for each line.
TEST(WorkBudgetHatch, ADashPatternSpendsItsStepsAcrossTheLinesAndTheBoundaryIsExact) {
  GSHGEOM_Hatch_Family f;
  const std::vector<double> dashes = {0.5, -0.5};  // the family points into it
  const GSHGEOM_Hatch h = one_family(f, 0.0, 250.0, dashes);
  GSHGEOM_Path * region = ring_path(8);
  const GSHGEOM_Limits max = with_work(UINT64_MAX);
  const Pieces ref = pieces_of(region, h, &max);
  ASSERT_EQ(ref.result, GSHGEOM_OK);
  EXPECT_EQ(ref.work, 11196u) << "3 units a step of the walk, after the first 256 units of each line";
  ASSERT_GT(ref.work, 0u);
  const GSHGEOM_Limits at = with_work(ref.work);
  EXPECT_EQ(pieces_of(region, h, &at).hash, ref.hash);
  const GSHGEOM_Limits below = with_work(ref.work - 1);
  EXPECT_EQ(pieces_of(region, h, &below).result, GSHGEOM_ERR_LIMIT);
  gshgeom_path_destroy(region);
}

// The two stages of the combined entry (the pieces, then their stroke) spend one budget: the least that passes is the sum of the
// least of each stage.
TEST(WorkBudgetHatch, TheCombinedEntrySpendsOneBudgetOverBothStages) {
  GSHGEOM_Hatch_Family f;
  const std::vector<double> dashes = {0.0, 6.0};
  const GSHGEOM_Hatch h = one_family(f, 0.0, 40.0, dashes);
  GSHGEOM_Path * region = ring_path(8000);
  GSHGEOM_Map chain = identity();
  const GSHGEOM_Stroke st = make_stroke(3, GSHGEOM_STROKE_JOIN_ROUND, GSHGEOM_STROKE_CAP_ROUND);
  const GSHGEOM_Limits max = with_work(UINT64_MAX);
  const Pieces made = pieces_of(region, h, &max);
  ASSERT_EQ(made.result, GSHGEOM_OK);
  GSHGEOM_Hatch_Pieces * pieces = nullptr;
  ASSERT_EQ(gshgeom_hatch_pieces(region, GSHGEOM_FILL_NONZERO, &h, &chain, 0, nullptr, &max, &pieces), GSHGEOM_OK);
  auto stage = [&](uint64_t budget) {
    const GSHGEOM_Limits l = with_work(budget);
    GSHGEOM_Stroke_Outline * o = nullptr;
    const GSHGEOM_Result r = gshgeom_hatch_expand_pieces(pieces, &st, &chain, GSHGEOM_FLATTEN_TOLERANCE, nullptr, &l, &o);
    gshgeom_stroke_outline_destroy(o);
    return r;
  };
  auto combined = [&](uint64_t budget) {
    const GSHGEOM_Limits l = with_work(budget);
    GSHGEOM_Stroke_Outline * o = nullptr;
    const GSHGEOM_Result r = gshgeom_hatch_expand(region, GSHGEOM_FILL_NONZERO, &h, &st, &chain, 0, nullptr, &l, &o);
    gshgeom_stroke_outline_destroy(o);
    return r;
  };
  const uint64_t expand_least = least_budget(stage);
  const uint64_t combined_least = least_budget(combined);
  EXPECT_GT(expand_least, 1u);
  EXPECT_GT(made.work, 0u) << "the first stage must spend for the second to be handed less";
  EXPECT_EQ(combined_least, made.work + expand_least);
  EXPECT_EQ(combined(made.work), GSHGEOM_ERR_LIMIT) << "a budget the first stage spends whole leaves the second nothing, which is not the default";
  EXPECT_EQ(combined(made.work + expand_least - 1), GSHGEOM_ERR_LIMIT);
  gshgeom_hatch_pieces_destroy(pieces);
  gshgeom_path_destroy(region);
}

namespace {

/** The three slow units of fuzz_hatch (pass B), as the harness decoded them. */
struct SlowHatch {
  const char * name;
  GSHGEOM_Hatch_Family family[3];
  size_t families;
  std::vector<double> dashes[3];
  int nesting;
  GSHGEOM_Fill_Rule rule;
  double width;
  bool omit;
  std::vector<std::vector<GSHGEOM_Point>> contours;
  std::vector<bool> closed;
  GSHGEOM_Map hmap;
};

GSHGEOM_Path * region_of(const SlowHatch & s) {
  GSHGEOM_Path * p = new_path();
  for (size_t c = 0; c < s.contours.size(); c++) {
    for (size_t i = 0; i < s.contours[c].size(); i++) {
      EXPECT_EQ(i == 0 ? gshgeom_path_move_to(p, s.contours[c][i]) : gshgeom_path_line_to(p, s.contours[c][i]), GSHGEOM_OK);
    }
    if (s.closed[c]) {
      EXPECT_EQ(gshgeom_path_close(p), GSHGEOM_OK);
    }
  }
  return p;
}

SlowHatch slow_dots() {
  // 0aabcde2: a dot hatch (a zero-length dash, a gap of 2) with spacing 1.5 over a 1000 by 960 rectangle, round dots 0.5 wide
  SlowHatch s{};
  s.name = "dots at 1.5 over 1000 by 960";
  s.families = 1;
  s.dashes[0] = {0.0, -2.0};
  s.family[0].angle = 0.0;
  s.family[0].offset_across = 1.5;
  s.nesting = 0;
  s.rule = GSHGEOM_FILL_NONZERO;
  s.width = 0.5;
  s.omit = true;
  s.contours = {{P(0x1p-1, 0x1p-2), P(0x1.f03p+9, 0x1p-2), P(0x1.f03p+9, 0x1.e02p+9), P(0x1p-1, 0x1.e02p+9)}};
  s.closed = {true};
  s.hmap = identity();
  return s;
}

SlowHatch slow_wide_dashes() {
  // 0ac26f78: lines 0.5 apart at 45 degrees with a dash of three entries, over a 2192 pixel square, stroked 2186.5 wide with round caps
  SlowHatch s{};
  s.name = "dashed lines 0.5 apart, stroked 2186.5 wide";
  s.families = 1;
  s.dashes[0] = {-0x1p-3, -0x1p-3, 0x1.45p+5};
  s.family[0].angle = 0x1.68p+5;
  s.family[0].offset_across = 0x1p-1;
  s.nesting = 0;
  s.rule = GSHGEOM_FILL_NONZERO;
  s.width = 0x1.1514p+11;
  s.contours = {{P(0x1.1514p+11, 0x1.1514p+11), P(0x1.1514p+12, 0x1.1514p+11), P(0x1.1514p+12, 0x1.1514p+12), P(0x1.1514p+11, 0x1.1514p+12)}};
  s.closed = {true};
  s.hmap = identity();
  return s;
}

SlowHatch slow_two_families() {
  // 41e010f4: two solid families 2 and 0.25 apart at 45 degrees (a third with no spacing), stroked 224.5 wide with round caps, the
  // hatch map a quarter turn, over three contours of an even-odd region
  SlowHatch s{};
  s.name = "solid lines 0.25 apart, stroked 224.5 wide";
  s.families = 3;
  s.family[0].angle = 0x1.68p+5;
  s.family[0].origin_x = -0x1p-1;
  s.family[0].origin_y = 0x1p-1;
  s.family[0].offset_across = -0x1p+1;
  s.family[1].angle = 0x1.68p+5;
  s.family[1].origin_x = -0x1p+1;
  s.family[1].origin_y = 0x1p-1;
  s.family[1].offset_across = 0x1p-2;
  s.family[2].angle = 0x1.1p+4;
  s.nesting = 1;
  s.rule = GSHGEOM_FILL_EVEN_ODD;
  s.width = 0x1.c04p+7;
  s.contours = {{P(-0x1.8p+1, 0.0), P(0x1.43acp+11, 0.0), P(0x1.434cp+10, 0x1.3p+2)},
      {P(0x1.8p+1, 0x1.2p+3), P(0.0, 0x1.2p+3), P(0.0, 0x1.2p+4), P(0x1.8p+1, 0x1.2p+4)},
      {P(0.0, 0x1.1e00p+5), P(0.0, 0x1.1e00p+5), P(0.0, 0x1.00e8p+10), P(0.0, 0x1.00e8p+10)}};
  s.closed = {false, true, true};
  const double lin[4] = {0.0, -1.0, 1.0, 0.0};  // rows (0, -1), (1, 0): a quarter turn
  const double t[2] = {0.0, 0.0};
  EXPECT_EQ(gshgeom_map_from_matrix(2, lin, t, &s.hmap), GSHGEOM_OK);
  return s;
}

/** Run one slow unit through the pieces and the stroke at the default budget; the result of the stage that ended it and its seconds. */
GSHGEOM_Result run_slow_hatch(SlowHatch & s, double * seconds, const GSHGEOM_Limits * limits, uint64_t * pieces_work) {
  GSHGEOM_Hatch h;
  std::memset(&h, 0, sizeof h);
  h.map = s.hmap;
  h.angle_unit = GSHGEOM_ANGLE_DEGREES;
  h.nesting = s.nesting;
  for (size_t i = 0; i < s.families; i++) {
    s.family[i].dash_count = s.dashes[i].size();
    s.family[i].dashes = s.dashes[i].empty() ? nullptr : s.dashes[i].data();
  }
  h.family_count = s.families;
  h.families = s.family;
  GSHGEOM_Stroke st = make_stroke(s.width, GSHGEOM_STROKE_JOIN_MITER, GSHGEOM_STROKE_CAP_ROUND);
  st.degenerate = s.omit ? GSHGEOM_STROKE_DEGENERATE_OMIT : GSHGEOM_STROKE_DEGENERATE_DRAW;
  GSHGEOM_Map chain = identity();
  GSHGEOM_Path * region = region_of(s);
  const auto t0 = Clock::now();
  GSHGEOM_Hatch_Pieces * pieces = nullptr;
  GSHGEOM_Result r = gshgeom_hatch_pieces(region, s.rule, &h, &chain, 0, nullptr, limits, &pieces);
  if (r == GSHGEOM_OK) {
    *pieces_work = gshgeom_hatch_pieces_work(pieces);
    GSHGEOM_Stroke_Outline * o = nullptr;
    r = gshgeom_hatch_expand_pieces(pieces, &st, &chain, GSHGEOM_FLATTEN_TOLERANCE, nullptr, limits, &o);
    gshgeom_stroke_outline_destroy(o);
    gshgeom_hatch_pieces_destroy(pieces);
  }
  *seconds = seconds_since(t0);
  gshgeom_path_destroy(region);
  return r;
}

}  // namespace

// The three slow units of fuzz_hatch (pass B) were 1.5 to 77 seconds: millions of dots or round caps, each through the discs' rings,
// until a count ended them.  At the default budget they are refused in a fraction of a second, by the budget.
TEST(WorkBudgetHatch, TheSlowUnitsOfThePassBFuzzerAreRefusedByTheBudgetAtOnce) {
  GSHGEOM_SKIP_UNDER_VALGRIND();
  SlowHatch cases[] = {slow_dots(), slow_wide_dashes(), slow_two_families()};
  for (SlowHatch & s : cases) {
    double seconds = 0.0;
    uint64_t pieces_work = 0;
    const GSHGEOM_Result r = run_slow_hatch(s, &seconds, nullptr, &pieces_work);
    EXPECT_EQ(r, GSHGEOM_ERR_LIMIT) << s.name;
    // a guard, not the requirement: 1.5 to 4 seconds before, a few tenths with the budget; sanitizers are several times slower
    EXPECT_LT(seconds, 8.0) << s.name;
    // and a smaller budget refuses the same unit the same way, sooner: the answer is the budget's, not a count's
    const GSHGEOM_Limits small = with_work(uint64_t(1) << 16);
    double seconds_small = 0.0;
    uint64_t pw = 0;
    EXPECT_EQ(run_slow_hatch(s, &seconds_small, &small, &pw), GSHGEOM_ERR_LIMIT) << s.name;
    EXPECT_LT(seconds_small, 8.0) << s.name;
  }
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
