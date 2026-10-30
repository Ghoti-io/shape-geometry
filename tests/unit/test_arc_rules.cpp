/**
 * @file
 *
 * The arc rules one at a time (AR2 to AR7, AR15): each rule broken alone, then in pairs where the order of the checks is the rule,
 * the arguments every entry point refuses, the fields equality, copy and readback keep, and the result codes that the
 * constructors report.  The cases are the ones the planted defects of tools/check-planted.py need: a rule checked in the
 * wrong place, a bound that is inclusive, a field that is dropped.
 *
 * Copyright 2026 by Corey Pennycuff
 */

#include "../test_helpers.h"

#include "../alloc_tracker.h"
extern "C" {
#include "../../src/exact/exact.h"
}

#include <ghoti.io/math/basic.h>
#include <ghoti.io/math/elementary.h>

#include <cmath>
#include <cstdint>
#include <cstring>

namespace {

GSHGEOM_Point P(double x, double y) {
  return GSHGEOM_Point{x, y, 0.0};
}

GSHGEOM_Path * make(int dim = 2) {
  GSHGEOM_Path * p = nullptr;
  EXPECT_EQ(gshgeom_path_create(nullptr, nullptr, dim, &p), GSHGEOM_OK);
  return p;
}

/** A record whose ends are the centre form's, for any conjugate pair. */
GSHGEOM_Arc consistent(GSHGEOM_Point c, GSHGEOM_Point u, GSHGEOM_Point v, double t0, double sweep,
    GSHGEOM_Angle_Unit unit = GSHGEOM_ANGLE_DEGREES, unsigned flags = 0u) {
  GSHGEOM_Arc a;
  std::memset(&a, 0, sizeof a);
  a.center = c;
  a.u = u;
  a.v = v;
  a.unit = unit;
  a.t0 = t0;
  a.sweep = sweep;
  a.flags = flags;
  GSHGEOM_Arc_Form f = {c, u, v, unit, t0, sweep, flags};
  gshgeom_arc_point(&f, t0, &a.start);
  gshgeom_arc_point(&f, t0 + sweep, &a.end);
  return a;
}

GSHGEOM_Result try_append(const GSHGEOM_Arc & a, GSHGEOM_Point current, size_t max_elements = 0) {
  GSHGEOM_Limits lim;
  gshgeom_limits_default(&lim);
  if (max_elements) {
    lim.max_elements = max_elements;
  }
  GSHGEOM_Path * p = nullptr;
  EXPECT_EQ(gshgeom_path_create(nullptr, &lim, 2, &p), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_move_to(p, current), GSHGEOM_OK);
  const GSHGEOM_Result r = gshgeom_path_append_arc(p, &a);
  gshgeom_path_destroy(p);
  return r;
}

GSHGEOM_Result try_append(const GSHGEOM_Arc & a) {
  return try_append(a, a.start);
}

}  // namespace

// AR2: every rule broken alone, on records whose ends agree with their centre forms, so that the rule and not the agreement is what is seen.
TEST(ArcRules, EachRuleOfValidityBrokenAloneOnAConsistentRecord) {
  const GSHGEOM_Arc ok = consistent(P(0, 0), P(10, 0), P(0, 10), 0, 90);
  EXPECT_EQ(try_append(ok), GSHGEOM_OK);
  GSHGEOM_Arc a = consistent(P(0, 0), P(10, 0), P(0, 10), 0, 0);
  EXPECT_EQ(try_append(a), GSHGEOM_ERR_INVALID) << "a zero sweep whose end is its start";
  a = consistent(P(0, 0), P(10, 0), P(0, 10), 0, -360.0 * (1 + 0x1p-52));
  EXPECT_EQ(try_append(a), GSHGEOM_ERR_INVALID) << "more than a turn backwards, the ends agreeing";
  a = consistent(P(0, 0), P(10, 0), P(0, 10), 0, 360.0 * (1 + 0x1p-52));
  EXPECT_EQ(try_append(a), GSHGEOM_ERR_INVALID);
  a = consistent(P(0, 0), P(10, 0), P(0, 10), 0, -360.0);
  EXPECT_EQ(try_append(a), GSHGEOM_OK) << "a full turn backwards closes";
  // an unknown unit alone, with a coordinate past the limit: invalid first
  a = ok;
  a.unit = GSHGEOM_ANGLE_UNIT_COUNT;
  a.center.x = 0x1.8p100;
  EXPECT_EQ(try_append(a, a.start), GSHGEOM_ERR_INVALID);
  // each of the twelve numbers not finite, with a coordinate past the limit: invalid first
  for (int k = 0; k < 12; k++) {
    GSHGEOM_Arc b = ok;
    double * fields[12] = {&b.start.x, &b.start.y, &b.end.x, &b.end.y, &b.center.x, &b.center.y, &b.u.x, &b.u.y, &b.v.x, &b.v.y,
        &b.t0, &b.sweep};
    *fields[k] = k % 2 ? NAN : INFINITY;
    GSHGEOM_Arc c = b;
    if (k == 4) {
      c.v.x = 0x1.8p100;
    } else {
      c.center.x = 0x1.8p100;
    }
    EXPECT_EQ(try_append(c, ok.start), GSHGEOM_ERR_INVALID) << k;
  }
  // a flat ellipse with u past the limit, and a zero sweep with v past it: invalid first
  a = consistent(P(0, 0), P(10, 0), P(0, 10), 0, 90);
  a.u = P(0x1.8p100, 0);
  a.v = P(0x1.8p100, 0);
  EXPECT_EQ(try_append(a, a.start), GSHGEOM_ERR_INVALID);
  a = ok;
  a.sweep = 0;
  a.v.x = 0x1.8p100;
  EXPECT_EQ(try_append(a, a.start), GSHGEOM_ERR_INVALID);
  // a start that is not the current point is the last check: a limit comes first
  a = ok;
  a.t0 = 90.0 * 0x1.8p40;
  EXPECT_EQ(try_append(a, P(5, 5)), GSHGEOM_ERR_LIMIT);
  // and the element limit comes before it, and after the invalid ones
  EXPECT_EQ(try_append(ok, ok.start, 1), GSHGEOM_ERR_LIMIT);
  a = ok;
  a.sweep = NAN;
  EXPECT_EQ(try_append(a, a.start, 1), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(try_append(ok, P(5, 5), 1), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(try_append(ok, ok.start, 2), GSHGEOM_OK);
  EXPECT_EQ(try_append(ok, P(5, 5), 2), GSHGEOM_ERR_INVALID);
  // negative zero is not the current point's zero
  a = consistent(P(0, 0), P(0, 10), P(-10, 0), 0, 90);
  a.start.x = -0.0;
  GSHGEOM_Path * p = make();
  EXPECT_EQ(gshgeom_path_move_to(p, P(0.0, a.start.y)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_append_arc(p, &a), GSHGEOM_ERR_INVALID);
  gshgeom_path_destroy(p);
}

TEST(ArcRules, TheLimitsAreInclusiveAndCoverEveryNumber) {
  // start and end past the limit while the centre and the half-diameters are within it
  const double L = 0x1p100;
  GSHGEOM_Arc a = consistent(P(0.75 * L, 0), P(0.75 * L, 0), P(0, 0.75 * L), 0, 90);
  EXPECT_EQ(try_append(a, P(0, 0)), GSHGEOM_ERR_LIMIT) << "a start past the limit cannot be a current point";
  a = consistent(P(0, 0.75 * L), P(0, 0.75 * L), P(-0.75 * L, 0), 0, 90);
  EXPECT_EQ(try_append(a, P(0, 0)), GSHGEOM_ERR_LIMIT) << "the start and end are checked, not only the centre form";
  a = consistent(P(0.75 * L, 0), P(-0.75 * L, 0), P(0, 0.75 * L), 0, 90);  // start (0, 0), end (0.75 L, 0.75 L): only the end is past
  a.end = P(1.5 * L, 0);
  EXPECT_EQ(try_append(a, a.start), GSHGEOM_ERR_LIMIT);
  a = consistent(P(L, 0), P(1, 0), P(0, 1), 90, 90);  // start (L, 1): exactly at the limit, accepted
  EXPECT_EQ(try_append(a), GSHGEOM_OK);
  // t0 at exactly the quarter-turn limit in degrees is accepted, one quarter more is not, and the sign does not matter
  a = consistent(P(0, 0), P(10, 0), P(0, 10), 90.0 * GSHGEOM_MAP_MAX_QUARTERS, -90);
  EXPECT_EQ(try_append(a), GSHGEOM_OK);
  a = consistent(P(0, 0), P(10, 0), P(0, 10), -90.0 * GSHGEOM_MAP_MAX_QUARTERS, 90);
  EXPECT_EQ(try_append(a), GSHGEOM_OK);
  a = consistent(P(0, 0), P(10, 0), P(0, 10), 0, 90);
  a.t0 = -90.0 * (GSHGEOM_MAP_MAX_QUARTERS + 4);
  EXPECT_EQ(try_append(a, a.start), GSHGEOM_ERR_LIMIT);
  // in radians the limit is measured in quarter turns of radians, not of degrees
  a = consistent(P(0, 0), P(10, 0), P(0, 10), 0, 1, GSHGEOM_ANGLE_RADIANS);
  a.t0 = 1.5707963267948966 * (GSHGEOM_MAP_MAX_QUARTERS * 1.01);
  EXPECT_EQ(try_append(a, a.start), GSHGEOM_ERR_LIMIT);
  a.t0 = 1.5707963267948966 * (GSHGEOM_MAP_MAX_QUARTERS * 0.5);
  GSHGEOM_Arc_Form f = {a.center, a.u, a.v, a.unit, a.t0, a.sweep, a.flags};
  gshgeom_arc_point(&f, a.t0, &a.start);
  gshgeom_arc_point(&f, a.t0 + a.sweep, &a.end);
  EXPECT_EQ(try_append(a), GSHGEOM_OK);
}

TEST(ArcRules, ARecordWhoseDecompositionLeavesTheCoordinateRangeIsALimit) {
  // every number of the record is within the limit, but the control point where the tangents at 45 and 135 degrees meet is past it
  const double h = 0x1p99;
  GSHGEOM_Arc a = consistent(P(0, h), P(0.9 * h, 0), P(0, 0.9 * h), 45, 90);
  EXPECT_EQ(try_append(a), GSHGEOM_ERR_LIMIT);
  // the same arc at a quarter of the size is fine
  a = consistent(P(0, h / 4), P(0.9 * h / 4, 0), P(0, 0.9 * h / 4), 45, 90);
  EXPECT_EQ(try_append(a), GSHGEOM_OK);
}

TEST(ArcRules, AnEllipseOfAreaThatARoundedCrossProductCallsFlatIsAccepted) {
  const double a1 = 1.0 + 0x1p-52, d1 = 1.0 - 0x1p-53;
  GSHGEOM_Arc a = consistent(P(0, 0), P(a1, 1.0), P(1.0, d1), 0, 90);
  EXPECT_EQ(try_append(a), GSHGEOM_OK) << "u x v = 2^-53 exactly, though fl(u.x v.y) = fl(u.y v.x)";
}

TEST(ArcRules, TheAgreementBoundIsInclusiveOnEveryEndAndAxis) {
  // radius 1024: scale 2048, bound 2^-29, offsets and values exact
  const GSHGEOM_Arc base = consistent(P(0, 0), P(1024, 0), P(0, 1024), 0, 90);
  const double tol = 0x1p-29;
  for (int which = 0; which < 4; which++) {
    GSHGEOM_Arc a = base;
    double * f[4] = {&a.start.x, &a.start.y, &a.end.x, &a.end.y};
    *f[which] += tol;
    EXPECT_EQ(try_append(a, a.start), GSHGEOM_OK) << "exactly at the bound " << which;
    *f[which] += tol * 0.5;
    EXPECT_EQ(try_append(a, a.start), GSHGEOM_ERR_INVALID) << "just past it " << which;
    a = base;
    *f[which] -= tol * 1.5;
    EXPECT_EQ(try_append(a, a.start), GSHGEOM_ERR_INVALID) << "below, past it " << which;
  }
  // the scale counts the ends and the centre: a far centre loosens the bound, a far end loosens it, and u alone does not
  GSHGEOM_Arc a = consistent(P(1000, 0), P(1, 0), P(0, 1), 0, 90);  // scale 1001 (the start), not 1000
  a.end.y += 0x1p-40 * 1000.5;
  EXPECT_EQ(try_append(a, a.start), GSHGEOM_OK);
  a = consistent(P(1000, 0), P(1, 0), P(0, 1), 0, 90);
  a.end.y += 0x1p-40 * 1001.5;
  EXPECT_EQ(try_append(a, a.start), GSHGEOM_ERR_INVALID);
  a = consistent(P(0, 3), P(1, 0), P(0, 1), 0, 90);  // the centre's y at 3 beats |u| + |v| = 2
  a.end.y += 0x1p-40 * 2.5;
  EXPECT_EQ(try_append(a, a.start), GSHGEOM_OK);
  // |u| + |v| are added: a circle of radius 1 has scale 2
  a = consistent(P(0, 0), P(1, 0), P(0, 1), 0, 90);
  a.end.y += 0x1p-40 * 1.75;
  EXPECT_EQ(try_append(a, a.start), GSHGEOM_OK);
  // a full turn closes to the bit, and -0 is not +0
  a = consistent(P(0, 0), P(10, 0), P(0, 10), 90, 360);
  a.start.x = 0.0;
  a.end.x = -0.0;
  GSHGEOM_Path * p = make();
  gshgeom_path_move_to(p, a.start);
  EXPECT_EQ(gshgeom_path_append_arc(p, &a), GSHGEOM_ERR_INVALID);
  gshgeom_path_destroy(p);
}

TEST(ArcRules, EveryEntryPointRefusesNullAndNeverWritesOnFailure) {
  GSHGEOM_Arc a = consistent(P(0, 0), P(10, 0), P(0, 10), 0, 90);
  GSHGEOM_Path * p = make();
  gshgeom_path_move_to(p, a.start);
  ASSERT_EQ(gshgeom_path_append_arc(p, &a), GSHGEOM_OK);
  GSHGEOM_Arc out;
  GSHGEOM_Curve pieces[GSHGEOM_ARC_MAX_PIECES];
  int n = 0;
  GSHGEOM_Path * q = nullptr;
  GSHGEOM_Map id;
  gshgeom_map_identity(2, &id);
  GSHGEOM_Point c;
  double r;
  GSHGEOM_Verb verb;
  EXPECT_EQ(gshgeom_path_arc(p, 1, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc(nullptr, 1, &out), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc(p, 2, &out), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_verb(p, 1, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_verb(p, 2, &verb), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_arc_pieces(nullptr, pieces, &n), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_arc_pieces(&a, nullptr, &n), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_arc_pieces(&a, pieces, nullptr), GSHGEOM_ERR_INVALID);
  GSHGEOM_Arc bad = a;
  bad.unit = GSHGEOM_ANGLE_UNIT_COUNT;
  EXPECT_EQ(gshgeom_arc_pieces(&bad, pieces, &n), GSHGEOM_ERR_INVALID);
  bad = a;
  bad.sweep = NAN;
  EXPECT_EQ(gshgeom_arc_pieces(&bad, pieces, &n), GSHGEOM_ERR_INVALID);
  bad = a;
  bad.t0 = INFINITY;
  EXPECT_EQ(gshgeom_arc_pieces(&bad, pieces, &n), GSHGEOM_ERR_INVALID);
  bad = a;
  bad.sweep = 0;
  EXPECT_EQ(gshgeom_arc_pieces(&bad, pieces, &n), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_lower_arcs(nullptr, nullptr, nullptr, &q), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_lower_arcs(p, nullptr, nullptr, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_copy(nullptr, nullptr, nullptr, &q), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_copy(p, nullptr, nullptr, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(q, nullptr);
  EXPECT_EQ(gshgeom_arc_map(nullptr, &a, &out), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_arc_map(&id, nullptr, &out), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_arc_map(&id, &a, nullptr), GSHGEOM_ERR_INVALID);
  GSHGEOM_Map m3;
  gshgeom_map_identity(3, &m3);
  EXPECT_EQ(gshgeom_arc_map(&m3, &a, &out), GSHGEOM_ERR_INVALID);
  GSHGEOM_Map broken = id;
  broken.m[0][0] = NAN;
  EXPECT_EQ(gshgeom_arc_map(&broken, &a, &out), GSHGEOM_ERR_INVALID);
  GSHGEOM_Map singular;
  gshgeom_map_scale(2, 0.0, 1.0, 1.0, &singular);
  EXPECT_EQ(gshgeom_arc_map(&singular, &a, &out), GSHGEOM_ERR_INVALID);
  // a map that throws a result past the limit is a limit, and the result is not written
  GSHGEOM_Map big;
  gshgeom_map_scale(2, 0x1p99, 0x1p99, 1.0, &big);
  std::memset(&out, 0x5a, sizeof out);
  EXPECT_EQ(gshgeom_arc_map(&big, &a, &out), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(out.unit, (GSHGEOM_Angle_Unit)0x5a5a5a5a);
  EXPECT_EQ(gshgeom_arc_circle(nullptr, &c, &r, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_arc_circle(&a, &c, &r, nullptr), GSHGEOM_ERR_INVALID) << "unflagged";
  a.flags = GSHGEOM_ARC_CIRCULAR;
  EXPECT_EQ(gshgeom_arc_circle(&a, nullptr, &r, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_arc_circle(&a, &c, nullptr, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_arc_circle(&a, &c, &r, nullptr), GSHGEOM_OK);
  EXPECT_EQ(c.x, 0.0);
  EXPECT_EQ(c.y, 0.0);
  EXPECT_EQ(r, 10.0) << "the radius is |u|";
  a.v = P(0, 4);
  a.u = P(3, 4);
  a.flags = GSHGEOM_ARC_CIRCULAR;
  EXPECT_EQ(gshgeom_arc_circle(&a, &c, &r, nullptr), GSHGEOM_OK);
  EXPECT_EQ(r, 5.0) << "hypot of both components of u";
  a.u = P(3, 4);
  a.v = P(6, 8);
  EXPECT_EQ(gshgeom_arc_circle(&a, &c, &r, nullptr), GSHGEOM_ERR_INVALID) << "a flagged arc whose ellipse is flat";
  a.center = P(7, -9);
  a.v = P(-4, 3);
  EXPECT_EQ(gshgeom_arc_circle(&a, &c, &r, nullptr), GSHGEOM_OK);
  EXPECT_EQ(c.x, 7.0);
  EXPECT_EQ(c.y, -9.0);
  gshgeom_path_destroy(p);
}

TEST(ArcRules, TheSenseIsTheSweepSignTimesTheOrientation) {
  for (int sw : {1, -1}) {
    for (int mirror : {0, 1}) {
      GSHGEOM_Arc a = consistent(P(0, 0), P(10, 0), P(0, mirror ? -10 : 10), 0, 90.0 * sw, GSHGEOM_ANGLE_DEGREES, GSHGEOM_ARC_CIRCULAR);
      GSHGEOM_Point c;
      double r;
      int ccw = -1;
      ASSERT_EQ(gshgeom_arc_circle(&a, &c, &r, &ccw), GSHGEOM_OK);
      EXPECT_EQ(ccw, (sw * (mirror ? -1 : 1)) > 0 ? 1 : 0) << sw << " " << mirror;
    }
  }
}

TEST(ArcRules, EqualityCopyAndReadbackKeepEveryField) {
  auto build = [](const GSHGEOM_Arc & a) {
    GSHGEOM_Path * p = make();
    gshgeom_path_move_to(p, a.start);
    EXPECT_EQ(gshgeom_path_append_arc(p, &a), GSHGEOM_OK);
    return p;
  };
  const GSHGEOM_Arc base = consistent(P(1, 2), P(10, 0), P(0, 10), 10, 100, GSHGEOM_ANGLE_DEGREES, GSHGEOM_ARC_CIRCULAR);
  GSHGEOM_Path * p = build(base);
  GSHGEOM_Path * same = build(base);
  EXPECT_TRUE(gshgeom_path_equal(p, same));
  EXPECT_TRUE(gshgeom_path_equal(nullptr, nullptr));
  EXPECT_FALSE(gshgeom_path_equal(p, nullptr));
  EXPECT_FALSE(gshgeom_path_equal(nullptr, p));
  // the same numbers in another unit, flag or sweep are another record
  GSHGEOM_Arc a = base;
  a.flags = 0u;
  GSHGEOM_Path * q = build(a);
  EXPECT_FALSE(gshgeom_path_equal(p, q));
  gshgeom_path_destroy(q);
  a = consistent(P(1, 2), P(10, 0), P(0, 10), 10.0 / 90.0 * 100.0, 100.0 / 90.0 * 100.0, GSHGEOM_ANGLE_GRADIANS, GSHGEOM_ARC_CIRCULAR);
  a.start = base.start;
  a.end = base.end;
  q = build(a);
  EXPECT_FALSE(gshgeom_path_equal(p, q)) << "the unit is part of the record";
  gshgeom_path_destroy(q);
  a = base;
  a.sweep = 100.0 + 0x1p-40;
  a.end = consistent(P(1, 2), P(10, 0), P(0, 10), 10, a.sweep).end;
  q = build(a);
  EXPECT_FALSE(gshgeom_path_equal(p, q));
  gshgeom_path_destroy(q);
  a = base;
  a.t0 = 10.0 + 0x1p-40;
  a.start = consistent(P(1, 2), P(10, 0), P(0, 10), a.t0, 100).start;
  q = build(a);
  EXPECT_FALSE(gshgeom_path_equal(p, q));
  gshgeom_path_destroy(q);
  a = base;
  a.u.y = 0x1p-45;
  q = build(a);
  EXPECT_FALSE(gshgeom_path_equal(p, q));
  gshgeom_path_destroy(q);
  a = base;
  a.center.z = 0.0;
  a.v.y = 10.0 + 0x1p-40;
  q = build(a);
  EXPECT_FALSE(gshgeom_path_equal(p, q));
  gshgeom_path_destroy(q);
  // an arc is not a line to the same end, and a 3D path is not a 2D one
  q = make();
  gshgeom_path_move_to(q, base.start);
  gshgeom_path_line_to(q, base.end);
  EXPECT_FALSE(gshgeom_path_equal(p, q));
  gshgeom_path_destroy(q);
  q = make(3);
  GSHGEOM_Path * q2 = make(2);
  EXPECT_FALSE(gshgeom_path_equal(q, q2));
  gshgeom_path_destroy(q);
  gshgeom_path_destroy(q2);
  // -0 and +0 are different in the points and in the weights (the angles)
  a = consistent(P(0, 0), P(10, 0), P(0, 10), 0, 90);
  GSHGEOM_Arc z = a;
  z.t0 = -0.0;
  GSHGEOM_Path * p0 = build(a);
  GSHGEOM_Path * pz = build(z);
  EXPECT_FALSE(gshgeom_path_equal(p0, pz));
  EXPECT_NE(gshgeom_path_hash(p0), gshgeom_path_hash(pz));
  z = a;
  z.u.y = -0.0;
  gshgeom_path_destroy(pz);
  pz = build(z);
  EXPECT_FALSE(gshgeom_path_equal(p0, pz));
  gshgeom_path_destroy(p0);
  gshgeom_path_destroy(pz);
  // a copy keeps the current point and the subpath start: the next verbs close and draw as they would on the original
  GSHGEOM_Path * cp = nullptr;
  ASSERT_EQ(gshgeom_path_copy(p, nullptr, nullptr, &cp), GSHGEOM_OK);
  EXPECT_TRUE(gshgeom_path_equal(p, cp));
  EXPECT_EQ(gshgeom_path_close(p), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_close(cp), GSHGEOM_OK);
  EXPECT_TRUE(gshgeom_path_equal(p, cp));
  EXPECT_EQ(gshgeom_path_arc_tangent_to(cp, P(50, 50), P(80, 50), 5), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_arc_tangent_to(p, P(50, 50), P(80, 50), 5), GSHGEOM_OK);
  EXPECT_TRUE(gshgeom_path_equal(p, cp));
  GSHGEOM_Path * open = build(base);
  GSHGEOM_Path * open2 = nullptr;
  ASSERT_EQ(gshgeom_path_copy(open, nullptr, nullptr, &open2), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_line_to(open, P(0, 0)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_line_to(open2, P(0, 0)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_close(open), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_close(open2), GSHGEOM_OK);
  EXPECT_TRUE(gshgeom_path_equal(open, open2));
  GSHGEOM_Path * fresh = nullptr;
  GSHGEOM_Path * none = make();
  ASSERT_EQ(gshgeom_path_copy(none, nullptr, nullptr, &fresh), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_line_to(fresh, P(1, 1)), GSHGEOM_ERR_INVALID) << "no current point, as the original";
  gshgeom_path_destroy(fresh);
  gshgeom_path_destroy(none);
  gshgeom_path_destroy(open);
  gshgeom_path_destroy(open2);
  gshgeom_path_destroy(cp);
  gshgeom_path_destroy(same);
  gshgeom_path_destroy(p);
}

TEST(ArcRules, TheBoxOfAnArcPathIsEveryPiecesAndEveryAxis) {
  // an ellipse arc through its maximum x, minimum y and others in each quadrant, over several pieces
  struct Case { double t0, sweep; };
  for (const Case k : {Case{-60, 300}, Case{100, 250}, Case{200, 320}, Case{-170, 340}}) {
    GSHGEOM_Arc a = consistent(P(10, -20), P(30, 5), P(-5, 15), k.t0, k.sweep);
    GSHGEOM_Path * p = make();
    gshgeom_path_move_to(p, a.start);
    ASSERT_EQ(gshgeom_path_append_arc(p, &a), GSHGEOM_OK);
    GSHGEOM_Bounds b;
    ASSERT_EQ(gshgeom_path_bounds(p, &b), GSHGEOM_OK);
    double lo[2] = {1e300, 1e300}, hi[2] = {-1e300, -1e300};
    for (int i = 0; i <= 20000; i++) {
      const double t = (k.t0 + k.sweep * i / 20000.0) * M_PI / 180.0;
      const double x = 10 + 30 * std::cos(t) - 5 * std::sin(t), y = -20 + 5 * std::cos(t) + 15 * std::sin(t);
      lo[0] = std::fmin(lo[0], x);
      lo[1] = std::fmin(lo[1], y);
      hi[0] = std::fmax(hi[0], x);
      hi[1] = std::fmax(hi[1], y);
    }
    const double tol = 1e-6;
    EXPECT_NEAR(b.min.x, lo[0], tol) << k.t0;
    EXPECT_NEAR(b.min.y, lo[1], tol) << k.t0;
    EXPECT_NEAR(b.max.x, hi[0], tol) << k.t0;
    EXPECT_NEAR(b.max.y, hi[1], tol) << k.t0;
    gshgeom_path_destroy(p);
  }
}

TEST(ArcRules, TheEndpointFormRefusesEachBadArgumentAndLeavesThePathAsItWas) {
  GSHGEOM_Path * p = make();
  gshgeom_path_move_to(p, P(0, 0));
  const uint64_t h0 = gshgeom_path_hash(p);
  GSHGEOM_Arc_Kind kind = GSHGEOM_ARC_KIND_ARC;
  auto call = [&](GSHGEOM_Point e, double rx, double ry, double phi, GSHGEOM_Angle_Unit unit, int large, int sweep) {
    kind = GSHGEOM_ARC_KIND_ARC;
    return gshgeom_path_arc_endpoint_to(p, e, rx, ry, phi, unit, large, sweep, &kind);
  };
  const GSHGEOM_Angle_Unit D = GSHGEOM_ANGLE_DEGREES;
  EXPECT_EQ(call(P(NAN, 1), 5, 5, 0, D, 0, 1), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(call(P(1, NAN), 5, 5, 0, D, 0, 1), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(call(P(1, 1), NAN, 5, 0, D, 0, 1), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(call(P(1, 1), 5, INFINITY, 0, D, 0, 1), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(call(P(1, 1), 5, 5, NAN, D, 0, 1), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(call(P(1, 1), 5, 5, 0, GSHGEOM_ANGLE_UNIT_COUNT, 0, 1), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(call(P(1, 1), 5, 5, 0, D, 2, 1), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(call(P(1, 1), 5, 5, 0, D, 0, 2), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(call(P(1, 1), 5, 5, 0, D, -1, 1), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(call(P(1, 1), 5, 5, 0, D, 0, -1), GSHGEOM_ERR_INVALID);
  // invalid arguments come before the equal-ends rule and the zero-radius rule
  EXPECT_EQ(call(P(0, 0), 0, 0, NAN, D, 0, 1), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(call(P(0, 0), 5, 5, 0, D, 3, 1), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(call(P(1, 1), 0, NAN, 0, D, 0, 1), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_hash(p), h0);
  // equal ends by bits: -0 differs from +0, so the arc is a line of a radius that is not zero or an arc
  EXPECT_EQ(call(P(0.0, 0), 5, 5, 0, D, 0, 1), GSHGEOM_OK) << "bit-equal ends add nothing";
  EXPECT_EQ(kind, GSHGEOM_ARC_KIND_NOTHING);
  EXPECT_EQ(call(P(-0.0, 0), 0, 5, 0, D, 0, 1), GSHGEOM_OK);
  EXPECT_EQ(kind, GSHGEOM_ARC_KIND_LINE) << "the end differs from the current point in the sign of a zero";
  // the kind is not written on a failure
  kind = GSHGEOM_ARC_KIND_NOTHING;
  EXPECT_EQ(gshgeom_path_arc_endpoint_to(p, P(1e30, 0), 0.0, 0.0, 0.0, D, 0, 1, &kind), GSHGEOM_OK);  // a line to 1e30
  EXPECT_EQ(kind, GSHGEOM_ARC_KIND_LINE);
  // a limit is the radii or the end past the coordinate limit
  GSHGEOM_Path * q = make();
  gshgeom_path_move_to(q, P(0, 0));
  kind = GSHGEOM_ARC_KIND_NOTHING;
  EXPECT_EQ(gshgeom_path_arc_endpoint_to(q, P(1, 1), 0x1.8p100, 5, 0, D, 0, 1, &kind), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(kind, GSHGEOM_ARC_KIND_NOTHING);
  EXPECT_EQ(gshgeom_path_arc_endpoint_to(q, P(1, 1), 5, 0x1.8p100, 0, D, 0, 1, nullptr), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(gshgeom_path_arc_endpoint_to(q, P(0x1.8p100, 1), 5, 5, 0, D, 0, 1, nullptr), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(gshgeom_path_arc_endpoint_to(q, P(1, 0x1.8p100), 5, 5, 0, D, 0, 1, nullptr), GSHGEOM_ERR_LIMIT);
  // a radius so large that its rotated components overflow is refused as a limit, not stored as a non-number
  EXPECT_EQ(gshgeom_path_arc_endpoint_to(q, P(1, 1), 1e308, 1e308, 45, D, 0, 1, nullptr), GSHGEOM_ERR_LIMIT);
  // no current point, a 3D path and no path
  GSHGEOM_Path * empty = make();
  EXPECT_EQ(gshgeom_path_arc_endpoint_to(empty, P(1, 1), 5, 5, 0, D, 0, 1, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_endpoint_to(nullptr, P(1, 1), 5, 5, 0, D, 0, 1, nullptr), GSHGEOM_ERR_INVALID);
  GSHGEOM_Path * three = make(3);
  gshgeom_path_move_to(three, P(0, 0));
  EXPECT_EQ(gshgeom_path_arc_endpoint_to(three, P(1, 1), 5, 5, 0, D, 0, 1, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_endpoint_to(three, P(1, 1), 0, 5, 0, D, 0, 1, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_center_to(three, P(1, 1), 5, 5, 0, D, 0, 90, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_tangent_to(three, P(1, 1), P(2, 3), 1), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_bulge_to(three, P(1, 1), 1), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_conjugate_to(three, nullptr, nullptr), GSHGEOM_ERR_INVALID);
  GSHGEOM_Arc_Form f = {P(0, 0), P(5, 0), P(0, 5), D, 0, 90, 0u};
  EXPECT_EQ(gshgeom_path_arc_conjugate_to(three, &f, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_conjugate_to(nullptr, &f, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_conjugate_to(q, nullptr, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_center_to(nullptr, P(1, 1), 5, 5, 0, D, 0, 90, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_tangent_to(nullptr, P(1, 1), P(2, 3), 1), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_bulge_to(nullptr, P(1, 1), 1), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_tangent_to(empty, P(1, 1), P(2, 3), 1), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_bulge_to(empty, P(1, 1), 1), GSHGEOM_ERR_INVALID);
  gshgeom_path_destroy(three);
  gshgeom_path_destroy(empty);
  gshgeom_path_destroy(q);
  gshgeom_path_destroy(p);
}

TEST(ArcRules, TheCentreFormRefusesEachBadArgument) {
  GSHGEOM_Path * p = make();
  const GSHGEOM_Angle_Unit D = GSHGEOM_ANGLE_DEGREES;
  auto c = [&](double rx, double ry, double phi, GSHGEOM_Angle_Unit u) {
    return gshgeom_path_arc_center_to(p, P(0, 0), rx, ry, phi, u, 0, 90, nullptr);
  };
  EXPECT_EQ(c(NAN, 5, 0, D), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(c(5, NAN, 0, D), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(c(INFINITY, 5, 0, D), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(c(5, INFINITY, 0, D), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(c(0, 5, 0, D), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(c(5, 0, 0, D), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(c(-5, 5, 0, D), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(c(5, -5, 0, D), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(c(5, 5, NAN, D), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(c(5, 5, 0, GSHGEOM_ANGLE_UNIT_COUNT), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_element_count(p), 0u);
  // nearly equal radii are not a circle; equal ones are, and only they
  for (const double ry : {5.0, std::nextafter(5.0, 6.0)}) {
    GSHGEOM_Path * q = make();
    ASSERT_EQ(gshgeom_path_arc_center_to(q, P(0, 0), 5, ry, 0, D, 0, 90, nullptr), GSHGEOM_OK);
    GSHGEOM_Arc a;
    ASSERT_EQ(gshgeom_path_arc(q, 1, &a), GSHGEOM_OK);
    EXPECT_EQ((a.flags & GSHGEOM_ARC_CIRCULAR) != 0u, ry == 5.0);
    gshgeom_path_destroy(q);
  }
  // the end override is kept to the bit, and the arc may begin a subpath with an end of its own
  GSHGEOM_Path * q = make();
  ASSERT_EQ(gshgeom_path_arc_center_to(q, P(0, 0), 5, 5, 0, D, 0, 90, nullptr), GSHGEOM_OK);
  GSHGEOM_Path * r = make();
  const GSHGEOM_Point e = P(0, 5.0 + 0x1p-45);
  ASSERT_EQ(gshgeom_path_arc_center_to(r, P(0, 0), 5, 5, 0, D, 0, 90, &e), GSHGEOM_OK);
  GSHGEOM_Arc a;
  ASSERT_EQ(gshgeom_path_arc(r, 1, &a), GSHGEOM_OK);
  EXPECT_EQ(std::memcmp(&a.end.y, &e.y, sizeof e.y), 0);
  EXPECT_EQ(a.end.z, 0.0);
  GSHGEOM_Arc b;
  ASSERT_EQ(gshgeom_path_arc(q, 1, &b), GSHGEOM_OK);
  EXPECT_EQ(b.center.z, 0.0);
  EXPECT_EQ(b.u.z, 0.0);
  EXPECT_EQ(b.v.z, 0.0);
  // a current point is the arc's start, not the form's own
  GSHGEOM_Path * s = make();
  gshgeom_path_move_to(s, P(5, 0));
  ASSERT_EQ(gshgeom_path_arc_center_to(s, P(0, 0), 5, 5, 0, D, 0, 90, nullptr), GSHGEOM_OK);
  GSHGEOM_Path * t = make();
  gshgeom_path_move_to(t, P(5, 1e-3));
  EXPECT_EQ(gshgeom_path_arc_center_to(t, P(0, 0), 5, 5, 0, D, 0, 90, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_element_count(t), 1u);
  // a form with a current point that does not match leaves no move behind either way
  GSHGEOM_Path * u = make();
  EXPECT_EQ(gshgeom_path_arc_center_to(u, P(0, 0), 5, 5, 0, D, 0, 0, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_element_count(u), 0u) << "a bad arc leaves no move behind";
  // the element limit counts the move that begins the subpath
  GSHGEOM_Limits lim;
  gshgeom_limits_default(&lim);
  lim.max_elements = 1;
  GSHGEOM_Path * w = nullptr;
  ASSERT_EQ(gshgeom_path_create(nullptr, &lim, 2, &w), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_arc_center_to(w, P(0, 0), 5, 5, 0, D, 0, 90, nullptr), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(gshgeom_path_element_count(w), 0u);
  lim.max_elements = 2;
  GSHGEOM_Path * x = nullptr;
  ASSERT_EQ(gshgeom_path_create(nullptr, &lim, 2, &x), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_arc_center_to(x, P(0, 0), 5, 5, 0, D, 0, 90, nullptr), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_element_count(x), 2u);
  for (GSHGEOM_Path * z : {p, q, r, s, t, u, w, x}) {
    gshgeom_path_destroy(z);
  }
}

TEST(ArcRules, TangentAndBulgeFormsRefuseEachBadArgument) {
  GSHGEOM_Path * p = make();
  gshgeom_path_move_to(p, P(0, 0));
  EXPECT_EQ(gshgeom_path_arc_tangent_to(p, P(10, 0), P(10, 10), NAN), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_tangent_to(p, P(10, 0), P(10, 10), INFINITY), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_tangent_to(p, P(10, 0), P(10, 10), -1), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_tangent_to(p, P(NAN, 0), P(10, 10), 1), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_tangent_to(p, P(10, NAN), P(10, 10), 1), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_tangent_to(p, P(10, 0), P(INFINITY, 10), 1), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_tangent_to(p, P(10, 0), P(10, NAN), 1), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_bulge_to(p, P(10, 0), NAN), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_bulge_to(p, P(10, 0), INFINITY), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_bulge_to(p, P(NAN, 0), 1), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_bulge_to(p, P(0, INFINITY), 1), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_element_count(p), 1u);
  // radius zero is a line (not refused), and so is a repeated point
  EXPECT_EQ(gshgeom_path_arc_tangent_to(p, P(10, 0), P(10, 10), 0), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_element_count(p), 2u);
  GSHGEOM_Element el;
  ASSERT_EQ(gshgeom_path_element(p, 1, &el), GSHGEOM_OK);
  EXPECT_EQ(el.verb, GSHGEOM_VERB_LINE);
  EXPECT_EQ(el.curve.p[1].x, 10.0);
  EXPECT_EQ(gshgeom_path_arc_tangent_to(p, P(10, 0), P(10, 10), 3), GSHGEOM_OK) << "P0 = P1: a line, not an arc";
  EXPECT_EQ(gshgeom_path_element_count(p), 3u);
  ASSERT_EQ(gshgeom_path_element(p, 2, &el), GSHGEOM_OK);
  EXPECT_EQ(el.verb, GSHGEOM_VERB_LINE);
  EXPECT_EQ(gshgeom_path_arc_tangent_to(p, P(30, 40), P(30, 40), 3), GSHGEOM_OK) << "P1 = P2: a line to P1";
  ASSERT_EQ(gshgeom_path_element(p, 3, &el), GSHGEOM_OK);
  EXPECT_EQ(el.verb, GSHGEOM_VERB_LINE);
  EXPECT_EQ(el.curve.p[1].y, 40.0);
  EXPECT_EQ(gshgeom_path_arc_tangent_to(p, P(60, 80), P(90, 120), 3), GSHGEOM_OK) << "collinear: a line";
  ASSERT_EQ(gshgeom_path_element(p, 4, &el), GSHGEOM_OK);
  EXPECT_EQ(el.verb, GSHGEOM_VERB_LINE);
  EXPECT_EQ(el.curve.p[1].x, 60.0) << "the line goes to P1";
  // a tangent point that is the current point adds no line; otherwise a line to it first
  GSHGEOM_Path * q = make();
  gshgeom_path_move_to(q, P(0, 0));
  ASSERT_EQ(gshgeom_path_arc_tangent_to(q, P(10, 0), P(10, 10), 10), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_element_count(q), 2u) << "T1 is the current point: only the arc";
  GSHGEOM_Verb verb;
  ASSERT_EQ(gshgeom_path_verb(q, 1, &verb), GSHGEOM_OK);
  EXPECT_EQ(verb, GSHGEOM_VERB_ARC);
  GSHGEOM_Arc a;
  ASSERT_EQ(gshgeom_path_arc(q, 1, &a), GSHGEOM_OK);
  EXPECT_NE(a.flags & GSHGEOM_ARC_CIRCULAR, 0u);
  gshgeom_path_destroy(q);
  q = make();
  gshgeom_path_move_to(q, P(0, 0));
  ASSERT_EQ(gshgeom_path_arc_tangent_to(q, P(100, 0), P(100, 100), 10), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_element_count(q), 3u) << "a line to T1, then the arc";
  // the bulge: equal ends add nothing (by bits), a zero bulge is a line, the circular flag is set, the sense follows the sign
  GSHGEOM_Path * b = make();
  gshgeom_path_move_to(b, P(0, 0));
  EXPECT_EQ(gshgeom_path_arc_bulge_to(b, P(0, 0), 1), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_element_count(b), 1u);
  EXPECT_EQ(gshgeom_path_arc_bulge_to(b, P(-0.0, 0), 0), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_element_count(b), 2u) << "a zero bulge to another point is a line, the sign of a zero counting";
  ASSERT_EQ(gshgeom_path_verb(b, 1, &verb), GSHGEOM_OK);
  EXPECT_EQ(verb, GSHGEOM_VERB_LINE);
  ASSERT_EQ(gshgeom_path_arc_bulge_to(b, P(10, 0), 1), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_arc(b, 2, &a), GSHGEOM_OK);
  EXPECT_NE(a.flags & GSHGEOM_ARC_CIRCULAR, 0u);
  EXPECT_NEAR(std::fabs(a.sweep), M_PI, 1e-15) << "a bulge of 1 is a half circle";
  EXPECT_GT(a.sweep, 0.0);
  ASSERT_EQ(gshgeom_path_arc_bulge_to(b, P(20, 0), -1), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_arc(b, 3, &a), GSHGEOM_OK);
  EXPECT_LT(a.sweep, 0.0);
  gshgeom_path_destroy(b);
  gshgeom_path_destroy(q);
  gshgeom_path_destroy(p);
}

TEST(ArcRules, UnitAndAnglesAreEachPartOfTheRecordWhenOnlyTheyDiffer) {
  // A sweep of 1e-14 in degrees or radians, or a t0 of 0 or 1e-14, leaves the ends where the agreement bound cannot tell them
  // apart, so the records differ in nothing but the unit or the angle.
  auto build = [](double t0, double sweep, GSHGEOM_Angle_Unit unit) {
    GSHGEOM_Arc a;
    std::memset(&a, 0, sizeof a);
    a.start = P(1, 0);
    a.end = P(1, 1.7e-16);
    a.center = P(0, 0);
    a.u = P(1, 0);
    a.v = P(0, 1);
    a.unit = unit;
    a.t0 = t0;
    a.sweep = sweep;
    GSHGEOM_Path * p = make();
    gshgeom_path_move_to(p, a.start);
    EXPECT_EQ(gshgeom_path_append_arc(p, &a), GSHGEOM_OK);
    return p;
  };
  GSHGEOM_Path * base = build(0, 1e-14, GSHGEOM_ANGLE_DEGREES);
  GSHGEOM_Path * same = build(0, 1e-14, GSHGEOM_ANGLE_DEGREES);
  GSHGEOM_Path * unit = build(0, 1e-14, GSHGEOM_ANGLE_RADIANS);
  GSHGEOM_Path * sweep = build(0, 2e-14, GSHGEOM_ANGLE_DEGREES);
  GSHGEOM_Path * t0 = build(1e-14, 1e-14, GSHGEOM_ANGLE_DEGREES);
  EXPECT_TRUE(gshgeom_path_equal(base, same));
  EXPECT_FALSE(gshgeom_path_equal(base, unit));
  EXPECT_FALSE(gshgeom_path_equal(base, sweep));
  EXPECT_FALSE(gshgeom_path_equal(base, t0));
  EXPECT_NE(gshgeom_path_hash(base), gshgeom_path_hash(unit));
  for (GSHGEOM_Path * z : {base, same, unit, sweep, t0}) {
    gshgeom_path_destroy(z);
  }
  // a 3D path's points are compared in z, and an arc is not a rational cubic that shares its numbers
  GSHGEOM_Path * a3 = make(3);
  GSHGEOM_Path * b3 = make(3);
  gshgeom_path_move_to(a3, GSHGEOM_Point{0, 0, 0});
  gshgeom_path_move_to(b3, GSHGEOM_Point{0, 0, 0});
  gshgeom_path_line_to(a3, GSHGEOM_Point{1, 1, 1});
  gshgeom_path_line_to(b3, GSHGEOM_Point{1, 1, 2});
  EXPECT_FALSE(gshgeom_path_equal(a3, b3));
  gshgeom_path_destroy(a3);
  gshgeom_path_destroy(b3);
  GSHGEOM_Path * ar = build(0, 1e-14, GSHGEOM_ANGLE_DEGREES);
  GSHGEOM_Path * rc = make();
  const double w[4] = {1, 1, 1, 1};
  gshgeom_path_move_to(rc, P(1, 0));
  gshgeom_path_rational_cubic_to(rc, P(0, 0), P(1, 0), P(0, 1.7e-16), w);
  EXPECT_FALSE(gshgeom_path_equal(ar, rc));
  gshgeom_path_destroy(ar);
  gshgeom_path_destroy(rc);
}

TEST(ArcRules, ACopyAndALoweringKeepTheStateTheNextVerbReads) {
  const GSHGEOM_Arc a = consistent(P(10, 10), P(10, 0), P(0, 10), -90, 230);
  GSHGEOM_Path * p = make();
  gshgeom_path_move_to(p, a.start);
  ASSERT_EQ(gshgeom_path_append_arc(p, &a), GSHGEOM_OK);
  GSHGEOM_Path * cp = nullptr;
  GSHGEOM_Path * lo = nullptr;
  ASSERT_EQ(gshgeom_path_copy(p, nullptr, nullptr, &cp), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_lower_arcs(p, nullptr, nullptr, &lo), GSHGEOM_OK);
  // the next verbs read the current point and the subpath start: arcTo from the arc's end, then a close
  for (GSHGEOM_Path * q : {p, cp, lo}) {
    EXPECT_EQ(gshgeom_path_arc_tangent_to(q, P(40, 40), P(40, 80), 4), GSHGEOM_OK);
    EXPECT_EQ(gshgeom_path_close(q), GSHGEOM_OK);
    EXPECT_EQ(gshgeom_path_line_to(q, P(1, 1)), GSHGEOM_OK);
  }
  EXPECT_TRUE(gshgeom_path_equal(p, cp));
  GSHGEOM_Path * lo2 = nullptr;
  ASSERT_EQ(gshgeom_path_lower_arcs(p, nullptr, nullptr, &lo2), GSHGEOM_OK);
  GSHGEOM_Bounds b1, b2;
  gshgeom_path_bounds(lo, &b1);
  gshgeom_path_bounds(lo2, &b2);
  EXPECT_EQ(std::memcmp(&b1, &b2, sizeof b1), 0);
  // the lowering is the pieces, in order, end to end
  GSHGEOM_Curve pieces[GSHGEOM_ARC_MAX_PIECES];
  int n = 0;
  ASSERT_EQ(gshgeom_arc_pieces(&a, pieces, &n), GSHGEOM_OK);
  ASSERT_EQ(n, 3);
  GSHGEOM_Path * only = nullptr;
  {
    GSHGEOM_Path * q = make();
    gshgeom_path_move_to(q, a.start);
    ASSERT_EQ(gshgeom_path_append_arc(q, &a), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_path_lower_arcs(q, nullptr, nullptr, &only), GSHGEOM_OK);
    gshgeom_path_destroy(q);
  }
  ASSERT_EQ(gshgeom_path_element_count(only), 4u);
  for (int k = 0; k < n; k++) {
    GSHGEOM_Element el;
    ASSERT_EQ(gshgeom_path_element(only, (size_t)k + 1, &el), GSHGEOM_OK);
    EXPECT_EQ(el.verb, GSHGEOM_VERB_RATIONAL_QUADRATIC);
    for (int i = 0; i < 3; i++) {
      EXPECT_EQ(std::memcmp(&el.curve.p[i], &pieces[k].p[i], sizeof el.curve.p[i]), 0) << k << " " << i;
      EXPECT_EQ(std::memcmp(&el.curve.w[i], &pieces[k].w[i], sizeof el.curve.w[i]), 0) << k << " " << i;
    }
  }
  // a failed lowering or copy releases everything it took
  for (long fail = 0; fail < 40; fail++) {
    tracker::Tracker tr;
    tr.fail_from = fail;
    GSHGEOM_Allocator al = tracker::allocator_for(&tr);
    GSHGEOM_Path * x = nullptr;
    const GSHGEOM_Result r1 = gshgeom_path_lower_arcs(p, &al, nullptr, &x);
    if (r1 == GSHGEOM_OK) {
      gshgeom_path_destroy(x);
    } else {
      EXPECT_EQ(x, nullptr);
    }
    GSHGEOM_Path * y = nullptr;
    const GSHGEOM_Result r2 = gshgeom_path_copy(p, &al, nullptr, &y);
    if (r2 == GSHGEOM_OK) {
      gshgeom_path_destroy(y);
    } else {
      EXPECT_EQ(y, nullptr);
    }
    EXPECT_TRUE(tr.live.empty()) << fail;
  }
  for (GSHGEOM_Path * z : {p, cp, lo, lo2, only}) {
    gshgeom_path_destroy(z);
  }
}

TEST(ArcRules, TheEvaluatorAndTheConstructorsReportWhatTheyCannotDo) {
  GSHGEOM_Arc_Form f = {P(1, 2), P(3, 0), P(0, 3), GSHGEOM_ANGLE_DEGREES, 0, 90, 0u};
  GSHGEOM_Point out = P(7, 7);
  EXPECT_EQ(gshgeom_arc_point(nullptr, 0, &out), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_arc_point(&f, 0, nullptr), GSHGEOM_ERR_INVALID);
  GSHGEOM_Arc_Form bad = f;
  bad.unit = GSHGEOM_ANGLE_UNIT_COUNT;
  EXPECT_EQ(gshgeom_arc_point(&bad, 0, &out), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_arc_point(&f, 1e300, &out), GSHGEOM_ERR_LIMIT) << "an angle past the quarter-turn limit";
  EXPECT_EQ(gshgeom_arc_point(&f, NAN, &out), GSHGEOM_ERR_INVALID);
  bad = f;
  bad.center.z = 5;
  bad.u.z = 6;
  ASSERT_EQ(gshgeom_arc_point(&bad, 0, &out), GSHGEOM_OK);
  EXPECT_EQ(out.z, 0.0);
  EXPECT_EQ(out.x, 4.0);
  EXPECT_EQ(out.y, 2.0);
  // a form with a z in it is stored with none
  GSHGEOM_Path * p = make();
  ASSERT_EQ(gshgeom_path_arc_conjugate_to(p, &bad, nullptr), GSHGEOM_OK);
  GSHGEOM_Arc a;
  ASSERT_EQ(gshgeom_path_arc(p, 1, &a), GSHGEOM_OK);
  EXPECT_EQ(a.center.z, 0.0);
  EXPECT_EQ(a.u.z, 0.0);
  EXPECT_EQ(a.start.z, 0.0);
  EXPECT_EQ(a.end.z, 0.0);
  gshgeom_path_destroy(p);
  // an angle past the limit in a constructor is a limit, and leaves nothing
  p = make();
  GSHGEOM_Arc_Form far = f;
  far.t0 = 1e300;
  EXPECT_EQ(gshgeom_path_arc_conjugate_to(p, &far, nullptr), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(gshgeom_path_arc_center_to(p, P(0, 0), 5, 5, 1e300, GSHGEOM_ANGLE_DEGREES, 0, 90, nullptr), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(gshgeom_path_element_count(p), 0u);
  gshgeom_path_move_to(p, P(0, 0));
  EXPECT_EQ(gshgeom_path_arc_endpoint_to(p, P(3, 4), 5, 5, 1e300, GSHGEOM_ANGLE_DEGREES, 0, 1, nullptr), GSHGEOM_ERR_LIMIT);
  // the result kind is written only by a call that succeeded
  GSHGEOM_Arc_Kind kind = GSHGEOM_ARC_KIND_NOTHING;
  GSHGEOM_Limits lim;
  gshgeom_limits_default(&lim);
  lim.max_elements = 1;
  GSHGEOM_Path * q = nullptr;
  ASSERT_EQ(gshgeom_path_create(nullptr, &lim, 2, &q), GSHGEOM_OK);
  gshgeom_path_move_to(q, P(0, 0));
  EXPECT_EQ(gshgeom_path_arc_endpoint_to(q, P(10, 0), 0, 5, 0, GSHGEOM_ANGLE_DEGREES, 0, 1, &kind), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(kind, GSHGEOM_ARC_KIND_NOTHING);
  EXPECT_EQ(gshgeom_path_arc_endpoint_to(q, P(10, 0), 5, 5, 0, GSHGEOM_ANGLE_DEGREES, 0, 1, &kind), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(kind, GSHGEOM_ARC_KIND_NOTHING);
  gshgeom_path_destroy(q);
  gshgeom_path_destroy(p);
}

TEST(ArcRules, EveryPieceBoundaryAndControlPointIsTheEvaluatorsAndTheTangentsOwn) {
  for (const double sweep : {100.0, 200.0, 250.0, 330.0, -200.0, -330.0}) {
    for (const GSHGEOM_Angle_Unit unit : {GSHGEOM_ANGLE_DEGREES, GSHGEOM_ANGLE_GRADIANS, GSHGEOM_ANGLE_TURNS}) {
      const double s = unit == GSHGEOM_ANGLE_DEGREES ? sweep : (unit == GSHGEOM_ANGLE_GRADIANS ? sweep / 0.9 : sweep / 360.0);
      const double t0 = unit == GSHGEOM_ANGLE_DEGREES ? 17.3 : (unit == GSHGEOM_ANGLE_GRADIANS ? 17.3 / 0.9 : 17.3 / 360.0);
      const GSHGEOM_Arc a = consistent(P(0.1, 0.3), P(1.3, 0.4), P(-0.9, 2.1), t0, s, unit);
      GSHGEOM_Curve pc[GSHGEOM_ARC_MAX_PIECES];
      int n = 0;
      ASSERT_EQ(gshgeom_arc_pieces(&a, pc, &n), GSHGEOM_OK);
      GSHGEOM_Arc_Form f = {a.center, a.u, a.v, a.unit, a.t0, a.sweep, a.flags};
      for (int k = 1; k < n; k++) {
        GSHGEOM_Point e;
        ASSERT_EQ(gshgeom_arc_point(&f, a.t0 + ((double)k * a.sweep) / (double)n, &e), GSHGEOM_OK);
        EXPECT_EQ(std::memcmp(&pc[k - 1].p[2], &e, 2 * sizeof(double)), 0) << sweep << " " << k;
        EXPECT_EQ(std::memcmp(&pc[k].p[0], &e, 2 * sizeof(double)), 0);
      }
      // the control point is the start plus tan(h) times the tangent there, in the order the decomposition writes it
      const double sg = a.sweep < 0 ? -1.0 : 1.0;
      const double q = unit == GSHGEOM_ANGLE_DEGREES ? 90.0 : (unit == GSHGEOM_ANGLE_GRADIANS ? 100.0 : 0.25);
      const double per = unit == GSHGEOM_ANGLE_DEGREES ? 0x1.1df46a2529d39p-6 : (unit == GSHGEOM_ANGLE_GRADIANS ? 0x1.015bf9217271ap-6 : 0x1.921fb54442d18p+2);
      const double h = (std::fabs(a.sweep) / (2.0 * n)) * per;
      for (int k = 0; k < n; k++) {
        double c, sn;
        ASSERT_EQ(gshgeom_sincos(a.unit, k == 0 ? a.t0 : a.t0 + ((double)k * a.sweep) / (double)n, &c, &sn), GSHGEOM_OK);
        const double tax = sg * ((a.v.x * c) - (a.u.x * sn));
        const double tay = sg * ((a.v.y * c) - (a.u.y * sn));
        const double th = gmath_tan(h);
        const double cx = pc[k].p[0].x + (th * tax), cy = pc[k].p[0].y + (th * tay);
        EXPECT_EQ(std::memcmp(&pc[k].p[1].x, &cx, sizeof cx), 0) << sweep << " " << k;
        EXPECT_EQ(std::memcmp(&pc[k].p[1].y, &cy, sizeof cy), 0) << sweep << " " << k;
      }
      (void)q;
    }
  }
  // a sweep just short of a quarter turn does not take the quarter turn's weight
  const GSHGEOM_Arc a = consistent(P(0, 0), P(1, 0), P(0, 1), 0, 90.0 * (1.0 - 0x1p-31));
  GSHGEOM_Curve pc[GSHGEOM_ARC_MAX_PIECES];
  int n = 0;
  ASSERT_EQ(gshgeom_arc_pieces(&a, pc, &n), GSHGEOM_OK);
  ASSERT_EQ(n, 1);
  EXPECT_GT(pc[0].w[1], gmath_sqrt(0.5));
  EXPECT_NEAR(pc[0].w[1], gmath_sqrt(0.5), 1e-8);
  // an off-grid quarter turn keeps its control point exact: (end.x, start.y) or (start.x, end.y)
  for (const double t0 : {0.0, 90.0, 180.0, 270.0}) {
    const GSHGEOM_Arc q = consistent(P(0.1, 0.3), P(0.7, 0), P(0, 0.7), t0, 90);
    GSHGEOM_Curve qc[GSHGEOM_ARC_MAX_PIECES];
    ASSERT_EQ(gshgeom_arc_pieces(&q, qc, &n), GSHGEOM_OK);
    ASSERT_EQ(n, 1);
    const bool horizontal_first = t0 == 90.0 || t0 == 270.0;
    EXPECT_EQ(qc[0].p[1].x, horizontal_first ? q.end.x : q.start.x) << t0;
    EXPECT_EQ(qc[0].p[1].y, horizontal_first ? q.start.y : q.end.y) << t0;
  }
}

TEST(ArcRules, AnExactSumWithABadTermConvertsToZero) {
  GSHGEOM_ExactSum s;
  gshgeom_exact_sum_init(&s);
  const double f[2] = {NAN, 1.0};
  gshgeom_exact_sum_add(&s, 0, 2, f);
  EXPECT_EQ(gshgeom_exact_sum_to_double(&s), 0.0);
  EXPECT_EQ(gshgeom_exact_sum_sign(&s), 0);
}

TEST(ArcRules, ThePathUnderAMapWithAToleranceRefusesOneBelowTheMinimumWhateverIsInIt) {
  GSHGEOM_Path * p = make();
  gshgeom_path_move_to(p, P(0, 0));
  gshgeom_path_line_to(p, P(5, 5));
  GSHGEOM_Map id;
  gshgeom_map_identity(2, &id);
  GSHGEOM_Device_Path * dp = nullptr;
  EXPECT_EQ(gshgeom_device_path_from_model_with_tolerance(p, &id, 1, nullptr, nullptr, &dp), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_device_path_from_model_with_tolerance(p, &id, 0, nullptr, nullptr, &dp), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_device_path_from_model_with_tolerance(p, &id, GSHGEOM_FLATTEN_MIN_TOLERANCE, nullptr, nullptr, &dp), GSHGEOM_OK);
  gshgeom_device_path_destroy(dp);
  gshgeom_path_destroy(p);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
