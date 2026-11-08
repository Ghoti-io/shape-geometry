/**
 * @file
 *
 * Arcs (arc.h): inputs built from the text of one rule each, for the rules the exhaustive arc grid could not tell from a
 * variant without them.  A refusal that a later check repeats is told apart by the code of a second fault that only the
 * later check sees (the earlier check decides which code is reported); a state is reached that the constructors do not
 * build (a map whose kind tag disagrees with its entries, a record stored raw); a constant is read at the one input that
 * uses it; and the rules that are repeated by a later check with the same code are pinned at the input that would show a
 * difference, which is how they are known to be equivalent rather than untested.
 *
 * Copyright 2026 by Corey Pennycuff
 */

#include "../test_helpers.h"

extern "C" {
#include "../../src/exact/exact.h"
#include "../../src/path/path_internal.h"
}

#include <ghoti.io/math/basic.h>

#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>

namespace {

GSHGEOM_Point P(double x, double y) {
  return GSHGEOM_Point{x, y, 0.0};
}

uint64_t bits(double v) {
  uint64_t b;
  std::memcpy(&b, &v, sizeof b);
  return b;
}

GSHGEOM_Path * make(const GSHGEOM_Limits * limits = nullptr) {
  GSHGEOM_Path * p = nullptr;
  EXPECT_EQ(gshgeom_path_create(nullptr, limits, 2, &p), GSHGEOM_OK);
  return p;
}

// a circle record in degrees about the origin: start (r, 0), a quarter turn
GSHGEOM_Arc unit_arc(double r = 1.0) {
  GSHGEOM_Arc a;
  std::memset(&a, 0, sizeof a);
  a.start = P(r, 0);
  a.end = P(0, r);
  a.center = P(0, 0);
  a.u = P(r, 0);
  a.v = P(0, r);
  a.unit = GSHGEOM_ANGLE_DEGREES;
  a.t0 = 0;
  a.sweep = 90;
  a.flags = GSHGEOM_ARC_CIRCULAR;
  return a;
}

GSHGEOM_Result append(const GSHGEOM_Arc & a, const GSHGEOM_Limits * limits = nullptr) {
  GSHGEOM_Path * p = make(limits);
  EXPECT_EQ(gshgeom_path_move_to(p, a.start), GSHGEOM_OK);
  const GSHGEOM_Result r = gshgeom_path_append_arc(p, &a);
  gshgeom_path_destroy(p);
  return r;
}

GSHGEOM_Map identity() {
  GSHGEOM_Map m;
  EXPECT_EQ(gshgeom_map_identity(2, &m), GSHGEOM_OK);
  return m;
}

const double kBeyond = 0x1p101;  // past GSHGEOM_COORD_LIMIT

}  // namespace

// ---------------------------------------------------------------------------
// validate: the order of the refusals
// ---------------------------------------------------------------------------

TEST(ArcAttack, AnAngleUnitPastTheLastIsRefusedBeforeACoordinateLimit) {
  GSHGEOM_Arc a = unit_arc();
  a.sweep = 1;  // within a full turn in whatever unit the record names, so the sweep check passes on its own
  a.unit = (GSHGEOM_Angle_Unit)GSHGEOM_ANGLE_UNIT_COUNT;
  a.center = P(kBeyond, 0);
  EXPECT_EQ(append(a), GSHGEOM_ERR_INVALID);
  a.unit = (GSHGEOM_Angle_Unit)77;
  EXPECT_EQ(append(a), GSHGEOM_ERR_INVALID);
  a = unit_arc();
  a.center = P(kBeyond, 0);
  EXPECT_EQ(append(a), GSHGEOM_ERR_LIMIT);  // the limit alone is a limit
}

TEST(ArcAttack, TheToleranceOfTheAgreementCountsTheCentre) {
  // The scale of the agreement bound is the largest of the end points, the centre and the half-diameters.  With the centre
  // the largest, a start that is off by more than 2^-40 of the next largest but less than 2^-40 of the centre is accepted.
  const double M = 1000.0;
  GSHGEOM_Arc a;
  std::memset(&a, 0, sizeof a);
  a.u = P(M, 0);
  a.v = P(0, M);
  a.center = P(2.5 * M, 0);  // the largest magnitude of the record
  a.unit = GSHGEOM_ANGLE_DEGREES;
  a.t0 = 180;
  a.sweep = -45;
  GSHGEOM_Arc_Form f = {a.center, a.u, a.v, a.unit, a.t0, a.sweep, 0u};
  ASSERT_EQ(gshgeom_arc_point(&f, a.t0, &a.start), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_arc_point(&f, a.t0 + a.sweep, &a.end), GSHGEOM_OK);
  ASSERT_LT(std::fabs(a.end.x), 2.5 * M);
  const double without_centre = 0x1p-40 * (2 * M);   // the scale the half-diameters make
  const double with_centre = 0x1p-40 * (2.5 * M);
  const double on = a.start.x;
  a.start.x = on + (without_centre + with_centre) / 2;  // off the exact start by a value between the two bounds
  EXPECT_EQ(append(a), GSHGEOM_OK);
  a.start.x = on + 2 * with_centre;                     // and past both
  EXPECT_EQ(append(a), GSHGEOM_ERR_INVALID);
}

TEST(ArcAttack, AnEndPastTheLimitOfQuarterTurnsIsALimitWhenTheStartIsNot) {
  // t0 is exactly at the limit of 2^40 quarter turns and the sweep carries the end one quarter past it: the point of the end cannot be
  // read, and the record is a limit (the code of the angle), not a disagreement of the end
  const double t0 = 90.0 * 0x1p40;
  GSHGEOM_Arc a = unit_arc();
  a.t0 = t0;
  a.sweep = 90;
  a.start = P(1, 0);   // t0 is a whole number of turns
  a.end = P(0, 1);
  EXPECT_EQ(append(a), GSHGEOM_ERR_LIMIT);
  GSHGEOM_Arc b = a;
  b.sweep = -90;       // the other way, the end is inside the limit and the record is as good as any
  b.end = P(0, -1);
  EXPECT_EQ(append(b), GSHGEOM_OK);
  // the same through the conjugate form, which reads the end itself
  GSHGEOM_Path * p = make();
  GSHGEOM_Arc_Form f = {P(0, 0), P(1, 0), P(0, 1), GSHGEOM_ANGLE_DEGREES, t0, 90, GSHGEOM_ARC_CIRCULAR};
  EXPECT_EQ(gshgeom_path_arc_conjugate_to(p, &f, nullptr), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(gshgeom_path_element_count(p), 0u);
  gshgeom_path_destroy(p);
}

// ---------------------------------------------------------------------------
// gshgeom_arc_pieces and the lowering
// ---------------------------------------------------------------------------

TEST(ArcAttack, PiecesRefuseANullOutputAndWriteNothingOnARefusal) {
  GSHGEOM_Arc a = unit_arc();
  int n = 5;
  EXPECT_EQ(gshgeom_arc_pieces(&a, nullptr, &n), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(n, 5);
  GSHGEOM_Curve out[GSHGEOM_ARC_MAX_PIECES];
  std::memset(out, 0x5a, sizeof out);
  GSHGEOM_Curve before[GSHGEOM_ARC_MAX_PIECES];
  std::memcpy(before, out, sizeof out);
  const double nan = std::nan("");
  const double inf = std::numeric_limits<double>::infinity();
  for (int which = 0; which < 5; which++) {
    GSHGEOM_Arc b = a;
    if (which == 0) b.unit = (GSHGEOM_Angle_Unit)GSHGEOM_ANGLE_UNIT_COUNT;
    if (which == 1) b.sweep = nan;
    if (which == 2) b.sweep = inf;
    if (which == 3) b.t0 = nan;
    if (which == 4) b.t0 = -inf;
    n = 5;
    EXPECT_EQ(gshgeom_arc_pieces(&b, out, &n), GSHGEOM_ERR_INVALID) << which;
    EXPECT_EQ(n, 5) << which;
    EXPECT_EQ(std::memcmp(out, before, sizeof out), 0) << which;
  }
}

TEST(ArcAttack, TheQuarterPieceWeightIsTheRoundedRootOfAHalf) {
  EXPECT_EQ(bits(gmath_sqrt(0.5)), bits(0.7071067811865476));
  GSHGEOM_Arc a = unit_arc();
  GSHGEOM_Curve out[GSHGEOM_ARC_MAX_PIECES];
  int n = 0;
  ASSERT_EQ(gshgeom_arc_pieces(&a, out, &n), GSHGEOM_OK);
  ASSERT_EQ(n, 1);
  EXPECT_EQ(bits(out[0].w[1]), bits(0.7071067811865476));
}

TEST(ArcAttack, ALoweringReportsAnArcThatCannotBeDecomposed) {
  // an arc record stored raw, as no constructor stores one (a NaN sweep): the lowering reads it back and must not skip it
  GSHGEOM_Path * p = make();
  ASSERT_EQ(gshgeom_path_move_to(p, P(1, 0)), GSHGEOM_OK);
  GSHGEOM_Point pts[4] = {P(0, 0), P(1, 0), P(0, 1), P(0, 1)};
  const double w[2] = {0.0, std::nan("")};
  ASSERT_EQ(gshgeom_path_append_raw(p, GSHGEOM_VERB_ARC, (unsigned)GSHGEOM_ANGLE_DEGREES, GSHGEOM_ARC_CIRCULAR, pts, w), GSHGEOM_OK);
  GSHGEOM_Path * lowered = nullptr;
  EXPECT_EQ(gshgeom_path_lower_arcs(p, nullptr, nullptr, &lowered), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(lowered, nullptr);
  gshgeom_path_destroy(p);
}

// ---------------------------------------------------------------------------
// gshgeom_arc_map
// ---------------------------------------------------------------------------

TEST(ArcAttack, AnInvalidMapIsReportedByTheMapCheckBeforeTheKindIsRead) {
  GSHGEOM_Map m = identity();
  m.m[0][0] = kBeyond;                  // past GSHGEOM_MAP_LIMIT: the map's check says LIMIT
  m.kind |= GSHGEOM_MAP_SINGULAR;       // and the tag says singular, which would say INVALID
  const GSHGEOM_Arc a = unit_arc();
  GSHGEOM_Arc out;
  EXPECT_EQ(gshgeom_arc_map(&m, &a, &out), GSHGEOM_ERR_LIMIT);
}

TEST(ArcAttack, ASingularMapIsRefusedWhateverTheRecord) {
  // the tag alone refuses: the entries are those of the identity, which would map the record to itself
  GSHGEOM_Map m = identity();
  m.kind |= GSHGEOM_MAP_SINGULAR;
  const GSHGEOM_Arc a = unit_arc();
  GSHGEOM_Arc out;
  EXPECT_EQ(gshgeom_arc_map(&m, &a, &out), GSHGEOM_ERR_INVALID);
}

TEST(ArcAttack, AFailedMappingOfAPointIsReportedWhateverTheLaterCheckSays) {
  // a centre that overflows under the map (the map's own result is not finite: LIMIT), with a record whose half-diameters are
  // collinear, which the validity check of the mapped record would call INVALID
  GSHGEOM_Map m = identity();
  m.m[0][0] = 0x1p100;
  m.m[1][1] = 0x1p100;
  m.kind = GSHGEOM_MAP_SCALES;
  for (int which = 0; which < 3; which++) {
    GSHGEOM_Arc a = unit_arc();
    a.v = P(2, 0);                       // collinear with u: no area
    if (which == 0) a.center = P(1e300, 0);
    if (which == 1) a.start = P(1e300, 0);
    if (which == 2) a.end = P(1e300, 0);
    GSHGEOM_Arc out;
    EXPECT_EQ(gshgeom_arc_map(&m, &a, &out), GSHGEOM_ERR_LIMIT) << which;
  }
  for (int which = 0; which < 2; which++) {
    GSHGEOM_Arc a = unit_arc();
    a.v = P(2, 0);
    if (which == 0) a.u = P(1e300, 0);
    if (which == 1) a.v = P(1e300, 0);
    GSHGEOM_Arc out;
    EXPECT_EQ(gshgeom_arc_map(&m, &a, &out), GSHGEOM_ERR_LIMIT) << "vector " << which;
  }
}

TEST(ArcAttack, AGeneralMapDoesNotKeepTheCircularFlag) {
  GSHGEOM_Map m = identity();
  m.kind = GSHGEOM_MAP_GENERAL;  // a tag of no known structure, over entries that scale one axis
  m.m[1][1] = 2.0;
  GSHGEOM_Arc a = unit_arc();
  GSHGEOM_Arc out;
  ASSERT_EQ(gshgeom_arc_map(&m, &a, &out), GSHGEOM_OK);
  EXPECT_EQ(out.flags & GSHGEOM_ARC_CIRCULAR, 0u);
  // the same tag on a map that only moves: also not similar by its tag
  GSHGEOM_Map t = identity();
  t.kind = GSHGEOM_MAP_GENERAL;
  ASSERT_EQ(gshgeom_arc_map(&t, &a, &out), GSHGEOM_OK);
  EXPECT_EQ(out.flags & GSHGEOM_ARC_CIRCULAR, 0u);
}

// ---------------------------------------------------------------------------
// The endpoint form
// ---------------------------------------------------------------------------

TEST(ArcAttack, AnEndpointArcWhoseRadiiAreFarBeyondTheBoundHasOneProductOverflowRefusedAsALimit) {
  // the product of a radius and the cosine of the angle is finite or infinite, never a NaN, and an infinity is past the limit
  GSHGEOM_Path * p = make();
  ASSERT_EQ(gshgeom_path_move_to(p, P(0, 0)), GSHGEOM_OK);
  GSHGEOM_Arc_Kind kind;
  EXPECT_EQ(gshgeom_path_arc_endpoint_to(p, P(1, 0), DBL_MAX, DBL_MAX, 0.0, GSHGEOM_ANGLE_RADIANS, 0, 1, &kind), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(gshgeom_path_arc_endpoint_to(p, P(1, 0), DBL_MAX, 1.0, 45.0, GSHGEOM_ANGLE_DEGREES, 0, 1, &kind), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(gshgeom_path_element_count(p), 1u);
  gshgeom_path_destroy(p);
}

TEST(ArcAttack, AnEndpointArcWhoseRadiiJustFitIsTheHalfTurnAboutTheMidpoint) {
  // N = 0 exactly: the bound's equality, which the scaled-radii branch and the fitting branch must give the same record for
  const double s[][4] = {{0, 0, 2, 0}, {0, 0, 0, 6}, {-3, 4, 3, 4}, {1, 1, 4, 5}};
  const double r[] = {1.0, 3.0, 3.0, 2.5};
  for (int i = 0; i < 4; i++) {
    for (int sweep = 0; sweep < 2; sweep++) {
      for (int large = 0; large < 2; large++) {
        GSHGEOM_Path * p = make();
        ASSERT_EQ(gshgeom_path_move_to(p, P(s[i][0], s[i][1])), GSHGEOM_OK);
        GSHGEOM_Arc_Kind kind;
        ASSERT_EQ(gshgeom_path_arc_endpoint_to(p, P(s[i][2], s[i][3]), r[i], r[i], 0.0, GSHGEOM_ANGLE_RADIANS, large, sweep, &kind),
            GSHGEOM_OK);
        GSHGEOM_Arc a;
        ASSERT_EQ(gshgeom_path_arc(p, 1, &a), GSHGEOM_OK);
        EXPECT_EQ(bits(a.center.x), bits((s[i][0] * 0.5) + (s[i][2] * 0.5)));
        EXPECT_EQ(bits(a.center.y), bits((s[i][1] * 0.5) + (s[i][3] * 0.5)));
        EXPECT_EQ(bits(a.u.x), bits(r[i])) << i;
        EXPECT_EQ(bits(std::fabs(a.sweep)), bits(3.141592653589793)) << i << " " << sweep << " " << large;
        EXPECT_EQ(a.sweep > 0, sweep == 1);
        gshgeom_path_destroy(p);
      }
    }
  }
}

TEST(ArcAttack, TheMidpointOfTheEndsIsTheSameFromTheSumOfHalvesAndHalfTheSum) {
  // The two spellings differ only when half of an end is subnormal and rounds; the other end is then the same size or the
  // sum is lost under it, and the chord's squared length underflows to zero in the exact sums' conversion, which the
  // endpoint form refuses.  So no accepted record can tell them apart: the tiny chords are refused, and across the whole
  // range of the ends the two spellings agree.
  const double u = 0x1p-1074;
  for (int k = 1; k < 40; k += 2) {
    GSHGEOM_Path * p = make();
    ASSERT_EQ(gshgeom_path_move_to(p, P(1 * u, 0)), GSHGEOM_OK);
    GSHGEOM_Arc_Kind kind;
    EXPECT_EQ(gshgeom_path_arc_endpoint_to(p, P((double)(k + 4) * u, 0), 1.0, 1.0, 0.0, GSHGEOM_ANGLE_RADIANS, 0, 1, &kind),
        GSHGEOM_ERR_INVALID);
    gshgeom_path_destroy(p);
  }
  uint64_t seed = 88172645463325252ull;
  auto next = [&]() {
    seed ^= seed << 13;
    seed ^= seed >> 7;
    seed ^= seed << 17;
    return seed;
  };
  for (int i = 0; i < 20000; i++) {
    const double a = std::ldexp((double)(next() >> 11) + 1.0, (int)(next() % 1200) - 1130);
    const double b = std::ldexp((double)(next() >> 11) + 1.0, (int)(next() % 1200) - 1130);
    volatile double halves = (a * 0.5) + (b * 0.5);
    volatile double sum = a + b;
    volatile double half_sum = sum * 0.5;
    if (std::fabs(a - b) >= 0x1p-500) {  // a chord the endpoint form accepts
      ASSERT_EQ(bits(halves), bits(half_sum)) << a << " " << b;
    }
  }
}

TEST(ArcAttack, ALargeArcWithAnEndAFewUlpsAwayKeepsTheSweepOneStepBelowTheTurn) {
  GSHGEOM_Path * p = make();
  ASSERT_EQ(gshgeom_path_move_to(p, P(0, 0)), GSHGEOM_OK);
  GSHGEOM_Arc_Kind kind;
  ASSERT_EQ(gshgeom_path_arc_endpoint_to(p, P(1e-17, 0), 1.0, 1.0, 0.0, GSHGEOM_ANGLE_RADIANS, 1, 1, &kind), GSHGEOM_OK);
  GSHGEOM_Arc a;
  ASSERT_EQ(gshgeom_path_arc(p, 1, &a), GSHGEOM_OK);
  EXPECT_EQ(bits(a.sweep), bits(0x1.921fb54442d17p+2));
  gshgeom_path_destroy(p);
}

// ---------------------------------------------------------------------------
// The centre and conjugate forms
// ---------------------------------------------------------------------------

TEST(ArcAttack, AFormsCentreHeightIsNotKept) {
  GSHGEOM_Path * p = make();
  GSHGEOM_Arc_Form f = {GSHGEOM_Point{1, 2, 7.0}, P(3, 0), P(0, 3), GSHGEOM_ANGLE_DEGREES, 0, 90, GSHGEOM_ARC_CIRCULAR};
  ASSERT_EQ(gshgeom_path_arc_conjugate_to(p, &f, nullptr), GSHGEOM_OK);
  GSHGEOM_Arc a;
  ASSERT_EQ(gshgeom_path_arc(p, 1, &a), GSHGEOM_OK);
  EXPECT_EQ(bits(a.center.z), bits(0.0));
  gshgeom_path_destroy(p);
}

TEST(ArcAttack, AnArcThatNeedsAMoveAlsoNeedsRoomForIt) {
  // with room for one more element only, an arc on a path with no current point needs two (the move and the arc): LIMIT,
  // whatever else is wrong with the arc (the end does not agree with the form: INVALID where the room is not looked at first)
  GSHGEOM_Limits lim;
  std::memset(&lim, 0, sizeof lim);
  lim.max_elements = 1;
  GSHGEOM_Path * p = make(&lim);
  GSHGEOM_Arc_Form f = {P(0, 0), P(1, 0), P(0, 1), GSHGEOM_ANGLE_DEGREES, 0, 90, GSHGEOM_ARC_CIRCULAR};
  const GSHGEOM_Point wrong_end = P(5, 5);
  EXPECT_EQ(gshgeom_path_arc_conjugate_to(p, &f, &wrong_end), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(gshgeom_path_element_count(p), 0u);
  EXPECT_EQ(gshgeom_path_arc_conjugate_to(p, &f, nullptr), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(gshgeom_path_element_count(p), 0u);
  gshgeom_path_destroy(p);
  lim.max_elements = 2;
  p = make(&lim);
  EXPECT_EQ(gshgeom_path_arc_conjugate_to(p, &f, &wrong_end), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_element_count(p), 0u);
  EXPECT_EQ(gshgeom_path_arc_conjugate_to(p, &f, nullptr), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_element_count(p), 2u);
  gshgeom_path_destroy(p);
}

TEST(ArcAttack, TheCentreFormRefusesWhatItsChecksCoverAndSaysINVALID) {
  GSHGEOM_Path * p = make();
  const double nan = std::nan("");
  const double inf = std::numeric_limits<double>::infinity();
  EXPECT_EQ(gshgeom_path_arc_center_to(p, P(0, 0), 1, 1, nan, GSHGEOM_ANGLE_DEGREES, 0, 90, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_center_to(p, P(0, 0), 1, 1, inf, GSHGEOM_ANGLE_DEGREES, 0, 90, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_center_to(p, P(0, 0), 1, 1, 0, (GSHGEOM_Angle_Unit)GSHGEOM_ANGLE_UNIT_COUNT, 0, 90, nullptr),
      GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_element_count(p), 0u);
  gshgeom_path_destroy(p);
}

// ---------------------------------------------------------------------------
// The tangent form
// ---------------------------------------------------------------------------

TEST(ArcAttack, TheTangentFormNeedsACurrentPointAndRefusesWithoutOne) {
  GSHGEOM_Path * p = make();
  EXPECT_EQ(gshgeom_path_arc_tangent_to(p, P(1, 0), P(1, 1), 0.5), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_tangent_to(p, P(1, 0), P(2, 0), 0.0), GSHGEOM_ERR_INVALID);   // the line case too
  EXPECT_EQ(gshgeom_path_arc_tangent_to(p, P(0, 0), P(0, 0), 0.5), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_element_count(p), 0u);
  gshgeom_path_destroy(p);
}

TEST(ArcAttack, TheTangentFormRefusesANonFinitePointAndAddsNothing) {
  const double nan = std::nan("");
  const double inf = std::numeric_limits<double>::infinity();
  GSHGEOM_Path * p = make();
  ASSERT_EQ(gshgeom_path_move_to(p, P(0, 0)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_arc_tangent_to(p, P(nan, 0), P(1, 1), 0.5), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_tangent_to(p, P(inf, 0), P(1, 1), 0.5), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_tangent_to(p, P(1, nan), P(1, 1), 0.5), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_tangent_to(p, P(1, 0), P(nan, 1), 0.5), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_tangent_to(p, P(1, 0), P(1, nan), 0.5), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_tangent_to(p, P(1, 0), P(1, -inf), 0.5), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_element_count(p), 1u);
  gshgeom_path_destroy(p);
}

TEST(ArcAttack, TheTangentFormReportsAnOrientationThatCannotBeDecided) {
  // the second point is past the coordinate limit, which the orientation test refuses; the first is fine, so a corner
  // taken as a straight line to it would succeed
  GSHGEOM_Path * p = make();
  ASSERT_EQ(gshgeom_path_move_to(p, P(0, 0)), GSHGEOM_OK);
  const GSHGEOM_Result r = gshgeom_path_arc_tangent_to(p, P(1, 0), P(1e300, 5), 0.5);
  EXPECT_NE(r, GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_element_count(p), 1u);
  gshgeom_path_destroy(p);
}

TEST(ArcAttack, TheTangentPointEqualToTheCurrentPointByValueButNotByBitsStillGetsItsLine) {
  // the current point is (-0, 0) and the tangent point comes out (+0, 0): the same value, different bits, so the line to it
  // is added (a degenerate one), as for any tangent point that is not the current point bit for bit
  GSHGEOM_Path * p = make();
  ASSERT_EQ(gshgeom_path_move_to(p, P(-0.0, 0.0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_arc_tangent_to(p, P(0.0, 1.0), P(1.0, 1.0), 1.0), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_element_count(p), 3u);
  gshgeom_path_destroy(p);
}

TEST(ArcAttack, TheTangentFormGivesALineWhenTheCornerIsDegenerate) {
  // the first point equal to the current point, or to the second, in value and in bits, and in value only
  GSHGEOM_Path * p = make();
  ASSERT_EQ(gshgeom_path_move_to(p, P(0, 0)), GSHGEOM_OK);
  GSHGEOM_Arc_Kind unused;
  (void)unused;
  EXPECT_EQ(gshgeom_path_arc_tangent_to(p, P(0, 0), P(3, 3), 1.0), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_arc_tangent_to(p, P(-0.0, -0.0), P(3, 3), 1.0), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_arc_tangent_to(p, P(5, 5), P(5, 5), 1.0), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_arc_tangent_to(p, P(5, 5), P(-5 + 10.0, 5), 1.0), GSHGEOM_OK);
  GSHGEOM_Verb v;
  for (size_t i = 1; i < gshgeom_path_element_count(p); i++) {
    ASSERT_EQ(gshgeom_path_verb(p, i, &v), GSHGEOM_OK);
    EXPECT_EQ(v, GSHGEOM_VERB_LINE) << i;
  }
  gshgeom_path_destroy(p);
}

// ---------------------------------------------------------------------------
// The bulge form
// ---------------------------------------------------------------------------

TEST(ArcAttack, TheBulgeFormRefusesANonFiniteEndForAStraightSegmentToo) {
  const double inf = std::numeric_limits<double>::infinity();
  const double nan = std::nan("");
  GSHGEOM_Path * p = make();
  ASSERT_EQ(gshgeom_path_move_to(p, P(0, 0)), GSHGEOM_OK);
  for (double b : {0.0, 1.0}) {
    EXPECT_EQ(gshgeom_path_arc_bulge_to(p, P(inf, 0), b), GSHGEOM_ERR_INVALID) << b;
    EXPECT_EQ(gshgeom_path_arc_bulge_to(p, P(0, -inf), b), GSHGEOM_ERR_INVALID) << b;
    EXPECT_EQ(gshgeom_path_arc_bulge_to(p, P(nan, 0), b), GSHGEOM_ERR_INVALID) << b;
    EXPECT_EQ(gshgeom_path_arc_bulge_to(p, P(0, nan), b), GSHGEOM_ERR_INVALID) << b;
  }
  EXPECT_EQ(gshgeom_path_element_count(p), 1u);
  gshgeom_path_destroy(p);
}

TEST(ArcAttack, TheBulgeCentreIsTheSumOfTheHalvesOfTheEnds) {
  const double u = 0x1p-1074;
  GSHGEOM_Path * p = make();
  ASSERT_EQ(gshgeom_path_move_to(p, P(1 * u, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_arc_bulge_to(p, P(5 * u, 0), 1e-300), GSHGEOM_OK);
  GSHGEOM_Arc a;
  ASSERT_EQ(gshgeom_path_arc(p, 1, &a), GSHGEOM_OK);
  volatile double h1 = (1 * u) * 0.5, h2 = (5 * u) * 0.5;
  EXPECT_EQ(bits(a.center.x), bits(h1 + h2));
  gshgeom_path_destroy(p);
}

TEST(ArcAttack, TheArcReadbackOfEveryFormIsInsensitiveToTheFormsFlagsField) {
  // the form's flags are read by the constructors that build an arc from a form and by nothing that reads one back
  GSHGEOM_Arc a = unit_arc();
  GSHGEOM_Arc_Form f = {a.center, a.u, a.v, a.unit, a.t0, a.sweep, 0u};
  GSHGEOM_Point p0, p1;
  GSHGEOM_Arc_Form g = f;
  g.flags = GSHGEOM_ARC_CIRCULAR;
  ASSERT_EQ(gshgeom_arc_point(&f, 30.0, &p0), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_arc_point(&g, 30.0, &p1), GSHGEOM_OK);
  EXPECT_EQ(bits(p0.x), bits(p1.x));
  GSHGEOM_Curve o1[GSHGEOM_ARC_MAX_PIECES], o2[GSHGEOM_ARC_MAX_PIECES];
  int n1 = 0, n2 = 0;
  GSHGEOM_Arc a2 = a;
  a2.flags = 0u;
  ASSERT_EQ(gshgeom_arc_pieces(&a, o1, &n1), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_arc_pieces(&a2, o2, &n2), GSHGEOM_OK);
  EXPECT_EQ(n1, n2);
  EXPECT_EQ(std::memcmp(o1, o2, sizeof(GSHGEOM_Curve) * (size_t)n1), 0);
}

TEST(ArcAttack, APointOfAFormWithAnUnitPastTheLastIsRefusedAndWritesNothing) {
  GSHGEOM_Arc_Form f = {P(0, 0), P(1, 0), P(0, 1), (GSHGEOM_Angle_Unit)GSHGEOM_ANGLE_UNIT_COUNT, 0, 90, 0u};
  GSHGEOM_Point out = GSHGEOM_Point{7, 8, 9};
  EXPECT_EQ(gshgeom_arc_point(&f, 0.0, &out), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(out.x, 7.0);
  EXPECT_EQ(out.y, 8.0);
  EXPECT_EQ(out.z, 9.0);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
