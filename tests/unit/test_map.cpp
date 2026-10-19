/**
 * @file
 *
 * Maps (CAP-2, CAP-4, AD-24): quarter turns that are exact, the sense of
 * rotation, the constructors and their kind tags, the one composer and its
 * order, tags that survive composition, application, the route predicate, the
 * hash and dump pinned by known answers, and the chiral fixture that pins the
 * orientation rules.
 *
 * Copyright 2026 by Corey Pennycuff
 */

#include "test_helpers.h"

#include "../known_answers.h"
#include "../dump_sink.h"

#include <ghoti.io/math/math.h>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

namespace {

const double kNaN = std::numeric_limits<double>::quiet_NaN();
const double kInf = std::numeric_limits<double>::infinity();

enum : unsigned {
  TR = GSHGEOM_MAP_TRANSLATES,
  SC = GSHGEOM_MAP_SCALES,
  NU = GSHGEOM_MAP_SCALES_NONUNIFORM,
  FY = GSHGEOM_MAP_FLIPS_Y,
  MI = GSHGEOM_MAP_MIRRORS,
  QT = GSHGEOM_MAP_QUARTER_TURNS,
  RO = GSHGEOM_MAP_ROTATES,
  SK = GSHGEOM_MAP_SKEWS,
  SG = GSHGEOM_MAP_SINGULAR,
  GE = GSHGEOM_MAP_GENERAL,
};

GSHGEOM_Point P(double x, double y, double z = 0.0) {
  return GSHGEOM_Point{x, y, z};
}

uint64_t bits(double v) {
  uint64_t b;
  std::memcpy(&b, &v, sizeof b);
  return b;
}

// Bit-exact equality, so that a negative zero is not a zero.
::testing::AssertionResult SameBits(double a, double b) {
  if (bits(a) == bits(b)) {
    return ::testing::AssertionSuccess();
  }
  return ::testing::AssertionFailure() << a << " (" << std::hex << bits(a)
                                       << ") vs " << b << " (" << bits(b) << ")";
}

GSHGEOM_Map Identity(int dim) {
  GSHGEOM_Map m;
  EXPECT_EQ(gshgeom_map_identity(dim, &m), GSHGEOM_OK);
  return m;
}

GSHGEOM_Map Translate(double x, double y, double z = 0.0, int dim = 2) {
  GSHGEOM_Map m;
  EXPECT_EQ(gshgeom_map_translate(dim, x, y, z, &m), GSHGEOM_OK);
  return m;
}

GSHGEOM_Map Scale(double x, double y, double z = 1.0, int dim = 2) {
  GSHGEOM_Map m;
  EXPECT_EQ(gshgeom_map_scale(dim, x, y, z, &m), GSHGEOM_OK);
  return m;
}

GSHGEOM_Map Rotate(GSHGEOM_Angle_Unit u, double a) {
  GSHGEOM_Map m;
  EXPECT_EQ(gshgeom_map_rotate(u, a, &m), GSHGEOM_OK);
  return m;
}

GSHGEOM_Map Matrix(int dim, std::vector<double> m, std::vector<double> t) {
  GSHGEOM_Map r;
  EXPECT_EQ(gshgeom_map_from_matrix(dim, m.data(), t.data(), &r), GSHGEOM_OK);
  return r;
}

GSHGEOM_Map Compose(const GSHGEOM_Map & outer, const GSHGEOM_Map & inner) {
  GSHGEOM_Map r;
  EXPECT_EQ(gshgeom_map_compose(&outer, &inner, &r), GSHGEOM_OK);
  return r;
}

GSHGEOM_Point Apply(const GSHGEOM_Map & m, GSHGEOM_Point p) {
  GSHGEOM_Point r{-1, -1, -1};
  EXPECT_EQ(gshgeom_map_apply_point(&m, p, &r), GSHGEOM_OK);
  return r;
}

bool Aligned(const GSHGEOM_Map & m, GSHGEOM_Grid * g = nullptr) {
  int a = -1;
  EXPECT_EQ(gshgeom_map_grid_aligned(&m, &a, g), GSHGEOM_OK);
  return a == 1;
}

::testing::AssertionResult SameMap(const GSHGEOM_Map & a, const GSHGEOM_Map & b) {
  if (a.dim != b.dim || a.kind != b.kind) {
    return ::testing::AssertionFailure() << "dim/kind " << a.dim << "/" << a.kind
                                         << " vs " << b.dim << "/" << b.kind;
  }
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (bits(a.m[i][j]) != bits(b.m[i][j])) {
        return ::testing::AssertionFailure() << "m[" << i << "][" << j << "] "
                                             << a.m[i][j] << " vs " << b.m[i][j];
      }
    }
    if (bits(a.t[i]) != bits(b.t[i])) {
      return ::testing::AssertionFailure() << "t[" << i << "]";
    }
  }
  return ::testing::AssertionSuccess();
}

std::string DumpOf(const GSHGEOM_Map & m) {
  FILE * f = std::tmpfile();
  EXPECT_NE(f, nullptr);
  EXPECT_EQ(gshgeom_map_dump(&m, f), GSHGEOM_OK);
  std::rewind(f);
  std::string s;
  char buf[256];
  size_t n;
  while ((n = std::fread(buf, 1, sizeof buf, f)) > 0) {
    s.append(buf, n);
  }
  std::fclose(f);
  return s;
}

uint64_t rng_state = UINT64_C(0x9E3779B97F4A7C15);

uint64_t next() {
  rng_state ^= rng_state >> 12;
  rng_state ^= rng_state << 25;
  rng_state ^= rng_state >> 27;
  return rng_state * UINT64_C(2685821657736338717);
}

// A number with a full mantissa in [-4, 4): products of these round.
double wild() {
  return ((double)(int64_t)(next() >> 11) / 9007199254740992.0 - 0.5) * 8.0 + 0.1;
}

GSHGEOM_Map WildMap(int dim) {
  std::vector<double> m, t;
  for (int i = 0; i < dim * dim; i++) m.push_back(wild());
  for (int i = 0; i < dim; i++) t.push_back(wild());
  return Matrix(dim, m, t);
}

}  // namespace

// ---------------------------------------------------------------------------
// Rotations
// ---------------------------------------------------------------------------

struct QuarterCase {
  GSHGEOM_Angle_Unit unit;
  double angle;
  double c, s;  // exact cosine and sine
};

TEST(MapRotate, QuarterTurnsGiveExactZeroAndOne) {
  const GSHGEOM_Angle_Unit D = GSHGEOM_ANGLE_DEGREES;
  const GSHGEOM_Angle_Unit G = GSHGEOM_ANGLE_GRADIANS;
  const GSHGEOM_Angle_Unit T = GSHGEOM_ANGLE_TURNS;
  const QuarterCase cases[] = {
      {D, 0, 1, 0}, {D, -0.0, 1, 0}, {D, 90, 0, 1}, {D, 180, -1, 0},
      {D, 270, 0, -1}, {D, -90, 0, -1}, {D, -180, -1, 0}, {D, -270, 0, 1},
      {D, 360, 1, 0}, {D, 450, 0, 1}, {D, -360, 1, 0}, {D, 720, 1, 0},
      {D, 90.0 * 1000000.0, 1, 0}, {D, 90.0 * 1000001.0, 0, 1},
      {D, 90.0 * 1000002.0, -1, 0}, {D, 90.0 * 1000003.0, 0, -1},
      {D, 90.0 * 1099511627776.0, 1, 0},
      {G, 100, 0, 1}, {G, 200, -1, 0}, {G, 300, 0, -1}, {G, 400, 1, 0},
      {G, -100, 0, -1},
      {T, 0.25, 0, 1}, {T, 0.5, -1, 0}, {T, 0.75, 0, -1}, {T, 1, 1, 0},
      {T, -0.25, 0, -1}, {T, 2.25, 0, 1}, {T, -1.5, -1, 0},
  };
  for (const QuarterCase & q : cases) {
    GSHGEOM_Map m = Rotate(q.unit, q.angle);
    SCOPED_TRACE(::testing::Message() << "unit " << q.unit << " angle " << q.angle);
    EXPECT_TRUE(SameBits(m.m[0][0], q.c));
    EXPECT_TRUE(SameBits(m.m[0][1], 0.0 - q.s));
    EXPECT_TRUE(SameBits(m.m[1][0], q.s));
    EXPECT_TRUE(SameBits(m.m[1][1], q.c));
    EXPECT_TRUE(SameBits(m.t[0], 0.0));
    EXPECT_TRUE(SameBits(m.t[1], 0.0));
    EXPECT_EQ(m.dim, 2);
    bool identity = q.c == 1.0;
    EXPECT_EQ(m.kind, identity ? 0u : (unsigned)QT);
  }
}

TEST(MapRotate, ThePositiveAngleTurnsXTowardY) {
  GSHGEOM_Map m = Rotate(GSHGEOM_ANGLE_DEGREES, 90);
  GSHGEOM_Point r = Apply(m, P(1, 0));
  EXPECT_TRUE(SameBits(r.x, 0.0));
  EXPECT_TRUE(SameBits(r.y, 1.0));
  r = Apply(m, P(0, 1));
  EXPECT_TRUE(SameBits(r.x, -1.0));
  EXPECT_TRUE(SameBits(r.y, 0.0));
  // On a screen with y down that is clockwise: (1, 0) is to the right and (0, 1)
  // is below it.
  r = Apply(Rotate(GSHGEOM_ANGLE_TURNS, 0.5), P(3, 4));
  EXPECT_TRUE(SameBits(r.x, -3.0));
  EXPECT_TRUE(SameBits(r.y, -4.0));
}

TEST(MapRotate, OtherAnglesAreTheTrigFunctionsOfTheAngle) {
  GSHGEOM_Map m = Rotate(GSHGEOM_ANGLE_DEGREES, 30);
  EXPECT_NEAR(m.m[0][0], 0.86602540378443864676, 4e-16);
  EXPECT_NEAR(m.m[1][0], 0.5, 4e-16);
  EXPECT_NEAR(m.m[0][1], -0.5, 4e-16);
  EXPECT_EQ(m.kind, (unsigned)RO);
  m = Rotate(GSHGEOM_ANGLE_DEGREES, 120);
  EXPECT_NEAR(m.m[0][0], -0.5, 4e-16);
  EXPECT_NEAR(m.m[1][0], 0.86602540378443864676, 4e-16);
  m = Rotate(GSHGEOM_ANGLE_GRADIANS, 50);
  EXPECT_NEAR(m.m[0][0], 0.70710678118654752440, 4e-16);
  m = Rotate(GSHGEOM_ANGLE_TURNS, 1.0 / 12.0);
  EXPECT_NEAR(m.m[0][0], 0.86602540378443864676, 4e-16);
  m = Rotate(GSHGEOM_ANGLE_RADIANS, 1.0);
  EXPECT_NEAR(m.m[0][0], 0.54030230586813971740, 2e-16);
  EXPECT_NEAR(m.m[1][0], 0.84147098480789650665, 2e-16);
  // Just past a quarter turn: cos(90 deg + 1e-9 deg) = -sin(1e-9 deg), which
  // a rotation not reduced by the exact quarter would get wrong in its
  // leading digits.
  const double just_past = 90.0 + 1e-9;
  m = Rotate(GSHGEOM_ANGLE_DEGREES, just_past);
  double tiny = just_past - 90.0;  // exact, as the call formed it
  EXPECT_NEAR(m.m[0][0], -tiny * 0.017453292519943295, 1e-25);
  EXPECT_NEAR(m.m[1][0], 1.0, 1e-15);
  // Just short of a quarter turn: the nearest quarter is the one above, and
  // the cosine is the sine of what is left, which a reduction to the quarter
  // below would lose in the leading digits.
  const double just_short = 90.0 - 1e-9;
  m = Rotate(GSHGEOM_ANGLE_DEGREES, just_short);
  double left = 90.0 - just_short;
  EXPECT_NEAR(m.m[0][0], left * 0.017453292519943295, 1e-25);
  EXPECT_NEAR(m.m[1][0], 1.0, 1e-15);
  m = Rotate(GSHGEOM_ANGLE_TURNS, 0.25 - 1e-12);
  EXPECT_NEAR(m.m[0][0], (0.25 - (0.25 - 1e-12)) * 6.283185307179586, 1e-26);
  m = Rotate(GSHGEOM_ANGLE_GRADIANS, 100.0 + 1e-9);
  EXPECT_NEAR(m.m[0][0], -((100.0 + 1e-9) - 100.0) * 0.015707963267948967, 1e-25);
  // The reduction is by exact quarters on either side of zero.
  GSHGEOM_Map a = Rotate(GSHGEOM_ANGLE_DEGREES, 400.0);
  GSHGEOM_Map b = Rotate(GSHGEOM_ANGLE_DEGREES, 40.0);
  EXPECT_NEAR(a.m[0][0], b.m[0][0], 4e-16);
  EXPECT_NEAR(a.m[1][0], b.m[1][0], 4e-16);
  // A rotation is orthogonal and has unit determinant.
  for (int u = 0; u < GSHGEOM_ANGLE_UNIT_COUNT; u++) {
    for (int i = 0; i < 300; i++) {
      double angle = wild() * 40.0;
      GSHGEOM_Map r = Rotate((GSHGEOM_Angle_Unit)u, angle);
      EXPECT_NEAR(r.m[0][0] * r.m[0][0] + r.m[1][0] * r.m[1][0], 1.0, 4e-16);
      EXPECT_EQ(r.m[0][0], r.m[1][1]);
      EXPECT_EQ(r.m[0][1], 0.0 - r.m[1][0]);
    }
  }
}

TEST(MapRotate, EveryUnitAgreesWithTheTrigFunctionsOfTheRadians) {
  // Angles in every quadrant and well past a whole turn, against gmath's sine
  // and cosine of the angle in radians (the straightforward way, which loses
  // the exactness of the quarters but not the value).
  const double per_unit[4] = {1.0, 0x1.1df46a2529d39p-6, 0x1.015bf9217271ap-6, 0x1.921fb54442d18p+2};
  for (int u = 0; u < GSHGEOM_ANGLE_UNIT_COUNT; u++) {
    double quarter_units[4] = {0.0, 90.0, 100.0, 0.25};
    for (int i = 0; i < 400; i++) {
      double angle = wild() * 40.0;
      if (u != 0 && i % 4 == 0) {
        // Land near a quarter: a multiple of it plus a small remainder.
        angle = quarter_units[u] * (double)((int)(next() % 17) - 8) + wild() * quarter_units[u] / 8.0;
      }
      GSHGEOM_Map m = Rotate((GSHGEOM_Angle_Unit)u, angle);
      double radians = angle * per_unit[u];
      // The straightforward way rounds the product angle * (radians per unit),
      // which is an error of about 1.1e-16 times the radians; the library does not.
      double tol = 1e-15 * (1.0 + (radians < 0 ? -radians : radians));
      EXPECT_NEAR(m.m[0][0], gmath_cos(radians), tol) << u << " " << angle;
      EXPECT_NEAR(m.m[1][0], gmath_sin(radians), tol) << u << " " << angle;
    }
  }
  // Each quadrant at a remainder that is not 0 or 45 degrees: 30, 120, 210, 300.
  const double quad[4] = {30, 120, 210, 300};
  const double cs[4] = {0.86602540378443864676, -0.5, -0.86602540378443864676, 0.5};
  const double sn[4] = {0.5, 0.86602540378443864676, -0.5, -0.86602540378443864676};
  for (int k = 0; k < 4; k++) {
    GSHGEOM_Map m = Rotate(GSHGEOM_ANGLE_DEGREES, quad[k]);
    EXPECT_NEAR(m.m[0][0], cs[k], 4e-16) << quad[k];
    EXPECT_NEAR(m.m[1][0], sn[k], 4e-16) << quad[k];
    m = Rotate(GSHGEOM_ANGLE_DEGREES, quad[k] - 360.0);
    EXPECT_NEAR(m.m[0][0], cs[k], 4e-16) << quad[k] - 360.0;
    EXPECT_NEAR(m.m[1][0], sn[k], 4e-16) << quad[k] - 360.0;
  }
}

TEST(MapRotate, ARadianQuarterIsNotSnapped) {
  GSHGEOM_Map m = Rotate(GSHGEOM_ANGLE_RADIANS, 0x1.921fb54442d18p+0);
  EXPECT_NE(m.m[0][0], 0.0) << "the double nearest pi/2 is not a quarter turn";
  EXPECT_NEAR(m.m[0][0], 6.123233995736766e-17, 1e-31);
  EXPECT_EQ(m.kind, (unsigned)RO);
  m = Rotate(GSHGEOM_ANGLE_RADIANS, 0.0);
  EXPECT_EQ(m.kind, 0u);
  EXPECT_TRUE(SameBits(m.m[0][1], 0.0));
  m = Rotate(GSHGEOM_ANGLE_RADIANS, -0.0);
  EXPECT_EQ(m.kind, 0u);
  // A huge angle is fine in radians: gmath reduces it.
  m = Rotate(GSHGEOM_ANGLE_RADIANS, 1e300);
  EXPECT_NEAR(m.m[0][0] * m.m[0][0] + m.m[1][0] * m.m[1][0], 1.0, 4e-16);
}

TEST(MapRotate, Refusals) {
  GSHGEOM_Map sentinel;
  std::memset(&sentinel, 0x5a, sizeof sentinel);
  GSHGEOM_Map m = sentinel;
  for (int u = 0; u < GSHGEOM_ANGLE_UNIT_COUNT; u++) {
    EXPECT_EQ(gshgeom_map_rotate((GSHGEOM_Angle_Unit)u, kNaN, &m), GSHGEOM_ERR_INVALID);
    EXPECT_EQ(gshgeom_map_rotate((GSHGEOM_Angle_Unit)u, kInf, &m), GSHGEOM_ERR_INVALID);
    EXPECT_EQ(gshgeom_map_rotate((GSHGEOM_Angle_Unit)u, -kInf, &m), GSHGEOM_ERR_INVALID);
  }
  // Values of the enumeration's range that name no unit (a value outside the
  // range would be undefined behaviour to form in C++).
  EXPECT_EQ(gshgeom_map_rotate((GSHGEOM_Angle_Unit)GSHGEOM_ANGLE_UNIT_COUNT, 1.0, &m), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_rotate((GSHGEOM_Angle_Unit)5, 1.0, &m), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_rotate((GSHGEOM_Angle_Unit)7, 1.0, &m), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_rotate(GSHGEOM_ANGLE_DEGREES, 1.0, nullptr), GSHGEOM_ERR_INVALID);
  // More than 2^40 quarter turns in a unit that is not radians.
  EXPECT_EQ(gshgeom_map_rotate(GSHGEOM_ANGLE_DEGREES, 1e300, &m), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(gshgeom_map_rotate(GSHGEOM_ANGLE_DEGREES, 90.0 * (1099511627776.0 + 1.0), &m), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(gshgeom_map_rotate(GSHGEOM_ANGLE_DEGREES, -90.0 * (1099511627776.0 + 1.0), &m), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(gshgeom_map_rotate(GSHGEOM_ANGLE_TURNS, 1e300, &m), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(gshgeom_map_rotate(GSHGEOM_ANGLE_GRADIANS, -1e300, &m), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(std::memcmp(&m, &sentinel, sizeof m), 0) << "untouched by a refusal";
}

TEST(MapRotate, AboutAnAxisTurnsTheNextAxisTowardTheOneAfter) {
  GSHGEOM_Map m;
  ASSERT_EQ(gshgeom_map_rotate_axis(GSHGEOM_AXIS_X, GSHGEOM_ANGLE_DEGREES, 90, &m), GSHGEOM_OK);
  EXPECT_EQ(m.dim, 3);
  EXPECT_EQ(m.kind, (unsigned)QT);
  GSHGEOM_Point r = Apply(m, P(0, 1, 0));
  EXPECT_TRUE(SameBits(r.x, 0.0));
  EXPECT_TRUE(SameBits(r.y, 0.0));
  EXPECT_TRUE(SameBits(r.z, 1.0));
  ASSERT_EQ(gshgeom_map_rotate_axis(GSHGEOM_AXIS_Y, GSHGEOM_ANGLE_DEGREES, 90, &m), GSHGEOM_OK);
  r = Apply(m, P(0, 0, 1));
  EXPECT_TRUE(SameBits(r.x, 1.0));
  EXPECT_TRUE(SameBits(r.y, 0.0));
  EXPECT_TRUE(SameBits(r.z, 0.0));
  ASSERT_EQ(gshgeom_map_rotate_axis(GSHGEOM_AXIS_Z, GSHGEOM_ANGLE_DEGREES, 90, &m), GSHGEOM_OK);
  r = Apply(m, P(1, 0, 0));
  EXPECT_TRUE(SameBits(r.x, 0.0));
  EXPECT_TRUE(SameBits(r.y, 1.0));
  // The axis itself stays put, and an arbitrary angle is a rotation.
  ASSERT_EQ(gshgeom_map_rotate_axis(GSHGEOM_AXIS_Z, GSHGEOM_ANGLE_DEGREES, 33, &m), GSHGEOM_OK);
  EXPECT_EQ(m.kind, (unsigned)GE) << "a 3D rotation of no axis-aligned structure";
  r = Apply(m, P(0, 0, 5));
  EXPECT_EQ(r.z, 5.0);
  GSHGEOM_Map a, z;
  ASSERT_EQ(gshgeom_map_rotate_axis(GSHGEOM_AXIS_Z, GSHGEOM_ANGLE_DEGREES, 90, &a), GSHGEOM_OK);
  GSHGEOM_Map b = Rotate(GSHGEOM_ANGLE_DEGREES, 90);
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      EXPECT_TRUE(SameBits(a.m[i][j], b.m[i][j])) << "the 2D rotation is about z";
    }
  }
  z = a;
  EXPECT_EQ(gshgeom_map_rotate_axis((GSHGEOM_Axis)3, GSHGEOM_ANGLE_DEGREES, 90, &z), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_rotate_axis(GSHGEOM_AXIS_X, GSHGEOM_ANGLE_DEGREES, kNaN, &z), GSHGEOM_ERR_INVALID);
  EXPECT_TRUE(SameMap(z, a));
  EXPECT_EQ(gshgeom_map_rotate_axis(GSHGEOM_AXIS_X, GSHGEOM_ANGLE_DEGREES, 90, nullptr), GSHGEOM_ERR_INVALID);
}

// ---------------------------------------------------------------------------
// Constructors and tags
// ---------------------------------------------------------------------------

TEST(MapKind, ConstructorsSetTheTagFromWhatTheyBuild) {
  EXPECT_EQ(Identity(2).kind, 0u);
  EXPECT_EQ(Identity(3).kind, 0u);
  EXPECT_EQ(Translate(1, 2).kind, (unsigned)TR);
  EXPECT_EQ(Translate(0, 0, 9).kind, 0u) << "z is ignored in 2D";
  EXPECT_EQ(Translate(0, 0, 9, 3).kind, (unsigned)TR);
  EXPECT_EQ(Translate(0, 0, 0, 3).kind, 0u);
  EXPECT_EQ(Scale(2, 2).kind, (unsigned)SC);
  EXPECT_EQ(Scale(0.5, 0.5).kind, (unsigned)SC);
  EXPECT_EQ(Scale(1, 1).kind, 0u);
  EXPECT_EQ(Scale(2, 3).kind, (unsigned)NU);
  EXPECT_EQ(Scale(0, 1).kind, (unsigned)SG);
  EXPECT_EQ(Scale(1, 0).kind, (unsigned)SG);
  EXPECT_EQ(Scale(-1, 1).kind, (unsigned)MI);
  EXPECT_EQ(Scale(1, -1).kind, (unsigned)FY);
  EXPECT_EQ(Scale(-1, -1).kind, (unsigned)QT) << "a half turn";
  EXPECT_EQ(Scale(-2, -2).kind, (unsigned)(QT | SC));
  EXPECT_EQ(Scale(2, -2).kind, (unsigned)(FY | SC));
  EXPECT_EQ(Scale(-2, 2).kind, (unsigned)(MI | SC));
  EXPECT_EQ(Scale(2, -3).kind, (unsigned)(FY | NU));
  EXPECT_EQ(Scale(1, 1, 1, 3).kind, 0u);
  EXPECT_EQ(Scale(2, 2, 2, 3).kind, (unsigned)SC);
  EXPECT_EQ(Scale(2, 2, 3, 3).kind, (unsigned)NU);
  EXPECT_EQ(Scale(1, -1, 1, 3).kind, (unsigned)FY);
  EXPECT_EQ(Scale(-1, 1, 1, 3).kind, (unsigned)MI);
  EXPECT_EQ(Scale(1, 1, -1, 3).kind, (unsigned)MI);
  EXPECT_EQ(Scale(-1, -1, 1, 3).kind, (unsigned)QT);
  EXPECT_EQ(Scale(1, -1, -1, 3).kind, (unsigned)QT);
  EXPECT_EQ(Scale(-1, -1, -1, 3).kind, (unsigned)MI) << "a point reflection";
  GSHGEOM_Map f;
  ASSERT_EQ(gshgeom_map_flip_y(2, &f), GSHGEOM_OK);
  EXPECT_EQ(f.kind, (unsigned)FY);
  EXPECT_TRUE(SameBits(f.m[1][1], -1.0));
  EXPECT_TRUE(SameBits(f.m[0][0], 1.0));
  ASSERT_EQ(gshgeom_map_flip_y(3, &f), GSHGEOM_OK);
  EXPECT_EQ(f.kind, (unsigned)FY);
  EXPECT_TRUE(SameBits(f.m[2][2], 1.0));
}

TEST(MapKind, ConstructorsFillTheRestOfTheMap) {
  GSHGEOM_Map m = Translate(1.5, -2.5, 7.0);
  EXPECT_EQ(m.dim, 2);
  EXPECT_TRUE(SameBits(m.m[2][2], 1.0));
  EXPECT_TRUE(SameBits(m.m[0][2], 0.0));
  EXPECT_TRUE(SameBits(m.m[2][0], 0.0));
  EXPECT_TRUE(SameBits(m.t[2], 0.0));
  m = Translate(-0.0, 0.0);
  EXPECT_TRUE(SameBits(m.t[0], 0.0)) << "no negative zero is stored";
  m = Scale(-0.0, 2.0);
  EXPECT_TRUE(SameBits(m.m[0][0], 0.0));
}

TEST(MapKind, ConstructorsRefuse) {
  GSHGEOM_Map sentinel;
  std::memset(&sentinel, 0x5a, sizeof sentinel);
  GSHGEOM_Map m = sentinel;
  EXPECT_EQ(gshgeom_map_identity(1, &m), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_identity(4, &m), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_identity(2, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_translate(5, 0, 0, 0, &m), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_translate(2, kNaN, 0, 0, &m), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_translate(2, 0, kInf, 0, &m), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_translate(2, 0, 0, kNaN, &m), GSHGEOM_OK) << "z is not read in 2D";
  m = sentinel;
  EXPECT_EQ(gshgeom_map_translate(3, 0, 0, kNaN, &m), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_translate(2, 0x1p101, 0, 0, &m), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(gshgeom_map_translate(2, -0x1p101, 0, 0, &m), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(gshgeom_map_translate(2, 0x1p100, -0x1p100, 0, &m), GSHGEOM_OK) << "the bound itself is allowed";
  m = sentinel;
  EXPECT_EQ(gshgeom_map_translate(2, 0x1p101, kNaN, 0, &m), GSHGEOM_ERR_INVALID) << "a NaN beats a large number";
  EXPECT_EQ(gshgeom_map_scale(2, kNaN, 1, 1, &m), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_scale(2, 1, 0x1p101, 1, &m), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(gshgeom_map_scale(3, 1, 1, kInf, &m), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_scale(0, 1, 1, 1, &m), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_flip_y(7, &m), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(std::memcmp(&m, &sentinel, sizeof m), 0);
}

TEST(MapKind, ARawMatrixIsClassifiedExactly) {
  size_t n = sizeof kMapFixtures / sizeof kMapFixtures[0];
  ASSERT_GE(n, 17u);
  for (size_t i = 0; i < n; i++) {
    const MapFixture & f = kMapFixtures[i];
    SCOPED_TRACE(f.name);
    GSHGEOM_Map m;
    ASSERT_EQ(gshgeom_map_from_matrix(f.dim, f.m, f.t, &m), GSHGEOM_OK);
    EXPECT_EQ(m.kind, f.kind);
    EXPECT_EQ(m.dim, f.dim);
    for (int r = 0; r < f.dim; r++) {
      for (int c = 0; c < f.dim; c++) {
        EXPECT_TRUE(SameBits(m.m[r][c], f.m[r * f.dim + c]));
      }
      EXPECT_TRUE(SameBits(m.t[r], f.t[r]));
    }
    if (f.dim == 2) {
      EXPECT_TRUE(SameBits(m.m[2][2], 1.0));
      EXPECT_TRUE(SameBits(m.t[2], 0.0));
    }
  }
}

TEST(MapKind, TheSignOfTheDeterminantIsExact) {
  // [[1, 1 + 2^-52], [1 - 2^-53, 1]]: the product (1 + 2^-52)(1 - 2^-53) is
  // 1 + 2^-53 - 2^-105, which a double rounds to 1, so the naive determinant
  // is zero; exactly it is not.
  const double a = 1.0 + 0x1p-52;
  const double b = 1.0 - 0x1p-53;
  GSHGEOM_Map m = Matrix(2, {1, a, b, 1}, {0, 0});
  EXPECT_EQ(m.kind & SG, 0u) << "(1 + 2^-52)(1 - 2^-53) is not 1";
  // [[a, b], [c, d]] with ad == bc exactly is singular.
  m = Matrix(2, {3, 6, 0.5, 1}, {0, 0});
  EXPECT_EQ(m.kind, (unsigned)SG);
  // A determinant below a double: [[2^-600, 0], [0, 2^-600]] has det 2^-1200,
  // zero as a product of doubles, positive exactly.
  m = Matrix(2, {0x1p-600, 0, 0, 0x1p-600}, {0, 0});
  EXPECT_EQ(m.kind & SG, 0u);
  EXPECT_EQ(m.kind, (unsigned)SC);
  // 3D: a skew whose rows are exactly dependent is singular.
  m = Matrix(3, {1, 2, 3, 2, 4, 6, 0, 0, 1}, {0, 0, 0});
  EXPECT_EQ(m.kind, (unsigned)SG);
  m = Matrix(3, {1, 2, 3, 4, 5, 6, 7, 8, 9}, {0, 0, 0});
  EXPECT_EQ(m.kind, (unsigned)SG);
  m = Matrix(3, {1, 2, 3, 4, 5, 6, 7, 8, 10}, {0, 0, 0});
  EXPECT_EQ(m.kind, (unsigned)GE);
}

TEST(MapKind, ARawMatrixRefuses) {
  GSHGEOM_Map sentinel;
  std::memset(&sentinel, 0x5a, sizeof sentinel);
  GSHGEOM_Map m = sentinel;
  double id[4] = {1, 0, 0, 1};
  double t[2] = {0, 0};
  EXPECT_EQ(gshgeom_map_from_matrix(2, nullptr, t, &m), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_from_matrix(2, id, nullptr, &m), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_from_matrix(2, id, t, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_from_matrix(5, id, t, &m), GSHGEOM_ERR_INVALID);
  double nan_m[4] = {1, kNaN, 0, 1};
  EXPECT_EQ(gshgeom_map_from_matrix(2, nan_m, t, &m), GSHGEOM_ERR_INVALID);
  double big_m[4] = {1, 0, 0x1p101, 1};
  EXPECT_EQ(gshgeom_map_from_matrix(2, big_m, t, &m), GSHGEOM_ERR_LIMIT);
  double nan_t[2] = {0, kInf};
  EXPECT_EQ(gshgeom_map_from_matrix(2, id, nan_t, &m), GSHGEOM_ERR_INVALID);
  double big_t[2] = {0x1p101, 0};
  EXPECT_EQ(gshgeom_map_from_matrix(2, id, big_t, &m), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(gshgeom_map_from_matrix(2, big_m, nan_t, &m), GSHGEOM_ERR_INVALID) << "a NaN beats a large number";
  EXPECT_EQ(std::memcmp(&m, &sentinel, sizeof m), 0);
}

TEST(MapKind, TheDesignToModelMap) {
  GSHGEOM_Map m;
  ASSERT_EQ(gshgeom_map_design_to_model(2.5, 10, 20, &m), GSHGEOM_OK);
  EXPECT_EQ(m.kind, (unsigned)(TR | SC | FY));
  EXPECT_TRUE(SameBits(m.m[0][0], 2.5));
  EXPECT_TRUE(SameBits(m.m[1][1], -2.5));
  EXPECT_TRUE(SameBits(m.m[0][1], 0.0));
  EXPECT_TRUE(SameBits(m.m[1][0], 0.0));
  EXPECT_TRUE(SameBits(m.t[0], 10.0));
  EXPECT_TRUE(SameBits(m.t[1], 20.0));
  EXPECT_TRUE(SameMap(m, Matrix(2, {2.5, 0, 0, -2.5}, {10, 20})));
  // A y-up point (3, 4) lands 2.5 * 4 above the translation's y.
  GSHGEOM_Point r = Apply(m, P(3, 4));
  EXPECT_EQ(r.x, 17.5);
  EXPECT_EQ(r.y, 10.0);
  GSHGEOM_Map sentinel = m;
  EXPECT_EQ(gshgeom_map_design_to_model(0.0, 0, 0, &m), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_design_to_model(-1.0, 0, 0, &m), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_design_to_model(kNaN, 0, 0, &m), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_design_to_model(kInf, 0, 0, &m), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_design_to_model(1.0, kNaN, 0, &m), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_design_to_model(0x1p101, 0, 0, &m), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(gshgeom_map_design_to_model(1.0, 0x1p101, 0, &m), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(gshgeom_map_design_to_model(1.0, 0, 0, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_TRUE(SameMap(m, sentinel));
}

TEST(MapKind, CheckRefusesAnInvalidMap) {
  GSHGEOM_Map m = Identity(2);
  EXPECT_EQ(gshgeom_map_check(&m), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_map_check(nullptr), GSHGEOM_ERR_INVALID);
  GSHGEOM_Map bad = m;
  bad.dim = 4;
  EXPECT_EQ(gshgeom_map_check(&bad), GSHGEOM_ERR_INVALID);
  bad = m;
  bad.dim = 1;
  EXPECT_EQ(gshgeom_map_check(&bad), GSHGEOM_ERR_INVALID);
  bad = m;
  bad.kind = GSHGEOM_MAP_KIND_ALL + 1;
  EXPECT_EQ(gshgeom_map_check(&bad), GSHGEOM_ERR_INVALID);
  bad = m;
  bad.kind = GSHGEOM_MAP_KIND_ALL;
  EXPECT_EQ(gshgeom_map_check(&bad), GSHGEOM_OK) << "every defined bit is valid";
  bad = m;
  bad.m[1][0] = kNaN;
  EXPECT_EQ(gshgeom_map_check(&bad), GSHGEOM_ERR_INVALID);
  bad = m;
  bad.t[1] = -kInf;
  EXPECT_EQ(gshgeom_map_check(&bad), GSHGEOM_ERR_INVALID);
  bad = m;
  bad.m[0][1] = 0x1p101;
  EXPECT_EQ(gshgeom_map_check(&bad), GSHGEOM_ERR_LIMIT);
  bad = m;
  bad.t[0] = -0x1p101;
  EXPECT_EQ(gshgeom_map_check(&bad), GSHGEOM_ERR_LIMIT);
  bad = m;
  bad.m[0][1] = 0x1p101;
  bad.t[0] = kNaN;
  EXPECT_EQ(gshgeom_map_check(&bad), GSHGEOM_ERR_INVALID) << "a NaN beats a large number";
  // What a 2D map does not use is not looked at.
  bad = m;
  bad.m[2][2] = kNaN;
  bad.t[2] = kInf;
  EXPECT_EQ(gshgeom_map_check(&bad), GSHGEOM_OK);
  // Every defined kind bit is valid; and the entry at the bound is.
  bad = m;
  bad.m[0][0] = 0x1p100;
  bad.t[1] = -0x1p100;
  EXPECT_EQ(gshgeom_map_check(&bad), GSHGEOM_OK);
}

// ---------------------------------------------------------------------------
// Composition
// ---------------------------------------------------------------------------

TEST(MapCompose, TheInnerMapActsFirst) {
  GSHGEOM_Map a = Translate(10, 0);
  GSHGEOM_Map b = Scale(2, 2);
  GSHGEOM_Point r = Apply(Compose(a, b), P(1, 0));
  EXPECT_EQ(r.x, 12.0) << "scale, then translate";
  r = Apply(Compose(b, a), P(1, 0));
  EXPECT_EQ(r.x, 22.0) << "translate, then scale";
  // Composing with the identity on either side changes nothing.
  GSHGEOM_Map w = WildMap(3);
  EXPECT_TRUE(SameMap(Compose(Identity(3), w), w));
  EXPECT_TRUE(SameMap(Compose(w, Identity(3)), w));
}

TEST(MapCompose, ChainIsRootToLeafThenTheTargetLast) {
  GSHGEOM_Map chain[3] = {Translate(100, 0), Scale(2, 2), Rotate(GSHGEOM_ANGLE_DEGREES, 90)};
  GSHGEOM_Map target = Translate(0, 5);
  GSHGEOM_Map out;
  ASSERT_EQ(gshgeom_map_compose_chain(chain, 3, &target, &out), GSHGEOM_OK);
  // The leaf turns (1, 0) to (0, 1), the middle doubles it, the root shifts
  // x by 100, the target shifts y by 5.
  GSHGEOM_Point r = Apply(out, P(1, 0));
  EXPECT_TRUE(SameBits(r.x, 100.0));
  EXPECT_TRUE(SameBits(r.y, 7.0));
  EXPECT_EQ(out.kind, (unsigned)(TR | SC | QT));
  // No target: the chain alone.
  ASSERT_EQ(gshgeom_map_compose_chain(chain, 3, nullptr, &out), GSHGEOM_OK);
  r = Apply(out, P(1, 0));
  EXPECT_TRUE(SameBits(r.y, 2.0));
  // One map; zero maps and a target; a chain aliasing its output.
  ASSERT_EQ(gshgeom_map_compose_chain(chain, 1, nullptr, &out), GSHGEOM_OK);
  EXPECT_TRUE(SameMap(out, chain[0]));
  ASSERT_EQ(gshgeom_map_compose_chain(nullptr, 0, &target, &out), GSHGEOM_OK);
  EXPECT_TRUE(SameMap(out, target));
  ASSERT_EQ(gshgeom_map_compose_chain(chain, 3, &target, &chain[1]), GSHGEOM_OK);
}

TEST(MapCompose, TheChainIsFormedInOneStatedOrder) {
  // Matrix products are associative on paper and not in a double; the chain
  // is root first, left to right, the target last, and every other order is a
  // different matrix for some inputs.
  int differs_right = 0, differs_target_first = 0;
  for (int trial = 0; trial < 300; trial++) {
    GSHGEOM_Map c[4] = {WildMap(3), WildMap(3), WildMap(3), WildMap(3)};
    GSHGEOM_Map target = WildMap(3);
    GSHGEOM_Map acc = c[0];
    for (int i = 1; i < 4; i++) {
      acc = Compose(acc, c[i]);
    }
    GSHGEOM_Map ref = Compose(target, acc);
    GSHGEOM_Map got;
    ASSERT_EQ(gshgeom_map_compose_chain(c, 4, &target, &got), GSHGEOM_OK);
    EXPECT_TRUE(SameMap(got, ref)) << trial;
    // From the leaf up.
    GSHGEOM_Map right = Compose(c[2], c[3]);
    right = Compose(c[1], right);
    right = Compose(c[0], right);
    right = Compose(target, right);
    if (!SameMap(right, ref)) differs_right++;
    // The target folded in first.
    GSHGEOM_Map first = Compose(target, c[0]);
    for (int i = 1; i < 4; i++) {
      first = Compose(first, c[i]);
    }
    if (!SameMap(first, ref)) differs_target_first++;
  }
  EXPECT_GT(differs_right, 100) << "leaf-up order is a different matrix, so the test can tell";
  EXPECT_GT(differs_target_first, 100) << "target-first order is a different matrix";
}

TEST(MapCompose, EachEntryIsFormedInOrderWithoutFusing) {
  int differs = 0;
  for (int trial = 0; trial < 500; trial++) {
    GSHGEOM_Map a = WildMap(3);
    GSHGEOM_Map b = WildMap(3);
    GSHGEOM_Map r = Compose(a, b);
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        volatile double p0 = a.m[i][0] * b.m[0][j];
        volatile double p1 = a.m[i][1] * b.m[1][j];
        volatile double p2 = a.m[i][2] * b.m[2][j];
        volatile double s = p0 + p1;
        double want = s + p2;
        EXPECT_TRUE(SameBits(r.m[i][j], want == 0.0 ? 0.0 : want));
        volatile double other = p1 + p2;
        double alt = p0 + other;
        if (bits(alt) != bits(want)) differs++;
      }
      volatile double q0 = a.m[i][0] * b.t[0];
      volatile double q1 = a.m[i][1] * b.t[1];
      volatile double q2 = a.m[i][2] * b.t[2];
      volatile double qs = q0 + q1;
      volatile double qt = qs + q2;
      double want = qt + a.t[i];
      EXPECT_TRUE(SameBits(r.t[i], want == 0.0 ? 0.0 : want));
    }
  }
  EXPECT_GT(differs, 100) << "another association gives other bits, so the test can tell";
  // In 2D the third term is not formed: a 2D entry is the two-term sum.
  for (int trial = 0; trial < 200; trial++) {
    GSHGEOM_Map a = WildMap(2);
    GSHGEOM_Map b = WildMap(2);
    GSHGEOM_Map r = Compose(a, b);
    volatile double p0 = a.m[0][0] * b.m[0][1];
    volatile double p1 = a.m[0][1] * b.m[1][1];
    EXPECT_TRUE(SameBits(r.m[0][1], p0 + p1));
    EXPECT_TRUE(SameBits(r.m[2][2], 1.0));
    EXPECT_TRUE(SameBits(r.t[2], 0.0));
  }
}

TEST(MapCompose, TheTagIsTheUnionOfTheFactorsAndSurvives) {
  std::vector<GSHGEOM_Map> two = {
      Identity(2), Translate(1, 2), Scale(2, 2), Scale(2, 3), Scale(1, -1),
      Scale(-1, 1), Scale(-1, -1), Rotate(GSHGEOM_ANGLE_DEGREES, 90),
      Rotate(GSHGEOM_ANGLE_DEGREES, 30), Matrix(2, {1, 0.5, 0, 1}, {0, 0}),
      Scale(0, 1), Matrix(2, {0, 1, 1, 0}, {0, 0}), Matrix(2, {0.6, 0.8, 0.8, -0.6}, {0, 0})};
  std::vector<GSHGEOM_Map> three = {
      Identity(3), Translate(1, 2, 3, 3), Scale(2, 2, 2, 3), Scale(1, -1, 1, 3),
      Matrix(3, {0, 1, 0, 0, 0, 1, 1, 0, 0}, {0, 0, 0}),
      Matrix(3, {1, 2, 3, 0, 1, 4, 5, 6, 0}, {0, 0, 0})};
  for (const auto * set : {&two, &three}) {
    for (const GSHGEOM_Map & a : *set) {
      for (const GSHGEOM_Map & b : *set) {
        GSHGEOM_Map r;
        GSHGEOM_Result res = gshgeom_map_compose(&a, &b, &r);
        ASSERT_EQ(res, GSHGEOM_OK);
        EXPECT_EQ(r.kind, a.kind | b.kind) << a.kind << " . " << b.kind;
      }
    }
  }
}

TEST(MapCompose, ATagIsNotLostWhenTheNetMapLooksPlain) {
  // A quarter turn undone by its inverse is the identity, and still says so.
  GSHGEOM_Map r = Compose(Rotate(GSHGEOM_ANGLE_DEGREES, -90), Rotate(GSHGEOM_ANGLE_DEGREES, 90));
  EXPECT_TRUE(SameBits(r.m[0][0], 1.0));
  EXPECT_TRUE(SameBits(r.m[0][1], 0.0));
  EXPECT_EQ(r.kind, (unsigned)QT);
  EXPECT_FALSE(Aligned(r)) << "the route predicate reads the tag";
  // A flip undone by a flip.
  GSHGEOM_Map f = Scale(1, -1);
  r = Compose(f, f);
  EXPECT_EQ(r.kind, (unsigned)FY);
  EXPECT_TRUE(SameBits(r.m[1][1], 1.0));
}

TEST(MapCompose, NoEntryIsANegativeZero) {
  // (-1)(0) + (-1)(0) is -0 + -0, which is -0: the entries are made +0.
  GSHGEOM_Map outer = Matrix(2, {-1, -1, 2, 3}, {0, 0});
  GSHGEOM_Map inner = Matrix(2, {1, 0, 1, 0}, {0, 0});
  GSHGEOM_Map r = Compose(outer, inner);
  EXPECT_TRUE(SameBits(r.m[0][1], 0.0));
  EXPECT_TRUE(SameBits(r.m[1][1], 0.0));
  // The same for the translation, with a hand-filled -0 in the outer map.
  outer = Matrix(2, {-1, -1, 2, 3}, {0, 0});
  outer.t[0] = -0.0;
  r = Compose(outer, Identity(2));
  EXPECT_TRUE(SameBits(r.t[0], 0.0));
  GSHGEOM_Map m3 = Matrix(3, {-1, -1, -1, 0, 1, 0, 0, 0, 1}, {0, 0, 0});
  GSHGEOM_Map zero3 = Matrix(3, {1, 0, 0, 1, 0, 0, 1, 0, 0}, {0, 0, 0});
  r = Compose(m3, zero3);
  EXPECT_TRUE(SameBits(r.m[0][1], 0.0)) << "three terms, all -0";
  EXPECT_TRUE(SameBits(r.m[0][2], 0.0));
}

TEST(MapCompose, AProductThatCollapsesIsSingular) {
  GSHGEOM_Map tiny = Scale(0x1p-600, 0x1p-600);
  EXPECT_EQ(tiny.kind, (unsigned)SC);
  GSHGEOM_Map r = Compose(tiny, tiny);
  EXPECT_EQ(r.m[0][0], 0.0) << "2^-1200 underflows";
  EXPECT_EQ(r.kind, (unsigned)(SC | SG));
}

TEST(MapCompose, Refusals) {
  GSHGEOM_Map m2 = Translate(1, 2);
  GSHGEOM_Map m3 = Identity(3);
  GSHGEOM_Map sentinel = Scale(7, 7);
  GSHGEOM_Map out = sentinel;
  EXPECT_EQ(gshgeom_map_compose(&m2, &m3, &out), GSHGEOM_ERR_INVALID) << "dimensions differ";
  EXPECT_EQ(gshgeom_map_compose(nullptr, &m2, &out), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_compose(&m2, nullptr, &out), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_compose(&m2, &m2, nullptr), GSHGEOM_ERR_INVALID);
  GSHGEOM_Map nan = m2;
  nan.m[0][0] = kNaN;
  EXPECT_EQ(gshgeom_map_compose(&nan, &m2, &out), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_compose(&m2, &nan, &out), GSHGEOM_ERR_INVALID);
  GSHGEOM_Map bigscale = Scale(0x1p60, 0x1p60);
  EXPECT_EQ(gshgeom_map_compose(&bigscale, &bigscale, &out), GSHGEOM_ERR_LIMIT) << "2^120 is past the bound";
  GSHGEOM_Map bigt = Translate(0x1p99, 0);
  GSHGEOM_Map scale4 = Scale(4, 4);
  EXPECT_EQ(gshgeom_map_compose(&scale4, &bigt, &out), GSHGEOM_ERR_LIMIT) << "a translation of 2^101";
  EXPECT_TRUE(SameMap(out, sentinel)) << "untouched by a refusal";
  GSHGEOM_Map bound = Translate(0x1p99, 0);
  EXPECT_EQ(gshgeom_map_compose(&bound, &bound, &out), GSHGEOM_OK) << "2^100 is allowed";
  EXPECT_EQ(out.t[0], 0x1p100);
  // The chain.
  out = sentinel;
  GSHGEOM_Map mixed[2] = {m2, m3};
  EXPECT_EQ(gshgeom_map_compose_chain(mixed, 2, nullptr, &out), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_compose_chain(&m2, 1, &m3, &out), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_compose_chain(nullptr, 0, nullptr, &out), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_compose_chain(nullptr, 1, nullptr, &out), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_compose_chain(&m2, 1, nullptr, nullptr), GSHGEOM_ERR_INVALID);
  GSHGEOM_Map badone[2] = {m2, nan};
  EXPECT_EQ(gshgeom_map_compose_chain(badone, 2, nullptr, &out), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_compose_chain(&nan, 1, nullptr, &out), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_compose_chain(nullptr, 0, &nan, &out), GSHGEOM_ERR_INVALID);
  GSHGEOM_Map grow[3] = {bigscale, bigscale, bigscale};
  EXPECT_EQ(gshgeom_map_compose_chain(grow, 3, nullptr, &out), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(gshgeom_map_compose_chain(&m2, 1, &grow[0], &out), GSHGEOM_OK);
  GSHGEOM_Map two[2] = {bigscale, scale4};
  EXPECT_EQ(gshgeom_map_compose_chain(two, 2, &bigscale, &out), GSHGEOM_ERR_LIMIT) << "the target's turn comes last";
}

// ---------------------------------------------------------------------------
// Application
// ---------------------------------------------------------------------------

TEST(MapApply, PointsAreFormedInOrder) {
  GSHGEOM_Map m = WildMap(3);
  for (int i = 0; i < 300; i++) {
    GSHGEOM_Point p = P(wild(), wild(), wild());
    GSHGEOM_Point r = Apply(m, p);
    for (int row = 0; row < 3; row++) {
      volatile double a = m.m[row][0] * p.x;
      volatile double b = m.m[row][1] * p.y;
      volatile double c = m.m[row][2] * p.z;
      volatile double s = a + b;
      volatile double u = s + c;
      double want = u + m.t[row];
      double got = row == 0 ? r.x : row == 1 ? r.y : r.z;
      EXPECT_TRUE(SameBits(got, want));
    }
    GSHGEOM_Point v;
    ASSERT_EQ(gshgeom_map_apply_vector(&m, p, &v), GSHGEOM_OK);
    volatile double a = m.m[0][0] * p.x;
    volatile double b = m.m[0][1] * p.y;
    volatile double c = m.m[0][2] * p.z;
    volatile double s = a + b;
    EXPECT_TRUE(SameBits(v.x, s + c)) << "a vector is not translated";
  }
}

TEST(MapApply, TwoDimensionalMapsIgnoreZAndGiveZero) {
  GSHGEOM_Map m = Compose(Translate(1, 2), Scale(2, 3));
  GSHGEOM_Point r = Apply(m, P(5, 7, 99));
  EXPECT_EQ(r.x, 11.0);
  EXPECT_EQ(r.y, 23.0);
  EXPECT_EQ(r.z, 0.0);
  GSHGEOM_Point v;
  ASSERT_EQ(gshgeom_map_apply_vector(&m, P(5, 7, 99), &v), GSHGEOM_OK);
  EXPECT_EQ(v.x, 10.0);
  EXPECT_EQ(v.y, 21.0);
  EXPECT_EQ(v.z, 0.0);
  // A NaN z in 2D is not read.
  GSHGEOM_Point q;
  EXPECT_EQ(gshgeom_map_apply_point(&m, P(5, 7, kNaN), &q), GSHGEOM_OK);
  GSHGEOM_Map m3 = Translate(1, 2, 3, 3);
  EXPECT_EQ(gshgeom_map_apply_point(&m3, P(5, 7, kNaN), &q), GSHGEOM_ERR_INVALID);
  r = Apply(m3, P(5, 7, 9));
  EXPECT_EQ(r.z, 12.0);
}

TEST(MapApply, Refusals) {
  GSHGEOM_Map m = Scale(2, 2);
  GSHGEOM_Point sentinel = P(-1, -2, -3);
  GSHGEOM_Point out = sentinel;
  EXPECT_EQ(gshgeom_map_apply_point(nullptr, P(1, 1), &out), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_apply_point(&m, P(1, 1), nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_apply_point(&m, P(kNaN, 1), &out), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_apply_point(&m, P(1, kInf), &out), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_apply_vector(&m, P(1, -kInf), &out), GSHGEOM_ERR_INVALID);
  GSHGEOM_Map bad = m;
  bad.m[0][0] = kNaN;
  EXPECT_EQ(gshgeom_map_apply_point(&bad, P(1, 1), &out), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_apply_vector(&bad, P(1, 1), &out), GSHGEOM_ERR_INVALID);
  // A finite point whose image is not finite.
  GSHGEOM_Map big = Scale(0x1p100, 0x1p100);
  EXPECT_EQ(gshgeom_map_apply_point(&big, P(1e300, 0), &out), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(gshgeom_map_apply_vector(&big, P(0, -1e300), &out), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(out.x, sentinel.x);
  EXPECT_EQ(out.y, sentinel.y);
  EXPECT_EQ(out.z, sentinel.z);
}

TEST(MapApply, ACurveMovesItsPointsAndKeepsItsWeights) {
  GSHGEOM_Curve c, r;
  const double w[4] = {1.0, 2.0, 0.5, 4.0};
  ASSERT_EQ(gshgeom_curve_rational_cubic(2, P(0, 0), P(1, 2), P(3, 2), P(4, 0), w, &c), GSHGEOM_OK);
  GSHGEOM_Map m = Compose(Translate(10, 20), Scale(2, 3));
  ASSERT_EQ(gshgeom_map_apply_curve(&m, &c, &r), GSHGEOM_OK);
  EXPECT_EQ(r.kind, GSHGEOM_CURVE_RATIONAL_CUBIC);
  EXPECT_EQ(r.dim, 2);
  for (int i = 0; i < 4; i++) {
    EXPECT_EQ(r.w[i], w[i]);
    GSHGEOM_Point want = Apply(m, c.p[i]);
    EXPECT_EQ(r.p[i].x, want.x);
    EXPECT_EQ(r.p[i].y, want.y);
  }
  // A point of the image is the image of the point (affine invariance).
  for (double t : {0.0, 0.25, 0.5, 0.9, 1.0}) {
    GSHGEOM_Point on, on_r;
    ASSERT_EQ(gshgeom_curve_eval(&c, t, &on), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_curve_eval(&r, t, &on_r), GSHGEOM_OK);
    GSHGEOM_Point image = Apply(m, on);
    EXPECT_NEAR(on_r.x, image.x, 1e-12);
    EXPECT_NEAR(on_r.y, image.y, 1e-12);
  }
  // A quarter turn takes the box of a curve to the turned box.
  GSHGEOM_Map q = Rotate(GSHGEOM_ANGLE_DEGREES, 90);
  GSHGEOM_Bounds b0, b1;
  ASSERT_EQ(gshgeom_curve_bounds(&c, &b0), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_map_apply_curve(&q, &c, &r), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_curve_bounds(&r, &b1), GSHGEOM_OK);
  EXPECT_NEAR(b1.min.x, -b0.max.y, 1e-12);
  EXPECT_NEAR(b1.max.x, -b0.min.y, 1e-12);
  EXPECT_NEAR(b1.min.y, b0.min.x, 1e-12);
  EXPECT_NEAR(b1.max.y, b0.max.x, 1e-12);
  // 3D.
  GSHGEOM_Curve c3, r3;
  ASSERT_EQ(gshgeom_curve_quadratic(3, P(0, 0, 0), P(1, 2, 3), P(2, 0, 1), &c3), GSHGEOM_OK);
  GSHGEOM_Map m3 = Translate(1, 1, 1, 3);
  ASSERT_EQ(gshgeom_map_apply_curve(&m3, &c3, &r3), GSHGEOM_OK);
  EXPECT_EQ(r3.p[1].z, 4.0);
  EXPECT_EQ(r3.p[2].x, 3.0);
}

TEST(MapApply, ACurveRefuses) {
  GSHGEOM_Curve c, r;
  ASSERT_EQ(gshgeom_curve_line(2, P(0, 0), P(0x1p99, 0), &c), GSHGEOM_OK);
  GSHGEOM_Map sentinel = Scale(1, 1);
  r = c;
  GSHGEOM_Curve before = r;
  GSHGEOM_Map m3 = Identity(3);
  EXPECT_EQ(gshgeom_map_apply_curve(&m3, &c, &r), GSHGEOM_ERR_INVALID) << "dimensions differ";
  GSHGEOM_Map big = Scale(4, 4);
  EXPECT_EQ(gshgeom_map_apply_curve(&big, &c, &r), GSHGEOM_ERR_LIMIT) << "the image passes the coordinate limit";
  EXPECT_EQ(gshgeom_map_apply_curve(nullptr, &c, &r), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_apply_curve(&sentinel, nullptr, &r), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_apply_curve(&sentinel, &c, nullptr), GSHGEOM_ERR_INVALID);
  GSHGEOM_Curve bad = c;
  bad.p[0].x = kNaN;
  EXPECT_EQ(gshgeom_map_apply_curve(&sentinel, &bad, &r), GSHGEOM_ERR_INVALID);
  GSHGEOM_Map badmap = sentinel;
  badmap.t[0] = kNaN;
  EXPECT_EQ(gshgeom_map_apply_curve(&badmap, &c, &r), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(std::memcmp(&r, &before, sizeof r), 0) << "untouched by a refusal";
}

// ---------------------------------------------------------------------------
// The route predicate
// ---------------------------------------------------------------------------

TEST(MapRoute, Accepts) {
  GSHGEOM_Grid g;
  std::memset(&g, 0, sizeof g);
  EXPECT_TRUE(Aligned(Identity(2), &g));
  EXPECT_EQ(g.scale, 1.0);
  EXPECT_EQ(g.flip_y, 0);
  EXPECT_TRUE(Aligned(Translate(3.5, -2), &g));
  EXPECT_EQ(g.tx, 3.5);
  EXPECT_EQ(g.ty, -2.0);
  EXPECT_TRUE(Aligned(Scale(2.5, 2.5), &g));
  EXPECT_EQ(g.scale, 2.5);
  EXPECT_EQ(g.flip_y, 0);
  GSHGEOM_Map d;
  ASSERT_EQ(gshgeom_map_design_to_model(2.5, 10, 20, &d), GSHGEOM_OK);
  EXPECT_TRUE(Aligned(d, &g));
  EXPECT_EQ(g.scale, 2.5);
  EXPECT_EQ(g.tx, 10.0);
  EXPECT_EQ(g.ty, 20.0);
  EXPECT_EQ(g.flip_y, 1);
  GSHGEOM_Map f;
  ASSERT_EQ(gshgeom_map_flip_y(2, &f), GSHGEOM_OK);
  EXPECT_TRUE(Aligned(f, &g));
  EXPECT_EQ(g.scale, 1.0);
  EXPECT_EQ(g.flip_y, 1);
  EXPECT_TRUE(Aligned(Matrix(2, {3, 0, 0, -3}, {1, 2}), &g));
  EXPECT_EQ(g.scale, 3.0);
  EXPECT_EQ(g.tx, 1.0);
  // Composed by the chain composer: the glyph node's map, the scene's, the target's.
  GSHGEOM_Map chain[3] = {Translate(5, 6), Scale(2, 2), f};
  GSHGEOM_Map target = Scale(1.5, 1.5);
  GSHGEOM_Map out;
  ASSERT_EQ(gshgeom_map_compose_chain(chain, 3, &target, &out), GSHGEOM_OK);
  EXPECT_TRUE(Aligned(out, &g));
  EXPECT_EQ(g.scale, 3.0);
  EXPECT_EQ(g.flip_y, 1);
  EXPECT_EQ(g.tx, 7.5);
  EXPECT_EQ(g.ty, 9.0);
  // A rotation by nothing is the identity and carries no tag.
  EXPECT_TRUE(Aligned(Compose(d, Rotate(GSHGEOM_ANGLE_DEGREES, 0))));
  EXPECT_TRUE(Aligned(Compose(d, Rotate(GSHGEOM_ANGLE_RADIANS, 0))));
  EXPECT_TRUE(Aligned(d, nullptr)) << "the grid is optional";
}

TEST(MapRoute, RejectsEachOtherClass) {
  const GSHGEOM_Map rejected[] = {
      Rotate(GSHGEOM_ANGLE_DEGREES, 30),
      Rotate(GSHGEOM_ANGLE_DEGREES, 90),
      Rotate(GSHGEOM_ANGLE_DEGREES, 180),
      Rotate(GSHGEOM_ANGLE_DEGREES, -90),
      Rotate(GSHGEOM_ANGLE_RADIANS, 1e-9),
      Matrix(2, {0, 1, 1, 0}, {0, 0}),
      Matrix(2, {1, 0.5, 0, 1}, {0, 0}),
      Matrix(2, {1, 0, 0.5, 1}, {0, 0}),
      Matrix(2, {0.6, 0.8, 0.8, -0.6}, {0, 0}),
      Scale(-1, 1),
      Scale(-2, 2),
      Scale(-2, -2),
      Scale(2, 3),
      Scale(2, -3),
      Scale(3, 2),
      Scale(0, 1),
      Scale(1, 0),
      Scale(0, 0),
      Matrix(2, {1, 0, 0, -1.0000000000000002}, {0, 0}),
      Matrix(2, {1.0000000000000002, 0, 0, 1}, {0, 0}),
      Identity(3),
      Translate(1, 2, 3, 3),
  };
  for (size_t i = 0; i < sizeof rejected / sizeof rejected[0]; i++) {
    SCOPED_TRACE(i);
    GSHGEOM_Grid g = {-9, -9, -9, -9};
    int a = -1;
    EXPECT_EQ(gshgeom_map_grid_aligned(&rejected[i], &a, &g), GSHGEOM_OK);
    EXPECT_EQ(a, 0);
    EXPECT_EQ(g.scale, -9.0) << "the grid is untouched when the map is rejected";
  }
  // A composition that mixes in a rejected class, even one that nets to nothing.
  GSHGEOM_Map d;
  ASSERT_EQ(gshgeom_map_design_to_model(2, 1, 1, &d), GSHGEOM_OK);
  EXPECT_FALSE(Aligned(Compose(d, Rotate(GSHGEOM_ANGLE_DEGREES, 90))));
  EXPECT_FALSE(Aligned(Compose(Rotate(GSHGEOM_ANGLE_DEGREES, 90), d)));
  EXPECT_FALSE(Aligned(Compose(d, Scale(1, 2))));
  EXPECT_FALSE(Aligned(Compose(d, Matrix(2, {1, 0.25, 0, 1}, {0, 0}))));
  EXPECT_FALSE(Aligned(Compose(d, Scale(-1, 1))));
}

TEST(MapRoute, ReadsTheTagAndTheValuesSeparately) {
  GSHGEOM_Map plain = Scale(2, -2);
  ASSERT_TRUE(Aligned(plain));
  // A tag that says rotation over values that are plain.
  GSHGEOM_Map stale = plain;
  stale.kind |= GSHGEOM_MAP_ROTATES;
  EXPECT_FALSE(Aligned(stale));
  for (unsigned bit : {NU, MI, QT, RO, SK, SG, GE}) {
    GSHGEOM_Map m = plain;
    m.kind |= bit;
    EXPECT_FALSE(Aligned(m)) << bit;
  }
  for (unsigned bit : {TR, SC, FY}) {
    GSHGEOM_Map m = plain;
    m.kind |= bit;
    EXPECT_TRUE(Aligned(m)) << bit;
  }
  // A zero tag over values that are not.
  const GSHGEOM_Map lying[] = {
      Rotate(GSHGEOM_ANGLE_DEGREES, 30), Rotate(GSHGEOM_ANGLE_DEGREES, 90),
      Scale(2, 3), Scale(-1, 1), Scale(2, -3), Scale(0, 1), Scale(-2, -2),
      Matrix(2, {1, 0.5, 0, 1}, {0, 0}),
      Matrix(2, {1, 0, 0.5, 1}, {0, 0}),
      Scale(0, 0),
      Matrix(2, {1, 0, 0, -1.0000000000000002}, {0, 0}),
      // Squares that underflow to the same zero: only the values tell them apart.
      Matrix(2, {0x1p-600, 0, 0, -(0x1p-600 + 0x1p-651)}, {0, 0}),
  };
  for (size_t i = 0; i < sizeof lying / sizeof lying[0]; i++) {
    GSHGEOM_Map m = lying[i];
    m.kind = 0;
    EXPECT_FALSE(Aligned(m)) << i;
    m.kind = TR | SC | FY;
    EXPECT_FALSE(Aligned(m)) << i;
  }
  // A map that is only a translation tagged as nothing is still aligned: the
  // values decide what the tag allows.
  GSHGEOM_Map t = Translate(1, 2);
  t.kind = 0;
  EXPECT_TRUE(Aligned(t));
  // A negative scale: the x scale must be positive.
  GSHGEOM_Map neg = Scale(-2, -2);
  neg.kind = 0;
  EXPECT_FALSE(Aligned(neg));
}

TEST(MapRoute, Refusals) {
  int a = 7;
  GSHGEOM_Map m = Identity(2);
  EXPECT_EQ(gshgeom_map_grid_aligned(nullptr, &a, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_grid_aligned(&m, nullptr, nullptr), GSHGEOM_ERR_INVALID);
  GSHGEOM_Map bad = m;
  bad.m[0][0] = kNaN;
  EXPECT_EQ(gshgeom_map_grid_aligned(&bad, &a, nullptr), GSHGEOM_ERR_INVALID);
  bad = m;
  bad.kind = 1u << 20;
  EXPECT_EQ(gshgeom_map_grid_aligned(&bad, &a, nullptr), GSHGEOM_ERR_INVALID);
  bad = m;
  bad.t[1] = 0x1p101;
  EXPECT_EQ(gshgeom_map_grid_aligned(&bad, &a, nullptr), GSHGEOM_ERR_LIMIT);
  bad = m;
  bad.dim = 5;
  EXPECT_EQ(gshgeom_map_grid_aligned(&bad, &a, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(a, 7) << "untouched by a refusal";
}

// ---------------------------------------------------------------------------
// Hash and dump
// ---------------------------------------------------------------------------

TEST(MapHash, KnownAnswersFromTheGenerator) {
  size_t n = sizeof kMapFixtures / sizeof kMapFixtures[0];
  for (size_t i = 0; i < n; i++) {
    const MapFixture & f = kMapFixtures[i];
    SCOPED_TRACE(f.name);
    GSHGEOM_Map m;
    ASSERT_EQ(gshgeom_map_from_matrix(f.dim, f.m, f.t, &m), GSHGEOM_OK);
    EXPECT_EQ(gshgeom_map_hash(&m), f.hash);
    EXPECT_EQ(DumpOf(m), std::string(f.dump));
  }
}

TEST(MapHash, EveryFieldIsHashed) {
  GSHGEOM_Map base = Matrix(3, {1, 2, 3, 0, 1, 4, 5, 6, 0}, {7, 8, 9});
  uint64_t h0 = gshgeom_map_hash(&base);
  std::vector<uint64_t> seen = {h0};
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      GSHGEOM_Map m = base;
      m.m[i][j] += 0.5;
      uint64_t h = gshgeom_map_hash(&m);
      for (uint64_t s : seen) EXPECT_NE(h, s) << i << j;
      seen.push_back(h);
    }
    GSHGEOM_Map m = base;
    m.t[i] += 0.5;
    uint64_t h = gshgeom_map_hash(&m);
    for (uint64_t s : seen) EXPECT_NE(h, s) << "t" << i;
    seen.push_back(h);
  }
  GSHGEOM_Map k = base;
  k.kind ^= GSHGEOM_MAP_SKEWS;
  EXPECT_NE(gshgeom_map_hash(&k), h0) << "the tag is hashed";
  GSHGEOM_Map z = base;
  z.m[0][0] = -0.0;
  GSHGEOM_Map p = base;
  p.m[0][0] = 0.0;
  EXPECT_NE(gshgeom_map_hash(&z), gshgeom_map_hash(&p)) << "-0.0 is not 0.0";
}

TEST(MapHash, ATwoDimensionalMapHashesOnlyItsTwoDimensions) {
  GSHGEOM_Map a = Compose(Translate(1, 2), Scale(2, 3));
  GSHGEOM_Map b = a;
  b.m[2][2] = 7;
  b.m[0][2] = 5;
  b.t[2] = 9;
  EXPECT_EQ(gshgeom_map_hash(&a), gshgeom_map_hash(&b));
  GSHGEOM_Map three = Identity(3);
  EXPECT_NE(gshgeom_map_hash(&a), gshgeom_map_hash(&three));
}

TEST(MapHash, TheDesignMapEqualsTheMatrixItIs) {
  GSHGEOM_Map d;
  ASSERT_EQ(gshgeom_map_design_to_model(2.5, 10, 20, &d), GSHGEOM_OK);
  GSHGEOM_Map m = Matrix(2, {2.5, 0, 0, -2.5}, {10, 20});
  EXPECT_EQ(gshgeom_map_hash(&d), gshgeom_map_hash(&m));
}

TEST(MapHash, AnInvalidMapHashesAsTheEmptyStream) {
  GSHGEOM_Hash h;
  gshgeom_hash_init(&h);
  uint64_t empty = gshgeom_hash_value(&h);
  EXPECT_EQ(gshgeom_map_hash(nullptr), empty);
  GSHGEOM_Map bad = Identity(2);
  bad.m[0][0] = kNaN;
  EXPECT_EQ(gshgeom_map_hash(&bad), empty);
  GSHGEOM_Map good = Identity(2);
  EXPECT_NE(gshgeom_map_hash(&good), empty);
}

TEST(MapDump, FailureAndRefusal) {
  GSHGEOM_Map m = Compose(Translate(1, 2), Scale(2, 3));
  FILE * full = std::fopen("/dev/full", "w");
  if (full) {
    EXPECT_EQ(gshgeom_map_dump(&m, full), GSHGEOM_ERR_IO);
    std::fclose(full);
  }
  FILE * ro = std::fopen("/dev/null", "r");
  ASSERT_NE(ro, nullptr);
  EXPECT_EQ(gshgeom_map_dump(&m, ro), GSHGEOM_ERR_IO);
  std::fclose(ro);
  EXPECT_EQ(gshgeom_map_dump(&m, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_map_dump(nullptr, stdout), GSHGEOM_ERR_INVALID);
  GSHGEOM_Map bad = m;
  bad.m[1][1] = kInf;
  FILE * f = std::tmpfile();
  ASSERT_NE(f, nullptr);
  EXPECT_EQ(gshgeom_map_dump(&bad, f), GSHGEOM_ERR_INVALID);
  std::fclose(f);
  // Zero and its sign are 16 digits too.
  GSHGEOM_Map z = Identity(2);
  z.t[1] = -0.0;
  std::string d = DumpOf(z);
  EXPECT_NE(d.find("\nt 0000000000000000 8000000000000000\n"), std::string::npos) << d;
}

TEST(MapDump, EveryWriteIsChecked) {
  GSHGEOM_Map m2 = Compose(Translate(1, 2), Scale(2, 3));
  GSHGEOM_Map m3 = Identity(3);
  GSHGEOM_EVERY_WRITE_CHECKED(gshgeom_map_dump, &m2, "2D map");
  GSHGEOM_EVERY_WRITE_CHECKED(gshgeom_map_dump, &m3, "3D map");
}

// ---------------------------------------------------------------------------
// The chiral fixture: the orientation rules, pinned
// ---------------------------------------------------------------------------

namespace {

// A letter F in y-down coordinates.  It has no symmetry, so a mirror image is
// not a rotation of it; walked as listed it goes +x, then +y, then -x: the
// positive sense.
const GSHGEOM_Point kF[] = {P(0, 0), P(3, 0), P(3, 1), P(1, 1), P(1, 2),
    P(2, 2), P(2, 3), P(1, 3), P(1, 5), P(0, 5)};
const size_t kFn = sizeof kF / sizeof kF[0];

double AreaOf(const std::vector<GSHGEOM_Point> & pts) {
  double a = 0;
  EXPECT_EQ(gshgeom_polygon_signed_area(pts.data(), pts.size(), &a), GSHGEOM_OK);
  return a;
}

std::vector<GSHGEOM_Point> Image(const GSHGEOM_Map & m) {
  std::vector<GSHGEOM_Point> out;
  for (size_t i = 0; i < kFn; i++) {
    out.push_back(Apply(m, kF[i]));
  }
  return out;
}

int Orient(GSHGEOM_Point a, GSHGEOM_Point b, GSHGEOM_Point c) {
  int s = 99;
  EXPECT_EQ(gshgeom_orient2d(a, b, c, &s), GSHGEOM_OK);
  return s;
}

}  // namespace

TEST(Chiral, TheFixtureWindsPositively) {
  std::vector<GSHGEOM_Point> f(kF, kF + kFn);
  EXPECT_EQ(AreaOf(f), 8.0);
  // Its first corner, (0,0) -> (3,0) -> (3,1), turns from +x toward +y.
  EXPECT_EQ(Orient(kF[0], kF[1], kF[2]), 1);
  // The notch corner (1,1) -> (1,2) -> (2,2) turns the other way.
  EXPECT_EQ(Orient(kF[3], kF[4], kF[5]), -1);
}

TEST(Chiral, ARotationKeepsTheSenseAndAMirrorReversesIt) {
  GSHGEOM_Map q = Rotate(GSHGEOM_ANGLE_DEGREES, 90);
  std::vector<GSHGEOM_Point> turned = Image(q);
  // The committed image: (x, y) -> (-y, x).
  const double want[10][2] = {{0, 0}, {0, 3}, {-1, 3}, {-1, 1}, {-2, 1},
      {-2, 2}, {-3, 2}, {-3, 1}, {-5, 1}, {-5, 0}};
  for (size_t i = 0; i < kFn; i++) {
    EXPECT_TRUE(SameBits(turned[i].x, want[i][0])) << i;
    EXPECT_TRUE(SameBits(turned[i].y, want[i][1])) << i;
  }
  EXPECT_EQ(AreaOf(turned), 8.0);
  EXPECT_EQ(Orient(turned[0], turned[1], turned[2]), 1);
  // Every quarter turn and every other rotation keeps +8 ...
  for (double deg : {180.0, 270.0, 33.0, -75.0}) {
    EXPECT_NEAR(AreaOf(Image(Rotate(GSHGEOM_ANGLE_DEGREES, deg))), 8.0, 1e-13) << deg;
  }
  // ... and the two mirrors reverse it, so no rotation takes F to its mirror.
  std::vector<GSHGEOM_Point> flipped = Image(Scale(1, -1));
  EXPECT_EQ(AreaOf(flipped), -8.0);
  EXPECT_EQ(Orient(flipped[0], flipped[1], flipped[2]), -1);
  std::vector<GSHGEOM_Point> mirrored = Image(Scale(-1, 1));
  EXPECT_EQ(AreaOf(mirrored), -8.0);
  EXPECT_EQ(Orient(mirrored[0], mirrored[1], mirrored[2]), -1);
  // Two mirrors are a rotation: a half turn.
  EXPECT_EQ(AreaOf(Image(Scale(-1, -1))), 8.0);
  // The canonical design-to-model map flips the sense (a y-up glyph walked
  // counterclockwise in design space is clockwise on a y-down screen).
  GSHGEOM_Map d;
  ASSERT_EQ(gshgeom_map_design_to_model(1.0, 0, 0, &d), GSHGEOM_OK);
  EXPECT_EQ(AreaOf(Image(d)), -8.0);
  // The inverse sense listed backwards.
  std::vector<GSHGEOM_Point> rev(kF, kF + kFn);
  std::vector<GSHGEOM_Point> r(rev.rbegin(), rev.rend());
  EXPECT_EQ(AreaOf(r), -8.0);
}

TEST(Chiral, TheLeftOfATangentIsATurnByANegativeAngle) {
  // Walking +x (to the right on the screen), the left hand is up the screen.
  GSHGEOM_Point t = P(1, 0);
  GSHGEOM_Point left = gshgeom_normal_left(t);
  GSHGEOM_Point turned = Apply(Rotate(GSHGEOM_ANGLE_DEGREES, -90), t);
  EXPECT_TRUE(SameBits(left.x, turned.x));
  EXPECT_TRUE(SameBits(left.y, turned.y));
  GSHGEOM_Point right = gshgeom_normal_right(t);
  turned = Apply(Rotate(GSHGEOM_ANGLE_DEGREES, 90), t);
  EXPECT_TRUE(SameBits(right.x, turned.x));
  EXPECT_TRUE(SameBits(right.y, turned.y));
  // For a positively wound boundary the inside is on the right, so the left
  // normal points out: a positive offset moves along it and grows the region.
  EXPECT_EQ(Orient(kF[0], kF[1], P(kF[0].x + 1.5, kF[0].y + 0.5)), 1) << "inside is on the positive side";
  GSHGEOM_Point edge = P(kF[1].x - kF[0].x, kF[1].y - kF[0].y);
  GSHGEOM_Point out = gshgeom_normal_left(edge);
  EXPECT_LT(out.y, 0.0) << "outward from the top edge is up the screen";
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
