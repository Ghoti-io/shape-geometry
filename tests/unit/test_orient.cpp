/**
 * @file
 *
 * The exact orientation predicates and the orientation rules (CAP-2, CAP-4,
 * AD-6): vectors whose signs were decided by rational arithmetic in Python
 * (tools/gen-hash-vectors.py), on which the naive double expression is shown
 * to be wrong; the filter against the exact path; the exact engine itself; the
 * normals and the signed area; and the refusals.
 *
 * Copyright 2026 by Corey Pennycuff
 */

#include "test_helpers.h"

#include "../known_answers.h"

extern "C" {
#include "../../src/exact/exact.h"
}

#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

namespace {

const double kNaN = std::numeric_limits<double>::quiet_NaN();
const double kInf = std::numeric_limits<double>::infinity();

GSHGEOM_Point P(double x, double y, double z = 0.0) {
  return GSHGEOM_Point{x, y, z};
}

GSHGEOM_Point P2(const double * v) {
  return P(v[0], v[1]);
}

GSHGEOM_Point P3(const double * v) {
  return P(v[0], v[1], v[2]);
}

int sgn(double v) {
  return (v > 0) - (v < 0);
}

// The expressions a caller would write, formed in the order the generator
// formed them (the tests are built with -ffp-contract=off).
int naive2(const Orient2dCase & c) {
  return sgn((c.b[0] - c.a[0]) * (c.c[1] - c.a[1]) -
      (c.b[1] - c.a[1]) * (c.c[0] - c.a[0]));
}

int naive3(const Orient3dCase & c) {
  double ux = c.b[0] - c.a[0], uy = c.b[1] - c.a[1], uz = c.b[2] - c.a[2];
  double vx = c.c[0] - c.a[0], vy = c.c[1] - c.a[1], vz = c.c[2] - c.a[2];
  double wx = c.d[0] - c.a[0], wy = c.d[1] - c.a[1], wz = c.d[2] - c.a[2];
  return sgn(ux * (vy * wz - vz * wy) - uy * (vx * wz - vz * wx) +
      uz * (vx * wy - vy * wx));
}

int o2(GSHGEOM_Point a, GSHGEOM_Point b, GSHGEOM_Point c) {
  int s = 99;
  EXPECT_EQ(gshgeom_orient2d(a, b, c, &s), GSHGEOM_OK);
  return s;
}

int o2x(GSHGEOM_Point a, GSHGEOM_Point b, GSHGEOM_Point c) {
  int s = 99;
  EXPECT_EQ(gshgeom_orient2d_exact(a, b, c, &s), GSHGEOM_OK);
  return s;
}

int o3(GSHGEOM_Point a, GSHGEOM_Point b, GSHGEOM_Point c, GSHGEOM_Point d) {
  int s = 99;
  EXPECT_EQ(gshgeom_orient3d(a, b, c, d, &s), GSHGEOM_OK);
  return s;
}

int o3x(GSHGEOM_Point a, GSHGEOM_Point b, GSHGEOM_Point c, GSHGEOM_Point d) {
  int s = 99;
  EXPECT_EQ(gshgeom_orient3d_exact(a, b, c, d, &s), GSHGEOM_OK);
  return s;
}

uint64_t rng_state = UINT64_C(0x9E3779B97F4A7C15);

uint64_t next() {
  rng_state ^= rng_state >> 12;
  rng_state ^= rng_state << 25;
  rng_state ^= rng_state >> 27;
  return rng_state * UINT64_C(2685821657736338717);
}

// A coordinate on a 2^-20 grid in [-16, 16).
double coord() {
  return (double)((int64_t)(next() >> 40) - (int64_t)(1 << 23)) / 1048576.0;
}

}  // namespace

// ---------------------------------------------------------------------------
// The vectors: exact signs from rational arithmetic
// ---------------------------------------------------------------------------

TEST(Orient2d, ExactSignOnTheGeneratedVectors) {
  size_t n = sizeof kOrient2dCases / sizeof kOrient2dCases[0];
  ASSERT_GT(n, 500u);
  unsigned wrong = 0;
  int seen[3] = {0, 0, 0};
  for (size_t i = 0; i < n; i++) {
    const Orient2dCase & c = kOrient2dCases[i];
    EXPECT_EQ(o2(P2(c.a), P2(c.b), P2(c.c)), c.sign) << "case " << i;
    EXPECT_EQ(o2x(P2(c.a), P2(c.b), P2(c.c)), c.sign) << "exact, case " << i;
    seen[c.sign + 1]++;
    if (naive2(c) != c.sign) {
      wrong++;
    }
  }
  // The vectors discriminate: the naive expression is wrong on as many as the
  // generator counted, and there are many, and all three signs occur.
  EXPECT_EQ(wrong, kNaiveWrong2d);
  EXPECT_GT(wrong, 100u);
  EXPECT_GT(seen[0], 50);
  EXPECT_GT(seen[1], 20);
  EXPECT_GT(seen[2], 50);
}

TEST(Orient3d, ExactSignOnTheGeneratedVectors) {
  size_t n = sizeof kOrient3dCases / sizeof kOrient3dCases[0];
  ASSERT_GT(n, 400u);
  unsigned wrong = 0;
  int seen[3] = {0, 0, 0};
  for (size_t i = 0; i < n; i++) {
    const Orient3dCase & c = kOrient3dCases[i];
    EXPECT_EQ(o3(P3(c.a), P3(c.b), P3(c.c), P3(c.d)), c.sign) << "case " << i;
    EXPECT_EQ(o3x(P3(c.a), P3(c.b), P3(c.c), P3(c.d)), c.sign)
        << "exact, case " << i;
    seen[c.sign + 1]++;
    if (naive3(c) != c.sign) {
      wrong++;
    }
  }
  EXPECT_EQ(wrong, kNaiveWrong3d);
  EXPECT_GT(wrong, 100u);
  EXPECT_GT(seen[0], 30);
  EXPECT_GT(seen[1], 10);
  EXPECT_GT(seen[2], 30);
}

TEST(Orient2d, ADoubleCannotSeeTheSmallestTriangle) {
  // The legs are 2^-1074; their product, 2^-2148, is zero in a double.
  const double h = std::ldexp(1.0, -1074);
  EXPECT_EQ(0.0, h * h);
  EXPECT_EQ(o2(P(0, 0), P(h, 0), P(0, h)), 1);
  EXPECT_EQ(o2(P(0, 0), P(h, 0), P(0, -h)), -1);
  EXPECT_EQ(o2(P(0, 0), P(0, h), P(h, 0)), -1);
  EXPECT_EQ(o2(P(0, 0), P(h, h), P(2 * h, 2 * h)), 0);
}

TEST(Orient2d, TheLargestCoordinates) {
  const double big = std::numeric_limits<double>::max();
  // (b-a) overflows in a naive double expression; exactly, the determinant of
  // (big, 0), (0, big), (-big, 0) is 2 big^2, positive, and swapping the ends
  // negates it.
  EXPECT_EQ(o2(P(big, 0), P(0, big), P(-big, 0)), 1);
  EXPECT_EQ(o2(P(-big, 0), P(0, big), P(big, 0)), -1);
  EXPECT_EQ(o2(P(big, big), P(big, big), P(big, big)), 0);
  EXPECT_EQ(o2(P(1e300, 1e300), P(2e300, 2e300), P(3e300, 3e300)), 0);
  EXPECT_EQ(o2x(P(1e300, 1e300), P(2e300, 2e300), P(3e300, 3e300)), 0);
}

TEST(Orient2d, ExactZeros) {
  EXPECT_EQ(o2(P(0, 0), P(1, 0), P(2, 0)), 0);
  EXPECT_EQ(o2(P(0, 0), P(0, 0), P(5, 7)), 0);
  EXPECT_EQ(o2(P(1, 1), P(1, 1), P(1, 1)), 0);
  EXPECT_EQ(o2(P(-0.0, 0), P(0, -0.0), P(0, 0)), 0);
  EXPECT_EQ(o2(P(0.1, 0.1), P(0.3, 0.3), P(0.7, 0.7)), 0);
}

TEST(Orient2d, SenseAndSymmetries) {
  // c at +y of the line +x is on the positive side (the rules in orient.h).
  EXPECT_EQ(o2(P(0, 0), P(1, 0), P(0, 1)), 1);
  EXPECT_EQ(o2(P(0, 0), P(1, 0), P(0, -1)), -1);
  // Antisymmetric in any two arguments, invariant under a cyclic shift.
  for (int i = 0; i < 2000; i++) {
    GSHGEOM_Point a = P(coord(), coord());
    GSHGEOM_Point b = P(coord(), coord());
    GSHGEOM_Point c = P(coord(), coord());
    int s = o2(a, b, c);
    EXPECT_EQ(o2(b, a, c), -s);
    EXPECT_EQ(o2(a, c, b), -s);
    EXPECT_EQ(o2(b, c, a), s);
    EXPECT_EQ(o2(c, a, b), s);
  }
}

TEST(Orient2d, ZIsIgnored) {
  EXPECT_EQ(o2(P(0, 0, 5), P(1, 0, -7), P(0, 1, 100)), 1);
}

TEST(Orient2d, TheFilterAgreesWithTheExactPath) {
  // Random triples, and triples made nearly collinear in the ways the generator
  // used, with magnitudes from tiny to near the filter's limit and past it.
  const double scales[] = {1.0, 0x1p-80, 0x1p-500, 0x1p-940, 0x1p60, 0x1p99, 0x1p120, 0x1p700};
  for (double sc : scales) {
    for (int i = 0; i < 4000; i++) {
      GSHGEOM_Point a = P(coord() * sc, coord() * sc);
      GSHGEOM_Point b = P(coord() * sc, coord() * sc);
      double t = (double)(next() % 17) / 8.0 - 1.0;
      GSHGEOM_Point c = P(a.x + t * (b.x - a.x), a.y + t * (b.y - a.y));
      if (i % 3 == 0) {
        c.x = std::nextafter(c.x, kInf);
      }
      if (i % 5 == 0) {
        c.y = std::nextafter(c.y, -kInf);
      }
      EXPECT_EQ(o2(a, b, c), o2x(a, b, c)) << sc << " " << i;
      GSHGEOM_Point r = P(coord() * sc, coord() * sc);
      EXPECT_EQ(o2(a, b, r), o2x(a, b, r)) << sc << " " << i;
    }
  }
}

TEST(Orient2d, TheLimitOfTheFilterIsNotAnEdgeOfTheAnswer) {
  // Coordinates at the filter's limit (2^100), one step past it either way, and coincident or collinear triples among them:
  // the answer is the exact one on both sides of the limit, and a point past it is not refused (only a non-finite one is).
  const double edge = 0x1p100;
  const double steps[] = {edge, std::nextafter(edge, kInf), std::nextafter(edge, 0.0), -edge, -std::nextafter(edge, kInf), 0.0,
      1.0, -0.0, 0x1p99};
  for (double x0 : steps) {
    for (double y0 : steps) {
      for (double x1 : steps) {
        for (double y1 : steps) {
          GSHGEOM_Point a = P(x0, y0);
          GSHGEOM_Point b = P(x1, y1);
          GSHGEOM_Point c = P(y0, x1);
          EXPECT_EQ(o2(a, b, c), o2x(a, b, c));
          EXPECT_EQ(o2(a, c, b), o2x(a, c, b));
          EXPECT_EQ(o2(a, a, b), 0);
          EXPECT_EQ(o2(a, b, b), 0);
        }
      }
    }
  }
}

TEST(Orient2d, RefusesNonFiniteInput) {
  int s = 7;
  const double bad[] = {kNaN, kInf, -kInf};
  for (double v : bad) {
    EXPECT_EQ(gshgeom_orient2d(P(v, 0), P(1, 0), P(0, 1), &s), GSHGEOM_ERR_INVALID);
    EXPECT_EQ(gshgeom_orient2d(P(0, 0), P(1, v), P(0, 1), &s), GSHGEOM_ERR_INVALID);
    EXPECT_EQ(gshgeom_orient2d(P(0, 0), P(1, 0), P(0, v), &s), GSHGEOM_ERR_INVALID);
    EXPECT_EQ(gshgeom_orient2d_exact(P(0, 0), P(1, 0), P(v, 1), &s), GSHGEOM_ERR_INVALID);
    // z is not a coordinate in 2D, and so is not looked at.
    EXPECT_EQ(gshgeom_orient2d(P(0, 0, v), P(1, 0), P(0, 1), &s), GSHGEOM_OK);
    s = 7;
  }
  EXPECT_EQ(s, 7) << "untouched by a refusal";
  EXPECT_EQ(gshgeom_orient2d(P(0, 0), P(1, 0), P(0, 1), nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_orient2d_exact(P(0, 0), P(1, 0), P(0, 1), nullptr), GSHGEOM_ERR_INVALID);
}

TEST(Orient3d, SenseAndSymmetries) {
  // d at +z above the counterclockwise (orient2d positive) triangle is positive.
  EXPECT_EQ(o3(P(0, 0, 0), P(1, 0, 0), P(0, 1, 0), P(0, 0, 1)), 1);
  EXPECT_EQ(o3(P(0, 0, 0), P(1, 0, 0), P(0, 1, 0), P(0, 0, -1)), -1);
  EXPECT_EQ(o3(P(0, 0, 0), P(1, 0, 0), P(0, 1, 0), P(5, 5, 0)), 0);
  // It reduces to orient2d of the projections when d is straight above a.
  for (int i = 0; i < 2000; i++) {
    GSHGEOM_Point a = P(coord(), coord(), 3.0);
    GSHGEOM_Point b = P(coord(), coord(), 3.0);
    GSHGEOM_Point c = P(coord(), coord(), 3.0);
    GSHGEOM_Point d = P(a.x, a.y, 4.0);
    EXPECT_EQ(o3(a, b, c, d), o2(a, b, c));
    // Odd permutations of the four points negate it, even ones keep it.
    GSHGEOM_Point e = P(coord(), coord(), coord());
    int s = o3(a, b, c, e);
    EXPECT_EQ(o3(b, a, c, e), -s);
    EXPECT_EQ(o3(a, c, b, e), -s);
    EXPECT_EQ(o3(a, b, e, c), -s);
    EXPECT_EQ(o3(b, c, a, e), s);
  }
}

TEST(Orient3d, ExactZerosAndTheSmallest) {
  EXPECT_EQ(o3(P(0, 0, 0), P(1, 1, 1), P(2, 2, 2), P(7, 8, 9)), 0);
  const double h = std::ldexp(1.0, -1074);
  // A volume of 2^-3222: far below a double, still positive.
  EXPECT_EQ(o3(P(0, 0, 0), P(h, 0, 0), P(0, h, 0), P(0, 0, h)), 1);
  EXPECT_EQ(o3(P(0, 0, 0), P(h, 0, 0), P(0, h, 0), P(0, 0, -h)), -1);
  const double big = 1e100;
  // b-a, c-a, d-a = (-B, B, 0), (-B, 0, B), (0, B, B): the determinant is 2 B^3.
  EXPECT_EQ(o3(P(big, 0, 0), P(0, big, 0), P(0, 0, big), P(big, big, big)), 1);
  EXPECT_EQ(o3x(P(big, 0, 0), P(0, big, 0), P(0, 0, big), P(big, big, big)), 1);
  EXPECT_EQ(o3x(P(0.1, 0.1, 0.1), P(0.2, 0.2, 0.2), P(0.3, 0.3, 0.3), P(0.4, 0.4, 0.4)), 0);
}

TEST(Orient3d, TheFilterAgreesWithTheExactPath) {
  const double scales[] = {1.0, 0x1p-80, 0x1p-300, 0x1p60, 0x1p99, 0x1p120, 0x1p300};
  for (double sc : scales) {
    for (int i = 0; i < 3000; i++) {
      GSHGEOM_Point a = P(coord() * sc, coord() * sc, coord() * sc);
      GSHGEOM_Point b = P(coord() * sc, coord() * sc, coord() * sc);
      GSHGEOM_Point c = P(coord() * sc, coord() * sc, coord() * sc);
      double s = (double)(next() % 17) / 8.0 - 1.0;
      double t = (double)(next() % 17) / 8.0 - 1.0;
      GSHGEOM_Point d = P(a.x + s * (b.x - a.x) + t * (c.x - a.x),
          a.y + s * (b.y - a.y) + t * (c.y - a.y),
          a.z + s * (b.z - a.z) + t * (c.z - a.z));
      if (i % 3 == 0) {
        d.z = std::nextafter(d.z, kInf);
      }
      EXPECT_EQ(o3(a, b, c, d), o3x(a, b, c, d)) << sc << " " << i;
      GSHGEOM_Point r = P(coord() * sc, coord() * sc, coord() * sc);
      EXPECT_EQ(o3(a, b, c, r), o3x(a, b, c, r)) << sc << " " << i;
    }
  }
}

TEST(Orient3d, RefusesNonFiniteInput) {
  int s = 7;
  const double bad[] = {kNaN, kInf, -kInf};
  for (double v : bad) {
    EXPECT_EQ(gshgeom_orient3d(P(0, 0, v), P(1, 0, 0), P(0, 1, 0), P(0, 0, 1), &s), GSHGEOM_ERR_INVALID);
    EXPECT_EQ(gshgeom_orient3d(P(0, 0, 0), P(v, 0, 0), P(0, 1, 0), P(0, 0, 1), &s), GSHGEOM_ERR_INVALID);
    EXPECT_EQ(gshgeom_orient3d(P(0, 0, 0), P(1, 0, 0), P(0, v, 0), P(0, 0, 1), &s), GSHGEOM_ERR_INVALID);
    EXPECT_EQ(gshgeom_orient3d(P(0, 0, 0), P(1, 0, 0), P(0, 1, 0), P(0, 0, v), &s), GSHGEOM_ERR_INVALID);
    EXPECT_EQ(gshgeom_orient3d_exact(P(0, 0, 0), P(1, 0, 0), P(0, 1, 0), P(0, 0, v), &s), GSHGEOM_ERR_INVALID);
  }
  EXPECT_EQ(s, 7);
  EXPECT_EQ(gshgeom_orient3d(P(0, 0, 0), P(1, 0, 0), P(0, 1, 0), P(0, 0, 1), nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_orient3d_exact(P(0, 0, 0), P(1, 0, 0), P(0, 1, 0), P(0, 0, 1), nullptr), GSHGEOM_ERR_INVALID);
}

// ---------------------------------------------------------------------------
// The exact engine
// ---------------------------------------------------------------------------

namespace {

int sum_sign(std::vector<GSHGEOM_ExactTerm> t) {
  return gshgeom_exact_sign(t.data(), t.size());
}

GSHGEOM_ExactTerm T(int neg, double a, double b = 1.0, double c = 1.0, int n = 3) {
  GSHGEOM_ExactTerm t;
  t.neg = neg;
  t.n = n;
  t.f[0] = a;
  t.f[1] = b;
  t.f[2] = c;
  return t;
}

}  // namespace

TEST(ExactSign, Sums) {
  EXPECT_EQ(sum_sign({}), 0);
  EXPECT_EQ(sum_sign({T(0, 3.0, 4.0, 5.0)}), 1);
  EXPECT_EQ(sum_sign({T(1, 3.0, 4.0, 5.0)}), -1);
  EXPECT_EQ(sum_sign({T(0, 0.0, 4.0, 5.0)}), 0);
  EXPECT_EQ(sum_sign({T(0, -3.0, 4.0, 5.0)}), -1);
  EXPECT_EQ(sum_sign({T(1, -3.0, -4.0, 5.0)}), -1);
  // A term cancelled by an equal one, exactly.
  EXPECT_EQ(sum_sign({T(0, 0.1, 0.2, 0.3), T(1, 0.2, 0.1, 0.3)}), 0);
  EXPECT_EQ(sum_sign({T(0, 0.1, 0.2, 0.3), T(1, 0.1, 0.2, 0.3), T(0, 1.0, 1.0, 1.0, 1)}), 1);
  // The smallest survivor of a long cancellation: 2^1000 + 1 - 2^1000 is 1,
  // and 2^1000 + 1 - 2^1000 - 1 is 0; the exponents differ by 1000 bits.
  const double e1000 = std::ldexp(1.0, 1000);
  EXPECT_EQ(sum_sign({T(0, e1000, 1.0, 1.0, 1), T(0, 1.0, 1.0, 1.0, 1), T(1, e1000, 1.0, 1.0, 1)}), 1);
  EXPECT_EQ(sum_sign({T(0, e1000, 1.0, 1.0, 1), T(0, 1.0, 1.0, 1.0, 1), T(1, e1000, 1.0, 1.0, 1), T(1, 1.0, 1.0, 1.0, 1)}), 0);
  // A subnormal factor against a normal number: the exponents of the two kinds
  // of double differ in how they are read, and a sum that mixes them is where
  // an off-by-one shows.  h * 2^100 is 2^-974, a normal number.
  {
    const double h = std::ldexp(1.0, -1074);
    const double p100 = std::ldexp(1.0, 100);
    EXPECT_EQ(sum_sign({T(0, h, p100, 1.0, 2), T(1, std::ldexp(1.0, -974), 1.0, 1.0, 1)}), 0);
    EXPECT_EQ(sum_sign({T(0, std::ldexp(1.0, -974), 1.0, 1.0, 1), T(1, h, p100, 1.0, 2)}), 0);
    EXPECT_EQ(sum_sign({T(0, h, p100, 1.0, 2), T(1, std::ldexp(1.0, -975), 1.0, 1.0, 1)}), 1);
    EXPECT_EQ(sum_sign({T(0, h, p100, 1.0, 2), T(1, std::ldexp(1.0, -973), 1.0, 1.0, 1)}), -1);
    // Two normal factors against one: a bias applied to normal exponents moves
    // the first by twice what it moves the second.
    EXPECT_EQ(sum_sign({T(0, std::ldexp(1.0, -500), std::ldexp(1.0, -474), 1.0, 2), T(1, h, p100, 1.0, 2)}), 0);
    EXPECT_EQ(sum_sign({T(0, std::ldexp(1.0, -500), std::ldexp(1.0, -473), 1.0, 2), T(1, h, p100, 1.0, 2)}), 1);
    EXPECT_EQ(sum_sign({T(0, std::ldexp(1.0, -500), std::ldexp(1.0, -475), 1.0, 2), T(1, h, p100, 1.0, 2)}), -1);
    const double h3 = 3.0 * h;
    EXPECT_EQ(sum_sign({T(0, h3, p100, 1.0, 2), T(1, 3.0 * std::ldexp(1.0, -974), 1.0, 1.0, 1)}), 0);
  }
  // A borrow that runs through every word above, and the carry that undoes it.
  EXPECT_EQ(sum_sign({T(1, 1.0, 1.0, 1.0, 1), T(0, e1000, 1.0, 1.0, 1)}), 1);
  EXPECT_EQ(sum_sign({T(0, 1.0, 1.0, 1.0, 1), T(1, e1000, 1.0, 1.0, 1)}), -1);
  EXPECT_EQ(sum_sign({T(1, 1.0, 1.0, 1.0, 1), T(0, e1000, 1.0, 1.0, 1), T(1, e1000, 1.0, 1.0, 1), T(0, 1.0, 1.0, 1.0, 1)}), 0);
  EXPECT_EQ(sum_sign({T(1, 1.0, 1.0, 1.0, 1), T(1, 1.0, 1.0, 1.0, 1), T(0, 2.0, 1.0, 1.0, 1)}), 0);
  // Many full mantissas added: the carries run up through the words.
  {
    const double mm = 9007199254740991.0;
    std::vector<GSHGEOM_ExactTerm> v;
    for (int i = 0; i < 16; i++) v.push_back(T(0, mm, mm, mm));
    v.push_back(T(1, 16.0 * mm, mm, mm));
    EXPECT_EQ(sum_sign(v), 0);
    v.push_back(T(1, 1.0, 1.0, 1.0, 1));
    EXPECT_EQ(sum_sign(v), -1);
  }
  // Products at both ends of the range: 2^-1074 * 2^-1074 * 2^-1074 against
  // 2^1023 * 2^1023 * 2^-1074 ... the sign is the larger's.
  const double h = std::ldexp(1.0, -1074);
  const double top = std::ldexp(1.0, 1023);
  EXPECT_EQ(sum_sign({T(0, top, top, h), T(1, h, h, h)}), 1);
  EXPECT_EQ(sum_sign({T(0, top, top, h), T(1, top, top, h)}), 0);
  EXPECT_EQ(sum_sign({T(1, top, top, top), T(0, top, top, top), T(0, h, h, h)}), 1);
  EXPECT_EQ(sum_sign({T(1, top, top, top), T(0, top, top, top), T(1, h, h, h)}), -1);
  // Full mantissas: (2^53-1)^3 against the same, shifted.
  const double m = 9007199254740991.0;
  EXPECT_EQ(sum_sign({T(0, m, m, m), T(1, m, m, m)}), 0);
  EXPECT_EQ(sum_sign({T(0, m, m, m), T(1, m, m, std::nextafter(m, 0.0))}), 1);
  // A caller bug: too many terms.
  std::vector<GSHGEOM_ExactTerm> many(GSHGEOM_EXACT_MAX_TERMS + 1, T(0, 1.0));
  EXPECT_EQ(sum_sign(many), 0);
  many.resize(GSHGEOM_EXACT_MAX_TERMS);
  EXPECT_EQ(sum_sign(many), 1);
}

// ---------------------------------------------------------------------------
// The normals and the signed area
// ---------------------------------------------------------------------------

TEST(Rules, NormalsAreTheTangentTurnedByAQuarter) {
  GSHGEOM_Point l = gshgeom_normal_left(P(1, 0));
  EXPECT_EQ(l.x, 0.0);
  EXPECT_EQ(l.y, -1.0) << "left of +x is up the screen: -y";
  GSHGEOM_Point r = gshgeom_normal_right(P(1, 0));
  EXPECT_EQ(r.x, 0.0);
  EXPECT_EQ(r.y, 1.0);
  GSHGEOM_Point l2 = gshgeom_normal_left(P(0, 1));
  EXPECT_EQ(l2.x, 1.0) << "travelling down the screen, left is +x";
  EXPECT_EQ(l2.y, 0.0);
  EXPECT_FALSE(std::signbit(l2.y)) << "no negative zero";
  EXPECT_FALSE(std::signbit(l.x));
  EXPECT_FALSE(std::signbit(r.x));
  GSHGEOM_Point z = gshgeom_normal_left(P(3, 4, 99));
  EXPECT_EQ(z.z, 0.0) << "z is dropped";
  z = gshgeom_normal_right(P(3, 4, 99));
  EXPECT_EQ(z.z, 0.0) << "z is dropped from the right normal too";
  // The left normal is on the negative side of the line, the right on the positive.
  for (int i = 0; i < 500; i++) {
    GSHGEOM_Point a = P(coord(), coord());
    GSHGEOM_Point t = P(coord(), coord());
    if (t.x == 0 && t.y == 0) {
      continue;
    }
    GSHGEOM_Point b = P(a.x + t.x, a.y + t.y);
    GSHGEOM_Point ln = gshgeom_normal_left(t);
    GSHGEOM_Point rn = gshgeom_normal_right(t);
    EXPECT_EQ(o2x(a, b, P(a.x + ln.x, a.y + ln.y)), -1);
    EXPECT_EQ(o2x(a, b, P(a.x + rn.x, a.y + rn.y)), 1);
    EXPECT_EQ(ln.x, 0.0 - rn.x);
    EXPECT_EQ(ln.y, 0.0 - rn.y);
  }
}

TEST(Rules, PolygonSignedArea) {
  const GSHGEOM_Point square[] = {P(0, 0), P(1, 0), P(1, 1), P(0, 1)};
  double a = 99;
  ASSERT_EQ(gshgeom_polygon_signed_area(square, 4, &a), GSHGEOM_OK);
  EXPECT_EQ(a, 1.0) << "+x then +y then -x is the positive sense";
  const GSHGEOM_Point rev[] = {P(0, 1), P(1, 1), P(1, 0), P(0, 0)};
  ASSERT_EQ(gshgeom_polygon_signed_area(rev, 4, &a), GSHGEOM_OK);
  EXPECT_EQ(a, -1.0);
  const GSHGEOM_Point tri[] = {P(0, 0), P(4, 0), P(0, 3)};
  ASSERT_EQ(gshgeom_polygon_signed_area(tri, 3, &a), GSHGEOM_OK);
  EXPECT_EQ(a, 6.0);
  // Where it is positive, the first corner's orient2d is positive too.
  EXPECT_EQ(o2(tri[0], tri[1], tri[2]), 1);
  ASSERT_EQ(gshgeom_polygon_signed_area(tri, 2, &a), GSHGEOM_OK);
  EXPECT_EQ(a, 0.0);
  ASSERT_EQ(gshgeom_polygon_signed_area(nullptr, 0, &a), GSHGEOM_OK);
  EXPECT_EQ(a, 0.0);
  // The area does not depend on where the polygon sits (the fan from the first point).
  const GSHGEOM_Point moved[] = {P(1e6, 1e6), P(1e6 + 4, 1e6), P(1e6, 1e6 + 3)};
  ASSERT_EQ(gshgeom_polygon_signed_area(moved, 3, &a), GSHGEOM_OK);
  EXPECT_EQ(a, 6.0);
}

TEST(Rules, PolygonSignedAreaRefusals) {
  double a = 5;
  EXPECT_EQ(gshgeom_polygon_signed_area(nullptr, 3, &a), GSHGEOM_ERR_INVALID);
  const GSHGEOM_Point tri[] = {P(0, 0), P(4, 0), P(0, 3)};
  EXPECT_EQ(gshgeom_polygon_signed_area(tri, 3, nullptr), GSHGEOM_ERR_INVALID);
  const GSHGEOM_Point nan[] = {P(0, 0), P(kNaN, 0), P(0, 3)};
  EXPECT_EQ(gshgeom_polygon_signed_area(nan, 3, &a), GSHGEOM_ERR_INVALID);
  const GSHGEOM_Point inf[] = {P(0, 0), P(4, kInf), P(0, 3)};
  EXPECT_EQ(gshgeom_polygon_signed_area(inf, 3, &a), GSHGEOM_ERR_INVALID);
  const double big = 1e300;
  const GSHGEOM_Point huge[] = {P(0, 0), P(big, 0), P(0, big)};
  EXPECT_EQ(gshgeom_polygon_signed_area(huge, 3, &a), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(a, 5.0) << "untouched by a refusal";
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
