/**
 * @file
 *
 * Curves (CAP-2): construction, refusal of what is not a number, evaluation
 * and bounds of every kind in 2D and 3D.  The bounds are checked three ways:
 * against closed-form extrema (the quadratic's and the cubic's roots written
 * out, and the circle for rational arcs), against dense sampling of the
 * curve, and across degrees (a curve raised to the next degree is the same
 * curve and must have the same box).
 *
 * Copyright 2026 by Corey Pennycuff
 */

#include "test_helpers.h"

#include <ghoti.io/math/math.h>

#include <algorithm>
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

double axis(const GSHGEOM_Point & p, int a) {
  return a == 0 ? p.x : (a == 1 ? p.y : p.z);
}

// A fixed generator: xorshift64*, so every run sees the same curves.
struct Rng {
  uint64_t s;
  explicit Rng(uint64_t seed) : s(seed) {}
  uint64_t next() {
    s ^= s >> 12;
    s ^= s << 25;
    s ^= s >> 27;
    return s * UINT64_C(2685821657736338717);
  }
  double unit() { return (double)(next() >> 11) * (1.0 / 9007199254740992.0); }
  double range(double lo, double hi) { return lo + (hi - lo) * unit(); }
};

GSHGEOM_Curve make(GSHGEOM_Curve_Kind kind, int dim, const GSHGEOM_Point * p,
    const double * w) {
  GSHGEOM_Curve c;
  GSHGEOM_Result r = GSHGEOM_ERR_INTERNAL;
  switch (kind) {
  case GSHGEOM_CURVE_LINE:
    r = gshgeom_curve_line(dim, p[0], p[1], &c);
    break;
  case GSHGEOM_CURVE_QUADRATIC:
    r = gshgeom_curve_quadratic(dim, p[0], p[1], p[2], &c);
    break;
  case GSHGEOM_CURVE_CUBIC:
    r = gshgeom_curve_cubic(dim, p[0], p[1], p[2], p[3], &c);
    break;
  case GSHGEOM_CURVE_RATIONAL_QUADRATIC:
    r = gshgeom_curve_rational_quadratic(dim, p[0], p[1], p[2], w, &c);
    break;
  case GSHGEOM_CURVE_RATIONAL_CUBIC:
    r = gshgeom_curve_rational_cubic(dim, p[0], p[1], p[2], p[3], w, &c);
    break;
  default:
    break;
  }
  EXPECT_EQ(r, GSHGEOM_OK);
  return c;
}

GSHGEOM_Curve make(GSHGEOM_Curve_Kind kind, int dim,
    const std::vector<GSHGEOM_Point> & p) {
  double w[4] = {1, 1, 1, 1};
  return make(kind, dim, p.data(), w);
}

GSHGEOM_Curve random_curve(Rng & rng, GSHGEOM_Curve_Kind kind, int dim) {
  GSHGEOM_Point p[4];
  double w[4];
  for (int i = 0; i < 4; i++) {
    p[i] = P(rng.range(-10, 10), rng.range(-10, 10),
        dim == 3 ? rng.range(-10, 10) : 0.0);
    w[i] = rng.range(0.2, 5.0);
  }
  return make(kind, dim, p, w);
}

const GSHGEOM_Curve_Kind kKinds[] = {GSHGEOM_CURVE_LINE,
    GSHGEOM_CURVE_QUADRATIC, GSHGEOM_CURVE_CUBIC,
    GSHGEOM_CURVE_RATIONAL_QUADRATIC, GSHGEOM_CURVE_RATIONAL_CUBIC};

GSHGEOM_Point eval(const GSHGEOM_Curve & c, double t) {
  GSHGEOM_Point out = P(-77, -77, -77);
  EXPECT_EQ(gshgeom_curve_eval(&c, t, &out), GSHGEOM_OK) << "t=" << t;
  return out;
}

GSHGEOM_Bounds bounds(const GSHGEOM_Curve & c) {
  GSHGEOM_Bounds b;
  EXPECT_EQ(gshgeom_curve_bounds(&c, &b), GSHGEOM_OK);
  EXPECT_EQ(b.empty, 0);
  return b;
}

// Size of the curve, the unit its tolerances are in.
double scale_of(const GSHGEOM_Curve & c) {
  double s = 1e-300;
  int n = gshgeom_curve_point_count(c.kind);
  for (int i = 0; i < n; i++) {
    for (int a = 0; a < 3; a++) {
      s = std::max(s, std::fabs(axis(c.p[i], a)));
    }
  }
  return s;
}

}  // namespace

// ---------------------------------------------------------------------------
// Construction and refusal
// ---------------------------------------------------------------------------

TEST(Curve, KindsAndPointCounts) {
  EXPECT_EQ(GSHGEOM_CURVE_LINE, 0);
  EXPECT_EQ(GSHGEOM_CURVE_QUADRATIC, 1);
  EXPECT_EQ(GSHGEOM_CURVE_CUBIC, 2);
  EXPECT_EQ(GSHGEOM_CURVE_RATIONAL_QUADRATIC, 3);
  EXPECT_EQ(GSHGEOM_CURVE_RATIONAL_CUBIC, 4);
  EXPECT_EQ(GSHGEOM_CURVE_KIND_COUNT, 5);
  EXPECT_EQ(gshgeom_curve_point_count(GSHGEOM_CURVE_LINE), 2);
  EXPECT_EQ(gshgeom_curve_point_count(GSHGEOM_CURVE_QUADRATIC), 3);
  EXPECT_EQ(gshgeom_curve_point_count(GSHGEOM_CURVE_CUBIC), 4);
  EXPECT_EQ(gshgeom_curve_point_count(GSHGEOM_CURVE_RATIONAL_QUADRATIC), 3);
  EXPECT_EQ(gshgeom_curve_point_count(GSHGEOM_CURVE_RATIONAL_CUBIC), 4);
  EXPECT_EQ(gshgeom_curve_point_count(GSHGEOM_CURVE_KIND_COUNT), 0);
  EXPECT_EQ(gshgeom_curve_point_count((GSHGEOM_Curve_Kind)-1), 0);
  EXPECT_EQ(gshgeom_curve_is_rational(GSHGEOM_CURVE_CUBIC), 0);
  EXPECT_EQ(gshgeom_curve_is_rational(GSHGEOM_CURVE_RATIONAL_QUADRATIC), 1);
  EXPECT_EQ(gshgeom_curve_is_rational(GSHGEOM_CURVE_RATIONAL_CUBIC), 1);
  EXPECT_EQ(gshgeom_curve_is_rational(GSHGEOM_CURVE_KIND_COUNT), 0);
}

TEST(Curve, ConstructorsZeroWhatTheyDoNotUse) {
  GSHGEOM_Curve c;
  ASSERT_EQ(gshgeom_curve_line(2, P(1, 2, 9), P(3, 4, 9), &c), GSHGEOM_OK);
  EXPECT_EQ(c.kind, GSHGEOM_CURVE_LINE);
  EXPECT_EQ(c.dim, 2);
  EXPECT_EQ(c.p[0].z, 0.0);  // z is ignored in 2D
  EXPECT_EQ(c.p[1].z, 0.0);
  for (int i = 2; i < 4; i++) {
    EXPECT_EQ(c.p[i].x, 0.0);
    EXPECT_EQ(c.p[i].y, 0.0);
    EXPECT_EQ(c.w[i], 0.0);
  }
  // A polynomial curve has no weights.
  ASSERT_EQ(gshgeom_curve_cubic(3, P(1, 2, 3), P(4, 5, 6), P(7, 8, 9),
                P(1, 1, 1), &c),
      GSHGEOM_OK);
  for (int i = 0; i < 4; i++) {
    EXPECT_EQ(c.w[i], 0.0);
  }
  EXPECT_EQ(c.p[2].z, 9.0);
  double w[3] = {1, 2, 3};
  ASSERT_EQ(gshgeom_curve_rational_quadratic(
                3, P(1, 2, 3), P(4, 5, 6), P(7, 8, 9), w, &c),
      GSHGEOM_OK);
  EXPECT_EQ(c.w[0], 1.0);
  EXPECT_EQ(c.w[2], 3.0);
  EXPECT_EQ(c.w[3], 0.0);
  EXPECT_EQ(c.p[3].x, 0.0);
}

TEST(Curve, ANaNInTheIgnoredZOf2DIsNotLookedAt) {
  GSHGEOM_Curve c;
  ASSERT_EQ(gshgeom_curve_line(2, P(1, 2, kNaN), P(3, 4, kInf), &c),
      GSHGEOM_OK);
  EXPECT_EQ(c.p[0].z, 0.0);
  EXPECT_EQ(c.p[1].z, 0.0);
  GSHGEOM_Point at = eval(c, 0.5);
  EXPECT_EQ(at.z, 0.0);
  // The same in 3D is refused.
  EXPECT_EQ(gshgeom_curve_line(3, P(1, 2, kNaN), P(3, 4, 0), &c),
      GSHGEOM_ERR_INVALID);
}

TEST(Curve, ConstructorsRefuseEveryNonFiniteComponent) {
  const double bad[] = {kNaN, kInf, -kInf};
  for (int dim = 2; dim <= 3; dim++) {
    for (int kind_i = 0; kind_i < 5; kind_i++) {
      GSHGEOM_Curve_Kind kind = kKinds[kind_i];
      int n = gshgeom_curve_point_count(kind);
      for (int pi = 0; pi < n; pi++) {
        for (int a = 0; a < dim; a++) {
          for (double v : bad) {
            GSHGEOM_Point pts[4] = {P(1, 2, 3), P(4, 5, 6), P(7, 8, 9),
                P(10, 11, 12)};
            double w[4] = {1, 1, 1, 1};
            if (a == 0) pts[pi].x = v;
            if (a == 1) pts[pi].y = v;
            if (a == 2) pts[pi].z = v;
            GSHGEOM_Curve c;
            std::memset(&c, 0x5a, sizeof c);
            GSHGEOM_Curve before = c;
            GSHGEOM_Result r = GSHGEOM_OK;
            switch (kind) {
            case GSHGEOM_CURVE_LINE:
              r = gshgeom_curve_line(dim, pts[0], pts[1], &c);
              break;
            case GSHGEOM_CURVE_QUADRATIC:
              r = gshgeom_curve_quadratic(dim, pts[0], pts[1], pts[2], &c);
              break;
            case GSHGEOM_CURVE_CUBIC:
              r = gshgeom_curve_cubic(
                  dim, pts[0], pts[1], pts[2], pts[3], &c);
              break;
            case GSHGEOM_CURVE_RATIONAL_QUADRATIC:
              r = gshgeom_curve_rational_quadratic(
                  dim, pts[0], pts[1], pts[2], w, &c);
              break;
            default:
              r = gshgeom_curve_rational_cubic(
                  dim, pts[0], pts[1], pts[2], pts[3], w, &c);
              break;
            }
            EXPECT_EQ(r, GSHGEOM_ERR_INVALID)
                << "kind " << kind << " dim " << dim << " point " << pi
                << " axis " << a << " value " << v;
            EXPECT_EQ(std::memcmp(&c, &before, sizeof c), 0)
                << "the output is written only on success";
          }
        }
      }
    }
  }
}

TEST(Curve, WeightsAreRefusedWhenNotPositiveAndFinite) {
  const double invalid[] = {kNaN, kInf, -kInf, 0.0, -0.0, -1.0, -1e-300};
  for (int dim = 2; dim <= 3; dim++) {
    for (int slot = 0; slot < 4; slot++) {
      for (double v : invalid) {
        double w[4] = {1, 2, 3, 4};
        w[slot] = v;
        GSHGEOM_Curve c;
        if (slot < 3) {
          EXPECT_EQ(gshgeom_curve_rational_quadratic(
                        dim, P(0, 0), P(1, 1), P(2, 0), w, &c),
              GSHGEOM_ERR_INVALID)
              << "slot " << slot << " value " << v;
        }
        EXPECT_EQ(gshgeom_curve_rational_cubic(dim, P(0, 0), P(1, 1), P(2, 1),
                      P(3, 0), w, &c),
            GSHGEOM_ERR_INVALID)
            << "slot " << slot << " value " << v;
      }
    }
  }
  // The weight a quadratic does not have is not looked at.
  double w[4] = {1, 2, 3, kNaN};
  GSHGEOM_Curve c;
  EXPECT_EQ(gshgeom_curve_rational_quadratic(
                2, P(0, 0), P(1, 1), P(2, 0), w, &c),
      GSHGEOM_OK);
  EXPECT_EQ(c.w[3], 0.0);
  EXPECT_EQ(gshgeom_curve_rational_quadratic(
                2, P(0, 0), P(1, 1), P(2, 0), nullptr, &c),
      GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_curve_rational_cubic(
                2, P(0, 0), P(1, 1), P(2, 0), P(3, 3), nullptr, &c),
      GSHGEOM_ERR_INVALID);
}

TEST(Curve, MagnitudeBoundsAreLimitsAndNotClamped) {
  GSHGEOM_Curve c;
  const double edge = GSHGEOM_COORD_LIMIT;
  const double over = std::nextafter(edge, kInf);
  EXPECT_EQ(gshgeom_curve_line(3, P(edge, -edge, edge), P(0, 0, 0), &c),
      GSHGEOM_OK);
  EXPECT_EQ(c.p[0].x, edge);
  for (int a = 0; a < 3; a++) {
    for (double sign : {1.0, -1.0}) {
      GSHGEOM_Point q = P(0, 0, 0);
      if (a == 0) q.x = sign * over;
      if (a == 1) q.y = sign * over;
      if (a == 2) q.z = sign * over;
      EXPECT_EQ(gshgeom_curve_line(3, P(0, 0, 0), q, &c), GSHGEOM_ERR_LIMIT);
    }
  }
  // A z out of range is not looked at in 2D, as a NaN is not.
  EXPECT_EQ(gshgeom_curve_line(2, P(0, 0, 0), P(0, 0, over * 4), &c),
      GSHGEOM_OK);

  double ok[3] = {GSHGEOM_WEIGHT_MIN, 1, GSHGEOM_WEIGHT_MAX};
  EXPECT_EQ(gshgeom_curve_rational_quadratic(
                2, P(0, 0), P(1, 1), P(2, 0), ok, &c),
      GSHGEOM_OK);
  for (int slot = 0; slot < 3; slot++) {
    double lo[3] = {1, 1, 1};
    double hi[3] = {1, 1, 1};
    lo[slot] = std::nextafter(GSHGEOM_WEIGHT_MIN, 0.0);
    hi[slot] = std::nextafter(GSHGEOM_WEIGHT_MAX, kInf);
    EXPECT_EQ(gshgeom_curve_rational_quadratic(
                  2, P(0, 0), P(1, 1), P(2, 0), lo, &c),
        GSHGEOM_ERR_LIMIT);
    EXPECT_EQ(gshgeom_curve_rational_quadratic(
                  2, P(0, 0), P(1, 1), P(2, 0), hi, &c),
        GSHGEOM_ERR_LIMIT);
  }
  // A NaN beats a limit: it is reported as what it is, wherever it sits.
  EXPECT_EQ(gshgeom_curve_line(2, P(over, 0), P(kNaN, 0), &c),
      GSHGEOM_ERR_INVALID);
  double mixed[3] = {GSHGEOM_WEIGHT_MAX * 2, 1, -1};
  EXPECT_EQ(gshgeom_curve_rational_quadratic(
                2, P(0, 0), P(1, 1), P(2, 0), mixed, &c),
      GSHGEOM_ERR_INVALID);
}

TEST(Curve, ArgumentsAreChecked) {
  GSHGEOM_Curve c;
  for (int dim : {-1, 0, 1, 4, 100}) {
    EXPECT_EQ(gshgeom_curve_line(dim, P(0, 0), P(1, 1), &c),
        GSHGEOM_ERR_INVALID);
  }
  EXPECT_EQ(gshgeom_curve_line(2, P(0, 0), P(1, 1), nullptr),
      GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_curve_quadratic(2, P(0, 0), P(1, 1), P(2, 0), nullptr),
      GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_curve_check(nullptr), GSHGEOM_ERR_INVALID);
  GSHGEOM_Point at;
  GSHGEOM_Bounds box;
  EXPECT_EQ(gshgeom_curve_eval(nullptr, 0.5, &at), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_curve_bounds(nullptr, &box), GSHGEOM_ERR_INVALID);
  ASSERT_EQ(gshgeom_curve_line(2, P(0, 0), P(1, 1), &c), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_curve_eval(&c, 0.5, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_curve_bounds(&c, nullptr), GSHGEOM_ERR_INVALID);
}

TEST(Curve, AHandFilledStructureIsCheckedLikeAConstructedOne) {
  GSHGEOM_Curve c;
  ASSERT_EQ(gshgeom_curve_cubic(3, P(0, 0, 0), P(1, 2, 3), P(4, 5, 6),
                P(7, 8, 9), &c),
      GSHGEOM_OK);
  EXPECT_EQ(gshgeom_curve_check(&c), GSHGEOM_OK);
  GSHGEOM_Point at;
  GSHGEOM_Bounds box;

  GSHGEOM_Curve bad = c;
  bad.kind = GSHGEOM_CURVE_KIND_COUNT;
  EXPECT_EQ(gshgeom_curve_check(&bad), GSHGEOM_ERR_INVALID);
  bad.kind = (GSHGEOM_Curve_Kind)-1;
  EXPECT_EQ(gshgeom_curve_check(&bad), GSHGEOM_ERR_INVALID);
  bad.kind = (GSHGEOM_Curve_Kind)1000;
  EXPECT_EQ(gshgeom_curve_eval(&bad, 0.5, &at), GSHGEOM_ERR_INVALID);
  bad = c;
  bad.dim = 4;
  EXPECT_EQ(gshgeom_curve_bounds(&bad, &box), GSHGEOM_ERR_INVALID);
  bad.dim = 1;
  EXPECT_EQ(gshgeom_curve_check(&bad), GSHGEOM_ERR_INVALID);
  bad = c;
  bad.p[3].y = kNaN;
  EXPECT_EQ(gshgeom_curve_eval(&bad, 0.5, &at), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_curve_bounds(&bad, &box), GSHGEOM_ERR_INVALID);
  bad = c;
  bad.p[1].z = 2 * GSHGEOM_COORD_LIMIT;
  EXPECT_EQ(gshgeom_curve_eval(&bad, 0.5, &at), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(gshgeom_curve_bounds(&bad, &box), GSHGEOM_ERR_LIMIT);
  // The z of a 2D curve is ignored, whatever it holds.
  GSHGEOM_Curve flat;
  ASSERT_EQ(gshgeom_curve_cubic(2, P(0, 0), P(1, 2), P(4, 5), P(7, 8), &flat),
      GSHGEOM_OK);
  for (int i = 0; i < 4; i++) {
    flat.p[i].z = (i & 1) ? kNaN : kInf;
  }
  EXPECT_EQ(gshgeom_curve_check(&flat), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_curve_eval(&flat, 0.5, &at), GSHGEOM_OK);
  EXPECT_EQ(at.z, 0.0);
  EXPECT_EQ(gshgeom_curve_bounds(&flat, &box), GSHGEOM_OK);
  EXPECT_EQ(box.max.z, 0.0);
  // The fourth point of a quadratic is not used.
  GSHGEOM_Curve q;
  ASSERT_EQ(gshgeom_curve_quadratic(2, P(0, 0), P(1, 1), P(2, 0), &q),
      GSHGEOM_OK);
  q.p[3].x = kNaN;
  EXPECT_EQ(gshgeom_curve_check(&q), GSHGEOM_OK);
}

// ---------------------------------------------------------------------------
// Evaluation
// ---------------------------------------------------------------------------

TEST(CurveEval, ParameterMustBeInUnitInterval) {
  GSHGEOM_Curve c = make(GSHGEOM_CURVE_QUADRATIC, 2,
      std::vector<GSHGEOM_Point>{P(0, 0), P(1, 1), P(2, 0)});
  GSHGEOM_Point out = P(-5, -5, -5);
  for (double t : {kNaN, kInf, -kInf, -1e-300, -0.5, 1.5,
           std::nextafter(1.0, 2.0), std::nextafter(0.0, -1.0), 1e300}) {
    EXPECT_EQ(gshgeom_curve_eval(&c, t, &out), GSHGEOM_ERR_INVALID) << t;
    EXPECT_EQ(out.x, -5.0) << "written only on success";
  }
  EXPECT_EQ(gshgeom_curve_eval(&c, 0.0, &out), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_curve_eval(&c, -0.0, &out), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_curve_eval(&c, 1.0, &out), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_curve_eval(&c, std::nextafter(1.0, 0.0), &out),
      GSHGEOM_OK);
}

TEST(CurveEval, EndsAreExactForEveryKindInBothDimensions) {
  Rng rng(1);
  for (int dim = 2; dim <= 3; dim++) {
    for (GSHGEOM_Curve_Kind kind : kKinds) {
      for (int i = 0; i < 50; i++) {
        GSHGEOM_Curve c = random_curve(rng, kind, dim);
        int n = gshgeom_curve_point_count(kind);
        GSHGEOM_Point a = eval(c, 0.0);
        GSHGEOM_Point b = eval(c, 1.0);
        EXPECT_EQ(std::memcmp(&a, &c.p[0], sizeof a), 0);
        EXPECT_EQ(std::memcmp(&b, &c.p[n - 1], sizeof b), 0);
      }
    }
  }
}

TEST(CurveEval, LineAndQuadraticAndCubicAgainstClosedForms) {
  for (int dim = 2; dim <= 3; dim++) {
    double z = dim == 3 ? 1.0 : 0.0;
    GSHGEOM_Curve l = make(GSHGEOM_CURVE_LINE, dim,
        std::vector<GSHGEOM_Point>{P(1, 2, 3 * z), P(5, -2, 7 * z)});
    GSHGEOM_Point at = eval(l, 0.25);
    EXPECT_EQ(at.x, 2.0);
    EXPECT_EQ(at.y, 1.0);
    EXPECT_EQ(at.z, 4.0 * z);

    GSHGEOM_Curve q = make(GSHGEOM_CURVE_QUADRATIC, dim,
        std::vector<GSHGEOM_Point>{P(0, 0, 0), P(4, 8, 4 * z), P(8, 0, 8 * z)});
    at = eval(q, 0.5);  // (p0 + 2 p1 + p2) / 4
    EXPECT_EQ(at.x, 4.0);
    EXPECT_EQ(at.y, 4.0);
    EXPECT_EQ(at.z, 4.0 * z);
    at = eval(q, 0.25);  // (9 p0 + 6 p1 + p2) / 16
    EXPECT_EQ(at.x, (0.0 + 24.0 + 8.0) / 16.0);
    EXPECT_EQ(at.y, 48.0 / 16.0);

    GSHGEOM_Curve c = make(GSHGEOM_CURVE_CUBIC, dim,
        std::vector<GSHGEOM_Point>{P(0, 0, 0), P(0, 8, 8 * z), P(8, 8, 0), P(8, 0, 8 * z)});
    at = eval(c, 0.5);  // (p0 + 3 p1 + 3 p2 + p3) / 8
    EXPECT_EQ(at.x, 4.0);
    EXPECT_EQ(at.y, 6.0);
    EXPECT_EQ(at.z, 4.0 * z);
  }
}

TEST(CurveEval, CubicAgainstTheBernsteinFormOnRandomCurves) {
  Rng rng(2);
  for (int dim = 2; dim <= 3; dim++) {
    for (int i = 0; i < 200; i++) {
      GSHGEOM_Curve c = random_curve(rng, GSHGEOM_CURVE_CUBIC, dim);
      double t = rng.unit();
      GSHGEOM_Point at = eval(c, t);
      double s = 1 - t;
      double b0 = s * s * s, b1 = 3 * s * s * t, b2 = 3 * s * t * t,
             b3 = t * t * t;
      for (int a = 0; a < dim; a++) {
        double want = b0 * axis(c.p[0], a) + b1 * axis(c.p[1], a) +
            b2 * axis(c.p[2], a) + b3 * axis(c.p[3], a);
        EXPECT_NEAR(axis(at, a), want, 1e-13 * 10);
      }
    }
  }
}

TEST(CurveEval, ARationalCurveWithEqualWeightsIsThePolynomialCurve) {
  Rng rng(3);
  for (int dim = 2; dim <= 3; dim++) {
    for (int i = 0; i < 100; i++) {
      GSHGEOM_Curve q = random_curve(rng, GSHGEOM_CURVE_QUADRATIC, dim);
      GSHGEOM_Curve c = random_curve(rng, GSHGEOM_CURVE_CUBIC, dim);
      double k = rng.range(0.25, 4.0);
      double w3[3] = {k, k, k};
      double w4[4] = {k, k, k, k};
      GSHGEOM_Curve rq, rc;
      ASSERT_EQ(gshgeom_curve_rational_quadratic(
                    dim, q.p[0], q.p[1], q.p[2], w3, &rq),
          GSHGEOM_OK);
      ASSERT_EQ(gshgeom_curve_rational_cubic(
                    dim, c.p[0], c.p[1], c.p[2], c.p[3], w4, &rc),
          GSHGEOM_OK);
      for (int j = 0; j <= 10; j++) {
        double t = j / 10.0;
        GSHGEOM_Point a = eval(q, t), b = eval(rq, t);
        GSHGEOM_Point c1 = eval(c, t), d = eval(rc, t);
        for (int ax = 0; ax < dim; ax++) {
          EXPECT_NEAR(axis(a, ax), axis(b, ax), 1e-13);
          EXPECT_NEAR(axis(c1, ax), axis(d, ax), 1e-13);
        }
      }
    }
  }
}

TEST(CurveEval, ScalingEveryWeightDoesNotMoveTheCurve) {
  Rng rng(4);
  for (int i = 0; i < 100; i++) {
    GSHGEOM_Curve c = random_curve(rng, GSHGEOM_CURVE_RATIONAL_CUBIC, 3);
    GSHGEOM_Curve d = c;
    for (int j = 0; j < 4; j++) {
      d.w[j] = c.w[j] * 8.0;  // a power of two: exact
    }
    for (int j = 0; j <= 8; j++) {
      double t = j / 8.0;
      GSHGEOM_Point a = eval(c, t), b = eval(d, t);
      for (int ax = 0; ax < 3; ax++) {
        EXPECT_NEAR(axis(a, ax), axis(b, ax), 1e-12);
      }
    }
  }
}

TEST(CurveEval, ARationalQuadraticWithTheRightWeightIsACircle) {
  const double w = gmath_sqrt(0.5);
  double ws[3] = {1.0, w, 1.0};
  // The quarter circle from (1, 0) to (0, 1) about the origin, in 2D and in
  // 3D at constant depth.
  for (int dim = 2; dim <= 3; dim++) {
    GSHGEOM_Curve c;
    ASSERT_EQ(gshgeom_curve_rational_quadratic(
                  dim, P(1, 0, 2), P(1, 1, 2), P(0, 1, 2), ws, &c),
        GSHGEOM_OK);
    for (int i = 0; i <= 100; i++) {
      GSHGEOM_Point p = eval(c, i / 100.0);
      EXPECT_NEAR(p.x * p.x + p.y * p.y, 1.0, 4e-16 * 4) << i;
      EXPECT_EQ(p.z, dim == 3 ? 2.0 : 0.0);
    }
  }
}

TEST(CurveEval, FiniteForTheLargestAndSmallestAllowedNumbers) {
  const double big = GSHGEOM_COORD_LIMIT;
  const double ext[2] = {GSHGEOM_WEIGHT_MIN, GSHGEOM_WEIGHT_MAX};
  for (int mask = 0; mask < 16; mask++) {
    double w[4];
    GSHGEOM_Point pts[4];
    for (int i = 0; i < 4; i++) {
      w[i] = ext[(mask >> i) & 1];
      pts[i] = P(((i & 1) ? big : -big), ((i & 2) ? -big : big),
          ((i + mask) & 1 ? big : -big));
    }
    GSHGEOM_Curve c = make(GSHGEOM_CURVE_RATIONAL_CUBIC, 3, pts, w);
    GSHGEOM_Bounds b = bounds(c);
    for (int i = 0; i <= 40; i++) {
      GSHGEOM_Point p = eval(c, i / 40.0);
      for (int a = 0; a < 3; a++) {
        EXPECT_TRUE(std::isfinite(axis(p, a)));
        EXPECT_LE(axis(b.min, a), axis(p, a));
        EXPECT_GE(axis(b.max, a), axis(p, a));
        EXPECT_LE(std::fabs(axis(p, a)), big * (1 + 1e-15));
      }
    }
    for (int a = 0; a < 3; a++) {
      EXPECT_TRUE(std::isfinite(axis(b.min, a)));
      EXPECT_TRUE(std::isfinite(axis(b.max, a)));
    }
  }
}

// ---------------------------------------------------------------------------
// Bounds
// ---------------------------------------------------------------------------

TEST(CurveBounds, ALineIsTheBoxOfItsEnds) {
  for (int dim = 2; dim <= 3; dim++) {
    GSHGEOM_Curve c = make(GSHGEOM_CURVE_LINE, dim,
        std::vector<GSHGEOM_Point>{P(3, -1, 5), P(-2, 4, -6)});
    GSHGEOM_Bounds b = bounds(c);
    EXPECT_EQ(b.min.x, -2.0);
    EXPECT_EQ(b.max.x, 3.0);
    EXPECT_EQ(b.min.y, -1.0);
    EXPECT_EQ(b.max.y, 4.0);
    EXPECT_EQ(b.min.z, dim == 3 ? -6.0 : 0.0);
    EXPECT_EQ(b.max.z, dim == 3 ? 5.0 : 0.0);
  }
}

TEST(CurveBounds, QuadraticAgainstTheClosedFormExtremum) {
  Rng rng(10);
  for (int dim = 2; dim <= 3; dim++) {
    for (int i = 0; i < 500; i++) {
      GSHGEOM_Curve c = random_curve(rng, GSHGEOM_CURVE_QUADRATIC, dim);
      GSHGEOM_Bounds b = bounds(c);
      double tol = 1e-13 * scale_of(c);
      for (int a = 0; a < dim; a++) {
        double p0 = axis(c.p[0], a), p1 = axis(c.p[1], a), p2 = axis(c.p[2], a);
        double lo = std::min(p0, p2), hi = std::max(p0, p2);
        double den = p0 - 2 * p1 + p2;
        if (den != 0.0) {
          double t = (p0 - p1) / den;
          if (t > 0 && t < 1) {
            double v = (1 - t) * (1 - t) * p0 + 2 * (1 - t) * t * p1 + t * t * p2;
            lo = std::min(lo, v);
            hi = std::max(hi, v);
          }
        }
        EXPECT_NEAR(axis(b.min, a), lo, tol) << "axis " << a;
        EXPECT_NEAR(axis(b.max, a), hi, tol) << "axis " << a;
      }
    }
  }
}

TEST(CurveBounds, ASymmetricQuadraticPeaksAtItsMiddle) {
  GSHGEOM_Curve c = make(GSHGEOM_CURVE_QUADRATIC, 2,
      std::vector<GSHGEOM_Point>{P(0, 0), P(1, 2), P(2, 0)});
  GSHGEOM_Bounds b = bounds(c);
  EXPECT_EQ(b.min.x, 0.0);
  EXPECT_EQ(b.max.x, 2.0);
  EXPECT_EQ(b.min.y, 0.0);
  EXPECT_NEAR(b.max.y, 1.0, 1e-16);
}

namespace {
// The roots in (0, 1) of the derivative of a cubic Bezier coordinate, from the
// quadratic formula (the test is allowed the C++ library's sqrt).
void cubic_expected(const GSHGEOM_Curve & c, int a, double * lo, double * hi) {
  double p0 = axis(c.p[0], a), p1 = axis(c.p[1], a), p2 = axis(c.p[2], a),
         p3 = axis(c.p[3], a);
  *lo = std::min(p0, p3);
  *hi = std::max(p0, p3);
  // B'(t) / 3 = A t^2 + B t + C
  double A = -p0 + 3 * p1 - 3 * p2 + p3;
  double B = 2 * (p0 - 2 * p1 + p2);
  double C = p1 - p0;
  double ts[2];
  int n = 0;
  if (A == 0.0) {
    if (B != 0.0) ts[n++] = -C / B;
  } else {
    double disc = B * B - 4 * A * C;
    if (disc >= 0) {
      double sq = std::sqrt(disc);
      double q = -0.5 * (B + (B >= 0 ? sq : -sq));
      ts[n++] = q / A;
      if (q != 0.0) ts[n++] = C / q;
    }
  }
  for (int i = 0; i < n; i++) {
    double t = ts[i];
    if (t > 0 && t < 1) {
      double s = 1 - t;
      double v = s * s * s * p0 + 3 * s * s * t * p1 + 3 * s * t * t * p2 +
          t * t * t * p3;
      *lo = std::min(*lo, v);
      *hi = std::max(*hi, v);
    }
  }
}
}  // namespace

TEST(CurveBounds, CubicAgainstTheClosedFormExtrema) {
  Rng rng(11);
  for (int dim = 2; dim <= 3; dim++) {
    for (int i = 0; i < 1000; i++) {
      GSHGEOM_Curve c = random_curve(rng, GSHGEOM_CURVE_CUBIC, dim);
      GSHGEOM_Bounds b = bounds(c);
      double tol = 1e-12 * scale_of(c);
      for (int a = 0; a < dim; a++) {
        double lo, hi;
        cubic_expected(c, a, &lo, &hi);
        EXPECT_NEAR(axis(b.min, a), lo, tol) << "axis " << a;
        EXPECT_NEAR(axis(b.max, a), hi, tol) << "axis " << a;
      }
    }
  }
}

TEST(CurveBounds, CubicWithTwoInteriorExtremaOnOneAxis) {
  // y(t) = 0, 3, -3, 0 control values: a doubly-bent curve with a maximum and
  // a minimum; the closed form is checked against the same helper and against
  // the values worked out by hand: y = 3 t (1-t)(1-2t)... evaluated at the
  // roots of 6 t^2 - 6 t + 1 = 0.
  GSHGEOM_Curve c = make(GSHGEOM_CURVE_CUBIC, 2,
      std::vector<GSHGEOM_Point>{P(0, 0), P(1, 3), P(2, -3), P(3, 0)});
  GSHGEOM_Bounds b = bounds(c);
  double t1 = 0.5 - std::sqrt(3.0) / 6.0;
  double t2 = 0.5 + std::sqrt(3.0) / 6.0;
  auto y = [](double t) {
    double s = 1 - t;
    return 3 * s * s * t * 3 + 3 * s * t * t * -3;
  };
  EXPECT_NEAR(b.max.y, std::max(y(t1), y(t2)), 1e-14);
  EXPECT_NEAR(b.min.y, std::min(y(t1), y(t2)), 1e-14);
  EXPECT_GT(b.max.y, 0.5);
  EXPECT_LT(b.min.y, -0.5);
  EXPECT_EQ(b.min.x, 0.0);
  EXPECT_EQ(b.max.x, 3.0);
}

TEST(CurveBounds, RationalArcsAgainstTheCircle) {
  // An arc of the unit circle from angle a0 to a1 (span under 180 degrees) is
  // a rational quadratic: ends at the circle, the control point where the
  // tangents meet, weights 1, cos(half span), 1.  Its box is the box of the
  // circle's arc, which is the box of its ends and of the axis points it
  // passes through.
  const double pi = 3.14159265358979323846;
  Rng rng(12);
  for (int i = 0; i < 2000; i++) {
    double a0 = rng.range(-pi, pi);
    double span = rng.range(0.05, pi * 0.98);
    double a1 = a0 + span;
    double m = a0 + span / 2;
    double h = span / 2;
    double w[3] = {1.0, std::cos(h), 1.0};
    GSHGEOM_Curve c;
    ASSERT_EQ(gshgeom_curve_rational_quadratic(2, P(std::cos(a0), std::sin(a0)),
                  P(std::cos(m) / std::cos(h), std::sin(m) / std::cos(h)),
                  P(std::cos(a1), std::sin(a1)), w, &c),
        GSHGEOM_OK);
    GSHGEOM_Bounds b = bounds(c);
    double lo[2] = {std::min(std::cos(a0), std::cos(a1)),
        std::min(std::sin(a0), std::sin(a1))};
    double hi[2] = {std::max(std::cos(a0), std::cos(a1)),
        std::max(std::sin(a0), std::sin(a1))};
    for (int k = -8; k <= 8; k++) {
      double q = k * pi / 2;
      if (q > a0 && q < a1) {
        double x = std::cos(q), y = std::sin(q);
        lo[0] = std::min(lo[0], x);
        hi[0] = std::max(hi[0], x);
        lo[1] = std::min(lo[1], y);
        hi[1] = std::max(hi[1], y);
      }
    }
    EXPECT_NEAR(b.min.x, lo[0], 1e-14) << a0 << " " << span;
    EXPECT_NEAR(b.max.x, hi[0], 1e-14) << a0 << " " << span;
    EXPECT_NEAR(b.min.y, lo[1], 1e-14) << a0 << " " << span;
    EXPECT_NEAR(b.max.y, hi[1], 1e-14) << a0 << " " << span;
  }
}

TEST(CurveBounds, TheQuarterCircleIsExactlyTheUnitSquare) {
  double ws[3] = {1.0, gmath_sqrt(0.5), 1.0};
  GSHGEOM_Curve c;
  ASSERT_EQ(gshgeom_curve_rational_quadratic(
                2, P(1, 0), P(1, 1), P(0, 1), ws, &c),
      GSHGEOM_OK);
  GSHGEOM_Bounds b = bounds(c);
  EXPECT_EQ(b.min.x, 0.0);
  EXPECT_EQ(b.max.x, 1.0);
  EXPECT_EQ(b.min.y, 0.0);
  EXPECT_EQ(b.max.y, 1.0);
}

TEST(CurveBounds, ADegreeRaisedCurveHasTheSameBox) {
  Rng rng(13);
  for (int dim = 2; dim <= 3; dim++) {
    for (int i = 0; i < 500; i++) {
      // Quadratic to cubic: Q1 = (P0 + 2 P1) / 3, Q2 = (2 P1 + P2) / 3.
      GSHGEOM_Curve q = random_curve(rng, GSHGEOM_CURVE_QUADRATIC, dim);
      GSHGEOM_Point e[4];
      e[0] = q.p[0];
      e[3] = q.p[2];
      for (int a = 0; a < 3; a++) {
        double v1 = (axis(q.p[0], a) + 2 * axis(q.p[1], a)) / 3.0;
        double v2 = (2 * axis(q.p[1], a) + axis(q.p[2], a)) / 3.0;
        (a == 0 ? e[1].x : a == 1 ? e[1].y : e[1].z) = v1;
        (a == 0 ? e[2].x : a == 1 ? e[2].y : e[2].z) = v2;
      }
      GSHGEOM_Curve c = make(GSHGEOM_CURVE_CUBIC, dim, e, nullptr);
      GSHGEOM_Bounds bq = bounds(q), bc = bounds(c);
      double tol = 1e-12 * scale_of(q);
      for (int a = 0; a < dim; a++) {
        EXPECT_NEAR(axis(bq.min, a), axis(bc.min, a), tol);
        EXPECT_NEAR(axis(bq.max, a), axis(bc.max, a), tol);
      }

      // Rational quadratic to rational cubic, in homogeneous coordinates:
      // H1 = (H0 + 2 H1) / 3, H2 = (2 H1 + H2) / 3.
      GSHGEOM_Curve rq = random_curve(rng, GSHGEOM_CURVE_RATIONAL_QUADRATIC, dim);
      double w[4];
      GSHGEOM_Point r[4];
      r[0] = rq.p[0];
      r[3] = rq.p[2];
      w[0] = rq.w[0];
      w[3] = rq.w[2];
      w[1] = (rq.w[0] + 2 * rq.w[1]) / 3.0;
      w[2] = (2 * rq.w[1] + rq.w[2]) / 3.0;
      for (int a = 0; a < 3; a++) {
        double h0 = rq.w[0] * axis(rq.p[0], a);
        double h1 = rq.w[1] * axis(rq.p[1], a);
        double h2 = rq.w[2] * axis(rq.p[2], a);
        double v1 = ((h0 + 2 * h1) / 3.0) / w[1];
        double v2 = ((2 * h1 + h2) / 3.0) / w[2];
        (a == 0 ? r[1].x : a == 1 ? r[1].y : r[1].z) = v1;
        (a == 0 ? r[2].x : a == 1 ? r[2].y : r[2].z) = v2;
      }
      GSHGEOM_Curve rc = make(GSHGEOM_CURVE_RATIONAL_CUBIC, dim, r, w);
      GSHGEOM_Bounds b1 = bounds(rq), b2 = bounds(rc);
      double tol2 = 1e-11 * scale_of(rq);
      for (int a = 0; a < dim; a++) {
        EXPECT_NEAR(axis(b1.min, a), axis(b2.min, a), tol2);
        EXPECT_NEAR(axis(b1.max, a), axis(b2.max, a), tol2);
      }
    }
  }
}

TEST(CurveBounds, DenseSamplingNeverEscapesAndComesCloseToTheBox) {
  // The box has to contain every point of the curve (to the rounding of the
  // evaluation) and has to be reached by the samples to within what a sample
  // spacing of 1/4000 allows near a smooth extremum.
  const int kSamples = 4000;
  Rng rng(14);
  for (int dim = 2; dim <= 3; dim++) {
    for (GSHGEOM_Curve_Kind kind : kKinds) {
      for (int i = 0; i < 150; i++) {
        GSHGEOM_Curve c = random_curve(rng, kind, dim);
        GSHGEOM_Bounds b = bounds(c);
        double s = scale_of(c);
        double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
        for (int k = 0; k <= kSamples; k++) {
          GSHGEOM_Point p = eval(c, (double)k / kSamples);
          for (int a = 0; a < dim; a++) {
            lo[a] = std::min(lo[a], axis(p, a));
            hi[a] = std::max(hi[a], axis(p, a));
          }
        }
        for (int a = 0; a < dim; a++) {
          EXPECT_LE(axis(b.min, a), lo[a] + 1e-13 * s)
              << "kind " << kind << " axis " << a;
          EXPECT_GE(axis(b.max, a), hi[a] - 1e-13 * s)
              << "kind " << kind << " axis " << a;
          EXPECT_GE(axis(b.min, a), lo[a] - 1e-5 * s)
              << "kind " << kind << " axis " << a << ": the box is too big";
          EXPECT_LE(axis(b.max, a), hi[a] + 1e-5 * s)
              << "kind " << kind << " axis " << a << ": the box is too big";
        }
        for (int a = dim; a < 3; a++) {
          EXPECT_EQ(axis(b.min, a), 0.0);
          EXPECT_EQ(axis(b.max, a), 0.0);
        }
      }
    }
  }
}

TEST(CurveBounds, TheBoxIsNotTheBoxOfTheControlPolygon) {
  // Control points far above the curve: the polygon's box is much bigger.
  GSHGEOM_Curve c = make(GSHGEOM_CURVE_CUBIC, 2,
      std::vector<GSHGEOM_Point>{P(0, 0), P(0, 10), P(10, 10), P(10, 0)});
  GSHGEOM_Bounds b = bounds(c);
  EXPECT_NEAR(b.max.y, 7.5, 1e-14);  // 3/4 of the control height
  EXPECT_EQ(b.min.y, 0.0);
  EXPECT_EQ(b.min.x, 0.0);
  EXPECT_EQ(b.max.x, 10.0);
}

TEST(CurveBounds, DegenerateCurves) {
  // All control points equal.
  for (GSHGEOM_Curve_Kind kind : kKinds) {
    GSHGEOM_Point pts[4] = {P(2, 3, 4), P(2, 3, 4), P(2, 3, 4), P(2, 3, 4)};
    double w[4] = {1, 2, 3, 4};
    GSHGEOM_Curve c = make(kind, 3, pts, w);
    GSHGEOM_Bounds b = bounds(c);
    EXPECT_EQ(b.min.x, 2.0);
    EXPECT_EQ(b.max.x, 2.0);
    EXPECT_EQ(b.min.y, 3.0);
    EXPECT_EQ(b.max.y, 3.0);
    EXPECT_EQ(b.min.z, 4.0);
    EXPECT_EQ(b.max.z, 4.0);
  }
  // Collinear control points out of order: the curve doubles back on a line.
  GSHGEOM_Curve c = make(GSHGEOM_CURVE_CUBIC, 2,
      std::vector<GSHGEOM_Point>{P(0, 0), P(3, 3), P(-1, -1), P(1, 1)});
  GSHGEOM_Bounds b = bounds(c);
  GSHGEOM_Point p;
  double lo = 1e9, hi = -1e9;
  for (int i = 0; i <= 10000; i++) {
    p = eval(c, i / 10000.0);
    lo = std::min(lo, p.x);
    hi = std::max(hi, p.x);
  }
  EXPECT_NEAR(b.min.x, lo, 1e-7);
  EXPECT_NEAR(b.max.x, hi, 1e-7);
  EXPECT_LE(b.min.x, lo + 1e-14);
  EXPECT_GE(b.max.x, hi - 1e-14);
  EXPECT_EQ(b.min.x, b.min.y);
  EXPECT_EQ(b.max.x, b.max.y);
  // A cubic whose coordinate is constant has a zero-width box on that axis,
  // exactly.
  c = make(GSHGEOM_CURVE_CUBIC, 2,
      std::vector<GSHGEOM_Point>{P(0, 5), P(3, 5), P(-1, 5), P(1, 5)});
  b = bounds(c);
  EXPECT_EQ(b.min.y, 5.0);
  EXPECT_EQ(b.max.y, 5.0);
  // A loop: the curve crosses itself, so the same point is reached twice.
  c = make(GSHGEOM_CURVE_CUBIC, 2,
      std::vector<GSHGEOM_Point>{P(0, 0), P(2, 2), P(-1, 2), P(1, 0)});
  b = bounds(c);
  lo = 1e9;
  hi = -1e9;
  double ylo = 1e9, yhi = -1e9;
  for (int i = 0; i <= 20000; i++) {
    p = eval(c, i / 20000.0);
    lo = std::min(lo, p.x);
    hi = std::max(hi, p.x);
    ylo = std::min(ylo, p.y);
    yhi = std::max(yhi, p.y);
  }
  EXPECT_NEAR(b.min.x, lo, 1e-6);
  EXPECT_NEAR(b.max.x, hi, 1e-6);
  EXPECT_NEAR(b.min.y, ylo, 1e-6);
  EXPECT_NEAR(b.max.y, yhi, 1e-6);
}

TEST(CurveBounds, SCurvesWithTwoExtremaPerAxisAreFoundInEveryKind) {
  // Control values that go up, down and up again give two interior extrema
  // per coordinate, so the root isolation has to cut [0, 1] correctly.  A
  // missed one is a large error, not a rounding.
  const int kSamples = 20000;
  Rng rng(16);
  for (int dim = 2; dim <= 3; dim++) {
    for (GSHGEOM_Curve_Kind kind : {GSHGEOM_CURVE_CUBIC,
             GSHGEOM_CURVE_RATIONAL_CUBIC}) {
      for (int i = 0; i < 300; i++) {
        GSHGEOM_Point p[4];
        double w[4];
        for (int k = 0; k < 4; k++) {
          double sgn = (k == 1) ? 1.0 : (k == 2 ? -1.0 : 0.0);
          p[k] = P(10 * sgn + rng.range(-1, 1), -10 * sgn + rng.range(-1, 1),
              dim == 3 ? 10 * sgn + rng.range(-1, 1) : 0.0);
          w[k] = rng.range(0.3, 3.0);
        }
        GSHGEOM_Curve c = make(kind, dim, p, w);
        GSHGEOM_Bounds b = bounds(c);
        double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
        for (int k = 0; k <= kSamples; k++) {
          GSHGEOM_Point q = eval(c, (double)k / kSamples);
          for (int a = 0; a < dim; a++) {
            lo[a] = std::min(lo[a], axis(q, a));
            hi[a] = std::max(hi[a], axis(q, a));
          }
        }
        for (int a = 0; a < dim; a++) {
          EXPECT_NEAR(axis(b.min, a), lo[a], 1e-6) << "kind " << kind;
          EXPECT_NEAR(axis(b.max, a), hi[a], 1e-6) << "kind " << kind;
        }
      }
    }
  }
}

TEST(CurveBounds, RationalCubicsWhereEveryLevelOfTheRootIsolationMatters) {
  // Found by searching: rational cubics (weights from 0.05 to 30) for which a
  // wrong root of the innermost, linear, polynomial changes the box.  The box
  // must still be that of the curve.
  struct Case {
    double x[4];
    double w[4];
  };
  const Case cases[] = {
      {{9.9, -6.6, -1.6, -1.8}, {0.2, 2.24, 11.67, 0.87}},
      {{-7.4, 9.4, -3.5, -1.9}, {0.06, 1.77, 26.64, 0.64}},
      {{-3.1, -1.5, -9.4, 7.1}, {0.44, 27.83, 3.26, 0.59}},
      {{5.8, 5.8, -7.1, 8.4}, {0.06, 1.79, 0.05, 14.49}},
      {{-0.9, -0.6, -9.9, 4.0}, {0.15, 2.9, 0.27, 0.89}},
      {{-1.9, -0.3, -9.1, 7.3}, {2.06, 14.74, 2.39, 0.14}},
  };
  const int kSamples = 400000;
  for (const Case & k : cases) {
    GSHGEOM_Point p[4];
    for (int i = 0; i < 4; i++) {
      p[i] = P(k.x[i], k.x[i]);
    }
    GSHGEOM_Curve c = make(GSHGEOM_CURVE_RATIONAL_CUBIC, 2, p, k.w);
    GSHGEOM_Bounds b = bounds(c);
    double lo = 1e300, hi = -1e300;
    for (int s = 0; s <= kSamples; s++) {
      GSHGEOM_Point q = eval(c, (double)s / kSamples);
      lo = std::min(lo, q.x);
      hi = std::max(hi, q.x);
    }
    EXPECT_NEAR(b.min.x, lo, 1e-7) << k.x[0] << " " << k.w[0];
    EXPECT_NEAR(b.max.x, hi, 1e-7) << k.x[0] << " " << k.w[0];
    EXPECT_NEAR(b.min.y, lo, 1e-7);
    EXPECT_NEAR(b.max.y, hi, 1e-7);
  }
}

TEST(CurveBounds, WildlyUnequalWeightsStillBoundTheCurve) {
  Rng rng(15);
  for (int i = 0; i < 300; i++) {
    GSHGEOM_Curve c = random_curve(rng, GSHGEOM_CURVE_RATIONAL_CUBIC, 3);
    for (int j = 0; j < 4; j++) {
      c.w[j] = std::ldexp(1.0, (int)(rng.next() % 61) - 30);  // 2^-30 .. 2^30
    }
    GSHGEOM_Bounds b = bounds(c);
    double s = scale_of(c);
    for (int k = 0; k <= 2000; k++) {
      GSHGEOM_Point p = eval(c, k / 2000.0);
      for (int a = 0; a < 3; a++) {
        EXPECT_LE(axis(b.min, a), axis(p, a) + 1e-12 * s);
        EXPECT_GE(axis(b.max, a), axis(p, a) - 1e-12 * s);
      }
    }
  }
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
