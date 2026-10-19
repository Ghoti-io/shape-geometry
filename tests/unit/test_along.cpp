/**
 * @file
 *
 * Along-path parametrization (CAP-5, CAP-13): lengths, points, tangents and
 * normals against closed forms (lines, circular arcs as rational curves) and
 * against an independent high-precision reference (cubics, rational curves,
 * a cusp, a loop, 3D), the tolerance and its refinement, monotonicity,
 * continuity and orthonormality, degenerate and hostile input, the limits,
 * allocation failure, and the hash and the dump.
 *
 * Copyright 2026 by Corey Pennycuff
 */

#include <cstdlib>
#include <set>

#include "test_helpers.h"

#include "../alloc_tracker.h"
#include "../along_vectors.h"
#include "../known_answers.h"
#include "../dump_sink.h"

#include <algorithm>
#include <clocale>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace {

const double kPi = 3.14159265358979323846;

GSHGEOM_Point P(double x, double y, double z = 0.0) {
  return GSHGEOM_Point{x, y, z};
}

GSHGEOM_Path * new_path(int dim = 2, const GSHGEOM_Allocator * a = nullptr) {
  GSHGEOM_Path * p = nullptr;
  EXPECT_EQ(gshgeom_path_create(a, nullptr, dim, &p), GSHGEOM_OK);
  return p;
}

GSHGEOM_Measured * measure(const GSHGEOM_Path * p, double tol = 0.0,
    const GSHGEOM_Limits * l = nullptr) {
  GSHGEOM_Measured * m = nullptr;
  EXPECT_EQ(gshgeom_measured_create(nullptr, l, p, tol, &m), GSHGEOM_OK);
  return m;
}

GSHGEOM_Measured_Frame at(const GSHGEOM_Measured * m, double s) {
  GSHGEOM_Measured_Frame f;
  std::memset(&f, 0, sizeof f);
  EXPECT_EQ(gshgeom_measured_at(m, s, &f), GSHGEOM_OK) << s;
  return f;
}

double tol_of(double requested) {
  return requested == 0.0 ? GSHGEOM_MEASURE_TOLERANCE_DEFAULT : requested;
}

double dist2(const GSHGEOM_Point & a, const GSHGEOM_Point & b) {
  return std::hypot(a.x - b.x, a.y - b.y);
}

double dist3(const GSHGEOM_Point & a, const GSHGEOM_Point & b) {
  return std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) +
      (a.z - b.z) * (a.z - b.z));
}

// The path of one reference curve: a move to its first point and the verb.
GSHGEOM_Path * path_of(const AlongRef & r) {
  GSHGEOM_Path * p = new_path(r.dim);
  auto pt = [&](int i) { return P(r.p[i][0], r.p[i][1], r.p[i][2]); };
  EXPECT_EQ(gshgeom_path_move_to(p, pt(0)), GSHGEOM_OK);
  switch (r.kind) {
  case 1:
    EXPECT_EQ(gshgeom_path_quadratic_to(p, pt(1), pt(2)), GSHGEOM_OK);
    break;
  case 2:
    EXPECT_EQ(gshgeom_path_cubic_to(p, pt(1), pt(2), pt(3)), GSHGEOM_OK);
    break;
  case 3:
    EXPECT_EQ(gshgeom_path_rational_quadratic_to(p, pt(1), pt(2), r.w), GSHGEOM_OK);
    break;
  default:
    EXPECT_EQ(gshgeom_path_rational_cubic_to(p, pt(1), pt(2), pt(3), r.w), GSHGEOM_OK);
    break;
  }
  return p;
}

const AlongRef * ref_named(const char * name) {
  for (const AlongRef & r : kAlongRefs) {
    if (std::strcmp(r.name, name) == 0) return &r;
  }
  return nullptr;
}

// An xorshift64* generator: the same numbers everywhere.
struct Rng {
  uint64_t s;
  uint64_t next() {
    s ^= s >> 12;
    s ^= s << 25;
    s ^= s >> 27;
    return s * UINT64_C(2685821657736338717);
  }
  double unit() { return (double)(next() >> 11) / 9007199254740992.0; }
  double range(double lo, double hi) { return lo + (hi - lo) * unit(); }
};

std::string dump_of(const GSHGEOM_Measured * m) {
  FILE * f = std::tmpfile();
  EXPECT_NE(f, nullptr);
  EXPECT_EQ(gshgeom_measured_dump(m, f), GSHGEOM_OK);
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

// ---------------------------------------------------------------------------
// Lines
// ---------------------------------------------------------------------------

TEST(AlongLine, TwoDimensionsExactly) {
  GSHGEOM_Path * p = new_path();
  ASSERT_EQ(gshgeom_path_move_to(p, P(0, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(p, P(3, 4)), GSHGEOM_OK);
  GSHGEOM_Measured * m = measure(p);
  EXPECT_EQ(gshgeom_measured_length(m), 5.0);
  EXPECT_EQ(gshgeom_measured_dimension(m), 2);
  EXPECT_EQ(gshgeom_measured_segment_count(m), 1u);
  EXPECT_EQ(gshgeom_measured_row_count(m), 1u) << "a line is one row";
  for (double s : {0.0, 0.5, 1.0, 2.5, 4.0, 5.0}) {
    GSHGEOM_Measured_Frame f = at(m, s);
    EXPECT_NEAR(f.point.x, 0.6 * s, 4e-16 * 5) << s;
    EXPECT_NEAR(f.point.y, 0.8 * s, 4e-16 * 5) << s;
    EXPECT_NEAR(f.tangent.x, 0.6, 2e-16);
    EXPECT_NEAR(f.tangent.y, 0.8, 2e-16);
    EXPECT_EQ(f.tangent.z, 0.0);
    EXPECT_TRUE(f.has_tangent);
    EXPECT_TRUE(f.has_normal);
    // The left normal is the tangent turned by a negative angle: (t.y, -t.x).
    EXPECT_EQ(f.normal.x, f.tangent.y);
    EXPECT_EQ(f.normal.y, -f.tangent.x);
    EXPECT_EQ(f.element, 1u);
    EXPECT_EQ(f.point.z, 0.0);
  }
  EXPECT_EQ(at(m, 0.0).t, 0.0);
  EXPECT_EQ(at(m, 5.0).t, 1.0);
  EXPECT_EQ(at(m, 5.0).point.x, 3.0);
  EXPECT_EQ(at(m, 5.0).point.y, 4.0);
  gshgeom_measured_destroy(m);
  gshgeom_path_destroy(p);
}

TEST(AlongLine, DyadicDistancesGiveExactPoints) {
  // s = 5 k / 8 is exact and s / 5 = k / 8 is exact, so the parameter and the
  // point are exact: a line is closed form, not a solved equation.
  GSHGEOM_Path * p = new_path();
  ASSERT_EQ(gshgeom_path_move_to(p, P(0, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(p, P(3, 4)), GSHGEOM_OK);
  GSHGEOM_Measured * m = measure(p);
  for (int k = 0; k <= 8; k++) {
    GSHGEOM_Measured_Frame f = at(m, 5.0 * k / 8.0);
    EXPECT_EQ(f.t, k / 8.0) << k;
    EXPECT_EQ(f.point.x, 3.0 * k / 8.0) << k;
    EXPECT_EQ(f.point.y, 4.0 * k / 8.0) << k;
  }
  gshgeom_measured_destroy(m);
  gshgeom_path_destroy(p);
}

TEST(AlongLine, ThreeDimensionsHaveNoNormal) {
  GSHGEOM_Path * p = new_path(3);
  ASSERT_EQ(gshgeom_path_move_to(p, P(1, 2, 3)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(p, P(4, 6, 15)), GSHGEOM_OK);
  GSHGEOM_Measured * m = measure(p);
  EXPECT_EQ(gshgeom_measured_length(m), 13.0);
  EXPECT_EQ(gshgeom_measured_dimension(m), 3);
  GSHGEOM_Measured_Frame f = at(m, 6.5);
  EXPECT_NEAR(f.point.x, 2.5, 1e-15);
  EXPECT_NEAR(f.point.y, 4.0, 1e-15);
  EXPECT_NEAR(f.point.z, 9.0, 1e-15);
  EXPECT_NEAR(f.tangent.x, 3.0 / 13, 1e-16);
  EXPECT_NEAR(f.tangent.y, 4.0 / 13, 1e-16);
  EXPECT_NEAR(f.tangent.z, 12.0 / 13, 1e-16);
  EXPECT_TRUE(f.has_tangent);
  EXPECT_FALSE(f.has_normal) << "left and right have no meaning in 3D";
  EXPECT_EQ(f.normal.x, 0.0);
  EXPECT_EQ(f.normal.y, 0.0);
  EXPECT_EQ(f.normal.z, 0.0);
  gshgeom_measured_destroy(m);
  gshgeom_path_destroy(p);
}

// ---------------------------------------------------------------------------
// Circular arcs, exactly representable as rational curves
// ---------------------------------------------------------------------------

struct Arc {
  GSHGEOM_Point c;
  double r, alpha, theta; // the start angle and the signed sweep
};

GSHGEOM_Point on_circle(const Arc & a, double phi, double radius) {
  return P(a.c.x + radius * std::cos(phi), a.c.y + radius * std::sin(phi));
}

// The arc as a rational quadratic (exact for |theta| < pi) or, with `elevated`,
// as the rational cubic of the same curve.
GSHGEOM_Path * arc_path(const Arc & a, bool elevated) {
  double half = std::fabs(a.theta) / 2;
  double mid = a.alpha + a.theta / 2;
  GSHGEOM_Point p0 = on_circle(a, a.alpha, a.r);
  GSHGEOM_Point p2 = on_circle(a, a.alpha + a.theta, a.r);
  GSHGEOM_Point p1 = on_circle(a, mid, a.r / std::cos(half));
  double w[3] = {1.0, std::cos(half), 1.0};
  GSHGEOM_Path * p = new_path();
  EXPECT_EQ(gshgeom_path_move_to(p, p0), GSHGEOM_OK);
  if (!elevated) {
    EXPECT_EQ(gshgeom_path_rational_quadratic_to(p, p1, p2, w), GSHGEOM_OK);
    return p;
  }
  // Degree elevation in homogeneous coordinates:
  // H1' = (H0 + 2 H1) / 3 and H2' = (2 H1 + H2) / 3.
  GSHGEOM_Point pts[3] = {p0, p1, p2};
  double hx[3], hy[3];
  for (int i = 0; i < 3; i++) {
    hx[i] = w[i] * pts[i].x;
    hy[i] = w[i] * pts[i].y;
  }
  double w1 = (w[0] + 2 * w[1]) / 3, w2 = (2 * w[1] + w[2]) / 3;
  GSHGEOM_Point c1 = P((hx[0] + 2 * hx[1]) / 3 / w1, (hy[0] + 2 * hy[1]) / 3 / w1);
  GSHGEOM_Point c2 = P((2 * hx[1] + hx[2]) / 3 / w2, (2 * hy[1] + hy[2]) / 3 / w2);
  double wc[4] = {w[0], w1, w2, w[2]};
  EXPECT_EQ(gshgeom_path_rational_cubic_to(p, c1, c2, p2, wc), GSHGEOM_OK);
  return p;
}

void check_arc(const Arc & a, bool elevated, double tol_request) {
  SCOPED_TRACE(testing::Message() << "r=" << a.r << " theta=" << a.theta << " alpha=" << a.alpha
                                  << " elevated=" << elevated << " tol=" << tol_request);
  GSHGEOM_Path * p = arc_path(a, elevated);
  GSHGEOM_Measured * m = measure(p, tol_request);
  double tol = tol_of(tol_request);
  double half = std::fabs(a.theta) / 2;
  double poly = 2 * a.r * std::tan(half); // the control polygon's length
  double scale = a.r + std::fabs(a.c.x) + std::fabs(a.c.y);
  double sgn = a.theta < 0 ? -1.0 : 1.0;
  double length = a.r * std::fabs(a.theta);
  // The estimated error is at most tol * poly; the control points are doubles
  // computed from cos and sin, which add a few units in the last place of the
  // largest coordinate.
  double len_bound = tol * poly + 4e-15 * scale;
  EXPECT_NEAR(gshgeom_measured_length(m), length, len_bound);
  double pos_bound = 3 * tol * poly + 1e-14 * scale;
  double dir_bound = pos_bound / a.r + 1e-14 * scale / a.r;
  for (int k = 0; k <= 16; k++) {
    double s = length * k / 16.0;
    GSHGEOM_Measured_Frame f;
    ASSERT_EQ(gshgeom_measured_at_clamped(m, s, &f), GSHGEOM_OK);
    double phi = a.alpha + sgn * s / a.r;
    GSHGEOM_Point want = on_circle(a, phi, a.r);
    EXPECT_NEAR(f.point.x, want.x, pos_bound) << k;
    EXPECT_NEAR(f.point.y, want.y, pos_bound) << k;
    // Positive angle turns +x toward +y; travelling that way the tangent is
    // (-sin, cos) and its left normal, (t.y, -t.x), points out.
    double tx = -sgn * std::sin(phi), ty = sgn * std::cos(phi);
    EXPECT_NEAR(f.tangent.x, tx, dir_bound) << k;
    EXPECT_NEAR(f.tangent.y, ty, dir_bound) << k;
    EXPECT_NEAR(f.normal.x, ty, dir_bound) << k;
    EXPECT_NEAR(f.normal.y, -tx, dir_bound) << k;
  }
  gshgeom_measured_destroy(m);
  gshgeom_path_destroy(p);
}

TEST(AlongArc, RationalQuadraticMatchesTheClosedForm) {
  const double thetas[] = {1.0, 15.0, 45.0, 90.0, 120.0, 170.0, 179.0, -30.0, -170.0};
  const double radii[] = {1.0, 0.001, 1.0e6};
  for (double r : radii) {
    for (double deg : thetas) {
      check_arc(Arc{P(3.0, -2.0), r, 0.3, deg * kPi / 180}, false, 0.0);
    }
  }
}

TEST(AlongArc, RationalCubicByDegreeElevationMatchesToo) {
  const double thetas[] = {5.0, 60.0, 90.0, 150.0, -80.0};
  for (double r : {1.0, 250.0}) {
    for (double deg : thetas) {
      check_arc(Arc{P(-7.0, 11.0), r, 2.0, deg * kPi / 180}, true, 0.0);
    }
  }
}

TEST(AlongArc, LooserAndTighterTolerancesStayInsideTheirBound) {
  for (double tol : {0x1p-10, 0x1p-20, 0x1p-30, 0x1p-46}) {
    check_arc(Arc{P(0.0, 0.0), 10.0, 0.0, 100.0 * kPi / 180}, false, tol);
    check_arc(Arc{P(0.0, 0.0), 10.0, 1.0, -160.0 * kPi / 180}, true, tol);
  }
}

TEST(AlongArc, AFullCircleOfFourQuartersIsClosed) {
  for (double r : {1.0, 1000.0}) {
    double w[3] = {1.0, std::sqrt(0.5), 1.0};
    GSHGEOM_Path * p = new_path();
    ASSERT_EQ(gshgeom_path_move_to(p, P(r, 0)), GSHGEOM_OK);
    const double q[4][4] = {{r, r, 0, r}, {-r, r, -r, 0}, {-r, -r, 0, -r}, {r, -r, r, 0}};
    for (auto & s : q) {
      ASSERT_EQ(gshgeom_path_rational_quadratic_to(p, P(s[0], s[1]), P(s[2], s[3]), w), GSHGEOM_OK);
    }
    ASSERT_EQ(gshgeom_path_close(p), GSHGEOM_OK);
    GSHGEOM_Measured * m = measure(p);
    double tol = GSHGEOM_MEASURE_TOLERANCE_DEFAULT;
    EXPECT_NEAR(gshgeom_measured_length(m), 2 * kPi * r, 4 * tol * 2 * r + 4e-15 * r);
    EXPECT_EQ(gshgeom_measured_segment_count(m), 4u) << "the closing line has length zero";
    ASSERT_EQ(gshgeom_measured_subpath_count(m), 1u);
    GSHGEOM_Measured_Subpath sp;
    ASSERT_EQ(gshgeom_measured_subpath(m, 0, &sp), GSHGEOM_OK);
    EXPECT_EQ(sp.begin, 0.0);
    EXPECT_EQ(sp.end, gshgeom_measured_length(m));
    EXPECT_TRUE(sp.closed);
    double length = gshgeom_measured_length(m);
    for (int k = 0; k <= 2000; k++) {
      double s = length * k / 2000.0;
      GSHGEOM_Measured_Frame f = at(m, s);
      double phi = s / r;
      EXPECT_NEAR(f.point.x, r * std::cos(phi), 2e-11 * r) << k;
      EXPECT_NEAR(f.point.y, r * std::sin(phi), 2e-11 * r) << k;
      EXPECT_NEAR(f.normal.x, std::cos(phi), 2e-11) << k;
      EXPECT_NEAR(f.normal.y, std::sin(phi), 2e-11) << k;
    }
    gshgeom_measured_destroy(m);
    gshgeom_path_destroy(p);
  }
}

// ---------------------------------------------------------------------------
// Curves against the independent reference
// ---------------------------------------------------------------------------

TEST(AlongReference, LengthIsWithinTheToleranceAtEveryToleranceAndRowsGrow) {
  const double tols[] = {0x1p-10, 0x1p-16, 0x1p-24, 0x1p-32, 0x1p-40, 0.0};
  ASSERT_GE(sizeof kAlongRefs / sizeof kAlongRefs[0], 14u);
  for (const AlongRef & r : kAlongRefs) {
    SCOPED_TRACE(r.name);
    GSHGEOM_Path * p = path_of(r);
    size_t rows_before = 0;
    double worst_coarse = 0, worst_fine = 0;
    for (double tol : tols) {
      GSHGEOM_Measured * m = measure(p, tol);
      double got = gshgeom_measured_length(m);
      double err = std::fabs(got - r.length);
      EXPECT_LE(err, tol_of(tol) * r.polygon + 4e-16 * r.length) << "tol " << tol;
      size_t rows = gshgeom_measured_row_count(m);
      EXPECT_GE(rows, rows_before) << "a tighter tolerance never has fewer rows";
      rows_before = rows;
      if (tol == 0x1p-10) worst_coarse = err;
      if (tol == 0x1p-40) worst_fine = err;
      gshgeom_measured_destroy(m);
    }
    EXPECT_LE(worst_fine, 1e-12 * r.length) << "converged";
    EXPECT_LE(worst_fine, worst_coarse + 4e-16 * r.length) << "refinement does not hurt";
    gshgeom_path_destroy(p);
  }
}

TEST(AlongReference, ARefinementActuallyRefinesACurveThatNeedsIt) {
  // The loop and the cusp are the curves on which the coarse tolerance is a real
  // restriction: more rows, and a smaller error, as the tolerance falls.
  for (const char * name : {"cubic_cusp_third", "cubic_loop", "rational_cubic_heavy"}) {
    const AlongRef * r = ref_named(name);
    ASSERT_NE(r, nullptr) << name;
    GSHGEOM_Path * p = path_of(*r);
    GSHGEOM_Measured * coarse = measure(p, 0x1p-10);
    GSHGEOM_Measured * fine = measure(p, 0x1p-40);
    EXPECT_GT(gshgeom_measured_row_count(fine), gshgeom_measured_row_count(coarse)) << name;
    gshgeom_measured_destroy(coarse);
    gshgeom_measured_destroy(fine);
    gshgeom_path_destroy(p);
  }
}

TEST(AlongReference, TheChordSumOfAPolylineConvergesToTheReference) {
  // An independent view of convergence: the length of the polyline through 2^k
  // equal parameter steps of the curve (gshgeom_curve_eval) rises toward the
  // reference with the error falling by about a factor of four per doubling.
  const AlongRef * r = ref_named("cubic_s");
  ASSERT_NE(r, nullptr);
  GSHGEOM_Curve c;
  ASSERT_EQ(gshgeom_curve_cubic(2, P(r->p[0][0], r->p[0][1]), P(r->p[1][0], r->p[1][1]),
                P(r->p[2][0], r->p[2][1]), P(r->p[3][0], r->p[3][1]), &c),
      GSHGEOM_OK);
  double prev_err = 0;
  for (int k = 2; k <= 12; k++) {
    int n = 1 << k;
    double sum = 0;
    GSHGEOM_Point last = c.p[0];
    for (int i = 1; i <= n; i++) {
      GSHGEOM_Point q;
      ASSERT_EQ(gshgeom_curve_eval(&c, (double)i / n, &q), GSHGEOM_OK);
      sum += dist2(last, q);
      last = q;
    }
    double err = r->length - sum;
    EXPECT_GT(err, 0.0) << "a chord is shorter than its arc, k=" << k;
    if (k >= 4) {
      EXPECT_GT(prev_err / err, 3.5) << k;
      EXPECT_LT(prev_err / err, 4.5) << k;
    }
    prev_err = err;
  }
  EXPECT_LT(prev_err, 1e-5 * r->length);
  GSHGEOM_Measured * m = nullptr;
  GSHGEOM_Path * p = path_of(*r);
  m = measure(p);
  EXPECT_LT(std::fabs(gshgeom_measured_length(m) - r->length), 1e-12 * r->length);
  EXPECT_GT(r->length - prev_err, 0.0);
  gshgeom_measured_destroy(m);
  gshgeom_path_destroy(p);
}

TEST(AlongReference, PositionAndTangentAtFractionsOfTheLength) {
  for (double tol : {0x1p-20, 0.0}) {
    for (const AlongRef & r : kAlongRefs) {
      SCOPED_TRACE(r.name);
      GSHGEOM_Path * p = path_of(r);
      GSHGEOM_Measured * m = measure(p, tol);
      double scale = r.length;
      for (int i = 0; i < 3 + 1; i++) {
        for (int a = 0; a < 3; a++) {
          scale = std::max(scale, std::fabs(r.pt[i][a]));
        }
      }
      for (int i = 0; i < 4; i++) {
        GSHGEOM_Measured_Frame f = at(m, r.s[i]);
        GSHGEOM_Point want = P(r.pt[i][0], r.pt[i][1], r.pt[i][2]);
        double err = r.dim == 3 ? dist3(f.point, want) : dist2(f.point, want);
        // The point is on the curve at a distance within the table's error
        // of the asked one, and a curve is 1-Lipschitz in arc length.
        EXPECT_LE(err, 4 * tol_of(tol) * r.polygon + 8e-15 * scale) << i;
        double te = std::sqrt(std::pow(f.tangent.x - r.tan[i][0], 2) +
            std::pow(f.tangent.y - r.tan[i][1], 2) + std::pow(f.tangent.z - r.tan[i][2], 2));
        // The direction changes by (curvature * the position error).
        EXPECT_LE(te, 1e-8 + 8 * tol_of(tol) * r.polygon * 100 / r.length) << i;
        EXPECT_NEAR(std::sqrt(f.tangent.x * f.tangent.x + f.tangent.y * f.tangent.y +
                        f.tangent.z * f.tangent.z),
            1.0, 4e-16);
      }
      gshgeom_measured_destroy(m);
      gshgeom_path_destroy(p);
    }
  }
}

// ---------------------------------------------------------------------------
// Structure of the answers
// ---------------------------------------------------------------------------

// A connected path of every kind that ends where it began.
GSHGEOM_Path * every_kind(int dim) {
  GSHGEOM_Path * p = new_path(dim);
  auto q = [&](double x, double y, double z) { return P(x, y, dim == 3 ? z : 0.0); };
  const double w3[3] = {1.0, 0.7, 1.0};
  const double w4[4] = {1.0, 1.5, 0.6, 1.0};
  EXPECT_EQ(gshgeom_path_move_to(p, q(0, 0, 0)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_line_to(p, q(4, 0, 1)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_quadratic_to(p, q(6, 0, 2), q(6, 3, 0)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_rational_quadratic_to(p, q(8, 6, -1), q(10, 3, 1), w3), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_cubic_to(p, q(12, 0, 2), q(14, -4, 0), q(9, -4, 3)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_rational_cubic_to(p, q(5, -4, 1), q(3, -2, 0), q(0, 0, 0), w4), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_close(p), GSHGEOM_OK);
  return p;
}

TEST(AlongProperties, MonotoneContinuousAndOrthonormalAcrossAMultiKindPath) {
  for (int dim : {2, 3}) {
    SCOPED_TRACE(dim);
    GSHGEOM_Path * p = every_kind(dim);
    GSHGEOM_Measured * m = measure(p);
    ASSERT_EQ(gshgeom_measured_segment_count(m), 5u) << "the close has length zero";
    double length = gshgeom_measured_length(m);
    const int n = 20000;
    GSHGEOM_Measured_Frame prev;
    std::memset(&prev, 0, sizeof prev);
    double chords = 0;
    for (int i = 0; i < n; i++) {
      double s = i == n - 1 ? length : length * i / (n - 1);
      GSHGEOM_Measured_Frame f = at(m, s);
      ASSERT_TRUE(f.has_tangent) << i;
      EXPECT_EQ(f.has_normal != 0, dim == 2);
      double tn = std::sqrt(f.tangent.x * f.tangent.x + f.tangent.y * f.tangent.y +
          f.tangent.z * f.tangent.z);
      EXPECT_NEAR(tn, 1.0, 4e-16) << i;
      if (dim == 2) {
        double nn = std::hypot(f.normal.x, f.normal.y);
        EXPECT_NEAR(nn, 1.0, 4e-16) << i;
        EXPECT_EQ(f.tangent.x * f.normal.x + f.tangent.y * f.normal.y, 0.0) << "exactly orthogonal";
        GSHGEOM_Point right = gshgeom_normal_right(f.tangent);
        EXPECT_EQ(right.x, -f.normal.x);
        EXPECT_EQ(right.y, -f.normal.y);
      }
      if (i > 0) {
        double ds = length / (n - 1);
        double d = dim == 3 ? dist3(prev.point, f.point) : dist2(prev.point, f.point);
        EXPECT_LE(d, ds * (1 + 1e-9) + 1e-11) << "continuous, 1-Lipschitz: " << i;
        chords += d;
        EXPECT_GE(f.element, prev.element) << i;
        if (f.element == prev.element) {
          EXPECT_GE(f.t, prev.t - 1e-9) << i;
        }
      }
      prev = f;
    }
    EXPECT_LE(chords, length * (1 + 1e-9));
    EXPECT_GE(chords, length * (1 - 2e-3));
    gshgeom_measured_destroy(m);
    gshgeom_path_destroy(p);
  }
}

TEST(AlongProperties, TheTableIsConsistentWithTheLength) {
  GSHGEOM_Path * p = every_kind(2);
  GSHGEOM_Measured * m = measure(p);
  // The dump lists every row: s0 + len of one is the s0 of the next, and the
  // last ends at the length.
  std::string d = dump_of(m);
  std::vector<std::pair<double, double>> rows;
  size_t pos = 0;
  auto bits_to_double = [](const std::string & hex) {
    uint64_t v = std::strtoull(hex.c_str(), nullptr, 16);
    double x;
    std::memcpy(&x, &v, sizeof x);
    return x;
  };
  while ((pos = d.find("\nR ", pos)) != std::string::npos) {
    pos += 3;
    char seg[32], t0[32], t1[32], s0[32], len[32];
    ASSERT_EQ(std::sscanf(d.c_str() + pos, "%31s %31s %31s %31s %31s", seg, t0, t1, s0, len), 5);
    rows.push_back({bits_to_double(s0), bits_to_double(len)});
  }
  ASSERT_EQ(rows.size(), gshgeom_measured_row_count(m));
  double sum = 0;
  for (size_t i = 0; i < rows.size(); i++) {
    EXPECT_EQ(rows[i].first, sum) << i;
    EXPECT_GE(rows[i].second, 0.0);
    sum += rows[i].second;
  }
  EXPECT_EQ(sum, gshgeom_measured_length(m));
  gshgeom_measured_destroy(m);
  gshgeom_path_destroy(p);
}

// ---------------------------------------------------------------------------
// Degenerate input: coincident control points, cusps, zero length
// ---------------------------------------------------------------------------

GSHGEOM_Measured_Frame frame_of_curve(GSHGEOM_Path * p, double s, GSHGEOM_Measured ** keep = nullptr) {
  GSHGEOM_Measured * m = measure(p);
  GSHGEOM_Measured_Frame f = at(m, s < 0 ? gshgeom_measured_length(m) : s);
  if (keep) {
    *keep = m;
  } else {
    gshgeom_measured_destroy(m);
  }
  return f;
}

void expect_direction(const GSHGEOM_Measured_Frame & f, double x, double y, const char * what) {
  double n = std::hypot(x, y);
  EXPECT_TRUE(f.has_tangent) << what;
  EXPECT_NEAR(f.tangent.x, x / n, 1e-12) << what;
  EXPECT_NEAR(f.tangent.y, y / n, 1e-12) << what;
}

TEST(AlongDegenerate, ACoincidentControlPointAtAnEndStillHasADirection) {
  {
    // P0 = P1: B'(t) is about 6 t (P2 - P0) near 0.
    GSHGEOM_Path * p = new_path();
    ASSERT_EQ(gshgeom_path_move_to(p, P(0, 0)), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_path_cubic_to(p, P(0, 0), P(3, 1), P(4, 4)), GSHGEOM_OK);
    expect_direction(frame_of_curve(p, 0.0), 3, 1, "cubic P0=P1 at the start");
    // At arc length s the parameter is about sqrt(s / 9.5) and the direction
    // has moved by about that much: 6 t (3, 1) + 3 t^2 (1, 3).
    GSHGEOM_Measured_Frame near_start = frame_of_curve(p, 1e-14);
    EXPECT_NEAR(near_start.tangent.x, 3 / std::hypot(3, 1), 2e-7);
    EXPECT_NEAR(near_start.tangent.y, 1 / std::hypot(3, 1), 2e-7);
    gshgeom_path_destroy(p);
  }
  {
    // P0 = P1 = P2: the first nonzero derivative is the third.
    GSHGEOM_Path * p = new_path();
    ASSERT_EQ(gshgeom_path_move_to(p, P(0, 0)), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_path_cubic_to(p, P(0, 0), P(0, 0), P(3, 4)), GSHGEOM_OK);
    expect_direction(frame_of_curve(p, 0.0), 3, 4, "cubic P0=P1=P2");
    gshgeom_path_destroy(p);
  }
  {
    GSHGEOM_Path * p = new_path();
    ASSERT_EQ(gshgeom_path_move_to(p, P(0, 0)), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_path_quadratic_to(p, P(0, 0), P(3, 4)), GSHGEOM_OK);
    expect_direction(frame_of_curve(p, 0.0), 3, 4, "quadratic P0=P1");
    gshgeom_path_destroy(p);
  }
  {
    // P2 = P3: the curve arrives along P2 - P1 (the left-hand limit at t = 1).
    GSHGEOM_Path * p = new_path();
    ASSERT_EQ(gshgeom_path_move_to(p, P(0, 0)), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_path_cubic_to(p, P(1, 2), P(3, 1), P(3, 1)), GSHGEOM_OK);
    expect_direction(frame_of_curve(p, -1.0), 2, -1, "cubic P2=P3 at the end");
    gshgeom_path_destroy(p);
  }
  {
    // The same, rational: B = t^2 P2 / D along (3, 4).
    GSHGEOM_Path * p = new_path();
    const double w[3] = {1.0, 2.0, 1.0};
    ASSERT_EQ(gshgeom_path_move_to(p, P(0, 0)), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_path_rational_quadratic_to(p, P(0, 0), P(3, 4), w), GSHGEOM_OK);
    expect_direction(frame_of_curve(p, 0.0), 3, 4, "rational P0=P1");
    gshgeom_path_destroy(p);
  }
}

TEST(AlongDegenerate, ARationalCurveWithThreeCoincidentPointsFollowsTheThirdDerivative) {
  // P0 = P1 = P2 with unequal weights: B(t) - B(0) is of order t^3 and E = N'D -
  // ND' vanishes to second order at 0, so the direction comes from the second
  // derivative of E, whose terms (Leibniz, with the binomial weights) mix P0 and
  // P3 until they cancel to the direction P3 - P0.
  GSHGEOM_Path * p = new_path();
  const double w[4] = {1.0, 2.0, 0.5, 1.5};
  ASSERT_EQ(gshgeom_path_move_to(p, P(1, 2)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_rational_cubic_to(p, P(1, 2), P(1, 2), P(4, 6), w), GSHGEOM_OK);
  expect_direction(frame_of_curve(p, 0.0), 3, 4, "rational cubic P0=P1=P2");
  gshgeom_path_destroy(p);
  // The mirror image at the other end.
  p = new_path();
  ASSERT_EQ(gshgeom_path_move_to(p, P(1, 2)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_rational_cubic_to(p, P(4, 6), P(4, 6), P(4, 6), w), GSHGEOM_OK);
  expect_direction(frame_of_curve(p, -1.0), 3, 4, "rational cubic P1=P2=P3 at the end");
  gshgeom_path_destroy(p);
}

TEST(AlongReference, TheEndsOfACurveAreExact) {
  for (const AlongRef & r : kAlongRefs) {
    SCOPED_TRACE(r.name);
    GSHGEOM_Path * p = path_of(r);
    GSHGEOM_Measured * m = measure(p);
    GSHGEOM_Measured_Frame start = at(m, 0.0);
    GSHGEOM_Measured_Frame end = at(m, gshgeom_measured_length(m));
    int n = gshgeom_curve_point_count((GSHGEOM_Curve_Kind)r.kind);
    EXPECT_EQ(start.t, 0.0);
    EXPECT_EQ(end.t, 1.0);
    EXPECT_EQ(start.point.x, r.p[0][0]);
    EXPECT_EQ(start.point.y, r.p[0][1]);
    EXPECT_EQ(start.point.z, r.p[0][2]);
    EXPECT_EQ(end.point.x, r.p[n - 1][0]);
    EXPECT_EQ(end.point.y, r.p[n - 1][1]);
    EXPECT_EQ(end.point.z, r.p[n - 1][2]);
    gshgeom_measured_destroy(m);
    gshgeom_path_destroy(p);
  }
}

TEST(AlongDegenerate, ACuspFlipsTheTangentAndNeverGivesNan) {
  const AlongRef * r = ref_named("cubic_cusp");
  ASSERT_NE(r, nullptr);
  ASSERT_GT(r->cusp_s, 0.0);
  GSHGEOM_Path * p = path_of(*r);
  GSHGEOM_Measured * m = measure(p);
  GSHGEOM_Measured_Frame before = at(m, r->cusp_s - 1e-6);
  GSHGEOM_Measured_Frame after = at(m, r->cusp_s + 1e-6);
  ASSERT_TRUE(before.has_tangent);
  ASSERT_TRUE(after.has_tangent);
  EXPECT_LT(before.tangent.x * after.tangent.x + before.tangent.y * after.tangent.y, -0.99)
      << "a cusp reverses the direction of travel";
  EXPECT_LT(dist2(before.point, after.point), 3e-6);
  for (double d : {-1e-12, -1e-15, 0.0, 1e-15, 1e-12}) {
    GSHGEOM_Measured_Frame f = at(m, r->cusp_s + d);
    ASSERT_TRUE(f.has_tangent) << d;
    EXPECT_NEAR(std::hypot(f.tangent.x, f.tangent.y), 1.0, 4e-16) << d;
    EXPECT_NEAR(std::hypot(f.normal.x, f.normal.y), 1.0, 4e-16) << d;
    EXPECT_TRUE(std::isfinite(f.point.x) && std::isfinite(f.point.y)) << d;
  }
  gshgeom_measured_destroy(m);
  gshgeom_path_destroy(p);
}

TEST(AlongDegenerate, ZeroLengthPathsHaveAPointAndNoDirection) {
  {
    GSHGEOM_Path * p = new_path();
    ASSERT_EQ(gshgeom_path_move_to(p, P(2, 5)), GSHGEOM_OK);
    GSHGEOM_Measured * m = measure(p);
    EXPECT_EQ(gshgeom_measured_length(m), 0.0);
    EXPECT_EQ(gshgeom_measured_row_count(m), 0u);
    EXPECT_EQ(gshgeom_measured_subpath_count(m), 0u) << "a lone move is not a subpath";
    GSHGEOM_Measured_Frame f = at(m, 0.0);
    EXPECT_EQ(f.point.x, 2.0);
    EXPECT_EQ(f.point.y, 5.0);
    EXPECT_FALSE(f.has_tangent);
    EXPECT_FALSE(f.has_normal);
    EXPECT_EQ(f.tangent.x, 0.0);
    GSHGEOM_Measured_Frame g;
    EXPECT_EQ(gshgeom_measured_at(m, 1.0, &g), GSHGEOM_ERR_INVALID);
    ASSERT_EQ(gshgeom_measured_at_clamped(m, 1.0, &g), GSHGEOM_OK);
    EXPECT_EQ(g.point.x, 2.0);
    gshgeom_measured_destroy(m);
    gshgeom_path_destroy(p);
  }
  {
    // M p L p Z: every segment has length zero, and the subpath exists.
    GSHGEOM_Path * p = new_path();
    ASSERT_EQ(gshgeom_path_move_to(p, P(7, 7)), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_path_line_to(p, P(7, 7)), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_path_close(p), GSHGEOM_OK);
    GSHGEOM_Measured * m = measure(p);
    EXPECT_EQ(gshgeom_measured_length(m), 0.0);
    EXPECT_EQ(gshgeom_measured_segment_count(m), 0u);
    ASSERT_EQ(gshgeom_measured_subpath_count(m), 1u);
    GSHGEOM_Measured_Subpath sp;
    ASSERT_EQ(gshgeom_measured_subpath(m, 0, &sp), GSHGEOM_OK);
    EXPECT_EQ(sp.begin, 0.0);
    EXPECT_EQ(sp.end, 0.0);
    EXPECT_TRUE(sp.closed);
    GSHGEOM_Measured_Frame f = at(m, 0.0);
    EXPECT_EQ(f.point.x, 7.0);
    EXPECT_FALSE(f.has_tangent);
    gshgeom_measured_destroy(m);
    gshgeom_path_destroy(p);
  }
  {
    // A curve whose control points all coincide has no extent either.
    GSHGEOM_Path * p = new_path(3);
    ASSERT_EQ(gshgeom_path_move_to(p, P(1, 2, 3)), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_path_cubic_to(p, P(1, 2, 3), P(1, 2, 3), P(1, 2, 3)), GSHGEOM_OK);
    const double w[3] = {1.0, 2.0, 3.0};
    ASSERT_EQ(gshgeom_path_rational_quadratic_to(p, P(1, 2, 3), P(1, 2, 3), w), GSHGEOM_OK);
    GSHGEOM_Measured * m = measure(p);
    EXPECT_EQ(gshgeom_measured_length(m), 0.0);
    EXPECT_EQ(gshgeom_measured_segment_count(m), 0u);
    EXPECT_EQ(gshgeom_measured_row_count(m), 0u);
    EXPECT_EQ(at(m, 0.0).point.z, 3.0);
    gshgeom_measured_destroy(m);
    gshgeom_path_destroy(p);
  }
  {
    // The empty path has no point at all.
    GSHGEOM_Path * p = new_path();
    GSHGEOM_Measured * m = measure(p);
    EXPECT_EQ(gshgeom_measured_length(m), 0.0);
    GSHGEOM_Measured_Frame f;
    EXPECT_EQ(gshgeom_measured_at(m, 0.0, &f), GSHGEOM_ERR_INVALID);
    EXPECT_EQ(gshgeom_measured_at_clamped(m, 0.0, &f), GSHGEOM_ERR_INVALID);
    EXPECT_EQ(gshgeom_measured_subpath_count(m), 0u);
    gshgeom_measured_destroy(m);
    gshgeom_path_destroy(p);
  }
}

TEST(AlongDegenerate, ZeroLengthSegmentsBetweenOthersAreNeverSelected) {
  GSHGEOM_Path * p = new_path();
  ASSERT_EQ(gshgeom_path_move_to(p, P(0, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(p, P(0, 0)), GSHGEOM_OK); // element 1: zero
  ASSERT_EQ(gshgeom_path_line_to(p, P(3, 4)), GSHGEOM_OK); // element 2
  ASSERT_EQ(gshgeom_path_cubic_to(p, P(3, 4), P(3, 4), P(3, 4)), GSHGEOM_OK); // element 3: zero
  ASSERT_EQ(gshgeom_path_line_to(p, P(3, 4)), GSHGEOM_OK); // element 4: zero
  ASSERT_EQ(gshgeom_path_line_to(p, P(3, 10)), GSHGEOM_OK); // element 5
  GSHGEOM_Measured * m = measure(p);
  EXPECT_EQ(gshgeom_measured_length(m), 11.0);
  EXPECT_EQ(gshgeom_measured_segment_count(m), 2u);
  EXPECT_EQ(at(m, 0.0).element, 2u);
  EXPECT_EQ(at(m, 4.9).element, 2u);
  EXPECT_EQ(at(m, 5.0).element, 5u) << "the boundary belongs to the later segment";
  EXPECT_EQ(at(m, 5.0).point.y, 4.0);
  EXPECT_EQ(at(m, 5.0).tangent.y, 1.0);
  EXPECT_EQ(at(m, 11.0).element, 5u);
  gshgeom_measured_destroy(m);
  gshgeom_path_destroy(p);
}

// ---------------------------------------------------------------------------
// Subpaths
// ---------------------------------------------------------------------------

// The fixture of the independent generator (tools/gen-hash-vectors.py): lines
// with exact integer lengths, an open subpath, a closed one with a zero-length
// line inside it, then a new subpath started by a drawing verb after the close.
GSHGEOM_Path * subpath_fixture2() {
  GSHGEOM_Path * p = new_path();
  EXPECT_EQ(gshgeom_path_move_to(p, P(0, 0)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_line_to(p, P(3, 4)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_move_to(p, P(10, 10)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_line_to(p, P(13, 10)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_line_to(p, P(13, 10)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_line_to(p, P(13, 14)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_close(p), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_line_to(p, P(10, 10)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_line_to(p, P(10, 16)), GSHGEOM_OK);
  return p;
}

GSHGEOM_Path * subpath_fixture3() {
  GSHGEOM_Path * p = new_path(3);
  EXPECT_EQ(gshgeom_path_move_to(p, P(1, 2, 3)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_line_to(p, P(4, 6, 15)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_move_to(p, P(0, 0, 0)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_line_to(p, P(0, 0, 0)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_line_to(p, P(2, 3, 6)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_close(p), GSHGEOM_OK);
  return p;
}

TEST(AlongSubpaths, ADistanceRunsAcrossThemInOrder) {
  GSHGEOM_Path * p = subpath_fixture2();
  GSHGEOM_Measured * m = measure(p);
  EXPECT_EQ(gshgeom_measured_length(m), 23.0);
  EXPECT_EQ(gshgeom_measured_segment_count(m), 5u);
  ASSERT_EQ(gshgeom_measured_subpath_count(m), 3u);
  const double want[3][3] = {{0, 5, 0}, {5, 17, 1}, {17, 23, 0}};
  for (size_t i = 0; i < 3; i++) {
    GSHGEOM_Measured_Subpath sp;
    ASSERT_EQ(gshgeom_measured_subpath(m, i, &sp), GSHGEOM_OK);
    EXPECT_EQ(sp.begin, want[i][0]) << i;
    EXPECT_EQ(sp.end, want[i][1]) << i;
    EXPECT_EQ(sp.closed != 0, want[i][2] != 0) << i;
  }
  GSHGEOM_Measured_Subpath sp;
  EXPECT_EQ(gshgeom_measured_subpath(m, 3, &sp), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_measured_subpath(nullptr, 0, &sp), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_measured_subpath(m, 0, nullptr), GSHGEOM_ERR_INVALID);

  // Just before and at the boundary between two subpaths: the point jumps and
  // the distance does not; the boundary is the later one's start.
  GSHGEOM_Measured_Frame f = at(m, 5.0 - 1e-9);
  EXPECT_NEAR(f.point.x, 3.0, 1e-9);
  EXPECT_NEAR(f.point.y, 4.0, 1e-9);
  EXPECT_EQ(f.element, 1u);
  f = at(m, 5.0);
  EXPECT_EQ(f.point.x, 10.0);
  EXPECT_EQ(f.point.y, 10.0);
  EXPECT_EQ(f.element, 3u);
  EXPECT_EQ(f.tangent.x, 1.0);
  // Inside the closed subpath: the line of 4 (element 5), then the close (6).
  f = at(m, 5.0 + 3.0 + 2.0);
  EXPECT_EQ(f.point.x, 13.0);
  EXPECT_EQ(f.point.y, 12.0);
  EXPECT_EQ(f.element, 5u);
  f = at(m, 12.0 + 2.5);
  EXPECT_EQ(f.element, 6u) << "the closing line is a segment";
  EXPECT_NEAR(f.point.x, 13.0 - 1.5, 1e-14);
  EXPECT_NEAR(f.point.y, 14.0 - 2.0, 1e-14);
  // The close ended at (10, 10), where a drawing verb starts a new subpath.
  f = at(m, 17.0);
  EXPECT_EQ(f.point.x, 10.0);
  EXPECT_EQ(f.point.y, 10.0);
  EXPECT_EQ(f.element, 8u);
  EXPECT_EQ(f.tangent.y, 1.0);
  f = at(m, 23.0);
  EXPECT_EQ(f.point.y, 16.0);
  gshgeom_measured_destroy(m);
  gshgeom_path_destroy(p);
}

// ---------------------------------------------------------------------------
// Out of range and arguments
// ---------------------------------------------------------------------------

TEST(AlongRange, ARefusedDistanceWritesNothing) {
  GSHGEOM_Path * p = new_path();
  ASSERT_EQ(gshgeom_path_move_to(p, P(0, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(p, P(3, 4)), GSHGEOM_OK);
  GSHGEOM_Measured * m = measure(p);
  GSHGEOM_Measured_Frame f, same;
  std::memset(&f, 0xAB, sizeof f);
  std::memset(&same, 0xAB, sizeof same);
  const double inf = HUGE_VAL;
  for (double s : {-1e-9, -1.0, 5.0 + 1e-9, 5.0000000000000009, 100.0, -inf, inf, std::nan("")}) {
    EXPECT_EQ(gshgeom_measured_at(m, s, &f), GSHGEOM_ERR_INVALID) << s;
    EXPECT_EQ(std::memcmp(&f, &same, sizeof f), 0) << "untouched";
  }
  EXPECT_EQ(gshgeom_measured_at(nullptr, 0.0, &f), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_measured_at(m, 0.0, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_measured_at_clamped(nullptr, 0.0, &f), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_measured_at_clamped(m, 0.0, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_measured_at(m, 0.0, &f), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_measured_at(m, 5.0, &f), GSHGEOM_OK) << "the length itself is in range";
  // Clamping: finite distances only.
  for (double s : {-1e-9, -1e300, -0.0}) {
    ASSERT_EQ(gshgeom_measured_at_clamped(m, s, &f), GSHGEOM_OK) << s;
    EXPECT_EQ(f.point.x, 0.0);
    EXPECT_EQ(f.point.y, 0.0);
    EXPECT_EQ(f.t, 0.0);
  }
  for (double s : {5.0 + 1e-9, 1e300, 5.0}) {
    ASSERT_EQ(gshgeom_measured_at_clamped(m, s, &f), GSHGEOM_OK) << s;
    EXPECT_EQ(f.point.x, 3.0);
    EXPECT_EQ(f.point.y, 4.0);
    EXPECT_EQ(f.t, 1.0);
  }
  std::memset(&f, 0xAB, sizeof f);
  for (double s : {-inf, inf, std::nan("")}) {
    EXPECT_EQ(gshgeom_measured_at_clamped(m, s, &f), GSHGEOM_ERR_INVALID) << s;
    EXPECT_EQ(std::memcmp(&f, &same, sizeof f), 0);
  }
  gshgeom_measured_destroy(m);
  gshgeom_path_destroy(p);
}

TEST(AlongRange, TheNullAndEmptyAccessors) {
  EXPECT_EQ(gshgeom_measured_length(nullptr), 0.0);
  EXPECT_EQ(gshgeom_measured_dimension(nullptr), 0);
  EXPECT_EQ(gshgeom_measured_segment_count(nullptr), 0u);
  EXPECT_EQ(gshgeom_measured_row_count(nullptr), 0u);
  EXPECT_EQ(gshgeom_measured_subpath_count(nullptr), 0u);
  gshgeom_measured_destroy(nullptr);
  GSHGEOM_Measured * m = nullptr;
  EXPECT_EQ(gshgeom_measured_create(nullptr, nullptr, nullptr, 0.0, &m), GSHGEOM_ERR_INVALID);
  GSHGEOM_Path * p = new_path();
  EXPECT_EQ(gshgeom_measured_create(nullptr, nullptr, p, 0.0, nullptr), GSHGEOM_ERR_INVALID);
  gshgeom_path_destroy(p);
}

TEST(AlongTolerance, TheArgumentAndItsBounds) {
  GSHGEOM_Path * p = new_path();
  ASSERT_EQ(gshgeom_path_move_to(p, P(0, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_quadratic_to(p, P(1, 2), P(3, 0)), GSHGEOM_OK);
  const double inf = HUGE_VAL;
  struct Case {
    double tol;
    GSHGEOM_Result want;
  };
  const Case cases[] = {
      {0.0, GSHGEOM_OK},
      {-0.0, GSHGEOM_OK},
      {GSHGEOM_MEASURE_TOLERANCE_MIN, GSHGEOM_OK},
      {std::nextafter(GSHGEOM_MEASURE_TOLERANCE_MIN, 0.0), GSHGEOM_ERR_LIMIT},
      {1e-300, GSHGEOM_ERR_LIMIT},
      {GSHGEOM_MEASURE_TOLERANCE_MAX, GSHGEOM_OK},
      {std::nextafter(GSHGEOM_MEASURE_TOLERANCE_MAX, 1.0), GSHGEOM_ERR_LIMIT},
      {1.0, GSHGEOM_ERR_LIMIT},
      {1e300, GSHGEOM_ERR_LIMIT},
      {-1e-12, GSHGEOM_ERR_INVALID},
      {-1.0, GSHGEOM_ERR_INVALID},
      {inf, GSHGEOM_ERR_INVALID},
      {-inf, GSHGEOM_ERR_INVALID},
      {std::nan(""), GSHGEOM_ERR_INVALID},
  };
  GSHGEOM_Curve curve;
  ASSERT_EQ(gshgeom_curve_quadratic(2, P(0, 0), P(1, 2), P(3, 0), &curve), GSHGEOM_OK);
  for (const Case & c : cases) {
    GSHGEOM_Measured * m = reinterpret_cast<GSHGEOM_Measured *>(0x1);
    GSHGEOM_Result r = gshgeom_measured_create(nullptr, nullptr, p, c.tol, &m);
    EXPECT_EQ(r, c.want) << c.tol;
    if (r == GSHGEOM_OK) {
      gshgeom_measured_destroy(m);
    } else {
      EXPECT_EQ(m, reinterpret_cast<GSHGEOM_Measured *>(0x1)) << "untouched";
    }
    double len = -5.0;
    EXPECT_EQ(gshgeom_measured_curve_length(&curve, c.tol, &len), c.want) << c.tol;
    if (c.want != GSHGEOM_OK) {
      EXPECT_EQ(len, -5.0);
    } else {
      EXPECT_GT(len, 0.0);
    }
  }
  // 0 means the default: the same table.
  GSHGEOM_Measured *a = measure(p, 0.0), *b = measure(p, GSHGEOM_MEASURE_TOLERANCE_DEFAULT);
  EXPECT_EQ(gshgeom_measured_hash(a), gshgeom_measured_hash(b));
  EXPECT_EQ(GSHGEOM_MEASURE_TOLERANCE_DEFAULT, 0x1p-40);
  EXPECT_EQ(GSHGEOM_MEASURE_TOLERANCE_MIN, 0x1p-46);
  EXPECT_EQ(GSHGEOM_MEASURE_TOLERANCE_MAX, 0x1p-4);
  gshgeom_measured_destroy(a);
  gshgeom_measured_destroy(b);
  gshgeom_path_destroy(p);
}

// ---------------------------------------------------------------------------
// The single-curve length
// ---------------------------------------------------------------------------

TEST(AlongCurveLength, MatchesTheReferenceAndTheMeasuredPath) {
  for (const AlongRef & r : kAlongRefs) {
    SCOPED_TRACE(r.name);
    GSHGEOM_Curve c;
    c.kind = (GSHGEOM_Curve_Kind)r.kind;
    c.dim = r.dim;
    std::memset(c.p, 0, sizeof c.p);
    std::memset(c.w, 0, sizeof c.w);
    int n = gshgeom_curve_point_count(c.kind);
    for (int i = 0; i < n; i++) {
      c.p[i] = P(r.p[i][0], r.p[i][1], r.p[i][2]);
      c.w[i] = r.w[i];
    }
    double len = 0;
    ASSERT_EQ(gshgeom_measured_curve_length(&c, 0.0, &len), GSHGEOM_OK);
    EXPECT_LE(std::fabs(len - r.length), GSHGEOM_MEASURE_TOLERANCE_DEFAULT * r.polygon + 4e-16 * r.length);
    GSHGEOM_Path * p = path_of(r);
    GSHGEOM_Measured * m = measure(p);
    EXPECT_EQ(len, gshgeom_measured_length(m)) << "one segment: the same rule, the same sum";
    gshgeom_measured_destroy(m);
    gshgeom_path_destroy(p);
  }
}

TEST(AlongCurveLength, CollinearCurvesAreIntegratedExactly) {
  // The speed of a curve along a line, monotone in its parameter, is a
  // polynomial of degree at most two: the eight-point rule integrates it
  // exactly, so the length is the distance between the ends.
  struct Case {
    GSHGEOM_Curve_Kind kind;
    double x[4];
  };
  const Case cases[] = {
      {GSHGEOM_CURVE_QUADRATIC, {0, 1, 4, 0}},
      {GSHGEOM_CURVE_QUADRATIC, {0, 2, 3, 0}},
      {GSHGEOM_CURVE_CUBIC, {0, 1, 3, 4}},
      {GSHGEOM_CURVE_CUBIC, {0, 0, 1, 4}},
      {GSHGEOM_CURVE_CUBIC, {0, 3, 3, 4}},
      {GSHGEOM_CURVE_CUBIC, {-2, -1, 5, 7}},
  };
  for (const Case & k : cases) {
    GSHGEOM_Curve c;
    std::memset(&c, 0, sizeof c);
    c.kind = k.kind;
    c.dim = 3;
    int n = gshgeom_curve_point_count(k.kind);
    for (int i = 0; i < n; i++) {
      c.p[i] = P(k.x[i], 2 * k.x[i], -2 * k.x[i]);
    }
    double len = 0;
    ASSERT_EQ(gshgeom_measured_curve_length(&c, 0.0, &len), GSHGEOM_OK);
    double want = 3.0 * std::fabs(k.x[n - 1] - k.x[0]);
    EXPECT_NEAR(len, want, 8e-15 * want) << (int)k.kind << " " << k.x[1];
  }
}

TEST(AlongCurveLength, RefusesWhatTheCurveCheckRefuses) {
  GSHGEOM_Curve c;
  ASSERT_EQ(gshgeom_curve_cubic(2, P(0, 0), P(1, 2), P(3, 2), P(4, 0), &c), GSHGEOM_OK);
  double len = -1;
  EXPECT_EQ(gshgeom_measured_curve_length(nullptr, 0.0, &len), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_measured_curve_length(&c, 0.0, nullptr), GSHGEOM_ERR_INVALID);
  GSHGEOM_Curve bad = c;
  bad.p[1].x = std::nan("");
  EXPECT_EQ(gshgeom_measured_curve_length(&bad, 0.0, &len), GSHGEOM_ERR_INVALID);
  bad = c;
  bad.p[2].y = HUGE_VAL;
  EXPECT_EQ(gshgeom_measured_curve_length(&bad, 0.0, &len), GSHGEOM_ERR_INVALID);
  bad = c;
  bad.p[2].y = 0x1p101;
  EXPECT_EQ(gshgeom_measured_curve_length(&bad, 0.0, &len), GSHGEOM_ERR_LIMIT);
  bad = c;
  bad.kind = (GSHGEOM_Curve_Kind)9;
  EXPECT_EQ(gshgeom_measured_curve_length(&bad, 0.0, &len), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(len, -1.0);
  // A curve with no extent has length zero.
  GSHGEOM_Curve point;
  ASSERT_EQ(gshgeom_curve_cubic(3, P(1, 1, 1), P(1, 1, 1), P(1, 1, 1), P(1, 1, 1), &point), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_measured_curve_length(&point, 0.0, &len), GSHGEOM_OK);
  EXPECT_EQ(len, 0.0);
}

// ---------------------------------------------------------------------------
// Hostile input
// ---------------------------------------------------------------------------

TEST(AlongHostile, HugeAndTinyCoordinatesStayFinite) {
  const double big = 0x1p100;
  {
    GSHGEOM_Path * p = new_path(3);
    ASSERT_EQ(gshgeom_path_move_to(p, P(-big, -big, -big)), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_path_line_to(p, P(big, big, big)), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_path_cubic_to(p, P(-big, big, big), P(big, -big, big), P(big, big, -big)), GSHGEOM_OK);
    GSHGEOM_Measured * m = measure(p);
    double length = gshgeom_measured_length(m);
    EXPECT_TRUE(std::isfinite(length));
    EXPECT_GE(length, std::sqrt(12.0) * big * 2) << "at least the line";
    GSHGEOM_Measured_Frame f = at(m, length / 2);
    EXPECT_TRUE(std::isfinite(f.point.x) && std::isfinite(f.point.y) && std::isfinite(f.point.z));
    EXPECT_NEAR(std::sqrt(f.tangent.x * f.tangent.x + f.tangent.y * f.tangent.y + f.tangent.z * f.tangent.z), 1.0, 4e-16);
    f = at(m, length);
    EXPECT_TRUE(std::isfinite(f.point.x));
    gshgeom_measured_destroy(m);
    gshgeom_path_destroy(p);
  }
  {
    // Subnormal-sized geometry: a length that is tiny or zero, never NaN, and
    // a query that does not trap.
    GSHGEOM_Path * p = new_path();
    const double tiny = 0x1p-1060;
    ASSERT_EQ(gshgeom_path_move_to(p, P(0, 0)), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_path_cubic_to(p, P(tiny, 0), P(tiny, tiny), P(0, 2 * tiny)), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_path_line_to(p, P(0, 3 * tiny)), GSHGEOM_OK);
    GSHGEOM_Measured * m = measure(p);
    double length = gshgeom_measured_length(m);
    EXPECT_TRUE(std::isfinite(length));
    GSHGEOM_Measured_Frame f;
    ASSERT_EQ(gshgeom_measured_at_clamped(m, length, &f), GSHGEOM_OK);
    EXPECT_TRUE(std::isfinite(f.point.x) && std::isfinite(f.point.y));
    gshgeom_measured_destroy(m);
    gshgeom_path_destroy(p);
  }
}

TEST(AlongHostile, ExtremeWeightsFinishFiniteWithinTheRowCap) {
  Rng rng{UINT64_C(0x1234567887654321)};
  const double weights[] = {0x1p-100, 0x1p-50, 1.0, 0x1p50, 0x1p100};
  GSHGEOM_Limits limits;
  gshgeom_limits_default(&limits);
  limits.max_points = 4000;
  int measured = 0, refused = 0;
  for (int i = 0; i < 120; i++) {
    int dim = (i & 1) ? 3 : 2;
    GSHGEOM_Path * p = new_path(dim);
    double scale = i % 3 == 0 ? 0x1p100 : (i % 3 == 1 ? 1.0 : 1e-30);
    auto pt = [&]() { return P(rng.range(-1, 1) * scale, rng.range(-1, 1) * scale, dim == 3 ? rng.range(-1, 1) * scale : 0.0); };
    ASSERT_EQ(gshgeom_path_move_to(p, pt()), GSHGEOM_OK);
    if (i % 2 == 0) {
      double w[3] = {weights[rng.next() % 5], weights[rng.next() % 5], weights[rng.next() % 5]};
      GSHGEOM_Point c = pt();
      ASSERT_EQ(gshgeom_path_rational_quadratic_to(p, c, pt(), w), GSHGEOM_OK);
    } else {
      double w[4] = {weights[rng.next() % 5], weights[rng.next() % 5], weights[rng.next() % 5], weights[rng.next() % 5]};
      GSHGEOM_Point c1 = pt(), c2 = pt();
      ASSERT_EQ(gshgeom_path_rational_cubic_to(p, c1, c2, pt(), w), GSHGEOM_OK);
    }
    GSHGEOM_Measured * m = nullptr;
    GSHGEOM_Result r = gshgeom_measured_create(nullptr, &limits, p, 0.0, &m);
    if (r == GSHGEOM_OK) {
      measured++;
      double length = gshgeom_measured_length(m);
      EXPECT_TRUE(std::isfinite(length)) << i;
      EXPECT_LE(gshgeom_measured_row_count(m), limits.max_points);
      for (int k = 0; k <= 8; k++) {
        GSHGEOM_Measured_Frame f;
        ASSERT_EQ(gshgeom_measured_at_clamped(m, length * k / 8.0, &f), GSHGEOM_OK) << i << " " << k;
        EXPECT_TRUE(std::isfinite(f.point.x) && std::isfinite(f.point.y) && std::isfinite(f.point.z)) << i;
        EXPECT_TRUE(std::isfinite(f.t)) << i;
        if (f.has_tangent) {
          EXPECT_NEAR(std::sqrt(f.tangent.x * f.tangent.x + f.tangent.y * f.tangent.y +
                          f.tangent.z * f.tangent.z),
              1.0, 1e-15) << i;
        }
      }
      gshgeom_measured_destroy(m);
    } else {
      refused++;
      EXPECT_EQ(r, GSHGEOM_ERR_LIMIT) << i;
      EXPECT_EQ(m, nullptr);
    }
    gshgeom_path_destroy(p);
  }
  EXPECT_GT(measured, 60) << "most extreme curves are measurable within the cap";
  (void)refused;
}

TEST(AlongHostile, RandomPathsOfEveryKindObeyTheInvariants) {
  Rng rng{UINT64_C(0xfeedfacecafebeef)};
  for (int i = 0; i < 150; i++) {
    int dim = (i % 3 == 0) ? 3 : 2;
    GSHGEOM_Path * p = new_path(dim);
    auto pt = [&]() { return P(rng.range(-50, 50), rng.range(-50, 50), dim == 3 ? rng.range(-50, 50) : 0.0); };
    ASSERT_EQ(gshgeom_path_move_to(p, pt()), GSHGEOM_OK);
    int n = 1 + (int)(rng.next() % 6);
    for (int k = 0; k < n; k++) {
      switch (rng.next() % 7) {
      case 0: ASSERT_EQ(gshgeom_path_line_to(p, pt()), GSHGEOM_OK); break;
      case 1: { GSHGEOM_Point c = pt(); ASSERT_EQ(gshgeom_path_quadratic_to(p, c, pt()), GSHGEOM_OK); break; }
      case 2: { GSHGEOM_Point a = pt(), b = pt(); ASSERT_EQ(gshgeom_path_cubic_to(p, a, b, pt()), GSHGEOM_OK); break; }
      case 3: {
        double w[3] = {rng.range(0.1, 5), rng.range(0.1, 5), rng.range(0.1, 5)};
        GSHGEOM_Point c = pt();
        ASSERT_EQ(gshgeom_path_rational_quadratic_to(p, c, pt(), w), GSHGEOM_OK);
        break;
      }
      case 4: {
        double w[4] = {rng.range(0.1, 5), rng.range(0.1, 5), rng.range(0.1, 5), rng.range(0.1, 5)};
        GSHGEOM_Point a = pt(), b = pt();
        ASSERT_EQ(gshgeom_path_rational_cubic_to(p, a, b, pt(), w), GSHGEOM_OK);
        break;
      }
      case 5: ASSERT_EQ(gshgeom_path_move_to(p, pt()), GSHGEOM_OK); break;
      default: ASSERT_EQ(gshgeom_path_close(p), GSHGEOM_OK); break;
      }
    }
    GSHGEOM_Measured * m = measure(p);
    double length = gshgeom_measured_length(m);
    EXPECT_GE(length, 0.0);
    double prev_s0 = 0.0;
    for (size_t k = 0; k < gshgeom_measured_subpath_count(m); k++) {
      GSHGEOM_Measured_Subpath sp;
      ASSERT_EQ(gshgeom_measured_subpath(m, k, &sp), GSHGEOM_OK);
      EXPECT_GE(sp.begin, prev_s0);
      EXPECT_GE(sp.end, sp.begin);
      EXPECT_LE(sp.end, length);
      prev_s0 = sp.end;
    }
    GSHGEOM_Measured_Frame last;
    std::memset(&last, 0, sizeof last);
    for (int k = 0; k <= 64; k++) {
      GSHGEOM_Measured_Frame f;
      ASSERT_EQ(gshgeom_measured_at_clamped(m, length * k / 64.0, &f), GSHGEOM_OK);
      EXPECT_TRUE(std::isfinite(f.point.x) && std::isfinite(f.point.y) && std::isfinite(f.point.z));
      if (f.has_tangent) {
        EXPECT_NEAR(std::sqrt(f.tangent.x * f.tangent.x + f.tangent.y * f.tangent.y + f.tangent.z * f.tangent.z), 1.0, 1e-15);
      }
      if (k > 0) {
        EXPECT_GE(f.element, last.element);
      }
      last = f;
    }
    gshgeom_measured_destroy(m);
    gshgeom_path_destroy(p);
  }
}

// ---------------------------------------------------------------------------
// Limits and allocation
// ---------------------------------------------------------------------------

GSHGEOM_Path * mixed_path() {
  GSHGEOM_Path * p = new_path();
  EXPECT_EQ(gshgeom_path_move_to(p, P(0, 0)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_line_to(p, P(4, 0)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_quadratic_to(p, P(6, 0), P(6, 3)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_cubic_to(p, P(8, 8), P(-3, 9), P(1, 0.5)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_move_to(p, P(20, 20)), GSHGEOM_OK);
  const double w[3] = {1.0, 2.0, 1.0};
  EXPECT_EQ(gshgeom_path_rational_quadratic_to(p, P(22, 25), P(30, 20), w), GSHGEOM_OK);
  return p;
}

TEST(AlongLimits, MaxPointsBoundsTheRowsExactly) {
  GSHGEOM_Path * p = mixed_path();
  GSHGEOM_Measured * all = measure(p);
  size_t rows = gshgeom_measured_row_count(all);
  ASSERT_GT(rows, 8u);
  gshgeom_measured_destroy(all);
  GSHGEOM_Limits l;
  gshgeom_limits_default(&l);
  l.max_points = rows;
  GSHGEOM_Measured * ok = measure(p, 0.0, &l);
  EXPECT_EQ(gshgeom_measured_row_count(ok), rows);
  gshgeom_measured_destroy(ok);
  for (size_t cap : {rows - 1, rows - 2, (size_t)1}) {
    l.max_points = cap;
    GSHGEOM_Measured * m = reinterpret_cast<GSHGEOM_Measured *>(0x1);
    EXPECT_EQ(gshgeom_measured_create(nullptr, &l, p, 0.0, &m), GSHGEOM_ERR_LIMIT) << cap;
    EXPECT_EQ(m, reinterpret_cast<GSHGEOM_Measured *>(0x1)) << "untouched";
  }
  l.max_points = 0;
  GSHGEOM_Measured * d = measure(p, 0.0, &l);
  EXPECT_EQ(gshgeom_measured_row_count(d), rows) << "0 means the default";
  gshgeom_measured_destroy(d);
  EXPECT_EQ(GSHGEOM_LIMITS_DEFAULT_MAX_POINTS, (size_t)1 << 22);
  gshgeom_path_destroy(p);
}

TEST(AlongAllocation, AFailedAllocationIsOomAndLeaksNothing) {
  GSHGEOM_Path * p = mixed_path();
  tracker::Tracker count;
  {
    GSHGEOM_Allocator a = tracker::allocator_for(&count);
    GSHGEOM_Measured * m = nullptr;
    ASSERT_EQ(gshgeom_measured_create(&a, nullptr, p, 0.0, &m), GSHGEOM_OK);
    gshgeom_measured_destroy(m);
    EXPECT_TRUE(count.live.empty());
  }
  ASSERT_GT(count.calls, 4);
  for (long k = 0; k < count.calls; k++) {
    for (int persistent = 0; persistent < 2; persistent++) {
      tracker::Tracker t;
      t.fail_from = k;
      t.persistent = persistent != 0;
      GSHGEOM_Allocator a = tracker::allocator_for(&t);
      GSHGEOM_Measured * m = nullptr;
      EXPECT_EQ(gshgeom_measured_create(&a, nullptr, p, 0.0, &m), GSHGEOM_ERR_OOM) << k;
      EXPECT_EQ(m, nullptr);
      EXPECT_TRUE(t.live.empty()) << "k=" << k << " leaked";
    }
  }
  gshgeom_path_destroy(p);
}

TEST(AlongAllocation, ALongPathGrowsEveryArrayThroughSeveralDoublings) {
  // More than 16 segments, rows and subpaths, so that each array grows twice or
  // more and a failure falls on a later doubling.
  GSHGEOM_Path * p = new_path();
  for (int i = 0; i < 40; i++) {
    ASSERT_EQ(gshgeom_path_move_to(p, P(i * 10.0, 0)), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_path_cubic_to(p, P(i * 10.0 + 3, 5), P(i * 10.0 + 6, -5), P(i * 10.0 + 9, 0)), GSHGEOM_OK);
  }
  tracker::Tracker count;
  {
    GSHGEOM_Allocator a = tracker::allocator_for(&count);
    GSHGEOM_Measured * m = nullptr;
    ASSERT_EQ(gshgeom_measured_create(&a, nullptr, p, 0.0, &m), GSHGEOM_OK);
    EXPECT_EQ(gshgeom_measured_subpath_count(m), 40u);
    EXPECT_GT(gshgeom_measured_row_count(m), 40u);
    gshgeom_measured_destroy(m);
  }
  for (long k = 0; k < count.calls; k++) {
    tracker::Tracker t;
    t.fail_from = k;
    GSHGEOM_Allocator a = tracker::allocator_for(&t);
    GSHGEOM_Measured * m = nullptr;
    EXPECT_EQ(gshgeom_measured_create(&a, nullptr, p, 0.0, &m), GSHGEOM_ERR_OOM) << k;
    EXPECT_TRUE(t.live.empty()) << k;
  }
  gshgeom_path_destroy(p);
}

// ---------------------------------------------------------------------------
// Hash and dump
// ---------------------------------------------------------------------------

TEST(AlongHash, LinesMatchTheIndependentImplementation) {
  GSHGEOM_Path * p2 = subpath_fixture2();
  GSHGEOM_Path * p3 = subpath_fixture3();
  GSHGEOM_Measured * m2 = measure(p2);
  GSHGEOM_Measured * m3 = measure(p3);
  EXPECT_EQ(gshgeom_measured_hash(m2), kMeasuredHash2);
  EXPECT_EQ(gshgeom_measured_hash(m3), kMeasuredHash3);
  EXPECT_EQ(dump_of(m2), kMeasuredDump2);
  EXPECT_EQ(dump_of(m3), kMeasuredDump3);
  EXPECT_NE(kMeasuredHash2, kMeasuredHash3);
  gshgeom_measured_destroy(m2);
  gshgeom_measured_destroy(m3);
  gshgeom_path_destroy(p2);
  gshgeom_path_destroy(p3);
  EXPECT_EQ(gshgeom_measured_hash(nullptr), UINT64_C(0xcbf29ce484222325));
}

// FNV-1a 64 as this file states it, to hash what a dump says.
struct Fnv {
  uint64_t h = UINT64_C(0xcbf29ce484222325);
  void byte(uint8_t b) {
    h ^= b;
    h *= UINT64_C(0x100000001b3);
  }
  void le(uint64_t v, int n) {
    for (int i = 0; i < n; i++) byte((uint8_t)((v >> (8 * i)) & 0xff));
  }
};

// Parse a dump and hash its fields by the spelling in along.h: what the dump
// says is what the hash covers.
uint64_t hash_of_dump(const std::string & dump, uint64_t * header_hash) {
  std::vector<std::string> lines;
  size_t pos = 0;
  while (pos < dump.size()) {
    size_t e = dump.find('\n', pos);
    lines.push_back(dump.substr(pos, e - pos));
    pos = e + 1;
  }
  unsigned version = 0;
  int dim = 0;
  unsigned long long nseg = 0, nrows = 0, nsub = 0, length = 0, hash = 0;
  EXPECT_EQ(std::sscanf(lines[0].c_str(),
                "gshgeom.measured %u dim=%d segments=%llu rows=%llu subpaths=%llu length=%llx hash=%llx",
                &version, &dim, &nseg, &nrows, &nsub, &length, &hash),
      7);
  *header_hash = hash;
  Fnv f;
  f.le(version, 4);
  f.le(0x5341454du, 4);
  f.le((uint64_t)dim, 4);
  f.le(nseg, 8);
  f.le(nrows, 8);
  f.le(nsub, 8);
  auto words = [](const std::string & s) {
    std::vector<std::string> w;
    size_t i = 0;
    while (i < s.size()) {
      size_t j = s.find(' ', i);
      if (j == std::string::npos) j = s.size();
      w.push_back(s.substr(i, j - i));
      i = j + 1;
    }
    return w;
  };
  auto hex64 = [&](const std::string & s) { f.le(std::strtoull(s.c_str(), nullptr, 16), 8); };
  std::vector<std::string> o = words(lines[1]);
  EXPECT_EQ(o[0], "O");
  for (size_t i = 1; i < o.size(); i++) hex64(o[i]);
  f.le(length, 8);
  size_t at = 2;
  for (unsigned long long i = 0; i < nseg; i++, at++) {
    std::vector<std::string> w = words(lines[at]);
    EXPECT_EQ(w[0], "G");
    f.le(std::strtoull(w[1].c_str(), nullptr, 10), 8);
    f.byte((uint8_t)std::strtoul(w[2].c_str(), nullptr, 10));
    for (size_t k = 3; k < w.size(); k++) hex64(w[k]);
  }
  for (unsigned long long i = 0; i < nrows; i++, at++) {
    std::vector<std::string> w = words(lines[at]);
    EXPECT_EQ(w[0], "R");
    f.le(std::strtoull(w[1].c_str(), nullptr, 10), 8);
    for (size_t k = 2; k < w.size(); k++) hex64(w[k]);
  }
  for (unsigned long long i = 0; i < nsub; i++, at++) {
    std::vector<std::string> w = words(lines[at]);
    EXPECT_EQ(w[0], "S");
    hex64(w[1]);
    hex64(w[2]);
    f.byte((uint8_t)std::strtoul(w[3].c_str(), nullptr, 10));
  }
  EXPECT_EQ(at, lines.size()) << "every line of the dump is accounted for";
  return f.h;
}

TEST(AlongHash, CurvesHashWhatTheirDumpSays) {
  std::vector<GSHGEOM_Path *> paths;
  for (const AlongRef & r : kAlongRefs) paths.push_back(path_of(r));
  paths.push_back(every_kind(2));
  paths.push_back(every_kind(3));
  paths.push_back(mixed_path());
  std::set<uint64_t> seen;
  for (GSHGEOM_Path * p : paths) {
    GSHGEOM_Measured * m = measure(p);
    std::string d = dump_of(m);
    uint64_t header = 0;
    uint64_t recomputed = hash_of_dump(d, &header);
    EXPECT_EQ(header, gshgeom_measured_hash(m));
    EXPECT_EQ(recomputed, gshgeom_measured_hash(m)) << d.substr(0, 80);
    EXPECT_TRUE(seen.insert(recomputed).second) << "no two fixtures share a hash";
    gshgeom_measured_destroy(m);
    gshgeom_path_destroy(p);
  }
}

TEST(AlongHash, EveryInputThatChangesTheTableChangesTheHash) {
  const AlongRef * r = ref_named("cubic_s");
  GSHGEOM_Path * p = path_of(*r);
  GSHGEOM_Measured * base = measure(p);
  uint64_t h = gshgeom_measured_hash(base);
  GSHGEOM_Measured * tight = measure(p, 0x1p-30);
  GSHGEOM_Measured * loose = measure(p, 0x1p-10);
  EXPECT_NE(gshgeom_measured_hash(tight), h) << "a different table is a different hash";
  EXPECT_NE(gshgeom_measured_hash(loose), h);
  AlongRef moved = *r;
  moved.p[2][1] = std::nextafter(moved.p[2][1], 10.0);
  GSHGEOM_Path * q = path_of(moved);
  GSHGEOM_Measured * other = measure(q);
  EXPECT_NE(gshgeom_measured_hash(other), h) << "one ulp of one control point";
  GSHGEOM_Measured * again = measure(p);
  EXPECT_EQ(gshgeom_measured_hash(again), h) << "deterministic";
  for (GSHGEOM_Measured * m : {base, tight, loose, other, again}) gshgeom_measured_destroy(m);
  gshgeom_path_destroy(p);
  gshgeom_path_destroy(q);
}

TEST(AlongHash, TheCurvedFixturesHaveTheirPinnedHashes) {
  // Regression pins: the values were produced by this library and are the same
  // on every target by the cross-target probe; the independent checks are the
  // dump round trip above and the high-precision reference.  The coarse
  // tolerances put the acceptance decision of the adaptive rule at stake: a
  // change to how a panel is judged changes which panels there are.
  struct Pin {
    const char * name;
    double tol;
    size_t rows;
    uint64_t hash;
  };
  const Pin pins[] = {
      {"cubic_s", 0.0, 20u, UINT64_C(0x4f9c5d7ba96f0bb0)},
      {"cubic_cusp", 0.0, 4u, UINT64_C(0xbade6eed950f5cd5)},
      {"rational_cubic_3d", 0.0, 16u, UINT64_C(0x4b22ed8b5249d9ba)},
      {"cubic_cusp_third", 0x1p-10, 6u, UINT64_C(0xefccf75eb5531336)},
      {"cubic_cusp_third", 0x1p-24, 34u, UINT64_C(0x6da33457daea2c33)},
      {"cubic_cusp_third", 0x1p-40, 66u, UINT64_C(0xa9cf6421460b8518)},
      {"cubic_loop", 0x1p-16, 4u, UINT64_C(0xa19c24fb84f818eb)},
      {"rational_cubic_heavy", 0x1p-14, 18u, UINT64_C(0x06acb62d9a57fd24)},
      {"rational_cubic_heavy", 0x1p-30, 36u, UINT64_C(0xea51aea8510d3e13)},
      {"quadratic", 0.0, 8u, UINT64_C(0x87d9aad4bb5332c5)},
      {"rational_quadratic", 0.0, 12u, UINT64_C(0xe6cdecc82077d44d)},
      {"quadratic_3d", 0x1p-10, 2u, UINT64_C(0xf4a189431b228976)},
      {"cubic_large", 0x1p-46, 20u, UINT64_C(0xfcedab83bb4a87d1)},
  };
  for (const Pin & pin : pins) {
    const AlongRef * r = ref_named(pin.name);
    ASSERT_NE(r, nullptr);
    GSHGEOM_Path * p = path_of(*r);
    GSHGEOM_Measured * m = measure(p, pin.tol);
    EXPECT_EQ(gshgeom_measured_row_count(m), pin.rows) << pin.name << " " << pin.tol;
    EXPECT_EQ(gshgeom_measured_hash(m), pin.hash) << pin.name << " " << pin.tol;
    gshgeom_measured_destroy(m);
    gshgeom_path_destroy(p);
  }
}

TEST(AlongDump, TheHashAndTheLengthAreAlwaysSixteenDigits) {
  bool leading_zero = false;
  for (int i = 0; i < 600; i++) {
    GSHGEOM_Path * p = new_path();
    ASSERT_EQ(gshgeom_path_move_to(p, P(i, 0)), GSHGEOM_OK);
    if (i % 2) {
      ASSERT_EQ(gshgeom_path_line_to(p, P(i, 3)), GSHGEOM_OK);
    }
    GSHGEOM_Measured * m = measure(p);
    leading_zero = leading_zero || gshgeom_measured_hash(m) < (UINT64_C(1) << 60);
    std::string d = dump_of(m);
    size_t at = d.find("hash=");
    ASSERT_NE(at, std::string::npos);
    EXPECT_EQ(d.find('\n', at) - (at + 5), 16u) << d;
    size_t len = d.find("length=");
    EXPECT_EQ(d.find(' ', len) - (len + 7), 16u) << d;
    gshgeom_measured_destroy(m);
    gshgeom_path_destroy(p);
  }
  EXPECT_TRUE(leading_zero) << "the search found no hash with a leading zero";
  GSHGEOM_Path * p = new_path();
  ASSERT_EQ(gshgeom_path_move_to(p, P(1, 1)), GSHGEOM_OK);
  GSHGEOM_Measured * m = measure(p);
  EXPECT_NE(dump_of(m).find("length=0000000000000000 "), std::string::npos);
  gshgeom_measured_destroy(m);
  gshgeom_path_destroy(p);
}

TEST(AlongDump, AFailingSinkIsReportedAndTheLocaleDoesNotMatter) {
  GSHGEOM_Path * p = subpath_fixture2();
  GSHGEOM_Measured * m = measure(p);
  FILE * full = std::fopen("/dev/full", "w");
  if (full) {
    EXPECT_EQ(gshgeom_measured_dump(m, full), GSHGEOM_ERR_IO);
    std::fclose(full);
  }
  FILE * ro = std::fopen("/dev/null", "r");
  ASSERT_NE(ro, nullptr);
  EXPECT_EQ(gshgeom_measured_dump(m, ro), GSHGEOM_ERR_IO);
  std::fclose(ro);
  EXPECT_EQ(gshgeom_measured_dump(m, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_measured_dump(nullptr, stdout), GSHGEOM_ERR_INVALID);
  const char * old = std::setlocale(LC_ALL, nullptr);
  std::string saved = old ? old : "C";
  for (const char * name : {"de_DE.UTF-8", "fr_FR.UTF-8", "de_DE", "C.UTF-8"}) {
    if (std::setlocale(LC_ALL, name)) {
      EXPECT_EQ(dump_of(m), kMeasuredDump2) << name;
    }
  }
  std::setlocale(LC_ALL, saved.c_str());
  gshgeom_measured_destroy(m);
  gshgeom_path_destroy(p);
}

TEST(AlongDump, EveryWriteIsChecked) {
  GSHGEOM_Path * p = subpath_fixture2();
  GSHGEOM_Measured * m = measure(p);
  GSHGEOM_EVERY_WRITE_CHECKED(gshgeom_measured_dump, m, "measured path");
  gshgeom_measured_destroy(m);
  gshgeom_path_destroy(p);
}

} // namespace

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
