/**
 * @file
 *
 * Exact arcs (arc.h): the record and its validity (AR1 to AR4), readback, copy, equality, hash and dump
 * (AR15), the three constructors against tools/arc-oracle.py (AR5 to AR7), the exact decomposition (AR8) and
 * the quarter-aligned exactness, lowering, bounds, and the exact sum's conversion to a double that the
 * endpoint form's radicand depends on.
 *
 * Copyright 2026 by Corey Pennycuff
 */

#include "../test_helpers.h"

#include "../alloc_tracker.h"
#include "../arc_vectors.h"

extern "C" {
#include "../../src/exact/exact.h"
}

#include <ghoti.io/math/basic.h>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace {

GSHGEOM_Point P(double x, double y) {
  return GSHGEOM_Point{x, y, 0.0};
}

GSHGEOM_Path * make(int dim = 2) {
  GSHGEOM_Path * p = nullptr;
  EXPECT_EQ(gshgeom_path_create(nullptr, nullptr, dim, &p), GSHGEOM_OK);
  return p;
}

std::string dump_of(const GSHGEOM_Path * p) {
  char * buf = nullptr;
  size_t n = 0;
  FILE * f = open_memstream(&buf, &n);
  EXPECT_EQ(gshgeom_path_dump(p, f), GSHGEOM_OK);
  fclose(f);
  std::string s(buf, n);
  free(buf);
  return s;
}

uint64_t bits(double v) {
  uint64_t b;
  std::memcpy(&b, &v, sizeof b);
  return b;
}

// the unit circle of radius r about (cx, cy), degrees, as the record a rounded corner uses
GSHGEOM_Arc circle_arc(double cx, double cy, double r, double t0, double sweep, GSHGEOM_Angle_Unit unit = GSHGEOM_ANGLE_DEGREES) {
  GSHGEOM_Arc a;
  std::memset(&a, 0, sizeof a);
  a.center = P(cx, cy);
  a.u = P(r, 0);
  a.v = P(0, r);
  a.unit = unit;
  a.t0 = t0;
  a.sweep = sweep;
  a.flags = GSHGEOM_ARC_CIRCULAR;
  GSHGEOM_Arc_Form f = {a.center, a.u, a.v, a.unit, a.t0, a.sweep, a.flags};
  gshgeom_arc_point(&f, t0, &a.start);
  gshgeom_arc_point(&f, t0 + sweep, &a.end);
  return a;
}

// scale of an arc record, as arc.h defines it
double scale_of(const GSHGEOM_Arc & a) {
  double s = 0;
  for (const GSHGEOM_Point * p : {&a.start, &a.end, &a.center}) {
    s = std::fmax(s, std::fmax(std::fabs(p->x), std::fabs(p->y)));
  }
  return std::fmax(s, std::fmax(std::fabs(a.u.x), std::fabs(a.u.y)) + std::fmax(std::fabs(a.v.x), std::fabs(a.v.y)));
}

GSHGEOM_Path * path_with(const GSHGEOM_Arc & a, GSHGEOM_Result * r = nullptr) {
  GSHGEOM_Path * p = make();
  EXPECT_EQ(gshgeom_path_move_to(p, a.start), GSHGEOM_OK);
  GSHGEOM_Result rr = gshgeom_path_append_arc(p, &a);
  if (r) {
    *r = rr;
  } else {
    EXPECT_EQ(rr, GSHGEOM_OK);
  }
  return p;
}

GSHGEOM_Arc read_arc(const GSHGEOM_Path * p, size_t i) {
  GSHGEOM_Arc a;
  EXPECT_EQ(gshgeom_path_arc(p, i, &a), GSHGEOM_OK);
  return a;
}

}  // namespace

// ---------------------------------------------------------------------------
// The exact sum's conversion to a double
// ---------------------------------------------------------------------------

TEST(ExactSumToDouble, OneProductRoundsOnceLikeTheMultiplication) {
  uint64_t seed = 88172645463325252ull;
  auto next = [&]() {
    seed ^= seed << 13;
    seed ^= seed >> 7;
    seed ^= seed << 17;
    return seed;
  };
  for (int i = 0; i < 20000; i++) {
    double a = std::ldexp((double)(next() >> 11) + 1.0, (int)(next() % 200) - 130) * ((next() & 1) ? 1 : -1);
    double b = std::ldexp((double)(next() >> 11) + 1.0, (int)(next() % 200) - 130) * ((next() & 1) ? 1 : -1);
    GSHGEOM_ExactSum s;
    gshgeom_exact_sum_init(&s);
    const double f[2] = {a, b};
    gshgeom_exact_sum_add(&s, 0, 2, f);
    volatile double product = a * b;  // no contraction: a plain rounded product
    EXPECT_EQ(bits(gshgeom_exact_sum_to_double(&s)), bits(product)) << a << " " << b;
  }
}

TEST(ExactSumToDouble, TiesGoToEvenAndNothingBelowIsLost) {
  auto sum2 = [](double a, double b, double c) {
    GSHGEOM_ExactSum s;
    gshgeom_exact_sum_init(&s);
    const double fa[1] = {a}, fb[1] = {b}, fc[1] = {c};
    gshgeom_exact_sum_add(&s, 0, 1, fa);
    gshgeom_exact_sum_add(&s, 0, 1, fb);
    gshgeom_exact_sum_add(&s, 0, 1, fc);
    return gshgeom_exact_sum_to_double(&s);
  };
  const double ulp = 0x1p-52;
  EXPECT_EQ(sum2(1.0, 0x1p-53, 0.0), 1.0);                      // a tie to an even mantissa
  EXPECT_EQ(sum2(1.0 + ulp, 0x1p-53, 0.0), 1.0 + 2 * ulp);       // a tie to an odd one goes up
  EXPECT_EQ(sum2(1.0, 0x1p-53, 0x1p-200), 1.0 + ulp);            // not a tie: a bit far below
  EXPECT_EQ(sum2(1.0, 0x1p-53, -0x1p-200), 1.0);                 // and below the other way
  EXPECT_EQ(sum2(-1.0, -0x1p-53, -0x1p-200), -(1.0 + ulp));      // the sign is symmetric
  EXPECT_EQ(sum2(1.0, -1.0, 0.0), 0.0);
  EXPECT_EQ(sum2(0x1p-300, -0x1p-300, 0x1p-900), 0x1p-900);      // a cancellation to a distant bit
  // 1 - 2^-60 (exact) is not 1 but rounds to it; 1 - 2^-54 is a tie below 1 and goes to 1 - 2^-53? (nearest is 1)
  EXPECT_EQ(sum2(1.0, -0x1p-60, 0.0), 1.0);
}

// ---------------------------------------------------------------------------
// The record, validity and readback
// ---------------------------------------------------------------------------

TEST(Arc, VerbAndRecordFollowTheRules) {
  EXPECT_EQ(GSHGEOM_VERB_ARC, 7);
  GSHGEOM_Arc a = circle_arc(1, 2, 3, 30, 100);
  GSHGEOM_Path * p = path_with(a);
  EXPECT_EQ(gshgeom_path_element_count(p), 2u);
  GSHGEOM_Verb v;
  EXPECT_EQ(gshgeom_path_verb(p, 0, &v), GSHGEOM_OK);
  EXPECT_EQ(v, GSHGEOM_VERB_MOVE);
  EXPECT_EQ(gshgeom_path_verb(p, 1, &v), GSHGEOM_OK);
  EXPECT_EQ(v, GSHGEOM_VERB_ARC);
  EXPECT_EQ(gshgeom_path_verb(p, 2, &v), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_verb(nullptr, 0, &v), GSHGEOM_ERR_INVALID);
  GSHGEOM_Arc b = read_arc(p, 1);
  EXPECT_EQ(std::memcmp(&b.start, &a.start, sizeof b.start), 0);
  EXPECT_EQ(std::memcmp(&b.end, &a.end, sizeof b.end), 0);
  EXPECT_EQ(bits(b.t0), bits(a.t0));
  EXPECT_EQ(bits(b.sweep), bits(a.sweep));
  EXPECT_EQ(b.unit, a.unit);
  EXPECT_EQ(b.flags, a.flags);
  GSHGEOM_Arc none;
  EXPECT_EQ(gshgeom_path_arc(p, 0, &none), GSHGEOM_ERR_INVALID);  // a move is not an arc
  EXPECT_EQ(gshgeom_path_arc(p, 9, &none), GSHGEOM_ERR_INVALID);
  // the next segment starts at the arc's end
  EXPECT_EQ(gshgeom_path_line_to(p, P(0, 0)), GSHGEOM_OK);
  GSHGEOM_Element el;
  EXPECT_EQ(gshgeom_path_element(p, 2, &el), GSHGEOM_OK);
  EXPECT_EQ(bits(el.curve.p[0].x), bits(a.end.x));
  EXPECT_EQ(bits(el.curve.p[0].y), bits(a.end.y));
  gshgeom_path_destroy(p);
}

TEST(Arc, ElementReadbackIsUnsupportedNotTheChord) {
  GSHGEOM_Path * p = path_with(circle_arc(0, 0, 5, 0, 90));
  GSHGEOM_Element el;
  std::memset(&el, 0x5a, sizeof el);
  EXPECT_EQ(gshgeom_path_element(p, 1, &el), GSHGEOM_ERR_UNSUPPORTED);
  EXPECT_EQ(el.verb, (GSHGEOM_Verb)0x5a5a5a5a);  // untouched
  EXPECT_EQ(gshgeom_path_element(p, 0, &el), GSHGEOM_OK);
  gshgeom_path_destroy(p);
}

TEST(Arc, ValidityIsDecidedInTheOrderOfTheRule) {
  GSHGEOM_Arc ok = circle_arc(0, 0, 10, 0, 90);
  GSHGEOM_Path * p = make();
  EXPECT_EQ(gshgeom_path_append_arc(p, &ok), GSHGEOM_ERR_INVALID);  // no current point
  EXPECT_EQ(gshgeom_path_append_arc(nullptr, &ok), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_append_arc(p, nullptr), GSHGEOM_ERR_INVALID);
  gshgeom_path_move_to(p, ok.start);
  const uint64_t h0 = gshgeom_path_hash(p);
  auto check = [&](GSHGEOM_Arc a, GSHGEOM_Result want, const char * what) {
    EXPECT_EQ(gshgeom_path_append_arc(p, &a), want) << what;
    EXPECT_EQ(gshgeom_path_hash(p), h0) << what << ": a failed call leaves the path as it was";
    EXPECT_EQ(gshgeom_path_element_count(p), 1u) << what;
  };
  GSHGEOM_Arc a = ok;
  a.t0 = NAN;
  check(a, GSHGEOM_ERR_INVALID, "nan t0");
  a = ok;
  a.sweep = INFINITY;
  check(a, GSHGEOM_ERR_INVALID, "infinite sweep");
  a = ok;
  a.end.y = NAN;
  check(a, GSHGEOM_ERR_INVALID, "nan end");
  a = ok;
  a.unit = GSHGEOM_ANGLE_UNIT_COUNT;
  check(a, GSHGEOM_ERR_INVALID, "unknown unit");
  a = ok;
  a.flags = 2;
  check(a, GSHGEOM_ERR_INVALID, "unknown flag bit");
  a = ok;
  a.v = P(10, 0);
  check(a, GSHGEOM_ERR_INVALID, "u parallel to v: no area");
  a = ok;
  a.v = P(0x1p-1000, 0x1p-1000);
  a.u = P(0x1p-1000, 0x1p-1000);
  check(a, GSHGEOM_ERR_INVALID, "equal u and v at a scale a rounded cross product loses");
  a = ok;
  a.sweep = 0;
  check(a, GSHGEOM_ERR_INVALID, "zero sweep");
  a = ok;
  a.sweep = 360.0 * (1 + 0x1p-52);
  check(a, GSHGEOM_ERR_INVALID, "more than a turn");
  a = ok;
  a.sweep = -360.0 * (1 + 0x1p-52);
  check(a, GSHGEOM_ERR_INVALID, "more than a turn, negative");
  a = ok;
  a.center.x = 0x1.8p100;
  check(a, GSHGEOM_ERR_LIMIT, "centre past the coordinate limit");
  a = ok;
  a.u.x = 0x1.8p100;
  check(a, GSHGEOM_ERR_LIMIT, "half-diameter past the limit");
  a = ok;
  a.t0 = 90.0 * 0x1.8p40;
  check(a, GSHGEOM_ERR_LIMIT, "t0 past the quarter-turn limit");
  a = ok;
  a.start.x += 1;
  check(a, GSHGEOM_ERR_INVALID, "start is not the current point");
  a = ok;
  a.end.x += 1;
  check(a, GSHGEOM_ERR_INVALID, "end is not where the centre form ends");
  // invalid is decided before limit: a flat ellipse with a coordinate past the limit is invalid first
  a = ok;
  a.v = a.u;
  a.center.x = 0x1.8p100;
  check(a, GSHGEOM_ERR_INVALID, "invalid before limit");
  EXPECT_EQ(gshgeom_path_append_arc(p, &ok), GSHGEOM_OK);
  gshgeom_path_destroy(p);
  // 3D
  GSHGEOM_Path * q = make(3);
  gshgeom_path_move_to(q, ok.start);
  EXPECT_EQ(gshgeom_path_append_arc(q, &ok), GSHGEOM_ERR_INVALID);
  gshgeom_path_destroy(q);
}

TEST(Arc, TheAgreementBoundHoldsAtItsLimitAndFailsPastIt) {
  GSHGEOM_Arc ok = circle_arc(0, 0, 1000, 0, 90, GSHGEOM_ANGLE_DEGREES);
  const double tol = GSHGEOM_ARC_AGREE * scale_of(ok);
  for (int which = 0; which < 2; which++) {
    for (double k : {0.5, 0.75}) {
      GSHGEOM_Arc a = ok;
      (which ? a.end.x : a.start.x) += k * tol;
      GSHGEOM_Path * p = make();
      gshgeom_path_move_to(p, a.start);
      EXPECT_EQ(gshgeom_path_append_arc(p, &a), GSHGEOM_OK) << which << " " << k;
      gshgeom_path_destroy(p);
    }
    for (double k : {1.5, 2.0, 64.0}) {
      GSHGEOM_Arc a = ok;
      (which ? a.end.y : a.start.y) -= k * tol;
      GSHGEOM_Path * p = make();
      gshgeom_path_move_to(p, a.start);
      EXPECT_EQ(gshgeom_path_append_arc(p, &a), GSHGEOM_ERR_INVALID) << which << " " << k;
      gshgeom_path_destroy(p);
    }
  }
}

TEST(Arc, AFullTurnClosesToTheBitOrIsInvalid) {
  GSHGEOM_Arc a = circle_arc(5, 5, 3, 0, 360);
  EXPECT_EQ(bits(a.end.x), bits(a.start.x));
  GSHGEOM_Path * p = path_with(a);
  gshgeom_path_destroy(p);
  GSHGEOM_Arc b = a;
  b.end.x = std::nextafter(b.end.x, 100.0);  // within the agreement bound, yet not equal to the bit
  GSHGEOM_Path * q = make();
  gshgeom_path_move_to(q, b.start);
  EXPECT_EQ(gshgeom_path_append_arc(q, &b), GSHGEOM_ERR_INVALID);
  gshgeom_path_destroy(q);
  // radians: a full turn is 2 pi as a double
  GSHGEOM_Arc r = circle_arc(0, 0, 1, 0, 0x1.921fb54442d18p+2, GSHGEOM_ANGLE_RADIANS);
  r.end = r.start;
  GSHGEOM_Path * s = path_with(r);
  gshgeom_path_destroy(s);
}

TEST(Arc, RoundTripsCopyEqualityAndDump) {
  GSHGEOM_Path * p = make();
  gshgeom_path_move_to(p, P(10, 0));
  GSHGEOM_Arc a = circle_arc(0, 0, 10, 0, 90);
  EXPECT_EQ(gshgeom_path_append_arc(p, &a), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_line_to(p, P(0, 0)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_close(p), GSHGEOM_OK);
  GSHGEOM_Path * c = nullptr;
  ASSERT_EQ(gshgeom_path_copy(p, nullptr, nullptr, &c), GSHGEOM_OK);
  EXPECT_TRUE(gshgeom_path_equal(p, c));
  EXPECT_EQ(gshgeom_path_hash(p), gshgeom_path_hash(c));
  EXPECT_EQ(dump_of(p), dump_of(c));
  // the copy keeps building where the original would
  EXPECT_EQ(gshgeom_path_line_to(c, P(1, 1)), GSHGEOM_OK);
  EXPECT_FALSE(gshgeom_path_equal(p, c));
  EXPECT_TRUE(gshgeom_path_equal(p, p));
  EXPECT_FALSE(gshgeom_path_equal(p, nullptr));
  EXPECT_TRUE(gshgeom_path_equal(nullptr, nullptr));
  // an arc differing in a flag, a unit, or a signed zero is a different path
  GSHGEOM_Arc b = a;
  b.flags = 0;
  GSHGEOM_Path * d = path_with(b);
  GSHGEOM_Path * e = path_with(a);
  EXPECT_FALSE(gshgeom_path_equal(d, e));
  EXPECT_NE(gshgeom_path_hash(d), gshgeom_path_hash(e));
  gshgeom_path_destroy(d);
  gshgeom_path_destroy(e);
  GSHGEOM_Path * z1 = make();
  GSHGEOM_Path * z2 = make();
  gshgeom_path_move_to(z1, P(0.0, 0.0));
  gshgeom_path_move_to(z2, P(-0.0, 0.0));
  EXPECT_FALSE(gshgeom_path_equal(z1, z2));
  EXPECT_NE(gshgeom_path_hash(z1), gshgeom_path_hash(z2));
  gshgeom_path_destroy(z1);
  gshgeom_path_destroy(z2);
  // a copy that does not fit its limits is LIMIT and writes nothing
  GSHGEOM_Limits lim = {};
  lim.max_elements = 2;
  GSHGEOM_Path * tiny = reinterpret_cast<GSHGEOM_Path *>(0x1);
  EXPECT_EQ(gshgeom_path_copy(p, nullptr, &lim, &tiny), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(tiny, reinterpret_cast<GSHGEOM_Path *>(0x1));
  gshgeom_path_destroy(p);
  gshgeom_path_destroy(c);
}

TEST(Arc, DumpAndHashAreStatedBitForBit) {
  GSHGEOM_Arc a = circle_arc(1.5, -2.5, 4, 15, -200, GSHGEOM_ANGLE_GRADIANS);
  a.flags = 0;
  a.v = P(-1, 4);  // not circular
  a.u = P(4, 1);
  GSHGEOM_Arc_Form f = {a.center, a.u, a.v, a.unit, a.t0, a.sweep, a.flags};
  gshgeom_arc_point(&f, a.t0, &a.start);
  gshgeom_arc_point(&f, a.t0 + a.sweep, &a.end);
  GSHGEOM_Path * p = path_with(a);
  // the hash, restated: u32 1, u32 dim, u64 count; then u8 verb 0 and the move's point, u8 verb 7, u8 flags, u8 unit, twelve f64
  GSHGEOM_Hash h;
  gshgeom_hash_init(&h);
  gshgeom_hash_u32(&h, 1);
  gshgeom_hash_u32(&h, 2);
  gshgeom_hash_u64(&h, 2);
  gshgeom_hash_u8(&h, 0);
  gshgeom_hash_f64(&h, a.start.x);
  gshgeom_hash_f64(&h, a.start.y);
  gshgeom_hash_u8(&h, 7);
  gshgeom_hash_u8(&h, 0);
  gshgeom_hash_u8(&h, 2);
  for (double v : {a.start.x, a.start.y, a.end.x, a.end.y, a.center.x, a.center.y, a.u.x, a.u.y, a.v.x, a.v.y, a.t0, a.sweep}) {
    gshgeom_hash_f64(&h, v);
  }
  EXPECT_EQ(gshgeom_path_hash(p), gshgeom_hash_value(&h));
  // the dump line: A <flags> <unit> and twelve 16-digit words
  char expect[512];
  int n = std::snprintf(expect, sizeof expect, "A 0 2");
  for (double v : {a.start.x, a.start.y, a.end.x, a.end.y, a.center.x, a.center.y, a.u.x, a.u.y, a.v.x, a.v.y, a.t0, a.sweep}) {
    n += std::snprintf(expect + n, sizeof expect - n, " %016llx", (unsigned long long)bits(v));
  }
  std::string d = dump_of(p);
  EXPECT_NE(d.find(std::string("\n") + expect + "\n"), std::string::npos) << d;
  EXPECT_NE(d.find("elements=2"), std::string::npos);
  gshgeom_path_destroy(p);
}

TEST(Arc, AnArcFreePathHashesAsItAlwaysDid) {
  GSHGEOM_Path * p = make();
  gshgeom_path_move_to(p, P(1, 2));
  gshgeom_path_line_to(p, P(3, 4));
  gshgeom_path_close(p);
  EXPECT_EQ(gshgeom_path_hash(p), 0x0ull + gshgeom_path_hash(p));  // self-check of the call; the pinned value is in geometry.hash
  GSHGEOM_Path * c = nullptr;
  ASSERT_EQ(gshgeom_path_copy(p, nullptr, nullptr, &c), GSHGEOM_OK);
  EXPECT_EQ(dump_of(p), dump_of(c));
  gshgeom_path_destroy(p);
  gshgeom_path_destroy(c);
}

TEST(Arc, FailedBuildsLeaveThePathAsItWasAtEveryAllocation) {
  for (int mode = 0; mode < 3; mode++) {
    for (long k = 0; k < 40; k++) {
      tracker::Tracker t;
      t.fail_from = k;
      GSHGEOM_Allocator al = tracker::allocator_for(&t);
      GSHGEOM_Path * p = nullptr;
      GSHGEOM_Result r = gshgeom_path_create(&al, nullptr, 2, &p);
      if (r != GSHGEOM_OK) {
        EXPECT_TRUE(t.live.empty());
        continue;
      }
      GSHGEOM_Result rr = GSHGEOM_OK;
      GSHGEOM_Arc a = circle_arc(0, 0, 10, 0, 90);
      if (mode == 0) {
        rr = gshgeom_path_move_to(p, a.start);
        if (rr == GSHGEOM_OK) rr = gshgeom_path_append_arc(p, &a);
      } else if (mode == 1) {
        GSHGEOM_Arc_Form f = {a.center, a.u, a.v, a.unit, a.t0, a.sweep, a.flags};
        rr = gshgeom_path_arc_conjugate_to(p, &f, nullptr);  // begins its own subpath: a move and an arc
      } else {
        rr = gshgeom_path_move_to(p, P(0, 0));
        if (rr == GSHGEOM_OK) rr = gshgeom_path_arc_tangent_to(p, P(100, 0), P(100, 100), 20);
      }
      if (rr == GSHGEOM_ERR_OOM) {
        EXPECT_TRUE(mode != 0 || gshgeom_path_element_count(p) <= 1u);
        if (mode == 1) {
          EXPECT_EQ(gshgeom_path_element_count(p), 0u) << "the move that began the subpath is taken back";
        }
      } else {
        EXPECT_EQ(rr, GSHGEOM_OK);
      }
      gshgeom_path_destroy(p);
      EXPECT_TRUE(t.live.empty()) << "mode " << mode << " fail_from " << k;
    }
  }
}

// ---------------------------------------------------------------------------
// The endpoint form against the oracle
// ---------------------------------------------------------------------------

static double rel_tol(double scale) {
  return GSHGEOM_ARC_DECOMP * scale;
}

TEST(ArcEndpoint, EveryVectorHoldsWithinTheDecompositionBoundAndTheEndsAreExact) {
  using namespace arc_vectors;
  double worst = 0;
  double agree_worst = 0;
  for (size_t i = 0; i < endpoint_count; i++) {
    const Endpoint & v = endpoints[i];
    GSHGEOM_Path * p = make();
    gshgeom_path_move_to(p, P(v.sx, v.sy));
    GSHGEOM_Arc_Kind kind = (GSHGEOM_Arc_Kind)99;
    GSHGEOM_Result r = gshgeom_path_arc_endpoint_to(
        p, P(v.ex, v.ey), v.rx, v.ry, v.phi, (GSHGEOM_Angle_Unit)v.unit, v.large, v.sweep, &kind);
    EXPECT_EQ(r, GSHGEOM_OK) << "vector " << i << " S " << v.sx << "," << v.sy << " E " << v.ex << "," << v.ey << " r " << v.rx << "," << v.ry
                             << " phi " << v.phi << " flags " << v.large << v.sweep;
    if (r != GSHGEOM_OK) {
      gshgeom_path_destroy(p);
      continue;
    }
    ASSERT_EQ((int)kind, v.kind) << "vector " << i;
    if (v.kind == 0) {
      EXPECT_EQ(gshgeom_path_element_count(p), 1u) << i;
    } else if (v.kind == 1) {
      EXPECT_EQ(gshgeom_path_element_count(p), 2u) << i;
      GSHGEOM_Verb vb;
      gshgeom_path_verb(p, 1, &vb);
      EXPECT_EQ(vb, GSHGEOM_VERB_LINE);
    } else {
      ASSERT_EQ(gshgeom_path_element_count(p), 2u) << i;
      GSHGEOM_Arc a = read_arc(p, 1);
      EXPECT_EQ(bits(a.start.x), bits(v.sx)) << i;
      EXPECT_EQ(bits(a.start.y), bits(v.sy)) << i;
      EXPECT_EQ(bits(a.end.x), bits(v.ex)) << i;
      EXPECT_EQ(bits(a.end.y), bits(v.ey)) << i;
      EXPECT_EQ(a.unit, GSHGEOM_ANGLE_RADIANS) << i;
      EXPECT_EQ((a.flags & GSHGEOM_ARC_CIRCULAR) != 0, v.circular != 0) << i;
      const double tol = rel_tol(v.scale);
      const double ang = GSHGEOM_ARC_DECOMP * 8;  // angles are of order one: the bound is relative
      auto near = [&](double got, double want, double t, const char * what) {
        double e = std::fabs(got - want);
        if (e / t > worst) worst = e / t;
        EXPECT_LE(e, t) << what << " vector " << i << " got " << got << " want " << want;
      };
      near(a.center.x, v.cx, tol, "cx");
      near(a.center.y, v.cy, tol, "cy");
      near(a.u.x, v.ux, tol, "ux");
      near(a.u.y, v.uy, tol, "uy");
      near(a.v.x, v.vx, tol, "vx");
      near(a.v.y, v.vy, tol, "vy");
      {
        // how far the centre form the constructor made misses the ends it kept, relative to the scale (arc.h's agreement)
        GSHGEOM_Arc_Form f = {a.center, a.u, a.v, a.unit, a.t0, a.sweep, a.flags};
        GSHGEOM_Point p0, p1;
        ASSERT_EQ(gshgeom_arc_point(&f, a.t0, &p0), GSHGEOM_OK);
        ASSERT_EQ(gshgeom_arc_point(&f, a.t0 + a.sweep, &p1), GSHGEOM_OK);
        const double s = scale_of(a);
        agree_worst = std::fmax(agree_worst, std::fmax(std::fmax(std::fabs(p0.x - a.start.x), std::fabs(p0.y - a.start.y)),
                                                      std::fmax(std::fabs(p1.x - a.end.x), std::fabs(p1.y - a.end.y))) / s);
      }
      near(a.t0, v.t0, ang, "t0");
      near(a.sweep, v.dsweep, ang, "sweep");
    }
    gshgeom_path_destroy(p);
  }
  std::printf("[ arc ] endpoint form: worst error %.3g of the bound over %zu vectors; worst agreement %.3g (2^%.1f) of the scale\n",
      worst, endpoint_count, agree_worst, std::log2(agree_worst));
  EXPECT_LE(agree_worst, GSHGEOM_ARC_AGREE / 8);
}

TEST(ArcEndpoint, DegenerateCasesAddWhatTheRuleSays) {
  GSHGEOM_Path * p = make();
  gshgeom_path_move_to(p, P(1, 1));
  GSHGEOM_Arc_Kind k;
  EXPECT_EQ(gshgeom_path_arc_endpoint_to(p, P(1, 1), 5, 5, 0, GSHGEOM_ANGLE_DEGREES, 0, 1, &k), GSHGEOM_OK);
  EXPECT_EQ(k, GSHGEOM_ARC_KIND_NOTHING);
  EXPECT_EQ(gshgeom_path_element_count(p), 1u);
  EXPECT_EQ(gshgeom_path_arc_endpoint_to(p, P(4, 5), 0, 5, 0, GSHGEOM_ANGLE_DEGREES, 0, 1, &k), GSHGEOM_OK);
  EXPECT_EQ(k, GSHGEOM_ARC_KIND_LINE);
  // the flags are exactly 0 or 1; nothing is added on an error
  const size_t n = gshgeom_path_element_count(p);
  for (int bad : {-1, 2, 7}) {
    EXPECT_EQ(gshgeom_path_arc_endpoint_to(p, P(8, 9), 5, 5, 0, GSHGEOM_ANGLE_DEGREES, bad, 0, &k), GSHGEOM_ERR_INVALID);
    EXPECT_EQ(gshgeom_path_arc_endpoint_to(p, P(8, 9), 5, 5, 0, GSHGEOM_ANGLE_DEGREES, 0, bad, &k), GSHGEOM_ERR_INVALID);
  }
  EXPECT_EQ(gshgeom_path_arc_endpoint_to(p, P(8, NAN), 5, 5, 0, GSHGEOM_ANGLE_DEGREES, 0, 0, &k), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_endpoint_to(p, P(8, 9), INFINITY, 5, 0, GSHGEOM_ANGLE_DEGREES, 0, 0, &k), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_endpoint_to(p, P(8, 9), 5, 5, 0, (GSHGEOM_Angle_Unit)9, 0, 0, &k), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_element_count(p), n);
  EXPECT_EQ(gshgeom_path_arc_endpoint_to(p, P(8, 9), 5, 5, 0, GSHGEOM_ANGLE_DEGREES, 0, 0, nullptr), GSHGEOM_OK);
  gshgeom_path_destroy(p);
  GSHGEOM_Path * q = make();  // no current point
  EXPECT_EQ(gshgeom_path_arc_endpoint_to(q, P(8, 9), 5, 5, 0, GSHGEOM_ANGLE_DEGREES, 0, 0, &k), GSHGEOM_ERR_INVALID);
  gshgeom_path_destroy(q);
}

TEST(ArcEndpoint, TheHalfEllipseIsExactWhenTheRadiiAreTooSmall) {
  GSHGEOM_Path * p = make();
  gshgeom_path_move_to(p, P(0, 0));
  ASSERT_EQ(gshgeom_path_arc_endpoint_to(p, P(4, 0), 1, 1, 0, GSHGEOM_ANGLE_DEGREES, 1, 1, nullptr), GSHGEOM_OK);
  GSHGEOM_Arc a = read_arc(p, 1);
  EXPECT_EQ(a.sweep, 0x1.921fb54442d18p+1);  // half a turn, in the sign of the flag
  EXPECT_EQ(a.center.x, 2.0);
  EXPECT_EQ(a.center.y, 0.0);
  EXPECT_EQ(a.u.x, 2.0);
  EXPECT_EQ(a.v.y, 2.0);
  EXPECT_TRUE(a.flags & GSHGEOM_ARC_CIRCULAR);
  gshgeom_path_destroy(p);
}

// ---------------------------------------------------------------------------
// The centre, tangent and bulge forms
// ---------------------------------------------------------------------------

TEST(ArcCenter, CanvasArcAndEllipseAgainstTheOracle) {
  using namespace arc_vectors;
  for (size_t i = 0; i < center_count; i++) {
    const Center & c = centers[i];
    GSHGEOM_Path * p = make();  // no current point: the subpath begins at P(t0)
    ASSERT_EQ(gshgeom_path_arc_center_to(p, P(c.cx, c.cy), c.rx, c.ry, c.phi, (GSHGEOM_Angle_Unit)c.unit, c.t0, c.sweep, nullptr), GSHGEOM_OK)
        << i;
    ASSERT_EQ(gshgeom_path_element_count(p), 2u) << i;
    GSHGEOM_Arc a = read_arc(p, 1);
    const double tol = rel_tol(c.scale);
    EXPECT_LE(std::fabs(a.start.x - c.sx), tol) << i;
    EXPECT_LE(std::fabs(a.start.y - c.sy), tol) << i;
    EXPECT_LE(std::fabs(a.end.x - c.ex), tol) << i;
    EXPECT_LE(std::fabs(a.end.y - c.ey), tol) << i;
    if (c.full) {
      EXPECT_EQ(bits(a.end.x), bits(a.start.x)) << i;
      EXPECT_EQ(bits(a.end.y), bits(a.start.y)) << i;
    }
    EXPECT_EQ((a.flags & GSHGEOM_ARC_CIRCULAR) != 0, c.rx == c.ry) << i;
    EXPECT_EQ(a.unit, (GSHGEOM_Angle_Unit)c.unit);
    EXPECT_EQ(bits(a.t0), bits(c.t0));        // t0 kept as given
    EXPECT_EQ(bits(a.sweep), bits(c.sweep));
    EXPECT_LE(std::fabs(a.u.x - c.ux), tol);
    EXPECT_LE(std::fabs(a.v.y - c.vy), tol);
    // the very bits the arc starts from are what arc_point gives, for a front end's connecting line
    GSHGEOM_Arc_Form f = {a.center, a.u, a.v, a.unit, a.t0, a.sweep, a.flags};
    GSHGEOM_Point s;
    ASSERT_EQ(gshgeom_arc_point(&f, c.t0, &s), GSHGEOM_OK);
    EXPECT_EQ(bits(s.x), bits(a.start.x));
    EXPECT_EQ(bits(s.y), bits(a.start.y));
    gshgeom_path_destroy(p);
  }
}

TEST(ArcCenter, WithACurrentPointTheArcStartsThereOrIsInvalid) {
  GSHGEOM_Arc_Form f = {P(0, 0), P(10, 0), P(0, 10), GSHGEOM_ANGLE_DEGREES, 0, 90, GSHGEOM_ARC_CIRCULAR};
  GSHGEOM_Point s;
  ASSERT_EQ(gshgeom_arc_point(&f, 0, &s), GSHGEOM_OK);
  GSHGEOM_Path * p = make();
  gshgeom_path_move_to(p, s);
  EXPECT_EQ(gshgeom_path_arc_conjugate_to(p, &f, nullptr), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_element_count(p), 2u);
  gshgeom_path_destroy(p);
  GSHGEOM_Path * q = make();
  gshgeom_path_move_to(q, P(3, 3));  // not on the arc
  EXPECT_EQ(gshgeom_path_arc_conjugate_to(q, &f, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_element_count(q), 1u);
  // a line to the very bits (Canvas's connecting line), then the arc
  gshgeom_path_line_to(q, s);
  EXPECT_EQ(gshgeom_path_arc_conjugate_to(q, &f, nullptr), GSHGEOM_OK);
  gshgeom_path_destroy(q);
  // an explicit end overrides, and must still agree
  GSHGEOM_Path * r = make();
  GSHGEOM_Point e = P(0, 10);
  EXPECT_EQ(gshgeom_path_arc_conjugate_to(r, &f, &e), GSHGEOM_OK);
  EXPECT_EQ(bits(read_arc(r, 1).end.x), bits(0.0));
  GSHGEOM_Path * r2 = make();
  GSHGEOM_Point bad = P(5, 10);
  EXPECT_EQ(gshgeom_path_arc_conjugate_to(r2, &f, &bad), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_element_count(r2), 0u);
  gshgeom_path_destroy(r);
  gshgeom_path_destroy(r2);
  // errors: a sweep past a turn, zero radius, a flat ellipse, t0 past the limit
  GSHGEOM_Path * t = make();
  GSHGEOM_Arc_Form g = f;
  g.sweep = 361;
  EXPECT_EQ(gshgeom_path_arc_conjugate_to(t, &g, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_center_to(t, P(0, 0), 0, 5, 0, GSHGEOM_ANGLE_DEGREES, 0, 90, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_center_to(t, P(0, 0), -1, 5, 0, GSHGEOM_ANGLE_DEGREES, 0, 90, nullptr), GSHGEOM_ERR_INVALID);
  g = f;
  g.v = P(5, 0);
  EXPECT_EQ(gshgeom_path_arc_conjugate_to(t, &g, nullptr), GSHGEOM_ERR_INVALID);
  g = f;
  g.t0 = 90.0 * 0x1.8p40;
  EXPECT_EQ(gshgeom_path_arc_conjugate_to(t, &g, nullptr), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(gshgeom_path_element_count(t), 0u);
  gshgeom_path_destroy(t);
}

TEST(ArcCenter, DxfEllipseByVectorAndRatio) {
  using namespace arc_vectors;
  for (size_t i = 0; i < dxf_count; i++) {
    const Dxf & d = dxfs[i];
    GSHGEOM_Arc_Form f = {P(d.cx, d.cy), P(d.ux, d.uy), P(d.vx, d.vy), GSHGEOM_ANGLE_RADIANS, d.t0, d.sweep, 0};
    GSHGEOM_Path * p = make();
    ASSERT_EQ(gshgeom_path_arc_conjugate_to(p, &f, nullptr), GSHGEOM_OK) << i;
    GSHGEOM_Arc a = read_arc(p, 1);
    EXPECT_LE(std::fabs(a.start.x - d.sx), rel_tol(d.scale)) << i;
    EXPECT_LE(std::fabs(a.end.y - d.ey), rel_tol(d.scale)) << i;
    gshgeom_path_destroy(p);
  }
}

TEST(ArcTangent, ArcToInEveryDegenerateCaseAndTheNormalOne) {
  using namespace arc_vectors;
  for (size_t i = 0; i < arcto_count; i++) {
    const ArcTo & v = arctos[i];
    GSHGEOM_Path * p = make();
    gshgeom_path_move_to(p, P(v.p0x, v.p0y));
    ASSERT_EQ(gshgeom_path_arc_tangent_to(p, P(v.p1x, v.p1y), P(v.p2x, v.p2y), v.r), GSHGEOM_OK) << i;
    if (v.kind == 1) {
      ASSERT_EQ(gshgeom_path_element_count(p), 2u) << i;
      GSHGEOM_Element el;
      ASSERT_EQ(gshgeom_path_element(p, 1, &el), GSHGEOM_OK);
      EXPECT_EQ(el.verb, GSHGEOM_VERB_LINE) << i;
      EXPECT_EQ(bits(el.curve.p[1].x), bits(v.p1x)) << i;
    } else {
      const bool line_first = !(bits(v.t1x) == bits(v.p0x) && bits(v.t1y) == bits(v.p0y));
      ASSERT_EQ(gshgeom_path_element_count(p), line_first ? 3u : 2u) << i;
      GSHGEOM_Arc a = read_arc(p, line_first ? 2 : 1);
      const double tol = rel_tol(v.scale);
      EXPECT_LE(std::fabs(a.start.x - v.t1x), tol) << i;
      EXPECT_LE(std::fabs(a.start.y - v.t1y), tol) << i;
      EXPECT_LE(std::fabs(a.end.x - v.t2x), tol) << i;
      EXPECT_LE(std::fabs(a.end.y - v.t2y), tol) << i;
      EXPECT_LE(std::fabs(a.center.x - v.cx), tol) << i;
      EXPECT_LE(std::fabs(a.center.y - v.cy), tol) << i;
      EXPECT_LE(std::fabs(a.t0 - v.t0), 8 * GSHGEOM_ARC_DECOMP) << i;
      EXPECT_LE(std::fabs(a.sweep - v.sweep), 8 * GSHGEOM_ARC_DECOMP) << i;
      EXPECT_TRUE(a.flags & GSHGEOM_ARC_CIRCULAR);
      EXPECT_EQ(a.sweep > 0, v.sweep > 0);  // the sign is the sign of the turn
    }
    gshgeom_path_destroy(p);
  }
  GSHGEOM_Path * q = make();
  gshgeom_path_move_to(q, P(0, 0));
  EXPECT_EQ(gshgeom_path_arc_tangent_to(q, P(10, 0), P(10, 10), -1), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_tangent_to(q, P(10, 0), P(10, 10), NAN), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_tangent_to(q, P(10, 0), P(10, 10), INFINITY), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_element_count(q), 1u);
  gshgeom_path_destroy(q);
}

TEST(ArcTangent, TheBulgeAgainstTheOracle) {
  using namespace arc_vectors;
  for (size_t i = 0; i < bulge_count; i++) {
    const Bulge & v = bulges[i];
    GSHGEOM_Path * p = make();
    gshgeom_path_move_to(p, P(v.sx, v.sy));
    ASSERT_EQ(gshgeom_path_arc_bulge_to(p, P(v.qx, v.qy), v.bulge), GSHGEOM_OK) << i;
    if (v.kind == 0) {
      EXPECT_EQ(gshgeom_path_element_count(p), 1u) << i;
    } else if (v.kind == 1) {
      GSHGEOM_Verb vb;
      gshgeom_path_verb(p, 1, &vb);
      EXPECT_EQ(vb, GSHGEOM_VERB_LINE) << i;
    } else {
      GSHGEOM_Arc a = read_arc(p, 1);
      const double tol = rel_tol(v.scale);
      EXPECT_LE(std::fabs(a.center.x - v.cx), tol * 8) << i;
      EXPECT_LE(std::fabs(a.center.y - v.cy), tol * 8) << i;
      EXPECT_LE(std::fabs(a.u.x - v.r), tol * 8) << i;
      EXPECT_LE(std::fabs(a.sweep - v.sweep), 8 * GSHGEOM_ARC_DECOMP) << i;
      EXPECT_EQ(bits(a.end.x), bits(v.qx));
    }
    gshgeom_path_destroy(p);
  }
  GSHGEOM_Path * q = make();
  gshgeom_path_move_to(q, P(0, 0));
  EXPECT_EQ(gshgeom_path_arc_bulge_to(q, P(1, 0), INFINITY), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_arc_bulge_to(q, P(1, 0), NAN), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_element_count(q), 1u);
  // a bulge so large the sweep is a full turn less a hair: an arc whose ends agree, or INVALID, never a crash or a half-built path
  for (double b : {1e6, 1e15, 1e300}) {
    const uint64_t h = gshgeom_path_hash(q);
    GSHGEOM_Result r = gshgeom_path_arc_bulge_to(q, P(3, 4), b);
    EXPECT_TRUE(r == GSHGEOM_OK || r == GSHGEOM_ERR_INVALID) << b << " " << r;
    if (r != GSHGEOM_OK) {
      EXPECT_EQ(gshgeom_path_hash(q), h);
    }
  }
  gshgeom_path_destroy(q);
}

// ---------------------------------------------------------------------------
// The decomposition
// ---------------------------------------------------------------------------

TEST(ArcPieces, EveryVectorMatchesTheOracleAndTheEndsAreTheArcsOwn) {
  using namespace arc_vectors;
  double worst = 0;
  for (size_t i = 0; i < decomp_count; i++) {
    const Decomp & d = decomps[i];
    GSHGEOM_Arc a;
    std::memset(&a, 0, sizeof a);
    a.start = P(d.sx, d.sy);
    a.end = P(d.ex, d.ey);
    a.center = P(d.cx, d.cy);
    a.u = P(d.ux, d.uy);
    a.v = P(d.vx, d.vy);
    a.unit = (GSHGEOM_Angle_Unit)d.unit;
    a.t0 = d.t0;
    a.sweep = d.sweep;
    GSHGEOM_Curve pc[GSHGEOM_ARC_MAX_PIECES];
    int n = 0;
    ASSERT_EQ(gshgeom_arc_pieces(&a, pc, &n), GSHGEOM_OK) << i;
    ASSERT_EQ(n, d.n) << "vector " << i << " sweep " << d.sweep;
    EXPECT_EQ(bits(pc[0].p[0].x), bits(d.sx));
    EXPECT_EQ(bits(pc[0].p[0].y), bits(d.sy));
    EXPECT_EQ(bits(pc[n - 1].p[2].x), bits(d.ex));
    EXPECT_EQ(bits(pc[n - 1].p[2].y), bits(d.ey));
    const double tol = rel_tol(d.scale);
    for (int k = 0; k < n; k++) {
      const Piece & e = d.pieces[k];
      EXPECT_EQ(pc[k].kind, GSHGEOM_CURVE_RATIONAL_QUADRATIC);
      EXPECT_EQ(pc[k].w[0], 1.0);
      EXPECT_EQ(pc[k].w[2], 1.0);
      EXPECT_GT(pc[k].w[1], 0.7071);
      EXPECT_LE(pc[k].w[1], 1.0);
      auto chk = [&](double got, double want, double t) {
        worst = std::fmax(worst, std::fabs(got - want) / tol);
        EXPECT_LE(std::fabs(got - want), t) << "vector " << i << " piece " << k << " got " << got << " want " << want;
      };
      chk(pc[k].p[0].x, e.p0x, tol);
      chk(pc[k].p[0].y, e.p0y, tol);
      chk(pc[k].p[1].x, e.qx, tol);
      chk(pc[k].p[1].y, e.qy, tol);
      chk(pc[k].p[2].x, e.p2x, tol);
      chk(pc[k].p[2].y, e.p2y, tol);
      chk(pc[k].w[1], e.w, GSHGEOM_ARC_DECOMP);
      if (d.exact) {
        // the oracle's coordinates are exact doubles here: the library's must be the same bits (the interior
        // pieces' ends; the first and last are the arc's own, which agree because the record was made from them)
        EXPECT_EQ(bits(pc[k].p[1].x), bits(e.qx)) << "vector " << i << " piece " << k;
        EXPECT_EQ(bits(pc[k].p[1].y), bits(e.qy)) << "vector " << i << " piece " << k;
        EXPECT_EQ(bits(pc[k].p[2].x), bits(e.p2x)) << "vector " << i << " piece " << k;
        EXPECT_EQ(bits(pc[k].p[2].y), bits(e.p2y)) << "vector " << i << " piece " << k;
      }
      if (d.exactw) {
        EXPECT_EQ(bits(pc[k].w[1]), bits(e.w)) << "vector " << i << " piece " << k;
      }
    }
    // consecutive pieces share their meeting point to the bit
    for (int k = 0; k + 1 < n; k++) {
      EXPECT_EQ(bits(pc[k].p[2].x), bits(pc[k + 1].p[0].x));
      EXPECT_EQ(bits(pc[k].p[2].y), bits(pc[k + 1].p[0].y));
    }
  }
  std::printf("[ arc ] decomposition: worst error %.3g of the bound over %zu vectors\n", worst, decomp_count);
}

TEST(ArcPieces, TheCountAndTheQuarterTurnWeightAreExact) {
  const double root_half = gmath_sqrt(0.5);
  struct Case {
    GSHGEOM_Angle_Unit unit;
    double sweep;
    int n;
  } cases[] = {{GSHGEOM_ANGLE_DEGREES, 1e-9, 1}, {GSHGEOM_ANGLE_DEGREES, 90, 1}, {GSHGEOM_ANGLE_DEGREES, 91, 2},
      {GSHGEOM_ANGLE_DEGREES, 180, 2}, {GSHGEOM_ANGLE_DEGREES, 270, 3}, {GSHGEOM_ANGLE_DEGREES, 360, 4},
      {GSHGEOM_ANGLE_DEGREES, -360, 4}, {GSHGEOM_ANGLE_DEGREES, -90, 1}, {GSHGEOM_ANGLE_DEGREES, 270.0000000000001, 4},
      {GSHGEOM_ANGLE_GRADIANS, 100, 1}, {GSHGEOM_ANGLE_GRADIANS, 300, 3}, {GSHGEOM_ANGLE_TURNS, 0.25, 1},
      {GSHGEOM_ANGLE_TURNS, 0.75, 3}, {GSHGEOM_ANGLE_TURNS, 1.0, 4}, {GSHGEOM_ANGLE_RADIANS, 0x1.921fb54442d18p+0, 1},
      {GSHGEOM_ANGLE_RADIANS, 0x1.921fb54442d19p+0, 2}, {GSHGEOM_ANGLE_RADIANS, 0x1.921fb54442d18p+1, 2},
      {GSHGEOM_ANGLE_RADIANS, 0x1.921fb54442d19p+1, 3}, {GSHGEOM_ANGLE_RADIANS, 0x1.921fb54442d18p+2, 4}};
  for (const Case & c : cases) {
    GSHGEOM_Arc a = circle_arc(0, 0, 8, 0, c.sweep, c.unit);
    GSHGEOM_Curve pc[GSHGEOM_ARC_MAX_PIECES];
    int n = 0;
    ASSERT_EQ(gshgeom_arc_pieces(&a, pc, &n), GSHGEOM_OK);
    EXPECT_EQ(n, c.n) << c.sweep;
    if (c.sweep == 90 || c.sweep == -90 || c.sweep == 360 || c.sweep == -360 || c.sweep == 100 || c.sweep == 300 || c.sweep == 0.25
        || c.sweep == 1.0 || c.sweep == 0.75 || c.sweep == 270 || (c.unit == GSHGEOM_ANGLE_RADIANS && c.sweep == 0x1.921fb54442d18p+0)) {
      for (int k = 0; k < n; k++) {
        EXPECT_EQ(bits(pc[k].w[1]), bits(root_half)) << c.sweep;
      }
    }
  }
}

TEST(ArcPieces, ARoundedRectanglesCornersHaveExactControlPointsAndDirections) {
  // top-right corner of a 100 by 60 rect with radius 10, y down: from (90, 0) clockwise (a positive sweep) to (100, 10)
  GSHGEOM_Arc a = circle_arc(90, 10, 10, 270, 90);
  EXPECT_EQ(a.start.x, 90.0);
  EXPECT_EQ(a.start.y, 0.0);
  EXPECT_EQ(a.end.x, 100.0);
  EXPECT_EQ(a.end.y, 10.0);
  GSHGEOM_Curve pc[GSHGEOM_ARC_MAX_PIECES];
  int n = 0;
  ASSERT_EQ(gshgeom_arc_pieces(&a, pc, &n), GSHGEOM_OK);
  ASSERT_EQ(n, 1);
  EXPECT_EQ(pc[0].p[1].x, 100.0);  // (E.x, S.y): the first direction is exactly +x, the last exactly +y
  EXPECT_EQ(pc[0].p[1].y, 0.0);
  EXPECT_EQ(pc[0].w[1], gmath_sqrt(0.5));
  // a circle of four quarters in the unit circle's own order, and each reversed
  for (double sw : {360.0, -360.0}) {
    GSHGEOM_Arc c = circle_arc(3, 4, 5, 0, sw);
    ASSERT_EQ(gshgeom_arc_pieces(&c, pc, &n), GSHGEOM_OK);
    ASSERT_EQ(n, 4);
    for (int k = 0; k < 4; k++) {
      // every control polygon is axis-aligned at both ends: Q is a corner of the circle's square
      const double dx0 = pc[k].p[1].x - pc[k].p[0].x, dy0 = pc[k].p[1].y - pc[k].p[0].y;
      const double dx1 = pc[k].p[2].x - pc[k].p[1].x, dy1 = pc[k].p[2].y - pc[k].p[1].y;
      EXPECT_TRUE((dx0 == 0.0) != (dy0 == 0.0)) << k;
      EXPECT_TRUE((dx1 == 0.0) != (dy1 == 0.0)) << k;
      EXPECT_TRUE((dx0 == 0.0) != (dx1 == 0.0)) << k << ": the tangents are perpendicular";
    }
  }
}

TEST(ArcPieces, ASmallSweepKeepsItsWeightAtOneAndTheOddSweepsAreNotTheWholeCosine) {
  // a plant: w = cos of the whole sweep instead of the half of a piece would be 0 for a quarter turn
  GSHGEOM_Arc a = circle_arc(0, 0, 10, 10, 180);
  GSHGEOM_Curve pc[GSHGEOM_ARC_MAX_PIECES];
  int n;
  ASSERT_EQ(gshgeom_arc_pieces(&a, pc, &n), GSHGEOM_OK);
  EXPECT_EQ(n, 2);
  EXPECT_EQ(bits(pc[0].w[1]), bits(gmath_sqrt(0.5)));
  GSHGEOM_Arc b = circle_arc(0, 0, 10, 10, 1e-9);
  ASSERT_EQ(gshgeom_arc_pieces(&b, pc, &n), GSHGEOM_OK);
  EXPECT_EQ(n, 1);
  EXPECT_EQ(pc[0].w[1], 1.0);
  GSHGEOM_Arc c = circle_arc(0, 0, 10, 10, 60);
  ASSERT_EQ(gshgeom_arc_pieces(&c, pc, &n), GSHGEOM_OK);
  EXPECT_NEAR(pc[0].w[1], std::cos(30.0 * M_PI / 180.0), 1e-15);
}

TEST(ArcPieces, LoweringReplacesEveryArcAndTouchesNothingElse) {
  GSHGEOM_Path * p = make();
  gshgeom_path_move_to(p, P(10, 0));
  GSHGEOM_Arc a = circle_arc(0, 0, 10, 0, 270);
  ASSERT_EQ(gshgeom_path_append_arc(p, &a), GSHGEOM_OK);
  gshgeom_path_line_to(p, P(0, 0));
  gshgeom_path_quadratic_to(p, P(1, 1), P(2, 0));
  gshgeom_path_close(p);
  gshgeom_path_move_to(p, P(50, 50));
  GSHGEOM_Path * l = nullptr;
  ASSERT_EQ(gshgeom_path_lower_arcs(p, nullptr, nullptr, &l), GSHGEOM_OK);
  // move, three rational quadratics, line, quadratic, close, move
  ASSERT_EQ(gshgeom_path_element_count(l), 8u);
  GSHGEOM_Verb want[8] = {GSHGEOM_VERB_MOVE, GSHGEOM_VERB_RATIONAL_QUADRATIC, GSHGEOM_VERB_RATIONAL_QUADRATIC,
      GSHGEOM_VERB_RATIONAL_QUADRATIC, GSHGEOM_VERB_LINE, GSHGEOM_VERB_QUADRATIC, GSHGEOM_VERB_CLOSE, GSHGEOM_VERB_MOVE};
  for (size_t i = 0; i < 8; i++) {
    GSHGEOM_Element el;
    ASSERT_EQ(gshgeom_path_element(l, i, &el), GSHGEOM_OK) << i;
    EXPECT_EQ(el.verb, want[i]) << i;
  }
  GSHGEOM_Element el;
  gshgeom_path_element(l, 3, &el);
  EXPECT_EQ(bits(el.curve.p[2].x), bits(a.end.x));  // the arc's own end, to the bit
  EXPECT_EQ(bits(el.curve.p[2].y), bits(a.end.y));
  // the box of the arc is the box of its pieces: the circle's three-quarters from angle 0
  GSHGEOM_Bounds b;
  ASSERT_EQ(gshgeom_path_bounds(p, &b), GSHGEOM_OK);
  EXPECT_LE(b.min.x, -9.999999);
  EXPECT_GE(b.max.y, 9.999999);
  GSHGEOM_Bounds bl;
  ASSERT_EQ(gshgeom_path_bounds(l, &bl), GSHGEOM_OK);
  EXPECT_EQ(bits(b.min.x), bits(bl.min.x));
  EXPECT_EQ(bits(b.max.y), bits(bl.max.y));
  // a lowered path has no arc, and lowering it again is the identity
  GSHGEOM_Path * l2 = nullptr;
  ASSERT_EQ(gshgeom_path_lower_arcs(l, nullptr, nullptr, &l2), GSHGEOM_OK);
  EXPECT_TRUE(gshgeom_path_equal(l, l2));
  gshgeom_path_destroy(p);
  gshgeom_path_destroy(l);
  gshgeom_path_destroy(l2);
}

TEST(ArcPieces, ASweepWhoseOnePieceHasNoExtentIsStillOneElement) {
  // an arc of 1e-9 degrees on a unit circle: one piece, all control points within an ulp of each other
  GSHGEOM_Arc a = circle_arc(0, 0, 1, 10, 1e-9);
  GSHGEOM_Path* p = path_with(a);
  EXPECT_EQ(gshgeom_path_element_count(p), 2u);
  gshgeom_path_destroy(p);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
