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
#include "../shapes_vectors.h"

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

GSHGEOM_Path * make(int dim = 2, const GSHGEOM_Limits * lim = nullptr, const GSHGEOM_Allocator * al = nullptr) {
  GSHGEOM_Path * p = nullptr;
  EXPECT_EQ(gshgeom_path_create(al, lim, dim, &p), GSHGEOM_OK);
  return p;
}

bool same(double a, double b) {
  return std::memcmp(&a, &b, sizeof a) == 0;
}

/* A piece's control point is compared as a value: A1's decomposition (AR8) does not promise the sign of a zero it derives. */
bool value_pt(const GSHGEOM_Point & p, const double e[2]) {
  return p.x == e[0] && p.y == e[1];
}

bool same_pt(const GSHGEOM_Point & p, const double e[2]) {
  return same(p.x, e[0]) && same(p.y, e[1]);
}

GSHGEOM_Result build(GSHGEOM_Path * p, const shapes_vectors::Case & c) {
  using namespace shapes_vectors;
  if (c.kind == 0) {
    return gshgeom_path_add_rect(p, c.a[0], c.a[1], c.a[2], c.a[3], c.fx, c.fy);
  }
  if (c.kind == 2) {
    return gshgeom_path_add_ellipse(p, c.a[0], c.a[1], c.a[2], c.a[3], c.fx, c.fy);
  }
  GSHGEOM_Corner r[4];
  for (int i = 0; i < 4; i++) {
    r[i].rx = c.r[2 * i];
    r[i].ry = c.r[2 * i + 1];
  }
  return gshgeom_path_add_rounded_rect(p, c.a[0], c.a[1], c.a[2], c.a[3], r);
}

/* Compare the elements from index @p at with the oracle's, bit for bit; returns a message or an empty string. */
std::string compare(const GSHGEOM_Path * p, size_t at, const shapes_vectors::Case & c) {
  using namespace shapes_vectors;
  char buf[200];
  if (gshgeom_path_element_count(p) != at + c.count) {
    std::snprintf(buf, sizeof buf, "element count %zu, expected %zu", gshgeom_path_element_count(p), at + c.count);
    return buf;
  }
  for (size_t i = 0; i < c.count; i++) {
    const Elem & e = elems[c.first + i];
    GSHGEOM_Verb verb;
    if (gshgeom_path_verb(p, at + i, &verb) != GSHGEOM_OK) {
      return "verb unreadable";
    }
    const GSHGEOM_Verb want = e.verb == 'M' ? GSHGEOM_VERB_MOVE : e.verb == 'L' ? GSHGEOM_VERB_LINE : e.verb == 'A' ? GSHGEOM_VERB_ARC : GSHGEOM_VERB_CLOSE;
    if (verb != want) {
      std::snprintf(buf, sizeof buf, "element %zu: verb %d, expected %d", i, (int)verb, (int)want);
      return buf;
    }
    if (e.verb == 'A') {
      GSHGEOM_Arc a;
      if (gshgeom_path_arc(p, at + i, &a) != GSHGEOM_OK) {
        return "arc unreadable";
      }
      const bool ok = same_pt(a.start, e.s) && same_pt(a.end, e.e) && same_pt(a.center, e.c) && same_pt(a.u, e.u) &&
                      same_pt(a.v, e.v) && same(a.t0, e.t0) && same(a.sweep, e.sweep) && a.unit == GSHGEOM_ANGLE_DEGREES &&
                      ((a.flags & GSHGEOM_ARC_CIRCULAR) != 0u) == (e.circ != 0);
      if (!ok) {
        std::snprintf(buf, sizeof buf, "element %zu: the arc differs from the oracle's", i);
        return buf;
      }
      GSHGEOM_Curve pc[GSHGEOM_ARC_MAX_PIECES];
      int n = 0;
      if (gshgeom_arc_pieces(&a, pc, &n) != GSHGEOM_OK || n != 1) {
        return "a quarter turn is not one piece";
      }
      if (!same_pt(pc[0].p[0], e.s) || !same_pt(pc[0].p[2], e.e) || !value_pt(pc[0].p[1], e.ctl)) {
        std::snprintf(buf, sizeof buf, "element %zu: the control polygon is not the exact one", i);
        return buf;
      }
    } else {
      GSHGEOM_Element el;
      if (gshgeom_path_element(p, at + i, &el) != GSHGEOM_OK) {
        return "element unreadable";
      }
      const GSHGEOM_Point end = el.curve.p[el.curve.kind == GSHGEOM_CURVE_LINE ? 1 : 0];
      if (!same_pt(end, e.e)) {
        std::snprintf(buf, sizeof buf, "element %zu (%c): the point differs from the oracle's (%a, %a)", i, e.verb, end.x, end.y);
        return buf;
      }
    }
  }
  return "";
}

}  // namespace

TEST(Shapes, EveryCaseEqualsTheOraclesExactPathToTheBit) {
  using namespace shapes_vectors;
  size_t ok = 0, refused = 0;
  for (size_t k = 0; k < case_count; k++) {
    const Case & c = cases[k];
    for (int prior = 0; prior < 2; prior++) {
      GSHGEOM_Path * p = make();
      if (prior) {
        ASSERT_EQ(gshgeom_path_move_to(p, P(-1, -2)), GSHGEOM_OK);
        ASSERT_EQ(gshgeom_path_line_to(p, P(-3, 4)), GSHGEOM_OK);
      }
      const size_t at = gshgeom_path_element_count(p);
      EXPECT_EQ(build(p, c), c.result) << "case " << k;
      if (c.result == GSHGEOM_OK) {
        const std::string m = compare(p, at, c);
        EXPECT_EQ(m, "") << "case " << k << " kind " << c.kind << " flips " << c.fx << c.fy;
        ok += prior == 0;
      } else {
        EXPECT_EQ(gshgeom_path_element_count(p), at) << "a refusal leaves the path as it was, case " << k;
        refused += prior == 0;
      }
      gshgeom_path_destroy(p);
    }
  }
  EXPECT_GT(ok, 600u);
  EXPECT_GT(refused, 60u);
}

TEST(Shapes, TheFourSignCasesOfCanvasRectAreTheFourFlipPairs) {
  /* Canvas rect(x, y, w, h): (x,y) (x+w,y) (x+w,y+h) (x,y+h) taken literally; the box is x=10, y=20, w=100, h=50.  Written out by hand. */
  struct Row {
    int fx, fy;
    double sx, sy;  // the start corner
    double nx, ny;  // the second point
    int clockwise;
  };
  const Row rows[] = {{0, 0, 10, 20, 110, 20, 1}, {1, 0, 110, 20, 10, 20, 0}, {0, 1, 10, 70, 110, 70, 0}, {1, 1, 110, 70, 10, 70, 1}};
  for (const Row & r : rows) {
    GSHGEOM_Path * p = make();
    ASSERT_EQ(gshgeom_path_add_rect(p, 10, 20, 100, 50, r.fx, r.fy), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_path_element_count(p), 5u);
    GSHGEOM_Element e0, e1, e2;
    ASSERT_EQ(gshgeom_path_element(p, 0, &e0), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_path_element(p, 1, &e1), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_path_element(p, 2, &e2), GSHGEOM_OK);
    EXPECT_EQ(e0.curve.p[0].x, r.sx) << r.fx << r.fy;
    EXPECT_EQ(e0.curve.p[0].y, r.sy) << r.fx << r.fy;
    EXPECT_EQ(e1.curve.p[1].x, r.nx) << r.fx << r.fy;
    EXPECT_EQ(e1.curve.p[1].y, r.ny) << r.fx << r.fy;
    /* the turn: the sign of the cross product at the second corner, in y-down screen terms positive is clockwise */
    const double cr = (e1.curve.p[1].x - e0.curve.p[0].x) * (e2.curve.p[1].y - e1.curve.p[1].y) -
                      (e1.curve.p[1].y - e0.curve.p[0].y) * (e2.curve.p[1].x - e1.curve.p[1].x);
    EXPECT_EQ(cr > 0.0, r.clockwise != 0) << r.fx << r.fy;
    gshgeom_path_destroy(p);
  }
}

TEST(Shapes, TheEllipseStartsAndTurnsAsCanvasEllipseDoes) {
  /* flips -> (start angle, positive turn?): none 0 yes; flip_y 0 no; flip_x pi no; both pi yes */
  struct Row {
    int fx, fy;
    double sx;
    double first_y;
    int positive;
  };
  const Row rows[] = {{0, 0, 150, 120, 1}, {0, 1, 150, 80, 0}, {1, 0, 50, 120, 0}, {1, 1, 50, 80, 1}};
  for (const Row & r : rows) {
    GSHGEOM_Path * p = make();
    ASSERT_EQ(gshgeom_path_add_ellipse(p, 100, 100, 50, 20, r.fx, r.fy), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_path_element_count(p), 6u);
    GSHGEOM_Arc a;
    ASSERT_EQ(gshgeom_path_arc(p, 1, &a), GSHGEOM_OK);
    EXPECT_EQ(a.start.x, r.sx);
    EXPECT_EQ(a.start.y, 100.0);
    EXPECT_EQ(a.end.x, 100.0);
    EXPECT_EQ(a.end.y, r.first_y);
    const double cross = a.u.x * a.v.y - a.u.y * a.v.x;
    EXPECT_EQ(cross > 0.0, r.positive != 0) << r.fx << r.fy;
    EXPECT_EQ(a.sweep, 90.0);
    EXPECT_EQ(a.t0, 0.0);
    gshgeom_path_destroy(p);
  }
}

TEST(Shapes, TheDirectionIsTheParityOfTheFlips) {
  for (int fx = 0; fx < 2; fx++) {
    for (int fy = 0; fy < 2; fy++) {
      GSHGEOM_Path * p = make();
      ASSERT_EQ(gshgeom_path_add_ellipse(p, 0, 0, 3, 3, fx, fy), GSHGEOM_OK);
      for (size_t i = 1; i <= 4; i++) {
        GSHGEOM_Arc a;
        ASSERT_EQ(gshgeom_path_arc(p, i, &a), GSHGEOM_OK);
        const double cross = a.u.x * a.v.y - a.u.y * a.v.x;
        EXPECT_EQ(cross < 0.0, fx != fy);
        EXPECT_EQ(cross > 0.0, fx == fy);
        EXPECT_EQ(std::fabs(a.u.x), 3.0);
        EXPECT_EQ(std::fabs(a.v.y), 3.0);
        EXPECT_EQ(a.u.x < 0.0, fx == 1);
        EXPECT_EQ(a.v.y < 0.0, fy == 1);
      }
      gshgeom_path_destroy(p);
    }
  }
}

TEST(Shapes, ARoundedRectStartsAfterTheTopLeftRadiusAndTurnsClockwise) {
  GSHGEOM_Path * p = make();
  const GSHGEOM_Corner r[4] = {{10, 12}, {6, 5}, {4, 4}, {8, 9}};
  ASSERT_EQ(gshgeom_path_add_rounded_rect(p, 0, 0, 100, 60, r), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_element_count(p), 10u);
  const GSHGEOM_Verb verbs[10] = {GSHGEOM_VERB_MOVE, GSHGEOM_VERB_LINE, GSHGEOM_VERB_ARC, GSHGEOM_VERB_LINE, GSHGEOM_VERB_ARC,
      GSHGEOM_VERB_LINE, GSHGEOM_VERB_ARC, GSHGEOM_VERB_LINE, GSHGEOM_VERB_ARC, GSHGEOM_VERB_CLOSE};
  for (size_t i = 0; i < 10; i++) {
    GSHGEOM_Verb v;
    ASSERT_EQ(gshgeom_path_verb(p, i, &v), GSHGEOM_OK);
    EXPECT_EQ(v, verbs[i]) << i;
  }
  GSHGEOM_Element m;
  ASSERT_EQ(gshgeom_path_element(p, 0, &m), GSHGEOM_OK);
  EXPECT_EQ(m.curve.p[0].x, 10.0);
  EXPECT_EQ(m.curve.p[0].y, 0.0);
  /* the arcs in order TR, BR, BL, TL, each with the starts the rule names */
  const double t0s[4] = {-90, 0, 90, 180};
  const GSHGEOM_Point starts[4] = {P(94, 0), P(100, 56), P(8, 60), P(0, 12)};
  const GSHGEOM_Point ends[4] = {P(100, 5), P(96, 60), P(0, 51), P(10, 0)};
  const GSHGEOM_Point centres[4] = {P(94, 5), P(96, 56), P(8, 51), P(10, 12)};
  for (size_t k = 0; k < 4; k++) {
    GSHGEOM_Arc a;
    ASSERT_EQ(gshgeom_path_arc(p, 2 + 2 * k, &a), GSHGEOM_OK);
    EXPECT_EQ(a.t0, t0s[k]);
    EXPECT_EQ(a.sweep, 90.0);
    EXPECT_EQ(a.start.x, starts[k].x);
    EXPECT_EQ(a.start.y, starts[k].y);
    EXPECT_EQ(a.end.x, ends[k].x);
    EXPECT_EQ(a.end.y, ends[k].y);
    EXPECT_EQ(a.center.x, centres[k].x);
    EXPECT_EQ(a.center.y, centres[k].y);
  }
  gshgeom_path_destroy(p);
}

TEST(Shapes, TheFourLinesAreEmittedEvenWhenTheRadiiSumToTheSide) {
  GSHGEOM_Path * p = make();
  const GSHGEOM_Corner r[4] = {{10, 5}, {10, 5}, {10, 5}, {10, 5}};  // a stadium: both sums are the side
  ASSERT_EQ(gshgeom_path_add_rounded_rect(p, 0, 0, 20, 10, r), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_element_count(p), 10u);
  GSHGEOM_Element e;
  ASSERT_EQ(gshgeom_path_element(p, 1, &e), GSHGEOM_OK);  // the top line, from (10, 0) to (10, 0)
  EXPECT_EQ(e.curve.p[0].x, e.curve.p[1].x);
  EXPECT_EQ(e.curve.p[0].y, e.curve.p[1].y);
  gshgeom_path_destroy(p);
}

TEST(Shapes, ACornerWithAZeroComponentIsSquareAndAllSquareIsTheRect) {
  GSHGEOM_Path * a = make();
  GSHGEOM_Path * b = make();
  GSHGEOM_Path * c = make();
  const GSHGEOM_Corner half[4] = {{5, 0}, {0, 5}, {0, 0}, {5, 0}};
  const GSHGEOM_Corner none[4] = {{0, 0}, {0, 0}, {0, 0}, {0, 0}};
  ASSERT_EQ(gshgeom_path_add_rounded_rect(a, 1, 2, 30, 40, half), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_add_rounded_rect(b, 1, 2, 30, 40, none), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_add_rect(c, 1, 2, 30, 40, 0, 0), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_element_count(a), 5u);
  EXPECT_TRUE(gshgeom_path_equal(a, b));
  EXPECT_TRUE(gshgeom_path_equal(b, c));
  gshgeom_path_destroy(a);
  gshgeom_path_destroy(b);
  gshgeom_path_destroy(c);
}

TEST(Shapes, OverlapIsDecidedOnTheExactSumNotTheRoundedOne) {
  /* 0.1 + 0.9 in the reals is a little over 1, and in double it rounds to exactly 1.0 */
  ASSERT_EQ(0.1 + 0.9, 1.0);
  const GSHGEOM_Corner over[4] = {{0.1, 0.125}, {0.9, 0.125}, {0, 0}, {0, 0}};
  EXPECT_EQ(gshgeom_rounded_rect_check(1.0, 10.0, over), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_rounded_rect_check(std::nextafter(1.0, 2.0), 10.0, over), GSHGEOM_OK);
  const GSHGEOM_Corner exact[4] = {{0.5, 0.25}, {0.5, 0.25}, {0, 0}, {0, 0}};
  EXPECT_EQ(gshgeom_rounded_rect_check(1.0, 10.0, exact), GSHGEOM_OK);
  const GSHGEOM_Corner one_over[4] = {{std::nextafter(0.5, 1.0), 0.25}, {0.5, 0.25}, {0, 0}, {0, 0}};
  EXPECT_EQ(gshgeom_rounded_rect_check(1.0, 10.0, one_over), GSHGEOM_ERR_INVALID);
  /* the same on each of the four sides (indices: top 0+1 in x, right 1+2 in y, bottom 2+3 in x, left 3+0 in y) */
  const GSHGEOM_Corner right[4] = {{0, 0}, {0.125, 0.1}, {0.125, 0.9}, {0, 0}};
  EXPECT_EQ(gshgeom_rounded_rect_check(10.0, 1.0, right), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_rounded_rect_check(10.0, std::nextafter(1.0, 2.0), right), GSHGEOM_OK);
  const GSHGEOM_Corner bottom[4] = {{0, 0}, {0, 0}, {0.1, 0.125}, {0.9, 0.125}};
  EXPECT_EQ(gshgeom_rounded_rect_check(1.0, 10.0, bottom), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_rounded_rect_check(std::nextafter(1.0, 2.0), 10.0, bottom), GSHGEOM_OK);
  const GSHGEOM_Corner left[4] = {{0.125, 0.9}, {0, 0}, {0, 0}, {0.125, 0.1}};
  EXPECT_EQ(gshgeom_rounded_rect_check(10.0, 1.0, left), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_rounded_rect_check(10.0, std::nextafter(1.0, 2.0), left), GSHGEOM_OK);
  /* a sum whose rounded value is under the side while the reals are over */
  const double big = 9007199254740992.0;  // 2^53
  const GSHGEOM_Corner huge[4] = {{big - 2.0, 1.0}, {3.0, 1.0}, {0, 0}, {0, 0}};
  EXPECT_EQ(gshgeom_rounded_rect_check(big, 8.0, huge), GSHGEOM_ERR_INVALID);
  const GSHGEOM_Corner huge2[4] = {{big - 2.0, 1.0}, {2.0, 1.0}, {0, 0}, {0, 0}};
  EXPECT_EQ(gshgeom_rounded_rect_check(big, 8.0, huge2), GSHGEOM_OK);
}

TEST(Shapes, ASquareCornerDoesNotCountTowardAnOverlap) {
  const GSHGEOM_Corner c[4] = {{30, 0}, {30, 5}, {0, 0}, {0, 0}};
  EXPECT_EQ(gshgeom_rounded_rect_check(40.0, 40.0, c), GSHGEOM_OK);
  const GSHGEOM_Corner d[4] = {{30, 1}, {30, 5}, {0, 0}, {0, 0}};
  EXPECT_EQ(gshgeom_rounded_rect_check(40.0, 40.0, d), GSHGEOM_ERR_INVALID);
}

TEST(Shapes, ArgumentErrorsAreInvalidInThePublishedOrder) {
  GSHGEOM_Path * p = make();
  GSHGEOM_Path * p3 = make(3);
  const GSHGEOM_Corner z[4] = {{0, 0}, {0, 0}, {0, 0}, {0, 0}};
  EXPECT_EQ(gshgeom_path_add_rect(nullptr, 0, 0, 1, 1, 0, 0), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_add_rect(p3, 0, 0, 1, 1, 0, 0), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_add_rect(p, 0, 0, 1, 1, 2, 0), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_add_rect(p, 0, 0, 1, 1, 0, -1), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_add_ellipse(p, 0, 0, 1, 1, 2, 0), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_add_ellipse(p, 0, 0, 1, 1, 0, 2), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_add_ellipse(p3, 0, 0, 1, 1, 0, 0), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_add_ellipse(nullptr, 0, 0, 1, 1, 0, 0), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_add_rounded_rect(nullptr, 0, 0, 1, 1, z), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_add_rounded_rect(p3, 0, 0, 1, 1, z), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_add_rounded_rect(p, 0, 0, 1, 1, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_rounded_rect_check(1, 1, nullptr), GSHGEOM_ERR_INVALID);
  /* an invalid number is reported before a limit: a NaN beside an out-of-range number */
  EXPECT_EQ(gshgeom_path_add_rect(p, std::nan(""), 0x1p101, 1, 1, 0, 0), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_add_rect(p, 0x1p101, 0, -1, 1, 0, 0), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_add_ellipse(p, 0x1p101, 0, 0, 1, 0, 0), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_add_ellipse(p, 0, 0, std::nan(""), 0x1p101, 0, 0), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_element_count(p), 0u);
  gshgeom_path_destroy(p);
  gshgeom_path_destroy(p3);
}

TEST(Shapes, TheCoordinateLimitIsAtTheDerivedCoordinateAndOneUlpPastIsRefused) {
  GSHGEOM_Path * p = make();
  const double lim = 0x1p100;
  EXPECT_EQ(gshgeom_path_add_rect(p, lim / 2, 0, lim / 2, 1, 0, 0), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_add_rect(p, lim / 2, 0, lim / 2 + 0x1p48, 1, 0, 0), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(gshgeom_path_add_rect(p, lim / 2, 0, lim / 2 + 0x1p47, 1, 0, 0), GSHGEOM_OK) << "x + w is the limit once rounded";
  EXPECT_EQ(gshgeom_path_add_rect(p, std::nextafter(lim, 0x1p200), 0, 1, 1, 0, 0), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(gshgeom_path_add_ellipse(p, 0, 0, lim, 1, 0, 0), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_add_ellipse(p, 0, 0, std::nextafter(lim, 0x1p200), 1, 0, 0), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(gshgeom_path_add_ellipse(p, lim, 0, 0x1p48, 1, 0, 0), GSHGEOM_ERR_LIMIT) << "cx + rx passes the limit";
  EXPECT_EQ(gshgeom_path_add_ellipse(p, -lim, 0, 0x1p48, 1, 1, 0), GSHGEOM_ERR_LIMIT) << "flip_x: cx - rx passes it the other way";
  EXPECT_EQ(gshgeom_path_element_count(p), 5u + 5u + 6u);
  gshgeom_path_destroy(p);
}

TEST(Shapes, TheElementLimitIsTheBoundAndOnePastTakesTheCallBack) {
  struct Row {
    int kind;
    size_t n;
  };
  const Row rows[] = {{0, GSHGEOM_SHAPE_RECT_ELEMENTS}, {1, GSHGEOM_SHAPE_ROUNDED_RECT_ELEMENTS}, {2, GSHGEOM_SHAPE_ELLIPSE_ELEMENTS}};
  const GSHGEOM_Corner r[4] = {{2, 2}, {2, 2}, {2, 2}, {2, 2}};
  for (const Row & row : rows) {
    for (size_t extra = 0; extra < 2; extra++) {
      GSHGEOM_Limits lim;
      gshgeom_limits_default(&lim);
      lim.max_elements = row.n - 1 + extra;
      GSHGEOM_Path * p = make(2, &lim);
      ASSERT_EQ(gshgeom_path_move_to(p, P(0, 0)), GSHGEOM_OK);
      const size_t before = gshgeom_path_element_count(p);
      lim.max_elements = before + row.n - 1 + extra;
      gshgeom_path_destroy(p);
      p = make(2, &lim);
      ASSERT_EQ(gshgeom_path_move_to(p, P(0, 0)), GSHGEOM_OK);
      GSHGEOM_Result res = row.kind == 0   ? gshgeom_path_add_rect(p, 0, 0, 10, 10, 0, 0)
                           : row.kind == 1 ? gshgeom_path_add_rounded_rect(p, 0, 0, 10, 10, r)
                                           : gshgeom_path_add_ellipse(p, 0, 0, 5, 5, 0, 0);
      if (extra == 0) {
        EXPECT_EQ(res, GSHGEOM_ERR_LIMIT) << row.kind;
        EXPECT_EQ(gshgeom_path_element_count(p), before) << "a refusal takes back what it added";
      } else {
        EXPECT_EQ(res, GSHGEOM_OK) << row.kind;
        EXPECT_EQ(gshgeom_path_element_count(p), before + row.n);
      }
      gshgeom_path_destroy(p);
    }
  }
}

TEST(Shapes, AFailureAtEveryAllocationLeavesThePathAsItWasAndLeaksNothing) {
  const GSHGEOM_Corner r[4] = {{2, 2}, {3, 3}, {2, 2}, {3, 3}};
  for (int kind = 0; kind < 3; kind++) {
    bool saw_oom = false;
    for (long k = 0; k < 60; k++) {
      tracker::Tracker t;
      t.fail_from = k;
      GSHGEOM_Allocator al = tracker::allocator_for(&t);
      GSHGEOM_Path * p = nullptr;
      if (gshgeom_path_create(&al, nullptr, 2, &p) != GSHGEOM_OK) {
        EXPECT_TRUE(t.live.empty());
        continue;
      }
      GSHGEOM_Result res = kind == 0   ? gshgeom_path_add_rect(p, 0, 0, 10, 10, 1, 1)
                           : kind == 1 ? gshgeom_path_add_rounded_rect(p, 0, 0, 20, 20, r)
                                       : gshgeom_path_add_ellipse(p, 0, 0, 5, 6, 0, 1);
      if (res == GSHGEOM_ERR_OOM) {
        saw_oom = true;
        EXPECT_EQ(gshgeom_path_element_count(p), 0u) << "kind " << kind << " fail_from " << k;
      } else {
        EXPECT_EQ(res, GSHGEOM_OK);
      }
      gshgeom_path_destroy(p);
      EXPECT_TRUE(t.live.empty()) << "kind " << kind << " fail_from " << k;
    }
    EXPECT_TRUE(saw_oom) << kind;
  }
}

TEST(Shapes, ThePathsAreClosedToTheBitAndAnArcFreeShapeIsJustLines) {
  GSHGEOM_Path * p = make();
  const GSHGEOM_Corner r[4] = {{0.1, 0.3}, {0.7, 0.2}, {0.3, 0.1}, {0.05, 0.15}};
  ASSERT_EQ(gshgeom_path_add_rounded_rect(p, 0.3, 0.7, 11.3, 5.9, r), GSHGEOM_OK);
  /* every arc's start is the end of the element before it, to the bit, and the last arc ends at the move */
  GSHGEOM_Point prev = {0, 0, 0};
  GSHGEOM_Point first = {0, 0, 0};
  for (size_t i = 0; i < gshgeom_path_element_count(p); i++) {
    GSHGEOM_Verb v;
    ASSERT_EQ(gshgeom_path_verb(p, i, &v), GSHGEOM_OK);
    if (v == GSHGEOM_VERB_ARC) {
      GSHGEOM_Arc a;
      ASSERT_EQ(gshgeom_path_arc(p, i, &a), GSHGEOM_OK);
      EXPECT_TRUE(same(a.start.x, prev.x) && same(a.start.y, prev.y));
      prev = a.end;
    } else {
      GSHGEOM_Element e;
      ASSERT_EQ(gshgeom_path_element(p, i, &e), GSHGEOM_OK);
      prev = e.curve.p[e.curve.kind == GSHGEOM_CURVE_LINE ? 1 : 0];
      if (i == 0) first = prev;
    }
  }
  EXPECT_TRUE(same(prev.x, first.x) && same(prev.y, first.y));
  gshgeom_path_destroy(p);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
