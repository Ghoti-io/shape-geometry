/**
 * @file
 *
 * The dash positions (stroke.h, SG21, SG29) and the DXF lowering (SG30).
 *
 * The positions are compared bit for bit with the exact model of
 * tools/gen-dash-vectors.py (SVG 2's text, transcribed with a running position in
 * rational arithmetic), over every row of its matrix; the phase's reduction is
 * checked against exact arithmetic on its own; the lowering is compared with the
 * vectors and with its definition (the pen-down intervals of the raw pattern).
 *
 * Copyright 2026 by Corey Pennycuff
 */

#include "test_helpers.h"

#include "../dash_vectors.h"

#include <cmath>
#include <vector>

namespace {

struct Pair {
  double a, b;
};

GSHGEOM_Result collect(void * ctx, double a, double b) {
  static_cast<std::vector<Pair> *>(ctx)->push_back({a, b});
  return GSHGEOM_OK;
}

std::vector<Pair> positions(const std::vector<double> & d, double offset, double length, GSHGEOM_Result * rr = nullptr) {
  std::vector<Pair> out;
  GSHGEOM_Result r = gshgeom_stroke_dash_each(d.empty() ? nullptr : d.data(), d.size(), offset, length,
      GSHGEOM_STROKE_DASH_STEPS_DEFAULT, collect, &out);
  if (rr) {
    *rr = r;
  } else {
    EXPECT_EQ(r, GSHGEOM_OK);
  }
  return out;
}

} // namespace

TEST(DashPositions, EveryRowOfTheExactModel) {
  EXPECT_GT(dash_vectors::positions_count, (size_t)100);
  for (size_t i = 0; i < dash_vectors::positions_count; i++) {
    const dash_vectors::Positions & v = dash_vectors::positions[i];
    std::vector<double> d(v.dashes, v.dashes + v.count);
    std::vector<Pair> got = positions(d, v.offset, v.length);
    ASSERT_EQ(got.size(), v.pair_count) << "row " << i;
    for (size_t k = 0; k < got.size(); k++) {
      EXPECT_EQ(got[k].a, v.pairs[k][0]) << "row " << i << " pair " << k;
      EXPECT_EQ(got[k].b, v.pairs[k][1]) << "row " << i << " pair " << k;
    }
  }
}

TEST(DashPositions, ChosenCasesByHand) {
  /* [5,3] on 20: 0-5, 8-13, 16-20 */
  auto p = positions({5, 3}, 0, 20);
  ASSERT_EQ(p.size(), 3u);
  EXPECT_EQ(p[0].a, 0); EXPECT_EQ(p[0].b, 5);
  EXPECT_EQ(p[1].a, 8); EXPECT_EQ(p[1].b, 13);
  EXPECT_EQ(p[2].a, 16); EXPECT_EQ(p[2].b, 20);
  /* an odd array repeats: [5,3,2] is [5,3,2,5,3,2]: dashes 0-5, 8-10, 15-18, 23-25 (gaps 3, 5, 5, 3, 5, ...) */
  p = positions({5, 3, 2}, 0, 30);
  ASSERT_GE(p.size(), 4u);
  EXPECT_EQ(p[1].a, 8); EXPECT_EQ(p[1].b, 10);
  EXPECT_EQ(p[2].a, 15); EXPECT_EQ(p[2].b, 18);
  /* a dotted line: dots at 0, 10, 20 and not at 30 (a dash may not begin at the end) */
  p = positions({0, 10}, 0, 30);
  ASSERT_EQ(p.size(), 3u);
  EXPECT_EQ(p[2].a, 20); EXPECT_EQ(p[2].b, 20);
  /* a negative offset counts back from the period's end: -3 on [5,3] is 5 into the period, exactly the end of the
   * dash, which the text's rule leaves as a dash of no length at the start; the gap runs 0-3 and the next dash
   * begins at 3 */
  p = positions({5, 3}, -3, 8);
  ASSERT_EQ(p.size(), 2u);
  EXPECT_EQ(p[0].a, 0);
  EXPECT_EQ(p[0].b, 0);
  EXPECT_EQ(p[1].a, 3);
  EXPECT_EQ(p[1].b, 8);
  /* all zero, a negative length, an empty array: solid */
  for (std::vector<double> d : {std::vector<double>{0, 0}, std::vector<double>{1, -1}, std::vector<double>{}}) {
    p = positions(d, 3, 9);
    ASSERT_EQ(p.size(), 1u);
    EXPECT_EQ(p[0].a, 0);
    EXPECT_EQ(p[0].b, 9);
  }
}

TEST(DashPositions, ThePhaseAtTheEndOfADash) {
  /* SVG 2's rule, read literally: the least index whose prefix sum reaches the phase, so a phase exactly at the end of a
   * dash is in that dash, with nothing left of it: a dash of no length at the start, then the gap 0-3 and the next
   * dash 3-8 */
  auto p = positions({5, 3}, 5, 20);
  ASSERT_GE(p.size(), 2u);
  EXPECT_EQ(p[0].a, 0);
  EXPECT_EQ(p[0].b, 0);
  EXPECT_EQ(p[1].a, 3);
  EXPECT_EQ(p[1].b, 8);
  /* a dotted line starts with a dot */
  p = positions({0, 10}, 0, 5);
  ASSERT_EQ(p.size(), 1u);
  EXPECT_EQ(p[0].a, 0);
  EXPECT_EQ(p[0].b, 0);
  /* the phase at the end of a gap is the start of the next dash: no dash of no length */
  p = positions({5, 3}, 8, 5);
  ASSERT_EQ(p.size(), 1u);
  EXPECT_EQ(p[0].a, 0);
  EXPECT_EQ(p[0].b, 5);
  /* inside a gap */
  p = positions({5, 3}, 6, 1);
  EXPECT_EQ(p.size(), 0u);
  /* a length of zero: the first piece only, when the pattern is on at the phase */
  p = positions({5, 3}, 0, 0);
  ASSERT_EQ(p.size(), 1u);
  EXPECT_EQ(p[0].b, 0);
  p = positions({5, 3}, 6, 0);
  EXPECT_EQ(p.size(), 0u);
  /* a dash that would begin exactly at the end of the subpath is not made (while position < pathlength) */
  p = positions({5, 5}, 0, 10);
  ASSERT_EQ(p.size(), 1u);
  EXPECT_EQ(p[0].b, 5);
}

TEST(DashPositions, TheReductionIsExact) {
  /* a huge phase against a small period: the position is the exact remainder, which a quotient in doubles loses */
  std::vector<double> d = {0.75, 0.25};
  double off = 0x1p51 + 0.5; /* representable: the ulp at 2^51 is 1/2 */
  auto p = positions(d, off, 4);
  /* (2^51 + 0.5) mod 1 = 0.5: the phase is 0.5, inside the dash [0, 0.75): 0.25 of it remains */
  ASSERT_GE(p.size(), 1u);
  EXPECT_EQ(p[0].a, 0);
  EXPECT_EQ(p[0].b, 0.25);
  p = positions(d, -0x1p51 - 0.5, 4); /* the floored modulo: 1 - 0.5 = 0.5 */
  ASSERT_GE(p.size(), 1u);
  EXPECT_EQ(p[0].b, 0.25);
  /* a period that is not a power of two: 0.1 + 0.2, and a phase that is a multiple of it in the doubles' sense */
  std::vector<double> e = {0.1, 0.2};
  double period = 0.1 + 0.2;
  p = positions(e, 40 * period, 1); /* 40 periods: the remainder is 0 up to rounding of the product */
  ASSERT_GE(p.size(), 1u);
}

TEST(DashPositions, TheErrorDoesNotGrowWithTheDashesBefore) {
  /* 1e6 periods of [0.5, 0.5]: the dash k begins at exactly k (the running sum would drift) */
  std::vector<Pair> out;
  std::vector<double> d = {0.5, 0.5};
  struct Ctx { std::vector<Pair> * v; size_t n; } c = {&out, 0};
  auto cb = [](void * ctx, double a, double b) -> GSHGEOM_Result {
    Ctx * k = static_cast<Ctx *>(ctx);
    if (k->n % 100000 == 0) {
      k->v->push_back({a, b});
    }
    k->n++;
    return GSHGEOM_OK;
  };
  ASSERT_EQ(gshgeom_stroke_dash_each(d.data(), 2, 0, 1000000, GSHGEOM_STROKE_DASH_STEPS_DEFAULT, cb, &c), GSHGEOM_OK);
  EXPECT_EQ(c.n, 1000000u);
  for (size_t i = 0; i < out.size(); i++) {
    EXPECT_EQ(out[i].a, (double)(i * 100000));
    EXPECT_EQ(out[i].b, (double)(i * 100000) + 0.5);
  }
}

TEST(DashPositions, ErrorsAndLimits) {
  std::vector<Pair> out;
  double d[2] = {5, 3};
  EXPECT_EQ(gshgeom_stroke_dash_each(d, 2, 0, 1, 10, nullptr, &out), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_stroke_dash_each(nullptr, 2, 0, 1, 10, collect, &out), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_stroke_dash_each(d, 2, NAN, 1, 10, collect, &out), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_stroke_dash_each(d, 2, 0, -1, 10, collect, &out), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_stroke_dash_each(d, 2, 0, INFINITY, 10, collect, &out), GSHGEOM_ERR_INVALID);
  double bad[2] = {INFINITY, 1};
  EXPECT_EQ(gshgeom_stroke_dash_each(bad, 2, 0, 1, 10, collect, &out), GSHGEOM_ERR_INVALID);
  /* SG2's bounds */
  double tiny[2] = {std::nextafter(0x1p-64, 0.0), 1};
  EXPECT_EQ(gshgeom_stroke_dash_each(tiny, 2, 0, 1, 10, collect, &out), GSHGEOM_ERR_LIMIT);
  double huge[2] = {std::nextafter(0x1p64, 1e300), 1};
  EXPECT_EQ(gshgeom_stroke_dash_each(huge, 2, 0, 1, 10, collect, &out), GSHGEOM_ERR_LIMIT);
  double ok[2] = {0x1p-64, 0x1p64};
  out.clear();
  EXPECT_EQ(gshgeom_stroke_dash_each(ok, 2, 0, 1, 10, collect, &out), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_stroke_dash_each(d, 2, std::nextafter(0x1p64, 1e300), 1, 10, collect, &out), GSHGEOM_ERR_LIMIT);
  /* too many entries */
  out.clear();
  EXPECT_EQ(gshgeom_stroke_dash_each(d, 2, 0, 1000, 5, collect, &out), GSHGEOM_ERR_LIMIT);
  /* more than the dash count can be */
  std::vector<double> many(65536, 1.0);
  EXPECT_EQ(gshgeom_stroke_dash_each(many.data(), many.size() / 2 * 2 + 1 - 1, 0, 1, 10, collect, &out), GSHGEOM_ERR_LIMIT);
  /* the callback's first error stops it and is returned */
  struct Stop { int n; } st = {0};
  auto cb = [](void * ctx, double, double) -> GSHGEOM_Result {
    return ++static_cast<Stop *>(ctx)->n == 2 ? GSHGEOM_ERR_IO : GSHGEOM_OK;
  };
  EXPECT_EQ(gshgeom_stroke_dash_each(d, 2, 0, 100, 1000, cb, &st), GSHGEOM_ERR_IO);
  EXPECT_EQ(st.n, 2);
}

TEST(DashDxf, EveryVector) {
  EXPECT_GT(dash_vectors::dxf_count, (size_t)10);
  for (size_t i = 0; i < dash_vectors::dxf_count; i++) {
    const dash_vectors::Dxf & v = dash_vectors::dxf[i];
    std::vector<double> out(v.count + 2);
    size_t n = 99;
    double off = -1;
    GSHGEOM_Result r = gshgeom_stroke_dashes_from_dxf(v.elements, v.count, v.scale, out.data(), out.size(), &n, &off);
    if (v.kind == 2) {
      EXPECT_EQ(r, GSHGEOM_ERR_UNSUPPORTED) << v.name;
      continue;
    }
    ASSERT_EQ(r, GSHGEOM_OK) << v.name;
    if (v.kind == 1) {
      EXPECT_EQ(n, 0u) << v.name;
      continue;
    }
    ASSERT_EQ(n, v.dash_count) << v.name;
    for (size_t k = 0; k < n; k++) {
      EXPECT_EQ(out[k], v.dashes[k]) << v.name << " entry " << k;
    }
    EXPECT_EQ(off, v.offset) << v.name;
  }
}

TEST(DashDxf, TheKnownPatterns) {
  /* the sign convention (dash positive, space negative, dot zero) and a few standard linetypes of acad.lin */
  double dotted[2] = {0.0, -0.2};
  double out[8];
  size_t n;
  double off;
  ASSERT_EQ(gshgeom_stroke_dashes_from_dxf(dotted, 2, 1.0, out, 8, &n, &off), GSHGEOM_OK);
  ASSERT_EQ(n, 2u);
  EXPECT_EQ(out[0], 0.0);
  EXPECT_EQ(out[1], 0.2);
  EXPECT_EQ(off, 0.0);
  /* DASHED: .5, -.25, scaled by 2 */
  double dashed[2] = {0.5, -0.25};
  ASSERT_EQ(gshgeom_stroke_dashes_from_dxf(dashed, 2, 2.0, out, 8, &n, &off), GSHGEOM_OK);
  ASSERT_EQ(n, 2u);
  EXPECT_EQ(out[0], 1.0);
  EXPECT_EQ(out[1], 0.5);
}

TEST(DashDxf, Errors) {
  double out[8], el[2] = {1, -1};
  size_t n;
  double off;
  EXPECT_EQ(gshgeom_stroke_dashes_from_dxf(nullptr, 2, 1, out, 8, &n, &off), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_stroke_dashes_from_dxf(el, 0, 1, out, 8, &n, &off), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_stroke_dashes_from_dxf(el, 2, 0, out, 8, &n, &off), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_stroke_dashes_from_dxf(el, 2, -1, out, 8, &n, &off), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_stroke_dashes_from_dxf(el, 2, NAN, out, 8, &n, &off), GSHGEOM_ERR_INVALID);
  double nan_el[2] = {NAN, -1};
  EXPECT_EQ(gshgeom_stroke_dashes_from_dxf(nan_el, 2, 1, out, 8, &n, &off), GSHGEOM_ERR_INVALID);
  double spaces[2] = {-1, -2};
  EXPECT_EQ(gshgeom_stroke_dashes_from_dxf(spaces, 2, 1, out, 8, &n, &off), GSHGEOM_ERR_INVALID); /* draws nothing */
  double dots[2] = {0, 0};
  EXPECT_EQ(gshgeom_stroke_dashes_from_dxf(dots, 2, 1, out, 8, &n, &off), GSHGEOM_ERR_INVALID);   /* no length */
  /* the capacity: three runs need four slots */
  double three[3] = {1, -1, 1};
  EXPECT_EQ(gshgeom_stroke_dashes_from_dxf(three, 3, 1, out, 2, &n, &off), GSHGEOM_ERR_LIMIT);
  /* SG2's bounds, after scaling */
  double big[2] = {0x1p63, -1};
  EXPECT_EQ(gshgeom_stroke_dashes_from_dxf(big, 2, 4.0, out, 8, &n, &off), GSHGEOM_ERR_LIMIT);
  double small[2] = {0x1p-60, -1};
  EXPECT_EQ(gshgeom_stroke_dashes_from_dxf(small, 2, 0x1p-10, out, 8, &n, &off), GSHGEOM_ERR_LIMIT);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
