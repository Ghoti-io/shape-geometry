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
   * dash, which is over (half open, F3): the gap runs 0-3 and the next dash begins at 3 */
  p = positions({5, 3}, -3, 8);
  ASSERT_EQ(p.size(), 1u);
  EXPECT_EQ(p[0].a, 3);
  EXPECT_EQ(p[0].b, 8);
  /* all zero, a negative length, an empty array: solid */
  for (std::vector<double> d : {std::vector<double>{0, 0}, std::vector<double>{1, -1}, std::vector<double>{}}) {
    p = positions(d, 3, 9);
    ASSERT_EQ(p.size(), 1u);
    EXPECT_EQ(p[0].a, 0);
    EXPECT_EQ(p[0].b, 9);
  }
}

TEST(DashPositions, ThePhaseAtTheEndOfADash) {
  /* F3: the CR takes the least index whose prefix sum reaches the phase, which leaves a dash of no length at the start
   * when the phase is exactly the end of a positive dash; every renderer draws nothing there, and so does this: the
   * dash is over, the gap runs 0-3 and the next dash is 3-8 */
  auto p = positions({5, 3}, 5, 20);
  ASSERT_GE(p.size(), 1u);
  EXPECT_EQ(p[0].a, 3);
  EXPECT_EQ(p[0].b, 8);
  /* a dash of no length exactly at the phase is kept: a dotted line starts with a dot */
  p = positions({0, 10}, 0, 5);
  ASSERT_EQ(p.size(), 1u);
  EXPECT_EQ(p[0].a, 0);
  EXPECT_EQ(p[0].b, 0);
  /* the phase at the end of a gap is the start of the next dash */
  p = positions({5, 3}, 8, 5);
  ASSERT_EQ(p.size(), 1u);
  EXPECT_EQ(p[0].a, 0);
  EXPECT_EQ(p[0].b, 5);
  /* a phase used up by the only dash of a short path: nothing */
  p = positions({1, 10}, 1, 5);
  EXPECT_EQ(p.size(), 0u);
  /* inside a gap */
  p = positions({5, 3}, 6, 1);
  EXPECT_EQ(p.size(), 0u);
  /* a length of zero: the first piece only, when the pattern is on at the phase */
  p = positions({5, 3}, 0, 0);
  ASSERT_EQ(p.size(), 1u);
  EXPECT_EQ(p[0].b, 0);
  p = positions({5, 3}, 6, 0);
  EXPECT_EQ(p.size(), 0u);
  p = positions({5, 3}, 5, 0);
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

TEST(DashDxf, ALineThatStartsJustAfterADot) {
  /* a space, a dash, a space and a dot, repeated: the line starts with the space, and the dot comes at the end of each
   * period, never at the start.  The positions of the lowered array must show no dot at 0 and the first dot at the
   * period (within the unit in the last place the lowering moves the start by) */
  const double el[4] = {-2, 3, -1, 0};
  double out[8];
  size_t n = 0;
  double off = 0;
  ASSERT_EQ(gshgeom_stroke_dashes_from_dxf(el, 4, 1.0, out, 8, &n, &off), GSHGEOM_OK);
  std::vector<double> d(out, out + n);
  std::vector<Pair> p = positions(d, off, 14);
  ASSERT_GE(p.size(), 3u);
  EXPECT_GT(p[0].a, 0.0); /* the first thing on the line is the dash, 2 in */
  EXPECT_NEAR(p[0].a, 2.0, 1e-9);
  EXPECT_NEAR(p[0].b, 5.0, 1e-9);
  EXPECT_NEAR(p[1].a, 6.0, 1e-9); /* the dot at the end of the period: 2 + 3 + 1 */
  EXPECT_NEAR(p[1].b, 6.0, 1e-9);
  EXPECT_EQ(p[1].a, p[1].b);
  EXPECT_NEAR(p[2].a, 8.0, 1e-9);
  /* the start is the next double above the dot's position, not the position itself */
  double at = 0;
  for (size_t i = 0; i + 2 < n; i++) {
    at += d[i];
  }
  EXPECT_EQ(off, std::nextafter(at, 1e300));
  /* the dotted line that starts with a dot is not touched: a dot first keeps it */
  const double dotted[2] = {0, -2};
  ASSERT_EQ(gshgeom_stroke_dashes_from_dxf(dotted, 2, 1.0, out, 8, &n, &off), GSHGEOM_OK);
  std::vector<double> dd(out, out + n);
  std::vector<Pair> q = positions(dd, off, 5);
  ASSERT_GE(q.size(), 1u);
  EXPECT_EQ(q[0].a, 0.0);
  EXPECT_EQ(q[0].b, 0.0);
  /* and one with a space first and a lone dot: no dot at 0, a dot at the period */
  const double lone[2] = {-1, 0};
  ASSERT_EQ(gshgeom_stroke_dashes_from_dxf(lone, 2, 1.0, out, 8, &n, &off), GSHGEOM_OK);
  std::vector<double> ld(out, out + n);
  std::vector<Pair> l = positions(ld, off, 3.5);
  ASSERT_EQ(l.size(), 3u);
  EXPECT_NEAR(l[0].a, 1.0, 1e-9);
  EXPECT_EQ(l[0].a, l[0].b);
  EXPECT_NEAR(l[2].a, 3.0, 1e-9);
  /* a pattern too small for the offset band to place the start after its lone dot */
  const double tiny[2] = {-0x1p-64, 0};
  EXPECT_EQ(gshgeom_stroke_dashes_from_dxf(tiny, 2, 1.0, out, 8, &n, &off), GSHGEOM_ERR_LIMIT);
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

namespace {

/** The winding number of an outline around (x, y) pixels. */
int winding_at(const GSHGEOM_Flat_Path * fp, double x, double y) {
  const GSHGEOM_Device_Point * pts = gshgeom_flat_path_points(fp);
  double px = x * 64.0, py = y * 64.0;
  int w = 0;
  for (size_t s = 0; s < gshgeom_flat_path_subpath_count(fp); s++) {
    GSHGEOM_Flat_Subpath sp;
    gshgeom_flat_path_subpath(fp, s, &sp);
    for (size_t i = 0; i < sp.count; i++) {
      GSHGEOM_Device_Point a = pts[sp.first + i], b = pts[sp.first + (i + 1) % sp.count];
      double c = (double)(b.x - a.x) * (py - a.y) - (px - a.x) * (double)(b.y - a.y);
      if (a.y <= py && b.y > py && c > 0) {
        w++;
      } else if (a.y > py && b.y <= py && c < 0) {
        w--;
      }
    }
  }
  return w;
}

/** A closed rectangle w by h from the origin, stroked at width 2 with miter joins and butt caps, dashed; whether the
 * outline covers the point (x, y). */
bool covered(double w, double h, std::vector<double> d, double offset, double x, double y) {
  GSHGEOM_Path * path = nullptr;
  EXPECT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &path), GSHGEOM_OK);
  GSHGEOM_Point p0 = {0, 0, 0}, p1 = {w, 0, 0}, p2 = {w, h, 0}, p3 = {0, h, 0};
  gshgeom_path_move_to(path, p0);
  gshgeom_path_line_to(path, p1);
  gshgeom_path_line_to(path, p2);
  gshgeom_path_line_to(path, p3);
  gshgeom_path_close(path);
  GSHGEOM_Stroke s;
  gshgeom_stroke_init(&s);
  s.width = 2;
  s.dashes = d.data();
  s.dash_count = d.size();
  s.dash_offset = offset;
  GSHGEOM_Map m;
  gshgeom_map_identity(2, &m);
  GSHGEOM_Stroke_Outline * o = nullptr;
  EXPECT_EQ(gshgeom_stroke_expand(path, &s, &m, nullptr, nullptr, &o), GSHGEOM_OK);
  bool in = winding_at(gshgeom_stroke_outline_body(o), x, y) != 0;
  gshgeom_stroke_outline_destroy(o);
  gshgeom_path_destroy(path);
  return in;
}

} // namespace

TEST(DashWrap, ADashAcrossTheStartIsOneDash) {
  /* a 20 by 20 square (perimeter 80), the start vertex (0, 0) joined by a miter whose outer corner is (-1, -1); the
   * point (-0.7, -0.7) is covered only by that join */
  /* [50, 10] at the phase 10: the first piece began before the path and the last dash (50 to 100) would go on past
   * the end: one dash, joined at the start vertex */
  EXPECT_TRUE(covered(20, 20, {50, 10}, 10, -0.7, -0.7));
  /* the same pattern with the phase 0: the first dash begins exactly at the start, so the two stay separate and the
   * start vertex has two butt ends and no join */
  EXPECT_FALSE(covered(20, 20, {50, 10}, 0, -0.7, -0.7));
  /* a dash that ends exactly at the start (20 by 10, perimeter 60, [35, 5] from the phase 15: dashes 0-20, 25-60) and
   * a first piece that began before it: two dashes */
  EXPECT_FALSE(covered(20, 10, {35, 5}, 15, -0.7, -0.7));
  /* a dash that covers the whole subpath is the closed stroke: joined at the start vertex, no caps */
  EXPECT_TRUE(covered(20, 20, {1000, 5}, 10, -0.7, -0.7));
  /* and a pattern that is on only inside, with gaps at the start, covers nothing there */
  EXPECT_FALSE(covered(20, 20, {10, 70}, 0, -0.7, -0.7));
  /* the body is still covered where the dash is, and the gap is a gap: (10, 0) is position 10 on the top edge, in the
   * first dash 0-40; (15, 20) is position 45 on the bottom edge, in the gap 40-50 */
  EXPECT_TRUE(covered(20, 20, {50, 10}, 10, 10, 0));
  EXPECT_FALSE(covered(20, 20, {50, 10}, 10, 15, 20));
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
