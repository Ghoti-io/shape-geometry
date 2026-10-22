/**
 * @file
 *
 * SG31, the canonical map (strokes story S4): the translation split exactly, the fraction in [0, 1), and the property
 * a cache of expansions rests on: the outline at the canonical map is the same bits at every whole-pixel translation of
 * the map, and moved by the whole pixels it is the stroke of the map itself.
 *
 * Copyright 2026 by Corey Pennycuff
 */

#include "../test_helpers.h"

#include "../stroke_matrix.h"

#include <cmath>
#include <cstdint>
#include <limits>

namespace sm = stroke_matrix;

namespace {

GSHGEOM_Map affine(double a, double b, double c, double d, double tx, double ty) {
  double m[4] = {a, b, c, d};
  double t[2] = {tx, ty};
  GSHGEOM_Map out;
  EXPECT_EQ(gshgeom_map_from_matrix(2, m, t, &out), GSHGEOM_OK);
  return out;
}

uint64_t hash_at(const GSHGEOM_Path * p, const GSHGEOM_Stroke & s, const GSHGEOM_Map & m) {
  GSHGEOM_Stroke_Outline * o = nullptr;
  EXPECT_EQ(gshgeom_stroke_expand(p, &s, &m, nullptr, nullptr, &o), GSHGEOM_OK);
  uint64_t h = gshgeom_stroke_outline_hash(o);
  gshgeom_stroke_outline_destroy(o);
  return h;
}

} // namespace

TEST(StrokeCanonical, TheTranslationIsSplitExactly) {
  struct Case {
    double tx, ty;
    double fx, fy;
    int64_t sx, sy;
  } cases[] = {
      {0.0, 0.0, 0.0, 0.0, 0, 0},
      {3.25, -1.5, 0.25, 0.5, 3, -2},
      {-0.0, 7.0, 0.0, 0.0, 0, 7},
      {-0.125, 0.125, 0.875, 0.125, -1, 0},
      {1099511627775.5, -1099511627775.5, 0.5, 0.5, 1099511627775, -1099511627776},
  };
  for (const Case & c : cases) {
    GSHGEOM_Map m = affine(1.5, 0.25, -0.5, 2.0, c.tx, c.ty);
    GSHGEOM_Map canon;
    int64_t sx = 99;
    int64_t sy = 99;
    ASSERT_EQ(gshgeom_stroke_canonical_map(&m, &canon, &sx, &sy), GSHGEOM_OK) << c.tx;
    EXPECT_EQ(canon.t[0], c.fx) << c.tx;
    EXPECT_EQ(canon.t[1], c.fy) << c.ty;
    EXPECT_FALSE(std::signbit(canon.t[0]));
    EXPECT_FALSE(std::signbit(canon.t[1]));
    EXPECT_EQ(sx, c.sx) << c.tx;
    EXPECT_EQ(sy, c.sy) << c.ty;
    EXPECT_GE(canon.t[0], 0.0);
    EXPECT_LT(canon.t[0], 1.0);
    EXPECT_GE(canon.t[1], 0.0);
    EXPECT_LT(canon.t[1], 1.0);
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        EXPECT_EQ(canon.m[i][j], m.m[i][j]) << "the linear part is untouched";
      }
    }
    EXPECT_EQ(gshgeom_map_check(&canon), GSHGEOM_OK);
  }
  // A tiny negative translation: 1 - 2^-60 rounds to 1, taken as 0 with the whole pixel one more.
  GSHGEOM_Map tiny = affine(1, 0, 0, 1, -0x1p-60, 0.0);
  GSHGEOM_Map canon;
  int64_t sx = 0;
  int64_t sy = 0;
  ASSERT_EQ(gshgeom_stroke_canonical_map(&tiny, &canon, &sx, &sy), GSHGEOM_OK);
  EXPECT_EQ(canon.t[0], 0.0);
  EXPECT_EQ(sx, 0);
}

TEST(StrokeCanonical, RefusesWhatItCannotSplit) {
  GSHGEOM_Map m = affine(1, 0, 0, 1, 0.5, 0.5);
  GSHGEOM_Map canon;
  int64_t sx = 0;
  int64_t sy = 0;
  EXPECT_EQ(gshgeom_stroke_canonical_map(nullptr, &canon, &sx, &sy), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_stroke_canonical_map(&m, nullptr, &sx, &sy), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_stroke_canonical_map(&m, &canon, nullptr, &sy), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_stroke_canonical_map(&m, &canon, &sx, nullptr), GSHGEOM_ERR_INVALID);
  GSHGEOM_Map three;
  ASSERT_EQ(gshgeom_map_identity(3, &three), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_stroke_canonical_map(&three, &canon, &sx, &sy), GSHGEOM_ERR_INVALID);
  GSHGEOM_Map bad = m;
  bad.t[0] = std::numeric_limits<double>::quiet_NaN();
  EXPECT_EQ(gshgeom_stroke_canonical_map(&bad, &canon, &sx, &sy), GSHGEOM_ERR_INVALID);
  GSHGEOM_Map far = affine(1, 0, 0, 1, 0x1p40, 0.0);
  EXPECT_EQ(gshgeom_stroke_canonical_map(&far, &canon, &sx, &sy), GSHGEOM_ERR_LIMIT);
  far = affine(1, 0, 0, 1, 0.0, -0x1p40);
  EXPECT_EQ(gshgeom_stroke_canonical_map(&far, &canon, &sx, &sy), GSHGEOM_ERR_LIMIT);
  GSHGEOM_Map edge = affine(1, 0, 0, 1, 0x1p40 - 1.0, 0.0);
  EXPECT_EQ(gshgeom_stroke_canonical_map(&edge, &canon, &sx, &sy), GSHGEOM_OK);
  EXPECT_EQ(sx, (int64_t)(0x1p40 - 1.0));
}

TEST(StrokeCanonical, TheOutlineIsTheSameBitsAtEveryWholePixelTranslation) {
  GSHGEOM_Path * p = sm::build({{'M', {0, 0}}, {'L', {10, 0}}, {'Q', {14, 0, 14, 6}}, {'L', {6, 9}}});
  GSHGEOM_Stroke s;
  gshgeom_stroke_init(&s);
  s.width = 2.5;
  s.join = GSHGEOM_STROKE_JOIN_ROUND;
  s.cap = GSHGEOM_STROKE_CAP_SQUARE;
  uint64_t first = 0;
  for (double tx : {0.3125, 1.3125, 17.3125, -4.6875, 1000.3125}) {
    for (double ty : {0.8125, -2.1875, 5.8125}) {
      GSHGEOM_Map m = affine(1.25, 0.5, -0.75, 1.5, tx, ty);
      GSHGEOM_Map canon;
      int64_t sx = 0;
      int64_t sy = 0;
      ASSERT_EQ(gshgeom_stroke_canonical_map(&m, &canon, &sx, &sy), GSHGEOM_OK);
      uint64_t h = hash_at(p, s, canon);
      // Dyadic translations: tx + 64 is exact, so the fraction is the same double (0.3 + 64 would round).
      GSHGEOM_Map shifted = affine(1.25, 0.5, -0.75, 1.5, tx + 64.0, ty - 3.0);
      GSHGEOM_Map canon2;
      int64_t sx2 = 0;
      int64_t sy2 = 0;
      ASSERT_EQ(gshgeom_stroke_canonical_map(&shifted, &canon2, &sx2, &sy2), GSHGEOM_OK);
      EXPECT_EQ(canon.t[0], canon2.t[0]);
      EXPECT_EQ(canon.t[1], canon2.t[1]);
      EXPECT_EQ(sx2, sx + 64);
      EXPECT_EQ(sy2, sy - 3);
      EXPECT_EQ(hash_at(p, s, canon2), h) << "tx " << tx << " ty " << ty;
      if (first == 0) {
        first = h;
      }
    }
  }
  EXPECT_NE(first, 0u);
  // A different fraction is a different expansion (26.6 rounds after the map).
  GSHGEOM_Map a = affine(1.25, 0.5, -0.75, 1.5, 0.0, 0.0);
  GSHGEOM_Map b = affine(1.25, 0.5, -0.75, 1.5, 0.3125, 0.0);
  EXPECT_NE(hash_at(p, s, a), hash_at(p, s, b));
  gshgeom_path_destroy(p);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
