/**
 * @file
 *
 * Two strokes, found by a differential search against a variant of the library (random polylines and curves on a small grid, every
 * join and cap, under the identity and under a shear), whose outline differs from the variant's: one with the miter-clip join of a
 * closed contour whose tip is reached through an error bound, one with an arcs join beside a cubic that doubles back.  The outline
 * hash is that of the library as it stands, taken before the variant was built, so that each test fails on the variant.
 *
 * Copyright 2026 by Corey Pennycuff
 */

#include "../test_helpers.h"

#include <cstdint>
#include <cstring>

namespace {

GSHGEOM_Point P(double x, double y) {
  return GSHGEOM_Point{x, y, 0.0};
}

uint64_t outline_hash(GSHGEOM_Path * p, const GSHGEOM_Stroke & st, const GSHGEOM_Map & m) {
  GSHGEOM_Stroke_Outline * o = nullptr;
  EXPECT_EQ(gshgeom_stroke_expand(p, &st, &m, nullptr, nullptr, &o), GSHGEOM_OK);
  const uint64_t h = o ? gshgeom_stroke_outline_hash(o) : 0;
  gshgeom_stroke_outline_destroy(o);
  return h;
}

}  // namespace

TEST(StrokeRegress, AClosedContourWithAMiterClipJoinAndACubicThatTurnsBack) {
  GSHGEOM_Path * p = nullptr;
  ASSERT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &p), GSHGEOM_OK);
  gshgeom_path_move_to(p, P(0x1p+1, 0));
  gshgeom_path_line_to(p, P(0x1.8p+2, 0));
  gshgeom_path_cubic_to(p, P(0x1p+3, 0), P(0x1p+3, 0), P(0x1.8p+2, 0x1p+0));
  gshgeom_path_quadratic_to(p, P(0x1.cp+2, 0x1p+0), P(0x1.cp+2, 0x1p+0));
  gshgeom_path_close(p);
  gshgeom_path_move_to(p, P(0x1.8p+2, 0x1p+0));
  gshgeom_path_line_to(p, P(0, 0));
  gshgeom_path_line_to(p, P(0, 0));
  const double w[3] = {1.0, 0x1.767c64401397cp+0, 1.0};
  gshgeom_path_rational_quadratic_to(p, P(0, 0x1p+0), P(0x1p+3, 0x1p+0), w);
  gshgeom_path_close(p);
  GSHGEOM_Stroke st;
  gshgeom_stroke_init(&st);
  st.width = 4.0;
  st.width_space = GSHGEOM_STROKE_WIDTH_LOCAL;
  st.join = GSHGEOM_STROKE_JOIN_MITER_CLIP;
  st.miter_limit = 3.0;
  st.cap = GSHGEOM_STROKE_CAP_BUTT;
  st.degenerate = GSHGEOM_STROKE_DEGENERATE_OMIT;
  st.alignment = GSHGEOM_STROKE_ALIGN_OUTSIDE;
  GSHGEOM_Map m;
  ASSERT_EQ(gshgeom_map_identity(2, &m), GSHGEOM_OK);
  EXPECT_EQ(outline_hash(p, st, m), UINT64_C(0x0841d210f837ce2b));
  gshgeom_path_destroy(p);
}

TEST(StrokeRegress, AnArcsJoinBesideACubicThatDoublesBackUnderAShear) {
  GSHGEOM_Path * p = nullptr;
  ASSERT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &p), GSHGEOM_OK);
  gshgeom_path_move_to(p, P(0, 1));
  gshgeom_path_line_to(p, P(8, 2));
  gshgeom_path_line_to(p, P(2, 1));
  gshgeom_path_cubic_to(p, P(8, 2), P(4, 1), P(8, 0));
  gshgeom_path_line_to(p, P(8, 0));
  gshgeom_path_line_to(p, P(0, 2));
  GSHGEOM_Stroke st;
  gshgeom_stroke_init(&st);
  st.width = 1.5;
  st.width_space = GSHGEOM_STROKE_WIDTH_DEVICE;
  st.join = GSHGEOM_STROKE_JOIN_ARCS;
  st.miter_limit = 11.0;
  st.cap = GSHGEOM_STROKE_CAP_ROUND;
  st.degenerate = GSHGEOM_STROKE_DEGENERATE_DRAW;
  GSHGEOM_Map m;
  ASSERT_EQ(gshgeom_map_identity(2, &m), GSHGEOM_OK);
  m.m[0][1] = 0.3;
  EXPECT_EQ(outline_hash(p, st, m), UINT64_C(0x16470a3a518282f7));
  gshgeom_path_destroy(p);
}

TEST(StrokeRegress, TheClipOfAnOffCenterStrokeFlattensConicsAtTheTolerance) {
  // a closed contour of four rational quadratics (a circle of radius 20) stroked inside and outside: the clip region F is the centerline
  // converted to the device at the tolerance asked for, so its points (and the outline's hash) follow the tolerance; the hashes are
  // those of the library as it stands, for the default tolerance and for 2, 8 and 64
  GSHGEOM_Path * p = nullptr;
  ASSERT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &p), GSHGEOM_OK);
  gshgeom_path_move_to(p, P(20, 0));
  const double w[3] = {1.0, 0.7071067811865476, 1.0};
  gshgeom_path_rational_quadratic_to(p, P(20, 20), P(0, 20), w);
  gshgeom_path_rational_quadratic_to(p, P(-20, 20), P(-20, 0), w);
  gshgeom_path_rational_quadratic_to(p, P(-20, -20), P(0, -20), w);
  gshgeom_path_rational_quadratic_to(p, P(20, -20), P(20, 0), w);
  gshgeom_path_close(p);
  GSHGEOM_Map m;
  ASSERT_EQ(gshgeom_map_identity(2, &m), GSHGEOM_OK);
  const int tols[4] = {0, 2, 8, 64};
  const uint64_t want[2][4] = {{UINT64_C(0xfe5b8d1309ac86b8), UINT64_C(0xbb9690f4df9129b1), UINT64_C(0x605341d33afae7ab), UINT64_C(0xda93f5ada516b717)},
      {UINT64_C(0x29ce73b75d95874e), UINT64_C(0xa754e57d942dabe3), UINT64_C(0x8e166163ed2a6c4d), UINT64_C(0xa7cdde8bc1d1a815)}};
  const size_t points[4] = {105, 129, 65, 17};
  for (int al = 1; al <= 2; al++) {
    for (int t = 0; t < 4; t++) {
      GSHGEOM_Stroke st;
      gshgeom_stroke_init(&st);
      st.width = 2.0;
      st.alignment = (uint32_t)al;
      GSHGEOM_Stroke_Outline * o = nullptr;
      const GSHGEOM_Result r = tols[t] == 0 ? gshgeom_stroke_expand(p, &st, &m, nullptr, nullptr, &o)
                                            : gshgeom_stroke_expand_with_tolerance(p, &st, &m, tols[t], nullptr, nullptr, &o);
      ASSERT_EQ(r, GSHGEOM_OK);
      EXPECT_EQ(gshgeom_stroke_outline_hash(o), want[al - 1][t]) << al << " " << tols[t];
      EXPECT_EQ(gshgeom_flat_path_point_count(gshgeom_stroke_outline_clip(o)), points[t]) << al << " " << tols[t];
      gshgeom_stroke_outline_destroy(o);
    }
  }
  gshgeom_path_destroy(p);
}

TEST(StrokeRegress, ATurnTheMapMakesInvisibleIsNotAJoin) {
  // a turn of 1e-300 at the second vertex, under a map that shrinks by 1e-30: the mapped points are exactly collinear and in line
  // (sigma 0, the turn's dot positive), so no join is made; a join taken for a reversal would add a round cap's worth of outline
  GSHGEOM_Path * p = nullptr;
  ASSERT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &p), GSHGEOM_OK);
  gshgeom_path_move_to(p, P(0, 0));
  gshgeom_path_line_to(p, P(0x1.8p+1, 0));
  gshgeom_path_line_to(p, P(0x1p+2, -0x1.56e1fc2f8f359p-997));
  GSHGEOM_Stroke st;
  gshgeom_stroke_init(&st);
  st.width = 1.5;
  st.width_space = GSHGEOM_STROKE_WIDTH_DEVICE;
  st.join = GSHGEOM_STROKE_JOIN_MITER_CLIP;
  st.miter_limit = 9.0;
  st.cap = GSHGEOM_STROKE_CAP_ROUND;
  st.degenerate = GSHGEOM_STROKE_DEGENERATE_DRAW;
  GSHGEOM_Map m;
  ASSERT_EQ(gshgeom_map_scale(2, 1e-30, 1e-30, 1.0, &m), GSHGEOM_OK);
  EXPECT_EQ(outline_hash(p, st, m), UINT64_C(0xa788bfaf0b61f577));
  gshgeom_path_destroy(p);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
