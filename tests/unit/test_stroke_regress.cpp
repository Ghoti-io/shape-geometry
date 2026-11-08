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

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
