/**
 * @file
 *
 * The work of the lattice's search and of the clip, counted through the internal seams (lattice_internal.h, clip_internal.h): the
 * guards that only save work (an interval that is already empty is not searched further, a side of zero length is not a pass over
 * the points) are told from their absence by the count, which nothing else reports.  The counts are those of the code as it stands:
 * a guard taken out raises one of them.
 *
 * Copyright 2026 by Corey Pennycuff
 */

#include "../test_helpers.h"

#include <ghoti.io/shape-geometry/clip.h>
#include <ghoti.io/shape-geometry/flatten.h>
#include <ghoti.io/shape-geometry/lattice.h>

extern "C" {
#include "../../src/clip/clip_internal.h"
#include "../../src/lattice/lattice_internal.h"
#include "../../src/path/path_internal.h"
}

#include <cmath>
#include <utility>
#include <cstdint>
#include <cstring>
#include <vector>

namespace {

uint64_t work(const GSHGEOM_Lattice & l, const GSHGEOM_Map & m, const int64_t rect[4], uint64_t * count = nullptr,
    GSHGEOM_Result * result = nullptr) {
  GSHGEOM_Lattice_Stats st = {12345};
  uint64_t n = 99;
  const GSHGEOM_Result r = gshgeom_lattice_each_counted(&m, &l, rect, 100000, nullptr, nullptr, &n, &st);
  if (count) {
    *count = n;
  }
  if (result) {
    *result = r;
  }
  return st.evaluations;
}

GSHGEOM_Map identity() {
  GSHGEOM_Map m;
  EXPECT_EQ(gshgeom_map_identity(2, &m), GSHGEOM_OK);
  return m;
}

}  // namespace

TEST(WorkCounts, ARectangleInAGapBetweenCellsEndsTheSearchAtTheFirstEmptyInterval) {
  const GSHGEOM_Map m = identity();
  // cells 1 wide every 1000: the rectangle [500, 600] holds none, which the first two constraints of the row already show
  GSHGEOM_Lattice l = {0.0, 0.0, 1.0, 1.0, 1000.0, 1.0, 1u, 0u};
  const int64_t gap[4] = {500, 0, 600, 1};
  uint64_t count = 7;
  GSHGEOM_Result r = GSHGEOM_OK;
  const uint64_t w = work(l, m, gap, &count, &r);
  EXPECT_EQ(r, GSHGEOM_OK);
  EXPECT_EQ(count, 0u);
  EXPECT_EQ(w, 4u);
  // the same lattice with a rectangle that holds a cell: every constraint is searched
  const int64_t hit[4] = {1000, 0, 1001, 1};
  EXPECT_EQ(work(l, m, hit, &count, &r), 26u);
  EXPECT_EQ(count, 1u);
}

TEST(WorkCounts, ACellSmallerThanTheSmallestMapEntryHasCornersThatCoincideAndTheirDirectionsAreSkipped) {
  GSHGEOM_Map m = identity();
  m.m[0][0] = 0x1p-600;
  m.m[1][1] = 0x1p-600;
  // the cell 2^-600 square under 2^-600: its four corners are the same point (the product underflows), so every pair gives the
  // direction zero
  GSHGEOM_Lattice l = {0.0, 0.0, 0x1p-600, 0x1p-600, 1.0, 1.0, 1u, 1u};
  const int64_t rect[4] = {0, 0, 0, 0};
  uint64_t count = 0;
  GSHGEOM_Result r = GSHGEOM_OK;
  const uint64_t w = work(l, m, rect, &count, &r);
  EXPECT_EQ(r, GSHGEOM_OK);
  EXPECT_EQ(w, 10u);
}

TEST(WorkCounts, TheSameInputsGiveTheSameWorkWhateverIsAskedFirst) {
  // a fixed set: the counts are a function of the input, and a plain call and a counted call report the same rows
  const GSHGEOM_Map m = identity();
  const struct { GSHGEOM_Lattice l; int64_t rect[4]; uint64_t work; } cases[] = {
      {{0.0, 0.0, 1.0, 1.0, 1000.0, 1.0, 1u, 0u}, {500, 0, 600, 1}, 4},
      {{0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1u, 1u}, {-3, -3, 3, 3}, 212},
      {{0.5, 0.25, 2.0, 1.0, 3.0, 2.0, 1u, 1u}, {-40, -40, 40, 40}, 1044},
      {{0.0, 0.0, 1.0, 1.0, 1e9, 1.0, 1u, 0u}, {5000, 0, 6000, 1}, 4},
      {{0.0, 0.0, 0x1p-30, 0x1p-30, 0x1p-30, 0x1p-30, 1u, 1u}, {0, 0, 0, 0}, 56},
      {{0.0, 0.0, 0x1p-600, 0x1p-600, 0x1p-600, 0x1p-600, 1u, 1u}, {0, 0, 0, 0}, 78},
      {{0.0, 0.0, 0x1p-600, 0x1p-600, 0x1p-600, 0x1p-600, 1u, 1u}, {5, 5, 6, 6}, 2},
      {{-0x1p-600, 0.0, 0x1p-600, 1.0, 0x1p-600, 1.0, 1u, 0u}, {0, 0, 0, 0}, 26},
      {{0.0, 0.0, 1.0, 1.0, 0x1p-20, 0x1p-20, 1u, 1u}, {-4, -4, 4, 4}, 30},
      {{0.0, 0.0, 3.0, 3.0, 7.0, 7.0, 0u, 1u}, {-100, -100, 100, 100}, 26},
  };
  for (const auto & c : cases) {
    EXPECT_EQ(work(c.l, m, c.rect), c.work);
  }
}

TEST(WorkCounts, AQuadrilateralWithARepeatedCornerCutsOncePerSideThatCuts) {
  const GSHGEOM_Device_Point pts[3] = {{-100, -100}, {300, 0}, {0, 300}};
  const GSHGEOM_Flat_Subpath sub = {0, 3, 1};
  GSHGEOM_Flat_Path * path = nullptr;
  ASSERT_EQ(gshgeom_flat_path_from_points(pts, 3, &sub, 1, nullptr, nullptr, &path), GSHGEOM_OK);
  const GSHGEOM_Device_Point four[4] = {{0, 0}, {200, 0}, {200, 200}, {0, 200}};
  const GSHGEOM_Device_Point three[4] = {{0, 0}, {0, 0}, {200, 0}, {100, 200}};
  GSHGEOM_Flat_Path * out = nullptr;
  GSHGEOM_Clip_Stats st = {0, 0, 0, 77};
  ASSERT_EQ(gshgeom_clip_convex_counted(path, four, nullptr, nullptr, &out, &st), GSHGEOM_OK);
  EXPECT_EQ(st.cuts, 4u);
  gshgeom_flat_path_destroy(out);
  st.cuts = 77;
  ASSERT_EQ(gshgeom_clip_convex_counted(path, three, nullptr, nullptr, &out, &st), GSHGEOM_OK);
  EXPECT_EQ(st.cuts, 3u);  // the side from the corner to itself is not a pass over the points
  gshgeom_flat_path_destroy(out);
  gshgeom_flat_path_destroy(path);
}

TEST(WorkCounts, AnArcThatCannotBeDecomposedIsRefusedBeforeAnyTrigonometry) {
  GSHGEOM_Arc a;
  std::memset(&a, 0, sizeof a);
  a.u = GSHGEOM_Point{1, 0, 0};
  a.v = GSHGEOM_Point{0, 1, 0};
  a.unit = GSHGEOM_ANGLE_DEGREES;
  a.sweep = 90;
  GSHGEOM_Curve out[GSHGEOM_ARC_MAX_PIECES];
  int np = 5;
  GSHGEOM_Arc_Stats st = {99};
  auto trig = [&](const GSHGEOM_Arc & arc, GSHGEOM_Curve * o, int * n, GSHGEOM_Result * r) {
    st.trig = 99;
    *r = gshgeom_arc_pieces_counted(&arc, o, n, &st);
    return st.trig;
  };
  GSHGEOM_Result r = GSHGEOM_OK;
  EXPECT_EQ(trig(a, nullptr, &np, &r), 0u);
  EXPECT_EQ(r, GSHGEOM_ERR_INVALID);
  GSHGEOM_Arc b = a;
  b.unit = (GSHGEOM_Angle_Unit)GSHGEOM_ANGLE_UNIT_COUNT;
  EXPECT_EQ(trig(b, out, &np, &r), 0u);
  EXPECT_EQ(r, GSHGEOM_ERR_INVALID);
  b = a;
  b.sweep = std::nan("");
  EXPECT_EQ(trig(b, out, &np, &r), 0u);
  EXPECT_EQ(r, GSHGEOM_ERR_INVALID);
  b.sweep = HUGE_VAL;
  EXPECT_EQ(trig(b, out, &np, &r), 0u);
  b = a;
  b.t0 = std::nan("");
  EXPECT_EQ(trig(b, out, &np, &r), 0u);
  b = a;
  b.sweep = 0;
  EXPECT_EQ(trig(b, out, &np, &r), 0u);
  EXPECT_EQ(trig(a, out, nullptr, &r), 0u);
  // the evaluations of arcs that are decomposed: a quarter turn needs no cosine or tangent (its weight is a constant), the others do
  EXPECT_EQ(trig(a, out, &np, &r), 2u);
  EXPECT_EQ(r, GSHGEOM_OK);
  b = a;
  b.sweep = 45;
  EXPECT_EQ(trig(b, out, &np, &r), 4u);
  b.sweep = 360;
  EXPECT_EQ(trig(b, out, &np, &r), 5u);
}

TEST(WorkCounts, TheTotalWorkOfAFixedSetOfThreeThousandLatticesIsPinned) {
  // lattices of ordinary, fractional, tiny and huge-index shapes from a fixed generator: the sum of the evaluations of the search (every
  // estimate, every step of a gallop and of a bisection) and of the copies reported is one number each, so a change to how the search
  // starts, steps or splits is a change of the sum
  uint64_t seed = 0x9e3779b97f4a7c15ull;
  auto next = [&]() {
    seed ^= seed << 13;
    seed ^= seed >> 7;
    seed ^= seed << 17;
    return seed;
  };
  auto range = [&](int64_t lo, int64_t hi) { return lo + (int64_t)(next() % (uint64_t)(hi - lo + 1)); };
  auto dyadic = [&](int64_t lo, int64_t hi, int shift) { return std::ldexp((double)range(lo, hi), shift); };
  uint64_t total_work = 0, total_count = 0;
  int limits = 0;
  for (int trial = 0; trial < 3000; trial++) {
    const int kind = (int)(next() % 4);
    GSHGEOM_Lattice l;
    std::memset(&l, 0, sizeof l);
    GSHGEOM_Map m = identity();
    int64_t rect[4];
    if (kind == 0) {  // small integers
      l.x = (double)range(-6, 6);
      l.y = (double)range(-6, 6);
      l.w = (double)range(1, 5);
      l.h = (double)range(1, 5);
      l.x_step = (double)range(1, 6);
      l.y_step = (double)range(1, 6);
      for (auto & row : m.m) for (int c = 0; c < 2; c++) row[c] = (double)range(-3, 3);
      m.t[0] = (double)range(-8, 8);
      m.t[1] = (double)range(-8, 8);
      for (int k = 0; k < 4; k++) rect[k] = range(-30, 30);
    } else if (kind == 1) {  // fractions
      l.x = dyadic(-40, 40, -3);
      l.y = dyadic(-40, 40, -3);
      l.w = dyadic(1, 40, -4);
      l.h = dyadic(1, 40, -4);
      l.x_step = dyadic(1, 60, -4);
      l.y_step = dyadic(1, 60, -4);
      m.m[0][0] = dyadic(-24, 24, -3);
      m.m[0][1] = dyadic(-24, 24, -4);
      m.m[1][0] = dyadic(-24, 24, -4);
      m.m[1][1] = dyadic(-24, 24, -3);
      m.t[0] = dyadic(-200, 200, -4);
      m.t[1] = dyadic(-200, 200, -4);
      for (int k = 0; k < 4; k++) rect[k] = range(-300, 300);
    } else if (kind == 2) {  // tiny cells and steps: the coefficients are tiny
      const int e = -(int)range(40, 560);
      l.x = dyadic(-8, 8, e);
      l.y = dyadic(-8, 8, e);
      l.w = dyadic(1, 8, e);
      l.h = dyadic(1, 8, e);
      l.x_step = dyadic(1, 8, e);
      l.y_step = dyadic(1, 8, e);
      m.m[0][0] = (double)range(-2, 2);
      m.m[0][1] = (double)range(-2, 2);
      m.m[1][0] = (double)range(-2, 2);
      m.m[1][1] = (double)range(-2, 2);
      m.t[0] = dyadic(-8, 8, e);
      m.t[1] = dyadic(-8, 8, e);
      for (int k = 0; k < 4; k++) rect[k] = range(-2, 2);
    } else {  // a long way out: indices past a billion, and the range of indices itself
      l.x = (double)range(-3, 3);
      l.y = (double)range(-3, 3);
      l.w = 1.0;
      l.h = 1.0;
      l.x_step = 1.0;
      l.y_step = (double)range(1, 3);
      m.t[0] = (double)range(-1000, 1000);
      const int64_t far = ((int64_t)1 << 53) - range(0, 40);
      rect[0] = (next() & 1) ? far - range(0, 20) : -far;
      rect[2] = (next() & 1) ? far : rect[0] + range(0, 20);
      rect[1] = range(-5, 5);
      rect[3] = rect[1] + range(0, 9);
    }
    l.repeat_x = (uint32_t)(next() % 8 != 0);
    l.repeat_y = (uint32_t)(next() % 3 != 0);
    if (rect[0] > rect[2]) std::swap(rect[0], rect[2]);
    if (rect[1] > rect[3]) std::swap(rect[1], rect[3]);
    GSHGEOM_Lattice_Stats st = {0};
    uint64_t count = 0;
    const GSHGEOM_Result r = gshgeom_lattice_each_counted(&m, &l, rect, 5000, nullptr, nullptr, &count, &st);
    total_work += st.evaluations;
    total_count += count;
    limits += r == GSHGEOM_ERR_LIMIT;
  }
  EXPECT_EQ(total_work, 21752323ull);
  EXPECT_EQ(total_count, 3877875ull);
  EXPECT_EQ(limits, 711);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
