/**
 * @file
 *
 * The flattened path built from explicit points (flatten.h, ::gshgeom_flat_path_from_points): the result is the
 * arguments and nothing else (points, subpaths, flags, hash), every way the subpaths can fail to tile the points
 * is refused before anything is allocated, the point cap is exact, and an allocation that fails at any call leaves
 * nothing behind.
 *
 * Copyright 2026 by Corey Pennycuff
 */

#include "test_helpers.h"

#include "../alloc_tracker.h"

#include <vector>

namespace {

using Dp = GSHGEOM_Device_Point;
using Sp = GSHGEOM_Flat_Subpath;

const std::vector<Dp> kPts = {{0, 0}, {64, 0}, {64, 64}, {0, 64}, {200, 5}, {-7, 9}, {300, 300}};
const std::vector<Sp> kSubs = {{0, 4, 1}, {4, 2, 0}, {6, 1, 1}};

GSHGEOM_Result build(const std::vector<Dp> & p, const std::vector<Sp> & s, GSHGEOM_Flat_Path ** out,
    const GSHGEOM_Allocator * a = nullptr, const GSHGEOM_Limits * l = nullptr) {
  return gshgeom_flat_path_from_points(p.empty() ? nullptr : p.data(), p.size(), s.empty() ? nullptr : s.data(),
      s.size(), a, l, out);
}

TEST(FlatPathFromPoints, TheResultIsTheArgumentsPointsSubpathsAndFlags) {
  GSHGEOM_Flat_Path * fp = nullptr;
  ASSERT_EQ(build(kPts, kSubs, &fp), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_flat_path_point_count(fp), kPts.size());
  ASSERT_EQ(gshgeom_flat_path_subpath_count(fp), kSubs.size());
  const Dp * got = gshgeom_flat_path_points(fp);
  for (size_t i = 0; i < kPts.size(); i++) {
    EXPECT_EQ(got[i].x, kPts[i].x);
    EXPECT_EQ(got[i].y, kPts[i].y);
  }
  for (size_t s = 0; s < kSubs.size(); s++) {
    Sp sub;
    ASSERT_EQ(gshgeom_flat_path_subpath(fp, s, &sub), GSHGEOM_OK);
    EXPECT_EQ(sub.first, kSubs[s].first);
    EXPECT_EQ(sub.count, kSubs[s].count);
    EXPECT_EQ(sub.closed, kSubs[s].closed);
  }
  gshgeom_flat_path_destroy(fp);
}

TEST(FlatPathFromPoints, ItHashesAsTheSamePathMadeByTheFlattener) {
  /* A unit square closed and a line left open, flattened from a model path. */
  GSHGEOM_Path * mp = nullptr;
  ASSERT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &mp), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_move_to(mp, GSHGEOM_Point{0, 0, 0}), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(mp, GSHGEOM_Point{1, 0, 0}), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(mp, GSHGEOM_Point{1, 1, 0}), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_close(mp), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_move_to(mp, GSHGEOM_Point{3, 3, 0}), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(mp, GSHGEOM_Point{4, 5, 0}), GSHGEOM_OK);
  GSHGEOM_Map map;
  ASSERT_EQ(gshgeom_map_identity(2, &map), GSHGEOM_OK);
  GSHGEOM_Flat_Path * flat = nullptr;
  ASSERT_EQ(gshgeom_path_flatten(mp, &map, nullptr, nullptr, &flat), GSHGEOM_OK);
  gshgeom_path_destroy(mp);
  std::vector<Dp> pts(gshgeom_flat_path_points(flat), gshgeom_flat_path_points(flat) + gshgeom_flat_path_point_count(flat));
  std::vector<Sp> subs;
  for (size_t s = 0; s < gshgeom_flat_path_subpath_count(flat); s++) {
    Sp sub;
    ASSERT_EQ(gshgeom_flat_path_subpath(flat, s, &sub), GSHGEOM_OK);
    subs.push_back(sub);
  }
  GSHGEOM_Flat_Path * built = nullptr;
  ASSERT_EQ(build(pts, subs, &built), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_flat_path_hash(built), gshgeom_flat_path_hash(flat));
  gshgeom_flat_path_destroy(built);
  gshgeom_flat_path_destroy(flat);
}

TEST(FlatPathFromPoints, NoPointsAndNoSubpathsMakeAnEmptyPath) {
  GSHGEOM_Flat_Path * fp = nullptr;
  ASSERT_EQ(gshgeom_flat_path_from_points(nullptr, 0, nullptr, 0, nullptr, nullptr, &fp), GSHGEOM_OK);
  ASSERT_NE(fp, nullptr);
  EXPECT_EQ(gshgeom_flat_path_point_count(fp), 0u);
  EXPECT_EQ(gshgeom_flat_path_subpath_count(fp), 0u);
  EXPECT_EQ(gshgeom_flat_path_points(fp), nullptr);
  gshgeom_flat_path_destroy(fp);
}

TEST(FlatPathFromPoints, EveryWayTheSubpathsFailToTileThePointsIsRefusedAndNothingIsWritten) {
  struct Case {
    const char * name;
    std::vector<Sp> subs;
  };
  const std::vector<Case> cases = {
      {"starts past 0", {{1, 6, 0}}},
      {"a gap", {{0, 3, 0}, {4, 3, 0}}},
      {"an overlap", {{0, 4, 0}, {3, 4, 0}}},
      {"an empty subpath", {{0, 0, 0}, {0, 7, 0}}},
      {"covers too few", {{0, 6, 0}}},
      {"covers too many", {{0, 8, 0}}},
      {"a count that wraps", {{0, (size_t)-3, 0}}},
      {"counts that wrap to the total", {{0, (size_t)-3, 0}, {(size_t)-3, 10, 0}}},
      {"out of order but covering", {{4, 3, 0}, {0, 4, 0}}},
      {"overlapping but summing right", {{0, 4, 0}, {2, 3, 0}}},
      {"a flag of 2", {{0, 7, 2}}},
      {"a flag of -1", {{0, 7, (int)-1}}},
      {"no subpaths for points", {}},
  };
  for (const Case & c : cases) {
    GSHGEOM_Flat_Path * fp = reinterpret_cast<GSHGEOM_Flat_Path *>(0x1);
    EXPECT_EQ(build(kPts, c.subs, &fp), GSHGEOM_ERR_INVALID) << c.name;
    EXPECT_EQ(fp, reinterpret_cast<GSHGEOM_Flat_Path *>(0x1)) << c.name << ": *out is written only on success";
  }
  GSHGEOM_Flat_Path * fp = nullptr;
  EXPECT_EQ(build({}, {{0, 1, 0}}, &fp), GSHGEOM_ERR_INVALID) << "subpaths for no points";
  EXPECT_EQ(gshgeom_flat_path_from_points(nullptr, 3, kSubs.data(), 1, nullptr, nullptr, &fp), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_flat_path_from_points(kPts.data(), 7, nullptr, 1, nullptr, nullptr, &fp), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_flat_path_from_points(kPts.data(), 7, kSubs.data(), 3, nullptr, nullptr, nullptr), GSHGEOM_ERR_INVALID);
}

TEST(FlatPathFromPoints, ThePointCapIsExactAndRefusesBeforeAllocating) {
  GSHGEOM_Limits l;
  gshgeom_limits_default(&l);
  l.max_points = kPts.size();
  GSHGEOM_Flat_Path * fp = nullptr;
  EXPECT_EQ(build(kPts, kSubs, &fp, nullptr, &l), GSHGEOM_OK);
  gshgeom_flat_path_destroy(fp);
  l.max_points = kPts.size() - 1;
  fp = nullptr;
  tracker::Tracker t;
  GSHGEOM_Allocator a = tracker::allocator_for(&t);
  EXPECT_EQ(build(kPts, kSubs, &fp, &a, &l), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(fp, nullptr);
  EXPECT_EQ(t.calls, 0) << "the cap is decided before anything is allocated";
}

TEST(FlatPathFromPoints, AnAllocationThatFailsAtAnyCallLeavesNothingBehind) {
  tracker::Tracker count;
  {
    GSHGEOM_Allocator a = tracker::allocator_for(&count);
    GSHGEOM_Flat_Path * fp = nullptr;
    ASSERT_EQ(build(kPts, kSubs, &fp, &a), GSHGEOM_OK);
    gshgeom_flat_path_destroy(fp);
    EXPECT_TRUE(count.live.empty());
  }
  ASSERT_GT(count.calls, 2);
  for (long k = 0; k < count.calls; k++) {
    for (int persistent = 0; persistent < 2; persistent++) {
      tracker::Tracker t;
      t.fail_from = k;
      t.persistent = persistent != 0;
      GSHGEOM_Allocator a = tracker::allocator_for(&t);
      GSHGEOM_Flat_Path * fp = nullptr;
      EXPECT_EQ(build(kPts, kSubs, &fp, &a), GSHGEOM_ERR_OOM) << k;
      EXPECT_EQ(fp, nullptr);
      EXPECT_TRUE(t.live.empty()) << "k=" << k << " leaked";
    }
  }
}

}  // namespace

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
