/**
 * @file
 *
 * The integer flattener and the flattened path (CAP-3, CAP-13):
 * bit-for-bit agreement with `font`'s flattener (at `font` commit 69d1860) on
 * recorded vectors, the guarantee (every point of the curve within the
 * tolerance less 0.12 units of the polyline) against dense exact evaluation on
 * every size class up to the whole `int32` range, vertex lists pinned by an
 * independent big-integer model, flattening after the map, the structure of
 * the result, the limits, and the hash and the dump pinned by known answers.
 *
 * Copyright 2026 by Corey Pennycuff
 */

#include "test_helpers.h"

#include "../alloc_tracker.h"

#include <chrono>
#include <cstdint>
#include <map>
#include <tuple>
#include <algorithm>
#include <cstdlib>
#include <vector>

namespace {

using Dp = GSHGEOM_Device_Point;

GSHGEOM_Point P(double x, double y) {
  return GSHGEOM_Point{x, y, 0.0};
}

/** A flattened path of the closed polygons @p polys, whose corners are in whole pixels (so in multiples of 64 in 26.6). */
GSHGEOM_Flat_Path * polygons(const std::vector<std::vector<std::pair<int, int>>> & polys) {
  GSHGEOM_Path * mp = nullptr;
  EXPECT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &mp), GSHGEOM_OK);
  for (const auto & poly : polys) {
    for (size_t i = 0; i < poly.size(); i++) {
      GSHGEOM_Point p = P(poly[i].first, poly[i].second);
      EXPECT_EQ(i == 0 ? gshgeom_path_move_to(mp, p) : gshgeom_path_line_to(mp, p), GSHGEOM_OK);
    }
    EXPECT_EQ(gshgeom_path_close(mp), GSHGEOM_OK);
  }
  GSHGEOM_Map map;
  EXPECT_EQ(gshgeom_map_identity(2, &map), GSHGEOM_OK);
  GSHGEOM_Flat_Path * fp = nullptr;
  EXPECT_EQ(gshgeom_path_flatten(mp, &map, nullptr, nullptr, &fp), GSHGEOM_OK);
  gshgeom_path_destroy(mp);
  return fp;
}

using Edge = std::tuple<int32_t, int32_t, int32_t, int32_t>;

/** The directed edges of @p fp, every subpath closed, zero-length ones left out, with their multiplicity. */
std::map<Edge, int> edges_of(const GSHGEOM_Flat_Path * fp) {
  std::map<Edge, int> edges;
  const Dp * pts = gshgeom_flat_path_points(fp);
  for (size_t s = 0; s < gshgeom_flat_path_subpath_count(fp); s++) {
    GSHGEOM_Flat_Subpath sub;
    EXPECT_EQ(gshgeom_flat_path_subpath(fp, s, &sub), GSHGEOM_OK);
    for (size_t k = 0; k < sub.count; k++) {
      Dp a = pts[sub.first + k];
      Dp b = pts[sub.first + (k + 1 == sub.count ? 0 : k + 1)];
      if (a.x != b.x || a.y != b.y) {
        edges[Edge(a.x, a.y, b.x, b.y)]++;
      }
    }
  }
  return edges;
}

/** The winding number of @p fp at (qx, qy) in 26.6 units, a point on no edge's line of a lattice point. */
int winding_at(const GSHGEOM_Flat_Path * fp, double qx, double qy) {
  int wn = 0;
  const Dp * pts = gshgeom_flat_path_points(fp);
  for (size_t s = 0; s < gshgeom_flat_path_subpath_count(fp); s++) {
    GSHGEOM_Flat_Subpath sub;
    EXPECT_EQ(gshgeom_flat_path_subpath(fp, s, &sub), GSHGEOM_OK);
    for (size_t k = 0; k < sub.count; k++) {
      Dp a = pts[sub.first + k];
      Dp b = pts[sub.first + (k + 1 == sub.count ? 0 : k + 1)];
      double left = (double)(b.x - a.x) * (qy - a.y) - (qx - a.x) * (double)(b.y - a.y);
      if (a.y <= qy) {
        if (b.y > qy && left > 0) {
          wn++;
        }
      } else if (b.y <= qy && left < 0) {
        wn--;
      }
    }
  }
  return wn;
}

/** The two paths have the same winding number at a grid of points between the lattice lines and at scattered ones. */
void expect_same_winding(const GSHGEOM_Flat_Path * a, const GSHGEOM_Flat_Path * b, int extent, uint32_t seed) {
  for (int i = 0; i < extent; i++) {
    for (int j = 0; j < extent; j++) {
      double qx = (i + 0.37) * 64.0;
      double qy = (j + 0.61) * 64.0;
      ASSERT_EQ(winding_at(a, qx, qy), winding_at(b, qx, qy)) << i << "," << j;
    }
  }
  uint32_t r = seed;
  for (int n = 0; n < 400; n++) {
    r = r * 1664525u + 1013904223u;
    double qx = (r >> 8) % (extent * 64 * 8) / 8.0 + 0.0625;
    r = r * 1664525u + 1013904223u;
    double qy = (r >> 8) % (extent * 64 * 8) / 8.0 + 0.0625;
    ASSERT_EQ(winding_at(a, qx, qy), winding_at(b, qx, qy)) << qx << "," << qy;
  }
}

/** No edge of @p fp has its reverse in it too. */
void expect_no_twins(const GSHGEOM_Flat_Path * fp) {
  auto edges = edges_of(fp);
  for (const auto & e : edges) {
    EXPECT_EQ(edges.count(Edge(std::get<2>(e.first), std::get<3>(e.first), std::get<0>(e.first), std::get<1>(e.first))), 0u);
  }
}

} // namespace

TEST(FlatPathSeams, TwoSquaresSharingASideLeaveTheOutlineOfBoth) {
  GSHGEOM_Flat_Path * fp = polygons({{{0, 0}, {10, 0}, {10, 10}, {0, 10}}, {{10, 0}, {20, 0}, {20, 10}, {10, 10}}});
  ASSERT_NE(fp, nullptr);
  GSHGEOM_Flat_Path * out = nullptr;
  ASSERT_EQ(gshgeom_flat_path_cancel_seams(fp, nullptr, nullptr, &out), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_flat_path_point_count(fp), 8u);
  EXPECT_EQ(gshgeom_flat_path_point_count(out), 6u) << "the shared side's two edges are gone";
  EXPECT_EQ(gshgeom_flat_path_subpath_count(out), 1u);
  expect_no_twins(out);
  expect_same_winding(fp, out, 24, 1);
  gshgeom_flat_path_destroy(out);
  gshgeom_flat_path_destroy(fp);
}

TEST(FlatPathSeams, APolygonAndItsReverseCancelWhollyAndLeaveNoSubpath) {
  GSHGEOM_Flat_Path * fp = polygons({{{0, 0}, {8, 0}, {4, 9}}, {{4, 9}, {8, 0}, {0, 0}}});
  GSHGEOM_Flat_Path * out = nullptr;
  ASSERT_EQ(gshgeom_flat_path_cancel_seams(fp, nullptr, nullptr, &out), GSHGEOM_OK);
  ASSERT_NE(out, nullptr);
  EXPECT_EQ(gshgeom_flat_path_subpath_count(out), 0u);
  EXPECT_EQ(gshgeom_flat_path_point_count(out), 0u);
  gshgeom_flat_path_destroy(out);
  gshgeom_flat_path_destroy(fp);
}

TEST(FlatPathSeams, AnEdgeRepeatedInOneDirectionKeepsItsMultiplicityAndWindingTwo) {
  GSHGEOM_Flat_Path * fp = polygons({{{0, 0}, {10, 0}, {10, 10}, {0, 10}}, {{0, 0}, {10, 0}, {10, 10}, {0, 10}}});
  GSHGEOM_Flat_Path * out = nullptr;
  ASSERT_EQ(gshgeom_flat_path_cancel_seams(fp, nullptr, nullptr, &out), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_flat_path_point_count(out), 8u);
  EXPECT_EQ(edges_of(out), edges_of(fp));
  EXPECT_EQ(std::abs(winding_at(out, 5.5 * 64, 5.5 * 64)), 2);
  gshgeom_flat_path_destroy(out);
  gshgeom_flat_path_destroy(fp);
}

TEST(FlatPathSeams, AZeroLengthEdgeIsDropped) {
  GSHGEOM_Flat_Path * fp = polygons({{{0, 0}, {10, 0}, {10, 0}, {10, 10}, {0, 10}}});
  GSHGEOM_Flat_Path * out = nullptr;
  ASSERT_EQ(gshgeom_flat_path_cancel_seams(fp, nullptr, nullptr, &out), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_flat_path_point_count(out), 4u);
  expect_same_winding(fp, out, 12, 2);
  gshgeom_flat_path_destroy(out);
  gshgeom_flat_path_destroy(fp);
}

TEST(FlatPathSeams, ManyOverlappingTrianglesKeepTheirWindingAndLoseTheirTwins) {
  uint32_t r = 12345;
  auto next = [&r](int n) {
    r = r * 1664525u + 1013904223u;
    return (int)((r >> 8) % (uint32_t)n);
  };
  size_t before = 0;
  size_t after = 0;
  for (int round = 0; round < 120; round++) {
    std::vector<std::vector<std::pair<int, int>>> polys;
    int count = 2 + next(40);
    for (int t = 0; t < count; t++) {
      polys.push_back({{next(7), next(7)}, {next(7), next(7)}, {next(7), next(7)}});
    }
    GSHGEOM_Flat_Path * fp = polygons(polys);
    GSHGEOM_Flat_Path * out = nullptr;
    ASSERT_EQ(gshgeom_flat_path_cancel_seams(fp, nullptr, nullptr, &out), GSHGEOM_OK);
    expect_no_twins(out);
    expect_same_winding(fp, out, 8, (uint32_t)round);
    EXPECT_LE(gshgeom_flat_path_point_count(out), gshgeom_flat_path_point_count(fp));
    before += gshgeom_flat_path_point_count(fp);
    after += gshgeom_flat_path_point_count(out);
    GSHGEOM_Flat_Path * again = nullptr;
    ASSERT_EQ(gshgeom_flat_path_cancel_seams(out, nullptr, nullptr, &again), GSHGEOM_OK);
    EXPECT_EQ(gshgeom_flat_path_hash(again), gshgeom_flat_path_hash(out)) << "cancelling what is cancelled changes nothing";
    gshgeom_flat_path_destroy(again);
    gshgeom_flat_path_destroy(out);
    gshgeom_flat_path_destroy(fp);
  }
  EXPECT_LT(after, before) << "random triangles on a small grid share sides";
}

TEST(FlatPathSeams, TheResultDoesNotDependOnTheOrderOfTheSubpaths) {
  std::vector<std::vector<std::pair<int, int>>> polys = {{{0, 0}, {6, 0}, {3, 5}}, {{3, 5}, {6, 0}, {9, 4}}, {{1, 1}, {5, 1}, {5, 5}}, {{0, 0}, {3, 5}, {0, 6}}};
  GSHGEOM_Flat_Path * a = polygons(polys);
  std::reverse(polys.begin(), polys.end());
  GSHGEOM_Flat_Path * b = polygons(polys);
  GSHGEOM_Flat_Path * ra = nullptr;
  GSHGEOM_Flat_Path * rb = nullptr;
  ASSERT_EQ(gshgeom_flat_path_cancel_seams(a, nullptr, nullptr, &ra), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_flat_path_cancel_seams(b, nullptr, nullptr, &rb), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_flat_path_hash(ra), gshgeom_flat_path_hash(rb));
  gshgeom_flat_path_destroy(ra);
  gshgeom_flat_path_destroy(rb);
  gshgeom_flat_path_destroy(a);
  gshgeom_flat_path_destroy(b);
}

TEST(FlatPathSeams, AStrokesOutlineLosesMostOfItsEdgesAndNoneOfItsCoverage) {
  GSHGEOM_Path * mp = nullptr;
  ASSERT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &mp), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_move_to(mp, P(6, 40)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_quadratic_to(mp, P(20, 4), P(34, 40)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(mp, P(50, 8)), GSHGEOM_OK);
  GSHGEOM_Stroke style;
  gshgeom_stroke_init(&style);
  style.width = 5.0;
  style.join = GSHGEOM_STROKE_JOIN_ROUND;
  style.cap = GSHGEOM_STROKE_CAP_ROUND;
  GSHGEOM_Map map;
  ASSERT_EQ(gshgeom_map_identity(2, &map), GSHGEOM_OK);
  GSHGEOM_Stroke_Outline * outline = nullptr;
  ASSERT_EQ(gshgeom_stroke_expand(mp, &style, &map, nullptr, nullptr, &outline), GSHGEOM_OK);
  const GSHGEOM_Flat_Path * body = gshgeom_stroke_outline_body(outline);
  ASSERT_NE(body, nullptr);
  GSHGEOM_Flat_Path * out = nullptr;
  ASSERT_EQ(gshgeom_flat_path_cancel_seams(body, nullptr, nullptr, &out), GSHGEOM_OK);
  EXPECT_LT(gshgeom_flat_path_point_count(out) * 10, gshgeom_flat_path_point_count(body) * 7) << "about two fifths of the edges go";
  expect_no_twins(out);
  expect_same_winding(body, out, 60, 7);
  gshgeom_flat_path_destroy(out);
  gshgeom_stroke_outline_destroy(outline);
  gshgeom_path_destroy(mp);
}

TEST(FlatPathSeams, ABendingStrokeOfStraightAndRoundPiecesCancelsToItsBoundaryAlone) {
  /* its pieces share their seams exactly, so only the outline of their union is left: one ring for the open stroke */
  GSHGEOM_Path * mp = nullptr;
  ASSERT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &mp), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_move_to(mp, P(6, 40)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(mp, P(30, 6)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(mp, P(54, 40)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(mp, P(70, 10)), GSHGEOM_OK);
  GSHGEOM_Stroke style;
  gshgeom_stroke_init(&style);
  style.width = 5.0;
  style.join = GSHGEOM_STROKE_JOIN_ROUND;
  style.cap = GSHGEOM_STROKE_CAP_SQUARE;
  GSHGEOM_Map map;
  ASSERT_EQ(gshgeom_map_identity(2, &map), GSHGEOM_OK);
  GSHGEOM_Stroke_Outline * outline = nullptr;
  ASSERT_EQ(gshgeom_stroke_expand(mp, &style, &map, nullptr, nullptr, &outline), GSHGEOM_OK);
  const GSHGEOM_Flat_Path * body = gshgeom_stroke_outline_body(outline);
  ASSERT_NE(body, nullptr);
  GSHGEOM_Flat_Path * out = nullptr;
  ASSERT_EQ(gshgeom_flat_path_cancel_seams(body, nullptr, nullptr, &out), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_flat_path_subpath_count(out), (size_t)1);
  expect_no_twins(out);
  expect_same_winding(body, out, 60, 11);
  gshgeom_flat_path_destroy(out);
  gshgeom_stroke_outline_destroy(outline);
  gshgeom_path_destroy(mp);
}

TEST(FlatPathSeams, TheArgumentsAndTheLimitAreChecked) {
  GSHGEOM_Flat_Path * fp = polygons({{{0, 0}, {10, 0}, {10, 10}, {0, 10}}});
  GSHGEOM_Flat_Path * out = nullptr;
  EXPECT_EQ(gshgeom_flat_path_cancel_seams(nullptr, nullptr, nullptr, &out), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_flat_path_cancel_seams(fp, nullptr, nullptr, nullptr), GSHGEOM_ERR_INVALID);
  GSHGEOM_Limits limits;
  gshgeom_limits_default(&limits);
  limits.max_points = 3;
  EXPECT_EQ(gshgeom_flat_path_cancel_seams(fp, nullptr, &limits, &out), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(out, nullptr);
  limits.max_points = 4;
  EXPECT_EQ(gshgeom_flat_path_cancel_seams(fp, nullptr, &limits, &out), GSHGEOM_OK);
  gshgeom_flat_path_destroy(out);
  gshgeom_flat_path_destroy(fp);
}

TEST(FlatPathSeams, EveryAllocationFailureIsReportedAndLeaksNothing) {
  GSHGEOM_Flat_Path * fp = polygons({{{0, 0}, {10, 0}, {10, 10}, {0, 10}}, {{10, 0}, {20, 0}, {20, 10}, {10, 10}}, {{3, 3}, {7, 3}, {5, 8}}});
  tracker::Tracker count;
  {
    GSHGEOM_Allocator a = tracker::allocator_for(&count);
    GSHGEOM_Flat_Path * out = nullptr;
    ASSERT_EQ(gshgeom_flat_path_cancel_seams(fp, &a, nullptr, &out), GSHGEOM_OK);
    gshgeom_flat_path_destroy(out);
    EXPECT_TRUE(count.live.empty());
  }
  ASSERT_GT(count.calls, 3);
  for (long k = 0; k < count.calls; k++) {
    for (int persistent = 0; persistent < 2; persistent++) {
      tracker::Tracker t;
      t.fail_from = k;
      t.persistent = persistent != 0;
      GSHGEOM_Allocator a = tracker::allocator_for(&t);
      GSHGEOM_Flat_Path * out = nullptr;
      EXPECT_EQ(gshgeom_flat_path_cancel_seams(fp, &a, nullptr, &out), GSHGEOM_ERR_OOM) << k;
      EXPECT_EQ(out, nullptr);
      EXPECT_TRUE(t.live.empty()) << "k=" << k << " leaked";
    }
  }
  gshgeom_flat_path_destroy(fp);
}

TEST(StrokeOutlineReduced, EachPartIsReducedBesideTheOutlineWhichIsNotChanged) {
  GSHGEOM_Path * mp = nullptr;
  ASSERT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &mp), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_move_to(mp, P(8, 8)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(mp, P(48, 10)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_quadratic_to(mp, P(60, 40), P(30, 50)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_close(mp), GSHGEOM_OK);
  for (uint32_t alignment : {(uint32_t)GSHGEOM_STROKE_ALIGN_CENTER, (uint32_t)GSHGEOM_STROKE_ALIGN_INSIDE}) {
    GSHGEOM_Stroke style;
    gshgeom_stroke_init(&style);
    style.width = 6.0;
    style.join = GSHGEOM_STROKE_JOIN_ROUND;
    style.alignment = alignment;
    GSHGEOM_Map map;
    ASSERT_EQ(gshgeom_map_identity(2, &map), GSHGEOM_OK);
    tracker::Tracker t;
    GSHGEOM_Allocator a = tracker::allocator_for(&t);
    GSHGEOM_Stroke_Outline * o = nullptr;
    ASSERT_EQ(gshgeom_stroke_expand(mp, &style, &map, &a, nullptr, &o), GSHGEOM_OK);
    uint64_t hash = gshgeom_stroke_outline_hash(o);
    EXPECT_EQ(gshgeom_stroke_outline_body_reduced(o), nullptr) << "not before it is asked for";
    ASSERT_EQ(gshgeom_stroke_outline_reduce(o, nullptr), GSHGEOM_OK);
    const GSHGEOM_Flat_Path * whole[3] = {gshgeom_stroke_outline_body(o), gshgeom_stroke_outline_aligned(o), gshgeom_stroke_outline_clip(o)};
    const GSHGEOM_Flat_Path * reduced[3] = {gshgeom_stroke_outline_body_reduced(o), gshgeom_stroke_outline_aligned_reduced(o),
        gshgeom_stroke_outline_clip_reduced(o)};
    int present = 0;
    for (int i = 0; i < 3; i++) {
      EXPECT_EQ(whole[i] == nullptr, reduced[i] == nullptr) << "a part is reduced if and only if it is there";
      if (whole[i]) {
        present++;
        EXPECT_LE(gshgeom_flat_path_point_count(reduced[i]), gshgeom_flat_path_point_count(whole[i]));
        expect_no_twins(reduced[i]);
        expect_same_winding(whole[i], reduced[i], 70, 3 + i);
      }
    }
    EXPECT_GE(present, 1);
    EXPECT_EQ(gshgeom_stroke_outline_hash(o), hash) << "the outline is the same outline";
    const GSHGEOM_Flat_Path * first = gshgeom_stroke_outline_body_reduced(o);
    EXPECT_EQ(gshgeom_stroke_outline_reduce(o, nullptr), GSHGEOM_OK);
    EXPECT_EQ(gshgeom_stroke_outline_body_reduced(o), first) << "a second call does nothing";
    gshgeom_stroke_outline_destroy(o);
    EXPECT_TRUE(t.live.empty()) << "the reduced parts are destroyed with the outline";
  }
  EXPECT_EQ(gshgeom_stroke_outline_reduce(nullptr, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_stroke_outline_body_reduced(nullptr), nullptr);
  gshgeom_path_destroy(mp);
}

TEST(StrokeOutlineReduced, AnAllocationFailureLeavesTheOutlineWithoutAReducedFormWhateverPartItFailsIn) {
  // An open stroke has one part; a closed one aligned inside has three, so a failure in the second or the third must free the
  // reduced parts already made.
  for (int closed = 0; closed < 2; closed++) {
    GSHGEOM_Path * mp = nullptr;
    ASSERT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &mp), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_path_move_to(mp, P(8, 8)), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_path_line_to(mp, P(48, 10)), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_path_line_to(mp, P(30, 50)), GSHGEOM_OK);
    if (closed) {
      ASSERT_EQ(gshgeom_path_close(mp), GSHGEOM_OK);
    }
    GSHGEOM_Stroke style;
    gshgeom_stroke_init(&style);
    style.width = 6.0;
    style.alignment = closed ? GSHGEOM_STROKE_ALIGN_INSIDE : GSHGEOM_STROKE_ALIGN_CENTER;
    GSHGEOM_Map map;
    ASSERT_EQ(gshgeom_map_identity(2, &map), GSHGEOM_OK);
    tracker::Tracker count;
    long made = 0;
    {
      GSHGEOM_Allocator a = tracker::allocator_for(&count);
      GSHGEOM_Stroke_Outline * o = nullptr;
      ASSERT_EQ(gshgeom_stroke_expand(mp, &style, &map, &a, nullptr, &o), GSHGEOM_OK);
      EXPECT_EQ(gshgeom_stroke_outline_aligned(o) != nullptr, closed != 0);
      made = count.calls;
      ASSERT_EQ(gshgeom_stroke_outline_reduce(o, nullptr), GSHGEOM_OK);
      gshgeom_stroke_outline_destroy(o);
    }
    ASSERT_GT(count.calls, made);
    for (long k = made; k < count.calls; k++) {
      tracker::Tracker t;
      t.fail_from = k;
      t.persistent = true;
      GSHGEOM_Allocator a = tracker::allocator_for(&t);
      GSHGEOM_Stroke_Outline * o = nullptr;
      ASSERT_EQ(gshgeom_stroke_expand(mp, &style, &map, &a, nullptr, &o), GSHGEOM_OK);
      EXPECT_EQ(gshgeom_stroke_outline_reduce(o, nullptr), GSHGEOM_ERR_OOM) << k;
      EXPECT_EQ(gshgeom_stroke_outline_body_reduced(o), nullptr);
      EXPECT_EQ(gshgeom_stroke_outline_aligned_reduced(o), nullptr);
      t.fail_from = -1;
      EXPECT_EQ(gshgeom_stroke_outline_reduce(o, nullptr), GSHGEOM_OK) << "it may be asked for again";
      gshgeom_stroke_outline_destroy(o);
      EXPECT_TRUE(t.live.empty()) << "closed=" << closed << " k=" << k << " leaked";
    }
    gshgeom_path_destroy(mp);
  }
}

namespace {

/** @p k triangles that each pass through the origin (their other corners are all different, so no edge is shared). */
GSHGEOM_Flat_Path * loops_through_one_vertex(size_t k) {
  std::vector<Dp> pts;
  std::vector<GSHGEOM_Flat_Subpath> subs;
  for (size_t i = 0; i < k; i++) {
    // the corner before the origin is the leftmost of its triangle, so the walk of each loop reaches the origin from outside
    pts.push_back(Dp{int32_t(-64 - int32_t(i)), int32_t(i % 97) * 64 + 64});
    pts.push_back(Dp{0, 0});
    pts.push_back(Dp{int32_t(64 + 2 * int32_t(i)), int32_t(i % 97) * 64 + 200});
    subs.push_back(GSHGEOM_Flat_Subpath{3 * i, 3, 1});
  }
  GSHGEOM_Flat_Path * fp = nullptr;
  EXPECT_EQ(gshgeom_flat_path_from_points(pts.data(), pts.size(), subs.data(), subs.size(), nullptr, nullptr, &fp), GSHGEOM_OK);
  return fp;
}

}  // namespace

// Many closed subpaths through one point: each walk past the point looked for an unused arc leaving it past every used one, which was
// the square of the number of subpaths (160 000 of them took 7 seconds, for an outline of half a million points).  The walk is
// linear now and its answer is the old one's: the walk takes the first unused arc leaving a point, so the loops through the
// vertex come back joined into one closed polyline of all their points, and its hash is the one the quadratic walk gave.
TEST(FlatSeams, ManyLoopsThroughOneVertexAreWalkedInLinearTimeAndGiveTheAnswerOfTheQuadraticWalk) {
  const struct {
    size_t k;
    uint64_t hash;
  } pinned[] = {{3, UINT64_C(0x2cbb33696175f33c)}, {500, UINT64_C(0x589f079c6ffb102c)}, {5000, UINT64_C(0x901ee5e6c86d53ae)},
      {20000, UINT64_C(0x126a2f8f0ceb5262)}};
  for (const auto & c : pinned) {
    GSHGEOM_Flat_Path * fp = loops_through_one_vertex(c.k);
    GSHGEOM_Flat_Path * out = nullptr;
    ASSERT_EQ(gshgeom_flat_path_cancel_seams(fp, nullptr, nullptr, &out), GSHGEOM_OK) << c.k;
    EXPECT_EQ(gshgeom_flat_path_subpath_count(out), 1u) << c.k;
    EXPECT_EQ(gshgeom_flat_path_point_count(out), 3 * c.k);
    EXPECT_EQ(gshgeom_flat_path_hash(out), c.hash) << c.k;
    gshgeom_flat_path_destroy(out);
    gshgeom_flat_path_destroy(fp);
  }
  GSHGEOM_Flat_Path * fp = loops_through_one_vertex(160000);
  GSHGEOM_Flat_Path * out = nullptr;
  const auto t0 = std::chrono::steady_clock::now();
  ASSERT_EQ(gshgeom_flat_path_cancel_seams(fp, nullptr, nullptr, &out), GSHGEOM_OK);
  const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  EXPECT_EQ(gshgeom_flat_path_subpath_count(out), 1u);
  EXPECT_EQ(gshgeom_flat_path_point_count(out), 480000u);
  // a guard, not the requirement: 7 seconds before, a tenth of a second now; a sanitizer is several times slower
  if (!under_valgrind()) {
    EXPECT_LT(seconds, 4.0);
  }
  gshgeom_flat_path_destroy(out);
  gshgeom_flat_path_destroy(fp);
}

// A failed allocation anywhere in the walk (the skip links are one) is an error and leaks nothing.
TEST(FlatSeams, AFailedAllocationInTheWalkOfLoopsThroughOneVertexLeaksNothing) {
  GSHGEOM_Flat_Path * fp = loops_through_one_vertex(40);
  for (long fail = 0; fail < 40; fail++) {
    tracker::Tracker t;
    t.fail_from = fail;
    const GSHGEOM_Allocator alloc = tracker::allocator_for(&t);
    GSHGEOM_Flat_Path * out = nullptr;
    const GSHGEOM_Result r = gshgeom_flat_path_cancel_seams(fp, &alloc, nullptr, &out);
    if (t.failed) {
      EXPECT_EQ(r, GSHGEOM_ERR_OOM) << fail;
      EXPECT_EQ(out, nullptr);
    } else {
      EXPECT_EQ(r, GSHGEOM_OK) << fail;
      gshgeom_flat_path_destroy(out);
    }
    EXPECT_TRUE(t.live.empty()) << fail;
  }
  gshgeom_flat_path_destroy(fp);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
