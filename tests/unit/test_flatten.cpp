/**
 * @file
 *
 * The integer flattener and the flattened path (CAP-3, CAP-13): bit-for-bit
 * agreement with `font`'s flattener on recorded vectors, the tolerance in
 * device pixels against dense exact evaluation, flattening after the map, the
 * structure of the result, the limits, and the hash and the dump pinned by
 * known answers.
 *
 * Copyright 2026 by Corey Pennycuff
 */

#include "test_helpers.h"

#include "../alloc_tracker.h"
#include "../flatten_vectors.h"
#include "../known_answers.h"

#include <algorithm>
#include <clocale>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <functional>
#include <string>
#include <vector>

namespace {

using Dp = GSHGEOM_Device_Point;

Dp D(int32_t x, int32_t y) {
  return Dp{x, y};
}

GSHGEOM_Point P(double x, double y) {
  return GSHGEOM_Point{x, y, 0.0};
}

std::string dump_of(const GSHGEOM_Flat_Path * p) {
  FILE * f = std::tmpfile();
  EXPECT_NE(f, nullptr);
  EXPECT_EQ(gshgeom_flat_path_dump(p, f), GSHGEOM_OK);
  std::rewind(f);
  std::string s;
  char buf[256];
  size_t n;
  while ((n = std::fread(buf, 1, sizeof buf, f)) > 0) {
    s.append(buf, n);
  }
  std::fclose(f);
  return s;
}

// FNV-1a 64 of the points as int32 little-endian pairs: this file's own
// statement of the spelling, the one the vectors were hashed with.
uint64_t points_hash(const Dp * p, size_t n) {
  uint64_t h = UINT64_C(0xcbf29ce484222325);
  auto mix = [&h](int32_t v) {
    uint32_t u = (uint32_t)v;
    for (int i = 0; i < 4; i++) {
      h ^= (u >> (8 * i)) & 0xff;
      h *= UINT64_C(0x100000001b3);
    }
  };
  for (size_t i = 0; i < n; i++) {
    mix(p[i].x);
    mix(p[i].y);
  }
  return h;
}

GSHGEOM_Device_Path * dev_path() {
  GSHGEOM_Device_Path * p = nullptr;
  EXPECT_EQ(gshgeom_device_path_create(nullptr, nullptr, &p), GSHGEOM_OK);
  return p;
}

GSHGEOM_Flat_Path * flatten(const GSHGEOM_Device_Path * dp, const GSHGEOM_Limits * l = nullptr) {
  GSHGEOM_Flat_Path * fp = nullptr;
  EXPECT_EQ(gshgeom_device_path_flatten(dp, nullptr, l, &fp), GSHGEOM_OK);
  return fp;
}

// One curve as a device path: move to p[0], then the quadratic or cubic.
GSHGEOM_Device_Path * curve_path(const FlattenVector & v) {
  GSHGEOM_Device_Path * dp = dev_path();
  EXPECT_EQ(gshgeom_device_path_move_to(dp, D(v.p[0][0], v.p[0][1])), GSHGEOM_OK);
  if (v.kind == 2) {
    EXPECT_EQ(gshgeom_device_path_quadratic_to(dp, D(v.p[1][0], v.p[1][1]), D(v.p[2][0], v.p[2][1])),
        GSHGEOM_OK);
  } else {
    EXPECT_EQ(gshgeom_device_path_cubic_to(dp, D(v.p[1][0], v.p[1][1]), D(v.p[2][0], v.p[2][1]),
                  D(v.p[3][0], v.p[3][1])),
        GSHGEOM_OK);
  }
  return dp;
}

// The exact numerators of a curve at t = j / 1024, in 64-bit integers (the vectors'
// coordinates are below 2^27, so a numerator is below 2^57); the division by a
// power of two is exact in a double.
void exact_curve_point(const FlattenVector & v, int j, double * x, double * y) {
  const int64_t N = 1024;
  int64_t s = N - j;
  int64_t t = j;
  int64_t nx, ny, den;
  if (v.kind == 2) {
    int64_t b0 = s * s, b1 = 2 * s * t, b2 = t * t;
    nx = b0 * v.p[0][0] + b1 * v.p[1][0] + b2 * v.p[2][0];
    ny = b0 * v.p[0][1] + b1 * v.p[1][1] + b2 * v.p[2][1];
    den = N * N;
  } else {
    int64_t b0 = s * s * s, b1 = 3 * s * s * t, b2 = 3 * s * t * t, b3 = t * t * t;
    nx = b0 * v.p[0][0] + b1 * v.p[1][0] + b2 * v.p[2][0] + b3 * v.p[3][0];
    ny = b0 * v.p[0][1] + b1 * v.p[1][1] + b2 * v.p[2][1] + b3 * v.p[3][1];
    den = N * N * N;
  }
  *x = (double)nx / (double)den;
  *y = (double)ny / (double)den;
}

double dist_to_segment(double px, double py, const Dp & a, const Dp & b) {
  double dx = (double)b.x - a.x, dy = (double)b.y - a.y;
  double len2 = dx * dx + dy * dy;
  double u = len2 > 0 ? ((px - a.x) * dx + (py - a.y) * dy) / len2 : 0.0;
  u = u < 0 ? 0 : (u > 1 ? 1 : u);
  double qx = a.x + u * dx - px, qy = a.y + u * dy - py;
  return std::sqrt(qx * qx + qy * qy);
}

// The largest distance from the exact curve to the flattened polyline.
double deviation(const FlattenVector & v, const Dp * pts, size_t n) {
  double worst = 0;
  for (int j = 0; j <= 1024; j++) {
    double x, y;
    exact_curve_point(v, j, &x, &y);
    double best = 1e300;
    for (size_t i = 0; i + 1 < n; i++) {
      best = std::min(best, dist_to_segment(x, y, pts[i], pts[i + 1]));
    }
    worst = std::max(worst, best);
  }
  return worst;
}

}  // namespace

// ---------------------------------------------------------------------------
// Against font
// ---------------------------------------------------------------------------

TEST(FlattenAgainstFont, EveryRecordedVectorMatchesBitForBit) {
  ASSERT_GE(kFlattenVectorCount, 900u);
  unsigned explicit_lists = 0, curves_with_many_lines = 0;
  for (unsigned i = 0; i < kFlattenVectorCount; i++) {
    const FlattenVector & v = kFlattenVectors[i];
    GSHGEOM_Device_Path * dp = curve_path(v);
    GSHGEOM_Flat_Path * fp = flatten(dp);
    ASSERT_NE(fp, nullptr) << "vector " << i;
    // The flattened path is the start point, then font's end points.
    ASSERT_EQ(gshgeom_flat_path_subpath_count(fp), 1u);
    ASSERT_EQ(gshgeom_flat_path_point_count(fp), (size_t)v.count + 1) << "vector " << i;
    const Dp * pts = gshgeom_flat_path_points(fp);
    EXPECT_EQ(pts[0].x, v.p[0][0]);
    EXPECT_EQ(pts[0].y, v.p[0][1]);
    EXPECT_EQ(points_hash(pts + 1, v.count), v.hash) << "vector " << i << " kind " << v.kind;
    if (v.offset != 0xffffffffu) {
      explicit_lists++;
      for (unsigned k = 0; k < v.count; k++) {
        ASSERT_EQ(pts[1 + k].x, kFlattenPoints[v.offset + 2 * k]) << "vector " << i << " point " << k;
        ASSERT_EQ(pts[1 + k].y, kFlattenPoints[v.offset + 2 * k + 1]) << "vector " << i << " point " << k;
      }
    }
    const int last = v.kind == 2 ? 2 : 3;
    EXPECT_EQ(pts[v.count].x, v.p[last][0]) << "the last line ends at the curve's end";
    EXPECT_EQ(pts[v.count].y, v.p[last][1]);
    if (v.count > 16) {
      curves_with_many_lines++;
    }
    gshgeom_flat_path_destroy(fp);
    gshgeom_device_path_destroy(dp);
  }
  EXPECT_GE(explicit_lists, 50u);
  EXPECT_GE(curves_with_many_lines, 100u) << "the vectors reach deep subdivisions";
}

TEST(FlattenAgainstFont, TheThresholdIsWhereFontPutsIt) {
  // d = |2 p1 - p0 - p2| summed over axes: 16 stays one line, 18 splits.
  struct Case {
    int32_t p1x, p1y;
    size_t points;
  };
  const Case cases[] = {{8, 8, 2}, {8, 9, 3}, {9, 8, 3}, {0, 0, 2}, {8, 7, 2}, {8, -8, 2}, {-9, -8, 3}};
  for (const Case & c : cases) {
    GSHGEOM_Device_Path * dp = dev_path();
    ASSERT_EQ(gshgeom_device_path_move_to(dp, D(0, 0)), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_device_path_quadratic_to(dp, D(c.p1x, c.p1y), D(16, 0)), GSHGEOM_OK);
    GSHGEOM_Flat_Path * fp = flatten(dp);
    EXPECT_EQ(gshgeom_flat_path_point_count(fp) == 2, c.points == 2) << c.p1x << "," << c.p1y;
    gshgeom_flat_path_destroy(fp);
    gshgeom_device_path_destroy(dp);
  }
  // The cubic needs both control points inside: d1 and d2.
  struct CCase {
    int32_t c1x, c2x;
    bool one_line;
  };
  // p0 = (0,0), p3 = (30, 0): d1 = |3 c1x - 30|, d2 = |3 c2x - 60|, y = 0.
  const CCase ccases[] = {{10, 20, true}, {15, 20, true}, {16, 20, false}, {10, 25, true},
      {10, 26, false}, {5, 20, true}, {4, 20, false}, {10, 15, true}, {10, 14, false}};
  for (const CCase & c : ccases) {
    GSHGEOM_Device_Path * dp = dev_path();
    ASSERT_EQ(gshgeom_device_path_move_to(dp, D(0, 0)), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_device_path_cubic_to(dp, D(c.c1x, 0), D(c.c2x, 0), D(30, 0)), GSHGEOM_OK);
    GSHGEOM_Flat_Path * fp = flatten(dp);
    EXPECT_EQ(gshgeom_flat_path_point_count(fp) == 2, c.one_line) << c.c1x << " " << c.c2x;
    gshgeom_flat_path_destroy(fp);
    gshgeom_device_path_destroy(dp);
  }
  EXPECT_EQ(GSHGEOM_FLATTEN_TOLERANCE, 4);
  EXPECT_EQ(GSHGEOM_FLATTEN_MAX_DEPTH, 16);
}

TEST(FlattenAgainstFont, CoordinatesPastWhatFontsInt32SumsHandleDoNotOverflow) {
  // 2 p1 - p0 - p2 overflows an int32 here; the 64-bit sum does not, and the
  // curve is split rather than treated as flat.
  GSHGEOM_Device_Path * dp = dev_path();
  ASSERT_EQ(gshgeom_device_path_move_to(dp, D(INT32_MIN, INT32_MIN)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_device_path_quadratic_to(dp, D(INT32_MAX, INT32_MAX), D(INT32_MIN, INT32_MIN)),
      GSHGEOM_OK);
  GSHGEOM_Flat_Path * fp = flatten(dp);
  ASSERT_NE(fp, nullptr);
  EXPECT_GT(gshgeom_flat_path_point_count(fp), 100u);
  const Dp * pts = gshgeom_flat_path_points(fp);
  for (size_t i = 0; i < gshgeom_flat_path_point_count(fp); i++) {
    ASSERT_GE(pts[i].x, INT32_MIN);
    ASSERT_LE(pts[i].x, INT32_MAX);
  }
  gshgeom_flat_path_destroy(fp);
  gshgeom_device_path_destroy(dp);
  // And a cubic across the whole range.
  dp = dev_path();
  ASSERT_EQ(gshgeom_device_path_move_to(dp, D(INT32_MIN, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_device_path_cubic_to(dp, D(INT32_MIN, INT32_MAX), D(INT32_MAX, INT32_MIN), D(INT32_MAX, 0)),
      GSHGEOM_OK);
  fp = flatten(dp);
  ASSERT_NE(fp, nullptr);
  EXPECT_GT(gshgeom_flat_path_point_count(fp), 100u);
  gshgeom_flat_path_destroy(fp);
  gshgeom_device_path_destroy(dp);
}

// ---------------------------------------------------------------------------
// The tolerance
// ---------------------------------------------------------------------------

TEST(FlattenTolerance, DeviationFromTheExactCurveIsTheThresholdPlusTheIntegerMidpointsDrift) {
  // What is asserted is what dense exact evaluation of every vector measures.
  // The recursion stops when a piece deviates by at most the tolerance (4 units,
  // 1/16 pixel) from its chord, and that holds for the curve through the
  // integer control points of the piece.  But each integer midpoint is a floor,
  // so the pieces' control points are not the true curve's, and the drift adds
  // up with the depth: a small curve stays within about the tolerance, a curve
  // spanning thousands of pixels strays by up to about three times it.  `font`
  // flattens exactly this way, bit for bit, and the pure move keeps it, so this
  // is the definition's real accuracy; the bound for each size is the measured
  // worst case plus a margin (story 7's notes give the table).
  struct Bucket {
    int bits;  // the largest coordinate is below 2^bits
    double bound;
  };
  const Bucket buckets[] = {{4, 4.0}, {8, 5.5}, {16, 8.5}, {28, 12.5}};
  double worst[4] = {0, 0, 0, 0};
  unsigned count[4] = {0, 0, 0, 0};
  for (unsigned i = 0; i < kFlattenVectorCount; i++) {
    const FlattenVector & v = kFlattenVectors[i];
    int64_t biggest = 0;
    for (int k = 0; k < (v.kind == 2 ? 3 : 4); k++) {
      biggest = std::max<int64_t>(biggest, std::llabs((int64_t)v.p[k][0]));
      biggest = std::max<int64_t>(biggest, std::llabs((int64_t)v.p[k][1]));
    }
    size_t b = 0;
    while (b < 3 && biggest >= ((int64_t)1 << buckets[b].bits)) {
      b++;
    }
    GSHGEOM_Device_Path * dp = curve_path(v);
    GSHGEOM_Flat_Path * fp = flatten(dp);
    double d = deviation(v, gshgeom_flat_path_points(fp), gshgeom_flat_path_point_count(fp));
    EXPECT_LE(d, buckets[b].bound) << "vector " << i << " coordinates below 2^" << buckets[b].bits;
    worst[b] = std::max(worst[b], d);
    count[b]++;
    gshgeom_flat_path_destroy(fp);
    gshgeom_device_path_destroy(dp);
  }
  for (size_t b = 0; b < 4; b++) {
    EXPECT_GT(count[b], 20u) << "bucket " << b;
    RecordProperty("worst_deviation_below_2_pow_" + std::to_string(buckets[b].bits), std::to_string(worst[b]));
  }
  // The smallest curves meet the tolerance itself.
  EXPECT_LE(worst[0], (double)GSHGEOM_FLATTEN_TOLERANCE);
}

TEST(FlattenTolerance, FlatteningHappensAfterTheMap) {
  // One model quadratic whose deviation from its chord is 0.25 model units, under
  // scales that make it 0.0025, 0.25 and 250 pixels.  After the map the error is
  // in device pixels, so the small one is one line and the large one is many.
  GSHGEOM_Path * mp = nullptr;
  ASSERT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &mp), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_move_to(mp, P(0, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_quadratic_to(mp, P(0.5, 0.5), P(1, 0)), GSHGEOM_OK);
  size_t prev = 0;
  for (double scale : {0.01, 1.0, 10.0, 1000.0}) {
    GSHGEOM_Map map;
    ASSERT_EQ(gshgeom_map_scale(2, scale, scale, 1.0, &map), GSHGEOM_OK);
    GSHGEOM_Flat_Path * fp = nullptr;
    ASSERT_EQ(gshgeom_path_flatten(mp, &map, nullptr, nullptr, &fp), GSHGEOM_OK);
    size_t n = gshgeom_flat_path_point_count(fp);
    if (scale == 0.01) {
      EXPECT_EQ(n, 2u) << "a deviation of 1/64 px needs no split";
    }
    EXPECT_GT(n, prev) << scale;
    prev = n;
    // Deviation of the polyline from the true curve in device units.
    const Dp * pts = gshgeom_flat_path_points(fp);
    double worst = 0;
    for (int j = 0; j <= 2000; j++) {
      double t = j / 2000.0;
      double x = 2 * t * (1 - t) * 0.5 + t * t * 1.0;
      double y = 2 * t * (1 - t) * 0.5;
      x *= scale * 64;
      y *= scale * 64;
      double best = 1e300;
      for (size_t i = 0; i + 1 < n; i++) {
        best = std::min(best, dist_to_segment(x, y, pts[i], pts[i + 1]));
      }
      worst = std::max(worst, best);
    }
    EXPECT_LE(worst, (double)GSHGEOM_FLATTEN_TOLERANCE + 0.5) << "scale " << scale;
    gshgeom_flat_path_destroy(fp);
  }
  gshgeom_path_destroy(mp);
}

// ---------------------------------------------------------------------------
// Structure
// ---------------------------------------------------------------------------

namespace {

// The device path whose flattening is FLAT_FIXTURE of tools/gen-hash-vectors.py.
GSHGEOM_Device_Path * flat_fixture_device() {
  GSHGEOM_Device_Path * dp = dev_path();
  EXPECT_EQ(gshgeom_device_path_move_to(dp, D(0, 0)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_device_path_line_to(dp, D(64, 0)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_device_path_quadratic_to(dp, D(96, 0), D(128, 0)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_device_path_cubic_to(dp, D(150, 1), D(170, 2), D(192, 3)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_device_path_move_to(dp, D(-1000, 1000)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_device_path_line_to(dp, D(1100, 1000)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_device_path_line_to(dp, D(-5, -7)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_device_path_close(dp), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_device_path_move_to(dp, D(5, 5)), GSHGEOM_OK);
  return dp;
}

}  // namespace

TEST(FlatPath, SubpathsPointsAndTheClosedFlag) {
  GSHGEOM_Device_Path * dp = flat_fixture_device();
  GSHGEOM_Flat_Path * fp = flatten(dp);
  ASSERT_EQ(gshgeom_flat_path_subpath_count(fp), 3u);
  ASSERT_EQ(gshgeom_flat_path_point_count(fp), 4u + 3u + 1u);
  GSHGEOM_Flat_Subpath s;
  ASSERT_EQ(gshgeom_flat_path_subpath(fp, 0, &s), GSHGEOM_OK);
  EXPECT_EQ(s.first, 0u);
  EXPECT_EQ(s.count, 4u);
  EXPECT_EQ(s.closed, 0);
  ASSERT_EQ(gshgeom_flat_path_subpath(fp, 1, &s), GSHGEOM_OK);
  EXPECT_EQ(s.first, 4u);
  EXPECT_EQ(s.count, 3u);
  EXPECT_EQ(s.closed, 1);
  ASSERT_EQ(gshgeom_flat_path_subpath(fp, 2, &s), GSHGEOM_OK) << "a lone move is a subpath";
  EXPECT_EQ(s.first, 7u);
  EXPECT_EQ(s.count, 1u);
  EXPECT_EQ(s.closed, 0);
  EXPECT_EQ(gshgeom_flat_path_subpath(fp, 3, &s), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_flat_path_subpath(nullptr, 0, &s), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_flat_path_subpath(fp, 0, nullptr), GSHGEOM_ERR_INVALID);
  const Dp * pts = gshgeom_flat_path_points(fp);
  EXPECT_EQ(pts[3].x, 192);
  EXPECT_EQ(pts[3].y, 3);
  EXPECT_EQ(pts[5].x, 1100);
  EXPECT_EQ(pts[7].x, 5);
  gshgeom_flat_path_destroy(fp);
  gshgeom_device_path_destroy(dp);
}

TEST(FlatPath, AVerbAfterACloseStartsANewSubpathAtTheClosedStart) {
  GSHGEOM_Device_Path * dp = dev_path();
  ASSERT_EQ(gshgeom_device_path_move_to(dp, D(10, 20)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_device_path_line_to(dp, D(30, 20)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_device_path_close(dp), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_device_path_line_to(dp, D(50, 60)), GSHGEOM_OK);
  GSHGEOM_Flat_Path * fp = flatten(dp);
  ASSERT_EQ(gshgeom_flat_path_subpath_count(fp), 2u);
  GSHGEOM_Flat_Subpath s;
  ASSERT_EQ(gshgeom_flat_path_subpath(fp, 1, &s), GSHGEOM_OK);
  EXPECT_EQ(s.count, 2u);
  EXPECT_EQ(s.closed, 0);
  const Dp * pts = gshgeom_flat_path_points(fp);
  EXPECT_EQ(pts[s.first].x, 10);
  EXPECT_EQ(pts[s.first].y, 20);
  EXPECT_EQ(pts[s.first + 1].x, 50);
  gshgeom_flat_path_destroy(fp);
  gshgeom_device_path_destroy(dp);
}

TEST(FlatPath, EmptyAndNull) {
  GSHGEOM_Device_Path * dp = dev_path();
  GSHGEOM_Flat_Path * fp = flatten(dp);
  EXPECT_EQ(gshgeom_flat_path_subpath_count(fp), 0u);
  EXPECT_EQ(gshgeom_flat_path_point_count(fp), 0u);
  EXPECT_EQ(gshgeom_flat_path_points(fp), nullptr);
  gshgeom_flat_path_destroy(fp);
  gshgeom_device_path_destroy(dp);
  EXPECT_EQ(gshgeom_flat_path_point_count(nullptr), 0u);
  EXPECT_EQ(gshgeom_flat_path_subpath_count(nullptr), 0u);
  EXPECT_EQ(gshgeom_flat_path_points(nullptr), nullptr);
  gshgeom_flat_path_destroy(nullptr);
  GSHGEOM_Flat_Path * out = nullptr;
  EXPECT_EQ(gshgeom_device_path_flatten(nullptr, nullptr, nullptr, &out), GSHGEOM_ERR_INVALID);
  dp = dev_path();
  EXPECT_EQ(gshgeom_device_path_flatten(dp, nullptr, nullptr, nullptr), GSHGEOM_ERR_INVALID);
  gshgeom_device_path_destroy(dp);
  GSHGEOM_Map map;
  ASSERT_EQ(gshgeom_map_identity(2, &map), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_flatten(nullptr, &map, nullptr, nullptr, &out), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_flatten(nullptr, &map, nullptr, nullptr, nullptr), GSHGEOM_ERR_INVALID);
  GSHGEOM_Path * mp = nullptr;
  ASSERT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &mp), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_move_to(mp, P(1, 1)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_flatten(mp, &map, nullptr, nullptr, nullptr), GSHGEOM_ERR_INVALID);
  gshgeom_path_destroy(mp);
  EXPECT_EQ(out, nullptr);
}

TEST(FlatPath, ACurveStartsWhereTheLastOneEnded) {
  // A cubic, then a quadratic, then a line: each starts at the previous end.
  GSHGEOM_Device_Path * dp = dev_path();
  ASSERT_EQ(gshgeom_device_path_move_to(dp, D(0, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_device_path_cubic_to(dp, D(0, 640), D(640, 640), D(640, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_device_path_quadratic_to(dp, D(900, -600), D(1200, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_device_path_cubic_to(dp, D(1300, 400), D(1500, 400), D(1600, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_device_path_line_to(dp, D(1700, 50)), GSHGEOM_OK);
  GSHGEOM_Flat_Path * fp = flatten(dp);
  size_t n = gshgeom_flat_path_point_count(fp);
  const Dp * pts = gshgeom_flat_path_points(fp);
  // The three curves flattened one by one from their own starts.
  size_t at = 0;
  const FlattenVector parts[] = {
      {3, {{0, 0}, {0, 640}, {640, 640}, {640, 0}}, 0, 0, 0},
      {2, {{640, 0}, {900, -600}, {1200, 0}, {0, 0}}, 0, 0, 0},
      {3, {{1200, 0}, {1300, 400}, {1500, 400}, {1600, 0}}, 0, 0, 0},
  };
  at = 1;
  for (const FlattenVector & v : parts) {
    GSHGEOM_Device_Path * one = curve_path(v);
    GSHGEOM_Flat_Path * f1 = flatten(one);
    size_t m = gshgeom_flat_path_point_count(f1);
    const Dp * q = gshgeom_flat_path_points(f1);
    for (size_t k = 1; k < m; k++) {
      ASSERT_LT(at, n);
      EXPECT_EQ(pts[at].x, q[k].x);
      EXPECT_EQ(pts[at].y, q[k].y);
      at++;
    }
    gshgeom_flat_path_destroy(f1);
    gshgeom_device_path_destroy(one);
  }
  EXPECT_EQ(at + 1, n) << "then the line's end point";
  EXPECT_EQ(pts[n - 1].x, 1700);
  gshgeom_flat_path_destroy(fp);
  gshgeom_device_path_destroy(dp);
}

TEST(FlatPath, ALineKeepsEveryPointEvenIfRepeated) {
  GSHGEOM_Device_Path * dp = dev_path();
  ASSERT_EQ(gshgeom_device_path_move_to(dp, D(1, 1)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_device_path_line_to(dp, D(1, 1)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_device_path_line_to(dp, D(1, 1)), GSHGEOM_OK);
  GSHGEOM_Flat_Path * fp = flatten(dp);
  EXPECT_EQ(gshgeom_flat_path_point_count(fp), 3u);
  gshgeom_flat_path_destroy(fp);
  gshgeom_device_path_destroy(dp);
}

// ---------------------------------------------------------------------------
// Limits and allocation
// ---------------------------------------------------------------------------

TEST(FlatPathLimits, TheCapOnPointsIsExact) {
  for (size_t cap : {1u, 2u, 7u, 100u}) {
    GSHGEOM_Device_Path * dp = dev_path();
    ASSERT_EQ(gshgeom_device_path_move_to(dp, D(0, 0)), GSHGEOM_OK);
    for (size_t i = 1; i < cap; i++) {
      ASSERT_EQ(gshgeom_device_path_line_to(dp, D((int32_t)i, 0)), GSHGEOM_OK);
    }
    GSHGEOM_Limits l;
    gshgeom_limits_default(&l);
    l.max_points = cap;
    GSHGEOM_Flat_Path * fp = nullptr;
    ASSERT_EQ(gshgeom_device_path_flatten(dp, nullptr, &l, &fp), GSHGEOM_OK) << cap;
    EXPECT_EQ(gshgeom_flat_path_point_count(fp), cap);
    gshgeom_flat_path_destroy(fp);
    ASSERT_EQ(gshgeom_device_path_line_to(dp, D(999, 0)), GSHGEOM_OK);
    fp = nullptr;
    EXPECT_EQ(gshgeom_device_path_flatten(dp, nullptr, &l, &fp), GSHGEOM_ERR_LIMIT) << cap;
    EXPECT_EQ(fp, nullptr);
    gshgeom_device_path_destroy(dp);
  }
}

TEST(FlatPathLimits, ACurveThatNeedsMoreIsRefusedAtTheBound) {
  const FlattenVector * big = nullptr;
  for (unsigned i = 0; i < kFlattenVectorCount; i++) {
    if (kFlattenVectors[i].count > 200 && kFlattenVectors[i].count < 400) {
      big = &kFlattenVectors[i];
      break;
    }
  }
  ASSERT_NE(big, nullptr);
  GSHGEOM_Device_Path * dp = curve_path(*big);
  GSHGEOM_Limits l;
  gshgeom_limits_default(&l);
  l.max_points = big->count + 1;  // the start point and every line
  GSHGEOM_Flat_Path * fp = nullptr;
  ASSERT_EQ(gshgeom_device_path_flatten(dp, nullptr, &l, &fp), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_flat_path_point_count(fp), (size_t)big->count + 1);
  gshgeom_flat_path_destroy(fp);
  l.max_points = big->count;
  fp = nullptr;
  EXPECT_EQ(gshgeom_device_path_flatten(dp, nullptr, &l, &fp), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(fp, nullptr);
  // Zero means the default.
  l.max_points = 0;
  ASSERT_EQ(gshgeom_device_path_flatten(dp, nullptr, &l, &fp), GSHGEOM_OK);
  gshgeom_flat_path_destroy(fp);
  gshgeom_device_path_destroy(dp);
}

TEST(FlatPathLimits, ThroughTheModelPathBothCapsApply) {
  GSHGEOM_Path * mp = nullptr;
  ASSERT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &mp), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_move_to(mp, P(0, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(mp, P(1, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(mp, P(1, 1)), GSHGEOM_OK);
  GSHGEOM_Map map;
  ASSERT_EQ(gshgeom_map_identity(2, &map), GSHGEOM_OK);
  GSHGEOM_Limits l;
  gshgeom_limits_default(&l);
  GSHGEOM_Flat_Path * fp = nullptr;
  l.max_elements = 3;
  l.max_points = 3;
  ASSERT_EQ(gshgeom_path_flatten(mp, &map, nullptr, &l, &fp), GSHGEOM_OK);
  gshgeom_flat_path_destroy(fp);
  fp = nullptr;
  l.max_elements = 2;
  EXPECT_EQ(gshgeom_path_flatten(mp, &map, nullptr, &l, &fp), GSHGEOM_ERR_LIMIT);
  l.max_elements = 3;
  l.max_points = 2;
  EXPECT_EQ(gshgeom_path_flatten(mp, &map, nullptr, &l, &fp), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(fp, nullptr);
  gshgeom_path_destroy(mp);
}

TEST(FlatPathAllocation, AFailedAllocationIsOomAndLeaksNothing) {
  GSHGEOM_Device_Path * dp = flat_fixture_device();
  // A quadratic and a cubic whose flattening makes the point array grow more
  // than once, so a failure falls inside the recursion, in a left half or a right.
  const FlattenVector * quad_v = nullptr;
  const FlattenVector * cubic_v = nullptr;
  for (unsigned i = 0; i < kFlattenVectorCount; i++) {
    const FlattenVector & v = kFlattenVectors[i];
    if (v.count > 100 && v.count < 300) {
      if (v.kind == 2 && !quad_v) quad_v = &v;
      if (v.kind == 3 && !cubic_v) cubic_v = &v;
    }
  }
  ASSERT_NE(quad_v, nullptr);
  ASSERT_NE(cubic_v, nullptr);
  GSHGEOM_Device_Path * quad_dp = curve_path(*quad_v);
  GSHGEOM_Device_Path * cubic_dp = curve_path(*cubic_v);
  for (GSHGEOM_Device_Path * path : {dp, quad_dp, cubic_dp}) {
    tracker::Tracker count;
    {
      GSHGEOM_Allocator a = tracker::allocator_for(&count);
      GSHGEOM_Flat_Path * fp = nullptr;
      ASSERT_EQ(gshgeom_device_path_flatten(path, &a, nullptr, &fp), GSHGEOM_OK);
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
        EXPECT_EQ(gshgeom_device_path_flatten(path, &a, nullptr, &fp), GSHGEOM_ERR_OOM) << k;
        EXPECT_EQ(fp, nullptr);
        EXPECT_TRUE(t.live.empty()) << "k=" << k << " leaked";
      }
    }
  }
  gshgeom_device_path_destroy(cubic_dp);
  gshgeom_device_path_destroy(quad_dp);
  gshgeom_device_path_destroy(dp);
}

TEST(FlatPathAllocation, TheOneCallLeaksNothingEither) {
  GSHGEOM_Path * mp = nullptr;
  ASSERT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &mp), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_move_to(mp, P(1, 1)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_cubic_to(mp, P(10, 90), P(80, 90), P(90, 1)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_close(mp), GSHGEOM_OK);
  GSHGEOM_Map map;
  ASSERT_EQ(gshgeom_map_scale(2, 4.0, 4.0, 1.0, &map), GSHGEOM_OK);
  tracker::Tracker count;
  {
    GSHGEOM_Allocator a = tracker::allocator_for(&count);
    GSHGEOM_Flat_Path * fp = nullptr;
    ASSERT_EQ(gshgeom_path_flatten(mp, &map, &a, nullptr, &fp), GSHGEOM_OK);
    gshgeom_flat_path_destroy(fp);
    EXPECT_TRUE(count.live.empty()) << "the intermediate device path is freed on success";
  }
  ASSERT_GT(count.calls, 4);
  for (long k = 0; k < count.calls; k++) {
    for (int persistent = 0; persistent < 2; persistent++) {
      tracker::Tracker t;
      t.fail_from = k;
      t.persistent = persistent != 0;
      GSHGEOM_Allocator a = tracker::allocator_for(&t);
      GSHGEOM_Flat_Path * fp = nullptr;
      EXPECT_EQ(gshgeom_path_flatten(mp, &map, &a, nullptr, &fp), GSHGEOM_ERR_OOM) << k;
      EXPECT_EQ(fp, nullptr);
      EXPECT_TRUE(t.live.empty()) << "k=" << k << " leaked";
    }
  }
  gshgeom_path_destroy(mp);
}

// ---------------------------------------------------------------------------
// Hash and dump
// ---------------------------------------------------------------------------

TEST(FlatPathHash, MatchesTheIndependentImplementation) {
  GSHGEOM_Device_Path * dp = flat_fixture_device();
  GSHGEOM_Flat_Path * fp = flatten(dp);
  EXPECT_EQ(gshgeom_flat_path_hash(fp), kFlatPathHash);
  EXPECT_EQ(dump_of(fp), kFlatPathDump);
  gshgeom_flat_path_destroy(fp);
  gshgeom_device_path_destroy(dp);
  EXPECT_EQ(gshgeom_flat_path_hash(nullptr), GSHGEOM_FNV_OFFSET_BASIS);
  // The empty path hashes the version, the kind and a zero count.
  dp = dev_path();
  fp = flatten(dp);
  GSHGEOM_Hash h;
  gshgeom_hash_init(&h);
  gshgeom_hash_u32(&h, 1);
  gshgeom_hash_u32(&h, 0x54414c46u);
  gshgeom_hash_u64(&h, 0);
  EXPECT_EQ(gshgeom_flat_path_hash(fp), gshgeom_hash_value(&h));
  gshgeom_flat_path_destroy(fp);
  gshgeom_device_path_destroy(dp);
}

TEST(FlatPathHash, EveryFieldMatters) {
  auto hash_of = [](int variant) {
    GSHGEOM_Device_Path * dp = dev_path();
    gshgeom_device_path_move_to(dp, D(1, 2));
    gshgeom_device_path_line_to(dp, D(3, 4));
    switch (variant) {
    case 1: gshgeom_device_path_line_to(dp, D(3, 4)); break;
    case 2: gshgeom_device_path_close(dp); break;
    case 3: gshgeom_device_path_move_to(dp, D(3, 4)); break;
    case 4: gshgeom_device_path_close(dp); gshgeom_device_path_close(dp); break;
    case 5: gshgeom_device_path_move_to(dp, D(5, 5)); gshgeom_device_path_line_to(dp, D(6, 6)); break;
    default: break;
    }
    GSHGEOM_Flat_Path * fp = flatten(dp);
    uint64_t h = gshgeom_flat_path_hash(fp);
    gshgeom_flat_path_destroy(fp);
    gshgeom_device_path_destroy(dp);
    return h;
  };
  std::vector<uint64_t> seen;
  for (int v = 0; v < 6; v++) {
    if (v == 4) {
      EXPECT_EQ(hash_of(4), hash_of(2)) << "a repeated close adds nothing to the flattened path";
      continue;
    }
    uint64_t h = hash_of(v);
    for (uint64_t o : seen) EXPECT_NE(h, o) << v;
    seen.push_back(h);
  }
}

TEST(FlatPathDump, TheHashIsAlwaysSixteenDigits) {
  bool leading_zero = false;
  for (int i = 0; i < 400; i++) {
    GSHGEOM_Device_Path * dp = dev_path();
    ASSERT_EQ(gshgeom_device_path_move_to(dp, D(i, 0)), GSHGEOM_OK);
    GSHGEOM_Flat_Path * fp = flatten(dp);
    leading_zero = leading_zero || gshgeom_flat_path_hash(fp) < (UINT64_C(1) << 60);
    std::string d = dump_of(fp);
    size_t at = d.find("hash=");
    ASSERT_NE(at, std::string::npos);
    size_t end = d.find('\n', at);
    EXPECT_EQ(end - (at + 5), 16u) << d;
    gshgeom_flat_path_destroy(fp);
    gshgeom_device_path_destroy(dp);
  }
  EXPECT_TRUE(leading_zero) << "the search found no hash with a leading zero";
}

TEST(FlatPathDump, AFailingSinkIsReportedAndTheLocaleDoesNotMatter) {
  GSHGEOM_Device_Path * dp = flat_fixture_device();
  GSHGEOM_Flat_Path * fp = flatten(dp);
  FILE * full = std::fopen("/dev/full", "w");
  if (full) {
    EXPECT_EQ(gshgeom_flat_path_dump(fp, full), GSHGEOM_ERR_IO);
    std::fclose(full);
  }
  FILE * ro = std::fopen("/dev/null", "r");
  ASSERT_NE(ro, nullptr);
  EXPECT_EQ(gshgeom_flat_path_dump(fp, ro), GSHGEOM_ERR_IO);
  std::fclose(ro);
  EXPECT_EQ(gshgeom_flat_path_dump(fp, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_flat_path_dump(nullptr, stdout), GSHGEOM_ERR_INVALID);
  const char * old = std::setlocale(LC_ALL, nullptr);
  std::string saved = old ? old : "C";
  for (const char * name : {"de_DE.UTF-8", "fr_FR.UTF-8", "de_DE", "C.UTF-8"}) {
    if (std::setlocale(LC_ALL, name)) {
      EXPECT_EQ(dump_of(fp), kFlatPathDump) << name;
    }
  }
  std::setlocale(LC_ALL, saved.c_str());
  gshgeom_flat_path_destroy(fp);
  gshgeom_device_path_destroy(dp);
}

TEST(FlatPathHash, AWholeModelPathGivesTheSameHashEverywhere) {
  // The fixture of the device tests, flattened through the model path: the value
  // is what the cross-target probe is asked to reproduce.
  GSHGEOM_Path * mp = nullptr;
  ASSERT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &mp), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_move_to(mp, P(1.5, 2.25)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_cubic_to(mp, P(40.5, -3), P(60, 55.75), P(99.25, 30)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_close(mp), GSHGEOM_OK);
  GSHGEOM_Map map;
  ASSERT_EQ(gshgeom_map_identity(2, &map), GSHGEOM_OK);
  GSHGEOM_Flat_Path *a = nullptr, *b = nullptr;
  ASSERT_EQ(gshgeom_path_flatten(mp, &map, nullptr, nullptr, &a), GSHGEOM_OK);
  GSHGEOM_Device_Path * dp = nullptr;
  ASSERT_EQ(gshgeom_device_path_from_model(mp, &map, nullptr, nullptr, &dp), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_device_path_flatten(dp, nullptr, nullptr, &b), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_flat_path_hash(a), gshgeom_flat_path_hash(b)) << "the one call is the two calls";
  gshgeom_flat_path_destroy(a);
  gshgeom_flat_path_destroy(b);
  gshgeom_device_path_destroy(dp);
  gshgeom_path_destroy(mp);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
