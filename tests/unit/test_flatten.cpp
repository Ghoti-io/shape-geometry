/**
 * @file
 *
 * The integer flattener, revision 2, and the flattened path (CAP-3, CAP-13):
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
#include "../flatten_vectors.h"
#include "../known_answers.h"
#include "../dump_sink.h"

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

double dist_to_segment(double px, double py, const Dp & a, const Dp & b) {
  double dx = (double)b.x - a.x, dy = (double)b.y - a.y;
  double len2 = dx * dx + dy * dy;
  double u = len2 > 0 ? ((px - a.x) * dx + (py - a.y) * dy) / len2 : 0.0;
  u = u < 0 ? 0 : (u > 1 ? 1 : u);
  double qx = a.x + u * dx - px, qy = a.y + u * dy - py;
  return std::sqrt(qx * qx + qy * qy);
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
  // d = |2 p1 - p0 - p2| summed over axes: 12 units (4 * (4 - 1)) stays one
  // line, 14 splits.
  struct Case {
    int32_t p1x, p1y;
    size_t points;
  };
  const Case cases[] = {{6, 6, 2}, {6, 7, 3}, {7, 6, 3}, {0, 0, 2}, {6, 5, 2}, {6, -6, 2}, {-7, -6, 3}};
  for (const Case & c : cases) {
    GSHGEOM_Device_Path * dp = dev_path();
    ASSERT_EQ(gshgeom_device_path_move_to(dp, D(0, 0)), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_device_path_quadratic_to(dp, D(c.p1x, c.p1y), D(12, 0)), GSHGEOM_OK);
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
  const CCase ccases[] = {{10, 20, true}, {14, 20, true}, {15, 20, false}, {10, 24, true},
      {10, 25, false}, {6, 20, true}, {5, 20, false}, {10, 16, true}, {10, 15, false}};
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
  EXPECT_EQ(GSHGEOM_FLATTEN_MAX_DEPTH, 20);
}

TEST(FlattenRevision, TheFunctionTheMacroAndTheVectorsAgree) {
  EXPECT_EQ(GSHGEOM_FLATTEN_REVISION, 3);
  EXPECT_EQ(gshgeom_flatten_revision(), GSHGEOM_FLATTEN_REVISION);
  EXPECT_EQ(kFlattenVectorsRevision, gshgeom_flatten_revision())
      << "tests/flatten_vectors.h was recorded from another revision";
}

TEST(FlattenRevision, TheRevisionIsInNoHash) {
  // A flattened path hashes its points, not the revision: the same points hash
  // the same, whatever produced them (AD-21; the format version stays 1).
  GSHGEOM_Device_Path * dp = dev_path();
  ASSERT_EQ(gshgeom_device_path_move_to(dp, D(0, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_device_path_line_to(dp, D(64, 0)), GSHGEOM_OK);
  GSHGEOM_Flat_Path * fp = flatten(dp);
  const Dp pts[2] = {D(0, 0), D(64, 0)};
  // The hash of a path of two points, spelled out here: version 1, kind FLAT,
  // one subpath, open, two points.
  uint64_t h = UINT64_C(0xcbf29ce484222325);
  auto byte = [&h](uint8_t b) { h = (h ^ b) * UINT64_C(0x100000001b3); };
  auto u32 = [&](uint32_t v) { for (int i = 0; i < 4; i++) byte((v >> (8 * i)) & 0xff); };
  auto u64 = [&](uint64_t v) { for (int i = 0; i < 8; i++) byte((v >> (8 * i)) & 0xff); };
  u32(1);
  u32(0x54414c46u);
  u64(1);
  byte(0);
  u64(2);
  for (const Dp & q : pts) {
    u32((uint32_t)q.x);
    u32((uint32_t)q.y);
  }
  EXPECT_EQ(gshgeom_flat_path_hash(fp), h);
  gshgeom_flat_path_destroy(fp);
  gshgeom_device_path_destroy(dp);
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

// ---------------------------------------------------------------------------
// The guarantee (flattener revision 2)
// ---------------------------------------------------------------------------
//
// src/flatten/flatten.c proves one property: for every curve with int32
// control points, every point of the curve is within 4 - 0.12 units of the
// polyline, and the depth cap never decides a piece.  Pixels cannot test that,
// so these tests evaluate the curve exactly (the Bernstein form in `long
// double` at dyadic parameters, which are exact in the type; the only error is
// the arithmetic of the sum, far below the 0.12 margin) at thirty-three
// parameters inside every piece the flattener made, and assert the distance to
// that piece's segment.  A piece-by-piece measure is used because a curve
// spanning the int32 range is 2^33 units long, so no affordable uniform
// sampling would see a unit of deviation between two samples; inside a piece
// the deviation is a smooth function of the parameter.
//
// Which parameter interval a piece covers is not part of the library's output,
// so a model of the recursion in this file (a second statement of the
// algorithm, in int64) recovers it; the model is first required to produce
// exactly the library's vertices, so a difference in either is a failure.

namespace {

struct WP {
  int64_t x, y;
};

struct Piece {
  Dp at;
  uint32_t leaf;
  int level;
};

int64_t fdiv(int64_t a, int64_t b) {
  int64_t q = a / b;
  if (a % b != 0 && a < 0) {
    q -= 1;
  }
  return q;
}

WP wmid(WP a, WP b) {
  return WP{fdiv(a.x + b.x + 1, 2), fdiv(a.y + b.y + 1, 2)};
}

WP widen(Dp p) {
  return WP{(int64_t)p.x * 256, (int64_t)p.y * 256};
}

Dp narrow(WP p) {
  return Dp{(int32_t)fdiv(p.x + 128, 256), (int32_t)fdiv(p.y + 128, 256)};
}

const int64_t kThreshold = 4 * 256 * (GSHGEOM_FLATTEN_TOLERANCE - 1);
const int kDepth = 20;

struct Model {
  std::vector<Piece> out;
  int forced = 0;
};

void model_quad(Model & m, WP p0, WP p1, WP p2, uint32_t leaf, int level) {
  int64_t sum = std::llabs(2 * p1.x - p0.x - p2.x) + std::llabs(2 * p1.y - p0.y - p2.y);
  if (sum <= kThreshold || level >= kDepth) {
    m.forced += sum > kThreshold;
    m.out.push_back({narrow(p2), leaf, level});
    return;
  }
  WP a = wmid(p0, p1), b = wmid(p1, p2), c = wmid(a, b);
  model_quad(m, p0, a, c, leaf * 2, level + 1);
  model_quad(m, c, b, p2, leaf * 2 + 1, level + 1);
}

void model_cubic(Model & m, WP p0, WP p1, WP p2, WP p3, uint32_t leaf, int level) {
  int64_t d1 = std::llabs(3 * p1.x - 2 * p0.x - p3.x) + std::llabs(3 * p1.y - 2 * p0.y - p3.y);
  int64_t d2 = std::llabs(3 * p2.x - p0.x - 2 * p3.x) + std::llabs(3 * p2.y - p0.y - 2 * p3.y);
  bool flat = d1 <= kThreshold && d2 <= kThreshold;
  if (flat || level >= kDepth) {
    m.forced += !flat;
    m.out.push_back({narrow(p3), leaf, level});
    return;
  }
  WP a = wmid(p0, p1), b = wmid(p1, p2), c = wmid(p2, p3);
  WP d = wmid(a, b), e = wmid(b, c), mm = wmid(d, e);
  model_cubic(m, p0, a, d, mm, leaf * 2, level + 1);
  model_cubic(m, mm, e, c, p3, leaf * 2 + 1, level + 1);
}

struct Curve {
  bool cubic;
  Dp p[4];
};

Model model_of(const Curve & c) {
  Model m;
  if (c.cubic) {
    model_cubic(m, widen(c.p[0]), widen(c.p[1]), widen(c.p[2]), widen(c.p[3]), 0, 0);
  } else {
    model_quad(m, widen(c.p[0]), widen(c.p[1]), widen(c.p[2]), 0, 0);
  }
  return m;
}

GSHGEOM_Device_Path * path_of(const Curve & c) {
  GSHGEOM_Device_Path * dp = dev_path();
  EXPECT_EQ(gshgeom_device_path_move_to(dp, c.p[0]), GSHGEOM_OK);
  if (c.cubic) {
    EXPECT_EQ(gshgeom_device_path_cubic_to(dp, c.p[1], c.p[2], c.p[3]), GSHGEOM_OK);
  } else {
    EXPECT_EQ(gshgeom_device_path_quadratic_to(dp, c.p[1], c.p[2]), GSHGEOM_OK);
  }
  return dp;
}

struct LVec {
  long double x, y;
};

LVec evaluate(const Curve & c, long double s) {
  const long double r = 1.0L - s;
  long double b[4];
  if (c.cubic) {
    b[0] = r * r * r;
    b[1] = 3 * s * r * r;
    b[2] = 3 * s * s * r;
    b[3] = s * s * s;
  } else {
    b[0] = r * r;
    b[1] = 2 * s * r;
    b[2] = s * s;
    b[3] = 0;
  }
  LVec out{0, 0};
  for (int i = 0; i < (c.cubic ? 4 : 3); ++i) {
    out.x += b[i] * c.p[i].x;
    out.y += b[i] * c.p[i].y;
  }
  return out;
}

long double seg_distance(LVec q, Dp a, Dp b) {
  const long double dx = (long double)b.x - a.x, dy = (long double)b.y - a.y;
  const long double len2 = dx * dx + dy * dy;
  long double along = 0;
  if (len2 > 0) {
    along = ((q.x - a.x) * dx + (q.y - a.y) * dy) / len2;
    along = along < 0 ? 0 : (along > 1 ? 1 : along);
  }
  const long double ex = q.x - (a.x + along * dx), ey = q.y - (a.y + along * dy);
  return std::sqrt(ex * ex + ey * ey);
}

const int kSamples = 32;

long double piece_deviation(const Curve & c, const std::vector<Piece> & pieces) {
  long double worst = 0;
  Dp from = c.p[0];
  for (const Piece & v : pieces) {
    const long double span = std::ldexp(1.0L, v.level);
    for (int i = 0; i <= kSamples; ++i) {
      const long double s = ((long double)v.leaf + (long double)i / kSamples) / span;
      worst = std::max(worst, seg_distance(evaluate(c, s), from, v.at));
    }
    from = v.at;
  }
  return worst;
}

// The pieces tile [0, 1]: each starts where the one before ended.
bool tiles(const std::vector<Piece> & pieces) {
  if (pieces.empty() || pieces[0].leaf != 0) {
    return false;
  }
  for (size_t i = 0; i + 1 < pieces.size(); i++) {
    const Piece &a = pieces[i], &b = pieces[i + 1];
    if (((uint64_t)a.leaf + 1) << b.level != (uint64_t)b.leaf << a.level) {
      return false;
    }
  }
  const Piece & last = pieces.back();
  return (uint64_t)last.leaf + 1 == (uint64_t)1 << last.level;
}

// 4 - 0.12, with the arithmetic's own slack.
const long double kBound = (long double)GSHGEOM_FLATTEN_TOLERANCE - 0.12L + 1e-4L;

// Flatten with the library, require the model to agree vertex for vertex,
// assert the bound; `worst` is the smallest margin seen (the proof says at
// least 0.12).
void check(const Curve & c, long double & worst, size_t & most) {
  GSHGEOM_Device_Path * dp = path_of(c);
  GSHGEOM_Flat_Path * fp = flatten(dp);
  ASSERT_NE(fp, nullptr);
  Model m = model_of(c);
  const size_t n = gshgeom_flat_path_point_count(fp);
  const Dp * pts = gshgeom_flat_path_points(fp);
  ASSERT_EQ(n, m.out.size() + 1) << "the model and the library disagree on the vertex count";
  for (size_t i = 0; i < m.out.size(); i++) {
    ASSERT_TRUE(pts[i + 1].x == m.out[i].at.x && pts[i + 1].y == m.out[i].at.y)
        << "vertex " << i << " differs from the model";
  }
  EXPECT_EQ(m.forced, 0) << "the depth cap decided a piece";
  EXPECT_TRUE(tiles(m.out)) << "the pieces do not tile the parameter interval";
  const long double d = piece_deviation(c, m.out);
  EXPECT_LE(d, kBound) << (c.cubic ? "cubic" : "quadratic") << " strays " << (double)d
                       << " units; first point " << c.p[0].x << "," << c.p[0].y;
  worst = std::min(worst, (long double)GSHGEOM_FLATTEN_TOLERANCE - d);
  most = std::max(most, m.out.size());
  gshgeom_flat_path_destroy(fp);
  gshgeom_device_path_destroy(dp);
}

struct Random {
  uint64_t state;
  uint64_t next() {
    state ^= state >> 12;
    state ^= state << 25;
    state ^= state >> 27;
    return state * UINT64_C(0x2545F4914F6CDD1D);
  }
  // A coordinate in [-2^bits, 2^bits), the whole int32 range for 31.
  int32_t coordinate(int bits) {
    const uint64_t span = (uint64_t)1 << (bits + 1);
    return (int32_t)((int64_t)(next() % span) - ((int64_t)1 << bits));
  }
};

Curve random_curve(Random & r, bool cubic, int bits) {
  Curve c{cubic, {}};
  for (int i = 0; i < (cubic ? 4 : 3); ++i) {
    c.p[i] = D(r.coordinate(bits), r.coordinate(bits));
  }
  return c;
}

Curve corner_curve(Random & r, bool cubic) {
  Curve c{cubic, {}};
  for (int i = 0; i < (cubic ? 4 : 3); ++i) {
    c.p[i] = D((r.next() & 1u) ? INT32_MAX : INT32_MIN, (r.next() & 1u) ? INT32_MAX : INT32_MIN);
  }
  return c;
}

}  // namespace

TEST(FlattenTolerance, EveryPieceOfEveryCurveIsWithinTheToleranceLessAFewHundredthsOfTheCurve) {
  Random random{UINT64_C(0x9E3779B97F4A7C15)};
  const int classes[] = {4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 31};
  for (const int bits : classes) {
    const int curves = bits <= 16 ? 40 : (bits <= 24 ? 10 : 4);
    long double worst = 1e9L;
    size_t most = 0;
    for (int i = 0; i < curves; ++i) {
      check(random_curve(random, false, bits), worst, most);
      check(random_curve(random, true, bits), worst, most);
    }
    RecordProperty("smallest_margin_below_2_pow_" + std::to_string(bits), std::to_string((double)worst));
    EXPECT_GE(worst, 0.12L) << "coordinates below 2^" << bits;
  }
}

TEST(FlattenTolerance, TheExtremeCurvesOfTheInt32SquareAreWithinTheTolerance) {
  Random random{UINT64_C(0x0123456789ABCDEF)};
  long double worst = 1e9L;
  size_t most = 0;
  for (int i = 0; i < 8; ++i) {
    check(corner_curve(random, false), worst, most);
    check(corner_curve(random, true), worst, most);
  }
  EXPECT_GE(worst, 0.12L);
  EXPECT_GE(most, 10000u) << "the corners make the longest polylines";
  // The two vectors the vertex lists below pin, by hand.
  const Curve q{false, {D(INT32_MIN, INT32_MAX), D(INT32_MAX, INT32_MIN), D(INT32_MIN, INT32_MIN)}};
  const Curve k{true, {D(INT32_MIN, INT32_MIN), D(INT32_MAX, INT32_MIN), D(INT32_MIN, INT32_MAX), D(INT32_MAX, INT32_MAX)}};
  check(q, worst, most);
  check(k, worst, most);
}

// A sum that wraps in 32 bits must not be taken for a small one: with p0 = p3
// and one control point 2^24 units away, 3 * 2^24 * 256 is a multiple of 2^32,
// so a 32-bit d1 (or d2) reads 0 and the curve would be called flat.
TEST(FlattenTolerance, ASumThatWrapsIn32BitsIsNotTakenForASmallOne) {
  const int32_t a = 1 << 24;
  const Curve curves[] = {
      {true, {D(0, 0), D(a, 0), D(0, 0), D(0, 0)}},
      {true, {D(0, 0), D(0, 0), D(a, 0), D(0, 0)}},
      {true, {D(0, 0), D(0, a), D(0, 0), D(0, 0)}},
      {true, {D(0, 0), D(0, 0), D(0, a), D(0, 0)}},
      {false, {D(0, 0), D(a, 0), D(0, 0), D(0, 0)}},
  };
  long double worst = 1e9L;
  size_t most = 0;
  for (const Curve & c : curves) {
    check(c, worst, most);
    EXPECT_GT(most, 100u);
  }
}

// Vertex lists worked out by an independent model: Python, unbounded integers,
// tools/gen-hash-vectors.py (flatten_model), sharing no code with the C.  Each
// is the vertex count, the FNV-1a 64 of the vertices, and the first and last
// vertex.  They pin the rounding rules - ties upward at both the midpoints and
// the narrowing, floors of negative values - that the distance test is too
// coarse to notice (a floor midpoint at eight extra bits stays inside the
// bound).  The first, second and fourth are also the lists `font`'s own test
// pins for its revision 2 at tolerance 4.
TEST(FlattenTolerance, TheVertexListsAreWhatTheIndependentModelSays) {
  struct Known {
    bool cubic;
    Dp p[4];
    size_t count;
    uint64_t hash;
    Dp first, last;
  };
  const Known known[] = {
      {false, {D(0, 0), D(3001, 4999), D(8003, 1), D(0, 0)}, 32, UINT64_C(0x5e9fe3e61234e30b), D(190, 303), D(8003, 1)},
      {false, {D(-1000, -777), D(-3001, 5001), D(777, -9999), D(0, 0)}, 64, UINT64_C(0x4bbe0f147b91fd7d), D(-1061, -602), D(777, -9999)},
      {false, {D(0, 0), D(400, -1200), D(800, 0), D(0, 0)}, 16, UINT64_C(0x90b0755175c5a1bc), D(50, -141), D(800, 0)},
      {true, {D(0, 0), D(1000, 3001), D(-2001, 4001), D(5000, -7)}, 59, UINT64_C(0xdd37fa8b0432a31f), D(44, 139), D(5000, -7)},
      {true, {D(-31, 17), D(-8003, 2999), D(9001, -6001), D(64, -64)}, 98, UINT64_C(0xc7393744da8497e9), D(-213, 85), D(64, -64)},
      {false, {D(INT32_MIN, INT32_MAX), D(INT32_MAX, INT32_MIN), D(INT32_MIN, INT32_MIN), D(0, 0)}, 32768,
          UINT64_C(0xf6d691a5fcc16620), D(-2147221512, 2147221507), D(INT32_MIN, INT32_MIN)},
      {true, {D(INT32_MIN, INT32_MIN), D(INT32_MAX, INT32_MIN), D(INT32_MIN, INT32_MAX), D(INT32_MAX, INT32_MAX)}, 53054,
          UINT64_C(0xce3c26488e495407), D(-2147287046, -2147483645), D(INT32_MAX, INT32_MAX)},
      {false, {D(-3, -3), D(-1, -1), D(1, 1), D(0, 0)}, 1, UINT64_C(0x29c7dd317360ac35), D(1, 1), D(1, 1)},
      {false, {D(0, 0), D(1, 13), D(2, 0), D(0, 0)}, 2, UINT64_C(0x05e8f83edd3d7f31), D(1, 7), D(2, 0)},
  };
  for (const Known & k : known) {
    Curve c{k.cubic, {k.p[0], k.p[1], k.p[2], k.p[3]}};
    GSHGEOM_Device_Path * dp = path_of(c);
    GSHGEOM_Flat_Path * fp = flatten(dp);
    ASSERT_EQ(gshgeom_flat_path_point_count(fp), k.count + 1);
    const Dp * pts = gshgeom_flat_path_points(fp);
    EXPECT_EQ(points_hash(pts + 1, k.count), k.hash);
    EXPECT_TRUE(pts[1].x == k.first.x && pts[1].y == k.first.y);
    EXPECT_TRUE(pts[k.count].x == k.last.x && pts[k.count].y == k.last.y);
    gshgeom_flat_path_destroy(fp);
    gshgeom_device_path_destroy(dp);
  }
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
    EXPECT_LE(worst, (double)GSHGEOM_FLATTEN_TOLERANCE - 0.12 + 1e-3) << "scale " << scale;
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

TEST(FlatPathDump, EveryWriteIsChecked) {
  GSHGEOM_Device_Path * dp = flat_fixture_device();
  GSHGEOM_Flat_Path * fp = flatten(dp);
  GSHGEOM_EVERY_WRITE_CHECKED(gshgeom_flat_path_dump, fp, "flat path");
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
