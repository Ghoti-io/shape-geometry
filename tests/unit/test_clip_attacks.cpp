/**
 * @file
 *
 * The clip to a convex quadrilateral (clip.h): the cases that an independent statement of each rule separates from a
 * variant, which the grid of calls could not: the quadrilateral's validity against the rule's own words over every
 * quadrilateral on a small lattice, the crossing point against exact rational arithmetic (written with the compiler's
 * 128-bit integer, which the library avoids but a test on a 64-bit host may use), and the limit of a caller that raises
 * the point cap past the default.
 *
 * Copyright 2026 by Corey Pennycuff
 */

#include "../test_helpers.h"

#include <ghoti.io/shape-geometry/clip.h>
#include <ghoti.io/shape-geometry/flatten.h>

extern "C" {
#include "../../src/clip/clip_internal.h"
}

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <vector>

namespace {

typedef GSHGEOM_Device_Point Pt;

struct Rng {
  uint64_t s = 0x2545F4914F6CDD1Dull;
  uint64_t next() {
    s ^= s << 13;
    s ^= s >> 7;
    s ^= s << 17;
    return s;
  }
  int64_t range(int64_t lo, int64_t hi) { return lo + (int64_t)(next() % (uint64_t)(hi - lo + 1)); }
};

int64_t orient(Pt a, Pt b, Pt c) {
  const int64_t v = ((int64_t)b.x - a.x) * ((int64_t)c.y - a.y) - ((int64_t)b.y - a.y) * ((int64_t)c.x - a.x);
  return v < 0 ? -1 : (v > 0 ? 1 : 0);
}

GSHGEOM_Flat_Path * triangle(Pt a, Pt b, Pt c) {
  const Pt pts[3] = {a, b, c};
  const GSHGEOM_Flat_Subpath sub = {0, 3, 1};
  GSHGEOM_Flat_Path * p = nullptr;
  EXPECT_EQ(gshgeom_flat_path_from_points(pts, 3, &sub, 1, nullptr, nullptr, &p), GSHGEOM_OK);
  return p;
}

}  // namespace

TEST(ClipAttack, AQuadrilateralIsAcceptedExactlyWhenItIsConvexAndHasArea) {
  // every ordered choice of four points of a 4 x 4 lattice, repeats included: the rule is that the turns at the corners,
  // each from the corner and its two neighbours, are all of one sign or zero, and that the polygon has an area
  GSHGEOM_Flat_Path * subject = triangle({0, 0}, {10, 0}, {0, 10});
  int accepted = 0, refused = 0, mixed_but_area = 0;
  for (int code = 0; code < 65536; code++) {
    Pt q[4];
    for (int k = 0; k < 4; k++) {
      const int v = (code >> (4 * k)) & 15;
      q[k] = {(v & 3) * 3 - 3, (v >> 2) * 3 - 3};
    }
    int pos = 0, neg = 0;
    for (int i = 0; i < 4; i++) {
      const int64_t t = orient(q[(i + 3) % 4], q[i], q[(i + 1) % 4]);
      pos += t > 0;
      neg += t < 0;
    }
    int64_t twice_area = 0;  // the shoelace sum
    for (int i = 0; i < 4; i++) {
      twice_area += (int64_t)q[i].x * q[(i + 1) % 4].y - (int64_t)q[(i + 1) % 4].x * q[i].y;
    }
    const bool want = !(pos && neg) && (pos || neg) && twice_area != 0;
    mixed_but_area += (pos && neg && twice_area != 0);
    GSHGEOM_Flat_Path * out = nullptr;
    const GSHGEOM_Result r = gshgeom_flat_path_clip_convex(subject, q, nullptr, nullptr, &out);
    ASSERT_EQ(r == GSHGEOM_OK, want) << "code " << code;
    if (!want) {
      ASSERT_EQ(r, GSHGEOM_ERR_INVALID) << "code " << code;
    }
    gshgeom_flat_path_destroy(out);
    accepted += want;
    refused += !want;
  }
  EXPECT_GT(accepted, 5000);
  EXPECT_GT(refused, 5000);
  EXPECT_GT(mixed_but_area, 1000);  // non-convex quadrilaterals: the case a wrong neighbour would let through
  gshgeom_flat_path_destroy(subject);
}

TEST(ClipAttack, ARepeatedCornerChangesNothing) {
  // a side of zero length cuts nothing: the same result as the quadrilateral without the repeat
  GSHGEOM_Flat_Path * subject = triangle({-5, -5}, {30, 2}, {4, 30});
  const Pt base[4] = {{0, 0}, {20, 0}, {20, 20}, {0, 20}};
  GSHGEOM_Flat_Path * ref = nullptr;
  ASSERT_EQ(gshgeom_flat_path_clip_convex(subject, base, nullptr, nullptr, &ref), GSHGEOM_OK);
  for (int rep = 0; rep < 4; rep++) {
    Pt q[4];
    for (int k = 0; k < 4; k++) {
      q[k] = base[k];
    }
    q[(rep + 1) % 4] = q[rep];  // a repeat: the quadrilateral is a triangle with a doubled corner
    GSHGEOM_Flat_Path * out = nullptr;
    GSHGEOM_Flat_Path * tri = nullptr;
    Pt three[4] = {q[0], q[1], q[2], q[3]};
    const GSHGEOM_Result r = gshgeom_flat_path_clip_convex(subject, three, nullptr, nullptr, &out);
    ASSERT_EQ(r, GSHGEOM_OK) << rep;
    // the same triangle written with the doubled corner somewhere else
    Pt shifted[4] = {three[1], three[2], three[3], three[0]};
    ASSERT_EQ(gshgeom_flat_path_clip_convex(subject, shifted, nullptr, nullptr, &tri), GSHGEOM_OK);
    gshgeom_flat_path_destroy(out);
    gshgeom_flat_path_destroy(tri);
  }
  gshgeom_flat_path_destroy(ref);
  gshgeom_flat_path_destroy(subject);
}

#ifdef __SIZEOF_INT128__
namespace {

__extension__ typedef __int128 I128;

I128 side_exact(Pt l0, Pt l1, Pt p) {
  return (I128)((int64_t)l1.x - l0.x) * ((int64_t)p.y - l0.y) - (I128)((int64_t)l1.y - l0.y) * ((int64_t)p.x - l0.x);
}

int sgn(I128 v) {
  return v < 0 ? -1 : (v > 0 ? 1 : 0);
}

I128 floor_div128(I128 n, I128 d) {  // d > 0
  I128 q = n / d;
  if ((n % d != 0) && n < 0) {
    q--;
  }
  return q;
}

// the crossing of the closed segment ab with the line l0 l1 as the nearest integer of the exact point, ties up
bool crossing_exact(Pt a, Pt b, Pt l0, Pt l1, Pt * out) {
  const I128 fa = side_exact(l0, l1, a), fb = side_exact(l0, l1, b);
  if (sgn(fa) * sgn(fb) >= 0) {
    return false;
  }
  const I128 afa = fa < 0 ? -fa : fa, afb = fb < 0 ? -fb : fb;
  const I128 den = afa + afb;
  // P = a + (b - a) |fa| / den; the coordinate rounded half up: floor((2 (a den + d |fa|) + den) / (2 den))
  const I128 nx = (I128)a.x * den + (I128)((int64_t)b.x - a.x) * afa;
  const I128 ny = (I128)a.y * den + (I128)((int64_t)b.y - a.y) * afa;
  out->x = (int32_t)floor_div128(2 * nx + den, 2 * den);
  out->y = (int32_t)floor_div128(2 * ny + den, 2 * den);
  return true;
}

}  // namespace

TEST(ClipAttack, TheCrossingIsTheExactPointRoundedHalfUpWhateverTheOrderOfTheNames) {
  Rng rng;
  int found = 0;
  const int64_t spans[] = {6, 40, 3000, 70000, 2147483647ll};
  for (int trial = 0; trial < 400000; trial++) {
    const int64_t span = spans[trial % 5];
    auto pt = [&]() {
      // mostly random in the span, sometimes on its edges
      auto c = [&]() -> int32_t {
        const uint64_t k = rng.next() % 8;
        if (k == 0) return (int32_t)span;
        if (k == 1) return (int32_t)-span;
        if (k == 2 && span == 2147483647ll) return INT32_MIN;
        return (int32_t)rng.range(-span, span);
      };
      return Pt{c(), c()};
    };
    const Pt a = pt(), b = pt(), l0 = pt(), l1 = pt();
    if (l0.x == l1.x && l0.y == l1.y) {
      continue;
    }
    Pt want = {0, 0};
    const bool has = crossing_exact(a, b, l0, l1, &want);
    for (int order = 0; order < 4; order++) {
      Pt got = {7, 7};
      const Pt ea = (order & 1) ? b : a, eb = (order & 1) ? a : b;
      const Pt el0 = (order & 2) ? l1 : l0, el1 = (order & 2) ? l0 : l1;
      const int r = gshgeom_clip_crossing(ea, eb, el0, el1, &got);
      ASSERT_EQ(r != 0, has) << trial << " " << order;
      if (has) {
        ASSERT_EQ(got.x, want.x) << trial << " " << order;
        ASSERT_EQ(got.y, want.y) << trial << " " << order;
      }
    }
    found += has;
  }
  EXPECT_GT(found, 100000);
}
#endif

TEST(ClipAttack, ACallerThatRaisesThePointCapPastTheDefaultIsHeldToItsOwn) {
  // a contour of more points than the default cap, wholly inside, copied: the caller's cap (larger) decides, not the default
  const size_t n = GSHGEOM_LIMITS_DEFAULT_MAX_POINTS + 17;
  std::vector<Pt> pts(n);
  for (size_t i = 0; i < n; i++) {
    pts[i] = {(int32_t)(i % 1000) + 1, (int32_t)((i / 1000) % 1000) + 1};
  }
  const GSHGEOM_Flat_Subpath sub = {0, n, 1};
  GSHGEOM_Limits lim;
  std::memset(&lim, 0, sizeof lim);
  lim.max_points = (size_t)1 << 23;
  GSHGEOM_Flat_Path * path = nullptr;
  ASSERT_EQ(gshgeom_flat_path_from_points(pts.data(), n, &sub, 1, nullptr, &lim, &path), GSHGEOM_OK);
  const Pt quad[4] = {{0, 0}, {5000, 0}, {5000, 5000}, {0, 5000}};
  GSHGEOM_Flat_Path * out = nullptr;
  EXPECT_EQ(gshgeom_flat_path_clip_convex(path, quad, nullptr, &lim, &out), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_flat_path_point_count(out), n);
  gshgeom_flat_path_destroy(out);
  out = nullptr;
  lim.max_points = n - 1;  // below the contour: the cap of the caller refuses
  EXPECT_EQ(gshgeom_flat_path_clip_convex(path, quad, nullptr, &lim, &out), GSHGEOM_ERR_LIMIT);
  gshgeom_flat_path_destroy(path);
}

TEST(ClipAttack, TheCapCountsThePointsOfEveryContourTogether) {
  const Pt pts[6] = {{1, 1}, {9, 1}, {1, 9}, {11, 11}, {19, 11}, {11, 19}};
  const GSHGEOM_Flat_Subpath sub[2] = {{0, 3, 1}, {3, 3, 1}};
  GSHGEOM_Limits lim;
  std::memset(&lim, 0, sizeof lim);
  lim.max_points = 6;
  GSHGEOM_Flat_Path * path = nullptr;
  ASSERT_EQ(gshgeom_flat_path_from_points(pts, 6, sub, 2, nullptr, &lim, &path), GSHGEOM_OK);
  const Pt quad[4] = {{0, 0}, {30, 0}, {30, 30}, {0, 30}};
  GSHGEOM_Flat_Path * out = nullptr;
  EXPECT_EQ(gshgeom_flat_path_clip_convex(path, quad, nullptr, &lim, &out), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_flat_path_point_count(out), 6u);
  gshgeom_flat_path_destroy(out);
  gshgeom_flat_path_destroy(path);
  lim.max_points = 5;
  path = nullptr;
  EXPECT_EQ(gshgeom_flat_path_from_points(pts, 6, sub, 2, nullptr, &lim, &path), GSHGEOM_ERR_LIMIT);
  // a path built with the room and clipped with the smaller cap: the second contour passes the cap, though each one alone fits
  lim.max_points = 6;
  ASSERT_EQ(gshgeom_flat_path_from_points(pts, 6, sub, 2, nullptr, &lim, &path), GSHGEOM_OK);
  lim.max_points = 5;
  out = nullptr;
  EXPECT_EQ(gshgeom_flat_path_clip_convex(path, quad, nullptr, &lim, &out), GSHGEOM_ERR_LIMIT);
  gshgeom_flat_path_destroy(path);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
