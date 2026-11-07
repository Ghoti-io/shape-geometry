/**
 * @file
 *
 * Convex clipping of flat paths (clip.h) over a grid of special and ordinary arguments.
 *
 * SPDX-License-Identifier: LGPL-3.0-only
 *
 * Copyright 2026 by Corey Pennycuff
 */

#include "../test_helpers.h"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

/**
 * Every entry point of arc.h is called on every combination of the values that decide a check (zero, one, a negative, a huge, the
 * coordinate limit and the first number past it, NaN, an infinity, a quarter turn and the full turn in each unit, a unit past the
 * last, a flag past 1) and what each returns (the result, and everything it wrote: the path's hash, the arc, the pieces) is folded
 * into one hash, which is pinned.  A check removed or weakened, a comparison changed, a constant moved or an error ignored changes
 * the answer for some combination, and the hash moves.  The values are the contract's: tests/unit/test_arc_rules.cpp holds the
 * assertion of each rule alone; this holds them all at once, and over the combinations in which one rule's refusal hides behind
 * another's.
 *
 * The hash is of the library as it stands (reviewed against arc.h for each refusal's code); a change that moves it changes what a
 * caller sees, and says so by failing here: regenerate it (GSHGEOM_ARC_GRID_PRINT=1) only for a change that is meant to.
 */

namespace {

constexpr uint64_t kGridHash = 0xcb4cb0658b59b01aull;

struct Fold {
  uint64_t h = 1469598103934665603ull;
  void word(uint64_t w) {
    for (int i = 0; i < 8; i++) {
      h ^= (w >> (8 * i)) & 0xffu;
      h *= 1099511628211ull;
    }
  }
  void dbl(double d) {
    uint64_t w;
    std::memcpy(&w, &d, 8);
    // all NaNs are one for the hash: the library does not promise a payload
    word(std::isnan(d) ? 0x7ff8000000000000ull : w);
  }
  void pt(const GSHGEOM_Point & p) {
    dbl(p.x);
    dbl(p.y);
    dbl(p.z);
  }
  void arc(const GSHGEOM_Arc & a) {
    pt(a.start);
    pt(a.end);
    pt(a.center);
    pt(a.u);
    pt(a.v);
    word((uint64_t)a.unit);
    dbl(a.t0);
    dbl(a.sweep);
    word(a.flags);
  }
};




GSHGEOM_Device_Point dp(int64_t x, int64_t y) {
  GSHGEOM_Device_Point p;
  p.x = (int32_t)x;
  p.y = (int32_t)y;
  return p;
}

static unsigned long g_stat[10];

void clip_grid(Fold & f) {
  uint64_t s = 0x853c49e6748fea9bull;
  auto rnd = [&]() {
    s ^= s << 13;
    s ^= s >> 7;
    s ^= s << 17;
    return s;
  };
  const int64_t ext[] = {4, 8, 16, 32, 64, 640, 6400, 1 << 20, 1 << 30, (int64_t)INT32_MAX - 1};
  auto coord = [&](int64_t range) { return (int64_t)(rnd() % (uint64_t)(2 * range + 1)) - range; };
  for (int it = 0; it < 30000; it++) {
    const int64_t range = ext[rnd() % 10];
    // a path: 1 to 3 subpaths of 2 to 9 points
    std::vector<GSHGEOM_Device_Point> pts;
    std::vector<GSHGEOM_Flat_Subpath> subs;
    const int ns = 1 + (int)(rnd() % 3);
    for (int k = 0; k < ns; k++) {
      GSHGEOM_Flat_Subpath sp;
      std::memset(&sp, 0, sizeof sp);
      sp.first = pts.size();
      sp.count = 2 + rnd() % 8;
      sp.closed = (int)(rnd() % 2);
      for (size_t q = 0; q < sp.count; q++) {
        if (rnd() % 6 == 0 && !pts.empty()) {
          pts.push_back(pts.back());  // a repeat
        } else {
          pts.push_back(dp(coord(range), coord(range)));
        }
      }
      subs.push_back(sp);
    }
    GSHGEOM_Flat_Path * fp = nullptr;
    if (gshgeom_flat_path_from_points(pts.data(), pts.size(), subs.data(), subs.size(), nullptr, nullptr, &fp) != GSHGEOM_OK) {
      continue;
    }
    // quads: a rectangle, a rotated square, a random quadrilateral (maybe not convex), degenerate ones
    GSHGEOM_Device_Point quads[8][4];
    const int64_t a = coord(range / 2 + 64), b = coord(range / 2 + 64), w = 64 + (int64_t)(rnd() % (uint64_t)(range + 64)), h = 64 + (int64_t)(rnd() % (uint64_t)(range + 64));
    quads[0][0] = dp(a, b); quads[0][1] = dp(a + w, b); quads[0][2] = dp(a + w, b + h); quads[0][3] = dp(a, b + h);
    quads[1][0] = dp(a, b + h); quads[1][1] = dp(a + w, b + h); quads[1][2] = dp(a + w, b); quads[1][3] = dp(a, b);  // reversed
    quads[2][0] = dp(a + w / 2, b); quads[2][1] = dp(a + w, b + h / 2); quads[2][2] = dp(a + w / 2, b + h); quads[2][3] = dp(a, b + h / 2);
    for (int c = 0; c < 4; c++) quads[3][c] = dp(coord(range), coord(range));
    quads[4][0] = dp(a, b); quads[4][1] = dp(a, b); quads[4][2] = dp(a + w, b); quads[4][3] = dp(a + w, b + h);  // a repeat
    quads[5][0] = dp(a, b); quads[5][1] = dp(a + w, b); quads[5][2] = dp(a + 2 * w, b); quads[5][3] = dp(a + w, b + h);  // collinear corners
    quads[6][0] = dp(a, b); quads[6][1] = dp(a + w, b + h); quads[6][2] = dp(a + w, b); quads[6][3] = dp(a, b + h);  // bow tie
    quads[7][0] = dp(a, b); quads[7][1] = dp(a + w, b); quads[7][2] = dp(a + w, b); quads[7][3] = dp(a, b);  // no area
    const size_t caps[] = {0, 3, 8, 40};
    for (int qi = 0; qi < 8; qi++) {
      if (((it + qi) % 2) && qi > 2) {
        continue;
      }
      const size_t cap = caps[(it + qi) % 4];
      GSHGEOM_Limits lim;
      gshgeom_limits_default(&lim);
      lim.max_points = cap;
      GSHGEOM_Flat_Path * res = reinterpret_cast<GSHGEOM_Flat_Path *>(0x1);
      const GSHGEOM_Result r = gshgeom_flat_path_clip_convex(fp, quads[qi], nullptr, cap ? &lim : nullptr, &res);
      f.word((uint64_t)r);
      g_stat[(int)r < 8 ? (int)r : 7]++;
      if (r == GSHGEOM_OK) {
        g_stat[8 + (gshgeom_flat_path_point_count(res) > 0)]++;
        f.word(gshgeom_flat_path_hash(res));
        f.word(gshgeom_flat_path_point_count(res));
        f.word(gshgeom_flat_path_subpath_count(res));
        gshgeom_flat_path_destroy(res);
      } else {
        f.word(res == reinterpret_cast<GSHGEOM_Flat_Path *>(0x1));
      }
    }
    gshgeom_flat_path_destroy(fp);
  }
  // NULL arguments
  GSHGEOM_Flat_Path * res = nullptr;
  GSHGEOM_Device_Point q[4] = {dp(0, 0), dp(64, 0), dp(64, 64), dp(0, 64)};
  f.word((uint64_t)gshgeom_flat_path_clip_convex(nullptr, q, nullptr, nullptr, &res));
}

}  // namespace

TEST(ClipGrid, EveryClipOfRandomPathsToQuadrilateralsOfEveryKindIsWhatItWas) {
  Fold f;
  clip_grid(f);
  if (std::getenv("GSHGEOM_CLIP_GRID_PRINT")) {
    std::printf("clip-grid: %016llx\n", (unsigned long long)f.h);
    for (int i = 0; i < 10; i++) std::printf("stat %d: %lu\n", i, g_stat[i]);
  }
  EXPECT_EQ(f.h, kGridHash);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
