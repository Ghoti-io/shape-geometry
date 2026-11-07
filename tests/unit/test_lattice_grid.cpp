/**
 * @file
 *
 * The lattice enumeration, canonical form and period over a grid of special and ordinary arguments.
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

constexpr uint64_t kGridHash = 0xdd239b9c4787ae3bull;

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



struct Rows {
  Fold * f;
  uint64_t rows = 0;
  uint64_t stop_after = 0;
};

int row_fn(void * ctx, int64_t j, int64_t lo, int64_t hi) {
  Rows * r = static_cast<Rows *>(ctx);
  r->f->word((uint64_t)j);
  r->f->word((uint64_t)lo);
  r->f->word((uint64_t)hi);
  r->rows++;
  return r->stop_after && r->rows >= r->stop_after ? 7 : 0;
}

const int64_t kIdx = (int64_t)1 << 53;

GSHGEOM_Map map_of(int k) {
  GSHGEOM_Map m;
  gshgeom_map_identity(2, &m);
  switch (k) {
    case 0: break;
    case 1: gshgeom_map_scale(2, 2.0, 2.0, 1.0, &m); break;
    case 2: gshgeom_map_rotate(GSHGEOM_ANGLE_DEGREES, 30.0, &m); break;
    case 3: gshgeom_map_rotate(GSHGEOM_ANGLE_DEGREES, 90.0, &m); break;
    case 4: gshgeom_map_scale(2, 0.0, 1.0, 1.0, &m); break;
    case 5: gshgeom_map_scale(2, 1e-30, 1e-30, 1.0, &m); break;
    case 6: gshgeom_map_scale(2, 1e30, 1e30, 1.0, &m); break;
    case 7: gshgeom_map_scale(2, -1.0, 1.0, 1.0, &m); break;
    case 8: gshgeom_map_translate(2, 1e15, -3.0, 0.0, &m); break;
    case 9: gshgeom_map_rotate(GSHGEOM_ANGLE_DEGREES, 30.0, &m); m.m[0][1] += 0.5; m.kind = 0xffffffffu; break;  // an invalid kind tag
    case 10: gshgeom_map_scale(2, 3.0, 0.5, 1.0, &m); break;
    case 11: gshgeom_map_identity(3, &m); break;
    case 12: gshgeom_map_rotate(GSHGEOM_ANGLE_DEGREES, 45.0, &m); break;
    default: gshgeom_map_scale(2, 0.25, 4.0, 1.0, &m); break;
  }
  return m;
}
const int kMaps = 14;

void lattice_grid(Fold & f) {
  const double xy[] = {0.0, -1.5, 1e15};
  const double wh[] = {1.0, 2.5, 0.0, -1.0, (double)NAN, 1e300, 1e-300};
  const double steps[] = {0.0, 1.0, 3.5, (double)NAN, (double)INFINITY, -1.0, 1e-300, 1e300, 0.1, 2.5};
  const uint32_t reps[] = {0, 1, 2};
  uint64_t s = 0x2545f4914f6cdd1dull;
  auto rnd = [&]() {
    s ^= s << 13;
    s ^= s >> 7;
    s ^= s << 17;
    return s;
  };
  const int64_t rects[][4] = {{0, 0, 3, 3}, {-5, -5, 5, 5}, {0, 0, 0, 0}, {3, 3, 0, 0}, {-100, -100, 100, 100}, {0, 0, kIdx, 1},
      {-kIdx, -kIdx, kIdx, kIdx}, {INT64_MIN, INT64_MIN, INT64_MAX, INT64_MAX}, {1000000, 1000000, 1000010, 1000010}, {0, 0, 1, 1}};
  const uint64_t caps[] = {0, 1, 10, 10000, UINT64_MAX};
  for (int it = 0; it < 2200; it++) {
    GSHGEOM_Lattice l;
    const bool valid = (rnd() % 4) != 0;
    l.x = xy[rnd() % 3];
    l.y = xy[rnd() % 3];
    l.w = valid ? (rnd() % 2 ? 1.0 : 2.5) : wh[rnd() % 7];
    l.h = valid ? (rnd() % 3 ? 1.0 : 0.75) : wh[rnd() % 7];
    l.x_step = valid ? steps[1 + rnd() % 3] * (rnd() % 2 ? 1.0 : 1.25) : steps[rnd() % 10];
    l.y_step = valid ? steps[1 + rnd() % 3] * (rnd() % 2 ? 1.0 : 0.5) : steps[rnd() % 10];
    l.repeat_x = valid ? (uint32_t)(rnd() % 2) : reps[rnd() % 3];
    l.repeat_y = valid ? (uint32_t)(rnd() % 2) : reps[rnd() % 3];
    if (valid && l.repeat_x == 0) { l.x_step = rnd() % 2 ? 0.0 : l.x_step; }
    for (int mk = 0; mk < kMaps; mk++) {
      if ((it + mk) % 3 == 2) {
        continue;
      }
      const GSHGEOM_Map m = map_of(mk);
      double out4[4] = {-7, -7, -7, -7};
      double out8[8];
      double out2[2] = {-7, -7};
      for (double & v : out8) {
        v = -7.0;
      }
      f.word((uint64_t)gshgeom_lattice_quad(&m, &l, out8));
      for (double v : out8) f.dbl(v);
      double a[2] = {-7, -7}, b[2] = {-7, -7};
      f.word((uint64_t)gshgeom_lattice_steps(&m, &l, a, b));
      f.dbl(a[0]); f.dbl(a[1]); f.dbl(b[0]); f.dbl(b[1]);
      int64_t ia[2] = {-7, -7}, ib[2] = {-7, -7};
      int integer = -7;
      f.word((uint64_t)gshgeom_lattice_integer_steps(&m, &l, ia, ib, &integer));
      f.word((uint64_t)ia[0]); f.word((uint64_t)ia[1]); f.word((uint64_t)ib[0]); f.word((uint64_t)ib[1]); f.word((uint64_t)(int64_t)integer);
      for (int64_t i : {(int64_t)0, (int64_t)1, (int64_t)-3, kIdx - 1, kIdx, -kIdx, ((int64_t)1 << 62)}) {
        for (int64_t j : {(int64_t)0, (int64_t)2, -kIdx + 1, kIdx}) {
          if (mk % 4 != 0 && (i == kIdx || j == kIdx)) {
            continue;
          }
          std::fill(out4, out4 + 4, -7.0);
          out2[0] = out2[1] = -7.0;
          f.word((uint64_t)gshgeom_lattice_cell(&l, i, j, out4));
          for (double v : out4) f.dbl(v);
          f.word((uint64_t)gshgeom_lattice_translation(&m, &l, i, j, out2));
          f.dbl(out2[0]); f.dbl(out2[1]);
          GSHGEOM_Map can;
          std::memset(&can, 0x42, sizeof can);
          int64_t sx = -7, sy = -7;
          f.word((uint64_t)gshgeom_lattice_canonical(&m, &l, i, j, &can, &sx, &sy));
          f.word((uint64_t)sx); f.word((uint64_t)sy);
          f.word(can.kind);
          for (int r = 0; r < 3; r++) { f.dbl(can.m[r][0]); f.dbl(can.m[r][1]); }
          f.dbl(can.t[0]); f.dbl(can.t[1]);
        }
      }
      // the enumeration, over the rectangles and the caps; stopped after a few rows
      const int nrect = (int)(sizeof rects / sizeof rects[0]);
      for (int r = 0; r < nrect; r++) {
        if (((it + mk + r) % 4) != 0 && r != 0 && r != 1) {
          continue;
        }
        for (uint64_t cap : caps) {
          if (((it + mk + r + (int)(cap % 7)) % 3) != 0) {
            continue;
          }
          Rows rows;
          rows.f = &f;
          rows.stop_after = (it % 5 == 0) ? 3 : 4000;
          uint64_t count = 12345;
          f.word((uint64_t)gshgeom_lattice_each(&m, &l, rects[r], cap, row_fn, &rows, &count));
          f.word(count);
          f.word(rows.rows);
        }
      }
    }
  }
  // each with no callback / no count / no rect
  {
    GSHGEOM_Lattice l = {0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 1, 1};
    GSHGEOM_Map m = map_of(0);
    const int64_t r0[4] = {0, 0, 3, 3};
    uint64_t count = 0;
    Rows rows;
    rows.f = &f;
    f.word((uint64_t)gshgeom_lattice_each(&m, &l, r0, 100, nullptr, &rows, &count));
    f.word((uint64_t)gshgeom_lattice_each(&m, &l, r0, 100, row_fn, &rows, nullptr));
    f.word((uint64_t)gshgeom_lattice_each(&m, &l, nullptr, 100, row_fn, &rows, &count));
    f.word((uint64_t)gshgeom_lattice_each(nullptr, &l, r0, 100, row_fn, &rows, &count));
    f.word((uint64_t)gshgeom_lattice_each(&m, nullptr, r0, 100, row_fn, &rows, &count));
  }
}

void period_grid(Fold & f) {
  const int64_t v[] = {0, 1, -1, 2, 3, -4, 7, 12, 100, -100, ((int64_t)1 << 31), ((int64_t)1 << 31) - 1, -((int64_t)1 << 31), ((int64_t)1 << 32),
      INT64_MAX, INT64_MIN};
  const int n = (int)(sizeof v / sizeof v[0]);
  for (int a0 = 0; a0 < n; a0++) {
    for (int a1 = 0; a1 < n; a1++) {
      for (int b0 = 0; b0 < n; b0++) {
        for (int b1 = 0; b1 < n; b1++) {
          if (((a0 * 7 + a1 * 5 + b0 * 3 + b1) % 3) == 0 || (a0 < 8 && a1 < 8 && b0 < 8 && b1 < 8)) {
            const int64_t a[2] = {v[a0], v[a1]};
            const int64_t b[2] = {v[b0], v[b1]};
            GSHGEOM_Lattice_Period pr;
            std::memset(&pr, 0x77, sizeof pr);
            const GSHGEOM_Result r = gshgeom_lattice_period(a, b, &pr);
            f.word((uint64_t)r);
            f.word((uint64_t)pr.p);
            f.word((uint64_t)pr.r);
            f.word((uint64_t)pr.s);
            if (r == GSHGEOM_OK) {
              for (int64_t x : {(int64_t)0, (int64_t)5, (int64_t)-9, (int64_t)1000003, -((int64_t)1 << 40)}) {
                for (int64_t y : {(int64_t)0, (int64_t)3, (int64_t)-17, ((int64_t)1 << 40)}) {
                  int64_t out[2] = {-7, -7};
                  gshgeom_lattice_reduce(&pr, x, y, out);
                  f.word((uint64_t)out[0]);
                  f.word((uint64_t)out[1]);
                }
              }
            }
          }
        }
      }
    }
  }
  GSHGEOM_Lattice_Period pr;
  const int64_t a[2] = {1, 0};
  f.word((uint64_t)gshgeom_lattice_period(nullptr, a, &pr));
  f.word((uint64_t)gshgeom_lattice_period(a, nullptr, &pr));
  f.word((uint64_t)gshgeom_lattice_period(a, a, nullptr));
}

}  // namespace

TEST(LatticeGrid, EveryEntryPointOverTheSpecialAndOrdinaryArgumentsIsWhatItWas) {
  Fold f;
  lattice_grid(f);
  const uint64_t after_lattice = f.h;
  period_grid(f);
  if (std::getenv("GSHGEOM_LATTICE_GRID_PRINT")) {
    std::printf("lattice-grid: lattice %016llx all %016llx\n", (unsigned long long)after_lattice, (unsigned long long)f.h);
  }
  EXPECT_EQ(f.h, kGridHash);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
