/**
 * @file
 *
 * The copies of a lattice (lattice.h, LT1) against a brute-force oracle written from the rule's own statement, and the
 * edges of the index range and of the canonical position that the exhaustive lattice grid did not separate from variants.
 *
 * The oracle takes integer parameters (a lattice, a map and a rectangle of small integers), so that every position is an
 * integer and the test of whether the hull of a copy's quadrilateral meets the closed rectangle is exact in 64-bit
 * arithmetic.  It then lists the copies of a window of indices that holds every copy that can meet the rectangle and
 * compares the rows with those the library reports.  The search the library runs (an estimate in doubles, fixed exactly
 * by galloping and bisection) shares no code with it.
 *
 * Copyright 2026 by Corey Pennycuff
 */

#include "../test_helpers.h"

#include <ghoti.io/shape-geometry/lattice.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <set>
#include <utility>
#include <vector>

namespace {

struct Rng {
  uint64_t s = 0x9e3779b97f4a7c15ull;
  uint64_t next() {
    s ^= s << 13;
    s ^= s >> 7;
    s ^= s << 17;
    return s;
  }
  int64_t range(int64_t lo, int64_t hi) { return lo + (int64_t)(next() % (uint64_t)(hi - lo + 1)); }
};

typedef std::pair<int64_t, int64_t> IJ;

int64_t iabs(int64_t v) {
  return v < 0 ? -v : v;
}

struct Pt {
  int64_t x, y;
};

int64_t orient(Pt a, Pt b, Pt c) {
  const int64_t v = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
  return v < 0 ? -1 : (v > 0 ? 1 : 0);
}

bool on_segment(Pt a, Pt b, Pt p) {
  return orient(a, b, p) == 0 && std::min(a.x, b.x) <= p.x && p.x <= std::max(a.x, b.x) && std::min(a.y, b.y) <= p.y &&
      p.y <= std::max(a.y, b.y);
}

// closed segments ab and cd share a point
bool segments_meet(Pt a, Pt b, Pt c, Pt d) {
  const int64_t o1 = orient(a, b, c), o2 = orient(a, b, d), o3 = orient(c, d, a), o4 = orient(c, d, b);
  if (o1 * o2 < 0 && o3 * o4 < 0) {
    return true;
  }
  return on_segment(a, b, c) || on_segment(a, b, d) || on_segment(c, d, a) || on_segment(c, d, b);
}

bool in_rect(Pt p, const int64_t r[4]) {
  return r[0] <= p.x && p.x <= r[2] && r[1] <= p.y && p.y <= r[3];
}

// p in the closed triangle abc, of nonzero area
bool in_triangle(Pt a, Pt b, Pt c, Pt p) {
  if (orient(a, b, c) == 0) {
    return false;
  }
  const int64_t s1 = orient(a, b, p), s2 = orient(b, c, p), s3 = orient(c, a, p);
  const bool neg = s1 < 0 || s2 < 0 || s3 < 0;
  const bool pos = s1 > 0 || s2 > 0 || s3 > 0;
  return !(neg && pos);
}

// the convex hull of the four points meets the closed rectangle
bool hull_meets(const Pt q[4], const int64_t r[4]) {
  for (int k = 0; k < 4; k++) {
    if (in_rect(q[k], r)) {
      return true;
    }
  }
  const Pt rc[4] = {{r[0], r[1]}, {r[2], r[1]}, {r[2], r[3]}, {r[0], r[3]}};
  for (int a = 0; a < 4; a++) {
    for (int b = a + 1; b < 4; b++) {
      for (int e = 0; e < 4; e++) {
        if (segments_meet(q[a], q[b], rc[e], rc[(e + 1) % 4])) {
          return true;
        }
      }
    }
  }
  for (int e = 0; e < 4; e++) {
    for (int skip = 0; skip < 4; skip++) {
      Pt t[3];
      int n = 0;
      for (int k = 0; k < 4; k++) {
        if (k != skip) {
          t[n++] = q[k];
        }
      }
      if (in_triangle(t[0], t[1], t[2], rc[e])) {
        return true;
      }
    }
  }
  return false;
}

struct Collect {
  std::set<IJ> got;
  std::vector<std::pair<int64_t, std::pair<int64_t, int64_t>>> rows;
};

int collect_row(void * ctx, int64_t j, int64_t lo, int64_t hi) {
  Collect * c = static_cast<Collect *>(ctx);
  c->rows.push_back({j, {lo, hi}});
  for (int64_t i = lo; i <= hi; i++) {
    c->got.insert({i, j});
  }
  return 0;
}

}  // namespace

TEST(LatticeOracle, TheRowsAreTheCopiesWhoseHullMeetsTheClosedRectangle) {
  Rng rng;
  int compared = 0;
  int nonempty = 0;
  for (int trial = 0; trial < 6000; trial++) {
    GSHGEOM_Lattice l;
    std::memset(&l, 0, sizeof l);
    l.x = (double)rng.range(-6, 6);
    l.y = (double)rng.range(-6, 6);
    l.w = (double)rng.range(1, 5);
    l.h = (double)rng.range(1, 5);
    l.x_step = (double)rng.range(1, 6);
    l.y_step = (double)rng.range(1, 6);
    const int mode = (int)(rng.next() % 4);
    l.repeat_x = (mode & 1) ? 1u : 0u;
    l.repeat_y = (mode & 2) ? 1u : 0u;
    if (rng.next() % 3 == 0) {
      l.x_step = l.w;  // abutting cells
    }
    if (rng.next() % 3 == 0) {
      l.y_step = l.h;
    }
    GSHGEOM_Map m;
    ASSERT_EQ(gshgeom_map_identity(2, &m), GSHGEOM_OK);
    int64_t mm[2][2];
    do {
      for (auto & row : mm) {
        for (auto & v : row) {
          v = rng.range(-3, 3);
        }
      }
    } while (mm[0][0] * mm[1][1] - mm[0][1] * mm[1][0] == 0);
    for (int a = 0; a < 2; a++) {
      for (int b = 0; b < 2; b++) {
        m.m[a][b] = (double)mm[a][b];
      }
    }
    const int64_t tau[2] = {rng.range(-8, 8), rng.range(-8, 8)};
    m.t[0] = (double)tau[0];
    m.t[1] = (double)tau[1];
    int64_t rect[4];
    rect[0] = rng.range(-14, 10);
    rect[1] = rng.range(-14, 10);
    rect[2] = rect[0] + rng.range(0, 14);
    rect[3] = rect[1] + rng.range(0, 14);
    const int64_t xs = (int64_t)l.x_step, ys = (int64_t)l.y_step;
    const int64_t A[2] = {mm[0][0] * xs, mm[1][0] * xs};
    const int64_t B[2] = {mm[0][1] * ys, mm[1][1] * ys};
    if (l.repeat_x && l.repeat_y && A[0] * B[1] - A[1] * B[0] == 0) {
      continue;  // parallel steps: a refusal
    }
    // the cell of copy (0, 0) under the linear part
    const int64_t cx[4] = {(int64_t)l.x, (int64_t)l.x + (int64_t)l.w, (int64_t)l.x + (int64_t)l.w, (int64_t)l.x};
    const int64_t cy[4] = {(int64_t)l.y, (int64_t)l.y, (int64_t)l.y + (int64_t)l.h, (int64_t)l.y + (int64_t)l.h};
    Pt q[4];
    for (int k = 0; k < 4; k++) {
      q[k] = {mm[0][0] * cx[k] + mm[0][1] * cy[k], mm[1][0] * cx[k] + mm[1][1] * cy[k]};
    }
    // the window of indices: the translation of a copy that meets the rectangle is within the rectangle grown by the hull
    int64_t reach = 0;
    for (int k = 0; k < 4; k++) {
      reach = std::max({reach, iabs(q[k].x), iabs(q[k].y)});
    }
    for (int k = 0; k < 4; k++) {
      reach = std::max({reach, iabs(rect[k]) + reach, iabs(tau[0]), iabs(tau[1])});
    }
    reach += 2 * 14 + 16;
    const int64_t det = A[0] * B[1] - A[1] * B[0];
    int64_t wi = 0, wj = 0;
    if (l.repeat_x) {
      wi = det != 0 ? (reach * (iabs(B[0]) + iabs(B[1])) * 2) / iabs(det) + 2 : 400;
      if (!l.repeat_y) {
        wi = reach / std::max<int64_t>(1, std::max(iabs(A[0]), iabs(A[1]))) + 2;
      }
    }
    if (l.repeat_y) {
      wj = det != 0 ? (reach * (iabs(A[0]) + iabs(A[1])) * 2) / iabs(det) + 2 : 400;
      if (!l.repeat_x) {
        wj = reach / std::max<int64_t>(1, std::max(iabs(B[0]), iabs(B[1]))) + 2;
      }
    }
    if ((2 * wi + 1) * (2 * wj + 1) > 60000) {
      continue;
    }
    std::set<IJ> want;
    for (int64_t j = -wj; j <= wj; j++) {
      for (int64_t i = -wi; i <= wi; i++) {
        Pt moved[4];
        for (int k = 0; k < 4; k++) {
          moved[k] = {q[k].x + tau[0] + i * A[0] + j * B[0], q[k].y + tau[1] + i * A[1] + j * B[1]};
        }
        if (hull_meets(moved, rect)) {
          want.insert({i, j});
        }
      }
    }
    Collect c;
    uint64_t count = 0;
    const GSHGEOM_Result r = gshgeom_lattice_each(&m, &l, rect, 1u << 30, collect_row, &c, &count);
    ASSERT_EQ(r, GSHGEOM_OK) << trial;
    EXPECT_EQ(count, c.got.size()) << trial;
    ASSERT_EQ(c.got.size(), want.size()) << "trial " << trial;
    ASSERT_TRUE(c.got == want) << "trial " << trial;
    // lattice order: rows increase in j and each row's copies in i
    for (size_t k = 1; k < c.rows.size(); k++) {
      ASSERT_LT(c.rows[k - 1].first, c.rows[k].first);
    }
    compared++;
    nonempty += want.empty() ? 0 : 1;
  }
  EXPECT_GT(compared, 3000);
  EXPECT_GT(nonempty, 1500);
}

TEST(LatticeOracle, TheIndexRangeIsClosedAtBothEnds) {
  // the cells abut at x_step = w = 1, so copy i spans [i, i + 1]; a rectangle at the far end of the range of indices
  GSHGEOM_Lattice l = {0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1u, 0u};
  GSHGEOM_Map m;
  ASSERT_EQ(gshgeom_map_identity(2, &m), GSHGEOM_OK);
  const int64_t big = (int64_t)1 << 53;
  const int64_t top = GSHGEOM_LATTICE_INDEX_MAX;
  {
    const int64_t rect[4] = {big - 4, 0, big, 1};
    Collect c;
    uint64_t count = 0;
    ASSERT_EQ(gshgeom_lattice_each(&m, &l, rect, 100, collect_row, &c, &count), GSHGEOM_OK);
    ASSERT_EQ(c.rows.size(), 1u);
    EXPECT_EQ(c.rows[0].second.first, big - 5);
    EXPECT_EQ(c.rows[0].second.second, top);  // the cell past the range is not a copy: the interval is clamped at INDEX_MAX
  }
  {
    const int64_t rect[4] = {-big, 0, -big + 4, 1};
    Collect c;
    uint64_t count = 0;
    ASSERT_EQ(gshgeom_lattice_each(&m, &l, rect, 100, collect_row, &c, &count), GSHGEOM_OK);
    ASSERT_EQ(c.rows.size(), 1u);
    EXPECT_EQ(c.rows[0].second.first, -top);  // clamped at -INDEX_MAX
    EXPECT_EQ(c.rows[0].second.second, -big + 4);
  }
  {
    // a rectangle that holds every index the lattice has: the whole range, counted and stopped at the limit
    const int64_t rect[4] = {-big, 0, big, 1};
    Collect c;
    uint64_t count = 0;
    EXPECT_EQ(gshgeom_lattice_each(&m, &l, rect, 7, collect_row, &c, &count), GSHGEOM_ERR_LIMIT);
    ASSERT_EQ(c.rows.size(), 1u);
    EXPECT_EQ(c.rows[0].second.first, -top);
    EXPECT_EQ(c.rows[0].second.second, -top + 6);
    EXPECT_EQ(count, 8u);
  }
}

TEST(LatticeOracle, ACopyPastTheCanonicalRangeIsALimitWhateverTheSignOfItsIndex) {
  // integer steps: the copy's whole pixels are floor(tau) + i a + j b, and a product past 2^40 is refused as it is formed, not
  // only when it is summed: a negative index whose product is past the range, with a translation that brings the sum back
  GSHGEOM_Lattice l = {0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1u, 1u};
  GSHGEOM_Map m;
  ASSERT_EQ(gshgeom_map_identity(2, &m), GSHGEOM_OK);
  m.t[0] = (double)(((int64_t)1 << 40) - 5);
  m.t[1] = 0.0;
  const int64_t over = ((int64_t)1 << 40) + 1;
  GSHGEOM_Map out;
  int64_t sx = 0, sy = 0;
  EXPECT_EQ(gshgeom_lattice_canonical(&m, &l, -over, 0, &out, &sx, &sy), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(gshgeom_lattice_canonical(&m, &l, 0, -over, &out, &sx, &sy), GSHGEOM_ERR_LIMIT);  // the second axis: b = (0, 1)
  m.t[0] = 0.0;
  m.t[1] = (double)(((int64_t)1 << 40) - 5);
  EXPECT_EQ(gshgeom_lattice_canonical(&m, &l, 0, -over, &out, &sx, &sy), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(gshgeom_lattice_canonical(&m, &l, over, 0, &out, &sx, &sy), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(gshgeom_lattice_canonical(&m, &l, -(((int64_t)1 << 40) - 1), 0, &out, &sx, &sy), GSHGEOM_OK);  // just inside the range
  EXPECT_EQ(sx, -(((int64_t)1 << 40) - 1));
}

TEST(LatticeOracle, ACellAtTheScaleOfTheSmallestDoublesIsDecidedExactly) {
  // Every coefficient of every half-plane is a sum of products of about 2^-530, so they are subnormal and carry a few bits: the
  // rounded evaluation has no error bound there and must not decide.  The copies that hold the origin (the rectangle is that
  // point) were found with an independent exact computation (rational arithmetic and a convex hull): the copies (87, j) for
  // j from -14 to -1, and no others.
  GSHGEOM_Lattice l;
  std::memset(&l, 0, sizeof l);
  l.x = -0x1.6e8e144954618p-527;
  l.y = -0x1.f77458b3bddc5p-530;
  l.w = 0x1.143cc0f3cb6e2p-533;
  l.h = 0x1.9ebfcd613f234p-528;
  l.x_step = 0x1.77219db76147p-533;
  l.y_step = 0x1.f11ed8255939ep-532;
  l.repeat_x = 1u;
  l.repeat_y = 1u;
  GSHGEOM_Map m;
  ASSERT_EQ(gshgeom_map_identity(2, &m), GSHGEOM_OK);
  m.m[0][0] = 2.0;
  m.m[0][1] = 1.0;
  m.m[1][0] = 1.0;
  m.m[1][1] = 2.0;
  m.t[0] = -0x1.a9b23a76503a4p-528;
  m.t[1] = 0x1.1e8385ce5b27p-533;
  const int64_t rect[4] = {0, 0, 0, 0};
  Collect c;
  uint64_t count = 0;
  ASSERT_EQ(gshgeom_lattice_each(&m, &l, rect, 2000, collect_row, &c, &count), GSHGEOM_OK);
  EXPECT_EQ(count, 14u);
  ASSERT_EQ(c.rows.size(), 14u);
  for (int64_t k = 0; k < 14; k++) {
    EXPECT_EQ(c.rows[(size_t)k].first, -14 + k);
    EXPECT_EQ(c.rows[(size_t)k].second.first, 87);
    EXPECT_EQ(c.rows[(size_t)k].second.second, 87);
  }
}

TEST(LatticeOracle, CellsFromTheSmallestDoublesUpAreCountedExactlyAtAPoint) {
  // Cells of side 2^e on steps of 2^e, from the smallest subnormal up to 2^-500, and a rectangle that is the origin: the four
  // cells that share the origin as a corner, at every scale.  (A rectangle of one unit holds about 2^-e cells per side, which is
  // the limit the call reports at any of these scales.)  For e below about -540 every coefficient of every half-plane underflows
  // to zero while its exact sign is not zero, so the estimate that starts the search is 0 / 0 and is not a number.
  const int64_t rect[4] = {0, 0, 0, 0};
  for (int e = -1074; e <= -500; e += 7) {
    const double t = std::ldexp(1.0, e);
    GSHGEOM_Lattice l = {0.0, 0.0, t, t, t, t, 1u, 1u};
    GSHGEOM_Map m;
    ASSERT_EQ(gshgeom_map_identity(2, &m), GSHGEOM_OK);
    Collect c;
    uint64_t count = 0;
    ASSERT_EQ(gshgeom_lattice_each(&m, &l, rect, 100, collect_row, &c, &count), GSHGEOM_OK) << e;
    EXPECT_EQ(count, 4u) << e;
    ASSERT_EQ(c.rows.size(), 2u) << e;
    EXPECT_EQ(c.rows[0].first, -1) << e;
    EXPECT_EQ(c.rows[0].second.first, -1) << e;
    EXPECT_EQ(c.rows[0].second.second, 0) << e;
    EXPECT_EQ(c.rows[1].first, 0) << e;
    EXPECT_EQ(c.rows[1].second.first, -1) << e;
    EXPECT_EQ(c.rows[1].second.second, 0) << e;
  }
  // and the one-unit rectangle is a limit at the scale where the cells are many
  GSHGEOM_Lattice l = {0.0, 0.0, 0x1p-600, 0x1p-600, 0x1p-600, 0x1p-600, 1u, 1u};
  GSHGEOM_Map m;
  ASSERT_EQ(gshgeom_map_identity(2, &m), GSHGEOM_OK);
  const int64_t one[4] = {0, 0, 1, 1};
  uint64_t count = 0;
  EXPECT_EQ(gshgeom_lattice_each(&m, &l, one, 100, nullptr, nullptr, &count), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(count, 101u);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
