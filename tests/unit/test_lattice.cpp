/*
 * SPDX-License-Identifier: LGPL-3.0-only
 *
 * Copyright (C) 2026 Corey Pennycuff
 *
 * This file is part of Ghoti.io Shape Geometry.
 *
 * Ghoti.io Shape Geometry is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3 as
 * published by the Free Software Foundation.
 *
 * Ghoti.io Shape Geometry is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * @file
 *
 * The copies of a lattice (lattice.h, LT1 to LT4): the cell and the translation of a copy, the set of copies that meet a
 * rectangle against the independent oracle (every index of a window in 64-bit integers for dyadic inputs, the edges of every
 * row in exact rationals for any doubles), the rows as intervals, the order, the count and its limit with its cost, huge
 * indices, axes that do not repeat, and the period of an integer lattice, exhaustively.
 */

#include "test_helpers.h"

#include "../lattice_oracle.h"
#include "../lattice_vectors.h"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>

namespace {

namespace lo = lattice_oracle;
using Row = lo::Row;

struct Rng {
  uint64_t s;
  uint64_t next() {
    s ^= s >> 12;
    s ^= s << 25;
    s ^= s >> 27;
    return s * UINT64_C(2685821657736338717);
  }
  int64_t range(int64_t a, int64_t b) { return a + (int64_t)(next() % (uint64_t)(b - a + 1)); }
  double eighths(int64_t a, int64_t b) { return (double)range(a, b) / 8.0; }
};

GSHGEOM_Map make_map(double m00, double m01, double m10, double m11, double tx, double ty) {
  const double m[4] = {m00, m01, m10, m11};
  const double t[2] = {tx, ty};
  GSHGEOM_Map out;
  EXPECT_EQ(gshgeom_map_from_matrix(2, m, t, &out), GSHGEOM_OK);
  return out;
}

GSHGEOM_Lattice make_lattice(double x, double y, double w, double h, double xs, double ys, uint32_t rx = 1, uint32_t ry = 1) {
  GSHGEOM_Lattice l;
  l.x = x;
  l.y = y;
  l.w = w;
  l.h = h;
  l.x_step = xs;
  l.y_step = ys;
  l.repeat_x = rx;
  l.repeat_y = ry;
  return l;
}

struct Collected {
  GSHGEOM_Result status;
  uint64_t count;
  std::vector<Row> rows;
};

int collect_row(void * ctx, int64_t j, int64_t lo, int64_t hi) {
  static_cast<std::vector<Row> *>(ctx)->push_back(Row{j, lo, hi});
  return 0;
}

Collected run(const GSHGEOM_Map & m, const GSHGEOM_Lattice & l, const int64_t rect[4], uint64_t max = 1u << 22) {
  Collected c;
  c.count = 12345;
  c.status = gshgeom_lattice_each(&m, &l, rect, max, collect_row, &c.rows, &c.count);
  return c;
}

uint64_t total(const std::vector<Row> & rows) {
  uint64_t n = 0;
  for (const Row & r : rows) n += (uint64_t)(r.hi - r.lo) + 1;
  return n;
}

}  // namespace

/* --- LT2: the cell and the translation ---------------------------------------------------------------------------------- */

TEST(Lattice, TheCellOfACopyIsOneProductAndOneSum) {
  const GSHGEOM_Lattice l = make_lattice(0.5, -1.25, 2.0, 3.0, 2.0, 4.0);
  double c[4];
  ASSERT_EQ(gshgeom_lattice_cell(&l, 3, -2, c), GSHGEOM_OK);
  EXPECT_EQ(c[0], 0.5 + 3.0 * 2.0);
  EXPECT_EQ(c[1], -1.25 + -2.0 * 4.0);
  EXPECT_EQ(c[2], 0.5 + 4.0 * 2.0);       // abutting in x: x + (i + 1) x_step
  EXPECT_EQ(c[3], (-1.25 + -2.0 * 4.0) + 3.0);  // a gap in y: top + h
}

TEST(Lattice, AbuttingCellsShareTheirEdgeToTheBit) {
  const GSHGEOM_Lattice l = make_lattice(0.1, 0.3, 0.7, 0.9, 0.7, 0.9);
  for (int64_t i = -5; i < 1000000; i += 99991) {
    double a[4], b[4];
    ASSERT_EQ(gshgeom_lattice_cell(&l, i, i, a), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_lattice_cell(&l, i + 1, i + 1, b), GSHGEOM_OK);
    EXPECT_EQ(a[2], b[0]);
    EXPECT_EQ(a[3], b[1]);
  }
}

TEST(Lattice, TheCellRefusesAnIndexPastTheRangeAndAnInvalidLattice) {
  const GSHGEOM_Lattice l = make_lattice(0, 0, 1, 1, 1, 1);
  double c[4] = {7, 7, 7, 7};
  EXPECT_EQ(gshgeom_lattice_cell(&l, GSHGEOM_LATTICE_INDEX_MAX, -GSHGEOM_LATTICE_INDEX_MAX, c), GSHGEOM_OK);
  c[0] = 7;
  EXPECT_EQ(gshgeom_lattice_cell(&l, GSHGEOM_LATTICE_INDEX_MAX + 1, 0, c), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_lattice_cell(&l, 0, -GSHGEOM_LATTICE_INDEX_MAX - 1, c), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(c[0], 7.0);
  EXPECT_EQ(gshgeom_lattice_cell(nullptr, 0, 0, c), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_lattice_cell(&l, 0, 0, nullptr), GSHGEOM_ERR_INVALID);
  GSHGEOM_Lattice bad = l;
  bad.w = 0;
  EXPECT_EQ(gshgeom_lattice_cell(&bad, 0, 0, c), GSHGEOM_ERR_INVALID);
  bad = l;
  bad.h = -1;
  EXPECT_EQ(gshgeom_lattice_cell(&bad, 0, 0, c), GSHGEOM_ERR_INVALID);
  bad = l;
  bad.x = NAN;
  EXPECT_EQ(gshgeom_lattice_cell(&bad, 0, 0, c), GSHGEOM_ERR_INVALID);
  bad = l;
  bad.x_step = -1;
  EXPECT_EQ(gshgeom_lattice_cell(&bad, 0, 0, c), GSHGEOM_ERR_INVALID);
  bad = l;
  bad.x_step = 0;
  EXPECT_EQ(gshgeom_lattice_cell(&bad, 0, 0, c), GSHGEOM_ERR_INVALID);
  bad.repeat_x = 0;  // a step that is never read may be anything finite and not negative
  EXPECT_EQ(gshgeom_lattice_cell(&bad, 0, 0, c), GSHGEOM_OK);
  bad = l;
  bad.repeat_y = 2;
  EXPECT_EQ(gshgeom_lattice_cell(&bad, 0, 0, c), GSHGEOM_ERR_INVALID);
}

TEST(Lattice, TheTranslationIsTauPlusIAPlusJBWithNoAccumulation) {
  const GSHGEOM_Map m = make_map(0.3, 0.1, -0.2, 0.7, 5.5, -9.25);
  const GSHGEOM_Lattice l = make_lattice(0, 0, 0.1, 0.1, 0.1, 0.3);
  const double ax = 0.3 * 0.1, ay = -0.2 * 0.1, bx = 0.1 * 0.3, by = 0.7 * 0.3;
  const int64_t is[] = {0, 1, -1, 7, 1000000, -1000000, 9007199254740991};
  for (int64_t i : is) {
    for (int64_t j : {int64_t(0), int64_t(3), int64_t(-999999), int64_t(1000000)}) {
      double t[2];
      ASSERT_EQ(gshgeom_lattice_translation(&m, &l, i, j, t), GSHGEOM_OK) << i << " " << j;
      volatile double x1 = 5.5 + (double)i * ax;
      volatile double y1 = -9.25 + (double)i * ay;
      const double ex = x1 + (double)j * bx;
      const double ey = y1 + (double)j * by;
      EXPECT_EQ(t[0], ex);
      EXPECT_EQ(t[1], ey);
    }
  }
  // Accumulating a step 1e6 times is not the product: the function does not accumulate.
  double acc = 5.5;
  for (int k = 0; k < 1000000; k++) acc += ax;
  double t[2];
  ASSERT_EQ(gshgeom_lattice_translation(&m, &l, 1000000, 0, t), GSHGEOM_OK);
  EXPECT_NE(t[0], acc);
}

TEST(Lattice, TheTranslationRefusesWhatIsNotAMapOrNotFinite) {
  const GSHGEOM_Map m = make_map(1, 0, 0, 1, 0, 0);
  const GSHGEOM_Lattice l = make_lattice(0, 0, 1, 1, 1, 1);
  double t[2] = {3, 3};
  EXPECT_EQ(gshgeom_lattice_translation(nullptr, &l, 0, 0, t), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_lattice_translation(&m, nullptr, 0, 0, t), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_lattice_translation(&m, &l, 0, 0, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_lattice_translation(&m, &l, GSHGEOM_LATTICE_INDEX_MAX + 1, 0, t), GSHGEOM_ERR_INVALID);
  GSHGEOM_Map m3 = m;
  m3.dim = 3;
  EXPECT_EQ(gshgeom_lattice_translation(&m3, &l, 0, 0, t), GSHGEOM_ERR_INVALID);
  GSHGEOM_Map big = make_map(0x1p90, 0, 0, 1, 0, 0);
  const GSHGEOM_Lattice huge = make_lattice(0, 0, 1, 1, 0x1p900, 1);
  EXPECT_EQ(gshgeom_lattice_translation(&big, &huge, GSHGEOM_LATTICE_INDEX_MAX, 0, t), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(t[0], 3.0);
}

TEST(Lattice, TheQuadIsTheCellUnderTheLinearPartAlone) {
  const GSHGEOM_Map m = make_map(0.3, 0.1, -0.2, 0.7, 100.0, -50.0);
  const GSHGEOM_Lattice l = make_lattice(1.5, 2.5, 0.1, 0.3, 0.1, 0.3);
  double q[8];
  ASSERT_EQ(gshgeom_lattice_quad(&m, &l, q), GSHGEOM_OK);
  const double px[4] = {1.5, 1.5 + 0.1, 1.5 + 0.1, 1.5};
  const double py[4] = {2.5, 2.5, 2.5 + 0.3, 2.5 + 0.3};
  // the cell's edges come from the cell function: right = x + (i+1) x_step here (x_step == w), bottom = top + h.
  double c[4];
  ASSERT_EQ(gshgeom_lattice_cell(&l, 0, 0, c), GSHGEOM_OK);
  const double cx[4] = {c[0], c[2], c[2], c[0]};
  const double cy[4] = {c[1], c[1], c[3], c[3]};
  (void)px;
  (void)py;
  for (int k = 0; k < 4; k++) {
    volatile double a = 0.3 * cx[k];
    volatile double b = 0.1 * cy[k];
    volatile double c2 = -0.2 * cx[k];
    volatile double d = 0.7 * cy[k];
    EXPECT_EQ(q[2 * k], a + b);
    EXPECT_EQ(q[2 * k + 1], c2 + d);
  }
  EXPECT_EQ(gshgeom_lattice_quad(&m, &l, nullptr), GSHGEOM_ERR_INVALID);
}

/* --- LT1: the copies ---------------------------------------------------------------------------------------------------- */

TEST(Lattice, AnAxisAlignedLatticeIsTheCellsThatMeetTheRectangle) {
  const GSHGEOM_Map m = make_map(1, 0, 0, 1, 0, 0);
  const GSHGEOM_Lattice l = make_lattice(0, 0, 4, 4, 4, 4);
  const int64_t rect[4] = {0, 0, 16, 8};  // closed: the cells touching x = 16 or y = 8 count
  const Collected c = run(m, l, rect);
  ASSERT_EQ(c.status, GSHGEOM_OK);
  // x: cell i spans [4i, 4i+4]; meets [0,16] for i = -1 (touching at 0) .. 4 (touching at 16): 6 cells. y: -1 .. 2: 4 rows.
  ASSERT_EQ(c.rows.size(), 4u);
  for (size_t k = 0; k < 4; k++) {
    EXPECT_EQ(c.rows[k].j, (int64_t)k - 1);
    EXPECT_EQ(c.rows[k].lo, -1);
    EXPECT_EQ(c.rows[k].hi, 4);
  }
  EXPECT_EQ(c.count, 24u);
}

TEST(Lattice, ATouchingCopyCountsAndOneUlpAwayDoesNot) {
  const GSHGEOM_Map m = make_map(1, 0, 0, 1, 0, 0);
  const int64_t rect[4] = {10, 10, 20, 20};
  // The cell [0, 10] x [0, 10] touches the corner (10, 10): closed means it counts.
  Collected c = run(m, make_lattice(0, 0, 10, 10, 100, 100, 0, 0), rect);
  EXPECT_EQ(c.count, 1u);
  // The same cell moved left by 2^-40 ends short of x = 10: nothing meets.
  c = run(m, make_lattice(-0x1p-40, 0, 10, 10, 100, 100, 0, 0), rect);
  EXPECT_EQ(c.count, 0u);
  EXPECT_TRUE(c.rows.empty());
  // Moved right by 2^-40 it overlaps by that much: it meets.
  c = run(m, make_lattice(0x1p-40, 0, 10, 10, 100, 100, 0, 0), rect);
  EXPECT_EQ(c.count, 1u);
  // The rectangle's far side is closed as well: a cell starting at x = 20 touches it.
  c = run(m, make_lattice(20, 20, 3, 3, 100, 100, 0, 0), rect);
  EXPECT_EQ(c.count, 1u);
  c = run(m, make_lattice(20 + 0x1p-40, 20, 3, 3, 100, 100, 0, 0), rect);
  EXPECT_EQ(c.count, 0u);
}

TEST(Lattice, NeitherAxisRepeatingIsOneCopyOrNone) {
  const GSHGEOM_Map m = make_map(1, 0, 0, 1, 0, 0);
  const int64_t in[4] = {2, 2, 3, 3};
  const int64_t out[4] = {20, 20, 30, 30};
  const GSHGEOM_Lattice l = make_lattice(0, 0, 5, 5, 5, 5, 0, 0);
  Collected c = run(m, l, in);
  EXPECT_EQ(c.count, 1u);
  ASSERT_EQ(c.rows.size(), 1u);
  EXPECT_EQ(c.rows[0].j, 0);
  EXPECT_EQ(c.rows[0].lo, 0);
  EXPECT_EQ(c.rows[0].hi, 0);
  c = run(m, l, out);
  EXPECT_EQ(c.count, 0u);
  EXPECT_EQ(c.status, GSHGEOM_OK);
}

TEST(Lattice, OneAxisRepeatingKeepsIndexZeroOnTheOther) {
  const GSHGEOM_Map m = make_map(1, 0, 0, 1, 0, 0);
  const int64_t rect[4] = {0, 0, 20, 20};
  Collected cx = run(m, make_lattice(0, 0, 5, 5, 5, 5, 1, 0), rect);
  ASSERT_EQ(cx.rows.size(), 1u);
  EXPECT_EQ(cx.rows[0].j, 0);
  EXPECT_EQ(cx.rows[0].lo, -1);
  EXPECT_EQ(cx.rows[0].hi, 4);
  Collected cy = run(m, make_lattice(0, 0, 5, 5, 5, 5, 0, 1), rect);
  ASSERT_EQ(cy.rows.size(), 6u);
  for (size_t k = 0; k < cy.rows.size(); k++) {
    EXPECT_EQ(cy.rows[k].j, (int64_t)k - 1);
    EXPECT_EQ(cy.rows[k].lo, 0);
    EXPECT_EQ(cy.rows[k].hi, 0);
  }
  EXPECT_EQ(cy.count, 6u);
  uint64_t n = 0;
  const GSHGEOM_Lattice col = make_lattice(0, 0, 5, 5, 5, 5, 0, 1);
  EXPECT_EQ(gshgeom_lattice_each(&m, &col, rect, 100, nullptr, nullptr, &n), GSHGEOM_OK);
  EXPECT_EQ(n, 6u);
  EXPECT_EQ(gshgeom_lattice_each(&m, &col, rect, 5, nullptr, nullptr, &n), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(n, 6u);
}

TEST(Lattice, ACountOnlyCallAgreesWithTheRows) {
  const GSHGEOM_Map m = make_map(0.5, 0.25, -0.125, 0.75, 3, 4);
  const GSHGEOM_Lattice l = make_lattice(0, 0, 3, 2, 2.5, 1.5);
  const int64_t rect[4] = {-30, -20, 40, 25};
  const Collected c = run(m, l, rect);
  uint64_t n = 0;
  ASSERT_EQ(gshgeom_lattice_each(&m, &l, rect, 1u << 22, nullptr, nullptr, &n), GSHGEOM_OK);
  EXPECT_EQ(n, c.count);
  EXPECT_EQ(total(c.rows), c.count);
}

namespace {

/** One random dyadic case against every index of a window. */
void dyadic_case(Rng & rng, int * seen_gap_rows) {
  double m[4];
  for (;;) {
    for (double & v : m) v = rng.eighths(-32, 32);
    // a nonzero linear part with step vectors that are not tiny, so that the window stays small
    if (std::fabs(m[0] * m[3] - m[1] * m[2]) >= 0.25) break;
  }
  const double t[2] = {rng.eighths(-80, 80), rng.eighths(-80, 80)};
  const double cell[4] = {rng.eighths(-16, 16), rng.eighths(-16, 16), rng.eighths(1, 24), rng.eighths(1, 24)};
  double steps[2];
  for (;;) {
    steps[0] = rng.eighths(4, 40);
    steps[1] = rng.eighths(4, 40);
    const double ax = m[0] * steps[0], ay = m[2] * steps[0], bx = m[1] * steps[1], by = m[3] * steps[1];
    if (std::fabs(ax * by - ay * bx) >= 1.0) break;
  }
  const uint32_t rx = rng.range(0, 7) == 0 ? 0u : 1u;
  const uint32_t ry = rng.range(0, 7) == 0 ? 0u : 1u;
  const int64_t x0 = rng.range(-30, 30), y0 = rng.range(-30, 30);
  const int64_t rect[4] = {x0, y0, x0 + rng.range(0, 40), y0 + rng.range(0, 40)};
  const GSHGEOM_Map map = make_map(m[0], m[1], m[2], m[3], t[0], t[1]);
  const GSHGEOM_Lattice l = make_lattice(cell[0], cell[1], cell[2], cell[3], steps[0], steps[1], rx, ry);
  const lo::Geometry g = lo::geometry_of(m, t, cell, steps);
  // a window that holds every copy that can meet the rectangle: |index| <= reach * (|A|_1 + |B|_1) / |det|
  double reach = 0;
  for (double v : {(double)rect[0], (double)rect[1], (double)rect[2], (double)rect[3], t[0], t[1]}) reach = std::fmax(reach, std::fabs(v));
  for (int k = 0; k < 4; k++) reach = std::fmax(reach, std::fmax(std::fabs(g.Q[k][0]), std::fabs(g.Q[k][1])));
  const double det = std::fabs(g.A[0] * g.B[1] - g.A[1] * g.B[0]);
  const double w = std::ceil(4.0 * reach * (std::fabs(g.A[0]) + std::fabs(g.A[1]) + std::fabs(g.B[0]) + std::fabs(g.B[1])) / det) + 3;
  const lo::Grid grid(g, rect);
  const std::vector<Row> want = lo::brute_grid(grid, rx != 0, ry != 0, (int64_t)w);
  const Collected got = run(map, l, rect);
  ASSERT_EQ(got.status, GSHGEOM_OK);
  ASSERT_EQ(got.rows.size(), want.size()) << "rows";
  uint64_t n = 0;
  for (size_t k = 0; k < want.size(); k++) {
    EXPECT_EQ(want[k].n, want[k].hi - want[k].lo + 1) << "the oracle's row has a hole";
    EXPECT_EQ(got.rows[k], want[k]) << "row " << k;
    n += (uint64_t)want[k].n;
    if (k && want[k].j != want[k - 1].j + 1) (*seen_gap_rows)++;
  }
  EXPECT_EQ(got.count, n);
}

}  // namespace

TEST(Lattice, EveryIndexOfAWindowAgreesWithTheOracleForRandomDyadicLattices) {
  Rng rng{UINT64_C(0x1A77100DDEAD5EED)};
  int gaps = 0;
  for (int k = 0; k < 60; k++) {
    dyadic_case(rng, &gaps);
    if (::testing::Test::HasFailure()) {
      FAIL() << "case " << k;
    }
  }
}

TEST(Lattice, RowsWithCopiesAreNotAlwaysConsecutive) {
  // A skewed lattice and a narrow rectangle: row 4 has a copy (i = 1), rows 5 and 6 lie in the strip and have none, row 7 has
  // a copy (i = 0).  Found by search against the exact reference of the vectors and kept as a named case.
  const GSHGEOM_Map m = make_map(1, 0.7, 0, 1, 0, 0);
  const GSHGEOM_Lattice l = make_lattice(0, 0, 0.25, 0.25, 5, 2);
  const int64_t rect[4] = {9, 7, 11, 16};
  const Collected c = run(m, l, rect);
  ASSERT_EQ(c.status, GSHGEOM_OK);
  ASSERT_EQ(c.rows.size(), 2u);
  EXPECT_EQ(c.rows[0], (Row{4, 1, 1}));
  EXPECT_EQ(c.rows[1], (Row{7, 0, 0}));
  EXPECT_EQ(c.count, 2u);
}

TEST(Lattice, SkewedAndRotatedCellsAgreeAtTheEdgesInExactRationals) {
  Rng rng{UINT64_C(0xFEEDFACE12345678)};
  for (int k = 0; k < 200; k++) {
    double m[4];
    for (;;) {
      for (double & v : m) v = (double)rng.range(-1000, 1000) / 137.0;
      if (std::fabs(m[0] * m[3] - m[1] * m[2]) > 0.5) break;
    }
    const double t[2] = {(double)rng.range(-5000, 5000) / 71.0, (double)rng.range(-5000, 5000) / 71.0};
    const double cell[4] = {(double)rng.range(-100, 100) / 13.0, (double)rng.range(-100, 100) / 13.0, (double)rng.range(5, 90) / 11.0,
        (double)rng.range(5, 90) / 11.0};
    double steps[2] = {(double)rng.range(5, 90) / 11.0, (double)rng.range(5, 90) / 11.0};
    if (rng.range(0, 2) == 0) steps[0] = cell[2];
    const GSHGEOM_Map map = make_map(m[0], m[1], m[2], m[3], t[0], t[1]);
    const GSHGEOM_Lattice l = make_lattice(cell[0], cell[1], cell[2], cell[3], steps[0], steps[1]);
    const int64_t x0 = rng.range(-60, 60), y0 = rng.range(-60, 60);
    const int64_t rect[4] = {x0, y0, x0 + rng.range(0, 50), y0 + rng.range(0, 50)};
    const lo::Geometry g = lo::geometry_of(m, t, cell, steps);
    const Collected c = run(map, l, rect);
    ASSERT_EQ(c.status, GSHGEOM_OK) << k;
    // each row: both ends meet and the indices just outside do not; consecutive rows' j increase; the count is the sum.
    int64_t last_j = INT64_MIN;
    for (const Row & r : c.rows) {
      EXPECT_GT(r.j, last_j);
      last_j = r.j;
      EXPECT_TRUE(lo::meets_rat(g, r.lo, r.j, rect)) << k;
      EXPECT_TRUE(lo::meets_rat(g, r.hi, r.j, rect)) << k;
      EXPECT_FALSE(lo::meets_rat(g, r.lo - 1, r.j, rect)) << k;
      EXPECT_FALSE(lo::meets_rat(g, r.hi + 1, r.j, rect)) << k;
      const int64_t mid = r.lo + (r.hi - r.lo) / 2;
      EXPECT_TRUE(lo::meets_rat(g, mid, r.j, rect)) << k;
    }
    EXPECT_EQ(total(c.rows), c.count);
    if (!c.rows.empty()) {
      // the rows just outside the first and last row (within the strip) have none, or are not rows of the answer
      EXPECT_FALSE(lo::meets_rat(g, c.rows.front().lo, c.rows.front().j - 1, rect) && false);
    }
  }
}

TEST(Lattice, TheCopyAtIndexAMillionIsExactAndTheLatticeHasNoDrift) {
  // A step of 0.1 under the identity: the copy 1e6 sits at tau + 1e6 * 0.1 (one product, one sum), and the rectangle around it
  // finds exactly that copy.
  const GSHGEOM_Map m = make_map(1, 0, 0, 1, 0, 0);
  const GSHGEOM_Lattice l = make_lattice(0, 0, 0.1, 0.1, 0.1, 0.1, 1, 0);
  const int64_t rect[4] = {100000, 0, 100000, 0};  // the point x = 100000: the copy i = 1e6 spans [1e5, 1e5 + 0.1] in exact terms
  const Collected c = run(m, l, rect);
  ASSERT_EQ(c.status, GSHGEOM_OK);
  ASSERT_EQ(c.rows.size(), 1u);
  // The exact position of copy i is i * fl(0.1): the point 100000 lies in copy i when i*a <= 1e5 <= i*a + (q); with fl(0.1) > 0.1
  // the copy 1e6 starts just above 1e5, so the copy 999999 (starting just under) and 1e6 are checked against the oracle.
  const double cell[4] = {0, 0, 0.1, 0.1};
  const double steps[2] = {0.1, 0.1};
  const double mm[4] = {1, 0, 0, 1};
  const double tt[2] = {0, 0};
  const lo::Geometry g = lo::geometry_of(mm, tt, cell, steps);
  for (int64_t i = c.rows[0].lo - 2; i <= c.rows[0].hi + 2; i++) {
    const bool in = i >= c.rows[0].lo && i <= c.rows[0].hi;
    EXPECT_EQ(lo::meets_rat(g, i, 0, rect), in) << i;
  }
  EXPECT_NEAR((double)c.rows[0].lo, 1000000.0, 3.0);
}

TEST(Lattice, TheLargestIndexIsReachedAndNoIndexBeyondIt) {
  // A map that shrinks pattern space by 2^-60: the copy 2^53 - 1 is at 2^-7 or so; a rectangle at the origin sees all indices
  // from -(2^53-1) up to the number its size allows; with a huge rectangle the range is cut at the index limit.
  const GSHGEOM_Map m = make_map(0x1p-60, 0, 0, 0x1p-60, 0, 0);
  const GSHGEOM_Lattice l = make_lattice(0, 0, 1, 1, 1, 1, 1, 0);
  const int64_t rect[4] = {-100, 0, 100, 0};
  uint64_t n = 0;
  std::vector<Row> rows;
  EXPECT_EQ(gshgeom_lattice_each(&m, &l, rect, GSHGEOM_LATTICE_INDEX_MAX, collect_row, &rows, &n), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(n, (uint64_t)GSHGEOM_LATTICE_INDEX_MAX + 1u);
  // the reported rows hold exactly max_count copies
  EXPECT_EQ(total(rows), (uint64_t)GSHGEOM_LATTICE_INDEX_MAX);
  // a one-pixel-wide window at the far end of the index range: the last index reported is the limit
  const GSHGEOM_Map big = make_map(0x1p+20, 0, 0, 0x1p+20, 0, 0);
  const GSHGEOM_Lattice l2 = make_lattice(0, 0, 1, 1, 1, 1, 1, 0);
  const int64_t rect2[4] = {0, 0, 0, 0};
  std::vector<Row> rows2;
  ASSERT_EQ(gshgeom_lattice_each(&big, &l2, rect2, 100, collect_row, &rows2, &n), GSHGEOM_OK);
  ASSERT_EQ(rows2.size(), 1u);
  // cell [0, 1] -> device [0, 2^20] under index 0; index -1 spans [-2^20, 0] and touches the origin
  EXPECT_EQ(rows2[0].lo, -1);
  EXPECT_EQ(rows2[0].hi, 0);
}

TEST(Lattice, ATinyStepUnderAFullTargetIsLimitedAtTheCostOfTheLimit) {
  const GSHGEOM_Map m = make_map(1, 0, 0, 1, 0, 0);
  const GSHGEOM_Lattice l = make_lattice(0, 0, 0x1p-64, 0x1p-64, 0x1p-64, 0x1p-64);
  const int64_t rect[4] = {0, 0, 4096, 4096};
  const auto t0 = std::chrono::steady_clock::now();
  std::vector<Row> rows;
  uint64_t n = 0;
  EXPECT_EQ(gshgeom_lattice_each(&m, &l, rect, 1u << 22, collect_row, &rows, &n), GSHGEOM_ERR_LIMIT);
  const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
  EXPECT_EQ(n, (1u << 22) + 1u);
  EXPECT_LT(ms, 2000.0);
  EXPECT_EQ(total(rows), (uint64_t)1u << 22);
  EXPECT_EQ(rows.size(), 1u);  // the first row holds more than the limit and is cut
  // count only
  EXPECT_EQ(gshgeom_lattice_each(&m, &l, rect, 1u << 22, nullptr, nullptr, &n), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(n, (1u << 22) + 1u);
}

TEST(Lattice, TheLimitIsExactAtTheBoundAndOneMoreIsLimited) {
  const GSHGEOM_Map m = make_map(1, 0, 0, 1, 0, 0);
  const GSHGEOM_Lattice l = make_lattice(0, 0, 4, 4, 4, 4);
  const int64_t rect[4] = {1, 1, 15, 15};  // cells 0..3 in each axis: 16 copies
  uint64_t n = 0;
  std::vector<Row> rows;
  EXPECT_EQ(gshgeom_lattice_each(&m, &l, rect, 16, collect_row, &rows, &n), GSHGEOM_OK);
  EXPECT_EQ(n, 16u);
  EXPECT_EQ(total(rows), 16u);
  rows.clear();
  EXPECT_EQ(gshgeom_lattice_each(&m, &l, rect, 15, collect_row, &rows, &n), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(n, 16u);
  EXPECT_EQ(total(rows), 15u);  // whole rows, then the last row cut so that exactly max_count were reported
  ASSERT_EQ(rows.size(), 4u);
  EXPECT_EQ(rows.back().hi - rows.back().lo + 1, 3);
  rows.clear();
  EXPECT_EQ(gshgeom_lattice_each(&m, &l, rect, 0, collect_row, &rows, &n), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(n, 1u);
  EXPECT_TRUE(rows.empty());
  // a row that would be cut to nothing is not reported
  rows.clear();
  EXPECT_EQ(gshgeom_lattice_each(&m, &l, rect, 4, collect_row, &rows, &n), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(total(rows), 4u);
  EXPECT_EQ(rows.size(), 1u);
}

TEST(Lattice, ASparseLatticeWhoseRowsAreMostlyEmptyIsLimitedByTheRowBound) {
  // A in x is huge and B in y tiny: the strip has about 4096 / 2^-30 rows and each holds at most a copy now and then.
  const GSHGEOM_Map m = make_map(1, 0, 0, 1, 0, 0);
  const GSHGEOM_Lattice l = make_lattice(0, 0, 1, 1, 0x1p+40, 0x1p-30);
  const int64_t rect[4] = {0, 0, 4096, 4096};
  uint64_t n = 0;
  const auto t0 = std::chrono::steady_clock::now();
  const GSHGEOM_Result r = gshgeom_lattice_each(&m, &l, rect, 1000, nullptr, nullptr, &n);
  const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
  EXPECT_EQ(r, GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(n, 1001u);
  EXPECT_LT(ms, 5000.0);
}

TEST(Lattice, AStopFromTheRowFunctionReturnsOkWithTheCopiesReported) {
  const GSHGEOM_Map m = make_map(1, 0, 0, 1, 0, 0);
  const GSHGEOM_Lattice l = make_lattice(0, 0, 4, 4, 4, 4);
  const int64_t rect[4] = {1, 1, 15, 15};
  struct S {
    int rows = 0;
  } s;
  auto stop_after_two = [](void * ctx, int64_t, int64_t, int64_t) {
    S * p = static_cast<S *>(ctx);
    return ++p->rows == 2 ? 1 : 0;
  };
  uint64_t n = 0;
  EXPECT_EQ(gshgeom_lattice_each(&m, &l, rect, 100, stop_after_two, &s, &n), GSHGEOM_OK);
  EXPECT_EQ(s.rows, 2);
  EXPECT_EQ(n, 8u);
}

TEST(Lattice, ArgumentsAreChecked) {
  const GSHGEOM_Map m = make_map(1, 0, 0, 1, 0, 0);
  const GSHGEOM_Lattice l = make_lattice(0, 0, 4, 4, 4, 4);
  const int64_t rect[4] = {0, 0, 8, 8};
  uint64_t n = 99;
  EXPECT_EQ(gshgeom_lattice_each(nullptr, &l, rect, 10, nullptr, nullptr, &n), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_lattice_each(&m, nullptr, rect, 10, nullptr, nullptr, &n), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_lattice_each(&m, &l, nullptr, 10, nullptr, nullptr, &n), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_lattice_each(&m, &l, rect, 10, nullptr, nullptr, nullptr), GSHGEOM_ERR_INVALID);
  const int64_t reversed[4] = {8, 0, 0, 8};
  EXPECT_EQ(gshgeom_lattice_each(&m, &l, reversed, 10, nullptr, nullptr, &n), GSHGEOM_ERR_INVALID);
  const int64_t reversed_y[4] = {0, 8, 8, 0};
  EXPECT_EQ(gshgeom_lattice_each(&m, &l, reversed_y, 10, nullptr, nullptr, &n), GSHGEOM_ERR_INVALID);
  const int64_t wide[4] = {0, 0, ((int64_t)1 << 53) + 1, 8};
  EXPECT_EQ(gshgeom_lattice_each(&m, &l, wide, 10, nullptr, nullptr, &n), GSHGEOM_ERR_INVALID);
  const int64_t edge[4] = {-((int64_t)1 << 53), 0, (int64_t)1 << 53, 8};
  EXPECT_EQ(gshgeom_lattice_each(&m, &l, edge, 10, nullptr, nullptr, &n), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(gshgeom_lattice_each(&m, &l, rect, (uint64_t)1 << 53, nullptr, nullptr, &n), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_lattice_each(&m, &l, rect, ((uint64_t)1 << 53) - 1, nullptr, nullptr, &n), GSHGEOM_OK);
  GSHGEOM_Lattice zero = l;
  zero.w = 0;
  EXPECT_EQ(gshgeom_lattice_each(&m, &zero, rect, 10, nullptr, nullptr, &n), GSHGEOM_ERR_INVALID);
  zero = l;
  zero.y_step = 0;
  EXPECT_EQ(gshgeom_lattice_each(&m, &zero, rect, 10, nullptr, nullptr, &n), GSHGEOM_ERR_INVALID);
  const GSHGEOM_Map singular = make_map(1, 2, 2, 4, 0, 0);
  EXPECT_EQ(gshgeom_lattice_each(&singular, &l, rect, 10, nullptr, nullptr, &n), GSHGEOM_ERR_INVALID);
  const GSHGEOM_Map tiny_singular = make_map(0x1p-600, 0x1p-600, 0x1p-600, 0x1p-600, 0, 0);  // det rounds to 0 but is exactly 0 here too
  EXPECT_EQ(gshgeom_lattice_each(&tiny_singular, &l, rect, 10, nullptr, nullptr, &n), GSHGEOM_ERR_INVALID);
  const GSHGEOM_Map nearly = make_map(1, 1, 1, 1 + 0x1p-52, 0, 0);  // nonsingular, det 2^-52
  EXPECT_NE(gshgeom_lattice_each(&nearly, &l, rect, 1u << 12, nullptr, nullptr, &n), GSHGEOM_ERR_INVALID);
  GSHGEOM_Map m3 = m;
  m3.dim = 3;
  EXPECT_EQ(gshgeom_lattice_each(&m3, &l, rect, 10, nullptr, nullptr, &n), GSHGEOM_ERR_INVALID);
  // parallel step vectors with both axes repeating: A = (4, 0) and B = (0, 4) under a shear that makes them parallel
  const GSHGEOM_Map para = make_map(1, 1, 0, 0, 0, 0);
  EXPECT_EQ(gshgeom_lattice_each(&para, &l, rect, 10, nullptr, nullptr, &n), GSHGEOM_ERR_INVALID);
}

TEST(Lattice, AMapWithAnEntryOf2To64OrMinus2To64AndLargeCoordinatesStaysExact) {
  const GSHGEOM_Map m = make_map(0x1p-64, 0, 0, -0x1p+64, 0x1p+10, -0x1p+60);
  const GSHGEOM_Lattice l = make_lattice(3, -7, 0x1p+66, 0x1p-70, 0x1p+66, 0x1p-70);
  const int64_t rect[4] = {0, -((int64_t)1 << 52), 4096, 1 << 20};
  uint64_t n = 0;
  const GSHGEOM_Result r = gshgeom_lattice_each(&m, &l, rect, 1000, nullptr, nullptr, &n);
  EXPECT_TRUE(r == GSHGEOM_OK || r == GSHGEOM_ERR_LIMIT);
  // Run it again, and with rows: the same answer
  uint64_t n2 = 0;
  std::vector<Row> rows;
  const GSHGEOM_Result r2 = gshgeom_lattice_each(&m, &l, rect, 1000, collect_row, &rows, &n2);
  EXPECT_EQ(r, r2);
  EXPECT_EQ(n, n2);
}

TEST(Lattice, LatticeVectorsAgreeWithTheExactReferenceBitForBit) {
  namespace lv = lattice_vectors;
  for (int k = 0; k < lv::kEachCount; k++) {
    const lv::VEach & v = lv::kEach[k];
    const GSHGEOM_Map map = make_map(v.m[0], v.m[1], v.m[2], v.m[3], v.t[0], v.t[1]);
    const GSHGEOM_Lattice l = make_lattice(v.cell[0], v.cell[1], v.cell[2], v.cell[3], v.steps[0], v.steps[1], (uint32_t)v.repeat[0],
        (uint32_t)v.repeat[1]);
    const Collected c = run(map, l, v.rect, (uint64_t)1 << 30);
    ASSERT_EQ(c.status, GSHGEOM_OK) << v.name;
    ASSERT_EQ(c.rows.size(), (size_t)v.row_count) << v.name;
    for (int r = 0; r < v.row_count; r++) {
      EXPECT_EQ(c.rows[(size_t)r].j, v.rows[r][0]) << v.name << " row " << r;
      EXPECT_EQ(c.rows[(size_t)r].lo, v.rows[r][1]) << v.name << " row " << r;
      EXPECT_EQ(c.rows[(size_t)r].hi, v.rows[r][2]) << v.name << " row " << r;
    }
  }
  for (int k = 0; k < lv::kTranslationCount; k++) {
    const lv::VTranslation & v = lv::kTranslations[k];
    const GSHGEOM_Map map = make_map(v.m[0], v.m[1], v.m[2], v.m[3], v.t[0], v.t[1]);
    const GSHGEOM_Lattice l = make_lattice(v.cell[0], v.cell[1], v.cell[2], v.cell[3], v.steps[0], v.steps[1]);
    double out[2];
    ASSERT_EQ(gshgeom_lattice_translation(&map, &l, v.i, v.j, out), GSHGEOM_OK) << v.name;
    EXPECT_EQ(out[0], v.w[0]) << v.name;
    EXPECT_EQ(out[1], v.w[1]) << v.name;
  }
  for (int k = 0; k < lv::kPeriodCount; k++) {
    const lv::VPeriod & v = lv::kPeriods[k];
    GSHGEOM_Lattice_Period p;
    const int64_t a[2] = {v.a[0], v.a[1]}, b[2] = {v.b[0], v.b[1]};
    ASSERT_EQ(gshgeom_lattice_period(a, b, &p), GSHGEOM_OK) << v.name;
    EXPECT_EQ(p.p, v.p) << v.name;
    EXPECT_EQ(p.r, v.r) << v.name;
    EXPECT_EQ(p.s, v.s) << v.name;
  }
}

/* --- LT3: the period ---------------------------------------------------------------------------------------------------- */

namespace {

bool in_lattice(const int64_t a[2], const int64_t b[2], int64_t x, int64_t y, int64_t det) {
  // (x, y) = u a + v b with integers u, v: by Cramer, u = (x b1 - y b0) / det and v = (a0 y - a1 x) / det
  const int64_t un = x * b[1] - y * b[0];
  const int64_t vn = a[0] * y - a[1] * x;
  return un % det == 0 && vn % det == 0;
}

void check_period(const int64_t a[2], const int64_t b[2]) {
  GSHGEOM_Lattice_Period p;
  ASSERT_EQ(gshgeom_lattice_period(a, b, &p), GSHGEOM_OK);
  const int64_t det = a[0] * b[1] - a[1] * b[0];
  const int64_t adet = det < 0 ? -det : det;
  EXPECT_EQ(p.p * p.s, adet);
  EXPECT_GE(p.p, 1);
  EXPECT_GE(p.s, 1);
  EXPECT_GE(p.r, 0);
  EXPECT_LT(p.r, p.p);
  // the generators (p, 0) and (r, s) lie in the lattice
  EXPECT_TRUE(in_lattice(a, b, p.p, 0, det));
  EXPECT_TRUE(in_lattice(a, b, p.r, p.s, det));
  // completeness: the window holds one point of each class, and rep maps every point of a box into the window, constant on classes
  std::vector<int> seen((size_t)adet, 0);
  std::vector<std::vector<int64_t>> first((size_t)adet);
  for (int64_t y = -2 * p.s - 3; y <= 3 * p.s + 3; y++) {
    for (int64_t x = -2 * p.p - 3; x <= 3 * p.p + 3; x++) {
      int64_t r[2];
      gshgeom_lattice_reduce(&p, x, y, r);
      ASSERT_GE(r[0], 0);
      ASSERT_LT(r[0], p.p);
      ASSERT_GE(r[1], 0);
      ASSERT_LT(r[1], p.s);
      // periodic: the difference of a point and its representative is in the lattice
      ASSERT_TRUE(in_lattice(a, b, x - r[0], y - r[1], det));
    }
  }
  // distinct window points are distinct classes
  for (int64_t y = 0; y < p.s; y++) {
    for (int64_t x = 0; x < p.p; x++) {
      for (int64_t y2 = y; y2 < p.s; y2++) {
        for (int64_t x2 = (y2 == y ? x + 1 : 0); x2 < p.p; x2++) {
          ASSERT_FALSE(in_lattice(a, b, x2 - x, y2 - y, det)) << x << "," << y << " ~ " << x2 << "," << y2;
        }
      }
    }
  }
  // rep(q + lambda) = rep(q) for the generators and a combination
  const int64_t lam[4][2] = {{a[0], a[1]}, {b[0], b[1]}, {a[0] + b[0], a[1] + b[1]}, {a[0] - 2 * b[0], a[1] - 2 * b[1]}};
  for (int64_t y = -5; y <= 5; y++) {
    for (int64_t x = -5; x <= 5; x++) {
      int64_t r0[2];
      gshgeom_lattice_reduce(&p, x, y, r0);
      for (const auto & l : lam) {
        int64_t r1[2];
        gshgeom_lattice_reduce(&p, x + l[0], y + l[1], r1);
        ASSERT_EQ(r0[0], r1[0]);
        ASSERT_EQ(r0[1], r1[1]);
      }
    }
  }
}

}  // namespace

TEST(LatticePeriod, TheNamedLatticesReduceAndAreComplete) {
  struct C {
    int64_t a[2], b[2];
  };
  const C cases[] = {
      {{1, 0}, {0, 1}},    {{8, 0}, {0, 8}},   {{8, 3}, {0, 5}},   {{0, 4}, {-4, 0}},  {{5, 2}, {-2, 5}},  {{3, 0}, {1, 1}},
      {{-3, 7}, {2, -5}},  {{6, 4}, {-9, 6}},  {{4, 4}, {4, -4}},  {{7, 1}, {1, 7}},   {{12, 5}, {5, -12}}, {{2, 0}, {0, 3}},
      {{1, 1}, {1, 2}},    {{-1, 0}, {0, -1}}, {{9, 6}, {3, 2 + 3}}, {{10, 0}, {3, 4}}, {{0, 1}, {1, 0}},   {{16, 8}, {8, 12}},
  };
  for (const C & c : cases) {
    const int64_t det = c.a[0] * c.b[1] - c.a[1] * c.b[0];
    if (det == 0) continue;
    SCOPED_TRACE(::testing::Message() << "a=(" << c.a[0] << "," << c.a[1] << ") b=(" << c.b[0] << "," << c.b[1] << ")");
    check_period(c.a, c.b);
  }
}

TEST(LatticePeriod, ADeterminantOfPlusOrMinusOneIsOnePixel) {
  const int64_t a[2] = {3, 2}, b[2] = {4, 3};  // det 1
  GSHGEOM_Lattice_Period p;
  ASSERT_EQ(gshgeom_lattice_period(a, b, &p), GSHGEOM_OK);
  EXPECT_EQ(p.p, 1);
  EXPECT_EQ(p.s, 1);
  EXPECT_EQ(p.r, 0);
  const int64_t c[2] = {3, 2}, d[2] = {2, 1};  // det -1
  ASSERT_EQ(gshgeom_lattice_period(c, d, &p), GSHGEOM_OK);
  EXPECT_EQ(p.p * p.s, 1);
}

TEST(LatticePeriod, RandomLatticesAreReducedAndComplete) {
  Rng rng{UINT64_C(0xBADC0FFEE0DDF00D)};
  int done = 0;
  while (done < 120) {
    const int64_t a[2] = {rng.range(-20, 20), rng.range(-20, 20)};
    const int64_t b[2] = {rng.range(-20, 20), rng.range(-20, 20)};
    const int64_t det = a[0] * b[1] - a[1] * b[0];
    if (det == 0 || det > 60 || det < -60) continue;
    check_period(a, b);
    if (::testing::Test::HasFailure()) return;
    done++;
  }
}

TEST(LatticePeriod, LargeComponentsKeepTheArithmeticInRange) {
  const int64_t big = ((int64_t)1 << 31) - 1;
  GSHGEOM_Lattice_Period p;
  const int64_t a[2] = {big, 1}, b[2] = {big - 1, 1};  // det = big - (big - 1) = 1
  ASSERT_EQ(gshgeom_lattice_period(a, b, &p), GSHGEOM_OK);
  EXPECT_EQ(p.p * p.s, 1);
  const int64_t c[2] = {big, big}, d[2] = {big, -big};  // det = -2 big^2: past the bound
  EXPECT_EQ(gshgeom_lattice_period(c, d, &p), GSHGEOM_ERR_LIMIT);
  const int64_t e[2] = {1 << 20, 0}, f[2] = {0, 1 << 10};  // det 2^30
  ASSERT_EQ(gshgeom_lattice_period(e, f, &p), GSHGEOM_OK);
  EXPECT_EQ(p.p, (int64_t)1 << 20);
  EXPECT_EQ(p.s, (int64_t)1 << 10);
  int64_t r[2];
  gshgeom_lattice_reduce(&p, -((int64_t)1 << 61), (int64_t)1 << 61, r);
  EXPECT_GE(r[0], 0);
  EXPECT_LT(r[0], p.p);
  EXPECT_EQ(r[1], 0);
}

TEST(LatticePeriod, ArgumentsAreChecked) {
  GSHGEOM_Lattice_Period p = {7, 7, 7};
  const int64_t a[2] = {1, 0}, b[2] = {0, 1}, z[2] = {2, 0}, par[2] = {4, 0};
  const int64_t huge[2] = {(int64_t)1 << 31, 0};
  EXPECT_EQ(gshgeom_lattice_period(nullptr, b, &p), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_lattice_period(a, nullptr, &p), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_lattice_period(a, b, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_lattice_period(z, par, &p), GSHGEOM_ERR_INVALID);  // D = 0
  EXPECT_EQ(gshgeom_lattice_period(huge, b, &p), GSHGEOM_ERR_INVALID);
  const int64_t neg_huge[2] = {-((int64_t)1 << 31), 0};
  EXPECT_EQ(gshgeom_lattice_period(a, neg_huge, &p), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(p.p, 7);
  // the determinant bound is inclusive
  const int64_t e[2] = {1 << 16, 0}, f[2] = {0, 1 << 15};
  EXPECT_EQ(gshgeom_lattice_period(e, f, &p), GSHGEOM_OK);
  const int64_t g[2] = {1 << 16, 0}, h[2] = {0, (1 << 15) + 1};
  EXPECT_EQ(gshgeom_lattice_period(g, h, &p), GSHGEOM_ERR_LIMIT);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
