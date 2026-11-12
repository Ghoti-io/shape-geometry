/**
 * @file
 *
 * The clip of a flattened path to a convex quadrilateral (clip.h, CK1 to CK6).
 *
 * What is measured, and against what:
 *
 * - CK2 against tests/clip_oracle.h, which knows nothing of how the clip works: the winding number of the input, at rational
 *   probe points at least three units from every edge of the input and of the quadrilateral (four cuts, each rounding a
 *   crossing by at most half a unit in each coordinate, move an edge by less than that), equals the result's inside the
 *   quadrilateral's interior and the result's is zero outside.  Over the handwritten shapes, every case of the generated
 *   vectors, and a seeded random sweep.  The sweep also gates CK4 and records the worst value seen.
 * - CK1, CK3, CK5 and CK6 against tests/clip_vectors.h (tools/gen-clip-vectors.py, exact rationals in Python): the points,
 *   their order and the flags, bit for bit, and the single crossing against the oracle's exact rounding.
 * - The refusals, the empty result, the limit, the bit copy, the half-unit ties and the allocation failures.
 *
 * Copyright 2026 by Corey Pennycuff
 */

#include "test_helpers.h"

#include "../alloc_tracker.h"
#include "../clip_oracle.h"
#include "../clip_vectors.h"
#include "../../src/clip/clip_internal.h"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

namespace {

using Dp = GSHGEOM_Device_Point;
namespace co = clip_oracle;

struct Contour {
  std::vector<Dp> pts;
  int closed = 1;
  bool operator==(const Contour & o) const {
    if (closed != o.closed || pts.size() != o.pts.size()) return false;
    for (size_t i = 0; i < pts.size(); i++) {
      if (pts[i].x != o.pts[i].x || pts[i].y != o.pts[i].y) return false;
    }
    return true;
  }
};
using Contours = std::vector<Contour>;

Dp D(int32_t x, int32_t y) {
  return Dp{x, y};
}

GSHGEOM_Flat_Path * make(const Contours & cs, const GSHGEOM_Allocator * a = nullptr) {
  std::vector<Dp> pts;
  std::vector<GSHGEOM_Flat_Subpath> subs;
  for (const Contour & c : cs) {
    subs.push_back(GSHGEOM_Flat_Subpath{pts.size(), c.pts.size(), c.closed});
    pts.insert(pts.end(), c.pts.begin(), c.pts.end());
  }
  GSHGEOM_Flat_Path * fp = nullptr;
  EXPECT_EQ(gshgeom_flat_path_from_points(pts.empty() ? nullptr : pts.data(), pts.size(), subs.empty() ? nullptr : subs.data(),
                subs.size(), a, nullptr, &fp),
      GSHGEOM_OK);
  return fp;
}

Contours read(const GSHGEOM_Flat_Path * fp) {
  Contours cs;
  const Dp * p = gshgeom_flat_path_points(fp);
  for (size_t s = 0; s < gshgeom_flat_path_subpath_count(fp); s++) {
    GSHGEOM_Flat_Subpath sub;
    EXPECT_EQ(gshgeom_flat_path_subpath(fp, s, &sub), GSHGEOM_OK);
    Contour c;
    c.closed = sub.closed;
    c.pts.assign(p + sub.first, p + sub.first + sub.count);
    cs.push_back(c);
  }
  return cs;
}

using Quad = std::vector<Dp>;

GSHGEOM_Result clip(const Contours & in, const Quad & q, Contours * out, GSHGEOM_Clip_Stats * stats = nullptr,
    const GSHGEOM_Limits * limits = nullptr) {
  GSHGEOM_Flat_Path * fp = make(in);
  GSHGEOM_Flat_Path * res = nullptr;
  GSHGEOM_Result r = gshgeom_clip_convex_counted(fp, q.data(), nullptr, limits, &res, stats);
  if (r == GSHGEOM_OK) {
    if (out) *out = read(res);
    gshgeom_flat_path_destroy(res);
  } else {
    EXPECT_EQ(res, nullptr);
  }
  gshgeom_flat_path_destroy(fp);
  return r;
}

co::Polygons polys(const Contours & cs) {
  co::Polygons ps;
  for (const Contour & c : cs) {
    co::Polygon p;
    for (const Dp & d : c.pts) p.push_back(co::IPt{d.x, d.y});
    ps.push_back(p);
  }
  return ps;
}

co::Polygon quad_poly(const Quad & q) {
  co::Polygon p;
  for (const Dp & d : q) p.push_back(co::IPt{d.x, d.y});
  return p;
}

/* The winding comparison of CK2 on a grid of probes.  Returns the number of probes compared and, through the references, how
 * many of them had a nonzero winding in the input (so a test that compares only zeros can be told from one that does not). */
struct ProbeStats {
  long compared = 0;
  long nonzero_in = 0;
  long nonzero_out = 0;
  long skipped_near = 0;
};

constexpr int kProbeN = 12;
constexpr int64_t kNearUnits = 3;  // four cuts at most 0.7072 each, rounded up

co::RPt probe_at(int i, int j, int64_t x0, int64_t y0, int64_t w, int64_t h) {
  // x0 + w (i + 0.37) / N, y0 + h (j + 0.61) / N, exactly.
  co::RPt p;
  p.x = co::Rat(co::Big(x0 * 100 * kProbeN + w * (100 * i + 37)), co::Big(100 * kProbeN));
  p.y = co::Rat(co::Big(y0 * 100 * kProbeN + h * (100 * j + 61)), co::Big(100 * kProbeN));
  return p;
}

void check_winding(const Contours & in, const Quad & q, const Contours & out, ProbeStats * st, const std::string & what) {
  const co::Polygons pin = polys(in);
  const co::Polygons pout = polys(out);
  const co::Polygon pq = quad_poly(q);
  co::Polygons edges = pin;
  edges.push_back(pq);
  int64_t lo_x = pq[0].x, hi_x = pq[0].x, lo_y = pq[0].y, hi_y = pq[0].y;
  for (const co::Polygon & p : edges) {
    for (const co::IPt & a : p) {
      lo_x = std::min(lo_x, a.x);
      hi_x = std::max(hi_x, a.x);
      lo_y = std::min(lo_y, a.y);
      hi_y = std::max(hi_y, a.y);
    }
  }
  lo_x -= 5;
  lo_y -= 5;
  hi_x += 5;
  hi_y += 5;
  for (int i = 0; i < kProbeN; i++) {
    for (int j = 0; j < kProbeN; j++) {
      const co::RPt p = probe_at(i, j, lo_x, lo_y, hi_x - lo_x, hi_y - lo_y);
      if (co::near_edges(edges, p, kNearUnits)) {
        st->skipped_near++;
        continue;
      }
      const int want = co::strictly_inside(pq, p) ? co::winding(pin, p) : 0;
      const int got = co::winding(pout, p);
      st->compared++;
      st->nonzero_in += want != 0;
      st->nonzero_out += got != 0;
      if (want != got) {
        ADD_FAILURE() << what << ": winding at probe (" << i << "," << j << ") is " << got << ", the exact clip's is " << want;
        return;
      }
    }
  }
}

/* --- the handwritten shapes -------------------------------------------------------------------------------------------- */

Contour rect(int32_t x0, int32_t y0, int32_t x1, int32_t y1, int closed = 1) {
  return Contour{{D(x0, y0), D(x1, y0), D(x1, y1), D(x0, y1)}, closed};
}
const Quad kCell = {D(0, 0), D(64, 0), D(64, 64), D(0, 64)};

TEST(Clip, TheVectorsAgreeBitForBit) {
  using namespace clip_vectors;
  ASSERT_GT(case_count, 250u);
  ProbeStats st;
  size_t empties = 0;
  for (size_t i = 0; i < case_count; i++) {
    const Case & c = cases[i];
    Contours in, want;
    for (size_t k = 0; k < c.in_count; k++) {
      const Sub & s = in_subs[c.in_first + k];
      Contour ct;
      ct.closed = s.closed;
      for (size_t j = 0; j < s.count; j++) ct.pts.push_back(D(in_points[s.first + j].x, in_points[s.first + j].y));
      in.push_back(ct);
    }
    for (size_t k = 0; k < c.out_count; k++) {
      const Sub & s = out_subs[c.out_first + k];
      Contour ct;
      ct.closed = s.closed;
      for (size_t j = 0; j < s.count; j++) ct.pts.push_back(D(out_points[s.first + j].x, out_points[s.first + j].y));
      want.push_back(ct);
    }
    Quad q;
    for (int k = 0; k < 4; k++) q.push_back(D(c.quad[k].x, c.quad[k].y));
    Contours got;
    ASSERT_EQ(clip(in, q, &got), GSHGEOM_OK) << c.name;
    ASSERT_EQ(got.size(), want.size()) << c.name;
    for (size_t k = 0; k < got.size(); k++) EXPECT_TRUE(got[k] == want[k]) << c.name << " contour " << k;
    empties += want.empty();
    check_winding(in, q, got, &st, c.name);
  }
  EXPECT_GT(empties, 20u) << "the empty result is exercised";
  EXPECT_GT(st.nonzero_in, 500) << "the winding comparison saw many probes inside the shapes";
  EXPECT_GT(st.nonzero_out, 400);
  std::printf("[clip] vectors: %zu cases, %ld probes compared, %ld nonzero, %ld skipped near an edge\n", case_count,
      st.compared, st.nonzero_in, st.skipped_near);
}

TEST(Clip, RectanglesAgainstASquareCellKeepTheirExactCorners) {
  Contours out;
  ASSERT_EQ(clip({rect(32, 32, 96, 96)}, kCell, &out), GSHGEOM_OK);
  ASSERT_EQ(out.size(), 1u);
  // The overlap is the square 32..64 in each axis: four corners, whatever the start.
  ASSERT_EQ(out[0].pts.size(), 4u);
  std::vector<std::pair<int, int>> got;
  for (const Dp & p : out[0].pts) got.emplace_back(p.x, p.y);
  std::sort(got.begin(), got.end());
  EXPECT_EQ(got, (std::vector<std::pair<int, int>>{{32, 32}, {32, 64}, {64, 32}, {64, 64}}));
  EXPECT_EQ(out[0].closed, 1);
}

TEST(Clip, AWindingTwoOverlapSurvivesAndTheNestedHoleStaysAHole) {
  ProbeStats st;
  {
    Contours in = {rect(-20, -20, 40, 40), rect(20, 20, 90, 90)};  // overlap 20..40 winds 2
    Contours out;
    ASSERT_EQ(clip(in, kCell, &out), GSHGEOM_OK);
    EXPECT_EQ(out.size(), 2u) << "neither contour is merged into the other";
    check_winding(in, kCell, out, &st, "overlap");
    const co::RPt p{co::Rat(30), co::Rat(30)};
    EXPECT_EQ(co::winding(polys(out), p), co::winding(polys(in), p));
    EXPECT_EQ(std::abs(co::winding(polys(out), p)), 2);
  }
  {
    // Outer clockwise as written by rect, inner reversed: a hole.  The hole (20..80) reaches past the cell on two sides.
    Contours in = {rect(-30, -30, 100, 100), Contour{{D(20, 20), D(20, 80), D(80, 80), D(80, 20)}, 1}};
    Contours out;
    ASSERT_EQ(clip(in, kCell, &out), GSHGEOM_OK);
    EXPECT_EQ(out.size(), 2u);
    check_winding(in, kCell, out, &st, "hole");
    const co::RPt hole{co::Rat(40), co::Rat(40)};
    const co::RPt ring{co::Rat(10), co::Rat(10)};
    EXPECT_EQ(co::winding(polys(out), hole), 0);
    EXPECT_NE(co::winding(polys(out), ring), 0);
  }
  EXPECT_GT(st.nonzero_in, 20);
}

TEST(Clip, TheBowtieKeepsBothLobesWithOppositeWinding) {
  Contours in = {Contour{{D(-20, -20), D(90, 60), D(90, -20), D(-20, 60)}, 1}};
  Contours out;
  ProbeStats st;
  ASSERT_EQ(clip(in, kCell, &out), GSHGEOM_OK);
  ASSERT_EQ(out.size(), 1u) << "one contour in, one contour out";
  check_winding(in, kCell, out, &st, "bowtie");
  const co::Polygons po = polys(out);
  const co::RPt left{co::Rat(10), co::Rat(20)}, right{co::Rat(60), co::Rat(20)}, above{co::Rat(35), co::Rat(5)};
  const int wl = co::winding(po, left), wr = co::winding(po, right);
  EXPECT_NE(wl, 0);
  EXPECT_NE(wr, 0);
  EXPECT_EQ(wl, -wr) << "the lobes of a bowtie wind opposite ways, and both survive";
  EXPECT_EQ(co::winding(po, above), 0) << "above the crossing the bowtie is empty";
  EXPECT_GT(st.nonzero_in, 5);
}

TEST(Clip, ATriangleFanKeepsEachTriangleAndTheirOverlap) {
  Contours in = {Contour{{D(32, 32), D(200, -40), D(200, 100)}, 1}, Contour{{D(32, 32), D(200, 100), D(-60, 200)}, 1},
      Contour{{D(32, 32), D(-60, 200), D(-100, 0)}, 1}, Contour{{D(32, 32), D(-100, 0), D(200, -40)}, 1}};
  Contours out;
  ProbeStats st;
  ASSERT_EQ(clip(in, kCell, &out), GSHGEOM_OK);
  EXPECT_EQ(out.size(), 4u);
  check_winding(in, kCell, out, &st, "fan");
  EXPECT_GT(st.nonzero_in, 5);
}

TEST(Clip, ContoursOnAndOutsideTheCellEdge) {
  Contours out;
  // On the boundary: the square of the cell itself is inside (closed) and copied bit for bit.
  Contours on = {rect(0, 0, 64, 64)};
  ASSERT_EQ(clip(on, kCell, &out), GSHGEOM_OK);
  ASSERT_EQ(out.size(), 1u);
  EXPECT_TRUE(out[0] == on[0]);
  // Along one side only, zero area: kept, as the closed quadrilateral contains it.
  Contours along = {Contour{{D(0, 0), D(30, 0), D(64, 0), D(30, 0)}, 1}};
  ASSERT_EQ(clip(along, kCell, &out), GSHGEOM_OK);
  ASSERT_EQ(out.size(), 1u);
  EXPECT_TRUE(out[0] == along[0]);
  // Outside, sharing only a side or a corner with the cell: nothing winds in the interior.
  ProbeStats st;
  for (const Contour & c : {rect(64, 0, 128, 64), rect(64, 64, 100, 100), rect(-40, 0, 0, 64), rect(100, 0, 140, 40)}) {
    ASSERT_EQ(clip({c}, kCell, &out), GSHGEOM_OK);
    check_winding({c}, kCell, out, &st, "outside");
    for (const Contour & o : out) {
      const co::RPt mid{co::Rat(32), co::Rat(32)};
      EXPECT_EQ(co::winding(polys({o}), mid), 0);
    }
  }
  // Wholly outside by a gap: no contour.
  ASSERT_EQ(clip({rect(100, 0, 140, 40)}, kCell, &out), GSHGEOM_OK);
  EXPECT_TRUE(out.empty());
  // A contour that touches the cell at one vertex on its side and goes away: winding zero inside either way.
  Contours touch = {Contour{{D(20, 20), D(64, 30), D(20, 40)}, 1}, Contour{{D(80, 10), D(64, 30), D(80, 50)}, 1}};
  ASSERT_EQ(clip(touch, kCell, &out), GSHGEOM_OK);
  check_winding(touch, kCell, out, &st, "touch");
  // A contour lying along a side and spilling out: the part on the side stays, the rest is cut off.
  Contours spill = {Contour{{D(0, -10), D(0, 0), D(0, 64), D(0, 80), D(-20, 40)}, 1}};
  ASSERT_EQ(clip(spill, kCell, &out), GSHGEOM_OK);
  check_winding(spill, kCell, out, &st, "spill");
}

TEST(Clip, ARotatedAndASkewedCellInEitherOrientation) {
  const Quad diamond = {D(64, 0), D(128, 64), D(64, 128), D(0, 64)};
  const Quad skew = {D(0, 0), D(100, 30), D(130, 90), D(30, 60)};
  const std::vector<Contours> shapes = {{rect(20, 20, 110, 110)}, {Contour{{D(10, 40), D(140, 70), D(60, 150)}, 1}},
      {rect(10, -20, 90, 70), Contour{{D(60, 100), D(140, 20), D(130, 130)}, 1}}};
  ProbeStats st;
  for (const Quad & q : {diamond, skew}) {
    Quad rev(q.rbegin(), q.rend());
    for (const Contours & in : shapes) {
      Contours a, b;
      ASSERT_EQ(clip(in, q, &a), GSHGEOM_OK);
      ASSERT_EQ(clip(in, rev, &b), GSHGEOM_OK);
      check_winding(in, q, a, &st, "forward");
      check_winding(in, rev, b, &st, "reversed");
      // The same region either way: the windings agree at every probe (the point lists need not).
      ProbeStats ignore;
      check_winding(in, rev, a, &ignore, "forward result against the reversed quad");
    }
  }
  EXPECT_GT(st.nonzero_in, 100);
}

TEST(Clip, RepeatedQuadPointsAreAcceptedAndATriangleClipsAsATriangle) {
  const std::vector<Quad> triangles = {
      {D(0, 0), D(100, 0), D(0, 100), D(0, 100)},  // the last point twice
      {D(0, 0), D(0, 0), D(100, 0), D(0, 100)},    // the first twice
      {D(0, 0), D(100, 0), D(100, 0), D(0, 100)},  // the middle twice
      {D(0, 0), D(50, 0), D(100, 0), D(0, 100)},   // a corner on a side (a straight angle)
      {D(0, 100), D(0, 100), D(0, 0), D(100, 0)}, {D(0, 0), D(0, 100), D(100, 0), D(100, 0)}};
  ProbeStats st;
  for (const Quad & q : triangles) {
    Contours in = {rect(-10, -10, 60, 60)};
    Contours out;
    ASSERT_EQ(clip(in, q, &out), GSHGEOM_OK);
    ASSERT_FALSE(out.empty());
    check_winding(in, q, out, &st, "triangle");
    // The point (45,45) is on the hypotenuse's far side (x + y = 90 < 100: inside); (70, 70) is outside the rectangle.
    EXPECT_NE(co::winding(polys(out), co::RPt{co::Rat(30), co::Rat(30)}), 0);
    EXPECT_EQ(co::winding(polys(out), co::RPt{co::Rat(55), co::Rat(55)}), 0) << "past the hypotenuse x + y = 100";
  }
  EXPECT_GT(st.nonzero_in, 30);
}

TEST(Clip, ANonConvexOrDegenerateQuadrilateralIsRefused) {
  const std::vector<std::pair<const char *, Quad>> bad = {
      {"all one point", {D(5, 5), D(5, 5), D(5, 5), D(5, 5)}},
      {"two points", {D(0, 0), D(9, 9), D(0, 0), D(9, 9)}},
      {"two points, adjacent pairs", {D(0, 0), D(0, 0), D(9, 9), D(9, 9)}},
      {"collinear", {D(0, 0), D(10, 10), D(20, 20), D(30, 30)}},
      {"collinear, out of order", {D(0, 0), D(30, 30), D(10, 10), D(20, 20)}},
      {"a bowtie", {D(0, 0), D(64, 64), D(64, 0), D(0, 64)}},
      {"a dart (one reflex corner)", {D(0, 0), D(64, 0), D(16, 16), D(0, 64)}},
      {"a reversal on a side", {D(0, 0), D(10, 0), D(5, 0), D(5, 5)}},
      {"a reversal with a repeated point", {D(0, 0), D(10, 0), D(0, 0), D(5, 5)}},
      {"three collinear and a point repeated", {D(0, 0), D(10, 0), D(20, 0), D(20, 0)}},
      {"a zig-zag", {D(0, 0), D(64, 64), D(0, 64), D(64, 0)}}};
  for (const auto & b : bad) {
    GSHGEOM_Flat_Path * fp = make({rect(0, 0, 64, 64)});
    GSHGEOM_Flat_Path * res = reinterpret_cast<GSHGEOM_Flat_Path *>(0x1);
    EXPECT_EQ(gshgeom_flat_path_clip_convex(fp, b.second.data(), nullptr, nullptr, &res), GSHGEOM_ERR_INVALID) << b.first;
    EXPECT_EQ(res, reinterpret_cast<GSHGEOM_Flat_Path *>(0x1)) << b.first << ": *out is written only on success";
    gshgeom_flat_path_destroy(fp);
  }
  // NULL arguments.
  GSHGEOM_Flat_Path * fp = make({rect(0, 0, 64, 64)});
  GSHGEOM_Flat_Path * res = nullptr;
  EXPECT_EQ(gshgeom_flat_path_clip_convex(nullptr, kCell.data(), nullptr, nullptr, &res), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_flat_path_clip_convex(fp, nullptr, nullptr, nullptr, &res), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_flat_path_clip_convex(fp, kCell.data(), nullptr, nullptr, nullptr), GSHGEOM_ERR_INVALID);
  gshgeom_flat_path_destroy(fp);
}

/**
 * The turns of the quadrilateral are decided by the filtered orientation (the exact path answers where the filter cannot certify a sign).
 * Whether a quadrilateral is accepted must be what the always-exact turns say, on corners that make the filter fail: the extremes of the
 * range, thin slivers off a long line by one unit, and corners that are collinear exactly.
 */
TEST(Clip, TheQuadrilateralVerdictIsTheExactTurnsOnCornersThatDefeatTheFilter) {
  const int32_t hi = INT32_MAX;
  const int32_t lo = INT32_MIN;
  const std::vector<Dp> corners = {D(0, 0), D(1, 0), D(0, 1), D(hi, hi), D(lo, lo), D(hi, lo), D(lo, hi), D(hi, hi - 1), D(hi - 1, hi),
      D(lo + 1, lo), D(1, hi), D(hi, 1), D(-1, lo), D(3, 5), D(hi - 2, hi - 3), D(lo / 2, hi / 2)};
  unsigned accepted = 0;
  unsigned refused = 0;
  GSHGEOM_Flat_Path * fp = make({rect(0, 0, 64, 64)});
  for (size_t a = 0; a < corners.size(); a++) {
    for (size_t b = 0; b < corners.size(); b++) {
      for (size_t c = 0; c < corners.size(); c++) {
        for (size_t d = 0; d < corners.size(); d += 3) {
          const Quad q = {corners[a], corners[b], corners[c], corners[d]};
          int sigma = 0;
          bool mixed = false;
          for (int i = 0; i < 4; i++) {
            const Dp & p0 = q[(i + 3) % 4];
            const Dp & p1 = q[i];
            const Dp & p2 = q[(i + 1) % 4];
            int turn = 0;
            ASSERT_EQ(gshgeom_orient2d_exact(GSHGEOM_Point{(double)p0.x, (double)p0.y, 0}, GSHGEOM_Point{(double)p1.x, (double)p1.y, 0},
                          GSHGEOM_Point{(double)p2.x, (double)p2.y, 0}, &turn),
                GSHGEOM_OK);
            if (turn != 0) {
              mixed = mixed || (sigma != 0 && turn != sigma);
              sigma = turn;
            }
          }
          GSHGEOM_Flat_Path * res = nullptr;
          const GSHGEOM_Result r = gshgeom_flat_path_clip_convex(fp, q.data(), nullptr, nullptr, &res);
          const bool expect_refused = mixed || sigma == 0;
          if (r == GSHGEOM_OK) {
            gshgeom_flat_path_destroy(res);
            accepted++;
          } else {
            refused++;
          }
          ASSERT_EQ(r, expect_refused ? GSHGEOM_ERR_INVALID : GSHGEOM_OK) << a << " " << b << " " << c << " " << d;
        }
      }
    }
  }
  gshgeom_flat_path_destroy(fp);
  EXPECT_GT(accepted, 0u) << "the sweep must reach the accepting side";
  EXPECT_GT(refused, 0u) << "the sweep must reach the refusing side";
}

TEST(Clip, TheQuadrilateralIsJudgedBeforeAnythingIsAllocatedOrCut) {
  tracker::Tracker t;
  GSHGEOM_Allocator a = tracker::allocator_for(&t);
  GSHGEOM_Flat_Path * fp = make({rect(0, 0, 64, 64)});
  GSHGEOM_Flat_Path * res = nullptr;
  const Quad bowtie = {D(0, 0), D(64, 64), D(64, 0), D(0, 64)};
  EXPECT_EQ(gshgeom_flat_path_clip_convex(fp, bowtie.data(), &a, nullptr, &res), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(t.calls, 0);
  gshgeom_flat_path_destroy(fp);
}

TEST(Clip, NothingInsideIsAnEmptyResultNotAnError) {
  GSHGEOM_Flat_Path * empty = make({});
  GSHGEOM_Flat_Path * res = nullptr;
  ASSERT_EQ(gshgeom_flat_path_clip_convex(empty, kCell.data(), nullptr, nullptr, &res), GSHGEOM_OK);
  ASSERT_NE(res, nullptr);
  EXPECT_EQ(gshgeom_flat_path_subpath_count(res), 0u);
  EXPECT_EQ(gshgeom_flat_path_point_count(res), 0u);
  gshgeom_flat_path_destroy(res);
  gshgeom_flat_path_destroy(empty);
  Contours out;
  ASSERT_EQ(clip({rect(200, 200, 300, 300), Contour{{D(-50, -50)}, 1}, Contour{{D(-50, 10), D(-10, 10)}, 0}}, kCell, &out),
      GSHGEOM_OK);
  EXPECT_TRUE(out.empty());
}

TEST(Clip, CrossingsReach4nOnAConstructedCaseAndNeverPassTheProvedBound) {
  // Four points, each edge crossing the cell's boundary on every cut: 16 crossings, four per input point.  The bound that
  // is proved is 8.125 per point; that four is the most there can be is not proved (the header says so), and this is the
  // case that shows it is reached.
  const Contours in = {Contour{{D(-197, 140), D(168, -40), D(-113, 85), D(208, -57)}, 1}};
  GSHGEOM_Clip_Stats st{};
  Contours out;
  ASSERT_EQ(clip(in, kCell, &out, &st), GSHGEOM_OK);
  EXPECT_EQ(st.crossings, 16u);
  EXPECT_LE(8 * st.crossings, 65u * 4u);
}

TEST(Clip, AContourWhollyInsideIsCopiedBitForBitFlagsAndTinyOnesIncluded) {
  Contours in = {Contour{{D(10, 10), D(50, 12), D(30, 50)}, 0}, Contour{{D(7, 7)}, 1}, Contour{{D(7, 7), D(9, 9)}, 0},
      Contour{{D(0, 0), D(64, 0), D(64, 64), D(0, 64)}, 1}, Contour{{D(5, 5), D(5, 5), D(5, 5)}, 1}};
  Contours out;
  GSHGEOM_Clip_Stats st{};
  ASSERT_EQ(clip(in, kCell, &out, &st), GSHGEOM_OK);
  ASSERT_EQ(out.size(), in.size());
  for (size_t i = 0; i < in.size(); i++) EXPECT_TRUE(out[i] == in[i]) << i;
  EXPECT_EQ(st.crossings, 0u);
  EXPECT_EQ(st.contours_copied, in.size());
  // Mixed with a cut one, the inside ones are still bit for bit and in order.
  Contours mixed = {in[0], rect(40, 40, 100, 100), in[1]};
  ASSERT_EQ(clip(mixed, kCell, &out), GSHGEOM_OK);
  ASSERT_EQ(out.size(), 3u);
  EXPECT_TRUE(out[0] == mixed[0]);
  EXPECT_TRUE(out[2] == mixed[2]);
  EXPECT_EQ(out[1].closed, 1) << "a cut contour is closed";
  // An open contour that is cut becomes closed; a one-point contour that is outside is gone.
  Contours cutopen = {Contour{{D(10, 10), D(90, 12), D(30, 90)}, 0}};
  ASSERT_EQ(clip(cutopen, kCell, &out), GSHGEOM_OK);
  ASSERT_EQ(out.size(), 1u);
  EXPECT_EQ(out[0].closed, 1);
}

TEST(Clip, TheOutputLimitIsExact) {
  Contours in = {rect(32, 32, 96, 96)};  // four points out
  GSHGEOM_Limits l;
  gshgeom_limits_default(&l);
  Contours out;
  l.max_points = 4;
  ASSERT_EQ(clip(in, kCell, &out, nullptr, &l), GSHGEOM_OK);
  EXPECT_EQ(out.size(), 1u);
  l.max_points = 3;
  EXPECT_EQ(clip(in, kCell, &out, nullptr, &l), GSHGEOM_ERR_LIMIT);
  // The limit is over all contours: two squares of four points need eight.
  Contours two = {rect(2, 2, 20, 20), rect(30, 30, 60, 60)};
  l.max_points = 8;
  ASSERT_EQ(clip(two, kCell, &out, nullptr, &l), GSHGEOM_OK);
  l.max_points = 7;
  EXPECT_EQ(clip(two, kCell, &out, nullptr, &l), GSHGEOM_ERR_LIMIT);
  // The input may be larger than the limit when the clip discards it.
  Contours far_away = {rect(1000, 1000, 2000, 2000), rect(3000, 0, 3100, 100)};
  l.max_points = 1;
  ASSERT_EQ(clip(far_away, kCell, &out, nullptr, &l), GSHGEOM_OK);
  EXPECT_TRUE(out.empty());
}

/* --- the crossing: CK3 -------------------------------------------------------------------------------------------------- */

TEST(ClipCrossing, OnlyAStrictCrossingIsACrossing) {
  Dp r{777, 777};
  // An end point on the line, both on one side, both on the line, and a segment that stops short of it: no crossing, *out untouched.
  EXPECT_EQ(gshgeom_clip_crossing(D(0, 0), D(5, 5), D(0, 0), D(10, 0), &r), 0);
  EXPECT_EQ(gshgeom_clip_crossing(D(5, 5), D(0, 0), D(0, 0), D(10, 0), &r), 0);
  EXPECT_EQ(gshgeom_clip_crossing(D(1, 1), D(5, 5), D(0, 0), D(10, 0), &r), 0);
  EXPECT_EQ(gshgeom_clip_crossing(D(1, 0), D(5, 0), D(0, 0), D(10, 0), &r), 0);
  EXPECT_EQ(gshgeom_clip_crossing(D(1, 1), D(2, 3), D(0, 5), D(10, 5), &r), 0) << "the line is infinite but the segment stops short: no sign change";
  EXPECT_EQ(r.x, 777);
  EXPECT_EQ(r.y, 777);
  EXPECT_EQ(gshgeom_clip_crossing(D(1, -1), D(1, 1), D(0, 0), D(10, 0), &r), 1);
  EXPECT_EQ(r.x, 1);
  EXPECT_EQ(r.y, 0);
  // The line is the whole line: the crossing may lie beyond the points that name it.
  EXPECT_EQ(gshgeom_clip_crossing(D(100, -1), D(100, 1), D(0, 0), D(10, 0), &r), 1);
  EXPECT_EQ(r.x, 100);
}

TEST(ClipCrossing, AHalfUnitTieRoundsUpAndIsTheSameEveryWayRound) {
  struct T {
    Dp a, b, l0, l1;
    int32_t x, y;
  };
  // Vertical line x = 1: (0,0)-(2,1) crosses at y = 1/2 -> 1; (0,0)-(2,-1) at y = -1/2 -> 0 (upward: toward +infinity).
  const std::vector<T> ties = {
      {D(0, 0), D(2, 1), D(1, 0), D(1, 64), 1, 1},
      {D(0, 0), D(2, -1), D(1, 0), D(1, 64), 1, 0},
      {D(0, 0), D(-2, 1), D(-1, 0), D(-1, 64), -1, 1},
      {D(0, 0), D(-2, -1), D(-1, 0), D(-1, 64), -1, 0},
      {D(0, 0), D(1, 2), D(0, 1), D(64, 1), 1, 1},   // horizontal line y = 1: x = 1/2 -> 1
      {D(0, 0), D(-1, 2), D(0, 1), D(64, 1), 0, 1},  // x = -1/2 -> 0
      {D(-1, 0), D(3, 0 + 4), D(0, 0), D(0, 64), 0, 1},  // x = 0: crosses at (0, 1)... exact, no tie
  };
  for (const T & t : ties) {
    Dp r{};
    for (int swap_edge = 0; swap_edge < 2; swap_edge++) {
      for (int swap_line = 0; swap_line < 2; swap_line++) {
        Dp r2{};
        ASSERT_EQ(gshgeom_clip_crossing(swap_edge ? t.b : t.a, swap_edge ? t.a : t.b, swap_line ? t.l1 : t.l0,
                      swap_line ? t.l0 : t.l1, &r2),
            1);
        if (swap_edge == 0 && swap_line == 0) r = r2;
        EXPECT_EQ(r2.x, r.x);
        EXPECT_EQ(r2.y, r.y);
      }
    }
    EXPECT_EQ(r.x, t.x) << t.a.x << "," << t.a.y;
    EXPECT_EQ(r.y, t.y) << t.a.x << "," << t.a.y;
  }
}

TEST(ClipCrossing, EveryTieInASmallBoxIsRoundedUpAndTheSameBothWays) {
  // All edges and lines with end points in a small box: find crossings whose exact coordinate is k + 1/2 and compare with the
  // oracle's exact rounding, in the four namings.
  long ties = 0, crossings = 0;
  std::vector<Dp> pts;
  for (int x = -2; x <= 2; x += 2) {
    for (int y = -3; y <= 3; y += 1) pts.push_back(D(x, y));
  }
  for (const Dp & a : pts) {
    for (const Dp & b : pts) {
      for (const Dp & l0 : pts) {
        for (const Dp & l1 : pts) {
          if (l0.x == l1.x && l0.y == l1.y) continue;
          Dp r;
          if (!gshgeom_clip_crossing(a, b, l0, l1, &r)) continue;
          crossings++;
          co::IPt want;
          ASSERT_TRUE(co::crossing_rounded(co::IPt{a.x, a.y}, co::IPt{b.x, b.y}, co::IPt{l0.x, l0.y}, co::IPt{l1.x, l1.y}, &want));
          ASSERT_EQ(r.x, want.x);
          ASSERT_EQ(r.y, want.y);
          Dp q;
          ASSERT_EQ(gshgeom_clip_crossing(b, a, l1, l0, &q), 1);
          ASSERT_EQ(q.x, r.x);
          ASSERT_EQ(q.y, r.y);
          // A tie: twice the exact coordinate is odd.
          const co::Rat fa = co::side_at(co::IPt{l0.x, l0.y}, co::IPt{l1.x, l1.y}, co::RPt{co::Rat(a.x), co::Rat(a.y)});
          const co::Rat fb = co::side_at(co::IPt{l0.x, l0.y}, co::IPt{l1.x, l1.y}, co::RPt{co::Rat(b.x), co::Rat(b.y)});
          const co::Rat t = fa / (fa - fb);
          for (int axis = 0; axis < 2; axis++) {
            const co::Rat v = co::Rat(axis ? a.y : a.x) + co::Rat((axis ? b.y : b.x) - (axis ? a.y : a.x)) * t;
            const co::Rat twice = v * co::Rat(2);
            // twice is an odd integer iff twice.n is an odd multiple of twice.d
            const co::Rat fl = co::Rat(co::round_half_up(twice));
            if (co::cmp(fl, twice) == 0 && (co::round_half_up(twice) % 2 != 0)) ties++;
          }
        }
      }
    }
  }
  EXPECT_GT(crossings, 10000);
  EXPECT_GT(ties, 2000) << "the box holds many exact ties";
  std::printf("[clip] ties: %ld tie coordinates among %ld crossings, all as the exact rounding\n", ties, crossings);
}

TEST(ClipCrossing, TheSideIsTheExactOrientationAndTheCrossingTheExactRoundingAcrossTheWholeInt32Range) {
  uint64_t s = UINT64_C(0xC0FFEE1234567);
  auto next = [&s]() {
    s ^= s >> 12;
    s ^= s << 25;
    s ^= s >> 27;
    return s * UINT64_C(2685821657736338717);
  };
  auto coord = [&](int mode) -> int32_t {
    const uint64_t v = next();
    switch (mode) {
      case 0: return (int32_t)(int64_t)(v % 2001) - 1000;
      case 1: return (int32_t)v;
      default: return (v & 1) ? INT32_MAX - (int32_t)(v >> 40 & 7) : INT32_MIN + (int32_t)(v >> 40 & 7);
    }
  };
  long agree = 0, crossed = 0;
  for (int i = 0; i < 6000; i++) {
    const int mode = i % 3;
    Dp a = D(coord(mode), coord(mode)), b = D(coord(mode), coord(mode)), l0 = D(coord(mode), coord(mode)),
       l1 = D(coord(mode), coord(mode));
    if (l0.x == l1.x && l0.y == l1.y) continue;
    int exact = 9;
    ASSERT_EQ(gshgeom_orient2d_exact(GSHGEOM_Point{(double)l0.x, (double)l0.y, 0}, GSHGEOM_Point{(double)l1.x, (double)l1.y, 0},
                  GSHGEOM_Point{(double)a.x, (double)a.y, 0}, &exact),
        GSHGEOM_OK);
    EXPECT_EQ(gshgeom_clip_side(l0, l1, a), exact) << i;
    agree++;
    Dp r;
    const int got = gshgeom_clip_crossing(a, b, l0, l1, &r);
    co::IPt want;
    const bool has = co::crossing_rounded(co::IPt{a.x, a.y}, co::IPt{b.x, b.y}, co::IPt{l0.x, l0.y}, co::IPt{l1.x, l1.y}, &want);
    ASSERT_EQ(got != 0, has) << i;
    if (has) {
      crossed++;
      ASSERT_EQ(r.x, want.x) << i;
      ASSERT_EQ(r.y, want.y) << i;
    }
  }
  EXPECT_GT(agree, 5000);
  EXPECT_GT(crossed, 1000);
}

/* --- allocation failure -------------------------------------------------------------------------------------------------- */

TEST(Clip, AnAllocationThatFailsAtAnyCallLeavesNothingBehind) {
  Contours in = {rect(-30, -30, 100, 100), Contour{{D(20, 20), D(20, 80), D(80, 80), D(80, 20)}, 1},
      Contour{{D(10, 10), D(50, 12), D(30, 50)}, 0}, rect(200, 0, 300, 50)};
  GSHGEOM_Flat_Path * fp = make(in);
  tracker::Tracker count;
  {
    GSHGEOM_Allocator a = tracker::allocator_for(&count);
    GSHGEOM_Flat_Path * out = nullptr;
    ASSERT_EQ(gshgeom_flat_path_clip_convex(fp, kCell.data(), &a, nullptr, &out), GSHGEOM_OK);
    gshgeom_flat_path_destroy(out);
    EXPECT_TRUE(count.live.empty());
  }
  ASSERT_GT(count.calls, 4);
  for (long k = 0; k < count.calls; k++) {
    for (int persistent = 0; persistent < 2; persistent++) {
      tracker::Tracker t;
      t.fail_from = k;
      t.persistent = persistent != 0;
      GSHGEOM_Allocator a = tracker::allocator_for(&t);
      GSHGEOM_Flat_Path * out = nullptr;
      EXPECT_EQ(gshgeom_flat_path_clip_convex(fp, kCell.data(), &a, nullptr, &out), GSHGEOM_ERR_OOM) << k;
      EXPECT_EQ(out, nullptr);
      EXPECT_TRUE(t.live.empty()) << "k=" << k << " leaked";
    }
  }
  gshgeom_flat_path_destroy(fp);
}

TEST(Clip, TheSameInputGivesTheSameBitsAndTheInputIsUntouched) {
  Contours in = {Contour{{D(-20, -20), D(90, 60), D(90, -20), D(-20, 60)}, 1}, rect(-30, -30, 100, 100)};
  GSHGEOM_Flat_Path * fp = make(in);
  const uint64_t before = gshgeom_flat_path_hash(fp);
  uint64_t first = 0;
  for (int i = 0; i < 3; i++) {
    GSHGEOM_Flat_Path * out = nullptr;
    ASSERT_EQ(gshgeom_flat_path_clip_convex(fp, kCell.data(), nullptr, nullptr, &out), GSHGEOM_OK);
    const uint64_t h = gshgeom_flat_path_hash(out);
    if (i == 0) first = h;
    EXPECT_EQ(h, first);
    gshgeom_flat_path_destroy(out);
  }
  EXPECT_EQ(gshgeom_flat_path_hash(fp), before);
  gshgeom_flat_path_destroy(fp);
}

/* --- the random sweep: CK2, CK4, CK5 -------------------------------------------------------------------------------------- */

struct Rng {
  uint64_t s;
  uint64_t next() {
    s ^= s >> 12;
    s ^= s << 25;
    s ^= s >> 27;
    return s * UINT64_C(2685821657736338717);
  }
  int64_t below(int64_t n) { return (int64_t)(next() % (uint64_t)n); }
};

bool convex_valid(const Quad & q) {
  int sigma = 0;
  int t[4];
  for (int i = 0; i < 4; i++) {
    const Dp a = q[(i + 3) % 4], b = q[i], c = q[(i + 1) % 4];
    const int64_t v = ((int64_t)b.x - a.x) * ((int64_t)c.y - a.y) - ((int64_t)b.y - a.y) * ((int64_t)c.x - a.x);
    t[i] = (v > 0) - (v < 0);
    if (t[i]) {
      if (sigma && sigma != t[i]) return false;
      sigma = t[i];
    }
  }
  return sigma != 0 && (t[1] != 0 || t[3] != 0);
}

/** The ratio of two big integers as a long double, by their leading limbs (the denominators here can pass 2^16384). */
long double ratio_ld(const co::Big & n, const co::Big & d) {
  if (n.m.empty()) return 0.0L;
  auto lead = [](const co::Big & b, long * exp) {
    long double v = 0;
    const size_t take = std::min<size_t>(b.m.size(), 4);
    for (size_t i = 0; i < take; i++) v = v * 4294967296.0L + (long double)b.m[b.m.size() - 1 - i];
    *exp = (long)(b.m.size() - take);
    return v;
  };
  long en, ed;
  const long double vn = lead(n, &en), vd = lead(d, &ed);
  const long double r = std::ldexp(vn / vd, (int)(32 * (en - ed)));
  return (n.neg != d.neg) ? -r : r;
}

long double to_ld(const co::Rat & a) {
  return ratio_ld(a.n, a.d);
}

/** The proved bound on the seam area, in square units: `c` crossings, each displaced by at most half a unit in each coordinate (sqrt(2)/2
 * in all), moves the two result edges at it, each at most the diameter of Q long inside Q, so the area changes by at most half the displacement
 * times the two lengths: sqrt(2)/2 * diameter per crossing.  (Documented in documentation/design.md.) */
long double proved_bound(size_t crossings, const Quad & q) {
  long double diam = 0;
  for (const Dp & a : q) {
    for (const Dp & b : q) diam = std::max(diam, std::hypot((long double)a.x - b.x, (long double)a.y - b.y));
  }
  return (long double)crossings * 0.70711L * diam;
}

struct Worst {
  long double literal = 0;      // area / (P / 128 px) with P the lengths of the sides crossed: recorded, not gated (it is not a bound)
  long double literal_cells = 0;  // the same over the cases whose quadrilateral is not thin (its sides within a factor of 3)
  long double proved = 0;       // area / (c * sqrt(2)/2 * diameter of Q): the bound that holds for every cell shape
  long cases = 0;
  long gated = 0;
};

TEST(ClipSweep, AThinCellShowsWhatTheLiteralBoundDoesNotCover) {
  // A cell 64 units wide and 6400 tall, and an edge that enters through its short side at a shallow angle and runs along
  // it: the crossing is rounded by up to half a unit perpendicular to an edge 6400 units long, which is a sliver of
  // up to 1600 square units, against P / 2 = 32 for a clip line of 64.
  const Quad thin = {D(0, 0), D(64, 0), D(64, 6400), D(0, 6400)};
  long double worst = 0;
  for (int dx = 1; dx <= 8; dx++) {
    // From (30, -6000) to (30 + dx, 6000): across y = 0 at x = 30 + dx / 2, a tie for odd dx.
    Contours in = {Contour{{D(30, -6000), D(30 + dx, 6000), D(60, 6000), D(60, -6000)}, 1}};
    Contours out;
    ASSERT_EQ(clip(in, thin, &out), GSHGEOM_OK);
    const co::AreaResult ar = co::symmetric_difference_area(polys(in), polys(out), quad_poly(thin));
    const long double area = to_ld(ar.area);
    worst = std::max(worst, area / (64.0L / 2));
    const long double bound = proved_bound(2, thin);  // the edge crosses y = 0 and, rounded, the contour's other edges do not
    EXPECT_LE(area, bound);
  }
  EXPECT_GT(worst, 10.0L) << "the literal P/128 px figure is exceeded by far in a thin cell, which is why the proved bound is the one";
  std::printf("[clip] thin cell: worst area / (P/128 px) = %.3Lf\n", worst);
}

TEST(ClipSweep, SeededRandomCasesAgreeWithTheOracleAndRecordTheSeamBound) {
  Rng rng{UINT64_C(0x5EA11CAFE0DDBA11)};
  Worst worst;
  ProbeStats st;
  long copied_checked = 0, max_c_over_n_x1000 = 0, empties = 0, cut_cases = 0;
  const int kCases = 360;
  for (int i = 0; i < kCases; i++) {
    // Draws, in order: scale class (1), centre (2), half sizes (2), shear (1), orientation (1), degenerate flag (1), four
    // corner jitters (4 when the general kind is drawn); then the contour count (1), and for each contour the point count (1),
    // each point's kind (1) and coordinates (2 or 1), and the closed flag (1).
    const int64_t span = (i % 6 == 5) ? (int64_t)1 << 20 : 300;
    const int64_t cx = rng.below(span) - span / 2, cy = rng.below(span) - span / 2;
    const int64_t hx = 8 + rng.below(span / 3), hy = 8 + rng.below(span / 3);
    const int64_t shear = rng.below(2 * (hx / 4) + 1) - hx / 4;
    const int flip = (int)rng.below(2);
    const int degenerate = (int)rng.below(6);
    Quad q = {D((int32_t)(cx - hx), (int32_t)(cy - hy)), D((int32_t)(cx + hx), (int32_t)(cy - hy + shear)),
        D((int32_t)(cx + hx + shear), (int32_t)(cy + hy)), D((int32_t)(cx - hx), (int32_t)(cy + hy - shear))};
    if (i % 4 == 3) {
      for (Dp & p : q) {
        p.x += (int32_t)(rng.below(hx / 2 + 1) - hx / 4);
        p.y += (int32_t)(rng.below(hy / 2 + 1) - hy / 4);
      }
    }
    if (flip) std::reverse(q.begin(), q.end());
    if (degenerate == 0) q[3] = q[2];
    if (!convex_valid(q)) {
      q = {D((int32_t)(cx - hx), (int32_t)(cy - hy)), D((int32_t)(cx + hx), (int32_t)(cy - hy)), D((int32_t)(cx + hx), (int32_t)(cy + hy)),
          D((int32_t)(cx - hx), (int32_t)(cy + hy))};
    }
    Contours in;
    const int nc = 1 + (int)rng.below(3);
    for (int c = 0; c < nc; c++) {
      const int np = 1 + (int)rng.below(8);
      Contour ct;
      for (int k = 0; k < np; k++) {
        const int kind = (int)rng.below(8);
        if (kind == 0) {  // a corner of the quadrilateral
          ct.pts.push_back(q[rng.below(4)]);
          continue;
        }
        if (kind == 1) {  // a point on a side of it, at a drawn fraction of the side, rounded to the grid
          const Dp a = q[rng.below(4)], b = q[(rng.below(4) + 1) % 4];
          const int64_t f = rng.below(9);
          ct.pts.push_back(D((int32_t)(a.x + (b.x - a.x) * f / 8), (int32_t)(a.y + (b.y - a.y) * f / 8)));
          continue;
        }
        const int64_t px = cx + rng.below(2 * span) - span;
        const int64_t py = cy + rng.below(2 * span) - span;
        ct.pts.push_back(D((int32_t)px, (int32_t)py));
      }
      ct.closed = rng.below(4) != 0;
      in.push_back(ct);
    }
    Contours out;
    GSHGEOM_Clip_Stats stats{};
    ASSERT_EQ(clip(in, q, &out, &stats), GSHGEOM_OK) << i;
    worst.cases++;
    empties += out.empty();
    check_winding(in, q, out, &st, "sweep case " + std::to_string(i));
    if (HasFailure()) {
      std::printf("[clip] failing case %d: quad", i);
      for (const Dp & p : q) std::printf(" (%d,%d)", p.x, p.y);
      std::printf("\n");
      return;
    }
    // CK5: points out <= points in + crossings.
    size_t n_in = 0, n_out = 0;
    for (const Contour & c : in) n_in += c.pts.size();
    for (const Contour & c : out) n_out += c.pts.size();
    EXPECT_LE(n_out, n_in + stats.crossings) << i;
    // CK5's proved bound on the crossings: at most 65/8 = 8.125 per input point (not the unproved 4).
    EXPECT_LE(8 * stats.crossings, 65 * n_in) << i;
    if (n_in) max_c_over_n_x1000 = std::max<long>(max_c_over_n_x1000, (long)(1000 * stats.crossings / n_in));
    // A contour entirely inside the closed quadrilateral is in the output, bit for bit.
    {
      const co::Polygon pq = quad_poly(q);
      int sigma = co::quad_sigma(pq);
      for (const Contour & c : in) {
        bool inside = true;
        for (const Dp & p : c.pts) {
          for (int s = 0; s < 4 && inside; s++) {
            const co::IPt a = pq[s], b = pq[(s + 1) % 4];
            if (a == b) continue;
            if (co::side_at(a, b, co::RPt{co::Rat(p.x), co::Rat(p.y)}).sign() * sigma < 0) inside = false;
          }
        }
        if (!inside) continue;
        bool found = false;
        for (const Contour & o : out) found = found || o == c;
        EXPECT_TRUE(found) << "case " << i << ": an inside contour is not copied";
        copied_checked++;
      }
    }
    // CK4: the area of the symmetric difference against the exact clip, and the bound it is gated by.
    {
      const co::Polygons pin = polys(in), pout = polys(out);
      const co::Polygon pq = quad_poly(q);
      const co::AreaResult ar = co::symmetric_difference_area(pin, pout, pq);
      // P: the lengths of the sides whose line an input edge crosses strictly, in units; the bound is P/2 square units
      // (P / 128 pixels with 64 units to the pixel is P / 8192 square pixels = P / 2 square units).
      long double P = 0;
      for (int s = 0; s < 4; s++) {
        const Dp a = q[s], b = q[(s + 1) % 4];
        if (a.x == b.x && a.y == b.y) continue;
        bool crossed = false;
        for (const Contour & c : in) {
          for (size_t k = 0; k < c.pts.size() && !crossed; k++) {
            const Dp u = c.pts[k], v = c.pts[(k + 1) % c.pts.size()];
            crossed = gshgeom_clip_side(a, b, u) * gshgeom_clip_side(a, b, v) < 0;
          }
        }
        if (crossed) P += std::hypot((long double)b.x - a.x, (long double)b.y - a.y);
      }
      const long double area = to_ld(ar.area);
      if (P > 0) {
        worst.gated++;
        worst.literal = std::max(worst.literal, area / (P / 2));
        if (std::max(hx, hy) <= 3 * std::min(hx, hy) && stats.crossings <= 8) {
          worst.literal_cells = std::max(worst.literal_cells, area / (P / 2));
        }
      } else {
        EXPECT_LE(area, 0.0L) << "case " << i << ": nothing crossed a side, so the cut changes nothing, and yet the area differs";
      }
      const long double bound = proved_bound(stats.crossings, q);
      EXPECT_LE(area, bound) << "case " << i << ": the seam area passes the proved bound";
      if (bound > 0) worst.proved = std::max(worst.proved, area / bound);
      cut_cases += stats.contours_in > 0;
    }
  }
  EXPECT_GT(st.compared, 20000);
  EXPECT_GT(st.nonzero_in, 800);
  EXPECT_GT(copied_checked, 20);
  EXPECT_GT(cut_cases, 100);
  EXPECT_GT(empties, 5);
  std::printf("[clip] sweep: %ld cases, %ld probes compared (%ld nonzero), %ld inside contours checked for the bit copy, "
              "%ld cases cut, worst crossings/points %.3f\n",
      worst.cases, st.compared, st.nonzero_in, copied_checked, cut_cases, (double)max_c_over_n_x1000 / 1000.0);
  std::printf("[clip] CK4 worst seam-bound ratio over %ld gated cases: area / (P/128 px) = %.4Lf (%.4Lf over the cases with a "
              "non-thin cell and at most 8 crossings); against the proved bound (crossings * sqrt(2)/2 * diameter) = %.5Lf\n",
      worst.gated, worst.literal, worst.literal_cells, worst.proved);
  RecordProperty("ck4_worst_ratio", std::to_string((double)worst.literal));
  RecordProperty("ck4_worst_ratio_proved", std::to_string((double)worst.proved));
}

}  // namespace

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
