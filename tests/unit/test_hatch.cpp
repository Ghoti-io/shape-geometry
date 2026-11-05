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
 * Hatches (hatch.h, HX2 and HX4 to HX12) against the independent known answers: exact rationals from the definition
 * (tools/gen-hatch-vectors.py), and the structural claims (the check's order, the limits, the differential against the
 * stroke of the same pieces).
 */

#include "test_helpers.h"

extern "C" {
#include "../../src/stroke/stroke_internal.h"
}
#include "../hatch_oracle.h"
#include "../hatch_vectors.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <deque>
#include <string>
#include <cstdio>
#include <cstring>
#include <vector>

namespace {

namespace hv = hatch_vectors;

GSHGEOM_Map map6(const double m[6]) {
  const double lin[4] = {m[0], m[1], m[2], m[3]};
  const double t[2] = {m[4], m[5]};
  GSHGEOM_Map out;
  EXPECT_EQ(gshgeom_map_from_matrix(2, lin, t, &out), GSHGEOM_OK);
  return out;
}

struct Flat {
  GSHGEOM_Flat_Path * p = nullptr;
  Flat(const std::vector<std::vector<GSHGEOM_Device_Point>> & contours) {
    std::vector<GSHGEOM_Device_Point> pts;
    std::vector<GSHGEOM_Flat_Subpath> subs;
    for (const auto & c : contours) {
      GSHGEOM_Flat_Subpath s;
      s.first = pts.size();
      s.count = c.size();
      s.closed = 1;
      subs.push_back(s);
      pts.insert(pts.end(), c.begin(), c.end());
    }
    EXPECT_EQ(gshgeom_flat_path_from_points(pts.data(), pts.size(), subs.data(), subs.size(), nullptr, nullptr, &p), GSHGEOM_OK);
  }
  ~Flat() { gshgeom_flat_path_destroy(p); }
};

std::vector<std::vector<GSHGEOM_Device_Point>> contours_of(const hv::VCase & c) {
  std::vector<std::vector<GSHGEOM_Device_Point>> out;
  const int32_t * q = c.points;
  for (int i = 0; i < c.contour_count; i++) {
    std::vector<GSHGEOM_Device_Point> one;
    for (int k = 0; k < c.contour_sizes[i]; k++) {
      GSHGEOM_Device_Point d;
      d.x = q[0];
      d.y = q[1];
      q += 2;
      one.push_back(d);
    }
    out.push_back(one);
  }
  return out;
}

bool near(double a, double b) {
  const double scale = std::fmax(1.0, std::fmax(std::fabs(a), std::fabs(b)));
  return std::fabs(a - b) <= std::ldexp(scale, -40);
}

struct HatchOwner {
  std::vector<GSHGEOM_Hatch_Family> fams;
  GSHGEOM_Hatch h;
  HatchOwner(const hv::VCase & c) {
    for (int i = 0; i < c.family_count; i++) {
      GSHGEOM_Hatch_Family f;
      f.angle = c.families[i].angle;
      f.origin_x = c.families[i].ox;
      f.origin_y = c.families[i].oy;
      f.offset_along = c.families[i].dx;
      f.offset_across = c.families[i].dy;
      f.dash_count = (size_t)c.families[i].dash_count;
      f.dashes = c.families[i].dashes;
      fams.push_back(f);
    }
    std::memset(&h, 0, sizeof h);
    h.angle_unit = GSHGEOM_ANGLE_DEGREES;
    h.nesting = (uint32_t)c.nesting;
    h.map = map6(c.hatch_map);
    h.family_count = fams.size();
    h.families = fams.data();
  }
};

}  // namespace

// The quarter turns and the 45, 30 and 60 degrees the vectors use are what the generator wrote.
TEST(HatchVectors, TheUnitVectorsAreTheGeneratorsOwn) {
  const double s2 = 0.70710678118654746, s2b = 0.70710678118654757, s3 = 0.86602540378443871, hf = 0.49999999999999994;
  const struct {
    double a, c, s;
  } t[] = {{0, 1, 0}, {90, 0, 1}, {180, -1, 0}, {270, 0, -1}, {45, s2, s2b}, {135, -s2b, s2}, {30, s3, hf}, {60, hf, s3}};
  for (const auto & e : t) {
    double c, s;
    ASSERT_EQ(gshgeom_sincos(GSHGEOM_ANGLE_DEGREES, e.a, &c, &s), GSHGEOM_OK);
    EXPECT_EQ(c, e.c) << e.a;
    EXPECT_EQ(s, e.s) << e.a;
  }
}

TEST(HatchVectors, EveryCaseMatchesTheExactAnswer) {
  for (size_t i = 0; i < hv::case_count; i++) {
    const hv::VCase & c = hv::cases[i];
    Flat region(contours_of(c));
    HatchOwner ho(c);
    GSHGEOM_Map chain = map6(c.chain);
    GSHGEOM_Hatch_Pieces * pieces = nullptr;
    ASSERT_EQ(gshgeom_hatch_pieces_from_flat(region.p, c.rule ? GSHGEOM_FILL_EVEN_ODD : GSHGEOM_FILL_NONZERO, &ho.h, &chain, 0,
                  nullptr, nullptr, &pieces), GSHGEOM_OK) << c.name;
    EXPECT_EQ(gshgeom_hatch_pieces_lines(pieces), c.lines) << c.name;
    ASSERT_EQ(gshgeom_hatch_pieces_count(pieces), (size_t)c.piece_count) << c.name;
    for (int k = 0; k < c.piece_count; k++) {
      double xy[4];
      int dot = 0;
      ASSERT_EQ(gshgeom_hatch_pieces_get(pieces, (size_t)k, xy, &dot), GSHGEOM_OK);
      EXPECT_EQ(dot, c.pieces[k].dot) << c.name << " piece " << k;
      for (int q = 0; q < 4; q++) {
        EXPECT_TRUE(near(xy[q], c.pieces[k].p[q])) << c.name << " piece " << k << " coord " << q << ": " << xy[q] << " vs "
                                                   << c.pieces[k].p[q];
      }
    }
    gshgeom_hatch_pieces_destroy(pieces);
  }
}

// --- helpers for the tests below --------------------------------------------------------------------------------------------

namespace {

namespace ho = hatch_oracle;

struct Rng {
  uint64_t s;
  uint64_t next() {
    s ^= s >> 12;
    s ^= s << 25;
    s ^= s >> 27;
    return s * UINT64_C(2685821657736338717);
  }
  int64_t range(int64_t a, int64_t b) { return a + (int64_t)(next() % (uint64_t)(b - a + 1)); }
};

GSHGEOM_Hatch one_family(GSHGEOM_Hatch_Family & f, double angle, double ox, double oy, double dx, double dy, const std::vector<double> & dashes,
                         const GSHGEOM_Map & map) {
  f.angle = angle;
  f.origin_x = ox;
  f.origin_y = oy;
  f.offset_along = dx;
  f.offset_across = dy;
  static std::deque<std::vector<double>> pool;  // the family borrows its elements, so they outlive the call
  pool.push_back(dashes);
  f.dash_count = dashes.size();
  f.dashes = dashes.empty() ? nullptr : pool.back().data();
  GSHGEOM_Hatch h;
  std::memset(&h, 0, sizeof h);
  h.angle_unit = GSHGEOM_ANGLE_DEGREES;
  h.nesting = GSHGEOM_HATCH_NESTING_RULE;
  h.map = map;
  h.family_count = 1;
  h.families = &f;
  return h;
}

GSHGEOM_Map identity() {
  GSHGEOM_Map m;
  EXPECT_EQ(gshgeom_map_identity(2, &m), GSHGEOM_OK);
  return m;
}

GSHGEOM_Path * rect_path(double x, double y, double w, double h) {
  GSHGEOM_Path * p = nullptr;
  EXPECT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &p), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_add_rect(p, x, y, w, h, 0, 0), GSHGEOM_OK);
  return p;
}

GSHGEOM_Stroke style(double width, uint32_t cap, uint32_t space = GSHGEOM_STROKE_WIDTH_LOCAL) {
  GSHGEOM_Stroke st;
  gshgeom_stroke_init(&st);
  st.width = width;
  st.cap = cap;
  st.width_space = space;
  return st;
}

/** The winding number of the 26.6 flat path at the point (x, y) in pixels (a ray in +x; the point is never on an edge in the tests). */
int winding(const GSHGEOM_Flat_Path * f, double x, double y) {
  const GSHGEOM_Device_Point * pts = gshgeom_flat_path_points(f);
  int w = 0;
  for (size_t si = 0; si < gshgeom_flat_path_subpath_count(f); si++) {
    GSHGEOM_Flat_Subpath s;
    EXPECT_EQ(gshgeom_flat_path_subpath(f, si, &s), GSHGEOM_OK);
    for (size_t i = 0; i < s.count; i++) {
      const GSHGEOM_Device_Point a = pts[s.first + i], b = pts[s.first + (i + 1) % s.count];
      const double ay = a.y / 64.0, by = b.y / 64.0;
      if ((ay <= y) != (by <= y)) {
        const double t = (y - ay) / (by - ay);
        if (a.x / 64.0 + t * ((b.x - a.x) / 64.0) > x) w += by > ay ? 1 : -1;
      }
    }
  }
  return w;
}

}  // namespace

// --- HX2: the geometric checks, in order ------------------------------------------------------------------------------------

TEST(HatchCheck, TheOrderIsInvalidThenUnsupportedThenLimit) {
  GSHGEOM_Hatch_Family f;
  std::vector<double> dashes = {1.0, -1.0};
  GSHGEOM_Hatch h = one_family(f, 45, 0, 0, 0, 1, dashes, identity());
  EXPECT_EQ(gshgeom_hatch_check(&h), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_hatch_check(nullptr), GSHGEOM_ERR_INVALID);
  // a limit and an invalid value together: invalid wins
  GSHGEOM_Hatch_Family g = f;
  g.offset_across = std::ldexp(1.0, 65);
  GSHGEOM_Hatch hg = h;
  hg.families = &g;
  EXPECT_EQ(gshgeom_hatch_check(&hg), GSHGEOM_ERR_LIMIT);
  hg.nesting = 3;
  EXPECT_EQ(gshgeom_hatch_check(&hg), GSHGEOM_ERR_INVALID);
  hg.nesting = 0;
  // a 3D map is unsupported, and an invalid map is invalid
  GSHGEOM_Hatch h3 = h;
  ASSERT_EQ(gshgeom_map_identity(3, &h3.map), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_hatch_check(&h3), GSHGEOM_ERR_UNSUPPORTED);
  GSHGEOM_Hatch h3l = h3;
  g.offset_across = std::ldexp(1.0, 65);
  h3l.families = &g;
  EXPECT_EQ(gshgeom_hatch_check(&h3l), GSHGEOM_ERR_UNSUPPORTED) << "unsupported is decided before the limits";
  GSHGEOM_Hatch hb = h;
  hb.map.m[0][0] = std::nan("");
  EXPECT_EQ(gshgeom_hatch_check(&hb), GSHGEOM_ERR_INVALID);
  GSHGEOM_Hatch hu = h;
  hu.angle_unit = GSHGEOM_ANGLE_UNIT_COUNT;
  EXPECT_EQ(gshgeom_hatch_check(&hu), GSHGEOM_ERR_INVALID);
  GSHGEOM_Hatch hz = h;
  hz.family_count = 0;
  EXPECT_EQ(gshgeom_hatch_check(&hz), GSHGEOM_ERR_INVALID);
  GSHGEOM_Hatch hn = h;
  hn.families = nullptr;
  EXPECT_EQ(gshgeom_hatch_check(&hn), GSHGEOM_ERR_INVALID);
}

TEST(HatchCheck, EachBoundHoldsAtTheBoundAndFailsOnePast) {
  GSHGEOM_Hatch_Family f;
  std::vector<double> dashes = {1.0, -1.0};
  GSHGEOM_Hatch h = one_family(f, 45, 0, 0, 0.5, 1, dashes, identity());
  struct Probe {
    const char * what;
    void (*set)(GSHGEOM_Hatch_Family &, double);
    double ok, bad;
  };
  const Probe probes[] = {
      {"origin x", [](GSHGEOM_Hatch_Family & g, double v) { g.origin_x = v; }, GSHGEOM_HATCH_VALUE_MAX, std::ldexp(1.0, 65)},
      {"origin y low", [](GSHGEOM_Hatch_Family & g, double v) { g.origin_y = v; }, GSHGEOM_HATCH_VALUE_MIN, std::ldexp(1.0, -65)},
      {"along", [](GSHGEOM_Hatch_Family & g, double v) { g.offset_along = v; }, -GSHGEOM_HATCH_VALUE_MAX, -std::ldexp(1.0, 65)},
      {"across", [](GSHGEOM_Hatch_Family & g, double v) { g.offset_across = v; }, GSHGEOM_HATCH_VALUE_MIN, std::ldexp(1.0, -65)},
      {"angle", [](GSHGEOM_Hatch_Family & g, double v) { g.angle = v; }, std::ldexp(1.0, 40), std::ldexp(1.0, 41)},
  };
  for (const Probe & p : probes) {
    GSHGEOM_Hatch_Family g = f;
    p.set(g, p.ok);
    GSHGEOM_Hatch hh = h;
    hh.families = &g;
    EXPECT_NE(gshgeom_hatch_check(&hh), GSHGEOM_ERR_INVALID) << p.what;
    if (std::string(p.what) != "angle") {
      EXPECT_EQ(gshgeom_hatch_check(&hh), GSHGEOM_OK) << p.what;
    }
    p.set(g, p.bad);
    EXPECT_EQ(gshgeom_hatch_check(&hh), GSHGEOM_ERR_LIMIT) << p.what;
  }
  // a zero value is allowed wherever the band allows zero, a dash too; a dash past the band is a limit
  GSHGEOM_Hatch_Family g = f;
  g.origin_x = 0.0;
  std::vector<double> big = {GSHGEOM_HATCH_VALUE_MAX, -1.0};
  g.dashes = big.data();
  GSHGEOM_Hatch hh = h;
  hh.families = &g;
  EXPECT_EQ(gshgeom_hatch_check(&hh), GSHGEOM_OK);
  big[0] = std::ldexp(1.0, 65);
  EXPECT_EQ(gshgeom_hatch_check(&hh), GSHGEOM_ERR_LIMIT);
  big[0] = std::ldexp(1.0, -65);
  EXPECT_EQ(gshgeom_hatch_check(&hh), GSHGEOM_ERR_LIMIT);
  // the map's entries are in the same band
  GSHGEOM_Hatch hm = h;
  hm.map.m[0][0] = std::ldexp(1.0, 65);
  hm.map.kind = 0;
  EXPECT_EQ(gshgeom_hatch_check(&hm), GSHGEOM_ERR_LIMIT);
  // families
  std::vector<GSHGEOM_Hatch_Family> many(GSHGEOM_HATCH_FAMILIES_MAX + 1, f);
  GSHGEOM_Hatch hf = h;
  hf.families = many.data();
  hf.family_count = GSHGEOM_HATCH_FAMILIES_MAX;
  EXPECT_EQ(gshgeom_hatch_check(&hf), GSHGEOM_OK);
  hf.family_count = GSHGEOM_HATCH_FAMILIES_MAX + 1;
  EXPECT_EQ(gshgeom_hatch_check(&hf), GSHGEOM_ERR_LIMIT);
  // dashes: 65535 elements pass, one more is a limit; elements with no pen-down element or no length are invalid
  std::vector<double> pat(GSHGEOM_HATCH_DASHES_MAX, 1.0);
  GSHGEOM_Hatch_Family fd = f;
  fd.dash_count = pat.size();
  fd.dashes = pat.data();
  GSHGEOM_Hatch hd = h;
  hd.families = &fd;
  EXPECT_EQ(gshgeom_hatch_check(&hd), GSHGEOM_OK);
  pat.push_back(1.0);
  fd.dash_count = pat.size();
  fd.dashes = pat.data();
  EXPECT_EQ(gshgeom_hatch_check(&hd), GSHGEOM_ERR_LIMIT);
  std::vector<double> spaces = {-1.0, -2.0};
  fd.dash_count = 2;
  fd.dashes = spaces.data();
  EXPECT_EQ(gshgeom_hatch_check(&hd), GSHGEOM_ERR_INVALID) << "no pen-down element";
  std::vector<double> zero = {0.0};
  fd.dash_count = 1;
  fd.dashes = zero.data();
  EXPECT_EQ(gshgeom_hatch_check(&hd), GSHGEOM_ERR_INVALID) << "no length";
  fd.dash_count = 1;
  fd.dashes = nullptr;
  EXPECT_EQ(gshgeom_hatch_check(&hd), GSHGEOM_ERR_INVALID) << "dashes counted and not given";
  // a dash element that is not finite is invalid even when a limit also holds
  std::vector<double> nan_dash = {std::nan(""), 1.0};
  fd.dash_count = 2;
  fd.dashes = nan_dash.data();
  EXPECT_EQ(gshgeom_hatch_check(&hd), GSHGEOM_ERR_INVALID);
}

// --- the lines (HX5) ----------------------------------------------------------------------------------------------------------

TEST(HatchLines, AQuarterTurnIsExactlyAxisAligned) {
  // the lines of a 90-degree family are exactly vertical: every piece has one x
  GSHGEOM_Hatch_Family f;
  GSHGEOM_Hatch h = one_family(f, 90, 0.0, 0.0, 0.0, 1.3, {}, identity());
  GSHGEOM_Path * region = rect_path(0.5, 0.5, 9.0, 7.0);
  GSHGEOM_Map chain = identity();
  GSHGEOM_Hatch_Pieces * p = nullptr;
  ASSERT_EQ(gshgeom_hatch_pieces(region, GSHGEOM_FILL_NONZERO, &h, &chain, 0, nullptr, nullptr, &p), GSHGEOM_OK);
  ASSERT_GT(gshgeom_hatch_pieces_count(p), 0u);
  for (size_t i = 0; i < gshgeom_hatch_pieces_count(p); i++) {
    double xy[4];
    int dot;
    ASSERT_EQ(gshgeom_hatch_pieces_get(p, i, xy, &dot), GSHGEOM_OK);
    EXPECT_EQ(xy[0], xy[2]) << i;
    EXPECT_EQ(xy[1], 0.5) << i;
    EXPECT_EQ(xy[3], 7.5) << i;
    EXPECT_FALSE(std::signbit(xy[0] * 0.0) && xy[0] == 0.0 && false);
  }
  gshgeom_hatch_pieces_destroy(p);
  gshgeom_path_destroy(region);
}

TEST(HatchLines, TheOriginIsAppliedInHatchSpaceBeforeTheMap) {
  // origin (3, 2), angle 0, spacing 4: lines at hatch y = 2 + 4k; the hatch map rotates a quarter turn about the origin, so
  // the lines in local space are vertical at x = -(2 + 4k): the origin moved with the lines, not after the rotation.
  GSHGEOM_Map rot;
  ASSERT_EQ(gshgeom_map_rotate(GSHGEOM_ANGLE_DEGREES, 90, &rot), GSHGEOM_OK);
  GSHGEOM_Hatch_Family f;
  GSHGEOM_Hatch h = one_family(f, 0, 3.0, 2.0, 0.0, 4.0, {}, rot);
  GSHGEOM_Path * region = rect_path(-20, -20, 40, 40);
  GSHGEOM_Map chain = identity();
  GSHGEOM_Hatch_Pieces * p = nullptr;
  ASSERT_EQ(gshgeom_hatch_pieces(region, GSHGEOM_FILL_NONZERO, &h, &chain, 0, nullptr, nullptr, &p), GSHGEOM_OK);
  size_t n = gshgeom_hatch_pieces_count(p);
  ASSERT_GT(n, 4u);
  for (size_t i = 0; i < n; i++) {
    double xy[4];
    int dot;
    ASSERT_EQ(gshgeom_hatch_pieces_get(p, i, xy, &dot), GSHGEOM_OK);
    const double r = std::fmod(std::fabs(xy[0]) - 2.0, 4.0);
    EXPECT_TRUE(r == 0.0 || std::fabs(std::fabs(r) - 4.0) < 1e-9) << xy[0];
    EXPECT_NEAR(xy[0], xy[2], 1e-9);
  }
  gshgeom_hatch_pieces_destroy(p);
  gshgeom_path_destroy(region);
}

TEST(HatchLines, TheLinesOfTheUseAreCountedAndLimited) {
  GSHGEOM_Hatch_Family f;
  GSHGEOM_Hatch h = one_family(f, 0, 0, 0, 0, 1.0, {}, identity());
  GSHGEOM_Path * region = rect_path(0, 0, 10, 10);  // lines y = 0 .. 10: eleven of them (the closed bounds)
  GSHGEOM_Map chain = identity();
  GSHGEOM_Hatch_Pieces * p = nullptr;
  ASSERT_EQ(gshgeom_hatch_pieces(region, GSHGEOM_FILL_NONZERO, &h, &chain, 11, nullptr, nullptr, &p), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_hatch_pieces_lines(p), 11u);
  EXPECT_EQ(gshgeom_hatch_pieces_count(p), 11u);
  gshgeom_hatch_pieces_destroy(p);
  p = nullptr;
  EXPECT_EQ(gshgeom_hatch_pieces(region, GSHGEOM_FILL_NONZERO, &h, &chain, 10, nullptr, nullptr, &p), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(p, nullptr);
  // a spacing of 2^-64 over a full target is a limit decided before a line is made
  GSHGEOM_Hatch_Family tiny;
  GSHGEOM_Hatch ht = one_family(tiny, 0, 0, 0, 0, std::ldexp(1.0, -64), {}, identity());
  EXPECT_EQ(gshgeom_hatch_pieces(region, GSHGEOM_FILL_NONZERO, &ht, &chain, 0, nullptr, nullptr, &p), GSHGEOM_ERR_LIMIT);
  // the lines of two families add up
  GSHGEOM_Hatch_Family two[2] = {f, f};
  GSHGEOM_Hatch h2 = h;
  h2.families = two;
  h2.family_count = 2;
  ASSERT_EQ(gshgeom_hatch_pieces(region, GSHGEOM_FILL_NONZERO, &h2, &chain, 22, nullptr, nullptr, &p), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_hatch_pieces_lines(p), 22u);
  gshgeom_hatch_pieces_destroy(p);
  p = nullptr;
  EXPECT_EQ(gshgeom_hatch_pieces(region, GSHGEOM_FILL_NONZERO, &h2, &chain, 21, nullptr, nullptr, &p), GSHGEOM_ERR_LIMIT);
  gshgeom_path_destroy(region);
}

TEST(HatchLines, AFamilyWithNoSpacingDrawsNothingAndTheOthersDo) {
  GSHGEOM_Hatch_Family fam[2];
  one_family(fam[0], 0, 0, 0, 0, 0.0, {}, identity());
  one_family(fam[1], 90, 0, 0, 0, 5.0, {}, identity());
  GSHGEOM_Hatch h = one_family(fam[0], 0, 0, 0, 0, 0.0, {}, identity());
  h.families = fam;
  h.family_count = 2;
  GSHGEOM_Path * region = rect_path(0, 0, 10, 10);
  GSHGEOM_Map chain = identity();
  GSHGEOM_Hatch_Pieces * p = nullptr;
  ASSERT_EQ(gshgeom_hatch_pieces(region, GSHGEOM_FILL_NONZERO, &h, &chain, 0, nullptr, nullptr, &p), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_hatch_pieces_lines(p), 3u);
  EXPECT_EQ(gshgeom_hatch_pieces_count(p), 3u);
  gshgeom_hatch_pieces_destroy(p);
  gshgeom_path_destroy(region);
}

// --- the stroke of the pieces (HX7, HX8) ----------------------------------------------------------------------------------------

TEST(HatchExpand, TheOutlineIsTheStrokeOfTheExplicitPieces) {
  // the differential of HX8: the hatch outline's hash equals the hash of the stroke expansion of the pieces as an explicit
  // path, one open two-point subpath each, with the same style and no dashes
  struct Spec {
    const char * name;
    double angle, dy;
    std::vector<double> dashes;
    uint32_t cap;
    uint32_t space;
    double width;
    bool nonuniform;
  };
  const Spec specs[] = {
      {"solid45 butt", 45, 1.5, {}, GSHGEOM_STROKE_CAP_BUTT, GSHGEOM_STROKE_WIDTH_LOCAL, 0.3, false},
      {"dashed round", 0, 2.0, {3.0, -1.0}, GSHGEOM_STROKE_CAP_ROUND, GSHGEOM_STROKE_WIDTH_LOCAL, 0.5, false},
      {"dashed square device", 30, 1.75, {2.0, -1.0, 0.5, -1.5}, GSHGEOM_STROKE_CAP_SQUARE, GSHGEOM_STROKE_WIDTH_DEVICE, 1.0, false},
      {"dotted round", 0, 2.0, {0.0, -1.5}, GSHGEOM_STROKE_CAP_ROUND, GSHGEOM_STROKE_WIDTH_LOCAL, 0.6, false},
      {"dotted butt", 0, 2.0, {0.0, -1.5}, GSHGEOM_STROKE_CAP_BUTT, GSHGEOM_STROKE_WIDTH_LOCAL, 0.6, false},
      {"nonuniform local", 90, 1.1, {1.0, -0.5}, GSHGEOM_STROKE_CAP_SQUARE, GSHGEOM_STROKE_WIDTH_LOCAL, 0.4, true},
  };
  for (const Spec & sp : specs) {
    GSHGEOM_Hatch_Family f;
    GSHGEOM_Hatch h = one_family(f, sp.angle, 0.25, 0.5, 0.0, sp.dy, sp.dashes, identity());
    GSHGEOM_Path * region = rect_path(0.5, 0.5, 9.0, 7.0);
    const double m[4] = {sp.nonuniform ? 3.0 : 2.0, 0.0, 0.0, 2.0}, t[2] = {1.5, 2.25};
    GSHGEOM_Map chain;
    ASSERT_EQ(gshgeom_map_from_matrix(2, m, t, &chain), GSHGEOM_OK);
    GSHGEOM_Stroke st = style(sp.width, sp.cap, sp.space);
    GSHGEOM_Stroke_Outline * a = nullptr;
    ASSERT_EQ(gshgeom_hatch_expand(region, GSHGEOM_FILL_NONZERO, &h, &st, &chain, 0, nullptr, nullptr, &a), GSHGEOM_OK) << sp.name;
    GSHGEOM_Hatch_Pieces * p = nullptr;
    ASSERT_EQ(gshgeom_hatch_pieces(region, GSHGEOM_FILL_NONZERO, &h, &chain, 0, nullptr, nullptr, &p), GSHGEOM_OK);
    GSHGEOM_Path * explicit_path = nullptr;
    ASSERT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &explicit_path), GSHGEOM_OK);
    bool dots = false;
    for (size_t i = 0; i < gshgeom_hatch_pieces_count(p); i++) {
      double xy[4];
      int dot;
      ASSERT_EQ(gshgeom_hatch_pieces_get(p, i, xy, &dot), GSHGEOM_OK);
      dots |= dot != 0;
      GSHGEOM_Point a0 = {xy[0], xy[1], 0}, a1 = {xy[2], xy[3], 0};
      ASSERT_EQ(gshgeom_path_move_to(explicit_path, a0), GSHGEOM_OK);
      ASSERT_EQ(gshgeom_path_line_to(explicit_path, a1), GSHGEOM_OK);
    }
    GSHGEOM_Stroke_Outline * b = nullptr;
    ASSERT_EQ(gshgeom_stroke_expand(explicit_path, &st, &chain, nullptr, nullptr, &b), GSHGEOM_OK);
    if (!dots || sp.cap != GSHGEOM_STROKE_CAP_SQUARE) {
      EXPECT_EQ(gshgeom_stroke_outline_hash(a), gshgeom_stroke_outline_hash(b)) << sp.name;
    }
    EXPECT_EQ(gshgeom_stroke_outline_status(a), GSHGEOM_STROKE_DRAWN);
    EXPECT_EQ(gshgeom_stroke_outline_aligned(a), nullptr);
    EXPECT_EQ(gshgeom_stroke_outline_clip(a), nullptr);
    gshgeom_stroke_outline_destroy(a);
    gshgeom_stroke_outline_destroy(b);
    gshgeom_path_destroy(explicit_path);
    gshgeom_hatch_pieces_destroy(p);
    gshgeom_path_destroy(region);
  }
}

TEST(HatchExpand, ADotIsDrawnAlongItsLineWhateverDegenerateSays) {
  // dots along a 45-degree line under a square cap: the square of side 2h is turned with the line; omit does not drop it
  GSHGEOM_Hatch_Family f;
  GSHGEOM_Hatch h = one_family(f, 45, 0, 0, 0, 1000.0, {0.0, -10.0}, identity());
  GSHGEOM_Path * region = rect_path(-3, -3, 30, 30);
  GSHGEOM_Map chain = identity();
  for (uint32_t degenerate : {(uint32_t)GSHGEOM_STROKE_DEGENERATE_DRAW, (uint32_t)GSHGEOM_STROKE_DEGENERATE_OMIT}) {
    GSHGEOM_Stroke st = style(2.0, GSHGEOM_STROKE_CAP_SQUARE);
    st.degenerate = degenerate;
    GSHGEOM_Stroke_Outline * o = nullptr;
    ASSERT_EQ(gshgeom_hatch_expand(region, GSHGEOM_FILL_NONZERO, &h, &st, &chain, 0, nullptr, nullptr, &o), GSHGEOM_OK);
    const GSHGEOM_Flat_Path * body = gshgeom_stroke_outline_body(o);
    ASSERT_NE(body, nullptr);
    // the dot at the line's base point is at the origin (0, 0): the square has corners at distance sqrt(2) along the axes of
    // the line: the point (1.3, 0) is inside an axis-aligned square of side 2 only if |x| <= 1, and a square turned 45 degrees
    // reaches (sqrt 2 - small, 0)
    EXPECT_NE(winding(body, 1.3, 0.01), 0) << "turned square reaches 1.3 along an axis";
    EXPECT_EQ(winding(body, 1.3, 0.9), 0) << "past the turned square's side";
    EXPECT_EQ(winding(body, 0.95, 0.95), 0) << "the turned square does not reach the corner an axis-aligned one has";
    gshgeom_stroke_outline_destroy(o);
  }
  gshgeom_path_destroy(region);
}

TEST(HatchExpand, StatusesAreAStrokesAndTheStyleIsChecked) {
  GSHGEOM_Hatch_Family f;
  GSHGEOM_Hatch h = one_family(f, 0, 0, 0, 0, 2.0, {}, identity());
  GSHGEOM_Path * region = rect_path(0, 0, 10, 10);
  GSHGEOM_Map chain = identity();
  GSHGEOM_Stroke_Outline * o = nullptr;
  GSHGEOM_Stroke st = style(0.0, GSHGEOM_STROKE_CAP_BUTT);
  ASSERT_EQ(gshgeom_hatch_expand(region, GSHGEOM_FILL_NONZERO, &h, &st, &chain, 0, nullptr, nullptr, &o), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_stroke_outline_status(o), GSHGEOM_STROKE_EMPTY);
  gshgeom_stroke_outline_destroy(o);
  // width 0 is decided before the lines: a limit that the lines would hit is not reached
  ASSERT_EQ(gshgeom_hatch_expand(region, GSHGEOM_FILL_NONZERO, &h, &st, &chain, 1, nullptr, nullptr, &o), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_stroke_outline_status(o), GSHGEOM_STROKE_EMPTY);
  gshgeom_stroke_outline_destroy(o);
  st = style(1.0, GSHGEOM_STROKE_CAP_BUTT);
  GSHGEOM_Map sing;
  ASSERT_EQ(gshgeom_map_scale(2, 0.0, 1.0, 1.0, &sing), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_hatch_expand(region, GSHGEOM_FILL_NONZERO, &h, &st, &sing, 0, nullptr, nullptr, &o), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_stroke_outline_status(o), GSHGEOM_STROKE_NOT_INVERTIBLE) << "a singular chain";
  gshgeom_stroke_outline_destroy(o);
  GSHGEOM_Hatch hs = h;
  hs.map = sing;
  ASSERT_EQ(gshgeom_hatch_expand(region, GSHGEOM_FILL_NONZERO, &hs, &st, &chain, 0, nullptr, nullptr, &o), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_stroke_outline_status(o), GSHGEOM_STROKE_NOT_INVERTIBLE) << "a singular hatch map";
  gshgeom_stroke_outline_destroy(o);
  // a style with dashes, or an alignment other than center, is not a hatch's
  std::vector<double> d = {1.0, 1.0};
  GSHGEOM_Stroke bad = st;
  bad.dash_count = 2;
  bad.dashes = d.data();
  EXPECT_EQ(gshgeom_hatch_expand(region, GSHGEOM_FILL_NONZERO, &h, &bad, &chain, 0, nullptr, nullptr, &o), GSHGEOM_ERR_INVALID);
  bad = st;
  bad.alignment = GSHGEOM_STROKE_ALIGN_INSIDE;
  EXPECT_EQ(gshgeom_hatch_expand(region, GSHGEOM_FILL_NONZERO, &h, &bad, &chain, 0, nullptr, nullptr, &o), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_hatch_expand(region, GSHGEOM_FILL_NONZERO, &h, nullptr, &chain, 0, nullptr, nullptr, &o), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_hatch_expand(region, (GSHGEOM_Fill_Rule)7, &h, &st, &chain, 0, nullptr, nullptr, &o), GSHGEOM_ERR_INVALID);
  gshgeom_path_destroy(region);
}

TEST(HatchExpand, ACurvedBoundaryEndsWithinTheToleranceAndAFinerOneMovesCloser) {
  GSHGEOM_Path * circle = nullptr;
  ASSERT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &circle), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_add_ellipse(circle, 20.0, 20.0, 15.0, 15.0, 0, 0), GSHGEOM_OK);
  GSHGEOM_Hatch_Family f;
  GSHGEOM_Hatch h = one_family(f, 90, 0.5, 0, 0, 2.5, {}, identity());
  GSHGEOM_Map chain = identity();
  double worst[2] = {0, 0};
  const int32_t taus[2] = {GSHGEOM_FLATTEN_TOLERANCE, GSHGEOM_FLATTEN_MIN_TOLERANCE};
  for (int k = 0; k < 2; k++) {
    GSHGEOM_Hatch_Pieces * p = nullptr;
    ASSERT_EQ(gshgeom_hatch_pieces_with_tolerance(circle, GSHGEOM_FILL_NONZERO, &h, &chain, taus[k], 0, nullptr, nullptr, &p), GSHGEOM_OK);
    ASSERT_GT(gshgeom_hatch_pieces_count(p), 5u);
    for (size_t i = 0; i < gshgeom_hatch_pieces_count(p); i++) {
      double xy[4];
      int dot;
      ASSERT_EQ(gshgeom_hatch_pieces_get(p, i, xy, &dot), GSHGEOM_OK);
      for (int e = 0; e < 2; e++) {
        const double dx = xy[2 * e] - 20.0, dy = xy[2 * e + 1] - 20.0;
        const double err = std::fabs(std::sqrt(dx * dx + dy * dy) - 15.0);
        worst[k] = std::fmax(worst[k], err);
      }
    }
    gshgeom_hatch_pieces_destroy(p);
  }
  EXPECT_LE(worst[0], taus[0] / 64.0) << "within the tolerance of the true circle";
  EXPECT_LE(worst[1], taus[1] / 64.0);
  EXPECT_LT(worst[1], worst[0]) << "a finer tolerance puts the ends nearer the circle";
  gshgeom_path_destroy(circle);
}

TEST(HatchExpand, RegionsFromAModelPathAndFromFlatContoursAgree) {
  GSHGEOM_Path * region = rect_path(1, 1, 8, 6);
  GSHGEOM_Hatch_Family f;
  GSHGEOM_Hatch h = one_family(f, 30, 0, 0, 0.5, 1.2, {2.0, -1.0}, identity());
  const double m[4] = {2, 0, 0, 2}, t[2] = {3, 4};
  GSHGEOM_Map chain;
  ASSERT_EQ(gshgeom_map_from_matrix(2, m, t, &chain), GSHGEOM_OK);
  GSHGEOM_Flat_Path * flat = nullptr;
  ASSERT_EQ(gshgeom_path_flatten(region, &chain, nullptr, nullptr, &flat), GSHGEOM_OK);
  GSHGEOM_Hatch_Pieces *a = nullptr, *b = nullptr;
  ASSERT_EQ(gshgeom_hatch_pieces(region, GSHGEOM_FILL_NONZERO, &h, &chain, 0, nullptr, nullptr, &a), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_hatch_pieces_from_flat(flat, GSHGEOM_FILL_NONZERO, &h, &chain, 0, nullptr, nullptr, &b), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_hatch_pieces_hash(a), gshgeom_hatch_pieces_hash(b));
  EXPECT_NE(gshgeom_hatch_pieces_count(a), 0u);
  gshgeom_hatch_pieces_destroy(a);
  gshgeom_hatch_pieces_destroy(b);
  gshgeom_flat_path_destroy(flat);
  gshgeom_path_destroy(region);
}

TEST(HatchPieces, TheHashAndTheDumpAreDeterministicAndAccessorsCheck) {
  GSHGEOM_Hatch_Family f;
  GSHGEOM_Hatch h = one_family(f, 45, 0, 0, 0, 2.0, {1.0, -1.0, 0.0, -1.0}, identity());
  GSHGEOM_Path * region = rect_path(0, 0, 9, 9);
  GSHGEOM_Map chain = identity();
  GSHGEOM_Hatch_Pieces *a = nullptr, *b = nullptr;
  ASSERT_EQ(gshgeom_hatch_pieces(region, GSHGEOM_FILL_NONZERO, &h, &chain, 0, nullptr, nullptr, &a), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_hatch_pieces(region, GSHGEOM_FILL_NONZERO, &h, &chain, 0, nullptr, nullptr, &b), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_hatch_pieces_hash(a), gshgeom_hatch_pieces_hash(b));
  EXPECT_NE(gshgeom_hatch_pieces_hash(a), gshgeom_hatch_pieces_hash(nullptr));
  double xy[4];
  int dot;
  EXPECT_EQ(gshgeom_hatch_pieces_get(a, gshgeom_hatch_pieces_count(a), xy, &dot), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_hatch_pieces_get(nullptr, 0, xy, &dot), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_hatch_pieces_count(nullptr), 0u);
  EXPECT_EQ(gshgeom_hatch_pieces_lines(nullptr), 0u);
  EXPECT_EQ(gshgeom_hatch_pieces_dump(nullptr, stdout), GSHGEOM_ERR_INVALID);
  FILE * f1 = tmpfile();
  ASSERT_NE(f1, nullptr);
  EXPECT_EQ(gshgeom_hatch_pieces_dump(a, f1), GSHGEOM_OK);
  rewind(f1);
  char line[256];
  ASSERT_NE(fgets(line, sizeof line, f1), nullptr);
  EXPECT_EQ(std::strncmp(line, "gshgeom.hatch-pieces 1 pieces=", 30), 0) << line;
  fclose(f1);
  EXPECT_EQ(gshgeom_hatch_revision(), 1);
  gshgeom_hatch_pieces_destroy(a);
  gshgeom_hatch_pieces_destroy(b);
  gshgeom_path_destroy(region);
}

// --- the random matrix against the exact oracle -------------------------------------------------------------------------------------

namespace {

struct OFam {
  double angle, ox, oy, dx, dy;
  std::vector<double> el;
};

struct OOut {
  uint64_t lines = 0;
  std::vector<std::array<double, 4>> xy;
  std::vector<int> dot;
  bool near_degenerate = false;
};

ho::RP rp(double x, double y) { return ho::RP{ho::from_double(x), ho::from_double(y)}; }

struct Aff {
  double m[4], t[2];
};

Aff compose(const Aff & a, const Aff & b) {
  Aff r;
  r.m[0] = a.m[0] * b.m[0] + a.m[1] * b.m[2];
  r.m[1] = a.m[0] * b.m[1] + a.m[1] * b.m[3];
  r.m[2] = a.m[2] * b.m[0] + a.m[3] * b.m[2];
  r.m[3] = a.m[2] * b.m[1] + a.m[3] * b.m[3];
  r.t[0] = (a.m[0] * b.t[0] + a.m[1] * b.t[1]) + a.t[0];
  r.t[1] = (a.m[2] * b.t[0] + a.m[3] * b.t[1]) + a.t[1];
  return r;
}

void apply(const Aff & a, double x, double y, double * ox, double * oy) {
  *ox = (a.m[0] * x + a.m[1] * y) + a.t[0];
  *oy = (a.m[2] * x + a.m[3] * y) + a.t[1];
}

/** HX4 to HX6 by the definition, in exact rationals (the library's doubles for the line's defining quantities, HX5). */
OOut oracle(const ho::Contours & cs, int nesting, int rule, const std::vector<OFam> & fams, const Aff & H, const Aff & chain) {
  OOut out;
  const Aff C = compose(chain, H);
  double lox = 1e300, hix = -1e300, loy = 1e300, hiy = -1e300;
  for (const auto & c : cs) {
    for (const ho::Pt & p : c) {
      lox = std::fmin(lox, p.x / 64.0);
      hix = std::fmax(hix, p.x / 64.0);
      loy = std::fmin(loy, p.y / 64.0);
      hiy = std::fmax(hiy, p.y / 64.0);
    }
  }
  const double cx[4] = {lox, hix, hix, lox}, cy[4] = {loy, loy, hiy, hiy};
  for (const OFam & f : fams) {
    if (f.dy == 0.0) continue;
    double ux, uy;
    EXPECT_EQ(gshgeom_sincos(GSHGEOM_ANGLE_DEGREES, f.angle, &ux, &uy), GSHGEOM_OK);
    const double nx = -uy, ny = ux;
    const double wx = C.m[0] * ux + C.m[1] * uy, wy = C.m[2] * ux + C.m[3] * uy;
    const double dvx = H.m[0] * ux + H.m[1] * uy, dvy = H.m[2] * ux + H.m[3] * uy;
    (void)dvx;
    (void)dvy;
    for (int k = -700; k <= 700; k++) {
      const double kdx = (double)k * f.dx, kdy = (double)k * f.dy;
      const double okx = (f.ox + kdx * ux) + kdy * nx, oky = (f.oy + kdx * uy) + kdy * ny;
      double ax, ay;
      apply(C, okx, oky, &ax, &ay);
      const ho::RP A = rp(ax, ay), w = rp(wx, wy);
      int pos = 0, neg = 0;
      for (int q = 0; q < 4; q++) {
        const ho::Rat v = w.x * (ho::from_double(cy[q]) - A.y) - w.y * (ho::from_double(cx[q]) - A.x);
        pos += v.sign() > 0;
        neg += v.sign() < 0;
      }
      if (pos == 4 || neg == 4) continue;
      out.lines++;
      std::vector<ho::Iv> ivs = ho::trim_intervals(cs, nesting, rule, A, w);
      // near-degenerate: a stretch (or a gap) of device length in (0, 2^-24) is outside what the library decides
      const double wl = std::sqrt(wx * wx + wy * wy);
      {
        auto short_pos = [&](const ho::Rat & d) { return d.sign() > 0 && ho::to_double(d) * wl < std::ldexp(1.0, -24); };
        for (size_t i = 0; i < ivs.size(); i++) {
          if (short_pos(ivs[i].b - ivs[i].a)) out.near_degenerate = true;
          if (i + 1 < ivs.size() && short_pos(ivs[i + 1].a - ivs[i].b)) out.near_degenerate = true;
        }
      }
      if (ivs.empty()) continue;
      std::vector<std::pair<ho::Rat, ho::Rat>> ps;
      double sum = 0;
      for (double e : f.el) sum += std::fabs(e);
      bool solid = true;
      if (!f.el.empty()) {
        const ho::Rat span = ho::from_double(sum) * ho::Rat(3);
        std::vector<ho::Iv> ds = ho::dash_set(f.el, ivs.front().a - span, ivs.back().b + span, &solid);
        if (!solid) {
          for (const ho::Iv & d : ds) {
            if (clip_oracle::cmp(d.a, d.b) == 0) {
              // HX6: a dot is kept at an end to within 2^-46 of the larger of the period, 1 and the largest trim position
              const double big = std::fmax(std::fmax(std::fabs(ho::to_double(ivs.front().a)), std::fabs(ho::to_double(ivs.back().b))), 1.0);
              const ho::Rat slack = ho::pow2(-46) * ho::from_double(std::fmax(sum, big));
              for (const ho::Iv & iv : ivs) {
                if (clip_oracle::cmp(iv.a - slack, d.a) <= 0 && clip_oracle::cmp(d.a, iv.b + slack) <= 0) {
                  const ho::Rat at = clip_oracle::cmp(d.a, iv.a) < 0 ? iv.a : (clip_oracle::cmp(d.a, iv.b) > 0 ? iv.b : d.a);
                  ps.emplace_back(at, at);
                  break;
                }
              }
            } else {
              for (const ho::Iv & iv : ivs) {
                const ho::Rat lo = clip_oracle::cmp(d.a, iv.a) > 0 ? d.a : iv.a;
                const ho::Rat hi = clip_oracle::cmp(d.b, iv.b) < 0 ? d.b : iv.b;
                const ho::Rat len = hi - lo;
                if (len.sign() > 0) {
                  if (ho::to_double(len) * wl < std::ldexp(1.0, -24)) out.near_degenerate = true;
                  ps.emplace_back(lo, hi);
                }
              }
            }
          }
        }
      }
      if (f.el.empty() || solid) {
        for (const ho::Iv & iv : ivs) ps.emplace_back(iv.a, iv.b);
      }
      for (const auto & pr : ps) {
        const bool dot = clip_oracle::cmp(pr.first, pr.second) == 0;
        const double s0 = ho::to_double(pr.first), s1 = ho::to_double(pr.second);
        double x0, y0, x1, y1;
        apply(H, okx + s0 * ux, oky + s0 * uy, &x0, &y0);
        apply(H, okx + s1 * ux, oky + s1 * uy, &x1, &y1);
        out.xy.push_back({x0, y0, dot ? x0 : x1, dot ? y0 : y1});
        out.dot.push_back(dot ? 1 : 0);
      }
    }
  }
  return out;
}

}  // namespace

namespace {

GSHGEOM_Map map_of(const Aff & a) {
  GSHGEOM_Map m;
  EXPECT_EQ(gshgeom_map_from_matrix(2, a.m, a.t, &m), GSHGEOM_OK);
  return m;
}

const Aff kMaps[] = {
    {{1, 0, 0, 1}, {0, 0}},   {{2, 0, 0, 2}, {5, 3}},   {{1, 0, 0, -1}, {0, 20}}, {{0, -1, 1, 0}, {20, 0}},
    {{2, 0, 0, 1}, {0, 0}},   {{0.5, 0, 0, 0.5}, {1, 1}}, {{1, 0, 0, 1}, {2.25, -3.5}}, {{0, 1, 1, 0}, {0, 0}},
};

ho::Contours random_region(Rng & r, int * rule) {
  ho::Contours cs;
  auto rect = [&](int64_t x0, int64_t y0, int64_t x1, int64_t y1, bool ccw) {
    ho::Contour c = {{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}};
    if (!ccw) std::reverse(c.begin(), c.end());
    return c;
  };
  auto q = [&](int64_t lo, int64_t hi) { return r.range(lo, hi) * 16; };  // quarter pixels
  *rule = (int)r.range(0, 1);
  switch (r.range(0, 5)) {
    case 0:
      cs.push_back(rect(q(0, 8), q(0, 8), q(40, 80), q(40, 80), r.range(0, 1)));
      break;
    case 1: {  // an L
      const int64_t x1 = q(40, 80), y1 = q(40, 80), xm = q(16, 30), ym = q(16, 30);
      cs.push_back({{0, 0}, {x1, 0}, {x1, ym}, {xm, ym}, {xm, y1}, {0, y1}});
      break;
    }
    case 2: {  // a rect with a hole, the hole the same or the other way round
      const int64_t x1 = q(60, 90), y1 = q(60, 90);
      cs.push_back(rect(0, 0, x1, y1, true));
      cs.push_back(rect(q(8, 20), q(8, 20), x1 - q(8, 20), y1 - q(8, 20), r.range(0, 1)));
      break;
    }
    case 3: {  // three levels
      cs.push_back(rect(0, 0, q(80, 90), q(80, 90), true));
      cs.push_back(rect(q(8, 14), q(8, 14), q(60, 70), q(60, 70), r.range(0, 1)));
      cs.push_back(rect(q(20, 26), q(20, 26), q(40, 50), q(40, 50), r.range(0, 1)));
      break;
    }
    case 4: {  // a star-ish polygon
      ho::Contour c;
      const int n = (int)r.range(5, 9);
      for (int i = 0; i < n; i++) c.push_back({q(0, 80), q(0, 80)});
      cs.push_back(c);
      break;
    }
    default:  // two apart, or touching
      cs.push_back(rect(0, 0, q(20, 30), q(20, 30), true));
      cs.push_back(rect(q(24, 36), q(10, 20), q(50, 60), q(30, 40), r.range(0, 1)));
      break;
  }
  return cs;
}

}  // namespace

TEST(HatchMatrix, RandomHatchesMatchTheExactOracle) {
  GSHGEOM_SKIP_UNDER_VALGRIND();
  Rng r{UINT64_C(0x9e3779b97f4a7c15)};
  const double angles[] = {0, 90, 180, 270, 30, 60, 45, 135, 17};
  const std::vector<std::vector<double>> patterns = {{}, {2, -1}, {1, -0.5, 0, -0.5}, {0, -1.25}, {-1, 0}, {0.5, -0.5, 1, -2}, {3, 0, -2}, {0, -1, 2}, {1.5}};
  const double dxs[] = {0, 0.5, -1.25, 0.75};
  const double dys[] = {0.75, 1, 1.5, 2, 2.5, 3, -1.25};
  int compared = 0, skipped = 0, with_dots = 0, with_dashes = 0, nested = 0;
  for (int iter = 0; iter < 260; iter++) {
    int rule = 0;
    const ho::Contours cs = random_region(r, &rule);
    const int nesting = iter % 5 == 4 ? (int)r.range(1, 2) : 0;
    std::vector<OFam> fams;
    const int nf = r.range(0, 3) == 0 ? 2 : 1;
    for (int i = 0; i < nf; i++) {
      OFam f;
      f.angle = angles[r.range(0, 8)];
      f.ox = (double)r.range(-8, 8) / 4.0;
      f.oy = (double)r.range(-8, 8) / 4.0;
      f.dx = dxs[r.range(0, 3)];
      f.dy = dys[r.range(0, 6)];
      f.el = patterns[(size_t)r.range(0, (int64_t)patterns.size() - 1)];
      fams.push_back(f);
    }
    const Aff & chain = kMaps[r.range(0, 7)];
    const Aff & H = kMaps[r.range(0, 7)];
    const OOut o = oracle(cs, nesting, rule, fams, H, chain);
    if (o.near_degenerate) {
      skipped++;
      continue;
    }
    // the library
    std::vector<std::vector<GSHGEOM_Device_Point>> lc;
    for (const auto & c : cs) {
      std::vector<GSHGEOM_Device_Point> one;
      for (const ho::Pt & p : c) one.push_back({(int32_t)p.x, (int32_t)p.y});
      lc.push_back(one);
    }
    Flat region(lc);
    std::vector<GSHGEOM_Hatch_Family> lf;
    for (const OFam & f : fams) {
      GSHGEOM_Hatch_Family g;
      g.angle = f.angle;
      g.origin_x = f.ox;
      g.origin_y = f.oy;
      g.offset_along = f.dx;
      g.offset_across = f.dy;
      g.dash_count = f.el.size();
      g.dashes = f.el.empty() ? nullptr : f.el.data();
      lf.push_back(g);
    }
    GSHGEOM_Hatch h;
    std::memset(&h, 0, sizeof h);
    h.angle_unit = GSHGEOM_ANGLE_DEGREES;
    h.nesting = (uint32_t)nesting;
    h.map = map_of(H);
    h.family_count = lf.size();
    h.families = lf.data();
    GSHGEOM_Map cm = map_of(chain);
    GSHGEOM_Hatch_Pieces * p = nullptr;
    ASSERT_EQ(gshgeom_hatch_pieces_from_flat(region.p, rule ? GSHGEOM_FILL_EVEN_ODD : GSHGEOM_FILL_NONZERO, &h, &cm, 1000, nullptr, nullptr, &p),
              GSHGEOM_OK) << iter;
    EXPECT_EQ(gshgeom_hatch_pieces_lines(p), o.lines) << "iteration " << iter;
    if (gshgeom_hatch_pieces_count(p) != o.xy.size()) {
      std::fprintf(stderr, "iteration %d nesting %d rule %d fams:", iter, nesting, rule);
      for (const OFam & f : fams) {
        std::fprintf(stderr, " [a%g o(%g,%g) dx%g dy%g el%zu]", f.angle, f.ox, f.oy, f.dx, f.dy, f.el.size());
      }
      std::fprintf(stderr, " chain m %g %g %g %g t %g %g H m %g %g %g %g t %g %g\ncontours:", chain.m[0], chain.m[1], chain.m[2], chain.m[3], chain.t[0], chain.t[1], H.m[0], H.m[1], H.m[2], H.m[3], H.t[0], H.t[1]);
      for (const auto & c : cs) {
        std::fprintf(stderr, " {");
        for (const ho::Pt & q : c) std::fprintf(stderr, "(%g,%g)", q.x / 64.0, q.y / 64.0);
        std::fprintf(stderr, "}");
      }
      std::fprintf(stderr, "\n");
      for (size_t i = 0; i < std::max(o.xy.size(), gshgeom_hatch_pieces_count(p)); i++) {
        double xy[4] = {0, 0, 0, 0};
        int dot = -1;
        if (i < gshgeom_hatch_pieces_count(p)) (void)gshgeom_hatch_pieces_get(p, i, xy, &dot);
        std::fprintf(stderr, "%3zu lib %d (%.17g,%.17g)-(%.17g,%.17g)", i, dot, xy[0], xy[1], xy[2], xy[3]);
        if (i < o.xy.size()) std::fprintf(stderr, "   orc %d (%.17g,%.17g)-(%.17g,%.17g)", o.dot[i], o.xy[i][0], o.xy[i][1], o.xy[i][2], o.xy[i][3]);
        std::fprintf(stderr, "\n");
      }
    }
    ASSERT_EQ(gshgeom_hatch_pieces_count(p), o.xy.size()) << "iteration " << iter;
    for (size_t i = 0; i < o.xy.size(); i++) {
      double xy[4];
      int dot;
      ASSERT_EQ(gshgeom_hatch_pieces_get(p, i, xy, &dot), GSHGEOM_OK);
      EXPECT_EQ(dot, o.dot[i]) << "iteration " << iter << " piece " << i;
      for (int q = 0; q < 4; q++) EXPECT_TRUE(near(xy[q], o.xy[i][q])) << "iteration " << iter << " piece " << i << " coord " << q << ": " << xy[q] << " vs " << o.xy[i][q];
      with_dots += dot;
    }
    with_dashes += !lf[0].dash_count ? 0 : 1;
    nested += nesting != 0;
    gshgeom_hatch_pieces_destroy(p);
    compared++;
  }
  EXPECT_GE(compared, 200) << "skipped near-degenerate: " << skipped;
  EXPECT_GT(with_dots, 20);
  EXPECT_GT(with_dashes, 50);
  EXPECT_GT(nested, 20);
}

// --- the outline against the definition of the stroke of the pieces (SG17) ---------------------------------------------------------------

TEST(HatchMatrix, TheOutlineCoversTheStrokeOfThePiecesAtEveryDecidedPoint) {
  Rng r{UINT64_C(0xd1b54a32d192ed03)};
  const uint32_t caps[] = {GSHGEOM_STROKE_CAP_BUTT, GSHGEOM_STROKE_CAP_ROUND, GSHGEOM_STROKE_CAP_SQUARE};
  const std::vector<std::vector<double>> patterns = {{}, {2, -1}, {0, -1.5}, {1, -0.5, 0, -0.5}};
  int decided = 0;
  for (int iter = 0; iter < 24; iter++) {
    GSHGEOM_Hatch_Family f;
    static const double kAngles[4] = {0, 90, 45, 30};
    const double angle = kAngles[iter % 4];
    GSHGEOM_Hatch h = one_family(f, angle, 0.25, 0.5, iter % 3 ? 0.0 : 0.5, 1.5 + (double)(iter % 3) * 0.25, patterns[(size_t)(iter / 4) % 4], identity());
    GSHGEOM_Path * region = rect_path(1.0, 1.0, 9.0, 8.0);
    GSHGEOM_Map chain = map_of(kMaps[iter % 2 == 0 ? 0 : 1]);
    const uint32_t cap = caps[iter % 3];
    GSHGEOM_Stroke st = style(0.6, cap);
    GSHGEOM_Hatch_Pieces * p = nullptr;
    ASSERT_EQ(gshgeom_hatch_pieces(region, GSHGEOM_FILL_NONZERO, &h, &chain, 0, nullptr, nullptr, &p), GSHGEOM_OK);
    GSHGEOM_Stroke_Outline * o = nullptr;
    ASSERT_EQ(gshgeom_hatch_expand_pieces(p, &st, &chain, GSHGEOM_FLATTEN_TOLERANCE, nullptr, nullptr, &o), GSHGEOM_OK);
    const GSHGEOM_Flat_Path * body = gshgeom_stroke_outline_body(o);
    ASSERT_NE(body, nullptr);
    const double hh = 0.3;
    double ux, uy;
    ASSERT_EQ(gshgeom_sincos(GSHGEOM_ANGLE_DEGREES, angle, &ux, &uy), GSHGEOM_OK);
    const double scale = chain.m[0][0];
    for (int t = 0; t < 400; t++) {
      // a random point in local space; the pixel point is its image
      const double lx = (double)r.range(0, 1200) / 100.0, ly = (double)r.range(0, 1200) / 100.0;
      // signed distance (local) to the nearest piece shape: negative inside; decided only outside a 0.12 pixel band
      double best = 1e9;
      for (size_t i = 0; i < gshgeom_hatch_pieces_count(p); i++) {
        double xy[4];
        int dot;
        ASSERT_EQ(gshgeom_hatch_pieces_get(p, i, xy, &dot), GSHGEOM_OK);
        double dxv = xy[2] - xy[0], dyv = xy[3] - xy[1];
        double L = std::sqrt(dxv * dxv + dyv * dyv);
        double ex, ey;
        if (dot) {
          ex = ux;
          ey = uy;
          L = 0;
        } else {
          ex = dxv / L;
          ey = dyv / L;
        }
        const double tt = (lx - xy[0]) * ex + (ly - xy[1]) * ey;
        const double vv = std::fabs((lx - xy[0]) * ey - (ly - xy[1]) * ex);
        double d;  // distance outside the shape (negative inside)
        if (cap == GSHGEOM_STROKE_CAP_BUTT) {
          if (dot) continue;
          d = std::fmax(std::fmax(-tt, tt - L), vv - hh);
        } else if (cap == GSHGEOM_STROKE_CAP_SQUARE) {
          d = std::fmax(std::fmax(-tt - hh, tt - L - hh), vv - hh);
        } else {
          const double cl = tt < 0 ? 0 : (tt > L ? L : tt);
          const double rx = (lx - xy[0]) - cl * ex, ry = (ly - xy[1]) - cl * ey;
          d = std::sqrt(rx * rx + ry * ry) - hh;
        }
        best = std::fmin(best, d);
      }
      const double band = 0.12 / scale;
      if (std::fabs(best) < band) continue;
      double px_, py_;
      px_ = chain.m[0][0] * lx + chain.m[0][1] * ly + chain.t[0];
      py_ = chain.m[1][0] * lx + chain.m[1][1] * ly + chain.t[1];
      const bool in_oracle = best < 0;
      const bool in_outline = winding(body, px_, py_) != 0;
      EXPECT_EQ(in_outline, in_oracle) << "iteration " << iter << " point " << lx << ", " << ly << " distance " << best;
      decided++;
    }
    gshgeom_stroke_outline_destroy(o);
    gshgeom_hatch_pieces_destroy(p);
    gshgeom_path_destroy(region);
  }
  EXPECT_GT(decided, 3000);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

// --- the rules the planted defects pinned down (check-planted.py, hx- plants) ------------------------------------------------------

namespace {

using DP = GSHGEOM_Device_Point;

DP dp(double x, double y) {
  DP d;
  d.x = (int32_t)std::lround(x * 64.0);
  d.y = (int32_t)std::lround(y * 64.0);
  return d;
}

GSHGEOM_Map scale_map(double s) {
  GSHGEOM_Map m;
  EXPECT_EQ(gshgeom_map_scale(2, s, s, 1.0, &m), GSHGEOM_OK);
  return m;
}

std::vector<std::array<double, 4>> pieces_of(const GSHGEOM_Hatch_Pieces * p, std::vector<int> * dots = nullptr) {
  std::vector<std::array<double, 4>> out;
  for (size_t i = 0; i < gshgeom_hatch_pieces_count(p); i++) {
    double xy[4];
    int dot = 0;
    EXPECT_EQ(gshgeom_hatch_pieces_get(p, i, xy, &dot), GSHGEOM_OK);
    out.push_back({xy[0], xy[1], xy[2], xy[3]});
    if (dots) dots->push_back(dot);
  }
  return out;
}

}  // namespace

TEST(HatchCheck, AnInvalidValueIsDecidedBeforeTheCountsAndTheAngleBoundIsInclusive) {
  GSHGEOM_Hatch_Family f;
  GSHGEOM_Hatch h = one_family(f, 45, 0, 0, 0, 1, {1.0, -1.0}, identity());
  std::vector<GSHGEOM_Hatch_Family> many(GSHGEOM_HATCH_FAMILIES_MAX + 1, f);
  GSHGEOM_Hatch hf = h;
  hf.families = many.data();
  hf.family_count = many.size();
  EXPECT_EQ(gshgeom_hatch_check(&hf), GSHGEOM_ERR_LIMIT);
  hf.map.m[0][0] = std::nan("");
  EXPECT_EQ(gshgeom_hatch_check(&hf), GSHGEOM_ERR_INVALID) << "an invalid map beats too many families";
  ASSERT_EQ(gshgeom_map_identity(3, &hf.map), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_hatch_check(&hf), GSHGEOM_ERR_UNSUPPORTED) << "so does a 3D map";
  // the angle bound is inclusive in every unit that can reach it
  GSHGEOM_Hatch_Family g = f;
  GSHGEOM_Hatch hr = h;
  hr.families = &g;
  hr.angle_unit = GSHGEOM_ANGLE_RADIANS;
  g.angle = GSHGEOM_HATCH_ANGLE_MAX;
  EXPECT_EQ(gshgeom_hatch_check(&hr), GSHGEOM_OK);
  g.angle = std::nextafter(GSHGEOM_HATCH_ANGLE_MAX, 1e300);
  EXPECT_EQ(gshgeom_hatch_check(&hr), GSHGEOM_ERR_LIMIT);
  g.angle = -GSHGEOM_HATCH_ANGLE_MAX;
  EXPECT_EQ(gshgeom_hatch_check(&hr), GSHGEOM_OK);
  // an angle whose quarter turns pass what ::gshgeom_sincos takes is a limit although the angle is in its own band
  hr.angle_unit = GSHGEOM_ANGLE_TURNS;
  g.angle = std::ldexp(1.0, 39);
  EXPECT_EQ(gshgeom_hatch_check(&hr), GSHGEOM_ERR_LIMIT);
  g.angle = std::ldexp(1.0, 38);
  EXPECT_EQ(gshgeom_hatch_check(&hr), GSHGEOM_OK);
  hr.angle_unit = GSHGEOM_ANGLE_DEGREES;
  g.angle = GSHGEOM_HATCH_ANGLE_MAX;
  EXPECT_EQ(gshgeom_hatch_check(&hr), GSHGEOM_OK);
}

TEST(HatchLines, ALineOnTheBoxEdgeIsCountedWhereTheDivisionRoundsItAway) {
  // spacing 0.1 from 0.1: line 319 is at y = 32 exactly (it touches the top of the box) and (32 - 0.1) / 0.1 is 318.99999999999994
  GSHGEOM_Hatch_Family f;
  GSHGEOM_Hatch h = one_family(f, 0, 0, 0.1, 0, 0.1, {}, identity());
  GSHGEOM_Path * region = rect_path(0, 0, 32, 32);
  GSHGEOM_Map chain = identity();
  GSHGEOM_Hatch_Pieces * p = nullptr;
  ASSERT_EQ(gshgeom_hatch_pieces(region, GSHGEOM_FILL_NONZERO, &h, &chain, 0, nullptr, nullptr, &p), GSHGEOM_OK);
  uint64_t expected = 0;
  for (int k = -5; k <= 330; k++) {
    const double y = (0.1 + 0.0) + ((double)k * 0.1) * 1.0;
    expected += y >= 0.0 && y <= 32.0;
  }
  EXPECT_GE(expected, 320u);
  EXPECT_EQ(gshgeom_hatch_pieces_lines(p), expected);
  gshgeom_hatch_pieces_destroy(p);
  gshgeom_path_destroy(region);
}

namespace {

/** The library's pieces against the exact oracle for one region, one family, and the maps. */
void expect_matches_oracle(const ho::Contours & cs, int rule, const OFam & fam, const Aff & H, const Aff & chain, const char * what) {
  const OOut o = oracle(cs, 0, rule, {fam}, H, chain);
  ASSERT_FALSE(o.near_degenerate) << what;
  std::vector<std::vector<DP>> lc;
  for (const auto & c : cs) {
    std::vector<DP> one;
    for (const ho::Pt & p : c) one.push_back({(int32_t)p.x, (int32_t)p.y});
    lc.push_back(one);
  }
  Flat region(lc);
  GSHGEOM_Hatch_Family g;
  g.angle = fam.angle;
  g.origin_x = fam.ox;
  g.origin_y = fam.oy;
  g.offset_along = fam.dx;
  g.offset_across = fam.dy;
  g.dash_count = fam.el.size();
  g.dashes = fam.el.empty() ? nullptr : fam.el.data();
  GSHGEOM_Hatch h;
  std::memset(&h, 0, sizeof h);
  h.angle_unit = GSHGEOM_ANGLE_DEGREES;
  h.map = map_of(H);
  h.family_count = 1;
  h.families = &g;
  GSHGEOM_Map cm = map_of(chain);
  GSHGEOM_Hatch_Pieces * p = nullptr;
  ASSERT_EQ(gshgeom_hatch_pieces_from_flat(region.p, rule ? GSHGEOM_FILL_EVEN_ODD : GSHGEOM_FILL_NONZERO, &h, &cm, 1000, nullptr, nullptr, &p),
            GSHGEOM_OK) << what;
  EXPECT_EQ(gshgeom_hatch_pieces_lines(p), o.lines) << what;
  const auto got = pieces_of(p);
  ASSERT_EQ(got.size(), o.xy.size()) << what;
  for (size_t i = 0; i < got.size(); i++) {
    for (int q = 0; q < 4; q++) EXPECT_TRUE(near(got[i][q], o.xy[i][q])) << what << " piece " << i << " coord " << q << ": " << got[i][q] << " vs " << o.xy[i][q];
  }
  gshgeom_hatch_pieces_destroy(p);
}

}  // namespace

TEST(HatchLines, AVertexAHairFromTheLineIsOnTheSideTheExactSumSays) {
  // a 45 degree line whose direction (cos, sin) differs in the last place runs along the diagonal of a triangle within 1e-15: the
  // double side value of the far corner has the wrong sign for these origins, the exact sum has the right one, and a wrong side
  // changes the piece by a crossing in the middle of the diagonal
  const Aff id = {{1, 0, 0, 1}, {0, 0}};
  // both corners of the diagonal are on the same side of the line exactly, so the edge is not crossed; the double value of the
  // far corner is the other sign for each of these
  const double origins[][2] = {{-9.95821614022951e-16, -2.658756415791537e-15}, {2.60150892154355e-15, 8.914977737327084e-16},
                               {8.298038105269666e-16, -7.869417060953046e-16}};
  int n = 0;
  for (const auto & o : origins) {
    for (int mirror = 0; mirror < 2; mirror++) {
      ho::Contours cs;
      if (mirror)
        cs.push_back({{0, 0}, {640, 640}, {640, 0}});
      else
        cs.push_back({{0, 0}, {640, 640}, {0, 640}});
      OFam f;
      f.angle = 45;
      f.ox = o[0];
      f.oy = o[1];
      f.dx = 0;
      f.dy = 1000;
      expect_matches_oracle(cs, 0, f, id, id, ("near-diagonal origin " + std::to_string(n++)).c_str());
    }
  }
}

namespace {

/** The pieces of one hatch line of a region, with the hatch at @p scale (the chain) and everything in device pixels as given. */
std::vector<std::array<double, 4>> pieces_at_scale(const std::vector<std::vector<DP>> & contours, double scale, double oy_device, double dy_device,
                                                   double angle, const std::vector<double> & dashes, uint64_t * lines) {
  Flat region(contours);
  GSHGEOM_Hatch_Family f;
  GSHGEOM_Hatch h = one_family(f, angle, 0, oy_device / scale, 0, dy_device / scale, dashes, identity());
  GSHGEOM_Map chain = scale_map(scale);
  GSHGEOM_Hatch_Pieces * p = nullptr;
  EXPECT_EQ(gshgeom_hatch_pieces_from_flat(region.p, GSHGEOM_FILL_NONZERO, &h, &chain, 0, nullptr, nullptr, &p), GSHGEOM_OK);
  if (lines) *lines = gshgeom_hatch_pieces_lines(p);
  auto out = pieces_of(p);
  gshgeom_hatch_pieces_destroy(p);
  return out;
}

}  // namespace

TEST(HatchBand, AnIntervalShorterThanTheBandIsDroppedAndOneLongerIsKept) {
  // a triangle with its apex at (50, 10): the line at depth delta below the apex crosses it in a chord 10 delta wide
  const std::vector<std::vector<DP>> tri = {{dp(0, 0), dp(100, 0), dp(50, 10)}};
  for (const double delta : {1e-9, 2e-9}) {
    uint64_t lines = 0;
    const auto pcs = pieces_at_scale(tri, 1.0, 10.0 - delta, 1.0, 0, {}, &lines);
    EXPECT_EQ(lines, 10u) << delta;
    const bool kept = delta > 1.5e-9;
    ASSERT_EQ(pcs.size(), kept ? 10u : 9u) << delta;
    if (kept) {
      const auto top = std::max_element(pcs.begin(), pcs.end(), [](const auto & a, const auto & b) { return a[1] < b[1]; });
      EXPECT_NEAR((*top)[2] - (*top)[0], 10.0 * delta, 1e-12);
    }
  }
}

TEST(HatchBand, CrossingsEachWithinTheBandOfTheLastNeverCutAnIntervalTheyLieIn) {
  // an outer square with three nested triangles, alternating in orientation, whose apexes meet on the point (0, 0): a line at depth
  // d below it crosses six edges at +-1.2 d, 0.8 d and 0.4 d, the winding goes 1 0 1 0 1 0 1 and the line is inside the region at
  // both ends.  At scale 1 and 1/4 the six are one position, so the line is one piece; at scale 4 and d = 5e-8 the steps are
  // beyond the band (in pixels), each is a position of its own, and the pieces are the four members of the pattern.
  struct Case {
    double scale, d;
    size_t pieces;
  };
  for (const Case c : {Case{1.0, 1e-8, 1}, Case{0.25, 1e-8, 1}, Case{4.0, 5e-8, 4}}) {
    std::vector<std::vector<DP>> cs;
    cs.push_back({dp(-50, -50), dp(50, -50), dp(50, 50), dp(-50, 50)});
    const double hw[3] = {6, 4, 2};
    for (int i = 0; i < 3; i++) {
      if (i % 2 == 0)
        cs.push_back({dp(0, 0), dp(hw[i], -5), dp(-hw[i], -5)});
      else
        cs.push_back({dp(0, 0), dp(-hw[i], -5), dp(hw[i], -5)});
    }
    uint64_t lines = 0;
    const auto pcs = pieces_at_scale(cs, c.scale, -c.d, 200.0, 0, {}, &lines);
    EXPECT_EQ(lines, 1u) << c.scale;
    ASSERT_EQ(pcs.size(), c.pieces) << "scale " << c.scale;
    EXPECT_NEAR(pcs.front()[0], -50.0 / c.scale, 1e-9);
    EXPECT_NEAR(pcs.back()[2], 50.0 / c.scale, 1e-9);
  }
}

TEST(HatchBand, AnIntervalBeginsAtTheLastCrossingOfTheGroupThatEntersTheRegion) {
  // three triangles of one orientation with a shared apex on (0, 0) and their lower edges at cot 0.4, 0.8 and 1.2: the line at depth
  // d meets them at 0.4 d, 0.8 d, 1.2 d (the winding goes 1, 2, 3) and the three right edges at x = 24.  When the three are one
  // position (a step within the band) the interval begins at the last, 1.2 d; when each is a position of its own it begins at the first.
  struct Case {
    double scale, d, x0;
  };
  for (const Case c : {Case{1.0, 1e-8, 1.2e-8}, Case{0.25, 1e-8, 1.2e-8 / 0.25}, Case{4.0, 5e-8, 0.4 * 5e-8 / 4.0}}) {
    std::vector<std::vector<DP>> cs;
    const double Y[3] = {20, 30, 60};
    for (int i = 0; i < 3; i++) cs.push_back({dp(0, 0), dp(24, -Y[i]), dp(24, 10)});
    uint64_t lines = 0;
    const auto pcs = pieces_at_scale(cs, c.scale, -c.d, 200.0, 0, {}, &lines);
    EXPECT_EQ(lines, 1u) << c.scale;
    ASSERT_EQ(pcs.size(), 1u) << "scale " << c.scale;
    EXPECT_NEAR(pcs[0][0], c.x0, 1e-12) << "scale " << c.scale;
    EXPECT_NEAR(pcs[0][2], 24.0 / c.scale, 1e-9);
  }
}

TEST(HatchDashes, ASliverOfADashAtTheEndOfAnIntervalIsDropped) {
  // a line 5e-9 short of the pattern: the last dash [6, 7] meets the interval [5e-9, 6 + 5e-9] in 5e-9, less than the band
  const std::vector<std::vector<DP>> rect = {{dp(0, 0), dp(6, 0), dp(6, 10), dp(0, 10)}};
  Flat region(rect);
  GSHGEOM_Hatch_Family f;
  GSHGEOM_Hatch h = one_family(f, 0, -5e-9, 5, 0, 100, {1.0, -1.0}, identity());
  GSHGEOM_Map chain = identity();
  GSHGEOM_Hatch_Pieces * p = nullptr;
  ASSERT_EQ(gshgeom_hatch_pieces_from_flat(region.p, GSHGEOM_FILL_NONZERO, &h, &chain, 0, nullptr, nullptr, &p), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_hatch_pieces_lines(p), 1u);
  const auto pcs = pieces_of(p);
  ASSERT_EQ(pcs.size(), 3u);
  EXPECT_NEAR(pcs[0][0], 0.0, 1e-12);
  EXPECT_NEAR(pcs[0][2], 1.0 - 5e-9, 1e-12);
  EXPECT_NEAR(pcs[2][2], 5.0 - 5e-9, 1e-12);
  gshgeom_hatch_pieces_destroy(p);
}

TEST(HatchDashes, ADotJustOutsideTheIntervalIsKeptAndPutOnItsEnd) {
  // dots at s = 0, 5, 10 ...; the interval is [-1e-14, 5 - 1e-14], so the dot at 5 is 1e-14 past its end, inside the slack
  const std::vector<std::vector<DP>> rect = {{dp(0, 0), dp(5, 0), dp(5, 10), dp(0, 10)}};
  Flat region(rect);
  GSHGEOM_Hatch_Family f;
  GSHGEOM_Hatch h = one_family(f, 0, 1e-14, 5, 0, 100, {0.0, -5.0}, identity());
  GSHGEOM_Map chain = identity();
  GSHGEOM_Hatch_Pieces * p = nullptr;
  ASSERT_EQ(gshgeom_hatch_pieces_from_flat(region.p, GSHGEOM_FILL_NONZERO, &h, &chain, 0, nullptr, nullptr, &p), GSHGEOM_OK);
  std::vector<int> dots;
  const auto pcs = pieces_of(p, &dots);
  ASSERT_EQ(pcs.size(), 2u);
  EXPECT_EQ(dots[0], 1);
  EXPECT_EQ(dots[1], 1);
  EXPECT_NEAR(pcs[0][0], 1e-14, 2e-15);
  EXPECT_LE(pcs[1][0], 5.0 + 2e-15) << "clamped onto the interval's end, not left 1e-14 past it";
  gshgeom_hatch_pieces_destroy(p);
}

TEST(HatchDashes, ALeadingSpaceShiftsThePatternByItsOffset) {
  const std::vector<std::vector<DP>> rect = {{dp(0, 0), dp(10, 0), dp(10, 10), dp(0, 10)}};
  for (const std::vector<double> & pat : std::vector<std::vector<double>>{{-1.0, 2.0}, {-1.0, 2.0, -1.0, 0.0}, {0.0, -3.0, 1.0}}) {
    Flat region(rect);
    GSHGEOM_Hatch_Family f;
    GSHGEOM_Hatch h = one_family(f, 0, 0, 5, 0, 100, pat, identity());
    GSHGEOM_Map chain = identity();
    GSHGEOM_Hatch_Pieces * p = nullptr;
    ASSERT_EQ(gshgeom_hatch_pieces_from_flat(region.p, GSHGEOM_FILL_NONZERO, &h, &chain, 0, nullptr, nullptr, &p), GSHGEOM_OK);
    const auto got = pieces_of(p);
    // the pattern repeats from s = 0: the elements' running sum gives the dashes, the dots and the spaces
    double period = 0;
    for (double e : pat) period += std::fabs(e);
    std::vector<std::pair<double, double>> want;
    for (int rep = -1; rep <= 11; rep++) {
      double at = rep * period;
      for (double e : pat) {
        if (e >= 0) {
          const double a = std::fmax(at, 0.0), b = std::fmin(at + e, 10.0);
          if (e == 0 ? (at >= 0 && at <= 10) : b > a) want.push_back({e == 0 ? at : a, e == 0 ? at : b});
        }
        at += std::fabs(e);
      }
    }
    ASSERT_EQ(got.size(), want.size()) << pat.size();
    for (size_t i = 0; i < got.size(); i++) {
      EXPECT_NEAR(got[i][0], want[i].first, 1e-12) << i;
      EXPECT_NEAR(got[i][2], want[i].second, 1e-12) << i;
    }
    gshgeom_hatch_pieces_destroy(p);
  }
}

TEST(HatchDashes, APatternThatLowersToSolidIsTheSolidFamilyBitForBit) {
  Rng r{UINT64_C(0x1234567)};
  for (int i = 0; i < 48; i++) {
    const double ox = (double)r.range(1, 999) / 997.0, oy = (double)r.range(1, 999) / 1009.0;
    const std::vector<std::vector<DP>> cs = {{dp(3.015625, 2.5), dp(21.2, 3.0), dp(17.75, 15.5), dp(2.5, 11.0)}, {dp(9, 6), dp(12, 6), dp(12, 9), dp(9, 9)}};
    uint64_t hashes[3];
    int k = 0;
    for (const std::vector<double> & pat : std::vector<std::vector<double>>{{}, {1.5}, {1.0, 0.5}}) {
      Flat region(cs);
      GSHGEOM_Hatch_Family f;
      GSHGEOM_Hatch h = one_family(f, 30, ox, oy, 0.37, 1.13, pat, identity());
      GSHGEOM_Map chain = identity();
      GSHGEOM_Hatch_Pieces * p = nullptr;
      ASSERT_EQ(gshgeom_hatch_pieces_from_flat(region.p, GSHGEOM_FILL_EVEN_ODD, &h, &chain, 0, nullptr, nullptr, &p), GSHGEOM_OK);
      hashes[k++] = gshgeom_hatch_pieces_hash(p);
      gshgeom_hatch_pieces_destroy(p);
    }
    EXPECT_EQ(hashes[0], hashes[1]) << i;
    EXPECT_EQ(hashes[0], hashes[2]) << i;
  }
}

TEST(HatchLimits, ThePiecesAndTheDashStepsAreBoundedByTheLimits) {
  GSHGEOM_Hatch_Family f;
  GSHGEOM_Hatch h = one_family(f, 0, 0, 0, 0, 1.0, {}, identity());
  GSHGEOM_Path * region = rect_path(0, 0, 10, 10);  // eleven lines, eleven pieces
  GSHGEOM_Map chain = identity();
  GSHGEOM_Limits limits;
  gshgeom_limits_default(&limits);
  GSHGEOM_Hatch_Pieces * p = nullptr;
  limits.max_elements = 22;  // half of it is the most pieces
  ASSERT_EQ(gshgeom_hatch_pieces(region, GSHGEOM_FILL_NONZERO, &h, &chain, 0, nullptr, &limits, &p), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_hatch_pieces_count(p), 11u);
  gshgeom_hatch_pieces_destroy(p);
  for (size_t cap : {size_t(21), size_t(20)}) {
    limits.max_elements = cap;
    p = nullptr;
    EXPECT_EQ(gshgeom_hatch_pieces(region, GSHGEOM_FILL_NONZERO, &h, &chain, 0, nullptr, &limits, &p), GSHGEOM_ERR_LIMIT) << cap;
    EXPECT_EQ(p, nullptr);
  }
  // a dashed line of a hundred pixels takes hundreds of steps: a cap of four steps is a limit and the default is not
  GSHGEOM_Hatch_Family g;
  GSHGEOM_Hatch hd = one_family(g, 0, 0, 5, 0, 100.0, {0.5, -0.5}, identity());
  GSHGEOM_Path * wide = rect_path(0, 0, 100, 10);
  gshgeom_limits_default(&limits);
  ASSERT_EQ(gshgeom_hatch_pieces(wide, GSHGEOM_FILL_NONZERO, &hd, &chain, 0, nullptr, &limits, &p), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_hatch_pieces_count(p), 100u);
  gshgeom_hatch_pieces_destroy(p);
  limits.max_dash_steps = 4;
  p = nullptr;
  EXPECT_EQ(gshgeom_hatch_pieces(wide, GSHGEOM_FILL_NONZERO, &hd, &chain, 0, nullptr, &limits, &p), GSHGEOM_ERR_LIMIT);
  gshgeom_path_destroy(wide);
  gshgeom_path_destroy(region);
}

TEST(HatchExpand, EveryEntryChecksTheStyleFirstAndDecidesWidthZeroBeforeTheLines) {
  GSHGEOM_Hatch_Family f;
  GSHGEOM_Hatch h = one_family(f, 0, 0, 0, 0, 2.0, {}, identity());
  GSHGEOM_Path * region = rect_path(0, 0, 10, 10);
  Flat flat({{dp(0, 0), dp(10, 0), dp(10, 10), dp(0, 10)}});
  GSHGEOM_Map chain = identity();
  std::vector<double> d = {1.0, 1.0};
  GSHGEOM_Stroke dashed = style(1.0, GSHGEOM_STROKE_CAP_BUTT);
  dashed.dash_count = 2;
  dashed.dashes = d.data();
  GSHGEOM_Stroke aligned = style(1.0, GSHGEOM_STROKE_CAP_BUTT);
  aligned.alignment = GSHGEOM_STROKE_ALIGN_INSIDE;
  GSHGEOM_Stroke zero = style(0.0, GSHGEOM_STROKE_CAP_BUTT);
  for (int entry = 0; entry < 3; entry++) {
    auto run = [&](const GSHGEOM_Stroke & st, uint64_t max_lines, GSHGEOM_Stroke_Outline ** o) {
      switch (entry) {
        case 0: return gshgeom_hatch_expand(region, GSHGEOM_FILL_NONZERO, &h, &st, &chain, max_lines, nullptr, nullptr, o);
        case 1: return gshgeom_hatch_expand_with_tolerance(region, GSHGEOM_FILL_NONZERO, &h, &st, &chain, GSHGEOM_FLATTEN_TOLERANCE, max_lines, nullptr, nullptr, o);
        default: return gshgeom_hatch_expand_from_flat(flat.p, GSHGEOM_FILL_NONZERO, &h, &st, &chain, max_lines, nullptr, nullptr, o);
      }
    };
    GSHGEOM_Stroke_Outline * o = nullptr;
    // six lines against a limit of one: the style is refused before the lines are made, so it is INVALID and not LIMIT
    EXPECT_EQ(run(dashed, 1, &o), GSHGEOM_ERR_INVALID) << entry;
    EXPECT_EQ(run(aligned, 1, &o), GSHGEOM_ERR_INVALID) << entry;
    EXPECT_EQ(run(style(1.0, GSHGEOM_STROKE_CAP_BUTT), 1, &o), GSHGEOM_ERR_LIMIT) << entry;
    ASSERT_EQ(run(zero, 1, &o), GSHGEOM_OK) << entry << ": width 0 is decided before the lines";
    EXPECT_EQ(gshgeom_stroke_outline_status(o), GSHGEOM_STROKE_EMPTY);
    gshgeom_stroke_outline_destroy(o);
  }
  // the pieces entry takes no dashes or alignment either
  GSHGEOM_Hatch_Pieces * p = nullptr;
  ASSERT_EQ(gshgeom_hatch_pieces(region, GSHGEOM_FILL_NONZERO, &h, &chain, 0, nullptr, nullptr, &p), GSHGEOM_OK);
  GSHGEOM_Stroke_Outline * o = nullptr;
  EXPECT_EQ(gshgeom_hatch_expand_pieces(p, &dashed, &chain, GSHGEOM_FLATTEN_TOLERANCE, nullptr, nullptr, &o), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_hatch_expand_pieces(p, &aligned, &chain, GSHGEOM_FLATTEN_TOLERANCE, nullptr, nullptr, &o), GSHGEOM_ERR_INVALID);
  // the tolerance reaches the stroke: a wide round cap at the finest tolerance has more vertices than at the default
  GSHGEOM_Stroke wide = style(6.0, GSHGEOM_STROKE_CAP_ROUND);
  GSHGEOM_Stroke_Outline *a = nullptr, *b = nullptr, *c = nullptr;
  ASSERT_EQ(gshgeom_hatch_expand_pieces(p, &wide, &chain, GSHGEOM_FLATTEN_TOLERANCE, nullptr, nullptr, &a), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_hatch_expand_pieces(p, &wide, &chain, GSHGEOM_FLATTEN_MIN_TOLERANCE, nullptr, nullptr, &b), GSHGEOM_OK);
  EXPECT_NE(gshgeom_stroke_outline_hash(a), gshgeom_stroke_outline_hash(b));
  ASSERT_EQ(gshgeom_hatch_expand_with_tolerance(region, GSHGEOM_FILL_NONZERO, &h, &wide, &chain, GSHGEOM_FLATTEN_MIN_TOLERANCE, 0, nullptr, nullptr, &c),
            GSHGEOM_OK);
  EXPECT_EQ(gshgeom_stroke_outline_hash(b), gshgeom_stroke_outline_hash(c));
  gshgeom_stroke_outline_destroy(a);
  gshgeom_stroke_outline_destroy(b);
  gshgeom_stroke_outline_destroy(c);
  gshgeom_hatch_pieces_destroy(p);
  gshgeom_path_destroy(region);
}

TEST(HatchExpand, APieceWhoseEndsCoincideInDoublesIsADotDrawnAlongTheMappedLineUnderOmit) {
  // the hatch map shears (x, y) to (x, x + y) and the line's origin is 2^28 in hatch space, so a chord 2.5e-8 wide, wider than the
  // band, has both ends on the same pair of doubles: the piece is a segment whose ends coincide, drawn as a dot along the line's
  // direction in local space, (1, 1), which a square cap turns by 45 degrees there and the chain's inverse shear takes back
  const std::vector<std::vector<DP>> tri = {{dp(0, 0), dp(100, 0), dp(50, 10)}};
  Flat region(tri);
  const double hlin[4] = {1, 0, 1, 1}, hzero[2] = {0, 0};
  GSHGEOM_Map hmap;
  ASSERT_EQ(gshgeom_map_from_matrix(2, hlin, hzero, &hmap), GSHGEOM_OK);
  GSHGEOM_Hatch_Family f;
  GSHGEOM_Hatch h = one_family(f, 0, 268435456.0, 10.0 - 2.5e-9, 0, 1000.0, {}, hmap);
  GSHGEOM_Map chain;
  const double clin[4] = {1, 0, -1, 1}, ct[2] = {-268435456.0, 0};
  ASSERT_EQ(gshgeom_map_from_matrix(2, clin, ct, &chain), GSHGEOM_OK);
  GSHGEOM_Hatch_Pieces * p = nullptr;
  ASSERT_EQ(gshgeom_hatch_pieces_from_flat(region.p, GSHGEOM_FILL_NONZERO, &h, &chain, 0, nullptr, nullptr, &p), GSHGEOM_OK);
  std::vector<int> dots;
  const auto pcs = pieces_of(p, &dots);
  ASSERT_EQ(pcs.size(), 1u);
  EXPECT_EQ(dots[0], 0) << "found as a segment";
  EXPECT_EQ(pcs[0][0], pcs[0][2]) << "whose ends are the same double";
  EXPECT_EQ(pcs[0][1], pcs[0][3]);
  GSHGEOM_Stroke st = style(2.0, GSHGEOM_STROKE_CAP_SQUARE);
  st.degenerate = GSHGEOM_STROKE_DEGENERATE_OMIT;
  GSHGEOM_Stroke_Outline * o = nullptr;
  ASSERT_EQ(gshgeom_hatch_expand_pieces(p, &st, &chain, GSHGEOM_FLATTEN_TOLERANCE, nullptr, nullptr, &o), GSHGEOM_OK);
  const GSHGEOM_Flat_Path * body = gshgeom_stroke_outline_body(o);
  ASSERT_NE(body, nullptr) << "drawn although degenerate is omit";
  EXPECT_NE(winding(body, 50.6, 10.0), 0) << "inside the turned square, which reaches 0.707 along the device x axis";
  EXPECT_EQ(winding(body, 50.8, 10.0), 0) << "past it (the square along the line's own axes would reach 1)";
  gshgeom_stroke_outline_destroy(o);
  gshgeom_hatch_pieces_destroy(p);
}

TEST(HatchExpand, EachDotIsDrawnAlongItsOwnDirectionMappedByTheMapUnderDeviceWidth) {
  GSHGEOM_Path * path = nullptr;
  ASSERT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &path), GSHGEOM_OK);
  for (double x : {10.0, 30.0}) {
    GSHGEOM_Point pt = {x, 0.0, 0.0};
    ASSERT_EQ(gshgeom_path_move_to(path, pt), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_path_line_to(path, pt), GSHGEOM_OK);
  }
  GSHGEOM_Stroke st = style(2.0, GSHGEOM_STROKE_CAP_SQUARE);
  st.degenerate = GSHGEOM_STROKE_DEGENERATE_OMIT;
  GSHGEOM_Map id = identity();
  const double dirs[4] = {1.0, 0.0, 1.0, 1.0};
  GSHGEOM_Stroke_Outline * o = nullptr;
  ASSERT_EQ(gshgeom_stroke_expand_dots(path, &st, &id, GSHGEOM_FLATTEN_TOLERANCE, nullptr, nullptr, dirs, 2, &o), GSHGEOM_OK);
  const GSHGEOM_Flat_Path * body = gshgeom_stroke_outline_body(o);
  ASSERT_NE(body, nullptr);
  EXPECT_EQ(winding(body, 11.3, 0.01), 0) << "the first dot is the axis-aligned square";
  EXPECT_NE(winding(body, 31.3, 0.01), 0) << "the second is turned by its own direction, not the first's";
  gshgeom_stroke_outline_destroy(o);
  // a degenerate subpath beyond the directions is an internal error, not a read past them
  EXPECT_EQ(gshgeom_stroke_expand_dots(path, &st, &id, GSHGEOM_FLATTEN_TOLERANCE, nullptr, nullptr, dirs, 1, &o), GSHGEOM_ERR_INTERNAL);
  // under a device width the direction goes through the map's linear part: (0, 1) under a shear is (1, 1)
  GSHGEOM_Stroke dev = style(2.0, GSHGEOM_STROKE_CAP_SQUARE, GSHGEOM_STROKE_WIDTH_DEVICE);
  dev.degenerate = GSHGEOM_STROKE_DEGENERATE_OMIT;
  const double shear[4] = {1, 1, 0, 1}, zero2[2] = {0, 0};
  GSHGEOM_Map sh;
  ASSERT_EQ(gshgeom_map_from_matrix(2, shear, zero2, &sh), GSHGEOM_OK);
  const double down[2] = {0.0, 1.0};
  GSHGEOM_Path * one = nullptr;
  ASSERT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &one), GSHGEOM_OK);
  GSHGEOM_Point pt = {10.0, 0.0, 0.0};
  ASSERT_EQ(gshgeom_path_move_to(one, pt), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(one, pt), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_stroke_expand_dots(one, &dev, &sh, GSHGEOM_FLATTEN_TOLERANCE, nullptr, nullptr, down, 1, &o), GSHGEOM_OK);
  body = gshgeom_stroke_outline_body(o);
  ASSERT_NE(body, nullptr);
  EXPECT_NE(winding(body, 11.3, 0.01), 0) << "turned a quarter of the way: reaches past the axis-aligned side";
  gshgeom_stroke_outline_destroy(o);
  gshgeom_path_destroy(one);
  gshgeom_path_destroy(path);
}
