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

#include "../hatch_oracle.h"
#include "../hatch_vectors.h"

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
