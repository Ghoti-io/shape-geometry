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
 * A stroke's parameters (SG1 to SG3), the statuses (SG14, SG15), the inputs and
 * their refusals, the parts and their order, the outline's hash and dump
 * (SG16, SG18), the limits, allocation failure, and the exact decisions of SG9.
 * Geometry against the definition is test_stroke_matrix.cpp's.
 */

#include "../test_helpers.h"

#include "../alloc_tracker.h"
#include "../dump_sink.h"
#include "../stroke_matrix.h"

extern "C" {
#include "../../src/exact/exact.h"
#include "../../src/stroke/stroke_internal.h"
}

#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

namespace sm = stroke_matrix;

namespace {

GSHGEOM_Map identity() {
  GSHGEOM_Map m;
  gshgeom_map_identity(2, &m);
  return m;
}

GSHGEOM_Map scale(double sx, double sy) {
  GSHGEOM_Map m;
  gshgeom_map_scale(2, sx, sy, 1.0, &m);
  return m;
}

GSHGEOM_Stroke plain() {
  GSHGEOM_Stroke s;
  gshgeom_stroke_init(&s);
  s.width = 2;
  return s;
}

GSHGEOM_Path * corner() {
  return sm::build({{'M', {0, 0}}, {'L', {10, 0}}, {'Q', {14, 0, 14, 6}}, {'L', {6, 9}}});
}

uint64_t outline_hash(const GSHGEOM_Path * p, const GSHGEOM_Stroke & s, const GSHGEOM_Map & m) {
  GSHGEOM_Stroke_Outline * o = nullptr;
  EXPECT_EQ(gshgeom_stroke_expand(p, &s, &m, nullptr, nullptr, &o), GSHGEOM_OK);
  uint64_t h = gshgeom_stroke_outline_hash(o);
  gshgeom_stroke_outline_destroy(o);
  return h;
}

} // namespace

/* --- SG1: the parameters ------------------------------------------------------- */

TEST(StrokeParams, TheEnumerationsCarryTheScenesPinnedValues) {
  EXPECT_EQ(GSHGEOM_STROKE_WIDTH_LOCAL, 0);
  EXPECT_EQ(GSHGEOM_STROKE_WIDTH_DEVICE, 1);
  EXPECT_EQ(GSHGEOM_STROKE_JOIN_MITER, 0);
  EXPECT_EQ(GSHGEOM_STROKE_JOIN_ROUND, 1);
  EXPECT_EQ(GSHGEOM_STROKE_JOIN_BEVEL, 2);
  EXPECT_EQ(GSHGEOM_STROKE_JOIN_MITER_CLIP, 3);
  EXPECT_EQ(GSHGEOM_STROKE_JOIN_ARCS, 4);
  EXPECT_EQ(GSHGEOM_STROKE_JOIN_COUNT, 5);
  EXPECT_EQ(GSHGEOM_STROKE_CAP_BUTT, 0);
  EXPECT_EQ(GSHGEOM_STROKE_CAP_ROUND, 1);
  EXPECT_EQ(GSHGEOM_STROKE_CAP_SQUARE, 2);
  EXPECT_EQ(GSHGEOM_STROKE_DEGENERATE_DRAW, 0);
  EXPECT_EQ(GSHGEOM_STROKE_DEGENERATE_OMIT, 1);
  EXPECT_EQ(GSHGEOM_STROKE_DASH_PER_SUBPATH, 0);
  EXPECT_EQ(GSHGEOM_STROKE_DASH_PATH, 1);
  EXPECT_EQ(GSHGEOM_STROKE_ALIGN_CENTER, 0);
  EXPECT_EQ(GSHGEOM_STROKE_ALIGN_INSIDE, 1);
  EXPECT_EQ(GSHGEOM_STROKE_ALIGN_OUTSIDE, 2);
  EXPECT_EQ(GSHGEOM_STROKE_DRAWN, 0);
  EXPECT_EQ(GSHGEOM_STROKE_EMPTY, 1);
  EXPECT_EQ(GSHGEOM_STROKE_NOT_INVERTIBLE, 2);
}

TEST(StrokeParams, InitWritesTheDefaults) {
  GSHGEOM_Stroke s;
  std::memset(&s, 0x5a, sizeof s);
  gshgeom_stroke_init(&s);
  EXPECT_EQ(s.width, 1.0);
  EXPECT_EQ(s.width_space, 0u);
  EXPECT_EQ(s.join, 0u);
  EXPECT_EQ(s.miter_limit, 4.0);
  EXPECT_EQ(s.cap, 0u);
  EXPECT_EQ(s.degenerate, 0u);
  EXPECT_EQ(s.dash_count, 0u);
  EXPECT_EQ(s.dashes, nullptr);
  EXPECT_EQ(s.dash_offset, 0.0);
  EXPECT_EQ(s.dash_scope, 0u);
  EXPECT_EQ(s.alignment, 0u);
  EXPECT_EQ(s.alignment_rule, GSHGEOM_FILL_NONZERO);
  EXPECT_EQ(s.width_profile, nullptr);
  for (uint8_t b : s._reserved) {
    EXPECT_EQ(b, 0);
  }
  EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_OK);
  gshgeom_stroke_init(nullptr);
}

TEST(StrokeParams, AZeroFilledStrokeIsNotValid) {
  GSHGEOM_Stroke s;
  std::memset(&s, 0, sizeof s);
  EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_ERR_LIMIT); /* miter limit 0 is out of range */
}

/* --- SG2: validity and bounds, in SK2's order -------------------------------------- */

TEST(StrokeParams, InvalidIsDecidedFirst) {
  GSHGEOM_Stroke s = plain();
  EXPECT_EQ(gshgeom_stroke_check(nullptr), GSHGEOM_ERR_INVALID);
  s.dash_count = 2;
  EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_ERR_INVALID);
  const double nan = std::numeric_limits<double>::quiet_NaN(), inf = std::numeric_limits<double>::infinity();
  for (double bad : {nan, inf, -inf}) {
    s = plain();
    s.width = bad;
    EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_ERR_INVALID);
    s = plain();
    s.miter_limit = bad;
    EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_ERR_INVALID);
    s = plain();
    s.dash_offset = bad;
    EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_ERR_INVALID);
    double d[2] = {1.0, bad};
    s = plain();
    s.dashes = d;
    s.dash_count = 2;
    EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_ERR_INVALID);
  }
  uint32_t * fields[] = {&s.width_space, &s.join, &s.cap, &s.degenerate, &s.dash_scope, &s.alignment};
  uint32_t counts[] = {2, 5, 3, 2, 2, 3};
  for (size_t i = 0; i < 6; i++) {
    s = plain();
    *fields[i] = counts[i];
    EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_ERR_INVALID) << i;
    s = plain();
    *fields[i] = counts[i] - 1;
    EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_OK) << i;
  }
  s = plain();
  s.alignment_rule = GSHGEOM_FILL_RULE_COUNT;
  EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_ERR_INVALID);
  s = plain();
  s.width = -1;
  EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_ERR_INVALID);
  s = plain();
  s._reserved[7] = 1;
  EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_ERR_INVALID);
  /* an invalid field wins over a limit and over the profile */
  s = plain();
  s.width = -1;
  s.miter_limit = 0.5;
  int profile = 0;
  s.width_profile = &profile;
  EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_ERR_INVALID);
}

TEST(StrokeParams, AWidthProfileIsUnsupportedBeforeAnyLimit) {
  GSHGEOM_Stroke s = plain();
  int profile = 0;
  s.width_profile = &profile;
  EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_ERR_UNSUPPORTED);
  s.miter_limit = 0.5;
  EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_ERR_UNSUPPORTED);
}

TEST(StrokeParams, EveryBoundAtAndPastItsEdge) {
  /* the constants equal the scene's (SK2): 2^-64, 2^64, 2^16, 65535 */
  EXPECT_EQ(GSHGEOM_STROKE_VALUE_MIN, 0x1p-64);
  EXPECT_EQ(GSHGEOM_STROKE_VALUE_MAX, 0x1p64);
  EXPECT_EQ(GSHGEOM_STROKE_MITER_LIMIT_MAX, 0x1p16);
  EXPECT_EQ(GSHGEOM_STROKE_DASHES_MAX, (size_t)65535);
  struct Edge {
    double ok, bad;
  };
  for (Edge e : {Edge{0x1p-64, std::nextafter(0x1p-64, 0.0)}, Edge{0x1p64, std::nextafter(0x1p64, 1e300)}}) {
    GSHGEOM_Stroke s = plain();
    s.width = e.ok;
    EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_OK);
    s.width = e.bad;
    EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_ERR_LIMIT);
    for (double sign : {1.0, -1.0}) {
      s = plain();
      s.dash_offset = sign * e.ok;
      EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_OK);
      s.dash_offset = sign * e.bad;
      EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_ERR_LIMIT);
      double d[2] = {1.0, sign * e.ok};
      s = plain();
      s.dashes = d;
      s.dash_count = 2;
      EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_OK);
      d[1] = sign * e.bad;
      EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_ERR_LIMIT);
    }
  }
  GSHGEOM_Stroke s = plain();
  s.width = 0.0;
  EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_OK);
  s = plain();
  s.miter_limit = 1.0;
  EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_OK);
  s.miter_limit = std::nextafter(1.0, 0.0);
  EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_ERR_LIMIT);
  s.miter_limit = 0x1p16;
  EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_OK);
  s.miter_limit = std::nextafter(0x1p16, 1e300);
  EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_ERR_LIMIT);
  /* the dash count after an odd array is repeated */
  std::vector<double> d(65535, 1.0);
  s = plain();
  s.dashes = d.data();
  s.dash_count = 65534; /* even: 65534 */
  EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_OK);
  s.dash_count = 65535; /* odd: 131070 after repeating */
  EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_ERR_LIMIT);
  d[3] = -1.0; /* solid: no dashes, so no count */
  EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_OK);
}

/* --- SG3: the canonical meaning ---------------------------------------------------- */

TEST(StrokeCanonical, AStyleAndItsCanonicalFormGiveTheSameBits) {
  GSHGEOM_Path * p = corner();
  GSHGEOM_Map m = identity();
  GSHGEOM_Stroke base = plain();
  base.join = GSHGEOM_STROKE_JOIN_ROUND;
  uint64_t h = outline_hash(p, base, m);
  GSHGEOM_Stroke s = base;
  s.miter_limit = 10; /* does not apply to a round join */
  EXPECT_EQ(outline_hash(p, s, m), h);
  double zeros[2] = {0, 0}, neg[2] = {1, -1};
  s = base;
  s.dashes = zeros;
  s.dash_count = 2;
  s.dash_offset = 5;
  s.dash_scope = GSHGEOM_STROKE_DASH_PATH;
  EXPECT_EQ(outline_hash(p, s, m), h);
  s.dashes = neg;
  EXPECT_EQ(outline_hash(p, s, m), h);
  s = base;
  s.alignment_rule = GSHGEOM_FILL_EVEN_ODD; /* does not apply when centered */
  EXPECT_EQ(outline_hash(p, s, m), h);
  /* -0 is +0: width -0 is width 0 */
  GSHGEOM_Stroke z = base, nz = base;
  z.width = 0.0;
  nz.width = -0.0;
  EXPECT_EQ(outline_hash(p, z, m), outline_hash(p, nz, m));
  gshgeom_path_destroy(p);
}

TEST(StrokeCanonical, DashesAreDrawnNeverSolid) {
  GSHGEOM_Path * p = corner();
  GSHGEOM_Map m = identity();
  GSHGEOM_Stroke solid = plain();
  uint64_t solid_hash = outline_hash(p, solid, m);
  double d1[2] = {2, 1}, d2[3] = {1, 2, 3}, d3[6] = {1, 2, 3, 1, 2, 3};
  GSHGEOM_Stroke a = plain(), b = plain(), c = plain();
  a.dashes = d1;
  a.dash_count = 2;
  b.dashes = d2;
  b.dash_count = 3;
  c.dashes = d3;
  c.dash_count = 6;
  EXPECT_NE(outline_hash(p, a, m), solid_hash);
  EXPECT_NE(outline_hash(p, b, m), solid_hash);
  EXPECT_NE(outline_hash(p, a, m), outline_hash(p, b, m));
  /* SG3: an odd array is the array repeated once */
  EXPECT_EQ(outline_hash(p, b, m), outline_hash(p, c, m));
  /* fields that do not apply to a solid style: the offset and the scope */
  double zeros[2] = {0, 0}, neg[2] = {1, -1};
  for (double * z : {zeros, neg}) {
    GSHGEOM_Stroke s = plain();
    s.dashes = z;
    s.dash_count = 2;
    s.dash_offset = 3;
    s.dash_scope = GSHGEOM_STROKE_DASH_PATH;
    EXPECT_EQ(outline_hash(p, s, m), solid_hash);
  }
  gshgeom_path_destroy(p);
}

/* --- SG4: inputs -------------------------------------------------------------------- */

TEST(StrokeInputs, RefusalsAndTolerances) {
  GSHGEOM_Path * p = corner();
  GSHGEOM_Map m = identity();
  GSHGEOM_Stroke s = plain();
  GSHGEOM_Stroke_Outline * o = nullptr;
  EXPECT_EQ(gshgeom_stroke_expand(nullptr, &s, &m, nullptr, nullptr, &o), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_stroke_expand(p, nullptr, &m, nullptr, nullptr, &o), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_stroke_expand(p, &s, nullptr, nullptr, nullptr, &o), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_stroke_expand(p, &s, &m, nullptr, nullptr, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_stroke_expand_with_tolerance(p, &s, &m, 1, nullptr, nullptr, &o), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(o, nullptr);
  GSHGEOM_Map m3;
  gshgeom_map_identity(3, &m3);
  EXPECT_EQ(gshgeom_stroke_expand(p, &s, &m3, nullptr, nullptr, &o), GSHGEOM_ERR_INVALID);
  GSHGEOM_Map bad = m;
  bad.m[0][0] = std::numeric_limits<double>::quiet_NaN();
  EXPECT_EQ(gshgeom_stroke_expand(p, &s, &bad, nullptr, nullptr, &o), GSHGEOM_ERR_INVALID);
  GSHGEOM_Path * p3 = nullptr;
  ASSERT_EQ(gshgeom_path_create(nullptr, nullptr, 3, &p3), GSHGEOM_OK);
  GSHGEOM_Point a = {0, 0, 0}, b = {1, 1, 1};
  gshgeom_path_move_to(p3, a);
  gshgeom_path_line_to(p3, b);
  EXPECT_EQ(gshgeom_stroke_expand(p3, &s, &m, nullptr, nullptr, &o), GSHGEOM_ERR_UNSUPPORTED);
  gshgeom_path_destroy(p3);
  /* tolerances 2, 4 and 64 are each expanded; 4 is the default exactly */
  uint64_t h4 = 0;
  for (int32_t tau : {2, 4, 64}) {
    ASSERT_EQ(gshgeom_stroke_expand_with_tolerance(p, &s, &m, tau, nullptr, nullptr, &o), GSHGEOM_OK) << tau;
    if (tau == 4) {
      h4 = gshgeom_stroke_outline_hash(o);
    }
    gshgeom_stroke_outline_destroy(o);
  }
  EXPECT_EQ(outline_hash(p, s, m), h4);
  gshgeom_path_destroy(p);
}

/* --- SG14, SG15: statuses ----------------------------------------------------------- */

TEST(StrokeStatus, ZeroWidthIsEmptyDegenerateDotsIncluded) {
  GSHGEOM_Path * p = sm::build({{'M', {0, 0}}, {'L', {5, 0}}, {'M', {9, 9}}, {'Z', {}}});
  GSHGEOM_Map m = identity();
  GSHGEOM_Stroke s = plain();
  s.width = 0.0;
  s.cap = GSHGEOM_STROKE_CAP_ROUND;
  GSHGEOM_Stroke_Outline * o = nullptr;
  ASSERT_EQ(gshgeom_stroke_expand(p, &s, &m, nullptr, nullptr, &o), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_stroke_outline_status(o), GSHGEOM_STROKE_EMPTY);
  EXPECT_EQ(gshgeom_stroke_outline_body(o), nullptr);
  EXPECT_EQ(gshgeom_stroke_outline_aligned(o), nullptr);
  EXPECT_EQ(gshgeom_stroke_outline_clip(o), nullptr);
  gshgeom_stroke_outline_destroy(o);
  gshgeom_path_destroy(p);
}

TEST(StrokeStatus, ASingularMapIsNotInvertibleAndDrawsNothing) {
  GSHGEOM_Path * p = corner();
  GSHGEOM_Stroke s = plain();
  double rank_one[4] = {1, 2, 2, 4}, t[2] = {0, 0};
  GSHGEOM_Map r1;
  ASSERT_EQ(gshgeom_map_from_matrix(2, rank_one, t, &r1), GSHGEOM_OK);
  /* a product of two invertible maps whose determinant underflows to exactly zero */
  GSHGEOM_Map a = scale(0x1p-600, 0x1p-600), b = scale(0x1p-600, 1.0), ab;
  ASSERT_EQ(gshgeom_map_compose(&a, &b, &ab), GSHGEOM_OK);
  for (const GSHGEOM_Map & m : {scale(0, 1), r1, ab}) {
    int inv = -1;
    ASSERT_EQ(gshgeom_stroke_map_invertible(&m, &inv), GSHGEOM_OK);
    EXPECT_EQ(inv, 0);
    for (uint32_t space : {0u, 1u}) {
      s.width_space = space;
      GSHGEOM_Stroke_Outline * o = nullptr;
      ASSERT_EQ(gshgeom_stroke_expand(p, &s, &m, nullptr, nullptr, &o), GSHGEOM_OK);
      EXPECT_EQ(gshgeom_stroke_outline_status(o), GSHGEOM_STROKE_NOT_INVERTIBLE);
      EXPECT_EQ(gshgeom_stroke_outline_body(o), nullptr);
      gshgeom_stroke_outline_destroy(o);
    }
  }
  /* decided before zero width */
  s.width = 0.0;
  GSHGEOM_Map z = scale(0, 1);
  GSHGEOM_Stroke_Outline * o = nullptr;
  ASSERT_EQ(gshgeom_stroke_expand(p, &s, &z, nullptr, nullptr, &o), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_stroke_outline_status(o), GSHGEOM_STROKE_NOT_INVERTIBLE);
  gshgeom_stroke_outline_destroy(o);
  int inv = -1;
  GSHGEOM_Map id = identity();
  EXPECT_EQ(gshgeom_stroke_map_invertible(&id, &inv), GSHGEOM_OK);
  EXPECT_EQ(inv, 1);
  EXPECT_EQ(gshgeom_stroke_map_invertible(nullptr, &inv), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_stroke_map_invertible(&id, nullptr), GSHGEOM_ERR_INVALID);
  GSHGEOM_Map m3;
  gshgeom_map_identity(3, &m3);
  EXPECT_EQ(gshgeom_stroke_map_invertible(&m3, &inv), GSHGEOM_ERR_INVALID);
  gshgeom_path_destroy(p);
}

/* --- SG12, SG16: parts ---------------------------------------------------------------- */

TEST(StrokeParts, AlignmentGivesBodyAlignedAndClip) {
  GSHGEOM_Path * p = sm::build({{'M', {0, 0}}, {'L', {8, 0}}, {'L', {8, 8}}, {'Z', {}}, {'M', {12, 0}}, {'L', {16, 4}}});
  GSHGEOM_Map m = identity();
  for (uint32_t align : {0u, 1u, 2u}) {
    for (int rule : {0, 1}) {
      GSHGEOM_Stroke s = plain();
      s.alignment = align;
      s.alignment_rule = (GSHGEOM_Fill_Rule)rule;
      GSHGEOM_Stroke_Outline * o = nullptr;
      ASSERT_EQ(gshgeom_stroke_expand(p, &s, &m, nullptr, nullptr, &o), GSHGEOM_OK);
      EXPECT_NE(gshgeom_stroke_outline_body(o), nullptr);
      EXPECT_EQ(gshgeom_stroke_outline_aligned(o) != nullptr, align != 0);
      EXPECT_EQ(gshgeom_stroke_outline_clip(o) != nullptr, align != 0);
      EXPECT_EQ(gshgeom_stroke_outline_clip_complement(o), align == 2 ? 1 : 0);
      EXPECT_EQ(gshgeom_stroke_outline_clip_rule(o), align == 0 ? GSHGEOM_FILL_NONZERO : (GSHGEOM_Fill_Rule)rule);
      if (align) {
        /* the clip is the closed subpath with extent through the device path and the flattener at the tolerance */
        GSHGEOM_Path * f = sm::build({{'M', {0, 0}}, {'L', {8, 0}}, {'L', {8, 8}}, {'Z', {}}});
        GSHGEOM_Flat_Path * want = nullptr;
        ASSERT_EQ(gshgeom_path_flatten(f, &m, nullptr, nullptr, &want), GSHGEOM_OK);
        EXPECT_EQ(gshgeom_flat_path_hash(gshgeom_stroke_outline_clip(o)), gshgeom_flat_path_hash(want));
        gshgeom_flat_path_destroy(want);
        gshgeom_path_destroy(f);
      }
      gshgeom_stroke_outline_destroy(o);
    }
  }
  gshgeom_path_destroy(p);
}

TEST(StrokeParts, NoClosedSubpathMeansNoAlignedPart) {
  GSHGEOM_Path * p = corner();
  GSHGEOM_Map m = identity();
  GSHGEOM_Stroke s = plain();
  s.alignment = GSHGEOM_STROKE_ALIGN_INSIDE;
  GSHGEOM_Stroke_Outline * o = nullptr;
  ASSERT_EQ(gshgeom_stroke_expand(p, &s, &m, nullptr, nullptr, &o), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_stroke_outline_aligned(o), nullptr);
  EXPECT_EQ(gshgeom_stroke_outline_clip(o), nullptr);
  gshgeom_stroke_outline_destroy(o);
  gshgeom_path_destroy(p);
}

TEST(StrokeParts, EveryRingIsClosedAndEveryVertexIsTwentySixSix) {
  GSHGEOM_Path * p = corner();
  GSHGEOM_Map m = scale(1.37, 2.11);
  GSHGEOM_Stroke s = plain();
  s.join = GSHGEOM_STROKE_JOIN_ROUND;
  GSHGEOM_Stroke_Outline * o = nullptr;
  ASSERT_EQ(gshgeom_stroke_expand(p, &s, &m, nullptr, nullptr, &o), GSHGEOM_OK);
  const GSHGEOM_Flat_Path * b = gshgeom_stroke_outline_body(o);
  ASSERT_GT(gshgeom_flat_path_subpath_count(b), 3u);
  for (size_t i = 0; i < gshgeom_flat_path_subpath_count(b); i++) {
    GSHGEOM_Flat_Subpath sp;
    ASSERT_EQ(gshgeom_flat_path_subpath(b, i, &sp), GSHGEOM_OK);
    EXPECT_EQ(sp.closed, 1);
    EXPECT_GE(sp.count, 3u);
  }
  gshgeom_stroke_outline_destroy(o);
  gshgeom_path_destroy(p);
}

/* --- SG18: determinism, revision, hash and dump --------------------------------------- */

TEST(StrokeHash, TheRevisionIsOne) {
  EXPECT_EQ(GSHGEOM_STROKE_REVISION, 1);
  EXPECT_EQ(gshgeom_stroke_revision(), 1);
}

TEST(StrokeHash, TheSpellingIsRestatedFieldByField) {
  GSHGEOM_Path * p = sm::build({{'M', {0, 0}}, {'L', {8, 0}}, {'L', {8, 8}}, {'Z', {}}, {'M', {12, 0}}, {'L', {16, 4}}});
  GSHGEOM_Map m = identity();
  for (uint32_t align : {0u, 2u}) {
    GSHGEOM_Stroke s = plain();
    s.alignment = align;
    s.alignment_rule = GSHGEOM_FILL_EVEN_ODD;
    GSHGEOM_Stroke_Outline * o = nullptr;
    ASSERT_EQ(gshgeom_stroke_expand(p, &s, &m, nullptr, nullptr, &o), GSHGEOM_OK);
    GSHGEOM_Hash h;
    gshgeom_hash_init(&h);
    gshgeom_hash_u32(&h, 1);
    gshgeom_hash_u32(&h, 0x4b525453u);
    gshgeom_hash_u8(&h, 0);
    gshgeom_hash_u8(&h, (uint8_t)align);
    gshgeom_hash_u8(&h, align ? 1 : 0); /* the canonical rule: even-odd off center, nonzero centered */
    gshgeom_hash_u8(&h, align == 2 ? 1 : 0);
    gshgeom_hash_u64(&h, gshgeom_flat_path_hash(gshgeom_stroke_outline_body(o)));
    gshgeom_hash_u64(&h, gshgeom_flat_path_hash(gshgeom_stroke_outline_aligned(o)));
    gshgeom_hash_u64(&h, gshgeom_flat_path_hash(gshgeom_stroke_outline_clip(o)));
    EXPECT_EQ(gshgeom_stroke_outline_hash(o), gshgeom_hash_value(&h));
    gshgeom_stroke_outline_destroy(o);
  }
  EXPECT_EQ(gshgeom_stroke_outline_hash(nullptr), GSHGEOM_FNV_OFFSET_BASIS);
  /* an absent part hashes as the offset basis (the empty stream) */
  GSHGEOM_Stroke s = plain();
  s.width = 0.0;
  GSHGEOM_Stroke_Outline * o = nullptr;
  ASSERT_EQ(gshgeom_stroke_expand(p, &s, &m, nullptr, nullptr, &o), GSHGEOM_OK);
  GSHGEOM_Hash h;
  gshgeom_hash_init(&h);
  gshgeom_hash_u32(&h, 1);
  gshgeom_hash_u32(&h, 0x4b525453u);
  gshgeom_hash_u8(&h, 1);
  gshgeom_hash_u8(&h, 0);
  gshgeom_hash_u8(&h, 0);
  gshgeom_hash_u8(&h, 0);
  for (int k = 0; k < 3; k++) {
    gshgeom_hash_u64(&h, GSHGEOM_FNV_OFFSET_BASIS);
  }
  EXPECT_EQ(gshgeom_stroke_outline_hash(o), gshgeom_hash_value(&h));
  gshgeom_stroke_outline_destroy(o);
  gshgeom_path_destroy(p);
}

TEST(StrokeHash, TheSameInputGivesTheSameBits) {
  GSHGEOM_Path * p = corner();
  GSHGEOM_Map m = scale(1.5, 0.75);
  GSHGEOM_Stroke s = plain();
  s.join = GSHGEOM_STROKE_JOIN_ARCS;
  s.cap = GSHGEOM_STROKE_CAP_ROUND;
  EXPECT_EQ(outline_hash(p, s, m), outline_hash(p, s, m));
  gshgeom_path_destroy(p);
}

TEST(StrokeDump, TheTextNamesEveryPart) {
  GSHGEOM_Path * p = sm::build({{'M', {0, 0}}, {'L', {8, 0}}, {'L', {8, 8}}, {'Z', {}}});
  GSHGEOM_Map m = identity();
  GSHGEOM_Stroke s = plain();
  s.alignment = GSHGEOM_STROKE_ALIGN_INSIDE;
  GSHGEOM_Stroke_Outline * o = nullptr;
  ASSERT_EQ(gshgeom_stroke_expand(p, &s, &m, nullptr, nullptr, &o), GSHGEOM_OK);
  char * buf = nullptr;
  size_t len = 0;
  FILE * f = open_memstream(&buf, &len);
  ASSERT_NE(f, nullptr);
  ASSERT_EQ(gshgeom_stroke_outline_dump(o, f), GSHGEOM_OK);
  std::fclose(f);
  std::string text(buf, len);
  std::free(buf);
  char head[160];
  std::snprintf(head, sizeof head, "gshgeom.stroke-outline 1 status=0 align=1 rule=0 complement=0 hash=%016llx\n",
      (unsigned long long)gshgeom_stroke_outline_hash(o));
  EXPECT_EQ(text.compare(0, std::strlen(head), head), 0) << text.substr(0, 200);
  EXPECT_NE(text.find("part body\ngshgeom.flat-path 1"), std::string::npos);
  EXPECT_NE(text.find("part aligned\ngshgeom.flat-path 1"), std::string::npos);
  EXPECT_NE(text.find("part clip\ngshgeom.flat-path 1"), std::string::npos);
  EXPECT_EQ(gshgeom_stroke_outline_dump(nullptr, stdout), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_stroke_outline_dump(o, nullptr), GSHGEOM_ERR_INVALID);
  gshgeom_stroke_outline_destroy(o);
  s.width = 0.0;
  ASSERT_EQ(gshgeom_stroke_expand(p, &s, &m, nullptr, nullptr, &o), GSHGEOM_OK);
  f = open_memstream(&buf, &len);
  ASSERT_EQ(gshgeom_stroke_outline_dump(o, f), GSHGEOM_OK);
  std::fclose(f);
  text.assign(buf, len);
  std::free(buf);
  EXPECT_NE(text.find("part body -\npart aligned -\npart clip -\n"), std::string::npos) << text;
  gshgeom_stroke_outline_destroy(o);
  gshgeom_path_destroy(p);
}

TEST(StrokeDump, EveryWriteIsChecked) {
  GSHGEOM_Path * p = sm::build({{'M', {0, 0}}, {'L', {8, 0}}, {'L', {8, 8}}, {'Z', {}}});
  GSHGEOM_Map m = identity();
  GSHGEOM_Stroke s = plain();
  s.alignment = GSHGEOM_STROKE_ALIGN_OUTSIDE;
  GSHGEOM_Stroke_Outline * o = nullptr;
  ASSERT_EQ(gshgeom_stroke_expand(p, &s, &m, nullptr, nullptr, &o), GSHGEOM_OK);
  GSHGEOM_EVERY_WRITE_CHECKED(gshgeom_stroke_outline_dump, o, "outline dump");
  gshgeom_stroke_outline_destroy(o);
  /* and a dump whose last lines are the outline's own ("part aligned -", "part clip -"), written by nobody else:
   * the final flush is the only check that they arrived */
  s.alignment = GSHGEOM_STROKE_ALIGN_CENTER;
  ASSERT_EQ(gshgeom_stroke_expand(p, &s, &m, nullptr, nullptr, &o), GSHGEOM_OK);
  GSHGEOM_EVERY_WRITE_CHECKED(gshgeom_stroke_outline_dump, o, "outline dump, absent parts last");
  gshgeom_stroke_outline_destroy(o);
  gshgeom_path_destroy(p);
}

/* --- limits and allocation ------------------------------------------------------------ */

TEST(StrokeLimits, AnOutlineOutsideTheRangeOrPastMaxPointsIsRefused) {
  GSHGEOM_Path * p = corner();
  GSHGEOM_Stroke s = plain();
  GSHGEOM_Stroke_Outline * o = nullptr;
  GSHGEOM_Map big = scale(0x1p22, 0x1p22); /* 14 * 2^22 px: past 2^25 */
  EXPECT_EQ(gshgeom_stroke_expand(p, &s, &big, nullptr, nullptr, &o), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(o, nullptr);
  GSHGEOM_Map m = identity();
  GSHGEOM_Limits lim;
  gshgeom_limits_default(&lim);
  lim.max_points = 20;
  EXPECT_EQ(gshgeom_stroke_expand(p, &s, &m, nullptr, &lim, &o), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(o, nullptr);
  /* a magnitude the arithmetic bound is not proved for: a huge translation cancelled by the model coordinates */
  GSHGEOM_Path * far = sm::build({{'M', {0x1p40, 0}}, {'L', {0x1p40 + 4, 0}}});
  double mm[4] = {1, 0, 0, 1}, t[2] = {-0x1p40, 0};
  GSHGEOM_Map shift;
  ASSERT_EQ(gshgeom_map_from_matrix(2, mm, t, &shift), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_stroke_expand(far, &s, &shift, nullptr, nullptr, &o), GSHGEOM_ERR_LIMIT);
  gshgeom_path_destroy(far);
  gshgeom_path_destroy(p);
}

TEST(StrokeAllocation, EveryFailureIsReportedAndNothingLeaks) {
  GSHGEOM_Path * p = sm::build({{'M', {0, 0}}, {'L', {8, 0}}, {'Q', {12, 0, 12, 4}}, {'L', {8, 8}}, {'Z', {}},
      {'M', {20, 20}}, {'Z', {}}});
  GSHGEOM_Map m = identity();
  GSHGEOM_Stroke s = plain();
  s.join = GSHGEOM_STROKE_JOIN_ROUND;
  s.cap = GSHGEOM_STROKE_CAP_ROUND;
  s.alignment = GSHGEOM_STROKE_ALIGN_INSIDE;
  tracker::Tracker count;
  {
    GSHGEOM_Allocator a = tracker::allocator_for(&count);
    GSHGEOM_Stroke_Outline * o = nullptr;
    ASSERT_EQ(gshgeom_stroke_expand(p, &s, &m, &a, nullptr, &o), GSHGEOM_OK);
    gshgeom_stroke_outline_destroy(o);
    EXPECT_TRUE(count.live.empty());
  }
  ASSERT_GT(count.calls, 5);
  for (long k = 0; k < count.calls; k++) {
    for (int persistent = 0; persistent < 2; persistent++) {
      tracker::Tracker t;
      t.fail_from = k;
      t.persistent = persistent != 0;
      GSHGEOM_Allocator a = tracker::allocator_for(&t);
      GSHGEOM_Stroke_Outline * o = nullptr;
      EXPECT_EQ(gshgeom_stroke_expand(p, &s, &m, &a, nullptr, &o), GSHGEOM_ERR_OOM) << k;
      EXPECT_EQ(o, nullptr);
      EXPECT_TRUE(t.live.empty()) << "k=" << k << " leaked";
    }
  }
  gshgeom_path_destroy(p);
}

/* --- SG9, exactly --------------------------------------------------------------------- */

TEST(StrokeMiterLimit, TheRatioFiveQuartersIsDecidedExactly) {
  /* (0,0) (4,3) (8,0): 1 / sin(theta/2) is exactly 5/4.  At L = 1.25 the miter is drawn, at the double below it is
   * not: the tip region (between the bevel chord and the tip) is covered or not. */
  GSHGEOM_Path * p = sm::build({{'M', {0, 0}}, {'L', {4, 3}}, {'L', {8, 0}}});
  GSHGEOM_Map m = identity();
  for (int below = 0; below < 2; below++) {
    GSHGEOM_Stroke s = plain();
    s.width = 8;
    s.miter_limit = below ? std::nextafter(1.25, 0.0) : 1.25;
    GSHGEOM_Stroke_Outline * o = nullptr;
    ASSERT_EQ(gshgeom_stroke_expand(p, &s, &m, nullptr, nullptr, &o), GSHGEOM_OK);
    stroke_oracle::Winding w;
    w.build(stroke_oracle::polygons_of(gshgeom_stroke_outline_body(o)));
    EXPECT_EQ(w.at({4.0, 7.6}) != 0, below == 0) << below;
    gshgeom_stroke_outline_destroy(o);
  }
  gshgeom_path_destroy(p);
}

TEST(StrokeMiterLimit, TheLimitIsDecidedInModelSpaceUnderLocalWidth) {
  /* Under scale(2, 1) the same corner is blunter in device space; local width decides in model space (SG10), so at
   * L = 1.25 the miter is drawn, and under device width (decided on the mapped points) it is decided there. */
  GSHGEOM_Path * p = sm::build({{'M', {0, 0}}, {'L', {4, 3}}, {'L', {8, 0}}});
  GSHGEOM_Map m = scale(1.0, 3.0); /* in device the corner is sharper: (4, 9), (4, -9) */
  GSHGEOM_Stroke s = plain();
  s.width = 8;
  s.miter_limit = 1.25;
  GSHGEOM_Stroke_Outline * o = nullptr;
  ASSERT_EQ(gshgeom_stroke_expand(p, &s, &m, nullptr, nullptr, &o), GSHGEOM_OK);
  stroke_oracle::Winding w;
  w.build(stroke_oracle::polygons_of(gshgeom_stroke_outline_body(o)));
  /* the model tip region (4, 7.6) maps to (4, 22.8): inside the miter */
  EXPECT_NE(w.at({4.0, 22.8}), 0);
  gshgeom_stroke_outline_destroy(o);
  s.width_space = GSHGEOM_STROKE_WIDTH_DEVICE;
  ASSERT_EQ(gshgeom_stroke_expand(p, &s, &m, nullptr, nullptr, &o), GSHGEOM_OK);
  stroke_oracle::Winding wd;
  wd.build(stroke_oracle::polygons_of(gshgeom_stroke_outline_body(o)));
  /* in device space the ratio is 1 / sin(atan(4/9)) > 1.25: a bevel, nothing at the tip's place */
  EXPECT_EQ(wd.at({4.0, 9.0 + 4.5}), 0);
  gshgeom_stroke_outline_destroy(o);
  gshgeom_path_destroy(p);
}

TEST(StrokeSeams, EveryPointBetweenPiecesIsCovered) {
  /* The pieces overlap, so no crack opens where they meet: scan finely across every vertex's normal segments and
   * every split of a curve, points deep inside the region (a crack of any width would show as a zero winding). */
  GSHGEOM_Path * p = sm::build({{'M', {0, 0}}, {'L', {10, 0}}, {'L', {14, 5}}, {'Q', {20, 12, 26, 4}}, {'L', {30, 9}}});
  for (int join : {GSHGEOM_STROKE_JOIN_BEVEL, GSHGEOM_STROKE_JOIN_MITER, GSHGEOM_STROKE_JOIN_ROUND}) {
    GSHGEOM_Stroke s = plain();
    s.width = 4;
    s.join = (uint32_t)join;
    GSHGEOM_Map m = scale(1.37, 1.37);
    GSHGEOM_Stroke_Outline * o = nullptr;
    ASSERT_EQ(gshgeom_stroke_expand(p, &s, &m, nullptr, nullptr, &o), GSHGEOM_OK);
    stroke_oracle::Winding w;
    w.build(stroke_oracle::polygons_of(gshgeom_stroke_outline_body(o)));
    stroke_oracle::Centerline cl = stroke_oracle::read_centerline(p);
    long bad = 0;
    for (const auto & seg : cl.subs[0].segs) {
      for (int k = 0; k <= 16; k++) {
        if ((&seg == &cl.subs[0].segs.front() && k == 0) || (&seg == &cl.subs[0].segs.back() && k == 16)) {
          continue; /* the butt ends are the region's boundary */
        }
        stroke_oracle::Poly poly = stroke_oracle::make_poly(seg);
        stroke_oracle::Real t = (stroke_oracle::Real)k / 16;
        stroke_oracle::V2 c = stroke_oracle::poly_point(poly, t);
        stroke_oracle::V2 n = stroke_oracle::left_of(stroke_oracle::unit(stroke_oracle::poly_dir(poly, t)));
        for (int j = -60; j <= 60; j++) {
          for (int i = -8; i <= 8; i++) {
            /* across the normal at 1/128 unit steps around it, within 0.9 h of the centerline */
            stroke_oracle::V2 q = c + ((stroke_oracle::Real)j / 60 * 0.9L * 2) * n;
            stroke_oracle::V2 d = {q.x * 1.37L + i / 8192.0L, q.y * 1.37L};
            if (w.at(d) == 0) {
              bad++;
            }
          }
        }
      }
    }
    EXPECT_EQ(bad, 0) << join;
    gshgeom_stroke_outline_destroy(o);
  }
  gshgeom_path_destroy(p);
}

TEST(StrokeLimits, AForcedPieceIsRefusedNeverDrawnCoarser) {
  /* with the depth cap at 2, a curve that needs more splits is ERR_LIMIT and returns nothing */
  GSHGEOM_Path * p = sm::build({{'M', {0, 0}}, {'C', {40, 60, 80, -60, 120, 0}}});
  GSHGEOM_Map m = identity();
  GSHGEOM_Stroke s = plain();
  s.width = 12;
  GSHGEOM_Stroke_Outline * o = nullptr;
  EXPECT_EQ(gshgeom_stroke_expand_capped(p, &s, &m, 4, nullptr, nullptr, 2, &o), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(o, nullptr);
  ASSERT_EQ(gshgeom_stroke_expand_capped(p, &s, &m, 4, nullptr, nullptr, GSHGEOM_STROKE_MAX_DEPTH, &o), GSHGEOM_OK);
  gshgeom_stroke_outline_destroy(o);
  gshgeom_path_destroy(p);
}

namespace {

struct Edges {
  std::vector<std::pair<GSHGEOM_SK_Vec, GSHGEOM_SK_Vec>> list;
};

GSHGEOM_Result collect_edge(void * ctx, GSHGEOM_SK_Vec from, GSHGEOM_SK_Vec to) {
  static_cast<Edges *>(ctx)->list.push_back({from, to});
  return GSHGEOM_OK;
}

GSHGEOM_SK_Vec at(double v) {
  return GSHGEOM_SK_Vec{0.0, v}; /* on the clip line x = 0 (its normal is +x, its direction +y): a position is its y */
}

} // namespace

TEST(StrokeNetChords, TwoChordsOverTheSameStretchInOppositeDirectionsCancel) {
  const GSHGEOM_SK_Vec y = {0.0, 0.0}, t = {1.0, 0.0};
  {
    /* the same stretch both ways: nothing remains */
    GSHGEOM_SK_Vec c[2][2] = {{at(1), at(3)}, {at(3), at(1)}};
    Edges e;
    EXPECT_EQ(gshgeom_sk_net_chords(c, 2, y, t, collect_edge, &e), GSHGEOM_OK);
    EXPECT_TRUE(e.list.empty());
  }
  {
    /* 1 -> 3 and 4 -> 2: the stretch 2 to 3 cancels; 1 -> 2 remains forward and 4 -> 3 remains backward */
    GSHGEOM_SK_Vec c[2][2] = {{at(1), at(3)}, {at(4), at(2)}};
    Edges e;
    EXPECT_EQ(gshgeom_sk_net_chords(c, 2, y, t, collect_edge, &e), GSHGEOM_OK);
    ASSERT_EQ(e.list.size(), (size_t)2);
    EXPECT_EQ(e.list[0].first.y, 1.0);
    EXPECT_EQ(e.list[0].second.y, 2.0);
    EXPECT_EQ(e.list[1].first.y, 4.0);
    EXPECT_EQ(e.list[1].second.y, 3.0);
  }
  {
    /* the same direction twice: the stretch counts twice, as winding does */
    GSHGEOM_SK_Vec c[2][2] = {{at(1), at(2)}, {at(1), at(2)}};
    Edges e;
    EXPECT_EQ(gshgeom_sk_net_chords(c, 2, y, t, collect_edge, &e), GSHGEOM_OK);
    EXPECT_EQ(e.list.size(), (size_t)2);
  }
  {
    /* one chord is itself, and a chord split by another's end stays whole as pieces that join */
    GSHGEOM_SK_Vec c[2][2] = {{at(0), at(5)}, {at(2), at(3)}};
    Edges e;
    EXPECT_EQ(gshgeom_sk_net_chords(c, 2, y, t, collect_edge, &e), GSHGEOM_OK);
    ASSERT_EQ(e.list.size(), (size_t)4); /* 0 -> 2, 2 -> 3 twice, 3 -> 5 */
    EXPECT_EQ(e.list[0].first.y, 0.0);
    EXPECT_EQ(e.list[3].second.y, 5.0);
  }
  EXPECT_EQ(gshgeom_sk_net_chords(nullptr, 0, y, t, collect_edge, nullptr), GSHGEOM_OK);
}

TEST(StrokeNorm, NeitherAHugeNorATinyVectorOverflowsOrVanishes) {
  EXPECT_NEAR(gshgeom_sk_norm(GSHGEOM_SK_Vec{3e-300, 4e-300}), 5e-300, 5e-300 * 1e-15);
  EXPECT_NEAR(gshgeom_sk_norm(GSHGEOM_SK_Vec{3e200, -4e200}), 5e200, 5e200 * 1e-15);
  EXPECT_EQ(gshgeom_sk_norm(GSHGEOM_SK_Vec{0.0, 0.0}), 0.0);
  GSHGEOM_SK_Vec u = gshgeom_sk_unit(GSHGEOM_SK_Vec{1e-300, 0.0});
  EXPECT_EQ(u.x, 1.0);
  EXPECT_EQ(u.y, 0.0);
}

TEST(StrokeCanonical, DegenerateHasNoEffectUnderAButtCap) {
  GSHGEOM_Stroke draw = plain(), omit = plain();
  draw.cap = omit.cap = GSHGEOM_STROKE_CAP_BUTT;
  draw.degenerate = GSHGEOM_STROKE_DEGENERATE_DRAW;
  omit.degenerate = GSHGEOM_STROKE_DEGENERATE_OMIT;
  GSHGEOM_Path * p = sm::build({{'M', {0, 0}}, {'L', {10, 0}}, {'M', {20, 20}}, {'Z', {}}});
  GSHGEOM_Map m = identity();
  EXPECT_EQ(outline_hash(p, draw, m), outline_hash(p, omit, m));
  GSHGEOM_Stroke c;
  ASSERT_EQ(gshgeom_stroke_canonical(&omit, &c), GSHGEOM_OK);
  EXPECT_EQ(c.degenerate, (uint32_t)GSHGEOM_STROKE_DEGENERATE_DRAW);
  omit.cap = GSHGEOM_STROKE_CAP_ROUND; /* under any other cap it still matters */
  ASSERT_EQ(gshgeom_stroke_canonical(&omit, &c), GSHGEOM_OK);
  EXPECT_EQ(c.degenerate, (uint32_t)GSHGEOM_STROKE_DEGENERATE_OMIT);
  gshgeom_path_destroy(p);
}

TEST(StrokeCanonical, AMinusZeroWidthIsAZeroWidth) {
  GSHGEOM_Stroke a = plain(), b = plain();
  a.width = 0.0;
  b.width = -0.0;
  GSHGEOM_Path * p = corner();
  GSHGEOM_Map m = identity();
  EXPECT_EQ(outline_hash(p, a, m), outline_hash(p, b, m));
  GSHGEOM_Stroke ca, cb;
  ASSERT_EQ(gshgeom_stroke_canonical(&a, &ca), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_stroke_canonical(&b, &cb), GSHGEOM_OK);
  EXPECT_EQ(std::signbit(cb.width), false);
  EXPECT_EQ(memcmp(&ca, &cb, sizeof ca), 0);
  gshgeom_path_destroy(p);
}

TEST(StrokeMiterLimit, AFilterWiderThanTheRoundingWouldDecideATieByNoise) {
  /* a = (x, 0), b = (7x, 24x): |a||b| = 25 x^2 and a . b = 7 x^2, so the miter is exactly at the limit 5/4. For an x
   * whose square is not representable the floating evaluation of the decision is off by rounding noise, far above
   * the width of a filter 2^-80 of the terms and well inside the real one: only the exact fallback gets every one
   * right. */
  const double xs[] = {1.0 + 0x1p-40, 1.0 - 0x1p-41, 1.5 + 0x1p-45, 1.0 + 0x1p-47, 0.75 + 0x1p-44, 1.0 + 0x1p-35,
      1.25 - 0x1p-43, 1.0 + 0x1p-38};
  for (double x : xs) {
    GSHGEOM_SK_Vec v = {0.0, 0.0}, q = {-x, 0.0}, s = {7.0 * x, 24.0 * x};
    EXPECT_EQ(gshgeom_sk_miter_within(q, v, s, 1.25), 1) << x;
    EXPECT_EQ(gshgeom_sk_miter_within(q, v, s, std::nextafter(1.25, 2.0)), 1) << x;
    EXPECT_EQ(gshgeom_sk_miter_within(q, v, s, std::nextafter(1.25, 0.0)), 0) << x;
  }
}

TEST(StrokeLimits, AMiterTipThatCannotBePlacedWithinTheToleranceIsRefused) {
  /* a turn of 4e-5 radians short of a reversal, width 2000, limit 65536: the tip is 5e7 pixels out, and the error of
   * the arithmetic that places it (2^-49 of the square of the ratio times the length) is hundreds of pixels */
  GSHGEOM_Path * p = sm::build({{'M', {0, 0}}, {'L', {1000, 0}}, {'L', {0, 0.04}}});
  GSHGEOM_Map m = identity();
  GSHGEOM_Stroke s = plain();
  s.width = 2000;
  s.miter_limit = 65536;
  GSHGEOM_Stroke_Outline * o = nullptr;
  EXPECT_EQ(gshgeom_stroke_expand(p, &s, &m, nullptr, nullptr, &o), GSHGEOM_ERR_LIMIT);
  /* and a turn that is merely sharp is drawn */
  GSHGEOM_Path * q = sm::build({{'M', {0, 0}}, {'L', {1000, 0}}, {'L', {200, 100}}});
  s.miter_limit = 10;
  EXPECT_EQ(gshgeom_stroke_expand(q, &s, &m, nullptr, nullptr, &o), GSHGEOM_OK);
  gshgeom_stroke_outline_destroy(o);
  gshgeom_path_destroy(p);
  gshgeom_path_destroy(q);
}

TEST(StrokeLimits, TheSmallestDifferenceBetweenControlPointsIsAStatedCapAtAndPastItsEdge) {
  /* a nonzero difference of neighbouring control points below 2^-200 is ERR_LIMIT: the build-space arithmetic
   * (cubes, fourth powers of a length) is only proved to stay in range above it */
  struct Row {
    double extent;
    double scale;
    GSHGEOM_Result want;
  };
  const Row rows[] = {{0x1p-200, 1.0, GSHGEOM_OK}, {0x1p-201, 1.0, GSHGEOM_ERR_LIMIT}, {1e-70, 1.0, GSHGEOM_ERR_LIMIT}};
  for (const Row & row : rows) {
    GSHGEOM_Path * p = sm::build({{'M', {0, 0}}, {'L', {row.extent, 0}}, {'L', {row.extent, row.extent}}});
    GSHGEOM_Map m = scale(row.scale, row.scale);
    GSHGEOM_Stroke s = plain();
    GSHGEOM_Stroke_Outline * o = nullptr;
    GSHGEOM_Result r = gshgeom_stroke_expand(p, &s, &m, nullptr, nullptr, &o);
    EXPECT_EQ(r, row.want) << "extent " << row.extent;
    if (r == GSHGEOM_OK) {
      gshgeom_stroke_outline_destroy(o);
    }
    gshgeom_path_destroy(p);
  }
}

TEST(StrokeDirections, ANearReversalWithATinyComponentIsDrawnOrRefusedNeverBroken) {
  /* the second segment's direction is (-10, 1e-300): its unit vector sums with the first's normal to a vector whose
   * squares underflow, which an unscaled norm turns into 0 / 0 */
  GSHGEOM_Path * p = sm::build({{'M', {0, 0}}, {'L', {10, 0}}, {'L', {0, 1e-300}}});
  GSHGEOM_Map m = identity();
  for (uint32_t join = 0; join < (uint32_t)GSHGEOM_STROKE_JOIN_COUNT; join++) {
    GSHGEOM_Stroke s = plain();
    s.join = join;
    GSHGEOM_Stroke_Outline * o = nullptr;
    GSHGEOM_Result r = gshgeom_stroke_expand(p, &s, &m, nullptr, nullptr, &o);
    EXPECT_TRUE(r == GSHGEOM_OK || r == GSHGEOM_ERR_LIMIT) << "join " << join << " gave " << r;
    if (r == GSHGEOM_OK) {
      EXPECT_EQ(gshgeom_stroke_outline_status(o), GSHGEOM_STROKE_DRAWN);
      gshgeom_stroke_outline_destroy(o);
    }
  }
  gshgeom_path_destroy(p);
}

TEST(StrokeExactSum, IdentitiesCancelExactlyAndATinyTermDecides) {
  /* (x + y)^2 - x^2 - 2xy - y^2 = 0 for doubles where x + y is exact; a product of two subnormals added or taken
   * away decides the sign; magnitudes from 2^-1074 to 2^1000 and twelve factors */
  const double xs[] = {1.0, 0x1p-1000, 0x1p1000, 3.0, -7.5, 0x1.fffffffffffffp+500, 0x1p-1074};
  for (double x : xs) {
    for (double y : xs) {
      double s = x + y;
      double bb = s - x;
      if ((x - (s - bb)) + (y - bb) != 0.0) {
        continue; /* x + y not exact: TwoSum's error term is not zero */
      }
      for (int tail = -1; tail <= 1; tail++) {
        GSHGEOM_ExactSum e;
        gshgeom_exact_sum_init(&e);
        double f[2];
        f[0] = s;
        f[1] = s;
        gshgeom_exact_sum_add(&e, 0, 2, f);
        f[0] = x;
        f[1] = x;
        gshgeom_exact_sum_add(&e, 1, 2, f);
        double g[3] = {2.0, x, y};
        gshgeom_exact_sum_add(&e, 1, 3, g);
        f[0] = y;
        f[1] = y;
        gshgeom_exact_sum_add(&e, 1, 2, f);
        if (tail) {
          double t2[2] = {0x1p-1074, 0x1p-1074};
          gshgeom_exact_sum_add(&e, tail < 0, 2, t2);
        }
        EXPECT_EQ(gshgeom_exact_sum_sign(&e), tail) << x << " " << y;
      }
    }
  }
  /* twelve factors: 2^1000 ^ 12 against itself */
  GSHGEOM_ExactSum e;
  gshgeom_exact_sum_init(&e);
  double big[12], tiny[12];
  for (int i = 0; i < 12; i++) {
    big[i] = 0x1p971;
    tiny[i] = 0x1p-1074;
  }
  gshgeom_exact_sum_add(&e, 0, 12, big);
  gshgeom_exact_sum_add(&e, 1, 12, big);
  gshgeom_exact_sum_add(&e, 0, 12, tiny);
  EXPECT_EQ(gshgeom_exact_sum_sign(&e), 1);
  gshgeom_exact_sum_add(&e, 1, 12, tiny);
  gshgeom_exact_sum_add(&e, 1, 12, tiny);
  EXPECT_EQ(gshgeom_exact_sum_sign(&e), -1);
  /* a non-finite factor is a caller bug: the sum says 0 */
  GSHGEOM_ExactSum bad;
  gshgeom_exact_sum_init(&bad);
  volatile double zero = 0.0;
  double nanv[1] = {zero / zero};
  gshgeom_exact_sum_add(&bad, 0, 1, nanv);
  EXPECT_EQ(gshgeom_exact_sum_sign(&bad), 0);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
