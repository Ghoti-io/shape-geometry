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
 * The curve outline (stroke_curves.h): the tolerance and its range (SG34), the statuses (SG14, SG15), the inputs and
 * their refusals, the forms, which elements each stroke style produces (lines only where the region is polygons, exact
 * rational-quadratic arcs where it is round, cubics for a curved offset), the clip part as the map applied to the
 * centerline's closed subpaths (SG32), the hash and the dump (SG41), the limits, allocation failure, and the magnitudes
 * the tolerance can be certified for (SG35).  Geometry against the definition is test_stroke_curves_oracle.cpp's.
 */

#include "../test_helpers.h"

#include "../alloc_tracker.h"
#include "../curve_eval.h"
#include "../dump_sink.h"
#include "../stroke_matrix.h"

extern "C" {
#include "../../src/stroke/stroke_internal.h"
}

#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace sm = stroke_matrix;
namespace ce = curve_eval;

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

GSHGEOM_Result expand(const GSHGEOM_Path * p, const GSHGEOM_Stroke & s, const GSHGEOM_Map & m, double t,
    GSHGEOM_Stroke_Curves_Form f, GSHGEOM_Stroke_Curves ** o) {
  return gshgeom_stroke_curves_expand(p, &s, &m, t, f, nullptr, nullptr, o);
}

uint64_t curves_hash(const GSHGEOM_Path * p, const GSHGEOM_Stroke & s, const GSHGEOM_Map & m,
    GSHGEOM_Stroke_Curves_Form f = GSHGEOM_STROKE_CURVES_RATIONAL, double t = GSHGEOM_STROKE_CURVES_TOLERANCE) {
  GSHGEOM_Stroke_Curves * o = nullptr;
  EXPECT_EQ(expand(p, s, m, t, f, &o), GSHGEOM_OK);
  uint64_t h = gshgeom_stroke_curves_hash(o);
  gshgeom_stroke_curves_destroy(o);
  return h;
}

/** How many elements of each verb a path holds. */
std::map<int, int> verbs(const GSHGEOM_Path * p) {
  std::map<int, int> out;
  for (size_t i = 0; p && i < gshgeom_path_element_count(p); i++) {
    GSHGEOM_Element e;
    gshgeom_path_element(p, i, &e);
    out[(int)e.verb]++;
  }
  return out;
}

} // namespace

/* --- SG34: the tolerance parameter ---------------------------------------------------- */

TEST(StrokeCurvesInputs, TheToleranceIsCheckedAtBothEdgesOfItsRange) {
  GSHGEOM_Path * p = corner();
  GSHGEOM_Map m = identity();
  GSHGEOM_Stroke s = plain();
  GSHGEOM_Stroke_Curves * o = nullptr;
  EXPECT_EQ(GSHGEOM_STROKE_CURVES_MIN_TOLERANCE, std::ldexp(1.0, -14));
  EXPECT_EQ(GSHGEOM_STROKE_CURVES_MAX_TOLERANCE, std::ldexp(1.0, 4));
  EXPECT_EQ(GSHGEOM_STROKE_CURVES_TOLERANCE, std::ldexp(1.0, -10));
  for (double t : {GSHGEOM_STROKE_CURVES_MIN_TOLERANCE, std::ldexp(1.0, -12), GSHGEOM_STROKE_CURVES_TOLERANCE, 1.0,
           GSHGEOM_STROKE_CURVES_MAX_TOLERANCE}) {
    ASSERT_EQ(expand(p, s, m, t, GSHGEOM_STROKE_CURVES_RATIONAL, &o), GSHGEOM_OK) << t;
    EXPECT_EQ(gshgeom_stroke_curves_status(o), GSHGEOM_STROKE_DRAWN);
    gshgeom_stroke_curves_destroy(o);
  }
  o = nullptr;
  double inf = std::numeric_limits<double>::infinity(), nan = std::numeric_limits<double>::quiet_NaN();
  for (double t : {std::nextafter(GSHGEOM_STROKE_CURVES_MIN_TOLERANCE, 0.0),
           std::nextafter(GSHGEOM_STROKE_CURVES_MAX_TOLERANCE, 100.0), 0.0, -1.0, nan, inf, -inf}) {
    EXPECT_EQ(expand(p, s, m, t, GSHGEOM_STROKE_CURVES_RATIONAL, &o), GSHGEOM_ERR_INVALID) << t;
    EXPECT_EQ(o, nullptr);
  }
  gshgeom_path_destroy(p);
}

TEST(StrokeCurvesInputs, RefusalsAndTheFormRange) {
  GSHGEOM_Path * p = corner();
  GSHGEOM_Map m = identity();
  GSHGEOM_Stroke s = plain();
  GSHGEOM_Stroke_Curves * o = nullptr;
  double t = GSHGEOM_STROKE_CURVES_TOLERANCE;
  auto F = GSHGEOM_STROKE_CURVES_RATIONAL;
  EXPECT_EQ(gshgeom_stroke_curves_expand(nullptr, &s, &m, t, F, nullptr, nullptr, &o), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_stroke_curves_expand(p, nullptr, &m, t, F, nullptr, nullptr, &o), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_stroke_curves_expand(p, &s, nullptr, t, F, nullptr, nullptr, &o), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_stroke_curves_expand(p, &s, &m, t, F, nullptr, nullptr, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(expand(p, s, m, t, GSHGEOM_STROKE_CURVES_FORM_COUNT, &o), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(expand(p, s, m, t, (GSHGEOM_Stroke_Curves_Form)-1, &o), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(o, nullptr);
  GSHGEOM_Map m3;
  gshgeom_map_identity(3, &m3);
  EXPECT_EQ(expand(p, s, m3, t, F, &o), GSHGEOM_ERR_INVALID);
  GSHGEOM_Map bad = m;
  bad.m[0][0] = std::numeric_limits<double>::quiet_NaN();
  EXPECT_EQ(expand(p, s, bad, t, F, &o), GSHGEOM_ERR_INVALID);
  GSHGEOM_Stroke neg = s;
  neg.width = -1;
  EXPECT_EQ(expand(p, neg, m, t, F, &o), GSHGEOM_ERR_INVALID);
  GSHGEOM_Stroke prof = s;
  int dummy;
  prof.width_profile = &dummy;
  EXPECT_EQ(expand(p, prof, m, t, F, &o), GSHGEOM_ERR_UNSUPPORTED);
  GSHGEOM_Path * p3 = nullptr;
  ASSERT_EQ(gshgeom_path_create(nullptr, nullptr, 3, &p3), GSHGEOM_OK);
  GSHGEOM_Point a = {0, 0, 0}, b = {1, 1, 1};
  gshgeom_path_move_to(p3, a);
  gshgeom_path_line_to(p3, b);
  EXPECT_EQ(expand(p3, s, m, t, F, &o), GSHGEOM_ERR_UNSUPPORTED);
  gshgeom_path_destroy(p3);
  EXPECT_EQ(o, nullptr);
  /* a NULL outline reads as absent everywhere */
  EXPECT_EQ(gshgeom_stroke_curves_status(nullptr), GSHGEOM_STROKE_STATUS_COUNT);
  EXPECT_EQ(gshgeom_stroke_curves_form(nullptr), GSHGEOM_STROKE_CURVES_FORM_COUNT);
  EXPECT_EQ(gshgeom_stroke_curves_body(nullptr), nullptr);
  EXPECT_EQ(gshgeom_stroke_curves_aligned(nullptr), nullptr);
  EXPECT_EQ(gshgeom_stroke_curves_clip(nullptr), nullptr);
  EXPECT_EQ(gshgeom_stroke_curves_hash(nullptr), GSHGEOM_FNV_OFFSET_BASIS);
  gshgeom_stroke_curves_destroy(nullptr);
  EXPECT_EQ(gshgeom_stroke_curves_revision(), GSHGEOM_STROKE_CURVES_REVISION);
  EXPECT_EQ(GSHGEOM_STROKE_CURVES_REVISION, 1);
  gshgeom_path_destroy(p);
}

/* --- SG14, SG15: statuses ------------------------------------------------------------- */

TEST(StrokeCurvesStatus, ZeroWidthIsEmptyAndASingularMapIsNotInvertible) {
  GSHGEOM_Path * p = sm::build({{'M', {0, 0}}, {'L', {5, 0}}, {'M', {9, 9}}, {'Z', {}}});
  GSHGEOM_Map m = identity();
  GSHGEOM_Stroke s = plain();
  s.width = 0.0;
  s.cap = GSHGEOM_STROKE_CAP_ROUND;
  GSHGEOM_Stroke_Curves * o = nullptr;
  ASSERT_EQ(expand(p, s, m, GSHGEOM_STROKE_CURVES_TOLERANCE, GSHGEOM_STROKE_CURVES_CUBIC, &o), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_stroke_curves_status(o), GSHGEOM_STROKE_EMPTY);
  EXPECT_EQ(gshgeom_stroke_curves_body(o), nullptr);
  EXPECT_EQ(gshgeom_stroke_curves_aligned(o), nullptr);
  EXPECT_EQ(gshgeom_stroke_curves_clip(o), nullptr);
  EXPECT_EQ(gshgeom_stroke_curves_form(o), GSHGEOM_STROKE_CURVES_CUBIC);
  gshgeom_stroke_curves_destroy(o);
  s.width = 2.0;
  for (GSHGEOM_Map sing : {scale(0.0, 1.0), scale(1.0, 0.0)}) {
    for (uint32_t space : {0u, 1u}) {
      s.width_space = space;
      ASSERT_EQ(expand(p, s, sing, GSHGEOM_STROKE_CURVES_TOLERANCE, GSHGEOM_STROKE_CURVES_RATIONAL, &o), GSHGEOM_OK);
      EXPECT_EQ(gshgeom_stroke_curves_status(o), GSHGEOM_STROKE_NOT_INVERTIBLE);
      EXPECT_EQ(gshgeom_stroke_curves_body(o), nullptr);
      gshgeom_stroke_curves_destroy(o);
    }
  }
  /* the singular map is decided before zero width (SG15) */
  s.width = 0.0;
  ASSERT_EQ(expand(p, s, scale(0.0, 1.0), GSHGEOM_STROKE_CURVES_TOLERANCE, GSHGEOM_STROKE_CURVES_RATIONAL, &o), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_stroke_curves_status(o), GSHGEOM_STROKE_NOT_INVERTIBLE);
  gshgeom_stroke_curves_destroy(o);
  /* and a rank-one skew and a product that underflows to a singular matrix, as stroke.h's predicate says */
  double mm[4] = {1, 2, 2, 4}, tt[2] = {0, 0};
  GSHGEOM_Map skew;
  ASSERT_EQ(gshgeom_map_from_matrix(2, mm, tt, &skew), GSHGEOM_OK);
  s.width = 2.0;
  ASSERT_EQ(expand(p, s, skew, GSHGEOM_STROKE_CURVES_TOLERANCE, GSHGEOM_STROKE_CURVES_RATIONAL, &o), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_stroke_curves_status(o), GSHGEOM_STROKE_NOT_INVERTIBLE);
  gshgeom_stroke_curves_destroy(o);
  gshgeom_path_destroy(p);
}

/* --- SG3: canonical meaning ----------------------------------------------------------- */

TEST(StrokeCurvesCanonical, AStyleAndItsCanonicalFormGiveTheSameBits) {
  GSHGEOM_Path * p = corner();
  GSHGEOM_Map m = identity();
  GSHGEOM_Stroke a = plain(), b = plain();
  a.join = b.join = GSHGEOM_STROKE_JOIN_ROUND;
  a.miter_limit = 10.0; /* does not apply to a round join */
  EXPECT_EQ(curves_hash(p, a, m), curves_hash(p, b, m));
  double d[2] = {1.0, -1.0}, z[2] = {0.0, 0.0};
  GSHGEOM_Stroke c = plain();
  c.dashes = d;
  c.dash_count = 2;
  EXPECT_EQ(curves_hash(p, c, m), curves_hash(p, plain(), m));
  c.dashes = z;
  EXPECT_EQ(curves_hash(p, c, m), curves_hash(p, plain(), m));
  GSHGEOM_Stroke e = plain(), f = plain();
  e.cap = f.cap = GSHGEOM_STROKE_CAP_BUTT;
  f.degenerate = GSHGEOM_STROKE_DEGENERATE_OMIT;
  EXPECT_EQ(curves_hash(p, e, m), curves_hash(p, f, m));
  GSHGEOM_Stroke g = plain();
  g.width = -0.0 + 2.0;
  EXPECT_EQ(curves_hash(p, g, m), curves_hash(p, plain(), m));
  /* the tolerance is in no hash: a polygon of lines is the same at every tolerance */
  GSHGEOM_Path * poly = sm::build({{'M', {0, 0}}, {'L', {10, 0}}, {'L', {10, 10}}});
  EXPECT_EQ(curves_hash(poly, plain(), m, GSHGEOM_STROKE_CURVES_RATIONAL, 0x1p-14),
      curves_hash(poly, plain(), m, GSHGEOM_STROKE_CURVES_RATIONAL, 8.0));
  gshgeom_path_destroy(poly);
  gshgeom_path_destroy(p);
}

/* --- SG37: which elements a style produces -------------------------------------------- */

TEST(StrokeCurvesElements, APolygonalStrokeIsLinesOnlyAndARoundOneHasExactArcs) {
  GSHGEOM_Map m = identity();
  GSHGEOM_Path * poly = sm::build({{'M', {0, 0}}, {'L', {10, 0}}, {'L', {10, 10}}, {'L', {2, 3}}});
  for (int join : {GSHGEOM_STROKE_JOIN_MITER, GSHGEOM_STROKE_JOIN_BEVEL, GSHGEOM_STROKE_JOIN_MITER_CLIP}) {
    for (int cap : {GSHGEOM_STROKE_CAP_BUTT, GSHGEOM_STROKE_CAP_SQUARE}) {
      GSHGEOM_Stroke s = plain();
      s.join = (uint32_t)join;
      s.cap = (uint32_t)cap;
      GSHGEOM_Stroke_Curves * o = nullptr;
      ASSERT_EQ(expand(poly, s, m, GSHGEOM_STROKE_CURVES_TOLERANCE, GSHGEOM_STROKE_CURVES_RATIONAL, &o), GSHGEOM_OK);
      auto v = verbs(gshgeom_stroke_curves_body(o));
      EXPECT_EQ(v.count(GSHGEOM_VERB_CUBIC), 0u) << join << " " << cap;
      EXPECT_EQ(v.count(GSHGEOM_VERB_RATIONAL_QUADRATIC), 0u);
      EXPECT_EQ(v.count(GSHGEOM_VERB_QUADRATIC), 0u);
      EXPECT_GT(v[GSHGEOM_VERB_LINE], 8);
      gshgeom_stroke_curves_destroy(o);
    }
  }
  /* round joins and caps: arcs of at most a quarter turn, weights in [sqrt(1/2), 1], and no cubic */
  GSHGEOM_Stroke s = plain();
  s.join = GSHGEOM_STROKE_JOIN_ROUND;
  s.cap = GSHGEOM_STROKE_CAP_ROUND;
  GSHGEOM_Stroke_Curves * o = nullptr;
  ASSERT_EQ(expand(poly, s, m, GSHGEOM_STROKE_CURVES_TOLERANCE, GSHGEOM_STROKE_CURVES_RATIONAL, &o), GSHGEOM_OK);
  const GSHGEOM_Path * body = gshgeom_stroke_curves_body(o);
  auto v = verbs(body);
  EXPECT_EQ(v.count(GSHGEOM_VERB_CUBIC), 0u);
  EXPECT_GE(v[GSHGEOM_VERB_RATIONAL_QUADRATIC], 4);
  for (size_t i = 0; i < gshgeom_path_element_count(body); i++) {
    GSHGEOM_Element e;
    gshgeom_path_element(body, i, &e);
    if (e.verb == GSHGEOM_VERB_RATIONAL_QUADRATIC) {
      EXPECT_EQ(e.curve.w[0], 1.0);
      EXPECT_EQ(e.curve.w[2], 1.0);
      EXPECT_GE(e.curve.w[1], std::sqrt(0.5) * (1 - 1e-15));
      EXPECT_LE(e.curve.w[1], 1.0);
    }
  }
  /* and the cubic form has the same region as lines and cubics only */
  GSHGEOM_Stroke_Curves * oc = nullptr;
  ASSERT_EQ(expand(poly, s, m, GSHGEOM_STROKE_CURVES_TOLERANCE, GSHGEOM_STROKE_CURVES_CUBIC, &oc), GSHGEOM_OK);
  auto vc = verbs(gshgeom_stroke_curves_body(oc));
  EXPECT_EQ(vc.count(GSHGEOM_VERB_RATIONAL_QUADRATIC), 0u);
  EXPECT_GE(vc[GSHGEOM_VERB_CUBIC], v[GSHGEOM_VERB_RATIONAL_QUADRATIC]);
  EXPECT_EQ(gshgeom_stroke_curves_form(oc), GSHGEOM_STROKE_CURVES_CUBIC);
  gshgeom_stroke_curves_destroy(oc);
  gshgeom_stroke_curves_destroy(o);
  gshgeom_path_destroy(poly);
}

TEST(StrokeCurvesElements, ACurvedSegmentGivesCubicOffsetsAndEveryRingIsClosed) {
  GSHGEOM_Path * p = sm::build({{'M', {0, 0}}, {'C', {4, 3, 6, -3, 10, 0}}});
  GSHGEOM_Map m = identity();
  GSHGEOM_Stroke s = plain();
  GSHGEOM_Stroke_Curves * o = nullptr;
  ASSERT_EQ(expand(p, s, m, GSHGEOM_STROKE_CURVES_TOLERANCE, GSHGEOM_STROKE_CURVES_RATIONAL, &o), GSHGEOM_OK);
  const GSHGEOM_Path * body = gshgeom_stroke_curves_body(o);
  auto v = verbs(body);
  EXPECT_GE(v[GSHGEOM_VERB_CUBIC], 2);
  EXPECT_EQ(v[GSHGEOM_VERB_MOVE], v[GSHGEOM_VERB_CLOSE]); /* every subpath is closed */
  /* a gentle cubic's offsets run as one contour: a handful of subpaths for the body, caps and nothing else */
  EXPECT_LE(v[GSHGEOM_VERB_MOVE], 3);
  gshgeom_stroke_curves_destroy(o);
  gshgeom_path_destroy(p);
}

TEST(StrokeCurvesElements, AnArcsJoinWhoseBoundaryIsAPolarGraphIsExactArcs) {
  /* two curved segments turning at a vertex with butt caps: the only rational quadratics are the join's extensions,
   * circles of curvature written as arcs, not as the polylines a polygon expansion needs */
  GSHGEOM_Path * p = sm::build({{'M', {0, 0}}, {'Q', {6, 0, 8, 5}}, {'Q', {14, 6, 16, 12}}});
  GSHGEOM_Map m = identity();
  GSHGEOM_Stroke s = plain();
  s.width = 1.0;
  s.join = GSHGEOM_STROKE_JOIN_ARCS;
  s.miter_limit = 10.0;
  GSHGEOM_Stroke_Curves * o = nullptr;
  ASSERT_EQ(expand(p, s, m, GSHGEOM_STROKE_CURVES_TOLERANCE, GSHGEOM_STROKE_CURVES_RATIONAL, &o), GSHGEOM_OK);
  auto v = verbs(gshgeom_stroke_curves_body(o));
  EXPECT_GE(v[GSHGEOM_VERB_RATIONAL_QUADRATIC], 1);
  /* a polyline boundary would need dozens of lines */
  EXPECT_LE(v[GSHGEOM_VERB_LINE], 24);
  gshgeom_stroke_curves_destroy(o);
  gshgeom_path_destroy(p);
}

TEST(StrokeCurvesElements, AnEllipticPenKeepsTheWeightsOfTheCircularOne) {
  /* local width under scale(8, 1): the round cap's arcs are ellipse arcs with the circle's weights (SG36) */
  GSHGEOM_Path * p = sm::build({{'M', {0, 0}}, {'L', {10, 0}}});
  GSHGEOM_Stroke s = plain();
  s.cap = GSHGEOM_STROKE_CAP_ROUND;
  std::vector<double> w[2];
  for (int k = 0; k < 2; k++) {
    GSHGEOM_Map m = k == 0 ? identity() : scale(8.0, 1.0);
    GSHGEOM_Stroke_Curves * o = nullptr;
    ASSERT_EQ(expand(p, s, m, GSHGEOM_STROKE_CURVES_TOLERANCE, GSHGEOM_STROKE_CURVES_RATIONAL, &o), GSHGEOM_OK);
    const GSHGEOM_Path * body = gshgeom_stroke_curves_body(o);
    for (size_t i = 0; i < gshgeom_path_element_count(body); i++) {
      GSHGEOM_Element e;
      gshgeom_path_element(body, i, &e);
      if (e.verb == GSHGEOM_VERB_RATIONAL_QUADRATIC) {
        w[k].push_back(e.curve.w[1]);
      }
    }
    gshgeom_stroke_curves_destroy(o);
  }
  ASSERT_EQ(w[0].size(), 4u);
  EXPECT_EQ(w[0], w[1]);
  gshgeom_path_destroy(p);
}

TEST(StrokeCurvesElements, AMirroringMapKeepsEveryRingPositive) {
  GSHGEOM_Path * p = sm::build({{'M', {0, 0}}, {'L', {6, 0}}, {'C', {9, 0, 9, 6, 6, 8}}, {'L', {1, 8}}});
  GSHGEOM_Stroke s = plain();
  s.join = GSHGEOM_STROKE_JOIN_ROUND;
  s.cap = GSHGEOM_STROKE_CAP_ROUND;
  for (GSHGEOM_Map m : {scale(-1.0, 1.0), scale(1.0, -1.0), scale(-2.0, -2.0), scale(2.0, 1.0)}) {
    for (uint32_t space : {0u, 1u}) {
      s.width_space = space;
      GSHGEOM_Stroke_Curves * o = nullptr;
      ASSERT_EQ(expand(p, s, m, GSHGEOM_STROKE_CURVES_TOLERANCE, GSHGEOM_STROKE_CURVES_RATIONAL, &o), GSHGEOM_OK);
      ce::Winding wd;
      wd.build(gshgeom_stroke_curves_body(o), 1e-5L);
      /* a point on the centerline is inside, and no ring is negative anywhere on a grid */
      ce::P2 on = {m.m[0][0] * 3.0, m.m[1][1] * 0.0};
      EXPECT_EQ(wd.at(on), 1);
      for (double x = -20; x <= 20; x += 0.7) {
        for (double y = -20; y <= 20; y += 0.7) {
          int least = 0;
          int w = wd.at({x, y}, &least);
          EXPECT_GE(w, 0);
          EXPECT_GE(least, 0);
        }
      }
      gshgeom_stroke_curves_destroy(o);
    }
  }
  gshgeom_path_destroy(p);
}

TEST(StrokeCurvesElements, TheBudgetIsDividedByTheMapsNormUnderLocalWidth) {
  /* the same model curve under a map sixteen times larger: its offsets are sixteen times larger in target units, the
   * tolerance is the same, so the build space (the model) must be fitted sixteen times finer, and the count of pieces,
   * which falls as the fourth root of the budget, doubles (SG36, E5) */
  GSHGEOM_Path * p = sm::build({{'M', {0, 0}}, {'C', {3, 5, 7, -5, 10, 0}}});
  GSHGEOM_Stroke s = plain();
  s.width = 1.0;
  int count[2] = {0, 0};
  for (int k = 0; k < 2; k++) {
    GSHGEOM_Map m = k == 0 ? identity() : scale(16.0, 16.0);
    GSHGEOM_Stroke_Curves * o = nullptr;
    ASSERT_EQ(expand(p, s, m, GSHGEOM_STROKE_CURVES_TOLERANCE, GSHGEOM_STROKE_CURVES_RATIONAL, &o), GSHGEOM_OK);
    count[k] = verbs(gshgeom_stroke_curves_body(o))[GSHGEOM_VERB_CUBIC];
    gshgeom_stroke_curves_destroy(o);
  }
  EXPECT_GT(count[0], 2);
  EXPECT_GE(count[1], count[0] * 3 / 2) << count[0] << " against " << count[1];
  gshgeom_path_destroy(p);
}

/* --- SG32: the parts, and the clip as the map applied --------------------------------- */

TEST(StrokeCurvesParts, AlignmentGivesBodyAlignedAndTheClipIsTheMapApplied) {
  /* a closed subpath with a curve and a rational arc, and an open one: the clip is the closed one, element for element */
  GSHGEOM_Path * p = nullptr;
  ASSERT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &p), GSHGEOM_OK);
  GSHGEOM_Point a = {0, 0, 0}, b = {8, 0, 0}, c = {8, 6, 0}, d = {4, 9, 0}, e = {0, 6, 0};
  double w3[3] = {1, 0.8, 1}, w4[4] = {1, 2, 0.5, 1};
  gshgeom_path_move_to(p, a);
  gshgeom_path_line_to(p, b);
  gshgeom_path_quadratic_to(p, GSHGEOM_Point{10, 3, 0}, c);
  gshgeom_path_rational_quadratic_to(p, GSHGEOM_Point{7, 8, 0}, d, w3);
  gshgeom_path_rational_cubic_to(p, GSHGEOM_Point{3, 10, 0}, GSHGEOM_Point{1, 9, 0}, e, w4);
  gshgeom_path_close(p);
  gshgeom_path_move_to(p, GSHGEOM_Point{20, 0, 0});
  gshgeom_path_line_to(p, GSHGEOM_Point{26, 4, 0});
  double mm[4] = {0.6, -0.8, 0.8, 0.6}, tt[2] = {3, 4};
  GSHGEOM_Map m;
  ASSERT_EQ(gshgeom_map_from_matrix(2, mm, tt, &m), GSHGEOM_OK);
  for (uint32_t align : {0u, 1u, 2u}) {
    for (int rule : {0, 1}) {
      GSHGEOM_Stroke s = plain();
      s.alignment = align;
      s.alignment_rule = (GSHGEOM_Fill_Rule)rule;
      GSHGEOM_Stroke_Curves * o = nullptr;
      ASSERT_EQ(expand(p, s, m, GSHGEOM_STROKE_CURVES_TOLERANCE, GSHGEOM_STROKE_CURVES_CUBIC, &o), GSHGEOM_OK);
      EXPECT_NE(gshgeom_stroke_curves_body(o), nullptr);
      EXPECT_EQ(gshgeom_stroke_curves_aligned(o) != nullptr, align != 0);
      EXPECT_EQ(gshgeom_stroke_curves_clip(o) != nullptr, align != 0);
      EXPECT_EQ(gshgeom_stroke_curves_clip_complement(o), align == 2 ? 1 : 0);
      EXPECT_EQ(gshgeom_stroke_curves_clip_rule(o), align == 0 ? GSHGEOM_FILL_NONZERO : (GSHGEOM_Fill_Rule)rule);
      if (align) {
        const GSHGEOM_Path * clip = gshgeom_stroke_curves_clip(o);
        /* the clip is not lowered even in the cubic form: it keeps the input's own element kinds */
        auto v = verbs(clip);
        EXPECT_EQ(v[GSHGEOM_VERB_QUADRATIC], 1);
        EXPECT_EQ(v[GSHGEOM_VERB_RATIONAL_QUADRATIC], 1);
        EXPECT_EQ(v[GSHGEOM_VERB_RATIONAL_CUBIC], 1);
        EXPECT_EQ(v[GSHGEOM_VERB_MOVE], 1);
        ASSERT_EQ(gshgeom_path_element_count(clip), 6u);
        /* every control point is the map applied to the input's, computed here independently */
        for (size_t i = 0; i < 6; i++) {
          GSHGEOM_Element ei, eo;
          gshgeom_path_element(p, i, &ei);
          gshgeom_path_element(clip, i, &eo);
          ASSERT_EQ(ei.verb, eo.verb);
          int np = gshgeom_curve_point_count(ei.curve.kind);
          for (int k = 0; k < np; k++) {
            double x = mm[0] * ei.curve.p[k].x + mm[1] * ei.curve.p[k].y + tt[0];
            double y = mm[2] * ei.curve.p[k].x + mm[3] * ei.curve.p[k].y + tt[1];
            EXPECT_NEAR(eo.curve.p[k].x, x, 1e-13) << i << " " << k;
            EXPECT_NEAR(eo.curve.p[k].y, y, 1e-13) << i << " " << k;
            if (gshgeom_curve_is_rational(ei.curve.kind)) {
              EXPECT_EQ(eo.curve.w[k], ei.curve.w[k]);
            }
          }
        }
      }
      gshgeom_stroke_curves_destroy(o);
    }
  }
  /* no closed subpath with extent: no aligned part, no clip */
  GSHGEOM_Path * open = corner();
  GSHGEOM_Stroke s = plain();
  s.alignment = GSHGEOM_STROKE_ALIGN_INSIDE;
  GSHGEOM_Stroke_Curves * o = nullptr;
  ASSERT_EQ(expand(open, s, m, GSHGEOM_STROKE_CURVES_TOLERANCE, GSHGEOM_STROKE_CURVES_RATIONAL, &o), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_stroke_curves_aligned(o), nullptr);
  EXPECT_EQ(gshgeom_stroke_curves_clip(o), nullptr);
  gshgeom_stroke_curves_destroy(o);
  gshgeom_path_destroy(open);
  gshgeom_path_destroy(p);
}

/* --- SG41: the hash and the dump ------------------------------------------------------ */

TEST(StrokeCurvesHash, TheSpellingIsRestatedFieldByField) {
  GSHGEOM_Path * p = sm::build({{'M', {0, 0}}, {'L', {8, 0}}, {'L', {8, 8}}, {'Z', {}}, {'M', {12, 0}}, {'L', {16, 4}}});
  GSHGEOM_Map m = identity();
  for (uint32_t align : {0u, 2u}) {
    GSHGEOM_Stroke s = plain();
    s.alignment = align;
    s.alignment_rule = GSHGEOM_FILL_EVEN_ODD;
    s.join = GSHGEOM_STROKE_JOIN_ROUND;
    GSHGEOM_Stroke_Curves * o = nullptr;
    ASSERT_EQ(expand(p, s, m, GSHGEOM_STROKE_CURVES_TOLERANCE, GSHGEOM_STROKE_CURVES_RATIONAL, &o), GSHGEOM_OK);
    GSHGEOM_Hash h;
    gshgeom_hash_init(&h);
    gshgeom_hash_u32(&h, 1);
    gshgeom_hash_u32(&h, 0x43525453u);
    gshgeom_hash_u8(&h, 0);
    gshgeom_hash_u8(&h, (uint8_t)align);
    gshgeom_hash_u8(&h, align ? 1 : 0);
    gshgeom_hash_u8(&h, align == 2 ? 1 : 0);
    gshgeom_hash_u64(&h, gshgeom_path_hash(gshgeom_stroke_curves_body(o)));
    gshgeom_hash_u64(&h, gshgeom_path_hash(gshgeom_stroke_curves_aligned(o)));
    gshgeom_hash_u64(&h, gshgeom_path_hash(gshgeom_stroke_curves_clip(o)));
    EXPECT_EQ(gshgeom_stroke_curves_hash(o), gshgeom_hash_value(&h));
    gshgeom_stroke_curves_destroy(o);
  }
  /* an absent part hashes as the offset basis */
  GSHGEOM_Stroke s = plain();
  s.width = 0.0;
  GSHGEOM_Stroke_Curves * o = nullptr;
  ASSERT_EQ(expand(p, s, m, GSHGEOM_STROKE_CURVES_TOLERANCE, GSHGEOM_STROKE_CURVES_RATIONAL, &o), GSHGEOM_OK);
  GSHGEOM_Hash h;
  gshgeom_hash_init(&h);
  gshgeom_hash_u32(&h, 1);
  gshgeom_hash_u32(&h, 0x43525453u);
  gshgeom_hash_u8(&h, 1);
  gshgeom_hash_u8(&h, 0);
  gshgeom_hash_u8(&h, 0);
  gshgeom_hash_u8(&h, 0);
  for (int k = 0; k < 3; k++) {
    gshgeom_hash_u64(&h, GSHGEOM_FNV_OFFSET_BASIS);
  }
  EXPECT_EQ(gshgeom_stroke_curves_hash(o), gshgeom_hash_value(&h));
  gshgeom_stroke_curves_destroy(o);
  /* the same input gives the same bits, and the polygon outline's hash is a different statement */
  GSHGEOM_Path * q = corner();
  GSHGEOM_Stroke t = plain();
  t.join = GSHGEOM_STROKE_JOIN_ARCS;
  t.cap = GSHGEOM_STROKE_CAP_ROUND;
  GSHGEOM_Map mm = scale(1.5, 0.75);
  EXPECT_EQ(curves_hash(q, t, mm), curves_hash(q, t, mm));
  gshgeom_path_destroy(q);
  gshgeom_path_destroy(p);
}

TEST(StrokeCurvesDump, TheTextNamesEveryPartAndEveryWriteIsChecked) {
  GSHGEOM_Path * p = sm::build({{'M', {0, 0}}, {'L', {8, 0}}, {'L', {8, 8}}, {'Z', {}}});
  GSHGEOM_Map m = identity();
  GSHGEOM_Stroke s = plain();
  s.alignment = GSHGEOM_STROKE_ALIGN_INSIDE;
  GSHGEOM_Stroke_Curves * o = nullptr;
  ASSERT_EQ(expand(p, s, m, GSHGEOM_STROKE_CURVES_TOLERANCE, GSHGEOM_STROKE_CURVES_CUBIC, &o), GSHGEOM_OK);
  char * buf = nullptr;
  size_t len = 0;
  FILE * f = open_memstream(&buf, &len);
  ASSERT_NE(f, nullptr);
  ASSERT_EQ(gshgeom_stroke_curves_dump(o, f), GSHGEOM_OK);
  std::fclose(f);
  std::string text(buf, len);
  std::free(buf);
  char head[160];
  std::snprintf(head, sizeof head, "gshgeom.stroke-curves 1 status=0 form=1 align=1 rule=0 complement=0 hash=%016llx\n",
      (unsigned long long)gshgeom_stroke_curves_hash(o));
  EXPECT_EQ(text.compare(0, std::strlen(head), head), 0) << text.substr(0, 200);
  EXPECT_NE(text.find("part body\ngshgeom.path 1"), std::string::npos);
  EXPECT_NE(text.find("part aligned\ngshgeom.path 1"), std::string::npos);
  EXPECT_NE(text.find("part clip\ngshgeom.path 1"), std::string::npos);
  EXPECT_EQ(gshgeom_stroke_curves_dump(nullptr, stdout), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_stroke_curves_dump(o, nullptr), GSHGEOM_ERR_INVALID);
  GSHGEOM_EVERY_WRITE_CHECKED(gshgeom_stroke_curves_dump, o, "curves dump");
  gshgeom_stroke_curves_destroy(o);
  s.alignment = GSHGEOM_STROKE_ALIGN_CENTER;
  ASSERT_EQ(expand(p, s, m, GSHGEOM_STROKE_CURVES_TOLERANCE, GSHGEOM_STROKE_CURVES_RATIONAL, &o), GSHGEOM_OK);
  GSHGEOM_EVERY_WRITE_CHECKED(gshgeom_stroke_curves_dump, o, "curves dump, absent parts last");
  gshgeom_stroke_curves_destroy(o);
  s.width = 0.0;
  ASSERT_EQ(expand(p, s, m, GSHGEOM_STROKE_CURVES_TOLERANCE, GSHGEOM_STROKE_CURVES_RATIONAL, &o), GSHGEOM_OK);
  f = open_memstream(&buf, &len);
  ASSERT_EQ(gshgeom_stroke_curves_dump(o, f), GSHGEOM_OK);
  std::fclose(f);
  text.assign(buf, len);
  std::free(buf);
  EXPECT_NE(text.find("part body -\npart aligned -\npart clip -\n"), std::string::npos) << text;
  gshgeom_stroke_curves_destroy(o);
  gshgeom_path_destroy(p);
}

/* --- limits, magnitudes and allocation ------------------------------------------------ */

TEST(StrokeCurvesLimits, AnElementCountPastMaxElementsAndAForcedPieceAreRefused) {
  GSHGEOM_Path * p = corner();
  GSHGEOM_Stroke s = plain();
  s.join = GSHGEOM_STROKE_JOIN_ROUND;
  GSHGEOM_Map m = identity();
  GSHGEOM_Stroke_Curves * o = nullptr;
  GSHGEOM_Limits lim;
  gshgeom_limits_default(&lim);
  lim.max_elements = 12;
  EXPECT_EQ(gshgeom_stroke_curves_expand(p, &s, &m, GSHGEOM_STROKE_CURVES_TOLERANCE, GSHGEOM_STROKE_CURVES_RATIONAL, nullptr,
                &lim, &o), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(o, nullptr);
  gshgeom_path_destroy(p);
  /* with the depth cap at 2 a curve that needs more splits is ERR_LIMIT and returns nothing (E11) */
  GSHGEOM_Path * q = sm::build({{'M', {0, 0}}, {'C', {40, 60, 80, -60, 120, 0}}});
  s.width = 12;
  s.join = GSHGEOM_STROKE_JOIN_MITER;
  EXPECT_EQ(gshgeom_stroke_curves_expand_capped(q, &s, &m, GSHGEOM_STROKE_CURVES_TOLERANCE,
                GSHGEOM_STROKE_CURVES_RATIONAL, nullptr, nullptr, 2, &o), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(o, nullptr);
  ASSERT_EQ(gshgeom_stroke_curves_expand_capped(q, &s, &m, GSHGEOM_STROKE_CURVES_TOLERANCE,
                GSHGEOM_STROKE_CURVES_RATIONAL, nullptr, nullptr, GSHGEOM_STROKE_MAX_DEPTH, &o), GSHGEOM_OK);
  gshgeom_stroke_curves_destroy(o);
  gshgeom_path_destroy(q);
}

TEST(StrokeCurvesLimits, TheMagnitudeTheToleranceCanBeCertifiedForIsAStatedBound) {
  /* SG35: the arithmetic's reach, 2^-43 of the largest magnitude, must be within t / 8: coordinates up to 2^30 at the
   * default tolerance (t 2^40 = 2^30), up to 2^26 at the finest */
  GSHGEOM_Stroke s = plain();
  GSHGEOM_Map m = identity();
  GSHGEOM_Stroke_Curves * o = nullptr;
  GSHGEOM_Path * at26 = sm::build({{'M', {0x1p26 - 8, 0}}, {'L', {0x1p26 - 4, 0}}});
  EXPECT_EQ(expand(at26, s, m, GSHGEOM_STROKE_CURVES_MIN_TOLERANCE, GSHGEOM_STROKE_CURVES_RATIONAL, &o), GSHGEOM_OK);
  gshgeom_stroke_curves_destroy(o);
  o = nullptr;
  GSHGEOM_Path * at27 = sm::build({{'M', {0x1p27, 0}}, {'L', {0x1p27 + 4, 0}}});
  EXPECT_EQ(expand(at27, s, m, GSHGEOM_STROKE_CURVES_MIN_TOLERANCE, GSHGEOM_STROKE_CURVES_RATIONAL, &o), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(o, nullptr);
  EXPECT_EQ(expand(at27, s, m, GSHGEOM_STROKE_CURVES_TOLERANCE, GSHGEOM_STROKE_CURVES_RATIONAL, &o), GSHGEOM_OK);
  gshgeom_stroke_curves_destroy(o);
  GSHGEOM_Path * at31 = sm::build({{'M', {0x1p31, 0}}, {'L', {0x1p31 + 4, 0}}});
  EXPECT_EQ(expand(at31, s, m, GSHGEOM_STROKE_CURVES_MAX_TOLERANCE, GSHGEOM_STROKE_CURVES_RATIONAL, &o), GSHGEOM_ERR_LIMIT);
  gshgeom_path_destroy(at26);
  gshgeom_path_destroy(at27);
  gshgeom_path_destroy(at31);
}

TEST(StrokeCurvesAllocation, EveryFailureIsReportedAndNothingLeaks) {
  GSHGEOM_Path * p = sm::build({{'M', {0, 0}}, {'L', {8, 0}}, {'Q', {12, 0, 12, 4}}, {'L', {8, 8}}, {'Z', {}},
      {'M', {20, 20}}, {'Z', {}}});
  GSHGEOM_Map m = identity();
  GSHGEOM_Stroke s = plain();
  s.join = GSHGEOM_STROKE_JOIN_ROUND;
  s.cap = GSHGEOM_STROKE_CAP_ROUND;
  s.alignment = GSHGEOM_STROKE_ALIGN_INSIDE;
  double dashes[2] = {5, 2};
  s.dashes = dashes;
  s.dash_count = 2;
  tracker::Tracker count;
  {
    GSHGEOM_Allocator a = tracker::allocator_for(&count);
    GSHGEOM_Stroke_Curves * o = nullptr;
    ASSERT_EQ(gshgeom_stroke_curves_expand(p, &s, &m, GSHGEOM_STROKE_CURVES_TOLERANCE, GSHGEOM_STROKE_CURVES_CUBIC, &a,
                  nullptr, &o), GSHGEOM_OK);
    gshgeom_stroke_curves_destroy(o);
    EXPECT_TRUE(count.live.empty());
  }
  ASSERT_GT(count.calls, 5);
  long step = count.calls > 400 ? count.calls / 400 : 1;
  for (long k = 0; k < count.calls; k += step) {
    for (int persistent = 0; persistent < 2; persistent++) {
      tracker::Tracker t;
      t.fail_from = k;
      t.persistent = persistent != 0;
      GSHGEOM_Allocator a = tracker::allocator_for(&t);
      GSHGEOM_Stroke_Curves * o = nullptr;
      EXPECT_EQ(gshgeom_stroke_curves_expand(p, &s, &m, GSHGEOM_STROKE_CURVES_TOLERANCE, GSHGEOM_STROKE_CURVES_CUBIC, &a,
                    nullptr, &o), GSHGEOM_ERR_OOM) << k;
      EXPECT_EQ(o, nullptr);
      EXPECT_TRUE(t.live.empty()) << "k=" << k << " leaked";
    }
  }
  gshgeom_path_destroy(p);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
