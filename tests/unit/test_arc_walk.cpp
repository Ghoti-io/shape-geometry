/**
 * @file
 *
 * What reads a path of arcs (arc.h): the bounds (AR10), the map and the circular flag (AR11), the measurement
 * (AR12), the stroke expansion and the dashes (AR13) and the path vertices with the directions at them (AR14),
 * against tools/arc-oracle.py's vectors and against the same geometry written as rational quadratics.
 *
 * Copyright 2026 by Corey Pennycuff
 */

#include "../test_helpers.h"

#include "../alloc_tracker.h"
#include "../arc_vectors.h"
#include "../stroke_oracle.h"

extern "C" {
#include "../../src/stroke/stroke_internal.h"
}

#include <ghoti.io/math/basic.h>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>

namespace {

GSHGEOM_Point P(double x, double y) {
  return GSHGEOM_Point{x, y, 0.0};
}

uint64_t bits(double v) {
  uint64_t b;
  std::memcpy(&b, &v, sizeof b);
  return b;
}

GSHGEOM_Path * make() {
  GSHGEOM_Path * p = nullptr;
  EXPECT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &p), GSHGEOM_OK);
  return p;
}

// the arc of the centre form, beginning a subpath at its start (or continuing from the current point)
void add_form(GSHGEOM_Path * p, GSHGEOM_Point c, GSHGEOM_Point u, GSHGEOM_Point v, int unit, double t0, double sweep,
    unsigned flags = 0) {
  GSHGEOM_Arc_Form f = {c, u, v, (GSHGEOM_Angle_Unit)unit, t0, sweep, flags};
  ASSERT_EQ(gshgeom_path_arc_conjugate_to(p, &f, nullptr), GSHGEOM_OK);
}

GSHGEOM_Path * form_path(const arc_vectors::Box & b) {
  GSHGEOM_Path * p = make();
  add_form(p, P(b.cx, b.cy), P(b.ux, b.uy), P(b.vx, b.vy), b.unit, b.t0, b.sweep);
  return p;
}

GSHGEOM_Arc arc_at(const GSHGEOM_Path * p, size_t i) {
  GSHGEOM_Arc a;
  EXPECT_EQ(gshgeom_path_arc(p, i, &a), GSHGEOM_OK);
  return a;
}

GSHGEOM_Arc_Form form_of(const GSHGEOM_Arc & a) {
  return GSHGEOM_Arc_Form{a.center, a.u, a.v, a.unit, a.t0, a.sweep, a.flags};
}

GSHGEOM_Stroke plain() {
  GSHGEOM_Stroke s;
  gshgeom_stroke_init(&s);
  s.width = 4.0;
  return s;
}

uint64_t stroke_hash(const GSHGEOM_Path * p, const GSHGEOM_Stroke & s, const GSHGEOM_Map & m, GSHGEOM_Result * rr = nullptr) {
  GSHGEOM_Stroke_Outline * o = nullptr;
  GSHGEOM_Result r = gshgeom_stroke_expand(p, &s, &m, nullptr, nullptr, &o);
  if (rr) {
    *rr = r;
  } else {
    EXPECT_EQ(r, GSHGEOM_OK);
  }
  uint64_t h = r == GSHGEOM_OK ? gshgeom_stroke_outline_hash(o) : 0;
  gshgeom_stroke_outline_destroy(o);
  return h;
}

GSHGEOM_Path * lowered(const GSHGEOM_Path * p) {
  GSHGEOM_Path * l = nullptr;
  EXPECT_EQ(gshgeom_path_lower_arcs(p, nullptr, nullptr, &l), GSHGEOM_OK);
  return l;
}

GSHGEOM_Map identity() {
  GSHGEOM_Map m;
  EXPECT_EQ(gshgeom_map_identity(2, &m), GSHGEOM_OK);
  return m;
}

GSHGEOM_Map matrix(double a, double b, double c, double d, double tx = 0, double ty = 0) {
  const double m[4] = {a, b, c, d};
  const double t[2] = {tx, ty};
  GSHGEOM_Map out;
  EXPECT_EQ(gshgeom_map_from_matrix(2, m, t, &out), GSHGEOM_OK);
  return out;
}

const int joins[] = {GSHGEOM_STROKE_JOIN_MITER, GSHGEOM_STROKE_JOIN_ROUND, GSHGEOM_STROKE_JOIN_BEVEL,
    GSHGEOM_STROKE_JOIN_MITER_CLIP, GSHGEOM_STROKE_JOIN_ARCS};
const int caps[] = {GSHGEOM_STROKE_CAP_BUTT, GSHGEOM_STROKE_CAP_ROUND, GSHGEOM_STROKE_CAP_SQUARE};

}  // namespace

// ---------------------------------------------------------------------------
// AR10 bounds
// ---------------------------------------------------------------------------

TEST(ArcBounds, EveryVectorsBoxHoldsWithinTheDecompositionBound) {
  for (size_t i = 0; i < arc_vectors::box_count; i++) {
    const arc_vectors::Box & b = arc_vectors::boxes[i];
    GSHGEOM_Path * p = form_path(b);
    GSHGEOM_Bounds got;
    ASSERT_EQ(gshgeom_path_bounds(p, &got), GSHGEOM_OK) << i;
    const double tol = GSHGEOM_ARC_DECOMP * b.scale;
    EXPECT_NEAR(got.min.x, b.minx, tol) << i;
    EXPECT_NEAR(got.min.y, b.miny, tol) << i;
    EXPECT_NEAR(got.max.x, b.maxx, tol) << i;
    EXPECT_NEAR(got.max.y, b.maxy, tol) << i;
    // the ends are always inside
    GSHGEOM_Arc a = arc_at(p, 1);
    EXPECT_LE(got.min.x, std::fmin(a.start.x, a.end.x)) << i;
    EXPECT_GE(got.max.x, std::fmax(a.start.x, a.end.x)) << i;
    gshgeom_path_destroy(p);
  }
}

TEST(ArcBounds, ACircleThroughAnExtremumAndATinyArc) {
  // the unit circle's quarter from 10 to 100 degrees reaches y = 1 between its ends
  GSHGEOM_Path * p = make();
  add_form(p, P(0, 0), P(1, 0), P(0, 1), GSHGEOM_ANGLE_DEGREES, 10, 90, GSHGEOM_ARC_CIRCULAR);
  GSHGEOM_Bounds b;
  ASSERT_EQ(gshgeom_path_bounds(p, &b), GSHGEOM_OK);
  EXPECT_NEAR(b.max.y, 1.0, 0x1p-44);
  EXPECT_GT(b.max.y, 1.0 - 0x1p-44);
  gshgeom_path_destroy(p);
}

// ---------------------------------------------------------------------------
// AR11 maps and the circular flag
// ---------------------------------------------------------------------------

TEST(ArcMap, EveryVectorsMappedEllipseHoldsWithinTheDecompositionBound) {
  for (size_t i = 0; i < arc_vectors::mapped_count; i++) {
    const arc_vectors::Mapped & v = arc_vectors::mappeds[i];
    GSHGEOM_Path * p = make();
    add_form(p, P(v.cx, v.cy), P(v.ux, v.uy), P(v.vx, v.vy), v.unit, v.t0, v.sweep);
    GSHGEOM_Arc a = arc_at(p, 1);
    GSHGEOM_Map m = matrix(v.m00, v.m01, v.m10, v.m11, v.tx, v.ty);
    GSHGEOM_Arc out;
    ASSERT_EQ(gshgeom_arc_map(&m, &a, &out), GSHGEOM_OK) << i;
    const GSHGEOM_Arc_Form f = form_of(out);
    for (int k = 0; k < 5; k++) {
      GSHGEOM_Point q;
      ASSERT_EQ(gshgeom_arc_point(&f, v.t0 + v.sweep * k / 4.0, &q), GSHGEOM_OK);
      EXPECT_NEAR(q.x, v.pts[k][0], GSHGEOM_ARC_DECOMP * v.scale) << i << " " << k;
      EXPECT_NEAR(q.y, v.pts[k][1], GSHGEOM_ARC_DECOMP * v.scale) << i << " " << k;
    }
    // the ends are the map's images of the ends, to the bit; the parameters do not change
    GSHGEOM_Point s, e;
    ASSERT_EQ(gshgeom_map_apply_point(&m, a.start, &s), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_map_apply_point(&m, a.end, &e), GSHGEOM_OK);
    EXPECT_EQ(bits(out.start.x), bits(s.x));
    EXPECT_EQ(bits(out.start.y), bits(s.y));
    EXPECT_EQ(bits(out.end.x), bits(e.x));
    EXPECT_EQ(bits(out.end.y), bits(e.y));
    EXPECT_EQ(bits(out.t0), bits(a.t0));
    EXPECT_EQ(bits(out.sweep), bits(a.sweep));
    EXPECT_EQ(out.unit, a.unit);
    gshgeom_path_destroy(p);
  }
}

TEST(ArcMap, SimilaritiesKeepTheCircularFlagAndTheRestClearIt) {
  GSHGEOM_Path * p = make();
  add_form(p, P(3, 4), P(5, 0), P(0, 5), GSHGEOM_ANGLE_DEGREES, 0, 100, GSHGEOM_ARC_CIRCULAR);
  GSHGEOM_Arc a = arc_at(p, 1);
  ASSERT_NE(a.flags & GSHGEOM_ARC_CIRCULAR, 0u);
  struct Case {
    const char * name;
    GSHGEOM_Map m;
    bool keeps;
  };
  GSHGEOM_Map t, sc, mir, flip, quarter, rot, nonu, skew, mirnonu;
  ASSERT_EQ(gshgeom_map_translate(2, 5, -2, 0, &t), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_map_scale(2, 2, 2, 1, &sc), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_map_scale(2, -1, 1, 1, &mir), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_map_flip_y(2, &flip), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_map_rotate(GSHGEOM_ANGLE_DEGREES, 90, &quarter), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_map_rotate(GSHGEOM_ANGLE_DEGREES, 30, &rot), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_map_scale(2, 3, 0.5, 1, &nonu), GSHGEOM_OK);
  skew = matrix(1, 0.75, 0, 1);
  ASSERT_EQ(gshgeom_map_scale(2, -3, 1, 1, &mirnonu), GSHGEOM_OK);
  const Case cases[] = {{"translate", t, true}, {"scale", sc, true}, {"mirror", mir, true}, {"flip y", flip, true},
      {"quarter turn", quarter, true}, {"rotation by 30", rot, true}, {"non-uniform scale", nonu, false},
      {"skew", skew, false}, {"mirrored non-uniform scale", mirnonu, false}};
  for (const Case & c : cases) {
    GSHGEOM_Arc out;
    ASSERT_EQ(gshgeom_arc_map(&c.m, &a, &out), GSHGEOM_OK) << c.name;
    EXPECT_EQ((out.flags & GSHGEOM_ARC_CIRCULAR) != 0u, c.keeps) << c.name;
    GSHGEOM_Point centre;
    double radius;
    int ccw = -1;
    GSHGEOM_Result r = gshgeom_arc_circle(&out, &centre, &radius, &ccw);
    EXPECT_EQ(r, c.keeps ? GSHGEOM_OK : GSHGEOM_ERR_INVALID) << c.name;
  }
  gshgeom_path_destroy(p);
}

TEST(ArcMap, AMirrorFlipsTheDirectionAndTheCircleSaysWhichWay) {
  GSHGEOM_Path * p = make();
  add_form(p, P(0, 0), P(5, 0), P(0, 5), GSHGEOM_ANGLE_DEGREES, 0, 90, GSHGEOM_ARC_CIRCULAR);
  GSHGEOM_Arc a = arc_at(p, 1);
  GSHGEOM_Point c;
  double r;
  int ccw = -1;
  ASSERT_EQ(gshgeom_arc_circle(&a, &c, &r, &ccw), GSHGEOM_OK);
  EXPECT_EQ(ccw, 1);
  EXPECT_EQ(r, 5.0);
  EXPECT_EQ(c.x, 0.0);
  GSHGEOM_Map flip;
  ASSERT_EQ(gshgeom_map_flip_y(2, &flip), GSHGEOM_OK);
  GSHGEOM_Arc m;
  ASSERT_EQ(gshgeom_arc_map(&flip, &a, &m), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_arc_circle(&m, &c, &r, &ccw), GSHGEOM_OK);
  EXPECT_EQ(ccw, 0);
  EXPECT_EQ(m.sweep, a.sweep);   // the parameters do not change
  EXPECT_EQ(m.v.y, -5.0);        // cross(u, v) is now negative
  // a negative sweep on a positively wound ellipse is the other way too
  GSHGEOM_Path * q = make();
  add_form(q, P(0, 0), P(5, 0), P(0, 5), GSHGEOM_ANGLE_DEGREES, 90, -90, GSHGEOM_ARC_CIRCULAR);
  GSHGEOM_Arc b = arc_at(q, 1);
  ASSERT_EQ(gshgeom_arc_circle(&b, &c, &r, &ccw), GSHGEOM_OK);
  EXPECT_EQ(ccw, 0);
  gshgeom_path_destroy(p);
  gshgeom_path_destroy(q);
}

TEST(ArcMap, ASingularMapIsInvalidAndAnUnflaggedArcIsNoCircle) {
  GSHGEOM_Path * p = make();
  add_form(p, P(0, 0), P(5, 0), P(0, 5), GSHGEOM_ANGLE_DEGREES, 0, 90, GSHGEOM_ARC_CIRCULAR);
  GSHGEOM_Arc a = arc_at(p, 1);
  GSHGEOM_Map zero;
  ASSERT_EQ(gshgeom_map_scale(2, 1, 0, 1, &zero), GSHGEOM_OK);
  GSHGEOM_Arc out;
  std::memset(&out, 0x5a, sizeof out);
  GSHGEOM_Arc before = out;
  EXPECT_EQ(gshgeom_arc_map(&zero, &a, &out), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(std::memcmp(&out, &before, sizeof out), 0);  // written only on success
  // the pieces of a singular image are still mapped exactly by the device route's rule: a segment
  GSHGEOM_Curve pieces[GSHGEOM_ARC_MAX_PIECES];
  int n = 0;
  ASSERT_EQ(gshgeom_arc_pieces(&a, pieces, &n), GSHGEOM_OK);
  GSHGEOM_Curve mapped;
  EXPECT_EQ(gshgeom_map_apply_curve(&zero, &pieces[0], &mapped), GSHGEOM_OK);
  // NULLs, a 3D map
  GSHGEOM_Map id = identity();
  EXPECT_EQ(gshgeom_arc_map(nullptr, &a, &out), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_arc_map(&id, nullptr, &out), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_arc_map(&id, &a, nullptr), GSHGEOM_ERR_INVALID);
  GSHGEOM_Map id3;
  ASSERT_EQ(gshgeom_map_identity(3, &id3), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_arc_map(&id3, &a, &out), GSHGEOM_ERR_INVALID);
  // an arc without the flag is no circle, however round its numbers
  GSHGEOM_Arc plain_arc = a;
  plain_arc.flags = 0;
  GSHGEOM_Point c;
  double r;
  EXPECT_EQ(gshgeom_arc_circle(&plain_arc, &c, &r, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_arc_circle(&a, &c, &r, nullptr), GSHGEOM_OK);
  // a coordinate past the limit after the map
  GSHGEOM_Map big;
  ASSERT_EQ(gshgeom_map_scale(2, GSHGEOM_MAP_LIMIT, GSHGEOM_MAP_LIMIT, 1, &big), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_arc_map(&big, &a, &out), GSHGEOM_ERR_LIMIT);
  gshgeom_path_destroy(p);
}

// ---------------------------------------------------------------------------
// AR12 measurement
// ---------------------------------------------------------------------------

TEST(ArcMeasure, EveryVectorsLengthHoldsWithinAlongsTolerance) {
  for (size_t i = 0; i < arc_vectors::length_count; i++) {
    const arc_vectors::Length & l = arc_vectors::lengths[i];
    GSHGEOM_Path * p = make();
    add_form(p, P(0, 0), P(l.ux, l.uy), P(l.vx, l.vy), l.unit, l.t0, l.sweep);
    GSHGEOM_Measured * m = nullptr;
    ASSERT_EQ(gshgeom_measured_create(nullptr, nullptr, p, 0, &m), GSHGEOM_OK) << i;
    // along.h: the estimated error of a segment is at most tol * S, S the length of its control polygon
    GSHGEOM_Path * low = lowered(p);
    double polygon = 0;
    for (size_t k = 1; k < gshgeom_path_element_count(low); k++) {
      GSHGEOM_Element el;
      ASSERT_EQ(gshgeom_path_element(low, k, &el), GSHGEOM_OK);
      for (int j = 0; j < 2; j++) {
        polygon += std::hypot(el.curve.p[j + 1].x - el.curve.p[j].x, el.curve.p[j + 1].y - el.curve.p[j].y);
      }
    }
    EXPECT_NEAR(gshgeom_measured_length(m), l.length, GSHGEOM_MEASURE_TOLERANCE_DEFAULT * polygon) << i;
    gshgeom_measured_destroy(m);
    gshgeom_path_destroy(low);
    gshgeom_path_destroy(p);
  }
}

TEST(ArcMeasure, ACircularArcIsRadiusTimesSweepAndPiecesAreSegments) {
  const double radius = 7.0;
  for (double sweep : {30.0, 90.0, 100.0, 200.0, 270.0, 360.0, -150.0}) {
    GSHGEOM_Path * p = make();
    add_form(p, P(1, 2), P(radius, 0), P(0, radius), GSHGEOM_ANGLE_DEGREES, 20, sweep, GSHGEOM_ARC_CIRCULAR);
    GSHGEOM_Measured * m = nullptr;
    ASSERT_EQ(gshgeom_measured_create(nullptr, nullptr, p, 0, &m), GSHGEOM_OK);
    const double want = radius * std::fabs(sweep) * (M_PI / 180.0);
    EXPECT_NEAR(gshgeom_measured_length(m), want, 4 * GSHGEOM_MEASURE_TOLERANCE_DEFAULT * want) << sweep;
    // segment_count counts pieces
    GSHGEOM_Arc a = arc_at(p, 1);
    GSHGEOM_Curve pieces[GSHGEOM_ARC_MAX_PIECES];
    int n = 0;
    ASSERT_EQ(gshgeom_arc_pieces(&a, pieces, &n), GSHGEOM_OK);
    EXPECT_EQ((int)gshgeom_measured_segment_count(m), n) << sweep;
    EXPECT_EQ(gshgeom_measured_subpath_count(m), 1u);
    gshgeom_measured_destroy(m);
    gshgeom_path_destroy(p);
  }
}

TEST(ArcMeasure, TheFrameNamesTheElementAndThePiece) {
  GSHGEOM_Path * p = make();
  ASSERT_EQ(gshgeom_path_move_to(p, P(10, 0)), GSHGEOM_OK);
  add_form(p, P(0, 0), P(10, 0), P(0, 10), GSHGEOM_ANGLE_DEGREES, 0, 360, GSHGEOM_ARC_CIRCULAR);
  ASSERT_EQ(gshgeom_path_line_to(p, P(20, 0)), GSHGEOM_OK);
  GSHGEOM_Measured * m = nullptr;
  ASSERT_EQ(gshgeom_measured_create(nullptr, nullptr, p, 0, &m), GSHGEOM_OK);
  const double circle = 2 * M_PI * 10;
  for (int k = 0; k < 4; k++) {
    GSHGEOM_Measured_Frame f;
    ASSERT_EQ(gshgeom_measured_at(m, circle * (k + 0.5) / 4, &f), GSHGEOM_OK);
    EXPECT_EQ(f.element, 1u) << k;
    EXPECT_EQ((int)f.piece, k);
    EXPECT_NEAR(f.t, 0.5, 1e-9);
  }
  GSHGEOM_Measured_Frame f;
  ASSERT_EQ(gshgeom_measured_at(m, circle + 3, &f), GSHGEOM_OK);
  EXPECT_EQ(f.element, 2u);
  EXPECT_EQ(f.piece, 0u);
  gshgeom_measured_destroy(m);
  gshgeom_path_destroy(p);
}

// ---------------------------------------------------------------------------
// AR13 stroke expansion and dashes
// ---------------------------------------------------------------------------

namespace {

// a rounded rectangle: lines and quarter arcs (every control point exact)
GSHGEOM_Path * rounded_rect(bool closed = true) {
  GSHGEOM_Path * p = make();
  EXPECT_EQ(gshgeom_path_move_to(p, P(10, 0)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_line_to(p, P(90, 0)), GSHGEOM_OK);
  add_form(p, P(90, 10), P(10, 0), P(0, 10), GSHGEOM_ANGLE_DEGREES, -90, 90, GSHGEOM_ARC_CIRCULAR);
  EXPECT_EQ(gshgeom_path_line_to(p, P(100, 90)), GSHGEOM_OK);
  add_form(p, P(90, 90), P(10, 0), P(0, 10), GSHGEOM_ANGLE_DEGREES, 0, 90, GSHGEOM_ARC_CIRCULAR);
  EXPECT_EQ(gshgeom_path_line_to(p, P(10, 100)), GSHGEOM_OK);
  add_form(p, P(10, 90), P(10, 0), P(0, 10), GSHGEOM_ANGLE_DEGREES, 90, 90, GSHGEOM_ARC_CIRCULAR);
  EXPECT_EQ(gshgeom_path_line_to(p, P(0, 10)), GSHGEOM_OK);
  add_form(p, P(10, 10), P(10, 0), P(0, 10), GSHGEOM_ANGLE_DEGREES, 180, 90, GSHGEOM_ARC_CIRCULAR);
  if (closed) {
    EXPECT_EQ(gshgeom_path_close(p), GSHGEOM_OK);
  }
  return p;
}

}  // namespace

TEST(ArcStroke, ARoundedRectangleStrokesAsItsConicsForEveryJoinAndCapSolidAndDashed) {
  GSHGEOM_Path * p = rounded_rect();
  GSHGEOM_Path * l = lowered(p);
  const double dashes[2] = {23.0, 9.0};
  GSHGEOM_Map maps[3] = {identity(), matrix(1.5, 0.3, -0.2, 1.1, 4, 5), matrix(0.5, 0, 0, 2)};
  for (int mi = 0; mi < 3; mi++) {
    for (uint32_t space : {GSHGEOM_STROKE_WIDTH_LOCAL, GSHGEOM_STROKE_WIDTH_DEVICE}) {
      for (int join : joins) {
        for (int cap : caps) {
          for (int dashed = 0; dashed < 2; dashed++) {
            GSHGEOM_Stroke s = plain();
            s.join = (uint32_t)join;
            s.cap = (uint32_t)cap;
            s.width_space = space;
            if (dashed) {
              s.dashes = dashes;
              s.dash_count = 2;
              s.dash_offset = 5;
            }
            EXPECT_EQ(stroke_hash(p, s, maps[mi]), stroke_hash(l, s, maps[mi]))
                << mi << " " << space << " " << join << " " << cap << " " << dashed;
          }
        }
      }
    }
  }
  gshgeom_path_destroy(p);
  gshgeom_path_destroy(l);
}

TEST(ArcStroke, AFullCircleOfFourPiecesJoinsNothingAtItsBoundaries) {
  // an open circle: caps at the start only (both ends there); every join type gives the same outline, and the conic path
  GSHGEOM_Path * p = make();
  add_form(p, P(50, 50), P(20, 0), P(0, 20), GSHGEOM_ANGLE_DEGREES, 0, 360, GSHGEOM_ARC_CIRCULAR);
  GSHGEOM_Path * l = lowered(p);
  const double dashes[2] = {17.0, 6.0};
  for (int dashed = 0; dashed < 2; dashed++) {
    GSHGEOM_Stroke s = plain();
    if (dashed) {
      s.dashes = dashes;
      s.dash_count = 2;
    }
    GSHGEOM_Map m = identity();
    const uint64_t first = stroke_hash(p, s, m);
    for (int join : joins) {
      s.join = (uint32_t)join;
      EXPECT_EQ(stroke_hash(p, s, m), first) << join;
      EXPECT_EQ(stroke_hash(l, s, m), first) << join;  // tangent-continuous to the bit: the conics join nothing either
    }
  }
  gshgeom_path_destroy(p);
  gshgeom_path_destroy(l);
}

TEST(ArcStroke, NoJoinIsMadeAtAPieceBoundaryInsideOneArc) {
  // 100 degrees of a 3:2 ellipse is two pieces whose tangents at the shared point are continuous but not to the bit
  GSHGEOM_Path * p = make();
  add_form(p, P(0, 0), P(45, 12), P(-8, 30), GSHGEOM_ANGLE_DEGREES, 17, 100);
  GSHGEOM_Arc a = arc_at(p, 1);
  GSHGEOM_Curve pieces[GSHGEOM_ARC_MAX_PIECES];
  int n = 0;
  ASSERT_EQ(gshgeom_arc_pieces(&a, pieces, &n), GSHGEOM_OK);
  ASSERT_EQ(n, 2);
  GSHGEOM_Path * l = lowered(p);
  const double dashes[2] = {31.0, 5.0};
  GSHGEOM_Map m = identity();
  bool conics_differ = false;
  for (int dashed = 0; dashed < 2; dashed++) {
    GSHGEOM_Stroke s = plain();
    s.width = 9;
    if (dashed) {
      s.dashes = dashes;
      s.dash_count = 2;
    }
    s.join = GSHGEOM_STROKE_JOIN_BEVEL;
    const uint64_t first = stroke_hash(p, s, m);
    for (int join : joins) {
      s.join = (uint32_t)join;
      EXPECT_EQ(stroke_hash(p, s, m), first) << join << " " << dashed;
      if (stroke_hash(l, s, m) != first) {
        conics_differ = true;  // the path of conics does make the join, for some join type
      }
    }
  }
  EXPECT_TRUE(conics_differ) << "the control: the same pieces written as conics are joined where the arc's are not";
  gshgeom_path_destroy(p);
  gshgeom_path_destroy(l);
}

TEST(ArcStroke, AJoinIsMadeAtAnArcsEndsFromItsEndDirections) {
  // a line then a 100 degree arc then a line: the same corners as the conic path (one or two pieces)
  GSHGEOM_Path* p = make();
  ASSERT_EQ(gshgeom_path_move_to(p, P(-20, 40)), GSHGEOM_OK);
  GSHGEOM_Point c = P(0, 0);
  GSHGEOM_Arc_Form f = {c, P(30, 0), P(0, 30), GSHGEOM_ANGLE_DEGREES, 0, 100, GSHGEOM_ARC_CIRCULAR};
  GSHGEOM_Point start;
  ASSERT_EQ(gshgeom_arc_point(&f, 0, &start), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(p, start), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_arc_conjugate_to(p, &f, nullptr), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(p, P(-40, 10)), GSHGEOM_OK);
  GSHGEOM_Path * l = lowered(p);
  GSHGEOM_Map m = identity();
  GSHGEOM_Stroke s = plain();
  s.join = GSHGEOM_STROKE_JOIN_BEVEL;
  // the joins at the arc's two ends are the conic path's; the conic path adds one more, at the piece boundary, a bevel of a
  // turn of a few ulps at most: the outlines differ by that sliver only, so compare the covered region at its corners
  GSHGEOM_Stroke_Outline * oa = nullptr;
  GSHGEOM_Stroke_Outline * ol = nullptr;
  ASSERT_EQ(gshgeom_stroke_expand(p, &s, &m, nullptr, nullptr, &oa), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_stroke_expand(l, &s, &m, nullptr, nullptr, &ol), GSHGEOM_OK);
  stroke_oracle::Winding wa, wl;
  wa.build(stroke_oracle::polygons_of(gshgeom_stroke_outline_body(oa)));
  wl.build(stroke_oracle::polygons_of(gshgeom_stroke_outline_body(ol)));
  // the arc's start corner: the line arrives along (1, -0.9), the arc leaves along +y; its outer side is covered by a bevel
  int checked = 0;
  for (double x = -30; x <= 40; x += 0.37) {
    for (double y = -10; y <= 50; y += 0.41) {
      const bool a = wa.at({x, y}) != 0;
      const bool b = wl.at({x, y}) != 0;
      if (a != b) {
        // the only disagreement allowed is the piece-boundary sliver, which lies within 0.1 pixel of the centerline's circle
        const double r = std::hypot(x, y);
        ASSERT_LT(std::fabs(r - 30.0) - 2.0, 0.1) << x << " " << y;
        EXPECT_LT(std::fabs(std::atan2(y, x) * 180 / M_PI - 50.0), 3.0) << x << " " << y;
      }
      checked++;
    }
  }
  EXPECT_GT(checked, 1000);
  gshgeom_stroke_outline_destroy(oa);
  gshgeom_stroke_outline_destroy(ol);
  gshgeom_path_destroy(p);
  gshgeom_path_destroy(l);
}

TEST(ArcStroke, ARegionOfAnEllipseArcIsTheBandAboutIt) {
  // the covered region under SG17's tolerance: the points within the half width of the arc and nowhere else
  GSHGEOM_Path * p = make();
  const GSHGEOM_Point c = P(60, 50), u = P(30, 0), v = P(0, 20);
  add_form(p, c, u, v, GSHGEOM_ANGLE_DEGREES, 20, 250);
  GSHGEOM_Stroke s = plain();
  s.width = 6;
  s.join = GSHGEOM_STROKE_JOIN_ROUND;
  s.cap = GSHGEOM_STROKE_CAP_BUTT;
  GSHGEOM_Map m = identity();
  GSHGEOM_Stroke_Outline * o = nullptr;
  ASSERT_EQ(gshgeom_stroke_expand(p, &s, &m, nullptr, nullptr, &o), GSHGEOM_OK);
  stroke_oracle::Winding w;
  w.build(stroke_oracle::polygons_of(gshgeom_stroke_outline_body(o)));
  const double h = 3.0;
  int inside = 0, outside = 0;
  for (int k = 1; k < 100; k++) {
    const double t = (20 + 250.0 * k / 100) * M_PI / 180;
    const double px = c.x + u.x * std::cos(t), py = c.y + v.y * std::sin(t);
    double tx = -u.x * std::sin(t), ty = v.y * std::cos(t);
    const double len = std::hypot(tx, ty);
    tx /= len;
    ty /= len;
    for (double d : {-h - 0.5, -h + 0.5, 0.0, h - 0.5, h + 0.5}) {
      const bool want = std::fabs(d) < h;
      const bool got = w.at({px + d * ty, py - d * tx}) != 0;
      EXPECT_EQ(got, want) << k << " " << d;
      (want ? inside : outside)++;
    }
  }
  EXPECT_GT(inside, 100);
  EXPECT_GT(outside, 100);
  gshgeom_stroke_outline_destroy(o);
  gshgeom_path_destroy(p);
}

TEST(ArcStroke, ADashCutInsideAnArcAndAcrossAPieceBoundaryIsStrokedAsTheCutConics) {
  // a half circle of two quarter pieces, dashed with a dash that straddles the boundary at its middle
  GSHGEOM_Path * p = make();
  add_form(p, P(0, 0), P(30, 0), P(0, 30), GSHGEOM_ANGLE_DEGREES, 0, 180, GSHGEOM_ARC_CIRCULAR);
  GSHGEOM_Path * l = lowered(p);
  const double quarter = 30 * M_PI / 2;
  for (double len : {10.0, 25.0, 60.0, 300.0}) {
    const double dashes[2] = {len, 7.0};
    for (double offset : {0.0, quarter - len / 2, quarter - 1, 3.0}) {
      for (uint32_t space : {GSHGEOM_STROKE_WIDTH_LOCAL, GSHGEOM_STROKE_WIDTH_DEVICE}) {
        GSHGEOM_Stroke s = plain();
        s.dashes = dashes;
        s.dash_count = 2;
        s.dash_offset = -offset;  // the pattern starts late, so a dash straddles the middle
        s.cap = GSHGEOM_STROKE_CAP_ROUND;
        s.width_space = space;
        GSHGEOM_Map m = matrix(1.3, 0.2, 0.1, 0.8, 10, 10);
        EXPECT_EQ(stroke_hash(p, s, m), stroke_hash(l, s, m)) << len << " " << offset << " " << space;
        GSHGEOM_Map id = identity();
        EXPECT_EQ(stroke_hash(p, s, id), stroke_hash(l, s, id)) << len << " " << offset << " " << space;
      }
    }
  }
  gshgeom_path_destroy(p);
  gshgeom_path_destroy(l);
}

// max_dash_steps (stroke.h) bounds the steps of a dashed stroke.  An arc is measured as its pieces (AR12), so a dashed arc and the same pieces
// written as conics are refused at the same cap: the least cap that passes is equal for the two, one less is refused for both, and for one
// open subpath it is the walk's own count over the measured length (the dash positions of that length, with the same bound).
TEST(ArcStroke, TheDashStepsOfAnArcAreThoseOfItsConics) {
  GSHGEOM_Path * p = make();
  ASSERT_EQ(gshgeom_path_move_to(p, P(300, 0)), GSHGEOM_OK);
  add_form(p, P(0, 0), P(300, 0), P(0, 300), GSHGEOM_ANGLE_DEGREES, 0, 270, GSHGEOM_ARC_CIRCULAR);
  GSHGEOM_Path * l = lowered(p);
  GSHGEOM_Measured * measured = nullptr;
  ASSERT_EQ(gshgeom_measured_create(nullptr, nullptr, p, GSHGEOM_MEASURE_TOLERANCE_DEFAULT, &measured), GSHGEOM_OK);
  const double length = gshgeom_measured_length(measured);
  gshgeom_measured_destroy(measured);
  EXPECT_NEAR(length, 300.0 * 1.5 * M_PI, 1e-3 * length);
  const double dashes[2] = {1.5, 1.5};
  GSHGEOM_Stroke s = plain();
  s.dashes = dashes;
  s.dash_count = 2;
  s.width = 1.0;
  const GSHGEOM_Map m = identity();
  auto passes = [&](const GSHGEOM_Path * path, uint64_t cap) {
    GSHGEOM_Limits lim;
    gshgeom_limits_default(&lim);
    lim.max_dash_steps = cap;
    GSHGEOM_Stroke_Outline * o = nullptr;
    const GSHGEOM_Result r = gshgeom_stroke_expand(path, &s, &m, nullptr, &lim, &o);
    gshgeom_stroke_outline_destroy(o);
    EXPECT_TRUE(r == GSHGEOM_OK || r == GSHGEOM_ERR_LIMIT) << r;
    return r == GSHGEOM_OK;
  };
  auto least = [&](const GSHGEOM_Path * path) {
    uint64_t lo = 1, hi = 1u << 22;
    EXPECT_TRUE(passes(path, hi));
    while (lo < hi) {
      const uint64_t mid = lo + (hi - lo) / 2;
      if (passes(path, mid)) {
        hi = mid;
      } else {
        lo = mid + 1;
      }
    }
    return lo;
  };
  const uint64_t arc_least = least(p);
  const uint64_t conic_least = least(l);
  EXPECT_EQ(arc_least, conic_least) << "an arc and its conics make the same steps";
  EXPECT_TRUE(passes(p, arc_least));
  EXPECT_FALSE(passes(p, arc_least - 1)) << "one past the bound is refused";
  EXPECT_FALSE(passes(l, conic_least - 1));
  auto nothing = [](void *, double, double) { return GSHGEOM_OK; };
  uint64_t a = 1, b = 1u << 22;
  while (a < b) {
    const uint64_t mid = a + (b - a) / 2;
    if (gshgeom_stroke_dash_each(dashes, 2, 0.0, length, (size_t)mid, nothing, nullptr) == GSHGEOM_OK) {
      b = mid;
    } else {
      a = mid + 1;
    }
  }
  EXPECT_EQ(arc_least, a) << "the walk's own count over the measured length";
  gshgeom_path_destroy(p);
  gshgeom_path_destroy(l);
}

TEST(ArcStroke, ZeroWidthAFullTurnAndADegenerateArcFollowTheStrokeRules) {
  GSHGEOM_Path * p = make();
  add_form(p, P(5, 5), P(10, 0), P(0, 10), GSHGEOM_ANGLE_DEGREES, 0, 360, GSHGEOM_ARC_CIRCULAR);
  GSHGEOM_Map m = identity();
  GSHGEOM_Stroke s = plain();
  s.width = 0;
  GSHGEOM_Stroke_Outline * o = nullptr;
  ASSERT_EQ(gshgeom_stroke_expand(p, &s, &m, nullptr, nullptr, &o), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_stroke_outline_status(o), GSHGEOM_STROKE_EMPTY);
  gshgeom_stroke_outline_destroy(o);
  // an arc with no extent is a degenerate subpath, as a point is
  GSHGEOM_Path * q = make();
  ASSERT_EQ(gshgeom_path_move_to(q, P(1e5, 1e5)), GSHGEOM_OK);
  GSHGEOM_Arc_Form f = {P(1e5, 1e5), P(1e-12, 0), P(0, 1e-12), GSHGEOM_ANGLE_DEGREES, 0, 90, 0};
  ASSERT_EQ(gshgeom_path_arc_conjugate_to(q, &f, nullptr), GSHGEOM_OK);
  GSHGEOM_Path * ql = lowered(q);
  GSHGEOM_Stroke r = plain();
  r.cap = GSHGEOM_STROKE_CAP_ROUND;
  EXPECT_EQ(stroke_hash(q, r, m), stroke_hash(ql, r, m));
  gshgeom_path_destroy(p);
  gshgeom_path_destroy(q);
  gshgeom_path_destroy(ql);
}

TEST(ArcStroke, TheCurveOutlineAndItsCubicLoweringTakeArcs) {
  GSHGEOM_Path * p = rounded_rect();
  GSHGEOM_Path * l = lowered(p);
  GSHGEOM_Map m = identity();
  GSHGEOM_Stroke s = plain();
  s.join = GSHGEOM_STROKE_JOIN_ROUND;
  for (GSHGEOM_Stroke_Curves_Form form : {GSHGEOM_STROKE_CURVES_RATIONAL, GSHGEOM_STROKE_CURVES_CUBIC}) {
    GSHGEOM_Stroke_Curves * a = nullptr;
    GSHGEOM_Stroke_Curves * b = nullptr;
    ASSERT_EQ(gshgeom_stroke_curves_expand(p, &s, &m, GSHGEOM_STROKE_CURVES_TOLERANCE, form, nullptr, nullptr, &a), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_stroke_curves_expand(l, &s, &m, GSHGEOM_STROKE_CURVES_TOLERANCE, form, nullptr, nullptr, &b), GSHGEOM_OK);
    EXPECT_EQ(gshgeom_stroke_curves_hash(a), gshgeom_stroke_curves_hash(b));
    gshgeom_stroke_curves_destroy(a);
    gshgeom_stroke_curves_destroy(b);
  }
  // the cubic lowering of a path of arcs is the lowering of its pieces
  GSHGEOM_Path * ca = nullptr;
  GSHGEOM_Path * cl = nullptr;
  ASSERT_EQ(gshgeom_stroke_curves_lower_cubic(p, 0.01, nullptr, nullptr, &ca), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_stroke_curves_lower_cubic(l, 0.01, nullptr, nullptr, &cl), GSHGEOM_OK);
  EXPECT_TRUE(gshgeom_path_equal(ca, cl));
  gshgeom_path_destroy(ca);
  gshgeom_path_destroy(cl);
  gshgeom_path_destroy(p);
  gshgeom_path_destroy(l);
}

TEST(ArcStroke, AlignmentClipsToTheArcsClosedRegion) {
  // a closed subpath of an arc and a line (a half disc), stroked inside and outside: the clip is built from the pieces
  GSHGEOM_Path * p = make();
  add_form(p, P(0, 0), P(20, 0), P(0, 20), GSHGEOM_ANGLE_DEGREES, 0, 180, GSHGEOM_ARC_CIRCULAR);
  ASSERT_EQ(gshgeom_path_close(p), GSHGEOM_OK);
  GSHGEOM_Path * l = lowered(p);
  GSHGEOM_Map m = identity();
  for (uint32_t al : {GSHGEOM_STROKE_ALIGN_INSIDE, GSHGEOM_STROKE_ALIGN_OUTSIDE}) {
    GSHGEOM_Stroke s = plain();
    s.alignment = al;
    EXPECT_EQ(stroke_hash(p, s, m), stroke_hash(l, s, m)) << al;
  }
  gshgeom_path_destroy(p);
  gshgeom_path_destroy(l);
}

// ---------------------------------------------------------------------------
// AR14 path vertices
// ---------------------------------------------------------------------------

namespace {

GSHGEOM_Result collect(const GSHGEOM_Vertex * v, void * ctx) {
  static_cast<std::vector<GSHGEOM_Vertex> *>(ctx)->push_back(*v);
  return GSHGEOM_OK;
}

std::vector<GSHGEOM_Vertex> vertices_of(const GSHGEOM_Path * p) {
  std::vector<GSHGEOM_Vertex> out;
  EXPECT_EQ(gshgeom_path_vertices(p, collect, &out), GSHGEOM_OK);
  size_t count = 0;
  EXPECT_EQ(gshgeom_path_vertex_count(p, &count), GSHGEOM_OK);
  EXPECT_EQ(count, out.size());
  EXPECT_EQ(count, gshgeom_path_element_count(p));
  return out;
}

bool same(GSHGEOM_Point a, GSHGEOM_Point b) {
  return bits(a.x) == bits(b.x) && bits(a.y) == bits(b.y);
}

GSHGEOM_SK_Seg sk_of(const GSHGEOM_Curve & c) {
  GSHGEOM_SK_Seg k;
  std::memset(&k, 0, sizeof k);
  k.n = gshgeom_curve_point_count(c.kind) - 1;
  for (int i = 0; i <= k.n; i++) {
    k.p[i] = gshgeom_sk_v(c.p[i].x, c.p[i].y);
  }
  return k;
}

}  // namespace

TEST(ArcVertices, ARoundedRectanglesCornersHaveBitEqualIncomingAndOutgoingDirections) {
  GSHGEOM_Path * p = rounded_rect();
  std::vector<GSHGEOM_Vertex> vs = vertices_of(p);
  ASSERT_EQ(vs.size(), 10u);  // move, 4 lines, 4 arcs, close
  const GSHGEOM_Point east = P(1, 0), south = P(0, 1), west = P(-1, 0), north = P(0, -1);
  // the line to the first corner: it ends going east and the arc leaves going east
  EXPECT_TRUE(same(vs[1].in, east));
  EXPECT_TRUE(same(vs[1].out, east));
  // the arc's end (100, 10): in and out are the same vector, south, to the bit, and the marker direction is exactly that
  EXPECT_EQ(vs[2].element, 2u);
  EXPECT_EQ(vs[2].point.x, 100.0);
  EXPECT_EQ(vs[2].point.y, 10.0);
  EXPECT_TRUE(same(vs[2].in, south));
  EXPECT_TRUE(same(vs[2].out, south));
  EXPECT_TRUE(same(vs[2].direction, south));
  EXPECT_TRUE(same(vs[4].in, west));
  EXPECT_TRUE(same(vs[4].out, west));
  EXPECT_TRUE(same(vs[4].direction, west));
  EXPECT_TRUE(same(vs[6].in, north));
  EXPECT_TRUE(same(vs[6].out, north));
  EXPECT_TRUE(same(vs[8].in, east));
  EXPECT_TRUE(same(vs[8].out, east));
  // an arc is one vertex: its piece boundaries add none
  for (size_t i = 0; i < vs.size(); i++) {
    EXPECT_EQ(vs[i].element, i);
  }
  gshgeom_path_destroy(p);
}

TEST(ArcVertices, TheDirectionsAreTheOnesTheStrokeCodeReadsOnTheSamePieces) {
  // arcs of 1 to 4 pieces, ellipses, both senses, a mirror: the unit directions are gshgeom_sk's, to the bit
  const double sweeps[] = {40, 91, 200, 330, 360, -75, -250};
  for (double sweep : sweeps) {
    for (int mirror = 0; mirror < 2; mirror++) {
      GSHGEOM_Path * p = make();
      GSHGEOM_Point c = P(3, -4), u = P(30, 6), v = mirror ? P(5, -20) : P(-5, 20);
      add_form(p, c, u, v, GSHGEOM_ANGLE_DEGREES, 23, sweep);
      ASSERT_EQ(gshgeom_path_line_to(p, P(-100, 100)), GSHGEOM_OK);
      GSHGEOM_Arc a = arc_at(p, 1);
      GSHGEOM_Curve pieces[GSHGEOM_ARC_MAX_PIECES];
      int n = 0;
      ASSERT_EQ(gshgeom_arc_pieces(&a, pieces, &n), GSHGEOM_OK);
      GSHGEOM_SK_Seg first = sk_of(pieces[0]), last = sk_of(pieces[n - 1]);
      GSHGEOM_SK_Vec ds, de;
      ASSERT_GE(gshgeom_sk_start_dir(&first, &ds), 0);
      ASSERT_GE(gshgeom_sk_end_dir(&last, &de), 0);
      std::vector<GSHGEOM_Vertex> vs = vertices_of(p);
      ASSERT_EQ(vs.size(), 3u);
      const GSHGEOM_SK_Vec us = gshgeom_sk_unit(ds), ue = gshgeom_sk_unit(de);
      EXPECT_EQ(bits(vs[1].in.x), bits(ue.x)) << sweep;
      EXPECT_EQ(bits(vs[1].in.y), bits(ue.y)) << sweep;
      EXPECT_EQ(bits(vs[0].out.x), bits(us.x)) << sweep;  // the move looks forward to the arc's start direction
      EXPECT_EQ(bits(vs[0].out.y), bits(us.y)) << sweep;
      EXPECT_EQ(vs[1].point.x, a.end.x);  // the vertex is the authoritative end
      EXPECT_EQ(vs[1].point.y, a.end.y);
      EXPECT_EQ(bits(vs[1].point.x), bits(a.end.x));
      gshgeom_path_destroy(p);
    }
  }
}

TEST(ArcVertices, ASmoothJoinBetweenArcsAReversalAndACorner) {
  // two quarter arcs of one circle: the second's start direction is the first's end direction, the bisector is exact
  GSHGEOM_Path * p = make();
  add_form(p, P(0, 0), P(10, 0), P(0, 10), GSHGEOM_ANGLE_DEGREES, 0, 90, GSHGEOM_ARC_CIRCULAR);
  add_form(p, P(0, 0), P(10, 0), P(0, 10), GSHGEOM_ANGLE_DEGREES, 90, 90, GSHGEOM_ARC_CIRCULAR);
  std::vector<GSHGEOM_Vertex> vs = vertices_of(p);
  ASSERT_EQ(vs.size(), 3u);
  EXPECT_TRUE(same(vs[1].in, P(-1, 0)));  // ends at (0, 10) going -x
  EXPECT_TRUE(same(vs[1].out, P(-1, 0)));
  EXPECT_TRUE(same(vs[1].direction, P(-1, 0)));
  gshgeom_path_destroy(p);
  // a reversal: an arc out and the same arc back (a negative sweep from its end): the bisector is the incoming turned
  p = make();
  add_form(p, P(0, 0), P(10, 0), P(0, 10), GSHGEOM_ANGLE_DEGREES, 0, 90, GSHGEOM_ARC_CIRCULAR);
  add_form(p, P(0, 0), P(10, 0), P(0, 10), GSHGEOM_ANGLE_DEGREES, 90, -90, GSHGEOM_ARC_CIRCULAR);
  vs = vertices_of(p);
  ASSERT_EQ(vs.size(), 3u);
  EXPECT_TRUE(same(vs[1].in, P(-1, 0)));
  EXPECT_TRUE(same(vs[1].out, P(1, 0)));
  EXPECT_TRUE(same(vs[1].direction, gshgeom_normal_right(P(-1, 0))));
  gshgeom_path_destroy(p);
  // a corner: a quarter arc then a line going back along its chord
  p = make();
  add_form(p, P(0, 0), P(10, 0), P(0, 10), GSHGEOM_ANGLE_DEGREES, 0, 90, GSHGEOM_ARC_CIRCULAR);
  ASSERT_EQ(gshgeom_path_line_to(p, P(0, 0)), GSHGEOM_OK);
  vs = vertices_of(p);
  ASSERT_EQ(vs.size(), 3u);
  EXPECT_TRUE(same(vs[1].in, P(-1, 0)));
  EXPECT_TRUE(same(vs[1].out, P(0, -1)));
  EXPECT_NEAR(vs[1].direction.x, -std::sqrt(0.5), 1e-15);
  EXPECT_NEAR(vs[1].direction.y, -std::sqrt(0.5), 1e-15);
  gshgeom_path_destroy(p);
}

TEST(ArcVertices, AClosedSubpathThatStartsWithAnArcIsReadRound) {
  GSHGEOM_Path * p = make();
  add_form(p, P(0, 0), P(10, 0), P(0, 10), GSHGEOM_ANGLE_DEGREES, 0, 90, GSHGEOM_ARC_CIRCULAR);  // (10, 0) to (0, 10)
  ASSERT_EQ(gshgeom_path_line_to(p, P(0, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_close(p), GSHGEOM_OK);  // back to (10, 0)
  std::vector<GSHGEOM_Vertex> vs = vertices_of(p);
  ASSERT_EQ(vs.size(), 4u);
  // the start vertex: incoming is the close's direction (+x), outgoing the arc's start direction (+y)
  EXPECT_TRUE(same(vs[0].in, P(1, 0)));
  EXPECT_TRUE(same(vs[0].out, P(0, 1)));
  EXPECT_TRUE(same(vs[3].in, P(1, 0)));
  EXPECT_TRUE(same(vs[3].out, P(0, 1)));
  EXPECT_NE(vs[3].flags & GSHGEOM_VERTEX_CLOSED, 0u);
  gshgeom_path_destroy(p);
}

TEST(ArcVertices, AFullTurnArcHasBothDirectionsAndAZeroExtentArcIsLookedPast) {
  GSHGEOM_Path * p = make();
  add_form(p, P(10, 0), P(10, 0), P(0, 10), GSHGEOM_ANGLE_DEGREES, 0, 360, GSHGEOM_ARC_CIRCULAR);
  std::vector<GSHGEOM_Vertex> vs = vertices_of(p);
  ASSERT_EQ(vs.size(), 2u);
  EXPECT_EQ(vs[1].flags & (GSHGEOM_VERTEX_HAS_IN | GSHGEOM_VERTEX_HAS_OUT), GSHGEOM_VERTEX_HAS_IN | 0u);  // the last vertex has no next
  EXPECT_TRUE(same(vs[1].in, P(0, 1)));
  EXPECT_TRUE(same(vs[0].out, P(0, 1)));
  gshgeom_path_destroy(p);
  // an arc too small to move any double of its position has no extent
  p = make();
  ASSERT_EQ(gshgeom_path_move_to(p, P(1e8 - 100, 1e8)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(p, P(1e8, 1e8)), GSHGEOM_OK);
  GSHGEOM_Arc_Form f = {P(1e8, 1e8), P(1e-9, 0), P(0, 1e-9), GSHGEOM_ANGLE_DEGREES, 0, 90, 0};
  ASSERT_EQ(gshgeom_path_arc_conjugate_to(p, &f, nullptr), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(p, P(1e8, 1e8 + 5)), GSHGEOM_OK);
  vs = vertices_of(p);
  ASSERT_EQ(vs.size(), 4u);
  EXPECT_TRUE(same(vs[2].in, P(1, 0)));   // the line before it
  EXPECT_TRUE(same(vs[2].out, P(0, 1)));  // the line after it
  gshgeom_path_destroy(p);
}

TEST(ArcVertices, AllocationFailuresInStrokeAndDashPathsLeaveNothing) {
  GSHGEOM_Path * p = rounded_rect();
  const double dashes[2] = {23.0, 9.0};
  GSHGEOM_Stroke s = plain();
  s.dashes = dashes;
  s.dash_count = 2;
  s.width_space = GSHGEOM_STROKE_WIDTH_DEVICE;  // the mapped centerline is built with the origins of its elements
  GSHGEOM_Map m = matrix(1.5, 0.3, -0.2, 1.1);
  long total = 0;
  {
    tracker::Tracker t;
    GSHGEOM_Allocator alloc = tracker::allocator_for(&t);
    GSHGEOM_Stroke_Outline * o = nullptr;
    ASSERT_EQ(gshgeom_stroke_expand(p, &s, &m, &alloc, nullptr, &o), GSHGEOM_OK);
    gshgeom_stroke_outline_destroy(o);
    total = t.calls;
    EXPECT_TRUE(t.live.empty());
  }
  ASSERT_GT(total, 0);
  for (long fail = 0; fail < total; fail++) {
    tracker::Tracker t;
    t.fail_from = fail;
    GSHGEOM_Allocator alloc = tracker::allocator_for(&t);
    GSHGEOM_Stroke_Outline * o = nullptr;
    GSHGEOM_Result r = gshgeom_stroke_expand(p, &s, &m, &alloc, nullptr, &o);
    EXPECT_EQ(r, GSHGEOM_ERR_OOM) << fail;
    EXPECT_EQ(o, nullptr);
    EXPECT_TRUE(t.live.empty()) << fail;
  }
  gshgeom_path_destroy(p);
}

// An aligned stroke clips to the region the centerline encloses: for an arc that region is built from its pieces after the map, and
// the curve outline carries it as conics.  Under maps that scale, skew and mirror, and for every join, the stroke and the curve outline
// of the arcs are the stroke and the curve outline of the pieces written as conics.
TEST(ArcStroke, TheClipOfAnAlignedStrokeIsTheMappedPiecesRegionForBothOutlines) {
  GSHGEOM_Path * p = rounded_rect();
  GSHGEOM_Path * l = lowered(p);
  GSHGEOM_Map maps[3] = {identity(), matrix(1.5, 0.3, -0.2, 1.1, 4, 5), matrix(0.5, 0, 0, -2)};
  for (int mi = 0; mi < 3; mi++) {
    for (uint32_t al : {GSHGEOM_STROKE_ALIGN_INSIDE, GSHGEOM_STROKE_ALIGN_OUTSIDE}) {
      for (int join : {GSHGEOM_STROKE_JOIN_MITER, GSHGEOM_STROKE_JOIN_ROUND}) {
        GSHGEOM_Stroke s = plain();
        s.alignment = al;
        s.join = (uint32_t)join;
        EXPECT_EQ(stroke_hash(p, s, maps[mi]), stroke_hash(l, s, maps[mi])) << mi << " " << al << " " << join;
        GSHGEOM_Stroke_Curves * a = nullptr;
        GSHGEOM_Stroke_Curves * b = nullptr;
        ASSERT_EQ(gshgeom_stroke_curves_expand(p, &s, &maps[mi], GSHGEOM_STROKE_CURVES_TOLERANCE, GSHGEOM_STROKE_CURVES_RATIONAL,
                      nullptr, nullptr, &a),
            GSHGEOM_OK);
        ASSERT_EQ(gshgeom_stroke_curves_expand(l, &s, &maps[mi], GSHGEOM_STROKE_CURVES_TOLERANCE, GSHGEOM_STROKE_CURVES_RATIONAL,
                      nullptr, nullptr, &b),
            GSHGEOM_OK);
        EXPECT_EQ(gshgeom_stroke_curves_hash(a), gshgeom_stroke_curves_hash(b)) << mi << " " << al << " " << join;
        gshgeom_stroke_curves_destroy(a);
        gshgeom_stroke_curves_destroy(b);
      }
    }
  }
  // the clip is flattened at the stroke's own tolerance: a finer one changes the region of an arc exactly as it changes the pieces'
  GSHGEOM_Stroke s = plain();
  s.alignment = GSHGEOM_STROKE_ALIGN_INSIDE;
  GSHGEOM_Stroke_Outline * fine_a = nullptr;
  GSHGEOM_Stroke_Outline * fine_b = nullptr;
  ASSERT_EQ(gshgeom_stroke_expand_with_tolerance(p, &s, &maps[1], 2, nullptr, nullptr, &fine_a), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_stroke_expand_with_tolerance(l, &s, &maps[1], 2, nullptr, nullptr, &fine_b), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_stroke_outline_hash(fine_a), gshgeom_stroke_outline_hash(fine_b));
  GSHGEOM_Stroke_Outline * dflt = nullptr;
  ASSERT_EQ(gshgeom_stroke_expand_with_tolerance(p, &s, &maps[1], GSHGEOM_FLATTEN_TOLERANCE, nullptr, nullptr, &dflt), GSHGEOM_OK);
  EXPECT_NE(gshgeom_stroke_outline_hash(fine_a), gshgeom_stroke_outline_hash(dflt));
  gshgeom_stroke_outline_destroy(fine_a);
  gshgeom_stroke_outline_destroy(fine_b);
  gshgeom_stroke_outline_destroy(dflt);
  gshgeom_path_destroy(p);
  gshgeom_path_destroy(l);
}

// The cubic lowering of a path of arcs is the lowering of its pieces, for an arc of several pieces (each lowered by its own weight).
TEST(ArcStroke, TheCubicLoweringOfAMultiPieceArcIsTheLoweringOfItsPieces) {
  for (double sweep : {100.0, 230.0, 330.0, -300.0}) {
    GSHGEOM_Path * p = make();
    add_form(p, P(10, 20), P(30, 5), P(-4, 12), GSHGEOM_ANGLE_DEGREES, 12.5, sweep);
    GSHGEOM_Path * l = lowered(p);
    GSHGEOM_Path * ca = nullptr;
    GSHGEOM_Path * cl = nullptr;
    ASSERT_EQ(gshgeom_stroke_curves_lower_cubic(p, 0.001, nullptr, nullptr, &ca), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_stroke_curves_lower_cubic(l, 0.001, nullptr, nullptr, &cl), GSHGEOM_OK);
    EXPECT_TRUE(gshgeom_path_equal(ca, cl)) << sweep;
    EXPECT_GE(gshgeom_path_element_count(ca), 3u);
    gshgeom_path_destroy(ca);
    gshgeom_path_destroy(cl);
    gshgeom_path_destroy(p);
    gshgeom_path_destroy(l);
  }
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
