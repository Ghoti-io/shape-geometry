/**
 * @file
 *
 * The walkers of a path that holds arcs (copy, equality, readback, bounds, dump, measure, stroke, dash, vertices, the lowering
 * to cubics), at the inputs where each one's error or field is the only thing that tells it from a variant: an arc record that
 * no constructor stores (stored raw, with a sweep that is not a number), a stream that fails at one write and no other, an
 * allocator that fails at each call in turn.
 *
 * Copyright 2026 by Corey Pennycuff
 */

#include "../test_helpers.h"

#include "../alloc_tracker.h"

extern "C" {
#include "../../src/path/path_internal.h"
}

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>

namespace {

GSHGEOM_Point P(double x, double y) {
  return GSHGEOM_Point{x, y, 0.0};
}

GSHGEOM_Path * make() {
  GSHGEOM_Path * p = nullptr;
  EXPECT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &p), GSHGEOM_OK);
  return p;
}

// a path with a move and a raw arc whose sweep is `sweep`
GSHGEOM_Path * raw_arc_path(double sweep, unsigned unit = GSHGEOM_ANGLE_DEGREES) {
  GSHGEOM_Path * p = make();
  EXPECT_EQ(gshgeom_path_move_to(p, P(1, 0)), GSHGEOM_OK);
  GSHGEOM_Point pts[4] = {P(0, 0), P(1, 0), P(0, 1), P(0, 1)};
  const double w[2] = {0.0, sweep};
  EXPECT_EQ(gshgeom_path_append_raw(p, GSHGEOM_VERB_ARC, unit, GSHGEOM_ARC_CIRCULAR, pts, w), GSHGEOM_OK);
  return p;
}

GSHGEOM_Path * good_arc_path() {
  GSHGEOM_Path * p = make();
  EXPECT_EQ(gshgeom_path_move_to(p, P(1, 0)), GSHGEOM_OK);
  GSHGEOM_Arc_Form f = {P(0, 0), P(1, 0), P(0, 1), GSHGEOM_ANGLE_DEGREES, 0, 90, GSHGEOM_ARC_CIRCULAR};
  EXPECT_EQ(gshgeom_path_arc_conjugate_to(p, &f, nullptr), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_line_to(p, P(0, 0)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_close(p), GSHGEOM_OK);
  return p;
}

GSHGEOM_Result count_vertex(const GSHGEOM_Vertex *, void * ctx) {
  ++*static_cast<int *>(ctx);
  return GSHGEOM_OK;
}

}  // namespace

namespace {

// a path with a move and one record stored raw that no constructor stores: 0 and 1, an arc whose sweep is not a number and one whose
// unit is past the last; 2, a line whose end is not a number; 3, a cubic whose end is not a number
GSHGEOM_Path * corrupt_path(int which) {
  const double nan = std::nan("");
  GSHGEOM_Path * p = make();
  EXPECT_EQ(gshgeom_path_move_to(p, P(1, 0)), GSHGEOM_OK);
  if (which < 2) {
    GSHGEOM_Point pts[4] = {P(0, 0), P(1, 0), P(0, 1), P(0, 1)};
    const double w[2] = {0.0, which == 0 ? nan : 90.0};
    EXPECT_EQ(gshgeom_path_append_raw(p, GSHGEOM_VERB_ARC, which == 0 ? (unsigned)GSHGEOM_ANGLE_DEGREES : (unsigned)GSHGEOM_ANGLE_UNIT_COUNT,
                  GSHGEOM_ARC_CIRCULAR, pts, w), GSHGEOM_OK);
  } else if (which == 2) {
    const GSHGEOM_Point q[1] = {P(nan, 1)};
    EXPECT_EQ(gshgeom_path_append_raw(p, GSHGEOM_VERB_LINE, 0, 0, q, nullptr), GSHGEOM_OK);
  } else {
    const GSHGEOM_Point q[3] = {P(2, 1), P(3, 1), P(4, nan)};
    EXPECT_EQ(gshgeom_path_append_raw(p, GSHGEOM_VERB_CUBIC, 0, 0, q, nullptr), GSHGEOM_OK);
  }
  return p;
}

}  // namespace

TEST(ArcWalkerEdges, EveryWalkerReportsARecordThatCannotBeRead) {
  // The code each walker reports for each corrupt path (the walkers read a stored element and every failure of that read is the
  // walker's own: an arc that cannot be decomposed, a curve that cannot be bounded or measured, an element that cannot be read).  A
  // walker that lost one of these errors would report success, or a different code.
  enum { BOUNDS, MEASURE, LOWER, STROKE, DASH, DASH_SCALED, DASH_DEVICE, VERTEX_COUNT, VERTICES, DEVICE, COLUMNS };
  const GSHGEOM_Result I = GSHGEOM_ERR_INTERNAL, V = GSHGEOM_ERR_INVALID, L = GSHGEOM_ERR_LIMIT, K = GSHGEOM_OK;
  const GSHGEOM_Result want[4][COLUMNS] = {
      // bounds measure lower stroke dash dash-scaled dash-device vertex-count vertices device
      {I, V, V, I, I, I, I, K, I, V},
      {I, V, V, I, I, I, I, K, I, V},
      {I, K, V, L, L, L, V, K, K, V},
      {I, L, V, L, L, L, V, K, K, V},
  };
  for (int which = 0; which < 4; which++) {
    GSHGEOM_Path * p = corrupt_path(which);
    GSHGEOM_Result got[COLUMNS];
    GSHGEOM_Map map, scaled, map0;
    ASSERT_EQ(gshgeom_map_identity(2, &map), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_map_identity(2, &map0), GSHGEOM_OK);
    GSHGEOM_Bounds b;
    got[BOUNDS] = gshgeom_path_bounds(p, &b);
    GSHGEOM_Measured * m = nullptr;
    got[MEASURE] = gshgeom_measured_create(nullptr, nullptr, p, 0.0, &m);  // the default tolerance: one outside the range is a limit of its own
    if (got[MEASURE] != GSHGEOM_OK) {
      EXPECT_EQ(m, nullptr);
    }
    gshgeom_measured_destroy(m);
    GSHGEOM_Device_Path * dev = nullptr;
    got[DEVICE] = gshgeom_device_path_from_model(p, &map0, nullptr, nullptr, &dev);
    gshgeom_device_path_destroy(dev);
    GSHGEOM_Path * lowered = nullptr;
    got[LOWER] = gshgeom_stroke_curves_lower_cubic(p, 0.25, nullptr, nullptr, &lowered);
    EXPECT_EQ(lowered, nullptr);
    GSHGEOM_Stroke s;
    gshgeom_stroke_init(&s);
    s.width = 2.0;
    ASSERT_EQ(gshgeom_map_scale(2, 2.0, 3.0, 1.0, &scaled), GSHGEOM_OK);
    GSHGEOM_Stroke_Outline * o = nullptr;
    got[STROKE] = gshgeom_stroke_expand(p, &s, &map, nullptr, nullptr, &o);
    EXPECT_EQ(o, nullptr);
    const double dashes[2] = {0.5, 0.5};
    s.dashes = dashes;
    s.dash_count = 2;
    got[DASH] = gshgeom_stroke_expand(p, &s, &map, nullptr, nullptr, &o);
    got[DASH_SCALED] = gshgeom_stroke_expand(p, &s, &scaled, nullptr, nullptr, &o);
    s.width_space = GSHGEOM_STROKE_WIDTH_DEVICE;
    got[DASH_DEVICE] = gshgeom_stroke_expand(p, &s, &scaled, nullptr, nullptr, &o);
    EXPECT_EQ(o, nullptr);
    size_t count = 0;
    got[VERTEX_COUNT] = gshgeom_path_vertex_count(p, &count);
    int n = 0;
    got[VERTICES] = gshgeom_path_vertices(p, count_vertex, &n);
    for (int k = 0; k < COLUMNS; k++) {
      EXPECT_EQ(got[k], want[which][k]) << "corrupt path " << which << ", walker " << k;
    }
    gshgeom_path_destroy(p);
  }
}

TEST(ArcWalkerEdges, ACopyCarriesTheCurrentPointAndTheStartOfTheSubpath) {
  // a subpath that ends in a move: the copy must close back to that move, and continue from it
  GSHGEOM_Path * a = make();
  ASSERT_EQ(gshgeom_path_move_to(a, P(1, 1)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(a, P(2, 2)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_move_to(a, P(5, 6)), GSHGEOM_OK);
  GSHGEOM_Path * c = nullptr;
  ASSERT_EQ(gshgeom_path_copy(a, nullptr, nullptr, &c), GSHGEOM_OK);
  ASSERT_TRUE(gshgeom_path_equal(a, c));
  for (int op = 0; op < 3; op++) {
    GSHGEOM_Result ra = GSHGEOM_OK, rc = GSHGEOM_OK;
    if (op == 0) {
      ra = gshgeom_path_line_to(a, P(7, 7));
      rc = gshgeom_path_line_to(c, P(7, 7));
      ra = gshgeom_path_close(a);
      rc = gshgeom_path_close(c);
    } else if (op == 1) {
      ra = gshgeom_path_close(a);
      rc = gshgeom_path_close(c);
    } else {
      GSHGEOM_Arc_Form f = {P(5, 7), P(1, 0), P(0, 1), GSHGEOM_ANGLE_DEGREES, 270, 90, GSHGEOM_ARC_CIRCULAR};
      ra = gshgeom_path_arc_conjugate_to(a, &f, nullptr);
      rc = gshgeom_path_arc_conjugate_to(c, &f, nullptr);
    }
    EXPECT_EQ(ra, rc) << op;
    EXPECT_TRUE(gshgeom_path_equal(a, c)) << op;
  }
  // the readback of the start of an element is the point the element starts from
  GSHGEOM_Path * g = good_arc_path();
  EXPECT_EQ(gshgeom_path_element_start(g, 0).x, 1.0);
  EXPECT_EQ(gshgeom_path_element_start(g, 1).x, 1.0);
  EXPECT_EQ(gshgeom_path_element_start(g, 1).y, 0.0);
  GSHGEOM_Path * gc = nullptr;
  ASSERT_EQ(gshgeom_path_copy(g, nullptr, nullptr, &gc), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_hash(g), gshgeom_path_hash(gc));
  gshgeom_path_destroy(gc);
  gshgeom_path_destroy(g);
  gshgeom_path_destroy(a);
  gshgeom_path_destroy(c);
}

TEST(ArcWalkerEdges, ACopyThatFailsLeavesNothingBehind) {
  GSHGEOM_Path * g = good_arc_path();
  for (long fail = 0; fail < 40; fail++) {
    tracker::Tracker t;
    t.fail_from = fail;
    t.persistent = true;
    GSHGEOM_Allocator al = tracker::allocator_for(&t);
    GSHGEOM_Path * c = nullptr;
    const GSHGEOM_Result r = gshgeom_path_copy(g, &al, nullptr, &c);
    if (r == GSHGEOM_OK) {
      gshgeom_path_destroy(c);
    } else {
      EXPECT_EQ(r, GSHGEOM_ERR_OOM) << fail;
      EXPECT_EQ(c, nullptr);
    }
    EXPECT_TRUE(t.live.empty()) << "a copy that failed at call " << fail << " left " << t.live.size() << " allocations";
  }
  gshgeom_path_destroy(g);
}

TEST(ArcWalkerEdges, PathsThatDifferOnlyInTheirVerbsAreNotEqual) {
  GSHGEOM_Path * a = make();
  GSHGEOM_Path * b = make();
  ASSERT_EQ(gshgeom_path_move_to(a, P(0, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(a, P(1, 1)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_move_to(b, P(0, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_move_to(b, P(1, 1)), GSHGEOM_OK);
  EXPECT_FALSE(gshgeom_path_equal(a, b));
  gshgeom_path_destroy(a);
  gshgeom_path_destroy(b);
}

namespace {

// a stream whose write number `fail_at` fails (and no other)
struct Spoil {
  long writes = 0;
  long fail_at = -1;
};

ssize_t spoil_write(void * cookie, const char *, size_t n) {
  Spoil * s = static_cast<Spoil *>(cookie);
  if (s->writes++ == s->fail_at) {
    return -1;
  }
  return (ssize_t)n;
}

}  // namespace

TEST(ArcWalkerEdges, ADumpReportsAWriteThatFailsWhereverItFails) {
  GSHGEOM_Path * g = good_arc_path();
  long total = 0;
  {
    Spoil s;
    cookie_io_functions_t io = {nullptr, spoil_write, nullptr, nullptr};
    FILE * f = fopencookie(&s, "w", io);
    ASSERT_NE(f, nullptr);
    char buf[16];
    setvbuf(f, buf, _IOFBF, sizeof buf);  // a buffer of 16 bytes: a write that fails surfaces in the fprintf that overflowed it
    ASSERT_EQ(gshgeom_path_dump(g, f), GSHGEOM_OK);
    fclose(f);
    total = s.writes;
  }
  ASSERT_GT(total, 10);
  for (long k = 0; k < total; k++) {
    Spoil s;
    s.fail_at = k;
    cookie_io_functions_t io = {nullptr, spoil_write, nullptr, nullptr};
    FILE * f = fopencookie(&s, "w", io);
    ASSERT_NE(f, nullptr);
    char buf[16];
    setvbuf(f, buf, _IOFBF, sizeof buf);  // a buffer of 16 bytes: a write that fails surfaces in the fprintf that overflowed it
    EXPECT_EQ(gshgeom_path_dump(g, f), GSHGEOM_ERR_IO) << "the write that fails is number " << k << " of " << total;
    fclose(f);
  }
  gshgeom_path_destroy(g);
}

TEST(ArcWalkerEdges, AShapeWithAnyPointPastTheLimitIsRefusedBeforeAnythingIsAllocated) {
  // the refusal is made from the numbers alone: the allocator is not called (a move appended and rolled back would call it)
  tracker::Tracker t;
  GSHGEOM_Allocator al = tracker::allocator_for(&t);
  GSHGEOM_Path * p = nullptr;
  ASSERT_EQ(gshgeom_path_create(&al, nullptr, 2, &p), GSHGEOM_OK);
  const long before = t.calls;
  const double big = 0x1p101, half = 0x1p50;
  const GSHGEOM_Corner round1[4] = {{0.25, 0.25}, {0.25, 0.25}, {0.25, 0.25}, {0.25, 0.25}};
  EXPECT_EQ(gshgeom_path_add_ellipse(p, 0, 0, 1, big, 0, 0), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(gshgeom_path_add_ellipse(p, -0x1p100, 0, half, 1, 0, 0), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(gshgeom_path_add_ellipse(p, 0, 0x1p100, 1, half, 0, 0), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(gshgeom_path_add_rounded_rect(p, 0, 0, big, 1, round1, 0, 0), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(gshgeom_path_add_rounded_rect(p, 0x1p100, 0, 0x1p49, 1, round1, 0, 0), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(gshgeom_path_add_rounded_rect(p, 0, -0x1p100 - 0x1p49, 1, 0x1p50, round1, 0, 0), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(t.calls, before);
  gshgeom_path_destroy(p);
  EXPECT_TRUE(t.live.empty());
}

TEST(ArcWalkerEdges, AShapeWithAnyPointPastTheLimitIsRefusedBeforeAnythingIsAppended) {
  const double big = 0x1p101;
  const double half = 0x1p50;
  GSHGEOM_Path * p = make();
  // ellipse: the ends are (cx + rx, cy), (cx, cy + ry), (cx - rx, cy), (cx, cy - ry); each in turn the only one past the limit
  struct E { double cx, cy, rx, ry; };
  const E cases[] = {
      {0x1p100 - 1, 0, 0x1p99, 1},   // e0 = cx + rx is past
      {0, 0, 1, big},                // e1.y and e3.y past (cy +- ry)
      {-0x1p100, 0, half, 1},        // only cx - rx past
      {0, -0x1p100, 1, half},        // only cy - ry past
      {0, 0x1p100, 1, half},         // only cy + ry past
  };
  for (const E & c : cases) {
    EXPECT_EQ(gshgeom_path_add_ellipse(p, c.cx, c.cy, c.rx, c.ry, 0, 0), GSHGEOM_ERR_LIMIT);
    EXPECT_EQ(gshgeom_path_element_count(p), 0u) << c.cx << " " << c.cy << " " << c.rx << " " << c.ry;
  }
  // rounded rectangle: each of the twelve points in turn the only one (or one of the few) past the limit
  const GSHGEOM_Corner round1[4] = {{0.25, 0.25}, {0.25, 0.25}, {0.25, 0.25}, {0.25, 0.25}};
  struct R { double x, y, w, h; };
  const R rects[] = {
      {0, 0, big, 1}, {0, 0, 1, big}, {-big, 0, 1, 1}, {0, -big, 1, 1},
      {0x1p100, 0, 0x1p49, 1},          // the left edge at the limit, the right past it
      {0, 0x1p100, 1, 0x1p49},
      {-0x1p100 - 0x1p49, 0, 0x1p50, 1},  // the left edge past, the right in range
      {0, -0x1p100 - 0x1p49, 1, 0x1p50},
  };
  for (const R & r : rects) {
    const GSHGEOM_Result e = gshgeom_path_add_rounded_rect(p, r.x, r.y, r.w, r.h, round1, 0, 0);
    EXPECT_EQ(e, GSHGEOM_ERR_LIMIT) << r.x << " " << r.y << " " << r.w << " " << r.h;
    EXPECT_EQ(gshgeom_path_element_count(p), 0u) << r.x << " " << r.y << " " << r.w << " " << r.h;
  }
  gshgeom_path_destroy(p);
}

TEST(ArcWalkerEdges, TheHalfDiametersOfAnEndpointArcAlwaysHaveAPositiveArea) {
  // The endpoint form builds u = (rx cos phi, rx sin phi) and v = (-ry sin phi, ry cos phi) from positive radii, and refuses
  // them when their cross product is not positive.  It is positive for every radii and angle (each term of the exact cross
  // product is a product of two like-signed roundings), including the radii at the bottom of the range, where the sum in
  // doubles underflows to zero: the refusal cannot be reached, so the three variants of it are the same program.
  uint64_t seed = 88172645463325252ull;
  auto next = [&]() {
    seed ^= seed << 13;
    seed ^= seed >> 7;
    seed ^= seed << 17;
    return seed;
  };
  for (int i = 0; i < 200000; i++) {
    const int e1 = (int)(next() % 2200) - 1100, e2 = (int)(next() % 2200) - 1100;
    const double rx = std::ldexp(1.0 + (double)(next() >> 12) / 4503599627370496.0, e1);
    const double ry = std::ldexp(1.0 + (double)(next() >> 12) / 4503599627370496.0, e2);
    if (!(rx > 0.0) || !(ry > 0.0) || !std::isfinite(rx) || !std::isfinite(ry)) {
      continue;
    }
    const GSHGEOM_Angle_Unit unit = (GSHGEOM_Angle_Unit)(next() % GSHGEOM_ANGLE_UNIT_COUNT);
    const double phi = (double)((int64_t)(next() % 2000001) - 1000000) / 977.0;
    double cp, sp;
    ASSERT_EQ(gshgeom_sincos(unit, phi, &cp, &sp), GSHGEOM_OK);
    const GSHGEOM_Point u = P(rx * cp, rx * sp), v = P(ry * (0.0 - sp), ry * cp);
    int sign = 0;
    ASSERT_EQ(gshgeom_orient2d(P(0, 0), u, v, &sign), GSHGEOM_OK);
    ASSERT_EQ(sign, 1) << rx << " " << ry << " " << phi;
  }
}

static_assert(GSHGEOM_ARC_DECOMP == 0x1p-44, "the decomposition's bound is 2^-44");
static_assert(GSHGEOM_ARC_MAX_PIECES == 4, "an arc is at most four pieces");

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
