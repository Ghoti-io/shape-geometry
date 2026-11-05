/**
 * @file
 *
 * The device box (A3, CU15): `gshgeom_device_box` is the lattice box of the device path
 * `gshgeom_device_path_from_model` makes, without making it, and Lemma CU15 (device.h) says every
 * point the sink draws lies in it.  The tests state it three ways: the box equals the extent of the
 * device path's own elements (control points included), every vertex of the flattened path lies in
 * it, and no point of the filled region, under either rule, lies outside it (sampled on a ring just
 * outside the box and at the polyline's own edge points), over every size class of curve, collinear
 * and near-collinear controls, rational curves of extreme weights, arcs, and random maps.
 *
 * Copyright 2026 by Corey Pennycuff
 */

#include "test_helpers.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

namespace {

__extension__ typedef __int128 i128;

GSHGEOM_Point P(double x, double y) {
  return GSHGEOM_Point{x, y, 0.0};
}

GSHGEOM_Map ident() {
  GSHGEOM_Map m;
  EXPECT_EQ(gshgeom_map_identity(2, &m), GSHGEOM_OK);
  return m;
}

struct Box {
  bool empty = true;
  int64_t x0 = 0, y0 = 0, x1 = 0, y1 = 0;
};

// The extent of every point the device path stores, control points included.
Box extent_of(const GSHGEOM_Device_Path * dp, bool drew_out[1]) {
  Box b;
  bool drew = false;
  size_t n = gshgeom_device_path_element_count(dp);
  for (size_t i = 0; i < n; i++) {
    GSHGEOM_Device_Element e;
    EXPECT_EQ(gshgeom_device_path_element(dp, i, &e), GSHGEOM_OK);
    drew |= e.verb != GSHGEOM_DEVICE_MOVE && e.verb != GSHGEOM_DEVICE_CLOSE;
    for (unsigned k = 0; k < e.count; k++) {
      if (e.verb == GSHGEOM_DEVICE_CLOSE) {
        continue;  // the start, already a point of the subpath
      }
      const GSHGEOM_Device_Point & p = e.p[k];
      if (b.empty) {
        b = Box{false, p.x, p.y, p.x, p.y};
      } else {
        b.x0 = std::min<int64_t>(b.x0, p.x);
        b.x1 = std::max<int64_t>(b.x1, p.x);
        b.y0 = std::min<int64_t>(b.y0, p.y);
        b.y1 = std::max<int64_t>(b.y1, p.y);
      }
    }
  }
  drew_out[0] = drew;
  return b;
}

struct PathHolder {
  GSHGEOM_Path * p = nullptr;
  PathHolder() { EXPECT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &p), GSHGEOM_OK); }
  ~PathHolder() { gshgeom_path_destroy(p); }
};

// The box must equal the device path's extent exactly, and the flat path's vertices must be in it.
void check_against_path(const GSHGEOM_Path * path, const GSHGEOM_Map & map, const char * what) {
  GSHGEOM_Device_Box box;
  GSHGEOM_Device_Path * dp = nullptr;
  GSHGEOM_Result rb = gshgeom_device_box(path, &map, &box);
  GSHGEOM_Result rp = gshgeom_device_path_from_model(path, &map, nullptr, nullptr, &dp);
  if (rp == GSHGEOM_ERR_LIMIT && rb == GSHGEOM_ERR_LIMIT) {
    return;
  }
  ASSERT_EQ(rp, GSHGEOM_OK) << what;
  ASSERT_EQ(rb, GSHGEOM_OK) << what;
  bool drew[1];
  Box e = extent_of(dp, drew);
  if (!drew[0]) {
    EXPECT_EQ(box.empty, 1) << what;
  } else {
    ASSERT_EQ(box.empty, 0) << what;
    EXPECT_EQ(box.x0, e.x0) << what;
    EXPECT_EQ(box.y0, e.y0) << what;
    EXPECT_EQ(box.x1, e.x1) << what;
    EXPECT_EQ(box.y1, e.y1) << what;
  }
  // Lemma CU15: every vertex of the polyline the integer flattener makes is in the box.
  GSHGEOM_Flat_Path * fp = nullptr;
  ASSERT_EQ(gshgeom_device_path_flatten(dp, nullptr, nullptr, &fp), GSHGEOM_OK) << what;
  const GSHGEOM_Device_Point * pts = gshgeom_flat_path_points(fp);
  size_t n = gshgeom_flat_path_point_count(fp);
  for (size_t i = 0; i < n && drew[0]; i++) {
    EXPECT_GE(pts[i].x, box.x0) << what;
    EXPECT_LE(pts[i].x, box.x1) << what;
    EXPECT_GE(pts[i].y, box.y0) << what;
    EXPECT_LE(pts[i].y, box.y1) << what;
  }
  gshgeom_flat_path_destroy(fp);
  gshgeom_device_path_destroy(dp);
}

}  // namespace

TEST(DeviceBox, EmptyPathsAndMovesOnlyAreEmpty) {
  PathHolder h;
  GSHGEOM_Map m = ident();
  GSHGEOM_Device_Box b = {7, 7, 7, 7, 7};
  ASSERT_EQ(gshgeom_device_box(h.p, &m, &b), GSHGEOM_OK);
  EXPECT_EQ(b.empty, 1);
  EXPECT_EQ(b.x0, 0);
  ASSERT_EQ(gshgeom_path_move_to(h.p, P(3, 4)), GSHGEOM_OK);
  b = {7, 7, 7, 7, 7};
  ASSERT_EQ(gshgeom_device_box(h.p, &m, &b), GSHGEOM_OK);
  EXPECT_EQ(b.empty, 1);
  EXPECT_EQ(b.x0, 0) << "an empty box has zero corners, not the point of the move";
  EXPECT_EQ(b.y0, 0);
  EXPECT_EQ(b.x1, 0);
  EXPECT_EQ(b.y1, 0);
  ASSERT_EQ(gshgeom_path_move_to(h.p, P(9, 9)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(h.p, P(9, 12)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_device_box(h.p, &m, &b), GSHGEOM_OK);
  EXPECT_EQ(b.empty, 0);
  // The lone move that began the first subpath is a point of the box (the box of the device path's
  // elements): the interior test ignores a box without width.
  EXPECT_EQ(b.x1, 9 * 64);
  EXPECT_EQ(b.y1, 12 * 64);
}

TEST(DeviceBox, RectangleUnderMapsIsTheLatticeBox) {
  PathHolder h;
  ASSERT_EQ(gshgeom_path_move_to(h.p, P(1, 2)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(h.p, P(4.5, 2)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(h.p, P(4.5, 7.25)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(h.p, P(1, 7.25)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_close(h.p), GSHGEOM_OK);
  GSHGEOM_Map m = ident();
  GSHGEOM_Device_Box b;
  ASSERT_EQ(gshgeom_device_box(h.p, &m, &b), GSHGEOM_OK);
  EXPECT_EQ(b.x0, 64);
  EXPECT_EQ(b.y0, 128);
  EXPECT_EQ(b.x1, 288);
  EXPECT_EQ(b.y1, 464);
  GSHGEOM_Map t;
  ASSERT_EQ(gshgeom_map_translate(2, 0.5, -1, 0, &t), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_device_box(h.p, &t, &b), GSHGEOM_OK);
  EXPECT_EQ(b.x0, 96);
  EXPECT_EQ(b.y0, 64);
  // a mirror reverses which corner is the minimum
  GSHGEOM_Map s;
  ASSERT_EQ(gshgeom_map_scale(2, -1, 1, 1, &s), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_device_box(h.p, &s, &b), GSHGEOM_OK);
  EXPECT_EQ(b.x0, -288);
  EXPECT_EQ(b.x1, -64);
}

TEST(DeviceBox, ControlPointsOfQuadraticsAndCubicsAreInTheBox) {
  PathHolder h;
  ASSERT_EQ(gshgeom_path_move_to(h.p, P(0, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_quadratic_to(h.p, P(5, -3), P(10, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_cubic_to(h.p, P(12, 6), P(-2, 8), P(0, 0)), GSHGEOM_OK);
  GSHGEOM_Map m = ident();
  GSHGEOM_Device_Box b;
  ASSERT_EQ(gshgeom_device_box(h.p, &m, &b), GSHGEOM_OK);
  EXPECT_EQ(b.x0, -2 * 64);
  EXPECT_EQ(b.x1, 12 * 64);
  EXPECT_EQ(b.y0, -3 * 64);
  EXPECT_EQ(b.y1, 8 * 64);
}

TEST(DeviceBox, RefusesWhatTheDevicePathRefuses) {
  PathHolder h;
  GSHGEOM_Map m = ident();
  GSHGEOM_Device_Box b = {5, 5, 5, 5, 5};
  EXPECT_EQ(gshgeom_device_box(nullptr, &m, &b), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_device_box(h.p, nullptr, &b), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_device_box(h.p, &m, nullptr), GSHGEOM_ERR_INVALID);
  GSHGEOM_Map m3;
  ASSERT_EQ(gshgeom_map_identity(3, &m3), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_device_box(h.p, &m3, &b), GSHGEOM_ERR_INVALID);
  GSHGEOM_Path * p3 = nullptr;
  ASSERT_EQ(gshgeom_path_create(nullptr, nullptr, 3, &p3), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_device_box(p3, &m, &b), GSHGEOM_ERR_INVALID);
  gshgeom_path_destroy(p3);
  ASSERT_EQ(gshgeom_path_move_to(h.p, P(0, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(h.p, P(1e12, 0)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_device_box(h.p, &m, &b), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(b.empty, 5);  // written only on success
  GSHGEOM_Map bad = m;
  bad.m[0][0] = std::nan("");
  EXPECT_EQ(gshgeom_device_box(h.p, &bad, &b), GSHGEOM_ERR_INVALID);
}

TEST(DeviceBox, ArcsAreBoxedByTheirDevicePathNotTheirChord) {
  PathHolder h;
  ASSERT_EQ(gshgeom_path_move_to(h.p, P(0, 0)), GSHGEOM_OK);
  GSHGEOM_Arc_Kind kind;
  // A half circle of radius 10 from (0,0) to (20,0): the chord has no height.
  ASSERT_EQ(gshgeom_path_arc_endpoint_to(h.p, P(20, 0), 10, 10, 0, GSHGEOM_ANGLE_RADIANS, 0, 1, &kind), GSHGEOM_OK);
  ASSERT_EQ(kind, GSHGEOM_ARC_KIND_ARC);
  GSHGEOM_Map m = ident();
  GSHGEOM_Device_Box b;
  ASSERT_EQ(gshgeom_device_box(h.p, &m, &b), GSHGEOM_OK);
  EXPECT_GE(b.y1 - b.y0, 9 * 64);  // the arc bulges by 10 units, within the lattice
  check_against_path(h.p, m, "half circle");
  // An arc whose bulge is inside a piece, not at its ends: 45 to 135 degrees on a circle of radius 10.  The
  // chord has no height to speak of; the device path's vertices are the conic route's.
  PathHolder q;
  ASSERT_EQ(gshgeom_path_move_to(q.p, P(7.0710678118654755, 7.0710678118654755)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_arc_endpoint_to(q.p, P(-7.0710678118654755, 7.0710678118654755), 10, 10, 0,
                GSHGEOM_ANGLE_RADIANS, 0, 1, &kind), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_device_box(q.p, &m, &b), GSHGEOM_OK);
  EXPECT_GE(b.y1, 10 * 64 - 8);  // the top of the arc, within the flattener's tolerance
  EXPECT_GE(b.y1 - b.y0, 2 * 64) << "the arc rises above its chord by about 2.9 units";
  GSHGEOM_Map r;
  ASSERT_EQ(gshgeom_map_rotate(GSHGEOM_ANGLE_DEGREES, 37, &r), GSHGEOM_OK);
  check_against_path(h.p, r, "rotated half circle");
}

// ---- the dense statement of the lemma -------------------------------------------------------------

namespace {

// The nonzero winding and the parity of the closed polygons of a flat path at a point given in 1/128
// lattice units (so a sample can lie between lattice points).  Crossing number with exact integers.
void windings(const GSHGEOM_Flat_Path * fp, int64_t px2, int64_t py2, int * wn, int * parity) {
  *wn = 0;
  *parity = 0;
  const GSHGEOM_Device_Point * pts = gshgeom_flat_path_points(fp);
  size_t ns = gshgeom_flat_path_subpath_count(fp);
  for (size_t s = 0; s < ns; s++) {
    GSHGEOM_Flat_Subpath sp;
    EXPECT_EQ(gshgeom_flat_path_subpath(fp, s, &sp), GSHGEOM_OK);
    for (size_t i = 0; i < sp.count; i++) {
      const GSHGEOM_Device_Point & a = pts[sp.first + i];
      const GSHGEOM_Device_Point & b = pts[sp.first + (i + 1) % sp.count];
      int64_t ay = (int64_t)a.y * 2, by = (int64_t)b.y * 2, ax = (int64_t)a.x * 2, bx = (int64_t)b.x * 2;
      if ((ay <= py2) == (by <= py2)) {
        continue;
      }
      // x of the edge at py2, compared with px2 without division: sign of cross product
      i128 cross = (i128)(bx - ax) * (py2 - ay) - (i128)(px2 - ax) * (by - ay);
      bool left = (by > ay) ? cross > 0 : cross < 0;  // the point is left of the edge: crossing to the right
      if (left) {
        *wn += (by > ay) ? 1 : -1;
        *parity ^= 1;
      }
    }
  }
}

struct Gen {
  std::mt19937_64 rng;
  explicit Gen(uint64_t seed) : rng(seed) {}
  double u(double lo, double hi) { return lo + (hi - lo) * (double)(rng() >> 11) / 9007199254740992.0; }
  int pick(int n) { return (int)(rng() % (uint64_t)n); }
  // A coordinate of a size class: tiny (sub lattice), pixel, glyph, page, huge.
  double coord(int cls) {
    static const double s[] = {0.01, 1.5, 40, 3000, 4.0e6};
    return u(-s[cls], s[cls]);
  }
};

void random_path(Gen & g, GSHGEOM_Path * path, int cls) {
  int nsub = 1 + g.pick(3);
  for (int s = 0; s < nsub; s++) {
    GSHGEOM_Point cur = P(g.coord(cls), g.coord(cls));
    ASSERT_EQ(gshgeom_path_move_to(path, cur), GSHGEOM_OK);
    int nseg = 1 + g.pick(5);
    for (int k = 0; k < nseg; k++) {
      GSHGEOM_Point a = P(g.coord(cls), g.coord(cls)), b = P(g.coord(cls), g.coord(cls)), c = P(g.coord(cls), g.coord(cls));
      int kind = g.pick(9);
      if (kind == 7) {  // collinear controls
        double t1 = g.u(-1, 2), t2 = g.u(-1, 2);
        a = P(cur.x + (c.x - cur.x) * t1, cur.y + (c.y - cur.y) * t1);
        b = P(cur.x + (c.x - cur.x) * t2, cur.y + (c.y - cur.y) * t2);
        kind = 3;
      }
      if (kind == 8) {  // near collinear
        double e = g.u(-1, 1) * 1e-3;
        a = P(cur.x + (c.x - cur.x) * 0.3 + e, cur.y + (c.y - cur.y) * 0.3);
        b = P(cur.x + (c.x - cur.x) * 0.7, cur.y + (c.y - cur.y) * 0.7 - e);
        kind = 3;
      }
      switch (kind) {
      case 0: case 1: ASSERT_EQ(gshgeom_path_line_to(path, a), GSHGEOM_OK); cur = a; break;
      case 2: ASSERT_EQ(gshgeom_path_quadratic_to(path, a, b), GSHGEOM_OK); cur = b; break;
      case 3: ASSERT_EQ(gshgeom_path_cubic_to(path, a, b, c), GSHGEOM_OK); cur = c; break;
      case 4: {
        double w[3] = {1, std::exp2(g.u(-30, 30)), std::exp2(g.u(-30, 30))};
        ASSERT_EQ(gshgeom_path_rational_quadratic_to(path, a, b, w), GSHGEOM_OK);
        cur = b;
        break;
      }
      case 5: {
        double w[4] = {1, std::exp2(g.u(-20, 20)), std::exp2(g.u(-20, 20)), std::exp2(g.u(-5, 5))};
        ASSERT_EQ(gshgeom_path_rational_cubic_to(path, a, b, c, w), GSHGEOM_OK);
        cur = c;
        break;
      }
      default: {
        GSHGEOM_Arc_Kind kind2;
        double rx = std::fabs(g.coord(cls)) + 1e-3, ry = std::fabs(g.coord(cls)) + 1e-3;
        ASSERT_EQ(gshgeom_path_arc_endpoint_to(path, a, rx, ry, g.u(-3, 3), GSHGEOM_ANGLE_RADIANS,
                      g.pick(2), g.pick(2), &kind2), GSHGEOM_OK);
        cur = a;
        break;
      }
      }
    }
    if (g.pick(2)) {
      ASSERT_EQ(gshgeom_path_close(path), GSHGEOM_OK);
    }
  }
}

GSHGEOM_Map random_map(Gen & g) {
  GSHGEOM_Map m = ident(), t, r, c;
  switch (g.pick(4)) {
  case 0: return m;
  case 1:
    EXPECT_EQ(gshgeom_map_scale(2, g.u(0.2, 3), g.u(-3, 3), 1, &t), GSHGEOM_OK);
    return t;
  case 2:
    EXPECT_EQ(gshgeom_map_rotate(GSHGEOM_ANGLE_DEGREES, g.u(-360, 360), &r), GSHGEOM_OK);
    EXPECT_EQ(gshgeom_map_translate(2, g.u(-50, 50), g.u(-50, 50), 0, &t), GSHGEOM_OK);
    EXPECT_EQ(gshgeom_map_compose(&t, &r, &c), GSHGEOM_OK);
    return c;
  default:
    EXPECT_EQ(gshgeom_map_scale(2, g.u(0.1, 8), g.u(0.1, 8), 1, &t), GSHGEOM_OK);
    EXPECT_EQ(gshgeom_map_rotate(GSHGEOM_ANGLE_DEGREES, g.u(-360, 360), &r), GSHGEOM_OK);
    EXPECT_EQ(gshgeom_map_compose(&r, &t, &c), GSHGEOM_OK);
    return c;
  }
}

}  // namespace

TEST(DeviceBox, DenseLemmaCU15EveryVertexAndEveryFilledPointIsInTheBox) {
  Gen g(0xA3C11B0F);
  int checked = 0, filled_outside_samples = 0;
  for (int iter = 0; iter < 1500; iter++) {
    PathHolder h;
    int cls = iter % 5;
    random_path(g, h.p, cls);
    GSHGEOM_Map m = random_map(g);
    check_against_path(h.p, m, "random");
    if (HasFatalFailure()) {
      return;
    }
    GSHGEOM_Device_Box box;
    if (gshgeom_device_box(h.p, &m, &box) != GSHGEOM_OK || box.empty) {
      continue;
    }
    GSHGEOM_Device_Path * dp = nullptr;
    ASSERT_EQ(gshgeom_device_path_from_model(h.p, &m, nullptr, nullptr, &dp), GSHGEOM_OK);
    GSHGEOM_Flat_Path * fp = nullptr;
    ASSERT_EQ(gshgeom_device_path_flatten(dp, nullptr, nullptr, &fp), GSHGEOM_OK);
    checked++;
    // Samples on a ring one half lattice unit outside the box (coordinates in half units): no
    // sample may be inside the region under either rule.  Also the four sides' lines at many points.
    int64_t X0 = (int64_t)box.x0 * 2 - 1, X1 = (int64_t)box.x1 * 2 + 1;
    int64_t Y0 = (int64_t)box.y0 * 2 - 1, Y1 = (int64_t)box.y1 * 2 + 1;
    for (int k = 0; k < 40; k++) {
      int64_t tx = X0 + (int64_t)(g.u(0, 1) * (double)(X1 - X0));
      int64_t ty = Y0 + (int64_t)(g.u(0, 1) * (double)(Y1 - Y0));
      int64_t pts[4][2] = {{tx, Y0}, {tx, Y1}, {X0, ty}, {X1, ty}};
      for (auto & p : pts) {
        int wn, par;
        windings(fp, p[0], p[1], &wn, &par);
        EXPECT_EQ(wn, 0) << "iter " << iter;
        EXPECT_EQ(par, 0) << "iter " << iter;
        filled_outside_samples++;
      }
    }
    gshgeom_flat_path_destroy(fp);
    gshgeom_device_path_destroy(dp);
  }
  EXPECT_GT(checked, 1000);
  EXPECT_GT(filled_outside_samples, 100000);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
