/**
 * @file
 *
 * Hit testing against a flattened path (hit.h, HT1 to HT6).  The hand cases pin each rule; the random cases compare every point
 * of a grid, on and off the edges, with a reference that shares no code with the library: Hormann and Agathos's quadrant count,
 * which sums the quadrant changes of the contour around the point and decides a half-turn by the sign of a cross product.
 */

#include "test_helpers.h"

#include <cstdint>
#include <random>
#include <utility>
#include <vector>

/**
 * @file
 *
 * The width profile (stroke.h, SW1 to SW9): the order of the refusals, the statuses, the other forms that refuse a profile, the outline
 * against an oracle that knows only the definition (a point is in the region when some point of the centerline has a disc of its
 * radius reaching it, tested by a convex minimization per piece, in model space), constant and tapering and swelling and jumping
 * profiles, closed and several subpaths, points, the maps, and every failing allocation.
 */

#include <gtest/gtest.h>

#include "../alloc_tracker.h"

#include <ghoti.io/shape-geometry/shape-geometry.h>

#include <cmath>
#include <cstring>
#include <random>
#include <vector>

namespace {

using Stop = GSHGEOM_Width_Stop;

GSHGEOM_Stroke round_stroke(double width, const GSHGEOM_Width_Profile * wp) {
  GSHGEOM_Stroke s;
  gshgeom_stroke_init(&s);
  s.width = width;
  s.join = GSHGEOM_STROKE_JOIN_ROUND;
  s.cap = GSHGEOM_STROKE_CAP_ROUND;
  s.width_profile = wp;
  return s;
}

GSHGEOM_Map identity() {
  GSHGEOM_Map m;
  gshgeom_map_identity(2, &m);
  return m;
}

struct Pt {
  double x, y;
};

GSHGEOM_Path * polyline(const std::vector<std::vector<Pt>> & subs, const std::vector<bool> & closed = {}) {
  GSHGEOM_Path * p = nullptr;
  gshgeom_path_create(nullptr, nullptr, 2, &p);
  for (size_t i = 0; i < subs.size(); i++) {
    for (size_t k = 0; k < subs[i].size(); k++) {
      GSHGEOM_Point q{subs[i][k].x, subs[i][k].y, 0};
      if (k == 0) {
        gshgeom_path_move_to(p, q);
      } else {
        gshgeom_path_line_to(p, q);
      }
    }
    if (i < closed.size() && closed[i]) {
      gshgeom_path_close(p);
    }
  }
  return p;
}

struct Outline {
  GSHGEOM_Stroke_Outline * o = nullptr;
  ~Outline() { gshgeom_stroke_outline_destroy(o); }
};

/** The winding of the body at a device point in pixels. */
int64_t winding(const GSHGEOM_Stroke_Outline * o, double x, double y) {
  const GSHGEOM_Flat_Path * body = gshgeom_stroke_outline_body(o);
  int64_t w = 0;
  int on = 0;
  GSHGEOM_Device_Point p{(int32_t)std::lround(x * 64.0), (int32_t)std::lround(y * 64.0)};
  EXPECT_EQ(gshgeom_flat_path_winding(body, p, &w, &on), GSHGEOM_OK);
  return on ? 99 : w;
}

/** f(s) of a profile, linear between stops, constant outside, the limit from @p side at a jump. */
double factor_of(const std::vector<Stop> & st, double s, int side) {
  const size_t n = st.size();
  if (s < st[0].position || (s == st[0].position && side < 0)) {
    return st[0].factor;
  }
  if (s > st[n - 1].position || (s == st[n - 1].position && side > 0)) {
    return st[n - 1].factor;
  }
  for (size_t i = 0; i + 1 < n; i++) {
    const double a = st[i].position, b = st[i + 1].position;
    if (a == b) {
      continue;
    }
    const bool inside = side > 0 ? (s >= a && s < b) : (s > a && s <= b);
    if (inside) {
      return st[i].factor + (st[i + 1].factor - st[i].factor) * ((s - a) / (b - a));
    }
  }
  return st[n - 1].factor;
}

/**
 * The oracle: the least over the pieces of a polyline (and over their cuts at the stops) of `|q - c| - r` along the piece, found by
 * ternary search of a convex function (the distance minus a linear radius).  Negative: inside by that much.
 */
double slack(const std::vector<std::vector<Pt>> & subs, const std::vector<bool> & closed, const std::vector<Stop> & st, double width, Pt q) {
  double best = 1e300;
  for (size_t si = 0; si < subs.size(); si++) {
    const std::vector<Pt> & v = subs[si];
    const bool cl = si < closed.size() && closed[si];
    if (v.size() == 1) {
      best = std::min(best, std::hypot(q.x - v[0].x, q.y - v[0].y) - width * std::max(factor_of(st, 0, -1), factor_of(st, 0, 1)) / 2.0);
      continue;
    }
    double s = 0;
    const size_t segs = cl ? v.size() : v.size() - 1;
    for (size_t i = 0; i < segs; i++) {
      const Pt a = v[i], b = v[(i + 1) % v.size()];
      const double len = std::hypot(b.x - a.x, b.y - a.y);
      // cut points: the segment's ends and every stop strictly inside it
      std::vector<double> cuts{s};
      for (const Stop & k : st) {
        if (k.position > s && k.position < s + len) {
          cuts.push_back(k.position);
        }
      }
      cuts.push_back(s + len);
      for (size_t c = 0; c + 1 < cuts.size(); c++) {
        const double s0 = cuts[c], s1 = cuts[c + 1];
        const double r0 = width * factor_of(st, s0, 1) / 2.0, r1 = width * factor_of(st, s1, -1) / 2.0;
        auto g = [&](double t) {
          const double ss = s0 + (s1 - s0) * t;
          const double u = len > 0 ? (ss - s) / len : 0.0;
          const double cx = a.x + (b.x - a.x) * u, cy = a.y + (b.y - a.y) * u;
          return std::hypot(q.x - cx, q.y - cy) - (r0 + (r1 - r0) * t);
        };
        double lo = 0, hi = 1;
        for (int it = 0; it < 80; it++) {
          const double m1 = lo + (hi - lo) / 3, m2 = hi - (hi - lo) / 3;
          if (g(m1) < g(m2)) {
            hi = m2;
          } else {
            lo = m1;
          }
        }
        best = std::min(best, std::min(g((lo + hi) / 2), std::min(g(0), g(1))));
      }
      // a jump exactly at the segment's start or end has both limits
      best = std::min(best, std::hypot(q.x - a.x, q.y - a.y) - width * factor_of(st, s, -1) / 2.0);
      s += len;
      best = std::min(best, std::hypot(q.x - b.x, q.y - b.y) - width * factor_of(st, s, 1) / 2.0);
    }
  }
  return best;
}

struct Case {
  std::vector<std::vector<Pt>> subs;
  std::vector<bool> closed;
  std::vector<Stop> stops;
  double width;
};

/** Compare the outline of @p c under @p map with the oracle at a grid of device points; @return the points compared. */
size_t compare(const Case & c, const GSHGEOM_Map & map, double margin, const char * what, double lo, double hi, double step) {
  GSHGEOM_Width_Profile wp{c.stops.size(), c.stops.data()};
  GSHGEOM_Stroke st = round_stroke(c.width, &wp);
  GSHGEOM_Path * path = polyline(c.subs, c.closed);
  Outline out;
  EXPECT_EQ(gshgeom_stroke_expand(path, &st, &map, nullptr, nullptr, &out.o), GSHGEOM_OK) << what;
  gshgeom_path_destroy(path);
  EXPECT_EQ(gshgeom_stroke_outline_status(out.o), GSHGEOM_STROKE_DRAWN) << what;
  // the inverse of the map's linear part takes a device point to model space
  const double det = map.m[0][0] * map.m[1][1] - map.m[0][1] * map.m[1][0];
  size_t compared = 0;
  for (double y = lo; y <= hi; y += step) {
    for (double x = lo; x <= hi; x += step) {
      const double dx = x - map.t[0], dy = y - map.t[1];
      const Pt q{(map.m[1][1] * dx - map.m[0][1] * dy) / det, (-map.m[1][0] * dx + map.m[0][0] * dy) / det};
      // a slack in model units is at most `norm` times that in pixels; the band is in pixels, so scale by the least stretch
      const double s = slack(c.subs, c.closed, c.stops, c.width, q);
      const double stretch = std::sqrt(std::fabs(det)) * 0.5;
      if (std::fabs(s) * stretch < margin) {
        continue;
      }
      const int64_t w = winding(out.o, x, y);
      if (w == 99) {
        continue;
      }
      EXPECT_EQ(w != 0, s < 0) << what << " at device (" << x << ", " << y << ") model (" << q.x << ", " << q.y << ") slack " << s;
      compared++;
    }
  }
  return compared;
}

}  // namespace

TEST(StrokeProfile, TheOrderOfTheRefusals) {
  const Stop ok[2] = {{0, 1}, {10, 2}};
  GSHGEOM_Width_Profile wp{2, ok};
  EXPECT_EQ(gshgeom_stroke_check(nullptr), GSHGEOM_ERR_INVALID);
  GSHGEOM_Stroke s = round_stroke(4, &wp);
  EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_OK);
  GSHGEOM_Width_Profile none{0, ok};
  s.width_profile = &none;
  EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_ERR_INVALID) << "no stops";
  GSHGEOM_Width_Profile null_stops{1, nullptr};
  s.width_profile = &null_stops;
  EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_ERR_INVALID);
  const Stop backwards[2] = {{5, 1}, {4, 1}};
  GSHGEOM_Width_Profile b{2, backwards};
  s.width_profile = &b;
  EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_ERR_INVALID) << "decreasing position";
  const Stop negative[1] = {{0, -1}};
  GSHGEOM_Width_Profile n{1, negative};
  s.width_profile = &n;
  EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_ERR_INVALID) << "negative factor";
  const Stop nan[1] = {{0, NAN}};
  GSHGEOM_Width_Profile nn{1, nan};
  s.width_profile = &nn;
  EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_ERR_INVALID);
  const Stop neg_pos[1] = {{-1, 1}};
  GSHGEOM_Width_Profile np{1, neg_pos};
  s.width_profile = &np;
  EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_ERR_INVALID);
  /* any join, cap, dash array and alignment is drawn (SW10 to SW13); a bad profile is invalid whatever the style is */
  s = round_stroke(4, &wp);
  s.join = GSHGEOM_STROKE_JOIN_MITER;
  s.cap = GSHGEOM_STROKE_CAP_BUTT;
  s.alignment = GSHGEOM_STROKE_ALIGN_INSIDE;
  EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_OK);
  s.width_profile = &b;
  EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_ERR_INVALID) << "a bad profile is invalid whatever the style";
  const Stop huge[1] = {{0, 0x1p20}};
  GSHGEOM_Width_Profile h{1, huge};
  s = round_stroke(4, &h);
  EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_ERR_LIMIT) << "a factor past 2^16";
  const Stop tiny[1] = {{1e-30, 1}};
  GSHGEOM_Width_Profile t{1, tiny};
  s = round_stroke(4, &t);
  EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_ERR_LIMIT) << "a position below 2^-64 that is not 0";
  s = round_stroke(4, &wp);
  s.join = GSHGEOM_STROKE_JOIN_MITER;
  s.width = 1e30;
  EXPECT_EQ(gshgeom_stroke_check(&s), GSHGEOM_ERR_LIMIT) << "a width past the band";
}

TEST(StrokeProfile, StatusesAndTheFormsThatDoNotDrawAProfile) {
  const Stop st[2] = {{0, 1}, {10, 0.2}};
  GSHGEOM_Width_Profile wp{2, st};
  GSHGEOM_Path * path = polyline({{{0, 0}, {10, 0}}});
  GSHGEOM_Map m = identity();
  GSHGEOM_Stroke s = round_stroke(0, &wp);
  Outline o;
  ASSERT_EQ(gshgeom_stroke_expand(path, &s, &m, nullptr, nullptr, &o.o), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_stroke_outline_status(o.o), GSHGEOM_STROKE_EMPTY) << "SG14";
  Outline sing;
  GSHGEOM_Map flat = identity();
  flat.m[1][1] = 0.0;
  s = round_stroke(4, &wp);
  ASSERT_EQ(gshgeom_stroke_expand(path, &s, &flat, nullptr, nullptr, &sing.o), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_stroke_outline_status(sing.o), GSHGEOM_STROKE_NOT_INVERTIBLE) << "SG15 before SG14";
  GSHGEOM_Stroke_Curves * cu = nullptr;
  ASSERT_EQ(gshgeom_stroke_curves_expand(path, &s, &m, 0.1, GSHGEOM_STROKE_CURVES_RATIONAL, nullptr, nullptr, &cu), GSHGEOM_OK) << "SW14";
  EXPECT_EQ(gshgeom_stroke_curves_status(cu), GSHGEOM_STROKE_DRAWN);
  EXPECT_GT(gshgeom_path_element_count(gshgeom_stroke_curves_body(cu)), 0u);
  gshgeom_stroke_curves_destroy(cu);
  cu = nullptr;
  EXPECT_EQ(gshgeom_stroke_curves_expand(path, &s, &m, 0.01, GSHGEOM_STROKE_CURVES_RATIONAL, nullptr, nullptr, &cu), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(cu, nullptr);
  Outline three;
  GSHGEOM_Map m3;
  gshgeom_map_identity(3, &m3);
  EXPECT_NE(gshgeom_stroke_expand(path, &s, &m3, nullptr, nullptr, &three.o), GSHGEOM_OK);
  gshgeom_path_destroy(path);
}

TEST(StrokeProfile, AConstantProfileIsARoundStrokeOfThatWidth) {
  Case c;
  c.subs = {{{10, 10}, {40, 12}, {30, 36}}};
  c.stops = {{0, 1}};
  c.width = 8;
  EXPECT_GT(compare(c, identity(), 0.12, "constant", 0, 50, 0.37), 1500u);
  /* and the same bytes as the round stroke, within the pen's polygon: the same points are in or out away from the edge */
  GSHGEOM_Width_Profile wp{1, c.stops.data()};
  GSHGEOM_Stroke a = round_stroke(8, &wp);
  GSHGEOM_Stroke b = round_stroke(8, nullptr);
  GSHGEOM_Path * path = polyline(c.subs);
  GSHGEOM_Map m = identity();
  Outline oa, ob;
  ASSERT_EQ(gshgeom_stroke_expand(path, &a, &m, nullptr, nullptr, &oa.o), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_stroke_expand(path, &b, &m, nullptr, nullptr, &ob.o), GSHGEOM_OK);
  gshgeom_path_destroy(path);
  size_t differ = 0, total = 0;
  for (double y = 0; y < 50; y += 0.41) {
    for (double x = 0; x < 50; x += 0.41) {
      differ += (winding(oa.o, x, y) != 0) != (winding(ob.o, x, y) != 0);
      total++;
    }
  }
  EXPECT_LT(differ, total / 150) << "only the points within the tolerance of the edge may differ";
}

TEST(StrokeProfile, ATaperASwellAndAJumpAgreeWithTheOracle) {
  Case taper;
  taper.subs = {{{5, 25}, {45, 25}}};
  taper.stops = {{0, 1}, {40, 0}};
  taper.width = 12;
  EXPECT_GT(compare(taper, identity(), 0.12, "taper", 0, 50, 0.31), 1500u);
  Case swell;
  swell.subs = {{{5, 10}, {20, 40}, {45, 10}}};
  swell.stops = {{0, 0.2}, {20, 1.5}, {40, 0.1}, {80, 1.0}};
  swell.width = 8;
  EXPECT_GT(compare(swell, identity(), 0.12, "swell", 0, 50, 0.31), 1500u);
  Case jump;
  jump.subs = {{{5, 25}, {45, 25}}};
  jump.stops = {{0, 1}, {20, 1}, {20, 3}, {40, 3}};
  jump.width = 4;
  EXPECT_GT(compare(jump, identity(), 0.12, "jump", 0, 50, 0.31), 1500u);
  Case jump_at_vertex;
  jump_at_vertex.subs = {{{5, 10}, {25, 10}, {25, 40}}};
  jump_at_vertex.stops = {{0, 0.5}, {20, 0.5}, {20, 2}, {50, 2}};
  jump_at_vertex.width = 6;
  EXPECT_GT(compare(jump_at_vertex, identity(), 0.12, "jump at a vertex", 0, 50, 0.31), 1500u);
  Case start_jump;
  start_jump.subs = {{{10, 25}, {40, 25}}};
  start_jump.stops = {{0, 3}, {0, 0.5}, {30, 0.5}, {30, 2}};
  start_jump.width = 4;
  EXPECT_GT(compare(start_jump, identity(), 0.12, "jumps at both ends", 0, 50, 0.31), 1500u);
  Case zero;
  zero.subs = {{{5, 25}, {45, 25}}};
  zero.stops = {{0, 0}, {20, 1}, {40, 0}};
  zero.width = 10;
  EXPECT_GT(compare(zero, identity(), 0.12, "a factor of 0 at both ends", 0, 50, 0.31), 1000u);
}

TEST(StrokeProfile, ClosedSeveralAndPointSubpaths) {
  Case c;
  c.subs = {{{10, 10}, {40, 10}, {40, 40}, {10, 40}}, {{25, 25}, {27, 25}}, {{45, 5}}};
  c.closed = {true, false, false};
  c.stops = {{0, 0.3}, {60, 1.2}, {120, 0.6}};
  c.width = 6;
  EXPECT_GT(compare(c, identity(), 0.12, "closed, two and a point", 0, 50, 0.29), 2000u);
  /* a point is a disc of f(0), and an omitted degenerate subpath is nothing */
  GSHGEOM_Width_Profile wp{c.stops.size(), c.stops.data()};
  GSHGEOM_Stroke s = round_stroke(6, &wp);
  s.degenerate = GSHGEOM_STROKE_DEGENERATE_OMIT;
  GSHGEOM_Path * dot = polyline({{{20, 20}}});
  GSHGEOM_Map m = identity();
  Outline o;
  ASSERT_EQ(gshgeom_stroke_expand(dot, &s, &m, nullptr, nullptr, &o.o), GSHGEOM_OK);
  EXPECT_EQ(winding(o.o, 20, 20), 0);
  s.degenerate = GSHGEOM_STROKE_DEGENERATE_DRAW;
  Outline d;
  ASSERT_EQ(gshgeom_stroke_expand(dot, &s, &m, nullptr, nullptr, &d.o), GSHGEOM_OK);
  EXPECT_NE(winding(d.o, 20, 20), 0);
  EXPECT_NE(winding(d.o, 20.8, 20), 0) << "radius 6 * 0.3 / 2 = 0.9";
  EXPECT_EQ(winding(d.o, 21.2, 20), 0);
  gshgeom_path_destroy(dot);
}

TEST(StrokeProfile, UnderAMapTheRegionIsTheImageOfTheModelRegion) {
  Case c;
  c.subs = {{{2, 2}, {12, 4}, {6, 12}}};
  c.stops = {{0, 1}, {10, 0.25}, {20, 1.0}};
  c.width = 3;
  GSHGEOM_Map scale;
  gshgeom_map_scale(2, 4.0, 4.0, 1.0, &scale);
  EXPECT_GT(compare(c, scale, 0.15, "scale", 0, 60, 0.37), 800u);
  GSHGEOM_Map aniso, rot, m, mv;
  gshgeom_map_scale(2, 3.0, 1.2, 1.0, &aniso);
  gshgeom_map_rotate(GSHGEOM_ANGLE_DEGREES, 33, &rot);
  gshgeom_map_compose(&rot, &aniso, &m);
  gshgeom_map_translate(2, 20.5, 10.25, 0.0, &mv);
  gshgeom_map_compose(&mv, &m, &m);
  EXPECT_GT(compare(c, m, 0.15, "anisotropic and rotated", 0, 70, 0.37), 800u);
  GSHGEOM_Map mirror;
  gshgeom_map_scale(2, -3.0, 3.0, 1.0, &mirror);
  mirror.t[0] = 60;
  EXPECT_GT(compare(c, mirror, 0.15, "mirrored", 0, 70, 0.37), 800u);
}

TEST(StrokeProfile, ACurvedCenterlineAgreesWithItsDenselySampledPolyline) {
  GSHGEOM_Path * p = nullptr;
  gshgeom_path_create(nullptr, nullptr, 2, &p);
  gshgeom_path_move_to(p, GSHGEOM_Point{5, 40, 0});
  gshgeom_path_cubic_to(p, GSHGEOM_Point{5, 5, 0}, GSHGEOM_Point{45, 5, 0}, GSHGEOM_Point{45, 40, 0});
  const Stop st[3] = {{0, 0.2}, {40, 1.4}, {90, 0.4}};
  GSHGEOM_Width_Profile wp{3, st};
  GSHGEOM_Stroke s = round_stroke(7, &wp);
  GSHGEOM_Map m = identity();
  Outline o;
  ASSERT_EQ(gshgeom_stroke_expand(p, &s, &m, nullptr, nullptr, &o.o), GSHGEOM_OK);
  gshgeom_path_destroy(p);
  Case c;
  std::vector<Pt> v;
  for (int i = 0; i <= 500; i++) {
    const double t = i / 500.0, u = 1 - t;
    v.push_back({u * u * u * 5 + 3 * u * u * t * 5 + 3 * u * t * t * 45 + t * t * t * 45,
        u * u * u * 40 + 3 * u * u * t * 5 + 3 * u * t * t * 5 + t * t * t * 40});
  }
  c.subs = {v};
  c.stops = {st[0], st[1], st[2]};
  c.width = 7;
  size_t compared = 0;
  for (double y = 0; y < 50; y += 1.1) {
    for (double x = 0; x < 50; x += 1.1) {
      const double sl = slack(c.subs, {}, c.stops, 7, Pt{x, y});
      if (std::fabs(sl) < 0.3) {
        continue; /* the flattener's and the pen polygon's tolerance */
      }
      EXPECT_EQ(winding(o.o, x, y) != 0, sl < 0) << x << "," << y << " slack " << sl;
      compared++;
    }
  }
  EXPECT_GT(compared, 800u);
}

TEST(StrokeProfile, RandomPolylinesProfilesAndMapsAgreeWithTheOracle) {
  std::mt19937 rng(20261113);
  auto real = [&](double lo, double hi) { return lo + (hi - lo) * ((double)(rng() % 100001) / 100000.0); };
  size_t compared = 0;
  for (int round = 0; round < 120; round++) {
    Case c;
    const int subs = 1 + (int)(rng() % 3);
    for (int i = 0; i < subs; i++) {
      std::vector<Pt> v;
      const int n = 1 + (int)(rng() % 5);
      for (int k = 0; k < n; k++) {
        v.push_back({std::round(real(0, 48) * 8) / 8, std::round(real(0, 48) * 8) / 8});
      }
      c.subs.push_back(v);
      c.closed.push_back(n >= 3 && rng() % 3 == 0);
    }
    const int stops = 1 + (int)(rng() % 5);
    double pos = 0;
    for (int k = 0; k < stops; k++) {
      c.stops.push_back({pos, rng() % 6 == 0 ? 0.0 : real(0.1, 2.0)});
      pos += rng() % 4 == 0 ? 0.0 : real(0.5, 40);
    }
    c.width = real(1, 9);
    GSHGEOM_Map m = identity();
    if (rng() % 2) {
      GSHGEOM_Map sc, rot;
      gshgeom_map_scale(2, real(0.7, 2.5), real(0.7, 2.5), 1.0, &sc);
      gshgeom_map_rotate(GSHGEOM_ANGLE_DEGREES, real(0, 360), &rot);
      gshgeom_map_compose(&rot, &sc, &m);
      m.t[0] = real(-5, 20);
      m.t[1] = real(-5, 20);
    }
    compared += compare(c, m, 0.2, "random", -10, 90, 1.7);
  }
  EXPECT_GT(compared, 20000u);
}

TEST(StrokeProfile, EveryFailingAllocationIsAnErrorAndLeaksNothing) {
  const Stop st[3] = {{0, 0.5}, {10, 1}, {30, 0.3}};
  GSHGEOM_Width_Profile wp{3, st};
  GSHGEOM_Stroke s = round_stroke(5, &wp);
  GSHGEOM_Path * path = polyline({{{5, 5}, {30, 8}, {20, 30}}});
  GSHGEOM_Map m = identity();
  tracker::Tracker t;
  GSHGEOM_Allocator a = tracker::allocator_for(&t);
  bool done = false;
  for (long k = 0; k < 400 && !done; k++) {
    t.persistent = false;
    t.failed = false;
    t.fail_from = t.calls + k;
    GSHGEOM_Stroke_Outline * o = nullptr;
    GSHGEOM_Result r = gshgeom_stroke_expand(path, &s, &m, &a, nullptr, &o);
    t.fail_from = -1;
    if (r == GSHGEOM_OK) {
      done = true;
      EXPECT_FALSE(t.failed);
      gshgeom_stroke_outline_destroy(o);
    } else {
      EXPECT_EQ(r, GSHGEOM_ERR_OOM) << "allocation " << k;
      EXPECT_EQ(o, nullptr);
    }
    EXPECT_TRUE(t.live.empty()) << "allocation " << k << ": " << t.live.size() << " blocks leaked";
  }
  EXPECT_TRUE(done);
  gshgeom_path_destroy(path);
}

TEST(StrokeProfile, LimitsAreErrorsNotTruncations) {
  const Stop st[1] = {{0, 1}};
  GSHGEOM_Width_Profile wp{1, st};
  GSHGEOM_Stroke s = round_stroke(5, &wp);
  GSHGEOM_Path * path = polyline({{{5, 5}, {30, 8}, {20, 30}, {5, 30}, {9, 2}}});
  GSHGEOM_Map m = identity();
  GSHGEOM_Limits lim;
  gshgeom_limits_default(&lim);
  lim.max_points = 20;
  GSHGEOM_Stroke_Outline * o = nullptr;
  EXPECT_EQ(gshgeom_stroke_expand(path, &s, &m, nullptr, &lim, &o), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(o, nullptr);
  gshgeom_limits_default(&lim);
  lim.max_work = 10;
  EXPECT_EQ(gshgeom_stroke_expand(path, &s, &m, nullptr, &lim, &o), GSHGEOM_ERR_LIMIT);
  GSHGEOM_Map far = identity();
  far.t[0] = 1e12;
  EXPECT_EQ(gshgeom_stroke_expand(path, &s, &far, nullptr, nullptr, &o), GSHGEOM_ERR_LIMIT);
  gshgeom_path_destroy(path);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
