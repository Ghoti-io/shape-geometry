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
 * A width profile with any join, any cap, dashes and an alignment (stroke.h, SW10 to SW13).  The oracle is the library itself without
 * a profile: a profile of one stop and factor 1 is the constant stroke, so the two outlines must cover the same points, whatever the
 * join, cap, miter limit, pattern, alignment, width space and map.  They are compared at a grid of device points; a point where they
 * differ is accepted only when the other outline reaches the first's answer within a quarter of a pixel (the two boundaries are the
 * same curve within the tolerance).  The tapering shapes (the trapezoid, the flat end, the step of a jump) are checked against
 * their definitions by hand.
 */

#include <gtest/gtest.h>

#include "../alloc_tracker.h"

#include <ghoti.io/shape-geometry/shape-geometry.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <vector>

namespace {

struct Pt {
  double x, y;
};

GSHGEOM_Path * polyline(const std::vector<std::vector<Pt>> & subs, const std::vector<bool> & closed) {
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

struct Out {
  GSHGEOM_Stroke_Outline * o = nullptr;
  ~Out() { gshgeom_stroke_outline_destroy(o); }
};

enum { OUT = 0, IN = 1, ON = 2 };

int wind(const GSHGEOM_Flat_Path * f, GSHGEOM_Fill_Rule rule, double x, double y) {
  if (!f) {
    return OUT;
  }
  int64_t w = 0;
  int on = 0;
  GSHGEOM_Device_Point p{(int32_t)std::lround(x * 64.0), (int32_t)std::lround(y * 64.0)};
  EXPECT_EQ(gshgeom_flat_path_winding(f, p, &w, &on), GSHGEOM_OK);
  if (on) {
    return ON;
  }
  return (rule == GSHGEOM_FILL_NONZERO ? w != 0 : (w & 1) != 0) ? IN : OUT;
}

/** Whether the outline of a stroke covers the device point: the body, or the aligned part within (or outside) the clip (SG12). */
int covers(const GSHGEOM_Stroke_Outline * o, double x, double y) {
  const int b = wind(gshgeom_stroke_outline_body(o), GSHGEOM_FILL_NONZERO, x, y);
  if (b != OUT) {
    return b;
  }
  const GSHGEOM_Flat_Path * al = gshgeom_stroke_outline_aligned(o);
  if (!al) {
    return OUT;
  }
  const int a = wind(al, GSHGEOM_FILL_NONZERO, x, y);
  const int c = wind(gshgeom_stroke_outline_clip(o), gshgeom_stroke_outline_clip_rule(o), x, y);
  if (a == ON || c == ON) {
    return ON;
  }
  if (a == OUT) {
    return OUT;
  }
  return (c == IN) != (gshgeom_stroke_outline_clip_complement(o) != 0) ? IN : OUT;
}

/** The two outlines agree at every point of the grid, up to a quarter pixel of boundary; @return the mismatches that remain. */
int disagreements(const GSHGEOM_Stroke_Outline * a, const GSHGEOM_Stroke_Outline * b, double lo, double hi, double step, size_t * inside_count) {
  int bad = 0;
  for (double y = lo; y <= hi; y += step) {
    for (double x = lo; x <= hi; x += step) {
      const int ca = covers(a, x, y), cb = covers(b, x, y);
      if (ca == ON || cb == ON) {
        continue;
      }
      if (ca == IN && inside_count) {
        (*inside_count)++;
      }
      if (ca == cb) {
        continue;
      }
      // the one that says OUT must reach an IN within a quarter of a pixel: the two boundaries are one curve within the tolerance
      const GSHGEOM_Stroke_Outline * out_one = ca == OUT ? a : b;
      bool near = false;
      for (int k = 0; k < 8 && !near; k++) {
        const double ang = 0.7853981633974483 * k;
        const int c2 = covers(out_one, x + 0.25 * std::cos(ang), y + 0.25 * std::sin(ang));
        near = c2 == ON || c2 == IN;
      }
      if (!near) {
        if (bad < 3) {
          fprintf(stderr, "mismatch at (%.3f, %.3f): profile %d plain %d\n", x, y, covers(a, x, y), covers(b, x, y));
        }
        bad++;
      }
    }
  }
  return bad;
}

GSHGEOM_Map map_of(int which) {
  GSHGEOM_Map m;
  gshgeom_map_identity(2, &m);
  switch (which) {
  case 1:
    m.m[0][0] = 1.6;
    m.m[1][1] = 0.8;
    m.t[0] = 3;
    m.t[1] = 2;
    break;
  case 2:
    m.m[0][0] = 0.9;
    m.m[0][1] = -0.5;
    m.m[1][0] = 0.4;
    m.m[1][1] = 1.1;
    break;
  case 3:
    m.m[0][0] = -1.2;
    m.m[1][1] = 1.0;
    m.t[0] = 60;
    break;
  default:
    break;
  }
  return m;
}

struct Spec {
  std::vector<std::vector<Pt>> subs;
  std::vector<bool> closed;
};

Spec random_spec(std::mt19937 & rng) {
  std::uniform_real_distribution<double> u(6.0, 54.0);
  Spec s;
  const int nsub = 1 + (int)(rng() % 3);
  for (int i = 0; i < nsub; i++) {
    const int n = 2 + (int)(rng() % 4);
    std::vector<Pt> v;
    for (int k = 0; k < n; k++) {
      v.push_back({u(rng), u(rng)});
    }
    s.subs.push_back(v);
    s.closed.push_back(n >= 3 && rng() % 3 == 0);
  }
  return s;
}

GSHGEOM_Stroke random_style(std::mt19937 & rng, bool dashes, bool align) {
  GSHGEOM_Stroke s;
  gshgeom_stroke_init(&s);
  static const double w[] = {1.5, 3.0, 5.5, 9.0};
  s.width = w[rng() % 4];
  s.join = (uint32_t)(rng() % GSHGEOM_STROKE_JOIN_COUNT);
  s.cap = (uint32_t)(rng() % GSHGEOM_STROKE_CAP_COUNT);
  s.miter_limit = 1.0 + (double)(rng() % 70) / 10.0;
  s.width_space = rng() % 3 == 0 ? GSHGEOM_STROKE_WIDTH_DEVICE : GSHGEOM_STROKE_WIDTH_LOCAL;
  if (dashes) {
    static const double pats[][4] = {{6, 4, 0, 0}, {10, 3, 2, 3}, {0, 5, 0, 0}, {14, 6, 0, 0}, {3, 2, 7, 4}};
    static const size_t counts[] = {2, 4, 2, 2, 4};
    const size_t k = rng() % 5;
    s.dashes = pats[k];
    s.dash_count = counts[k];
    s.dash_offset = (double)(rng() % 120) / 10.0 - 3.0;
    s.dash_scope = rng() % 2 ? GSHGEOM_STROKE_DASH_PATH : GSHGEOM_STROKE_DASH_PER_SUBPATH;
  }
  if (align) {
    s.alignment = rng() % 2 ? GSHGEOM_STROKE_ALIGN_INSIDE : GSHGEOM_STROKE_ALIGN_OUTSIDE;
    s.alignment_rule = rng() % 2 ? GSHGEOM_FILL_NONZERO : GSHGEOM_FILL_EVEN_ODD;
  }
  return s;
}

/**
 * Whether a dash boundary falls within 0.1 pixel of a vertex.  The library measures the pattern along the flattened centerline, whose
 * points are rounded to 1/64 pixel, and the plain stroke measures the curve itself: a boundary that close to a vertex is a join in one
 * and not in the other, which is the rounding of the measure and not a defect.  Such a case is not compared.
 */
bool boundary_near_vertex(const Spec & sp, const GSHGEOM_Stroke & st, const GSHGEOM_Map & map) {
  if (!st.dash_count) {
    return false;
  }
  double begin = 0.0;
  for (size_t i = 0; i < sp.subs.size(); i++) {
    const std::vector<Pt> & v = sp.subs[i];
    const bool cl = i < sp.closed.size() && sp.closed[i];
    std::vector<double> verts{0.0};
    double d = 0.0;
    const size_t segs = cl ? v.size() : v.size() - 1;
    for (size_t k = 0; k < segs; k++) {
      const Pt a = v[k], b = v[(k + 1) % v.size()];
      double dx = b.x - a.x, dy = b.y - a.y;
      if (st.width_space == GSHGEOM_STROKE_WIDTH_DEVICE) {
        const double mx = map.m[0][0] * dx + map.m[0][1] * dy, my = map.m[1][0] * dx + map.m[1][1] * dy;
        dx = mx;
        dy = my;
      }
      d += std::hypot(dx, dy);
      verts.push_back(d);
    }
    std::vector<double> bounds;
    auto fn = [](void * c, double a, double b) -> GSHGEOM_Result {
      auto * out = static_cast<std::vector<double> *>(c);
      out->push_back(a);
      out->push_back(b);
      return GSHGEOM_OK;
    };
    const double off = st.dash_offset + (st.dash_scope == GSHGEOM_STROKE_DASH_PATH ? begin : 0.0);
    gshgeom_stroke_dash_each(st.dashes, st.dash_count, off, d, 1u << 20, fn, &bounds);
    for (double b : bounds) {
      if (b < 1e-9 || b > d - 1e-9) {
        continue; // the subpath's own ends are exact, not measured
      }
      for (double vd : verts) {
        if (std::fabs(b - vd) < 0.1) {
          return true;
        }
      }
    }
    begin += d;
  }
  return false;
}

/** The profile stroke and the plain stroke of one style, compared. */
void differential(const Spec & sp, GSHGEOM_Stroke style, const GSHGEOM_Map & map, const char * what) {
  if (boundary_near_vertex(sp, style, map)) {
    return;
  }
  const GSHGEOM_Width_Stop one[1] = {{0, 1}};
  GSHGEOM_Width_Profile wp{1, one};
  GSHGEOM_Stroke plain = style;
  GSHGEOM_Stroke prof = style;
  prof.width_profile = &wp;
  GSHGEOM_Path * path = polyline(sp.subs, sp.closed);
  Out a, b;
  const GSHGEOM_Result ra = gshgeom_stroke_expand(path, &prof, &map, nullptr, nullptr, &a.o);
  const GSHGEOM_Result rb = gshgeom_stroke_expand(path, &plain, &map, nullptr, nullptr, &b.o);
  gshgeom_path_destroy(path);
  ASSERT_EQ(rb, GSHGEOM_OK) << what;
  ASSERT_EQ(ra, GSHGEOM_OK) << what << " join " << style.join << " cap " << style.cap;
  ASSERT_EQ(gshgeom_stroke_outline_status(a.o), gshgeom_stroke_outline_status(b.o)) << what;
  if (gshgeom_stroke_outline_status(a.o) != GSHGEOM_STROKE_DRAWN) {
    return;
  }
  size_t inside = 0;
  const int bad = disagreements(a.o, b.o, -2, 66, 0.37, &inside);
  if (bad) {
    for (size_t i = 0; i < sp.subs.size(); i++) {
      fprintf(stderr, "sub %zu closed %d:", i, i < sp.closed.size() ? (int)sp.closed[i] : 0);
      for (const Pt & q : sp.subs[i]) {
        fprintf(stderr, " (%.4f,%.4f)", q.x, q.y);
      }
      fprintf(stderr, "\n");
    }
    fprintf(stderr, "map %.3f %.3f %.3f %.3f t %.2f %.2f\n", map.m[0][0], map.m[0][1], map.m[1][0], map.m[1][1], map.t[0], map.t[1]);
  }
  EXPECT_EQ(bad, 0) << what << ": join " << style.join << " cap " << style.cap << " miter " << style.miter_limit << " width " << style.width << " space "
                    << style.width_space << " dashes " << style.dash_count << " align " << style.alignment << " pat " << (style.dash_count ? style.dashes[0] : 0) << "," << (style.dash_count ? style.dashes[1] : 0) << " off " << style.dash_offset << " scope " << style.dash_scope << " closed0 " << (sp.closed.empty() ? 0 : (int)sp.closed[0]) << " n0 " << sp.subs[0].size() << " nsub " << sp.subs.size();
}

}  // namespace

TEST(StrokeProfileForms, AConstantProfileIsTheConstantStrokeForEveryJoinAndCap) {
  std::mt19937 rng(7001);
  for (int it = 0; it < 300; it++) {
    const Spec sp = random_spec(rng);
    differential(sp, random_style(rng, false, false), map_of((int)(rng() % 4)), "solid");
  }
}

TEST(StrokeProfileForms, AConstantProfileIsTheConstantStrokeUnderDashes) {
  std::mt19937 rng(7002);
  for (int it = 0; it < 300; it++) {
    const Spec sp = random_spec(rng);
    differential(sp, random_style(rng, true, false), map_of((int)(rng() % 4)), "dashed");
  }
}

TEST(StrokeProfileForms, AConstantProfileIsTheConstantStrokeUnderAnAlignment) {
  std::mt19937 rng(7003);
  for (int it = 0; it < 300; it++) {
    Spec sp = random_spec(rng);
    // make sure a closed subpath is there
    if (it % 2 == 0) {
      sp.subs.push_back({{12, 12}, {50, 14}, {44, 50}, {16, 40}});
      sp.closed.push_back(true);
    }
    differential(sp, random_style(rng, it % 3 == 0, true), map_of((int)(rng() % 4)), "aligned");
  }
}

TEST(StrokeProfileForms, ADegenerateSubpathAndADotFollowTheConstantStroke) {
  std::mt19937 rng(7004);
  for (int it = 0; it < 120; it++) {
    Spec sp;
    sp.subs = {{{20, 20}, {20, 20}}, {{40, 30}, {40, 30}}, {{10, 50}, {50, 50}}};
    sp.closed = {false, true, false};
    GSHGEOM_Stroke st = random_style(rng, it % 2 == 0, false);
    st.degenerate = rng() % 2 ? GSHGEOM_STROKE_DEGENERATE_DRAW : GSHGEOM_STROKE_DEGENERATE_OMIT;
    differential(sp, st, map_of((int)(rng() % 4)), "degenerate");
  }
}

TEST(StrokeProfileForms, ATaperWithFlatEndsIsTheTrapezoid) {
  // a line from (10, 20) to (50, 20), width 8, factor 1 at 0 down to 0.25 at 40: the radius at x is 4 (1 - 0.75 (x - 10) / 40)
  const GSHGEOM_Width_Stop st[2] = {{0, 1}, {40, 0.25}};
  GSHGEOM_Width_Profile wp{2, st};
  GSHGEOM_Stroke s;
  gshgeom_stroke_init(&s);
  s.width = 8;
  s.cap = GSHGEOM_STROKE_CAP_BUTT;
  s.join = GSHGEOM_STROKE_JOIN_MITER;
  s.width_profile = &wp;
  GSHGEOM_Map m = map_of(0);
  GSHGEOM_Path * path = polyline({{{10, 20}, {50, 20}}}, {});
  Out o;
  ASSERT_EQ(gshgeom_stroke_expand(path, &s, &m, nullptr, nullptr, &o.o), GSHGEOM_OK);
  gshgeom_path_destroy(path);
  auto r = [](double x) { return 4.0 * (1.0 - 0.75 * (x - 10.0) / 40.0); };
  for (double x : {10.5, 20.0, 30.0, 49.5}) {
    EXPECT_EQ(covers(o.o, x, 20.0 + r(x) - 0.12), IN) << x;
    EXPECT_EQ(covers(o.o, x, 20.0 - r(x) + 0.12), IN) << x;
    EXPECT_EQ(covers(o.o, x, 20.0 + r(x) + 0.12), OUT) << x;
    EXPECT_EQ(covers(o.o, x, 20.0 - r(x) - 0.12), OUT) << x;
  }
  EXPECT_EQ(covers(o.o, 9.8, 20.0), OUT) << "a butt cap is flat";
  EXPECT_EQ(covers(o.o, 50.2, 20.0), OUT);
  EXPECT_EQ(covers(o.o, 10.2, 23.6), IN) << "the full width at the start, to the corner";
  // a square cap reaches the radius beyond the end, a round one a half disc
  s.cap = GSHGEOM_STROKE_CAP_SQUARE;
  Out q;
  path = polyline({{{10, 20}, {50, 20}}}, {});
  ASSERT_EQ(gshgeom_stroke_expand(path, &s, &m, nullptr, nullptr, &q.o), GSHGEOM_OK);
  gshgeom_path_destroy(path);
  EXPECT_EQ(covers(q.o, 10.0 - 3.9, 23.8), IN);
  EXPECT_EQ(covers(q.o, 10.0 - 4.2, 20.0), OUT);
  EXPECT_EQ(covers(q.o, 50.0 + 0.9, 20.0 + 0.9), IN) << "the end radius is 1";
  EXPECT_EQ(covers(q.o, 50.0 + 1.2, 20.0), OUT);
  s.cap = GSHGEOM_STROKE_CAP_ROUND;
  Out rd;
  path = polyline({{{10, 20}, {50, 20}}}, {});
  ASSERT_EQ(gshgeom_stroke_expand(path, &s, &m, nullptr, nullptr, &rd.o), GSHGEOM_OK);
  gshgeom_path_destroy(path);
  EXPECT_EQ(covers(rd.o, 10.0 - 3.9, 20.0), IN);
  EXPECT_EQ(covers(rd.o, 10.0 - 3.0, 22.8), OUT) << "a half disc, not a square";
}

TEST(StrokeProfileForms, AJumpWithFlatEndsIsAStep) {
  // two stops at 20 (the line starts at x = 10, so the jump is at x = 30): the width is 2 before it and 8 after it
  const GSHGEOM_Width_Stop st[3] = {{0, 0.25}, {20, 0.25}, {20, 1}};
  GSHGEOM_Width_Profile wp{3, st};
  GSHGEOM_Stroke s;
  gshgeom_stroke_init(&s);
  s.width = 8;
  s.cap = GSHGEOM_STROKE_CAP_BUTT;
  s.width_profile = &wp;
  GSHGEOM_Map m = map_of(0);
  GSHGEOM_Path * path = polyline({{{10, 20}, {50, 20}}}, {});
  Out o;
  ASSERT_EQ(gshgeom_stroke_expand(path, &s, &m, nullptr, nullptr, &o.o), GSHGEOM_OK);
  gshgeom_path_destroy(path);
  EXPECT_EQ(covers(o.o, 15.0, 20.0 + 0.9), IN);
  EXPECT_EQ(covers(o.o, 15.0, 20.0 + 1.2), OUT);
  EXPECT_EQ(covers(o.o, 35.0, 20.0 + 3.8), IN);
  EXPECT_EQ(covers(o.o, 35.0, 20.0 + 4.2), OUT);
  EXPECT_EQ(covers(o.o, 29.9, 20.0 + 3.0), OUT);
  EXPECT_EQ(covers(o.o, 30.1, 20.0 + 3.0), IN);
}

TEST(StrokeProfileForms, AProfileAndADashEachMeasureFromTheSubpathStart) {
  // dashes 10 on, 10 off along a line of 60, the width tapering from 8 to 0 over 60: each dash has the radius of its own place
  const GSHGEOM_Width_Stop st[2] = {{0, 1}, {60, 0}};
  GSHGEOM_Width_Profile wp{2, st};
  const double dashes[2] = {10, 10};
  GSHGEOM_Stroke s;
  gshgeom_stroke_init(&s);
  s.width = 8;
  s.cap = GSHGEOM_STROKE_CAP_BUTT;
  s.dashes = dashes;
  s.dash_count = 2;
  s.width_profile = &wp;
  GSHGEOM_Map m = map_of(0);
  GSHGEOM_Path * path = polyline({{{0, 30}, {60, 30}}}, {});
  Out o;
  ASSERT_EQ(gshgeom_stroke_expand(path, &s, &m, nullptr, nullptr, &o.o), GSHGEOM_OK);
  gshgeom_path_destroy(path);
  auto r = [](double x) { return 4.0 * (1.0 - x / 60.0); };
  for (double x : {2.0, 8.0, 22.0, 28.0, 42.0, 48.0}) {
    EXPECT_EQ(covers(o.o, x, 30.0 + r(x) - 0.15), IN) << x;
    EXPECT_EQ(covers(o.o, x, 30.0 + r(x) + 0.15), OUT) << x;
  }
  for (double x : {12.0, 18.0, 32.0, 38.0, 52.0, 58.0}) {
    EXPECT_EQ(covers(o.o, x, 30.0), OUT) << "a gap at " << x;
  }
}

TEST(StrokeProfileForms, EveryFailingAllocationIsAnErrorAndLeaksNothing) {
  const GSHGEOM_Width_Stop st[3] = {{0, 1}, {30, 0.3}, {60, 1.5}};
  GSHGEOM_Width_Profile wp{3, st};
  const double dashes[2] = {9, 5};
  for (int variant = 0; variant < 4; variant++) {
    GSHGEOM_Stroke s;
    gshgeom_stroke_init(&s);
    s.width = 6;
    s.width_profile = &wp;
    s.join = variant == 0 ? GSHGEOM_STROKE_JOIN_MITER : GSHGEOM_STROKE_JOIN_ROUND;
    s.cap = variant == 1 ? GSHGEOM_STROKE_CAP_SQUARE : (variant == 2 ? GSHGEOM_STROKE_CAP_BUTT : GSHGEOM_STROKE_CAP_ROUND);
    if (variant >= 2) {
      s.dashes = dashes;
      s.dash_count = 2;
    }
    if (variant == 3) {
      s.alignment = GSHGEOM_STROKE_ALIGN_INSIDE;
    }
    GSHGEOM_Map m = map_of(2);
    GSHGEOM_Path * path = polyline({{{10, 10}, {50, 14}, {44, 50}, {16, 40}}}, {true});
    tracker::Tracker t;
    GSHGEOM_Allocator a = tracker::allocator_for(&t);
    bool done = false;
    for (long k = 0; k < 4000 && !done; k++) {
      t.persistent = false;
      t.failed = false;
      t.fail_from = t.calls + k;
      GSHGEOM_Stroke_Outline * o = nullptr;
      const GSHGEOM_Result r = gshgeom_stroke_expand(path, &s, &m, &a, nullptr, &o);
      t.fail_from = -1;
      if (r == GSHGEOM_OK) {
        done = true;
        EXPECT_FALSE(t.failed);
        gshgeom_stroke_outline_destroy(o);
      } else {
        EXPECT_EQ(r, GSHGEOM_ERR_OOM) << "variant " << variant << " allocation " << k;
        EXPECT_EQ(o, nullptr);
      }
      EXPECT_TRUE(t.live.empty()) << "variant " << variant << " allocation " << k << ": " << t.live.size() << " blocks leaked";
    }
    EXPECT_TRUE(done) << "variant " << variant;
    gshgeom_path_destroy(path);
  }
}

TEST(StrokeProfileForms, WorkAndPointLimitsAreErrorsNotTruncations) {
  const GSHGEOM_Width_Stop st[2] = {{0, 1}, {200, 3}};
  GSHGEOM_Width_Profile wp{2, st};
  GSHGEOM_Stroke s;
  gshgeom_stroke_init(&s);
  s.width = 4;
  s.join = GSHGEOM_STROKE_JOIN_MITER;
  s.cap = GSHGEOM_STROKE_CAP_SQUARE;
  s.width_profile = &wp;
  GSHGEOM_Map m = map_of(0);
  GSHGEOM_Path * path = polyline({{{10, 10}, {50, 14}, {44, 50}, {16, 40}, {30, 5}, {58, 58}}}, {});
  GSHGEOM_Limits lim;
  gshgeom_limits_default(&lim);
  lim.max_work = 8;
  GSHGEOM_Stroke_Outline * o = nullptr;
  EXPECT_EQ(gshgeom_stroke_expand(path, &s, &m, nullptr, &lim, &o), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(o, nullptr);
  const double dashes[2] = {0.5, 0.5};
  s.dashes = dashes;
  s.dash_count = 2;
  gshgeom_limits_default(&lim);
  lim.max_dash_steps = 5;
  EXPECT_EQ(gshgeom_stroke_expand(path, &s, &m, nullptr, &lim, &o), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(o, nullptr);
  gshgeom_path_destroy(path);
}

namespace {

/** A hatch of two families (a solid one at 45 degrees and a dashed one at 100) over a box, and the map that carries it. */
void hatch_outline(const GSHGEOM_Stroke & style, Out * out, GSHGEOM_Result * result) {
  static const double dashes[2] = {5, -3};
  GSHGEOM_Hatch_Family fam[2];
  std::memset(fam, 0, sizeof fam);
  fam[0].angle = 45;
  fam[0].offset_across = 7;
  fam[1].angle = 100;
  fam[1].origin_x = 3;
  fam[1].offset_across = 11;
  fam[1].dash_count = 2;
  fam[1].dashes = dashes;
  GSHGEOM_Hatch h;
  std::memset(&h, 0, sizeof h);
  h.angle_unit = GSHGEOM_ANGLE_DEGREES;
  h.nesting = GSHGEOM_HATCH_NESTING_RULE;
  gshgeom_map_identity(2, &h.map);
  h.family_count = 2;
  h.families = fam;
  GSHGEOM_Path * region = polyline({{{8, 8}, {52, 10}, {50, 50}, {12, 54}}}, {true});
  GSHGEOM_Map chain = map_of(2);
  *result = gshgeom_hatch_expand(region, GSHGEOM_FILL_NONZERO, &h, &style, &chain, 100000, nullptr, nullptr, &out->o);
  gshgeom_path_destroy(region);
}

}  // namespace

TEST(StrokeProfileForms, AConstantProfileOnAHatchPenIsTheHatchOfThePlainPen) {
  std::mt19937 rng(7006);
  const GSHGEOM_Width_Stop one[1] = {{0, 1}};
  GSHGEOM_Width_Profile wp{1, one};
  for (int it = 0; it < 24; it++) {
    GSHGEOM_Stroke plain = random_style(rng, false, false);
    plain.width = 1.0 + (double)(rng() % 20) / 4.0;
    plain.degenerate = GSHGEOM_STROKE_DEGENERATE_DRAW;
    GSHGEOM_Stroke prof = plain;
    prof.width_profile = &wp;
    Out a, b;
    GSHGEOM_Result ra, rb;
    hatch_outline(prof, &a, &ra);
    hatch_outline(plain, &b, &rb);
    ASSERT_EQ(rb, GSHGEOM_OK);
    ASSERT_EQ(ra, GSHGEOM_OK) << "join " << plain.join << " cap " << plain.cap;
    EXPECT_EQ(disagreements(a.o, b.o, -1.9269, 66, 0.3713, nullptr), 0) << "join " << plain.join << " cap " << plain.cap << " width " << plain.width;
  }
}

TEST(StrokeProfileForms, AHatchPenThatTapersTapersEachPieceFromItsOwnStart) {
  // a single family of horizontal lines: each piece starts at the region's left edge and the profile falls from 1 to 0.2 over 40 units
  const GSHGEOM_Width_Stop st[2] = {{0, 1}, {40, 0.2}};
  GSHGEOM_Width_Profile wp{2, st};
  GSHGEOM_Stroke s;
  gshgeom_stroke_init(&s);
  s.width = 6;
  s.cap = GSHGEOM_STROKE_CAP_BUTT;
  s.join = GSHGEOM_STROKE_JOIN_MITER;
  s.width_profile = &wp;
  GSHGEOM_Hatch_Family fam;
  std::memset(&fam, 0, sizeof fam);
  fam.angle = 0;
  fam.origin_y = 20;
  fam.offset_across = 100; // one line, at y = 20
  GSHGEOM_Hatch h;
  std::memset(&h, 0, sizeof h);
  h.angle_unit = GSHGEOM_ANGLE_DEGREES;
  h.nesting = GSHGEOM_HATCH_NESTING_RULE;
  gshgeom_map_identity(2, &h.map);
  h.family_count = 1;
  h.families = &fam;
  GSHGEOM_Path * region = polyline({{{10, 5}, {60, 5}, {60, 35}, {10, 35}}}, {true});
  GSHGEOM_Map chain = map_of(0);
  Out o;
  ASSERT_EQ(gshgeom_hatch_expand(region, GSHGEOM_FILL_NONZERO, &h, &s, &chain, 1000, nullptr, nullptr, &o.o), GSHGEOM_OK);
  gshgeom_path_destroy(region);
  // the one piece runs along y = 20 from x = 10 to 60: the radius at x is 3 (1 - 0.8 (x - 10) / 40), constant past x = 50
  auto r = [](double x) { return 3.0 * (1.0 - 0.8 * std::min(x - 10.0, 40.0) / 40.0); };
  for (double x : {11.0, 25.0, 40.0, 58.0}) {
    EXPECT_EQ(covers(o.o, x, 20.0 + r(x) - 0.12), IN) << x;
    EXPECT_EQ(covers(o.o, x, 20.0 - r(x) + 0.12), IN) << x;
    EXPECT_EQ(covers(o.o, x, 20.0 + r(x) + 0.12), OUT) << x;
    EXPECT_EQ(covers(o.o, x, 20.0 - r(x) - 0.12), OUT) << x;
  }
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
