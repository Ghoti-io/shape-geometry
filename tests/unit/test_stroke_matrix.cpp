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
 * The expansion against the geometric oracle over the stroke's case matrix
 * (stroke.h, SG17): lines and every join, the miter boundary, joins against caps,
 * every curve kind (cusps, loops, circles with h below, at and above the
 * radius), zero extent, degenerate subpaths, lone moves, both width spaces under
 * every kind of map, alignment, overlaps, and every tolerance.  For each the
 * outline's nonzero fill agrees with the region at every decided point, no
 * subpath winds negatively, every vertex is inside the conservative box plus the
 * tolerance, and the deviation profile's maximum is recorded as a test property
 * and is at most the tolerance; and the exact vectors of tools/gen-stroke-vectors.py
 * are classified as the exact model classifies them.
 */

#include "../test_helpers.h"

#include "../stroke_matrix.h"
#include "../stroke_vectors.h"

#include <string>
#include <vector>

namespace sm = stroke_matrix;
namespace so = stroke_oracle;

namespace {

typedef sm::Ops Ops;

GSHGEOM_Map map_id() {
  GSHGEOM_Map m;
  gshgeom_map_identity(2, &m);
  return m;
}

GSHGEOM_Map map_of(double a, double b, double c, double d, double tx, double ty) {
  double mm[4] = {a, b, c, d}, t[2] = {tx, ty};
  GSHGEOM_Map m;
  EXPECT_EQ(gshgeom_map_from_matrix(2, mm, t, &m), GSHGEOM_OK);
  return m;
}

GSHGEOM_Stroke style(double width, int join, double limit, int cap) {
  GSHGEOM_Stroke s;
  gshgeom_stroke_init(&s);
  s.width = width;
  s.join = (uint32_t)join;
  s.miter_limit = limit;
  s.cap = (uint32_t)cap;
  return s;
}

struct Case {
  std::string name;
  Ops ops;
  GSHGEOM_Stroke st;
  GSHGEOM_Map map;
  int32_t tau;
  std::vector<double> dashes = {}; ///< the dash array (SG19), empty for a solid stroke
};

/** The worst deviation seen over every case, in 26.6 units. */
double g_worst = 0.0;
double g_worst_margin_tau = 0.0;

void run(const Case & c) {
  SCOPED_TRACE(c.name);
  GSHGEOM_Stroke st = c.st;
  if (!c.dashes.empty()) {
    st.dashes = c.dashes.data();
    st.dash_count = c.dashes.size();
  }
  GSHGEOM_Path * path = sm::build(c.ops);
  ASSERT_NE(path, nullptr);
  GSHGEOM_Stroke_Outline * o = nullptr;
  ASSERT_EQ(gshgeom_stroke_expand_with_tolerance(path, &st, &c.map, c.tau, nullptr, nullptr, &o), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_stroke_outline_status(o), GSHGEOM_STROKE_DRAWN);
  so::Centerline cl = so::read_centerline(path);
  so::Stroked sk = so::build_stroke(cl, sm::oracle_style(st), so::affine_of(c.map), (int)st.alignment);
  sm::Box box = sm::stroke_box(path, st, c.map);
  struct Part {
    const GSHGEOM_Flat_Path * fp;
    const so::Region * r;
    const char * name;
  } parts[2] = {{gshgeom_stroke_outline_body(o), &sk.body, "body"},
      {gshgeom_stroke_outline_aligned(o), &sk.aligned, "aligned"}};
  EXPECT_EQ(gshgeom_stroke_outline_aligned(o) != nullptr, sk.has_aligned);
  for (const Part & p : parts) {
    if (!p.fp) {
      continue;
    }
    sm::PartResult pr = sm::compare_part(p.fp, *p.r, box, c.tau, 2500, 0x51ULL);
    EXPECT_GT(pr.decided, 0) << p.name;
    EXPECT_EQ(pr.wrong, 0) << p.name << ": " << pr.first;
    EXPECT_EQ(pr.negative, 0) << p.name;
    EXPECT_EQ(pr.outside_box, 0) << p.name;
    double dev_units = (double)(pr.worst_outward * 64);
    EXPECT_LE(dev_units, (double)c.tau) << p.name;
    if (dev_units > g_worst) {
      g_worst = dev_units;
    }
    if (dev_units / c.tau > g_worst_margin_tau) {
      g_worst_margin_tau = dev_units / c.tau;
    }
  }
  gshgeom_stroke_outline_destroy(o);
  gshgeom_path_destroy(path);
}

const int MITER = GSHGEOM_STROKE_JOIN_MITER, ROUND = GSHGEOM_STROKE_JOIN_ROUND, BEVEL = GSHGEOM_STROKE_JOIN_BEVEL,
          MCLIP = GSHGEOM_STROKE_JOIN_MITER_CLIP, ARCS = GSHGEOM_STROKE_JOIN_ARCS;
const int BUTT = GSHGEOM_STROKE_CAP_BUTT, RCAP = GSHGEOM_STROKE_CAP_ROUND, SQUARE = GSHGEOM_STROKE_CAP_SQUARE;

Ops circle(double cx, double cy, double r) {
  double w = 0.70710678118654752440;
  return Ops{{'M', {cx + r, cy}}, {'R', {cx + r, cy + r, cx, cy + r, 1, w, 1}}, {'R', {cx - r, cy + r, cx - r, cy, 1, w, 1}},
      {'R', {cx - r, cy - r, cx, cy - r, 1, w, 1}}, {'R', {cx + r, cy - r, cx + r, cy, 1, w, 1}}, {'Z', {}}};
}

} // namespace

TEST(StrokeMatrix, LinesEveryJoinAndCap) {
  std::vector<Ops> paths = {
      {{'M', {0, 0}}, {'L', {10, 0}}, {'L', {10, 10}}},                 // right angle
      {{'M', {0, 0}}, {'L', {10, 0}}, {'L', {1, 3}}},                   // acute
      {{'M', {0, 0}}, {'L', {10, 0}}, {'L', {20, 4}}},                  // obtuse
      {{'M', {0, 0}}, {'L', {5, 0}}, {'L', {10, 0}}},                   // collinear, same direction (smooth)
      {{'M', {0, 0}}, {'L', {10, 0}}, {'L', {4, 0}}},                   // reversal
      {{'M', {0, 0}}, {'L', {6, 0}}, {'L', {6, 6}}, {'L', {0, 6}}, {'Z', {}}}, // closed
      {{'M', {0, 0}}, {'L', {8, 0}}, {'L', {8.2, 0.3}}, {'L', {0, 6}}},    // a very short segment
  };
  /* every join on every path; the caps in turn (every cap on the first path) */
  int turn = 0;
  for (const Ops & p : paths) {
    for (int join : {MITER, ROUND, BEVEL, MCLIP, ARCS}) {
      for (int cap : {BUTT, RCAP, SQUARE}) {
        if (&p != &paths[0] && cap != (turn % 3)) {
          continue;
        }
        Case c{"lines join " + std::to_string(join) + " cap " + std::to_string(cap), p, style(3, join, 4, cap),
            map_id(), GSHGEOM_FLATTEN_TOLERANCE};
        run(c);
      }
      turn++;
    }
  }
}

TEST(StrokeMatrix, TheMiterBoundaryIsDecidedExactly) {
  Ops tri = {{'M', {0, 0}}, {'L', {4, 3}}, {'L', {8, 0}}};
  for (double L : {1.25, std::nextafter(1.25, 0.0), 1.1, 3.0}) {
    for (int join : {MITER, MCLIP}) {
      run({"miter boundary", tri, style(8, join, L, BUTT), map_id(), GSHGEOM_FLATTEN_TOLERANCE});
    }
  }
}

TEST(StrokeMatrix, EveryCurveKind) {
  std::vector<std::pair<std::string, Ops>> curves = {
      {"quadratic", {{'M', {0, 0}}, {'Q', {5, 9, 10, 0}}}},
      {"cubic inflection", {{'M', {0, 0}}, {'C', {4, 8, 6, -8, 10, 0}}}},
      {"cubic loop", {{'M', {0, 0}}, {'C', {12, 8, -2, 8, 10, 0}}}},
      {"cubic cusp", {{'M', {0, 0}}, {'C', {8, 4, 4, 4, 4, 0}}}},
      {"cubic near cusp", {{'M', {0, 0}}, {'C', {8, 4, 4.01, 4, 4, 0}}}},
      {"cubic p0 is p1", {{'M', {0, 0}}, {'C', {0, 0, 6, 2, 8, 6}}}},
      {"rational cubic", {{'M', {0, 0}}, {'K', {3, 6, 7, 6, 10, 0, 1, 2, 0.5, 1.5}}}},
      {"extreme weights", {{'M', {0, 0}}, {'R', {5, 8, 10, 0, 1, 1e6, 1}}}},
      {"light middle", {{'M', {0, 0}}, {'R', {5, 8, 10, 0, 1, 1e-6, 1}}}},
      {"circle h below r", circle(0, 0, 6)},
      {"circle h is r", circle(0, 0, 2)},
      {"circle h above r", circle(0, 0, 1)},
      {"collinear back and forth", {{'M', {0, 0}}, {'C', {10, 0, -5, 0, 5, 0}}}},
  };
  const int pairs[3][2] = {{MITER, BUTT}, {ROUND, RCAP}, {MITER, SQUARE}};
  for (const auto & cv : curves) {
    for (const auto & jc : pairs) {
      run({cv.first, cv.second, style(4, jc[0], 4, jc[1]), map_id(), GSHGEOM_FLATTEN_TOLERANCE});
    }
  }
}

TEST(StrokeMatrix, ZeroExtentAndDegenerateSubpaths) {
  std::vector<Ops> paths = {
      {{'M', {2, 2}}, {'L', {2, 2}}},
      {{'M', {2, 2}}, {'Z', {}}},
      {{'M', {2, 2}}, {'C', {2, 2, 2, 2, 2, 2}}},
      {{'M', {0, 0}}, {'L', {6, 0}}, {'L', {6, 0}}, {'L', {6, 6}}},
      {{'M', {0, 0}}, {'L', {3, 4}}, {'M', {10, 10}}, {'Z', {}}},
      {{'M', {10, 10}}, {'L', {10, 10}}, {'M', {0, 0}}, {'L', {4, -3}}},
      {{'M', {0, 0}}, {'L', {5, 5}}, {'M', {3, 3}}, {'M', {10, 10}}, {'L', {10, 10}}},
      {{'M', {0, 0}}, {'L', {6, 0}}, {'L', {6, 6}}, {'Z', {}}, {'M', {20, 20}}, {'Z', {}}},
      {{'M', {4, 4}}},
  };
  std::vector<GSHGEOM_Map> maps = {map_id(), map_of(0.6, -0.8, 0.8, 0.6, 1, 2), map_of(1, 1, 0, 1, 0, 0)};
  for (const Ops & p : paths) {
    for (int cap : {BUTT, RCAP, SQUARE}) {
      for (uint32_t deg : {0u, 1u}) {
        for (uint32_t space : {0u, 1u}) {
          for (const GSHGEOM_Map & m : maps) {
            if ((cap != SQUARE || deg != 0) && (&m != &maps[0] || space != 0)) {
              continue; /* every cap and both settings under the identity; the square under every map and space */
            }
            GSHGEOM_Stroke s = style(3, MITER, 4, cap);
            s.degenerate = deg;
            s.width_space = space;
            GSHGEOM_Path * pp = sm::build(p);
            so::Centerline cl = so::read_centerline(pp);
            gshgeom_path_destroy(pp);
            bool draws = false;
            for (const so::Sub & sb : cl.subs) {
              draws = draws || (sb.draws && (!sb.segs.empty() || (deg == 0 && cap != BUTT)));
            }
            if (!draws) {
              /* nothing is drawn: the body is empty */
              GSHGEOM_Path * q = sm::build(p);
              GSHGEOM_Stroke_Outline * o = nullptr;
              ASSERT_EQ(gshgeom_stroke_expand(q, &s, &m, nullptr, nullptr, &o), GSHGEOM_OK);
              EXPECT_EQ(gshgeom_flat_path_subpath_count(gshgeom_stroke_outline_body(o)), 0u);
              gshgeom_stroke_outline_destroy(o);
              gshgeom_path_destroy(q);
              continue;
            }
            run({"degenerate", p, s, m, GSHGEOM_FLATTEN_TOLERANCE});
          }
        }
      }
    }
  }
}

TEST(StrokeMatrix, BothWidthSpacesUnderEveryMap) {
  Ops p = {{'M', {0, 0}}, {'L', {6, 0}}, {'Q', {9, 0, 9, 4}}, {'L', {2, 6}}};
  std::vector<std::pair<std::string, GSHGEOM_Map>> maps = {
      {"identity", map_id()},
      {"uniform", map_of(2, 0, 0, 2, 3, -1)},
      {"scale 2 1", map_of(2, 0, 0, 1, 0, 0)},
      {"rotation", map_of(0.6, -0.8, 0.8, 0.6, 5, 5)},
      {"skew", map_of(1, 0.75, 0, 1, 0, 0)},
      {"mirror", map_of(-1, 0, 0, 1, 20, 0)},
      {"tiny", map_of(0x1p-40, 0, 0, 0x1p-40, 0, 0)},
  };
  for (const auto & mp : maps) {
    for (uint32_t space : {0u, 1u}) {
      for (int join : {MITER, ROUND, BEVEL}) {
        GSHGEOM_Stroke s = style(space ? 2 : 1.5, join, 4, RCAP);
        s.width_space = space;
        Ops q = p;
        if (mp.first == "tiny" && space == 0) {
          /* model coordinates large enough that the device stroke is pixels wide */
          for (auto & op : q) {
            for (double & v : op.second) {
              v *= 0x1p40;
            }
          }
          s.width = 1.5 * 0x1p40;
        } else if (mp.first == "tiny") {
          for (auto & op : q) {
            for (double & v : op.second) {
              v *= 0x1p40;
            }
          }
        }
        run({mp.first, q, s, mp.second, GSHGEOM_FLATTEN_TOLERANCE});
      }
    }
  }
  /* the hairline: device width 1 */
  GSHGEOM_Stroke hair = style(1, MITER, 4, BUTT);
  hair.width_space = GSHGEOM_STROKE_WIDTH_DEVICE;
  run({"hairline", p, hair, map_of(3, 0, 0, 3, 0, 0), GSHGEOM_FLATTEN_TOLERANCE});
}

TEST(StrokeMatrix, Alignment) {
  std::vector<std::pair<std::string, Ops>> paths = {
      {"square", {{'M', {0, 0}}, {'L', {8, 0}}, {'L', {8, 8}}, {'L', {0, 8}}, {'Z', {}}}},
      {"star", {{'M', {5, 0}}, {'L', {7, 4}}, {'L', {11, 4}}, {'L', {8, 7}}, {'L', {9, 11}}, {'L', {5, 8}}, {'L', {1, 11}},
                   {'L', {2, 7}}, {'L', {-1, 4}}, {'L', {3, 4}}, {'Z', {}}}},
      {"self intersecting", {{'M', {0, 0}}, {'L', {10, 10}}, {'L', {10, 0}}, {'L', {0, 10}}, {'Z', {}}}},
      {"open and closed", {{'M', {0, 0}}, {'L', {8, 0}}, {'L', {8, 8}}, {'Z', {}}, {'M', {12, 0}}, {'L', {16, 4}}}},
  };
  for (const auto & p : paths) {
    for (uint32_t align : {1u, 2u}) {
      for (int rule : {GSHGEOM_FILL_NONZERO, GSHGEOM_FILL_EVEN_ODD}) {
        GSHGEOM_Stroke s = style(2, MITER, 4, BUTT);
        s.alignment = align;
        s.alignment_rule = (GSHGEOM_Fill_Rule)rule;
        run({p.first, p.second, s, map_id(), GSHGEOM_FLATTEN_TOLERANCE});
      }
    }
  }
}

TEST(StrokeMatrix, OverlapsAreOneUnion) {
  std::vector<std::pair<std::string, Ops>> paths = {
      {"figure eight", {{'M', {0, 0}}, {'C', {10, 10, 10, -10, 0, 0}}, {'C', {-10, 10, -10, -10, 0, 0}}, {'Z', {}}}},
      {"retrace", {{'M', {0, 0}}, {'L', {10, 0}}, {'L', {0, 0}}, {'L', {10, 0}}}},
      {"overlapping subpaths", {{'M', {0, 0}}, {'L', {10, 10}}, {'M', {0, 10}}, {'Q', {5, -6, 10, 0}}}},
      {"loop", {{'M', {0, 0}}, {'C', {14, 10, -4, 10, 10, 0}}}},
      {"coincident", {{'M', {0, 0}}, {'L', {8, 0}}, {'M', {0, 0}}, {'L', {8, 0}}, {'M', {8, 0}}, {'L', {0, 0}}}},
  };
  for (const auto & p : paths) {
    for (int join : {MITER, ROUND, BEVEL}) {
      run({p.first, p.second, style(3, join, 4, RCAP), map_id(), GSHGEOM_FLATTEN_TOLERANCE});
    }
  }
}

TEST(StrokeMatrix, EveryTolerance) {
  Ops p = {{'M', {0, 0}}, {'C', {8, 12, 14, -6, 20, 4}}, {'L', {12, 10}}};
  for (int32_t tau : {2, 3, 4, 16, 64}) {
    for (int join : {MITER, ROUND, ARCS}) {
      run({"tolerance " + std::to_string(tau), p, style(5, join, 4, RCAP), map_id(), tau});
    }
  }
}

TEST(StrokeMatrix, ArcsJoins) {
  /* the exact vectors' arcs cases, through the expansion, at every tolerance the matrix uses */
  for (size_t i = 0; i < stroke_vectors::case_count; i++) {
    const stroke_vectors::Case & vc = stroke_vectors::cases[i];
    if (vc.join != ARCS) {
      continue;
    }
    Ops ops;
    for (size_t k = 0; k < vc.op_count; k++) {
      ops.push_back({vc.ops[k].verb, std::vector<double>(vc.ops[k].v, vc.ops[k].v + vc.ops[k].n)});
    }
    GSHGEOM_Stroke s = style(vc.width, vc.join, vc.limit, vc.cap);
    s.width_space = (uint32_t)vc.space;
    for (int32_t tau : {2, 4}) {
      run({vc.name, ops, s, map_of(vc.map[0], vc.map[1], vc.map[2], vc.map[3], vc.map[4], vc.map[5]), tau});
    }
  }
}

TEST(StrokeMatrix, TheExactVectorsThroughTheExpansion) {
  /* every exact verdict, through the outline: the points are 0.2 px from every boundary the exact model drew, beyond
   * the default tolerance's reach */
  size_t checked = 0;
  for (size_t i = 0; i < stroke_vectors::case_count; i++) {
    const stroke_vectors::Case & vc = stroke_vectors::cases[i];
    Ops ops;
    for (size_t k = 0; k < vc.op_count; k++) {
      ops.push_back({vc.ops[k].verb, std::vector<double>(vc.ops[k].v, vc.ops[k].v + vc.ops[k].n)});
    }
    GSHGEOM_Stroke s = style(vc.width, vc.join, vc.limit, vc.cap);
    s.width_space = (uint32_t)vc.space;
    s.degenerate = (uint32_t)vc.degenerate;
    s.alignment = (uint32_t)vc.alignment;
    GSHGEOM_Map m = map_of(vc.map[0], vc.map[1], vc.map[2], vc.map[3], vc.map[4], vc.map[5]);
    GSHGEOM_Path * path = sm::build(ops);
    ASSERT_NE(path, nullptr);
    GSHGEOM_Stroke_Outline * o = nullptr;
    ASSERT_EQ(gshgeom_stroke_expand(path, &s, &m, nullptr, nullptr, &o), GSHGEOM_OK) << vc.name;
    so::Winding body, aligned;
    body.build(so::polygons_of(gshgeom_stroke_outline_body(o)));
    aligned.build(so::polygons_of(gshgeom_stroke_outline_aligned(o)));
    for (size_t k = 0; k < vc.point_count; k++) {
      const stroke_vectors::Point & pt = vc.points[k];
      int w = (pt.part == 0 ? body : aligned).at({pt.x, pt.y});
      EXPECT_EQ(w != 0, pt.inside != 0) << vc.name << " point " << k << " (" << pt.x << ", " << pt.y << ")";
      checked++;
    }
    gshgeom_stroke_outline_destroy(o);
    gshgeom_path_destroy(path);
  }
  EXPECT_GT(checked, (size_t)1000);
}

namespace {

Case dashed(const std::string & name, const Ops & ops, GSHGEOM_Stroke st, std::vector<double> d, double offset,
    uint32_t scope = 0, GSHGEOM_Map m = map_id()) {
  st.dash_offset = offset;
  st.dash_scope = scope;
  Case c{name, ops, st, m, GSHGEOM_FLATTEN_TOLERANCE};
  c.dashes = d;
  return c;
}

} // namespace

TEST(StrokeMatrix, DashesOnLines) {
  std::vector<std::pair<std::string, Ops>> paths = {
      {"line", {{'M', {0, 0}}, {'L', {40, 0}}}},
      {"slanted", {{'M', {0, 0}}, {'L', {24, 32}}}},
      {"corners", {{'M', {0, 0}}, {'L', {20, 0}}, {'L', {20, 15}}, {'L', {5, 18}}}},
      {"closed", {{'M', {0, 0}}, {'L', {16, 0}}, {'L', {16, 12}}, {'L', {0, 12}}, {'Z', {}}}},
      {"smooth and reversal", {{'M', {0, 0}}, {'L', {10, 0}}, {'L', {20, 0}}, {'L', {6, 0}}}},
      {"two subpaths", {{'M', {0, 0}}, {'L', {14, 0}}, {'M', {0, 8}}, {'L', {9, 8}}, {'L', {9, 20}}}},
  };
  std::vector<std::vector<double>> arrays = {{5, 3}, {4, 2, 1}, {0, 6}, {6, 0, 2, 3}, {20, 4}};
  int turn = 0;
  for (const auto & p : paths) {
    for (const auto & d : arrays) {
      for (int cap : {BUTT, RCAP, SQUARE}) {
        if ((turn++ % 2) && cap == BUTT) {
          continue;
        }
        int join = (turn % 5);
        GSHGEOM_Stroke st = style(3, join, 4, cap);
        for (double off : {0.0, 2.5, -1.0}) {
          if (off != 0.0 && (turn % 3)) {
            continue;
          }
          for (uint32_t scope : {0u, 1u}) {
            if (scope == 1u && p.first != "two subpaths" && (turn % 4)) {
              continue;
            }
            run(dashed("dashes " + p.first + " cap " + std::to_string(cap), p.second, st, d, off, scope));
          }
        }
      }
    }
  }
}

TEST(StrokeMatrix, DashesAtVerticesAndDots) {
  /* dashes whose cuts fall exactly on vertices (axis-aligned lengths a pattern divides), dots at corners, at the
   * ends and at the start of a closed subpath */
  Ops sq = {{'M', {0, 0}}, {'L', {10, 0}}, {'L', {10, 10}}, {'L', {0, 10}}, {'Z', {}}};
  Ops open = {{'M', {0, 0}}, {'L', {10, 0}}, {'L', {10, 10}}};
  for (int cap : {BUTT, RCAP, SQUARE}) {
    for (int join : {MITER, ROUND, BEVEL}) {
      GSHGEOM_Stroke st = style(2, join, 4, cap);
      run(dashed("on vertices", sq, st, {5, 5}, 0));
      run(dashed("on vertices, offset 5", sq, st, {10, 10}, 5));
      run(dashed("dots on a square", sq, st, {0, 5}, 0));
      run(dashed("dots on a path", open, st, {0, 10}, 0));
      run(dashed("a dash across a corner", open, st, {6, 4}, 7));
      run(dashed("a dash across the start", sq, st, {12, 6}, 3));
      run(dashed("one dash for the whole subpath", sq, st, {1000, 1}, 0));
    }
  }
  /* degenerate subpaths with dashes: the pattern is on at the point or it is not; degenerate omit draws nothing */
  Ops deg = {{'M', {2, 2}}, {'L', {2, 2}}};
  for (int cap : {RCAP, SQUARE}) {
    GSHGEOM_Stroke st = style(3, MITER, 4, cap);
    run(dashed("degenerate, on", Ops{{'M', {0, 0}}, {'L', {8, 0}}, {'M', {20, 5}}, {'L', {20, 5}}}, st, {5, 3}, 0));
    run(dashed("degenerate path scope", Ops{{'M', {0, 0}}, {'L', {8, 0}}, {'M', {20, 5}}, {'L', {20, 5}}}, st, {5, 3}, 0, 1));
    (void)deg;
  }
  /* a zero-length dash is drawn whatever `degenerate` says (SG24), and a degenerate subpath is omitted by it */
  for (int cap : {RCAP, SQUARE}) {
    GSHGEOM_Stroke st = style(3, MITER, 4, cap);
    st.degenerate = GSHGEOM_STROKE_DEGENERATE_OMIT;
    run(dashed("dots under omit", open, st, {0, 5}, 0));
    run(dashed("dots under omit on a corner", sq, st, {0, 5}, 0));
  }
}

TEST(StrokeMatrix, DashesAcrossTheStartOfAClosedSubpath) {
  /* SG23: a dash that crosses the start is one dash joined at the start vertex; one that ends or begins exactly there
   * is two; one that covers the whole subpath is its closed stroke.  The square's perimeter is 56 */
  Ops sq = {{'M', {0, 0}}, {'L', {16, 0}}, {'L', {16, 12}}, {'L', {0, 12}}, {'Z', {}}};
  Ops rect = {{'M', {0, 0}}, {'L', {20, 0}}, {'L', {20, 10}}, {'L', {0, 10}}, {'Z', {}}};
  Ops tri = {{'M', {0, 0}}, {'L', {12, 3}}, {'L', {5, 14}}, {'Z', {}}};
  Ops round = circle(0, 0, 9);
  for (int join : {MITER, ROUND, BEVEL}) {
    for (int cap : {BUTT, RCAP, SQUARE}) {
      GSHGEOM_Stroke st = style(3, join, 4, cap);
      run(dashed("wraps", sq, st, {30, 10}, 10));
      run(dashed("starts exactly at the start", sq, st, {30, 10}, 0));
      run(dashed("ends exactly at the end", rect, st, {35, 5}, 15));
      run(dashed("covers it all", sq, st, {1000, 5}, 10));
      run(dashed("wraps, a corner inside the wrapped part", tri, st, {20, 7}, 4));
      run(dashed("wraps on a circle", round, st, {30, 9}, 12));
      run(dashed("wraps, one pattern for the path", sq, st, {30, 10}, 10, 1));
    }
  }
}

TEST(StrokeMatrix, DashesOnCurves) {
  std::vector<std::pair<std::string, Ops>> curves = {
      {"quadratic", {{'M', {0, 0}}, {'Q', {10, 18, 20, 0}}}},
      {"cubic inflection", {{'M', {0, 0}}, {'C', {8, 16, 12, -16, 20, 0}}}},
      {"cubic loop", {{'M', {0, 0}}, {'C', {24, 16, -4, 16, 20, 0}}}},
      {"cubic cusp", {{'M', {0, 0}}, {'C', {16, 8, 8, 8, 8, 0}}}},
      {"rational cubic", {{'M', {0, 0}}, {'K', {6, 12, 14, 12, 20, 0, 1, 2, 0.5, 1.5}}}},
      {"circle", circle(0, 0, 10)},
      {"extreme weights", {{'M', {0, 0}}, {'R', {10, 16, 20, 0, 1, 40, 1}}}},
  };
  const int pairs[3][2] = {{MITER, BUTT}, {ROUND, RCAP}, {BEVEL, SQUARE}};
  for (const auto & cv : curves) {
    for (const auto & jc : pairs) {
      GSHGEOM_Stroke st = style(2, jc[0], 4, jc[1]);
      run(dashed(cv.first, cv.second, st, {7, 4}, 0));
      run(dashed(cv.first + " dots", cv.second, st, {0, 5}, 1.5));
    }
  }
}

TEST(StrokeMatrix, DashesUnderBothWidthSpaces) {
  Ops p = {{'M', {0, 0}}, {'L', {12, 0}}, {'Q', {18, 0, 18, 8}}, {'L', {4, 12}}};
  std::vector<std::pair<std::string, GSHGEOM_Map>> maps = {
      {"identity", map_id()},
      {"scale 2 1", map_of(2, 0, 0, 1, 0, 0)},
      {"rotation", map_of(0.6, -0.8, 0.8, 0.6, 5, 5)},
      {"skew", map_of(1, 0.75, 0, 1, 0, 0)},
      {"mirror", map_of(-1, 0, 0, 1, 40, 0)},
  };
  for (const auto & mp : maps) {
    for (uint32_t space : {0u, 1u}) {
      GSHGEOM_Stroke st = style(space ? 2 : 1.5, ROUND, 4, SQUARE);
      st.width_space = space;
      run(dashed("dashes " + mp.first, p, st, {5, 3, 1}, 1.0, 0, mp.second));
      run(dashed("dots " + mp.first, p, st, {0, 6}, 0.0, 0, mp.second));
    }
  }
}

TEST(StrokeMatrix, DashedAlignment) {
  Ops sq = {{'M', {0, 0}}, {'L', {16, 0}}, {'L', {16, 16}}, {'L', {0, 16}}, {'Z', {}}};
  Ops mixed = {{'M', {0, 0}}, {'L', {12, 0}}, {'L', {12, 12}}, {'Z', {}}, {'M', {20, 0}}, {'L', {34, 5}}};
  for (uint32_t align : {1u, 2u}) {
    for (const Ops * p : {&sq, &mixed}) {
      GSHGEOM_Stroke st = style(2, MITER, 4, BUTT);
      st.alignment = align;
      run(dashed("dashed alignment", *p, st, {6, 4}, 1));
    }
  }
}

TEST(StrokeMatrix, DashEdgesOnACircleAreAtTheExactArcLength) {
  /* a circle of radius 20 stroked at width 2 with butt caps and the pattern [10, 5]: the dash k covers the arc
   * [15 k, 15 k + 10], so the centerline point at the arc length s is inside exactly when s mod 15 < 10; probed 0.2
   * beyond every edge (past the tolerance) at the centerline's radius */
  Ops c = circle(0, 0, 20);
  GSHGEOM_Path * path = sm::build(c);
  GSHGEOM_Stroke st = style(2, MITER, 4, BUTT);
  double d[2] = {10, 5};
  st.dashes = d;
  st.dash_count = 2;
  GSHGEOM_Map m = map_id();
  GSHGEOM_Stroke_Outline * o = nullptr;
  ASSERT_EQ(gshgeom_stroke_expand(path, &st, &m, nullptr, nullptr, &o), GSHGEOM_OK);
  so::Winding w;
  w.build(so::polygons_of(gshgeom_stroke_outline_body(o)));
  double total = 2 * 3.14159265358979323846 * 20;
  int checked = 0;
  for (double s0 = 0; s0 < total - 0.5; s0 += 5.0) {
    for (double delta : {-0.2, 0.2}) {
      double edge = s0;
      double s = edge + delta;
      if (s < 0 || s > total - 0.3) {
        continue;
      }
      double mod = std::fmod(s, 15.0);
      /* only points 0.2 from an edge: edges are at multiples of 5 */
      bool inside = mod < 10.0;
      double th = s / 20.0;
      so::V2 pt = {20 * std::cos(th), 20 * std::sin(th)};
      EXPECT_EQ(w.at(pt) != 0, inside) << "arc length " << s;
      checked++;
    }
  }
  EXPECT_GT(checked, 20);
  gshgeom_stroke_outline_destroy(o);
  gshgeom_path_destroy(path);
}

TEST(StrokeMatrix, DashesAreRefusedPastTheErrorGate) {
  /* SG26: a cut error past the gate is ERR_LIMIT, never a coarser cut: 2^-46 of a polygon 2^29 long is 2^-17 */
  GSHGEOM_Path * path = sm::build(Ops{{'M', {0, 0}}, {'L', {0x1p29, 0}}});
  GSHGEOM_Stroke st = style(2, MITER, 4, BUTT);
  double d[2] = {0x1p22, 0x1p22}; /* a few dozen dashes: only the cut position's error bound refuses this */
  st.dashes = d;
  st.dash_count = 2;
  GSHGEOM_Map m = map_id();
  GSHGEOM_Stroke_Outline * o = nullptr;
  EXPECT_EQ(gshgeom_stroke_expand(path, &st, &m, nullptr, nullptr, &o), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(o, nullptr);
  gshgeom_path_destroy(path);
}

TEST(StrokeMatrix, TheDeviationProfileIsRecorded) {
  /* runs last in this file (gtest runs tests in order of definition): the worst outward deviation of any outline
   * vertex from its region over every case above, as a property and against the proof's margin */
  RecordProperty("worst_outward_deviation_units", std::to_string(g_worst));
  RecordProperty("worst_deviation_over_tolerance", std::to_string(g_worst_margin_tau));
  std::printf("stroke matrix: worst outward deviation %.4f units (%.4f of the tolerance)\n", g_worst,
      g_worst_margin_tau);
  EXPECT_GT(g_worst, 0.0) << "nothing was measured";
  EXPECT_LE(g_worst_margin_tau, 1.0);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
