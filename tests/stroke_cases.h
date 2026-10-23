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
 * The stroke case matrix (tests/unit/test_stroke_matrix.cpp): lines and every join and cap, the miter boundary, every
 * curve kind, both width spaces under every map, alignment, overlaps, every tolerance and the dashed families, as
 * lists of cases that the polygon expansion's matrix and the curve outline's both run.
 */

#ifndef GHOTI_IO_GSHGEOM_TESTS_STROKE_CASES_H
#define GHOTI_IO_GSHGEOM_TESTS_STROKE_CASES_H

#include "stroke_matrix.h"

#include <string>
#include <vector>

namespace stroke_cases {

namespace sm = stroke_matrix;

typedef sm::Ops Ops;


inline GSHGEOM_Map map_id() {
  GSHGEOM_Map m;
  gshgeom_map_identity(2, &m);
  return m;
}

inline GSHGEOM_Map map_of(double a, double b, double c, double d, double tx, double ty) {
  double mm[4] = {a, b, c, d}, t[2] = {tx, ty};
  GSHGEOM_Map m;
  EXPECT_EQ(gshgeom_map_from_matrix(2, mm, t, &m), GSHGEOM_OK);
  return m;
}

inline GSHGEOM_Stroke style(double width, int join, double limit, int cap) {
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


const int MITER = GSHGEOM_STROKE_JOIN_MITER, ROUND = GSHGEOM_STROKE_JOIN_ROUND, BEVEL = GSHGEOM_STROKE_JOIN_BEVEL,
          MCLIP = GSHGEOM_STROKE_JOIN_MITER_CLIP, ARCS = GSHGEOM_STROKE_JOIN_ARCS;
const int BUTT = GSHGEOM_STROKE_CAP_BUTT, RCAP = GSHGEOM_STROKE_CAP_ROUND, SQUARE = GSHGEOM_STROKE_CAP_SQUARE;

inline Ops circle(double cx, double cy, double r) {
  double w = 0.70710678118654752440;
  return Ops{{'M', {cx + r, cy}}, {'R', {cx + r, cy + r, cx, cy + r, 1, w, 1}}, {'R', {cx - r, cy + r, cx - r, cy, 1, w, 1}},
      {'R', {cx - r, cy - r, cx, cy - r, 1, w, 1}}, {'R', {cx + r, cy - r, cx + r, cy, 1, w, 1}}, {'Z', {}}};
}


inline Case dashed(const std::string & name, const Ops & ops, GSHGEOM_Stroke st, std::vector<double> d, double offset,
    uint32_t scope = 0, GSHGEOM_Map m = map_id()) {
  st.dash_offset = offset;
  st.dash_scope = scope;
  Case c{name, ops, st, m, GSHGEOM_FLATTEN_TOLERANCE};
  c.dashes = d;
  return c;
}


inline void add(std::vector<Case> & out, const Case & c) {
  out.push_back(c);
}

inline std::vector<Case> linesEveryJoinAndCap() {
  std::vector<Case> out;
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
        add(out, c);
      }
      turn++;
    }
  }
  return out;
}

inline std::vector<Case> theMiterBoundaryIsDecidedExactly() {
  std::vector<Case> out;
  Ops tri = {{'M', {0, 0}}, {'L', {4, 3}}, {'L', {8, 0}}};
  for (double L : {1.25, std::nextafter(1.25, 0.0), 1.1, 3.0}) {
    for (int join : {MITER, MCLIP}) {
      add(out, {"miter boundary", tri, style(8, join, L, BUTT), map_id(), GSHGEOM_FLATTEN_TOLERANCE});
    }
  }
  return out;
}

inline std::vector<Case> everyCurveKind() {
  std::vector<Case> out;
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
      add(out, {cv.first, cv.second, style(4, jc[0], 4, jc[1]), map_id(), GSHGEOM_FLATTEN_TOLERANCE});
    }
  }
  return out;
}

inline std::vector<Case> bothWidthSpacesUnderEveryMap() {
  std::vector<Case> out;
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
        add(out, {mp.first, q, s, mp.second, GSHGEOM_FLATTEN_TOLERANCE});
      }
    }
  }
  /* the hairline: device width 1 */
  GSHGEOM_Stroke hair = style(1, MITER, 4, BUTT);
  hair.width_space = GSHGEOM_STROKE_WIDTH_DEVICE;
  add(out, {"hairline", p, hair, map_of(3, 0, 0, 3, 0, 0), GSHGEOM_FLATTEN_TOLERANCE});
  return out;
}

inline std::vector<Case> alignment() {
  std::vector<Case> out;
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
        add(out, {p.first, p.second, s, map_id(), GSHGEOM_FLATTEN_TOLERANCE});
      }
    }
  }
  return out;
}

inline std::vector<Case> overlapsAreOneUnion() {
  std::vector<Case> out;
  std::vector<std::pair<std::string, Ops>> paths = {
      {"figure eight", {{'M', {0, 0}}, {'C', {10, 10, 10, -10, 0, 0}}, {'C', {-10, 10, -10, -10, 0, 0}}, {'Z', {}}}},
      {"retrace", {{'M', {0, 0}}, {'L', {10, 0}}, {'L', {0, 0}}, {'L', {10, 0}}}},
      {"overlapping subpaths", {{'M', {0, 0}}, {'L', {10, 10}}, {'M', {0, 10}}, {'Q', {5, -6, 10, 0}}}},
      {"loop", {{'M', {0, 0}}, {'C', {14, 10, -4, 10, 10, 0}}}},
      {"coincident", {{'M', {0, 0}}, {'L', {8, 0}}, {'M', {0, 0}}, {'L', {8, 0}}, {'M', {8, 0}}, {'L', {0, 0}}}},
  };
  for (const auto & p : paths) {
    for (int join : {MITER, ROUND, BEVEL}) {
      add(out, {p.first, p.second, style(3, join, 4, RCAP), map_id(), GSHGEOM_FLATTEN_TOLERANCE});
    }
  }
  return out;
}

inline std::vector<Case> everyTolerance() {
  std::vector<Case> out;
  Ops p = {{'M', {0, 0}}, {'C', {8, 12, 14, -6, 20, 4}}, {'L', {12, 10}}};
  for (int32_t tau : {2, 3, 4, 16, 64}) {
    for (int join : {MITER, ROUND, ARCS}) {
      add(out, {"tolerance " + std::to_string(tau), p, style(5, join, 4, RCAP), map_id(), tau});
    }
  }
  return out;
}

inline std::vector<Case> dashesOnLines() {
  std::vector<Case> out;
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
            add(out, dashed("dashes " + p.first + " cap " + std::to_string(cap), p.second, st, d, off, scope));
          }
        }
      }
    }
  }
  return out;
}

inline std::vector<Case> dashesAtVerticesAndDots() {
  std::vector<Case> out;
  /* dashes whose cuts fall exactly on vertices (axis-aligned lengths a pattern divides), dots at corners, at the
   * ends and at the start of a closed subpath */
  Ops sq = {{'M', {0, 0}}, {'L', {10, 0}}, {'L', {10, 10}}, {'L', {0, 10}}, {'Z', {}}};
  Ops open = {{'M', {0, 0}}, {'L', {10, 0}}, {'L', {10, 10}}};
  for (int cap : {BUTT, RCAP, SQUARE}) {
    for (int join : {MITER, ROUND, BEVEL}) {
      GSHGEOM_Stroke st = style(2, join, 4, cap);
      add(out, dashed("on vertices", sq, st, {5, 5}, 0));
      add(out, dashed("on vertices, offset 5", sq, st, {10, 10}, 5));
      add(out, dashed("dots on a square", sq, st, {0, 5}, 0));
      add(out, dashed("dots on a path", open, st, {0, 10}, 0));
      add(out, dashed("a dash across a corner", open, st, {6, 4}, 7));
      add(out, dashed("a dash across the start", sq, st, {12, 6}, 3));
      add(out, dashed("one dash for the whole subpath", sq, st, {1000, 1}, 0));
    }
  }
  /* degenerate subpaths with dashes: the pattern is on at the point or it is not; degenerate omit draws nothing */
  Ops deg = {{'M', {2, 2}}, {'L', {2, 2}}};
  for (int cap : {RCAP, SQUARE}) {
    GSHGEOM_Stroke st = style(3, MITER, 4, cap);
    add(out, dashed("degenerate, on", Ops{{'M', {0, 0}}, {'L', {8, 0}}, {'M', {20, 5}}, {'L', {20, 5}}}, st, {5, 3}, 0));
    add(out, dashed("degenerate path scope", Ops{{'M', {0, 0}}, {'L', {8, 0}}, {'M', {20, 5}}, {'L', {20, 5}}}, st, {5, 3}, 0, 1));
    (void)deg;
  }
  /* a zero-length dash is drawn whatever `degenerate` says (SG24), and a degenerate subpath is omitted by it */
  for (int cap : {RCAP, SQUARE}) {
    GSHGEOM_Stroke st = style(3, MITER, 4, cap);
    st.degenerate = GSHGEOM_STROKE_DEGENERATE_OMIT;
    add(out, dashed("dots under omit", open, st, {0, 5}, 0));
    add(out, dashed("dots under omit on a corner", sq, st, {0, 5}, 0));
  }
  return out;
}

inline std::vector<Case> dashesAcrossTheStartOfAClosedSubpath() {
  std::vector<Case> out;
  /* SG23: a dash that crosses the start is one dash joined at the start vertex; one that ends or begins exactly there
   * is two; one that covers the whole subpath is its closed stroke.  The square's perimeter is 56 */
  Ops sq = {{'M', {0, 0}}, {'L', {16, 0}}, {'L', {16, 12}}, {'L', {0, 12}}, {'Z', {}}};
  Ops rect = {{'M', {0, 0}}, {'L', {20, 0}}, {'L', {20, 10}}, {'L', {0, 10}}, {'Z', {}}};
  Ops tri = {{'M', {0, 0}}, {'L', {12, 3}}, {'L', {5, 14}}, {'Z', {}}};
  Ops round = circle(0, 0, 9);
  for (int join : {MITER, ROUND, BEVEL}) {
    for (int cap : {BUTT, RCAP, SQUARE}) {
      GSHGEOM_Stroke st = style(3, join, 4, cap);
      add(out, dashed("wraps", sq, st, {30, 10}, 10));
      add(out, dashed("starts exactly at the start", sq, st, {30, 10}, 0));
      add(out, dashed("ends exactly at the end", rect, st, {35, 5}, 15));
      add(out, dashed("covers it all", sq, st, {1000, 5}, 10));
      add(out, dashed("wraps, a corner inside the wrapped part", tri, st, {20, 7}, 4));
      add(out, dashed("wraps on a circle", round, st, {30, 9}, 12));
      add(out, dashed("wraps, one pattern for the path", sq, st, {30, 10}, 10, 1));
    }
  }
  return out;
}

inline std::vector<Case> dashesOnCurves() {
  std::vector<Case> out;
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
      add(out, dashed(cv.first, cv.second, st, {7, 4}, 0));
      add(out, dashed(cv.first + " dots", cv.second, st, {0, 5}, 1.5));
    }
  }
  return out;
}

inline std::vector<Case> dashesUnderBothWidthSpaces() {
  std::vector<Case> out;
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
      add(out, dashed("dashes " + mp.first, p, st, {5, 3, 1}, 1.0, 0, mp.second));
      add(out, dashed("dots " + mp.first, p, st, {0, 6}, 0.0, 0, mp.second));
    }
  }
  return out;
}

inline std::vector<Case> dashedAlignment() {
  std::vector<Case> out;
  Ops sq = {{'M', {0, 0}}, {'L', {16, 0}}, {'L', {16, 16}}, {'L', {0, 16}}, {'Z', {}}};
  Ops mixed = {{'M', {0, 0}}, {'L', {12, 0}}, {'L', {12, 12}}, {'Z', {}}, {'M', {20, 0}}, {'L', {34, 5}}};
  for (uint32_t align : {1u, 2u}) {
    for (const Ops * p : {&sq, &mixed}) {
      GSHGEOM_Stroke st = style(2, MITER, 4, BUTT);
      st.alignment = align;
      add(out, dashed("dashed alignment", *p, st, {6, 4}, 1));
    }
  }
  return out;
}

/** Every case of the families above, in the matrix's order. */
inline std::vector<Case> everything() {
  std::vector<Case> all;
  for (auto f : {&linesEveryJoinAndCap, &theMiterBoundaryIsDecidedExactly, &everyCurveKind, &bothWidthSpacesUnderEveryMap,
           &alignment, &overlapsAreOneUnion, &everyTolerance, &dashesOnLines, &dashesAtVerticesAndDots,
           &dashesAcrossTheStartOfAClosedSubpath, &dashesOnCurves, &dashesUnderBothWidthSpaces, &dashedAlignment}) {
    std::vector<Case> v = f();
    all.insert(all.end(), v.begin(), v.end());
  }
  return all;
}

} // namespace stroke_cases

#endif /* GHOTI_IO_GSHGEOM_TESTS_STROKE_CASES_H */
