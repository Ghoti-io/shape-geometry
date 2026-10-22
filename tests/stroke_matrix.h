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
 * The expansion checked against the oracle (tests/stroke_oracle.h) on one case:
 * a centerline, parameters, a map and a tolerance.  Every point the oracle
 * decides (farther than `tau + m` from every boundary sample, `m` the sample
 * spacing and a margin for the oracle's own arithmetic) must be classified by
 * the outline's nonzero fill as the region classifies it; no subpath of the
 * outline may wind negatively at any tested point; every outline vertex must lie
 * inside the conservative box of the scene's stroke bounds (restated here) plus
 * `tau`; and the deviation profile (how far an outline vertex lies outside the
 * region) is measured.  The points are a grid over that box, two probes per
 * boundary sample at `tau + m + 1/256` pixel on each side along both axes, and
 * random points.
 */

#ifndef GHOTI_IO_GSHGEOM_TESTS_STROKE_MATRIX_H
#define GHOTI_IO_GSHGEOM_TESTS_STROKE_MATRIX_H

#include "stroke_oracle.h"

#include <cstdio>
#include <random>
#include <string>
#include <vector>

namespace stroke_matrix {

namespace so = stroke_oracle;

typedef std::vector<std::pair<char, std::vector<double>>> Ops;

inline GSHGEOM_Point pt(double x, double y) {
  GSHGEOM_Point p = {x, y, 0.0};
  return p;
}

/** A path from the compact description of test_stroke_oracle.cpp. */
inline GSHGEOM_Path * build(const Ops & ops) {
  GSHGEOM_Path * p = nullptr;
  if (gshgeom_path_create(nullptr, nullptr, 2, &p) != GSHGEOM_OK) {
    return nullptr;
  }
  for (const auto & op : ops) {
    const std::vector<double> & a = op.second;
    GSHGEOM_Result r = GSHGEOM_OK;
    switch (op.first) {
    case 'M':
      r = gshgeom_path_move_to(p, pt(a[0], a[1]));
      break;
    case 'L':
      r = gshgeom_path_line_to(p, pt(a[0], a[1]));
      break;
    case 'Q':
      r = gshgeom_path_quadratic_to(p, pt(a[0], a[1]), pt(a[2], a[3]));
      break;
    case 'C':
      r = gshgeom_path_cubic_to(p, pt(a[0], a[1]), pt(a[2], a[3]), pt(a[4], a[5]));
      break;
    case 'R': {
      double w[3] = {a[4], a[5], a[6]};
      r = gshgeom_path_rational_quadratic_to(p, pt(a[0], a[1]), pt(a[2], a[3]), w);
      break;
    }
    case 'K': {
      double w[4] = {a[6], a[7], a[8], a[9]};
      r = gshgeom_path_rational_cubic_to(p, pt(a[0], a[1]), pt(a[2], a[3]), pt(a[4], a[5]), w);
      break;
    }
    case 'Z':
      r = gshgeom_path_close(p);
      break;
    }
    if (r != GSHGEOM_OK) {
      gshgeom_path_destroy(p);
      return nullptr;
    }
  }
  return p;
}

inline so::Style oracle_style(const GSHGEOM_Stroke & s) {
  so::Style o;
  o.width = s.width;
  o.space = (int)s.width_space;
  o.join = (int)s.join;
  o.limit = s.miter_limit;
  o.cap = (int)s.cap;
  o.degenerate = (int)s.degenerate;
  return o;
}

/** The scene's conservative stroke box, restated: the outset half width and the box in device pixels. */
struct Box {
  so::Real x0 = 0, y0 = 0, x1 = 0, y1 = 0;
  bool empty = true;
};

inline Box stroke_box(const GSHGEOM_Path * path, const GSHGEOM_Stroke & s, const GSHGEOM_Map & m) {
  Box b;
  GSHGEOM_Bounds mb;
  if (gshgeom_path_bounds(path, &mb) != GSHGEOM_OK) {
    return b;
  }
  so::Centerline cl = so::read_centerline(path);
  bool any = !mb.empty;
  so::Real x0 = mb.min.x, y0 = mb.min.y, x1 = mb.max.x, y1 = mb.max.y;
  for (const so::Sub & sb : cl.subs) {
    if (!sb.draws) {
      continue; /* a lone move draws nothing; a degenerate subpath is its point */
    }
    if (!any) {
      x0 = x1 = sb.start.x;
      y0 = y1 = sb.start.y;
      any = true;
    }
    x0 = std::min(x0, sb.start.x);
    y0 = std::min(y0, sb.start.y);
    x1 = std::max(x1, sb.start.x);
    y1 = std::max(y1, sb.start.y);
  }
  if (!any) {
    return b;
  }
  bool miterish = s.join == GSHGEOM_STROKE_JOIN_MITER || s.join == GSHGEOM_STROKE_JOIN_MITER_CLIP ||
      s.join == GSHGEOM_STROKE_JOIN_ARCS;
  so::Real f = 1;
  if (miterish) {
    f = std::max(f, (so::Real)s.miter_limit);
  }
  if (s.cap == GSHGEOM_STROKE_CAP_SQUARE) {
    f = std::max(f, (so::Real)1.5);
  }
  so::Real a = s.alignment == GSHGEOM_STROKE_ALIGN_CENTER ? 1 : 2;
  so::Real hb = s.width / 2 * f * a;
  so::Affine af = so::affine_of(m);
  so::Real out = 0;
  if (s.width_space == GSHGEOM_STROKE_WIDTH_LOCAL) {
    x0 -= hb;
    y0 -= hb;
    x1 += hb;
    y1 += hb;
  } else {
    out = hb;
  }
  so::V2 cs[4] = {{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}};
  b.x0 = b.y0 = 1e300L;
  b.x1 = b.y1 = -1e300L;
  for (so::V2 c : cs) {
    so::V2 d = af.apply(c);
    b.x0 = std::min(b.x0, d.x - out);
    b.y0 = std::min(b.y0, d.y - out);
    b.x1 = std::max(b.x1, d.x + out);
    b.y1 = std::max(b.y1, d.y + out);
  }
  b.empty = false;
  return b;
}

/** What one part's comparison found. */
struct PartResult {
  long tested = 0;
  long decided = 0;
  long wrong = 0;
  long negative = 0;
  long outside_box = 0;
  long ambiguous = 0;
  so::Real worst_outward = 0; ///< device pixels: the deviation profile's maximum
  std::string first;
};

/**
 * Compare one part of the outline (polygons in 26.6) with one region of the oracle, over the box outset by tau.
 * `budget` caps the number of grid points; probes and random points are added.
 */
inline PartResult compare_part(const GSHGEOM_Flat_Path * fp, const so::Region & region, const Box & box, int32_t tau,
    long budget, uint64_t seed) {
  PartResult pr;
  so::Real tau_px = (so::Real)tau / 64;
  so::Real sigma = 1.0L / 256;
  so::Real m = sigma + 1e-6L;
  so::Samples samples;
  samples.cell = tau_px + m + 2 * sigma;
  so::sample_region(samples, region, sigma);
  so::Winding wd;
  wd.build(so::polygons_of(fp));
  so::Real bx0 = box.x0 - tau_px, by0 = box.y0 - tau_px, bx1 = box.x1 + tau_px, by1 = box.y1 + tau_px;
  /* an arcs join reaches its meeting point and the arcs on the way, past what the path and the width bound: its own
   * boundary samples extend the box */
  for (const so::Shape & sh : region.shapes) {
    if (sh.kind != so::Shape::ARCS) {
      continue;
    }
    so::Region one;
    one.local = region.local;
    one.to_device = region.to_device;
    one.h = region.h;
    one.shapes.push_back(sh);
    so::Samples ts;
    ts.cell = 1;
    so::sample_region(ts, one, 1.0L / 16);
    for (const auto & cell : ts.grid) {
      for (const so::V2 & q : cell.second) {
        bx0 = std::min(bx0, q.x - tau_px);
        by0 = std::min(by0, q.y - tau_px);
        bx1 = std::max(bx1, q.x + tau_px);
        by1 = std::max(by1, q.y + tau_px);
      }
    }
  }
  /* no subpath winds negatively anywhere: each ring's fan from its first vertex has no negative triangle (then its
   * winding is a sum of triangles' 0s and 1s), decided exactly on the 26.6 integers */
  for (const auto & ring : wd.polys.rings) {
    for (size_t i = 1; i + 1 < ring.size(); i++) {
      GSHGEOM_Point a = {(double)ring[0].first, (double)ring[0].second, 0.0};
      GSHGEOM_Point b = {(double)ring[i].first, (double)ring[i].second, 0.0};
      GSHGEOM_Point c = {(double)ring[i + 1].first, (double)ring[i + 1].second, 0.0};
      int s = 0;
      gshgeom_orient2d(a, b, c, &s);
      if (s < 0) {
        pr.negative++;
        break;
      }
    }
  }
  /* every outline vertex inside the box plus tau; the deviation of those outside the region */
  for (const auto & ring : wd.polys.rings) {
    for (const auto & q : ring) {
      so::V2 v = {(so::Real)q.first / 64, (so::Real)q.second / 64};
      if (v.x < bx0 - 1e-9L || v.x > bx1 + 1e-9L || v.y < by0 - 1e-9L || v.y > by1 + 1e-9L) {
        pr.outside_box++;
      }
      so::Verdict vd = so::region_contains(region, v);
      if (!vd.inside) {
        so::Real d = samples.nearest(v, samples.cell);
        pr.worst_outward = std::max(pr.worst_outward, d);
      }
    }
  }
  auto test = [&](so::V2 x) {
    pr.tested++;
    int w = wd.at(x);
    if (w < 0) {
      pr.negative++;
    }
    if (samples.nearest(x, samples.cell) <= tau_px + m) {
      return;
    }
    so::Verdict v = so::region_contains(region, x);
    if (v.ambiguous) {
      pr.ambiguous++;
      return;
    }
    pr.decided++;
    if (v.inside != (w != 0)) {
      if (pr.wrong++ == 0) {
        char buf[256];
        std::snprintf(buf, sizeof buf, "(%.6Lf, %.6Lf): the region says %s, the outline winds %d", x.x, x.y,
            v.inside ? "inside" : "outside", w);
        pr.first = buf;
        if (v.inside) {
          /* which shape of the region holds the point: the first thing to look at */
          so::V2 xl = region.local ? region.to_device.inverse(x) : x;
          for (size_t si = 0; si < region.shapes.size(); si++) {
            if (so::shape_contains(region.shapes[si], xl)) {
              char b2[160];
              std::snprintf(b2, sizeof b2, " [shape %zu kind %d, X (%.4Lf, %.4Lf)]", si, (int)region.shapes[si].kind,
                  region.shapes[si].X.x, region.shapes[si].X.y);
              pr.first += b2;
            }
          }
        }
      }
    }
  };
  so::Real step = 1.0L / 8;
  so::Real area = (bx1 - bx0) * (by1 - by0);
  if (area / (step * step) > (so::Real)budget) {
    step = std::sqrt(area / (so::Real)budget);
  }
  for (so::Real y = by0; y <= by1; y += step) {
    for (so::Real x = bx0; x <= bx1; x += step) {
      test({x, y});
    }
  }
  /* probes beside the boundary samples (a subsample, so a long boundary does not dominate) */
  size_t stride = samples.count / 600 + 1;
  size_t idx = 0;
  so::Real off = tau_px + m + 1.0L / 256;
  for (const auto & cell : samples.grid) {
    for (const so::V2 & s : cell.second) {
      if (idx++ % stride) {
        continue;
      }
      test({s.x + off, s.y});
      test({s.x - off, s.y});
      test({s.x, s.y + off});
      test({s.x, s.y - off});
    }
  }
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> ux((double)bx0, (double)bx1), uy((double)by0, (double)by1);
  for (int k = 0; k < 512; k++) {
    test({ux(rng), uy(rng)});
  }
  return pr;
}

} // namespace stroke_matrix

#endif /* GHOTI_IO_GSHGEOM_TESTS_STROKE_MATRIX_H */
