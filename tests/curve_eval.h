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
 * The curve-winding evaluator of the curve outline's tests (stroke_curves.h): the
 * nonzero fill of a path of lines, quadratics, cubics, rational quadratics and
 * rational cubics, read through the public readback and shared with nothing in
 * src/.
 *
 * Each subpath is a closed ring.  A curve is flattened in `long double` by
 * adaptive bisection (de Casteljau in homogeneous coordinates) until every control
 * point lies within `flat` of its chord, so the polyline is within `flat` of the
 * curve (a curve with positive weights lies in the hull of its control points);
 * the winding number at a point is the signed crossing count of the horizontal ray
 * to +x, by the polyline of each ring.  The evaluator is first checked where the
 * answer is exact (tests/unit/test_stroke_curves_eval.cpp): circles and ellipses
 * built from rational quadratics, and cubic approximations of them.
 */

#ifndef GHOTI_IO_GSHGEOM_TESTS_CURVE_EVAL_H
#define GHOTI_IO_GSHGEOM_TESTS_CURVE_EVAL_H

#include <ghoti.io/shape-geometry/shape-geometry.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace curve_eval {

typedef long double Real;

/** A defect planted in the evaluator itself, which its self-check must see: 1 flips the sign of a rational term
 * (the middle weight `w` becomes `2 - w`). */
inline int plant = 0;

struct P2 {
  Real x = 0, y = 0;
};

/** One curve: degree 1 to 3, control points, weights (1 for a polynomial kind). */
struct Elem {
  int n = 1;
  P2 p[4];
  Real w[4] = {1, 1, 1, 1};
};

inline Elem elem_of(const GSHGEOM_Element & e) {
  Elem r;
  r.n = (int)gshgeom_curve_point_count(e.curve.kind) - 1;
  bool rat = gshgeom_curve_is_rational(e.curve.kind) != 0;
  for (int k = 0; k <= r.n; k++) {
    r.p[k] = {(Real)e.curve.p[k].x, (Real)e.curve.p[k].y};
    r.w[k] = rat ? (Real)e.curve.w[k] : 1;
  }
  if (plant == 1 && rat && r.n == 2) {
    r.w[1] = 2 - r.w[1];
  }
  return r;
}

/** The point at t, by de Casteljau on the homogeneous coordinates. */
inline P2 eval(const Elem & e, Real t) {
  Real x[4], y[4], w[4];
  for (int k = 0; k <= e.n; k++) {
    x[k] = e.w[k] * e.p[k].x;
    y[k] = e.w[k] * e.p[k].y;
    w[k] = e.w[k];
  }
  for (int lev = 1; lev <= e.n; lev++) {
    for (int i = 0; i + lev <= e.n; i++) {
      x[i] = (1 - t) * x[i] + t * x[i + 1];
      y[i] = (1 - t) * y[i] + t * y[i + 1];
      w[i] = (1 - t) * w[i] + t * w[i + 1];
    }
  }
  return {x[0] / w[0], y[0] / w[0]};
}

inline void split(const Elem & e, Elem & l, Elem & r) {
  Real x[4], y[4], w[4];
  for (int k = 0; k <= e.n; k++) {
    x[k] = e.w[k] * e.p[k].x;
    y[k] = e.w[k] * e.p[k].y;
    w[k] = e.w[k];
  }
  l.n = r.n = e.n;
  auto put = [&](Elem & d, int k, int src) {
    d.w[k] = w[src];
    d.p[k] = {x[src] / w[src], y[src] / w[src]};
  };
  put(l, 0, 0);
  put(r, e.n, e.n);
  for (int lev = 1; lev <= e.n; lev++) {
    for (int i = 0; i + lev <= e.n; i++) {
      x[i] = (x[i] + x[i + 1]) / 2;
      y[i] = (y[i] + y[i + 1]) / 2;
      w[i] = (w[i] + w[i + 1]) / 2;
    }
    put(l, lev, 0);
    put(r, e.n - lev, e.n - lev);
  }
}

inline Real dist_to_segment(P2 p, P2 a, P2 b) {
  Real dx = b.x - a.x, dy = b.y - a.y;
  Real l2 = dx * dx + dy * dy;
  Real t = l2 > 0 ? ((p.x - a.x) * dx + (p.y - a.y) * dy) / l2 : 0;
  t = std::max<Real>(0, std::min<Real>(1, t));
  Real ex = p.x - (a.x + t * dx), ey = p.y - (a.y + t * dy);
  return std::sqrt(ex * ex + ey * ey);
}

/** Append the polyline of e (excluding its first point) to out, within flat of the curve. */
inline void flatten(const Elem & e, Real flat, std::vector<P2> & out, int depth = 0) {
  Real worst = 0;
  for (int k = 1; k < e.n; k++) {
    worst = std::max(worst, dist_to_segment(e.p[k], e.p[0], e.p[e.n]));
  }
  if (e.n == 1 || worst <= flat || depth >= 40) {
    out.push_back(e.p[e.n]);
    return;
  }
  Elem l, r;
  split(e, l, r);
  flatten(l, flat, out, depth + 1);
  flatten(r, flat, out, depth + 1);
}

/** The rings of a path, flattened. */
struct Rings {
  std::vector<std::vector<P2>> ring;
};

inline Rings rings_of(const GSHGEOM_Path * path, Real flat) {
  Rings out;
  if (!path) {
    return out;
  }
  size_t n = gshgeom_path_element_count(path);
  std::vector<P2> cur;
  auto flush = [&]() {
    if (cur.size() >= 3) {
      out.ring.push_back(cur);
    }
    cur.clear();
  };
  for (size_t i = 0; i < n; i++) {
    GSHGEOM_Element el;
    gshgeom_path_element(path, i, &el);
    if (el.verb == GSHGEOM_VERB_MOVE) {
      flush();
      cur.push_back({(Real)el.curve.p[0].x, (Real)el.curve.p[0].y});
    } else if (el.verb == GSHGEOM_VERB_CLOSE) {
      flush();
    } else {
      if (cur.empty()) {
        cur.push_back({(Real)el.curve.p[0].x, (Real)el.curve.p[0].y});
      }
      flatten(elem_of(el), flat, cur);
    }
  }
  flush();
  return out;
}

/** The winding number of one ring at x: signed crossings of the ray to +x. */
inline int ring_winding(const std::vector<P2> & ring, P2 x) {
  int w = 0;
  size_t n = ring.size();
  for (size_t i = 0; i < n; i++) {
    P2 a = ring[i], b = ring[(i + 1) % n];
    if (a.y <= x.y) {
      if (b.y > x.y && (b.x - a.x) * (x.y - a.y) - (x.x - a.x) * (b.y - a.y) > 0) {
        w++;
      }
    } else if (b.y <= x.y && (b.x - a.x) * (x.y - a.y) - (x.x - a.x) * (b.y - a.y) < 0) {
      w--;
    }
  }
  return w;
}

/** A path's rings bucketed by row, for many queries. */
struct Winding {
  Rings rings;
  struct Edge {
    Real ax, ay, bx, by;
  };
  Real band = 1;
  Real lo = 0;
  std::vector<std::vector<Edge>> rows;

  void build(const GSHGEOM_Path * path, Real flat) {
    rings = rings_of(path, flat);
    Real mn = 1e300L, mx = -1e300L;
    for (auto & r : rings.ring) {
      for (auto & q : r) {
        mn = std::min(mn, q.y);
        mx = std::max(mx, q.y);
      }
    }
    rows.clear();
    if (mn > mx) {
      return;
    }
    band = std::max<Real>((mx - mn) / 512, 1e-3L);
    lo = std::floor(mn / band);
    int64_t nrow = (int64_t)(std::floor(mx / band) - lo) + 1;
    rows.assign((size_t)nrow, {});
    for (auto & r : rings.ring) {
      size_t n = r.size();
      for (size_t i = 0; i < n; i++) {
        Edge e = {r[i].x, r[i].y, r[(i + 1) % n].x, r[(i + 1) % n].y};
        int64_t a = (int64_t)(std::floor(std::min(e.ay, e.by) / band) - lo);
        int64_t b = (int64_t)(std::floor(std::max(e.ay, e.by) / band) - lo);
        for (int64_t k = a; k <= b; k++) {
          rows[(size_t)k].push_back(e);
        }
      }
    }
  }

  int at(P2 x) const {
    int64_t k = (int64_t)(std::floor(x.y / band) - lo);
    if (k < 0 || (size_t)k >= rows.size()) {
      return 0;
    }
    int w = 0;
    for (const Edge & e : rows[(size_t)k]) {
      Real c = (e.bx - e.ax) * (x.y - e.ay) - (x.x - e.ax) * (e.by - e.ay);
      if (e.ay <= x.y) {
        if (e.by > x.y && c > 0) {
          w++;
        }
      } else if (e.by <= x.y && c < 0) {
        w--;
      }
    }
    return w;
  }

  /** The least winding any single ring has at x. */
  int least_ring(P2 x) const {
    int least = 0;
    for (auto & r : rings.ring) {
      least = std::min(least, ring_winding(r, x));
    }
    return least;
  }
};

/** Points along every curve of a path (for the deviation profile): `per` samples per curve element. */
inline std::vector<P2> sample_path(const GSHGEOM_Path * path, int per) {
  std::vector<P2> out;
  size_t n = path ? gshgeom_path_element_count(path) : 0;
  for (size_t i = 0; i < n; i++) {
    GSHGEOM_Element el;
    gshgeom_path_element(path, i, &el);
    if (el.verb == GSHGEOM_VERB_MOVE || el.verb == GSHGEOM_VERB_CLOSE) {
      continue;
    }
    Elem e = elem_of(el);
    for (int k = 0; k <= per; k++) {
      out.push_back(eval(e, (Real)k / per));
    }
  }
  return out;
}

} // namespace curve_eval

#endif /* GHOTI_IO_GSHGEOM_TESTS_CURVE_EVAL_H */
