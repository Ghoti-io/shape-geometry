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
 * An independent oracle for the copies of a lattice (lattice.h, LT1), by the definition and not by the library's method.
 *
 * The definition: copy `(i, j)` is a copy of the use when the convex hull of the four exact points `Q_k + tau + i A + j B`
 * meets the closed rectangle.  The library decides it with sixteen half-plane constraints, an estimate and a gallop.  This
 * file decides it directly, by geometry: the hull meets the rectangle exactly when a point of the four lies in the rectangle,
 * or one of the six segments between the points meets one of the rectangle's four sides (or lies inside it), or a corner of
 * the rectangle lies in a triangle of three of the points.  Two decisions are made here, one in exact rationals
 * (::lattice_oracle::meets_rat, any doubles) and one in 64-bit integers for inputs that are multiples of 2^-6 and small
 * (::lattice_oracle::meets_grid, fast enough to visit every index of a window).
 */

#ifndef GHOTI_IO_GSHGEOM_TEST_LATTICE_ORACLE_H
#define GHOTI_IO_GSHGEOM_TEST_LATTICE_ORACLE_H

#include "clip_oracle.h"

#include <cmath>
#include <cstdint>
#include <vector>

namespace lattice_oracle {

using clip_oracle::Big;
using clip_oracle::Rat;

/** The exact rational equal to the finite double @p v. */
inline Rat from_double(double v) {
  if (v == 0.0) {
    return Rat(0);
  }
  int e = 0;
  const double f = std::frexp(v, &e);  // v = f 2^e, 0.5 <= |f| < 1
  const int64_t m = (int64_t)std::ldexp(f, 53);  // exact: |m| < 2^53
  e -= 53;
  Big p(m);
  Big two30((int64_t)1 << 30);
  Big pow2(1);
  for (int k = 0; k < (e < 0 ? -e : e) / 30; k++) pow2 = pow2 * two30;
  for (int k = 0; k < (e < 0 ? -e : e) % 30; k++) pow2 = pow2 * Big(2);
  return e >= 0 ? Rat(p * pow2, Big(1)) : Rat(p, pow2);
}

struct RPt {
  Rat x, y;
};

inline int orient(const RPt & a, const RPt & b, const RPt & c) {
  return ((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)).sign();
}

inline bool between(const Rat & v, const Rat & a, const Rat & b) {
  return (clip_oracle::cmp(a, v) <= 0 && clip_oracle::cmp(v, b) <= 0) || (clip_oracle::cmp(b, v) <= 0 && clip_oracle::cmp(v, a) <= 0);
}

inline bool on_segment(const RPt & a, const RPt & b, const RPt & p) {
  return orient(a, b, p) == 0 && between(p.x, a.x, b.x) && between(p.y, a.y, b.y);
}

/** Whether the closed segments ab and cd share a point. */
inline bool segments_meet(const RPt & a, const RPt & b, const RPt & c, const RPt & d) {
  const int o1 = orient(a, b, c);
  const int o2 = orient(a, b, d);
  const int o3 = orient(c, d, a);
  const int o4 = orient(c, d, b);
  if (o1 * o2 < 0 && o3 * o4 < 0) return true;
  return on_segment(a, b, c) || on_segment(a, b, d) || on_segment(c, d, a) || on_segment(c, d, b);
}

inline bool in_rect(const RPt & p, const RPt r[4]) {
  return clip_oracle::cmp(r[0].x, p.x) <= 0 && clip_oracle::cmp(p.x, r[2].x) <= 0 && clip_oracle::cmp(r[0].y, p.y) <= 0 &&
      clip_oracle::cmp(p.y, r[2].y) <= 0;
}

inline bool in_triangle(const RPt & a, const RPt & b, const RPt & c, const RPt & p) {
  const int s = orient(a, b, c);
  if (s == 0) return false;
  const int o1 = orient(a, b, p) * s;
  const int o2 = orient(b, c, p) * s;
  const int o3 = orient(c, a, p) * s;
  return o1 >= 0 && o2 >= 0 && o3 >= 0;
}

/** Whether the convex hull of the four points @p p meets the closed rectangle @p r (corners in order, (x0,y0) first). */
inline bool hull_meets(const RPt p[4], const RPt r[4]) {
  for (int k = 0; k < 4; k++) {
    if (in_rect(p[k], r)) return true;
  }
  for (int a = 0; a < 4; a++) {
    for (int b = a + 1; b < 4; b++) {
      for (int s = 0; s < 4; s++) {
        if (segments_meet(p[a], p[b], r[s], r[(s + 1) % 4])) return true;
      }
    }
  }
  for (int c = 0; c < 4; c++) {
    for (int a = 0; a < 4; a++) {
      for (int b = a + 1; b < 4; b++) {
        for (int d = b + 1; d < 4; d++) {
          if (in_triangle(p[a], p[b], p[d], r[c])) return true;
        }
      }
    }
  }
  return false;
}

/** The inputs of a copy's position, all doubles: the quad of copy (0,0), tau, and the step vectors. */
struct Geometry {
  double Q[4][2];
  double tau[2];
  double A[2];
  double B[2];
};

/** Q, tau, A, B from a map's linear part and translation and a lattice, by the formulas of the header, in plain doubles. */
inline Geometry geometry_of(const double m[4], const double t[2], const double cell[4], const double steps[2]) {
  Geometry g;
  const double px[4] = {cell[0], cell[0] + cell[2], cell[0] + cell[2], cell[0]};
  const double py[4] = {cell[1], cell[1], cell[1] + cell[3], cell[1] + cell[3]};
  for (int k = 0; k < 4; k++) {
    volatile double a = m[0] * px[k];
    volatile double b = m[1] * py[k];
    g.Q[k][0] = a + b;
    volatile double c = m[2] * px[k];
    volatile double d = m[3] * py[k];
    g.Q[k][1] = c + d;
  }
  g.tau[0] = t[0];
  g.tau[1] = t[1];
  g.A[0] = m[0] * steps[0];
  g.A[1] = m[2] * steps[0];
  g.B[0] = m[1] * steps[1];
  g.B[1] = m[3] * steps[1];
  return g;
}

/** Whether copy (i, j) meets the closed rectangle {x0, y0, x1, y1}, in exact rationals. */
inline bool meets_rat(const Geometry & g, int64_t i, int64_t j, const int64_t rect[4]) {
  RPt p[4];
  const Rat ri(i);
  const Rat rj(j);
  for (int k = 0; k < 4; k++) {
    p[k].x = from_double(g.Q[k][0]) + from_double(g.tau[0]) + ri * from_double(g.A[0]) + rj * from_double(g.B[0]);
    p[k].y = from_double(g.Q[k][1]) + from_double(g.tau[1]) + ri * from_double(g.A[1]) + rj * from_double(g.B[1]);
  }
  const RPt r[4] = {{Rat(rect[0]), Rat(rect[1])}, {Rat(rect[2]), Rat(rect[1])}, {Rat(rect[2]), Rat(rect[3])}, {Rat(rect[0]), Rat(rect[3])}};
  return hull_meets(p, r);
}

/* --- the same in 64-bit integers, for values that are multiples of 2^-6 ----------------------------------------------- */

struct GPt {
  int64_t x, y;
};

inline int gorient(const GPt & a, const GPt & b, const GPt & c) {
  const int64_t v = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
  return (v > 0) - (v < 0);
}
inline bool gbetween(int64_t v, int64_t a, int64_t b) {
  return (a <= v && v <= b) || (b <= v && v <= a);
}
inline bool gon(const GPt & a, const GPt & b, const GPt & p) {
  return gorient(a, b, p) == 0 && gbetween(p.x, a.x, b.x) && gbetween(p.y, a.y, b.y);
}
inline bool gsegs(const GPt & a, const GPt & b, const GPt & c, const GPt & d) {
  const int o1 = gorient(a, b, c), o2 = gorient(a, b, d), o3 = gorient(c, d, a), o4 = gorient(c, d, b);
  if (o1 * o2 < 0 && o3 * o4 < 0) return true;
  return gon(a, b, c) || gon(a, b, d) || gon(c, d, a) || gon(c, d, b);
}

/** The scaled integer of a multiple of 2^-6 (the caller guarantees it is one). */
inline int64_t to_grid(double v) {
  return (int64_t)std::llround(v * 64.0);
}

/** Geometry scaled by 64, as integers, and the decision for any index. */
struct Grid {
  int64_t Q[4][2], tau[2], A[2], B[2];
  GPt r[4];
  explicit Grid(const Geometry & g, const int64_t rect[4]) {
    for (int k = 0; k < 4; k++) {
      Q[k][0] = to_grid(g.Q[k][0]);
      Q[k][1] = to_grid(g.Q[k][1]);
    }
    for (int d = 0; d < 2; d++) {
      tau[d] = to_grid(g.tau[d]);
      A[d] = to_grid(g.A[d]);
      B[d] = to_grid(g.B[d]);
    }
    r[0] = {rect[0] * 64, rect[1] * 64};
    r[1] = {rect[2] * 64, rect[1] * 64};
    r[2] = {rect[2] * 64, rect[3] * 64};
    r[3] = {rect[0] * 64, rect[3] * 64};
  }
  bool meets(int64_t i, int64_t j) const {
    GPt p[4];
    for (int k = 0; k < 4; k++) {
      p[k] = {Q[k][0] + tau[0] + i * A[0] + j * B[0], Q[k][1] + tau[1] + i * A[1] + j * B[1]};
    }
    for (int k = 0; k < 4; k++) {
      if (r[0].x <= p[k].x && p[k].x <= r[2].x && r[0].y <= p[k].y && p[k].y <= r[2].y) return true;
    }
    for (int a = 0; a < 4; a++)
      for (int b = a + 1; b < 4; b++)
        for (int s = 0; s < 4; s++)
          if (gsegs(p[a], p[b], r[s], r[(s + 1) % 4])) return true;
    for (int c = 0; c < 4; c++) {
      for (int a = 0; a < 4; a++)
        for (int b = a + 1; b < 4; b++)
          for (int d = b + 1; d < 4; d++) {
            const int s = gorient(p[a], p[b], p[d]);
            if (s == 0) continue;
            if (gorient(p[a], p[b], r[c]) * s >= 0 && gorient(p[b], p[d], r[c]) * s >= 0 && gorient(p[d], p[a], r[c]) * s >= 0) return true;
          }
    }
    return false;
  }
};

struct Row {
  int64_t j, lo, hi;
  int64_t n = 0;  // the copies of the span that meet: hi - lo + 1 when the row is an interval
  bool operator==(const Row & o) const { return j == o.j && lo == o.lo && hi == o.hi; }
};

/** Every copy of the window [-w, w]^2 that meets the rectangle, as rows in lattice order (nonrepeating axes restricted to 0). */
inline std::vector<Row> brute_grid(const Grid & gd, bool rx, bool ry, int64_t w) {
  std::vector<Row> rows;
  for (int64_t j = ry ? -w : 0; j <= (ry ? w : 0); j++) {
    int64_t lo = 0, hi = -1, n = 0;
    bool open = false;
    for (int64_t i = rx ? -w : 0; i <= (rx ? w : 0); i++) {
      if (gd.meets(i, j)) {
        if (!open) {
          lo = i;
          open = true;
        }
        hi = i;
        n++;
      }
    }
    if (open) {
      // report the whole span; the caller checks that the span has no holes
      Row row{j, lo, hi};
      row.n = n;
      rows.push_back(row);
    }
  }
  return rows;
}

}  // namespace lattice_oracle

#endif /* GHOTI_IO_GSHGEOM_TEST_LATTICE_ORACLE_H */
