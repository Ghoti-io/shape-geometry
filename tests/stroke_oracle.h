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
 * The geometric oracle of a stroke: which points a stroke covers, decided
 * from the definition in stroke.h (SG5 to SG13) and from nothing else.
 *
 * It shares no code with the expansion (src/stroke/): it reads a path through
 * the public readback, works in `long double`, and decides membership the way
 * the rules are written.  A point is in a segment's body when it lies on one of
 * the segment's normal segments: the parameters where `(x - c(t)) . c'(t) = 0`
 * are the roots of a polynomial (for `c = N/D`, of `(x D - N) . E` with
 * `E = N'D - ND'`), isolated in the Bernstein basis by sign variations and
 * bisected; the root is a member when `|x - c(t)| <= h`.  A cusp (`E(t*) = 0`
 * inside the segment) adds its disc (SG6).  Joins (SG8, the 2018 SVG 2 CR's arcs
 * included), caps (SG7) and degenerate shapes (SG11) are tested as the shapes
 * they are.  Local width maps the point back into model space; device width
 * maps the centerline forward.
 *
 * The oracle also samples a superset of the region's boundary (rails, end
 * normals, the evolute within reach, cusp circles, every join, cap and
 * degenerate edge), at most `sigma` apart, so that a point farther than
 * `tau + m` from every sample is decided: the expansion's outline must classify
 * it as the region does.
 *
 * It is first checked where it can be checked exactly: with round joins and
 * round caps the region is the set of points within `h` of the centerline
 * (tests/unit/test_stroke_oracle.cpp, `OracleSelfCheck`).
 */

#ifndef GHOTI_IO_GSHGEOM_TESTS_STROKE_ORACLE_H
#define GHOTI_IO_GSHGEOM_TESTS_STROKE_ORACLE_H

#include <ghoti.io/shape-geometry/shape-geometry.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace stroke_oracle {

typedef long double Real;

/** pi, to the precision of a long double. */
const Real kPi = 3.14159265358979323846264338327950288L;

/** A point or a vector. */
struct V2 {
  Real x = 0, y = 0;
};

inline V2 operator+(V2 a, V2 b) { return {a.x + b.x, a.y + b.y}; }
inline V2 operator-(V2 a, V2 b) { return {a.x - b.x, a.y - b.y}; }
inline V2 operator*(Real s, V2 a) { return {s * a.x, s * a.y}; }
inline Real dot(V2 a, V2 b) { return a.x * b.x + a.y * b.y; }
inline Real cross(V2 a, V2 b) { return a.x * b.y - a.y * b.x; }
inline Real norm(V2 a) { return std::sqrt(dot(a, a)); }
inline V2 unit(V2 a) {
  Real n = norm(a);
  return n > 0 ? (1 / n) * a : V2{0, 0};
}
/** orient.h's left normal: the tangent turned by a negative quarter turn. */
inline V2 left_of(V2 t) { return {t.y, -t.x}; }
inline V2 right_of(V2 t) { return {-t.y, t.x}; }

/** The parameters, as the oracle reads them (stroke.h's values; it does not include stroke.h's types). */
struct Style {
  Real width = 1;
  int space = 0;     ///< 0 local, 1 device
  int join = 0;      ///< 0 miter, 1 round, 2 bevel, 3 miter-clip, 4 arcs
  Real limit = 4;    ///< the miter limit
  int cap = 0;       ///< 0 butt, 1 round, 2 square
  int degenerate = 0; ///< 0 draw, 1 omit
  /** A defect planted in the oracle itself, which its self-check must see: 1 puts the round join's sector on
   * the inner side. */
  int plant = 0;
  /** The dash array as given (SG3 is the oracle's own to apply), the offset and the scope (0 per subpath, 1 path). */
  std::vector<Real> dashes;
  Real dash_offset = 0;
  int dash_scope = 0;
};

enum { JOIN_MITER = 0, JOIN_ROUND = 1, JOIN_BEVEL = 2, JOIN_MITER_CLIP = 3, JOIN_ARCS = 4 };
enum { CAP_BUTT = 0, CAP_ROUND = 1, CAP_SQUARE = 2 };

/** An affine map: `x' = a x + b y + tx`, `y' = c x + d y + ty`. */
struct Affine {
  Real a = 1, b = 0, c = 0, d = 1, tx = 0, ty = 0;
  V2 apply(V2 p) const { return {a * p.x + b * p.y + tx, c * p.x + d * p.y + ty}; }
  V2 linear(V2 v) const { return {a * v.x + b * v.y, c * v.x + d * v.y}; }
  Real det() const { return a * d - b * c; }
  V2 inverse(V2 q) const {
    Real dt = det();
    V2 r = {q.x - tx, q.y - ty};
    return {(d * r.x - b * r.y) / dt, (-c * r.x + a * r.y) / dt};
  }
  /** An upper bound of the spectral norm: the Frobenius norm. */
  Real frobenius() const { return std::sqrt(a * a + b * b + c * c + d * d); }
};

inline Affine affine_of(const GSHGEOM_Map & m) {
  Affine f;
  f.a = m.m[0][0];
  f.b = m.m[0][1];
  f.c = m.m[1][0];
  f.d = m.m[1][1];
  f.tx = m.t[0];
  f.ty = m.t[1];
  return f;
}

/** A segment: degree `n` (1 to 3), control points, weights (1 for a polynomial kind). */
struct Seg {
  int n = 1;
  V2 p[4];
  Real w[4] = {1, 1, 1, 1};
  bool rational = false;
};

/** SG5: a segment has extent when a control point differs, exactly, from the first. */
inline bool has_extent(const Seg & s) {
  for (int k = 1; k <= s.n; k++) {
    if (s.p[k].x != s.p[0].x || s.p[k].y != s.p[0].y) {
      return true;
    }
  }
  return false;
}

/** SG5's start direction `P_k - P_0` (least k with P_k != P_0) and the control point that defines it. */
inline V2 start_dir(const Seg & s, V2 * defining = nullptr) {
  for (int k = 1; k <= s.n; k++) {
    if (s.p[k].x != s.p[0].x || s.p[k].y != s.p[0].y) {
      if (defining) {
        *defining = s.p[k];
      }
      return s.p[k] - s.p[0];
    }
  }
  return {0, 0};
}

/** SG5's end direction `P_n - P_k` (greatest k with P_k != P_n). */
inline V2 end_dir(const Seg & s, V2 * defining = nullptr) {
  for (int k = s.n - 1; k >= 0; k--) {
    if (s.p[k].x != s.p[s.n].x || s.p[k].y != s.p[s.n].y) {
      if (defining) {
        *defining = s.p[k];
      }
      return s.p[s.n] - s.p[k];
    }
  }
  return {0, 0};
}

/* --- polynomials in the Bernstein basis -------------------------------------------- */

inline Real binom(int n, int k) {
  Real r = 1;
  for (int i = 1; i <= k; i++) {
    r = r * (n - k + i) / i;
  }
  return r;
}

/** Bernstein polynomial of degree size()-1 on [0, 1]. */
typedef std::vector<Real> Bern;

inline Real bern_eval(const Bern & b, Real t) {
  std::vector<Real> q(b);
  for (size_t r = 1; r < q.size(); r++) {
    for (size_t i = 0; i + r < q.size(); i++) {
      q[i] = (1 - t) * q[i] + t * q[i + 1];
    }
  }
  return q.empty() ? 0 : q[0];
}

/** The product of two Bernstein polynomials. */
inline Bern bern_mul(const Bern & a, const Bern & b) {
  int m = (int)a.size() - 1, n = (int)b.size() - 1;
  Bern out((size_t)(m + n + 1), 0);
  for (int i = 0; i <= m; i++) {
    for (int j = 0; j <= n; j++) {
      out[(size_t)(i + j)] += a[(size_t)i] * b[(size_t)j] * binom(m, i) * binom(n, j) / binom(m + n, i + j);
    }
  }
  return out;
}

/** Raise the degree by one. */
inline Bern bern_elevate(const Bern & a) {
  int n = (int)a.size() - 1;
  Bern out((size_t)(n + 2), 0);
  out[0] = a[0];
  out[(size_t)(n + 1)] = a[(size_t)n];
  for (int i = 1; i <= n; i++) {
    out[(size_t)i] = ((Real)i / (n + 1)) * a[(size_t)(i - 1)] + (1 - (Real)i / (n + 1)) * a[(size_t)i];
  }
  return out;
}

inline Bern bern_add(Bern a, Bern b) {
  while (a.size() < b.size()) {
    a = bern_elevate(a);
  }
  while (b.size() < a.size()) {
    b = bern_elevate(b);
  }
  for (size_t i = 0; i < a.size(); i++) {
    a[i] += b[i];
  }
  return a;
}

inline Bern bern_scale(Bern a, Real s) {
  for (Real & v : a) {
    v *= s;
  }
  return a;
}

inline Bern bern_deriv(const Bern & a) {
  int n = (int)a.size() - 1;
  if (n == 0) {
    return Bern{0};
  }
  Bern out((size_t)n);
  for (int i = 0; i < n; i++) {
    out[(size_t)i] = n * (a[(size_t)(i + 1)] - a[(size_t)i]);
  }
  return out;
}

/** de Casteljau split at s: the left part over [0, s] and the right over [s, 1]. */
inline void bern_split(const Bern & a, Real s, Bern & left, Bern & right) {
  size_t n = a.size();
  std::vector<Real> q(a);
  left.assign(n, 0);
  right.assign(n, 0);
  left[0] = q[0];
  right[n - 1] = q[n - 1];
  for (size_t r = 1; r < n; r++) {
    for (size_t i = 0; i + r < n; i++) {
      q[i] = (1 - s) * q[i] + s * q[i + 1];
    }
    left[r] = q[0];
    right[n - 1 - r] = q[n - 1 - r];
  }
}

inline int sign_changes(const Bern & a, Real eps) {
  int changes = 0, last = 0;
  for (Real v : a) {
    int s = v > eps ? 1 : (v < -eps ? -1 : 0);
    if (s != 0) {
      if (last != 0 && s != last) {
        changes++;
      }
      last = s;
    }
  }
  return changes;
}

/**
 * The roots of `a` in [lo, hi] (a's own [0, 1] mapped there), isolated by sign
 * variations and bisected.  Sets `*ambiguous` when an interval still has two or
 * more variations at the depth limit (a double root: a point on the evolute).
 */
inline void bern_roots_rec(const Bern & a, Real lo, Real hi, int depth, std::vector<Real> & out, bool * ambiguous) {
  Real scale = 0;
  for (Real v : a) {
    scale = std::max(scale, std::fabs(v));
  }
  if (scale == 0) {
    *ambiguous = true; /* identically zero on the interval */
    return;
  }
  Real eps = scale * 1e-30L;
  int v = sign_changes(a, eps);
  if (v == 0) {
    /* a root at either end of the interval (a point on an end rung) */
    if (std::fabs(a.front()) <= eps) {
      out.push_back(lo);
    }
    if (std::fabs(a.back()) <= eps) {
      out.push_back(hi);
    }
    return;
  }
  if (v == 1) {
    /* exactly one root: bisect the polynomial itself */
    Real l = 0, h = 1;
    Real fl = a.front();
    for (int it = 0; it < 90; it++) {
      Real m = (l + h) / 2;
      Real fm = bern_eval(a, m);
      if (fm == 0) {
        l = h = m;
        break;
      }
      if ((fm < 0) == (fl < 0)) {
        l = m;
        fl = fm;
      } else {
        h = m;
      }
    }
    out.push_back(lo + (hi - lo) * (l + h) / 2);
    return;
  }
  if (depth >= 60) {
    *ambiguous = true;
    out.push_back((lo + hi) / 2);
    return;
  }
  Bern l, r;
  bern_split(a, 0.5L, l, r);
  Real mid = (lo + hi) / 2;
  bern_roots_rec(l, lo, mid, depth + 1, out, ambiguous);
  bern_roots_rec(r, mid, hi, depth + 1, out, ambiguous);
}

/* --- a curve as polynomials -------------------------------------------------------- */

/** A segment as Bernstein polynomials: c = (Nx, Ny) / D, and E = N'D - ND'. */
struct Poly {
  Seg seg;
  Bern nx, ny, d;
  Bern ex, ey;   ///< E: the direction of c'(t)
  Bern dex, dey; ///< E'
  Bern ax, ay;   ///< D E
  Bern b;        ///< N . E
  std::vector<Real> cusps;
};

inline Poly make_poly(const Seg & s) {
  Poly p;
  p.seg = s;
  int n = s.n;
  p.nx.resize((size_t)n + 1);
  p.ny.resize((size_t)n + 1);
  p.d.resize((size_t)n + 1);
  for (int i = 0; i <= n; i++) {
    p.nx[(size_t)i] = s.w[i] * s.p[i].x;
    p.ny[(size_t)i] = s.w[i] * s.p[i].y;
    p.d[(size_t)i] = s.w[i];
  }
  Bern dnx = bern_deriv(p.nx), dny = bern_deriv(p.ny), dd = bern_deriv(p.d);
  p.ex = bern_add(bern_mul(dnx, p.d), bern_scale(bern_mul(p.nx, dd), -1));
  p.ey = bern_add(bern_mul(dny, p.d), bern_scale(bern_mul(p.ny, dd), -1));
  p.dex = bern_deriv(p.ex);
  p.dey = bern_deriv(p.ey);
  p.ax = bern_mul(p.d, p.ex);
  p.ay = bern_mul(p.d, p.ey);
  p.b = bern_add(bern_mul(p.nx, p.ex), bern_mul(p.ny, p.ey));
  /* cusps: common roots of Ex and Ey inside (0, 1) */
  Real scale = 0;
  for (Real v : p.ex) {
    scale = std::max(scale, std::fabs(v));
  }
  for (Real v : p.ey) {
    scale = std::max(scale, std::fabs(v));
  }
  bool amb = false;
  std::vector<Real> rx, ry;
  bern_roots_rec(p.ex, 0, 1, 0, rx, &amb);
  amb = false;
  bern_roots_rec(p.ey, 0, 1, 0, ry, &amb);
  const std::vector<Real> & first = rx.empty() ? ry : rx;
  for (Real t : first) {
    if (t <= 0 || t >= 1) {
      continue;
    }
    if (std::fabs(bern_eval(p.ex, t)) <= scale * 1e-12L && std::fabs(bern_eval(p.ey, t)) <= scale * 1e-12L) {
      p.cusps.push_back(t);
    }
  }
  return p;
}

inline V2 poly_point(const Poly & p, Real t) {
  Real dv = bern_eval(p.d, t);
  return {bern_eval(p.nx, t) / dv, bern_eval(p.ny, t) / dv};
}

inline V2 poly_dir(const Poly & p, Real t) { return {bern_eval(p.ex, t), bern_eval(p.ey, t)}; }

/** The segment mapped by an affine map (weights unchanged). */
inline Seg map_seg(const Seg & s, const Affine & f) {
  Seg o = s;
  for (int i = 0; i <= s.n; i++) {
    o.p[i] = f.apply(s.p[i]);
  }
  return o;
}

/* --- the shapes --------------------------------------------------------------------- */

/** A closed convex polygon given counter... in any orientation; membership by sides. */
inline bool in_convex(const std::vector<V2> & poly, V2 x) {
  int pos = 0, neg = 0;
  size_t n = poly.size();
  for (size_t i = 0; i < n; i++) {
    Real c = cross(poly[(i + 1) % n] - poly[i], x - poly[i]);
    if (c > 0) {
      pos++;
    } else if (c < 0) {
      neg++;
    }
  }
  return pos == 0 || neg == 0;
}

/** The angular wedge from unit u0 to unit u1 turning by the short way (angle in [0, pi]), inclusive. */
inline bool in_wedge(V2 u0, V2 u1, V2 v) {
  Real c01 = cross(u0, u1);
  if (c01 >= 0) {
    return cross(u0, v) >= 0 && cross(v, u1) >= 0;
  }
  return cross(u0, v) <= 0 && cross(v, u1) <= 0;
}

/** The region of one join, cap, disc or square, in one space. */
struct Shape {
  enum Kind { POLY, SECTOR, DISC, ARCS } kind = POLY;
  std::vector<V2> poly; ///< POLY: convex polygon
  V2 c;                 ///< SECTOR, DISC: centre
  Real r = 0;           ///< radius
  V2 u0, u1;            ///< SECTOR: the wedge's ends (unit), short way
  /* ARCS: the CR's region */
  V2 V, A, B, X;
  bool line_a = false, line_b = false;
  V2 za, zb;    ///< circle centres
  Real ra = 0, rb = 0;
  V2 ahat, bhat;
  bool clip = false;
  V2 clip_y, clip_t;
  Real span = 0; ///< a size bound, for sampling
};

/** The angle the segment from a to b sweeps as seen from x. */
inline Real seen_angle(V2 x, V2 a, V2 b) {
  V2 p = a - x, q = b - x;
  return std::atan2(cross(p, q), dot(p, q));
}

/**
 * The angle the arc of circle (z, r) sweeps as seen from x, the arc starting at angle a0 and turning by `sweep`
 * (>= 0) in the sense `sense` (+1 counterclockwise in atan2's sense).  Each piece of at most an eighth of a turn
 * sweeps its chord's angle, plus a whole turn in the arc's sense when x lies between the chord and the arc.
 */
inline Real arc_seen_angle(V2 x, V2 z, Real r, Real a0, Real sense, Real sweep) {
  int n = (int)std::ceil(sweep / (kPi / 4));
  if (n < 1) {
    n = 1;
  }
  Real total = 0;
  for (int k = 0; k < n; k++) {
    Real t0 = a0 + sense * sweep * k / n, t1 = a0 + sense * sweep * (k + 1) / n;
    V2 p0 = z + r * V2{std::cos(t0), std::sin(t0)}, p1 = z + r * V2{std::cos(t1), std::sin(t1)};
    total += seen_angle(x, p0, p1);
    if (norm(x - z) < r) {
      Real sx = cross(p1 - p0, x - p0), sz = cross(p1 - p0, z - p0);
      if ((sx > 0 && sz < 0) || (sx < 0 && sz > 0)) {
        total += sense * 2 * kPi;
      }
    }
  }
  return total;
}

/** The start angle, sense and sweep of the arc of a circle centred at z from s, leaving in direction dir, to e. */
inline void arc_extent(V2 z, V2 s, V2 dir, V2 e, Real * a0, Real * sense, Real * sweep) {
  *sense = cross(s - z, dir) >= 0 ? 1 : -1;
  *a0 = std::atan2(s.y - z.y, s.x - z.x);
  Real a1 = std::atan2(e.y - z.y, e.x - z.x);
  Real sw = *sense * (a1 - *a0);
  while (sw < 0) {
    sw += 2 * kPi;
  }
  while (sw >= 2 * kPi) {
    sw -= 2 * kPi;
  }
  *sweep = sw;
}

/**
 * Whether x is in the CR's arcs region before its clip: the region the closed boundary V, A, A's arc (or line)
 * to X, B's arc (or line) back to B, V encloses (nonzero winding).  The boundary can curl past V (a reversal's
 * two arcs run on round a half circle to their far meeting point), so this is a winding number, not a ray from
 * V.
 */
inline bool arcs_region_contains(const Shape & s, V2 x) {
  Real total = seen_angle(x, s.V, s.A);
  Real a0, sense, sweep;
  if (s.line_a) {
    total += seen_angle(x, s.A, s.X);
  } else {
    arc_extent(s.za, s.A, s.ahat, s.X, &a0, &sense, &sweep);
    total += arc_seen_angle(x, s.za, s.ra, a0, sense, sweep);
  }
  if (s.line_b) {
    total += seen_angle(x, s.X, s.B);
  } else {
    arc_extent(s.zb, s.B, (Real)-1 * s.bhat, s.X, &a0, &sense, &sweep);
    /* traversed from X back to B: the same arc the other way */
    Real ax = std::atan2(s.X.y - s.zb.y, s.X.x - s.zb.x);
    total += arc_seen_angle(x, s.zb, s.rb, ax, -sense, sweep);
  }
  total += seen_angle(x, s.B, s.V);
  return std::fabs(total) > kPi;
}

inline bool shape_contains(const Shape & s, V2 x) {
  switch (s.kind) {
  case Shape::POLY:
    return in_convex(s.poly, x);
  case Shape::DISC:
    return norm(x - s.c) <= s.r;
  case Shape::SECTOR:
    return norm(x - s.c) <= s.r && in_wedge(s.u0, s.u1, x - s.c);
  case Shape::ARCS: {
    if (s.clip && dot(x - s.clip_y, s.clip_t) > 0) {
      return false;
    }
    /* the bevel (a triangle with no area at a reversal contains nothing: every point of its line would pass the
     * sign test) */
    if (cross(s.A - s.V, s.B - s.V) != 0 && in_convex({s.V, s.A, s.B}, x)) {
      return true;
    }
    return arcs_region_contains(s, x);
  }
  }
  return false;
}

/* --- the CR's arcs join (SG8) ------------------------------------------------------ */

/** Intersections of two circles. */
inline int circle_circle(V2 z1, Real r1, V2 z2, Real r2, V2 out[2]) {
  V2 dz = z2 - z1;
  Real d = norm(dz);
  if (d == 0 || d > r1 + r2 || d < std::fabs(r1 - r2)) {
    return 0;
  }
  Real a = (r1 * r1 - r2 * r2 + d * d) / (2 * d);
  Real hh = r1 * r1 - a * a;
  if (hh < 0) {
    hh = 0;
  }
  Real hl = std::sqrt(hh);
  V2 e = (1 / d) * dz;
  V2 m = z1 + a * e;
  out[0] = m + hl * V2{-e.y, e.x};
  out[1] = m - hl * V2{-e.y, e.x};
  return hl == 0 ? 1 : 2;
}

inline int circle_line(V2 z, Real r, V2 p, V2 dir, V2 out[2]) {
  V2 d = unit(dir);
  V2 m = p - z;
  Real bq = dot(m, d), cq = dot(m, m) - r * r;
  Real disc = bq * bq - cq;
  if (disc < 0) {
    return 0;
  }
  Real sq = std::sqrt(disc);
  out[0] = p + (-bq - sq) * d;
  out[1] = p + (-bq + sq) * d;
  return sq == 0 ? 1 : 2;
}

/** Least positive root of q2 d^2 + q1 d + q0 = 0; -1 if none. */
inline Real least_positive(Real q2, Real q1, Real q0) {
  Real best = -1;
  auto take = [&](Real r) {
    if (r > 0 && (best < 0 || r < best)) {
      best = r;
    }
  };
  if (std::fabs(q2) < 1e-300L) {
    if (q1 != 0) {
      take(-q0 / q1);
    }
    return best;
  }
  Real disc = q1 * q1 - 4 * q2 * q0;
  if (disc < 0) {
    if (disc > -1e-24L * (q1 * q1 + std::fabs(4 * q2 * q0))) {
      disc = 0;
    } else {
      return -1;
    }
  }
  /* the stable pair: one root from the sum without cancellation, the other from the product */
  Real sq = std::sqrt(disc);
  Real qq = -(q1 + (q1 < 0 ? -sq : sq)) / 2;
  if (qq != 0) {
    take(qq / q2);
    take(q0 / qq);
  } else {
    take(0);
  }
  return best;
}

/**
 * The curvature at the end of a segment that touches the vertex, per the 2018 CR and SG8: `((n-1)/n) (w_0 w_2 /
 * w_1^2) |d1 x d2| / |d1|^3` over the three control points nearest that end, signed by the turning sense (positive:
 * turning toward +y, the centre of curvature on the right of travel).  `*line` when the end's first difference
 * vanishes or the segment is a line.
 */
inline Real end_curvature(const Seg & s, bool at_end, bool * line) {
  *line = false;
  if (s.n == 1) {
    *line = true;
    return 0;
  }
  V2 q0, q1, q2;
  Real w0, w1, w2;
  if (at_end) {
    q0 = s.p[s.n];
    q1 = s.p[s.n - 1];
    q2 = s.p[s.n - 2];
    w0 = s.w[s.n];
    w1 = s.w[s.n - 1];
    w2 = s.w[s.n - 2];
  } else {
    q0 = s.p[0];
    q1 = s.p[1];
    q2 = s.p[2];
    w0 = s.w[0];
    w1 = s.w[1];
    w2 = s.w[2];
  }
  V2 d1 = q1 - q0;
  if (d1.x == 0 && d1.y == 0) {
    *line = true;
    return 0;
  }
  /* turning sense along travel: at the start, travel is q0 -> q1 -> q2; at the end, q2 -> q1 -> q0 */
  Real x = at_end ? cross(q1 - q2, q0 - q1) : cross(q1 - q0, q2 - q1);
  Real k = ((Real)(s.n - 1) / s.n) * (w0 * w2 / (w1 * w1)) * std::fabs(x) / std::pow(norm(d1), 3);
  if (x == 0) {
    return 0;
  }
  return x > 0 ? k : -k;
}

struct JoinInput {
  V2 V;
  V2 a, b;      ///< raw end direction of the incoming segment, start direction of the outgoing
  V2 Q, S;      ///< the control points defining them (a = V - Q, b = S - V)
  const Seg * in = nullptr;
  const Seg * out = nullptr;
};

/** One arcs join with the outer side given (outer_left), per SG8 (the CR). Appends shapes. */
inline void arcs_join(const JoinInput & j, Real h, Real L, bool outer_left, std::vector<Shape> & shapes);

/** A miter-clip join (also the arcs fallback), outer side given; the miter within-limit flag decided by the caller. */
inline void miter_clip_join(V2 V, V2 ah, V2 bh, V2 oa, V2 ob, Real h, Real L, bool within, std::vector<Shape> & shapes) {
  V2 A = V + h * oa, B = V + h * ob;
  V2 u = unit(oa + ob);
  Real cosphi = dot(oa, ob);
  if (within) {
    V2 T = V + (h / (1 + cosphi)) * (oa + ob);
    shapes.push_back(Shape{});
    shapes.back().poly = {V, A, T, B};
    return;
  }
  /* A' and B' on the tangent lines A + l a^, B - m b^ at (x - V) . u = L h */
  Real la = (L * h - dot(A - V, u)) / dot(ah, u);
  Real lb = (L * h - dot(B - V, u)) / dot((Real)-1 * bh, u);
  V2 A2 = A + la * ah, B2 = B - lb * bh;
  shapes.push_back(Shape{});
  shapes.back().poly = {V, A, A2, B2, B};
}

/** The miter limit test of SG9 in the oracle's own form: 1 / sin(theta/2) <= L. */
inline bool miter_within(V2 a, V2 b, Real L) {
  Real na = norm(a), nb = norm(b);
  return (na * nb + dot(a, b)) * L * L >= 2 * na * nb;
}

inline void arcs_join(const JoinInput & j, Real h, Real L, bool outer_left, std::vector<Shape> & shapes) {
  V2 ah = unit(j.a), bh = unit(j.b);
  V2 oa = outer_left ? left_of(ah) : right_of(ah);
  V2 ob = outer_left ? left_of(bh) : right_of(bh);
  V2 V = j.V, A = V + h * oa, B = V + h * ob;
  bool la = false, lb = false;
  Real ka = end_curvature(*j.in, true, &la);
  Real kb = end_curvature(*j.out, false, &lb);
  /* signed positive toward the inner side: with the outer side on the left, the inner side is the right, where a
   * positive turning sense puts the centre */
  if (!outer_left) {
    ka = -ka;
    kb = -kb;
  }
  bool reversal = cross(j.a, j.b) == 0 && dot(j.a, j.b) < 0;
  if (ka == 0 && kb == 0) {
    miter_clip_join(V, ah, bh, oa, ob, h, L, !reversal && miter_within(j.a, j.b, L), shapes);
    if (reversal) {
      shapes.back().poly = {V + h * oa, V + h * oa + (L * h) * ah, V - h * oa + (L * h) * ah, V - h * oa};
    }
    return;
  }
  if (std::fabs(ka) > 1 / h || std::fabs(kb) > 1 / h) {
    Shape s;
    s.kind = reversal ? Shape::POLY : Shape::SECTOR;
    if (reversal) {
      /* the round join's limit: the half disc ahead */
      s.kind = Shape::SECTOR;
      s.c = V;
      s.r = h;
      s.u0 = oa;
      s.u1 = ah;
      shapes.push_back(s);
      s.u0 = ah;
      s.u1 = (Real)-1 * oa;
      shapes.push_back(s);
      return;
    }
    s.c = V;
    s.r = h;
    s.u0 = oa;
    s.u1 = ob;
    shapes.push_back(s);
    return;
  }
  Shape s;
  s.kind = Shape::ARCS;
  s.V = V;
  s.A = A;
  s.B = B;
  s.ahat = ah;
  s.bhat = bh;
  s.line_a = ka == 0;
  s.line_b = kb == 0;
  if (!s.line_a) {
    s.za = V - (1 / ka) * oa;
    s.ra = std::fabs(1 / ka + h);
  }
  if (!s.line_b) {
    s.zb = V - (1 / kb) * ob;
    s.rb = std::fabs(1 / kb + h);
  }
  /* the meeting point closest to V */
  V2 cand[2];
  int nc = 0;
  if (s.line_a && s.line_b) {
    Real den = cross(ah, bh);
    if (den != 0) {
      Real t = cross(B - A, bh) / den;
      cand[0] = A + t * ah;
      nc = 1;
    }
  } else if (s.line_a) {
    nc = circle_line(s.zb, s.rb, A, ah, cand);
  } else if (s.line_b) {
    nc = circle_line(s.za, s.ra, B, bh, cand);
  } else {
    nc = circle_circle(s.za, s.ra, s.zb, s.rb, cand);
  }
  auto nu_of = [&](Real k, V2 o) { return (1 / k + h) > 0 ? (Real)-1 * o : o; };
  if (nc == 2 && s.line_a != s.line_b) {
    /* a line's extension is a ray: forward from A along a, backward from B along -b (SG8) */
    V2 lp = s.line_a ? A : B, ld = s.line_a ? ah : (Real)-1 * bh;
    bool ok0 = dot(cand[0] - lp, ld) >= 0, ok1 = dot(cand[1] - lp, ld) >= 0;
    if (!ok0 && !ok1) {
      nc = 0;
    } else if (!ok0) {
      cand[0] = cand[1];
      nc = 1;
    } else if (!ok1) {
      nc = 1;
    }
  }
  if (nc == 0) {
    /* adjust (the CR's 13.5.9) */
    Real delta = -1;
    if (!s.line_a && !s.line_b) {
      V2 nua = nu_of(ka, oa), nub = nu_of(kb, ob);
      Real dz = norm(s.za - s.zb);
      bool enclosing = dz < std::fabs(s.ra - s.rb);
      bool a_big = s.ra > s.rb;
      if (enclosing) {
        /* big shrinks, small grows */
        V2 g = a_big ? (Real)-1 * (nua + nub) : nua + nub; /* d(Za - Zb)/d delta */
        Real c0 = std::fabs(s.ra - s.rb);
        V2 w0 = s.za - s.zb;
        delta = least_positive(dot(g, g) - 4, 2 * dot(w0, g) + 4 * c0, dot(w0, w0) - c0 * c0);
        if (delta > 0) {
          if (a_big) {
            s.za = s.za - delta * nua;
            s.ra -= delta;
            s.zb = s.zb + delta * nub;
            s.rb += delta;
          } else {
            s.za = s.za + delta * nua;
            s.ra += delta;
            s.zb = s.zb - delta * nub;
            s.rb -= delta;
          }
        }
      } else {
        V2 g = nua - nub;
        Real c0 = s.ra + s.rb;
        V2 w0 = s.za - s.zb;
        delta = least_positive(dot(g, g) - 4, 2 * dot(w0, g) - 4 * c0, dot(w0, w0) - c0 * c0);
        if (delta > 0) {
          s.za = s.za + delta * nua;
          s.ra += delta;
          s.zb = s.zb + delta * nub;
          s.rb += delta;
        }
      }
      if (delta > 0) {
        /* the point of touching: on the line of centres, at the larger (or, apart, either) circle's radius */
        if (enclosing && !a_big) {
          cand[0] = s.zb + s.rb * unit(s.za - s.zb);
        } else {
          cand[0] = s.za + s.ra * unit(s.zb - s.za);
        }
        nc = 1;
      }
    } else if (s.line_a != s.line_b) {
      /* the circle alone grows until it touches the line */
      bool ca = !s.line_a;
      V2 z = ca ? s.za : s.zb;
      Real r = ca ? s.ra : s.rb;
      V2 nu = ca ? nu_of(ka, oa) : nu_of(kb, ob);
      V2 lp = ca ? B : A, ld = ca ? bh : ah;
      V2 ray = ca ? (Real)-1 * bh : ah; /* the extension's direction of travel */
      V2 m = left_of(ld);
      Real d0 = dot(z - lp, m), g = dot(nu, m);
      Real best = -1;
      for (int sgn = -1; sgn <= 1; sgn += 2) {
        Real den = g - sgn;
        if (den != 0) {
          Real dl = (sgn * r - d0) / den;
          if (dl > 0 && (best < 0 || dl < best)) {
            /* the touching point must be on the ray */
            V2 zz = z + dl * nu;
            V2 foot = zz - dot(zz - lp, m) * m;
            if (dot(foot - lp, ray) >= 0) {
              best = dl;
            }
          }
        }
      }
      if (best > 0) {
        z = z + best * nu;
        r += best;
        if (ca) {
          s.za = z;
          s.ra = r;
        } else {
          s.zb = z;
          s.rb = r;
        }
        cand[0] = z - dot(z - lp, m) * m;
        nc = 1;
      }
    }
    if (nc == 0) {
      /* the CR's rectangle: parallel tangents that never meet */
      shapes.push_back(Shape{});
      shapes.back().poly = {V + h * oa, V + h * oa + (L * h) * ah, V - h * oa + (L * h) * ah, V - h * oa};
      return;
    }
  }
  int pick = 0;
  if (nc == 2) {
    /* the closer candidate is the one on V's side of the line of centres (the two are mirror images across it); they
     * tie exactly when V is on it, which is a reversal; a tie goes to the one reached first along A's extension */
    Real side;
    if (!s.line_a && !s.line_b) {
      side = -cross(s.zb - s.za, V - s.za);
    } else {
      V2 z = s.line_a ? s.zb : s.za;
      V2 lp = s.line_a ? A : B, ld = s.line_a ? ah : bh;
      Real bq = dot(lp - z, ld);
      side = bq - dot(lp, ld) + dot(V, ld); /* V's place along the line minus the circle's foot */
    }
    auto along = [&](V2 p) {
      if (s.line_a) {
        return dot(p - A, ah);
      }
      Real sense = cross(A - s.za, ah) >= 0 ? 1 : -1;
      Real ang = sense * (std::atan2(p.y - s.za.y, p.x - s.za.x) - std::atan2(A.y - s.za.y, A.x - s.za.x));
      while (ang < 0) {
        ang += 2 * kPi;
      }
      return ang;
    };
    if (reversal || side == 0) {
      pick = along(cand[1]) < along(cand[0]) ? 1 : 0;
    } else {
      pick = side > 0 ? 1 : 0;
    }
  }
  s.X = cand[pick];
  /* Gamma, the miter length along it, and the clip */
  V2 u = unit(oa + ob);
  if (reversal) {
    u = ah;
  }
  V2 d = s.X - V;
  Real dn = norm(d);
  s.span = std::max(dn, h) + L * h;
  if (dn > 0) {
    Real cr = cross(u, d), dt = dot(u, d);
    Real Lh = L * h;
    if (cr == 0) {
      if (dt > 0) {
        if (dn > Lh) {
          s.clip = true;
          s.clip_y = V + Lh * u;
          s.clip_t = u;
        }
      } else {
        s.clip = true;
        s.clip_y = V + Lh * u;
        s.clip_t = u;
      }
    } else {
      V2 w = cr > 0 ? right_of(u) : left_of(u);
      /* w is the unit normal of u on the side of d */
      if (dot(w, d) < 0) {
        w = (Real)-1 * w;
      }
      Real rho = dn * dn / (2 * dot(d, w));
      Real beta = std::atan2(std::fabs(cr), dt);
      Real sl = 2 * beta * rho;
      if (sl > Lh) {
        Real psi = Lh / rho;
        V2 O = V + rho * w;
        s.clip = true;
        s.clip_y = O + rho * ((Real)-1 * std::cos(psi) * w + std::sin(psi) * u);
        s.clip_t = std::sin(psi) * w + std::cos(psi) * u;
      }
    }
  }
  shapes.push_back(s);
}

/* --- a stroked path ----------------------------------------------------------------- */

/** One subpath as SG5 reads it. */
struct Sub {
  std::vector<Seg> segs; ///< with extent only, in order
  bool closed = false;
  bool draws = false; ///< a drawing verb or a close
  V2 start;
  /** A dash piece (stroke.h, SG22): the joins are given, between the original segments, not between the cut ones. */
  bool use_joins = false;
  std::vector<std::pair<Seg, Seg>> joins;
  /** A zero-length dash (SG24): a point and a tangent; drawn whatever `degenerate` is. */
  bool dot = false;
  V2 dot_dir;
  /** For alignment: the subpath this was cut from is a closed one with extent. */
  bool aligned_closed = false;
};

/** The centerline read from a path (model space), every subpath in order. */
struct Centerline {
  std::vector<Sub> subs;
};

inline Seg seg_of_curve(const GSHGEOM_Curve & c) {
  Seg s;
  s.n = gshgeom_curve_point_count(c.kind) - 1;
  s.rational = gshgeom_curve_is_rational(c.kind) != 0;
  for (int i = 0; i <= s.n; i++) {
    s.p[i] = {c.p[i].x, c.p[i].y};
    s.w[i] = s.rational ? c.w[i] : 1;
  }
  return s;
}

inline Centerline read_centerline(const GSHGEOM_Path * path) {
  Centerline cl;
  size_t n = gshgeom_path_element_count(path);
  bool open = false;
  for (size_t i = 0; i < n; i++) {
    GSHGEOM_Element el;
    gshgeom_path_element(path, i, &el);
    if (el.verb == GSHGEOM_VERB_MOVE) {
      cl.subs.push_back(Sub{});
      cl.subs.back().start = {el.curve.p[0].x, el.curve.p[0].y};
      open = true;
      continue;
    }
    if (!open) {
      /* a drawing verb after a close starts a subpath at the closed one's start */
      V2 st = cl.subs.back().start;
      cl.subs.push_back(Sub{});
      cl.subs.back().start = st;
      open = true;
    }
    Sub & sb = cl.subs.back();
    sb.draws = true;
    Seg s = seg_of_curve(el.curve);
    if (has_extent(s)) {
      sb.segs.push_back(s);
    }
    if (el.verb == GSHGEOM_VERB_CLOSE) {
      sb.closed = true;
      open = false;
    }
  }
  return cl;
}

/** Everything one part covers, in one space (model under local width, device under device width). */
struct Region {
  Affine to_device;   ///< identity under device width
  bool local = true;
  Real h = 0;
  std::vector<Poly> bodies;
  std::vector<Shape> shapes;
  Real device_scale = 1; ///< an upper bound of how far a unit of this space reaches in device pixels
};

inline void add_cap(Region & r, V2 E, V2 d, int cap) {
  Real h = r.h;
  V2 dh = unit(d), m = left_of(dh);
  if (cap == CAP_SQUARE) {
    Shape s;
    s.poly = {E + h * m, E + h * m + h * dh, E - h * m + h * dh, E - h * m};
    r.shapes.push_back(s);
  } else if (cap == CAP_ROUND) {
    Shape s;
    s.kind = Shape::SECTOR;
    s.c = E;
    s.r = h;
    s.u0 = m;
    s.u1 = dh;
    r.shapes.push_back(s);
    s.u0 = dh;
    s.u1 = (Real)-1 * m;
    r.shapes.push_back(s);
  }
}

inline void add_join(Region & r, const Seg & in, const Seg & out, int join, Real L, int plant) {
  JoinInput j;
  j.in = &in;
  j.out = &out;
  j.V = in.p[in.n];
  j.a = end_dir(in, &j.Q);
  j.b = start_dir(out, &j.S);
  Real h = r.h;
  Real cr = cross(j.a, j.b), dt = dot(j.a, j.b);
  if (cr == 0 && dt > 0) {
    return; /* smooth */
  }
  V2 ah = unit(j.a), bh = unit(j.b);
  V2 V = j.V;
  if (cr == 0) {
    /* reversal: each join's limit */
    V2 n = left_of(ah);
    if (join == JOIN_ROUND) {
      Shape s;
      s.kind = Shape::SECTOR;
      s.c = V;
      s.r = h;
      s.u0 = n;
      s.u1 = ah;
      r.shapes.push_back(s);
      s.u0 = ah;
      s.u1 = (Real)-1 * n;
      r.shapes.push_back(s);
    } else if (join == JOIN_MITER_CLIP) {
      Shape s;
      s.poly = {V + h * n, V + h * n + (L * h) * ah, V - h * n + (L * h) * ah, V - h * n};
      r.shapes.push_back(s);
    } else if (join == JOIN_ARCS) {
      arcs_join(j, h, L, true, r.shapes);
      arcs_join(j, h, L, false, r.shapes);
    }
    return;
  }
  bool outer_left = cr > 0; /* sigma > 0: the outer side is the left */
  V2 oa = outer_left ? left_of(ah) : right_of(ah);
  V2 ob = outer_left ? left_of(bh) : right_of(bh);
  V2 A = V + h * oa, B = V + h * ob;
  Shape bevel;
  bevel.poly = {V, A, B};
  r.shapes.push_back(bevel);
  switch (join) {
  case JOIN_BEVEL:
    break;
  case JOIN_MITER:
    if (miter_within(j.a, j.b, L)) {
      miter_clip_join(V, ah, bh, oa, ob, h, L, true, r.shapes);
    }
    break;
  case JOIN_MITER_CLIP:
    miter_clip_join(V, ah, bh, oa, ob, h, L, miter_within(j.a, j.b, L), r.shapes);
    break;
  case JOIN_ROUND: {
    Shape s;
    s.kind = Shape::SECTOR;
    s.c = V;
    s.r = h;
    s.u0 = plant == 1 ? (Real)-1 * oa : oa;
    s.u1 = plant == 1 ? (Real)-1 * ob : ob;
    r.shapes.push_back(s);
    break;
  }
  case JOIN_ARCS:
    arcs_join(j, h, L, outer_left, r.shapes);
    break;
  }
}

/**
 * Build the region of every subpath selected by `take` at half width h.  `d_tangent` gives the degenerate square's
 * axis for subpath i (in the region's space).
 */
template <typename Take>
inline void build_region(Region & r, const std::vector<Sub> & subs, const Style & st, Real h, Take take,
    const std::vector<V2> & d_tangent) {
  r.h = h;
  for (size_t i = 0; i < subs.size(); i++) {
    if (!take(i)) {
      continue;
    }
    const Sub & sb = subs[i];
    if (sb.dot) {
      V2 P = sb.start;
      V2 d = unit(sb.dot_dir);
      V2 m = left_of(d);
      if (st.cap == CAP_ROUND) {
        Shape s;
        s.kind = Shape::DISC;
        s.c = P;
        s.r = h;
        r.shapes.push_back(s);
      } else if (st.cap == CAP_SQUARE) {
        Shape s;
        s.poly = {P + h * m + h * d, P - h * m + h * d, P - h * m - h * d, P + h * m - h * d};
        r.shapes.push_back(s);
      }
      continue;
    }
    if (sb.segs.empty()) {
      if (!sb.draws || st.degenerate == 1) {
        continue;
      }
      V2 P = sb.start;
      V2 d = unit(d_tangent[i]);
      V2 m = left_of(d);
      if (st.cap == CAP_ROUND) {
        Shape s;
        s.kind = Shape::DISC;
        s.c = P;
        s.r = h;
        r.shapes.push_back(s);
      } else if (st.cap == CAP_SQUARE) {
        Shape s;
        s.poly = {P + h * m + h * d, P - h * m + h * d, P - h * m - h * d, P + h * m - h * d};
        r.shapes.push_back(s);
      }
      continue;
    }
    for (const Seg & s : sb.segs) {
      r.bodies.push_back(make_poly(s));
    }
    size_t n = sb.segs.size();
    if (sb.use_joins) {
      for (const auto & jp : sb.joins) {
        add_join(r, jp.first, jp.second, st.join, st.limit, st.plant);
      }
    }
    for (size_t k = 0; !sb.use_joins && k + 1 < n; k++) {
      add_join(r, sb.segs[k], sb.segs[k + 1], st.join, st.limit, st.plant);
    }
    if (sb.closed) {
      add_join(r, sb.segs[n - 1], sb.segs[0], st.join, st.limit, st.plant);
    } else {
      V2 sd = start_dir(sb.segs[0]);
      add_cap(r, sb.segs[0].p[0], (Real)-1 * sd, st.cap);
      add_cap(r, sb.segs[n - 1].p[sb.segs[n - 1].n], end_dir(sb.segs[n - 1]), st.cap);
    }
  }
}

/** The result of one membership query. */
struct Verdict {
  bool inside = false;
  bool ambiguous = false; ///< a root the isolation could not separate (a point on the evolute)
};

/** Whether x (in the region's space) is on one of the body's normal segments, or in one of its cusp discs. */
inline bool body_contains(const Poly & p, Real h, V2 x, bool * ambiguous) {
  /* quick reject: farther than h from the control hull's box */
  Real mnx = p.seg.p[0].x, mxx = mnx, mny = p.seg.p[0].y, mxy = mny;
  for (int i = 1; i <= p.seg.n; i++) {
    mnx = std::min(mnx, p.seg.p[i].x);
    mxx = std::max(mxx, p.seg.p[i].x);
    mny = std::min(mny, p.seg.p[i].y);
    mxy = std::max(mxy, p.seg.p[i].y);
  }
  if (x.x < mnx - h || x.x > mxx + h || x.y < mny - h || x.y > mxy + h) {
    return false;
  }
  for (Real t : p.cusps) {
    if (norm(x - poly_point(p, t)) <= h) {
      return true;
    }
  }
  /* f(t) = x . (D E) - N . E */
  Bern f = bern_add(bern_add(bern_scale(p.ax, x.x), bern_scale(p.ay, x.y)), bern_scale(p.b, -1));
  std::vector<Real> roots;
  bool amb = false;
  bern_roots_rec(f, 0, 1, 0, roots, &amb);
  if (amb) {
    *ambiguous = true;
  }
  for (Real t : roots) {
    V2 c = poly_point(p, t);
    if (norm(x - c) <= h) {
      V2 e = poly_dir(p, t);
      if (e.x != 0 || e.y != 0) {
        return true;
      }
    }
  }
  return false;
}

inline Verdict region_contains(const Region & r, V2 x_device) {
  Verdict v;
  V2 x = r.local ? r.to_device.inverse(x_device) : x_device;
  for (const Shape & s : r.shapes) {
    if (shape_contains(s, x)) {
      v.inside = true;
      return v;
    }
  }
  for (const Poly & p : r.bodies) {
    if (body_contains(p, r.h, x, &v.ambiguous)) {
      v.inside = true;
      return v;
    }
  }
  return v;
}

/* --- dashes: the definition of a dashed stroke (stroke.h, SG19 to SG24) ------------------- */

/** A segment's arc length by composite Gauss-Legendre (8 points) over 1024 panels: independent of along.h. */
struct Arc {
  Poly p;
  std::vector<Real> cum; ///< the length at each panel boundary
  Real total = 0;
};

inline Real arc_speed(const Poly & p, Real t) {
  Real dv = bern_eval(p.d, t);
  return norm(poly_dir(p, t)) / (dv * dv);
}

inline Real gl8(const Poly & p, Real a, Real b) {
  static const Real x[4] = {0.18343464249564980493947614236018L, 0.52553240991632898581773904918925L,
      0.79666647741362673959155393647583L, 0.96028985649753623168356086856947L};
  static const Real w[4] = {0.36268378337836198296515044927720L, 0.31370664587788728733796220198660L,
      0.22238103445337447054435599442624L, 0.10122853629037625915253135430996L};
  Real c = (a + b) / 2, h = (b - a) / 2, sum = 0;
  for (int i = 0; i < 4; i++) {
    sum += w[i] * (arc_speed(p, c - h * x[i]) + arc_speed(p, c + h * x[i]));
  }
  return sum * h;
}

inline Arc make_arc(const Seg & s) {
  Arc a;
  a.p = make_poly(s);
  const int panels = 1024;
  a.cum.assign(panels + 1, 0);
  for (int i = 0; i < panels; i++) {
    a.cum[(size_t)i + 1] = a.cum[(size_t)i] + gl8(a.p, (Real)i / panels, (Real)(i + 1) / panels);
  }
  a.total = a.cum[panels];
  return a;
}

/** The length from 0 to t. */
inline Real arc_to(const Arc & a, Real t) {
  const int panels = 1024;
  int i = (int)std::floor(t * panels);
  if (i >= panels) {
    return a.total;
  }
  return a.cum[(size_t)i] + gl8(a.p, (Real)i / panels, t);
}

/** The parameter at length s (0 <= s <= total), by bisection. */
inline Real arc_param(const Arc & a, Real s) {
  Real lo = 0, hi = 1;
  for (int i = 0; i < 80; i++) {
    Real mid = (lo + hi) / 2;
    if (arc_to(a, mid) < s) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  return (lo + hi) / 2;
}

/** The piece of a segment on [a, b] (homogeneous de Casteljau, twice). */
inline void seg_halves(const Seg & s, Real t, Seg & l, Seg & r) {
  int n = s.n;
  Real X[4], Y[4], W[4];
  for (int i = 0; i <= n; i++) {
    W[i] = s.w[i];
    X[i] = s.w[i] * s.p[i].x;
    Y[i] = s.w[i] * s.p[i].y;
  }
  l = s;
  r = s;
  for (int lev = 1; lev <= n; lev++) {
    for (int i = 0; i + lev <= n; i++) {
      X[i] = (1 - t) * X[i] + t * X[i + 1];
      Y[i] = (1 - t) * Y[i] + t * Y[i + 1];
      W[i] = (1 - t) * W[i] + t * W[i + 1];
    }
    l.p[lev] = {X[0] / W[0], Y[0] / W[0]};
    l.w[lev] = W[0];
    r.p[n - lev] = {X[n - lev] / W[n - lev], Y[n - lev] / W[n - lev]};
    r.w[n - lev] = W[n - lev];
  }
}

inline Seg seg_between(const Seg & s, Real a, Real b) {
  Seg cur = s, l, r;
  if (a > 0) {
    seg_halves(cur, a, l, r);
    cur = r;
  }
  if (b < 1) {
    Real u = a > 0 ? (b - a) / (1 - a) : b;
    seg_halves(cur, u, l, r);
    cur = l;
  }
  return cur;
}

/** The unit tangent at t: SG5's directions at the ends, the hodograph inside. */
inline V2 seg_tangent(const Seg & s, Real t) {
  if (t <= 0) {
    return start_dir(s);
  }
  if (t >= 1) {
    return end_dir(s);
  }
  Poly p = make_poly(s);
  V2 e = poly_dir(p, t);
  if (e.x != 0 || e.y != 0) {
    return e;
  }
  return start_dir(seg_between(s, t, 1));
}

/**
 * SVG 2's dash positions, transcribed from the text with a running position (13.5.7), in `long double`.
 */
inline std::vector<std::pair<Real, Real>> dash_positions(const std::vector<Real> & raw, Real offset, Real length) {
  bool solid = false, any = false;
  for (Real v : raw) {
    solid = solid || v < 0;
    any = any || v != 0;
  }
  if (raw.empty() || solid || !any) {
    return {{0, length}};
  }
  std::vector<Real> d = raw;
  if (d.size() % 2) {
    d.insert(d.end(), raw.begin(), raw.end());
  }
  Real sum = 0;
  for (Real v : d) {
    sum += v;
  }
  if (offset < 0) {
    offset = sum - std::fabs(offset);
  }
  offset = std::fmod(offset, sum);
  if (offset < 0) {
    offset += sum;
  }
  std::vector<std::pair<Real, Real>> out;
  size_t index = 0;
  Real acc = d[0];
  while (acc < offset) {
    index++;
    acc += d[index];
  }
  Real remaining = acc - offset;
  Real dashlength = std::min(remaining, length);
  if (index % 2 == 0) {
    out.push_back({0, dashlength});
  }
  Real position = dashlength;
  while (position < length) {
    index = (index + 1) % d.size();
    dashlength = std::min(d[index], length - position);
    if (index % 2 == 0) {
      out.push_back({position, position + dashlength});
    }
    position += dashlength;
  }
  return out;
}

/** A subpath's segments located by length: the segment and parameter at length s (a boundary is the later segment's). */
struct Located {
  size_t k = 0;
  Real t = 0;
};

/** The dashed pieces of the subpaths (stroke.h, SG21 to SG24), as subpaths of their own. */
inline std::vector<Sub> dash_subs(const std::vector<Sub> & subs, const Style & st, const std::vector<V2> & tangents) {
  std::vector<Sub> out;
  Real path_pos = 0;
  for (size_t i = 0; i < subs.size(); i++) {
    const Sub & sb = subs[i];
    if (!sb.draws) {
      continue;
    }
    std::vector<Arc> arcs;
    std::vector<Real> start;
    Real total = 0;
    for (const Seg & g : sb.segs) {
      arcs.push_back(make_arc(g));
      start.push_back(total);
      total += arcs.back().total;
    }
    Real phase = st.dash_scope == 1 ? st.dash_offset + path_pos : st.dash_offset;
    std::vector<std::pair<Real, Real>> pairs = dash_positions(st.dashes, phase, total);
    bool closed_ext = sb.closed && !sb.segs.empty();
    auto locate = [&](Real s, bool end_side) {
      Located l;
      if (s >= total) {
        l.k = sb.segs.size() - 1;
        l.t = 1;
        return l;
      }
      size_t k = 0;
      while (k + 1 < sb.segs.size() && start[k + 1] <= s) {
        k++;
      }
      Real local = s - start[k];
      l.k = k;
      l.t = local <= 0 ? 0 : arc_param(arcs[k], local);
      if (end_side && l.t <= 0 && k > 0) {
        l.k = k - 1;
        l.t = 1;
      }
      return l;
    };
    if (sb.segs.empty()) {
      /* a degenerate subpath: the first piece only, a dot when it exists (unless degenerate omit) */
      if (st.degenerate != 1 && !pairs.empty()) {
        Sub d;
        d.draws = true;
        d.dot = true;
        d.start = sb.start;
        d.dot_dir = tangents[i];
        out.push_back(d);
      }
      path_pos += total;
      continue;
    }
    for (const auto & pr : pairs) {
      Real a = std::min(pr.first, total), b = std::min(pr.second, total);
      Located la = locate(a, false);
      Sub d;
      d.draws = true;
      d.aligned_closed = closed_ext;
      if (!(b > a)) {
        d.dot = true;
        d.start = seg_between(sb.segs[la.k], la.t, 1).p[0];
        d.dot_dir = seg_tangent(sb.segs[la.k], la.t);
        out.push_back(d);
        continue;
      }
      Located lb = locate(b, true);
      if (lb.k < la.k || (lb.k == la.k && lb.t <= la.t)) {
        d.dot = true;
        d.start = seg_between(sb.segs[la.k], la.t, 1).p[0];
        d.dot_dir = seg_tangent(sb.segs[la.k], la.t);
        out.push_back(d);
        continue;
      }
      for (size_t k = la.k; k <= lb.k; k++) {
        Real ta = k == la.k ? la.t : 0, tb = k == lb.k ? lb.t : 1;
        if (!(tb > ta)) {
          continue;
        }
        Seg piece = seg_between(sb.segs[k], ta, tb);
        if (!has_extent(piece)) {
          continue;
        }
        if (!d.segs.empty() && k > 0) {
          d.joins.push_back({sb.segs[k - 1], sb.segs[k]});
        }
        d.segs.push_back(piece);
      }
      if (d.segs.empty()) {
        d.dot = true;
        d.start = seg_between(sb.segs[la.k], la.t, 1).p[0];
        d.dot_dir = seg_tangent(sb.segs[la.k], la.t);
      } else {
        d.use_joins = true;
        /* the caps' directions: the centerline's at the two cut points */
        d.start = d.segs.front().p[0];
      }
      out.push_back(d);
    }
    path_pos += total;
  }
  return out;
}

/** The two regions of SG12: the body (R, or R_o off center) and the aligned part (R_c at twice the width). */
struct Stroked {
  Region body, aligned;
  bool has_aligned = false;
};

/** SG11's effective tangent of every degenerate subpath (in the space the subpaths are given in). */
inline std::vector<V2> degenerate_tangents(const std::vector<Sub> & subs, V2 fallback) {
  std::vector<V2> out(subs.size(), fallback);
  for (size_t i = 0; i < subs.size(); i++) {
    if (!subs[i].segs.empty()) {
      continue;
    }
    bool found = false;
    for (size_t j = i; j-- > 0;) {
      if (!subs[j].segs.empty()) {
        out[i] = end_dir(subs[j].segs.back());
        found = true;
        break;
      }
    }
    for (size_t j = i + 1; !found && j < subs.size(); j++) {
      if (!subs[j].segs.empty()) {
        out[i] = start_dir(subs[j].segs.front());
        found = true;
      }
    }
  }
  return out;
}

/**
 * The regions of a stroke (SG4 to SG13).  `alignment` 0 center, 1 inside, 2 outside.  Under local width the
 * regions are built in model space with `to_device` the map; under device width the centerline is mapped first.
 */
inline Stroked build_stroke(const Centerline & cl, const Style & st, const Affine & map, int alignment) {
  Stroked out;
  std::vector<Sub> subs = cl.subs;
  bool local = st.space == 0;
  V2 fallback = {1, 0};
  if (!local) {
    for (Sub & s : subs) {
      s.start = map.apply(s.start);
      std::vector<Seg> kept;
      for (Seg & g : s.segs) {
        Seg m = map_seg(g, map);
        if (has_extent(m)) {
          kept.push_back(m);
        }
      }
      s.segs = kept;
    }
    fallback = map.linear({1, 0});
  }
  std::vector<V2> tangents = degenerate_tangents(subs, fallback);
  bool any_closed_extent = false;
  for (const Sub & sb : subs) {
    any_closed_extent = any_closed_extent || (sb.closed && !sb.segs.empty());
  }
  bool dashed = false;
  {
    bool neg = false, nz = false;
    for (Real v : st.dashes) {
      neg = neg || v < 0;
      nz = nz || v != 0;
    }
    dashed = !st.dashes.empty() && !neg && nz;
  }
  if (dashed) {
    subs = dash_subs(subs, st, tangents);
    tangents.assign(subs.size(), fallback);
  }
  Real h = st.width / 2;
  for (Region * r : {&out.body, &out.aligned}) {
    r->local = local;
    r->to_device = local ? map : Affine{};
    r->device_scale = local ? map.frobenius() : 1;
  }
  auto closed_with_extent = [&](size_t i) {
    return subs[i].aligned_closed || (!dashed && subs[i].closed && !subs[i].segs.empty());
  };
  if (alignment == 0) {
    build_region(out.body, subs, st, h, [](size_t) { return true; }, tangents);
  } else {
    build_region(out.body, subs, st, h, [&](size_t i) { return !closed_with_extent(i); }, tangents);
    out.has_aligned = any_closed_extent;
    build_region(out.aligned, subs, st, 2 * h, closed_with_extent, tangents);
  }
  return out;
}

/** The distance from x (in the curve's space) to a segment: the self-check's other formulation. */
inline Real distance_to(const Poly & p, V2 x) {
  Real best = std::min(norm(x - p.seg.p[0]), norm(x - p.seg.p[p.seg.n]));
  Bern f = bern_add(bern_add(bern_scale(p.ax, x.x), bern_scale(p.ay, x.y)), bern_scale(p.b, -1));
  std::vector<Real> roots;
  bool amb = false;
  bern_roots_rec(f, 0, 1, 0, roots, &amb);
  for (Real t : roots) {
    best = std::min(best, norm(x - poly_point(p, t)));
  }
  /* a double root the isolation reported once: also try a dense scan near it */
  if (amb) {
    for (int i = 0; i <= 4096; i++) {
      best = std::min(best, norm(x - poly_point(p, (Real)i / 4096)));
    }
  }
  return best;
}

/* --- the boundary samples ------------------------------------------------------------ */

/** A spatial hash of device-space boundary samples, for "farther than d from every sample". */
struct Samples {
  Real cell = 1;
  std::unordered_map<int64_t, std::vector<V2>> grid;
  size_t count = 0;
  static int64_t key(int64_t i, int64_t j) { return (i << 32) ^ (j & 0xffffffffLL); }
  void add(V2 p) {
    int64_t i = (int64_t)std::floor(p.x / cell), j = (int64_t)std::floor(p.y / cell);
    grid[key(i, j)].push_back(p);
    count++;
  }
  /** The distance to the nearest sample, or `limit` if none is nearer than `limit` (limit <= cell). */
  Real nearest(V2 p, Real limit) const {
    int64_t i = (int64_t)std::floor(p.x / cell), j = (int64_t)std::floor(p.y / cell);
    Real best = limit;
    for (int64_t di = -1; di <= 1; di++) {
      for (int64_t dj = -1; dj <= 1; dj++) {
        auto it = grid.find(key(i + di, j + dj));
        if (it == grid.end()) {
          continue;
        }
        for (const V2 & q : it->second) {
          best = std::min(best, norm(p - q));
        }
      }
    }
    return best;
  }
};

/** Sample a curve given by f on [t0, t1] so that consecutive device samples are at most sigma apart. */
template <typename F>
inline void sample_curve(Samples & out, F f, Real t0, Real t1, Real sigma, int depth = 0) {
  V2 a = f(t0), b = f(t1);
  if (!std::isfinite((double)a.x) || !std::isfinite((double)a.y) || !std::isfinite((double)b.x) ||
      !std::isfinite((double)b.y)) {
    abort(); /* a sampled curve must be finite everywhere: a bug of the sampler, not a long run */
  }
  if (depth == 0) {
    out.add(a);
  }
  V2 m = f((t0 + t1) / 2);
  if ((norm(b - a) <= sigma && norm(m - a) <= sigma && norm(m - b) <= sigma) || depth > 40) {
    out.add(b);
    return;
  }
  sample_curve(out, f, t0, (t0 + t1) / 2, sigma, depth + 1);
  sample_curve(out, f, (t0 + t1) / 2, t1, sigma, depth + 1);
}

inline void sample_segment(Samples & out, V2 a, V2 b, Real sigma) {
  int n = (int)std::ceil(norm(b - a) / sigma) + 1;
  for (int i = 0; i <= n; i++) {
    out.add(a + ((Real)i / n) * (b - a));
  }
}

/** The unit direction of a segment's tangent at t: at a zero (a stopped end, a cusp) the one-sided limit. */
inline V2 tangent_at(const Poly & p, Real t) {
  V2 d = poly_dir(p, t);
  if (d.x != 0 || d.y != 0) {
    return unit(d);
  }
  for (Real step : {(Real)1e-9L, (Real)1e-6L, (Real)1e-3L}) {
    Real u = t + step <= 1 ? t + step : t - step;
    V2 e = poly_dir(p, u);
    if (e.x != 0 || e.y != 0) {
      return unit(e);
    }
  }
  return V2{1, 0};
}

inline void sample_region(Samples & out, const Region & r, Real sigma_device) {
  auto dev = [&](V2 p) { return r.local ? r.to_device.apply(p) : p; };
  Real h = r.h;
  for (const Poly & p : r.bodies) {
    for (int side = -1; side <= 1; side += 2) {
      sample_curve(out, [&](Real t) {
        V2 n = left_of(tangent_at(p, t));
        return dev(poly_point(p, t) + ((Real)side * h) * n);
      }, 0, 1, sigma_device);
    }
    for (Real t : {(Real)0, (Real)1}) {
      V2 c = poly_point(p, t), n = left_of(tangent_at(p, t));
      if (t == 0) {
        V2 sd = unit(start_dir(p.seg));
        n = left_of(sd);
      } else {
        n = left_of(unit(end_dir(p.seg)));
      }
      sample_curve(out, [&](Real v) { return dev(c + v * n); }, -h, h, sigma_device);
    }
    /* the evolute within reach: c + rho n where |rho| <= h; sampled where defined */
    sample_curve(out, [&](Real t) {
      V2 e = poly_dir(p, t), de = {bern_eval(p.dex, t), bern_eval(p.dey, t)};
      Real dv = bern_eval(p.d, t);
      Real num = cross(e, de);
      Real ne = norm(e);
      V2 c = poly_point(p, t);
      /* no curvature to speak of (a straight stretch, or a stopped end where the tangent is still): the evolute is
       * beyond any reach, and a rounding residue must not pick a side at random */
      Real nde = norm(de);
      if (ne == 0 || nde == 0 || std::fabs(num) <= 1e-9L * ne * nde) {
        return dev(c);
      }
      /* kappa = (E x E') D^2 / |E|^3, signed toward the left of travel */
      Real kappa = num * dv * dv / (ne * ne * ne);
      Real rho = -1 / kappa; /* the centre lies on the left normal's side when kappa < 0 */
      /* beyond reach the clamped point runs along a rail, so the sampled curve stays continuous */
      rho = std::max(-h, std::min(h, rho));
      return dev(c + rho * left_of(unit(e)));
    }, 0, 1, sigma_device);
    for (Real t : p.cusps) {
      V2 c = poly_point(p, t);
      sample_curve(out, [&](Real a) { return dev(c + h * V2{std::cos(a), std::sin(a)}); }, 0, 2 * kPi,
          sigma_device);
    }
  }
  for (const Shape & s : r.shapes) {
    switch (s.kind) {
    case Shape::POLY:
      for (size_t i = 0; i < s.poly.size(); i++) {
        V2 a = s.poly[i], b = s.poly[(i + 1) % s.poly.size()];
        sample_curve(out, [&](Real t) { return dev(a + t * (b - a)); }, 0, 1, sigma_device);
      }
      break;
    case Shape::DISC:
      sample_curve(out, [&](Real a) { return dev(s.c + s.r * V2{std::cos(a), std::sin(a)}); }, 0, 2 * kPi,
          sigma_device);
      break;
    case Shape::SECTOR: {
      Real a0 = std::atan2(s.u0.y, s.u0.x), a1 = std::atan2(s.u1.y, s.u1.x);
      Real sweep = a1 - a0;
      if (cross(s.u0, s.u1) >= 0) {
        while (sweep < 0) {
          sweep += 2 * kPi;
        }
      } else {
        while (sweep > 0) {
          sweep -= 2 * kPi;
        }
      }
      sample_curve(out, [&](Real t) { return dev(s.c + s.r * V2{std::cos(a0 + t * sweep), std::sin(a0 + t * sweep)}); },
          0, 1, sigma_device);
      sample_curve(out, [&](Real t) { return dev(s.c + t * s.r * s.u0); }, 0, 1, sigma_device);
      sample_curve(out, [&](Real t) { return dev(s.c + t * s.r * s.u1); }, 0, 1, sigma_device);
      break;
    }
    case Shape::ARCS: {
      auto arc = [&](V2 z, Real rr, V2 from, V2 dir, V2 to) {
        Real sense = cross(from - z, dir) >= 0 ? 1 : -1;
        Real a0 = std::atan2(from.y - z.y, from.x - z.x);
        Real a1 = std::atan2(to.y - z.y, to.x - z.x);
        Real sweep = sense * (a1 - a0);
        while (sweep < 0) {
          sweep += 2 * kPi;
        }
        sample_curve(out, [&](Real t) {
          Real a = a0 + sense * t * sweep;
          return dev(z + rr * V2{std::cos(a), std::sin(a)});
        }, 0, 1, sigma_device);
      };
      if (s.line_a) {
        sample_curve(out, [&](Real t) { return dev(s.A + t * (s.X - s.A)); }, 0, 1, sigma_device);
      } else {
        arc(s.za, s.ra, s.A, s.ahat, s.X);
      }
      if (s.line_b) {
        sample_curve(out, [&](Real t) { return dev(s.B + t * (s.X - s.B)); }, 0, 1, sigma_device);
      } else {
        arc(s.zb, s.rb, s.B, (Real)-1 * s.bhat, s.X);
      }
      /* the boundary closes along V -> A and B -> V (the bevel's chord and the rays to X are inside the region) */
      for (V2 e : {s.A, s.B}) {
        sample_curve(out, [&](Real t) { return dev(s.V + t * (e - s.V)); }, 0, 1, sigma_device);
      }
      if (s.clip) {
        V2 perp = left_of(s.clip_t);
        Real reach = 2 * s.span + 2 * h;
        sample_curve(out, [&](Real t) { return dev(s.clip_y + t * perp); }, -reach, reach, sigma_device);
      }
      break;
    }
    }
  }
}

/* --- the outline's nonzero fill ----------------------------------------------------- */

/** The polygons of a flattened path, every subpath closed, in 26.6 units. */
struct Polygons {
  std::vector<std::vector<std::pair<int64_t, int64_t>>> rings;
};

inline Polygons polygons_of(const GSHGEOM_Flat_Path * fp) {
  Polygons out;
  if (!fp) {
    return out;
  }
  const GSHGEOM_Device_Point * pts = gshgeom_flat_path_points(fp);
  size_t ns = gshgeom_flat_path_subpath_count(fp);
  for (size_t i = 0; i < ns; i++) {
    GSHGEOM_Flat_Subpath sp;
    gshgeom_flat_path_subpath(fp, i, &sp);
    std::vector<std::pair<int64_t, int64_t>> ring;
    for (size_t k = 0; k < sp.count; k++) {
      ring.push_back({pts[sp.first + k].x, pts[sp.first + k].y});
    }
    out.rings.push_back(ring);
  }
  return out;
}

/** The winding number of one ring around x (pixels), by crossings of the horizontal ray to +x. */
inline int ring_winding(const std::vector<std::pair<int64_t, int64_t>> & ring, V2 x) {
  int w = 0;
  size_t n = ring.size();
  Real px = x.x * 64, py = x.y * 64;
  for (size_t i = 0; i < n; i++) {
    Real ax = (Real)ring[i].first, ay = (Real)ring[i].second;
    Real bx = (Real)ring[(i + 1) % n].first, by = (Real)ring[(i + 1) % n].second;
    if (ay <= py) {
      if (by > py) {
        if ((bx - ax) * (py - ay) - (px - ax) * (by - ay) > 0) {
          w++;
        }
      }
    } else if (by <= py) {
      if ((bx - ax) * (py - ay) - (px - ax) * (by - ay) < 0) {
        w--;
      }
    }
  }
  return w;
}

/** Polygons bucketed by row band, for many queries. */
struct Winding {
  Polygons polys;
  struct Edge {
    Real ax, ay, bx, by;
  };
  Real band = 64; /* units */
  int64_t lo = 0;
  std::vector<std::vector<Edge>> rows;
  void build(const Polygons & p) {
    polys = p;
    int64_t mn = INT64_MAX, mx = INT64_MIN;
    for (auto & r : p.rings) {
      for (auto & q : r) {
        mn = std::min(mn, q.second);
        mx = std::max(mx, q.second);
      }
    }
    if (mn > mx) {
      return;
    }
    lo = (int64_t)std::floor(mn / band);
    int64_t hi = (int64_t)std::floor(mx / band);
    rows.assign((size_t)(hi - lo + 1), {});
    for (auto & r : p.rings) {
      size_t n = r.size();
      for (size_t i = 0; i < n; i++) {
        Edge e = {(Real)r[i].first, (Real)r[i].second, (Real)r[(i + 1) % n].first, (Real)r[(i + 1) % n].second};
        int64_t a = (int64_t)std::floor(std::min(e.ay, e.by) / band) - lo;
        int64_t b = (int64_t)std::floor(std::max(e.ay, e.by) / band) - lo;
        for (int64_t k = a; k <= b; k++) {
          rows[(size_t)k].push_back(e);
        }
      }
    }
  }
  int at(V2 x) const {
    Real px = x.x * 64, py = x.y * 64;
    int64_t k = (int64_t)std::floor(py / band) - lo;
    if (k < 0 || (size_t)k >= rows.size()) {
      return 0;
    }
    int w = 0;
    for (const Edge & e : rows[(size_t)k]) {
      if (e.ay <= py) {
        if (e.by > py && (e.bx - e.ax) * (py - e.ay) - (px - e.ax) * (e.by - e.ay) > 0) {
          w++;
        }
      } else if (e.by <= py && (e.bx - e.ax) * (py - e.ay) - (px - e.ax) * (e.by - e.ay) < 0) {
        w--;
      }
    }
    return w;
  }
  /** The least winding any single ring has at x: a positively wound piece never goes below 0. */
  int least_ring(V2 x) const {
    int least = 0;
    for (auto & r : polys.rings) {
      least = std::min(least, ring_winding(r, x));
    }
    return least;
  }
};

} // namespace stroke_oracle

#endif /* GHOTI_IO_GSHGEOM_TESTS_STROKE_ORACLE_H */
