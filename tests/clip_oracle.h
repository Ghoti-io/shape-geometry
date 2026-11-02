/**
 * @file
 *
 * An independent oracle for the clip (clip.h, CK1 to CK6), by the definition and not by an algorithm.
 *
 * What the clip must do is a statement about winding numbers: at a point away from every edge, the result winds as the input
 * does inside the interior of the quadrilateral and not at all outside it.  This file holds nothing of the clip's method.
 * It has exact rational arithmetic on arbitrary-size integers (no `double`, no `__int128`: the library avoids both and so do
 * the oracles), the winding number of a set of closed polygons at a rational point by counting the crossings of a ray, the
 * distance from a point to a segment (to keep probes away from edges, CK2's exclusion), the exact rounding of the crossing of
 * an edge with a line (CK3) from the rational, and the area of the symmetric difference of two windings (CK4) by a sweep of
 * horizontal slabs: between two consecutive heights at which any edge ends or two edges meet, every edge crosses the slab
 * whole, the order in which they cross a horizontal line does not change, and the length of the stretches where the two
 * windings differ is linear in the height, so the area of a slab is that length at its middle times its height, exactly.
 *
 * Copyright 2026 by Corey Pennycuff
 */

#ifndef GHOTI_IO_GSHGEOM_TEST_CLIP_ORACLE_H
#define GHOTI_IO_GSHGEOM_TEST_CLIP_ORACLE_H

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <utility>
#include <vector>

namespace clip_oracle {

/* --- integers of any size ----------------------------------------------------------------------------------------- */

struct Big {
  std::vector<uint32_t> m;  // least significant first, no leading zero limb; zero is empty
  bool neg = false;

  Big() = default;
  explicit Big(int64_t v) {
    neg = v < 0;
    uint64_t a = v < 0 ? (uint64_t)0 - (uint64_t)v : (uint64_t)v;
    while (a) {
      m.push_back((uint32_t)(a & 0xFFFFFFFFu));
      a >>= 32;
    }
  }
  void trim() {
    while (!m.empty() && m.back() == 0) m.pop_back();
    if (m.empty()) neg = false;
  }
  int sign() const { return m.empty() ? 0 : (neg ? -1 : 1); }
};

inline int cmp_mag(const Big & a, const Big & b) {
  if (a.m.size() != b.m.size()) return a.m.size() < b.m.size() ? -1 : 1;
  for (size_t i = a.m.size(); i-- > 0;) {
    if (a.m[i] != b.m[i]) return a.m[i] < b.m[i] ? -1 : 1;
  }
  return 0;
}

inline std::vector<uint32_t> add_mag(const std::vector<uint32_t> & a, const std::vector<uint32_t> & b) {
  std::vector<uint32_t> r;
  uint64_t carry = 0;
  for (size_t i = 0; i < std::max(a.size(), b.size()) || carry; i++) {
    uint64_t s = carry + (i < a.size() ? a[i] : 0) + (i < b.size() ? b[i] : 0);
    r.push_back((uint32_t)(s & 0xFFFFFFFFu));
    carry = s >> 32;
  }
  return r;
}

/** a - b for |a| >= |b|. */
inline std::vector<uint32_t> sub_mag(const std::vector<uint32_t> & a, const std::vector<uint32_t> & b) {
  std::vector<uint32_t> r;
  int64_t borrow = 0;
  for (size_t i = 0; i < a.size(); i++) {
    int64_t d = (int64_t)a[i] - borrow - (i < b.size() ? (int64_t)b[i] : 0);
    borrow = d < 0;
    if (d < 0) d += (int64_t)1 << 32;
    r.push_back((uint32_t)d);
  }
  return r;
}

inline Big operator-(Big a) {
  if (!a.m.empty()) a.neg = !a.neg;
  return a;
}

inline Big operator+(const Big & a, const Big & b) {
  Big r;
  if (a.neg == b.neg) {
    r.m = add_mag(a.m, b.m);
    r.neg = a.neg;
  } else if (cmp_mag(a, b) >= 0) {
    r.m = sub_mag(a.m, b.m);
    r.neg = a.neg;
  } else {
    r.m = sub_mag(b.m, a.m);
    r.neg = b.neg;
  }
  r.trim();
  return r;
}

inline Big operator-(const Big & a, const Big & b) {
  return a + (-b);
}

inline Big operator*(const Big & a, const Big & b) {
  Big r;
  if (a.m.empty() || b.m.empty()) return r;
  r.m.assign(a.m.size() + b.m.size(), 0);
  for (size_t i = 0; i < a.m.size(); i++) {
    uint64_t carry = 0;
    for (size_t j = 0; j < b.m.size() || carry; j++) {
      uint64_t cur = r.m[i + j] + carry + (j < b.m.size() ? (uint64_t)a.m[i] * b.m[j] : 0);
      r.m[i + j] = (uint32_t)(cur & 0xFFFFFFFFu);
      carry = cur >> 32;
    }
  }
  r.neg = a.neg != b.neg;
  r.trim();
  return r;
}

inline int cmp(const Big & a, const Big & b) {
  const int d = (a - b).sign();
  return d;
}

/* --- rationals (not reduced: the sizes stay small enough for the cases) ------------------------------------------------ */

struct Rat {
  Big n;
  Big d{1};  // positive

  Rat() = default;
  Rat(int64_t v) : n(v) {}
  Rat(Big num, Big den) : n(std::move(num)), d(std::move(den)) {
    if (d.neg) {
      n = -n;
      d = -d;
    }
  }
  int sign() const { return n.sign(); }
};

inline Rat operator+(const Rat & a, const Rat & b) { return Rat(a.n * b.d + b.n * a.d, a.d * b.d); }
inline Rat operator-(const Rat & a, const Rat & b) { return Rat(a.n * b.d - b.n * a.d, a.d * b.d); }
inline Rat operator*(const Rat & a, const Rat & b) { return Rat(a.n * b.n, a.d * b.d); }
inline Rat operator/(const Rat & a, const Rat & b) { return Rat(a.n * b.d, a.d * b.n); }
inline int cmp(const Rat & a, const Rat & b) { return (a.n * b.d - b.n * a.d).sign(); }
inline bool operator<(const Rat & a, const Rat & b) { return cmp(a, b) < 0; }
inline bool operator==(const Rat & a, const Rat & b) { return cmp(a, b) == 0; }

/** floor(v + 1/2): the nearest integer, ties upward, by bisection on exact comparisons (no division). */
inline int64_t round_half_up(const Rat & v) {
  // The largest k with k - 1/2 <= v, i.e. (2k - 1) / 2 <= v.
  int64_t lo = -((int64_t)1 << 40);
  int64_t hi = (int64_t)1 << 40;
  while (lo < hi) {
    const int64_t mid = lo + (hi - lo + 1) / 2;
    if (cmp(Rat(Big(2 * mid - 1), Big(2)), v) <= 0) {
      lo = mid;
    } else {
      hi = mid - 1;
    }
  }
  return lo;
}

/* --- points and polygons ------------------------------------------------------------------------------------------------ */

struct IPt {
  int64_t x, y;
};
struct RPt {
  Rat x, y;
};
using Polygon = std::vector<IPt>;
using Polygons = std::vector<Polygon>;

inline bool operator==(const IPt & a, const IPt & b) { return a.x == b.x && a.y == b.y; }

/** The exact side value (l1 - l0) x (p - l0) at a rational point. */
inline Rat side_at(const IPt & l0, const IPt & l1, const RPt & p) {
  return Rat(l1.x - l0.x) * (p.y - Rat(l0.y)) - Rat(l1.y - l0.y) * (p.x - Rat(l0.x));
}

/** The exact crossing of the segment a-b with the line l0-l1 rounded half up, or false if they are not strictly on opposite sides. */
inline bool crossing_rounded(const IPt & a, const IPt & b, const IPt & l0, const IPt & l1, IPt * out) {
  const Rat fa = side_at(l0, l1, RPt{Rat(a.x), Rat(a.y)});
  const Rat fb = side_at(l0, l1, RPt{Rat(b.x), Rat(b.y)});
  if (fa.sign() * fb.sign() >= 0) return false;
  const Rat t = fa / (fa - fb);
  out->x = round_half_up(Rat(a.x) + Rat(b.x - a.x) * t);
  out->y = round_half_up(Rat(a.y) + Rat(b.y - a.y) * t);
  return true;
}

/**
 * The winding number of the closed polygons about a rational point not on any edge: the signed number of times the edges cross
 * the horizontal ray to the right, a downward edge counting +1.  (A point on an edge has no winding number; the callers keep
 * their probes away from the edges.)
 */
inline int winding(const Polygons & polys, const RPt & q) {
  int w = 0;
  for (const Polygon & poly : polys) {
    const size_t n = poly.size();
    for (size_t k = 0; k < n; k++) {
      const IPt a = poly[k];
      const IPt b = poly[(k + 1) % n];
      const bool a_above = Rat(a.y) < q.y || Rat(a.y) == q.y;  // a.y <= q.y
      const bool b_above = Rat(b.y) < q.y || Rat(b.y) == q.y;
      if (a_above == b_above) continue;  // does not cross the line y = q.y (half-open)
      // x of the edge at y = q.y, compared with q.x.
      const Rat x = Rat(a.x) + Rat(b.x - a.x) * (q.y - Rat(a.y)) / Rat(b.y - a.y);
      if (q.x < x) w += b.y > a.y ? 1 : -1;
    }
  }
  return w;
}

/** The squared distance from the point to the segment a-b, exactly. */
inline Rat dist2_to_segment(const IPt & a, const IPt & b, const RPt & q) {
  const Rat ex(b.x - a.x), ey(b.y - a.y);
  const Rat fx = q.x - Rat(a.x), fy = q.y - Rat(a.y);
  const Rat len2 = ex * ex + ey * ey;
  if (len2.sign() == 0) return fx * fx + fy * fy;
  const Rat dot = fx * ex + fy * ey;
  if (dot.sign() <= 0) return fx * fx + fy * fy;
  if (cmp(dot, len2) >= 0) {
    const Rat gx = q.x - Rat(b.x), gy = q.y - Rat(b.y);
    return gx * gx + gy * gy;
  }
  const Rat cross = ex * fy - ey * fx;
  return cross * cross / len2;
}

/** Whether the point is within `r` (an integer number of units) of an edge of any polygon (closed). */
inline bool near_edges(const Polygons & polys, const RPt & q, int64_t r) {
  const Rat r2(r * r);
  for (const Polygon & poly : polys) {
    for (size_t k = 0; k < poly.size(); k++) {
      if (cmp(dist2_to_segment(poly[k], poly[(k + 1) % poly.size()], q), r2) <= 0) return true;
    }
  }
  return false;
}

/** The orientation of a convex quadrilateral: +1 or -1 by the sign of its turns, 0 if it has none. */
inline int quad_sigma(const Polygon & quad) {
  int sigma = 0;
  for (size_t i = 0; i < 4; i++) {
    const IPt a = quad[(i + 3) % 4], b = quad[i], c = quad[(i + 1) % 4];
    const int t = (Rat((b.x - a.x)) * Rat(c.y - a.y) - Rat(b.y - a.y) * Rat(c.x - a.x)).sign();
    if (t) sigma = t;
  }
  return sigma;
}

/** Whether the point is strictly inside the convex quadrilateral (a side of zero length does not count). */
inline bool strictly_inside(const Polygon & quad, const RPt & q) {
  const int sigma = quad_sigma(quad);
  for (size_t i = 0; i < 4; i++) {
    const IPt a = quad[i], b = quad[(i + 1) % 4];
    if (a == b) continue;
    if (side_at(a, b, q).sign() * sigma <= 0) return false;
  }
  return true;
}

/** The result of the area measure: exact area of the symmetric difference, in square units, as a Rat. */
struct AreaResult {
  Rat area;
  size_t slabs = 0;
};

namespace detail {

struct Seg {
  IPt a, b;
  int kind;  // 0 input, 1 output, 2 the quadrilateral
};

/** The height at which the two segments meet, if they properly cross or touch in a single point. */
inline bool meet_height(const Seg & s, const Seg & t, Rat * y) {
  const Rat d1x(s.b.x - s.a.x), d1y(s.b.y - s.a.y);
  const Rat d2x(t.b.x - t.a.x), d2y(t.b.y - t.a.y);
  const Rat den = d1x * d2y - d1y * d2x;
  if (den.sign() == 0) return false;
  const Rat wx(t.a.x - s.a.x), wy(t.a.y - s.a.y);
  const Rat u = (wx * d2y - wy * d2x) / den;  // along s
  const Rat v = (wx * d1y - wy * d1x) / den;  // along t
  const Rat zero(0), one(1);
  if (cmp(u, zero) < 0 || cmp(u, one) > 0 || cmp(v, zero) < 0 || cmp(v, one) > 0) return false;
  *y = Rat(s.a.y) + d1y * u;
  return true;
}

}  // namespace detail

/**
 * The area of the symmetric difference, weighted by the difference of the windings, between `out` and the exact clip of `in`
 * (which winds as `in` does strictly inside `quad` and not at all outside).
 */
inline AreaResult symmetric_difference_area(const Polygons & in, const Polygons & out, const Polygon & quad) {
  using detail::Seg;
  std::vector<Seg> segs;
  auto take = [&](const Polygons & ps, int kind) {
    for (const Polygon & p : ps) {
      for (size_t k = 0; k < p.size(); k++) {
        const IPt a = p[k], b = p[(k + 1) % p.size()];
        if (a.y != b.y) segs.push_back(Seg{a, b, kind});  // a horizontal edge is crossed by no horizontal line
      }
    }
  };
  take(in, 0);
  take(out, 1);
  {
    Polygons q{quad};
    take(q, 2);
  }
  std::vector<Rat> ys;
  for (const Seg & s : segs) {
    ys.push_back(Rat(s.a.y));
    ys.push_back(Rat(s.b.y));
  }
  // Horizontal edges end at heights of their own; their end points are vertices of other edges too, or the polygon degenerates.
  for (size_t i = 0; i < segs.size(); i++) {
    for (size_t j = i + 1; j < segs.size(); j++) {
      Rat y;
      if (detail::meet_height(segs[i], segs[j], &y)) ys.push_back(y);
    }
  }
  std::sort(ys.begin(), ys.end(), [](const Rat & a, const Rat & b) { return a < b; });
  ys.erase(std::unique(ys.begin(), ys.end(), [](const Rat & a, const Rat & b) { return a == b; }), ys.end());
  AreaResult res;
  res.area = Rat(0);
  const int sigma = quad_sigma(quad);
  (void)sigma;
  for (size_t s = 0; s + 1 < ys.size(); s++) {
    const Rat ym = (ys[s] + ys[s + 1]) / Rat(2);
    struct Hit {
      Rat x;
      int kind;
      int dir;
    };
    std::vector<Hit> hits;
    for (const Seg & g : segs) {
      const int64_t lo = std::min(g.a.y, g.b.y), hi = std::max(g.a.y, g.b.y);
      if (!(Rat(lo) < ym && ym < Rat(hi))) continue;
      const Rat x = Rat(g.a.x) + Rat(g.b.x - g.a.x) * (ym - Rat(g.a.y)) / Rat(g.b.y - g.a.y);
      hits.push_back(Hit{x, g.kind, g.b.y > g.a.y ? 1 : -1});
    }
    std::sort(hits.begin(), hits.end(), [](const Hit & a, const Hit & b) { return a.x < b.x; });
    int w_in = 0, w_out = 0, q_in = 0;
    Rat length(0);
    for (size_t k = 0; k + 1 < hits.size(); k++) {
      const Hit & h = hits[k];
      if (h.kind == 0) w_in += h.dir;
      if (h.kind == 1) w_out += h.dir;
      if (h.kind == 2) q_in += h.dir;
      // The interval between hit k and hit k+1: inside the convex quadrilateral iff its winding is nonzero.
      const int w_exact = q_in != 0 ? w_in : 0;
      const int diff = std::abs(w_out - w_exact);
      if (diff) length = length + Rat(diff) * (hits[k + 1].x - h.x);
    }
    res.area = res.area + length * (ys[s + 1] - ys[s]);
    res.slabs++;
  }
  return res;
}

}  // namespace clip_oracle

#endif /* GHOTI_IO_GSHGEOM_TEST_CLIP_ORACLE_H */
