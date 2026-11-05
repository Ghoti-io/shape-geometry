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
 * An independent oracle for the pieces of a hatch (hatch.h, HX4 to HX6), by the definition and not by the library's method,
 * in exact rationals.  It shares no code with the library (its arithmetic is clip_oracle.h's rationals) and none with the
 * generator of the committed vectors (tools/gen-hatch-vectors.py), which states the same definition in Python.
 *
 * The trim set of a line is found as the definition says: every position where a contour edge meets the line is a
 * breakpoint (the projection of both ends for an edge that lies along it); the open stretch between two consecutive
 * breakpoints is tested at two exact points, 2^-100 of the line's direction to the left of the midpoint and to the right
 * of it (a distance smaller than any feature of the inputs the tests use), with the winding number and the even-odd
 * parity of every contour taken by a ray in exact rationals.  The dash positions are those of the DXF elements laid on the
 * infinite line from the base point (a run of pen-down elements is one dash, a dot first with a dash last is kept apart).
 */

#ifndef GHOTI_IO_GSHGEOM_TEST_HATCH_ORACLE_H
#define GHOTI_IO_GSHGEOM_TEST_HATCH_ORACLE_H

#include "lattice_oracle.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace hatch_oracle {

using clip_oracle::Big;
using clip_oracle::Rat;
using lattice_oracle::from_double;

inline Rat pow2(int e) {
  Big p(1);
  Big two(2);
  for (int i = 0; i < (e < 0 ? -e : e); i++) p = p * two;
  return e >= 0 ? Rat(p, Big(1)) : Rat(Big(1), p);
}

/** The nearest long double of a rational (for comparing against a double within 2^-40). */
inline long double to_ld(const Big & b) {
  long double r = 0;
  for (size_t i = b.m.size(); i-- > 0;) r = r * 4294967296.0L + (long double)b.m[i];
  return b.neg ? -r : r;
}

inline double to_double(const Rat & r) {
  if (r.n.m.empty()) return 0.0;
  // scale both by a power of two so that neither overflows long double
  size_t ln = r.n.m.size(), ld = r.d.m.size();
  const size_t drop = std::min(ln, ld) > 400 ? std::min(ln, ld) - 400 : 0;
  Big n = r.n, d = r.d;
  n.m.erase(n.m.begin(), n.m.begin() + (std::ptrdiff_t)std::min(drop, n.m.size()));
  d.m.erase(d.m.begin(), d.m.begin() + (std::ptrdiff_t)std::min(drop, d.m.size()));
  if (d.m.empty()) return 0.0;
  return (double)(to_ld(n) / to_ld(d));
}

struct Pt {
  int64_t x, y;  // 26.6
};
using Contour = std::vector<Pt>;
using Contours = std::vector<Contour>;

struct RP {
  Rat x, y;
};

inline Rat px(int64_t v) { return Rat(Big(v), Big(64)); }

/** Membership of the exact point @p p in the region of the nesting style (HX4), by a ray in +x. */
inline bool in_region(const Contours & cs, int nesting, int rule, const RP & p) {
  int64_t wind = 0, odd = 0;
  for (const Contour & c : cs) {
    int64_t w = 0, cr = 0;
    for (size_t i = 0; i < c.size(); i++) {
      const Pt a = c[i], b = c[(i + 1) % c.size()];
      const Rat ay = px(a.y), by = px(b.y);
      if ((clip_oracle::cmp(ay, p.y) <= 0) != (clip_oracle::cmp(by, p.y) <= 0)) {
        const Rat t = (p.y - ay) / (by - ay);
        const Rat x = px(a.x) + t * (px(b.x) - px(a.x));
        if (clip_oracle::cmp(x, p.x) > 0) {
          cr++;
          w += b.y > a.y ? 1 : -1;
        }
      }
    }
    wind += w;
    odd += cr & 1;
  }
  if (nesting == 1) return odd == 1;
  if (nesting == 2) return odd >= 1;
  return rule == 0 ? wind != 0 : (wind % 2) != 0;
}

struct Iv {
  Rat a, b;
};

/** The trim intervals of the line `A + s w`. */
inline std::vector<Iv> trim_intervals(const Contours & cs, int nesting, int rule, const RP & A, const RP & w) {
  const Rat ww = w.x * w.x + w.y * w.y;
  auto proj = [&](const RP & p) { return ((p.x - A.x) * w.x + (p.y - A.y) * w.y) / ww; };
  auto side = [&](const RP & p) { return w.x * (p.y - A.y) - w.y * (p.x - A.x); };
  std::vector<Rat> brk;
  for (const Contour & c : cs) {
    for (size_t i = 0; i < c.size(); i++) {
      const RP a{px(c[i].x), px(c[i].y)};
      const RP b{px(c[(i + 1) % c.size()].x), px(c[(i + 1) % c.size()].y)};
      const Rat sa = side(a), sb = side(b);
      const int za = sa.sign(), zb = sb.sign();
      if (za == 0 && zb == 0) {
        brk.push_back(proj(a));
        brk.push_back(proj(b));
      } else if (za == 0) {
        brk.push_back(proj(a));
      } else if (zb == 0) {
        brk.push_back(proj(b));
      } else if ((za > 0) != (zb > 0)) {
        const Rat t = sa / (sa - sb);
        brk.push_back(proj(RP{a.x + t * (b.x - a.x), a.y + t * (b.y - a.y)}));
      }
    }
  }
  std::sort(brk.begin(), brk.end(), [](const Rat & p, const Rat & q) { return clip_oracle::cmp(p, q) < 0; });
  std::vector<Rat> bs;
  for (const Rat & r : brk) {
    if (bs.empty() || clip_oracle::cmp(bs.back(), r) != 0) bs.push_back(r);
  }
  const Rat eps = pow2(-100);
  std::vector<Iv> out;
  for (size_t i = 0; i + 1 < bs.size(); i++) {
    const Rat m = (bs[i] + bs[i + 1]) / Rat(2);
    const RP M{A.x + m * w.x, A.y + m * w.y};
    bool ok = false;
    for (int sg = -1; sg <= 1; sg += 2) {
      // eps along the left normal (-w.y, w.x), in units of the line's own length scale (exact, never zero)
      const Rat e = eps * Rat(sg);
      const RP q{M.x - e * w.y, M.y + e * w.x};
      if (in_region(cs, nesting, rule, q)) ok = true;
    }
    if (!ok) continue;
    if (!out.empty() && clip_oracle::cmp(out.back().b, bs[i]) == 0) {
      out.back().b = bs[i + 1];
    } else {
      out.push_back(Iv{bs[i], bs[i + 1]});
    }
  }
  return out;
}

/** The pen-down pieces (a, b; a == b a dot) of the pattern of DXF elements on the infinite line, covering [lo, hi] and a period beyond. */
inline std::vector<Iv> dash_set(const std::vector<double> & elems, const Rat & lo, const Rat & hi, bool * solid) {
  struct Run {
    char kind;
    Rat st, ln;
  };
  std::vector<Run> runs;
  Rat pos(0);
  for (double e : elems) {
    const char kind = e >= 0 ? 'd' : 'g';
    const Rat ln = from_double(e < 0 ? -e : e);
    if (!runs.empty() && runs.back().kind == kind) {
      runs.back().ln = runs.back().ln + ln;
    } else {
      runs.push_back(Run{kind, pos, ln});
    }
    pos = pos + ln;
  }
  const Rat period = pos;
  *solid = true;
  for (const Run & r : runs) {
    if (r.kind == 'g') *solid = false;
  }
  std::vector<Iv> out;
  if (*solid) {
    out.push_back(Iv{lo, hi});
    return out;
  }
  const bool merge = runs.front().kind == 'd' && runs.back().kind == 'd' && runs.size() > 1 && runs.front().ln.sign() != 0;
  std::vector<Iv> items;
  for (size_t i = 0; i < runs.size(); i++) {
    if (runs[i].kind != 'd') continue;
    if (merge && i == 0) continue;
    if (merge && i + 1 == runs.size()) {
      items.push_back(Iv{runs[i].st, period + runs.front().ln});
    } else {
      items.push_back(Iv{runs[i].st, runs[i].st + runs[i].ln});
    }
  }
  // the periods that can reach [lo, hi]: q from floor(lo / period) - 2
  const double pd = to_double(period);
  const int64_t q0 = (int64_t)std::floor(to_double(lo) / pd) - 2;
  const int64_t q1 = (int64_t)std::floor(to_double(hi) / pd) + 2;
  for (int64_t q = q0; q <= q1; q++) {
    for (const Iv & it : items) {
      out.push_back(Iv{it.a + period * Rat(Big(q), Big(1)), it.b + period * Rat(Big(q), Big(1))});
    }
  }
  std::sort(out.begin(), out.end(), [](const Iv & p, const Iv & q) { return clip_oracle::cmp(p.a, q.a) < 0; });
  return out;
}

struct Piece {
  bool dot;
  Rat s0, s1;  // positions on the line
  RP a;        // the hatch-space base point O_k, as exact rationals of the doubles the library forms
  double ox, oy;  // the same, as the doubles
};

struct Line {
  uint64_t k;
  std::vector<Piece> pieces;
};

}  // namespace hatch_oracle

#endif
