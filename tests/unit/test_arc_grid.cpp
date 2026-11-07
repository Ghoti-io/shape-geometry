/**
 * @file
 *
 * The arc constructors, readers and decompositions over a grid of special and ordinary arguments.
 *
 * SPDX-License-Identifier: LGPL-3.0-only
 *
 * Copyright 2026 by Corey Pennycuff
 */

#include "../test_helpers.h"

extern "C" {
#include "../../src/exact/exact.h"
}

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

/**
 * Every entry point of arc.h is called on every combination of the values that decide a check (zero, one, a negative, a huge, the
 * coordinate limit and the first number past it, NaN, an infinity, a quarter turn and the full turn in each unit, a unit past the
 * last, a flag past 1) and what each returns (the result, and everything it wrote: the path's hash, the arc, the pieces) is folded
 * into one hash, which is pinned.  A check removed or weakened, a comparison changed, a constant moved or an error ignored changes
 * the answer for some combination, and the hash moves.  The values are the contract's: tests/unit/test_arc_rules.cpp holds the
 * assertion of each rule alone; this holds them all at once, and over the combinations in which one rule's refusal hides behind
 * another's.
 *
 * The hash is of the library as it stands (reviewed against arc.h for each refusal's code); a change that moves it changes what a
 * caller sees, and says so by failing here: regenerate it (GSHGEOM_ARC_GRID_PRINT=1) only for a change that is meant to.
 */

namespace {

constexpr uint64_t kGridHash = 0xf64c369796dded10ull;

struct Fold {
  uint64_t h = 1469598103934665603ull;
  void word(uint64_t w) {
    for (int i = 0; i < 8; i++) {
      h ^= (w >> (8 * i)) & 0xffu;
      h *= 1099511628211ull;
    }
  }
  void dbl(double d) {
    uint64_t w;
    std::memcpy(&w, &d, 8);
    // all NaNs are one for the hash: the library does not promise a payload
    word(std::isnan(d) ? 0x7ff8000000000000ull : w);
  }
  void pt(const GSHGEOM_Point & p) {
    dbl(p.x);
    dbl(p.y);
    dbl(p.z);
  }
  void arc(const GSHGEOM_Arc & a) {
    pt(a.start);
    pt(a.end);
    pt(a.center);
    pt(a.u);
    pt(a.v);
    unsigned unit_bits;  // the field's bytes: a record a failed call left alone may hold a value that is no unit
    std::memcpy(&unit_bits, &a.unit, sizeof unit_bits);
    word(unit_bits);
    dbl(a.t0);
    dbl(a.sweep);
    word(a.flags);
  }
};

GSHGEOM_Point P(double x, double y) {
  return GSHGEOM_Point{x, y, 0.0};
}

GSHGEOM_Path * make(const GSHGEOM_Limits * lim = nullptr) {
  GSHGEOM_Path * p = nullptr;
  gshgeom_path_create(nullptr, lim, 2, &p);
  return p;
}

const double kBig = 0x1p100;
const double kPastBig = 0x1.0000000000001p100;
const double kTwoPiBelow = 0x1.921fb54442d17p+2;
const double kTwoPi = 6.283185307179586;

/** One call's record: the result, then the path (its hash and element count) so a write after a failure shows. */
void fold_path(Fold & f, GSHGEOM_Result r, GSHGEOM_Path * p) {
  f.word((uint64_t)r);
  f.word(gshgeom_path_hash(p));
  f.word(gshgeom_path_element_count(p));
}

void endpoint_grid(Fold & f) {
  const double ends[] = {0.0, 1.0, 7.5, NAN, INFINITY, kBig, kPastBig};
  const double radii[] = {0.0, -0.0, 1.0, 5.0, -1.0, NAN, INFINITY, 0x1p-1000, kPastBig};
  const double phis[] = {0.0, 90.0, NAN, 1e300};
  const int units[] = {GSHGEOM_ANGLE_RADIANS, GSHGEOM_ANGLE_DEGREES, GSHGEOM_ANGLE_GRADIANS, GSHGEOM_ANGLE_TURNS,
      (int)GSHGEOM_ANGLE_UNIT_COUNT, -1};
  const int flags[] = {0, 1, 2, -1};
  for (int current = 0; current < 2; current++) {
    for (double ex : ends) {
      for (double ey : ends) {
        for (double rx : radii) {
          for (double ry : radii) {
            for (double phi : phis) {
              for (int u : units) {
                for (int large : flags) {
                  for (int sweep : flags) {
                    GSHGEOM_Path * p = make();
                    if (current) {
                      gshgeom_path_move_to(p, P(0.0, 0.0));
                    }
                    GSHGEOM_Arc_Kind kind = GSHGEOM_ARC_KIND_ARC;  // a value the call can leave alone (it writes the kind on success only)
                    GSHGEOM_Result r = gshgeom_path_arc_endpoint_to(p, P(ex, ey), rx, ry, phi, (GSHGEOM_Angle_Unit)u, large, sweep, &kind);
                    fold_path(f, r, p);
                    f.word((uint64_t)kind);
                    gshgeom_path_destroy(p);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

void center_grid(Fold & f) {
  const double cs[] = {0.0, NAN, kBig};
  const double radii[] = {0.0, 1.0, 5.0, -1.0, NAN, INFINITY, kPastBig};
  const double phis[] = {0.0, 30.0, NAN};
  const int units[] = {GSHGEOM_ANGLE_RADIANS, GSHGEOM_ANGLE_DEGREES, GSHGEOM_ANGLE_GRADIANS, GSHGEOM_ANGLE_TURNS,
      (int)GSHGEOM_ANGLE_UNIT_COUNT};
  const double t0s[] = {0.0, 90.0, NAN, 720.0, 1e300, -0.0, kTwoPiBelow};
  const double sweeps[] = {0.0, 90.0, 360.0, 361.0, -360.0, NAN, -0.0, kTwoPi, kTwoPiBelow, 1.0, 400.0};
  const GSHGEOM_Point ends[] = {P(0.0, 0.0), P(1000.0, 0.0), P(NAN, 0.0), P(0.0, INFINITY)};
  for (int current = 0; current < 2; current++) {
    for (double cx : cs) {
      for (double cy : cs) {
        for (double rx : radii) {
          for (double ry : radii) {
            for (double phi : phis) {
              for (int u : units) {
                for (double t0 : t0s) {
                  for (double sw : sweeps) {
                    for (int e = -1; e < 4; e++) {
                      GSHGEOM_Path * p = make();
                      if (current) {
                        gshgeom_path_move_to(p, P(1.0, 1.0));
                      }
                      GSHGEOM_Result r = gshgeom_path_arc_center_to(p, P(cx, cy), rx, ry, phi, (GSHGEOM_Angle_Unit)u, t0, sw,
                          e < 0 ? nullptr : &ends[e]);
                      fold_path(f, r, p);
                      gshgeom_path_destroy(p);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

void tangent_bulge_grid(Fold & f) {
  const double v[] = {0.0, 1.0, 10.0, NAN, INFINITY, -10.0, kBig};
  const double rs[] = {0.0, 1.0, 3.0, NAN, -1.0, INFINITY, 0x1p-1000, kBig};
  for (int current = 0; current < 2; current++) {
    for (double a : v) {
      for (double b : v) {
        for (double c : v) {
          for (double d : v) {
            for (double r : rs) {
              GSHGEOM_Path * p = make();
              if (current) {
                gshgeom_path_move_to(p, P(0.0, 0.0));
              }
              GSHGEOM_Result res = gshgeom_path_arc_tangent_to(p, P(a, b), P(c, d), r);
              fold_path(f, res, p);
              gshgeom_path_destroy(p);
            }
          }
        }
      }
    }
    for (double qx : v) {
      for (double qy : v) {
        for (double b : std::vector<double>{0.0, 1.0, -1.0, 0.5, (double)NAN, (double)INFINITY, 1e300, 1e-300, -1e-300, 4.0}) {
          GSHGEOM_Path * p = make();
          if (current) {
            gshgeom_path_move_to(p, P(0.0, 0.0));
          }
          GSHGEOM_Result res = gshgeom_path_arc_bulge_to(p, P(qx, qy), b);
          fold_path(f, res, p);
          gshgeom_path_destroy(p);
        }
      }
    }
  }
}

/** A deterministic family of arc records: consistent ones, and each field of one disturbed. */
std::vector<GSHGEOM_Arc> arcs() {
  std::vector<GSHGEOM_Arc> out;
  const GSHGEOM_Angle_Unit units[] = {GSHGEOM_ANGLE_RADIANS, GSHGEOM_ANGLE_DEGREES, GSHGEOM_ANGLE_GRADIANS, GSHGEOM_ANGLE_TURNS};
  const double turns[] = {kTwoPi, 360.0, 400.0, 1.0};
  const double quarters[] = {kTwoPi / 4.0, 90.0, 100.0, 0.25};
  uint64_t s = 88172645463325252ull;
  auto rnd = [&]() {
    s ^= s << 13;
    s ^= s >> 7;
    s ^= s << 17;
    return (double)(s >> 11) / 9007199254740992.0;
  };
  for (int i = 0; i < 4; i++) {
    for (int k = 0; k < 90; k++) {
      GSHGEOM_Arc a;
      std::memset(&a, 0, sizeof a);
      a.unit = units[i];
      const double mult[] = {1e-9, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, -0.25, -1.0, -2.0, -3.5, 4.0 / 3.0};
      const double t0s[] = {0.0, 1.0, -0.5, 4.5, 2.0};
      double t0 = t0s[k % 5] * quarters[i];
      double sw = mult[(k / 5) % 13] * quarters[i];
      if (k >= 65) {  // random bases
        t0 = (rnd() - 0.5) * 3.0 * turns[i];
        sw = (rnd() - 0.5) * 2.0 * turns[i];
      }
      a.center = P(rnd() * 20.0 - 10.0, rnd() * 20.0 - 10.0);
      a.u = P(rnd() * 8.0 + 1.0, rnd() * 4.0);
      const double ratio = 0.3 + rnd();
      a.v = P(-a.u.y * ratio, a.u.x * ratio);
      if (k % 7 == 0) {  // a circle
        a.v = P(-a.u.y, a.u.x);
        a.flags = GSHGEOM_ARC_CIRCULAR;
      }
      a.t0 = t0;
      a.sweep = sw;
      GSHGEOM_Arc_Form fm = {a.center, a.u, a.v, a.unit, a.t0, a.sweep, a.flags};
      gshgeom_arc_point(&fm, a.t0, &a.start);
      gshgeom_arc_point(&fm, a.t0 + a.sweep, &a.end);
      out.push_back(a);
      if (k < 40) {  // each field disturbed in turn, from this base
        const double bad[] = {NAN, INFINITY, -INFINITY, 0.0, kBig, kPastBig, -kPastBig, 1e-300};
        for (double b : bad) {
          for (int field = 0; field < 16; field++) {
            GSHGEOM_Arc d = a;
            double * targets[] = {&d.start.x, &d.start.y, &d.end.x, &d.end.y, &d.center.x, &d.center.y, &d.u.x, &d.u.y, &d.v.x,
                &d.v.y, &d.t0, &d.sweep, &d.start.z, &d.end.z, &d.center.z, &d.u.z};
            *targets[field] = b;
            out.push_back(d);
          }
        }
        GSHGEOM_Arc d = a;
        d.unit = GSHGEOM_ANGLE_UNIT_COUNT;
        out.push_back(d);
        d = a;
        d.flags = 2u;
        out.push_back(d);
        d = a;
        d.flags = a.flags ^ 1u;
        out.push_back(d);
        d = a;
        d.sweep = a.sweep >= 0 ? turns[i] : -turns[i];
        out.push_back(d);
        d = a;
        d.sweep = (a.sweep >= 0 ? 1.0 : -1.0) * std::nextafter(turns[i], 2.0 * turns[i]);
        out.push_back(d);
        d = a;
        d.v = d.u;
        out.push_back(d);
        d = a;
        d.v = P(-a.u.y * 0x1p-1000, a.u.x * 0x1p-1000);
        out.push_back(d);
        d = a;
        d.end.x = std::nextafter(a.end.x, 1e300);
        out.push_back(d);
        d = a;
        d.end.x = a.end.x + 1e-9 * (std::fabs(a.u.x) + 1.0);
        out.push_back(d);
        d = a;
        d.start.y = a.start.y + 1e-9 * (std::fabs(a.u.x) + 1.0);
        out.push_back(d);
      }
    }
  }
  return out;
}

void record_grid(Fold & f) {
  const std::vector<GSHGEOM_Arc> all = arcs();
  GSHGEOM_Map ident;
  gshgeom_map_identity(2, &ident);
  GSHGEOM_Map maps[6];
  gshgeom_map_scale(2, 2.0, 2.0, 1.0, &maps[0]);
  gshgeom_map_scale(2, 2.0, 0.5, 1.0, &maps[1]);
  gshgeom_map_scale(2, 0.0, 1.0, 1.0, &maps[2]);
  gshgeom_map_rotate(GSHGEOM_ANGLE_DEGREES, 30.0, &maps[3]);
  gshgeom_map_scale(2, -1.0, 1.0, 1.0, &maps[4]);
  gshgeom_map_translate(2, 1e30, -3.0, 0.0, &maps[5]);
  GSHGEOM_Map skew = maps[3];
  skew.m[0][1] += 0.5;
  for (const GSHGEOM_Arc & a : all) {
    // append at a path whose current point is the arc's start, at another, at none
    for (int mode = 0; mode < 3; mode++) {
      GSHGEOM_Path * p = make();
      if (mode == 0) {
        gshgeom_path_move_to(p, a.start);
      } else if (mode == 1) {
        gshgeom_path_move_to(p, P(a.start.x + 1.0, a.start.y));
      }
      GSHGEOM_Result r = gshgeom_path_append_arc(p, &a);
      fold_path(f, r, p);
      if (r == GSHGEOM_OK) {
        GSHGEOM_Arc back;
        std::memset(&back, 0x5a, sizeof back);
        GSHGEOM_Result rr = gshgeom_path_arc(p, 1, &back);
        f.word((uint64_t)rr);
        f.arc(back);
        f.word((uint64_t)gshgeom_path_arc(p, 0, &back));
        f.word((uint64_t)gshgeom_path_arc(p, 2, &back));
        GSHGEOM_Path * low = nullptr;
        GSHGEOM_Result lr = gshgeom_path_lower_arcs(p, nullptr, nullptr, &low);
        f.word((uint64_t)lr);
        if (lr == GSHGEOM_OK && low) {
          f.word(gshgeom_path_hash(low));
          f.word(gshgeom_path_element_count(low));
        }
        gshgeom_path_destroy(low);
      }
      gshgeom_path_destroy(p);
    }
    // the pieces, the circle, the point
    GSHGEOM_Curve pieces[GSHGEOM_ARC_MAX_PIECES];
    std::memset(pieces, 0, sizeof pieces);
    int n = -5;
    GSHGEOM_Result pr = gshgeom_arc_pieces(&a, pieces, &n);
    f.word((uint64_t)pr);
    f.word((uint64_t)(int64_t)n);
    if (pr == GSHGEOM_OK) {
      for (int k = 0; k < n; k++) {
        f.word((uint64_t)pieces[k].kind);
        f.word((uint64_t)pieces[k].dim);
        for (int q = 0; q < 4; q++) {
          f.pt(pieces[k].p[q]);
          f.dbl(pieces[k].w[q]);
        }
      }
    }
    GSHGEOM_Point c = P(-1.0, -1.0);
    double radius = -1.0;
    int ccw = -7;
    GSHGEOM_Result cr = gshgeom_arc_circle(&a, &c, &radius, &ccw);
    f.word((uint64_t)cr);
    f.pt(c);
    f.dbl(radius);
    f.word((uint64_t)(int64_t)ccw);
    GSHGEOM_Arc_Form fm = {a.center, a.u, a.v, a.unit, a.t0, a.sweep, a.flags};
    for (double t : std::vector<double>{a.t0, a.t0 + a.sweep, 0.0, 1e300, (double)NAN}) {
      GSHGEOM_Point pt = P(-9.0, -9.0);
      GSHGEOM_Result tr = gshgeom_arc_point(&fm, t, &pt);
      f.word((uint64_t)tr);
      f.pt(pt);
    }
    // under maps
    for (const GSHGEOM_Map * m : std::vector<const GSHGEOM_Map *>{&ident, &maps[0], &maps[1], &maps[2], &maps[3], &maps[4], &maps[5], &skew}) {
      GSHGEOM_Arc out;
      std::memset(&out, 0x3c, sizeof out);
      GSHGEOM_Result mr = gshgeom_arc_map(m, &a, &out);
      f.word((uint64_t)mr);
      f.arc(out);
    }
  }
}


/** The entry points on what they must refuse whatever else is right: no path, no record, no output, an index past the end, a path of three dimensions. */
void refusal_grid(Fold & f) {
  GSHGEOM_Arc a;
  std::memset(&a, 0, sizeof a);
  a.center = P(2.0, 3.0);
  a.u = P(4.0, 0.0);
  a.v = P(0.0, 4.0);
  a.unit = GSHGEOM_ANGLE_DEGREES;
  a.t0 = 10.0;
  a.sweep = 100.0;
  a.flags = GSHGEOM_ARC_CIRCULAR;
  GSHGEOM_Arc_Form fm = {a.center, a.u, a.v, a.unit, a.t0, a.sweep, a.flags};
  gshgeom_arc_point(&fm, a.t0, &a.start);
  gshgeom_arc_point(&fm, a.t0 + a.sweep, &a.end);
  GSHGEOM_Curve pieces[GSHGEOM_ARC_MAX_PIECES];
  int n = 0;
  GSHGEOM_Point c;
  double radius;
  int ccw;
  GSHGEOM_Arc back;
  GSHGEOM_Map ident;
  gshgeom_map_identity(2, &ident);
  GSHGEOM_Point endpt = P(9.0, 9.0);
  for (int dim = 2; dim <= 3; dim++) {
    for (int withcur = 0; withcur < 2; withcur++) {
      GSHGEOM_Path * p = nullptr;
      gshgeom_path_create(nullptr, nullptr, dim, &p);
      if (withcur) {
        gshgeom_path_move_to(p, a.start);
      }
      f.word((uint64_t)gshgeom_path_append_arc(p, &a));
      f.word((uint64_t)gshgeom_path_append_arc(p, nullptr));
      f.word((uint64_t)gshgeom_path_append_arc(nullptr, &a));
      for (size_t idx : {(size_t)0, (size_t)1, (size_t)2, (size_t)3, (size_t)100, (size_t)-1}) {
        std::memset(&back, 0x11, sizeof back);
        f.word((uint64_t)gshgeom_path_arc(p, idx, &back));
        f.arc(back);
        f.word((uint64_t)gshgeom_path_arc(p, idx, nullptr));
      }
      f.word((uint64_t)gshgeom_path_arc(nullptr, 0, &back));
      GSHGEOM_Arc_Kind kind;
      f.word((uint64_t)gshgeom_path_arc_endpoint_to(p, P(9.0, 1.0), 5.0, 5.0, 0.0, GSHGEOM_ANGLE_DEGREES, 0, 1, &kind));
      f.word((uint64_t)gshgeom_path_arc_endpoint_to(p, P(9.0, 1.0), 5.0, 5.0, 0.0, GSHGEOM_ANGLE_DEGREES, 0, 1, nullptr));
      f.word((uint64_t)gshgeom_path_arc_endpoint_to(nullptr, P(9.0, 1.0), 5.0, 5.0, 0.0, GSHGEOM_ANGLE_DEGREES, 0, 1, &kind));
      f.word((uint64_t)gshgeom_path_arc_conjugate_to(p, &fm, &endpt));
      f.word((uint64_t)gshgeom_path_arc_conjugate_to(p, &fm, nullptr));
      f.word((uint64_t)gshgeom_path_arc_conjugate_to(p, nullptr, nullptr));
      f.word((uint64_t)gshgeom_path_arc_conjugate_to(nullptr, &fm, nullptr));
      f.word((uint64_t)gshgeom_path_arc_center_to(nullptr, P(0, 0), 1.0, 1.0, 0.0, GSHGEOM_ANGLE_DEGREES, 0.0, 90.0, nullptr));
      f.word((uint64_t)gshgeom_path_arc_tangent_to(nullptr, P(1, 1), P(2, 0), 1.0));
      f.word((uint64_t)gshgeom_path_arc_bulge_to(nullptr, P(1, 1), 1.0));
      fold_path(f, GSHGEOM_OK, p);
      GSHGEOM_Path * low = nullptr;
      f.word((uint64_t)gshgeom_path_lower_arcs(p, nullptr, nullptr, &low));
      gshgeom_path_destroy(low);
      low = nullptr;
      f.word((uint64_t)gshgeom_path_lower_arcs(nullptr, nullptr, nullptr, &low));
      f.word((uint64_t)gshgeom_path_lower_arcs(p, nullptr, nullptr, nullptr));
      gshgeom_path_destroy(p);
    }
  }
  f.word((uint64_t)gshgeom_arc_pieces(nullptr, pieces, &n));
  f.word((uint64_t)gshgeom_arc_pieces(&a, nullptr, &n));
  f.word((uint64_t)gshgeom_arc_pieces(&a, pieces, nullptr));
  f.word((uint64_t)gshgeom_arc_circle(nullptr, &c, &radius, &ccw));
  f.word((uint64_t)gshgeom_arc_circle(&a, nullptr, &radius, &ccw));
  f.word((uint64_t)gshgeom_arc_circle(&a, &c, nullptr, &ccw));
  f.word((uint64_t)gshgeom_arc_circle(&a, &c, &radius, nullptr));
  f.word((uint64_t)gshgeom_arc_map(nullptr, &a, &back));
  f.word((uint64_t)gshgeom_arc_map(&ident, nullptr, &back));
  f.word((uint64_t)gshgeom_arc_map(&ident, &a, nullptr));
  f.word((uint64_t)gshgeom_arc_point(nullptr, 0.0, &c));
  f.word((uint64_t)gshgeom_arc_point(&fm, 0.0, nullptr));
  GSHGEOM_Arc_Form bad = fm;
  bad.unit = (GSHGEOM_Angle_Unit)-1;
  f.word((uint64_t)gshgeom_arc_point(&bad, 0.0, &c));
}

/** The exact sums the arc routines decide with, converted to a double: ties, sticky bits, cancellation, the ends of the exponent range. */
void sum_grid(Fold & f) {
  uint64_t s = 0x9e3779b97f4a7c15ull;
  auto rnd = [&]() {
    s ^= s << 13;
    s ^= s >> 7;
    s ^= s << 17;
    return s;
  };
  auto num = [&]() {
    const uint64_t r = rnd();
    const int kind = (int)(r % 9);
    const int e = (int)((r >> 8) % 120) - 60;
    double m = 1.0 + (double)((rnd() >> 12) & ((1ull << 52) - 1)) * 0x1p-52;
    switch (kind) {
      case 0: return std::ldexp(m, e);
      case 1: return -std::ldexp(m, e);
      case 2: return std::ldexp(1.0, e);
      case 3: return std::ldexp(1.0 + 0x1p-52, e);
      case 4: return std::ldexp(1.0 + 0x1p-26, e);
      case 5: return std::ldexp(m, e - 900);
      case 6: return std::ldexp(m, e + 900);
      case 7: return 0.0;
      default: return std::ldexp(1.0 - 0x1p-53, e);
    }
  };
  for (int i = 0; i < 40000; i++) {
    GSHGEOM_ExactSum sum;
    gshgeom_exact_sum_init(&sum);
    const int terms = 1 + (int)(rnd() % 4);
    for (int t = 0; t < terms; t++) {
      double fac[3] = {num(), num(), num()};
      const int n = 1 + (int)(rnd() % 3);
      gshgeom_exact_sum_add(&sum, (int)(rnd() & 1), n, fac);
      if ((rnd() & 7) == 0) {  // the term again, negated, with a small one left over: cancellation
        gshgeom_exact_sum_add(&sum, 1, n, fac);
        double tiny[1] = {std::ldexp(1.0, (int)(rnd() % 200) - 100)};
        gshgeom_exact_sum_add(&sum, (int)(rnd() & 1), 1, tiny);
      }
    }
    f.dbl(gshgeom_exact_sum_to_double(&sum));
    f.word((uint64_t)(int64_t)gshgeom_exact_sum_sign(&sum));
  }
  // ties: 2^53 + 1, 2^53 + 3, 2^54 + 2 and the sticky bit far below
  for (int k = 0; k < 140; k++) {
    GSHGEOM_ExactSum sum;
    gshgeom_exact_sum_init(&sum);
    double big[1] = {std::ldexp(1.0, 53 + k % 3)};
    double odd[1] = {(double)(1 + 2 * (k % 5))};
    double sticky[1] = {std::ldexp(1.0, -(k + 1))};
    gshgeom_exact_sum_add(&sum, k & 1, 1, big);
    gshgeom_exact_sum_add(&sum, k & 1, 1, odd);
    if (k % 4 >= 2) {
      gshgeom_exact_sum_add(&sum, k & 1, 1, sticky);
    }
    f.dbl(gshgeom_exact_sum_to_double(&sum));
  }
}

}  // namespace

TEST(ArcGrid, EveryEntryPointOverTheSpecialAndOrdinaryArgumentsIsWhatItWas) {
  Fold f;
  endpoint_grid(f);
  const uint64_t after_endpoint = f.h;
  center_grid(f);
  const uint64_t after_center = f.h;
  tangent_bulge_grid(f);
  const uint64_t after_tangent = f.h;
  record_grid(f);
  const uint64_t after_records = f.h;
  refusal_grid(f);
  sum_grid(f);
  if (std::getenv("GSHGEOM_ARC_GRID_PRINT")) {
    std::printf("arc-grid: endpoint %016llx center %016llx tangent+bulge %016llx records %016llx all %016llx\n", (unsigned long long)after_endpoint,
        (unsigned long long)after_center, (unsigned long long)after_tangent, (unsigned long long)after_records, (unsigned long long)f.h);
  }
  EXPECT_EQ(f.h, kGridHash);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
