/**
 * @file
 *
 * The conic route (flatten.h, revision 3): the device path of a rational quadratic, a rational
 * cubic and an arc, and the proof above src/device/device.c that every point of the exact curve is
 * within the tolerance less 0.29 units of the emitted polyline, and that the depth cap never ends a
 * piece that is not flat.
 *
 * Pixels cannot test that, so these tests evaluate the curve exactly (the homogeneous Bernstein
 * form in `long double` at dyadic parameters, which are exact in the type; the only error is the
 * arithmetic of the sum, below a millionth of a unit even where `long double` is a plain `double`,
 * against the 0.0018 of margin the proof leaves) at thirty-three parameters inside every piece the
 * route made, on every size class up to the whole `int32` range, and assert the distance to that
 * piece's segment.  A piece-by-piece measure is used because a curve spanning the `int32` range is
 * 2^33 units long, so no affordable uniform sampling would see a unit of deviation between two
 * samples.  A curve whose weights make a piece deeper than 32 levels (a weight ratio of 2^200) has
 * no exact parameter interval for its pieces, and is measured against the whole polyline at a
 * dense set of parameters that includes 2^-k and 1 - 2^-k for every depth the proof allows.
 *
 * The depth cap is shown to matter in the other direction: the same curve under a cap a curve does
 * reach ends pieces that are not flat, the count says so, and the bound is lost.
 *
 * Copyright 2026 by Corey Pennycuff
 */

#include "test_helpers.h"

#include "../alloc_tracker.h"
#include "../curve_eval.h"

#include "../../src/device/conic_internal.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

namespace {

using Dp = GSHGEOM_Device_Point;
using curve_eval::Real;
using curve_eval::P2;

GSHGEOM_Point P(double x, double y) {
  return GSHGEOM_Point{x, y, 0.0};
}

struct Vertex {
  Dp at;
  uint32_t leaf;
  int level;
};

GSHGEOM_Result collect(void * user, Dp to, uint32_t leaf, int level) {
  static_cast<std::vector<Vertex> *>(user)->push_back({to, leaf, level});
  return GSHGEOM_OK;
}

/** A conic in device units: n control points, weights. */
struct Conic {
  int n = 3;
  double x[4] = {0, 0, 0, 0};
  double y[4] = {0, 0, 0, 0};
  double w[4] = {1, 1, 1, 1};
};

curve_eval::Elem elem_of(const Conic & c) {
  curve_eval::Elem e;
  e.n = c.n - 1;
  for (int i = 0; i < c.n; i++) {
    e.p[i] = {(Real)c.x[i], (Real)c.y[i]};
    e.w[i] = (Real)c.w[i];
  }
  return e;
}

/** The rule, floor(v + 1/2), in long double (exact for |v| < 2^31 with a double's mantissa). */
int32_t rounded(double v) {
  return (int32_t)std::floor((Real)v + 0.5L);
}

Real distance(P2 q, Dp a, Dp b) {
  return curve_eval::dist_to_segment(q, {(Real)a.x, (Real)a.y}, {(Real)b.x, (Real)b.y});
}

Real polyline_distance(P2 q, const std::vector<Dp> & pl) {
  Real best = 1e30L;
  for (size_t i = 0; i + 1 < pl.size(); i++) {
    best = std::min(best, distance(q, pl[i], pl[i + 1]));
  }
  return best;
}

struct Measure {
  Real worst = 0;       ///< the largest distance of a sampled point of the curve from the polyline
  Real at = -1;         ///< the parameter where it was seen
  size_t pieces = 0;
  int deepest = 0;
  int forced = 0;
  bool tiled = true;    ///< the pieces tile [0, 1] in order (where their intervals are exact)
  bool exact = true;    ///< every piece was measured against its own segment
  GSHGEOM_Result result = GSHGEOM_OK;
};

constexpr int kSamples = 32;

/** Route @p c at @p tolerance with the cap at @p levels and measure the exact curve against the polyline. */
Measure measure(const Conic & c, int32_t tolerance, int levels = GSHGEOM_CONIC_MAX_DEPTH) {
  Measure m;
  std::vector<Vertex> out;
  m.result = gshgeom_conic_route_capped(c.n, c.x, c.y, c.w, tolerance, levels, collect, &out, &m.forced, &m.deepest);
  if (m.result != GSHGEOM_OK) {
    return m;
  }
  m.pieces = out.size();
  const curve_eval::Elem e = elem_of(c);
  int deepest_level = 0;
  for (const Vertex & v : out) {
    deepest_level = std::max(deepest_level, v.level);
  }
  Dp from = {rounded(c.x[0]), rounded(c.y[0])};
  std::vector<Dp> pl;
  pl.push_back(from);
  for (const Vertex & v : out) {
    pl.push_back(v.at);
  }
  double wmin = c.w[0], wmax = c.w[0];
  for (int i = 1; i < c.n; i++) {
    wmin = std::min(wmin, c.w[i]);
    wmax = std::max(wmax, c.w[i]);
  }
  if (deepest_level <= 32) {
    // The pieces tile the parameter interval: in order, each starting where the last ended.
    Real expect = 0;
    for (const Vertex & v : out) {
      const Real span = std::ldexp(1.0L, v.level);
      const Real lo = (Real)v.leaf / span;
      if (lo != expect) {
        m.tiled = false;
      }
      expect = ((Real)v.leaf + 1.0L) / span;
      for (int i = 0; i <= kSamples; i++) {
        const Real s = ((Real)v.leaf + (Real)i / kSamples) / span;
        const Real d = distance(curve_eval::eval(e, s), from, v.at);
        if (d > m.worst) {
          m.worst = d;
          m.at = s;
        }
      }
      from = v.at;
    }
    m.tiled = m.tiled && expect == 1.0L;
  } else {
    m.exact = false;
  }
  if (deepest_level <= 32 && wmax / wmin <= 4096.0) {
    return m;
  }
  // Weights that are far apart put the curve's length into a sliver of the parameter, where a piece's own samples do not look (and
  // a piece deeper than 32 levels has no exact interval): the whole polyline at a dense set of parameters, among them 2^-k and
  // 1 - 2^-k at every depth the proof allows, a few multiples of each.
  curve_eval::Elem rev = e;
  for (int i = 0; i < e.n + 1; i++) {
    rev.p[i] = e.p[e.n - i];
    rev.w[i] = e.w[e.n - i];
  }
  auto see = [&](const curve_eval::Elem & el, Real s, bool reversed) {
    const Real d = polyline_distance(curve_eval::eval(el, s), pl);
    if (d > m.worst) {
      m.worst = d;
      m.at = reversed ? 1.0L - s : s;
    }
  };
  if (pl.size() <= 4000) {
    for (int i = 0; i <= 4096; i++) {
      see(e, (Real)i / 4096.0L, false);
    }
  }
  for (int k = 1; k <= 160; k++) {
    const Real s = std::ldexp(1.0L, -k);
    for (const Real mult : {1.0L, 1.25L, 1.5L, 1.75L}) {
      if (mult * s <= 0.5L) {
        see(e, mult * s, false);
        see(rev, mult * s, true);
      }
    }
    if (k >= 2) {
      see(e, 3 * s, false);
      see(rev, 3 * s, true);
    }
  }
  return m;
}

std::string describe(const Conic & c) {
  char buf[512];
  std::string out;
  for (int i = 0; i < c.n; i++) {
    std::snprintf(buf, sizeof buf, " (%.17g, %.17g; w %.17g)", c.x[i], c.y[i], c.w[i]);
    out += buf;
  }
  return out;
}

/** The proof's bound for the tolerance, with a little for the evaluation's own arithmetic. */
Real bound(int32_t tolerance) {
  return (Real)tolerance - 0.29L + 1e-6L;
}

struct Random {
  uint64_t state;
  uint64_t next() {
    state ^= state >> 12;
    state ^= state << 25;
    state ^= state >> 27;
    return state * UINT64_C(0x2545F4914F6CDD1D);
  }
  /** A device coordinate in [-2^bits, 2^bits) with a fraction, that rounds into int32. */
  double coordinate(int bits) {
    const uint64_t span = (uint64_t)1 << (bits + 1);
    double v = (double)((int64_t)(next() % span) - ((int64_t)1 << bits)) + (double)(next() % 1024u) / 1024.0;
    if (v >= 2147483647.0) {
      v = 2147483646.5;
    }
    return v;
  }
  /** A weight of the class: 0 about 1, 1 from 2^-10 to 2^10, 2 from 2^-40 to 2^40, 3 the extremes of what a path accepts. */
  double weight(int klass) {
    switch (klass) {
    case 0:
      return 0.25 + (double)(next() % 4096u) / 1024.0;
    case 1:
      return std::ldexp(1.0 + (double)(next() % 1024u) / 1024.0, (int)(next() % 21u) - 10);
    case 2:
      return std::ldexp(1.0 + (double)(next() % 1024u) / 1024.0, (int)(next() % 81u) - 40);
    default: {
      const uint64_t pick = next() % 4u;
      return pick == 0 ? 0x1p100 : pick == 1 ? 0x1p-100 : pick == 2 ? 1.0 : std::ldexp(1.0, (int)(next() % 201u) - 100);
    }
    }
  }
};

Conic random_conic(Random & r, int n, int bits, int wclass) {
  Conic c;
  c.n = n;
  for (int i = 0; i < n; i++) {
    c.x[i] = r.coordinate(bits);
    c.y[i] = r.coordinate(bits);
    c.w[i] = r.weight(wclass);
  }
  return c;
}

}  // namespace

// The theorem, at every size a conic can have, at three tolerances and four classes of weights.
TEST(ConicBound, EveryPieceOfEveryConicIsWithinTheToleranceLessTwoTenthsOfTheCurve) {
  GSHGEOM_SKIP_UNDER_VALGRIND();
  Random random{UINT64_C(0x9E3779B97F4A7C15)};
  const int classes[] = {2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 31};
  const int32_t tolerances[] = {2, 4, 9};
  size_t most = 0;
  int deepest_ordinary = 0;
  for (const int bits : classes) {
    const int curves = bits <= 14 ? 6 : (bits <= 24 ? 2 : 1);
    Real worst_margin = 1e9L;
    for (int wclass = 0; wclass < 4; wclass++) {
      for (const int32_t t : tolerances) {
        for (int i = 0; i < curves; i++) {
          for (int n = 3; n <= 4; n++) {
            const Conic c = random_conic(random, n, bits, wclass);
            const Measure m = measure(c, t);
            ASSERT_EQ(m.result, GSHGEOM_OK);
            EXPECT_EQ(m.forced, 0) << "the depth cap ended a piece that was not flat";
            EXPECT_TRUE(m.tiled) << "the pieces do not tile the parameter interval";
            EXPECT_LE(m.worst, bound(t)) << "n=" << n << " bits=" << bits << " wclass=" << wclass << " t=" << t
                                         << " strays " << (double)m.worst << " units at s=" << (double)m.at << " pieces=" << m.pieces
                                         << " deepest=" << m.deepest << " " << describe(c);
            EXPECT_LE(m.deepest, GSHGEOM_CONIC_MAX_DEPTH);
            if (wclass < 2) {
              deepest_ordinary = std::max(deepest_ordinary, m.deepest);
            }
            worst_margin = std::min(worst_margin, (Real)t - m.worst);
            most = std::max(most, m.pieces);
          }
        }
      }
    }
    RecordProperty("smallest_margin_below_2_pow_" + std::to_string(bits), std::to_string((double)worst_margin));
  }
  EXPECT_GE(most, 10000u) << "the large classes make the long polylines";
  EXPECT_LE(deepest_ordinary, 40) << "an ordinary curve is flat by about twenty levels";
  RecordProperty("most_pieces_of_one_conic", std::to_string(most));
}

// The corners of the int32 square, with the weights of every extreme: a piece count that is bounded by what the route does, and a bound that holds.
TEST(ConicBound, TheCornersOfTheInt32SquareAndTheExtremeWeights) {
  GSHGEOM_SKIP_UNDER_VALGRIND();
  Random random{UINT64_C(0x0123456789ABCDEF)};
  size_t most = 0;
  for (int i = 0; i < 12; i++) {
    Conic c;
    c.n = 3 + (i & 1);
    for (int k = 0; k < c.n; k++) {
      c.x[k] = (random.next() & 1u) ? 2147483647.25 : -2147483648.25;
      c.y[k] = (random.next() & 1u) ? 2147483647.25 : -2147483648.25;
      c.w[k] = random.weight(3);
    }
    for (const int32_t t : {2, 4}) {
      const Measure m = measure(c, t);
      ASSERT_EQ(m.result, GSHGEOM_OK);
      EXPECT_EQ(m.forced, 0);
      EXPECT_LE(m.worst, bound(t)) << i << " t=" << t;
      most = std::max(most, m.pieces);
    }
  }
  // The route emits one line per leaf and the device path holds at most max_elements of them: a curve is never more than
  // that, and the work is bounded by it (the next test).
  EXPECT_LE(most, (size_t)GSHGEOM_LIMITS_DEFAULT_MAX_ELEMENTS);
  RecordProperty("most_pieces_of_a_corner_conic", std::to_string(most));
}

// A weight ratio of 2^200 is what makes the proof need 236 levels: a uniform split of such a curve spends about a hundred levels before
// it divides the curve's length at all.  The real cap holds it, a cap of twenty (the integer flattener's) does not, and the bound is then lost.
TEST(ConicDepth, AWeightRatioOfTwoToThe200NeedsMoreThanAHundredLevelsAndTheCapHoldsIt) {
  GSHGEOM_SKIP_UNDER_VALGRIND();
  Conic c;
  c.n = 3;
  c.x[0] = 0;       c.y[0] = 0;
  c.x[1] = 64000;   c.y[1] = 64000;
  c.x[2] = 128000;  c.y[2] = 0;
  c.w[0] = 0x1p-100;
  c.w[1] = 1.0;
  c.w[2] = 0x1p100;
  const int32_t t = GSHGEOM_FLATTEN_TOLERANCE;
  const Measure whole = measure(c, t);
  ASSERT_EQ(whole.result, GSHGEOM_OK);
  EXPECT_EQ(whole.forced, 0);
  EXPECT_GT(whole.deepest, 100) << "the split spends its levels where the weight is";
  EXPECT_LE(whole.deepest, 236) << "and the proof needs no more than 236";
  EXPECT_LE(whole.worst, bound(t));
  // The mirror image of the same curve.
  Conic r = c;
  std::swap(r.x[0], r.x[2]);
  std::swap(r.w[0], r.w[2]);
  const Measure mirror = measure(r, t);
  ASSERT_EQ(mirror.result, GSHGEOM_OK);
  EXPECT_EQ(mirror.forced, 0);
  EXPECT_GT(mirror.deepest, 100);
  EXPECT_LE(mirror.worst, bound(t));
  // The cap lowered to twenty ends pieces that are not flat, says so, and the bound is lost: the check sees a violation.
  const Measure capped = measure(c, t, 20);
  ASSERT_EQ(capped.result, GSHGEOM_OK);
  EXPECT_GT(capped.forced, 0);
  EXPECT_GT(capped.worst, bound(t)) << "a cap that binds loses the guarantee, which the bound check sees";
  // And the public call, at the real cap, never reports a forced piece (it would be an internal error).
  std::vector<Vertex> out;
  EXPECT_EQ(gshgeom_conic_route(c.n, c.x, c.y, c.w, t, collect, &out), GSHGEOM_OK);
}

// The cap that never binds for an ordinary curve is a long way from where it would: a curve whose weights are within 2^10 of each other
// is flat by the integer flattener's twenty-odd levels, so the sweep of caps below ends pieces and the real one does not.
TEST(ConicDepth, AnOrdinaryCurveNeedsAboutTwentyLevelsAndAnyCapAboveIsTheSameRoute) {
  Random random{UINT64_C(0xD1B54A32D192ED03)};
  const Conic c = random_conic(random, 4, 28, 1);
  const Measure real = measure(c, 4);
  ASSERT_EQ(real.result, GSHGEOM_OK);
  EXPECT_EQ(real.forced, 0);
  const int needed = real.deepest;
  EXPECT_LE(needed, 40);
  for (int cap = needed; cap <= GSHGEOM_CONIC_MAX_DEPTH; cap += 37) {
    const Measure m = measure(c, 4, cap);
    EXPECT_EQ(m.forced, 0) << cap;
    EXPECT_EQ(m.pieces, real.pieces) << "a cap above the depth reached changes nothing, " << cap;
  }
  const Measure low = measure(c, 4, needed - 1);
  EXPECT_GT(low.forced, 0) << "one level less ends a piece that is not flat";
}

TEST(ConicRoute, RefusesWhatItCannotRoute) {
  Conic c;
  c.n = 3;
  c.x[1] = 10;
  c.y[1] = 10;
  c.x[2] = 20;
  std::vector<Vertex> out;
  int forced = 0;
  EXPECT_EQ(gshgeom_conic_route_capped(2, c.x, c.y, c.w, 4, 10, collect, &out, &forced, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_conic_route_capped(5, c.x, c.y, c.w, 4, 10, collect, &out, &forced, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_conic_route_capped(3, c.x, c.y, c.w, 1, 10, collect, &out, &forced, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_conic_route_capped(3, c.x, c.y, c.w, 4, -1, collect, &out, &forced, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_conic_route_capped(3, c.x, c.y, c.w, 4, GSHGEOM_CONIC_MAX_DEPTH + 1, collect, &out, &forced, nullptr),
      GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_conic_route_capped(3, nullptr, c.y, c.w, 4, 10, collect, &out, &forced, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_conic_route_capped(3, c.x, c.y, c.w, 4, 10, nullptr, &out, &forced, nullptr), GSHGEOM_ERR_INVALID);
  // A vertex that does not round into int32 is the caller's ERR_LIMIT, the route's own end point included.
  c.x[2] = 2147483648.0;
  EXPECT_EQ(gshgeom_conic_route(3, c.x, c.y, c.w, 4, collect, &out), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(GSHGEOM_CONIC_MAX_DEPTH, 240);
}

// ---------------------------------------------------------------------------
// The device path of a model path: rational curves and arcs, through the public calls
// ---------------------------------------------------------------------------

namespace {

GSHGEOM_Path * model() {
  GSHGEOM_Path * p = nullptr;
  EXPECT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &p), GSHGEOM_OK);
  return p;
}

GSHGEOM_Map matrix(double a, double b, double c, double d, double tx = 0, double ty = 0) {
  const double m[4] = {a, b, c, d};
  const double t[2] = {tx, ty};
  GSHGEOM_Map map;
  EXPECT_EQ(gshgeom_map_from_matrix(2, m, t, &map), GSHGEOM_OK);
  return map;
}

/** The points of a device path of a move and lines. */
std::vector<Dp> polyline_of(const GSHGEOM_Device_Path * dp) {
  std::vector<Dp> pts;
  const size_t n = gshgeom_device_path_element_count(dp);
  for (size_t i = 0; i < n; i++) {
    GSHGEOM_Device_Element e;
    EXPECT_EQ(gshgeom_device_path_element(dp, i, &e), GSHGEOM_OK);
    EXPECT_TRUE(e.verb == GSHGEOM_DEVICE_MOVE || e.verb == GSHGEOM_DEVICE_LINE) << "verb " << e.verb;
    pts.push_back(e.p[0]);
  }
  return pts;
}

/** The exact point of an arc's centre form at fraction tau of its sweep, mapped to device units, in long double. */
P2 ellipse_device_point(const GSHGEOM_Arc & a, const GSHGEOM_Map & map, Real tau) {
  const P2 q = curve_eval::arc_point(a, tau);
  return {((Real)map.m[0][0] * q.x + (Real)map.m[0][1] * q.y + (Real)map.t[0]) * 64.0L,
      ((Real)map.m[1][0] * q.x + (Real)map.m[1][1] * q.y + (Real)map.t[1]) * 64.0L};
}

/** The largest distance of the exact ellipse arc, mapped, from the polyline, over @p samples parameters. */
Real arc_deviation(const GSHGEOM_Arc & a, const GSHGEOM_Map & map, const std::vector<Dp> & pl, int samples) {
  Real worst = 0;
  size_t hint = 0;
  for (int i = 0; i <= samples; i++) {
    const P2 q = ellipse_device_point(a, map, (Real)i / samples);
    // The nearest segment is usually near the last one's: look there first, and everywhere only when that is not within the bound.
    Real best = 1e30L;
    size_t where = hint;
    const size_t lo = hint > 8 ? hint - 8 : 0;
    const size_t hi = std::min(pl.size() - 1, hint + 96);
    for (size_t k = lo; k < hi; k++) {
      const Real d = distance(q, pl[k], pl[k + 1]);
      if (d < best) {
        best = d;
        where = k;
      }
    }
    if (best > 3.0L) {
      for (size_t k = 0; k + 1 < pl.size(); k++) {
        const Real d = distance(q, pl[k], pl[k + 1]);
        if (d < best) {
          best = d;
          where = k;
        }
      }
    }
    hint = where;
    worst = std::max(worst, best);
  }
  return worst;
}

}  // namespace

// Arcs under maps that rotate, skew, mirror, stretch unevenly and collapse, of every sweep and at sizes from a few pixels to the int32
// range, against the exact ellipse (not its pieces: the decomposition's own error is in what is asserted).
TEST(ConicArc, AnArcUnderAnyMapIsWithinTheToleranceOfTheExactEllipseAndEndsWhereItShould) {
  GSHGEOM_SKIP_UNDER_VALGRIND();
  Random random{UINT64_C(0xA24BAED4963EE407)};
  const double sweeps[] = {5.0, 45.0, 90.0, 91.0, 180.0, 270.0, -200.0, 360.0};
  const int shifts[] = {-6, 0, 6, 12, 18};
  size_t most = 0;
  for (const int shift : shifts) {
    const double k = std::ldexp(1.0, shift);
    const GSHGEOM_Map maps[] = {
        matrix(k, 0, 0, k, 3 * k, -2 * k),
        matrix(0.8660254037844386 * k, -0.5 * k, 0.5 * k, 0.8660254037844386 * k, 1 * k, 1 * k),
        matrix(k, 0.5 * k, 0, k),
        matrix(-k, 0, 0, k),
        matrix(3 * k, 0, 0, 0.25 * k),
        matrix(k, 2 * k, 2 * k, 4 * k),  // singular: the ellipse is a segment
    };
    for (size_t mi = 0; mi < sizeof maps / sizeof maps[0]; mi++) {
      const bool big = shift >= 18;
      if (big && mi != 1 && mi != 2) {
        continue;
      }
      for (const double sweep : sweeps) {
        if (big && sweep != 360.0 && sweep != 91.0) {
          continue;
        }
        GSHGEOM_Path * p = model();
        const double ra = 2.0 + (double)(random.next() % 11u);
        const double rb = (random.next() % 3u == 0u) ? ra : 1.0 + (double)(random.next() % 11u);
        const double phi = (double)(random.next() % 360u) * 3.14159265358979323846 / 180.0;
        GSHGEOM_Arc_Form f;
        f.center = P(0.0625 * (double)(random.next() % 64u), -0.125 * (double)(random.next() % 32u));
        f.u = P(ra * std::cos(phi), ra * std::sin(phi));
        f.v = P(-rb * std::sin(phi), rb * std::cos(phi));
        f.unit = GSHGEOM_ANGLE_DEGREES;
        f.t0 = (double)(random.next() % 360u);
        f.sweep = sweep;
        f.flags = 0;
        ASSERT_EQ(gshgeom_path_arc_conjugate_to(p, &f, nullptr), GSHGEOM_OK);
        GSHGEOM_Arc a;
        ASSERT_EQ(gshgeom_path_arc(p, 1, &a), GSHGEOM_OK);
        GSHGEOM_Device_Path * dp = nullptr;
        const GSHGEOM_Result r = gshgeom_device_path_from_model(p, &maps[mi], nullptr, nullptr, &dp);
        ASSERT_EQ(r, GSHGEOM_OK) << shift << " " << mi << " " << sweep;
        const std::vector<Dp> pl = polyline_of(dp);
        // The ends are the arc's own ends, mapped and rounded once.
        Dp s0, e0;
        ASSERT_EQ(gshgeom_device_point(&maps[mi], a.start, &s0), GSHGEOM_OK);
        ASSERT_EQ(gshgeom_device_point(&maps[mi], a.end, &e0), GSHGEOM_OK);
        ASSERT_GE(pl.size(), 2u);
        EXPECT_TRUE(pl.front().x == s0.x && pl.front().y == s0.y) << "the start";
        EXPECT_TRUE(pl.back().x == e0.x && pl.back().y == e0.y) << "the end";
        if (sweep == 360.0) {
          EXPECT_TRUE(pl.front().x == pl.back().x && pl.front().y == pl.back().y) << "a full turn closes in device space";
        }
        const Real d = arc_deviation(a, maps[mi], pl, shift >= 12 ? 1500 : 3000);
        EXPECT_LE(d, 4.0L - 0.12L) << "shift " << shift << " map " << mi << " sweep " << sweep;
        most = std::max(most, pl.size());
        gshgeom_device_path_destroy(dp);
        gshgeom_path_destroy(p);
      }
    }
  }
  EXPECT_GE(most, 5000u) << "the largest arcs are many lines";
}

// The same arc through the tolerance stated: a smaller one is more lines and closer, and below the minimum is refused.
TEST(ConicArc, TheToleranceIsTheCallersAndTheMinimumIsTwo) {
  GSHGEOM_SKIP_UNDER_VALGRIND();
  GSHGEOM_Path * p = model();
  ASSERT_EQ(gshgeom_path_move_to(p, P(500, 0)), GSHGEOM_OK);
  GSHGEOM_Arc_Form f = {P(0, 0), P(500, 0), P(0, 500), GSHGEOM_ANGLE_DEGREES, 0, 270, GSHGEOM_ARC_CIRCULAR};
  ASSERT_EQ(gshgeom_path_arc_conjugate_to(p, &f, nullptr), GSHGEOM_OK);
  GSHGEOM_Arc a;
  ASSERT_EQ(gshgeom_path_arc(p, 1, &a), GSHGEOM_OK);
  const GSHGEOM_Map map = matrix(1, 0, 0, 1);
  size_t previous = 0;
  for (const int32_t t : {64, 16, 4, 2}) {
    GSHGEOM_Device_Path * dp = nullptr;
    ASSERT_EQ(gshgeom_device_path_from_model_with_tolerance(p, &map, t, nullptr, nullptr, &dp), GSHGEOM_OK);
    const std::vector<Dp> pl = polyline_of(dp);
    EXPECT_GT(pl.size(), previous) << t;
    previous = pl.size();
    EXPECT_LE(arc_deviation(a, map, pl, 4000), (Real)t - 0.29L + 0.001L) << t;
    gshgeom_device_path_destroy(dp);
  }
  GSHGEOM_Device_Path * dp = reinterpret_cast<GSHGEOM_Device_Path *>(0x1);
  EXPECT_EQ(gshgeom_device_path_from_model_with_tolerance(p, &map, 1, nullptr, nullptr, &dp), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_device_path_from_model_with_tolerance(p, &map, 0, nullptr, nullptr, &dp), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_device_path_from_model_with_tolerance(p, &map, INT32_MIN, nullptr, nullptr, &dp), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(dp, reinterpret_cast<GSHGEOM_Device_Path *>(0x1)) << "written only on success";
  // The call without a tolerance is the call at the default.
  GSHGEOM_Device_Path * one = nullptr;
  GSHGEOM_Device_Path * two = nullptr;
  ASSERT_EQ(gshgeom_device_path_from_model(p, &map, nullptr, nullptr, &one), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_device_path_from_model_with_tolerance(p, &map, GSHGEOM_FLATTEN_TOLERANCE, nullptr, nullptr, &two), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_device_path_hash(one), gshgeom_device_path_hash(two));
  gshgeom_device_path_destroy(one);
  gshgeom_device_path_destroy(two);
  gshgeom_path_destroy(p);
}

// A rounded rectangle is four lines and four quarter arcs; closed, it closes in device space under any map, and its device path is the
// lines of the same path with the arcs lowered to their pieces.
TEST(ConicArc, ARoundedRectangleClosesAndIsTheLoweredPathsLines) {
  GSHGEOM_Path * p = model();
  const double x0 = 3, y0 = 5, w = 40, h = 24, r = 6;
  ASSERT_EQ(gshgeom_path_move_to(p, P(x0 + r, y0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(p, P(x0 + w - r, y0)), GSHGEOM_OK);
  GSHGEOM_Arc_Form f = {P(x0 + w - r, y0 + r), P(r, 0), P(0, r), GSHGEOM_ANGLE_DEGREES, -90, 90, GSHGEOM_ARC_CIRCULAR};
  ASSERT_EQ(gshgeom_path_arc_conjugate_to(p, &f, nullptr), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(p, P(x0 + w, y0 + h - r)), GSHGEOM_OK);
  f.center = P(x0 + w - r, y0 + h - r);
  f.t0 = 0;
  ASSERT_EQ(gshgeom_path_arc_conjugate_to(p, &f, nullptr), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(p, P(x0 + r, y0 + h)), GSHGEOM_OK);
  f.center = P(x0 + r, y0 + h - r);
  f.t0 = 90;
  ASSERT_EQ(gshgeom_path_arc_conjugate_to(p, &f, nullptr), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(p, P(x0, y0 + r)), GSHGEOM_OK);
  f.center = P(x0 + r, y0 + r);
  f.t0 = 180;
  ASSERT_EQ(gshgeom_path_arc_conjugate_to(p, &f, nullptr), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_close(p), GSHGEOM_OK);
  GSHGEOM_Path * lowered = nullptr;
  ASSERT_EQ(gshgeom_path_lower_arcs(p, nullptr, nullptr, &lowered), GSHGEOM_OK);
  const GSHGEOM_Map maps[] = {matrix(2, 0, 0, 2), matrix(1.5, 0.75, -0.25, 3, 10, 10), matrix(0.8660254037844386 * 9, -4.5, 4.5, 0.8660254037844386 * 9)};
  for (const GSHGEOM_Map & map : maps) {
    GSHGEOM_Device_Path * dp = nullptr;
    GSHGEOM_Device_Path * dl = nullptr;
    ASSERT_EQ(gshgeom_device_path_from_model(p, &map, nullptr, nullptr, &dp), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_device_path_from_model(lowered, &map, nullptr, nullptr, &dl), GSHGEOM_OK);
    EXPECT_EQ(gshgeom_device_path_hash(dp), gshgeom_device_path_hash(dl));
    const size_t n = gshgeom_device_path_element_count(dp);
    GSHGEOM_Device_Element first, last_line, close;
    ASSERT_EQ(gshgeom_device_path_element(dp, 0, &first), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_device_path_element(dp, n - 2, &last_line), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_device_path_element(dp, n - 1, &close), GSHGEOM_OK);
    EXPECT_EQ(close.verb, GSHGEOM_DEVICE_CLOSE);
    EXPECT_TRUE(last_line.p[0].x == first.p[0].x && last_line.p[0].y == first.p[0].y)
        << "the last arc ends where the subpath began, in device space";
    gshgeom_device_path_destroy(dp);
    gshgeom_device_path_destroy(dl);
  }
  gshgeom_path_destroy(p);
  gshgeom_path_destroy(lowered);
}

// CV6: hostile input.  Never a hang, never a depth that binds, ERR_LIMIT where a coordinate does not round into int32.
TEST(ConicHostile, LimitsAreLimitsAndTheWorkIsBounded) {
  const GSHGEOM_Map id = matrix(1, 0, 0, 1);
  // A control point at the map's limit, and a weight at either extreme of what the path accepts.
  {
    GSHGEOM_Path * p = model();
    ASSERT_EQ(gshgeom_path_move_to(p, P(0, 0)), GSHGEOM_OK);
    const double w[3] = {0x1p100, 1.0, 0x1p-100};
    ASSERT_EQ(gshgeom_path_rational_quadratic_to(p, P(0x1p99, 0), P(10, 10), w), GSHGEOM_OK);
    GSHGEOM_Device_Path * dp = nullptr;
    EXPECT_EQ(gshgeom_device_path_from_model(p, &id, nullptr, nullptr, &dp), GSHGEOM_ERR_LIMIT);
    EXPECT_EQ(dp, nullptr);
    // A map that brings the same curve into range: the limit is the device range, not the model's.
    const GSHGEOM_Map small = matrix(0x1p-80, 0, 0, 0x1p-80);
    ASSERT_EQ(gshgeom_device_path_from_model(p, &small, nullptr, nullptr, &dp), GSHGEOM_OK);
    EXPECT_GE(gshgeom_device_path_element_count(dp), 3u);
    gshgeom_device_path_destroy(dp);
    gshgeom_path_destroy(p);
  }
  // A curve across the int32 range under a small element limit is refused, and quickly: the work is bounded by the limit and the cap.
  {
    GSHGEOM_Path * p = model();
    ASSERT_EQ(gshgeom_path_move_to(p, P(-3e7, -3e7)), GSHGEOM_OK);
    const double w[4] = {1.0, 4.0, 0.25, 1.0};
    ASSERT_EQ(gshgeom_path_rational_cubic_to(p, P(3e7, -3e7), P(-3e7, 3e7), P(3e7, 3e7), w), GSHGEOM_OK);
    const GSHGEOM_Map big = matrix(64, 0, 0, 64);  // 64 model units are 64 * 64 device units
    GSHGEOM_Limits l;
    gshgeom_limits_default(&l);
    l.max_elements = 100;
    GSHGEOM_Device_Path * dp = nullptr;
    EXPECT_EQ(gshgeom_device_path_from_model(p, &big, nullptr, &l, &dp), GSHGEOM_ERR_LIMIT) << "3e7 * 64 is past the device range";
    const GSHGEOM_Map fit = matrix(0.5, 0, 0, 0.5);  // 3e7 * 0.5 * 64 = 9.6e8 units, in range
    EXPECT_EQ(gshgeom_device_path_from_model(p, &fit, nullptr, &l, &dp), GSHGEOM_ERR_LIMIT) << "100 elements are too few";
    EXPECT_EQ(dp, nullptr);
    l.max_elements = 0;
    ASSERT_EQ(gshgeom_device_path_from_model(p, &fit, nullptr, &l, &dp), GSHGEOM_OK);
    EXPECT_LT(gshgeom_device_path_element_count(dp), (size_t)GSHGEOM_LIMITS_DEFAULT_MAX_ELEMENTS);
    gshgeom_device_path_destroy(dp);
    gshgeom_path_destroy(p);
  }
}

// An allocation that fails at any point of the device path of arcs and conics leaves nothing behind.
TEST(ConicHostile, AllocationFailureAtEveryCallOfAnArcAndConicPath) {
  GSHGEOM_Path * p = model();
  ASSERT_EQ(gshgeom_path_move_to(p, P(10, 0)), GSHGEOM_OK);
  GSHGEOM_Arc_Form f = {P(0, 0), P(10, 0), P(0, 10), GSHGEOM_ANGLE_DEGREES, 0, 200, GSHGEOM_ARC_CIRCULAR};
  ASSERT_EQ(gshgeom_path_arc_conjugate_to(p, &f, nullptr), GSHGEOM_OK);
  const double w[3] = {1.0, 2.0, 1.0};
  ASSERT_EQ(gshgeom_path_rational_quadratic_to(p, P(-30, 30), P(-20, -20), w), GSHGEOM_OK);
  const GSHGEOM_Map map = matrix(8, 0, 0, 8);
  tracker::Tracker count;
  {
    GSHGEOM_Allocator a = tracker::allocator_for(&count);
    GSHGEOM_Device_Path * dp = nullptr;
    ASSERT_EQ(gshgeom_device_path_from_model(p, &map, &a, nullptr, &dp), GSHGEOM_OK);
    gshgeom_device_path_destroy(dp);
  }
  ASSERT_GT(count.calls, 2);
  for (long k = 0; k < count.calls; k++) {
    tracker::Tracker tr;
    tr.fail_from = k;
    tr.persistent = true;
    GSHGEOM_Allocator a = tracker::allocator_for(&tr);
    GSHGEOM_Device_Path * dp = nullptr;
    EXPECT_EQ(gshgeom_device_path_from_model(p, &map, &a, nullptr, &dp), GSHGEOM_ERR_OOM) << k;
    EXPECT_EQ(dp, nullptr);
    EXPECT_TRUE(tr.live.empty()) << "k=" << k << " leaked";
  }
  gshgeom_path_destroy(p);
}

// Determinism: the route is a function of the control points, the weights and the tolerance.
TEST(ConicRoute, TheSameInputGivesTheSameLinesAndTheLinesDoNotDependOnTheCall) {
  Random random{UINT64_C(0x5851F42D4C957F2D)};
  for (int i = 0; i < 20; i++) {
    const Conic c = random_conic(random, 3 + (i & 1), 12 + i, i % 4);
    std::vector<Vertex> a, b;
    ASSERT_EQ(gshgeom_conic_route(c.n, c.x, c.y, c.w, 4, collect, &a), GSHGEOM_OK);
    int forced = -1, deepest = -1;
    ASSERT_EQ(gshgeom_conic_route_capped(c.n, c.x, c.y, c.w, 4, GSHGEOM_CONIC_MAX_DEPTH, collect, &b, &forced, &deepest), GSHGEOM_OK);
    ASSERT_EQ(a.size(), b.size());
    for (size_t k = 0; k < a.size(); k++) {
      EXPECT_TRUE(a[k].at.x == b[k].at.x && a[k].at.y == b[k].at.y && a[k].leaf == b[k].leaf && a[k].level == b[k].level);
    }
    EXPECT_EQ(forced, 0);
    // The last vertex is the curve's end point, rounded once.
    EXPECT_EQ(a.back().at.x, rounded(c.x[c.n - 1]));
    EXPECT_EQ(a.back().at.y, rounded(c.y[c.n - 1]));
  }
}

// A piece whose interior control point is exactly delta = t - 1 units from its chord is flat: the criterion is "at most", so the
// route emits the curve's end and nothing else, at every tolerance; a control point a hair past it splits.
TEST(ConicRoute, ABoundaryPieceIsFlatAndOneHairPastItIsNot) {
  for (const int32_t t : {2, 3, 4, 7, 16}) {
    Conic c;
    c.n = 3;
    c.x[1] = 64;
    c.x[2] = 128;
    c.y[1] = (double)(t - 1);
    std::vector<Vertex> out;
    ASSERT_EQ(gshgeom_conic_route(3, c.x, c.y, c.w, t, collect, &out), GSHGEOM_OK);
    ASSERT_EQ(out.size(), 1u) << t;
    EXPECT_EQ(out[0].at.x, 128);
    EXPECT_EQ(out[0].at.y, 0);
    EXPECT_EQ(out[0].level, 0);
    c.y[1] = (double)(t - 1) + 0x1p-40;
    out.clear();
    ASSERT_EQ(gshgeom_conic_route(3, c.x, c.y, c.w, t, collect, &out), GSHGEOM_OK);
    EXPECT_GT(out.size(), 1u) << t;
  }
}

// The arguments the route reads, and the out parameters it writes on every call.
TEST(ConicRoute, EveryArgumentIsCheckedAndTheCountsAreWrittenOnEveryCall) {
  Conic c;
  c.n = 3;
  c.x[1] = 10;
  c.y[1] = 10;
  c.x[2] = 20;
  std::vector<Vertex> out;
  EXPECT_EQ(gshgeom_conic_route(3, c.x, c.y, nullptr, 4, collect, &out), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_conic_route(3, c.x, nullptr, c.w, 4, collect, &out), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_conic_route(3, c.x, c.y, c.w, 4, nullptr, &out), GSHGEOM_ERR_INVALID);
  EXPECT_TRUE(out.empty());
  // Each coordinate of the end point is range checked on its own (a vertex that does not round into int32 is the caller's limit).
  Conic big = c;
  big.y[2] = 2147483648.0;
  EXPECT_EQ(gshgeom_conic_route(3, big.x, big.y, big.w, 4, collect, &out), GSHGEOM_ERR_LIMIT);
  big = c;
  big.x[2] = -2147483649.0;
  EXPECT_EQ(gshgeom_conic_route(3, big.x, big.y, big.w, 4, collect, &out), GSHGEOM_ERR_LIMIT);
  big = c;
  big.y[2] = -2147483649.0;
  EXPECT_EQ(gshgeom_conic_route(3, big.x, big.y, big.w, 4, collect, &out), GSHGEOM_ERR_LIMIT);
  // The counts are reset by every call, not accumulated or left as they were.
  int forced = 7, deepest = 9;
  ASSERT_EQ(gshgeom_conic_route_capped(3, c.x, c.y, c.w, 4, GSHGEOM_CONIC_MAX_DEPTH, collect, &out, &forced, &deepest), GSHGEOM_OK);
  EXPECT_EQ(forced, 0);
  EXPECT_LE(deepest, 3) << "a gentle piece is flat after a few halvings";
  forced = 7;
  deepest = 9;
  Conic flat;
  flat.n = 3;
  flat.x[1] = 64;
  flat.x[2] = 128;
  ASSERT_EQ(gshgeom_conic_route_capped(3, flat.x, flat.y, flat.w, 4, GSHGEOM_CONIC_MAX_DEPTH, collect, &out, &forced, &deepest), GSHGEOM_OK);
  EXPECT_EQ(forced, 0);
  EXPECT_EQ(deepest, 0);
  // A call that is refused for its arguments leaves nothing half done.
  forced = 7;
  deepest = 9;
  EXPECT_EQ(gshgeom_conic_route_capped(2, c.x, c.y, c.w, 4, 10, collect, &out, &forced, &deepest), GSHGEOM_ERR_INVALID);
  // The cap is the caller's to lower down to zero: a piece at depth zero is ended as it stands, and counted when it is not flat.
  EXPECT_EQ(gshgeom_conic_route_capped(3, c.x, c.y, c.w, 4, 0, collect, &out, &forced, &deepest), GSHGEOM_OK);
  EXPECT_EQ(deepest, 0);
  c.y[1] = 1e6;
  EXPECT_EQ(gshgeom_conic_route_capped(3, c.x, c.y, c.w, 4, 0, collect, &out, &forced, &deepest), GSHGEOM_OK);
  EXPECT_EQ(forced, 1);
}

// A split point that rounds outside the box of the control points is brought back into it.  A curve whose control points all lie at
// one coordinate that is exactly half way between two device units (1.5) must emit that coordinate rounded half up (2) at every
// vertex however its weights make the quotient round, and one at the double just below it (1.4999999999999998) must emit 1: without
// the clamp's low side some quotient comes out one ulp low and the first rounds down, without its high side one comes out one ulp
// high and the second rounds up.
TEST(ConicRoute, ASplitPointNeverLeavesTheBoxOfTheControlPoints) {
  Random random{UINT64_C(0x2545F4914F6CDD1D)};
  int checked = 0;
  struct Constant {
    double at;
    int32_t rounds_to;
  };
  for (const Constant k : {Constant{1.5, 2}, Constant{0x1.7ffffffffffffp+0, 1}}) {
    for (int axis = 0; axis < 2; axis++) {
      for (int i = 0; i < 4000; i++) {
        Conic c;
        c.n = 3 + (i & 1);
        for (int j = 0; j < c.n; j++) {
          const double along = (double)(random.next() % 100000);  // not monotone: a control point past the end splits the piece
          (axis == 0 ? c.x : c.y)[j] = k.at;
          (axis == 0 ? c.y : c.x)[j] = along;
          c.w[j] = 1.0 + (double)(random.next() % 1000000) / 77777.0;
        }
        std::vector<Vertex> out;
        ASSERT_EQ(gshgeom_conic_route(c.n, c.x, c.y, c.w, 2, collect, &out), GSHGEOM_OK);
        for (const Vertex & v : out) {
          EXPECT_EQ(axis == 0 ? v.at.x : v.at.y, k.rounds_to) << k.at << " " << axis << " " << i;
          checked++;
        }
      }
    }
  }
  EXPECT_GT(checked, 80000);
}

// The device path of a model path is the route's lines at the tolerance the caller named, exactly: not the route at a neighbouring
// tolerance, and the same for a rational quadratic, a rational cubic and the pieces of an arc.
TEST(ConicArc, TheDevicePathIsTheRoutesLinesAtTheCallersToleranceExactly) {
  const double wq[3] = {1.0, 0.6, 1.0};
  const double wc[4] = {1.0, 2.5, 0.4, 1.0};
  for (const int32_t t : {2, 3, 4, 5, 9, 33}) {
    GSHGEOM_Path * p = model();
    ASSERT_EQ(gshgeom_path_move_to(p, P(0, 0)), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_path_rational_quadratic_to(p, P(30, 40), P(60, 0), wq), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_path_rational_cubic_to(p, P(80, -50), P(100, 70), P(120, 0), wc), GSHGEOM_OK);
    GSHGEOM_Arc_Form f = {P(60, 0), P(60, 0), P(0, 30), GSHGEOM_ANGLE_DEGREES, 0, 200, 0u};
    ASSERT_EQ(gshgeom_path_arc_conjugate_to(p, &f, nullptr), GSHGEOM_OK);
    const GSHGEOM_Map id = matrix(1, 0, 0, 1);
    GSHGEOM_Device_Path * dp = nullptr;
    ASSERT_EQ(gshgeom_device_path_from_model_with_tolerance(p, &id, t, nullptr, nullptr, &dp), GSHGEOM_OK);
    const std::vector<Dp> got = polyline_of(dp);
    // The same lines by the route itself, curve by curve.
    std::vector<Dp> want;
    want.push_back({0, 0});
    auto route = [&](int n, const GSHGEOM_Point * pts, const double * w) {
      Conic c;
      c.n = n;
      for (int k = 0; k < n; k++) {
        c.x[k] = pts[k].x * 64.0;
        c.y[k] = pts[k].y * 64.0;
        c.w[k] = w[k];
      }
      std::vector<Vertex> out;
      ASSERT_EQ(gshgeom_conic_route(n, c.x, c.y, c.w, t, collect, &out), GSHGEOM_OK);
      for (const Vertex & v : out) {
        want.push_back(v.at);
      }
    };
    const GSHGEOM_Point q[3] = {P(0, 0), P(30, 40), P(60, 0)};
    route(3, q, wq);
    const GSHGEOM_Point cu[4] = {P(60, 0), P(80, -50), P(100, 70), P(120, 0)};
    route(4, cu, wc);
    GSHGEOM_Arc a;
    ASSERT_EQ(gshgeom_path_arc(p, 3, &a), GSHGEOM_OK);
    GSHGEOM_Curve pieces[GSHGEOM_ARC_MAX_PIECES];
    int n = 0;
    ASSERT_EQ(gshgeom_arc_pieces(&a, pieces, &n), GSHGEOM_OK);
    for (int k = 0; k < n; k++) {
      route(3, pieces[k].p, pieces[k].w);
    }
    ASSERT_EQ(got.size(), want.size()) << t;
    for (size_t k = 0; k < got.size(); k++) {
      EXPECT_TRUE(got[k].x == want[k].x && got[k].y == want[k].y) << t << " " << k;
    }
    gshgeom_device_path_destroy(dp);
    gshgeom_path_destroy(p);
  }
}

int main(int argc, char ** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
