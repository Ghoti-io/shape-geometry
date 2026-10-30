/**
 * @file
 *
 * The flattener at the tolerances a caller may ask for, checked against the
 * curve.
 *
 * `src/flatten/flatten.c` states and proves one property: for every curve with
 * int32 control points and every tolerance of 2 units and up, the polyline is
 * within `tolerance - 0.12` units of the curve, and the depth cap never forces
 * a piece.  Pixels cannot test that - a hash says "different", never "by how
 * much" - so this file measures it directly: the curve is evaluated exactly at
 * a dense set of parameters inside every piece the flattener made, and the
 * distance to that piece's segment is the number asserted.
 *
 * The pieces are the real ones: the emit callback is told which parameter
 * interval each vertex's piece covers (`leaf` and `level`), so the flattener
 * is measured, not a model of it.  (test_flatten.cpp measures a model of the
 * recursion at the default tolerance and requires the model to give the
 * library's vertices; this file is the same measure on the library's own
 * pieces, at every tolerance.)  The tests here began as `font`'s, when `font`
 * had its own copy of this flattener.
 *
 * **The evaluation shares nothing with the flattener.** It is the Bernstein
 * form in `long double` at dyadic parameters, which are exact in the type, so
 * the only error is the arithmetic of the sum: below 2^-33 of a unit where
 * `long double` is the x87 type and below a millionth of a unit where it is a
 * plain `double`, against a margin of 0.12 unit.
 *
 * Why a piece-by-piece measure and not a uniform one: a curve spanning the
 * int32 range is 2^33 units long, so no affordable uniform sampling is dense
 * enough to see a deviation of a unit between two samples.  Inside a piece the
 * deviation from the chord is a smooth function of the parameter (a cubic
 * polynomial), so thirty-three samples find its maximum to within a thousandth
 * of it.
 *
 * Copyright 2026 by Corey Pennycuff
 */

#include "test_helpers.h"

#include "../alloc_tracker.h"

#include "../../src/flatten/flatten_internal.h"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

namespace {

using Dp = GSHGEOM_Device_Point;

constexpr int32_t kMin = INT32_MIN;

/** One emitted vertex and the piece of the curve it ends. */
struct Vertex {
  Dp at;
  uint32_t leaf;
  int level;
};

GSHGEOM_Result collect(void * user, Dp to, uint32_t leaf, int level) {
  static_cast<std::vector<Vertex> *>(user)->push_back({to, leaf, level});
  return GSHGEOM_OK;
}

/** A curve of either degree. */
struct Curve {
  bool cubic;
  Dp p[4];
};

/** Flatten through the public functions; the pieces the cap forced come back as -1 if refused. */
GSHGEOM_Result flatten(const Curve & curve, int32_t tolerance, std::vector<Vertex> & out) {
  out.clear();
  return curve.cubic
      ? gshgeom_flatten_cubic(tolerance, curve.p[0], curve.p[1], curve.p[2], curve.p[3], collect, &out)
      : gshgeom_flatten_quadratic(tolerance, curve.p[0], curve.p[1], curve.p[2], collect, &out);
}

/** The same, with a depth cap, reporting the pieces the cap forced. */
GSHGEOM_Result flatten_capped(const Curve & curve, int32_t tolerance, int levels,
    std::vector<Vertex> & out, int * forced) {
  out.clear();
  return curve.cubic
      ? gshgeom_flatten_cubic_capped(tolerance, levels, curve.p[0], curve.p[1], curve.p[2],
            curve.p[3], collect, &out, forced)
      : gshgeom_flatten_quadratic_capped(tolerance, levels, curve.p[0], curve.p[1], curve.p[2],
            collect, &out, forced);
}

struct Vec {
  long double x;
  long double y;
};

/** The curve at parameter `s`, from the Bernstein form. */
Vec evaluate(const Curve & curve, long double s) {
  const long double r = 1.0L - s;
  long double b[4];
  Vec out{0, 0};

  if (curve.cubic) {
    b[0] = r * r * r;
    b[1] = 3 * s * r * r;
    b[2] = 3 * s * s * r;
    b[3] = s * s * s;
  } else {
    b[0] = r * r;
    b[1] = 2 * s * r;
    b[2] = s * s;
    b[3] = 0;
  }
  for (int i = 0; i < (curve.cubic ? 4 : 3); ++i) {
    out.x += b[i] * curve.p[i].x;
    out.y += b[i] * curve.p[i].y;
  }
  return out;
}

/** Distance from `q` to the segment `a`-`b`. */
long double distance(Vec q, Dp a, Dp b) {
  const long double ax = a.x;
  const long double ay = a.y;
  const long double dx = static_cast<long double>(b.x) - ax;
  const long double dy = static_cast<long double>(b.y) - ay;
  const long double length2 = dx * dx + dy * dy;
  long double along = 0;

  if (length2 > 0) {
    along = ((q.x - ax) * dx + (q.y - ay) * dy) / length2;
    along = along < 0 ? 0 : (along > 1 ? 1 : along);
  }
  const long double ex = q.x - (ax + along * dx);
  const long double ey = q.y - (ay + along * dy);
  return std::sqrt(ex * ex + ey * ey);
}

constexpr int kSamples = 32;

/**
 * The largest distance from the curve to the polyline, taken over every piece:
 * thirty-three exact points in each, measured to that piece's segment.
 */
long double deviation(const Curve & curve, const std::vector<Vertex> & v) {
  long double worst = 0;
  Dp from = curve.p[0];

  for (const Vertex & vertex : v) {
    const long double span = std::ldexp(1.0L, vertex.level);

    for (int i = 0; i <= kSamples; ++i) {
      const long double s = (static_cast<long double>(vertex.leaf)
          + static_cast<long double>(i) / kSamples) / span;
      const long double d = distance(evaluate(curve, s), from, vertex.at);

      worst = d > worst ? d : worst;
    }
    from = vertex.at;
  }
  return worst;
}

/** xorshift64*: the same stream on every host. */
struct Random {
  uint64_t state;

  uint64_t next() {
    state ^= state >> 12;
    state ^= state << 25;
    state ^= state >> 27;
    return state * 0x2545F4914F6CDD1DULL;
  }

  /** A coordinate in `[-2^bits, 2^bits)`, or the whole int32 range for 31. */
  int32_t coordinate(int bits) {
    const uint64_t span = uint64_t{1} << (bits + 1);
    const int64_t value = static_cast<int64_t>(next() % span) - (int64_t{1} << bits);

    return static_cast<int32_t>(value);
  }
};

Curve random_curve(Random & random, bool cubic, int bits) {
  Curve curve{cubic, {}};

  for (int i = 0; i < (cubic ? 4 : 3); ++i) {
    curve.p[i] = {random.coordinate(bits), random.coordinate(bits)};
  }
  return curve;
}

/** Every control point at one of the four corners of the int32 square. */
Curve corner_curve(Random & random, bool cubic) {
  Curve curve{cubic, {}};

  for (int i = 0; i < (cubic ? 4 : 3); ++i) {
    curve.p[i] = {(random.next() & 1u) ? INT32_MAX : kMin, (random.next() & 1u) ? INT32_MAX : kMin};
  }
  return curve;
}

/** `t - 0.12`, the bound flatten.c proves, with the arithmetic's own slack. */
long double bound(int32_t tolerance) {
  return static_cast<long double>(tolerance) - 0.12L + 1e-4L;
}

/**
 * Flatten, measure, assert.  `worst` is the smallest margin seen, the tolerance
 * minus the deviation: the proof says it is at least 0.12, so it is the number
 * to watch.
 */
void check(const Curve & curve, int32_t tolerance, long double & worst, size_t & most_vertices) {
  std::vector<Vertex> vertices;
  int forced = -1;
  ASSERT_EQ(flatten_capped(curve, tolerance, GSHGEOM_FLATTEN_MAX_DEPTH, vertices, &forced), GSHGEOM_OK);
  ASSERT_FALSE(vertices.empty());
  const long double d = deviation(curve, vertices);

  EXPECT_EQ(forced, 0) << "the depth cap forced a piece, tolerance " << tolerance;
  EXPECT_LE(d, bound(tolerance)) << (curve.cubic ? "cubic" : "quadratic") << " at tolerance "
      << tolerance << " strays " << static_cast<double>(d) << " units; first point "
      << curve.p[0].x << "," << curve.p[0].y << " last " << vertices.back().at.x << ","
      << vertices.back().at.y;
  const long double margin = static_cast<long double>(tolerance) - d;

  worst = margin < worst ? margin : worst;
  most_vertices = vertices.size() > most_vertices ? vertices.size() : most_vertices;
}

const int32_t kTolerances[] = {2, 3, 4, 5, 8, 16, 64, 1000};

}  // namespace

TEST(Tolerance, TheMinimumIsWhatTheHeaderSaysAndTheDefaultIsAboveIt) {
  EXPECT_EQ(GSHGEOM_FLATTEN_MIN_TOLERANCE, 2);
  EXPECT_EQ(GSHGEOM_FLATTEN_TOLERANCE, 4);
  EXPECT_EQ(GSHGEOM_FLATTEN_MAX_DEPTH, 20);
  EXPECT_EQ(GSHGEOM_FLATTEN_REVISION, 1);
}

/**
 * The property itself, at every size a curve can have.  The classes are the
 * coordinate range `[-2^b, 2^b)`: 4 is a curve inside a quarter of a pixel's
 * cell, 16 is a few hundred pixels, 31 is the whole int32 range, which is what
 * a hostile outline can ask for and what flatten.c proves for.
 */
TEST(Tolerance, EveryPieceOfEveryCurveIsWithinTheToleranceOfTheCurve) {
  Random random{0x9E3779B97F4A7C15ULL};
  const int classes[] = {4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 31};

  for (const int bits : classes) {
    // Few big curves: a curve of the whole range is tens of thousands of
    // pieces at the smallest tolerance, each sampled thirty-three times.
    const int curves = bits <= 16 ? 24 : (bits <= 24 ? 6 : 3);
    long double worst = 1e9L;
    size_t most = 0;

    for (const int32_t tolerance : kTolerances) {
      for (int i = 0; i < curves; ++i) {
        check(random_curve(random, false, bits), tolerance, worst, most);
        check(random_curve(random, true, bits), tolerance, worst, most);
      }
    }
    // The smallest margin seen, per class, so a regression that eats the
    // margin shows in the log before it crosses the bound.
    printf("flatten: coordinates below 2^%d, smallest margin %.3f units, most vertices %zu\n",
        bits, static_cast<double>(worst), most);
  }
}

/** The curves that make the second differences as large as an int32 pair allows. */
TEST(Tolerance, TheExtremeCurvesOfTheInt32SquareAreWithinTheTolerance) {
  Random random{0x0123456789ABCDEFULL};
  long double worst = 1e9L;
  size_t most = 0;

  for (const int32_t tolerance : {2, 4, 64}) {
    for (int i = 0; i < 4; ++i) {
      check(corner_curve(random, false), tolerance, worst, most);
      check(corner_curve(random, true), tolerance, worst, most);
    }
  }
  printf("flatten: int32 corners, smallest margin %.3f units, most vertices %zu\n",
      static_cast<double>(worst), most);
}

/** Small tolerances are where a constant threshold would fail; see flatten.c. */
TEST(Tolerance, EveryToleranceFromTheMinimumIsHonouredOnTheSameCurves) {
  Random random{0xDEADBEEFCAFEF00DULL};

  for (int32_t tolerance = GSHGEOM_FLATTEN_MIN_TOLERANCE; tolerance <= 40; ++tolerance) {
    long double worst = 1e9L;
    size_t most = 0;

    for (int i = 0; i < 6; ++i) {
      check(random_curve(random, false, 16), tolerance, worst, most);
      check(random_curve(random, true, 16), tolerance, worst, most);
    }
  }
}

/** The largest tolerance there is still gets a polyline, a coarse one, and the bound. */
TEST(Tolerance, TheLargestToleranceIsAcceptedAndStillWithinTheBound) {
  Random random{0x0F0F0F0F12345678ULL};
  long double worst = 1e9L;
  size_t most = 0;

  for (int i = 0; i < 4; ++i) {
    check(random_curve(random, false, 31), INT32_MAX, worst, most);
    check(random_curve(random, true, 31), INT32_MAX, worst, most);
  }
}

/**
 * Answers worked out by arithmetic independent of this library: a Python model
 * of the flattener using unbounded integers (`tools/gen-hash-vectors.py`, and
 * `font`'s own model when it had this flattener), recorded as a vertex count,
 * the first and last vertex and a hash of all of them.  They pin the rounding
 * rules - ties upward at both the midpoints and the narrowing, floors of
 * negative values - which the distance test is too coarse to notice (a floor
 * midpoint stays within the bound on most curves), at tolerances other than
 * the default, where the threshold's derivation from the tolerance is what is
 * being pinned.
 */
TEST(Tolerance, TheVertexListsAreWhatTheIndependentModelSays) {
  struct Known {
    bool cubic;
    int32_t tolerance;
    Dp p[4];
    size_t count;
    uint64_t hash;
    Dp first;
    Dp last;
  };
  const Known known[] = {
      {false, 4, {{0, 0}, {3001, 4999}, {8003, 1}, {0, 0}}, 32, 0x5e9fe3e61234e30bULL, {190, 303}, {8003, 1}},
      {false, 4, {{-1000, -777}, {-3001, 5001}, {777, -9999}, {0, 0}}, 64, 0x4bbe0f147b91fd7dULL, {-1061, -602}, {777, -9999}},
      {false, 2, {{0, 0}, {400, -1200}, {800, 0}, {0, 0}}, 32, 0xc71abb9947fcae38ULL, {25, -73}, {800, 0}},
      {true, 4, {{0, 0}, {1000, 3001}, {-2001, 4001}, {5000, -7}}, 59, 0xdd37fa8b0432a31fULL, {44, 139}, {5000, -7}},
      {true, 3, {{-31, 17}, {-8003, 2999}, {9001, -6001}, {64, -64}}, 109, 0x29ecf687a1bbf682ULL, {-213, 85}, {64, -64}},
      {false, 64, {{INT32_MIN, INT32_MAX}, {INT32_MAX, INT32_MIN}, {INT32_MIN, INT32_MIN}, {0, 0}}, 8192, 0x1d73f31c84978d9dULL, {-2146435200, 2146435135}, {INT32_MIN, INT32_MIN}},
      {true, 2, {{INT32_MIN, INT32_MIN}, {INT32_MAX, INT32_MIN}, {INT32_MIN, INT32_MAX}, {INT32_MAX, INT32_MAX}}, 97787, 0xea23931ffec1d46bULL, {-2147385345, -INT32_MAX}, {INT32_MAX, INT32_MAX}},
  };

  for (const Known & k : known) {
    Curve curve{k.cubic, {k.p[0], k.p[1], k.p[2], k.p[3]}};
    std::vector<Vertex> vertices;
    uint64_t hash = 0xCBF29CE484222325ULL;

    ASSERT_EQ(flatten(curve, k.tolerance, vertices), GSHGEOM_OK);
    ASSERT_EQ(vertices.size(), k.count);
    for (const Vertex & vertex : vertices) {
      for (const int32_t value : {vertex.at.x, vertex.at.y}) {
        const uint32_t bits = static_cast<uint32_t>(value);

        for (int shift = 0; shift < 32; shift += 8) {
          hash = (hash ^ ((bits >> shift) & 0xFFu)) * 0x100000001B3ULL;
        }
      }
    }
    EXPECT_EQ(hash, k.hash);
    EXPECT_EQ(vertices.front().at.x, k.first.x);
    EXPECT_EQ(vertices.front().at.y, k.first.y);
    EXPECT_EQ(vertices.back().at.x, k.last.x);
    EXPECT_EQ(vertices.back().at.y, k.last.y);
  }
}

/**
 * The depth cap, which the proof says no input reaches, is exercised with a
 * smaller one: a piece it forces is counted, still emitted, and the pieces
 * still tile the curve.  Without this the count is a number nothing has ever
 * seen move.
 */
TEST(Tolerance, APieceTheCapForcesIsCountedAndStillTilesTheCurve) {
  const Curve quad{false, {{0, 0}, {1 << 20, 1 << 22}, {1 << 21, 0}, {0, 0}}};
  const Curve cubic{true, {{0, 0}, {1 << 20, 1 << 22}, {-(1 << 21), 1 << 21}, {1 << 21, 0}}};

  for (const Curve * curve : {&quad, &cubic}) {
    for (const int levels : {1, 3, 6}) {
      std::vector<Vertex> vertices;
      int forced = -1;
      uint64_t next = 0;

      ASSERT_EQ(flatten_capped(*curve, 4, levels, vertices, &forced), GSHGEOM_OK);
      EXPECT_GT(forced, 0) << "levels " << levels;
      EXPECT_EQ(static_cast<size_t>(forced), vertices.size())
          << "every piece of a capped curve this large is forced";
      for (const Vertex & vertex : vertices) {
        EXPECT_EQ(vertex.level, levels);
        EXPECT_EQ(static_cast<uint64_t>(vertex.leaf), next);
        next += 1;
        // The vertex a forced piece emits is where the curve is at the end of
        // the piece, to the one rounding of a vertex (at most 0.71 units) and
        // the noise of the halvings (a few thousandths): not the piece's start
        // or its control point, which a flattener that ended a piece early
        // with the wrong point would give.
        const long double s = static_cast<long double>(vertex.leaf + 1)
            / std::ldexp(1.0L, vertex.level);
        const Vec at = evaluate(*curve, s);
        const long double dx = at.x - vertex.at.x;
        const long double dy = at.y - vertex.at.y;

        EXPECT_LE(std::sqrt(dx * dx + dy * dy), 1.0L) << "levels " << levels << " piece " << vertex.leaf;
      }
      EXPECT_EQ(next, uint64_t{1} << levels);
    }
  }
  // The same cap applied to a curve it does not bind counts nothing.
  std::vector<Vertex> vertices;
  int forced = -1;
  const Curve flat{false, {{0, 0}, {320, 160}, {640, 320}, {0, 0}}};

  ASSERT_EQ(flatten_capped(flat, 4, 1, vertices, &forced), GSHGEOM_OK);
  EXPECT_EQ(forced, 0);
}

/** A cap outside 1..GSHGEOM_FLATTEN_MAX_DEPTH is brought inside it, never past it. */
TEST(Tolerance, TheCapIsClampedToWhatTheProofCovers) {
  const Curve curve{true, {{0, 0}, {1 << 20, 1 << 22}, {-(1 << 21), 1 << 21}, {1 << 21, 0}}};
  std::vector<Vertex> low, one, high, real;
  int forced = -1;

  ASSERT_EQ(flatten_capped(curve, 4, -3, low, &forced), GSHGEOM_OK);
  ASSERT_EQ(flatten_capped(curve, 4, 1, one, &forced), GSHGEOM_OK);
  ASSERT_EQ(low.size(), one.size());
  ASSERT_EQ(flatten_capped(curve, 4, 1000, high, &forced), GSHGEOM_OK);
  EXPECT_EQ(forced, 0);
  ASSERT_EQ(flatten(curve, 4, real), GSHGEOM_OK);
  ASSERT_EQ(high.size(), real.size());
  for (size_t i = 0; i < high.size(); ++i) {
    EXPECT_EQ(high[i].at.x, real[i].at.x);
    EXPECT_EQ(high[i].at.y, real[i].at.y);
    EXPECT_LE(high[i].level, GSHGEOM_FLATTEN_MAX_DEPTH);
  }
}

/** The pieces tile the parameter interval, in order, with no gap or overlap. */
TEST(Tolerance, ThePiecesTileTheCurveInOrderAndEndAtTheLastPoint) {
  Random random{0x5555AAAA5555AAAAULL};

  for (const int bits : {6, 14, 22, 31}) {
    for (const bool cubic : {false, true}) {
      const Curve curve = random_curve(random, cubic, bits);
      std::vector<Vertex> vertices;
      uint64_t next = 0;

      ASSERT_EQ(flatten(curve, 4, vertices), GSHGEOM_OK);
      ASSERT_FALSE(vertices.empty());
      for (const Vertex & vertex : vertices) {
        ASSERT_LE(vertex.level, GSHGEOM_FLATTEN_MAX_DEPTH);
        const uint64_t start = static_cast<uint64_t>(vertex.leaf)
            << (GSHGEOM_FLATTEN_MAX_DEPTH - vertex.level);
        const uint64_t width = uint64_t{1} << (GSHGEOM_FLATTEN_MAX_DEPTH - vertex.level);

        ASSERT_EQ(start, next);
        next = start + width;
      }
      EXPECT_EQ(next, uint64_t{1} << GSHGEOM_FLATTEN_MAX_DEPTH);
      const Dp last = curve.p[cubic ? 3 : 2];

      EXPECT_EQ(vertices.back().at.x, last.x);
      EXPECT_EQ(vertices.back().at.y, last.y);
    }
  }
}

TEST(Tolerance, AFlatCurveIsOnePieceAndADegenerateOneIsToo) {
  std::vector<Vertex> vertices;
  const Curve straight{false, {{0, 0}, {320, 160}, {640, 320}, {0, 0}}};
  const Curve point{true, {{-7, 9}, {-7, 9}, {-7, 9}, {-7, 9}}};

  ASSERT_EQ(flatten(straight, 4, vertices), GSHGEOM_OK);
  EXPECT_EQ(vertices.size(), 1u);
  ASSERT_EQ(flatten(point, 4, vertices), GSHGEOM_OK);
  ASSERT_EQ(vertices.size(), 1u);
  EXPECT_EQ(vertices[0].at.x, -7);
  EXPECT_EQ(vertices[0].at.y, 9);
}

/**
 * A tolerance below the minimum is refused, not raised to it: the caller asked
 * for a promise the proof cannot keep, and nothing is emitted.
 */
TEST(Tolerance, ATooSmallToleranceIsRefusedAndEmitsNothing) {
  const Curve curve{true, {{0, 0}, {900, 3000}, {2000, -100}, {64, 64}}};
  const Curve quad{false, {{0, 0}, {900, 3000}, {2000, -100}, {0, 0}}};

  for (const int32_t tolerance : {INT32_MIN, -5, -1, 0, 1}) {
    std::vector<Vertex> vertices;
    int forced = 7;

    EXPECT_EQ(flatten(curve, tolerance, vertices), GSHGEOM_ERR_INVALID) << tolerance;
    EXPECT_TRUE(vertices.empty());
    EXPECT_EQ(flatten(quad, tolerance, vertices), GSHGEOM_ERR_INVALID) << tolerance;
    EXPECT_TRUE(vertices.empty());
    EXPECT_EQ(flatten_capped(quad, tolerance, 3, vertices, &forced), GSHGEOM_ERR_INVALID);
    EXPECT_EQ(forced, 0) << "a refused call reports no forced pieces";
  }
  std::vector<Vertex> vertices;

  EXPECT_EQ(flatten(curve, GSHGEOM_FLATTEN_MIN_TOLERANCE, vertices), GSHGEOM_OK);
  EXPECT_FALSE(vertices.empty());
}

TEST(Tolerance, ANullSinkIsRefused) {
  const Dp a{0, 0}, b{100, 300}, c{200, 0}, d{300, 300};

  EXPECT_EQ(gshgeom_flatten_quadratic(4, a, b, c, nullptr, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_flatten_cubic(4, a, b, c, d, nullptr, nullptr), GSHGEOM_ERR_INVALID);
}

/**
 * A sink that cannot take another vertex says so and the flattening stops
 * there, with that result, instead of running the curve to its end.
 */
TEST(Tolerance, ASinkThatRefusesStopsTheFlatteningAndItsResultIsReturned) {
  struct Stop {
    size_t taken = 0;
    size_t room;
  };
  const auto sink = [](void * user, Dp, uint32_t, int) -> GSHGEOM_Result {
    Stop * stop = static_cast<Stop *>(user);

    if (stop->taken >= stop->room) {
      return GSHGEOM_ERR_LIMIT;
    }
    stop->taken += 1;
    return GSHGEOM_OK;
  };
  const Dp a{0, 0}, b{4000, 9000}, c{-3000, 9000}, d{8000, 0};
  std::vector<Vertex> all;

  ASSERT_EQ(gshgeom_flatten_cubic(4, a, b, c, d, collect, &all), GSHGEOM_OK);
  ASSERT_GT(all.size(), 10u);
  for (const size_t room : {size_t{0}, size_t{1}, all.size() / 2, all.size() - 1}) {
    Stop stop{0, room};

    EXPECT_EQ(gshgeom_flatten_cubic(4, a, b, c, d, sink, &stop), GSHGEOM_ERR_LIMIT) << room;
    EXPECT_EQ(stop.taken, room) << "nothing was emitted after the refusal";
  }
  Stop exact{0, all.size()};

  EXPECT_EQ(gshgeom_flatten_cubic(4, a, b, c, d, sink, &exact), GSHGEOM_OK);
  EXPECT_EQ(exact.taken, all.size());
  for (const size_t room : {size_t{0}, size_t{3}}) {
    Stop stop{0, room};

    EXPECT_EQ(gshgeom_flatten_quadratic(4, a, b, d, sink, &stop), GSHGEOM_ERR_LIMIT) << room;
    EXPECT_EQ(stop.taken, room);
  }
}

// ---------------------------------------------------------------------------
// The path functions
// ---------------------------------------------------------------------------

namespace {

GSHGEOM_Device_Path * path_of(const Curve & curve) {
  GSHGEOM_Device_Path * path = nullptr;

  EXPECT_EQ(gshgeom_device_path_create(nullptr, nullptr, &path), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_device_path_move_to(path, curve.p[0]), GSHGEOM_OK);
  if (curve.cubic) {
    EXPECT_EQ(gshgeom_device_path_cubic_to(path, curve.p[1], curve.p[2], curve.p[3]), GSHGEOM_OK);
  } else {
    EXPECT_EQ(gshgeom_device_path_quadratic_to(path, curve.p[1], curve.p[2]), GSHGEOM_OK);
  }
  return path;
}

/** The flattened path's points after the move, as vertices of the curve. */
std::vector<Dp> points_of(const GSHGEOM_Flat_Path * flat) {
  const Dp * points = gshgeom_flat_path_points(flat);
  const size_t count = gshgeom_flat_path_point_count(flat);

  return std::vector<Dp>(points, points + count);
}

}  // namespace

TEST(TolerancePath, ThePathAtTheDefaultToleranceIsTheOneWithoutAToleranceArgument) {
  Random random{0x13579BDF02468ACEULL};

  for (const bool cubic : {false, true}) {
    for (const int bits : {8, 16, 31}) {
      const Curve curve = random_curve(random, cubic, bits);
      GSHGEOM_Device_Path * path = path_of(curve);
      GSHGEOM_Flat_Path * plain = nullptr;
      GSHGEOM_Flat_Path * with = nullptr;

      ASSERT_EQ(gshgeom_device_path_flatten(path, nullptr, nullptr, &plain), GSHGEOM_OK);
      ASSERT_EQ(gshgeom_device_path_flatten_with_tolerance(
                    path, GSHGEOM_FLATTEN_TOLERANCE, nullptr, nullptr, &with), GSHGEOM_OK);
      EXPECT_EQ(gshgeom_flat_path_hash(plain), gshgeom_flat_path_hash(with));
      gshgeom_flat_path_destroy(plain);
      gshgeom_flat_path_destroy(with);
      gshgeom_device_path_destroy(path);
    }
  }
}

/** The path function is the curve function with the start point in front. */
TEST(TolerancePath, ThePathAtEveryToleranceIsTheCurveFunctionsVertices) {
  Random random{0x2468ACE013579BDFULL};

  for (const int32_t tolerance : {2, 3, 4, 7, 64, 1000, INT32_MAX}) {
    for (const bool cubic : {false, true}) {
      const Curve curve = random_curve(random, cubic, 20);
      GSHGEOM_Device_Path * path = path_of(curve);
      GSHGEOM_Flat_Path * flat = nullptr;
      std::vector<Vertex> vertices;

      ASSERT_EQ(flatten(curve, tolerance, vertices), GSHGEOM_OK);
      ASSERT_EQ(gshgeom_device_path_flatten_with_tolerance(path, tolerance, nullptr, nullptr, &flat),
          GSHGEOM_OK);
      const std::vector<Dp> points = points_of(flat);

      ASSERT_EQ(points.size(), vertices.size() + 1) << "the move, then one point per vertex";
      EXPECT_EQ(points[0].x, curve.p[0].x);
      EXPECT_EQ(points[0].y, curve.p[0].y);
      for (size_t i = 0; i < vertices.size(); ++i) {
        EXPECT_EQ(points[i + 1].x, vertices[i].at.x) << "tolerance " << tolerance << " vertex " << i;
        EXPECT_EQ(points[i + 1].y, vertices[i].at.y);
      }
      gshgeom_flat_path_destroy(flat);
      gshgeom_device_path_destroy(path);
    }
  }
}

TEST(TolerancePath, ALargerToleranceGivesAShorterPolyline) {
  const Curve curve{true, {{0, 0}, {6400, 12800}, {-3200, 12800}, {9600, 0}}};
  GSHGEOM_Device_Path * path = path_of(curve);
  size_t before = SIZE_MAX;

  for (const int32_t tolerance : {2, 4, 16, 128, 4096}) {
    GSHGEOM_Flat_Path * flat = nullptr;

    ASSERT_EQ(gshgeom_device_path_flatten_with_tolerance(path, tolerance, nullptr, nullptr, &flat),
        GSHGEOM_OK);
    EXPECT_LT(gshgeom_flat_path_point_count(flat), before) << tolerance;
    before = gshgeom_flat_path_point_count(flat);
    gshgeom_flat_path_destroy(flat);
  }
  gshgeom_device_path_destroy(path);
}

TEST(TolerancePath, ATooSmallToleranceIsRefusedBeforeAnythingIsAllocatedAndNullsAreRefused) {
  const Curve curve{false, {{0, 0}, {100, 300}, {200, 0}, {0, 0}}};
  GSHGEOM_Device_Path * path = path_of(curve);
  tracker::Tracker count;
  GSHGEOM_Allocator allocator = tracker::allocator_for(&count);

  for (const int32_t tolerance : {INT32_MIN, -1, 0, 1}) {
    GSHGEOM_Flat_Path * flat = reinterpret_cast<GSHGEOM_Flat_Path *>(0x1);

    EXPECT_EQ(gshgeom_device_path_flatten_with_tolerance(path, tolerance, &allocator, nullptr, &flat),
        GSHGEOM_ERR_INVALID) << tolerance;
    EXPECT_EQ(flat, reinterpret_cast<GSHGEOM_Flat_Path *>(0x1)) << "the output is written only on success";
  }
  EXPECT_EQ(count.calls, 0) << "a refused tolerance allocates nothing";
  GSHGEOM_Flat_Path * flat = nullptr;

  EXPECT_EQ(gshgeom_device_path_flatten_with_tolerance(nullptr, 4, nullptr, nullptr, &flat),
      GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_device_path_flatten_with_tolerance(path, 4, nullptr, nullptr, nullptr),
      GSHGEOM_ERR_INVALID);
  gshgeom_device_path_destroy(path);
}

/** The limit on points applies at any tolerance, and the result is the limit's. */
TEST(TolerancePath, TheLimitOnPointsBindsAtEveryTolerance) {
  const Curve curve{true, {{0, 0}, {64000, 128000}, {-32000, 128000}, {96000, 0}}};
  GSHGEOM_Device_Path * path = path_of(curve);

  for (const int32_t tolerance : {2, 4, 64}) {
    GSHGEOM_Flat_Path * flat = nullptr;
    GSHGEOM_Limits limits;

    ASSERT_EQ(gshgeom_device_path_flatten_with_tolerance(path, tolerance, nullptr, nullptr, &flat),
        GSHGEOM_OK);
    const size_t needed = gshgeom_flat_path_point_count(flat);

    gshgeom_flat_path_destroy(flat);
    gshgeom_limits_default(&limits);
    limits.max_points = needed;
    ASSERT_EQ(gshgeom_device_path_flatten_with_tolerance(path, tolerance, nullptr, &limits, &flat),
        GSHGEOM_OK);
    gshgeom_flat_path_destroy(flat);
    limits.max_points = needed - 1;
    flat = nullptr;
    EXPECT_EQ(gshgeom_device_path_flatten_with_tolerance(path, tolerance, nullptr, &limits, &flat),
        GSHGEOM_ERR_LIMIT) << tolerance;
    EXPECT_EQ(flat, nullptr);
  }
  gshgeom_device_path_destroy(path);
}

int main(int argc, char ** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
