/**
 * @file
 *
 * The 26.6 rounding rule and the device path (CAP-3, CAP-13): the rule against
 * an independent exact reference at every half-way point and at the range
 * limit, the builder's rules and limits, a model path under a map, the hash and
 * the dump pinned by known answers, and allocation failure at every call.
 *
 * Copyright 2026 by Corey Pennycuff
 */

#include "test_helpers.h"

#include "../alloc_tracker.h"
#include "../known_answers.h"
#include "../dump_sink.h"

#include <algorithm>
#include <cfloat>
#include <clocale>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <functional>
#include <limits>
#include <set>
#include <string>
#include <vector>

namespace {

const double kNaN = std::numeric_limits<double>::quiet_NaN();
const double kInf = std::numeric_limits<double>::infinity();
const int32_t kSentinel = 123456789;

using Dp = GSHGEOM_Device_Point;

Dp D(int32_t x, int32_t y) {
  return Dp{x, y};
}

GSHGEOM_Point P(double x, double y) {
  return GSHGEOM_Point{x, y, 0.0};
}

// A second, independent statement of the rule: exact integer arithmetic on the
// mantissa and exponent of x (no division, no floor of a double).  Returns false
// when the result does not fit an int32.
bool reference_round(double x, int32_t * out) {
  if (x == 0.0) {
    *out = 0;
    return true;
  }
  int e;
  double frac = std::frexp(std::fabs(x), &e);  // |x| = frac * 2^e, frac in [0.5, 1)
  int64_t m = (int64_t)std::ldexp(frac, 53);   // exact, 2^52 <= m < 2^53: |x| = m * 2^(e - 53)
  int shift = (e - 53) + 6;                    // |x| * 64 = m * 2^shift
  if (shift >= 0) {
    return false;  // |x * 64| >= 2^52: far outside an int32
  }
  int s = -shift;
  if (s >= 55) {
    *out = 0;  // |x * 64| < 1/4: floor(v + 1/2) is 0 on both sides of zero
    return true;
  }
  int64_t n = x < 0.0 ? -m : m;
  int64_t d = (int64_t)1 << s;
  // floor((2 n + d) / (2 d)) = floor(n / d + 1/2), with an arithmetic shift as
  // the floor of a signed quotient by a power of two.
  int64_t value = (2 * n + d) >> (s + 1);
  if (value < INT32_MIN || value > INT32_MAX) {
    return false;
  }
  *out = (int32_t)value;
  return true;
}

std::string dump_of(const GSHGEOM_Device_Path * p) {
  FILE * f = std::tmpfile();
  EXPECT_NE(f, nullptr);
  EXPECT_EQ(gshgeom_device_path_dump(p, f), GSHGEOM_OK);
  std::rewind(f);
  std::string s;
  char buf[256];
  size_t n;
  while ((n = std::fread(buf, 1, sizeof buf, f)) > 0) {
    s.append(buf, n);
  }
  std::fclose(f);
  return s;
}

struct Snapshot {
  size_t count;
  uint64_t hash;
  std::string dump;
  explicit Snapshot(const GSHGEOM_Device_Path * p)
      : count(gshgeom_device_path_element_count(p)),
        hash(gshgeom_device_path_hash(p)),
        dump(dump_of(p)) {}
  bool operator==(const Snapshot & o) const {
    return count == o.count && hash == o.hash && dump == o.dump;
  }
};

GSHGEOM_Device_Path * make_dp(const GSHGEOM_Allocator * a = nullptr,
    const GSHGEOM_Limits * l = nullptr) {
  GSHGEOM_Device_Path * p = nullptr;
  EXPECT_EQ(gshgeom_device_path_create(a, l, &p), GSHGEOM_OK);
  return p;
}

// The sequence of tools/gen-hash-vectors.py's DEV_FIXTURE.
using DStep = std::function<GSHGEOM_Result(GSHGEOM_Device_Path *)>;

std::vector<DStep> device_steps() {
  return {
      [](GSHGEOM_Device_Path * p) { return gshgeom_device_path_move_to(p, D(64, -128)); },
      [](GSHGEOM_Device_Path * p) { return gshgeom_device_path_line_to(p, D(192, 256)); },
      [](GSHGEOM_Device_Path * p) {
        return gshgeom_device_path_quadratic_to(p, D(320, 416), D(448, 512));
      },
      [](GSHGEOM_Device_Path * p) {
        return gshgeom_device_path_cubic_to(p, D(576, 640), D(720, 768), D(832, -896));
      },
      [](GSHGEOM_Device_Path * p) { return gshgeom_device_path_close(p); },
      [](GSHGEOM_Device_Path * p) {
        return gshgeom_device_path_line_to(p, D(INT32_MIN, INT32_MAX));
      },
      [](GSHGEOM_Device_Path * p) { return gshgeom_device_path_close(p); },
  };
}

GSHGEOM_Device_Path * fixture_dp(const GSHGEOM_Allocator * a = nullptr) {
  GSHGEOM_Device_Path * p = make_dp(a);
  for (const DStep & s : device_steps()) {
    EXPECT_EQ(s(p), GSHGEOM_OK);
  }
  return p;
}

GSHGEOM_Map matrix_map(const double m[4], const double t[2]) {
  GSHGEOM_Map map;
  EXPECT_EQ(gshgeom_map_from_matrix(2, m, t, &map), GSHGEOM_OK);
  return map;
}

GSHGEOM_Map identity2() {
  GSHGEOM_Map map;
  EXPECT_EQ(gshgeom_map_identity(2, &map), GSHGEOM_OK);
  return map;
}

GSHGEOM_Path * make_model(int dim = 2, const GSHGEOM_Allocator * a = nullptr,
    const GSHGEOM_Limits * l = nullptr) {
  GSHGEOM_Path * p = nullptr;
  EXPECT_EQ(gshgeom_path_create(a, l, dim, &p), GSHGEOM_OK);
  return p;
}

using MStep = std::function<GSHGEOM_Result(GSHGEOM_Path *)>;

// The model path of tools/gen-hash-vectors.py's MODEL_DEV_PATH.
std::vector<MStep> model_steps() {
  return {
      [](GSHGEOM_Path * p) { return gshgeom_path_move_to(p, P(0.0078125, -0.0078125)); },
      [](GSHGEOM_Path * p) { return gshgeom_path_line_to(p, P(1.5, -2.25)); },
      [](GSHGEOM_Path * p) { return gshgeom_path_quadratic_to(p, P(3.0, 4.5), P(5.0, 1.0)); },
      [](GSHGEOM_Path * p) {
        return gshgeom_path_cubic_to(p, P(6.0, 7.0), P(8.25, 9.0), P(10.0, 11.5));
      },
      [](GSHGEOM_Path * p) { return gshgeom_path_close(p); },
      [](GSHGEOM_Path * p) { return gshgeom_path_line_to(p, P(2.5, 3.5)); },
  };
}

GSHGEOM_Path * model_fixture() {
  GSHGEOM_Path * p = make_model();
  for (const MStep & s : model_steps()) {
    EXPECT_EQ(s(p), GSHGEOM_OK);
  }
  return p;
}

}  // namespace

// ---------------------------------------------------------------------------
// The rounding rule
// ---------------------------------------------------------------------------

TEST(Round, KnownAnswersFromTheIndependentReference) {
  size_t n = sizeof kRoundCases / sizeof kRoundCases[0];
  ASSERT_GT(n, 900u);
  size_t limits = 0;
  for (size_t i = 0; i < n; i++) {
    const RoundCase & c = kRoundCases[i];
    int32_t v = kSentinel;
    GSHGEOM_Result r = gshgeom_round_26_6(c.x, &v);
    if (c.limit) {
      limits++;
      EXPECT_EQ(r, GSHGEOM_ERR_LIMIT) << std::hexfloat << c.x;
      EXPECT_EQ(v, kSentinel) << "written only on success: " << std::hexfloat << c.x;
    } else {
      ASSERT_EQ(r, GSHGEOM_OK) << std::hexfloat << c.x;
      EXPECT_EQ(v, c.value) << std::hexfloat << c.x;
    }
  }
  EXPECT_GT(limits, 20u);  // the vectors do reach the limit
}

TEST(Round, TheMathematicalRuleIsNotTheDoubleExpression) {
  // floor(t + 0.5) in double rounds 0.5 - 2^-54 + 0.5 up to 1; the rule says 0.
  double x = kNaiveRoundWrongInput;
  EXPECT_EQ(std::floor(x * 64.0 + 0.5), 1.0) << "the double expression is wrong here";
  int32_t v = kSentinel;
  ASSERT_EQ(gshgeom_round_26_6(x, &v), GSHGEOM_OK);
  EXPECT_EQ(v, 0);
  ASSERT_EQ(gshgeom_round_26_6(-x, &v), GSHGEOM_OK);
  EXPECT_EQ(v, 0);
  // The vectors contain the inputs on which the double expression is wrong, and
  // the library is right on every one of them.
  unsigned wrong = 0;
  for (const RoundCase & c : kRoundCases) {
    if (!c.limit && (int32_t)std::floor(c.x * 64.0 + 0.5) != c.value) {
      wrong++;
    }
  }
  EXPECT_EQ(wrong, kNaiveRoundWrong);
  EXPECT_GE(wrong, 2u);
}

TEST(Round, HalfWayPointsGoUp) {
  for (int k = -300; k <= 300; k++) {
    int32_t v = kSentinel;
    double half = (k + 0.5) / 64.0;
    ASSERT_EQ(gshgeom_round_26_6(half, &v), GSHGEOM_OK);
    EXPECT_EQ(v, k + 1) << "k=" << k;
    ASSERT_EQ(gshgeom_round_26_6(std::nextafter(half, -kInf), &v), GSHGEOM_OK);
    EXPECT_EQ(v, k) << "just below, k=" << k;
    ASSERT_EQ(gshgeom_round_26_6(std::nextafter(half, kInf), &v), GSHGEOM_OK);
    EXPECT_EQ(v, k + 1) << "just above, k=" << k;
    ASSERT_EQ(gshgeom_round_26_6(k / 64.0, &v), GSHGEOM_OK);
    EXPECT_EQ(v, k) << "an integer is itself, k=" << k;
  }
  int32_t v;
  ASSERT_EQ(gshgeom_round_26_6(-0.5 / 64.0, &v), GSHGEOM_OK);
  EXPECT_EQ(v, 0) << "-0.5 goes up to 0, not away from zero to -1";
  ASSERT_EQ(gshgeom_round_26_6(-1.5 / 64.0, &v), GSHGEOM_OK);
  EXPECT_EQ(v, -1);
  ASSERT_EQ(gshgeom_round_26_6(-0.0, &v), GSHGEOM_OK);
  EXPECT_EQ(v, 0);
  ASSERT_EQ(gshgeom_round_26_6(std::numeric_limits<double>::denorm_min(), &v), GSHGEOM_OK);
  EXPECT_EQ(v, 0);
  ASSERT_EQ(gshgeom_round_26_6(-std::numeric_limits<double>::denorm_min(), &v), GSHGEOM_OK);
  EXPECT_EQ(v, 0);
}

TEST(Round, AgreesWithAnExactIntegerReferenceOnASweep) {
  uint64_t s = UINT64_C(0x9E3779B97F4A7C15);
  auto next = [&s]() {
    s ^= s >> 12;
    s ^= s << 25;
    s ^= s >> 27;
    return s * UINT64_C(2685821657736338717);
  };
  unsigned ok = 0, limit = 0;
  for (int i = 0; i < 300000; i++) {
    double x;
    uint64_t r = next();
    switch (i % 3) {
    case 0: {  // random magnitude from 2^-30 to 2^28
      int e = (int)(r % 59) - 30;
      x = std::ldexp((double)((next() >> 11) | (UINT64_C(1) << 52)), e - 52);
      break;
    }
    case 1: {  // within a few ulps of a half-way point
      int64_t k = (int64_t)(next() >> 33) - (int64_t(1) << 30);
      x = ((double)k + 0.5) / 64.0;
      for (unsigned j = r % 5; j > 0; j--) {
        x = std::nextafter(x, (next() & 1) ? kInf : -kInf);
      }
      break;
    }
    default: {  // any bit pattern with a modest exponent
      int e = (int)(r % 80) - 50;
      x = std::ldexp((double)(next() >> 11), e);
      break;
    }
    }
    if (next() & 1) {
      x = -x;
    }
    int32_t want = 0, got = kSentinel;
    bool fits = reference_round(x, &want);
    GSHGEOM_Result res = gshgeom_round_26_6(x, &got);
    if (fits) {
      ASSERT_EQ(res, GSHGEOM_OK) << std::hexfloat << x;
      ASSERT_EQ(got, want) << std::hexfloat << x;
      ok++;
    } else {
      ASSERT_EQ(res, GSHGEOM_ERR_LIMIT) << std::hexfloat << x;
      ASSERT_EQ(got, kSentinel);
      limit++;
    }
  }
  EXPECT_GT(ok, 100000u);
  EXPECT_GT(limit, 1000u);
}

TEST(Round, NearTheRangeLimitNeverWraps) {
  struct Case {
    double x;
    bool fits;
    int32_t value;
  };
  const double hi = 2147483647.0;   // t values
  const double lo = -2147483648.0;
  const Case cases[] = {
      {hi / 64.0, true, INT32_MAX},
      {std::nextafter(hi / 64.0, kInf), true, INT32_MAX},
      {std::nextafter((hi + 0.5) / 64.0, 0.0), true, INT32_MAX},
      {(hi + 0.5) / 64.0, false, 0},
      {std::nextafter((hi + 0.5) / 64.0, kInf), false, 0},
      {(hi + 1.0) / 64.0, false, 0},
      {lo / 64.0, true, INT32_MIN},
      {std::nextafter(lo / 64.0, kInf), true, INT32_MIN},
      {(lo - 0.5) / 64.0, true, INT32_MIN},  // half up: -2147483648.5 + 0.5
      {std::nextafter((lo - 0.5) / 64.0, -kInf), false, 0},
      {std::nextafter((lo - 0.5) / 64.0, 0.0), true, INT32_MIN},
      {(lo - 1.0) / 64.0, false, 0},
      // Where a wrapped 32-bit value would look plausible.
      {4294967296.0 / 64.0, false, 0},
      {-4294967296.0 / 64.0, false, 0},
      {(4294967296.0 + 64.0) / 64.0, false, 0},
      {1e18, false, 0},
      {-1e18, false, 0},
      {DBL_MAX, false, 0},
      {-DBL_MAX, false, 0},
  };
  for (const Case & c : cases) {
    int32_t v = kSentinel;
    GSHGEOM_Result r = gshgeom_round_26_6(c.x, &v);
    if (c.fits) {
      ASSERT_EQ(r, GSHGEOM_OK) << std::hexfloat << c.x;
      EXPECT_EQ(v, c.value) << std::hexfloat << c.x;
    } else {
      EXPECT_EQ(r, GSHGEOM_ERR_LIMIT) << std::hexfloat << c.x;
      EXPECT_EQ(v, kSentinel) << "the output is untouched: " << std::hexfloat << c.x;
    }
  }
}

TEST(Round, RefusesWhatIsNotANumber) {
  int32_t v = kSentinel;
  EXPECT_EQ(gshgeom_round_26_6(kNaN, &v), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_round_26_6(kInf, &v), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_round_26_6(-kInf, &v), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(v, kSentinel);
  EXPECT_EQ(gshgeom_round_26_6(1.0, nullptr), GSHGEOM_ERR_INVALID);
}

// ---------------------------------------------------------------------------
// The device path
// ---------------------------------------------------------------------------

TEST(DeviceVerb, TheValuesAreFixed) {
  EXPECT_EQ(GSHGEOM_DEVICE_MOVE, 0);
  EXPECT_EQ(GSHGEOM_DEVICE_LINE, 1);
  EXPECT_EQ(GSHGEOM_DEVICE_QUADRATIC, 2);
  EXPECT_EQ(GSHGEOM_DEVICE_CUBIC, 3);
  EXPECT_EQ(GSHGEOM_DEVICE_CLOSE, 4);
  EXPECT_EQ(GSHGEOM_DEVICE_VERB_COUNT, 5);
  EXPECT_EQ(GSHGEOM_DEVICE_ONE, 64);
}

TEST(DevicePath, ReadbackOfEveryVerb) {
  GSHGEOM_Device_Path * p = fixture_dp();
  // M L Q C Z M L Z: the second M is stored by the builder.
  ASSERT_EQ(gshgeom_device_path_element_count(p), 8u);
  GSHGEOM_Device_Element e;
  ASSERT_EQ(gshgeom_device_path_element(p, 0, &e), GSHGEOM_OK);
  EXPECT_EQ(e.verb, GSHGEOM_DEVICE_MOVE);
  EXPECT_EQ(e.count, 1u);
  EXPECT_EQ(e.p[0].x, 64);
  EXPECT_EQ(e.p[0].y, -128);
  EXPECT_EQ(e.p[1].x, 0);
  ASSERT_EQ(gshgeom_device_path_element(p, 2, &e), GSHGEOM_OK);
  EXPECT_EQ(e.verb, GSHGEOM_DEVICE_QUADRATIC);
  EXPECT_EQ(e.count, 2u);
  EXPECT_EQ(e.p[0].x, 320);
  EXPECT_EQ(e.p[0].y, 416);
  EXPECT_EQ(e.p[1].x, 448);
  EXPECT_EQ(e.p[1].y, 512);
  ASSERT_EQ(gshgeom_device_path_element(p, 3, &e), GSHGEOM_OK);
  EXPECT_EQ(e.verb, GSHGEOM_DEVICE_CUBIC);
  EXPECT_EQ(e.count, 3u);
  EXPECT_EQ(e.p[0].x, 576);
  EXPECT_EQ(e.p[1].y, 768);
  EXPECT_EQ(e.p[2].x, 832);
  EXPECT_EQ(e.p[2].y, -896);
  ASSERT_EQ(gshgeom_device_path_element(p, 4, &e), GSHGEOM_OK);
  EXPECT_EQ(e.verb, GSHGEOM_DEVICE_CLOSE);
  EXPECT_EQ(e.count, 1u);
  EXPECT_EQ(e.p[0].x, 64) << "a close stores the subpath's start";
  EXPECT_EQ(e.p[0].y, -128);
  ASSERT_EQ(gshgeom_device_path_element(p, 5, &e), GSHGEOM_OK);
  EXPECT_EQ(e.verb, GSHGEOM_DEVICE_MOVE) << "stored for the line after the close";
  EXPECT_EQ(e.p[0].x, 64);
  EXPECT_EQ(e.p[0].y, -128);
  ASSERT_EQ(gshgeom_device_path_element(p, 6, &e), GSHGEOM_OK);
  EXPECT_EQ(e.verb, GSHGEOM_DEVICE_LINE);
  EXPECT_EQ(e.p[0].x, INT32_MIN);
  EXPECT_EQ(e.p[0].y, INT32_MAX);
  EXPECT_EQ(gshgeom_device_path_element(p, 8, &e), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_device_path_element(nullptr, 0, &e), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_device_path_element(p, 0, nullptr), GSHGEOM_ERR_INVALID);
  gshgeom_device_path_destroy(p);
}

TEST(DevicePath, TheBuildersRefuseWhatHasNoSubpath) {
  GSHGEOM_Device_Path * p = make_dp();
  Snapshot empty(p);
  EXPECT_EQ(gshgeom_device_path_line_to(p, D(1, 2)), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_device_path_quadratic_to(p, D(1, 2), D(3, 4)), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_device_path_cubic_to(p, D(1, 2), D(3, 4), D(5, 6)), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_device_path_close(p), GSHGEOM_ERR_INVALID);
  EXPECT_TRUE(Snapshot(p) == empty);
  EXPECT_EQ(gshgeom_device_path_move_to(nullptr, D(0, 0)), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_device_path_line_to(nullptr, D(0, 0)), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_device_path_quadratic_to(nullptr, D(0, 0), D(0, 0)), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_device_path_cubic_to(nullptr, D(0, 0), D(0, 0), D(0, 0)), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_device_path_close(nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_device_path_create(nullptr, nullptr, nullptr), GSHGEOM_ERR_INVALID);
  gshgeom_device_path_destroy(nullptr);
  gshgeom_device_path_destroy(p);
  EXPECT_EQ(gshgeom_device_path_element_count(nullptr), 0u);
}

TEST(DevicePath, EveryDrawingVerbAfterACloseStoresAMove) {
  for (int verb = 0; verb < 3; verb++) {
    GSHGEOM_Device_Path * p = make_dp();
    ASSERT_EQ(gshgeom_device_path_move_to(p, D(7, 9)), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_device_path_line_to(p, D(70, 90)), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_device_path_close(p), GSHGEOM_OK);
    GSHGEOM_Result r = verb == 0   ? gshgeom_device_path_line_to(p, D(1, 1))
        : verb == 1                ? gshgeom_device_path_quadratic_to(p, D(1, 1), D(2, 2))
                                   : gshgeom_device_path_cubic_to(p, D(1, 1), D(2, 2), D(3, 3));
    ASSERT_EQ(r, GSHGEOM_OK);
    ASSERT_EQ(gshgeom_device_path_element_count(p), 5u);
    GSHGEOM_Device_Element e;
    ASSERT_EQ(gshgeom_device_path_element(p, 3, &e), GSHGEOM_OK);
    EXPECT_EQ(e.verb, GSHGEOM_DEVICE_MOVE);
    EXPECT_EQ(e.p[0].x, 7);
    EXPECT_EQ(e.p[0].y, 9);
    // A second drawing verb stores no further move.
    ASSERT_EQ(gshgeom_device_path_line_to(p, D(5, 5)), GSHGEOM_OK);
    EXPECT_EQ(gshgeom_device_path_element_count(p), 6u);
    // A move is explicit, and a new start replaces the old.
    ASSERT_EQ(gshgeom_device_path_close(p), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_device_path_move_to(p, D(-3, -4)), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_device_path_close(p), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_device_path_line_to(p, D(0, 0)), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_device_path_element(p, 9, &e), GSHGEOM_OK);
    EXPECT_EQ(e.verb, GSHGEOM_DEVICE_MOVE);
    EXPECT_EQ(e.p[0].x, -3);
    EXPECT_EQ(e.p[0].y, -4);
    gshgeom_device_path_destroy(p);
  }
  // A second close in a row is a second close, with no move between.
  GSHGEOM_Device_Path * p = make_dp();
  ASSERT_EQ(gshgeom_device_path_move_to(p, D(1, 1)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_device_path_close(p), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_device_path_close(p), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_device_path_element_count(p), 3u);
  gshgeom_device_path_destroy(p);
}

TEST(DevicePathLimits, TheElementCapIsExact) {
  for (size_t cap : {1u, 2u, 5u, 64u}) {
    GSHGEOM_Limits l;
    gshgeom_limits_default(&l);
    l.max_elements = cap;
    GSHGEOM_Device_Path * p = make_dp(nullptr, &l);
    ASSERT_EQ(gshgeom_device_path_move_to(p, D(0, 0)), GSHGEOM_OK);
    for (size_t i = 1; i < cap; i++) {
      ASSERT_EQ(gshgeom_device_path_line_to(p, D((int32_t)i, 0)), GSHGEOM_OK) << cap << " " << i;
    }
    EXPECT_EQ(gshgeom_device_path_element_count(p), cap);
    Snapshot full(p);
    EXPECT_EQ(gshgeom_device_path_line_to(p, D(9, 9)), GSHGEOM_ERR_LIMIT) << cap;
    EXPECT_EQ(gshgeom_device_path_move_to(p, D(9, 9)), GSHGEOM_ERR_LIMIT);
    EXPECT_EQ(gshgeom_device_path_close(p), GSHGEOM_ERR_LIMIT);
    EXPECT_TRUE(Snapshot(p) == full);
    gshgeom_device_path_destroy(p);
  }
  // A drawing verb after a close needs two elements and is all or nothing.
  GSHGEOM_Limits l;
  gshgeom_limits_default(&l);
  l.max_elements = 3;
  GSHGEOM_Device_Path * p = make_dp(nullptr, &l);
  ASSERT_EQ(gshgeom_device_path_move_to(p, D(0, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_device_path_line_to(p, D(1, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_device_path_close(p), GSHGEOM_OK);
  Snapshot closed(p);
  EXPECT_EQ(gshgeom_device_path_line_to(p, D(2, 2)), GSHGEOM_ERR_LIMIT);
  EXPECT_TRUE(Snapshot(p) == closed) << "no move was stored on its own";
  gshgeom_device_path_destroy(p);
  l.max_elements = 4;
  p = make_dp(nullptr, &l);
  ASSERT_EQ(gshgeom_device_path_move_to(p, D(0, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_device_path_line_to(p, D(1, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_device_path_close(p), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_device_path_line_to(p, D(2, 2)), GSHGEOM_ERR_LIMIT) << "needs two, one is left";
  gshgeom_device_path_destroy(p);
  l.max_elements = 5;
  p = make_dp(nullptr, &l);
  ASSERT_EQ(gshgeom_device_path_move_to(p, D(0, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_device_path_line_to(p, D(1, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_device_path_close(p), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_device_path_line_to(p, D(2, 2)), GSHGEOM_OK) << "exactly two are left";
  EXPECT_EQ(gshgeom_device_path_element_count(p), 5u);
  gshgeom_device_path_destroy(p);
  // Zero means the default.
  l.max_elements = 0;
  p = make_dp(nullptr, &l);
  for (int i = 0; i < 100; i++) {
    ASSERT_EQ(gshgeom_device_path_move_to(p, D(i, i)), GSHGEOM_OK);
  }
  gshgeom_device_path_destroy(p);
}

TEST(DevicePathHash, MatchesTheIndependentImplementation) {
  GSHGEOM_Device_Path * p = fixture_dp();
  EXPECT_EQ(gshgeom_device_path_hash(p), kDevicePathHash);
  EXPECT_EQ(dump_of(p), kDevicePathDump);
  gshgeom_device_path_destroy(p);
  // An empty path and NULL.
  GSHGEOM_Device_Path * e = make_dp();
  GSHGEOM_Hash h;
  gshgeom_hash_init(&h);
  gshgeom_hash_u32(&h, 1);
  gshgeom_hash_u32(&h, 0x50564544u);
  gshgeom_hash_u64(&h, 0);
  EXPECT_EQ(gshgeom_device_path_hash(e), gshgeom_hash_value(&h));
  gshgeom_device_path_destroy(e);
  EXPECT_EQ(gshgeom_device_path_hash(nullptr), GSHGEOM_FNV_OFFSET_BASIS);
}

TEST(DevicePathHash, EveryFieldMatters) {
  // Two paths that differ in one field each hash differently.
  auto hash_of = [](int variant) {
    GSHGEOM_Device_Path * p = make_dp();
    gshgeom_device_path_move_to(p, D(1, 2));
    switch (variant) {
    case 0: gshgeom_device_path_line_to(p, D(3, 4)); break;
    case 1: gshgeom_device_path_line_to(p, D(3, 5)); break;
    case 2: gshgeom_device_path_line_to(p, D(4, 3)); break;
    case 3: gshgeom_device_path_quadratic_to(p, D(3, 4), D(3, 4)); break;
    case 4: gshgeom_device_path_cubic_to(p, D(3, 4), D(3, 4), D(3, 4)); break;
    case 5: gshgeom_device_path_close(p); break;
    case 6: gshgeom_device_path_line_to(p, D(-3, 4)); break;
    default: break;
    }
    uint64_t h = gshgeom_device_path_hash(p);
    gshgeom_device_path_destroy(p);
    return h;
  };
  std::set<uint64_t> seen;
  for (int v = 0; v < 8; v++) {
    EXPECT_TRUE(seen.insert(hash_of(v)).second) << v;
  }
}

TEST(DevicePathDump, AFailingSinkIsReported) {
  GSHGEOM_Device_Path * p = fixture_dp();
  FILE * full = std::fopen("/dev/full", "w");
  if (!full) {
    gshgeom_device_path_destroy(p);
    GTEST_SKIP() << "no /dev/full on this host";
  }
  EXPECT_EQ(gshgeom_device_path_dump(p, full), GSHGEOM_ERR_IO);
  std::fclose(full);
  FILE * ro = std::fopen("/dev/null", "r");
  ASSERT_NE(ro, nullptr);
  EXPECT_EQ(gshgeom_device_path_dump(p, ro), GSHGEOM_ERR_IO);
  std::fclose(ro);
  EXPECT_EQ(gshgeom_device_path_dump(p, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_device_path_dump(nullptr, stdout), GSHGEOM_ERR_INVALID);
  gshgeom_device_path_destroy(p);
}

TEST(DevicePathDump, EveryWriteIsChecked) {
  GSHGEOM_Device_Path * p = fixture_dp();
  GSHGEOM_EVERY_WRITE_CHECKED(gshgeom_device_path_dump, p, "device path");
  gshgeom_device_path_destroy(p);
}

TEST(DevicePathDump, TheHashIsAlwaysSixteenDigits) {
  // Some hash among these begins with a zero digit, and prints with it.
  bool leading_zero = false;
  for (int i = 0; i < 400; i++) {
    GSHGEOM_Device_Path * p = make_dp();
    ASSERT_EQ(gshgeom_device_path_move_to(p, D(i, 0)), GSHGEOM_OK);
    uint64_t h = gshgeom_device_path_hash(p);
    leading_zero = leading_zero || h < (UINT64_C(1) << 60);
    std::string d = dump_of(p);
    size_t at = d.find("hash=");
    ASSERT_NE(at, std::string::npos);
    size_t end = d.find('\n', at);
    EXPECT_EQ(end - (at + 5), 16u) << d;
    gshgeom_device_path_destroy(p);
  }
  EXPECT_TRUE(leading_zero) << "the search found no hash with a leading zero";
}

TEST(DevicePathDump, DoesNotDependOnTheLocale) {
  GSHGEOM_Device_Path * p = fixture_dp();
  const char * old = std::setlocale(LC_ALL, nullptr);
  std::string saved = old ? old : "C";
  for (const char * name : {"de_DE.UTF-8", "fr_FR.UTF-8", "de_DE", "C.UTF-8"}) {
    if (std::setlocale(LC_ALL, name)) {
      EXPECT_EQ(dump_of(p), kDevicePathDump) << name;
    }
  }
  std::setlocale(LC_ALL, saved.c_str());
  gshgeom_device_path_destroy(p);
}

namespace {

// The fixture, then forty more lines, so that both arrays grow more than once.
std::vector<DStep> long_steps() {
  std::vector<DStep> steps = device_steps();
  for (int i = 0; i < 40; i++) {
    steps.push_back([i](GSHGEOM_Device_Path * p) {
      return gshgeom_device_path_line_to(p, D(1000 + i, -i));
    });
  }
  return steps;
}

}  // namespace

TEST(DevicePathAllocation, AFailedAllocationIsOomAndChangesNothing) {
  uint64_t want_hash = 0;
  tracker::Tracker count;
  {
    GSHGEOM_Allocator a = tracker::allocator_for(&count);
    GSHGEOM_Device_Path * p = make_dp(&a);
    for (const DStep & s : long_steps()) ASSERT_EQ(s(p), GSHGEOM_OK);
    want_hash = gshgeom_device_path_hash(p);
    gshgeom_device_path_destroy(p);
    EXPECT_TRUE(count.live.empty());
  }
  ASSERT_GE(count.calls, 6);
  for (int persistent = 0; persistent < 2; persistent++) {
    for (long k = 0; k < count.calls; k++) {
      tracker::Tracker t;
      t.fail_from = k;
      t.persistent = persistent != 0;
      GSHGEOM_Allocator a = tracker::allocator_for(&t);
      GSHGEOM_Device_Path * p = nullptr;
      GSHGEOM_Result r = gshgeom_device_path_create(&a, nullptr, &p);
      if (r != GSHGEOM_OK) {
        EXPECT_EQ(r, GSHGEOM_ERR_OOM);
        EXPECT_EQ(p, nullptr);
        EXPECT_TRUE(t.live.empty());
        continue;
      }
      std::vector<DStep> steps = long_steps();
      size_t done = 0;
      for (const DStep & s : steps) {
        Snapshot before(p);
        GSHGEOM_Result sr = s(p);
        if (sr == GSHGEOM_OK) {
          done++;
          continue;
        }
        ASSERT_EQ(sr, GSHGEOM_ERR_OOM) << "k=" << k;
        EXPECT_TRUE(Snapshot(p) == before) << "k=" << k << " step " << done;
        if (!t.persistent) {
          ASSERT_EQ(s(p), GSHGEOM_OK) << "retry, k=" << k;
          done++;
        } else {
          break;
        }
      }
      if (!t.persistent) {
        EXPECT_EQ(done, steps.size());
        EXPECT_EQ(gshgeom_device_path_hash(p), want_hash) << "k=" << k;
      }
      gshgeom_device_path_destroy(p);
      EXPECT_TRUE(t.live.empty()) << "k=" << k << " leaked";
    }
  }
}

// ---------------------------------------------------------------------------
// A model path under a map
// ---------------------------------------------------------------------------

TEST(DeviceFromModel, MapsThenRoundsEachPoint) {
  GSHGEOM_Path * mp = model_fixture();
  const double m[4] = {2.0, 0.0, 0.0, 4.0};
  const double t[2] = {0.2421875, 0.0234375};
  GSHGEOM_Map map = matrix_map(m, t);
  GSHGEOM_Device_Path * dp = nullptr;
  ASSERT_EQ(gshgeom_device_path_from_model(mp, &map, nullptr, nullptr, &dp), GSHGEOM_OK);
  ASSERT_NE(dp, nullptr);
  EXPECT_EQ(gshgeom_device_path_hash(dp), kModelDeviceHash);
  EXPECT_EQ(dump_of(dp), kModelDeviceDump);
  gshgeom_device_path_destroy(dp);
  gshgeom_path_destroy(mp);
}

TEST(DeviceFromModel, TheMapIsAppliedNotJustAScale) {
  // A y flip and a translation: (x, y) -> (x + 10, -y + 20), in pixels.
  GSHGEOM_Map map;
  ASSERT_EQ(gshgeom_map_design_to_model(1.0, 10.0, 20.0, &map), GSHGEOM_OK);
  GSHGEOM_Path * mp = make_model();
  ASSERT_EQ(gshgeom_path_move_to(mp, P(1.0, 2.0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(mp, P(-0.5, 3.0)), GSHGEOM_OK);
  GSHGEOM_Device_Path * dp = nullptr;
  ASSERT_EQ(gshgeom_device_path_from_model(mp, &map, nullptr, nullptr, &dp), GSHGEOM_OK);
  GSHGEOM_Device_Element e;
  ASSERT_EQ(gshgeom_device_path_element(dp, 0, &e), GSHGEOM_OK);
  EXPECT_EQ(e.p[0].x, 11 * 64);
  EXPECT_EQ(e.p[0].y, 18 * 64);
  ASSERT_EQ(gshgeom_device_path_element(dp, 1, &e), GSHGEOM_OK);
  EXPECT_EQ(e.p[0].x, 9 * 64 + 32);
  EXPECT_EQ(e.p[0].y, 17 * 64);
  gshgeom_device_path_destroy(dp);
  gshgeom_path_destroy(mp);
  // The same point through the point function.
  Dp d;
  ASSERT_EQ(gshgeom_device_point(&map, P(1.0, 2.0), &d), GSHGEOM_OK);
  EXPECT_EQ(d.x, 11 * 64);
  EXPECT_EQ(d.y, 18 * 64);
}

TEST(DeviceFromModel, RefusesWhatIsNot2DOrNull) {
  GSHGEOM_Path * mp = model_fixture();
  GSHGEOM_Map map = identity2();
  GSHGEOM_Device_Path * dp = nullptr;
  GSHGEOM_Map map3;
  ASSERT_EQ(gshgeom_map_identity(3, &map3), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_device_path_from_model(mp, &map3, nullptr, nullptr, &dp), GSHGEOM_ERR_INVALID);
  GSHGEOM_Path * p3 = make_model(3);
  ASSERT_EQ(gshgeom_path_move_to(p3, GSHGEOM_Point{1, 2, 3}), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_device_path_from_model(p3, &map, nullptr, nullptr, &dp), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_device_path_from_model(p3, &map3, nullptr, nullptr, &dp), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_device_path_from_model(nullptr, &map, nullptr, nullptr, &dp), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_device_path_from_model(mp, nullptr, nullptr, nullptr, &dp), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_device_path_from_model(mp, &map, nullptr, nullptr, nullptr), GSHGEOM_ERR_INVALID);
  GSHGEOM_Map bad = map;
  bad.m[0][0] = kNaN;
  EXPECT_EQ(gshgeom_device_path_from_model(mp, &bad, nullptr, nullptr, &dp), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(dp, nullptr);
  Dp d = D(5, 6);
  EXPECT_EQ(gshgeom_device_point(&map3, P(1, 1), &d), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_device_point(nullptr, P(1, 1), &d), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_device_point(&map, P(kNaN, 1), &d), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_device_point(&map, P(1, 1), nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(d.x, 5);
  EXPECT_EQ(d.y, 6);
  gshgeom_path_destroy(p3);
  gshgeom_path_destroy(mp);
}

TEST(DeviceFromModel, AnEmptyPathGivesAnEmptyDevicePath) {
  GSHGEOM_Path * mp = make_model();
  GSHGEOM_Map map = identity2();
  GSHGEOM_Device_Path * dp = nullptr;
  ASSERT_EQ(gshgeom_device_path_from_model(mp, &map, nullptr, nullptr, &dp), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_device_path_element_count(dp), 0u);
  gshgeom_device_path_destroy(dp);
  gshgeom_path_destroy(mp);
}

TEST(DeviceFromModel, OutOfRangeIsRefusedNotWrappedOrClipped) {
  GSHGEOM_Map map = identity2();
  const double top = 2147483647.0 / 64.0;  // the largest coordinate that fits
  struct Case {
    double x, y;
    bool fits;
  };
  const Case cases[] = {
      {top, 0.0, true},
      {0.0, top, true},
      {-2147483648.0 / 64.0, 0.0, true},
      {std::nextafter(top, kInf), 0.0, true},  // still rounds to INT32_MAX
      {(2147483647.5) / 64.0, 0.0, false},
      {0.0, (2147483647.5) / 64.0, false},
      {-2147483649.0 / 64.0, 0.0, false},
      {0.0, -2147483649.0 / 64.0, false},
      {4294967296.0 / 64.0, 0.0, false},  // would wrap to 0
      {1e30, 0.0, false},
      {0.0, -1e30, false},
  };
  for (const Case & c : cases) {
    GSHGEOM_Path * mp = make_model();
    ASSERT_EQ(gshgeom_path_move_to(mp, P(0, 0)), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_path_line_to(mp, P(c.x, c.y)), GSHGEOM_OK);
    GSHGEOM_Device_Path * dp = nullptr;
    GSHGEOM_Result r = gshgeom_device_path_from_model(mp, &map, nullptr, nullptr, &dp);
    if (c.fits) {
      ASSERT_EQ(r, GSHGEOM_OK) << c.x << " " << c.y;
      GSHGEOM_Device_Element e;
      ASSERT_EQ(gshgeom_device_path_element(dp, 1, &e), GSHGEOM_OK);
      int32_t want;
      ASSERT_EQ(gshgeom_round_26_6(c.x, &want), GSHGEOM_OK);
      EXPECT_EQ(e.p[0].x, want);
      gshgeom_device_path_destroy(dp);
    } else {
      EXPECT_EQ(r, GSHGEOM_ERR_LIMIT) << c.x << " " << c.y;
      EXPECT_EQ(dp, nullptr) << "nothing returned";
    }
    gshgeom_path_destroy(mp);
  }
  // The range is the range after the map: a small path under a large scale.
  GSHGEOM_Map big;
  ASSERT_EQ(gshgeom_map_scale(2, 1e9, 1e9, 1.0, &big), GSHGEOM_OK);
  GSHGEOM_Path * mp = make_model();
  ASSERT_EQ(gshgeom_path_move_to(mp, P(0, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(mp, P(1.0, 1.0)), GSHGEOM_OK);
  GSHGEOM_Device_Path * dp = nullptr;
  EXPECT_EQ(gshgeom_device_path_from_model(mp, &big, nullptr, nullptr, &dp), GSHGEOM_ERR_LIMIT);
  gshgeom_path_destroy(mp);
  Dp d = D(5, 6);
  EXPECT_EQ(gshgeom_device_point(&big, P(1.0, 1.0), &d), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(d.x, 5) << "written only on success";
  EXPECT_EQ(d.y, 6);
}

TEST(DeviceFromModel, AControlPointOutOfRangeIsRefusedEvenIfTheCurveIsInside) {
  GSHGEOM_Map map = identity2();
  GSHGEOM_Path * mp = make_model();
  ASSERT_EQ(gshgeom_path_move_to(mp, P(0, 0)), GSHGEOM_OK);
  // The curve peaks at x = 2.5e7 pixels, inside the range; the control point is
  // at 5e7 pixels, outside it.
  ASSERT_EQ(gshgeom_path_quadratic_to(mp, P(5e7, 0), P(0, 1)), GSHGEOM_OK);
  GSHGEOM_Device_Path * dp = nullptr;
  EXPECT_EQ(gshgeom_device_path_from_model(mp, &map, nullptr, nullptr, &dp), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(dp, nullptr);
  gshgeom_path_destroy(mp);
  // The same for a cubic's second control point and a rational curve's.
  mp = make_model();
  ASSERT_EQ(gshgeom_path_move_to(mp, P(0, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_cubic_to(mp, P(1, 1), P(0, -5e7), P(2, 2)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_device_path_from_model(mp, &map, nullptr, nullptr, &dp), GSHGEOM_ERR_LIMIT);
  gshgeom_path_destroy(mp);
  mp = make_model();
  ASSERT_EQ(gshgeom_path_move_to(mp, P(0, 0)), GSHGEOM_OK);
  const double w[3] = {1.0, 1.0, 1.0};
  ASSERT_EQ(gshgeom_path_rational_quadratic_to(mp, P(5e7, 0), P(0, 1), w), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_device_path_from_model(mp, &map, nullptr, nullptr, &dp), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(dp, nullptr);
  gshgeom_path_destroy(mp);
  mp = make_model();
  ASSERT_EQ(gshgeom_path_move_to(mp, P(0, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_rational_quadratic_to(mp, P(0, -5e7), P(0, 1), w), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_device_path_from_model(mp, &map, nullptr, nullptr, &dp), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(dp, nullptr);
  gshgeom_path_destroy(mp);
}

TEST(DeviceFromModel, TheElementLimitApplies) {
  GSHGEOM_Path * mp = model_fixture();  // 7 device elements: M L Q C Z M L
  GSHGEOM_Map map = identity2();
  GSHGEOM_Limits l;
  gshgeom_limits_default(&l);
  GSHGEOM_Device_Path * dp = nullptr;
  l.max_elements = 7;
  ASSERT_EQ(gshgeom_device_path_from_model(mp, &map, nullptr, &l, &dp), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_device_path_element_count(dp), 7u);
  gshgeom_device_path_destroy(dp);
  dp = nullptr;
  l.max_elements = 6;
  EXPECT_EQ(gshgeom_device_path_from_model(mp, &map, nullptr, &l, &dp), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(dp, nullptr);
  gshgeom_path_destroy(mp);
}

TEST(DeviceFromModel, AllocationFailureAtEveryCall) {
  GSHGEOM_Path * mp = model_fixture();
  const double m[4] = {2.0, 0.0, 0.0, 4.0};
  const double t[2] = {0.2421875, 0.0234375};
  GSHGEOM_Map map = matrix_map(m, t);
  tracker::Tracker count;
  {
    GSHGEOM_Allocator a = tracker::allocator_for(&count);
    GSHGEOM_Device_Path * dp = nullptr;
    ASSERT_EQ(gshgeom_device_path_from_model(mp, &map, &a, nullptr, &dp), GSHGEOM_OK);
    gshgeom_device_path_destroy(dp);
  }
  ASSERT_GT(count.calls, 2);
  for (long k = 0; k < count.calls; k++) {
    tracker::Tracker tr;
    tr.fail_from = k;
    tr.persistent = true;
    GSHGEOM_Allocator a = tracker::allocator_for(&tr);
    GSHGEOM_Device_Path * dp = nullptr;
    EXPECT_EQ(gshgeom_device_path_from_model(mp, &map, &a, nullptr, &dp), GSHGEOM_ERR_OOM) << k;
    EXPECT_EQ(dp, nullptr);
    EXPECT_TRUE(tr.live.empty()) << "k=" << k << " leaked";
  }
  gshgeom_path_destroy(mp);
}

// ---------------------------------------------------------------------------
// Rational curves
// ---------------------------------------------------------------------------

namespace {

// All the points of a device path that is only a move and lines, in 26.6.
std::vector<std::pair<double, double>> polyline_of(const GSHGEOM_Device_Path * dp) {
  std::vector<std::pair<double, double>> pts;
  size_t n = gshgeom_device_path_element_count(dp);
  for (size_t i = 0; i < n; i++) {
    GSHGEOM_Device_Element e;
    EXPECT_EQ(gshgeom_device_path_element(dp, i, &e), GSHGEOM_OK);
    EXPECT_TRUE(e.verb == GSHGEOM_DEVICE_MOVE || e.verb == GSHGEOM_DEVICE_LINE) << "verb " << e.verb;
    pts.push_back({(double)e.p[0].x, (double)e.p[0].y});
  }
  return pts;
}

double dist_to_segment(double px, double py, std::pair<double, double> a, std::pair<double, double> b) {
  double dx = b.first - a.first, dy = b.second - a.second;
  double len2 = dx * dx + dy * dy;
  double u = len2 > 0 ? ((px - a.first) * dx + (py - a.second) * dy) / len2 : 0.0;
  u = u < 0 ? 0 : (u > 1 ? 1 : u);
  double qx = a.first + u * dx - px, qy = a.second + u * dy - py;
  return std::sqrt(qx * qx + qy * qy);
}

double dist_to_polyline(double px, double py, const std::vector<std::pair<double, double>> & pl) {
  double best = 1e300;
  for (size_t i = 0; i + 1 < pl.size(); i++) {
    best = std::min(best, dist_to_segment(px, py, pl[i], pl[i + 1]));
  }
  return best;
}

// The largest distance from the model curve (scaled to device units) to the
// polyline, over many parameters of the curve itself.
double curve_to_polyline(const GSHGEOM_Curve & c, double scale, const std::vector<std::pair<double, double>> & pl) {
  double worst = 0;
  const int n = 4000;
  for (int j = 0; j <= n; j++) {
    GSHGEOM_Point q;
    EXPECT_EQ(gshgeom_curve_eval(&c, (double)j / n, &q), GSHGEOM_OK);
    worst = std::max(worst, dist_to_polyline(q.x * scale * 64.0, q.y * scale * 64.0, pl));
  }
  return worst;
}

}  // namespace

TEST(DeviceRational, AQuarterCircleStaysWithinTheToleranceAndEndsWhereItShould) {
  for (double radius : {1.0, 10.0, 100.0, 1000.0}) {
    const double w[3] = {1.0, std::sqrt(0.5), 1.0};
    GSHGEOM_Path * mp = make_model();
    ASSERT_EQ(gshgeom_path_move_to(mp, P(radius, 0)), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_path_rational_quadratic_to(mp, P(radius, radius), P(0, radius), w), GSHGEOM_OK);
    GSHGEOM_Map map = identity2();
    GSHGEOM_Device_Path * dp = nullptr;
    ASSERT_EQ(gshgeom_device_path_from_model(mp, &map, nullptr, nullptr, &dp), GSHGEOM_OK);
    auto pl = polyline_of(dp);
    ASSERT_GE(pl.size(), 2u);
    EXPECT_EQ(pl.front().first, radius * 64);
    EXPECT_EQ(pl.front().second, 0.0);
    EXPECT_EQ(pl.back().first, 0.0);
    EXPECT_EQ(pl.back().second, radius * 64) << "the ends are the mapped end points";
    if (radius >= 10) {
      EXPECT_GT(pl.size(), 3u) << "a real arc needs several lines, radius " << radius;
    }
    if (radius == 1000.0) {
      // A quarter circle of 64,000 units at 2 units of deviation needs about a
      // hundred pieces; far more means the subdivision is not stopping.
      EXPECT_GE(pl.size(), 64u);
      EXPECT_LE(pl.size(), 300u);
    }
    double R = radius * 64;
    for (auto & v : pl) {
      EXPECT_LE(std::fabs(std::hypot(v.first, v.second) - R), 1.0) << "a vertex is on the circle, to rounding";
    }
    // True circle points, not the library's own evaluation.
    double worst = 0;
    for (int j = 0; j <= 4000; j++) {
      double a = (M_PI / 2) * j / 4000;
      worst = std::max(worst, dist_to_polyline(R * std::cos(a), R * std::sin(a), pl));
    }
    EXPECT_LE(worst, 4.0) << "within 1/16 px, radius " << radius;
    gshgeom_device_path_destroy(dp);
    gshgeom_path_destroy(mp);
  }
}

TEST(DeviceRational, TheToleranceIsInDevicePixels) {
  // The same unit-radius arc under a small and a large scale: the large one
  // takes many lines, the small one a few, both within 1/16 px of the arc.
  const double w[3] = {1.0, std::sqrt(0.5), 1.0};
  GSHGEOM_Path * mp = make_model();
  ASSERT_EQ(gshgeom_path_move_to(mp, P(1, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_rational_quadratic_to(mp, P(1, 1), P(0, 1), w), GSHGEOM_OK);
  size_t previous = 0;
  for (double scale : {0.5, 8.0, 100.0, 2000.0}) {
    GSHGEOM_Map map;
    ASSERT_EQ(gshgeom_map_scale(2, scale, scale, 1.0, &map), GSHGEOM_OK);
    GSHGEOM_Device_Path * dp = nullptr;
    ASSERT_EQ(gshgeom_device_path_from_model(mp, &map, nullptr, nullptr, &dp), GSHGEOM_OK);
    size_t n = gshgeom_device_path_element_count(dp);
    EXPECT_GT(n, previous) << scale;
    previous = n;
    double R = scale * 64;
    double worst = 0;
    auto pl = polyline_of(dp);
    for (int j = 0; j <= 3000; j++) {
      double a = (M_PI / 2) * j / 3000;
      worst = std::max(worst, dist_to_polyline(R * std::cos(a), R * std::sin(a), pl));
    }
    EXPECT_LE(worst, 4.0) << "scale " << scale;
    gshgeom_device_path_destroy(dp);
  }
  gshgeom_path_destroy(mp);
}

TEST(DeviceRational, EqualWeightsAgreeWithThePolynomialCurve) {
  GSHGEOM_Map map = identity2();
  const double w[3] = {3.0, 3.0, 3.0};
  const double wc[4] = {0.25, 0.25, 0.25, 0.25};
  GSHGEOM_Path * mp = make_model();
  ASSERT_EQ(gshgeom_path_move_to(mp, P(0, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_rational_quadratic_to(mp, P(40, 90), P(80, 0), w), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_rational_cubic_to(mp, P(100, -50), P(120, 50), P(160, 10), wc), GSHGEOM_OK);
  GSHGEOM_Device_Path * dp = nullptr;
  ASSERT_EQ(gshgeom_device_path_from_model(mp, &map, nullptr, nullptr, &dp), GSHGEOM_OK);
  auto pl = polyline_of(dp);
  GSHGEOM_Curve q, c;
  ASSERT_EQ(gshgeom_curve_quadratic(2, P(0, 0), P(40, 90), P(80, 0), &q), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_curve_cubic(2, P(80, 0), P(100, -50), P(120, 50), P(160, 10), &c), GSHGEOM_OK);
  EXPECT_LE(curve_to_polyline(q, 1.0, pl), 4.0);
  EXPECT_LE(curve_to_polyline(c, 1.0, pl), 4.0);
  gshgeom_device_path_destroy(dp);
  gshgeom_path_destroy(mp);
}

TEST(DeviceRational, VeryUnequalWeightsStayWithinTheTolerance) {
  GSHGEOM_Map map = identity2();
  const double wq[3] = {1.0, 1000.0, 1.0};
  const double wc[4] = {1.0, 0.001, 500.0, 1.0};
  GSHGEOM_Path * mp = make_model();
  ASSERT_EQ(gshgeom_path_move_to(mp, P(0, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_rational_quadratic_to(mp, P(40, 90), P(80, 0), wq), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_rational_cubic_to(mp, P(100, -50), P(120, 50), P(160, 10), wc), GSHGEOM_OK);
  GSHGEOM_Device_Path * dp = nullptr;
  ASSERT_EQ(gshgeom_device_path_from_model(mp, &map, nullptr, nullptr, &dp), GSHGEOM_OK);
  auto pl = polyline_of(dp);
  GSHGEOM_Curve q, c;
  ASSERT_EQ(gshgeom_curve_rational_quadratic(2, P(0, 0), P(40, 90), P(80, 0), wq, &q), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_curve_rational_cubic(2, P(80, 0), P(100, -50), P(120, 50), P(160, 10), wc, &c), GSHGEOM_OK);
  EXPECT_LE(curve_to_polyline(q, 1.0, pl), 4.0);
  EXPECT_LE(curve_to_polyline(c, 1.0, pl), 4.0);
  EXPECT_EQ(pl.back().first, 160.0 * 64);
  gshgeom_device_path_destroy(dp);
  gshgeom_path_destroy(mp);
}

TEST(DeviceRational, ACollinearOvershootIsNotFlat) {
  // The control point lies on the chord's line but beyond its end: the curve runs
  // out to x = 133 and comes back to 100, so its distance from the chord's
  // *line* is zero and from the chord is 33 pixels.
  GSHGEOM_Map map = identity2();
  const double w[3] = {1.0, 1.0, 1.0};
  GSHGEOM_Path * mp = make_model();
  ASSERT_EQ(gshgeom_path_move_to(mp, P(0, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_rational_quadratic_to(mp, P(200, 0), P(100, 0), w), GSHGEOM_OK);
  GSHGEOM_Device_Path * dp = nullptr;
  ASSERT_EQ(gshgeom_device_path_from_model(mp, &map, nullptr, nullptr, &dp), GSHGEOM_OK);
  auto pl = polyline_of(dp);
  EXPECT_GT(pl.size(), 2u);
  double far = 0;
  for (auto & v : pl) {
    far = std::max(far, v.first);
  }
  EXPECT_GE(far, 130.0 * 64) << "the polyline reaches the overshoot";
  gshgeom_device_path_destroy(dp);
  gshgeom_path_destroy(mp);
}

TEST(DeviceRational, EveryControlPointOfACubicCounts) {
  // One interior control point on the chord, the other far from it; and the
  // reverse.
  GSHGEOM_Map map = identity2();
  const double w[4] = {1.0, 2.0, 1.5, 1.0};
  for (int which = 0; which < 2; which++) {
    GSHGEOM_Point c1 = which ? P(30, 0) : P(30, 60);
    GSHGEOM_Point c2 = which ? P(70, 60) : P(70, 0);
    GSHGEOM_Path * mp = make_model();
    ASSERT_EQ(gshgeom_path_move_to(mp, P(0, 0)), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_path_rational_cubic_to(mp, c1, c2, P(100, 0), w), GSHGEOM_OK);
    GSHGEOM_Device_Path * dp = nullptr;
    ASSERT_EQ(gshgeom_device_path_from_model(mp, &map, nullptr, nullptr, &dp), GSHGEOM_OK);
    auto pl = polyline_of(dp);
    GSHGEOM_Curve c;
    ASSERT_EQ(gshgeom_curve_rational_cubic(2, P(0, 0), c1, c2, P(100, 0), w, &c), GSHGEOM_OK);
    EXPECT_GT(pl.size(), 4u) << which;
    EXPECT_LE(curve_to_polyline(c, 1.0, pl), 4.0) << which;
    gshgeom_device_path_destroy(dp);
    gshgeom_path_destroy(mp);
  }
}

TEST(DeviceRational, ADegreeElevatedCircleIsACubicOfTheSameShape) {
  // The quarter circle as a rational quadratic, elevated to a rational cubic
  // (homogeneous control points (1/3) H0 + (2/3) H1 and (2/3) H1 + (1/3) H2).
  const double s = std::sqrt(0.5);
  const double R = 1000.0;
  const double wc = (1.0 + 2.0 * s) / 3.0;
  const double w[4] = {1.0, wc, wc, 1.0};
  GSHGEOM_Point p0 = P(R, 0), p1 = P(R, R), p2 = P(0, R);
  GSHGEOM_Point c1 = P(((p0.x + 2.0 * s * p1.x) / 3.0) / wc, ((p0.y + 2.0 * s * p1.y) / 3.0) / wc);
  GSHGEOM_Point c2 = P(((2.0 * s * p1.x + p2.x) / 3.0) / wc, ((2.0 * s * p1.y + p2.y) / 3.0) / wc);
  GSHGEOM_Map map = identity2();
  GSHGEOM_Path * mp = make_model();
  ASSERT_EQ(gshgeom_path_move_to(mp, p0), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_rational_cubic_to(mp, c1, c2, p2, w), GSHGEOM_OK);
  GSHGEOM_Device_Path * dp = nullptr;
  ASSERT_EQ(gshgeom_device_path_from_model(mp, &map, nullptr, nullptr, &dp), GSHGEOM_OK);
  auto pl = polyline_of(dp);
  EXPECT_GE(pl.size(), 64u);
  EXPECT_LE(pl.size(), 300u);
  double worst = 0;
  for (int j = 0; j <= 4000; j++) {
    double a = (M_PI / 2) * j / 4000;
    worst = std::max(worst, dist_to_polyline(R * 64 * std::cos(a), R * 64 * std::sin(a), pl));
  }
  EXPECT_LE(worst, 4.0);
  gshgeom_device_path_destroy(dp);
  gshgeom_path_destroy(mp);
}

TEST(DeviceRational, AllocationFailureAtEveryCallLeavesNoPartialResult) {
  const double w[3] = {1.0, std::sqrt(0.5), 1.0};
  const double wc[4] = {1.0, 3.0, 0.5, 1.0};
  GSHGEOM_Path * mp = make_model();
  ASSERT_EQ(gshgeom_path_move_to(mp, P(1000, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_rational_quadratic_to(mp, P(1000, 1000), P(0, 1000), w), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_rational_cubic_to(mp, P(-500, 2000), P(-1000, 500), P(-1000, 0), wc), GSHGEOM_OK);
  GSHGEOM_Map map = identity2();
  tracker::Tracker count;
  size_t elements = 0;
  {
    GSHGEOM_Allocator a = tracker::allocator_for(&count);
    GSHGEOM_Device_Path * dp = nullptr;
    ASSERT_EQ(gshgeom_device_path_from_model(mp, &map, &a, nullptr, &dp), GSHGEOM_OK);
    elements = gshgeom_device_path_element_count(dp);
    gshgeom_device_path_destroy(dp);
  }
  ASSERT_GT(elements, 300u) << "enough lines for the arrays to grow more than once";
  ASSERT_GT(count.calls, 4);
  for (long k = 0; k < count.calls; k++) {
    for (int persistent = 0; persistent < 2; persistent++) {
      tracker::Tracker t;
      t.fail_from = k;
      t.persistent = persistent != 0;
      GSHGEOM_Allocator a = tracker::allocator_for(&t);
      GSHGEOM_Device_Path * dp = nullptr;
      EXPECT_EQ(gshgeom_device_path_from_model(mp, &map, &a, nullptr, &dp), GSHGEOM_ERR_OOM) << k;
      EXPECT_EQ(dp, nullptr);
      EXPECT_TRUE(t.live.empty()) << "k=" << k << " leaked";
    }
  }
  gshgeom_path_destroy(mp);
}

TEST(DeviceRational, ADegenerateCurveIsNotAnError) {
  GSHGEOM_Map map = identity2();
  const double w[3] = {1.0, 2.0, 1.0};
  GSHGEOM_Path * mp = make_model();
  ASSERT_EQ(gshgeom_path_move_to(mp, P(5, 5)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_rational_quadratic_to(mp, P(5, 5), P(5, 5), w), GSHGEOM_OK);
  // A closed loop: the chord has zero length.
  ASSERT_EQ(gshgeom_path_rational_quadratic_to(mp, P(105, 5), P(5, 5), w), GSHGEOM_OK);
  GSHGEOM_Device_Path * dp = nullptr;
  ASSERT_EQ(gshgeom_device_path_from_model(mp, &map, nullptr, nullptr, &dp), GSHGEOM_OK);
  auto pl = polyline_of(dp);
  ASSERT_GE(pl.size(), 3u);
  GSHGEOM_Curve loop;
  ASSERT_EQ(gshgeom_curve_rational_quadratic(2, P(5, 5), P(105, 5), P(5, 5), w, &loop), GSHGEOM_OK);
  EXPECT_LE(curve_to_polyline(loop, 1.0, std::vector<std::pair<double, double>>(pl.begin() + 1, pl.end())), 4.0);
  gshgeom_device_path_destroy(dp);
  gshgeom_path_destroy(mp);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
