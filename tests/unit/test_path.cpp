/**
 * @file
 *
 * The model path (CAP-2, CAP-13): the builder's rules, readback, bounds, the
 * limits, refusal that leaves the path unchanged, allocation failure at every
 * call, and the hash and the dump pinned by known answers.
 *
 * Copyright 2026 by Corey Pennycuff
 */

#include "test_helpers.h"

#include "../known_answers.h"

#include <ghoti.io/math/math.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <clocale>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <limits>
#include <set>
#include <string>
#include <vector>

namespace {

const double kNaN = std::numeric_limits<double>::quiet_NaN();
const double kInf = std::numeric_limits<double>::infinity();

GSHGEOM_Point P(double x, double y, double z = 0.0) {
  return GSHGEOM_Point{x, y, z};
}

// The fixture of tools/gen-hash-vectors.py: z is (x + y) / 4.
GSHGEOM_Point FP(double x, double y) {
  return P(x, y, (x + y) / 4.0);
}

using Step = std::function<GSHGEOM_Result(GSHGEOM_Path *)>;

std::vector<Step> fixture_steps() {
  static const double wq[3] = {1.0, 0.5, 2.0};
  static const double wc[4] = {1.0, 2.0, 0.5, 4.0};
  return {
      [](GSHGEOM_Path * p) { return gshgeom_path_move_to(p, FP(1.5, -2.25)); },
      [](GSHGEOM_Path * p) { return gshgeom_path_line_to(p, FP(3.0, 4.0)); },
      [](GSHGEOM_Path * p) {
        return gshgeom_path_quadratic_to(p, FP(5.0, 6.5), FP(7.0, 8.0));
      },
      [](GSHGEOM_Path * p) {
        return gshgeom_path_cubic_to(
            p, FP(9.0, 10.0), FP(11.25, 12.0), FP(13.0, 14.0));
      },
      [](GSHGEOM_Path * p) {
        return gshgeom_path_rational_quadratic_to(
            p, FP(15.0, 16.0), FP(17.0, 18.0), wq);
      },
      [](GSHGEOM_Path * p) {
        return gshgeom_path_rational_cubic_to(
            p, FP(19.0, 20.0), FP(21.0, 22.0), FP(23.0, 24.1), wc);
      },
      [](GSHGEOM_Path * p) { return gshgeom_path_close(p); },
  };
}

GSHGEOM_Path * make_path(int dim, const GSHGEOM_Allocator * alloc = nullptr,
    const GSHGEOM_Limits * limits = nullptr) {
  GSHGEOM_Path * p = nullptr;
  EXPECT_EQ(gshgeom_path_create(alloc, limits, dim, &p), GSHGEOM_OK);
  return p;
}

GSHGEOM_Path * fixture_path(int dim) {
  GSHGEOM_Path * p = make_path(dim);
  for (const Step & s : fixture_steps()) {
    EXPECT_EQ(s(p), GSHGEOM_OK);
  }
  return p;
}

std::string dump_of(const GSHGEOM_Path * p) {
  FILE * f = std::tmpfile();
  EXPECT_NE(f, nullptr);
  EXPECT_EQ(gshgeom_path_dump(p, f), GSHGEOM_OK);
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

// A path's identity for "unchanged": its size, its hash and its dump.
struct Snapshot {
  size_t count;
  uint64_t hash;
  std::string dump;
  explicit Snapshot(const GSHGEOM_Path * p)
      : count(gshgeom_path_element_count(p)),
        hash(gshgeom_path_hash(p)),
        dump(dump_of(p)) {}
  bool operator==(const Snapshot & o) const {
    return count == o.count && hash == o.hash && dump == o.dump;
  }
};

}  // namespace

// ---------------------------------------------------------------------------
// Core
// ---------------------------------------------------------------------------

TEST(Core, EveryResultHasADistinctString) {
  std::set<std::string> seen;
  for (int r = 0; r < GSHGEOM_RESULT_COUNT; r++) {
    const char * s = gshgeom_result_string((GSHGEOM_Result)r);
    ASSERT_NE(s, nullptr);
    EXPECT_GT(std::strlen(s), 0u);
    EXPECT_TRUE(seen.insert(s).second) << s;
    EXPECT_STRNE(s, "unknown result") << r;
  }
  EXPECT_STREQ(gshgeom_result_string(GSHGEOM_RESULT_COUNT), "unknown result");
  EXPECT_STREQ(gshgeom_result_string((GSHGEOM_Result)-1), "unknown result");
  EXPECT_STREQ(gshgeom_result_string((GSHGEOM_Result)1000), "unknown result");
  EXPECT_EQ(GSHGEOM_OK, 0);
}

TEST(Core, LimitsDefaultsAndZeroMeansDefault) {
  GSHGEOM_Limits l;
  std::memset(&l, 0xff, sizeof l);
  gshgeom_limits_default(&l);
  EXPECT_EQ(l.max_elements, GSHGEOM_LIMITS_DEFAULT_MAX_ELEMENTS);
  EXPECT_EQ(GSHGEOM_LIMITS_DEFAULT_MAX_ELEMENTS, (size_t)4194304);
  for (unsigned char b : l._reserved) {
    EXPECT_EQ(b, 0);
  }
  gshgeom_limits_default(nullptr);  // ignored
  EXPECT_EQ(gshgeom_allocator_default(), gcu_allocator_default());
}

// ---------------------------------------------------------------------------
// Building
// ---------------------------------------------------------------------------

TEST(Path, CreateChecksItsArguments) {
  GSHGEOM_Path * p = (GSHGEOM_Path *)0x1;
  for (int dim : {-1, 0, 1, 4, 99}) {
    p = (GSHGEOM_Path *)0x1;
    EXPECT_EQ(gshgeom_path_create(nullptr, nullptr, dim, &p),
        GSHGEOM_ERR_INVALID);
    EXPECT_EQ(p, (GSHGEOM_Path *)0x1) << "written only on success";
  }
  EXPECT_EQ(gshgeom_path_create(nullptr, nullptr, 2, nullptr),
      GSHGEOM_ERR_INVALID);
  gshgeom_path_destroy(nullptr);
  for (int dim = 2; dim <= 3; dim++) {
    GSHGEOM_Path * q = make_path(dim);
    EXPECT_EQ(gshgeom_path_dimension(q), dim);
    EXPECT_EQ(gshgeom_path_element_count(q), 0u);
    gshgeom_path_destroy(q);
  }
  EXPECT_EQ(gshgeom_path_dimension(nullptr), 0);
  EXPECT_EQ(gshgeom_path_element_count(nullptr), 0u);
}

TEST(Path, VerbValuesArePinned) {
  EXPECT_EQ(GSHGEOM_VERB_MOVE, 0);
  EXPECT_EQ(GSHGEOM_VERB_LINE, 1);
  EXPECT_EQ(GSHGEOM_VERB_QUADRATIC, 2);
  EXPECT_EQ(GSHGEOM_VERB_CUBIC, 3);
  EXPECT_EQ(GSHGEOM_VERB_RATIONAL_QUADRATIC, 4);
  EXPECT_EQ(GSHGEOM_VERB_RATIONAL_CUBIC, 5);
  EXPECT_EQ(GSHGEOM_VERB_CLOSE, 6);
  EXPECT_EQ(GSHGEOM_VERB_COUNT, 7);
}

TEST(Path, ElementsReadBackAsTheyWereGiven) {
  for (int dim = 2; dim <= 3; dim++) {
    GSHGEOM_Path * p = fixture_path(dim);
    ASSERT_EQ(gshgeom_path_element_count(p), 7u);
    double z = dim == 3 ? 1.0 : 0.0;
    auto zof = [&](double x, double y) { return z * (x + y) / 4.0; };
    GSHGEOM_Element e;

    ASSERT_EQ(gshgeom_path_element(p, 0, &e), GSHGEOM_OK);
    EXPECT_EQ(e.verb, GSHGEOM_VERB_MOVE);
    EXPECT_EQ(e.curve.kind, GSHGEOM_CURVE_LINE);
    EXPECT_EQ(e.curve.dim, dim);
    EXPECT_EQ(e.curve.p[0].x, 1.5);
    EXPECT_EQ(e.curve.p[0].y, -2.25);
    EXPECT_EQ(e.curve.p[0].z, zof(1.5, -2.25));
    EXPECT_EQ(std::memcmp(&e.curve.p[0], &e.curve.p[1], sizeof e.curve.p[0]), 0);

    ASSERT_EQ(gshgeom_path_element(p, 1, &e), GSHGEOM_OK);
    EXPECT_EQ(e.verb, GSHGEOM_VERB_LINE);
    EXPECT_EQ(e.curve.kind, GSHGEOM_CURVE_LINE);
    EXPECT_EQ(e.curve.p[0].x, 1.5);  // the start is the end of the element before
    EXPECT_EQ(e.curve.p[1].x, 3.0);
    EXPECT_EQ(e.curve.p[1].y, 4.0);
    EXPECT_EQ(e.curve.p[1].z, zof(3.0, 4.0));

    ASSERT_EQ(gshgeom_path_element(p, 2, &e), GSHGEOM_OK);
    EXPECT_EQ(e.verb, GSHGEOM_VERB_QUADRATIC);
    EXPECT_EQ(e.curve.kind, GSHGEOM_CURVE_QUADRATIC);
    EXPECT_EQ(e.curve.p[0].x, 3.0);
    EXPECT_EQ(e.curve.p[1].x, 5.0);
    EXPECT_EQ(e.curve.p[1].y, 6.5);
    EXPECT_EQ(e.curve.p[2].x, 7.0);

    ASSERT_EQ(gshgeom_path_element(p, 3, &e), GSHGEOM_OK);
    EXPECT_EQ(e.verb, GSHGEOM_VERB_CUBIC);
    EXPECT_EQ(e.curve.kind, GSHGEOM_CURVE_CUBIC);
    EXPECT_EQ(e.curve.p[0].x, 7.0);
    EXPECT_EQ(e.curve.p[2].x, 11.25);
    EXPECT_EQ(e.curve.p[3].x, 13.0);
    EXPECT_EQ(e.curve.w[0], 0.0);

    ASSERT_EQ(gshgeom_path_element(p, 4, &e), GSHGEOM_OK);
    EXPECT_EQ(e.verb, GSHGEOM_VERB_RATIONAL_QUADRATIC);
    EXPECT_EQ(e.curve.kind, GSHGEOM_CURVE_RATIONAL_QUADRATIC);
    EXPECT_EQ(e.curve.p[0].x, 13.0);
    EXPECT_EQ(e.curve.p[2].x, 17.0);
    EXPECT_EQ(e.curve.w[0], 1.0);
    EXPECT_EQ(e.curve.w[1], 0.5);
    EXPECT_EQ(e.curve.w[2], 2.0);

    ASSERT_EQ(gshgeom_path_element(p, 5, &e), GSHGEOM_OK);
    EXPECT_EQ(e.verb, GSHGEOM_VERB_RATIONAL_CUBIC);
    EXPECT_EQ(e.curve.kind, GSHGEOM_CURVE_RATIONAL_CUBIC);
    EXPECT_EQ(e.curve.p[0].x, 17.0);
    EXPECT_EQ(e.curve.p[3].y, 24.1);
    EXPECT_EQ(e.curve.w[0], 1.0);
    EXPECT_EQ(e.curve.w[1], 2.0);
    EXPECT_EQ(e.curve.w[2], 0.5);
    EXPECT_EQ(e.curve.w[3], 4.0);

    // The close is the line back to the start of the subpath.
    ASSERT_EQ(gshgeom_path_element(p, 6, &e), GSHGEOM_OK);
    EXPECT_EQ(e.verb, GSHGEOM_VERB_CLOSE);
    EXPECT_EQ(e.curve.kind, GSHGEOM_CURVE_LINE);
    EXPECT_EQ(e.curve.p[0].x, 23.0);
    EXPECT_EQ(e.curve.p[0].y, 24.1);
    EXPECT_EQ(e.curve.p[1].x, 1.5);
    EXPECT_EQ(e.curve.p[1].y, -2.25);

    // Every element read back is a valid curve.
    for (size_t i = 0; i < 7; i++) {
      ASSERT_EQ(gshgeom_path_element(p, i, &e), GSHGEOM_OK);
      EXPECT_EQ(gshgeom_curve_check(&e.curve), GSHGEOM_OK) << i;
    }
    gshgeom_path_destroy(p);
  }
}

TEST(Path, ElementReadbackChecksItsArguments) {
  GSHGEOM_Path * p = fixture_path(2);
  GSHGEOM_Element e;
  std::memset(&e, 0x77, sizeof e);
  GSHGEOM_Element before = e;
  EXPECT_EQ(gshgeom_path_element(p, 7, &e), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_element(p, (size_t)-1, &e), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_element(nullptr, 0, &e), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_element(p, 0, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(std::memcmp(&e, &before, sizeof e), 0);
  EXPECT_EQ(gshgeom_path_element(p, 6, &e), GSHGEOM_OK);
  gshgeom_path_destroy(p);
}

TEST(Path, TheSubpathRules) {
  for (int dim = 2; dim <= 3; dim++) {
    GSHGEOM_Path * p = make_path(dim);
    double w3[3] = {1, 1, 1};
    double w4[4] = {1, 1, 1, 1};
    // Nothing draws, and nothing closes, before a move.
    EXPECT_EQ(gshgeom_path_line_to(p, P(1, 1)), GSHGEOM_ERR_INVALID);
    EXPECT_EQ(gshgeom_path_quadratic_to(p, P(1, 1), P(2, 2)), GSHGEOM_ERR_INVALID);
    EXPECT_EQ(gshgeom_path_cubic_to(p, P(1, 1), P(2, 2), P(3, 3)),
        GSHGEOM_ERR_INVALID);
    EXPECT_EQ(gshgeom_path_rational_quadratic_to(p, P(1, 1), P(2, 2), w3),
        GSHGEOM_ERR_INVALID);
    EXPECT_EQ(gshgeom_path_rational_cubic_to(p, P(1, 1), P(2, 2), P(3, 3), w4),
        GSHGEOM_ERR_INVALID);
    EXPECT_EQ(gshgeom_path_close(p), GSHGEOM_ERR_INVALID);
    EXPECT_EQ(gshgeom_path_element_count(p), 0u);

    // Two moves are two elements.
    EXPECT_EQ(gshgeom_path_move_to(p, P(1, 1)), GSHGEOM_OK);
    EXPECT_EQ(gshgeom_path_move_to(p, P(2, 2)), GSHGEOM_OK);
    EXPECT_EQ(gshgeom_path_element_count(p), 2u);
    EXPECT_EQ(gshgeom_path_line_to(p, P(5, 2)), GSHGEOM_OK);
    GSHGEOM_Element e;
    ASSERT_EQ(gshgeom_path_element(p, 2, &e), GSHGEOM_OK);
    EXPECT_EQ(e.curve.p[0].x, 2.0) << "starts at the second move";

    // After a close the current point is the start, and a drawing verb
    // begins a new subpath there.
    EXPECT_EQ(gshgeom_path_close(p), GSHGEOM_OK);
    EXPECT_EQ(gshgeom_path_line_to(p, P(9, 9)), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_path_element(p, 4, &e), GSHGEOM_OK);
    EXPECT_EQ(e.verb, GSHGEOM_VERB_LINE);
    EXPECT_EQ(e.curve.p[0].x, 2.0);
    EXPECT_EQ(e.curve.p[0].y, 2.0);
    // The new subpath closes to the same start.
    EXPECT_EQ(gshgeom_path_close(p), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_path_element(p, 5, &e), GSHGEOM_OK);
    EXPECT_EQ(e.curve.p[0].x, 9.0);
    EXPECT_EQ(e.curve.p[1].x, 2.0);
    // A move starts a different subpath, and a close goes back to it.
    EXPECT_EQ(gshgeom_path_move_to(p, P(40, 40)), GSHGEOM_OK);
    EXPECT_EQ(gshgeom_path_line_to(p, P(41, 40)), GSHGEOM_OK);
    EXPECT_EQ(gshgeom_path_close(p), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_path_element(p, 8, &e), GSHGEOM_OK);
    EXPECT_EQ(e.curve.p[1].x, 40.0);
    // A close right after a move is a line of length zero.
    EXPECT_EQ(gshgeom_path_move_to(p, P(7, 7)), GSHGEOM_OK);
    EXPECT_EQ(gshgeom_path_close(p), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_path_element(p, 10, &e), GSHGEOM_OK);
    EXPECT_EQ(e.curve.p[0].x, 7.0);
    EXPECT_EQ(e.curve.p[1].x, 7.0);
    gshgeom_path_destroy(p);
  }
}

TEST(Path, EveryVerbRefusesWhatIsNotANumberAndLeavesThePathAlone) {
  const double bad[] = {kNaN, kInf, -kInf};
  for (int dim = 2; dim <= 3; dim++) {
    GSHGEOM_Path * p = make_path(dim);
    ASSERT_EQ(gshgeom_path_move_to(p, P(1, 2, 3)), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_path_line_to(p, P(4, 5, 6)), GSHGEOM_OK);
    Snapshot before(p);
    for (double v : bad) {
      for (int a = 0; a < dim; a++) {
        auto q = [&](double x, double y, double z) {
          GSHGEOM_Point r = P(x, y, z);
          (a == 0 ? r.x : a == 1 ? r.y : r.z) = v;
          return r;
        };
        double w3[3] = {1, 2, 3}, w4[4] = {1, 2, 3, 4};
        EXPECT_EQ(gshgeom_path_move_to(p, q(0, 0, 0)), GSHGEOM_ERR_INVALID);
        EXPECT_EQ(gshgeom_path_line_to(p, q(0, 0, 0)), GSHGEOM_ERR_INVALID);
        EXPECT_EQ(gshgeom_path_quadratic_to(p, q(0, 0, 0), P(1, 1)),
            GSHGEOM_ERR_INVALID);
        EXPECT_EQ(gshgeom_path_quadratic_to(p, P(1, 1), q(0, 0, 0)),
            GSHGEOM_ERR_INVALID);
        EXPECT_EQ(gshgeom_path_cubic_to(p, q(0, 0, 0), P(1, 1), P(2, 2)),
            GSHGEOM_ERR_INVALID);
        EXPECT_EQ(gshgeom_path_cubic_to(p, P(1, 1), q(0, 0, 0), P(2, 2)),
            GSHGEOM_ERR_INVALID);
        EXPECT_EQ(gshgeom_path_cubic_to(p, P(1, 1), P(2, 2), q(0, 0, 0)),
            GSHGEOM_ERR_INVALID);
        EXPECT_EQ(gshgeom_path_rational_quadratic_to(
                      p, q(0, 0, 0), P(1, 1), w3),
            GSHGEOM_ERR_INVALID);
        EXPECT_EQ(gshgeom_path_rational_quadratic_to(
                      p, P(1, 1), q(0, 0, 0), w3),
            GSHGEOM_ERR_INVALID);
        EXPECT_EQ(gshgeom_path_rational_cubic_to(
                      p, q(0, 0, 0), P(1, 1), P(2, 2), w4),
            GSHGEOM_ERR_INVALID);
        EXPECT_EQ(gshgeom_path_rational_cubic_to(
                      p, P(1, 1), P(2, 2), q(0, 0, 0), w4),
            GSHGEOM_ERR_INVALID);
      }
      for (int slot = 0; slot < 4; slot++) {
        double w3[3] = {1, 2, 3}, w4[4] = {1, 2, 3, 4};
        w3[slot % 3] = v;
        w4[slot] = v;
        EXPECT_EQ(gshgeom_path_rational_quadratic_to(p, P(1, 1), P(2, 2), w3),
            GSHGEOM_ERR_INVALID);
        EXPECT_EQ(gshgeom_path_rational_cubic_to(
                      p, P(1, 1), P(2, 2), P(3, 3), w4),
            GSHGEOM_ERR_INVALID);
      }
    }
    for (double v : {0.0, -1.0}) {
      double w3[3] = {1, v, 3};
      EXPECT_EQ(gshgeom_path_rational_quadratic_to(p, P(1, 1), P(2, 2), w3),
          GSHGEOM_ERR_INVALID);
    }
    EXPECT_EQ(gshgeom_path_rational_quadratic_to(p, P(1, 1), P(2, 2), nullptr),
        GSHGEOM_ERR_INVALID);
    EXPECT_EQ(gshgeom_path_rational_cubic_to(
                  p, P(1, 1), P(2, 2), P(3, 3), nullptr),
        GSHGEOM_ERR_INVALID);
    // Too big is a limit, not an invalid argument, and equally harmless.
    double big = 2 * GSHGEOM_COORD_LIMIT;
    EXPECT_EQ(gshgeom_path_line_to(p, P(big, 0)), GSHGEOM_ERR_LIMIT);
    EXPECT_EQ(gshgeom_path_move_to(p, P(0, -big)), GSHGEOM_ERR_LIMIT);
    double wbig[3] = {1, 2 * GSHGEOM_WEIGHT_MAX, 1};
    EXPECT_EQ(gshgeom_path_rational_quadratic_to(p, P(1, 1), P(2, 2), wbig),
        GSHGEOM_ERR_LIMIT);
    EXPECT_TRUE(Snapshot(p) == before);
    // And the path is still usable.
    EXPECT_EQ(gshgeom_path_line_to(p, P(7, 8)), GSHGEOM_OK);
    EXPECT_EQ(gshgeom_path_element_count(p), 3u);
    gshgeom_path_destroy(p);
  }
}

TEST(Path, TheIgnoredZOfA2DPathIsNotLookedAt) {
  GSHGEOM_Path * a = make_path(2);
  GSHGEOM_Path * b = make_path(2);
  EXPECT_EQ(gshgeom_path_move_to(a, P(1, 2, kNaN)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_line_to(a, P(3, 4, kInf)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_move_to(b, P(1, 2, 5)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_line_to(b, P(3, 4, -9)), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_hash(a), gshgeom_path_hash(b));
  EXPECT_EQ(dump_of(a), dump_of(b));
  gshgeom_path_destroy(a);
  gshgeom_path_destroy(b);
}

TEST(Path, NullPathIsInvalidEverywhere) {
  double w3[3] = {1, 1, 1}, w4[4] = {1, 1, 1, 1};
  EXPECT_EQ(gshgeom_path_move_to(nullptr, P(0, 0)), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_line_to(nullptr, P(0, 0)), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_quadratic_to(nullptr, P(0, 0), P(1, 1)),
      GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_cubic_to(nullptr, P(0, 0), P(1, 1), P(2, 2)),
      GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_rational_quadratic_to(nullptr, P(0, 0), P(1, 1), w3),
      GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_rational_cubic_to(
                nullptr, P(0, 0), P(1, 1), P(2, 2), w4),
      GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_close(nullptr), GSHGEOM_ERR_INVALID);
  GSHGEOM_Bounds b;
  EXPECT_EQ(gshgeom_path_bounds(nullptr, &b), GSHGEOM_ERR_INVALID);
  GSHGEOM_Path * p = make_path(2);
  EXPECT_EQ(gshgeom_path_bounds(p, nullptr), GSHGEOM_ERR_INVALID);
  gshgeom_path_destroy(p);
  EXPECT_EQ(gshgeom_path_hash(nullptr), GSHGEOM_FNV_OFFSET_BASIS);
  EXPECT_EQ(gshgeom_path_dump(nullptr, stdout), GSHGEOM_ERR_INVALID);
}

// ---------------------------------------------------------------------------
// Limits
// ---------------------------------------------------------------------------

TEST(PathLimits, TheBoundIsAssertedNotTheStatus) {
  for (size_t cap : {(size_t)1, (size_t)2, (size_t)5, (size_t)64, (size_t)65}) {
    GSHGEOM_Limits l;
    gshgeom_limits_default(&l);
    l.max_elements = cap;
    GSHGEOM_Path * p = make_path(2, nullptr, &l);
    ASSERT_EQ(gshgeom_path_move_to(p, P(0, 0)), GSHGEOM_OK);
    for (size_t i = 1; i < cap; i++) {
      ASSERT_EQ(gshgeom_path_line_to(p, P((double)i, 0)), GSHGEOM_OK) << i;
    }
    EXPECT_EQ(gshgeom_path_element_count(p), cap) << "element N is accepted";
    Snapshot before(p);
    EXPECT_EQ(gshgeom_path_line_to(p, P(-1, -1)), GSHGEOM_ERR_LIMIT)
        << "and N + 1 is not";
    EXPECT_EQ(gshgeom_path_move_to(p, P(-1, -1)), GSHGEOM_ERR_LIMIT);
    EXPECT_EQ(gshgeom_path_close(p), GSHGEOM_ERR_LIMIT);
    double w[3] = {1, 1, 1};
    EXPECT_EQ(gshgeom_path_rational_quadratic_to(p, P(1, 1), P(2, 2), w),
        GSHGEOM_ERR_LIMIT);
    EXPECT_EQ(gshgeom_path_element_count(p), cap) << "no truncation either way";
    EXPECT_TRUE(Snapshot(p) == before);
    gshgeom_path_destroy(p);
  }
}

TEST(PathLimits, ZeroMeansDefaultAndNullMeansDefaults) {
  GSHGEOM_Limits zero;
  std::memset(&zero, 0, sizeof zero);
  GSHGEOM_Path * p = make_path(2, nullptr, &zero);
  for (int i = 0; i < 1000; i++) {
    ASSERT_EQ(gshgeom_path_move_to(p, P(i, i)), GSHGEOM_OK);
  }
  gshgeom_path_destroy(p);
}

TEST(PathLimits, ALimitBeatsAnInvalidArgumentOnlyWhenTheArgumentIsValid) {
  GSHGEOM_Limits l;
  gshgeom_limits_default(&l);
  l.max_elements = 1;
  GSHGEOM_Path * p = make_path(2, nullptr, &l);
  ASSERT_EQ(gshgeom_path_move_to(p, P(0, 0)), GSHGEOM_OK);
  // A NaN is refused as what it is, and is not hidden by being over the limit.
  EXPECT_EQ(gshgeom_path_line_to(p, P(kNaN, 0)), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_line_to(p, P(1, 0)), GSHGEOM_ERR_LIMIT);
  gshgeom_path_destroy(p);
}

// ---------------------------------------------------------------------------
// Bounds
// ---------------------------------------------------------------------------

TEST(PathBounds, EmptyWhenNothingIsDrawn) {
  for (int dim = 2; dim <= 3; dim++) {
    GSHGEOM_Path * p = make_path(dim);
    GSHGEOM_Bounds b;
    std::memset(&b, 0x77, sizeof b);
    ASSERT_EQ(gshgeom_path_bounds(p, &b), GSHGEOM_OK);
    EXPECT_EQ(b.empty, 1);
    EXPECT_EQ(b.min.x, 0.0);
    EXPECT_EQ(b.max.z, 0.0);
    ASSERT_EQ(gshgeom_path_move_to(p, P(5, 5, 5)), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_path_move_to(p, P(9, 9, 9)), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_path_close(p), GSHGEOM_OK);
    ASSERT_EQ(gshgeom_path_bounds(p, &b), GSHGEOM_OK);
    EXPECT_EQ(b.empty, 1) << "moves and a close of nothing add nothing";
    gshgeom_path_destroy(p);
  }
}

TEST(PathBounds, AMoveThatStartsNoSegmentAddsNothing) {
  GSHGEOM_Path * p = make_path(2);
  ASSERT_EQ(gshgeom_path_move_to(p, P(0, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(p, P(2, 1)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_move_to(p, P(100, -100)), GSHGEOM_OK);
  GSHGEOM_Bounds b;
  ASSERT_EQ(gshgeom_path_bounds(p, &b), GSHGEOM_OK);
  EXPECT_EQ(b.empty, 0);
  EXPECT_EQ(b.min.x, 0.0);
  EXPECT_EQ(b.max.x, 2.0);
  EXPECT_EQ(b.min.y, 0.0);
  EXPECT_EQ(b.max.y, 1.0);
  gshgeom_path_destroy(p);
}

TEST(PathBounds, UnionOfSegmentsOfAClosedShape2D) {
  GSHGEOM_Path * p = make_path(2);
  ASSERT_EQ(gshgeom_path_move_to(p, P(0, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(p, P(10, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_quadratic_to(p, P(15, 10), P(20, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_close(p), GSHGEOM_OK);
  GSHGEOM_Bounds b;
  ASSERT_EQ(gshgeom_path_bounds(p, &b), GSHGEOM_OK);
  EXPECT_EQ(b.min.x, 0.0);
  EXPECT_EQ(b.max.x, 20.0);
  EXPECT_EQ(b.min.y, 0.0);
  EXPECT_NEAR(b.max.y, 5.0, 1e-15);
  EXPECT_EQ(b.min.z, 0.0);
  EXPECT_EQ(b.max.z, 0.0);
  gshgeom_path_destroy(p);
}

TEST(PathBounds, ARationalCircleOfFourArcsIsTheBoxOfTheCircle) {
  // The unit circle about (3, -2): four rational quadratic quarter arcs.
  const double w = gmath_sqrt(0.5);
  double ws[3] = {1.0, w, 1.0};
  GSHGEOM_Path * p = make_path(2);
  ASSERT_EQ(gshgeom_path_move_to(p, P(4, -2)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_rational_quadratic_to(p, P(4, -1), P(3, -1), ws),
      GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_rational_quadratic_to(p, P(2, -1), P(2, -2), ws),
      GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_rational_quadratic_to(p, P(2, -3), P(3, -3), ws),
      GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_rational_quadratic_to(p, P(4, -3), P(4, -2), ws),
      GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_close(p), GSHGEOM_OK);
  GSHGEOM_Bounds b;
  ASSERT_EQ(gshgeom_path_bounds(p, &b), GSHGEOM_OK);
  EXPECT_EQ(b.min.x, 2.0);
  EXPECT_EQ(b.max.x, 4.0);
  EXPECT_EQ(b.min.y, -3.0);
  EXPECT_EQ(b.max.y, -1.0);
  gshgeom_path_destroy(p);
}

TEST(PathBounds, TheBoxOfThePathIsTheUnionOfTheBoxesOfItsElements) {
  for (int dim = 2; dim <= 3; dim++) {
    GSHGEOM_Path * p = fixture_path(dim);
    GSHGEOM_Bounds want;
    want.empty = 1;
    for (size_t i = 0; i < gshgeom_path_element_count(p); i++) {
      GSHGEOM_Element e;
      ASSERT_EQ(gshgeom_path_element(p, i, &e), GSHGEOM_OK);
      if (e.verb == GSHGEOM_VERB_MOVE || e.verb == GSHGEOM_VERB_CLOSE) {
        continue;
      }
      GSHGEOM_Bounds b;
      ASSERT_EQ(gshgeom_curve_bounds(&e.curve, &b), GSHGEOM_OK);
      if (want.empty) {
        want = b;
      } else {
        want.min.x = std::min(want.min.x, b.min.x);
        want.min.y = std::min(want.min.y, b.min.y);
        want.min.z = std::min(want.min.z, b.min.z);
        want.max.x = std::max(want.max.x, b.max.x);
        want.max.y = std::max(want.max.y, b.max.y);
        want.max.z = std::max(want.max.z, b.max.z);
      }
    }
    GSHGEOM_Bounds got;
    ASSERT_EQ(gshgeom_path_bounds(p, &got), GSHGEOM_OK);
    EXPECT_EQ(got.empty, 0);
    EXPECT_EQ(std::memcmp(&got.min, &want.min, sizeof got.min), 0);
    EXPECT_EQ(std::memcmp(&got.max, &want.max, sizeof got.max), 0);
    EXPECT_GT(got.max.x, 20.0);
    if (dim == 3) {
      EXPECT_GT(got.max.z, got.min.z);
    }
    gshgeom_path_destroy(p);
  }
}

TEST(PathBounds, TwoSubpathsAreOneBox) {
  GSHGEOM_Path * p = make_path(3);
  ASSERT_EQ(gshgeom_path_move_to(p, P(0, 0, 0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_cubic_to(p, P(0, 1, 2), P(1, 1, 2), P(1, 0, 0)),
      GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_move_to(p, P(-5, -6, -7)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(p, P(-4, -5, -6)), GSHGEOM_OK);
  GSHGEOM_Bounds b;
  ASSERT_EQ(gshgeom_path_bounds(p, &b), GSHGEOM_OK);
  EXPECT_EQ(b.min.x, -5.0);
  EXPECT_EQ(b.min.y, -6.0);
  EXPECT_EQ(b.min.z, -7.0);
  EXPECT_EQ(b.max.x, 1.0);
  EXPECT_NEAR(b.max.y, 0.75, 1e-15);
  EXPECT_NEAR(b.max.z, 1.5, 1e-15);
  gshgeom_path_destroy(p);
}

// ---------------------------------------------------------------------------
// Hash and dump
// ---------------------------------------------------------------------------

TEST(PathHash, KnownAnswers) {
  GSHGEOM_Path * p2 = fixture_path(2);
  GSHGEOM_Path * p3 = fixture_path(3);
  EXPECT_EQ(gshgeom_path_hash(p2), kPathHash2);
  EXPECT_EQ(gshgeom_path_hash(p3), kPathHash3);
  gshgeom_path_destroy(p2);
  gshgeom_path_destroy(p3);
  GSHGEOM_Path * e2 = make_path(2);
  GSHGEOM_Path * e3 = make_path(3);
  EXPECT_EQ(gshgeom_path_hash(e2), kPathHashEmpty2);
  EXPECT_EQ(gshgeom_path_hash(e3), kPathHashEmpty3);
  gshgeom_path_destroy(e2);
  gshgeom_path_destroy(e3);
}

TEST(PathHash, TheSpellingIsWrittenOutHere) {
  // The same hash computed field by field with the public hasher, in the
  // order path.h states, for the empty and for a one-segment path.
  GSHGEOM_Path * p = make_path(3);
  ASSERT_EQ(gshgeom_path_move_to(p, P(1, 2, 3)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(p, P(4, 5, 6)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_close(p), GSHGEOM_OK);
  GSHGEOM_Hash h;
  gshgeom_hash_init(&h);
  gshgeom_hash_u32(&h, 1);
  gshgeom_hash_u32(&h, 3);
  gshgeom_hash_u64(&h, 3);
  gshgeom_hash_u8(&h, 0);
  for (double v : {1.0, 2.0, 3.0}) gshgeom_hash_f64(&h, v);
  gshgeom_hash_u8(&h, 1);
  for (double v : {4.0, 5.0, 6.0}) gshgeom_hash_f64(&h, v);
  gshgeom_hash_u8(&h, 6);
  EXPECT_EQ(gshgeom_path_hash(p), gshgeom_hash_value(&h));
  gshgeom_path_destroy(p);
}

TEST(PathHash, EveryPartOfAPathChangesIt) {
  auto build = [](int dim, double x, double zz, bool close, int verb,
                   double w1) {
    GSHGEOM_Path * p = make_path(dim);
    EXPECT_EQ(gshgeom_path_move_to(p, P(0, 0, 0)), GSHGEOM_OK);
    double w[3] = {1.0, w1, 1.0};
    switch (verb) {
    case 0:
      EXPECT_EQ(gshgeom_path_line_to(p, P(x, 1, zz)), GSHGEOM_OK);
      break;
    case 1:
      EXPECT_EQ(gshgeom_path_quadratic_to(p, P(0.5, 0.5, 0.5), P(x, 1, zz)),
          GSHGEOM_OK);
      break;
    default:
      EXPECT_EQ(gshgeom_path_rational_quadratic_to(
                    p, P(0.5, 0.5, 0.5), P(x, 1, zz), w),
          GSHGEOM_OK);
      break;
    }
    if (close) {
      EXPECT_EQ(gshgeom_path_close(p), GSHGEOM_OK);
    }
    uint64_t h = gshgeom_path_hash(p);
    gshgeom_path_destroy(p);
    return h;
  };
  uint64_t base = build(3, 2.0, 3.0, true, 2, 0.5);
  EXPECT_EQ(base, build(3, 2.0, 3.0, true, 2, 0.5));
  EXPECT_NE(base, build(2, 2.0, 3.0, true, 2, 0.5)) << "the dimension";
  EXPECT_NE(base, build(3, 2.5, 3.0, true, 2, 0.5)) << "a coordinate";
  EXPECT_NE(base, build(3, 2.0, 3.5, true, 2, 0.5)) << "z in 3D";
  EXPECT_EQ(build(2, 2.0, 3.0, true, 2, 0.5), build(2, 2.0, -8.0, true, 2, 0.5))
      << "but z in 2D is ignored";
  EXPECT_NE(base, build(3, 2.0, 3.0, false, 2, 0.5)) << "the close";
  EXPECT_NE(base, build(3, 2.0, 3.0, true, 1, 0.5)) << "the verb";
  EXPECT_NE(base, build(3, 2.0, 3.0, true, 2, 0.75)) << "a weight";
  EXPECT_NE(build(3, 0.0, 0.0, true, 0, 1.0), build(3, -0.0, 0.0, true, 0, 1.0))
      << "-0.0 is not 0.0";
}

TEST(PathDump, KnownText) {
  GSHGEOM_Path * p2 = fixture_path(2);
  GSHGEOM_Path * p3 = fixture_path(3);
  EXPECT_EQ(dump_of(p2), std::string(kPathDump2));
  EXPECT_EQ(dump_of(p3), std::string(kPathDump3));
  gshgeom_path_destroy(p2);
  gshgeom_path_destroy(p3);
}

TEST(PathDump, ZeroIsSixteenDigitsToo) {
  GSHGEOM_Path * p = make_path(2);
  ASSERT_EQ(gshgeom_path_move_to(p, P(0, -0.0)), GSHGEOM_OK);
  ASSERT_EQ(gshgeom_path_line_to(p, P(1, 0)), GSHGEOM_OK);
  std::string d = dump_of(p);
  EXPECT_NE(d.find("\nM 0000000000000000 8000000000000000\n"), std::string::npos)
      << d;
  EXPECT_NE(d.find("\nL 3ff0000000000000 0000000000000000\n"), std::string::npos)
      << d;
  gshgeom_path_destroy(p);
}

TEST(PathDump, AFailingSinkIsReported) {
  GSHGEOM_Path * p = fixture_path(2);
  FILE * full = std::fopen("/dev/full", "w");
  if (!full) {
    gshgeom_path_destroy(p);
    GTEST_SKIP() << "no /dev/full on this host";
  }
  EXPECT_EQ(gshgeom_path_dump(p, full), GSHGEOM_ERR_IO);
  std::fclose(full);
  // A stream open only for reading refuses the first write.
  FILE * ro = std::fopen("/dev/null", "r");
  ASSERT_NE(ro, nullptr);
  EXPECT_EQ(gshgeom_path_dump(p, ro), GSHGEOM_ERR_IO);
  std::fclose(ro);
  EXPECT_EQ(gshgeom_path_dump(p, nullptr), GSHGEOM_ERR_INVALID);
  gshgeom_path_destroy(p);
}

TEST(PathDump, DoesNotDependOnTheLocale) {
  // The dump prints bit patterns as integers; a comma locale changes nothing.
  GSHGEOM_Path * p = fixture_path(2);
  std::string before = dump_of(p);
  const char * old = std::setlocale(LC_ALL, nullptr);
  std::string saved = old ? old : "C";
  for (const char * name : {"de_DE.UTF-8", "fr_FR.UTF-8", "de_DE", "C.UTF-8"}) {
    if (std::setlocale(LC_ALL, name)) {
      EXPECT_EQ(dump_of(p), before) << name;
    }
  }
  std::setlocale(LC_ALL, saved.c_str());
  gshgeom_path_destroy(p);
}

// ---------------------------------------------------------------------------
// Allocation failure
// ---------------------------------------------------------------------------

namespace {

struct Tracker {
  long fail_from = -1;    // the call number that first fails; -1 never
  bool persistent = false;  // every call from fail_from on fails
  long calls = 0;
  bool failed = false;
  std::set<void *> live;
  bool fails() {
    long n = calls++;
    if (fail_from >= 0 && (persistent ? n >= fail_from : n == fail_from)) {
      failed = true;
      return true;
    }
    return false;
  }
};

void * t_malloc(void * ctx, size_t n) {
  Tracker * t = (Tracker *)ctx;
  if (t->fails()) return nullptr;
  void * p = std::malloc(n ? n : 1);
  if (p) {
    t->live.insert(p);
  }
  return p;
}
void * t_calloc(void * ctx, size_t c, size_t n) {
  Tracker * t = (Tracker *)ctx;
  if (t->fails()) return nullptr;
  void * p = std::calloc(c ? c : 1, n ? n : 1);
  if (p) {
    t->live.insert(p);
  }
  return p;
}
void * t_realloc(void * ctx, void * old, size_t n) {
  Tracker * t = (Tracker *)ctx;
  if (t->fails()) return nullptr;
  void * p = std::realloc(old, n ? n : 1);
  if (p) {
    t->live.erase(old);
    t->live.insert(p);
  }
  return p;
}
void t_free(void * ctx, void * p) {
  Tracker * t = (Tracker *)ctx;
  if (p) {
    EXPECT_EQ(t->live.erase(p), 1u) << "freed what was not allocated here";
  }
  std::free(p);
}

GSHGEOM_Allocator allocator_for(Tracker * t) {
  GSHGEOM_Allocator a;
  a.ctx = t;
  a.malloc_fn = t_malloc;
  a.calloc_fn = t_calloc;
  a.realloc_fn = t_realloc;
  a.free_fn = t_free;
  return a;
}

}  // namespace

TEST(PathAllocation, ARefusedAllocationIsOomAndChangesNothing) {
  // How many allocator calls does building the fixture make?
  Tracker count;
  {
    GSHGEOM_Allocator a = allocator_for(&count);
    GSHGEOM_Path * p = make_path(3, &a);
    for (const Step & s : fixture_steps()) ASSERT_EQ(s(p), GSHGEOM_OK);
    gshgeom_path_destroy(p);
    EXPECT_TRUE(count.live.empty());
  }
  ASSERT_GT(count.calls, 3);

  for (int persistent = 0; persistent < 2; persistent++) {
    for (long k = 0; k < count.calls; k++) {
      Tracker t;
      t.fail_from = k;
      t.persistent = persistent != 0;
      GSHGEOM_Allocator a = allocator_for(&t);
      GSHGEOM_Path * p = nullptr;
      GSHGEOM_Result r = gshgeom_path_create(&a, nullptr, 3, &p);
      if (r != GSHGEOM_OK) {
        EXPECT_EQ(r, GSHGEOM_ERR_OOM);
        EXPECT_EQ(p, nullptr) << "written only on success";
        EXPECT_TRUE(t.live.empty());
        continue;
      }
      std::vector<Step> steps = fixture_steps();
      size_t done = 0;
      for (const Step & s : steps) {
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
        EXPECT_EQ(gshgeom_path_hash(p), kPathHash3) << "k=" << k;
      }
      gshgeom_path_destroy(p);
      EXPECT_TRUE(t.live.empty()) << "k=" << k << " leaked";
    }
  }
}

TEST(PathAllocation, GrowthKeepsEverythingAcrossManyReallocations) {
  GSHGEOM_Path * p = make_path(3);
  const int n = 5000;
  ASSERT_EQ(gshgeom_path_move_to(p, P(0, 0, 0)), GSHGEOM_OK);
  double w[4] = {1, 2, 3, 4};
  for (int i = 1; i < n; i++) {
    double x = i;
    switch (i % 5) {
    case 0:
      ASSERT_EQ(gshgeom_path_line_to(p, P(x, 1, 2)), GSHGEOM_OK);
      break;
    case 1:
      ASSERT_EQ(gshgeom_path_quadratic_to(p, P(x, 1, 2), P(x + 1, 2, 3)),
          GSHGEOM_OK);
      break;
    case 2:
      ASSERT_EQ(gshgeom_path_cubic_to(p, P(x, 1, 2), P(x, 2, 3), P(x + 1, 3, 4)),
          GSHGEOM_OK);
      break;
    case 3:
      ASSERT_EQ(gshgeom_path_rational_quadratic_to(
                    p, P(x, 1, 2), P(x + 1, 2, 3), w),
          GSHGEOM_OK);
      break;
    default:
      ASSERT_EQ(gshgeom_path_rational_cubic_to(
                    p, P(x, 1, 2), P(x, 2, 3), P(x + 1, 3, 4), w),
          GSHGEOM_OK);
      break;
    }
  }
  ASSERT_EQ(gshgeom_path_element_count(p), (size_t)n);
  // Each element reads back with its own weights.
  for (int i = 1; i < n; i++) {
    GSHGEOM_Element e;
    ASSERT_EQ(gshgeom_path_element(p, (size_t)i, &e), GSHGEOM_OK);
    int kind = i % 5;
    EXPECT_EQ((int)e.verb, kind == 0 ? 1 : kind + 1);
    if (kind == 3) {
      EXPECT_EQ(e.curve.w[2], 3.0);
    }
    if (kind == 4) {
      EXPECT_EQ(e.curve.w[3], 4.0);
    }
  }
  gshgeom_path_destroy(p);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
