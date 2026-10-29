/**
 * @file
 *
 * The vertices of a path and the directions at them (vertices.h): the positions, the incoming and outgoing
 * directions, the closed subpath's round reading, the rule for segments without extent, the marker direction and
 * its bisector, and the flags; against a second reading of the definition (the cyclic scan, written out the long
 * way in long double) over random paths with deliberate repeats.
 *
 * Copyright 2026 by Corey Pennycuff
 */

#include "test_helpers.h"

#include <ghoti.io/math/math.h>

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>

namespace {

enum Verb : char { M = 'M', L = 'L', Q = 'Q', C = 'C', Z = 'Z' };

struct Op {
  char verb;
  double v[6];
};

Op op(char verb, double a = 0, double b = 0, double c = 0, double d = 0, double e = 0, double f = 0) {
  return Op{verb, {a, b, c, d, e, f}};
}

GSHGEOM_Point P(double x, double y) {
  return GSHGEOM_Point{x, y, 0.0};
}

GSHGEOM_Path * build(const std::vector<Op> & ops, int dim = 2) {
  GSHGEOM_Path * p = nullptr;
  EXPECT_EQ(gshgeom_path_create(nullptr, nullptr, dim, &p), GSHGEOM_OK);
  for (const Op & o : ops) {
    GSHGEOM_Result r = GSHGEOM_OK;
    switch (o.verb) {
    case M: r = gshgeom_path_move_to(p, P(o.v[0], o.v[1])); break;
    case L: r = gshgeom_path_line_to(p, P(o.v[0], o.v[1])); break;
    case Q: r = gshgeom_path_quadratic_to(p, P(o.v[0], o.v[1]), P(o.v[2], o.v[3])); break;
    case C: r = gshgeom_path_cubic_to(p, P(o.v[0], o.v[1]), P(o.v[2], o.v[3]), P(o.v[4], o.v[5])); break;
    case Z: r = gshgeom_path_close(p); break;
    }
    EXPECT_EQ(r, GSHGEOM_OK);
  }
  return p;
}

GSHGEOM_Result collect_one(const GSHGEOM_Vertex * v, void * ctx) {
  static_cast<std::vector<GSHGEOM_Vertex> *>(ctx)->push_back(*v);
  return GSHGEOM_OK;
}

std::vector<GSHGEOM_Vertex> vertices_of(const std::vector<Op> & ops) {
  GSHGEOM_Path * p = build(ops);
  std::vector<GSHGEOM_Vertex> out;
  EXPECT_EQ(gshgeom_path_vertices(p, collect_one, &out), GSHGEOM_OK);
  size_t count = 0;
  EXPECT_EQ(gshgeom_path_vertex_count(p, &count), GSHGEOM_OK);
  EXPECT_EQ(count, out.size());
  gshgeom_path_destroy(p);
  return out;
}

const double kRoot2 = 1.0 / std::sqrt(2.0); // (1, 1) / |(1, 1)| in the library's own operations: sqrt(1 + 1), then x / n

void expect_pt(const GSHGEOM_Point & got, double x, double y, const char * what, size_t at) {
  EXPECT_EQ(got.x, x) << what << " x, vertex " << at;
  EXPECT_EQ(got.y, y) << what << " y, vertex " << at;
  EXPECT_EQ(got.z, 0.0) << what << " z, vertex " << at;
}

constexpr uint32_t kIn = GSHGEOM_VERTEX_HAS_IN, kOut = GSHGEOM_VERTEX_HAS_OUT, kPS = GSHGEOM_VERTEX_PATH_START,
                   kPE = GSHGEOM_VERTEX_PATH_END, kSS = GSHGEOM_VERTEX_SUBPATH_START, kSE = GSHGEOM_VERTEX_SUBPATH_END,
                   kCl = GSHGEOM_VERTEX_CLOSED;

} // namespace

TEST(Vertices, AnOpenPolylineHasItsEndsAndOneCorner) {
  auto v = vertices_of({op(M, 0, 0), op(L, 10, 0), op(L, 10, 10)});
  ASSERT_EQ(v.size(), 3u);
  EXPECT_EQ(v[0].flags, kOut | kPS | kSS);
  expect_pt(v[0].point, 0, 0, "point", 0);
  expect_pt(v[0].in, 0, 0, "in", 0);
  expect_pt(v[0].out, 1, 0, "out", 0);
  expect_pt(v[0].direction, 1, 0, "direction", 0);
  EXPECT_EQ(v[1].flags, kIn | kOut);
  expect_pt(v[1].point, 10, 0, "point", 1);
  expect_pt(v[1].in, 1, 0, "in", 1);
  expect_pt(v[1].out, 0, 1, "out", 1);
  expect_pt(v[1].direction, kRoot2, kRoot2, "direction", 1);
  EXPECT_EQ(v[2].flags, kIn | kPE | kSE);
  expect_pt(v[2].in, 0, 1, "in", 2);
  expect_pt(v[2].out, 0, 0, "out", 2);
  expect_pt(v[2].direction, 0, 1, "direction", 2);
  for (size_t i = 0; i < v.size(); i++) {
    EXPECT_EQ(v[i].element, i);
    EXPECT_EQ(v[i].subpath, 0u);
  }
}

TEST(Vertices, AClosedSubpathIsReadRoundAtItsStartAndItsClose) {
  auto v = vertices_of({op(M, 0, 0), op(L, 1, 0), op(L, 1, 1), op(L, 0, 1), op(Z)});
  ASSERT_EQ(v.size(), 5u);
  // The move and the close are one point and one corner: the closing line comes in going -y, the first segment goes +x.
  for (size_t i : {size_t{0}, size_t{4}}) {
    expect_pt(v[i].point, 0, 0, "point", i);
    expect_pt(v[i].in, 0, -1, "in", i);
    expect_pt(v[i].out, 1, 0, "out", i);
    expect_pt(v[i].direction, kRoot2, -kRoot2, "direction", i);
    EXPECT_TRUE(v[i].flags & kCl);
    EXPECT_TRUE(v[i].flags & kIn);
    EXPECT_TRUE(v[i].flags & kOut);
  }
  EXPECT_EQ(v[0].flags, kIn | kOut | kPS | kSS | kCl);
  EXPECT_EQ(v[4].flags, kIn | kOut | kPE | kSE | kCl);
  // The corners between are as in an open subpath.
  expect_pt(v[1].direction, kRoot2, kRoot2, "direction", 1);
  expect_pt(v[2].direction, -kRoot2, kRoot2, "direction", 2);
  expect_pt(v[3].direction, -kRoot2, -kRoot2, "direction", 3);
}

TEST(Vertices, ASegmentWithoutExtentHasNoDirectionAndTheVertexLooksPastIt) {
  // A zero-length line at the start, in the middle and at the end.
  auto v = vertices_of({op(M, 0, 0), op(L, 0, 0), op(L, 5, 0), op(L, 5, 0), op(L, 5, 5), op(L, 5, 5)});
  ASSERT_EQ(v.size(), 6u);
  // v1 is the zero line: nothing before it with extent, the real line after it.
  EXPECT_EQ(v[1].flags & (kIn | kOut), kOut);
  expect_pt(v[1].out, 1, 0, "out", 1);
  expect_pt(v[1].direction, 1, 0, "direction", 1);
  // v2 and v3 are both at (5, 0): the corner is the same whichever of them the marker is on.
  for (size_t i : {size_t{2}, size_t{3}}) {
    expect_pt(v[i].point, 5, 0, "point", i);
    expect_pt(v[i].in, 1, 0, "in", i);
    expect_pt(v[i].out, 0, 1, "out", i);
    expect_pt(v[i].direction, kRoot2, kRoot2, "direction", i);
  }
  // v5 is the zero line at the end: it comes in as the last real line did and goes nowhere.
  EXPECT_EQ(v[5].flags & (kIn | kOut), kIn);
  expect_pt(v[5].direction, 0, 1, "direction", 5);
}

TEST(Vertices, AClosingLineOfNoLengthIsPassedOverAndTheCornerIsStillTheStarts) {
  // The last point is the start, so the close has no extent; the corner at the start is between the last line and the first.
  auto v = vertices_of({op(M, 0, 0), op(L, 4, 0), op(L, 4, 4), op(L, 0, 0), op(Z)});
  ASSERT_EQ(v.size(), 5u);
  const double d = kRoot2; // (-4, -4) normalised
  for (size_t i : {size_t{0}, size_t{3}, size_t{4}}) {
    EXPECT_TRUE(v[i].flags & kIn) << i;
    EXPECT_TRUE(v[i].flags & kOut) << i;
    expect_pt(v[i].in, -d, -d, "in", i);
    expect_pt(v[i].out, 1, 0, "out", i);
  }
}

TEST(Vertices, AReversalTurnsTheIncomingDirectionAQuarterTurn) {
  auto v = vertices_of({op(M, 0, 0), op(L, 3, 0), op(L, 0, 0)});
  ASSERT_EQ(v.size(), 3u);
  expect_pt(v[1].in, 1, 0, "in", 1);
  expect_pt(v[1].out, -1, 0, "out", 1);
  EXPECT_EQ(v[1].direction.x, 0.0);
  EXPECT_EQ(v[1].direction.y, 1.0);
  // The same turn for a reversal going up: in (0, -1) turns to (1, 0).
  v = vertices_of({op(M, 0, 0), op(L, 0, -3), op(L, 0, 0)});
  EXPECT_EQ(v[1].direction.x, 1.0);
  EXPECT_EQ(v[1].direction.y, 0.0);
}

TEST(Vertices, ACurveIsReadByTheControlPointThatDiffersFromItsEnd) {
  // The first control point is the start, so the start direction is toward the second; the last is the end, so the
  // end direction is from the second.
  auto v = vertices_of({op(M, 0, 0), op(C, 0, 0, 4, 0, 4, 4), op(C, 4, 4, 4, 9, 4, 9)});
  ASSERT_EQ(v.size(), 3u);
  expect_pt(v[0].out, 1, 0, "out", 0);
  expect_pt(v[1].in, 0, 1, "in", 1);
  expect_pt(v[1].out, 0, 1, "out", 1);
  expect_pt(v[1].direction, 0, 1, "direction", 1);
  expect_pt(v[2].in, 0, 1, "in", 2);
  // A quadratic whose control point is its start, and a cubic of four equal points (no extent at all).
  v = vertices_of({op(M, 1, 1), op(Q, 1, 1, 2, 1), op(C, 2, 1, 2, 1, 2, 1), op(L, 2, 3)});
  ASSERT_EQ(v.size(), 4u);
  expect_pt(v[1].in, 1, 0, "in", 1);
  expect_pt(v[2].in, 1, 0, "in", 2);
  expect_pt(v[2].out, 0, 1, "out", 2);
  expect_pt(v[2].direction, kRoot2, kRoot2, "direction", 2);
}

TEST(Vertices, EachSubpathIsOwnAndEachMoveStartsOne) {
  auto v = vertices_of({op(M, 0, 0), op(L, 1, 0), op(M, 5, 5), op(L, 6, 5), op(L, 6, 6)});
  ASSERT_EQ(v.size(), 5u);
  EXPECT_EQ(v[0].subpath, 0u);
  EXPECT_EQ(v[1].subpath, 0u);
  EXPECT_EQ(v[2].subpath, 1u);
  EXPECT_EQ(v[1].flags, kIn | kSE);
  EXPECT_EQ(v[2].flags, kOut | kSS);
  EXPECT_EQ(v[4].flags, kIn | kPE | kSE);
  // The second subpath's first vertex does not see the first subpath.
  EXPECT_FALSE(v[2].flags & kIn);
  expect_pt(v[2].direction, 1, 0, "direction", 2);
}

TEST(Vertices, ADrawingVerbAfterACloseStartsASubpathWithNoMoveAndNoVertexOfItsOwn) {
  auto v = vertices_of({op(M, 0, 0), op(L, 2, 0), op(L, 2, 2), op(Z), op(L, 5, 5)});
  ASSERT_EQ(v.size(), 5u);
  EXPECT_EQ(v[3].flags, kIn | kOut | kSE | kCl); // the close: round to the first segment
  EXPECT_EQ(v[4].subpath, 1u);
  EXPECT_EQ(v[4].flags, kIn | kPE | kSE); // the line from (0, 0): open, no start move
  expect_pt(v[4].point, 5, 5, "point", 4);
  expect_pt(v[4].in, kRoot2, kRoot2, "in", 4);
  EXPECT_FALSE(v[4].flags & kSS);
}

TEST(Vertices, ALoneMoveAndAMoveAndCloseHaveAVertexAndTheDirectionZero) {
  auto v = vertices_of({op(M, 3, 4)});
  ASSERT_EQ(v.size(), 1u);
  EXPECT_EQ(v[0].flags, kPS | kPE | kSS | kSE);
  expect_pt(v[0].point, 3, 4, "point", 0);
  expect_pt(v[0].direction, 1, 0, "direction", 0);
  v = vertices_of({op(M, 3, 4), op(Z)});
  ASSERT_EQ(v.size(), 2u);
  EXPECT_EQ(v[1].flags, kPE | kSE | kCl);
  expect_pt(v[1].direction, 1, 0, "direction", 1);
  // Two lone moves in turn.
  v = vertices_of({op(M, 1, 1), op(M, 2, 2)});
  ASSERT_EQ(v.size(), 2u);
  EXPECT_EQ(v[0].flags, kPS | kSS | kSE);
  EXPECT_EQ(v[1].flags, kPE | kSS | kSE);
}

TEST(Vertices, AnEmptyPathHasNoVertexAndTheArgumentsAreChecked) {
  GSHGEOM_Path * p = build({});
  std::vector<GSHGEOM_Vertex> out;
  EXPECT_EQ(gshgeom_path_vertices(p, collect_one, &out), GSHGEOM_OK);
  EXPECT_TRUE(out.empty());
  size_t count = 7;
  EXPECT_EQ(gshgeom_path_vertex_count(p, &count), GSHGEOM_OK);
  EXPECT_EQ(count, 0u);
  EXPECT_EQ(gshgeom_path_vertices(nullptr, collect_one, &out), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_vertices(p, nullptr, &out), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_vertex_count(p, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_path_vertex_count(nullptr, &count), GSHGEOM_ERR_INVALID);
  gshgeom_path_destroy(p);
  GSHGEOM_Path * p3 = build({}, 3);
  EXPECT_EQ(gshgeom_path_vertices(p3, collect_one, &out), GSHGEOM_ERR_INVALID);
  count = 7;
  EXPECT_EQ(gshgeom_path_vertex_count(p3, &count), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(count, 7u) << "written only on success";
  gshgeom_path_destroy(p3);
}

TEST(Vertices, ACallbackThatRefusesStopsTheWalkAndItsResultIsReturned) {
  GSHGEOM_Path * p = build({op(M, 0, 0), op(L, 1, 0), op(L, 2, 0), op(L, 3, 0)});
  struct Stop {
    int calls = 0;
  } stop;
  auto fn = [](const GSHGEOM_Vertex *, void * ctx) -> GSHGEOM_Result {
    Stop * s = static_cast<Stop *>(ctx);
    return ++s->calls == 2 ? GSHGEOM_ERR_LIMIT : GSHGEOM_OK;
  };
  EXPECT_EQ(gshgeom_path_vertices(p, fn, &stop), GSHGEOM_ERR_LIMIT);
  EXPECT_EQ(stop.calls, 2);
  gshgeom_path_destroy(p);
}

TEST(Vertices, TheDirectionsAreUnitAtAnyMagnitude) {
  // A hair-thin segment and a huge one: the length is computed after scaling by the larger coordinate.
  auto v = vertices_of({op(M, 0, 0), op(L, 1e-300, 0), op(L, 1e-300, 1e-300)});
  expect_pt(v[1].in, 1, 0, "in", 1);
  expect_pt(v[1].out, 0, 1, "out", 1);
  v = vertices_of({op(M, -1e29, -1e29), op(L, 1e29, -1e29), op(L, 1e29, 1e29)});
  expect_pt(v[1].in, 1, 0, "in", 1);
  expect_pt(v[1].out, 0, 1, "out", 1);
  for (const auto & x : v) {
    EXPECT_EQ(std::hypot(x.direction.x, x.direction.y) > 0.9999999999999, true);
  }
}

TEST(Vertices, AHugeRunOfZeroLengthSegmentsIsPassedOverOnceNotOncePerVertex) {
  const size_t n = 200000;
  std::vector<Op> ops{op(M, 0, 0)};
  for (size_t i = 0; i < n; i++) {
    ops.push_back(op(L, 0, 0));
  }
  ops.push_back(op(L, 0, 7));
  GSHGEOM_Path * p = build(ops);
  struct Count {
    size_t seen = 0, wrong = 0;
  } c;
  auto fn = [](const GSHGEOM_Vertex * v, void * ctx) -> GSHGEOM_Result {
    Count * k = static_cast<Count *>(ctx);
    k->seen++;
    if ((v->flags & GSHGEOM_VERTEX_HAS_OUT) && v->out.y != 1.0) {
      k->wrong++;
    }
    return GSHGEOM_OK;
  };
  auto t0 = std::chrono::steady_clock::now();
  EXPECT_EQ(gshgeom_path_vertices(p, fn, &c), GSHGEOM_OK);
  auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - t0).count();
  EXPECT_EQ(c.seen, n + 2);
  EXPECT_EQ(c.wrong, 0u);
  EXPECT_LT(ms, 3000) << "a vertex rescanned the run of segments before the next one with extent";
  gshgeom_path_destroy(p);
}

// ---- the second reading ---------------------------------------------------------------------------------------

namespace {

struct Ref {
  bool has_in = false, has_out = false;
  long double in[2] = {0, 0}, out[2] = {0, 0};
};

// A segment's control points, from the path's own element readback.
struct RSeg {
  int n;
  double p[4][2];
  bool extent;
  long double s[2], e[2];
};

long double lnorm(long double x, long double y) {
  return std::sqrt(x * x + y * y);
}

RSeg rseg(const GSHGEOM_Element & el) {
  RSeg s{};
  s.n = gshgeom_curve_point_count(el.curve.kind) - 1;
  for (int i = 0; i <= s.n; i++) {
    s.p[i][0] = el.curve.p[i].x;
    s.p[i][1] = el.curve.p[i].y;
  }
  int a = -1, b = -1;
  for (int k = 1; k <= s.n; k++) {
    if (s.p[k][0] != s.p[0][0] || s.p[k][1] != s.p[0][1]) {
      a = k;
      break;
    }
  }
  for (int k = s.n - 1; k >= 0; k--) {
    if (s.p[k][0] != s.p[s.n][0] || s.p[k][1] != s.p[s.n][1]) {
      b = k;
      break;
    }
  }
  s.extent = a >= 0 && b >= 0;
  if (s.extent) {
    long double sx = (long double)s.p[a][0] - s.p[0][0], sy = (long double)s.p[a][1] - s.p[0][1];
    long double ex = (long double)s.p[s.n][0] - s.p[b][0], ey = (long double)s.p[s.n][1] - s.p[b][1];
    s.s[0] = sx / lnorm(sx, sy);
    s.s[1] = sy / lnorm(sx, sy);
    s.e[0] = ex / lnorm(ex, ey);
    s.e[1] = ey / lnorm(ex, ey);
  }
  return s;
}

} // namespace

TEST(Vertices, RandomPathsAgreeWithTheDefinitionReadTheLongWay) {
  uint64_t state = 0x1234abcd5678ef01ull;
  auto next = [&]() {
    state = state * 6364136223846793005ull + 1442695040888963407ull;
    return state >> 33;
  };
  auto coord = [&]() { return (double)((int)(next() % 7u) - 3); }; // few values: repeats and zero-length segments are common
  size_t compared = 0, bisected = 0, reversals = 0;
  for (int round = 0; round < 4000; round++) {
    std::vector<Op> ops;
    const size_t len = 1 + next() % 14;
    for (size_t i = 0; i < len; i++) {
      unsigned pick = (unsigned)(next() % 9u);
      if (i == 0 || pick == 0) {
        ops.push_back(op(M, coord(), coord()));
      } else if (pick <= 3) {
        ops.push_back(op(L, coord(), coord()));
      } else if (pick == 4) {
        ops.push_back(op(Q, coord(), coord(), coord(), coord()));
      } else if (pick == 5) {
        ops.push_back(op(C, coord(), coord(), coord(), coord(), coord(), coord()));
      } else if (pick == 6) {
        ops.push_back(op(Z));
      } else {
        ops.push_back(op(L, ops.back().verb == Z ? 0 : ops.back().v[0], 0)); // repeat an x: often a zero-length line
      }
    }
    GSHGEOM_Path * p = build(ops);
    const size_t n = gshgeom_path_element_count(p);
    std::vector<GSHGEOM_Element> els(n);
    for (size_t i = 0; i < n; i++) {
      ASSERT_EQ(gshgeom_path_element(p, i, &els[i]), GSHGEOM_OK);
    }
    // The subpaths: a move starts one; so does a drawing verb after a close; a close ends one.
    std::vector<size_t> sub_of(n);
    std::vector<size_t> first_of, last_of; // per subpath, [first, last]
    std::vector<char> closed_of;
    for (size_t i = 0; i < n; i++) {
      bool start = i == 0 || els[i].verb == GSHGEOM_VERB_MOVE || els[i - 1].verb == GSHGEOM_VERB_CLOSE;
      if (start) {
        first_of.push_back(i);
        last_of.push_back(i);
        closed_of.push_back(0);
      }
      sub_of[i] = first_of.size() - 1;
      last_of.back() = i;
      if (els[i].verb == GSHGEOM_VERB_CLOSE) {
        closed_of.back() = 1;
      }
    }
    std::vector<GSHGEOM_Vertex> got;
    ASSERT_EQ(gshgeom_path_vertices(p, collect_one, &got), GSHGEOM_OK);
    ASSERT_EQ(got.size(), n);
    for (size_t i = 0; i < n; i++) {
      const size_t s = sub_of[i], a = first_of[s], b = last_of[s];
      const bool closed = closed_of[s];
      // Nearest segment with extent at or before i, and strictly after i, then round the closed subpath.
      auto seg_at = [&](size_t j) {
        return els[j].verb == GSHGEOM_VERB_MOVE ? RSeg{0, {}, false, {}, {}} : rseg(els[j]);
      };
      Ref ref;
      long double in[2], out[2];
      bool found = false;
      for (size_t j = i + 1; j-- > a && !found;) {
        RSeg g = seg_at(j);
        if (g.extent) {
          found = true;
          in[0] = g.e[0];
          in[1] = g.e[1];
        }
      }
      for (size_t j = b + 1; j-- > a && !found && closed;) {
        RSeg g = seg_at(j);
        if (g.extent) {
          found = true;
          in[0] = g.e[0];
          in[1] = g.e[1];
        }
      }
      ref.has_in = found;
      if (found) {
        ref.in[0] = in[0];
        ref.in[1] = in[1];
      }
      found = false;
      for (size_t j = i + 1; j <= b && !found; j++) {
        RSeg g = seg_at(j);
        if (g.extent) {
          found = true;
          out[0] = g.s[0];
          out[1] = g.s[1];
        }
      }
      for (size_t j = a; j <= b && !found && closed; j++) {
        RSeg g = seg_at(j);
        if (g.extent) {
          found = true;
          out[0] = g.s[0];
          out[1] = g.s[1];
        }
      }
      ref.has_out = found;
      if (found) {
        ref.out[0] = out[0];
        ref.out[1] = out[1];
      }
      const GSHGEOM_Vertex & v = got[i];
      EXPECT_EQ(v.element, i);
      EXPECT_EQ(v.subpath, s);
      EXPECT_EQ(!!(v.flags & kIn), ref.has_in) << round << ":" << i;
      EXPECT_EQ(!!(v.flags & kOut), ref.has_out) << round << ":" << i;
      EXPECT_EQ(!!(v.flags & kCl), closed) << round << ":" << i;
      EXPECT_EQ(!!(v.flags & kPS), i == 0);
      EXPECT_EQ(!!(v.flags & kPE), i + 1 == n);
      EXPECT_EQ(!!(v.flags & kSS), els[i].verb == GSHGEOM_VERB_MOVE);
      EXPECT_EQ(!!(v.flags & kSE), i == b);
      const GSHGEOM_Point & end = els[i].curve.p[gshgeom_curve_point_count(els[i].curve.kind) - 1];
      EXPECT_EQ(v.point.x, end.x);
      EXPECT_EQ(v.point.y, end.y);
      const long double tol = 4e-16L;
      if (ref.has_in) {
        EXPECT_NEAR((double)v.in.x, (double)ref.in[0], tol);
        EXPECT_NEAR((double)v.in.y, (double)ref.in[1], tol);
      } else {
        EXPECT_EQ(v.in.x, 0.0);
        EXPECT_EQ(v.in.y, 0.0);
      }
      if (ref.has_out) {
        EXPECT_NEAR((double)v.out.x, (double)ref.out[0], tol);
        EXPECT_NEAR((double)v.out.y, (double)ref.out[1], tol);
      }
      // The marker direction, by the angle average the specification states (the shorter arc).
      long double dx, dy;
      if (ref.has_in && ref.has_out) {
        long double ai = std::atan2(ref.in[1], ref.in[0]), ao = std::atan2(ref.out[1], ref.out[0]);
        long double diff = ao - ai;
        const long double pi = 3.14159265358979323846264338327950288L;
        if (std::fabs((double)(std::fabs(diff) - pi)) < 1e-9L) {
          reversals++; // a reversal: the choice of side is the library's own (documented); the tests above pin it
          dx = v.direction.x;
          dy = v.direction.y;
        } else {
          if (diff > pi) diff -= 2 * pi;
          if (diff < -pi) diff += 2 * pi;
          long double am = ai + diff / 2;
          dx = std::cos(am);
          dy = std::sin(am);
          bisected++;
        }
      } else if (ref.has_in) {
        dx = ref.in[0];
        dy = ref.in[1];
      } else if (ref.has_out) {
        dx = ref.out[0];
        dy = ref.out[1];
      } else {
        dx = 1;
        dy = 0;
      }
      // The unit vector of a sum has the sum's relative error: near a reversal the sum is small and the error large.
      long double dtol = 2e-15L;
      if (ref.has_in && ref.has_out) {
        const long double sn = lnorm(ref.in[0] + ref.out[0], ref.in[1] + ref.out[1]);
        if (sn > 1e-3L) {
          dtol = 8e-16L / sn + 2e-15L;
        } else {
          dtol = 1.0L; // too near a reversal for the angle average to say which side; the unit length below still holds
        }
      }
      EXPECT_NEAR((double)v.direction.x, (double)dx, (double)dtol) << round << ":" << i;
      EXPECT_NEAR((double)v.direction.y, (double)dy, (double)dtol) << round << ":" << i;
      EXPECT_NEAR(std::hypot(v.direction.x, v.direction.y), 1.0, 4e-16) << round << ":" << i;
      compared++;
    }
    gshgeom_path_destroy(p);
  }
  EXPECT_GT(compared, 20000u);
  EXPECT_GT(bisected, 3000u);
  EXPECT_GT(reversals, 50u) << "the population holds reversals";
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
