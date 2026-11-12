/**
 * @file
 *
 * Hit testing against a flattened path (hit.h, HT1 to HT6).  The hand cases pin each rule; the random cases compare every point
 * of a grid, on and off the edges, with a reference that shares no code with the library: Hormann and Agathos's quadrant count,
 * which sums the quadrant changes of the contour around the point and decides a half-turn by the sign of a cross product.
 */

#include "test_helpers.h"

#include <cstdint>
#include <random>
#include <utility>
#include <vector>

namespace {

using Dp = GSHGEOM_Device_Point;
using Poly = std::vector<Dp>;

/** A flat path of the subpaths @p polys (each as points; @p open[i] true leaves it unclosed in the flag, which HT2 ignores). */
GSHGEOM_Flat_Path * flat(const std::vector<Poly> & polys, const std::vector<bool> & open = {}) {
  std::vector<Dp> pts;
  std::vector<GSHGEOM_Flat_Subpath> subs;
  for (size_t i = 0; i < polys.size(); i++) {
    GSHGEOM_Flat_Subpath s;
    s.first = pts.size();
    s.count = polys[i].size();
    s.closed = (i < open.size() && open[i]) ? 0 : 1;
    pts.insert(pts.end(), polys[i].begin(), polys[i].end());
    subs.push_back(s);
  }
  GSHGEOM_Flat_Path * fp = nullptr;
  EXPECT_EQ(gshgeom_flat_path_from_points(pts.data(), pts.size(), subs.data(), subs.size(), nullptr, nullptr, &fp), GSHGEOM_OK);
  return fp;
}

struct Answer {
  int64_t w;
  int on;
};

Answer winding(const GSHGEOM_Flat_Path * fp, int32_t x, int32_t y) {
  Answer a{12345, 12345};
  EXPECT_EQ(gshgeom_flat_path_winding(fp, Dp{x, y}, &a.w, &a.on), GSHGEOM_OK);
  return a;
}

GSHGEOM_Hit hit(const GSHGEOM_Flat_Path * fp, int32_t x, int32_t y, GSHGEOM_Fill_Rule rule) {
  GSHGEOM_Hit h = GSHGEOM_HIT_COUNT;
  EXPECT_EQ(gshgeom_flat_path_hit(fp, Dp{x, y}, rule, &h), GSHGEOM_OK);
  return h;
}

/** The square (0,0) (10,0) (10,10) (0,10): on a y-down screen it runs right, down, left, up. */
Poly square(int32_t x0, int32_t y0, int32_t x1, int32_t y1) {
  return Poly{{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}};
}

Poly reversed(Poly p) {
  Poly r(p.rbegin(), p.rend());
  return r;
}

TEST(Hit, ASquareHasWindingMinusOneInsideAndPlusOneReversedAndTheSignIsTheDocumentedOne) {
  GSHGEOM_Flat_Path * fp = flat({square(0, 0, 10, 10)});
  Answer a = winding(fp, 5, 5);
  EXPECT_EQ(a.w, -1) << "the edge up the left side has the point to its right: an edge directed up counts -1 (HT4)";
  EXPECT_EQ(a.on, 0);
  gshgeom_flat_path_destroy(fp);
  fp = flat({reversed(square(0, 0, 10, 10))});
  EXPECT_EQ(winding(fp, 5, 5).w, 1);
  gshgeom_flat_path_destroy(fp);
}

TEST(Hit, OutsideInsideAndBothRules) {
  GSHGEOM_Flat_Path * fp = flat({square(0, 0, 10, 10)});
  for (GSHGEOM_Fill_Rule rule : {GSHGEOM_FILL_NONZERO, GSHGEOM_FILL_EVEN_ODD}) {
    EXPECT_EQ(hit(fp, 5, 5, rule), GSHGEOM_HIT_INSIDE);
    EXPECT_EQ(hit(fp, 11, 5, rule), GSHGEOM_HIT_OUTSIDE);
    EXPECT_EQ(hit(fp, -1, 5, rule), GSHGEOM_HIT_OUTSIDE);
    EXPECT_EQ(hit(fp, 5, -1, rule), GSHGEOM_HIT_OUTSIDE);
    EXPECT_EQ(hit(fp, 5, 11, rule), GSHGEOM_HIT_OUTSIDE);
  }
  gshgeom_flat_path_destroy(fp);
}

TEST(Hit, EveryPointOnAnEdgeOrAVertexIsOnTheBoundaryUnderBothRulesAndWhicheverWayTheContourRuns) {
  for (int flip = 0; flip < 2; flip++) {
    Poly sq = square(0, 0, 10, 10);
    GSHGEOM_Flat_Path * fp = flat({flip ? reversed(sq) : sq});
    const Dp on[] = {{0, 0}, {10, 0}, {10, 10}, {0, 10}, {5, 0}, {10, 5}, {5, 10}, {0, 5}, {1, 0}, {10, 9}};
    for (const Dp & p : on) {
      Answer a = winding(fp, p.x, p.y);
      EXPECT_EQ(a.on, 1) << p.x << "," << p.y;
      EXPECT_EQ(a.w, 0) << "a boundary point reports winding 0";
      EXPECT_EQ(hit(fp, p.x, p.y, GSHGEOM_FILL_NONZERO), GSHGEOM_HIT_BOUNDARY);
      EXPECT_EQ(hit(fp, p.x, p.y, GSHGEOM_FILL_EVEN_ODD), GSHGEOM_HIT_BOUNDARY);
    }
    const Dp off[] = {{-1, 0}, {11, 0}, {5, -1}, {5, 11}, {1, 1}, {9, 9}};
    for (const Dp & p : off) {
      EXPECT_EQ(winding(fp, p.x, p.y).on, 0) << p.x << "," << p.y;
    }
    gshgeom_flat_path_destroy(fp);
  }
}

TEST(Hit, ADiagonalEdgeIsOnTheBoundaryExactlyWhereTheLatticeHasAPointOnIt) {
  GSHGEOM_Flat_Path * fp = flat({Poly{{0, 0}, {6, 4}, {0, 4}}});
  EXPECT_EQ(winding(fp, 3, 2).on, 1) << "the midpoint of (0,0) to (6,4)";
  EXPECT_EQ(winding(fp, 6, 4).on, 1);
  EXPECT_EQ(winding(fp, 3, 3).on, 0) << "above the midpoint's height on the line, not on it";
  EXPECT_EQ(winding(fp, 4, 2).on, 0);
  EXPECT_EQ(winding(fp, 2, 1).on, 0) << "the line passes (1.5, 1), not (2, 1)";
  EXPECT_EQ(hit(fp, 2, 3, GSHGEOM_FILL_NONZERO), GSHGEOM_HIT_INSIDE);
  gshgeom_flat_path_destroy(fp);
}

TEST(Hit, AVertexOnTheRayIsCountedOnceWhenTheContourCrossesAndNotAtAllWhenItOnlyTouches) {
  /* A diamond with its left vertex at the height of the probe: the probe to the right of it, at the vertex height. */
  GSHGEOM_Flat_Path * fp = flat({Poly{{0, 5}, {5, 0}, {10, 5}, {5, 10}}});
  EXPECT_EQ(hit(fp, 5, 5, GSHGEOM_FILL_NONZERO), GSHGEOM_HIT_INSIDE) << "the ray to the left meets the vertex (0,5) where the contour crosses it";
  EXPECT_EQ(winding(fp, 5, 5).w, -1);
  EXPECT_EQ(hit(fp, 11, 5, GSHGEOM_FILL_NONZERO), GSHGEOM_HIT_OUTSIDE) << "past the right vertex: both vertices crossed, in and out";
  EXPECT_EQ(winding(fp, 11, 5).w, 0);
  EXPECT_EQ(hit(fp, -1, 5, GSHGEOM_FILL_NONZERO), GSHGEOM_HIT_OUTSIDE);
  EXPECT_EQ(winding(fp, 20, 0).w, 0) << "the top vertex is touched, not crossed";
  EXPECT_EQ(winding(fp, 20, 10).w, 0) << "so is the bottom one";
  EXPECT_EQ(winding(fp, 20, 5).w, 0);
  gshgeom_flat_path_destroy(fp);
}

TEST(Hit, AnEdgeAlongTheRayCountsNothingAndNeitherDoesItsNeighbourPairTwice) {
  /* An L whose horizontal edge lies at the probe's height. */
  GSHGEOM_Flat_Path * fp = flat({Poly{{0, 0}, {10, 0}, {10, 5}, {20, 5}, {20, 10}, {0, 10}}});
  EXPECT_EQ(winding(fp, 30, 5).w, 0) << "to the right of the whole shape at the height of its step";
  EXPECT_EQ(winding(fp, 15, 5).on, 1);
  EXPECT_EQ(winding(fp, 15, 7).w, -1) << "below the step, inside";
  EXPECT_EQ(winding(fp, 15, 3).w, 0) << "above the step, outside";
  EXPECT_EQ(winding(fp, 5, 5).w, -1);
  gshgeom_flat_path_destroy(fp);
}

TEST(Hit, EvenOddAndNonzeroDifferWhereAContourWindsTwice) {
  /* The same square twice in the same direction (winding -2) and a square with an inner square the other way (a hole). */
  GSHGEOM_Flat_Path * twice = flat({square(0, 0, 10, 10), square(0, 0, 10, 10)});
  EXPECT_EQ(winding(twice, 5, 5).w, -2);
  EXPECT_EQ(hit(twice, 5, 5, GSHGEOM_FILL_NONZERO), GSHGEOM_HIT_INSIDE);
  EXPECT_EQ(hit(twice, 5, 5, GSHGEOM_FILL_EVEN_ODD), GSHGEOM_HIT_OUTSIDE);
  gshgeom_flat_path_destroy(twice);
  GSHGEOM_Flat_Path * hole = flat({square(0, 0, 10, 10), reversed(square(3, 3, 7, 7))});
  EXPECT_EQ(winding(hole, 5, 5).w, 0);
  EXPECT_EQ(hit(hole, 5, 5, GSHGEOM_FILL_NONZERO), GSHGEOM_HIT_OUTSIDE);
  EXPECT_EQ(hit(hole, 5, 5, GSHGEOM_FILL_EVEN_ODD), GSHGEOM_HIT_OUTSIDE);
  EXPECT_EQ(hit(hole, 1, 5, GSHGEOM_FILL_NONZERO), GSHGEOM_HIT_INSIDE);
  EXPECT_EQ(hit(hole, 3, 5, GSHGEOM_FILL_NONZERO), GSHGEOM_HIT_BOUNDARY);
  gshgeom_flat_path_destroy(hole);
  GSHGEOM_Flat_Path * same = flat({square(0, 0, 10, 10), square(3, 3, 7, 7)});
  EXPECT_EQ(hit(same, 5, 5, GSHGEOM_FILL_NONZERO), GSHGEOM_HIT_INSIDE);
  EXPECT_EQ(hit(same, 5, 5, GSHGEOM_FILL_EVEN_ODD), GSHGEOM_HIT_OUTSIDE) << "two turns: even";
  gshgeom_flat_path_destroy(same);
}

TEST(Hit, AnOpenSubpathIsClosedByTheLineBackToItsStartAndTheFlagIsIgnored) {
  Poly tri{{0, 0}, {10, 0}, {0, 10}};
  GSHGEOM_Flat_Path * closed = flat({tri});
  GSHGEOM_Flat_Path * open = flat({tri}, {true});
  for (int x = -2; x <= 12; x++) {
    for (int y = -2; y <= 12; y++) {
      Answer a = winding(closed, x, y);
      Answer b = winding(open, x, y);
      EXPECT_EQ(a.w, b.w) << x << "," << y;
      EXPECT_EQ(a.on, b.on) << x << "," << y;
    }
  }
  EXPECT_EQ(winding(open, 5, 5).on, 1) << "on the closing line from (0,10) to (10,0)... the hypotenuse x+y=10";
  EXPECT_EQ(winding(open, 2, 8).on, 1);
  gshgeom_flat_path_destroy(closed);
  gshgeom_flat_path_destroy(open);
}

TEST(Hit, ASubpathWithNoEdgeTakesNoPartNotEvenAsABoundary) {
  GSHGEOM_Flat_Path * one = flat({Poly{{4, 4}}});
  EXPECT_EQ(winding(one, 4, 4).on, 0);
  EXPECT_EQ(hit(one, 4, 4, GSHGEOM_FILL_NONZERO), GSHGEOM_HIT_OUTSIDE);
  gshgeom_flat_path_destroy(one);
  GSHGEOM_Flat_Path * same = flat({Poly{{4, 4}, {4, 4}, {4, 4}}});
  EXPECT_EQ(winding(same, 4, 4).on, 0);
  gshgeom_flat_path_destroy(same);
  GSHGEOM_Flat_Path * seg = flat({Poly{{0, 0}, {8, 0}}});
  EXPECT_EQ(winding(seg, 4, 0).on, 1) << "a two-point subpath is a line out and back: its points are on the boundary";
  EXPECT_EQ(winding(seg, 4, 1).w, 0);
  gshgeom_flat_path_destroy(seg);
  GSHGEOM_Flat_Path * dup = flat({Poly{{0, 0}, {0, 0}, {8, 0}, {8, 8}, {8, 8}, {0, 8}}});
  EXPECT_EQ(winding(dup, 4, 4).w, -1) << "repeated vertices add no edge";
  gshgeom_flat_path_destroy(dup);
}

TEST(Hit, TheEmptyAndTheNullPathEncloseNothing) {
  for (int i = 0; i < 2; i++) {
    GSHGEOM_Flat_Path * fp = i ? nullptr : flat({});
    Answer a = winding(fp, 0, 0);
    EXPECT_EQ(a.w, 0);
    EXPECT_EQ(a.on, 0);
    EXPECT_EQ(hit(fp, 3, 3, GSHGEOM_FILL_NONZERO), GSHGEOM_HIT_OUTSIDE);
    gshgeom_flat_path_destroy(fp);
  }
}

TEST(Hit, TheArgumentsAreChecked) {
  GSHGEOM_Flat_Path * fp = flat({square(0, 0, 10, 10)});
  int64_t w = 77;
  int on = 77;
  EXPECT_EQ(gshgeom_flat_path_winding(fp, Dp{1, 1}, nullptr, &on), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(on, 77) << "nothing is written on an error";
  EXPECT_EQ(gshgeom_flat_path_winding(fp, Dp{1, 1}, &w, nullptr), GSHGEOM_OK) << "the boundary flag is optional";
  EXPECT_EQ(w, -1);
  GSHGEOM_Hit h = GSHGEOM_HIT_COUNT;
  EXPECT_EQ(gshgeom_flat_path_hit(fp, Dp{1, 1}, GSHGEOM_FILL_NONZERO, nullptr), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_flat_path_hit(fp, Dp{1, 1}, GSHGEOM_FILL_RULE_COUNT, &h), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(gshgeom_flat_path_hit(fp, Dp{1, 1}, static_cast<GSHGEOM_Fill_Rule>(-1), &h), GSHGEOM_ERR_INVALID);
  EXPECT_EQ(h, GSHGEOM_HIT_COUNT) << "nothing is written on an error";
  gshgeom_flat_path_destroy(fp);
}

TEST(Hit, TheValuesAreNamedAndPinned) {
  EXPECT_EQ(static_cast<int>(GSHGEOM_HIT_OUTSIDE), 0);
  EXPECT_EQ(static_cast<int>(GSHGEOM_HIT_INSIDE), 1);
  EXPECT_EQ(static_cast<int>(GSHGEOM_HIT_BOUNDARY), 2);
  EXPECT_STREQ(gshgeom_hit_string(GSHGEOM_HIT_OUTSIDE), "outside");
  EXPECT_STREQ(gshgeom_hit_string(GSHGEOM_HIT_INSIDE), "inside");
  EXPECT_STREQ(gshgeom_hit_string(GSHGEOM_HIT_BOUNDARY), "boundary");
  EXPECT_STREQ(gshgeom_hit_string(GSHGEOM_HIT_COUNT), "unknown hit");
  EXPECT_STREQ(gshgeom_hit_string(static_cast<GSHGEOM_Hit>(-5)), "unknown hit");
}

TEST(Hit, TheLimitsOfTheLatticeDoNotOverflowTheSideTest) {
  const int32_t lo = INT32_MIN;
  const int32_t hi = INT32_MAX;
  /* A triangle spanning the whole lattice: the products of the differences are near 2^64. */
  GSHGEOM_Flat_Path * fp = flat({Poly{{lo, lo}, {hi, lo}, {hi, hi}}});
  /* The triangle is the points with x > y (y down): top left, top right, bottom right. */
  EXPECT_EQ(hit(fp, hi - 1, lo + 1, GSHGEOM_FILL_NONZERO), GSHGEOM_HIT_INSIDE);
  EXPECT_EQ(hit(fp, lo + 1, hi - 1, GSHGEOM_FILL_NONZERO), GSHGEOM_HIT_OUTSIDE);
  EXPECT_EQ(hit(fp, 0, 0, GSHGEOM_FILL_NONZERO), GSHGEOM_HIT_BOUNDARY) << "the origin is on the diagonal y = x from (lo,lo) to (hi,hi)";
  gshgeom_flat_path_destroy(fp);
}

/** The reference: the quadrant count of Hormann and Agathos over @p Wide arithmetic; the winding it returns counts a counterclockwise
 *  contour of the mathematical plane (x right, y up, read straight off the coordinates) +1, so the library's number is its negation. */
template <typename Wide>
Answer reference(const std::vector<Poly> & polys, int32_t px, int32_t py) {
  int64_t total = 0;
  auto quadrant = [](int64_t dx, int64_t dy) {
    if (dx > 0 && dy >= 0) return 0;
    if (dx <= 0 && dy > 0) return 1;
    if (dx < 0 && dy <= 0) return 2;
    return 3;
  };
  for (const Poly & poly : polys) {
    const size_t n = poly.size();
    for (size_t i = 0; i < n; i++) {
      const Dp a = poly[i];
      const Dp b = poly[(i + 1) % n];
      if (a.x == b.x && a.y == b.y) {
        continue;
      }
      const int64_t ax = (int64_t)a.x - px, ay = (int64_t)a.y - py, bx = (int64_t)b.x - px, by = (int64_t)b.y - py;
      const Wide cross = (Wide)ax * by - (Wide)ay * bx;
      const Wide dot = (Wide)ax * bx + (Wide)ay * by;
      if (cross == 0 && dot <= 0) {
        return Answer{0, 1}; /* the point is on the segment: the vectors to its ends are opposite or one is zero */
      }
      int d = quadrant(bx, by) - quadrant(ax, ay);
      if (d == 3) {
        d = -1;
      } else if (d == -3) {
        d = 1;
      } else if (d == 2 || d == -2) {
        d = cross > 0 ? 2 : -2;
      }
      total += d;
    }
  }
  return Answer{-(total / 4), 0};
}

Poly random_poly(std::mt19937_64 & rng, int32_t lo, int32_t hi, size_t n) {
  std::uniform_int_distribution<int32_t> d(lo, hi);
  Poly p;
  for (size_t i = 0; i < n; i++) {
    p.push_back(Dp{d(rng), d(rng)});
  }
  return p;
}

TEST(Hit, EveryPointOfAGridAgreesWithTheQuadrantCountOnSmallRandomPathsFullOfDegeneracies) {
  std::mt19937_64 rng(0x6869742dULL);
  size_t on_seen = 0, nonzero_seen = 0, multi_seen = 0;
  for (int trial = 0; trial < 400; trial++) {
    std::vector<Poly> polys;
    const size_t subs = 1 + rng() % 3;
    for (size_t s = 0; s < subs; s++) {
      polys.push_back(random_poly(rng, -6, 6, 1 + rng() % 9));
    }
    GSHGEOM_Flat_Path * fp = flat(polys);
    for (int32_t y = -8; y <= 8; y++) {
      for (int32_t x = -8; x <= 8; x++) {
        const Answer got = winding(fp, x, y);
        const Answer want = reference<int64_t>(polys, x, y);
        ASSERT_EQ(got.on, want.on) << "trial " << trial << " at " << x << "," << y;
        ASSERT_EQ(got.w, want.w) << "trial " << trial << " at " << x << "," << y;
        on_seen += got.on;
        nonzero_seen += got.w != 0;
        multi_seen += (got.w > 1 || got.w < -1);
        const GSHGEOM_Hit hn = hit(fp, x, y, GSHGEOM_FILL_NONZERO);
        const GSHGEOM_Hit he = hit(fp, x, y, GSHGEOM_FILL_EVEN_ODD);
        if (got.on) {
          ASSERT_EQ(hn, GSHGEOM_HIT_BOUNDARY);
          ASSERT_EQ(he, GSHGEOM_HIT_BOUNDARY);
        } else {
          ASSERT_EQ(hn, got.w != 0 ? GSHGEOM_HIT_INSIDE : GSHGEOM_HIT_OUTSIDE);
          ASSERT_EQ(he, (got.w & 1) ? GSHGEOM_HIT_INSIDE : GSHGEOM_HIT_OUTSIDE);
        }
      }
    }
    gshgeom_flat_path_destroy(fp);
  }
  EXPECT_GT(on_seen, 2000u) << "the grid must reach the edges and vertices";
  EXPECT_GT(nonzero_seen, 5000u);
  EXPECT_GT(multi_seen, 100u) << "self-overlapping contours must be in the population";
}

TEST(Hit, LargeCoordinatesAgreeWithTheQuadrantCountWhereTheProductsFitSixtyFourBits) {
  std::mt19937_64 rng(0x6269672dULL);
  const int32_t m = 1 << 29;
  std::uniform_int_distribution<int32_t> d(-m, m);
  for (int trial = 0; trial < 3000; trial++) {
    std::vector<Poly> polys{random_poly(rng, -m, m, 3 + rng() % 8)};
    GSHGEOM_Flat_Path * fp = flat(polys);
    for (int k = 0; k < 20; k++) {
      const int32_t x = d(rng), y = d(rng);
      const Answer got = winding(fp, x, y);
      const Answer want = reference<int64_t>(polys, x, y);
      ASSERT_EQ(got.on, want.on);
      ASSERT_EQ(got.w, want.w) << "trial " << trial;
    }
    gshgeom_flat_path_destroy(fp);
  }
}

#ifdef __SIZEOF_INT128__
__extension__ typedef __int128 Int128;

TEST(Hit, TheWholeInt32LatticeAgreesWithTheQuadrantCountInSixtyFourBitsOfProductAndBeyond) {
  std::mt19937_64 rng(0x77696465ULL);
  size_t on_seen = 0;
  for (int trial = 0; trial < 6000; trial++) {
    std::vector<Poly> polys{random_poly(rng, INT32_MIN, INT32_MAX, 3 + rng() % 6)};
    if (trial % 3 == 0) {
      /* A repeated vertex, which adds no edge. */
      Poly & p = polys[0];
      p[1].x = p[0].x;
      p[1].y = p[0].y;
    }
    GSHGEOM_Flat_Path * fp = flat(polys);
    std::uniform_int_distribution<int32_t> d(INT32_MIN, INT32_MAX);
    for (int k = 0; k < 12; k++) {
      int32_t x = d(rng), y = d(rng);
      if (k == 0) {
        x = polys[0][0].x;
        y = polys[0][0].y;
      } else if (k == 1) {
        const Dp a = polys[0][1], b = polys[0][2]; /* the midpoint of an edge, on it when the span is even in both axes */
        x = (int32_t)(((int64_t)a.x + b.x) / 2);
        y = (int32_t)(((int64_t)a.y + b.y) / 2);
      }
      const Answer got = winding(fp, x, y);
      const Answer want = reference<Int128>(polys, x, y);
      ASSERT_EQ(got.on, want.on) << "trial " << trial;
      ASSERT_EQ(got.w, want.w) << "trial " << trial;
      on_seen += got.on;
    }
    gshgeom_flat_path_destroy(fp);
  }
  EXPECT_GT(on_seen, 1000u);
}
#endif

}  // namespace

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
