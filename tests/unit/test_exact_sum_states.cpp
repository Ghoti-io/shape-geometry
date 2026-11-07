/**
 * @file
 *
 * The exact sum's conversion to a double on the states the public adders cannot reach (the accumulator's top word, a
 * sum already marked bad) and the ranges that make some of its steps unobservable.  Each test is a separation of a
 * plant that the arc grid could not tell from the real code.
 *
 * Copyright 2026 by Corey Pennycuff
 */

#include "../test_helpers.h"

extern "C" {
#include "../../src/exact/exact.h"
}

#include <cmath>
#include <cstdint>
#include <cstring>

namespace {

uint64_t bits(double x) {
  uint64_t b;
  std::memcpy(&b, &x, sizeof b);
  return b;
}

}  // namespace

TEST(ExactSumStates, ABadSumConvertsToZeroWhateverItHolds) {
  GSHGEOM_ExactSum s;
  gshgeom_exact_sum_init(&s);
  const double one[1] = {3.0};
  gshgeom_exact_sum_add(&s, 0, 1, one);
  EXPECT_EQ(gshgeom_exact_sum_to_double(&s), 3.0);
  const double nan_term[1] = {std::nan("")};
  gshgeom_exact_sum_add(&s, 0, 1, nan_term);
  EXPECT_EQ(s.bad, 1);
  EXPECT_EQ(bits(gshgeom_exact_sum_to_double(&s)), bits(0.0));  // the 3 it holds is not read
  EXPECT_EQ(gshgeom_exact_sum_sign(&s), 0);
}

TEST(ExactSumStates, TheTopWordDecidesTheSignAndTheMagnitude) {
  // The adders reach the top word only past 2^12680, a sum no twelve doubles and no feasible count of terms make, so the
  // words are written directly: the sign is bit 31 of the top word, and a value in the top word is read.
  GSHGEOM_ExactSum s;
  gshgeom_exact_sum_init(&s);
  s.w[GSHGEOM_EXACT_SUM_WORDS - 1] = 0x40000000u;  // large and positive: bit 30 is a magnitude bit, not the sign
  EXPECT_EQ(bits(gshgeom_exact_sum_to_double(&s)), bits(HUGE_VAL));
  gshgeom_exact_sum_init(&s);
  s.w[GSHGEOM_EXACT_SUM_WORDS - 1] = 0xc0000000u;  // the same magnitude, negative (two's complement of 0x40000000)
  EXPECT_EQ(bits(gshgeom_exact_sum_to_double(&s)), bits(-HUGE_VAL));
  gshgeom_exact_sum_init(&s);
  s.w[GSHGEOM_EXACT_SUM_WORDS - 1] = 1u;  // a one in the top word and nothing under it
  EXPECT_EQ(bits(gshgeom_exact_sum_to_double(&s)), bits(HUGE_VAL));
  EXPECT_EQ(gshgeom_exact_sum_sign(&s), 1);
}

TEST(ExactSumStates, NothingBelowTheFirstBitSurvivesTheScale) {
  // The conversion reads below bit 0 only when the leading bit is under bit 63, and a value that small is far under the
  // smallest double: every such accumulator converts to zero, with bit 0 set or not.
  for (unsigned top = 0; top < 64; top++) {
    for (int lowbit = 0; lowbit < 2; lowbit++) {
      GSHGEOM_ExactSum s;
      gshgeom_exact_sum_init(&s);
      s.w[1] = top >= 32 ? (1u << (top - 32)) : 0u;
      s.w[0] = (top < 32 ? (1u << top) : 0u) | (unsigned)lowbit;
      if (top == 0 && lowbit == 0) {
        continue;
      }
      EXPECT_EQ(bits(gshgeom_exact_sum_to_double(&s)), bits(0.0)) << top << " " << lowbit;
    }
  }
}

TEST(ExactSumStates, ARoundUpToThePowerOfTwoAndTheLargeScaleStepsAgreeWithTheHardwareSum) {
  // 2^53 - 1 plus a half rounds (a tie to even) to 2^53, the case that renormalizes; the result must equal the
  // hardware sum at every scale across the steps of the scale loop, including those that overflow or underflow.
  for (int e = -1100; e <= 1023; e += 1) {
    const double a = std::ldexp(9007199254740991.0, e - 53);  // (2^53 - 1) 2^(e-53)
    const double b = std::ldexp(0.5, e - 53);
    GSHGEOM_ExactSum s;
    gshgeom_exact_sum_init(&s);
    const double fa[1] = {a}, fb[1] = {b};
    gshgeom_exact_sum_add(&s, 0, 1, fa);
    gshgeom_exact_sum_add(&s, 0, 1, fb);
    volatile double want = a + b;
    // the conversion rounds the 53-bit value and then scales it, so it agrees with the hardware sum wherever the sum is
    // a normal number (a subnormal result is rounded twice, once to 53 bits and once at the scale)
    if (std::fabs((double)want) >= 0x1p-1022 || (double)want == 0.0) {
      EXPECT_EQ(bits(gshgeom_exact_sum_to_double(&s)), bits(want)) << e;
    }
  }
  // 1000 and the step below it: the value just across 2^1000 and 2^999, and the scale 2^-1000 and 2^-999
  for (int e = -1010; e <= -990; e++) {
    const double a = std::ldexp(1.0, e);
    GSHGEOM_ExactSum s;
    gshgeom_exact_sum_init(&s);
    const double fa[1] = {a};
    gshgeom_exact_sum_add(&s, 0, 1, fa);
    EXPECT_EQ(bits(gshgeom_exact_sum_to_double(&s)), bits(a)) << e;
  }
  for (int e = 990; e <= 1023; e++) {
    const double a = std::ldexp(1.0, e);
    GSHGEOM_ExactSum s;
    gshgeom_exact_sum_init(&s);
    const double fa[1] = {a};
    gshgeom_exact_sum_add(&s, 0, 1, fa);
    EXPECT_EQ(bits(gshgeom_exact_sum_to_double(&s)), bits(a)) << e;
  }
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
