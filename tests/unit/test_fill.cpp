/**
 * @file
 *
 * The fill rules (AD-4): the pinned values, the names, and the closed count.
 *
 * Copyright 2026 by Corey Pennycuff
 */

#include "test_helpers.h"

#include <ghoti.io/shape-geometry/fill.h>

#include <cstring>

TEST(Fill, TheValuesArePinned) {
  // AD-4: 0 nonzero (the default), 1 even-odd, anything further appended.
  EXPECT_EQ(static_cast<int>(GSHGEOM_FILL_NONZERO), 0);
  EXPECT_EQ(static_cast<int>(GSHGEOM_FILL_EVEN_ODD), 1);
  EXPECT_EQ(static_cast<int>(GSHGEOM_FILL_RULE_COUNT), 2);
}

TEST(Fill, EveryRuleHasItsName) {
  EXPECT_STREQ(gshgeom_fill_rule_string(GSHGEOM_FILL_NONZERO), "non-zero");
  EXPECT_STREQ(gshgeom_fill_rule_string(GSHGEOM_FILL_EVEN_ODD), "even-odd");
}

TEST(Fill, AValueOutsideTheEnumIsNamedAndNeverNull) {
  const char * unknown = "unknown fill rule";

  EXPECT_STREQ(gshgeom_fill_rule_string(GSHGEOM_FILL_RULE_COUNT), unknown);
  EXPECT_STREQ(gshgeom_fill_rule_string(static_cast<GSHGEOM_Fill_Rule>(7)),
      unknown);
  EXPECT_STREQ(gshgeom_fill_rule_string(static_cast<GSHGEOM_Fill_Rule>(-1)),
      unknown);
}

TEST(Fill, TheNamesAreDistinct) {
  for (int a = 0; a < GSHGEOM_FILL_RULE_COUNT; ++a) {
    for (int b = a + 1; b < GSHGEOM_FILL_RULE_COUNT; ++b) {
      EXPECT_STRNE(gshgeom_fill_rule_string(static_cast<GSHGEOM_Fill_Rule>(a)),
          gshgeom_fill_rule_string(static_cast<GSHGEOM_Fill_Rule>(b)));
    }
  }
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
