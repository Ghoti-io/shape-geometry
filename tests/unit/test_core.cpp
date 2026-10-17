/**
 * @file
 *
 * The version: the placeholder test the scaffold ships with.
 *
 * Copyright 2026 by Corey Pennycuff
 */

#include "test_helpers.h"

#include <ghoti.io/shape-geometry/libver.h>

#include <cstring>

TEST(Version, StringMatchesTheGeneratedHeader) {
  EXPECT_STREQ(gshgeom_version_string(), GSHGEOM_VERSION_STRING);
  EXPECT_EQ(gshgeom_version_number(), GSHGEOM_VERSION_NUMBER);
  EXPECT_EQ(GSHGEOM_VERSION_NUMBER,
      GSHGEOM_MAKE_VERSION(GSHGEOM_VERSION_MAJOR, GSHGEOM_VERSION_MINOR,
          GSHGEOM_VERSION_PATCH));
}

TEST(Version, StringIsNeverNullOrEmpty) {
  ASSERT_NE(gshgeom_version_string(), nullptr);
  EXPECT_GT(std::strlen(gshgeom_version_string()), 0u);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
