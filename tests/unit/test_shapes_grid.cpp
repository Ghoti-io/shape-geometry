/**
 * @file
 *
 * The shapes (rect, rounded rect, ellipse) and flat paths from points over a grid of special and ordinary arguments (tests/grid_probe.h), against the hash committed in tests/corpus/geometry.hash.
 *
 * SPDX-License-Identifier: LGPL-3.0-only
 *
 * Copyright 2026 by Corey Pennycuff
 */

#include "../test_helpers.h"

#include "../grid_probe.h"

#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <string>

/**
 * The grid itself is written once, in tests/grid_probe.h, in the subset of C that C++ also reads, so that the cross-target probe
 * (tools/corpus-hash.c) and this test run the same text: this test compares what the grid hashes to on the host with the committed
 * lines `shapes-grid-hash` and `shapes-grid-count`, and the cross-target gates compare the probe's output on every target with the same lines.  A
 * check removed or weakened, a comparison changed, a constant moved or an error ignored moves the hash.  The values are the
 * contract's, reviewed against the header for each refusal's code; regenerate the lines (the probe prints them) only for a change that
 * is meant to move what a caller sees.
 */

namespace {

bool committed(const char * key, std::string * hash, unsigned long * count) {
  FILE * f = std::fopen("tests/corpus/geometry.hash", "r");
  if (!f) {
    return false;
  }
  char name[128];
  char value[128];
  const std::string want_hash = std::string(key) + "-hash";
  const std::string want_count = std::string(key) + "-count";
  bool got_hash = false;
  bool got_count = false;
  while (std::fscanf(f, "%127s %127s", name, value) == 2) {
    if (want_hash == name) {
      *hash = value;
      got_hash = true;
    } else if (want_count == name) {
      *count = std::strtoul(value, nullptr, 10);
      got_count = true;
    }
  }
  std::fclose(f);
  return got_hash && got_count;
}

}  // namespace

TEST(ShapesGrid, EveryEntryPointOverTheSpecialAndOrdinaryArgumentsIsWhatItWas) {
  uint64_t words = 0;
  const uint64_t h = grid_shapes_hash(&words);
  char text[32];
  std::snprintf(text, sizeof text, "%016" PRIx64, h);
  std::string want_hash;
  unsigned long want_count = 0;
  ASSERT_TRUE(committed("shapes-grid", &want_hash, &want_count)) << "no shapes-grid lines in tests/corpus/geometry.hash (run from the library's root)";
  EXPECT_EQ(std::string(text), want_hash) << "the grid's hash over " << words << " words";
  EXPECT_EQ((unsigned long)words, want_count);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
