/**
 * @file
 *
 * Shared helpers for the Ghoti.io Shape Geometry unit tests.
 *
 * Copyright 2026 by Corey Pennycuff
 */

#ifndef GHOTI_IO_GSHGEOM_TEST_HELPERS_H
#define GHOTI_IO_GSHGEOM_TEST_HELPERS_H

#include <gtest/gtest.h>

#include <ghoti.io/shape-geometry/shape-geometry.h>

#include <cstdlib>
#include <cstring>

/**
 * True under Valgrind, whose x87 is a double.  The oracles (stroke_oracle.h, curve_eval.h) compute in `long double` and decide
 * membership at ties, so under Valgrind they report disagreements the library does not have; the library holds no `long
 * double`.  A test that compares against them skips there and says why.
 */
inline bool under_valgrind() {
  const char * preload = std::getenv("LD_PRELOAD");
  return std::getenv("VALGRIND_LAUNCHER") != nullptr || (preload && std::strstr(preload, "vgpreload"));
}

#define GSHGEOM_SKIP_UNDER_VALGRIND() \
  do { \
    if (under_valgrind()) { \
      GTEST_SKIP() << "the oracle computes in long double, which Valgrind computes as double"; \
    } \
  } while (0)

#endif /* GHOTI_IO_GSHGEOM_TEST_HELPERS_H */
