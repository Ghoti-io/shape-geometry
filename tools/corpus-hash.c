/*
 * SPDX-License-Identifier: LGPL-3.0-only
 *
 * Copyright (C) 2026 Corey Pennycuff
 *
 * This file is part of Ghoti.io Shape Geometry.
 *
 * Ghoti.io Shape Geometry is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3 as
 * published by the Free Software Foundation.
 *
 * Ghoti.io Shape Geometry is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * @file
 *
 * The cross-target probe (CAP-2, CAP-13, AD-17): a fixed set of curves and
 * paths, evaluated, bounded and hashed in the library's own spelling, printed
 * as one hash.  Every target must print the committed value in
 * tests/corpus/geometry.hash: x86-64, i686 (SSE), arm64, big-endian and
 * Windows.  It is built by `make check-corpus-hash` on the host and by
 * suite/tools/xarch/shape.sh and suite/tools/xwin/shape.sh for the others.
 *
 * What it can see: x87 excess precision and fused multiply-adds (de
 * Casteljau's s*a + t*b is exactly the shape a compiler fuses), byte order
 * (through the hash), and anything in the root isolation that depends on the
 * order of operations.  The inputs are generated from an integer generator and
 * are exact in double, so a difference cannot come from parsing a decimal.
 *
 * -DSHAPE_GEOM_PLANT_BIT flips one bit of one result before it is hashed: the
 * planted differing bit that the gate has to catch on the target it is built
 * for.
 */

#include <ghoti.io/shape-geometry/shape-geometry.h>

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>

static GSHGEOM_Hash total;
static unsigned long count;

static uint64_t rng_state = UINT64_C(0x9E3779B97F4A7C15);

/* xorshift64*: the same integers everywhere. */
static uint64_t next(void) {
  rng_state ^= rng_state >> 12;
  rng_state ^= rng_state << 25;
  rng_state ^= rng_state >> 27;
  return rng_state * UINT64_C(2685821657736338717);
}

/* A coordinate in [-16, 16) with 2^-19 granularity: exact in double. */
static double coord(void) {
  return (double)((int64_t)(next() >> 40) - (int64_t)(1 << 23)) / 524288.0;
}

/* A weight in [0.25, 4.25) with 2^-12 granularity: exact in double. */
static double weight(void) {
  return 0.25 + (double)(next() >> 52) / 1024.0;
}

static void feed(double v) {
  uint64_t bits = gshgeom_f64_bits(v);
#ifdef SHAPE_GEOM_PLANT_BIT
  if (count == 1000) {
    bits ^= 1u;
  }
#endif
  count++;
  gshgeom_hash_u64(&total, bits);
}

static void feed_point(GSHGEOM_Point p, int dim) {
  feed(p.x);
  feed(p.y);
  if (dim == 3) {
    feed(p.z);
  }
}

static int fail(const char * what) {
  fprintf(stderr, "corpus-hash: %s\n", what);
  return 1;
}

static GSHGEOM_Point random_point(int dim) {
  /* One statement each: the order of the elements of an initializer list is
   * unspecified in C, as is the order of the arguments of a call. */
  GSHGEOM_Point p = {0.0, 0.0, 0.0};
  p.x = coord();
  p.y = coord();
  if (dim == 3) {
    p.z = coord();
  }
  return p;
}

int main(void) {
  gshgeom_hash_init(&total);

  for (int dim = 2; dim <= 3; dim++) {
    for (int kind = 0; kind < GSHGEOM_CURVE_KIND_COUNT; kind++) {
      for (int i = 0; i < 300; i++) {
        GSHGEOM_Point p[4];
        double w[4];
        for (int k = 0; k < 4; k++) {
          p[k] = random_point(dim);
          w[k] = weight();
        }
        GSHGEOM_Curve c;
        GSHGEOM_Result r = GSHGEOM_ERR_INTERNAL;
        switch (kind) {
        case GSHGEOM_CURVE_LINE:
          r = gshgeom_curve_line(dim, p[0], p[1], &c);
          break;
        case GSHGEOM_CURVE_QUADRATIC:
          r = gshgeom_curve_quadratic(dim, p[0], p[1], p[2], &c);
          break;
        case GSHGEOM_CURVE_CUBIC:
          r = gshgeom_curve_cubic(dim, p[0], p[1], p[2], p[3], &c);
          break;
        case GSHGEOM_CURVE_RATIONAL_QUADRATIC:
          r = gshgeom_curve_rational_quadratic(dim, p[0], p[1], p[2], w, &c);
          break;
        default:
          r = gshgeom_curve_rational_cubic(dim, p[0], p[1], p[2], p[3], w, &c);
          break;
        }
        if (r != GSHGEOM_OK) {
          return fail("a curve was refused");
        }
        gshgeom_hash_u32(&total, (uint32_t)kind);
        GSHGEOM_Bounds b;
        if (gshgeom_curve_bounds(&c, &b) != GSHGEOM_OK) {
          return fail("bounds failed");
        }
        feed_point(b.min, dim);
        feed_point(b.max, dim);
        for (int k = 0; k <= 16; k++) {
          GSHGEOM_Point at;
          if (gshgeom_curve_eval(&c, (double)k / 16.0, &at) != GSHGEOM_OK) {
            return fail("eval failed");
          }
          feed_point(at, dim);
        }
      }
    }
  }

  /* Paths of every verb, hashed whole and bounded. */
  for (int dim = 2; dim <= 3; dim++) {
    for (int i = 0; i < 40; i++) {
      GSHGEOM_Path * path = NULL;
      if (gshgeom_path_create(NULL, NULL, dim, &path) != GSHGEOM_OK) {
        return fail("path create failed");
      }
      GSHGEOM_Result r = gshgeom_path_move_to(path, random_point(dim));
      int steps = 4 + (int)(next() % 12);
      for (int s = 0; s < steps && r == GSHGEOM_OK; s++) {
        double w3[3];
        double w4[4];
        for (int k = 0; k < 4; k++) {
          w4[k] = weight();
        }
        for (int k = 0; k < 3; k++) {
          w3[k] = weight();
        }
        /* The arguments are drawn into variables first: the order in which a
         * C compiler evaluates the arguments of a call is unspecified (GCC
         * goes right to left on x86 and left to right elsewhere), and the
         * first version of this probe printed two different hashes for it. */
        GSHGEOM_Point q1 = random_point(dim);
        GSHGEOM_Point q2 = random_point(dim);
        GSHGEOM_Point q3 = random_point(dim);
        switch (next() % 7) {
        case 0:
          r = gshgeom_path_line_to(path, q1);
          break;
        case 1:
          r = gshgeom_path_quadratic_to(path, q1, q2);
          break;
        case 2:
          r = gshgeom_path_cubic_to(path, q1, q2, q3);
          break;
        case 3:
          r = gshgeom_path_rational_quadratic_to(path, q1, q2, w3);
          break;
        case 4:
          r = gshgeom_path_rational_cubic_to(path, q1, q2, q3, w4);
          break;
        case 5:
          r = gshgeom_path_close(path);
          break;
        default:
          r = gshgeom_path_move_to(path, q1);
          break;
        }
      }
      GSHGEOM_Bounds b;
      if (r != GSHGEOM_OK || gshgeom_path_bounds(path, &b) != GSHGEOM_OK) {
        gshgeom_path_destroy(path);
        return fail("a path was refused");
      }
      gshgeom_hash_u64(&total, gshgeom_path_hash(path));
      feed_point(b.min, dim);
      feed_point(b.max, dim);
      gshgeom_path_destroy(path);
    }
  }

  printf("shape-geometry.hash %016" PRIx64 "\n", gshgeom_hash_value(&total));
  printf("shape-geometry.count %lu\n", count);
  return 0;
}
