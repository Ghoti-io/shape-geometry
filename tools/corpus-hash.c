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

static void feed_int(int v) {
  count++;
  gshgeom_hash_u32(&total, (uint32_t)v);
}

static void feed_map(const GSHGEOM_Map * m) {
  count++;
  gshgeom_hash_u64(&total, gshgeom_map_hash(m));
  for (int i = 0; i < m->dim; i++) {
    for (int j = 0; j < m->dim; j++) {
      feed(m->m[i][j]);
    }
    feed(m->t[i]);
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


  /* Maps (story 6): the constructors of every angle unit, the composer and the
   * chain, application, the route predicate, and the exact predicates.  The
   * matrices are built from the generator's grid so that they are exact in
   * double and every product rounds the same way on every target; what the
   * probe can see is a fused multiply-add in the composer or in the filter's
   * floating-point part, a different association, and gmath's sine and
   * cosine (the angles are in every unit). */
  for (int dim = 2; dim <= 3; dim++) {
    for (int i = 0; i < 60; i++) {
      GSHGEOM_Map a, b, c, d, e;
      double ma[9], mb[9], mc[9], md[9], ta[3], tb[3], tc[3], td[3];
      for (int k = 0; k < 9; k++) {
        ma[k] = coord() / 8.0;
        mb[k] = coord() / 8.0;
        mc[k] = coord() / 8.0;
        md[k] = coord() / 8.0;
      }
      for (int k = 0; k < 3; k++) {
        ta[k] = coord();
        tb[k] = coord();
        tc[k] = coord();
        td[k] = coord();
      }
      if (gshgeom_map_from_matrix(dim, ma, ta, &a) != GSHGEOM_OK ||
          gshgeom_map_from_matrix(dim, mb, tb, &b) != GSHGEOM_OK ||
          gshgeom_map_from_matrix(dim, mc, tc, &c) != GSHGEOM_OK ||
          gshgeom_map_from_matrix(dim, md, td, &d) != GSHGEOM_OK) {
        return fail("a map was refused");
      }
      GSHGEOM_Map chain[3];
      chain[0] = a;
      chain[1] = b;
      chain[2] = c;
      if (gshgeom_map_compose(&a, &b, &e) != GSHGEOM_OK) {
        return fail("compose failed");
      }
      feed_map(&e);
      if (gshgeom_map_compose_chain(chain, 3, &d, &e) != GSHGEOM_OK) {
        return fail("compose_chain failed");
      }
      feed_map(&e);
      GSHGEOM_Point p = random_point(dim);
      GSHGEOM_Point q;
      if (gshgeom_map_apply_point(&e, p, &q) != GSHGEOM_OK) {
        return fail("apply failed");
      }
      feed_point(q, dim);
      if (gshgeom_map_apply_vector(&e, p, &q) != GSHGEOM_OK) {
        return fail("apply_vector failed");
      }
      feed_point(q, dim);
    }
  }
  {
    static const GSHGEOM_Angle_Unit units[4] = {GSHGEOM_ANGLE_RADIANS,
        GSHGEOM_ANGLE_DEGREES, GSHGEOM_ANGLE_GRADIANS, GSHGEOM_ANGLE_TURNS};
    static const double quarters[12] = {0.0, 90.0, 180.0, 270.0, -90.0, 360.0,
        450.0, 100.0, 0.25, -0.75, 30.0, 1.0};
    for (int u = 0; u < 4; u++) {
      for (int i = 0; i < 12; i++) {
        GSHGEOM_Map m;
        if (gshgeom_map_rotate(units[u], quarters[i], &m) != GSHGEOM_OK) {
          return fail("a quarter-turn rotation was refused");
        }
        feed_map(&m);
      }
      for (int i = 0; i < 80; i++) {
        double angle = coord() * 12.0;
        GSHGEOM_Map m;
        if (gshgeom_map_rotate(units[u], angle, &m) != GSHGEOM_OK) {
          return fail("a rotation was refused");
        }
        feed_map(&m);
        GSHGEOM_Map m3;
        if (gshgeom_map_rotate_axis((GSHGEOM_Axis)(i % 3), units[u], angle, &m3) != GSHGEOM_OK) {
          return fail("an axis rotation was refused");
        }
        feed_map(&m3);
      }
    }
    for (int i = 0; i < 40; i++) {
      double s = weight();
      double tx = coord();
      double ty = coord();
      GSHGEOM_Map dm, rot, comp;
      if (gshgeom_map_design_to_model(s, tx, ty, &dm) != GSHGEOM_OK ||
          gshgeom_map_rotate(GSHGEOM_ANGLE_DEGREES, (i % 3) ? 0.0 : 90.0, &rot) != GSHGEOM_OK ||
          gshgeom_map_compose(&dm, &rot, &comp) != GSHGEOM_OK) {
        return fail("a design map failed");
      }
      int aligned = 0;
      GSHGEOM_Grid g;
      if (gshgeom_map_grid_aligned(&dm, &aligned, &g) != GSHGEOM_OK || !aligned) {
        return fail("the design map is not grid aligned");
      }
      feed_int(aligned);
      feed(g.scale);
      feed_int(g.flip_y);
      if (gshgeom_map_grid_aligned(&comp, &aligned, NULL) != GSHGEOM_OK) {
        return fail("the route predicate failed");
      }
      feed_int(aligned);
    }
    /* Exact predicates: exactly and nearly degenerate triples and quadruples. */
    for (int i = 0; i < 400; i++) {
      GSHGEOM_Point a = random_point(3);
      GSHGEOM_Point b = random_point(3);
      GSHGEOM_Point c = random_point(3);
      int k = (int)(next() % 17) - 8;
      int j = (int)(next() % 9) - 4;
      int perturb = (int)(next() % 5) - 2;
      GSHGEOM_Point on;
      on.x = a.x + ((double)k / 8.0) * (b.x - a.x);
      on.y = a.y + ((double)k / 8.0) * (b.y - a.y);
      on.z = a.z + ((double)k / 8.0) * (b.z - a.z) + ((double)j / 8.0) * (c.z - a.z);
      GSHGEOM_Point on3 = on;
      on3.x += ((double)j / 8.0) * (c.x - a.x);
      on3.y += ((double)j / 8.0) * (c.y - a.y);
      on.z = a.z + ((double)k / 8.0) * (b.z - a.z);
      GSHGEOM_Point near2 = on;
      near2.z = 0.0;
      near2.x += (double)perturb * 0x1p-40;
      GSHGEOM_Point near3 = on3;
      near3.z += (double)perturb * 0x1p-40;
      int s2 = 9, s3 = 9, s2x = 9, s3x = 9;
      if (gshgeom_orient2d(a, b, near2, &s2) != GSHGEOM_OK ||
          gshgeom_orient2d_exact(a, b, near2, &s2x) != GSHGEOM_OK ||
          gshgeom_orient3d(a, b, c, near3, &s3) != GSHGEOM_OK ||
          gshgeom_orient3d_exact(a, b, c, near3, &s3x) != GSHGEOM_OK) {
        return fail("a predicate was refused");
      }
      if (s2 != s2x || s3 != s3x) {
        return fail("the filtered predicate disagrees with the exact one");
      }
      feed_int(s2);
      feed_int(s3);
    }
    /* The smallest triangle, and one at the top of the range. */
    {
      GSHGEOM_Point o = {0.0, 0.0, 0.0};
      GSHGEOM_Point h1 = {0x1p-1074, 0.0, 0.0};
      GSHGEOM_Point h2 = {0.0, 0x1p-1074, 0.0};
      GSHGEOM_Point h3 = {0.0, 0.0, 0x1p-1074};
      GSHGEOM_Point big1 = {0x1p1000, 0.0, 0.0};
      GSHGEOM_Point big2 = {0.0, 0x1p1000, 0.0};
      GSHGEOM_Point big3 = {0.0, 0.0, 0x1p1000};
      int s;
      if (gshgeom_orient2d(o, h1, h2, &s) != GSHGEOM_OK) return fail("orient2d failed");
      feed_int(s);
      if (gshgeom_orient3d(o, h1, h2, h3, &s) != GSHGEOM_OK) return fail("orient3d failed");
      feed_int(s);
      if (gshgeom_orient3d(o, big1, big2, big3, &s) != GSHGEOM_OK) return fail("orient3d failed");
      feed_int(s);
    }
  }

  printf("shape-geometry.hash %016" PRIx64 "\n", gshgeom_hash_value(&total));
  printf("shape-geometry.count %lu\n", count);
  return 0;
}
