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
#include <string.h>

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
#ifdef SHAPE_GEOM_PLANT_BIT_INT
  /* A late result of story 7's rounding, device path or flattener. */
  static int flipped;
  if (count >= 160000 && !flipped) {
    flipped = 1;
    v ^= 1;
  }
#endif
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

  /* Story 7: the rounding rule, the device path and the integer flattener. */
  {
    /* The rule, on values around the half-way points and the range limit. */
    for (int i = 0; i < 6000; i++) {
      double x;
      switch (i % 3) {
      case 0:
        x = coord() * (double)(1 << (i % 24));
        break;
      case 1:
        x = ((double)((int64_t)(next() >> 32) - ((int64_t)1 << 31)) + 0.5) / 64.0;
        break;
      default:
        x = ((double)(int64_t)(next() >> 40) - 8388608.0) / 64.0 + 1.0 / 128.0;
        break;
      }
      int32_t v = 0;
      GSHGEOM_Result r = gshgeom_round_26_6(x, &v);
      feed_int((int)r);
      if (r == GSHGEOM_OK) {
        feed_int((int)v);
      }
    }
    /* The edges: half-way points, the range limit and the neighbours of each,
     * one unit in the last place either side. */
    {
      static const double edges[] = {2147483647.0, 2147483647.5, 2147483646.5,
          -2147483648.0, -2147483648.5, -2147483649.0, 0.5, -0.5, 1.5, -1.5,
          0.49999999999999994, -0.49999999999999994, 4294967296.0, 0.0};
      for (unsigned e = 0; e < sizeof edges / sizeof edges[0]; e++) {
        for (int d = -1; d <= 1; d++) {
          double x = edges[e] / 64.0;
          /* One unit in the last place of the magnitude, either sign: the
           * sign bit is separate, so adding to the pattern changes |x|. */
          uint64_t bits = gshgeom_f64_bits(x);
          if (x == 0.0) {
            bits = (d < 0) ? UINT64_C(0x8000000000000001) : (d > 0 ? 1u : 0u);
          } else {
            bits += (uint64_t)(int64_t)d;
          }
          double y;
          memcpy(&y, &bits, sizeof y);
          int32_t v = 0;
          GSHGEOM_Result r = gshgeom_round_26_6(y, &v);
          feed_int((int)r);
          if (r == GSHGEOM_OK) {
            feed_int((int)v);
          }
        }
      }
    }
    /* Random 2D paths of every verb under random maps, mapped, rounded and
     * flattened: the device path's and the flattened path's hashes. */
    for (int i = 0; i < 300; i++) {
      GSHGEOM_Path * path = NULL;
      if (gshgeom_path_create(NULL, NULL, 2, &path) != GSHGEOM_OK) {
        return fail("path create failed");
      }
      GSHGEOM_Point m0 = random_point(2);
      GSHGEOM_Result r = gshgeom_path_move_to(path, m0);
      int steps = 3 + (int)(next() % 8);
      for (int s = 0; s < steps && r == GSHGEOM_OK; s++) {
        double w3[3];
        double w4[4];
        for (int k = 0; k < 4; k++) {
          w4[k] = weight();
        }
        for (int k = 0; k < 3; k++) {
          w3[k] = weight();
        }
        GSHGEOM_Point q1 = random_point(2);
        GSHGEOM_Point q2 = random_point(2);
        GSHGEOM_Point q3 = random_point(2);
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
      if (r != GSHGEOM_OK) {
        return fail("a path step was refused");
      }
      GSHGEOM_Map map;
      switch (i % 4) {
      case 0: {
        /* Drawn into variables: the order of a call's arguments is unspecified
         * (the first run of this section printed one hash on x86 and another on
         * every other target). */
        double sx = weight() * 16.0;
        double sy = weight() * 16.0;
        r = gshgeom_map_scale(2, sx, sy, 1.0, &map);
        break;
      }
      case 1: {
        GSHGEOM_Map rot, sc;
        double angle = coord() * 20.0;
        double sx = weight() * 8.0;
        double sy = weight() * 8.0;
        r = gshgeom_map_rotate(GSHGEOM_ANGLE_DEGREES, angle, &rot);
        if (r == GSHGEOM_OK) {
          r = gshgeom_map_scale(2, sx, sy, 1.0, &sc);
        }
        if (r == GSHGEOM_OK) {
          r = gshgeom_map_compose(&sc, &rot, &map);
        }
        break;
      }
      case 2: {
        double scale = weight() * 12.0;
        double tx = coord();
        double ty = coord();
        r = gshgeom_map_design_to_model(scale, tx, ty, &map);
        break;
      }
      default:
        r = gshgeom_map_scale(2, 0x1p20, 0x1p20, 1.0, &map);
        break;
      }
      if (r != GSHGEOM_OK) {
        return fail("a map was refused");
      }
      GSHGEOM_Device_Path * dp = NULL;
      r = gshgeom_device_path_from_model(path, &map, NULL, NULL, &dp);
      feed_int((int)r);
      if (r == GSHGEOM_OK) {
        count++;
        gshgeom_hash_u64(&total, gshgeom_device_path_hash(dp));
        GSHGEOM_Flat_Path * fp = NULL;
        r = gshgeom_device_path_flatten(dp, NULL, NULL, &fp);
        feed_int((int)r);
        if (r == GSHGEOM_OK) {
          count++;
          gshgeom_hash_u64(&total, gshgeom_flat_path_hash(fp));
          gshgeom_flat_path_destroy(fp);
        }
        gshgeom_device_path_destroy(dp);
      }
      gshgeom_path_destroy(path);
    }
    /* The flattener on integer curves across scales, as font would see them. */
    for (int i = 0; i < 400; i++) {
      GSHGEOM_Device_Path * dp = NULL;
      if (gshgeom_device_path_create(NULL, NULL, &dp) != GSHGEOM_OK) {
        return fail("device path create failed");
      }
      int bits = 4 + (int)(next() % 22);
      int64_t span = (int64_t)1 << (bits + 1);
      GSHGEOM_Device_Point q[4];
      for (int k = 0; k < 4; k++) {
        q[k].x = (int32_t)((int64_t)(next() % (uint64_t)span) - ((int64_t)1 << bits));
        q[k].y = (int32_t)((int64_t)(next() % (uint64_t)span) - ((int64_t)1 << bits));
      }
      GSHGEOM_Result r = gshgeom_device_path_move_to(dp, q[0]);
      if (r == GSHGEOM_OK) {
        r = (i & 1) ? gshgeom_device_path_quadratic_to(dp, q[1], q[2])
                    : gshgeom_device_path_cubic_to(dp, q[1], q[2], q[3]);
      }
      GSHGEOM_Flat_Path * fp = NULL;
      if (r == GSHGEOM_OK) {
        r = gshgeom_device_path_flatten(dp, NULL, NULL, &fp);
      }
      if (r != GSHGEOM_OK) {
        return fail("a curve was not flattened");
      }
      count++;
      gshgeom_hash_u64(&total, gshgeom_flat_path_hash(fp));
      gshgeom_flat_path_destroy(fp);
      gshgeom_device_path_destroy(dp);
    }
  }

  /* Story 8: measured paths.  Random paths of every kind, 2D and 3D, measured at
   * the default tolerance and at two others, and queried at 13 distances each:
   * the table's hash, the length, and every number of every frame.  The adaptive
   * rule's decisions are comparisons of sums of products, so a fused
   * multiply-add or an excess of precision changes which panels are accepted and
   * the whole table with them. */
  {
    static const double tolerances[3] = {0.0, 0x1p-20, 0x1p-34};
    for (int i = 0; i < 90; i++) {
      int dim = (i % 3 == 2) ? 3 : 2;
      GSHGEOM_Path * path = NULL;
      if (gshgeom_path_create(NULL, NULL, dim, &path) != GSHGEOM_OK) {
        return fail("path create failed");
      }
      GSHGEOM_Point start = random_point(dim);
      GSHGEOM_Result r = gshgeom_path_move_to(path, start);
      int steps = 1 + (int)(next() % 5);
      for (int k = 0; k < steps && r == GSHGEOM_OK; k++) {
        double w3[3], w4[4];
        for (int j = 0; j < 3; j++) {
          w3[j] = weight();
        }
        for (int j = 0; j < 4; j++) {
          w4[j] = weight();
        }
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
      if (r != GSHGEOM_OK) {
        return fail("a path step was refused");
      }
      GSHGEOM_Measured * m = NULL;
      r = gshgeom_measured_create(NULL, NULL, path, tolerances[i % 3], &m);
      if (r != GSHGEOM_OK) {
        return fail("a path was not measured");
      }
      count++;
      gshgeom_hash_u64(&total, gshgeom_measured_hash(m));
      double length = gshgeom_measured_length(m);
      feed(length);
      feed_int((int)gshgeom_measured_row_count(m));
      for (int k = 0; k <= 12; k++) {
        GSHGEOM_Measured_Frame f;
        r = gshgeom_measured_at_clamped(m, length * (double)k / 12.0, &f);
        feed_int((int)r);
        if (r != GSHGEOM_OK) {
          continue;
        }
        feed_point(f.point, dim);
        feed_point(f.tangent, dim);
        feed_point(f.normal, dim);
        feed(f.t);
        feed_int((int)f.element);
        feed_int(f.has_tangent);
      }
      gshgeom_measured_destroy(m);
      gshgeom_path_destroy(path);
    }
  }

  /* Strokes (stroke.h, SG18): random 2D paths of every kind, stroked with every
   * join and cap, both width spaces and every alignment, under a map that
   * rotates and scales unequally, at the default tolerance and at 2: each
   * outline's hash, which covers every vertex of every part.  The expansion's
   * decisions are exact predicates and its constructions sums of products, the
   * arc length of the arcs join math's trigonometry, so a fused multiply-add, an
   * excess of precision or a libm call changes vertices and with them the hash. */
  {
    GSHGEOM_Map rot, sc, m;
    if (gshgeom_map_rotate(GSHGEOM_ANGLE_DEGREES, 17.0, &rot) != GSHGEOM_OK ||
        gshgeom_map_scale(2, 6.0, 4.5, 1.0, &sc) != GSHGEOM_OK || gshgeom_map_compose(&rot, &sc, &m) != GSHGEOM_OK) {
      return fail("stroke map");
    }
    for (int i = 0; i < 60; i++) {
      GSHGEOM_Path * path = NULL;
      if (gshgeom_path_create(NULL, NULL, 2, &path) != GSHGEOM_OK) {
        return fail("path create failed");
      }
      GSHGEOM_Result r = gshgeom_path_move_to(path, random_point(2));
      int steps = 1 + (int)(next() % 4);
      for (int k = 0; k < steps && r == GSHGEOM_OK; k++) {
        double w3[3], w4[4];
        for (int j = 0; j < 3; j++) {
          w3[j] = weight();
        }
        for (int j = 0; j < 4; j++) {
          w4[j] = weight();
        }
        GSHGEOM_Point q1 = random_point(2);
        GSHGEOM_Point q2 = random_point(2);
        GSHGEOM_Point q3 = random_point(2);
        switch (next() % 6) {
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
        default:
          r = gshgeom_path_close(path);
          break;
        }
      }
      if (r != GSHGEOM_OK) {
        return fail("a stroke path step was refused");
      }
      GSHGEOM_Stroke s;
      gshgeom_stroke_init(&s);
      s.width = 0.25 + (double)(next() % 64) / 32.0;
      s.join = (uint32_t)(i % 5);
      s.cap = (uint32_t)((i / 5) % 3);
      s.width_space = (uint32_t)((i / 15) % 2);
      s.alignment = (uint32_t)((i / 30) % 3);
      s.miter_limit = 1.0 + (double)(next() % 16) / 2.0;
      if (s.width_space == GSHGEOM_STROKE_WIDTH_DEVICE) {
        s.width = s.width * 4.0;
      }
      GSHGEOM_Stroke_Outline * o = NULL;
      r = gshgeom_stroke_expand_with_tolerance(path, &s, &m, (i % 4 == 3) ? 2 : GSHGEOM_FLATTEN_TOLERANCE, NULL, NULL,
          &o);
      feed_int((int)r);
      count++;
      if (r == GSHGEOM_OK) {
        gshgeom_hash_u64(&total, gshgeom_stroke_outline_hash(o));
        feed_int((int)gshgeom_flat_path_point_count(gshgeom_stroke_outline_body(o)));
      }
      gshgeom_stroke_outline_destroy(o);
      gshgeom_path_destroy(path);
    }
  }

  printf("shape-geometry.hash %016" PRIx64 "\n", gshgeom_hash_value(&total));
  printf("shape-geometry.count %lu\n", count);
  return 0;
}
