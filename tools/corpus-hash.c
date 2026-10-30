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

/* The vertex section's accumulator. */
typedef struct VertexSink {
  GSHGEOM_Hash total;
  unsigned long count;
} VertexSink;

/* A coordinate for the vertex section: a multiple of 2 on a small grid, or the finer ones of coord(). */
static GSHGEOM_Point vertex_point(unsigned fine) {
  GSHGEOM_Point p;
  if (fine) {
    p.x = coord();
    p.y = coord();
  } else {
    p.x = (double)((int)(next() % 5u) - 2) * 2.0;
    p.y = (double)((int)(next() % 5u) - 2) * 2.0;
  }
  p.z = 0.0;
  return p;
}

static void vertex_f64(VertexSink * vs, double v) {
  uint64_t bits = gshgeom_f64_bits(v);
#ifdef SHAPE_GEOM_PLANT_BIT_VERTEX
  if (vs->count == 300) {
    bits ^= 1u;
  }
#endif
  vs->count++;
  gshgeom_hash_u64(&vs->total, bits);
}

static GSHGEOM_Result vertex_feed(const GSHGEOM_Vertex * v, void * ctx) {
  VertexSink * vs = (VertexSink *)ctx;
  gshgeom_hash_u64(&vs->total, (uint64_t)v->element);
  gshgeom_hash_u64(&vs->total, (uint64_t)v->subpath);
  gshgeom_hash_u32(&vs->total, v->flags);
  vertex_f64(vs, v->point.x);
  vertex_f64(vs, v->point.y);
  vertex_f64(vs, v->in.x);
  vertex_f64(vs, v->in.y);
  vertex_f64(vs, v->out.x);
  vertex_f64(vs, v->out.y);
  vertex_f64(vs, v->direction.x);
  vertex_f64(vs, v->direction.y);
  return GSHGEOM_OK;
}

/* 2^e for a small integer e, by doubling: exact, and no libm. */
static double pow2(int e) {
  double v = 1.0;
  for (; e > 0; e--) {
    v *= 2.0;
  }
  for (; e < 0; e++) {
    v *= 0.5;
  }
  return v;
}

static void feed_arc_int(int v, GSHGEOM_Hash * h, unsigned long * n) {
  (*n)++;
  gshgeom_hash_u32(h, (uint32_t)v);
}

static void feed_arc_u64(uint64_t v, GSHGEOM_Hash * h, unsigned long * n) {
  (*n)++;
  gshgeom_hash_u64(h, v);
}

/* The end point of the last element of @p path (an arc's authoritative end, a curve's last control point). */
static GSHGEOM_Point path_end(const GSHGEOM_Path * path) {
  const size_t n = gshgeom_path_element_count(path);
  GSHGEOM_Verb verb;
  GSHGEOM_Point p = {0.0, 0.0, 0.0};
  if (n == 0 || gshgeom_path_verb(path, n - 1, &verb) != GSHGEOM_OK) {
    return p;
  }
  if (verb == GSHGEOM_VERB_ARC) {
    GSHGEOM_Arc a;
    if (gshgeom_path_arc(path, n - 1, &a) == GSHGEOM_OK) {
      p = a.end;
    }
    return p;
  }
  GSHGEOM_Element el;
  if (gshgeom_path_element(path, n - 1, &el) == GSHGEOM_OK) {
    p = el.curve.p[gshgeom_curve_point_count(el.curve.kind) - 1];
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

  /* Dashed strokes (stroke.h, SG19 to SG30), appended after the solid ones so that every earlier result is unchanged:
   * the same kind of random paths under five dash patterns (a dot pattern, an odd array and a zero gap among them),
   * offsets either side of zero, both scopes, every join and cap, both width spaces and every alignment.  The cut
   * positions come from along.h's quadrature and the pattern's arithmetic, the dots' directions from the
   * hodograph: a fused multiply-add, an excess of precision or a libm call moves a dash end and with it the hash. */
  {
    static const double pat0[2] = {2.5, 1.5}, pat1[3] = {1.0, 2.0, 3.0}, pat2[2] = {0.0, 1.75},
                        pat3[4] = {4.0, 1.0, 0.5, 1.0}, pat4[2] = {0.75, 0.25};
    static const double * const pats[5] = {pat0, pat1, pat2, pat3, pat4};
    static const size_t counts[5] = {2, 3, 2, 4, 2};
    GSHGEOM_Map rot, sc, m;
    if (gshgeom_map_rotate(GSHGEOM_ANGLE_DEGREES, 17.0, &rot) != GSHGEOM_OK ||
        gshgeom_map_scale(2, 6.0, 4.5, 1.0, &sc) != GSHGEOM_OK || gshgeom_map_compose(&rot, &sc, &m) != GSHGEOM_OK) {
      return fail("dash map");
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
          r = gshgeom_path_move_to(path, q1);
          break;
        default:
          r = gshgeom_path_close(path);
          break;
        }
      }
      if (r != GSHGEOM_OK) {
        return fail("a dashed stroke path step was refused");
      }
      GSHGEOM_Stroke s;
      gshgeom_stroke_init(&s);
      s.width = 0.25 + (double)(next() % 32) / 32.0;
      s.join = (uint32_t)(i % 5);
      s.cap = (uint32_t)((i / 5) % 3);
      s.width_space = (uint32_t)((i / 15) % 2);
      s.alignment = (uint32_t)((i / 30) % 3);
      s.miter_limit = 1.0 + (double)(next() % 16) / 2.0;
      s.dashes = pats[i % 5];
      s.dash_count = counts[i % 5];
      s.dash_offset = (double)((int)(next() % 9) - 3) * 0.75;
      s.dash_scope = (uint32_t)((i / 5) % 2);
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
    /* the DXF lowering of a few patterns, and the positions of a pattern on a few lengths */
    {
      static const double dxf[6] = {0.5, -0.25, 0.0, -0.25, 1.25, -0.125};
      double out[8];
      size_t n = 0;
      double off = 0.0;
      GSHGEOM_Result rr = gshgeom_stroke_dashes_from_dxf(dxf, 6, 3.0, out, 8, &n, &off);
      feed_int((int)rr);
      feed_int((int)n);
      for (size_t k = 0; k < n; k++) {
        feed(out[k]);
      }
      feed(off);
      count++;
    }
  }

  /* Curve outlines (stroke_curves.h, SG41), appended after the dashed strokes so that every earlier result is
   * unchanged: random paths of every kind, solid and under three dash patterns, in both forms, at two tolerances, every
   * join and cap, both width spaces and every alignment, under the same unequal map: each outline's hash covers every
   * control point and weight of every part.  The offsets are Hermite fits from certified bounds, the arcs are
   * square roots, the lowering of the cubic form bisection: all sums of products in a fixed order, so a fused
   * multiply-add or an excess of precision moves a control point and with it the hash. */
  {
    static const double pat0[2] = {2.5, 1.5}, pat1[3] = {1.0, 2.0, 3.0}, pat2[2] = {0.0, 1.75};
    static const double * const pats[3] = {pat0, pat1, pat2};
    static const size_t counts[3] = {2, 3, 2};
    GSHGEOM_Map rot, sc, m;
    if (gshgeom_map_rotate(GSHGEOM_ANGLE_DEGREES, 17.0, &rot) != GSHGEOM_OK ||
        gshgeom_map_scale(2, 6.0, 4.5, 1.0, &sc) != GSHGEOM_OK || gshgeom_map_compose(&rot, &sc, &m) != GSHGEOM_OK) {
      return fail("curve map");
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
        return fail("a curve outline path step was refused");
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
      if (i % 3 != 0) {
        s.dashes = pats[(i % 3) - 1 + (i / 3) % 2];
        s.dash_count = counts[(i % 3) - 1 + (i / 3) % 2];
        s.dash_offset = (double)((int)(next() % 9) - 3) * 0.75;
        s.dash_scope = (uint32_t)((i / 5) % 2);
      }
      GSHGEOM_Stroke_Curves * o = NULL;
      r = gshgeom_stroke_curves_expand(path, &s, &m, (i % 4 == 3) ? 0x1p-8 : GSHGEOM_STROKE_CURVES_TOLERANCE,
          (i / 2) % 2 ? GSHGEOM_STROKE_CURVES_CUBIC : GSHGEOM_STROKE_CURVES_RATIONAL, NULL, NULL, &o);
      feed_int((int)r);
      count++;
      if (r == GSHGEOM_OK) {
        gshgeom_hash_u64(&total, gshgeom_stroke_curves_hash(o));
        feed_int((int)gshgeom_path_element_count(gshgeom_stroke_curves_body(o)));
      }
      gshgeom_stroke_curves_destroy(o);
      gshgeom_path_destroy(path);
    }
  }

  printf("shape-geometry.hash %016" PRIx64 "\n", gshgeom_hash_value(&total));
  printf("shape-geometry.count %lu\n", count);

  /* The predicate section (masks-images-patterns story I2, IP8): whether a map stretches every direction, over maps
   * that are exact in binary and over the unit scale and its neighbours by a few ulps.  Its own hash, so the first
   * pair above is unchanged.  Pure integer arithmetic inside, so every target must print the same. */
  GSHGEOM_Hash predicate_total;
  unsigned long predicate_count = 0;
  gshgeom_hash_init(&predicate_total);
  for (unsigned i = 0; i < 600; i++) {
    double m[4];
    for (int k = 0; k < 4; k++) {
      unsigned pick = (unsigned)(next() % 5u);
      uint64_t draw = next();
      if (pick == 0) {
        m[k] = 1.0 + (double)((int)(draw % 9u) - 4) * 0x1p-52;
      } else if (pick == 1) {
        m[k] = (double)((int)(draw % 17u) - 8) / 4.0;
      } else if (pick == 2) {
        m[k] = coord();
      } else if (pick == 3) {
        m[k] = coord() / 16.0;
      } else {
        m[k] = (double)((int64_t)(draw % 9u) - 4);
      }
    }
    const double t[2] = {0.0, 0.0};
    GSHGEOM_Map map;
    if (gshgeom_map_from_matrix(2, m, t, &map) != GSHGEOM_OK) {
      continue;
    }
    gshgeom_hash_u64(&predicate_total, (uint64_t)gshgeom_map_stretches_every_direction(&map));
    predicate_count++;
  }
  printf("shape-geometry.predicate-hash %016" PRIx64 "\n", gshgeom_hash_value(&predicate_total));
  printf("shape-geometry.predicate-count %lu\n", predicate_count);

  /* The vertex section (vertices.h): the vertices of random paths and the directions at them, bit for bit. The paths
   * mix coordinates on a coarse grid (so that segments without extent, reversals and closed corners are common) with
   * the fine ones above (so that the unit vectors and the bisectors are irrational, and a square root, a quotient or a
   * contraction that differs between targets shows).  Its own pair, so the first two never move. */
  {
    VertexSink vs;
    gshgeom_hash_init(&vs.total);
    vs.count = 0;
    for (unsigned round = 0; round < 500; round++) {
      GSHGEOM_Path * path = NULL;
      if (gshgeom_path_create(NULL, NULL, 2, &path) != GSHGEOM_OK) {
        return fail("cannot make a path for the vertex section");
      }
      const unsigned fine = round % 2u;
      const unsigned len = 1u + (unsigned)(next() % 12u);
      for (unsigned i = 0; i < len; i++) {
        const unsigned pick = i == 0 ? 0u : (unsigned)(next() % 9u);
        GSHGEOM_Point a = vertex_point(fine), b = vertex_point(fine), c = vertex_point(fine);
        GSHGEOM_Result r = GSHGEOM_OK;
        if (pick == 0) {
          r = gshgeom_path_move_to(path, a);
        } else if (pick <= 3) {
          r = gshgeom_path_line_to(path, a);
        } else if (pick == 4) {
          r = gshgeom_path_quadratic_to(path, a, b);
        } else if (pick == 5 || pick == 6) {
          r = gshgeom_path_cubic_to(path, a, b, c);
        } else if (pick == 7) {
          r = gshgeom_path_close(path);
        } else {
          GSHGEOM_Element last;
          size_t n = gshgeom_path_element_count(path);
          (void)gshgeom_path_element(path, n - 1, &last);
          GSHGEOM_Point q = last.curve.p[gshgeom_curve_point_count(last.curve.kind) - 1];
          r = gshgeom_path_line_to(path, q); /* a line of no extent */
        }
        if (r != GSHGEOM_OK) {
          return fail("cannot build a path for the vertex section");
        }
      }
      if (gshgeom_path_vertices(path, vertex_feed, &vs) != GSHGEOM_OK) {
        return fail("the vertex walk failed");
      }
      gshgeom_path_destroy(path);
    }
    printf("shape-geometry.vertex-hash %016" PRIx64 "\n", gshgeom_hash_value(&vs.total));
    printf("shape-geometry.vertex-count %lu\n", vs.count);
  }

  /* The arc section (arc.h): random arcs from every constructor and their pieces, bit for bit.  The inputs are doubles
   * with short fractions (so that quarter-aligned and near-degenerate cases are common) and fine ones (so that the
   * sines, cosines, arctangents and square roots of math, the exact sums and the rounding of the radicand all show a
   * difference between targets).  Its own pair after the others, so they never move. */
  {
    GSHGEOM_Hash arc_total;
    unsigned long arc_count = 0;
    gshgeom_hash_init(&arc_total);
    static const GSHGEOM_Angle_Unit units[4] = {
        GSHGEOM_ANGLE_RADIANS, GSHGEOM_ANGLE_DEGREES, GSHGEOM_ANGLE_GRADIANS, GSHGEOM_ANGLE_TURNS};
    static const double phis[8] = {0.0, 90.0, 45.5, -30.0, 180.0, 360.0, 720.0, 12.25};
    for (unsigned round = 0; round < 500; round++) {
      GSHGEOM_Path * path = NULL;
      if (gshgeom_path_create(NULL, NULL, 2, &path) != GSHGEOM_OK) {
        return fail("cannot make a path for the arc section");
      }
      const unsigned fine = round % 2u;
      GSHGEOM_Point s = vertex_point(fine);
      if (gshgeom_path_move_to(path, s) != GSHGEOM_OK) {
        return fail("cannot begin a path for the arc section");
      }
      for (unsigned step = 0; step < 4; step++) {
        const unsigned form = (unsigned)(next() % 6u);
        GSHGEOM_Result r = GSHGEOM_OK;
        GSHGEOM_Point e = vertex_point(fine);
        if (form <= 2) {
          const double rx = fine ? 0.25 + (double)(next() % 4096u) / 256.0 : (double)(1u + next() % 12u);
          const double ry = (next() % 3u == 0u) ? rx : (fine ? 0.25 + (double)(next() % 4096u) / 256.0 : (double)(1u + next() % 12u));
          GSHGEOM_Arc_Kind kind = GSHGEOM_ARC_KIND_NOTHING;
          r = gshgeom_path_arc_endpoint_to(path, e, rx, ry, phis[next() % 8u], GSHGEOM_ANGLE_DEGREES,
              (int)(next() % 2u), (int)(next() % 2u), &kind);
          feed_arc_int((int)kind, &arc_total, &arc_count);
        } else if (form == 3) {
          const double t0 = (double)((int)(next() % 721u) - 360);
          const double sw = (double)((int)(next() % 719u) - 359) + (fine ? (double)(next() % 8u) / 8.0 : 0.0);
          const unsigned u = (unsigned)(next() % 4u);
          const double scale = u == 1u ? 1.0 : u == 2u ? 400.0 / 360.0 : u == 3u ? 1.0 / 360.0 : 3.14159265358979323846 / 180.0;
          /* the centre form begins its own subpath (no current point), so it goes into a path of its own */
          GSHGEOM_Path * own = NULL;
          if (gshgeom_path_create(NULL, NULL, 2, &own) != GSHGEOM_OK) {
            return fail("cannot make a path for a centre-form arc");
          }
          r = gshgeom_path_arc_center_to(own, e, 2.0 + (double)(next() % 9u), 1.0 + (double)(next() % 9u),
              phis[next() % 8u] * scale, units[u], t0 * scale, sw == 0.0 ? scale : sw * scale, NULL);
          feed_arc_u64(gshgeom_path_hash(own), &arc_total, &arc_count);
          gshgeom_path_destroy(own);
          r = r == GSHGEOM_OK ? GSHGEOM_OK : r;
        } else if (form == 4) {
          r = gshgeom_path_arc_tangent_to(path, e, vertex_point(fine), (double)(next() % 5u) * 0.75);
        } else {
          r = gshgeom_path_arc_bulge_to(path, e, (double)((int)(next() % 41u) - 20) / 8.0);
        }
        feed_arc_int((int)r, &arc_total, &arc_count);
      }
      feed_arc_u64(gshgeom_path_hash(path), &arc_total, &arc_count);
      GSHGEOM_Path * lowered = NULL;
      if (gshgeom_path_lower_arcs(path, NULL, NULL, &lowered) != GSHGEOM_OK) {
        return fail("cannot lower the paths of the arc section");
      }
      feed_arc_u64(gshgeom_path_hash(lowered), &arc_total, &arc_count);
      GSHGEOM_Bounds b;
      if (gshgeom_path_bounds(path, &b) != GSHGEOM_OK) {
        return fail("the bounds of an arc path failed");
      }
      feed_arc_u64(gshgeom_f64_bits(b.min.x), &arc_total, &arc_count);
      feed_arc_u64(gshgeom_f64_bits(b.max.y), &arc_total, &arc_count);
      gshgeom_path_destroy(lowered);
      gshgeom_path_destroy(path);
    }
    printf("shape-geometry.arc-hash %016" PRIx64 "\n", gshgeom_hash_value(&arc_total));
    printf("shape-geometry.arc-count %lu\n", arc_count);
  }

  /* The arc-vertex section (vertices.h, AR14): the vertices of random paths that hold arcs, and the directions at them,
   * bit for bit.  Lines, arcs of every constructor (a circle or ellipse that begins at the current point, SVG's endpoint
   * form, arcTo, a bulge), closes and lines of no extent; coordinates on a coarse grid (axis-aligned tangents,
   * reversals, smooth joins) and fine ones.  Its own pair after the others, so they never move. */
  {
    VertexSink vs;
    gshgeom_hash_init(&vs.total);
    vs.count = 0;
    static const double quarter_cos[4] = {1.0, 0.0, -1.0, 0.0};
    static const double quarter_sin[4] = {0.0, 1.0, 0.0, -1.0};
    for (unsigned round = 0; round < 500; round++) {
      GSHGEOM_Path * path = NULL;
      if (gshgeom_path_create(NULL, NULL, 2, &path) != GSHGEOM_OK) {
        return fail("cannot make a path for the arc-vertex section");
      }
      const unsigned fine = round % 2u;
      GSHGEOM_Point cur = vertex_point(fine);
      if (gshgeom_path_move_to(path, cur) != GSHGEOM_OK) {
        return fail("cannot begin a path for the arc-vertex section");
      }
      const unsigned len = 2u + (unsigned)(next() % 8u);
      for (unsigned i = 0; i < len; i++) {
        const unsigned pick = (unsigned)(next() % 9u);
        GSHGEOM_Point e = vertex_point(fine);
        GSHGEOM_Result r = GSHGEOM_OK;
        if (pick == 0) {
          r = gshgeom_path_line_to(path, e);
        } else if (pick <= 2) {
          /* a circle or an ellipse that begins at the current point, at a quarter turn of its parameter */
          const unsigned k = (unsigned)(next() % 4u);
          const double rx = 1.0 + (double)(next() % 6u), ry = (next() % 2u == 0u) ? rx : 1.0 + (double)(next() % 6u);
          const double sweep = (double)(1 + (int)(next() % 8u)) * 45.0 * ((next() % 2u == 0u) ? 1.0 : -1.0);
          GSHGEOM_Arc_Form f;
          f.center.x = cur.x - (rx * quarter_cos[k]);
          f.center.y = cur.y - (ry * quarter_sin[k]);
          f.center.z = 0.0;
          f.u.x = rx;
          f.u.y = 0.0;
          f.u.z = 0.0;
          f.v.x = 0.0;
          f.v.y = ry;
          f.v.z = 0.0;
          f.unit = GSHGEOM_ANGLE_DEGREES;
          f.t0 = (double)k * 90.0;
          f.sweep = sweep;
          f.flags = rx == ry ? GSHGEOM_ARC_CIRCULAR : 0u;
          r = gshgeom_path_arc_conjugate_to(path, &f, NULL);
        } else if (pick <= 4) {
          r = gshgeom_path_arc_endpoint_to(path, e, 1.0 + (double)(next() % 8u), 1.0 + (double)(next() % 8u),
              (double)(next() % 4u) * 30.0, GSHGEOM_ANGLE_DEGREES, (int)(next() % 2u), (int)(next() % 2u), NULL);
        } else if (pick == 5) {
          r = gshgeom_path_arc_tangent_to(path, e, vertex_point(fine), (double)(1u + next() % 4u));
        } else if (pick == 6) {
          r = gshgeom_path_arc_bulge_to(path, e, (double)((int)(next() % 17u) - 8) / 4.0);
        } else if (pick == 7) {
          r = gshgeom_path_close(path);
        } else {
          r = gshgeom_path_line_to(path, cur); /* a line of no extent */
        }
        if (r != GSHGEOM_OK && r != GSHGEOM_ERR_LIMIT && r != GSHGEOM_ERR_INVALID) {
          return fail("cannot build a path for the arc-vertex section");
        }
        cur = path_end(path);
      }
      if (gshgeom_path_vertices(path, vertex_feed, &vs) != GSHGEOM_OK) {
        return fail("the arc-vertex walk failed");
      }
      gshgeom_path_destroy(path);
    }
    printf("shape-geometry.arc-vertex-hash %016" PRIx64 "\n", gshgeom_hash_value(&vs.total));
    printf("shape-geometry.arc-vertex-count %lu\n", vs.count);
  }

  /* The conic section (flatten.h, revision 3; CV5): the device path of random rational curves and arcs, bit for bit, under maps
   * that rotate, skew, mirror, stretch unevenly, collapse (a singular map) and scale from 2^-4 to 2^20, at tolerances 2, 4 and 37.
   * The weights are short fractions (so the halving is exact), powers of two to 2^+-40 and, in a few rounds, the 2^+-100 the path
   * accepts; the arcs are every sweep from a few degrees to a full turn, in the endpoint form and the centre form.  What is hashed is
   * the result code, the element count and the device path's hash: every line of every conic is in it.  Its own pair after the
   * others, so they never move. */
  {
    GSHGEOM_Hash conic_total;
    unsigned long conic_count = 0;
    gshgeom_hash_init(&conic_total);
    static const int32_t tolerances[4] = {2, 4, 4, 37};
    static const double sweeps[8] = {5.0, 45.0, 90.0, 91.0, 180.0, 270.0, -200.0, 360.0};
    for (unsigned round = 0; round < 480; round++) {
      GSHGEOM_Path * path = NULL;
      if (gshgeom_path_create(NULL, NULL, 2, &path) != GSHGEOM_OK) {
        return fail("cannot make a path for the conic section");
      }
      const unsigned wclass = round % 4u;
      GSHGEOM_Point s = vertex_point(1u);
      if (gshgeom_path_move_to(path, s) != GSHGEOM_OK) {
        return fail("cannot begin a path for the conic section");
      }
      GSHGEOM_Result r = GSHGEOM_OK;
      for (unsigned step = 0; step < 3; step++) {
        double w4[4];
        for (unsigned k = 0; k < 4; k++) {
          const uint64_t pick = next();
          if (wclass == 0u) {
            w4[k] = weight();
          } else if (wclass == 1u) {
            w4[k] = pow2((int)(pick % 81u) - 40);
          } else if (wclass == 2u) {
            w4[k] = (round % 24u == 2u) ? ((pick & 1u) ? 0x1p100 : 0x1p-100) : weight();
          } else {
            w4[k] = (pick & 1u) ? 1.0 : 0.5;
          }
        }
        GSHGEOM_Point q1 = random_point(2);
        GSHGEOM_Point q2 = random_point(2);
        GSHGEOM_Point q3 = random_point(2);
        if (step == 0u) {
          r = gshgeom_path_rational_quadratic_to(path, q1, q2, w4);
        } else if (step == 1u) {
          r = gshgeom_path_rational_cubic_to(path, q1, q2, q3, w4);
        } else {
          const double rx = 1.0 + (double)(next() % 9u);
          const double ry = 1.0 + (double)(next() % 9u);
          const double sweep = sweeps[next() % 8u];
          GSHGEOM_Arc_Form f;
          f.center = q1;
          f.u.x = rx;
          f.u.y = 0.0;
          f.u.z = 0.0;
          f.v.x = 0.0;
          f.v.y = ry;
          f.v.z = 0.0;
          f.unit = GSHGEOM_ANGLE_DEGREES;
          f.t0 = (double)(next() % 8u) * 45.0;
          f.sweep = sweep;
          f.flags = rx == ry ? GSHGEOM_ARC_CIRCULAR : 0u;
          GSHGEOM_Path * own = NULL;
          if (gshgeom_path_create(NULL, NULL, 2, &own) != GSHGEOM_OK) {
            return fail("cannot make a path for a centre-form arc of the conic section");
          }
          /* the centre form begins its own subpath: lowered into this path it is a move and its pieces */
          r = gshgeom_path_arc_conjugate_to(own, &f, NULL);
          if (r == GSHGEOM_OK) {
            GSHGEOM_Arc a;
            r = gshgeom_path_arc(own, 1, &a);
            if (r == GSHGEOM_OK) {
              r = gshgeom_path_line_to(path, a.start);
            }
            if (r == GSHGEOM_OK) {
              r = gshgeom_path_append_arc(path, &a);
            }
          }
          gshgeom_path_destroy(own);
        }
        feed_arc_int((int)r, &conic_total, &conic_count);
        if (r != GSHGEOM_OK && r != GSHGEOM_ERR_LIMIT && r != GSHGEOM_ERR_INVALID) {
          return fail("cannot build a path for the conic section");
        }
      }
      if (round % 3u == 0u) {
        (void)gshgeom_path_close(path);
      }
      feed_arc_u64(gshgeom_path_hash(path), &conic_total, &conic_count);
      GSHGEOM_Map map;
      const unsigned mclass = round % 8u;
      const int e = (int)(round / 8u % 7u) * 4 - 4; /* 2^-4 to 2^20 */
      const double scale = pow2(e);
      double m[4] = {scale, 0.0, 0.0, scale};
      double t[2] = {coord() * scale, coord() * scale};
      if (mclass == 1u) {
        m[0] = scale; m[1] = scale * 0.5; m[2] = 0.0; m[3] = scale; /* a skew */
      } else if (mclass == 2u) {
        m[0] = -scale; /* a mirror */
      } else if (mclass == 3u) {
        m[0] = scale * 3.0; m[3] = scale * 0.25; /* uneven */
      } else if (mclass == 4u) {
        m[0] = scale; m[1] = 2.0 * scale; m[2] = 2.0 * scale; m[3] = 4.0 * scale; /* singular */
      } else if (mclass == 5u) {
        m[0] = scale * 0.8660254037844386; m[1] = -scale * 0.5; m[2] = scale * 0.5; m[3] = scale * 0.8660254037844386; /* a rotation */
      }
      if (gshgeom_map_from_matrix(2, m, t, &map) != GSHGEOM_OK) {
        return fail("a map of the conic section was refused");
      }
      GSHGEOM_Device_Path * dp = NULL;
      r = gshgeom_device_path_from_model_with_tolerance(path, &map, tolerances[round % 4u], NULL, NULL, &dp);
      feed_arc_int((int)r, &conic_total, &conic_count);
      if (r == GSHGEOM_OK) {
        feed_arc_u64((uint64_t)gshgeom_device_path_element_count(dp), &conic_total, &conic_count);
        feed_arc_u64(gshgeom_device_path_hash(dp), &conic_total, &conic_count);
        gshgeom_device_path_destroy(dp);
      }
      gshgeom_path_destroy(path);
    }
    printf("shape-geometry.conic-hash %016" PRIx64 "\n", gshgeom_hash_value(&conic_total));
    printf("shape-geometry.conic-count %lu\n", conic_count);
  }

  /* The named primitives (shapes.h; PR2 to PR5, PR13): random boxes, radii and flips through the three builders, the result code and
   * the path's hash (every coordinate, every arc record).  Coordinates are sums that round (a 2^-19 grid plus offsets from 2^-3 to
   * 2^40), a share of the radii are zero (square corners) and a share overlap or are negative, so the refusals are in it too.  Its own
   * pair after the others, so they never move. */
  {
    GSHGEOM_Hash shapes_total;
    unsigned long shapes_count = 0;
    gshgeom_hash_init(&shapes_total);
    for (unsigned round = 0; round < 1200; round++) {
      GSHGEOM_Path * path = NULL;
      if (gshgeom_path_create(NULL, NULL, 2, &path) != GSHGEOM_OK) {
        return fail("cannot make a path for the shapes section");
      }
      const uint64_t pick = next();
      const double scale = (double)(UINT64_C(1) << (pick & 31u)) / 8.0;
      const double x = coord() * scale;
      const double y = coord() * scale;
      double w = (coord() + 16.0) * scale / 32.0;
      double h = (coord() + 16.0) * scale / 32.0;
      if ((pick >> 8) % 23u == 0u) {
        w = 0.0;
      }
      if ((pick >> 12) % 29u == 0u) {
        h = -h;
      }
      const int fx = (int)((pick >> 20) & 1u);
      const int fy = (int)((pick >> 21) & 1u);
      GSHGEOM_Result r;
      const unsigned kind = round % 3u;
      if (kind == 0u) {
        r = gshgeom_path_add_rect(path, x, y, w, h, fx, fy);
      } else if (kind == 1u) {
        GSHGEOM_Corner c[4];
        for (unsigned k = 0; k < 4; k++) {
          const uint64_t q = next();
          c[k].rx = ((q & 7u) == 0u) ? 0.0 : w * (double)((q >> 8) & 255u) / 1024.0;
          c[k].ry = ((q & 7u) == 1u) ? 0.0 : h * (double)((q >> 20) & 255u) / 1024.0;
        }
        r = gshgeom_path_add_rounded_rect(path, x, y, w, h, c);
      } else {
        r = gshgeom_path_add_ellipse(path, x, y, w, h < 0.0 ? -h : h, fx, fy);
      }
      feed_arc_int((int)r, &shapes_total, &shapes_count);
      feed_arc_u64((uint64_t)gshgeom_path_element_count(path), &shapes_total, &shapes_count);
      feed_arc_u64(gshgeom_path_hash(path), &shapes_total, &shapes_count);
      gshgeom_path_destroy(path);
    }
    printf("shape-geometry.shapes-hash %016" PRIx64 "\n", gshgeom_hash_value(&shapes_total));
    printf("shape-geometry.shapes-count %lu\n", shapes_count);
  }
  return 0;
}
