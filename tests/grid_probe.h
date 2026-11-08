/**
 * @file
 *
 * The arc, lattice and clip entry points over grids of special and ordinary arguments, folded into one hash each.
 *
 * SPDX-License-Identifier: LGPL-3.0-only
 *
 * Copyright 2026 by Corey Pennycuff
 *
 * Written in the subset of C that C++ also reads, and header-only, so that three programs share one text: the cross-target probe
 * (tools/corpus-hash.c, which every target of the cross gates builds and whose three pairs `arc-grid`, `lattice-grid` and `clip-grid`
 * are committed in tests/corpus/geometry.hash), and the unit tests of the grids (tests/unit/test_arc_grid.cpp and its two
 * siblings), which compare the same hashes with that file.
 *
 * Every entry point is called on every combination of the values that decide a check (zero, one, a negative, a huge, the coordinate
 * limit and the first number past it, NaN, an infinity, a quarter turn and the full turn in each unit, a unit past the last, a flag
 * past 1, NULL arguments, indices past the end, paths of three dimensions, lattices and maps of every kind, random paths clipped to
 * quadrilaterals of every kind) and what each returns, and everything it wrote, is folded into the hash.  A check removed or weakened,
 * a comparison changed, a constant moved or an error ignored changes the answer for some combination, and the hash moves.
 *
 * What the draws of the pseudo-random generator are: one draw to a statement, from an xorshift of 64 bits and the integer-to-double
 * conversions below, which are exact.  The order in which the arguments of a call or the operands of an operator are evaluated is
 * unspecified, so no statement draws twice; no library function of the C library that is not exact (sin, cos, pow) is used here.
 */

#ifndef GHOTI_IO_GSHGEOM_TESTS_GRID_PROBE_H
#define GHOTI_IO_GSHGEOM_TESTS_GRID_PROBE_H

#include <ghoti.io/shape-geometry/shape-geometry.h>

#include <ghoti.io/math/basic.h>

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif
#include "../src/exact/exact.h"
#ifdef __cplusplus
}
#endif

#define GRID_BIG 0x1p100
#define GRID_PAST_BIG 0x1.0000000000001p100
#define GRID_TWO_PI_BELOW 0x1.921fb54442d17p+2
#define GRID_TWO_PI 6.283185307179586
#define GRID_ARRAY_COUNT(a) ((int)(sizeof(a) / sizeof((a)[0])))

/* No function of the C library's math is called here (a cross probe is built without it): the scaling is math's own `gmath_scalbn`, correctly
 * rounded, and the two below are bit moves. */
static inline int grid_isnan(double d) {
  return d != d;
}

static inline double grid_fabs(double d) {
  return d < 0.0 ? -d : d;
}

/** The next double from @p x toward @p toward (both finite, @p x not zero in the uses here). */
static inline double grid_nextafter(double x, double toward) {
  uint64_t b;
  memcpy(&b, &x, 8);
  if (x == toward) {
    return toward;
  }
  if (x == 0.0) {
    b = 1; /* the least subnormal, with the sign of the direction */
    if (toward < 0.0) {
      b |= 0x8000000000000000ull;
    }
  } else if ((x < toward) == (x > 0.0)) {
    b++;
  } else {
    b--;
  }
  memcpy(&x, &b, 8);
  return x;
}

/** Added to every seed: 0 in the pinned grids; the probe takes it from GSHGEOM_GRID_SALT, to try other inputs (the audit of the argued plants). */
static uint64_t grid_seed_salt = 0;

typedef struct Grid_Fold {
  uint64_t h;
  uint64_t n; /* the words folded */
} Grid_Fold;

static inline void grid_fold_init(Grid_Fold * f) {
  f->h = 1469598103934665603ull;
  f->n = 0;
}

static inline void grid_word(Grid_Fold * f, uint64_t w) {
  int i;
  for (i = 0; i < 8; i++) {
    f->h ^= (w >> (8 * i)) & 0xffu;
    f->h *= 1099511628211ull;
  }
  f->n++;
}

static inline void grid_dbl(Grid_Fold * f, double d) {
  uint64_t w;
  memcpy(&w, &d, 8);
  /* all NaNs are one for the hash: the library does not promise a payload */
  grid_word(f, grid_isnan(d) ? 0x7ff8000000000000ull : w);
}

static inline void grid_pt(Grid_Fold * f, const GSHGEOM_Point * p) {
  grid_dbl(f, p->x);
  grid_dbl(f, p->y);
  grid_dbl(f, p->z);
}

static inline void grid_arc(Grid_Fold * f, const GSHGEOM_Arc * a) {
  unsigned unit_bits; /* the field's bytes: a record a failed call left alone may hold a value that is no unit */
  grid_pt(f, &a->start);
  grid_pt(f, &a->end);
  grid_pt(f, &a->center);
  grid_pt(f, &a->u);
  grid_pt(f, &a->v);
  memcpy(&unit_bits, &a->unit, sizeof unit_bits);
  grid_word(f, unit_bits);
  grid_dbl(f, a->t0);
  grid_dbl(f, a->sweep);
  grid_word(f, a->flags);
}

typedef struct Grid_Rng {
  uint64_t s;
} Grid_Rng;

static inline uint64_t grid_next(Grid_Rng * r) {
  r->s ^= r->s << 13;
  r->s ^= r->s >> 7;
  r->s ^= r->s << 17;
  return r->s;
}

/** A draw in [0, 1): the top 53 bits over 2^53, which is exact. */
static inline double grid_unit(Grid_Rng * r) {
  return (double)(grid_next(r) >> 11) / 9007199254740992.0;
}

static inline GSHGEOM_Point grid_P(double x, double y) {
  GSHGEOM_Point p;
  p.x = x;
  p.y = y;
  p.z = 0.0;
  return p;
}

static inline GSHGEOM_Path * grid_make(void) {
  GSHGEOM_Path * p = NULL;
  gshgeom_path_create(NULL, NULL, 2, &p);
  return p;
}

/** One call's record: the result, then the path (its hash and element count) so a write after a failure shows. */
static inline void grid_fold_path(Grid_Fold * f, GSHGEOM_Result r, GSHGEOM_Path * p) {
  grid_word(f, (uint64_t)r);
  grid_word(f, gshgeom_path_hash(p));
  grid_word(f, gshgeom_path_element_count(p));
}

/* ---- arcs ---------------------------------------------------------------------------------------------------------------- */

static inline void grid_arc_endpoint(Grid_Fold * f) {
  static const double ends[] = {0.0, 1.0, 7.5, (double)NAN, (double)INFINITY, GRID_BIG, GRID_PAST_BIG};
  static const double radii[] = {0.0, -0.0, 1.0, 5.0, -1.0, (double)NAN, (double)INFINITY, 0x1p-1000, GRID_PAST_BIG};
  static const double phis[] = {0.0, 90.0, (double)NAN, 1e300};
  static const int units[] = {GSHGEOM_ANGLE_RADIANS, GSHGEOM_ANGLE_DEGREES, GSHGEOM_ANGLE_GRADIANS, GSHGEOM_ANGLE_TURNS,
      (int)GSHGEOM_ANGLE_UNIT_COUNT, -1};
  static const int flags[] = {0, 1, 2, -1};
  int current, ix, iy, irx, iry, iphi, iu, il, is;
  for (current = 0; current < 2; current++) {
    for (ix = 0; ix < GRID_ARRAY_COUNT(ends); ix++) {
      for (iy = 0; iy < GRID_ARRAY_COUNT(ends); iy++) {
        for (irx = 0; irx < GRID_ARRAY_COUNT(radii); irx++) {
          for (iry = 0; iry < GRID_ARRAY_COUNT(radii); iry++) {
            for (iphi = 0; iphi < GRID_ARRAY_COUNT(phis); iphi++) {
              for (iu = 0; iu < GRID_ARRAY_COUNT(units); iu++) {
                for (il = 0; il < GRID_ARRAY_COUNT(flags); il++) {
                  for (is = 0; is < GRID_ARRAY_COUNT(flags); is++) {
                    GSHGEOM_Path * p = grid_make();
                    GSHGEOM_Arc_Kind kind = GSHGEOM_ARC_KIND_ARC; /* a value the call can leave alone: it writes the kind on success only */
                    GSHGEOM_Result r;
                    if (current) {
                      gshgeom_path_move_to(p, grid_P(0.0, 0.0));
                    }
                    r = gshgeom_path_arc_endpoint_to(p, grid_P(ends[ix], ends[iy]), radii[irx], radii[iry], phis[iphi],
                        (GSHGEOM_Angle_Unit)units[iu], flags[il], flags[is], &kind);
                    grid_fold_path(f, r, p);
                    grid_word(f, (uint64_t)kind);
                    gshgeom_path_destroy(p);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

static inline void grid_arc_center(Grid_Fold * f) {
  static const double cs[] = {0.0, (double)NAN, GRID_BIG};
  static const double radii[] = {0.0, 1.0, 5.0, -1.0, (double)NAN, (double)INFINITY, GRID_PAST_BIG};
  static const double phis[] = {0.0, 30.0, (double)NAN};
  static const int units[] = {GSHGEOM_ANGLE_RADIANS, GSHGEOM_ANGLE_DEGREES, GSHGEOM_ANGLE_GRADIANS, GSHGEOM_ANGLE_TURNS,
      (int)GSHGEOM_ANGLE_UNIT_COUNT};
  static const double t0s[] = {0.0, 90.0, (double)NAN, 720.0, 1e300, -0.0, GRID_TWO_PI_BELOW};
  static const double sweeps[] = {0.0, 90.0, 360.0, 361.0, -360.0, (double)NAN, -0.0, GRID_TWO_PI, GRID_TWO_PI_BELOW, 1.0, 400.0};
  GSHGEOM_Point ends[4];
  int current, icx, icy, irx, iry, iphi, iu, it, isw, e;
  ends[0] = grid_P(0.0, 0.0);
  ends[1] = grid_P(1000.0, 0.0);
  ends[2] = grid_P((double)NAN, 0.0);
  ends[3] = grid_P(0.0, (double)INFINITY);
  for (current = 0; current < 2; current++) {
    for (icx = 0; icx < GRID_ARRAY_COUNT(cs); icx++) {
      for (icy = 0; icy < GRID_ARRAY_COUNT(cs); icy++) {
        for (irx = 0; irx < GRID_ARRAY_COUNT(radii); irx++) {
          for (iry = 0; iry < GRID_ARRAY_COUNT(radii); iry++) {
            for (iphi = 0; iphi < GRID_ARRAY_COUNT(phis); iphi++) {
              for (iu = 0; iu < GRID_ARRAY_COUNT(units); iu++) {
                for (it = 0; it < GRID_ARRAY_COUNT(t0s); it++) {
                  for (isw = 0; isw < GRID_ARRAY_COUNT(sweeps); isw++) {
                    for (e = -1; e < 4; e++) {
                      GSHGEOM_Path * p = grid_make();
                      GSHGEOM_Result r;
                      if (current) {
                        gshgeom_path_move_to(p, grid_P(1.0, 1.0));
                      }
                      r = gshgeom_path_arc_center_to(p, grid_P(cs[icx], cs[icy]), radii[irx], radii[iry], phis[iphi],
                          (GSHGEOM_Angle_Unit)units[iu], t0s[it], sweeps[isw], e < 0 ? NULL : &ends[e]);
                      grid_fold_path(f, r, p);
                      gshgeom_path_destroy(p);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

static inline void grid_arc_tangent_bulge(Grid_Fold * f) {
  static const double v[] = {0.0, 1.0, 10.0, (double)NAN, (double)INFINITY, -10.0, GRID_BIG};
  static const double rs[] = {0.0, 1.0, 3.0, (double)NAN, -1.0, (double)INFINITY, 0x1p-1000, GRID_BIG};
  static const double bulges[] = {0.0, 1.0, -1.0, 0.5, (double)NAN, (double)INFINITY, 1e300, 1e-300, -1e-300, 4.0};
  int current, a, b, c, d, r, qx, qy, bi;
  for (current = 0; current < 2; current++) {
    for (a = 0; a < GRID_ARRAY_COUNT(v); a++) {
      for (b = 0; b < GRID_ARRAY_COUNT(v); b++) {
        for (c = 0; c < GRID_ARRAY_COUNT(v); c++) {
          for (d = 0; d < GRID_ARRAY_COUNT(v); d++) {
            for (r = 0; r < GRID_ARRAY_COUNT(rs); r++) {
              GSHGEOM_Path * p = grid_make();
              GSHGEOM_Result res;
              if (current) {
                gshgeom_path_move_to(p, grid_P(0.0, 0.0));
              }
              res = gshgeom_path_arc_tangent_to(p, grid_P(v[a], v[b]), grid_P(v[c], v[d]), rs[r]);
              grid_fold_path(f, res, p);
              gshgeom_path_destroy(p);
            }
          }
        }
      }
    }
    for (qx = 0; qx < GRID_ARRAY_COUNT(v); qx++) {
      for (qy = 0; qy < GRID_ARRAY_COUNT(v); qy++) {
        for (bi = 0; bi < GRID_ARRAY_COUNT(bulges); bi++) {
          GSHGEOM_Path * p = grid_make();
          GSHGEOM_Result res;
          if (current) {
            gshgeom_path_move_to(p, grid_P(0.0, 0.0));
          }
          res = gshgeom_path_arc_bulge_to(p, grid_P(v[qx], v[qy]), bulges[bi]);
          grid_fold_path(f, res, p);
          gshgeom_path_destroy(p);
        }
      }
    }
  }
}

/** A deterministic family of arc records: consistent ones, and each field of one disturbed.  The caller frees the array. */
static inline GSHGEOM_Arc * grid_arc_family(size_t * count_out) {
  static const GSHGEOM_Angle_Unit units[4] = {GSHGEOM_ANGLE_RADIANS, GSHGEOM_ANGLE_DEGREES, GSHGEOM_ANGLE_GRADIANS, GSHGEOM_ANGLE_TURNS};
  static const double turns[4] = {GRID_TWO_PI, 360.0, 400.0, 1.0};
  static const double quarters[4] = {GRID_TWO_PI / 4.0, 90.0, 100.0, 0.25};
  static const double mult[13] = {1e-9, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, -0.25, -1.0, -2.0, -3.5, 4.0 / 3.0};
  static const double t0s[5] = {0.0, 1.0, -0.5, 4.5, 2.0};
  static const double bad[8] = {(double)NAN, (double)INFINITY, -(double)INFINITY, 0.0, GRID_BIG, GRID_PAST_BIG, -GRID_PAST_BIG, 1e-300};
  Grid_Rng rng;
  GSHGEOM_Arc * out = NULL;
  size_t n = 0;
  size_t cap = 0;
  int i, k;
  rng.s = 88172645463325252ull ^ (grid_seed_salt * 0x9e3779b97f4a7c15ull);
#define GRID_PUSH(arcval) \
  do { \
    if (n == cap) { \
      cap = cap ? cap * 2 : 1024; \
      out = (GSHGEOM_Arc *)realloc(out, cap * sizeof *out); \
    } \
    out[n++] = (arcval); \
  } while (0)
  for (i = 0; i < 4; i++) {
    for (k = 0; k < 90; k++) {
      GSHGEOM_Arc a;
      GSHGEOM_Arc_Form fm;
      double t0, sw, cx, cy, ux, uy, ratio;
      memset(&a, 0, sizeof a);
      a.unit = units[i];
      t0 = t0s[k % 5] * quarters[i];
      sw = mult[(k / 5) % 13] * quarters[i];
      if (k >= 65) { /* random bases */
        double d0 = grid_unit(&rng);
        double d1 = grid_unit(&rng);
        t0 = (d0 - 0.5) * 3.0 * turns[i];
        sw = (d1 - 0.5) * 2.0 * turns[i];
      }
      cx = grid_unit(&rng) * 20.0 - 10.0;
      cy = grid_unit(&rng) * 20.0 - 10.0;
      a.center = grid_P(cx, cy);
      ux = grid_unit(&rng) * 8.0 + 1.0;
      uy = grid_unit(&rng) * 4.0;
      a.u = grid_P(ux, uy);
      ratio = 0.3 + grid_unit(&rng);
      a.v = grid_P(-a.u.y * ratio, a.u.x * ratio);
      if (k % 7 == 0) { /* a circle */
        a.v = grid_P(-a.u.y, a.u.x);
        a.flags = GSHGEOM_ARC_CIRCULAR;
      }
      a.t0 = t0;
      a.sweep = sw;
      fm.center = a.center;
      fm.u = a.u;
      fm.v = a.v;
      fm.unit = a.unit;
      fm.t0 = a.t0;
      fm.sweep = a.sweep;
      fm.flags = a.flags;
      gshgeom_arc_point(&fm, a.t0, &a.start);
      gshgeom_arc_point(&fm, a.t0 + a.sweep, &a.end);
      GRID_PUSH(a);
      if (k < 40) { /* each field disturbed in turn, from this base */
        int bi, field;
        GSHGEOM_Arc d;
        for (bi = 0; bi < 8; bi++) {
          for (field = 0; field < 16; field++) {
            double * targets[16];
            d = a;
            targets[0] = &d.start.x;
            targets[1] = &d.start.y;
            targets[2] = &d.end.x;
            targets[3] = &d.end.y;
            targets[4] = &d.center.x;
            targets[5] = &d.center.y;
            targets[6] = &d.u.x;
            targets[7] = &d.u.y;
            targets[8] = &d.v.x;
            targets[9] = &d.v.y;
            targets[10] = &d.t0;
            targets[11] = &d.sweep;
            targets[12] = &d.start.z;
            targets[13] = &d.end.z;
            targets[14] = &d.center.z;
            targets[15] = &d.u.z;
            *targets[field] = bad[bi];
            GRID_PUSH(d);
          }
        }
        d = a;
        d.unit = GSHGEOM_ANGLE_UNIT_COUNT;
        GRID_PUSH(d);
        d = a;
        d.flags = 2u;
        GRID_PUSH(d);
        d = a;
        d.flags = a.flags ^ 1u;
        GRID_PUSH(d);
        d = a;
        d.sweep = a.sweep >= 0 ? turns[i] : -turns[i];
        GRID_PUSH(d);
        d = a;
        d.sweep = (a.sweep >= 0 ? 1.0 : -1.0) * grid_nextafter(turns[i], 2.0 * turns[i]);
        GRID_PUSH(d);
        d = a;
        d.v = d.u;
        GRID_PUSH(d);
        d = a;
        d.v = grid_P(-a.u.y * 0x1p-1000, a.u.x * 0x1p-1000);
        GRID_PUSH(d);
        d = a;
        d.end.x = grid_nextafter(a.end.x, 1e300);
        GRID_PUSH(d);
        d = a;
        d.end.x = a.end.x + 1e-9 * (grid_fabs(a.u.x) + 1.0);
        GRID_PUSH(d);
        d = a;
        d.start.y = a.start.y + 1e-9 * (grid_fabs(a.u.x) + 1.0);
        GRID_PUSH(d);
      }
    }
  }
#undef GRID_PUSH
  *count_out = n;
  return out;
}

static inline void grid_arc_records(Grid_Fold * f) {
  size_t count = 0;
  GSHGEOM_Arc * all = grid_arc_family(&count);
  GSHGEOM_Map ident;
  GSHGEOM_Map maps[6];
  GSHGEOM_Map skew;
  const GSHGEOM_Map * under[8];
  size_t ai;
  gshgeom_map_identity(2, &ident);
  gshgeom_map_scale(2, 2.0, 2.0, 1.0, &maps[0]);
  gshgeom_map_scale(2, 2.0, 0.5, 1.0, &maps[1]);
  gshgeom_map_scale(2, 0.0, 1.0, 1.0, &maps[2]);
  gshgeom_map_rotate(GSHGEOM_ANGLE_DEGREES, 30.0, &maps[3]);
  gshgeom_map_scale(2, -1.0, 1.0, 1.0, &maps[4]);
  gshgeom_map_translate(2, 1e30, -3.0, 0.0, &maps[5]);
  skew = maps[3];
  skew.m[0][1] += 0.5;
  under[0] = &ident;
  under[1] = &maps[0];
  under[2] = &maps[1];
  under[3] = &maps[2];
  under[4] = &maps[3];
  under[5] = &maps[4];
  under[6] = &maps[5];
  under[7] = &skew;
  for (ai = 0; ai < count; ai++) {
    const GSHGEOM_Arc * a = &all[ai];
    int mode, k, q, t, mi;
    GSHGEOM_Curve pieces[GSHGEOM_ARC_MAX_PIECES];
    int n = -5;
    GSHGEOM_Result pr;
    GSHGEOM_Point c;
    double radius = -1.0;
    int ccw = -7;
    GSHGEOM_Arc_Form fm;
    double ts[5];
    /* append at a path whose current point is the arc's start, at another, at none */
    for (mode = 0; mode < 3; mode++) {
      GSHGEOM_Path * p = grid_make();
      GSHGEOM_Result r;
      if (mode == 0) {
        gshgeom_path_move_to(p, a->start);
      } else if (mode == 1) {
        gshgeom_path_move_to(p, grid_P(a->start.x + 1.0, a->start.y));
      }
      r = gshgeom_path_append_arc(p, a);
      grid_fold_path(f, r, p);
      if (r == GSHGEOM_OK) {
        GSHGEOM_Arc back;
        GSHGEOM_Result rr;
        GSHGEOM_Path * low = NULL;
        GSHGEOM_Result lr;
        memset(&back, 0x11, sizeof back);
        rr = gshgeom_path_arc(p, 1, &back);
        grid_word(f, (uint64_t)rr);
        grid_arc(f, &back);
        grid_word(f, (uint64_t)gshgeom_path_arc(p, 0, &back));
        grid_word(f, (uint64_t)gshgeom_path_arc(p, 2, &back));
        lr = gshgeom_path_lower_arcs(p, NULL, NULL, &low);
        grid_word(f, (uint64_t)lr);
        if (lr == GSHGEOM_OK && low) {
          grid_word(f, gshgeom_path_hash(low));
          grid_word(f, gshgeom_path_element_count(low));
        }
        gshgeom_path_destroy(low);
      }
      gshgeom_path_destroy(p);
    }
    /* the pieces, the circle, the point */
    memset(pieces, 0, sizeof pieces);
    pr = gshgeom_arc_pieces(a, pieces, &n);
    grid_word(f, (uint64_t)pr);
    grid_word(f, (uint64_t)(int64_t)n);
    if (pr == GSHGEOM_OK) {
      for (k = 0; k < n; k++) {
        grid_word(f, (uint64_t)pieces[k].kind);
        grid_word(f, (uint64_t)pieces[k].dim);
        for (q = 0; q < 4; q++) {
          grid_pt(f, &pieces[k].p[q]);
          grid_dbl(f, pieces[k].w[q]);
        }
      }
    }
    c = grid_P(-1.0, -1.0);
    grid_word(f, (uint64_t)gshgeom_arc_circle(a, &c, &radius, &ccw));
    grid_pt(f, &c);
    grid_dbl(f, radius);
    grid_word(f, (uint64_t)(int64_t)ccw);
    fm.center = a->center;
    fm.u = a->u;
    fm.v = a->v;
    fm.unit = a->unit;
    fm.t0 = a->t0;
    fm.sweep = a->sweep;
    fm.flags = a->flags;
    ts[0] = a->t0;
    ts[1] = a->t0 + a->sweep;
    ts[2] = 0.0;
    ts[3] = 1e300;
    ts[4] = (double)NAN;
    for (t = 0; t < 5; t++) {
      GSHGEOM_Point pt = grid_P(-9.0, -9.0);
      GSHGEOM_Result tr = gshgeom_arc_point(&fm, ts[t], &pt);
      grid_word(f, (uint64_t)tr);
      grid_pt(f, &pt);
    }
    /* under maps */
    for (mi = 0; mi < 8; mi++) {
      GSHGEOM_Arc out;
      GSHGEOM_Result mr;
      memset(&out, 0x3c, sizeof out);
      mr = gshgeom_arc_map(under[mi], a, &out);
      grid_word(f, (uint64_t)mr);
      grid_arc(f, &out);
    }
  }
  free(all);
}

/** The entry points on what they must refuse whatever else is right: no path, no record, no output, an index past the end, a path of three dimensions. */
static inline void grid_arc_refusals(Grid_Fold * f) {
  GSHGEOM_Arc a;
  GSHGEOM_Arc_Form fm;
  GSHGEOM_Curve pieces[GSHGEOM_ARC_MAX_PIECES];
  int n = 0;
  GSHGEOM_Point c;
  double radius;
  int ccw;
  GSHGEOM_Arc back;
  GSHGEOM_Map ident;
  GSHGEOM_Point endpt = grid_P(9.0, 9.0);
  GSHGEOM_Arc_Form bad;
  static const size_t indices[6] = {0, 1, 2, 3, 100, (size_t)-1};
  int dim, withcur, ii;
  memset(&a, 0, sizeof a);
  a.center = grid_P(2.0, 3.0);
  a.u = grid_P(4.0, 0.0);
  a.v = grid_P(0.0, 4.0);
  a.unit = GSHGEOM_ANGLE_DEGREES;
  a.t0 = 10.0;
  a.sweep = 100.0;
  a.flags = GSHGEOM_ARC_CIRCULAR;
  fm.center = a.center;
  fm.u = a.u;
  fm.v = a.v;
  fm.unit = a.unit;
  fm.t0 = a.t0;
  fm.sweep = a.sweep;
  fm.flags = a.flags;
  gshgeom_arc_point(&fm, a.t0, &a.start);
  gshgeom_arc_point(&fm, a.t0 + a.sweep, &a.end);
  memset(&back, 0, sizeof back);
  memset(pieces, 0, sizeof pieces);
  gshgeom_map_identity(2, &ident);
  for (dim = 2; dim <= 3; dim++) {
    for (withcur = 0; withcur < 2; withcur++) {
      GSHGEOM_Path * p = NULL;
      GSHGEOM_Path * low = NULL;
      GSHGEOM_Arc_Kind kind = GSHGEOM_ARC_KIND_NOTHING;
      gshgeom_path_create(NULL, NULL, dim, &p);
      if (withcur) {
        gshgeom_path_move_to(p, a.start);
      }
      grid_word(f, (uint64_t)gshgeom_path_append_arc(p, &a));
      grid_word(f, (uint64_t)gshgeom_path_append_arc(p, NULL));
      grid_word(f, (uint64_t)gshgeom_path_append_arc(NULL, &a));
      for (ii = 0; ii < 6; ii++) {
        memset(&back, 0x11, sizeof back);
        grid_word(f, (uint64_t)gshgeom_path_arc(p, indices[ii], &back));
        grid_arc(f, &back);
        grid_word(f, (uint64_t)gshgeom_path_arc(p, indices[ii], NULL));
      }
      grid_word(f, (uint64_t)gshgeom_path_arc(NULL, 0, &back));
      grid_word(f, (uint64_t)gshgeom_path_arc_endpoint_to(p, grid_P(9.0, 1.0), 5.0, 5.0, 0.0, GSHGEOM_ANGLE_DEGREES, 0, 1, &kind));
      grid_word(f, (uint64_t)gshgeom_path_arc_endpoint_to(p, grid_P(9.0, 1.0), 5.0, 5.0, 0.0, GSHGEOM_ANGLE_DEGREES, 0, 1, NULL));
      grid_word(f, (uint64_t)gshgeom_path_arc_endpoint_to(NULL, grid_P(9.0, 1.0), 5.0, 5.0, 0.0, GSHGEOM_ANGLE_DEGREES, 0, 1, &kind));
      grid_word(f, (uint64_t)gshgeom_path_arc_conjugate_to(p, &fm, &endpt));
      grid_word(f, (uint64_t)gshgeom_path_arc_conjugate_to(p, &fm, NULL));
      grid_word(f, (uint64_t)gshgeom_path_arc_conjugate_to(p, NULL, NULL));
      grid_word(f, (uint64_t)gshgeom_path_arc_conjugate_to(NULL, &fm, NULL));
      grid_word(f, (uint64_t)gshgeom_path_arc_center_to(NULL, grid_P(0, 0), 1.0, 1.0, 0.0, GSHGEOM_ANGLE_DEGREES, 0.0, 90.0, NULL));
      grid_word(f, (uint64_t)gshgeom_path_arc_tangent_to(NULL, grid_P(1, 1), grid_P(2, 0), 1.0));
      grid_word(f, (uint64_t)gshgeom_path_arc_bulge_to(NULL, grid_P(1, 1), 1.0));
      grid_fold_path(f, GSHGEOM_OK, p);
      grid_word(f, (uint64_t)gshgeom_path_lower_arcs(p, NULL, NULL, &low));
      gshgeom_path_destroy(low);
      low = NULL;
      grid_word(f, (uint64_t)gshgeom_path_lower_arcs(NULL, NULL, NULL, &low));
      grid_word(f, (uint64_t)gshgeom_path_lower_arcs(p, NULL, NULL, NULL));
      gshgeom_path_destroy(p);
    }
  }
  grid_word(f, (uint64_t)gshgeom_arc_pieces(NULL, pieces, &n));
  grid_word(f, (uint64_t)gshgeom_arc_pieces(&a, NULL, &n));
  grid_word(f, (uint64_t)gshgeom_arc_pieces(&a, pieces, NULL));
  grid_word(f, (uint64_t)gshgeom_arc_circle(NULL, &c, &radius, &ccw));
  grid_word(f, (uint64_t)gshgeom_arc_circle(&a, NULL, &radius, &ccw));
  grid_word(f, (uint64_t)gshgeom_arc_circle(&a, &c, NULL, &ccw));
  grid_word(f, (uint64_t)gshgeom_arc_circle(&a, &c, &radius, NULL));
  grid_word(f, (uint64_t)gshgeom_arc_map(NULL, &a, &back));
  grid_word(f, (uint64_t)gshgeom_arc_map(&ident, NULL, &back));
  grid_word(f, (uint64_t)gshgeom_arc_map(&ident, &a, NULL));
  grid_word(f, (uint64_t)gshgeom_arc_point(NULL, 0.0, &c));
  grid_word(f, (uint64_t)gshgeom_arc_point(&fm, 0.0, NULL));
  bad = fm;
  memset(&bad.unit, 0xff, sizeof bad.unit); /* a unit that is no unit, made without naming a value the enum does not have */
  grid_word(f, (uint64_t)gshgeom_arc_point(&bad, 0.0, &c));
}

/** The exact sums the arc routines decide with, converted to a double: ties, sticky bits, cancellation, the ends of the exponent range. */
static inline double grid_sum_number(Grid_Rng * rng) {
  uint64_t r = grid_next(rng);
  int kind = (int)(r % 9);
  int e = (int)((r >> 8) % 120) - 60;
  uint64_t bitsrc = grid_next(rng);
  double m = 1.0 + (double)((bitsrc >> 12) & ((1ull << 52) - 1)) * 0x1p-52;
  switch (kind) {
    case 0: return gmath_scalbn(m, e);
    case 1: return -gmath_scalbn(m, e);
    case 2: return gmath_scalbn(1.0, e);
    case 3: return gmath_scalbn(1.0 + 0x1p-52, e);
    case 4: return gmath_scalbn(1.0 + 0x1p-26, e);
    case 5: return gmath_scalbn(m, e - 900);
    case 6: return gmath_scalbn(m, e + 900);
    case 7: return 0.0;
    default: return gmath_scalbn(1.0 - 0x1p-53, e);
  }
}

static inline void grid_arc_sums(Grid_Fold * f) {
  Grid_Rng rng;
  int i, k;
  rng.s = 0x9e3779b97f4a7c15ull ^ (grid_seed_salt * 0x9e3779b97f4a7c15ull);
  for (i = 0; i < 40000; i++) {
    GSHGEOM_ExactSum sum;
    int terms, t;
    gshgeom_exact_sum_init(&sum);
    terms = 1 + (int)(grid_next(&rng) % 4);
    for (t = 0; t < terms; t++) {
      double fac[3];
      int n, neg;
      fac[0] = grid_sum_number(&rng);
      fac[1] = grid_sum_number(&rng);
      fac[2] = grid_sum_number(&rng);
      n = 1 + (int)(grid_next(&rng) % 3);
      neg = (int)(grid_next(&rng) & 1);
      gshgeom_exact_sum_add(&sum, neg, n, fac);
      if ((grid_next(&rng) & 7) == 0) { /* the term again, negated, with a small one left over: cancellation */
        double tiny[1];
        int e2 = (int)(grid_next(&rng) % 200) - 100;
        int neg2 = (int)(grid_next(&rng) & 1);
        gshgeom_exact_sum_add(&sum, 1, n, fac);
        tiny[0] = gmath_scalbn(1.0, e2);
        gshgeom_exact_sum_add(&sum, neg2, 1, tiny);
      }
    }
    grid_dbl(f, gshgeom_exact_sum_to_double(&sum));
    grid_word(f, (uint64_t)(int64_t)gshgeom_exact_sum_sign(&sum));
  }
  /* ties: 2^53 + 1, 2^53 + 3, 2^54 + 2 and the sticky bit far below */
  for (k = 0; k < 140; k++) {
    GSHGEOM_ExactSum sum;
    double big[1];
    double odd[1];
    double sticky[1];
    gshgeom_exact_sum_init(&sum);
    big[0] = gmath_scalbn(1.0, 53 + k % 3);
    odd[0] = (double)(1 + 2 * (k % 5));
    sticky[0] = gmath_scalbn(1.0, -(k + 1));
    gshgeom_exact_sum_add(&sum, k & 1, 1, big);
    gshgeom_exact_sum_add(&sum, k & 1, 1, odd);
    if (k % 4 >= 2) {
      gshgeom_exact_sum_add(&sum, k & 1, 1, sticky);
    }
    grid_dbl(f, gshgeom_exact_sum_to_double(&sum));
  }
}

/** The arcs' grid: the hash of every call and the number of words folded. */
static inline uint64_t grid_arc_hash(uint64_t * words) {
  Grid_Fold f;
  grid_fold_init(&f);
  grid_arc_endpoint(&f);
  grid_arc_center(&f);
  grid_arc_tangent_bulge(&f);
  grid_arc_records(&f);
  grid_arc_refusals(&f);
  grid_arc_sums(&f);
  *words = f.n;
  return f.h;
}

/* ---- lattices ------------------------------------------------------------------------------------------------------------- */

typedef struct Grid_Rows {
  Grid_Fold * f;
  uint64_t rows;
  uint64_t stop_after;
} Grid_Rows;

static inline int grid_row_fn(void * ctx, int64_t j, int64_t lo, int64_t hi) {
  Grid_Rows * r = (Grid_Rows *)ctx;
  grid_word(r->f, (uint64_t)j);
  grid_word(r->f, (uint64_t)lo);
  grid_word(r->f, (uint64_t)hi);
  r->rows++;
  return r->stop_after && r->rows >= r->stop_after ? 7 : 0;
}

static inline GSHGEOM_Map grid_lattice_map(int k) {
  GSHGEOM_Map m;
  gshgeom_map_identity(2, &m);
  switch (k) {
    case 0: break;
    case 1: gshgeom_map_scale(2, 2.0, 2.0, 1.0, &m); break;
    case 2: gshgeom_map_rotate(GSHGEOM_ANGLE_DEGREES, 30.0, &m); break;
    case 3: gshgeom_map_rotate(GSHGEOM_ANGLE_DEGREES, 90.0, &m); break;
    case 4: gshgeom_map_scale(2, 0.0, 1.0, 1.0, &m); break;
    case 5: gshgeom_map_scale(2, 1e-30, 1e-30, 1.0, &m); break;
    case 6: gshgeom_map_scale(2, 1e30, 1e30, 1.0, &m); break;
    case 7: gshgeom_map_scale(2, -1.0, 1.0, 1.0, &m); break;
    case 8: gshgeom_map_translate(2, 1e15, -3.0, 0.0, &m); break;
    case 9:
      gshgeom_map_rotate(GSHGEOM_ANGLE_DEGREES, 30.0, &m);
      m.m[0][1] += 0.5;
      m.kind = 0xffffffffu; /* an invalid kind tag */
      break;
    case 10: gshgeom_map_scale(2, 3.0, 0.5, 1.0, &m); break;
    case 11: gshgeom_map_identity(3, &m); break;
    case 12: gshgeom_map_rotate(GSHGEOM_ANGLE_DEGREES, 45.0, &m); break;
    default: gshgeom_map_scale(2, 0.25, 4.0, 1.0, &m); break;
  }
  return m;
}

#define GRID_IDX ((int64_t)1 << 53)
#define GRID_LATTICE_MAPS 14

static inline void grid_lattice_run(Grid_Fold * f) {
  static const double xy[3] = {0.0, -1.5, 1e15};
  static const double wh[7] = {1.0, 2.5, 0.0, -1.0, (double)NAN, 1e300, 1e-300};
  static const double steps[10] = {0.0, 1.0, 3.5, (double)NAN, (double)INFINITY, -1.0, 1e-300, 1e300, 0.1, 2.5};
  static const uint32_t reps[3] = {0, 1, 2};
  static const int64_t rects[10][4] = {{0, 0, 3, 3}, {-5, -5, 5, 5}, {0, 0, 0, 0}, {3, 3, 0, 0}, {-100, -100, 100, 100},
      {0, 0, GRID_IDX, 1}, {-GRID_IDX, -GRID_IDX, GRID_IDX, GRID_IDX}, {INT64_MIN, INT64_MIN, INT64_MAX, INT64_MAX},
      {1000000, 1000000, 1000010, 1000010}, {0, 0, 1, 1}};
  static const uint64_t caps[5] = {0, 1, 10, 10000, UINT64_MAX};
  static const int64_t iis[7] = {0, 1, -3, GRID_IDX - 1, GRID_IDX, -GRID_IDX, ((int64_t)1 << 62)};
  static const int64_t jjs[4] = {0, 2, -GRID_IDX + 1, GRID_IDX};
  Grid_Rng rng;
  int it;
  rng.s = 0x2545f4914f6cdd1dull ^ (grid_seed_salt * 0x9e3779b97f4a7c15ull);
  for (it = 0; it < 2200; it++) {
    GSHGEOM_Lattice l;
    int valid, mk;
    uint64_t d;
    valid = (grid_next(&rng) % 4) != 0;
    l.x = xy[grid_next(&rng) % 3];
    l.y = xy[grid_next(&rng) % 3];
    d = grid_next(&rng);
    {
      uint64_t d2 = grid_next(&rng);
      l.w = valid ? (d % 2 ? 1.0 : 2.5) : wh[d2 % 7];
    }
    d = grid_next(&rng);
    {
      uint64_t d2 = grid_next(&rng);
      l.h = valid ? (d % 3 ? 1.0 : 0.75) : wh[d2 % 7];
    }
    {
      /* one draw to a statement: the order of the operands of a product is unspecified */
      double sx_base = steps[1 + grid_next(&rng) % 3];
      double sx_mul = grid_next(&rng) % 2 ? 1.0 : 1.25;
      double sx_any = steps[grid_next(&rng) % 10];
      double sy_base, sy_mul, sy_any;
      l.x_step = valid ? sx_base * sx_mul : sx_any;
      sy_base = steps[1 + grid_next(&rng) % 3];
      sy_mul = grid_next(&rng) % 2 ? 1.0 : 0.5;
      sy_any = steps[grid_next(&rng) % 10];
      l.y_step = valid ? sy_base * sy_mul : sy_any;
    }
    {
      uint64_t a1 = grid_next(&rng);
      uint64_t a2 = grid_next(&rng);
      uint64_t b1 = grid_next(&rng);
      uint64_t b2 = grid_next(&rng);
      l.repeat_x = valid ? (uint32_t)(a1 % 2) : reps[a2 % 3];
      l.repeat_y = valid ? (uint32_t)(b1 % 2) : reps[b2 % 3];
    }
    if (valid && l.repeat_x == 0) {
      uint64_t z = grid_next(&rng);
      l.x_step = z % 2 ? 0.0 : l.x_step;
    }
    for (mk = 0; mk < GRID_LATTICE_MAPS; mk++) {
      GSHGEOM_Map m;
      double out4[4];
      double out8[8];
      double out2[2];
      double a[2], b[2];
      int64_t ia[2], ib[2];
      int integer = -7;
      int ii, jj, r, ci, w;
      if ((it + mk) % 3 == 2) {
        continue;
      }
      m = grid_lattice_map(mk);
      for (w = 0; w < 8; w++) {
        out8[w] = -7.0;
      }
      grid_word(f, (uint64_t)gshgeom_lattice_quad(&m, &l, out8));
      for (w = 0; w < 8; w++) {
        grid_dbl(f, out8[w]);
      }
      a[0] = a[1] = b[0] = b[1] = -7;
      grid_word(f, (uint64_t)gshgeom_lattice_steps(&m, &l, a, b));
      grid_dbl(f, a[0]);
      grid_dbl(f, a[1]);
      grid_dbl(f, b[0]);
      grid_dbl(f, b[1]);
      ia[0] = ia[1] = ib[0] = ib[1] = -7;
      grid_word(f, (uint64_t)gshgeom_lattice_integer_steps(&m, &l, ia, ib, &integer));
      grid_word(f, (uint64_t)ia[0]);
      grid_word(f, (uint64_t)ia[1]);
      grid_word(f, (uint64_t)ib[0]);
      grid_word(f, (uint64_t)ib[1]);
      grid_word(f, (uint64_t)(int64_t)integer);
      for (ii = 0; ii < 7; ii++) {
        for (jj = 0; jj < 4; jj++) {
          GSHGEOM_Map can;
          int64_t sx = -7, sy = -7;
          int row;
          if (mk % 4 != 0 && (iis[ii] == GRID_IDX || jjs[jj] == GRID_IDX)) {
            continue;
          }
          out4[0] = out4[1] = out4[2] = out4[3] = -7.0;
          out2[0] = out2[1] = -7.0;
          grid_word(f, (uint64_t)gshgeom_lattice_cell(&l, iis[ii], jjs[jj], out4));
          for (w = 0; w < 4; w++) {
            grid_dbl(f, out4[w]);
          }
          grid_word(f, (uint64_t)gshgeom_lattice_translation(&m, &l, iis[ii], jjs[jj], out2));
          grid_dbl(f, out2[0]);
          grid_dbl(f, out2[1]);
          memset(&can, 0x42, sizeof can);
          grid_word(f, (uint64_t)gshgeom_lattice_canonical(&m, &l, iis[ii], jjs[jj], &can, &sx, &sy));
          grid_word(f, (uint64_t)sx);
          grid_word(f, (uint64_t)sy);
          grid_word(f, can.kind);
          for (row = 0; row < 3; row++) {
            grid_dbl(f, can.m[row][0]);
            grid_dbl(f, can.m[row][1]);
          }
          grid_dbl(f, can.t[0]);
          grid_dbl(f, can.t[1]);
        }
      }
      /* the enumeration, over the rectangles and the caps; stopped after a few rows */
      for (r = 0; r < 10; r++) {
        if (((it + mk + r) % 4) != 0 && r != 0 && r != 1) {
          continue;
        }
        for (ci = 0; ci < 5; ci++) {
          Grid_Rows rows;
          uint64_t count = 12345;
          if (((it + mk + r + (int)(caps[ci] % 7)) % 3) != 0) {
            continue;
          }
          rows.f = f;
          rows.rows = 0;
          rows.stop_after = (it % 5 == 0) ? 3 : 4000;
          grid_word(f, (uint64_t)gshgeom_lattice_each(&m, &l, rects[r], caps[ci], grid_row_fn, &rows, &count));
          grid_word(f, count);
          grid_word(f, rows.rows);
        }
      }
    }
  }
  /* each with no callback / no count / no rect */
  {
    GSHGEOM_Lattice l;
    GSHGEOM_Map m = grid_lattice_map(0);
    static const int64_t r0[4] = {0, 0, 3, 3};
    uint64_t count = 0;
    Grid_Rows rows;
    l.x = 0.0;
    l.y = 0.0;
    l.w = 1.0;
    l.h = 1.0;
    l.x_step = 2.0;
    l.y_step = 2.0;
    l.repeat_x = 1;
    l.repeat_y = 1;
    rows.f = f;
    rows.rows = 0;
    rows.stop_after = 0;
    grid_word(f, (uint64_t)gshgeom_lattice_each(&m, &l, r0, 100, NULL, &rows, &count));
    grid_word(f, (uint64_t)gshgeom_lattice_each(&m, &l, r0, 100, grid_row_fn, &rows, NULL));
    grid_word(f, (uint64_t)gshgeom_lattice_each(&m, &l, NULL, 100, grid_row_fn, &rows, &count));
    grid_word(f, (uint64_t)gshgeom_lattice_each(NULL, &l, r0, 100, grid_row_fn, &rows, &count));
    grid_word(f, (uint64_t)gshgeom_lattice_each(&m, NULL, r0, 100, grid_row_fn, &rows, &count));
  }
}

static inline void grid_lattice_period(Grid_Fold * f) {
  static const int64_t v[16] = {0, 1, -1, 2, 3, -4, 7, 12, 100, -100, ((int64_t)1 << 31), ((int64_t)1 << 31) - 1,
      -((int64_t)1 << 31), ((int64_t)1 << 32), INT64_MAX, INT64_MIN};
  static const int64_t xs[5] = {0, 5, -9, 1000003, -((int64_t)1 << 40)};
  static const int64_t ys[4] = {0, 3, -17, ((int64_t)1 << 40)};
  int a0, a1, b0, b1;
  GSHGEOM_Lattice_Period pr;
  int64_t a[2];
  for (a0 = 0; a0 < 16; a0++) {
    for (a1 = 0; a1 < 16; a1++) {
      for (b0 = 0; b0 < 16; b0++) {
        for (b1 = 0; b1 < 16; b1++) {
          if (((a0 * 7 + a1 * 5 + b0 * 3 + b1) % 3) == 0 || (a0 < 8 && a1 < 8 && b0 < 8 && b1 < 8)) {
            int64_t aa[2];
            int64_t bb[2];
            GSHGEOM_Result r;
            aa[0] = v[a0];
            aa[1] = v[a1];
            bb[0] = v[b0];
            bb[1] = v[b1];
            memset(&pr, 0x77, sizeof pr);
            r = gshgeom_lattice_period(aa, bb, &pr);
            grid_word(f, (uint64_t)r);
            grid_word(f, (uint64_t)pr.p);
            grid_word(f, (uint64_t)pr.r);
            grid_word(f, (uint64_t)pr.s);
            if (r == GSHGEOM_OK) {
              int xi, yi;
              for (xi = 0; xi < 5; xi++) {
                for (yi = 0; yi < 4; yi++) {
                  int64_t out[2];
                  out[0] = out[1] = -7;
                  gshgeom_lattice_reduce(&pr, xs[xi], ys[yi], out);
                  grid_word(f, (uint64_t)out[0]);
                  grid_word(f, (uint64_t)out[1]);
                }
              }
            }
          }
        }
      }
    }
  }
  a[0] = 1;
  a[1] = 0;
  grid_word(f, (uint64_t)gshgeom_lattice_period(NULL, a, &pr));
  grid_word(f, (uint64_t)gshgeom_lattice_period(a, NULL, &pr));
  grid_word(f, (uint64_t)gshgeom_lattice_period(a, a, NULL));
}

/** The lattices' grid: the hash of every call and the number of words folded. */
static inline uint64_t grid_lattice_hash(uint64_t * words) {
  Grid_Fold f;
  grid_fold_init(&f);
  grid_lattice_run(&f);
  grid_lattice_period(&f);
  *words = f.n;
  return f.h;
}

/* ---- clipping -------------------------------------------------------------------------------------------------------------- */

static inline GSHGEOM_Device_Point grid_dp(int64_t x, int64_t y) {
  GSHGEOM_Device_Point p;
  p.x = (int32_t)x;
  p.y = (int32_t)y;
  return p;
}

static inline int64_t grid_coord(Grid_Rng * rng, int64_t range) {
  return (int64_t)(grid_next(rng) % (uint64_t)(2 * range + 1)) - range;
}

/** The clips' grid: the hash of every call and the number of words folded. */
static inline uint64_t grid_clip_hash(uint64_t * words) {
  static const int64_t ext[10] = {4, 8, 16, 32, 64, 640, 6400, 1 << 20, 1 << 30, (int64_t)INT32_MAX - 1};
  static const size_t caps[4] = {0, 3, 8, 40};
  Grid_Fold fold;
  Grid_Fold * f = &fold;
  Grid_Rng rng;
  int it;
  GSHGEOM_Flat_Path * res = NULL;
  GSHGEOM_Device_Point q[4];
  grid_fold_init(f);
  rng.s = 0x853c49e6748fea9bull ^ (grid_seed_salt * 0x9e3779b97f4a7c15ull);
  for (it = 0; it < 30000; it++) {
    int64_t range = ext[grid_next(&rng) % 10];
    GSHGEOM_Device_Point pts[32];
    GSHGEOM_Flat_Subpath subs[3];
    size_t npts = 0;
    int ns = 1 + (int)(grid_next(&rng) % 3);
    int k;
    size_t qq;
    GSHGEOM_Flat_Path * fp = NULL;
    GSHGEOM_Device_Point quads[8][4];
    int64_t a, b, w, h;
    int qi, c;
    for (k = 0; k < ns; k++) {
      GSHGEOM_Flat_Subpath sp;
      memset(&sp, 0, sizeof sp);
      sp.first = npts;
      sp.count = 2 + grid_next(&rng) % 8;
      sp.closed = (int)(grid_next(&rng) % 2);
      for (qq = 0; qq < sp.count; qq++) {
        uint64_t rep = grid_next(&rng);
        if (rep % 6 == 0 && npts > 0) {
          pts[npts] = pts[npts - 1]; /* a repeat */
        } else {
          int64_t px = grid_coord(&rng, range); /* one draw to a statement */
          int64_t py = grid_coord(&rng, range);
          pts[npts] = grid_dp(px, py);
        }
        npts++;
      }
      subs[k] = sp;
    }
    if (gshgeom_flat_path_from_points(pts, npts, subs, (size_t)ns, NULL, NULL, &fp) != GSHGEOM_OK) {
      continue;
    }
    /* quads: a rectangle, a rotated square, a random quadrilateral (maybe not convex), degenerate ones */
    a = grid_coord(&rng, range / 2 + 64);
    b = grid_coord(&rng, range / 2 + 64);
    w = 64 + (int64_t)(grid_next(&rng) % (uint64_t)(range + 64));
    h = 64 + (int64_t)(grid_next(&rng) % (uint64_t)(range + 64));
    quads[0][0] = grid_dp(a, b);
    quads[0][1] = grid_dp(a + w, b);
    quads[0][2] = grid_dp(a + w, b + h);
    quads[0][3] = grid_dp(a, b + h);
    quads[1][0] = grid_dp(a, b + h);
    quads[1][1] = grid_dp(a + w, b + h);
    quads[1][2] = grid_dp(a + w, b);
    quads[1][3] = grid_dp(a, b); /* reversed */
    quads[2][0] = grid_dp(a + w / 2, b);
    quads[2][1] = grid_dp(a + w, b + h / 2);
    quads[2][2] = grid_dp(a + w / 2, b + h);
    quads[2][3] = grid_dp(a, b + h / 2);
    for (c = 0; c < 4; c++) {
      int64_t qx = grid_coord(&rng, range);
      int64_t qy = grid_coord(&rng, range);
      quads[3][c] = grid_dp(qx, qy);
    }
    quads[4][0] = grid_dp(a, b);
    quads[4][1] = grid_dp(a, b);
    quads[4][2] = grid_dp(a + w, b);
    quads[4][3] = grid_dp(a + w, b + h); /* a repeat */
    quads[5][0] = grid_dp(a, b);
    quads[5][1] = grid_dp(a + w, b);
    quads[5][2] = grid_dp(a + 2 * w, b);
    quads[5][3] = grid_dp(a + w, b + h); /* collinear corners */
    quads[6][0] = grid_dp(a, b);
    quads[6][1] = grid_dp(a + w, b + h);
    quads[6][2] = grid_dp(a + w, b);
    quads[6][3] = grid_dp(a, b + h); /* bow tie */
    quads[7][0] = grid_dp(a, b);
    quads[7][1] = grid_dp(a + w, b);
    quads[7][2] = grid_dp(a + w, b);
    quads[7][3] = grid_dp(a, b); /* no area */
    for (qi = 0; qi < 8; qi++) {
      size_t cap;
      GSHGEOM_Limits lim;
      GSHGEOM_Result r;
      if (((it + qi) % 2) && qi > 2) {
        continue;
      }
      cap = caps[(it + qi) % 4];
      gshgeom_limits_default(&lim);
      lim.max_points = cap;
      res = (GSHGEOM_Flat_Path *)(uintptr_t)1;
      r = gshgeom_flat_path_clip_convex(fp, quads[qi], NULL, cap ? &lim : NULL, &res);
      grid_word(f, (uint64_t)r);
      if (r == GSHGEOM_OK) {
        grid_word(f, gshgeom_flat_path_hash(res));
        grid_word(f, gshgeom_flat_path_point_count(res));
        grid_word(f, gshgeom_flat_path_subpath_count(res));
        gshgeom_flat_path_destroy(res);
      } else {
        grid_word(f, res == (GSHGEOM_Flat_Path *)(uintptr_t)1);
      }
    }
    gshgeom_flat_path_destroy(fp);
  }
  /* NULL arguments */
  res = NULL;
  q[0] = grid_dp(0, 0);
  q[1] = grid_dp(64, 0);
  q[2] = grid_dp(64, 64);
  q[3] = grid_dp(0, 64);
  grid_word(f, (uint64_t)gshgeom_flat_path_clip_convex(NULL, q, NULL, NULL, &res));
  *words = f->n;
  return f->h;
}


/* ---- walkers, shapes and strokes ----------------------------------------------------------------------------------------- */

/** A coordinate on the grid of an eighth, in [-range, range]: one draw. */
static inline double grid_q8(Grid_Rng * rng, int range) {
  uint64_t d = grid_next(rng);
  return (double)((int64_t)(d % (uint64_t)(range * 16 + 1)) - (int64_t)(range * 8)) / 8.0;
}

/** A positive length on the grid of an eighth, in [1/8, range]: one draw. */
static inline double grid_len8(Grid_Rng * rng, int range) {
  uint64_t d = grid_next(rng);
  return (double)(1 + (int64_t)(d % (uint64_t)(range * 8))) / 8.0;
}

/** A path of a few random elements of every verb, arcs among them (every constructor of arc.h). */
static inline GSHGEOM_Path * grid_walk_path(Grid_Rng * rng) {
  GSHGEOM_Path * p = grid_make();
  int n, i;
  double x = grid_q8(rng, 20);
  double y = grid_q8(rng, 20);
  gshgeom_path_move_to(p, grid_P(x, y));
  n = 3 + (int)(grid_next(rng) % 6);
  for (i = 0; i < n; i++) {
    int kind = (int)(grid_next(rng) % 10);
    double ax = grid_q8(rng, 20);
    double ay = grid_q8(rng, 20);
    double bx = grid_q8(rng, 20);
    double by = grid_q8(rng, 20);
    double cx = grid_q8(rng, 20);
    double cy = grid_q8(rng, 20);
    switch (kind) {
      case 0:
      case 1:
        gshgeom_path_line_to(p, grid_P(ax, ay));
        break;
      case 2:
        gshgeom_path_quadratic_to(p, grid_P(ax, ay), grid_P(bx, by));
        break;
      case 3:
        gshgeom_path_cubic_to(p, grid_P(ax, ay), grid_P(bx, by), grid_P(cx, cy));
        break;
      case 4:
      case 5: {
        double rx = grid_len8(rng, 15);
        double ry = grid_len8(rng, 15);
        double phi = grid_q8(rng, 45);
        int large = (int)(grid_next(rng) % 2);
        int sweep = (int)(grid_next(rng) % 2);
        GSHGEOM_Arc_Kind ak;
        gshgeom_path_arc_endpoint_to(p, grid_P(ax, ay), rx, ry, phi, GSHGEOM_ANGLE_DEGREES, large, sweep, &ak);
        break;
      }
      case 6: {
        double rx = grid_len8(rng, 15);
        double ry = grid_len8(rng, 15);
        double t0 = grid_q8(rng, 200);
        double sw = grid_q8(rng, 45) * 8.0;
        gshgeom_path_arc_center_to(p, grid_P(ax, ay), rx, ry, 0.0, GSHGEOM_ANGLE_DEGREES, t0, sw, NULL);
        break;
      }
      case 7:
        gshgeom_path_arc_tangent_to(p, grid_P(ax, ay), grid_P(bx, by), grid_len8(rng, 10));
        break;
      case 8:
        gshgeom_path_arc_bulge_to(p, grid_P(ax, ay), grid_q8(rng, 3));
        break;
      default:
        gshgeom_path_close(p);
        gshgeom_path_move_to(p, grid_P(cx, cy));
        break;
    }
  }
  return p;
}

typedef struct Grid_Dump {
  Grid_Fold * f;
} Grid_Dump;

/** The text a dump writes, folded: the stream is a temporary file, the only portable sink of a FILE *. */
static inline void grid_fold_text(Grid_Fold * f, FILE * t) {
  long n;
  long i;
  if (!t) {
    grid_word(f, 0xdeadull);
    return;
  }
  n = ftell(t);
  rewind(t);
  grid_word(f, (uint64_t)n);
  for (i = 0; i < n; i++) {
    int c = fgetc(t);
    f->h ^= (uint64_t)(unsigned char)c;
    f->h *= 1099511628211ull;
  }
  f->n++;
  fclose(t);
}

static inline GSHGEOM_Result grid_vertex_fn(const GSHGEOM_Vertex * v, void * ctx) {
  Grid_Fold * f = (Grid_Fold *)ctx;
  grid_word(f, v->element);
  grid_word(f, v->subpath);
  grid_word(f, v->flags);
  grid_pt(f, &v->point);
  grid_pt(f, &v->in);
  grid_pt(f, &v->out);
  grid_pt(f, &v->direction);
  return GSHGEOM_OK;
}

static inline void grid_fold_bounds(Grid_Fold * f, GSHGEOM_Result r, const GSHGEOM_Bounds * b) {
  grid_word(f, (uint64_t)r);
  grid_word(f, (uint64_t)(int64_t)b->empty);
  grid_pt(f, &b->min);
  grid_pt(f, &b->max);
}

static inline void grid_walk_one(Grid_Fold * f, Grid_Rng * rng, GSHGEOM_Path * p, GSHGEOM_Path * other) {
  static const double dash_a[2] = {3.0, 2.0};
  static const double dash_b[1] = {1.5};
  static const double dash_c[4] = {0.25, 0.5, 4.0, 1.0};
  size_t n = gshgeom_path_element_count(p);
  size_t i;
  GSHGEOM_Bounds b;
  GSHGEOM_Path * copy = NULL;
  GSHGEOM_Map maps[3];
  int mi;
  FILE * t;
  size_t vcount = 0;
  GSHGEOM_Measured * m = NULL;
  int ti;
  gshgeom_map_identity(2, &maps[0]);
  gshgeom_map_scale(2, 3.0, 0.5, 1.0, &maps[1]);
  gshgeom_map_rotate(GSHGEOM_ANGLE_DEGREES, 30.0, &maps[2]);
  grid_word(f, gshgeom_path_hash(p));
  grid_word(f, n);
  for (i = 0; i <= n; i++) {
    GSHGEOM_Verb verb = GSHGEOM_VERB_MOVE;
    GSHGEOM_Element el;
    GSHGEOM_Arc arc;
    GSHGEOM_Result vr = gshgeom_path_verb(p, i, &verb);
    GSHGEOM_Result er;
    unsigned vb = 0;
    grid_word(f, (uint64_t)vr);
    memcpy(&vb, &verb, sizeof vb < sizeof verb ? sizeof vb : sizeof verb);
    grid_word(f, vb);
    memset(&el, 0x22, sizeof el);
    er = gshgeom_path_element(p, i, &el);
    grid_word(f, (uint64_t)er);
    if (er == GSHGEOM_OK) {
      int q;
      unsigned verb_bits;
      memcpy(&verb_bits, &el.verb, sizeof verb_bits < sizeof el.verb ? sizeof verb_bits : sizeof el.verb);
      grid_word(f, verb_bits);
      grid_word(f, (uint64_t)el.curve.kind);
      grid_word(f, (uint64_t)(int64_t)el.curve.dim);
      for (q = 0; q < 4; q++) {
        grid_pt(f, &el.curve.p[q]);
        grid_dbl(f, el.curve.w[q]);
      }
    }
    memset(&arc, 0, sizeof arc);
    grid_word(f, (uint64_t)gshgeom_path_arc(p, i, &arc));
    grid_arc(f, &arc);
  }
  if (gshgeom_path_copy(p, NULL, NULL, &copy) == GSHGEOM_OK) {
    grid_word(f, gshgeom_path_hash(copy));
    grid_word(f, (uint64_t)(int64_t)gshgeom_path_equal(p, copy));
    grid_word(f, (uint64_t)(int64_t)gshgeom_path_equal(copy, other));
    gshgeom_path_destroy(copy);
  }
  grid_word(f, (uint64_t)(int64_t)gshgeom_path_equal(p, other));
  memset(&b, 0, sizeof b);
  grid_fold_bounds(f, gshgeom_path_bounds(p, &b), &b);
  t = tmpfile();
  if (t) {
    grid_word(f, (uint64_t)gshgeom_path_dump(p, t));
  }
  grid_fold_text(f, t);
  grid_word(f, (uint64_t)gshgeom_path_vertex_count(p, &vcount));
  grid_word(f, vcount);
  grid_word(f, (uint64_t)gshgeom_path_vertices(p, grid_vertex_fn, f));
  for (mi = 0; mi < 3; mi++) {
    GSHGEOM_Device_Path * dp = NULL;
    GSHGEOM_Result dr = gshgeom_device_path_from_model(p, &maps[mi], NULL, NULL, &dp);
    grid_word(f, (uint64_t)dr);
    if (dr == GSHGEOM_OK) {
      grid_word(f, gshgeom_device_path_hash(dp));
      grid_word(f, gshgeom_device_path_element_count(dp));
    }
    gshgeom_device_path_destroy(dp);
  }
  {
    GSHGEOM_Path * low = NULL;
    GSHGEOM_Result lr = gshgeom_stroke_curves_lower_cubic(p, 0.25, NULL, NULL, &low);
    grid_word(f, (uint64_t)lr);
    if (lr == GSHGEOM_OK && low) {
      grid_word(f, gshgeom_path_hash(low));
    }
    gshgeom_path_destroy(low);
  }
  for (ti = 0; ti < 2; ti++) {
    GSHGEOM_Result mr = gshgeom_measured_create(NULL, NULL, p, ti ? 0.0625 : 0.5, &m);
    grid_word(f, (uint64_t)mr);
    if (mr == GSHGEOM_OK && m) {
      double len = gshgeom_measured_length(m);
      int q;
      grid_dbl(f, len);
      grid_word(f, gshgeom_measured_segment_count(m));
      grid_word(f, gshgeom_measured_row_count(m));
      grid_word(f, gshgeom_measured_subpath_count(m));
      grid_word(f, gshgeom_measured_hash(m));
      for (q = -1; q <= 5; q++) {
        GSHGEOM_Measured_Frame fr;
        double s = len * (double)q / 4.0;
        GSHGEOM_Result ar;
        memset(&fr, 0, sizeof fr);
        ar = gshgeom_measured_at(m, s, &fr);
        grid_word(f, (uint64_t)ar);
        grid_pt(f, &fr.point);
        grid_pt(f, &fr.tangent);
        grid_pt(f, &fr.normal);
        grid_word(f, fr.element);
        grid_word(f, fr.piece);
        grid_dbl(f, fr.t);
        memset(&fr, 0, sizeof fr);
        grid_word(f, (uint64_t)gshgeom_measured_at_clamped(m, s, &fr));
        grid_pt(f, &fr.point);
      }
      t = tmpfile();
      if (t) {
        grid_word(f, (uint64_t)gshgeom_measured_dump(m, t));
      }
      grid_fold_text(f, t);
    }
    gshgeom_measured_destroy(m);
    m = NULL;
  }
  /* strokes: four of the joins, caps, dash arrays and maps, drawn at random for each path (the whole product is 45 outlines a path and is slow) */
  for (ti = 0; ti < 4; ti++) {
    GSHGEOM_Stroke st;
    GSHGEOM_Stroke_Outline * o = NULL;
    GSHGEOM_Result sr;
    uint64_t dd = grid_next(rng);
    int dsel = (int)(dd % 4);
    int jn = (int)(dd >> 4 & 7) % 5;
    int cp = (int)(dd >> 8 & 3) % 3;
    mi = (int)(dd >> 12 & 3) % 3;
    memset(&st, 0, sizeof st);
    st.width = 1.0 + (double)(dd >> 16 & 7) * 0.5;
    st.width_space = GSHGEOM_STROKE_WIDTH_LOCAL;
    st.join = (uint32_t)jn;
    st.miter_limit = 4.0;
    st.cap = (uint32_t)cp;
    st.degenerate = GSHGEOM_STROKE_DEGENERATE_DRAW;
    st.dash_count = dsel == 1 ? 2 : dsel == 2 ? 1 : dsel == 3 ? 4 : 0;
    st.dashes = dsel == 1 ? dash_a : dsel == 2 ? dash_b : dsel == 3 ? dash_c : NULL;
    st.dash_offset = dsel ? (double)(dd >> 20 & 3) * 0.75 : 0.0;
    st.dash_scope = (dd >> 24 & 1) ? GSHGEOM_STROKE_DASH_PATH : GSHGEOM_STROKE_DASH_PER_SUBPATH;
    st.alignment = GSHGEOM_STROKE_ALIGN_CENTER;
    sr = gshgeom_stroke_expand(p, &st, &maps[mi], NULL, NULL, &o);
    grid_word(f, (uint64_t)sr);
    if (sr == GSHGEOM_OK && o) {
      grid_word(f, gshgeom_stroke_outline_hash(o));
    }
    gshgeom_stroke_outline_destroy(o);
  }
}

/** The grid of the walkers: random paths of every verb, arcs among them, through every reader and the strokes. */
static inline uint64_t grid_walk_hash(uint64_t * words) {
  Grid_Fold f;
  Grid_Rng rng;
  int it;
  GSHGEOM_Path * prev = NULL;
  grid_fold_init(&f);
  rng.s = 0x9b05688c2b3e6c1full ^ (grid_seed_salt * 0x9e3779b97f4a7c15ull);
  prev = grid_walk_path(&rng);
  for (it = 0; it < 250; it++) {
    GSHGEOM_Path * p = grid_walk_path(&rng);
    grid_walk_one(&f, &rng, p, prev);
    gshgeom_path_destroy(prev);
    prev = p;
  }
  gshgeom_path_destroy(prev);
  *words = f.n;
  return f.h;
}

/** The shapes (rect, rounded rect, ellipse) over special and ordinary arguments, and flat paths from points. */
static inline uint64_t grid_shapes_hash(uint64_t * words) {
  static const double v[8] = {0.0, 1.0, 4.0, -3.0, (double)NAN, (double)INFINITY, GRID_BIG, 0x1p-1000};
  static const double rr[6] = {0.0, 1.0, 2.0, 5.0, -1.0, (double)NAN};
  Grid_Fold f;
  int fx, fy, a, b, c, d, i, j;
  grid_fold_init(&f);
  for (fx = 0; fx < 2; fx++) {
    for (fy = 0; fy < 2; fy++) {
      for (a = 0; a < 8; a++) {
        for (b = 0; b < 8; b++) {
          for (c = 0; c < 8; c++) {
            for (d = 0; d < 8; d++) {
              GSHGEOM_Path * p = grid_make();
              grid_word(&f, (uint64_t)gshgeom_path_add_rect(p, v[a], v[b], v[c], v[d], fx, fy));
              grid_word(&f, gshgeom_path_hash(p));
              gshgeom_path_destroy(p);
              p = grid_make();
              grid_word(&f, (uint64_t)gshgeom_path_add_ellipse(p, v[a], v[b], v[c], v[d], fx, fy));
              grid_word(&f, gshgeom_path_hash(p));
              gshgeom_path_destroy(p);
            }
          }
        }
      }
      for (a = 0; a < 6; a++) {
        for (b = 0; b < 6; b++) {
          for (c = 0; c < 6; c++) {
            for (d = 0; d < 6; d++) {
              GSHGEOM_Corner r[4];
              GSHGEOM_Path * p = grid_make();
              double w = (a % 3) == 0 ? 10.0 : (a % 3) == 1 ? 3.0 : 0.0;
              double h = (b % 3) == 0 ? 8.0 : (b % 3) == 1 ? 2.0 : (double)NAN;
              r[0].rx = rr[a];
              r[0].ry = rr[b];
              r[1].rx = rr[c];
              r[1].ry = rr[d];
              r[2].rx = rr[b];
              r[2].ry = rr[c];
              r[3].rx = rr[d];
              r[3].ry = rr[a];
              grid_word(&f, (uint64_t)gshgeom_rounded_rect_check(w, h, r));
              grid_word(&f, (uint64_t)gshgeom_path_add_rounded_rect(p, 1.0, 2.0, w, h, r, fx, fy));
              grid_word(&f, gshgeom_path_hash(p));
              gshgeom_path_destroy(p);
            }
          }
        }
      }
    }
  }
  /* flat paths from points: subpaths that do not tile the points, none, NULL arrays, counts past the limit */
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      GSHGEOM_Device_Point pts[4];
      GSHGEOM_Flat_Subpath subs[3];
      GSHGEOM_Flat_Path * fp = NULL;
      GSHGEOM_Limits lim;
      size_t np = (size_t)i;
      size_t ns = (size_t)j;
      pts[0] = grid_dp(0, 0);
      pts[1] = grid_dp(64, 0);
      pts[2] = grid_dp(64, 64);
      pts[3] = grid_dp(0, 64);
      subs[0].first = 0;
      subs[0].count = (size_t)(i > 1 ? 2 : 1);
      subs[0].closed = j % 2;
      subs[1].first = subs[0].count;
      subs[1].count = (size_t)(i > 2 ? i - 2 : 1);
      subs[1].closed = 1;
      /* the count j reaches 3: the third entry is read when the first two tile, and is one that never tiles (its first is past every
       * count here), so the answer is invalid, as it was when this entry was whatever the stack held */
      subs[2].first = SIZE_MAX;
      subs[2].count = 1;
      subs[2].closed = 0;
      gshgeom_limits_default(&lim);
      lim.max_points = (size_t)(1 + (i + j) % 3);
      grid_word(&f, (uint64_t)gshgeom_flat_path_from_points(np ? pts : NULL, np, ns ? subs : NULL, ns, NULL, NULL, &fp));
      if (fp) {
        grid_word(&f, gshgeom_flat_path_hash(fp));
        gshgeom_flat_path_destroy(fp);
        fp = NULL;
      }
      grid_word(&f, (uint64_t)gshgeom_flat_path_from_points(pts, np, subs, ns, NULL, &lim, &fp));
      if (fp) {
        grid_word(&f, gshgeom_flat_path_hash(fp));
        gshgeom_flat_path_destroy(fp);
      }
      grid_word(&f, (uint64_t)gshgeom_flat_path_from_points(pts, np, subs, ns, NULL, NULL, NULL));
    }
  }
  *words = f.n;
  return f.h;
}

#endif /* GHOTI_IO_GSHGEOM_TESTS_GRID_PROBE_H */
