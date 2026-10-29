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
 * The rounding rule, the device path (builder, readback, hash and dump) and the
 * conversion of a model path under a map into it (device.h).
 *
 * Elements store their points in one array, in readback order: control points
 * then the end point, and a close stores the subpath's start.  The builders
 * reserve everything an operation needs before they store anything, so a
 * failure changes capacity and nothing else.
 */

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/math/basic.h>
#include <ghoti.io/shape-geometry/device.h>
#include <ghoti.io/shape-geometry/hash.h>

#include "../core/grow.h"

#include <string.h>

/* One pixel in device units: the factor of the rule (an exact power of two). */
#define SCALE ((double)GSHGEOM_DEVICE_ONE)

/* The bounds of the rule, on t = x * 64 (exact): floor(t + 0.5) fits an int32
 * exactly when -2^31 - 0.5 <= t < 2^31 - 0.5. */
#define ROUND_T_MAX 2147483647.5
#define ROUND_T_MIN (-2147483648.5)

/**
 * floor(t + 0.5) for a finite, already scaled t, without ever adding 0.5: the
 * floor and the fraction are both exact, and the half-way point is a
 * comparison.  Infinite t (a finite x so large that x * 64 overflowed) falls to
 * the range check.
 */
static GSHGEOM_Result round_scaled(double t, int32_t * out) {
  if (!(t < ROUND_T_MAX) || t < ROUND_T_MIN) {
    return GSHGEOM_ERR_LIMIT;
  }
  double f = gmath_floor(t);
  double r = (t - f) >= 0.5 ? f + 1.0 : f;
  *out = (int32_t)r;
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_round_26_6(double x, int32_t * out) {
  if (!out || !((x - x) == 0.0)) {
    return GSHGEOM_ERR_INVALID;
  }
  return round_scaled(x * SCALE, out);
}

GSHGEOM_Result gshgeom_scale_26_6(
    int32_t p, int32_t scale, int flip, int32_t * out) {
  if (!out || (flip != 0 && flip != 1)) {
    return GSHGEOM_ERR_INVALID;
  }
  /* |p * scale| <= 2^62: the product and its negation fit an int64_t. */
  int64_t v = (int64_t)p * (int64_t)scale;
  if (flip) {
    v = -v;
  }
  /* floor((v + 2^15) / 2^16) by shifts of a magnitude, never of a negative
   * number (whose right shift C leaves to the implementation). */
  int64_t n = v + 32768;
  int64_t q = n >= 0 ? (int64_t)((uint64_t)n >> 16)
                     : -(int64_t)(((uint64_t)(-n) + 65535u) >> 16);
  if (q > INT32_MAX || q < INT32_MIN) {
    return GSHGEOM_ERR_LIMIT;
  }
  *out = (int32_t)q;
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_device_point(const GSHGEOM_Map * map, GSHGEOM_Point p,
    GSHGEOM_Device_Point * out) {
  if (!map || !out) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Result r = gshgeom_map_check(map);
  if (r != GSHGEOM_OK) {
    return r;
  }
  if (map->dim != 2) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Point m;
  r = gshgeom_map_apply_point(map, p, &m);
  if (r != GSHGEOM_OK) {
    return r;
  }
  GSHGEOM_Device_Point d;
  r = gshgeom_round_26_6(m.x, &d.x);
  if (r == GSHGEOM_OK) {
    r = gshgeom_round_26_6(m.y, &d.y);
  }
  if (r == GSHGEOM_OK) {
    *out = d;
  }
  return r;
}

/* ---------------------------------------------------------------------- */
/* The device path                                                         */
/* ---------------------------------------------------------------------- */

typedef struct Elem {
  uint8_t verb;
  size_t point; /* index of its first stored point */
} Elem;

struct GSHGEOM_Device_Path {
  GSHGEOM_Allocator alloc;
  size_t max_elements;
  Elem * el;
  size_t el_n, el_cap;
  GSHGEOM_Device_Point * pt;
  size_t pt_n, pt_cap;
  int has_current; /* a move has been made */
  int closed;      /* the last element is a close */
  GSHGEOM_Device_Point start;
};

static unsigned verb_points(int verb) {
  switch (verb) {
  case GSHGEOM_DEVICE_MOVE:
  case GSHGEOM_DEVICE_LINE:
  case GSHGEOM_DEVICE_CLOSE:
    return 1;
  case GSHGEOM_DEVICE_QUADRATIC:
    return 2;
  case GSHGEOM_DEVICE_CUBIC:
    return 3;
  default:
    return 0;
  }
}

GSHGEOM_Result gshgeom_device_path_create(const GSHGEOM_Allocator * alloc,
    const GSHGEOM_Limits * limits, GSHGEOM_Device_Path ** out) {
  if (!out) {
    return GSHGEOM_ERR_INVALID;
  }
  if (!alloc) {
    alloc = gshgeom_allocator_default();
  }
  GSHGEOM_Device_Path * path = gcu_allocator_calloc(alloc, 1, sizeof *path);
  if (!path) {
    return GSHGEOM_ERR_OOM;
  }
  path->alloc = *alloc;
  path->max_elements = (limits && limits->max_elements)
      ? limits->max_elements
      : GSHGEOM_LIMITS_DEFAULT_MAX_ELEMENTS;
  *out = path;
  return GSHGEOM_OK;
}

void gshgeom_device_path_destroy(GSHGEOM_Device_Path * path) {
  if (!path) {
    return;
  }
  GSHGEOM_Allocator alloc = path->alloc;
  gcu_allocator_free(&alloc, path->el);
  gcu_allocator_free(&alloc, path->pt);
  gcu_allocator_free(&alloc, path);
}

size_t gshgeom_device_path_element_count(const GSHGEOM_Device_Path * path) {
  return path ? path->el_n : 0;
}

/* Add `verb` with its points, preceded by a move to the closed subpath's start
 * when `with_move` (a drawing verb after a close).  All or nothing. */
static GSHGEOM_Result append(GSHGEOM_Device_Path * path, int verb,
    const GSHGEOM_Device_Point * pts, int with_move) {
  size_t np = verb_points(verb);
  size_t add_el = with_move ? 2 : 1;
  size_t add_pt = np + (with_move ? 1 : 0);
  if (path->el_n > path->max_elements ||
      add_el > path->max_elements - path->el_n) {
    return GSHGEOM_ERR_LIMIT;
  }
  void * grown = path->el;
  GSHGEOM_Result r = gshgeom_grow(&path->alloc, &grown, &path->el_cap,
      path->el_n + add_el, sizeof *path->el);
  if (r != GSHGEOM_OK) {
    return r;
  }
  path->el = grown;
  grown = path->pt;
  r = gshgeom_grow(&path->alloc, &grown, &path->pt_cap, path->pt_n + add_pt,
      sizeof *path->pt);
  if (r != GSHGEOM_OK) {
    return r;
  }
  path->pt = grown;
  if (with_move) {
    path->el[path->el_n].verb = GSHGEOM_DEVICE_MOVE;
    path->el[path->el_n].point = path->pt_n;
    path->pt[path->pt_n] = path->start;
    path->el_n++;
    path->pt_n++;
  }
  path->el[path->el_n].verb = (uint8_t)verb;
  path->el[path->el_n].point = path->pt_n;
  for (size_t i = 0; i < np; i++) {
    path->pt[path->pt_n + i] = pts[i];
  }
  path->el_n++;
  path->pt_n += np;
  path->closed = verb == GSHGEOM_DEVICE_CLOSE;
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_device_path_move_to(
    GSHGEOM_Device_Path * path, GSHGEOM_Device_Point p) {
  if (!path) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Result r = append(path, GSHGEOM_DEVICE_MOVE, &p, 0);
  if (r == GSHGEOM_OK) {
    path->has_current = 1;
    path->start = p;
  }
  return r;
}

GSHGEOM_Result gshgeom_device_path_line_to(
    GSHGEOM_Device_Path * path, GSHGEOM_Device_Point p) {
  if (!path || !path->has_current) {
    return GSHGEOM_ERR_INVALID;
  }
  return append(path, GSHGEOM_DEVICE_LINE, &p, path->closed);
}

GSHGEOM_Result gshgeom_device_path_quadratic_to(GSHGEOM_Device_Path * path,
    GSHGEOM_Device_Point c, GSHGEOM_Device_Point p) {
  if (!path || !path->has_current) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Device_Point pts[2] = {c, p};
  return append(path, GSHGEOM_DEVICE_QUADRATIC, pts, path->closed);
}

GSHGEOM_Result gshgeom_device_path_cubic_to(GSHGEOM_Device_Path * path,
    GSHGEOM_Device_Point c1, GSHGEOM_Device_Point c2, GSHGEOM_Device_Point p) {
  if (!path || !path->has_current) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Device_Point pts[3] = {c1, c2, p};
  return append(path, GSHGEOM_DEVICE_CUBIC, pts, path->closed);
}

GSHGEOM_Result gshgeom_device_path_close(GSHGEOM_Device_Path * path) {
  if (!path || !path->has_current) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Device_Point start = path->start;
  return append(path, GSHGEOM_DEVICE_CLOSE, &start, 0);
}

GSHGEOM_Result gshgeom_device_path_element(const GSHGEOM_Device_Path * path,
    size_t index, GSHGEOM_Device_Element * out) {
  if (!path || !out || index >= path->el_n) {
    return GSHGEOM_ERR_INVALID;
  }
  const Elem * e = &path->el[index];
  GSHGEOM_Device_Element el;
  memset(&el, 0, sizeof el);
  el.verb = (GSHGEOM_Device_Verb)e->verb;
  el.count = verb_points(e->verb);
  for (unsigned i = 0; i < el.count; i++) {
    el.p[i] = path->pt[e->point + i];
  }
  *out = el;
  return GSHGEOM_OK;
}

/* ---------------------------------------------------------------------- */
/* A model path under a map                                                */
/* ---------------------------------------------------------------------- */

#define RATIONAL_MAX_DEPTH 16
/* A piece of a rational curve is flat when its control points are within 2
 * units (1/32 pixel) of its chord; squared. */
#define RATIONAL_TOLERANCE_SQ 4.0

/* A point in homogeneous coordinates: (w X, w Y, w). */
typedef struct Hp {
  double x, y, w;
} Hp;

static Hp hp_mid(Hp a, Hp b) {
  Hp m;
  m.x = 0.5 * (a.x + b.x);
  m.y = 0.5 * (a.y + b.y);
  m.w = 0.5 * (a.w + b.w);
  return m;
}

/* The squared distance from (px, py) to the segment (sx, sy) - (ex, ey). */
static double dist2_to_segment(
    double px, double py, double sx, double sy, double ex, double ey) {
  double dx = ex - sx;
  double dy = ey - sy;
  double ax = px - sx;
  double ay = py - sy;
  double len2 = dx * dx + dy * dy;
  if (!(len2 > 0.0)) {
    return ax * ax + ay * ay;
  }
  double u = (ax * dx + ay * dy) / len2;
  if (u < 0.0) {
    u = 0.0;
  } else if (u > 1.0) {
    u = 1.0;
  }
  double qx = ax - u * dx;
  double qy = ay - u * dy;
  return qx * qx + qy * qy;
}

/* Lines for the rational piece h[0..n-1] (n is 3 or 4) whose projected end
 * points are (sx, sy) and (ex, ey): the caller's exact values. */
static GSHGEOM_Result rational_lines(GSHGEOM_Device_Path * dp, int n,
    const Hp * h, double sx, double sy, double ex, double ey, int depth) {
  int flat = 1;
  for (int i = 1; i < n - 1; i++) {
    double px = h[i].x / h[i].w;
    double py = h[i].y / h[i].w;
    if (dist2_to_segment(px, py, sx, sy, ex, ey) > RATIONAL_TOLERANCE_SQ) {
      flat = 0;
      break;
    }
  }
  if (depth <= 0 || flat) {
    GSHGEOM_Device_Point d;
    GSHGEOM_Result r = round_scaled(ex, &d.x);
    if (r == GSHGEOM_OK) {
      r = round_scaled(ey, &d.y);
    }
    return r != GSHGEOM_OK ? r : gshgeom_device_path_line_to(dp, d);
  }
  Hp left[4], right[4];
  if (n == 3) {
    Hp a = hp_mid(h[0], h[1]);
    Hp b = hp_mid(h[1], h[2]);
    Hp m = hp_mid(a, b);
    left[0] = h[0]; left[1] = a; left[2] = m;
    right[0] = m; right[1] = b; right[2] = h[2];
    double mx = m.x / m.w;
    double my = m.y / m.w;
    GSHGEOM_Result r = rational_lines(dp, 3, left, sx, sy, mx, my, depth - 1);
    return r != GSHGEOM_OK
        ? r
        : rational_lines(dp, 3, right, mx, my, ex, ey, depth - 1);
  }
  Hp a = hp_mid(h[0], h[1]);
  Hp b = hp_mid(h[1], h[2]);
  Hp c = hp_mid(h[2], h[3]);
  Hp d = hp_mid(a, b);
  Hp e = hp_mid(b, c);
  Hp m = hp_mid(d, e);
  left[0] = h[0]; left[1] = a; left[2] = d; left[3] = m;
  right[0] = m; right[1] = e; right[2] = c; right[3] = h[3];
  double mx = m.x / m.w;
  double my = m.y / m.w;
  GSHGEOM_Result r = rational_lines(dp, 4, left, sx, sy, mx, my, depth - 1);
  return r != GSHGEOM_OK
      ? r
      : rational_lines(dp, 4, right, mx, my, ex, ey, depth - 1);
}

/* Map a model point to device units, scaled but not rounded (x * 64 is exact),
 * and check that it rounds into range. */
static GSHGEOM_Result map_scaled(const GSHGEOM_Map * map, GSHGEOM_Point p,
    double * x, double * y) {
  GSHGEOM_Point m;
  GSHGEOM_Result r = gshgeom_map_apply_point(map, p, &m);
  if (r != GSHGEOM_OK) {
    return r;
  }
  double sx = m.x * SCALE;
  double sy = m.y * SCALE;
  int32_t probe;
  r = round_scaled(sx, &probe);
  if (r == GSHGEOM_OK) {
    r = round_scaled(sy, &probe);
  }
  if (r == GSHGEOM_OK) {
    *x = sx;
    *y = sy;
  }
  return r;
}

static GSHGEOM_Result map_rounded(const GSHGEOM_Map * map, GSHGEOM_Point p,
    GSHGEOM_Device_Point * out) {
  double x, y;
  GSHGEOM_Result r = map_scaled(map, p, &x, &y);
  if (r == GSHGEOM_OK) {
    r = round_scaled(x, &out->x);
  }
  if (r == GSHGEOM_OK) {
    r = round_scaled(y, &out->y);
  }
  return r;
}

static GSHGEOM_Result convert_element(GSHGEOM_Device_Path * dp,
    const GSHGEOM_Map * map, const GSHGEOM_Element * el) {
  const GSHGEOM_Curve * c = &el->curve;
  GSHGEOM_Device_Point q[3];
  GSHGEOM_Result r = GSHGEOM_OK;
  switch (el->verb) {
  case GSHGEOM_VERB_MOVE:
    r = map_rounded(map, c->p[1], &q[0]);
    return r != GSHGEOM_OK ? r : gshgeom_device_path_move_to(dp, q[0]);
  case GSHGEOM_VERB_LINE:
    r = map_rounded(map, c->p[1], &q[0]);
    return r != GSHGEOM_OK ? r : gshgeom_device_path_line_to(dp, q[0]);
  case GSHGEOM_VERB_QUADRATIC:
    r = map_rounded(map, c->p[1], &q[0]);
    if (r == GSHGEOM_OK) {
      r = map_rounded(map, c->p[2], &q[1]);
    }
    return r != GSHGEOM_OK ? r : gshgeom_device_path_quadratic_to(dp, q[0], q[1]);
  case GSHGEOM_VERB_CUBIC:
    r = map_rounded(map, c->p[1], &q[0]);
    if (r == GSHGEOM_OK) {
      r = map_rounded(map, c->p[2], &q[1]);
    }
    if (r == GSHGEOM_OK) {
      r = map_rounded(map, c->p[3], &q[2]);
    }
    return r != GSHGEOM_OK
        ? r
        : gshgeom_device_path_cubic_to(dp, q[0], q[1], q[2]);
  case GSHGEOM_VERB_RATIONAL_QUADRATIC:
  case GSHGEOM_VERB_RATIONAL_CUBIC: {
    int n = el->verb == GSHGEOM_VERB_RATIONAL_QUADRATIC ? 3 : 4;
    Hp h[4];
    double sx[4], sy[4];
    for (int i = 0; i < n && r == GSHGEOM_OK; i++) {
      r = map_scaled(map, c->p[i], &sx[i], &sy[i]);
      h[i].x = c->w[i] * sx[i];
      h[i].y = c->w[i] * sy[i];
      h[i].w = c->w[i];
    }
    if (r != GSHGEOM_OK) {
      return r;
    }
    return rational_lines(dp, n, h, sx[0], sy[0], sx[n - 1], sy[n - 1],
        RATIONAL_MAX_DEPTH);
  }
  case GSHGEOM_VERB_CLOSE:
    return gshgeom_device_path_close(dp);
  case GSHGEOM_VERB_ARC:
    return GSHGEOM_ERR_UNSUPPORTED; /* the device route of an arc is not defined yet; never its chord */
  case GSHGEOM_VERB_COUNT:
    break;
  }
  return GSHGEOM_ERR_INTERNAL;
}

GSHGEOM_Result gshgeom_device_path_from_model(const GSHGEOM_Path * path,
    const GSHGEOM_Map * map, const GSHGEOM_Allocator * alloc,
    const GSHGEOM_Limits * limits, GSHGEOM_Device_Path ** out) {
  if (!path || !map || !out || gshgeom_path_dimension(path) != 2) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Result r = gshgeom_map_check(map);
  if (r != GSHGEOM_OK) {
    return r;
  }
  if (map->dim != 2) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Device_Path * dp;
  r = gshgeom_device_path_create(alloc, limits, &dp);
  if (r != GSHGEOM_OK) {
    return r;
  }
  size_t n = gshgeom_path_element_count(path);
  for (size_t i = 0; i < n && r == GSHGEOM_OK; i++) {
    GSHGEOM_Element el;
    r = gshgeom_path_element(path, i, &el);
    if (r != GSHGEOM_OK) {
      r = r == GSHGEOM_ERR_UNSUPPORTED ? r : GSHGEOM_ERR_INTERNAL; /* an arc is refused, never drawn as its chord */
      break;
    }
    r = convert_element(dp, map, &el);
  }
  if (r != GSHGEOM_OK) {
    gshgeom_device_path_destroy(dp);
    return r;
  }
  *out = dp;
  return GSHGEOM_OK;
}

/* ---------------------------------------------------------------------- */
/* Hash and dump                                                           */
/* ---------------------------------------------------------------------- */

/** The format version of the hash and the dump. */
#define DEVICE_FORMAT_VERSION 1u
/** The kind field: "DEVP" as little-endian bytes. */
#define DEVICE_KIND 0x50564544u

uint64_t gshgeom_device_path_hash(const GSHGEOM_Device_Path * path) {
  GSHGEOM_Hash h;
  gshgeom_hash_init(&h);
  if (!path) {
    return gshgeom_hash_value(&h);
  }
  gshgeom_hash_u32(&h, DEVICE_FORMAT_VERSION);
  gshgeom_hash_u32(&h, DEVICE_KIND);
  gshgeom_hash_u64(&h, (uint64_t)path->el_n);
  for (size_t i = 0; i < path->el_n; i++) {
    const Elem * e = &path->el[i];
    gshgeom_hash_u8(&h, e->verb);
    if (e->verb == GSHGEOM_DEVICE_CLOSE) {
      continue;
    }
    for (unsigned k = 0; k < verb_points(e->verb); k++) {
      gshgeom_hash_i32(&h, path->pt[e->point + k].x);
      gshgeom_hash_i32(&h, path->pt[e->point + k].y);
    }
  }
  return gshgeom_hash_value(&h);
}

static int put_i32(FILE * out, int32_t v) {
  return fprintf(out, " %08x", (unsigned)(uint32_t)v);
}

GSHGEOM_Result gshgeom_device_path_dump(
    const GSHGEOM_Device_Path * path, FILE * out) {
  static const char names[GSHGEOM_DEVICE_VERB_COUNT] = {'M', 'L', 'Q', 'C', 'Z'};
  if (!path || !out) {
    return GSHGEOM_ERR_INVALID;
  }
  /* The count goes through an unsigned long long: the width of size_t is a
   * property of the host. */
  if (fprintf(out, "gshgeom.device-path %u elements=%llu hash=%016llx\n",
          DEVICE_FORMAT_VERSION, (unsigned long long)path->el_n,
          (unsigned long long)gshgeom_device_path_hash(path)) < 0) {
    return GSHGEOM_ERR_IO;
  }
  for (size_t i = 0; i < path->el_n; i++) {
    const Elem * e = &path->el[i];
    if (fputc(names[e->verb], out) == EOF) {
      return GSHGEOM_ERR_IO;
    }
    if (e->verb != GSHGEOM_DEVICE_CLOSE) {
      for (unsigned k = 0; k < verb_points(e->verb); k++) {
        if (put_i32(out, path->pt[e->point + k].x) < 0 ||
            put_i32(out, path->pt[e->point + k].y) < 0) {
          return GSHGEOM_ERR_IO;
        }
      }
    }
    if (fputc('\n', out) == EOF) {
      return GSHGEOM_ERR_IO;
    }
  }
  /* A buffered stream reports a write error at the flush, and a stream's error flag is sticky:
   * a failed write whose return value was lost (an unbuffered stream's fprintf) still sets it. */
  return (fflush(out) == 0 && !ferror(out)) ? GSHGEOM_OK : GSHGEOM_ERR_IO;
}
