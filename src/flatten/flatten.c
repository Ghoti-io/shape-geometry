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
 * The integer flattener (flatten.h) and the flattened path.
 *
 * Ported from `font`'s scan converter (`gfnt_raster_quad`, `gfnt_raster_cubic`
 * and the midpoint they use), whose recursion it keeps exactly: the tests
 * compare its output with the points `font`'s functions produced for the same
 * input, bit for bit.  Two things differ, neither visible in the output: the
 * deviation sums are 64-bit (`font`'s 32-bit ones overflow for a coordinate past
 * about 2^29), and a line is appended to the flattened path instead of being
 * deposited into scan converter cells.
 *
 * A failed call frees what it made and writes nothing to its output.
 */

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/shape-geometry/flatten.h>
#include <ghoti.io/shape-geometry/hash.h>

#include "../core/grow.h"

#include <string.h>

typedef struct Sub {
  size_t first;
  size_t count;
  uint8_t closed;
} Sub;

struct GSHGEOM_Flat_Path {
  GSHGEOM_Allocator alloc;
  size_t max_points;
  GSHGEOM_Device_Point * pt;
  size_t pt_n, pt_cap;
  Sub * sub;
  size_t sub_n, sub_cap;
};

typedef GSHGEOM_Device_Point Pt;

static GSHGEOM_Result add_point(GSHGEOM_Flat_Path * fp, Pt p) {
  if (fp->pt_n >= fp->max_points) {
    return GSHGEOM_ERR_LIMIT;
  }
  void * grown = fp->pt;
  GSHGEOM_Result r = gshgeom_grow(
      &fp->alloc, &grown, &fp->pt_cap, fp->pt_n + 1, sizeof *fp->pt);
  if (r != GSHGEOM_OK) {
    return r;
  }
  fp->pt = grown;
  fp->pt[fp->pt_n++] = p;
  fp->sub[fp->sub_n - 1].count++;
  return GSHGEOM_OK;
}

static GSHGEOM_Result begin_subpath(GSHGEOM_Flat_Path * fp, Pt p) {
  void * grown = fp->sub;
  GSHGEOM_Result r = gshgeom_grow(
      &fp->alloc, &grown, &fp->sub_cap, fp->sub_n + 1, sizeof *fp->sub);
  if (r != GSHGEOM_OK) {
    return r;
  }
  fp->sub = grown;
  fp->sub[fp->sub_n].first = fp->pt_n;
  fp->sub[fp->sub_n].count = 0;
  fp->sub[fp->sub_n].closed = 0;
  fp->sub_n++;
  r = add_point(fp, p);
  if (r != GSHGEOM_OK) {
    fp->sub_n--;
  }
  return r;
}

/* The midpoint of two points, floored (font's gfnt_raster_mid). */
static int32_t mid(int32_t a, int32_t b) {
  int64_t sum = (int64_t)a + b;
  int64_t q = sum / 2;
  if (sum % 2 != 0 && sum < 0) {
    q -= 1;
  }
  return (int32_t)q;
}

static Pt mid_pt(Pt a, Pt b) {
  Pt m;
  m.x = mid(a.x, b.x);
  m.y = mid(a.y, b.y);
  return m;
}

static int64_t abs64(int64_t v) {
  return v < 0 ? -v : v;
}

/* A quadratic's deviation from its chord is |2 p1 - p0 - p2| / 4 at its
 * midpoint, so the test is that sum against four tolerances. */
static GSHGEOM_Result quad(
    GSHGEOM_Flat_Path * fp, Pt p0, Pt p1, Pt p2, int depth) {
  int64_t dx = abs64(2 * (int64_t)p1.x - p0.x - p2.x);
  int64_t dy = abs64(2 * (int64_t)p1.y - p0.y - p2.y);
  if (depth <= 0 || dx + dy <= 4 * GSHGEOM_FLATTEN_TOLERANCE) {
    return add_point(fp, p2);
  }
  Pt a = mid_pt(p0, p1);
  Pt b = mid_pt(p1, p2);
  Pt m = mid_pt(a, b);
  GSHGEOM_Result r = quad(fp, p0, a, m, depth - 1);
  return r != GSHGEOM_OK ? r : quad(fp, m, b, p2, depth - 1);
}

/* A cubic's deviation is bounded by three quarters of the larger of its two
 * control-point deviations, so both are measured. */
static GSHGEOM_Result cubic(
    GSHGEOM_Flat_Path * fp, Pt p0, Pt p1, Pt p2, Pt p3, int depth) {
  int64_t d1 = abs64(3 * (int64_t)p1.x - 2 * (int64_t)p0.x - p3.x) +
      abs64(3 * (int64_t)p1.y - 2 * (int64_t)p0.y - p3.y);
  int64_t d2 = abs64(3 * (int64_t)p2.x - p0.x - 2 * (int64_t)p3.x) +
      abs64(3 * (int64_t)p2.y - p0.y - 2 * (int64_t)p3.y);
  if (depth <= 0 ||
      (d1 <= 4 * GSHGEOM_FLATTEN_TOLERANCE &&
          d2 <= 4 * GSHGEOM_FLATTEN_TOLERANCE)) {
    return add_point(fp, p3);
  }
  Pt a = mid_pt(p0, p1);
  Pt b = mid_pt(p1, p2);
  Pt c = mid_pt(p2, p3);
  Pt d = mid_pt(a, b);
  Pt e = mid_pt(b, c);
  Pt m = mid_pt(d, e);
  GSHGEOM_Result r = cubic(fp, p0, a, d, m, depth - 1);
  return r != GSHGEOM_OK ? r : cubic(fp, m, e, c, p3, depth - 1);
}

static void flat_destroy_parts(GSHGEOM_Flat_Path * fp) {
  GSHGEOM_Allocator alloc = fp->alloc;
  gcu_allocator_free(&alloc, fp->pt);
  gcu_allocator_free(&alloc, fp->sub);
  gcu_allocator_free(&alloc, fp);
}

void gshgeom_flat_path_destroy(GSHGEOM_Flat_Path * path) {
  if (path) {
    flat_destroy_parts(path);
  }
}

GSHGEOM_Result gshgeom_device_path_flatten(const GSHGEOM_Device_Path * path,
    const GSHGEOM_Allocator * alloc, const GSHGEOM_Limits * limits,
    GSHGEOM_Flat_Path ** out) {
  if (!path || !out) {
    return GSHGEOM_ERR_INVALID;
  }
  if (!alloc) {
    alloc = gshgeom_allocator_default();
  }
  GSHGEOM_Flat_Path * fp = gcu_allocator_calloc(alloc, 1, sizeof *fp);
  if (!fp) {
    return GSHGEOM_ERR_OOM;
  }
  fp->alloc = *alloc;
  fp->max_points = (limits && limits->max_points)
      ? limits->max_points
      : GSHGEOM_LIMITS_DEFAULT_MAX_POINTS;
  GSHGEOM_Result r = GSHGEOM_OK;
  Pt current = {0, 0};
  size_t n = gshgeom_device_path_element_count(path);
  for (size_t i = 0; i < n && r == GSHGEOM_OK; i++) {
    GSHGEOM_Device_Element el;
    if (gshgeom_device_path_element(path, i, &el) != GSHGEOM_OK) {
      r = GSHGEOM_ERR_INTERNAL;
      break;
    }
    /* The builder guarantees a move first and after every close. */
    if (el.verb != GSHGEOM_DEVICE_MOVE && fp->sub_n == 0) {
      r = GSHGEOM_ERR_INTERNAL;
      break;
    }
    switch (el.verb) {
    case GSHGEOM_DEVICE_MOVE:
      r = begin_subpath(fp, el.p[0]);
      current = el.p[0];
      break;
    case GSHGEOM_DEVICE_LINE:
      r = add_point(fp, el.p[0]);
      current = el.p[0];
      break;
    case GSHGEOM_DEVICE_QUADRATIC:
      r = quad(fp, current, el.p[0], el.p[1], GSHGEOM_FLATTEN_MAX_DEPTH);
      current = el.p[1];
      break;
    case GSHGEOM_DEVICE_CUBIC:
      r = cubic(fp, current, el.p[0], el.p[1], el.p[2], GSHGEOM_FLATTEN_MAX_DEPTH);
      current = el.p[2];
      break;
    case GSHGEOM_DEVICE_CLOSE:
      fp->sub[fp->sub_n - 1].closed = 1;
      break;
    case GSHGEOM_DEVICE_VERB_COUNT:
      r = GSHGEOM_ERR_INTERNAL;
      break;
    }
  }
  if (r != GSHGEOM_OK) {
    flat_destroy_parts(fp);
    return r;
  }
  *out = fp;
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_path_flatten(const GSHGEOM_Path * path,
    const GSHGEOM_Map * map, const GSHGEOM_Allocator * alloc,
    const GSHGEOM_Limits * limits, GSHGEOM_Flat_Path ** out) {
  if (!out) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Device_Path * dp;
  GSHGEOM_Result r = gshgeom_device_path_from_model(path, map, alloc, limits, &dp);
  if (r != GSHGEOM_OK) {
    return r;
  }
  r = gshgeom_device_path_flatten(dp, alloc, limits, out);
  gshgeom_device_path_destroy(dp);
  return r;
}

/* ---------------------------------------------------------------------- */
/* Readback                                                                */
/* ---------------------------------------------------------------------- */

size_t gshgeom_flat_path_point_count(const GSHGEOM_Flat_Path * path) {
  return path ? path->pt_n : 0;
}

size_t gshgeom_flat_path_subpath_count(const GSHGEOM_Flat_Path * path) {
  return path ? path->sub_n : 0;
}

const GSHGEOM_Device_Point * gshgeom_flat_path_points(
    const GSHGEOM_Flat_Path * path) {
  return (path && path->pt_n) ? path->pt : NULL;
}

GSHGEOM_Result gshgeom_flat_path_subpath(
    const GSHGEOM_Flat_Path * path, size_t index, GSHGEOM_Flat_Subpath * out) {
  if (!path || !out || index >= path->sub_n) {
    return GSHGEOM_ERR_INVALID;
  }
  out->first = path->sub[index].first;
  out->count = path->sub[index].count;
  out->closed = path->sub[index].closed;
  return GSHGEOM_OK;
}

/* ---------------------------------------------------------------------- */
/* Hash and dump                                                           */
/* ---------------------------------------------------------------------- */

/** The format version of the hash and the dump. */
#define FLAT_FORMAT_VERSION 1u
/** The kind field: "FLAT" as little-endian bytes. */
#define FLAT_KIND 0x54414c46u

uint64_t gshgeom_flat_path_hash(const GSHGEOM_Flat_Path * path) {
  GSHGEOM_Hash h;
  gshgeom_hash_init(&h);
  if (!path) {
    return gshgeom_hash_value(&h);
  }
  gshgeom_hash_u32(&h, FLAT_FORMAT_VERSION);
  gshgeom_hash_u32(&h, FLAT_KIND);
  gshgeom_hash_u64(&h, (uint64_t)path->sub_n);
  for (size_t i = 0; i < path->sub_n; i++) {
    const Sub * s = &path->sub[i];
    gshgeom_hash_u8(&h, s->closed);
    gshgeom_hash_u64(&h, (uint64_t)s->count);
    for (size_t k = 0; k < s->count; k++) {
      gshgeom_hash_i32(&h, path->pt[s->first + k].x);
      gshgeom_hash_i32(&h, path->pt[s->first + k].y);
    }
  }
  return gshgeom_hash_value(&h);
}

GSHGEOM_Result gshgeom_flat_path_dump(
    const GSHGEOM_Flat_Path * path, FILE * out) {
  if (!path || !out) {
    return GSHGEOM_ERR_INVALID;
  }
  if (fprintf(out, "gshgeom.flat-path %u subpaths=%llu points=%llu hash=%016llx\n",
          FLAT_FORMAT_VERSION, (unsigned long long)path->sub_n,
          (unsigned long long)path->pt_n,
          (unsigned long long)gshgeom_flat_path_hash(path)) < 0) {
    return GSHGEOM_ERR_IO;
  }
  for (size_t i = 0; i < path->sub_n; i++) {
    const Sub * s = &path->sub[i];
    if (fprintf(out, "S closed=%d points=%llu\n", s->closed ? 1 : 0,
            (unsigned long long)s->count) < 0) {
      return GSHGEOM_ERR_IO;
    }
    for (size_t k = 0; k < s->count; k++) {
      if (fprintf(out, "P %08x %08x\n",
              (unsigned)(uint32_t)path->pt[s->first + k].x,
              (unsigned)(uint32_t)path->pt[s->first + k].y) < 0) {
        return GSHGEOM_ERR_IO;
      }
    }
  }
  return fflush(out) == 0 ? GSHGEOM_OK : GSHGEOM_ERR_IO;
}
