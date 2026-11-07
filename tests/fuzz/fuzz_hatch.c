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
 * Hatches (hatch.h, HX2 and HX4 to HX12) driven by bytes: one entry point for libFuzzer and for the replay gate
 * (`make fuzz-replay`).
 *
 * Header bytes choose the limits (a small `max_points`, a small line limit), the allocation that fails, the style (width
 * space, cap, `degenerate`) and the fill rule.  Then the hatch map and the chain (identity, scale, rotation, a singular map or
 * raw values), one to three families (each field from a table of degenerate values or from the bytes, the elements of up to
 * five dashes), the nesting and the angle unit, and one to three contours (rectangles and polygons).  After each call the
 * library's stated properties are checked, and a violation is an `abort`:
 *
 * - the result is one the header lists, a refusal returns nothing and the counting allocator holds nothing live at the
 *   end, and ::gshgeom_hatch_check says the same as the call when the arguments are the cause (HX2);
 * - the lines of the use are at most the limit, and the limit is the answer when one more would be made (HX12);
 * - every piece is finite, a dot has both points equal, and the pieces are the same on a second call (the hash);
 * - the expansion is a stroke's: the same status codes, no aligned or clip part, and, when no piece is a dot under a square cap
 *   or under `degenerate` omit (a dot is a dash of the pattern and is drawn along its line whatever the style says), the hash of
 *   the stroke expansion of the pieces as an explicit path (HX8).
 *
 * The replay fails if the corpus never reaches an outcome: tools/gen-fuzz-corpus.py writes it.
 */

#include <ghoti.io/shape-geometry/shape-geometry.h>

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(cond) \
  do { \
    if (!(cond)) { \
      fprintf(stderr, "fuzz_hatch: invariant failed: %s (line %d)\n", #cond, __LINE__); \
      abort(); \
    } \
  } while (0)

enum {
  OUT_PIECES = 0, OUT_DOT, OUT_DASHED, OUT_NESTED, OUT_TWO_FAMILIES, OUT_LINE_LIMIT, OUT_POINT_LIMIT, OUT_OOM, OUT_INVALID,
  OUT_UNSUPPORTED, OUT_NOT_INVERTIBLE, OUT_EMPTY_WIDTH, OUT_NO_PIECES, OUT_DIFFERENTIAL, OUT_DEVICE_WIDTH, OUT_COUNT
};

static const char * const outcome_names[OUT_COUNT] = {"pieces made", "a dot", "a dashed family", "outer or ignore nesting",
    "two families", "the lines limit", "a part past max_points", "a failed allocation", "an invalid hatch",
    "an unsupported (3D) map", "a singular map", "width 0", "a region with no pieces", "the differential compared",
    "device width"};

static unsigned long outcomes[OUT_COUNT];

int gshgeom_fuzz_unreached(void);

int gshgeom_fuzz_unreached(void) {
  int missing = 0;
  for (int i = 0; i < OUT_COUNT; i++) {
    if (outcomes[i] == 0) {
      fprintf(stderr, "fuzz_hatch: the corpus never reached: %s\n", outcome_names[i]);
      missing++;
    }
  }
  return missing;
}

typedef struct Counter {
  long live;
  long calls;
  long fail_at;
  int failed;
} Counter;

static int fails(Counter * c) {
  if (c->fail_at >= 0 && c->calls++ == c->fail_at) {
    c->failed = 1;
    return 1;
  }
  return 0;
}

static void * c_malloc(void * ctx, size_t n) {
  Counter * c = ctx;
  void * p = fails(c) ? NULL : malloc(n ? n : 1);
  c->live += p != NULL;
  return p;
}

static void * c_calloc(void * ctx, size_t k, size_t n) {
  Counter * c = ctx;
  void * p = fails(c) ? NULL : calloc(k ? k : 1, n ? n : 1);
  c->live += p != NULL;
  return p;
}

static void * c_realloc(void * ctx, void * old, size_t n) {
  Counter * c = ctx;
  void * p = fails(c) ? NULL : realloc(old, n ? n : 1);
  c->live += p != NULL && old == NULL;
  return p;
}

static void c_free(void * ctx, void * p) {
  ((Counter *)ctx)->live -= p != NULL;
  free(p);
}

typedef struct Input {
  const uint8_t * data;
  size_t size;
  size_t at;
} Input;

static uint8_t next(Input * in) {
  return in->at < in->size ? in->data[in->at++] : 0;
}

/* A value: a degenerate one from the table, or a 16-bit value over 8 from two bytes. */
static double value(Input * in) {
  static const double table[] = {0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 7.5, 10.0, -1.0, -2.0, 0x1p-70, 0x1p70, 0x1p-64, 0x1p64,
                                 90.0, 45.0, 30.0, 0.0, 1e300, -0.0, 5.0, 12.0};
  const uint8_t pick = next(in);
  if (pick < 24) {
    return table[pick];
  }
  if (pick == 24) {
    return NAN;
  }
  if (pick == 25) {
    return INFINITY;
  }
  const unsigned hi = next(in);
  const unsigned lo = next(in);
  return (double)(int16_t)(uint16_t)((hi << 8) | lo) / 8.0;
}

static GSHGEOM_Map map_of(Input * in) {
  GSHGEOM_Map m;
  const uint8_t kind = next(in) % 7u;
  GSHGEOM_Result r = GSHGEOM_OK;
  switch (kind) {
    case 0:
      r = gshgeom_map_identity(2, &m);
      break;
    case 1:
      r = gshgeom_map_scale(2, 2.0, 2.0, 1.0, &m);
      break;
    case 2:
      r = gshgeom_map_rotate(GSHGEOM_ANGLE_DEGREES, 90.0, &m);
      break;
    case 3:
      r = gshgeom_map_scale(2, 0.0, 1.0, 1.0, &m); /* singular */
      break;
    case 4: {
      const double lin[4] = {value(in), value(in), value(in), value(in)};
      const double t[2] = {value(in), value(in)};
      r = gshgeom_map_from_matrix(2, lin, t, &m);
      break;
    }
    case 5:
      r = gshgeom_map_identity(3, &m); /* 3D */
      break;
    default:
      r = gshgeom_map_flip_y(2, &m);
      break;
  }
  if (r != GSHGEOM_OK) {
    gshgeom_map_identity(2, &m);
  }
  return m;
}

static int is_square(const GSHGEOM_Stroke * st) {
  return st->cap == GSHGEOM_STROKE_CAP_SQUARE;
}

int LLVMFuzzerTestOneInput(const uint8_t * data, size_t size);

int LLVMFuzzerTestOneInput(const uint8_t * data, size_t size) {
  Input in = {data, size, 0};
  const uint8_t opts = next(&in);
  const uint8_t fail = next(&in);
  const int small_points = (opts & 1) != 0;
  const int fail_alloc = (opts & 2) != 0;
  const int small_lines = (opts & 4) != 0;
  const int device = (opts & 8) != 0;
  const unsigned cap = (opts >> 4) & 3u;
  const int omit = (opts & 64) != 0;
  const GSHGEOM_Fill_Rule rule = (opts & 128) ? GSHGEOM_FILL_EVEN_ODD : GSHGEOM_FILL_NONZERO;

  GSHGEOM_Map hmap = map_of(&in);
  GSHGEOM_Map chain = map_of(&in);
  GSHGEOM_Hatch hatch;
  memset(&hatch, 0, sizeof hatch);
  hatch.map = hmap;
  hatch.angle_unit = next(&in) % 5u;
  hatch.nesting = next(&in) % 4u;
  GSHGEOM_Hatch_Family fam[3];
  double dashes[3][5];
  const size_t nf = next(&in) % 4u;
  for (size_t i = 0; i < nf; i++) {
    GSHGEOM_Hatch_Family * f = &fam[i];
    memset(f, 0, sizeof *f);
    f->angle = value(&in);
    f->origin_x = value(&in);
    f->origin_y = value(&in);
    f->offset_along = value(&in);
    f->offset_across = value(&in);
    const size_t nd = next(&in) % 6u;
    for (size_t k = 0; k < nd; k++) {
      dashes[i][k] = value(&in);
    }
    f->dash_count = nd;
    f->dashes = nd ? dashes[i] : NULL;
  }
  hatch.family_count = nf;
  hatch.families = fam;

  GSHGEOM_Stroke st;
  gshgeom_stroke_init(&st);
  st.width = value(&in);
  st.width_space = device ? GSHGEOM_STROKE_WIDTH_DEVICE : GSHGEOM_STROKE_WIDTH_LOCAL;
  st.cap = cap;
  st.degenerate = omit ? GSHGEOM_STROKE_DEGENERATE_OMIT : GSHGEOM_STROKE_DEGENERATE_DRAW;

  GSHGEOM_Path * region = NULL;
  CHECK(gshgeom_path_create(NULL, NULL, 2, &region) == GSHGEOM_OK);
  const size_t nc = 1u + next(&in) % 3u;
  for (size_t c = 0; c < nc; c++) {
    const uint8_t kind = next(&in) % 3u;
    const double x = value(&in), y = value(&in), w = value(&in), h = value(&in);
    if (kind == 0) {
      (void)gshgeom_path_add_rect(region, x, y, w, h, c & 1, 0);
    } else {
      GSHGEOM_Point p = {x, y, 0.0};
      if (gshgeom_path_move_to(region, p) == GSHGEOM_OK) {
        GSHGEOM_Point q = {x + w, y, 0.0}, s = {x + w / 2.0, y + h, 0.0};
        (void)gshgeom_path_line_to(region, q);
        (void)gshgeom_path_line_to(region, s);
        if (kind == 2) {
          (void)gshgeom_path_close(region);
        }
      }
    }
  }

  Counter counter = {0, 0, fail_alloc ? (long)(fail % 40u) : -1, 0};
  GSHGEOM_Allocator alloc = {&counter, c_malloc, c_calloc, c_realloc, c_free};
  GSHGEOM_Limits limits;
  gshgeom_limits_default(&limits);
  if (small_points) {
    limits.max_points = 8u + fail % 40u;
    limits.max_elements = 16u + fail % 90u;
  }
  const uint64_t max_lines = small_lines ? 1u + fail % 12u : 0;

  const GSHGEOM_Result hc = gshgeom_hatch_check(&hatch);
  GSHGEOM_Hatch_Pieces * pieces = NULL;
  const GSHGEOM_Result r = gshgeom_hatch_pieces(region, rule, &hatch, &chain, max_lines, &alloc, &limits, &pieces);
  CHECK(r == GSHGEOM_OK || r == GSHGEOM_ERR_INVALID || r == GSHGEOM_ERR_UNSUPPORTED || r == GSHGEOM_ERR_LIMIT || r == GSHGEOM_ERR_OOM);
  if (hc != GSHGEOM_OK && r != GSHGEOM_ERR_OOM) {
    /* HX2's order for the whole call: an invalid hatch or chain, then a 3D path, chain or hatch map, then a bound of the hatch. */
    const GSHGEOM_Result want = hc == GSHGEOM_ERR_INVALID || gshgeom_map_check(&chain) != GSHGEOM_OK
                                    ? GSHGEOM_ERR_INVALID
                                    : (gshgeom_path_dimension(region) != 2 || chain.dim != 2) ? GSHGEOM_ERR_UNSUPPORTED : hc;
    CHECK(r == want);
  }
  if (r != GSHGEOM_OK) {
    CHECK(pieces == NULL);
    CHECK(counter.live == 0);
    if (r == GSHGEOM_ERR_INVALID) {
      outcomes[OUT_INVALID]++;
    } else if (r == GSHGEOM_ERR_UNSUPPORTED) {
      outcomes[OUT_UNSUPPORTED]++;
    } else if (r == GSHGEOM_ERR_LIMIT) {
      outcomes[small_lines ? OUT_LINE_LIMIT : OUT_POINT_LIMIT]++;
    } else {
      CHECK(fail_alloc && counter.failed);
      outcomes[OUT_OOM]++;
    }
    gshgeom_path_destroy(region);
    return 0;
  }
  CHECK(pieces != NULL);
  CHECK(max_lines == 0 || gshgeom_hatch_pieces_lines(pieces) <= max_lines);
  const size_t n = gshgeom_hatch_pieces_count(pieces);
  int dots = 0;
  for (size_t i = 0; i < n; i++) {
    double xy[4];
    int dot = 0;
    CHECK(gshgeom_hatch_pieces_get(pieces, i, xy, &dot) == GSHGEOM_OK);
    for (int q = 0; q < 4; q++) {
      CHECK(isfinite(xy[q]));
    }
    if (dot) {
      CHECK(xy[0] == xy[2] && xy[1] == xy[3]);
      dots++;
    }
  }
  outcomes[n ? OUT_PIECES : OUT_NO_PIECES]++;
  outcomes[OUT_DOT] += dots != 0;
  for (size_t i = 0; i < nf; i++) {
    outcomes[OUT_DASHED] += fam[i].dash_count != 0 && n != 0;
  }
  outcomes[OUT_NESTED] += hatch.nesting != 0 && n != 0;
  outcomes[OUT_TWO_FAMILIES] += nf >= 2 && n != 0;
  outcomes[OUT_DEVICE_WIDTH] += device && n != 0;
  {
    /* the same call twice gives the same pieces */
    GSHGEOM_Hatch_Pieces * again = NULL;
    if (!fail_alloc) {
      CHECK(gshgeom_hatch_pieces(region, rule, &hatch, &chain, max_lines, NULL, &limits, &again) == GSHGEOM_OK);
      CHECK(gshgeom_hatch_pieces_hash(again) == gshgeom_hatch_pieces_hash(pieces));
      CHECK(gshgeom_hatch_pieces_lines(again) == gshgeom_hatch_pieces_lines(pieces));
      gshgeom_hatch_pieces_destroy(again);
    }
  }

  /* the expansion */
  GSHGEOM_Stroke_Outline * o = NULL;
  const GSHGEOM_Result er = gshgeom_hatch_expand_pieces(pieces, &st, &chain, GSHGEOM_FLATTEN_TOLERANCE, &alloc, &limits, &o);
  CHECK(er == GSHGEOM_OK || er == GSHGEOM_ERR_INVALID || er == GSHGEOM_ERR_LIMIT || er == GSHGEOM_ERR_OOM);
  if (er != GSHGEOM_OK) {
    CHECK(o == NULL);
    if (er == GSHGEOM_ERR_INVALID) {
      CHECK(!(st.width >= 0.0 && st.width <= 1.7976931348623157e308) || st.cap >= (unsigned)GSHGEOM_STROKE_CAP_COUNT);
      outcomes[OUT_INVALID]++;
    } else if (er == GSHGEOM_ERR_LIMIT) {
      outcomes[OUT_POINT_LIMIT]++;
    } else {
      outcomes[OUT_OOM]++;
    }
  } else {
    CHECK(o != NULL);
    CHECK(gshgeom_stroke_outline_aligned(o) == NULL && gshgeom_stroke_outline_clip(o) == NULL);
    const GSHGEOM_Stroke_Status status = gshgeom_stroke_outline_status(o);
    if (status == GSHGEOM_STROKE_NOT_INVERTIBLE) {
      outcomes[OUT_NOT_INVERTIBLE]++;
    } else if (status == GSHGEOM_STROKE_EMPTY) {
      CHECK(st.width == 0.0);
      outcomes[OUT_EMPTY_WIDTH]++;
    } else if (!fail_alloc && !(dots && (is_square(&st) || st.degenerate == GSHGEOM_STROKE_DEGENERATE_OMIT)) && n < 400) {
      /* HX8: the stroke expansion of the pieces as an explicit path */
      GSHGEOM_Path * path = NULL;
      if (gshgeom_path_create(NULL, &limits, 2, &path) == GSHGEOM_OK) {
        int ok = 1;
        for (size_t i = 0; i < n && ok; i++) {
          double xy[4];
          int dot;
          (void)gshgeom_hatch_pieces_get(pieces, i, xy, &dot);
          GSHGEOM_Point a = {xy[0], xy[1], 0.0}, b = {xy[2], xy[3], 0.0};
          ok = gshgeom_path_move_to(path, a) == GSHGEOM_OK && gshgeom_path_line_to(path, b) == GSHGEOM_OK;
        }
        GSHGEOM_Stroke_Outline * ref = NULL;
        if (ok && gshgeom_stroke_expand(path, &st, &chain, NULL, &limits, &ref) == GSHGEOM_OK) {
          CHECK(gshgeom_stroke_outline_hash(ref) == gshgeom_stroke_outline_hash(o));
          outcomes[OUT_DIFFERENTIAL]++;
        }
        gshgeom_stroke_outline_destroy(ref);
        gshgeom_path_destroy(path);
      }
    }
  }
#ifdef GSHGEOM_FUZZ_PLANT
  CHECK(n == (size_t)-1); /* the replay gate's plant: a claim no pieces list can meet, which the replay must see */
#endif
  gshgeom_stroke_outline_destroy(o);
  gshgeom_hatch_pieces_destroy(pieces);
  gshgeom_path_destroy(region);
  CHECK(counter.live == 0);
  return 0;
}
