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
 * The curve outline's readback, hash and dump (stroke_curves.h, SG32 and SG41; format 1).
 */

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/shape-geometry/hash.h>

#include "../stroke/stroke_internal.h"

/** The format version of the hash and the dump. */
#define CURVES_FORMAT_VERSION 1u
/** The kind field: "STRC" as little-endian bytes. */
#define CURVES_KIND 0x43525453u

void gshgeom_stroke_curves_destroy(GSHGEOM_Stroke_Curves * o) {
  if (!o) {
    return;
  }
  gshgeom_path_destroy(o->body);
  gshgeom_path_destroy(o->aligned);
  gshgeom_path_destroy(o->clip);
  GSHGEOM_Allocator alloc = o->alloc;
  gcu_allocator_free(&alloc, o);
}

GSHGEOM_Stroke_Status gshgeom_stroke_curves_status(const GSHGEOM_Stroke_Curves * o) {
  return o ? o->status : GSHGEOM_STROKE_STATUS_COUNT;
}

GSHGEOM_Stroke_Curves_Form gshgeom_stroke_curves_form(const GSHGEOM_Stroke_Curves * o) {
  return o ? o->form : GSHGEOM_STROKE_CURVES_FORM_COUNT;
}

const GSHGEOM_Path * gshgeom_stroke_curves_body(const GSHGEOM_Stroke_Curves * o) {
  return o ? o->body : NULL;
}

const GSHGEOM_Path * gshgeom_stroke_curves_aligned(const GSHGEOM_Stroke_Curves * o) {
  return o ? o->aligned : NULL;
}

const GSHGEOM_Path * gshgeom_stroke_curves_clip(const GSHGEOM_Stroke_Curves * o) {
  return o ? o->clip : NULL;
}

GSHGEOM_Fill_Rule gshgeom_stroke_curves_clip_rule(const GSHGEOM_Stroke_Curves * o) {
  return o ? o->rule : GSHGEOM_FILL_NONZERO;
}

int gshgeom_stroke_curves_clip_complement(const GSHGEOM_Stroke_Curves * o) {
  return o ? o->complement : 0;
}

int gshgeom_stroke_curves_revision(void) {
  return GSHGEOM_STROKE_CURVES_REVISION;
}

uint64_t gshgeom_stroke_curves_hash(const GSHGEOM_Stroke_Curves * o) {
  GSHGEOM_Hash h;
  gshgeom_hash_init(&h);
  if (!o) {
    return gshgeom_hash_value(&h);
  }
  gshgeom_hash_u32(&h, CURVES_FORMAT_VERSION);
  gshgeom_hash_u32(&h, CURVES_KIND);
  gshgeom_hash_u8(&h, (uint8_t)o->status);
  gshgeom_hash_u8(&h, (uint8_t)o->alignment);
  gshgeom_hash_u8(&h, (uint8_t)o->rule);
  gshgeom_hash_u8(&h, (uint8_t)(o->complement ? 1 : 0));
  /* an absent part hashes as the empty stream: the offset basis */
  gshgeom_hash_u64(&h, gshgeom_path_hash(o->body));
  gshgeom_hash_u64(&h, gshgeom_path_hash(o->aligned));
  gshgeom_hash_u64(&h, gshgeom_path_hash(o->clip));
  return gshgeom_hash_value(&h);
}

static GSHGEOM_Result dump_part(const char * name, const GSHGEOM_Path * p, FILE * out) {
  if (!p) {
    return fprintf(out, "part %s -\n", name) < 0 ? GSHGEOM_ERR_IO : GSHGEOM_OK;
  }
  if (fprintf(out, "part %s\n", name) < 0) {
    return GSHGEOM_ERR_IO;
  }
  return gshgeom_path_dump(p, out);
}

GSHGEOM_Result gshgeom_stroke_curves_dump(const GSHGEOM_Stroke_Curves * o, FILE * out) {
  if (!o || !out) {
    return GSHGEOM_ERR_INVALID;
  }
  if (fprintf(out, "gshgeom.stroke-curves %u status=%u form=%u align=%u rule=%u complement=%d hash=%016llx\n",
          CURVES_FORMAT_VERSION, (unsigned)o->status, (unsigned)o->form, (unsigned)o->alignment, (unsigned)o->rule,
          o->complement ? 1 : 0, (unsigned long long)gshgeom_stroke_curves_hash(o)) < 0) {
    return GSHGEOM_ERR_IO;
  }
  GSHGEOM_Result r = dump_part("body", o->body, out);
  if (r == GSHGEOM_OK) {
    r = dump_part("aligned", o->aligned, out);
  }
  if (r == GSHGEOM_OK) {
    r = dump_part("clip", o->clip, out);
  }
  if (r != GSHGEOM_OK) {
    return r;
  }
  return (fflush(out) == 0 && !ferror(out)) ? GSHGEOM_OK : GSHGEOM_ERR_IO;
}
