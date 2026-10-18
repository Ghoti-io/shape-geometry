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
 * The hash spelling (see hash.h): FNV-1a 64 over little-endian fields.
 */

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/shape-geometry/hash.h>

#include <string.h>

/* Every multi-byte field is cut into bytes by shifting, least significant
 * first, so that the stream does not depend on the host's byte order. */
static void mix_byte(GSHGEOM_Hash * hash, uint8_t byte) {
  hash->state = (hash->state ^ byte) * GSHGEOM_FNV_PRIME;
}

static void mix_le(GSHGEOM_Hash * hash, uint64_t value, int bytes) {
  for (int i = 0; i < bytes; i++) {
    mix_byte(hash, (uint8_t)((value >> (8 * i)) & 0xffu));
  }
}

void gshgeom_hash_init(GSHGEOM_Hash * hash) {
  if (hash) {
    hash->state = GSHGEOM_FNV_OFFSET_BASIS;
  }
}

void gshgeom_hash_u8(GSHGEOM_Hash * hash, uint8_t value) {
  if (hash) {
    mix_byte(hash, value);
  }
}

void gshgeom_hash_u16(GSHGEOM_Hash * hash, uint16_t value) {
  if (hash) {
    mix_le(hash, value, 2);
  }
}

void gshgeom_hash_u32(GSHGEOM_Hash * hash, uint32_t value) {
  if (hash) {
    mix_le(hash, value, 4);
  }
}

void gshgeom_hash_u64(GSHGEOM_Hash * hash, uint64_t value) {
  if (hash) {
    mix_le(hash, value, 8);
  }
}

void gshgeom_hash_i32(GSHGEOM_Hash * hash, int32_t value) {
  gshgeom_hash_u32(hash, (uint32_t)value);
}

void gshgeom_hash_i64(GSHGEOM_Hash * hash, int64_t value) {
  gshgeom_hash_u64(hash, (uint64_t)value);
}

uint64_t gshgeom_f64_bits(double value) {
  uint64_t bits;
  /* Same-width reinterpretation through memcpy: no aliasing cast. */
  memcpy(&bits, &value, sizeof bits);
  return bits;
}

void gshgeom_hash_f64(GSHGEOM_Hash * hash, double value) {
  gshgeom_hash_u64(hash, gshgeom_f64_bits(value));
}

void gshgeom_hash_bytes(GSHGEOM_Hash * hash, const void * bytes, size_t size) {
  if (!hash || !bytes) {
    return;
  }
  const uint8_t * p = (const uint8_t *)bytes;
  for (size_t i = 0; i < size; i++) {
    mix_byte(hash, p[i]);
  }
}

void gshgeom_hash_string(GSHGEOM_Hash * hash, const char * string) {
  if (!hash || !string) {
    return;
  }
  size_t length = strlen(string);
  gshgeom_hash_u32(hash, (uint32_t)(length & 0xffffffffu));
  gshgeom_hash_bytes(hash, string, length);
}

uint64_t gshgeom_hash_value(const GSHGEOM_Hash * hash) {
  return hash ? hash->state : GSHGEOM_FNV_OFFSET_BASIS;
}
