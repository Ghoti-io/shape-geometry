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
 * @file hash.h
 *
 * The hash spelling: the one definition, for every intermediate of the shape
 * family (AD-21, CAP-13).
 *
 * A hash of an intermediate is FNV-1a 64 over a byte stream, and the byte
 * stream is a sequence of fields written in a platform-independent spelling.
 * This header is where that spelling is stated; `shape-scene`, `shape-cpu` and
 * every sink hash through ::GSHGEOM_Hash and do not restate it.
 *
 * ### The spelling
 *
 * - **The function.**  FNV-1a, 64 bits.  The state starts at the offset basis
 *   ::GSHGEOM_FNV_OFFSET_BASIS (`0xcbf29ce484222325`).  For each byte `b` in
 *   order, `state = (state XOR b) * ::GSHGEOM_FNV_PRIME (0x100000001b3)`,
 *   modulo 2^64.  The hash is the final state.  An empty stream hashes to the
 *   offset basis, not to zero.
 * - **Integers** are written as their two's-complement value in a fixed number
 *   of bytes, least significant byte first: 1 byte for `u8`, 2 for `u16`, 4 for
 *   `u32` and `i32`, 8 for `u64` and `i64`.  Each byte is produced by shifting
 *   the value, never by reading the bytes of a variable, so the stream is the
 *   same on a big-endian host.
 * - **A `double`** is written as the 64-bit pattern of its IEEE 754 binary64
 *   representation, as a `u64`.  Nothing is canonicalised: `-0.0` and `+0.0`
 *   differ, and so would two NaNs with different payloads (the shape
 *   libraries refuse non-finite values, so the case does not arise in their
 *   intermediates).  A `float` is converted to `double` first (exactly).
 * - **A string** is its length in bytes as a `u32`, then its bytes, with no
 *   terminator.  A byte sequence of known length is written with
 *   ::gshgeom_hash_bytes and is not prefixed.
 * - **Booleans and enumerations** are written as `u8` or `u32` as the owner of
 *   the intermediate states, with the enumeration's pinned value.
 * - **Each intermediate states its own fields in order** in the documentation
 *   of its `_hash`, and starts with a `u32` format version where it can change.
 *
 * `font`'s `gfnt_coverage_hash` is this spelling (`u32` fields, then bytes).
 * The known answers that pin it are in `tests/unit/test_hash.cpp`, produced by
 * `tools/gen-hash-vectors.py`, which shares no code with this library.
 */

#ifndef GHOTI_IO_GSHGEOM_HASH_H
#define GHOTI_IO_GSHGEOM_HASH_H

#include <ghoti.io/shape-geometry/macros.h>

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/** FNV-1a 64's offset basis. */
#define GSHGEOM_FNV_OFFSET_BASIS UINT64_C(0xcbf29ce484222325)

/** FNV-1a 64's prime. */
#define GSHGEOM_FNV_PRIME UINT64_C(0x00000100000001b3)

/**
 * @brief A hash in progress.
 *
 * A value: copy it to fork a hash.  Fill it with ::gshgeom_hash_init.
 */
typedef struct GSHGEOM_Hash {
  uint64_t state; ///< The FNV-1a state.  Read it with ::gshgeom_hash_value.
} GSHGEOM_Hash;

/**
 * @brief Start a hash: the state becomes ::GSHGEOM_FNV_OFFSET_BASIS.
 *
 * @param hash The hash to start; NULL is ignored (as are NULL arguments to the
 *   functions below).
 */
GSHGEOM_API void gshgeom_hash_init(GSHGEOM_Hash * hash);

/** @brief Mix one byte. */
GSHGEOM_API void gshgeom_hash_u8(GSHGEOM_Hash * hash, uint8_t value);

/** @brief Mix a 16-bit unsigned integer, least significant byte first. */
GSHGEOM_API void gshgeom_hash_u16(GSHGEOM_Hash * hash, uint16_t value);

/** @brief Mix a 32-bit unsigned integer, least significant byte first. */
GSHGEOM_API void gshgeom_hash_u32(GSHGEOM_Hash * hash, uint32_t value);

/** @brief Mix a 64-bit unsigned integer, least significant byte first. */
GSHGEOM_API void gshgeom_hash_u64(GSHGEOM_Hash * hash, uint64_t value);

/** @brief Mix a 32-bit signed integer as its two's complement `u32`. */
GSHGEOM_API void gshgeom_hash_i32(GSHGEOM_Hash * hash, int32_t value);

/** @brief Mix a 64-bit signed integer as its two's complement `u64`. */
GSHGEOM_API void gshgeom_hash_i64(GSHGEOM_Hash * hash, int64_t value);

/**
 * @brief Mix a `double` as the `u64` of its IEEE 754 binary64 bit pattern.
 *
 * Nothing is canonicalised (see the spelling above).
 */
GSHGEOM_API void gshgeom_hash_f64(GSHGEOM_Hash * hash, double value);

/**
 * @brief Mix @p size bytes, in order, with no length prefix.
 *
 * @param bytes The bytes; may be NULL only when @p size is 0.
 * @param size How many.
 */
GSHGEOM_API void gshgeom_hash_bytes(
    GSHGEOM_Hash * hash, const void * bytes, size_t size);

/**
 * @brief Mix a string: its length in bytes as a `u32`, then the bytes.
 *
 * @param string A NUL-terminated string.  A string longer than 2^32-1 bytes
 *   mixes its length modulo 2^32 and all its bytes.
 */
GSHGEOM_API void gshgeom_hash_string(GSHGEOM_Hash * hash, const char * string);

/**
 * @brief The hash of everything mixed so far.
 *
 * @return The state; the offset basis for a hash with nothing mixed in, and
 *   for a NULL @p hash.
 */
GSHGEOM_API uint64_t gshgeom_hash_value(const GSHGEOM_Hash * hash);

/**
 * @brief The 64-bit pattern of a `double`, as ::gshgeom_hash_f64 mixes it.
 *
 * What a dump prints, so that a dump and a hash agree on what a number is.
 */
GSHGEOM_API uint64_t gshgeom_f64_bits(double value);

#ifdef __cplusplus
}
#endif

#endif /* GHOTI_IO_GSHGEOM_HASH_H */
