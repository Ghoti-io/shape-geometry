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
 * The exact sign of a long sum of products of up to twelve doubles (see
 * exact.h): the degree the stroke's miter limit (SG9) and the arcs join's
 * curvature test (SG8) need, beyond the three-factor sums of exact.c.
 *
 * Every product of `n <= 12` finite doubles is `m 2^e` with `m` an integer of at
 * most `12 * 53 = 636` bits and `e >= -1074 n >= -12888`.  The accumulator holds
 * the sum scaled by `2^12888`, so every product is added at a non-negative bit
 * offset, at most `12 * 971 + 12888 = 24540`; with the mantissa and the carries
 * of a few thousand terms the sum stays below `2^25180`, inside the 25,600-bit
 * two's complement accumulator.  Nothing is rounded: the sign is exact for every
 * finite input, subnormals and values near `DBL_MAX` included.  The method is
 * the plain long-integer sum, the same as exact.c's at a wider size; the tests
 * check it against Python's rationals.
 */

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/shape-geometry/hash.h>

#include "exact.h"

#include <stdint.h>
#include <string.h>

/** The scale's exponent: every product is added at (its exponent + SUM_BASE) bits. */
#define SUM_BASE 12888L

void gshgeom_exact_sum_init(GSHGEOM_ExactSum * s) {
  memset(s->w, 0, sizeof s->w);
  s->bad = 0;
}

/** Cut a finite double into `m * 2^e` and a sign. */
static void cut(double x, uint64_t * m, long * e, int * neg) {
  uint64_t bits = gshgeom_f64_bits(x);
  unsigned ex = (unsigned)((bits >> 52) & 0x7ffu);
  uint64_t frac = bits & ((UINT64_C(1) << 52) - 1u);
  *neg = (int)(bits >> 63);
  if (ex == 0) {
    *m = frac;
    *e = -1074;
  } else {
    *m = frac | (UINT64_C(1) << 52);
    *e = (long)ex - 1075;
  }
}

void gshgeom_exact_sum_add(GSHGEOM_ExactSum * s, int neg, int n, const double * f) {
  if (n < 1 || n > GSHGEOM_EXACT_SUM_FACTORS) {
    s->bad = 1;
    return;
  }
  uint32_t v[GSHGEOM_EXACT_SUM_PRODUCT_WORDS];
  unsigned len = 0;
  long e = 0;
  int sign = neg ? 1 : 0;
  memset(v, 0, sizeof v);
  for (int k = 0; k < n; k++) {
    uint64_t m;
    long ek;
    int nk;
    if (f[k] - f[k] != 0.0) {
      s->bad = 1; /* not finite: a caller bug */
      return;
    }
    cut(f[k], &m, &ek, &nk);
    if (m == 0) {
      return; /* a zero product adds nothing */
    }
    sign ^= nk;
    e += ek;
    if (k == 0) {
      v[0] = (uint32_t)(m & 0xffffffffu);
      v[1] = (uint32_t)(m >> 32);
      len = 2;
      continue;
    }
    uint32_t out[GSHGEOM_EXACT_SUM_PRODUCT_WORDS];
    memset(out, 0, sizeof out);
    for (unsigned j = 0; j < 2; j++) {
      uint64_t mj = (m >> (32u * j)) & 0xffffffffu;
      uint64_t carry = 0;
      for (unsigned i = 0; i < len; i++) {
        uint64_t cur = (uint64_t)v[i] * mj + out[i + j] + carry;
        out[i + j] = (uint32_t)(cur & 0xffffffffu);
        carry = cur >> 32;
      }
      out[len + j] = (uint32_t)carry;
    }
    len += 2;
    memcpy(v, out, sizeof v);
  }
  long off = e + SUM_BASE;
  unsigned long q = (unsigned long)off / 32u;
  unsigned r = (unsigned)((unsigned long)off % 32u);
  /* shifted product: len + 1 words */
  uint32_t sh[GSHGEOM_EXACT_SUM_PRODUCT_WORDS + 1];
  if (r == 0) {
    for (unsigned i = 0; i < len; i++) {
      sh[i] = v[i];
    }
    sh[len] = 0;
  } else {
    sh[0] = v[0] << r;
    for (unsigned i = 1; i < len; i++) {
      sh[i] = (v[i] << r) | (v[i - 1] >> (32u - r));
    }
    sh[len] = v[len - 1] >> (32u - r);
  }
  if (!sign) {
    uint64_t carry = 0;
    unsigned long i = 0;
    for (; i <= len; i++) {
      uint64_t cur = (uint64_t)s->w[q + i] + sh[i] + carry;
      s->w[q + i] = (uint32_t)(cur & 0xffffffffu);
      carry = cur >> 32;
    }
    for (i = q + len + 1; carry && i < GSHGEOM_EXACT_SUM_WORDS; i++) {
      uint64_t cur = (uint64_t)s->w[i] + carry;
      s->w[i] = (uint32_t)(cur & 0xffffffffu);
      carry = cur >> 32;
    }
  } else {
    uint64_t borrow = 0;
    unsigned long i = 0;
    for (; i <= len; i++) {
      uint64_t cur = (uint64_t)s->w[q + i] - sh[i] - borrow;
      s->w[q + i] = (uint32_t)(cur & 0xffffffffu);
      borrow = (cur >> 63) & 1u;
    }
    for (i = q + len + 1; borrow && i < GSHGEOM_EXACT_SUM_WORDS; i++) {
      uint64_t cur = (uint64_t)s->w[i] - borrow;
      s->w[i] = (uint32_t)(cur & 0xffffffffu);
      borrow = (cur >> 63) & 1u;
    }
  }
}

int gshgeom_exact_sum_sign(const GSHGEOM_ExactSum * s) {
  if (s->bad) {
    return 0;
  }
  if (s->w[GSHGEOM_EXACT_SUM_WORDS - 1] & 0x80000000u) {
    return -1;
  }
  for (unsigned i = 0; i < GSHGEOM_EXACT_SUM_WORDS; i++) {
    if (s->w[i]) {
      return 1;
    }
  }
  return 0;
}

double gshgeom_exact_sum_to_double(const GSHGEOM_ExactSum * s) {
  if (s->bad) {
    return 0.0;
  }
  uint32_t mag[GSHGEOM_EXACT_SUM_WORDS];
  int neg = (s->w[GSHGEOM_EXACT_SUM_WORDS - 1] & 0x80000000u) != 0;
  if (neg) {
    /* two's complement negation: invert and add one */
    uint64_t carry = 1;
    for (unsigned i = 0; i < GSHGEOM_EXACT_SUM_WORDS; i++) {
      uint64_t cur = (uint64_t)(uint32_t)~s->w[i] + carry;
      mag[i] = (uint32_t)(cur & 0xffffffffu);
      carry = cur >> 32;
    }
  } else {
    memcpy(mag, s->w, sizeof mag);
  }
  long hi = -1; /* index of the highest nonzero word */
  for (long i = GSHGEOM_EXACT_SUM_WORDS - 1; i >= 0; i--) {
    if (mag[i]) {
      hi = i;
      break;
    }
  }
  if (hi < 0) {
    return 0.0;
  }
  unsigned lead = 31u;
  while (!((mag[hi] >> lead) & 1u)) {
    lead--;
  }
  long top = hi * 32 + (long)lead; /* the position of the leading bit */
  /* the 64 bits from the leading bit down (zero-filled below bit 0), and whether anything lies below them */
  uint64_t m64 = 0;
  int sticky = 0;
  for (int k = 0; k < 64; k++) {
    long b = top - k;
    m64 <<= 1;
    if (b >= 0 && ((mag[b / 32] >> (b % 32)) & 1u)) {
      m64 |= 1u;
    }
  }
  for (long b = top - 64; b >= 0 && !sticky; b--) {
    if ((mag[b / 32] >> (b % 32)) & 1u) {
      sticky = 1;
    }
  }
  uint64_t m = m64 >> 11;
  uint64_t rest = m64 & 0x7ffu;
  if (rest > 0x400u || (rest == 0x400u && (sticky || (m & 1u)))) {
    m++;
  }
  long e2 = top - 52 - SUM_BASE; /* the value is m * 2^e2 */
  if (m == (UINT64_C(1) << 53)) {
    m >>= 1;
    e2++;
  }
  /* m < 2^53, so each of these terms and every partial sum is an exact integer.  Not `(double)(int64_t)m`: a 32-bit
   * x86 build converts a 64-bit integer with the x87 `fild`, and this library's code holds no x87 instruction. */
  double v = ((double)(int32_t)(uint32_t)(m >> 32) * 4294967296.0 + (double)(int32_t)(uint32_t)((m >> 16) & 0xffffu) * 65536.0) +
      (double)(int32_t)(uint32_t)(m & 0xffffu);
  /* scale by 2^e2 in steps that are themselves exact powers of two */
  while (e2 > 0) {
    long st = e2 > 1000 ? 1000 : e2;
    uint64_t bits = (uint64_t)(st + 1023) << 52;
    double p;
    memcpy(&p, &bits, sizeof p);
    v *= p;
    e2 -= st;
  }
  while (e2 < 0) {
    long st = e2 < -1000 ? -1000 : e2;
    uint64_t bits = (uint64_t)(st + 1023) << 52;
    double p;
    memcpy(&p, &bits, sizeof p);
    v *= p;
    e2 -= st;
  }
  return neg ? -v : v;
}
