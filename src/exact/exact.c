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
 * The exact sign of a sum of products (see exact.h).
 */

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/shape-geometry/hash.h>

#include "exact.h"

#include <stdint.h>

/** Accumulator width in 32-bit words: 7,168 bits.  The widest sum is three
 * factors at the full exponent range, 6,135 bits of shift plus 159 of mantissa
 * product and a few of carry: under 6,400. */
#define ACC_WORDS 224u

/** A product's mantissa: up to six words (three 53-bit factors fit in five). */
typedef struct Prod {
  uint32_t v[6]; ///< The integer, least significant word first.
  int neg;       ///< The sign.
  long e;        ///< The power of two.
  int zero;      ///< Nonzero if the product is zero.
} Prod;

/** Cut a finite double into `m * 2^e` and a sign. */
static void decompose(double x, uint64_t * m, int * e, int * neg) {
  uint64_t bits = gshgeom_f64_bits(x);
  unsigned ex = (unsigned)((bits >> 52) & 0x7ffu);
  uint64_t frac = bits & ((UINT64_C(1) << 52) - 1u);
  *neg = (int)(bits >> 63);
  if (ex == 0) {
    *m = frac;
    *e = -1074;
  } else {
    *m = frac | (UINT64_C(1) << 52);
    *e = (int)ex - 1075;
  }
}

/** Form the product of a term's factors exactly. */
static void make_prod(const GSHGEOM_ExactTerm * t, Prod * p) {
  uint32_t v[6] = {0, 0, 0, 0, 0, 0};
  unsigned len = 0;
  for (unsigned i = 0; i < 6; i++) {
    p->v[i] = 0;
  }
  p->neg = t->neg ? 1 : 0;
  p->e = 0;
  p->zero = 0;
  for (int k = 0; k < t->n; k++) {
    uint64_t m;
    int e;
    int neg;
    decompose(t->f[k], &m, &e, &neg);
    if (m == 0) {
      p->zero = 1;
      return;
    }
    p->neg ^= neg;
    p->e += e;
    if (k == 0) {
      v[0] = (uint32_t)(m & 0xffffffffu);
      v[1] = (uint32_t)(m >> 32);
      len = 2;
      continue;
    }
    uint32_t out[6] = {0, 0, 0, 0, 0, 0};
    for (unsigned j = 0; j < 2; j++) {
      uint64_t mj = (m >> (32u * j)) & 0xffffffffu;
      uint64_t carry = 0;
      for (unsigned i = 0; i < len; i++) {
        /* At most (2^32-1)^2 + 2 (2^32-1) = 2^64 - 1: no overflow. */
        uint64_t cur = (uint64_t)v[i] * mj + out[i + j] + carry;
        out[i + j] = (uint32_t)(cur & 0xffffffffu);
        carry = cur >> 32;
      }
      out[len + j] = (uint32_t)carry;
    }
    len += 2;
    for (unsigned i = 0; i < 6; i++) {
      v[i] = out[i];
    }
  }
  for (unsigned i = 0; i < 6; i++) {
    p->v[i] = v[i];
  }
}

/** acc += (or -=) p << off bits, in two's complement. */
static void accumulate(uint32_t * acc, const Prod * p, unsigned long off) {
  unsigned long q = off / 32u;
  unsigned r = (unsigned)(off % 32u);
  uint32_t s[7];
  if (r == 0) {
    for (unsigned i = 0; i < 6; i++) {
      s[i] = p->v[i];
    }
    s[6] = 0;
  } else {
    s[0] = p->v[0] << r;
    for (unsigned i = 1; i < 6; i++) {
      s[i] = (p->v[i] << r) | (p->v[i - 1] >> (32u - r));
    }
    s[6] = p->v[5] >> (32u - r);
  }
  if (!p->neg) {
    uint64_t carry = 0;
    for (unsigned i = 0; i < 7; i++) {
      uint64_t cur = (uint64_t)acc[q + i] + s[i] + carry;
      acc[q + i] = (uint32_t)(cur & 0xffffffffu);
      carry = cur >> 32;
    }
    for (unsigned long i = q + 7; carry && i < ACC_WORDS; i++) {
      uint64_t cur = (uint64_t)acc[i] + carry;
      acc[i] = (uint32_t)(cur & 0xffffffffu);
      carry = cur >> 32;
    }
  } else {
    uint64_t borrow = 0;
    for (unsigned i = 0; i < 7; i++) {
      uint64_t cur = (uint64_t)acc[q + i] - s[i] - borrow;
      acc[q + i] = (uint32_t)(cur & 0xffffffffu);
      borrow = (cur >> 63) & 1u;
    }
    for (unsigned long i = q + 7; borrow && i < ACC_WORDS; i++) {
      uint64_t cur = (uint64_t)acc[i] - borrow;
      acc[i] = (uint32_t)(cur & 0xffffffffu);
      borrow = (cur >> 63) & 1u;
    }
  }
}

int gshgeom_exact_sign(const GSHGEOM_ExactTerm * terms, size_t count) {
  if (count > GSHGEOM_EXACT_MAX_TERMS || (count && !terms)) {
    return 0;
  }
  Prod ps[GSHGEOM_EXACT_MAX_TERMS];
  long emin = 0;
  int have = 0;
  for (size_t i = 0; i < count; i++) {
    make_prod(&terms[i], &ps[i]);
    if (!ps[i].zero && (!have || ps[i].e < emin)) {
      emin = ps[i].e;
      have = 1;
    }
  }
  if (!have) {
    return 0;
  }
  uint32_t acc[ACC_WORDS];
  for (unsigned i = 0; i < ACC_WORDS; i++) {
    acc[i] = 0;
  }
  for (size_t i = 0; i < count; i++) {
    if (!ps[i].zero) {
      accumulate(acc, &ps[i], (unsigned long)(ps[i].e - emin));
    }
  }
  if (acc[ACC_WORDS - 1] & 0x80000000u) {
    return -1;
  }
  for (unsigned i = 0; i < ACC_WORDS; i++) {
    if (acc[i]) {
      return 1;
    }
  }
  return 0;
}

void gshgeom_exact_put(GSHGEOM_ExactTerm * t, size_t * n, int neg, int factors,
    double f0, double f1, double f2) {
  t[*n].neg = neg;
  t[*n].n = factors;
  t[*n].f[0] = f0;
  t[*n].f[1] = f1;
  t[*n].f[2] = f2;
  (*n)++;
}

void gshgeom_exact_put_det3(GSHGEOM_ExactTerm * t, size_t * n, int neg,
    const double * r1, const double * r2, const double * r3) {
  /* r1x (r2y r3z - r2z r3y) - r1y (r2x r3z - r2z r3x) + r1z (r2x r3y - r2y r3x) */
  gshgeom_exact_put(t, n, neg, 3, r1[0], r2[1], r3[2]);
  gshgeom_exact_put(t, n, !neg, 3, r1[0], r2[2], r3[1]);
  gshgeom_exact_put(t, n, !neg, 3, r1[1], r2[0], r3[2]);
  gshgeom_exact_put(t, n, neg, 3, r1[1], r2[2], r3[0]);
  gshgeom_exact_put(t, n, neg, 3, r1[2], r2[0], r3[1]);
  gshgeom_exact_put(t, n, !neg, 3, r1[2], r2[1], r3[0]);
}
