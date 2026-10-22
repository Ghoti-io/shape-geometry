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
 * The exact sign of a sum of products of doubles (internal; AD-6).
 *
 * A double is `m * 2^e` with `m` an integer below `2^53`.  A product of up to
 * three doubles is then an integer of up to 159 bits times a power of two, and
 * a sum of such terms is an integer once every term is shifted to the smallest
 * exponent among them.  The shift can be as large as the exponent range of
 * three doubles (about 6,100 bits), so the sum is formed in a two's complement
 * accumulator of 7,168 bits: wide enough that nothing overflows, and nothing
 * underflows, because no term is ever rounded.  The sign of the accumulator is
 * the sign of the sum, for every finite input.
 *
 * The method is the plain one (a long-integer sum of exact products); it
 * shares nothing with any implementation of the predicates and is checked
 * against Python's rational arithmetic in the tests.
 */

#ifndef GHOTI_IO_GSHGEOM_EXACT_H
#define GHOTI_IO_GSHGEOM_EXACT_H

#include <ghoti.io/shape-geometry/macros.h>

#include <stddef.h>
#include <stdint.h>

/** The most terms one call takes. */
#define GSHGEOM_EXACT_MAX_TERMS 32

/** One term: `(-1)^neg * f[0] * ... * f[n-1]`, with `n` from 1 to 3. */
typedef struct GSHGEOM_ExactTerm {
  int neg;      ///< Nonzero if the term is subtracted.
  int n;        ///< The number of factors, 1 to 3.
  double f[3];  ///< The factors.  Every one is finite.
} GSHGEOM_ExactTerm;

/** Append the term `(-1)^neg * f0 [* f1 [* f2]]` (`factors` of them) at `t[*n]`. */
void gshgeom_exact_put(GSHGEOM_ExactTerm * t, size_t * n, int neg, int factors,
    double f0, double f1, double f2);

/**
 * Append the six terms of `(-1)^neg * det[r1; r2; r3]`, the 3 by 3 determinant
 * of the rows `r1`, `r2`, `r3` (each three doubles: x, y, z).
 */
void gshgeom_exact_put_det3(GSHGEOM_ExactTerm * t, size_t * n, int neg,
    const double * r1, const double * r2, const double * r3);

/**
 * The sign of the sum of the terms: -1, 0 or +1, exactly.
 *
 * Every factor must be finite and `count` at most ::GSHGEOM_EXACT_MAX_TERMS
 * (a caller bug otherwise; zero is returned).
 */
int gshgeom_exact_sign(const GSHGEOM_ExactTerm * terms, size_t count);

/* --- the wide sum (sum.c) ------------------------------------------------------- */

/** The most factors one term of the wide sum takes. */
#define GSHGEOM_EXACT_SUM_FACTORS 12
/** Words of one product of GSHGEOM_EXACT_SUM_FACTORS mantissas (53 bits each, two words each). */
#define GSHGEOM_EXACT_SUM_PRODUCT_WORDS 24
/** The wide accumulator: 25,600 bits of two's complement (sum.c says why that is enough). */
#define GSHGEOM_EXACT_SUM_WORDS 800

/**
 * An exact sum of products of up to twelve doubles, built term by term.  Its
 * sign is exact for any number of terms below a few thousand.
 */
typedef struct GSHGEOM_ExactSum {
  uint32_t w[GSHGEOM_EXACT_SUM_WORDS]; ///< The sum times 2^12888, two's complement, least significant first.
  int bad;                             ///< Set by a term with a non-finite factor or a bad count (a caller bug).
} GSHGEOM_ExactSum;

/** Start an empty sum. */
void gshgeom_exact_sum_init(GSHGEOM_ExactSum * s);

/** Add `(-1)^neg * f[0] * ... * f[n-1]`, `n` from 1 to GSHGEOM_EXACT_SUM_FACTORS, every factor finite. */
void gshgeom_exact_sum_add(GSHGEOM_ExactSum * s, int neg, int n, const double * f);

/** The sign of the sum: -1, 0 or +1, exactly (0 after a bad term). */
int gshgeom_exact_sum_sign(const GSHGEOM_ExactSum * s);

#endif /* GHOTI_IO_GSHGEOM_EXACT_H */
