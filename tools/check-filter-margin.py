#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-only
#
# Copyright (C) 2026 Corey Pennycuff
#
# This file is part of Ghoti.io Shape Geometry.
#
# Ghoti.io Shape Geometry is free software: you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License version 3
# as published by the Free Software Foundation.
"""The sign filter of src/lattice/lattice.c at its boundary, against exact rational arithmetic.

con_sign decides `S = c i + e j + g` from `s = (fl(c i) + fl(e j)) + g` when `|s| > bound`, `bound = 2^-48 (|ci| + |ej| + |g|)`, and
leaves the rest to the exact sum.  c, e and g are doubles that were each rounded once from an exact sum, so the S that must be decided
is `c i (1 + dc) + e j (1 + de) + g (1 + dg)` with `|d| <= 2^-53`.  This checks, for a deterministic table of values at the edges of
that margin (mantissas 1, 1 + 2^-52, 1.5, 2 - 2^-52; exponents within a few of each other; i and j from 1 to 2^53 - 1; g the rounding
of -(c i + e j) and its neighbours from -64 to 64 ulps, which is across the bound (about 16 to 32 ulps of g), so that S is as near zero as the bound lets a decision be made; the three variants the plants name,
`ignores g`, `ignores j` and the real bound, and a control with a bound of 2^-56 of the sum that must be seen to decide wrongly; and the roundings of c, e and g taken as the worst case in both directions), that no
decision a variant makes disagrees with the exact sign.  The bound's derivation, which the table exercises and does not replace:

  |s - S| <= 2^-53 (3 (|ci| + |ej|) + |g| + |s|)
  real bound 32 * 2^-53 (|ci| + |ej| + |g|);   without |g|: if |g| <= 29 (|ci| + |ej|) the bound is at least the error, else |s| >= 0.96 |g|;
  without |ej|: if |ej| <= 4 (|ci| + |g|) the error is at most 15/32 of the bound, else |s| >= 0.75 |ej| against an error of 3.1 * 2^-53 |ej|.

Run: python3 tools/check-filter-margin.py   (exit 0 when no variant ever decides wrongly; prints the count of decisions it checked)
"""
import itertools
import sys
from fractions import Fraction as F
import math

B = 2.0 ** -48
mants = [1.0, 1.0 + 2.0 ** -52, 1.5, 2.0 - 2.0 ** -52]
exps = range(-2, 3)
ij = [1, 3, (1 << 20) + 1, (1 << 53) - 1]


def decide(variant, c, e, g, di, dj):
    ci = c * di
    ej = e * dj
    s = (ci + ej) + g
    a = abs(ci) + (0.0 if variant == "ignores_j" else abs(ej)) + (0.0 if variant == "ignores_g" else abs(g))
    bound = a * (B if variant != "control" else 2.0 ** -56)
    if s > bound:
        return 1
    if s < -bound:
        return -1
    return 0


def sign(x):
    return (x > 0) - (x < 0)


decided = wrong = control_wrong = 0
for mc, me, ec, ee, si, sj in itertools.product(mants, mants, exps, exps, (1, -1), (1, -1)):
    c = math.ldexp(mc, ec)
    e = math.ldexp(me, ee)
    for i in ij:
        for j in ij:
            di, dj = float(si * i), float(sj * j)
            t = F(c) * int(di) + F(e) * int(dj)
            g0 = -float(t)
            ulp = math.ulp(g0) if g0 != 0 else 5e-324
            for k in (-64, -40, -33, -32, -31, -24, -17, -16, -15, -8, -2, 0, 2, 8, 15, 16, 17, 24, 31, 32, 33, 40, 64):
                g = g0 + k * ulp
                # the exact sum is decided for the doubles c, e, g, and for each of their roundings at the worst end of 2^-53
                for dc, de, dg in ((0, 0, 0), (1, 1, 1), (-1, -1, -1), (1, -1, 1), (-1, 1, -1), (1, 1, -1), (-1, -1, 1)):
                    S = F(c) * (1 + dc * F(1, 1 << 53)) * int(di) + F(e) * (1 + de * F(1, 1 << 53)) * int(dj) + F(g) * (1 + dg * F(1, 1 << 53))
                    for variant in ("real", "ignores_g", "ignores_j", "control"):
                        d = decide(variant, c, e, g, di, dj)
                        if d != 0 and d != sign(S):
                            if variant == "control":
                                control_wrong += 1
                            else:
                                wrong += 1
                                if wrong <= 5:
                                    print("WRONG", variant, c, e, g, di, dj, d, sign(S))
                        if d != 0 and variant != "control":
                            decided += 1
# the third plant: the filter trusts a coefficient that underflowed to zero (c = 0.0 while the exact c is not zero).  The exact c is below
# 2^-1075, so c i is below 2^-1022 in magnitude; any decision has |s| > bound >= 2^-48 * 2^-900 (e and g are in the filter's range or
# zero), far above it.  The table: e and g at the edges of the range, c_exact below the underflow, i and j to 2^53.
under = wrong_under = 0
edge = [0.0, 2.0 ** -900, -(2.0 ** -900), 2.0 ** -900 * (1 + 2.0 ** -52), 2.0 ** 900, -(2.0 ** 900), 2.0 ** -899, 1.0]
for e_, g_ in itertools.product(edge, edge):
    for ce in (F(1, 1 << 1076), F(-1, 1 << 1076), F(1, 1 << 1100), F(-3, 1 << 1075)):
        for i in ij:
            for j in ij + [0]:
                for si, sj in itertools.product((1, -1), (1, -1)):
                    di, dj = float(si * i), float(sj * j)
                    for k in (0, 1, -1, 16, -16, 40, -40):
                        g = g_ + k * (math.ulp(g_) if g_ else 5e-324)
                        if g != 0.0 and not (2.0 ** -900 <= abs(g) <= 2.0 ** 900):
                            continue  # the filter is off for a coefficient outside its range (in_filter_range)
                        S = ce * int(di) + F(e_) * int(dj) + F(g)
                        d = decide("real", 0.0, e_, g, di, dj)
                        if d != 0:
                            under += 1
                            if d != sign(S):
                                wrong_under += 1
                                print("WRONG underflow", e_, g, di, dj, d, sign(S))
print("underflowed c: decisions checked:", under, " wrong:", wrong_under)
wrong += wrong_under
print("decisions checked:", decided, " wrong:", wrong, " (control, a bound 2^-56 of the sum: wrong", control_wrong, "times, as it must)")
sys.exit(1 if wrong or not control_wrong else 0)
