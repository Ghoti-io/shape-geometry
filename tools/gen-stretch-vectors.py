#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-only
#
# Copyright (C) 2026 Corey Pennycuff
#
# This file is part of Ghoti.io Shape Geometry.
#
# Ghoti.io Shape Geometry is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3 as
# published by the Free Software Foundation.
"""Print the known answers of gshgeom_map_stretches_every_direction.

The predicate is true exactly when the smallest singular value of a map's 2x2 linear part is greater than 1.  With the
columns a = (m00, m10) and b = (m01, m11) and the exact reals A = a.a, B = a.b and C = b.b, that is A > 1, C > 1 and
(A - 1)(C - 1) > B^2.  This script decides every case that way, in exact rationals (`fractions.Fraction` of each
double's own value), through the eigenvalues' defining inequalities and never through a square root or a determinant, so
that a test comparing it with the library compares two readings of the definition.  It writes tests/stretch_vectors.h,
which is committed:

    python3 tools/gen-stretch-vectors.py > tests/stretch_vectors.h
    python3 tools/gen-stretch-vectors.py --check      # fails if the file differs

The cases are the ones that tell the predicate from its plausible wrong forms: a determinant (an area ratio), the larger
singular value, a comparison with >= and a floating square root; each is a map the wrong form decides differently.
"""
from fractions import Fraction as F
import os
import struct
import sys

M64 = (1 << 64) - 1


def ulp_up(x, n=1):
    bits = struct.unpack("<Q", struct.pack("<d", x))[0]
    return struct.unpack("<d", struct.pack("<Q", bits + n))[0]


def ulp_down(x, n=1):
    bits = struct.unpack("<Q", struct.pack("<d", x))[0]
    return struct.unpack("<d", struct.pack("<Q", bits - n))[0]


def stretches(m00, m01, m10, m11):
    p, q, r, s = F(m00), F(m10), F(m01), F(m11)
    a = p * p + q * q
    b = p * r + q * s
    c = r * r + s * s
    return a > 1 and c > 1 and (a - 1) * (c - 1) > b * b


def case(name, m00, m01, m10, m11):
    return (name, (m00, m01, m10, m11), int(stretches(m00, m01, m10, m11)))


def hand():
    import math
    out = []
    one = 1.0
    out.append(case("identity", 1.0, 0.0, 0.0, 1.0))
    out.append(case("scale 7", 7.0, 0.0, 0.0, 7.0))
    out.append(case("scale 0.3", 0.3, 0.0, 0.0, 0.3))
    out.append(case("scale 7 by 0.3: area up, one direction shrunk", 7.0, 0.0, 0.0, 0.3))
    out.append(case("scale 0.3 by 7", 0.3, 0.0, 0.0, 7.0))
    out.append(case("scale 1 by 2: the smallest singular value is exactly 1", 1.0, 0.0, 0.0, 2.0))
    out.append(case("scale 2 by 1", 2.0, 0.0, 0.0, 1.0))
    out.append(case("scale 2: every direction", 2.0, 0.0, 0.0, 2.0))
    out.append(case("determinant 1, singular values 2 and 0.5", 2.0, 0.0, 0.0, 0.5))
    out.append(case("determinant 4, singular values 4 and 1", 4.0, 0.0, 0.0, 1.0))
    out.append(case("determinant 4.0000001 and one value below 1", 4.0000001, 0.0, 0.0, 0.9999999))
    out.append(case("one ulp above unit scale", ulp_up(one), 0.0, 0.0, ulp_up(one)))
    out.append(case("one ulp below unit scale", ulp_down(one), 0.0, 0.0, ulp_down(one)))
    out.append(case("one axis one ulp above 1, the other 2", ulp_up(one), 0.0, 0.0, 2.0))
    out.append(case("one axis exactly 1, the other one ulp above", 1.0, 0.0, 0.0, ulp_up(one)))
    out.append(case("a mirror at scale 2", -2.0, 0.0, 0.0, 2.0))
    out.append(case("a mirror at scale 1", -1.0, 0.0, 0.0, 1.0))
    out.append(case("a quarter turn at scale 1", 0.0, -1.0, 1.0, 0.0))
    out.append(case("a quarter turn at scale 2", 0.0, -2.0, 2.0, 0.0))
    out.append(case("a quarter turn at scale 0.9", 0.0, -0.9, 0.9, 0.0))
    # 30 degrees at scale 2 and at 0.9 (the doubles nearest cos and sin, so the map is the product to rounding only)
    c30, s30 = math.cos(math.pi / 6), math.sin(math.pi / 6)
    out.append(case("a 30 degree turn at scale 2", 2 * c30, -2 * s30, 2 * s30, 2 * c30))
    out.append(case("a 30 degree turn at scale 0.9", 0.9 * c30, -0.9 * s30, 0.9 * s30, 0.9 * c30))
    out.append(case("a 30 degree turn at scale 1", c30, -s30, s30, c30))
    out.append(case("a shear: singular values 1.618 and 0.618", 1.0, 1.0, 0.0, 1.0))
    out.append(case("a shear with both values above 1", 2.0, 1.0, 0.0, 2.0))
    out.append(case("a shear whose smaller value is just under 1", 1.5, 1.0, 0.0, 1.5))
    out.append(case("a shear whose smaller value is just over 1", 1.6, 1.0, 0.0, 1.6))
    out.append(case("a skew with the columns not orthogonal", 3.0, 2.0, 0.0, 3.0))
    out.append(case("singular, rank one, large", 2.0, 4.0, 1.0, 2.0))
    out.append(case("singular, a zero column", 5.0, 0.0, 5.0, 0.0))
    out.append(case("all zero", 0.0, 0.0, 0.0, 0.0))
    out.append(case("a huge scale", 2.0 ** 60, 0.0, 0.0, 2.0 ** 60))
    out.append(case("a tiny scale", 2.0 ** -60, 0.0, 0.0, 2.0 ** -60))
    out.append(case("a huge scale by a tiny one", 2.0 ** 60, 0.0, 0.0, 2.0 ** -60))
    out.append(case("negative zero entries at scale 2", 2.0, -0.0, -0.0, 2.0))
    # Off the axes, a smallest singular value of exactly 1 with A and C both above 1: the columns (1, 1/2) and (1, -1) have
    # A = 5/4, C = 2, B = 1/2, so (A - 1)(C - 1) = 1/4 = B^2.  Every entry is exact in binary.
    out.append(case("an exact singular value of 1 off the axes", 1.0, 1.0, 0.5, -1.0))
    out.append(case("the same with one entry one ulp larger", ulp_up(1.0), 1.0, 0.5, -1.0))
    out.append(case("the same with one entry one ulp smaller", ulp_down(1.0), 1.0, 0.5, -1.0))
    out.append(case("the same with the other column one ulp larger", 1.0, ulp_up(1.0), 0.5, -1.0))
    return out


def random_cases(seed, count):
    s = seed
    out = []

    def nxt():
        nonlocal s
        s ^= s >> 12
        s ^= (s << 25) & M64
        s ^= s >> 27
        return (s * 2685821657736338717) & M64

    def val():
        k = nxt() % 6
        if k == 0:
            return float((nxt() % 17) - 8) / 4.0
        if k == 1:
            return float((nxt() % 4001) - 2000) / 1024.0
        if k == 2:
            return ulp_up(1.0, (nxt() % 5)) * (1 if nxt() & 1 else -1)
        if k == 3:
            return float((nxt() % 41) - 20) / 16.0
        if k == 4:
            return (float(nxt() % 100001) / 100000.0) * (1 if nxt() & 1 else -1)
        return float((nxt() % 7) - 3)

    for i in range(count):
        m = (val(), val(), val(), val())
        out.append(case("random %d" % i, *m))
    return out


def render():
    cases = hand() + random_cases(0x9E3779B97F4A7C15, 400)
    o = []
    o.append("// Generated by tools/gen-stretch-vectors.py: gshgeom_map_stretches_every_direction decided in exact rationals.\n")
    o.append("// Do not edit; run the generator.  Independent of the library's code.\n")
    o.append("#ifndef GHOTI_IO_GSHGEOM_TEST_STRETCH_VECTORS_H\n#define GHOTI_IO_GSHGEOM_TEST_STRETCH_VECTORS_H\n\nnamespace stretchvec {\n\n")
    o.append("struct Case {\n  const char * name;\n  double m00, m01, m10, m11;  // the linear part: columns (m00, m10) and (m01, m11)\n  int expected;\n};\n\n")
    o.append("static const Case kCases[] = {\n")
    for name, (a, b, c, d), exp in cases:
        o.append('    {"%s", %s, %s, %s, %s, %d},\n' % (name, a.hex(), b.hex(), c.hex(), d.hex(), exp))
    o.append("};\n\n")
    o.append("static const unsigned kCaseCount = sizeof kCases / sizeof kCases[0];\n\n")
    n_true = sum(1 for _, _, e in cases if e)
    o.append("static const unsigned kStretching = %d;  // how many of them stretch every direction\n\n" % n_true)
    o.append("}  // namespace stretchvec\n\n#endif\n")
    return "".join(o)


if __name__ == "__main__":
    text = render()
    if len(sys.argv) > 1 and sys.argv[1] == "--check":
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tests", "stretch_vectors.h")
        with open(path) as f:
            if f.read() != text:
                sys.exit("gen-stretch-vectors: tests/stretch_vectors.h is not what the generator writes")
        print("gen-stretch-vectors: tests/stretch_vectors.h is current")
    else:
        sys.stdout.write(text)
