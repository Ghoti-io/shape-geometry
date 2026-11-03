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
"""Print tests/lattice_vectors.h: the known answers of lattice.h (LT1 to LT3), in exact rationals.

    python3 tools/gen-lattice-vectors.py > tests/lattice_vectors.h

Shares no code with the library.  The copies of a use are found by the definition: copy (i, j) is a copy when the convex hull
of the four exact points ``Q_k + tau + i A + j B`` (Fractions built from the doubles) meets the closed rectangle.  That is
decided by geometry and not by the library's half-plane method: a point of the four in the rectangle, or one of the six
segments between them meeting a side of the rectangle, or a corner of the rectangle in a triangle of three of them.  Every
index of a window that holds every copy that can meet (computed with a margin from the inverse of the step matrix) is
tested, with a float bounding-box test to skip the clear misses.  ``Q``, ``A`` and ``B`` are formed in doubles by the
header's formulas (Python floats are IEEE doubles and do not fuse).  The translation vectors restate
``((tau + i a) + j b)``; the period vectors find the Hermite normal form by search (``s`` the gcd of the y components,
``p = |D| / s``, ``r`` the first residue with ``(r, s)`` in the lattice).
"""
import math
import random
import sys
from fractions import Fraction as F


def hexf(v):
    return float(v).hex()


def quad_and_steps(m, t, cell, steps):
    x, y, w, h = cell
    xs, ys = steps
    # PT5 / LT2 cell of copy (0, 0)
    left, top = x + 0.0 * xs, y + 0.0 * ys
    right = x + 1.0 * xs if xs == w else left + w
    bottom = y + 1.0 * ys if ys == h else top + h
    px = [left, right, right, left]
    py = [top, top, bottom, bottom]
    Q = [(m[0] * px[k] + m[1] * py[k], m[2] * px[k] + m[3] * py[k]) for k in range(4)]
    A = (m[0] * xs, m[2] * xs)
    B = (m[1] * ys, m[3] * ys)
    return Q, t, A, B


def orient(a, b, c):
    v = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    return (v > 0) - (v < 0)


def between(v, a, b):
    return min(a, b) <= v <= max(a, b)


def on_seg(a, b, p):
    return orient(a, b, p) == 0 and between(p[0], a[0], b[0]) and between(p[1], a[1], b[1])


def segs_meet(a, b, c, d):
    o1, o2, o3, o4 = orient(a, b, c), orient(a, b, d), orient(c, d, a), orient(c, d, b)
    if o1 * o2 < 0 and o3 * o4 < 0:
        return True
    return on_seg(a, b, c) or on_seg(a, b, d) or on_seg(c, d, a) or on_seg(c, d, b)


def hull_meets(P, rect):
    x0, y0, x1, y1 = rect
    R = [(F(x0), F(y0)), (F(x1), F(y0)), (F(x1), F(y1)), (F(x0), F(y1))]
    for p in P:
        if x0 <= p[0] <= x1 and y0 <= p[1] <= y1:
            return True
    for a in range(4):
        for b in range(a + 1, 4):
            for s in range(4):
                if segs_meet(P[a], P[b], R[s], R[(s + 1) % 4]):
                    return True
    for c in R:
        for a in range(4):
            for b in range(a + 1, 4):
                for d in range(b + 1, 4):
                    s = orient(P[a], P[b], P[d])
                    if s == 0:
                        continue
                    if orient(P[a], P[b], c) * s >= 0 and orient(P[b], P[d], c) * s >= 0 and orient(P[d], P[a], c) * s >= 0:
                        return True
    return False


def copy_meets(Q, tau, A, B, i, j, rect):
    # float rejection with a generous margin (the exact test decides the rest)
    fx = [q[0] + tau[0] + i * A[0] + j * B[0] for q in Q]
    fy = [q[1] + tau[1] + i * A[1] + j * B[1] for q in Q]
    mag = max(1.0, max(abs(v) for v in fx + fy))
    eps = mag * 1e-9
    if max(fx) < rect[0] - eps or min(fx) > rect[2] + eps or max(fy) < rect[1] - eps or min(fy) > rect[3] + eps:
        return False
    P = [(F(q[0]) + F(tau[0]) + i * F(A[0]) + j * F(B[0]), F(q[1]) + F(tau[1]) + i * F(A[1]) + j * F(B[1])) for q in Q]
    return hull_meets(P, rect)


def rows_of(m, t, cell, steps, repeat, rect):
    Q, tau, A, B = quad_and_steps(m, t, cell, steps)
    det = A[0] * B[1] - A[1] * B[0]
    # window: indices of the points of K, from the inverse step matrix, with a margin
    pts = [(rect[0], rect[1]), (rect[2], rect[1]), (rect[2], rect[3]), (rect[0], rect[3])]
    lo_i = lo_j = 10 ** 18
    hi_i = hi_j = -10 ** 18
    for px, py in pts:
        for q in Q:
            dx, dy = px - q[0] - tau[0], py - q[1] - tau[1]
            ii = (dx * B[1] - dy * B[0]) / det
            jj = (A[0] * dy - A[1] * dx) / det
            lo_i, hi_i = min(lo_i, ii), max(hi_i, ii)
            lo_j, hi_j = min(lo_j, jj), max(hi_j, jj)
    ir = range(math.floor(lo_i) - 3, math.ceil(hi_i) + 4) if repeat[0] else range(0, 1)
    jr = range(math.floor(lo_j) - 3, math.ceil(hi_j) + 4) if repeat[1] else range(0, 1)
    rows = []
    for j in jr:
        hits = [i for i in ir if copy_meets(Q, tau, A, B, i, j, rect)]
        if hits:
            assert hits == list(range(hits[0], hits[-1] + 1)), "a row with a hole"
            rows.append((j, hits[0], hits[-1]))
    return rows


def rot(deg, sx=1.0, sy=1.0):
    r = math.radians(deg)
    c, s = math.cos(r), math.sin(r)
    return [c * sx, -s * sy, s * sx, c * sy]


def each_cases():
    rnd = random.Random(20261102)
    cases = []

    def add(name, m, t, cell, steps, repeat, rect):
        cases.append((name, m, t, cell, steps, repeat, rect))

    add("abut-identity", [1, 0, 0, 1], [0, 0], [0, 0, 4, 4], [4, 4], [1, 1], [3, 3, 21, 17])
    add("gap-identity", [1, 0, 0, 1], [0.5, 0.25], [0, 0, 3, 3], [5, 7], [1, 1], [0, 0, 40, 40])
    add("overlap-identity", [1, 0, 0, 1], [0, 0], [0, 0, 6, 6], [4, 4], [1, 1], [10, 10, 30, 30])
    add("rot30", rot(30), [3.3, 4.4], [0, 0, 5, 3], [5, 3], [1, 1], [0, 0, 60, 45])
    add("rot45-gap", rot(45, 1.5, 0.75), [-7.1, 2.2], [0.5, 0.5, 2.5, 4.5], [3.7, 6.1], [1, 1], [-20, -10, 25, 33])
    add("skew", [1, 0.6, 0.2, 1], [0.1, 0.2], [0, 0, 7, 5], [7, 5], [1, 1], [0, 0, 64, 64])
    add("flip-y", [1, 0, 0, -1], [0, 100], [0, 0, 8, 8], [8, 8], [1, 1], [0, 0, 100, 100])
    add("tiny-step-small-rect", [1, 0, 0, 1], [0, 0], [0, 0, 0.3, 0.2], [0.3, 0.2], [1, 1], [0, 0, 6, 4])
    add("one-axis-x", rot(17), [2, 3], [0, 0, 4, 4], [4, 4], [1, 0], [-50, -50, 50, 50])
    add("one-axis-y", rot(17), [2, 3], [0, 0, 4, 4], [4, 4], [0, 1], [-50, -50, 50, 50])
    add("no-repeat", rot(17), [2, 3], [0, 0, 4, 4], [4, 4], [0, 0], [0, 0, 8, 8])
    add("point-rect", rot(11), [0.3, 0.7], [0, 0, 3, 2], [3, 2], [1, 1], [5, 5, 5, 5])
    add("far-offset", [1, 0, 0, 1], [1000000.5, -250000.25], [0, 0, 0.1, 0.1], [0.1, 0.1], [1, 0], [1000000, -250001, 1000003, -250000])
    add("far-index-1e6", [1, 0, 0, 1], [0, 0], [0, 0, 0.1, 0.1], [0.1, 0.1], [1, 0], [100000, 0, 100001, 1])
    add("negative-scale", [-0.5, 0, 0, -0.5], [10, 10], [0, 0, 3, 3], [3, 3], [1, 1], [-5, -5, 12, 12])
    add("thin-sliver", [1, 0.001, 0, 1], [0, 0], [0, 0, 100, 0.01], [100, 0.5], [1, 1], [0, 0, 300, 20])
    for k in range(24):
        m = [rnd.uniform(-3, 3) for _ in range(4)]
        while abs(m[0] * m[3] - m[1] * m[2]) < 0.4:
            m = [rnd.uniform(-3, 3) for _ in range(4)]
        t = [rnd.uniform(-50, 50), rnd.uniform(-50, 50)]
        cell = [rnd.uniform(-5, 5), rnd.uniform(-5, 5), rnd.uniform(0.5, 6), rnd.uniform(0.5, 6)]
        steps = [cell[2] if rnd.random() < 0.4 else rnd.uniform(0.7, 7), cell[3] if rnd.random() < 0.4 else rnd.uniform(0.7, 7)]
        # keep the step vectors well conditioned and not tiny so that the window is small
        Q, tau, A, B = quad_and_steps(m, t, cell, steps)
        if abs(A[0] * B[1] - A[1] * B[0]) < 0.5:
            continue
        x0, y0 = rnd.randint(-40, 40), rnd.randint(-40, 40)
        rect = [x0, y0, x0 + rnd.randint(0, 45), y0 + rnd.randint(0, 45)]
        add("random-%02d" % k, m, t, cell, steps, [1, 1] if rnd.random() < 0.85 else [rnd.randint(0, 1), rnd.randint(0, 1)], rect)
    return cases


def translation_cases():
    rnd = random.Random(7)
    out = []
    for k in range(40):
        m = [rnd.uniform(-2, 2) for _ in range(4)]
        t = [rnd.uniform(-1e4, 1e4), rnd.uniform(-1e4, 1e4)]
        cell = [0.0, 0.0, rnd.uniform(0.05, 3), rnd.uniform(0.05, 3)]
        steps = [rnd.choice([0.1, 0.3, 1 / 3, cell[2]]), rnd.choice([0.1, 0.7, 1 / 7, cell[3]])]
        i = rnd.choice([0, 1, -1, 17, 10 ** 6, -10 ** 6, 12345678901, 2 ** 53 - 1, -(2 ** 53 - 1)])
        j = rnd.choice([0, 1, -3, 999999, 10 ** 6, -10 ** 7, 2 ** 40])
        a = (m[0] * steps[0], m[2] * steps[0])
        b = (m[1] * steps[1], m[3] * steps[1])
        w = ((t[0] + float(i) * a[0]) + float(j) * b[0], (t[1] + float(i) * a[1]) + float(j) * b[1])
        # sanity: within a few ulps of the exact value
        ex = F(t[0]) + i * F(a[0]) + j * F(b[0])
        if abs(w[0]) < 1e300:
            assert abs(F(w[0]) - ex) <= F(abs(w[0])) * F(1, 2 ** 50) + F(abs(i * a[0]) + abs(j * b[0])) * F(1, 2 ** 50)
        out.append(("translation-%02d" % k, m, t, cell, steps, i, j, w))
    return out


def period_cases():
    out = []
    pairs = [((1, 0), (0, 1)), ((8, 0), (0, 8)), ((8, 3), (0, 5)), ((0, 4), (-4, 0)), ((5, 2), (-2, 5)), ((3, 0), (1, 1)), ((-3, 7), (2, -5)),
             ((6, 4), (-9, 6)), ((4, 4), (4, -4)), ((7, 1), (1, 7)), ((12, 5), (5, -12)), ((2, 0), (0, 3)), ((1, 1), (1, 2)), ((16, 8), (8, 12)),
             ((2147483647, 1), (2147483646, 1)), ((1 << 20, 0), (0, 1 << 10)), ((100, 37), (-41, 90))]
    for a, b in pairs:
        D = a[0] * b[1] - a[1] * b[0]
        if D == 0 or abs(D) > 2 ** 31:
            continue
        s = math.gcd(abs(a[1]), abs(b[1]))
        p = abs(D) // s

        def in_lat(x, y):
            return (x * b[1] - y * b[0]) % D == 0 and (a[0] * y - a[1] * x) % D == 0

        if p * s <= 2 ** 24:
            r = next(r for r in range(p) if in_lat(r, s))
        else:
            # the search is too long: solve u a1 + v b1 = s by brute force over u in the small range of residues
            r = None
            for u in range(-abs(b[1]) // s - 1, abs(b[1]) // s + 2):
                num = s - u * a[1]
                if num % b[1] == 0 if b[1] else False:
                    v = num // b[1]
                    r = (u * a[0] + v * b[0]) % p
                    break
            assert r is not None
        out.append(("period-%d" % len(out), a, b, p, r, s))
    return out


def main():
    cases = each_cases()
    print("// Generated by tools/gen-lattice-vectors.py (exact rationals).  Do not edit; run the generator.")
    print("// Independent of the library's code.")
    print("#ifndef GHOTI_IO_GSHGEOM_TEST_LATTICE_VECTORS_H")
    print("#define GHOTI_IO_GSHGEOM_TEST_LATTICE_VECTORS_H")
    print("")
    print("#include <cstdint>")
    print("")
    print("namespace lattice_vectors {")
    print("")
    print("// The copies of a use: the map's linear part {m00, m01, m10, m11} and translation, the lattice, the rect {x0, y0, x1, y1},")
    print("// and the rows {j, i_lo, i_hi} in lattice order.")
    print("struct VEach {")
    print("  const char * name;")
    print("  double m[4];")
    print("  double t[2];")
    print("  double cell[4];")
    print("  double steps[2];")
    print("  int repeat[2];")
    print("  int64_t rect[4];")
    print("  int row_count;")
    print("  const int64_t (*rows)[3];")
    print("};")
    print("")
    print("// LT2: the translation ((tau + i a) + j b) of copy (i, j), as the header spells it.")
    print("struct VTranslation {")
    print("  const char * name;")
    print("  double m[4];")
    print("  double t[2];")
    print("  double cell[4];")
    print("  double steps[2];")
    print("  int64_t i;")
    print("  int64_t j;")
    print("  double w[2];")
    print("};")
    print("")
    print("// LT3: the Hermite normal form of Z a + Z b.")
    print("struct VPeriod {")
    print("  const char * name;")
    print("  int64_t a[2];")
    print("  int64_t b[2];")
    print("  int64_t p;")
    print("  int64_t r;")
    print("  int64_t s;")
    print("};")
    print("")
    total_rows = 0
    entries = []
    for idx, (name, m, t, cell, steps, repeat, rect) in enumerate(cases):
        rows = rows_of(m, t, cell, steps, repeat, rect)
        total_rows += len(rows)
        if rows:
            print("static const int64_t e%02d_rows[][3] = {%s};" % (idx, ", ".join("{%d, %d, %d}" % r for r in rows)))
        entries.append((idx, name, m, t, cell, steps, repeat, rect, rows))
    print("")
    print("static const VEach kEach[] = {")
    for idx, name, m, t, cell, steps, repeat, rect, rows in entries:
        print("    {\"%s\", {%s}, {%s}, {%s}, {%s}, {%d, %d}, {%s}, %d, %s}," % (
            name, ", ".join(hexf(v) for v in m), ", ".join(hexf(v) for v in t), ", ".join(hexf(v) for v in cell),
            ", ".join(hexf(v) for v in steps), repeat[0], repeat[1], ", ".join("%d" % v for v in rect), len(rows),
            ("e%02d_rows" % idx) if rows else "nullptr"))
    print("};")
    print("static const int kEachCount = %d;" % len(entries))
    print("")
    tcs = translation_cases()
    print("static const VTranslation kTranslations[] = {")
    for name, m, t, cell, steps, i, j, w in tcs:
        print("    {\"%s\", {%s}, {%s}, {%s}, {%s}, %d, %d, {%s}}," % (name, ", ".join(hexf(v) for v in m), ", ".join(hexf(v) for v in t),
              ", ".join(hexf(v) for v in cell), ", ".join(hexf(v) for v in steps), i, j, ", ".join(hexf(v) for v in w)))
    print("};")
    print("static const int kTranslationCount = %d;" % len(tcs))
    print("")
    pcs = period_cases()
    print("static const VPeriod kPeriods[] = {")
    for name, a, b, p, r, s in pcs:
        print("    {\"%s\", {%d, %d}, {%d, %d}, %d, %d, %d}," % (name, a[0], a[1], b[0], b[1], p, r, s))
    print("};")
    print("static const int kPeriodCount = %d;" % len(pcs))
    print("")
    print("// %d cases, %d rows, %d translations, %d periods." % (len(entries), total_rows, len(tcs), len(pcs)))
    print("")
    print("}  // namespace lattice_vectors")
    print("")
    print("#endif")


if __name__ == "__main__":
    main()
