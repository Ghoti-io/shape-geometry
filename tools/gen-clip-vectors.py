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
"""Print tests/clip_vectors.h: the known answers of clip.h (CK1 to CK6), in exact rationals.

    python3 tools/gen-clip-vectors.py > tests/clip_vectors.h

Shares no code with the library.  Every number is an integer (a 26.6 coordinate).  The side of a point is the sign of an exact
integer determinant (Python integers); a crossing is the exact rational ``Fraction`` intersection of the edge with the line
through two points of the quadrilateral, rounded once to the nearest integer with ties upward (``floor(v + 1/2)``).  The
rational is formed from the edge's end points as the edge is given (not in the library's lexicographic order), and from the
line's points in the quadrilateral's order, so that the library's canonical order is checked against an ordinary
formulation: the two agree only because the exact rational does not depend on the order.

The cut is the one clip.h states: one side of the quadrilateral at a time, in order, a side of zero length skipped; a
contour with every point inside or on the side's line is unchanged, one with none is dropped, otherwise the closed polygon
is walked from its last point to its first and on, a crossing emitted when the two ends are strictly on opposite sides and
the end kept when it is not strictly outside.  A contour left with fewer than three points after a cut is dropped.  The
flag of a contour no cut changed is its input's; a cut contour is closed.

The cases: the handwritten shapes of the test (a square cell and rectangles, a fan, a bowtie, a nested pair, contours on
and outside the cell edge, rotated and skewed cells, both orientations, a triangle given as four points), ties at a half
unit in each coordinate and each direction, extremes of the int32 range, and seeded random polygons against random
convex quadrilaterals (seed and order of draws below).
"""
import random
from fractions import Fraction as F

I32 = 2 ** 31


def side(a, b, p):
    return (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0])


def sgn(v):
    return (v > 0) - (v < 0)


def half_up(v):
    """floor(v + 1/2) for a Fraction."""
    w = v + F(1, 2)
    return w.numerator // w.denominator


def crossing(a, b, l0, l1):
    fa, fb = side(l0, l1, a), side(l0, l1, b)
    assert sgn(fa) * sgn(fb) < 0
    t = F(fa, fa - fb)
    return (half_up(a[0] + (b[0] - a[0]) * t), half_up(a[1] + (b[1] - a[1]) * t))


def quad_sides(q):
    """The cutting sides with the interior sign, or None for a quad the library must refuse."""
    turns = [sgn(side(q[(i + 3) % 4], q[i], q[(i + 1) % 4])) for i in range(4)]
    signs = {t for t in turns if t}
    if len(signs) != 1 or (turns[1] == 0 and turns[3] == 0):
        return None
    sigma = signs.pop()
    return [(q[i], q[(i + 1) % 4], sigma) for i in range(4) if q[i] != q[(i + 1) % 4]]


def cut(poly, s):
    a, b, sigma = s
    cls = [sgn(side(a, b, p)) * sigma for p in poly]
    if all(c >= 0 for c in cls):
        return poly, False
    if all(c < 0 for c in cls):
        return None, True
    out = []
    n = len(poly)
    for k in range(n):
        s_, e = poly[k - 1], poly[k]
        cs, ce = cls[k - 1], cls[k]
        if cs * ce < 0:
            out.append(crossing(s_, e, a, b))
        if ce >= 0:
            out.append(e)
    return out, True


def clip(contours, q):
    """contours: list of (points, closed).  Returns the clipped list."""
    sides = quad_sides(q)
    out = []
    for pts, closed in contours:
        cur, changed = list(pts), False
        for s in sides:
            nxt, ch = cut(cur, s)
            if nxt is None:
                cur = None
                break
            if ch:
                cur, changed = nxt, True
        if cur is None:
            continue
        if changed:
            if len(cur) < 3:
                continue
            out.append((cur, 1))
        else:
            out.append((cur, closed))
    return out


# --- the cases -------------------------------------------------------------------------------------------------------

CASES = []


def add(name, q, contours):
    CASES.append((name, q, contours))


def rect(x0, y0, x1, y1):
    return ([(x0, y0), (x1, y0), (x1, y1), (x0, y1)], 1)


S64 = [(0, 0), (64, 0), (64, 64), (0, 64)]
add("cell-rect-overlapping-corner", S64, [rect(32, 32, 96, 96)])
add("cell-rect-crossing-two-sides", S64, [rect(-32, 16, 96, 48)])
add("cell-rect-covering", S64, [rect(-10, -10, 100, 100)])
add("cell-rect-inside", S64, [rect(8, 8, 40, 40)])
add("cell-rect-outside", S64, [rect(100, 0, 140, 40)])
add("cell-rect-touching-side", S64, [rect(64, 0, 128, 64)])
add("cell-rect-on-boundary", S64, [rect(0, 0, 64, 64)])
add("cell-rect-touching-corner", S64, [rect(64, 64, 100, 100)])
add("cell-rect-reversed", S64, [([(0, 0), (0, 64), (64, 64), (64, 0)], 1), rect(-5, -5, 30, 30)[0:1] + (0,)])
add("cell-fan", S64, [([(32, 32), (200, -40), (200, 100)], 1), ([(32, 32), (200, 100), (-60, 200)], 1),
                      ([(32, 32), (-60, 200), (-100, 0)], 1), ([(32, 32), (-100, 0), (200, -40)], 1)])
add("cell-bowtie", S64, [([(-20, -20), (90, 60), (90, -20), (-20, 60)], 1)])
add("cell-bowtie-within", S64, [([(8, 8), (56, 40), (56, 8), (8, 40)], 1)])
add("cell-nested-hole", S64, [rect(-30, -30, 100, 100), ([(20, 20), (20, 80), (80, 80), (80, 20)], 1)])
add("cell-nested-hole-inside", S64, [rect(4, 4, 60, 60), ([(16, 16), (16, 48), (48, 48), (48, 16)], 1)])
add("cell-overlapping-same-winding", S64, [rect(-20, -20, 40, 40), rect(20, 20, 90, 90)])
add("cell-open-subpath-inside", S64, [([(10, 10), (50, 12), (30, 50)], 0)])
add("cell-open-subpath-cut", S64, [([(10, 10), (90, 12), (30, 90)], 0)])
add("cell-one-two-point-contours", S64, [([(10, 10)], 0), ([(10, 10), (20, 20)], 1), ([(100, 100)], 1),
                                         ([(10, 10), (100, 100)], 1)])
add("cell-on-the-line-detour", S64, [([(0, -10), (0, 0), (0, 64), (0, 80), (-20, 40)], 1)])
add("cell-vertex-on-corner", S64, [([(0, 0), (-30, 40), (40, 90), (64, 64), (90, 20)], 1)])
add("cell-touches-side-at-point", S64, [([(20, 20), (64, 30), (20, 40)], 1), ([(80, 10), (64, 30), (80, 50)], 1)])
add("cell-zero-area-along-edge", S64, [([(0, 0), (30, 0), (64, 0), (30, 0)], 1)])
add("cell-all-collinear-inside", S64, [([(8, 8), (16, 16), (24, 24)], 1)])
# rotated and skewed cells, both orientations, a triangle as four points, a repeated point in the middle
DIAMOND = [(64, 0), (128, 64), (64, 128), (0, 64)]
add("diamond-square", DIAMOND, [rect(20, 20, 110, 110)])
add("diamond-triangle", DIAMOND, [([(10, 40), (140, 70), (60, 150)], 1)])
add("diamond-reversed-quad", DIAMOND[::-1], [rect(20, 20, 110, 110)])
add("skew-cell", [(0, 0), (100, 30), (130, 90), (30, 60)], [rect(10, -20, 90, 70), ([(60, 100), (140, 20), (130, 130)], 1)])
add("skew-cell-reversed", [(0, 0), (30, 60), (130, 90), (100, 30)], [rect(10, -20, 90, 70)])
add("triangle-as-four-points", [(0, 0), (100, 0), (0, 100), (0, 100)], [rect(-10, -10, 60, 60)])
add("triangle-repeat-first", [(0, 0), (0, 0), (100, 0), (0, 100)], [rect(-10, -10, 60, 60)])
add("triangle-repeat-middle", [(0, 0), (100, 0), (100, 0), (0, 100)], [rect(40, -10, 160, 60)])
add("triangle-with-corner-on-side", [(0, 0), (50, 0), (100, 0), (0, 100)], [rect(-10, -10, 60, 60)])
# half-unit ties: a vertical line x = 1 through (1, 0)-(1, 64); the edge (0,0)-(2,1) crosses it at y = 1/2
TIE = [(1, 0), (1, 64), (-60, 64), (-60, 0)]
add("tie-y-up", TIE, [([(0, 0), (2, 1), (2, 10), (0, 10)], 1)])
add("tie-y-down", TIE, [([(0, 0), (2, -1), (2, -10), (0, -10)], 1)])
add("tie-y-up-reversed-edge", TIE, [([(2, 1), (0, 0), (0, 10), (2, 10)], 1)])
add("tie-y-down-reversed-edge", TIE, [([(2, -1), (0, 0), (0, -10), (2, -10)], 1)])
add("tie-x", [(0, 1), (64, 1), (64, 60), (0, 60)], [([(0, 0), (1, 2), (10, 2), (10, 0)], 1)])
add("tie-x-negative", [(0, 1), (64, 1), (64, 60), (0, 60)], [([(0, 0), (-1, 2), (-10, 2), (-10, 0)], 1)])
add("tie-both", [(1, 1), (60, 1), (60, 60), (1, 60)], [([(0, 0), (2, 2), (-3, 3)], 1), ([(0, 0), (2, 1), (4, 1)], 1)])
add("tie-line-reversed", [(1, 64), (1, 0), (60, 0), (60, 64)], [([(0, 0), (2, 1), (2, 10), (0, 10)], 1)])
add("tie-negative-coordinates", [(-1, 0), (-1, 64), (-60, 64), (-60, 0)], [([(0, 0), (-2, 1), (-2, 10), (0, 10)], 1)])
add("tie-diagonal-line", [(0, 0), (64, 64), (0, 64), (0, 0)], [([(-1, 0), (3, 0), (3, 4), (-1, 4)], 1)])
# the extremes of int32
M = I32 - 1
add("extreme-full-range", [(-M - 1, -M - 1), (M, -M - 1), (M, M), (-M - 1, M)], [rect(-M, -M, M - 5, M - 7)])
add("extreme-edge-across", [(-100, -100), (100, -100), (100, 100), (-100, 100)],
    [([(-M - 1, -M - 1), (M, M), (M, -M - 1)], 1)])
add("extreme-sliver", [(0, 0), (1000, 0), (1000, 1000), (0, 1000)], [([(-M - 1, 1), (M, 2), (M, 3)], 1)])
add("extreme-thin-quad", [(0, 0), (M, 0), (M, 1), (0, 1)], [([(-1000, -M - 1), (M, M), (5, 5)], 1)])
add("extreme-diagonal-cell", [(-M, -M), (M, M), (M, M - 1), (-M, -M - 1)], [rect(-1000, -1000, 1000, 1000)])
add("empty-path", S64, [])

# seeded random polygons against random convex quadrilaterals.  The generator is xorshift64*, drawn in a stated order per
# case: the quadrilateral's centre (2 draws), its half-sizes (2), its rotation (1), the number of contours (1), and for
# each contour the number of points (1) and each point's two coordinates (2 each), all reduced by modulo to the ranges
# below.  The same 64-bit integers on every machine; the Python and C++ sides never need to agree on a draw because the
# cases are written out in full in the header.


class Rng:
    def __init__(self, seed):
        self.s = seed

    def next(self):
        s = self.s
        s ^= s >> 12
        s ^= (s << 25) & (2 ** 64 - 1)
        s ^= s >> 27
        self.s = s
        return (s * 2685821657736338717) & (2 ** 64 - 1)

    def below(self, n):
        return self.next() % n


def random_cases(count, seed, span, big):
    rng = Rng(seed)
    for i in range(count):
        cx, cy = rng.below(span) - span // 2, rng.below(span) - span // 2
        hx, hy = 1 + rng.below(span // 3), 1 + rng.below(span // 3)
        shear = rng.below(7) - 3
        q = [(cx - hx, cy - hy), (cx + hx, cy - hy + shear), (cx + hx + shear, cy + hy), (cx - hx, cy + hy - shear)]
        if rng.below(2):
            q = q[::-1]
        if rng.below(5) == 0:
            q[3] = q[2]
        if quad_sides(q) is None:
            q = [(cx - hx, cy - hy), (cx + hx, cy - hy), (cx + hx, cy + hy), (cx - hx, cy + hy)]
        contours = []
        for _ in range(1 + rng.below(3)):
            n = 1 + rng.below(9)
            pts = [(cx + rng.below(2 * span) - span, cy + rng.below(2 * span) - span) for _ in range(n)]
            contours.append((pts, rng.below(4) != 0))
        add("random-%s-%d" % (big, i), q, contours)


random_cases(150, 0x9E3779B97F4A7C15, 300, "small")
random_cases(60, 0xD1B54A32D192ED03, 1 << 24, "large")


def random_extreme(count, seed):
    rng = Rng(seed)
    for i in range(count):
        def c():
            return rng.below(1 << 32) - I32
        q = None
        while q is None or quad_sides(q) is None:
            x0, y0 = c(), c()
            x1, y1 = c(), c()
            q = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
            if x0 == x1 or y0 == y1:
                q = None
        contours = [([(c(), c()) for _ in range(1 + rng.below(6))], 1) for _ in range(1 + rng.below(2))]
        add("random-extreme-%d" % i, q, contours)


random_extreme(40, 0xA24BAED4963EE407)

# --- output ----------------------------------------------------------------------------------------------------------


def main():
    pts_in, subs_in, pts_out, subs_out, cases = [], [], [], [], []
    for name, q, contours in CASES:
        want = clip(contours, q)
        first_in = len(subs_in)
        for pts, closed in contours:
            subs_in.append((len(pts_in), len(pts), int(bool(closed))))
            pts_in.extend(pts)
        first_out = len(subs_out)
        for pts, closed in want:
            subs_out.append((len(pts_out), len(pts), int(bool(closed))))
            pts_out.extend(pts)
        cases.append((name, q, first_in, len(contours), first_out, len(want)))
    print("/* Generated by tools/gen-clip-vectors.py; do not edit. */")
    print("/* SPDX-License-Identifier: LGPL-3.0-only */")
    print("#ifndef GHOTI_IO_GSHGEOM_TESTS_CLIP_VECTORS_H")
    print("#define GHOTI_IO_GSHGEOM_TESTS_CLIP_VECTORS_H")
    print()
    print("#include <cstddef>")
    print("#include <cstdint>")
    print()
    print("namespace clip_vectors {")
    print()
    print("struct Pt {")
    print("  int32_t x, y;")
    print("};")
    print("/* A contour: its first point in the points array, the number of points, its closed flag. */")
    print("struct Sub {")
    print("  size_t first, count;")
    print("  int closed;")
    print("};")
    print("/* A case: a name, the quadrilateral, and the input and expected contours as ranges of the contour arrays. */")
    print("struct Case {")
    print("  const char * name;")
    print("  Pt quad[4];")
    print("  size_t in_first, in_count, out_first, out_count;")
    print("};")
    print()

    def pts_arr(name, pts):
        print("static const Pt %s[] = {" % name)
        for i in range(0, len(pts), 6):
            print("    " + " ".join("{%d, %d}," % p for p in pts[i:i + 6]))
        print("    {0, 0}};  /* a sentinel, so the array is never empty */")

    def subs_arr(name, subs):
        print("static const Sub %s[] = {" % name)
        for i in range(0, len(subs), 4):
            print("    " + " ".join("{%d, %d, %d}," % s for s in subs[i:i + 4]))
        print("    {0, 0, 0}};")

    pts_arr("in_points", pts_in)
    subs_arr("in_subs", subs_in)
    pts_arr("out_points", pts_out)
    subs_arr("out_subs", subs_out)
    print()
    print("static const Case cases[] = {")
    for name, q, fi, ci, fo, co in cases:
        print('    {"%s", {%s}, %d, %d, %d, %d},' % (name, ", ".join("{%d, %d}" % p for p in q), fi, ci, fo, co))
    print("};")
    print()
    print("static const size_t case_count = %d;" % len(cases))
    print()
    print("}  // namespace clip_vectors")
    print()
    print("#endif")


if __name__ == "__main__":
    main()
