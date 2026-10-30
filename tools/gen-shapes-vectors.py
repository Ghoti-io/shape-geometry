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
"""Print tests/shapes_vectors.h: the known answers of shapes.h (PR2 to PR5, PR13), in exact rationals.

    python3 tools/gen-shapes-vectors.py > tests/shapes_vectors.h

Shares no code with the library.  Every input is an exact double.  A coordinate is the correctly rounded value of the real
expression the rule names, computed in ``fractions.Fraction`` and rounded once (``float(Fraction)`` rounds to nearest even); a
sum with a zero term is the other term bit for bit, and a sum that is exactly zero is +0.  A rule states a shape once, unflipped,
in the box's terms (a distance from the left or the right edge, from the top or the bottom edge); the flipped forms of PR13 are
that shape under ONE exchange of edges, so the three flipped pairs cannot disagree with the unflipped one.

The rounded rectangle has no flipped form here (see the header: one engine mirrors the radii and not the start).

The cases: boxes (zero width, zero height, negative zero, a point where x + w rounds, huge, the coordinate limit and one past) by the
four flips for the rect and the ellipse; rounded rectangles (one to four round corners, elliptical corners, radii summing exactly to a
side, a corner with one zero component) and the overlap table (each side, one ulp over, equal, one ulp under, and sums that round the
other way in double); the refusals.
"""
import random
import struct
from fractions import Fraction as F

LIMIT = F(2) ** 100


def hexf(x):
    return float(x).hex()


def nextafter(x, up=True):
    b = struct.unpack("<q", struct.pack("<d", x))[0]
    if x == 0:
        return 5e-324 if up else -5e-324
    b += 1 if (x > 0) == up else -1
    return struct.unpack("<d", struct.pack("<q", b))[0]


def rnd(fr):
    """The nearest double to an exact value; inf when it overflows; +0 for zero."""
    try:
        return float(fr)
    except OverflowError:
        return float("inf") if fr > 0 else float("-inf")


def S(*terms):
    """A sum of doubles (each +/- given as a signed double), rounded once."""
    return rnd(sum((F(t) for t in terms), F(0)))


class Box:
    """x, y, w, h: a coordinate is a distance from the near (left, top) or the far (right, bottom) edge, the edge exchanged by a flip."""

    def __init__(self, x, y, w, h, fx, fy):
        self.x, self.y, self.w, self.h, self.fx, self.fy = x, y, w, h, fx, fy

    def X(self, far, a):
        far ^= self.fx
        if far:
            return S(self.x, self.w, -a)
        return self.x if a == 0 else S(self.x, a)

    def Y(self, far, b):
        far ^= self.fy
        if far:
            return S(self.y, self.h, -b)
        return self.y if b == 0 else S(self.y, b)


def too_big(v):
    return not (abs(v) <= LIMIT) or v != v


def finite(v):
    return v == v and abs(v) != float("inf")


def pt(x, y):
    return (x, y)


def rect_elems(b):
    p = [(b.X(0, 0), b.Y(0, 0)), (b.X(1, 0), b.Y(0, 0)), (b.X(1, 0), b.Y(1, 0)), (b.X(0, 0), b.Y(1, 0))]
    return [("M", p[0])] + [("L", q) for q in p[1:]] + [("Z", p[0])], p


def arc(start, end, c, u, v, t0, circ):
    return ("A", start, end, c, u, v, t0, 90.0, circ)


def ellipse_elems(cx, cy, rx, ry, fx, fy):
    sx = -rx if fx else rx
    sy = -ry if fy else ry
    c = pt(cx, cy)
    u = pt(sx, 0.0)
    v = pt(0.0, sy)
    ends = [pt(S(cx, sx), cy), pt(cx, S(cy, sy)), pt(S(cx, -sx), cy), pt(cx, S(cy, -sy))]
    circ = 1 if rx == ry else 0
    el = [("M", ends[0])]
    for k in range(4):
        el.append(arc(ends[k], ends[(k + 1) % 4], c, u, v, 90.0 * k, circ))
    el.append(("Z", ends[0]))
    return el, ends


def rrect_elems(x, y, w, h, r):
    """r: four (rx, ry), top-left, top-right, bottom-right, bottom-left; unflipped."""
    b = Box(x, y, w, h, 0, 0)
    c = [(0.0, 0.0) if (a == 0 or d == 0) else (a, d) for a, d in r]
    if all(k == (0.0, 0.0) for k in c):
        return rect_elems(b)[0], None
    (tlx, tly), (trx, try_), (brx, bry), (blx, bly) = c
    start = pt(b.X(0, tlx), b.Y(0, 0))
    el = [("M", start), ("L", pt(b.X(1, trx), b.Y(0, 0)))]
    if trx:
        e = pt(b.X(1, 0), b.Y(0, try_))
        el.append(arc(el[-1][1], e, pt(b.X(1, trx), b.Y(0, try_)), pt(trx, 0.0), pt(0.0, try_), -90.0, 1 if trx == try_ else 0))
    el.append(("L", pt(b.X(1, 0), b.Y(1, bry))))
    if brx:
        e = pt(b.X(1, brx), b.Y(1, 0))
        el.append(arc(el[-1][1], e, pt(b.X(1, brx), b.Y(1, bry)), pt(brx, 0.0), pt(0.0, bry), 0.0, 1 if brx == bry else 0))
    el.append(("L", pt(b.X(0, blx), b.Y(1, 0))))
    if blx:
        e = pt(b.X(0, 0), b.Y(1, bly))
        el.append(arc(el[-1][1], e, pt(b.X(0, blx), b.Y(1, bly)), pt(blx, 0.0), pt(0.0, bly), 90.0, 1 if blx == bly else 0))
    el.append(("L", pt(b.X(0, 0), b.Y(0, tly))))
    if tlx:
        el.append(arc(el[-1][1], start, pt(b.X(0, tlx), b.Y(0, tly)), pt(tlx, 0.0), pt(0.0, tly), 180.0, 1 if tlx == tly else 0))
    el.append(("Z", start))
    return el, c


def overlap(w, h, r):
    c = [(0.0, 0.0) if (a == 0 or d == 0) else (a, d) for a, d in r]
    return (F(c[0][0]) + F(c[1][0]) > F(w) or F(c[2][0]) + F(c[3][0]) > F(w)
            or F(c[3][1]) + F(c[0][1]) > F(h) or F(c[1][1]) + F(c[2][1]) > F(h))


def control(a):
    """The control point of a quarter piece: the meeting of the tangents at its ends.  Both tangents are parallel to an axis (the
    half-diameters are), so it is the corner of the box of the two named ends: a piece whose start tangent is vertical (t0 = 0 or 180)
    takes the start's x and the end's y, one whose start tangent is horizontal (t0 = 90 or -90) the end's x and the start's y."""
    _, s, e, c, u, v, t0, sw, circ = a
    k = int(t0) // 90 % 4
    return pt(s[0], e[1]) if k % 2 == 0 else pt(e[0], s[1])


cases = []  # (kind, a[4], r[8], fx, fy, result, elements)


def add_rect(x, y, w, h, fx, fy):
    res = "GSHGEOM_OK"
    el = []
    if not all(finite(v) for v in (x, y, w, h)) or w < 0 or h < 0:
        res = "GSHGEOM_ERR_INVALID"
    else:
        b = Box(x, y, w, h, fx, fy)
        el, p = rect_elems(b)
        if too_big(w) or too_big(h) or any(too_big(c) for q in p for c in q):
            res, el = "GSHGEOM_ERR_LIMIT", []
    cases.append((0, [x, y, w, h], [0.0] * 8, fx, fy, res, el))


def add_rrect(x, y, w, h, r):
    res = "GSHGEOM_OK"
    el = []
    rr = [v for k in r for v in k]
    if not all(finite(v) for v in (x, y, w, h) + tuple(rr)) or w < 0 or h < 0 or any(v < 0 for v in rr) or overlap(w, h, r):
        res = "GSHGEOM_ERR_INVALID"
    else:
        el, _ = rrect_elems(x, y, w, h, r)
        pts = [q for e in el for q in e[1:4] if isinstance(q, tuple)]
        radii = [v for k in r for v in k]
        if too_big(w) or too_big(h) or any(too_big(v) for v in radii) or any(too_big(c) for q in pts for c in q):
            res, el = "GSHGEOM_ERR_LIMIT", []
    cases.append((1, [x, y, w, h], rr, 0, 0, res, el))


def add_ellipse(cx, cy, rx, ry, fx, fy):
    res = "GSHGEOM_OK"
    el = []
    if not all(finite(v) for v in (cx, cy, rx, ry)) or not rx > 0 or not ry > 0:
        res = "GSHGEOM_ERR_INVALID"
    else:
        el, ends = ellipse_elems(cx, cy, rx, ry, fx, fy)
        if any(too_big(v) for v in (cx, cy, rx, ry)) or any(too_big(c) for q in ends for c in q):
            res, el = "GSHGEOM_ERR_LIMIT", []
    cases.append((2, [cx, cy, rx, ry], [0.0] * 8, fx, fy, res, el))


BIG = 2.0 ** 100
boxes = [(10.0, 20.0, 100.0, 50.0), (0.0, 0.0, 1.0, 1.0), (-0.0, -0.0, 0.0, 0.0), (5.0, 7.0, 0.0, 3.0), (5.0, 7.0, 3.0, 0.0),
         (0.1, 0.2, 0.3, 0.7), (-3.5, -4.25, 8.75, 2.5), (1e16, 1e16, 3.0, 3.0), (1e16, -1e16, 1.0, 1.0), (-0.0, 4.0, 0.0, 2.0),
         (2.0 ** 99, -(2.0 ** 99), 2.0 ** 99, 2.0 ** 99), (2.0 ** 99, 0.0, 2.0 ** 99 + 2.0 ** 47, 1.0), (2.0 ** 99, 0.0, 2.0 ** 99 + 2.0 ** 48, 1.0), (BIG, 0.0, 1.0, 1.0),
         (-BIG, -BIG, BIG, BIG), (-BIG, 0.0, nextafter(BIG), 1.0), (1e-300, 1e-300, 1e-300, 1e-300),
         (5e-324, -5e-324, 5e-324, 5e-324), (1.0, 1.0, 1e308, 1.0), (0.0, 0.0, -1.0, 1.0), (0.0, 0.0, 1.0, -1.0),
         (float("nan"), 0.0, 1.0, 1.0), (0.0, float("inf"), 1.0, 1.0), (0.0, 0.0, float("inf"), 1.0)]
for b in boxes:
    for fx in (0, 1):
        for fy in (0, 1):
            add_rect(*b, fx, fy)

R = lambda *a: [(a[i], a[i + 1]) for i in range(0, 8, 2)]
rr_boxes = [(10.0, 20.0, 100.0, 60.0), (0.0, 0.0, 40.0, 40.0), (-7.5, 3.25, 31.0, 17.0), (1e16, 1e16, 64.0, 64.0)]
radii = [R(0, 0, 0, 0, 0, 0, 0, 0), R(10, 10, 0, 0, 0, 0, 0, 0), R(0, 0, 8, 6, 0, 0, 0, 0), R(0, 0, 0, 0, 5, 5, 0, 0),
         R(0, 0, 0, 0, 0, 0, 12, 4), R(5, 5, 7, 3, 0, 0, 9, 9), R(4, 4, 4, 4, 4, 4, 4, 4), R(10, 5, 6, 12, 3, 8, 14, 2),
         R(3, 0, 5, 5, 0, 7, 4, 4), R(0.25, 0.5, 0.125, 0.75, 1.5, 2.5, 0.1, 0.3), R(5, 5, 5, 5, 5, 5, 5, 5)]
for x, y, w, h in rr_boxes:
    for r in radii:
        add_rrect(x, y, w, h, r)
# radii that sum exactly to a side: a stadium (zero-length lines), a circle, and a one-side-only case
add_rrect(0.0, 0.0, 20.0, 10.0, R(10, 5, 10, 5, 10, 5, 10, 5))
add_rrect(0.0, 0.0, 10.0, 10.0, R(5, 5, 5, 5, 5, 5, 5, 5))
add_rrect(0.0, 0.0, 20.0, 10.0, R(8, 5, 12, 5, 12, 5, 8, 5))
add_rrect(3.0, 4.0, 0.0, 0.0, R(0, 0, 0, 0, 0, 0, 0, 0))
add_rrect(3.0, 4.0, 0.0, 6.0, R(0, 3, 0, 3, 0, 3, 0, 3))
add_rrect(-0.0, -0.0, 8.0, 8.0, R(0, 0, 4, 4, 0, 0, 4, 4))
# the overlap table, each side, w = 1 with 0.1 + 0.9: the doubles sum above 1 in the reals and round to 1
a, b = 0.1, 0.9
for one in (1.0, nextafter(1.0), nextafter(1.0, False)):
    for side in range(4):
        r = [(0.0, 0.0)] * 4
        r = [list(k) for k in r]
        # top: r0.x + r1.x; right: r1.y + r2.y; bottom: r2.x + r3.x; left: r3.y + r0.y
        slots = [((0, 0), (1, 0)), ((1, 1), (2, 1)), ((2, 0), (3, 0)), ((3, 1), (0, 1))][side]
        for (i, j), val in zip(slots, (a, b)):
            r[i][j] = val
        # the partners of a one-component corner must be nonzero to be round: give them the other component 1/8 of the side
        for (i, j) in slots:
            r[i][1 - j] = 0.125
        rt = [tuple(k) for k in r]
        add_rrect(0.0, 0.0, one if side in (0, 2) else 10.0, one if side in (1, 3) else 10.0, rt)
# a rounded sum that is exactly the side while the reals are not: 2^53 + 1 style
big = 2.0 ** 53
add_rrect(0.0, 0.0, big, 100.0, R(big - 1, 1, 1.0 + 2.0 ** -52, 1, 0, 0, 0, 0))
add_rrect(0.0, 0.0, big, 100.0, R(big - 2, 1, 2.0, 1, 0, 0, 0, 0))
add_rrect(0.0, 0.0, 1.0, 1.0, R(nextafter(0.5), 0.25, 0.5, 0.25, 0, 0, 0, 0))
add_rrect(0.0, 0.0, 1.0, 1.0, R(0.5, 0.25, 0.5, 0.25, 0, 0, 0, 0))
# a corner with one zero component is square and does not count toward an overlap
add_rrect(0.0, 0.0, 40.0, 40.0, R(30, 0, 30, 5, 0, 0, 0, 0))
add_rrect(0.0, 0.0, 40.0, 40.0, R(30, 5, 30, 5, 0, 0, 0, 0))
# refusals
add_rrect(0.0, 0.0, -1.0, 5.0, R(0, 0, 0, 0, 0, 0, 0, 0))
add_rrect(0.0, 0.0, 5.0, -1.0, R(0, 0, 0, 0, 0, 0, 0, 0))
add_rrect(0.0, 0.0, 5.0, 5.0, R(-1, 1, 0, 0, 0, 0, 0, 0))
add_rrect(0.0, 0.0, 5.0, 5.0, R(1, -1, 0, 0, 0, 0, 0, 0))
add_rrect(0.0, 0.0, 5.0, 5.0, R(float("nan"), 1, 0, 0, 0, 0, 0, 0))
add_rrect(0.0, 0.0, 5.0, 5.0, R(1, float("inf"), 0, 0, 0, 0, 0, 0))
add_rrect(float("nan"), 0.0, 5.0, 5.0, R(1, 1, 0, 0, 0, 0, 0, 0))
add_rrect(0.0, 0.0, BIG, 5.0, R(1, 1, 0, 0, 0, 0, 0, 0))
add_rrect(BIG, 0.0, 1.0, 5.0, R(1, 1, 0, 0, 0, 0, 0, 0))
add_rrect(BIG, 0.0, 1.0, 5.0, R(0, 0, 0, 0, 0, 0, 0, 0))
add_rrect(-BIG, -BIG, BIG, BIG, R(BIG / 2, BIG / 2, 0, 0, 0, 0, 0, 0))
add_rrect(1e-300, 1e-300, 1e-300, 1e-300, R(1e-301, 1e-301, 0, 0, 0, 0, 0, 0))

ells = [(0.0, 0.0, 1.0, 1.0), (10.0, 20.0, 5.0, 3.0), (-4.5, 8.25, 0.5, 40.0), (1e16, 1e16, 3.0, 3.0), (-0.0, -0.0, 2.0, 2.0),
        (0.1, 0.2, 0.3, 0.7), (1.0, 1.0, 1.0, nextafter(1.0)), (0.0, 0.0, 1e-300, 1e-300), (0.0, 0.0, 5e-324, 5e-324),
        (0.0, 0.0, 2.0 ** 99, 2.0 ** 99), (2.0 ** 99, 0.0, 2.0 ** 99, 1.0), (BIG, 0.0, 1.0, 1.0), (BIG, 0.0, 2.0 ** 48, 1.0), (-BIG, 0.0, 2.0 ** 48, 1.0), (0.0, 0.0, nextafter(BIG), 1.0),
        (0.0, 0.0, 0.0, 1.0), (0.0, 0.0, 1.0, 0.0), (0.0, 0.0, -1.0, 1.0), (0.0, 0.0, 1.0, -1.0), (float("nan"), 0.0, 1.0, 1.0),
        (0.0, 0.0, float("inf"), 1.0), (0.0, 0.0, 1e308, 1.0), (1.0, 1.0, 1e300, 1e-300)]
for e in ells:
    for fx in (0, 1):
        for fy in (0, 1):
            add_ellipse(*e, fx, fy)

# a random sweep, to put a coordinate that rounds in every position
rng = random.Random(0x5A2)


def rd(scale):
    return rng.choice([1, -1]) * (rng.random() * scale)


for _ in range(160):
    k = 10.0 ** rng.randint(-3, 9)
    add_rect(rd(k), rd(k), rng.random() * k, rng.random() * k, rng.randint(0, 1), rng.randint(0, 1))
    add_ellipse(rd(k), rd(k), rng.random() * k + 1e-9, rng.random() * k + 1e-9, rng.randint(0, 1), rng.randint(0, 1))
    w, h = rng.random() * k + 1.0, rng.random() * k + 1.0
    fr = [rng.random() * 0.25 for _ in range(8)]
    rr = [(w * fr[0], h * fr[1]), (w * fr[2], h * fr[3]), (w * fr[4], h * fr[5]), (w * fr[6], h * fr[7])]
    rr = [(a if rng.random() > 0.15 else 0.0, d if rng.random() > 0.15 else 0.0) for a, d in rr]
    add_rrect(rd(k), rd(k), w, h, rr)

out = []
w = out.append
w("/* Generated by tools/gen-shapes-vectors.py; do not edit. */")
w("/* SPDX-License-Identifier: LGPL-3.0-only */")
w("#ifndef GHOTI_IO_GSHGEOM_TESTS_SHAPES_VECTORS_H")
w("#define GHOTI_IO_GSHGEOM_TESTS_SHAPES_VECTORS_H")
w("")
w("#include <ghoti.io/shape-geometry/core.h>")
w("")
w("#include <cstddef>")
w("")
w("namespace shapes_vectors {")
w("")
w("/* verb: 'M', 'L', 'A' (an arc: s start, e end, c centre, u and v half-diameters, t0 and sweep in degrees, circ, and ctl the control point")
w("   of its one piece: the corner of the box of the two ends), 'Z' (e is the subpath's start).  Every number is the oracle's, rounded once. */")
w("struct Elem {")
w("  char verb;")
w("  double e[2], s[2], c[2], u[2], v[2], t0, sweep, ctl[2];")
w("  int circ;")
w("};")
w("/* kind: 0 rect, 1 rounded rect, 2 ellipse.  a: x y w h (or cx cy rx ry); r: the eight radius components, TL TR BR BL, rx then ry. */")
w("struct Case {")
w("  int kind;")
w("  double a[4], r[8];")
w("  int fx, fy;")
w("  GSHGEOM_Result result;")
w("  size_t first, count;")
w("};")
w("")
flat = []
rows = []
for kind, a, r, fx, fy, res, el in cases:
    rows.append((kind, a, r, fx, fy, res, len(flat), len(el)))
    for e in el:
        if e[0] in "MLZ":
            flat.append("{'%s', {%s, %s}, {0,0}, {0,0}, {0,0}, {0,0}, 0, 0, {0,0}, 0}" % (e[0], hexf(e[1][0]), hexf(e[1][1])))
        else:
            ct = control(e)
            flat.append("{'A', {%s, %s}, {%s, %s}, {%s, %s}, {%s, %s}, {%s, %s}, %s, %s, {%s, %s}, %d}" % (
                hexf(e[2][0]), hexf(e[2][1]), hexf(e[1][0]), hexf(e[1][1]), hexf(e[3][0]), hexf(e[3][1]), hexf(e[4][0]),
                hexf(e[4][1]), hexf(e[5][0]), hexf(e[5][1]), hexf(e[6]), hexf(e[7]), hexf(ct[0]), hexf(ct[1]), e[8]))
w("static const Elem elems[] = {")
for f in flat:
    w("    " + f + ",")
w("};")
w("static const Case cases[] = {")
for kind, a, r, fx, fy, res, first, n in rows:
    def hx(v):
        if v != v:
            return "__builtin_nan(\"\")"
        if v in (float("inf"), float("-inf")):
            return ("" if v > 0 else "-") + "__builtin_inf()"
        return hexf(v)
    w("    {%d, {%s}, {%s}, %d, %d, %s, %d, %d}," % (kind, ", ".join(hx(v) for v in a), ", ".join(hx(v) for v in r), fx, fy, res, first, n))
w("};")
w("static const size_t case_count = sizeof cases / sizeof cases[0];")
w("")
w("} // namespace shapes_vectors")
w("")
w("#endif")
print("\n".join(out))
