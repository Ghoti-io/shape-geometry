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
"""Print the known answers of the hash spelling and the predicates (CAP-2, CAP-4, CAP-13, AD-21).

This is a second statement of the spelling in include/ghoti.io/shape-geometry/
hash.h and path.h, written with Python's `struct` and none of the library's
code, so that a test comparing the two compares two readings of the text.  It
writes tests/known_answers.h, which is committed; run it after changing the
spelling, with the fixtures below (the same numbers as test_hash.cpp builds
its paths from).

    python3 tools/gen-hash-vectors.py > tests/known_answers.h

Story 6 adds the map fixtures (their hash and dump, and the kind tag each is
expected to be classified as, stated here by hand from the rules in map.h) and
the orientation vectors: points chosen to be exactly or nearly degenerate, with
the sign decided by Python's `fractions` (exact rational arithmetic), and the
number of them on which the naive double expression is wrong, so that a test
can show the vectors are discriminating.
"""
import math
import random
import struct
import sys
from fractions import Fraction

BASIS = 0xCBF29CE484222325
PRIME = 0x100000001B3
M64 = (1 << 64) - 1


def fnv(data, state=BASIS):
    for b in data:
        state = ((state ^ b) * PRIME) & M64
    return state


# FNV's published test vectors (the reference suite at isthe.com/chongo/tech/
# comp/fnv/): "" is the basis, and these two are in the reference listing.
assert fnv(b"") == 0xCBF29CE484222325
assert fnv(b"a") == 0xAF63DC4C8601EC8C
assert fnv(b"foobar") == 0x85944171F73967E8


def u8(v):
    return struct.pack("<B", v)


def u16(v):
    return struct.pack("<H", v)


def u32(v):
    return struct.pack("<I", v)


def u64(v):
    return struct.pack("<Q", v)


def i32(v):
    return struct.pack("<i", v)


def i64(v):
    return struct.pack("<q", v)


def f64(v):
    return struct.pack("<d", v)


def bits(v):
    return struct.unpack("<Q", struct.pack("<d", v))[0]


def string(s):
    b = s.encode("ascii")
    return u32(len(b)) + b


# The fixture paths.  ("M", [points]) and so on; points are 3-tuples whose z
# is used only in 3D.
def fixture(dim):
    def P(x, y):
        return (x, y, (x + y) / 4.0)
    return [
        ("M", [P(1.5, -2.25)], None),
        ("L", [P(3.0, 4.0)], None),
        ("Q", [P(5.0, 6.5), P(7.0, 8.0)], None),
        ("C", [P(9.0, 10.0), P(11.25, 12.0), P(13.0, 14.0)], None),
        ("RQ", [P(15.0, 16.0), P(17.0, 18.0)], [1.0, 0.5, 2.0]),
        ("RC", [P(19.0, 20.0), P(21.0, 22.0), P(23.0, 24.1)], [1.0, 2.0, 0.5, 4.0]),
        ("Z", [], None),
    ]


VERB = {"M": 0, "L": 1, "Q": 2, "C": 3, "RQ": 4, "RC": 5, "Z": 6}


def path_bytes(dim):
    els = fixture(dim)
    out = u32(1) + u32(dim) + u64(len(els))
    for name, pts, w in els:
        out += u8(VERB[name])
        for p in pts:
            out += f64(p[0]) + f64(p[1]) + (f64(p[2]) if dim == 3 else b"")
        for x in (w or []):
            out += f64(x)
    return out


def dump_text(dim):
    els = fixture(dim)
    h = fnv(path_bytes(dim))
    lines = ["gshgeom.path 1 dim=%d elements=%d hash=%016x" % (dim, len(els), h)]
    for name, pts, w in els:
        s = name
        for x in (w or []):
            s += " %016x" % bits(x)
        for p in pts:
            s += " %016x %016x" % (bits(p[0]), bits(p[1]))
            if dim == 3:
                s += " %016x" % bits(p[2])
        lines.append(s)
    return "\n".join(lines) + "\n"


# --- maps --------------------------------------------------------------------
TR, SC, NU, FY, MI, QT, RO, SK, SG, GE = (1 << i for i in range(10))

# (name, dim, row-major matrix, translation, expected kind tag)
MAP_FIXTURES = [
    ("identity2", 2, [1, 0, 0, 1], [0, 0], 0),
    ("translate2", 2, [1, 0, 0, 1], [1.5, -2.25], TR),
    ("scale2", 2, [2, 0, 0, 2], [0, 0], SC),
    ("design2", 2, [2.5, 0, 0, -2.5], [10, 20], TR | SC | FY),
    ("quarter2", 2, [0, -1, 1, 0], [0, 0], QT),
    ("half2", 2, [-1, 0, 0, -1], [0, 0], QT),
    ("nonuniform2", 2, [2, 0, 0, 3], [0, 0], NU),
    ("mirrorx2", 2, [-1, 0, 0, 1], [0, 0], MI),
    ("rotate2", 2, [0.8, -0.6, 0.6, 0.8], [0.5, 0], RO | TR),
    ("mirrorrot2", 2, [0.6, 0.8, 0.8, -0.6], [0, 0], RO | MI),
    ("skew2", 2, [1, 0.5, 0, 1], [0, 0], SK),
    ("singular2", 2, [1, 2, 2, 4], [0, 0], SG),
    ("swap2", 2, [0, 1, 1, 0], [0, 0], QT | MI),
    ("cycle3", 3, [0, 1, 0, 0, 0, 1, 1, 0, 0], [0, 0, 0], QT),
    ("flip3", 3, [1, 0, 0, 0, -1, 0, 0, 0, 1], [1, 2, 3], FY | TR),
    ("general3", 3, [1, 2, 3, 0, 1, 4, 5, 6, 0], [0, 0, 0], GE),
    ("identity3", 3, [1, 0, 0, 0, 1, 0, 0, 0, 1], [0, 0, 0], 0),
]


def map_bytes(dim, m, t, kind):
    out = u32(1) + u32(dim) + u32(kind)
    for v in m:
        out += f64(float(v))
    for v in t:
        out += f64(float(v))
    return out


def map_dump(dim, m, t, kind):
    h = fnv(map_bytes(dim, m, t, kind))
    lines = ["gshgeom.map 1 dim=%d kind=%08x hash=%016x" % (dim, kind, h)]
    for i in range(dim):
        lines.append("m" + "".join(" %016x" % bits(float(m[i * dim + j])) for j in range(dim)))
    lines.append("t" + "".join(" %016x" % bits(float(t[i])) for i in range(dim)))
    return "\n".join(lines) + "\n"


def hexf(v):
    return float(v).hex()


# --- orientation -------------------------------------------------------------
def sgn(x):
    return (x > 0) - (x < 0)


def exact2(a, b, c):
    A, B, C = [[Fraction(v) for v in p] for p in (a, b, c)]
    return sgn((B[0] - A[0]) * (C[1] - A[1]) - (B[1] - A[1]) * (C[0] - A[0]))


def naive2(a, b, c):
    return sgn((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))


def exact3(a, b, c, d):
    A, B, C, D = [[Fraction(v) for v in p] for p in (a, b, c, d)]
    u = [B[i] - A[i] for i in range(3)]
    v = [C[i] - A[i] for i in range(3)]
    w = [D[i] - A[i] for i in range(3)]
    return sgn(u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0])
               + u[2] * (v[0] * w[1] - v[1] * w[0]))


def naive3(a, b, c, d):
    ux, uy, uz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
    vx, vy, vz = c[0] - a[0], c[1] - a[1], c[2] - a[2]
    wx, wy, wz = d[0] - a[0], d[1] - a[1], d[2] - a[2]
    return sgn(ux * (vy * wz - vz * wy) - uy * (vx * wz - vz * wx) + uz * (vx * wy - vy * wx))


def scaled(p, e):
    return tuple(math.ldexp(v, e) for v in p)


def orient2d_cases():
    rng = random.Random(20261018)
    cases = []
    u = 2.0 ** -53
    for i in range(16):
        for j in range(16):
            cases.append(((0.5 + i * u, 0.5 + j * u), (12.0, 12.0), (24.0, 24.0)))

    def grid():
        return tuple(rng.randrange(-(1 << 24), 1 << 24) / float(1 << 20) for _ in range(2))

    base = []
    for _ in range(300):
        a, b = grid(), grid()
        s = Fraction(rng.randrange(-16, 25), 8)
        c = [float(Fraction(a[k]) + s * (Fraction(b[k]) - Fraction(a[k]))) for k in range(2)]
        for _ in range(rng.randrange(0, 3)):
            k = rng.randrange(2)
            step = rng.choice((-1, 1))
            for _ in range(rng.randrange(1, 3)):
                c[k] = math.nextafter(c[k], math.inf * step)
        base.append((a, b, tuple(c)))
    cases += base
    for e in (-1040, -700, 600, 970):
        for case in base[:40]:
            cases.append(tuple(scaled(p, e) for p in case))
    h = 2.0 ** -1074
    cases += [
        ((0.0, 0.0), (h, 0.0), (0.0, h)),
        ((0.0, 0.0), (h, 0.0), (0.0, -h)),
        ((0.0, 0.0), (h, h), (2 * h, 2 * h)),
        ((0.0, 0.0), (1.0, 0.0), (0.0, 0.0)),
        ((1e300, 1e300), (2e300, 2e300), (3e300, 3e300)),
        ((1e300, 1e300), (2e300, 2e300), (3e300, math.nextafter(3e300, math.inf))),
        ((-0.0, 0.0), (0.0, -0.0), (0.0, 0.0)),
        ((1.7976931348623157e308, 0.0), (0.0, 1.7976931348623157e308), (-1.7976931348623157e308, 0.0)),
    ]
    return cases


def orient3d_cases():
    rng = random.Random(20261019)
    cases = []
    u = 2.0 ** -53
    for i in range(6):
        for j in range(6):
            for k in range(6):
                cases.append(((0.5 + i * u, 0.5 + j * u, 0.5 + k * u),
                              (24.0, 0.0, 0.0), (0.0, 24.0, 0.0), (0.0, 0.0, 24.0)))

    def grid():
        return tuple(rng.randrange(-(1 << 24), 1 << 24) / float(1 << 20) for _ in range(3))

    base = []
    for _ in range(300):
        a, b, c = grid(), grid(), grid()
        s = Fraction(rng.randrange(-16, 25), 8)
        t = Fraction(rng.randrange(-16, 25), 8)
        d = [float(Fraction(a[k]) + s * (Fraction(b[k]) - Fraction(a[k]))
                   + t * (Fraction(c[k]) - Fraction(a[k]))) for k in range(3)]
        for _ in range(rng.randrange(0, 3)):
            k = rng.randrange(3)
            step = rng.choice((-1, 1))
            d[k] = math.nextafter(d[k], math.inf * step)
        base.append((a, b, c, tuple(d)))
    cases += base
    for e in (-1040, -700, 600, 960):
        for case in base[:40]:
            cases.append(tuple(scaled(p, e) for p in case))
    # Three points nearly collinear in the plane (the Kettner family) and the
    # fourth far above it, with the other two points almost at the plane's
    # height: the determinant is one z-difference times a cancelling 2D
    # determinant, so one term of the filter's permanent dominates, and the
    # bound is wrong if that term is missing.  One family for each term.
    u = 2.0 ** -53
    for i in range(9):
        for j in range(9):
            p = (0.5 + i * u, 0.5 + j * u)
            q = (12.0, 12.0)
            r = (24.0, 24.0)
            far = (0.3, 7.7, 5.0)
            tiny = 2.0 ** -60
            # |adz| small, |bdz| small, |cdz| large: the third term.
            cases.append(((p[0], p[1], tiny), (q[0], q[1], tiny / 2), far, (r[0], r[1], 0.0)))
            # |adz| large: the first term (a above, b c d nearly collinear).
            cases.append((far, (p[0], p[1], tiny), (q[0], q[1], tiny / 2), (r[0], r[1], 0.0)))
            # |bdz| large: the second term (b above, a c d nearly collinear).
            cases.append(((p[0], p[1], tiny), far, (q[0], q[1], tiny / 2), (r[0], r[1], 0.0)))
    h = 2.0 ** -1074
    cases += [
        ((0.0, 0.0, 0.0), (h, 0.0, 0.0), (0.0, h, 0.0), (0.0, 0.0, h)),
        ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (1.0, 1.0, 0.0)),
        ((1e100, 0.0, 0.0), (0.0, 1e100, 0.0), (0.0, 0.0, 1e100), (1e100, 1e100, 1e100)),
    ]
    return cases


# --- story 7: the 26.6 rounding, the device path, the flattened path ----------
INT32_MIN, INT32_MAX = -(1 << 31), (1 << 31) - 1


def round_ref(x):
    """floor(x * 64 + 1/2) in exact rational arithmetic; None outside int32."""
    r = math.floor(Fraction(x) * 64 + Fraction(1, 2))
    return r if INT32_MIN <= r <= INT32_MAX else None


def round_naive(x):
    """The same expression evaluated in double: floor(fl(fl(x * 64) + 0.5))."""
    return math.floor(x * 64.0 + 0.5)


def round_cases():
    cases = []
    half = 0.5 - 2.0 ** -54
    ks = [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 100, -100, 1000, -1000, (1 << 20) + 1, -(1 << 20) - 1,
          (1 << 30) - 1, -(1 << 30), INT32_MAX - 1, INT32_MAX, INT32_MIN, INT32_MIN + 1,
          INT32_MAX + 1, INT32_MIN - 1]
    for k in ks:
        x = (k + 0.5) / 64.0
        cases += [math.nextafter(x, -math.inf), x, math.nextafter(x, math.inf)]
        y = k / 64.0
        cases += [math.nextafter(y, -math.inf), y, math.nextafter(y, math.inf)]
    cases += [half / 64.0, -half / 64.0, (1 + half) / 64.0, (-1 + half) / 64.0, (-1 - half) / 64.0,
              0.5 / 64.0, -0.5 / 64.0, 1.5 / 64.0, -1.5 / 64.0, 0.25 / 64.0, 0.75 / 64.0, -0.75 / 64.0]
    cases += [0.0, -0.0, 2.0 ** -1074, -(2.0 ** -1074), 2.0 ** -1022, 2.0 ** -200, 2.0 ** -8, -(2.0 ** -8),
              1.0, -1.0, 0.1, -0.1, 33554432.0, -33554432.0, 33554431.984375, -33554432.015625,
              1e10, -1e10, 1.7976931348623157e308, -1.7976931348623157e308, 2.0 ** 1000, -(2.0 ** 1000)]
    rng = random.Random(20261020)
    for _ in range(600):
        e = rng.randrange(-30, 27)
        x = math.ldexp(rng.randrange(1 << 52, 1 << 53), e - 52)
        cases.append(-x if rng.random() < 0.5 else x)
    for _ in range(150):
        # within a few units of a half-way point, in the last bits of the double
        k = rng.randrange(-(1 << 24), 1 << 24)
        x = (k + 0.5) / 64.0
        for _ in range(rng.randrange(0, 4)):
            x = math.nextafter(x, math.inf if rng.random() < 0.5 else -math.inf)
        cases.append(x)
    return cases


# Device path: the same sequence the test makes through the builder.  Verbs:
# M 0, L 1, Q 2, C 3, Z 4; a verb after a close stores a move first.
DEV_VERB = {"M": 0, "L": 1, "Q": 2, "C": 3, "Z": 4}
DEV_FIXTURE = [
    ("M", [(64, -128)]),
    ("L", [(192, 256)]),
    ("Q", [(320, 416), (448, 512)]),
    ("C", [(576, 640), (720, 768), (832, -896)]),
    ("Z", []),
    ("M", [(64, -128)]),  # the move the builder stores for the line below
    ("L", [(-2147483648, 2147483647)]),
    ("Z", []),
]


def dev_bytes(els):
    out = u32(1) + u32(0x50564544) + u64(len(els))
    for name, pts in els:
        out += u8(DEV_VERB[name])
        for x, y in pts:
            out += i32(x) + i32(y)
    return out


def hex32(v):
    return "%08x" % (v & 0xFFFFFFFF)


def dev_dump(els):
    lines = ["gshgeom.device-path 1 elements=%d hash=%016x" % (len(els), fnv(dev_bytes(els)))]
    for name, pts in els:
        lines.append(name + "".join(" %s %s" % (hex32(x), hex32(y)) for x, y in pts))
    return "\n".join(lines) + "\n"


# Flattened path: (closed, [points]) per subpath.
FLAT_FIXTURE = [
    (0, [(0, 0), (64, 0), (128, 0), (192, 3)]),
    (1, [(-1000, 1000), (1100, 1000), (-5, -7)]),
    (0, [(5, 5)]),
]


def flat_bytes(subs):
    out = u32(1) + u32(0x54414C46) + u64(len(subs))
    for closed, pts in subs:
        out += u8(closed) + u64(len(pts))
        for x, y in pts:
            out += i32(x) + i32(y)
    return out


def flat_dump(subs):
    npts = sum(len(p) for _, p in subs)
    lines = ["gshgeom.flat-path 1 subpaths=%d points=%d hash=%016x" % (len(subs), npts, fnv(flat_bytes(subs)))]
    for closed, pts in subs:
        lines.append("S closed=%d points=%d" % (closed, len(pts)))
        for x, y in pts:
            lines.append("P %s %s" % (hex32(x), hex32(y)))
    return "\n".join(lines) + "\n"


# A model path under a map to the device: all numbers are dyadic, so the map's
# arithmetic is exact in double and the reference is exact rational arithmetic.
MODEL_DEV_PATH = [
    ("M", [(0.0078125, -0.0078125)]),
    ("L", [(1.5, -2.25)]),
    ("Q", [(3.0, 4.5), (5.0, 1.0)]),
    ("C", [(6.0, 7.0), (8.25, 9.0), (10.0, 11.5)]),
    ("Z", []),
    ("L", [(2.5, 3.5)]),
]
MODEL_DEV_MAP = ([2.0, 0.0, 0.0, 4.0], [0.2421875, 0.0234375])


def model_to_device():
    m, t = MODEL_DEV_MAP
    els = []
    closed = False
    start = None
    for name, pts in MODEL_DEV_PATH:
        dev = []
        for x, y in pts:
            fx = Fraction(m[0]) * Fraction(x) + Fraction(m[1]) * Fraction(y) + Fraction(t[0])
            fy = Fraction(m[2]) * Fraction(x) + Fraction(m[3]) * Fraction(y) + Fraction(t[1])
            rx, ry = round_ref(float(fx)), round_ref(float(fy))
            assert float(fx) == fx and float(fy) == fy and rx is not None and ry is not None
            dev.append((rx, ry))
        if name == "M":
            start = dev[0]
        if closed and name in ("L", "Q", "C"):
            els.append(("M", [start]))
        closed = name == "Z"
        els.append((name, dev))
    return els


# --- story 8: the measured path of lines (lengths are exact integers) ---------

# (verb, points) in element order; every segment is a straight line whose
# length is an exact integer (3-4-5 and 5-12-13 triangles), so the whole table
# is computable here without any quadrature.
MEASURED_FIXTURES = {
    2: [("M", [(0.0, 0.0)]), ("L", [(3.0, 4.0)]),
        ("M", [(10.0, 10.0)]), ("L", [(13.0, 10.0)]), ("L", [(13.0, 10.0)]),
        ("L", [(13.0, 14.0)]), ("Z", []), ("L", [(10.0, 10.0)]), ("L", [(10.0, 16.0)])],
    3: [("M", [(1.0, 2.0, 3.0)]), ("L", [(4.0, 6.0, 15.0)]),
        ("M", [(0.0, 0.0, 0.0)]), ("L", [(0.0, 0.0, 0.0)]), ("L", [(2.0, 3.0, 6.0)]), ("Z", [])],
}


def meas_build(dim, els):
    """Segments, rows and subpaths by the rules in along.h, from scratch."""
    import math
    segs, rows, subs = [], [], []
    cur = None
    start = None
    total = 0.0
    need_sub = True
    origin = None
    for i, (verb, pts) in enumerate(els):
        if verb == "M":
            if origin is None:
                origin = pts[0]
            cur = start = pts[0]
            need_sub = True
            continue
        end = start if verb == "Z" else pts[0]
        if need_sub:
            subs.append([total, total, 0])
            need_sub = False
        d2 = sum((e - c) ** 2 for c, e in zip(cur, end))
        length = float(math.isqrt(int(d2)))
        assert length * length == d2, "not an exact length"
        if length != 0.0:
            segs.append((i, cur, end))
            rows.append((len(segs) - 1, 0.0, 1.0, total, length))
            total += length
        subs[-1][1] = total
        if verb == "Z":
            subs[-1][2] = 1
            need_sub = True
        cur = end
    return origin, segs, rows, subs, total


def meas_bytes(dim, els):
    origin, segs, rows, subs, total = meas_build(dim, els)
    b = u32(1) + u32(0x5341454D) + u32(dim) + u64(len(segs)) + u64(len(rows)) + u64(len(subs))
    for c in origin[:dim]:
        b += f64(c)
    b += f64(total)
    for el, p0, p1 in segs:
        b += u64(el) + u8(0)
        for c in p0[:dim] + p1[:dim]:
            b += f64(c)
    for seg, t0, t1, s0, ln in rows:
        b += u64(seg) + f64(t0) + f64(t1) + f64(s0) + f64(ln)
    for begin, end, closed in subs:
        b += f64(begin) + f64(end) + u8(closed)
    return b


def meas_dump(dim, els):
    origin, segs, rows, subs, total = meas_build(dim, els)
    h = fnv(meas_bytes(dim, els))
    out = ["gshgeom.measured 1 dim=%d segments=%d rows=%d subpaths=%d length=%016x hash=%016x" % (
        dim, len(segs), len(rows), len(subs), bits(total), h)]
    out.append("O " + " ".join("%016x" % bits(c) for c in origin[:dim]))
    for el, p0, p1 in segs:
        out.append("G %d 0 " % el + " ".join("%016x" % bits(c) for c in p0[:dim] + p1[:dim]))
    for seg, t0, t1, s0, ln in rows:
        out.append("R %d " % seg + " ".join("%016x" % bits(c) for c in (t0, t1, s0, ln)))
    for begin, end, closed in subs:
        out.append("S %016x %016x %d" % (bits(begin), bits(end), closed))
    return "\n".join(out) + "\n"


# --- story 7a: the flattener, revision 2 -------------------------------------
# A second statement of the algorithm in src/flatten/flatten.c, in Python's
# unbounded integers and sharing no code with it.  `python3 tools/gen-hash-vectors.py
# --check-flatten tests/flatten_vectors.h` runs it over every recorded vector
# (recorded from font's flattener, tools/vectors-from-font) and fails on the
# first difference in count, hash or explicit vertex list.
FLATTEN_REVISION = 2
FLATTEN_TOLERANCE = 4
FLATTEN_DEPTH = 20


def _wmid(a, b):
    return ((a[0] + b[0] + 1) // 2, (a[1] + b[1] + 1) // 2)


def _narrow(p):
    return ((p[0] + 128) // 256, (p[1] + 128) // 256)


def flatten_model(ctrl, tol=FLATTEN_TOLERANCE, depth=FLATTEN_DEPTH):
    """The vertices (not the start point) of a quadratic (3 points) or cubic (4)."""
    assert tol >= 2
    limit = 4 * 256 * (tol - 1)
    out = []

    def rec(P, level):
        if len(P) == 3:
            p0, p1, p2 = P
            flat = abs(2 * p1[0] - p0[0] - p2[0]) + abs(2 * p1[1] - p0[1] - p2[1]) <= limit
        else:
            p0, p1, p2, p3 = P
            d1 = abs(3 * p1[0] - 2 * p0[0] - p3[0]) + abs(3 * p1[1] - 2 * p0[1] - p3[1])
            d2 = abs(3 * p2[0] - p0[0] - 2 * p3[0]) + abs(3 * p2[1] - p0[1] - 2 * p3[1])
            flat = d1 <= limit and d2 <= limit
        if flat or level >= depth:
            assert flat, "the depth cap forced a piece"
            out.append(_narrow(P[-1]))
            return
        if len(P) == 3:
            a, b = _wmid(p0, p1), _wmid(p1, p2)
            m = _wmid(a, b)
            rec((p0, a, m), level + 1)
            rec((m, b, p2), level + 1)
        else:
            a, b, c = _wmid(p0, p1), _wmid(p1, p2), _wmid(p2, p3)
            d, e = _wmid(a, b), _wmid(b, c)
            m = _wmid(d, e)
            rec((p0, a, d, m), level + 1)
            rec((m, e, c, p3), level + 1)

    rec(tuple((x * 256, y * 256) for x, y in ctrl), 0)
    return out


def check_flatten_vectors(path):
    import re
    text = open(path).read()
    assert re.search(r"kFlattenVectorsRevision = (\d+);", text).group(1) == str(FLATTEN_REVISION)
    pts_body = re.search(r"kFlattenPoints\[\] = \{(.*?)\};", text, re.S).group(1)
    pts = [int(v) for v in re.findall(r"-?\d+", pts_body)]
    vec_body = re.search(r"kFlattenVectors\[\] = \{(.*?)\n\};", text, re.S).group(1)
    rows = re.findall(
        r"\{(\d), \{\{(-?\d+), (-?\d+)\}, \{(-?\d+), (-?\d+)\}, \{(-?\d+), (-?\d+)\}, \{(-?\d+), (-?\d+)\}\}, "
        r"(\d+), UINT64_C\(0x([0-9a-f]+)\), (\d+)u?\}", vec_body)
    n = int(re.search(r"kFlattenVectorCount = (\d+);", text).group(1))
    assert len(rows) == n, (len(rows), n)
    explicit = 0
    for i, r in enumerate(rows):
        kind = int(r[0])
        v = [int(x) for x in r[1:9]]
        count, h, off = int(r[9]), int(r[10], 16), int(r[11])
        ctrl = [(v[2 * k], v[2 * k + 1]) for k in range(kind + 1)]
        got = flatten_model(ctrl)
        assert len(got) == count, "vector %d: %d vertices, recorded %d" % (i, len(got), count)
        assert fnv(b"".join(i32(x) + i32(y) for x, y in got)) == h, "vector %d: hash" % i
        assert got[-1] == ctrl[-1], "vector %d: the last vertex is the curve's end" % i
        if off != 0xFFFFFFFF:
            explicit += 1
            want = [(pts[off + 2 * k], pts[off + 2 * k + 1]) for k in range(count)]
            assert got == want, "vector %d: vertex list" % i
    print("flatten model (revision %d): %d vectors agree, %d with explicit lists" % (FLATTEN_REVISION, n, explicit))


def main():
    if len(sys.argv) == 3 and sys.argv[1] == "--check-flatten":
        check_flatten_vectors(sys.argv[2])
        return
    o = sys.stdout
    o.write("// Generated by tools/gen-hash-vectors.py.  Do not edit.\n")
    o.write("// Known answers of the hash spelling, from an implementation that\n")
    o.write("// shares no code with the library.\n\n")
    o.write("#ifndef GHOTI_IO_GSHGEOM_KNOWN_ANSWERS_H\n#define GHOTI_IO_GSHGEOM_KNOWN_ANSWERS_H\n\n")
    o.write("#include <cstdint>\n\n")

    def c(name, value):
        o.write("static const std::uint64_t %s = UINT64_C(0x%016x);\n" % (name, value))

    c("kFnvEmpty", fnv(b""))
    c("kFnvA", fnv(b"a"))
    c("kFnvFoobar", fnv(b"foobar"))
    c("kHashU8", fnv(u8(0xA5)))
    c("kHashU16", fnv(u16(0x1234)))
    c("kHashU32", fnv(u32(0x01020304)))
    c("kHashU64", fnv(u64(0x0102030405060708)))
    c("kHashI32", fnv(i32(-2)))
    c("kHashI64", fnv(i64(-3)))
    c("kHashF64One", fnv(f64(1.0)))
    c("kHashF64NegZero", fnv(f64(-0.0)))
    c("kHashF64Tenth", fnv(f64(0.1)))
    c("kHashString", fnv(string("abc")))
    c("kHashEmptyString", fnv(string("")))
    c("kHashMixed", fnv(u8(1) + u16(2) + u32(3) + u64(4) + f64(5.5) + string("xy") + b"zw"))
    c("kPathHash2", fnv(path_bytes(2)))
    c("kPathHash3", fnv(path_bytes(3)))
    c("kPathHashEmpty2", fnv(u32(1) + u32(2) + u64(0)))
    c("kPathHashEmpty3", fnv(u32(1) + u32(3) + u64(0)))
    o.write("\n")
    for dim in (2, 3):
        o.write("static const char kPathDump%d[] =\n" % dim)
        for line in dump_text(dim).splitlines():
            o.write('    "%s\\n"\n' % line)
        o.write("    ;\n")

    o.write("struct MapFixture {\n  const char * name;\n  int dim;\n  double m[9];\n  double t[3];\n  unsigned kind;\n  std::uint64_t hash;\n  const char * dump;\n};\n")
    o.write("static const MapFixture kMapFixtures[] = {\n")
    for name, dim, m, t, kind in MAP_FIXTURES:
        mm = [hexf(v) for v in m] + ["0.0"] * (9 - len(m))
        tt = [hexf(v) for v in t] + ["0.0"] * (3 - len(t))
        o.write('    {"%s", %d, {%s}, {%s}, %du, UINT64_C(0x%016x),\n' % (
            name, dim, ", ".join(mm), ", ".join(tt), kind, fnv(map_bytes(dim, m, t, kind))))
        for line in map_dump(dim, m, t, kind).splitlines():
            o.write('     "%s\\n"\n' % line)
        o.write("    },\n")
    o.write("};\n\n")

    c2 = orient2d_cases()
    c3 = orient3d_cases()
    o.write("struct Orient2dCase {\n  double a[2], b[2], c[2];\n  int sign;\n};\n")
    o.write("struct Orient3dCase {\n  double a[3], b[3], c[3], d[3];\n  int sign;\n};\n")
    wrong2 = sum(1 for a, b, c in c2 if naive2(a, b, c) != exact2(a, b, c))
    wrong3 = sum(1 for a, b, c, d in c3 if naive3(a, b, c, d) != exact3(a, b, c, d))
    o.write("static const unsigned kNaiveWrong2d = %d;\n" % wrong2)
    o.write("static const unsigned kNaiveWrong3d = %d;\n" % wrong3)
    o.write("static const Orient2dCase kOrient2dCases[] = {\n")
    for a, b, c in c2:
        o.write("    {{%s}, {%s}, {%s}, %d},\n" % tuple(
            [", ".join(hexf(v) for v in p) for p in (a, b, c)] + [exact2(a, b, c)]))
    o.write("};\n")
    o.write("static const Orient3dCase kOrient3dCases[] = {\n")
    for a, b, c, d in c3:
        o.write("    {{%s}, {%s}, {%s}, {%s}, %d},\n" % tuple(
            [", ".join(hexf(v) for v in p) for p in (a, b, c, d)] + [exact3(a, b, c, d)]))
    o.write("};\n")

    # --- story 7 ---
    cases = round_cases()
    o.write("\nstruct RoundCase {\n  double x;\n  int limit;  // 1 when the result does not fit an int32\n  std::int32_t value;\n};\n")
    o.write("static const unsigned kNaiveRoundWrong = %d;\n" % sum(
        1 for x in cases if round_ref(x) is not None and round_naive(x) != round_ref(x)))
    o.write("static const RoundCase kRoundCases[] = {\n")
    for x in cases:
        r = round_ref(x)
        o.write("    {%s, %d, %d},\n" % (hexf(x), 1 if r is None else 0, 0 if r is None else r))
    o.write("};\n")
    o.write("static const double kNaiveRoundWrongInput = %s;\n" % hexf((0.5 - 2.0 ** -54) / 64.0))
    assert round_naive((0.5 - 2.0 ** -54) / 64.0) == 1 and round_ref((0.5 - 2.0 ** -54) / 64.0) == 0
    def c64(name, value):
        o.write("static const std::uint64_t %s = UINT64_C(0x%016x);\n" % (name, value))

    c64("kDevicePathHash", fnv(dev_bytes(DEV_FIXTURE)))
    c64("kFlatPathHash", fnv(flat_bytes(FLAT_FIXTURE)))
    md = model_to_device()
    c64("kModelDeviceHash", fnv(dev_bytes(md)))
    for name, text in (("kDevicePathDump", dev_dump(DEV_FIXTURE)), ("kFlatPathDump", flat_dump(FLAT_FIXTURE)),
                       ("kModelDeviceDump", dev_dump(md))):
        o.write("static const char %s[] =\n" % name)
        for line in text.splitlines():
            o.write('    "%s\\n"\n' % line)
        o.write("    ;\n")
    # --- story 8 ---
    for dim in (2, 3):
        c64("kMeasuredHash%d" % dim, fnv(meas_bytes(dim, MEASURED_FIXTURES[dim])))
        o.write("static const char kMeasuredDump%d[] =\n" % dim)
        for line in meas_dump(dim, MEASURED_FIXTURES[dim]).splitlines():
            o.write('    "%s\\n"\n' % line)
        o.write("    ;\n")
    o.write("\n#endif\n")


main()
