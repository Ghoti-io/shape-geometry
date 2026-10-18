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
"""Print the reference vectors of along-path parametrization (CAP-5, AD-1).

An independent high-precision statement of the arc length of a curve and of the
point and tangent at a distance along it, for the curves listed in CURVES.  It
uses mpmath (MPFR-grade, 50 digits here) and none of the library's code: the
curve is expanded into power-basis polynomials, the velocity is the exact
derivative of that expansion (checked against numerical differentiation of the
curve itself), the length is mpmath's tanh-sinh/Gauss-Legendre quadrature with
the integrand split at the speed's known kinks, and the position at a distance
is found by root finding on the arc-length function.

It also regenerates the Gauss-Legendre nodes and weights (400 bits) and checks
that the literals in src/along/along.c are the nearest doubles.

Run it in the pinned oracle image (mpmath is not installed on the host):

    podman run --rm --memory=3g --memory-swap=3g --cpus=2 \\
        -v "$PWD:/w:Z" -w /w localhost/ghoti-math-oracle-mpfr:2.3.2 \\
        python3 tools/gen-along-vectors.py > tests/along_vectors.h

The image is the one libs/math/tools/oracle builds (gmpy2 2.3.2, mpmath 1.4.1).
Every number is written as a hexadecimal floating literal, so nothing is lost
in a decimal parse.  The inputs are exact doubles.
"""
import re
import sys
from pathlib import Path

import mpmath
from mpmath import mp, mpf

mp.dps = 50

ROOT = Path(__file__).resolve().parent.parent

# kind: 0 line 1 quadratic 2 cubic 3 rational quadratic 4 rational cubic.
# points: control points (3 coordinates, z = 0 in 2D); weights or None;
# kinks: parameters at which the speed is not smooth.
CURVES = [
    ("quadratic", 1, 2, [(0, 0), (2, 3), (5, -1)], None, []),
    ("quadratic_3d", 1, 3, [(0, 0, 0), (2, 3, 1), (5, -1, 4)], None, []),
    ("cubic_s", 2, 2, [(0, 0), (1, 3), (2, -3), (3, 0)], None, []),
    ("cubic_loop", 2, 2, [(0, 0), (3, 3), (-2, 3), (1, 0)], None, []),
    # B'(1/2) = 0 and B''(1/2) != 0: a cusp at t = 1/2
    ("cubic_cusp", 2, 2, [(0, 0), (2, 1), (1, 1), (1, 0)], None, [mpf(1) / 2]),
    # B'(1/3) = 0: the cusp is not at a panel boundary of the first split
    ("cubic_cusp_third", 2, 2, [(0, 0), (2, 1), (3, 2), (-9, -6)], None, [mpf(1) / 3]),
    ("cubic_nearline", 2, 2, [(0, 0), (1, 2.0**-10), (2, -2.0**-10), (3, 0)], None, []),
    ("cubic_3d", 2, 3, [(0, 0, 0), (1, 2, 3), (-1, 1, 2), (2, -1, 4)], None, []),
    ("cubic_large", 2, 2, [(0, 0), (1e6, 2e6), (-5e5, 1e6), (3e5, -1e5)], None, []),
    ("cubic_p0_is_p1", 2, 2, [(0, 0), (0, 0), (3, 1), (4, 4)], None, []),
    ("rational_quadratic", 3, 2, [(0, 0), (1, 2), (3, 0)], [1, 3, 1], []),
    ("rational_quadratic_light", 3, 2, [(0, 0), (1, 2), (3, 0)], [1, 0.125, 2], []),
    ("rational_cubic", 4, 2, [(0, 0), (1, 2), (3, 2), (4, 0)], [1, 2, 0.5, 1.5], []),
    ("rational_cubic_3d", 4, 3, [(0, 0, 0), (1, 2, 1), (3, 2, -1), (4, 0, 2)], [1, 0.75, 3, 1], []),
    ("rational_cubic_heavy", 4, 2, [(0, 0), (1, 2), (3, 2), (4, 0)], [1, 50, 0.0625, 1], []),
]

FRACTIONS = [0.1, 0.25, 0.5, 0.9]


def hexf(v):
    return float(v).hex()


def poly_mul(a, b):
    out = [mpf(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def poly_add(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]


def poly_scale(a, k):
    return [k * x for x in a]


def binom(n, k):
    return mpmath.binomial(n, k)


def bernstein_poly(n, i):
    # C(n,i) t^i (1-t)^(n-i) in the power basis
    p = [mpf(1)]
    for _ in range(i):
        p = poly_mul(p, [mpf(0), mpf(1)])
    for _ in range(n - i):
        p = poly_mul(p, [mpf(1), mpf(-1)])
    return poly_scale(p, binom(n, i))


def poly_eval(a, t):
    r = mpf(0)
    for c in reversed(a):
        r = r * t + c
    return r


def poly_diff(a):
    return [i * a[i] for i in range(1, len(a))] or [mpf(0)]


class Curve:
    def __init__(self, name, kind, dim, pts, w, kinks):
        self.name, self.kind, self.dim, self.kinks = name, kind, dim, kinks
        self.pts = [tuple(mpf(float(c)) for c in (list(p) + [0] * (3 - len(p)))) for p in pts]
        self.w = [mpf(float(x)) for x in w] if w else None
        n = len(pts) - 1
        self.n = n
        wt = self.w or [mpf(1)] * (n + 1)
        self.N = []
        for a in range(3):
            poly = [mpf(0)]
            for i in range(n + 1):
                poly = poly_add(poly, poly_scale(bernstein_poly(n, i), wt[i] * self.pts[i][a]))
            self.N.append(poly)
        self.D = [mpf(0)]
        for i in range(n + 1):
            self.D = poly_add(self.D, poly_scale(bernstein_poly(n, i), wt[i]))
        self.dN = [poly_diff(p) for p in self.N]
        self.dD = poly_diff(self.D)

    def point(self, t):
        d = poly_eval(self.D, t)
        return [poly_eval(self.N[a], t) / d for a in range(3)]

    def e_vec(self, t):
        d, dd = poly_eval(self.D, t), poly_eval(self.dD, t)
        return [poly_eval(self.dN[a], t) * d - poly_eval(self.N[a], t) * dd for a in range(3)]

    def velocity(self, t):
        d = poly_eval(self.D, t)
        return [x / (d * d) for x in self.e_vec(t)]

    def speed(self, t):
        v = self.velocity(t)
        return mpmath.sqrt(sum(x * x for x in v))

    def length(self, a, b, method="tanh-sinh"):
        # the kinks and a grid of 1/32, so that a peaked integrand (a heavy
        # weight) is resolved by either method
        grid = [mpf(i) / 32 for i in range(1, 32)]
        pts = sorted(set([a] + [k for k in self.kinks if a < k < b] + [g for g in grid if a < g < b] + [b]))
        return mpmath.quad(self.speed, pts, method=method, maxdegree=10)

    def polygon(self):
        total = mpf(0)
        for i in range(self.n):
            total += mpmath.sqrt(sum((self.pts[i + 1][a] - self.pts[i][a]) ** 2 for a in range(3)))
        return total

    def at_distance(self, s):
        f = lambda u: self.length(mpf(0), u) - s
        u = mpmath.findroot(f, (mpf("0.001"), mpf("0.999")), solver="illinois", tol=mpf(10) ** -40, maxsteps=200)
        return u

    def tangent(self, t):
        e = self.e_vec(t)
        m = mpmath.sqrt(sum(x * x for x in e))
        return [x / m for x in e]


def check_velocity(c):
    # the exact derivative of the expansion against numerical differentiation of the curve
    for t in (mpf("0.13"), mpf("0.5") + mpf("0.07"), mpf("0.91")):
        for a in range(3):
            num = mpmath.diff(lambda x: c.point(x)[a], t)
            assert abs(num - c.velocity(t)[a]) < mpf(10) ** -30 * (1 + abs(num)), (c.name, a)


def gl_check():
    src = (ROOT / "src/along/along.c").read_text()
    nodes = [float.fromhex(x) for x in re.findall(r"0x1\.[0-9a-f]+p-\d", src.split("gl_node")[1].split("};")[0])]
    weights = [float.fromhex(x) for x in re.findall(r"0x1\.[0-9a-f]+p-\d", src.split("gl_weight")[1].split("};")[0])]
    mp.prec = 400
    n = 8

    def legendre(x, m):
        p0, p1 = mpf(1), x
        for k in range(2, m + 1):
            p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
        return p1

    want = []
    for i in range(1, n // 2 + 1):
        x = mpmath.cos(mpmath.pi * (i - mpf("0.25")) / (n + mpf("0.5")))
        for _ in range(100):
            p, pm = legendre(x, n), legendre(x, n - 1)
            dp = n * (x * p - pm) / (x * x - 1)
            x = x - p / dp
        p, pm = legendre(x, n), legendre(x, n - 1)
        dp = n * (x * p - pm) / (x * x - 1)
        want.append((float(x), float(2 / ((1 - x * x) * dp * dp))))
    want.sort()
    mp.dps = 50
    assert nodes == [x for x, _ in want], (nodes, want)
    assert weights == [w for _, w in want], (weights, want)


def main():
    gl_check()
    out = []
    out.append("/* Generated by tools/gen-along-vectors.py (mpmath, 50 digits). Do not edit. */")
    out.append("#ifndef GHOTI_IO_GSHGEOM_TEST_ALONG_VECTORS_H")
    out.append("#define GHOTI_IO_GSHGEOM_TEST_ALONG_VECTORS_H")
    out.append("")
    out.append("struct AlongRef {")
    out.append("  const char * name;")
    out.append("  int kind;        // GSHGEOM_Curve_Kind")
    out.append("  int dim;")
    out.append("  double p[4][3];")
    out.append("  double w[4];")
    out.append("  double length;   // the exact length, rounded to nearest")
    out.append("  double polygon;  // the length of the control polygon, rounded")
    out.append("  double cusp_s;   // the distance of the cusp, or -1")
    out.append("  double s[4];     // distances: FRACTIONS of the length, rounded to double")
    out.append("  double pt[4][3]; // the exact point at the double s")
    out.append("  double tan[4][3];// the exact unit tangent there")
    out.append("};")
    out.append("")
    out.append("static const AlongRef kAlongRefs[] = {")
    for spec in CURVES:
        c = Curve(*spec)
        check_velocity(c)
        L = c.length(mpf(0), mpf(1))
        # a second quadrature method must agree
        L2 = c.length(mpf(0), mpf(1), method="gauss-legendre")
        assert abs(L - L2) < mpf(10) ** -30 * L, (c.name, L, L2)
        Lf = float(L)
        row = []
        row.append('    "%s", %d, %d,' % (c.name, c.kind, c.dim))
        pts = ", ".join("{%s}" % ", ".join(hexf(x) for x in p) for p in c.pts + [(0, 0, 0)] * (4 - len(c.pts)))
        row.append("    {%s}," % pts)
        ws = [hexf(x) for x in c.w] if c.w else []
        ws += ["0.0"] * (4 - len(ws))
        row.append("    {%s}," % ", ".join(ws))
        row.append("    %s, %s," % (hexf(L), hexf(c.polygon())))
        cusp = hexf(c.length(mpf(0), c.kinks[0])) if c.kinks else "-1.0"
        row.append("    %s," % cusp)
        ss, pp, tt = [], [], []
        for f in FRACTIONS:
            sd = float(f) * Lf
            u = c.at_distance(mpf(sd))
            ss.append(hexf(sd))
            pp.append("{%s}" % ", ".join(hexf(x) for x in c.point(u)))
            tt.append("{%s}" % ", ".join(hexf(x) for x in c.tangent(u)))
        row.append("    {%s}," % ", ".join(ss))
        row.append("    {%s}," % ", ".join(pp))
        row.append("    {%s}" % ", ".join(tt))
        out.append("  {\n" + "\n".join(row) + "\n  },")
    out.append("};")
    out.append("")
    out.append("#endif")
    print("\n".join(out))


if __name__ == "__main__":
    main()
