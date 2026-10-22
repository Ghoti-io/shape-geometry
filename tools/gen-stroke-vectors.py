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
"""Print the exact verdicts of stroke membership (SG5 to SG13) for tests/stroke_vectors.h.

An independent statement of what a stroke covers, in mpmath at 30 digits, sharing
no code with the library or with the C++ oracle (tests/stroke_oracle.h).  For each
case below (a centerline, a style and a map) it decides points inside or outside
the region by the rules as written in stroke.h: a point is in a segment's body
when it lies on a normal segment of length h either side (the parameters where
(x - c(t)) . c'(t) = 0, found by bracketing and bisection), in a cusp disc, or in
a join, cap or degenerate shape.  The miter limit (SG9) is decided in exact
rational arithmetic.  The arcs join follows the SVG 2 Candidate Recommendation of
2018-10-04 (sections 13.5.4 to 13.5.9), with the gaps settled as SG8 states, and
its constructions (the meeting point X, the clip point Y and the clip's direction)
are printed too.

A point is kept only when it is decided: every point at 0.2 px from it in sixteen
directions has the same verdict, so the expansion's 1/16 px tolerance cannot
reach it.  Explicit probes name the point each case exists for (the miter tip,
the cusp disc, the far side of a degenerate square); the rest come from a grid.

Run it in the pinned oracle image (mpmath is not installed on the host):

    podman run --rm --memory=3g --memory-swap=3g --cpus=2 \\
        -v "$PWD:/w:Z" -w /w localhost/ghoti-math-oracle-mpfr:2.3.2 \\
        python3 tools/gen-stroke-vectors.py > tests/stroke_vectors.h

The image is the one libs/math/tools/oracle builds (gmpy2 2.3.2, mpmath 1.4.1).
Every input is an exact double; every number printed is a hexadecimal floating
literal.
"""
import math
import sys
from fractions import Fraction

import mpmath
from mpmath import mp, mpf

mp.dps = 30

JOIN = {"miter": 0, "round": 1, "bevel": 2, "miter-clip": 3, "arcs": 4}
CAP = {"butt": 0, "round": 1, "square": 2}


def V(x, y):
    return (mpf(x), mpf(y))


def add(a, b):
    return (a[0] + b[0], a[1] + b[1])


def sub(a, b):
    return (a[0] - b[0], a[1] - b[1])


def mul(s, a):
    return (s * a[0], s * a[1])


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1]


def crs(a, b):
    return a[0] * b[1] - a[1] * b[0]


def nrm(a):
    return mpmath.sqrt(dot(a, a))


def unit(a):
    n = nrm(a)
    return (a[0] / n, a[1] / n)


def left(t):
    # orient.h: the tangent turned by a negative quarter turn
    return (t[1], -t[0])


def right(t):
    return (-t[1], t[0])


# --- the centerline ---------------------------------------------------------------


class Seg:
    def __init__(self, pts, w=None):
        self.n = len(pts) - 1
        self.p = [V(*q) for q in pts]
        self.w = [mpf(x) for x in (w or [1] * len(pts))]
        self.rational = w is not None
        self.raw = [tuple(q) for q in pts]

    def extent(self):
        return any(q != self.raw[0] for q in self.raw[1:])

    def start_dir(self):
        for k in range(1, self.n + 1):
            if self.raw[k] != self.raw[0]:
                return sub(self.p[k], self.p[0]), k
        return None, None

    def end_dir(self):
        for k in range(self.n - 1, -1, -1):
            if self.raw[k] != self.raw[self.n]:
                return sub(self.p[self.n], self.p[k]), k
        return None, None

    def homog(self, t):
        n = self.n
        num = [mul(self.w[i], self.p[i]) for i in range(n + 1)]
        den = list(self.w)
        for r in range(1, n + 1):
            for i in range(n + 1 - r):
                num[i] = add(mul(1 - t, num[i]), mul(t, num[i + 1]))
                den[i] = (1 - t) * den[i] + t * den[i + 1]
        return num[0], den[0]

    def point(self, t):
        num, den = self.homog(t)
        return (num[0] / den, num[1] / den)

    def deriv(self, t):
        """The direction of c'(t): N'D - ND' by finite exact derivatives of the Bernstein form."""
        n = self.n
        num = [mul(self.w[i], self.p[i]) for i in range(n + 1)]
        den = list(self.w)

        def bern(c, t):
            c = list(c)
            for r in range(1, len(c)):
                for i in range(len(c) - r):
                    c[i] = (1 - t) * c[i] + t * c[i + 1]
            return c[0]

        dnum = [mul(n, sub(num[i + 1], num[i])) for i in range(n)]
        dden = [n * (den[i + 1] - den[i]) for i in range(n)]
        N = (bern([q[0] for q in num], t), bern([q[1] for q in num], t))
        D = bern(den, t)
        dN = (bern([q[0] for q in dnum], t), bern([q[1] for q in dnum], t))
        dD = bern(dden, t)
        return (dN[0] * D - N[0] * dD, dN[1] * D - N[1] * dD)


def read(ops):
    """SG5: subpaths in order, each with its segments with extent, closed flag, whether it draws, its start."""
    subs = []
    cur = None
    is_open = False
    for verb, a in ops:
        if verb == "M":
            subs.append({"segs": [], "closed": False, "draws": False, "start": tuple(a)})
            cur = tuple(a)
            is_open = True
            continue
        if not is_open:
            subs.append({"segs": [], "closed": False, "draws": False, "start": subs[-1]["start"]})
            is_open = True
        sb = subs[-1]
        sb["draws"] = True
        if verb == "L":
            s = Seg([cur, (a[0], a[1])])
            cur = (a[0], a[1])
        elif verb == "Q":
            s = Seg([cur, (a[0], a[1]), (a[2], a[3])])
            cur = (a[2], a[3])
        elif verb == "C":
            s = Seg([cur, (a[0], a[1]), (a[2], a[3]), (a[4], a[5])])
            cur = (a[4], a[5])
        elif verb == "R":
            s = Seg([cur, (a[0], a[1]), (a[2], a[3])], [a[4], a[5], a[6]])
            cur = (a[2], a[3])
        elif verb == "K":
            s = Seg([cur, (a[0], a[1]), (a[2], a[3]), (a[4], a[5])], [a[6], a[7], a[8], a[9]])
            cur = (a[4], a[5])
        elif verb == "Z":
            s = Seg([cur, sb["start"]])
            cur = sb["start"]
            sb["closed"] = True
            is_open = False
        if s.extent():
            sb["segs"].append(s)
    return subs


def mapped(seg, M):
    pts = []
    for q in seg.raw:
        x = M[0] * q[0] + M[1] * q[1] + M[4]
        y = M[2] * q[0] + M[3] * q[1] + M[5]
        pts.append((x, y))
    return Seg(pts, [float(w) for w in seg.w] if seg.rational else None)


# --- SG9, exactly ---------------------------------------------------------------------


def miter_within(Q, Vx, S, L):
    """1 / sin(theta/2) <= L, as (L^2 - 2) |a||b| >= -L^2 (a . b), in rationals."""
    Q = [Fraction(c) for c in Q]
    Vx = [Fraction(c) for c in Vx]
    S = [Fraction(c) for c in S]
    L = Fraction(L)
    a = (Vx[0] - Q[0], Vx[1] - Q[1])
    b = (S[0] - Vx[0], S[1] - Vx[1])
    ab = a[0] * b[0] + a[1] * b[1]
    na2 = a[0] ** 2 + a[1] ** 2
    nb2 = b[0] ** 2 + b[1] ** 2
    P = L * L - 2
    rhs = -L * L * ab
    # compare P sqrt(na2 nb2) with rhs
    if P >= 0 and rhs <= 0:
        return True
    if P <= 0 and rhs > 0:
        return False
    if P < 0 and rhs <= 0:
        # both <= 0: P r >= rhs  <=>  |P| r <= |rhs|
        return P * P * na2 * nb2 <= rhs * rhs
    # P > 0, rhs > 0
    return P * P * na2 * nb2 >= rhs * rhs


# --- shapes ---------------------------------------------------------------------------


def in_convex(poly, x):
    pos = neg = 0
    n = len(poly)
    for i in range(n):
        c = crs(sub(poly[(i + 1) % n], poly[i]), sub(x, poly[i]))
        if c > 0:
            pos += 1
        elif c < 0:
            neg += 1
    return pos == 0 or neg == 0


def in_wedge(u0, u1, v):
    if crs(u0, u1) >= 0:
        return crs(u0, v) >= 0 and crs(v, u1) >= 0
    return crs(u0, v) <= 0 and crs(v, u1) <= 0


class Sector:
    def __init__(self, c, r, u0, u1):
        self.c, self.r, self.u0, self.u1 = c, r, u0, u1

    def has(self, x):
        d = sub(x, self.c)
        return nrm(d) <= self.r and in_wedge(self.u0, self.u1, d)


class Poly:
    def __init__(self, pts):
        self.pts = pts

    def has(self, x):
        return in_convex(self.pts, x)


class Disc:
    def __init__(self, c, r):
        self.c, self.r = c, r

    def has(self, x):
        return nrm(sub(x, self.c)) <= self.r


def half_disc(c, r, n, d):
    return [Sector(c, r, n, d), Sector(c, r, d, mul(-1, n))]


def end_curv(seg, at_end):
    """The CR's end curvature generalized (SG8): returns (kappa >= 0, turning sign, is_line)."""
    if seg.n == 1:
        return mpf(0), 0, True
    if at_end:
        idx = [seg.n, seg.n - 1, seg.n - 2]
    else:
        idx = [0, 1, 2]
    q0, q1, q2 = (seg.p[i] for i in idx)
    w0, w1, w2 = (seg.w[i] for i in idx)
    d1 = sub(q1, q0)
    if seg.raw[idx[1]] == seg.raw[idx[0]]:
        return mpf(0), 0, True
    if at_end:
        x = crs(sub(q1, q2), sub(q0, q1))
    else:
        x = crs(d1, sub(q2, q1))
    k = mpf(seg.n - 1) / seg.n * (w0 * w2 / (w1 * w1)) * abs(x) / nrm(d1) ** 3
    return k, (1 if x > 0 else (-1 if x < 0 else 0)), False


def exceeds_inverse_h(seg, at_end, h):
    """|kappa| > 1/h in rationals: (n-1)^2 h^2 w0^2 w2^2 X^2 > n^2 w1^4 |d1|^6."""
    if seg.n == 1:
        return False
    idx = [seg.n, seg.n - 1, seg.n - 2] if at_end else [0, 1, 2]
    q = [tuple(Fraction(c) for c in seg.raw[i]) for i in idx]
    w = [Fraction(float(seg.w[i])) for i in idx]
    if q[1] == q[0]:
        return False
    d1 = (q[1][0] - q[0][0], q[1][1] - q[0][1])
    d2 = (q[2][0] - q[1][0], q[2][1] - q[1][1])
    X = d1[0] * d2[1] - d1[1] * d2[0]
    h = Fraction(h)
    n = seg.n
    lhs = (n - 1) ** 2 * h * h * w[0] ** 2 * w[2] ** 2 * X * X
    rhs = n * n * w[1] ** 4 * (d1[0] ** 2 + d1[1] ** 2) ** 3
    return lhs > rhs


class ArcsRegion:
    """The 2018 CR's arcs join, as SG8 settles its gaps."""

    def __init__(self, Vv, A, B, ah, bh, ca, cb, X, clip):
        self.V, self.A, self.B, self.ah, self.bh = Vv, A, B, ah, bh
        self.ca, self.cb = ca, cb  # None for a line, else (centre, radius)
        self.X = X
        self.clip = clip  # None or (Y, T)

    @staticmethod
    def _seen(x, a, b):
        p, q = sub(a, x), sub(b, x)
        return mpmath.atan2(crs(p, q), dot(p, q))

    @staticmethod
    def _extent(z, s, d, e):
        sense = 1 if crs(sub(s, z), d) >= 0 else -1
        a0 = mpmath.atan2(s[1] - z[1], s[0] - z[0])
        a1 = mpmath.atan2(e[1] - z[1], e[0] - z[0])
        sweep = sense * (a1 - a0)
        while sweep < 0:
            sweep += 2 * mpmath.pi
        while sweep >= 2 * mpmath.pi:
            sweep -= 2 * mpmath.pi
        return a0, sense, sweep

    def _arc_seen(self, x, z, r, a0, sense, sweep):
        """The angle an arc sweeps seen from x: each piece of at most an eighth of a turn sweeps its chord's
        angle, plus a whole turn in the arc's sense when x lies between the chord and the arc."""
        n = max(1, int(mpmath.ceil(sweep / (mpmath.pi / 4))))
        total = mpf(0)
        inside = nrm(sub(x, z)) < r
        for k in range(n):
            t0 = a0 + sense * sweep * k / n
            t1 = a0 + sense * sweep * (k + 1) / n
            p0 = add(z, mul(r, (mpmath.cos(t0), mpmath.sin(t0))))
            p1 = add(z, mul(r, (mpmath.cos(t1), mpmath.sin(t1))))
            total += self._seen(x, p0, p1)
            if inside:
                sx, sz = crs(sub(p1, p0), sub(x, p0)), crs(sub(p1, p0), sub(z, p0))
                if (sx > 0 and sz < 0) or (sx < 0 and sz > 0):
                    total += sense * 2 * mpmath.pi
        return total

    def has(self, x):
        """The region the closed boundary V, A, A's arc (or line) to X, B's arc (or line) back to B, V
        encloses (nonzero winding), before the clip.  A winding number, not a ray from V: a reversal's arcs
        can run half way round to their far meeting point and curl past V."""
        if self.clip is not None:
            Y, T = self.clip
            if dot(sub(x, Y), T) > 0:
                return False
        if crs(sub(self.A, self.V), sub(self.B, self.V)) != 0 and in_convex([self.V, self.A, self.B], x):
            return True
        total = self._seen(x, self.V, self.A)
        if self.ca is None:
            total += self._seen(x, self.A, self.X)
        else:
            z, r = self.ca
            a0, sense, sweep = self._extent(z, self.A, self.ah, self.X)
            total += self._arc_seen(x, z, r, a0, sense, sweep)
        if self.cb is None:
            total += self._seen(x, self.X, self.B)
        else:
            z, r = self.cb
            a0, sense, sweep = self._extent(z, self.B, mul(-1, self.bh), self.X)
            ax = mpmath.atan2(self.X[1] - z[1], self.X[0] - z[0])
            total += self._arc_seen(x, z, r, ax, -sense, sweep)
        total += self._seen(x, self.B, self.V)
        return abs(total) > mpmath.pi


def rectangle(Vv, n, ah, h, L):
    return Poly([add(Vv, mul(h, n)), add(add(Vv, mul(h, n)), mul(L * h, ah)),
                 add(sub(Vv, mul(h, n)), mul(L * h, ah)), sub(Vv, mul(h, n))])


def miter_clip_shape(Vv, ah, bh, oa, ob, h, L, within):
    A = add(Vv, mul(h, oa))
    B = add(Vv, mul(h, ob))
    u = unit(add(oa, ob))
    if within:
        T = add(Vv, mul(h / (1 + dot(oa, ob)), add(oa, ob)))
        return Poly([Vv, A, T, B])
    la = (L * h - dot(sub(A, Vv), u)) / dot(ah, u)
    lb = (L * h - dot(sub(B, Vv), u)) / dot(mul(-1, bh), u)
    return Poly([Vv, A, add(A, mul(la, ah)), sub(B, mul(lb, bh)), B])


def arcs_join(segin, segout, Vv, a, b, h, L, outer_left, record):
    """Returns a list of shapes; records the constructions when `record` is a list."""
    ah, bh = unit(a), unit(b)
    oa = left(ah) if outer_left else right(ah)
    ob = left(bh) if outer_left else right(bh)
    A, B = add(Vv, mul(h, oa)), add(Vv, mul(h, ob))
    reversal = crs(a, b) == 0 and dot(a, b) < 0
    ka, sa, la = end_curv(segin, True)
    kb, sb, lb = end_curv(segout, False)
    # signed positive toward the inner side: the inner side is the right when the outer is the left
    ka = ka * sa * (1 if outer_left else -1)
    kb = kb * sb * (1 if outer_left else -1)
    if ka == 0 and kb == 0:
        if reversal:
            return [rectangle(Vv, oa, ah, h, L)], "miter-clip"
        within = miter_within_vec(a, b, L)
        return [miter_clip_shape(Vv, ah, bh, oa, ob, h, L, within)], "miter-clip"
    if exceeds_inverse_h(segin, True, h) or exceeds_inverse_h(segout, False, h):
        if reversal:
            return half_disc(Vv, h, oa, ah), "round"
        return [Sector(Vv, h, oa, ob)], "round"
    ca = None if ka == 0 else (sub(Vv, mul(1 / ka, oa)), abs(1 / ka + h))
    cb = None if kb == 0 else (sub(Vv, mul(1 / kb, ob)), abs(1 / kb + h))

    def nu(k, o):
        return mul(-1, o) if (1 / k + h) > 0 else o

    cands = []
    branch = "meet"
    if ca is None and cb is None:
        den = crs(ah, bh)
        if den != 0:
            t = crs(sub(B, A), bh) / den
            cands = [add(A, mul(t, ah))]
    elif ca is None or cb is None:
        z, r = cb if ca is None else ca
        p0, d0 = (A, ah) if ca is None else (B, bh)
        m = sub(p0, z)
        bq = dot(m, d0)
        cq = dot(m, m) - r * r
        disc = bq * bq - cq
        if disc >= 0:
            sq = mpmath.sqrt(disc)
            cands = [add(p0, mul(-bq - sq, d0)), add(p0, mul(-bq + sq, d0))]
    else:
        (z1, r1), (z2, r2) = ca, cb
        dz = sub(z2, z1)
        d = nrm(dz)
        if d != 0 and abs(r1 - r2) <= d <= r1 + r2:
            aa = (r1 * r1 - r2 * r2 + d * d) / (2 * d)
            hh = mpmath.sqrt(max(r1 * r1 - aa * aa, 0))
            e = unit(dz)
            mm = add(z1, mul(aa, e))
            cands = [add(mm, mul(hh, (-e[1], e[0]))), sub(mm, mul(hh, (-e[1], e[0])))]
    if len(cands) == 2 and (ca is None) != (cb is None):
        # a line's extension is a ray: forward from A along a, backward from B along -b (SG8)
        p0, d0 = (A, ah) if ca is None else (B, mul(-1, bh))
        cands = [c for c in cands if dot(sub(c, p0), d0) >= 0]
    if not cands:
        if ca is not None and cb is not None:
            (z1, r1), (z2, r2) = ca, cb
            n1, n2 = nu(ka, oa), nu(kb, ob)
            dz = nrm(sub(z1, z2))
            if dz < abs(r1 - r2):
                branch = "enclosing"
                big_a = r1 > r2
                s1, s2 = (-1, 1) if big_a else (1, -1)
            else:
                branch = "apart"
                s1, s2 = 1, 1

            def gap(dl):
                za = add(z1, mul(s1 * dl, n1))
                zb = add(z2, mul(s2 * dl, n2))
                ra, rb = r1 + s1 * dl, r2 + s2 * dl
                dd = nrm(sub(za, zb))
                return (dd - abs(ra - rb)) if branch == "enclosing" else (dd - (ra + rb))

            dl = first_zero(gap)
            if dl is not None:
                za = add(z1, mul(s1 * dl, n1))
                zb = add(z2, mul(s2 * dl, n2))
                ra, rb = r1 + s1 * dl, r2 + s2 * dl
                ca, cb = (za, ra), (zb, rb)
                if branch == "enclosing" and not big_a:
                    cands = [add(zb, mul(rb, unit(sub(za, zb))))]
                else:
                    cands = [add(za, mul(ra, unit(sub(zb, za))))]
        elif ca is not None or cb is not None:
            branch = "grow"
            on_a = ca is not None
            z, r = ca if on_a else cb
            n_ = nu(ka, oa) if on_a else nu(kb, ob)
            p0, d0 = (B, bh) if on_a else (A, ah)
            m = left(d0)

            def gap(dl):
                zz = add(z, mul(dl, n_))
                return abs(dot(sub(zz, p0), m)) - (r + dl)

            ray = mul(-1, bh) if on_a else ah  # the line extension's direction of travel (a: forward; b: backward)

            def on_ray(dl):
                zz = add(z, mul(dl, n_))
                foot = sub(zz, mul(dot(sub(zz, p0), m), m))
                return dot(sub(foot, p0), ray) >= 0

            dl = first_zero(gap, on_ray)
            if dl is not None:
                zz = add(z, mul(dl, n_))
                if on_a:
                    ca = (zz, r + dl)
                else:
                    cb = (zz, r + dl)
                cands = [sub(zz, mul(dot(sub(zz, p0), m), m))]
        if not cands:
            if record is not None:
                record.append({"X": None, "clip": None, "branch": "rectangle"})
            return [rectangle(Vv, oa, ah, h, L)], "rectangle"
    def along_a(p):
        """How far along A's extension, from A in its direction of travel, a point lies (SG8's tie rule)."""
        if ca is None:
            return dot(sub(p, A), ah)
        z = ca[0]
        sense = 1 if crs(sub(A, z), ah) >= 0 else -1
        ang = sense * (mpmath.atan2(p[1] - z[1], p[0] - z[0]) - mpmath.atan2(A[1] - z[1], A[0] - z[0]))
        while ang < 0:
            ang += 2 * mpmath.pi
        return ang

    if len(cands) == 2:
        # The closer candidate is the one on V's side of the line of centres (two circles) or nearer V's foot along the
        # line (a circle and a line): the two are mirror images across it.  They tie exactly when V is on it, which is
        # a reversal; the tie goes to the one reached first along A's extension.
        if ca is not None and cb is not None:
            side = -crs(sub(cb[0], ca[0]), sub(Vv, ca[0]))
        else:
            z = cb[0] if ca is None else ca[0]
            p0, d0 = (A, ah) if ca is None else (B, bh)
            side = dot(sub(p0, z), d0) + dot(sub(Vv, p0), d0)
        if reversal or side == 0:
            cands.sort(key=along_a)
        else:
            cands = [cands[1] if side > 0 else cands[0]]
    X = cands[0]
    u = ah if reversal else unit(add(oa, ob))
    d = sub(X, Vv)
    dn = nrm(d)
    clip = None
    Lh = L * h
    if dn > 0:
        cr, dt = crs(u, d), dot(u, d)
        if cr == 0:
            if dt <= 0 or dn > Lh:
                clip = (add(Vv, mul(Lh, u)), u)
        else:
            w = right(u) if cr > 0 else left(u)
            if dot(w, d) < 0:
                w = mul(-1, w)
            rho = dn * dn / (2 * dot(d, w))
            beta = mpmath.atan2(abs(cr), dt)
            if 2 * beta * rho > Lh:
                psi = Lh / rho
                O = add(Vv, mul(rho, w))
                Y = add(O, mul(rho, add(mul(-mpmath.cos(psi), w), mul(mpmath.sin(psi), u))))
                T = add(mul(mpmath.sin(psi), w), mul(mpmath.cos(psi), u))
                clip = (Y, T)
    if record is not None:
        record.append({"X": X, "clip": clip, "branch": branch})
    return [ArcsRegion(Vv, A, B, ah, bh, ca, cb, X, clip)], branch


def first_zero(g, valid=None):
    """The least positive root of g (with valid(root), if given), by stepping outward then bisecting."""
    start = mpf(0)
    for _ in range(8):
        r = first_zero_from(g, start)
        if r is None:
            return None
        if valid is None or valid(r):
            return r
        start = r + mpf("1e-9")
    return None


def first_zero_from(g, start):
    lo = start
    g0 = g(lo)
    step = mpf("1e-6")
    hi = lo + step
    # up to 1e-6 * 2^60 (about 1e12): far past any scale of the cases, and short of where the 30-digit arithmetic
    # would lose the difference of the two sides
    for _ in range(60):
        gh = g(hi)
        if (gh > 0) != (g0 > 0) or gh == 0:
            break
        lo = hi
        hi = lo + 2 * (hi - start)
    else:
        return None
    for _ in range(300):
        mid = (lo + hi) / 2
        gm = g(mid)
        if (gm > 0) == (g(lo) > 0):
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def miter_within_vec(a, b, L):
    na, nb = nrm(a), nrm(b)
    return (na * nb + dot(a, b)) * L * L >= 2 * na * nb


# --- a stroke --------------------------------------------------------------------------


class Body:
    def __init__(self, seg, h):
        self.seg, self.h = seg, h
        # a box the body cannot leave: the control points' box outset by h
        xs = [q[0] for q in seg.p]
        ys = [q[1] for q in seg.p]
        self.box = (min(xs) - h, min(ys) - h, max(xs) + h, max(ys) + h)
        # cusps: parameters in (0, 1) where the derivative vanishes, at minima of |E| refined by golden section
        self.cusps = []
        N = 512
        vals = [nrm(seg.deriv(mpf(i) / N)) for i in range(N + 1)]
        scale = max(vals)
        for i in range(1, N):
            if vals[i] <= vals[i - 1] and vals[i] <= vals[i + 1]:
                a, b = mpf(i - 1) / N, mpf(i + 1) / N
                for _ in range(160):
                    m1 = a + (b - a) / 3
                    m2 = b - (b - a) / 3
                    if nrm(seg.deriv(m1)) < nrm(seg.deriv(m2)):
                        b = m2
                    else:
                        a = m1
                t = (a + b) / 2
                if nrm(seg.deriv(t)) <= scale * mpf("1e-20") and 0 < t < 1:
                    self.cusps.append(t)
        self.ts = [mpf(i) / 96 for i in range(97)]

    def has(self, x):
        seg, h = self.seg, self.h
        if not (self.box[0] <= x[0] <= self.box[2] and self.box[1] <= x[1] <= self.box[3]):
            return False
        for t in self.cusps:
            if nrm(sub(x, seg.point(t))) <= h:
                return True

        def f(t):
            return dot(sub(x, seg.point(t)), seg.deriv(t))

        ts = self.ts
        fs = [f(t) for t in ts]
        for i in range(len(ts)):
            if fs[i] == 0 and nrm(sub(x, seg.point(ts[i]))) <= h and nrm(seg.deriv(ts[i])) > 0:
                return True
        for i in range(len(ts) - 1):
            if (fs[i] > 0) != (fs[i + 1] > 0) and fs[i] != 0 and fs[i + 1] != 0:
                a, b = ts[i], ts[i + 1]
                fa = fs[i]
                for _ in range(90):
                    m = (a + b) / 2
                    fm = f(m)
                    if (fm > 0) == (fa > 0):
                        a, fa = m, fm
                    else:
                        b = m
                t = (a + b) / 2
                if nrm(sub(x, seg.point(t))) <= h and nrm(seg.deriv(t)) > 0:
                    return True
        return False


def tangents(subs, fallback):
    """SG11: the effective tangent of each degenerate subpath."""
    out = {}
    for i, s in enumerate(subs):
        if s["segs"]:
            continue
        d = None
        for j in range(i - 1, -1, -1):
            if subs[j]["segs"]:
                d = subs[j]["segs"][-1].end_dir()[0]
                break
        if d is None:
            for j in range(i + 1, len(subs)):
                if subs[j]["segs"]:
                    d = subs[j]["segs"][0].start_dir()[0]
                    break
        out[i] = d if d is not None else fallback
    return out


def build(subs, style, h, take, tang, record):
    shapes = []
    jn, L, cap = style["join"], mpf(style["limit"]), style["cap"]
    for i, sb in enumerate(subs):
        if not take(i):
            continue
        segs = sb["segs"]
        if not segs:
            if not sb["draws"] or style["degenerate"] == 1:
                continue
            P = V(*sb["start"])
            d = unit(tang[i])
            m = left(d)
            if cap == 1:
                shapes.append(Disc(P, h))
            elif cap == 2:
                shapes.append(Poly([add(add(P, mul(h, m)), mul(h, d)), add(sub(P, mul(h, m)), mul(h, d)),
                                    sub(sub(P, mul(h, m)), mul(h, d)), sub(add(P, mul(h, m)), mul(h, d))]))
            continue
        for s in segs:
            shapes.append(Body(s, h))
        pairs = list(zip(segs, segs[1:]))
        if sb["closed"]:
            pairs.append((segs[-1], segs[0]))
        for si, so_ in pairs:
            a, qi = si.end_dir()
            b, _ = so_.start_dir()
            Vv = si.p[si.n]
            Q = si.raw[qi]
            S = so_.raw[[k for k in range(1, so_.n + 1) if so_.raw[k] != so_.raw[0]][0]]
            cr, dt = crs(a, b), dot(a, b)
            if cr == 0 and dt > 0:
                continue
            ah, bh = unit(a), unit(b)
            if cr == 0:
                n = left(ah)
                if jn == 1:
                    shapes.extend(half_disc(Vv, h, n, ah))
                elif jn == 3:
                    shapes.append(rectangle(Vv, n, ah, h, L))
                elif jn == 4:
                    for ol in (True, False):
                        sh, _ = arcs_join(si, so_, Vv, a, b, h, L, ol, record)
                        shapes.extend(sh)
                continue
            outer_left = cr > 0
            oa = left(ah) if outer_left else right(ah)
            ob = left(bh) if outer_left else right(bh)
            A, B = add(Vv, mul(h, oa)), add(Vv, mul(h, ob))
            shapes.append(Poly([Vv, A, B]))
            if jn == 0:
                if miter_within(Q, si.raw[si.n], S, style["limit"]):
                    shapes.append(miter_clip_shape(Vv, ah, bh, oa, ob, h, L, True))
            elif jn == 3:
                shapes.append(miter_clip_shape(Vv, ah, bh, oa, ob, h, L,
                                               miter_within(Q, si.raw[si.n], S, style["limit"])))
            elif jn == 1:
                shapes.append(Sector(Vv, h, oa, ob))
            elif jn == 4:
                sh, _ = arcs_join(si, so_, Vv, a, b, h, L, outer_left, record)
                shapes.extend(sh)
        if not sb["closed"]:
            for E, d in ((segs[0].p[0], mul(-1, segs[0].start_dir()[0])),
                         (segs[-1].p[segs[-1].n], segs[-1].end_dir()[0])):
                dh = unit(d)
                m = left(dh)
                if cap == 2:
                    shapes.append(Poly([add(E, mul(h, m)), add(add(E, mul(h, m)), mul(h, dh)),
                                        add(sub(E, mul(h, m)), mul(h, dh)), sub(E, mul(h, m))]))
                elif cap == 1:
                    shapes.extend(half_disc(E, h, m, dh))
    return shapes


def stroke(case, record=None):
    """(body shapes, aligned shapes, model-to-device map, local?)."""
    subs = read(case["ops"])
    M = [float(v) for v in case["map"]]
    local = case["space"] == 0
    if not local:
        for s in subs:
            s["segs"] = [m for m in (mapped(g, M) for g in s["segs"]) if m.extent()]
            x, y = s["start"]
            s["start"] = (M[0] * x + M[1] * y + M[4], M[2] * x + M[3] * y + M[5])
        fallback = V(M[0], M[2])
    else:
        fallback = V(1, 0)
    tang = tangents(subs, fallback)
    h = mpf(case["width"]) / 2
    style = case

    def cwe(i):
        return subs[i]["closed"] and bool(subs[i]["segs"])

    if case.get("alignment", 0) == 0:
        body = build(subs, style, h, lambda i: True, tang, record)
        aligned = []
    else:
        body = build(subs, style, h, lambda i: not cwe(i), tang, record)
        aligned = build(subs, style, 2 * h, cwe, tang, record)
    return body, aligned, M, local


def inside(shapes, M, local, x):
    if local:
        det = mpf(M[0]) * M[3] - mpf(M[1]) * M[2]
        rx, ry = x[0] - M[4], x[1] - M[5]
        x = ((M[3] * rx - M[1] * ry) / det, (-M[2] * rx + M[0] * ry) / det)
    return any(s.has(x) for s in shapes)


def decided(shapes, M, local, x, r=mpf("0.2")):
    v = inside(shapes, M, local, x)
    for k in range(16):
        a = 2 * mpmath.pi * k / 16
        y = (x[0] + r * mpmath.cos(a), x[1] + r * mpmath.sin(a))
        if inside(shapes, M, local, y) != v:
            return None
    return v


# --- the cases ------------------------------------------------------------------------

ID = [1, 0, 0, 1, 0, 0]


def case(name, ops, width, join="miter", limit=4, cap="butt", degenerate=0, space=0, mp_=ID, alignment=0,
         probes=(), box=None, grid=8, arcs=False, empty_ok=False):
    return {"empty_ok": empty_ok, "name": name, "ops": ops, "width": width, "join": JOIN[join], "limit": limit, "cap": CAP[cap],
            "degenerate": degenerate, "space": space, "map": mp_, "alignment": alignment, "probes": probes,
            "box": box, "grid": grid, "arcs": arcs}


BELOW_125 = math.nextafter(1.25, 0)
TRI = [("M", [0, 0]), ("L", [4, 3]), ("L", [8, 0])]
REV = [("M", [0, 0]), ("L", [10, 0]), ("L", [4, 0])]

CASES = [
    # SG9: the ratio is exactly 5/4 here
    case("miter_at_limit", TRI, 8, "miter", 1.25, probes=[(4, 7.2), (3.6, 6.8)], box=(0, 0, 8, 9)),
    case("miter_below_limit", TRI, 8, "miter", BELOW_125, probes=[(4, 7.2), (3.6, 6.8)], box=(0, 0, 8, 9)),
    case("miter_clip_below_limit", TRI, 8, "miter-clip", BELOW_125, probes=[(4, 7.2)], box=(0, 0, 8, 9)),
    case("miter_clip_one_one", TRI, 8, "miter-clip", 1.1, probes=[(4, 7.2), (4, 6.8)], box=(0, 0, 8, 9)),
    case("bevel_right_angle", [("M", [0, 0]), ("L", [6, 0]), ("L", [6, 6])], 2, "bevel", box=(3, -3, 9, 3)),
    case("miter_acute", [("M", [0, 0]), ("L", [10, 0]), ("L", [1, 2])], 1, "miter", 10, box=(6, -3, 16, 3)),
    case("round_right_angle", [("M", [0, 0]), ("L", [6, 0]), ("L", [6, 6])], 3, "round", box=(3, -3, 9, 3)),
    # a round join is a sector, not a disc: beside a short segment with a butt end
    case("round_not_disc", [("M", [0, 0]), ("L", [0.5, 0]), ("L", [0.5, 8])], 4, "round", cap="butt",
         probes=[(-1.2, 1.2), (-1.2, -1.2)], box=(-3, -3, 3, 3)),
    case("reversal_round", REV, 2, "round", probes=[(10.6, 0.3)], box=(7, -2, 12, 2)),
    case("reversal_miter_clip", REV, 2, "miter-clip", 3, probes=[(12.5, 0.5), (13.4, 0)], box=(8, -2, 15, 2)),
    case("reversal_miter", REV, 2, "miter", 10, probes=[(10.5, 0.5)], box=(8, -2, 13, 2)),
    case("reversal_bevel", REV, 2, "bevel", probes=[(10.5, 0.5)], box=(8, -2, 13, 2)),
    case("cusp_disc", [("M", [0, 0]), ("C", [8, 4, 4, 4, 4, 0])], 3, "miter", probes=[(5.2, 3.0), (6.8, 3.0)],
         box=(1, 0, 9, 6)),
    case("smooth_vertex", [("M", [0, 0]), ("Q", [4, 0, 8, 4]), ("L", [12, 8])], 2, "miter", 4, box=(5, -1, 11, 5)),
    case("square_caps", [("M", [0, 0]), ("L", [6, 2])], 2, "miter", cap="square", box=(-3, -3, 9, 5)),
    case("round_caps_curve", [("M", [0, 0]), ("Q", [3, 6, 6, 0])], 2, "miter", cap="round", box=(-2, -2, 8, 5)),
    # SG11: the effective tangent of a degenerate square (h = 3: a probe inside the turned square and outside the
    # axis-aligned one, and one the other way round)
    case("degenerate_alone", [("M", [2, 2]), ("L", [2, 2])], 6, cap="square", probes=[(4.5, 4.5)], box=(-3, -3, 7, 7)),
    case("degenerate_after", [("M", [0, 0]), ("L", [3, 4]), ("M", [10, 10]), ("Z", [])], 6, cap="square",
         probes=[(13.5, 10.5), (12.7, 12.7)], box=(5, 5, 15, 15)),
    case("degenerate_before", [("M", [10, 10]), ("L", [10, 10]), ("M", [0, 0]), ("L", [4, -3])], 6, cap="square",
         probes=[(10.5, 6.5), (12.7, 12.7)], box=(5, 5, 15, 15)),
    case("degenerate_between", [("M", [0, 0]), ("L", [3, 4]), ("M", [10, 10]), ("L", [10, 10]), ("M", [20, 0]),
                                ("L", [25, 0])], 6, cap="square", probes=[(13.5, 10.5), (12.7, 12.7)],
         box=(5, 5, 15, 15)),
    case("degenerate_after_lone_move", [("M", [0, 0]), ("L", [5, 5]), ("M", [3, 3]), ("M", [10, 10]),
                                        ("L", [10, 10])], 6, cap="square", probes=[(13.8, 10.0), (12.7, 12.7)],
         box=(5, 5, 15, 15)),
    case("degenerate_after_closed", [("M", [0, 0]), ("L", [6, 0]), ("L", [6, 6]), ("Z", []), ("M", [20, 20]),
                                     ("Z", [])], 6, cap="square", probes=[(23.8, 20.0), (22.7, 22.7)],
         box=(15, 15, 25, 25)),
    case("degenerate_round", [("M", [2, 2]), ("Z", [])], 6, cap="round", box=(-3, -3, 7, 7)),
    case("degenerate_butt", [("M", [2, 2]), ("Z", [])], 6, cap="butt", box=(-3, -3, 7, 7), empty_ok=True),
    case("degenerate_omit", [("M", [2, 2]), ("Z", [])], 6, cap="square", degenerate=1, box=(-3, -3, 7, 7),
         empty_ok=True),
    case("degenerate_local_skew", [("M", [0, 0]), ("L", [4, 0]), ("M", [6, 0]), ("L", [6, 0])], 2, cap="square",
         mp_=[1, 1, 0, 1, 0, 0], box=(4, -2, 10, 2)),
    case("degenerate_device_skew", [("M", [0, 0]), ("L", [4, 0]), ("M", [6, 0]), ("L", [6, 0])], 2, cap="square",
         space=1, mp_=[1, 1, 0, 1, 0, 0], box=(4, -2, 10, 2)),
    case("degenerate_device_alone", [("M", [3, 3]), ("Z", [])], 2, cap="square", space=1, mp_=[0, -2, 1, 0, 0, 0],
         box=(-9, 0, -3, 6)),
    case("local_scale", [("M", [0, 0]), ("L", [4, 0]), ("L", [4, 4])], 1, "miter", 4, mp_=[3, 0, 0, 1, 0, 0],
         box=(8, -2, 16, 6)),
    case("device_scale", [("M", [0, 0]), ("L", [4, 0]), ("L", [4, 4])], 1, "miter", 4, space=1,
         mp_=[3, 0, 0, 1, 0, 0], box=(8, -2, 16, 6)),
    case("inside_alignment", [("M", [0, 0]), ("L", [8, 0]), ("L", [8, 8]), ("L", [0, 8]), ("Z", []),
                              ("M", [12, 0]), ("L", [16, 0])], 2, "miter", alignment=1, box=(-3, -3, 18, 11)),
    # the 2018 CR's arcs join
    case("arcs_two_lines", TRI, 8, "arcs", 1.1, box=(0, 0, 8, 9), arcs=True),
    case("arcs_curve_line", [("M", [0, 0]), ("Q", [6, 0, 8, 4]), ("L", [6, 10])], 2, "arcs", 8, box=(5, 0, 12, 8),
         arcs=True),
    case("arcs_meet_clipped", [("M", [0, 0]), ("Q", [5, 3, 10, 4]), ("Q", [12, 6, 1, 17])], 2, "arcs", 1,
         box=(6, 0, 16, 10), grid=12, arcs=True),
    case("arcs_meet", [("M", [0, 0]), ("Q", [11, -6, 10, 4]), ("Q", [11, 8, 4, 11])], 1, "arcs", 30,
         box=(6, 0, 16, 10), grid=12, arcs=True),
    case("arcs_apart_clipped", [("M", [0, 0]), ("Q", [14, 0, 10, 4]), ("Q", [13, 9, -4, 14])], 4, "arcs", 1.5,
         box=(4, -2, 18, 12), grid=12, arcs=True),
    case("arcs_apart", [("M", [0, 0]), ("Q", [12, 4, 10, 4]), ("Q", [15, 14, -2, 10])], 2, "arcs", 30,
         box=(4, -2, 18, 12), grid=12, arcs=True),
    case("arcs_enclosing_clipped", [("M", [0, 0]), ("Q", [5, 0, 10, 4]), ("Q", [6, 4, 1, 13])], 4, "arcs", 4,
         box=(4, -2, 18, 12), grid=12, arcs=True),
    case("arcs_enclosing", [("M", [0, 0]), ("Q", [11, 6, 10, 4]), ("Q", [15, 5, -2, 19])], 1, "arcs", 30,
         box=(4, -2, 18, 12), grid=12, arcs=True),
    case("arcs_grow_clipped", [("M", [0, 0]), ("Q", [10, 4, 10, 4]), ("Q", [5, 6, 5, 15])], 4, "arcs", 1,
         box=(4, -2, 18, 12), grid=12, arcs=True),
    case("arcs_first_difference_vanishes", [("M", [0, 0]), ("C", [0, 0, 6, 2, 8, 4]), ("Q", [6, 8, 0, 8])], 2, "arcs",
         2, box=(5, 0, 12, 8), arcs=True),
    case("arcs_reversal_enclosing", [("M", [0, 0]), ("Q", [6, -1, 10, 0]), ("Q", [2, -2, 8, 8])], 1, "arcs", 4,
         box=(6, -4, 16, 4), grid=12, arcs=True),
    case("arcs_reversal_retrace", [("M", [0, 0]), ("Q", [5, 3, 10, 0]), ("Q", [5, 3, 0, 0])], 2, "arcs", 4,
         box=(6, -5, 17, 5), grid=12, arcs=True),
    case("arcs_reversal_rectangle", [("M", [0, 0]), ("L", [10, 0]), ("Q", [6, 0, 2, 5])], 2, "arcs", 4,
         box=(6, -4, 16, 4), grid=12, arcs=True),
    case("arcs_local_nonuniform", [("M", [0, 0]), ("Q", [5, 3, 10, 4]), ("Q", [12, 6, 1, 17])], 2, "arcs", 1,
         mp_=[1.5, 0, 0, 1, 0, 0], box=(9, 0, 24, 10), grid=12, arcs=True),
    case("arcs_device", [("M", [0, 0]), ("Q", [5, 3, 10, 4]), ("Q", [12, 6, 1, 17])], 2, "arcs", 1, space=1,
         mp_=[1.5, 0, 0, 1, 0, 0], box=(9, 0, 24, 10), grid=12, arcs=True),
]


def kappa_exactly_inverse_h_case():
    """A quadratic whose end curvature is exactly 1/h: kappa(1) = (1/2) |X| / |d1|^3 with d1 = (2, 0), X = 16:
    kappa = 1, so h = 1 (width 2) is the boundary; width 2 + an ulp is past it."""
    ops = [("M", [0, 0]), ("Q", [0, 8, 2, 8]), ("L", [2, 12])]
    return [case("arcs_kappa_at_inverse_h", ops, 2, "arcs", 6, box=(-1, 5, 6, 12), arcs=True),
            case("arcs_kappa_past_inverse_h", ops, math.nextafter(2.0, 3.0), "arcs", 6, box=(-1, 5, 6, 12),
                 arcs=True)]


CASES += kappa_exactly_inverse_h_case()


def hexf(v):
    return float(v).hex()


def compute(ci):
    """The decided points of case ci, and its arcs constructions."""
    c = CASES[ci]
    rec = [] if c["arcs"] else None
    body, aligned, M, local = stroke(c, rec)
    pts = []
    problems = []
    for pr in c["probes"]:
        x = (mpf(pr[0]), mpf(pr[1]))
        v = decided(body, M, local, x)
        if v is None:
            problems.append("probe %r is not decided" % (pr,))
            continue
        pts.append((x, v, 0))
    x0, y0, x1, y1 = c["box"]
    g = c["grid"]
    for i in range(g + 1):
        for j in range(g + 1):
            x = (mpf(x0) + (mpf(x1) - x0) * i / g, mpf(y0) + (mpf(y1) - y0) * j / g)
            for part, shapes in ((0, body), (1, aligned)):
                if part == 1 and not aligned:
                    continue
                v = decided(shapes, M, local, x)
                if v is not None:
                    pts.append((x, v, part))
    if (not any(p[1] for p in pts) or not any(not p[1] for p in pts)) and not (c["empty_ok"] and not any(p[1] for p in pts)):
        problems.append("every point has one verdict; the case measures nothing")
    if problems:
        return ("%s: %s" % (c["name"], "; ".join(problems)), None, None)
    rows = []
    if rec is not None:
        # numbered among the arcs regions only (a rectangle is a polygon, not an arcs region)
        for r in rec:
            if r["X"] is not None:
                rows.append((len(rows), r["X"], r["clip"], r["branch"]))
    return (None, pts, rows)


def main():
    import multiprocessing
    with multiprocessing.Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 2) as pool:
        results = pool.map(compute, range(len(CASES)))
    errors = [err for err, _, _ in results if err]
    if errors:
        raise SystemExit("\n".join(errors))
    out = sys.stdout
    out.write("/* Generated by tools/gen-stroke-vectors.py; do not edit. */\n")
    out.write("/* SPDX-License-Identifier: LGPL-3.0-only */\n")
    out.write("#ifndef GHOTI_IO_GSHGEOM_TESTS_STROKE_VECTORS_H\n#define GHOTI_IO_GSHGEOM_TESTS_STROKE_VECTORS_H\n\n")
    out.write("#include <cstddef>\n\nnamespace stroke_vectors {\n\n")
    out.write("struct Op {\n  char verb;\n  size_t n;\n  double v[10];\n};\n\n")
    out.write("struct Point {\n  double x, y;\n  int inside;\n  int part; ///< 0 body, 1 aligned\n};\n\n")
    out.write("struct Case {\n  const char * name;\n  const Op * ops;\n  size_t op_count;\n  double width;\n"
              "  int space, join;\n  double limit;\n  int cap, degenerate, alignment;\n  double map[6];\n"
              "  const Point * points;\n  size_t point_count;\n};\n\n")
    out.write("struct Arcs {\n  size_t case_index;\n  size_t join_index;\n  double x[2];\n  int clipped;\n"
              "  double y[2], t[2];\n  const char * branch;\n};\n\n")
    arcs_rows = []
    total = 0
    for ci, (c, (_, pts, rows)) in enumerate(zip(CASES, results)):
        out.write("static const Op ops_%d[] = {\n" % ci)
        for verb, a in c["ops"]:
            vals = ", ".join(hexf(v) for v in a) if a else "0"
            out.write("    {'%s', %d, {%s}},\n" % (verb, len(a), vals))
        out.write("};\n")
        out.write("static const Point points_%d[] = {\n" % ci)
        for x, v, part in pts:
            out.write("    {%s, %s, %d, %d},\n" % (hexf(x[0]), hexf(x[1]), 1 if v else 0, part))
        out.write("};\n\n")
        total += len(pts)
        for ji, X, clip, branch in rows:
            arcs_rows.append((ci, ji, X, clip, branch))
    out.write("static const Case cases[] = {\n")
    for ci, (c, (_, pts, _)) in enumerate(zip(CASES, results)):
        out.write("    {\"%s\", ops_%d, %d, %s, %d, %d, %s, %d, %d, %d, {%s}, points_%d, %d},\n" % (
            c["name"], ci, len(c["ops"]), hexf(c["width"]), c["space"], c["join"], hexf(c["limit"]), c["cap"],
            c["degenerate"], c["alignment"], ", ".join(hexf(v) for v in c["map"]), ci, len(pts)))
    out.write("};\nstatic const size_t case_count = sizeof cases / sizeof cases[0];\n\n")
    out.write("static const Arcs arcs[] = {\n")
    for ci, ji, X, clip, branch in arcs_rows:
        if clip is None:
            out.write("    {%d, %d, {%s, %s}, 0, {0, 0}, {0, 0}, \"%s\"},\n" % (ci, ji, hexf(X[0]), hexf(X[1]), branch))
        else:
            Y, T = clip
            out.write("    {%d, %d, {%s, %s}, 1, {%s, %s}, {%s, %s}, \"%s\"},\n" % (
                ci, ji, hexf(X[0]), hexf(X[1]), hexf(Y[0]), hexf(Y[1]), hexf(T[0]), hexf(T[1]), branch))
    out.write("};\nstatic const size_t arcs_count = sizeof arcs / sizeof arcs[0];\n\n")
    out.write("} // namespace stroke_vectors\n\n#endif /* GHOTI_IO_GSHGEOM_TESTS_STROKE_VECTORS_H */\n")
    branches = sorted(set(r[4] for r in arcs_rows))
    sys.stderr.write("cases %d, points %d, arcs constructions %d, branches %s\n" % (
        len(CASES), total, len(arcs_rows), ", ".join(branches)))
    for ci, ji, X, clip, branch in arcs_rows:
        sys.stderr.write("  %-34s join %d: %s%s\n" % (CASES[ci]["name"], ji, branch, " clipped" if clip else ""))


if __name__ == "__main__":
    main()
