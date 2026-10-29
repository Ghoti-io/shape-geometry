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
"""An independent statement of the arcs of arc.h, in exact rationals and 130-digit decimals.

It shares no code with the library.  Geometry that is rational (the midpoint, the quarter-turn
sines and cosines, the radius check of a quarter-turn rotation) is done in `Fraction`;
`sin`, `cos`, `atan2` and `sqrt` are series and Newton in `Decimal` at 130 digits, with a
function defined here (not libm, not the library's `math`).  An input is a Python float, taken
exactly (`Fraction(x)`, `Decimal(x)`); an answer is a `Decimal`, which a caller turns into the
nearest double with `float()`.

What it states (the numbers in the names are arc.h's rules):

* the endpoint form: SVG 2's conversion from endpoint to centre parameterization and
  the correction of out-of-range radii, written from the text of the Candidate
  Recommendation as it stands in its section on implementation notes (numbered B.2.4 and
  B.2.5 in the text read on 2026-10-29; the formulas are the same), each step as the text
  has it: `x1'`, `y1'` by the rotation of half the difference, `Lambda`, the replacement of
  the radii, the centre by the sign rule on the two flags, `theta1` and `dtheta` as signed
  angles between vectors, and the adjustment by the sweep flag;
* the centre form, the tangent form (Canvas `arcTo`) and the bulge (DXF);
* the decomposition into rational quadratics (the meeting of the tangents, `w = cos h`).

Run `python3 tools/arc-oracle.py` to check the oracle against itself: that the ends lie on the
ellipse it returns to 1e-100, that the flags select the four arcs, and the identities of the
pieces.
"""
from decimal import Decimal, getcontext
from fractions import Fraction

getcontext().prec = 130
D = Decimal
ZERO = D(0)
ONE = D(1)
TWO = D(2)


# --- elementary functions ---------------------------------------------------

def _atan_series(x):
    # |x| <= 0.01
    x2 = x * x
    term = x
    total = x
    n = 1
    eps = D(10) ** (-(getcontext().prec + 5))
    while abs(term) > eps:
        term = -term * x2
        n += 2
        total += term / n
    return total


def _compute_pi():
    # Machin: pi = 16 atan(1/5) - 4 atan(1/239), each by the plain series
    def atan_inv(n):
        x = ONE / n
        x2 = x * x
        term = x
        total = x
        k = 1
        eps = D(10) ** (-(getcontext().prec + 5))
        while abs(term) > eps:
            term = -term * x2
            k += 2
            total += term / k
        return total
    return 16 * atan_inv(5) - 4 * atan_inv(239)


PI = _compute_pi()


def atan(x):
    x = D(x)
    if x < 0:
        return -atan(-x)
    if x > 1:
        return PI / 2 - atan(ONE / x)
    k = 0
    while x > D("0.01"):
        x = x / (ONE + (ONE + x * x).sqrt())
        k += 1
    return _atan_series(x) * (2 ** k)


def atan2(y, x):
    y, x = D(y), D(x)
    if x == 0:
        if y == 0:
            return ZERO
        return PI / 2 if y > 0 else -PI / 2
    a = atan(y / x)
    if x > 0:
        return a
    return a + PI if y >= 0 else a - PI


def _sincos_small(t):
    # |t| <= pi/4 or so: Taylor
    t2 = t * t
    s = t
    c = ONE
    ts = t
    tc = ONE
    eps = D(10) ** (-(getcontext().prec + 5))
    n = 0
    while True:
        n += 1
        tc = -tc * t2 / ((2 * n - 1) * (2 * n))
        ts = -ts * t2 / ((2 * n) * (2 * n + 1))
        c += tc
        s += ts
        if abs(tc) < eps and abs(ts) < eps:
            break
    return s, c


def sincos_rad(t):
    """sin and cos of a Decimal in radians (reduced by quarter turns of pi/2 first)."""
    t = D(t)
    q = (t / (PI / 2)).to_integral_value()
    r = t - q * (PI / 2)
    s, c = _sincos_small(r)
    k = int(q) % 4
    return [(s, c), (c, -s), (-s, -c), (-c, s)][k]


def tan_rad(t):
    s, c = sincos_rad(t)
    return s / c


UNIT_QUARTER = {"radians": None, "degrees": Fraction(90), "gradians": Fraction(100), "turns": Fraction(1, 4)}
UNITS = ["radians", "degrees", "gradians", "turns"]  # the library's enum order


def to_radians(angle, unit):
    """A float angle in a unit to radians (a Decimal)."""
    a = D(angle)
    if unit == "radians":
        return a
    q = UNIT_QUARTER[unit]
    return a * (PI / 2) / D(q.numerator) * D(q.denominator)


def sincos_unit(angle, unit):
    """cos and sin of a float angle in a unit; exact 0 and +-1 at quarter turns (Fraction test)."""
    if unit != "radians":
        q = Fraction(angle) / UNIT_QUARTER[unit]
        if q.denominator == 1:
            k = int(q) % 4
            return [(ONE, ZERO), (ZERO, ONE), (-ONE, ZERO), (ZERO, -ONE)][k]
    s, c = sincos_rad(to_radians(angle, unit))
    return c, s


def turn_radians(unit):
    return PI * 2


# --- the ellipse -------------------------------------------------------------

def ellipse_point(c, u, v, unit, t):
    """P(t) = c + u cos t + v sin t, t a float in the unit."""
    cs, sn = sincos_unit(t, unit)
    return (D(c[0]) + D(u[0]) * cs + D(v[0]) * sn, D(c[1]) + D(u[1]) * cs + D(v[1]) * sn)


def uv_from_radii(rx, ry, phi, unit):
    """u = rx (cos phi, sin phi), v = ry (-sin phi, cos phi) as exact-as-possible decimals."""
    c, s = sincos_unit(phi, unit)
    return (D(rx) * c, D(rx) * s), (-D(ry) * s, D(ry) * c)


# --- the endpoint form: SVG 2, conversion from endpoint to center parameterization ----

def endpoint_arc(S, E, rx, ry, phi, unit, fA, fS):
    """The arc the endpoint form names.

    Returns {'kind': 'nothing'|'line'} or {'kind': 'arc', 'center', 'u', 'v', 'unit': 'radians',
    't0', 'sweep', 'circular', 'lambda'} with Decimal numbers.
    """
    if S[0] == E[0] and S[1] == E[1]:
        return {"kind": "nothing"}
    if rx == 0 or ry == 0:
        return {"kind": "line"}
    rx, ry = abs(rx), abs(ry)
    circular = rx == ry
    x1, y1, x2, y2 = D(S[0]), D(S[1]), D(E[0]), D(E[1])
    cphi, sphi = sincos_unit(phi, unit)
    # step 1: (x1', y1') = [cos sin; -sin cos] ((x1-x2)/2, (y1-y2)/2)
    hx, hy = (x1 - x2) / 2, (y1 - y2) / 2
    x1p = cphi * hx + sphi * hy
    y1p = -sphi * hx + cphi * hy
    # the correction of out-of-range radii: Lambda = x1'^2/rx^2 + y1'^2/ry^2
    rxd, ryd = D(rx), D(ry)
    lam = x1p * x1p / (rxd * rxd) + y1p * y1p / (ryd * ryd)
    corrected = False
    if lam >= ONE - D(10) ** -100:
        # at 1 or above: the replacement rx <- sqrt(Lambda) rx (a no-op at 1), centre at the midpoint
        if lam > ONE:
            rxd, ryd = lam.sqrt() * rxd, lam.sqrt() * ryd
        corrected = True
    if corrected:
        cxp, cyp = ZERO, ZERO
    else:
        num = rxd ** 2 * ryd ** 2 - rxd ** 2 * y1p ** 2 - ryd ** 2 * x1p ** 2
        den = rxd ** 2 * y1p ** 2 + ryd ** 2 * x1p ** 2
        k = (num / den).sqrt()
        if fA == fS:
            k = -k  # the + sign is chosen if fA != fS
        cxp = k * rxd * y1p / ryd
        cyp = -k * ryd * x1p / rxd
    cx = cphi * cxp - sphi * cyp + (x1 + x2) / 2
    cy = sphi * cxp + cphi * cyp + (y1 + y2) / 2
    # theta1 = angle((1, 0), ((x1'-cx')/rx, (y1'-cy')/ry)); dtheta = angle(v1, v2) mod 2 pi, adjusted by fS
    v1 = ((x1p - cxp) / rxd, (y1p - cyp) / ryd)
    v2 = ((-x1p - cxp) / rxd, (-y1p - cyp) / ryd)
    theta1 = atan2(v1[1], v1[0])
    if corrected:
        dtheta = PI if fS else -PI
    else:
        dtheta = atan2(v1[0] * v2[1] - v1[1] * v2[0], v1[0] * v2[0] + v1[1] * v2[1])
        if not fS and dtheta > 0:
            dtheta -= 2 * PI
        elif fS and dtheta < 0:
            dtheta += 2 * PI
    ux, uy = rxd * cphi, rxd * sphi
    vx, vy = -ryd * sphi, ryd * cphi
    return {"kind": "arc", "center": (cx, cy), "u": (ux, uy), "v": (vx, vy), "unit": "radians",
            "t0": theta1, "sweep": dtheta, "circular": circular, "lambda": lam}


# --- the tangent form and the bulge ------------------------------------------

def arc_to(P0, P1, P2, r):
    """Canvas arcTo: ('line',) or ('arcto', T1, C, T2, t0, sweep) with Decimals."""
    if r == 0 or P0 == P1 or P1 == P2:
        return {"kind": "line"}
    p0, p1, p2 = [(D(p[0]), D(p[1])) for p in (P0, P1, P2)]
    cross = (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p1[1] - p0[1]) * (p2[0] - p0[0])
    if cross == 0:
        return {"kind": "line"}
    e0 = (p0[0] - p1[0], p0[1] - p1[1])
    e2 = (p2[0] - p1[0], p2[1] - p1[1])
    l0 = (e0[0] ** 2 + e0[1] ** 2).sqrt()
    l2 = (e2[0] ** 2 + e2[1] ** 2).sqrt()
    e0 = (e0[0] / l0, e0[1] / l0)
    e2 = (e2[0] / l2, e2[1] / l2)
    cosq = e0[0] * e2[0] + e0[1] * e2[1]
    sinq = abs(e0[0] * e2[1] - e0[1] * e2[0])
    dist = D(r) * (ONE + cosq) / sinq  # r / tan(theta / 2)
    t1 = (p1[0] + dist * e0[0], p1[1] + dist * e0[1])
    t2 = (p1[0] + dist * e2[0], p1[1] + dist * e2[1])
    # the centre lies on the bisector at r / sin(theta / 2) from P1
    bis = (e0[0] + e2[0], e0[1] + e2[1])
    bl = (bis[0] ** 2 + bis[1] ** 2).sqrt()
    half_sin = ((ONE - cosq) / 2).sqrt()
    dc = D(r) / half_sin
    c = (p1[0] + dc * bis[0] / bl, p1[1] + dc * bis[1] / bl)
    a = (t1[0] - c[0], t1[1] - c[1])
    b = (t2[0] - c[0], t2[1] - c[1])
    t0 = atan2(a[1], a[0])
    sweep = atan2(a[0] * b[1] - a[1] * b[0], a[0] * b[0] + a[1] * b[1])
    return {"kind": "arcto", "t1": t1, "t2": t2, "center": c, "radius": D(r), "t0": t0, "sweep": sweep,
            "positive": cross > 0}


def bulge_arc(S, Q, bulge):
    """A DXF bulge: {'kind': 'nothing'|'line'|'arc', ...}; bulge = tan(sweep / 4)."""
    if S == Q:
        return {"kind": "nothing"}
    if bulge == 0:
        return {"kind": "line"}
    s, q = (D(S[0]), D(S[1])), (D(Q[0]), D(Q[1]))
    sweep = 4 * atan(D(bulge))
    half = sweep / 2
    chord = ((q[0] - s[0]) ** 2 + (q[1] - s[1]) ** 2).sqrt()
    # the circle through the ends with that sweep; its centre is to the left of the chord by (chord/2)/tan(half)
    sh, ch = sincos_rad(half)
    off = chord / 2 * ch / sh
    nrm = (-(q[1] - s[1]) / chord, (q[0] - s[0]) / chord)
    c = ((s[0] + q[0]) / 2 + off * nrm[0], (s[1] + q[1]) / 2 + off * nrm[1])
    radius = abs(chord / 2 / sh)
    t0 = atan2(s[1] - c[1], s[0] - c[0])
    return {"kind": "arc", "center": c, "radius": radius, "t0": t0, "sweep": sweep}


# --- the decomposition ---------------------------------------------------------

def piece_count(sweep, unit):
    q = UNIT_QUARTER[unit]
    a = abs(Fraction(sweep))
    if unit == "radians":
        qq = Fraction(float.fromhex("0x1.921fb54442d18p+0"))
    else:
        qq = q
    for n in (1, 2, 3):
        if a <= n * qq:
            return n
    return 4


def decompose(c, u, v, unit, t0, sweep):
    """The exact pieces of the arc: a list of (P0, Q, P2, w), Decimals, in real arithmetic."""
    n = piece_count(sweep, unit)
    sg = ONE if sweep > 0 else -ONE
    t0d, swd = D(t0), D(sweep)
    pieces = []
    if unit == "radians":
        scale = ONE
    else:
        scale = (PI / 2) / D(UNIT_QUARTER[unit].numerator) * D(UNIT_QUARTER[unit].denominator)
    h = abs(swd) / (2 * n) * scale
    w = sincos_rad(h)[1]

    def at(k):
        t = t0d + swd * k / n
        cs, sn = sincos_unit_dec(t, unit)
        p = (D(c[0]) + D(u[0]) * cs + D(v[0]) * sn, D(c[1]) + D(u[1]) * cs + D(v[1]) * sn)
        tg = (sg * (-D(u[0]) * sn + D(v[0]) * cs), sg * (-D(u[1]) * sn + D(v[1]) * cs))
        return p, tg
    for k in range(n):
        pa, ta = at(k)
        pb, tb = at(k + 1)
        # the control point is the meeting of the tangent lines at the two ends: pa + s ta = pb + r tb
        den = ta[0] * tb[1] - ta[1] * tb[0]
        s = ((pb[0] - pa[0]) * tb[1] - (pb[1] - pa[1]) * tb[0]) / den
        q = (pa[0] + s * ta[0], pa[1] + s * ta[1])
        pieces.append((pa, q, pb, w))
    return pieces


def sincos_unit_dec(t, unit):
    """cos and sin of a Decimal parameter in a unit (exact at quarter turns)."""
    if unit != "radians":
        q = Fraction(t) / UNIT_QUARTER[unit] if not isinstance(t, Decimal) else Fraction(t) / UNIT_QUARTER[unit]
        if q.denominator == 1:
            return [(ONE, ZERO), (ZERO, ONE), (-ONE, ZERO), (ZERO, -ONE)][int(q) % 4]
        s, c = sincos_rad(D(t) * (PI / 2) / D(UNIT_QUARTER[unit].numerator) * D(UNIT_QUARTER[unit].denominator))
        return c, s
    s, c = sincos_rad(D(t))
    return c, s


# --- self check -----------------------------------------------------------------

def _close(a, b, tol=D(10) ** -100):
    return abs(a - b) <= tol


def _self_check():
    # constants and functions
    assert _close(PI, D("3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460955058223172535940812848111745028410270193852110555964462294895493038196"), D(10) ** -120)
    s, c = sincos_rad(D(1))
    assert _close(s * s + c * c, ONE)
    assert _close(atan2(sincos_rad(D("0.7"))[0], sincos_rad(D("0.7"))[1]), D("0.7"))
    assert _close(atan(ONE), PI / 4)
    # the four arcs of one chord: ends on the ellipse, flags choose the sense and the size
    S, E = (1.0, 0.0), (0.0, 1.0)
    seen = set()
    for phi, unit in ((0.0, "degrees"), (45.5, "degrees"), (-30.0, "degrees"), (90.0, "degrees")):
        for fA in (0, 1):
            for fS in (0, 1):
                a = endpoint_arc(S, E, 1.5, 0.75, phi, unit, fA, fS)
                assert a["kind"] == "arc"
                for t, end in ((a["t0"], S), (a["t0"] + a["sweep"], E)):
                    sn, cs = sincos_rad(t)
                    px = a["center"][0] + a["u"][0] * cs + a["v"][0] * sn
                    py = a["center"][1] + a["u"][1] * cs + a["v"][1] * sn
                    assert _close(px, D(end[0])) and _close(py, D(end[1])), (phi, fA, fS)
                assert (a["sweep"] > 0) == bool(fS)
                assert (abs(a["sweep"]) > PI) == bool(fA)
                seen.add((phi, fA, fS, a["center"]))
    # radius correction: a chord longer than the diameter gives the half ellipse
    a = endpoint_arc((0.0, 0.0), (4.0, 0.0), 1.0, 1.0, 0.0, "degrees", 1, 1)
    assert abs(a["sweep"]) == PI and a["center"] == (D(2), D(0)) and _close(a["u"][0], D(2))
    # pieces: a quarter piece has w = sqrt(1/2) and the end points on the circle
    ps = decompose((0.0, 0.0), (2.0, 0.0), (0.0, 2.0), "degrees", 0.0, 360.0)
    assert len(ps) == 4 and all(_close(p[3], D("0.5").sqrt()) for p in ps)
    assert ps[0][0] == (D(2), D(0)) and ps[0][2] == (D(0), D(2))
    assert _close(ps[0][1][0], D(2)) and _close(ps[0][1][1], D(2))
    print("arc-oracle: self check ok")


if __name__ == "__main__":
    _self_check()
