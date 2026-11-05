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
"""Print tests/hatch_vectors.h: the known answers of hatch.h (HX4 to HX6), in exact rationals.

    python3 tools/gen-hatch-vectors.py > tests/hatch_vectors.h

Shares no code with the library.  The definition is applied, not the library's method:

* The lines of a use are those ``k`` for which the line meets the closed bounding box of the contours (an exact
  test of the four corners).  The line is the one through the double point ``A_k`` with the double direction ``w``
  (formed by the header's formulas, in Python floats, which are IEEE doubles and do not fuse), all of it exact
  rationals after that.
* The trim set is found by the definition of HX4: every position where a contour edge meets the line is a
  breakpoint (the projection of both ends for an edge that lies along the line); the open stretch between two
  consecutive breakpoints is tested at two exact points, one at 2^-100 to the left of its midpoint and one to the
  right, with the winding number and the even-odd parity of every contour taken by casting a ray in exact
  rationals; a stretch is in the set when either point is in the region ``G`` of the nesting style, and
  consecutive such stretches are one interval.  An isolated breakpoint is no interval.
* The dashes are the pattern of the DXF elements (groups 49) laid on the infinite line from the base point
  ``s = 0``, in both directions: a maximal run of pen-down elements (a dash, or a dot adding zero) is one dash, a run
  of spaces one gap, a run of zero length a dot, and a dot first with a dash last keeps the dot apart (HX6, stroke.h SG30).
* The pieces are the dashes cut by the intervals; a dash cut to a point is dropped; a dot inside an interval, ends
  included, is kept.  Local coordinates are ``H(O_k + s u)`` rounded once to the nearest double.
"""
import math
import sys
from fractions import Fraction as F

EPS = F(1, 2 ** 100)


def hexf(v):
    return float(v).hex()


def fl(x):
    """Nearest double of an exact rational."""
    return float(x)


# --- exact sincos for the angles the vectors use (the library's gshgeom_sincos agrees; a test checks it) ----------------
# the doubles the library's ::gshgeom_sincos returns for these degrees (the angle is multiplied by pi/180 once, in
# double, and a test checks that they are these bits): the lines are defined by these doubles, not by the true values
S2 = 0.70710678118654746  # cos 45 and sin 135
S2B = 0.70710678118654757  # sin 45 and -cos 135
S3 = 0.86602540378443871
HF = 0.49999999999999994


def unit_vector(deg):
    t = {0: (1.0, 0.0), 90: (0.0, 1.0), 180: (-1.0, 0.0), 270: (0.0, -1.0),
         45: (S2, S2B), 135: (-S2B, S2),
         30: (S3, HF), 60: (HF, S3)}
    return t[deg % 360]


# --- contours and the region --------------------------------------------------------------------------------------------
def ray_winding(poly, px, py):
    """Winding number and crossing-count parity of closed polygon poly (integer pixels*64) around exact point (px, py)."""
    wind = 0
    cross = 0
    n = len(poly)
    for i in range(n):
        ax, ay = poly[i]
        bx, by = poly[(i + 1) % n]
        if (ay <= py) != (by <= py):
            # x of the crossing at height py
            t = (py - ay) / F(by - ay)
            x = ax + t * (bx - ax)
            if x > px:
                cross += 1
                wind += 1 if by > ay else -1
    return wind, cross & 1


def in_region(contours, nesting, rule, p):
    wind = 0
    odd = 0
    for c in contours:
        w, par = ray_winding(c, p[0], p[1])
        wind += w
        odd += par
    if nesting == 1:
        return odd == 1
    if nesting == 2:
        return odd >= 1
    if rule == 0:
        return wind != 0
    return wind % 2 != 0


def trim_intervals(contours, nesting, rule, A, w):
    """Trim intervals of the line A + s w (exact rationals; contours in 26.6 integers, scaled to pixels here)."""
    ww = w[0] * w[0] + w[1] * w[1]

    def proj(p):
        return ((p[0] - A[0]) * w[0] + (p[1] - A[1]) * w[1]) / ww

    def side(p):
        return w[0] * (p[1] - A[1]) - w[1] * (p[0] - A[0])

    brk = set()
    polys = [[(F(x, 64), F(y, 64)) for (x, y) in c] for c in contours]
    for c in polys:
        n = len(c)
        for i in range(n):
            a, b = c[i], c[(i + 1) % n]
            sa, sb = side(a), side(b)
            if sa == 0 and sb == 0:
                brk.add(proj(a))
                brk.add(proj(b))
            elif sa == 0:
                brk.add(proj(a))
            elif sb == 0:
                brk.add(proj(b))
            elif (sa > 0) != (sb > 0):
                t = sa / (sa - sb)
                x = (a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]))
                brk.add(proj(x))
    bs = sorted(brk)
    nrm = (-w[1], w[0])  # left normal (not unit)
    nl = math.sqrt(float(ww))  # only to scale EPS; the test point need only be closer than any feature
    scale = F(EPS) / F(nl).limit_denominator(10 ** 9)
    segs = []
    for i in range(len(bs) - 1):
        m = (bs[i] + bs[i + 1]) / 2
        M = (A[0] + m * w[0], A[1] + m * w[1])
        ok = False
        for sg in (1, -1):
            q = (M[0] + sg * scale * nrm[0], M[1] + sg * scale * nrm[1])
            if in_region_pix(polys, nesting, rule, q):
                ok = True
        if ok:
            segs.append((bs[i], bs[i + 1]))
    out = []
    for s in segs:
        if out and out[-1][1] == s[0]:
            out[-1] = (out[-1][0], s[1])
        else:
            out.append(s)
    return out


def in_region_pix(polys, nesting, rule, p):
    wind = 0
    odd = 0
    for c in polys:
        n = len(c)
        w = 0
        cr = 0
        for i in range(n):
            ax, ay = c[i]
            bx, by = c[(i + 1) % n]
            if (ay <= p[1]) != (by <= p[1]):
                t = (p[1] - ay) / (by - ay)
                x = ax + t * (bx - ax)
                if x > p[0]:
                    cr += 1
                    w += 1 if by > ay else -1
        wind += w
        odd += cr & 1
    if nesting == 1:
        return odd == 1
    if nesting == 2:
        return odd >= 1
    if rule == 0:
        return wind != 0
    return wind % 2 != 0


# --- dashes ---------------------------------------------------------------------------------------------------------------
def pattern_runs(elems):
    """The runs of the DXF elements: [(kind, start, length)] for one period, kind 'd' (pen down), 'g' (gap); plus the period."""
    runs = []
    pos = F(0)
    for e in elems:
        e = F(e)
        kind = 'd' if e >= 0 else 'g'
        ln = abs(e)
        if runs and runs[-1][0] == kind:
            runs[-1] = (kind, runs[-1][1], runs[-1][2] + ln)
        else:
            runs.append((kind, pos, ln))
        pos += ln
    return runs, pos


def dash_set(elems, lo, hi):
    """Pen-down pieces of the infinite pattern covering [lo, hi] (exact), as (a, b) with a == b a dot; and whether it is solid."""
    runs, period = pattern_runs(elems)
    if not any(r[0] == 'g' for r in runs):
        return [(lo, hi)], True
    # The first and last runs are one run across the repeat when both are pen-down, except that a dot first (a run of
    # no length) is kept apart from a dash last (HX6, SG30).  A dot last next to a dash first is absorbed in it.
    merge = runs[0][0] == 'd' and runs[-1][0] == 'd' and len(runs) > 1 and runs[0][2] != 0
    items = []
    for i, (kind, st, ln) in enumerate(runs):
        if kind != 'd':
            continue
        if merge and i == 0:
            continue
        if merge and i == len(runs) - 1:
            items.append((st, period + runs[0][2]))
        else:
            items.append((st, st + ln))
    out = []
    q0 = math.floor(lo / period) - 2
    q1 = math.floor(hi / period) + 2
    for q in range(int(q0), int(q1) + 1):
        for a, b in items:
            out.append((a + q * period, b + q * period))
    out.sort()
    return out, False


# --- one case -------------------------------------------------------------------------------------------------------------
def apply(m, p):
    # the library's order: (m00 x + m01 y) + t
    return (m[0] * p[0] + m[1] * p[1]) + m[4], (m[2] * p[0] + m[3] * p[1]) + m[5]


def compose(outer, inner):
    # outer . inner as 6-tuples (m00 m01 m10 m11 tx ty), in doubles (integers and halves in the vectors: exact)
    a = outer
    b = inner
    return (a[0] * b[0] + a[1] * b[2], a[0] * b[1] + a[1] * b[3], a[2] * b[0] + a[3] * b[2], a[2] * b[1] + a[3] * b[3],
            (a[0] * b[4] + a[1] * b[5]) + a[4], (a[2] * b[4] + a[3] * b[5]) + a[5])


def case_pieces(contours, nesting, rule, deg, families, H, chain):
    C = compose(chain, H)
    pts = [p for c in contours for p in c]
    lox, hix = min(p[0] for p in pts) / 64.0, max(p[0] for p in pts) / 64.0
    loy, hiy = min(p[1] for p in pts) / 64.0, max(p[1] for p in pts) / 64.0
    corners = [(lox, loy), (hix, loy), (hix, hiy), (lox, hiy)]
    pieces = []
    lines = 0
    for fam in families:
        ang, ox, oy, dx, dy, el = fam
        if dy == 0:
            continue
        ux, uy = unit_vector(int(ang))
        nx, ny = -uy, ux
        wv = (C[0] * ux + C[1] * uy, C[2] * ux + C[3] * uy)
        dvec = (H[0] * ux + H[1] * uy, H[2] * ux + H[3] * uy)
        for k in range(-400, 401):
            kdx = float(k) * dx
            kdy = float(k) * dy
            okx = (ox + kdx * ux) + kdy * nx
            oky = (oy + kdx * uy) + kdy * ny
            A = apply(C, (okx, oky))
            Af = (F(A[0]), F(A[1]))
            wf = (F(wv[0]), F(wv[1]))
            sg = []
            for c in corners:
                v = wf[0] * (F(c[1]) - Af[1]) - wf[1] * (F(c[0]) - Af[0])
                sg.append((v > 0) - (v < 0))
            if all(s > 0 for s in sg) or all(s < 0 for s in sg):
                continue
            lines += 1
            ivs = trim_intervals(contours, nesting, rule, Af, wf)
            if not ivs:
                continue
            ps = []
            if not el or all(e >= 0 for e in el):
                ps = [(a, b) for (a, b) in ivs if b > a]
            else:
                dsh, _ = dash_set(el, ivs[0][0] - 3 * sum(abs(F(e)) for e in el), ivs[-1][1] + 3 * sum(abs(F(e)) for e in el))
                for (da, db) in dsh:
                    if da == db:
                        # a dot is kept at an end to within 2^-46 of the larger of the period, 1 and the largest trim position
                        big = max(float(abs(ivs[0][0])), float(abs(ivs[-1][1])), 1.0, float(sum(abs(F(e)) for e in el)))
                        slack = F(big) / F(2 ** 46)
                        for (a, b) in ivs:
                            if a - slack <= da <= b + slack:
                                at = a if da < a else (b if da > b else da)
                                ps.append((at, at))
                                break
                    else:
                        for (a, b) in ivs:
                            lo, hi = max(da, a), min(db, b)
                            if hi > lo:
                                ps.append((lo, hi))
                ps.sort()
            for (s0, s1) in ps:
                p0 = apply(H, (fl(F(okx) + s0 * F(ux)), fl(F(oky) + s0 * F(uy))))
                p1 = apply(H, (fl(F(okx) + s1 * F(ux)), fl(F(oky) + s1 * F(uy))))
                pieces.append((1 if s0 == s1 else 0, p0, p1 if s1 != s0 else p0))
    return lines, pieces


# --- the matrix -----------------------------------------------------------------------------------------------------------
def sq(x0, y0, x1, y1, ccw=True):
    p = [(x0 * 64, y0 * 64), (x1 * 64, y0 * 64), (x1 * 64, y1 * 64), (x0 * 64, y1 * 64)]
    return p if ccw else p[::-1]


REGIONS = {
    'square': ([sq(0, 0, 10, 10)], 0),
    'square_eo': ([sq(0, 0, 10, 10)], 1),
    'lshape': ([[(0, 0), (12 * 64, 0), (12 * 64, 4 * 64), (4 * 64, 4 * 64), (4 * 64, 12 * 64), (0, 12 * 64)]], 0),
    'star': ([[(int(64 * x), int(64 * y)) for (x, y) in [(5, 0), (6.5, 3.5), (10, 4), (7.5, 6.5), (8, 10), (5, 8), (2, 10), (2.5, 6.5), (0, 4), (3.5, 3.5)]]], 0),
    'ring3': ([sq(0, 0, 16, 16), sq(3, 3, 13, 13, False), sq(6, 6, 10, 10)], 0),
    'ring3_eo': ([sq(0, 0, 16, 16), sq(3, 3, 13, 13), sq(6, 6, 10, 10)], 1),
    'ring3_same': ([sq(0, 0, 16, 16), sq(3, 3, 13, 13), sq(6, 6, 10, 10)], 0),
    'bowtie': ([[(0, 0), (10 * 64, 10 * 64), (10 * 64, 0), (0, 10 * 64)]], 0),
    'bowtie_eo': ([[(0, 0), (10 * 64, 10 * 64), (10 * 64, 0), (0, 10 * 64)]], 1),
    'tri': ([[(0, 0), (12 * 64, 0), (6 * 64, 9 * 64)]], 0),
    'two': ([sq(0, 0, 6, 6), sq(9, 2, 15, 8)], 0),
    'touch': ([sq(0, 0, 5, 5), sq(5, 5, 10, 10)], 0),
    'slit': ([sq(0, 0, 8, 8), [(8 * 64, 0), (8 * 64, 8 * 64), (0, 8 * 64), (0, 0)]], 0),
}
FAMILIES = {
    'solid0': [(0, 0, 0, 0, 2.0, [])],
    'solid45': [(45, 0, 0, 0, 1.5, [])],
    'solid90_o': [(90, 3.0, 2.0, 0, 2.5, [])],
    'solid30': [(30, 0.5, 0.25, 0, 1.0, [])],
    'edge0': [(0, 0, 0, 0, 5.0, [])],
    'dashed': [(0, 0, 0, 0, 2.0, [3.0, -1.0])],
    'dashed45': [(45, 0, 0, 0, 1.5, [2.0, -1.0, 0.5, -1.5])],
    'dotted': [(0, 0, 0, 0, 1.5, [0.0, -2.0])],
    'dotted_sd': [(90, 0, 0, 0, 2.0, [-2.0, 0.0])],
    'dotdash': [(0, 0, 0, 0, 2.0, [0.0, -1.0, 2.0])],
    'brick': [(0, 0, 0, 1.5, 2.0, [3.0, -1.0])],
    'brick_neg': [(0, 0, 0, -1.5, -2.0, [3.0, -1.0])],
    'cross': [(45, 0, 0, 0, 2.0, []), (135, 0, 0, 0, 2.0, [])],
    'cross_dash': [(0, 0, 0, 0, 2.5, [2.0, -1.0]), (90, 0, 0, 0, 2.5, [1.0, -1.0, 0.0, -1.0])],
    'degen_dy': [(0, 0, 0, 0, 0, []), (90, 0, 0, 0, 3.0, [])],
    'dashed270': [(270, 1.0, 1.0, 0.5, 1.75, [1.5, -0.5, 0.0, -0.5])],
}
MAPS = {
    'id': ((1.0, 0.0, 0.0, 1.0, 0.0, 0.0), (1.0, 0.0, 0.0, 1.0, 0.0, 0.0)),
    'scale2': ((2.0, 0.0, 0.0, 2.0, 5.0, 3.0), (1.0, 0.0, 0.0, 1.0, 0.0, 0.0)),
    'hscale': ((1.0, 0.0, 0.0, 1.0, 0.0, 0.0), (0.5, 0.0, 0.0, 0.5, 1.0, 1.0)),
    'flipy': ((1.0, 0.0, 0.0, -1.0, 0.0, 20.0), (1.0, 0.0, 0.0, 1.0, 0.0, 0.0)),
    'rot90': ((0.0, -1.0, 1.0, 0.0, 20.0, 0.0), (1.0, 0.0, 0.0, 1.0, 0.0, 0.0)),
    'nonuni': ((2.0, 0.0, 0.0, 1.0, 0.0, 0.0), (1.0, 0.0, 0.0, 2.0, 0.0, 0.0)),
}

PLAN = [
    # (region, family, nesting, map)
    ('square', 'solid0', 0, 'id'), ('square', 'solid45', 0, 'id'), ('square', 'solid90_o', 0, 'id'), ('square', 'solid30', 0, 'id'),
    ('square', 'edge0', 0, 'id'), ('square', 'dashed', 0, 'id'), ('square', 'dashed45', 0, 'id'), ('square', 'dotted', 0, 'id'),
    ('square', 'dotted_sd', 0, 'id'), ('square', 'dotdash', 0, 'id'), ('square', 'brick', 0, 'id'), ('square', 'brick_neg', 0, 'id'),
    ('square', 'cross', 0, 'id'), ('square', 'cross_dash', 0, 'id'), ('square', 'degen_dy', 0, 'id'), ('square', 'dashed270', 0, 'id'),
    ('lshape', 'solid45', 0, 'id'), ('lshape', 'dashed', 0, 'id'), ('lshape', 'edge0', 0, 'id'), ('lshape', 'brick', 0, 'id'),
    ('star', 'solid0', 0, 'id'), ('star', 'solid45', 0, 'id'), ('star', 'dashed45', 0, 'id'), ('star', 'dotted', 0, 'id'),
    ('star', 'solid90_o', 1, 'id'), ('star', 'solid30', 0, 'scale2'),
    ('ring3', 'solid0', 0, 'id'), ('ring3', 'dashed', 0, 'id'), ('ring3', 'dotted', 0, 'id'), ('ring3', 'solid45', 0, 'id'),
    ('ring3_eo', 'solid0', 0, 'id'), ('ring3_eo', 'dashed', 0, 'id'), ('ring3_same', 'solid0', 0, 'id'), ('ring3_same', 'solid0', 1, 'id'),
    ('ring3', 'solid0', 1, 'id'), ('ring3', 'solid0', 2, 'id'), ('ring3', 'solid90_o', 1, 'id'), ('ring3', 'solid90_o', 2, 'id'),
    ('ring3', 'dashed', 1, 'id'), ('ring3', 'dashed', 2, 'id'), ('ring3_same', 'dashed', 1, 'id'), ('ring3_same', 'dashed', 2, 'id'),
    ('bowtie', 'solid0', 0, 'id'), ('bowtie_eo', 'solid0', 0, 'id'), ('bowtie', 'dashed', 0, 'id'), ('bowtie', 'solid45', 0, 'id'),
    ('tri', 'solid0', 0, 'id'), ('tri', 'solid45', 0, 'id'), ('tri', 'dashed45', 0, 'id'), ('tri', 'cross', 0, 'id'),
    ('two', 'solid0', 0, 'id'), ('two', 'dashed', 0, 'id'), ('two', 'dotted', 0, 'id'), ('two', 'solid45', 0, 'id'),
    ('touch', 'solid0', 0, 'id'), ('touch', 'edge0', 0, 'id'), ('touch', 'solid90_o', 0, 'id'),
    ('slit', 'solid0', 0, 'id'), ('slit', 'solid90_o', 0, 'id'),
    ('square', 'solid0', 0, 'scale2'), ('square', 'dashed', 0, 'scale2'), ('square', 'solid45', 0, 'scale2'), ('square', 'dotted', 0, 'scale2'),
    ('square', 'solid0', 0, 'hscale'), ('square', 'dashed', 0, 'hscale'), ('square', 'brick', 0, 'hscale'),
    ('square', 'solid0', 0, 'flipy'), ('square', 'dashed45', 0, 'flipy'), ('square', 'brick', 0, 'flipy'),
    ('square', 'solid0', 0, 'rot90'), ('square', 'dashed', 0, 'rot90'),
    ('square', 'solid0', 0, 'nonuni'), ('square', 'dashed', 0, 'nonuni'), ('square', 'solid45', 0, 'nonuni'),
    ('star', 'dashed', 0, 'flipy'), ('ring3', 'dashed45', 0, 'scale2'), ('lshape', 'dotted', 0, 'hscale'),
]


def fmt_arr(vals):
    return '{' + ', '.join(hexf(v) for v in vals) + '}'


def main():
    out = []
    out.append('// Generated by tools/gen-hatch-vectors.py (exact rationals).  Do not edit; run the generator.')
    out.append('// Independent of the library\'s code.')
    out.append('#ifndef GHOTI_IO_GSHGEOM_TEST_HATCH_VECTORS_H')
    out.append('#define GHOTI_IO_GSHGEOM_TEST_HATCH_VECTORS_H')
    out.append('')
    out.append('#include <cstddef>')
    out.append('#include <cstdint>')
    out.append('')
    out.append('namespace hatch_vectors {')
    out.append('')
    out.append('// One piece: the dot flag and the two points in local space (the first twice for a dot).')
    out.append('struct VPiece {')
    out.append('  int dot;')
    out.append('  double p[4];')
    out.append('};')
    out.append('// One family: angle in degrees, origin, offset along, offset across, and the elements.')
    out.append('struct VFamily {')
    out.append('  double angle, ox, oy, dx, dy;')
    out.append('  int dash_count;')
    out.append('  const double * dashes;')
    out.append('};')
    out.append('// One case: the contours (26.6), the nesting, the rule (0 nonzero, 1 even-odd), the families, the hatch map and')
    out.append('// the chain {m00, m01, m10, m11, tx, ty}, the lines of the use and the pieces.')
    out.append('struct VCase {')
    out.append('  const char * name;')
    out.append('  int contour_count;')
    out.append('  const int * contour_sizes;')
    out.append('  const int32_t * points;')
    out.append('  int nesting, rule;')
    out.append('  int family_count;')
    out.append('  const VFamily * families;')
    out.append('  double hatch_map[6], chain[6];')
    out.append('  uint64_t lines;')
    out.append('  int piece_count;')
    out.append('  const VPiece * pieces;')
    out.append('};')
    out.append('')
    cases = []
    for idx, (rn, fn, nest, mn) in enumerate(PLAN):
        contours, rule = REGIONS[rn]
        fams = FAMILIES[fn]
        chain, H = MAPS[mn]
        # the vectors state H as the hatch map and chain as the chain
        lines, pieces = case_pieces(contours, nest, rule, 1, fams, H, chain)
        cases.append((idx, rn, fn, nest, mn, contours, rule, fams, H, chain, lines, pieces))
    for (idx, rn, fn, nest, mn, contours, rule, fams, H, chain, lines, pieces) in cases:
        out.append('static const int c%d_sizes[] = {%s};' % (idx, ', '.join(str(len(c)) for c in contours)))
        out.append('static const int32_t c%d_points[] = {%s};' % (idx, ', '.join('%d, %d' % p for c in contours for p in c)))
        for fi, f in enumerate(fams):
            if f[5]:
                out.append('static const double c%d_f%d_dashes[] = %s;' % (idx, fi, fmt_arr(f[5])))
        out.append('static const VFamily c%d_families[] = {' % idx)
        for fi, f in enumerate(fams):
            out.append('    {%s, %s, %s, %s, %s, %d, %s},' % (hexf(f[0]), hexf(f[1]), hexf(f[2]), hexf(f[3]), hexf(f[4]), len(f[5]),
                                                          'c%d_f%d_dashes' % (idx, fi) if f[5] else 'nullptr'))
        out.append('};')
        if pieces:
            out.append('static const VPiece c%d_pieces[] = {' % idx)
            for (dot, p0, p1) in pieces:
                out.append('    {%d, {%s, %s, %s, %s}},' % (dot, hexf(p0[0]), hexf(p0[1]), hexf(p1[0]), hexf(p1[1])))
            out.append('};')
    out.append('')
    out.append('static const VCase cases[] = {')
    for (idx, rn, fn, nest, mn, contours, rule, fams, H, chain, lines, pieces) in cases:
        out.append('    {"%s/%s/n%d/%s", %d, c%d_sizes, c%d_points, %d, %d, %d, c%d_families, %s, %s, %d, %d, %s},' % (
            rn, fn, nest, mn, len(contours), idx, idx, nest, rule, len(fams), idx, fmt_arr(H), fmt_arr(chain), lines, len(pieces),
            'c%d_pieces' % idx if pieces else 'nullptr'))
    out.append('};')
    out.append('static const size_t case_count = sizeof cases / sizeof cases[0];')
    out.append('')
    out.append('}  // namespace hatch_vectors')
    out.append('')
    out.append('#endif')
    sys.stdout.write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
