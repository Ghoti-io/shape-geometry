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
"""Print the exact vectors of the dash rules (stroke.h, SG21 and SG30).

An independent statement in exact rational arithmetic (Python's Fraction, no
library code): the dash positions of the SVG 2 Candidate Recommendation of
2018-10-04 (painting chapter, "dash positions"), transcribed literally with a
running position, and the lowering of a DXF simple linetype, checked by the
definition of the lowering rather than by a second copy of its algorithm: the
pen-down intervals of the raw pattern (dash positive, space negative, dot zero,
repeated from the start of the line) must equal those of the lowered array and
offset under the positions model.

Every input is a dyadic rational small enough that the library's double
arithmetic is exact, so the vectors are compared bit for bit.  Needs only the
standard library:

    python3 tools/gen-dash-vectors.py > tests/dash_vectors.h
"""
from fractions import Fraction as F
import sys


def hexd(x):
    """A Fraction that is exactly a double, as a C hexadecimal literal."""
    v = float(x)
    assert F(v) == x, "not a double: %s" % x
    return v.hex()


def positions(dashes, offset, length):
    """SVG 2's dash positions, literally.  dashes: the canonical, even-length, positive-sum list."""
    count = len(dashes)
    total = sum(dashes)
    if total == 0:
        return [(F(0), length)]
    out = []
    if offset < 0:
        offset = total - abs(offset)
    offset = offset % total  # Python's % is the floored modulo
    index = 0
    while sum(dashes[: index + 1]) < offset:
        index += 1
    remaining = sum(dashes[: index + 1]) - offset
    dashlength = min(remaining, length)
    # half open (stroke.h SG21, F3): a positive dash the phase has used up is over; the CR's loop would make a dash of
    # no length here.  A dash of no length exactly at the phase is kept.
    if index % 2 == 0 and (remaining > 0 or dashes[index] == 0):
        out.append((F(0), dashlength))
    position = dashlength
    while position < length:
        index = (index + 1) % count
        dashlength = min(dashes[index], length - position)
        if index % 2 == 0:
            out.append((position, position + dashlength))
        position += dashlength
    return out


def canonical(dashes):
    if any(d < 0 for d in dashes) or all(d == 0 for d in dashes):
        return []
    return list(dashes) * 2 if len(dashes) % 2 else list(dashes)


ARRAYS = [
    [5, 3],
    [5, 3, 2],
    [0, 10],
    [10, 0],
    [3, 0, 2, 5],
    [1, 1],
    [F(1, 2), F(1, 4), F(1, 8)],
    [4, 2, 0, 2],
    [7],
    [0, 0],
    [1, -1],
    [],
    [2, 0, 0, 1],
]
OFFSETS = [0, 1, 5, 8, -3, -13, 10, 26, F(1, 2), F(-1, 4), -8, -10, 32, 2**40 + 3]
LENGTHS = [0, 1, 7, 10, 30, 100, F(21, 2)]


def pos_cases():
    cases = []
    for arr in ARRAYS:
        arr = [F(a) for a in arr]
        canon = canonical(arr)
        for off in OFFSETS:
            for length in LENGTHS:
                off, length = F(off), F(length)
                if canon:
                    pairs = positions(canon, off, length)
                else:
                    pairs = [(F(0), length)]
                cases.append((arr, off, length, pairs))
    return cases


# --- DXF ------------------------------------------------------------------------------

DXF = [
    # (name, elements, scale)
    ("dotted", [F(0), F(-1, 4)], F(1)),
    ("dotted_double", [F(0), F(-1, 2)], F(1)),
    ("dashed", [F(1, 2), F(-1, 4)], F(1)),
    ("dashed_scaled", [F(1, 2), F(-1, 4)], F(4)),
    ("dashdot", [F(1, 2), F(-1, 4), F(0), F(-1, 4)], F(1)),
    ("divide", [F(1, 2), F(-1, 4), F(0), F(-1, 4), F(0), F(-1, 4)], F(1)),
    ("center", [F(5, 4), F(-1, 4), F(1, 4), F(-1, 4)], F(1)),
    ("phantom", [F(5, 4), F(-1, 4), F(1, 4), F(-1, 4), F(1, 4), F(-1, 4)], F(1)),
    ("border", [F(1, 2), F(-1, 4), F(1, 2), F(-1, 4), F(0), F(-1, 4)], F(1)),
    ("space_first", [F(-1), F(2)], F(1)),
    ("space_first_last", [F(-2), F(5), F(-3)], F(1)),
    ("dash_last", [F(5), F(-2), F(1)], F(1)),
    ("dot_first_dash_last", [F(0), F(-2), F(3)], F(1)),
    ("dash_dash", [F(2), F(1), F(-1), F(2)], F(1)),
    ("space_space", [F(2), F(-1), F(-1), F(2), F(-2)], F(1)),
    ("dot_dash", [F(0), F(2), F(-1)], F(1)),
    ("solid", [F(3), F(1)], F(1)),
    ("space_first_dot_last", [F(-1), F(2), F(-3), F(0)], F(1)),
    ("space_first_dash_last", [F(-1), F(2), F(1)], F(1)),
    ("spaces_merge_dot_last", [F(-1), F(2), F(-3), F(0), F(-2)], F(1)),
    ("dots_only_spaced", [F(0), F(-2), F(0), F(-4)], F(1)),
]


def normalize(intervals, span):
    """Merge touching intervals, clip to [0, span]; a zero-length interval absorbed by a neighbour disappears."""
    iv = sorted((a, min(b, span)) for a, b in intervals if a <= span)
    out = []
    for a, b in iv:
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def lower(vals):
    """The lowering (stroke.h SG30), in Fractions: (dash list, offset), None for solid, "unsupported" for a pattern
    that begins with a space no trailing space continues and ends with a dot."""
    down = [v >= 0 for v in vals]
    if all(down):
        return None
    runs = []  # [is_down, length]
    for v, d in zip(vals, down):
        if runs and runs[-1][0] == d:
            runs[-1][1] += abs(v)
        else:
            runs.append([d, abs(v)])
    start_off = F(0)
    if runs[0][0] and runs[0][1] == 0 and runs[-1][0]:
        runs.append([False, F(0)])  # a dot first, a dash last: kept apart so the dot is not lost
    elif len(runs) > 1 and runs[0][0] == runs[-1][0]:
        last = runs.pop()
        runs[0][1] += last[1]
        start_off = last[1]
    offset = start_off
    if not runs[0][0]:
        if start_off == 0 and runs[-1][1] == 0:
            return "unsupported"  # the phase would stand exactly at the dot that ends the period (SG21 keeps it)
        first = runs.pop(0)
        runs.append(first)
        period = sum(r[1] for r in runs)
        offset = start_off - first[1]
        if offset < 0:
            offset += period
    return [r[1] for r in runs], offset


def dxf_cases():
    cases = []
    for name, elements, scale in DXF:
        vals = [e * scale for e in elements]
        res = lower(vals)
        period = sum(abs(v) for v in vals)
        span = period * 3 + F(1, 64)
        raw = []
        x = F(0)
        while x < span:
            for v in vals:
                if v >= 0:
                    raw.append((x, x + v))
                x += abs(v)
        want = normalize(raw, span)
        if res == "unsupported":
            cases.append((name, elements, scale, "unsupported", F(0)))
            continue
        if res is None:
            got = [(F(0), span)]
            # a solid pattern has no spaces: every interval touches the next
            assert len(want) == 1 and want[0][0] == 0 and want[0][1] >= span - 0, name
            cases.append((name, elements, scale, None, F(0)))
            continue
        dashes, offset = res
        got = normalize(positions(dashes, offset, span), span)
        # the raw pattern's last interval may run past span: compare clipped
        assert got == want, (name, got, want)
        cases.append((name, elements, scale, dashes, offset))
    return cases


def arr(name, values):
    if not values:
        return "static const double %s[1] = {0.0};" % name
    return "static const double %s[%d] = {%s};" % (name, len(values), ", ".join(hexd(F(v)) for v in values))


def main():
    out = []
    out.append("/* Generated by tools/gen-dash-vectors.py (exact rationals). Do not edit. */")
    out.append("#ifndef GHOTI_IO_GSHGEOM_TEST_DASH_VECTORS_H")
    out.append("#define GHOTI_IO_GSHGEOM_TEST_DASH_VECTORS_H")
    out.append("")
    out.append("#include <cstddef>")
    out.append("")
    out.append("namespace dash_vectors {")
    out.append("")
    out.append("struct Positions {")
    out.append("  const double * dashes; size_t count; double offset; double length;")
    out.append("  const double (*pairs)[2]; size_t pair_count;")
    out.append("};")
    out.append("")
    cases = pos_cases()
    for i, (arr_, off, length, pairs) in enumerate(cases):
        out.append(arr("d%d" % i, arr_))
        if pairs:
            out.append("static const double p%d[%d][2] = {%s};" % (
                i, len(pairs), ", ".join("{%s, %s}" % (hexd(a), hexd(b)) for a, b in pairs)))
        else:
            out.append("static const double p%d[1][2] = {{0.0, 0.0}};" % i)
    out.append("static const Positions positions[] = {")
    for i, (arr_, off, length, pairs) in enumerate(cases):
        out.append("  {d%d, %d, %s, %s, p%d, %d}," % (i, len(arr_), hexd(off), hexd(length), i, len(pairs)))
    out.append("};")
    out.append("static const size_t positions_count = %d;" % len(cases))
    out.append("")
    out.append("struct Dxf {")
    out.append("  const char * name; const double * elements; size_t count; double scale;")
    out.append("  const double * dashes; size_t dash_count; double offset; int kind; /* 0 lowered, 1 solid, 2 unsupported */")
    out.append("};")
    dc = dxf_cases()
    for i, (name, elements, scale, dashes, offset) in enumerate(dc):
        out.append(arr("xe%d" % i, elements))
        out.append(arr("xd%d" % i, dashes if isinstance(dashes, list) else []))
    out.append("static const Dxf dxf[] = {")
    for i, (name, elements, scale, dashes, offset) in enumerate(dc):
        kind = 2 if dashes == "unsupported" else (1 if dashes is None else 0)
        out.append('  {"%s", xe%d, %d, %s, xd%d, %d, %s, %d},' % (
            name, i, len(elements), hexd(scale), i, len(dashes) if kind == 0 else 0, hexd(offset), kind))
    out.append("};")
    out.append("static const size_t dxf_count = %d;" % len(dc))
    out.append("")
    out.append("} // namespace dash_vectors")
    out.append("")
    out.append("#endif")
    sys.stdout.write("\n".join(out) + "\n")


main()
