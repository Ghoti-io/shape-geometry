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
"""Write the seed corpora of the path, stroke, clip and hatch fuzzers, deterministically (story 16; CAP-15, AD-19).

    python3 tools/gen-fuzz-corpus.py            # write tests/fuzz/corpus/path/, stroke/ and clip/
    python3 tools/gen-fuzz-corpus.py --check    # fail if the committed files differ

The stroke corpus (tests/fuzz/fuzz_stroke.c) is described at stroke_files() below, and the clip corpus
(tests/fuzz/fuzz_clip.c) at clip_files().

Every file is an input of tests/fuzz/fuzz_path.c: an option byte (the dimension, the
limits, whether an allocation fails), the allocation to fail, then a stream of ops,
each an op byte and argument bytes.  The **biased** files choose ops by a weight
table that favours curves and the passes that read them (the device path, the
flattener, the measured path) and argument bytes that lean on the table of degenerate
coordinates; the **plain** files are random bytes, which the fuzzer reads as valid
input by construction; the **scripted** files reach, by construction, each outcome
the replay requires.  The replay gate (`make fuzz-replay`) fails if these files,
together, never reach one of those outcomes: that is what keeps the corpus from
measuring nothing, and what to read if an op number below drifts from the C file.
"""
import base64
import os
import sys

# Keep in step with the op enum in tests/fuzz/fuzz_path.c.
OPS = ["MOVE", "LINE", "QUAD", "CUBIC", "RQUAD", "RCUBIC", "CLOSE", "MOVE_FIXED", "LINE_FIXED", "QUAD_FIXED",
       "CUBIC_FIXED", "DEVICE", "MEASURE", "ORIENT", "BOUNDS", "ARC_ENDPOINT", "ARC_CENTER", "ARC_CONJUGATE", "ARC_TANGENT",
       "ARC_BULGE", "ARC_RAW", "ARC_PASS", "ARC_MAP", "SHAPE"]
OP = {name: i for i, name in enumerate(OPS)}

BIASED = {"MOVE": 3, "LINE": 5, "QUAD": 5, "CUBIC": 6, "RQUAD": 4, "RCUBIC": 3, "CLOSE": 3, "MOVE_FIXED": 1,
          "LINE_FIXED": 2, "QUAD_FIXED": 2, "CUBIC_FIXED": 2, "DEVICE": 6, "MEASURE": 4, "ORIENT": 3, "BOUNDS": 1,
          "ARC_ENDPOINT": 4, "ARC_CENTER": 3, "ARC_CONJUGATE": 2, "ARC_TANGENT": 3, "ARC_BULGE": 3, "ARC_RAW": 3,
          "ARC_PASS": 3, "ARC_MAP": 3, "SHAPE": 4}

BIASED_FILES = 120
PLAIN_FILES = 40

# A coordinate byte below 24 picks the table in fuzz_path.c (21 NaN, 22 +inf, 20 and 23 -inf, 9 is 2^25: past
# the 26.6 range under the identity, 2 is 1.0); 24 and above takes the next two bytes as a 16-bit value in 1/16.
NAN, ONE, PAST_26_6 = 21, 2, 9


class Rng:
    def __init__(self, seed):
        self.s = (seed * 0x9E3779B97F4A7C15 + 54321) & (2**64 - 1)

    def next(self):
        s = self.s
        s ^= s >> 12
        s ^= (s << 25) & (2**64 - 1)
        s ^= s >> 27
        self.s = s
        return (s * 2685821657736338717) & (2**64 - 1)

    def below(self, n):
        return self.next() % n

    def byte(self):
        return self.next() & 255


def raw(value16):
    """A coordinate of value16 / 16 (value16 in [-32768, 32767])."""
    v = value16 & 0xFFFF
    return bytes([24, v >> 8, v & 255])


def biased(seed):
    rng = Rng(seed)
    table = []
    for name, weight in BIASED.items():
        table += [OP[name]] * weight
    # Mostly 2D, the default limits and no failing allocation; sometimes each of the others.
    options = 0
    if rng.below(4) == 0:
        options |= 1
    if rng.below(6) == 0:
        options |= 2 | (rng.below(16) << 2)
    if rng.below(6) == 0:
        options |= 64
    if rng.below(8) == 0:
        options |= 128
    out = bytearray([options, rng.byte()])
    for _ in range(10 + rng.below(60)):
        out.append(table[rng.below(len(table))])
        for _ in range(2 + rng.below(20)):
            # lean on the degenerate table one time in four, else ordinary coordinates
            out += bytes([rng.below(24)]) if rng.below(4) == 0 else raw(rng.below(4096) - 2048)
    return bytes(out)


def plain(seed):
    rng = Rng(1000 + seed)
    return bytes(rng.byte() for _ in range(2 + 50 + rng.below(400)))


def identity_device(which=0):
    # DEVICE: no map steps (the identity), then the flattener's pick (below 128: the default tolerance)
    return bytes([OP["DEVICE"], 0, which])


def cubic():
    return (bytes([OP["MOVE"]]) + raw(0) + raw(0) + bytes([OP["CUBIC"]]) + raw(160) + raw(0) + raw(320) + raw(480)
            + raw(0) + raw(640))


# Arc ops (tests/fuzz/fuzz_path.c).  A radius is a byte below 16 (the table) or a byte of 16 and up and two more; an angle a
# byte below 24 or one of 24 and up and two more; raw() writes 24, which both read as a value.
A_ZERO, A_DEG = 0, 1


def arc_endpoint(e, rx, ry, phi=bytes([0]), unit=A_DEG, flags=2):
    return bytes([OP["ARC_ENDPOINT"]]) + e + rx + ry + phi + bytes([unit, flags])


def arc_center(c, rx, ry, phi, unit, t0, sweep, end=None):
    return (bytes([OP["ARC_CENTER"]]) + c + rx + ry + phi + bytes([unit]) + t0 + sweep + bytes([0 if end is None else 1])
            + (end or raw(0) + raw(0)))


def arc_conj(c, u, v, unit, t0, sweep, flags=0, end=None):
    return (bytes([OP["ARC_CONJUGATE"]]) + c + u + v + bytes([unit]) + t0 + sweep + bytes([flags, 0 if end is None else 1])
            + (end or raw(0) + raw(0)))


def pt2(x, y):
    return raw(x) + raw(y)


def move0():
    return bytes([OP["MOVE"]]) + pt2(0, 0)


ARC_PASS = bytes([OP["ARC_PASS"]])
SEMI = move0() + arc_endpoint(pt2(320, 0), raw(160), raw(160))


def arc_scripted():
    """Arc inputs: each constructor, each kind of result, the agreement bound, the readers, the map, the limits."""
    out = {}
    out["scripted-arc-endpoint.bin"] = bytes([0, 0]) + SEMI + ARC_PASS
    out["scripted-arc-nothing.bin"] = bytes([0, 0]) + move0() + arc_endpoint(pt2(0, 0), raw(160), raw(160)) + ARC_PASS
    out["scripted-arc-line.bin"] = bytes([0, 0]) + move0() + arc_endpoint(pt2(320, 0), bytes([A_ZERO]), raw(160)) + ARC_PASS
    out["scripted-arc-endpoint-flag.bin"] = bytes([0, 0]) + move0() + arc_endpoint(pt2(320, 0), raw(160), raw(160), flags=0xF0)
    out["scripted-arc-ellipse.bin"] = bytes([0, 0]) + move0() + arc_endpoint(pt2(160, 160), raw(320), raw(96), raw(30 * 16), flags=1) + ARC_PASS
    # a full-turn circle by the centre form, in degrees (angle table index 4 is 360), beginning a subpath of its own
    out["scripted-arc-full-turn.bin"] = bytes([0, 0]) + arc_center(pt2(0, 0), raw(160), raw(160), bytes([0]), A_DEG, bytes([0]),
                                                                    bytes([4])) + ARC_PASS + identity_device()
    # a conjugate form with a negative sweep and an explicit end, and an elliptical arc in turns
    out["scripted-arc-conjugate.bin"] = bytes([0, 0]) + arc_conj(pt2(0, 0), pt2(160, 0), pt2(0, 80), A_DEG, bytes([0]), bytes([5])) + ARC_PASS
    out["scripted-arc-conjugate-flag.bin"] = bytes([0, 0]) + arc_conj(pt2(0, 0), pt2(160, 0), pt2(0, 160), A_DEG, bytes([0]), bytes([2]), flags=1) + ARC_PASS
    out["scripted-arc-conjugate-end.bin"] = bytes([0, 0]) + arc_conj(pt2(0, 0), pt2(160, 0), pt2(0, 160), 3, bytes([0]), bytes([13]), end=pt2(0, 160)) + ARC_PASS
    out["scripted-arc-tangent.bin"] = bytes([0, 0]) + move0() + bytes([OP["ARC_TANGENT"]]) + pt2(160, 0) + pt2(160, 160) + raw(32) + ARC_PASS
    out["scripted-arc-tangent-line.bin"] = bytes([0, 0]) + move0() + bytes([OP["ARC_TANGENT"]]) + pt2(160, 0) + pt2(320, 0) + raw(32)
    out["scripted-arc-bulge.bin"] = bytes([0, 0]) + move0() + bytes([OP["ARC_BULGE"]]) + pt2(320, 0) + bytes([1]) + ARC_PASS
    out["scripted-arc-bulge-infinite.bin"] = bytes([0, 0]) + move0() + bytes([OP["ARC_BULGE"]]) + pt2(320, 0) + bytes([10])
    # a raw record: the end as the form gives it, a nudge inside the agreement bound (nudge 2), and past it (nudge 1: refused)
    def raw_arc(nudge):
        return (bytes([OP["ARC_RAW"]]) + pt2(0, 0) + pt2(160, 0) + pt2(0, 160) + bytes([A_DEG]) + bytes([0]) + bytes([1])
                + bytes([1, nudge]))
    out["scripted-arc-raw.bin"] = bytes([0, 0]) + raw_arc(0) + ARC_PASS
    out["scripted-arc-raw-nudged.bin"] = bytes([0, 0]) + raw_arc(2) + ARC_PASS
    out["scripted-arc-raw-disagree.bin"] = bytes([0, 0]) + raw_arc(1) + ARC_PASS
    # the map: the identity, a scale by 0 and 1 (singular: refused), a rotation by 30 degrees, a non-uniform scale
    out["scripted-arc-map.bin"] = bytes([0, 0]) + SEMI + bytes([OP["ARC_MAP"], 0])
    out["scripted-arc-map-singular.bin"] = bytes([0, 0]) + SEMI + bytes([OP["ARC_MAP"], 1, 1, 0, 2])
    out["scripted-arc-map-rotate.bin"] = bytes([0, 0]) + SEMI + bytes([OP["ARC_MAP"], 1, 2, 1]) + raw(480)
    out["scripted-arc-map-nonuniform.bin"] = bytes([0, 0]) + SEMI + bytes([OP["ARC_MAP"], 1, 1]) + raw(32) + raw(16)
    # an arc in a 3D path, the element cap, a point cap, a failed allocation at several places
    out["scripted-arc-3d.bin"] = bytes([1, 0]) + move0() + arc_endpoint(pt2(320, 0), raw(160), raw(160))
    out["scripted-arc-element-cap.bin"] = bytes([2 | (1 << 2), 0]) + SEMI + arc_endpoint(pt2(640, 0), raw(160), raw(160)) + ARC_PASS
    for k, fail in enumerate((1, 3, 5, 7, 9, 12)):
        out["scripted-arc-oom-%d.bin" % k] = bytes([128, fail]) + SEMI + ARC_PASS + identity_device() + bytes([OP["MEASURE"], 0])
    out["scripted-arc-device.bin"] = bytes([0, 0]) + SEMI + identity_device() + bytes([OP["MEASURE"], 0, 10, 100, 200, 255, ONE])
    return out


def shape(kind, flags, x, y, w, h, radii=()):
    """OP_SHAPE: the kind (0 rect, 1 rounded rect, 2 ellipse), the flip byte (bit 0 flip_x, bit 1 flip_y, 0x10 a flip that is not 0 or 1), the
    box and, for a rounded rect, the eight radii (a byte below 16 is the table: 0 is a zero radius)."""
    out = bytes([OP["SHAPE"], kind, flags]) + raw(x) + raw(y) + raw(w) + raw(h)
    for r in radii:
        out += raw(r) if r else bytes([0])
    return out


def shape_scripted():
    """Primitive inputs: each kind, each flip pair, an overlap, a bad flip, a 3D path, the element cap and failed allocations."""
    out = {}
    for fl in range(4):
        out["scripted-shape-rect-%d.bin" % fl] = bytes([0, 0]) + shape(0, fl, 160, 320, 1600, 960)
        out["scripted-shape-ellipse-%d.bin" % fl] = bytes([0, 0]) + shape(2, fl, 160, 320, 1600, 960)
    rr = (160, 96, 160, 96, 160, 96, 160, 96)
    out["scripted-shape-rounded.bin"] = bytes([0, 0]) + shape(1, 0, 0, 0, 1600, 960, rr) + ARC_PASS + identity_device()
    out["scripted-shape-rounded-two.bin"] = bytes([0, 0]) + shape(1, 0, 0, 0, 1600, 960, (160, 96, 0, 0, 0, 0, 96, 160)) + ARC_PASS
    for fl in (1, 2, 3):
        out["scripted-shape-rounded-flipped-%d.bin" % fl] = bytes([0, 0]) + shape(1, fl, 0, 0, 1600, 960, (160, 96, 0, 0, 0, 0, 96, 160)) + ARC_PASS
    out["scripted-shape-overlap.bin"] = bytes([0, 0]) + shape(1, 0, 0, 0, 1600, 960, (1000, 96, 1000, 96, 0, 0, 0, 0))
    out["scripted-shape-bad-flip.bin"] = bytes([0, 0]) + shape(0, 0x10, 0, 0, 1600, 960) + shape(2, 0x30, 0, 0, 1600, 960)
    out["scripted-shape-3d.bin"] = bytes([1, 0]) + shape(0, 0, 0, 0, 1600, 960)
    out["scripted-shape-element-cap.bin"] = bytes([2 | (3 << 2), 0]) + shape(0, 0, 0, 0, 1600, 960) + shape(2, 0, 0, 0, 1600, 960)
    for k, fail in enumerate((1, 2, 4, 6, 9)):
        out["scripted-shape-oom-%d.bin" % k] = bytes([128, fail]) + shape(1, 0, 0, 0, 1600, 960, rr) + shape(2, 3, 0, 0, 800, 800) + ARC_PASS
    return out


def scripted():
    """Inputs that reach, by construction, the outcomes that random bytes reach only by luck."""
    out = {}
    out.update(arc_scripted())
    out.update(shape_scripted())
    out["scripted-nan.bin"] = bytes([0, 0, OP["MOVE"], NAN, 0])
    # options 2: max_elements 1, so the line after the move is refused by the cap.
    out["scripted-element-cap.bin"] = bytes([2, 0, OP["MOVE"], 0, 0, OP["LINE"], ONE, ONE])
    out["scripted-curve.bin"] = bytes([0, 0]) + cubic() + identity_device() + bytes([OP["MEASURE"], 0, 10, 100, 200,
                                                                                       255, ONE])
    out["scripted-device-range.bin"] = bytes([0, 0, OP["MOVE"], PAST_26_6, 0, OP["LINE"], 0, 0]) + identity_device()
    # options 64: max_points 2, which the cubic's flattening passes.
    out["scripted-point-cap.bin"] = bytes([64, 0]) + cubic() + identity_device()
    # A flattening pick of 128 asks for a tolerance of -2 units: below the minimum, refused.
    out["scripted-tolerance.bin"] = bytes([0, 0]) + cubic() + identity_device(128)
    # MEASURE's tolerance picks 4 (2^-47, below the smallest) and 6 (negative) are refused.
    out["scripted-measure-refused.bin"] = bytes([0, 0]) + cubic() + bytes([OP["MEASURE"], 4, OP["MEASURE"], 6])
    # options 128: the allocation numbered by the second byte fails; 0 is the path's own.
    out["scripted-oom.bin"] = bytes([128, 0, OP["MOVE"], 0, 0])
    out["scripted-oom-later.bin"] = bytes([128, 3]) + cubic() + identity_device() + bytes([OP["MEASURE"], 0])
    # A design-to-model map (scale 1, no translation) is the one the grid predicate accepts.
    out["scripted-grid.bin"] = bytes([0, 0]) + cubic() + bytes([OP["DEVICE"], 1, 4, ONE, 0, 0, 0])
    # A rotation in an angle unit that does not exist (4) is an invalid map.
    out["scripted-map-invalid.bin"] = bytes([0, 0]) + cubic() + bytes([OP["DEVICE"], 1, 2, 4, 0])
    # Four points at the origin: every orientation is 0.
    out["scripted-orient-zero.bin"] = bytes([0, 0, OP["ORIENT"]] + [0] * 12)
    return out


def files():
    out = {}
    for i in range(BIASED_FILES):
        out["biased-%03d.bin" % i] = biased(i + 1)
    for i in range(PLAIN_FILES):
        out["plain-%03d.bin" % i] = plain(i + 1)
    out.update(scripted())
    return out


# --- the stroke fuzzer (tests/fuzz/fuzz_stroke.c) -------------------------------------------------------------
# An input: options (bits 0-1 the limits, 2 a failing allocation, 3 the tolerance from the third byte), the
# allocation to fail, the tolerance; the parameters (width, width space, join, miter limit, cap, degenerate,
# alignment, alignment rule, a dash mode); the map (a step count and the steps); a byte that makes the map 3D when
# 9; the path (a byte that makes it 3D when 7, an op count, the ops: 0 move, 1 line, 2 quadratic, 3 cubic,
# 4 rational quadratic, 5 rational cubic, 6 close).  A coordinate byte below 20 picks the C file's table
# (0 is 0, 2 is 1, 16 is 3, 17 is 8, 18 NaN); 24 and up takes two more bytes, as raw() writes.

S_ZERO, S_ONE = 0, 2


def s_params(width=raw(32), space=0, join=0, limit=raw(64), cap=0, degenerate=0, align=0, rule=0, dash=b"\x00"):
    return width + bytes([space, join]) + limit + bytes([cap, degenerate, align, rule]) + dash


def s_path(ops):
    out = bytearray([0, len(ops)])
    for op in ops:
        out += op
    return bytes(out)


def s_pt(x, y):
    return raw(x) + raw(y)


S_SQUARE = [bytes([0]) + s_pt(0, 0), bytes([1]) + s_pt(128, 0), bytes([1]) + s_pt(128, 128), bytes([1]) + s_pt(0, 128),
            bytes([6])]
S_POLY = [bytes([0]) + s_pt(0, 0), bytes([1]) + s_pt(160, 0), bytes([1]) + s_pt(80, 96)]
S_CURVE = [bytes([0]) + s_pt(0, 0), bytes([2]) + s_pt(80, 160) + s_pt(160, 0), bytes([1]) + s_pt(240, 64)]


def s_input(params, ops, opts=0, fail=0, tol=4, map_bytes=b"\x00", map3d=0, tail=b""):
    return bytes([opts, fail, tol]) + params + map_bytes + bytes([map3d]) + s_path(ops) + tail


def s_arc_endpoint(e, rx, ry, phi=bytes([0]), flags=2):
    return bytes([7]) + e + rx + ry + phi + bytes([flags])


S_ARC = [bytes([0]) + s_pt(0, 0), s_arc_endpoint(s_pt(320, 0), raw(160), raw(160))]
S_ARC_ELLIPSE = [bytes([0]) + s_pt(0, 0), s_arc_endpoint(s_pt(160, 160), raw(320), raw(96), raw(30 * 16), 1)]
S_ARC_CIRCLE = [bytes([8]) + s_pt(160, 0) + raw(160) + raw(160) + bytes([0, 1, 0, 8])]
S_ARC_ROUNDED = [bytes([0]) + s_pt(0, 0), bytes([1]) + s_pt(240, 0), s_arc_endpoint(s_pt(320, 80), raw(80), raw(80)),
                 bytes([1]) + s_pt(320, 160), s_arc_endpoint(s_pt(240, 240), raw(80), raw(80)), bytes([6])]
S_ARC_TANGENT = [bytes([0]) + s_pt(0, 0), bytes([9]) + s_pt(160, 0) + s_pt(160, 160) + raw(32)]
S_ARC_BULGE = [bytes([0]) + s_pt(0, 0), bytes([10]) + s_pt(320, 0) + raw(64)]


def stroke_arc_scripted():
    out = {}
    out["s-arc.bin"] = s_input(s_params(), S_ARC)
    out["s-arc-ellipse.bin"] = s_input(s_params(join=1, cap=1), S_ARC_ELLIPSE)
    out["s-arc-circle.bin"] = s_input(s_params(), S_ARC_CIRCLE)
    out["s-arc-rounded.bin"] = s_input(s_params(join=1), S_ARC_ROUNDED)
    out["s-arc-tangent.bin"] = s_input(s_params(join=3, limit=raw(64)), S_ARC_TANGENT)
    out["s-arc-bulge.bin"] = s_input(s_params(cap=2), S_ARC_BULGE)
    out["s-arc-arcs-join.bin"] = s_input(s_params(join=4, limit=raw(160)), S_ARC_ROUNDED)
    out["s-arc-dashed.bin"] = s_input(s_params(cap=1, dash=bytes([6, 0, 5, 3, 0, 0])), S_ARC)
    out["s-arc-dashed-closed.bin"] = s_input(s_params(join=1, dash=bytes([6, 1, 4, 3, 5, 0, 0])), S_ARC_ROUNDED)
    out["s-arc-dashed-device.bin"] = s_input(s_params(space=1, cap=2, dash=bytes([6, 0, 5, 3, 3, 0])), S_ARC_ELLIPSE,
                                             map_bytes=bytes([1, 1]) + raw(48) + raw(16))
    out["s-arc-device.bin"] = s_input(s_params(space=1), S_ARC_CIRCLE, map_bytes=bytes([1, 1]) + raw(48) + raw(16))
    out["s-arc-inside.bin"] = s_input(s_params(align=1, rule=1), S_ARC_ROUNDED)
    out["s-arc-outside.bin"] = s_input(s_params(align=2), S_ARC_CIRCLE)
    out["s-arc-curves.bin"] = s_input(s_params(join=1, cap=1), S_ARC_ELLIPSE, tail=bytes([0, 0]))
    out["s-arc-curves-cubic.bin"] = s_input(s_params(join=1, cap=1), S_ARC_ROUNDED, tail=bytes([0, 1]))
    out["s-arc-curves-aligned.bin"] = s_input(s_params(align=1, rule=1), S_ARC_ROUNDED, tail=bytes([1, 0]))
    out["s-arc-rotated.bin"] = s_input(s_params(), S_ARC_ELLIPSE, map_bytes=bytes([1, 2, 1]) + raw(480))
    out["s-arc-oom.bin"] = s_input(s_params(join=1, cap=1), S_ARC_ELLIPSE, opts=4, fail=9)
    out["s-arc-limit.bin"] = s_input(s_params(join=1), S_ARC_ROUNDED, opts=2)
    return out


def stroke_scripted():
    out = stroke_arc_scripted()
    out["s-drawn.bin"] = s_input(s_params(), S_POLY)
    out["s-empty.bin"] = s_input(s_params(width=bytes([S_ZERO])), S_POLY)
    # a map step of kind 1 (a scale) by 0 and 1: singular
    out["s-singular.bin"] = s_input(s_params(), S_POLY, map_bytes=bytes([1, 1, S_ZERO, S_ONE]))
    out["s-invalid.bin"] = s_input(s_params(space=2), S_POLY)
    out["s-map-3d.bin"] = s_input(s_params(), S_POLY, map3d=9)
    out["s-dashes.bin"] = s_input(s_params(dash=bytes([1]) + raw(32) + raw(16)), S_POLY)
    # dash mode 6: a count byte (2 + b % 3), the entries and the offset as indexes into the C file's table (0 is 0, 2 is
    # 1, 3 is 2, 4 is 3, 5 is 4, 7 is 8), a scope byte
    out["s-dashed.bin"] = s_input(s_params(dash=bytes([6, 0, 5, 3, 0, 0])), S_POLY)
    out["s-dashed-path-scope.bin"] = s_input(
        s_params(cap=1, dash=bytes([6, 0, 5, 3, 2, 1])),
        S_POLY + [bytes([0]) + s_pt(0, 96), bytes([1]) + s_pt(100, 96)])
    out["s-dashed-dots.bin"] = s_input(s_params(cap=1, dash=bytes([6, 0, 0, 5, 0, 0])), S_POLY)
    out["s-dashed-closed.bin"] = s_input(s_params(join=1, dash=bytes([6, 1, 4, 3, 5, 0, 0])), S_SQUARE)
    out["s-dashed-device.bin"] = s_input(s_params(space=1, cap=2, dash=bytes([6, 0, 5, 3, 3, 0])), S_CURVE,
                                         map_bytes=bytes([1, 1]) + raw(48) + raw(16))
    out["s-dashed-inside.bin"] = s_input(s_params(align=1, dash=bytes([6, 0, 5, 3, 0, 0])), S_SQUARE)
    out["s-dashed-bad-scope.bin"] = s_input(s_params(dash=bytes([6, 0, 5, 3, 0, 2])), S_POLY)
    out["s-profile.bin"] = s_input(s_params(dash=bytes([4])), S_POLY)
    out["s-limit.bin"] = s_input(s_params(limit=raw(8)), S_POLY)
    out["s-tolerance.bin"] = s_input(s_params(), S_POLY, opts=8, tol=1)
    out["s-point-cap.bin"] = s_input(s_params(join=1), S_CURVE, opts=2)
    out["s-oom.bin"] = s_input(s_params(), S_POLY, opts=4, fail=0)
    out["s-oom-later.bin"] = s_input(s_params(join=1, cap=1), S_CURVE, opts=4, fail=5)
    out["s-inside.bin"] = s_input(s_params(align=1, rule=1), S_SQUARE)
    out["s-outside.bin"] = s_input(s_params(align=2), S_SQUARE + S_POLY)
    out["s-curve-round.bin"] = s_input(s_params(join=1, cap=1), S_CURVE)
    out["s-arcs.bin"] = s_input(s_params(join=4, limit=raw(160)), S_CURVE)
    out["s-miter-clip.bin"] = s_input(s_params(join=3, limit=raw(20), cap=2), S_POLY)
    out["s-device.bin"] = s_input(s_params(space=1), S_CURVE, map_bytes=bytes([1, 1]) + raw(48) + raw(16))
    out["s-degenerate.bin"] = s_input(s_params(cap=1), [bytes([0]) + s_pt(32, 32), bytes([6])])
    out["s-degenerate-square.bin"] = s_input(s_params(cap=2), [bytes([0]) + s_pt(32, 32), bytes([1]) + s_pt(32, 32)])
    out["s-cusp.bin"] = s_input(s_params(join=1, cap=1, width=raw(48)),
                                [bytes([0]) + s_pt(0, 0), bytes([3]) + s_pt(128, 64) + s_pt(64, 64) + s_pt(64, 0)])
    out["s-rational.bin"] = s_input(s_params(cap=1),
                                    [bytes([0]) + s_pt(0, 0), bytes([4]) + s_pt(80, 80) + s_pt(160, 0) + bytes([0, 1, 0])])
    # the curve outline (stroke_curves.h): two bytes after the path, the tolerance (an index into the C file's table: 0 is 2^-10,
    # 1 the finest, 5 below the finest, 6 NaN, 7 infinity) and the form (0 rational, 1 cubic, 2 out of range)
    out["s-curves-cubic.bin"] = s_input(s_params(join=1, cap=1), S_CURVE, tail=bytes([0, 1]))
    out["s-curves-rational.bin"] = s_input(s_params(join=1, cap=1), S_CURVE, tail=bytes([0, 0]))
    out["s-curves-aligned.bin"] = s_input(s_params(align=1, rule=1), S_SQUARE, tail=bytes([1, 0]))
    out["s-curves-dashed.bin"] = s_input(s_params(cap=1, dash=bytes([6, 0, 5, 3, 0, 0])), S_CURVE, tail=bytes([2, 1]))
    out["s-curves-bad-tolerance.bin"] = s_input(s_params(), S_POLY, tail=bytes([5, 0]))
    out["s-curves-nan-tolerance.bin"] = s_input(s_params(), S_POLY, tail=bytes([6, 0]))
    out["s-curves-infinite-tolerance.bin"] = s_input(s_params(), S_POLY, tail=bytes([7, 0]))
    out["s-curves-bad-form.bin"] = s_input(s_params(), S_POLY, tail=bytes([0, 2]))
    out["s-curves-limit.bin"] = s_input(s_params(join=1, cap=1), S_CURVE, opts=2, tail=bytes([0, 0]))
    out["s-curves-oom.bin"] = s_input(s_params(join=1, cap=1), S_CURVE, opts=4, fail=40, tail=bytes([0, 1]))
    out["s-curves-cusp.bin"] = s_input(s_params(join=1, cap=1, width=raw(48)),
                                       [bytes([0]) + s_pt(0, 0), bytes([3]) + s_pt(128, 64) + s_pt(64, 64) + s_pt(64, 0)],
                                       tail=bytes([3, 0]))
    out["s-curves-device.bin"] = s_input(s_params(space=1, join=1), S_CURVE, map_bytes=bytes([1, 1]) + raw(48) + raw(16),
                                         tail=bytes([2, 0]))
    # found by the fuzzer: a join between a segment of a few ulps' extent and its neighbour, whose direction needs a
    # norm that neither underflows nor overflows
    out["s-regress-norm.bin"] = base64.b64decode(
        "AFEEGAAyAAAYAJACAQABAAABsj9aAAAAAAAALAYBBgQYAE4Y/60Y/9MYALwHAwcHARgAhwQYAHYY/8cY/9cY/20HAooYAKQBAQUHABgAMBgAZA==")
    # found by the fuzzer: the clip part of a closed subpath that goes out and back is a ring of two vertices
    out["s-regress-clip-ring.bin"] = base64.b64decode(
        "AAcEGAAxAAIYALwBAQEAAAEBGAAkGAAnAAAFABgAMhj/+QEY/18YaGhoaGhoaGho/zkFGAAzGP/KGP/4GP+kGACMGP+3BgUABwAY//8vfRif/68DGP93GP95GP/YGP9tAwAAAG0YAH0YAH0YAFY=")
    # found by the fuzzer: a miter-clip join whose clipped corners lie a little past L h from the vertex
    out["s-regress-miter-clip-box.bin"] = base64.b64decode(
        "AAgEGDMQAQMYABoAAQAAAAAAAAUAGP9aGP8BAAAANhgAlxgAIhj/rwIYAEIY/+gY/8UAQAAAAP+PGP+yABgArRj/SA==")
    # found by the fuzzer: a coordinate that is not a number reached the snap table of the conforming seams, whose key is
    # the coordinate cast to an integer; the rounding that refuses it now comes first
    out["s-regress-snap-nan.bin"] = base64.b64decode(
        "AAAEGAAgQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQECKioqKQEBAQEBAQAABGABAAgAAAABL4AAGAQQDBQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIABBgAIAABGABAAQAAAAAAAAEDABcAABi4uLgAAAIYAFAYAKAYAAABGAAFABgAABgAAAEYAIAAQBgAAQAAAA==")
    return out


def stroke_biased(seed):
    rng = Rng(5000 + seed)
    opts = rng.below(4) if rng.below(3) == 0 else 0
    if rng.below(8) == 0:
        opts |= 4
    params = s_params(width=raw(8 + rng.below(80)), space=rng.below(2), join=rng.below(5), limit=raw(16 + rng.below(200)),
                      cap=rng.below(3), degenerate=rng.below(2), align=rng.below(3) if rng.below(4) == 0 else 0,
                      rule=rng.below(2))
    ops = [bytes([0]) + s_pt(rng.below(400) - 200, rng.below(400) - 200)]
    for _ in range(1 + rng.below(8)):
        k = rng.below(7)
        npts = {0: 1, 1: 1, 2: 2, 3: 3, 4: 2, 5: 3, 6: 0}[k]
        op = bytes([k]) + b"".join(s_pt(rng.below(400) - 200, rng.below(400) - 200) for _ in range(npts))
        if k in (4, 5):
            op += bytes(rng.below(8) for _ in range(3 if k == 4 else 4))
        ops.append(op)
    mp = bytes([0])
    if rng.below(3) == 0:
        mp = bytes([1, 1]) + raw(8 + rng.below(40)) + raw(8 + rng.below(40))
    return s_input(params, ops, opts=opts, fail=rng.below(12), map_bytes=mp)


STROKE_BIASED_FILES = 60
STROKE_PLAIN_FILES = 20


def stroke_files():
    out = {}
    for i in range(STROKE_BIASED_FILES):
        out["biased-%03d.bin" % i] = stroke_biased(i + 1)
    for i in range(STROKE_PLAIN_FILES):
        rng = Rng(7000 + i)
        out["plain-%03d.bin" % i] = bytes(rng.byte() for _ in range(3 + 20 + rng.below(160)))
    out.update(stroke_scripted())
    return out


# --- the clip corpus (tests/fuzz/fuzz_clip.c) -------------------------------------------------------------------------
#
# A file is: an option byte (bit 0 the whole int32 range, bit 1 a small max_points, bit 2 a failing allocation), a byte for the
# limit and the allocation to fail, a quadrilateral kind (0 a rectangle, 1 reversed, 2 a repeated point, 3 a bowtie, 4 a
# sheared rectangle, 5 four raw points) with its centre, half sizes and shear (or the four points), the number of contours and
# for each the number of points, each point's kind byte (0 a corner, 1 a point on a side, else a coordinate pair) and a closed
# flag.  A coordinate is a byte below 12 (a table of degenerate values; below 21 in the wide range) or 200 and two bytes (a
# 16-bit value, over four in the small range, times 32768 in the wide one).


def cc(v):
    """A small-range coordinate of value v (|v| <= 8191)."""
    r = (v * 4) & 0xFFFF
    return bytes([200, r >> 8, r & 255])


def cw(v16):
    """A wide-range coordinate of v16 * 32768."""
    r = v16 & 0xFFFF
    return bytes([200, r >> 8, r & 255])


def clip_quad(kind, cx, cy, hx, hy, shear=0, wide=False):
    c = cw if wide else cc
    scale = 1 if not wide else 1
    return (bytes([kind]) + c(cx) + c(cy) + c(2 * (hx - 1)) + c(2 * (hy - 1)) + bytes([shear & 255]))


def clip_raw_quad(points, wide=False):
    c = cw if wide else cc
    return bytes([5]) + b"".join(c(x) + c(y) for x, y in points)


def clip_poly(points, closed=1, wide=False):
    """points: ('p', x, y), ('c', corner) or ('s', a, b, f)."""
    c = cw if wide else cc
    out = bytes([len(points)])
    for pt in points:
        if pt[0] == "p":
            out += bytes([2]) + c(pt[1]) + c(pt[2])
        elif pt[0] == "c":
            out += bytes([0, pt[1]])
        else:
            out += bytes([1, pt[1], pt[2], pt[3]])
    return out + bytes([closed])


def clip_input(quad, polys, opts=0, fail=0):
    return bytes([opts, fail]) + quad + bytes([len(polys) - 1]) + b"".join(polys)


def clip_scripted():
    out = {}
    q = clip_quad(0, 20, 20, 30, 30)
    inside = clip_poly([("p", 10, 10), ("p", 40, 12), ("p", 20, 40)])
    cover = clip_poly([("p", -60, -60), ("p", 90, -60), ("p", 90, 90), ("p", -60, 90)])
    hole = clip_poly([("p", 5, 5), ("p", 5, 35), ("p", 35, 35), ("p", 35, 5)])
    cutting = clip_poly([("p", 0, 0), ("p", 120, 10), ("p", 30, 120)])
    far = clip_poly([("p", 400, 400), ("p", 500, 400), ("p", 450, 500)])
    out["invalid-collinear.bin"] = clip_input(clip_raw_quad([(0, 0), (10, 10), (20, 20), (30, 30)]), [inside])
    out["invalid-bowtie.bin"] = clip_input(clip_quad(3, 20, 20, 30, 30), [inside])
    out["invalid-one-point.bin"] = clip_input(clip_raw_quad([(5, 5)] * 4), [inside])
    out["limit.bin"] = clip_input(q, [cutting], opts=2, fail=1)
    out["oom-first.bin"] = clip_input(q, [cutting], opts=4, fail=0)
    out["oom-later.bin"] = clip_input(q, [cutting, hole], opts=4, fail=3)
    out["empty.bin"] = clip_input(q, [far])
    out["inside-open.bin"] = clip_input(q, [clip_poly([("p", 10, 10), ("p", 40, 12), ("p", 20, 40)], closed=0)])
    out["inside-closed.bin"] = clip_input(q, [inside, hole])
    out["cut.bin"] = clip_input(q, [cutting, inside])
    out["cover-and-hole.bin"] = clip_input(q, [cover, hole])
    out["clockwise.bin"] = clip_input(clip_quad(1, 20, 20, 30, 30), [cover, inside])
    out["repeated-point.bin"] = clip_input(clip_quad(2, 20, 20, 30, 30), [cutting, inside])
    out["sheared.bin"] = clip_input(clip_quad(4, 20, 20, 30, 30, shear=24), [cutting, cover])
    out["on-the-line.bin"] = clip_input(q, [clip_poly([("c", 0), ("c", 1), ("p", 60, 60), ("s", 1, 2, 4)])])
    out["corners-and-sides.bin"] = clip_input(q, [clip_poly([("c", 0), ("s", 0, 1, 3), ("c", 2), ("s", 3, 0, 5)]),
                                                  clip_poly([("s", 1, 2, 8), ("c", 3), ("p", 100, 100)], closed=0)])
    wq = clip_quad(0, 100, 50, 200, 150, wide=True)
    wp = clip_poly([("p", -300, -300), ("p", 300, -200), ("p", 100, 300)], wide=True)
    out["wide.bin"] = clip_input(wq, [wp, wp], opts=1)
    out["wide-extremes.bin"] = bytes([1, 0, 0]) + bytes([12, 13, 12, 14, 0, 0, 0, 0, 0]) + b"\x00\x01\x03\x02" + bytes([12, 13, 14, 15, 16]) + b"\x01"
    out["wide-invalid.bin"] = clip_input(clip_raw_quad([(0, 0), (5, 5), (10, 10), (20, 20)], wide=True), [wp], opts=1)
    out["wide-limit.bin"] = clip_input(wq, [wp], opts=3, fail=0)
    return out


def clip_biased(seed):
    rng = Rng(seed)
    opts = (1 if rng.below(5) == 0 else 0) | (2 if rng.below(8) == 0 else 0) | (4 if rng.below(8) == 0 else 0)
    wide = opts & 1
    c = cw if wide else cc
    span = 120 if wide else 400
    if rng.below(10) == 0:
        quad = clip_raw_quad([(rng.below(2 * span) - span, rng.below(2 * span) - span) for _ in range(4)], bool(wide))
    else:
        quad = clip_quad(rng.below(5), rng.below(span) - span // 2, rng.below(span) - span // 2, 4 + rng.below(span // 4),
                         4 + rng.below(span // 4), rng.below(40) - 20, bool(wide))
    polys = []
    for _ in range(1 + rng.below(3)):
        pts = []
        for _ in range(1 + rng.below(9)):
            k = rng.below(6)
            if k == 0:
                pts.append(("c", rng.below(4)))
            elif k == 1:
                pts.append(("s", rng.below(4), rng.below(4), rng.below(9)))
            else:
                pts.append(("p", rng.below(2 * span) - span, rng.below(2 * span) - span))
        polys.append(clip_poly(pts, rng.below(4) != 0, bool(wide)))
    return clip_input(quad, polys, opts=opts, fail=rng.below(24))


CLIP_BIASED_FILES = 100
CLIP_PLAIN_FILES = 30


def clip_files():
    out = {}
    for i in range(CLIP_BIASED_FILES):
        out["biased-%03d.bin" % i] = clip_biased(i + 1)
    for i in range(CLIP_PLAIN_FILES):
        rng = Rng(9000 + i)
        out["plain-%03d.bin" % i] = bytes(rng.byte() for _ in range(3 + 20 + rng.below(200)))
    out.update(clip_scripted())
    return out


# --- the hatch corpus (tests/fuzz/fuzz_hatch.c) -----------------------------------------------------------------------
#
# A file is: an option byte (bit 0 a small max_points, bit 1 a failing allocation, bit 2 a small line limit, bit 3 device
# width, bits 4 and 5 the cap, bit 6 degenerate omit, bit 7 even-odd), a byte for the limit and the allocation to fail, the hatch
# map and the chain (a kind byte below 7, then four values and two for kind 4), the angle unit, the nesting, the family count (mod
# 4) with each family's angle, origin x and y, offsets along and across and its dash count with the dashes, the style's width,
# and the contours (a count, then for each a kind byte, x, y, w and h).  A value is a byte below 24 (a table of degenerate and
# round numbers, in fuzz_hatch.c), 24 for NaN, 25 for infinity, or 30 and a 16-bit value over 8.

HV = {0: 0.0, 0.25: 1, 0.5: 2, 1.0: 3, 1.5: 4, 2.0: 5, 3.0: 6, 4.0: 7, 7.5: 8, 10.0: 9, -1.0: 10, -2.0: 11, 90.0: 16, 45.0: 17,
      30.0: 18, 5.0: 22, 12.0: 23}


def hv(x):
    """A value: a table byte if x is in the table, else the three-byte form."""
    if x == 0:
        return bytes([0])
    if x in HV:
        return bytes([HV[x]])
    r = int(round(x * 8)) & 0xFFFF
    return bytes([30, r >> 8, r & 255])


def hmap(kind, *vals):
    return bytes([kind]) + b"".join(hv(v) for v in vals)


H_ID = hmap(0)
H_SCALE = hmap(1)
H_ROT = hmap(2)
H_SING = hmap(3)
H_3D = hmap(5)
H_FLIP = hmap(6)


def hfam(angle, ox, oy, dx, dy, dashes=()):
    return hv(angle) + hv(ox) + hv(oy) + hv(dx) + hv(dy) + bytes([len(dashes)]) + b"".join(hv(d) for d in dashes)


def hcontour(kind, x, y, w, h):
    return bytes([kind]) + hv(x) + hv(y) + hv(w) + hv(h)


def hatch_input(fams, contours, opts=0, fail=0, hm=H_ID, cm=H_ID, unit=1, nesting=0, width=0.5):
    return (bytes([opts, fail]) + hm + cm + bytes([unit, nesting, len(fams)]) + b"".join(fams) + hv(width) +
            bytes([len(contours) - 1]) + b"".join(contours))


def hatch_scripted():
    out = {}
    sq = hcontour(0, 0, 0, 10, 10)
    hole = hcontour(0, 3, 3, 4, 4)
    tri = hcontour(2, 0, 0, 12, 9)
    solid = hfam(0, 0, 0, 0, 2.0)
    diag = hfam(45, 0, 0, 0, 1.5)
    dashed = hfam(0, 0, 0, 0, 2.0, (3.0, -1.0))
    dotted = hfam(0, 0, 0, 0, 1.5, (0, -2.0))
    out["solid.bin"] = hatch_input([solid], [sq])
    out["diag-round.bin"] = hatch_input([diag], [sq, hole], opts=1 << 4 | 0)
    out["dashed-square.bin"] = hatch_input([dashed, diag], [sq], opts=2 << 4)
    out["dots-square-cap.bin"] = hatch_input([dotted], [sq], opts=2 << 4)
    out["dots-round.bin"] = hatch_input([dotted], [tri], opts=1 << 4)
    out["dots-omit.bin"] = hatch_input([dotted], [sq], opts=(1 << 4) | 64)
    out["outer.bin"] = hatch_input([solid], [sq, hole], nesting=1)
    out["ignore.bin"] = hatch_input([solid, dashed], [sq, hole], nesting=2)
    out["device.bin"] = hatch_input([dashed], [sq], opts=8 | (1 << 4), width=1.0)
    out["even-odd.bin"] = hatch_input([solid], [sq, hole], opts=128)
    out["scaled.bin"] = hatch_input([dashed], [tri], hm=H_SCALE, cm=H_ROT)
    out["flipped.bin"] = hatch_input([diag], [sq], hm=H_FLIP, cm=H_SCALE)
    out["zero-width.bin"] = hatch_input([solid], [sq], width=0)
    out["singular-hatch.bin"] = hatch_input([solid], [sq], hm=H_SING)
    out["singular-chain.bin"] = hatch_input([solid], [sq], cm=H_SING)
    out["3d.bin"] = hatch_input([solid], [sq], hm=H_3D)
    out["no-family.bin"] = hatch_input([], [sq])
    out["bad-nesting.bin"] = hatch_input([solid], [sq], nesting=3)
    out["bad-unit.bin"] = hatch_input([solid], [sq], unit=4)
    out["zero-spacing.bin"] = hatch_input([hfam(0, 0, 0, 0, 0), hfam(90, 0, 0, 0, 2.0)], [sq])
    out["no-lines.bin"] = hatch_input([hfam(0, 100, 100, 0, 2.0)], [hcontour(0, 0, 0, 0, 0)])
    out["line-limit.bin"] = hatch_input([solid], [sq], opts=4, fail=1)
    out["point-limit.bin"] = hatch_input([dashed], [sq], opts=1, fail=0)
    out["oom-first.bin"] = hatch_input([solid], [sq], opts=2, fail=0)
    out["oom-later.bin"] = hatch_input([dashed], [sq, hole], opts=2, fail=7)
    out["bad-values.bin"] = hatch_input([bytes([24, 0, 0, 0, 3, 0])], [sq])
    out["infinite-dash.bin"] = hatch_input([bytes([0, 0, 0, 0, 3, 1, 25])], [sq])
    out["spaces-only.bin"] = hatch_input([hfam(0, 0, 0, 0, 2.0, (-1.0, -2.0))], [sq])
    return out


def hatch_random(seed):
    rng = Rng(seed)
    vals = [0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 7.5, 10.0, -1.0, -2.0, 5.0, 12.0]
    angles = [0, 90, 45, 30, 17.0, 135.0, 180.0, 270.0]
    fams = []
    for _ in range(1 + rng.below(3)):
        nd = rng.below(6)
        dashes = tuple(vals[rng.below(len(vals))] * (1 if rng.below(3) else 0) for _ in range(nd))
        fams.append(hfam(angles[rng.below(len(angles))] if angles[rng.below(len(angles))] in HV else 45, vals[rng.below(len(vals))] - 1.0,
                         vals[rng.below(len(vals))] - 1.0, [0, 0.5, 1.5, -1.0][rng.below(4)], vals[rng.below(len(vals))], dashes))
    contours = [hcontour(rng.below(3), vals[rng.below(len(vals))] - 1.0, vals[rng.below(len(vals))] - 1.0, vals[rng.below(len(vals))] + 4.0,
                         vals[rng.below(len(vals))] + 4.0) for _ in range(1 + rng.below(3))]
    maps = [H_ID, H_SCALE, H_ROT, H_FLIP, H_ID]
    opts = rng.below(256) & ~2 if rng.below(4) else rng.below(256)
    return hatch_input(fams, contours, opts=opts, fail=rng.below(40), hm=maps[rng.below(5)], cm=maps[rng.below(5)],
                       unit=1 if rng.below(5) else rng.below(5), nesting=rng.below(3), width=[0.25, 0.5, 1.0, 1.5][rng.below(4)])


HATCH_RANDOM_FILES = 90
HATCH_PLAIN_FILES = 30


def hatch_files():
    out = hatch_scripted()
    # Found by pass B's libFuzzer run: a hatch with a 3D map and a region past the point limit, which answered with the region's
    # limit where the hatch check's unsupported comes first.
    out["regression-refused-hatch-and-a-region-past-the-limit.bin"] = base64.b64decode("AQAFAAEAAQAAAAAFAAICAAUACQk=")
    for i in range(HATCH_RANDOM_FILES):
        out["random-%03d.bin" % i] = hatch_random(i + 1)
    for i in range(HATCH_PLAIN_FILES):
        rng = Rng(19000 + i)
        out["plain-%03d.bin" % i] = bytes(rng.byte() for _ in range(40 + rng.below(160)))
    return out


def check_or_write(root, want, check):
    if check:
        have = sorted(os.listdir(root)) if os.path.isdir(root) else []
        if have != sorted(want):
            sys.exit("gen-fuzz-corpus: the committed corpus in %s is not the generator's set of files" % root)
        for name, data in want.items():
            with open(os.path.join(root, name), "rb") as f:
                if f.read() != data:
                    sys.exit("gen-fuzz-corpus: %s differs from the generator's" % name)
        print("gen-fuzz-corpus: the committed corpus (%d files) in %s is current" % (len(want), os.path.basename(root)))
        return
    os.makedirs(root, exist_ok=True)
    for name, data in want.items():
        with open(os.path.join(root, name), "wb") as f:
            f.write(data)
    print("gen-fuzz-corpus: wrote %d files" % len(want))


def main():
    base = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tests", "fuzz", "corpus")
    check = "--check" in sys.argv
    check_or_write(os.path.join(base, "path"), files(), check)
    check_or_write(os.path.join(base, "stroke"), stroke_files(), check)
    check_or_write(os.path.join(base, "clip"), clip_files(), check)
    check_or_write(os.path.join(base, "hatch"), hatch_files(), check)


if __name__ == "__main__":
    main()
