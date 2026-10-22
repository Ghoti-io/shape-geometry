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
"""Write the seed corpora of the path and stroke fuzzers, deterministically (story 16; CAP-15, AD-19).

    python3 tools/gen-fuzz-corpus.py            # write tests/fuzz/corpus/path/ and tests/fuzz/corpus/stroke/
    python3 tools/gen-fuzz-corpus.py --check    # fail if the committed files differ

The stroke corpus (tests/fuzz/fuzz_stroke.c) is described at stroke_files() below.

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
       "CUBIC_FIXED", "DEVICE", "MEASURE", "ORIENT", "BOUNDS"]
OP = {name: i for i, name in enumerate(OPS)}

BIASED = {"MOVE": 3, "LINE": 5, "QUAD": 5, "CUBIC": 6, "RQUAD": 4, "RCUBIC": 3, "CLOSE": 3, "MOVE_FIXED": 1,
          "LINE_FIXED": 2, "QUAD_FIXED": 2, "CUBIC_FIXED": 2, "DEVICE": 6, "MEASURE": 4, "ORIENT": 3, "BOUNDS": 1}

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


def scripted():
    """Inputs that reach, by construction, the outcomes that random bytes reach only by luck."""
    out = {}
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


def s_input(params, ops, opts=0, fail=0, tol=4, map_bytes=b"\x00", map3d=0):
    return bytes([opts, fail, tol]) + params + map_bytes + bytes([map3d]) + s_path(ops)


def stroke_scripted():
    out = {}
    out["s-drawn.bin"] = s_input(s_params(), S_POLY)
    out["s-empty.bin"] = s_input(s_params(width=bytes([S_ZERO])), S_POLY)
    # a map step of kind 1 (a scale) by 0 and 1: singular
    out["s-singular.bin"] = s_input(s_params(), S_POLY, map_bytes=bytes([1, 1, S_ZERO, S_ONE]))
    out["s-invalid.bin"] = s_input(s_params(space=2), S_POLY)
    out["s-map-3d.bin"] = s_input(s_params(), S_POLY, map3d=9)
    out["s-dashes.bin"] = s_input(s_params(dash=bytes([1]) + raw(32) + raw(16)), S_POLY)
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


if __name__ == "__main__":
    main()
