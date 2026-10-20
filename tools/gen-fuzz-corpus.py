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
"""Write the seed corpus of the path fuzzer, deterministically (story 16; CAP-15, AD-19).

    python3 tools/gen-fuzz-corpus.py            # write tests/fuzz/corpus/path/
    python3 tools/gen-fuzz-corpus.py --check    # fail if the committed files differ

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


def main():
    root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tests", "fuzz", "corpus", "path")
    want = files()
    if "--check" in sys.argv:
        have = sorted(os.listdir(root)) if os.path.isdir(root) else []
        if have != sorted(want):
            sys.exit("gen-fuzz-corpus: the committed corpus is not the generator's set of files")
        for name, data in want.items():
            with open(os.path.join(root, name), "rb") as f:
                if f.read() != data:
                    sys.exit("gen-fuzz-corpus: %s differs from the generator's" % name)
        print("gen-fuzz-corpus: the committed corpus (%d files) is current" % len(want))
        return
    os.makedirs(root, exist_ok=True)
    for name, data in want.items():
        with open(os.path.join(root, name), "wb") as f:
            f.write(data)
    print("gen-fuzz-corpus: wrote %d files" % len(want))


if __name__ == "__main__":
    main()
