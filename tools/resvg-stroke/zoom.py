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
"""The zoom gate of the curve outline (stroke_curves.h, SG43).

A stroke written as a fill is exact only to the tolerance of the file it is written in, and a viewer zooms.  For each
case of tools/stroke-curves-svg.c this renders, with the pinned resvg and aliased, three documents at 1x, 16x and 64x:
the stroke as a stroke (the renderer's own reading), the curve outline in form CUBIC written as a fill, and the control,
the polygon the polygon expansion makes for the 1x pixel, written as a fill.  At each zoom `z` the fill must agree with
the stroke:

  * from 16x up, the **mean displacement** of the boundary, `XOR / B` (the pixels on which the two renderings differ over
    the pixels on the stroke rendering's boundary), is at most `t z + M` with `t` the file's tolerance (2^-10) and
    `M` = 0.10 pixel, the renderer's own stroker's disagreement with the exact region plus the quantization of an aliased
    edge (at 1x a dot of a dozen boundary pixels flips by a quarter of them, which is why the mean gate starts at 16x); and
  * no pixel on which they differ is farther than `t z + W` (W = 1.5 pixels) from the stroke rendering's boundary.

And the control must be seen to fail the same gate on every curved case at 16x and 64x: a gate that cannot fail on the
facets this story exists to remove measures nothing.

Run it in the pinned image (tools/resvg-stroke/fetch.sh builds it) over the directory tools/stroke-curves-svg wrote:

    podman run --rm -v "$PWD:/w:Z" -w /w localhost/ghoti-shape-oracle-resvg:0.45.1 \\
        python3 tools/resvg-stroke/zoom.py build/linux/release/zoom

Exit status 0 when the fills pass and the controls fail, 1 otherwise.
"""
import os
import struct
import subprocess
import sys
import zlib

T = 2.0 ** -10
MEAN_M = 0.10
WORST_W = 1.5
ZOOMS = (1, 16, 64)


def alpha_of(path):
    """The alpha channel of an 8-bit RGBA non-interlaced PNG, as a list of bytearrays (rows)."""
    data = open(path, "rb").read()
    pos = 8
    idat = b""
    width = height = None
    while pos < len(data):
        length = struct.unpack(">I", data[pos:pos + 4])[0]
        kind = data[pos + 4:pos + 8]
        body = data[pos + 8:pos + 8 + length]
        pos += 12 + length
        if kind == b"IHDR":
            width, height, depth, ctype, _, _, interlace = struct.unpack(">IIBBBBB", body)
            if depth != 8 or ctype != 6 or interlace != 0:
                raise SystemExit("zoom: %s is not 8-bit RGBA non-interlaced" % path)
        elif kind == b"IDAT":
            idat += body
    raw = zlib.decompress(idat)
    stride = width * 4
    rows = []
    prev = bytearray(stride)
    i = 0
    for _ in range(height):
        f = raw[i]
        line = bytearray(raw[i + 1:i + 1 + stride])
        i += 1 + stride
        if f == 1:
            for x in range(4, stride):
                line[x] = (line[x] + line[x - 4]) & 255
        elif f == 2:
            line = bytearray(((a + b) & 255) for a, b in zip(line, prev))
        elif f == 3:
            for x in range(stride):
                a = line[x - 4] if x >= 4 else 0
                line[x] = (line[x] + ((a + prev[x]) >> 1)) & 255
        elif f == 4:
            for x in range(stride):
                a = line[x - 4] if x >= 4 else 0
                b = prev[x]
                c = prev[x - 4] if x >= 4 else 0
                p = a + b - c
                pa, pb, pc = abs(p - a), abs(p - b), abs(p - c)
                pr = a if (pa <= pb and pa <= pc) else (b if pb <= pc else c)
                line[x] = (line[x] + pr) & 255
        rows.append(bytearray(1 if v >= 128 else 0 for v in line[3::4]))
        prev = line
    return width, height, rows


def render(svg, z, out):
    subprocess.run(["resvg", "-z", str(z), svg, out], check=True)
    return alpha_of(out)


def compare(ref, got):
    """(mean displacement, worst displacement) of `got` against `ref`, in pixels."""
    h = len(ref)
    w = len(ref[0])
    boundary = []
    for y in range(h):
        row = ref[y]
        for x in range(w):
            v = row[x]
            if (x > 0 and row[x - 1] != v) or (x + 1 < w and row[x + 1] != v) or (y > 0 and ref[y - 1][x] != v) or \
                    (y + 1 < h and ref[y + 1][x] != v):
                boundary.append((x, y))
    diff = []
    for y in range(h):
        a, b = ref[y], got[y]
        if a != b:
            for x in range(w):
                if a[x] != b[x]:
                    diff.append((x, y))
    nb = max(len(boundary), 1)
    mean = len(diff) / nb
    # the worst: each differing pixel's distance to the nearest boundary pixel, by buckets of 16 pixels
    cell = 16
    buckets = {}
    for (x, y) in boundary:
        buckets.setdefault((x // cell, y // cell), []).append((x, y))
    worst = 0.0
    for (x, y) in diff:
        best = 1e9
        r = 1
        while True:
            for gy in range(y // cell - r, y // cell + r + 1):
                for gx in range(x // cell - r, x // cell + r + 1):
                    for (bx, by) in buckets.get((gx, gy), ()):
                        d = ((bx - x) ** 2 + (by - y) ** 2) ** 0.5
                        if d < best:
                            best = d
            if best <= (r - 1) * cell + 1 or r > 8 or (best < 1e9 and best <= r * cell):
                break
            r += 1
        worst = max(worst, best if best < 1e9 else 1e9)
    return mean, worst


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: zoom.py DIR")
    d = sys.argv[1]
    cases = []
    for line in open(os.path.join(d, "cases.txt")):
        name, curved = line.split()
        cases.append((name, int(curved)))
    tmp = os.path.join(os.environ.get("TMPDIR", "/tmp"), "zoom-%d" % os.getpid())
    os.makedirs(tmp, exist_ok=True)
    failed = False
    for name, curved in cases:
        for z in ZOOMS:
            ref = render(os.path.join(d, name + ".stroke.svg"), z, os.path.join(tmp, "r.png"))[2]
            fill = render(os.path.join(d, name + ".fill.svg"), z, os.path.join(tmp, "f.png"))[2]
            poly = render(os.path.join(d, name + ".poly.svg"), z, os.path.join(tmp, "p.png"))[2]
            fm, fw = compare(ref, fill)
            pm, pw = compare(ref, poly)
            limit_m, limit_w = T * z + MEAN_M, T * z + WORST_W
            use_mean = z >= 16
            fill_ok = (not use_mean or fm <= limit_m) and fw <= limit_w
            poly_fails = (use_mean and pm > limit_m) or pw > limit_w
            need_fail = curved and z >= 16
            verdict = []
            if not fill_ok:
                verdict.append("FILL FAILS")
                failed = True
            if need_fail and not poly_fails:
                verdict.append("CONTROL DID NOT FAIL")
                failed = True
            print("%-16s %2dx fill mean %.3f worst %.2f | control mean %.3f worst %.2f | limit %.3f %.2f %s" % (
                name, z, fm, fw, pm, pw, limit_m, limit_w, " ".join(verdict) or "ok"))
            sys.stdout.flush()
    for f in os.listdir(tmp):
        os.remove(os.path.join(tmp, f))
    os.rmdir(tmp)
    if failed:
        sys.exit(1)
    print("zoom: the fills agree with the strokes at 1x, 16x and 64x, and the 1x polygon fails at 16x and 64x on every curved case")


if __name__ == "__main__":
    main()
