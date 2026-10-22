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
"""Render SVG strokes with resvg and record what it covers, for tests/stroke_resvg.h.

The behavior check of the stroke rules' reading (stroke.h): resvg is a second,
external reading of SVG 2's stroke shape, so a misreading that the expansion and
the oracle share (both written from the same text) shows as a disagreement.  It
supplies behavior only: nothing of it is read.

Each case is an SVG path with an SVG-expressible style (miter, round and bevel
joins, the three caps, a miter limit, curves, closed and overlapping subpaths, a
zero-length subpath with round caps), rendered aliased (`shape-rendering:
crispEdges`) at 16 pixels per user unit.  A sample point is the centre of every
other pixel on both axes (every 1/8 unit); its verdict is the pixel's alpha at
least half.  tests/unit/test_stroke_resvg.cpp compares the expansion's outline
with these verdicts at every point the geometric oracle decides with a margin
for resvg's own flattening and rasterization (a pixel at 16x and its flattening
tolerance).

Run it in the pinned image (tools/resvg-stroke/fetch.sh builds it):

    podman run --rm -v "$PWD:/w:Z" -w /w localhost/ghoti-shape-oracle-resvg:0.45.1 \\
        python3 tools/resvg-stroke/gen.py > tests/stroke_resvg.h

`--check` regenerates into memory and fails unless the committed header is the
same, so a changed renderer is seen.
"""
import os
import struct
import subprocess
import sys
import zlib

SCALE = 16
STEP = 2  # pixels between samples

JOIN = {"miter": 0, "round": 1, "bevel": 2}
CAP = {"butt": 0, "round": 1, "square": 2}

CASES = [
    ("polyline_miter", [("M", [2, 2]), ("L", [14, 3]), ("L", [5, 12]), ("L", [18, 16])], 2.5, "miter", 4, "butt"),
    ("polyline_miter_limited", [("M", [2, 2]), ("L", [14, 3]), ("L", [5, 12]), ("L", [18, 16])], 2.5, "miter", 1.5,
     "square"),
    ("polyline_round", [("M", [2, 2]), ("L", [14, 3]), ("L", [5, 12]), ("L", [18, 16])], 2.5, "round", 4, "round"),
    ("polyline_bevel", [("M", [2, 2]), ("L", [14, 3]), ("L", [5, 12]), ("L", [18, 16])], 2.5, "bevel", 4, "butt"),
    ("acute_miter", [("M", [2, 4]), ("L", [18, 6]), ("L", [3, 9])], 2, "miter", 10, "butt"),
    ("quadratic_round", [("M", [2, 14]), ("Q", [10, -4, 18, 14])], 3, "round", 4, "round"),
    ("cubic_square", [("M", [2, 10]), ("C", [6, 0, 14, 20, 18, 10])], 2, "miter", 4, "square"),
    ("cubic_loop", [("M", [3, 15]), ("C", [22, 2, -2, 2, 17, 15])], 1.5, "round", 4, "butt"),
    ("closed_square_miter", [("M", [4, 4]), ("L", [16, 4]), ("L", [16, 16]), ("L", [4, 16]), ("Z", [])], 3, "miter",
     4, "butt"),
    ("closed_triangle_bevel", [("M", [3, 16]), ("L", [10, 3]), ("L", [17, 16]), ("Z", [])], 2.5, "bevel", 4, "butt"),
    ("figure_eight", [("M", [10, 10]), ("C", [20, 20, 20, 0, 10, 10]), ("C", [0, 20, 0, 0, 10, 10]), ("Z", [])], 2,
     "round", 4, "round"),
    ("overlapping_subpaths", [("M", [2, 2]), ("L", [18, 18]), ("M", [2, 18]), ("Q", [10, -2, 18, 10])], 3, "miter", 4,
     "round"),
    ("dot_round", [("M", [10, 10]), ("Z", [])], 6, "miter", 4, "round"),
    ("cusp_round", [("M", [3, 15]), ("C", [19, 5, 11, 5, 11, 15])], 2, "round", 4, "round"),
]

BOX = (0.0, 0.0, 20.0, 20.0)


def svg_d(ops):
    parts = []
    for verb, a in ops:
        parts.append(verb + " " + " ".join(repr(float(v)) for v in a))
    return " ".join(parts)


def svg(case):
    name, ops, width, join, limit, cap = case
    x0, y0, w, h = BOX
    return ('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="%r %r %r %r" '
            'shape-rendering="crispEdges"><path d="%s" fill="none" stroke="#000000" stroke-width="%r" '
            'stroke-linejoin="%s" stroke-linecap="%s" stroke-miterlimit="%r"/></svg>\n' % (
                int(w * SCALE), int(h * SCALE), x0, y0, w, h, svg_d(ops), float(width), join, cap, float(limit)))


def read_png(path):
    """The alpha channel of an 8-bit RGBA, non-interlaced PNG, as rows of bytes."""
    data = open(path, "rb").read()
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        raise SystemExit("resvg-stroke: %s is not a PNG" % path)
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
                raise SystemExit("resvg-stroke: %s is not 8-bit RGBA non-interlaced" % path)
        elif kind == b"IDAT":
            idat += body
    raw = zlib.decompress(idat)
    bpp, stride = 4, width * 4
    rows = []
    prev = bytearray(stride)
    i = 0
    for _ in range(height):
        f = raw[i]
        line = bytearray(raw[i + 1:i + 1 + stride])
        i += 1 + stride
        for x in range(stride):
            a = line[x - bpp] if x >= bpp else 0
            b = prev[x]
            c = prev[x - bpp] if x >= bpp else 0
            if f == 1:
                line[x] = (line[x] + a) & 255
            elif f == 2:
                line[x] = (line[x] + b) & 255
            elif f == 3:
                line[x] = (line[x] + ((a + b) >> 1)) & 255
            elif f == 4:
                p = a + b - c
                pa, pb, pc = abs(p - a), abs(p - b), abs(p - c)
                pr = a if (pa <= pb and pa <= pc) else (b if pb <= pc else c)
                line[x] = (line[x] + pr) & 255
        rows.append(bytes(line[3::4]))
        prev = line
    return width, height, rows


def rle(bits):
    """Run lengths of a 0/1 string, starting with a run of 0s (possibly empty), comma-separated."""
    out = []
    cur = "0"
    n = 0
    for b in bits:
        if b == cur:
            n += 1
        else:
            out.append(n)
            cur = b
            n = 1
    out.append(n)
    return ",".join(str(v) for v in out)


def render(case, tmp):
    name = case[0]
    sp = os.path.join(tmp, name + ".svg")
    pp = os.path.join(tmp, name + ".png")
    with open(sp, "w") as f:
        f.write(svg(case))
    subprocess.run(["resvg", sp, pp], check=True)
    width, height, alpha = read_png(pp)
    cols, rows = width // STEP, height // STEP
    bits = []
    for j in range(rows):
        for i in range(cols):
            bits.append("1" if alpha[j * STEP][i * STEP] >= 128 else "0")
    return cols, rows, rle("".join(bits))


def hexf(v):
    return float(v).hex()


def generate():
    version = subprocess.run(["resvg", "--version"], capture_output=True, text=True, check=True).stdout.strip()
    out = []
    out.append("/* Generated by tools/resvg-stroke/gen.py from resvg %s; do not edit. */\n" % version)
    out.append("/* SPDX-License-Identifier: LGPL-3.0-only */\n")
    out.append("#ifndef GHOTI_IO_GSHGEOM_TESTS_STROKE_RESVG_H\n#define GHOTI_IO_GSHGEOM_TESTS_STROKE_RESVG_H\n\n")
    out.append("#include <cstddef>\n\n#include \"stroke_vectors.h\"\n\nnamespace stroke_resvg {\n\n")
    out.append("static const char * const resvg_version = \"%s\";\n" % version)
    out.append("static const int scale = %d; ///< pixels per user unit\nstatic const int step = %d; ///< pixels between samples\n\n" % (SCALE, STEP))
    out.append("struct Case {\n  const char * name;\n  const stroke_vectors::Op * ops;\n  size_t op_count;\n  double width;\n"
               "  int join;\n  double limit;\n  int cap;\n  double box[4];\n  int cols, rows;\n"
               "  const char * runs; ///< run lengths of the row-major verdicts, starting with outside\n};\n\n")
    tmp = os.path.join(os.environ.get("TMPDIR", "/tmp"), "resvg-stroke-%d" % os.getpid())
    os.makedirs(tmp, exist_ok=True)
    if True:
        rows_out = []
        for ci, case in enumerate(CASES):
            name, ops, width, join, limit, cap = case
            cols, rows, runs = render(case, tmp)
            out.append("static const stroke_vectors::Op resvg_ops_%d[] = {\n" % ci)
            for verb, a in ops:
                vals = ", ".join(hexf(v) for v in a) if a else "0"
                out.append("    {'%s', %d, {%s}},\n" % (verb, len(a), vals))
            out.append("};\n")
            rows_out.append('    {"%s", resvg_ops_%d, %d, %s, %d, %s, %d, {%s}, %d, %d,\n     "%s"},\n' % (
                name, ci, len(ops), hexf(width), JOIN[join], hexf(limit), CAP[cap],
                ", ".join(hexf(v) for v in BOX), cols, rows, runs))
    for f in os.listdir(tmp):
        os.remove(os.path.join(tmp, f))
    os.rmdir(tmp)
    out.append("\nstatic const Case cases[] = {\n")
    out.extend(rows_out)
    out.append("};\nstatic const size_t case_count = sizeof cases / sizeof cases[0];\n\n")
    out.append("} // namespace stroke_resvg\n\n#endif /* GHOTI_IO_GSHGEOM_TESTS_STROKE_RESVG_H */\n")
    return "".join(out)


def main():
    text = generate()
    if "--check" in sys.argv:
        here = os.path.dirname(os.path.abspath(__file__))
        committed = open(os.path.join(here, "..", "..", "tests", "stroke_resvg.h")).read()
        if committed != text:
            sys.exit("resvg-stroke: tests/stroke_resvg.h is not what this resvg renders")
        print("resvg-stroke: tests/stroke_resvg.h is what resvg renders")
        return
    sys.stdout.write(text)


if __name__ == "__main__":
    main()
