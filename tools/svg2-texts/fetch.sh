#!/bin/sh
# SPDX-License-Identifier: LGPL-3.0-only
#
# Copyright (C) 2026 Corey Pennycuff
#
# Fetch the pinned texts the stroke rules (stroke.h, SG5 to SG13) were checked
# against, and verify each by its SHA-256.  Nothing is built from them: they are
# the specification, kept so that a reader can check a rule against the exact
# edition it cites.
#
#   - SVG 2, W3C Candidate Recommendation 2018-10-04: the painting chapter
#     (stroke shape, line caps, line joins including miter-clip and arcs, the
#     miter limit, the arcs circles and their adjustment) and the paths chapter
#     (the direction at the start and end of a path segment, zero-length
#     segments).  The arcs and miter-clip joins are defined only here: the
#     editor's draft below has withdrawn both values.
#   - SVG 2 editor's draft, w3c/svgwg at commit f846f12e3f50: the same two
#     chapters in their source form.
#   - The SVG Strokes module draft at the same commit (stroke-alignment).
#   - The WHATWG HTML standard's source at commit c6c3055c1287: the canvas
#     line styles and the "trace a path" algorithm.
#
# Usage: tools/svg2-texts/fetch.sh [DIR]   (default build/svg2-texts)
set -eu
dir="${1:-build/svg2-texts}"
mkdir -p "$dir"

svgwg=https://raw.githubusercontent.com/w3c/svgwg/f846f12e3f504aec18f6ba3dacfb8db7ab7e0e36
whatwg=https://raw.githubusercontent.com/whatwg/html/c6c3055c1287bf50070260d53555d8d6ab530771

fetch() {
  name="$1"
  url="$2"
  want="$3"
  if [ ! -f "$dir/$name" ]; then
    curl -sSfL -o "$dir/$name.part" "$url"
    mv "$dir/$name.part" "$dir/$name"
  fi
  got=$(sha256sum "$dir/$name" | awk '{print $1}')
  if [ "$got" != "$want" ]; then
    printf 'svg2-texts: %s has SHA-256 %s, pinned %s\n' "$name" "$got" "$want" >&2
    exit 1
  fi
  printf 'svg2-texts: %s ok\n' "$name"
}

fetch svg2-cr-20181004-painting.html https://www.w3.org/TR/2018/CR-SVG2-20181004/painting.html \
  7d6def7e874c400920e09c534ed52307c66ec7f4275c15a8c7f08c54bf73dca7
fetch svg2-cr-20181004-paths.html https://www.w3.org/TR/2018/CR-SVG2-20181004/paths.html \
  6213e9b47b068443841fc840dde69afdd188a25aa6ea066e284b612b28f44a8d
fetch svg2-ed-painting.bs "$svgwg/master/painting.bs" \
  ce4d9c621dc95aa12e0f0a76dc41fb941c0a1bebb87a93297715af8d803caffc
fetch svg2-ed-paths.bs "$svgwg/master/paths.bs" \
  e3efdc40baab30f99bc3d2ca7dc968395aded3ed5b6f91dd987366e8c25e7838
fetch svg-strokes-ed.html "$svgwg/specs/strokes/master/Overview.html" \
  fe5458059120c80609f1a9d8108c448058b3c51f33887744c7200c3ce413889f
fetch whatwg-html-source "$whatwg/source" \
  915736f7cceded6c0b2f3a886aead8a1b058d7d01f512602731789f00abda9fa
