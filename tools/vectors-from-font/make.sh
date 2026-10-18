#!/bin/sh
# SPDX-License-Identifier: LGPL-3.0-only
#
# A RECORD, not part of the build, and the only place in this library that names
# `font`: how tests/flatten_vectors.h was made (story 7).  shape-geometry neither
# includes nor links font; this builds a throwaway program, outside any
# repository's build, that runs font's own gfnt_raster_quad and gfnt_raster_cubic.
#
#   tools/vectors-from-font/make.sh <font checkout> <workdir> > tests/flatten_vectors.h
#
# It copies font's src/raster/raster.c into <workdir> with two lines inserted at
# the top of gfnt_raster_line(), which record the end point of each line the
# flattener produces and return (so no scan converter state is needed), links
# the rest of font's src/core, src/outline and src/sfnt by symlink, compiles
# gen.c (which #includes that copy), and runs it.  The vectors were recorded at
# font commit 588437c.  The libfont it links is an installed build of the same
# tree (GFNT symbols the copy needs but does not define); pass FONT_LIB_DIR,
# FONT_GEN_DIR (the font build's generated headers), FONT_LIB (e.g.
# ghoti.io-font-8) and CUTIL_CFLAGS to match the machine.
set -e
FONT=$1
WORK=$2
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$WORK/src/raster"
python3 - "$FONT" "$WORK" <<'PY'
import sys
font, work = sys.argv[1:3]
s = open(font + "/src/raster/raster.c").read()
s = s.replace("#define GFNT_PIXEL_BITS 6", "void scratch_record(GFNT_Point to);\n#define GFNT_PIXEL_BITS 6", 1)
old = "  GFNT_Point from = raster->current;\n  int32_t ey1;"
assert s.count(old) == 1
s = s.replace(old, "  scratch_record(to);\n  if (to.x == to.x) return;\n" + old, 1)
open(work + "/src/raster/raster.c", "w").write(s)
PY
for d in core outline sfnt; do ln -sfn "$FONT/src/$d" "$WORK/src/$d"; done
cp "$HERE/gen.c" "$WORK/gen.c"
cc -std=gnu17 -O1 -w -I"$FONT/include" -I"$FONT_GEN_DIR" $CUTIL_CFLAGS -I"$WORK" "$WORK/gen.c" -o "$WORK/gen" \
  -L"$FONT_LIB_DIR" -Wl,-rpath,"$FONT_LIB_DIR" -l"$FONT_LIB"
"$WORK/gen"
