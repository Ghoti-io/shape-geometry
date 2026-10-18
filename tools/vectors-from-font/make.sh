#!/bin/sh
# SPDX-License-Identifier: LGPL-3.0-only
#
# A RECORD, not part of the build, and the only place in this library that names
# `font`: how tests/flatten_vectors.h was made (story 7, remade for flattener
# revision 2 in story 7a).  shape-geometry neither includes nor links font; this
# builds a throwaway program, outside any repository's build, that #includes
# font's own src/raster/flatten.c (the unmodified file) and runs its
# gfnt_flatten_quad and gfnt_flatten_cubic at tolerance 4, recording the
# vertices they emit.
#
#   FONT_GEN_DIR=<font build's generated headers> CUTIL_CFLAGS=... \
#     tools/vectors-from-font/make.sh <font checkout> <workdir> > tests/flatten_vectors.h
#
# The vectors were recorded at font commit 69d1860 (flattener revision 2).  Nothing
# of font is linked: flatten.c needs only headers and the inline helpers of
# src/core/fixed.h.  FONT_GEN_DIR is the directory holding
# ghoti.io/font/libver_gen.h (e.g. libs/font/build/linux/release/generated).
set -e
FONT=$1
WORK=$2
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$WORK/src"
for d in core raster; do ln -sfn "$FONT/src/$d" "$WORK/src/$d"; done
cp "$HERE/gen.c" "$WORK/gen.c"
cc -std=gnu17 -O1 -w -I"$FONT/include" -I"$FONT_GEN_DIR" $CUTIL_CFLAGS -I"$WORK" "$WORK/gen.c" -o "$WORK/gen"
"$WORK/gen"
