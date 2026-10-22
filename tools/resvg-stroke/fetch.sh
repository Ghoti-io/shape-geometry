#!/bin/sh
# SPDX-License-Identifier: LGPL-3.0-only
#
# Copyright (C) 2026 Corey Pennycuff
#
# Build the pinned resvg image (tools/resvg-stroke/Containerfile) when it is not there, and check its version.
#
# Usage: tools/resvg-stroke/fetch.sh
# The engine is $GHOTI_CONTAINER_ENGINE, else podman, else docker.
set -eu
here=$(cd "$(dirname "$0")" && pwd)
image=localhost/ghoti-shape-oracle-resvg:0.45.1
pin=0.45.1
engine="${GHOTI_CONTAINER_ENGINE:-}"
if [ -z "$engine" ]; then
  if command -v podman >/dev/null 2>&1; then engine=podman; elif command -v docker >/dev/null 2>&1; then engine=docker; fi
fi
if [ -z "$engine" ]; then
  echo "resvg-stroke: no container engine (podman or docker)" >&2
  exit 3
fi
if ! $engine image inspect "$image" >/dev/null 2>&1; then
  $engine build -t "$image" -f "$here/Containerfile" "$here"
fi
got=$($engine run --rm "$image" resvg --version | tr -d '\r')
if [ "$got" != "$pin" ]; then
  echo "resvg-stroke: the image says resvg $got, the pin is $pin" >&2
  exit 1
fi
echo "resvg-stroke: $image is resvg $got"
