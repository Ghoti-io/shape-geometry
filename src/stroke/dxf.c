/*
 * SPDX-License-Identifier: LGPL-3.0-only
 *
 * Copyright (C) 2026 Corey Pennycuff
 *
 * This file is part of Ghoti.io Shape Geometry.
 *
 * Ghoti.io Shape Geometry is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3 as
 * published by the Free Software Foundation.
 *
 * Ghoti.io Shape Geometry is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * @file
 *
 * A DXF simple linetype lowered to a dash array and an offset (stroke.h, SG30).
 *
 * The runs of the pattern are written into the caller's array, merged in place:
 * consecutive pen-down elements (dashes and dots) are one run, consecutive spaces
 * one run.  The pattern repeats, so when the first and last runs are of one kind
 * they are one dash (or one space) across the repetition, and the line starts
 * `lead` into it, where `lead` is the length of the last run: the offset.  The
 * runs are rotated to start with a pen-down run.  Two cases keep the dash
 * positions' ties out of the way: a dot first with a dash last keeps the dot apart
 * (a zero-length space between), and a pattern that begins with a space which no
 * trailing space continues and ends with a dot cannot be expressed (the phase would
 * stand exactly at that dot, which the positions keep).
 */

#include <ghoti.io/shape-geometry/macros.h>

#include "stroke_internal.h"

static int finite(double x) {
  return x - x == 0.0;
}

static int in_band(double x) {
  double m = x < 0.0 ? -x : x;
  return m == 0.0 || (m >= GSHGEOM_STROKE_VALUE_MIN && m <= GSHGEOM_STROKE_VALUE_MAX);
}

GSHGEOM_Result gshgeom_stroke_dashes_from_dxf(const double * elements, size_t count, double scale, double * dashes,
    size_t capacity, size_t * dash_count, double * offset) {
  if (!elements || count == 0 || !dashes || !dash_count || !offset || !finite(scale) || !(scale > 0.0)) {
    return GSHGEOM_ERR_INVALID;
  }
  int any_down = 0, any_up = 0;
  double total = 0.0;
  for (size_t i = 0; i < count; i++) {
    if (!finite(elements[i])) {
      return GSHGEOM_ERR_INVALID;
    }
  }
  for (size_t i = 0; i < count; i++) {
    double v = elements[i] * scale;
    if (!finite(v) || !in_band(v)) {
      return GSHGEOM_ERR_LIMIT;
    }
    if (v >= 0.0) {
      any_down = 1;
    } else {
      any_up = 1;
    }
    total += v < 0.0 ? -v : v;
  }
  if (!any_down || !(total > 0.0)) {
    return GSHGEOM_ERR_INVALID; /* draws nothing, or has no length */
  }
  if (!any_up) {
    *dash_count = 0; /* no space: solid */
    *offset = 0.0;
    return GSHGEOM_OK;
  }
  /* how many runs */
  size_t m = 0;
  int prev = -1;
  for (size_t i = 0; i < count; i++) {
    int down = elements[i] * scale >= 0.0;
    if (down != prev) {
      m++;
      prev = down;
    }
  }
  size_t need = m + 1; /* a zero-length space may be appended below */
  if (need > capacity || need > GSHGEOM_STROKE_DASHES_MAX) {
    return GSHGEOM_ERR_LIMIT;
  }
  size_t k = 0;
  int first_down = elements[0] * scale >= 0.0;
  prev = -1;
  for (size_t i = 0; i < count; i++) {
    double v = elements[i] * scale;
    int down = v >= 0.0;
    double len = v < 0.0 ? -v : v;
    if (down != prev) {
      dashes[k++] = len;
      prev = down;
    } else {
      dashes[k - 1] += len;
    }
  }
  int last_down = prev;
  double lead = 0.0;
  if (first_down && dashes[0] == 0.0 && last_down) {
    dashes[k++] = 0.0; /* a dot first, a dash last: the dot is kept apart */
  } else if (k > 1 && first_down == last_down) {
    lead = dashes[k - 1];
    dashes[0] += dashes[k - 1];
    k--;
  }
  double off = lead;
  if (!first_down) {
    /* rotate to start with a pen-down run */
    /* a space first with nothing merged into it and a dot last: the phase would stand exactly at that dot, which
     * SG21 keeps (a dash of no length exactly at the phase), a dot at the start of the line that the pattern does not
     * have */
    if (lead == 0.0 && dashes[k - 1] == 0.0) {
      return GSHGEOM_ERR_UNSUPPORTED;
    }
    double first = dashes[0];
    for (size_t i = 1; i < k; i++) {
      dashes[i - 1] = dashes[i];
    }
    dashes[k - 1] = first;
    double period = 0.0;
    for (size_t i = 0; i < k; i++) {
      period += dashes[i];
    }
    off = lead - first;
    if (off < 0.0) {
      off += period;
    }
  }
  for (size_t i = 0; i < k; i++) {
    if (!in_band(dashes[i])) {
      return GSHGEOM_ERR_LIMIT;
    }
  }
  if (!in_band(off)) {
    return GSHGEOM_ERR_LIMIT;
  }
  *dash_count = k;
  *offset = off;
  return GSHGEOM_OK;
}
