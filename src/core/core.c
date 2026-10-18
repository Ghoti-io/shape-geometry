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
 * The result strings, the default limits and the version this build reports.
 */

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/shape-geometry/core.h>
#include <ghoti.io/shape-geometry/libver.h>

const char * gshgeom_result_string(GSHGEOM_Result result) {
  switch (result) {
  case GSHGEOM_OK:
    return "ok";
  case GSHGEOM_ERR_IO:
    return "input/output error";
  case GSHGEOM_ERR_FORMAT:
    return "not a format handled here";
  case GSHGEOM_ERR_UNSUPPORTED:
    return "feature not supported";
  case GSHGEOM_ERR_LIMIT:
    return "a limit was exceeded";
  case GSHGEOM_ERR_CORRUPT:
    return "corrupt data";
  case GSHGEOM_ERR_OOM:
    return "out of memory";
  case GSHGEOM_ERR_INVALID:
    return "invalid argument";
  case GSHGEOM_ERR_INTERNAL:
    return "internal error";
  case GSHGEOM_RESULT_COUNT:
    break;
  }
  return "unknown result";
}

void gshgeom_limits_default(GSHGEOM_Limits * limits) {
  if (!limits) {
    return;
  }
  *limits = (GSHGEOM_Limits){0};
  limits->max_elements = GSHGEOM_LIMITS_DEFAULT_MAX_ELEMENTS;
}

const char * gshgeom_version_string(void) {
  return GSHGEOM_VERSION_STRING;
}

unsigned gshgeom_version_number(void) {
  return GSHGEOM_VERSION_NUMBER;
}
