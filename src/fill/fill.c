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
 * The names of the fill rules.
 */

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/shape-geometry/fill.h>

const char * gshgeom_fill_rule_string(GSHGEOM_Fill_Rule rule) {
  switch (rule) {
  case GSHGEOM_FILL_NONZERO:
    return "non-zero";
  case GSHGEOM_FILL_EVEN_ODD:
    return "even-odd";
  case GSHGEOM_FILL_RULE_COUNT:
    break;
  }
  return "unknown fill rule";
}
