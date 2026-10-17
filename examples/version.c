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
 * Prints what this build is.
 *
 * The smallest program that links the library, so that `make examples` fails
 * if the installed headers and the shared library disagree.
 */

#include <ghoti.io/shape-geometry/shape-geometry.h>

#include <stdio.h>

int main(void) {
  printf("Ghoti.io Shape Geometry %s (0x%06x)\n", gshgeom_version_string(),
      gshgeom_version_number());
  return 0;
}
