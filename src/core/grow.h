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
 * Capacity growth shared by the device path and the flattened path.  Internal:
 * not installed, not part of the interface.
 */

#ifndef GHOTI_IO_GSHGEOM_CORE_GROW_H
#define GHOTI_IO_GSHGEOM_CORE_GROW_H

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/shape-geometry/allocator.h>
#include <ghoti.io/shape-geometry/core.h>

#include <stddef.h>
#include <stdint.h>

/**
 * Make room for @p need items of @p size bytes in *array, whose capacity is
 * *cap.  The capacity doubles from 16.  Only the capacity changes: on failure
 * *array and *cap are as they were.
 */
static inline GSHGEOM_Result gshgeom_grow(const GSHGEOM_Allocator * alloc,
    void ** array, size_t * cap, size_t need, size_t size) {
  if (need <= *cap) {
    return GSHGEOM_OK;
  }
  size_t grown = *cap ? *cap : 16;
  while (grown < need) {
    if (grown > SIZE_MAX / 2) {
      return GSHGEOM_ERR_OOM;
    }
    grown *= 2;
  }
  if (grown > SIZE_MAX / size) {
    return GSHGEOM_ERR_OOM;
  }
  void * bigger = gcu_allocator_realloc(alloc, *array, grown * size);
  if (!bigger) {
    return GSHGEOM_ERR_OOM;
  }
  *array = bigger;
  *cap = grown;
  return GSHGEOM_OK;
}

#endif /* GHOTI_IO_GSHGEOM_CORE_GROW_H */
