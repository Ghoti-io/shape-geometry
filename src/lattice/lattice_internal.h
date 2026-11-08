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
 * The lattice's seam for the tests (internal): ::gshgeom_lattice_each, also counting the evaluations of a half-plane's sign that the
 * search made.  The count is the work of the search (an estimate far from the answer, an interval already empty and an index that is
 * none are all work), which nothing else reports.  Not installed and not part of the interface; the tests link the archive.
 */

#ifndef GHOTI_IO_GSHGEOM_LATTICE_INTERNAL_H
#define GHOTI_IO_GSHGEOM_LATTICE_INTERNAL_H

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/shape-geometry/lattice.h>

#ifdef __cplusplus
extern "C" {
#endif

/** What the search of ::gshgeom_lattice_each did. */
typedef struct GSHGEOM_Lattice_Stats {
  uint64_t evaluations; /**< Evaluations of a half-plane's sign (a probe of the gallop or the bisection, or a row's whole test). */
} GSHGEOM_Lattice_Stats;

/** ::gshgeom_lattice_each, also filling @p stats (NULL to leave it). */
GSHGEOM_Result gshgeom_lattice_each_counted(const GSHGEOM_Map * map, const GSHGEOM_Lattice * lattice, const int64_t rect[4],
    uint64_t max_count, GSHGEOM_Lattice_Row_Fn row_fn, void * ctx, uint64_t * count, GSHGEOM_Lattice_Stats * stats);

#ifdef __cplusplus
}
#endif

#endif /* GHOTI_IO_GSHGEOM_LATTICE_INTERNAL_H */
