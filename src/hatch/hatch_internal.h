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
 * The hatch's seam for the tests (internal): the units of the work budget (the limits' `max_work`) that making a set of pieces
 * spent.  Not installed and not part of the interface; the tests link the archive.
 */

#ifndef GHOTI_IO_GSHGEOM_HATCH_INTERNAL_H
#define GHOTI_IO_GSHGEOM_HATCH_INTERNAL_H

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/shape-geometry/hatch.h>

#ifdef __cplusplus
extern "C" {
#endif

/** The work units that making @p pieces spent (the sweeps of its lines and the steps of its dash walks); 0 for NULL. */
uint64_t gshgeom_hatch_pieces_work(const GSHGEOM_Hatch_Pieces * pieces);

#ifdef __cplusplus
}
#endif

#endif /* GHOTI_IO_GSHGEOM_HATCH_INTERNAL_H */
