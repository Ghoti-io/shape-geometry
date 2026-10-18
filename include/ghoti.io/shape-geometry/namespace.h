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
 * @file namespace.h
 *
 * Maps every public name of this library into its version namespace.
 *
 * Kept in one file rather than beside each declaration: a type rename has to
 * be in effect before any struct tag that uses the name, and an internal
 * header may define such a tag without including the public header that
 * declares the typedef.
 *
 * `make check-symbols` fails if an exported symbol is missing from this list.
 *
 * See CONVENTIONS.md section 4.
 */

#ifndef GHOTI_IO_GSHGEOM_NAMESPACE_H
#define GHOTI_IO_GSHGEOM_NAMESPACE_H

#include <ghoti.io/shape-geometry/libver.h>

/// @cond HIDDEN_SYMBOLS

/* Public types. GCU_* names are cutil's; cutil has already renamed them. */
#define GSHGEOM_Allocator GHOTIIO_SHAPE_GEOMETRY(GSHGEOM_Allocator)
#define GSHGEOM_Result GHOTIIO_SHAPE_GEOMETRY(GSHGEOM_Result)
#define GSHGEOM_Limits GHOTIIO_SHAPE_GEOMETRY(GSHGEOM_Limits)
#define GSHGEOM_Hash GHOTIIO_SHAPE_GEOMETRY(GSHGEOM_Hash)
#define GSHGEOM_Curve_Kind GHOTIIO_SHAPE_GEOMETRY(GSHGEOM_Curve_Kind)
#define GSHGEOM_Point GHOTIIO_SHAPE_GEOMETRY(GSHGEOM_Point)
#define GSHGEOM_Curve GHOTIIO_SHAPE_GEOMETRY(GSHGEOM_Curve)
#define GSHGEOM_Bounds GHOTIIO_SHAPE_GEOMETRY(GSHGEOM_Bounds)
#define GSHGEOM_Path GHOTIIO_SHAPE_GEOMETRY(GSHGEOM_Path)
#define GSHGEOM_Verb GHOTIIO_SHAPE_GEOMETRY(GSHGEOM_Verb)
#define GSHGEOM_Element GHOTIIO_SHAPE_GEOMETRY(GSHGEOM_Element)

/* Public functions. */
#define gshgeom_allocator_default GHOTIIO_SHAPE_GEOMETRY(gshgeom_allocator_default)
#define gshgeom_result_string GHOTIIO_SHAPE_GEOMETRY(gshgeom_result_string)
#define gshgeom_limits_default GHOTIIO_SHAPE_GEOMETRY(gshgeom_limits_default)
#define gshgeom_hash_init GHOTIIO_SHAPE_GEOMETRY(gshgeom_hash_init)
#define gshgeom_hash_u8 GHOTIIO_SHAPE_GEOMETRY(gshgeom_hash_u8)
#define gshgeom_hash_u16 GHOTIIO_SHAPE_GEOMETRY(gshgeom_hash_u16)
#define gshgeom_hash_u32 GHOTIIO_SHAPE_GEOMETRY(gshgeom_hash_u32)
#define gshgeom_hash_u64 GHOTIIO_SHAPE_GEOMETRY(gshgeom_hash_u64)
#define gshgeom_hash_i32 GHOTIIO_SHAPE_GEOMETRY(gshgeom_hash_i32)
#define gshgeom_hash_i64 GHOTIIO_SHAPE_GEOMETRY(gshgeom_hash_i64)
#define gshgeom_hash_f64 GHOTIIO_SHAPE_GEOMETRY(gshgeom_hash_f64)
#define gshgeom_hash_bytes GHOTIIO_SHAPE_GEOMETRY(gshgeom_hash_bytes)
#define gshgeom_hash_string GHOTIIO_SHAPE_GEOMETRY(gshgeom_hash_string)
#define gshgeom_hash_value GHOTIIO_SHAPE_GEOMETRY(gshgeom_hash_value)
#define gshgeom_f64_bits GHOTIIO_SHAPE_GEOMETRY(gshgeom_f64_bits)
#define gshgeom_curve_point_count GHOTIIO_SHAPE_GEOMETRY(gshgeom_curve_point_count)
#define gshgeom_curve_is_rational GHOTIIO_SHAPE_GEOMETRY(gshgeom_curve_is_rational)
#define gshgeom_curve_line GHOTIIO_SHAPE_GEOMETRY(gshgeom_curve_line)
#define gshgeom_curve_quadratic GHOTIIO_SHAPE_GEOMETRY(gshgeom_curve_quadratic)
#define gshgeom_curve_cubic GHOTIIO_SHAPE_GEOMETRY(gshgeom_curve_cubic)
#define gshgeom_curve_rational_quadratic GHOTIIO_SHAPE_GEOMETRY(gshgeom_curve_rational_quadratic)
#define gshgeom_curve_rational_cubic GHOTIIO_SHAPE_GEOMETRY(gshgeom_curve_rational_cubic)
#define gshgeom_curve_check GHOTIIO_SHAPE_GEOMETRY(gshgeom_curve_check)
#define gshgeom_curve_eval GHOTIIO_SHAPE_GEOMETRY(gshgeom_curve_eval)
#define gshgeom_curve_bounds GHOTIIO_SHAPE_GEOMETRY(gshgeom_curve_bounds)
#define gshgeom_path_create GHOTIIO_SHAPE_GEOMETRY(gshgeom_path_create)
#define gshgeom_path_destroy GHOTIIO_SHAPE_GEOMETRY(gshgeom_path_destroy)
#define gshgeom_path_dimension GHOTIIO_SHAPE_GEOMETRY(gshgeom_path_dimension)
#define gshgeom_path_element_count GHOTIIO_SHAPE_GEOMETRY(gshgeom_path_element_count)
#define gshgeom_path_move_to GHOTIIO_SHAPE_GEOMETRY(gshgeom_path_move_to)
#define gshgeom_path_line_to GHOTIIO_SHAPE_GEOMETRY(gshgeom_path_line_to)
#define gshgeom_path_quadratic_to GHOTIIO_SHAPE_GEOMETRY(gshgeom_path_quadratic_to)
#define gshgeom_path_cubic_to GHOTIIO_SHAPE_GEOMETRY(gshgeom_path_cubic_to)
#define gshgeom_path_rational_quadratic_to GHOTIIO_SHAPE_GEOMETRY(gshgeom_path_rational_quadratic_to)
#define gshgeom_path_rational_cubic_to GHOTIIO_SHAPE_GEOMETRY(gshgeom_path_rational_cubic_to)
#define gshgeom_path_close GHOTIIO_SHAPE_GEOMETRY(gshgeom_path_close)
#define gshgeom_path_element GHOTIIO_SHAPE_GEOMETRY(gshgeom_path_element)
#define gshgeom_path_bounds GHOTIIO_SHAPE_GEOMETRY(gshgeom_path_bounds)
#define gshgeom_path_hash GHOTIIO_SHAPE_GEOMETRY(gshgeom_path_hash)
#define gshgeom_path_dump GHOTIIO_SHAPE_GEOMETRY(gshgeom_path_dump)
#define gshgeom_version_number GHOTIIO_SHAPE_GEOMETRY(gshgeom_version_number)
#define gshgeom_version_string GHOTIIO_SHAPE_GEOMETRY(gshgeom_version_string)

/// @endcond

#endif // GHOTI_IO_GSHGEOM_NAMESPACE_H
