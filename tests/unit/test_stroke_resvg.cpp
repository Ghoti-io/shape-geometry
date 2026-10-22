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
 * The behavior check of the stroke rules' reading (stroke.h): an external
 * renderer's verdicts (tests/stroke_resvg.h, rendered by the pinned resvg with
 * tools/resvg-stroke/gen.py, aliased at 16 pixels per unit) against both the
 * geometric oracle and the expansion's outline, at every sample point the oracle
 * decides with a margin for resvg's own rasterization (two of its pixels) and
 * flattening, plus the expansion's tolerance.  The oracle and the expansion are
 * written from the same text; resvg is not, so a misreading they share shows
 * here.
 */

#include "../test_helpers.h"

#include "../stroke_matrix.h"
#include "../stroke_resvg.h"

#include <string>
#include <vector>

namespace sm = stroke_matrix;
namespace so = stroke_oracle;

namespace {

std::vector<int> verdicts(const stroke_resvg::Case & c) {
  std::vector<int> out;
  int v = 0;
  const char * p = c.runs;
  while (*p) {
    long n = std::strtol(p, const_cast<char **>(&p), 10);
    out.insert(out.end(), (size_t)n, v);
    v = !v;
    if (*p == ',') {
      p++;
    }
  }
  return out;
}

} // namespace

TEST(StrokeResvg, TheRenderersVerdictsAgreeAtEveryDecidedPoint) {
  long decided = 0;
  for (size_t i = 0; i < stroke_resvg::case_count; i++) {
    const stroke_resvg::Case & c = stroke_resvg::cases[i];
    SCOPED_TRACE(c.name);
    sm::Ops ops;
    for (size_t k = 0; k < c.op_count; k++) {
      ops.push_back({c.ops[k].verb, std::vector<double>(c.ops[k].v, c.ops[k].v + c.ops[k].n)});
    }
    GSHGEOM_Path * path = sm::build(ops);
    ASSERT_NE(path, nullptr);
    GSHGEOM_Stroke s;
    gshgeom_stroke_init(&s);
    s.width = c.width;
    s.join = (uint32_t)c.join;
    s.miter_limit = c.limit;
    s.cap = (uint32_t)c.cap;
    GSHGEOM_Map m;
    gshgeom_map_identity(2, &m);
    GSHGEOM_Stroke_Outline * o = nullptr;
    ASSERT_EQ(gshgeom_stroke_expand(path, &s, &m, nullptr, nullptr, &o), GSHGEOM_OK);
    so::Winding w;
    w.build(so::polygons_of(gshgeom_stroke_outline_body(o)));
    so::Stroked sk = so::build_stroke(so::read_centerline(path), sm::oracle_style(s), so::affine_of(m), 0);
    so::Samples samples;
    const so::Real margin = 0.25L; /* our 1/16 px, two resvg pixels at 16x, its flattening, the sampling */
    samples.cell = margin * 2;
    so::sample_region(samples, sk.body, 1.0L / 256);
    std::vector<int> v = verdicts(c);
    ASSERT_EQ(v.size(), (size_t)c.cols * (size_t)c.rows);
    long wrong_ours = 0, wrong_oracle = 0;
    std::string first;
    for (int j = 0; j < c.rows; j++) {
      for (int k = 0; k < c.cols; k++) {
        /* the centre of pixel (k step, j step) at `scale` pixels per unit */
        so::V2 x = {c.box[0] + ((so::Real)(k * stroke_resvg::step) + 0.5L) / stroke_resvg::scale,
            c.box[1] + ((so::Real)(j * stroke_resvg::step) + 0.5L) / stroke_resvg::scale};
        if (samples.nearest(x, samples.cell) <= margin) {
          continue;
        }
        so::Verdict ov = so::region_contains(sk.body, x);
        if (ov.ambiguous) {
          continue;
        }
        decided++;
        int rv = v[(size_t)j * (size_t)c.cols + (size_t)k];
        if ((int)ov.inside != rv) {
          if (wrong_oracle++ == 0) {
            first = "oracle at (" + std::to_string((double)x.x) + ", " + std::to_string((double)x.y) + ")";
          }
        }
        if ((w.at(x) != 0) != (rv != 0)) {
          if (wrong_ours++ == 0 && first.empty()) {
            first = "outline at (" + std::to_string((double)x.x) + ", " + std::to_string((double)x.y) + ")";
          }
        }
      }
    }
    EXPECT_EQ(wrong_oracle, 0) << first;
    EXPECT_EQ(wrong_ours, 0) << first;
    gshgeom_stroke_outline_destroy(o);
    gshgeom_path_destroy(path);
  }
  RecordProperty("resvg_version", stroke_resvg::resvg_version);
  EXPECT_GT(decided, 10000);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
