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
 * One part of a curve outline (stroke_curves.h) checked against the geometric oracle (tests/stroke_oracle.h), as
 * tests/stroke_matrix.h checks a polygon outline: every point the oracle decides (farther than `t + m` from every
 * boundary sample, the samples at most `sigma = t * sigma_over_t` apart and `m` that spacing and a margin for the oracle's own
 * arithmetic) must be classified by the part's nonzero fill as the region classifies it; no ring of the part may
 * wind negatively at any tested point; every sample of the part's curves must lie inside the stroke box plus `t`; and
 * the deviation profile (how far a sample of the part lies outside the region) is measured.  The part's winding is
 * tests/curve_eval.h's, in `long double`, its curves flattened to `t / 64`.
 */

#ifndef GHOTI_IO_GSHGEOM_TESTS_STROKE_CURVES_MATRIX_H
#define GHOTI_IO_GSHGEOM_TESTS_STROKE_CURVES_MATRIX_H

#include "curve_eval.h"
#include "stroke_matrix.h"

#include <chrono>

namespace stroke_curves_matrix {

namespace so = stroke_oracle;
namespace ce = curve_eval;
namespace sm = stroke_matrix;

/** The comparison of one part: see stroke_matrix.h's PartResult. */
typedef sm::PartResult PartResult;

inline PartResult compare_curves(const GSHGEOM_Path * part, const so::Region & region, const sm::Box & box, double t,
    long budget, uint64_t seed, double sigma_over_t = 1.0) {
  PartResult pr;
  auto T0 = std::chrono::steady_clock::now();
  auto lap = [&](const char * what) { if (std::getenv("SCM_TIME")) { auto T1 = std::chrono::steady_clock::now(); std::fprintf(stderr, "  %s %.3f\n", what, std::chrono::duration<double>(T1 - T0).count()); T0 = T1; } };
  so::Real tp = (so::Real)t;
  so::Real sigma = tp * (so::Real)sigma_over_t;
  so::Real m = sigma + 1e-9L;
  so::Samples samples;
  samples.cell = tp + m + 2 * sigma;
  so::sample_region(samples, region, sigma);
  lap("samples");
  ce::Winding wd;
  wd.build(part, tp / 64);
  lap("winding build");
  so::Real bx0 = box.x0 - tp, by0 = box.y0 - tp, bx1 = box.x1 + tp, by1 = box.y1 + tp;
  for (const so::Shape & sh : region.shapes) {
    if (sh.kind != so::Shape::ARCS) {
      continue;
    }
    so::Region one;
    one.local = region.local;
    one.to_device = region.to_device;
    one.h = region.h;
    one.shapes.push_back(sh);
    so::Samples ts;
    ts.cell = 1;
    so::sample_region(ts, one, 1.0L / 16);
    for (const auto & cell : ts.grid) {
      for (const so::V2 & q : cell.second) {
        bx0 = std::min(bx0, q.x - tp);
        by0 = std::min(by0, q.y - tp);
        bx1 = std::max(bx1, q.x + tp);
        by1 = std::max(by1, q.y + tp);
      }
    }
  }
  /* every sample of the part's curves inside the box plus t; the deviation of those outside the region */
  /* at most about 20000 samples of the part's curves, evenly over its elements */
  size_t nel = part ? gshgeom_path_element_count(part) : 0;
  std::vector<ce::P2> cs = ce::sample_path(part, 6, nel > 3000 ? nel / 3000 : 1);
  for (const ce::P2 & q : cs) {
    so::V2 v = {q.x, q.y};
    if (v.x < bx0 - 1e-9L || v.x > bx1 + 1e-9L || v.y < by0 - 1e-9L || v.y > by1 + 1e-9L) {
      if (pr.outside_box++ == 0) {
        std::fprintf(stderr, "outside the box: (%.9Lf, %.9Lf) box [%.6Lf %.6Lf] x [%.6Lf %.6Lf]\n", v.x, v.y, bx0, bx1, by0, by1);
      }
    }
    so::Verdict vd = so::region_contains(region, v);
    if (!vd.inside) {
      so::Real d = samples.nearest(v, samples.cell);
      pr.worst_outward = std::max(pr.worst_outward, d);
    }
  }
  lap("outline samples");
  auto test = [&](so::V2 x) {
    pr.tested++;
    int least = 0;
    int w = wd.at({x.x, x.y}, &least);
    if (w < 0 || least < 0) {
      pr.negative++;
    }
    if (samples.nearest(x, samples.cell) <= tp + m) {
      return;
    }
    so::Verdict v = so::region_contains(region, x);
    if (v.ambiguous) {
      pr.ambiguous++;
      return;
    }
    pr.decided++;
    if (v.inside != (w != 0)) {
      if (pr.wrong++ == 0) {
        char buf[256];
        std::snprintf(buf, sizeof buf, "(%.6Lf, %.6Lf): the region says %s, the outline winds %d", x.x, x.y,
            v.inside ? "inside" : "outside", w);
        pr.first = buf;
      }
    }
  };
  so::Real step = 1.0L / 8;
  so::Real area = (bx1 - bx0) * (by1 - by0);
  if (area / (step * step) > (so::Real)budget) {
    step = std::sqrt(area / (so::Real)budget);
  }
  for (so::Real y = by0; y <= by1; y += step) {
    for (so::Real x = bx0; x <= bx1; x += step) {
      test({x, y});
    }
  }
  /* probes beside the boundary samples (a subsample, so a long boundary does not dominate) */
  size_t stride = samples.count / 600 + 1;
  size_t idx = 0;
  so::Real off = tp + m + sigma;
  for (const auto & cell : samples.grid) {
    for (const so::V2 & s : cell.second) {
      if (idx++ % stride) {
        continue;
      }
      test({s.x + off, s.y});
      test({s.x - off, s.y});
      test({s.x, s.y + off});
      test({s.x, s.y - off});
    }
  }
  lap("grid+probes");
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> ux((double)bx0, (double)bx1), uy((double)by0, (double)by1);
  for (int k = 0; k < 512; k++) {
    test({ux(rng), uy(rng)});
  }
  lap("random");
  return pr;
}

} // namespace stroke_curves_matrix

#endif /* GHOTI_IO_GSHGEOM_TESTS_STROKE_CURVES_MATRIX_H */
