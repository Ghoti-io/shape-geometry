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
 * The body of a segment (SG6) as pieces (emit.c), and the bounds that make each
 * piece's `d_in` and `d_out`.
 *
 * A line's body is the rectangle of its corners `P0 +- h n`, `P1 +- h n`: one
 * piece, exact.  A curve is split by de Casteljau, in homogeneous coordinates
 * relative to its start, at the parameter 7/16 (so that no split point lands on
 * 1/2, 1/4 or any other parameter of the form `odd / 2^k`, where a symmetric
 * construction puts a cusp: every split point has an odd numerator over a power
 * of 16), until each piece is accepted as a **band** or as a **fan**, at most
 * ::GSHGEOM_STROKE_MAX_DEPTH levels deep; a piece still refused there is
 * ::GSHGEOM_ERR_LIMIT, never drawn coarser.
 *
 * ## What a piece knows
 *
 * For a piece with homogeneous control points `N_i = w_i (P_i - o)`, `D_i = w_i`
 * (each with an error bound carried through the splits, emit.c (a)), the
 * direction of `c'(u)` is that of `E(u) = N'D - ND'`, a polynomial of degree
 * `2n - 1` whose Bernstein coefficients `E_k` are formed from the `N_i`, `D_i`
 * (each with its error bound `e_k`).  `E(u)` is a positive combination of the
 * `E_k`, so its direction lies in their cone.  A coefficient that is exactly
 * zero with a zero error bound adds no direction; one within four times its
 * error of zero makes the cone unknown.  The cone's width is bounded from the
 * pairwise angles of the coefficients, each widened by `asin(e_k / |E_k|)`.
 * The curvature `|E x E'| D^2 / |E|^3` is at most
 * `max|E_k| max|E'_k| (max D_i)^2 / m^3`, `m` the least projection of a
 * coefficient on a coefficient's direction (all of it in build space).
 *
 * ## A band
 *
 * Accepted when the build-space cone and the device cone (the coefficients
 * mapped by the map's linear part) are each narrower than a quarter turn, the
 * curvature bound `kappa` satisfies `kappa h < cos(alpha)` (`alpha` the
 * build-space cone's width), and the four device corners `C0 +- H0`, `C1 +- H1`
 * (`H = h lin n`, `n` the unit left normal of the end coefficients) are a strictly
 * convex quadrilateral.  Then:
 *  - **No two normal segments of the piece meet.**  Two normals of a curve whose
 *    tangent turns within a cone of width `alpha < pi/2` meet at a distance from
 *    each foot of at least `cos(alpha) / kappa > h` (the distance is
 *    `int cos(theta - theta_2) ds / sin(theta_2 - theta_1)`, and
 *    `|sin(theta_2 - theta_1)| <= kappa s`).  So the band is the image of a
 *    rectangle under an injective map: a disc bounded by the two rails, the two
 *    end normals.
 *  - **Each rail is near its chord.**  A rail `c +- h n` has the tangent
 *    `c' (1 -+ h kappa_signed)`, the factor positive by the bound above, so its
 *    tangents lie in the device cone; a curve whose tangents lie in a cone of
 *    width `alpha_d < pi/2` containing its chord's direction is a graph over the
 *    chord within `(lambda / 2) sin(alpha_d)` of it, `lambda` its length, and
 *    `lambda <= (1 + h kappa) Lpoly`, `Lpoly` the device control polygon's length
 *    (a curve with positive weights is no longer than its control polygon).
 *  - So the band and the quadrilateral have boundaries within
 *    `delta = (lambda / 2) sin(alpha_d)` plus the corners' error of each other,
 *    point for point, and each region is within that of the other (a point of
 *    either far from both boundaries has the same winding in both):
 *    `d_in = d_out = delta + corners`.
 *
 * ## A fan
 *
 * When a band is not accepted and the piece is short: `dm` bounds the device
 * distance from the piece's start point to its control points (the piece lies
 * in their hull), and `n0`, `n1` are the normals at the ends (the first and the
 * last coefficient that is not exactly zero: one-sided limits where `c'`
 * vanishes).  The normal segments of the piece are within `dm` of the segments
 * from the start point along the same directions:
 *  - **a double sector**, when the cone is narrower than a half turn: the
 *    normals lie in the cone (widened), and the true normal turns continuously
 *    from `n0` to `n1` inside it, sweeping the arc between them, so the sector
 *    over the cone is within `dm` plus `h` times the chord of its overshoot past
 *    `n0` and `n1` of the region;
 *  - **a disc**, when `n0` and `n1` are nearly opposite: either `c'` vanishes
 *    inside (a cusp, whose disc is SG6's) or the normal turns continuously by at
 *    least the angle between them, and the full normal segments then sweep the
 *    disc but for a double wedge whose chord is `|n0 + n1|`.
 * In both `d_in = dm` (and the ends' uncertainty) and `d_out` adds the overshoot.
 *
 * ## Curve mode (stroke_curves.h, SG37, SG38)
 *
 * With the context's `cv` set the same recursion writes curves instead of polygons, in target units.  A piece is a
 * **band** when band_analyze() says its normals do not meet: it is then bounded by two offset cubics, the cubic
 * Hermite interpolants of `c +- h n` at the piece's ends, whose error (the supremum of the fourth derivative of the
 * offset over 384) is certified by fourth_bound() from Bernstein bounds of the derivatives of `c` (for a rational
 * piece through `D c = N`) and the recurrences for the unit tangent's derivatives.  Consecutive bands of a segment are
 * one contour (the shared end normals cancel and are not written), their shared end points and tangents computed once
 * where the piece was split.  A piece whose band is not simple only because the offset folds (the curvature bound alone
 * refuses) is written as polygons by chord_curve(): the quadrilateral of its two end normals, or the two lobes when
 * they cross, within an eighth of the supremum of the second derivative of `c + lambda n` of the region its normals
 * sweep, second order in its extent however the offset folds.  A piece the certificate cannot reach at all (a cusp, a
 * stopped end) is covered by the fan below, in curve mode measured in build space and inflated by its own displacement
 * so that neighbouring fans overlap, down to the budget's extent.  Everything else (lines, joins, caps) is the polygon
 * mode's, written through gshgeom_sk_emit() as lines and contour.c's exact arcs.
 */

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/math/basic.h>
#include <ghoti.io/shape-geometry/orient.h>

#include "stroke_internal.h"

#include <string.h>

/** The split parameter and its complement: exact doubles. */
#define SPLIT 0.4375
#define SPLIT_C 0.5625

/** The unit roundoff. */
#define U 0x1p-53

/** How much a computed bound is raised for the rounding of its own computation. */
#define UP (1.0 + 0x1p-40)

/** A piece: homogeneous control points relative to the segment's start, with per-point error bounds. */
typedef struct Piece {
  double nx[4], ny[4], d[4];
  double en[4], ed[4]; ///< bounds of the error of each coordinate of N_i, and of D_i
  int z0;              ///< the first control point truly equals the second (point and weight): E's terms there vanish
  int zn;              ///< the last control point truly equals the one before it
  /* curve mode (SG37): the tangent directions at the piece's two ends, each computed once where the piece's end is a
   * split point and shared by the piece on either side of it, with a bound of the angle they may be off by */
  GSHGEOM_SK_Vec dir0, dir1;
  double bd0, bd1;
  int has0, has1;      ///< the direction is known to within its bound
} Piece;

static double fabs_(double x) {
  return x < 0 ? -x : x;
}

static double max2(double a, double b) {
  return a > b ? a : b;
}

static double mag(double x, double y) {
  return max2(fabs_(x), fabs_(y));
}

/** Split p at SPLIT into left and right (degree n). */
static void split(const Piece * p, int n, Piece * l, Piece * r) {
  double nx[4], ny[4], d[4], en[4], ed[4];
  memcpy(nx, p->nx, sizeof nx);
  memcpy(ny, p->ny, sizeof ny);
  memcpy(d, p->d, sizeof d);
  memcpy(en, p->en, sizeof en);
  memcpy(ed, p->ed, sizeof ed);
  memset(l, 0, sizeof *l);
  memset(r, 0, sizeof *r);
  /* a coincident end pair stays coincident in the half that keeps that end (de Casteljau of equal points is exact in
   * the true values), and nowhere else */
  l->z0 = p->z0;
  r->zn = p->zn;
  l->nx[0] = nx[0];
  l->ny[0] = ny[0];
  l->d[0] = d[0];
  l->en[0] = en[0];
  l->ed[0] = ed[0];
  r->nx[n] = nx[n];
  r->ny[n] = ny[n];
  r->d[n] = d[n];
  r->en[n] = en[n];
  r->ed[n] = ed[n];
  l->dir0 = p->dir0;
  l->bd0 = p->bd0;
  l->has0 = p->has0;
  r->dir1 = p->dir1;
  r->bd1 = p->bd1;
  r->has1 = p->has1;
  for (int lev = 1; lev <= n; lev++) {
    for (int i = 0; i + lev <= n; i++) {
      double ax = nx[i], bx = nx[i + 1], ay = ny[i], by = ny[i + 1], ad = d[i], bd = d[i + 1];
      double ma = mag(ax, ay), mb = mag(bx, by);
      nx[i] = SPLIT_C * ax + SPLIT * bx;
      ny[i] = SPLIT_C * ay + SPLIT * by;
      d[i] = SPLIT_C * ad + SPLIT * bd;
      /* a convex combination of two values with their errors, and three roundings of the result's size */
      en[i] = SPLIT_C * en[i] + SPLIT * en[i + 1] + 3.0 * U * (SPLIT_C * ma + SPLIT * mb);
      ed[i] = SPLIT_C * ed[i] + SPLIT * ed[i + 1] + 3.0 * U * (SPLIT_C * fabs_(ad) + SPLIT * fabs_(bd));
    }
    if (lev == n - 1) {
      /* the last two points before the split point: the curve's tangent there points from the first to the second
       * (their projections), the one direction both halves use at their shared end */
      double p0x = nx[0] / d[0], p0y = ny[0] / d[0], p1x = nx[1] / d[1], p1y = ny[1] / d[1];
      GSHGEOM_SK_Vec dv = gshgeom_sk_v(p1x - p0x, p1y - p0y);
      double e0 = (en[0] + mag(p0x, p0y) * ed[0]) / fabs_(d[0]) + 2.0 * U * mag(p0x, p0y);
      double e1 = (en[1] + mag(p1x, p1y) * ed[1]) / fabs_(d[1]) + 2.0 * U * mag(p1x, p1y);
      double len = gshgeom_sk_norm(dv);
      double beta = len > 0.0 ? 4.5 * (e0 + e1) / len * UP : 1.0;
      int ok = len > 0.0 && beta <= 0.0625 && dv.x - dv.x == 0.0 && dv.y - dv.y == 0.0;
      l->dir1 = r->dir0 = dv;
      l->bd1 = r->bd0 = ok ? (beta > 0x1p-40 ? beta : 0x1p-40) : 1.0;
      l->has1 = r->has0 = ok;
    }
    l->nx[lev] = nx[0];
    l->ny[lev] = ny[0];
    l->d[lev] = d[0];
    l->en[lev] = en[0];
    l->ed[lev] = ed[0];
    r->nx[n - lev] = nx[n - lev];
    r->ny[n - lev] = ny[n - lev];
    r->d[n - lev] = d[n - lev];
    r->en[n - lev] = en[n - lev];
    r->ed[n - lev] = ed[n - lev];
  }
}

/** n C(n-1, i) C(n, j) / C(2n-1, i+j): the weights of E's coefficients. */
static double ecoef(int n, int i, int j) {
  static const double binom[6][6] = {{1, 0, 0, 0, 0, 0}, {1, 1, 0, 0, 0, 0}, {1, 2, 1, 0, 0, 0}, {1, 3, 3, 1, 0, 0},
      {1, 4, 6, 4, 1, 0}, {1, 5, 10, 10, 5, 1}};
  return (double)n * binom[n - 1][i] * binom[n][j] / binom[2 * n - 1][i + j];
}

/** What a piece's coefficients tell: the cone and the bounds (build space unless named device). */
typedef struct Know {
  int m;                 ///< number of coefficients (2n)
  GSHGEOM_SK_Vec e[6];   ///< E_k
  double err[6];         ///< bound of |E_k - true E_k| (a vector bound)
  int zero[6];           ///< exactly zero, error zero: no direction
  int unknown;           ///< a coefficient too near zero to know its direction
} Know;

static void know(const Piece * p, int n, Know * k) {
  memset(k, 0, sizeof *k);
  k->m = 2 * n;
  for (int i = 0; i < n; i++) {
    if ((i == 0 && p->z0) || (i == n - 1 && p->zn)) {
      continue; /* the difference of two truly equal points is exactly zero, with no error: so are its terms */
    }
    for (int j = 0; j <= n; j++) {
      double cf = ecoef(n, i, j);
      double dnx = p->nx[i + 1] - p->nx[i], dny = p->ny[i + 1] - p->ny[i];
      double dd = p->d[i + 1] - p->d[i];
      double tx = dnx * p->d[j] - dd * p->nx[j];
      double ty = dny * p->d[j] - dd * p->ny[j];
      int kk = i + j;
      k->e[kk].x += cf * tx;
      k->e[kk].y += cf * ty;
      double sn = mag(p->nx[i + 1], p->ny[i + 1]) + mag(p->nx[i], p->ny[i]);
      double sd = fabs_(p->d[i + 1]) + fabs_(p->d[i]);
      double en2 = p->en[i + 1] + p->en[i], ed2 = p->ed[i + 1] + p->ed[i];
      double t = en2 * (fabs_(p->d[j]) + p->ed[j]) + sn * p->ed[j] + ed2 * (mag(p->nx[j], p->ny[j]) + p->en[j]) +
          sd * p->en[j] + 8.0 * U * (sn * fabs_(p->d[j]) + sd * mag(p->nx[j], p->ny[j]));
      k->err[kk] += cf * t;
    }
  }
  for (int i = 0; i < k->m; i++) {
    k->err[i] = k->err[i] * 1.5 * UP; /* a coordinate bound to a vector bound (sqrt 2), and the sums' rounding */
    double len = gshgeom_sk_norm(k->e[i]);
    if (len == 0.0 && k->err[i] == 0.0) {
      k->zero[i] = 1;
    } else if (k->err[i] * 4.0 > len) {
      k->unknown = 1;
    }
  }
}

/** The angular widening of a known coefficient: asin(err / |E|) <= 1.6 err / |E| for err / |E| <= 1/4. */
static double widen(const Know * k, int i) {
  double b = 1.6 * k->err[i] / gshgeom_sk_norm(k->e[i]) * UP;
  return b > 0x1p-40 ? b : 0x1p-40; /* never exactly zero: a cone always has a width to turn by */
}

/**
 * The cone of a set of vectors with angular uncertainties beta: lower bound of the cosine and upper bound of the sine
 * of its width, from the pairwise angles; ok when every widened pair is within a quarter turn.
 */
static int quarter_cone(const GSHGEOM_SK_Vec * v, const double * beta, int m, double * cosw, double * sinw) {
  double cmin = 1.0, smax = 0.0;
  int any = 0;
  for (int i = 0; i < m; i++) {
    for (int j = i; j < m; j++) {
      double c = gshgeom_sk_dot(v[i], v[j]) - (beta[i] + beta[j]);
      double s = fabs_(gshgeom_sk_cross(v[i], v[j])) + beta[i] + beta[j];
      cmin = c < cmin ? c : cmin;
      smax = s > smax ? s : smax;
      any = 1;
    }
  }
  if (!any || !(cmin > 0.0)) {
    return 0;
  }
  *cosw = cmin;
  *sinw = smax < 1.0 ? smax : 1.0;
  return 1;
}

/** The recursion's state for one segment. */
typedef struct Seg {
  GSHGEOM_SK_Ctx * c;
  int n;
  GSHGEOM_SK_Vec anchor; ///< the device position of the segment's start
} Seg;

/** The unit vectors (build space) and widenings of the known nonzero coefficients. */
static int units(const Know * k, GSHGEOM_SK_Vec * u, double * beta) {
  int m = 0;
  for (int i = 0; i < k->m; i++) {
    if (k->zero[i]) {
      continue;
    }
    u[m] = gshgeom_sk_unit(k->e[i]);
    beta[m] = widen(k, i);
    m++;
  }
  return m;
}


/** What the band's tests learn about a piece (SG6, the quarter cones and the curvature bound; see the file comment). */
typedef struct BandInfo {
  GSHGEOM_SK_Vec u[6];     ///< the unit directions of the known coefficients (build space)
  double beta[6];          ///< their widenings
  int m;                   ///< how many
  double cosb, sinb;       ///< the build-space cone: a lower bound of the cosine, an upper bound of the sine of its width
  double cosd, sind;       ///< the same for the device cone
  double kappa;            ///< an upper bound of the curvature
  double minproj;          ///< the least projection of a coefficient on u[0]
  double maxd;             ///< an upper bound of the weight function
  GSHGEOM_SK_Vec P[4];     ///< the piece's control points, relative to the segment's start (build space)
  double ep;               ///< an upper bound of the error of each coordinate of P
} BandInfo;

/**
 * The tests of a band: 1 when the piece's normals do not meet and its cones are a quarter turn, with the numbers; 2 when
 * only the curvature bound refuses (the normals may meet within h: a fold, the numbers set all the same); 0 otherwise.
 */
static int band_analyze(const Seg * sg, const Piece * p, const Know * k, BandInfo * bi) {
  GSHGEOM_SK_Ctx * c = sg->c;
  int n = sg->n;
  if (k->unknown) {
    return 0;
  }
  for (int i = 0; i < k->m; i++) {
    if (k->zero[i]) {
      return 0; /* |E| reaches zero: the curvature is not bounded */
    }
  }
  GSHGEOM_SK_Vec ud[6];
  double betad[6];
  int m = units(k, bi->u, bi->beta);
  bi->m = m;
  if (!quarter_cone(bi->u, bi->beta, m, &bi->cosb, &bi->sinb)) {
    return 0;
  }
  /* the device cone: the coefficients mapped, each with its error times the norm and the mapping's rounding */
  for (int i = 0; i < m; i++) {
    GSHGEOM_SK_Vec ed = gshgeom_sk_lin(c, k->e[i]);
    double len = gshgeom_sk_norm(ed);
    double er = (k->err[i] * c->norm + 4.0 * U * c->norm * gshgeom_sk_norm(k->e[i])) * UP;
    if (!(len > 4.0 * er)) {
      return 0;
    }
    ud[i] = gshgeom_sk_unit(ed);
    betad[i] = 1.6 * er / len * UP;
  }
  if (!quarter_cone(ud, betad, m, &bi->cosd, &bi->sind)) {
    return 0;
  }
  /* the least projection on u[0], and the curvature bound */
  double minproj = -1.0, maxe = 0.0, maxde = 0.0, maxd = 0.0;
  for (int i = 0; i < k->m; i++) {
    double pr = gshgeom_sk_dot(k->e[i], bi->u[0]) - k->err[i];
    minproj = (minproj < 0.0 || pr < minproj) ? pr : minproj;
    maxe = max2(maxe, gshgeom_sk_norm(k->e[i]) + k->err[i]);
  }
  if (!(minproj > 0.0)) {
    return 0;
  }
  for (int i = 0; i + 1 < k->m; i++) {
    GSHGEOM_SK_Vec de = gshgeom_sk_sub(k->e[i + 1], k->e[i]);
    double v = (double)(k->m - 1) * (gshgeom_sk_norm(de) + k->err[i] + k->err[i + 1]);
    maxde = max2(maxde, v);
  }
  for (int i = 0; i <= n; i++) {
    maxd = max2(maxd, p->d[i] + p->ed[i]);
  }
  double kappa = maxe * maxde * maxd * maxd / (minproj * minproj * minproj) * UP;
  int fold = !(kappa * c->h < bi->cosb); /* the normals may meet within h: not a band */
  bi->kappa = kappa;
  bi->minproj = minproj;
  bi->maxd = maxd;
  /* the projected control points (build space, relative to the segment's start) and their error */
  double ep = 0.0;
  for (int i = 0; i <= n; i++) {
    bi->P[i] = gshgeom_sk_v(p->nx[i] / p->d[i], p->ny[i] / p->d[i]);
    double e = (p->en[i] + mag(bi->P[i].x, bi->P[i].y) * p->ed[i]) / p->d[i] + 2.0 * U * mag(bi->P[i].x, bi->P[i].y);
    ep = max2(ep, e * 1.5);
  }
  bi->ep = ep;
  return fold ? 2 : 1;
}

/** Try the band as a quadrilateral (polygon mode); returns 1 when emitted (or on an error, in *r), 0 to try something else. */
static int band(const Seg * sg, const Piece * p, const Know * k, GSHGEOM_Result * r) {
  GSHGEOM_SK_Ctx * c = sg->c;
  int n = sg->n;
  BandInfo bi;
  if (band_analyze(sg, p, k, &bi) != 1) {
    return 0;
  }
  GSHGEOM_SK_Vec * P = bi.P;
  double ep = bi.ep, kappa = bi.kappa, sind = bi.sind;
  GSHGEOM_SK_Vec * u = bi.u;
  double * beta = bi.beta;
  int m = bi.m;
  double lpoly = 0.0;
  for (int i = 0; i < n; i++) {
    lpoly += gshgeom_sk_norm(gshgeom_sk_lin(c, gshgeom_sk_sub(P[i + 1], P[i])));
  }
  lpoly = (lpoly + 2.0 * (double)n * c->norm * ep) * UP;
  double lambda = (1.0 + c->h * kappa) * lpoly * UP;
  double delta = 0.5 * lambda * sind * UP;
  /* the corners */
  GSHGEOM_SK_Vec n0 = gshgeom_sk_left(u[0]), n1 = gshgeom_sk_left(u[m - 1]);
  double corner = c->norm * (ep + c->h * max2(beta[0], beta[m - 1]) * 1.1) * UP;
  GSHGEOM_SK_Vec q[4] = {gshgeom_sk_add(P[0], gshgeom_sk_scale(c->h, n0)), gshgeom_sk_add(P[n], gshgeom_sk_scale(c->h, n1)),
      gshgeom_sk_sub(P[n], gshgeom_sk_scale(c->h, n1)), gshgeom_sk_sub(P[0], gshgeom_sk_scale(c->h, n0))};
  /* strictly convex in device space, exactly */
  GSHGEOM_Point dq[4];
  for (int i = 0; i < 4; i++) {
    GSHGEOM_SK_Vec dv = gshgeom_sk_add(sg->anchor, gshgeom_sk_lin(c, q[i]));
    dq[i].x = dv.x;
    dq[i].y = dv.y;
    dq[i].z = 0.0;
  }
  int first = 0;
  for (int i = 0; i < 4; i++) {
    int s = 0;
    (void)gshgeom_orient2d(dq[i], dq[(i + 1) % 4], dq[(i + 2) % 4], &s);
    if (s == 0 || (first != 0 && s != first)) {
      return 0;
    }
    first = s;
  }
  double d = delta + corner;
  if (!gshgeom_sk_fits_conform(c, d, d)) {
    return 0;
  }
  /* the centerline points at the middle of the two end edges, so that the band beside this one, or the join or cap, shares
   * each end edge with it exactly */
  GSHGEOM_SK_Vec piv[2] = {P[0], P[n]};
  *r = gshgeom_sk_emit_conform(c, sg->anchor, q, 4, d, d, piv, 2);
  return 1;
}

/* --- curve mode (SG37, SG38) ----------------------------------------------------------- */

static double binom_c(int n, int k) {
  static const double b[6][6] = {{1, 0, 0, 0, 0, 0}, {1, 1, 0, 0, 0, 0}, {1, 2, 1, 0, 0, 0}, {1, 3, 3, 1, 0, 0},
      {1, 4, 6, 4, 1, 0}, {1, 5, 10, 10, 5, 1}};
  return b[n][k];
}

/** An upper bound of the norm of the k-th derivative of a Bernstein polynomial of degree n, from its control values. */
static double deriv_bound(int n, int k, const double * vx, const double * vy, const double * errc) {
  if (k > n) {
    return 0.0;
  }
  double x[4], y[4];
  memcpy(x, vx, sizeof x);
  memcpy(y, vy, sizeof y);
  double emax = 0.0;
  for (int i = 0; i <= n; i++) {
    emax = max2(emax, errc[i]);
  }
  for (int lev = 0; lev < k; lev++) {
    for (int i = 0; i + 1 <= n - lev; i++) {
      x[i] = x[i + 1] - x[i];
      y[i] = y[i + 1] - y[i];
    }
  }
  double scale = 1.0;
  for (int i = 0; i < k; i++) {
    scale *= (double)(n - i);
  }
  double best = 0.0;
  for (int i = 0; i + k <= n; i++) {
    best = max2(best, gshgeom_sk_norm(gshgeom_sk_v(x[i], y[i])));
  }
  /* the error of each control value, k differences deep: 2^k of them */
  double ek = emax * (double)(1 << k) * 1.5;
  return (best + ek) * scale * UP;
}

/** sup |c^(j)|, j = 0..5, by D c = N, and the lower bound m of |c'|; 0 when the piece cannot be certified. */
static int bounds_g(const Piece * p, int n, const BandInfo * bi, double * G, double * m) {
  double dmin = 1e300;
  for (int i = 0; i <= n; i++) {
    dmin = p->d[i] - p->ed[i] < dmin ? p->d[i] - p->ed[i] : dmin;
  }
  *m = bi->minproj / (bi->maxd * bi->maxd) / UP;
  if (!(dmin > 0.0) || !(*m > 0.0)) {
    return 0;
  }
  double BN[6], BD[6];
  double enc[4], edc[4];
  for (int i = 0; i <= n; i++) {
    enc[i] = p->en[i];
    edc[i] = p->ed[i];
  }
  for (int j = 0; j <= 5; j++) {
    BN[j] = deriv_bound(n, j, p->nx, p->ny, enc);
    double zero[4] = {0, 0, 0, 0};
    BD[j] = deriv_bound(n, j, p->d, zero, edc);
  }
  double g0 = 0.0;
  for (int i = 0; i <= n; i++) {
    g0 = max2(g0, gshgeom_sk_norm(bi->P[i]));
  }
  G[0] = g0 + bi->ep;
  for (int j = 1; j <= 5; j++) {
    double sum = BN[j];
    for (int i = 0; i < j; i++) {
      sum += binom_c(j, i) * BD[j - i] * G[i];
    }
    G[j] = sum / dmin * UP;
  }
  return 1;
}

/*
 * The certified bound of SG38: with G[j] an upper bound of sup |c^(j)| on the piece (G[0] of the distance from its
 * start) and m a lower bound of |c'|, an upper bound of sup |O''''| for the offsets c +- h n, namely G[4] + h F[4],
 * F[k] bounding the k-th derivative of the unit tangent T = c' / |c'|.  With s = |c'|, g = c', g_j its j-th derivative
 * (G[j+1]): s T = g gives s T^(k) = g^(k) - sum_{i<k} C(k,i) s^(k-i) T^(i), and s^2 = g . g gives
 * 2 s s^(k) = sum_{i<=k} C(k,i) g^(i) . g^(k-i) - sum_{0<i<k} C(k,i) s^(i) s^(k-i), so |s^(k)| <= S[k] and
 * |T^(k)| <= F[k] follow term by term (S[1] = g_1, since |s'| = |g . g'| / |g| <= |g'|; F[1] = g_1 / m, since
 * T' = g'_perp / |g|).  Every quantity is a sum of products of nonnegative bounds, raised by UP at each step.
 */
static void tangent_bounds(const double * G, double m, double * F) {
  double g[5], S[5];
  for (int j = 0; j < 5; j++) {
    g[j] = G[j + 1];
  }
  S[1] = g[1];
  F[0] = 1.0;
  F[1] = g[1] / m * UP;
  for (int k = 2; k <= 4; k++) {
    double sum = 0.0;
    for (int i = 0; i <= k; i++) {
      sum += binom_c(k, i) * g[i] * g[k - i];
    }
    for (int i = 1; i < k; i++) {
      sum += binom_c(k, i) * S[i] * S[k - i];
    }
    S[k] = sum / (2.0 * m) * UP;
    double f = g[k];
    for (int i = 0; i < k; i++) {
      f += binom_c(k, i) * S[k - i] * F[i];
    }
    F[k] = f / m * UP;
  }
}

/** An upper bound of sup of the fourth derivative of the offsets c +- h n (SG38). */
static double fourth_bound(const double * G, double m, double h) {
  double F[5];
  tangent_bounds(G, m, F);
  return (G[4] + h * F[4]) * UP;
}

/** An upper bound of sup |F_uu| for F(u, lambda) = c(u) + lambda n(u), |lambda| <= h: `c'' + lambda n''` (SG37). */
static double second_bound(const double * G, double m, double h) {
  double F[5];
  tangent_bounds(G, m, F);
  return (G[2] + h * F[2]) * UP;
}

/** The first and second derivatives of the piece's curve at the start (end 0) or the end (end 1), relative to its start. */
static void end_derivs(const Piece * p, int n, int end, GSHGEOM_SK_Vec * c1, GSHGEOM_SK_Vec * c2) {
  double nx[4], ny[4], d[4];
  for (int i = 0; i <= n; i++) {
    int j = end ? n - i : i;
    nx[i] = p->nx[j];
    ny[i] = p->ny[j];
    d[i] = p->d[j];
  }
  double nn = (double)n, n1 = (double)(n * (n - 1));
  GSHGEOM_SK_Vec c0 = gshgeom_sk_v(nx[0] / d[0], ny[0] / d[0]);
  GSHGEOM_SK_Vec N1 = gshgeom_sk_v(nn * (nx[1] - nx[0]), nn * (ny[1] - ny[0]));
  double D1 = nn * (d[1] - d[0]);
  GSHGEOM_SK_Vec v1 = gshgeom_sk_scale(1.0 / d[0], gshgeom_sk_sub(N1, gshgeom_sk_scale(D1, c0)));
  GSHGEOM_SK_Vec N2 = gshgeom_sk_v(n1 * (nx[2] - 2.0 * nx[1] + nx[0]), n1 * (ny[2] - 2.0 * ny[1] + ny[0]));
  double D2 = n1 * (d[2] - 2.0 * d[1] + d[0]);
  GSHGEOM_SK_Vec v2 = gshgeom_sk_scale(1.0 / d[0],
      gshgeom_sk_sub(gshgeom_sk_sub(N2, gshgeom_sk_scale(2.0 * D1, v1)), gshgeom_sk_scale(D2, c0)));
  *c1 = end ? gshgeom_sk_scale(-1.0, v1) : v1;
  *c2 = v2;
}

/** Whether two convex hulls of four points each are apart on one of the axes their points define, by a margin. */
static int hulls_apart(const GSHGEOM_SK_Vec * a, const GSHGEOM_SK_Vec * b, double margin) {
  const GSHGEOM_SK_Vec * sets[2] = {a, b};
  for (int s = 0; s < 2; s++) {
    for (int i = 0; i < 4; i++) {
      for (int j = i + 1; j < 4; j++) {
        GSHGEOM_SK_Vec e = gshgeom_sk_sub(sets[s][j], sets[s][i]);
        GSHGEOM_SK_Vec ax = gshgeom_sk_left(e);
        double len = gshgeom_sk_norm(ax);
        if (!(len > 0.0)) {
          continue;
        }
        ax = gshgeom_sk_scale(1.0 / len, ax);
        double amin = 1e300, amax = -1e300, bmin = 1e300, bmax = -1e300;
        for (int k = 0; k < 4; k++) {
          double pa = gshgeom_sk_dot(a[k], ax), pb = gshgeom_sk_dot(b[k], ax);
          amin = pa < amin ? pa : amin;
          amax = pa > amax ? pa : amax;
          bmin = pb < bmin ? pb : bmin;
          bmax = pb > bmax ? pb : bmax;
        }
        if (amax + margin < bmin || bmax + margin < amin) {
          return 1;
        }
      }
    }
  }
  return 0;
}

/*
 * The band of a piece as two offset cubics (curve mode).  Accepted when the piece's band is simple (band_analyze), the
 * fit error of SG38 is within the budget and a quarter of h, the two fitted rails run forward along both end tangents
 * at every control vector, and the hulls of the rails are apart.  Then the fitted loop (left rail, end normal, right
 * rail backward, start normal) is a simple closed curve: each rail is monotone along each end tangent, so it never
 * returns to an end normal's line, and the two rails, in disjoint hulls, never meet.  Its orientation is the true
 * band's: the fit moves the boundary by at most eps <= h / 4, and the band's signed area 2 h length(c) cannot change
 * sign under a move of eps along a boundary of length at most 4 length(c).  Returns 1 when the piece was added to
 * the run, 0 to try something else.
 */
static int band_curve(const Seg * sg, const Piece * p, const BandInfo * bip, GSHGEOM_Result * r) {
  GSHGEOM_SK_Ctx * c = sg->c;
  int n = sg->n;
  double h = c->h;
  BandInfo bi = *bip;
  double G[6], m;
  if (!bounds_g(p, n, &bi, G, &m)) {
    return 0;
  }
  double o4 = fourth_bound(G, m, h);
  double eps = o4 / 384.0 * UP;
  double bmax = max2(max2(bi.beta[0], bi.beta[bi.m - 1]), max2(p->bd0, p->bd1));
  double corner = c->norm * (bi.ep + h * bmax * 1.1) * UP;
  double d = eps * c->norm + corner;
  if (!(eps <= 0.25 * h) || !(bi.kappa * (h + eps) < bi.cosb) || !gshgeom_sk_fits(c, d, d)) {
    return 0;
  }
  /* the end data: the shared directions, and the curve's derivatives along them */
  GSHGEOM_SK_Vec t0 = gshgeom_sk_unit(p->dir0), t1 = gshgeom_sk_unit(p->dir1);
  GSHGEOM_SK_Vec c1s, c2s, c1e, c2e;
  end_derivs(p, n, 0, &c1s, &c2s);
  end_derivs(p, n, 1, &c1e, &c2e);
  double sp0 = gshgeom_sk_dot(c1s, t0), sp1 = gshgeom_sk_dot(c1e, t1);
  if (!(sp0 > 0.0) || !(sp1 > 0.0)) {
    return 0;
  }
  /* the directions the coefficients and the shared tangents give must agree, within their bounds */
  if (!(fabs_(gshgeom_sk_cross(t0, bi.u[0])) <= 4.0 * (bi.beta[0] + p->bd0) + 0x1p-30) ||
      !(fabs_(gshgeom_sk_cross(t1, bi.u[bi.m - 1])) <= 4.0 * (bi.beta[bi.m - 1] + p->bd1) + 0x1p-30)) {
    return 0;
  }
  GSHGEOM_SK_Vec tp0 = gshgeom_sk_scale(1.0 / sp0, gshgeom_sk_sub(c2s, gshgeom_sk_scale(gshgeom_sk_dot(c2s, t0), t0)));
  GSHGEOM_SK_Vec tp1 = gshgeom_sk_scale(1.0 / sp1, gshgeom_sk_sub(c2e, gshgeom_sk_scale(gshgeom_sk_dot(c2e, t1), t1)));
  GSHGEOM_SK_Vec P0 = bi.P[0], P3 = bi.P[n];
  GSHGEOM_SK_Vec rail[2][4]; /* [0] the left offset (s = +1), [1] the right (s = -1), build space */
  for (int side = 0; side < 2; side++) {
    double sgn = side == 0 ? 1.0 : -1.0;
    GSHGEOM_SK_Vec b0 = gshgeom_sk_add(P0, gshgeom_sk_scale(sgn * h, gshgeom_sk_left(t0)));
    GSHGEOM_SK_Vec b3 = gshgeom_sk_add(P3, gshgeom_sk_scale(sgn * h, gshgeom_sk_left(t1)));
    GSHGEOM_SK_Vec o0 = gshgeom_sk_add(c1s, gshgeom_sk_scale(sgn * h, gshgeom_sk_left(tp0)));
    GSHGEOM_SK_Vec o1 = gshgeom_sk_add(c1e, gshgeom_sk_scale(sgn * h, gshgeom_sk_left(tp1)));
    rail[side][0] = b0;
    rail[side][1] = gshgeom_sk_add(b0, gshgeom_sk_scale(1.0 / 3.0, o0));
    rail[side][2] = gshgeom_sk_sub(b3, gshgeom_sk_scale(1.0 / 3.0, o1));
    rail[side][3] = b3;
    /* the rail runs forward along both end tangents at every control vector: injective, never back to an end normal */
    for (int i = 0; i < 3; i++) {
      GSHGEOM_SK_Vec hv = gshgeom_sk_sub(rail[side][i + 1], rail[side][i]);
      double hl = gshgeom_sk_norm(hv);
      if (!(gshgeom_sk_dot(hv, t0) > 0x1p-20 * hl) || !(gshgeom_sk_dot(hv, t1) > 0x1p-20 * hl)) {
        return 0;
      }
    }
  }
  double scale = 0.0;
  for (int i = 0; i < 4; i++) {
    scale = max2(scale, max2(mag(rail[0][i].x, rail[0][i].y), mag(rail[1][i].x, rail[1][i].y)));
  }
  if (!hulls_apart(rail[0], rail[1], 0x1p-30 * scale + 4.0 * eps)) {
    return 0;
  }
  GSHGEOM_SK_Band b;
  for (int i = 0; i < 4; i++) {
    b.l[i] = gshgeom_sk_add(sg->anchor, gshgeom_sk_lin(c, rail[0][i]));
    b.r[i] = gshgeom_sk_add(sg->anchor, gshgeom_sk_lin(c, rail[1][i]));
  }
  *r = gshgeom_skc_add_band(c, &b);
  return 1;
}

/** A triangle of target points, written positively wound (decided exactly) and dropped when it has no area. */
static GSHGEOM_Result triangle(GSHGEOM_SK_Ctx * c, GSHGEOM_SK_Vec a, GSHGEOM_SK_Vec b, GSHGEOM_SK_Vec d) {
  GSHGEOM_Point pa = {a.x, a.y, 0.0}, pb = {b.x, b.y, 0.0}, pd = {d.x, d.y, 0.0};
  int s = 0;
  (void)gshgeom_orient2d(pa, pb, pd, &s);
  if (s == 0) {
    return GSHGEOM_OK;
  }
  GSHGEOM_SK_Vec ring[3] = {a, s > 0 ? b : d, s > 0 ? d : b};
  return gshgeom_skc_ring(c, ring, 3);
}

static int orient_of(GSHGEOM_SK_Vec a, GSHGEOM_SK_Vec b, GSHGEOM_SK_Vec d) {
  GSHGEOM_Point pa = {a.x, a.y, 0.0}, pb = {b.x, b.y, 0.0}, pd = {d.x, d.y, 0.0};
  int s = 0;
  (void)gshgeom_orient2d(pa, pb, pd, &s);
  return s;
}

/*
 * The piece of a fold as polygons (curve mode; SG37), where the offset folds and no band is simple.  The region the
 * normals of the piece sweep is the image of `F(u, lambda) = c(u) + lambda n(u)` over `[0, 1] x [-h, h]`; the bilinear
 * patch `B` through the four corners `L0 R0 L1 R1` (the ends of the two end normals) is the union of the segments
 * from `(1-s) L0 + s L1` to `(1-s) R0 + s R1`: the quadrilateral `L0 L1 R1 R0` when the end normals do not cross, and the
 * two lobes `L0 L1 X` and `R1 R0 X` (`X` where they cross) when they do.  For each `lambda`, `F(u, lambda)` is within
 * `sup |F_uu| / 8` of its chord, so the two images are within `(sup |c''| + h sup |n''|) / 8` of each other,
 * point for point: the bound of second_bound(), second order in the piece's extent however the offset folds.  The
 * piece is written as positively wound triangles, so the nonzero fill is their union.
 */
static int chord_curve(const Seg * sg, const Piece * p, const BandInfo * bip, GSHGEOM_Result * r) {
  GSHGEOM_SK_Ctx * c = sg->c;
  int n = sg->n;
  double h = c->h;
  BandInfo bi = *bip;
  double G[6], m;
  if (!bounds_g(p, n, &bi, G, &m)) {
    return 0;
  }
  double d2 = second_bound(G, m, h) / 8.0 * UP;
  double bmax = max2(max2(bi.beta[0], bi.beta[bi.m - 1]), max2(p->bd0, p->bd1));
  double corner = c->norm * (bi.ep + h * bmax * 1.1) * UP;
  double d = d2 * c->norm + corner;
  if (!gshgeom_sk_fits(c, d, d)) {
    return 0;
  }
  GSHGEOM_SK_Vec t0 = gshgeom_sk_unit(p->dir0), t1 = gshgeom_sk_unit(p->dir1);
  GSHGEOM_SK_Vec hn0 = gshgeom_sk_scale(h, gshgeom_sk_left(t0)), hn1 = gshgeom_sk_scale(h, gshgeom_sk_left(t1));
  GSHGEOM_SK_Vec P0 = bi.P[0], P1 = bi.P[n];
  GSHGEOM_SK_Vec L0 = gshgeom_sk_add(sg->anchor, gshgeom_sk_lin(c, gshgeom_sk_add(P0, hn0)));
  GSHGEOM_SK_Vec R0 = gshgeom_sk_add(sg->anchor, gshgeom_sk_lin(c, gshgeom_sk_sub(P0, hn0)));
  GSHGEOM_SK_Vec L1 = gshgeom_sk_add(sg->anchor, gshgeom_sk_lin(c, gshgeom_sk_add(P1, hn1)));
  GSHGEOM_SK_Vec R1 = gshgeom_sk_add(sg->anchor, gshgeom_sk_lin(c, gshgeom_sk_sub(P1, hn1)));
  if ((*r = gshgeom_skc_flush_run(c)) != GSHGEOM_OK) {
    return 1;
  }
  int s1 = orient_of(L0, R0, L1), s2 = orient_of(L0, R0, R1), s3 = orient_of(L1, R1, L0), s4 = orient_of(L1, R1, R0);
  if (s1 * s2 < 0 && s3 * s4 < 0) {
    /* the end normals cross at X: the lobes */
    GSHGEOM_SK_Vec e0 = gshgeom_sk_sub(R0, L0), e1 = gshgeom_sk_sub(R1, L1);
    double den = gshgeom_sk_cross(e0, e1);
    GSHGEOM_SK_Vec X = gshgeom_sk_add(L0, gshgeom_sk_scale(gshgeom_sk_cross(gshgeom_sk_sub(L1, L0), e1) / den, e0));
    *r = triangle(c, X, L0, L1);
    if (*r == GSHGEOM_OK) {
      *r = triangle(c, X, R1, R0);
    }
    return 1;
  }
  /* the quadrilateral, as the two triangles of the diagonal that lies inside it */
  int a1 = orient_of(L0, L1, R1), a2 = orient_of(L0, R1, R0);
  if (a1 * a2 > 0 || (a1 == 0 && a2 == 0)) {
    *r = triangle(c, L0, L1, R1);
    if (*r == GSHGEOM_OK) {
      *r = triangle(c, L0, R1, R0);
    }
  } else {
    *r = triangle(c, L1, R1, R0);
    if (*r == GSHGEOM_OK) {
      *r = triangle(c, L1, R0, L0);
    }
  }
  return 1;
}

/** Rotate a unit vector by at least `beta` (tan approximated from above) in the sense s (+1 positive). */
static GSHGEOM_SK_Vec turn_by(GSHGEOM_SK_Vec u, double beta, double s) {
  GSHGEOM_SK_Vec perp = s > 0 ? gshgeom_sk_right(u) : gshgeom_sk_left(u);
  return gshgeom_sk_unit(gshgeom_sk_add(u, gshgeom_sk_scale(1.1 * beta, perp)));
}

/**
 * Try the fan; returns 1 when emitted (or on an error, in *r), 0 to split.  In curve mode the fan is measured in build
 * space (the displacement of the piece's points from its start, the overshoot of the cone) and inflated there by that
 * displacement, so that it contains the piece instead of coming within the displacement of it: neighbouring fans
 * overlap and no crack is left between them.  The cost in target units is the build-space bound times the norm.
 */
static int fan(const Seg * sg, const Piece * p, const Know * k, GSHGEOM_Result * r) {
  GSHGEOM_SK_Ctx * c = sg->c;
  int n = sg->n;
  int cvm = c->cv != NULL;
  /* the ends' normals: the first and the last coefficient that is not exactly zero, which must be known */
  int i0 = -1, i1 = -1;
  for (int i = 0; i < k->m && i0 < 0; i++) {
    if (!k->zero[i]) {
      i0 = i;
    }
  }
  for (int i = k->m - 1; i >= 0 && i1 < 0; i--) {
    if (!k->zero[i]) {
      i1 = i;
    }
  }
  if (i0 < 0 || i1 < 0) {
    return 0;
  }
  double l0 = gshgeom_sk_norm(k->e[i0]), l1 = gshgeom_sk_norm(k->e[i1]);
  if (!(l0 > 4.0 * k->err[i0]) || !(l1 > 4.0 * k->err[i1])) {
    return 0;
  }
  GSHGEOM_SK_Vec t0 = gshgeom_sk_unit(k->e[i0]), t1 = gshgeom_sk_unit(k->e[i1]);
  double b0 = widen(k, i0), b1 = widen(k, i1);
  /* dm: the piece's control points' distance from its start (device pixels; build units in curve mode), with their error */
  GSHGEOM_SK_Vec P[4];
  double ep = 0.0, dm = 0.0;
  for (int i = 0; i <= n; i++) {
    P[i] = gshgeom_sk_v(p->nx[i] / p->d[i], p->ny[i] / p->d[i]);
    double e = (p->en[i] + mag(P[i].x, P[i].y) * p->ed[i]) / p->d[i] + 2.0 * U * mag(P[i].x, P[i].y);
    ep = max2(ep, e * 1.5);
  }
  double nrm = cvm ? 1.0 : c->norm;
  for (int i = 1; i <= n; i++) {
    GSHGEOM_SK_Vec v = gshgeom_sk_sub(P[i], P[0]);
    dm = max2(dm, gshgeom_sk_norm(cvm ? v : gshgeom_sk_lin(c, v)));
  }
  dm = (dm + 2.0 * nrm * ep) * UP;
  GSHGEOM_SK_Vec centre;
  centre = gshgeom_sk_add(sg->anchor, gshgeom_sk_lin(c, P[0]));
  double d_in = dm + nrm * ep;
  double hh = cvm ? c->h : c->h_dev;
  /* a double sector: every known direction within an arc below a half turn */
  if (!k->unknown) {
    GSHGEOM_SK_Vec u[6];
    double beta[6];
    int m = units(k, u, beta);
    for (int lo = 0; lo < m; lo++) {
      GSHGEOM_SK_Vec a = turn_by(u[lo], beta[lo], -1.0);
      for (int hi = 0; hi < m; hi++) {
        GSHGEOM_SK_Vec b = turn_by(u[hi], beta[hi], 1.0);
        if (!(gshgeom_sk_cross(a, b) > 0.0)) {
          continue; /* not an arc below a half turn from a to b */
        }
        /* every other direction, widened, inside: on the positive side of a and the negative side of b by at least
         * sin(beta), which 1.1 beta exceeds; a and b were turned past lo and hi by at least their beta */
        int inside = 1;
        for (int i = 0; i < m && inside; i++) {
          if (i == lo || i == hi) {
            continue;
          }
          double s = 1.1 * beta[i];
          if (!(gshgeom_sk_cross(a, u[i]) >= s && gshgeom_sk_cross(u[i], b) >= s)) {
            inside = 0;
          }
        }
        if (!inside) {
          continue;
        }
        /* the covered arc runs between t0 and t1 (both inside [a, b]); the overshoot on each side */
        int fwd = gshgeom_sk_cross(t0, t1) >= 0.0;
        GSHGEOM_SK_Vec first = fwd ? t0 : t1, last = fwd ? t1 : t0;
        double bf = fwd ? b0 : b1, bl = fwd ? b1 : b0;
        double over = max2(gshgeom_sk_norm(gshgeom_sk_sub(a, first)) + 1.1 * bf,
            gshgeom_sk_norm(gshgeom_sk_sub(last, b)) + 1.1 * bl);
        double extra = (dm + hh * over) * UP;
        if (cvm ? !(c->norm * (extra + 1.4143 * d_in) * UP <= c->cv->budget) :
                  !(gshgeom_sk_room(c, d_in) - extra > 0.0)) {
          return 0;
        }
        GSHGEOM_SK_Vec na = gshgeom_sk_left(a), nb = gshgeom_sk_left(b);
        if (cvm) {
          if ((*r = gshgeom_skc_flush_run(c)) != GSHGEOM_OK) {
            return 1;
          }
          *r = gshgeom_skc_fan_sector(c, centre, na, nb, d_in * UP);
          if (*r == GSHGEOM_OK) {
            *r = gshgeom_skc_fan_sector(c, centre, gshgeom_sk_scale(-1.0, na), gshgeom_sk_scale(-1.0, nb), d_in * UP);
          }
          return 1;
        }
        *r = gshgeom_sk_sector(c, centre, na, nb, d_in, extra);
        if (*r == GSHGEOM_OK) {
          *r = gshgeom_sk_sector(c, centre, gshgeom_sk_scale(-1.0, na), gshgeom_sk_scale(-1.0, nb), d_in, extra);
        }
        return 1;
      }
    }
  }
  /* a disc: the ends' tangents nearly opposite */
  if (gshgeom_sk_dot(t0, t1) < 0.0) {
    double deficit = gshgeom_sk_norm(gshgeom_sk_add(t0, t1)) + 1.1 * (b0 + b1);
    double extra = (dm + hh * deficit) * UP;
    if (cvm ? (c->norm * (extra + d_in) * UP <= c->cv->budget) : (gshgeom_sk_room(c, d_in) - extra > 0.0)) {
      if (cvm) {
        if ((*r = gshgeom_skc_flush_run(c)) != GSHGEOM_OK) {
          return 1;
        }
        *r = gshgeom_skc_disc(c, centre, c->h + d_in * UP);
        return 1;
      }
      *r = gshgeom_sk_disc(c, centre, d_in, extra);
      return 1;
    }
  }
  return 0;
}

static GSHGEOM_Result piece(const Seg * sg, const Piece * p, int depth) {
  /* a bound on the work: a piece is either emitted (at least three points of max_points) or split, and a curve that
   * needs more pieces than the part can hold is refused before it is searched without end */
  if (++sg->c->work > sg->c->work_cap) {
    return GSHGEOM_ERR_LIMIT;
  }
  Know k;
  know(p, sg->n, &k);
  GSHGEOM_Result r = GSHGEOM_OK;
  if (sg->c->cv) {
    BandInfo bi;
    int st = band_analyze(sg, p, &k, &bi);
    if (st == 1 && p->has0 && p->has1 && band_curve(sg, p, &bi, &r)) {
      return r;
    }
    if (st == 2 && p->has0 && p->has1 && chord_curve(sg, p, &bi, &r)) {
      return r;
    }
  } else if (band(sg, p, &k, &r)) {
    return r;
  }
  if (fan(sg, p, &k, &r)) {
    return r;
  }
  if (depth >= sg->c->depth_cap) {
    return GSHGEOM_ERR_LIMIT; /* a forced piece is refused, never drawn coarser */
  }
  Piece l, rt;
  split(p, sg->n, &l, &rt);
  r = piece(sg, &l, depth + 1);
  if (r == GSHGEOM_OK) {
    r = piece(sg, &rt, depth + 1);
  }
  return r;
}

GSHGEOM_Result gshgeom_sk_body(GSHGEOM_SK_Ctx * c, const GSHGEOM_SK_Seg * s) {
  GSHGEOM_SK_Vec o = s->p[0];
  Seg sg;
  sg.c = c;
  sg.n = s->n;
  GSHGEOM_Result r = gshgeom_sk_anchor(c, o, &sg.anchor);
  if (r != GSHGEOM_OK) {
    return r;
  }
  if (s->n == 1) {
    GSHGEOM_SK_Vec a = gshgeom_sk_sub(s->p[1], o);
    GSHGEOM_SK_Vec hn = gshgeom_sk_scale(c->h, gshgeom_sk_left(gshgeom_sk_unit(a)));
    GSHGEOM_SK_Vec q[4] = {hn, gshgeom_sk_add(a, hn), gshgeom_sk_sub(a, hn), gshgeom_sk_scale(-1.0, hn)};
    GSHGEOM_SK_Vec piv[2] = {gshgeom_sk_v(0.0, 0.0), a};
    return gshgeom_sk_emit_conform(c, sg.anchor, q, 4, 0.0, 0.0, piv, 2);
  }
  Piece p;
  memset(&p, 0, sizeof p);
  p.z0 = s->p[1].x == s->p[0].x && s->p[1].y == s->p[0].y && s->w[1] == s->w[0];
  p.zn = s->p[s->n].x == s->p[s->n - 1].x && s->p[s->n].y == s->p[s->n - 1].y && s->w[s->n] == s->w[s->n - 1];
  for (int i = 0; i <= s->n; i++) {
    double dx = s->p[i].x - o.x, dy = s->p[i].y - o.y;
    p.nx[i] = s->w[i] * dx;
    p.ny[i] = s->w[i] * dy;
    p.d[i] = s->w[i];
    /* the difference and the product, each rounded once */
    p.en[i] = 2.0 * U * mag(p.nx[i], p.ny[i]) * UP;
    p.ed[i] = 0.0;
  }
  if (c->cv) {
    /* the segment's own end directions (SG5): raw differences of the control points, shared with the joins and caps */
    p.has0 = gshgeom_sk_start_dir(s, &p.dir0) >= 0;
    p.has1 = gshgeom_sk_end_dir(s, &p.dir1) >= 0;
    p.bd0 = p.bd1 = 0x1p-40;
  }
  GSHGEOM_Result res = piece(&sg, &p, 0);
  if (c->cv) {
    GSHGEOM_Result fr = gshgeom_skc_flush_run(c);
    res = res != GSHGEOM_OK ? res : fr;
  }
  return res;
}
