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

/** Try the band; returns 1 when emitted (or on an error, in *r), 0 to try something else. */
static int band(const Seg * sg, const Piece * p, const Know * k, GSHGEOM_Result * r) {
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
  GSHGEOM_SK_Vec u[6], ud[6];
  double beta[6], betad[6];
  int m = units(k, u, beta);
  double cosb, sinb, cosd, sind;
  if (!quarter_cone(u, beta, m, &cosb, &sinb)) {
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
  if (!quarter_cone(ud, betad, m, &cosd, &sind)) {
    return 0;
  }
  /* the least projection on u[0], and the curvature bound */
  double minproj = -1.0, maxe = 0.0, maxde = 0.0, maxd = 0.0;
  for (int i = 0; i < k->m; i++) {
    double pr = gshgeom_sk_dot(k->e[i], u[0]) - k->err[i];
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
  if (!(kappa * c->h < cosb)) {
    return 0; /* the normals may meet within h: not a band */
  }
  /* the projected control points (build space, relative to the segment's start) and their error */
  GSHGEOM_SK_Vec P[4];
  double ep = 0.0;
  for (int i = 0; i <= n; i++) {
    P[i] = gshgeom_sk_v(p->nx[i] / p->d[i], p->ny[i] / p->d[i]);
    double e = (p->en[i] + mag(P[i].x, P[i].y) * p->ed[i]) / p->d[i] + 2.0 * U * mag(P[i].x, P[i].y);
    ep = max2(ep, e * 1.5);
  }
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
  if (!gshgeom_sk_fits(c, d, d)) {
    return 0;
  }
  *r = gshgeom_sk_emit(c, sg->anchor, q, 4, d, d);
  return 1;
}

/** Rotate a unit vector by at least `beta` (tan approximated from above) in the sense s (+1 positive). */
static GSHGEOM_SK_Vec turn_by(GSHGEOM_SK_Vec u, double beta, double s) {
  GSHGEOM_SK_Vec perp = s > 0 ? gshgeom_sk_right(u) : gshgeom_sk_left(u);
  return gshgeom_sk_unit(gshgeom_sk_add(u, gshgeom_sk_scale(1.1 * beta, perp)));
}

/** Try the fan; returns 1 when emitted (or on an error, in *r), 0 to split. */
static int fan(const Seg * sg, const Piece * p, const Know * k, GSHGEOM_Result * r) {
  GSHGEOM_SK_Ctx * c = sg->c;
  int n = sg->n;
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
  /* dm: the piece's control points' device distance from its start, with their error */
  GSHGEOM_SK_Vec P[4];
  double ep = 0.0, dm = 0.0;
  for (int i = 0; i <= n; i++) {
    P[i] = gshgeom_sk_v(p->nx[i] / p->d[i], p->ny[i] / p->d[i]);
    double e = (p->en[i] + mag(P[i].x, P[i].y) * p->ed[i]) / p->d[i] + 2.0 * U * mag(P[i].x, P[i].y);
    ep = max2(ep, e * 1.5);
  }
  for (int i = 1; i <= n; i++) {
    dm = max2(dm, gshgeom_sk_norm(gshgeom_sk_lin(c, gshgeom_sk_sub(P[i], P[0]))));
  }
  dm = (dm + 2.0 * c->norm * ep) * UP;
  GSHGEOM_SK_Vec centre;
  centre = gshgeom_sk_add(sg->anchor, gshgeom_sk_lin(c, P[0]));
  double d_in = dm + c->norm * ep;
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
        double extra = (dm + c->h_dev * over) * UP;
        if (!(gshgeom_sk_room(c, d_in) - extra > 0.0)) {
          return 0;
        }
        GSHGEOM_SK_Vec na = gshgeom_sk_left(a), nb = gshgeom_sk_left(b);
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
    double extra = (dm + c->h_dev * deficit) * UP;
    if (gshgeom_sk_room(c, d_in) - extra > 0.0) {
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
  if (band(sg, p, &k, &r) || fan(sg, p, &k, &r)) {
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
    return gshgeom_sk_emit(c, sg.anchor, q, 4, 0.0, 0.0);
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
  return piece(&sg, &p, 0);
}
