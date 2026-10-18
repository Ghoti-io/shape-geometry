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
 * Along-path parametrization (along.h): the measured path, its table of
 * quadrature panels, and the queries by distance.
 *
 * Everything here is `+ - * /`, `gmath_hypot` and a fixed order of operations,
 * so the table and every answer are the same bits on every target.  Nothing
 * comes from libm and no operation is fused (-ffp-contract=off).
 *
 * A curve is held in homogeneous form: numerator coefficients `N_i = w_i P_i`
 * and denominator coefficients `w_i` for a rational kind, the control points
 * and no denominator for a polynomial one (so that a polynomial curve never
 * multiplies by a `1.0` that de Casteljau's recurrence would have rounded).  The
 * velocity of a rational curve is `(N'D - ND') / D^2`, and `E = N'D - ND'`
 * is its direction; `E` is a polynomial, so its derivatives at a parameter,
 * which the tangent at a cusp needs, come from the same recurrences.
 */

#include <ghoti.io/shape-geometry/macros.h>

#include <ghoti.io/math/elementary.h>
#include <ghoti.io/shape-geometry/along.h>
#include <ghoti.io/shape-geometry/hash.h>

#include "../core/grow.h"

#include <string.h>

/* ---------------------------------------------------------------------- */
/* The Gauss-Legendre rule                                                 */
/* ---------------------------------------------------------------------- */

/* Eight points, as four symmetric pairs, the smallest node first.  Generated
 * by tools/gen-along-vectors.py (mpmath at 400 bits, rounded to nearest); a
 * test integrates every monomial through degree 15 exactly. */
#define GL_HALF 4
static const double gl_node[GL_HALF] = {
    0x1.77ac94f3c7345p-3, 0x1.0d129583284b4p-1,
    0x1.97e4ab249f41ep-1, 0x1.ebab1cb0acc67p-1};
static const double gl_weight[GL_HALF] = {
    0x1.736360b199343p-2, 0x1.413c50a255615p-2,
    0x1.c76fb531d2b96p-3, 0x1.9ea1d04ca0374p-4};

/* ---------------------------------------------------------------------- */
/* Curves in homogeneous form                                              */
/* ---------------------------------------------------------------------- */

typedef struct {
  int n;          /* degree, 1 to 3 */
  int dim;        /* 2 or 3 */
  int rational;
  double c[3][4]; /* c[axis][i]: the numerator coefficients */
  double d[4];    /* the denominator coefficients (rational only) */
} Hom;

/* A finite double: v - v is zero for every finite v. */
static int is_finite(double v) {
  return (v - v) == 0.0;
}

static double dabs(double v) {
  return v < 0.0 ? -v : v;
}

static double axis_of(const GSHGEOM_Point * p, int axis) {
  return axis == 0 ? p->x : (axis == 1 ? p->y : p->z);
}

static void load(const GSHGEOM_Curve * curve, Hom * h) {
  memset(h, 0, sizeof *h);
  h->n = gshgeom_curve_point_count(curve->kind) - 1;
  h->dim = curve->dim;
  h->rational = gshgeom_curve_is_rational(curve->kind);
  for (int i = 0; i <= h->n; i++) {
    for (int a = 0; a < h->dim; a++) {
      double v = axis_of(&curve->p[i], a);
      h->c[a][i] = h->rational ? curve->w[i] * v : v;
    }
    h->d[i] = h->rational ? curve->w[i] : 1.0;
  }
}

/* The length of a vector in 2 or 3 dimensions. */
static double mag(int dim, const double * v) {
  double h = gmath_hypot(v[0], v[1]);
  return dim == 3 ? gmath_hypot(h, v[2]) : h;
}

/* de Casteljau's algorithm for a Bezier polynomial of degree m (m + 1
 * coefficients) at t.  t = 0 and t = 1 give the end coefficients exactly. */
static double bez(const double * c, int m, double t) {
  double h[4];
  for (int i = 0; i <= m; i++) {
    h[i] = c[i];
  }
  double s = 1.0 - t;
  for (int level = m; level > 0; level--) {
    for (int i = 0; i < level; i++) {
      h[i] = s * h[i] + t * h[i + 1];
    }
  }
  return h[0];
}

/* The k-th derivative of a Bezier polynomial of degree m at t: the k-th
 * difference of the coefficients, times m!/(m-k)!, as a Bezier polynomial of
 * degree m - k. */
static double deriv(const double * c, int m, int k, double t) {
  if (k > m) {
    return 0.0;
  }
  double h[4];
  for (int i = 0; i <= m; i++) {
    h[i] = c[i];
  }
  for (int j = 0; j < k; j++) {
    for (int i = 0; i < m - j; i++) {
      h[i] = h[i + 1] - h[i];
    }
  }
  double factor = 1.0;
  for (int j = 0; j < k; j++) {
    factor *= (double)(m - j);
  }
  return factor * bez(h, m - k, t);
}

static const int binom[5][5] = {
    {1, 0, 0, 0, 0}, {1, 1, 0, 0, 0}, {1, 2, 1, 0, 0}, {1, 3, 3, 1, 0},
    {1, 4, 6, 4, 1}};

/* The k-th derivative of E = N'D - ND' at t, by Leibniz's rule:
 * E^(k) = sum over j of C(k, j) (N^(j+1) D^(k-j) - N^(j) D^(k-j+1)).  For a
 * polynomial curve (D = 1) it is N^(k+1). */
static void e_derivative(const Hom * h, int k, double t, double e[3]) {
  for (int a = 0; a < h->dim; a++) {
    if (!h->rational) {
      e[a] = deriv(h->c[a], h->n, k + 1, t);
      continue;
    }
    double sum = 0.0;
    for (int j = 0; j <= k; j++) {
      double term = deriv(h->c[a], h->n, j + 1, t) * deriv(h->d, h->n, k - j, t) -
          deriv(h->c[a], h->n, j, t) * deriv(h->d, h->n, k - j + 1, t);
      sum += (double)binom[k][j] * term;
    }
    e[a] = sum;
  }
  if (h->dim == 2) {
    e[2] = 0.0;
  }
}

/* The speed |B'(t)|. */
static double speed(const Hom * h, double t) {
  double e[3];
  e_derivative(h, 0, t, e);
  double m = mag(h->dim, e);
  if (!h->rational) {
    return m;
  }
  double d = bez(h->d, h->n, t);
  return m / (d * d);
}

/* The integral of the speed over [a, b] by the eight-point rule. */
static double gl(const Hom * h, double a, double b) {
  double c = (a + b) * 0.5;
  double half = (b - a) * 0.5;
  double sum = 0.0;
  for (int i = 0; i < GL_HALF; i++) {
    double off = half * gl_node[i];
    sum += gl_weight[i] * (speed(h, c - off) + speed(h, c + off));
  }
  return half * sum;
}

/* The length of the control polygon of a curve: the scale of its tolerance.
 * Zero exactly when the curve has no extent. */
static double polygon_length(const GSHGEOM_Curve * curve) {
  int n = gshgeom_curve_point_count(curve->kind) - 1;
  double total = 0.0;
  for (int i = 0; i < n; i++) {
    double d[3] = {curve->p[i + 1].x - curve->p[i].x,
        curve->p[i + 1].y - curve->p[i].y, curve->p[i + 1].z - curve->p[i].z};
    total += mag(curve->dim, d);
  }
  return total;
}

/* ---------------------------------------------------------------------- */
/* The table                                                               */
/* ---------------------------------------------------------------------- */

typedef struct {
  size_t seg;
  double t0, t1;
  double s0;
  double len;
} Row;

/* Rows being made.  With `store` clear they are only counted and summed (the
 * single-curve length). */
typedef struct {
  GSHGEOM_Allocator alloc;
  int store;
  Row * row;
  size_t n, cap, max;
  double sum; /* the cumulative length so far */
  size_t cur; /* the segment the rows belong to */
} Table;

static GSHGEOM_Result emit(Table * tab, double t0, double t1, double len) {
  if (tab->n >= tab->max) {
    return GSHGEOM_ERR_LIMIT;
  }
  if (tab->store) {
    void * grown = tab->row;
    GSHGEOM_Result r = gshgeom_grow(
        &tab->alloc, &grown, &tab->cap, tab->n + 1, sizeof *tab->row);
    if (r != GSHGEOM_OK) {
      return r;
    }
    tab->row = grown;
    Row * row = &tab->row[tab->n];
    row->seg = tab->cur;
    row->t0 = t0;
    row->t1 = t1;
    row->s0 = tab->sum;
    row->len = len;
  }
  tab->n++;
  tab->sum += len;
  return GSHGEOM_OK;
}

/* Accept [a, b] when its value differs from its two halves by at most
 * `scaled * (b - a)`, where `scaled` is the tolerance times the polygon length;
 * the halves are then the rows.  Otherwise refine each half. */
static GSHGEOM_Result refine(const Hom * h, Table * tab, double a, double b,
    double whole, double scaled, int depth) {
  double m = (a + b) * 0.5;
  double left = gl(h, a, m);
  double right = gl(h, m, b);
  double diff = dabs(whole - (left + right));
  if (depth >= GSHGEOM_MEASURE_MAX_DEPTH || diff <= scaled * (b - a)) {
    GSHGEOM_Result r = emit(tab, a, m, left);
    return r != GSHGEOM_OK ? r : emit(tab, m, b, right);
  }
  GSHGEOM_Result r = refine(h, tab, a, m, left, scaled, depth + 1);
  return r != GSHGEOM_OK ? r : refine(h, tab, m, b, right, scaled, depth + 1);
}

/* Measure one curve that has extent (its polygon length @p poly is not zero)
 * into @p tab. */
static GSHGEOM_Result measure(const GSHGEOM_Curve * curve, double poly,
    double tolerance, Table * tab) {
  if (curve->kind == GSHGEOM_CURVE_LINE) {
    return emit(tab, 0.0, 1.0, poly);
  }
  Hom h;
  load(curve, &h);
  double whole = gl(&h, 0.0, 1.0);
  return refine(&h, tab, 0.0, 1.0, whole, tolerance * poly, 0);
}

static GSHGEOM_Result check_tolerance(double tolerance, double * out) {
  if (!is_finite(tolerance) || tolerance < 0.0) {
    return GSHGEOM_ERR_INVALID;
  }
  if (tolerance == 0.0) {
    *out = GSHGEOM_MEASURE_TOLERANCE_DEFAULT;
    return GSHGEOM_OK;
  }
  if (tolerance < GSHGEOM_MEASURE_TOLERANCE_MIN ||
      tolerance > GSHGEOM_MEASURE_TOLERANCE_MAX) {
    return GSHGEOM_ERR_LIMIT;
  }
  *out = tolerance;
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_measured_curve_length(
    const GSHGEOM_Curve * curve, double tolerance, double * out) {
  if (!curve || !out) {
    return GSHGEOM_ERR_INVALID;
  }
  GSHGEOM_Result r = gshgeom_curve_check(curve);
  if (r != GSHGEOM_OK) {
    return r;
  }
  double tol = 0.0;
  r = check_tolerance(tolerance, &tol);
  if (r != GSHGEOM_OK) {
    return r;
  }
  double poly = polygon_length(curve);
  if (poly == 0.0) {
    *out = 0.0;
    return GSHGEOM_OK;
  }
  Table tab;
  memset(&tab, 0, sizeof tab);
  tab.max = GSHGEOM_LIMITS_DEFAULT_MAX_POINTS;
  r = measure(curve, poly, tol, &tab);
  if (r != GSHGEOM_OK) {
    return r;
  }
  *out = tab.sum;
  return GSHGEOM_OK;
}

/* ---------------------------------------------------------------------- */
/* The measured path                                                       */
/* ---------------------------------------------------------------------- */

typedef struct {
  GSHGEOM_Curve curve;
  size_t element;
} Seg;

struct GSHGEOM_Measured {
  GSHGEOM_Allocator alloc;
  int dim;
  int has_origin;
  GSHGEOM_Point origin;
  Seg * seg;
  size_t seg_n, seg_cap;
  Table tab;
  GSHGEOM_Measured_Subpath * sub;
  size_t sub_n, sub_cap;
};

void gshgeom_measured_destroy(GSHGEOM_Measured * m) {
  if (!m) {
    return;
  }
  GSHGEOM_Allocator alloc = m->alloc;
  gcu_allocator_free(&alloc, m->seg);
  gcu_allocator_free(&alloc, m->tab.row);
  gcu_allocator_free(&alloc, m->sub);
  gcu_allocator_free(&alloc, m);
}

static GSHGEOM_Result open_subpath(GSHGEOM_Measured * m) {
  void * grown = m->sub;
  GSHGEOM_Result r = gshgeom_grow(
      &m->alloc, &grown, &m->sub_cap, m->sub_n + 1, sizeof *m->sub);
  if (r != GSHGEOM_OK) {
    return r;
  }
  m->sub = grown;
  m->sub[m->sub_n].begin = m->tab.sum;
  m->sub[m->sub_n].end = m->tab.sum;
  m->sub[m->sub_n].closed = 0;
  m->sub_n++;
  return GSHGEOM_OK;
}

static GSHGEOM_Result add_segment(GSHGEOM_Measured * m, size_t element,
    const GSHGEOM_Curve * curve, double tolerance) {
  double poly = polygon_length(curve);
  if (poly == 0.0) {
    return GSHGEOM_OK; /* no extent: nothing selects it */
  }
  void * grown = m->seg;
  GSHGEOM_Result r = gshgeom_grow(
      &m->alloc, &grown, &m->seg_cap, m->seg_n + 1, sizeof *m->seg);
  if (r != GSHGEOM_OK) {
    return r;
  }
  m->seg = grown;
  m->seg[m->seg_n].curve = *curve;
  m->seg[m->seg_n].element = element;
  m->tab.cur = m->seg_n;
  m->seg_n++;
  return measure(curve, poly, tolerance, &m->tab);
}

GSHGEOM_Result gshgeom_measured_create(const GSHGEOM_Allocator * alloc,
    const GSHGEOM_Limits * limits, const GSHGEOM_Path * path, double tolerance,
    GSHGEOM_Measured ** out) {
  if (!path || !out) {
    return GSHGEOM_ERR_INVALID;
  }
  double tol = 0.0;
  GSHGEOM_Result r = check_tolerance(tolerance, &tol);
  if (r != GSHGEOM_OK) {
    return r;
  }
  if (!alloc) {
    alloc = gshgeom_allocator_default();
  }
  GSHGEOM_Measured * m = gcu_allocator_calloc(alloc, 1, sizeof *m);
  if (!m) {
    return GSHGEOM_ERR_OOM;
  }
  m->alloc = *alloc;
  m->dim = gshgeom_path_dimension(path);
  m->tab.alloc = *alloc;
  m->tab.store = 1;
  m->tab.max = (limits && limits->max_points) ? limits->max_points
                                               : GSHGEOM_LIMITS_DEFAULT_MAX_POINTS;
  int need_sub = 1;
  size_t count = gshgeom_path_element_count(path);
  for (size_t i = 0; i < count && r == GSHGEOM_OK; i++) {
    GSHGEOM_Element el;
    r = gshgeom_path_element(path, i, &el);
    if (r != GSHGEOM_OK) {
      break;
    }
    if (el.verb == GSHGEOM_VERB_MOVE) {
      if (!m->has_origin) {
        m->origin = el.curve.p[0];
        m->has_origin = 1;
      }
      need_sub = 1;
      continue;
    }
    if (need_sub) {
      r = open_subpath(m);
      need_sub = 0;
      if (r != GSHGEOM_OK) {
        break;
      }
    }
    r = add_segment(m, i, &el.curve, tol);
    if (r != GSHGEOM_OK) {
      break;
    }
    m->sub[m->sub_n - 1].end = m->tab.sum;
    if (el.verb == GSHGEOM_VERB_CLOSE) {
      m->sub[m->sub_n - 1].closed = 1;
      need_sub = 1;
    }
  }
  if (r != GSHGEOM_OK) {
    gshgeom_measured_destroy(m);
    return r;
  }
  *out = m;
  return GSHGEOM_OK;
}

int gshgeom_measured_dimension(const GSHGEOM_Measured * m) {
  return m ? m->dim : 0;
}

double gshgeom_measured_length(const GSHGEOM_Measured * m) {
  return m ? m->tab.sum : 0.0;
}

size_t gshgeom_measured_segment_count(const GSHGEOM_Measured * m) {
  return m ? m->seg_n : 0;
}

size_t gshgeom_measured_row_count(const GSHGEOM_Measured * m) {
  return m ? m->tab.n : 0;
}

size_t gshgeom_measured_subpath_count(const GSHGEOM_Measured * m) {
  return m ? m->sub_n : 0;
}

GSHGEOM_Result gshgeom_measured_subpath(const GSHGEOM_Measured * m, size_t index,
    GSHGEOM_Measured_Subpath * out) {
  if (!m || !out || index >= m->sub_n) {
    return GSHGEOM_ERR_INVALID;
  }
  *out = m->sub[index];
  return GSHGEOM_OK;
}

/* ---------------------------------------------------------------------- */
/* Queries                                                                 */
/* ---------------------------------------------------------------------- */

/* Solve the integral of the speed over [row.t0, u] = target for u, for a
 * target strictly inside (0, row.len): Newton's method inside a bracket, the
 * midpoint whenever the step leaves it or the speed vanishes, to 2^-52 of the
 * row's length or 64 steps. */
static double solve(const Hom * h, const Row * row, double target) {
  double lo = row->t0;
  double hi = row->t1;
  double u = lo + (hi - lo) * (target / row->len);
  if (!(u > lo && u < hi)) {
    u = lo + (hi - lo) * 0.5;
  }
  double stop = row->len * 0x1p-52;
  for (int iteration = 0; iteration < 64; iteration++) {
    double f = gl(h, row->t0, u) - target;
    if (dabs(f) <= stop) {
      return u;
    }
    if (f > 0.0) {
      hi = u;
    } else {
      lo = u;
    }
    double next = u - f / speed(h, u);
    if (!(next > lo && next < hi)) {
      next = lo + (hi - lo) * 0.5;
    }
    if (next == u) {
      return u;
    }
    u = next;
  }
  return u;
}

/* The largest term of E in a segment, for the test of "vanishes". */
static double e_scale(const Hom * h) {
  double dc = 0.0, cc = 0.0, dd = 0.0, dmax = 0.0;
  for (int i = 0; i <= h->n; i++) {
    for (int a = 0; a < h->dim; a++) {
      double v = dabs(h->c[a][i]);
      cc = v > cc ? v : cc;
      if (i < h->n) {
        double d = dabs(h->c[a][i + 1] - h->c[a][i]);
        dc = d > dc ? d : dc;
      }
    }
    dmax = h->d[i] > dmax ? h->d[i] : dmax;
    if (i < h->n) {
      double d = dabs(h->d[i + 1] - h->d[i]);
      dd = d > dd ? d : dd;
    }
  }
  if (!h->rational) {
    return dc;
  }
  return dc * dmax + cc * dd;
}

/* The unit tangent at u; 0 if the segment has no direction there. */
static int tangent(const Hom * h, double u, double out[3]) {
  double thr = e_scale(h) * 0x1p-48;
  for (int k = 0; k <= 4; k++) {
    double e[3];
    e_derivative(h, k, u, e);
    double m = mag(h->dim, e);
    if (m > thr) {
      /* The right-hand limit has the sign of the first nonzero derivative; at
       * u = 1, where there is no right-hand side, the left-hand one has the
       * sign of (-1)^k. */
      double sign = (u >= 1.0 && (k & 1)) ? -1.0 : 1.0;
      for (int a = 0; a < 3; a++) {
        out[a] = (sign * e[a]) / m;
      }
      return 1;
    }
  }
  return 0;
}

static GSHGEOM_Result query(const GSHGEOM_Measured * m, double s, int clamp,
    GSHGEOM_Measured_Frame * out) {
  if (!m || !out || !is_finite(s) || !m->has_origin) {
    return GSHGEOM_ERR_INVALID;
  }
  double length = m->tab.sum;
  if (s < 0.0 || s > length) {
    if (!clamp) {
      return GSHGEOM_ERR_INVALID;
    }
    s = s < 0.0 ? 0.0 : length;
  }
  GSHGEOM_Measured_Frame f;
  memset(&f, 0, sizeof f);
  if (m->tab.n == 0) {
    f.point = m->origin;
    *out = f;
    return GSHGEOM_OK;
  }
  /* The last row that starts at or before s: a boundary belongs to the later
   * segment. */
  size_t lo = 0;
  size_t hi = m->tab.n;
  while (lo < hi) {
    size_t mid = lo + (hi - lo) / 2;
    if (m->tab.row[mid].s0 <= s) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }
  const Row * row = &m->tab.row[lo - 1];
  const Seg * seg = &m->seg[row->seg];
  Hom h;
  load(&seg->curve, &h);
  double target = s - row->s0;
  double u;
  if (target <= 0.0) {
    u = row->t0;
  } else if (s >= row->s0 + row->len) {
    /* The end of the row, by the very sum that made the next row's start (or
     * the length): s - s0 can round to less than len, and then the end of the
     * path would not be exactly the end of its last curve. */
    u = row->t1;
  } else if (seg->curve.kind == GSHGEOM_CURVE_LINE) {
    u = target / row->len;
  } else {
    u = solve(&h, row, target);
  }
  if (gshgeom_curve_eval(&seg->curve, u, &f.point) != GSHGEOM_OK) {
    return GSHGEOM_ERR_INTERNAL;
  }
  double tan[3];
  if (tangent(&h, u, tan)) {
    f.has_tangent = 1;
    f.tangent.x = tan[0];
    f.tangent.y = tan[1];
    f.tangent.z = h.dim == 3 ? tan[2] : 0.0;
    if (h.dim == 2) {
      f.normal = gshgeom_normal_left(f.tangent);
      f.has_normal = 1;
    }
  }
  f.element = seg->element;
  f.t = u;
  *out = f;
  return GSHGEOM_OK;
}

GSHGEOM_Result gshgeom_measured_at(
    const GSHGEOM_Measured * m, double s, GSHGEOM_Measured_Frame * out) {
  return query(m, s, 0, out);
}

GSHGEOM_Result gshgeom_measured_at_clamped(
    const GSHGEOM_Measured * m, double s, GSHGEOM_Measured_Frame * out) {
  return query(m, s, 1, out);
}

/* ---------------------------------------------------------------------- */
/* Hash and dump                                                           */
/* ---------------------------------------------------------------------- */

/** The format version of the hash and the dump. */
#define MEASURED_FORMAT_VERSION 1u
/** The kind field: "MEAS" as little-endian bytes. */
#define MEASURED_KIND 0x5341454du

static void hash_point(GSHGEOM_Hash * h, const GSHGEOM_Point * p, int dim) {
  gshgeom_hash_f64(h, p->x);
  gshgeom_hash_f64(h, p->y);
  if (dim == 3) {
    gshgeom_hash_f64(h, p->z);
  }
}

uint64_t gshgeom_measured_hash(const GSHGEOM_Measured * m) {
  GSHGEOM_Hash h;
  gshgeom_hash_init(&h);
  if (!m) {
    return gshgeom_hash_value(&h);
  }
  gshgeom_hash_u32(&h, MEASURED_FORMAT_VERSION);
  gshgeom_hash_u32(&h, MEASURED_KIND);
  gshgeom_hash_u32(&h, (uint32_t)m->dim);
  gshgeom_hash_u64(&h, (uint64_t)m->seg_n);
  gshgeom_hash_u64(&h, (uint64_t)m->tab.n);
  gshgeom_hash_u64(&h, (uint64_t)m->sub_n);
  hash_point(&h, &m->origin, m->dim);
  gshgeom_hash_f64(&h, m->tab.sum);
  for (size_t i = 0; i < m->seg_n; i++) {
    const Seg * s = &m->seg[i];
    int n = gshgeom_curve_point_count(s->curve.kind);
    gshgeom_hash_u64(&h, (uint64_t)s->element);
    gshgeom_hash_u8(&h, (uint8_t)s->curve.kind);
    for (int k = 0; k < n; k++) {
      hash_point(&h, &s->curve.p[k], m->dim);
    }
    if (gshgeom_curve_is_rational(s->curve.kind)) {
      for (int k = 0; k < n; k++) {
        gshgeom_hash_f64(&h, s->curve.w[k]);
      }
    }
  }
  for (size_t i = 0; i < m->tab.n; i++) {
    const Row * r = &m->tab.row[i];
    gshgeom_hash_u64(&h, (uint64_t)r->seg);
    gshgeom_hash_f64(&h, r->t0);
    gshgeom_hash_f64(&h, r->t1);
    gshgeom_hash_f64(&h, r->s0);
    gshgeom_hash_f64(&h, r->len);
  }
  for (size_t i = 0; i < m->sub_n; i++) {
    gshgeom_hash_f64(&h, m->sub[i].begin);
    gshgeom_hash_f64(&h, m->sub[i].end);
    gshgeom_hash_u8(&h, m->sub[i].closed ? 1 : 0);
  }
  return gshgeom_hash_value(&h);
}

static int put_f64(FILE * out, double v) {
  return fprintf(out, " %016llx", (unsigned long long)gshgeom_f64_bits(v));
}

static int put_point(FILE * out, const GSHGEOM_Point * p, int dim) {
  int n = put_f64(out, p->x);
  if (n >= 0) n = put_f64(out, p->y);
  if (n >= 0 && dim == 3) n = put_f64(out, p->z);
  return n;
}

GSHGEOM_Result gshgeom_measured_dump(const GSHGEOM_Measured * m, FILE * out) {
  if (!m || !out) {
    return GSHGEOM_ERR_INVALID;
  }
  if (fprintf(out,
          "gshgeom.measured %u dim=%d segments=%llu rows=%llu subpaths=%llu "
          "length=%016llx hash=%016llx\n",
          MEASURED_FORMAT_VERSION, m->dim, (unsigned long long)m->seg_n,
          (unsigned long long)m->tab.n, (unsigned long long)m->sub_n,
          (unsigned long long)gshgeom_f64_bits(m->tab.sum),
          (unsigned long long)gshgeom_measured_hash(m)) < 0) {
    return GSHGEOM_ERR_IO;
  }
  if (fputc('O', out) == EOF || put_point(out, &m->origin, m->dim) < 0 ||
      fputc('\n', out) == EOF) {
    return GSHGEOM_ERR_IO;
  }
  for (size_t i = 0; i < m->seg_n; i++) {
    const Seg * s = &m->seg[i];
    int n = gshgeom_curve_point_count(s->curve.kind);
    if (fprintf(out, "G %llu %d", (unsigned long long)s->element,
            (int)s->curve.kind) < 0) {
      return GSHGEOM_ERR_IO;
    }
    for (int k = 0; k < n; k++) {
      if (put_point(out, &s->curve.p[k], m->dim) < 0) {
        return GSHGEOM_ERR_IO;
      }
    }
    if (gshgeom_curve_is_rational(s->curve.kind)) {
      for (int k = 0; k < n; k++) {
        if (put_f64(out, s->curve.w[k]) < 0) {
          return GSHGEOM_ERR_IO;
        }
      }
    }
    if (fputc('\n', out) == EOF) {
      return GSHGEOM_ERR_IO;
    }
  }
  for (size_t i = 0; i < m->tab.n; i++) {
    const Row * r = &m->tab.row[i];
    if (fprintf(out, "R %llu", (unsigned long long)r->seg) < 0 ||
        put_f64(out, r->t0) < 0 || put_f64(out, r->t1) < 0 ||
        put_f64(out, r->s0) < 0 || put_f64(out, r->len) < 0 ||
        fputc('\n', out) == EOF) {
      return GSHGEOM_ERR_IO;
    }
  }
  for (size_t i = 0; i < m->sub_n; i++) {
    if (fputc('S', out) == EOF || put_f64(out, m->sub[i].begin) < 0 ||
        put_f64(out, m->sub[i].end) < 0 ||
        fprintf(out, " %d\n", m->sub[i].closed ? 1 : 0) < 0) {
      return GSHGEOM_ERR_IO;
    }
  }
  /* A buffered stream reports a write error at the flush. */
  return fflush(out) == 0 ? GSHGEOM_OK : GSHGEOM_ERR_IO;
}
