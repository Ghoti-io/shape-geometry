/* Throwaway generator of known answers from font's own flattener.
 *
 * Built outside every repository: it #includes font's src/raster/flatten.c itself
 * (the unmodified file; see make.sh) and calls gfnt_flatten_quad and
 * gfnt_flatten_cubic at tolerance 4 with an emit callback that records every
 * vertex.  Prints tests/flatten_vectors.h of shape-geometry.
 */
#include "src/raster/flatten.c"
#include <stdio.h>
#include <stdlib.h>

static GFNT_Point rec[1 << 20];
static unsigned rec_n;
static void record(void * user, GFNT_Point to, uint32_t leaf, int level) {
  (void)user; (void)leaf; (void)level;
  rec[rec_n++] = to;
}

static uint64_t rs = UINT64_C(0x2545F4914F6CDD1D);
static uint64_t nx(void) { rs ^= rs >> 12; rs ^= rs << 25; rs ^= rs >> 27; return rs * UINT64_C(2685821657736338717); }
static int32_t rr(int bits) { /* uniform in (-2^bits, 2^bits) */
  int64_t span = (int64_t)1 << (bits + 1);
  return (int32_t)((int64_t)(nx() % (uint64_t)span) - ((int64_t)1 << bits));
}

static uint64_t fnv_i32(uint64_t h, int32_t v) {
  uint32_t u = (uint32_t)v;
  for (int i = 0; i < 4; i++) { h ^= (u >> (8 * i)) & 0xff; h *= UINT64_C(0x100000001b3); }
  return h;
}

typedef struct { int kind; int32_t p[4][2]; } Vec;

static uint64_t run(const Vec *v, unsigned *count) {
  rec_n = 0;
  GFNT_Point p[4];
  for (int i = 0; i < 4; i++) { p[i].x = v->p[i][0]; p[i].y = v->p[i][1]; }
  int forced;
  if (v->kind == 2) forced = gfnt_flatten_quad(4, p[0], p[1], p[2], record, 0);
  else forced = gfnt_flatten_cubic(4, p[0], p[1], p[2], p[3], record, 0);
  if (forced) { fprintf(stderr, "the depth cap forced a piece\n"); exit(1); }
  uint64_t h = UINT64_C(0xcbf29ce484222325);
  for (unsigned i = 0; i < rec_n; i++) { h = fnv_i32(h, rec[i].x); h = fnv_i32(h, rec[i].y); }
  *count = rec_n;
  return h;
}

static Vec mk2(int32_t a, int32_t b, int32_t c, int32_t d, int32_t e, int32_t f) {
  Vec v = {2, {{a, b}, {c, d}, {e, f}, {0, 0}}}; return v;
}
static Vec mk3(int32_t a, int32_t b, int32_t c, int32_t d, int32_t e, int32_t f, int32_t g, int32_t h) {
  Vec v = {3, {{a, b}, {c, d}, {e, f}, {g, h}}}; return v;
}

int main(void) {
  Vec ex[256]; int nex = 0;
  /* hand-picked, explicit lists (their vertices are font's; the thresholds
   in the comments are in units: a piece is flat at 4 * (t - 1) = 12) */
  ex[nex++] = mk2(0, 0, 32, 0, 64, 0);          /* straight, flat */
  ex[nex++] = mk2(0, 0, 0, 0, 0, 0);            /* a point */
  ex[nex++] = mk2(10, 10, 10, 10, 10, 10);
  ex[nex++] = mk2(0, 0, 64, 64, 128, 0);        /* a hump of 32 units */
  ex[nex++] = mk2(0, 0, 640, 640, 1280, 0);
  ex[nex++] = mk2(0, 0, 100, -3000, 200, 0);
  ex[nex++] = mk2(0, 0, 5000, 0, 0, 0);         /* out and back */
  ex[nex++] = mk2(-640, -640, 0, 640, 640, -640);
  ex[nex++] = mk2(-1, -1, 0, 0, 1, 1);
  ex[nex++] = mk2(7, 7, 7, 7, 700, 7);          /* p1 = p0 */
  ex[nex++] = mk2(0, 0, 33, 33, 64, 0);         /* d = 66: well over the threshold */
  ex[nex++] = mk2(0, 0, 6, 6, 12, 0);           /* d = 12: at the threshold */
  ex[nex++] = mk2(0, 0, 6, 7, 12, 0);           /* d = 14: just over */
  ex[nex++] = mk2(0, 0, 1 << 24, 1 << 24, 1 << 25, 0);
  ex[nex++] = mk2(-(1 << 26), -(1 << 26), 0, 1 << 26, 1 << 26, -(1 << 26));
  ex[nex++] = mk3(0, 0, 0, 0, 0, 0, 0, 0);
  ex[nex++] = mk3(0, 0, 21, 0, 43, 0, 64, 0);   /* straight cubic */
  ex[nex++] = mk3(0, 0, 0, 640, 640, 640, 640, 0);
  ex[nex++] = mk3(0, 0, 640, 640, 0, 640, 640, 0);     /* loop-free S */
  ex[nex++] = mk3(0, 0, 1280, 640, -640, 640, 640, 0); /* a loop */
  ex[nex++] = mk3(0, 0, 640, 640, 0, 640, 640, 0);
  ex[nex++] = mk3(0, 0, 640, 0, 0, 0, 640, 0);         /* a cusp on a line */
  ex[nex++] = mk3(0, 0, 2000, 2000, -2000, 2000, 0, 0);   /* closed: p0 = p3 */
  ex[nex++] = mk3(5, 5, 5, 5, 200, 200, 200, 200);   /* p1 = p0, p2 = p3 */
  ex[nex++] = mk3(-300, 120, -1000, 5000, 4000, -5000, 300, -120);
  ex[nex++] = mk3(0, 0, 12, 12, 24, 12, 36, 0);        /* d1, d2 both small */
  ex[nex++] = mk3(0, 0, 1 << 22, 1 << 22, -(1 << 22), 1 << 22, 1 << 21, 0);
  ex[nex++] = mk3(-(1 << 26), 0, -(1 << 26), 1 << 26, 1 << 26, 1 << 26, 1 << 26, 0);
  while (nex < 60) {
    int bits = 4 + (int)(nex % 8);
    if (nex & 1) ex[nex] = mk2(rr(bits), rr(bits), rr(bits), rr(bits), rr(bits), rr(bits));
    else ex[nex] = mk3(rr(bits), rr(bits), rr(bits), rr(bits), rr(bits), rr(bits), rr(bits), rr(bits));
    nex++;
  }
  printf("// Generated: known answers of font's flattener (libs/font/src/raster/\n"
         "// flatten.c, gfnt_flatten_quad and gfnt_flatten_cubic at tolerance 4, font commit\n"
         "// 69d1860), recorded by tools/vectors-from-font (story 7a).  Do not edit.\n"
         "// Each vector is a quadratic (kind 2) or cubic (kind 3) in 26.6, the number of\n"
         "// lines font emits, the FNV-1a 64 of their end points (x then y, int32 little-endian),\n"
         "// and, for those with an offset, the end points themselves.\n\n"
         "#ifndef GHOTI_IO_GSHGEOM_FLATTEN_VECTORS_H\n#define GHOTI_IO_GSHGEOM_FLATTEN_VECTORS_H\n\n#include <cstdint>\n\n"
         "// The flattener revision these vertices are those of.\nstatic const int kFlattenVectorsRevision = 1;\n\n"
         "struct FlattenVector {\n  int kind;\n  std::int32_t p[4][2];\n  unsigned count;\n  std::uint64_t hash;\n  unsigned offset;  // into kFlattenPoints; UINT32_MAX when only the hash is kept\n};\n\n");
  /* explicit points first */
  static int32_t pts[1 << 20]; unsigned npts = 0;
  static Vec all[8192]; static unsigned off[8192]; static unsigned cnt[8192]; static uint64_t hs[8192]; unsigned n = 0;
  for (int i = 0; i < nex; i++) {
    unsigned c; uint64_t h = run(&ex[i], &c);
    all[n] = ex[i]; cnt[n] = c; hs[n] = h;
    if (c > 150) { off[n] = UINT32_MAX; n++; continue; }
    off[n] = npts; n++;
    for (unsigned k = 0; k < c; k++) { pts[npts++] = rec[k].x; pts[npts++] = rec[k].y; }
  }
  
  unsigned maxc = 0;
  for (int i = 0; i < 900; i++) {
    int bits = 4 + (int)(nx() % 24);          /* up to 2^27 */
    Vec v = (i & 1) ? mk2(rr(bits), rr(bits), rr(bits), rr(bits), rr(bits), rr(bits))
                    : mk3(rr(bits), rr(bits), rr(bits), rr(bits), rr(bits), rr(bits), rr(bits), rr(bits));
    unsigned c; uint64_t h = run(&v, &c);
    if (c > maxc) maxc = c;
    all[n] = v; off[n] = UINT32_MAX; cnt[n] = c; hs[n] = h; n++;
  }
  printf("static const unsigned kFlattenVectorCount = %u;\n\n", n);
  printf("static const std::int32_t kFlattenPoints[] = {");
  for (unsigned i = 0; i < npts; i++) printf("%s%d,", (i % 12) ? " " : "\n    ", pts[i]);
  printf("\n};\n\n");
  printf("static const FlattenVector kFlattenVectors[] = {\n");
  for (unsigned i = 0; i < n; i++) {
    printf("    {%d, {{%d, %d}, {%d, %d}, {%d, %d}, {%d, %d}}, %u, UINT64_C(0x%016llx), %s%u%s},\n", all[i].kind,
      all[i].p[0][0], all[i].p[0][1], all[i].p[1][0], all[i].p[1][1], all[i].p[2][0], all[i].p[2][1], all[i].p[3][0], all[i].p[3][1],
      cnt[i], (unsigned long long)hs[i], off[i] == UINT32_MAX ? "" : "", off[i] == UINT32_MAX ? 0xffffffffu : off[i], off[i] == UINT32_MAX ? "u" : "");
  }
  printf("};\n\n#endif\n");
  fprintf(stderr, "vectors %u points %u max lines %u\n", n, npts / 2, maxc);
  return 0;
}
