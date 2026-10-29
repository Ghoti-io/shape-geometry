/**
 * @file
 *
 * One case per audited site that reads a path by verb (tools/verb-walkers.audit): each is fed a path that
 * holds an arc.  A site that handles an arc returns its answer; a site that does not know arcs yet receives
 * ::GSHGEOM_ERR_UNSUPPORTED and passes it on, writes nothing, and never draws the chord.  A site that cannot be
 * reached by an arc is shown unreachable by building the path that would reach it.
 *
 * Copyright 2026 by Corey Pennycuff
 */

#include "../test_helpers.h"

#include <cstdint>
#include <cstring>
#include <string>

namespace {

GSHGEOM_Point P(double x, double y) {
  return GSHGEOM_Point{x, y, 0.0};
}

// move to (10, 0), a quarter circle about the origin to (0, 10), then a line back to the origin
GSHGEOM_Path * arc_path() {
  GSHGEOM_Path * p = nullptr;
  EXPECT_EQ(gshgeom_path_create(nullptr, nullptr, 2, &p), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_move_to(p, P(10, 0)), GSHGEOM_OK);
  GSHGEOM_Arc_Form f = {P(0, 0), P(10, 0), P(0, 10), GSHGEOM_ANGLE_DEGREES, 0, 90, GSHGEOM_ARC_CIRCULAR};
  EXPECT_EQ(gshgeom_path_arc_conjugate_to(p, &f, nullptr), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_line_to(p, P(0, 0)), GSHGEOM_OK);
  return p;
}

}  // namespace

// path.c: gshgeom_path_element, _bounds, _hash, _dump, _copy, _equal, _verb and the tables (updated, AR1, AR10, AR15)
TEST(Walkers, PathElementBoundsHashDumpCopyEqualAndVerbKnowTheArc) {
  GSHGEOM_Path * p = arc_path();
  GSHGEOM_Element el;
  EXPECT_EQ(gshgeom_path_element(p, 1, &el), GSHGEOM_ERR_UNSUPPORTED);
  GSHGEOM_Verb v;
  ASSERT_EQ(gshgeom_path_verb(p, 1, &v), GSHGEOM_OK);
  EXPECT_EQ(v, GSHGEOM_VERB_ARC);
  GSHGEOM_Bounds b;
  ASSERT_EQ(gshgeom_path_bounds(p, &b), GSHGEOM_OK);
  EXPECT_LE(b.min.x, 0.0);
  EXPECT_GE(b.max.x, 10.0);
  EXPECT_GE(b.max.y, 10.0);
  GSHGEOM_Path * c = nullptr;
  ASSERT_EQ(gshgeom_path_copy(p, nullptr, nullptr, &c), GSHGEOM_OK);
  EXPECT_TRUE(gshgeom_path_equal(p, c));
  EXPECT_EQ(gshgeom_path_hash(p), gshgeom_path_hash(c));
  char * buf = nullptr;
  size_t n = 0;
  FILE * f = open_memstream(&buf, &n);
  ASSERT_EQ(gshgeom_path_dump(p, f), GSHGEOM_OK);
  fclose(f);
  EXPECT_NE(std::string(buf, n).find("\nA "), std::string::npos);
  free(buf);
  gshgeom_path_destroy(p);
  gshgeom_path_destroy(c);
}

// path.c's builders, shape-cpu's render/sink.c and render/stroke.c (element counts): an arc is one element (AR16)
TEST(Walkers, ElementCountCountsAnArcOnce) {
  GSHGEOM_Path * p = arc_path();
  EXPECT_EQ(gshgeom_path_element_count(p), 3u);
  GSHGEOM_Path * lowered = nullptr;
  ASSERT_EQ(gshgeom_path_lower_arcs(p, nullptr, nullptr, &lowered), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_element_count(lowered), 3u);  // a quarter turn is one piece
  size_t count = 0;
  EXPECT_EQ(gshgeom_path_vertex_count(p, &count), GSHGEOM_OK);
  EXPECT_EQ(count, 3u);
  gshgeom_path_destroy(p);
  gshgeom_path_destroy(lowered);
}

static GSHGEOM_Result ignore_vertex(const GSHGEOM_Vertex *, void *) {
  return GSHGEOM_OK;
}

// vertices.c (read_seg, element, gshgeom_path_vertices) and vector's marker lowering that calls it: refused loudly
TEST(Walkers, VerticesRefuseAnArcUntilItReadsOne) {
  GSHGEOM_Path * p = arc_path();
  EXPECT_EQ(gshgeom_path_vertices(p, ignore_vertex, nullptr), GSHGEOM_ERR_UNSUPPORTED);
  GSHGEOM_Path * lowered = nullptr;
  ASSERT_EQ(gshgeom_path_lower_arcs(p, nullptr, nullptr, &lowered), GSHGEOM_OK);
  EXPECT_EQ(gshgeom_path_vertices(lowered, ignore_vertex, nullptr), GSHGEOM_OK);  // the lowered path is read
  gshgeom_path_destroy(p);
  gshgeom_path_destroy(lowered);
}

// along.c: gshgeom_measured_create
TEST(Walkers, MeasuredCreateRefusesAnArcAndMeasuresTheLoweredPath) {
  GSHGEOM_Path * p = arc_path();
  GSHGEOM_Measured * m = reinterpret_cast<GSHGEOM_Measured *>(0x1);
  EXPECT_EQ(gshgeom_measured_create(nullptr, nullptr, p, 0.001, &m), GSHGEOM_ERR_UNSUPPORTED);
  EXPECT_EQ(m, reinterpret_cast<GSHGEOM_Measured *>(0x1));
  GSHGEOM_Path * lowered = nullptr;
  ASSERT_EQ(gshgeom_path_lower_arcs(p, nullptr, nullptr, &lowered), GSHGEOM_OK);
  GSHGEOM_Measured * ml = nullptr;
  ASSERT_EQ(gshgeom_measured_create(nullptr, nullptr, lowered, 0.001, &ml), GSHGEOM_OK);
  gshgeom_measured_destroy(ml);
  gshgeom_path_destroy(p);
  gshgeom_path_destroy(lowered);
}

// stroke/expand.c (read_centerline, build_clip, build_clip_curves) and stroke/dash.c (gshgeom_sk_dash_begin)
TEST(Walkers, StrokeExpansionRefusesAnArcSolidOrDashed) {
  GSHGEOM_Path * p = arc_path();
  GSHGEOM_Map m;
  ASSERT_EQ(gshgeom_map_identity(2, &m), GSHGEOM_OK);
  const double dashes[2] = {3.0, 2.0};
  for (int dashed = 0; dashed < 2; dashed++) {
    GSHGEOM_Stroke s;
    gshgeom_stroke_init(&s);
    s.width = 2.0;
    if (dashed) {
      s.dashes = dashes;
      s.dash_count = 2;
    }
    GSHGEOM_Stroke_Outline * out = reinterpret_cast<GSHGEOM_Stroke_Outline *>(0x1);
    EXPECT_EQ(gshgeom_stroke_expand(p, &s, &m, nullptr, nullptr, &out), GSHGEOM_ERR_UNSUPPORTED) << dashed;
    EXPECT_EQ(out, reinterpret_cast<GSHGEOM_Stroke_Outline *>(0x1));
  }
  GSHGEOM_Path * lowered = nullptr;
  ASSERT_EQ(gshgeom_path_lower_arcs(p, nullptr, nullptr, &lowered), GSHGEOM_OK);
  GSHGEOM_Stroke s;
  gshgeom_stroke_init(&s);
  s.width = 2.0;
  GSHGEOM_Stroke_Outline * out = nullptr;
  ASSERT_EQ(gshgeom_stroke_expand(lowered, &s, &m, nullptr, nullptr, &out), GSHGEOM_OK);
  gshgeom_stroke_outline_destroy(out);
  gshgeom_path_destroy(p);
  gshgeom_path_destroy(lowered);
}

// stroke_curves/lower.c: gshgeom_stroke_curves_lower_cubic
TEST(Walkers, StrokeCurvesLowerRefusesAnArc) {
  GSHGEOM_Path * p = arc_path();
  GSHGEOM_Path * out = reinterpret_cast<GSHGEOM_Path *>(0x1);
  EXPECT_EQ(gshgeom_stroke_curves_lower_cubic(p, 0.01, nullptr, nullptr, &out), GSHGEOM_ERR_UNSUPPORTED);
  EXPECT_EQ(out, reinterpret_cast<GSHGEOM_Path *>(0x1));
  gshgeom_path_destroy(p);
}

// device/device.c (convert_element, gshgeom_device_path_from_model), flatten/flatten.c, shape-cpu's render/geometry.c
TEST(Walkers, DevicePathFromModelAndFlattenRefuseAnArc) {
  GSHGEOM_Path * p = arc_path();
  GSHGEOM_Map m;
  ASSERT_EQ(gshgeom_map_identity(2, &m), GSHGEOM_OK);
  GSHGEOM_Device_Path * dp = reinterpret_cast<GSHGEOM_Device_Path *>(0x1);
  EXPECT_EQ(gshgeom_device_path_from_model(p, &m, nullptr, nullptr, &dp), GSHGEOM_ERR_UNSUPPORTED);
  EXPECT_EQ(dp, reinterpret_cast<GSHGEOM_Device_Path *>(0x1));
  GSHGEOM_Path * lowered = nullptr;
  ASSERT_EQ(gshgeom_path_lower_arcs(p, nullptr, nullptr, &lowered), GSHGEOM_OK);
  GSHGEOM_Device_Path * ok = nullptr;
  EXPECT_EQ(gshgeom_device_path_from_model(lowered, &m, nullptr, nullptr, &ok), GSHGEOM_OK);
  gshgeom_device_path_destroy(ok);
  gshgeom_path_destroy(p);
  gshgeom_path_destroy(lowered);
}

// shape-scene (resource.c, glyph.c, image_paint.c) read a path through gshgeom_path_element; A2 edits them.  Until then
// they receive this error from the one readback and so fail rather than copy, intern, map or paint a chord.
TEST(Walkers, TheReadbackTheSceneLibraryUsesRefusesAnArc) {
  GSHGEOM_Path * p = arc_path();
  for (size_t i = 0; i < gshgeom_path_element_count(p); i++) {
    GSHGEOM_Element el;
    GSHGEOM_Result r = gshgeom_path_element(p, i, &el);
    GSHGEOM_Verb v;
    ASSERT_EQ(gshgeom_path_verb(p, i, &v), GSHGEOM_OK);
    EXPECT_EQ(r, v == GSHGEOM_VERB_ARC ? GSHGEOM_ERR_UNSUPPORTED : GSHGEOM_OK) << i;
  }
  gshgeom_path_destroy(p);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
