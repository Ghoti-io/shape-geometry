# Ghoti.io Shape Geometry

Paths, transforms and flattening for the shape family, bit-identical on every host. It is part of the
[Ghoti.io suite](https://github.com/Ghoti-io/suite).

## What is implemented

- Curves (line, quadratic, cubic, rational quadratic and cubic, 2D and 3D) with
  evaluation and exact-extremum bounds, and the model path (`double`, y down).
- Maps with kind tags, the one composer, the orientation rules, the exact
  `orient2d` and `orient3d`, and the grid-aligned route predicate.
- The 26.6 device path: the one rounding rule `floor(x * 64 + 0.5)`, a model
  path under a map converted to it, and the integer flattener (1/16 pixel by
  default and any tolerance of 2 units and up on request, revision 2, proved
  within the tolerance for every `int32` curve; ported from `font`'s scan
  converter, which now calls it) producing a flattened path, or the vertices of
  one curve.
- Along-path parametrization: a measured path built once from a model path (any
  curve kind, 2D or 3D, any number of subpaths) gives its arc length and, at any
  distance along it, the point, the unit tangent and the left normal. The arc
  length of a curve is adaptive eight-point Gauss-Legendre to a stated
  tolerance; the table and every answer are bit-identical on every target.
- Path vertices (`vertices.h`): each vertex of a 2D path with its position, the unit incoming and outgoing directions
  (segments without extent are looked past; a closed subpath is read round) and the direction an `orient="auto"` marker
  takes (the bisector, bit-exact on every target).
- Exact arcs (`arc.h`): an elliptical arc is one element of a 2D path, its start and end kept to the bit and its centre form
  a centre and two conjugate half-diameters (closed under every affine map), built by the SVG endpoint form (with the
  radius correction and the near-semicircle decided exactly), the Canvas and DXF centre forms, `arcTo` and the DXF bulge.
  An arc decomposes exactly into at most four rational quadratics (a quarter turn's weight is exactly `sqrt(1/2)`, a
  rounded rectangle's corner has exact control points); bounds, maps (with the circular flag), measurement, strokes,
  dashes, vertices and the device path all take arcs, and the element readback refuses one so that a reader that does not
  know arcs fails instead of drawing the chord. `examples/arc.c` builds a rounded rectangle and a circle.
- Named primitives (`shapes.h`): a rectangle, a rounded rectangle and an ellipse as lines and quarter-turn arcs with SVG 2's start
  point and direction, the four Canvas `rect` sign cases as two flips (and `roundRect` with a negative size as Chromium draws it: the radii to the mirrored corners, the path reversed for one flip), radii that overlap refused on an exact sum, and every coordinate
  the correctly rounded value of its expression. `examples/shapes.c` builds them.
- The conic route: the device path of a rational curve or an arc, proved within the tolerance less 0.29 units (flattener
  revision 3).
- The fill rules (`fill.h`: nonzero 0, even-odd 1, pinned), which `shape-cpu`
  fills under and `font`'s `GFNT_FillRule` is the type of.
- Strokes (`stroke.h`): the region a stroke covers, stated once as rules SG1 to
  SG30 after SVG 2 and the WHATWG canvas text (joins miter, round, bevel,
  miter-clip and arcs, caps, the miter limit decided exactly, local and device
  width, degenerate subpaths, inside and outside alignment as separate parts),
  and `gshgeom_stroke_expand`, which turns a centerline, the parameters and a map
  into 26.6 polygons whose nonzero fill contains that region and lies within the
  tolerance of it (proved; checked against an independent geometric oracle).
  Dashes cut the centerline first (SVG 2's dash positions, measured with the along-path
  parametrization), a DXF simple linetype lowers to a dash array, and `examples/stroke.c` shows a stroke.
- Stroke curves (`stroke_curves.h`): the same region as an outline of lines, cubics and exact
  rational-quadratic arcs in target units, with a tolerance stated in those units (rules SG32 to SG43), for a writer that
  must expand a stroke its format cannot express and write it as a fill that stays exact under zoom;
  `gshgeom_stroke_curves_lower_cubic` turns a path of curves into lines and cubics within a stated bound, and
  `examples/stroke_curves.c` writes a stroke as SVG path data.
- `_dump` and `_hash` for each (FNV-1a 64 over little-endian fields, the same
  value on every target).

`documentation/design.md` has the model.

## Building

```bash
make test PREFIX="$PWD/.local"
make install PREFIX="$PWD/.local"
```

Dependencies are found by pkg-config and nothing else; a missing one is an
error that names the fix. The libraries this one needs: cutil, math.

Determinism is a requirement, not a hope: every `double` stage builds with
`-ffp-contract=off` and no x87, and elementary functions come from `math`, not
libm. `make check-gates` shows each of those rules failing on a planted
violation.

## License

LGPL-3.0-only. See `COPYING` and `COPYING.LESSER`. Patches are not being
accepted at this time; see `CONTRIBUTING.md`.
