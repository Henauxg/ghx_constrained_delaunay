# Changelog

## Version 0.3.0 (2026-11-16)

### Fixes

- fix: a vertex close to the line through the finite edge of a nearly flat semi-infinite triangle, but outside that edge, was inserted on the edge. This created folded triangles and broken neighbor links, and the CDT could then fail with "aborting cycle when searching for the first intersection of an edge" (#6).
- fix: the orientation and in-circle tests now use exact predicates (`robust` crate), so their sign is always correct. This fixes most failures on nearly colinear inputs, such as points a few ULPs off a straight constrained border.
- fix: a vertex inside the convex hull, but closer than `Float::EPSILON` to a hull edge, could stay on the hull.
- fix: a vertex inside the convex hull could also stay on the boundary of the output when it was the most extreme vertex in two of the four infinite directions, so a few triangles along the hull were missing. Triangles with two infinite vertices are now swapped when a new vertex is further than their finite vertex in one of these directions.

### Changes

- changed: a vertex is only placed on an existing edge when it is exactly on that edge's line, instead of within `Float::EPSILON` of it. Vertices too close to an existing vertex are still merged.
- changed: the output triangles come in a different order than in 0.2.0, because the insertion order changed. Where 4 or more vertices lie on the same circle, the triangulation can also contain a different, equally valid, set of triangles.
- changed: `glam` versions up to 0.30 are accepted. `hashbrown` is updated to 0.16 and `thiserror` to 2.0 (#3).
- added: `triangle_orientation` field in `TriangulationConfiguration` and `ConstrainedTriangulationConfiguration`. It sets the vertex order of the output triangles: `TriangleOrientation::Clockwise` (default) or `TriangleOrientation::CounterClockwise` (#4).
- removed: `bin_vertex_density_power` from `TriangulationConfiguration` and `ConstrainedTriangulationConfiguration`, and `DEFAULT_BIN_VERTEX_DENSITY_POWER`. The insertion order no longer uses bins.
- removed: the geometry helpers in `utils` are no longer public: `egdes_intersect`, `EdgesIntersectionResult`, `triplet_orientation`, `Orientation`, `test_point_edge_side`, `EdgeSideTestResult`, `is_point_on_right_side_of_edge`, `is_point_strictly_on_right_side_of_edge`, `line_slope` and `on_segment`. `Quad::diagonals_intersection_test` is no longer public either.

### Optimizations

- improved: in CDT, a constrained edge that already exists in the triangulation is found by walking around its first vertex, without geometric tests.
- improved: in CDT, the final filtering reads constrained edges from a flag per triangle edge instead of a hash set. Together with the previous change, CDT is 8 to 27% faster on the coastline benchmarks.
- improved: in CDT, constrained edges are marked while the constraints are applied, instead of in a second pass over all of them. About 8% faster on the coastline benchmarks.
- improved: vertices are inserted in BRIO order (random rounds of doubling size, each following a Hilbert curve) instead of in bins. The Delaunay triangulation is 5 to 20% faster on the benchmarks, with no setting to tune.
- improved: the insertion order is computed about 40% faster. The radix sort counts are taken while computing the Hilbert keys, and the key computation uses vector instructions. The Delaunay triangulation is 3 to 5% faster.
- improved: from 50,000 vertices, the vertices are inserted from a copy sorted in insertion order, so that the vertices of nearby triangles are close in memory. On inputs whose vertices come in random order, DT is about 30% faster and CDT 12 to 20% faster. Inputs that are already in spatial order, such as coastlines, are unchanged.

### Examples, benchmarks & tests

- added: a fuzzer (`fuzzer` crate). It runs the CDT on random near-degenerate boxes, and checks the DT of random points against their exact convex hull and the Delaunay property.
- added: the `cdt_random_points_with_loops` benchmark example: CDT on random points with constrained loops, where most constrained edges are missing from the Delaunay triangulation.
- updated: benchmark results, with CDT tables and new DT charts.

## Version 0.2.0 (2025-03-13)

### Fixes 

- fix: In CDT, update internal function `remove_crossed_edges` to be faillible, to avoid a possible infinite loop.
  
### Changes

- changed: defaults to using `f64` instead of `f32` inside the algorithm for more stable results (benchmarks were already using `f64`). `f32` can still be enabled via the `f32` cargo feature.
- changed: `tracing` dependency is now optional. Error messages (and `DebugContext` when the `debug_context` feature is enabled) are returned within a `TriangulationError`.
- changed: `glam` dependency now accepts version of glam up to 0.29
- changed: now accept slices as inputs instead of `Vec` in most functions
  
- vertices format:
  - added: `Vertex2d` and `Vertex3d` traits defining the algorithm vertex input format.
  - changed: `constrained_triangulation_from_3d_planar_vertices`, `constrained_triangulation_from_2d_vertices`, `triangulation_from_3d_planar_vertices`, `triangulation_from_2d_vertices` and `transform_to_2d_planar_coordinate_system` to be generic over the input vertex format. This does not affect the vertex format used during the algorithm. This allows to decouple the vertex format used in the algorithm (`f32` or `f64`) from the vertex input format.
    - As an example: it is now possible to give vertices with f32 coordinates as an input, but compute the triangulation using f64.
  - removed: `Vector3A` type. It was not used except for the `plane_normal` input format, and it was already being converted immediately to a `Vector3`.
  - removed: the `Vector3` type now superseded by the `Vertex3d` trait.
  
### Optimizations

- improved: find_vertex_placement hot path for edge side tests, and move the vertex merging to the end
- improved: some triangles neigbor checks

### Examples, benchmarks & profiling

- changed: the visual debugger in the examples now always uses `glam::f32::Vec3` for vertices storage.
- updated: profiling & benchmarks

## Version 0.1.2 (2024-07-26)

- fix: properly handle the case where q4 is infinite in `should_swap_diagonals` for CDT (via a new `constrained_should_swap_diagonals`)
- optimizations: more small optimizations for `should_swap_diagonals` in DT and for the `swaps_history` in `restore_delaunay_triangulation_constrained` in CDT

## Version 0.1.1 (2024-07-25)

- fix: properly handle the case where q4 is infinite in `should_swap_diagonals`

## Version 0.1.0 (2024-07-24)

- Initial release