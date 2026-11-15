# Changelog

## Unreleased

### Fixes

- fix: a vertex close to the line through the finite edge of a nearly flat semi-infinite triangle, but outside that edge, was inserted on the edge. This created folded triangles and broken neighbor links, and the CDT could then fail with "aborting cycle when searching for the first intersection of an edge" (#6).
- fix: the orientation and in-circle tests now use exact predicates (`robust` crate), so their sign is always correct. This fixes most failures on nearly colinear inputs, such as points a few ULPs off a straight constrained border.
- fix: a vertex inside the convex hull, but closer than `Float::EPSILON` to a hull edge, could stay on the hull.

### Changes

- changed: a vertex is only placed on an existing edge when it is exactly on that edge's line, instead of within `Float::EPSILON` of it. Vertices too close to an existing vertex are still merged.
- added: `triangle_orientation` field in `TriangulationConfiguration` and `ConstrainedTriangulationConfiguration`. It sets the vertex order of the output triangles: `TriangleOrientation::Clockwise` (default) or `TriangleOrientation::CounterClockwise` (#4).
- removed: the geometry helpers in `utils` are no longer public: `egdes_intersect`, `EdgesIntersectionResult`, `triplet_orientation`, `Orientation`, `test_point_edge_side`, `EdgeSideTestResult`, `is_point_on_right_side_of_edge`, `is_point_strictly_on_right_side_of_edge`, `line_slope` and `on_segment`. `Quad::diagonals_intersection_test` is no longer public either.

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