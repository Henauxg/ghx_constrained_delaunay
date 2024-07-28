# Changelog

## Version 0.2.0 (TBD)

- features: Defaults to using `f64` instead of `f32` during the triangulation for more stable results (benchmarks were already using `f64`). `f32` can still be enabled via the `f32` cargo feature.
- changes: `constrained_triangulation_from_3d_planar_vertices`, `triangulation_from_3d_planar_vertices` and `transform_to_2d_planar_coordinate_system` now take a `Vector3` for `plane_normal` instead of a `Vector3A`.
  - `Vector3A` type has been removed. It was not used except for the `plane_normal` input format, and it was already being converted immediately to a `Vector3`.

- examples: The visual debugger in the examples now always uses `glam::f32::Vec3` for vertices storage.

## Version 0.1.2 (2024-07-26)

- fix: properly handle the case where q4 is infinite in `should_swap_diagonals` for CDT (via a new `constrained_should_swap_diagonals`)
- optimizations: more small optimizations for `should_swap_diagonals` in DT and for the `swaps_history` in `restore_delaunay_triangulation_constrained` in CDT

## Version 0.1.1 (2024-07-25)

- fix: properly handle the case where q4 is infinite in `should_swap_diagonals`

## Version 0.1.0 (2024-07-24)

- Initial release