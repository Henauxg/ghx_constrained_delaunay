use crate::infinite::{
    collect_infinite_triangle_vertices, is_finite, is_infinite, is_vertex_beyond_infinite_corner,
    is_vertex_in_half_plane_1, vertex_placement_1_infinite_vertex,
    vertex_placement_2_infinite_vertex, INFINITE_V0_ID, INFINITE_V1_ID, INFINITE_V2_ID,
    INFINITE_V3_ID,
};
use crate::types::{
    next_clockwise_edge_index, next_counter_clockwise_edge_index, opposite_vertex_index,
    vertex_next_ccw_edge_index, vertex_next_cw_edge_index, Edge, Float, Neighbor, Quad,
    TriangleData, TriangleEdgeIndex, TriangleId, Triangles, Vertex, Vertex2d, Vertex3d, VertexId,
    EDGE_12, EDGE_23, EDGE_31, VERT_1, VERT_2,
};
use crate::utils::{is_vertex_in_triangle_circumcircle, test_point_edge_side};

#[cfg(feature = "parallel_filtering")]
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};

#[cfg(feature = "progress_log")]
use tracing::info;

#[cfg(feature = "debug_context")]
use crate::debug::{DebugConfiguration, DebugContext, EventInfo, Phase};

#[cfg(feature = "profile_traces")]
use tracing::{span, Level};

/// Good default value for `filter_parallel_tri_count_threshold` in [`TriangulationConfiguration`]
pub const DEFAULT_ILTER_PARALLEL_TRI_COUNT_THRESHOLD: usize = 100_000;
/// Good default value for `filter_parallel_min_batch_len` in [`TriangulationConfiguration`]
pub const DEFAULT_FILTER_PARALLEL_MIN_BATCH_LEN: usize = 1000;

/// Order of the 3 vertices of each triangle in a [`Triangulation`].
///
/// The orientation is measured with the y axis pointing up. With the 3d functions, it is seen from the side the plane normal points to.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum TriangleOrientation {
    #[default]
    Clockwise,
    CounterClockwise,
}
impl TriangleOrientation {
    /// `verts` MUST be in clockwise order
    #[inline(always)]
    pub(crate) fn orient(self, verts: [VertexId; 3]) -> [VertexId; 3] {
        match self {
            TriangleOrientation::Clockwise => verts,
            TriangleOrientation::CounterClockwise => [verts[0], verts[2], verts[1]],
        }
    }
}

#[derive(Clone, Debug)]
pub struct TriangulationConfiguration {
    /// Order of the vertices of each triangle in the result.
    pub triangle_orientation: TriangleOrientation,

    /// Minimum count of triangles that needs to be reached in order to enable parallel filtering.
    #[cfg(feature = "parallel_filtering")]
    pub filter_parallel_tri_count_threshold: usize,
    /// Minimum size of the work batches during parallel triangles filtering.
    #[cfg(feature = "parallel_filtering")]
    pub filter_parallel_min_batch_len: usize,

    #[cfg(feature = "debug_context")]
    pub debug_config: DebugConfiguration,
}
// Not derivable with the `parallel_filtering` feature
#[allow(clippy::derivable_impls)]
impl Default for TriangulationConfiguration {
    fn default() -> Self {
        Self {
            triangle_orientation: TriangleOrientation::default(),

            #[cfg(feature = "parallel_filtering")]
            filter_parallel_tri_count_threshold: DEFAULT_ILTER_PARALLEL_TRI_COUNT_THRESHOLD,
            #[cfg(feature = "parallel_filtering")]
            filter_parallel_min_batch_len: DEFAULT_FILTER_PARALLEL_MIN_BATCH_LEN,

            #[cfg(feature = "debug_context")]
            debug_config: DebugConfiguration::default(),
        }
    }
}

#[derive(Default)]
pub struct Triangulation {
    /// Indices of the original vertices by groups of 3 to from triangles.
    pub triangles: Vec<[VertexId; 3]>,

    #[cfg(feature = "debug_context")]
    pub debug_context: DebugContext,
}

/// An internal error occured, which could be due to an invalid input, or possibly an error in the algorithm.
#[derive(thiserror::Error, Debug)]
#[error("triangulation error")]
pub struct TriangulationError {
    pub msg: String,
    #[cfg(feature = "debug_context")]
    pub debug_context: DebugContext,
}
impl TriangulationError {
    pub(crate) fn new(
        msg: String,
        #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
    ) -> Self {
        Self {
            msg,
            #[cfg(feature = "debug_context")]
            debug_context: debug_context.clone(),
        }
    }
}

/// Same as [triangulation_from_2d_vertices] but input vertices are in 3d and will be transformed to 2d before the triangulation.
///
/// Additional requirements:
/// - All the vertices are expected to belong to the same 2d plane, with the provided `plane_normal`.
/// - `plane_normal` must be normalized
pub fn triangulation_from_3d_planar_vertices<T: Vertex3d>(
    vertices: &[T],
    plane_normal: T,
    config: TriangulationConfiguration,
) -> Result<Triangulation, TriangulationError> {
    if vertices.len() < 3 {
        return Ok(Triangulation::default());
    }
    let mut planar_vertices = transform_to_2d_planar_coordinate_system(vertices, plane_normal);
    triangulation_from_2d_vertices(&mut planar_vertices, config)
}

/// Creates a Delaunay triangulation with the input vertices.
///
/// Vertices that are identical (or extremely close to one another) will me "merged" together. This means that only one of those will appear in the final triangulation.
///
/// Vertices requirements:
/// - Vertices are expected to be valid floating points values. You can use [crate::utils::validate_vertices] to check your vertices beforehand for NaN or infinity.
pub fn triangulation_from_2d_vertices<T: Vertex2d>(
    vertices: &[T],
    config: TriangulationConfiguration,
) -> Result<Triangulation, TriangulationError> {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "triangulation_from_2d_vertices").entered();

    if vertices.len() < 3 {
        return Ok(Triangulation::default());
    }

    // Uniformly scale the coordinates of the points so that they all lie between 0 and 1.
    let (normalized_vertices, _scale_factor, _x_min, _y_min) =
        normalize_vertices_coordinates(vertices);

    #[cfg(feature = "debug_context")]
    let mut debug_context =
        DebugContext::new(config.debug_config.clone(), _scale_factor, _x_min, _y_min);

    let triangles = wrap_and_triangulate_2d_normalized_vertices(
        &normalized_vertices,
        &mut None,
        #[cfg(feature = "debug_context")]
        &mut debug_context,
    )?;

    let vert_indices = filter_triangles(
        &triangles,
        &config,
        #[cfg(feature = "debug_context")]
        &mut debug_context,
    );

    Ok(Triangulation {
        triangles: vert_indices,
        #[cfg(feature = "debug_context")]
        debug_context,
    })
}

/// Transforms 3d coordinates of all vertices into 2d coordinates on a plane defined by the given normal and vertices.
/// - Input vertices need to all belong to the same 3d plan
/// - There must be at least two vertices
pub fn transform_to_2d_planar_coordinate_system<T: Vertex3d>(
    vertices: &[T],
    plane_normal: T,
) -> Vec<Vertex> {
    // Create a base B for the plan, using the first two vertices as the first base vector
    let basis_1 = (vertices[0] - vertices[1]).normalize();
    let basis_2 = basis_1.cross(-plane_normal);
    // basis_3 would be plane_normal

    // Transform every vertices into base B 2D coordinates
    let mut vertices_2d = Vec::with_capacity(vertices.len());
    for vertex in vertices {
        vertices_2d.push(Vertex::new(vertex.dot(basis_1), vertex.dot(basis_2)));
    }
    vertices_2d
}

/// This scaling ensures that all of the coordinates are between 0 and 1 but does not modify the relative positions of the points in the x-y plane.
/// The use of normalized coordinates, although not essential, reduces the effects of roundoff error and is also convenient from a computational point of view.
pub(crate) fn normalize_vertices_coordinates<T: Vertex2d>(
    vertices: &[T],
) -> (Vec<Vertex>, Float, Float, Float) {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "normalize_vertices_coordinates").entered();

    let mut normalized_vertices = Vec::with_capacity(vertices.len());
    let (mut x_min, mut y_min, mut x_max, mut y_max) =
        (Float::MAX, Float::MAX, Float::MIN, Float::MIN);

    for vertex in vertices.iter() {
        if vertex.x() < x_min {
            x_min = vertex.x();
        }
        if vertex.x() > x_max {
            x_max = vertex.x();
        }
        if vertex.y() < y_min {
            y_min = vertex.y();
        }
        if vertex.y() > y_max {
            y_max = vertex.y();
        }
    }

    let scale_factor = (x_max - x_min).max(y_max - y_min);

    for vertex in vertices.iter() {
        normalized_vertices.push(Vertex {
            x: (vertex.x() - x_min) / scale_factor,
            y: (vertex.y() - y_min) / scale_factor,
        });
    }

    (normalized_vertices, scale_factor, x_min, y_min)
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub(crate) enum VertexPlacement {
    InsideTriangle(TriangleId),
    OnTriangleEdge(TriangleId, TriangleEdgeIndex, Edge),
    NotFound,
}

fn find_vertex_placement(
    vertex: Vertex,
    from: TriangleId,
    triangles: &Triangles,
    vertices: &[Vertex],
    #[cfg(feature = "debug_context")] _debug_context: &mut DebugContext,
) -> VertexPlacement {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "find_vertex_placement").entered();

    let mut current_triangle: Neighbor = from.into();
    let mut previous_triangle = current_triangle;

    // We use `triangles.len()` as an upper bound on the number of triangles
    for _index in 0..triangles.count() {
        if !current_triangle.exists() {
            break;
        }

        let triangle_id = current_triangle.id;
        let triangle = triangles.get(triangle_id);

        let infinite_verts = collect_infinite_triangle_vertices(&triangle.verts);

        // This is the hot branch
        if infinite_verts.is_empty() {
            let (v1, v2, v3) = triangle.to_vertices(vertices);

            // TODO Optim: We could in theory avoid checking 1 edge. See the "previous_edge" branch. In practice, it seems to have too many cases disjunctions (initial search versus general iteration).
            // We don't need to check the edge we come from
            // These 3 checks are the hot ones. Most iterations of the loop will go from one triangle to another here.
            // This is why we don't do the side tests until we need them. (As a side effect, for the last triangle, we will do 2 of those tests twice)
            if previous_triangle != triangle.neighbor12() {
                let edge_12_test = test_point_edge_side((v1, v2), vertex);
                if edge_12_test.is_on_left_side() {
                    previous_triangle = current_triangle;
                    current_triangle = triangle.neighbor12();
                    continue;
                }
            }
            if previous_triangle != triangle.neighbor23() {
                let edge_23_test = test_point_edge_side((v2, v3), vertex);
                if edge_23_test.is_on_left_side() {
                    previous_triangle = current_triangle;
                    current_triangle = triangle.neighbor23();
                    continue;
                }
            }
            if previous_triangle != triangle.neighbor31() {
                let edge_31_test = test_point_edge_side((v3, v1), vertex);
                if edge_31_test.is_on_left_side() {
                    previous_triangle = current_triangle;
                    current_triangle = triangle.neighbor31();
                    continue;
                }
            }

            let edge_test = test_point_edge_side((v1, v2), vertex);
            if edge_test.is_colinear() {
                return VertexPlacement::OnTriangleEdge(triangle_id, EDGE_12, triangle.edge12());
            }
            let edge_test = test_point_edge_side((v2, v3), vertex);
            if edge_test.is_colinear() {
                return VertexPlacement::OnTriangleEdge(triangle_id, EDGE_23, triangle.edge23());
            }
            let edge_test = test_point_edge_side((v3, v1), vertex);
            if edge_test.is_colinear() {
                return VertexPlacement::OnTriangleEdge(triangle_id, EDGE_31, triangle.edge31());
            }

            return VertexPlacement::InsideTriangle(triangle_id);
        } else if infinite_verts.len() == 1 {
            let res = vertex_placement_1_infinite_vertex(
                vertices,
                vertex,
                triangle,
                triangle_id,
                infinite_verts[0],
                &mut previous_triangle,
                &mut current_triangle,
            );
            match res {
                Some(placement) => return placement,
                None => continue,
            }
        } else {
            let res = vertex_placement_2_infinite_vertex(
                vertices,
                vertex,
                triangle,
                triangle_id,
                infinite_verts,
                &mut previous_triangle,
                &mut current_triangle,
            );
            match res {
                Some(placement) => return placement,
                None => continue,
            }
        }
    }

    VertexPlacement::NotFound
}

/// Splits the infinite quad into 4 triangles at p, the first inserted vertex.
///
/// ```text
///          1
///          |
///     t1   |    t2
///          |
/// 0--------p---------2
///          |
///     t4   |    t3
///          |
///          3
/// ```
fn create_initial_triangles(vertices_count: usize, first_vertex_id: VertexId) -> Triangles {
    let mut triangles = Triangles::with_capacity(vertices_count * 2 + 2);

    let (t1, t2, t3, t4) = (0.into(), 1.into(), 2.into(), 3.into());
    triangles.buffer_mut().push(TriangleData {
        verts: [first_vertex_id, INFINITE_V0_ID, INFINITE_V1_ID],
        neighbors: [t4, Neighbor::NONE, t2],
    });
    triangles.buffer_mut().push(TriangleData {
        verts: [first_vertex_id, INFINITE_V1_ID, INFINITE_V2_ID],
        neighbors: [t1, Neighbor::NONE, t3],
    });
    triangles.buffer_mut().push(TriangleData {
        verts: [first_vertex_id, INFINITE_V2_ID, INFINITE_V3_ID],
        neighbors: [t2, Neighbor::NONE, t4],
    });
    triangles.buffer_mut().push(TriangleData {
        verts: [first_vertex_id, INFINITE_V3_ID, INFINITE_V0_ID],
        neighbors: [t3, Neighbor::NONE, t1],
    });
    triangles
}

#[inline]
pub(crate) fn is_vertex_pair_too_close(a: Vertex, b: Vertex) -> bool {
    let dist = a - b;
    dist.x.abs() < Float::EPSILON && dist.y.abs() < Float::EPSILON
}

fn verts_count_to_restoration_stack_initial_capacity(verts_count: usize) -> usize {
    if verts_count < 9_000 {
        16
    } else if verts_count < 150_000 {
        32
    } else if verts_count < 500_000 {
        64
    } else {
        128
    }
}

/// - `vertices` should be normalized with their cooridnates in [0,1]
pub(crate) fn wrap_and_triangulate_2d_normalized_vertices(
    vertices: &[Vertex],
    vertex_merge_mapping: &mut Option<Vec<VertexId>>,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) -> Result<Triangles, TriangulationError> {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "wrap_and_triangulate_2d_normalized_vertices").entered();

    let ordered_vertices = insertion_order(vertices);

    let mut vertices_iterator = ordered_vertices.iter().enumerate();
    // There is always at least 3 vertices
    let (_, first_vertex_id) = vertices_iterator.next().unwrap();
    let mut triangles = create_initial_triangles(vertices.len(), *first_vertex_id);
    // Id of the triangle we are looking at
    let mut triangle_id = triangles.last_id();

    #[cfg(feature = "debug_context")]
    debug_context.push_snapshot(Phase::FirstVertexInsertion, &triangles, &[0], &[]);

    // This buffer is used in the loop. Filled and cleared each time.
    // We create if here to share the allocation between all those calls as an optimization.
    let mut quads_to_check = Vec::with_capacity(verts_count_to_restoration_stack_initial_capacity(
        vertices.len(),
    ));

    // Loop over all the input vertices
    for (_index, &vertex_id) in vertices_iterator {
        #[cfg(feature = "debug_context")]
        {
            let force_end = debug_context.advance_step();
            if force_end {
                break;
            }
        }

        // Find an existing triangle which encloses P
        let vertex = vertices[vertex_id as usize];
        let vertex_place = find_vertex_placement(
            vertex,
            triangle_id,
            &triangles,
            &vertices,
            #[cfg(feature = "debug_context")]
            debug_context,
        );
        match vertex_place {
            // Cases sorted by frequency
            VertexPlacement::InsideTriangle(enclosing_triangle_id) => {
                let close_vertex =
                    triangles
                        .get(enclosing_triangle_id)
                        .verts
                        .into_iter()
                        .find(|v| {
                            is_finite(*v) && is_vertex_pair_too_close(vertices[*v as usize], vertex)
                        });
                if let Some(close_vertex) = close_vertex {
                    if let Some(vertex_merge_mapping) = vertex_merge_mapping {
                        vertex_merge_mapping[vertex_id as usize] = close_vertex;
                    }
                    continue;
                }
                split_triangle_into_three_triangles(
                    &mut triangles,
                    enclosing_triangle_id,
                    vertex_id,
                    &mut quads_to_check,
                    #[cfg(feature = "debug_context")]
                    debug_context,
                )
            }
            VertexPlacement::OnTriangleEdge(enclosing_triangle_id, edge_index, edge) => {
                // Merge close vertices. We do it here, once, to avoid doing it multiple times during the vertex placement research.
                if is_finite(edge.from)
                    && is_vertex_pair_too_close(vertices[edge.from as usize], vertex)
                {
                    if let Some(vertex_merge_mapping) = vertex_merge_mapping {
                        vertex_merge_mapping[vertex_id as usize] = edge.from;
                    }
                    continue;
                }
                if is_finite(edge.to)
                    && is_vertex_pair_too_close(vertices[edge.to as usize], vertex)
                {
                    if let Some(vertex_merge_mapping) = vertex_merge_mapping {
                        vertex_merge_mapping[vertex_id as usize] = edge.to;
                    }
                    continue;
                }

                split_quad_into_four_triangles(
                    &mut triangles,
                    enclosing_triangle_id,
                    edge_index,
                    edge,
                    vertex_id,
                    &mut quads_to_check,
                    #[cfg(feature = "debug_context")]
                    debug_context,
                )
            }
            VertexPlacement::NotFound => {
                return Err(TriangulationError::new(
                    format!(
                        "Internal error, failed to find placement for vertex {:?}, step {}",
                        vertex_id, _index
                    ),
                    #[cfg(feature = "debug_context")]
                    debug_context,
                ));
            }
        };

        restore_delaunay_triangulation(
            &mut triangles,
            vertices,
            vertex_id,
            &mut quads_to_check,
            #[cfg(feature = "debug_context")]
            debug_context,
        );

        // We'll start the search for the next enclosing triangle from the last created triangle.
        // This is a pretty good heuristic since the vertices were spatially partitionned
        triangle_id = triangles.last_id();

        #[cfg(feature = "progress_log")]
        {
            if _index % ((ordered_vertices.len() / 50) + 1) == 0 {
                let progress = 100. * _index as f32 / ordered_vertices.len() as f32;
                info!(
                    "Triangulation progress, step n°{}, {}%: {}/{}",
                    debug_context.current_step,
                    progress,
                    _index,
                    ordered_vertices.len()
                );
            }
        }
    }

    Ok(triangles)
}

pub(crate) fn filter_triangles(
    triangles: &Triangles,
    config: &TriangulationConfiguration,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) -> Vec<[VertexId; 3]> {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "filter_triangles").entered();

    #[cfg(feature = "debug_context")]
    let mut filtered_debug_triangles = Triangles::new();

    #[cfg(feature = "parallel_filtering")]
    let indices = _multithreaded_filter_triangles(
        triangles,
        config,
        #[cfg(feature = "debug_context")]
        &mut filtered_debug_triangles,
    );
    #[cfg(not(feature = "parallel_filtering"))]
    let indices = _singlethreaded_filter_triangles(
        triangles,
        config.triangle_orientation,
        #[cfg(feature = "debug_context")]
        &mut filtered_debug_triangles,
    );

    #[cfg(feature = "debug_context")]
    debug_context.push_snapshot(Phase::FilterTriangles, &filtered_debug_triangles, &[], &[]);

    indices
}

#[cfg(feature = "parallel_filtering")]
pub(crate) fn _multithreaded_filter_triangles(
    triangles: &Triangles,
    config: &TriangulationConfiguration,
    #[cfg(feature = "debug_context")] filtered_debug_triangles: &mut Triangles,
) -> Vec<[VertexId; 3]> {
    let orientation = config.triangle_orientation;
    if triangles.count() > config.filter_parallel_tri_count_threshold {
        // Debug loop out of the // loop because filtered_debug_triangles cannot be accessed as mut from multiple tasks.
        #[cfg(feature = "debug_context")]
        {
            filtered_debug_triangles.buffer = triangles
                .buffer()
                .iter()
                .filter(|t| t.is_finite())
                .map(|t| t.clone())
                .collect();
        }

        triangles
            .buffer()
            .par_iter()
            .with_min_len(config.filter_parallel_min_batch_len)
            .filter_map(|t| match t.is_finite() {
                true => Some(orientation.orient(t.verts)),
                false => None,
            })
            // TODO Is it possible to pre-set the capacity of indices ?
            // .collect_into_vec(&mut indices);
            .collect()
    } else {
        _singlethreaded_filter_triangles(
            triangles,
            orientation,
            #[cfg(feature = "debug_context")]
            filtered_debug_triangles,
        )
    }
}

pub(crate) fn _singlethreaded_filter_triangles(
    triangles: &Triangles,
    orientation: TriangleOrientation,
    #[cfg(feature = "debug_context")] filtered_debug_triangles: &mut Triangles,
) -> Vec<[VertexId; 3]> {
    let mut indices = Vec::with_capacity(triangles.count());
    for t in triangles.buffer().iter() {
        if t.is_finite() {
            indices.push(orientation.orient(t.verts));
            #[cfg(feature = "debug_context")]
            filtered_debug_triangles.push(t.clone());
        }
    }
    indices
}

/// Order in which to insert the vertices: BRIO, from "Incremental constructions con BRIO" (Amenta, Choi, Rote, 2003).
///
/// Vertices go into random rounds of doubling size, and each round follows a Hilbert curve. Consecutive vertices are
/// then close to each other, so the search for their triangle is short. The random rounds avoid inserting long runs
/// of aligned vertices, which create thin triangles that later vertices keep flipping.
///
/// - `vertices` should be normalized with their coordinates in [0,1]
fn insertion_order(vertices: &[Vertex]) -> Vec<VertexId> {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "insertion_order").entered();

    // About 16 grid cells per vertex
    let order = (((vertices.len() as f64).log2() / 2.).ceil() as u32 + 2).clamp(1, 16);
    let max_cell = ((1u32 << order) - 1) as Float;
    let mut keys: Vec<(u32, VertexId)> = vertices
        .iter()
        .enumerate()
        .map(|(id, v)| {
            let x = (v.x * max_cell) as u32;
            let y = (v.y * max_cell) as u32;
            (hilbert_index(order, x, y), id as VertexId)
        })
        .collect();

    // Radix sort of the keys, which have 2 * order bits: two passes of `order` bits
    let buckets = 1usize << order;
    let mask = (buckets - 1) as u32;
    let mut sorted = vec![(0, 0); keys.len()];
    let mut counts = vec![0usize; buckets + 1];
    for shift in [0, order] {
        counts.fill(0);
        for key in &keys {
            counts[((key.0 >> shift) & mask) as usize + 1] += 1;
        }
        for bucket in 0..buckets {
            counts[bucket + 1] += counts[bucket];
        }
        for key in &keys {
            let bucket = ((key.0 >> shift) & mask) as usize;
            sorted[counts[bucket]] = *key;
            counts[bucket] += 1;
        }
        std::mem::swap(&mut keys, &mut sorted);
    }

    // Round LAST_ROUND - k holds about 1 / 2^(k+1) of the vertices. The seed is fixed so that results are reproducible.
    const LAST_ROUND: usize = 40;
    let mut state = 0x9E37_79B9_7F4A_7C15u64;
    let rounds: Vec<u8> = (0..vertices.len())
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (LAST_ROUND - (state.trailing_zeros() as usize).min(LAST_ROUND)) as u8
        })
        .collect();

    // Stable counting sort by round, which keeps the curve order within each round
    let mut round_starts = [0usize; LAST_ROUND + 2];
    for &(_, id) in &keys {
        round_starts[rounds[id as usize] as usize + 1] += 1;
    }
    for round in 0..=LAST_ROUND {
        round_starts[round + 1] += round_starts[round];
    }
    let mut ordered = vec![0; keys.len()];
    for &(_, id) in &keys {
        let round = rounds[id as usize] as usize;
        ordered[round_starts[round]] = id;
        round_starts[round] += 1;
    }
    ordered
}

/// Position of the cell (x, y) along a Hilbert curve over a 2^order x 2^order grid, with order <= 16.
///
/// Branchless version from <https://github.com/rawrunprotected/hilbert_curves> (public domain).
fn hilbert_index(order: u32, x: u32, y: u32) -> u32 {
    let x = x << (16 - order);
    let y = y << (16 - order);

    let (mut a_, mut b_, mut c_, mut d_);
    {
        let a = x ^ y;
        let b = 0xFFFF ^ a;
        let c = 0xFFFF ^ (x | y);
        let d = x & (y ^ 0xFFFF);
        a_ = a | (b >> 1);
        b_ = (a >> 1) ^ a;
        c_ = ((c >> 1) ^ (b & (d >> 1))) ^ c;
        d_ = ((a & (c >> 1)) ^ (d >> 1)) ^ d;
    }
    for shift in [2, 4] {
        let (a, b, c, d) = (a_, b_, c_, d_);
        a_ = (a & (a >> shift)) ^ (b & (b >> shift));
        b_ = (a & (b >> shift)) ^ (b & ((a ^ b) >> shift));
        c_ ^= (a & (c >> shift)) ^ (b & (d >> shift));
        d_ ^= (b & (c >> shift)) ^ ((a ^ b) & (d >> shift));
    }
    {
        let (a, b, c, d) = (a_, b_, c_, d_);
        c_ ^= (a & (c >> 8)) ^ (b & (d >> 8));
        d_ ^= (b & (c >> 8)) ^ ((a ^ b) & (d >> 8));
    }
    let a = c_ ^ (c_ >> 1);
    let b = d_ ^ (d_ >> 1);
    let i0 = x ^ y;
    let i1 = b | (0xFFFF ^ (i0 | a));
    ((interleave(i1) << 1) | interleave(i0)) >> (32 - 2 * order)
}

/// Spreads the 16 low bits of `x` to the even bits
fn interleave(x: u32) -> u32 {
    let x = (x | (x << 8)) & 0x00FF00FF;
    let x = (x | (x << 4)) & 0x0F0F0F0F;
    let x = (x | (x << 2)) & 0x33333333;
    (x | (x << 1)) & 0x55555555
}

/// Splits a quad (triangle pair) into 4 triangles
///```text
///      from
///      /|\
///     / | \
/// n2 /  |  \ n3
///   /   |   \
///  /    |    \
/// /  t1 | t2  \
/// \     |     /
///  \    |    /
///   \   |   /
/// n1 \  |  / n4
///     \ | /
///      \|/
///       to
///
/// Becomes
///
///      /|\
///     /3|2\
/// n2 /  |  \ n3
///   / t1|t2 \
///  /    |    \
/// /2___1|1___3\
/// \3   1|1   2/
///  \ t3 | t4 /
///   \   |   /
/// n1 \  |  / n4
///     \2|3/
///      \|/
/// ```
pub(crate) fn split_quad_into_four_triangles(
    triangles: &mut Triangles,
    triangle_id: TriangleId,
    edge_index: TriangleEdgeIndex,
    edge: Edge,
    vertex_id: VertexId,
    quads_to_check: &mut Vec<(TriangleId, TriangleId)>,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "split_quad_into_four_triangles").entered();

    // Re-use the existing triangle id for the first triangle
    let t1 = triangle_id;
    // TODO Comment is not true anymore in the general sense, but still true in this function
    // The triangle is always guaranteed to have a neighbor on the edge since a vertex
    // cannot land on a container triangle's edge.
    let t2 = triangles.get(triangle_id).neighbor(edge_index).id;
    // Create two new triangles for the other two
    let t3 = triangles.next_id();
    let t4 = triangles.next_id() + 1;

    let t1_before = triangles.get(t1).clone();
    let t2_before = triangles.get(t2).clone();

    let n1 = t1_before.neighbor(next_clockwise_edge_index(edge_index));
    let n2 = t1_before.neighbor(next_counter_clockwise_edge_index(edge_index));

    // t3
    let t3_v3 = t1_before.v(opposite_vertex_index(edge_index));
    triangles.create([vertex_id, edge.to, t3_v3], [t4.into(), n1, t1.into()]);

    // TODO Optim: Can do a tiny bit better
    let t2_opposite_v_id = t2_before.get_opposite_vertex_id(&edge);
    let t2_opposite_v_index = t2_before.vertex_index(t2_opposite_v_id);

    let n3 = t2_before.neighbor(vertex_next_ccw_edge_index(t2_opposite_v_index));
    let n4 = t2_before.neighbor(vertex_next_cw_edge_index(t2_opposite_v_index));

    // t4
    triangles.create(
        [vertex_id, t2_opposite_v_id, edge.to],
        [t2.into(), n4, t3.into()],
    );

    // Update t1 verts
    let verts = [vertex_id, t3_v3, edge.from];
    triangles.get_mut(t1).verts = verts;
    // Update t1 neighbors
    let neighbors = [t3.into(), n2, t2.into()];
    triangles.get_mut(t1).neighbors = neighbors;

    // Update t2 verts
    let verts = [vertex_id, edge.from, t2_opposite_v_id];
    triangles.get_mut(t2).verts = verts;
    // Update t2 neighbors
    let neighbors = [t1.into(), n3, t4.into()];
    triangles.get_mut(t2).neighbors = neighbors;

    // Update neighbors and fill in quads_to_check directly for later. Way better for performances.
    if n1.exists() {
        update_triangle_neighbor(triangles.get_mut(n1.id), t1.into(), t3.into());
        quads_to_check.push((t3, n1.id));
    }
    if n4.exists() {
        update_triangle_neighbor(triangles.get_mut(n4.id), t2.into(), t4.into());
        quads_to_check.push((t4, n4.id));
    }
    if n2.exists() {
        quads_to_check.push((t1, n2.id));
    }
    if n3.exists() {
        quads_to_check.push((t2, n3.id));
    }

    #[cfg(feature = "debug_context")]
    debug_context.push_snapshot_event(
        Phase::SplitQuad,
        EventInfo::Split(vertex_id),
        &triangles,
        &[t1, t2, t3, t4],
        &[],
    );
}

/// Splits `triangle_id` into 3 triangles (re-using the existing triangle id)
///
/// All the resulting triangles will share `vertex_id` as their first vertex, and will be oriented in a CW order
///
/// ```text
///                  v1
///                / | \
///               / 3|2 \
///              /   |   \
///             /    |    \
///       n6   / t1  |  t3 \   n4
///           /     1|1     \
///          /      /1\      \
///         /     /     \     \
///        /    /         \    \
///       /   /             \   \
///      /2 /        t2       \ 3\
///     / / 3                 2 \ \
///   v3 ------------------------- v2
///                  n5
/// ```
pub(crate) fn split_triangle_into_three_triangles(
    triangles: &mut Triangles,
    triangle_id: TriangleId,
    vertex_id: VertexId,
    quads_to_check: &mut Vec<(TriangleId, TriangleId)>,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "split_triangle_into_three_triangles").entered();

    // Re-use the existing triangle id for the first triangle
    let t1 = triangle_id;
    // Create two new triangles for the other two
    let t2 = triangles.next_id();
    let t3 = triangles.next_id() + 1;
    let n4 = triangles.get(t1).neighbor12();
    let n5 = triangles.get(t1).neighbor23();
    let n6 = triangles.get(t1).neighbor31();

    // t2
    triangles.create(
        [vertex_id, triangles.get(t1).v2(), triangles.get(t1).v3()],
        [t3.into(), n5, t1.into()],
    );
    // t3
    triangles.create(
        [vertex_id, triangles.get(t1).v1(), triangles.get(t1).v2()],
        [t1.into(), n4, t2.into()],
    );

    // Update t1 verts
    let verts = [vertex_id, triangles.get(t1).v3(), triangles.get(t1).v1()];
    triangles.get_mut(t1).verts = verts;
    // Update t1 neighbors
    let neighbors = [t2.into(), n6, t3.into()];
    triangles.get_mut(t1).neighbors = neighbors;

    // Update neighbors and fill in quads_to_check directly for later. Way better for performances.
    if n4.exists() {
        update_triangle_neighbor(triangles.get_mut(n4.id), t1.into(), t3.into());
        quads_to_check.push((t3, n4.id));
    }
    if n5.exists() {
        update_triangle_neighbor(triangles.get_mut(n5.id), t1.into(), t2.into());
        quads_to_check.push((t2, n5.id));
    }
    if n6.exists() {
        quads_to_check.push((t1, n6.id));
    }

    #[cfg(feature = "debug_context")]
    debug_context.push_snapshot_event(
        Phase::SplitTriangle,
        EventInfo::Split(vertex_id),
        &triangles,
        &[t1, t2, t3],
        // Neighbor data is out of date here and is not that interesting
        &[],
    );
}

pub(crate) fn update_neighbor_neighbor(
    triangle: Neighbor,
    old_neighbour_id: Neighbor,
    new_neighbour_id: Neighbor,
    triangles: &mut Triangles,
) {
    #[cfg(feature = "more_profile_traces")]
    let _span = span!(Level::TRACE, "update_neighbor_neighbor").entered();

    if triangle.exists() {
        let t = triangles.get_mut(triangle.id);
        if t.neighbor12() == old_neighbour_id {
            *t.neighbor12_mut() = new_neighbour_id;
        } else if t.neighbor23() == old_neighbour_id {
            *t.neighbor23_mut() = new_neighbour_id;
        } else if t.neighbor31() == old_neighbour_id {
            *t.neighbor31_mut() = new_neighbour_id;
        }
    }
}

pub(crate) fn update_triangle_neighbor(
    t: &mut TriangleData,
    old_neighbour_id: Neighbor,
    new_neighbour_id: Neighbor,
) {
    #[cfg(feature = "more_profile_traces")]
    let _span = span!(Level::TRACE, "update_triangle_neighbor").entered();

    if t.neighbor12() == old_neighbour_id {
        *t.neighbor12_mut() = new_neighbour_id;
    } else if t.neighbor23() == old_neighbour_id {
        *t.neighbor23_mut() = new_neighbour_id;
    } else if t.neighbor31() == old_neighbour_id {
        *t.neighbor31_mut() = new_neighbour_id;
    }
}

/// `quads_to_check` used the shared pre-allocated buffer for efficiency.
/// - It does not need to be cleared since it is fully emptied by each call to restore_delaunay_triangulation
fn restore_delaunay_triangulation(
    triangles: &mut Triangles,
    vertices: &[Vertex],
    from_vertex_id: VertexId,
    quads_to_check: &mut Vec<(TriangleId, TriangleId)>,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "restore_delaunay_triangulation").entered();

    while let Some((from_triangle_id, opposite_triangle_id)) = quads_to_check.pop() {
        check_and_swap_quad_diagonal(
            triangles,
            vertices,
            from_vertex_id,
            from_triangle_id,
            opposite_triangle_id,
            quads_to_check,
            #[cfg(feature = "debug_context")]
            debug_context,
        );
    }
}

#[inline(always)]
pub(crate) fn should_swap_diagonals(quad: &Quad, vertices: &[Vertex]) -> bool {
    #[cfg(feature = "more_profile_traces")]
    let _span = span!(Level::TRACE, "should_swap_diagonals").entered();

    // q4 is the inserted vertex, so it is never infinite here.
    // An infinite q3 means q1 q2 q3 has no circumcircle: it covers a strip or a corner of the plane, which
    // is_vertex_beyond_infinite_corner tests instead.
    if is_infinite(quad.v3()) {
        return is_vertex_beyond_infinite_corner(vertices, quad);
    }

    let is_q1_finite = is_finite(quad.v1());
    let is_q2_finite = is_finite(quad.v2());
    if is_q1_finite && is_q2_finite {
        is_vertex_in_triangle_circumcircle(
            vertices[quad.v1() as usize],
            vertices[quad.v2() as usize],
            vertices[quad.v3() as usize],
            vertices[quad.v4() as usize],
        )
    } else if !is_q1_finite && !is_q2_finite {
        true
    } else if !is_q1_finite {
        is_vertex_in_half_plane_1(vertices, quad, VERT_1)
    } else {
        is_vertex_in_half_plane_1(vertices, quad, VERT_2)
    }
}

/// ```text
///                q3
///         n3   /    \   n4
///            /   t2   \
///          /            \
///         q1 ---------- q2
///          \ 2        3 /
///            \   t1   /
///         n2   \ 1  /   n1
///                q4
/// ```
///
/// If q4 is in the circumcircle of the tirangle q1q2q3, the two triangles form a convex quadrilateral
/// whose diagonal is drawn in the wrong direction. We swap this diagonal to form two new triangles so
/// that the structure of the Delaunay triangulation is locally restored.
///
/// The quad becomes
///
/// ```text
///               q3
///        n3   / 3|2 \   n4
///           /    |    \
///         /      |      \
///       q1 2  t1 | t2  3 q2
///         \      |      /
///           \    |    /
///        n2   \ 1|1 /   n1
///               q4
/// ```
pub(crate) fn check_and_swap_quad_diagonal(
    triangles: &mut Triangles,
    vertices: &[Vertex],
    from_vertex_id: VertexId,
    t1_id: TriangleId,
    t2_id: TriangleId,
    quads_to_check: &mut Vec<(TriangleId, TriangleId)>,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) {
    #[cfg(feature = "more_profile_traces")]
    let _span = span!(Level::TRACE, "check_and_swap_quad_diagonal").entered();

    let t2 = triangles.get(t2_id);

    let (quad, n3, n4) =
    // No need to check if neighbor exists, handled by the == check since `t1_id` exists
        if t2.neighbor12().id == t1_id {
            (
                Quad::new([
                    t2.v2(),
                    t2.v1(),
                    t2.v3(),
                    from_vertex_id,
                ]),
                t2.neighbor23(),
                t2.neighbor31(),
            )
        } else if t2.neighbor23().id == t1_id {
            (
                Quad::new([
                    t2.v3(),
                    t2.v2(),
                    t2.v1(),
                    from_vertex_id,
                ]),
                t2.neighbor31(),
                t2.neighbor12(),
            )
        } else {
            (
                Quad::new([
                    t2.v1(),
                    t2.v3(),
                    t2.v2(),
                    from_vertex_id,
                ]),
                t2.neighbor12(),
                t2.neighbor23(),
            )
        };

    if should_swap_diagonals(&quad, vertices) {
        let t2_as_neighbor = t2_id.into();
        let t1_as_neighbor = t1_id.into();

        let n1 = triangles.get(t1_id).neighbor31();
        let n2 = triangles.get(t1_id).neighbor12();

        // Tmp
        *triangles.get_mut(t1_id) = TriangleData {
            verts: [quad.v4(), quad.v1(), quad.v3()],
            neighbors: [n2, n3, t2_as_neighbor],
        };
        *triangles.get_mut(t2_id) = TriangleData {
            verts: [quad.v4(), quad.v3(), quad.v2()],
            neighbors: [t1_as_neighbor, n4, n1],
        };

        // Update neighbors and place any new triangles pairs which are now opposite to `from_vertex_id` on the stack, to be checked
        if n1.exists() {
            update_triangle_neighbor(triangles.get_mut(n1.id), t1_as_neighbor, t2_as_neighbor);
        }
        if n3.exists() {
            update_triangle_neighbor(triangles.get_mut(n3.id), t2_as_neighbor, t1_as_neighbor);
            quads_to_check.push((t1_id, n3.id));
        }
        if n4.exists() {
            quads_to_check.push((t2_id, n4.id));
        }

        #[cfg(feature = "debug_context")]
        debug_context.push_snapshot(
            Phase::DelaunayRestoreSwapQuadDiagonals,
            &triangles,
            &[t1_id, t2_id],
            &[n3, n1],
        );
    }
}

///////////////////////////////////////////////////////////
///                                                     ///
///                        Tests                        ///
///                                                     ///
///////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use glam::Vec3;

    use crate::{
        triangulation::{
            hilbert_index, insertion_order, normalize_vertices_coordinates,
            transform_to_2d_planar_coordinate_system,
        },
        types::{Float, Vertex, VertexId},
    };

    #[test]
    fn hilbert_curve_visits_each_cell_once_through_adjacent_cells() {
        for order in [1, 2, 6] {
            let side = 1u32 << order;
            let mut cells: Vec<(u32, u32, u32)> = (0..side)
                .flat_map(|x| (0..side).map(move |y| (hilbert_index(order, x, y), x, y)))
                .collect();
            cells.sort();
            assert!(cells.iter().enumerate().all(|(i, cell)| cell.0 == i as u32));
            for pair in cells.windows(2) {
                assert_eq!(
                    1,
                    pair[0].1.abs_diff(pair[1].1) + pair[0].2.abs_diff(pair[1].2)
                );
            }
        }
    }

    #[test]
    fn insertion_order_has_each_vertex_once() {
        let vertices: Vec<Vertex> = (0..1000)
            .map(|i| {
                Vertex::new(
                    ((i * 37) % 1000) as Float / 1000.,
                    ((i * 91) % 1000) as Float / 1000.,
                )
            })
            .collect();
        let mut order = insertion_order(&vertices);
        order.sort_unstable();
        assert_eq!((0..1000).collect::<Vec<VertexId>>(), order);
    }

    #[test]
    fn triangulation_normalize_set_of_vertices() {
        let vertices = vec![
            Vertex::new(3., 2.),
            Vertex::new(-1., 2.),
            Vertex::new(-1., -2.),
            Vertex::new(3., -2.),
        ];

        let normalized_vertices = normalize_vertices_coordinates(&vertices);

        assert_eq!(
            vec![
                Vertex::new(1., 1.),
                Vertex::new(0., 1.),
                Vertex::new(0., 0.),
                Vertex::new(1., 0.),
            ],
            normalized_vertices.0
        );
    }

    #[test]
    fn triangulation_set_to_2d_plane_vertices() {
        let vertices = vec![
            Vec3::new(-3., 2., 0.),
            Vec3::new(1., 2., 0.),
            Vec3::new(11., -2., 0.),
            Vec3::new(-3., -2., 0.),
        ];

        let planar_vertices = transform_to_2d_planar_coordinate_system(&vertices, Vec3::Z);

        assert_eq!(
            vec![
                -Vertex::new(-3., 2.),
                -Vertex::new(1., 2.),
                -Vertex::new(11., -2.),
                -Vertex::new(-3., -2.),
            ],
            planar_vertices
        );
    }
}
