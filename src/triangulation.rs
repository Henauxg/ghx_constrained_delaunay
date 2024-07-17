use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};
use tracing::error;

use crate::types::{
    infinite_vertex_local_quad_index, is_infinite, next_clockwise_edge_index,
    next_clockwise_vertex_index, next_counter_clockwise_edge_index,
    next_counter_clockwise_vertex_index, opposite_edge_index, opposite_vertex_index,
    opposite_vertex_index_from_edge, vertex_next_ccw_edge_index, vertex_next_cw_edge_index,
    EdgeVertices, Float, Neighbor, Quad, QuadVertexIndex, TriangleData, TriangleEdgeIndex,
    TriangleId, TriangleVertexIndex, Triangles, Vector3, Vector3A, Vertex, VertexId, EDGE_12,
    EDGE_23, EDGE_31, EDGE_TO_VERTS, INFINITE_V0_ID, INFINITE_V1_ID, INFINITE_V2_ID,
    INFINITE_V3_ID, NEXT_CCW_VERTEX_INDEX, NEXT_CW_VERTEX_INDEX, QUAD_1, QUAD_2, QUAD_3, QUAD_4,
    VERT_1, VERT_2, VERT_3,
};
use crate::utils::{
    is_point_strictly_on_right_side_of_edge, is_vertex_in_triangle_circumcircle,
    test_point_edge_side,
};

#[cfg(feature = "progress_log")]
use tracing::info;

#[cfg(feature = "debug_context")]
use crate::debug::{DebugConfiguration, DebugContext, EventInfo, Phase};

#[cfg(feature = "profile_traces")]
use tracing::{span, Level};

/// Binsort will cover the region to be triangulated by a rectangular grid so that each bin contains roughly N^(density_power) points.
pub const DEFAULT_BIN_VERTEX_DENSITY_POWER: f64 = 0.5;

pub const DEFAULT_ILTER_PARALLEL_TRI_COUNT_THRESHOLD: usize = 100_000;
pub const DEFAULT_FILTER_PARALLEL_MIN_BATCH_LEN: usize = 1000;

#[derive(Clone, Debug)]
pub struct TriangulationConfiguration {
    /// Binsort will cover the region to be triangulated by a rectangular grid so that each bin contains roughly N^(density_power) points.
    pub bin_vertex_density_power: f64,
    // TODO Doc
    pub filter_parallel_tri_count_threshold: usize,
    pub filter_parallel_min_batch_len: usize,
    #[cfg(feature = "debug_context")]
    pub debug_config: DebugConfiguration,
}
impl Default for TriangulationConfiguration {
    fn default() -> Self {
        Self {
            bin_vertex_density_power: DEFAULT_BIN_VERTEX_DENSITY_POWER,
            filter_parallel_tri_count_threshold: DEFAULT_ILTER_PARALLEL_TRI_COUNT_THRESHOLD,
            filter_parallel_min_batch_len: DEFAULT_FILTER_PARALLEL_MIN_BATCH_LEN,
            #[cfg(feature = "debug_context")]
            debug_config: DebugConfiguration::default(),
        }
    }
}

/// TODO A doc comment about infinite vertices positions

// Slightly higher values than 1.0 to be safe, which is enough since all our vertices coordinates are normalized inside the unit square.
pub const DELTA_VALUE: Float = 1.42;

/// plane_normal must be normalized
/// vertices must all belong to a 3d plane
pub fn triangulation_from_3d_planar_vertices(
    vertices: &Vec<Vector3>,
    plane_normal: Vector3A,
    config: TriangulationConfiguration,
) -> Triangulation {
    if vertices.len() < 3 {
        return Triangulation::default();
    }
    let mut planar_vertices = transform_to_2d_planar_coordinate_system(vertices, plane_normal);
    triangulation_from_2d_vertices(&mut planar_vertices, config)
}

#[derive(Default)]
pub struct Triangulation {
    /// Indices of the original vertices by groups of 3 to from triangles.
    pub triangles: Vec<[VertexId; 3]>,

    #[cfg(feature = "debug_context")]
    pub debug_context: DebugContext,
}

pub fn triangulation_from_2d_vertices(
    vertices: &Vec<Vertex>,
    config: TriangulationConfiguration,
) -> Triangulation {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "triangulation_from_2d_vertices").entered();

    if vertices.len() < 3 {
        return Triangulation::default();
    }

    // Uniformly scale the coordinates of the points so that they all lie between 0 and 1.
    let (mut normalized_vertices, _scale_factor, _x_min, _y_min) =
        normalize_vertices_coordinates(vertices);

    #[cfg(feature = "debug_context")]
    let mut debug_context =
        DebugContext::new(config.debug_config.clone(), _scale_factor, _x_min, _y_min);

    let triangles = wrap_and_triangulate_2d_normalized_vertices(
        &mut normalized_vertices,
        config.bin_vertex_density_power,
        &mut None,
        #[cfg(feature = "debug_context")]
        &mut debug_context,
    );

    let vert_indices = remove_wrapping(
        &triangles,
        &config,
        #[cfg(feature = "debug_context")]
        &mut debug_context,
    );

    Triangulation {
        triangles: vert_indices,
        #[cfg(feature = "debug_context")]
        debug_context,
    }
}

/// Transforms 3d coordinates of all vertices into 2d coordinates on a plane defined by the given normal and vertices.
/// - Input vertices need to all belong to the same 3d plan
/// - There must be at least two vertices
pub fn transform_to_2d_planar_coordinate_system(
    vertices: &Vec<Vector3>,
    plane_normal: Vector3A,
) -> Vec<Vertex> {
    // Create a base B for the plan, using the first two vertices as the first base vector
    let basis_1 = (vertices[0] - vertices[1]).normalize();
    let basis_2 = basis_1.cross((-plane_normal).into());
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
pub(crate) fn normalize_vertices_coordinates(
    vertices: &Vec<Vertex>,
) -> (Vec<Vertex>, Float, Float, Float) {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "normalize_vertices_coordinates").entered();

    let mut normalized_vertices = Vec::with_capacity(vertices.len());
    let (mut x_min, mut y_min, mut x_max, mut y_max) =
        (Float::MAX, Float::MAX, Float::MIN, Float::MIN);

    for vertex in vertices.iter() {
        if vertex.x < x_min {
            x_min = vertex.x;
        }
        if vertex.x > x_max {
            x_max = vertex.x;
        }
        if vertex.y < y_min {
            y_min = vertex.y;
        }
        if vertex.y > y_max {
            y_max = vertex.y;
        }
    }

    let scale_factor = (x_max - x_min).max(y_max - y_min);

    for vertex in vertices.iter() {
        normalized_vertices.push(Vertex {
            x: (vertex.x - x_min) / scale_factor,
            y: (vertex.y - y_min) / scale_factor,
        });
    }

    (normalized_vertices, scale_factor, x_min, y_min)
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
enum VertexPlacement {
    InsideTriangle(TriangleId),
    OnTriangleEdge(TriangleId, TriangleEdgeIndex),
    OnVertex(VertexId),
}

// TODO Details
pub struct TriangulationError;

// TODO Clean: this function should be shorter
fn find_vertex_placement(
    vertex: Vertex,
    from: TriangleId,
    triangles: &Triangles,
    vertices: &Vec<Vertex>,
    #[cfg(feature = "debug_context")] _debug_context: &mut DebugContext,
) -> Result<VertexPlacement, TriangulationError> {
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

        let mut infinite_verts = Vec::new();
        if is_infinite(triangle.v1()) {
            infinite_verts.push(VERT_1);
        }
        if is_infinite(triangle.v2()) {
            infinite_verts.push(VERT_2);
        }
        if is_infinite(triangle.v3()) {
            infinite_verts.push(VERT_3);
        }

        if infinite_verts.len() == 2 {
            let finite_v_index =
                opposite_vertex_index_from_edge(infinite_verts[0], infinite_verts[1]);
            let finite_vert_id = triangle.v(finite_v_index);
            let finite_vertex = vertices[finite_vert_id as usize];

            // No need to check if vertex is too close to an infinite vertex here
            if is_vertex_pair_too_close(finite_vertex, vertex) {
                return Ok(VertexPlacement::OnVertex(finite_vert_id));
            }

            // Check if the point is more towards this triangle's neighbor
            // Avoid looping between two triangles
            let infinite_vert_a_local_index = infinite_vertex_local_quad_index(
                triangle.v(next_clockwise_vertex_index(finite_v_index)),
            );
            let edge_a = edge_from_semi_infinite_edge(finite_vertex, infinite_vert_a_local_index);
            let edge_a_index = vertex_next_cw_edge_index(finite_v_index);
            let edge_a_test = test_point_edge_side(edge_a, vertex);
            if edge_a_test.is_on_left_side() && previous_triangle != triangle.neighbor(edge_a_index)
            {
                previous_triangle = current_triangle;
                current_triangle = triangle.neighbor(edge_a_index);
                continue;
            }

            let infinite_vert_b_local_index = infinite_vertex_local_quad_index(
                triangle.v(next_counter_clockwise_vertex_index(finite_v_index)),
            );
            let edge_b = edge_from_semi_infinite_edge(finite_vertex, infinite_vert_b_local_index);
            let edge_b_index = vertex_next_ccw_edge_index(finite_v_index);
            let edge_b_test = test_point_edge_side(edge_b, vertex);
            // Test for right side because edge_b is oriented the wrong way. We could also reverse edge_b
            if edge_b_test.is_on_right_side()
                && previous_triangle != triangle.neighbor(edge_b_index)
            {
                previous_triangle = current_triangle;
                current_triangle = triangle.neighbor(edge_b_index);
                continue;
            }

            if edge_a_test.is_near_edge() {
                return Ok(VertexPlacement::OnTriangleEdge(triangle_id, edge_a_index));
            } else if edge_b_test.is_near_edge() {
                return Ok(VertexPlacement::OnTriangleEdge(triangle_id, edge_b_index));
            } else {
                return Ok(VertexPlacement::InsideTriangle(triangle_id));
            }
        } else if infinite_verts.len() == 1 {
            let finite_vert_a_index = NEXT_CW_VERTEX_INDEX[infinite_verts[0] as usize];
            let finite_vert_a_id = triangle.v(finite_vert_a_index);
            let finite_vertex_a = vertices[finite_vert_a_id as usize];

            if is_vertex_pair_too_close(finite_vertex_a, vertex) {
                return Ok(VertexPlacement::OnVertex(finite_vert_a_id));
            }

            let finite_vert_b_index = NEXT_CCW_VERTEX_INDEX[infinite_verts[0] as usize];
            let finite_vert_b_id = triangle.v(finite_vert_b_index);
            let finite_vertex_b = vertices[finite_vert_b_id as usize];

            if is_vertex_pair_too_close(finite_vertex_b, vertex) {
                return Ok(VertexPlacement::OnVertex(finite_vert_b_id));
            }

            // Check if the point is more towards this triangle's neighbor
            // Avoid looping between two triangles
            let edge_ab_test = test_point_edge_side((finite_vertex_a, finite_vertex_b), vertex);
            let edge_ab_index = vertex_next_cw_edge_index(finite_vert_a_index);
            if edge_ab_test.is_on_left_side()
                && previous_triangle != triangle.neighbor(edge_ab_index)
            {
                previous_triangle = current_triangle;
                current_triangle = triangle.neighbor(edge_ab_index);
                continue;
            }

            let infinite_vert_local_index =
                infinite_vertex_local_quad_index(triangle.v(infinite_verts[0]));
            let edge_a = edge_from_semi_infinite_edge(finite_vertex_a, infinite_vert_local_index);
            let edge_a_index = vertex_next_ccw_edge_index(finite_vert_a_index);
            let edge_a_test = test_point_edge_side(edge_a, vertex);
            // Test for right side because edge_a is oriented the wrong way. We could also reverse edge_a
            if edge_a_test.is_on_right_side()
                && previous_triangle != triangle.neighbor(edge_a_index)
            {
                previous_triangle = current_triangle;
                current_triangle = triangle.neighbor(edge_a_index);
                continue;
            }

            let edge_b = edge_from_semi_infinite_edge(finite_vertex_b, infinite_vert_local_index);
            let edge_b_index = vertex_next_cw_edge_index(finite_vert_b_index);
            let edge_b_test = test_point_edge_side(edge_b, vertex);
            if edge_b_test.is_on_left_side() && previous_triangle != triangle.neighbor(edge_b_index)
            {
                previous_triangle = current_triangle;
                current_triangle = triangle.neighbor(edge_b_index);
                continue;
            }

            if edge_ab_test.is_near_edge() {
                return Ok(VertexPlacement::OnTriangleEdge(triangle_id, edge_ab_index));
            } else if edge_a_test.is_near_edge() {
                return Ok(VertexPlacement::OnTriangleEdge(triangle_id, edge_a_index));
            } else if edge_b_test.is_near_edge() {
                return Ok(VertexPlacement::OnTriangleEdge(triangle_id, edge_b_index));
            } else {
                // Imposssible to be on super triangle edge, would be colinear
                return Ok(VertexPlacement::InsideTriangle(triangle_id));
            }
        } else {
            let (v1, v2, v3) = triangle.to_vertices(vertices);

            // TODO Profile this check's position
            if is_vertex_pair_too_close(v1, vertex) {
                return Ok(VertexPlacement::OnVertex(triangle.v1()));
            }
            if is_vertex_pair_too_close(v2, vertex) {
                return Ok(VertexPlacement::OnVertex(triangle.v2()));
            }
            if is_vertex_pair_too_close(v3, vertex) {
                return Ok(VertexPlacement::OnVertex(triangle.v3()));
            }

            // Check if the point is more towards this triangle's neighbor
            // Avoid looping between two triangles
            let edge_12_test = test_point_edge_side((v1, v2), vertex);
            if edge_12_test.is_on_left_side() && previous_triangle != triangle.neighbor12() {
                previous_triangle = current_triangle;
                current_triangle = triangle.neighbor12();
                continue;
            }
            let edge_23_test = test_point_edge_side((v2, v3), vertex);
            if edge_23_test.is_on_left_side() && previous_triangle != triangle.neighbor23() {
                previous_triangle = current_triangle;
                current_triangle = triangle.neighbor23();
                continue;
            }
            let edge_31_test = test_point_edge_side((v3, v1), vertex);
            if edge_31_test.is_on_left_side() && previous_triangle != triangle.neighbor31() {
                previous_triangle = current_triangle;
                current_triangle = triangle.neighbor31();
                continue;
            }

            if edge_12_test.is_near_edge() {
                return Ok(VertexPlacement::OnTriangleEdge(triangle_id, EDGE_12));
            } else if edge_23_test.is_near_edge() {
                return Ok(VertexPlacement::OnTriangleEdge(triangle_id, EDGE_23));
            } else if edge_31_test.is_near_edge() {
                return Ok(VertexPlacement::OnTriangleEdge(triangle_id, EDGE_31));
            } else {
                return Ok(VertexPlacement::InsideTriangle(triangle_id));
            }
        }
    }

    Err(TriangulationError)
}

pub const INFINITE_VERTS_Y_DELTAS: [Float; 4] = [0., DELTA_VALUE, 0., -DELTA_VALUE];
pub const INFINITE_VERTS_X_DELTAS: [Float; 4] = [-DELTA_VALUE, 0., DELTA_VALUE, 0.];

/// Returns a finite segment from an edge between a finite vertex and an infinite vertex
///  - infinite_vert_local_index is the local index of the infinite vertex
#[inline]
pub fn edge_from_semi_infinite_edge(
    finite_vertex: Vertex,
    infinite_vert_local_index: QuadVertexIndex,
) -> EdgeVertices {
    (
        finite_vertex,
        Vertex::new(
            // We care about the delta sign here since we create a segment from an infinite line,
            // starting from the finite point and aimed towards the infinite point.
            finite_vertex.x + INFINITE_VERTS_X_DELTAS[infinite_vert_local_index as usize],
            finite_vertex.y + INFINITE_VERTS_Y_DELTAS[infinite_vert_local_index as usize],
        ),
    )
}

/// TODO Ascii art
fn create_initial_triangles(vertices_count: usize, first_vertex_id: u32) -> Triangles {
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

/// - `vertices` should be normalized with their cooridnates in [0,1]
pub(crate) fn wrap_and_triangulate_2d_normalized_vertices(
    vertices: &mut Vec<Vertex>,
    bin_vertex_density_power: f64,
    vertex_merge_mapping: &mut Option<Vec<VertexId>>,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) -> Triangles {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "wrap_and_triangulate_2d_normalized_vertices").entered();

    // Sort points into bins. Cover the region to be triangulated by a rectangular grid so that each bin contains roughly N^(1/2) points.
    // Label the bins so that consecutive bins are adjacent to one another, and then allocate each point to its appropriate bin.
    // Sort the list of points in ascending sequence of their bin numbers so that consecutive points are grouped together in the x-y plane.
    let partitioned_vertices = VertexBinSort::sort(&vertices, bin_vertex_density_power);

    let mut vertices_iterator = partitioned_vertices.iter().enumerate();
    // There is always at least 3 vertices
    let (_, first_vertex_id) = vertices_iterator.next().unwrap();
    let mut triangles = create_initial_triangles(vertices.len(), *first_vertex_id);
    // Id of the triangle we are looking at
    let mut triangle_id = triangles.last_id();

    #[cfg(feature = "debug_context")]
    debug_context.push_snapshot(Phase::FirstVertexInsertion, &triangles, &[0], &[]);

    // This buffer is used by all calls to `restore_delaunay_triangulation`.
    // We create if here to share the allocation between all those calls as an optimization.
    let mut quads_to_check = Vec::new();

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
        let Ok(vertex_place) = find_vertex_placement(
            vertex,
            triangle_id,
            &triangles,
            &vertices,
            #[cfg(feature = "debug_context")]
            debug_context,
        ) else {
            // TODO Internal error
            error!(
                "Internal error, found no triangle containing vertex {:?}, step {}",
                vertex_id, _index
            );
            break;
        };

        let new_triangles = match vertex_place {
            VertexPlacement::InsideTriangle(enclosing_triangle_id) => {
                // Form three new triangles by connecting P to each of the enclosing triangle's vertices.
                split_triangle_into_three_triangles(
                    &mut triangles,
                    enclosing_triangle_id,
                    vertex_id,
                    #[cfg(feature = "debug_context")]
                    debug_context,
                )
            }
            VertexPlacement::OnTriangleEdge(enclosing_triangle_id, edge_index) => {
                split_quad_into_four_triangles(
                    &mut triangles,
                    enclosing_triangle_id,
                    edge_index,
                    vertex_id,
                    #[cfg(feature = "debug_context")]
                    debug_context,
                )
            }
            VertexPlacement::OnVertex(existing_vertex_id) => {
                if let Some(vertex_merge_mapping) = vertex_merge_mapping {
                    vertex_merge_mapping[vertex_id as usize] = existing_vertex_id;
                }
                continue;
            }
        };

        restore_delaunay_triangulation(
            &mut triangles,
            vertices,
            vertex_id,
            new_triangles,
            &mut quads_to_check,
            #[cfg(feature = "debug_context")]
            debug_context,
        );

        // We'll start the search for the next enclosing triangle from the last created triangle.
        // This is a pretty good heuristic since the vertices were spatially partitionned
        triangle_id = triangles.last_id();

        #[cfg(feature = "progress_log")]
        {
            if _index % ((partitioned_vertices.len() / 50) + 1) == 0 {
                let progress = 100. * _index as f32 / partitioned_vertices.len() as f32;
                info!(
                    "Triangulation progress, step n°{}, {}%: {}/{}",
                    debug_context.current_step,
                    progress,
                    _index,
                    partitioned_vertices.len()
                );
            }
        }
    }

    triangles
}

pub(crate) fn remove_wrapping(
    triangles: &Triangles,
    config: &TriangulationConfiguration,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) -> Vec<[VertexId; 3]> {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "remove_wrapping").entered();

    #[cfg(feature = "debug_context")]
    let mut filtered_debug_triangles = Triangles::new();

    let indices = if triangles.count() > config.filter_parallel_tri_count_threshold {
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
                true => Some(t.verts),
                false => None,
            })
            // TODO Is it possible to pre-set the capacity of indices ?
            // .collect_into_vec(&mut indices);
            .collect()
    } else {
        let mut indices = Vec::with_capacity(triangles.count());
        for t in triangles.buffer().iter() {
            if t.is_finite() {
                indices.push(t.verts);
                #[cfg(feature = "debug_context")]
                filtered_debug_triangles.push(t.clone());
            }
        }
        indices
    };

    #[cfg(feature = "debug_context")]
    debug_context.push_snapshot(Phase::FilterTriangles, &filtered_debug_triangles, &[], &[]);

    indices
}

pub(crate) struct VertexBinSort {
    bins_per_row: usize,
    bins_count: usize,
}

impl VertexBinSort {
    // Each bin will contain roughly vertices.len()^(vertex_density_power) vertices
    pub fn sort(vertices: &Vec<Vertex>, vertex_density_power: f64) -> Vec<VertexId> {
        #[cfg(feature = "profile_traces")]
        let _span = span!(Level::TRACE, "sort").entered();

        let bins_per_row = (vertices.len() as f64)
            .powf(vertex_density_power / 2.)
            .round() as usize;

        let bin_sort = Self {
            bins_per_row,
            bins_count: bins_per_row * bins_per_row,
        };

        // Indexes of the bin corresponding to each vertex
        let mut vertices_bin_indexes = Vec::with_capacity(vertices.len());
        // Will be used to represent the index of a sorted vertex
        let mut bins_counters = vec![0; bin_sort.bins_count];

        for vertex in vertices {
            let bin_index = bin_sort.bin_index_from_vertex(*vertex);
            vertices_bin_indexes.push(bin_index);
            // Start by counting vertices in each bin
            bins_counters[bin_index] += 1;
        }

        // Add the amount of vertices in previous bins to each bin
        for bin_index in 1..bin_sort.bins_count {
            bins_counters[bin_index] += bins_counters[bin_index - 1];
        }

        let mut sorted = vec![0; vertices.len()];
        for vertex_id in 0..vertices.len() {
            let bin_index = vertices_bin_indexes[vertex_id];
            bins_counters[bin_index] -= 1;
            sorted[bins_counters[bin_index]] = vertex_id as VertexId;
        }
        sorted
    }

    fn bin_index_from_vertex(&self, vertex: Vertex) -> usize {
        // Compute a bin index from a vertex position which is in [0.,1]
        let bin_x = (0.99 * self.bins_per_row as Float * vertex.x) as usize;
        let bin_y: usize = (0.99 * self.bins_per_row as Float * vertex.y) as usize;
        self.bin_index_from_bin_position(bin_x, bin_y)
    }

    // Label the bins so that bins with consecutive indexes are spatially adjacent to one another
    fn bin_index_from_bin_position(&self, x: usize, y: usize) -> usize {
        if y % 2 == 0 {
            (y * self.bins_per_row) + x
        } else {
            (y + 1) * self.bins_per_row - x - 1
        }
    }
}

/// Splits a quad (triangle pair) into 4 triangles
///```text
///      from
///      /|\
///     / | \
///    /  |  \
///   /   |   \
///  /    |    \
/// /  t1 | t2  \
/// \     |     /
///  \    |    /
///   \   |   /
///    \  |  /
///     \ | /
///      \|/
///       to
///
/// Becomes
///
///      /|\
///     /3|2\
///    /  |  \
///   / t1|t2 \
///  /    |    \
/// /2___1|1___3\
/// \3   1|1   2/
///  \ t3 | t4 /
///   \   |   /
///    \  |  /
///     \2|3/
///      \|/
/// ```
pub(crate) fn split_quad_into_four_triangles(
    triangles: &mut Triangles,
    triangle_id: TriangleId,
    edge_index: TriangleEdgeIndex,
    vertex_id: VertexId,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) -> Vec<TriangleId> {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "split_quad_into_four_triangles").entered();

    // Re-use the existing triangle id for the first triangle
    let t1 = triangle_id;
    // TODO Comment is not true anymore in the general sense, but still true in this function
    // The triangle is always guaranteed to have a neighbor on the edge since a vertex
    // cannot land on a container triangle's edge.
    let t2 = triangles.get(triangle_id).neighbor(edge_index);
    // Create two new triangles for the other two
    let t3 = triangles.next_id();
    let t4 = triangles.next_id() + 1;

    let edge = triangles.get(t1).edge(edge_index);

    // t3
    let t3_neighbor_23 = triangles
        .get(t1)
        .neighbor(next_clockwise_edge_index(edge_index));
    let t3_v3 = triangles.get(t1).v(opposite_vertex_index(edge_index));
    triangles.create(
        [vertex_id, edge.to, t3_v3],
        [t4.into(), t3_neighbor_23, t1.into()],
    );

    let t2_opposite_v_id = triangles.get(t2.id).get_opposite_vertex_id(&edge);
    let t2_opposite_v_index = triangles.get(t2.id).vertex_index(t2_opposite_v_id);

    // t4
    let t4_neighbor_23 = triangles
        .get(t2.id)
        .neighbor(vertex_next_cw_edge_index(t2_opposite_v_index));
    triangles.create(
        [vertex_id, t2_opposite_v_id, edge.to],
        [t2.into(), t4_neighbor_23, t3.into()],
    );

    // Update triangle indexes
    update_triangle_neighbor(t3_neighbor_23, t1.into(), t3.into(), triangles);
    update_triangle_neighbor(t4_neighbor_23, t2.into(), t4.into(), triangles);

    // Update t1 verts
    let verts = [vertex_id, t3_v3, edge.from];
    triangles.get_mut(t1).verts = verts;
    // Update t1 neighbors
    let neighbors = [
        t3.into(),
        triangles
            .get(t1)
            .neighbor(next_counter_clockwise_edge_index(edge_index)),
        t2.into(),
    ];
    triangles.get_mut(t1).neighbors = neighbors;

    // Update t2 verts
    let verts = [vertex_id, edge.from, t2_opposite_v_id];
    triangles.get_mut(t2.id).verts = verts;
    // Update t2 neighbors
    let neighbors = [
        t1.into(),
        triangles
            .get(t2.id)
            .neighbor(vertex_next_ccw_edge_index(t2_opposite_v_index)),
        t4.into(),
    ];
    triangles.get_mut(t2.id).neighbors = neighbors;

    #[cfg(feature = "debug_context")]
    debug_context.push_snapshot_event(
        Phase::SplitQuad,
        EventInfo::Split(vertex_id),
        &triangles,
        &[t1, t2.id, t3, t4],
        // Neighbor data is out of date here and is not that interesting
        &[],
    );

    vec![t1, t2.id, t3, t4]
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
///            / t1  |  t3 \
///           /     1|1     \
///          /      /1\      \
///         /     /     \     \
///        /    /         \    \
///       /   /             \   \
///      /2 /        t2       \ 3\
///     / / 3                 2 \ \
///   v3 ------------------------- v2
/// ```
pub(crate) fn split_triangle_into_three_triangles(
    triangles: &mut Triangles,
    triangle_id: TriangleId,
    vertex_id: VertexId,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) -> Vec<TriangleId> {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "split_triangle_in_three_at_vertex").entered();

    // Re-use the existing triangle id for the first triangle
    let t1 = triangle_id;
    // Create two new triangles for the other two
    let t2 = triangles.next_id();
    let t3 = triangles.next_id() + 1;

    // t2
    triangles.create(
        [vertex_id, triangles.get(t1).v2(), triangles.get(t1).v3()],
        [t3.into(), triangles.get(t1).neighbor23(), t1.into()],
    );
    // t3
    triangles.create(
        [vertex_id, triangles.get(t1).v1(), triangles.get(t1).v2()],
        [t1.into(), triangles.get(t1).neighbor12(), t2.into()],
    );

    // Update triangle indexes
    update_triangle_neighbor(
        triangles.get(t1).neighbor12(),
        t1.into(),
        t3.into(),
        triangles,
    );
    update_triangle_neighbor(
        triangles.get(t1).neighbor23(),
        t1.into(),
        t2.into(),
        triangles,
    );

    // Update t1 verts
    let verts = [vertex_id, triangles.get(t1).v3(), triangles.get(t1).v1()];
    triangles.get_mut(t1).verts = verts;
    // Update t1 neighbors
    let neighbors = [t2.into(), triangles.get(t1).neighbor31(), t3.into()];
    triangles.get_mut(t1).neighbors = neighbors;

    #[cfg(feature = "debug_context")]
    debug_context.push_snapshot_event(
        Phase::SplitTriangle,
        EventInfo::Split(vertex_id),
        &triangles,
        &[t1, t2, t3],
        // Neighbor data is out of date here and is not that interesting
        &[],
    );

    vec![t1, t2, t3]
}

pub(crate) fn update_triangle_neighbor(
    triangle: Neighbor,
    old_neighbour_id: Neighbor,
    new_neighbour_id: Neighbor,
    triangles: &mut Triangles,
) {
    if triangle.exists() {
        for neighbor in triangles.get_mut(triangle.id).neighbors.iter_mut() {
            if *neighbor == old_neighbour_id {
                *neighbor = new_neighbour_id;
                break;
            }
        }
    }
}

/// `quads_to_check` used the shared pre-allocated buffer for efficiency.
/// - It does not need to be cleared since it is fully emptied by each call to restore_delaunay_triangulation
fn restore_delaunay_triangulation(
    triangles: &mut Triangles,
    vertices: &Vec<Vertex>,
    from_vertex_id: VertexId,
    new_triangles: Vec<TriangleId>,
    quads_to_check: &mut Vec<(TriangleId, TriangleId)>,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "restore_delaunay_triangulation").entered();

    for &from_triangle_id in new_triangles.iter() {
        // EDGE_23 is the opposite edge of `from_vertex_id` in the 3 new triangles
        let neighbor = triangles.get(from_triangle_id).neighbor23();
        if neighbor.exists() {
            quads_to_check.push((from_triangle_id, neighbor.id));
        }
    }

    while let Some((from_triangle_id, opposite_triangle_id)) = quads_to_check.pop() {
        match check_and_swap_quad_diagonal(
            triangles,
            vertices,
            from_vertex_id,
            from_triangle_id,
            opposite_triangle_id,
            #[cfg(feature = "debug_context")]
            debug_context,
        ) {
            QuadSwapResult::Swapped(quad_1, quad_2) => {
                // Place any new triangles pairs which are now opposite to `from_vertex_id` on the stack, to be checked
                // Unrolled loop for performances
                if quad_1.1.exists() {
                    quads_to_check.push((quad_1.0, quad_1.1.id));
                }
                if quad_2.1.exists() {
                    quads_to_check.push((quad_2.0, quad_2.1.id));
                }
            }
            QuadSwapResult::NotSwapped => (),
        }
    }
}

#[cold]
fn is_vertex_in_half_plane_1(
    vertices: &Vec<Vertex>,
    quad: &Quad,
    infinite_vert: TriangleVertexIndex,
) -> bool {
    // Test if q4 is inside the circle with 1 infinite point (half-plane defined by the 2 finite points)
    // TODO Clean: utils functions in Quad/Triangle
    let edge_index = opposite_edge_index(infinite_vert);
    let finite_vert_indexes = EDGE_TO_VERTS[edge_index as usize];
    // q1q2q3 is in a CCW order, so we reverse the edge
    let edge_vertices = (
        vertices[quad.v(finite_vert_indexes[1]) as usize],
        vertices[quad.v(finite_vert_indexes[0]) as usize],
    );

    // Strict seems to be necessary here, in order to avoid creating flat triangles
    is_point_strictly_on_right_side_of_edge(edge_vertices, vertices[quad.v4() as usize])
}

/// Test if q4 is inside the circle with 2 infinite points (half-plane defined by the finite point and the slope between the 2 infinite points)
#[cold]
fn is_vertex_in_half_plane_2(
    vertices: &Vec<Vertex>,
    quad: &Quad,
    infinite_v1: TriangleVertexIndex,
    infinite_v2: TriangleVertexIndex,
) -> bool {
    const INFINITE_VERTS_SLOPES: [Float; 3] = [1., -1., 1.];

    let infinite_vert_1_local_index = infinite_vertex_local_quad_index(quad.v(infinite_v1));
    let infinite_vert_2_local_index = infinite_vertex_local_quad_index(quad.v(infinite_v2));

    // Index of the finite vertex in q1q2q3
    let finite_vert_index = opposite_vertex_index_from_edge(infinite_v1, infinite_v2);

    let line_point_1 = vertices[quad.v(finite_vert_index) as usize];
    // Sum/2 defines a convenient surjection to (0..2)
    let a = INFINITE_VERTS_SLOPES
        [((infinite_vert_1_local_index + infinite_vert_2_local_index) / 2) as usize];
    let b = line_point_1.y - a * line_point_1.x;
    // Another point on the line, its position depends of which infinite vertices are considered. We need delta > 0 for Q1-Q2 and Q2-Q3, negative otherwise.
    let delta_x = if infinite_vert_1_local_index == QUAD_4 || infinite_vert_2_local_index == QUAD_4
    {
        DELTA_VALUE
    } else {
        -DELTA_VALUE
    };
    let line_point_2_x = line_point_1.x + delta_x;
    let line_point_2 = Vertex::new(line_point_2_x, a * line_point_2_x + b);
    // Strict seems to be necessary here, in order to avoid creating flat triangles
    is_point_strictly_on_right_side_of_edge(
        (line_point_1, line_point_2),
        vertices[quad.v4() as usize],
    )
}

#[inline(always)]
pub(crate) fn should_swap_diagonals(quad: &Quad, vertices: &Vec<Vertex>) -> bool {
    // TODO Performance: try stack/pre-alloc
    let mut infinite_verts = Vec::new();
    if is_infinite(quad.v1()) {
        infinite_verts.push(QUAD_1);
    }
    if is_infinite(quad.v2()) {
        infinite_verts.push(QUAD_2);
    }
    if is_infinite(quad.v3()) {
        infinite_verts.push(QUAD_3);
    }
    // TODO Fix: Remove this debug code
    // if is_infinite(quad.v4(), min_container_vertex_id) {
    //     error!("Temporary debug log, infinite q4, should not be possible");
    // }

    // TODO Performance: try to play with the #cold attribute for the infinite cases.
    if infinite_verts.is_empty() {
        // General case: no infinite vertices
        // Test if `from_vertex_id` is inside the circumcircle of `opposite_triangle`
        let quad_vertices = quad.to_vertices(vertices);
        is_vertex_in_triangle_circumcircle(&quad_vertices.verts[0..=2], quad_vertices.q4())
    } else if infinite_verts.len() == 1 {
        is_vertex_in_half_plane_1(vertices, quad, infinite_verts[0])
    } else {
        is_vertex_in_half_plane_2(vertices, quad, infinite_verts[0], infinite_verts[1])
    }
    // 3 infinite vertices is not possible by construction, the container triangle is split into 3 triangles as soon as the first point is inserted.
}

#[derive(PartialEq, Eq, Debug)]
pub enum QuadSwapResult {
    /// Contains the new triangle pairs to check
    Swapped((TriangleId, Neighbor), (TriangleId, Neighbor)),
    NotSwapped,
}

/// ```text
///                q3
///         t3   /    \   t4
///            /   To   \
///          /            \
///         q1 ---------- q2
///          \ 2        3 /
///            \   Tf   /
///              \ 1  /
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
///         t3  / 3|2 \   t4
///           /    |    \
///         /      |      \
///        q1 2  Tf|To   3 q2
///         \      |      /
///           \    |    /
///             \ 1|1 /
///               q4
/// ```
pub(crate) fn check_and_swap_quad_diagonal(
    triangles: &mut Triangles,
    vertices: &Vec<Vertex>,
    from_vertex_id: VertexId,
    from_triangle_id: TriangleId,
    opposite_triangle_id: TriangleId,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) -> QuadSwapResult {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "check_and_swap_quad_diagonal").entered();

    let opposite_triangle = triangles.get(opposite_triangle_id);

    let (quad, triangle_3, triangle_4) =
    // No need to check if neighbor exists, handled by the == check since `from_triangle_id` exists
        if opposite_triangle.neighbor12().id == from_triangle_id {
            (
                Quad::new([
                    opposite_triangle.v2(),
                    opposite_triangle.v1(),
                    opposite_triangle.v3(),
                    from_vertex_id,
                ]),
                opposite_triangle.neighbor23(),
                opposite_triangle.neighbor31(),
            )
        } else if opposite_triangle.neighbor23().id == from_triangle_id {
            (
                Quad::new([
                    opposite_triangle.v3(),
                    opposite_triangle.v2(),
                    opposite_triangle.v1(),
                    from_vertex_id,
                ]),
                opposite_triangle.neighbor31(),
                opposite_triangle.neighbor12(),
            )
        } else {
            (
                Quad::new([
                    opposite_triangle.v1(),
                    opposite_triangle.v3(),
                    opposite_triangle.v2(),
                    from_vertex_id,
                ]),
                opposite_triangle.neighbor12(),
                opposite_triangle.neighbor23(),
            )
        };

    if should_swap_diagonals(&quad, vertices) {
        let opposite_neighbor = opposite_triangle_id.into();
        let from_neighbor = from_triangle_id.into();

        update_triangle_neighbor(triangle_3, opposite_neighbor, from_neighbor, triangles);
        update_triangle_neighbor(
            triangles.get(from_triangle_id).neighbor31(),
            from_neighbor,
            opposite_neighbor,
            triangles,
        );

        triangles.get_mut(from_triangle_id).verts = [quad.v4(), quad.v1(), quad.v3()];
        triangles.get_mut(opposite_triangle_id).verts = [quad.v4(), quad.v3(), quad.v2()];

        triangles.get_mut(opposite_triangle_id).neighbors = [
            from_neighbor,
            triangle_4,
            triangles.get(from_triangle_id).neighbor31(),
        ];
        *triangles.get_mut(from_triangle_id).neighbor23_mut() = triangle_3;
        *triangles.get_mut(from_triangle_id).neighbor31_mut() = opposite_neighbor;

        #[cfg(feature = "debug_context")]
        debug_context.push_snapshot(
            Phase::DelaunayRestoreSwapQuadDiagonals,
            &triangles,
            &[from_triangle_id, opposite_triangle_id],
            &[triangle_3, triangles.get(from_triangle_id).neighbor31()],
        );

        QuadSwapResult::Swapped(
            (from_triangle_id, triangle_3),
            (opposite_triangle_id, triangle_4),
        )
    } else {
        QuadSwapResult::NotSwapped
    }
}

///////////////////////////////////////////////////////////
///                                                     ///
///                        Tests                        ///
///                                                     ///
///////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    #[cfg(feature = "debug_context")]
    use crate::debug::{DebugConfiguration, DebugContext};
    use crate::{
        triangulation::{
            check_and_swap_quad_diagonal, normalize_vertices_coordinates,
            transform_to_2d_planar_coordinate_system, QuadSwapResult,
        },
        types::{Float, Neighbor, TriangleData, Triangles, Vector3A, Vertex},
    };

    #[test]
    fn triangulation_normalize_set_of_vertices() {
        let mut vertices = Vec::<[Float; 2]>::new();

        vertices.push([3.0, 2.0]);
        vertices.push([-1.0, 2.0]);
        vertices.push([-1.0, -2.0]);
        vertices.push([3.0, -2.0]);

        let mut planar_vertices = Vec::with_capacity(vertices.len());
        for v in &vertices {
            planar_vertices.push(Vertex::from_array(*v));
        }

        normalize_vertices_coordinates(&mut planar_vertices);

        assert_eq!(
            Vec::from([
                Vertex::from([1., 1.]),
                Vertex::from([0., 1.]),
                Vertex::from([0., 0.]),
                Vertex::from([1., 0.])
            ]),
            planar_vertices
        );
    }

    #[test]
    fn triangulation_set_to_2d_plane_vertices() {
        let mut vertices = Vec::<[Float; 3]>::new();

        vertices.push([-3., 2., 0.]);
        vertices.push([1., 2., 0.]);
        vertices.push([1., -2., 0.]);
        vertices.push([-3., -2., 0.]);

        let plane_normal = Vector3A::Z;

        let mut vertices_data = Vec::with_capacity(vertices.len());
        for v in &vertices {
            vertices_data.push(Vector3A::from_array(*v));
        }

        let planar_vertices =
            transform_to_2d_planar_coordinate_system(&mut vertices_data, plane_normal);

        assert_eq!(
            Vec::from([
                Vertex::from([3.0, 2.0]),
                Vertex::from([-1.0, 2.0]),
                Vertex::from([-1.0, -2.0]),
                Vertex::from([3.0, -2.0])
            ]),
            planar_vertices
        );
    }

    #[test]
    fn no_swap() {
        let mut vertices = Vec::<Vertex>::new();
        vertices.push(Vertex::new(0.5, 3.));
        vertices.push(Vertex::new(-2., -2.));
        vertices.push(Vertex::new(1., -4.));
        vertices.push(Vertex::new(3., -2.));

        let triangle_1 = TriangleData {
            verts: [3, 1, 0],
            neighbors: [Neighbor::NONE, Neighbor::NONE, Neighbor::new(1)],
        };

        let triangle_2 = TriangleData {
            verts: [1, 2, 3],
            neighbors: [Neighbor::NONE, Neighbor::NONE, Neighbor::new(0)],
        };

        let mut triangles = Triangles::with_capacity(2);
        triangles.push(triangle_1);
        triangles.push(triangle_2);

        #[cfg(feature = "debug_context")]
        let mut debug_context = DebugContext::new(DebugConfiguration::default(), 0., 0., 0.);
        let quad_swap = check_and_swap_quad_diagonal(
            &mut triangles,
            &vertices,
            1,
            0,
            1,
            #[cfg(feature = "debug_context")]
            &mut debug_context,
        );

        assert_eq!(QuadSwapResult::NotSwapped, quad_swap);
        assert_eq!(2, triangles.count());
    }

    #[test]
    fn swap() {
        let mut vertices = Vec::<Vertex>::new();
        vertices.push(Vertex::new(0.5, 3.));
        vertices.push(Vertex::new(-2., -2.));
        vertices.push(Vertex::new(1., -4.));
        vertices.push(Vertex::new(3., -2.));

        let triangle_1 = TriangleData {
            verts: [0, 1, 2],
            neighbors: [Neighbor::NONE, Neighbor::NONE, Neighbor::new(1)],
        };

        let triangle_2 = TriangleData {
            verts: [2, 3, 0],
            neighbors: [Neighbor::NONE, Neighbor::NONE, Neighbor::new(0)],
        };

        let mut triangles = Triangles::with_capacity(2);
        triangles.push(triangle_1);
        triangles.push(triangle_2);

        #[cfg(feature = "debug_context")]
        let mut debug_context = DebugContext::new(DebugConfiguration::default(), 0., 0., 0.);
        let quad_swap = check_and_swap_quad_diagonal(
            &mut triangles,
            &vertices,
            1,
            0,
            1,
            #[cfg(feature = "debug_context")]
            &mut debug_context,
        );

        assert_ne!(QuadSwapResult::NotSwapped, quad_swap);
        assert_eq!(2, triangles.count());
    }
}
