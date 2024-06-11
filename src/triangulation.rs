use log::error;
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};

use crate::types::{
    opposite_edge_index, Float, Neighbor, Quad, QuadVertices, TriangleData, TriangleId,
    TriangleVertexIndex, Triangles, Vector3A, Vertex, VertexId, EDGE_TO_VERTS, VERT_1, VERT_2,
    VERT_3,
};
use crate::utils::{
    is_point_on_right_side_of_edge, is_vertex_in_triangle_circumcircle, line_slope,
};

#[cfg(feature = "progress_log")]
use log::info;

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

pub const CONTAINER_TRIANGLE_COORDINATE: Float = 5.;

pub const CONTAINER_TRIANGLE_VERTICES: [Vertex; 3] = [
    Vertex::new(
        -CONTAINER_TRIANGLE_COORDINATE,
        -CONTAINER_TRIANGLE_COORDINATE,
    ),
    Vertex::new(0., CONTAINER_TRIANGLE_COORDINATE),
    Vertex::new(
        CONTAINER_TRIANGLE_COORDINATE,
        -CONTAINER_TRIANGLE_COORDINATE,
    ),
];

/// plane_normal must be normalized
/// vertices must all belong to a 3d plane
pub fn triangulation_from_3d_planar_vertices(
    vertices: &Vec<[Float; 3]>,
    plane_normal: Vector3A,
    config: TriangulationConfiguration,
) -> Triangulation {
    // TODO Clean: See what we need for input data format of `triangulate`
    let mut vertices_data = Vec::with_capacity(vertices.len());
    for v in vertices {
        vertices_data.push(Vector3A::from_array(*v));
    }

    let mut planar_vertices =
        transform_to_2d_planar_coordinate_system(&mut vertices_data, plane_normal);

    triangulation_from_2d_vertices(&mut planar_vertices, config)
}
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

    // Uniformly scale the coordinates of the points so that they all lie between 0 and 1.
    let (mut normalized_vertices, _scale_factor, _x_min, _y_min) =
        normalize_vertices_coordinates(vertices);

    #[cfg(feature = "debug_context")]
    let mut debug_context =
        DebugContext::new(config.debug_config.clone(), _scale_factor, _x_min, _y_min);

    let (triangles, min_container_vertex_id) = wrap_and_triangulate_2d_normalized_vertices(
        &mut normalized_vertices,
        config.bin_vertex_density_power,
        &mut None,
        #[cfg(feature = "debug_context")]
        &mut debug_context,
    );

    let vert_indices = remove_wrapping(
        &triangles,
        min_container_vertex_id,
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
pub(crate) fn transform_to_2d_planar_coordinate_system(
    vertices: &mut Vec<Vector3A>,
    plane_normal: Vector3A,
) -> Vec<Vertex> {
    // Create a base, using the first two vertices as the first base vector and plane_normal as the second
    let basis_1 = (vertices[0] - vertices[1]).normalize();
    // basis_3 is already normalized since basis_1 and plane_normal are normalized and orthogonal
    let basis_3 = basis_1.cross(plane_normal);

    // Project every vertices into the base B
    let mut vertices_2d = Vec::with_capacity(vertices.len());
    for vertex in vertices {
        vertices_2d.push(Vertex::new(vertex.dot(basis_1), vertex.dot(basis_3)));
    }
    vertices_2d
}

/// This scaling ensures that all of the coordinates are between 0 and 1 but does not modify the relative positions of the points in the x-y plane.
/// The use of normalized coordinates, although not essential, reduces the effects of  roundoff error and is also convenient from a computational point of view.
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
enum SearchResult {
    EnlosingTriangle(TriangleId),
    NotFound,
}

fn search_enclosing_triangle(
    vertex: Vertex,
    from: Neighbor,
    triangles: &Triangles,
    vertices: &Vec<Vertex>,
) -> SearchResult {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "search_enclosing_triangle").entered();

    let mut neighbor = from;

    let mut search_result = SearchResult::NotFound;
    // We use `triangles.len()` as an upper bound on the number of triangles
    for _ in 0..triangles.count() {
        if !neighbor.exists() {
            break;
        }
        let triangle = triangles.get(neighbor.id);
        let (v1, v2, v3) = triangle.to_vertices(vertices);

        // Check if the point is inside the triangle, if not check the neighbours
        if !is_point_on_right_side_of_edge((v1, v2), vertex) {
            neighbor = triangle.neighbor12();
        } else if !is_point_on_right_side_of_edge((v2, v3), vertex) {
            neighbor = triangle.neighbor23();
        } else if !is_point_on_right_side_of_edge((v3, v1), vertex) {
            neighbor = triangle.neighbor31();
        } else {
            search_result = SearchResult::EnlosingTriangle(neighbor.id);
            break;
        }
    }

    search_result
}

/// Select three dummy points to form a supertriangle that completely encompasses all of the points to be triangulated.
///  This supertriangle initially defines a Delaunay triangulation which is comprised of a single triangle.
///  Its vertices are defined in terms of normalized coordinates and are usually located at a considerable distance from the window which encloses the set of points.
pub(crate) fn add_container_triangle_vertices(
    vertices: &mut Vec<Vertex>,
) -> (TriangleData, VertexId) {
    let min_container_triangle_vertex_id = vertices.len() as VertexId;
    let container_triangle = TriangleData::new_container_triangle(min_container_triangle_vertex_id);
    vertices.extend(CONTAINER_TRIANGLE_VERTICES.clone());
    (container_triangle, min_container_triangle_vertex_id)
}

pub(crate) fn find_existing_close_vertex(
    triangle: &TriangleData,
    triangle_verts: &[Vertex; 3],
    vertex: Vertex,
) -> Option<VertexId> {
    for (vertex_index, triangle_vertex) in triangle_verts.iter().enumerate() {
        let dist = *triangle_vertex - vertex;
        if dist.x.abs() < Float::EPSILON && dist.y.abs() < Float::EPSILON {
            return Some(triangle.verts[vertex_index]);
        }
    }
    None
}

/// - `vertices` should be normalized with their cooridnates in [0,1]
pub(crate) fn wrap_and_triangulate_2d_normalized_vertices(
    vertices: &mut Vec<Vertex>,
    bin_vertex_density_power: f64,
    vertex_merge_mapping: &mut Option<Vec<VertexId>>,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) -> (Triangles, VertexId) {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "wrap_and_triangulate_2d_normalized_vertices").entered();

    // Sort points into bins. Cover the region to be triangulated by a rectangular grid so that each bin contains roughly N^(1/2) points.
    // Label the bins so that consecutive bins are adjacent to one another, and then allocate each point to its appropriate bin.
    // Sort the list of points in ascending sequence of their bin numbers so that consecutive points are grouped together in the x-y plane.
    let partitioned_vertices = VertexBinSort::sort(&vertices, bin_vertex_density_power);

    let (container_triangle, min_container_vertex_id) = add_container_triangle_vertices(vertices);

    let mut triangles = Triangles::with_capacity(vertices.len() * 2 + 1);
    triangles.buffer_mut().push(container_triangle);

    // Id of the triangle we are looking at
    let mut triangle_id = Neighbor::new(0);

    #[cfg(feature = "debug_context")]
    debug_context.push_snapshot(Phase::ContainerVerticesInsertion, &triangles, &[0], &[]);

    // This buffer is used by all calls to `restore_delaunay_triangulation`.
    // We create if here to share the allocation between all those calls as an optimization.
    let mut quads_to_check = Vec::new();

    // Loop over all the input vertices
    for (_index, &vertex_id) in partitioned_vertices.iter().enumerate() {
        #[cfg(feature = "debug_context")]
        {
            let force_end = debug_context.advance_step();
            if force_end {
                break;
            }
        }

        // Find an existing triangle which encloses P
        let vertex = vertices[vertex_id as usize];
        match search_enclosing_triangle(vertex, triangle_id, &triangles, &vertices) {
            SearchResult::EnlosingTriangle(enclosing_triangle_id) => {
                // Compare to the points in the triangle, if too close to one, continue
                let triangle = &triangles.get(enclosing_triangle_id);
                let triangle_verts = triangle.to_vertices_array(vertices);
                if let Some(existing_vertex_id) =
                    find_existing_close_vertex(&triangle, &triangle_verts, vertex)
                {
                    if let Some(vertex_merge_mapping) = vertex_merge_mapping {
                        vertex_merge_mapping[vertex_id as usize] = existing_vertex_id;
                    }
                    continue;
                }

                // Form three new triangles by connecting P to each of the enclosing triangle's vertices.
                let new_triangles = split_triangle_in_three_at_vertex(
                    &mut triangles,
                    enclosing_triangle_id,
                    vertex_id,
                    #[cfg(feature = "debug_context")]
                    debug_context,
                );

                restore_delaunay_triangulation(
                    &mut triangles,
                    vertices,
                    min_container_vertex_id,
                    vertex_id,
                    new_triangles,
                    &mut quads_to_check,
                    #[cfg(feature = "debug_context")]
                    debug_context,
                );

                // We'll start the search for the next enclosing triangle from the last created triangle.
                // This is a pretty good heuristic since the vertices were spatially partitionned
                triangle_id = Neighbor::new(triangles.last_id());

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
            SearchResult::NotFound => {
                // TODO Error
                error!(
                    "Found no triangle enclosing vertex {:?}, step {}",
                    vertex_id, _index
                );
                break;
            }
        }
    }

    (triangles, min_container_vertex_id)
}

pub(crate) fn remove_wrapping(
    triangles: &Triangles,
    min_container_vertex_id: VertexId,
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
        for triangle in triangles.buffer().iter() {
            let mut filtered = false;
            for &vert in triangle.verts.iter() {
                if vert >= min_container_vertex_id {
                    filtered = true;
                    break;
                }
            }
            if !filtered {
                filtered_debug_triangles.push(triangle.clone());
            }
        }

        triangles
            .buffer()
            .par_iter()
            .with_min_len(config.filter_parallel_min_batch_len)
            .filter_map(|t| {
                let mut filtered = false;
                for &vert in t.verts.iter() {
                    if vert >= min_container_vertex_id {
                        filtered = true;
                        break;
                    }
                }
                match filtered {
                    true => None,
                    false => Some(t.verts),
                }
            })
            // TODO Is it possible to pre-set the capacity of indices ?
            // .collect_into_vec(&mut indices);
            .collect()
    } else {
        let mut indices = Vec::with_capacity(triangles.count());
        for t in triangles.buffer().iter() {
            let mut filtered = false;
            for &vert in t.verts.iter() {
                if vert >= min_container_vertex_id {
                    filtered = true;
                    break;
                }
            }
            if !filtered {
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

/// Splits `triangle_id` into 3 triangles (re-using the existing triangle id)
///
/// All the resulting triangles will share `vertex_id` as their frist vertex, and will be oriented in a CW order
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
pub(crate) fn split_triangle_in_three_at_vertex(
    triangles: &mut Triangles,
    triangle_id: TriangleId,
    vertex_id: VertexId,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) -> [TriangleId; 3] {
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

    let verts = [vertex_id, triangles.get(t1).v3(), triangles.get(t1).v1()];
    triangles.get_mut(t1).verts = verts;

    let neighbors = [t2.into(), triangles.get(t1).neighbor31(), t3.into()];
    triangles.get_mut(t1).neighbors = neighbors;

    #[cfg(feature = "debug_context")]
    debug_context.push_snapshot_event(
        Phase::SplitTriangle,
        EventInfo::SplitTriangle(vertex_id),
        &triangles,
        &[t1, t2, t3],
        // Neighbor data is out of date here and is not that interesting
        &[],
    );

    [t1, t2, t3]
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
    min_container_vertex_id: VertexId,
    from_vertex_id: VertexId,
    new_triangles: [TriangleId; 3],
    quads_to_check: &mut Vec<(TriangleId, TriangleId)>,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "restore_delaunay_triangulation").entered();

    for &from_triangle_id in &new_triangles {
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
            min_container_vertex_id,
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
    infinite_vert: TriangleVertexIndex,
    quad_vertices: &QuadVertices,
) -> bool {
    // Test if q4 is inside the circle with 1 inifite point (half-plane defined by the 2 finite points)
    // TODO Clean: utils functions in Quad/Triangle
    let infinite_vert = infinite_vert;
    let edge = opposite_edge_index(infinite_vert);
    let finite_vert_indexes = EDGE_TO_VERTS[edge as usize];
    // q1q2q3 is in a CCW order, so we reverse the edge
    let edge_vertices = (
        quad_vertices.verts[finite_vert_indexes[1] as usize],
        quad_vertices.verts[finite_vert_indexes[0] as usize],
    );
    is_point_on_right_side_of_edge(edge_vertices, quad_vertices.q4())
}

#[cold]
fn is_vertex_in_half_plane_2(
    infinite_vert_1: TriangleVertexIndex,
    infinite_vert_2: TriangleVertexIndex,
    quad_vertices: &QuadVertices,
) -> bool {
    // Test if q4 is inside the circle with 2 infinite points (half-plane defined by the finite point and the slope between the 2 infinite points)
    // Index of the finite vertex in q1q2q3
    let finite_vert_index = (3 - infinite_vert_1 + infinite_vert_2) as usize;
    let line_point = quad_vertices.verts[finite_vert_index];
    // TODO Improvement: Could use a slope LUT since we know container vertices as const
    let a = line_slope(
        quad_vertices.verts[infinite_vert_1 as usize],
        quad_vertices.verts[infinite_vert_2 as usize],
    );
    let b = line_point.y - a * line_point.x;
    // q1q2q3 is CCW
    // q1 q2 ou q3 q1 => y > a.x + b
    // q2 q3 => y < a.x + b
    if infinite_vert_1 == VERT_2 && infinite_vert_2 == VERT_3 {
        quad_vertices.q4().y < a * quad_vertices.q4().x + b
    } else {
        quad_vertices.q4().y > a * quad_vertices.q4().x + b
    }
}

#[inline(always)]
pub(crate) fn should_swap_diagonals(
    quad: &Quad,
    vertices: &Vec<Vertex>,
    min_container_vertex_id: VertexId,
) -> bool {
    let quad_vertices = quad.to_vertices(vertices);
    // TODO Performance: try stack/pre-alloc
    let mut infinite_verts = Vec::new();
    if quad.v1() >= min_container_vertex_id {
        infinite_verts.push(VERT_1);
    }
    if quad.v2() >= min_container_vertex_id {
        infinite_verts.push(VERT_2);
    }
    if quad.v3() >= min_container_vertex_id {
        infinite_verts.push(VERT_3);
    }

    // TODO Performance: try to play with the #cold attribute for the infinite cases.
    if infinite_verts.is_empty() {
        // General case: no infinite vertices
        // Test if `from_vertex_id` is inside the circumcircle of `opposite_triangle`
        is_vertex_in_triangle_circumcircle(&quad_vertices.verts[0..=2], quad_vertices.q4())
    } else if infinite_verts.len() == 1 {
        is_vertex_in_half_plane_1(infinite_verts[0], &quad_vertices)
    } else {
        is_vertex_in_half_plane_2(infinite_verts[0], infinite_verts[1], &quad_vertices)
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
    min_container_vertex_id: VertexId,
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

    if should_swap_diagonals(&quad, vertices, min_container_vertex_id) {
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
            split_triangle_in_three_at_vertex, transform_to_2d_planar_coordinate_system,
            QuadSwapResult,
        },
        types::{Float, Neighbor, TriangleData, TriangleId, Triangles, Vector3A, Vertex},
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
    fn split_in_three_triangle() {
        let mut vertices = Vec::<Vertex>::new();
        vertices.push(Vertex::new(0., 0.)); // vertex to be added

        // container triangle to be splited by the vertex
        let container_triangle = TriangleData::new_container_triangle(vertices.len() as TriangleId);

        // vertices of the container triangle
        vertices.extend([
            Vertex::new(1., 1.),
            Vertex::new(1., -2.),
            Vertex::new(-3., 2.),
        ]);

        let mut triangles = Triangles::new();
        triangles.push(container_triangle);

        #[cfg(feature = "debug_context")]
        let mut debug_context = DebugContext::new(DebugConfiguration::default(), 0., 0., 0.);
        let _new_triangles = split_triangle_in_three_at_vertex(
            &mut triangles,
            0,
            0,
            #[cfg(feature = "debug_context")]
            &mut debug_context,
        );

        assert_eq!(3, triangles.count());
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
            6,
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
            6,
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
