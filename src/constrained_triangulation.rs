use std::collections::VecDeque;

use hashbrown::HashSet;
use tracing::error;

use crate::triangulation::{
    edge_from_semi_infinite_edge, normalize_vertices_coordinates, should_swap_diagonals,
    Triangulation, DEFAULT_BIN_VERTEX_DENSITY_POWER,
};
use crate::types::{
    is_infinite, next_counter_clockwise_edge_index, QuadVertices, TriangleEdgeIndex, Triangles,
    Vector3, Vector3A, Vertex, ADJACENT_QUAD_VERTICES_INDEXES, QUAD_1, QUAD_2, QUAD_3, QUAD_4,
};
use crate::utils::{egdes_intersect, EdgesIntersectionResult};

#[cfg(feature = "debug_context")]
use crate::debug::{DebugConfiguration, DebugContext, Phase};

#[cfg(feature = "profile_traces")]
use tracing::{span, Level};

#[cfg(feature = "progress_log")]
use tracing::info;

use crate::{
    triangulation::{self, wrap_and_triangulate_2d_normalized_vertices},
    types::{
        next_clockwise_edge_index, opposite_edge_index, outermost_clockwise_edge_index_around,
        outermost_counter_clockwise_edge_index_around, Edge, EdgeVertices, Neighbor, Quad,
        TriangleData, TriangleId, TriangleVertexIndex, TriangleVertices, VertexId,
    },
};

#[derive(Clone, Debug)]
pub struct ConstrainedTriangulationConfiguration {
    /// Binsort will cover the region to be triangulated by a rectangular grid so that each bin contains roughly N^(density_power) points.
    pub bin_vertex_density_power: f64,

    #[cfg(feature = "debug_context")]
    pub debug_config: DebugConfiguration,
}

impl Default for ConstrainedTriangulationConfiguration {
    fn default() -> Self {
        Self {
            bin_vertex_density_power: DEFAULT_BIN_VERTEX_DENSITY_POWER,
            #[cfg(feature = "debug_context")]
            debug_config: DebugConfiguration::default(),
        }
    }
}

/// - `plane_normal` **MUST** be normalized
/// - `vertices` **MUST** all belong to a 3d plane
/// - `constrained_edges` **MUST**:
///     - be oriented: vi -> vj
///     - not contain edges with v0==v1
pub fn constrained_triangulation_from_3d_planar_vertices(
    vertices: &Vec<Vector3>,
    plane_normal: Vector3A,
    constrained_edges: &Vec<Edge>,
    config: ConstrainedTriangulationConfiguration,
) -> Triangulation {
    let mut planar_vertices =
        triangulation::transform_to_2d_planar_coordinate_system(vertices, plane_normal);
    constrained_triangulation_from_2d_vertices(&mut planar_vertices, constrained_edges, config)
}

/// - Constrained edges **MUST** not contain edges with from == to
///
/// As visible in the figure below,
/// - Contours (set of constrained edges) forming a CW cycle indicate a domain to triangulate.
/// - Contours forming a CCW cycle indicate a domain that should be discarded from the triangulation.
///
/// ```text
/// ------------->-----------------------------------
/// |                   keep (CW)                   |
/// | -----------<--------------------------------- |
/// | |                remove (CCW)               | |
/// | |   ------->--------     ------->--------   | |
/// | |   |              |     |              |   | |
/// ^ v   ^   keep (CW)  v     ^   keep (CW)  v   ^ v
/// | |   |              |     |              |   | |
/// | |   -------<--------     -------<--------   | |
/// | ----------->--------------------------------- |
/// -------------<-----------------------------------
/// ```
pub fn constrained_triangulation_from_2d_vertices(
    vertices: &Vec<Vertex>,
    constrained_edges: &Vec<Edge>,
    config: ConstrainedTriangulationConfiguration,
) -> Triangulation {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "constrained_triangulation_from_2d_vertices").entered();

    // Uniformly scale the coordinates of the points so that they all lie between 0 and 1.
    let (mut normalized_vertices, _scale_factor, _x_min, _y_min) =
        normalize_vertices_coordinates(vertices);

    #[cfg(feature = "debug_context")]
    let mut debug_context =
        DebugContext::new(config.debug_config.clone(), _scale_factor, _x_min, _y_min);

    let mut vertex_merge_mapping = Some((0..vertices.len() as VertexId).collect());
    let (mut triangles, min_container_vertex_id) = wrap_and_triangulate_2d_normalized_vertices(
        &mut normalized_vertices,
        config.bin_vertex_density_power,
        &mut vertex_merge_mapping,
        #[cfg(feature = "debug_context")]
        &mut debug_context,
    );

    #[cfg(feature = "debug_context")]
    debug_context.push_snapshot(Phase::BeforeConstraints, &triangles, &[], &[]);

    // TODO Debug: consider adding debug support to this function
    let constrained_edges_set = apply_constraints(
        &normalized_vertices,
        &mut triangles,
        min_container_vertex_id,
        constrained_edges,
        vertex_merge_mapping.unwrap(),
        #[cfg(feature = "debug_context")]
        &mut debug_context,
    );

    #[cfg(feature = "debug_context")]
    debug_context.push_snapshot(Phase::AfterConstraints, &triangles, &[], &[]);

    let vert_indices = remove_wrapping_and_unconstrained_domains(
        &triangles,
        min_container_vertex_id,
        constrained_edges_set,
        #[cfg(feature = "debug_context")]
        &mut debug_context,
    );

    Triangulation {
        triangles: vert_indices,
        #[cfg(feature = "debug_context")]
        debug_context,
    }
}

/// Filter triangles that should be removed due to the input domains constraints or due to being part of
/// the container triangle.
fn remove_wrapping_and_unconstrained_domains(
    triangles: &Triangles,
    min_container_vertex_id: VertexId,
    constrained_edges: HashSet<Edge>,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) -> Vec<[VertexId; 3]> {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "remove_wrapping_and_unconstrained_domains").entered();

    let mut visited_triangles = vec![false; triangles.count()];

    // TODO Clean: Size approx
    let mut indices = Vec::with_capacity(3 * triangles.count());
    let mut triangles_to_explore = Vec::new();

    #[cfg(feature = "debug_context")]
    let mut filtered_debug_triangles = Triangles::new();

    // Loop over all triangles
    for (index, triangle) in triangles.buffer().iter().enumerate() {
        if visited_triangles[index] {
            continue;
        }

        let triangle_edge_is_constrained = [
            constrained_edges.contains(&triangle.edge12()),
            constrained_edges.contains(&triangle.edge23()),
            constrained_edges.contains(&triangle.edge31()),
        ];

        // Stop at the first non-visited triangle with a constrained edge: we found a new unvisited domain
        if triangle_edge_is_constrained.contains(&true) {
            triangles_to_explore.clear();

            if triangle.has_no_container_vertex(min_container_vertex_id) {
                indices.push(triangle.verts);
                #[cfg(feature = "debug_context")]
                filtered_debug_triangles.push(triangle.clone());
            }
            for (edge_index, is_constrained) in triangle_edge_is_constrained.iter().enumerate() {
                if !is_constrained {
                    triangles_to_explore.push(triangle.neighbors[edge_index]);
                }
            }
            visited_triangles[index] = true;

            // Triangle-walk in the domain to register the domain's triangles
            while let Some(neighbor) = triangles_to_explore.pop() {
                if !neighbor.exists() {
                    continue;
                }
                let triangle_id = neighbor.id;

                if visited_triangles[triangle_id as usize] {
                    continue;
                } else {
                    let t = triangles.get(triangle_id);
                    if t.has_no_container_vertex(min_container_vertex_id) {
                        indices.push(t.verts);
                        #[cfg(feature = "debug_context")]
                        filtered_debug_triangles.push(t.clone());
                    }

                    visited_triangles[triangle_id as usize] = true;

                    for (edge_index, edge) in t.edges().iter().enumerate() {
                        if !constrained_edges.contains(edge) {
                            // TODO Clean: Would it be quicker to also check for visited before adding to the stack
                            triangles_to_explore.push(t.neighbors[edge_index]);
                        }
                    }
                }
            }
        }
    }

    #[cfg(feature = "debug_context")]
    debug_context.push_snapshot(Phase::FilterTriangles, &filtered_debug_triangles, &[], &[]);

    indices
}

fn apply_constraints(
    vertices: &Vec<Vertex>,
    triangles: &mut Triangles,
    min_container_vertex_id: VertexId,
    constrained_edges: &Vec<Edge>,
    vertex_merge_mapping: Vec<VertexId>,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) -> HashSet<Edge> {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "apply_constraints").entered();

    let mut constrained_edges_set = HashSet::with_capacity(constrained_edges.len());
    // Map each vertex to one of the triangles (the last) that contains it
    let mut vertex_to_triangle = vec![0; vertices.len()];
    for (index, t) in triangles.buffer().iter().enumerate() {
        vertex_to_triangle[t.v1() as usize] = index as TriangleId;
        vertex_to_triangle[t.v2() as usize] = index as TriangleId;
        vertex_to_triangle[t.v3() as usize] = index as TriangleId;
    }

    // Those buffers are used by all calls to `register_intersected_edges`, `remove_crossed_edges` and `restore_delaunay_triangulation_constrained`.
    // We create them here to share the allocation between all those calls as an optimization.
    let mut intersections = VecDeque::new();
    let mut new_diagonals_created = VecDeque::new();

    for (_edge_index, constrained_edge) in constrained_edges.iter().enumerate() {
        #[cfg(feature = "debug_context")]
        {
            let force_end = debug_context.advance_step();
            if force_end {
                break;
            }
        }

        #[cfg(feature = "progress_log")]
        {
            if _edge_index % ((constrained_edges.len() / 50) + 1) == 0 {
                let progress = 100. * _edge_index as f32 / constrained_edges.len() as f32;
                info!(
                    "Constraints progress, step n°{}, {}%: {}/{}",
                    debug_context.current_step,
                    progress,
                    _edge_index,
                    constrained_edges.len()
                );
            }
        }

        // Transform the edge to use the merged vertices indices
        let constrained_edge = Edge::new(
            vertex_merge_mapping[constrained_edge.from as usize],
            vertex_merge_mapping[constrained_edge.to as usize],
        );

        // Check for duplicate edges
        if constrained_edges_set.insert(constrained_edge) == false {
            // Skipping duplicate constrained edge
            continue;
        }

        // Stores all of the edges that cross the constrained edge
        let constrained_edge_vertices = &constrained_edge.to_vertices(vertices);
        register_intersected_edges(
            triangles,
            vertices,
            min_container_vertex_id,
            constrained_edge,
            constrained_edge_vertices,
            &vertex_to_triangle,
            &mut intersections,
            #[cfg(feature = "debug_context")]
            debug_context,
        );
        if intersections.is_empty() {
            // Skipping constrained edge already in triangulation
            continue;
        }

        // Remove intersecting edges
        remove_crossed_edges(
            triangles,
            vertices,
            min_container_vertex_id,
            &mut vertex_to_triangle,
            constrained_edge_vertices,
            &mut intersections,
            &mut new_diagonals_created,
            #[cfg(feature = "debug_context")]
            debug_context,
        );

        // Restore Delaunay triangulation
        restore_delaunay_triangulation_constrained(
            triangles,
            vertices,
            min_container_vertex_id,
            &mut vertex_to_triangle,
            constrained_edge,
            &mut new_diagonals_created,
            #[cfg(feature = "debug_context")]
            debug_context,
        );
    }
    constrained_edges_set
}

#[derive(Debug)]
enum EdgeFirstIntersection {
    /// The constrained edge is already present in the triangulation
    EdgeAlreadyInTriangulation,
    /// Intersection found
    Intersection(EdgeData),
    /// No intersection found
    NotFound,
}

fn edge_and_constrained_edge_intersection(
    edge: &Edge,
    edge_vertices: &EdgeVertices,
    constrained_edge_verts: &EdgeVertices,
    min_container_vertex_id: VertexId,
    #[cfg(feature = "debug_context")] _debug_context: &mut DebugContext,
) -> EdgesIntersectionResult {
    // Handle infinite vertices
    // We know that constrained_edge_verts are finite here
    let mut infinite_verts = Vec::new();
    if is_infinite(edge.from, min_container_vertex_id) {
        infinite_verts.push((
            (edge.from - min_container_vertex_id) as TriangleVertexIndex,
            edge_vertices.1,
        ));
    }
    if is_infinite(edge.to, min_container_vertex_id) {
        infinite_verts.push((
            (edge.to - min_container_vertex_id) as TriangleVertexIndex,
            edge_vertices.0,
        ));
    }

    if infinite_verts.is_empty() {
        egdes_intersect(&constrained_edge_verts, &edge_vertices)
    } else if infinite_verts.len() == 1 {
        // TODO Cold function  ?
        // Test intersection between edge and infinite edge segment
        let (infinite_vert_id, finite_vert) = infinite_verts[0];
        let extrapolated_segment = edge_from_semi_infinite_edge(finite_vert, infinite_vert_id);
        egdes_intersect(&constrained_edge_verts, &extrapolated_segment)
    } else {
        EdgesIntersectionResult::None
    }
}

fn loop_around_vertex_and_search_intersection(
    triangles: &Triangles,
    vertices: &Vec<Vertex>,
    min_container_vertex_id: VertexId,
    from_triangle_id: TriangleId,
    constrained_edge: Edge,
    constrained_edge_verts: &EdgeVertices,
    next_edge_index: fn(TriangleVertexIndex) -> TriangleEdgeIndex,
    #[cfg(feature = "debug_context")] _debug_context: &mut DebugContext,
) -> EdgeFirstIntersection {
    let mut triangle_id = from_triangle_id;

    // We search by circling around the first vertex of the constrained edge
    // We use `triangles.len()` as an upper bound on the number of triangles around a vertex
    for _i in 0..triangles.count() {
        if _i > 0 && triangle_id == from_triangle_id {
            error!("Internal error, aborting cycle when searching for the first intersection of an edge");
            break;
        }

        let triangle = triangles.get(triangle_id);

        if triangle.verts.contains(&constrained_edge.to) {
            return EdgeFirstIntersection::EdgeAlreadyInTriangulation;
        }

        let vert_index = triangle.vertex_index(constrained_edge.from);
        let edge_index = opposite_edge_index(vert_index);
        let edge = triangle.edge(edge_index);
        let edge_vertices = edge.to_vertices(vertices);

        let intersection = edge_and_constrained_edge_intersection(
            &edge,
            &edge_vertices,
            constrained_edge_verts,
            min_container_vertex_id,
            #[cfg(feature = "debug_context")]
            _debug_context,
        );
        if intersection == EdgesIntersectionResult::Crossing {
            let neighbor_triangle = triangle.neighbor(edge_index);
            if neighbor_triangle.exists() {
                return EdgeFirstIntersection::Intersection(EdgeData {
                    from_triangle_id: triangle_id,
                    to_triangle_id: neighbor_triangle.id,
                    edge,
                });
            } else {
                error!("Internal error, loop_around_vertex_and_search_intersection found a triangle with a crossed edge but no neihgbor, triangle {} {:?}, starting_vertex {}", triangle_id, triangle, constrained_edge.from);
                break;
            }
        }

        let neighbor = triangle.neighbor(next_edge_index(vert_index));
        if neighbor.exists() {
            triangle_id = neighbor.id;
        } else {
            break;
        }
    }
    EdgeFirstIntersection::NotFound
}

fn search_first_interstected_quad(
    triangles: &Triangles,
    vertices: &Vec<Vertex>,
    min_container_vertex_id: VertexId,
    constrained_edge: Edge,
    constrained_edge_vertices: &EdgeVertices,
    start_triangle: TriangleId,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) -> EdgeFirstIntersection {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "search_first_interstected_quad").entered();

    // Search clockwise
    let intersection = loop_around_vertex_and_search_intersection(
        triangles,
        vertices,
        min_container_vertex_id,
        start_triangle,
        constrained_edge,
        constrained_edge_vertices,
        outermost_clockwise_edge_index_around,
        #[cfg(feature = "debug_context")]
        debug_context,
    );
    match &intersection {
        EdgeFirstIntersection::NotFound => (),
        _ => return intersection,
    }

    // Get the counter-clockwise neighbor of the starting triangle
    let triangle = triangles.get(start_triangle);
    // Not having the constrained edge first vertex in the starting triangle is not possible
    let vert_index = triangle.vertex_index(constrained_edge.from);
    let edge_index = outermost_counter_clockwise_edge_index_around(vert_index);
    // It is impossible to have no neighbor here. It would mean that there is no triangle with the constrained edge intersection
    let neighbor_triangle = triangle.neighbor(edge_index).id;

    // Search counterclockwise
    let intersection = loop_around_vertex_and_search_intersection(
        triangles,
        vertices,
        min_container_vertex_id,
        neighbor_triangle,
        constrained_edge,
        constrained_edge_vertices,
        outermost_counter_clockwise_edge_index_around,
        #[cfg(feature = "debug_context")]
        debug_context,
    );
    match &intersection {
        EdgeFirstIntersection::NotFound => {
            // TODO Return an internal error ?
            error!("Internal error, search_first_interstected_quad found no triangle with the constrained edge intersection, starting_vertex {}", constrained_edge.from);
        }
        _ => (),
    }
    return intersection;
}

fn register_intersected_edges(
    triangles: &Triangles,
    vertices: &Vec<Vertex>,
    min_container_vertex_id: VertexId,
    constrained_edge: Edge,
    constrained_edge_vertices: &EdgeVertices,
    vertex_to_triangle: &Vec<TriangleId>,
    intersections: &mut VecDeque<EdgeData>,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "register_intersected_edges").entered();

    let search_result = search_first_interstected_quad(
        triangles,
        vertices,
        min_container_vertex_id,
        constrained_edge,
        constrained_edge_vertices,
        vertex_to_triangle[constrained_edge.from as usize],
        #[cfg(feature = "debug_context")]
        debug_context,
    );

    let mut intersection = match search_result {
        EdgeFirstIntersection::Intersection(edge_data) => {
            intersections.push_back(edge_data.clone());
            edge_data
        }
        _ => return,
    };

    // Search all the edges intersected by this constrained edge
    // We use `triangles.len()` as an upper bound on the number of triangles intersecting a constrained edge
    for _ in 0..triangles.count() {
        let triangle_id = intersection.to_triangle_id;
        let triangle = triangles.get(triangle_id);

        if triangle.verts.contains(&constrained_edge.to) {
            // Stop if we reached the end of the constrained edge
            break;
        }

        match get_next_triangle_edge_intersection(
            min_container_vertex_id,
            &triangle.to_vertices(vertices),
            constrained_edge_vertices,
            triangle,
            triangle_id,
            &intersection.edge,
            #[cfg(feature = "debug_context")]
            debug_context,
        ) {
            Some(next_intersection) => {
                intersection = next_intersection;
                intersections.push_back(intersection.clone());
            }
            None => {
                // Not possible, the edge must intersects at least one other side of the triangle here
                error!("Internal error, triangle {} {:?} should be intersecting the constrained edge {:?}", triangle_id, triangle, constrained_edge);
                break;
            }
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct EdgeData {
    from_triangle_id: TriangleId,
    to_triangle_id: TriangleId,
    edge: Edge,
}
impl EdgeData {
    fn new(from_triangle_id: TriangleId, to_triangle_id: TriangleId, crossed_edge: Edge) -> Self {
        Self {
            from_triangle_id,
            to_triangle_id,
            edge: crossed_edge,
        }
    }

    #[inline]
    fn from(&self) -> TriangleId {
        self.from_triangle_id
    }
    #[inline]
    fn to(&self) -> TriangleId {
        self.to_triangle_id
    }
}

/// Represents an edge between two triangles
///
/// ```text
///               q3
///              /  \
///            /  Tt  \
///          /          \
///      q1=e1 -------- q2=e2
///          \          /
///            \  Tf  /
///              \  /
///               q4
/// ```
///
/// where:
/// Tf: from triangle
/// Tt: to triangle
/// e1, e2: edge vertices
/// q1, q2, q3, q4: quad coords
impl EdgeData {
    pub fn to_quad(&self, triangles: &Triangles) -> Quad {
        Quad::new([
            self.edge.from,
            self.edge.to,
            triangles
                .get(self.to_triangle_id)
                .get_opposite_vertex_id(&self.edge),
            triangles
                .get(self.from_triangle_id)
                .get_opposite_vertex_id(&self.edge),
        ])
    }
}

fn get_next_triangle_edge_intersection(
    min_container_vertex_id: VertexId,
    triangle_verts: &TriangleVertices,
    constrained_edge: &EdgeVertices,
    triangle: &TriangleData,
    triangle_id: TriangleId,
    from_crossed_edge: &Edge,
    #[cfg(feature = "debug_context")] _debug_context: &mut DebugContext,
) -> Option<EdgeData> {
    for (edge_index, (edge_vertices, edge)) in vec![
        (triangle_verts.0, triangle_verts.1),
        (triangle_verts.1, triangle_verts.2),
        (triangle_verts.2, triangle_verts.0),
    ]
    .iter()
    .zip(triangle.edges())
    .enumerate()
    {
        if edge.undirected_equals(&from_crossed_edge) {
            continue;
        }
        // TODO Doc: Add a note about why we only care about the `Crossing` result
        if edge_and_constrained_edge_intersection(
            &edge,
            &edge_vertices,
            constrained_edge,
            min_container_vertex_id,
            #[cfg(feature = "debug_context")]
            _debug_context,
        ) == EdgesIntersectionResult::Crossing
        {
            let neighbor_triangle = triangle.neighbors[edge_index];
            if neighbor_triangle.exists() {
                return Some(EdgeData {
                    // crossed_edge_index: edge_index,
                    from_triangle_id: triangle_id,
                    to_triangle_id: neighbor_triangle.id,
                    edge,
                });
            }
        }
    }
    None
}

/// Swap diagonals of a quad like in the following diagram:
///
/// ```text
///               q3
///              /  \
///      Ttl   /  Tt  \    Ttr
///          /          \
///      q1=c1 -------- q2=c2
///          \          /
///            \  Tf  /
///      Tfl     \  /      Tfr
///               q4
///
///             Becomes
///
///               q3
///              / | \
///      Ttl   /   |   \   Ttr
///          /     |     \
///      q1=c1  Tf | Tt  q2=c2
///          \     |     /
///            \   |   /
///      Tfl     \ | /     Tfr
///               q4
/// ```
///
/// where:
/// Tf: triangle from
/// Tt: triangle to
/// c1, c2: crossed edge
/// q1, q2, q3, q4: quad coords
pub(crate) fn swap_quad_diagonal(
    triangles: &mut Triangles,
    vertex_to_triangle: &mut Vec<TriangleId>,
    edges_data_collections: &mut [&mut VecDeque<EdgeData>],
    edge_data: &EdgeData,
    quad: &Quad,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "swap_quad_diagonal").entered();

    let from = edge_data.from();
    let to = edge_data.to();

    // TODO Design: Memorize from_crossed_edge_index in Intersection ?
    // `tt`: triangle_to
    // `tf`: triangle_from
    let tt_left_edge_index = triangles.get(to).opposite_edge_index_from_vertex(quad.v2());
    let tt_left_neighbor = triangles.get(to).neighbor(tt_left_edge_index);
    let tt_right_neighbor = triangles
        .get(to)
        .neighbor(next_clockwise_edge_index(tt_left_edge_index));

    let tf_left_edge_index = triangles
        .get(from)
        .opposite_edge_index_from_vertex(quad.v2());
    let tf_left_neighbor = triangles.get(from).neighbor(tf_left_edge_index);
    let tf_right_neighbor = triangles
        .get(from)
        .neighbor(next_counter_clockwise_edge_index(tf_left_edge_index));

    triangles.get_mut(from).verts = [quad.v3(), quad.v4(), quad.v1()];
    triangles.get_mut(to).verts = [quad.v3(), quad.v2(), quad.v4()];

    triangles.get_mut(from).neighbors = [to.into(), tf_left_neighbor, tt_left_neighbor];
    triangles.get_mut(to).neighbors = [tt_right_neighbor, tf_right_neighbor, from.into()];

    triangulation::update_triangle_neighbor(tt_left_neighbor, to.into(), from.into(), triangles);
    triangulation::update_triangle_neighbor(tf_right_neighbor, from.into(), to.into(), triangles);

    vertex_to_triangle[quad.v1() as usize] = from;
    vertex_to_triangle[quad.v2() as usize] = to;

    #[cfg(feature = "debug_context")]
    debug_context.push_snapshot(
        Phase::ConstrainedSwapQuadDiagonals,
        &triangles,
        &[from, to],
        &[tt_left_neighbor, tf_right_neighbor],
    );

    // Update in memory collectiosn that may contain now out of date data
    for edges_data_collection in edges_data_collections.iter_mut() {
        update_edges_data(
            *edges_data_collection,
            tt_left_neighbor,
            tf_right_neighbor,
            from,
            to,
        );
    }
}

fn update_edges_data(
    edges_data: &mut VecDeque<EdgeData>,
    tt_left_neighbor: Neighbor,
    tf_right_neighbor: Neighbor,
    from: TriangleId,
    to: TriangleId,
) {
    for edge_data in edges_data.iter_mut() {
        // We do not need to check for updates of triangle_from<->triangle_to pairs since there should only be one in the collection and we are currently treating it.
        if tt_left_neighbor.exists() {
            if edge_data.from_triangle_id == tt_left_neighbor.id && edge_data.to_triangle_id == to {
                edge_data.to_triangle_id = from;
            } else if edge_data.from_triangle_id == to
                && edge_data.to_triangle_id == tt_left_neighbor.id
            {
                edge_data.from_triangle_id = from;
                // TODO Design: inter.from_edge_index becomes EDGE_12
            }
        }
        if tf_right_neighbor.exists() {
            if edge_data.from_triangle_id == tf_right_neighbor.id
                && edge_data.to_triangle_id == from
            {
                edge_data.to_triangle_id = to;
            } else if edge_data.from_triangle_id == from
                && edge_data.to_triangle_id == tf_right_neighbor.id
            {
                edge_data.from_triangle_id = to;
                // TODO Design: inter.from_edge_index becomes EDGE_23
            }
        }
        // TODO Design: If we decide to store from_edge_index, we also need to update tfl<->tf and ttr<->tt pairs
    }
}

fn quad_diagonals_intersection(
    quad: &Quad,
    quad_vertices: &QuadVertices,
    min_container_vertex_id: VertexId,
) -> EdgesIntersectionResult {
    // TODO Share alloc ?
    let mut infinite_verts = Vec::new();
    if is_infinite(quad.v1(), min_container_vertex_id) {
        infinite_verts.push((
            (quad.v1() - min_container_vertex_id) as TriangleVertexIndex,
            quad_vertices.q2(),
            QUAD_1,
        ));
    }
    if is_infinite(quad.v2(), min_container_vertex_id) {
        infinite_verts.push((
            (quad.v2() - min_container_vertex_id) as TriangleVertexIndex,
            quad_vertices.q1(),
            QUAD_2,
        ));
    }
    if is_infinite(quad.v3(), min_container_vertex_id) {
        infinite_verts.push((
            (quad.v3() - min_container_vertex_id) as TriangleVertexIndex,
            quad_vertices.q4(),
            QUAD_3,
        ));
    }
    if is_infinite(quad.v4(), min_container_vertex_id) {
        infinite_verts.push((
            (quad.v4() - min_container_vertex_id) as TriangleVertexIndex,
            quad_vertices.q3(),
            QUAD_4,
        ));
    }

    if infinite_verts.is_empty() {
        quad_vertices.diagonals_intersection_test()
    }
    // For the cases with infite vertices, we know that a quad's diagonal has at least 1 finite vertex
    // TODO Those cases do not seem to be present in tested data sets. Try to have a dataset/test making use of it
    else if infinite_verts.len() == 1 {
        // TODO Cold function
        let (infinite_vert_index, finite_vertex, finite_vert_index) = infinite_verts[0];
        let other_diagonal_verts = ADJACENT_QUAD_VERTICES_INDEXES[finite_vert_index as usize];
        let other_diagonal_edge = (
            quad_vertices.verts[other_diagonal_verts[0] as usize],
            quad_vertices.verts[other_diagonal_verts[1] as usize],
        );

        egdes_intersect(
            &other_diagonal_edge,
            &edge_from_semi_infinite_edge(finite_vertex, infinite_vert_index),
        )
    } else {
        // TODO Cold function
        let (infinite_vert_id_0, finite_vert_0, _) = infinite_verts[0];
        let infinite_edge_segment_0 =
            edge_from_semi_infinite_edge(finite_vert_0, infinite_vert_id_0);

        let (infinite_vert_id_1, finite_vert_1, _) = infinite_verts[1];
        let infinite_edge_segment_1 =
            edge_from_semi_infinite_edge(finite_vert_1, infinite_vert_id_1);

        egdes_intersect(&infinite_edge_segment_0, &infinite_edge_segment_1)
    }
}

fn remove_crossed_edges(
    triangles: &mut Triangles,
    vertices: &Vec<Vertex>,
    min_container_vertex_id: VertexId,
    vertex_to_triangle: &mut Vec<TriangleId>,
    constrained_edge_vertices: &EdgeVertices,
    intersections: &mut VecDeque<EdgeData>,
    new_diagonals_created: &mut VecDeque<EdgeData>,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "remove_crossed_edges").entered();

    while let Some(intersection) = intersections.pop_front() {
        let quad = intersection.to_quad(triangles);
        let quad_vertices = quad.to_vertices(vertices);
        // If the quad is not convex (diagonals do not cross) or if an edge tip lie on the other edge,
        // we skip this edge and put it on the stack to re-process it later
        if quad_diagonals_intersection(&quad, &quad_vertices, min_container_vertex_id)
            != EdgesIntersectionResult::Crossing
        {
            intersections.push_back(intersection);
        }
        // Swap the diagonal of this strictly convex quadrilateral if the two diagonals cross normaly
        else {
            swap_quad_diagonal(
                triangles,
                vertex_to_triangle,
                &mut [intersections, new_diagonals_created],
                &intersection,
                &quad,
                #[cfg(feature = "debug_context")]
                debug_context,
            );
            // TODO Clean: this could be returned by swap_quad_diagonal ?
            let edge_data = EdgeData::new(
                intersection.from_triangle_id,
                intersection.to_triangle_id,
                Edge::new(quad.v3(), quad.v4()),
            );
            // If the new diagonal still intersects the constrained edge,
            // then place it on the list of intersecting edges
            if edge_and_constrained_edge_intersection(
                &edge_data.edge,
                &(quad_vertices.q3(), quad_vertices.q4()),
                constrained_edge_vertices,
                min_container_vertex_id,
                #[cfg(feature = "debug_context")]
                debug_context,
            ) == EdgesIntersectionResult::Crossing
            {
                intersections.push_back(edge_data);
            } else {
                new_diagonals_created.push_back(edge_data);
            }
        }
    }
}

fn restore_delaunay_triangulation_constrained(
    triangles: &mut Triangles,
    vertices: &Vec<Vertex>,
    min_container_vertex_id: VertexId,
    vertex_to_triangle: &mut Vec<TriangleId>,
    constrained_edge: Edge,
    new_diagonals_created: &mut VecDeque<EdgeData>,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "restore_delaunay_triangulation_constrained").entered();

    let mut swaps_history = HashSet::new();
    while let Some(new_edge) = new_diagonals_created.pop_front() {
        if new_edge.edge.undirected_equals(&constrained_edge) {
            continue;
        } else {
            let quad = new_edge.to_quad(triangles);

            let should_swap = should_swap_diagonals(&quad, vertices, min_container_vertex_id);

            // Use a small swap history to detect infinite swap loops & cycles.
            // Ideally we would like to prevent them instead of detecting them but this seems to be harder than it looks.
            // Simply being a bit more more strict on the circumcircle test causes issues on some datasets.
            if should_swap && swaps_history.insert(new_edge.clone()) {
                swap_quad_diagonal(
                    triangles,
                    vertex_to_triangle,
                    &mut [new_diagonals_created],
                    &new_edge,
                    &quad,
                    #[cfg(feature = "debug_context")]
                    debug_context,
                );

                new_diagonals_created.push_back(EdgeData::new(
                    new_edge.from_triangle_id,
                    new_edge.to_triangle_id,
                    Edge::new(quad.v3(), quad.v4()),
                ));
            }
        }
    }
}

///////////////////////////////////////////////////////////
///                                                     ///
///                        Tests                        ///
///                                                     ///
///////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use std::collections::VecDeque;

    use crate::types::{Edge, Neighbor, Quad, TriangleData, TriangleId, Triangles, Vertex};

    use super::{swap_quad_diagonal, EdgeData};

    #[cfg(feature = "debug_context")]
    use crate::debug::{DebugConfiguration, DebugContext};

    #[test]
    fn swap_quad_diag_constrained() {
        let mut vertices = Vec::<Vertex>::new();
        vertices.push(Vertex::new(0.5, 3.));
        vertices.push(Vertex::new(-2., -2.));
        vertices.push(Vertex::new(1., -4.));
        vertices.push(Vertex::new(3., -2.));

        let triangle_1 = TriangleData {
            verts: [0, 1, 3],
            neighbors: [Neighbor::NONE, Neighbor::new(1), Neighbor::NONE],
        };

        let triangle_2 = TriangleData {
            verts: [1, 2, 3],
            neighbors: [Neighbor::NONE, Neighbor::NONE, Neighbor::new(0)],
        };

        let mut triangles = Triangles::new();
        triangles.buffer.extend([triangle_1, triangle_2]);

        let mut vertex_to_triangle = vec![0; vertices.len()];
        for (index, triangle) in triangles.buffer().iter().enumerate() {
            vertex_to_triangle[triangle.v1() as usize] = index as TriangleId;
            vertex_to_triangle[triangle.v2() as usize] = index as TriangleId;
            vertex_to_triangle[triangle.v3() as usize] = index as TriangleId;
        }

        let mut new_diagonals_created = VecDeque::new();

        let intersection = EdgeData::new(0, 1, Edge::new(1, 3));

        let mut edge_data_collection = [&mut new_diagonals_created];

        let quad = Quad::new([1, 3, 0, 2]);

        #[cfg(feature = "debug_context")]
        let mut debug_context = DebugContext::new(DebugConfiguration::default(), 0., 0., 0.);
        swap_quad_diagonal(
            &mut triangles,
            &mut vertex_to_triangle,
            &mut edge_data_collection,
            &intersection,
            &quad,
            #[cfg(feature = "debug_context")]
            &mut debug_context,
        );

        assert_eq!(2, triangles.count());
        assert_eq!([2, 3, 0], triangles.get(0).verts);
        assert_eq!([2, 0, 1], triangles.get(1).verts);
    }
}
