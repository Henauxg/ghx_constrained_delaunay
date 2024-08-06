use std::collections::VecDeque;

use hashbrown::HashSet;

use crate::infinite::{
    collect_infinite_quad_vertices, edge_from_semi_infinite_edge, infinite_vertex_local_quad_index,
    is_finite, is_infinite, is_vertex_in_half_plane_1, quad_diagonals_intersection_1_infinite,
    quad_diagonals_intersection_2_infinite,
};
use crate::triangulation::{
    normalize_vertices_coordinates, update_neighbor_neighbor, Triangulation, TriangulationError,
    DEFAULT_BIN_VERTEX_DENSITY_POWER,
};
use crate::types::{
    next_counter_clockwise_edge_index, TriangleEdgeIndex, Triangles, Vertex, Vertex2d, Vertex3d,
    VERT_1, VERT_2,
};
use crate::utils::{egdes_intersect, is_vertex_in_triangle_circumcircle, EdgesIntersectionResult};

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
        TriangleData, TriangleId, TriangleVertexIndex, VertexId,
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

/// Same as [constrained_triangulation_from_2d_vertices] but input vertices are in 3d and will be transformed to 2d before the triangulation.
///
/// Additional requirements:
/// - All the vertices are expected to belong to the same 2d plane, with the provided `plane_normal`.
/// - `plane_normal` must be normalized
pub fn constrained_triangulation_from_3d_planar_vertices<T: Vertex3d>(
    vertices: &[T],
    plane_normal: T,
    constrained_edges: &[Edge],
    config: ConstrainedTriangulationConfiguration,
) -> Result<Triangulation, TriangulationError> {
    if vertices.len() < 3 {
        return Ok(Triangulation::default());
    }

    let mut planar_vertices =
        triangulation::transform_to_2d_planar_coordinate_system(vertices, plane_normal);
    constrained_triangulation_from_2d_vertices(&mut planar_vertices, constrained_edges, config)
}

/// Creates a constrained Delaunay triangulation with the input vertices and edges.
///
/// Vertices that are identical (or extremely close to one another) will me "merged" together. This means that only one of those will appear in the final triangulation.
///
/// Vertices requirements:
/// - Vertices are expected to be valid floating points values. You can use [crate::utils::validate_vertices] to check your vertices beforehand for NaN or infinity.
///
/// Edges requirements:
/// - Two constrained egdes should NOT cross
/// - A constrained edge (two vertices) should NOT cross another vertex of the input data
/// - Vertices should contain valid floating point values
///
/// Constrained edges orientation:
/// - Contours (set of constrained edges) forming a CW cycle indicate a domain to triangulate.
/// - Contours forming a CCW cycle indicate a domain that should be discarded from the triangulation.
///
/// Example figure below:
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
pub fn constrained_triangulation_from_2d_vertices<T: Vertex2d>(
    vertices: &[T],
    constrained_edges: &[Edge],
    config: ConstrainedTriangulationConfiguration,
) -> Result<Triangulation, TriangulationError> {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "constrained_triangulation_from_2d_vertices").entered();

    if vertices.len() < 3 {
        return Ok(Triangulation::default());
    }

    // Uniformly scale the coordinates of the points so that they all lie between 0 and 1.
    let (normalized_vertices, _scale_factor, _x_min, _y_min) =
        normalize_vertices_coordinates(vertices);

    #[cfg(feature = "debug_context")]
    let mut debug_context =
        DebugContext::new(config.debug_config.clone(), _scale_factor, _x_min, _y_min);

    let mut vertex_merge_mapping = Some((0..vertices.len() as VertexId).collect());
    let mut triangles = wrap_and_triangulate_2d_normalized_vertices(
        &normalized_vertices,
        config.bin_vertex_density_power,
        &mut vertex_merge_mapping,
        #[cfg(feature = "debug_context")]
        &mut debug_context,
    )?;

    #[cfg(feature = "debug_context")]
    debug_context.push_snapshot(Phase::BeforeConstraints, &triangles, &[], &[]);

    // TODO Debug: consider adding debug support to this function
    let constrained_edges_set = apply_constraints(
        &normalized_vertices,
        &mut triangles,
        constrained_edges,
        vertex_merge_mapping.unwrap(),
        #[cfg(feature = "debug_context")]
        &mut debug_context,
    )?;

    #[cfg(feature = "debug_context")]
    debug_context.push_snapshot(Phase::AfterConstraints, &triangles, &[], &[]);

    let vert_indices = cdt_filter_triangles(
        &triangles,
        constrained_edges_set,
        #[cfg(feature = "debug_context")]
        &mut debug_context,
    );

    Ok(Triangulation {
        triangles: vert_indices,
        #[cfg(feature = "debug_context")]
        debug_context,
    })
}

/// Filter triangles that should be removed due to the input domains constraints or due to being part of
/// the container triangle.
fn cdt_filter_triangles(
    triangles: &Triangles,
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

            if triangle.is_finite() {
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
                    if t.is_finite() {
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
    vertices: &[Vertex],
    triangles: &mut Triangles,
    constrained_edges: &[Edge],
    vertex_merge_mapping: Vec<VertexId>,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) -> Result<HashSet<Edge>, TriangulationError> {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "apply_constraints").entered();

    let mut constrained_edges_set = HashSet::with_capacity(constrained_edges.len());
    // Map each finite vertex to one of the triangles (the last) that contains it
    let mut vertex_to_triangle = vec![0; vertices.len()];
    for (index, t) in triangles.buffer().iter().enumerate() {
        if is_finite(t.v1()) {
            vertex_to_triangle[t.v1() as usize] = index as TriangleId;
        }
        if is_finite(t.v2()) {
            vertex_to_triangle[t.v2() as usize] = index as TriangleId;
        }
        if is_finite(t.v3()) {
            vertex_to_triangle[t.v3() as usize] = index as TriangleId;
        }
    }

    // Those buffers are used by all calls to `register_intersected_edges`, `remove_crossed_edges` and `restore_delaunay_triangulation_constrained`.
    // We create them here to share the allocation between all those calls as an optimization.
    let mut intersections = VecDeque::new();
    let mut new_diagonals_created = VecDeque::new();
    // Shared alloc for `restore_delaunay_triangulation_constrained`
    let mut swaps_history = HashSet::new();

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

        // Check for invalid or duplicate edges
        // We will use this Set built iteratively for domains removal at the end
        if constrained_edge.from == constrained_edge.to
            || constrained_edges_set.insert(constrained_edge) == false
        {
            // Skipping
            continue;
        }

        let constrained_edge_vertices = &constrained_edge.to_vertices(vertices);
        // 'intersections' will store all the edges crossed by the constrained edge
        register_intersected_edges(
            triangles,
            vertices,
            constrained_edge,
            constrained_edge_vertices,
            &vertex_to_triangle,
            &mut intersections,
            #[cfg(feature = "debug_context")]
            debug_context,
        )?;
        if intersections.is_empty() {
            // Skipping, constrained edge is already an edge of the current triangulation
            continue;
        }

        // Remove intersecting edges
        remove_crossed_edges(
            triangles,
            vertices,
            &mut vertex_to_triangle,
            constrained_edge_vertices,
            &mut intersections,
            &mut new_diagonals_created,
            #[cfg(feature = "debug_context")]
            debug_context,
        )?;

        // Restore Delaunay triangulation
        restore_delaunay_triangulation_constrained(
            triangles,
            vertices,
            &mut vertex_to_triangle,
            constrained_edge,
            &mut new_diagonals_created,
            &mut swaps_history,
            #[cfg(feature = "debug_context")]
            debug_context,
        );
    }
    Ok(constrained_edges_set)
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
    vertices: &[Vertex],
    edge: &Edge,
    constrained_edge_verts: &EdgeVertices,
    #[cfg(feature = "debug_context")] _debug_context: &mut DebugContext,
) -> EdgesIntersectionResult {
    // Handle infinite vertices
    // We know that constrained_edge_verts are finite here
    let mut infinite_verts = Vec::new();
    if is_infinite(edge.from) {
        infinite_verts.push((infinite_vertex_local_quad_index(edge.from), edge.to));
    }
    if is_infinite(edge.to) {
        infinite_verts.push((infinite_vertex_local_quad_index(edge.to), edge.from));
    }

    if infinite_verts.is_empty() {
        egdes_intersect(&constrained_edge_verts, &edge.to_vertices(vertices))
    } else if infinite_verts.len() == 1 {
        // TODO Cold function  ?
        // Test intersection between edge and infinite edge segment
        let (infinite_vertex_index, finite_vert_id) = infinite_verts[0];
        let extrapolated_segment =
            edge_from_semi_infinite_edge(vertices[finite_vert_id as usize], infinite_vertex_index);
        egdes_intersect(&constrained_edge_verts, &extrapolated_segment)
    } else {
        EdgesIntersectionResult::None
    }
}

fn loop_around_vertex_and_search_intersection(
    triangles: &Triangles,
    vertices: &[Vertex],
    from_triangle_id: TriangleId,
    constrained_edge: Edge,
    constrained_edge_verts: &EdgeVertices,
    next_edge_index: fn(TriangleVertexIndex) -> TriangleEdgeIndex,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) -> Result<EdgeFirstIntersection, TriangulationError> {
    let mut triangle_id = from_triangle_id;

    // We search by circling around the first vertex of the constrained edge
    // We use `triangles.len()` as an upper bound on the number of triangles around a vertex
    for _i in 0..triangles.count() {
        if _i > 0 && triangle_id == from_triangle_id {
            return Err(TriangulationError::new(
                format!("Internal error, aborting cycle when searching for the first intersection of an edge"),
                #[cfg(feature = "debug_context")]
                debug_context,
            ));
        }

        let triangle = triangles.get(triangle_id);

        if triangle.verts.contains(&constrained_edge.to) {
            return Ok(EdgeFirstIntersection::EdgeAlreadyInTriangulation);
        }

        let vert_index = triangle.vertex_index(constrained_edge.from);
        let edge_index = opposite_edge_index(vert_index);
        let edge = triangle.edge(edge_index);

        let intersection = edge_and_constrained_edge_intersection(
            vertices,
            &edge,
            constrained_edge_verts,
            #[cfg(feature = "debug_context")]
            debug_context,
        );
        if intersection == EdgesIntersectionResult::Crossing {
            let neighbor_triangle = triangle.neighbor(edge_index);
            if neighbor_triangle.exists() {
                return Ok(EdgeFirstIntersection::Intersection(EdgeData {
                    from_triangle_id: triangle_id,
                    to_triangle_id: neighbor_triangle.id,
                    edge,
                }));
            } else {
                return Err(TriangulationError::new(
                    format!("Internal error, loop_around_vertex_and_search_intersection found a triangle with a crossed edge but no neihgbor, triangle {} {:?}, starting_vertex {}", triangle_id, triangle, constrained_edge.from),
                    #[cfg(feature = "debug_context")]
                    debug_context,
                ));
            }
        }

        let neighbor = triangle.neighbor(next_edge_index(vert_index));
        if neighbor.exists() {
            triangle_id = neighbor.id;
        } else {
            break;
        }
    }
    Ok(EdgeFirstIntersection::NotFound)
}

fn search_first_interstected_quad(
    triangles: &Triangles,
    vertices: &[Vertex],
    constrained_edge: Edge,
    constrained_edge_vertices: &EdgeVertices,
    start_triangle: TriangleId,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) -> Result<EdgeFirstIntersection, TriangulationError> {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "search_first_interstected_quad").entered();

    // Search clockwise. We may not find it.
    let intersection = loop_around_vertex_and_search_intersection(
        triangles,
        vertices,
        start_triangle,
        constrained_edge,
        constrained_edge_vertices,
        outermost_clockwise_edge_index_around,
        #[cfg(feature = "debug_context")]
        debug_context,
    )?;
    match &intersection {
        EdgeFirstIntersection::NotFound => (),
        _ => return Ok(intersection),
    }

    // Get the counter-clockwise neighbor of the starting triangle
    let triangle = triangles.get(start_triangle);
    // Not having the constrained edge first vertex in the starting triangle is not possible
    let vert_index = triangle.vertex_index(constrained_edge.from);
    let edge_index = outermost_counter_clockwise_edge_index_around(vert_index);
    // It is impossible to have no neighbor here. It would mean that there is no triangle with the constrained edge intersection
    let neighbor_triangle = triangle.neighbor(edge_index).id;

    // Search counterclockwise. We must find it.
    let intersection = loop_around_vertex_and_search_intersection(
        triangles,
        vertices,
        neighbor_triangle,
        constrained_edge,
        constrained_edge_vertices,
        outermost_counter_clockwise_edge_index_around,
        #[cfg(feature = "debug_context")]
        debug_context,
    )?;
    Ok(intersection)
}

fn register_intersected_edges(
    triangles: &Triangles,
    vertices: &[Vertex],
    constrained_edge: Edge,
    constrained_edge_vertices: &EdgeVertices,
    vertex_to_triangle: &[TriangleId],
    intersections: &mut VecDeque<EdgeData>,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) -> Result<(), TriangulationError> {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "register_intersected_edges").entered();

    let search_result = search_first_interstected_quad(
        triangles,
        vertices,
        constrained_edge,
        constrained_edge_vertices,
        vertex_to_triangle[constrained_edge.from as usize],
        #[cfg(feature = "debug_context")]
        debug_context,
    )?;

    let mut intersection = match search_result {
        EdgeFirstIntersection::Intersection(edge_data) => {
            intersections.push_back(edge_data.clone());
            edge_data
        }
        EdgeFirstIntersection::EdgeAlreadyInTriangulation => return Ok(()),
        EdgeFirstIntersection::NotFound => return Err(TriangulationError::new(
            format!("Internal error, the search for the first intersected quad found no triangle with the constrained edge intersection, from starting_vertex {}", constrained_edge.from
           ),
           #[cfg(feature = "debug_context")]
            debug_context,
        )),
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
            vertices,
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
                return Err(TriangulationError::new(
                    format!(
                    "Internal error, triangle {} {:?} should have at least 1 intersection with the constrained edge {:?}", triangle_id, triangle, constrained_edge
                    ),
                    #[cfg(feature = "debug_context")]
                    debug_context,
                ));
            }
        }
    }
    Ok(())
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
    vertices: &[Vertex],
    constrained_edge: &EdgeVertices,
    triangle: &TriangleData,
    triangle_id: TriangleId,
    from_crossed_edge: &Edge,
    #[cfg(feature = "debug_context")] _debug_context: &mut DebugContext,
) -> Option<EdgeData> {
    for (edge_index, edge) in triangle.edges().iter().enumerate() {
        if edge.undirected_equals(&from_crossed_edge) {
            continue;
        }
        // TODO Doc: Add a note about why we only care about the `Crossing` result
        if edge_and_constrained_edge_intersection(
            vertices,
            &edge,
            constrained_edge,
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
                    edge: *edge,
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

    update_neighbor_neighbor(tt_left_neighbor, to.into(), from.into(), triangles);
    update_neighbor_neighbor(tf_right_neighbor, from.into(), to.into(), triangles);

    if is_finite(quad.v1()) {
        vertex_to_triangle[quad.v1() as usize] = from;
    }
    if is_finite(quad.v2()) {
        vertex_to_triangle[quad.v2() as usize] = to;
    }

    #[cfg(feature = "debug_context")]
    debug_context.push_snapshot(
        Phase::ConstrainedSwapQuadDiagonals,
        &triangles,
        &[from, to],
        &[tt_left_neighbor, tf_right_neighbor],
    );

    // Update in memory collections that may contain now out of date data
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

fn quad_diagonals_intersection(vertices: &[Vertex], quad: &Quad) -> EdgesIntersectionResult {
    let infinite_verts = collect_infinite_quad_vertices(&quad.verts);

    if infinite_verts.is_empty() {
        quad.to_vertices(vertices).diagonals_intersection_test()
    }
    // For the cases with infinite vertices, we know that a quad's diagonal has at least 1 finite vertex.
    // TODO Those cases do not seem to be present in the tested data sets. Try to have a dataset/test making use of it
    else if infinite_verts.len() == 1 {
        quad_diagonals_intersection_1_infinite(vertices, quad, infinite_verts[0])
    } else {
        quad_diagonals_intersection_2_infinite(vertices, quad, infinite_verts[0], infinite_verts[1])
    }
}

fn remove_crossed_edges(
    triangles: &mut Triangles,
    vertices: &[Vertex],
    vertex_to_triangle: &mut Vec<TriangleId>,
    constrained_edge_vertices: &EdgeVertices,
    intersections: &mut VecDeque<EdgeData>,
    new_diagonals_created: &mut VecDeque<EdgeData>,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) -> Result<(), TriangulationError> {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "remove_crossed_edges").entered();

    // TODO Optimization: could try to store quads directly in the intersections queue (small gain to avoid multiple conversions when some edges are pushed back for re-processing)
    let mut non_convex_counter = 0;
    while let Some(intersection) = intersections.pop_front() {
        // This is a safety against a possible infinite loop.
        if non_convex_counter > intersections.len() {
            return Err(TriangulationError::new(
                format!("Failed to remove crossed edges for constrained egde {:?}, no convex quads in the intersection queue", constrained_edge_vertices),
                #[cfg(feature = "debug_context")]
                debug_context,
            ));
        }
        let quad = intersection.to_quad(triangles);
        // If the quad is not convex (diagonals do not cross) or if an edge tip lie on the other edge,
        // we skip this edge and put it on the stack to re-process it later
        if quad_diagonals_intersection(vertices, &quad) != EdgesIntersectionResult::Crossing {
            intersections.push_back(intersection);
            non_convex_counter += 1;
        }
        // Swap the diagonal of this strictly convex quadrilateral if the two diagonals cross normaly
        else {
            non_convex_counter = 0;
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
                vertices,
                &edge_data.edge,
                constrained_edge_vertices,
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

    Ok(())
}

fn restore_delaunay_triangulation_constrained(
    triangles: &mut Triangles,
    vertices: &[Vertex],
    vertex_to_triangle: &mut Vec<TriangleId>,
    constrained_edge: Edge,
    new_diagonals_created: &mut VecDeque<EdgeData>,
    swaps_history: &mut HashSet<(TriangleId, TriangleId)>,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "restore_delaunay_triangulation_constrained").entered();

    swaps_history.clear();
    while let Some(new_edge) = new_diagonals_created.pop_front() {
        if new_edge.edge.undirected_equals(&constrained_edge) {
            continue;
        } else {
            let quad = new_edge.to_quad(triangles);

            let should_swap = constrained_should_swap_diagonals(&quad, vertices);

            // Use a small swap history to detect infinite swap loops & cycles.
            // Ideally we would like to prevent them instead of detecting them but this seems to be harder than it looks.
            // Simply being a bit more more strict on the circumcircle test causes issues on some datasets.
            if should_swap && swaps_history.insert((new_edge.from(), new_edge.to())) {
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

/// This is almost the exact copy of [triangulation::should_swap_diagonals] but this needs to also check q4 for infinite.
/// TODO Optimization: If we could ensure that we build all [EdgeData] with a finite q4 we could forego this chekck here too and use the same implem (and optimize a check away by the same occasion)
#[inline(always)]
pub fn constrained_should_swap_diagonals(quad: &Quad, vertices: &[Vertex]) -> bool {
    #[cfg(feature = "more_profile_traces")]
    let _span = span!(Level::TRACE, "should_swap_diagonals").entered();

    // Get rid of this case early: an infinite vertex cannot be in the circumcircle of the other triangle
    // Note: infinite q4 is not possible in DT (we always place a finite vertex in q4). This can only occur in CDT
    if is_infinite(quad.v3()) || is_infinite(quad.v4()) {
        return false;
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
