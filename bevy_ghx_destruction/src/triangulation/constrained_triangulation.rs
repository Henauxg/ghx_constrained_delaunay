use std::collections::VecDeque;

use bevy::{
    log::error,
    math::{Vec2, Vec3A},
    utils::hashbrown::HashSet,
};

use crate::utils::{egdes_intersect, is_vertex_in_triangle_circumcircle, EdgesIntersectionResult};

use super::{
    next_clockwise_edge_index, next_clockwise_edge_index_around, next_counter_clockwise_edge_index,
    next_counter_clockwise_edge_index_around, opposite_edge_index,
    triangulation::{self, wrap_and_triangulate_2d_vertices},
    Edge, EdgeVertices, Neighbor, Quad, TriangleData, TriangleId, TriangleVertexIndex,
    TriangleVertices, VertexId,
};

/// - `plane_normal` **MUST** be normalized
/// - `vertices` **MUST** all belong to a 3d plane
/// - `constrained_edges` **MUST**:
///     - be oriented: vi -> vj
///     - not contain edges with v0==v1
pub fn constrained_triangulation_from_3d_planar_vertices(
    vertices: &Vec<[f32; 3]>,
    plane_normal: Vec3A,
    constrained_edges: &HashSet<Edge>,
) -> (Vec<VertexId>, Vec<Vec<TriangleData>>) {
    // TODO Clean: See what we need for input data format of `triangulate`
    let mut vertices_data = Vec::with_capacity(vertices.len());
    for v in vertices {
        vertices_data.push(Vec3A::from_array(*v));
    }

    let mut planar_vertices =
        triangulation::transform_to_2d_planar_coordinate_system(&mut vertices_data, plane_normal);

    constrained_triangulation_from_2d_vertices(&mut planar_vertices, constrained_edges)
}

/// - `constrained_edges` **MUST**:
///     - be oriented: vi -> vj
///     - not contain edges with v0==v1
pub fn constrained_triangulation_from_2d_vertices(
    vertices: &mut Vec<Vec2>,
    constrained_edges: &HashSet<Edge>,
) -> (Vec<VertexId>, Vec<Vec<TriangleData>>) {
    let (mut triangles, container_triangle, mut debug_data) =
        wrap_and_triangulate_2d_vertices(vertices);

    // TODO Debug: consider adding debug support to this function
    apply_constraints(vertices, &mut triangles, constrained_edges);

    debug_data.push(triangles.clone());

    let indices = unholy_triangles_to_be_purify(
        &triangles,
        &container_triangle,
        constrained_edges,
        &mut debug_data,
    );

    (indices, debug_data)
}

/// Mark triangles that should be removed due to the input domains constraints
/// This supposes that all groups of constrained edges are cycles, and oriented, as shown below:
///
/// ---------------<-------------------------------------
/// |                        keep                       |
/// |   ----------->---------------------------------   |
/// |   |                   remove                  |   |
/// |   |   -------<--------     -------<--------   |   |
/// |   |   |              |     |              |   |   |
/// v   ∧   v     keep     ∧     v     keep     ∧   v   ∧
/// |   |   |              |     |              |   |   |
/// |   |   ------->--------     ------->--------   |   |
/// |   |                                           |   |
/// |   -----------<---------------------------------   |
/// |                                                   |
/// --------------->-------------------------------------
///
/// - visited_triangles: List of every triangles already checked. Pre-allocated by the caller. Reinitialized at the beginning of this function.
///
/// remove_wrapping_and_unconstrained_domains
fn unholy_triangles_to_be_purify(
    triangles: &Vec<TriangleData>,
    container_triangle: &TriangleData,
    constrained_edges: &HashSet<Edge>,
    debug_data: &mut Vec<Vec<TriangleData>>,
) -> Vec<VertexId> {
    let mut visited_triangles = vec![false; triangles.len()];

    // TODO Clean: Size approx
    let mut indices = Vec::with_capacity(3 * triangles.len());
    let mut triangles_to_explore = Vec::new();
    let container_verts: HashSet<VertexId> = HashSet::from(container_triangle.verts);

    // TODO: Debug-only
    let mut filtered_debug_triangles = Vec::new();

    // Loop over all triangles
    for (index, triangle) in triangles.iter().enumerate() {
        if visited_triangles[index] {
            continue;
        }

        let mut triangle_edge_is_constrained = [false; 3];
        for (edge_index, edge) in triangle.edges().iter().enumerate() {
            triangle_edge_is_constrained[edge_index] = constrained_edges.contains(edge);
        }

        // Stop at the first non-visited triangle with a constrained edge: we found a new unvisited domain
        if triangle_edge_is_constrained.contains(&true) {
            triangles_to_explore.clear();
            register_triangle(
                &mut indices,
                triangle,
                &container_verts,
                &mut filtered_debug_triangles,
            );
            for (edge_index, is_constrained) in triangle_edge_is_constrained.iter().enumerate() {
                if !is_constrained {
                    triangles_to_explore.push(triangle.neighbors[edge_index]);
                }
            }
            visited_triangles[index] = true;

            // Triangle-walk in the domain to register the domain's triangles
            while let Some(neighbor) = triangles_to_explore.pop() {
                match neighbor {
                    Some(triangle_index) => {
                        if visited_triangles[triangle_index] {
                            continue;
                        } else {
                            let triangle = &triangles[triangle_index];
                            register_triangle(
                                &mut indices,
                                triangle,
                                &container_verts,
                                &mut filtered_debug_triangles,
                            );
                            visited_triangles[triangle_index] = true;

                            for (edge_index, edge) in triangle.edges().iter().enumerate() {
                                if !constrained_edges.contains(edge) {
                                    // TODO Clean: Would it be quicker to also check for visited before adding to the stack
                                    triangles_to_explore.push(triangle.neighbors[edge_index]);
                                }
                            }
                        }
                    }
                    None => {
                        continue;
                    }
                }
            }
        }
    }

    debug_data.push(filtered_debug_triangles);

    indices
}

fn register_triangle(
    indices: &mut Vec<VertexId>,
    triangle: &TriangleData,
    container_verts: &HashSet<VertexId>,
    filtered_debug_triangles: &mut Vec<TriangleData>,
) {
    let mut filtered = false;
    for vert in triangle.verts.iter() {
        if container_verts.contains(vert) {
            filtered = true;
        }
    }
    if !filtered {
        indices.push(triangle.v1());
        indices.push(triangle.v2());
        indices.push(triangle.v3());
        filtered_debug_triangles.push(triangle.clone());
    }
}

fn apply_constraints(
    vertices: &Vec<Vec2>,
    triangles: &mut Vec<TriangleData>,
    constrained_edges: &HashSet<Edge>,
) {
    // Map each verticex to one of the triangles (the last) that contains it
    let mut vertex_to_triangle = vec![0; vertices.len()];
    // TODO Performance:
    // let mut vertex_to_triangle = Vec::with_capacity(vertices.len());
    // unsafe {
    //     vertex_to_triangle.set_len(vertices.len());
    // }
    for (index, triangle) in triangles.iter().enumerate() {
        vertex_to_triangle[triangle.v1()] = index;
        vertex_to_triangle[triangle.v2()] = index;
        vertex_to_triangle[triangle.v3()] = index;
    }

    for constrained_edge in constrained_edges {
        // TODO Performance: Use vertex_to_triangle to quicken this search
        // -> this will probably be done by handling the SharedEdge case when testing intersection
        // from the start of a crossed edge in register_intersected_edges
        if is_constrained_edge_inside_triangulation(triangles, constrained_edge) {
            // Nothing to do for this edge
            continue;
        };
        let constrained_edge_vertices = &constrained_edge.to_vertices(vertices);

        // Stores all of the edges that cross the constrained edge
        let intersections = register_intersected_edges(
            triangles,
            vertices,
            constrained_edge,
            constrained_edge_vertices,
            &vertex_to_triangle,
        );

        // Remove intersecting edges
        let mut new_diagonals_created = remove_crossed_edges(
            triangles,
            vertices,
            &mut vertex_to_triangle,
            constrained_edge_vertices,
            intersections,
        );

        // Restore Delaunay triangulation
        restore_delaunay_triangulation_constrained(
            triangles,
            vertices,
            &mut vertex_to_triangle,
            constrained_edge,
            &mut new_diagonals_created,
        );
    }
}

fn is_constrained_edge_inside_triangulation(
    triangles: &Vec<TriangleData>,
    constrained_edge: &Edge,
) -> bool {
    for triangle in triangles {
        if triangle.verts.contains(&constrained_edge.from)
            && triangle.verts.contains(&constrained_edge.to)
        {
            return true;
        }
    }
    false
}

fn loop_around_vertex_and_search_intersection(
    triangles: &Vec<TriangleData>,
    vertices: &Vec<Vec2>,
    from_triangle_id: TriangleId,
    constrained_edge_vertex: VertexId,
    constrained_edge_verts: &EdgeVertices,
    next_edge_index: fn(VertexId) -> TriangleVertexIndex,
) -> Option<EdgeData> {
    let mut triangle_id = from_triangle_id;
    let mut triangle = &triangles[from_triangle_id];

    // We search by circling around the first vertex of the constrained edge
    // We use `triangles.len()` as an upper bound on the number of triangles around a vertex
    for _ in 0..triangles.len() {
        let Some(vert_index) = triangle.vertex_index(constrained_edge_vertex) else {
            // Not possible
            error!("Internal error, search_first_interstected_quad found a triangle without the constrained edge starting vertex, triangle {} {:?}, starting_vertex {}", triangle_id, triangle, constrained_edge_vertex);
            return None;
        };

        let edge_index = opposite_edge_index(vert_index);
        let edge = triangle.edge(edge_index);
        let edge_vertices = edge.to_vertices(vertices);

        if egdes_intersect(&constrained_edge_verts, &edge_vertices)
            == EdgesIntersectionResult::Crossing
        {
            if let Some(neighbor_triangle_id) = triangle.neighbors[edge_index] {
                return Some(EdgeData {
                    from_triangle_id: triangle_id,
                    to_triangle_id: neighbor_triangle_id,
                    edge,
                });
            } else {
                // Not possible
                error!("Internal error, search_first_interstected_quad found a triangle with a crossed edge but no neihgbor, triangle {} {:?}, starting_vertex {}", triangle_id, triangle, constrained_edge_vertex);
                break;
            }
        }

        match triangle.neighbors[next_edge_index(vert_index)] {
            Some(neighbor_id) => {
                triangle_id = neighbor_id;
                triangle = &triangles[triangle_id];
            }
            None => break,
        }
    }
    None
}

/// - visited_triangles: List of every triangles already checked. Pre-allocated and initialized by the caller
fn search_first_interstected_quad(
    triangles: &Vec<TriangleData>,
    vertices: &Vec<Vec2>,
    constrained_edge_first_vertex: VertexId,
    constrained_edge: &EdgeVertices,
    start_triangle: TriangleId,
) -> EdgeData {
    // Search clockwise
    if let Some(intersection) = loop_around_vertex_and_search_intersection(
        triangles,
        vertices,
        start_triangle,
        constrained_edge_first_vertex,
        constrained_edge,
        next_clockwise_edge_index_around,
    ) {
        return intersection;
    }

    // Get the counter-clockwise neighbor of the starting triangle
    let triangle = &triangles[start_triangle];
    // Not having the constrained edge first vertex in the starting triangle is not possible
    let vert_index = triangle
        .vertex_index(constrained_edge_first_vertex)
        .unwrap();
    let edge_index = next_counter_clockwise_edge_index_around(vert_index);
    let neighbor_triangle = triangle.neighbors[edge_index].expect(&format!("Internal error, search_first_interstected_quad found no triangle with the constrained edge intersection, starting_vertex {}", constrained_edge_first_vertex));

    // Search counterclockwise
    loop_around_vertex_and_search_intersection(
        triangles,
        vertices,
        neighbor_triangle,
        constrained_edge_first_vertex,
        constrained_edge,
        next_counter_clockwise_edge_index_around,
    ).expect(&format!("Internal error, search_first_interstected_quad found no triangle with the constrained edge intersection, starting_vertex {}", constrained_edge_first_vertex))
}

/// - visited_triangles: List of every triangles already checked. Pre-allocated by the caller. Reinitialized at the beginning of this function.
fn register_intersected_edges(
    triangles: &Vec<TriangleData>,
    vertices: &Vec<Vec2>,
    constrained_edge: &Edge,
    constrained_edge_vertices: &EdgeVertices,
    vertex_to_triangle: &Vec<TriangleId>,
) -> VecDeque<EdgeData> {
    let mut intersection = search_first_interstected_quad(
        triangles,
        vertices,
        constrained_edge.from,
        constrained_edge_vertices,
        vertex_to_triangle[constrained_edge.from],
    );

    let mut intersections = VecDeque::new();
    intersections.push_back(intersection.clone());

    // Search all the edges intersected by this constrained edge
    // We use `triangles.len()` as an upper bound on the number of triangles intersecting a constrained edge
    for _ in 0..triangles.len() {
        let triangle_id = intersection.to_triangle_id;
        let triangle = &triangles[triangle_id];

        // Stop if we reached the end of the constrained edge
        if triangle.verts.contains(&constrained_edge.to) {
            break;
        }

        match get_next_triangle_edge_intersection(
            &triangle.to_vertices(vertices),
            constrained_edge_vertices,
            triangle,
            triangle_id,
            &intersection.edge,
        ) {
            Some(next_intersection) => {
                intersection = next_intersection;
                intersections.push_back(intersection.clone());
            }
            None => {
                // TODO Doc: Justify why it is not possible
                error!("Internal error, triangle {} {:?} should be intersecting the constrained edge {:?}", triangle_id, triangle, constrained_edge);
                break;
            }
        }
    }

    intersections
}

#[derive(Clone, Debug)]
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

///               q4
///              /  \
///            /  Tt  \
///          /          \
///      q2=c2 -------- q1=c1
///          \          /
///            \  Tf  /
///              \  /
///               q3
/// where:
/// Tf: from triangle
/// Tt: to triangle
/// c1, c2: crossed edge
/// q1, q2, q3, q4: quad coords
impl EdgeData {
    pub fn to_quad(&self, triangles: &Vec<TriangleData>) -> Quad {
        Quad::new([
            self.edge.from,
            self.edge.to,
            triangles[self.from_triangle_id].get_opposite_vertex_index(&self.edge),
            triangles[self.to_triangle_id].get_opposite_vertex_index(&self.edge),
        ])
    }
}

/// - visited_triangles: List of every triangles already checked. Pre-allocated and initialized by the caller
fn get_next_triangle_edge_intersection(
    triangle_verts: &TriangleVertices,
    constrained_edge: &EdgeVertices,
    triangle: &TriangleData,
    triangle_id: TriangleId,
    from_crossed_edge: &Edge,
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
        if egdes_intersect(constrained_edge, edge_vertices) == EdgesIntersectionResult::Crossing {
            if let Some(neighbor_triangle_id) = triangle.neighbors[edge_index] {
                return Some(EdgeData {
                    // crossed_edge_index: edge_index,
                    from_triangle_id: triangle_id,
                    to_triangle_id: neighbor_triangle_id,
                    edge,
                });
            }
        }
    }
    None
}

///               q4
///              /  \
///      Ttl   /  Tt  \    Ttr
///          /          \
///      q2=c2 -------- q1=c1
///          \          /
///            \  Tf  /
///      Tfl     \  /      Tfr
///               q3
///
///             Becomes
///
///               q4
///              / | \
///      Ttl   /   |   \   Ttr
///          /     |     \
///      q2=c2  Tf | Tt  q1=c1
///          \     |     /
///            \   |   /
///      Tfl     \ | /     Tfr
///               q3
/// where:
/// Tf: triangle from
/// Tt: triangle to
/// c1, c2: crossed edge
/// q1, q2, q3, q4: quad coords
pub(crate) fn swap_quad_diagonal(
    triangles: &mut Vec<TriangleData>,
    vertex_to_triangle: &mut Vec<TriangleId>,
    edges_data_collections: &mut [&mut VecDeque<EdgeData>],
    edge_data: &EdgeData,
    quad: &Quad,
) {
    let from = edge_data.from();
    let to = edge_data.to();

    // TODO Design: Memorize from_crossed_edge_index in Intersection ?
    // `tt`: triangle_to
    // `tf`: triangle_from
    let tt_left_edge_index = triangles[to].opposite_edge_index_from_vertex(quad.v1());
    let tt_left_neighbor = triangles[to].neighbors[tt_left_edge_index];
    let tt_right_neighbor = triangles[to].neighbors[next_clockwise_edge_index(tt_left_edge_index)];

    let tf_left_edge_index = triangles[from].opposite_edge_index_from_vertex(quad.v1());
    let tf_left_neighbor = triangles[from].neighbors[tf_left_edge_index];
    let tf_right_neighbor =
        triangles[from].neighbors[next_counter_clockwise_edge_index(tf_left_edge_index)];

    triangles[from].verts = [quad.v4(), quad.v2(), quad.v3()];
    triangles[to].verts = [quad.v4(), quad.v3(), quad.v1()];

    triangles[from].neighbors = [tt_left_neighbor, tf_left_neighbor, Some(to)];
    triangles[to].neighbors = [Some(from), tf_right_neighbor, tt_right_neighbor];

    triangulation::update_triangle_neighbour(tt_left_neighbor, Some(to), Some(from), triangles);
    triangulation::update_triangle_neighbour(tf_right_neighbor, Some(from), Some(to), triangles);

    vertex_to_triangle[quad.v1()] = to;
    vertex_to_triangle[quad.v2()] = from;

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
        if let Some(tt_left_neighbor) = tt_left_neighbor {
            if edge_data.from_triangle_id == tt_left_neighbor && edge_data.to_triangle_id == to {
                edge_data.to_triangle_id = from;
            } else if edge_data.from_triangle_id == to
                && edge_data.to_triangle_id == tt_left_neighbor
            {
                edge_data.from_triangle_id = from;
                // TODO Design: inter.from_edge_index becomes EDGE_12
            }
        }
        if let Some(tf_right_neighbor) = tf_right_neighbor {
            if edge_data.from_triangle_id == tf_right_neighbor && edge_data.to_triangle_id == from {
                edge_data.to_triangle_id = to;
            } else if edge_data.from_triangle_id == from
                && edge_data.to_triangle_id == tf_right_neighbor
            {
                edge_data.from_triangle_id = to;
                // TODO Design: inter.from_edge_index becomes EDGE_23
            }
        }
        // TODO Design: IF we decide to store from_edge_index, we also need to update tfl<->tf and ttr<->tt pairs
    }
}

fn remove_crossed_edges(
    triangles: &mut Vec<TriangleData>,
    vertices: &Vec<Vec2>,
    vertex_to_triangle: &mut Vec<TriangleId>,
    constrained_edge_vertices: &EdgeVertices,
    mut intersections: VecDeque<EdgeData>,
) -> VecDeque<EdgeData> {
    let mut new_diagonals_created = VecDeque::new();

    while let Some(intersection) = intersections.pop_front() {
        let quad = intersection.to_quad(triangles);
        let quad_vertices = quad.to_vertices(vertices);
        // If the quad is not convex (diagonals do not cross) or if an edge tip lie on the other edge,
        // we skip this edge and put it on the stack to re-process it later
        if quad_vertices.diagonals_intersection_test() != EdgesIntersectionResult::Crossing {
            intersections.push_back(intersection);
        }
        // Swap the diagonal of this strictly convex quadrilateral if the two diagonals cross normaly
        else {
            swap_quad_diagonal(
                triangles,
                vertex_to_triangle,
                &mut [&mut intersections, &mut new_diagonals_created],
                &intersection,
                &quad,
            );
            // TODO Clean: this could be returned by swap_quad_diagonal ?
            let edge_data = EdgeData::new(
                intersection.from_triangle_id,
                intersection.to_triangle_id,
                Edge::new(quad.v3(), quad.v4()),
            );
            // If the new diagonal still intersects the constrained edge,
            // then place it on the list of intersecting edges
            if egdes_intersect(
                &(quad_vertices.q3(), quad_vertices.q4()),
                constrained_edge_vertices,
            ) == EdgesIntersectionResult::Crossing
            {
                intersections.push_back(edge_data);
            } else {
                new_diagonals_created.push_back(edge_data);
            }
        }
    }
    new_diagonals_created
}

fn restore_delaunay_triangulation_constrained(
    triangles: &mut Vec<TriangleData>,
    vertices: &Vec<Vec2>,
    vertex_to_triangle: &mut Vec<TriangleId>,
    constrained_edge: &Edge,
    new_diagonals_created: &mut VecDeque<EdgeData>,
) {
    while let Some(new_edge) = new_diagonals_created.pop_front() {
        if new_edge.edge.undirected_equals(constrained_edge) {
            continue;
        } else {
            let quad = new_edge.to_quad(triangles);
            let quad_vertices = quad.to_vertices(vertices);
            if is_vertex_in_triangle_circumcircle(&quad_vertices.0[0..=3], quad_vertices.q4()) {
                swap_quad_diagonal(
                    triangles,
                    vertex_to_triangle,
                    &mut [new_diagonals_created],
                    &new_edge,
                    &quad,
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

#[cfg(test)]
mod tests {
    use std::collections::VecDeque;

    use bevy::math::Vec2;

    use crate::triangulation::{Edge, Quad, TriangleData};

    use super::{swap_quad_diagonal, EdgeData};

    #[test]
    fn swap_quad_diag_constrained() {
        let mut vertices = Vec::<Vec2>::new();
        vertices.push(Vec2::new(0.5, 3.));
        vertices.push(Vec2::new(-2., -2.));
        vertices.push(Vec2::new(1., -4.));
        vertices.push(Vec2::new(3., -2.));

        let triangle_1 = TriangleData {
            verts: [0, 1, 3],
            neighbors: [None, Some(1), None],
        };

        let triangle_2 = TriangleData {
            verts: [1, 2, 3],
            neighbors: [None, None, Some(0)],
        };

        let mut triangles = vec![triangle_1, triangle_2];

        let mut vertex_to_triangle = vec![0; vertices.len()];
        for (index, triangle) in triangles.iter().enumerate() {
            vertex_to_triangle[triangle.v1()] = index;
            vertex_to_triangle[triangle.v2()] = index;
            vertex_to_triangle[triangle.v3()] = index;
        }

        let mut new_diagonals_created = VecDeque::new();

        let intersection = EdgeData::new(0, 1, Edge::new(1, 3));

        let mut edge_data_collection = [&mut new_diagonals_created];

        let quad = Quad::new([1, 3, 0, 2]);

        swap_quad_diagonal(
            &mut triangles,
            &mut vertex_to_triangle,
            &mut edge_data_collection,
            &intersection,
            &quad,
        );

        assert_eq!(2, triangles.len());
        assert_eq!([2, 3, 0], triangles[0].verts);
        assert_eq!([2, 0, 1], triangles[1].verts);
    }
}
