use std::collections::VecDeque;

use bevy::{
    log::error,
    math::{Vec2, Vec3A},
    utils::hashbrown::HashSet,
};

use crate::utils::{egdes_intersect, EdgesIntersectionResult};

use super::{
    triangulation::{self, wrap_and_triangulate_2d_vertices},
    Edge, EdgeVertices, Quad, TriangleData, TriangleId, TriangleVertexIndex, TriangleVertices,
    VertexId, EDGE_23, EDGE_31,
};

/// plane_normal must be normalized
/// vertices must all belong to a 3d plane
/// constrained edges must be oriented: vi -> vj
pub fn constrained_triangulation_from_3d_planar_vertices(
    vertices: &Vec<[f32; 3]>,
    plane_normal: Vec3A,
    constrained_edges: &HashSet<Edge>,
) -> (Vec<VertexId>, Vec<Vec<TriangleData>>) {
    // TODO See what we need for input data format of `triangulate`
    let mut vertices_data = Vec::with_capacity(vertices.len());
    for v in vertices {
        vertices_data.push(Vec3A::from_array(*v));
    }

    let mut planar_vertices =
        triangulation::transform_to_2d_planar_coordinate_system(&mut vertices_data, plane_normal);

    constrained_triangulation_from_2d_vertices(&mut planar_vertices, constrained_edges)
}

/// constrained_edges should not contain edges with v0==v1
pub fn constrained_triangulation_from_2d_vertices(
    vertices: &mut Vec<Vec2>,
    constrained_edges: &HashSet<Edge>,
) -> (Vec<VertexId>, Vec<Vec<TriangleData>>) {
    let (mut triangles, container_triangle, mut debug_data) =
        wrap_and_triangulate_2d_vertices(vertices);

    // TODO Debug: consider adding debug support to this function

    apply_constraints(vertices, &mut triangles, constrained_edges);

    let indices = remove_wrapping_and_unconstrained_domains(
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
fn remove_wrapping_and_unconstrained_domains(
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
    // TODO Check performances
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
        // TODO Use vertex_to_triangle to quicken this search
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
        let mut new_diagonals_created = remove_crossing_edges(
            triangles,
            vertices,
            &mut vertex_to_triangle,
            constrained_edge,
            constrained_edge_vertices,
            intersections,
        );

        // Restore Delaunay triangulation
        restore_delaunay_triangulation_constrained(
            triangles,
            vertices,
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
) -> Option<TriangleEdgeIntersection> {
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

        let edge_index = TriangleData::get_opposite_edge_index(vert_index);
        let edge = triangle.edge(edge_index);
        let edge_vertices = edge.to_vertices(vertices);

        if egdes_intersect(&constrained_edge_verts, &edge_vertices)
            == EdgesIntersectionResult::Crossing
        {
            if let Some(neighbor_triangle_id) = triangle.neighbors[edge_index] {
                return Some(TriangleEdgeIntersection {
                    from_triangle_id: triangle_id,
                    to_triangle_id: neighbor_triangle_id,
                    crossed_edge: edge,
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
) -> TriangleEdgeIntersection {
    // Search left
    if let Some(intersection) = loop_around_vertex_and_search_intersection(
        triangles,
        vertices,
        start_triangle,
        constrained_edge_first_vertex,
        constrained_edge,
        TriangleData::get_left_edge_index_around,
    ) {
        return intersection;
    }

    // Get the right neighbor of the starting triangle
    let triangle = &triangles[start_triangle];
    // Not having the constrained edge first vertex in the starting triangle is not possible
    let vert_index = triangle
        .vertex_index(constrained_edge_first_vertex)
        .unwrap();
    let edge_index = TriangleData::get_right_edge_index_around(vert_index);
    let neighbor_triangle = triangle.neighbors[edge_index].expect(&format!("Internal error, search_first_interstected_quad found no triangle with the constrained edge intersection, starting_vertex {}", constrained_edge_first_vertex));

    // Search right
    loop_around_vertex_and_search_intersection(
        triangles,
        vertices,
        neighbor_triangle,
        constrained_edge_first_vertex,
        constrained_edge,
        TriangleData::get_right_edge_index_around,
    ).expect(&format!("Internal error, search_first_interstected_quad found no triangle with the constrained edge intersection, starting_vertex {}", constrained_edge_first_vertex))
}

/// - visited_triangles: List of every triangles already checked. Pre-allocated by the caller. Reinitialized at the beginning of this function.
fn register_intersected_edges(
    triangles: &Vec<TriangleData>,
    vertices: &Vec<Vec2>,
    constrained_edge: &Edge,
    constrained_edge_vertices: &EdgeVertices,
    vertex_to_triangle: &Vec<TriangleId>,
) -> VecDeque<TriangleEdgeIntersection> {
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
            &intersection.crossed_edge,
        ) {
            Some(next_intersection) => {
                intersection = next_intersection;
                intersections.push_back(intersection.clone());
            }
            None => {
                // TODO Justify why it is not possible
                error!("Internal error, triangle {} {:?} should be intersecting the constrained edge {:?}", triangle_id, triangle, constrained_edge);
                break;
            }
        }
    }

    intersections
}

#[derive(Clone)]
struct TriangleEdgeIntersection {
    from_triangle_id: TriangleId,
    to_triangle_id: TriangleId,
    crossed_edge: Edge,
}
impl TriangleEdgeIntersection {
    fn new(from_triangle_id: TriangleId, to_triangle_id: TriangleId, crossed_edge: Edge) -> Self {
        Self {
            from_triangle_id,
            to_triangle_id,
            crossed_edge,
        }
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
impl TriangleEdgeIntersection {
    pub fn to_quad(&self, triangles: &Vec<TriangleData>) -> Quad {
        Quad::new([
            self.crossed_edge.from,
            self.crossed_edge.to,
            triangles[self.from_triangle_id].get_opposite_vertex_index(&self.crossed_edge),
            triangles[self.to_triangle_id].get_opposite_vertex_index(&self.crossed_edge),
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
) -> Option<TriangleEdgeIntersection> {
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
        // TODO: Add a note about why we only care about the `Crossing` result
        if egdes_intersect(constrained_edge, edge_vertices) == EdgesIntersectionResult::Crossing {
            if let Some(neighbor_triangle_id) = triangle.neighbors[edge_index] {
                return Some(TriangleEdgeIntersection {
                    // crossed_edge_index: edge_index,
                    from_triangle_id: triangle_id,
                    to_triangle_id: neighbor_triangle_id,
                    crossed_edge: edge,
                });
            }
        }
    }
    None
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
///
///             Becomes
///
///               q4
///              / | \
///            /   |   \
///          /     |     \
///      q2=c2  Tf | Tt  q1=c1
///          \     |     /
///            \   |   /
///              \ | /
///               q3
/// where:
/// Tf: from triangle
/// Tt: to triangle
/// c1, c2: crossed edge
/// q1, q2, q3, q4: quad coords
fn new_swap_quad_diagonal(
    triangles: &mut Vec<TriangleData>,
    vertex_to_triangle: &mut Vec<TriangleId>,
    intersection: &TriangleEdgeIntersection,
    quad: &Quad,
) {
    let from = intersection.from_triangle_id;
    let to = intersection.to_triangle_id;

    // TODO Memorize from_crossed_edge_index in Intersection
    // TODO Compute to_crossed_edge_index (here?)
    // TODO TtLeft : get it form somewhere or just compute it from triangleData
    // TODO Simplify
    let tt_left_edge_index = triangles[to].get_opposite_edge_index_from_vertex(quad.v1());
    let tt_left_neighbor = triangles[to].neighbors[tt_left_edge_index];
    let tt_right_neighbor = triangles[to].neighbors[(tt_left_edge_index + 2) % 3];

    let tf_left_edge_index = triangles[from].get_opposite_edge_index_from_vertex(quad.v1());
    let tf_left_neighbor = triangles[from].neighbors[tf_left_edge_index];
    let tf_right_neighbor = triangles[from].neighbors[(tf_left_edge_index + 2) % 3];

    triangles[from].verts = [quad.v4(), quad.v2(), quad.v3()];
    triangles[to].verts = [quad.v4(), quad.v3(), quad.v1()];

    triangles[from].neighbors = [
        tt_left_neighbor,
        tf_left_neighbor,
        Some(intersection.to_triangle_id),
    ];
    triangles[to].neighbors = [
        Some(intersection.from_triangle_id),
        tf_right_neighbor,
        tt_right_neighbor,
    ];

    // TODO Should exist on a non optional triangle
    triangulation::update_triangle_neighbour(
        tt_left_neighbor,
        Some(intersection.to_triangle_id),
        Some(intersection.from_triangle_id),
        triangles,
    );
    triangulation::update_triangle_neighbour(
        tf_right_neighbor,
        Some(intersection.from_triangle_id),
        Some(intersection.to_triangle_id),
        triangles,
    );

    vertex_to_triangle[quad.v1()] = intersection.to_triangle_id;
    vertex_to_triangle[quad.v2()] = intersection.from_triangle_id;
}

fn remove_crossing_edges(
    triangles: &mut Vec<TriangleData>,
    vertices: &Vec<Vec2>,
    vertex_to_triangle: &mut Vec<TriangleId>,
    constrained_edge: &Edge,
    constrained_edge_vertices: &EdgeVertices,
    mut intersections: VecDeque<TriangleEdgeIntersection>,
) -> VecDeque<QuadDiag> {
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
            // let (t3, t4) = legacy_swap_quad_diagonal(
            //     triangles,
            //     quad.v4(),
            //     intersection.from_triangle_id,
            //     intersection.to_triangle_id,
            // );
            // let (t3, t4) =
            new_swap_quad_diagonal(triangles, vertex_to_triangle, &intersection, &quad);

            // TODO Update data in intersections
            // TODO Update data in new_diagonals_created

            // If the new diagonal still intersects the constrained edge,
            // then place it on the list of intersecting edges
            if egdes_intersect(
                &(quad_vertices.q3(), quad_vertices.q4()),
                constrained_edge_vertices,
            ) == EdgesIntersectionResult::Crossing
            {
                intersections.push_back(TriangleEdgeIntersection::new(
                    intersection.from_triangle_id,
                    intersection.to_triangle_id,
                    Edge::new(quad.v3(), quad.v4()),
                ));
            } else {
                new_diagonals_created.push_back(QuadDiag::new(
                    intersection.from_triangle_id,
                    intersection.to_triangle_id,
                    Edge::new(quad.v3(), quad.v4()),
                ));
            }
        }
    }
    new_diagonals_created
}

struct QuadDiag {
    from_triangle: TriangleId,
    to_triangle: TriangleId,
    edge: Edge,
}
impl QuadDiag {
    #[inline]
    fn new(from: TriangleId, to: TriangleId, edge: Edge) -> Self {
        Self {
            from_triangle: from,
            to_triangle: to,
            edge,
        }
    }
}

fn restore_delaunay_triangulation_constrained(
    triangles: &mut Vec<TriangleData>,
    vertices: &Vec<Vec2>,
    constrained_edge: &Edge,
    new_diagonals_created: &mut VecDeque<QuadDiag>,
) {
    while let Some(new_diag) = new_diagonals_created.pop_front() {
        if new_diag.edge.undirected_equals(constrained_edge) {
            continue;
        } else {
            // let (quad_diag_swapped, quad_diag) = swap_quad_diagonal_to_edge(
            //     triangles,
            //     vertices,
            //     new_diag.edge.from,
            //     new_diag.from_triangle,
            //     new_diag.to_triangle,
            // );
            // if quad_diag_swapped {
            //     new_diagonals_created.push_back(quad_diag);
            // }
        }
    }
}

// fn swap_quad_diagonal_to_edge(
//     triangles: &mut Vec<TriangleData>,
//     vertices: &Vec<Vec2>,
//     vertex_id: VertexId,
//     triangle_id: TriangleId,
//     adjacent_triangle_id: TriangleId,
// ) -> (bool, QuadDiag) {
//     let adjacent_triangle = &triangles[adjacent_triangle_id];
//     let (quad, triangle_3_id, triangle_4_id) =
//         if adjacent_triangle.neighbor12() == Some(triangle_id) {
//             (
//                 Quad {
//                     v1: adjacent_triangle.v2(),
//                     v2: adjacent_triangle.v1(),
//                     v3: adjacent_triangle.v3(),
//                     v4: vertex_id,
//                 },
//                 adjacent_triangle.neighbor23(),
//                 adjacent_triangle.neighbor31(),
//             )
//         } else if adjacent_triangle.neighbor23() == Some(triangle_id) {
//             (
//                 Quad {
//                     v1: adjacent_triangle.v3(),
//                     v2: adjacent_triangle.v2(),
//                     v3: adjacent_triangle.v1(),
//                     v4: vertex_id,
//                 },
//                 adjacent_triangle.neighbor31(),
//                 adjacent_triangle.neighbor12(),
//             )
//         } else {
//             (
//                 Quad {
//                     v1: adjacent_triangle.v1(),
//                     v2: adjacent_triangle.v3(),
//                     v3: adjacent_triangle.v2(),
//                     v4: vertex_id,
//                 },
//                 adjacent_triangle.neighbor12(),
//                 adjacent_triangle.neighbor23(),
//             )
//         };

//     // Check if the vertex is on the circumcircle of the adjacent triangle:
//     let mut swapped_quad_diagonal = false;
//     if triangulation::is_vertex_in_triangle_circumcircle(
//         vertices[quad.v1],
//         vertices[quad.v2],
//         vertices[quad.v3],
//         vertices[vertex_id],
//     ) {
//         // The triangle containing P as a vertex and the unstacked triangle form a convex quadrilateral whose diagonal is drawn in the wrong direction.
//         // Swap this diagonal so that two old triangles are replaced by two new triangles and the structure of the Delaunay triangulation is locally restored.
//         triangulation::update_triangle_neighbour(
//             triangle_3_id,
//             Some(adjacent_triangle_id),
//             Some(triangle_id),
//             triangles,
//         );
//         triangulation::update_triangle_neighbour(
//             triangles[triangle_id].neighbor31(),
//             Some(triangle_id),
//             Some(adjacent_triangle_id),
//             triangles,
//         );

//         triangles[triangle_id].verts = [quad.v4, quad.v1, quad.v3];

//         triangles[adjacent_triangle_id].verts = [quad.v4, quad.v3, quad.v2];

//         triangles[adjacent_triangle_id].neighbors = [
//             Some(triangle_id),
//             triangle_4_id,
//             triangles[triangle_id].neighbor31(),
//         ];

//         triangles[triangle_id].neighbors[EDGE_23] = triangle_3_id;
//         triangles[triangle_id].neighbors[EDGE_31] = Some(adjacent_triangle_id);

//         swapped_quad_diagonal = true;
//     }
//     (
//         swapped_quad_diagonal,
//         // TODO Justify unwraps
//         QuadDiag::new(
//             triangle_3_id.unwrap(),
//             triangle_4_id.unwrap(),
//             Edge::new(quad.v1, quad.v2),
//         ),
//     )
// }
