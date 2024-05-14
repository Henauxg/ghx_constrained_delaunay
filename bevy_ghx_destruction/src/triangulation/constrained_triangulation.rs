use std::collections::VecDeque;

use bevy::{
    math::{Vec2, Vec3A},
    utils::hashbrown::HashSet,
};

use crate::utils::{egdes_intersect, EdgesIntersectionResult};

use super::{
    triangulation::{self, wrap_and_triangulate_2d_vertices},
    Edge, Neighbor, Quad, TriangleData, TriangleId, VertexId, EDGE_23, EDGE_31,
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
fn remove_wrapping_and_unconstrained_domains(
    triangles: &Vec<TriangleData>,
    container_triangle: &TriangleData,
    constrained_edges: &HashSet<Edge>,
    debug_data: &mut Vec<Vec<TriangleData>>,
) -> Vec<VertexId> {
    // TODO Clean: Size approx
    let mut indices = Vec::with_capacity(3 * triangles.len());
    let mut visited_triangles = vec![false; triangles.len()];
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
    vertices: &mut Vec<Vec2>,
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

        // Stores all of the edges that cross the constrained edge
        let intersected_edges =
            register_intersected_edges(triangles, vertices, constrained_edge, &vertex_to_triangle);

        // Remove intersecting edges
        let mut new_diagonals_created =
            remove_crossing_edges(triangles, vertices, constrained_edge, intersected_edges);

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

fn register_intersected_edges(
    triangles: &Vec<TriangleData>,
    vertices: &mut Vec<Vec2>,
    constrained_edge: &Edge,
    vertex_to_triangle: &Vec<TriangleId>,
) -> VecDeque<(TriangleId, Edge)> {
    // List of every triangles already checked
    //TODO: change data structure
    let mut visited_triangles = vec![false; triangles.len()];

    let mut intersected_edges = VecDeque::new();

    // providing a starting triangle to begin the search for the edges which cross the constrained edge
    let mut current_triangle_id = vertex_to_triangle[constrained_edge.from];
    visited_triangles[current_triangle_id] = true;

    let mut crossing = false;

    // we circle constrained_edge.0
    while !crossing {
        // Vertices of the current triangle
        let triangle_vertex_1 = vertices[triangles[current_triangle_id].v1()];
        let triangle_vertex_2 = vertices[triangles[current_triangle_id].v2()];
        let triangle_vertex_3 = vertices[triangles[current_triangle_id].v3()];

        (crossing, _, _) = constrained_edge_intersect_triangle(
            triangles,
            vertices,
            constrained_edge,
            current_triangle_id,
            &visited_triangles,
        );

        if crossing {
            break;
        }

        // If the current triangle has already been visited, we go the right
        if visited_triangles[current_triangle_id] {
            if triangle_vertex_1 == vertices[constrained_edge.from] {
                match triangles[current_triangle_id].neighbor12() {
                    Some(neighbour_triangle_id) => current_triangle_id = neighbour_triangle_id,
                    None => todo!(), // TODO: show error,
                }
            } else if triangle_vertex_2 == vertices[constrained_edge.from] {
                match triangles[current_triangle_id].neighbor23() {
                    Some(neighbour_triangle_id) => current_triangle_id = neighbour_triangle_id,
                    None => todo!(), // TODO: show error,
                }
            } else if triangle_vertex_3 == vertices[constrained_edge.from] {
                match triangles[current_triangle_id].neighbor31() {
                    Some(neighbour_triangle_id) => current_triangle_id = neighbour_triangle_id,
                    None => todo!(), // TODO: show error,
                }
            }
        } else {
            //march from one triangle to the next if there are no crossing edges
            if triangle_vertex_1 == vertices[constrained_edge.from] {
                match triangles[current_triangle_id].neighbor31() {
                    // By default, we turn to the left
                    Some(neighbour_triangle_id) => {
                        // if the left triangle exist
                        if !visited_triangles[neighbour_triangle_id] {
                            // if the left triangle hasn't been visited yet
                            current_triangle_id = neighbour_triangle_id; // we got to this triangle
                            visited_triangles[current_triangle_id] = true;
                        } else {
                            // if the left triangle has been visited
                            current_triangle_id =
                                triangles[current_triangle_id].neighbor12().unwrap(); // we turn to the right triangle
                            visited_triangles[current_triangle_id] = true;
                        }
                    }
                    None => {
                        // If there is no triangle in the left direction, we got to the right direction
                        current_triangle_id = triangles[current_triangle_id].neighbor12().unwrap();
                        visited_triangles[current_triangle_id] = true;
                    }
                }
            } else if triangle_vertex_2 == vertices[constrained_edge.from] {
                match triangles[current_triangle_id].neighbor12() {
                    // By default, we turn to the left
                    Some(neighbour_triangle_id) => {
                        // if the left triangle exist
                        if !visited_triangles[neighbour_triangle_id] {
                            // if the left triangle hasn't been visited yet
                            current_triangle_id = neighbour_triangle_id; // we got to this triangle
                            visited_triangles[current_triangle_id] = true;
                        } else {
                            // if the left triangle has been visited
                            current_triangle_id =
                                triangles[current_triangle_id].neighbor23().unwrap(); // we turn to the right triangle
                            visited_triangles[current_triangle_id] = true;
                        }
                    }
                    None => {
                        // If there is no triangle in the left direction, we got to the right direction
                        current_triangle_id = triangles[current_triangle_id].neighbor23().unwrap();
                        visited_triangles[current_triangle_id] = true;
                    }
                }
            } else if triangle_vertex_3 == vertices[constrained_edge.from] {
                match triangles[current_triangle_id].neighbor23() {
                    // By default, we turn to the left
                    Some(neighbour_triangle_id) => {
                        // if the left triangle exist
                        if !visited_triangles[neighbour_triangle_id] {
                            // if the left triangle hasn't been visited yet
                            current_triangle_id = neighbour_triangle_id; // we got to this triangle
                            visited_triangles[current_triangle_id] = true;
                        } else {
                            // if the left triangle has been visited
                            current_triangle_id =
                                triangles[current_triangle_id].neighbor31().unwrap(); // we turn to the right triangle
                            visited_triangles[current_triangle_id] = true;
                        }
                    }
                    None => {
                        // If there is no triangle in the left direction, we got to the right direction
                        current_triangle_id = triangles[current_triangle_id].neighbor31().unwrap();
                        visited_triangles[current_triangle_id] = true;
                    }
                }
            }
        }
    }

    // Flag to check whenever we are at the end of the constrained edge
    let mut arrived_to_constrained_edge_end = false;

    while !arrived_to_constrained_edge_end {
        visited_triangles[current_triangle_id] = true;

        // Vertices of the current triangle
        let triangle_vertex_1 = vertices[triangles[current_triangle_id].v1()];
        let triangle_vertex_2 = vertices[triangles[current_triangle_id].v2()];
        let triangle_vertex_3 = vertices[triangles[current_triangle_id].v3()];

        let (_, crossing_edge, crossed_triangle) = constrained_edge_intersect_triangle(
            triangles,
            vertices,
            constrained_edge,
            current_triangle_id,
            &visited_triangles,
        );

        // Check if we are at constrained_edge.1
        if triangle_vertex_1 == vertices[constrained_edge.to]
            || triangle_vertex_2 == vertices[constrained_edge.to]
            || triangle_vertex_3 == vertices[constrained_edge.to]
        {
            arrived_to_constrained_edge_end = true;
        } else {
            current_triangle_id = crossed_triangle.unwrap();
            intersected_edges.push_back((current_triangle_id, crossing_edge.unwrap()));
        }
    }

    intersected_edges
}

fn remove_crossing_edges(
    triangles: &mut Vec<TriangleData>,
    vertices: &mut Vec<Vec2>,
    constrained_edge: &Edge,
    mut intersected_edges: VecDeque<(TriangleId, Edge)>,
) -> VecDeque<QuadDiag> {
    let mut new_diagonals_created = VecDeque::new();

    while let Some((current_triangle_id, current_crossing_edge)) = intersected_edges.pop_front() {
        // Construct a quad from the current triangle and its edge by looking at the neighbour triangle of this edge
        let (quad, triangle_id, adjacent_triangle_id) = find_quad_from_triangle_and_crossing_edge(
            current_triangle_id,
            &current_crossing_edge,
            triangles,
        );

        // If the quad is not convex or if an edge tip lie on the other edge, we skip this edge
        if quad_diagonals_intersection_test(&quad, vertices) != EdgesIntersectionResult::Crossing {
            intersected_edges.push_back((current_triangle_id, current_crossing_edge));
        }
        // Swap the diagonal of this strictly convex quadrilateral if the two diagonals cross normaly
        else {
            let (t3, t4) =
                swap_quad_diagonal(quad.v4, adjacent_triangle_id, triangle_id, triangles);

            // If the new diagonal still intersects the constrained edge,
            // then place it on the list of intersecting edges
            let quad_vertex_2 = vertices[quad.v2];
            let quad_vertex_4 = vertices[quad.v4];
            let constrained_edge_vertex_1 = vertices[constrained_edge.from];
            let constrained_edge_vertex_2 = vertices[constrained_edge.to];

            if egdes_intersect(
                quad_vertex_2,
                quad_vertex_4,
                constrained_edge_vertex_1,
                constrained_edge_vertex_2,
            ) == EdgesIntersectionResult::Crossing
            {
                intersected_edges.push_back((t3.unwrap(), Edge::new(quad.v2, quad.v4)));
            } else {
                // TODO Justfiy unwraps
                new_diagonals_created.push_back(QuadDiag::new(
                    t3.unwrap(),
                    t4.unwrap(),
                    Edge::new(quad.v2, quad.v4),
                ));
            }
        }
    }
    new_diagonals_created
}

fn swap_quad_diagonal(
    vertex_id: VertexId,
    triangle_id: TriangleId,
    adjacent_triangle_id: TriangleId,
    triangles: &mut Vec<TriangleData>,
) -> (Option<TriangleId>, Option<TriangleId>) {
    let adjacent_triangle = &triangles[adjacent_triangle_id];
    let (quad, triangle_3_id, triangle_4_id) =
        if adjacent_triangle.neighbor12() == Some(triangle_id) {
            (
                Quad {
                    v1: adjacent_triangle.v2(),
                    v2: adjacent_triangle.v1(),
                    v3: adjacent_triangle.v3(),
                    v4: vertex_id,
                },
                adjacent_triangle.neighbor23(),
                adjacent_triangle.neighbor31(),
            )
        } else if adjacent_triangle.neighbor23() == Some(triangle_id) {
            (
                Quad {
                    v1: adjacent_triangle.v3(),
                    v2: adjacent_triangle.v2(),
                    v3: adjacent_triangle.v1(),
                    v4: vertex_id,
                },
                adjacent_triangle.neighbor31(),
                adjacent_triangle.neighbor12(),
            )
        } else {
            (
                Quad {
                    v1: adjacent_triangle.v1(),
                    v2: adjacent_triangle.v3(),
                    v3: adjacent_triangle.v2(),
                    v4: vertex_id,
                },
                adjacent_triangle.neighbor12(),
                adjacent_triangle.neighbor23(),
            )
        };

    // The triangle containing P as a vertex and the unstacked triangle form a convex quadrilateral whose diagonal is drawn in the wrong direction.
    // Swap this diagonal so that two old triangles are replaced by two new triangles and the structure of the Delaunay triangulation is locally restored.
    triangulation::update_triangle_neighbour(
        triangle_3_id,
        Some(adjacent_triangle_id),
        Some(triangle_id),
        triangles,
    );
    triangulation::update_triangle_neighbour(
        triangles[triangle_id].neighbor31(),
        Some(triangle_id),
        Some(adjacent_triangle_id),
        triangles,
    );

    triangles[triangle_id].verts = [quad.v4, quad.v1, quad.v3];

    triangles[adjacent_triangle_id].verts = [quad.v4, quad.v3, quad.v2];

    triangles[adjacent_triangle_id].neighbors = [
        Some(triangle_id),
        triangle_4_id,
        triangles[triangle_id].neighbor31(),
    ];

    triangles[triangle_id].neighbors[EDGE_23] = triangle_3_id;
    triangles[triangle_id].neighbors[EDGE_31] = Some(adjacent_triangle_id);

    (triangle_3_id, triangle_4_id)
}

struct QuadDiag {
    from_triangle: TriangleId,
    to_triangle: TriangleId,
    edge: Edge,
}
impl QuadDiag {
    #[inline]
    fn new(from: usize, to: usize, edge: Edge) -> Self {
        Self {
            from_triangle: from,
            to_triangle: to,
            edge,
        }
    }
}

fn restore_delaunay_triangulation_constrained(
    triangles: &mut Vec<TriangleData>,
    vertices: &mut Vec<Vec2>,
    constrained_edge: &Edge,
    new_diagonals_created: &mut VecDeque<QuadDiag>,
) {
    while let Some(new_diag) = new_diagonals_created.pop_front() {
        if new_diag.edge.undirected_equals(constrained_edge) {
            continue;
        } else {
            let (quad_diag_swapped, quad_diag) = swap_quad_diagonal_to_edge(
                triangles,
                vertices,
                new_diag.edge.from,
                new_diag.from_triangle,
                new_diag.to_triangle,
            );
            if quad_diag_swapped {
                new_diagonals_created.push_back(quad_diag);
            }
        }
    }
}

fn swap_quad_diagonal_to_edge(
    triangles: &mut Vec<TriangleData>,
    vertices: &Vec<Vec2>,
    vertex_id: VertexId,
    triangle_id: TriangleId,
    adjacent_triangle_id: TriangleId,
) -> (bool, QuadDiag) {
    let adjacent_triangle = &triangles[adjacent_triangle_id];
    let (quad, triangle_3_id, triangle_4_id) =
        if adjacent_triangle.neighbor12() == Some(triangle_id) {
            (
                Quad {
                    v1: adjacent_triangle.v2(),
                    v2: adjacent_triangle.v1(),
                    v3: adjacent_triangle.v3(),
                    v4: vertex_id,
                },
                adjacent_triangle.neighbor23(),
                adjacent_triangle.neighbor31(),
            )
        } else if adjacent_triangle.neighbor23() == Some(triangle_id) {
            (
                Quad {
                    v1: adjacent_triangle.v3(),
                    v2: adjacent_triangle.v2(),
                    v3: adjacent_triangle.v1(),
                    v4: vertex_id,
                },
                adjacent_triangle.neighbor31(),
                adjacent_triangle.neighbor12(),
            )
        } else {
            (
                Quad {
                    v1: adjacent_triangle.v1(),
                    v2: adjacent_triangle.v3(),
                    v3: adjacent_triangle.v2(),
                    v4: vertex_id,
                },
                adjacent_triangle.neighbor12(),
                adjacent_triangle.neighbor23(),
            )
        };

    // Check if the vertex is on the circumcircle of the adjacent triangle:
    let mut swapped_quad_diagonal = false;
    if triangulation::is_vertex_in_triangle_circumcircle(
        vertices[quad.v1],
        vertices[quad.v2],
        vertices[quad.v3],
        vertices[vertex_id],
    ) {
        // The triangle containing P as a vertex and the unstacked triangle form a convex quadrilateral whose diagonal is drawn in the wrong direction.
        // Swap this diagonal so that two old triangles are replaced by two new triangles and the structure of the Delaunay triangulation is locally restored.
        triangulation::update_triangle_neighbour(
            triangle_3_id,
            Some(adjacent_triangle_id),
            Some(triangle_id),
            triangles,
        );
        triangulation::update_triangle_neighbour(
            triangles[triangle_id].neighbor31(),
            Some(triangle_id),
            Some(adjacent_triangle_id),
            triangles,
        );

        triangles[triangle_id].verts = [quad.v4, quad.v1, quad.v3];

        triangles[adjacent_triangle_id].verts = [quad.v4, quad.v3, quad.v2];

        triangles[adjacent_triangle_id].neighbors = [
            Some(triangle_id),
            triangle_4_id,
            triangles[triangle_id].neighbor31(),
        ];

        triangles[triangle_id].neighbors[EDGE_23] = triangle_3_id;
        triangles[triangle_id].neighbors[EDGE_31] = Some(adjacent_triangle_id);

        swapped_quad_diagonal = true;
    }
    (
        swapped_quad_diagonal,
        // TODO Justify unwraps
        QuadDiag::new(
            triangle_3_id.unwrap(),
            triangle_4_id.unwrap(),
            Edge::new(quad.v1, quad.v2),
        ),
    )
}

///               q4
///              /  \
///            /  Tn  \
///          /          \
///      q1=e1 -------- e2=q3
///          \          /
///            \  Tc  /
///              \  /
///               q2
/// where:
/// Tn: triangle neighbour id
/// Tc: triangle current id
/// q1,q2,q3,q4: quad coords
fn find_quad_from_triangle_and_crossing_edge(
    triangle_id: TriangleId,
    crossing_edge: &Edge,
    triangles: &Vec<TriangleData>,
) -> (Quad, TriangleId, TriangleId) {
    // Crossing edge vertices
    let crossing_edge_v1 = crossing_edge.from;
    let crossing_edge_v2 = crossing_edge.to;

    // Triangle vertices
    let vertex_1 = triangles[triangle_id].v1();
    let vertex_2 = triangles[triangle_id].v2();
    let vertex_3 = triangles[triangle_id].v3();

    // Quad coords
    let q1 = crossing_edge_v1;
    let q3 = crossing_edge_v2;

    // Shared triangle
    // Crossing edge starts at v1
    let (neighbour_triangle_id, q2) = if crossing_edge_v1 == vertex_1 {
        // Crossing edge ends at v2
        if crossing_edge_v2 == vertex_2 {
            (triangles[triangle_id].neighbor12().unwrap(), vertex_3)
        }
        // Crossing edge ends at v3
        else {
            (triangles[triangle_id].neighbor31().unwrap(), vertex_2)
        }
    }
    // Crossing edge starts at v2
    else if crossing_edge_v1 == vertex_2 {
        // Crossing edge ends at v1
        if crossing_edge_v2 == vertex_1 {
            (triangles[triangle_id].neighbor12().unwrap(), vertex_3)
        }
        // Crossing edge ends at v3
        else {
            (triangles[triangle_id].neighbor31().unwrap(), vertex_1)
        }
    }
    // Crossing edge starts at v3
    else {
        // Crossing edge ends at v1
        if crossing_edge_v2 == vertex_1 {
            (triangles[triangle_id].neighbor31().unwrap(), vertex_2)
        }
        // Crossing edge ends at v2
        else {
            (triangles[triangle_id].neighbor23().unwrap(), vertex_1)
        }
    };

    let v1 = triangles[neighbour_triangle_id].v1();
    let v2 = triangles[neighbour_triangle_id].v2();
    let v3 = triangles[neighbour_triangle_id].v3();
    let q4 = if !(v1 == crossing_edge_v1 || v1 == crossing_edge_v2) {
        v1
    } else if !(v2 == crossing_edge_v1 || v2 == crossing_edge_v2) {
        v2
    } else {
        v3
    };

    //Form the quad
    (
        Quad {
            v1: q1,
            v2: q2,
            v3: q3,
            v4: q4,
        },
        triangle_id,
        neighbour_triangle_id,
    )
}

fn quad_diagonals_intersection_test(
    quad: &Quad,
    vertices: &mut Vec<Vec2>,
) -> EdgesIntersectionResult {
    let v1 = vertices[quad.v1];
    let v2 = vertices[quad.v2];
    let v3 = vertices[quad.v3];
    let v4 = vertices[quad.v4];

    // A quad is convex if diagonals intersect
    egdes_intersect(v1, v3, v2, v4)
}

fn constrained_edge_intersect_triangle(
    triangles: &Vec<TriangleData>,
    vertices: &mut Vec<Vec2>,
    constrained_edge: &Edge,
    current_triangle_id: TriangleId,
    visited_triangles: &Vec<bool>,
) -> (bool, Option<Edge>, Neighbor) {
    // Set the constrained edge vertices and the triangle vertices
    let constrained_edge_start = vertices[constrained_edge.from];
    let constrained_edge_end = vertices[constrained_edge.to];

    let triangle_vertex_1 = vertices[triangles[current_triangle_id].v1()];
    let triangle_vertex_2 = vertices[triangles[current_triangle_id].v2()];
    let triangle_vertex_3 = vertices[triangles[current_triangle_id].v3()];

    // Check if constrained edge cross v1-v2
    if egdes_intersect(
        constrained_edge_start,
        constrained_edge_end,
        triangle_vertex_1,
        triangle_vertex_2,
    ) == EdgesIntersectionResult::Crossing
    {
        let neighbour_triangle_id = triangles[current_triangle_id].neighbor12();
        match neighbour_triangle_id {
            Some(neighbour_triangle_id) => {
                if !visited_triangles[neighbour_triangle_id] {
                    return (
                        true,
                        Some(Edge::new(
                            triangles[current_triangle_id].v1(),
                            triangles[current_triangle_id].v2(),
                        )),
                        triangles[current_triangle_id].neighbor12(),
                    );
                }
            }
            None => return (false, None, None), // TODO: not possible error,
        }
    }
    // Check if constrained edge cross v2-v3
    if egdes_intersect(
        constrained_edge_start,
        constrained_edge_end,
        triangle_vertex_2,
        triangle_vertex_3,
    ) == EdgesIntersectionResult::Crossing
    {
        let neighbour_triangle_id = triangles[current_triangle_id].neighbor23();
        match neighbour_triangle_id {
            Some(neighbour_triangle_id) => {
                if !visited_triangles[neighbour_triangle_id] {
                    return (
                        true,
                        Some(Edge::new(
                            triangles[current_triangle_id].v2(),
                            triangles[current_triangle_id].v3(),
                        )),
                        triangles[current_triangle_id].neighbor23(),
                    );
                }
            }
            None => return (false, None, None), // TODO: not possible error,
        }
    }

    // Check if constrained edge cross v3-v1
    if egdes_intersect(
        constrained_edge_start,
        constrained_edge_end,
        triangle_vertex_3,
        triangle_vertex_1,
    ) == EdgesIntersectionResult::Crossing
    {
        let neighbour_triangle_id = triangles[current_triangle_id].neighbor31();
        match neighbour_triangle_id {
            Some(neighbour_triangle_id) => {
                if !visited_triangles[neighbour_triangle_id] {
                    return (
                        true,
                        Some(Edge::new(
                            triangles[current_triangle_id].v3(),
                            triangles[current_triangle_id].v1(),
                        )),
                        triangles[current_triangle_id].neighbor31(),
                    );
                }
            }
            None => return (false, None, None), // TODO: not possible error,
        }
    }

    // The constrained edge don't cross the triangle
    (false, None, None)
}
