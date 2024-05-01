use bevy::{
    log::info,
    math::{Vec2, Vec3A},
};

use super::triangulation::{self, Quad, TriangleData, TriangleId, VertexId};

/// plane_normal must be normalized
/// vertices must all belong to a 3d plane
pub fn triangulate_3d_planar_vertices_constrained(
    vertices: &Vec<[f32; 3]>,
    plane_normal: Vec3A,
    mut constrained_edges: &Vec<(usize, usize)>,
) -> Vec<VertexId> {
    // TODO See what we need for input data format of `triangulate`
    let mut vertices_data = Vec::with_capacity(vertices.len());
    for v in vertices {
        vertices_data.push(Vec3A::from_array(*v));
    }

    let mut planar_vertices =
        triangulation::transform_to_2d_planar_coordinate_system(&mut vertices_data, plane_normal);

    // Delaunay triangulation
    triangulate_2d_vertices_constrained(&mut planar_vertices, &mut constrained_edges)
}

fn triangulate_2d_vertices_constrained(
    vertices: &mut Vec<Vec2>,
    constrained_edges: &Vec<(usize, usize)>,
) -> Vec<VertexId> {
    // Uniformly scale the coordinates of the points so that they all lie between 0 and 1.
    triangulation::normalize_vertices_coordinates(vertices);

    // Sort points into bins. Cover the region to be triangulated by a rectangular grid so that each bin contains roughly N^(1/2) points. Label the bins so that consecutive bins are adjacent to one another, and then allocate each point to its appropriate bin. Sort the list of points in ascending sequence of their bin numbers so that consecutive points are grouped together in the x-y plane.
    let partitioned_vertices = triangulation::VertexBinSort::sort(&vertices, 0.5);

    // Select three dummy points to form a supertriangle that completely encompasses all of the points to be triangulated. This supertriangle initially defines a Delaunay triangulation which is comprised of a single triangle. Its vertices are defined in terms of normalized coordinates and are usually located at a considerable distance from the window which encloses the set of points.
    // TODO Clean: constants + comments
    // TODO See if we could avoid merging fake data
    let mut triangles = Vec::<TriangleData>::new();
    let container_triangle = TriangleData {
        v1: vertices.len(),
        v2: vertices.len() + 1,
        v3: vertices.len() + 2,
        edge12: None,
        edge23: None,
        edge31: None,
    };
    triangles.push(container_triangle.clone());
    vertices.push(Vec2::new(-100., -100.));
    vertices.push(Vec2::new(0., 100.));
    vertices.push(Vec2::new(100., -100.));

    // Id of the triangle we are looking at
    let mut triangle_id = Some(0);

    // Loop over all the input vertices
    for sorted_vertex in partitioned_vertices.iter() {
        // Find an existing triangle which encloses P
        match triangulation::search_enclosing_triangle(
            sorted_vertex.vertex,
            triangle_id,
            &triangles,
            &vertices,
        ) {
            Some(enclosing_triangle_id) => {
                // Delete this triangle and form three new triangles by connecting P to each of its vertices.
                triangulation::split_triangle_in_three_at_vertex(
                    enclosing_triangle_id,
                    sorted_vertex.original_id,
                    &mut triangles,
                    &vertices,
                );
                triangle_id = Some(triangles.len() - 1);
            }
            None => (),
        }
    }

    //CONSTRAINED: --------------------------------------------------------------------------------------------------------------------------------

    //Map each vertices to a triangle that contains it
    let mut triangle_vertices = vec![0; vertices.len()];
    for (indexe, triangle) in triangles.iter().enumerate() {
        triangle_vertices[triangle.v1] = indexe;
        triangle_vertices[triangle.v2] = indexe;
        triangle_vertices[triangle.v3] = indexe;
    }

    // Loop over each constrained edge
    for constrained_edge in constrained_edges {
        //If the constrained is already present in the triangulation, continue
        if is_constrained_edge_inside_triangulation(&triangles, constrained_edge) {
            continue;
        };

        //Store all of the edges that cross the constraine edge
        let intersected_edges =
            intersected_edges(&triangles, vertices, constrained_edge, &triangle_vertices);

        // Remove intersecting edges
        let mut new_edges_created = remove_crossing_edges(
            &mut triangles,
            vertices,
            constrained_edge,
            intersected_edges,
        );

        // Restore Delaunay triangulation
        restore_delaunay_triangulation_constrained(
            &mut triangles,
            vertices,
            constrained_edge,
            &mut new_edges_created,
        );
    }

    //CONSTRAINED: --------------------------------------------------------------------------------------------------------------------------------

    // TODO Clean: Size approx
    let mut indices = Vec::with_capacity(3 * triangles.len());
    for triangle in triangles.iter() {
        if triangle.v1 != container_triangle.v1
            && triangle.v2 != container_triangle.v1
            && triangle.v3 != container_triangle.v1
            && triangle.v1 != container_triangle.v2
            && triangle.v2 != container_triangle.v2
            && triangle.v3 != container_triangle.v2
            && triangle.v1 != container_triangle.v3
            && triangle.v2 != container_triangle.v3
            && triangle.v3 != container_triangle.v3
        {
            indices.push(triangle.v1);
            indices.push(triangle.v2);
            indices.push(triangle.v3);
        }
    }
    indices
}

fn is_constrained_edge_inside_triangulation(
    triangles: &Vec<TriangleData>,
    constrained_edge: &(usize, usize),
) -> bool {
    for triangle in triangles {
        if (triangle.v1 == constrained_edge.0
            || triangle.v2 == constrained_edge.0
            || triangle.v3 == constrained_edge.0)
            && (triangle.v1 == constrained_edge.1
                || triangle.v2 == constrained_edge.1
                || triangle.v3 == constrained_edge.1)
        {
            return true;
        }
    }
    false
}

fn intersected_edges(
    triangles: &Vec<TriangleData>,
    vertices: &mut Vec<Vec2>,
    constrained_edge: &(usize, usize),
    triangle_vertices: &Vec<usize>,
) -> Vec<(TriangleId, (VertexId, VertexId))> {
    // List of every triangles already checked
    //TODO: change data structure
    let mut visited_triangles = vec![false; triangles.len()];

    let mut intersected_edges = Vec::new();

    // providing a starting triangle to begin the search for the edges which cross the constrained edge
    let mut current_triangle_id = triangle_vertices[constrained_edge.0];
    visited_triangles[current_triangle_id] = true;

    // Flag to check whenever we are at the end of the constrained edge
    let mut arrived_to_constrained_edge_end = false;

    // we circle constrained_edge.0
    while !arrived_to_constrained_edge_end {
        info!("intersected_edges");
        let (crossing, crossing_edge, crossed_triangle) = constrained_edge_intersect_triangle(
            triangles,
            vertices,
            constrained_edge,
            current_triangle_id,
            &visited_triangles,
        );

        // Vertices of the current triangle
        let triangle_vertex_1 = vertices[triangles[current_triangle_id].v1];
        let triangle_vertex_2 = vertices[triangles[current_triangle_id].v2];
        let triangle_vertex_3 = vertices[triangles[current_triangle_id].v3];

        //march from one triangle to the next if there are no crossing edges
        //TODO: cut in two different algorithms:research and follow
        if !crossing {
            if triangle_vertex_1 == vertices[constrained_edge.0] {
                match triangles[current_triangle_id].edge31 {
                    // By default, we turn to the left
                    Some(neighbour_triangle_id) => {
                        // if the left triangle exist
                        if !visited_triangles[neighbour_triangle_id] {
                            // if the left triangle hasn't been visited yet
                            current_triangle_id = neighbour_triangle_id; // we got to this triangle
                            visited_triangles[current_triangle_id] = true;
                        } else {
                            // if the left triangle has been visited
                            current_triangle_id = triangles[current_triangle_id].edge12.unwrap(); // we turn to the right triangle
                            visited_triangles[current_triangle_id] = true;
                        }
                    }
                    None => {
                        // If there is no triangle in the left direction, we got to the right direction
                        current_triangle_id = triangles[current_triangle_id].edge12.unwrap();
                        visited_triangles[current_triangle_id] = true;
                    }
                }
            } else if triangle_vertex_2 == vertices[constrained_edge.0] {
                match triangles[current_triangle_id].edge12 {
                    // By default, we turn to the left
                    Some(neighbour_triangle_id) => {
                        // if the left triangle exist
                        if !visited_triangles[neighbour_triangle_id] {
                            // if the left triangle hasn't been visited yet
                            current_triangle_id = neighbour_triangle_id; // we got to this triangle
                            visited_triangles[current_triangle_id] = true;
                        } else {
                            // if the left triangle has been visited
                            current_triangle_id = triangles[current_triangle_id].edge23.unwrap(); // we turn to the right triangle
                            visited_triangles[current_triangle_id] = true;
                        }
                    }
                    None => {
                        // If there is no triangle in the left direction, we got to the right direction
                        current_triangle_id = triangles[current_triangle_id].edge23.unwrap();
                        visited_triangles[current_triangle_id] = true;
                    }
                }
            } else if triangle_vertex_3 == vertices[constrained_edge.0] {
                match triangles[current_triangle_id].edge23 {
                    // By default, we turn to the left
                    Some(neighbour_triangle_id) => {
                        // if the left triangle exist
                        if !visited_triangles[neighbour_triangle_id] {
                            // if the left triangle hasn't been visited yet
                            current_triangle_id = neighbour_triangle_id; // we got to this triangle
                            visited_triangles[current_triangle_id] = true;
                        } else {
                            // if the left triangle has been visited
                            current_triangle_id = triangles[current_triangle_id].edge31.unwrap(); // we turn to the right triangle
                            visited_triangles[current_triangle_id] = true;
                        }
                    }
                    None => {
                        // If there is no triangle in the left direction, we got to the right direction
                        current_triangle_id = triangles[current_triangle_id].edge31.unwrap();
                        visited_triangles[current_triangle_id] = true;
                    }
                }
            }
        }
        // Add edge if we find one crossing edge and march from one triangle to the next in the general direction of constrained_edge.1,
        else if crossing {
            current_triangle_id = crossed_triangle.unwrap();
            visited_triangles[current_triangle_id] = true; // not necessary since we will laways go to a new triangle here
            intersected_edges.push((current_triangle_id, crossing_edge.unwrap()));

            // Check if we are at constrained_edge.1
            if crossing_edge.unwrap().0 == constrained_edge.1
                || crossing_edge.unwrap().1 == constrained_edge.1
            {
                arrived_to_constrained_edge_end = true;
            }
        }
    }

    intersected_edges
}

fn remove_crossing_edges(
    triangles: &mut Vec<TriangleData>,
    vertices: &mut Vec<Vec2>,
    constrained_edge: &(usize, usize),
    mut intersected_edges: Vec<(usize, (usize, usize))>,
) -> Vec<(usize, usize, (usize, usize))> {
    let mut new_edges_created: Vec<(usize, usize, (usize, usize))> = Vec::new();

    while intersected_edges.len() != 0 {
        info!("remove_crossing_edges");
        let (current_triangle_id, current_crossing_edge) = intersected_edges.pop().unwrap();

        // Construct a quad from the current triangle and its edge by looking at the neighbour triangle of this edge
        let (quad, triangle_id, adjacent_triangle_id) = find_quad_from_triangle_and_crossing_edge(
            current_triangle_id,
            &current_crossing_edge,
            triangles,
        );

        // If the quad is not convex, we skip this edge
        if !is_quad_convex(&quad, vertices) {
            continue;
        }
        // swap the diagonal of this strictly convex quadrilateral
        else {
            let (t3, t4) =
                swap_quad_diagonal(quad.v4, triangle_id, adjacent_triangle_id, triangles);

            //  If the new diagonal still intersects the constrained edge, then place it on the list of intersecting edges
            let quad_vertex_1 = vertices[quad.v1];
            let quad_vertex_4 = vertices[quad.v4];
            let constrained_edge_vertex_1 = vertices[constrained_edge.0];
            let constrained_edge_vertex_2 = vertices[constrained_edge.1];

            if egdes_intersect(
                quad_vertex_1,
                quad_vertex_4,
                constrained_edge_vertex_1,
                constrained_edge_vertex_2,
            ) {
                intersected_edges.push((t3.unwrap(), (quad.v1, quad.v4)))
            }
            // If the edge does not intersect the constrained edge, then place it on a list of newly created edges
            else {
                new_edges_created.push((t3.unwrap(), t4.unwrap(), (quad.v1, quad.v4)));
            }
        }
    }

    new_edges_created
}

fn swap_quad_diagonal(
    vertex_id: VertexId,
    triangle_id: TriangleId,
    adjacent_triangle_id: TriangleId,
    triangles: &mut Vec<TriangleData>,
) -> (Option<TriangleId>, Option<TriangleId>) {
    let adjacent_triangle = &triangles[adjacent_triangle_id];
    let (quad, triangle_3_id, triangle_4_id) = if adjacent_triangle.edge12 == Some(triangle_id) {
        (
            Quad {
                v1: adjacent_triangle.v2,
                v2: adjacent_triangle.v1,
                v3: adjacent_triangle.v3,
                v4: vertex_id,
            },
            adjacent_triangle.edge23,
            adjacent_triangle.edge31,
        )
    } else if adjacent_triangle.edge23 == Some(triangle_id) {
        (
            Quad {
                v1: adjacent_triangle.v3,
                v2: adjacent_triangle.v2,
                v3: adjacent_triangle.v1,
                v4: vertex_id,
            },
            adjacent_triangle.edge31,
            adjacent_triangle.edge12,
        )
    } else {
        (
            Quad {
                v1: adjacent_triangle.v1,
                v2: adjacent_triangle.v3,
                v3: adjacent_triangle.v2,
                v4: vertex_id,
            },
            adjacent_triangle.edge12,
            adjacent_triangle.edge23,
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
        triangles[triangle_id].edge31,
        Some(triangle_id),
        Some(adjacent_triangle_id),
        triangles,
    );

    triangles[triangle_id].v1 = quad.v4;
    triangles[triangle_id].v2 = quad.v1;
    triangles[triangle_id].v3 = quad.v3;

    triangles[adjacent_triangle_id].v1 = quad.v4;
    triangles[adjacent_triangle_id].v2 = quad.v3;
    triangles[adjacent_triangle_id].v3 = quad.v2;

    triangles[adjacent_triangle_id].edge12 = Some(triangle_id);
    triangles[adjacent_triangle_id].edge23 = triangle_4_id;
    triangles[adjacent_triangle_id].edge31 = triangles[triangle_id].edge31;

    triangles[triangle_id].edge23 = triangle_3_id;
    triangles[triangle_id].edge31 = Some(adjacent_triangle_id);

    (triangle_3_id, triangle_4_id)
}

fn restore_delaunay_triangulation_constrained(
    triangles: &mut Vec<TriangleData>,
    vertices: &mut Vec<Vec2>,
    constrained_edge: &(usize, usize),
    new_edges_created: &mut Vec<(usize, usize, (usize, usize))>,
) {
    while new_edges_created.len() != 0 {
        let (current_triangle_id_1, current_triangle_id_2, current_edge) =
            new_edges_created.pop().unwrap();
        // If the edge is equal to the constrained edge, then skip
        if (current_edge.0 == constrained_edge.0 && current_edge.1 == constrained_edge.1)
            || (current_edge.0 == constrained_edge.1 && current_edge.1 == constrained_edge.0)
        {
            break;
        } else {
            let (quad_diag_swapped, t3, t4, edge_1, edge_2) = swap_quad_diagonal_to_edge(
                current_edge.0,
                current_triangle_id_1,
                current_triangle_id_2,
                triangles,
                vertices,
            );

            if quad_diag_swapped {
                new_edges_created.push((t3.unwrap(), t4.unwrap(), (edge_1, edge_2)));
            }
        }
    }
}

fn swap_quad_diagonal_to_edge(
    vertex_id: VertexId,
    triangle_id: TriangleId,
    adjacent_triangle_id: TriangleId,
    triangles: &mut Vec<TriangleData>,
    vertices: &Vec<Vec2>,
) -> (bool, Option<TriangleId>, Option<TriangleId>, usize, usize) {
    let adjacent_triangle = &triangles[adjacent_triangle_id];
    let (quad, triangle_3_id, triangle_4_id) = if adjacent_triangle.edge12 == Some(triangle_id) {
        (
            Quad {
                v1: adjacent_triangle.v2,
                v2: adjacent_triangle.v1,
                v3: adjacent_triangle.v3,
                v4: vertex_id,
            },
            adjacent_triangle.edge23,
            adjacent_triangle.edge31,
        )
    } else if adjacent_triangle.edge23 == Some(triangle_id) {
        (
            Quad {
                v1: adjacent_triangle.v3,
                v2: adjacent_triangle.v2,
                v3: adjacent_triangle.v1,
                v4: vertex_id,
            },
            adjacent_triangle.edge31,
            adjacent_triangle.edge12,
        )
    } else {
        (
            Quad {
                v1: adjacent_triangle.v1,
                v2: adjacent_triangle.v3,
                v3: adjacent_triangle.v2,
                v4: vertex_id,
            },
            adjacent_triangle.edge12,
            adjacent_triangle.edge23,
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
            triangles[triangle_id].edge31,
            Some(triangle_id),
            Some(adjacent_triangle_id),
            triangles,
        );

        triangles[triangle_id].v1 = quad.v4;
        triangles[triangle_id].v2 = quad.v1;
        triangles[triangle_id].v3 = quad.v3;

        triangles[adjacent_triangle_id].v1 = quad.v4;
        triangles[adjacent_triangle_id].v2 = quad.v3;
        triangles[adjacent_triangle_id].v3 = quad.v2;

        triangles[adjacent_triangle_id].edge12 = Some(triangle_id);
        triangles[adjacent_triangle_id].edge23 = triangle_4_id;
        triangles[adjacent_triangle_id].edge31 = triangles[triangle_id].edge31;

        triangles[triangle_id].edge23 = triangle_3_id;
        triangles[triangle_id].edge31 = Some(adjacent_triangle_id);

        swapped_quad_diagonal = true;
    }
    (
        swapped_quad_diagonal,
        triangle_3_id,
        triangle_4_id,
        quad.v1,
        quad.v2,
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
///               q1
/// where:
/// Tn: triangle neighbour id
/// Tc: triangle current id
/// q1,q2,q3,q4: quad coords
fn find_quad_from_triangle_and_crossing_edge(
    triangle_id: usize,
    crossing_edge: &(usize, usize),
    triangles: &Vec<TriangleData>,
) -> (Quad, usize, usize) {
    // Crossing edge vertices
    let crossing_edge_v1 = crossing_edge.0;
    let crossing_edge_v2 = crossing_edge.1;

    // Triangle vertices
    let vertex_1 = triangles[triangle_id].v1;
    let vertex_2 = triangles[triangle_id].v2;
    let vertex_3 = triangles[triangle_id].v3;

    // Quad coords
    let q1 = crossing_edge_v1;
    let q2 = crossing_edge_v2;

    // Shared triangle
    // Crossing edge starts at v1
    let (neighbour_triangle_id, q3) = if crossing_edge_v1 == vertex_1 {
        // Crossing edge ends at v2
        if crossing_edge_v2 == vertex_2 {
            (triangles[triangle_id].edge12.unwrap(), vertex_3)
        }
        // Crossing edge ends at v3
        else {
            (triangles[triangle_id].edge31.unwrap(), vertex_2)
        }
    }
    // Crossing edge starts at v2
    else if crossing_edge_v1 == vertex_2 {
        // Crossing edge ends at v1
        if crossing_edge_v2 == vertex_1 {
            (triangles[triangle_id].edge12.unwrap(), vertex_3)
        }
        // Crossing edge ends at v3
        else {
            (triangles[triangle_id].edge31.unwrap(), vertex_1)
        }
    }
    // Crossing edge starts at v3
    else {
        // Crossing edge ends at v1
        if crossing_edge_v2 == vertex_1 {
            (triangles[triangle_id].edge31.unwrap(), vertex_2)
        }
        // Crossing edge ends at v2
        else {
            (triangles[triangle_id].edge23.unwrap(), vertex_1)
        }
    };

    // Neighbour triangle vertices
    // let neighbour_triangle_vertex_1 = triangles[neighbour_triangle_id].v1;
    // let neighbour_triangle_vertex_2 = triangles[neighbour_triangle_id].v2;
    // let neighbour_triangle_vertex_3 = triangles[neighbour_triangle_id].v3;

    // // Looking for the last quad corrd which is not the vertex from the shared edge of the two triangles
    // let mut neighbour_vertices_id = Vec::new();
    // neighbour_vertices_id.push(neighbour_triangle_vertex_1);
    // neighbour_vertices_id.push(neighbour_triangle_vertex_2);
    // neighbour_vertices_id.push(neighbour_triangle_vertex_3);

    // let mut q4 = triangles[neighbour_triangle_id].v1;
    // for vertex_id in vec![
    //     triangles[neighbour_triangle_id].v1,
    //     triangles[neighbour_triangle_id].v2,
    //     triangles[neighbour_triangle_id].v3,
    // ] {
    //     if !(vertex_id == crossing_edge_v1 || vertex_id == crossing_edge_v2) {
    //         q4 = vertex_id;
    //         break;
    //     }
    // }

    let v1 = triangles[neighbour_triangle_id].v1;
    let v2 = triangles[neighbour_triangle_id].v2;
    let v3 = triangles[neighbour_triangle_id].v3;
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

fn is_quad_convex(quad: &Quad, vertices: &mut Vec<Vec2>) -> bool {
    let v1 = vertices[quad.v1];
    let v2 = vertices[quad.v2];
    let v3 = vertices[quad.v3];
    let v4 = vertices[quad.v4];

    // A quad is convex if diagonals intersect
    egdes_intersect(v1, v4, v2, v3)
}

fn constrained_edge_intersect_triangle(
    triangles: &Vec<TriangleData>,
    vertices: &mut Vec<Vec2>,
    constrained_edge: &(usize, usize),
    current_triangle_id: usize,
    visited_triangles: &Vec<bool>,
) -> (bool, Option<(usize, usize)>, Option<usize>) {
    // Set the constrained edge vertices and the triangle vertices
    let constrained_edge_start = vertices[constrained_edge.0];
    let constrained_edge_end = vertices[constrained_edge.1];
    let triangle_vertex_1 = vertices[triangles[current_triangle_id].v1];
    let triangle_vertex_2 = vertices[triangles[current_triangle_id].v2];
    let triangle_vertex_3 = vertices[triangles[current_triangle_id].v3];

    // Check if constrained edge cross v1-v2
    if egdes_intersect(
        constrained_edge_start,
        constrained_edge_end,
        triangle_vertex_1,
        triangle_vertex_2,
    ) {
        let neighbour_triangle_id = triangles[current_triangle_id].edge12;
        match neighbour_triangle_id {
            Some(neighbour_triangle_id) => {
                if !visited_triangles[neighbour_triangle_id] {
                    return (
                        true,
                        Some((
                            triangles[current_triangle_id].v1,
                            triangles[current_triangle_id].v2,
                        )),
                        triangles[current_triangle_id].edge12,
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
    ) {
        let neighbour_triangle_id = triangles[current_triangle_id].edge23;
        match neighbour_triangle_id {
            Some(neighbour_triangle_id) => {
                if !visited_triangles[neighbour_triangle_id] {
                    return (
                        true,
                        Some((
                            triangles[current_triangle_id].v2,
                            triangles[current_triangle_id].v3,
                        )),
                        triangles[current_triangle_id].edge23,
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
    ) {
        let neighbour_triangle_id = triangles[current_triangle_id].edge31;
        match neighbour_triangle_id {
            Some(neighbour_triangle_id) => {
                if !visited_triangles[neighbour_triangle_id] {
                    return (
                        true,
                        Some((
                            triangles[current_triangle_id].v3,
                            triangles[current_triangle_id].v1,
                        )),
                        triangles[current_triangle_id].edge31,
                    );
                }
            }
            None => return (false, None, None), // TODO: not possible error,
        }
    }
    // The constrained edge don't cross the triangle

    (false, None, None)
}

fn egdes_intersect(e1_0: Vec2, e1_1: Vec2, e2_0: Vec2, e2_1: Vec2) -> bool {
    edges_intersect_internal(e1_0, e1_1, e2_0, e2_1, false)
}

fn edges_intersect_internal(
    e1_0: Vec2,
    e1_1: Vec2,
    e2_0: Vec2,
    e2_1: Vec2,
    include_shared_end_points: bool,
) -> bool {
    let e1_12 = Vec2::new(e1_1.x - e1_0.x, e1_1.y - e1_0.y);
    let e2_12 = Vec2::new(e2_1.x - e2_0.x, e2_1.y - e2_0.y);

    // If any of the vertices are shared between the two diagonals,
    // the quad collapses into a triangle and is convex by default.
    if e1_0 == e2_0 || e1_0 == e2_1 || e1_1 == e2_0 || e1_1 == e2_1 {
        include_shared_end_points
    } else {
        // Compute cross product between each point and the opposite diagonal
        // Look at sign of the Z component to see which side of line point is on
        let e1_0_cross_e2 = (e1_0.x - e2_0.x) * e2_12.y - (e1_0.y - e2_0.y) * e2_12.x;
        let e1_1_cross_e2 = (e1_1.x - e2_0.x) * e2_12.y - (e1_1.y - e2_0.y) * e2_12.x;
        let e2_0_cross_e1 = (e2_0.x - e1_0.x) * e1_12.y - (e2_0.y - e1_0.y) * e1_12.x;
        let e2_1_cross_e1 = (e2_1.x - e1_0.x) * e1_12.y - (e2_1.y - e1_0.y) * e1_12.x;

        // Check that the points for each diagonal lie on opposite sides of the other
        // diagonal. Quad is also convex if a1/a2 lie on b1->b2 (and vice versa) since
        // the shape collapses into a triangle (hence >= instead of >)
        return ((e1_0_cross_e2 >= 0. && e1_1_cross_e2 <= 0.)
            || (e1_0_cross_e2 <= 0. && e1_1_cross_e2 >= 0.))
            && ((e2_0_cross_e1 >= 0. && e2_1_cross_e1 <= 0.)
                || (e2_0_cross_e1 <= 0. && e2_1_cross_e1 >= 0.));
    }
}
