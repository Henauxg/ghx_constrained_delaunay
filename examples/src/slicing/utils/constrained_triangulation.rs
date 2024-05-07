use std::collections::{HashMap, HashSet, VecDeque};

use bevy::{
    log::info,
    math::{Vec2, Vec3A},
};

use super::triangulation::{
    self, wrap_and_triangulate_2d_vertices, Quad, TriangleData, TriangleId, VertexId,
};

/// plane_normal must be normalized
/// vertices must all belong to a 3d plane
/// constrained edges must be oriented: vi -> vj
pub fn triangulate_3d_planar_vertices_constrained(
    vertices: &Vec<[f32; 3]>,
    plane_normal: Vec3A,
    constrained_edges: &mut VecDeque<(usize, usize)>,
) -> (Vec<VertexId>, Vec<Vec<TriangleData>>) {
    // TODO See what we need for input data format of `triangulate`
    let mut vertices_data = Vec::with_capacity(vertices.len());
    for v in vertices {
        vertices_data.push(Vec3A::from_array(*v));
    }

    let mut planar_vertices =
        triangulation::transform_to_2d_planar_coordinate_system(&mut vertices_data, plane_normal);

    // Delaunay triangulation
    triangulate_2d_vertices_constrained(&mut planar_vertices, constrained_edges)
}

fn triangulate_2d_vertices_constrained(
    vertices: &mut Vec<Vec2>,
    constrained_edges: &mut VecDeque<(usize, usize)>,
) -> (Vec<VertexId>, Vec<Vec<TriangleData>>) {
    let (mut triangles, container_triangle, mut debugger) =
        wrap_and_triangulate_2d_vertices(vertices);

    apply_constrains(vertices, constrained_edges, &mut triangles, &mut debugger);

    let killed_triangles = marked_triangles_to_be_killed(&triangles, constrained_edges);

    info!("killed_triangles {:?}", killed_triangles);

    let indices = remove_constrained_wrapping(
        &triangles,
        &container_triangle,
        &killed_triangles,
        &mut debugger,
    );

    (indices, debugger)
}

fn remove_constrained_wrapping(
    triangles: &Vec<TriangleData>,
    container_triangle: &TriangleData,
    killed_triangles: &Vec<bool>,
    debugger: &mut Vec<Vec<TriangleData>>,
) -> Vec<VertexId> {
    // TODO Clean: Size approx
    let mut indices = Vec::with_capacity(3 * triangles.len());
    let mut filtered_debug_triangles = Vec::new();
    for (index, triangle) in triangles.iter().enumerate() {
        if triangle.v1 != container_triangle.v1
            && triangle.v2 != container_triangle.v1
            && triangle.v3 != container_triangle.v1
            && triangle.v1 != container_triangle.v2
            && triangle.v2 != container_triangle.v2
            && triangle.v3 != container_triangle.v2
            && triangle.v1 != container_triangle.v3
            && triangle.v2 != container_triangle.v3
            && triangle.v3 != container_triangle.v3
            && !killed_triangles[index]
        {
            indices.push(triangle.v1);
            indices.push(triangle.v2);
            indices.push(triangle.v3);
            filtered_debug_triangles.push(triangle.clone());
        }
    }
    debugger.push(filtered_debug_triangles);

    indices
}

/// Mark as repmove all triangles that are not suppose tostay in theoriginal configuration
/// This supposes that all group of constrained edges are oriented in different orientations, as shown:
/// ---------------<----------------
/// |                              |
/// |   ----------->------------   |
/// |   |  zone to be removed  |   |
/// |   |   -------<--------   |   |
/// |   |   |              |   |   |
/// v   ∧   v              ∧   v   ∧
/// |   |   |              |   |   |
/// |   |   ------->--------   |   |
/// |   |                      |   |
/// |   -----------<------------   |
/// |                              |
/// --------------->----------------
fn marked_triangles_to_be_killed(
    triangles: &Vec<TriangleData>,
    constrained_edges: &mut VecDeque<(usize, usize)>,
) -> Vec<bool> {
    let mut killed_triangles = vec![true; triangles.len()];
    let mut visited_triangles = vec![false; triangles.len()];

    visited_triangles[0] = true;

    let mut triangle_buffer = Vec::new();

    for (index, triangle) in triangles.iter().enumerate() {
        if visited_triangles[index] {
            continue;
        }

        let v1 = triangle.v1;
        let v2 = triangle.v2;
        let v3 = triangle.v3;

        let e12: bool = constrained_edges.contains(&(v1, v2));
        let e23 = constrained_edges.contains(&(v2, v3));
        let e31 = constrained_edges.contains(&(v3, v1));

        if e12 || e23 || e31 {
            killed_triangles[index] = false;

            triangle_buffer.clear();

            if !e12 {
                triangle_buffer.push(triangle.edge12);
            }
            if !e23 {
                triangle_buffer.push(triangle.edge23);
            }
            if !e31 {
                triangle_buffer.push(triangle.edge31);
            }

            while triangle_buffer.len() > 0 {
                let triangle_index = triangle_buffer.pop().unwrap();

                match triangle_index {
                    Some(neighbor) => {
                        if visited_triangles[neighbor] {
                            continue;
                        } else {
                            killed_triangles[neighbor] = false;
                            visited_triangles[neighbor] = true;

                            let v1 = triangles[neighbor].v1;
                            let v2 = triangles[neighbor].v2;
                            let v3 = triangles[neighbor].v3;

                            if !constrained_edges.contains(&(v1, v2)) {
                                triangle_buffer.push(triangles[neighbor].edge12);
                            }
                            if !constrained_edges.contains(&(v2, v3)) {
                                triangle_buffer.push(triangles[neighbor].edge23);
                            }
                            if !constrained_edges.contains(&(v3, v1)) {
                                triangle_buffer.push(triangles[neighbor].edge31);
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

    killed_triangles
}

// type FamilyId = usize;

// fn get_sets_of_triangles(
//     constrained_edges: &mut VecDeque<(usize, usize)>,
//     vertices: &mut Vec<Vec2>,
// ) -> Vec<(FamilyId, Vec<(usize, usize)>)> {
//     // buffer of all the constrained cycles
//     let mut sets = VecDeque::new();

//     // fill the buffer with all the cycles
//     while constrained_edges.len() != 0 {
//         sets.push_back(determined_set(constrained_edges));
//     }

//     let mut polygons: Vec<(FamilyId, Vec<(usize, usize)>)> = Vec::new();

//     // if there is only one cycle, the parent tree is only the set present in the buffer
//     if sets.len() == 1 {
//         polygons.push((usize::MAX, sets.pop_front().unwrap()));
//         return polygons;
//     } else if sets.len() == 2 {
//         let polygon_1: Vec<(usize, usize)> = sets[0].clone();
//         let polygon_2 = sets[1].clone();

//         let family_trait = get_polygons_relationship(&polygon_1, &polygon_2, vertices);
//         if family_trait == Family::PARENT {
//             polygons.push((usize::MAX, polygon_1));
//             polygons.push((1, polygon_2));
//         } else {
//             polygons.push((usize::MAX, polygon_1));
//             polygons.push((usize::MAX, polygon_2));
//         }
//     } else {
//         // we compare every sets 2 by 2
//         for index in 0..sets.len() - 1 {
//             // polygons to be compare
//             let polygon_1 = sets[index].clone();
//             let polygon_2 = sets[index + 1].clone();

//             // we check rather the polygons are parent-child related or twins (not inside of each others but in the same level in thee family tree)
//             let family_trait = get_polygons_relationship(&polygon_1, &polygon_2, vertices);

//             //add the polygon into the polygon buffer with the family index
//             if index == 0 {
//                 info!("index start");
//                 // specific case for start
//                 if family_trait == Family::PARENT {
//                     polygons.push((usize::MAX, polygon_1));
//                 } else {
//                     polygons.push((usize::MAX, polygon_1));
//                     polygons.push((usize::MAX, polygon_2));
//                     //todo
//                 }
//             } else if index == sets.len() - 2 {
//                 info!("index end");
//                 // specific case for end
//                 if family_trait == Family::PARENT {
//                     polygons.push((index, polygon_1));
//                     polygons.push((index + 1, polygon_2));
//                 } else {
//                     polygons.push((index, polygon_1));
//                     polygons.push((index, polygon_2));
//                 }
//                 break;
//             } else {
//                 if family_trait == Family::PARENT {
//                     polygons.push((index, polygon_1));
//                 } else {
//                     polygons.push((index, polygon_1));
//                     polygons.push((index, polygon_2));
//                     //todo
//                 }
//             }
//         }
//     }

//     polygons
// }

// #[derive(PartialEq)]
// pub enum Family {
//     PARENT,
//     TWINS,
// }

// fn get_polygons_relationship(
//     polygon_1: &Vec<(usize, usize)>,
//     polygon_2: &Vec<(usize, usize)>,
//     vertices: &mut Vec<Vec2>,
// ) -> Family {
//     //source: https://stackoverflow.com/questions/4833802/check-if-polygon-is-inside-a-polygon
//     // we need to check rather a set is inside another or if the two set are all children of the upper set
//     // if one of the two sets is not convex, we need to check every edge intersection
//     if !is_polygon_convex(&polygon_1, vertices) || !is_polygon_convex(&polygon_2, vertices) {
//         // If really necessary, you can speed up the line intersection tests using the sweep line algorithm.
//         if !are_polygons_intersecting(&polygon_1, &polygon_2, vertices) {
//             // polygon 2 is inside polygon 1
//             if is_polygon_inside_polygon(&polygon_2, &polygon_1, vertices) {
//                 return Family::PARENT;
//             }
//             // both polygons are children
//             else {
//                 return Family::TWINS;
//             }
//         }
//     }

//     // if the polygons are convex or there are no edge intersection, we need to check if a polygon is inside another by looking at every points of it
//     // polygon 2 is inside polygon 1
//     if is_polygon_inside_polygon(&polygon_2, &polygon_1, vertices) {
//         return Family::PARENT;
//     }
//     // both polygons are children
//     else {
//         return Family::TWINS;
//     }
// }

// ///source: https://stackoverflow.com/questions/217578/how-can-i-determine-whether-a-2d-point-is-within-a-polygon
// fn is_polygon_inside_polygon(
//     tested_polygon: &Vec<(usize, usize)>,
//     encompass_polygon: &Vec<(usize, usize)>,
//     vertices: &mut Vec<Vec2>,
// ) -> bool {
//     let mut c = false;
//     for edge_index_1 in tested_polygon {
//         let p = vertices[edge_index_1.0];

//         for edge_index_2 in encompass_polygon.iter() {
//             let q = vertices[edge_index_2.0];
//             let r = vertices[edge_index_2.1];

//             if ((r.y > p.y) != (q.y > p.y)) && (p.x < (q.x - r.x) * (p.y - r.y) / (q.y - r.y) + r.x)
//             {
//                 c = !c;
//             }
//         }
//     }

//     true
// }

// fn are_polygons_intersecting(
//     polygon_1: &Vec<(usize, usize)>,
//     polygon_2: &Vec<(usize, usize)>,
//     vertices: &mut Vec<Vec2>,
// ) -> bool {
//     // for all edges of the first polygon
//     for edge_index_1 in polygon_1 {
//         for edge_index_2 in polygon_2 {
//             let e1_0 = vertices[edge_index_1.0];
//             let e1_1 = vertices[edge_index_1.1];
//             let e2_0 = vertices[edge_index_2.0];
//             let e2_1 = vertices[edge_index_2.1];

//             if egdes_intersect(e1_0, e1_1, e2_0, e2_1) != EdgesIntersectionResult::None {
//                 return true;
//             }
//         }
//     }

//     false
// }

// fn determined_set(constrained_edges: &mut VecDeque<(usize, usize)>) -> Vec<(usize, usize)> {
//     //buffer of a cycle
//     let mut set = Vec::new();

//     // get the initial constrained edge to start with
//     let initial_edge = constrained_edges.pop_front().unwrap();
//     set.push(initial_edge);

//     // check every constrained edges until the constrained edge list is empty
//     while constrained_edges.len() != 0 {
//         // get the next edge and add it
//         let current_edge = constrained_edges.pop_front().unwrap();
//         set.push(current_edge);

//         // if we end a cycle,  we stop. We will start a new cycle from were we ended up
//         if current_edge.1 == initial_edge.0 {
//             break;
//         }
//     }

//     set
// }

// ///source:https://stackoverflow.com/questions/471962/how-do-i-efficiently-determine-if-a-polygon-is-convex-non-convex-or-complex
// ///For each consecutive pair of edges of the polygon (each triplet of points), compute the z-component of the cross product of the vectors defined by the edges pointing towards the points in increasing order. Take the cross product of these vectors.
// /// The polygon is convex if the z-components of the cross products are either all positive or all negative. Otherwise the polygon is nonconvex.
// fn is_polygon_convex(polygon: &Vec<(usize, usize)>, vertices: &mut Vec<Vec2>) -> bool {
//     for index in 0..polygon.len() - 1 {
//         let p0 = vertices[polygon[index].0];
//         let p1 = vertices[polygon[index].1];
//         let p2 = vertices[polygon[index + 1].1];

//         let dx1 = p1.x - p0.x;
//         let dy1 = p1.y - p0.y;
//         let dx2 = p2.x - p1.x;
//         let dy2 = p2.y - p1.y;

//         if dx1 * dy2 - dy1 * dx2 == 0. {
//             return false;
//         }
//     }

//     true
// }

fn apply_constrains(
    vertices: &mut Vec<Vec2>,
    constrained_edges: &VecDeque<(usize, usize)>,
    mut triangles: &mut Vec<TriangleData>,
    mut debugger: &mut Vec<Vec<TriangleData>>,
) {
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
        let mut new_edges_created =
            remove_crossing_edges(triangles, vertices, constrained_edge, intersected_edges);

        // Restore Delaunay triangulation
        restore_delaunay_triangulation_constrained(
            &mut triangles,
            vertices,
            constrained_edge,
            &mut new_edges_created,
        );
    }
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
) -> VecDeque<(TriangleId, (VertexId, VertexId))> {
    // List of every triangles already checked
    //TODO: change data structure
    let mut visited_triangles = vec![false; triangles.len()];

    let mut intersected_edges = VecDeque::new();

    // providing a starting triangle to begin the search for the edges which cross the constrained edge
    let mut current_triangle_id = triangle_vertices[constrained_edge.0];
    visited_triangles[current_triangle_id] = true;

    let mut crossing = false;

    // we circle constrained_edge.0
    while !crossing {
        // Vertices of the current triangle
        let triangle_vertex_1 = vertices[triangles[current_triangle_id].v1];
        let triangle_vertex_2 = vertices[triangles[current_triangle_id].v2];
        let triangle_vertex_3 = vertices[triangles[current_triangle_id].v3];

        #[cfg(feature = "debug_constrained_traces")]
        info!("circle around constrained edge {:?}-----------------------------------------------------------", constrained_edge);
        (crossing, _, _) = constrained_edge_intersect_triangle(
            triangles,
            vertices,
            constrained_edge,
            current_triangle_id,
            &visited_triangles,
        );

        if crossing {
            #[cfg(feature = "debug_constrained_traces")]
            info!("end of circle around constrained edge {:?}-------------------------------------------------", constrained_edge);
            break;
        }

        // If the current triangle has already been visited, we go the right
        if visited_triangles[current_triangle_id] {
            #[cfg(feature = "debug_constrained_traces")]
            info!(
                "The triangle {:?} has already been visited",
                current_triangle_id
            );
            if triangle_vertex_1 == vertices[constrained_edge.0] {
                match triangles[current_triangle_id].edge12 {
                    Some(neighbour_triangle_id) => current_triangle_id = neighbour_triangle_id,
                    None => todo!(), // TODO: show error,
                }
            } else if triangle_vertex_2 == vertices[constrained_edge.0] {
                match triangles[current_triangle_id].edge23 {
                    Some(neighbour_triangle_id) => current_triangle_id = neighbour_triangle_id,
                    None => todo!(), // TODO: show error,
                }
            } else if triangle_vertex_3 == vertices[constrained_edge.0] {
                match triangles[current_triangle_id].edge31 {
                    Some(neighbour_triangle_id) => current_triangle_id = neighbour_triangle_id,
                    None => todo!(), // TODO: show error,
                }
            }
        } else {
            //march from one triangle to the next if there are no crossing edges
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
    }

    // Flag to check whenever we are at the end of the constrained edge
    let mut arrived_to_constrained_edge_end = false;

    while !arrived_to_constrained_edge_end {
        #[cfg(feature = "debug_constrained_traces")]
        info!("Check for corssing edges {:?} ------------------------------------------------------------------------", constrained_edge);

        visited_triangles[current_triangle_id] = true;

        // Vertices of the current triangle
        let triangle_vertex_1 = vertices[triangles[current_triangle_id].v1];
        let triangle_vertex_2 = vertices[triangles[current_triangle_id].v2];
        let triangle_vertex_3 = vertices[triangles[current_triangle_id].v3];

        let (_, crossing_edge, crossed_triangle) = constrained_edge_intersect_triangle(
            triangles,
            vertices,
            constrained_edge,
            current_triangle_id,
            &visited_triangles,
        );

        // Check if we are at constrained_edge.1
        if triangle_vertex_1 == vertices[constrained_edge.1]
            || triangle_vertex_2 == vertices[constrained_edge.1]
            || triangle_vertex_3 == vertices[constrained_edge.1]
        {
            #[cfg(feature = "debug_constrained_traces")]
            info!(
                "Triangle {:?} contains constrained edge final vertex {:?}",
                current_triangle_id, constrained_edge
            );
            #[cfg(feature = "debug_constrained_traces")]
            info!("End of checking for corssing edges {:?}----------------------------------------------------------------", constrained_edge);
            arrived_to_constrained_edge_end = true;
        } else {
            current_triangle_id = crossed_triangle.unwrap();
            intersected_edges.push_back((current_triangle_id, crossing_edge.unwrap()));
        }
    }

    #[cfg(feature = "debug_constrained_traces")]
    info!("intersected edges: {:?}", intersected_edges);
    intersected_edges
}

fn remove_crossing_edges(
    triangles: &mut Vec<TriangleData>,
    vertices: &mut Vec<Vec2>,
    constrained_edge: &(usize, usize),
    mut intersected_edges: VecDeque<(usize, (usize, usize))>,
) -> VecDeque<(usize, usize, (usize, usize))> {
    let mut new_edges_created = VecDeque::new();

    #[cfg(feature = "debug_constrained_traces")]
    info!("removing crossing edges for constrained edge {:?}------------------------------------------------------", constrained_edge);

    while intersected_edges.len() != 0 {
        let (current_triangle_id, current_crossing_edge) = intersected_edges.pop_front().unwrap();

        #[cfg(feature = "debug_constrained_traces")]
        info!(
            "crossing edge {:?}, inside triangle {:?}",
            current_crossing_edge, current_triangle_id
        );

        // Construct a quad from the current triangle and its edge by looking at the neighbour triangle of this edge
        let (quad, triangle_id, adjacent_triangle_id) = find_quad_from_triangle_and_crossing_edge(
            current_triangle_id,
            &current_crossing_edge,
            triangles,
        );

        #[cfg(feature = "debug_constrained_traces")]
        info!("quad created: {:?}", quad);

        let quad_diagonal_test = quad_diagonals_intersection_test(&quad, vertices);

        #[cfg(feature = "debug_constrained_traces")]
        info!("quad test: {:?}", quad_diagonal_test);

        // If the quad is not convex or if an edge tip lie on the other edge, we skip this edge
        if quad_diagonal_test == EdgesIntersectionResult::None
            || quad_diagonal_test == EdgesIntersectionResult::OnEdgeTip
            || quad_diagonal_test == EdgesIntersectionResult::SharedEdges
        {
            #[cfg(feature = "debug_constrained_traces")]
            info!("quad not convex");
            intersected_edges.push_back((current_triangle_id, current_crossing_edge));
        }
        // swap the diagonal of this strictly convex quadrilateral if the two diagonals cross normaly
        else {
            #[cfg(feature = "debug_constrained_traces")]
            info!("quad merged");
            let (t3, t4) =
                swap_quad_diagonal(quad.v4, adjacent_triangle_id, triangle_id, triangles); //quad.v2, triangle_id, adjacent_triangle_id, triangles pour 5.0, quad.v4, adjacent_triangle_id,triangle_id, triangles pour 1.3 mais pk ??? => edge orientée !! (1,3), (0,5)

            //  If the new diagonal still intersects the constrained edge, then place it on the list of intersecting edges
            let quad_vertex_2 = vertices[quad.v2];
            let quad_vertex_4 = vertices[quad.v4];
            let constrained_edge_vertex_1 = vertices[constrained_edge.0];
            let constrained_edge_vertex_2 = vertices[constrained_edge.1];

            if egdes_intersect(
                quad_vertex_2,
                quad_vertex_4,
                constrained_edge_vertex_1,
                constrained_edge_vertex_2,
            ) == EdgesIntersectionResult::Crossing
            {
                intersected_edges.push_back((t3.unwrap(), (quad.v2, quad.v4)));
            } else {
                #[cfg(feature = "debug_constrained_traces")]
                info!("edge added: {:?}", (quad.v2, quad.v4));
                new_edges_created.push_back((t3.unwrap(), t4.unwrap(), (quad.v2, quad.v4)));
            }
        }
    }

    #[cfg(feature = "debug_constrained_traces")]
    info!("end of removing crossing edges for constrained edge {:?}------------------------------------------", constrained_edge);
    #[cfg(feature = "debug_constrained_traces")]
    info!("new_edges_created {:?}", new_edges_created);
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
    new_edges_created: &mut VecDeque<(usize, usize, (usize, usize))>,
) {
    while new_edges_created.len() != 0 {
        #[cfg(feature = "debug_constrained_traces")]
        info!("restauring delaunay triangulation constrained--------------------");
        let (current_triangle_id_1, current_triangle_id_2, current_edge) =
            new_edges_created.pop_front().unwrap();

        #[cfg(feature = "debug_constrained_traces")]
        info!(
            "triangle 1: {:?}, triangle 2 {:?}, current edge: {:?}",
            current_triangle_id_1, current_triangle_id_2, current_edge
        );
        // If the edge is equal to the constrained edge, then skip
        if (current_edge.0 == constrained_edge.0 && current_edge.1 == constrained_edge.1)
            || (current_edge.0 == constrained_edge.1 && current_edge.1 == constrained_edge.0)
        {
            #[cfg(feature = "debug_constrained_traces")]
            info!("edge skipped");
            continue;
        } else {
            let (quad_diag_swapped, t3, t4, edge_1, edge_2) = swap_quad_diagonal_to_edge(
                current_edge.0,
                current_triangle_id_1,
                current_triangle_id_2,
                triangles,
                vertices,
            );

            if quad_diag_swapped {
                new_edges_created.push_back((t3.unwrap(), t4.unwrap(), (edge_1, edge_2)));
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
///               q2
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

    #[cfg(feature = "debug_constrained_traces")]
    info!("v1 :{:?},v2 :{:?},v3 :{:?}", vertex_1, vertex_2, vertex_3);

    // Quad coords
    let q1 = crossing_edge_v1;
    let q3 = crossing_edge_v2;

    // Shared triangle
    // Crossing edge starts at v1
    let (neighbour_triangle_id, q2) = if crossing_edge_v1 == vertex_1 {
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
    constrained_edge: &(usize, usize),
    current_triangle_id: TriangleId,
    visited_triangles: &Vec<bool>,
) -> (bool, Option<(usize, usize)>, Option<usize>) {
    // Set the constrained edge vertices and the triangle vertices
    let constrained_edge_start = vertices[constrained_edge.0];
    let constrained_edge_end = vertices[constrained_edge.1];

    let triangle_vertex_1 = vertices[triangles[current_triangle_id].v1];
    let triangle_vertex_2 = vertices[triangles[current_triangle_id].v2];
    let triangle_vertex_3 = vertices[triangles[current_triangle_id].v3];

    #[cfg(feature = "debug_constrained_traces")]
    info!(
        "check for intersection inside triangle {}",
        current_triangle_id
    );
    #[cfg(feature = "debug_constrained_traces")]
    info!("check edge 12 inside triangle {}", current_triangle_id);
    // Check if constrained edge cross v1-v2
    if egdes_intersect(
        constrained_edge_start,
        constrained_edge_end,
        triangle_vertex_1,
        triangle_vertex_2,
    ) == EdgesIntersectionResult::Crossing
    {
        #[cfg(feature = "debug_constrained_traces")]
        info!("edge 12 is interected");
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
    #[cfg(feature = "debug_constrained_traces")]
    info!("check edge 23 inside triangle {}", current_triangle_id);
    // Check if constrained edge cross v2-v3
    if egdes_intersect(
        constrained_edge_start,
        constrained_edge_end,
        triangle_vertex_2,
        triangle_vertex_3,
    ) == EdgesIntersectionResult::Crossing
    {
        #[cfg(feature = "debug_constrained_traces")]
        info!("edge 23 is interected");
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
                    #[cfg(feature = "debug_constrained_traces")]
                    info!("SHOULD NOT PRINT")
                }
            }
            None => return (false, None, None), // TODO: not possible error,
        }
    }

    #[cfg(feature = "debug_constrained_traces")]
    info!("check edge 31 inside triangle {}", current_triangle_id);

    // Check if constrained edge cross v3-v1
    if egdes_intersect(
        constrained_edge_start,
        constrained_edge_end,
        triangle_vertex_3,
        triangle_vertex_1,
    ) == EdgesIntersectionResult::Crossing
    {
        #[cfg(feature = "debug_constrained_traces")]
        info!("edge 31 is interected");
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

#[derive(PartialEq, Eq, Debug)]
enum EdgesIntersectionResult {
    None,
    Crossing,
    OnEdgeTip,
    SharedEdges,
}

// source: https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
fn egdes_intersect(p1: Vec2, q1: Vec2, p2: Vec2, q2: Vec2) -> EdgesIntersectionResult {
    if p1 == p2 || p1 == q2 || q1 == p2 || q1 == q2 {
        return EdgesIntersectionResult::SharedEdges;
    }

    let orientation_1 = orientation(p1, q1, p2);
    let orientation_2 = orientation(p1, q1, q2);
    let orientation_3 = orientation(p2, q2, p1);
    let orientation_4 = orientation(p2, q2, q1);

    // Special Cases
    // p1, q1 and p2 are collinear and p2 lies on segment p1q1
    if orientation_1 == 0 && on_segment(p1, p2, q1) {
        return EdgesIntersectionResult::OnEdgeTip;
    }

    // p1, q1 and q2 are collinear and q2 lies on segment p1q1
    if orientation_2 == 0 && on_segment(p1, q2, q1) {
        return EdgesIntersectionResult::OnEdgeTip;
    }

    // p2, q2 and p1 are collinear and p1 lies on segment p2q2
    if orientation_3 == 0 && on_segment(p2, p1, q2) {
        return EdgesIntersectionResult::OnEdgeTip;
    }

    // p2, q2 and q1 are collinear and q1 lies on segment p2q2
    if orientation_4 == 0 && on_segment(p2, q1, q2) {
        return EdgesIntersectionResult::OnEdgeTip;
    }

    // General case
    if orientation_1 != orientation_2 && orientation_3 != orientation_4 {
        return EdgesIntersectionResult::Crossing;
    }

    EdgesIntersectionResult::None // Doesn't fall in any of the above cases
}

// Given three collinear points p, q, r, the function checks if
// point q lies on line segment 'pr'
fn on_segment(p: Vec2, q: Vec2, r: Vec2) -> bool {
    q.x <= p.x.max(r.x) && q.x >= p.x.min(r.x) && q.y <= p.y.max(r.y) && q.y >= p.y.min(r.y)
}

// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are collinear
// 1 --> Clockwise
// 2 --> Counterclockwise
fn orientation(p: Vec2, q: Vec2, r: Vec2) -> usize {
    let val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);

    if val == 0. {
        0
    } else if val > 0. {
        1
    } else {
        2
    }
}
