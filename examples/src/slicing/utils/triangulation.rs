use bevy::{
    log::info,
    math::{Vec2, Vec3A, Vec3Swizzles},
};

/// plane_normal must be normalized
/// vertices must all belong to a 3d plane
pub fn triangulate_3d_planar_vertices(
    vertices: &Vec<[f32; 3]>,
    plane_normal: Vec3A,
) -> (Vec<VertexId>, Vec<Vec<TriangleData>>) {
    // TODO See what we need for input data format of `triangulate`
    let mut vertices_data = Vec::with_capacity(vertices.len());
    for v in vertices {
        vertices_data.push(Vec3A::from_array(*v));
    }

    let mut planar_vertices =
        transform_to_2d_planar_coordinate_system(&mut vertices_data, plane_normal);

    info!(
        "transform_to_2d_planar_coordinate_system {:?}",
        planar_vertices
    );

    // Delaunay triangulation
    triangulate_2d_vertices(&mut planar_vertices)
}

/// Transforms 3d coordinates of all vertices into 2d coordinates on a plane defined by the given normal and vertices.
/// - Input vertices need to all belong to the same 3d plan
/// - There must be at least two vertices
pub(crate) fn transform_to_2d_planar_coordinate_system(
    vertices: &mut Vec<Vec3A>,
    plane_normal: Vec3A,
) -> Vec<Vec2> {
    // Create a base, using the first two vertices as the first base vector and plane_normal as the second
    let basis_1 = (vertices[0] - vertices[1]).normalize();
    // basis_3 is already normalized since basis_1 and plane_normal are normalized and orthogonal
    let basis_3 = basis_1.cross(plane_normal);

    // Project every vertices into the base B
    let mut vertices_2d = Vec::with_capacity(vertices.len());
    for vertex in vertices {
        vertices_2d.push(Vec2::new(vertex.dot(basis_1), vertex.dot(basis_3)));

        //vertices_2d.push(vertex.xy());
    }
    vertices_2d
}

/// This scaling ensures that all of the coordinates are between 0 and 1 but does not modify the relative positions of the points in the x-y plane.
/// The use of normalized coordinates, although not essential, reduces the effects of  roundoff error and is also convenient from a computational point of view.
pub(crate) fn normalize_vertices_coordinates(vertices: &mut Vec<Vec2>) {
    let (mut x_min, mut y_min, mut x_max, mut y_max) = (f32::MAX, f32::MAX, f32::MIN, f32::MIN);

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

    for vertex in vertices.iter_mut() {
        vertex.x = (vertex.x - x_min) / scale_factor;
        vertex.y = (vertex.y - y_min) / scale_factor;
    }
}

pub(crate) type VertexId = usize;
pub(crate) type TriangleId = usize;

#[derive(Debug, Clone)]
pub struct TriangleData {
    // Vertices ids
    pub v1: VertexId,
    pub v2: VertexId,
    pub v3: VertexId,
    // Neighbours
    pub edge12: Option<TriangleId>,
    pub edge23: Option<TriangleId>,
    pub edge31: Option<TriangleId>,
    //pub edges: [TriangleId;3],
}

const EDGE_12: usize = 0; 
const EDGE_23: usize = 1; 
const EDGE_31: usize = 2; 
const EDGES: [usize;3] = [EDGE_12, EDGE_23, EDGE_31];

pub(crate) fn search_enclosing_triangle(
    vertex: Vec2,
    from: Option<TriangleId>,
    triangles: &Vec<TriangleData>,
    vertices: &Vec<Vec2>,
) -> Option<TriangleId> {
    let mut triangle_id = from;
    #[cfg(feature = "debug_traces")]
    info!("search_enclosing_triangle----------------------------------------------------------------------------");
    for iter in 0..triangles.len() {
        #[cfg(feature = "debug_traces")]
        info!("inter{}", iter);
        let Some(current_triangle_id) = triangle_id else {
            break;
        };

        // Vertices of the actual triangle
        let v1 = vertices[triangles[current_triangle_id].v1];
        let v2 = vertices[triangles[current_triangle_id].v2];
        let v3 = vertices[triangles[current_triangle_id].v3];

        #[cfg(feature = "debug_traces")]
        info!("point to be placed {}", vertex);

        // Check if the point is inside the triangle, if not check the neighbours
        if !is_vertex_on_right_side_of_edge(v1, v2, vertex) {
            #[cfg(feature = "debug_traces")]
            info!(
                "triangle{},v1{}, v2{}, vertex{} ",
                current_triangle_id, v1, v2, vertex
            );
            triangle_id = triangles[current_triangle_id].edge12;
        } else if !is_vertex_on_right_side_of_edge(v2, v3, vertex) {
            #[cfg(feature = "debug_traces")]
            info!(
                "triangle{},v1{}, v2{}, vertex{} ",
                current_triangle_id, v1, v2, vertex
            );
            triangle_id = triangles[current_triangle_id].edge23;
        } else if !is_vertex_on_right_side_of_edge(v3, v1, vertex) {
            #[cfg(feature = "debug_traces")]
            info!(
                "triangle{},v1{}, v2{}, vertex{} ",
                current_triangle_id, v1, v2, vertex
            );
            triangle_id = triangles[current_triangle_id].edge31;
        } else {
            #[cfg(feature = "debug_traces")]
            info!("current triangle{}", current_triangle_id);
            return Some(current_triangle_id);
        }
    }
    None
}

fn triangulate_2d_vertices(vertices: &mut Vec<Vec2>) -> (Vec<VertexId>, Vec<Vec<TriangleData>>) {
    let (triangles, container_triangle, mut debugger) = wrap_and_triangulate_2d_vertices(vertices);

    let indices = remove_wrapping(&triangles, &container_triangle, &mut debugger);

    (indices, debugger)
}

pub(crate) fn wrap_and_triangulate_2d_vertices(
    vertices: &mut Vec<Vec2>,
) -> (Vec<TriangleData>, TriangleData, Vec<Vec<TriangleData>>) {
    //Debug
    let mut debugger: Vec<Vec<TriangleData>> = Vec::new();

    // Uniformly scale the coordinates of the points so that they all lie between 0 and 1.
    normalize_vertices_coordinates(vertices);

    // Sort points into bins. Cover the region to be triangulated by a rectangular grid so that each bin contains roughly N^(1/2) points.
    // Label the bins so that consecutive bins are adjacent to one another, and then allocate each point to its appropriate bin.
    // Sort the list of points in ascending sequence of their bin numbers so that consecutive points are grouped together in the x-y plane.
    let partitioned_vertices = VertexBinSort::sort(&vertices, 0.5);

    // Select three dummy points to form a supertriangle that completely encompasses all of the points to be triangulated.
    // This supertriangle initially defines a Delaunay triangulation which is comprised of a single triangle.
    // Its vertices are defined in terms of normalized coordinates and are usually located at a considerable distance from the window which encloses the set of points.
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

    debugger.push(triangles.clone());

    // Loop over all the input vertices
    for sorted_vertex in partitioned_vertices.iter() {
        // Find an existing triangle which encloses P
        match search_enclosing_triangle(sorted_vertex.vertex, triangle_id, &triangles, &vertices) {
            Some(enclosing_triangle_id) => {
                // Delete this triangle and form three new triangles by connecting P to each of its vertices.
                split_triangle_in_three_at_vertex(
                    enclosing_triangle_id,
                    sorted_vertex.original_id,
                    &mut triangles,
                    &vertices,
                );
                triangle_id = Some(triangles.len() - 1);
            }
            None => (),
        }
        debugger.push(triangles.clone());
    }

    (triangles, container_triangle, debugger)
}

pub(crate) fn remove_wrapping(
    triangles: &Vec<TriangleData>,
    container_triangle: &TriangleData,
    mut debugger: &mut Vec<Vec<TriangleData>>,
) -> Vec<VertexId> {
    // TODO Clean: Size approx
    let mut indices = Vec::with_capacity(3 * triangles.len());
    let mut filtered_debug_triangles = Vec::new();
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
            filtered_debug_triangles.push(triangle.clone());
        }
    }
    debugger.push(filtered_debug_triangles);

    indices
}

pub(crate) struct VertexBinSort {
    bins_per_row: usize,
    bins_count: usize,
}

#[derive(Debug, Default, Copy, Clone)]
pub(crate) struct SortedVertex {
    pub(crate) original_id: VertexId,
    pub(crate) vertex: Vec2,
}

impl VertexBinSort {
    // Each bin will contain roughly vertices.len()^(vertex_density_power) vertices
    pub fn sort(vertices: &Vec<Vec2>, vertex_density_power: f32) -> Vec<SortedVertex> {
        let bins_per_row = (vertices.len() as f32)
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

        let mut sorted = vec![SortedVertex::default(); vertices.len()];
        for (vertex_id, vertex) in vertices.iter().enumerate() {
            let bin_index = vertices_bin_indexes[vertex_id];
            bins_counters[bin_index] -= 1;
            sorted[bins_counters[bin_index]] = SortedVertex {
                original_id: vertex_id,
                vertex: *vertex,
            };
        }
        sorted
    }

    fn bin_index_from_vertex(&self, vertex: Vec2) -> usize {
        // Compute a bin index from a vertox position which is in [0.,1]
        let bin_x = (0.99 * self.bins_per_row as f32 * vertex.x) as usize;
        let bin_y: usize = (0.99 * self.bins_per_row as f32 * vertex.y) as usize;
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

fn is_vertex_on_right_side_of_edge(v1: Vec2, v2: Vec2, p: Vec2) -> bool {
    // Cross product of vectors v1-v2 and v1-p
    //((v2.x - v1.x) * (p.y - v1.y) - (v2.y - v1.y) * (p.x - v1.x)) <= 0.
    ((p.x - v1.x) * (v2.y - v1.y) - (p.y - v1.y) * (v2.x - v1.x)) >= 0.
}

pub(crate) fn split_triangle_in_three_at_vertex(
    triangle_id: TriangleId,
    vertex_id: VertexId,
    triangles: &mut Vec<TriangleData>,
    vertices: &Vec<Vec2>,
) {
    #[cfg(feature = "debug_traces")]
    info!(
        "split_triangle_in_three_at_vertex, vertex {}",
        vertices[vertex_id]
    );

    // Re-use the existing triangle id for the first triangle
    let triangle_1_id = triangle_id;
    // Create two new triangles for the other two
    let triangle_2_id = triangles.len();
    let triangle_3_id = triangles.len() + 1;

    triangles.push(TriangleData {
        v1: vertex_id,
        v2: triangles[triangle_id].v2,
        v3: triangles[triangle_id].v3,
        edge12: Some(triangle_3_id),
        edge23: triangles[triangle_id].edge23,
        edge31: Some(triangle_1_id),
    });
    triangles.push(TriangleData {
        v1: vertex_id,
        v2: triangles[triangle_id].v1,
        v3: triangles[triangle_id].v2,
        edge12: Some(triangle_1_id),
        edge23: triangles[triangle_id].edge12,
        edge31: Some(triangle_2_id),
    });

    // Update triangle indexes
    update_triangle_neighbour(
        triangles[triangle_id].edge12,
        Some(triangle_id),
        Some(triangle_3_id),
        triangles,
    );
    update_triangle_neighbour(
        triangles[triangle_id].edge23,
        Some(triangle_id),
        Some(triangle_2_id),
        triangles,
    );

    // Replace id_triangle with id_triangle_1
    triangles[triangle_1_id].v2 = triangles[triangle_id].v3;
    triangles[triangle_1_id].v3 = triangles[triangle_id].v1;
    triangles[triangle_1_id].v1 = vertex_id;
    triangles[triangle_1_id].edge23 = triangles[triangle_id].edge31;
    triangles[triangle_1_id].edge12 = Some(triangle_2_id);
    triangles[triangle_1_id].edge31 = Some(triangle_3_id);

    restore_delaunay_triangulation(
        vertex_id,
        triangle_1_id,
        triangle_2_id,
        triangle_3_id,
        triangles,
        vertices,
    );
}

pub(crate) fn update_triangle_neighbour(
    triangle_id: Option<TriangleId>,
    old_neighbour_id: Option<TriangleId>,
    new_neighbour_id: Option<TriangleId>,
    triangles: &mut Vec<TriangleData>,
) {
    match triangle_id {
        Some(triangle_id) => {
            let triangle = &mut triangles[triangle_id];
            if triangle.edge12 == old_neighbour_id {
                triangle.edge12 = new_neighbour_id;
            } else if triangle.edge23 == old_neighbour_id {
                triangle.edge23 = new_neighbour_id;
            } else if triangle.edge31 == old_neighbour_id {
                triangle.edge31 = new_neighbour_id;
            }
        }
        None => (),
    };
}

fn restore_delaunay_triangulation(
    vertex_id: VertexId,
    triangle_id_1: TriangleId,
    triangle_id_2: TriangleId,
    triangle_id_3: TriangleId,
    triangles: &mut Vec<TriangleData>,
    vertices: &Vec<Vec2>,
) {
    let mut stack = Vec::<(TriangleId, Option<TriangleId>)>::new();

    stack.push((triangle_id_1, triangles[triangle_id_1].edge23));
    stack.push((triangle_id_2, triangles[triangle_id_2].edge23));
    stack.push((triangle_id_3, triangles[triangle_id_3].edge23));

    while let Some((triangle_id, adjacent_triangle_id)) = stack.pop() {
        let Some(adjacent_triangle_id) = adjacent_triangle_id else {
            continue;
        };

        let (quad_diag_swapped, t3, t4) = check_and_swap_quad_diagonal(
            vertex_id,
            triangle_id,
            adjacent_triangle_id,
            triangles,
            vertices,
        );
        if quad_diag_swapped {
            // Place any triangles which are now opposite P on the stack.
            stack.push((triangle_id, t3));
            stack.push((adjacent_triangle_id, t4));
        }
    }
}

#[derive(Debug)]
pub(crate) struct Quad {
    pub(crate) v1: VertexId,
    pub(crate) v2: VertexId,
    pub(crate) v3: VertexId,
    pub(crate) v4: VertexId,
}

pub(crate) fn check_and_swap_quad_diagonal(
    vertex_id: VertexId,
    triangle_id: TriangleId,
    adjacent_triangle_id: TriangleId,
    triangles: &mut Vec<TriangleData>,
    vertices: &Vec<Vec2>,
) -> (bool, Option<TriangleId>, Option<TriangleId>) {
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
    if is_vertex_in_triangle_circumcircle(
        vertices[quad.v1],
        vertices[quad.v2],
        vertices[quad.v3],
        vertices[vertex_id],
    ) {
        // The triangle containing P as a vertex and the unstacked triangle form a convex quadrilateral whose diagonal is drawn in the wrong direction.
        // Swap this diagonal so that two old triangles are replaced by two new triangles and the structure of the Delaunay triangulation is locally restored.
        update_triangle_neighbour(
            triangle_3_id,
            Some(adjacent_triangle_id),
            Some(triangle_id),
            triangles,
        );
        update_triangle_neighbour(
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
    (swapped_quad_diagonal, triangle_3_id, triangle_4_id)
}

/// See ref: Cline and Renka [5]
/// The algo is sensitive to the order of the given vertices. In order:
///
/// n3 --------------------- n2
/// |                       /|
/// |                     /  |
/// |                   /    |
/// |                 /      |
/// |               /        |
/// |             /          |
/// |           /            |
/// |         /              |
/// |       /                |
/// |     /                  |
/// |   /                    |
/// | /                      |
/// n1 ---------------------- p
///
/// where n1, n2 and n3 are the given vertices of the given triangle and p the vertex to be check
///
pub(crate) fn is_vertex_in_triangle_circumcircle(
    triangle_v1: Vec2,
    triangle_v2: Vec2,
    triangle_v3: Vec2,
    vertex: Vec2,
) -> bool {
    let x13 = triangle_v1.x - triangle_v3.x;
    let x23 = triangle_v2.x - triangle_v3.x;
    let y13 = triangle_v1.y - triangle_v3.y;
    let y23 = triangle_v2.y - triangle_v3.y;
    let x14 = triangle_v1.x - vertex.x;
    let x24 = triangle_v2.x - vertex.x;
    let y14 = triangle_v1.y - vertex.y;
    let y24 = triangle_v2.y - vertex.y;

    let cos_a = x13 * x23 + y13 * y23;
    let cos_b = x24 * x14 + y24 * y14;

    if cos_a >= 0. && cos_b >= 0. {
        false
    } else if cos_a < 0. && cos_b < 0. {
        true
    } else {
        let sin_a = x13 * y23 - x23 * y13;
        let sin_b = x24 * y14 - x14 * y24;
        let sin_ab = sin_a * cos_b + sin_b * cos_a;
        sin_ab < 0.
    }
}
