use bevy::math::{Vec2, Vec3A};

use super::{Quad, TriangleData, TriangleId, VertexId, TRIANGLE_EDGES};

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

pub(crate) fn search_enclosing_triangle(
    vertex: Vec2,
    from: Option<TriangleId>,
    triangles: &Vec<TriangleData>,
    vertices: &Vec<Vec2>,
) -> Option<TriangleId> {
    let mut triangle_id = from;

    for _ in 0..triangles.len() {
        let Some(current_triangle_id) = triangle_id else {
            break;
        };
        let triangle = &triangles[current_triangle_id];

        // Check if the point is inside the triangle, if not check the neighbours
        let mut inside_triangle = true;
        for (triangle_edge_index, triangle_edge) in TRIANGLE_EDGES.iter().enumerate() {
            let edge = triangle.edge_from_triangle_edge(triangle_edge);
            if !is_vertex_on_right_side_of_edge(vertices[edge.0], vertices[edge.1], vertex) {
                triangle_id = triangle.neighbors[triangle_edge_index];
                inside_triangle = false;
                break;
            }
        }
        if inside_triangle {
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
        verts: [vertices.len(), vertices.len() + 1, vertices.len() + 2],
        neighbors: [None, None, None],
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
    debugger: &mut Vec<Vec<TriangleData>>,
) -> Vec<VertexId> {
    // TODO Clean: Size approx
    let mut indices = Vec::with_capacity(3 * triangles.len());
    let mut filtered_debug_triangles = Vec::new();

    for triangle in triangles.iter() {
        let mut filtered = false;
        for &vert in triangle.verts.iter() {
            for &container_vert in container_triangle.verts.iter() {
                if vert == container_vert {
                    filtered = true;
                }
            }
        }
        if !filtered {
            indices.push(triangle.v1());
            indices.push(triangle.v2());
            indices.push(triangle.v3());
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
    // Re-use the existing triangle id for the first triangle
    let triangle_1_id = triangle_id;
    // Create two new triangles for the other two
    let triangle_2_id = triangles.len();
    let triangle_3_id = triangles.len() + 1;

    triangles.push(TriangleData {
        verts: [
            vertex_id,
            triangles[triangle_id].v2(),
            triangles[triangle_id].v3(),
        ],
        neighbors: [
            Some(triangle_3_id),
            triangles[triangle_id].neighbor23(),
            Some(triangle_1_id),
        ],
    });
    triangles.push(TriangleData {
        verts: [
            vertex_id,
            triangles[triangle_id].v1(),
            triangles[triangle_id].v2(),
        ],
        neighbors: [
            Some(triangle_1_id),
            triangles[triangle_id].neighbor12(),
            Some(triangle_2_id),
        ],
    });

    // Update triangle indexes
    update_triangle_neighbour(
        triangles[triangle_id].neighbor12(),
        Some(triangle_id),
        Some(triangle_3_id),
        triangles,
    );
    update_triangle_neighbour(
        triangles[triangle_id].neighbor23(),
        Some(triangle_id),
        Some(triangle_2_id),
        triangles,
    );

    // Replace id_triangle with id_triangle_1
    triangles[triangle_1_id].verts = [
        vertex_id,
        triangles[triangle_id].v3(),
        triangles[triangle_id].v1(),
    ];

    triangles[triangle_1_id].neighbors = [
        Some(triangle_2_id),
        triangles[triangle_id].neighbor31(),
        Some(triangle_3_id),
    ];

    restore_delaunay_triangulation(
        triangles,
        vertex_id,
        [triangle_1_id, triangle_2_id, triangle_3_id],
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
            for neighbor in triangle.neighbors.iter_mut() {
                if *neighbor == old_neighbour_id {
                    *neighbor = new_neighbour_id;
                    break;
                }
            }
        }
        None => (),
    };
}

fn restore_delaunay_triangulation(
    triangles: &mut Vec<TriangleData>,
    vertex_id: VertexId,
    triangles_to_check: [TriangleId; 3],
    vertices: &Vec<Vec2>,
) {
    let mut stack = Vec::<(TriangleId, Option<TriangleId>)>::new();

    for &triangle_id in &triangles_to_check {
        stack.push((triangle_id, triangles[triangle_id].neighbor23()));
    }

    while let Some((triangle_id, adjacent_triangle_id)) = stack.pop() {
        let Some(adjacent_triangle_id) = adjacent_triangle_id else {
            continue;
        };

        let (quad_diag_swapped, t3, t4) = check_and_swap_quad_diagonal(
            triangles,
            vertices,
            vertex_id,
            triangle_id,
            adjacent_triangle_id,
        );
        if quad_diag_swapped == QuadSwapResult::Swapped {
            // Place any triangles which are now opposite P on the stack.
            stack.push((triangle_id, t3));
            stack.push((adjacent_triangle_id, t4));
        }
    }
}

#[derive(PartialEq, Eq)]
pub enum QuadSwapResult {
    NotSwapped,
    Swapped,
}

pub fn check_and_swap_quad_diagonal(
    triangles: &mut Vec<TriangleData>,
    vertices: &Vec<Vec2>,
    vertex_id: VertexId,
    triangle_id: TriangleId,
    adjacent_triangle_id: TriangleId,
) -> (QuadSwapResult, Option<TriangleId>, Option<TriangleId>) {
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
    let swapped_quad_diagonal = if is_vertex_in_triangle_circumcircle(
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

        triangles[triangle_id].set_neighbor23(triangle_3_id);
        triangles[triangle_id].set_neighbor31(Some(adjacent_triangle_id));

        QuadSwapResult::Swapped
    } else {
        QuadSwapResult::NotSwapped
    };
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
