use bevy::math::Vec3A;

use super::utils::{cross_from_vec3a, dot_from_array, normalize_vertexe, substract_two_vertices};

pub struct Triangulation {}

struct TriangleData {
    vert_indices: [usize; 3],
    neighbours: [usize; 3],
}

// impl Triangulation {
pub fn triangulate(vertices_to_be_triangulate: &Vec<[f32; 3]>, normal_vec: Vec3A) -> Vec<[f32; 2]> {
    // Vec that contains the id of every triangles and there vertices + id of the three neighbour triangles:
    let mut triangle_data = Vec::<Vec<usize>>::new();

    //Transform every 3D coords of all vertices into the same 2D plane
    let mut projected_vertices =
        project_vertices_coords_in_2d_plane(vertices_to_be_triangulate, normal_vec);

    // add big triangle into the triangle vertices list
    add_big_triangle(&mut projected_vertices, &mut triangle_data);

    // make sure that every vertices of all triangles are into the big triangle, by normalize coords (not big one !)
    let mut normalize_vertices = get_projected_vertices_in_big_triangle(&mut projected_vertices);

    // Delaunay triangulation
    triangulate_projected_vertices(&mut normalize_vertices, &mut triangle_data);

    // Remove extra triangles (by selecting them from big triangle vertices)
    let triangulate_vertices = remove_extra_triangles(&mut normalize_vertices, &mut triangle_data);

    triangulate_vertices
}

fn project_vertices_coords_in_2d_plane(
    vertices_to_be_triangulate: &Vec<[f32; 3]>,
    normal_vec: Vec3A,
) -> Vec<[f32; 2]> {
    // First need to create the base B(e1,e2)
    let base_vec_1 = normalize_vertexe(substract_two_vertices(
        vertices_to_be_triangulate[0],
        vertices_to_be_triangulate[1],
    ));

    let base_vec_2 = normalize_vertexe(cross_from_vec3a(base_vec_1, normal_vec.normalize()));

    // Project every vertices into the base B
    let mut projected_vertices = Vec::new();

    for vertexe in vertices_to_be_triangulate {
        let x = dot_from_array(vertexe, &base_vec_1);
        let y = dot_from_array(vertexe, &base_vec_2);
        let projected_vertexe = [x, y];
        projected_vertices.push(projected_vertexe);
    }

    projected_vertices
}

fn add_big_triangle(projected_vertices: &mut Vec<[f32; 2]>, triangle_data: &mut Vec<Vec<usize>>) {
    let triangles = projected_vertices.len();

    projected_vertices.push([-100., -100.]);
    projected_vertices.push([100., 0.]);
    projected_vertices.push([100., -100.]);

    let mut big_triangle = Vec::new();
    big_triangle.push(triangles);
    big_triangle.push(triangles + 1);
    big_triangle.push(triangles + 2);

    big_triangle.push(100000);
    big_triangle.push(100000);
    big_triangle.push(100000);

    triangle_data.push(big_triangle);
}

fn get_projected_vertices_in_big_triangle(projected_vertices: &mut Vec<[f32; 2]>) -> Vec<[f32; 2]> {
    // the max/min of the coords of vertices
    let mut x_min = 0.;
    let mut y_min = 0.;
    let mut x_max = 0.;
    let mut y_max = 0.;

    for vertexe in projected_vertices.iter() {
        let x = vertexe[0];
        let y = vertexe[1];
        if x < x_min {
            x_min = x;
        }
        if x > x_max {
            x_max = x;
        }
        if y < y_min {
            y_min = y;
        }
        if y > y_max {
            y_max = y;
        }
    }

    // Scale number for normalization
    let scale_factor = (x_max - x_min).max(y_max - y_min);

    let mut normalize_vertices = Vec::new();

    // Create the new normalized vertices
    for vertexe in projected_vertices.iter() {
        let x = vertexe[0];
        let y = vertexe[1];

        let normalize_x = (x - x_min) / scale_factor;
        let normalize_y = (y - y_min) / scale_factor;

        normalize_vertices.push([normalize_x, normalize_y])
    }

    normalize_vertices
}

fn triangulate_projected_vertices(
    normalize_vertices: &mut Vec<[f32; 2]>,
    triangle_data: &mut Vec<Vec<usize>>,
    // triangle_data: &mut Vec<TriangleData>,
) {
    // Total number of triangles
    let mut number_of_triangles = normalize_vertices.len() / 3;

    // Id of the triangle we are looking at
    let mut id_triangle = 0;

    // The id of the last triangle created
    let id_last_triangle_formed = 0;

    let mut bin_vertices = sort_vertices_into_bins(normalize_vertices);

    for (id_vertex, vertexe) in bin_vertices.iter().enumerate() {
        let mut is_vertexe_inserted = false;
        let mut iteration_counter = 0;

        while !is_vertexe_inserted {
            iteration_counter += 1;
            // We can't go past the max number of triangles
            if iteration_counter > number_of_triangles || id_triangle == 100000 {
                break;
            }

            // vertices of the actual triangle
            let v1 = bin_vertices[triangle_data[id_triangle][0]];
            let v2 = bin_vertices[triangle_data[id_triangle][1]];
            let v3 = bin_vertices[triangle_data[id_triangle][2]];

            // Check if the point is inside the triangle, if not check the neighbours
            if !is_vertexe_on_right_side_of_an_edge(v1, v2, *vertexe) {
                id_triangle = triangle_data[id_triangle][3];
            }
            if !is_vertexe_on_right_side_of_an_edge(v2, v3, *vertexe) {
                id_triangle = triangle_data[id_triangle][4];
            }
            if !is_vertexe_on_right_side_of_an_edge(v3, v1, *vertexe) {
                id_triangle = triangle_data[id_triangle][5];
            } else {
                add_vertexe_in_triangle(
                    id_triangle,
                    id_last_triangle_formed,
                    triangle_data,
                    id_vertex,
                    &bin_vertices,
                );
                number_of_triangles += 2;
                id_triangle = id_last_triangle_formed;
                is_vertexe_inserted = true;
            }
        }
    }
}

struct VerticesBins {
    bins_per_row: usize,
    bins_count: usize,
}

impl VerticesBins {
    // Each bin will contain roughly N^(vertex_density_power) points
    pub fn sort(vertices: &[[f32; 2]], vertex_density_power: f32) {
        let bins_per_row = (vertices.len() as f32)
            .powf(vertex_density_power / 2.)
            .round() as usize;
        let bins = Self {
            bins_per_row,
            bins_count: bins_per_row * bins_per_row,
        };

        // Indexes of the bin corresponding to each vertex
        let mut vertices_bin_indexes = Vec::with_capacity(vertices.len());

        for vertex in vertices {
            vertices_bin_indexes.push(bins.bin_index_from_vertex(vertex));
        }
    }

    fn bin_index_from_vertex(&self, vertex: &[f32; 2]) -> usize {
        let bin_x = (0.99 * self.bins_per_row as f32 * vertex[0]) as usize;
        let bin_y: usize = (0.99 * self.bins_per_row as f32 * vertex[1]) as usize;
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

fn sort_vertices_into_bins(normalize_vertices: &mut Vec<[f32; 2]>) -> Vec<[f32; 2]> {
    // Ignore the vertices of the initial wrapping triangle
    let vertex_count = normalize_vertices.len() - 3;

    // Cover the region to be triangulated by a rectangular grid so that each bin contains roughly N^(1/2) points.
    let partitioned_vertices = VerticesBins::sort(&normalize_vertices[0..vertex_count], 0.5);

    // let bins_per_row: usize = (vertex_count as f32).powf(0.25).round() as usize;

    todo!()
}

// fn bin_index_from_vertex(vertex: &[f32; 2], bins_per_row: usize) -> usize {
//     let bin_x = (0.99 * bins_per_row as f32 * vertex[0]).round() as usize;
//     let bin_y: usize = (0.99 * bins_per_row as f32 * vertex[1]).round() as usize;
//     bin_index_from_bin_position(bin_x, bin_y, bins_per_row)
// }

// // Label the bins so that bins with consecutive indexes are spatially adjacent to one another
// fn bin_index_from_bin_position(x: usize, y: usize, bins_per_row: usize) -> usize {
//     if y % 2 == 0 {
//         (y * bins_per_row) + x
//     } else {
//         (y + 1) * bins_per_row - x - 1
//     }
// }

// fn sort_bin_vertices ()

fn remove_extra_triangles(
    delaunay_triangulate_vertices: &mut Vec<[f32; 2]>,
    triangle_data: &mut Vec<Vec<usize>>,
) -> Vec<[f32; 2]> {
    let triangles_to_be_removed =
        remove_extra_triangles_vertices(delaunay_triangulate_vertices, triangle_data);

    let mut triangles = Vec::new();

    for id_triangle in (0..delaunay_triangulate_vertices.len()).step_by(3) {
        if !triangles_to_be_removed[id_triangle] {
            triangles.push(delaunay_triangulate_vertices[triangle_data[id_triangle][0]]);
            triangles.push(delaunay_triangulate_vertices[triangle_data[id_triangle][1]]);
            triangles.push(delaunay_triangulate_vertices[triangle_data[id_triangle][2]]);
        }
    }

    triangles
}

fn remove_extra_triangles_vertices(
    delaunay_triangulate_vertices: &Vec<[f32; 2]>,
    triangle_data: &Vec<Vec<usize>>,
) -> Vec<bool> {
    let big_id_triangle_vertexe = delaunay_triangulate_vertices.len() - 3;
    let mut triangles_to_be_removed = Vec::new();
    for id_triangle in (0..delaunay_triangulate_vertices.len()).step_by(3) {
        if contains_big_triangle_vertexe(&id_triangle, &big_id_triangle_vertexe, triangle_data)
            || contains_big_triangle_vertexe(
                &id_triangle,
                &(big_id_triangle_vertexe - 1),
                triangle_data,
            )
            || contains_big_triangle_vertexe(
                &id_triangle,
                &(big_id_triangle_vertexe - 2),
                triangle_data,
            )
        {
            triangles_to_be_removed[id_triangle] = true;
        } else {
            triangles_to_be_removed[id_triangle] = false;
        }
    }

    triangles_to_be_removed
}

fn contains_big_triangle_vertexe(
    id_triangle: &usize,
    big_id_triangle_vertexe: &usize,
    triangle_data: &Vec<Vec<usize>>,
) -> bool {
    &triangle_data[*id_triangle][0] == big_id_triangle_vertexe
        || &triangle_data[*id_triangle][1] == big_id_triangle_vertexe
        || &triangle_data[*id_triangle][2] == big_id_triangle_vertexe
}

fn is_vertexe_on_right_side_of_an_edge(v1: [f32; 2], v2: [f32; 2], v3: [f32; 2]) -> bool {
    ((v2[0] - v1[0]) * (v3[1] - v1[1]) - (v2[1] - v1[1]) * (v3[0] - v1[0])) <= 0.
}

fn add_vertexe_in_triangle(
    id_triangle: usize,
    id_last_triangle_formed: usize,
    triangle_data: &mut Vec<Vec<usize>>,
    id_vertex: usize,
    normalize_vertices: &Vec<[f32; 2]>,
) {
    let triangle_id_1 = id_triangle;
    let triangle_id_2 = id_last_triangle_formed + 1;
    let triangle_id_3 = id_last_triangle_formed + 2;

    // Add all the inforamtions into the triangle_data
    let mut triangle_2 = Vec::new();
    // Vertices
    triangle_2.push(id_vertex);
    triangle_2.push(triangle_data[0][1]);
    triangle_2.push(triangle_data[0][2]);

    // Neighbours id
    triangle_2.push(triangle_id_3);
    triangle_2.push(triangle_data[0][4]);
    triangle_2.push(triangle_id_1);

    let mut triangle_3 = Vec::new();
    // Vertices
    triangle_3.push(id_vertex);
    triangle_3.push(triangle_data[0][0]);
    triangle_3.push(triangle_data[0][1]);

    // Neighbours id
    triangle_3.push(triangle_id_1);
    triangle_3.push(triangle_data[0][3]);
    triangle_3.push(triangle_id_2);

    triangle_data.push(triangle_2);
    triangle_data.push(triangle_3);

    // Update triangle indexes
    update_triangle_adjancy(
        triangle_data[id_triangle][3],
        id_triangle,
        triangle_id_3,
        triangle_data,
    );
    update_triangle_adjancy(
        triangle_data[id_triangle][4],
        id_triangle,
        triangle_id_2,
        triangle_data,
    );

    // Replace id_triangle with id_triangle_1
    triangle_data[triangle_id_1][1] = triangle_data[id_triangle][2];
    triangle_data[triangle_id_1][2] = triangle_data[id_triangle][1];
    triangle_data[triangle_id_1][0] = id_vertex;

    triangle_data[triangle_id_1][4] = triangle_data[id_triangle][5];
    triangle_data[triangle_id_1][3] = triangle_id_2;
    triangle_data[triangle_id_1][5] = triangle_id_3;

    // Make sure to still be in a Delaunay triangulation:
    udpate_delaunay_triangulation(
        id_vertex,
        triangle_id_1,
        triangle_id_2,
        triangle_id_3,
        triangle_data,
        normalize_vertices,
    );
}

fn update_triangle_adjancy(
    id_triangle: usize,
    id_triangle_to_be_raplced: usize,
    id_new_triangle: usize,
    triangle_data: &mut Vec<Vec<usize>>,
) {
    let (edge, is_edge_shared) =
        is_edge_shared(id_triangle, id_triangle_to_be_raplced, triangle_data);

    if id_triangle == 100000 {
        return;
    } else if is_edge_shared {
        triangle_data[id_triangle][edge] = id_new_triangle;
    }
}

fn is_edge_shared(
    id_triangle: usize,
    id_adjacent_triangle: usize,
    triangle_data: &Vec<Vec<usize>>,
) -> (usize, bool) {
    if id_triangle == 100000 {
        (0, false)
    } else if triangle_data[id_triangle][3] == id_adjacent_triangle {
        (3, true)
    } else if triangle_data[id_triangle][4] == id_adjacent_triangle {
        (4, true)
    } else if triangle_data[id_triangle][5] == id_adjacent_triangle {
        (5, true)
    } else {
        (0, false)
    }
}

fn udpate_delaunay_triangulation(
    id_vertex: usize,
    triangle_id_1: usize,
    triangle_id_2: usize,
    triangle_id_3: usize,
    triangle_data: &mut Vec<Vec<usize>>,
    normalize_vertices: &Vec<[f32; 2]>,
) {
    let mut triangles = Vec::<(usize, usize)>::new();

    triangles.push((triangle_id_1, triangle_data[triangle_id_1][4]));
    triangles.push((triangle_id_2, triangle_data[triangle_id_2][4]));
    triangles.push((triangle_id_3, triangle_data[triangle_id_3][4]));

    while triangles.len() > 0 {
        let (id_triangle, id_adjacent_triangle) = triangles.pop().unwrap();

        let (is_swap_quad, t3, t4) = swap_quad_diagonal_of_triangle(
            id_vertex,
            id_triangle,
            id_adjacent_triangle,
            triangle_data,
            normalize_vertices,
        );

        if id_adjacent_triangle == 100000 {
            continue;
        } else if is_swap_quad {
            triangles.push((id_triangle, t3));
            triangles.push((id_adjacent_triangle, t4));
        }
    }
}

fn swap_quad_diagonal_of_triangle(
    id_vertex: usize,
    id_triangle: usize,
    id_adjacent_triangle: usize,
    triangle_data: &mut Vec<Vec<usize>>,
    normalize_vertices: &Vec<[f32; 2]>,
) -> (bool, usize, usize) {
    let id_quad_1;
    let id_quad_2;
    let id_quad_3;
    let id_quad_4 = id_vertex;

    let triangle_id_3;
    let triangle_id_4;

    let mut is_swap_quad = false;

    // Get orientation
    if triangle_data[id_adjacent_triangle][3] == id_triangle {
        id_quad_1 = triangle_data[id_adjacent_triangle][1];
        id_quad_2 = triangle_data[id_adjacent_triangle][0];
        id_quad_3 = triangle_data[id_adjacent_triangle][2];

        triangle_id_3 = triangle_data[id_adjacent_triangle][4];
        triangle_id_4 = triangle_data[id_adjacent_triangle][5];
    } else if triangle_data[id_adjacent_triangle][4] == id_triangle {
        id_quad_1 = triangle_data[id_adjacent_triangle][2];
        id_quad_2 = triangle_data[id_adjacent_triangle][1];
        id_quad_3 = triangle_data[id_adjacent_triangle][0];

        triangle_id_3 = triangle_data[id_adjacent_triangle][5];
        triangle_id_4 = triangle_data[id_adjacent_triangle][3];
    } else {
        id_quad_1 = triangle_data[id_adjacent_triangle][0];
        id_quad_2 = triangle_data[id_adjacent_triangle][2];
        id_quad_3 = triangle_data[id_adjacent_triangle][1];

        triangle_id_3 = triangle_data[id_adjacent_triangle][3];
        triangle_id_4 = triangle_data[id_adjacent_triangle][4];
    }

    // Check if the vertexe is on the circumcircle of the adjacent triangle:
    if is_triangle_circumscribing_point(
        normalize_vertices[id_quad_1],
        normalize_vertices[id_quad_2],
        normalize_vertices[id_quad_3],
        normalize_vertices[id_quad_4],
    ) {
        update_triangle_adjancy(
            triangle_id_3,
            id_adjacent_triangle,
            id_triangle,
            triangle_data,
        );
        update_triangle_adjancy(
            triangle_data[id_triangle][5],
            id_triangle,
            id_adjacent_triangle,
            triangle_data,
        );

        triangle_data[id_triangle][0] = id_quad_4;
        triangle_data[id_triangle][1] = id_quad_1;
        triangle_data[id_triangle][2] = id_quad_3;

        triangle_data[id_adjacent_triangle][0] = id_quad_4;
        triangle_data[id_adjacent_triangle][1] = id_quad_3;
        triangle_data[id_adjacent_triangle][2] = id_quad_2;

        triangle_data[id_adjacent_triangle][3] = id_triangle;
        triangle_data[id_adjacent_triangle][4] = triangle_id_4;
        triangle_data[id_adjacent_triangle][5] = triangle_data[id_triangle][5];

        triangle_data[id_triangle][4] = triangle_id_3;
        triangle_data[id_triangle][5] = id_adjacent_triangle;

        is_swap_quad = true;
    }
    (is_swap_quad, triangle_id_3, triangle_id_4)
}

fn is_triangle_circumscribing_point(
    vertexe_1: [f32; 2],
    vertexe_2: [f32; 2],
    vertexe_3: [f32; 2],
    point: [f32; 2],
) -> bool {
    let x13 = vertexe_1[0] - vertexe_3[0];
    let x23 = vertexe_2[0] - vertexe_3[0];
    let y13 = vertexe_1[1] - vertexe_3[1];
    let y23 = vertexe_2[1] - vertexe_3[1];
    let x14 = vertexe_1[0] - point[0];
    let x24 = vertexe_2[0] - point[0];
    let y14 = vertexe_1[1] - point[1];
    let y24 = vertexe_2[1] - point[1];

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
// }
