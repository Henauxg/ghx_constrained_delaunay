use bevy::math::Vec3A;

use super::utils;

pub struct Triangulation {}

impl Triangulation {
    pub fn triangulate(
        vertices_to_be_triangulate: Vec<[f32; 3]>,
        normal_vec: Vec3A,
    ) -> Vec<[f32; 2]> {
        // Vec that contains the id of every triangles and there vertices + id of the three neighbour triangles:
        let mut triangle_data = Vec::<Vec<usize>>::new();

        //Transform every 3D coords of all vertices into the same 2D plane
        let mut projected_vertices = Triangulation::project_vertices_coords_in_2d_plane(
            vertices_to_be_triangulate,
            normal_vec,
        );

        // add big triangle into the triangle vertices list
        Triangulation::add_big_triangle(projected_vertices, triangle_data);

        // make sure that every vertices of all triangles are into the big triangle, by normalize coords (not big one !)
        let mut normalize_vertices =
            Triangulation::get_projected_vertices_in_big_triangle(projected_vertices);

        // Delaunay triangulation
        let mut delaunay_triangulate_vertices =
            Triangulation::triangulate_projected_vertices(normalize_vertices, triangle_data);

        // Remove extra triangles (by selecting them from big triangle vertices)
        let mut triangulate_vertices =
            Triangulation::remove_extra_triangles(delaunay_triangulate_vertices);
        triangulate_vertices
    }

    fn project_vertices_coords_in_2d_plane(
        vertices_to_be_triangulate: Vec<[f32; 3]>,
        normal_vec: Vec3A,
    ) -> Vec<[f32; 2]> {
        // First need to create the base B(e1,e2)
        let mut base_vec_1 = utils::normalize_vertexe(utils::substract_two_vertices(
            vertices_to_be_triangulate[0],
            vertices_to_be_triangulate[1],
        ));

        let mut base_vec_2 =
            utils::normalize_vertexe(utils::cross_from_vec3a(base_vec_1, normal_vec.normalize()));

        // Project every vertices into the base B
        let mut projected_vertices = Vec::new();

        for vertexe in vertices_to_be_triangulate {
            let mut x = utils::dot_from_array(vertexe, base_vec_1);
            let mut y = utils::dot_from_array(vertexe, base_vec_2);
            let mut projected_vertexe = [x, y];
            projected_vertices.push(projected_vertexe);
        }

        projected_vertices
    }

    fn add_big_triangle(mut projected_vertices: Vec<[f32; 2]>, mut triangle_data: Vec<Vec<usize>>) {
        let triangles = projected_vertices.len();

        projected_vertices.push([0., 100.]); // top vertexe
        projected_vertices.push([-100., 100.]); // bottom right vertexe
        projected_vertices.push([-100., -100.]); // bottom left vertexe

        let mut big_triangle = Vec::new();
        big_triangle.push(triangles);
        big_triangle.push(triangles + 1);
        big_triangle.push(triangles + 2);

        big_triangle.push(100000);
        big_triangle.push(100000);
        big_triangle.push(100000);

        triangle_data.push(big_triangle);
    }

    fn get_projected_vertices_in_big_triangle(
        mut projected_vertices: Vec<[f32; 2]>,
    ) -> Vec<[f32; 2]> {
        // the max/min of the coords of vertices
        let mut x_min = 0.;
        let mut y_min = 0.;
        let mut x_max = 0.;
        let mut y_max = 0.;

        for vertexe in projected_vertices {
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
        for vertexe in projected_vertices {
            let x = vertexe[0];
            let y = vertexe[1];

            let normalize_x = (x - x_min) / scale_factor;
            let normalize_y = (y - y_min) / scale_factor;

            normalize_vertices.push([normalize_x, normalize_y])
        }

        normalize_vertices
    }

    fn triangulate_projected_vertices(
        mut normalize_vertices: Vec<[f32; 2]>,
        triangle_data: Vec<Vec<usize>>,
    ) -> Vec<[f32; 2]> {
        // Total number of triangles
        let mut number_of_triangles = normalize_vertices.len() / 3;

        // Id of the triangle we are looking at
        let mut id_triangle = 0;

        // The id of the last triangle created
        let mut id_last_triangle_formed = 0;

        for (id_vertex, vertexe) in normalize_vertices.iter().enumerate() {
            let mut is_vertexe_inserted = false;
            let mut counter = 0;

            while !is_vertexe_inserted {
                counter += 1;
                // We can't outreach the limit of triangles
                if counter > number_of_triangles || id_triangle == 100000 {
                    break;
                }

                // vertices of the actual triangle
                let v1 = normalize_vertices[triangle_data[id_triangle][0]];
                let v2 = normalize_vertices[triangle_data[id_triangle][1]];
                let v3 = normalize_vertices[triangle_data[id_triangle][2]];

                // Check if the point is inside the triangle, if not check the neighbours
                if !Triangulation::is_vertexe_on_right_side_of_an_edge(v1, v2, *vertexe) {
                    id_triangle = triangle_data[id_triangle][3];
                }
                if !Triangulation::is_vertexe_on_right_side_of_an_edge(v2, v3, *vertexe) {
                    id_triangle = triangle_data[id_triangle][4];
                }
                if !Triangulation::is_vertexe_on_right_side_of_an_edge(v3, v1, *vertexe) {
                    id_triangle = triangle_data[id_triangle][5];
                } else {
                    Triangulation::add_vertexe_in_triangle(
                        *vertexe,
                        id_triangle,
                        id_last_triangle_formed,
                        triangle_data,
                        id_vertex,
                    );
                    number_of_triangles += 2;
                    id_triangle = id_last_triangle_formed;
                    is_vertexe_inserted = true;
                }
            }
        }

        todo!()
    }

    fn remove_extra_triangles(mut delaunay_triangulate_vertices: Vec<[f32; 2]>) -> Vec<[f32; 2]> {
        todo!()
    }

    fn is_vertexe_on_right_side_of_an_edge(v1: [f32; 2], v2: [f32; 2], v3: [f32; 2]) -> bool {
        ((v2[0] - v1[0]) * (v3[1] - v1[1]) - (v2[1] - v1[1]) * (v3[0] - v1[0])) <= 0.
    }

    fn add_vertexe_in_triangle(
        vertexe: [f32; 2],
        id_triangle: usize,
        id_last_triangle_formed: usize,
        triangle_data: Vec<Vec<usize>>,
        id_vertex: usize,
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
    }
}
