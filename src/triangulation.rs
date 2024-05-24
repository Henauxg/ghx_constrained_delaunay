use hashbrown::HashSet;
use log::error;

use crate::types::{
    Float, Neighbor, Quad, TriangleData, TriangleId, Vector3A, VertexId, Vertice, EDGE_12, EDGE_23,
    EDGE_31, VERT_1, VERT_2, VERT_3,
};
use crate::utils::{is_point_on_right_side_of_edge, is_vertex_in_triangle_circumcircle};

#[cfg(feature = "progress_log")]
use log::info;

#[cfg(feature = "debug_context")]
use crate::debug::{DebugContext, TriangulationPhase};

pub const CONTAINER_TRIANGLE_COORDINATE: Float = 5.;

pub const CONTAINER_TRIANGLE_VERTICES: [Vertice; 3] = [
    Vertice::new(
        -CONTAINER_TRIANGLE_COORDINATE,
        -CONTAINER_TRIANGLE_COORDINATE,
    ),
    Vertice::new(0., CONTAINER_TRIANGLE_COORDINATE),
    Vertice::new(
        CONTAINER_TRIANGLE_COORDINATE,
        -CONTAINER_TRIANGLE_COORDINATE,
    ),
];

/// plane_normal must be normalized
/// vertices must all belong to a 3d plane
pub fn triangulation_from_3d_planar_vertices(
    vertices: &Vec<[Float; 3]>,
    plane_normal: Vector3A,
) -> Triangulation {
    // TODO Clean: See what we need for input data format of `triangulate`
    let mut vertices_data = Vec::with_capacity(vertices.len());
    for v in vertices {
        vertices_data.push(Vector3A::from_array(*v));
    }

    let mut planar_vertices =
        transform_to_2d_planar_coordinate_system(&mut vertices_data, plane_normal);

    triangulation_from_2d_vertices(&mut planar_vertices)
}
pub struct Triangulation {
    pub vert_indices: Vec<VertexId>,
    #[cfg(feature = "debug_context")]
    pub debug_context: DebugContext,
}

pub fn triangulation_from_2d_vertices(vertices: &Vec<Vertice>) -> Triangulation {
    // Uniformly scale the coordinates of the points so that they all lie between 0 and 1.
    let (mut normalized_vertices, _scale_factor, _x_min, _y_min) =
        normalize_vertices_coordinates(vertices);

    #[cfg(feature = "debug_context")]
    let mut debug_context = DebugContext::new(_scale_factor, _x_min, _y_min);

    let (triangles, container_triangle) = wrap_and_triangulate_2d_normalized_vertices(
        &mut normalized_vertices,
        #[cfg(feature = "debug_context")]
        &mut debug_context,
    );

    let vert_indices = remove_wrapping(
        &triangles,
        &container_triangle,
        #[cfg(feature = "debug_context")]
        &mut debug_context,
    );

    Triangulation {
        vert_indices,
        #[cfg(feature = "debug_context")]
        debug_context,
    }
}

/// Transforms 3d coordinates of all vertices into 2d coordinates on a plane defined by the given normal and vertices.
/// - Input vertices need to all belong to the same 3d plan
/// - There must be at least two vertices
pub(crate) fn transform_to_2d_planar_coordinate_system(
    vertices: &mut Vec<Vector3A>,
    plane_normal: Vector3A,
) -> Vec<Vertice> {
    // Create a base, using the first two vertices as the first base vector and plane_normal as the second
    let basis_1 = (vertices[0] - vertices[1]).normalize();
    // basis_3 is already normalized since basis_1 and plane_normal are normalized and orthogonal
    let basis_3 = basis_1.cross(plane_normal);

    // Project every vertices into the base B
    let mut vertices_2d = Vec::with_capacity(vertices.len());
    for vertex in vertices {
        vertices_2d.push(Vertice::new(vertex.dot(basis_1), vertex.dot(basis_3)));
    }
    vertices_2d
}

/// This scaling ensures that all of the coordinates are between 0 and 1 but does not modify the relative positions of the points in the x-y plane.
/// The use of normalized coordinates, although not essential, reduces the effects of  roundoff error and is also convenient from a computational point of view.
pub(crate) fn normalize_vertices_coordinates(
    vertices: &Vec<Vertice>,
) -> (Vec<Vertice>, Float, Float, Float) {
    let mut normalized_vertices = Vec::with_capacity(vertices.len());
    let (mut x_min, mut y_min, mut x_max, mut y_max) =
        (Float::MAX, Float::MAX, Float::MIN, Float::MIN);

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

    for vertex in vertices.iter() {
        normalized_vertices.push(Vertice {
            x: (vertex.x - x_min) / scale_factor,
            y: (vertex.y - y_min) / scale_factor,
        });
    }

    (normalized_vertices, scale_factor, x_min, y_min)
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
enum SearchResult {
    EnlosingTriangle(TriangleId),
    NotFound,
}

fn search_enclosing_triangle(
    sorted_vertex: &SortedVertex,
    from: Neighbor,
    triangles: &Vec<TriangleData>,
    vertices: &Vec<Vertice>,
) -> SearchResult {
    let mut triangle_id = from;

    let mut search_result = SearchResult::NotFound;
    // We use `triangles.len()` as an upper bound on the number of triangles
    for _ in 0..triangles.len() {
        let Some(current_triangle_id) = triangle_id else {
            break;
        };
        let triangle = &triangles[current_triangle_id];
        let (v1, v2, v3) = triangle.to_vertices(vertices);

        // Check if the point is inside the triangle, if not check the neighbours
        if !is_point_on_right_side_of_edge((v1, v2), sorted_vertex.vertex) {
            triangle_id = triangle.neighbors[EDGE_12];
        } else if !is_point_on_right_side_of_edge((v2, v3), sorted_vertex.vertex) {
            triangle_id = triangle.neighbors[EDGE_23];
        } else if !is_point_on_right_side_of_edge((v3, v1), sorted_vertex.vertex) {
            triangle_id = triangle.neighbors[EDGE_31];
        } else {
            search_result = SearchResult::EnlosingTriangle(current_triangle_id);
            break;
        }
    }

    search_result
}

/// Select three dummy points to form a supertriangle that completely encompasses all of the points to be triangulated.
///  This supertriangle initially defines a Delaunay triangulation which is comprised of a single triangle.
///  Its vertices are defined in terms of normalized coordinates and are usually located at a considerable distance from the window which encloses the set of points.
pub(crate) fn add_container_triangle_vertices(vertices: &mut Vec<Vertice>) -> TriangleData {
    let container_triangle = TriangleData {
        verts: [vertices.len(), vertices.len() + 1, vertices.len() + 2],
        neighbors: [None, None, None],
    };
    vertices.extend(CONTAINER_TRIANGLE_VERTICES.clone());
    container_triangle
}

/// - `vertices` should be normalized with their cooridnates in [0,1]
pub(crate) fn wrap_and_triangulate_2d_normalized_vertices(
    vertices: &mut Vec<Vertice>,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) -> (Vec<TriangleData>, TriangleData) {
    // Sort points into bins. Cover the region to be triangulated by a rectangular grid so that each bin contains roughly N^(1/2) points.
    // Label the bins so that consecutive bins are adjacent to one another, and then allocate each point to its appropriate bin.
    // Sort the list of points in ascending sequence of their bin numbers so that consecutive points are grouped together in the x-y plane.
    let partitioned_vertices = VertexBinSort::sort(&vertices, 0.5);

    let container_triangle = add_container_triangle_vertices(vertices);

    let mut triangles = Vec::<TriangleData>::new();
    triangles.push(container_triangle.clone());

    // Id of the triangle we are looking at
    let mut triangle_id = Some(0);

    #[cfg(feature = "debug_context")]
    debug_context.push_snapshot(
        TriangulationPhase::ContainerVerticesInsertion,
        &triangles,
        &[0],
        &[],
    );

    // Loop over all the input vertices
    for (_step, sorted_vertex) in partitioned_vertices.iter().enumerate() {
        #[cfg(feature = "debug_context")]
        debug_context.set_step(_step);

        // Find an existing triangle which encloses P
        match search_enclosing_triangle(sorted_vertex, triangle_id, &triangles, &vertices) {
            SearchResult::EnlosingTriangle(enclosing_triangle_id) => {
                let vertex_id = sorted_vertex.original_id;

                // Form three new triangles by connecting P to each of the enclosing triangle's vertices.
                let new_triangles = split_triangle_in_three_at_vertex(
                    &mut triangles,
                    enclosing_triangle_id,
                    vertex_id,
                    #[cfg(feature = "debug_context")]
                    debug_context,
                );

                restore_delaunay_triangulation(
                    &mut triangles,
                    vertices,
                    vertex_id,
                    new_triangles,
                    #[cfg(feature = "debug_context")]
                    debug_context,
                );

                // We'll start the search for the next enclosing triangle from the last created triangle.
                // This is a pretty good heuristic since the vertices were spatially partitionned
                triangle_id = Some(triangles.len() - 1);

                #[cfg(feature = "progress_log")]
                {
                    if _step % ((partitioned_vertices.len() / 50) + 1) == 0 {
                        let progress = 100. * _step as f32 / partitioned_vertices.len() as f32;
                        info!(
                            "Triangulation progress {}%: {}/{}",
                            progress,
                            _step,
                            partitioned_vertices.len()
                        );
                    }
                }
            }
            SearchResult::NotFound => {
                // TODO Error
                error!(
                    "Found no triangle enclosing vertex {:?}, step {}",
                    sorted_vertex, _step
                );
                break;
            }
        }
    }

    (triangles, container_triangle)
}

pub(crate) fn remove_wrapping(
    triangles: &Vec<TriangleData>,
    container_triangle: &TriangleData,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) -> Vec<VertexId> {
    // TODO Clean: Size approx
    let mut indices = Vec::with_capacity(3 * triangles.len());

    #[cfg(feature = "debug_context")]
    let mut filtered_debug_triangles = Vec::new();

    let container_verts: HashSet<VertexId> = HashSet::from(container_triangle.verts);
    for triangle in triangles.iter() {
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
            #[cfg(feature = "debug_context")]
            filtered_debug_triangles.push(triangle.clone());
        }
    }

    #[cfg(feature = "debug_context")]
    debug_context.push_snapshot(
        TriangulationPhase::RemoveWrapping,
        &filtered_debug_triangles,
        &[],
        &[],
    );

    indices
}

pub(crate) struct VertexBinSort {
    bins_per_row: usize,
    bins_count: usize,
}

#[derive(Debug, Default, Copy, Clone)]
pub(crate) struct SortedVertex {
    pub(crate) original_id: VertexId,
    pub(crate) vertex: Vertice,
}

impl VertexBinSort {
    // Each bin will contain roughly vertices.len()^(vertex_density_power) vertices
    pub fn sort(vertices: &Vec<Vertice>, vertex_density_power: Float) -> Vec<SortedVertex> {
        let bins_per_row = (vertices.len() as Float)
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

    fn bin_index_from_vertex(&self, vertex: Vertice) -> usize {
        // Compute a bin index from a vertox position which is in [0.,1]
        let bin_x = (0.99 * self.bins_per_row as Float * vertex.x) as usize;
        let bin_y: usize = (0.99 * self.bins_per_row as Float * vertex.y) as usize;
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

/// Splits `triangle_id` into 3 triangles (re-using the existing triangle id)
///
/// All the resulting triangles will share `vertex_id` as their frist vertex, and will be oriented in a CW order
///
/// ```text
///                  v1
///                / | \
///               / 3|2 \
///              /   |   \
///             /    |    \
///            / t1  |  t3 \
///           /     1|1     \
///          /      /1\      \
///         /     /     \     \
///        /    /         \    \
///       /   /             \   \
///      /2 /        t2       \ 3\
///     / / 3                 2 \ \
///   v3 ------------------------- v2
/// ```
pub(crate) fn split_triangle_in_three_at_vertex(
    triangles: &mut Vec<TriangleData>,
    triangle_id: TriangleId,
    vertex_id: VertexId,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) -> [TriangleId; 3] {
    // Re-use the existing triangle id for the first triangle
    let t1 = triangle_id;
    // Create two new triangles for the other two
    let t2 = triangles.len();
    let t3 = triangles.len() + 1;

    // t2
    triangles.push(TriangleData {
        verts: [vertex_id, triangles[t1].v2(), triangles[t1].v3()],
        neighbors: [Some(t3), triangles[t1].neighbor23(), Some(t1)],
    });
    // t3
    triangles.push(TriangleData {
        verts: [vertex_id, triangles[t1].v1(), triangles[t1].v2()],
        neighbors: [Some(t1), triangles[t1].neighbor12(), Some(t2)],
    });

    // Update triangle indexes
    update_triangle_neighbor(triangles[t1].neighbor12(), Some(t1), Some(t3), triangles);
    update_triangle_neighbor(triangles[t1].neighbor23(), Some(t1), Some(t2), triangles);

    triangles[t1].verts[VERT_2] = triangles[t1].v3();
    triangles[t1].verts[VERT_3] = triangles[t1].v1();
    triangles[t1].verts[VERT_1] = vertex_id;

    triangles[t1].neighbors[EDGE_23] = triangles[t1].neighbor31();
    triangles[t1].neighbors[EDGE_12] = Some(t2);
    triangles[t1].neighbors[EDGE_31] = Some(t3);

    #[cfg(feature = "debug_context")]
    debug_context.push_snapshot(
        TriangulationPhase::SplitTriangle(vertex_id),
        &triangles,
        &[t1, t2, t3],
        // Neighbor data is out of date here and is nto that interesting
        &[],
    );

    [t1, t2, t3]
}

pub(crate) fn update_triangle_neighbor(
    triangle_id: Neighbor,
    old_neighbour_id: Neighbor,
    new_neighbour_id: Neighbor,
    triangles: &mut Vec<TriangleData>,
) {
    match triangle_id {
        Some(triangle_id) => {
            for neighbor in triangles[triangle_id].neighbors.iter_mut() {
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
    vertices: &Vec<Vertice>,
    from_vertex_id: VertexId,
    new_triangles: [TriangleId; 3],
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) {
    let mut quads_to_check = Vec::<(TriangleId, Neighbor)>::new();

    for &from_triangle_id in &new_triangles {
        // EDGE_23 is the opposite edge of `from_vertex_id` in the 3 new triangles
        quads_to_check.push((from_triangle_id, triangles[from_triangle_id].neighbor23()));
    }

    while let Some((from_triangle_id, opposite_triangle_id)) = quads_to_check.pop() {
        let Some(opposite_triangle_id) = opposite_triangle_id else {
            continue;
        };

        match check_and_swap_quad_diagonal(
            triangles,
            vertices,
            from_vertex_id,
            from_triangle_id,
            opposite_triangle_id,
            #[cfg(feature = "debug_context")]
            debug_context,
        ) {
            QuadSwapResult::Swapped(pairs) => {
                // Place any new triangles pairs which are now opposite `vertex_id` on the stack, to be checked
                quads_to_check.extend(pairs);
            }
            QuadSwapResult::NotSwapped => (),
        }
    }
}

#[derive(PartialEq, Eq, Debug)]
pub enum QuadSwapResult {
    /// Contains the new triangle pairs to check
    Swapped([(TriangleId, Neighbor); 2]),
    NotSwapped,
}

pub fn check_and_swap_quad_diagonal(
    triangles: &mut Vec<TriangleData>,
    vertices: &Vec<Vertice>,
    from_vertex_id: VertexId,
    from_triangle_id: TriangleId,
    opposite_triangle_id: TriangleId,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) -> QuadSwapResult {
    let opposite_triangle = &triangles[opposite_triangle_id];

    // ```text
    //                q3
    //         t3   /    \   t4
    //            /   To   \
    //          /            \
    //         q1 ---------- q2
    //          \ 2        3 /
    //            \   Tf   /
    //              \ 1  /
    //                q4
    // ```
    let (quad, triangle_3_id, triangle_4_id) =
        if opposite_triangle.neighbor12() == Some(from_triangle_id) {
            (
                Quad::new([
                    opposite_triangle.v2(),
                    opposite_triangle.v1(),
                    opposite_triangle.v3(),
                    from_vertex_id,
                ]),
                opposite_triangle.neighbor23(),
                opposite_triangle.neighbor31(),
            )
        } else if opposite_triangle.neighbor23() == Some(from_triangle_id) {
            (
                Quad::new([
                    opposite_triangle.v3(),
                    opposite_triangle.v2(),
                    opposite_triangle.v1(),
                    from_vertex_id,
                ]),
                opposite_triangle.neighbor31(),
                opposite_triangle.neighbor12(),
            )
        } else {
            (
                Quad::new([
                    opposite_triangle.v1(),
                    opposite_triangle.v3(),
                    opposite_triangle.v2(),
                    from_vertex_id,
                ]),
                opposite_triangle.neighbor12(),
                opposite_triangle.neighbor23(),
            )
        };

    let quad_vertices = quad.to_vertices(vertices);

    // Check if `from_vertex_id` is on the circumcircle of `opposite_triangle`:
    if is_vertex_in_triangle_circumcircle(&quad_vertices.0[0..=3], quad_vertices.q4()) {
        // The two triangles form a convex quadrilateral whose diagonal is drawn in the wrong direction.
        // Swap this diagonal to form two new triangles so that the structure of the Delaunay triangulation
        // is locally restored.
        //
        // The quad becomes
        // ```text
        //               q3
        //         t3  / 3|2 \   t4
        //           /    |    \
        //         /      |      \
        //        q1 2  Tf|To   3 q2
        //         \      |      /
        //           \    |    /
        //             \ 1|1 /
        //               q4
        // ```
        update_triangle_neighbor(
            triangle_3_id,
            Some(opposite_triangle_id),
            Some(from_triangle_id),
            triangles,
        );
        update_triangle_neighbor(
            triangles[from_triangle_id].neighbor31(),
            Some(from_triangle_id),
            Some(opposite_triangle_id),
            triangles,
        );

        triangles[from_triangle_id].verts = [quad.v4(), quad.v1(), quad.v3()];
        triangles[opposite_triangle_id].verts = [quad.v4(), quad.v3(), quad.v2()];
        triangles[opposite_triangle_id].neighbors = [
            Some(from_triangle_id),
            triangle_4_id,
            triangles[from_triangle_id].neighbor31(),
        ];
        triangles[from_triangle_id].neighbors[EDGE_23] = triangle_3_id;
        triangles[from_triangle_id].neighbors[EDGE_31] = Some(opposite_triangle_id);

        #[cfg(feature = "debug_context")]
        debug_context.push_snapshot(
            TriangulationPhase::SwapQuadDiagonal,
            &triangles,
            &[from_triangle_id, opposite_triangle_id],
            &[triangle_3_id, triangles[from_triangle_id].neighbor31()],
        );

        QuadSwapResult::Swapped([
            (from_triangle_id, triangle_3_id),
            (opposite_triangle_id, triangle_4_id),
        ])
    } else {
        QuadSwapResult::NotSwapped
    }
}

///////////////////////////////////////////////////////////
///                                                     ///
///                        Tests                        ///
///                                                     ///
///////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::{
        triangulation::{
            check_and_swap_quad_diagonal, normalize_vertices_coordinates,
            split_triangle_in_three_at_vertex, transform_to_2d_planar_coordinate_system,
            DebugContext, QuadSwapResult,
        },
        types::{Float, TriangleData, Vector3A, Vertice},
    };

    #[test]
    fn triangulation_normalize_set_of_vertices() {
        let mut vertices = Vec::<[Float; 2]>::new();

        vertices.push([3.0, 2.0]);
        vertices.push([-1.0, 2.0]);
        vertices.push([-1.0, -2.0]);
        vertices.push([3.0, -2.0]);

        let mut planar_vertices = Vec::with_capacity(vertices.len());
        for v in &vertices {
            planar_vertices.push(Vertice::from_array(*v));
        }

        normalize_vertices_coordinates(&mut planar_vertices);

        assert_eq!(
            Vec::from([
                Vertice::from([1., 1.]),
                Vertice::from([0., 1.]),
                Vertice::from([0., 0.]),
                Vertice::from([1., 0.])
            ]),
            planar_vertices
        );
    }

    #[test]
    fn triangulation_set_to_2d_plane_vertices() {
        let mut vertices = Vec::<[Float; 3]>::new();

        vertices.push([-3., 2., 0.]);
        vertices.push([1., 2., 0.]);
        vertices.push([1., -2., 0.]);
        vertices.push([-3., -2., 0.]);

        let plane_normal = Vector3A::Z;

        let mut vertices_data = Vec::with_capacity(vertices.len());
        for v in &vertices {
            vertices_data.push(Vector3A::from_array(*v));
        }

        let planar_vertices =
            transform_to_2d_planar_coordinate_system(&mut vertices_data, plane_normal);

        assert_eq!(
            Vec::from([
                Vertice::from([3.0, 2.0]),
                Vertice::from([-1.0, 2.0]),
                Vertice::from([-1.0, -2.0]),
                Vertice::from([3.0, -2.0])
            ]),
            planar_vertices
        );
    }

    #[test]
    fn split_in_three_triangle() {
        let mut vertices = Vec::<Vertice>::new();
        vertices.push(Vertice::new(0., 0.)); // vertex to be added

        // container triangle to be splited by the vertex
        let container_triangle = TriangleData {
            verts: [vertices.len(), vertices.len() + 1, vertices.len() + 2],
            neighbors: [None, None, None],
        };

        // vertices of the container triangle
        vertices.push(Vertice::new(1., 1.));
        vertices.push(Vertice::new(1., -2.));
        vertices.push(Vertice::new(-3., 2.));

        let mut triangles = Vec::<TriangleData>::new();
        triangles.push(container_triangle);

        #[cfg(feature = "debug_context")]
        let mut debug_context = DebugContext::new(0., 0., 0.);
        let _new_triangles = split_triangle_in_three_at_vertex(
            &mut triangles,
            0,
            0,
            #[cfg(feature = "debug_context")]
            &mut debug_context,
        );

        assert_eq!(3, triangles.len());
    }

    #[test]
    fn no_swap() {
        let mut vertices = Vec::<Vertice>::new();
        vertices.push(Vertice::new(0.5, 3.));
        vertices.push(Vertice::new(-2., -2.));
        vertices.push(Vertice::new(1., -4.));
        vertices.push(Vertice::new(3., -2.));

        let triangle_1 = TriangleData {
            verts: [3, 1, 0],
            neighbors: [None, None, Some(1)],
        };

        let triangle_2 = TriangleData {
            verts: [1, 2, 3],
            neighbors: [None, None, Some(0)],
        };

        let mut triangles = Vec::<TriangleData>::new();
        triangles.push(triangle_1);
        triangles.push(triangle_2);

        #[cfg(feature = "debug_context")]
        let mut debug_context = DebugContext::new(0., 0., 0.);
        let quad_swap = check_and_swap_quad_diagonal(
            &mut triangles,
            &vertices,
            1,
            0,
            1,
            #[cfg(feature = "debug_context")]
            &mut debug_context,
        );

        assert_eq!(QuadSwapResult::NotSwapped, quad_swap);
        assert_eq!(2, triangles.len());
    }

    #[test]
    fn swap() {
        let mut vertices = Vec::<Vertice>::new();
        vertices.push(Vertice::new(0.5, 3.));
        vertices.push(Vertice::new(-2., -2.));
        vertices.push(Vertice::new(1., -4.));
        vertices.push(Vertice::new(3., -2.));

        let triangle_1 = TriangleData {
            verts: [0, 1, 2],
            neighbors: [None, None, Some(1)],
        };

        let triangle_2 = TriangleData {
            verts: [2, 3, 0],
            neighbors: [None, None, Some(0)],
        };

        let mut triangles = Vec::<TriangleData>::new();
        triangles.push(triangle_1);
        triangles.push(triangle_2);

        #[cfg(feature = "debug_context")]
        let mut debug_context = DebugContext::new(0., 0., 0.);
        let quad_swap = check_and_swap_quad_diagonal(
            &mut triangles,
            &vertices,
            1,
            0,
            1,
            #[cfg(feature = "debug_context")]
            &mut debug_context,
        );

        assert_ne!(QuadSwapResult::NotSwapped, quad_swap);
        assert_eq!(2, triangles.len());
    }
}
