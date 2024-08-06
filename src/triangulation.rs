use crate::infinite::{
    collect_infinite_triangle_vertices, is_finite, is_infinite, is_vertex_in_half_plane_1,
    vertex_placement_1_infinite_vertex, vertex_placement_2_infinite_vertex, INFINITE_V0_ID,
    INFINITE_V1_ID, INFINITE_V2_ID, INFINITE_V3_ID,
};
use crate::types::{
    next_clockwise_edge_index, next_counter_clockwise_edge_index, opposite_vertex_index,
    vertex_next_ccw_edge_index, vertex_next_cw_edge_index, Float, Neighbor, Quad, TriangleData,
    TriangleEdgeIndex, TriangleId, Triangles, Vertex, Vertex2d, Vertex3d, VertexId, EDGE_12,
    EDGE_23, EDGE_31, VERT_1, VERT_2,
};
use crate::utils::{is_vertex_in_triangle_circumcircle, test_point_edge_side};

#[cfg(feature = "parallel_filtering")]
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};

#[cfg(feature = "progress_log")]
use tracing::info;

#[cfg(feature = "debug_context")]
use crate::debug::{DebugConfiguration, DebugContext, EventInfo, Phase};

#[cfg(feature = "profile_traces")]
use tracing::{span, Level};

/// Binsort will cover the region to be triangulated by a rectangular grid so that each bin contains roughly N^(density_power) points.
pub const DEFAULT_BIN_VERTEX_DENSITY_POWER: f64 = 0.5;

/// Good default value for `filter_parallel_tri_count_threshold` in [`TriangulationConfiguration`]
pub const DEFAULT_ILTER_PARALLEL_TRI_COUNT_THRESHOLD: usize = 100_000;
/// Good default value for `filter_parallel_min_batch_len` in [`TriangulationConfiguration`]
pub const DEFAULT_FILTER_PARALLEL_MIN_BATCH_LEN: usize = 1000;

#[derive(Clone, Debug)]
pub struct TriangulationConfiguration {
    /// Binsort will cover the region to be triangulated by a rectangular grid so that each bin contains roughly N^(density_power) points.
    pub bin_vertex_density_power: f64,

    /// Minimum count of triangles that needs to be reached in order to enable parallel filtering.
    #[cfg(feature = "parallel_filtering")]
    pub filter_parallel_tri_count_threshold: usize,
    /// Minimum size of the work batches during parallel triangles filtering.
    #[cfg(feature = "parallel_filtering")]
    pub filter_parallel_min_batch_len: usize,

    #[cfg(feature = "debug_context")]
    pub debug_config: DebugConfiguration,
}
impl Default for TriangulationConfiguration {
    fn default() -> Self {
        Self {
            bin_vertex_density_power: DEFAULT_BIN_VERTEX_DENSITY_POWER,

            #[cfg(feature = "parallel_filtering")]
            filter_parallel_tri_count_threshold: DEFAULT_ILTER_PARALLEL_TRI_COUNT_THRESHOLD,
            #[cfg(feature = "parallel_filtering")]
            filter_parallel_min_batch_len: DEFAULT_FILTER_PARALLEL_MIN_BATCH_LEN,

            #[cfg(feature = "debug_context")]
            debug_config: DebugConfiguration::default(),
        }
    }
}

#[derive(Default)]
pub struct Triangulation {
    /// Indices of the original vertices by groups of 3 to from triangles.
    pub triangles: Vec<[VertexId; 3]>,

    #[cfg(feature = "debug_context")]
    pub debug_context: DebugContext,
}

/// An internal error occured, which could be due to an invalid input, or possibly an error in the algorithm.
#[derive(thiserror::Error, Debug)]
#[error("triangulation error")]
pub struct TriangulationError {
    pub msg: String,
    #[cfg(feature = "debug_context")]
    pub debug_context: DebugContext,
}
impl TriangulationError {
    pub(crate) fn new(
        msg: String,
        #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
    ) -> Self {
        Self {
            msg,
            #[cfg(feature = "debug_context")]
            debug_context: debug_context.clone(),
        }
    }
}

/// Same as [triangulation_from_2d_vertices] but input vertices are in 3d and will be transformed to 2d before the triangulation.
///
/// Additional requirements:
/// - All the vertices are expected to belong to the same 2d plane, with the provided `plane_normal`.
/// - `plane_normal` must be normalized
pub fn triangulation_from_3d_planar_vertices<T: Vertex3d>(
    vertices: &Vec<T>,
    plane_normal: T,
    config: TriangulationConfiguration,
) -> Result<Triangulation, TriangulationError> {
    if vertices.len() < 3 {
        return Ok(Triangulation::default());
    }
    let mut planar_vertices = transform_to_2d_planar_coordinate_system(vertices, plane_normal);
    triangulation_from_2d_vertices(&mut planar_vertices, config)
}

/// Creates a Delaunay triangulation with the input vertices.
///
/// Vertices that are identical (or extremely close to one another) will me "merged" together. This means that only one of those will appear in the final triangulation.
///
/// Vertices requirements:
/// - Vertices are expected to be valid floating points values. You can use [crate::utils::validate_vertices] to check your vertices beforehand for NaN or infinity.
pub fn triangulation_from_2d_vertices<T: Vertex2d>(
    vertices: &Vec<T>,
    config: TriangulationConfiguration,
) -> Result<Triangulation, TriangulationError> {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "triangulation_from_2d_vertices").entered();

    if vertices.len() < 3 {
        return Ok(Triangulation::default());
    }

    // Uniformly scale the coordinates of the points so that they all lie between 0 and 1.
    let (mut normalized_vertices, _scale_factor, _x_min, _y_min) =
        normalize_vertices_coordinates(vertices);

    #[cfg(feature = "debug_context")]
    let mut debug_context =
        DebugContext::new(config.debug_config.clone(), _scale_factor, _x_min, _y_min);

    let triangles = wrap_and_triangulate_2d_normalized_vertices(
        &mut normalized_vertices,
        config.bin_vertex_density_power,
        &mut None,
        #[cfg(feature = "debug_context")]
        &mut debug_context,
    )?;

    let vert_indices = filter_triangles(
        &triangles,
        &config,
        #[cfg(feature = "debug_context")]
        &mut debug_context,
    );

    Ok(Triangulation {
        triangles: vert_indices,
        #[cfg(feature = "debug_context")]
        debug_context,
    })
}

/// Transforms 3d coordinates of all vertices into 2d coordinates on a plane defined by the given normal and vertices.
/// - Input vertices need to all belong to the same 3d plan
/// - There must be at least two vertices
pub fn transform_to_2d_planar_coordinate_system<T: Vertex3d>(
    vertices: &Vec<T>,
    plane_normal: T,
) -> Vec<Vertex> {
    // Create a base B for the plan, using the first two vertices as the first base vector
    let basis_1 = (vertices[0] - vertices[1]).normalize();
    let basis_2 = basis_1.cross(-plane_normal);
    // basis_3 would be plane_normal

    // Transform every vertices into base B 2D coordinates
    let mut vertices_2d = Vec::with_capacity(vertices.len());
    for vertex in vertices {
        vertices_2d.push(Vertex::new(vertex.dot(basis_1), vertex.dot(basis_2)));
    }
    vertices_2d
}

/// This scaling ensures that all of the coordinates are between 0 and 1 but does not modify the relative positions of the points in the x-y plane.
/// The use of normalized coordinates, although not essential, reduces the effects of roundoff error and is also convenient from a computational point of view.
pub(crate) fn normalize_vertices_coordinates<T: Vertex2d>(
    vertices: &Vec<T>,
) -> (Vec<Vertex>, Float, Float, Float) {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "normalize_vertices_coordinates").entered();

    let mut normalized_vertices = Vec::with_capacity(vertices.len());
    let (mut x_min, mut y_min, mut x_max, mut y_max) =
        (Float::MAX, Float::MAX, Float::MIN, Float::MIN);

    for vertex in vertices.iter() {
        if vertex.x() < x_min {
            x_min = vertex.x();
        }
        if vertex.x() > x_max {
            x_max = vertex.x();
        }
        if vertex.y() < y_min {
            y_min = vertex.y();
        }
        if vertex.y() > y_max {
            y_max = vertex.y();
        }
    }

    let scale_factor = (x_max - x_min).max(y_max - y_min);

    for vertex in vertices.iter() {
        normalized_vertices.push(Vertex {
            x: (vertex.x() - x_min) / scale_factor,
            y: (vertex.y() - y_min) / scale_factor,
        });
    }

    (normalized_vertices, scale_factor, x_min, y_min)
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub(crate) enum VertexPlacement {
    InsideTriangle(TriangleId),
    OnTriangleEdge(TriangleId, TriangleEdgeIndex),
    OnVertex(VertexId),
}

fn find_vertex_placement(
    vertex: Vertex,
    from: TriangleId,
    triangles: &Triangles,
    vertices: &Vec<Vertex>,
    // #[cfg(feature = "debug_context")] _debug_context: &mut DebugContext,
) -> Option<VertexPlacement> {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "find_vertex_placement").entered();

    let mut current_triangle: Neighbor = from.into();
    let mut previous_triangle = current_triangle;

    // We use `triangles.len()` as an upper bound on the number of triangles
    for _index in 0..triangles.count() {
        if !current_triangle.exists() {
            break;
        }

        let triangle_id = current_triangle.id;
        let triangle = triangles.get(triangle_id);

        let infinite_verts = collect_infinite_triangle_vertices(&triangle.verts);

        if infinite_verts.is_empty() {
            let (v1, v2, v3) = triangle.to_vertices(vertices);

            // TODO Profile this check's position
            if is_vertex_pair_too_close(v1, vertex) {
                return Some(VertexPlacement::OnVertex(triangle.v1()));
            }
            if is_vertex_pair_too_close(v2, vertex) {
                return Some(VertexPlacement::OnVertex(triangle.v2()));
            }
            if is_vertex_pair_too_close(v3, vertex) {
                return Some(VertexPlacement::OnVertex(triangle.v3()));
            }

            // Check if the point is towards this triangle's neighbor
            // Check previous triangle to avoid looping between two triangles
            let edge_12_test = test_point_edge_side((v1, v2), vertex);
            if edge_12_test.is_on_left_side() && previous_triangle != triangle.neighbor12() {
                previous_triangle = current_triangle;
                current_triangle = triangle.neighbor12();
                continue;
            }
            let edge_23_test = test_point_edge_side((v2, v3), vertex);
            if edge_23_test.is_on_left_side() && previous_triangle != triangle.neighbor23() {
                previous_triangle = current_triangle;
                current_triangle = triangle.neighbor23();
                continue;
            }
            let edge_31_test = test_point_edge_side((v3, v1), vertex);
            if edge_31_test.is_on_left_side() && previous_triangle != triangle.neighbor31() {
                previous_triangle = current_triangle;
                current_triangle = triangle.neighbor31();
                continue;
            }

            if edge_12_test.is_near_edge() {
                return Some(VertexPlacement::OnTriangleEdge(triangle_id, EDGE_12));
            } else if edge_23_test.is_near_edge() {
                return Some(VertexPlacement::OnTriangleEdge(triangle_id, EDGE_23));
            } else if edge_31_test.is_near_edge() {
                return Some(VertexPlacement::OnTriangleEdge(triangle_id, EDGE_31));
            } else {
                return Some(VertexPlacement::InsideTriangle(triangle_id));
            }
        } else if infinite_verts.len() == 1 {
            let res = vertex_placement_1_infinite_vertex(
                vertices,
                vertex,
                triangle,
                triangle_id,
                infinite_verts,
                &mut previous_triangle,
                &mut current_triangle,
            );
            match res {
                Some(_) => return res,
                None => continue,
            }
        } else if infinite_verts.len() == 2 {
            let res = vertex_placement_2_infinite_vertex(
                vertices,
                vertex,
                triangle,
                triangle_id,
                infinite_verts,
                &mut previous_triangle,
                &mut current_triangle,
            );
            match res {
                Some(_) => return res,
                None => continue,
            }
        }
    }

    None
}

/// Splits the infinite quad into 4 triangles at p, the first inserted vertex.
///
/// ```text
///          1
///          |
///     t1   |    t2
///          |
/// 0--------p---------2
///          |
///     t4   |    t3
///          |
///          3
/// ```
fn create_initial_triangles(vertices_count: usize, first_vertex_id: VertexId) -> Triangles {
    let mut triangles = Triangles::with_capacity(vertices_count * 2 + 2);

    let (t1, t2, t3, t4) = (0.into(), 1.into(), 2.into(), 3.into());
    triangles.buffer_mut().push(TriangleData {
        verts: [first_vertex_id, INFINITE_V0_ID, INFINITE_V1_ID],
        neighbors: [t4, Neighbor::NONE, t2],
    });
    triangles.buffer_mut().push(TriangleData {
        verts: [first_vertex_id, INFINITE_V1_ID, INFINITE_V2_ID],
        neighbors: [t1, Neighbor::NONE, t3],
    });
    triangles.buffer_mut().push(TriangleData {
        verts: [first_vertex_id, INFINITE_V2_ID, INFINITE_V3_ID],
        neighbors: [t2, Neighbor::NONE, t4],
    });
    triangles.buffer_mut().push(TriangleData {
        verts: [first_vertex_id, INFINITE_V3_ID, INFINITE_V0_ID],
        neighbors: [t3, Neighbor::NONE, t1],
    });
    triangles
}

#[inline]
pub(crate) fn is_vertex_pair_too_close(a: Vertex, b: Vertex) -> bool {
    let dist = a - b;
    dist.x.abs() < Float::EPSILON && dist.y.abs() < Float::EPSILON
}

fn verts_count_to_restoration_stack_initial_capacity(verts_count: usize) -> usize {
    if verts_count < 9_000 {
        16
    } else if verts_count < 150_000 {
        32
    } else if verts_count < 500_000 {
        64
    } else {
        128
    }
}

/// - `vertices` should be normalized with their cooridnates in [0,1]
pub(crate) fn wrap_and_triangulate_2d_normalized_vertices(
    vertices: &Vec<Vertex>,
    bin_vertex_density_power: f64,
    vertex_merge_mapping: &mut Option<Vec<VertexId>>,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) -> Result<Triangles, TriangulationError> {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "wrap_and_triangulate_2d_normalized_vertices").entered();

    // Sort points into bins. Cover the region to be triangulated by a rectangular grid so that each bin contains roughly N^(1/2) points.
    // Label the bins so that consecutive bins are adjacent to one another, and then allocate each point to its appropriate bin.
    // Sort the list of points in ascending sequence of their bin numbers so that consecutive points are grouped together in the x-y plane.
    let partitioned_vertices = VertexBinSort::sort(&vertices, bin_vertex_density_power);

    let mut vertices_iterator = partitioned_vertices.iter().enumerate().peekable();
    // There is always at least 3 vertices
    let (_, first_vertex_id) = vertices_iterator.next().unwrap();
    let mut triangles = create_initial_triangles(vertices.len(), *first_vertex_id);
    // Id of the triangle we are looking at
    let mut triangle_id = triangles.last_id();

    #[cfg(feature = "debug_context")]
    debug_context.push_snapshot(Phase::FirstVertexInsertion, &triangles, &[0], &[]);

    // This buffer is used in the loop. Filled and cleared each time.
    // We create if here to share the allocation between all those calls as an optimization.
    let mut quads_to_check = Vec::with_capacity(verts_count_to_restoration_stack_initial_capacity(
        vertices.len(),
    ));

    // let triangles_double_buffer = TrianglesDoubleBuffer {
    //     b1: triangles.clone(),
    //     b2: triangles,
    // };

    // let (to_predictor_send, to_predictor_recv) =
    //     mpsc::channel::<(VertexId, TriangleId, Vec<(TriangleId,TriangleData)>)>();
    // let (from_predictor_send, from_predictor_recv) =
    //     mpsc::channel::<Option<VertexPlacement>>();

    // // Allocate/copy once.
    // // TODO Could we avoid this clone with lifetimes ? (read-only here)
    // let predictor_vertices = vertices.clone();

    // // TODO Some measures: when should we start predicting vertex placement? Right at the start, or from a certain number of triangles ?
    // {
    //     let (_, &next_vertex_id) = vertices_iterator.peek().unwrap();
    //     match to_predictor_send.send((next_vertex_id, triangle_id, Vec::new())) {
    //         Ok(_) => (),
    //         // TODO Handle error
    //         Err(_) => println!("Failed to send vertex placement request to the predictor"),
    //     };

    //     // Allocate once, and then memcpy to refresh
    //     let mut triangles_mirror = Triangles::with_capacity(vertices.len() * 2 + 2);
    //     triangles_mirror
    //         .buffer
    //         .extend_from_slice(&triangles.buffer);

    //     thread::spawn(move ||
    //         // triangles_snapshot is a snapshot of the triangles at some point in the past.
    //         while let Ok((vertex_id, triangle_id, triangles_deltas)) = to_predictor_recv.recv() {
    //             for (triangle_id, triangle) in triangles_deltas.into_iter(){
    //                 // TODO Constraint, new ids are in order in triangles_deltas (so that i can just push and assume it is ok)
    //                 if triangle_id as usize >= triangles_mirror.count() {
    //                     triangles_mirror.push(triangle);
    //                 } else {
    //                     *triangles_mirror.get_mut(triangle_id) = triangle;
    //                 }
    //             }

    //             let vertex = predictor_vertices[vertex_id as usize];
    //             let predicted_vertex_place = find_vertex_placement(
    //                 vertex,
    //                 triangle_id,
    //                 &triangles_mirror,
    //                 &predictor_vertices,
    //             );
    //             if let Err(_) = from_predictor_send.send(predicted_vertex_place) {
    //                 // TODO Handle error
    //                 println!("Failed to send vertex placement response from the predictor")
    //             }
    //         }
    //     );
    // }

    // TODO Performance: triangle_deltas could be double buffered to avoid allocations (since we clear it each iteration)
    // let mut triangle_deltas = Vec::with_capacity(25);
    // let mut modified_triangles = BTreeSet::new();

    // Loop over all the input vertices
    while let Some((_index, &vertex_id)) = vertices_iterator.next() {
        // for (_index, &vertex_id) in vertices_iterator {
        #[cfg(feature = "debug_context")]
        {
            let force_end = debug_context.advance_step();
            if force_end {
                break;
            }
        }

        // TODO Handle the first loop (no prediction pending? or launch one early ? or ?)
        // TODO Wait for the previous search job to be finished + Check its results
        // let Ok(predicted_vertex_place) = from_predictor_recv.recv() else {
        //     // TODO Handle this case
        //     println!("Error on prediction recv");
        //    break;
        // };
        // if predicted_vertex_place.is_none()  {
        //     println!("Internal error, failed to predict placement for vertex {:?}, step {}",
        //     vertex_id, _index);
        // };

        // If there is a vertex after the current one, queue a prediction
        // if let Some((_, &next_vertex_id)) = vertices_iterator.peek() {
        //     // TODO Better
        //     match to_predictor_send.send((next_vertex_id, triangle_id, mem::take(&mut triangle_deltas))) {
        //         Ok(_) => (),
        //         // TODO Handle error
        //         Err(_) => println!("Failed to send vertex placement request to the predictor"),
        //     };
        // }

        // TODO Prediciton verifier: check if the result is valid
        // OnVertex => no need to check, will be true.
        // OnEdge: TODO
        // InsideTriangle: TODO
        // If the prediciton is wrong: still, do a vertex placement search from the predicted triangle (which should be better than just the last inserted triangle)

        // Find an existing triangle which encloses P
        let vertex = vertices[vertex_id as usize];
        let Some(vertex_place) = find_vertex_placement(
            vertex,
            triangle_id,
            &triangles,
            &vertices,
            // #[cfg(feature = "debug_context")]
            // debug_context,
        ) else {
            return Err(TriangulationError::new(
                format!(
                    "Internal error, failed to find placement for vertex {:?}, step {}",
                    vertex_id, _index
                ),
                #[cfg(feature = "debug_context")]
                debug_context,
            ));
        };

        // TODO Tmp
        // if vertex_place == predicted_vertex_place {
        //     println!("Accurate vertex prediction, step {}",_index);
        // } else {
        //     println!("Failed vertex prediction, step {}, vertex_place {:?}, predicted_vertex_place {:?}",_index,vertex_place,predicted_vertex_place);
        // }

        match vertex_place {
            VertexPlacement::InsideTriangle(enclosing_triangle_id) => {
                split_triangle_into_three_triangles(
                    &mut triangles,
                    enclosing_triangle_id,
                    vertex_id,
                    &mut quads_to_check,
                    // &mut modified_triangles,
                    // &mut triangle_deltas,
                    #[cfg(feature = "debug_context")]
                    debug_context,
                )
            }
            VertexPlacement::OnTriangleEdge(enclosing_triangle_id, edge_index) => {
                split_quad_into_four_triangles(
                    &mut triangles,
                    enclosing_triangle_id,
                    edge_index,
                    vertex_id,
                    &mut quads_to_check,
                    // &mut triangle_deltas,
                    #[cfg(feature = "debug_context")]
                    debug_context,
                )
            }
            VertexPlacement::OnVertex(existing_vertex_id) => {
                if let Some(vertex_merge_mapping) = vertex_merge_mapping {
                    vertex_merge_mapping[vertex_id as usize] = existing_vertex_id;
                }
                continue;
            }
        };

        restore_delaunay_triangulation(
            &mut triangles,
            vertices,
            vertex_id,
            &mut quads_to_check,
            // &mut triangle_deltas,
            #[cfg(feature = "debug_context")]
            debug_context,
        );

        // We'll start the search for the next enclosing triangle from the last created triangle.
        // This is a pretty good heuristic since the vertices were spatially partitionned
        triangle_id = triangles.last_id();

        #[cfg(feature = "progress_log")]
        {
            if _index % ((partitioned_vertices.len() / 50) + 1) == 0 {
                let progress = 100. * _index as f32 / partitioned_vertices.len() as f32;
                info!(
                    "Triangulation progress, step n°{}, {}%: {}/{}",
                    debug_context.current_step,
                    progress,
                    _index,
                    partitioned_vertices.len()
                );
            }
        }
    }

    Ok(triangles)
}

pub(crate) fn filter_triangles(
    triangles: &Triangles,
    _config: &TriangulationConfiguration,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) -> Vec<[VertexId; 3]> {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "filter_triangles").entered();

    #[cfg(feature = "debug_context")]
    let mut filtered_debug_triangles = Triangles::new();

    #[cfg(feature = "parallel_filtering")]
    let indices = _multithreaded_filter_triangles(
        triangles,
        _config,
        #[cfg(feature = "debug_context")]
        &mut filtered_debug_triangles,
    );
    #[cfg(not(feature = "parallel_filtering"))]
    let indices = _singlethreaded_filter_triangles(
        triangles,
        #[cfg(feature = "debug_context")]
        &mut filtered_debug_triangles,
    );

    #[cfg(feature = "debug_context")]
    debug_context.push_snapshot(Phase::FilterTriangles, &filtered_debug_triangles, &[], &[]);

    indices
}

#[cfg(feature = "parallel_filtering")]
pub(crate) fn _multithreaded_filter_triangles(
    triangles: &Triangles,
    config: &TriangulationConfiguration,
    #[cfg(feature = "debug_context")] filtered_debug_triangles: &mut Triangles,
) -> Vec<[VertexId; 3]> {
    if triangles.count() > config.filter_parallel_tri_count_threshold {
        // Debug loop out of the // loop because filtered_debug_triangles cannot be accessed as mut from multiple tasks.
        #[cfg(feature = "debug_context")]
        {
            filtered_debug_triangles.buffer = triangles
                .buffer()
                .iter()
                .filter(|t| t.is_finite())
                .map(|t| t.clone())
                .collect();
        }

        triangles
            .buffer()
            .par_iter()
            .with_min_len(config.filter_parallel_min_batch_len)
            .filter_map(|t| match t.is_finite() {
                true => Some(t.verts),
                false => None,
            })
            // TODO Is it possible to pre-set the capacity of indices ?
            // .collect_into_vec(&mut indices);
            .collect()
    } else {
        let mut indices = Vec::with_capacity(triangles.count());
        for t in triangles.buffer().iter() {
            if t.is_finite() {
                indices.push(t.verts);
                #[cfg(feature = "debug_context")]
                filtered_debug_triangles.push(t.clone());
            }
        }
        indices
    }
}

pub(crate) fn _singlethreaded_filter_triangles(
    triangles: &Triangles,
    #[cfg(feature = "debug_context")] filtered_debug_triangles: &mut Triangles,
) -> Vec<[VertexId; 3]> {
    let mut indices = Vec::with_capacity(triangles.count());
    for t in triangles.buffer().iter() {
        if t.is_finite() {
            indices.push(t.verts);
            #[cfg(feature = "debug_context")]
            filtered_debug_triangles.push(t.clone());
        }
    }
    indices
}

pub(crate) struct VertexBinSort {
    bins_per_row: usize,
    bins_count: usize,
}

impl VertexBinSort {
    // Each bin will contain roughly vertices.len()^(vertex_density_power) vertices
    pub fn sort(vertices: &Vec<Vertex>, vertex_density_power: f64) -> Vec<VertexId> {
        #[cfg(feature = "profile_traces")]
        let _span = span!(Level::TRACE, "sort").entered();

        let bins_per_row = (vertices.len() as f64)
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

        let mut sorted = vec![0; vertices.len()];
        for vertex_id in 0..vertices.len() {
            let bin_index = vertices_bin_indexes[vertex_id];
            bins_counters[bin_index] -= 1;
            sorted[bins_counters[bin_index]] = vertex_id as VertexId;
        }
        sorted
    }

    fn bin_index_from_vertex(&self, vertex: Vertex) -> usize {
        // Compute a bin index from a vertex position which is in [0.,1]
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

/// Splits a quad (triangle pair) into 4 triangles
///```text
///      from
///      /|\
///     / | \
/// n2 /  |  \ n3
///   /   |   \
///  /    |    \
/// /  t1 | t2  \
/// \     |     /
///  \    |    /
///   \   |   /
/// n1 \  |  / n4
///     \ | /
///      \|/
///       to
///
/// Becomes
///
///      /|\
///     /3|2\
/// n2 /  |  \ n3
///   / t1|t2 \
///  /    |    \
/// /2___1|1___3\
/// \3   1|1   2/
///  \ t3 | t4 /
///   \   |   /
/// n1 \  |  / n4
///     \2|3/
///      \|/
/// ```
pub(crate) fn split_quad_into_four_triangles(
    triangles: &mut Triangles,
    triangle_id: TriangleId,
    edge_index: TriangleEdgeIndex,
    vertex_id: VertexId,
    quads_to_check: &mut Vec<(TriangleId, TriangleId)>,
    // triangle_deltas: &mut Vec<(TriangleId, TriangleData)>,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "split_quad_into_four_triangles").entered();

    // Re-use the existing triangle id for the first triangle
    let t1 = triangle_id;
    // TODO Comment is not true anymore in the general sense, but still true in this function
    // The triangle is always guaranteed to have a neighbor on the edge since a vertex
    // cannot land on a container triangle's edge.
    let t2 = triangles.get(triangle_id).neighbor(edge_index).id;
    // Create two new triangles for the other two
    let t3 = triangles.next_id();
    let t4 = triangles.next_id() + 1;

    let n1 = triangles
        .get(t1)
        .neighbor(next_clockwise_edge_index(edge_index));
    let n2 = triangles
        .get(t1)
        .neighbor(next_counter_clockwise_edge_index(edge_index));

    let edge = triangles.get(t1).edge(edge_index);

    // t3
    let t3_v3 = triangles.get(t1).v(opposite_vertex_index(edge_index));
    triangles.create([vertex_id, edge.to, t3_v3], [t4.into(), n1, t1.into()]);

    // TODO Optim: Can do a tiny bit better
    let t2_opposite_v_id = triangles.get(t2).get_opposite_vertex_id(&edge);
    let t2_opposite_v_index = triangles.get(t2).vertex_index(t2_opposite_v_id);

    let n4 = triangles
        .get(t2)
        .neighbor(vertex_next_cw_edge_index(t2_opposite_v_index));
    let n3 = triangles
        .get(t2)
        .neighbor(vertex_next_ccw_edge_index(t2_opposite_v_index));

    // t4
    triangles.create(
        [vertex_id, t2_opposite_v_id, edge.to],
        [t2.into(), n4, t3.into()],
    );

    // Update t1 verts
    let verts = [vertex_id, t3_v3, edge.from];
    triangles.get_mut(t1).verts = verts;
    // Update t1 neighbors
    let neighbors = [t3.into(), n2, t2.into()];
    triangles.get_mut(t1).neighbors = neighbors;

    // Update t2 verts
    let verts = [vertex_id, edge.from, t2_opposite_v_id];
    triangles.get_mut(t2).verts = verts;
    // Update t2 neighbors
    let neighbors = [t1.into(), n3, t4.into()];
    triangles.get_mut(t2).neighbors = neighbors;

    // triangle_deltas.push((t1, triangles.get(t1).clone()));
    // triangle_deltas.push((t2, triangles.get(t2).clone()));
    // triangle_deltas.push((t3, triangles.get(t3).clone()));
    // triangle_deltas.push((t4, triangles.get(t4).clone()));

    // Update neighbors and fill in quads_to_check directly for later. Way better for performances.
    if n1.exists() {
        update_triangle_neighbor(n1.id, t1.into(), t3.into(), triangles);
        quads_to_check.push((t3, n1.id));
        // triangle_deltas.push((n1.id, triangles.get(n1.id).clone()));
    }
    if n4.exists() {
        update_triangle_neighbor(n4.id, t2.into(), t4.into(), triangles);
        quads_to_check.push((t4, n4.id));
        // triangle_deltas.push((n4.id, triangles.get(n4.id).clone()));
    }
    if n2.exists() {
        quads_to_check.push((t1, n2.id));
    }
    if n3.exists() {
        quads_to_check.push((t2, n3.id));
    }

    #[cfg(feature = "debug_context")]
    debug_context.push_snapshot_event(
        Phase::SplitQuad,
        EventInfo::Split(vertex_id),
        &triangles,
        &[t1, t2, t3, t4],
        // Neighbor data is out of date here and is not that interesting
        &[],
    );
}

/// Splits `triangle_id` into 3 triangles (re-using the existing triangle id)
///
/// All the resulting triangles will share `vertex_id` as their first vertex, and will be oriented in a CW order
///
/// ```text
///                  v1
///                / | \
///               / 3|2 \
///              /   |   \
///             /    |    \
///       n6   / t1  |  t3 \   n4
///           /     1|1     \
///          /      /1\      \
///         /     /     \     \
///        /    /         \    \
///       /   /             \   \
///      /2 /        t2       \ 3\
///     / / 3                 2 \ \
///   v3 ------------------------- v2
///                  n5
/// ```
pub(crate) fn split_triangle_into_three_triangles(
    triangles: &mut Triangles,
    triangle_id: TriangleId,
    vertex_id: VertexId,
    quads_to_check: &mut Vec<(TriangleId, TriangleId)>,
    // triangle_deltas: &mut Vec<(TriangleId, TriangleData)>,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "split_triangle_into_three_triangles").entered();

    // Re-use the existing triangle id for the first triangle
    let t1 = triangle_id;
    // Create two new triangles for the other two
    let t2 = triangles.next_id();
    let t3 = triangles.next_id() + 1;
    let n4 = triangles.get(t1).neighbor12();
    let n5 = triangles.get(t1).neighbor23();
    let n6 = triangles.get(t1).neighbor31();

    // t2
    triangles.create(
        [vertex_id, triangles.get(t1).v2(), triangles.get(t1).v3()],
        [t3.into(), triangles.get(t1).neighbor23(), t1.into()],
    );
    // t3
    triangles.create(
        [vertex_id, triangles.get(t1).v1(), triangles.get(t1).v2()],
        [t1.into(), n4, t2.into()],
    );

    // Update t1 verts
    let verts = [vertex_id, triangles.get(t1).v3(), triangles.get(t1).v1()];
    triangles.get_mut(t1).verts = verts;
    // Update t1 neighbors
    let neighbors = [t2.into(), n6, t3.into()];
    triangles.get_mut(t1).neighbors = neighbors;

    // triangle_deltas.push((t1, triangles.get(t1).clone()));
    // triangle_deltas.push((t2, triangles.get(t2).clone()));
    // triangle_deltas.push((t3, triangles.get(t3).clone()));

    // Update neighbors and fill in quads_to_check directly for later. Way better for performances.
    if n4.exists() {
        update_triangle_neighbor(n4.id, t1.into(), t3.into(), triangles);
        quads_to_check.push((t3, n4.id));
        // triangle_deltas.push((n4.id, triangles.get(n4.id).clone()));
    }
    if n5.exists() {
        update_triangle_neighbor(n5.id, t1.into(), t2.into(), triangles);
        quads_to_check.push((t2, n5.id));
        // triangle_deltas.push((n5.id, triangles.get(n5.id).clone()));
    }
    if n6.exists() {
        quads_to_check.push((t1, n6.id));
    }

    #[cfg(feature = "debug_context")]
    debug_context.push_snapshot_event(
        Phase::SplitTriangle,
        EventInfo::Split(vertex_id),
        &triangles,
        &[t1, t2, t3],
        // Neighbor data is out of date here and is not that interesting
        &[],
    );
}

pub(crate) fn update_neighbor_neighbor(
    triangle: Neighbor,
    old_neighbour_id: Neighbor,
    new_neighbour_id: Neighbor,
    triangles: &mut Triangles,
) {
    #[cfg(feature = "more_profile_traces")]
    let _span = span!(Level::TRACE, "update_neighbor_neighbor").entered();

    if triangle.exists() {
        let t = triangles.get_mut(triangle.id);
        if t.neighbor12() == old_neighbour_id {
            *t.neighbor12_mut() = new_neighbour_id;
        } else if t.neighbor23() == old_neighbour_id {
            *t.neighbor23_mut() = new_neighbour_id;
        } else if t.neighbor31() == old_neighbour_id {
            *t.neighbor31_mut() = new_neighbour_id;
        }
    }
}

pub(crate) fn update_triangle_neighbor(
    triangle_id: TriangleId,
    old_neighbour_id: Neighbor,
    new_neighbour_id: Neighbor,
    triangles: &mut Triangles,
) {
    #[cfg(feature = "more_profile_traces")]
    let _span = span!(Level::TRACE, "update_triangle_neighbor").entered();

    let t = triangles.get_mut(triangle_id);
    if t.neighbor12() == old_neighbour_id {
        *t.neighbor12_mut() = new_neighbour_id;
    } else if t.neighbor23() == old_neighbour_id {
        *t.neighbor23_mut() = new_neighbour_id;
    } else if t.neighbor31() == old_neighbour_id {
        *t.neighbor31_mut() = new_neighbour_id;
    }
}

/// `quads_to_check` used the shared pre-allocated buffer for efficiency.
/// - It does not need to be cleared since it is fully emptied by each call to restore_delaunay_triangulation
fn restore_delaunay_triangulation(
    triangles: &mut Triangles,
    vertices: &Vec<Vertex>,
    from_vertex_id: VertexId,
    quads_to_check: &mut Vec<(TriangleId, TriangleId)>,
    // triangle_deltas: &mut Vec<(TriangleId, TriangleData)>,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "restore_delaunay_triangulation").entered();

    while let Some((from_triangle_id, opposite_triangle_id)) = quads_to_check.pop() {
        check_and_swap_quad_diagonal(
            triangles,
            vertices,
            from_vertex_id,
            from_triangle_id,
            opposite_triangle_id,
            quads_to_check,
            // triangle_deltas,
            #[cfg(feature = "debug_context")]
            debug_context,
        );
    }
}

#[inline(always)]
pub(crate) fn should_swap_diagonals(quad: &Quad, vertices: &Vec<Vertex>) -> bool {
    #[cfg(feature = "more_profile_traces")]
    let _span = span!(Level::TRACE, "should_swap_diagonals").entered();

    // Get rid of this case early: an infinite vertex cannot be in the circumcircle of the other triangle
    // Note: infinite q4 is not possible in DT (we always place a finite vertex in q4).
    if is_infinite(quad.v3()) {
        return false;
    }

    let is_q1_finite = is_finite(quad.v1());
    let is_q2_finite = is_finite(quad.v2());
    if is_q1_finite && is_q2_finite {
        is_vertex_in_triangle_circumcircle(
            vertices[quad.v1() as usize],
            vertices[quad.v2() as usize],
            vertices[quad.v3() as usize],
            vertices[quad.v4() as usize],
        )
    } else if !is_q1_finite && !is_q2_finite {
        true
    } else if !is_q1_finite {
        is_vertex_in_half_plane_1(vertices, quad, VERT_1)
    } else {
        is_vertex_in_half_plane_1(vertices, quad, VERT_2)
    }
}

/// ```text
///                q3
///         n3   /    \   n4
///            /   To   \
///          /            \
///         q1 ---------- q2
///          \ 2        3 /
///            \   Tf   /
///         n2   \ 1  /   n1
///                q4
/// ```
///
/// If q4 is in the circumcircle of the tirangle q1q2q3, the two triangles form a convex quadrilateral
/// whose diagonal is drawn in the wrong direction. We swap this diagonal to form two new triangles so
/// that the structure of the Delaunay triangulation is locally restored.
///
/// The quad becomes
///
/// ```text
///               q3
///        n3   / 3|2 \   n4
///           /    |    \
///         /      |      \
///        q1 2  Tf|To   3 q2
///         \      |      /
///           \    |    /
///        n2   \ 1|1 /   n1
///               q4
/// ```
pub(crate) fn check_and_swap_quad_diagonal(
    triangles: &mut Triangles,
    vertices: &Vec<Vertex>,
    from_vertex_id: VertexId,
    from_triangle_id: TriangleId,
    opposite_triangle_id: TriangleId,
    quads_to_check: &mut Vec<(TriangleId, TriangleId)>,
    // triangle_deltas: &mut Vec<(TriangleId, TriangleData)>,
    #[cfg(feature = "debug_context")] debug_context: &mut DebugContext,
) {
    #[cfg(feature = "more_profile_traces")]
    let _span = span!(Level::TRACE, "check_and_swap_quad_diagonal").entered();

    let opposite_triangle = triangles.get(opposite_triangle_id);

    let (quad, n3, n4) =
    // No need to check if neighbor exists, handled by the == check since `from_triangle_id` exists
        if opposite_triangle.neighbor12().id == from_triangle_id {
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
        } else if opposite_triangle.neighbor23().id == from_triangle_id {
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

    if should_swap_diagonals(&quad, vertices) {
        let opposite_neighbor = opposite_triangle_id.into();
        let from_neighbor = from_triangle_id.into();
        let n1 = triangles.get(from_triangle_id).neighbor31();

        triangles.get_mut(from_triangle_id).verts = [quad.v4(), quad.v1(), quad.v3()];
        triangles.get_mut(opposite_triangle_id).verts = [quad.v4(), quad.v3(), quad.v2()];

        triangles.get_mut(opposite_triangle_id).neighbors = [from_neighbor, n4, n1];
        *triangles.get_mut(from_triangle_id).neighbor23_mut() = n3;
        *triangles.get_mut(from_triangle_id).neighbor31_mut() = opposite_neighbor;

        // triangle_deltas.push((from_triangle_id, triangles.get(from_triangle_id).clone()));
        // triangle_deltas.push((
        //     opposite_triangle_id,
        //     triangles.get(opposite_triangle_id).clone(),
        // ));

        // Update neighbors and place any new triangles pairs which are now opposite to `from_vertex_id` on the stack, to be checked
        // TODO Optim: could organize it better to regroup n1 and n4 ?
        if n1.exists() {
            update_triangle_neighbor(n1.id, from_neighbor, opposite_neighbor, triangles);
            // triangle_deltas.push((n1.id, triangles.get(n1.id).clone()));
        }
        if n3.exists() {
            update_triangle_neighbor(n3.id, opposite_neighbor, from_neighbor, triangles);
            quads_to_check.push((from_triangle_id, n3.id));
            // triangle_deltas.push((n3.id, triangles.get(n3.id).clone()));
        }
        if n4.exists() {
            quads_to_check.push((opposite_triangle_id, n4.id));
        }

        #[cfg(feature = "debug_context")]
        debug_context.push_snapshot(
            Phase::DelaunayRestoreSwapQuadDiagonals,
            &triangles,
            &[from_triangle_id, opposite_triangle_id],
            &[n3, triangles.get(from_triangle_id).neighbor31()],
        );
    }
}

///////////////////////////////////////////////////////////
///                                                     ///
///                        Tests                        ///
///                                                     ///
///////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use glam::Vec3;

    use crate::{
        triangulation::{normalize_vertices_coordinates, transform_to_2d_planar_coordinate_system},
        types::Vertex,
    };

    #[test]
    fn triangulation_normalize_set_of_vertices() {
        let vertices = vec![
            Vertex::new(3., 2.),
            Vertex::new(-1., 2.),
            Vertex::new(-1., -2.),
            Vertex::new(3., -2.),
        ];

        let normalized_vertices = normalize_vertices_coordinates(&vertices);

        assert_eq!(
            vec![
                Vertex::new(1., 1.),
                Vertex::new(0., 1.),
                Vertex::new(0., 0.),
                Vertex::new(1., 0.),
            ],
            normalized_vertices.0
        );
    }

    #[test]
    fn triangulation_set_to_2d_plane_vertices() {
        let vertices = vec![
            Vec3::new(-3., 2., 0.),
            Vec3::new(1., 2., 0.),
            Vec3::new(11., -2., 0.),
            Vec3::new(-3., -2., 0.),
        ];

        let planar_vertices = transform_to_2d_planar_coordinate_system(&vertices, Vec3::Z);

        assert_eq!(
            vec![
                -Vertex::new(-3., 2.),
                -Vertex::new(1., 2.),
                -Vertex::new(11., -2.),
                -Vertex::new(-3., -2.),
            ],
            planar_vertices
        );
    }
}
