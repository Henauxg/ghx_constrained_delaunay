use crate::{
    triangulation::{is_vertex_pair_too_close, VertexPlacement, DELTA_VALUE},
    types::{
        next_clockwise_vertex_index, next_counter_clockwise_vertex_index, opposite_edge_index,
        opposite_vertex_index_from_edge, vertex_next_ccw_edge_index, vertex_next_cw_edge_index,
        EdgeVertices, Float, Neighbor, Quad, QuadVertexIndex, TriangleData, TriangleId,
        TriangleVertexIndex, Vertex, VertexId, EDGE_TO_VERTS, NEXT_CCW_VERTEX_INDEX,
        NEXT_CW_VERTEX_INDEX, QUAD_4, VERT_1, VERT_2, VERT_3,
    },
    utils::{is_point_strictly_on_right_side_of_edge, test_point_edge_side},
};

use arrayvec::ArrayVec;
#[cfg(feature = "profile_traces")]
use tracing::{span, Level};

pub const INFINITE_V0_ID: VertexId = VertexId::MAX - 3;
pub const INFINITE_V1_ID: VertexId = VertexId::MAX - 2;
pub const INFINITE_V2_ID: VertexId = VertexId::MAX - 1;
pub const INFINITE_V3_ID: VertexId = VertexId::MAX;

#[inline]
/// For an infinite vertex id, returns its local vertex index in the infinite quad
///
/// INVALID for a finite vertex
pub fn infinite_vertex_local_quad_index(vertex_id: VertexId) -> QuadVertexIndex {
    (vertex_id - INFINITE_V0_ID) as QuadVertexIndex
}

#[inline]
/// From a vertex id, returns true if it represents an infinite vertex
pub fn is_infinite(vert_id: VertexId) -> bool {
    vert_id >= INFINITE_V0_ID
}

#[inline]
/// From a vertex id, returns true if it represents an infinite vertex
pub fn is_finite(vert_id: VertexId) -> bool {
    vert_id < INFINITE_V0_ID
}

pub const INFINITE_VERTS_Y_DELTAS: [Float; 4] = [0., DELTA_VALUE, 0., -DELTA_VALUE];
pub const INFINITE_VERTS_X_DELTAS: [Float; 4] = [-DELTA_VALUE, 0., DELTA_VALUE, 0.];

/// Returns a finite segment from an edge between a finite vertex and an infinite vertex
///  - infinite_vert_local_index is the local index of the infinite vertex
#[inline]
pub fn edge_from_semi_infinite_edge(
    finite_vertex: Vertex,
    infinite_vert_local_index: QuadVertexIndex,
) -> EdgeVertices {
    (
        finite_vertex,
        Vertex::new(
            // We care about the delta sign here since we create a segment from an infinite line,
            // starting from the finite point and aimed towards the infinite point.
            finite_vertex.x + INFINITE_VERTS_X_DELTAS[infinite_vert_local_index as usize],
            finite_vertex.y + INFINITE_VERTS_Y_DELTAS[infinite_vert_local_index as usize],
        ),
    )
}

#[inline(always)]
/// The slice MUST have at least 3 vertices
pub(crate) fn collect_infinite_triangle_vertices(
    tri: &[VertexId],
) -> ArrayVec<TriangleVertexIndex, 2> {
    #[cfg(feature = "more_profile_traces")]
    let _span = span!(Level::TRACE, "collect_infinite_triangle_vertices").entered();

    let mut infinite_verts = ArrayVec::new();
    if is_infinite(tri[VERT_1 as usize]) {
        infinite_verts.push(VERT_1);
    }
    if is_infinite(tri[VERT_2 as usize]) {
        infinite_verts.push(VERT_2);
    }
    if is_infinite(tri[VERT_3 as usize]) {
        infinite_verts.push(VERT_3);
    }
    infinite_verts
}

#[cold]
pub(crate) fn is_vertex_in_half_plane_1(
    vertices: &Vec<Vertex>,
    quad: &Quad,
    infinite_vert: TriangleVertexIndex,
) -> bool {
    #[cfg(feature = "more_profile_traces")]
    let _span = span!(Level::TRACE, "is_vertex_in_half_plane_1").entered();

    // Test if q4 is inside the circle with 1 infinite point (half-plane defined by the 2 finite points)
    // TODO Clean: utils functions in Quad/Triangle
    let edge_index = opposite_edge_index(infinite_vert);
    let finite_vert_indexes = EDGE_TO_VERTS[edge_index as usize];
    // q1q2q3 is in a CCW order, so we reverse the edge
    let edge_vertices = (
        vertices[quad.v(finite_vert_indexes[1]) as usize],
        vertices[quad.v(finite_vert_indexes[0]) as usize],
    );

    // Strict seems to be necessary here, in order to avoid creating flat triangles
    is_point_strictly_on_right_side_of_edge(edge_vertices, vertices[quad.v4() as usize])
}

/// Test if q4 is inside the circle with 2 infinite points (half-plane defined by the finite point and the slope between the 2 infinite points)
#[cold]
pub(crate) fn is_vertex_in_half_plane_2(
    vertices: &Vec<Vertex>,
    quad: &Quad,
    infinite_v1: TriangleVertexIndex,
    infinite_v2: TriangleVertexIndex,
) -> bool {
    #[cfg(feature = "more_profile_traces")]
    let _span = span!(Level::TRACE, "is_vertex_in_half_plane_2").entered();

    const INFINITE_VERTS_SLOPES: [Float; 3] = [1., -1., 1.];

    let infinite_vert_1_local_index = infinite_vertex_local_quad_index(quad.v(infinite_v1));
    let infinite_vert_2_local_index = infinite_vertex_local_quad_index(quad.v(infinite_v2));

    // Index of the finite vertex in q1q2q3
    let finite_vert_index = opposite_vertex_index_from_edge(infinite_v1, infinite_v2);

    let line_point_1 = vertices[quad.v(finite_vert_index) as usize];
    // Sum/2 defines a convenient surjection to (0..2)
    let a = INFINITE_VERTS_SLOPES
        [((infinite_vert_1_local_index + infinite_vert_2_local_index) / 2) as usize];
    let b = line_point_1.y - a * line_point_1.x;
    // Another point on the line, its position depends of which infinite vertices are considered. We need delta < 0 for edges Q1-Q2 and Q2-Q3, positive otherwise.
    let delta_x = if infinite_vert_1_local_index == QUAD_4 || infinite_vert_2_local_index == QUAD_4
    {
        DELTA_VALUE
    } else {
        -DELTA_VALUE
    };
    let line_point_2_x = line_point_1.x + delta_x;
    let line_point_2 = Vertex::new(line_point_2_x, a * line_point_2_x + b);
    // Strict seems to be necessary here, in order to avoid creating flat triangles
    is_point_strictly_on_right_side_of_edge(
        (line_point_1, line_point_2),
        vertices[quad.v4() as usize],
    )
}

#[cold]
pub(crate) fn vertex_placement_1_infinite_vertex(
    vertices: &Vec<Vertex>,
    vertex: Vertex,
    triangle: &TriangleData,
    triangle_id: TriangleId,
    infinite_verts: ArrayVec<TriangleVertexIndex, 2>,
    previous_triangle: &mut Neighbor,
    current_triangle: &mut Neighbor,
) -> Option<VertexPlacement> {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "vertex_placement_1_infinite_vertex").entered();

    let finite_vert_a_index = NEXT_CW_VERTEX_INDEX[infinite_verts[0] as usize];
    let finite_vert_a_id = triangle.v(finite_vert_a_index);
    let finite_vertex_a = vertices[finite_vert_a_id as usize];

    if is_vertex_pair_too_close(finite_vertex_a, vertex) {
        return Some(VertexPlacement::OnVertex(finite_vert_a_id));
    }

    let finite_vert_b_index = NEXT_CCW_VERTEX_INDEX[infinite_verts[0] as usize];
    let finite_vert_b_id = triangle.v(finite_vert_b_index);
    let finite_vertex_b = vertices[finite_vert_b_id as usize];

    if is_vertex_pair_too_close(finite_vertex_b, vertex) {
        return Some(VertexPlacement::OnVertex(finite_vert_b_id));
    }

    // Check if the point is towards this triangle's neighbor
    // Check previous triangle to avoid looping between two triangles
    let edge_ab_test = test_point_edge_side((finite_vertex_a, finite_vertex_b), vertex);
    let edge_ab_index = vertex_next_cw_edge_index(finite_vert_a_index);
    if edge_ab_test.is_on_left_side()
           // Avoid looping between two triangles
        && *previous_triangle != triangle.neighbor(edge_ab_index)
    {
        *previous_triangle = *current_triangle;
        *current_triangle = triangle.neighbor(edge_ab_index);
        return None;
    }

    let infinite_vert_local_index = infinite_vertex_local_quad_index(triangle.v(infinite_verts[0]));
    let edge_a = edge_from_semi_infinite_edge(finite_vertex_a, infinite_vert_local_index);
    let edge_a_index = vertex_next_ccw_edge_index(finite_vert_a_index);
    let edge_a_test = test_point_edge_side(edge_a, vertex);
    // Test for right side because edge_a is oriented the wrong way. We could also reverse edge_a
    if edge_a_test.is_on_right_side() && *previous_triangle != triangle.neighbor(edge_a_index) {
        *previous_triangle = *current_triangle;
        *current_triangle = triangle.neighbor(edge_a_index);
        return None;
    }

    let edge_b = edge_from_semi_infinite_edge(finite_vertex_b, infinite_vert_local_index);
    let edge_b_index = vertex_next_cw_edge_index(finite_vert_b_index);
    let edge_b_test = test_point_edge_side(edge_b, vertex);
    if edge_b_test.is_on_left_side() && *previous_triangle != triangle.neighbor(edge_b_index) {
        *previous_triangle = *current_triangle;
        *current_triangle = triangle.neighbor(edge_b_index);
        return None;
    }

    if edge_ab_test.is_near_edge() {
        return Some(VertexPlacement::OnTriangleEdge(triangle_id, edge_ab_index));
    } else if edge_a_test.is_near_edge() {
        return Some(VertexPlacement::OnTriangleEdge(triangle_id, edge_a_index));
    } else if edge_b_test.is_near_edge() {
        return Some(VertexPlacement::OnTriangleEdge(triangle_id, edge_b_index));
    } else {
        return Some(VertexPlacement::InsideTriangle(triangle_id));
    }
}

#[cold]
pub(crate) fn vertex_placement_2_infinite_vertex(
    vertices: &Vec<Vertex>,
    vertex: Vertex,
    triangle: &TriangleData,
    triangle_id: TriangleId,
    infinite_verts: ArrayVec<TriangleVertexIndex, 2>,
    previous_triangle: &mut Neighbor,
    current_triangle: &mut Neighbor,
) -> Option<VertexPlacement> {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "vertex_placement_2_infinite_vertex").entered();

    let finite_v_index = opposite_vertex_index_from_edge(infinite_verts[0], infinite_verts[1]);
    let finite_vert_id = triangle.v(finite_v_index);
    let finite_vertex = vertices[finite_vert_id as usize];

    // No need to check if vertex is too close to an infinite vertex here
    if is_vertex_pair_too_close(finite_vertex, vertex) {
        return Some(VertexPlacement::OnVertex(finite_vert_id));
    }

    // Check if the point is towards this triangle's neighbor
    // Check previous triangle to avoid looping between two triangles
    let infinite_vert_a_local_index =
        infinite_vertex_local_quad_index(triangle.v(next_clockwise_vertex_index(finite_v_index)));
    let edge_a = edge_from_semi_infinite_edge(finite_vertex, infinite_vert_a_local_index);
    let edge_a_index = vertex_next_cw_edge_index(finite_v_index);
    let edge_a_test = test_point_edge_side(edge_a, vertex);
    if edge_a_test.is_on_left_side() && *previous_triangle != triangle.neighbor(edge_a_index) {
        *previous_triangle = *current_triangle;
        *current_triangle = triangle.neighbor(edge_a_index);
        return None;
    }

    let infinite_vert_b_local_index = infinite_vertex_local_quad_index(
        triangle.v(next_counter_clockwise_vertex_index(finite_v_index)),
    );
    let edge_b = edge_from_semi_infinite_edge(finite_vertex, infinite_vert_b_local_index);
    let edge_b_index = vertex_next_ccw_edge_index(finite_v_index);
    let edge_b_test = test_point_edge_side(edge_b, vertex);
    // Test for right side because edge_b is oriented the wrong way. We could also reverse edge_b
    if edge_b_test.is_on_right_side() && *previous_triangle != triangle.neighbor(edge_b_index) {
        *previous_triangle = *current_triangle;
        *current_triangle = triangle.neighbor(edge_b_index);
        return None;
    }

    if edge_a_test.is_near_edge() {
        return Some(VertexPlacement::OnTriangleEdge(triangle_id, edge_a_index));
    } else if edge_b_test.is_near_edge() {
        return Some(VertexPlacement::OnTriangleEdge(triangle_id, edge_b_index));
    } else {
        return Some(VertexPlacement::InsideTriangle(triangle_id));
    }
}
