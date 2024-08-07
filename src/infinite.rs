use crate::{
    triangulation::VertexPlacement,
    types::{
        next_ccw_vertex_index, next_cw_vertex_index, opposite_edge_index,
        opposite_vertex_index_from_edge, vertex_next_ccw_edge_index, vertex_next_cw_edge_index,
        Edge, EdgeVertices, Float, Neighbor, Quad, QuadVertexIndex, TriangleData, TriangleId,
        TriangleVertexIndex, Vertex, VertexId, ADJACENT_QUAD_VERTICES_INDEXES, EDGE_TO_VERTS,
        NEXT_CCW_VERTEX_INDEX, NEXT_CW_VERTEX_INDEX, OPPOSITE_QUAD_VERTEX_INDEX, QUAD_1, QUAD_2,
        QUAD_3, QUAD_4, VERT_1, VERT_2, VERT_3,
    },
    utils::{
        egdes_intersect, is_point_strictly_on_right_side_of_edge, test_point_edge_side,
        EdgesIntersectionResult,
    },
};

use arrayvec::ArrayVec;
#[cfg(feature = "profile_traces")]
use tracing::{span, Level};

pub const INFINITE_V0_ID: VertexId = VertexId::MAX - 3;
pub const INFINITE_V1_ID: VertexId = VertexId::MAX - 2;
pub const INFINITE_V2_ID: VertexId = VertexId::MAX - 1;
pub const INFINITE_V3_ID: VertexId = VertexId::MAX;

// Slightly higher values than 1.0 to be safe, which is enough since all our vertices coordinates are normalized inside the unit square.
pub const EXTRAPOLATION_DELTA_VALUE: Float = 1.42;

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

pub const INFINITE_VERTS_Y_DELTAS: [Float; 4] = [
    0.,
    EXTRAPOLATION_DELTA_VALUE,
    0.,
    -EXTRAPOLATION_DELTA_VALUE,
];
pub const INFINITE_VERTS_X_DELTAS: [Float; 4] = [
    -EXTRAPOLATION_DELTA_VALUE,
    0.,
    EXTRAPOLATION_DELTA_VALUE,
    0.,
];

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
pub fn collect_infinite_triangle_vertices(tri: &[VertexId]) -> ArrayVec<TriangleVertexIndex, 2> {
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

#[inline(always)]
/// The slice MUST have at least 3 vertices
pub fn collect_infinite_quad_vertices(quad: &[VertexId]) -> ArrayVec<QuadVertexIndex, 2> {
    #[cfg(feature = "more_profile_traces")]
    let _span = span!(Level::TRACE, "collect_infinite_quad_vertices").entered();

    let mut infinite_verts = ArrayVec::new();
    if is_infinite(quad[QUAD_1 as usize]) {
        infinite_verts.push(QUAD_1);
    }
    if is_infinite(quad[QUAD_2 as usize]) {
        infinite_verts.push(QUAD_2);
    }
    if is_infinite(quad[QUAD_3 as usize]) {
        infinite_verts.push(QUAD_3);
    }
    if is_infinite(quad[QUAD_4 as usize]) {
        infinite_verts.push(QUAD_3);
    }
    infinite_verts
}

#[cold]
pub(crate) fn is_vertex_in_half_plane_1(
    vertices: &[Vertex],
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

#[cold]
pub(crate) fn vertex_placement_1_infinite_vertex(
    vertices: &[Vertex],
    vertex: Vertex,
    triangle: &TriangleData,
    triangle_id: TriangleId,
    infinite_vertex_index: TriangleVertexIndex,
    previous_triangle: &mut Neighbor,
    current_triangle: &mut Neighbor,
) -> Option<VertexPlacement> {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "vertex_placement_1_infinite_vertex").entered();

    let finite_vert_a_index = NEXT_CW_VERTEX_INDEX[infinite_vertex_index as usize];
    let finite_vert_a_id = triangle.v(finite_vert_a_index);
    let finite_vertex_a = vertices[finite_vert_a_id as usize];

    let finite_vert_b_index = NEXT_CCW_VERTEX_INDEX[infinite_vertex_index as usize];
    let finite_vert_b_id = triangle.v(finite_vert_b_index);
    let finite_vertex_b = vertices[finite_vert_b_id as usize];

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

    let infinite_vertex_id = triangle.v(infinite_vertex_index);
    let infinite_vert_local_index = infinite_vertex_local_quad_index(infinite_vertex_id);
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
        return Some(VertexPlacement::OnTriangleEdge(
            triangle_id,
            edge_ab_index,
            Edge::new(finite_vert_a_id, finite_vert_b_id),
        ));
    } else if edge_a_test.is_near_edge() {
        return Some(VertexPlacement::OnTriangleEdge(
            triangle_id,
            edge_a_index,
            Edge::new(infinite_vertex_id, finite_vert_a_id),
        ));
    } else if edge_b_test.is_near_edge() {
        return Some(VertexPlacement::OnTriangleEdge(
            triangle_id,
            edge_b_index,
            Edge::new(finite_vert_b_id, infinite_vertex_id),
        ));
    } else {
        return Some(VertexPlacement::InsideTriangle(triangle_id));
    }
}

#[cold]
pub(crate) fn vertex_placement_2_infinite_vertex(
    vertices: &[Vertex],
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

    // Check if the point is towards this triangle's neighbors

    let infinite_vert_a_id = triangle.v(next_cw_vertex_index(finite_v_index));
    let infinite_vert_a_local_index = infinite_vertex_local_quad_index(infinite_vert_a_id);
    let edge_a = edge_from_semi_infinite_edge(finite_vertex, infinite_vert_a_local_index);
    let edge_a_index = vertex_next_cw_edge_index(finite_v_index);
    let edge_a_test = test_point_edge_side(edge_a, vertex);
    // Check previous triangle to avoid looping between two triangles
    if edge_a_test.is_on_left_side() && *previous_triangle != triangle.neighbor(edge_a_index) {
        *previous_triangle = *current_triangle;
        *current_triangle = triangle.neighbor(edge_a_index);
        return None;
    }

    let infinite_vert_b_id = triangle.v(next_ccw_vertex_index(finite_v_index));
    let infinite_vert_b_local_index = infinite_vertex_local_quad_index(infinite_vert_b_id);
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
        return Some(VertexPlacement::OnTriangleEdge(
            triangle_id,
            edge_a_index,
            Edge::new(finite_vert_id, infinite_vert_a_id),
        ));
    } else if edge_b_test.is_near_edge() {
        return Some(VertexPlacement::OnTriangleEdge(
            triangle_id,
            edge_b_index,
            Edge::new(infinite_vert_b_id, finite_vert_id),
        ));
    } else {
        return Some(VertexPlacement::InsideTriangle(triangle_id));
    }
}

#[cold]
pub(crate) fn quad_diagonals_intersection_1_infinite(
    vertices: &[Vertex],
    quad: &Quad,
    infinite_vertex_index: QuadVertexIndex,
) -> EdgesIntersectionResult {
    let infinite_vert_local_index = infinite_vertex_local_quad_index(quad.v(infinite_vertex_index));
    let finite_vertex =
        vertices[quad.v(OPPOSITE_QUAD_VERTEX_INDEX[infinite_vertex_index as usize]) as usize];
    let other_diagonal_verts = ADJACENT_QUAD_VERTICES_INDEXES[infinite_vertex_index as usize];
    let other_diagonal_edge = (
        vertices[quad.v(other_diagonal_verts[0]) as usize],
        vertices[quad.v(other_diagonal_verts[1]) as usize],
    );

    egdes_intersect(
        &other_diagonal_edge,
        &edge_from_semi_infinite_edge(finite_vertex, infinite_vert_local_index),
    )
}

#[cold]
pub(crate) fn quad_diagonals_intersection_2_infinite(
    vertices: &[Vertex],
    quad: &Quad,
    infinite_vert_index_0: QuadVertexIndex,
    infinite_vert_index_1: QuadVertexIndex,
) -> EdgesIntersectionResult {
    let infinite_vert_local_index_0 =
        infinite_vertex_local_quad_index(quad.v(infinite_vert_index_0));
    let finite_vert_0 =
        vertices[quad.v(OPPOSITE_QUAD_VERTEX_INDEX[infinite_vert_index_0 as usize]) as usize];
    let infinite_edge_segment_0 =
        edge_from_semi_infinite_edge(finite_vert_0, infinite_vert_local_index_0);

    let infinite_vert_local_index_1 =
        infinite_vertex_local_quad_index(quad.v(infinite_vert_index_1));
    let finite_vert_1 =
        vertices[quad.v(OPPOSITE_QUAD_VERTEX_INDEX[infinite_vert_index_1 as usize]) as usize];
    let infinite_edge_segment_1 =
        edge_from_semi_infinite_edge(finite_vert_1, infinite_vert_local_index_1);

    egdes_intersect(&infinite_edge_segment_0, &infinite_edge_segment_1)
}
