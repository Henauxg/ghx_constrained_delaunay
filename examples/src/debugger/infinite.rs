use bevy::math::Vec3;
use ghx_constrained_delaunay::{
    infinite::{
        infinite_vertex_local_quad_index, INFINITE_VERTS_X_DELTAS, INFINITE_VERTS_Y_DELTAS,
    },
    types::{
        next_cw_vertex_index, next_ccw_vertex_index, Float, TriangleData,
        TriangleVertexIndex, NEXT_CCW_VERTEX_INDEX, NEXT_CW_VERTEX_INDEX,
    },
};

use super::Basis;

// TODO Change to a factor dependent on the triangle/triangulation
pub const SEMI_INFINITE_EDGE_VISUAL_LEN_FACTOR: Float = 50.;

#[inline]
pub fn extrapolated_vertex_from_semi_infinite_edge(
    finite_vertex: Vec3,
    infinite_vert_local_index: TriangleVertexIndex,
    basis: &Basis,
) -> Vec3 {
    finite_vertex
        + (SEMI_INFINITE_EDGE_VISUAL_LEN_FACTOR
            * INFINITE_VERTS_X_DELTAS[infinite_vert_local_index as usize]) as f32
            * basis.e1
        + (SEMI_INFINITE_EDGE_VISUAL_LEN_FACTOR
            * INFINITE_VERTS_Y_DELTAS[infinite_vert_local_index as usize]) as f32
            * basis.e2
}

pub fn get_2_extrapolated_vertices_from_2_finite_vertices(
    vertices: &[Vec3],
    triangle: &TriangleData,
    basis: &Basis,
    infinite_vert_index: TriangleVertexIndex,
) -> (Vec3, Vec3, Vec3, Vec3) {
    let infinite_vert_local_index =
        infinite_vertex_local_quad_index(triangle.v(infinite_vert_index));

    let finite_vert_a_index = NEXT_CW_VERTEX_INDEX[infinite_vert_index as usize];
    let finite_vert_a_id = triangle.v(finite_vert_a_index);
    let finite_vertex_a = vertices[finite_vert_a_id as usize];
    let finite_tip_a = extrapolated_vertex_from_semi_infinite_edge(
        finite_vertex_a,
        infinite_vert_local_index,
        basis,
    );

    let finite_vert_b_index = NEXT_CCW_VERTEX_INDEX[infinite_vert_index as usize];
    let finite_vert_b_id = triangle.v(finite_vert_b_index);
    let finite_vertex_b = vertices[finite_vert_b_id as usize];
    let finite_tip_b = extrapolated_vertex_from_semi_infinite_edge(
        finite_vertex_b,
        infinite_vert_local_index,
        basis,
    );
    (finite_vertex_a, finite_tip_a, finite_vertex_b, finite_tip_b)
}

pub fn get_2_extrapolated_vertices_from_1_finite_vertex(
    vertices: &[Vec3],
    triangle: &TriangleData,
    basis: &Basis,
    finite_vertex_index: TriangleVertexIndex,
) -> (Vec3, Vec3, Vec3) {
    let finite_vert_id = triangle.v(finite_vertex_index);
    let finite_vertex = vertices[finite_vert_id as usize];

    let infinite_vert_a_local_index = infinite_vertex_local_quad_index(
        triangle.v(next_cw_vertex_index(finite_vertex_index)),
    );
    let finite_tip_a = extrapolated_vertex_from_semi_infinite_edge(
        finite_vertex,
        infinite_vert_a_local_index,
        basis,
    );

    let infinite_vert_b_local_index = infinite_vertex_local_quad_index(
        triangle.v(next_ccw_vertex_index(finite_vertex_index)),
    );
    let finite_tip_b = extrapolated_vertex_from_semi_infinite_edge(
        finite_vertex,
        infinite_vert_b_local_index,
        basis,
    );

    (finite_tip_a, finite_vertex, finite_tip_b)
}
