use bevy::{
    color::{palettes::css::RED, Color},
    ecs::system::Res,
    gizmos::gizmos::Gizmos,
    math::Dir3,
};
use ghx_constrained_delaunay::glam::Quat;
use ghx_constrained_delaunay::{
    debug::EventInfo,
    infinite::infinite_vertex_local_quad_index,
    types::{
        next_clockwise_vertex_index, next_counter_clockwise_vertex_index, TriangleData, TriangleId,
        Vector3, NEXT_CCW_VERTEX_INDEX, NEXT_CW_VERTEX_INDEX, VERT_1, VERT_2, VERT_3,
    },
};

use crate::{finite_vertex_from_semi_infinite_edge, COLORS};

use super::{Basis, TrianglesDebugData, TrianglesDebugViewConfig, TrianglesDrawMode};

pub fn draw_triangles_debug_data_gizmos(
    mut gizmos: Gizmos,
    debug_data: Res<TrianglesDebugData>,
    view_config: Res<TrianglesDebugViewConfig>,
) {
    let Some(snapshot) = debug_data.current_snapshot() else {
        return;
    };
    let vertices = &debug_data.vertices;

    // Draw triangles
    match view_config.triangles_draw_mode {
        TrianglesDrawMode::AllAsGizmos => {
            for (triangle_id, triangle) in snapshot.triangles.buffer().iter().enumerate() {
                draw_triangle_gizmo(
                    &mut gizmos,
                    vertices,
                    &debug_data.basis,
                    triangle_id as TriangleId,
                    triangle,
                );
            }
        }
        TrianglesDrawMode::ChangedAsGizmos => {
            for triangle_id in snapshot.changed_ids.iter() {
                let triangle = &snapshot.triangles.get(*triangle_id);
                draw_triangle_gizmo(
                    &mut gizmos,
                    vertices,
                    &debug_data.basis,
                    *triangle_id,
                    triangle,
                );
            }
        }
        _ => (),
    }

    // Draw last changes boundaries
    gizmos.rect(
        debug_data.current_changes_bounds.0,
        Quat::IDENTITY,
        debug_data.current_changes_bounds.1,
        Color::WHITE,
    );

    // Draw specific event info
    if let EventInfo::Split(vertex_id) = snapshot.event {
        // Set circle radius as proportional to the changes bounding box
        let circle_radius = debug_data
            .current_changes_bounds
            .1
            .x
            .min(debug_data.current_changes_bounds.1.y)
            / 10.;
        gizmos.circle(
            vertices[vertex_id as usize].as_vec3(),
            Dir3::Z,
            circle_radius,
            RED,
        );
    }
}

pub fn draw_triangle_gizmo(
    gizmos: &mut Gizmos,
    vertices: &Vec<Vector3>,
    _basis: &Option<Basis>,
    triangle_id: TriangleId,
    triangle: &TriangleData,
) {
    let color = COLORS[triangle_id as usize % COLORS.len()];

    let mut infinite_verts = Vec::new();
    if triangle.v1() as usize >= vertices.len() {
        infinite_verts.push(VERT_1);
    }
    if triangle.v2() as usize >= vertices.len() {
        infinite_verts.push(VERT_2);
    }
    if triangle.v3() as usize >= vertices.len() {
        infinite_verts.push(VERT_3);
    }

    let linestrip = if infinite_verts.is_empty() {
        let (v1, v2, v3) = (
            vertices[triangle.v1() as usize].as_vec3(),
            vertices[triangle.v2() as usize].as_vec3(),
            vertices[triangle.v3() as usize].as_vec3(),
        );
        vec![v1, v2, v3, v1]
    } else if infinite_verts.len() == 1 {
        let infinite_vert_local_index =
            infinite_vertex_local_quad_index(triangle.v(infinite_verts[0]));

        let finite_vert_a_index = NEXT_CW_VERTEX_INDEX[infinite_verts[0] as usize];
        let finite_vert_a_id = triangle.v(finite_vert_a_index);
        let finite_vertex_a = vertices[finite_vert_a_id as usize];
        let finite_tip_a =
            finite_vertex_from_semi_infinite_edge(finite_vertex_a, infinite_vert_local_index);

        let finite_vert_b_index = NEXT_CCW_VERTEX_INDEX[infinite_verts[0] as usize];
        let finite_vert_b_id = triangle.v(finite_vert_b_index);
        let finite_vertex_b = vertices[finite_vert_b_id as usize];
        let finite_tip_b =
            finite_vertex_from_semi_infinite_edge(finite_vertex_b, infinite_vert_local_index);

        // if let Some(basis) = basis {
        //     finite_tip_a = detransfrom(finite_tip_a, basis);
        //     finite_tip_b = detransfrom(finite_tip_b, basis);
        // }
        // info!("draw gizmos, infinite_verts.len() == 1");

        vec![
            finite_tip_a,
            finite_vertex_a.as_vec3(),
            finite_vertex_b.as_vec3(),
            finite_tip_b,
        ]
    } else if infinite_verts.len() == 2 {
        //TODO
        let finite_v_index = 3 - (infinite_verts[0] + infinite_verts[1]);
        let finite_vert_id = triangle.v(finite_v_index);
        let finite_vertex = vertices[finite_vert_id as usize];

        let infinite_vert_a_local_index = infinite_vertex_local_quad_index(
            triangle.v(next_clockwise_vertex_index(finite_v_index)),
        );
        let finite_tip_a =
            finite_vertex_from_semi_infinite_edge(finite_vertex, infinite_vert_a_local_index);

        let infinite_vert_b_local_index = infinite_vertex_local_quad_index(
            triangle.v(next_counter_clockwise_vertex_index(finite_v_index)),
        );
        let finite_tip_b =
            finite_vertex_from_semi_infinite_edge(finite_vertex, infinite_vert_b_local_index);

        // if let Some(basis) = basis {
        //     finite_tip_a = detransfrom(finite_tip_a, basis);
        //     finite_tip_b = detransfrom(finite_tip_b, basis);
        // }

        vec![finite_tip_a, finite_vertex.as_vec3(), finite_tip_b]
    } else {
        vec![]
    };
    gizmos.linestrip(linestrip, color);
    // TODO Could draw something for 3 infinite vertices
}
