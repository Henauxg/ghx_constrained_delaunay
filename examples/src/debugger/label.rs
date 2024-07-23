use bevy::{
    ecs::system::Commands,
    log::error,
    prelude::default,
    text::{Text, TextSection, TextStyle},
    transform::components::Transform,
};
use bevy_mod_billboard::BillboardTextBundle;
use ghx_constrained_delaunay::glam::Vec3;
use ghx_constrained_delaunay::{
    infinite::infinite_vertex_local_quad_index,
    types::{
        next_clockwise_vertex_index, next_counter_clockwise_vertex_index,
        opposite_vertex_index_from_edge, TriangleData, TriangleId, Vector3, NEXT_CCW_VERTEX_INDEX,
        NEXT_CW_VERTEX_INDEX, VERT_1, VERT_2, VERT_3,
    },
};

use crate::COLORS;

pub const DEBUG_LABEL_FONT_SIZE: f32 = 60.0;

pub const BILLBOARD_DEFAULT_SCALE: Vec3 = Vec3::splat(0.0005);

use super::{finite_vertex_from_semi_infinite_edge, Basis, TriangleDebugEntity, VertexLabelMode};

pub fn spawn_label(
    commands: &mut Commands,
    vertex_label_mode: VertexLabelMode,
    vertices: &Vec<Vector3>,
    _basis: &Option<Basis>,
    triangle_id: TriangleId,
    triangle: &TriangleData,
) {
    let color = COLORS[(triangle_id as usize) % COLORS.len()];

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

    let (v1, v2, v3) = if infinite_verts.is_empty() {
        (
            vertices[triangle.v1() as usize].as_vec3(),
            vertices[triangle.v2() as usize].as_vec3(),
            vertices[triangle.v3() as usize].as_vec3(),
        )
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

        if infinite_verts[0] == VERT_1 {
            (
                (finite_tip_a + finite_tip_b) / 2.,
                vertices[triangle.v2() as usize].as_vec3(),
                vertices[triangle.v3() as usize].as_vec3(),
            )
        } else if infinite_verts[0] == VERT_2 {
            (
                vertices[triangle.v1() as usize].as_vec3(),
                (finite_tip_a + finite_tip_b) / 2.,
                vertices[triangle.v3() as usize].as_vec3(),
            )
        } else {
            (
                vertices[triangle.v1() as usize].as_vec3(),
                vertices[triangle.v2() as usize].as_vec3(),
                (finite_tip_a + finite_tip_b) / 2.,
            )
        }
    } else if infinite_verts.len() == 2 {
        let finite_v_index = opposite_vertex_index_from_edge(infinite_verts[0], infinite_verts[1]);

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

        if finite_v_index == VERT_1 {
            (finite_vertex.as_vec3(), finite_tip_a, finite_tip_b)
        } else if finite_v_index == VERT_2 {
            (finite_tip_b, finite_vertex.as_vec3(), finite_tip_a)
        } else {
            (finite_tip_a, finite_tip_b, finite_vertex.as_vec3())
        }
    } else {
        error!("Handle 2 infinite vertices for labels");
        return;
    };

    let center = (v1 + v2 + v3) / 3.;

    let v1v2 = v2 - v1;
    let v1v3 = v3 - v1;
    let v2v3 = v3 - v2;
    // We could also use triangle_area = 0.5 * v1v2.cross(v1v3).length();
    let max_side_length = v1v2.length().max(v1v3.length()).max(v2v3.length()) as f32;
    let billboard_scale = BILLBOARD_DEFAULT_SCALE * max_side_length * Vec3::ONE;

    commands.spawn((
        BillboardTextBundle {
            transform: Transform::from_translation(center).with_scale(billboard_scale),
            text: Text::from_sections([TextSection {
                value: format!("t{}", triangle_id.to_string()),
                style: TextStyle {
                    font_size: DEBUG_LABEL_FONT_SIZE,
                    color,
                    ..default()
                },
            }]),
            ..default()
        },
        TriangleDebugEntity,
    ));

    for (v, v_id, label) in vec![
        (v1, triangle.v1(), "v1"),
        (v2, triangle.v2(), "v2"),
        (v3, triangle.v3(), "v3"),
    ] {
        let v_display_pos = v + (center - v) * 0.15;
        let vertex_label = match vertex_label_mode {
            VertexLabelMode::LocalIndex => label.to_string(),
            VertexLabelMode::GlobalIndex => {
                if (v_id as usize) < vertices.len() {
                    v_id.to_string()
                } else {
                    format!("inf_{}", infinite_vertex_local_quad_index(v_id).to_string())
                }
            }
        };
        commands.spawn((
            BillboardTextBundle {
                transform: Transform::from_translation(v_display_pos).with_scale(billboard_scale),
                text: Text::from_sections([TextSection {
                    value: vertex_label,
                    style: TextStyle {
                        font_size: DEBUG_LABEL_FONT_SIZE,
                        color,
                        ..default()
                    },
                }]),
                ..default()
            },
            TriangleDebugEntity,
        ));
    }
}
