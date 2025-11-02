use bevy::{
    ecs::system::Commands,
    log::error,
    math::Vec3,
    prelude::default,
    text::{Text, TextSection, TextStyle},
    transform::components::Transform,
};
use bevy_mod_billboard::BillboardTextBundle;
use ghx_constrained_delaunay::{infinite::collect_infinite_triangle_vertices};
use ghx_constrained_delaunay::{
    infinite::infinite_vertex_local_quad_index,
    types::{opposite_vertex_index_from_edge, TriangleData, TriangleId, VERT_1, VERT_2},
};

use crate::COLORS;

pub const DEBUG_LABEL_FONT_SIZE: f32 = 60.0;

pub const BILLBOARD_DEFAULT_SCALE: Vec3 = Vec3::splat(0.0005);

use super::{
    infinite::{
        get_2_extrapolated_vertices_from_1_finite_vertex,
        get_2_extrapolated_vertices_from_2_finite_vertices,
    },
    Basis, TriangleDebugEntity, VertexLabelMode,
};

pub fn spawn_label(
    commands: &mut Commands,
    vertex_label_mode: VertexLabelMode,
    vertices: &[Vec3],
    triangle_id: TriangleId,
    triangle: &TriangleData,
    basis: &Basis,
) {
    let color = COLORS[(triangle_id as usize) % COLORS.len()];

    let infinite_verts = collect_infinite_triangle_vertices(&triangle.verts);

    let (v1, v2, v3) = if infinite_verts.is_empty() {
        (
            vertices[triangle.v1() as usize],
            vertices[triangle.v2() as usize],
            vertices[triangle.v3() as usize],
        )
    } else if infinite_verts.len() == 1 {
        let (_a, finite_tip_a, _b, finite_tip_b) =
            get_2_extrapolated_vertices_from_2_finite_vertices(
                vertices,
                triangle,
                basis,
                infinite_verts[0],
            );

        if infinite_verts[0] == VERT_1 {
            (
                (finite_tip_a + finite_tip_b) / 2.,
                vertices[triangle.v2() as usize],
                vertices[triangle.v3() as usize],
            )
        } else if infinite_verts[0] == VERT_2 {
            (
                vertices[triangle.v1() as usize],
                (finite_tip_a + finite_tip_b) / 2.,
                vertices[triangle.v3() as usize],
            )
        } else {
            (
                vertices[triangle.v1() as usize],
                vertices[triangle.v2() as usize],
                (finite_tip_a + finite_tip_b) / 2.,
            )
        }
    } else if infinite_verts.len() == 2 {
        let finite_vertex_index =
            opposite_vertex_index_from_edge(infinite_verts[0], infinite_verts[1]);
        let (finite_tip_a, finite_vertex, finite_tip_b) =
            get_2_extrapolated_vertices_from_1_finite_vertex(
                vertices,
                triangle,
                basis,
                finite_vertex_index,
            );

        if finite_vertex_index == VERT_1 {
            (finite_vertex, finite_tip_a, finite_tip_b)
        } else if finite_vertex_index == VERT_2 {
            (finite_tip_b, finite_vertex, finite_tip_a)
        } else {
            (finite_tip_a, finite_tip_b, finite_vertex)
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
