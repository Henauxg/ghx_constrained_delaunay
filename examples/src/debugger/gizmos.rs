use bevy::{
    color::{palettes::css::RED, Color},
    ecs::system::Res,
    gizmos::gizmos::Gizmos,
    math::{Dir3, Vec3},
};
use ghx_constrained_delaunay::{
    debug::EventInfo,
    types::{opposite_vertex_index_from_edge, TriangleData, TriangleId},
};
use ghx_constrained_delaunay::{glam::Quat, infinite::collect_infinite_triangle_vertices};

use crate::{
    infinite::{
        get_2_extrapolated_vertices_from_1_finite_vertex,
        get_2_extrapolated_vertices_from_2_finite_vertices,
    },
    COLORS,
};

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
                    triangle_id as TriangleId,
                    triangle,
                    &debug_data.basis,
                );
            }
        }
        TrianglesDrawMode::ChangedAsGizmos => {
            for triangle_id in snapshot.changed_ids.iter() {
                let triangle = &snapshot.triangles.get(*triangle_id);
                draw_triangle_gizmo(
                    &mut gizmos,
                    vertices,
                    *triangle_id,
                    triangle,
                    &debug_data.basis,
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
        gizmos.circle(vertices[vertex_id as usize], Dir3::Z, circle_radius, RED);
    }
}

pub fn draw_triangle_gizmo(
    gizmos: &mut Gizmos,
    vertices: &Vec<Vec3>,
    triangle_id: TriangleId,
    triangle: &TriangleData,
    basis: &Basis,
) {
    let color = COLORS[triangle_id as usize % COLORS.len()];

    let infinite_verts = collect_infinite_triangle_vertices(&triangle.verts);

    let linestrip = if infinite_verts.is_empty() {
        let (v1, v2, v3) = (
            vertices[triangle.v1() as usize],
            vertices[triangle.v2() as usize],
            vertices[triangle.v3() as usize],
        );
        vec![v1, v2, v3, v1]
    } else if infinite_verts.len() == 1 {
        let (finite_vertex_a, finite_tip_a, finite_vertex_b, finite_tip_b) =
            get_2_extrapolated_vertices_from_2_finite_vertices(
                vertices,
                triangle,
                basis,
                infinite_verts[0],
            );
        vec![finite_tip_a, finite_vertex_a, finite_vertex_b, finite_tip_b]
    } else if infinite_verts.len() == 2 {
        let (finite_tip_a, finite_vertex, finite_tip_b) =
            get_2_extrapolated_vertices_from_1_finite_vertex(
                vertices,
                triangle,
                basis,
                opposite_vertex_index_from_edge(infinite_verts[0], infinite_verts[1]),
            );
        vec![finite_tip_a, finite_vertex, finite_tip_b]
    } else {
        vec![]
        // TODO Could draw something for 3 infinite vertices
    };
    gizmos.linestrip(linestrip, color);
}
