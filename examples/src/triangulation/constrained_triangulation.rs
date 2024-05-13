use bevy::{prelude::*, utils::hashbrown::HashSet};

use bevy_ghx_destruction::triangulation::constrained_triangulation::constrained_triangulation_from_3d_planar_vertices;
use bevy_mod_billboard::plugin::BillboardPlugin;
use examples::{
    debug_utils::{
        create_displayed_vertices, draw_triangles_debug_data, update_triangles_debug_index,
        update_triangles_debugs_labels, TriangleDebugDataUpdated, TrianglesDebugData,
    },
    plugin::ExamplesPlugin,
};

fn main() {
    App::new()
        .add_plugins((DefaultPlugins, ExamplesPlugin, BillboardPlugin))
        .add_event::<TriangleDebugDataUpdated>()
        .add_systems(Startup, setup)
        .add_systems(
            Update,
            (
                update_triangles_debug_index,
                update_triangles_debugs_labels,
                draw_triangles_debug_data,
            )
                .chain(),
        )
        .run();
}

fn setup(mut commands: Commands) {
    let mut vertices = Vec::<[f32; 3]>::new();
    vertices.push([0., 0., 0.]);
    vertices.push([0., 6., 0.]);
    vertices.push([6., 6., 0.]);
    vertices.push([6., 0., 0.]);

    vertices.push([1., 1., 0.]);
    vertices.push([5., 1., 0.]);
    vertices.push([5., 5., 0.]);
    vertices.push([1., 5., 0.]);

    vertices.push([2., 2., 0.]);
    vertices.push([4., 2., 0.]);
    vertices.push([4., 4., 0.]);
    vertices.push([2., 4., 0.]);

    let mut constrained_edges = HashSet::new();
    constrained_edges.insert((0, 1));
    constrained_edges.insert((1, 2));
    constrained_edges.insert((2, 3));
    constrained_edges.insert((3, 0));

    constrained_edges.insert((4, 7));
    constrained_edges.insert((7, 6));
    constrained_edges.insert((6, 5));
    constrained_edges.insert((5, 4));

    constrained_edges.insert((8, 9));
    constrained_edges.insert((9, 10));
    constrained_edges.insert((10, 11));
    constrained_edges.insert((11, 8));

    let plane_normal = Vec3::Z;
    let (_, debug_data) = constrained_triangulation_from_3d_planar_vertices(
        &vertices,
        plane_normal.into(),
        &constrained_edges,
    );

    let displayed_vertices = create_displayed_vertices(vertices, plane_normal);

    commands.insert_resource(TrianglesDebugData {
        vertices: displayed_vertices,
        triangles_buffers: debug_data,
        current_buffer_index: 0,
    })
}
