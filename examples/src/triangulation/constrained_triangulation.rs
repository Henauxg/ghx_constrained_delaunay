use bevy::{prelude::*, utils::hashbrown::HashSet};

use bevy_ghx_destruction::triangulation::{
    constrained_triangulation::constrained_triangulation_from_3d_planar_vertices, Edge,
};
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
    let vertices = vec![
        [0., 0., 0.],
        [2., 0., 0.],
        [2., 2., 0.],
        [2., 4., 0.],
        [11., 4., 0.],
        [11., 1., 0.],
        [7., -1., 0.],
        [1., -2., 0.],
        [-5., 0., 0.],
        [-7., 4., 0.],
        [-4., 6., 0.],
        [0., 6., 0.],
        [6., 1., 0.],
        [8., 3., 0.],
        [-3., 1., 0.],
        [-4., 4., 0.],
    ];

    let mut constrained_edges: HashSet<Edge> = HashSet::new();

    constrained_edges.insert(Edge::new(0, 11));
    constrained_edges.insert(Edge::new(11, 10));
    constrained_edges.insert(Edge::new(10, 9));
    constrained_edges.insert(Edge::new(9, 8));
    constrained_edges.insert(Edge::new(8, 7));
    constrained_edges.insert(Edge::new(7, 6));
    constrained_edges.insert(Edge::new(6, 5));
    constrained_edges.insert(Edge::new(5, 4));
    constrained_edges.insert(Edge::new(4, 3));
    constrained_edges.insert(Edge::new(3, 2));
    constrained_edges.insert(Edge::new(2, 1));
    constrained_edges.insert(Edge::new(1, 0));

    constrained_edges.insert(Edge::new(12, 13));

    constrained_edges.insert(Edge::new(14, 15));

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
