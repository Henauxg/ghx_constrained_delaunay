use bevy::{prelude::*, utils::hashbrown::HashSet};

use ghx_constrained_delaunay::{
    constrained_triangulation::constrained_triangulation_from_3d_planar_vertices, types::Edge,
};

use utils::{create_displayed_vertices, ExamplesPlugin, TriangleDebugPlugin, TrianglesDebugData};

mod utils;

fn main() {
    App::new()
        .add_plugins((DefaultPlugins, ExamplesPlugin, TriangleDebugPlugin))
        .add_systems(Startup, setup)
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

    let constrained_edges: HashSet<Edge> = HashSet::from([
        Edge::new(0, 11),
        Edge::new(11, 10),
        Edge::new(10, 9),
        Edge::new(9, 8),
        Edge::new(8, 7),
        Edge::new(7, 6),
        Edge::new(6, 5),
        Edge::new(5, 4),
        Edge::new(4, 3),
        Edge::new(3, 2),
        Edge::new(2, 1),
        Edge::new(1, 0),
        Edge::new(12, 13),
        Edge::new(14, 15),
    ]);

    let plane_normal = Vec3::Z;
    let (_, debug_data) = constrained_triangulation_from_3d_planar_vertices(
        &vertices,
        plane_normal.into(),
        &constrained_edges,
    );

    let displayed_vertices = create_displayed_vertices(vertices, plane_normal);
    commands.insert_resource(TrianglesDebugData::new(displayed_vertices, debug_data))
}
