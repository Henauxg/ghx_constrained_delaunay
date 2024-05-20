use bevy::prelude::*;

use ghx_constrained_delaunay::triangulation::triangulation_from_3d_planar_vertices;
use utils::{create_displayed_vertices, ExamplesPlugin, TriangleDebugPlugin, TrianglesDebugData};

mod utils;

fn main() {
    App::new()
        .add_plugins((DefaultPlugins, ExamplesPlugin, TriangleDebugPlugin))
        .add_systems(Startup, setup)
        .run();
}

fn setup(mut commands: Commands) {
    let vertices = vec![[0., 0., 0.], [0., 5., 0.], [5., 5., 0.], [5., 0., 0.]];
    let plane_normal = Vec3::Z;

    let (indices, debug_data) =
        triangulation_from_3d_planar_vertices(&vertices, plane_normal.into());

    let mut constrained_vertices = Vec::new();

    for vertex_id in indices {
        constrained_vertices.push(vertices[vertex_id]);
    }

    let displayed_vertices = create_displayed_vertices(vertices, plane_normal);
    commands.insert_resource(TrianglesDebugData::new(displayed_vertices, debug_data))
}
