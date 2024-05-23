use bevy::{
    app::{App, Startup, Update},
    ecs::system::Commands,
    gizmos::gizmos::Gizmos,
    log::info,
    math::{primitives::Direction3d, Vec3},
    render::color::Color,
    DefaultPlugins,
};
use ghx_constrained_delaunay::triangulation::triangulation_from_3d_planar_vertices;
use glam::{Quat, Vec2};
use utils::{
    extend_displayed_vertices_with_container_vertice, ExamplesPlugin, TriangleDebugPlugin,
    TrianglesDebugData,
};

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

    let triangulation = triangulation_from_3d_planar_vertices(&vertices, plane_normal.into());

    let mut displayed_vertices = vertices.iter().map(|v| Vec3::from_slice(v)).collect();
    extend_displayed_vertices_with_container_vertice(
        &mut displayed_vertices,
        plane_normal,
        &triangulation.debug_context,
        true,
    );
    commands.insert_resource(TrianglesDebugData::new(
        displayed_vertices,
        triangulation.debug_context,
        true,
    ));
}
