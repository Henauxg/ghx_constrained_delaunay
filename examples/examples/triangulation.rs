use bevy::{
    app::{App, Startup},
    ecs::system::Commands,
    math::Vec3,
    DefaultPlugins,
};
use examples::{
    extend_displayed_vertices_with_container_vertice, ExamplesPlugin, TriangleDebugPlugin,
    TrianglesDebugData,
};
use ghx_constrained_delaunay::triangulation::triangulation_from_3d_planar_vertices;

fn main() {
    App::new()
        .add_plugins((
            DefaultPlugins,
            ExamplesPlugin,
            TriangleDebugPlugin::default(),
        ))
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
