use bevy::{
    app::{App, Startup},
    ecs::system::Commands,
    log::info,
    DefaultPlugins,
};
use examples::{
    extend_displayed_vertices_with_container_vertice, ExamplesPlugin, LabelMode,
    TriangleDebugPlugin, TrianglesDebugData, TrianglesDebugViewConfig, TrianglesDrawMode,
    VertexLabelMode,
};
use ghx_constrained_delaunay::{
    triangulation::{triangulation_from_3d_planar_vertices, TriangulationConfiguration},
    types::{Vector3, Vertex},
    utils::check_delaunay_optimal,
};

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
    let plane_normal = Vector3::Z;

    let vertices = vec![
        Vector3::new(0., 0., 0.),
        Vector3::new(0., 5., 0.),
        Vector3::new(5., 5., 0.),
        Vector3::new(5., 0., 0.),
    ];
    let triangulation = triangulation_from_3d_planar_vertices(
        &vertices,
        plane_normal,
        TriangulationConfiguration::default(),
    );

    let q = check_delaunay_optimal(
        triangulation.triangles.iter().copied(),
        &vertices.iter().map(|v| Vertex::new(v.x, v.y)).collect(),
        false,
    );
    info!("DT quality: {:?}", q);

    let mut displayed_vertices = vertices
        .iter()
        .map(|v| Vector3::new(v.x, v.y, 0.))
        .collect();
    extend_displayed_vertices_with_container_vertice(
        &mut displayed_vertices,
        plane_normal,
        &triangulation.debug_context,
        true, // Change depending on planar transformation
    );
    commands.insert_resource(TrianglesDebugData::new(
        displayed_vertices,
        triangulation.debug_context,
    ));
    commands.insert_resource(TrianglesDebugViewConfig::new(
        LabelMode::All,
        VertexLabelMode::LocalIndex,
        TrianglesDrawMode::AllAsGizmos,
    ));
}
