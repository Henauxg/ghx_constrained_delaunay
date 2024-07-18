use bevy::{
    app::{App, Startup},
    ecs::system::Commands,
    log::info,
    DefaultPlugins,
};
use examples::{
    ExamplesPlugin, LabelMode, TriangleDebugPlugin, TrianglesDebugData, TrianglesDebugViewConfig,
    TrianglesDrawMode, VertexLabelMode,
};
use ghx_constrained_delaunay::{
    triangulation::TriangulationConfiguration,
    triangulation_from_2d_vertices,
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
    let vertices = vec![
        Vertex::new(0., 0.),
        Vertex::new(0., 5.),
        Vertex::new(5., 5.),
        Vertex::new(5., 0.),
    ];
    let triangulation =
        triangulation_from_2d_vertices(&vertices, TriangulationConfiguration::default());

    let q = check_delaunay_optimal(
        triangulation.triangles.iter().copied(),
        &vertices.iter().map(|v| Vertex::new(v.x, v.y)).collect(),
        false,
    );
    info!("DT quality: {:?}", q);

    let displayed_vertices = vertices
        .iter()
        .map(|v| Vector3::new(v.x, v.y, 0.))
        .collect();
    let debug_data = TrianglesDebugData::new(displayed_vertices, triangulation.debug_context);
    commands.insert_resource(debug_data);
    commands.insert_resource(TrianglesDebugViewConfig::new(
        LabelMode::All,
        VertexLabelMode::LocalIndex,
        TrianglesDrawMode::AllAsGizmos,
        true,
    ));
}
