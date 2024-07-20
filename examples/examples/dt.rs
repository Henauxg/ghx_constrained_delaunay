use bevy::{
    app::{App, Startup},
    ecs::system::Commands,
    log::info,
    prelude::{EventWriter, IntoSystemConfigs, ResMut},
    DefaultPlugins,
};
use examples::{
    ExamplesPlugin, LabelMode, TriangleDebugCursorUpdate, TriangleDebugPlugin, TrianglesDebugData,
    TrianglesDebugViewConfig, TrianglesDrawMode, VertexLabelMode,
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
        .add_systems(Startup, (setup, display_last_snapshot).chain())
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
    info!("DT quality info: {:?}", q);

    commands.insert_resource(TrianglesDebugData::new(
        vertices
            .iter()
            .map(|v| Vector3::new(v.x, v.y, 0.))
            .collect(),
        triangulation.debug_context,
    ));
    commands.insert_resource(TrianglesDebugViewConfig::new(
        LabelMode::All,
        VertexLabelMode::LocalIndex,
        TrianglesDrawMode::AllAsGizmos,
        true,
    ));
}

pub fn display_last_snapshot(
    mut debug_data: ResMut<TrianglesDebugData>,
    mut debug_data_updates_events: EventWriter<TriangleDebugCursorUpdate>,
) {
    debug_data.set_cursor_to_last_snapshot();
    debug_data_updates_events.send(TriangleDebugCursorUpdate);
}
