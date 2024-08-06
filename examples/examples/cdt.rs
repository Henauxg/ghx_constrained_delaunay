use bevy::{
    app::{App, Startup},
    ecs::system::Commands,
    log::info,
    math::Vec3,
    prelude::{EventWriter, IntoSystemConfigs, ResMut},
    DefaultPlugins,
};

use examples::{
    ExamplesPlugin, LabelMode, TriangleDebugCursorUpdate, TriangleDebugPlugin, TrianglesDebugData,
    TrianglesDebugViewConfig, TrianglesDrawMode, VertexLabelMode,
};
use ghx_constrained_delaunay::{
    constrained_triangulation::ConstrainedTriangulationConfiguration,
    constrained_triangulation_from_2d_vertices,
    types::{Edge, Vertex},
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
        Vertex::new(2., 0.),
        Vertex::new(2., 2.),
        Vertex::new(2., 4.),
        Vertex::new(11., 4.),
        Vertex::new(11., 1.),
        Vertex::new(7., -1.),
        Vertex::new(1., -2.),
        Vertex::new(-5., 0.),
        Vertex::new(-7., 4.),
        Vertex::new(-4., 6.),
        Vertex::new(0., 6.),
        Vertex::new(6., 1.),
        Vertex::new(8., 3.),
        Vertex::new(-3., 1.),
        Vertex::new(-4., 4.),
    ];

    let constrained_edges = vec![
        Edge::new(11, 0),
        Edge::new(10, 11),
        Edge::new(9, 10),
        Edge::new(8, 9),
        Edge::new(7, 8),
        Edge::new(6, 7),
        Edge::new(5, 6),
        Edge::new(4, 5),
        Edge::new(3, 4),
        Edge::new(2, 3),
        Edge::new(1, 2),
        Edge::new(0, 1),
        Edge::new(13, 12),
        Edge::new(15, 14),
    ];

    let triangulation = constrained_triangulation_from_2d_vertices(
        &vertices,
        &constrained_edges,
        ConstrainedTriangulationConfiguration::default(),
    )
    .unwrap();

    let delaunay_quality = check_delaunay_optimal(
        triangulation.triangles,
        &vertices
            .iter()
            .map(|v| Vertex::new(v[0], v[1]))
            .collect::<Vec<Vertex>>(),
        false,
    );
    info!("CDT quality info: {:?}", delaunay_quality);

    commands.insert_resource(TrianglesDebugData::new_with_constrained_edges(
        vertices
            .iter()
            .map(|v| Vec3::new(v.x as f32, v.y as f32, 0.))
            .collect(),
        &constrained_edges,
        triangulation.debug_context,
    ));
    commands.insert_resource(TrianglesDebugViewConfig::new(
        LabelMode::All,
        VertexLabelMode::LocalIndex,
        TrianglesDrawMode::ContoursAndInteriorAsMeshes,
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
