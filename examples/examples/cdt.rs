use bevy::{
    app::{App, Startup},
    ecs::system::Commands,
    log::info,
    DefaultPlugins,
};

use examples::{
    extend_displayed_vertices_with_container_vertice, ExamplesPlugin, LabelMode,
    TriangleDebugPlugin, TrianglesDebugData, TrianglesDebugViewConfig, TrianglesDrawMode,
};
use ghx_constrained_delaunay::{
    constrained_triangulation::{
        constrained_triangulation_from_3d_planar_vertices, ConstrainedTriangulationConfiguration,
    },
    types::{Edge, Vector3, Vertex},
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
        Vector3::new(0., 0., 0.),
        Vector3::new(2., 0., 0.),
        Vector3::new(2., 2., 0.),
        Vector3::new(2., 4., 0.),
        Vector3::new(11., 4., 0.),
        Vector3::new(11., 1., 0.),
        Vector3::new(7., -1., 0.),
        Vector3::new(1., -2., 0.),
        Vector3::new(-5., 0., 0.),
        Vector3::new(-7., 4., 0.),
        Vector3::new(-4., 6., 0.),
        Vector3::new(0., 6., 0.),
        Vector3::new(6., 1., 0.),
        Vector3::new(8., 3., 0.),
        Vector3::new(-3., 1., 0.),
        Vector3::new(-4., 4., 0.),
    ];

    let constrained_edges = vec![
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
    ];

    let plane_normal = Vector3::Z;
    let triangulation = constrained_triangulation_from_3d_planar_vertices(
        &vertices,
        plane_normal.into(),
        &constrained_edges,
        ConstrainedTriangulationConfiguration::default(),
    );

    let delaunay_quality = check_delaunay_optimal(
        triangulation.triangles,
        &vertices.iter().map(|v| Vertex::new(v[0], v[1])).collect(),
        false,
    );
    info!("Delaunay quality info:: {:?}", delaunay_quality);

    let mut displayed_vertices = vertices.clone();
    extend_displayed_vertices_with_container_vertice(
        &mut displayed_vertices,
        plane_normal,
        &triangulation.debug_context,
        true,
    );
    commands.insert_resource(TrianglesDebugData::new(
        displayed_vertices,
        triangulation.debug_context,
    ));
    commands.insert_resource(TrianglesDebugViewConfig::new(
        LabelMode::All,
        TrianglesDrawMode::AllAsGizmos,
    ));
}
