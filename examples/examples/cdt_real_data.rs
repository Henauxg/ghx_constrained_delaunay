use std::time::Instant;

use bevy::{
    app::{App, Startup, Update},
    color::palettes::css::ALICE_BLUE,
    ecs::system::{Commands, Res},
    gizmos::gizmos::Gizmos,
    log::{error, info},
    math::{Dir3, Vec3},
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
    debug::{DebugConfiguration, Phase, PhaseRecord},
    hashbrown::HashSet,
    types::{Edge, Vertex, VertexId},
    utils::check_degenerate_triangles,
};

const SHP_FILE_NAME: &str = "light/ne_50m_coastline";
// const SHP_FILE_NAME: &str = "heavy/ne_10m_coastline";
// const SHP_FILE_NAME: &str = "heavy/Europe_coastline";

// Modify to change the display scale
const DISPLAY_SCALE: f32 = 1.;

const SHP_FILES_PATH: &str = "../assets";

fn main() {
    App::new()
        .add_plugins((DefaultPlugins, ExamplesPlugin, TriangleDebugPlugin))
        .add_systems(Startup, (setup, display_last_snapshot).chain())
        .add_systems(Update, draw_origin_circle)
        .run();
}

fn setup(mut commands: Commands) {
    let shape_file_path = format!("{SHP_FILES_PATH}/{SHP_FILE_NAME}.shp");
    info!("Loading {} ...", shape_file_path);
    let mut reader = shapefile::Reader::from_path(shape_file_path).unwrap();

    let mut vertices = Vec::new();
    let mut constrained_edges = Vec::new();

    for shape_record in reader.iter_shapes_and_records() {
        let (shape, _) = shape_record.unwrap();

        match shape {
            shapefile::Shape::Polyline(line) => {
                for part in line.parts() {
                    let first_vertex = vertices.len();
                    vertices.extend(part.iter().map(|p| Vertex::new(p.x, p.y)));
                    let last_vertex = vertices.len() - 1;
                    constrained_edges.extend(
                        (first_vertex..last_vertex)
                            .map(|i| Edge::new(i as VertexId, (i + 1) as VertexId)),
                    );
                }
            }
            _ => unimplemented!(),
        }
    }

    info!("{} vertices", vertices.len());
    info!("{} constraint edges", constrained_edges.len());

    vertices.shrink_to_fit();
    constrained_edges.shrink_to_fit();
    let config = ConstrainedTriangulationConfiguration {
        debug_config: DebugConfiguration {
            phase_record: PhaseRecord::InAny(HashSet::from([
                Phase::BeforeConstraints,
                Phase::AfterConstraints,
                Phase::FilterTriangles,
            ])),
            ..Default::default()
        },
        ..Default::default()
    };

    info!("Loading cdt (ghx_cdt crate)");
    let now = Instant::now();
    let triangulation_result =
        constrained_triangulation_from_2d_vertices(&vertices, &constrained_edges, config);
    let debug_context = match triangulation_result {
        Ok(triangulation) => {
            let delaunay_quality =
                check_degenerate_triangles(triangulation.triangles.iter().copied(), &vertices);
            info!("Delaunay quality info:: {:?}", delaunay_quality);
            info!(
                "loading time (ghx_cdt crate with constraints): {}ms",
                now.elapsed().as_millis()
            );
            triangulation.debug_context
        }
        Err(err) => {
            error!("Failed triangulation: {:?}", err.msg);
            err.debug_context
        }
    };

    let displayed_vertices = vertices
        .iter()
        .map(|v| DISPLAY_SCALE * Vec3::new(v.x as f32, v.y as f32, 0.))
        .collect();
    commands.insert_resource(TrianglesDebugData::new_with_constrained_edges(
        displayed_vertices,
        &constrained_edges,
        debug_context,
    ));
    commands.insert_resource(TrianglesDebugViewConfig::new(
        LabelMode::Changed,
        VertexLabelMode::GlobalIndex,
        TrianglesDrawMode::ContoursAndInteriorAsMeshes,
        true,
    ));
    // TODO Center camera on data
}

fn draw_origin_circle(mut gizmos: Gizmos, triangle_debug_data: Res<TrianglesDebugData>) {
    gizmos.circle(
        Vec3::new(0., 0., 0.),
        Dir3::Z,
        0.01 * triangle_debug_data.context.scale_factor as f32,
        ALICE_BLUE,
    );
}

pub fn display_last_snapshot(
    mut debug_data: ResMut<TrianglesDebugData>,
    mut debug_data_updates_events: EventWriter<TriangleDebugCursorUpdate>,
) {
    debug_data.set_cursor_to_last_snapshot();
    debug_data_updates_events.send(TriangleDebugCursorUpdate);
}
