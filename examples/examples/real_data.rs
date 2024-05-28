use std::time::Instant;

use bevy::{
    app::{App, Startup, Update},
    ecs::system::Commands,
    gizmos::gizmos::Gizmos,
    math::{primitives::Direction3d, Vec3},
    render::color::Color,
    DefaultPlugins,
};
use examples::{
    extend_displayed_vertices_with_container_vertice, ExamplesPlugin, LabelMode,
    TriangleDebugPlugin, TrianglesDebugData, TrianglesDebugViewConfig, TrianglesDrawMode,
};
use ghx_constrained_delaunay::{
    constrained_triangulation::ConstrainedTriangulationConfiguration,
    debug::{DebugConfiguration, PhaseRecord, TriangulationPhase},
    hashbrown::HashSet,
    triangulation::TriangulationConfiguration,
    types::{Edge, Float, Vector3, VertexId, Vertex},
    Triangulation,
};
use ordered_float::OrderedFloat;

fn main() {
    App::new()
        .add_plugins((DefaultPlugins, ExamplesPlugin, TriangleDebugPlugin))
        .add_systems(Startup, setup)
        .add_systems(Update, draw_debug_circle)
        .run();
}

const RESIZE_FACTOR: Float = 1.0;

fn setup(mut commands: Commands) {
    let shape_file_path = "../delaunay_compare/examples/Europe_coastline.shp";
    println!("Loading {} ...", shape_file_path);
    let mut reader = shapefile::Reader::from_path(shape_file_path).unwrap();
    println!("Loaded!");

    println!("Extracting data...");
    let mut vertices = Vec::new();
    let mut edges = Vec::new();

    for shape_record in reader.iter_shapes_and_records() {
        let (shape, _) = shape_record.unwrap();

        let mut uniques: HashSet<[OrderedFloat<Float>; 2]> = HashSet::new();
        match shape {
            shapefile::Shape::Polyline(line) => {
                for part in line.parts() {
                    let first_vertex = vertices.len();
                    for p in part.iter() {
                        let x = RESIZE_FACTOR * p.x as Float;
                        let y = RESIZE_FACTOR * p.y as Float;
                        match uniques.insert([OrderedFloat(x), OrderedFloat(y)]) {
                            true => vertices.push(Vertex::new(x, y)),
                            false => (),
                        }
                    }
                    let last_vertex = vertices.len() - 1;
                    edges.extend((first_vertex..last_vertex).map(|i| [i, i + 1]));
                }
            }
            _ => unimplemented!(),
        }
    }

    println!("Done!");
    println!("{} vertices", vertices.len());
    println!("{} constraint edges", edges.len());

    vertices.shrink_to_fit();
    edges.shrink_to_fit();

    let triangulation = load_with_ghx_cdt_crate(&vertices, &edges);

    let plane_normal = Vector3::Z;
    let mut displayed_vertices = vertices
        .iter()
        .map(|v| Vector3::new(v.x, v.y, 0.))
        .collect();
    extend_displayed_vertices_with_container_vertice(
        &mut displayed_vertices,
        plane_normal,
        &triangulation.debug_context,
        false,
    );

    commands.insert_resource(TrianglesDebugData::new(
        displayed_vertices,
        triangulation.debug_context,
    ));
    commands.insert_resource(TrianglesDebugViewConfig::new(
        LabelMode::Changed,
        TrianglesDrawMode::AllAsMeshBatches { batch_size: 150 },
    ));
    // TODO Center camera on data
}

fn load_with_ghx_cdt_crate(vertices: &[Vertex], edges: &[[usize; 2]]) -> Triangulation {
    let vertices_clone = vertices.iter().map(|p| p.clone()).collect::<Vec<_>>();

    let config = ConstrainedTriangulationConfiguration {
        triangulation: TriangulationConfiguration {
            debug_config: DebugConfiguration {
                phase_record: PhaseRecord::InAny(HashSet::from([
                    TriangulationPhase::BeforeConstraints,
                    TriangulationPhase::AfterConstraints,
                    TriangulationPhase::RemoveWrapping,
                ])),
                ..Default::default()
            },
            ..Default::default()
        },
    };

    println!("Loading cdt (ghx_cdt crate)");
    let edges = edges
        .iter()
        // TODO into()
        .map(|[from, to]| Edge::new(*from as VertexId, *to as VertexId))
        .collect::<Vec<_>>();
    let now = Instant::now();
    let triangulation = ghx_constrained_delaunay::constrained_triangulation_from_2d_vertices(
        &vertices_clone,
        &edges,
        config,
    );
    println!("Done!");
    println!(
        "loading time (ghx_cdt crate with constraints): {}ms",
        now.elapsed().as_millis()
    );

    triangulation
}

fn draw_debug_circle(mut gizmos: Gizmos) {
    gizmos.circle(
        Vec3::new(0., 0., 0.),
        Direction3d::Z,
        100.,
        Color::ALICE_BLUE,
    );
}
