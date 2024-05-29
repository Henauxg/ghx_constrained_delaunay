use std::time::Instant;

use examples::extend_displayed_vertices_with_container_vertice;
use ghx_constrained_delaunay::{
    constrained_triangulation::ConstrainedTriangulationConfiguration,
    debug::{DebugConfiguration, PhaseRecord, TriangulationPhase},
    hashbrown::HashSet,
    triangulation::TriangulationConfiguration,
    types::{Edge, Float, Vector3, Vertex, VertexId},
    Triangulation,
};
use ordered_float::OrderedFloat;
use tracing_subscriber::{layer::SubscriberExt, Registry};
use tracing_tracy::TracyLayer;

fn main() {
    let subscriber = Registry::default().with(TracyLayer::new());
    tracing::subscriber::set_global_default(subscriber).expect("Failed to set subscriber");

    let shape_file_path = "../delaunay_compare/examples/Europe_coastline.shp";
    println!("Loading {} ...", shape_file_path);
    let mut reader = shapefile::Reader::from_path(shape_file_path).unwrap();

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
                        let x = p.x as Float;
                        let y = p.y as Float;
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
    println!(
        "loading time (ghx_cdt crate with constraints): {}ms",
        now.elapsed().as_millis()
    );

    triangulation
}
