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
    extend_displayed_vertices_with_container_vertice, ExamplesPlugin, TriangleDebugPlugin,
    TrianglesDebugData,
};
use ghx_constrained_delaunay::Triangulation;
use glam::Vec2;

fn main() {
    App::new()
        .add_plugins((
            DefaultPlugins,
            ExamplesPlugin,
            TriangleDebugPlugin {
                draw_mode: examples::DrawMode::MeshBatches { batch_size: 15 },
            },
        ))
        .add_systems(Startup, setup)
        .add_systems(Update, draw_debug_circle)
        .run();
}

const RESIZE_FACTOR: f32 = 0.001;

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

        match shape {
            shapefile::Shape::Polyline(line) => {
                for part in line.parts() {
                    let first_vertex = vertices.len();
                    vertices.extend(part.iter().map(|p| {
                        Vec2::new(RESIZE_FACTOR * p.x as f32, RESIZE_FACTOR * p.y as f32)
                    }));
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

    // let plane_normal = Vec3::Z;
    // let displayed_vertices = custom_create_displayed_vertices(
    //     vertices.iter().map(|v| Vec3::new(v.x, v.y, 0.)).collect(),
    //     plane_normal,
    //     &triangulation.debug_context,
    // );

    let plane_normal = Vec3::Z;
    let mut displayed_vertices = vertices.iter().map(|v| Vec3::new(v.x, v.y, 0.)).collect();
    extend_displayed_vertices_with_container_vertice(
        &mut displayed_vertices,
        plane_normal,
        &triangulation.debug_context,
        false,
    );

    commands.insert_resource(TrianglesDebugData::new(
        displayed_vertices,
        triangulation.debug_context,
        false,
    ));
}

fn load_with_ghx_cdt_crate(vertices: &[Vec2], _edges: &[[usize; 2]]) -> Triangulation {
    let vertices_clone = vertices
        .iter()
        .map(|p| Vec2::new(p.x as f32, p.y as f32))
        .collect::<Vec<_>>();

    // println!("Loading cdt (ghx_cdt crate)");
    // let edges = edges
    //     .iter()
    //     .map(|[from, to]| Edge::new(*from, *to))
    //     .collect::<HashSet<_>>();
    // let now = Instant::now();
    // ghx_constrained_delaunay::constrained_triangulation_from_2d_vertices(&vertices_clone, &edges);
    // println!("Done!");
    // println!(
    //     "loading time (ghx_cdt crate with constraints): {}ms",
    //     now.elapsed().as_millis()
    // );

    let now = Instant::now();
    let triangulation = ghx_constrained_delaunay::triangulation_from_2d_vertices(&vertices_clone);
    println!("Done!");
    println!(
        "loading time (ghx_cdt crate without constraints): {}ms",
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
