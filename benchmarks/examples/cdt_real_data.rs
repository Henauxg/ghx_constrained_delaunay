use std::time::Instant;

use anyhow::{Context, Ok};
use env_logger::Env;
use ghx_constrained_delaunay::{
    constrained_triangulation::ConstrainedTriangulationConfiguration,
    hashbrown::HashSet,
    triangulation::TriangulationConfiguration,
    types::{Edge, Float, Vertex, VertexId},
};
use ordered_float::OrderedFloat;
use spade::{ConstrainedDelaunayTriangulation, Point2, Triangulation};
use tiny_skia::{Paint, PathBuilder, Pixmap, Stroke, Transform};

const SHP_FILE_NAME: &str = "Europe_coastline";
// const SHP_FILE_NAME: &str = "ne_10m_coastline";
// const SHP_FILE_NAME: &str = "ne_50m_coastline";

const SHP_FILES_PATH: &str = "../assets";

fn main() -> anyhow::Result<()> {
    env_logger::Builder::from_env(Env::default().default_filter_or("debug")).init();

    let shape_file = format!("{SHP_FILES_PATH}/{SHP_FILE_NAME}.shp");
    println!("Loading {} ...", shape_file);
    let mut reader = shapefile::Reader::from_path(shape_file)?;

    let mut vertices = Vec::new();
    let mut edges = Vec::new();

    for shape_record in reader.iter_shapes_and_records() {
        let (shape, _) = shape_record?;

        let mut uniques: HashSet<[OrderedFloat<f64>; 2]> = HashSet::new();
        match shape {
            shapefile::Shape::Polyline(line) => {
                for part in line.parts() {
                    let first_vertex = vertices.len();
                    for p in part.iter() {
                        match uniques.insert([OrderedFloat(p.x), OrderedFloat(p.y)]) {
                            true => vertices.push(Point2::new(p.x, p.y)),
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
    println!();

    vertices.shrink_to_fit();
    edges.shrink_to_fit();

    load_with_spade(&vertices, &edges)?;
    println!();

    load_with_ghx_cdt_crate(&vertices, &edges)?;
    println!();

    load_with_cdt_crate(&vertices, &edges)?;

    Ok(())
}

fn load_with_spade(vertices: &Vec<Point2<f64>>, edges: &Vec<[usize; 2]>) -> anyhow::Result<()> {
    let vertices_clone = vertices.clone();
    let edges_clone = edges.clone();

    println!("Loading triangulation (spade)...");
    let now = Instant::now();
    let cdt =
        spade::ConstrainedDelaunayTriangulation::<_>::bulk_load_cdt(vertices_clone, edges_clone)?;
    let elapsed = now.elapsed().as_millis();
    println!("{} vertices (without duplicates)", cdt.num_vertices());
    println!("{} undirected edges", cdt.num_undirected_edges());
    println!("{} constraint edges", cdt.num_constraints());
    println!("{} triangles", cdt.num_inner_faces());
    println!("{} convex hull edges", cdt.convex_hull_size());
    println!();
    println!("loading time (spade with constraints): {}ms", elapsed);

    let vertices_clone = vertices.clone();
    let edges_clone = edges.clone();

    let now = Instant::now();
    spade::ConstrainedDelaunayTriangulation::<_>::bulk_load_cdt_stable(
        vertices_clone,
        edges_clone,
    )?;

    println!(
        "loading time (spade with constraints, stable): {}ms",
        now.elapsed().as_millis()
    );

    let vertices_clone = vertices.clone();

    let now = Instant::now();
    spade::ConstrainedDelaunayTriangulation::<_>::bulk_load(vertices_clone)?;
    println!(
        "loading time (spade without constraints): {}ms",
        now.elapsed().as_millis()
    );

    println!("Creating and saving output image...");
    _draw_to_pixmap(cdt)?.save_png(format!("{SHP_FILES_PATH}/{SHP_FILE_NAME}.png"))?;
    println!("Done!");

    Ok(())
}

fn load_with_cdt_crate(vertices: &[Point2<f64>], edges: &[[usize; 2]]) -> anyhow::Result<()> {
    let vertices_clone = vertices.iter().map(|p| (p.x, p.y)).collect::<Vec<_>>();

    println!("Loading cdt (cdt crate)");
    let edges = edges
        .iter()
        .map(|[from, to]| (*from, *to))
        .collect::<Vec<_>>();
    let now = Instant::now();
    cdt::triangulate_with_edges(&vertices_clone, &edges)?;
    println!(
        "loading time (cdt crate with constraints): {}ms",
        now.elapsed().as_millis()
    );

    let now = Instant::now();
    cdt::triangulate_points(&vertices_clone)?;
    println!(
        "loading time (cdt crate without constraints): {}ms",
        now.elapsed().as_millis()
    );

    Ok(())
}

fn load_with_ghx_cdt_crate(vertices: &[Point2<f64>], edges: &[[usize; 2]]) -> anyhow::Result<()> {
    let vertices_clone = vertices
        .iter()
        .map(|p| Vertex::new(p.x as Float, p.y as Float))
        .collect::<Vec<_>>();

    println!("Loading cdt (ghx_cdt crate)");
    let edges = edges
        .iter()
        // TODO into()
        .map(|[from, to]| Edge::new(*from as VertexId, *to as VertexId))
        .collect::<Vec<_>>();
    let now = Instant::now();
    ghx_constrained_delaunay::constrained_triangulation_from_2d_vertices(
        &vertices_clone,
        &edges,
        ConstrainedTriangulationConfiguration::default(),
    );
    println!(
        "loading time (ghx_cdt crate with constraints): {}ms",
        now.elapsed().as_millis()
    );

    let now = Instant::now();
    ghx_constrained_delaunay::triangulation_from_2d_vertices(
        &vertices_clone,
        TriangulationConfiguration::default(),
    );
    println!(
        "loading time (ghx_cdt crate without constraints): {}ms",
        now.elapsed().as_millis()
    );

    Ok(())
}

fn _draw_to_pixmap(cdt: ConstrainedDelaunayTriangulation<Point2<f64>>) -> anyhow::Result<Pixmap> {
    let mut min_x = f64::MAX;
    let mut min_y = f64::MAX;
    let mut max_x = f64::MIN;
    let mut max_y = f64::MIN;

    for position in cdt.convex_hull().map(|edge| edge.from().position()) {
        if position.x < min_x {
            min_x = position.x;
        }
        if position.y < min_y {
            min_y = position.y;
        }
        if position.x > max_x {
            max_x = position.x;
        }
        if position.y > max_y {
            max_y = position.y;
        }
    }

    let mut constraints = PathBuilder::new();
    let mut edges = PathBuilder::new();
    for ([from, to], is_constraint_edge) in cdt
        .undirected_edges()
        .map(|edge| (edge.positions(), edge.is_constraint_edge()))
    {
        if is_constraint_edge {
            constraints.move_to(from.x as f32, from.y as f32);
            constraints.line_to(to.x as f32, to.y as f32);
        } else {
            edges.move_to(from.x as f32, from.y as f32);
            edges.line_to(to.x as f32, to.y as f32);
        }
    }
    let constraints = constraints
        .finish()
        .context("Failed to finish constraint path")?;
    let edges = edges.finish().context("Failed to finish edge path")?;

    let res_x = 1024.0;
    let res_y = 1024.0;

    let mut pixmap = Pixmap::new(res_x as u32, res_y as u32).unwrap();

    let scale_x = res_x / (max_x - min_x);
    let scale_y = res_y / (max_y - min_y);

    let transform = Transform::from_translate(-min_x as f32, -min_y as f32)
        .post_scale(scale_x as f32, -scale_y as f32)
        .post_translate(0.0, res_y as f32);

    let mut stroke = Stroke::default();
    stroke.width = 0.1 / scale_x as f32;

    let mut paint = Paint::default();

    paint.set_color_rgba8(50, 127, 150, 255);

    pixmap.stroke_path(&edges, &paint, &stroke, transform, None);
    paint.set_color_rgba8(200, 126, 150, 255);
    stroke.width = 2.0 / scale_x as f32;
    pixmap.stroke_path(&constraints, &paint, &stroke, transform, None);

    Ok(pixmap)
}
