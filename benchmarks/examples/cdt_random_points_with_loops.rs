//! CDT on random points with constrained loops, to vary how many constrained edges already exist after the
//! Delaunay step. In the coastline datasets of `cdt_benchmark`, almost all of them do.
//!
//! Run with `cargo run --release --example cdt_random_points_with_loops` in the `benchmarks` directory.

use std::{collections::HashSet, time::Instant};

use ghx_constrained_delaunay::{
    constrained_triangulation::ConstrainedTriangulationConfiguration,
    constrained_triangulation_from_2d_vertices,
    triangulation::TriangulationConfiguration,
    triangulation_from_2d_vertices,
    types::{Edge, Vertex, VertexId},
};
use rand::{rngs::StdRng, Rng, SeedableRng};

const RUNS: usize = 10;

struct Case {
    name: &'static str,
    random_points: usize,
    loops: usize,
    loop_vertices: usize,
}

const CASES: [Case; 2] = [
    // Few constrained edges, each much longer than the distance between points
    Case {
        name: "long constrained edges",
        random_points: 1_000_000,
        loops: 100,
        loop_vertices: 100,
    },
    // As many constrained edges as points, each 1 to 3 times the distance between points
    Case {
        name: "short constrained edges",
        random_points: 1_000_000,
        loops: 1000,
        loop_vertices: 1000,
    },
];

fn main() {
    for case in CASES {
        let (vertices, edges) = generate(&case);
        let already_present = already_present_ratio(&vertices, &edges);

        let mut times = Vec::with_capacity(RUNS);
        for _ in 0..RUNS {
            let now = Instant::now();
            constrained_triangulation_from_2d_vertices(
                &vertices,
                &edges,
                ConstrainedTriangulationConfiguration::default(),
            )
            .unwrap();
            times.push(now.elapsed().as_secs_f64() * 1000.);
        }
        times.sort_by(|a, b| a.partial_cmp(b).unwrap());

        println!(
            "{}: {} vertices, {} constrained edges, {:.1}% already in the Delaunay triangulation, CDT median {:.1} ms",
            case.name,
            vertices.len(),
            edges.len(),
            100. * already_present,
            times[RUNS / 2]
        );
    }
}

/// Concentric polygons as constrained loops, which never cross each other, plus random points in the unit square
fn generate(case: &Case) -> (Vec<Vertex>, Vec<Edge>) {
    let mut rng = StdRng::seed_from_u64(42);
    let mut vertices = Vec::new();
    let mut edges = Vec::new();
    for l in 0..case.loops {
        let radius = 0.5 * (l + 1) as f64 / case.loops as f64;
        let first = vertices.len();
        for i in 0..case.loop_vertices {
            // Clockwise, so that the loops delimit domains to keep
            let angle = -std::f64::consts::TAU * i as f64 / case.loop_vertices as f64;
            vertices.push(Vertex::new(
                0.5 + radius * angle.cos(),
                0.5 + radius * angle.sin(),
            ));
        }
        for i in 0..case.loop_vertices {
            edges.push(Edge::new(
                (first + i) as VertexId,
                (first + (i + 1) % case.loop_vertices) as VertexId,
            ));
        }
    }
    for _ in 0..case.random_points {
        vertices.push(Vertex::new(rng.gen(), rng.gen()));
    }
    (vertices, edges)
}

/// Share of the constrained edges that are edges of the Delaunay triangulation without constraints
fn already_present_ratio(vertices: &[Vertex], edges: &[Edge]) -> f64 {
    let triangulation =
        triangulation_from_2d_vertices(vertices, TriangulationConfiguration::default()).unwrap();
    let undirected = |a: VertexId, b: VertexId| (a.min(b), a.max(b));
    let delaunay_edges: HashSet<(VertexId, VertexId)> = triangulation
        .triangles
        .iter()
        .flat_map(|t| {
            [
                undirected(t[0], t[1]),
                undirected(t[1], t[2]),
                undirected(t[2], t[0]),
            ]
        })
        .collect();
    let present = edges
        .iter()
        .filter(|e| delaunay_edges.contains(&undirected(e.from, e.to)))
        .count();
    present as f64 / edges.len() as f64
}
