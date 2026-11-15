//! Runs the CDT on random inputs where vertices sit a few ULPs off the border of a box,
//! and checks each output.
//!
//! Usage:
//! - `cargo run --release -p fuzzer -- [runs] [first_seed]` prints each failing seed.
//! - `cargo run --release -p fuzzer -- --seed <seed>` prints the input of one seed as Rust code, to paste into a test.

use ghx_constrained_delaunay::{
    constrained_triangulation::ConstrainedTriangulationConfiguration,
    constrained_triangulation_from_2d_vertices,
    types::{Edge, Float, Vertex, VertexId},
};
use rand::{rngs::StdRng, Rng, SeedableRng};

struct Input {
    vertices: Vec<Vertex>,
    edges: Vec<Edge>,
    expected_area: Float,
}

enum Outcome {
    Ok,
    Error(String),
    MissingEdge(Edge),
    InvertedTriangle([VertexId; 3]),
    WrongArea { area: Float, expected: Float },
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.first().map(String::as_str) == Some("--seed") {
        let seed = args
            .get(1)
            .and_then(|a| a.parse().ok())
            .expect("--seed <seed>");
        print_input(seed);
        return;
    }
    let runs: u64 = args.first().map_or(10_000, |a| a.parse().expect("runs"));
    let first_seed: u64 = args.get(1).map_or(0, |a| a.parse().expect("first_seed"));

    let mut failures = 0;
    for seed in first_seed..first_seed + runs {
        if let Some(message) = describe(check(&generate(seed))) {
            failures += 1;
            println!("seed {seed}: {message}");
        }
    }
    println!("{failures} failures out of {runs} runs");
}

/// Returns `None` if the outcome is a success
fn describe(outcome: Outcome) -> Option<String> {
    Some(match outcome {
        Outcome::Ok => return None,
        Outcome::Error(msg) => format!("error: {msg}"),
        Outcome::MissingEdge(edge) => format!("missing constrained edge {edge:?}"),
        Outcome::InvertedTriangle(t) => format!("inverted triangle {t:?}"),
        Outcome::WrongArea { area, expected } => format!("area {area}, expected {expected}"),
    })
}

fn print_input(seed: u64) {
    let input = generate(seed);
    let outcome = describe(check(&input)).unwrap_or_else(|| "ok".to_string());
    println!("// Fuzzer seed {seed}: {outcome}");
    println!("let vertices = vec![");
    for v in &input.vertices {
        println!("    Vertex::new({:?}, {:?}),", v.x, v.y);
    }
    println!("];");
    println!("let edges = vec![");
    for e in &input.edges {
        println!("    Edge::new({}, {}),", e.from, e.to);
    }
    println!("];");
}

/// A box whose border is the constrained CW loop. Extra vertices on each side are moved
/// a few ULPs off the side, and some inner vertices sit very close to the right side.
fn generate(seed: u64) -> Input {
    let mut rng = StdRng::seed_from_u64(seed);
    let w: Float = rng.gen_range(0.05..1.0);
    let h: Float = rng.gen_range(0.5..3.0);

    let corners = [(-w, 0.), (-w, h), (w, h), (w, 0.)];
    let mut vertices = Vec::new();
    for side in 0..4 {
        let (a, b) = (corners[side], corners[(side + 1) % 4]);
        vertices.push(Vertex::new(a.0, a.1));
        let mut ts: Vec<Float> = (0..rng.gen_range(0..4))
            .map(|_| rng.gen_range(0.05..0.95))
            .collect();
        ts.sort_by(|x, y| x.partial_cmp(y).unwrap());
        for t in ts {
            let mut x = a.0 + (b.0 - a.0) * t;
            let mut y = a.1 + (b.1 - a.1) * t;
            let offset = rng.gen_range(-3..=3);
            if a.0 == b.0 {
                x = shift_ulps(x, offset);
            } else {
                y = shift_ulps(y, offset);
            }
            vertices.push(Vertex::new(x, y));
        }
    }

    let loop_len = vertices.len();
    let edges = (0..loop_len)
        .map(|i| Edge::new(i as VertexId, ((i + 1) % loop_len) as VertexId))
        .collect();
    // The loop is CW, so its shoelace sum is negative.
    let expected_area = -0.5
        * (0..loop_len)
            .map(|i| vertices[i].perp_dot(vertices[(i + 1) % loop_len]))
            .sum::<Float>();

    for _ in 0..rng.gen_range(0..10) {
        let x = if rng.gen_bool(0.5) {
            shift_ulps(w * 0.999999, rng.gen_range(-3..3))
        } else {
            rng.gen_range(-w * 0.9..w * 0.9)
        };
        vertices.push(Vertex::new(x, rng.gen_range(0.05 * h..0.95 * h)));
    }

    Input {
        vertices,
        edges,
        expected_area,
    }
}

// `Float` is `f32` with the `f32` feature.
#[allow(clippy::unnecessary_cast)]
fn coord(v: Vertex) -> robust::Coord<f64> {
    robust::Coord {
        x: v.x as f64,
        y: v.y as f64,
    }
}

fn check(input: &Input) -> Outcome {
    let triangulation = match constrained_triangulation_from_2d_vertices(
        &input.vertices,
        &input.edges,
        ConstrainedTriangulationConfiguration::default(),
    ) {
        Ok(t) => t,
        Err(e) => return Outcome::Error(e.msg),
    };

    for edge in &input.edges {
        let found = triangulation
            .triangles
            .iter()
            .any(|t| t.contains(&edge.from) && t.contains(&edge.to));
        if !found {
            return Outcome::MissingEdge(*edge);
        }
    }

    let mut area = 0.;
    for t in &triangulation.triangles {
        let [a, b, c] = t.map(|v| input.vertices[v as usize]);
        let signed_area = 0.5 * (b - a).perp_dot(c - a);
        if robust::orient2d(coord(a), coord(b), coord(c)) > 0. {
            return Outcome::InvertedTriangle(*t);
        }
        area -= signed_area;
    }
    if (area - input.expected_area).abs() > 1e-6 * input.expected_area {
        return Outcome::WrongArea {
            area,
            expected: input.expected_area,
        };
    }
    Outcome::Ok
}

/// Moves `x` by `n` representable values. Zero moves by a tiny fixed step instead.
fn shift_ulps(x: Float, n: i32) -> Float {
    if x == 0. {
        return n as Float * Float::EPSILON * 1e-3;
    }
    let mut x = x;
    for _ in 0..n.unsigned_abs() {
        x = if (n > 0) == (x > 0.) {
            Float::from_bits(x.to_bits() + 1)
        } else {
            Float::from_bits(x.to_bits() - 1)
        };
    }
    x
}
