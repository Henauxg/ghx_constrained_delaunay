//! Runs the triangulations on random near-degenerate inputs and checks each output.
//!
//! Generators:
//! - `box`: CDT of a box whose border is the constrained loop, with vertices a few ULPs off its sides.
//! - `points`: Delaunay triangulation of clusters of points at very different scales, with points on the
//!   horizontal and vertical lines of other points.
//!
//! Usage:
//! - `cargo run --release -p fuzzer -- [runs] [first_seed]` runs all generators and prints each failing seed.
//! - `cargo run --release -p fuzzer -- --seed <generator> <seed>` prints the input of one seed as Rust code, to
//!   paste into a test.

use std::collections::HashMap;

use ghx_constrained_delaunay::{
    constrained_triangulation::ConstrainedTriangulationConfiguration,
    constrained_triangulation_from_2d_vertices,
    triangulation::TriangulationConfiguration,
    triangulation_from_2d_vertices,
    types::{Edge, Float, Vertex, VertexId},
};
use rand::{rngs::StdRng, Rng, SeedableRng};

struct Input {
    vertices: Vec<Vertex>,
    edges: Vec<Edge>,
    /// Area of the domain kept by the CDT
    expected_area: Float,
}

enum Outcome {
    Ok,
    Error(String),
    MissingEdge(Edge),
    InvertedTriangle([VertexId; 3]),
    FlatTriangle([VertexId; 3]),
    WrongArea {
        area: Float,
        expected: Float,
    },
    WrongTriangleCount {
        count: usize,
        expected: usize,
    },
    NotDelaunay {
        triangle: [VertexId; 3],
        vertex: VertexId,
    },
}

struct Generator {
    name: &'static str,
    generate: fn(u64) -> Input,
    check: fn(&Input) -> Outcome,
}

const GENERATORS: [Generator; 2] = [
    Generator {
        name: "box",
        generate: box_input,
        check: check_cdt,
    },
    Generator {
        name: "points",
        generate: points_input,
        check: check_delaunay,
    },
];

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.first().map(String::as_str) == Some("--seed") {
        let generator = args
            .get(1)
            .and_then(|name| GENERATORS.iter().find(|g| g.name == name))
            .expect("--seed <generator> <seed>");
        let seed = args
            .get(2)
            .and_then(|a| a.parse().ok())
            .expect("--seed <generator> <seed>");
        print_input(generator, seed);
        return;
    }
    let runs: u64 = args.first().map_or(10_000, |a| a.parse().expect("runs"));
    let first_seed: u64 = args.get(1).map_or(0, |a| a.parse().expect("first_seed"));

    for generator in &GENERATORS {
        let mut failures = 0;
        for seed in first_seed..first_seed + runs {
            if let Some(message) = describe((generator.check)(&(generator.generate)(seed))) {
                failures += 1;
                println!("{} seed {seed}: {message}", generator.name);
            }
        }
        println!("{}: {failures} failures out of {runs} runs", generator.name);
    }
}

/// Returns `None` if the outcome is a success
fn describe(outcome: Outcome) -> Option<String> {
    Some(match outcome {
        Outcome::Ok => return None,
        Outcome::Error(msg) => format!("error: {msg}"),
        Outcome::MissingEdge(edge) => format!("missing constrained edge {edge:?}"),
        Outcome::InvertedTriangle(t) => format!("inverted triangle {t:?}"),
        Outcome::FlatTriangle(t) => format!("flat triangle {t:?}"),
        Outcome::WrongArea { area, expected } => format!("area {area}, expected {expected}"),
        Outcome::WrongTriangleCount { count, expected } => {
            format!("{count} triangles, expected 2V - 2 - h = {expected}")
        }
        Outcome::NotDelaunay { triangle, vertex } => {
            format!("vertex {vertex} is inside the circumcircle of triangle {triangle:?}")
        }
    })
}

fn print_input(generator: &Generator, seed: u64) {
    let input = (generator.generate)(seed);
    let outcome = describe((generator.check)(&input)).unwrap_or_else(|| "ok".to_string());
    println!("// Fuzzer {} seed {seed}: {outcome}", generator.name);
    println!("let vertices = vec![");
    for v in &input.vertices {
        println!("    Vertex::new({:?}, {:?}),", v.x, v.y);
    }
    println!("];");
    if !input.edges.is_empty() {
        println!("let edges = vec![");
        for e in &input.edges {
            println!("    Edge::new({}, {}),", e.from, e.to);
        }
        println!("];");
    }
}

/// A box whose border is the constrained CW loop. Extra vertices on each side are moved
/// a few ULPs off the side, and some inner vertices sit very close to the right side.
fn box_input(seed: u64) -> Input {
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

/// Clusters of points at very different scales, some far from each other, plus points that share the x or y
/// of another point within a few ULPs: the infinite vertices are in the 4 axis directions.
fn points_input(seed: u64) -> Input {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut vertices = Vec::new();
    for _ in 0..rng.gen_range(1..5) {
        let scale = (10 as Float).powf(rng.gen_range(-4.0..2.0));
        let center = Vertex::new(rng.gen_range(-100.0..100.0), rng.gen_range(-100.0..100.0));
        for _ in 0..rng.gen_range(1..25) {
            let offset = Vertex::new(rng.gen_range(-1.0..1.0), rng.gen_range(-1.0..1.0));
            vertices.push(center + offset * scale);
        }
    }
    for _ in 0..rng.gen_range(0..8) {
        let base = vertices[rng.gen_range(0..vertices.len())];
        let other: Float = rng.gen_range(-100.0..100.0);
        let offset = rng.gen_range(-2..=2);
        vertices.push(if rng.gen_bool(0.5) {
            Vertex::new(shift_ulps(base.x, offset), other)
        } else {
            Vertex::new(other, shift_ulps(base.y, offset))
        });
    }
    Input {
        vertices,
        edges: Vec::new(),
        expected_area: 0.,
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

fn check_cdt(input: &Input) -> Outcome {
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

/// Checks the Delaunay triangulation in the normalized coordinates the library works with: no flat or inverted
/// triangle, a triangle count that matches the exact convex hull, and every shared edge locally Delaunay.
fn check_delaunay(input: &Input) -> Outcome {
    let triangulation = match triangulation_from_2d_vertices(
        &input.vertices,
        TriangulationConfiguration::default(),
    ) {
        Ok(t) => t,
        Err(e) => return Outcome::Error(e.msg),
    };
    let vertices = normalize(&input.vertices);
    if hull(&vertices, false).len() < 3 {
        // All the points are on one line
        return Outcome::Ok;
    }

    for t in &triangulation.triangles {
        let [a, b, c] = t.map(|v| coord(vertices[v as usize]));
        let orientation = robust::orient2d(a, b, c);
        if orientation == 0. {
            return Outcome::FlatTriangle(*t);
        }
        if orientation > 0. {
            return Outcome::InvertedTriangle(*t);
        }
    }

    let mut used: Vec<VertexId> = triangulation.triangles.iter().flatten().copied().collect();
    used.sort_unstable();
    used.dedup();
    let used_points: Vec<Vertex> = used.iter().map(|&v| vertices[v as usize]).collect();
    let expected = 2 * used.len() - 2 - hull(&used_points, true).len();
    if triangulation.triangles.len() != expected {
        return Outcome::WrongTriangleCount {
            count: triangulation.triangles.len(),
            expected,
        };
    }

    // For each edge, the first triangle that has it and the vertex of that triangle opposite to it
    let mut edges: HashMap<(VertexId, VertexId), ([VertexId; 3], VertexId)> = HashMap::new();
    for t in &triangulation.triangles {
        for k in 0..3 {
            let (a, b, opposite) = (t[k], t[(k + 1) % 3], t[(k + 2) % 3]);
            let key = (a.min(b), a.max(b));
            let Some((other, other_opposite)) = edges.insert(key, (*t, opposite)) else {
                continue;
            };
            if is_in_circumcircle(&vertices, other, opposite) {
                return Outcome::NotDelaunay {
                    triangle: other,
                    vertex: opposite,
                };
            }
            if is_in_circumcircle(&vertices, *t, other_opposite) {
                return Outcome::NotDelaunay {
                    triangle: *t,
                    vertex: other_opposite,
                };
            }
        }
    }
    Outcome::Ok
}

/// `triangle` is clockwise
fn is_in_circumcircle(vertices: &[Vertex], triangle: [VertexId; 3], vertex: VertexId) -> bool {
    let [a, b, c] = triangle.map(|v| coord(vertices[v as usize]));
    robust::incircle(c, b, a, coord(vertices[vertex as usize])) > 0.
}

/// Same operations as the library, so that the coordinates are identical
fn normalize(vertices: &[Vertex]) -> Vec<Vertex> {
    let (mut x_min, mut y_min, mut x_max, mut y_max) =
        (Float::MAX, Float::MAX, Float::MIN, Float::MIN);
    for v in vertices {
        x_min = x_min.min(v.x);
        x_max = x_max.max(v.x);
        y_min = y_min.min(v.y);
        y_max = y_max.max(v.y);
    }
    let scale = (x_max - x_min).max(y_max - y_min);
    vertices
        .iter()
        .map(|v| Vertex::new((v.x - x_min) / scale, (v.y - y_min) / scale))
        .collect()
}

/// Exact convex hull. With `keep_colinear`, also keeps the points on its edges.
fn hull(points: &[Vertex], keep_colinear: bool) -> Vec<Vertex> {
    let is_not_on_hull = |a: Vertex, b: Vertex, c: Vertex| {
        let orientation = robust::orient2d(coord(a), coord(b), coord(c));
        if keep_colinear {
            orientation < 0.
        } else {
            orientation <= 0.
        }
    };
    let mut sorted = points.to_vec();
    sorted.sort_by(|a, b| {
        a.x.partial_cmp(&b.x)
            .unwrap()
            .then(a.y.partial_cmp(&b.y).unwrap())
    });
    sorted.dedup();
    let half = |points: &mut dyn Iterator<Item = &Vertex>| {
        let mut chain: Vec<Vertex> = Vec::new();
        for &p in points {
            while chain.len() >= 2
                && is_not_on_hull(chain[chain.len() - 2], chain[chain.len() - 1], p)
            {
                chain.pop();
            }
            chain.push(p);
        }
        chain.pop();
        chain
    };
    let mut result = half(&mut sorted.iter());
    result.extend(half(&mut sorted.iter().rev()));
    result
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
