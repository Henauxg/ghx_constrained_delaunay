use ghx_constrained_delaunay::triangulation_from_2d_vertices;
use ghx_constrained_delaunay::{triangulation::TriangulationConfiguration, types::Vertex};
use rand::rngs::StdRng;
use rand::{distributions::uniform::SampleUniform, Rng, SeedableRng};
use tracing_subscriber::{layer::SubscriberExt, Registry};
use tracing_tracy::TracyLayer;

pub const SEED: &[u8; 32] = b"\xfb\xdc\x4e\xa0\x30\xde\x82\xba\x69\x97\x3c\x52\x49\x4d\x00\xca
\x5c\x21\xa3\x8d\x5c\xf2\x34\x4e\x58\x7d\x80\x16\x66\x23\x30";

pub fn walk_f64() -> impl Iterator<Item = [f64; 2]> {
    random_walk_distribution(1.0, *SEED)
}

pub fn random_walk_distribution<S>(step_size: S, seed: [u8; 32]) -> impl Iterator<Item = [S; 2]>
where
    S: SampleUniform + num_traits::Float,
    S::Sampler: Copy,
{
    let range = rand::distributions::Uniform::new_inclusive(-step_size, step_size);
    let mut last_x = S::zero();
    let mut last_y = S::one();

    let mut rng = StdRng::from_seed(seed);
    let step_fn = move || {
        last_x = last_x + rng.sample(range);
        last_y = last_y + rng.sample(range);

        Some([last_x, last_y])
    };
    core::iter::from_fn(step_fn)
}

const VERTICES_COUNT: usize = 2000;

fn main() {
    let subscriber = Registry::default().with(TracyLayer::default());
    tracing::subscriber::set_global_default(subscriber).expect("Failed to set subscriber");

    let vertices: Vec<Vertex> = walk_f64()
        .take(VERTICES_COUNT)
        .map(|vertex| Vertex::new(vertex[0], vertex[1]))
        .collect();

    let _triangulation = triangulation_from_2d_vertices(
        &vertices,
        TriangulationConfiguration {
            ..Default::default()
        },
    );
}
