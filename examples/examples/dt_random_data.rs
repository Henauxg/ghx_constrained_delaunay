use bevy::app::{App, Startup};
use bevy::prelude::{Commands, IntoSystemConfigs};
use bevy::DefaultPlugins;
use examples::{ExamplesPlugin, TriangleDebugPlugin};
use ghx_constrained_delaunay::debug::{DebugConfiguration, PhaseRecord, StepsRecord};
use ghx_constrained_delaunay::triangulation_from_2d_vertices;
use ghx_constrained_delaunay::{triangulation::TriangulationConfiguration, types::Vertex};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
// use tracing_subscriber::{layer::SubscriberExt, Registry};
// use tracing_tracy::TracyLayer;

pub const SEED: &[u8; 32] = b"\xfb\xdc\x4e\xa0\x30\xde\x82\xba\x69\x97\x3c\x52\x49\x4d\x00\xca
\x5c\x21\xa3\x8d\x5c\xf2\x34\x4e\x58\x7d\x80\x16\x66\x23\x30";

pub fn walk_f64() -> impl Iterator<Item = [f64; 2]> {
    random_walk_distribution(1.0, *SEED)
}

fn main() {
    App::new()
        .add_plugins((
            DefaultPlugins,
            ExamplesPlugin,
            TriangleDebugPlugin::default(),
        ))
        .add_systems(Startup, (setup).chain())
        .run();
}

pub fn random_walk_distribution(step_size: f64, seed: [u8; 32]) -> impl Iterator<Item = [f64; 2]> {
    let range = rand::distributions::Uniform::new_inclusive(-step_size, step_size);
    let mut last_x = 0.;
    let mut last_y = 1.;

    let mut rng = StdRng::from_seed(seed);
    let step_fn = move || {
        last_x = last_x + rng.sample(range);
        last_y = last_y + rng.sample(range);

        Some([last_x, last_y])
    };
    core::iter::from_fn(step_fn)
}

const VERTICES_COUNT: usize = 550_000;

fn setup(mut commands: Commands) {
    let vertices: Vec<Vertex> = walk_f64()
        .take(VERTICES_COUNT)
        .map(|vertex| Vertex::new(vertex[0], vertex[1]))
        .collect();

    let _triangulation = triangulation_from_2d_vertices(
        &vertices,
        TriangulationConfiguration {
            debug_config: DebugConfiguration {
                phase_record: PhaseRecord::All,
                steps_record: StepsRecord::None,
                force_end_at_step: None,
            },
            ..Default::default()
        },
    )
    .unwrap();
}
