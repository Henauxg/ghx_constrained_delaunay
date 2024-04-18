use bevy::prelude::*;

use examples::plugin::ExamplesPlugin;
use utils::triangulation::triangulate_3d_planar_vertices;

pub mod utils;

fn main() {
    App::new()
        .add_plugins((DefaultPlugins, ExamplesPlugin))
        .add_systems(Startup, setup)
        .add_systems(Update, draw_triangulation)
        .run();
}

fn setup(mut commands: Commands) {
    let mut vertices = Vec::<[f32; 3]>::new();
    vertices.push([0., 5., 0.]);
    vertices.push([5., 5., 0.]);
    vertices.push([5., 0., 0.]);
    vertices.push([0., 0., 0.]);

    let indices = triangulate_3d_planar_vertices(&vertices, Vec3::Z.into());

    commands.insert_resource(DebugVertData { vertices, indices })
}

#[derive(Resource)]
struct DebugVertData {
    vertices: Vec<[f32; 3]>,
    indices: Vec<usize>,
}

fn draw_triangulation(mut gizmos: Gizmos, debug_vert_data: Res<DebugVertData>) {
    for (triangle_id, triangle_indices) in debug_vert_data.indices.chunks_exact(3).enumerate() {
        let (v1, v2, v3) = (
            debug_vert_data.vertices[triangle_indices[0]],
            debug_vert_data.vertices[triangle_indices[1]],
            debug_vert_data.vertices[triangle_indices[2]],
        );
        let color = COLORS[triangle_id % COLORS.len()];
        gizmos.linestrip(
            vec![
                Vec3::from_array(v1),
                Vec3::from_array(v2),
                Vec3::from_array(v3),
                Vec3::from_array(v1),
            ],
            color,
        );
    }
}

const COLORS: &'static [Color] = &[
    Color::GREEN,
    Color::BLUE,
    Color::BLACK,
    Color::RED,
    Color::YELLOW,
    Color::MAROON,
    Color::PURPLE,
    Color::SALMON,
    Color::ORANGE,
    Color::CYAN,
    Color::NAVY,
    Color::OLIVE,
    Color::PINK,
    Color::ALICE_BLUE,
    Color::CRIMSON,
    Color::TURQUOISE,
    Color::YELLOW_GREEN,
    Color::TEAL,
];
