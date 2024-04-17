use bevy::{math::Vec3A, prelude::*};

use rand::Rng;
use utils::triangulation::triangulate;

pub mod utils;

fn main() {
    App::new()
        .add_plugins(DefaultPlugins)
        .init_gizmo_group::<MyRoundGizmos>()
        .add_systems(Startup, setup)
        .add_systems(Update, rotate_camera)
        .add_systems(Update, (draw_example_collection, update_config))
        .run();
}

// We can create our own gizmo config group!
#[derive(Default, Reflect, GizmoConfigGroup)]
struct MyRoundGizmos {}

fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    commands.spawn(Camera3dBundle {
        transform: Transform::from_xyz(0., 1.5, 6.).looking_at(Vec3::ZERO, Vec3::Y),
        ..default()
    });
    // plane
    commands.spawn(PbrBundle {
        mesh: meshes.add(Plane3d::default().mesh().size(5.0, 5.0)),
        material: materials.add(Color::rgb(0.3, 0.5, 0.3)),
        ..default()
    });

    // light
    commands.spawn(PointLightBundle {
        point_light: PointLight {
            shadows_enabled: true,
            ..default()
        },
        transform: Transform::from_xyz(4.0, 8.0, 4.0),
        ..default()
    });
}

fn rotate_camera(mut query: Query<&mut Transform, With<Camera>>, time: Res<Time>) {
    let mut transform = query.single_mut();

    transform.rotate_around(Vec3::ZERO, Quat::from_rotation_y(time.delta_seconds() / 2.));
}

fn draw_example_collection(mut gizmos: Gizmos) {
    let mut vertices = Vec::<[f32; 3]>::new();
    vertices.push([0., 3., 0.]);
    vertices.push([3., 3., 0.]);
    vertices.push([0., 0., 0.]);
    vertices.push([3., 0., 0.]);

    let triangulate_vertices = triangulate(&vertices, get_random_normalized_vec());

    for (triagle_id, triangle) in triangulate_vertices.chunks_exact(3).enumerate() {
        // let (v1, v2, v3) = (&triangle[0], &triangle[1], &triangle[2]);
        let color = COLORS[triagle_id % COLORS.len()];
        gizmos.linestrip_2d(
            vec![
                Vec2::from_array(triangle[0]),
                Vec2::from_array(triangle[1]),
                Vec2::from_array(triangle[2]),
                Vec2::from_array(triangle[0]),
            ],
            color,
        );

        // gizmos.line(
        //     triangulate_vertices[vertex_id],
        //     triangulate_vertices[vertex_id + 1],
        //     color,
        // )
    }
}

fn get_random_normalized_vec() -> Vec3A {
    let mut rng = rand::thread_rng();
    Vec3A::new(rng.gen::<f32>(), rng.gen::<f32>(), rng.gen::<f32>()).normalize()
}

fn update_config(
    mut config_store: ResMut<GizmoConfigStore>,
    keyboard: Res<ButtonInput<KeyCode>>,
    time: Res<Time>,
) {
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
