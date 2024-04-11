use std::f32::consts::PI;

use bevy::{
    app::{App, Startup, Update},
    asset::Assets,
    core_pipeline::core_3d::Camera3dBundle,
    ecs::{
        schedule::IntoSystemConfigs,
        system::{Commands, ResMut},
    },
    input::{common_conditions::input_just_pressed, keyboard::KeyCode},
    math::{primitives::Plane3d, EulerRot, Quat, Vec3},
    pbr::{
        AmbientLight, CascadeShadowConfigBuilder, DirectionalLight, DirectionalLightBundle,
        PbrBundle, StandardMaterial,
    },
    render::{
        color::Color,
        mesh::{Mesh, Meshable},
    },
    transform::components::Transform,
    utils::default,
    DefaultPlugins,
};
use bevy_ghx_utils::camera::{pan_orbit_camera, toggle_auto_orbit, PanOrbitCamera};

fn main() {
    let mut app = App::new();

    app.insert_resource(AmbientLight {
        color: Color::WHITE,
        brightness: 2000.,
    })
    .add_plugins(DefaultPlugins);

    app.add_systems(Startup, (setup_camera, setup_sandbox));
    app.add_systems(
        Update,
        (
            toggle_auto_orbit.run_if(input_just_pressed(KeyCode::F5)),
            (pan_orbit_camera,),
        ),
    );

    app.run();
}

pub fn setup_camera(mut commands: Commands) {
    // Camera
    let camera_position = Vec3::new(0., 1.5, 2.5);
    let look_target = Vec3::ZERO;
    commands.spawn((
        Camera3dBundle {
            transform: Transform::from_translation(camera_position)
                .looking_at(look_target, Vec3::Y),
            ..default()
        },
        PanOrbitCamera {
            radius: (look_target - camera_position).length(),
            ..Default::default()
        },
    ));
}

pub fn setup_sandbox(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    // // Plane
    commands.spawn(PbrBundle {
        mesh: meshes.add(Plane3d::default().mesh().size(500000.0, 500000.0)),
        transform: Transform::from_xyz(0.0, -0.5, 0.0),
        material: materials.add(Color::rgb(0.3, 0.5, 0.3)),
        ..default()
    });

    // Light
    commands.spawn(DirectionalLightBundle {
        transform: Transform::from_rotation(Quat::from_euler(EulerRot::ZYX, 0.0, 1.0, -PI / 4.)),
        directional_light: DirectionalLight {
            shadows_enabled: true,
            ..default()
        },
        cascade_shadow_config: CascadeShadowConfigBuilder {
            first_cascade_far_bound: 200.0,
            maximum_distance: 400.0,
            ..default()
        }
        .into(),
        ..default()
    });
}
