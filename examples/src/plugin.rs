use bevy::{
    app::{App, Plugin, Startup, Update},
    core_pipeline::core_3d::Camera3dBundle,
    ecs::system::Commands,
    input::ButtonInput,
    log::info,
    math::Vec3,
    prelude::{KeyCode, Query, Res},
    transform::components::Transform,
    utils::default,
};

use crate::camera::{update_pan_orbit_camera, PanOrbitCamera};

pub struct ExamplesPlugin;

impl Plugin for ExamplesPlugin {
    fn build(&self, app: &mut App) {
        app.add_systems(Startup, setup_camera);
        app.add_systems(Update, (update_pan_orbit_camera, camera_keyboard_control));
    }
}

pub fn setup_camera(mut commands: Commands) {
    let camera_position = Vec3::new(0., 0., 120.);
    let look_target = Vec3::ZERO;
    commands.spawn((
        Camera3dBundle {
            transform: Transform::from_translation(camera_position)
                .looking_at(look_target, Vec3::Y),
            ..default()
        },
        PanOrbitCamera {
            radius: (look_target - camera_position).length(),
            auto_orbit: false,
            ..Default::default()
        },
    ));
}

pub fn camera_keyboard_control(
    keyboard_input: Res<ButtonInput<KeyCode>>,
    cam_query: Query<&PanOrbitCamera>,
) {
    if keyboard_input.just_pressed(KeyCode::KeyC) {
        info!("Camera info: {:?}", cam_query.single());
    }
}
