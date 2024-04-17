use std::f32::consts::PI;

use bevy::{
    app::{App, Plugin, Startup, Update},
    asset::Assets,
    core_pipeline::core_3d::Camera3dBundle,
    diagnostic::FrameTimeDiagnosticsPlugin,
    ecs::{
        component::Component,
        schedule::IntoSystemConfigs,
        system::{Commands, ResMut},
    },
    hierarchy::BuildChildren,
    input::{
        common_conditions::{input_just_pressed, input_toggle_active},
        keyboard::KeyCode,
    },
    math::{primitives::Cylinder, EulerRot, Quat, Vec3},
    pbr::{AmbientLight, DirectionalLight, DirectionalLightBundle, PbrBundle, StandardMaterial},
    render::{color::Color, mesh::Mesh},
    text::{BreakLineOn, Text, TextSection, TextStyle},
    transform::components::Transform,
    ui::{
        node_bundles::{NodeBundle, TextBundle},
        PositionType, Style, UiRect, Val,
    },
    utils::default,
};
use bevy_ghx_utils::{
    camera::{toggle_auto_orbit, update_pan_orbit_camera, PanOrbitCamera},
    systems::toggle_visibility,
};
use bevy_inspector_egui::quick::WorldInspectorPlugin;
use bevy_rapier3d::geometry::{ActiveCollisionTypes, Collider, Friction, Restitution};

use crate::{
    ball::BallPlugin,
    fps::{FpsDisplayPlugin, FpsRoot},
    DEFAULT_EXAMPLES_FONT_SIZE,
};

pub struct ExamplesPlugin;

impl Plugin for ExamplesPlugin {
    fn build(&self, app: &mut App) {
        app.add_plugins((
            FrameTimeDiagnosticsPlugin::default(),
            FpsDisplayPlugin,
            WorldInspectorPlugin::default().run_if(input_toggle_active(false, KeyCode::F3)),
            BallPlugin,
        ));
        app.add_systems(Startup, (setup_camera, setup_ui, setup_environment));
        app.add_systems(
            Update,
            (
                toggle_visibility::<ExamplesUiRoot>.run_if(input_just_pressed(KeyCode::F1)),
                toggle_visibility::<FpsRoot>.run_if(input_just_pressed(KeyCode::F2)),
                toggle_auto_orbit.run_if(input_just_pressed(KeyCode::F5)),
                update_pan_orbit_camera,
            ),
        );
    }
}

/// Marker to find the container entity so we can show/hide the UI node
#[derive(Component)]
pub struct ExamplesUiRoot;

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

pub fn setup_ui(mut commands: Commands) {
    let ui_root = commands
        .spawn((
            ExamplesUiRoot,
            NodeBundle {
                style: Style {
                    left: Val::Percent(1.),
                    height: Val::Vh(100.),
                    ..default()
                },
                ..default()
            },
            // Pickable::IGNORE, // TODO
        ))
        .id();
    let keybindings_text = "Toggles:\n\
        'F1' Show/hide UI\n\
        'F2' Show/hide fps\n\
        'F3' Show/hide inspector\n\
        'F5' Enable/disable camera rotation\n\
        "
    .to_string();

    let keybindings_ui_background = commands
        .spawn((
            // Pickable::IGNORE,
            NodeBundle {
                background_color: Color::BLACK.with_a(0.6).into(),
                style: Style {
                    position_type: PositionType::Absolute,
                    top: Val::Percent(1.),
                    padding: UiRect {
                        left: Val::Px(6.),
                        right: Val::Px(6.),
                        top: Val::Px(6.),
                        bottom: Val::Px(6.),
                    },
                    ..default()
                },
                ..default()
            },
        ))
        .id();
    let keybindings_ui = commands
        .spawn((
            // Pickable::IGNORE // TODO,
            TextBundle {
                style: Style {
                    position_type: PositionType::Relative,
                    ..default()
                },
                text: Text {
                    sections: vec![TextSection::new(
                        keybindings_text,
                        TextStyle {
                            font_size: DEFAULT_EXAMPLES_FONT_SIZE,
                            ..Default::default()
                        },
                    )],
                    linebreak_behavior: BreakLineOn::NoWrap,
                    ..default()
                },
                ..default()
            },
        ))
        .id();

    commands
        .entity(ui_root)
        .add_child(keybindings_ui_background);
    commands
        .entity(keybindings_ui_background)
        .add_child(keybindings_ui);
}

fn setup_environment(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    // Plane
    let radius = 2000.;
    let height = 200.;
    commands.spawn((
        PbrBundle {
            mesh: meshes.add(Cylinder::new(radius, height)),
            material: materials.add(Color::WHITE),
            transform: Transform::from_xyz(0.0, -height / 2., 0.0),
            ..default()
        },
        Collider::cylinder(height / 2., radius),
        (ActiveCollisionTypes::default()),
        Friction::coefficient(0.7),
        Restitution::coefficient(0.0),
    ));

    // Lights
    commands.insert_resource(AmbientLight {
        color: Color::WHITE,
        brightness: 50.,
    });
    commands.spawn(DirectionalLightBundle {
        transform: Transform::from_rotation(Quat::from_euler(
            EulerRot::ZYX,
            0.,
            -PI / 5.,
            -PI / 3.,
        )),
        directional_light: DirectionalLight {
            shadows_enabled: true,
            illuminance: 4000.,
            color: Color::WHITE,
            ..default()
        },
        ..default()
    });
    commands.spawn(DirectionalLightBundle {
        directional_light: DirectionalLight {
            shadows_enabled: false,
            illuminance: 1000.,
            color: Color::WHITE,
            ..default()
        },
        transform: Transform::from_rotation(Quat::from_euler(
            EulerRot::ZYX,
            0.,
            PI * 4. / 5.,
            -PI / 3.,
        )),
        ..default()
    });
}
