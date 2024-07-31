use std::fs;

use bevy::{
    app::{App, Startup, Update},
    core_pipeline::{bloom::BloomSettings, tonemapping::Tonemapping},
    ecs::system::Commands,
    input::ButtonInput,
    log::{error, info},
    math::Vec3,
    prelude::{
        Camera, Camera3dBundle, EventWriter, IntoSystemConfigs, KeyCode, Query, Res, ResMut,
        Transform,
    },
    utils::default,
    DefaultPlugins,
};

use examples::{
    camera::{update_pan_orbit_camera, PanOrbitCamera},
    LabelMode, TriangleDebugCursorUpdate, TriangleDebugPlugin, TrianglesDebugData,
    TrianglesDebugViewConfig, TrianglesDrawMode, VertexLabelMode,
};
use ghx_constrained_delaunay::{
    constrained_triangulation::ConstrainedTriangulationConfiguration,
    constrained_triangulation_from_2d_vertices,
    debug::{DebugConfiguration, PhaseRecord},
    types::{Edge, Float, Vertex, VertexId},
};
use serde::Deserialize;

pub const FRAME_INDEX: usize = 497;

// Low quality settings
// First 500 frames, as low quality
const FRAMES_FILE: &str = "../assets/light/bad_apple_frames_part.msgpack";
// Whole video, around 6500 frames, low quality
// const FRAMES_FILE: &str = "../assets/heavy/bad_apple_frames_lq.msgpack";
const CAMERA_FOCUS: Vec3 = Vec3::new(240.02797, -182.32127, 0.0);
const CAMERA_RADIUS: f32 = 475.5253;

// High quality settings
// Whole video, around 13000 frames, high quality
// const FRAMES_FILE: &str = "../assets/heavy/bad_apple_frames_hq.msgpack";
// const CAMERA_FOCUS: Vec3 = Vec3::new(1284.7368, -960.08685, 0.0);
// const CAMERA_RADIUS: f32 = 2327.6062;
// const DEFAULT_FRAME_DISPLAY_PERIOD_MICROS: u64 = 16500;

#[derive(Debug, Deserialize)]
pub struct Frame {
    pub vertices: Vec<(i32, i32)>,
    pub edges: Vec<(usize, usize)>,
}

fn main() {
    App::new()
        .add_plugins((DefaultPlugins, TriangleDebugPlugin::default()))
        .add_systems(
            Startup,
            (setup, setup_camera, display_first_snapshot).chain(),
        )
        .add_systems(Update, (update_pan_orbit_camera, camera_keyboard_control))
        .run();
}

pub fn setup_camera(mut commands: Commands) {
    let camera_position = Vec3::new(0., 0., 100.);
    let look_target = Vec3::ZERO;
    commands.spawn((
        Camera3dBundle {
            transform: Transform::from_translation(camera_position)
                .looking_at(look_target, Vec3::Y),
            camera: Camera {
                hdr: true,
                ..default()
            },
            tonemapping: Tonemapping::TonyMcMapface,
            ..default()
        },
        BloomSettings::NATURAL,
        PanOrbitCamera {
            focus: CAMERA_FOCUS,
            radius: CAMERA_RADIUS,
            upside_down: false,
            auto_orbit: false,
            auto_orbit_factor: 0.3,
            needs_transform_refresh: true,
        },
    ));
}

fn setup(mut commands: Commands) {
    info!("Deserializing frames from file...");
    let frames_bytes = fs::read(FRAMES_FILE).expect("Should have been able to load frames file");
    let frames: Vec<Frame> = rmp_serde::from_slice(&frames_bytes).unwrap();

    let frame = &frames[FRAME_INDEX];

    let config = ConstrainedTriangulationConfiguration {
        debug_config: DebugConfiguration {
            phase_record: PhaseRecord::All,
            ..Default::default()
        },
        ..Default::default()
    };

    let edges = frame
        .edges
        .iter()
        .map(|e| Edge::new(e.0 as VertexId, e.1 as VertexId))
        .collect();
    let triangulation_result = constrained_triangulation_from_2d_vertices(
        &frame
            .vertices
            .iter()
            .map(|v| Vertex::new(v.0 as Float, v.1 as Float))
            .collect(),
        &edges,
        config,
    );
    let debug_context = match triangulation_result {
        Ok(triangulation) => triangulation.debug_context,
        Err(err) => {
            error!("Failed triangulation: {:?}", err.msg);
            err.debug_context
        }
    };

    let displayed_vertices = frame
        .vertices
        .iter()
        .map(|v| Vec3::new(v.0 as f32, v.1 as f32, 0.))
        .collect();
    commands.insert_resource(TrianglesDebugData::new_with_constrained_edges(
        displayed_vertices,
        &edges,
        debug_context,
    ));
    commands.insert_resource(TrianglesDebugViewConfig::new(
        LabelMode::All,
        VertexLabelMode::GlobalIndex,
        TrianglesDrawMode::AllAsGizmos,
        true,
    ));
}

pub fn display_first_snapshot(
    mut debug_data: ResMut<TrianglesDebugData>,
    mut debug_data_updates_events: EventWriter<TriangleDebugCursorUpdate>,
) {
    debug_data.set_cursor_to_first_snapshot();
    debug_data_updates_events.send(TriangleDebugCursorUpdate);
}

pub fn camera_keyboard_control(
    keyboard_input: Res<ButtonInput<KeyCode>>,
    mut cam_query: Query<&mut PanOrbitCamera>,
) {
    if keyboard_input.just_pressed(KeyCode::KeyC) {
        info!("Camera info: {:?}", cam_query.single());
    }
    if keyboard_input.pressed(KeyCode::NumpadAdd) {
        let mut cam = cam_query.single_mut();
        cam.radius -= 1.;
        cam.needs_transform_refresh = true;
    }
    if keyboard_input.pressed(KeyCode::NumpadSubtract) {
        let mut cam = cam_query.single_mut();
        cam.radius += 1.;
        cam.needs_transform_refresh = true;
    }
}
