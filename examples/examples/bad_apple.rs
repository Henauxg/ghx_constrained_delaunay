use std::{fs, time::Duration};

use bevy::{
    app::{App, Startup, Update},
    core_pipeline::{bloom::BloomSettings, tonemapping::Tonemapping},
    ecs::system::{Commands, Res},
    input::{common_conditions::input_just_pressed, ButtonInput},
    log::info,
    math::Vec3,
    prelude::{Camera3dBundle, EventWriter, IntoSystemConfigs, KeyCode, ResMut, Resource},
    render::camera::Camera,
    time::{Time, Timer, TimerMode},
    transform::components::Transform,
    utils::default,
    DefaultPlugins,
};
use examples::{
    camera::{update_pan_orbit_camera, PanOrbitCamera},
    update_triangles_debug_entities, LabelMode, TriangleDebugCursorUpdate, TriangleDebugPlugin,
    TrianglesDebugData, TrianglesDebugViewConfig, TrianglesDrawMode, VertexLabelMode,
};
use ghx_constrained_delaunay::{
    constrained_triangulation::ConstrainedTriangulationConfiguration,
    debug::{DebugConfiguration, Phase, PhaseRecord},
    types::{Edge, Float, Vector3, Vertex, VertexId},
    Triangulation,
};
use serde::Deserialize;

const FRAMES_FILE: &str = "./assets/bad_apple_frames.msgpack";
const DEFAULT_FRAME_DISPLAY_PERIOD_MICROS: u64 = 29850;
const DEFAULT_FRAME_SWITCH_MODE: FrameSwitchMode = FrameSwitchMode::Auto;

#[derive(Debug)]
pub enum FrameSwitchMode {
    Manual,
    Auto,
}

fn main() {
    App::new()
        .add_plugins((DefaultPlugins, TriangleDebugPlugin))
        .add_systems(Startup, (setup, setup_camera))
        .add_systems(
            Update,
            (
                update_pan_orbit_camera,
                switch_frame_mode.run_if(input_just_pressed(KeyCode::Space)),
                next_frame.before(update_triangles_debug_entities),
            ),
        )
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
            focus: Vec3::new(240.02797, -182.32127, 0.0),
            radius: 475.5253,
            upside_down: false,
            auto_orbit: false,
            auto_orbit_factor: 0.3,
            needs_transform_refresh: true,
        },
    ));
}

#[derive(Debug, Deserialize)]
pub struct Frame {
    pub vertices: Vec<(i32, i32)>,
    pub edges: Vec<(usize, usize)>,
}

pub struct TriangulatedFrame {
    pub triangulation: Triangulation,
    pub edges: Vec<Edge>,
    pub displayed_vertices: Vec<Vector3>,
}

#[derive(Resource)]
pub struct Frames {
    pub frame_index: usize,
    pub next_frame_timer: Timer,
    pub triangulated_frames: Vec<TriangulatedFrame>,
    pub switch_mode: FrameSwitchMode,
}

fn setup(mut commands: Commands) {
    info!("Deserializing frames from file...");
    let frames_bytes = fs::read(FRAMES_FILE).expect("Should have been able to load frames file");
    let frames: Vec<Frame> = rmp_serde::from_slice(&frames_bytes).unwrap();

    info!("Triangulating {} frames...", frames.len());

    let config = ConstrainedTriangulationConfiguration {
        debug_config: DebugConfiguration {
            phase_record: PhaseRecord::In(Phase::FilterTriangles),
            ..Default::default()
        },
        ..Default::default()
    };

    let mut triangulated_frames = Vec::new();
    for (_i, frame) in frames.iter().enumerate() {
        // if i != 897 {
        //     continue;
        // }
        // info!("Loading frame n°{}", _i);
        let edges = frame
            .edges
            .iter()
            .map(|e| Edge::new(e.0 as VertexId, e.1 as VertexId))
            .collect();
        let triangulation = ghx_constrained_delaunay::constrained_triangulation_from_2d_vertices(
            &frame
                .vertices
                .iter()
                .map(|v| Vertex::new(v.0 as Float, v.1 as Float))
                .collect(),
            &edges,
            config.clone(),
        );
        let displayed_vertices = frame
            .vertices
            .iter()
            .map(|v| Vector3::new(v.0 as Float, v.1 as Float, 0.))
            .collect();
        triangulated_frames.push(TriangulatedFrame {
            triangulation,
            edges,
            displayed_vertices,
        });
    }

    commands.insert_resource(TrianglesDebugData::new_with_constraints(
        triangulated_frames[0].displayed_vertices.clone(),
        &triangulated_frames[0].edges,
        triangulated_frames[0].triangulation.debug_context.clone(),
    ));

    commands.insert_resource(TrianglesDebugViewConfig::new(
        LabelMode::None,
        VertexLabelMode::GlobalIndex,
        TrianglesDrawMode::AllAsContourAndInteriorMeshes,
        false,
    ));

    commands.insert_resource(Frames {
        frame_index: 0,
        next_frame_timer: Timer::new(
            Duration::from_micros(DEFAULT_FRAME_DISPLAY_PERIOD_MICROS),
            TimerMode::Repeating,
        ),
        switch_mode: DEFAULT_FRAME_SWITCH_MODE,
        triangulated_frames,
    });

    info!("All {} frames are loaded and ready", frames.len());
}

fn next_frame(
    mut commands: Commands,
    time: Res<Time>,
    mut frames: ResMut<Frames>,
    keyboard_input: Res<ButtonInput<KeyCode>>,
    mut debug_data_updates_events: EventWriter<TriangleDebugCursorUpdate>,
) {
    let mut updated = false;
    match frames.switch_mode {
        FrameSwitchMode::Manual => {
            if keyboard_input.just_pressed(KeyCode::ArrowLeft)
                || keyboard_input.pressed(KeyCode::ArrowDown)
            {
                if frames.frame_index == 0 {
                    frames.frame_index = frames.triangulated_frames.len() - 1;
                } else {
                    frames.frame_index -= 1;
                }
                updated = true;
                info!("Frame n°{}", frames.frame_index);
            } else if keyboard_input.just_pressed(KeyCode::ArrowRight)
                || keyboard_input.pressed(KeyCode::ArrowUp)
            {
                frames.frame_index = (frames.frame_index + 1) % frames.triangulated_frames.len();
                updated = true;
                info!("Frame n°{}", frames.frame_index);
            }
        }
        FrameSwitchMode::Auto => {
            frames.next_frame_timer.tick(time.delta());
            if frames.next_frame_timer.finished() {
                frames.frame_index = (frames.frame_index + 1) % frames.triangulated_frames.len();
                updated = true;
            }
        }
    }

    if updated {
        let frame_data = &frames.triangulated_frames[frames.frame_index];
        commands.insert_resource(TrianglesDebugData::new_with_constraints(
            frame_data.displayed_vertices.clone(),
            &frame_data.edges,
            frame_data.triangulation.debug_context.clone(),
        ));
        debug_data_updates_events.send(TriangleDebugCursorUpdate);
    }
}

pub fn switch_frame_mode(mut frames: ResMut<Frames>) {
    frames.switch_mode = match frames.switch_mode {
        FrameSwitchMode::Manual => FrameSwitchMode::Auto,
        FrameSwitchMode::Auto => FrameSwitchMode::Manual,
    };
    info!("Frame switch mode set to {:?}", frames.switch_mode);
}
