use std::{fs, time::Duration};

use bevy::{
    app::{App, Startup, Update},
    core_pipeline::{bloom::BloomSettings, tonemapping::Tonemapping},
    ecs::system::{Commands, Res},
    gizmos::gizmos::Gizmos,
    log::info,
    math::{primitives::Direction3d, Vec3},
    prelude::{Camera3dBundle, EventWriter, ResMut, Resource},
    render::{camera::Camera, color::Color},
    time::{Time, Timer, TimerMode},
    transform::components::Transform,
    utils::default,
    DefaultPlugins,
};
use bevy_ghx_utils::camera::{update_pan_orbit_camera, PanOrbitCamera};
use examples::{
    extend_displayed_vertices_with_container_vertice, LabelMode, TriangleDebugCursorUpdate,
    TriangleDebugPlugin, TrianglesDebugData, TrianglesDebugViewConfig, TrianglesDrawMode,
};
use ghx_constrained_delaunay::{
    constrained_triangulation::ConstrainedTriangulationConfiguration,
    debug::{DebugConfiguration, Phase, PhaseRecord},
    types::{Edge, Float, Vector3, Vertex, VertexId},
    Triangulation,
};
use serde::Deserialize;

const FRAMES_FILE: &str = "./assets/bad_apple_frames.msgpack";
const DEFAULT_FRAME_DISPLAY_PERIOD_MS: u64 = 60;
const DEFAULT_FRAME_SWITCH_MODE: FrameSwitchMode = FrameSwitchMode::Auto;

pub enum FrameSwitchMode {
    Manual,
    Auto,
}

fn main() {
    App::new()
        .add_plugins((DefaultPlugins, TriangleDebugPlugin))
        .add_systems(Startup, (setup, setup_camera))
        .add_systems(Update, (update_pan_orbit_camera, next_frame))
        // .add_systems(Update, draw_origin_circle)
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
            radius: (look_target - camera_position).length(),
            auto_orbit: false,
            ..Default::default()
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
    let frames_bytes = fs::read(FRAMES_FILE).expect("Should have been able to load frames file");
    let frames: Vec<Frame> = rmp_serde::from_slice(&frames_bytes).unwrap();

    let config = ConstrainedTriangulationConfiguration {
        debug_config: DebugConfiguration {
            phase_record: PhaseRecord::In(Phase::FilterTriangles),
            // phase_record: PhaseRecord::All,
            // force_end_at_step: Some(160),
            ..Default::default()
        },
        ..Default::default()
    };

    let mut triangulated_frames = Vec::new();
    for (i, frame) in frames.iter().enumerate() {
        info!("Loading frame n°{}", i);
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
        let mut displayed_vertices = frame
            .vertices
            .iter()
            .map(|v| Vector3::new(v.0 as Float, v.1 as Float, 0.))
            .collect();

        let plane_normal = Vector3::Z;
        // TODO Tmp debug, remove
        extend_displayed_vertices_with_container_vertice(
            &mut displayed_vertices,
            plane_normal,
            &triangulation.debug_context,
            false,
        );

        triangulated_frames.push(TriangulatedFrame {
            triangulation,
            edges,
            displayed_vertices,
        });
    }

    commands.insert_resource(TrianglesDebugData::new_with_constraintss(
        triangulated_frames[0].displayed_vertices.clone(),
        &triangulated_frames[0].edges,
        triangulated_frames[0].triangulation.debug_context.clone(),
    ));

    commands.insert_resource(TrianglesDebugViewConfig::new(
        LabelMode::None,
        TrianglesDrawMode::AllAsContourAndInteriorMeshes,
    ));

    commands.insert_resource(Frames {
        frame_index: 0,
        next_frame_timer: Timer::new(
            Duration::from_millis(DEFAULT_FRAME_DISPLAY_PERIOD_MS),
            TimerMode::Repeating,
        ),
        switch_mode: DEFAULT_FRAME_SWITCH_MODE,
        triangulated_frames,
    });
}

fn draw_origin_circle(mut gizmos: Gizmos, triangle_debug_data: Res<TrianglesDebugData>) {
    gizmos.circle(
        Vec3::new(0., 0., 0.),
        Direction3d::Z,
        0.01 * triangle_debug_data.context.scale_factor as f32,
        Color::ALICE_BLUE,
    );
}

fn next_frame(
    mut commands: Commands,
    time: Res<Time>,
    mut frames: ResMut<Frames>,
    mut debug_data_updates_events: EventWriter<TriangleDebugCursorUpdate>,
) {
    match frames.switch_mode {
        FrameSwitchMode::Manual => (),
        FrameSwitchMode::Auto => {
            frames.next_frame_timer.tick(time.delta());
            if frames.next_frame_timer.finished() {
                frames.frame_index = (frames.frame_index + 1) % frames.triangulated_frames.len();
                // info!("Swapping to frame n°{}", frames.frame_index);

                let frame_data = &frames.triangulated_frames[frames.frame_index];
                commands.insert_resource(TrianglesDebugData::new_with_constraintss(
                    frame_data.displayed_vertices.clone(),
                    &frame_data.edges,
                    frame_data.triangulation.debug_context.clone(),
                ));
                debug_data_updates_events.send(TriangleDebugCursorUpdate);
            }
        }
    }
}
