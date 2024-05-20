use bevy::prelude::*;

use bevy_ghx_utils::camera::{update_pan_orbit_camera, PanOrbitCamera};
use bevy_mod_billboard::plugin::BillboardPlugin;

use ghx_constrained_delaunay::triangulation::triangulation_from_3d_planar_vertices;
use utils::{
    create_displayed_vertices, draw_triangles_debug_data, update_triangles_debug_index,
    update_triangles_debugs_labels, TriangleDebugDataUpdated, TrianglesDebugData,
};

mod utils;

pub struct ExamplesPlugin;

impl Plugin for ExamplesPlugin {
    fn build(&self, app: &mut App) {
        app.add_systems(Startup, setup_camera);
        app.add_systems(Update, update_pan_orbit_camera);
    }
}

pub fn setup_camera(mut commands: Commands) {
    // Camera
    let camera_position = Vec3::new(0., 0., 120.5);
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

fn main() {
    App::new()
        .add_plugins((DefaultPlugins, ExamplesPlugin, BillboardPlugin))
        .add_event::<TriangleDebugDataUpdated>()
        .add_systems(Startup, setup)
        .add_systems(
            Update,
            (
                update_triangles_debug_index,
                update_triangles_debugs_labels,
                draw_triangles_debug_data,
            )
                .chain(),
        )
        .run();
}

fn setup(mut commands: Commands) {
    let mut vertices = Vec::<[f32; 3]>::new();
    vertices.push([0., 0., 0.]);
    vertices.push([0., 5., 0.]);
    vertices.push([5., 5., 0.]);
    vertices.push([5., 0., 0.]);

    // vertices.push([0., 0., 0.]);
    // vertices.push([0., 10., 0.]);
    // vertices.push([0., 5., 0.]);
    // vertices.push([5., 0., 0.]);
    // vertices.push([5., 5., 0.]);
    // vertices.push([5., 10., 0.]);

    // vertices.push([5., 5., 0.]);
    // vertices.push([0., 6., 0.]);
    // vertices.push([-4., 4., 0.]);
    // vertices.push([-5., 0., 0.]);
    // vertices.push([-3., -2., 0.]);
    // vertices.push([0., -4., 0.]);
    // vertices.push([4., -5., 0.]);
    // vertices.push([6., -3., 0.]);
    // vertices.push([7., 0., 0.]);
    // vertices.push([5., 1., 0.]);
    // vertices.push([0., -6., 0.]);
    // vertices.push([4., -6., 0.]);
    // vertices.push([0., -8., 0.]);
    // vertices.push([4., -8., 0.]);

    let plane_normal = Vec3::Z;
    let (vertex_id, debug_data) =
        triangulation_from_3d_planar_vertices(&vertices, plane_normal.into());

    let mut constrained_vertices = Vec::new();

    for vertexid in vertex_id {
        constrained_vertices.push(vertices[vertexid]);
    }

    let displayed_vertices = create_displayed_vertices(vertices, plane_normal);

    info!("vertices {:?}", constrained_vertices);

    commands.insert_resource(TrianglesDebugData {
        vertices: displayed_vertices,
        triangles_buffers: debug_data,
        current_buffer_index: 0,
    })
}
