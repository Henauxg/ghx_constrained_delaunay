use std::f32::consts::PI;

use bevy::prelude::*;

use bevy::{
    app::{App, Startup, Update},
    asset::Assets,
    core_pipeline::core_3d::Camera3dBundle,
    ecs::system::{Commands, ResMut},
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
use bevy_ghx_utils::camera::PanOrbitCamera;

use examples::plugin::ExamplesPlugin;
use examples::COLORS;
use simple_delaunay_lib::delaunay_3d::delaunay_struct_3d::DelaunayStructure3D;
use simple_delaunay_lib::delaunay_3d::simplicial_struct_3d;

fn main() {
    let mut app = App::new();

    app.insert_resource(AmbientLight {
        color: Color::WHITE,
        brightness: 2000.,
    })
    .add_plugins((DefaultPlugins, ExamplesPlugin));

    app.add_systems(Startup, setup_sandbox);
    app.add_systems(Update, generate_tera);

    app.run();
}

fn generate_tera(mut gizmos: Gizmos) {
    let cube_vertices = generate_cube();

    let mut del_struct = DelaunayStructure3D::new();
    del_struct.insert_vertices(&cube_vertices, true).unwrap();

    info!(
        "There are {} tetras",
        del_struct.get_simplicial().get_nb_tetrahedra()
    );

    for tetra_index in 0..del_struct.get_simplicial().get_nb_tetrahedra() {
        let tetra_nodes = del_struct
            .get_simplicial()
            .get_tetrahedron(tetra_index)
            .unwrap()
            .nodes();

        let mut cube_vert_indexes = vec![0; 4];
        for vert_index in 0..4 {
            cube_vert_indexes[vert_index] = match tetra_nodes[vert_index] {
                simplicial_struct_3d::Node::Infinity => {
                    error!("Infinity detected");
                    continue;
                }
                simplicial_struct_3d::Node::Value(cube_vert_index) => cube_vert_index,
            };
        }
        let color = COLORS[tetra_index % COLORS.len()];

        let mut draw_gizmo_line =
            |color: Color, cube_vert_index_1: usize, cube_vert_index_2: usize| {
                gizmos.line(
                    Vec3::new(
                        cube_vertices[cube_vert_index_1][0] as f32,
                        cube_vertices[cube_vert_index_1][1] as f32,
                        cube_vertices[cube_vert_index_1][2] as f32,
                    ),
                    Vec3::new(
                        cube_vertices[cube_vert_index_2][0] as f32,
                        cube_vertices[cube_vert_index_2][1] as f32,
                        cube_vertices[cube_vert_index_2][2] as f32,
                    ),
                    color,
                );
            };

        draw_gizmo_line(color, cube_vert_indexes[0], cube_vert_indexes[1]);
        draw_gizmo_line(color, cube_vert_indexes[0], cube_vert_indexes[2]);
        draw_gizmo_line(color, cube_vert_indexes[0], cube_vert_indexes[3]);
        draw_gizmo_line(color, cube_vert_indexes[1], cube_vert_indexes[3]);
        draw_gizmo_line(color, cube_vert_indexes[1], cube_vert_indexes[2]);
        draw_gizmo_line(color, cube_vert_indexes[3], cube_vert_indexes[2]);
    }
}

fn generate_cube() -> Vec<[f64; 3]> {
    let mut vec_pts: Vec<[f64; 3]> = Vec::new();
    vec_pts.push([0.0, 0.0, 0.0]);
    vec_pts.push([0.0, 1.0, 1.0]);
    vec_pts.push([0.0, 0.0, 1.0]);
    vec_pts.push([1.0, 0.0, 1.0]);
    vec_pts.push([1.0, 1.0, 0.0]);
    vec_pts.push([1.0, 1.0, 1.0]);
    vec_pts.push([1.0, 0.0, 0.0]);
    vec_pts.push([0.0, 1.0, 0.0]);
    vec_pts
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

    // Create and save a handle to the mesh.
    // let cube_mesh_handle: Handle<Mesh> = meshes.add(create_cube_mesh());

    // Render the mesh with the custom texture using a PbrBundle, add the marker.
    // commands.spawn((
    //     PbrBundle {
    //         mesh: cube_mesh_handle,
    //         material: materials.add(Color::ORANGE_RED),
    //         transform: Transform::from_xyz(0., 1., 0.),
    //         ..default()
    //     },
    //     // CustomUV,
    // ));
}
