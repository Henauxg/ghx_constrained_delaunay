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
        mesh::{Indices, Mesh, Meshable, PrimitiveTopology},
        render_asset::RenderAssetUsages,
    },
    transform::components::Transform,
    utils::default,
    DefaultPlugins,
};
use bevy_ghx_utils::camera::PanOrbitCamera;

use examples::plugin::ExamplesPlugin;
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

#[rustfmt::skip]
fn create_cube_mesh() -> Mesh {
    // Keep the mesh data accessible in future frames to be able to mutate it in toggle_texture.
    Mesh::new(PrimitiveTopology::TriangleList, RenderAssetUsages::MAIN_WORLD | RenderAssetUsages::RENDER_WORLD)
    .with_inserted_attribute(
        Mesh::ATTRIBUTE_POSITION,
        // Each array is an [x, y, z] coordinate in local space.
        // Meshes always rotate around their local [0, 0, 0] when a rotation is applied to their Transform.
        // By centering our mesh around the origin, rotating the mesh preserves its center of mass.
        vec![
            // top (facing towards +y)
            [-0.5, 0.5, -0.5], // vertex with index 0
            [0.5, 0.5, -0.5], // vertex with index 1
            [0.5, 0.5, 0.5], // etc. until 23
            [-0.5, 0.5, 0.5],
            // bottom   (-y)
            [-0.5, -0.5, -0.5],
            [0.5, -0.5, -0.5],
            [0.5, -0.5, 0.5],
            [-0.5, -0.5, 0.5],
            // right    (+x)
            [0.5, -0.5, -0.5],
            [0.5, -0.5, 0.5],
            [0.5, 0.5, 0.5], // This vertex is at the same position as vertex with index 2, but they'll have different UV and normal
            [0.5, 0.5, -0.5],
            // left     (-x)
            [-0.5, -0.5, -0.5],
            [-0.5, -0.5, 0.5],
            [-0.5, 0.5, 0.5],
            [-0.5, 0.5, -0.5],
            // back     (+z)
            [-0.5, -0.5, 0.5],
            [-0.5, 0.5, 0.5],
            [0.5, 0.5, 0.5],
            [0.5, -0.5, 0.5],
            // forward  (-z)
            [-0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [0.5, 0.5, -0.5],
            [0.5, -0.5, -0.5],
        ],
    )
    // Set-up UV coordinates to point to the upper (V < 0.5), "dirt+grass" part of the texture.
    // Take a look at the custom image (assets/textures/array_texture.png)
    // so the UV coords will make more sense
    // Note: (0.0, 0.0) = Top-Left in UV mapping, (1.0, 1.0) = Bottom-Right in UV mapping
    // .with_inserted_attribute(
    //     Mesh::ATTRIBUTE_UV_0,
    //     vec![
    //         // Assigning the UV coords for the top side.
    //         [0.0, 0.2], [0.0, 0.0], [1.0, 0.0], [1.0, 0.25],
    //         // Assigning the UV coords for the bottom side.
    //         [0.0, 0.45], [0.0, 0.25], [1.0, 0.25], [1.0, 0.45],
    //         // Assigning the UV coords for the right side.
    //         [1.0, 0.45], [0.0, 0.45], [0.0, 0.2], [1.0, 0.2],
    //         // Assigning the UV coords for the left side.
    //         [1.0, 0.45], [0.0, 0.45], [0.0, 0.2], [1.0, 0.2],
    //         // Assigning the UV coords for the back side.
    //         [0.0, 0.45], [0.0, 0.2], [1.0, 0.2], [1.0, 0.45],
    //         // Assigning the UV coords for the forward side.
    //         [0.0, 0.45], [0.0, 0.2], [1.0, 0.2], [1.0, 0.45],
    //     ],
    // )
    // For meshes with flat shading, normals are orthogonal (pointing out) from the direction of
    // the surface.
    // Normals are required for correct lighting calculations.
    // Each array represents a normalized vector, which length should be equal to 1.0.
    .with_inserted_attribute(
        Mesh::ATTRIBUTE_NORMAL,
        vec![
            // Normals for the top side (towards +y)
            [0.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            // Normals for the bottom side (towards -y)
            [0.0, -1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, -1.0, 0.0],
            // Normals for the right side (towards +x)
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            // Normals for the left side (towards -x)
            [-1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            // Normals for the back side (towards +z)
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 1.0],
            // Normals for the forward side (towards -z)
            [0.0, 0.0, -1.0],
            [0.0, 0.0, -1.0],
            [0.0, 0.0, -1.0],
            [0.0, 0.0, -1.0],
        ],
    )
    // Create the triangles out of the 24 vertices we created.
    // To construct a square, we need 2 triangles, therefore 12 triangles in total.
    // To construct a triangle, we need the indices of its 3 defined vertices, adding them one
    // by one, in a counter-clockwise order (relative to the position of the viewer, the order
    // should appear counter-clockwise from the front of the triangle, in this case from outside the cube).
    // Read more about how to correctly build a mesh manually in the Bevy documentation of a Mesh,
    // further examples and the implementation of the built-in shapes.
    .with_inserted_indices(Indices::U32(vec![
        0,3,1 , 1,3,2, // triangles making up the top (+y) facing side.
        4,5,7 , 5,6,7, // bottom (-y)
        8,11,9 , 9,11,10, // right (+x)
        12,13,15 , 13,14,15, // left (-x)
        16,19,17 , 17,19,18, // back (+z)
        20,21,23 , 21,22,23, // forward (-z)
    ]))
}
