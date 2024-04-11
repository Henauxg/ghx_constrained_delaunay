use std::f32::consts::PI;

use bevy::{
    app::{App, Startup},
    asset::{Assets, Handle},
    ecs::system::{Commands, ResMut},
    math::{primitives::Plane3d, EulerRot, Quat},
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
use examples::plugin::ExamplesPlugin;
use tritet::Tetgen;

fn main() {
    let mut app = App::new();

    app.insert_resource(AmbientLight {
        color: Color::WHITE,
        brightness: 2000.,
    })
    .add_plugins((DefaultPlugins, ExamplesPlugin));

    app.add_systems(Startup, setup_sandbox);

    app.run();
}

pub fn setup_sandbox(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    // mut destruction_mesh: ResMut<>,
) {
    // // Plane
    commands.spawn(PbrBundle {
        mesh: meshes.add(Plane3d::default().mesh().size(500000.0, 500000.0)),
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
    let cube_mesh_handle: Handle<Mesh> = meshes.add(create_cube_mesh());

    // Render the mesh with the custom texture using a PbrBundle, add the marker.
    commands.spawn((
        PbrBundle {
            mesh: cube_mesh_handle,
            // material: materials.add(StandardMaterial {
            //      base_color_texture: Some(custom_texture_handle),
            //     ..default()
            // }),
            material: materials.add(Color::ORANGE_RED),
            transform: Transform::from_xyz(0., 1., 0.),
            ..default()
        },
        // CustomUV,
    ));

    // allocate data for 8 points
    let mut delaunay = Tetgen::new(8, None, None, None).unwrap();

    // set points
    delaunay
        .set_point(0, 0, 0.0, 0.0, 0.0)
        .unwrap()
        .set_point(1, 0, 1.0, 0.0, 0.0)
        .unwrap()
        .set_point(2, 0, 1.0, 1.0, 0.0)
        .unwrap()
        .set_point(3, 0, 0.0, 1.0, 0.0)
        .unwrap()
        .set_point(4, 0, 0.0, 0.0, 1.0)
        .unwrap()
        .set_point(5, 0, 1.0, 0.0, 1.0)
        .unwrap()
        .set_point(6, 0, 1.0, 1.0, 1.0)
        .unwrap()
        .set_point(7, 0, 0.0, 1.0, 1.0)
        .unwrap();

    // generate Delaunay triangulation
    delaunay.generate_delaunay(false).unwrap();

    // let ntet = tetgen_data.out_ncell();
    const EDGES: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
    let tetras_count = delaunay.out_ncell();
    // let attributes = Vec::with_capacity(tetras_count);
    for tetra_idx in 0..tetras_count {
        // attributes.push(delaunay.out_cell_attribute(tetra_idx));
        for vert_local_index in 0..4 {
            let vert_id = delaunay.out_cell_point(tetra_idx, vert_local_index);
            let xyz = vec![
                delaunay.out_point(vert_id, 0),
                delaunay.out_point(vert_id, 1),
                delaunay.out_point(vert_id, 2),
            ];
            // for dim in 0..3 {
            //     x[dim] = tetgen_data.out_point(vert_id, dim);
            //     min[dim] = f64::min(min[dim], x[dim]);
            //     max[dim] = f64::max(max[dim], x[dim]);
            //     xcen[dim] += x[dim] / 4.0;
            // }
        }
        // for (vert_local_idx_a, vert_local_idx_b) in &EDGES {
        //     let a = delaunay.out_cell_point(tetra_idx, *vert_local_idx_a);
        //     let b = delaunay.out_cell_point(tetra_idx, *vert_local_idx_b);
        //     for dim in 0..3 {
        //         xa[dim] = self.out_point(a, dim);
        //         xb[dim] = self.out_point(b, dim);
        //     }
        //     canvas.polyline_3d_begin();
        //     canvas.polyline_3d_add(xa[0], xa[1], xa[2]);
        //     canvas.polyline_3d_add(xb[0], xb[1], xb[2]);
        //     canvas.polyline_3d_end();
        // }
    }
    let ntet = delaunay.out_ncell();

    // commands.insert_resource(DestructionMesh { delaunay });

    // tetgen_data.gene
}

// #[derive(Resource)]
// pub struct DestructionMesh {
//     delaunay: Arc<Tetgen>,
// }

// impl Send for Arc<Mutex<Tetgen>> {}

// pub fn draw_destrtuction_mesh(destr_mesh: Res<DestructionMesh>, mut gizmos: Gizmos) {}

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
